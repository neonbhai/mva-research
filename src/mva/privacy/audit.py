"""The privacy audit engine.

Every check answers one question about one concrete failure mode, and every check
reports **paths, line numbers, span lengths and rule IDs only** (GP-41). The audit
is run by an agent and by CI; its output lands in a model context window and a
build log, so a check that echoed what it found would be the largest leak in the
system.

The check set is chosen around how patient data actually escapes a research repo,
which is almost never "someone published a VCF on purpose":

* it is already **tracked**, so ``.gitignore`` is irrelevant to it
  (``git_tracked_sensitive``);
* it is **staged right now**, and the worktree copy has since been cleaned
  (``git_staged_sensitive`` reads the staged blob, not the file on disk);
* the ignore rule that was supposed to catch it is **cancelled by a negation**
  further down the file (``gitignore_effectiveness``, ``gitignore_negation_safety``);
* the workspace is a **symlink** into the repo, or a symlink in the repo points
  out into the workspace (``workspace_containment``, ``symlink_escape``);
* it is in a file whose extension says nothing (``content_scan``);
* it is in the FILE NAME rather than the file — a sequencing filename is
  routinely the MRN, the NHS number or the accession, so paths are scanned as
  well as scanned-for, and every path is redacted before it reaches the report
  (``_path_findings``, :func:`redact_path`);
* it is inside a gzip member, where every plaintext rule is blind
  (``_scan_gzip_member``);
* it is under an allowlisted prefix, where softening a finding used to be
  automatic. It is not: see :data:`PATH_DOWNGRADABLE_RULES` and
  :data:`DECLARED_SYNTHETIC_DOWNGRADABLE_RULES`;
* a real file was dropped into the synthetic-fixture directory, where the
  ``.gitignore`` negations deliberately re-admit files
  (``synthetic_fixtures_marked``);
* a real file was dropped into a *public reference* fixture directory, where the
  negations re-admit ``.vcf.gz`` and ``.tbi`` outright. That category is granted
  only against declared, hash-pinned provenance AND a verification that the exact
  buffer being scanned is a single sites-only VCF -- "sites-only" is not
  "not derived from an individual", and a sites-only export of a proband is still
  that child's complete variant list (``reference_fixtures_declared``,
  :func:`buffer_is_sites_only_vcf`, ADR 0012);
* it is in a notebook output cell (``notebook_output_purity``);
* it went to a log (``log_redaction_probe``, which is a live runtime probe, not a
  static check — the only way to know redaction works is to try to defeat it);
* the workspace is outside the repo but inside iCloud Drive
  (``cloud_sync_location``).
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import zlib
from bisect import bisect_right
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final, cast

import yaml

from mva.config import (
    CLOUD_SYNCED_HOME_DIRS,
    CLOUD_SYNCED_MARKERS,
    path_is_within,
)
from mva.determinism import hash_file
from mva.errors import PrivacyViolationError
from mva.privacy.classify import is_sensitive_extension
from mva.privacy.patterns import (
    BINARY_SEVERITY,
    HEAD_BYTES,
    HPO_DISTINCT_FAIL_THRESHOLD,
    MAX_SCAN_BYTES,
    RULES,
    Severity,
    correlation_id,
    decode_lossy,
    decode_scrubbed,
    gunzip_capped,
    read_capped,
    sniff_binary,
)
from mva.privacy.redact import (
    redact_text,
    redaction_installed,
    unfiltered_handlers,
)

# ---------------------------------------------------------------------------
# Policy constants
# ---------------------------------------------------------------------------

#: Tracked files here may carry a patient-data EXTENSION without failing
#: ``git_tracked_sensitive``. Only the synthetic fixture directory qualifies: it is
#: the one place ``.gitignore`` deliberately re-admits ``.vcf``/``.ped``/``.bed``,
#: and it is the one place with a marker check (``synthetic_fixtures_marked``)
#: standing behind the exemption. ``knowledge/public/`` used to be here too and is
#: not any more — it holds curated public TSVs, none of which has a sensitive
#: extension, so the exemption bought nothing and would have waved through a
#: ``knowledge/public/proband.vcf``.
TRACKED_EXEMPT_PREFIXES: Final[tuple[str, ...]] = ("tests/fixtures/synthetic/",)

#: Paths where genomic-looking CONTENT can be legitimate. Being on this list is
#: not by itself permission to look like patient data — see
#: :func:`_resolve_context` for exactly which rules it can soften and under what
#: additional condition.
CONTENT_ALLOWLIST_PREFIXES: Final[tuple[str, ...]] = (
    "tests/fixtures/synthetic/",
    "knowledge/public/",
    "knowledge/manifests/",
    "tests/golden/",
)

#: The ONLY rules an allowlisted path may soften on the strength of its path
#: alone. Both are the documented false-positive-prone rules: a public coordinate
#: table really is VCF-shaped, and a public ontology table really is full of HPO
#: identifiers. Nothing else is on this list, and the omissions are the point.
#:
#: Until this existed, ``_resolve_context`` downgraded EVERY ``fail`` to ``warn``
#: for any path under the four allowlisted prefixes. That is an unconditional
#: bypass: a complete, real VCF committed to ``knowledge/public/variants.tsv`` or
#: ``tests/golden/expected.csv`` produced twelve rule hits, every one of them
#: downgraded, and the audit reported ``passed=True``. ``vcf_header``,
#: ``vcf_chrom_line``, ``genotype_field``, ``fastq_record``, ``sam_rg_sample``,
#: ``fasta_record`` and ``plink_ped_line`` have near-zero false-positive risk and
#: are precisely what a real file trips.
PATH_DOWNGRADABLE_RULES: Final[frozenset[str]] = frozenset({"vcf_data_line", "hpo_term"})

#: Softened for a PUBLIC REFERENCE SLICE that structurally carries no sample
#: columns (ADR 0012). A sites-only VCF necessarily looks like a VCF — that is
#: what it is — and it has no individual in it to disclose. ``genotype_field`` is
#: included for a specific reason: when the ``#CHROM`` line declares no sample
#: columns, a genotype-shaped match CANNOT be a genotype, so it is a pattern
#: coincidence in INFO text rather than a finding.
#:
#: The ``magic:*`` container rules are deliberately absent, keeping full strength
#: everywhere, per the reasoning in :func:`_resolve_context`.
SITES_ONLY_DOWNGRADABLE_RULES: Final[frozenset[str]] = frozenset(
    {"vcf_header", "vcf_chrom_line", "vcf_data_line", "genotype_field"}
)

#: Rules an allowlisted path may soften ONLY when the file also carries a verified
#: synthetic declaration (see :func:`declares_synthetic`). This is what keeps
#: ``tests/fixtures/synthetic/synthetic_case.vcf`` — a deliberately VCF-shaped
#: fixture — from failing the audit forever, without extending the same courtesy
#: to an undeclared file that merely happens to sit in the same directory.
DECLARED_SYNTHETIC_DOWNGRADABLE_RULES: Final[frozenset[str]] = frozenset(
    {
        "vcf_header",
        "vcf_chrom_line",
        "genotype_field",
        "fastq_record",
        "sam_rg_sample",
        "fasta_record",
        "plink_ped_line",
    }
)

#: Byte markers a file uses to declare itself synthetic. Accepted only in the head
#: of the file (:data:`SYNTHETIC_MARKER_BYTES`), never buried in the body.
SYNTHETIC_MARKERS: Final[tuple[bytes, ...]] = (b"mva_synthetic=true", b"SYNTH_", b"SYNTH-")

#: How far into a file a synthetic marker is looked for. A declaration is a header,
#: not a needle: allowing it anywhere in 8 MiB meant any file containing the
#: substring ``SYNTH_`` once, at any depth, counted as declared.
SYNTHETIC_MARKER_BYTES: Final[int] = 4096

#: Prefixes where an HPO identifier is a CONSTANT rather than a phenotype profile:
#: reviewed source, its tests, docs and prompt briefs. This downgrade is specific
#: to ``hpo_term`` — every other fail-severity rule keeps full strength here.
#:
#: The reasoning: an HPO ID is a public ontology term. It becomes disclosive only
#: as part of a patient's *profile*, and by GP-40 a patient profile lives in the
#: external workspace and never in the repository at all — which is enforced by
#: ``workspace_containment`` and ``git_tracked_sensitive``, not by counting terms
#: in a ``.py`` file. Without this, a docstring showing the canonical ``HP:``
#: format and a test fixture listing three phenotypes both fail forever, and an
#: audit that is permanently red is an audit that gets switched off.
#: ``knowledge/``, ``tools/`` and ``docs/`` are added for the same reason and
#: under the same rule as ``src/``: they hold gene→phenotype panel tables, the HPO
#: frequency subontology (``HP:0040280`` to ``HP:0040284``, which is ontology
#: *metadata*), and reviewed prose about both. What makes a cluster of distinct
#: HPO terms identifying is that it is **tied to a subject**; a table keyed by gene,
#: by disease or by ontology term has no subject to identify. That danger model is
#: NOT weakened here — see :func:`has_phenotype_profile_shape`, which is the other
#: half of the condition and is checked against the bytes.
HPO_CODE_PREFIXES: Final[tuple[str, ...]] = (
    "src/",
    "tests/",
    "docs/",
    "prompts/",
    "knowledge/",
    "tools/",
)

#: Clinical observation statuses (`mva.models.phenotype.PhenotypeStatus`) as a
#: whole DELIMITED FIELD. A term carrying one of these has been *asserted about
#: somebody*, which a gene- or ontology-keyed reference table never does.
#:
#: Anchored to field boundaries, not to ``\b``: a bare word match hit every
#: docstring containing "observed" and every prose line containing "excluded",
#: which is most of the phenotype package.
_PHENOTYPE_STATUS_FIELD: Final = re.compile(
    rb"(?:^|[\t,])(?:observed|excluded|not_assessed)(?:[\t,]|$)"
)

#: The header row of a loadable phenotype profile
#: (`mva.phenotype.loader.REQUIRED_COLUMNS` = ``hpo_id``, ``label``, ``status``).
#: Anchored to the START of a line and followed by a delimiter, so the mere
#: MENTION of ``hpo_id`` — in the loader, in a builder, in a test, in a docstring —
#: is not a profile. That distinction is the whole check: code that reads the
#: format names the column; only the data itself has it as a header.
_PHENOTYPE_HEADER_ROW: Final = re.compile(rb"^[\"']?hpo_id[\"']?[\t,]")

_HPO_TERM_BYTES: Final = re.compile(rb"HP:\d{7}")

#: How many HPO terms must sit on status-bearing DATA ROWS before the file is
#: treated as a profile. One is a code constant or a doc example; a cluster of
#: status-annotated terms is a person.
_PROFILE_SHAPE_MIN_TERMS: Final[int] = 3

#: Minimum delimiters for a line to count as a data row rather than prose. A
#: profile record is ``hpo_id``, ``label``, ``status`` at least.
_PROFILE_ROW_MIN_DELIMITERS: Final[int] = 2


def has_phenotype_profile_shape(data: bytes) -> bool:
    """Whether these bytes look like a phenotype profile keyed to ONE individual.

    This is the half of the ``hpo_term`` path downgrade that does the real work, and
    it is deliberately a property of the bytes rather than of the path — the same
    reasoning ADR 0012 applies to reference fixtures: provenance is a claim in a
    comment, shape is checkable.

    A phenotype profile is not "a file with HPO terms in it". It is a table of terms
    each carrying a **clinical observation status** — ``observed`` / ``excluded`` /
    ``not_assessed`` — which is an assertion about a person. A gene→phenotype panel
    table, the HPO frequency subontology and a docstring showing the ``HP:`` format
    all lack that, and none of them has a subject to identify.

    Two independent signals, either sufficient:

    * a header ROW declaring the loader's required columns; or
    * at least :data:`_PROFILE_SHAPE_MIN_TERMS` HPO terms on delimited data rows
      that also carry an observation status as a field.

    Both are anchored to line and field structure rather than to substrings. The
    unanchored version of this check flagged eleven files whose only offence was
    naming the format — ``src/mva/phenotype/hpo.py``, the HPO fixture builders,
    every phenotype test — which is the false-positive spiral that gets a scanner
    switched off. Measured over the whole tree, the anchored version flags exactly
    one file: ``tests/fixtures/synthetic/synthetic_phenotype.tsv``, which really is
    a phenotype profile (a declared-synthetic one, downgraded by the allowlist).

    Note what this does NOT rely on: a subject identifier in the file. The real
    proband profile does not contain one — ``subject_id`` is passed to the loader as
    an argument — so keying the check on a visible patient ID would have been a
    check that fails open on the exact file it exists for.
    """
    lines = data.splitlines()
    for line in lines:
        if _PHENOTYPE_HEADER_ROW.match(line) and b"label" in line and b"status" in line:
            return True

    seen = 0
    for line in lines:
        delimiters = line.count(b"\t") + line.count(b",")
        if delimiters < _PROFILE_ROW_MIN_DELIMITERS:
            continue
        if _HPO_TERM_BYTES.search(line) is None or _PHENOTYPE_STATUS_FIELD.search(line) is None:
            continue
        seen += len(_HPO_TERM_BYTES.findall(line))
        if seen >= _PROFILE_SHAPE_MIN_TERMS:
            return True
    return False


#: Slices of PUBLIC reference releases committed as adapter test fixtures.
#: These hold real data — that is the point, since a fixture of invented records
#: would not test the parser against what the source actually emits — but real
#: *public reference* data, which has no individual to identify. See ADR 0012.
#:
#: Membership here is NOT by itself permission, and is not even half of it. Three
#: further conditions are ANDed with it at every use site:
#:
#: * the file is DECLARED in its directory's :data:`FIXTURE_PROVENANCE_FILENAME`,
#:   naming the manifest resource it was cut from and the regions that were cut
#:   (ADR 0012 conditions 2 and 4);
#: * that resource is registered in :data:`RESOURCE_MANIFEST_RELPATH`, is
#:   ``status: fetched`` and carries a real sha256 (condition 2);
#: * the exact BUFFER being scanned is a single sites-only VCF (condition 3),
#:   which is checked by :func:`buffer_is_sites_only_vcf` against those bytes and
#:   never against a path reopened independently of them.
PUBLIC_REFERENCE_FIXTURE_PREFIXES: Final[tuple[str, ...]] = (
    "tests/fixtures/clinvar/",
    "tests/fixtures/gnomad/",
    "tests/fixtures/mane/",
    "tests/fixtures/hpo/",
)

#: A conforming VCF header line is EXACTLY these eight columns, in this order,
#: with nothing after ``INFO``. Compared as a whole tuple rather than counted: a
#: length test accepts ``#CHROM POS ID REF ALT QUAL FILTER SAMPLE01`` and a
#: startswith test accepts ``#CHROMOSOME``.
_VCF_FIXED_COLUMNS: Final[tuple[bytes, ...]] = (
    b"#CHROM",
    b"POS",
    b"ID",
    b"REF",
    b"ALT",
    b"QUAL",
    b"FILTER",
    b"INFO",
)

#: A conforming VCF header line has 8 fixed columns before FORMAT/samples.
_VCF_FIXED_FIELDS: Final[int] = len(_VCF_FIXED_COLUMNS)

#: Bound on the gzip member chain :func:`inflate_fully` will walk. BGZF is a chain
#: of ~64 KiB members, so a few hundred covers any committable fixture; the bound
#: exists so a crafted chain cannot spin. Hitting it is NOT truncation-with-a-
#: shrug: :func:`inflate_fully` returns ``None``, which denies the exemption.
_MAX_INFLATE_MEMBERS: Final[int] = 4096


def inflate_fully(data: bytes, limit: int = MAX_SCAN_BYTES) -> bytes | None:
    """Inflate a gzip/BGZF byte string COMPLETELY, or return ``None``.

    Deliberately not :func:`~mva.privacy.patterns.gunzip_capped`. That function is
    built for scanning, where a truncated prefix is still worth matching rules
    against, so it returns whatever it managed to inflate and says nothing about
    whether it reached the end. For an exemption decision that is exactly wrong: a
    prefix that contains no ``FORMAT`` column proves nothing about the member
    nobody inflated. So this returns ``None`` — deny — for a stream that is
    truncated, corrupt, longer than ``limit``, or made of more members than
    :data:`_MAX_INFLATE_MEMBERS`.

    Concatenated members are followed to the end, because a second gzip member is
    the cheapest way to hide a genotyped VCF behind a sites-only one.
    """
    if data[:2] != b"\x1f\x8b":
        return None
    out = bytearray()
    remaining = data
    for _ in range(_MAX_INFLATE_MEMBERS):
        obj = zlib.decompressobj(16 + zlib.MAX_WBITS)
        try:
            # One byte of headroom past `limit`, so an over-long stream is
            # detectable rather than silently trimmed to the cap.
            out += obj.decompress(remaining, max(0, limit - len(out)) + 1)
        except zlib.error:
            return None
        if len(out) > limit or not obj.eof:
            return None
        if not obj.unused_data:
            return bytes(out)
        remaining = obj.unused_data
    return None


def _plaintext_is_sites_only_vcf(text: bytes) -> bool:
    """Whether this DECOMPRESSED buffer is one complete, sites-only VCF.

    ADR 0012 condition 3, applied to the whole buffer rather than to its first
    header line. The predicate this replaces returned at the first ``#CHROM``,
    which a file could satisfy with an eight-column decoy and then follow with a
    second header carrying ``FORMAT`` and a patient sample. Every rule here exists
    to close a specific way of doing that:

    * **exactly one** ``#CHROM`` line in the whole buffer — a second one is a
      concatenation, whatever the first one said;
    * that line is exactly :data:`_VCF_FIXED_COLUMNS` — no ``FORMAT``, no sample
      names, no trailing tab;
    * no ``##`` meta line appears AFTER the header — that is the shape of two
      VCFs stapled together;
    * no other ``#``-prefixed line anywhere — unrecognised is not benign;
    * every data row has exactly eight tab-separated fields, so a row cannot be
      wider than the header that describes it;
    * no data row precedes the header.

    Blank lines are skipped; nothing else is. Returns ``False`` for a buffer with
    no ``#CHROM`` line at all, which is how every non-VCF fixture (the MANE GTF,
    the HPO OBO, a ``.tbi``) correctly fails to earn the exemption.
    """
    headers = 0
    for raw in text.split(b"\n"):
        line = raw.rstrip(b"\r")
        if not line:
            continue
        if line.startswith(b"##"):
            if headers:
                return False
            continue
        if line.startswith(b"#CHROM"):
            if headers or tuple(line.split(b"\t")) != _VCF_FIXED_COLUMNS:
                return False
            headers = 1
            continue
        if line.startswith(b"#"):
            return False
        if headers != 1 or line.count(b"\t") != _VCF_FIXED_FIELDS - 1:
            return False
    return headers == 1


def buffer_is_sites_only_vcf(data: bytes) -> bool:
    """ADR 0012 condition 3, decided against **the exact bytes being scanned**.

    This is the whole guarantee, and taking it from the right buffer is half of
    it. The privacy risk being managed is not "real data is in the repo" — ClinVar
    and gnomAD are public and already sit on every clinical genomics machine. It
    is **"a file carrying per-sample genotypes got committed"**, because that is
    what discloses an individual.

    Callers pass the buffer they are about to scan or about to commit — a staged
    blob, an index blob, a decompressed member — never a path. A path reopened
    beside a scan is a second file: staging a genotyped fixture and then cleaning
    the worktree copy made the audit probe one and report on the other.

    **Fails closed** on everything else: a buffer at or over
    :data:`~mva.privacy.patterns.MAX_SCAN_BYTES` may have been truncated by the
    read cap and cannot be cleared; a gzip stream that will not fully inflate
    cannot be cleared; a gzip inside a gzip is ambiguous and is not cleared.
    """
    if len(data) >= MAX_SCAN_BYTES:
        return False
    if data[:2] == b"\x1f\x8b":
        plain = inflate_fully(data)
        if plain is None or plain[:2] == b"\x1f\x8b":
            return False
        return _plaintext_is_sites_only_vcf(plain)
    return _plaintext_is_sites_only_vcf(data)


def vcf_declares_sample_columns(path: Path) -> bool:
    """Whether the file at ``path`` fails to be a verified sites-only VCF.

    A path-shaped convenience over :func:`buffer_is_sites_only_vcf`, kept for
    callers that genuinely hold nothing but a path. The audit's own checks do NOT
    use it: they hold the buffer that is about to be committed or scanned, and
    reopening the path beside that buffer is the bug this module was fixed for.

    **Fails closed**: unreadable, missing, malformed or not-a-VCF all report
    ``True``.
    """
    try:
        data = read_capped(path, MAX_SCAN_BYTES)
    except OSError:
        return True
    return not buffer_is_sites_only_vcf(data)


# ---------------------------------------------------------------------------
# ADR 0012 conditions 1, 2 and 4: verifiable provenance
# ---------------------------------------------------------------------------

#: The committed, hash-pinned index of public releases this project fetched.
#: Pinned here as a CONSTANT and never read out of a fixture's own metadata: a
#: directory allowed to name its own manifest could ship that manifest too, which
#: is precisely the general-purpose escape hatch condition 2 exists to close.
RESOURCE_MANIFEST_RELPATH: Final[str] = "knowledge/manifests/resources.yaml"

#: One provenance record per fixture directory, naming — per committed file — the
#: manifest resource it was cut from and the regions that were cut.
#:
#: Why a per-directory YAML sidecar rather than an in-band marker (the shape the
#: synthetic category uses): there is no in-band place to put it. The MANE GTF has
#: no comment lines at all, a ``.tbi`` is binary, and ``hp.obo`` is a third
#: format — an in-band claim would need a different, forgeable mechanism per file
#: type. A sidecar is also the only form that shows up in review as a diff line:
#: admitting a new file to the category means writing its name and its claimed
#: descent where a reviewer reads them.
#:
#: What this does and does not buy is worth stating plainly. Naming a resource is
#: a CLAIM, and a claim can be false — nothing offline can prove a given slice was
#: cut from ClinVar. What the requirement removes is the silent path: a file can
#: no longer enter the category by being dropped in a directory. It has to be
#: named, beside a generator, against a resource the manifest records fetching.
#: The property that actually protects the proband is still condition 3, which is
#: a fact about the bytes.
FIXTURE_PROVENANCE_FILENAME: Final[str] = "fixture_provenance.yaml"

#: A manifest entry qualifies only with a REAL digest. ``sha256: null`` is what
#: the manifest carries for a resource that was never successfully fetched, and
#: a placeholder string is what a hand-edit produces.
_SHA256_HEX: Final = re.compile(r"\A[0-9a-f]{64}\Z")


@dataclass(frozen=True, slots=True)
class FixtureProvenance:
    """One committed fixture's declared descent.

    Exactly one of ``resource`` (a data fixture cut from a manifest resource) and
    ``indexes`` (a sidecar inheriting from a sibling data fixture) is set.

    An index holds byte offsets, never records, so it cannot disclose anything its
    data file does not; it therefore inherits that file's verdict rather than
    being judged on its own bytes, where the sites-only predicate would
    (correctly, but uselessly) refuse it. That inheritance is DECLARED rather than
    inferred from the file name -- ``x.vcf.gz.tbi`` is not evidence that
    ``x.vcf.gz`` is what it indexes, and a sidecar whose named companion is absent
    from the source under audit inherits nothing.
    """

    resource: str | None
    #: Bare file name of the sibling data fixture this sidecar indexes.
    indexes: str | None
    #: What was cut from the release, for a human: intervals, or the selection
    #: rule where the release has no intervals. ADR 0012 condition 4's record.
    regions: str


@dataclass(frozen=True, slots=True)
class FixtureProvenanceDoc:
    """A validated ``fixture_provenance.yaml``."""

    generator: str
    fixtures: dict[str, FixtureProvenance]


@dataclass(frozen=True, slots=True)
class ReferenceFixtureVerdict:
    """Whether ADR 0012 grants a path the public-reference exemption, and why not.

    The reason is a fixed structural string, never an echo of a name, a resource
    key or a matched byte: it is rendered into audit findings, and GP-41 says a
    finding carries paths, counts and rule IDs only.
    """

    exempt: bool
    reason: str


def _load_yaml_mapping(data: bytes | None) -> dict[str, object] | None:
    """Parse one small YAML document into a mapping, or ``None`` for anything else."""
    if data is None:
        return None
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return None
    try:
        doc: object = yaml.safe_load(text)
    except yaml.YAMLError:
        return None
    if not isinstance(doc, dict):
        return None
    return {str(key): value for key, value in cast(dict[object, object], doc).items()}


#: Memo for the parsed resource manifest, keyed on a digest of the BYTES rather
#: than on a path: the manifest now arrives through the same loader as everything
#: else, so two sources (worktree, index) can legitimately hold two versions of it
#: in one audit run. A 50 KiB YAML parse per fixture file would otherwise dominate
#: the check.
_MANIFEST_CACHE: dict[bytes, dict[str, dict[str, object]] | None] = {}
_MANIFEST_CACHE_MAX: Final[int] = 8


def _manifest_resources(data: bytes | None) -> dict[str, dict[str, object]] | None:
    """The ``resources:`` block of the committed manifest, or ``None`` if unusable."""
    if data is None:
        return None
    key = hashlib.blake2b(data, digest_size=16).digest()
    if key in _MANIFEST_CACHE:
        return _MANIFEST_CACHE[key]
    parsed: dict[str, dict[str, object]] | None = None
    doc = _load_yaml_mapping(data)
    if doc is not None:
        raw = doc.get("resources")
        if isinstance(raw, dict):
            parsed = {}
            for name, entry in cast(dict[object, object], raw).items():
                if isinstance(entry, dict):
                    parsed[str(name)] = {
                        str(k): v for k, v in cast(dict[object, object], entry).items()
                    }
    if len(_MANIFEST_CACHE) >= _MANIFEST_CACHE_MAX:
        _MANIFEST_CACHE.clear()
    _MANIFEST_CACHE[key] = parsed
    return parsed


def _is_bare_name(value: object) -> bool:
    """A non-empty file name with no path structure in it."""
    return (
        isinstance(value, str)
        and bool(value)
        and "/" not in value
        and "\\" not in value
        and value not in {".", ".."}
    )


def _regions_text(value: object) -> str:
    """Flatten a ``regions:`` declaration to one line, or ``""`` if unusable."""
    if isinstance(value, str):
        return value.strip()
    if isinstance(value, list):
        parts = [str(item).strip() for item in cast(list[object], value)]
        return "; ".join(part for part in parts if part)
    return ""


def _parse_fixture_entry(raw_entry: object) -> tuple[FixtureProvenance | None, str]:
    """Validate one ``fixtures:`` entry. ``(None, reason)`` for anything malformed."""
    if not isinstance(raw_entry, dict):
        return None, "a fixture entry is not a mapping"
    entry = {str(k): v for k, v in cast(dict[object, object], raw_entry).items()}
    resource = entry.get("resource")
    indexes = entry.get("indexes")
    if (resource is None) == (indexes is None):
        return None, "a fixture entry does not declare exactly one of `resource:` or `indexes:`"
    if resource is not None:
        if not isinstance(resource, str) or not resource.strip():
            return None, "a `resource:` is not a manifest resource name"
        regions = _regions_text(entry.get("regions"))
        if not regions:
            return None, "a data fixture does not record `regions:` (condition 4)"
        return FixtureProvenance(resource, None, regions), "declared"
    if not _is_bare_name(indexes):
        return None, "an `indexes:` does not name a sibling fixture by bare file name"
    return FixtureProvenance(None, cast(str, indexes), ""), "declared"


def _load_fixture_provenance(
    prefix: str, *, load: Callable[[str], bytes | None]
) -> tuple[FixtureProvenanceDoc | None, str]:
    """Parse and VALIDATE one fixture directory's provenance record.

    Read through ``load`` like every other input, and for the same reason: the
    record is part of the exemption decision, so taking it from the worktree while
    judging an index blob would recreate the staged-versus-worktree divergence one
    level up. A directory can legitimately have one record on disk and a different
    one staged.

    Every deviation is fatal for the whole directory and returns ``(None, reason)``.
    A record that is partly wrong is a record nobody is maintaining, and the
    failure direction that matters here is refusing an exemption, not granting one.
    The reason is a fixed string — see :class:`ReferenceFixtureVerdict`.
    """
    doc = _load_yaml_mapping(load(prefix + FIXTURE_PROVENANCE_FILENAME))
    if doc is None:
        return None, f"no readable {FIXTURE_PROVENANCE_FILENAME} beside it"
    if doc.get("manifest") != RESOURCE_MANIFEST_RELPATH:
        return None, f"the record does not point `manifest:` at {RESOURCE_MANIFEST_RELPATH}"
    generator = doc.get("generator")
    if not _is_bare_name(generator):
        return None, "`generator:` is missing or is not a bare file name"
    generator = cast(str, generator)
    if load(prefix + generator) is None:
        return None, "the declared generator is not present beside the fixtures (condition 4)"
    raw_fixtures = doc.get("fixtures")
    if not isinstance(raw_fixtures, dict):
        return None, "`fixtures:` is missing or is not a mapping"

    parsed: dict[str, FixtureProvenance] = {}
    for raw_name, raw_entry in cast(dict[object, object], raw_fixtures).items():
        if not _is_bare_name(raw_name):
            return None, "a fixture key is not a bare file name"
        entry, reason = _parse_fixture_entry(raw_entry)
        if entry is None:
            return None, reason
        parsed[str(raw_name)] = entry

    for entry_value in parsed.values():
        if entry_value.indexes is None:
            continue
        target = parsed.get(entry_value.indexes)
        if target is None or target.resource is None:
            return None, "an `indexes:` names a file that is not itself a declared data fixture"
    return FixtureProvenanceDoc(generator, parsed), "declared"


@dataclass(frozen=True, slots=True)
class PinnedResource:
    """One manifest entry that is actually usable as a fixture's ancestor."""

    name: str
    #: Path relative to the external resource root, as the manifest records it.
    path: str
    sha256: str


def pinned_resource(repo_root: Path, resource: str) -> PinnedResource:
    """The manifest entry ``resource`` names, if it is a real hash pin.

    Raises ``LookupError`` when the manifest is unreadable, the resource is not
    registered, it is not ``status: fetched``, or its digest is not a real
    sha256 — the same four conditions the audit refuses an exemption for, so the
    generator and the audit cannot drift into disagreeing about what "pinned"
    means.
    """
    resources = _manifest_resources(worktree_loader(repo_root)(RESOURCE_MANIFEST_RELPATH))
    if resources is None:
        raise LookupError(f"cannot read {RESOURCE_MANIFEST_RELPATH} under {repo_root}")
    record = resources.get(resource)
    if record is None:
        raise LookupError(f"{resource!r} is not registered in {RESOURCE_MANIFEST_RELPATH}")
    if record.get("status") != "fetched":
        raise LookupError(f"{resource!r} is not recorded as fetched; acquire it first")
    digest = record.get("sha256")
    if not isinstance(digest, str) or _SHA256_HEX.match(digest) is None:
        raise LookupError(f"{resource!r} carries no real sha256 in the manifest")
    path = record.get("path")
    if not isinstance(path, str) or not path:
        raise LookupError(f"{resource!r} records no path in the manifest")
    return PinnedResource(resource, path, digest)


def pinned_source(
    repo_root: Path,
    resource: str,
    *,
    source: Path | None = None,
    resource_root: Path | None = None,
) -> Path:
    """Locate and HASH-VERIFY the release a fixture generator is about to cut from.

    ADR 0012 condition 2 is a claim the privacy audit *checks* and this function
    is where a generator *makes it true*. The audit can only see that a fixture
    names a resource the manifest records fetching; nothing in the committed bytes
    proves the slice came from that release. Verifying the input here is what puts
    something real behind the claim.

    It closes the concrete hole in generators that took ``--source`` and read
    whatever was there: pointed at the proband's VCF, such a generator produced a
    file the audit would then clear on the strength of a ``resource: clinvar_vcf``
    line, and nothing in the repository would disagree. ``--source`` still exists,
    because a release may legitimately live outside the resource root — but it now
    has to hash to the pinned digest like everything else.

    ``resource_root`` is a parameter rather than an environment lookup because
    where the acquisition tool puts its downloads is that tool's business
    (``tools.acquire.fetch.resolve_resource_root``), and this package must not
    grow a second opinion about it.
    """
    pin = pinned_resource(repo_root, resource)
    candidate = source if source is not None else None
    if candidate is None:
        if resource_root is None:
            raise LookupError(
                f"pass --source, or --resource-root / $MVA_RESOURCES so {pin.path!r} can be found"
            )
        candidate = resource_root / pin.path
    if not candidate.is_file():
        raise LookupError(f"the pinned release for {resource!r} is not at the path given")
    actual = hash_file(candidate)
    if actual != pin.sha256:
        raise LookupError(
            f"sha256 mismatch for {resource!r}: the manifest pins {pin.sha256}, "
            f"the file given hashes to {actual}. Refusing to cut a public-reference "
            "fixture from bytes that are not the pinned release."
        )
    return candidate


def worktree_loader(repo_root: Path) -> Callable[[str], bytes | None]:
    """Read repo-relative paths from the WORKTREE, capped and memoised.

    ``None`` for anything that is not a readable regular file — a symlink
    included, since following one out of the repository is how the containment
    checks get defeated.
    """
    seen: dict[str, bytes | None] = {}

    def load(path: str) -> bytes | None:
        if path not in seen:
            candidate = repo_root / path
            try:
                data = (
                    None
                    if candidate.is_symlink() or not candidate.is_file()
                    else read_capped(candidate, MAX_SCAN_BYTES)
                )
            except OSError:
                data = None
            seen[path] = data
        return seen[path]

    return load


def index_loader(repo_root: Path) -> Callable[[str], bytes | None]:
    """Read repo-relative paths from the git INDEX, capped and memoised.

    The index — not the worktree — is what ``git ls-files`` enumerates and what a
    commit will contain, so it is the only source an exemption for a tracked or
    staged path may be decided from. Memoised because resolving one fixture reads
    up to four blobs (the record, the generator, the manifest, the data file) and
    each is two subprocesses.
    """
    seen: dict[str, bytes | None] = {}

    def load(path: str) -> bytes | None:
        if path not in seen:
            data: bytes | None = None
            code, out = _git(repo_root, ["rev-parse", f":{path}"])
            if code == 0:
                sha = out.decode("ascii", errors="replace").strip()
                blob_code, blob = _git(repo_root, ["cat-file", "blob", sha])
                if blob_code == 0:
                    data = blob[:MAX_SCAN_BYTES]
            seen[path] = data
        return seen[path]

    return load


def _reference_fixture_declaration(
    path: str, *, load: Callable[[str], bytes | None]
) -> tuple[ReferenceFixtureVerdict, str]:
    """ADR 0012 conditions 1, 2 and 4 for one repo-relative path.

    Returns the verdict and the repo-relative path of the DATA file the entry
    stands for — itself, or the fixture an index sidecar declares it indexes.
    Condition 3 is **not** decided here: it is a property of bytes, and this
    function reads none of the fixture's own.
    """
    prefix = next((p for p in PUBLIC_REFERENCE_FIXTURE_PREFIXES if path.startswith(p)), None)
    if prefix is None:
        return ReferenceFixtureVerdict(False, "not under a public-reference fixture prefix"), path
    name = path[len(prefix) :]
    if "/" in name:
        return (
            ReferenceFixtureVerdict(
                False, "reference fixtures are declared per directory; this one is nested"
            ),
            path,
        )
    doc, reason = _load_fixture_provenance(prefix, load=load)
    if doc is None:
        return ReferenceFixtureVerdict(False, reason), path
    entry = doc.fixtures.get(name)
    if entry is None:
        return (
            ReferenceFixtureVerdict(False, f"not declared in {FIXTURE_PROVENANCE_FILENAME}"),
            path,
        )
    resource_name = entry.resource
    data_path = path
    if resource_name is None:
        indexed = entry.indexes or ""
        target = doc.fixtures.get(indexed)
        if target is None or target.resource is None:
            # Unreachable through _load_fixture_provenance, which rejects such a
            # document outright; kept so this function does not depend on that.
            return ReferenceFixtureVerdict(False, "index sidecar names no data fixture"), path
        resource_name = target.resource
        data_path = prefix + indexed
    resources = _manifest_resources(load(RESOURCE_MANIFEST_RELPATH))
    if resources is None:
        return (
            ReferenceFixtureVerdict(False, "the resource manifest is missing or unreadable"),
            data_path,
        )
    record = resources.get(resource_name)
    if record is None:
        return (
            ReferenceFixtureVerdict(
                False, "the declared resource is not registered in the manifest"
            ),
            data_path,
        )
    if record.get("status") != "fetched":
        return (
            ReferenceFixtureVerdict(False, "the declared resource is not recorded as fetched"),
            data_path,
        )
    digest = record.get("sha256")
    if not isinstance(digest, str) or _SHA256_HEX.match(digest) is None:
        return (
            ReferenceFixtureVerdict(False, "the declared resource carries no real sha256"),
            data_path,
        )
    return (
        ReferenceFixtureVerdict(
            True, "descent declared from a fetched, hash-pinned public release"
        ),
        data_path,
    )


def reference_fixture_provenance(
    path: str, *, load: Callable[[str], bytes | None]
) -> ReferenceFixtureVerdict:
    """ADR 0012 conditions 1, 2 and 4 for one repo-relative path — NOT condition 3.

    This is the half that can be decided without the fixture's own bytes. The
    other half is :func:`buffer_is_sites_only_vcf`, applied by the caller to
    whichever buffer it is actually handling.
    """
    return _reference_fixture_declaration(path, load=load)[0]


def reference_fixture_exemption(
    path: str, *, load: Callable[[str], bytes | None]
) -> ReferenceFixtureVerdict:
    """All four ADR 0012 conditions, for a caller that does not hold the bytes.

    ``load`` returns the bytes of a repo-relative path **from the source under
    audit** — the git index for the tracked and staged checks, the worktree for a
    filesystem walk. It is a parameter rather than an ``open()`` because the defect
    this replaces probed the worktree while the check reported on the staged blob:
    two different files, one verdict, and the genotypes were in the one nobody
    looked at.

    An index sidecar is resolved to the fixture it DECLARES it indexes and judged
    on that file's bytes from the same source. A sidecar whose companion is absent
    from that source earns nothing.
    """
    verdict, data_path = _reference_fixture_declaration(path, load=load)
    if not verdict.exempt:
        return verdict
    data = load(data_path)
    if data is None:
        return ReferenceFixtureVerdict(
            False, "the data file it stands for is absent or unreadable in the source under audit"
        )
    if not buffer_is_sites_only_vcf(data):
        return ReferenceFixtureVerdict(
            False, "the bytes are not a single sites-only VCF (condition 3)"
        )
    return ReferenceFixtureVerdict(
        True, "sites-only VCF with declared descent from a hash-pinned public release"
    )


#: The only places a ``!`` negation may re-admit a FILE.
NEGATION_ALLOWED_PREFIXES: Final[tuple[str, ...]] = (
    "tests/fixtures/synthetic/",
    "knowledge/",
    "config/",
    "templates/",
    "tests/golden/",
    *PUBLIC_REFERENCE_FIXTURE_PREFIXES,
)

#: Directories never walked: build caches and the object store itself.
SKIP_DIR_NAMES: Final[frozenset[str]] = frozenset(
    {
        ".git",
        ".venv",
        "venv",
        "env",
        "__pycache__",
        ".pytest_cache",
        ".ruff_cache",
        ".mypy_cache",
        ".hypothesis",
        "node_modules",
        ".tox",
        ".nox",
        "htmlcov",
    }
)

#: Synthetic probe paths for ``gitignore_effectiveness``. These files are never
#: created: ``git check-ignore --no-index`` evaluates the patterns against a
#: hypothetical path, so the probe cannot itself write a sensitive-looking file
#: into the repo it is auditing.
GITIGNORE_PROBES: Final[tuple[str, ...]] = (
    "workspace/proband.vcf",
    "workspace/runs/r1/calls.vcf.gz",
    "runs/case01/out.bam",
    "runs/case01/out.cram",
    "a/b/sample_R1.fastq.gz",
    "a/b/sample_R2.fastq",
    "evidence.duckdb",
    "variants.parquet",
    "data/cohort.sqlite",
    "notes/family.ped",
    "raw/aligned.sam",
    "case.phenopacket.json",
    "secrets.yaml",
    "pipeline.log",
    "analysis.ipynb",
    ".env",
)

_MAX_FINDINGS_PER_RULE: Final[int] = 3
_MAX_MATCHES_PER_RULE: Final[int] = 500
_MAX_FINDINGS_PER_CHECK: Final[int] = 200

_EMPTY_TREE: Final[str] = "4b825dc642cb6eb9a060e54bf8d69288fbee4904"
_GIT: Final[str] = shutil.which("git") or "git"

#: Stand-in for the absolute workspace path in every emitted finding. The location
#: and naming of the patient directory is itself disclosive (``.../smith_family/``),
#: and this report is an artifact that may be committed.
WORKSPACE_LABEL: Final[str] = "$MVA_WORKSPACE"

#: Path labels that are structural, not filesystem paths, and are emitted verbatim.
_LITERAL_PATH_LABELS: Final[frozenset[str]] = frozenset({".", "<logging>", WORKSPACE_LABEL})

#: A digit run long enough to be an identifier rather than a version or a year.
#: Sequencing filenames routinely ARE the identifier — an MRN, an NHS number, a
#: hospital accession — and until now the report printed them verbatim.
_PATH_IDENTIFIER_RUN: Final = re.compile(r"(?<!\d)\d{7,}(?!\d)")

#: File-name suffixes are shapes, not values, so they survive redaction. Anything
#: that does not look like a plain extension does not.
_SAFE_SUFFIX: Final = re.compile(r"\.[A-Za-z0-9]{1,8}\Z")


class PathRedactor:
    """Redacts paths for one rendered report, labelling by first appearance.

    GP-41 says the audit emits "paths and counts, never matched content". That is
    sound right up to the point where the path IS the content. Sequencing
    filenames carry MRNs, NHS numbers, accessions and surnames as a matter of
    routine — ``data/NHS9999999999_Faketon.vcf`` is an ordinary filename, not a
    contrived one — and this report is ``cat``-ed by the justfile, written into
    every run as a ``DERIVED_SAFE`` artifact, and read into a model context.
    Emitting ``finding.path`` verbatim defeated every other control in the package.

    A component is rewritten when it carries anything rule-detectable or an
    identifier-length digit run. The whole component goes, not just the matched
    span, because the surviving remainder of ``NHS9999999999_Faketon`` is the half
    that names a person. What is kept is the directory chain (structure) and the
    extension chain (shape).

    **Why a counter and not a hash.** The label used to be
    :func:`~mva.privacy.patterns.correlation_id` — an HMAC under a
    ``secrets.token_bytes(16)`` module salt regenerated at every interpreter start.
    That made ``privacy/privacy_audit.md`` a random function of a secret, and that
    file is a registered artifact inside the GP-30 byte-identity claim and is not
    in ``verify_determinism``'s skip set. Two processes rendering the same
    ``AuditReport`` produced different bytes. It stayed green only because a clean
    tree has no redacted paths to label, and ``just demo-determinism`` runs both
    passes in one interpreter, where the salt is constant — so the check could
    never observe it.

    A counter fixes determinism without weakening anything, because it is
    *strictly less* invertible than the HMAC was: ``001`` is a function of position
    in the report, not of the path at all. There is no key to leak and no
    keyspace to search. It keeps the one property the tag existed for — two
    findings about the same file share a label, two different files do not — and
    that property is what a reader needs to act ("the same file trips three
    rules").

    The random salt stays in :func:`~mva.privacy.patterns.correlation_id` for
    in-memory distinct-counting (``scan_bytes`` counting distinct HPO terms), which
    is the use it was designed for and where nothing is ever rendered.

    One instance per rendered report: labels are meaningful only within the
    document that assigned them, and saying so in the type is cheaper than saying
    so in a comment.
    """

    __slots__ = ("_labels",)

    def __init__(self) -> None:
        self._labels: dict[str, str] = {}

    def redact(self, path: str) -> str:
        """Rewrite every disclosive component of ``path``."""
        if path in _LITERAL_PATH_LABELS or not path:
            return path
        return "/".join(self._redact_component(part) for part in path.split("/"))

    def _redact_component(self, component: str) -> str:
        if not component or component in {".", ".."}:
            return component
        cleaned = _PATH_IDENTIFIER_RUN.sub("0", redact_text(component))
        if cleaned == component:
            return component
        label = self._labels.get(component)
        if label is None:
            # Keyed by the raw component so the SAME file gets the SAME label
            # throughout one report. The dict is report-local and dies with it;
            # it is never rendered, only its values are.
            label = f"{len(self._labels) + 1:03d}"
            self._labels[component] = label
        return f"<REDACTED:path:{label}>{''.join(_safe_suffixes(component))}"


def redact_path(path: str) -> str:
    """Redact a single path with a throwaway :class:`PathRedactor`.

    For callers that need one path and no cross-path distinctness: the
    ``did this path change?`` test in :func:`_path_findings`, and the ``privacy
    classify`` style of one-shot use. A report renders through ONE
    :class:`PathRedactor` instead, or every redacted component in it would be
    labelled ``001``.
    """
    return PathRedactor().redact(path)


def _safe_suffixes(component: str) -> list[str]:
    """The trailing extension chain, at most two deep, and only if it looks like one."""
    suffixes = [s for s in Path(component).suffixes[-2:] if _SAFE_SUFFIX.fullmatch(s)]
    return suffixes


# ---------------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Finding:
    """One observation. ``detail`` is prose about the RULE, never about the match."""

    check: str
    path: str
    line: int | None
    rule_id: str | None
    span_len: int | None
    severity: Severity
    detail: str


@dataclass(frozen=True, slots=True)
class CheckResult:
    name: str
    passed: bool
    severity: str
    findings: tuple[Finding, ...]
    summary: str


@dataclass(frozen=True, slots=True)
class AuditReport:
    results: tuple[CheckResult, ...]
    passed: bool
    failed_checks: tuple[str, ...]

    @property
    def findings(self) -> tuple[Finding, ...]:
        return tuple(f for result in self.results for f in result.findings)

    def to_markdown(self) -> str:
        """Render for a human or an agent. Contains no matched content by construction.

        Byte-identical across processes for the same report: the only non-obvious
        ingredient is the redacted-path label, and that is a first-appearance
        counter over the findings in their existing deterministic order (GP-30).
        See :class:`PathRedactor`.
        """
        redactor = PathRedactor()
        status = "PASS" if self.passed else "FAIL"
        lines: list[str] = [
            "# Privacy audit",
            "",
            f"**Result: {status}** — {len(self.results)} checks, "
            f"{len(self.failed_checks)} failed, {len(self.findings)} findings.",
            "",
            "> This report intentionally contains no matched content: only paths, "
            "line numbers, span lengths and rule IDs (GP-41).",
            "",
            "| Check | Status | Findings | Summary |",
            "| --- | --- | --- | --- |",
        ]
        for result in self.results:
            mark = "pass" if result.passed else "**FAIL**"
            summary = result.summary.replace("|", "\\|")
            lines.append(f"| `{result.name}` | {mark} | {len(result.findings)} | {summary} |")

        for result in self.results:
            if not result.findings:
                continue
            lines += ["", f"## {result.name}", ""]
            for finding in result.findings:
                location = redactor.redact(finding.path)
                if finding.line is not None:
                    location = f"{location}:{finding.line}"
                bits = [f"- `{finding.severity}` `{location}`"]
                if finding.rule_id is not None:
                    bits.append(f"rule=`{finding.rule_id}`")
                if finding.span_len is not None:
                    bits.append(f"span_len={finding.span_len}")
                bits.append(f"— {finding.detail}")
                lines.append(" ".join(bits))

        if not self.passed:
            lines += [
                "",
                "## Remediation",
                "",
                "Failing checks are listed above. Do NOT weaken a check to make it "
                "pass: fix the repository, or record a decision in `docs/decisions/` "
                "explaining why a specific finding is acceptable.",
            ]
        return "\n".join(lines) + "\n"

    def to_dict(self) -> dict[str, object]:
        """As :meth:`to_markdown`, and equally process-stable. One redactor for the
        whole document, so labels are consistent within it."""
        redactor = PathRedactor()
        return {
            "passed": self.passed,
            "failed_checks": list(self.failed_checks),
            "results": [
                {
                    "name": result.name,
                    "passed": result.passed,
                    "severity": result.severity,
                    "summary": result.summary,
                    "findings": [
                        {
                            "check": f.check,
                            "path": redactor.redact(f.path),
                            "line": f.line,
                            "rule_id": f.rule_id,
                            "span_len": f.span_len,
                            "severity": f.severity,
                            "detail": f.detail,
                        }
                        for f in result.findings
                    ],
                }
                for result in self.results
            ],
        }


def _result(name: str, findings: Sequence[Finding], summary: str, *, strict: bool) -> CheckResult:
    """Derive pass/fail. Under ``strict`` a warning is also disqualifying."""
    trimmed = tuple(findings[:_MAX_FINDINGS_PER_CHECK])
    if len(findings) > _MAX_FINDINGS_PER_CHECK:
        summary = f"{summary} (showing {_MAX_FINDINGS_PER_CHECK} of {len(findings)} findings)"
    has_fail = any(f.severity == "fail" for f in findings)
    has_warn = any(f.severity == "warn" for f in findings)
    passed = not has_fail and not (strict and has_warn)
    severity = "fail" if has_fail else ("warn" if has_warn else "pass")
    return CheckResult(
        name=name, passed=passed, severity=severity, findings=trimmed, summary=summary
    )


# ---------------------------------------------------------------------------
# Git plumbing
# ---------------------------------------------------------------------------


def _git(repo_root: Path, args: Sequence[str], *, stdin: bytes | None = None) -> tuple[int, bytes]:
    """Run a git plumbing command. Never raises on non-zero; callers decide."""
    try:
        proc = subprocess.run(  # noqa: S603
            [_GIT, "-C", str(repo_root), *args],
            capture_output=True,
            input=stdin,
            check=False,
        )
    except OSError:
        return 127, b""
    return proc.returncode, proc.stdout


def _split_nul(data: bytes) -> list[str]:
    return [chunk.decode("utf-8", errors="replace") for chunk in data.split(b"\0") if chunk]


def _git_available(repo_root: Path) -> bool:
    code, _ = _git(repo_root, ["rev-parse", "--git-dir"])
    return code == 0


def _exempt(path: str, prefixes: Sequence[str]) -> bool:
    return any(path.startswith(prefix) for prefix in prefixes)


def _unavailable(name: str, reason: str, *, strict: bool) -> CheckResult:
    finding = Finding(
        check=name,
        path=".",
        line=None,
        rule_id=None,
        span_len=None,
        severity="warn",
        detail=reason,
    )
    return _result(name, [finding], f"Skipped: {reason}", strict=strict)


# ---------------------------------------------------------------------------
# Content scanning
# ---------------------------------------------------------------------------


def _line_offsets(data: bytes) -> list[int]:
    offsets: list[int] = []
    start = 0
    while True:
        index = data.find(b"\n", start)
        if index == -1:
            return offsets
        offsets.append(index)
        start = index + 1


def _line_of(offsets: Sequence[int], position: int) -> int:
    return bisect_right(offsets, position - 1) + 1


#: Sample columns of a VCF ``#CHROM`` line, and ``@RG SM:`` values, must ALL carry
#: a synthetic prefix for a file to count as declared synthetic. A marker comment
#: is a claim; the sample names are the part a real export would have got wrong.
_VCF_CHROM_SAMPLES: Final = re.compile(
    rb"^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO(?:\tFORMAT((?:\t[^\t\r\n]+)+))?",
    re.MULTILINE,
)
_RG_SAMPLE_VALUE: Final = re.compile(rb"^@RG[ \t](?:[!-~]+[ \t])*SM:([!-~]+)", re.MULTILINE)

#: Sample-name prefixes that mark a fabricated subject.
_SYNTHETIC_SAMPLE_PREFIXES: Final[tuple[bytes, ...]] = (b"SYNTH_", b"SYNTH-", b"SYNTHETIC")


def synthetic_marker_present(data: bytes) -> bool:
    """Does the HEAD of these bytes carry a synthetic declaration?"""
    head = data[:SYNTHETIC_MARKER_BYTES]
    return any(marker in head for marker in SYNTHETIC_MARKERS)


def unmarked_sample_names(data: bytes) -> int:
    """How many VCF sample columns / ``@RG SM:`` values are NOT synthetic-prefixed.

    A declaration in a comment is unverifiable prose. A sample name is not: it is
    the field a real sequencing centre fills in with an accession or a patient's
    initials, and it is the field a genuinely fabricated fixture controls. So the
    marker is only believed when every subject named in the file is named as a
    fabricated one.
    """
    unmarked = 0
    for match in _VCF_CHROM_SAMPLES.finditer(data):
        columns = match.group(1)
        if not columns:
            continue
        unmarked += sum(
            1
            for sample in columns.split(b"\t")
            if sample and not sample.startswith(_SYNTHETIC_SAMPLE_PREFIXES)
        )
    unmarked += sum(
        1
        for match in _RG_SAMPLE_VALUE.finditer(data)
        if not match.group(1).startswith(_SYNTHETIC_SAMPLE_PREFIXES)
    )
    return unmarked


def declares_synthetic(data: bytes) -> bool:
    """Whether a file credibly declares itself fabricated.

    Both halves are required, and each closes a different hole:

    * a marker in the first :data:`SYNTHETIC_MARKER_BYTES` bytes — a declaration
      is a header. Accepting ``SYNTH_`` anywhere in 8 MiB meant a real call set
      with one incidental ``SYNTH_`` token, at any depth, counted as declared;
    * every named subject carries a synthetic sample prefix — which is what a real
      exported VCF or BAM header cannot satisfy without someone deliberately
      rewriting its sample columns.
    """
    return synthetic_marker_present(data) and unmarked_sample_names(data) == 0


def _resolve_context(
    *,
    rule_id: str,
    base_severity: Severity,
    description: str,
    matched_ids: set[str],
    distinct_hpo: int,
    allowlisted: bool,
    synthetic_declared: bool,
    hpo_is_constant: bool,
    sites_only_reference: bool,
) -> tuple[Severity, str]:
    """Refine a rule's isolated severity using whole-file context.

    Kept separate from :func:`scan_bytes` because this is the part reviewers should
    argue about: it is where a `warn` becomes a `fail`, and where a documented
    false-positive class is held down. Every branch states its reason in the
    returned note, so the report explains its own severities.

    The allowlist is deliberately NOT a blanket downgrade. It softens
    :data:`PATH_DOWNGRADABLE_RULES` on the strength of the path, and
    :data:`DECLARED_SYNTHETIC_DOWNGRADABLE_RULES` only when the file also passes
    :func:`declares_synthetic`. Everything else — every ``magic:*`` container hit
    included — keeps full strength wherever it is found, because a rule with
    near-zero false-positive risk firing inside an allowlisted directory is not
    noise, it is the exact event the audit exists to catch.
    """
    severity: Severity = base_severity
    note = description
    if rule_id == "vcf_data_line":
        promoted = bool({"vcf_header", "genotype_field"} & matched_ids)
        severity = "fail" if promoted else "warn"
        note = (
            f"{note} Promoted to fail: the file also matches a VCF header or a genotype field."
            if promoted
            else f"{note} Held at warn: no VCF header or genotype field in this file."
        )
    elif rule_id == "hpo_term":
        over_threshold = distinct_hpo >= HPO_DISTINCT_FAIL_THRESHOLD
        severity = "fail" if over_threshold and not hpo_is_constant else "warn"
        note = f"{note} {distinct_hpo} distinct term(s) in this file."
        if over_threshold and hpo_is_constant:
            note = (
                f"{note} Held at warn: a reviewed source/reference path AND the bytes "
                "do not have per-subject phenotype-profile shape."
            )

    if severity != "fail":
        return severity, note
    # Checked BEFORE the allowlist gate, and deliberately: this condition does not
    # rest on the path allowlist at all. By the time it is true, all four ADR 0012
    # conditions hold -- an allowlisted prefix, a declared descent from a manifest
    # resource recorded as fetched with a real sha256, a generator beside it, and
    # a verification against THESE bytes that they form one sites-only VCF. That
    # is a strictly stronger claim than allowlist membership, so it stands alone.
    if sites_only_reference and rule_id in SITES_ONLY_DOWNGRADABLE_RULES:
        return "warn", (
            f"{note} Downgraded: a declared public reference slice, and this exact "
            "buffer is a single sites-only VCF -- one #CHROM line of exactly 8 "
            "columns, no FORMAT, no wider row (ADR 0012). It holds real public "
            "records but no individual."
        )
    if not allowlisted:
        return severity, note
    if rule_id in PATH_DOWNGRADABLE_RULES:
        return "warn", f"{note} Downgraded: path is on the audited public/synthetic allowlist."
    if synthetic_declared and rule_id in DECLARED_SYNTHETIC_DOWNGRADABLE_RULES:
        return "warn", (
            f"{note} Downgraded: allowlisted path AND the file declares itself synthetic "
            "(marker in the head, all sample names synthetic-prefixed)."
        )
    return severity, (
        f"{note} NOT downgraded despite the allowlisted path: this rule is only "
        "softened for a file that declares itself synthetic, and this one does not."
    )


def scan_bytes(
    data: bytes,
    *,
    check: str,
    path_label: str,
    allowlisted: bool = False,
    hpo_is_constant: bool = False,
    reference_fixture: bool = False,
    _gzip_depth: int = 0,
) -> list[Finding]:
    """Apply the rule battery to one buffer and return content-free findings.

    Two rules are resolved here rather than in :mod:`mva.privacy.patterns`, because
    only this function has whole-file context:

    * ``vcf_data_line`` is promoted from ``warn`` to ``fail`` when the same file
      also matches ``vcf_header`` or ``genotype_field``. A tab-separated
      coordinate table is a public resource; the same table next to a VCF header
      is a patient's call set.
    * ``hpo_term`` fails only at :data:`HPO_DISTINCT_FAIL_THRESHOLD` DISTINCT terms
      in a non-allowlisted file. Distinctness is counted over
      :func:`~mva.privacy.patterns.correlation_id` values, so the terms themselves
      are never accumulated in a data structure that could later be printed.
      ``hpo_is_constant`` — a claim about the path — is additionally verified here
      against the buffer with :func:`has_phenotype_profile_shape`.
    """
    offsets = _line_offsets(data)
    spans_by_rule: dict[str, list[tuple[int, int]]] = {}
    distinct_hpo = 0
    synthetic_declared = allowlisted and declares_synthetic(data)
    # `hpo_is_constant` arrives as a claim about the PATH. This is its verification
    # against the BYTES, and both halves must hold: a reviewed-source prefix does
    # not license a file that is shaped like one patient's phenotype profile. The
    # two are ANDed here rather than at the call sites because only this function
    # has the buffer -- the same reason `vcf_data_line` and `hpo_term` are resolved
    # here at all. Deliberately NOT routed through `synthetic_declared`: per ADR
    # 0012 the exemptions do not share a code path, so widening one cannot widen
    # the other by accident.
    hpo_constant = hpo_is_constant and not has_phenotype_profile_shape(data)
    # ADR 0012 condition 3, decided against THIS buffer and no other. The caller
    # supplies only the provenance half (`reference_fixture`, a claim about the
    # path); the sites-only half is never taken from a path reopened beside the
    # scan, because a staged blob and its worktree file are routinely different
    # files -- which is exactly how a genotyped fixture used to be staged and then
    # cleaned on disk, leaving the audit probing one file and reporting on another.
    sites_only_reference = reference_fixture and buffer_is_sites_only_vcf(data)

    for rule in RULES:
        spans: list[tuple[int, int]] = []
        seen: set[str] = set()
        for match in rule.pattern.finditer(data):
            spans.append(match.span())
            if rule.rule_id == "hpo_term":
                # The matched bytes exist only as an argument to the salted HMAC;
                # they are never bound to a name that outlives this expression.
                seen.add(correlation_id(match.group(0)))
            if len(spans) >= _MAX_MATCHES_PER_RULE:
                break
        if spans:
            spans_by_rule[rule.rule_id] = spans
            if rule.rule_id == "hpo_term":
                distinct_hpo = len(seen)

    matched_ids = set(spans_by_rule)
    findings: list[Finding] = []

    for rule in RULES:
        spans = spans_by_rule.get(rule.rule_id, [])
        if not spans:
            continue
        severity, note = _resolve_context(
            rule_id=rule.rule_id,
            base_severity=rule.severity,
            description=rule.description,
            matched_ids=matched_ids,
            distinct_hpo=distinct_hpo,
            allowlisted=allowlisted,
            synthetic_declared=synthetic_declared,
            hpo_is_constant=hpo_constant,
            sites_only_reference=sites_only_reference,
        )
        total = len(spans)
        for start, end in spans[:_MAX_FINDINGS_PER_RULE]:
            findings.append(
                Finding(
                    check=check,
                    path=path_label,
                    line=_line_of(offsets, start),
                    rule_id=rule.rule_id,
                    span_len=end - start,
                    severity=severity,
                    detail=f"{note} ({total} occurrence(s) in this file.)",
                )
            )

    kind = sniff_binary(data[:HEAD_BYTES])
    if kind is not None:
        findings.append(
            Finding(
                check=check,
                path=path_label,
                line=None,
                rule_id=f"magic:{kind}",
                span_len=None,
                severity=BINARY_SEVERITY[kind],
                # A magic-byte hit is NEVER downgraded by the allowlist. There is
                # no false-positive story for "these four bytes are a BAM": an
                # aligned read set inside tests/fixtures/synthetic/ is a leak in
                # exactly the way it is anywhere else.
                detail=(
                    f"Container identified from magic bytes as {kind!r}. "
                    "Extension checks are defeated by a rename; this is not. "
                    "Container hits are never softened by the path allowlist."
                ),
            )
        )

    findings.extend(
        _scan_gzip_member(
            data,
            check=check,
            path_label=path_label,
            allowlisted=allowlisted,
            hpo_is_constant=hpo_is_constant,
            reference_fixture=reference_fixture,
            gzip_depth=_gzip_depth,
        )
    )
    return findings


#: How deep a gzip-in-gzip chain is followed before giving up. One level covers
#: every real case (``sample.vcf.gz``); the bound exists so a nested bomb cannot
#: recurse.
_MAX_GZIP_DEPTH: Final[int] = 2


def _scan_gzip_member(
    data: bytes,
    *,
    check: str,
    path_label: str,
    allowlisted: bool,
    hpo_is_constant: bool,
    reference_fixture: bool,
    gzip_depth: int,
) -> list[Finding]:
    """Inflate a gzip member and run the whole battery over the plaintext.

    Every rule matches plaintext, so ``sample.vcf.gz`` — the form variant callers
    actually write — matched nothing at all: not the header, not the ``#CHROM``
    line, not one genotype. BGZF at least tripped the container sniffer (at
    ``warn``, because every bgzipped public resource has that header); plain gzip
    tripped nothing whatsoever. Compression is not a privacy control, and treating
    it as an opaque blob made it one.

    Line numbers are reported in the DECOMPRESSED stream, which is the only place
    a line number means anything, and the detail says so.
    """
    if gzip_depth >= _MAX_GZIP_DEPTH:
        return []
    plain = gunzip_capped(data, MAX_SCAN_BYTES)
    if plain is None:
        return []
    return [
        Finding(
            check=finding.check,
            path=finding.path,
            line=finding.line,
            rule_id=finding.rule_id,
            span_len=finding.span_len,
            severity=finding.severity,
            detail=(
                f"{finding.detail} Found INSIDE a gzip member; the line number is "
                "an offset into the decompressed stream."
            ),
        )
        for finding in scan_bytes(
            plain,
            check=check,
            path_label=path_label,
            allowlisted=allowlisted,
            hpo_is_constant=hpo_is_constant,
            reference_fixture=reference_fixture,
            _gzip_depth=gzip_depth + 1,
        )
    ]


def scan_file(
    path: Path,
    *,
    check: str,
    path_label: str,
    allowlisted: bool,
    hpo_is_constant: bool = False,
    reference_fixture: bool = False,
) -> list[Finding]:
    """Scan one file on disk, capped, in bytes, never decoding for detection."""
    try:
        size = path.stat().st_size
    except OSError:
        return []
    try:
        data = read_capped(path, MAX_SCAN_BYTES)
    except OSError as exc:
        return [
            Finding(
                check=check,
                path=path_label,
                line=None,
                rule_id=None,
                span_len=None,
                severity="warn",
                detail=f"Unreadable ({type(exc).__name__}); not scanned.",
            )
        ]
    findings = scan_bytes(
        data,
        check=check,
        path_label=path_label,
        allowlisted=allowlisted,
        hpo_is_constant=hpo_is_constant,
        reference_fixture=reference_fixture,
    )
    if size > MAX_SCAN_BYTES:
        findings.append(
            Finding(
                check=check,
                path=path_label,
                line=None,
                rule_id=None,
                span_len=None,
                severity="warn",
                detail=(
                    f"File is {size} bytes; only the first {MAX_SCAN_BYTES} were scanned. "
                    "Magic-byte identification still applied to the head."
                ),
            )
        )
    return findings


# ---------------------------------------------------------------------------
# Workspace resolution shared by two checks
# ---------------------------------------------------------------------------


def _configured_workspace(workspace: Path | None) -> Path | None:
    raw: str | Path | None = workspace if workspace is not None else os.environ.get("MVA_WORKSPACE")
    if raw is None:
        return None
    try:
        return Path(raw).expanduser().resolve()
    except OSError:
        return None


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------


def _reference_exemption_note(verdict: ReferenceFixtureVerdict, path: str) -> str:
    """Why a file under a reference-fixture prefix did NOT earn the exemption.

    Attached to the extension findings so a refusal explains itself: without it a
    typo in ``fixture_provenance.yaml`` and a genotyped fixture produce the same
    opaque failure, and the first one gets "fixed" by deleting the check.
    """
    if not _exempt(path, PUBLIC_REFERENCE_FIXTURE_PREFIXES):
        return ""
    return (
        " It is under a public-reference fixture prefix but is NOT exempt "
        f"(ADR 0012): {verdict.reason}."
    )


def check_git_tracked_sensitive(repo_root: Path, *, strict: bool = False) -> CheckResult:
    """Tracked files with a sensitive extension.

    This is the check that matters most, and the one people assume ``.gitignore``
    already covers. It does not: ignore rules are consulted only for *untracked*
    files. Once a path is in the index, git tracks it forever regardless of any
    pattern, and `git rm --cached` leaves it recoverable in history.

    The ADR 0012 exemption is decided from the INDEX blob, never from the file on
    disk. ``git ls-files`` reports what the index holds, and the index is what a
    commit will carry; clearing a path on the strength of a worktree copy that
    was cleaned after staging is the same defect twice.
    """
    name = "git_tracked_sensitive"
    if not _git_available(repo_root):
        return _unavailable(name, "not a git repository (or git unavailable)", strict=strict)
    code, out = _git(repo_root, ["ls-files", "-z"])
    if code != 0:
        return _unavailable(name, "git ls-files failed", strict=strict)
    paths = _split_nul(out)
    load = index_loader(repo_root)
    findings: list[Finding] = []
    for path in sorted(paths):
        if not is_sensitive_extension(Path(path)) or _exempt(path, TRACKED_EXEMPT_PREFIXES):
            continue
        verdict = reference_fixture_exemption(path, load=load)
        if verdict.exempt:
            continue
        findings.append(
            Finding(
                check=name,
                path=path,
                line=None,
                rule_id=None,
                span_len=None,
                severity="fail",
                detail=(
                    "Tracked file has a patient-data extension. .gitignore does not "
                    "apply to tracked paths. Remove it from the index and purge history."
                    + _reference_exemption_note(verdict, path)
                ),
            )
        )
    return _result(
        name,
        findings,
        f"{len(paths)} tracked paths inspected; {len(findings)} with sensitive extensions.",
        strict=strict,
    )


def check_git_staged_sensitive(repo_root: Path, *, strict: bool = False) -> CheckResult:
    """Sensitive extensions or content in the STAGED blobs.

    Deliberately reads ``git cat-file blob`` rather than the worktree. The staged
    content is what a commit will contain; a file can be staged dirty and then
    cleaned on disk, and only the index still holds the payload. This is the
    pre-commit hook's check.
    """
    name = "git_staged_sensitive"
    if not _git_available(repo_root):
        return _unavailable(name, "not a git repository (or git unavailable)", strict=strict)

    args = ["diff", "--cached", "--name-only", "-z", "--diff-filter=ACMR"]
    code, out = _git(repo_root, args)
    if code != 0:
        # No HEAD yet: diff against the empty tree so the first commit is audited too.
        code, out = _git(repo_root, [*args, _EMPTY_TREE])
        if code != 0:
            return _unavailable(name, "git diff --cached failed", strict=strict)

    paths = sorted(_split_nul(out))
    findings: list[Finding] = []
    load = index_loader(repo_root)
    for path in paths:
        # The pre-commit hook runs this subset, so the path rule has to be here
        # too: a filename that IS the identifier is committable without the file
        # containing a single byte of genomic content.
        findings.extend(_path_findings(path, check=name))
        allowlisted = _exempt(path, CONTENT_ALLOWLIST_PREFIXES)
        verdict = reference_fixture_exemption(path, load=load)
        if (
            is_sensitive_extension(Path(path))
            and not _exempt(path, TRACKED_EXEMPT_PREFIXES)
            and not verdict.exempt
        ):
            findings.append(
                Finding(
                    check=name,
                    path=path,
                    line=None,
                    rule_id=None,
                    span_len=None,
                    severity="fail",
                    detail=(
                        "Staged file has a patient-data extension. Unstage it."
                        + _reference_exemption_note(verdict, path)
                    ),
                )
            )
        blob = load(path)
        if blob is None:
            continue
        # `reference_fixture` carries ONLY the provenance half of ADR 0012.
        # scan_bytes decides the sites-only half against `blob` itself, which is
        # the buffer about to be committed -- the worktree file of the same name
        # may have been cleaned since it was staged, and used to be what decided
        # this.
        findings.extend(
            scan_bytes(
                blob,
                check=name,
                path_label=path,
                allowlisted=allowlisted,
                hpo_is_constant=_exempt(path, HPO_CODE_PREFIXES),
                reference_fixture=reference_fixture_provenance(path, load=load).exempt,
            )
        )
    return _result(name, findings, f"{len(paths)} staged blob(s) inspected.", strict=strict)


def _parse_check_ignore(line: str) -> tuple[str, str] | None:
    """Parse one ``git check-ignore -v`` line into (source_location, pattern)."""
    if "\t" not in line:
        return None
    left, _, _ = line.partition("\t")
    parts = left.split(":", 2)
    if len(parts) != 3:
        return None
    source, lineno, pattern = parts
    return f"{source}:{lineno}", pattern


def check_gitignore_effectiveness(repo_root: Path, *, strict: bool = False) -> CheckResult:
    """Does ``.gitignore`` actually ignore the paths it is supposed to?

    The subtlety that makes this check worth writing: ``git check-ignore`` exits 0
    when *a rule matched*, and a NEGATION (``!pattern``) is a rule. So exit 0 alone
    means "some rule had an opinion", not "the file is ignored". A pass therefore
    requires exit 0 **and** a reported pattern that does not begin with ``!``.

    ``--no-index`` evaluates the patterns against a hypothetical path, so no probe
    file is ever written — which matters, because writing ``workspace/proband.vcf``
    to test whether it would be ignored is exactly the accident being guarded
    against.
    """
    name = "gitignore_effectiveness"
    if not _git_available(repo_root):
        return _unavailable(name, "not a git repository (or git unavailable)", strict=strict)

    findings: list[Finding] = []
    for probe in GITIGNORE_PROBES:
        code, out = _git(repo_root, ["check-ignore", "-v", "--no-index", "--", probe])
        text = out.decode("utf-8", errors="replace").strip()
        if code != 0 or not text:
            findings.append(
                Finding(
                    check=name,
                    path=probe,
                    line=None,
                    rule_id=None,
                    span_len=None,
                    severity="fail",
                    detail=(
                        "No .gitignore rule matches this probe path; a real file here "
                        "would be committable."
                    ),
                )
            )
            continue
        parsed = _parse_check_ignore(text.splitlines()[0])
        if parsed is None:
            findings.append(
                Finding(
                    check=name,
                    path=probe,
                    line=None,
                    rule_id=None,
                    span_len=None,
                    severity="warn",
                    detail="Unparseable `git check-ignore -v` output.",
                )
            )
            continue
        source, pattern = parsed
        if pattern.startswith("!"):
            findings.append(
                Finding(
                    check=name,
                    path=probe,
                    line=None,
                    rule_id=None,
                    span_len=None,
                    severity="fail",
                    detail=(
                        f"Matched by a NEGATION rule at {source} — the path is "
                        "RE-INCLUDED, not ignored. Exit code 0 from check-ignore does "
                        "not mean ignored."
                    ),
                )
            )
    return _result(
        name,
        findings,
        f"{len(GITIGNORE_PROBES)} synthetic probe paths evaluated with --no-index.",
        strict=strict,
    )


def check_gitignore_negation_safety(repo_root: Path, *, strict: bool = False) -> CheckResult:
    """No ``!`` rule may re-admit a file outside the audited public directories.

    A negation is the only way to punch a hole in a deny-by-default ignore file,
    so the blast radius of each one is reviewed here rather than trusted.

    One carve-out, stated explicitly because it is a real weakening: a *pure
    directory* re-inclusion (``!dir/`` or ``!dir/**/``) that is an ancestor of an
    allowed prefix is permitted. Git will not descend into an excluded directory,
    so those lines are structurally required for the narrow file negations beneath
    them to have any effect — and a directory negation cannot by itself re-admit a
    file, because the file patterns still apply inside.
    """
    name = "gitignore_negation_safety"
    gitignore = repo_root / ".gitignore"
    if not gitignore.is_file():
        return _unavailable(name, "no .gitignore at the repository root", strict=strict)

    findings: list[Finding] = []
    text = decode_lossy(read_capped(gitignore, MAX_SCAN_BYTES))
    negations = 0
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line.startswith("!"):
            continue
        negations += 1
        target = line[1:].strip().lstrip("/")
        if _exempt(target, NEGATION_ALLOWED_PREFIXES):
            continue
        directory = target.removesuffix("**/").removesuffix("**").rstrip("/")
        is_directory_rule = target.endswith("/")
        if is_directory_rule and any(
            allowed.startswith(f"{directory}/") for allowed in NEGATION_ALLOWED_PREFIXES
        ):
            continue
        findings.append(
            Finding(
                check=name,
                path=".gitignore",
                line=lineno,
                rule_id=None,
                span_len=None,
                severity="fail",
                detail=(
                    f"Negation `{line}` re-admits a path outside the audited set "
                    f"({', '.join(NEGATION_ALLOWED_PREFIXES)}). Every negation is a "
                    "hole in a deny-by-default policy and needs a decision record."
                ),
            )
        )
    return _result(name, findings, f"{negations} negation rule(s) reviewed.", strict=strict)


def check_workspace_containment(
    repo_root: Path, *, workspace: Path | None = None, strict: bool = False
) -> CheckResult:
    """GP-40: the patient workspace must resolve outside the repository.

    Containment is decided by :func:`mva.config.path_is_within`, which compares
    filesystem identity rather than resolved strings. A workspace that is a
    symlink into the repo tree passes a naive string comparison and fails this
    one — and so, now, does a workspace spelled in a different CASE from the
    repository, which on APFS is the same directory and which ``Path.resolve()``
    reports as unrelated.

    The absolute path is never echoed — the workspace directory name is itself
    disclosive — so findings report the relationship, not the location.
    """
    name = "workspace_containment"
    resolved = _configured_workspace(workspace)
    if resolved is None:
        return _unavailable(
            name,
            "no workspace configured (pass workspace= or set MVA_WORKSPACE)",
            strict=strict,
        )
    repo = repo_root.resolve()
    findings: list[Finding] = []
    if path_is_within(resolved, repo):
        findings.append(
            Finding(
                check=name,
                path=WORKSPACE_LABEL,
                line=None,
                rule_id=None,
                span_len=None,
                severity="fail",
                detail=(
                    "Workspace resolves INSIDE the repository (symlinks followed). "
                    "Move it out; data written there is one `git add -A` from the "
                    "index and stays recoverable from history afterwards."
                ),
            )
        )
    if not resolved.exists():
        findings.append(
            Finding(
                check=name,
                path=WORKSPACE_LABEL,
                line=None,
                rule_id=None,
                span_len=None,
                severity="warn",
                detail="Configured workspace does not exist yet.",
            )
        )
    return _result(name, findings, "Workspace containment evaluated.", strict=strict)


def check_symlink_escape(
    repo_root: Path, *, workspace: Path | None = None, strict: bool = False
) -> CheckResult:
    """Symlinks in the repo that reach into the workspace or at sensitive data.

    Walked with ``followlinks=False``: following them would both loop and,
    ironically, read the patient data this check exists to keep at arm's length.
    A link is a leak because git stores the *link*, but every tool that reads the
    repo — including an agent — reads through it.
    """
    name = "symlink_escape"
    repo = repo_root.resolve()
    resolved_workspace = _configured_workspace(workspace)
    findings: list[Finding] = []

    for dirpath, dirnames, filenames in os.walk(repo, followlinks=False):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIR_NAMES]
        base = Path(dirpath)
        for entry in sorted(dirnames) + sorted(filenames):
            link = base / entry
            if not link.is_symlink():
                continue
            label = link.relative_to(repo).as_posix()
            target = Path(os.path.realpath(link))
            if resolved_workspace is not None and (
                target == resolved_workspace or target.is_relative_to(resolved_workspace)
            ):
                findings.append(
                    Finding(
                        check=name,
                        path=label,
                        line=None,
                        rule_id=None,
                        span_len=None,
                        severity="fail",
                        detail=(
                            f"Symlink resolves into {WORKSPACE_LABEL}. The repo now "
                            "reads patient data through a tracked path."
                        ),
                    )
                )
            elif is_sensitive_extension(target):
                findings.append(
                    Finding(
                        check=name,
                        path=label,
                        line=None,
                        rule_id=None,
                        span_len=None,
                        severity="fail",
                        detail="Symlink resolves to a path with a patient-data extension.",
                    )
                )
            elif not target.is_relative_to(repo):
                findings.append(
                    Finding(
                        check=name,
                        path=label,
                        line=None,
                        rule_id=None,
                        span_len=None,
                        severity="warn",
                        detail="Symlink resolves outside the repository.",
                    )
                )
    return _result(name, findings, "Repository walked for symlinks.", strict=strict)


def _scannable_paths(repo_root: Path) -> list[str]:
    """Tracked plus untracked-but-not-ignored paths, deduplicated and sorted.

    Ignored files are excluded on purpose: an ignored patient VCF sitting in a
    developer's worktree is expected and is not this check's business. The checks
    that DO care about it are ``workspace_containment`` and ``gitignore_effectiveness``.
    """
    seen: set[str] = set()
    for args in (["ls-files", "-z"], ["ls-files", "-z", "--others", "--exclude-standard"]):
        code, out = _git(repo_root, args)
        if code == 0:
            seen.update(_split_nul(out))
    return sorted(seen)


def _path_findings(path: str, *, check: str) -> list[Finding]:
    """Run the battery over the PATH ITSELF, not only over the bytes it names.

    No rule had ever been applied to a path component, which left the commonest
    real-world carrier untouched: the identifier is in the filename long before it
    is in the file. A committed ``data/<accession>_<surname>.txt`` carries no
    genomic content and passed every content rule.

    The emitted finding is redacted like every other path (:func:`redact_path`),
    so reporting the problem does not reproduce it.
    """
    detected = redact_path(path) != path
    if not detected:
        return []
    return [
        Finding(
            check=check,
            path=path,
            line=None,
            rule_id="path_identifier",
            span_len=None,
            severity="fail",
            detail=(
                "The PATH carries an identifier: a keyword-anchored PHI token or a "
                "run of >=7 digits in a name component. Sequencing filenames are "
                "routinely the MRN, the NHS number or the accession. Rename the "
                "file to a case-local alias before it is committed."
            ),
        )
    ]


def check_content_scan(repo_root: Path, *, strict: bool = False) -> CheckResult:
    """The rule battery over every file git would let you commit."""
    name = "content_scan"
    if not _git_available(repo_root):
        return _unavailable(name, "not a git repository (or git unavailable)", strict=strict)
    paths = _scannable_paths(repo_root)
    load = worktree_loader(repo_root)
    findings: list[Finding] = []
    scanned = 0
    for path in paths:
        candidate = repo_root / path
        if not candidate.is_file() or candidate.is_symlink():
            continue
        if any(part in SKIP_DIR_NAMES for part in Path(path).parts):
            continue
        scanned += 1
        findings.extend(_path_findings(path, check=name))
        findings.extend(
            scan_file(
                candidate,
                check=name,
                path_label=path,
                allowlisted=_exempt(path, CONTENT_ALLOWLIST_PREFIXES),
                hpo_is_constant=_exempt(path, HPO_CODE_PREFIXES),
                reference_fixture=reference_fixture_provenance(path, load=load).exempt,
            )
        )
    fails = sum(1 for f in findings if f.severity == "fail")
    return _result(
        name,
        findings,
        f"{scanned} file(s) scanned (cap {MAX_SCAN_BYTES} bytes); {fails} fail-severity hit(s).",
        strict=strict,
    )


def check_synthetic_fixtures_marked(repo_root: Path, *, strict: bool = False) -> CheckResult:
    """Everything under ``tests/fixtures/synthetic/`` must declare itself synthetic.

    This directory is the one place ``.gitignore`` deliberately re-admits ``.vcf``,
    ``.ped`` and ``.bed`` files. Without this check, the safest-looking way to leak
    a real patient VCF is to save it here, because every other control assumes the
    directory is synthetic. The marker requirement inverts that: an unmarked file
    FAILS rather than sliding in.
    """
    name = "synthetic_fixtures_marked"
    root = repo_root / "tests" / "fixtures" / "synthetic"
    if not root.is_dir():
        return _result(name, [], "No synthetic fixture directory present.", strict=strict)

    findings: list[Finding] = []
    checked = 0
    for candidate in sorted(root.rglob("*")):
        if not candidate.is_file() or candidate.name == ".DS_Store":
            continue
        checked += 1
        label = candidate.relative_to(repo_root).as_posix()
        try:
            data = read_capped(candidate, MAX_SCAN_BYTES)
        except OSError:
            data = b""
        plain = gunzip_capped(data, MAX_SCAN_BYTES)
        if plain is not None:
            data = plain
        if not synthetic_marker_present(data):
            findings.append(
                Finding(
                    check=name,
                    path=label,
                    line=None,
                    rule_id=None,
                    span_len=None,
                    severity="fail",
                    detail=(
                        "No synthetic marker in the first "
                        f"{SYNTHETIC_MARKER_BYTES} bytes. Add `mva_synthetic=true` as "
                        "a VCF header line or a leading comment. A declaration is a "
                        "header, not a needle: accepting the marker anywhere in the "
                        "file meant one incidental occurrence at any depth counted. "
                        "If the file is not synthetic it must not be here at all."
                    ),
                )
            )
            continue
        unmarked = unmarked_sample_names(data)
        if unmarked:
            findings.append(
                Finding(
                    check=name,
                    path=label,
                    line=None,
                    rule_id=None,
                    span_len=None,
                    severity="fail",
                    detail=(
                        f"Declares itself synthetic, but {unmarked} named subject(s) "
                        "(VCF sample column or @RG SM: value) do not carry a SYNTH_ "
                        "prefix. The marker comment is a claim; the sample names are "
                        "the field a real export cannot satisfy by accident. Sample "
                        "names are counted, never emitted."
                    ),
                )
            )
    return _result(name, findings, f"{checked} fixture file(s) checked.", strict=strict)


def check_reference_fixtures_declared(repo_root: Path, *, strict: bool = False) -> CheckResult:
    """Every file under a public-reference fixture prefix declares its descent.

    ADR 0012 conditions 1, 2 and 4, enforced the way ``synthetic_fixtures_marked``
    enforces the synthetic marker: an UNDECLARED file **fails** rather than merely
    losing an exemption it never asked for.

    The distinction is the whole point of the check. Without it, provenance is
    only ever consulted when something else already wants to soften a finding, so
    a file that trips no other rule -- a sites-only export with a harmless
    extension, a stray index, a second copy of a slice -- enters these directories
    in silence, and the ``.gitignore`` negations that re-admit them are the widest
    doors in the repository. With it, joining the category costs a named entry in
    a committed file, beside a committed generator, against a resource the
    manifest records having fetched with a real digest.

    What it does NOT establish, stated plainly because the ADR's condition 2 is
    easy to over-read: naming a resource is a CLAIM. Nothing available offline
    proves a given slice was cut from ClinVar rather than from a patient. What the
    requirement removes is the *silent* path in. The property that actually
    protects an individual is condition 3 -- no sample columns, no genotypes --
    which is a fact about the bytes and is checked against them
    (:func:`buffer_is_sites_only_vcf`), not here.
    """
    name = "reference_fixtures_declared"
    load = worktree_loader(repo_root)
    findings: list[Finding] = []
    checked = 0
    for prefix in PUBLIC_REFERENCE_FIXTURE_PREFIXES:
        directory = repo_root / prefix.rstrip("/")
        if not directory.is_dir():
            continue
        doc, reason = _load_fixture_provenance(prefix, load=load)
        if doc is None:
            findings.append(
                Finding(
                    check=name,
                    path=prefix,
                    line=None,
                    rule_id=None,
                    span_len=None,
                    severity="fail",
                    detail=(
                        "This directory re-admits files that .gitignore otherwise "
                        f"denies, and it has no usable provenance record: {reason}. "
                        f"Add {FIXTURE_PROVENANCE_FILENAME} naming, per committed "
                        "file, the manifest resource it was cut from and the regions "
                        "that were cut (ADR 0012 conditions 2 and 4). Until then no "
                        "file here can earn the public-reference exemption."
                    ),
                )
            )
            continue
        exempt_names = {FIXTURE_PROVENANCE_FILENAME, doc.generator, ".DS_Store"}
        for candidate in sorted(directory.rglob("*")):
            if not candidate.is_file() or candidate.name in exempt_names:
                continue
            relative = candidate.relative_to(directory)
            if any(part in SKIP_DIR_NAMES for part in relative.parts):
                continue
            checked += 1
            label = candidate.relative_to(repo_root).as_posix()
            entry = doc.fixtures.get(candidate.name) if len(relative.parts) == 1 else None
            if entry is not None and entry.indexes and not (directory / entry.indexes).is_file():
                findings.append(
                    Finding(
                        check=name,
                        path=label,
                        line=None,
                        rule_id=None,
                        span_len=None,
                        severity="fail",
                        detail=(
                            "Orphan index: it declares the fixture it indexes, and "
                            "that fixture is not here. An index inherits its data "
                            "file's verdict, so an index with no data file inherits "
                            "nothing and must not be committed alone."
                        ),
                    )
                )
                continue
            verdict = reference_fixture_provenance(label, load=load)
            if verdict.exempt:
                continue
            findings.append(
                Finding(
                    check=name,
                    path=label,
                    line=None,
                    rule_id=None,
                    span_len=None,
                    severity="fail",
                    detail=(
                        "File in a public-reference fixture directory without "
                        f"verifiable provenance: {verdict.reason}. Declare it in "
                        f"{FIXTURE_PROVENANCE_FILENAME} against a resource the "
                        "manifest records as fetched with a real sha256, or remove "
                        "it. A directory whose .gitignore negations re-admit genomic "
                        "formats is exactly where an undeclared file must fail."
                    ),
                )
            )
    return _result(
        name,
        findings,
        f"{checked} reference fixture file(s) checked for declared provenance.",
        strict=strict,
    )


def check_notebook_output_purity(repo_root: Path, *, strict: bool = False) -> CheckResult:
    """No tracked notebook may carry executed outputs.

    A notebook output cell is a verbatim dump of whatever the cell printed: a
    DataFrame head of variant rows, a rendered pedigree, a stack trace with frame
    locals. It is also invisible in a normal diff review.
    """
    name = "notebook_output_purity"
    if not _git_available(repo_root):
        return _unavailable(name, "not a git repository (or git unavailable)", strict=strict)
    code, out = _git(repo_root, ["ls-files", "-z", "--", "*.ipynb"])
    if code != 0:
        return _unavailable(name, "git ls-files failed", strict=strict)

    findings: list[Finding] = []
    paths = sorted(_split_nul(out))
    for path in paths:
        candidate = repo_root / path
        if not candidate.is_file():
            continue
        try:
            # decode_scrubbed, not bytes-into-json: json.loads on bytes raises a
            # UnicodeDecodeError whose str() embeds the offending bytes, and that
            # string travels to the terminal, the CI log and a model context. The
            # guard only protects anything if real decodes are routed through it.
            text = decode_scrubbed(read_capped(candidate, MAX_SCAN_BYTES), path=candidate)
            document = cast(dict[str, Any], json.loads(text))
        except (ValueError, PrivacyViolationError) as exc:
            findings.append(
                Finding(
                    check=name,
                    path=path,
                    line=None,
                    rule_id=None,
                    span_len=None,
                    severity="warn",
                    detail=f"Not parseable as a notebook ({type(exc).__name__}).",
                )
            )
            continue
        cells = document.get("cells")
        if not isinstance(cells, list):
            continue
        for index, cell in enumerate(cast(list[Any], cells)):
            if not isinstance(cell, dict):
                continue
            outputs = cast(dict[str, Any], cell).get("outputs")
            if isinstance(outputs, list) and outputs:
                findings.append(
                    Finding(
                        check=name,
                        path=path,
                        line=None,
                        rule_id=None,
                        span_len=None,
                        severity="fail",
                        detail=(
                            f"Cell {index} carries {len(cast(list[Any], outputs))} stored "
                            "output(s). Strip outputs before committing."
                        ),
                    )
                )
    return _result(name, findings, f"{len(paths)} tracked notebook(s) checked.", strict=strict)


# ---------------------------------------------------------------------------
# Runtime probe (GP-42)
#
# The canary payloads below are ASSEMBLED AT RUNTIME from fragments. Written as
# plain literals they would make this very file trip `content_scan`, and lowering
# the rules to accommodate the audit's own source would be exactly the kind of
# self-serving weakening the audit exists to prevent. Assembly also documents the
# rules' shape: a bare position or a bare digit run is intentionally not
# detectable, so the fragments alone are inert.
# ---------------------------------------------------------------------------

_CANARY_POS: Final[str] = "40200000"
_CANARY_ID: Final[str] = "4457821"
_CANARY_HPO_DIGITS: Final[str] = "0001250"
_CANARY_SEQ: Final[str] = "ACGT" * 12


def canary_payloads() -> tuple[tuple[str, str, tuple[str, ...]], ...]:
    """(label, payload, substrings that must NOT survive redaction)."""
    vcf_line = "\t".join(
        ("chr15", _CANARY_POS, ".", "C", "T", "820.5", "PASS", "SYNTH", "GT:DP", "0/1:45")
    )
    fastq = "\n".join(("@synthetic_read_1", _CANARY_SEQ, "+", "I" * len(_CANARY_SEQ)))
    hpo = f"HP:{_CANARY_HPO_DIGITS}"
    mrn = f"MRN: {_CANARY_ID}"
    return (
        ("vcf_data_line", vcf_line, (f"chr15\t{_CANARY_POS}", _CANARY_POS)),
        ("fastq_record", fastq, (_CANARY_SEQ,)),
        ("hpo_term", hpo, (hpo,)),
        ("mrn", mrn, (_CANARY_ID,)),
    )


class _CaptureHandler(logging.Handler):
    """Collects what a handler would have written, without writing it."""

    def __init__(self, sink: list[str]) -> None:
        super().__init__(level=logging.DEBUG)
        self._sink = sink

    def emit(self, record: logging.LogRecord) -> None:
        try:
            self._sink.append(self.format(record))
        except Exception:  # a formatting failure must not abort the audit
            self._sink.append("<format-error>")


def _registry_handlers() -> list[logging.Handler]:
    """Every handler attached ANYWHERE in the logger registry, root's first.

    The probe swaps each one's ``emit`` for a capture, so containment of the
    canaries is structural rather than a consequence of the control being
    measured. Looking only at ``root.handlers`` left a handler attached to a
    library's own logger free to write a canary to the terminal for real, and the
    probe covered that gap by arming redaction first -- which is precisely the
    side effect that made this check report a different answer on its second run
    in the same process.
    """
    handlers: list[logging.Handler] = list(logging.getLogger().handlers)
    seen = {id(handler) for handler in handlers}
    for entry in list(logging.Logger.manager.loggerDict.values()):
        if not isinstance(entry, logging.Logger):
            continue  # a PlaceHolder, which holds no handlers
        for handler in entry.handlers:
            if id(handler) not in seen:
                seen.add(id(handler))
                handlers.append(handler)
    return handlers


def check_log_redaction_probe(*, strict: bool = False) -> CheckResult:
    """Push real canaries through the CONFIGURED logging stack and look for survivors.

    This is a live probe rather than a static assertion, because every static
    argument for "our logs are clean" has failed in practice: the filter was
    attached to a logger instead of a handler, a library added its own handler
    after configuration, a formatter cached ``exc_text``. The only reliable test is
    to log something that must not appear and then read every handler's output.

    Each handler's ``emit`` is temporarily swapped for a capture that formats the
    record and discards it. Filters still run — they are applied in
    ``Handler.handle`` before ``emit`` — so this measures the real pipeline while
    keeping the canaries off the terminal, out of the log file, and out of any
    model context reading this run.

    The state of the pipeline is recorded BEFORE anything is installed, across the
    whole logger registry rather than the root logger alone. Both were defects in
    their own right: a probe that installs the control it measures can only report
    success, and a probe that inspects only ``root.handlers`` reports "clean" for
    the library-attached handler that is the actual leak.
    """
    name = "log_redaction_probe"
    canaries = canary_payloads()
    captured: list[str] = []
    root = logging.getLogger()
    probe = _CaptureHandler(captured)

    original_level = root.level
    swapped: list[tuple[logging.Handler, bool, Any]] = []

    root.addHandler(probe)
    for handler in _registry_handlers():
        if handler is probe:
            continue
        had_own = "emit" in handler.__dict__
        swapped.append((handler, had_own, handler.__dict__.get("emit")))

        def _capture_for(bound: logging.Handler) -> Callable[[logging.LogRecord], None]:
            def _emit(record: logging.LogRecord) -> None:
                try:
                    captured.append(bound.format(record))
                except Exception:  # never let a formatter abort the probe
                    captured.append("<format-error>")

            return _emit

        cast(Any, handler).emit = _capture_for(handler)

    # OBSERVE, AND CHANGE NOTHING. The probe used to call install_redaction()
    # here, which cost it both of the properties it needs:
    #
    # * it could not report the one failure it exists to report -- that the
    #   application never armed GP-42 -- because by the time it looked it had
    #   repaired the stack itself. That half was fixed by reading the state
    #   first; the install stayed, and left the second half broken:
    # * it MUTATED what it measures, so its answer depended on execution history
    #   rather than on inputs. The first audit in a process reported "NOT armed"
    #   and armed it; every later audit in that process reported "armed". The
    #   report is a registered artifact (`privacy/privacy_audit.md`), so two
    #   pipeline runs in one process produced two different files and GP-30's
    #   byte-identity claim was false -- in exactly the situation where someone
    #   re-runs a submission to check reproducibility.
    #
    # Containment of the canaries does not need the install: every handler in the
    # registry has had its `emit` swapped for a capture above. Arming GP-42 is
    # the composition root's job (`mva.cli._install_privacy_guards`), which is
    # what this check exists to verify it did.
    armed_before = redaction_installed()
    unfiltered = unfiltered_handlers()

    logger = logging.getLogger("mva.privacy.audit.canary")
    logger.setLevel(logging.NOTSET)
    logger.propagate = True
    root.setLevel(logging.DEBUG)
    try:
        for _, payload, _ in canaries:
            # Both paths matter: %-args are scrubbed separately from record.msg.
            logger.debug("canary via args: %s", payload)
            logger.debug(payload)
    finally:
        root.setLevel(original_level)
        root.removeHandler(probe)
        for handler, had_own, previous in swapped:
            if had_own:
                cast(Any, handler).emit = previous
            else:
                handler.__dict__.pop("emit", None)

    haystack = "\n".join(captured)
    findings: list[Finding] = [
        Finding(
            check=name,
            path="<logging>",
            line=None,
            rule_id=rule_id,
            span_len=len(needle),
            severity="fail",
            detail=(
                "A canary substring survived the logging pipeline. Redaction is not "
                "in force for at least one handler."
            ),
        )
        for rule_id, _, needles in canaries
        for needle in needles
        if needle in haystack
    ]
    if not captured:
        findings.append(
            Finding(
                check=name,
                path="<logging>",
                line=None,
                rule_id=None,
                span_len=None,
                severity="warn",
                detail="No handler output captured; the probe was inconclusive.",
            )
        )
    if not armed_before:
        findings.append(
            Finding(
                check=name,
                path="<logging>",
                line=None,
                rule_id=None,
                span_len=None,
                severity="fail",
                detail=(
                    "Log redaction was NOT armed when the audit started. The "
                    "composition root must call mva.privacy.redact.install_redaction() "
                    "before any stage runs; every record logged before that point "
                    "left the process unscrubbed. The probe has armed it now, which "
                    "does nothing for the records already written."
                ),
            )
        )
    findings.extend(
        Finding(
            check=name,
            path="<logging>",
            line=None,
            rule_id=None,
            span_len=None,
            severity="warn",
            detail=(
                f"Handler {handler_name} carried no GenomicRedactionFilter when the "
                "audit started. Counted across the whole logger registry, not just "
                "the root: a handler a library attached to its own logger is never "
                "consulted by the root and was previously invisible here."
            ),
        )
        for handler_name in unfiltered
    )
    return _result(
        name,
        findings,
        f"{len(canaries)} canary payload(s) pushed at DEBUG; "
        f"{len(captured)} handler emission(s) captured.",
        strict=strict,
    )


def check_cloud_sync_location(
    repo_root: Path, *, workspace: Path | None = None, strict: bool = False
) -> CheckResult:
    """The workspace must not sit under a cloud-synced root.

    ``~/Desktop`` and ``~/Documents`` are included because macOS syncs both to
    iCloud Drive by default ("Desktop & Documents Folders"). A VCF written there is
    uploaded within seconds, is outside the researcher's control from that moment,
    and cannot be recalled. Findings name the matched marker, never the path.
    """
    name = "cloud_sync_location"
    findings: list[Finding] = []
    home = Path.home().resolve()

    resolved = _configured_workspace(workspace)
    if resolved is not None:
        posix = resolved.as_posix()
        findings.extend(
            Finding(
                check=name,
                path=WORKSPACE_LABEL,
                line=None,
                rule_id=None,
                span_len=None,
                severity="fail",
                detail=(
                    f"Workspace is under a cloud-synced root ({marker!r}). Patient "
                    "data placed there is uploaded to a third party automatically."
                ),
            )
            for marker in CLOUD_SYNCED_MARKERS
            if marker in posix
        )
        for dirname in CLOUD_SYNCED_HOME_DIRS:
            synced = home / dirname
            if resolved == synced or resolved.is_relative_to(synced):
                findings.append(
                    Finding(
                        check=name,
                        path=WORKSPACE_LABEL,
                        line=None,
                        rule_id=None,
                        span_len=None,
                        severity="fail",
                        detail=(
                            f"Workspace is under ~/{dirname}, which macOS syncs to "
                            "iCloud Drive by default. Use a non-synced directory."
                        ),
                    )
                )

    repo_posix = repo_root.resolve().as_posix()
    findings.extend(
        Finding(
            check=name,
            path=".",
            line=None,
            rule_id=None,
            span_len=None,
            severity="warn",
            detail=(
                f"The repository itself is under a cloud-synced root ({marker!r}). "
                "Audit reports and logs written here leave the machine."
            ),
        )
        for marker in CLOUD_SYNCED_MARKERS
        if marker in repo_posix
    )
    summary = (
        "Workspace and repository checked against cloud-sync markers."
        if resolved is not None
        else "No workspace configured; only the repository was checked."
    )
    return _result(name, findings, summary, strict=strict)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

#: Checks that inspect the index/worktree only — the pre-commit subset.
STAGED_CHECKS: Final[tuple[str, ...]] = (
    "git_staged_sensitive",
    "gitignore_effectiveness",
    "gitignore_negation_safety",
    "synthetic_fixtures_marked",
    "reference_fixtures_declared",
    "notebook_output_purity",
    "log_redaction_probe",
)


def run_audit(
    repo_root: Path,
    *,
    workspace: Path | None = None,
    staged_only: bool = False,
    strict: bool = False,
) -> AuditReport:
    """Run the privacy audit.

    ``staged_only`` runs the fast pre-commit subset (:data:`STAGED_CHECKS`): what is
    about to be committed, plus the ignore-file and logging invariants. It skips the
    whole-tree content scan and the symlink walk, which are the slow checks and
    which say nothing about *this* commit.

    ``strict`` promotes warnings to failures. Off by default because the
    highest-value warnings (``vcf_data_line`` in a public coordinate table, ISO
    dates everywhere) are dominated by legitimate matches; a strict run is for CI
    on a release branch, not for the inner loop.

    ``workspace`` defaults to ``$MVA_WORKSPACE``. When neither is set, the two
    workspace checks report a warning and are skipped rather than failing: an
    audit of the repository alone is a legitimate thing to want.
    """
    root = repo_root.resolve()
    results: list[CheckResult] = []
    if staged_only:
        results = [
            check_git_staged_sensitive(root, strict=strict),
            check_gitignore_effectiveness(root, strict=strict),
            check_gitignore_negation_safety(root, strict=strict),
            check_synthetic_fixtures_marked(root, strict=strict),
            check_reference_fixtures_declared(root, strict=strict),
            check_notebook_output_purity(root, strict=strict),
            check_log_redaction_probe(strict=strict),
        ]
    else:
        results = [
            check_git_tracked_sensitive(root, strict=strict),
            check_git_staged_sensitive(root, strict=strict),
            check_gitignore_effectiveness(root, strict=strict),
            check_gitignore_negation_safety(root, strict=strict),
            check_workspace_containment(root, workspace=workspace, strict=strict),
            check_symlink_escape(root, workspace=workspace, strict=strict),
            check_content_scan(root, strict=strict),
            check_synthetic_fixtures_marked(root, strict=strict),
            check_reference_fixtures_declared(root, strict=strict),
            check_notebook_output_purity(root, strict=strict),
            check_log_redaction_probe(strict=strict),
            check_cloud_sync_location(root, workspace=workspace, strict=strict),
        ]

    failed = tuple(result.name for result in results if not result.passed)
    return AuditReport(results=tuple(results), passed=not failed, failed_checks=failed)
