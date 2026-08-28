"""Population frequencies from local gnomAD sites VCFs — a REAL frequency adapter.

This replaces :class:`~mva.annotation.local_tables.LocalFrequencyAdapter` (invented
numbers, ``synthetic-v0.0``) with tabix region queries against the gnomAD sites
VCFs that the offline acquisition step downloaded. It is the highest-value half of
TD-01: allele frequency is the strongest rarity signal the ranking has, and until
now it was fabricated.

Seven things here are load-bearing, and each exists because the obvious shortcut is
a silent, plausible-looking wrong answer.

**Absence stays absent (GP-14).** A variant with no gnomAD record is *omitted* from
the returned mapping. It is never returned with ``allele_frequency=0.0``. gnomAD
draws the same distinction the pipeline does: ``AC=0`` over ``AN=403,848`` at
``chr21:5031905`` is a real observation of zero carriers in 400,000 chromosomes,
while a site gnomAD never called is no evidence at all. Defaulting absence to zero
converts "the reference cohort covers this ancestry badly" into "ultra-rare
candidate", and it does so most often exactly where reference panels are thinnest.
The same rule applies one level down: a genetic ancestry group with ``AN=0`` at a
site has no ``AF_<grp>`` key in the record at all (verified against v4.1 exomes —
see ``chr21:6086421 G>T``, where six of nine groups are absent and three report a
genuine zero), and that group is omitted rather than emitted at zero.

**"We could not look" is a third fact, and it raises.** The rule above governs
variants gnomAD was *asked about*. A variant on a contig with no complete shard
was never asked about at all, and returning it as a missing key would collapse the
two: a common allele on a chromosome whose 3 GB shard is still downloading would
look identical to a variant gnomAD has never seen, evade the common-frequency
down-rank, and be ranked as an ultra-rare candidate. :meth:`frequencies
<GnomadSitesFrequencyAdapter.frequencies>` therefore raises on a coverage hole.
Running against a partial release is possible but must be *asked for*, via
:meth:`lookup_partial <GnomadSitesFrequencyAdapter.lookup_partial>`, which returns
the gap inside a :class:`FrequencyLookup` rather than as an adapter property a
caller can forget to read.

**Error paths leak too.** cyvcf2 puts the region string — a proband coordinate —
into the message of anything it raises, and an unwrapped backend failure therefore
prints ``chr15:40200239-40200239`` to the terminal, the log, a crash report and an
agent's context. Region queries are wrapped, the replacement exception carries a
one-way :func:`~mva.models.base.error_token` handle instead of the coordinate, and
the original is raised ``from None`` so chaining cannot print it again (PRIV-09).

**Filtered sites are carried, not dropped.** gnomAD's ``AC0`` / ``AS_VQSR`` /
``InbreedingCoeff`` flags mean "we saw this and distrust the call". Dropping those
records would turn that into "gnomAD never saw it" — the requirement above, failed
by another route, and on this dataset it is not a hypothetical: a full scan of
chr21 (2,188,842 records) found **no** ``PASS`` record with ``AC=0``, because zero
post-QC alleles is precisely what the ``AC0`` filter marks. Every genuine
"observed in 400,000 chromosomes, carried by none of them" record — the strongest
rarity evidence the dataset can give — is a filtered record. The verbatim FILTER
string travels in :attr:`~mva.models.variant.PopulationFrequency.filter_status`,
where a reader and a downstream rule can weigh it.

**PASS and "." are different facts.** cyvcf2 reports both a passing record's
FILTER and an unfiltered record's ``.`` as ``None``, so this module reads
``FILTERS`` (the list) instead: ``["PASS"]`` becomes the literal string ``"PASS"``,
an empty list becomes ``None``. ``None`` in ``filter_status`` therefore means "the
file recorded no FILTER opinion", which is a different claim from "the file said
this call passed", and writing the second where the first is true would be
inventing a quality judgement the source never made.

**Per-population, not just global.** The ranking takes a population *maximum*
under an allele-number guard (ADR 0010). That guard reads ``allele_number``, so
per-group ``AC``/``AN`` must be present and correct or it silently does nothing and
a founder allele in an under-sampled cohort gets destroyed. ``chr21:5036078 A>C``
is the case in real data: global AF 0.0037, ``afr`` AF 0.170 over AN 15,066. The
set of genetic ancestry groups is read out of the VCF header
(``##INFO=<ID=AF_afr,...``), not hard-coded, so a release that adds or renames a
group is followed rather than silently truncated.

**One canonicalisation rule, never a second copy (ADR 0018).** Both sides of the
join — the caller's variant IDs and every ALT read out of the release — are put
through :func:`mva.alleles.canonicalise_allele`, the same function
:mod:`mva.ingestion.normalise` and
:class:`~mva.annotation.clinvar_vcf.ClinvarVcfAdapter` call. This module keeps no
trimming or shifting logic of its own. It used to: a private
``minimal_representation`` that trimmed but could not left-align, which meant a
left-aligned query key and a gnomAD record spelled elsewhere in the same repeat
tract simply did not join — and a variant with no frequency record is scored as
novel and ultra-rare, the strongest promoting signal the ranker has. Two
implementations of one rule was the defect; deleting the second one is the fix.
The adapter therefore accepts an optional :class:`~mva.alleles.ReferenceLookup`,
left-aligns when it has one, and says so through
:attr:`~GnomadSitesFrequencyAdapter.representation_status` when it does not
(GP-14) rather than letting the miss read as "gnomAD has no record".

**Never the network (PRIV-05).** Both cyvcf2 and pysam will happily open an
``https://`` tabix URL and range-request a remote index. That capability is not
used and must not be: a proband coordinate in an outbound request is an
irreversible disclosure of patient genetic data. Every source path is checked for
a URI scheme *before* it is touched, and the check runs on the raw text because
``Path`` silently collapses ``https://host`` to ``https:/host`` — a guard written
against the two-slash form never fires. ``mva.annotation`` is additionally
forbidden from importing any network client by
``tests/unit/test_architecture.py::test_no_network_clients_in_sensitive_stages``.

**Never a half-written file, and never a half-written index.** The sites VCFs are
~250 GB and arrive over hours. Reading a truncated BGZF stream yields corrupt or —
worse — silently partial results, which look exactly like "this variant is novel".
:func:`check_source_complete` must pass before a file is opened: index present,
BGZF end-of-file marker present, the index actually reaching the end of the data,
and (optionally) size stable across a probe interval. Incomplete shards are
excluded and named in :attr:`GnomadSitesFrequencyAdapter.incomplete_sources` rather
than being read.

The index check is the one that is easy to leave out, because the EOF marker looks
like it has already answered the question. It has not: it proves the *data* stream
is whole and says nothing about the *index*, and htslib will region-query a
complete data file through an index built from a shorter one without complaint,
answering "no record" for everything past the point the index reaches. htslib does
warn when an index is older than its data — and on the real release it warns on all
24 shards, because the ``.tbi`` files finished downloading before the multi-gigabyte
shards did. Every one of those indexes is nevertheless complete, so mtime is not
the signal; :func:`index_covers_data` compares the index's furthest reachable block
offset against the file size instead, which is exact and which the real release
passes. See :func:`tabix_index_reach`.
"""

from __future__ import annotations

import re
import time
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Final, Protocol, cast

from mva.alleles import (
    CanonicalAllele,
    LeftAlignmentReport,
    LeftAlignmentStatus,
    ReferenceLookup,
    ReferenceStatus,
    canonicalise_allele,
    is_sequence_allele,
    rightmost_equivalent_bound,
    summarise_left_alignment,
)
from mva.annotation.bgzf import (
    has_bgzf_eof,
    index_covers_data,
    index_path_for,
)
from mva.errors import AdapterUnavailableError, GenomeBuildMismatchError, NetworkDeniedError
from mva.models.base import error_token
from mva.models.genome import GenomeBuild, contig_sort_key, normalise_contig
from mva.models.variant import PopulationFrequency

# --------------------------------------------------------------------------- names

#: Adapter identity. Names the *mechanism* (local sites VCF + tabix), because the
#: dataset and release it is pointed at are properties of the files, not of the code.
ADAPTER_NAME: Final = "gnomad-sites-vcf"

#: Label used for gnomAD's un-stratified ``AC``/``AN``/``AF``. gnomAD's header calls
#: it nothing at all (the fields are simply unsuffixed); ``global`` is this
#: repository's canonical name for that cohort and matches
#: ``knowledge/public/frequencies.tsv``.
GLOBAL_POPULATION: Final = "global"

#: Written into ``filter_status`` for a record whose FILTER column is ``PASS``.
PASS_FILTER: Final = "PASS"  # noqa: S105 - VCF FILTER value, not a credential

#: Sites-VCF extensions the acquisition step may produce.
_SOURCE_SUFFIXES: Final[tuple[str, ...]] = (".vcf.bgz", ".vcf.gz")

#: Positions within this many bases are fetched in one tabix region query rather
#: than one each. Purely a query-count optimisation: every fetched record is still
#: matched on the exact ``(pos, ref, alt)`` key, so a wider window can only cause
#: records to be read and discarded, never joined to the wrong variant.
#:
#: 1000 is tuned for this dataset rather than copied. gnomAD v4.1 exomes carries
#: ~47 records per kb in exons and ~8 KB per record, so an un-merged batch spends
#: one BGZF block decompression per variant; merging a clustered candidate set
#: into a handful of spans is what takes a 5,000-variant lookup from ~56 s to
#: under a second. Spans grow only by chaining, and a chained span of length L
#: necessarily contains at least L/window queried positions, so the window cannot
#: be tricked into reading a whole chromosome for two distant variants.
DEFAULT_MERGE_WINDOW_BP: Final = 1000

#: ``AF_<group>`` INFO IDs that are not genetic ancestry groups. ``grpmax`` is a
#: derived maximum (it would double-count the group that produced it) and ``raw``
#: is the pre-QC callset. Sex strata (``AF_XX``) and the UKB-excluded subset
#: (``AF_non_ukb``, ``AF_non_ukb_nfe``) are already excluded by the pattern, which
#: requires a single all-lowercase token terminated by the field separator.
_NON_ANCESTRY_AF_SUFFIXES: Final[frozenset[str]] = frozenset({"grpmax", "raw"})

#: Facts a gnomAD sites record carries that :class:`PopulationFrequency` has no
#: field for, recorded here rather than left implicit — an undocumented drop is
#: indistinguishable from a parser that never looked.
#:
#: The consequential one is the first. A record whose ``AN`` is 0 in every cohort
#: (33,596 of chr21's 2,188,842 records, all ``AC0``-filtered) is *dropped*, not
#: emitted, because ``allele_frequency`` is a required non-optional float and
#: there is no number to put in it. The dropped fact is "gnomAD attempted this
#: site and retained no genotype", which is weaker than "gnomAD has no record"
#: but is not evidence of rarity either, so collapsing the two is the
#: conservative direction. Emitting 0.0 would not be.
UNREPRESENTABLE_GNOMAD_FACTS: Final[tuple[str, ...]] = (
    "a site observed but with AN=0 in every cohort (no allele_frequency to record)",
    "faf95 / faf99 (ACMG filtering allele frequency, the CI lower bound)",
    "AC_grpmax / AF_grpmax / grpmax (gnomAD's own guarded population maximum)",
    "age_hist / GQ / DP / QUALapprox quality histograms",
    "AC_XX / AC_XY sex strata",
    "AC_non_ukb* (the UK-Biobank-excluded subset)",
    "vep, revel_max, spliceai_ds_max and the other per-transcript annotations",
)

_CONTIG_META_RE: Final = re.compile(r"^##contig=<ID=([^,>]+)[^>]*?assembly=([^,>]+)", re.MULTILINE)
_AF_GROUP_RE: Final = re.compile(r"^##INFO=<ID=AF_([a-z]+),", re.MULTILINE)
_FILTER_ID_RE: Final = re.compile(r"^##FILTER=<ID=([^,>]+)", re.MULTILINE)
_SIMPLE_META_RE: Final = re.compile(r"^##([A-Za-z0-9_]+)=(.*)$", re.MULTILINE)

#: A ``v``-prefixed dotted release token as it appears in a gnomAD filename
#: (``gnomad.exomes.v4.1.sites.chr21.vcf.bgz`` -> ``v4.1``). Used to check the
#: declared release against the filename by whole token: a plain substring test
#: would accept ``release="v4"`` for a ``v4.1`` download and ``release="v4.1"``
#: for a ``v4.10`` one.
_FILENAME_RELEASE_RE: Final = re.compile(r"v[0-9]+(?:\.[0-9]+)*")

#: Anything of the form ``scheme:/`` at the start of a path. Deliberately matches
#: the single-slash form: ``Path("https://host/x")`` stringifies as
#: ``"https:/host/x"``, so a guard testing for ``"://"`` never fires on a caller
#: who wrapped a URL in a ``Path``. A single leading letter is required so that a
#: Windows drive letter is not mistaken for a scheme.
_URI_SCHEME_RE: Final = re.compile(r"^[A-Za-z][A-Za-z0-9+.\-]{1,}:/")

#: Header keys that pin which pipeline run produced the release. gnomAD's VCF
#: header carries **no** release string — verified against the v4.1 exomes chr21
#: file, whose only version lines are the tool versions below plus
#: ``mane_select``/``pangolin``/``phylop``/``polyphen``/``seqrepo``/``sift``/``vrs``
#: — so these are the closest thing the data holds to a build identity: two files
#: that agree on all of them came out of the same production run.
_FINGERPRINT_KEYS: Final[tuple[str, ...]] = (
    "hailversion",
    "vep_version",
    "dbsnp_version",
    "gencode_version",
    "cadd_version",
    "revel_version",
    "spliceai_version",
)


# ------------------------------------------------------------------ reader surface


class _SitesRecord(Protocol):
    """The slice of a cyvcf2 ``Variant`` this adapter reads.

    ``FILTERS`` rather than ``FILTER``: the scalar renders both ``PASS`` and ``.``
    as ``None``, and those are different claims (see the module docstring).
    """

    @property
    def POS(self) -> int: ...  # noqa: N802 - cyvcf2's attribute name

    @property
    def REF(self) -> str: ...  # noqa: N802 - cyvcf2's attribute name

    @property
    def ALT(self) -> list[str]: ...  # noqa: N802 - cyvcf2's attribute name

    @property
    def FILTERS(self) -> list[str]: ...  # noqa: N802 - cyvcf2's attribute name

    @property
    def INFO(self) -> object: ...  # noqa: N802 - cyvcf2's attribute name


class _SitesReader(Protocol):
    """The slice of cyvcf2's ``VCF`` surface this adapter uses.

    cyvcf2 ships no type stubs, so without this the reader is ``Unknown`` and leaks
    into every signature that touches it. Naming the three members actually used
    also documents the dependency: header text, region query, close. Nothing here
    can open a URL — but see :func:`_reject_remote`, which is what enforces that.
    """

    @property
    def raw_header(self) -> str: ...

    def __call__(self, region: str) -> Iterator[_SitesRecord]: ...

    def close(self) -> None: ...


class _TabixHandle(Protocol):
    """The pysam surface used to read which contigs an index actually holds."""

    @property
    def contigs(self) -> list[str]: ...

    def close(self) -> None: ...


def _open_reader(path: Path) -> _SitesReader:
    """Open a sites VCF for region queries, or explain why we cannot.

    cyvcf2 is imported here rather than at module scope so that importing this
    module — which the architecture tests do — does not require the optional
    ``genomics`` extra to be installed.
    """
    try:
        from cyvcf2 import VCF  # noqa: PLC0415 - native backend, imported on demand
    except ImportError as exc:  # pragma: no cover - guarded by the genomics extra
        msg = (
            "The cyvcf2 backend is not installed; install the 'genomics' extra. Tabix "
            "region queries are how this adapter avoids scanning a 2 GB shard per variant."
        )
        raise AdapterUnavailableError(msg) from exc
    return cast(_SitesReader, VCF(str(path)))


def _open_tabix(path: Path) -> _TabixHandle:
    """Open the tabix index of ``path``. Same lazy-import rationale as above."""
    try:
        import pysam  # noqa: PLC0415 - native backend, imported on demand
    except ImportError as exc:  # pragma: no cover - guarded by the genomics extra
        msg = "The pysam backend is not installed; install the 'genomics' extra."
        raise AdapterUnavailableError(msg) from exc
    return cast(_TabixHandle, pysam.TabixFile(str(path)))


# --------------------------------------------------------------------------- paths


def _reject_remote(raw: object) -> None:
    """Refuse anything that names a URI scheme. PRIV-05, structurally.

    htslib accepts ``https://`` and ``s3://`` URLs and will range-request a remote
    tabix index; handing it a proband's coordinates would disclose patient genetic
    data to a third party irreversibly. The capability is removed here rather than
    governed by convention.

    Raises :class:`~mva.errors.NetworkDeniedError` rather than
    :class:`~mva.errors.AdapterUnavailableError` on purpose: this is a privacy
    control, not a missing resource, and a caller that catches
    ``AdapterUnavailableError`` to fall back to a synthetic adapter must not
    swallow it.
    """
    text = str(raw)
    if _URI_SCHEME_RE.match(text) is None:
        return
    scheme = text.split(":", 1)[0]
    msg = (
        f"Refusing a non-local gnomAD source with URI scheme {scheme!r}. This adapter opens "
        "local files only: htslib would range-request the remote index and send proband "
        "coordinates to a third party (PRIV-05). Download the release with the offline "
        "acquisition step first. The path is not echoed."
    )
    raise NetworkDeniedError(msg)


def _local_source_path(path: Path) -> Path:
    """Resolve a path and refuse anything that is not an existing local file."""
    _reject_remote(path)
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        msg = (
            f"gnomAD sites file {resolved.as_posix()!r} is not a readable local file. "
            "The adapter never fetches anything; point it at a downloaded release."
        )
        raise AdapterUnavailableError(msg)
    return resolved


# --------------------------------------------------------------- completeness check


@dataclass(frozen=True, slots=True)
class SourceCompleteness:
    """Whether one sites file is safe to open, and if not, why not.

    Reasons are file *names* and structural facts only — never record content
    (PRIV-09). This is a value object rather than a bare bool so that a caller can
    report which shard is still arriving instead of silently annotating against a
    dataset with a hole in it.
    """

    path: Path
    exists: bool
    has_index: bool
    has_bgzf_eof: bool
    size_bytes: int
    #: ``None`` when no stability probe was requested — unknown, not stable.
    size_stable: bool | None
    #: Whether the index reaches the end of the data. ``None`` when it could not be
    #: measured (a CSI index, or an unparseable one) — unknown, not broken. See
    #: :func:`index_covers_data`, and note that this is a content check: the real
    #: release's indexes are all *older* than their data and all complete.
    index_covers_data: bool | None = None

    @property
    def is_complete(self) -> bool:
        """True only when every check that ran passed.

        Two of the four signals are tri-state, and ``None`` means "not measured"
        rather than "failed" in both. ``size_stable is None`` (not probed) does not
        block: the BGZF end-of-file marker is the authoritative signal for the data
        stream, and a stream that carries it is a finished stream. Nor does
        ``index_covers_data is None``, which is what a CSI index reports: refusing a
        legitimate release we have no parser for would be a false negative of our
        own making. ``False`` from either one blocks.
        """
        return (
            self.exists
            and self.has_index
            and self.has_bgzf_eof
            and self.size_stable is not False
            and self.index_covers_data is not False
        )

    @property
    def reasons(self) -> tuple[str, ...]:
        """Human-readable failures, in a fixed order (GP-30).

        A file that does not exist reports only that: appending "no index" and
        "no BGZF marker" to it would bury the one fact the reader can act on.
        """
        if not self.exists:
            return ("file does not exist",)
        reasons: list[str] = []
        if not self.has_index:
            reasons.append("no .tbi/.csi index beside it")
        if not self.has_bgzf_eof:
            reasons.append("BGZF end-of-file marker missing (file is truncated or still arriving)")
        if self.size_stable is False:
            reasons.append("size changed during the stability probe (still being written)")
        if self.index_covers_data is False:
            reasons.append(
                "the index does not reach the end of the data (built from a shorter file; "
                "every record past its reach would report as absent)"
            )
        return tuple(reasons)


def check_source_complete(
    path: Path,
    *,
    stability_probe_seconds: float = 0.0,
    sleeper: Callable[[float], None] = time.sleep,
) -> SourceCompleteness:
    """Decide whether ``path`` may be opened as data.

    Three independent signals, because each misses a case the others catch:

    * an index beside the file — htslib cannot region-query without one, and the
      acquisition step fetches ``.tbi`` after the data;
    * the BGZF end-of-file marker — proves the compressed stream is whole;
    * optionally, size stability across ``stability_probe_seconds`` — catches a
      writer that has produced a valid tail and is still appending.

    The probe is **opt-in and defaults to off** so that construction does not sleep,
    and it takes an injected ``sleeper`` so a test can advance it without waiting.
    No wall-clock value is read (GP-30): the probe compares two ``st_size`` samples,
    it does not timestamp anything.
    """
    exists = path.is_file()
    if not exists:
        return SourceCompleteness(
            path=path,
            exists=False,
            has_index=False,
            has_bgzf_eof=False,
            size_bytes=0,
            size_stable=None,
        )

    size_before = path.stat().st_size
    stable: bool | None = None
    if stability_probe_seconds > 0.0:
        sleeper(stability_probe_seconds)
        stable = path.is_file() and path.stat().st_size == size_before

    index = index_path_for(path)
    return SourceCompleteness(
        path=path,
        exists=True,
        has_index=index is not None,
        has_bgzf_eof=has_bgzf_eof(path),
        size_bytes=size_before,
        size_stable=stable,
        index_covers_data=index_covers_data(path, index),
    )


# ------------------------------------------------------------------- header facts


@dataclass(frozen=True, slots=True)
class GnomadHeaderFacts:
    """What the sites VCF header actually says about itself.

    Deliberately not "what the filename says". The one thing this header does *not*
    carry is a gnomAD release string — v4.1 exome headers name Hail, VEP, dbSNP,
    GENCODE, CADD, REVEL and SpliceAI versions, and nothing else. See
    :class:`GnomadSitesFrequencyAdapter` for how the release is supplied and checked.
    """

    dataset: str
    """Dataset name from ``##contig ... assembly=`` (``gnomAD_GRCh38`` -> ``gnomAD``)."""

    build: GenomeBuild
    contigs: tuple[str, ...]
    """Contig names verbatim from the header, in header order."""

    ancestry_groups: tuple[str, ...]
    """Genetic ancestry groups the release reports, sorted; from ``AF_<grp>`` INFO IDs."""

    filter_ids: tuple[str, ...]
    tool_versions: tuple[tuple[str, str], ...]
    """``(key, value)`` pairs from :data:`_FINGERPRINT_KEYS`, sorted by key."""

    @property
    def fingerprint(self) -> str:
        """A stable identity for the production run that emitted these files.

        Two shards of one release agree on this string; a directory that mixes
        releases does not, which is the only way this adapter can detect a mixture
        (the header has no release field to compare).
        """
        tools = ",".join(f"{key}={value}" for key, value in self.tool_versions)
        return f"{self.dataset}/{self.build.value}/{tools}"


def _parse_header_facts(raw_header: str, *, path: Path) -> GnomadHeaderFacts:
    contig_meta = _CONTIG_META_RE.findall(raw_header)
    if not contig_meta:
        msg = (
            f"gnomAD sites file {path.name} has no '##contig=<ID=...,assembly=...>' header "
            "lines. The adapter reads the dataset name and genome build from that "
            "attribute rather than from the filename (GP-11); a header without it cannot "
            "be attributed to a build, and guessing GRCh38 would mis-locate every variant."
        )
        raise AdapterUnavailableError(msg)

    assemblies = {assembly for _, assembly in contig_meta}
    if len(assemblies) != 1:
        msg = (
            f"gnomAD sites file {path.name} declares {len(assemblies)} different contig "
            "assemblies. A single sites file must describe one assembly."
        )
        raise AdapterUnavailableError(msg)
    assembly = next(iter(assemblies))
    dataset, separator, build_token = assembly.rpartition("_")
    if not separator or not dataset:
        msg = (
            f"gnomAD sites file {path.name} declares assembly {assembly!r}, which carries no "
            "'<dataset>_<build>' form (expected e.g. 'gnomAD_GRCh38'). The dataset name "
            "reaches every PopulationFrequency.source and may not be invented by the adapter "
            "(GP-18)."
        )
        raise AdapterUnavailableError(msg)
    try:
        build = GenomeBuild.parse(build_token)
    except ValueError as exc:
        msg = f"gnomAD sites file {path.name} declares an unrecognised genome build: {exc}"
        raise AdapterUnavailableError(msg) from exc

    groups = sorted(set(_AF_GROUP_RE.findall(raw_header)) - _NON_ANCESTRY_AF_SUFFIXES)
    if not groups:
        msg = (
            f"gnomAD sites file {path.name} declares no 'AF_<group>' INFO fields, so it "
            "carries no per-population frequencies. The ranking takes a population maximum "
            "under an allele-number guard (ADR 0010); a global-only source would make that "
            "guard silently inert."
        )
        raise AdapterUnavailableError(msg)

    meta = dict(_SIMPLE_META_RE.findall(raw_header))
    tools = tuple(
        (key, meta[key].strip()) for key in sorted(_FINGERPRINT_KEYS) if meta.get(key, "").strip()
    )
    return GnomadHeaderFacts(
        dataset=dataset,
        build=build,
        contigs=tuple(contig for contig, _ in contig_meta),
        ancestry_groups=tuple(groups),
        filter_ids=tuple(sorted(set(_FILTER_ID_RE.findall(raw_header)))),
        tool_versions=tools,
    )


# --------------------------------------------------------------- query planning


def merge_query_regions(
    spans: Sequence[tuple[int, int]], window: int
) -> tuple[tuple[int, int], ...]:
    """Collapse 1-based inclusive spans into the fewest region queries.

    Two spans merge when the gap between them is at most ``window``. Purely a
    query-count optimisation: every fetched record is still matched on the exact
    ``(pos, ref, alt)`` key, so a wider window can only cause records to be read
    and discarded, never joined to the wrong variant.

    Takes spans rather than positions because a deletion's REF covers several
    bases and the region query has to reach the far end of it.
    """
    if not spans:
        return ()
    ordered = sorted(set(spans))
    merged: list[tuple[int, int]] = []
    low, high = ordered[0]
    for start, end in ordered[1:]:
        if start - high <= window:
            high = max(high, end)
        else:
            merged.append((low, high))
            low, high = start, end
    merged.append((low, high))
    return tuple(merged)


@dataclass(frozen=True, slots=True)
class _Query:
    """One caller-supplied variant ID, decomposed into a joinable coordinate.

    ``position``/``ref``/``alt`` are the **canonical** form, not the raw text of
    the ID: they have been through :func:`mva.alleles.canonicalise_allele`, which
    is the same function every gnomAD record read out of the release is reduced
    by. Both sides of the join therefore agree by construction rather than by the
    caller having happened to normalise first.
    """

    variant_id: str
    contig: str
    position: int
    ref: str
    alt: str

    search_end: int
    """Right-most POS at which an equivalent spelling of this event could sit.

    Equal to ``position`` unless a reference is configured and the variant sits
    in a repeat tract. See :func:`mva.alleles.rightmost_equivalent_bound` and
    :attr:`span`, which is where it is actually spent."""

    reference_status: ReferenceStatus
    """How far the reference could be trusted while building *this* query.

    ``UNUSABLE`` when either the join key or ``search_end`` needed a base the
    reference could not supply, ``SHIFT_LIMIT_REACHED`` when every base was
    supplied and the search ran out of budget instead. Carried rather than
    discarded because the query side is where the reference is actually read: a
    gnomAD release is already left-aligned, so a release record almost never asks
    for a base, and an adapter that counted only the release side would report
    ``APPLIED`` over a FASTA that failed on every caller lookup."""

    is_indel: bool
    """Whether left-alignment was even applicable to this query.

    Carried on the query rather than recomputed from ``ref``/``alt`` at the
    counting site so that "is this an indel" is decided once, by the same rule
    that decided whether to shift it. It is what makes ``NOT_REQUIRED``
    distinguishable from ``APPLIED``: a batch of SNVs had nothing to left-align,
    which is a different statement from having left-aligned everything."""

    left_aligned: bool
    """Whether canonicalisation actually moved this query's POS leftwards."""

    @property
    def span(self) -> tuple[int, int]:
        """The 1-based inclusive window a region query must cover to be complete.

        Fixing the join key is only half the fix. An index query finds records by
        the bases they occupy, so a key that is now correct still loses its record
        if the fetch window never reaches the place the release spells it — and
        that miss is indistinguishable from "gnomAD has no record", which is the
        exact failure this adapter exists to not have.

        The window is ``[canonical POS, max(end of the canonical REF span,
        right-most equivalent POS)]``, and it is complete in both directions:

        * **Nothing is missed to the left.** With a reference, ``position`` is the
          *left-most* legal spelling of the event, so no record of it can sit
          further left. Without one, both sides are trimmed only and the window is
          byte-for-byte the one this adapter used before ADR 0018 — the degraded
          state is unchanged, not newly degraded.
        * **Non-minimal records to the left are already caught.** Trimming moves a
          POS rightwards but never outside the raw REF span, so a record whose
          canonical POS lands in this window necessarily *overlaps* this window
          with its raw span, and htslib returns every overlapping record.
          ``chr21:100 AT>AG`` is fetched by a query at 101.
        * **Right-shifted records need the right edge.** Left-alignment moves a
          record's POS leftwards, out of its raw span, so a release record
          spelling this insertion further along the tract occupies a span disjoint
          from ``position`` and would never be fetched at all.
          ``rightmost_equivalent_bound`` bounds that from the reference instead
          of from a guessed padding constant, and the bound is provably sufficient:
          a record's raw POS is never greater than its own trimmed POS, which is
          never greater than the right-most legal spelling of the event.

        The REF-span end is retained in the ``max`` because a deletion's REF covers
        several bases and, with no reference configured, ``search_end`` is just
        ``position``. Over-reaching costs bases of index query that are already
        inside a merged span; under-reaching re-opens the silent miss.
        """
        return (self.position, max(self.position + len(self.ref) - 1, self.search_end))


def _parse_variant_id(
    variant_id: str, *, build: GenomeBuild, reference: ReferenceLookup | None
) -> _Query:
    """Split a canonical variant ID, refusing anything from another assembly (GP-11).

    A GRCh37 coordinate looked up in a GRCh38 dataset does not return nothing — it
    returns whatever unrelated variant happens to sit at that offset, or nothing,
    and both answers are confidently wrong. Cross-build lookup raises rather than
    guessing, exactly as ``GenomicCoordinate.assert_same_build`` does.

    The caller's ID is canonicalised here rather than trusted. Ingestion already
    normalises proband records, so inside the pipeline this is usually a no-op —
    but the adapter is also called with hand-written and third-party IDs, and an
    adapter that joins correctly only when its caller happened to normalise first
    is an adapter whose correctness lives somewhere else.
    """
    parts = variant_id.split(":")
    if len(parts) != 5:
        msg = (
            f"Malformed variant ID <variant:{error_token(variant_id)}>: expected "
            "'{build}:{contig}:{pos}:{ref}:{alt}'. The value is tokenised rather than "
            "echoed (PRIV-09)."
        )
        raise ValueError(msg)
    build_token, contig, position_token, ref, alt = parts
    try:
        record_build = GenomeBuild.parse(build_token)
    except ValueError as exc:
        msg = (
            f"Variant ID <variant:{error_token(variant_id)}> names an unrecognised genome "
            f"build: {exc}"
        )
        raise ValueError(msg) from exc
    if record_build is not build:
        msg = (
            f"Variant ID <variant:{error_token(variant_id)}> is {record_build.value} but this "
            f"gnomAD release is {build.value}. Cross-build lookup is refused: the same offset "
            "names a different locus in each assembly, so the answer would be confidently "
            "wrong rather than missing (GP-11). Lift-over must be an explicit, "
            "provenance-tracked stage."
        )
        raise GenomeBuildMismatchError(msg)
    try:
        position = int(position_token)
    except ValueError as exc:
        msg = f"Variant ID <variant:{error_token(variant_id)}> has a non-integer position."
        raise ValueError(msg) from exc
    canonical_contig = normalise_contig(contig)
    canonical = canonicalise_allele(
        contig=canonical_contig,
        position=position,
        ref=ref.strip().upper(),
        alt=alt.strip().upper(),
        reference=reference,
    )
    search_end = canonical.position
    status = canonical.reference_status
    if reference is not None:
        # `.proven` is False when the reference could not be read to the right of
        # the record. The window is then the furthest point the search reached,
        # which may be short — but the key is trim-only for the same reason, so
        # the limitation is the one already reported by representation_status
        # rather than a second, undeclared one.
        bound = rightmost_equivalent_bound(
            contig=canonical_contig,
            position=canonical.position,
            ref=canonical.ref,
            alt=canonical.alt,
            reference=reference,
        )
        search_end = max(search_end, bound.position)
        # An unproven bound degrades the query, but *how* matters: an unreadable
        # reference and an exhausted shift budget send an operator to different
        # places, and only the first is a reason to go and look at the FASTA. The
        # worse of the two wins when the key and the bound disagree.
        if bound.reference_status is ReferenceStatus.UNUSABLE:
            status = ReferenceStatus.UNUSABLE
        elif (
            bound.reference_status is ReferenceStatus.SHIFT_LIMIT_REACHED
            and status is ReferenceStatus.USABLE
        ):
            status = ReferenceStatus.SHIFT_LIMIT_REACHED
    return _Query(
        variant_id=variant_id,
        contig=canonical_contig,
        position=canonical.position,
        ref=canonical.ref,
        alt=canonical.alt,
        search_end=search_end,
        reference_status=status,
        is_indel=(
            is_sequence_allele(canonical.ref)
            and is_sequence_allele(canonical.alt)
            and canonical.is_indel
        ),
        left_aligned=canonical.left_aligned,
    )


# ------------------------------------------------------------------ INFO accessors


def _as_int(value: object) -> int | None:
    """An INFO integer, or ``None`` when the key was absent. Never a default of 0."""
    if isinstance(value, bool) or not isinstance(value, int):
        return None
    return value


def _as_float(value: object) -> float | None:
    """An INFO float, or ``None`` when the key was absent. Never a default of 0.0."""
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return float(value)
    if isinstance(value, float):
        return value
    return None


def _per_allele(value: object, index: int) -> object:
    """Select this ALT's element from a ``Number=A`` INFO value.

    cyvcf2 returns a bare scalar when the record has one ALT and a tuple when it
    has several. gnomAD sites VCFs are split — a full scan of the v4.1 exomes
    chr21 shard found 0 multi-ALT records in 2,188,842 — so ``index`` is always 0
    in practice; handling the sequence case keeps the parser correct for any
    un-split source rather than silently attributing allele 1's frequency to
    allele 2.
    """
    if isinstance(value, (tuple, list)):
        sequence: Sequence[object] = value
        return sequence[index] if index < len(sequence) else None
    return value if index == 0 else None


def _get(info: object, key: str) -> object:
    """``INFO.get`` that returns ``None`` for an absent key rather than raising.

    cyvcf2 raises ``KeyError`` from ``INFO[key]`` and returns ``None`` from
    ``INFO.get(key)``; this wrapper exists so the untyped third-party surface is
    touched in exactly one place.
    """
    getter = getattr(info, "get", None)
    if getter is None:  # pragma: no cover - defensive
        return None
    return getter(key)


def _filter_status(filters: Sequence[str]) -> str | None:
    """The record's FILTER column, verbatim, or ``None`` when it held ``.``.

    cyvcf2 renders both ``PASS`` and ``.`` as ``FILTER = None``; ``FILTERS`` keeps
    them apart (``["PASS"]`` versus ``[]``). Writing ``"PASS"`` for an unfiltered
    record would attribute to the source a quality judgement it never made.
    Several filters are rejoined with ``;`` in file order, so
    ``AC0;AS_VQSR`` reaches the reader as it was written.
    """
    if not filters:
        return None
    return ";".join(str(item) for item in filters)


# --------------------------------------------------------------------- the adapter


@dataclass(frozen=True, slots=True)
class _ContigSource:
    """One contig, the file that holds it, and the contig name that file uses."""

    canonical: str
    raw: str
    path: Path


@dataclass(frozen=True, slots=True)
class FrequencyLookup:
    """A partial-coverage lookup: what was answered, and what could not be asked.

    Returned only by :meth:`GnomadSitesFrequencyAdapter.lookup_partial`, which is
    the *explicit opt-in* to running against an incomplete release. The gap is a
    field of the result rather than a property of the adapter on purpose: a
    property is something a caller can forget to read, and forgetting here means
    "we could not look" is silently rendered as "we looked and found nothing",
    which is GP-14's core failure mode. A common variant on a chromosome whose
    shard has not finished downloading would evade the common-frequency
    down-rank and be scored as rare.
    """

    frequencies: Mapping[str, tuple[PopulationFrequency, ...]]
    """Answers for the variants on contigs that *do* have a complete shard."""

    unqueryable_contigs: tuple[str, ...]
    """Queried contigs with no complete, indexed shard, in karyotype order."""

    unqueryable_variant_ids: tuple[str, ...]
    """The caller's IDs that could not be looked up at all, in caller order.

    Absent from :attr:`frequencies` and absent for a *different reason* from the
    variants gnomAD genuinely has no record for. Nothing may score these as rare.
    """

    incomplete_sources: tuple[str, ...]
    """Shards excluded as still-arriving, truncated or unreadable, with the reason."""

    @property
    def is_complete(self) -> bool:
        """True when every queried contig was backed by a complete shard."""
        return not self.unqueryable_contigs

    def describe_gap(self) -> str:
        """One sentence naming the coverage hole, for a warning or a report footer."""
        if self.is_complete:
            return "Every queried contig was backed by a complete gnomAD shard."
        return (
            f"{len(self.unqueryable_variant_ids)} variant(s) on "
            f"{', '.join(self.unqueryable_contigs)} could not be looked up: no complete "
            "gnomAD shard covers those contigs. This is absence of a resource, not "
            "absence of the variant from gnomAD, and must not be scored as rarity "
            "(GP-14)."
        )


class GnomadSitesFrequencyAdapter:
    """Allele frequencies read from local, tabix-indexed gnomAD sites VCFs.

    Constructed over a directory of shards (``gnomad.exomes.v4.1.sites.chr21.vcf.bgz``
    and friends) or over a single file. Each shard is completeness-checked before it
    is opened, its header is parsed, and the contigs it actually indexes are read
    from the tabix index — not guessed from the filename, and not taken from the
    header, which declares all 24 contigs in every shard.

    ``release`` is a required keyword argument and has no default. That is not
    laziness: **the gnomAD v4.1 sites VCF header contains no release string.** It
    names the Hail, VEP, dbSNP, GENCODE, CADD, REVEL and SpliceAI versions used to
    build it and nothing more (verified by
    ``test_the_release_header_line_does_not_exist``). Since ``version`` reaches
    every ``Citation`` and ``EvidenceItem`` refuses a ``DATABASE_ASSERTION`` whose
    citation has no version, the release has to come from the acquisition step that
    fetched the files — exactly as ``knowledge/manifests/knowledge.yaml`` records a
    table's version today. It is then *checked*, not trusted: it must appear as a
    whole version token in every source filename, and every shard's header
    fingerprint must agree, so pointing a ``v4.1`` adapter at a v2.1.1 download, or
    at a directory that mixes releases, fails loudly at construction instead of
    mislabelling every frequency it emits.

    There is deliberately **no** ``expected_sha256`` pin, which is where this
    adapter parts company with :class:`~mva.annotation.clinvar_vcf.ClinvarVcfAdapter`.
    ClinVar is one 185 MB file and hashing it costs a few seconds; a gnomAD exomes
    release is ~250 GB across 24 shards, so hashing at construction would add
    minutes to every run and would be skipped in practice — a pin nobody can
    afford is a pin nobody uses. The integrity story here is structural instead:
    :func:`check_source_complete` proves each BGZF stream is whole,
    :meth:`_check_filename` proves the release label matches the file, and
    :attr:`header_fingerprint` proves the shards came from one production run.
    A content hash belongs in the acquisition step's resource manifest, computed
    once at download time.
    """

    def __init__(
        self,
        sites: Path,
        *,
        release: str,
        subset: str = "exomes",
        require_contigs: Sequence[str] | None = None,
        merge_window_bp: int = DEFAULT_MERGE_WINDOW_BP,
        reference: ReferenceLookup | None = None,
        stability_probe_seconds: float = 0.0,
        sleeper: Callable[[float], None] = time.sleep,
    ) -> None:
        """Open a gnomAD sites release for lookup.

        Args:
            sites: directory of sites VCFs, or one sites VCF. Must be a local path;
                a URL is refused before it is touched.
            release: release identifier recorded by the acquisition step, e.g.
                ``"v4.1"``. Must appear as a whole token in every source filename.
            subset: which gnomAD callset these files are (``"exomes"`` /
                ``"genomes"``). Must also appear in every source filename. It becomes
                the second half of ``PopulationFrequency.source``; the first half is
                the dataset name read out of the header.
            require_contigs: contigs the caller cannot proceed without. Any that is
                missing or incomplete raises here, at construction, rather than at
                the first lookup that touches it. ``frequencies`` fails closed on a
                coverage hole either way; this makes the failure happen before a
                run has done any work.
            merge_window_bp: query-coalescing window; see :func:`merge_query_regions`.
            reference: optional 1-based inclusive reference accessor, the same
                :class:`~mva.alleles.ReferenceLookup` ingestion and the ClinVar
                adapter take. Trimming to a minimal representation never needs one
                and is always applied, so ``5031991 GAA>GA`` versus ``GA>G`` joins
                with or without it. Left-alignment does need one: gnomAD stores
                left-aligned alleles, so without a reference an indel spelled
                anywhere else in the same repeat tract cannot be reconciled, is
                reported as having no frequency record, and is then scored as
                novel and ultra-rare. :attr:`representation_status` states which
                of the two this adapter is in rather than leaving a reader to
                infer it from an absent frequency (GP-14).
            stability_probe_seconds: when > 0, re-stat each file after this delay and
                reject it if the size changed. Off by default so construction does
                not sleep.
            sleeper: injected for tests; never called when the probe is off.

        Raises:
            NetworkDeniedError: ``sites`` names a URI scheme (PRIV-05).
            AdapterUnavailableError: no usable shard, a required contig missing or
                incomplete, a filename that disagrees with ``release``/``subset``, a
                directory mixing releases, or two shards claiming the same contig.
        """
        _reject_remote(sites)
        if merge_window_bp < 0:
            msg = "merge_window_bp must not be negative."
            raise AdapterUnavailableError(msg)
        self._release = release.strip()
        self._subset = subset.strip()
        self._merge_window_bp = merge_window_bp
        self._reference = reference
        # Per-run counts over every allele this adapter canonicalised, on either
        # side of the join. Counts, not booleans: the report states how many
        # records are affected, and "one unreadable base" and "the FASTA is gone"
        # are different operator problems — as are "the FASTA is broken" and "this
        # repeat tract is longer than the shift budget", which is why the last two
        # are tallied apart.
        self._unusable_reference_alleles = 0
        self._shift_limited_alleles = 0
        self._indel_alleles = 0
        self._shifted_alleles = 0
        if not self._release or not self._subset:
            msg = "gnomAD adapter requires a non-empty 'release' and 'subset' (GP-18)."
            raise AdapterUnavailableError(msg)

        candidates = self._discover(sites)
        complete: list[Path] = []
        incomplete: list[str] = []
        for candidate in candidates:
            status = check_source_complete(
                candidate,
                stability_probe_seconds=stability_probe_seconds,
                sleeper=sleeper,
            )
            if status.is_complete:
                complete.append(candidate)
            else:
                incomplete.append(f"{candidate.name}: {'; '.join(status.reasons)}")

        self._readers: dict[Path, _SitesReader] = {}
        self._closed = False
        try:
            facts, contig_sources = self._open_all(complete, incomplete)
        except BaseException:
            self.close()
            raise
        self._incomplete_sources: tuple[str, ...] = tuple(sorted(incomplete))

        if facts is None or not contig_sources:
            detail = "; ".join(self._incomplete_sources) or "no *.vcf.bgz files found"
            self.close()
            msg = (
                f"No complete, readable gnomAD sites file under {sites.as_posix()!r} yielded an "
                f"indexed canonical contig. Refusing to open a partially written BGZF stream: it "
                f"decompresses cleanly up to the truncation point, so every variant past it "
                f"would report as absent and read downstream as novel. Details: {detail}"
            )
            raise AdapterUnavailableError(msg)

        self._facts = facts
        self._contig_sources = contig_sources
        self._source = f"{self._facts.dataset}_{self._subset}"

        if require_contigs is not None:
            missing = sorted(
                {normalise_contig(contig) for contig in require_contigs} - set(self._contig_sources)
            )
            if missing:
                detail = "; ".join(self._incomplete_sources) or "none"
                self.close()
                msg = (
                    f"gnomAD release {self._release} has no complete shard for required "
                    f"contig(s) {missing}. Proceeding would return 'no frequency data' for "
                    f"every variant on them, which is scored as absence rather than as a "
                    f"missing resource. Incomplete sources: {detail}"
                )
                raise AdapterUnavailableError(msg)

    # ------------------------------------------------------------------ discovery

    def _discover(self, sites: Path) -> tuple[Path, ...]:
        """Sites files under ``sites``, name-sorted, each checked to be a local file."""
        expanded = sites.expanduser()
        if expanded.is_file():
            return (_local_source_path(expanded),)
        if not expanded.is_dir():
            msg = (
                f"gnomAD sites path {expanded.as_posix()!r} is neither a file nor a directory. "
                "The adapter reads a pre-downloaded release from disk and never fetches."
            )
            raise AdapterUnavailableError(msg)
        found = sorted(
            {
                path
                for suffix in _SOURCE_SUFFIXES
                for path in expanded.glob(f"*{suffix}")
                if path.is_file()
            },
            key=lambda path: path.name,
        )
        return tuple(found)

    def _open_all(
        self, sources: Sequence[Path], incomplete: list[str]
    ) -> tuple[GnomadHeaderFacts | None, dict[str, _ContigSource]]:
        """Validate every shard, then map canonical contig -> shard.

        A shard whose index cannot be opened is added to ``incomplete`` rather than
        raising: the ``.tbi`` may itself be mid-download. That is deliberately not
        silent — the shard is named in :attr:`incomplete_sources`, and
        ``require_contigs`` is how a caller turns a missing contig into a hard
        failure rather than into "gnomAD has nothing there".
        """
        facts: GnomadHeaderFacts | None = None
        contig_sources: dict[str, _ContigSource] = {}

        for candidate in sources:
            path = _local_source_path(candidate)
            self._check_filename(path)
            try:
                raw_contigs = self._indexed_contigs(path)
            except (OSError, ValueError) as exc:
                incomplete.append(f"{path.name}: index unreadable ({type(exc).__name__})")
                continue
            reader = _open_reader(path)
            self._readers[path] = reader
            shard_facts = _parse_header_facts(str(reader.raw_header), path=path)
            if facts is None:
                facts = shard_facts
            elif shard_facts.fingerprint != facts.fingerprint:
                msg = (
                    f"gnomAD sites file {path.name} was produced by a different pipeline run "
                    f"than the other shards in this directory ({shard_facts.fingerprint!r} vs "
                    f"{facts.fingerprint!r}). The header carries no release string, so this "
                    "fingerprint is the only signal that a directory mixes releases; mixing "
                    "them would attribute one release's frequencies to another (GP-18)."
                )
                raise AdapterUnavailableError(msg)

            for raw_contig in raw_contigs:
                canonical = self._canonicalise(raw_contig)
                if canonical is None:
                    continue
                existing = contig_sources.get(canonical)
                if existing is not None:
                    msg = (
                        f"Contig {canonical} is indexed by both {existing.path.name} and "
                        f"{path.name}. Two shards claiming one contig makes the lookup "
                        "order-dependent, so the release layout is refused rather than "
                        "silently resolved (GP-30)."
                    )
                    raise AdapterUnavailableError(msg)
                contig_sources[canonical] = _ContigSource(
                    canonical=canonical, raw=raw_contig, path=path
                )

        return facts, contig_sources

    def _check_filename(self, path: Path) -> None:
        """Cross-check the declared release and subset against the filename.

        The release is matched as a whole ``v``-dotted token rather than as a
        substring: ``"v4" in "gnomad.exomes.v4.1..."`` is true, and so is
        ``"v4.1" in "...v4.10..."``, so a substring test would happily stamp the
        wrong release onto every frequency and every citation (GP-18). A filename
        that carries no version token at all is refused for the same reason — it
        offers no cross-check, and this is the only one available.
        """
        name = path.name
        tokens = set(_FILENAME_RELEASE_RE.findall(name))
        if self._release not in tokens:
            msg = (
                f"gnomAD sites file {name} does not name the declared release "
                f"{self._release!r} as a whole version token (found {sorted(tokens) or 'none'}). "
                "The VCF header carries no release string, so the filename is the only "
                "cross-check available; a mislabelled release would be stamped onto every "
                "frequency and every citation (GP-18)."
            )
            raise AdapterUnavailableError(msg)
        if self._subset.lower() not in name.lower():
            msg = (
                f"gnomAD sites file {name} does not name the declared subset "
                f"{self._subset!r}. Exome and genome callsets have different cohort sizes and "
                "different ancestry coverage; labelling one as the other misstates every "
                "allele number the ranking guards on (ADR 0010)."
            )
            raise AdapterUnavailableError(msg)

    @staticmethod
    def _indexed_contigs(path: Path) -> tuple[str, ...]:
        """Contigs the tabix index actually holds, verbatim.

        Read from the index rather than the header: a gnomAD shard's header declares
        all 24 contigs while the file holds one. Reading the header instead would
        make every shard claim every contig and the first shard opened would answer
        for the whole genome.
        """
        tabix = _open_tabix(path)
        try:
            return tuple(str(contig) for contig in tabix.contigs)
        finally:
            tabix.close()

    @staticmethod
    def _canonicalise(raw_contig: str) -> str | None:
        """UCSC-style contig, or ``None`` for alts/decoys this pipeline cannot reason about."""
        try:
            return normalise_contig(raw_contig)
        except ValueError:
            return None

    # -------------------------------------------------------------------- identity

    @property
    def name(self) -> str:
        return ADAPTER_NAME

    @property
    def version(self) -> str:
        """The gnomAD release, as declared to the constructor and checked against the files.

        Not read from the header, because there is nothing there to read: see the
        class docstring and :attr:`header_fingerprint`, which is the identity the
        *data* carries.
        """
        return self._release

    @property
    def synthetic(self) -> bool:
        """False, deliberately (GP-20).

        ``is_synthetic()`` fails closed, so this line is what opts a real adapter out
        of the "NOT biologically valid" disclosure. It is true here: every number
        this adapter emits was read out of a gnomAD sites VCF.
        """
        return False

    @property
    def source_label(self) -> str:
        """``PopulationFrequency.source`` for every record: dataset from the header + subset."""
        return self._source

    @property
    def header_facts(self) -> GnomadHeaderFacts:
        return self._facts

    @property
    def header_fingerprint(self) -> str:
        """The identity the data itself carries, as opposed to the declared release."""
        return self._facts.fingerprint

    @property
    def build(self) -> GenomeBuild:
        return self._facts.build

    @property
    def populations(self) -> tuple[str, ...]:
        """Populations this adapter can emit: ``global`` plus the header's ancestry groups."""
        return (GLOBAL_POPULATION, *self._facts.ancestry_groups)

    @property
    def available_contigs(self) -> tuple[str, ...]:
        """Canonical contigs backed by a complete, indexed shard, in karyotype order."""
        return tuple(sorted(self._contig_sources, key=contig_sort_key))

    @property
    def contig_map(self) -> Mapping[str, str]:
        """Canonical UCSC contig -> the name this release's index actually uses.

        Exposed so the mapping is assertable directly rather than inferred from a
        successful lookup. A wrong map fails by finding nothing, which reads
        exactly like "gnomAD has no record" — for an entire chromosome.
        """
        return MappingProxyType(
            {canonical: source.raw for canonical, source in self._contig_sources.items()}
        )

    @property
    def incomplete_sources(self) -> tuple[str, ...]:
        """Shards excluded as still-arriving, truncated or unreadable, with the reason."""
        return self._incomplete_sources

    @property
    def left_alignment(self) -> LeftAlignmentReport:
        """The full typed report, derived through the one shared rule.

        Derived, never asserted. This adapter used to label its own batch::

            if self._reference is None:
                return LeftAlignmentStatus.UNAVAILABLE_NO_REFERENCE
            if self._unusable_reference_alleles:
                return LeftAlignmentStatus.INCOMPLETE_REFERENCE_UNUSABLE
            return LeftAlignmentStatus.APPLIED

        which has no case for "there were no indels" and so answered ``APPLIED``
        — rendered by ``describe()`` as "left-alignment applied to all 0 indel
        records" — over a batch that had nothing to left-align.
        :func:`~mva.alleles.summarise_left_alignment` answers ``NOT_REQUIRED``
        there, and the two must not share a value: as
        ``local_tables._left_alignment_report`` puts it, "could not left-align" and
        "had nothing to left-align" are opposite claims about how far to trust the
        rarity of every indel in the run. Two implementations of one status rule is
        the same defect as two implementations of the representation rule
        (ADR 0018), one level up.
        """
        return summarise_left_alignment(
            indel_count=self._indel_alleles,
            shifted_count=self._shifted_alleles,
            unaligned_indel_count=self._unusable_reference_alleles,
            shift_limited_count=self._shift_limited_alleles,
            reference_available=self._reference is not None,
        )

    @property
    def representation_status(self) -> LeftAlignmentStatus:
        """Whether this adapter can reconcile a *shifted* indel spelling, typed.

        ``APPLIED`` only when a reference was supplied, **every read the rule
        needed from it succeeded**, and there was at least one indel to place. A
        reference object that is merely non-``None`` proves nothing: a FASTA that
        raises on every read yields trim-only keys, and reporting those as
        reference-backed is a provenance lie no downstream consumer can detect.
        Trimming is unconditional either way, so the non-minimal class of mismatch
        joins in both states; what the degraded state costs is the repeat-tract
        class, where gnomAD and the caller place one insertion at different
        positions. Typed rather than logged so a report can state the limitation
        instead of a reader having to infer it from a missing key — which is
        precisely the inference GP-14 forbids.

        Read off :attr:`left_alignment` rather than re-derived here.
        """
        return self.left_alignment.status

    @property
    def representation_limitation(self) -> str | None:
        """One sentence for a report footer, or ``None`` when nothing is degraded."""
        status = self.representation_status
        if status in (LeftAlignmentStatus.APPLIED, LeftAlignmentStatus.NOT_REQUIRED):
            return None
        if status is LeftAlignmentStatus.INCOMPLETE_SHIFT_LIMIT:
            return (
                f"The reference FASTA was readable throughout, but "
                f"{self._shift_limited_alleles} gnomAD allele(s) sit in a repeat tract "
                "longer than the shift budget and stopped short of their left-most "
                "position. Their keys are reproducible but not left-most, so they may "
                "fail to join against gnomAD's left-aligned alleles and would then be "
                "reported as having no frequency data — absence of information, not "
                "evidence of rarity. The FASTA is not at fault here."
            )
        if self._reference is not None:
            return (
                f"A reference FASTA was supplied but could not be read for "
                f"{self._unusable_reference_alleles} gnomAD lookup(s), which were "
                "canonicalised by trimming only. Those variants may fail to join "
                "against gnomAD's left-aligned alleles and would then be reported as "
                "having no frequency data — absence of information, not evidence of "
                "rarity. Check that the FASTA is the same assembly and patch release "
                "as the callset, and that its .fai index matches the file."
            )
        return (
            "gnomAD lookups were canonicalised by trimming only: no reference FASTA "
            "was supplied to the adapter, so an indel that gnomAD places at a "
            "different position within the same repeat tract could not be matched. "
            "Such a variant is reported as having no frequency data, which is absence "
            "of information and not evidence of rarity."
        )

    def close(self) -> None:
        """Release every open shard handle. Idempotent."""
        for reader in self._readers.values():
            reader.close()
        self._readers = {}
        self._closed = True

    # ------------------------------------------------------------ canonicalisation

    def canonicalise(self, contig: str, position: int, ref: str, alt: str) -> CanonicalAllele:
        """The single entry point both sides of the join go through.

        Delegates to :func:`mva.alleles.canonicalise_allele` — the same function
        :mod:`mva.ingestion.normalise` and
        :class:`~mva.annotation.clinvar_vcf.ClinvarVcfAdapter` call. There is
        deliberately no trimming or shifting logic in this module: a second
        implementation of the rule is what made a left-aligned proband indel and
        the gnomAD record for the same event fail to join, silently, in the
        highest-weight signal the ranker has (ADR 0018).

        Public so that a test can compare this adapter's representation against
        ingestion's and ClinVar's *directly*, rather than inferring that all three
        agree from a join that happened to succeed. Agreement inferred from a
        passing join is exactly the evidence that was available while they
        disagreed.
        """
        canonical = canonicalise_allele(
            contig=contig,
            position=position,
            ref=ref,
            alt=alt,
            reference=self._reference,
        )
        return canonical

    def _count(self, *, is_indel: bool, left_aligned: bool, status: ReferenceStatus) -> None:
        """Fold one canonicalised **query** allele into the left-alignment counts.

        The query side, not the release side, and the distinction is the same one
        ``local_tables._left_alignment_report`` gets wrong when it counts indels in
        the lookup table. This adapter canonicalises every release record it reads
        as well, and a window queried for one SNV routinely holds release indels —
        counting those would report left-alignment work over records the caller
        never asked about, and would put an SNV-only batch back on ``APPLIED`` by a
        different route.

        A release record's key can only ever answer a query with the same key, so
        an SNV-only batch is unaffected by how the release's indels were
        represented. "This run's indel joins may be wrong" is a statement about the
        indels this run looked up.

        Symbolic alleles are excluded from the indel count: ``<DEL>`` against ``A``
        has different lengths and is not a shiftable indel, and counting it would
        claim left-alignment work over records the rule declines to touch.
        """
        if is_indel:
            self._indel_alleles += 1
            if left_aligned:
                self._shifted_alleles += 1
        if status is ReferenceStatus.UNUSABLE:
            self._unusable_reference_alleles += 1
        elif status is ReferenceStatus.SHIFT_LIMIT_REACHED:
            self._shift_limited_alleles += 1

    # --------------------------------------------------------------------- lookup

    def frequencies(
        self, variant_ids: Sequence[str]
    ) -> Mapping[str, tuple[PopulationFrequency, ...]]:
        """Look up frequencies for canonical variant IDs. Fails closed on a coverage hole.

        A variant with no gnomAD record is **absent from the result**. It is never
        present with ``allele_frequency=0.0`` (GP-14). Callers must handle the
        missing key; ``VariantRecord.has_frequency_data`` exists for exactly that.

        A variant on a contig with **no complete shard** is a different fact, and
        this method raises rather than omitting it. Omitting would render "we could
        not look" as "we looked and found nothing": a common variant on a
        chromosome whose 3 GB shard is still downloading would evade the
        common-frequency down-rank and be scored as rare. That is GP-14's core
        failure mode, and it is a live one — a 250 GB release arrives over hours,
        so a run started mid-download hits this path for most of the genome. Use
        :meth:`lookup_partial` to run against an incomplete release deliberately;
        it returns the gap in the result instead of hiding it.

        Work is grouped by contig and coalesced into as few region queries as the
        merge window allows, so each shard is opened once and read forward, and so
        repeat runs issue a byte-identical query sequence (GP-30). The returned
        mapping preserves the caller's variant order.

        Raises:
            AdapterUnavailableError: a queried contig has no complete shard, or the
                adapter is closed, or the backend failed on a region query.
            GenomeBuildMismatchError: a variant ID names another assembly (GP-11).
            ValueError: a variant ID is not the canonical five-field form, or names
                a contig this pipeline does not reason about.
        """
        result = self.lookup_partial(variant_ids)
        if not result.is_complete:
            msg = (
                f"gnomAD release {self._release} has no complete shard for queried "
                f"contig(s) {list(result.unqueryable_contigs)}, so "
                f"{len(result.unqueryable_variant_ids)} variant(s) could not be looked up "
                "at all. Refusing to return a partial answer: a missing key is read "
                "downstream as 'gnomAD has no record', which scores a common variant as "
                "rare (GP-14). Wait for the acquisition step to finish, pass "
                "require_contigs at construction to fail earlier, or call lookup_partial() "
                "to accept the gap explicitly. Coordinates are not echoed (PRIV-09). "
                f"Incomplete sources: {'; '.join(result.incomplete_sources) or 'none'}"
            )
            raise AdapterUnavailableError(msg)
        return result.frequencies

    def lookup_partial(self, variant_ids: Sequence[str]) -> FrequencyLookup:
        """Answer what can be answered, and **name** what could not be asked.

        The explicit opt-in to running against an incomplete release. Everything
        :meth:`frequencies` guarantees still holds for the contigs that are backed
        by a complete shard; the contigs that are not are reported in
        :attr:`FrequencyLookup.unqueryable_contigs` and their variants in
        :attr:`FrequencyLookup.unqueryable_variant_ids`, so a caller that ignores
        them is making a visible choice rather than an invisible one.
        """
        if self._closed:
            msg = (
                "This gnomAD adapter has been closed; its shard handles are released. "
                "Construct a new adapter rather than reusing a closed one — a lookup "
                "against a closed handle would return nothing, which reads as 'gnomAD "
                "has no record' for every variant in the batch."
            )
            raise AdapterUnavailableError(msg)

        queries = [
            _parse_variant_id(variant_id, build=self.build, reference=self._reference)
            for variant_id in dict.fromkeys(variant_ids)
        ]
        for query in queries:
            self._count(
                is_indel=query.is_indel,
                left_aligned=query.left_aligned,
                status=query.reference_status,
            )
        by_contig: dict[str, list[_Query]] = {}
        for query in queries:
            by_contig.setdefault(query.contig, []).append(query)

        found: dict[str, tuple[PopulationFrequency, ...]] = {}
        missing_contigs: list[str] = []
        for contig in sorted(by_contig, key=contig_sort_key):
            source = self._contig_sources.get(contig)
            if source is None:
                missing_contigs.append(contig)
                continue
            self._lookup_contig(source, by_contig[contig], found)

        unqueryable = frozenset(missing_contigs)
        return FrequencyLookup(
            frequencies={
                query.variant_id: found[query.variant_id]
                for query in queries
                if query.variant_id in found
            },
            unqueryable_contigs=tuple(missing_contigs),
            unqueryable_variant_ids=tuple(
                query.variant_id for query in queries if query.contig in unqueryable
            ),
            incomplete_sources=self._incomplete_sources,
        )

    def _lookup_contig(
        self,
        source: _ContigSource,
        queries: Sequence[_Query],
        found: dict[str, tuple[PopulationFrequency, ...]],
    ) -> None:
        """Resolve every query on one contig with as few region queries as possible."""
        reader = self._readers[source.path]
        # The match key carries the position as well as the alleles. A region query
        # spans several bases and returns neighbouring records; matching on
        # (ref, alt) alone would let a different variant a base away answer for
        # this one.
        #
        # A *list* of variant IDs per key, not one: two different caller spellings
        # of the same event canonicalise to the same key (``chr21:5031991 GAA>GA``
        # and ``GA>G`` are one deletion; with a reference configured, so are two
        # spellings from opposite ends of one repeat tract), and a plain dict would
        # let the second silently evict the first, so one of the two callers would
        # be told gnomAD has no record. Assignment into ``found`` rather than
        # accumulation, so a record returned by two overlapping regions cannot be
        # counted twice.
        #
        # Records are matched on the canonical key, never on their raw POS. The
        # spans below no longer partition the queried positions — a query that was
        # left-aligned reaches out to the right-most spelling of its own event — so
        # one release record can be handed back by two fetches, and a record whose
        # raw POS falls outside the span it answers is the very case being fixed.
        wanted: dict[tuple[int, str, str], list[str]] = {}
        for query in queries:
            wanted.setdefault((query.position, query.ref, query.alt), []).append(query.variant_id)
        spans = merge_query_regions([query.span for query in queries], self._merge_window_bp)
        for start, end in spans:
            # Both the region call and the iteration are inside the guard. cyvcf2
            # puts the region string — a proband coordinate — into the message of
            # anything it raises, and htslib defers work to the first `next()`, so
            # guarding only the construction would let the leak through on
            # iteration instead. Exception text and traceback frames reach
            # terminals, log files, crash reports and agent context (PRIV-09).
            try:
                for record in reader(f"{source.raw}:{start}-{end}"):
                    for allele_key, frequencies in self._record_frequencies(
                        record, contig=source.canonical
                    ):
                        if not frequencies:
                            continue
                        for variant_id in wanted.get(allele_key, ()):
                            found[variant_id] = frequencies
            except Exception as exc:
                # Deliberately broad, and deliberately re-raised with the context
                # SUPPRESSED. `raise ... from exc` would chain the original, whose
                # message and frame are exactly what must not be printed; a bare
                # `raise` inside `except` re-exposes it the same way. The token is
                # a one-way handle so two messages about one region still correlate
                # within a run.
                #
                # The shard's *filename* is withheld too, even though it is public
                # reference data: a gnomAD sites file is named after its contig, so
                # printing it would restate the half of the coordinate the token
                # exists to hide. The actionable diagnostic does not need it —
                # `check_source_complete` over the release directory names the bad
                # shard without any reference to what was being looked up.
                handle = error_token((source.canonical, start, end))
                msg = (
                    f"The gnomAD backend failed on a region query against the "
                    f"{self._release} {self._subset} release ({type(exc).__name__}). The "
                    f"region is a proband coordinate and is not echoed, nor is the shard "
                    f"it names; the correlation handle is <region:{handle}> (PRIV-09). A "
                    "shard may be truncated, its index may not match its bytes, or a file "
                    "may have changed under an open handle. Run check_source_complete over "
                    "the release directory to find which one."
                )
                raise AdapterUnavailableError(msg) from None

    def _record_frequencies(
        self, record: _SitesRecord, *, contig: str
    ) -> Iterable[tuple[tuple[int, str, str], tuple[PopulationFrequency, ...]]]:
        """Yield ``((pos, ref, alt), frequencies)`` for each ALT of one gnomAD record.

        Emitted per ALT rather than per record so that an un-split source would still
        attribute ``Number=A`` values to the right allele instead of to the first one.

        The key is built from the **canonicalised** allele, through
        :meth:`canonicalise`, which is the same call the caller's variant ID went
        through — one rule, not two (ADR 0018). ``contig`` is threaded in only so
        the reference can be read; it is not part of the key, which is scoped to
        one contig's shard by construction.

        On a release that is already left-aligned and minimal — which gnomAD's is —
        this costs nothing: :func:`~mva.alleles.canonicalise_allele` exits its shift
        loop on the first comparison for an allele that cannot move, so no
        reference base is fetched for the overwhelming majority of records.
        """
        position = record.POS
        ref = record.REF.upper()
        filter_status = _filter_status(record.FILTERS)
        info = record.INFO

        for index, raw_alt in enumerate(record.ALT):
            canonical = self.canonicalise(contig, position, ref, str(raw_alt).upper())
            yield (
                (canonical.position, canonical.ref, canonical.alt),
                self._build_frequencies(info, index=index, filter_status=filter_status),
            )

    def _build_frequencies(
        self, info: object, *, index: int, filter_status: str | None
    ) -> tuple[PopulationFrequency, ...]:
        """One :class:`PopulationFrequency` per population that this record reports.

        A population whose ``AF_<grp>`` key is absent is **skipped**, not emitted at
        zero. gnomAD omits the key entirely when the group has no called alleles at
        the site (``AN_<grp>=0``, verified on v4.1 exomes), and that is absence of
        information, not an observation of zero carriers. A group with ``AN>0`` and
        ``AC=0`` *is* an observation of zero carriers and is emitted with
        ``allele_frequency=0.0``, which is the strongest rarity evidence gnomAD can
        give — the distinction this whole adapter exists to preserve (GP-14).
        ``chr21:6086421 G>T`` carries both shapes in one record.

        Order is fixed — global, then ancestry groups A-Z — and comes from
        :meth:`_population_fields` rather than from iterating a set (GP-30).
        """
        records: list[PopulationFrequency] = []
        for population, suffix in self._population_fields():
            allele_frequency = _as_float(_per_allele(_get(info, f"AF{suffix}"), index))
            allele_number = _as_int(_get(info, f"AN{suffix}"))
            if allele_frequency is None or allele_number == 0:
                continue
            records.append(
                PopulationFrequency(
                    # GP-18: source, version and population describe the DATA. `source`
                    # is the dataset name from the header plus the declared subset;
                    # `population` is the header's own ancestry-group token. None of
                    # the three is the adapter's own identity.
                    source=self._source,
                    version=self._release,
                    population=population,
                    allele_frequency=allele_frequency,
                    allele_count=_as_int(_per_allele(_get(info, f"AC{suffix}"), index)),
                    allele_number=allele_number,
                    homozygote_count=_as_int(_per_allele(_get(info, f"nhomalt{suffix}"), index)),
                    filter_status=filter_status,
                )
            )
        return tuple(records)

    def _population_fields(self) -> tuple[tuple[str, str], ...]:
        """``(population label, INFO suffix)`` in a fixed order: global, then groups A-Z."""
        return (
            (GLOBAL_POPULATION, ""),
            *((group, f"_{group}") for group in self._facts.ancestry_groups),
        )

    # ------------------------------------------------------------ context manager

    def __enter__(self) -> GnomadSitesFrequencyAdapter:
        return self

    def __exit__(self, *_exc: object) -> None:
        self.close()
