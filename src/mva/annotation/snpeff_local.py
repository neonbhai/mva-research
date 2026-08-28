"""Molecular consequences from a locally installed SnpEff and an offline database.

This is the **real** consequence adapter: it shells out to SnpEff, which is a
published, independently validated variant-effect predictor, and reads back the
``ANN`` field it writes into the VCF. Nothing in this module predicts anything
itself. Transcript models, splice boundaries, frameshift-versus-inframe at exon
edges, NMD and strand handling are exactly the decades of edge cases that a
hand-rolled caller gets subtly and invisibly wrong, so they are delegated whole.

Why SnpEff rather than Ensembl VEP
----------------------------------
Both are acceptable tools. SnpEff was chosen for reasons of *what can actually be
run offline on this machine today*: it is a single self-contained JAR needing
only a JRE, and its GRCh38 database is ~1 GB rather than VEP's ~26 GB cache.
VEP's official distribution here would require starting Docker Desktop — a GUI
application on the operator's machine — which is not a decision an automated
agent may take. See the module tests for the offline guarantees that follow.

The offline guarantee (PRIV-05)
-------------------------------
This module imports no network client, and the subprocess is pinned offline by
argv rather than by hope. Three flags carry that weight and every one of them is
asserted in :data:`SNPEFF_OFFLINE_FLAGS`:

* ``-noLog`` — SnpEff *by default* posts an anonymous usage record to its own
  server on every run. That is an outbound connection from the annotation stage
  and it is switched off explicitly.
* ``-nodownload`` — SnpEff *by default* downloads a missing genome database over
  the network. Without this flag a typo in the database name turns a local run
  into a silent remote fetch.
* ``-noStats`` — SnpEff otherwise writes a summary HTML/TXT file stamped with the
  wall-clock time of the run, next to the working directory. That is a
  determinism hazard (GP-30) before it is anything else.

The proband's coordinates never leave this process: the input VCF is handed to
SnpEff on **stdin**, so it is never written to a file, never named in a command
line, and never recorded in the ``##SnpEffCmd`` header of the output.

Locating the JVM, explicitly
----------------------------
The Java runtime is named by an absolute path -- a constructor argument, or
``$JAVA_HOME/bin/java`` -- and **never** by the bare name ``java`` resolved
through ``PATH``. That is not defensive style; it is the observed behaviour of
the machine this was built on. macOS ships ``/usr/bin/java`` as a *shim* that is
present, executable, and not a runtime: it prints "Unable to locate a Java
Runtime" and exits non-zero. ``shutil.which("java")`` finds it, ``Path.is_file``
confirms it, and every check short of running it passes. A resolver that fell
back to ``PATH`` would therefore find the shim and fail inside SnpEff's version
probe, where the message would blame SnpEff for a missing JDK.
:func:`resolve_java_binary` refuses that fallback, and :meth:`
SnpEffConsequenceAdapter._verify_jvm` runs the candidate before anything else, so
a non-runtime is diagnosed as a non-runtime at construction time.

What comes back, and what does not
----------------------------------
Every ``ANN`` entry SnpEff emits becomes one :class:`ConsequenceAnnotation`.
Nothing is collapsed to the canonical or MANE-Select transcript: that is a
data-loss bug, because a variant can be benign on MANE-Select and
splice-disrupting on the tissue-relevant isoform. Ordering is presentation only
(see :func:`consequence_sort_key`), never selection.

Two fields SnpEff's ``ANN`` format simply does not carry:

* ``is_canonical`` — not reported. It stays ``False``, which the model cannot
  distinguish from "known not to be canonical". This is a limitation of
  ``ConsequenceAnnotation``, recorded here rather than papered over.
* ``is_mane_select`` — not reported either, but it *is* recoverable from the NCBI
  MANE summary release, a small public reference file. Pass ``mane_summary`` and
  the flag is set from that file's ``MANE Select`` rows; omit it and the flag
  stays ``False`` for every transcript.

One thing SnpEff reports that this adapter currently has nowhere to put: the
per-annotation diagnostic codes in the last ``ANN`` sub-field.
``WARNING_REF_DOES_NOT_MATCH_GENOME`` in particular says the REF allele the
pipeline supplied disagrees with the reference the database was built on, which
means that entry's HGVS and impact were computed from the wrong base. The
annotation is still returned — dropping it would be data loss, and ingestion
raises its own ``ref_allele_mismatch`` flag for the same condition — but nothing
in ``ConsequenceAnnotation`` carries the warning, so a reader of the evidence
trail cannot see it. Closing that gap needs a field on the model
(``tool_warnings: tuple[str, ...]``), which is outside this module.

``splice_ai_delta_max`` and ``pathogenicity_scores`` stay empty: base SnpEff
computes neither. They belong in this adapter (per the adapter README, SpliceAI
is a ``ConsequenceAdapter`` concern and not a fourth Protocol) and would arrive
through a SnpSift/dbNSFP annotation pass, which is not wired here.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import subprocess
from collections.abc import Generator, Iterable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass, replace
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Final, NoReturn, cast

from mva.determinism import hash_file
from mva.errors import AdapterUnavailableError, GenomeBuildMismatchError
from mva.models.base import error_token
from mva.models.genome import ContigStyle, GenomeBuild, contig_sort_key, normalise_contig
from mva.models.variant import ConsequenceAnnotation, ImpactSeverity

__all__ = [
    "JAVA_HOME_ENV",
    "SNPEFF_ADAPTER_NAME",
    "SNPEFF_OFFLINE_FLAGS",
    "AnnEntries",
    "SnpEffArtifactPins",
    "SnpEffConsequenceAdapter",
    "SnpEffRunReport",
    "SnpEffSite",
    "consequence_sort_key",
    "genome_is_declared",
    "load_mane_select_ids",
    "parse_ann_entries",
    "plan_batch",
    "render_input_vcf",
    "resolve_java_binary",
]

#: Stamped onto every EvidenceItem this adapter justifies (GP-20).
SNPEFF_ADAPTER_NAME: Final = "snpeff"

#: Flags without which SnpEff may reach the network or emit a timestamp.
#:
#: Kept as a named constant so the offline guarantee is a value a test can assert
#: against the constructed argv, rather than a comment somebody has to believe.
SNPEFF_OFFLINE_FLAGS: Final[tuple[str, ...]] = ("-noLog", "-nodownload", "-noStats")

#: How long one SnpEff invocation may run before it is killed, in seconds.
DEFAULT_TIMEOUT_SECONDS: Final = 1800.0

#: JVM heap. GRCh38 predictors need well over the 1 GB default.
DEFAULT_JAVA_HEAP: Final = "6g"

#: How long the construction-time JVM and SnpEff probes may run, in seconds.
#: Short: both are expected to answer immediately, and a hang here is a broken
#: install rather than a big workload.
PROBE_TIMEOUT_SECONDS: Final = 120.0

#: The environment variable consulted when no ``java_binary`` is passed. The only
#: fallback there is -- see :func:`resolve_java_binary`.
JAVA_HOME_ENV: Final = "JAVA_HOME"


def resolve_java_binary(explicit: Path | None = None) -> Path:
    """Locate a Java runtime by absolute path, never through ``PATH``.

    Resolution order, and there is deliberately no third step:

    1. ``explicit``, when the caller passed one.
    2. ``$JAVA_HOME/bin/java``.

    A bare ``java`` found on ``PATH`` is **not** accepted, because on macOS that
    name resolves to a system shim which exists, is executable, and is not a
    runtime. Accepting it would convert a clear "no JDK is installed" into a
    SnpEff failure several steps later. Returning a path here is not a claim that
    the file works: :meth:`SnpEffConsequenceAdapter._verify_jvm` runs it.

    Raises:
        AdapterUnavailableError: no candidate, or the candidate is not a file.
    """
    if explicit is not None:
        candidate = explicit.expanduser()
        if not candidate.is_file():
            msg = (
                f"Java runtime not found at {candidate}. Pass java_binary=<path to a real "
                "JRE's bin/java>, or run tools/setup/install_snpeff.sh, which symlinks one "
                "into the install root."
            )
            raise AdapterUnavailableError(msg)
        return candidate.resolve()

    java_home = os.environ.get(JAVA_HOME_ENV, "").strip()
    if not java_home:
        msg = (
            f"No Java runtime was given and {JAVA_HOME_ENV} is unset. SnpEff needs a JVM. "
            "Pass java_binary=<path>/bin/java explicitly, or export "
            f"{JAVA_HOME_ENV}. A bare 'java' on PATH is deliberately NOT used as a "
            "fallback: on macOS that name resolves to a system shim which is present and "
            "executable but is not a runtime, so the fallback would look like it worked "
            "and then fail inside SnpEff."
        )
        raise AdapterUnavailableError(msg)
    candidate = Path(java_home).expanduser() / "bin" / "java"
    if not candidate.is_file():
        msg = (
            f"{JAVA_HOME_ENV} is set but {candidate} is not a file, so it does not point at "
            "a Java runtime. Fix the variable or pass java_binary explicitly; PATH is "
            "deliberately not consulted as a fallback."
        )
        raise AdapterUnavailableError(msg)
    return candidate.resolve()


def genome_is_declared(config_path: Path, genome_database: str) -> bool:
    """Whether ``snpEff.config`` declares ``<genome_database>.genome``.

    SnpEff resolves a genome name through its config file, not through the data
    directory, so a database whose files are present but whose name is undeclared
    fails at run time with a message about downloading it -- which, with
    ``-nodownload`` set, is a confusing way to say "typo". Checking here turns
    that into a construction-time error, offline, before a JVM is started.

    Scanned line by line and abandoned on the first match: the shipped 5.4c config
    is ~19 MB and declares tens of thousands of genomes.
    """
    prefix = f"{genome_database}.genome"
    try:
        with config_path.open("r", encoding="utf-8", errors="replace") as handle:
            return any(line.lstrip().startswith(prefix) for line in handle)
    except OSError:
        return False


#: Variants per SnpEff invocation. The JVM start-up and the GRCh38.115 database
#: load are the entire cost and are paid once per chunk, so this is set high; the
#: reason it is bounded at all is memory, not speed.
#:
#: Measured on this machine: 35.5 s for ONE variant, 31.6 s for 5,000. The
#: marginal per-variant cost is indistinguishable from zero, so the number of
#: launches *is* the runtime. At 4,962,060 variants that is the difference between
#: 20 launches (~25-30 min) and 993 (~9.8 h of reloading the same database).
#:
#: It must equal :data:`mva.annotation.service.DEFAULT_ANNOTATION_BATCH_SIZE`.
#: Whichever is smaller decides the launch count, so a mismatch silently keeps the
#: reloads: 25,000 under a 250,000 service batch means ten launches per batch and
#: 199 of the original 993 survive.
#:
#: The old 25,000 ceiling existed because ``_invoke`` used ``capture_output=True``
#: and held the whole annotated VCF in memory as one ``bytes`` object. It no longer
#: does — stdout is spooled to a file in the run's own scratch directory and read
#: back a line at a time — so the ceiling that number encoded is gone. What still
#: bounds it is the *input* side and the ``VariantRecord`` batch upstream; 250,000
#: is ~918 MB of records against 24 GB of RAM and a 6 GB JVM heap, and 500,000 is
#: deliberately not taken (the scale sweep puts batch-1M at 7.79 GiB against a
#: 12 GiB watchdog, and there is no throughput left to buy).
DEFAULT_BATCH_SIZE: Final = 250_000

#: Bytes of a spooled stderr hashed into the ``<stderr:...>`` correlation handle.
#: A tail rather than the whole file: SnpEff's stderr is unbounded (it reports
#: progress per chromosome), and the handle only has to be stable and one-way.
_STDERR_TOKEN_BYTES: Final = 4096


def _tail_bytes(path: Path, limit: int = _STDERR_TOKEN_BYTES) -> bytes:
    """The last ``limit`` bytes of ``path``, or ``b""`` if it cannot be read.

    Seeks rather than reading the file in, because this runs on the failure path
    against a stream whose size is the child's business. Returns bytes, never text:
    the caller hashes them into a one-way handle, and decoding is both unnecessary
    and a way for undecodable output to raise inside an error handler.
    """
    try:
        with path.open("rb") as handle:
            size = handle.seek(0, os.SEEK_END)
            handle.seek(max(0, size - limit))
            return handle.read()
    except OSError:  # pragma: no cover - the spool file is ours and was just written
        return b""


_VCF_HEADER: Final = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

# ANN sub-field order, fixed by the SnpEff "ANN" specification (VCF annotation
# format v1.0). Named rather than indexed inline so a spec change is one edit.
_ANN_ALLELE: Final = 0
_ANN_EFFECT: Final = 1
_ANN_IMPACT: Final = 2
_ANN_GENE_NAME: Final = 3
_ANN_GENE_ID: Final = 4
_ANN_FEATURE_TYPE: Final = 5
_ANN_FEATURE_ID: Final = 6
_ANN_BIOTYPE: Final = 7
_ANN_RANK: Final = 8
_ANN_HGVS_C: Final = 9
_ANN_HGVS_P: Final = 10
_ANN_AA_POS: Final = 13
_ANN_ERRORS: Final = 15
_ANN_MIN_FIELDS: Final = 16

_IMPACT_BY_TOKEN: Final[Mapping[str, ImpactSeverity]] = {
    "HIGH": ImpactSeverity.HIGH,
    "MODERATE": ImpactSeverity.MODERATE,
    "LOW": ImpactSeverity.LOW,
    "MODIFIER": ImpactSeverity.MODIFIER,
}

_IMPACT_RANK: Final[Mapping[ImpactSeverity, int]] = {
    ImpactSeverity.HIGH: 0,
    ImpactSeverity.MODERATE: 1,
    ImpactSeverity.LOW: 2,
    ImpactSeverity.MODIFIER: 3,
}

#: Sort rank for ``impact=None`` -- "not assessed". Ranked *after* MODIFIER and
#: given its own value rather than folded into MODIFIER's 3, because the two are
#: different claims: MODIFIER is a positive prediction of negligible effect, and
#: None is the absence of any prediction (GP-14). This adapter never emits None --
#: :func:`_impact` refuses an ANN entry whose impact token it cannot recognise --
#: but :func:`consequence_sort_key` is public and orders annotations from any
#: source, including ``local_tables``, which can.
_NOT_ASSESSED_RANK: Final = 4


def _impact_rank(impact: ImpactSeverity | None) -> int:
    """Ordinal for sorting. ``None`` sorts last and never collides with MODIFIER."""
    return _NOT_ASSESSED_RANK if impact is None else _IMPACT_RANK[impact]


#: SnpEff's way of saying "that contig is not in my database".
#:
#: Measured against SnpEff 5.4c rather than assumed: SnpEff strips a leading
#: ``chr`` before looking a chromosome up, so sending ``chr15`` to an
#: Ensembl-named database annotates correctly rather than failing. The explicit
#: contig mapping in :func:`plan_batch` is therefore defence in depth, not the
#: only thing standing between this pipeline and silence — but it stays, because
#: the tolerance is SnpEff's behaviour and not a guarantee of the ANN format, and
#: because a RefSeq-accession database (``NC_000015.10``) is not covered by it.
#:
#: What this constant does catch is a contig the database genuinely lacks, where
#: SnpEff emits a placeholder entry with no gene, no feature and this code. A run
#: in which *every* variant comes back that way is a configuration failure that
#: looks exactly like a clean run with nothing to report, so it is raised.
_CHROMOSOME_NOT_FOUND: Final = "ERROR_CHROMOSOME_NOT_FOUND"

#: SnpEff diagnostics meaning "this variant is not on the map I was given". An
#: entry carrying one of these is a failure notice, not an annotation: there is no
#: transcript, no gene and no consequence behind it, only a placeholder row. It is
#: dropped rather than handed downstream as a typed, confident-looking record.
#:
#: Every OTHER diagnostic SnpEff can attach — notably
#: ``WARNING_REF_DOES_NOT_MATCH_GENOME`` and the incomplete-transcript warnings —
#: is left alone, because those entries do describe a real transcript and dropping
#: them would be the data loss this adapter exists to avoid. They are, however,
#: currently *invisible* downstream: ``ConsequenceAnnotation`` has no field for a
#: tool warning. See the module docstring.
_UNPLACEABLE_CODES: Final[frozenset[str]] = frozenset(
    {_CHROMOSOME_NOT_FOUND, "ERROR_OUT_OF_CHROMOSOME_RANGE"}
)

#: VCF INFO percent-escapes, per the VCF 4.3 specification.
_VCF_UNESCAPE: Final[tuple[tuple[str, str], ...]] = (
    ("%3A", ":"),
    ("%3B", ";"),
    ("%3D", "="),
    ("%2C", ","),
    ("%09", "\t"),
    ("%0D", "\r"),
    ("%0A", "\n"),
    ("%25", "%"),
)

#: Alleles that are legal VCF but are not independently annotatable: the spanning
#: deletion and the missing-allele placeholder. Variants carrying them are omitted
#: from the result, because "SnpEff was not asked" is not "SnpEff found nothing".
_UNANNOTATABLE_ALLELES: Final[frozenset[str]] = frozenset({"*", "."})

_VARIANT_ID_PARTS: Final = 5

#: Fixed VCF column positions in SnpEff's output.
_VCF_ID_COLUMN: Final = 2
_VCF_INFO_COLUMN: Final = 7
_VCF_MIN_COLUMNS: Final = 8


# --------------------------------------------------------------------------- pinning

#: Artifact roles a pins manifest must carry. Required rather than best-effort:
#: a manifest that names three of the four leaves the fourth free to change while
#: the run still reports itself as pinned.
_PIN_REQUIRED_ROLES: Final[tuple[str, ...]] = ("jar", "config", "predictor")

#: Roles that are pinned only when the corresponding artifact is configured.
_PIN_OPTIONAL_ROLES: Final[tuple[str, ...]] = ("mane_summary",)

#: A sha256 as ``sha256sum`` and :func:`mva.determinism.hash_file` render one.
#: Matched whole, so a truncated or upper-cased digest is refused rather than
#: silently compared against something that can never equal it.
_SHA256_RE: Final = re.compile(r"[0-9a-f]{64}")


@dataclass(frozen=True, slots=True)
class SnpEffArtifactPins:
    """Expected sha256 of every installed file whose bytes can change an answer.

    Required, and verified before the first variant is annotated, for a reason
    that is specific to how this tool is distributed. The core archive is
    published as ``snpEff_latest_core.zip``: the name is a promise that the bytes
    behind it *will* change. The database archive is served under a stable name
    too, and SnpEff rebuilds a genome under the same ``GRCh38.115`` label when the
    underlying Ensembl annotation is corrected.

    Without a pin, two runs a month apart can produce different transcripts,
    different HGVS and different impacts while stamping the *identical* provenance
    string onto every ``EvidenceItem`` -- which is precisely the failure that
    stamping a version is supposed to make impossible. A version that cannot
    distinguish two different answers is decoration.

    Four files, because all four are inputs to the answer:

    * ``jar`` -- the predictor code itself.
    * ``config`` -- the genome declarations SnpEff resolves a database name through.
    * ``predictor`` -- ``snpEffectPredictor.bin``, the gene model: every transcript,
      exon boundary and CDS offset the annotation is computed from.
    * ``mane_summary`` -- optional, and pinned only when supplied; it decides which
      transcripts are flagged MANE Select.

    The reference sequence files (``sequence.*.bin``) are deliberately not hashed:
    they are ~550 MB, they are re-read on every construction, and a change to them
    without a change to the predictor is not a case SnpEff's packaging produces.
    That is a stated limit of this pin, not an oversight.

    Two further fields, ``snpeff_release`` and ``genome_database``, are *declared
    identity* rather than artifact bytes. They say which release and which genome
    the reviewed digests above were reviewed *for*, and the adapter checks the
    installation it is actually about to run against them
    (:meth:`SnpEffConsequenceAdapter._verify_pins` and
    :meth:`~SnpEffConsequenceAdapter._probe_version`). Four matching digests
    handed to an adapter configured for a different genome would verify bytes
    that answer a different question, and stamp one database's provenance onto
    the other's answers. They are deliberately **not** part of :meth:`as_rows`
    or :attr:`composite_digest`, which are over bytes.
    """

    jar: str
    config: str
    predictor: str
    mane_summary: str | None = None
    snpeff_release: str | None = None
    """The SnpEff release these digests were reviewed for, e.g. ``5.4c``. Checked
    against the release the JAR reports about itself; ``None`` skips the check."""
    genome_database: str | None = None
    """The genome the ``predictor`` digest belongs to, e.g. ``GRCh38.115``. Checked
    against the adapter's configured database; ``None`` skips the check."""

    @classmethod
    def from_manifest(cls, path: Path, *, genome_database: str | None = None) -> SnpEffArtifactPins:
        """Load the reviewed pins ``tools/setup/install_snpeff.sh`` wrote.

        This is the counterpart to :meth:`measure`, and the difference between
        them is the whole point. ``measure`` hashes whatever is installed, so
        feeding its output back into the constructor pins the installation to
        itself and verifies nothing. The installer writes ``snpeff_pins.json``
        only *after* checking every digest against the reviewed ``EXPECT_*``
        constants in the script and refusing to finish otherwise, so the manifest
        is an assertion about which bytes were reviewed rather than a recording of
        which bytes are present.

        Every field is validated. A manifest missing a role would leave that
        artifact unpinned while the run reported itself as pinned, which is worse
        than no pin at all: a version string that cannot distinguish two different
        answers is decoration.

        Args:
            path: the installer's ``snpeff_pins.json``.
            genome_database: when given, the database the caller intends to run.
                A manifest describing a different genome is refused here rather
                than at the first annotation.

        Raises:
            AdapterUnavailableError: the file is missing, is not a JSON object, is
                missing a required key, carries a digest that is not a sha256, or
                describes a different genome from ``genome_database``.
        """
        if not path.is_file():
            msg = (
                f"SnpEff pins manifest {path.as_posix()!r} not found. Run "
                "tools/setup/install_snpeff.sh, which verifies every artifact against "
                "the reviewed digests in the script and then writes this file. Running "
                "unpinned is refused: see SnpEffArtifactPins."
            )
            raise AdapterUnavailableError(msg)
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            msg = (
                f"SnpEff pins manifest {path.name!r} could not be read as JSON "
                f"({type(exc).__name__}). It is written by tools/setup/install_snpeff.sh; "
                "re-run it rather than hand-editing the file."
            )
            raise AdapterUnavailableError(msg) from exc
        if not isinstance(payload, dict):
            msg = (
                f"SnpEff pins manifest {path.name!r} is not a JSON object, so it names no "
                "artifact digests at all."
            )
            raise AdapterUnavailableError(msg)
        entries = cast("dict[str, object]", payload)

        digests: dict[str, str] = {}
        for role in (*_PIN_REQUIRED_ROLES, *_PIN_OPTIONAL_ROLES):
            raw = entries.get(role)
            if raw is None:
                if role in _PIN_REQUIRED_ROLES:
                    msg = (
                        f"SnpEff pins manifest {path.name!r} carries no {role!r} digest. "
                        "Every artifact that can change an answer must be pinned: a "
                        "manifest that pins three of four would verify the run and leave "
                        "the fourth free to change while the provenance string stayed "
                        "identical."
                    )
                    raise AdapterUnavailableError(msg)
                continue
            digest = str(raw).strip()
            if _SHA256_RE.fullmatch(digest) is None:
                msg = (
                    f"SnpEff pins manifest {path.name!r} gives {role!r} a value that is not "
                    "a sha256 digest (64 lowercase hex characters). Refusing to treat an "
                    "unreadable pin as a pin."
                )
                raise AdapterUnavailableError(msg)
            digests[role] = digest

        declared_database = entries.get("genome_database")
        if declared_database is None or not str(declared_database).strip():
            msg = (
                f"SnpEff pins manifest {path.name!r} names no 'genome_database'. The "
                "predictor digest belongs to one genome and means nothing without it: "
                "the same four digests handed to an adapter configured for another "
                "database would verify successfully and mislabel every annotation."
            )
            raise AdapterUnavailableError(msg)
        database = str(declared_database).strip()
        if genome_database is not None and database != genome_database:
            msg = (
                f"SnpEff pins manifest {path.name!r} describes genome {database!r}, but "
                f"{genome_database!r} was requested. Refusing to verify one genome's "
                "bytes and annotate against another's."
            )
            raise AdapterUnavailableError(msg)

        release = entries.get("snpeff_release")
        return cls(
            jar=digests["jar"],
            config=digests["config"],
            predictor=digests["predictor"],
            mane_summary=digests.get("mane_summary"),
            snpeff_release=str(release).strip() if release is not None else None,
            genome_database=database,
        )

    @classmethod
    def measure(
        cls,
        *,
        jar_path: Path,
        config_path: Path,
        predictor_path: Path,
        mane_summary: Path | None = None,
    ) -> SnpEffArtifactPins:
        """Hash what is on disk **now**, to write into a manifest.

        This is a recording step, not a verification step. Feeding its output
        straight back into the constructor pins the adapter to whatever happens to
        be installed and checks nothing; the digests belong in a reviewed manifest
        entry, and the constructor should be given *those*.
        """
        return cls(
            jar=hash_file(jar_path),
            config=hash_file(config_path),
            predictor=hash_file(predictor_path),
            mane_summary=hash_file(mane_summary) if mane_summary is not None else None,
        )

    def as_rows(self) -> tuple[tuple[str, str], ...]:
        """``(role, digest)`` pairs in a fixed order, for hashing and for display."""
        rows = [("jar", self.jar), ("config", self.config), ("predictor", self.predictor)]
        if self.mane_summary is not None:
            rows.append(("mane_summary", self.mane_summary))
        return tuple(rows)

    @property
    def composite_digest(self) -> str:
        """A single short digest over all pinned artifacts, for the version string.

        Twelve hex characters of a sha256 over ``role=digest`` lines. Short enough
        to sit in a report footer, wide enough that two distinct installations
        colliding is not a practical concern.
        """
        rendered = "\n".join(f"{role}={digest}" for role, digest in self.as_rows())
        return hashlib.sha256(rendered.encode("utf-8")).hexdigest()[:12]


# --------------------------------------------------------------------------- input


@dataclass(frozen=True, slots=True)
class SnpEffSite:
    """One variant, rendered into the contig style the local database expects.

    ``key`` is the handle written into the VCF ``ID`` column and read back out of
    SnpEff's output. Joining on it rather than on the coordinate means the result
    survives any renaming SnpEff does to a chromosome on the way through, which is
    precisely the failure mode that silently produces zero matches.
    """

    variant_id: str
    key: str
    contig: str
    position: int
    ref: str
    alt: str

    def sort_key(self) -> tuple[int, int, str, str]:
        return (contig_sort_key(self.contig), self.position, self.ref, self.alt)


def _unique(variant_ids: Sequence[str]) -> tuple[str, ...]:
    """Deduplicate while preserving caller order."""
    seen: dict[str, None] = {}
    for variant_id in variant_ids:
        seen.setdefault(variant_id, None)
    return tuple(seen)


def plan_batch(
    variant_ids: Sequence[str],
    *,
    build: GenomeBuild,
    contig_style: ContigStyle,
) -> tuple[SnpEffSite, ...]:
    """Turn canonical variant IDs into annotatable sites, in a stable order.

    Sorted by coordinate so that the bytes handed to SnpEff depend only on the
    *set* of variants, not on the order the caller happened to iterate in (GP-30).

    Raises :class:`GenomeBuildMismatchError` when an ID names a different assembly
    from the one the local database was built on. Annotating GRCh37 coordinates
    against a GRCh38 database does not fail — it returns confident nonsense — so
    the mismatch is refused rather than tolerated.
    """
    sites: list[SnpEffSite] = []
    for variant_id in _unique(variant_ids):
        parts = variant_id.split(":", _VARIANT_ID_PARTS - 1)
        if len(parts) != _VARIANT_ID_PARTS:
            msg = (
                f"Malformed variant ID <variant:{error_token(variant_id)}>: expected "
                "'{build}:{contig}:{pos}:{ref}:{alt}'. The value is tokenised rather "
                "than echoed (PRIV-09)."
            )
            raise AdapterUnavailableError(msg)
        raw_build, raw_contig, raw_position, ref, alt = parts
        if raw_build != build.value:
            msg = (
                f"Genome build mismatch: <variant:{error_token(variant_id)}> is "
                f"{raw_build!r} but the SnpEff database is {build.value}. Cross-build "
                "annotation is refused; lift-over must be an explicit, "
                "provenance-tracked stage."
            )
            raise GenomeBuildMismatchError(msg)
        if ref.upper() in _UNANNOTATABLE_ALLELES or alt.upper() in _UNANNOTATABLE_ALLELES:
            continue
        try:
            position = int(raw_position)
            contig = normalise_contig(raw_contig, contig_style)
        except ValueError as exc:
            msg = (
                f"Unusable coordinate in <variant:{error_token(variant_id)}>; the value "
                "is tokenised rather than echoed (PRIV-09)."
            )
            raise AdapterUnavailableError(msg) from exc
        sites.append(
            SnpEffSite(
                variant_id=variant_id,
                key="",
                contig=contig,
                position=position,
                ref=ref.upper(),
                alt=alt.upper(),
            )
        )
    # Keys are assigned *after* sorting, so the rendered VCF -- and therefore every
    # byte the subprocess sees -- depends only on the set of variants, never on the
    # order the caller happened to iterate in (GP-30).
    return tuple(
        replace(site, key=f"mva{index}")
        for index, site in enumerate(sorted(sites, key=SnpEffSite.sort_key))
    )


def render_input_vcf(sites: Iterable[SnpEffSite]) -> str:
    """A minimal, header-stable VCF for SnpEff's stdin.

    Deliberately carries no ``##fileDate`` and no sample columns: a timestamp in
    the input is a timestamp in every downstream hash (GP-30), and a genotype is
    not something SnpEff needs in order to predict a consequence (PRIV-09).
    """
    rows = "".join(
        f"{site.contig}\t{site.position}\t{site.key}\t{site.ref}\t{site.alt}\t.\t.\t.\n"
        for site in sites
    )
    return _VCF_HEADER + rows


# --------------------------------------------------------------------------- output


def _vcf_unescape(value: str) -> str:
    for escape, literal in _VCF_UNESCAPE:
        value = value.replace(escape, literal)
    return value


def _optional(value: str) -> str | None:
    """Empty ANN cells are absence, never an empty string (GP-14)."""
    stripped = value.strip()
    return stripped or None


def _leading_int(value: str) -> int | None:
    """``'412/1863'`` -> ``412``; anything unparseable or non-positive -> ``None``."""
    head = value.split("/", 1)[0].strip()
    if not head.isdigit():
        return None
    parsed = int(head)
    return parsed if parsed >= 1 else None


def _impact(token: str, *, ann_index: int) -> ImpactSeverity:
    impact = _IMPACT_BY_TOKEN.get(token.strip().upper())
    if impact is None:
        msg = (
            f"SnpEff emitted an unrecognised impact class {token.strip()!r} in ANN entry "
            f"{ann_index}. Expected one of {sorted(_IMPACT_BY_TOKEN)}."
        )
        raise AdapterUnavailableError(msg)
    return impact


@dataclass(frozen=True, slots=True)
class AnnEntries:
    """The outcome of parsing one VCF record's ``ANN`` field.

    Richer than a list of annotations because an empty list has four distinct
    causes, and telling them apart is the difference between "SnpEff says this
    variant does nothing" and "SnpEff never saw this variant". Only the caller can
    act on that difference, so the parser reports it rather than flattening it.
    """

    annotations: tuple[ConsequenceAnnotation, ...]
    """Every transcript annotation SnpEff emitted, order as emitted."""

    present: bool
    """Whether an ``ANN=`` field existed at all. ``False`` is anomalous for a
    record SnpEff actually processed: it annotates even an intergenic site."""

    contig_unknown: bool
    """SnpEff reported ``ERROR_CHROMOSOME_NOT_FOUND`` for this record."""

    unplaceable: int
    """Entries dropped because SnpEff could not place the variant on the map."""

    incomplete: int
    """Entries dropped for carrying no effect term, gene or feature to hang an
    annotation on. Distinct from ``unplaceable``: the variant was located, but the
    entry describes nothing that ``ConsequenceAnnotation`` could truthfully hold."""


def parse_ann_entries(
    info: str,
    *,
    tool: str,
    tool_version: str,
    mane_select_ids: frozenset[str] = frozenset(),
) -> AnnEntries:
    """Parse one record's ``ANN=`` INFO field into typed annotations.

    Every ANN entry is retained. Where SnpEff describes a non-transcript feature
    (an intergenic region, a motif) it leaves the biotype column empty and names
    the feature in ``Feature_ID``; those become the ``transcript_id`` and the
    ``Feature_Type`` becomes the ``transcript_biotype``, so the tool's answer is
    preserved verbatim instead of being dropped into a false absence.

    Verified against SnpEff 5.4c / GRCh38.115: every one of the 42 ANN entries the
    two ClinVar variants in the integration test produce carries exactly the 16
    sub-fields the format specifies, including the trailing empty ERRORS cell.
    """
    raw = _extract_ann(info)
    if raw is None:
        return AnnEntries((), present=False, contig_unknown=False, unplaceable=0, incomplete=0)
    annotations: list[ConsequenceAnnotation] = []
    contig_unknown = False
    unplaceable = 0
    incomplete = 0
    for ann_index, entry in enumerate(raw.split(",")):
        fields = entry.split("|")
        if len(fields) < _ANN_MIN_FIELDS:
            msg = (
                f"SnpEff ANN entry {ann_index} has {len(fields)} sub-fields, expected at "
                f"least {_ANN_MIN_FIELDS}. The adapter refuses to guess at a changed "
                "annotation format."
            )
            raise AdapterUnavailableError(msg)
        codes = {code.strip() for code in fields[_ANN_ERRORS].split("&") if code.strip()}
        if _CHROMOSOME_NOT_FOUND in codes:
            contig_unknown = True
        if codes & _UNPLACEABLE_CODES:
            unplaceable += 1
            continue
        terms = tuple(term.strip() for term in fields[_ANN_EFFECT].split("&") if term.strip())
        if not terms:
            incomplete += 1
            continue
        gene_symbol = _optional(_vcf_unescape(fields[_ANN_GENE_NAME]))
        feature_id = _optional(_vcf_unescape(fields[_ANN_FEATURE_ID]))
        if gene_symbol is None or feature_id is None:
            # No feature to hang the annotation on. Silently keeping it would put a
            # blank transcript ID into the evidence trail. Counted, not swallowed:
            # a record whose every entry lands here is reported as such.
            incomplete += 1
            continue
        biotype = _optional(fields[_ANN_BIOTYPE]) or _optional(fields[_ANN_FEATURE_TYPE])
        rank = _optional(fields[_ANN_RANK])
        in_intron = "intron_variant" in terms
        annotations.append(
            ConsequenceAnnotation(
                gene_symbol=gene_symbol,
                gene_id=_optional(_vcf_unescape(fields[_ANN_GENE_ID])),
                transcript_id=feature_id,
                transcript_biotype=biotype if biotype is not None else "protein_coding",
                # SnpEff's ANN format carries no canonical flag; see the module
                # docstring. False here means "not reported", not "not canonical".
                is_canonical=False,
                is_mane_select=_strip_version(feature_id) in mane_select_ids,
                consequence_terms=terms,
                impact=_impact(fields[_ANN_IMPACT], ann_index=ann_index),
                hgvs_c=_optional(_vcf_unescape(fields[_ANN_HGVS_C])),
                hgvs_p=_optional(_vcf_unescape(fields[_ANN_HGVS_P])),
                exon=None if in_intron else rank,
                intron=rank if in_intron else None,
                protein_position=_leading_int(fields[_ANN_AA_POS]),
                # SnpEff reports HGVS.p but not VEP's separate reference/alternate
                # amino-acid pair. Deriving it by re-parsing three-letter codes would
                # be inventing a field, so it stays absent (GP-14).
                amino_acids=None,
                splice_ai_delta_max=None,
                pathogenicity_scores={},
                source_tool=tool,
                source_tool_version=tool_version,
            )
        )
    return AnnEntries(
        tuple(annotations),
        present=True,
        contig_unknown=contig_unknown,
        unplaceable=unplaceable,
        incomplete=incomplete,
    )


def _extract_ann(info: str) -> str | None:
    for field in info.split(";"):
        if field.startswith("ANN="):
            return field[len("ANN=") :]
    return None


def _strip_version(feature_id: str) -> str:
    """``ENST00000456328.2`` -> ``ENST00000456328``.

    SnpEff's Ensembl databases carry unversioned transcript IDs while the MANE
    summary carries versioned ones; matching on the accession is what makes the
    two joinable at all.
    """
    return feature_id.split(".", 1)[0]


def consequence_sort_key(
    annotation: ConsequenceAnnotation,
) -> tuple[int, int, str, str, int, str, str, str, str]:
    """MANE, canonical, transcript, gene, impact, terms, HGVS.c, HGVS.p, rank.

    A *total* order over what SnpEff can emit, which is what makes a repeat run
    byte-identical (GP-30). The first two components mirror
    ``local_tables._consequence_sort_key`` so that swapping adapters does not
    reshuffle a report. All transcripts are kept; this is presentation only.

    The trailing three components exist for totality, not for readability. One
    variant can produce two ANN entries agreeing on transcript, gene, impact and
    effect terms and differing only in the HGVS or the exon rank -- SnpEff emits
    exactly that for a variant spanning an exon boundary. Ties there would leave
    the order decided by ``sorted``'s stability over SnpEff's own emission order,
    which is one more thing a repeat run would have to be trusted to reproduce.
    """
    return (
        0 if annotation.is_mane_select else 1,
        0 if annotation.is_canonical else 1,
        annotation.transcript_id,
        annotation.gene_symbol,
        _impact_rank(annotation.impact),
        ",".join(annotation.consequence_terms),
        annotation.hgvs_c or "",
        annotation.hgvs_p or "",
        annotation.exon or annotation.intron or "",
    )


# --------------------------------------------------------------------------- report


@dataclass(frozen=True, slots=True)
class SnpEffRunReport:
    """What happened to every variant in one ``annotate`` call.

    This type exists because of a specific, verified failure. ``annotate``
    previously skipped any output row it could not interpret -- a truncated
    record, an ID it did not recognise -- and returned the mapping it had built so
    far. A SnpEff run that died halfway therefore returned ``{}`` and *looked
    exactly like a clean run over variants that have no consequence*.

    That is not a lost annotation, it is a deleted variant.
    ``VariantRecord.gene_symbols`` is derived purely from ``consequences``, and
    ``prioritization.pairing`` groups candidates by gene symbol, so a variant with
    no consequence is in no gene, is in no group, and is never considered for a
    compound-heterozygous pair. A truncated subprocess would silently delete
    candidate pairs, possibly including the right one, with no error anywhere.

    So the adapter now fails closed on anything it cannot account for, and every
    variant that legitimately produces no annotation lands in exactly one named
    bucket below. The buckets partition the request; :meth:`assert_partitions`
    checks that, because a bucket that quietly failed to cover a case would
    reintroduce the same silence one level up.
    """

    requested: tuple[str, ...]
    """Unique caller-supplied variant IDs, in caller order."""

    skipped_unannotatable: tuple[str, ...]
    """Never sent: a spanning deletion (``*``) or missing allele (``.``), which is
    legal VCF but is not an independently annotatable variant."""

    annotated: tuple[str, ...]
    """SnpEff returned at least one usable transcript annotation."""

    without_ann: tuple[str, ...]
    """A record came back carrying no ``ANN`` field at all. Anomalous for a
    variant SnpEff placed -- it annotates intergenic sites too -- so this is
    reported rather than treated as 'no consequence'."""

    unplaceable: tuple[str, ...]
    """SnpEff could not place the variant: an unknown contig, or a position past
    the end of one. A naming or build problem, never a biological finding."""

    incomplete: tuple[str, ...]
    """Every ANN entry was dropped for describing no effect, gene or feature that
    ``ConsequenceAnnotation`` could hold without inventing a value."""

    records_returned: int
    """VCF data rows read back from the subprocess."""

    @property
    def unannotated(self) -> tuple[str, ...]:
        """Every requested variant that produced no annotation, for any reason."""
        return self.skipped_unannotatable + self.without_ann + self.unplaceable + self.incomplete

    def assert_partitions(self) -> None:
        """Every requested variant is in exactly one bucket, or this is a bug.

        Raises:
            AdapterUnavailableError: the buckets do not partition the request.
        """
        buckets = (
            self.skipped_unannotatable,
            self.annotated,
            self.without_ann,
            self.unplaceable,
            self.incomplete,
        )
        placed = [variant_id for bucket in buckets for variant_id in bucket]
        if len(placed) != len(set(placed)) or set(placed) != set(self.requested):
            msg = (
                f"SnpEff run accounting is inconsistent: {len(self.requested)} variants "
                f"requested, {len(placed)} placed into outcome buckets "
                f"({len(set(placed))} distinct). This is an adapter bug, and it is raised "
                "rather than ignored because an unaccounted variant is exactly the silent "
                "deletion this report exists to prevent. No coordinate is echoed (PRIV-09)."
            )
            raise AdapterUnavailableError(msg)


# --------------------------------------------------------------------------- MANE


def load_mane_select_ids(summary_path: Path) -> frozenset[str]:
    """Unversioned transcript accessions flagged ``MANE Select`` by NCBI.

    Reads the MANE summary release (plain or gzipped) that
    ``tools/setup/install_snpeff.sh`` downloads. Both the Ensembl and the RefSeq
    accession of each row are returned, so the same set works against a SnpEff
    Ensembl database and a RefSeq one.
    """
    try:
        if summary_path.suffix == ".gz":
            with gzip.open(summary_path, "rt", encoding="utf-8") as handle:
                text = handle.read()
        else:
            text = summary_path.read_text(encoding="utf-8")
    except OSError as exc:
        msg = (
            f"MANE summary {summary_path.name} could not be read. Run "
            "tools/setup/install_snpeff.sh, or construct the adapter without "
            "mane_summary to leave every is_mane_select flag unset."
        )
        raise AdapterUnavailableError(msg) from exc
    lines = text.splitlines()
    if not lines:
        msg = f"MANE summary {summary_path.name} is empty."
        raise AdapterUnavailableError(msg)
    # The release's own header line starts with '#', and so may any commentary
    # above it, so the header is found by content rather than by position.
    header: list[str] = []
    first_data_row = 0
    for index, line in enumerate(lines):
        columns = [column.strip().lstrip("#") for column in line.split("\t")]
        if "MANE_status" in columns:
            header = columns
            first_data_row = index + 1
            break
    try:
        ensembl_at = header.index("Ensembl_nuc")
        refseq_at = header.index("RefSeq_nuc")
        status_at = header.index("MANE_status")
    except ValueError as exc:
        msg = (
            f"MANE summary {summary_path.name} lacks the expected columns "
            "('Ensembl_nuc', 'RefSeq_nuc', 'MANE_status'); refusing to guess at the layout."
        )
        raise AdapterUnavailableError(msg) from exc
    selected: set[str] = set()
    width = max(ensembl_at, refseq_at, status_at) + 1
    for line in lines[first_data_row:]:
        if line.startswith("#"):
            continue
        cells = line.split("\t")
        if len(cells) < width or cells[status_at].strip() != "MANE Select":
            continue
        for column in (ensembl_at, refseq_at):
            accession = cells[column].strip()
            if accession:
                selected.add(_strip_version(accession))
    return frozenset(selected)


# --------------------------------------------------------------------------- adapter


class SnpEffConsequenceAdapter:
    """Consequences from a locally installed SnpEff run against a local database.

    Construct with absolute paths to the SnpEff JAR and the *parent* of the
    database directory (SnpEff's ``data.dir``), plus the database name as SnpEff
    knows it, e.g. ``GRCh38.115``. The JVM is named explicitly too -- by
    ``java_binary``, or by ``$JAVA_HOME`` -- and never resolved through ``PATH``;
    see :func:`resolve_java_binary`.

    Everything is validated at construction, in cheapest-first order: paths, then
    the genome's declaration in ``snpEff.config``, then a JVM probe, then SnpEff's
    own version. A missing JAR, an undeclared genome or a ``java`` that is not a
    runtime therefore fails at wiring time with a message naming that specific
    cause, rather than partway through a patient run with a message about
    something else.
    """

    def __init__(
        self,
        *,
        jar_path: Path,
        data_dir: Path,
        genome_database: str,
        pins: SnpEffArtifactPins | None = None,
        java_binary: Path | None = None,
        config_path: Path | None = None,
        build: GenomeBuild = GenomeBuild.GRCH38,
        contig_style: ContigStyle = ContigStyle.ENSEMBL,
        mane_summary: Path | None = None,
        java_heap: str = DEFAULT_JAVA_HEAP,
        batch_size: int = DEFAULT_BATCH_SIZE,
        timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
    ) -> None:
        self._java_binary = resolve_java_binary(java_binary)
        self._jar_path = jar_path.resolve()
        self._data_dir = data_dir.resolve()
        self._genome_database = genome_database
        self._config_path = (
            config_path.resolve()
            if config_path is not None
            else self._jar_path.parent / "snpEff.config"
        )
        self._build = build
        self._contig_style = contig_style
        self._java_heap = java_heap
        self._timeout_seconds = timeout_seconds
        if batch_size < 1:
            msg = f"batch_size must be at least 1, got {batch_size}."
            raise AdapterUnavailableError(msg)
        self._batch_size = batch_size
        self._mane_summary = mane_summary
        self._last_run_report: SnpEffRunReport | None = None
        self._require_installed()
        self._pins = self._verify_pins(pins)
        self._verify_jvm()
        self._mane_select_ids = (
            load_mane_select_ids(mane_summary) if mane_summary is not None else frozenset()
        )
        self._version = self._probe_version()

    # ------------------------------------------------------------------ identity

    @property
    def name(self) -> str:
        return SNPEFF_ADAPTER_NAME

    @property
    def version(self) -> str:
        """``<release>/<database>+<digest>``, e.g. ``5.4c/GRCh38.115+3f9a1c07b2de``.

        Every part is read off the installed artifacts. The release comes from the
        JAR's own ``-version`` output, the database name from the directory that
        was verified and loaded, and the digest from
        :attr:`SnpEffArtifactPins.composite_digest` over the bytes of the JAR, the
        config and the gene model.

        The digest is the part that makes this a version rather than a label.
        ``5.4c/GRCh38.115`` is not unique: SnpEff rebuilds a genome under an
        unchanged name, and the core archive is published under a name that
        guarantees its bytes change. Two runs that disagree about a transcript
        would otherwise cite identical provenance. With the digest they do not.
        """
        return self._version

    @property
    def pins(self) -> SnpEffArtifactPins:
        """The verified artifact digests behind every answer this adapter gives."""
        return self._pins

    @property
    def synthetic(self) -> bool:
        """False, deliberately (GP-20).

        ``is_synthetic()`` fails closed, so this property is the only thing that
        stops a "SYNTHETIC STAND-IN" disclosure being stamped on every evidence
        item. It is set here because the adapter genuinely executes SnpEff against
        a real gene model, and for no other reason.
        """
        return False

    @property
    def database_name(self) -> str:
        return self._genome_database

    @property
    def data_dir(self) -> Path:
        return self._data_dir

    @property
    def jar_path(self) -> Path:
        return self._jar_path

    @property
    def java_binary(self) -> Path:
        """The JVM this adapter will run, resolved to an absolute path.

        Exposed so a caller -- and a test -- can see *which* runtime was chosen
        rather than infer it from a failure. There is no spelling of this that
        means "whatever `java` PATH finds"; see :func:`resolve_java_binary`.
        """
        return self._java_binary

    @property
    def config_path(self) -> Path:
        return self._config_path

    @property
    def last_run_report(self) -> SnpEffRunReport | None:
        """Accounting for the most recent :meth:`annotate` call, or ``None``.

        Every variant that produced no annotation is in exactly one named bucket
        here. That is what makes an omission a classified outcome rather than an
        inference from an empty result.
        """
        return self._last_run_report

    @property
    def mane_select_count(self) -> int:
        """How many MANE Select accessions were loaded; ``0`` means the flag is unset."""
        return len(self._mane_select_ids)

    # ------------------------------------------------------------------ execution

    def build_argv(self) -> tuple[str, ...]:
        """The exact argv used for annotation. Public so tests can assert on it.

        No shell, no string interpolation into a command line, and every offline
        flag present by construction.
        """
        return (
            str(self._java_binary),
            f"-Xmx{self._java_heap}",
            # A fixed locale and timezone inside the JVM as well as outside it:
            # SnpEff formats numbers and dates, and a comma decimal separator would
            # change the bytes of an otherwise identical run (GP-30).
            "-Duser.language=en",
            "-Duser.country=US",
            "-Duser.timezone=UTC",
            "-Dfile.encoding=UTF-8",
            "-jar",
            str(self._jar_path),
            "ann",
            *SNPEFF_OFFLINE_FLAGS,
            "-config",
            str(self._config_path),
            "-dataDir",
            str(self._data_dir),
            "-i",
            "vcf",
            "-o",
            "vcf",
            self._genome_database,
        )

    def annotate(
        self, variant_ids: Sequence[str]
    ) -> Mapping[str, tuple[ConsequenceAnnotation, ...]]:
        """Annotate canonical variant IDs against the local database.

        Variants SnpEff returns nothing for are **omitted** from the mapping rather
        than present with an empty tuple: this adapter cannot tell "no consequence"
        from "never reached the tool", and it is not entitled to imply the former
        (GP-14). Which of those it was, for each omitted variant, is recorded in
        :attr:`last_run_report` -- omission is classified, never inferred from
        silence.

        **This method additionally fails closed on output it could not read.** It
        is the one the ``ConsequenceAdapter`` Protocol exposes, and therefore the
        only one ``annotation.service`` calls, so a caller here has no way to reach
        the report: a missing key is all it sees, and a missing key is read as "this
        variant has no consequence". A record that came back with no ``ANN`` field
        at all is not that. SnpEff annotates even an intergenic site, so a record
        without one is a record whose output this adapter cannot interpret -- a
        truncated writer, a full disk, a JVM that died between the columns -- and
        returning it as absence deletes the variant from ``gene_symbols``, from
        gene grouping and from pairing, silently. :meth:`annotate_with_report` is
        the deliberate way through and hands the classification back in the result
        rather than leaving it on a property a caller can forget to read.

        The other omission buckets are *interpretable* answers and stay classified
        rather than fatal: ``unplaceable`` is SnpEff explicitly reporting
        ``ERROR_CHROMOSOME_NOT_FOUND`` (escalated separately when it is the whole
        batch, which is a naming or build misconfiguration), ``incomplete`` is an
        entry this repository's model cannot hold, and ``skipped_unannotatable``
        was never sent because the caller supplied ``*`` or ``.``.

        Raises:
            AdapterUnavailableError: SnpEff's output could not be fully accounted
                for (see :meth:`annotate_with_report`), or a returned record
                carried no ``ANN`` field.
        """
        annotations, report = self.annotate_with_report(variant_ids)
        if report.without_ann:
            msg = (
                f"SnpEff {self._genome_database} returned {len(report.without_ann)} of "
                f"{report.records_returned} records with no ANN field at all. SnpEff "
                "annotates every site it can place, including intergenic ones, so a "
                "record without one is output this adapter cannot interpret -- not a "
                "finding that the variant has no consequence. Returning it as an absent "
                "key would remove the variant from gene grouping and therefore from "
                "pairing, deleting candidate pairs with nothing in the run saying so. "
                "Call annotate_with_report() to accept the gap explicitly; it returns the "
                "classification in the result. No record content is echoed (PRIV-09)."
            )
            raise AdapterUnavailableError(msg)
        return annotations

    def annotate_with_report(
        self, variant_ids: Sequence[str]
    ) -> tuple[Mapping[str, tuple[ConsequenceAnnotation, ...]], SnpEffRunReport]:
        """:meth:`annotate`, plus the accounting of what happened to every variant.

        **Fails closed on an output it cannot account for.** SnpEff can exit zero
        having written a truncated final record -- a killed pipe, a full disk, an
        OOM in the JVM's writer thread -- and the previous shape of this method
        skipped any row it could not parse. That turned a half-finished run into
        an empty mapping, which downstream reads as "these variants have no
        consequence", which (via ``VariantRecord.gene_symbols``) removes them from
        gene grouping and therefore from pairing altogether. A silently truncated
        run would delete candidate pairs and report success.

        So every planned site must come back exactly once, in a well-formed row,
        with the ID it was sent with. Four conditions are refused outright:

        * a data row with fewer than the eight mandatory VCF columns;
        * a row whose ID was never sent (the join key is not what we think);
        * a row whose ID appears twice (one variant's answer would overwrite another);
        * a planned site with no row at all (the run stopped early).

        Only *after* the record is accounted for is its ``ANN`` interpreted.
        """
        requested = _unique(variant_ids)
        sites = plan_batch(variant_ids, build=self._build, contig_style=self._contig_style)
        planned_ids = {site.variant_id for site in sites}
        skipped = tuple(v for v in requested if v not in planned_ids)

        annotated: dict[str, tuple[ConsequenceAnnotation, ...]] = {}
        without_ann: list[str] = []
        unplaceable: list[str] = []
        incomplete: list[str] = []
        records_returned = 0

        # Chunked because the rendered input is held in memory: a whole genome in
        # one call is tens of GB of VCF text. The *output* no longer is -- `_run`
        # spools it to a file and hands back a line iterator -- which is what
        # allowed DEFAULT_BATCH_SIZE to rise from 25,000 to 250,000 and the JVM
        # launch count over the 4.9 M-variant callset to fall from 993 to 20.
        # Chunk boundaries fall after the coordinate sort and the key assignment, so
        # they are a function of the variant set alone and do not affect the result
        # (GP-30) -- proven by test_batch_size_does_not_change_the_answer.
        for start in range(0, len(sites), self._batch_size):
            chunk = sites[start : start + self._batch_size]
            by_key = {site.key: site.variant_id for site in chunk}
            seen: set[str] = set()
            with self._run(render_input_vcf(chunk)) as lines:
                for line in lines:
                    if not line or line.startswith("#"):
                        continue
                    columns = line.split("\t")
                    if len(columns) < _VCF_MIN_COLUMNS:
                        self._refuse_output(
                            f"a data row carried {len(columns)} tab-separated columns, fewer "
                            f"than the {_VCF_MIN_COLUMNS} mandatory VCF columns, so the run "
                            "was truncated mid-record"
                        )
                    key = columns[_VCF_ID_COLUMN]
                    variant_id = by_key.get(key)
                    if variant_id is None:
                        self._refuse_output(
                            "a record came back under an ID that was never sent, so the "
                            "output cannot be joined to the input"
                        )
                    if key in seen:
                        self._refuse_output(
                            "a record ID came back twice; one variant's annotations would "
                            "silently overwrite another's"
                        )
                    seen.add(key)
                    records_returned += 1
                    parsed = parse_ann_entries(
                        columns[_VCF_INFO_COLUMN],
                        tool=SNPEFF_ADAPTER_NAME,
                        tool_version=self._version,
                        mane_select_ids=self._mane_select_ids,
                    )
                    if parsed.annotations:
                        annotated[variant_id] = tuple(
                            sorted(parsed.annotations, key=consequence_sort_key)
                        )
                    elif not parsed.present:
                        without_ann.append(variant_id)
                    elif parsed.unplaceable:
                        unplaceable.append(variant_id)
                    else:
                        incomplete.append(variant_id)
            if len(seen) != len(chunk):
                self._refuse_output(
                    f"{len(chunk) - len(seen)} of {len(chunk)} variants in a chunk came "
                    "back with no record at all, so SnpEff stopped before finishing"
                )

        # Counted over records SnpEff actually returned. Every one unplaceable is a
        # naming or build failure that would otherwise look like a clean run with
        # nothing to report.
        if records_returned > 0 and len(unplaceable) == records_returned:
            msg = (
                f"SnpEff database {self._genome_database} recognised none of the "
                f"{records_returned} contigs supplied. The adapter is sending "
                f"{self._contig_style.value}-style names; the database expects the other "
                "convention. Construct the adapter with the matching contig_style. No "
                "coordinate is echoed here (PRIV-09)."
            )
            raise AdapterUnavailableError(msg)

        report = SnpEffRunReport(
            requested=requested,
            skipped_unannotatable=skipped,
            annotated=tuple(v for v in requested if v in annotated),
            without_ann=tuple(without_ann),
            unplaceable=tuple(unplaceable),
            incomplete=tuple(incomplete),
            records_returned=records_returned,
        )
        report.assert_partitions()
        self._last_run_report = report
        # Restore the caller's order over the variants that were annotated, so the
        # mapping's iteration order depends on the request rather than on the
        # coordinate sort used for the subprocess input.
        return (
            {
                variant_id: annotated[variant_id]
                for variant_id in requested
                if variant_id in annotated
            },
            report,
        )

    def _refuse_output(self, what: str) -> NoReturn:
        """Raise on subprocess output that cannot be accounted for (PRIV-09).

        The message describes the *shape* of the failure and never the record: a
        malformed row still contains a proband coordinate, and an exception
        message reaches terminals, logs and agent context.
        """
        msg = (
            f"SnpEff {self._genome_database} exited successfully but its output is "
            f"unusable: {what}. Refusing to return a partial result, because an "
            "annotation missing from this mapping is read downstream as 'this variant "
            "has no consequence' -- which removes it from gene grouping and from "
            "pairing entirely. No record content is echoed (PRIV-09)."
        )
        raise AdapterUnavailableError(msg)

    # ------------------------------------------------------------------ internals

    def _require_installed(self) -> None:
        """Every precondition that can be checked without starting a JVM."""
        for label, path in (
            ("Java runtime", self._java_binary),
            ("SnpEff JAR", self._jar_path),
            ("SnpEff config", self._config_path),
        ):
            if not path.is_file():
                msg = (
                    f"{label} not found at {path}. Run tools/setup/install_snpeff.sh to "
                    "install SnpEff and its GRCh38 database."
                )
                raise AdapterUnavailableError(msg)
        predictor = self._data_dir / self._genome_database / "snpEffectPredictor.bin"
        if not predictor.is_file():
            msg = (
                f"SnpEff database {self._genome_database} is not installed under "
                f"{self._data_dir} (expected {predictor.name}). This adapter never "
                "downloads it: run tools/setup/install_snpeff.sh first."
            )
            raise AdapterUnavailableError(msg)
        if not genome_is_declared(self._config_path, self._genome_database):
            msg = (
                f"SnpEff database {self._genome_database} has its files under "
                f"{self._data_dir} but is not declared in {self._config_path.name} "
                f"(no '{self._genome_database}.genome' entry). SnpEff resolves a genome "
                "through the config, so it would report the database as missing and -- "
                "with -nodownload set -- refuse to fetch it, which reads as a network "
                "problem rather than a naming one."
            )
            raise AdapterUnavailableError(msg)

    @property
    def predictor_path(self) -> Path:
        """``snpEffectPredictor.bin`` for the configured database: the gene model."""
        return self._data_dir / self._genome_database / "snpEffectPredictor.bin"

    def _verify_pins(self, pins: SnpEffArtifactPins | None) -> SnpEffArtifactPins:
        """Check every artifact's bytes before one variant is annotated.

        Fails closed when no pins are supplied. An unpinned SnpEff installation is
        not merely unverified -- because the core archive is published under a
        name that guarantees its bytes change, and because a genome can be rebuilt
        under an unchanged name, it is an installation whose provenance string is
        actively misleading.

        Messages name the file and both digests. A hash of public reference data is
        not patient-derived, and file contents are never echoed (PRIV-09).
        """
        if pins is None:
            msg = (
                "Refusing to run SnpEff without artifact pins. Pass "
                "pins=SnpEffArtifactPins(jar=..., config=..., predictor=..., "
                "mane_summary=...) with digests from the resource manifest. SnpEff's "
                "core is published as 'snpEff_latest_core.zip' -- the name guarantees "
                "the bytes change -- and a genome can be rebuilt under an unchanged "
                "name, so an unpinned install lets two different answers carry the "
                "identical provenance string. Record today's digests with "
                "SnpEffArtifactPins.measure(...) and review them into the manifest."
            )
            raise AdapterUnavailableError(msg)
        if pins.genome_database is not None and pins.genome_database != self._genome_database:
            msg = (
                f"These pins were reviewed for genome {pins.genome_database!r} but this "
                f"adapter is configured for {self._genome_database!r}. Refusing to verify "
                "one genome's bytes and annotate against another's: the predictor digest "
                "belongs to a specific gene model, so four matching digests would report a "
                "verified run while the transcripts, HGVS and impacts came from a "
                "different database than the provenance string names."
            )
            raise AdapterUnavailableError(msg)
        if (pins.mane_summary is None) != (self._mane_summary is None):
            msg = (
                "MANE summary pin and MANE summary path must be supplied together: a "
                "file that decides which transcripts are flagged MANE Select is an "
                "input to the answer, so it is either pinned or not used."
            )
            raise AdapterUnavailableError(msg)
        checks: list[tuple[str, Path, str]] = [
            ("SnpEff JAR", self._jar_path, pins.jar),
            ("SnpEff config", self._config_path, pins.config),
            ("SnpEff gene model", self.predictor_path, pins.predictor),
        ]
        if self._mane_summary is not None and pins.mane_summary is not None:
            checks.append(("MANE summary", self._mane_summary, pins.mane_summary))
        for label, path, expected in checks:
            actual = hash_file(path)
            if actual != expected.strip().lower():
                msg = (
                    f"{label} {path.name!r} failed its sha256 pin: expected "
                    f"{expected.strip().lower()}, found {actual}. Refusing to annotate "
                    "against an artifact that is not the reviewed one -- a changed gene "
                    "model produces different transcripts and different HGVS while "
                    "reporting the same database name."
                )
                raise AdapterUnavailableError(msg)
        return pins

    def _verify_jvm(self) -> None:
        """Prove the named ``java`` is actually a runtime, before SnpEff needs it.

        Run first, and separately from the SnpEff probe, because the failure it
        catches is otherwise misattributed. macOS's ``/usr/bin/java`` shim exists
        and is executable, so every path check passes; it exits non-zero saying
        "Unable to locate a Java Runtime". Without this step the first symptom
        would be "SnpEff version probe failed with exit code 1", pointing the
        operator at SnpEff instead of at a missing JDK.

        ``java -version`` is fed no input and its argv contains only paths this
        adapter owns, so its output cannot contain patient data. The first line is
        therefore quoted back, bounded, because it is the whole diagnosis.
        """
        argv = (str(self._java_binary), "-version")
        try:
            completed = subprocess.run(  # noqa: S603 - fixed argv, no shell, no input
                argv,
                input=b"",
                capture_output=True,
                timeout=PROBE_TIMEOUT_SECONDS,
                check=False,
                shell=False,
            )
        except OSError as exc:
            msg = (
                f"{self._java_binary} could not be executed as a Java runtime "
                f"({type(exc).__name__}). Pass java_binary=<path>/bin/java for a real JRE, "
                f"or set {JAVA_HOME_ENV}; PATH is deliberately not consulted."
            )
            raise AdapterUnavailableError(msg) from exc
        except subprocess.TimeoutExpired as exc:
            msg = (
                f"{self._java_binary} -version did not answer within "
                f"{PROBE_TIMEOUT_SECONDS:.0f}s, so it is not a usable Java runtime."
            )
            raise AdapterUnavailableError(msg) from exc
        if completed.returncode != 0:
            detail = (completed.stderr or completed.stdout).decode("utf-8", errors="replace")
            first_line = detail.strip().splitlines()[0][:200] if detail.strip() else "(no output)"
            msg = (
                f"{self._java_binary} is not a working Java runtime: 'java -version' exited "
                f"{completed.returncode} saying {first_line!r}. On macOS /usr/bin/java is a "
                "stub that is present and executable but is not a JRE, which is why this "
                "adapter never resolves 'java' through PATH. Install a JDK (the setup "
                f"script uses Homebrew openjdk@21) and pass java_binary, or set "
                f"{JAVA_HOME_ENV} to its home. No patient data reaches this probe."
            )
            raise AdapterUnavailableError(msg)

    def _probe_version(self) -> str:
        """Read the SnpEff release out of the JAR, then pair it with the database."""
        argv = (
            str(self._java_binary),
            f"-Xmx{self._java_heap}",
            "-jar",
            str(self._jar_path),
            "-version",
        )
        # SnpEff prints `-version` on **stderr**, not stdout, so both streams are
        # scanned. There is no patient data in either: the probe is fed no input.
        stdout, stderr = self._invoke(
            argv, stdin_text="", what="version probe", timeout=PROBE_TIMEOUT_SECONDS
        )
        for line in (stdout + "\n" + stderr).splitlines():
            fields = [field.strip() for field in line.replace("\t", " ").split() if field.strip()]
            if len(fields) >= 2 and fields[0].lower().startswith("snpeff"):
                release = fields[1]
                declared = self._pins.snpeff_release
                if declared is not None and release != declared:
                    msg = (
                        f"The JAR at {self._jar_path.name!r} reports SnpEff {release}, but "
                        f"its pins were reviewed for {declared}. The digests are over "
                        "bytes and the release is over behaviour; they disagree, so one of "
                        "the two is stale. Refusing to annotate: a run must be able to say "
                        "which predictor produced its answers."
                    )
                    raise AdapterUnavailableError(msg)
                return f"{release}/{self._genome_database}+{self._pins.composite_digest}"
        msg = (
            f"Could not read a release string from {self._jar_path.name}. The adapter "
            "refuses to report a version it invented: `version` is what makes an "
            "annotation reproducible, and 'unknown' is not a version."
        )
        raise AdapterUnavailableError(msg)

    @contextmanager
    def _run(self, vcf_text: str) -> Generator[Iterator[str]]:
        """Annotate one chunk, yielding its output **lines** rather than one string.

        A context manager because the lines are read from a spool file that lives
        in the run's scratch directory: the directory has to outlive the read and
        be removed afterwards, and tying that to a ``with`` block is the only way
        to guarantee it on the error paths too.

        Streaming rather than returning ``str`` is what lifted the batch ceiling.
        The previous shape held the entire annotated VCF in memory twice — once as
        the ``bytes`` ``capture_output=True`` accumulates, once as the ``str`` that
        ``.splitlines()`` was called on — and at ~700 bytes of annotated record
        that is ~350 MB per 250,000 variants before the list of lines. Spooling to
        disk and iterating makes the resident cost a line at a time, so the chunk
        size can be chosen on how many JVM launches the run can afford instead.
        """
        with (
            self._spool(
                self.build_argv(),
                stdin_text=vcf_text,
                what="annotation",
                timeout=self._timeout_seconds,
            ) as (stdout_path, _stderr_path),
            stdout_path.open("r", encoding="utf-8", errors="replace") as handle,
        ):
            # `rstrip("\r\n")` rather than `splitlines()`, which also splits on
            # \x0b, \x1c-\x1e and U+2028. Those cannot appear in a VCF data row,
            # but the two must agree exactly: batch size may not change the answer
            # (GP-30), and the whole-input path reads the same lines.
            yield (line.rstrip("\r\n") for line in handle)

    def _invoke(
        self, argv: tuple[str, ...], *, stdin_text: str, what: str, timeout: float
    ) -> tuple[str, str]:
        """Run SnpEff and return both streams as strings.

        The small-output form, for the version probe. :meth:`_run` is the form
        annotation uses, because an annotated whole-genome chunk is not something
        to hold in a ``str``.
        """
        with self._spool(argv, stdin_text=stdin_text, what=what, timeout=timeout) as (
            stdout_path,
            stderr_path,
        ):
            return (
                stdout_path.read_text(encoding="utf-8", errors="replace"),
                stderr_path.read_text(encoding="utf-8", errors="replace"),
            )

    @contextmanager
    def _spool(
        self, argv: tuple[str, ...], *, stdin_text: str, what: str, timeout: float
    ) -> Generator[tuple[Path, Path]]:
        """Run SnpEff in a scratch directory, with both streams spooled to files.

        The environment is rebuilt rather than inherited: an operator's ``LANG`` or
        ``TZ`` must not be able to change the bytes of an annotation run (GP-30),
        and ``HOME``/``TMPDIR`` are pointed at the scratch directory so nothing the
        JVM decides to write lands anywhere permanent. The spool files go in that
        same directory, so they are inside the run's own workspace and are unlinked
        with it — they carry annotated proband coordinates and must not outlive the
        call or land in a shared temp directory (GP-40).

        Redirecting the child's streams to files rather than pipes also removes the
        deadlock the old shape only avoided by accident: SnpEff echoes offending
        records to stderr, and a child filling a stderr pipe while the parent is
        still writing a 250,000-record VCF to its stdin would block both ends.

        Nothing from the child's stderr, and no coordinate from its stdin, is ever
        put into the exception message: an exception message reaches terminals, log
        files, crash reports and agent context (PRIV-09).
        """
        with TemporaryDirectory(prefix="mva-snpeff-") as scratch:
            scratch_path = Path(scratch)
            stdout_path = scratch_path / "snpeff.stdout"
            stderr_path = scratch_path / "snpeff.stderr"
            env = {
                "PATH": "/usr/bin:/bin",
                "LC_ALL": "C",
                "LANG": "C",
                "TZ": "UTC",
                "HOME": scratch,
                "TMPDIR": scratch,
            }
            try:
                with (
                    stdout_path.open("wb") as out_handle,
                    stderr_path.open("wb") as err_handle,
                ):
                    completed = subprocess.run(  # noqa: S603 - fixed argv, no shell
                        argv,
                        input=stdin_text.encode("utf-8"),
                        stdout=out_handle,
                        stderr=err_handle,
                        cwd=scratch,
                        env=env,
                        timeout=timeout,
                        check=False,
                        shell=False,
                    )
            except FileNotFoundError as exc:
                msg = (
                    f"SnpEff {what} could not start: {self._java_binary} is not an "
                    "executable Java runtime. Install one with "
                    "tools/setup/install_snpeff.sh."
                )
                raise AdapterUnavailableError(msg) from exc
            except subprocess.TimeoutExpired as exc:
                msg = (
                    f"SnpEff {what} exceeded its {timeout:.0f}s timeout and "
                    "was killed. Any partial output was discarded with the scratch "
                    "directory; input is not echoed (PRIV-09)."
                )
                raise AdapterUnavailableError(msg) from exc
            if completed.returncode != 0:
                msg = (
                    f"SnpEff {what} failed with exit code {completed.returncode} "
                    f"(database {self._genome_database}). Diagnostics are withheld: SnpEff "
                    "echoes offending VCF records on failure and an exception message "
                    "reaches terminals, logs and agent context (PRIV-09). Reproduce with a "
                    f"non-patient VCF to see them. stderr handle "
                    f"<stderr:{error_token(_tail_bytes(stderr_path))}>."
                )
                raise AdapterUnavailableError(msg)
            yield stdout_path, stderr_path
