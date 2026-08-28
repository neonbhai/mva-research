"""Variant normalisation: parsimony trimming, left-alignment, REF validation.

Two representations of the same deletion — ``chr1:100 AAT>A`` and
``chr1:98 TAA>T`` — are the same biological event but different join keys. Every
downstream stage in this pipeline joins on ``GenomicCoordinate.variant_id``, so an
un-normalised record silently fails to match ClinVar, gnomAD and the challenge
answer key while looking perfectly reasonable to a reader.

**The rule itself is not here.** It lives in :mod:`mva.alleles`, which is a
foundation module both this stage and :mod:`mva.annotation.clinvar_vcf` import.
That is deliberate and it is the fix for a real defect: this module trimmed, the
ClinVar adapter did not, and equivalent variants therefore failed to join — a
failure that presents as "the database has no record", which the ranker reads as
novel and ultra-rare. Two implementations of one rule is the bug. See the citation
of Tan, Abecasis & Kang (2015) in :mod:`mva.alleles`.

What this module will and will not claim:

* **Trimming needs no reference.** Stripping a shared suffix then a shared prefix
  (always keeping at least one base) is pure string surgery and is always applied.
* **Left-alignment needs a reference.** Shifting an indel leftwards through a
  repeat requires reading the bases to its left. Without a ``ReferenceLookup`` the
  shift is not performed, ``left_align`` is *not* recorded in
  ``normalisation_ops`` — recording an operation that did not happen would make the
  provenance trail actively misleading — and the result carries an explicit,
  typed :class:`~mva.alleles.LeftAlignmentReport` saying the run's indel joins are
  unreliable (GP-14). :func:`open_reference_fasta` is how a caller supplies one.
* **A REF disagreement is flagged, not deleted.** GP-13: a record whose REF does not
  match the reference is usually a build or reference-FASTA problem, not a bad
  variant. It is marked ``ref_allele_mismatch`` so :mod:`mva.ingestion.qc` can raise
  a contradicting evidence item, and it only becomes fatal under
  ``strict_reference=True``.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final, Protocol, cast

from pydantic import ValidationError

from mva.alleles import (
    MAX_SHIFT_BP,
    CanonicalAllele,
    LeftAlignmentReport,
    LeftAlignmentStatus,
    ReferenceLookup,
    ReferenceStatus,
    canonicalise_allele,
    is_sequence_allele,
    summarise_left_alignment,
    trim_parsimoniously,
)
from mva.errors import AdapterUnavailableError, ReferenceMismatchError
from mva.models.base import error_token
from mva.models.genome import GenomicCoordinate
from mva.models.variant import (
    FLAG_REF_ALLELE_MISMATCH,
    OP_LEFT_ALIGN,
    OP_SPLIT_MULTIALLELIC,
    OP_TRIM,
    VariantRecord,
    Zygosity,
)

__all__ = [
    "MAX_SHIFT_BP",
    "OP_LEFT_ALIGN",
    "OP_SPLIT_MULTIALLELIC",
    "OP_TRIM",
    "REF_ALLELE_MISMATCH_FLAG",
    "WARN_NO_REFERENCE",
    "WARN_REFERENCE_LOOKUP_FAILED",
    "WARN_REF_ALLELE_MISMATCH",
    "CanonicalAllele",
    "FastaReference",
    "LeftAlignmentReport",
    "LeftAlignmentStatus",
    "NormalisationResult",
    "ReferenceLookup",
    "ReferenceStatus",
    "canonicalise_allele",
    "normalise_variants",
    "open_reference_fasta",
    "split_multiallelic",
    "trim_and_left_align",
    "trim_and_left_align_with_status",
]

# ---------------------------------------------------------------------------
# Operation and flag vocabulary
# ---------------------------------------------------------------------------

#: ``OP_SPLIT_MULTIALLELIC``, ``OP_TRIM`` and ``OP_LEFT_ALIGN`` are imported from
#: :mod:`mva.models.variant`, which owns the vocabulary this stage shares with
#: prioritisation. GP-01/GP-03 keep those two packages from importing each other,
#: so a literal spelled out in both is a literal that drifts.
#:
#: QC flag raised when the record's REF disagrees with the reference sequence.
#: Consumed by :mod:`mva.ingestion.qc` to emit a CONTRADICTS evidence item.
REF_ALLELE_MISMATCH_FLAG: Final = FLAG_REF_ALLELE_MISMATCH

WARN_NO_REFERENCE: Final = "no_reference_lookup_left_alignment_skipped"
WARN_REF_ALLELE_MISMATCH: Final = "ref_allele_mismatch"
WARN_REFERENCE_LOOKUP_FAILED: Final = "reference_lookup_failed"

#: Bases fetched per FASTA read in :class:`FastaReference`. Left-alignment walks
#: one base at a time; a per-base htslib call for every step of every indel in a
#: whole-genome VCF is the difference between minutes and hours. The block is a
#: pure read cache over immutable reference data, so it cannot change a result.
_FASTA_BLOCK_BP: Final = 4096

#: Blocks retained. Normalisation visits records in coordinate order, so a handful
#: of blocks covers the working set; eviction is oldest-first and therefore
#: deterministic (GP-30).
_FASTA_BLOCK_CACHE: Final = 8


# ---------------------------------------------------------------------------
# Public types
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class NormalisationResult:
    """Normalised records plus an auditable account of what was done to them.

    ``operations_applied`` counts *records carrying each operation in the output*,
    which is what a provenance manifest wants to record ("2 records came from a
    multiallelic split, 3 were trimmed"). It therefore includes operations applied
    upstream — a multiallelic site is necessarily decomposed at parse time, since a
    ``VariantRecord`` cannot hold a comma-separated ALT.

    ``left_alignment`` is the typed degraded state (GP-14). A run with no reference
    FASTA cannot left-align, and the consequence is not cosmetic: every indel whose
    source spells it differently from ClinVar or gnomAD will miss its join and be
    scored as novel. That must be a value a caller receives and writes into a QC
    artifact, not a log line. ``warnings`` carries the same statement in prose for
    the run manifest, but a reader who needs to branch on it should branch on
    ``left_alignment.status``.
    """

    variants: tuple[VariantRecord, ...]
    operations_applied: dict[str, int]
    warnings: tuple[str, ...]
    left_alignment: LeftAlignmentReport


@dataclass(frozen=True, slots=True)
class _Outcome:
    """One record's normalisation, plus what the batch report needs to count."""

    record: VariantRecord
    reference_consulted: bool
    """False when no reference was supplied, when the one supplied could not be
    read or disagreed with this record's REF, or when the REF span read cleanly but
    a base the *shift* needed did not. Such a record was trimmed only."""


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def normalise_variants(
    variants: Sequence[VariantRecord],
    *,
    reference: ReferenceLookup | None = None,
    strict_reference: bool = False,
) -> NormalisationResult:
    """Normalise records to a canonical, joinable representation.

    Order of work per record: resolve any multi-allele genotype, validate REF
    against the reference (when one is supplied), then trim and — only with a
    reference — left-align.

    Args:
        variants: records to normalise. Not mutated; new records are returned.
        reference: optional 1-based inclusive reference accessor, e.g. from
            :func:`open_reference_fasta`. Without it, left-alignment is impossible;
            the result then carries a degraded
            :class:`~mva.alleles.LeftAlignmentReport` and a warning, rather than
            skipping quietly.
        strict_reference: when ``True``, a REF/reference disagreement raises
            :class:`~mva.errors.ReferenceMismatchError` instead of flagging. Default
            ``False`` because a mismatch usually indicts the reference or the build,
            and deleting the record would destroy the evidence of that (GP-13).

    Raises:
        ReferenceMismatchError: only under ``strict_reference=True``.
    """
    warnings: Counter[str] = Counter()
    outcomes: list[_Outcome] = []

    for record in variants:
        for split in split_multiallelic(record):
            outcomes.append(
                _normalise_one(
                    split,
                    reference=reference,
                    strict_reference=strict_reference,
                    warnings=warnings,
                )
            )

    outcomes.sort(key=lambda outcome: outcome.record.sort_key())
    output = [outcome.record for outcome in outcomes]

    indels = [outcome for outcome in outcomes if outcome.record.coordinate.is_indel]
    unaligned = [outcome for outcome in indels if not outcome.reference_consulted]
    left_alignment = summarise_left_alignment(
        indel_count=len(indels),
        shifted_count=sum(
            1 for outcome in indels if OP_LEFT_ALIGN in outcome.record.normalisation_ops
        ),
        unaligned_indel_count=len(unaligned),
        reference_available=reference is not None,
    )
    if reference is None and unaligned:
        warnings[WARN_NO_REFERENCE] += len(unaligned)

    operations: Counter[str] = Counter()
    for record in output:
        for operation in record.normalisation_ops:
            operations[operation] += 1

    # Codes first (sorted, so the tuple is stable), then the prose statement of the
    # degraded state last. A reader scanning `warnings` in the run manifest sees the
    # sentence that tells them what it costs them, not only a token.
    messages = [f"{code} (n={warnings[code]})" for code in sorted(warnings)]
    if left_alignment.is_degraded:
        messages.append(left_alignment.describe())

    return NormalisationResult(
        variants=tuple(output),
        operations_applied={key: operations[key] for key in sorted(operations)},
        warnings=tuple(messages),
        left_alignment=left_alignment,
    )


def _normalise_one(
    record: VariantRecord,
    *,
    reference: ReferenceLookup | None,
    strict_reference: bool,
    warnings: Counter[str],
) -> _Outcome:
    current = record
    aligner = reference
    if reference is not None:
        matches = _reference_matches(current, reference)
        if matches is None:
            warnings[WARN_REFERENCE_LOOKUP_FAILED] += 1
            aligner = None
        elif not matches:
            if strict_reference:
                msg = (
                    "REF allele disagrees with the supplied reference sequence for at "
                    "least one record and strict_reference=True. This normally means "
                    "the VCF and the reference FASTA are different assemblies or "
                    "different patch releases; confirm the build before continuing, or "
                    "run with strict_reference=False to flag the record instead of "
                    "refusing the file."
                )
                raise ReferenceMismatchError(msg)
            warnings[WARN_REF_ALLELE_MISMATCH] += 1
            current = current.with_qc_flags(REF_ALLELE_MISMATCH_FLAG)
            # Shifting an indel against a reference that already disagrees would
            # produce a confidently wrong coordinate, so alignment is withheld.
            aligner = None
    normalised, status = trim_and_left_align_with_status(current, aligner)
    if status is ReferenceStatus.UNUSABLE:
        # The REF span read cleanly but a base the shift needed did not. Without
        # this the batch reports APPLIED for a record that was trimmed only.
        warnings[WARN_REFERENCE_LOOKUP_FAILED] += 1
    return _Outcome(
        record=normalised,
        reference_consulted=aligner is not None and status is not ReferenceStatus.UNUSABLE,
    )


# ---------------------------------------------------------------------------
# Multiallelic resolution
# ---------------------------------------------------------------------------


def split_multiallelic(record: VariantRecord) -> tuple[VariantRecord, ...]:
    """Resolve a record whose GT still describes a multi-ALT *site*.

    A ``VariantRecord`` is biallelic by construction — ``GenomicCoordinate.alt``
    accepts one allele and rejects a comma — so the decomposition of ``A>G,GT`` into
    two records happens where the full ALT list still exists, at parse time in
    :func:`mva.ingestion.reader.read_vcf`. That is the only place it *can* happen
    without inventing a second, parallel variant representation.

    This function is the record-level half of the same contract, and it is what
    makes the stage safe to re-run: given a record built from a multiallelic line
    whose genotype was never resolved to its own allele (a hand-built record, or one
    from an adapter that split coordinates but not calls), it re-derives the
    genotype for this record's own allele and marks the operation. An
    already-decomposed or genuinely biallelic record is returned untouched, so the
    function is idempotent.

    Returns a tuple rather than a bare record so callers treat "one in, N out"
    uniformly.
    """
    if OP_SPLIT_MULTIALLELIC in record.normalisation_ops:
        return (record,)

    alleles = _genotype_allele_indices(record.genotype.genotype_string)
    distinct_alts = {allele for allele in alleles if allele is not None and allele > 0}
    if len(distinct_alts) <= 1:
        return (record,)

    # The GT names several ALT alleles; this record stands for exactly one of them,
    # so relative to its own allele the call is heterozygous.
    genotype = record.genotype.model_copy(update={"zygosity": Zygosity.HET})
    updated = record.model_copy(
        update={
            "genotype": genotype,
            "normalisation_ops": _with_ops(record.normalisation_ops, OP_SPLIT_MULTIALLELIC),
        }
    )
    return (updated,)


def _genotype_allele_indices(genotype_string: str) -> tuple[int | None, ...]:
    tokens = genotype_string.strip().replace("|", "/").split("/")
    indices: list[int | None] = []
    for token in tokens:
        try:
            indices.append(int(token))
        except ValueError:
            indices.append(None)
    return tuple(indices)


# ---------------------------------------------------------------------------
# Trimming and left-alignment
# ---------------------------------------------------------------------------


def trim_and_left_align(record: VariantRecord, reference: ReferenceLookup | None) -> VariantRecord:
    """Reduce a record to its most parsimonious, left-most representation.

    A thin record-level wrapper over :func:`mva.alleles.canonicalise_allele`, which
    is the single definition of the rule. The wrapper's only jobs are to unwrap the
    coordinate, to re-validate the result through ``GenomicCoordinate``, and to
    append the operations the shared rule reports to ``normalisation_ops`` — so a
    POS that moved always leaves a trace (GP-31).

    Trimming is always applied. Left-alignment is applied only when ``reference`` is
    not ``None``; when it is ``None`` the returned record records ``trim`` alone and
    never ``left_align``, because claiming an alignment that was never computed
    would be a provenance lie that no downstream consumer could detect.

    The operation is idempotent: normalising an already-normalised record returns an
    equal record with the same ``normalisation_ops``.
    """
    return trim_and_left_align_with_status(record, reference)[0]


def trim_and_left_align_with_status(
    record: VariantRecord, reference: ReferenceLookup | None
) -> tuple[VariantRecord, ReferenceStatus]:
    """:func:`trim_and_left_align`, plus how far the reference could be trusted.

    Separate from the record-only form because ``VariantRecord`` has nowhere to
    put the status: ``normalisation_ops`` records what happened, and a record that
    did not move is ambiguous between "already left-most" and "the reference could
    not be read". The batch report needs the second, so it has to be returned.
    """
    coordinate = record.coordinate
    canonical = canonicalise_allele(
        contig=coordinate.contig,
        position=coordinate.position,
        ref=coordinate.ref,
        alt=coordinate.alt,
        reference=reference,
    )
    return _apply_canonical(record, canonical), canonical.reference_status


def _apply_canonical(record: VariantRecord, canonical: CanonicalAllele) -> VariantRecord:
    """Write a :class:`~mva.alleles.CanonicalAllele` back onto its record.

    The single implementation shared by both public spellings above: ADR 0018
    forbids a second copy of the write-back, for the same reason it forbids a
    second copy of the rule itself.
    """
    if not canonical.changed:
        return record

    coordinate = record.coordinate
    operations = record.normalisation_ops
    for operation in canonical.operations:
        operations = _with_ops(operations, operation)

    try:
        updated = GenomicCoordinate(
            build=coordinate.build,
            contig=coordinate.contig,
            position=canonical.position,
            ref=canonical.ref,
            alt=canonical.alt,
        )
    except (ValidationError, ValueError):  # pragma: no cover - defensive
        return record
    return record.model_copy(update={"coordinate": updated, "normalisation_ops": operations})


#: Retained spelling of the shared trimming rule. Kept because it is imported by
#: name elsewhere in the test suite; the implementation is
#: :func:`mva.alleles.trim_parsimoniously` and there is no second copy of it.
_parsimony_trim = trim_parsimoniously


def _reference_matches(record: VariantRecord, reference: ReferenceLookup) -> bool | None:
    """``True``/``False`` for match/mismatch, ``None`` when the lookup was unusable.

    ``None`` is distinct from ``False`` on purpose: a lookup that failed is absence
    of information, not evidence of a mismatch (GP-14).
    """
    coordinate = record.coordinate
    if not is_sequence_allele(coordinate.ref):
        return None
    try:
        sequence = reference.fetch(coordinate.contig, coordinate.position, coordinate.end)
    except Exception:  # a caller-supplied lookup may raise anything
        return None
    observed = sequence.strip().upper()
    if len(observed) != len(coordinate.ref):
        return None
    return observed == coordinate.ref


def _with_ops(operations: tuple[str, ...], operation: str) -> tuple[str, ...]:
    """Append an operation, order-stable and never duplicated."""
    if operation in operations:
        return operations
    return (*operations, operation)


# ---------------------------------------------------------------------------
# The reference FASTA
# ---------------------------------------------------------------------------


class _FastaHandle(Protocol):
    """Exactly the pysam surface :class:`FastaReference` uses.

    pysam ships no ``py.typed`` marker, so everything it returns is Unknown to
    pyright. Narrowing it to a Protocol at the single point of construction keeps
    the rest of the module typed against a contract written down here.
    """

    @property
    def references(self) -> Sequence[str]: ...

    def fetch(self, reference: str, start: int, end: int) -> str: ...

    def close(self) -> None: ...


class FastaReference:
    """A :class:`~mva.alleles.ReferenceLookup` over an indexed reference FASTA.

    Two boundary translations happen here, once each, so no other module has to
    know about them (GP-02):

    * **Coordinates.** VCF and this pipeline are 1-based inclusive; pysam's
      ``fetch`` is 0-based half-open. Getting this wrong shifts every left-aligned
      indel by one base, which is a silent, systematic, whole-run error.
    * **Contig names.** The GRCh38 analysis-set FASTA is ``chr``-prefixed and the
      proband VCF from GATK is not, while ``GenomicCoordinate`` always holds the
      UCSC spelling. The map from UCSC name to the name this FASTA actually uses is
      resolved once, from the FASTA's own index, and exposed as :attr:`contig_map`
      so a test can assert it instead of inferring it from a miss. An unmapped
      contig raises ``KeyError``, which the shift loop treats as "cannot read here"
      and stops — never as "the bases matched".

    Sequence is upper-cased on the way out: reference FASTAs may be soft-masked,
    and a lower-case ``a`` compared against a VCF's ``A`` is a REF mismatch that
    exists only in the string.
    """

    def __init__(self, handle: _FastaHandle, *, path: Path) -> None:
        self._handle = handle
        self._path = path
        self._contig_map = _resolve_fasta_contigs(handle.references)
        self._blocks: dict[tuple[str, int], str] = {}

    @property
    def path(self) -> Path:
        return self._path

    @property
    def contig_map(self) -> dict[str, str]:
        """UCSC contig -> the name this FASTA uses, e.g. ``chr1 -> 1``."""
        return dict(self._contig_map)

    def fetch(self, contig: str, start: int, end: int) -> str:
        """Reference bases in the 1-based inclusive range ``[start, end]``."""
        if start < 1 or end < start:
            msg = (
                "Reference fetch received an empty or negative 1-based inclusive range. "
                f"Range handle <range:{error_token((contig, start, end))}>; the "
                "coordinate is tokenised rather than echoed (PRIV-09)."
            )
            raise ValueError(msg)
        name = self._contig_map.get(contig)
        if name is None:
            msg = (
                f"Contig <contig:{error_token(contig)}> is not present in reference FASTA "
                f"{self._path.name!r}. The name is tokenised rather than echoed (PRIV-09)."
            )
            raise KeyError(msg)

        first = (start - 1) // _FASTA_BLOCK_BP
        last = (end - 1) // _FASTA_BLOCK_BP
        if last - first > 1:
            # Longer than the cache is for (a large REF span); read it directly.
            return self._handle.fetch(name, start - 1, end).upper()
        assembled = "".join(self._block(name, index) for index in range(first, last + 1))
        offset = start - 1 - first * _FASTA_BLOCK_BP
        return assembled[offset : offset + (end - start + 1)]

    def close(self) -> None:
        """Release the htslib handle."""
        self._handle.close()

    def _block(self, name: str, index: int) -> str:
        cached = self._blocks.get((name, index))
        if cached is not None:
            return cached
        block = self._handle.fetch(
            name, index * _FASTA_BLOCK_BP, (index + 1) * _FASTA_BLOCK_BP
        ).upper()
        if len(self._blocks) >= _FASTA_BLOCK_CACHE:
            # Oldest-first eviction: dicts preserve insertion order, so which entry
            # goes is a function of the input order alone (GP-30).
            del self._blocks[next(iter(self._blocks))]
        self._blocks[(name, index)] = block
        return block


def _resolve_fasta_contigs(references: Sequence[str]) -> dict[str, str]:
    """Map each UCSC contig this pipeline uses to the FASTA's own spelling.

    Resolved from the index rather than assumed. The failure mode of guessing is
    silence: ``fetch("chr1", ...)`` against a bare-contig FASTA raises, left-shift
    stops on the first base, and every indel keeps its input representation while
    the run reports that left-alignment was applied.
    """
    available = set(references)
    resolved: dict[str, str] = {}
    for index in range(1, 23):
        _map_contig(f"chr{index}", (f"chr{index}", str(index)), available, resolved)
    _map_contig("chrX", ("chrX", "X"), available, resolved)
    _map_contig("chrY", ("chrY", "Y"), available, resolved)
    _map_contig("chrM", ("chrM", "chrMT", "MT", "M"), available, resolved)
    return resolved


def _map_contig(
    ucsc: str, candidates: Sequence[str], available: set[str], resolved: dict[str, str]
) -> None:
    for candidate in candidates:
        if candidate in available:
            resolved[ucsc] = candidate
            return


def open_reference_fasta(path: Path) -> FastaReference:
    """Open an indexed reference FASTA as a :class:`~mva.alleles.ReferenceLookup`.

    Creates the ``.fai`` index if it is missing, because there is no ``samtools``
    binary on the target machine and an un-indexed 3 GB FASTA would otherwise be
    read linearly for every base.

    pysam is imported inside the function so that importing this module — which the
    composition root and the architecture tests both do — does not require the
    optional ``genomics`` extra.

    Raises:
        AdapterUnavailableError: the file is missing, pysam is not installed, or
            htslib cannot index or open it. Fails loudly: silently continuing
            without a reference is exactly the degraded run this function exists to
            prevent, and the caller must choose it explicitly by passing ``None``.
    """
    if not path.is_file():
        msg = (
            f"Reference FASTA {path.as_posix()!r} not found. Left-alignment needs the "
            "reference bases to the left of an indel; without them indel join keys "
            "cannot be made to agree with ClinVar or gnomAD."
        )
        raise AdapterUnavailableError(msg)
    try:
        import pysam  # noqa: PLC0415 - native backend, imported on demand
    except ImportError as exc:  # pragma: no cover - guarded by the genomics extra
        msg = (
            "The pysam backend is not installed; install the 'genomics' extra. It is "
            "what reads the reference FASTA that left-alignment requires."
        )
        raise AdapterUnavailableError(msg) from exc

    index_path = path.with_name(path.name + ".fai")
    if not index_path.is_file():
        try:
            pysam.faidx(str(path))
        except Exception as exc:  # htslib raises its own exception types
            msg = (
                f"Could not build a .fai index for reference FASTA {path.name!r}. The file "
                "may be truncated, still downloading, or not bgzip/plain FASTA."
            )
            raise AdapterUnavailableError(msg) from exc
    try:
        return FastaReference(cast(_FastaHandle, pysam.FastaFile(str(path))), path=path)
    except (OSError, ValueError) as exc:
        msg = (
            f"htslib could not open reference FASTA {path.name!r} with its index; the "
            "file and the index may describe different bytes."
        )
        raise AdapterUnavailableError(msg) from exc
