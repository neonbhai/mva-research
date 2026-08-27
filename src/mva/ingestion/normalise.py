"""Variant normalisation: parsimony trimming, left-alignment, REF validation.

Two representations of the same deletion — ``chr1:100 AAT>A`` and
``chr1:98 TAA>T`` — are the same biological event but different join keys. Every
downstream stage in this pipeline joins on ``GenomicCoordinate.variant_id``, so an
un-normalised record silently fails to match ClinVar, gnomAD and the challenge
answer key while looking perfectly reasonable to a reader.

What this module will and will not claim:

* **Trimming needs no reference.** Stripping a shared suffix then a shared prefix
  (always keeping at least one base) is pure string surgery and is always applied.
* **Left-alignment needs a reference.** Shifting an indel leftwards through a
  repeat requires reading the bases to its left. Without a ``ReferenceLookup`` the
  shift is simply not performed, and — critically — ``left_align`` is *not* recorded
  in ``normalisation_ops``. Recording an operation that did not happen would make
  the provenance trail actively misleading.
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
from typing import Final, Protocol

from pydantic import ValidationError

from mva.errors import ReferenceMismatchError
from mva.models.genome import GenomicCoordinate
from mva.models.variant import (
    FLAG_REF_ALLELE_MISMATCH,
    OP_LEFT_ALIGN,
    OP_SPLIT_MULTIALLELIC,
    OP_TRIM,
    VariantRecord,
    Zygosity,
)

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

_NUCLEOTIDES: Final[frozenset[str]] = frozenset("ACGTN")

#: Ceiling on left-shift iterations. A pathological or wrapping ``ReferenceLookup``
#: must not be able to spin this loop forever; 1 kb is far beyond any indel this
#: pipeline is entitled to reason about.
_MAX_LEFT_SHIFT: Final = 1000


# ---------------------------------------------------------------------------
# Public types
# ---------------------------------------------------------------------------


class ReferenceLookup(Protocol):
    """Minimal reference-sequence accessor.

    ``start`` and ``end`` are **1-based and inclusive**, matching VCF convention
    rather than pysam's 0-based half-open ``fetch``. The adapter that wraps a FASTA
    is responsible for that translation, once, at the boundary (GP-02).
    """

    def fetch(self, contig: str, start: int, end: int) -> str:
        """Return the reference bases in ``[start, end]`` on ``contig``."""
        ...


@dataclass(frozen=True, slots=True)
class NormalisationResult:
    """Normalised records plus an auditable account of what was done to them.

    ``operations_applied`` counts *records carrying each operation in the output*,
    which is what a provenance manifest wants to record ("2 records came from a
    multiallelic split, 3 were trimmed"). It therefore includes operations applied
    upstream — a multiallelic site is necessarily decomposed at parse time, since a
    ``VariantRecord`` cannot hold a comma-separated ALT.
    """

    variants: tuple[VariantRecord, ...]
    operations_applied: dict[str, int]
    warnings: tuple[str, ...]


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
        reference: optional 1-based inclusive reference accessor. Without it,
            left-alignment is skipped and a warning is emitted.
        strict_reference: when ``True``, a REF/reference disagreement raises
            :class:`~mva.errors.ReferenceMismatchError` instead of flagging. Default
            ``False`` because a mismatch usually indicts the reference or the build,
            and deleting the record would destroy the evidence of that (GP-13).

    Raises:
        ReferenceMismatchError: only under ``strict_reference=True``.
    """
    warnings: Counter[str] = Counter()
    output: list[VariantRecord] = []

    for record in variants:
        for split in split_multiallelic(record):
            output.append(
                _normalise_one(
                    split,
                    reference=reference,
                    strict_reference=strict_reference,
                    warnings=warnings,
                )
            )

    output.sort(key=lambda record: record.sort_key())

    if reference is None:
        unaligned = sum(1 for record in output if record.coordinate.is_indel)
        if unaligned:
            warnings[WARN_NO_REFERENCE] += unaligned

    operations: Counter[str] = Counter()
    for record in output:
        for operation in record.normalisation_ops:
            operations[operation] += 1

    return NormalisationResult(
        variants=tuple(output),
        operations_applied={key: operations[key] for key in sorted(operations)},
        warnings=tuple(f"{code} (n={warnings[code]})" for code in sorted(warnings)),
    )


def _normalise_one(
    record: VariantRecord,
    *,
    reference: ReferenceLookup | None,
    strict_reference: bool,
    warnings: Counter[str],
) -> VariantRecord:
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
    return trim_and_left_align(current, aligner)


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

    Trimming is always applied. Left-alignment is applied only when ``reference`` is
    not ``None``; when it is ``None`` the returned record records ``trim`` alone and
    never ``left_align``, because claiming an alignment that was never computed
    would be a provenance lie that no downstream consumer could detect.

    The operation is idempotent: normalising an already-normalised record returns an
    equal record with the same ``normalisation_ops``.
    """
    coordinate = record.coordinate
    ref, alt = coordinate.ref, coordinate.alt
    if not _is_simple_allele(ref) or not _is_simple_allele(alt):
        return record

    position, trimmed_ref, trimmed_alt = _parsimony_trim(coordinate.position, ref, alt)
    trimmed = (position, trimmed_ref, trimmed_alt) != (coordinate.position, ref, alt)

    shifted = False
    if reference is not None:
        shift_position, shift_ref, shift_alt = _left_shift(
            coordinate.contig, position, trimmed_ref, trimmed_alt, reference
        )
        shifted = (shift_position, shift_ref, shift_alt) != (position, trimmed_ref, trimmed_alt)
        position, trimmed_ref, trimmed_alt = _parsimony_trim(shift_position, shift_ref, shift_alt)

    if not trimmed and not shifted:
        return record

    operations = record.normalisation_ops
    if trimmed:
        operations = _with_ops(operations, OP_TRIM)
    if shifted:
        operations = _with_ops(operations, OP_LEFT_ALIGN)

    try:
        updated = GenomicCoordinate(
            build=coordinate.build,
            contig=coordinate.contig,
            position=position,
            ref=trimmed_ref,
            alt=trimmed_alt,
        )
    except (ValidationError, ValueError):  # pragma: no cover - defensive
        return record
    return record.model_copy(update={"coordinate": updated, "normalisation_ops": operations})


def _parsimony_trim(position: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Strip the shared suffix, then the shared prefix, always keeping >= 1 base.

    This is the VCF convention: an indel retains one anchoring base, so ``G>GAT`` is
    already minimal while ``AAT>AAG`` reduces to ``T>G`` three bases to the right.
    """
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]

    limit = min(len(ref), len(alt)) - 1
    shared = 0
    while shared < limit and ref[shared] == alt[shared]:
        shared += 1
    if shared:
        ref, alt = ref[shared:], alt[shared:]
        position += shared
    return position, ref, alt


def _left_shift(
    contig: str, position: int, ref: str, alt: str, reference: ReferenceLookup
) -> tuple[int, str, str]:
    """Roll an indel as far left as the reference allows (the vt/bcftools loop)."""
    for _ in range(_MAX_LEFT_SHIFT):
        while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        if len(ref) > 1 and len(alt) > 1:
            break  # a complex substitution, not a shiftable indel
        if ref[-1] != alt[-1]:
            break  # already left-most
        if position <= 1:
            break
        base = _fetch_base(reference, contig, position - 1)
        if base is None:
            break
        ref, alt = base + ref, base + alt
        position -= 1
    return position, ref, alt


def _fetch_base(reference: ReferenceLookup, contig: str, position: int) -> str | None:
    try:
        sequence = reference.fetch(contig, position, position)
    except Exception:  # a caller-supplied lookup may raise anything
        return None
    base = sequence.strip().upper()
    return base if len(base) == 1 and base in _NUCLEOTIDES else None


def _reference_matches(record: VariantRecord, reference: ReferenceLookup) -> bool | None:
    """``True``/``False`` for match/mismatch, ``None`` when the lookup was unusable.

    ``None`` is distinct from ``False`` on purpose: a lookup that failed is absence
    of information, not evidence of a mismatch (GP-14).
    """
    coordinate = record.coordinate
    if not _is_simple_allele(coordinate.ref):
        return None
    try:
        sequence = reference.fetch(coordinate.contig, coordinate.position, coordinate.end)
    except Exception:  # a caller-supplied lookup may raise anything
        return None
    observed = sequence.strip().upper()
    if len(observed) != len(coordinate.ref):
        return None
    return observed == coordinate.ref


def _is_simple_allele(allele: str) -> bool:
    return bool(allele) and set(allele) <= _NUCLEOTIDES


def _with_ops(operations: tuple[str, ...], operation: str) -> tuple[str, ...]:
    """Append an operation, order-stable and never duplicated."""
    if operation in operations:
        return operations
    return (*operations, operation)
