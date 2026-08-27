"""The one representation rule, shared by every stage that builds a join key.

Two spellings of the same event — ``chr1:100 AT>AG`` and ``chr1:101 T>G``,
``chr1:101 A>AA`` and ``chr1:100 A>AA`` — are the same biology and different
strings. Every join in this pipeline is a string comparison on
``GenomicCoordinate.variant_id``, so a disagreement about representation does not
surface as an error. It surfaces as *absence*: the frequency is not found, the
ClinVar assertion is not found, and the variant is scored as novel and ultra-rare,
which is the strongest promoting signal the ranker has. A representation bug
therefore manufactures false positives and deletes true pathogenic assertions at
the same time, silently, in both directions.

This module exists because that had already happened twice, in two places, with
two implementations. The ingestion normaliser trimmed; the ClinVar adapter did
not. Neither was wrong on its own terms and together they could not join. **The
rule lives here once and both sides call it**; a second implementation is the
defect, not the cure.

The algorithm is not invented here. It is the one specified in:

* Tan, A., Abecasis, G. R. & Kang, H. M. (2015) "Unified representation of genetic
  variants", *Bioinformatics* **31**(13):2202-4 — Algorithm 1, "variant
  normalization": left-align while the alleles end in the same base and one of
  them is length 1, then trim the shared suffix and the shared prefix, always
  keeping at least one base.
* The VCF 4.2 specification, §1.4.1 and §5, which fixes the anchor-base convention
  an indel is written in and requires POS to be the left-most position at which
  the event can be placed.

Two properties of that standard decide the shape of everything below:

* **Trimming needs no reference.** It is pure string surgery on REF and ALT and is
  always available, everywhere, to every caller.
* **Left-alignment needs the reference.** Rolling an indel leftwards through a
  homopolymer or a tandem repeat means reading the bases to its left. A caller
  without a :class:`ReferenceLookup` cannot do it, must not pretend it did, and
  must say so — see :class:`LeftAlignmentReport`.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Final, Protocol

from mva.models.variant import OP_LEFT_ALIGN, OP_TRIM

#: Bases a plain, non-symbolic allele may contain. ``<DEL>``, ``*`` and ``.`` are
#: legal VCF and are deliberately excluded: none of them can be trimmed or shifted
#: as a sequence, so they are returned untouched rather than mangled.
NUCLEOTIDES: Final[frozenset[str]] = frozenset("ACGTN")

#: Ceiling on shift iterations, in either direction. A pathological or wrapping
#: :class:`ReferenceLookup` must not be able to spin the loop forever, and 1 kb is
#: far beyond any indel this pipeline is entitled to reason about. Reaching the
#: ceiling stops the shift where it is; the representation stays valid, it is just
#: not proven left-most.
MAX_SHIFT_BP: Final = 1000


class ReferenceLookup(Protocol):
    """Minimal reference-sequence accessor.

    ``start`` and ``end`` are **1-based and inclusive**, matching VCF convention
    rather than pysam's 0-based half-open ``fetch``. The adapter that wraps a FASTA
    is responsible for that translation, once, at the boundary (GP-02) — see
    :class:`mva.ingestion.normalise.FastaReference`.
    """

    def fetch(self, contig: str, start: int, end: int) -> str:
        """Return the reference bases in ``[start, end]`` on ``contig``."""
        ...


@dataclass(frozen=True, slots=True)
class CanonicalAllele:
    """A canonicalised ``(position, ref, alt)`` and the account of how it got there.

    ``operations`` is the audit trail required by GP-31: left-alignment *changes a
    variant's POS*, and the submission this pipeline feeds is scored on exact
    coordinates. A coordinate that moved with no record of the move is
    unauditable, so every mutation is named here and the caller copies the names
    into ``VariantRecord.normalisation_ops``.

    It names only what actually happened. ``left_align`` is absent when no
    reference was supplied, because recording an operation that was never
    performed is a provenance lie no downstream consumer can detect.
    """

    position: int
    ref: str
    alt: str
    operations: tuple[str, ...]

    @property
    def trimmed(self) -> bool:
        return OP_TRIM in self.operations

    @property
    def left_aligned(self) -> bool:
        """True only when the position actually moved leftwards."""
        return OP_LEFT_ALIGN in self.operations

    @property
    def changed(self) -> bool:
        return bool(self.operations)

    @property
    def is_indel(self) -> bool:
        return len(self.ref) != len(self.alt)


# ---------------------------------------------------------------------------
# The rule
# ---------------------------------------------------------------------------


def canonicalise_allele(
    *,
    contig: str,
    position: int,
    ref: str,
    alt: str,
    reference: ReferenceLookup | None = None,
) -> CanonicalAllele:
    """Reduce one allele pair to its minimal, left-most representation.

    This is the function every join key in the pipeline must be built through.
    Ingestion enters here via :func:`mva.ingestion.normalise.trim_and_left_align`;
    the ClinVar adapter enters here for both the query it is given and the record
    it reads out of the release.

    Args:
        contig: needed only to read the reference; unused when ``reference`` is None.
        position: 1-based VCF POS.
        ref: REF allele, already upper-cased.
        alt: ALT allele, already upper-cased.
        reference: 1-based inclusive accessor. Without it the result is trimmed but
            **not** left-aligned, and ``operations`` will not claim otherwise.

    Returns:
        The canonical allele. Non-sequence alleles (``*``, ``.``, ``<DEL>``, an
        empty string) are returned verbatim with no operations, because there is no
        defined trimming or shifting of a symbolic allele and inventing one would
        move a coordinate on a guess.
    """
    if not is_sequence_allele(ref) or not is_sequence_allele(alt):
        return CanonicalAllele(position=position, ref=ref, alt=alt, operations=())

    trimmed_position, trimmed_ref, trimmed_alt = trim_parsimoniously(position, ref, alt)
    trimmed = (trimmed_position, trimmed_ref, trimmed_alt) != (position, ref, alt)

    shifted = False
    if reference is not None:
        shift_position, shift_ref, shift_alt = _left_shift(
            contig, trimmed_position, trimmed_ref, trimmed_alt, reference
        )
        shifted = (shift_position, shift_ref, shift_alt) != (
            trimmed_position,
            trimmed_ref,
            trimmed_alt,
        )
        # Re-trim: the shift loop grows the alleles leftwards by one anchor base at
        # a time, and the standard requires the final form to be parsimonious.
        after = trim_parsimoniously(shift_position, shift_ref, shift_alt)
        trimmed = trimmed or after != (shift_position, shift_ref, shift_alt)
        trimmed_position, trimmed_ref, trimmed_alt = after

    operations: list[str] = []
    if shifted:
        operations.append(OP_LEFT_ALIGN)
    if trimmed:
        operations.append(OP_TRIM)
    return CanonicalAllele(
        position=trimmed_position,
        ref=trimmed_ref,
        alt=trimmed_alt,
        operations=tuple(operations),
    )


def trim_parsimoniously(position: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Strip the shared suffix, then the shared prefix, always keeping >= 1 base.

    Tan et al. (2015) Algorithm 1, steps 2-3, and the VCF 4.2 anchor-base
    convention: an indel retains one anchoring base, so ``G>GAT`` is already
    minimal while ``AAT>AAG`` reduces to ``T>G`` two bases to the right.

    Suffix first, then prefix: doing it the other way round turns ``AT>AG`` into a
    prefix-trimmed ``T>G`` correctly, but turns ``TA>TAA`` into ``A>AA`` with the
    anchor on the wrong side. Order is part of the standard, not a preference.
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
    """Roll an indel as far left as the reference allows (Tan et al. Algorithm 1).

    While the alleles end in the same base and either is length 1, prepend the
    preceding reference base to both and step POS back by one. A complex
    substitution — both alleles longer than one base after trimming — is not a
    shiftable indel and is left where it is.
    """
    for _ in range(MAX_SHIFT_BP):
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


def rightmost_equivalent_position(
    *,
    contig: str,
    position: int,
    ref: str,
    alt: str,
    reference: ReferenceLookup,
) -> int:
    """POS of the right-most legal spelling of this same event.

    The mirror of :func:`_left_shift`, and it exists for one job: deciding how far
    to the *right* a region query must reach to be sure it has seen every source
    record that canonicalises back to this variant.

    An index query finds records by the span they occupy. A left-aligned query at
    POS 100 and a source record that spells the identical insertion at POS 106
    occupy disjoint spans, so the record is never fetched, never canonicalised, and
    the join fails exactly as if the source had no record — the failure mode this
    whole module exists to remove. Querying out to the right-most equivalent
    position closes that hole with a bound computed from the reference rather than
    a guessed padding constant. For a variant not in a repeat this returns
    ``position`` unchanged and costs nothing.

    The bound is deliberately not tight: the loop stops only when growing rightwards
    buys no further shift, which can leave it one base past the right-most
    *VCF-conventional* spelling (VCF puts the padding base before the event; this
    loop will also accept it after). Over-reaching costs one base of index query.
    Under-reaching re-opens the silent miss, so the asymmetry is the right way round.

    Returns ``position`` for symbolic alleles and for complex substitutions, which
    do not shift.
    """
    if not is_sequence_allele(ref) or not is_sequence_allele(alt):
        return position
    position, ref, alt = trim_parsimoniously(position, ref, alt)
    for _ in range(MAX_SHIFT_BP):
        if len(ref) > 1 and len(alt) > 1:
            break  # a complex substitution, not a shiftable indel
        base = _fetch_base(reference, contig, position + len(ref))
        if base is None:
            break
        grown_ref, grown_alt = ref + base, alt + base
        moved, next_ref, next_alt = _trim_shared_prefix(position, grown_ref, grown_alt)
        if moved == position:
            break  # growing rightwards bought no shift; this is the right-most form
        position, ref, alt = moved, next_ref, next_alt
    return position


def _trim_shared_prefix(position: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Prefix-only trim. Used by the right shift, where a suffix trim would undo it."""
    limit = min(len(ref), len(alt)) - 1
    shared = 0
    while shared < limit and ref[shared] == alt[shared]:
        shared += 1
    if shared:
        ref, alt = ref[shared:], alt[shared:]
        position += shared
    return position, ref, alt


def _fetch_base(reference: ReferenceLookup, contig: str, position: int) -> str | None:
    """One reference base, or ``None`` when the lookup cannot supply it.

    ``None`` stops the shift rather than raising: an unreadable base is absence of
    information (GP-14), and the un-shifted representation is still a valid one.
    A caller-supplied lookup may raise anything, so nothing narrower is caught.
    """
    if position < 1:
        return None
    try:
        sequence = reference.fetch(contig, position, position)
    except Exception:
        return None
    base = sequence.strip().upper()
    return base if len(base) == 1 and base in NUCLEOTIDES else None


def is_sequence_allele(allele: str) -> bool:
    """True for a plain ACGTN allele: the only kind that may be trimmed or shifted."""
    return bool(allele) and set(allele) <= NUCLEOTIDES


# ---------------------------------------------------------------------------
# The degraded state (GP-14)
# ---------------------------------------------------------------------------


class LeftAlignmentStatus(StrEnum):
    """Whether a batch's indels were actually placed at their left-most position.

    Typed rather than a boolean or a log line, because "we could not left-align"
    and "there was nothing to left-align" have opposite meanings for how much a
    reader should trust the rarity of every indel in the run, and a warning string
    cannot be branched on.
    """

    APPLIED = "applied"
    """A reference was supplied and every indel was placed against it."""

    NOT_REQUIRED = "not_required"
    """The batch contains no indels. SNVs are unaffected by left-alignment."""

    UNAVAILABLE_NO_REFERENCE = "unavailable_no_reference"
    """No reference was supplied. **Indel joins are unreliable in both directions.**"""

    INCOMPLETE_REFERENCE_UNUSABLE = "incomplete_reference_unusable"
    """A reference was supplied but could not be read, or disagreed with REF, for
    some records. Those records were trimmed only."""


@dataclass(frozen=True, slots=True)
class LeftAlignmentReport:
    """What a run may and may not claim about its indel representations.

    GP-14: without a reference FASTA, left-alignment is impossible, and the honest
    output is not "left-alignment skipped" buried in a log. It is a first-class
    statement that every indel frequency and every indel clinical assertion in this
    run may be missing for representational reasons, carried in the result object
    so a caller has to receive it.

    Counts only. No coordinate, allele or genotype appears in any field or in
    :meth:`describe` (PRIV-09).
    """

    status: LeftAlignmentStatus
    indel_count: int
    shifted_count: int
    unaligned_indel_count: int
    reference_available: bool

    @property
    def is_degraded(self) -> bool:
        """True when at least one indel was never checked against a reference."""
        return self.unaligned_indel_count > 0

    def describe(self) -> str:
        """One sentence a reader who has never seen this code can act on."""
        if self.status is LeftAlignmentStatus.NOT_REQUIRED:
            return (
                "Left-alignment was not required: this run contains no indel records, "
                "and SNV representation is unaffected by it."
            )
        if self.status is LeftAlignmentStatus.APPLIED:
            return (
                f"Left-alignment applied against the configured reference to all "
                f"{self.indel_count} indel records; {self.shifted_count} moved to a "
                "left-most position. Indel join keys agree with left-aligned sources."
            )
        if self.status is LeftAlignmentStatus.UNAVAILABLE_NO_REFERENCE:
            return (
                f"DEGRADED: no reference FASTA was supplied, so none of the "
                f"{self.indel_count} indel records in this run could be left-aligned. "
                "ClinVar and gnomAD store left-aligned alleles, so any indel this VCF "
                "spells at a different position simply will not join. A failed join is "
                "indistinguishable from 'no record found', which this pipeline scores as "
                "novel and ultra-rare — the strongest promoting signal it has. Set "
                "inputs.reference_fasta in the case config to remove this limitation."
            )
        return (
            f"DEGRADED: a reference was supplied but was unusable for "
            f"{self.unaligned_indel_count} of {self.indel_count} indel records, which "
            "were trimmed only. Those records may fail to join against left-aligned "
            "sources and would then be scored as novel. Check that the FASTA is the "
            "same assembly and patch release as the VCF."
        )

    def as_dict(self) -> dict[str, object]:
        """JSON-safe, fixed key order, for a QC artifact a reader actually opens."""
        return {
            "status": self.status.value,
            "reference_available": self.reference_available,
            "indel_count": self.indel_count,
            "shifted_count": self.shifted_count,
            "unaligned_indel_count": self.unaligned_indel_count,
            "degraded": self.is_degraded,
            "summary": self.describe(),
        }


def summarise_left_alignment(
    *,
    indel_count: int,
    shifted_count: int,
    unaligned_indel_count: int,
    reference_available: bool,
) -> LeftAlignmentReport:
    """Derive the status from the counts, so no caller can label a batch by hand."""
    if indel_count == 0:
        status = LeftAlignmentStatus.NOT_REQUIRED
    elif not reference_available:
        status = LeftAlignmentStatus.UNAVAILABLE_NO_REFERENCE
    elif unaligned_indel_count:
        status = LeftAlignmentStatus.INCOMPLETE_REFERENCE_UNUSABLE
    else:
        status = LeftAlignmentStatus.APPLIED
    return LeftAlignmentReport(
        status=status,
        indel_count=indel_count,
        shifted_count=shifted_count,
        unaligned_indel_count=unaligned_indel_count,
        reference_available=reference_available,
    )
