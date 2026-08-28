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

And one property of *failure* decides the rest:

* **A reference that cannot be read is not the same as a position with no base.**
  Position 0 is before the start of a contig and the base past the end of a contig
  does not exist; stopping there is the reference being complete. A FASTA that
  raises, that is missing the contig, that returns an empty string mid-contig or
  that returns something which is not a nucleotide is *broken*, and nothing about
  the record is known. Both used to arrive at the same silent ``None``, which
  turned an I/O error into a quality downgrade no caller could observe — while the
  adapters went on reporting that left-alignment had been applied. The distinction
  is now made in :func:`_read_one_base` and carried out of here in
  :class:`ReferenceStatus`, on the same principle as
  :attr:`CanonicalAllele.operations`: never claim what did not happen. ADR 0026
  records the reasoning and the limits that remain.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Final, Protocol

from mva.errors import ReferenceUnusableError
from mva.models.base import error_token
from mva.models.variant import OP_LEFT_ALIGN, OP_TRIM

#: Bases a plain, non-symbolic allele may contain. ``<DEL>``, ``*`` and ``.`` are
#: legal VCF and are deliberately excluded: none of them can be trimmed or shifted
#: as a sequence, so they are returned untouched rather than mangled.
NUCLEOTIDES: Final[frozenset[str]] = frozenset("ACGTN")

#: Ceiling on shift iterations, in either direction. A pathological or wrapping
#: :class:`ReferenceLookup` must not be able to spin the loop forever, and 1 kb is
#: far beyond any indel this pipeline is entitled to reason about. Reaching the
#: ceiling stops the shift where it is; the representation stays valid, it is just
#: not proven left-most — which is reported as
#: :attr:`ReferenceStatus.SHIFT_LIMIT_REACHED` rather than passed off as success.
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


class ReferenceStatus(StrEnum):
    """What the reference was able to tell the rule about *this one* allele.

    Per-allele, and deliberately not the same axis as
    :class:`LeftAlignmentStatus`, which is a per-batch summary. A caller building
    a batch report derives the batch status from these; a caller building a single
    join key branches on this one directly.
    """

    NOT_SUPPLIED = "not_supplied"
    """No :class:`ReferenceLookup` was passed. The result is trimmed only, and an
    indel may be spelled differently from a left-aligned source."""

    USABLE = "usable"
    """A reference was supplied and every base the rule needed from it was read.

    Vacuously true when the rule needed no base at all — a SNV, an allele already
    proven left-most by its own last bases, or a variant at position 1. The claim
    being made is "nothing the rule asked for was denied", not "a read happened",
    and for representation purposes those have the same consequence: the result
    *is* the left-most spelling."""

    UNUSABLE = "unusable"
    """A base the rule needed could not be read: the accessor raised, the contig is
    absent, the sequence came back empty mid-contig, or it was not a nucleotide.

    **The result is trimmed only and is not proven left-most.** It must never be
    reported as reference-backed canonicalisation. This is the state that used to
    be indistinguishable from success, and the one that silently turned a common
    gnomAD allele into "frequency unknown"."""

    SHIFT_LIMIT_REACHED = "shift_limit_reached"
    """The shift ran for :data:`MAX_SHIFT_BP` iterations and was still moving.

    **Deliberately not** :attr:`UNUSABLE`. Nothing is wrong with the reference —
    every base it was asked for was supplied — and nothing is wrong with the
    record. What ran out is this module's own iteration budget, which is a
    property of the *code*, and an operator who reads "the reference is unusable"
    will go and check a FASTA that is perfectly healthy.

    The consequence for the join key is nevertheless the same as ``UNUSABLE``'s:
    the allele sits somewhere inside a repeat tract at neither its input position
    nor its proven left-most one, so :attr:`CanonicalAllele.left_alignment_proven`
    is False and the key must not be reported as reference-backed. It used to
    return ``USABLE`` here — a partial join key marked proven, uncounted by
    ``unaligned_indel_count``, in a run whose report went on to state that every
    indel had been placed against the reference.

    Reachable only in a tract longer than 1 kb of the same repeat unit. Rare, and
    exactly the kind of locus where a mis-joined indel is most likely to be
    scored novel and ultra-rare."""


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

    ``reference_status`` is that same discipline applied to the failure case.
    ``operations`` answers "did the position move"; a position that did *not* move
    is ambiguous on its own — it may be already left-most, or the reference may
    have been unreadable. Only ``reference_status`` separates those, and the
    separation is the difference between a trustworthy join key and a silent miss.
    """

    position: int
    ref: str
    alt: str
    operations: tuple[str, ...]
    reference_status: ReferenceStatus = ReferenceStatus.NOT_SUPPLIED

    @property
    def trimmed(self) -> bool:
        return OP_TRIM in self.operations

    @property
    def left_aligned(self) -> bool:
        """True only when the position actually moved leftwards."""
        return OP_LEFT_ALIGN in self.operations

    @property
    def left_alignment_proven(self) -> bool:
        """True when this spelling is the left-most one, and that was *established*.

        The property an adapter must branch on before claiming
        :attr:`LeftAlignmentStatus.APPLIED`. False both when no reference was
        supplied and when the one supplied could not be read — two different
        operator problems (configure a FASTA / fix the FASTA you configured) with
        the same consequence for this allele's join key.
        """
        return self.reference_status is ReferenceStatus.USABLE

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
        The canonical allele, carrying both what was done to it (``operations``)
        and how far the reference could be trusted (``reference_status``).
        Non-sequence alleles (``*``, ``.``, ``<DEL>``, an empty string) are returned
        verbatim with no operations, because there is no defined trimming or
        shifting of a symbolic allele and inventing one would move a coordinate on
        a guess.

    Does not raise on a broken reference. This function's return type has room for
    the degraded state, so it reports it rather than raising; a per-record
    exception in a whole-callset loop would turn one bad base into a dead run,
    which GP-13 forbids for a record-level problem. The one thing it may not do is
    stay silent — see :class:`ReferenceStatus`.
    """
    if not is_sequence_allele(ref) or not is_sequence_allele(alt):
        return CanonicalAllele(
            position=position,
            ref=ref,
            alt=alt,
            operations=(),
            reference_status=(
                ReferenceStatus.NOT_SUPPLIED if reference is None else ReferenceStatus.USABLE
            ),
        )

    trimmed_position, trimmed_ref, trimmed_alt = trim_parsimoniously(position, ref, alt)
    trimmed = (trimmed_position, trimmed_ref, trimmed_alt) != (position, ref, alt)

    shifted = False
    status = ReferenceStatus.NOT_SUPPLIED
    if reference is not None:
        status = ReferenceStatus.USABLE
        try:
            shift_position, shift_ref, shift_alt, limit_reached = _left_shift(
                contig, trimmed_position, trimmed_ref, trimmed_alt, reference
            )
        except ReferenceUnusableError:
            # The single place in this module that converts the loud failure into a
            # declared state, and the only legitimate thing to do with it is report
            # it. The *trimmed* form is kept rather than the partially shifted one:
            # an allele abandoned half way through a repeat tract sits at neither
            # the input position nor the left-most one, so it would join against
            # neither the source nor another run of this pipeline over the same
            # input. Trim-only is at least a defined, reproducible representation.
            status = ReferenceStatus.UNUSABLE
        else:
            if limit_reached:
                # The budget ran out mid-tract. The reference answered every read,
                # so this is not UNUSABLE — but the allele is not left-most either,
                # and saying USABLE here is the claim that used to be made and was
                # never true. The partially shifted form is kept: unlike the
                # unreadable-reference case there is no doubt about the bases, and
                # every run over this input stops at the same place, so the key
                # stays reproducible even though it is not proven left-most.
                status = ReferenceStatus.SHIFT_LIMIT_REACHED
            shifted = (shift_position, shift_ref, shift_alt) != (
                trimmed_position,
                trimmed_ref,
                trimmed_alt,
            )
            # Re-trim: the shift loop grows the alleles leftwards by one anchor base
            # at a time, and the standard requires the final form to be parsimonious.
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
        reference_status=status,
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
) -> tuple[int, str, str, bool]:
    """Roll an indel as far left as the reference allows (Tan et al. Algorithm 1).

    While the alleles end in the same base and either is length 1, prepend the
    preceding reference base to both and step POS back by one. A complex
    substitution — both alleles longer than one base after trimming — is not a
    shiftable indel and is left where it is.

    Every base this loop asks for is one it has already established the shift needs
    and that lies at ``position - 1`` with ``position > 1``, i.e. inside a contig
    the record claims to be on. A reference that cannot supply it is broken, not
    exhausted, so this raises rather than stopping — stopping would be
    indistinguishable from "already left-most".

    Returns the shifted allele and **whether the iteration budget ran out while it
    was still moving**. The flag is not optional bookkeeping. Every ``break`` above
    is a *reason* the shift finished — left-most, complex substitution, start of
    contig — and falling out of the ``for`` is the one case that is not a reason
    but a stop, leaving the allele at neither its input position nor its left-most
    one. Returning only the tuple made those indistinguishable, and the caller
    then labelled the stop ``USABLE``.

    Raises:
        ReferenceUnusableError: the reference could not supply a base the shift
            required.
    """
    for _ in range(MAX_SHIFT_BP):
        while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        if len(ref) > 1 and len(alt) > 1:
            return position, ref, alt, False  # a complex substitution, not shiftable
        if ref[-1] != alt[-1]:
            return position, ref, alt, False  # already left-most
        if position <= 1:
            return position, ref, alt, False  # nothing to the left of a contig's first base
        base = _required_base(reference, contig, position - 1)
        ref, alt = base + ref, base + alt
        position -= 1
    return position, ref, alt, True


@dataclass(frozen=True, slots=True)
class QueryBound:
    """How far right a region query must reach, and whether that bound was proven.

    The two halves have to travel together. A bound alone cannot distinguish "the
    reference proved there is no equivalent spelling further right" from "the
    reference could not be read, so this is merely the furthest point I got to" —
    and those have opposite consequences for whether an empty query result means
    the source has no record.
    """

    position: int
    reference_status: ReferenceStatus

    @property
    def proven(self) -> bool:
        """True when the bound is the reference's answer rather than a stopping point."""
        return self.reference_status is ReferenceStatus.USABLE


def rightmost_equivalent_bound(
    *,
    contig: str,
    position: int,
    ref: str,
    alt: str,
    reference: ReferenceLookup,
) -> QueryBound:
    """POS of the right-most legal spelling of this same event, plus its provenance.

    **This is the form callers should use.** :func:`rightmost_equivalent_position`
    returns the same number with the provenance discarded, which is the shape that
    made an unreadable reference indistinguishable from a proven bound.

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

    Two ways of running out of reference, and they are not the same:

    * **Past the end of the contig.** This loop reads *beyond* the record's own
      span, so an empty read there is the reference correctly reporting where the
      chromosome stops. The bound is complete and ``reference_status`` stays
      ``USABLE``.
    * **The reference is broken.** It raised, the contig is absent, or it returned
      something that is not a nucleotide. The loop stops where it is and the bound
      is reported ``UNUSABLE`` — because a bound that stops early makes the fetch
      window too small, and a source record lost to the *fetch* is
      indistinguishable, from outside, from one lost to the key, which is in turn
      indistinguishable from the source genuinely having nothing.

    An under-reaching bound is not, on its own, a new failure: an unreadable
    reference has already left the *query key* trimmed-only, so a right-shifted
    source record would not have joined even if it had been fetched. What matters
    is that the caller can tell, which is what this return type is for.
    """
    if not is_sequence_allele(ref) or not is_sequence_allele(alt):
        return QueryBound(position=position, reference_status=ReferenceStatus.USABLE)
    position, ref, alt = trim_parsimoniously(position, ref, alt)
    if len(ref) == len(alt):
        # A substitution occupies fixed bases and cannot shift.
        return QueryBound(position=position, reference_status=ReferenceStatus.USABLE)
    for _ in range(MAX_SHIFT_BP):
        if len(ref) > 1 and len(alt) > 1:
            break  # a complex substitution, not a shiftable indel
        try:
            base = _read_one_base(reference, contig, position + len(ref))
        except ReferenceUnusableError:
            return QueryBound(position=position, reference_status=ReferenceStatus.UNUSABLE)
        if base is None:
            break  # past the end of the contig: absence, not breakage
        grown_ref, grown_alt = ref + base, alt + base
        moved, next_ref, next_alt = _trim_shared_prefix(position, grown_ref, grown_alt)
        if moved == position:
            break  # growing rightwards bought no shift; this is the right-most form
        position, ref, alt = moved, next_ref, next_alt
    else:
        # Fell out of the loop still moving. Every ``break`` above is a reason the
        # search finished; this is the one exit that is not. The bound is therefore
        # a stopping point, not the reference's answer, and reporting it USABLE
        # told the caller the query window provably covers every equivalent
        # spelling when it demonstrably does not. Not UNUSABLE either: the
        # reference supplied every base it was asked for.
        return QueryBound(position=position, reference_status=ReferenceStatus.SHIFT_LIMIT_REACHED)
    return QueryBound(position=position, reference_status=ReferenceStatus.USABLE)


def rightmost_equivalent_position(
    *,
    contig: str,
    position: int,
    ref: str,
    alt: str,
    reference: ReferenceLookup,
) -> int:
    """The bound from :func:`rightmost_equivalent_bound`, with its provenance dropped.

    **Retained for the two annotation adapters that already call it, and not the
    form to write new code against.** An ``int`` has nowhere to say "I could not
    compute this", so a caller of this function cannot tell a proven right-most
    position from the point at which a broken reference stopped the search — and
    that is precisely the class of silence ADR 0018 exists to remove.

    Migrating a caller is one line: take :class:`QueryBound` and branch on
    :attr:`QueryBound.proven` before treating an empty query result as evidence
    that the source holds no record (GP-14). The exact diff for both adapters is in
    ``docs/handoff-integrity.md``.
    """
    return rightmost_equivalent_bound(
        contig=contig, position=position, ref=ref, alt=alt, reference=reference
    ).position


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


def _read_one_base(reference: ReferenceLookup, contig: str, position: int) -> str | None:
    """One reference base; ``None`` when the reference *has* no base there.

    The whole point of this function is that those two sentences are different
    things, and that only one of them is a property of the data:

    * ``None`` — the position lies outside the reference. Before the first base of
      a contig, or past its last. The reference is working and is telling the
      truth; a caller that reads beyond a record's own span must accept this.
    * :class:`~mva.errors.ReferenceUnusableError` — the accessor raised, the contig
      is not in the index, the sequence came back empty where a base must exist, or
      it came back as something that is not a nucleotide. Nothing is known.

    A caller-supplied lookup may raise anything, so everything is caught — but it
    is caught in order to be *re-raised as a named condition*, not to be turned
    into a return value that reads like ordinary absence. The original exception is
    suppressed with ``from None`` rather than chained: a genomics backend routinely
    puts the region it was asked for into its own message, and a chained traceback
    would carry a proband coordinate into every terminal, log and crash report the
    traceback reaches (PRIV-09). The exception *type* is safe and is kept, because
    ``OSError`` and ``KeyError`` mean very different things to whoever debugs this.
    """
    if position < 1:
        return None
    try:
        sequence = reference.fetch(contig, position, position)
    except Exception as exc:
        raise ReferenceUnusableError(
            _unusable_message(contig, position, reason=f"the accessor raised {type(exc).__name__}")
        ) from None
    base = sequence.strip().upper()
    if not base:
        return None
    if len(base) != 1 or base not in NUCLEOTIDES:
        raise ReferenceUnusableError(
            _unusable_message(
                contig,
                position,
                reason=(
                    f"the accessor returned {len(base)} character(s) that are not a single "
                    f"nucleotide from {''.join(sorted(NUCLEOTIDES))}"
                ),
            )
        )
    return base


def _required_base(reference: ReferenceLookup, contig: str, position: int) -> str:
    """A base the rule has already established must exist, or the reference is broken.

    Used by the left shift, which only ever reads at ``position - 1`` for a
    ``position`` greater than 1 on a contig the record claims to sit on. There is
    no legitimate way for that base to be missing, so an empty read here is
    breakage — a truncated FASTA, a mis-built ``.fai``, or the wrong assembly —
    and not the end of the contig.
    """
    base = _read_one_base(reference, contig, position)
    if base is None:
        raise ReferenceUnusableError(
            _unusable_message(
                contig,
                position,
                reason=(
                    "the accessor returned no sequence for a position inside the contig, "
                    "which usually means a truncated FASTA, a stale .fai index, or a "
                    "reference shorter than the assembly the records were called against"
                ),
            )
        )
    return base


def _unusable_message(contig: str, position: int, *, reason: str) -> str:
    """PRIV-09-safe diagnostic: the problem and the field, never the coordinate."""
    return (
        "The reference could not supply a base that left-alignment required: "
        f"{reason}. Locus handle <locus:{error_token((contig, position))}>; the contig "
        "and position are tokenised rather than echoed (PRIV-09), and the handle is "
        "stable within this run so repeated failures at one locus correlate. This is a "
        "broken reference, NOT a position that legitimately has no base — records "
        "affected by it are trimmed only and must never be reported as left-aligned "
        "(ADR 0018, GP-14)."
    )


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

    INCOMPLETE_SHIFT_LIMIT = "incomplete_shift_limit"
    """The reference answered every read, but some indels sat in a repeat tract
    longer than :data:`MAX_SHIFT_BP` and the shift stopped before reaching the
    left-most position.

    A separate value from ``INCOMPLETE_REFERENCE_UNUSABLE`` because the operator
    action is different — nothing is wrong with the FASTA and re-downloading it
    changes nothing — while the consequence for those records' join keys is the
    same: they are not left-most and may not match a left-aligned source. Kept
    apart from ``APPLIED`` for the reason this enum exists at all: "every indel
    was placed against the reference" was being stated over a batch in which some
    had not been."""


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

    shift_limited_count: int = 0
    """Indels whose shift hit :data:`MAX_SHIFT_BP` while still moving.

    Disjoint from ``unaligned_indel_count``, which counts records the reference
    could not be *read* for. These were read successfully and are still not
    left-most, and before this field existed they were counted as neither — so a
    run could report ``APPLIED`` and "every indel was placed against it" over a
    batch holding partially shifted keys."""

    @property
    def is_degraded(self) -> bool:
        """True when at least one indel's key is not proven left-most.

        Both ways of failing to reach the left-most position count, because the
        consequence downstream is identical: a key that may not match a
        left-aligned source, whose miss is then read as "no record" and scored
        novel and ultra-rare. The two causes are told apart by
        :attr:`status`, not by whether the run is degraded at all.
        """
        return self.unaligned_indel_count > 0 or self.shift_limited_count > 0

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
        if self.status is LeftAlignmentStatus.INCOMPLETE_SHIFT_LIMIT:
            return (
                f"DEGRADED: the reference was readable throughout, but "
                f"{self.shift_limited_count} of {self.indel_count} indel records sit in a "
                f"repeat tract longer than the {MAX_SHIFT_BP} bp shift budget and stopped "
                "short of their left-most position. Their join keys are reproducible but "
                "not left-most, so they may fail to match ClinVar or gnomAD and would then "
                "be scored as novel. The FASTA is not at fault and re-downloading it will "
                "not change this; raising MAX_SHIFT_BP is the only remedy."
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
            "shift_limited_count": self.shift_limited_count,
            "degraded": self.is_degraded,
            "summary": self.describe(),
        }


def summarise_left_alignment(
    *,
    indel_count: int,
    shifted_count: int,
    unaligned_indel_count: int,
    reference_available: bool,
    shift_limited_count: int = 0,
) -> LeftAlignmentReport:
    """Derive the status from the counts, so no caller can label a batch by hand.

    Precedence when more than one degradation is present is worst-first:
    no reference at all, then records the reference could not be read for, then
    records that ran out of shift budget. The report carries every count either
    way, so nothing is lost by the status naming only the most severe.
    """
    if indel_count == 0:
        status = LeftAlignmentStatus.NOT_REQUIRED
    elif not reference_available:
        status = LeftAlignmentStatus.UNAVAILABLE_NO_REFERENCE
    elif unaligned_indel_count:
        status = LeftAlignmentStatus.INCOMPLETE_REFERENCE_UNUSABLE
    elif shift_limited_count:
        status = LeftAlignmentStatus.INCOMPLETE_SHIFT_LIMIT
    else:
        status = LeftAlignmentStatus.APPLIED
    return LeftAlignmentReport(
        status=status,
        indel_count=indel_count,
        shifted_count=shifted_count,
        unaligned_indel_count=unaligned_indel_count,
        reference_available=reference_available,
        shift_limited_count=shift_limited_count,
    )
