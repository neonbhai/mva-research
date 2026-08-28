"""A broken reference must not be reported as a successful left-alignment.

`mva.alleles` used to catch *every* exception from `ReferenceLookup.fetch` and
return ``None``, which the shift loop read as "no base here, stop". Two very
different things arrived at that same ``None``:

* **There is genuinely no reference base at this position.** Position 0 is before
  the start of a contig; the base one past the end of a contig does not exist.
  Stopping is correct, and the un-shifted representation is a valid one.
* **The reference is broken.** The FASTA raised ``OSError``, the contig is not in
  the index, the accessor returned an empty string in the middle of a contig, or
  it returned something that is not a nucleotide. Nothing about the variant is
  known, and the result is *not* proven left-most.

Collapsing the second into the first is a provenance lie, and it is a worse
failure than the one it hides. Both annotation adapters report left-alignment as
``APPLIED`` whenever their reference object is merely non-``None``, so a FASTA
that raises on every read produced trim-only join keys while the run claimed
reference-backed canonicalisation. A left-aligned gnomAD record and a
right-shifted proband spelling then fail to join, the variant is scored as having
no frequency data, and GP-14 is violated in the most expensive direction: absence
of information is read as evidence of rarity, which is the strongest promoting
signal the ranker has.

`CanonicalAllele.operations` already exists so that an operation which did not
happen is never claimed. `CanonicalAllele.reference_status` is the same
discipline applied to the failure case: it says whether the reference was
supplied, and whether everything the rule needed to read from it could actually
be read.
"""

from __future__ import annotations

import traceback
from collections.abc import Mapping

import pytest

from mva.alleles import (
    CanonicalAllele,
    LeftAlignmentStatus,
    ReferenceStatus,
    _left_shift,
    canonicalise_allele,
    rightmost_equivalent_bound,
    rightmost_equivalent_position,
    summarise_left_alignment,
)
from mva.errors import ReferenceUnusableError
from mva.models.variant import OP_LEFT_ALIGN

pytestmark = [pytest.mark.unit]

CONTIG = "chr" + "21"

#: A ten-base homopolymer at 100-109 with distinct flanks. An insertion anywhere
#: inside it has ten legal spellings, so it is the case left-alignment exists for.
#: 1-based, inclusive.
REPEAT_TRACT: Mapping[int, str] = {
    **dict.fromkeys(range(90, 100), "C"),
    **dict.fromkeys(range(100, 110), "A"),
    110: "G",
}

#: The right-shifted spelling a caller hands in, inside the tract.
SHIFTED_POSITION = 105
SHIFTED_REF = "A"
SHIFTED_ALT = "AA"


class WorkingReference:
    """The reference the rule is entitled to assume."""

    def __init__(self, bases: Mapping[int, str] = REPEAT_TRACT) -> None:
        self._bases = bases
        self.reads = 0

    def fetch(self, contig: str, start: int, end: int) -> str:
        self.reads += 1
        assert contig == CONTIG
        # Empty beyond the mapped region: that is what a real FASTA returns past
        # the end of a contig, and it is genuine absence, not breakage.
        return "".join(self._bases.get(position, "") for position in range(start, end + 1))


class RaisingReference:
    """A FASTA whose handle has gone away. htslib raises ``OSError`` for this.

    The message deliberately carries a coordinate, because a real backend's does:
    the point of the test is that the coordinate must not survive into the
    exception this module raises (PRIV-09).
    """

    def fetch(self, contig: str, start: int, end: int) -> str:
        msg = f"{contig}:{start}-{end}: htslib read failed"
        raise OSError(msg)


class MissingContigReference:
    """The FASTA is fine but does not contain this contig — a wrong reference."""

    def fetch(self, contig: str, start: int, end: int) -> str:
        raise KeyError(contig)


class EmptyReference:
    """Returns nothing, in the middle of a contig. A truncated or mis-indexed FASTA."""

    def fetch(self, contig: str, start: int, end: int) -> str:
        return ""


class InvalidBaseReference:
    """Returns something that is not a nucleotide — an IUPAC code, or noise."""

    def fetch(self, contig: str, start: int, end: int) -> str:
        return "M" * (end - start + 1)


BROKEN_REFERENCES = [
    pytest.param(RaisingReference(), id="raises-oserror"),
    pytest.param(MissingContigReference(), id="missing-contig"),
    pytest.param(EmptyReference(), id="empty-sequence"),
    pytest.param(InvalidBaseReference(), id="invalid-base"),
]

#: The subset a *rightward* read can still recognise as breakage. An empty read to
#: the right of a record is genuinely ambiguous — see
#: :func:`test_an_empty_read_past_the_end_of_a_record_cannot_be_called_breakage`.
RIGHTWARD_BREAKAGE = [
    pytest.param(RaisingReference(), id="raises-oserror"),
    pytest.param(MissingContigReference(), id="missing-contig"),
    pytest.param(InvalidBaseReference(), id="invalid-base"),
]


def canonicalise_shifted(reference: object | None) -> CanonicalAllele:
    return canonicalise_allele(
        contig=CONTIG,
        position=SHIFTED_POSITION,
        ref=SHIFTED_REF,
        alt=SHIFTED_ALT,
        reference=reference,  # type: ignore[arg-type]
    )


# ---------------------------------------------------------------------------
# The provenance lie
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("reference", BROKEN_REFERENCES)
def test_a_broken_reference_is_never_reported_as_a_successful_left_alignment(
    reference: object,
) -> None:
    """The defect, reproduced four ways.

    Each of these is a distinct way for a FASTA to be unusable, and every one of
    them used to arrive at the same silent ``None`` as "there is no base here".
    """
    canonical = canonicalise_shifted(reference)

    assert canonical.reference_status is ReferenceStatus.UNUSABLE, (
        "a reference that could not be read was not reported as unusable; the caller "
        "cannot tell this apart from a successful left-alignment"
    )
    assert not canonical.left_alignment_proven, (
        "the result claims to be proven left-most although no reference base was read"
    )
    assert OP_LEFT_ALIGN not in canonical.operations, (
        "operations claims an alignment that never happened (GP-31)"
    )


@pytest.mark.parametrize("reference", BROKEN_REFERENCES)
def test_a_broken_reference_yields_exactly_the_no_reference_representation(
    reference: object,
) -> None:
    """The degraded result must be the trim-only one, not a third spelling.

    A shift that failed half way through a repeat tract would leave the allele at
    neither the input position nor the left-most one. That is strictly worse than
    not shifting: it cannot join against the source *or* against another run of
    this pipeline over the same input.
    """
    degraded = canonicalise_shifted(reference)
    trim_only = canonicalise_shifted(None)

    assert (degraded.position, degraded.ref, degraded.alt) == (
        trim_only.position,
        trim_only.ref,
        trim_only.alt,
    )


def test_no_reference_and_an_unusable_reference_are_different_states() -> None:
    """GP-14 cuts both ways: "not supplied" and "supplied but broken" are not one thing.

    They call for different operator actions — configure a FASTA, versus fix the
    FASTA you configured — so a single degraded flag would be a worse report.
    """
    assert canonicalise_shifted(None).reference_status is ReferenceStatus.NOT_SUPPLIED
    assert canonicalise_shifted(RaisingReference()).reference_status is ReferenceStatus.UNUSABLE


def test_a_working_reference_still_left_aligns_and_says_so() -> None:
    """The fix must not turn a healthy reference into a degraded one."""
    reference = WorkingReference()
    canonical = canonicalise_allele(
        contig=CONTIG,
        position=SHIFTED_POSITION,
        ref=SHIFTED_REF,
        alt=SHIFTED_ALT,
        reference=reference,
    )
    assert canonical.reference_status is ReferenceStatus.USABLE
    assert canonical.left_alignment_proven
    assert canonical.left_aligned
    assert canonical.position < SHIFTED_POSITION
    assert reference.reads > 0


def test_a_position_with_no_base_to_its_left_is_absence_not_breakage() -> None:
    """Position 1 has nothing to its left. That is the reference being complete.

    The rule must not report a healthy reference as unusable merely because the
    variant sits at the very start of a contig; over-stating a limitation is its
    own dishonesty (ADR 0018).
    """
    canonical = canonicalise_allele(
        contig=CONTIG,
        position=1,
        ref="A",
        alt="AA",
        reference=WorkingReference({1: "A"}),
    )
    assert canonical.reference_status is ReferenceStatus.USABLE


def test_an_allele_the_rule_never_needs_to_read_for_is_not_degraded() -> None:
    """A SNV is already left-most by construction; no read is required to prove it.

    Reporting it as degraded because some *other* record could not be read would
    make the status useless for deciding whether an indel join can be trusted.
    """
    canonical = canonicalise_allele(
        contig=CONTIG, position=105, ref="A", alt="G", reference=RaisingReference()
    )
    assert canonical.reference_status is ReferenceStatus.USABLE
    assert canonical.left_alignment_proven


# ---------------------------------------------------------------------------
# The query window
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("reference", RIGHTWARD_BREAKAGE)
def test_the_query_bound_says_when_it_is_not_a_proven_bound(reference: object) -> None:
    """A bare ``int`` has nowhere to carry "I could not compute this".

    This bound decides how far right an index query reaches. Stopping early when
    the reference cannot be read makes the fetch window too small, so a
    right-shifted source record is never even retrieved — a join failure
    indistinguishable, from outside, from the source having no record, which is
    the inference GP-14 forbids.
    """
    bound = rightmost_equivalent_bound(
        contig=CONTIG,
        position=SHIFTED_POSITION,
        ref=SHIFTED_REF,
        alt=SHIFTED_ALT,
        reference=reference,  # type: ignore[arg-type]
    )
    assert bound.reference_status is ReferenceStatus.UNUSABLE
    assert not bound.proven


def test_the_query_bound_is_proven_when_it_stops_at_the_end_of_a_contig() -> None:
    """Running off the end of a contig is absence, and does not degrade the bound.

    The right-shift reads *beyond* the record's own span, so an empty read there
    is the reference telling the truth about where the chromosome ends. Reporting
    that as a broken reference would overstate the limitation, which ADR 0018 is
    explicit is its own kind of dishonesty.
    """
    bound = rightmost_equivalent_bound(
        contig=CONTIG,
        position=SHIFTED_POSITION,
        ref=SHIFTED_REF,
        alt=SHIFTED_ALT,
        reference=WorkingReference(),
    )
    assert bound.proven
    assert bound.position >= SHIFTED_POSITION


def test_an_empty_read_past_the_end_of_a_record_cannot_be_called_breakage() -> None:
    """The one direction where the distinction genuinely cannot be made, stated openly.

    The right shift reads *beyond* the record, where "no sequence" is the correct
    answer at the end of a chromosome. A FASTA that returns nothing everywhere is
    therefore indistinguishable, on a rightward read alone, from a variant sitting
    on the last base of a contig. It is reported as a proven bound, and it is the
    *leftward* read — where the base must exist — that catches the same FASTA.

    Documented as a test rather than left implicit, because the alternative is a
    reader later assuming the bound proves more than it does.
    """
    bound = rightmost_equivalent_bound(
        contig=CONTIG,
        position=SHIFTED_POSITION,
        ref=SHIFTED_REF,
        alt=SHIFTED_ALT,
        reference=EmptyReference(),
    )
    assert bound.proven, "a rightward empty read is absence; only the left read can tell"
    assert canonicalise_shifted(EmptyReference()).reference_status is ReferenceStatus.UNUSABLE, (
        "the same reference must still be caught on the leftward read"
    )


def test_the_int_only_bound_still_agrees_with_the_typed_one() -> None:
    """The retained ``-> int`` shim must not be a second implementation (ADR 0018)."""
    for reference in (WorkingReference(), RaisingReference()):
        typed = rightmost_equivalent_bound(
            contig=CONTIG,
            position=SHIFTED_POSITION,
            ref=SHIFTED_REF,
            alt=SHIFTED_ALT,
            reference=reference,  # type: ignore[arg-type]
        )
        plain = rightmost_equivalent_position(
            contig=CONTIG,
            position=SHIFTED_POSITION,
            ref=SHIFTED_REF,
            alt=SHIFTED_ALT,
            reference=reference,  # type: ignore[arg-type]
        )
        assert plain == typed.position


# ---------------------------------------------------------------------------
# PRIV-09
# ---------------------------------------------------------------------------


def test_the_unusable_reference_error_does_not_carry_the_coordinate() -> None:
    """The loud failure must not become a new disclosure route.

    A backend's own exception routinely embeds the region it was asked for, and
    chaining would keep that message and frame visible even behind a clean one.
    The whole traceback is rendered here, because that is what a reader sees.
    """
    position = 100 + 5
    with pytest.raises(ReferenceUnusableError) as excinfo:
        _left_shift(CONTIG, position, SHIFTED_REF, SHIFTED_ALT, RaisingReference())
    rendered = "".join(
        traceback.format_exception(type(excinfo.value), excinfo.value, excinfo.value.__traceback__)
    )
    assert "chr21" not in rendered, "the contig reached the traceback"
    assert str(position) not in rendered, "the position reached the traceback"
    assert "htslib read failed" not in rendered, "the backend's own message was chained through"
    # Still debuggable: the failing accessor's exception type and a correlation
    # handle for the locus survive.
    assert "OSError" in rendered
    assert "<locus:" in rendered


def test_the_unusable_reference_error_names_the_problem_not_the_value() -> None:
    """A message that only says "reference error" moves the debugging cost, it does not remove it.

    Reasons are distinguishable, and each of them points at a different fix: a
    stale handle, a wrong assembly, a truncated FASTA, a non-nucleotide sequence.
    """
    reasons: list[str] = []
    for reference in (RaisingReference(), MissingContigReference(), EmptyReference()):
        with pytest.raises(ReferenceUnusableError) as excinfo:
            _left_shift(CONTIG, SHIFTED_POSITION, SHIFTED_REF, SHIFTED_ALT, reference)
        reasons.append(str(excinfo.value))
    assert len(set(reasons)) == len(reasons), "every failure mode reads identically"
    assert "OSError" in reasons[0]
    assert "KeyError" in reasons[1]
    assert "truncated" in reasons[2]


# ---------------------------------------------------------------------------
# Running out of shift budget is not the same as succeeding
# ---------------------------------------------------------------------------
#
# Both shift loops ran `for _ in range(MAX_SHIFT_BP)` and, on falling out of the
# loop still moving, returned as if they had finished. `canonicalise_allele`
# labelled that `ReferenceStatus.USABLE` and `rightmost_equivalent_bound` returned
# a `QueryBound` marked `proven`. A repeat tract longer than 1 kb therefore
# produced a partial join key stamped as reference-backed, `unaligned_indel_count`
# never counted it, and the batch report went on to say that every indel had been
# placed against the reference.
#
# The state is real: >1 kb homopolymer and short-tandem-repeat tracts exist in
# GRCh38, and they are precisely the loci where a mis-placed indel is most likely
# to be scored novel and ultra-rare.


class UnboundedHomopolymerReference:
    """A contig that is nothing but `A`. Every read succeeds; the tract never ends.

    Deliberately NOT a broken reference: it answers every question it is asked,
    with a valid nucleotide, for ever. It is the shift *budget* that runs out, not
    the FASTA, which is why the resulting status must not be `UNUSABLE` — an
    operator told the reference is unusable would go and check a file that is
    perfectly healthy.
    """

    def __init__(self) -> None:
        self.reads = 0

    def fetch(self, contig: str, start: int, end: int) -> str:
        self.reads += 1
        return "A" * (end - start + 1)


#: A one-base insertion inside the unbounded tract. Shiftable in both directions
#: without limit, which is the whole point.
TRACT_POSITION = 500_000
TRACT_REF = "A"
TRACT_ALT = "AA"


def test_a_left_shift_that_exhausts_its_budget_says_so() -> None:
    """`_left_shift` reports the stop, so its caller can tell it from a finish."""
    reference = UnboundedHomopolymerReference()
    *_, limit_reached = _left_shift(CONTIG, TRACT_POSITION, TRACT_REF, TRACT_ALT, reference)
    assert limit_reached is True
    assert reference.reads > 0, "the tract was never walked at all"


def test_a_left_shift_that_finishes_normally_does_not_claim_the_limit() -> None:
    """The healthy path must stay healthy: every ordinary stop is a reason, not a limit."""
    reference = WorkingReference()
    *_, limit_reached = _left_shift(CONTIG, SHIFTED_POSITION, SHIFTED_REF, SHIFTED_ALT, reference)
    assert limit_reached is False


def test_exhausting_the_shift_budget_is_not_reported_as_a_proven_alignment() -> None:
    """The defect. A partial shift must not come back stamped reference-backed.

    Before the fix this asserted `USABLE` and `left_alignment_proven` was True, so
    the key looked identical to one that really had reached its left-most position.
    """
    canonical = canonicalise_allele(
        contig=CONTIG,
        position=TRACT_POSITION,
        ref=TRACT_REF,
        alt=TRACT_ALT,
        reference=UnboundedHomopolymerReference(),
    )
    assert canonical.reference_status is ReferenceStatus.SHIFT_LIMIT_REACHED
    assert not canonical.left_alignment_proven, (
        "a key that stopped mid-tract claims to be the proven left-most spelling"
    )


def test_the_shift_limit_is_a_distinct_state_from_a_broken_reference() -> None:
    """Different causes, different fixes, so they must not share a value.

    `UNUSABLE` sends an operator to the FASTA. Here the FASTA is fine and the only
    remedy is a larger `MAX_SHIFT_BP`; conflating them would send someone to
    re-verify a file that was never at fault.
    """
    limited = canonicalise_allele(
        contig=CONTIG,
        position=TRACT_POSITION,
        ref=TRACT_REF,
        alt=TRACT_ALT,
        reference=UnboundedHomopolymerReference(),
    )
    broken = canonicalise_shifted(RaisingReference())
    assert limited.reference_status is not ReferenceStatus.UNUSABLE
    assert limited.reference_status is not broken.reference_status
    assert limited.reference_status is not ReferenceStatus.NOT_SUPPLIED


def test_a_rightward_bound_that_exhausts_its_budget_is_not_proven() -> None:
    """The mirror loop had the same fall-through, with the same consequence.

    An unproven bound means the fetch window may be short, so a source record
    beyond it is never even read — a miss indistinguishable from "the source holds
    no record", which is what `QueryBound` carries its provenance to prevent.
    """
    bound = rightmost_equivalent_bound(
        contig=CONTIG,
        position=TRACT_POSITION,
        ref=TRACT_REF,
        alt=TRACT_ALT,
        reference=UnboundedHomopolymerReference(),
    )
    assert bound.reference_status is ReferenceStatus.SHIFT_LIMIT_REACHED
    assert not bound.proven, "a bound that ran out of budget claims to be the proven right-most"
    assert bound.position != TRACT_POSITION, "the search did not move at all"


def test_a_rightward_bound_that_finishes_normally_is_still_proven() -> None:
    """The healthy path is unchanged; the new state must not leak into it."""
    bound = rightmost_equivalent_bound(
        contig=CONTIG,
        position=SHIFTED_POSITION,
        ref=SHIFTED_REF,
        alt=SHIFTED_ALT,
        reference=WorkingReference(),
    )
    assert bound.reference_status is ReferenceStatus.USABLE
    assert bound.proven


def test_the_batch_report_surfaces_a_shift_limited_indel() -> None:
    """`alleles.py` must stop reporting "every indel was placed against it".

    The count is the whole point: before the fix a shift-limited indel was in
    neither `unaligned_indel_count` nor any other tally, so the batch summed to
    APPLIED and `describe()` stated that all of them had been left-aligned.
    """
    report = summarise_left_alignment(
        indel_count=3,
        shifted_count=3,
        unaligned_indel_count=0,
        shift_limited_count=1,
        reference_available=True,
    )
    assert report.status is LeftAlignmentStatus.INCOMPLETE_SHIFT_LIMIT
    assert report.is_degraded, "an indel that never reached its left-most position is not healthy"
    assert report.shift_limited_count == 1
    assert report.as_dict()["shift_limited_count"] == 1

    described = report.describe()
    assert "1 of 3" in described
    assert "shift budget" in described
    assert "not at fault" in described
    assert "Indel join keys agree with left-aligned sources" not in described


def test_an_unreadable_reference_still_outranks_a_shift_limit_in_the_batch_status() -> None:
    """Worst-first. Both counts survive in the report either way, so nothing is lost."""
    report = summarise_left_alignment(
        indel_count=4,
        shifted_count=2,
        unaligned_indel_count=1,
        shift_limited_count=1,
        reference_available=True,
    )
    assert report.status is LeftAlignmentStatus.INCOMPLETE_REFERENCE_UNUSABLE
    assert report.shift_limited_count == 1
    assert report.unaligned_indel_count == 1
