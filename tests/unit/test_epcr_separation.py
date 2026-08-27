"""EPCR ties are a scoring defect, and these tests are what stop them recurring.

The challenge scorer sees exactly two properties of the EPCR column: the ORDER of
the rows, and which rows are TIED. Rank comes from
``sorted(enumerate(rows), key=(-epcr, file_index))``; F-max sweeps thresholds over
the unique emitted values with ``row.epcr >= t`` and unions
``predicted_variants |= row.variants``. Executed against the real scorer source
(see `docs/references/track1-submission-contract.md`, re-verified 2026-08-28), a
true pair tied with a wrong pair scores F-max 0.6667 instead of 1.0000, and loses
50 rank points on top if file order puts the wrong row first.

So two things have to hold, and they pull in opposite directions:

* the emitted values must be **strictly decreasing and well separated**, and
* the **relative order must be exactly what it was before separation** — the pass
  is a formatting fix, and a formatting fix that reordered two hypotheses would
  be a scientific change smuggled in as one.

The second is the one worth testing hardest, so it is asserted against an
independently recomputed pre-fix ordering rather than against the pass's own
output.
"""

from __future__ import annotations

import csv
import io
from dataclasses import replace

import pytest

from mva.models.genome import GenomeBuild, GenomicCoordinate
from mva.models.pair import (
    CandidatePair,
    ComponentScores,
    InheritanceModel,
    PhaseEvidence,
    PhaseStatus,
    make_pair_id,
)
from mva.models.variant import (
    ConsequenceAnnotation,
    FilterStatus,
    Genotype,
    ImpactSeverity,
    VariantRecord,
    Zygosity,
)
from mva.reporting.track1 import (
    ACCEPTED_PROBAND_ID,
    EPCR_CEILING,
    EPCR_DECIMALS,
    EPCR_FLOOR,
    MAX_SUBMISSION_ROWS,
    MIN_EPCR_SEPARATION,
    SubmissionRow,
    build_submission_rows,
    composite_to_epcr,
    render_submission_csv,
    validate_submission,
)

pytestmark = pytest.mark.unit

GENE = "SYNTHKIN1"


def _variant(position: int, *, contig: str = "chr15", alt: str = "T") -> VariantRecord:
    return VariantRecord(
        coordinate=GenomicCoordinate(
            build=GenomeBuild.GRCH38, contig=contig, position=position, ref="C", alt=alt
        ),
        genotype=Genotype(
            zygosity=Zygosity.HET,
            genotype_string="0/1",
            depth=45,
            ref_reads=23,
            alt_reads=22,
            genotype_quality=99,
        ),
        filter_status=FilterStatus.PASS,
        consequences=(
            ConsequenceAnnotation(
                gene_symbol=GENE,
                transcript_id="SYNTHT0001.1",
                consequence_terms=("stop_gained",),
                impact=ImpactSeverity.HIGH,
            ),
        ),
        source_artifact="test",
    )


def _pair(
    variant_a: VariantRecord, variant_b: VariantRecord | None, *, composite: float
) -> CandidatePair:
    variant_ids = (
        (variant_a.variant_id,)
        if variant_b is None
        else (variant_a.variant_id, variant_b.variant_id)
    )
    return CandidatePair(
        pair_id=make_pair_id(GENE, variant_ids),
        gene_symbol=GENE,
        variant_a=variant_a,
        variant_b=variant_b,
        inheritance_model=(
            InheritanceModel.COMPOUND_HETEROZYGOUS
            if variant_b is not None
            else InheritanceModel.UNKNOWN
        ),
        phase=PhaseEvidence(status=PhaseStatus.UNKNOWN, method="none"),
        scores=ComponentScores(
            analytical_validity=0.9,
            rarity=0.9,
            molecular_consequence=0.9,
            inheritance_consistency=0.5,
            phenotype_similarity=0.8,
            mechanistic_relevance=0.7,
            evidence_quality=0.6,
            contradiction_penalty=0.0,
        ),
        composite_score=composite,
        recommended_next_test="Parental segregation testing.",
    )


def _colliding_pairs() -> list[CandidatePair]:
    """Ten candidates reproducing the tie structure of the golden demo submission.

    Two collisions at four decimals (0.4420 x2, 0.3900 x2) plus a third pairing
    that lands 0.0038 below its predecessor — separated by rounding, but by less
    than one hundredth.
    """
    composites = [
        0.877899,
        0.667470,
        0.617133,
        0.610899,
        0.590040,
        0.522515,
        0.436315,
        0.436315,
        0.383881,
        0.383881,
    ]
    return [
        _pair(
            _variant(40_200_000 + index * 1_000),
            _variant(40_500_000 + index * 1_000, alt="A"),
            composite=composite,
        )
        for index, composite in enumerate(composites)
    ]


def _reference_order(pairs: list[CandidatePair]) -> list[str]:
    """The submitted order the PRE-FIX code produced, recomputed independently.

    Deliberately not derived from ``build_submission_rows``: this is the baseline
    the separation pass must not disturb, so it is rebuilt from the ranked list
    the same way the pre-fix implementation did — candidate sort key, then a
    stable sort by rounded EPCR descending.
    """
    ordered = sorted(pairs, key=lambda pair: pair.sort_key())
    ordered.sort(key=lambda pair: -composite_to_epcr(pair.composite_score))
    return [pair.pair_id for pair in ordered[:MAX_SUBMISSION_ROWS]]


def _row_ids(rows: tuple[SubmissionRow, ...]) -> list[str]:
    """Row identity by coordinates, so it can be matched back to a candidate."""
    return [f"{row.chrom_1}:{row.pos_1}:{row.chrom_2}:{row.pos_2}" for row in rows]


def _pair_ids(pairs: list[CandidatePair]) -> dict[str, str]:
    keys: dict[str, str] = {}
    for pair in pairs:
        second = pair.variant_b
        first = pair.variant_a.coordinate
        key = (
            f"{first.contig}:{first.position}:"
            f"{second.coordinate.contig if second else ''}:"
            f"{second.coordinate.position if second else ''}"
        )
        keys[key] = pair.pair_id
    return keys


# ---------------------------------------------------------------------------
# The defect, and the property that fixes it.
# ---------------------------------------------------------------------------


def test_the_golden_tie_structure_no_longer_produces_ties() -> None:
    """The live defect: four tied EPCRs in the demo submission (0.4420, 0.3900).

    Without the separation pass these ten composites round to eight distinct
    values. A tie at the answer's own threshold costs up to 0.333 F-max, and up to
    a further 50 rank points if file order puts the wrong row first.
    """
    pairs = _colliding_pairs()
    unseparated = [composite_to_epcr(pair.composite_score) for pair in pairs]
    assert len(set(unseparated)) == 8, "fixture no longer reproduces the tie defect"

    rows = build_submission_rows(pairs, proband_id=ACCEPTED_PROBAND_ID)
    epcrs = [row.epcr for row in rows]
    assert len(set(epcrs)) == len(epcrs), f"emitted EPCRs still contain a tie: {epcrs}"


def test_emitted_epcrs_are_strictly_decreasing_with_the_minimum_separation() -> None:
    pairs = _colliding_pairs()
    rows = build_submission_rows(pairs, proband_id=ACCEPTED_PROBAND_ID)
    epcrs = [row.epcr for row in rows]

    assert len(epcrs) == MAX_SUBMISSION_ROWS
    for index in range(1, len(epcrs)):
        gap = epcrs[index - 1] - epcrs[index]
        assert gap >= MIN_EPCR_SEPARATION - 1e-9, (
            f"rows {index} and {index + 1} are {gap:.4f} apart, below the "
            f"{MIN_EPCR_SEPARATION} minimum; the F-max sweep cannot resolve them "
            "if a threshold grid is coarser than the gap."
        )


def test_separation_preserves_the_pre_fix_relative_order_exactly() -> None:
    """The load-bearing property: this is a formatting pass, not a ranking one.

    If separation could reorder two rows it would be changing which hypothesis we
    rank above which — a scientific claim — while presenting itself as a rendering
    detail. Asserted against an independently recomputed pre-fix ordering.
    """
    pairs = _colliding_pairs()
    expected = _reference_order(pairs)
    rows = build_submission_rows(pairs, proband_id=ACCEPTED_PROBAND_ID)

    keys = _pair_ids(pairs)
    emitted = [keys[key] for key in _row_ids(rows)]
    assert emitted == expected, (
        "the EPCR separation pass reordered the submission. It may only turn a "
        "collision into a separation; it may never change which row outranks which."
    )


def test_separation_never_raises_a_row_above_its_own_calibrated_value() -> None:
    """Downward-only, except where the floor forces a lift.

    Lowering a colliding row is conservative — it under-claims confidence. Raising
    one would over-claim it. With ten rows and 0.09 of range needed, the floor is
    never reached, so no row is ever lifted here.
    """
    pairs = _colliding_pairs()
    rows = build_submission_rows(pairs, proband_id=ACCEPTED_PROBAND_ID)
    keys = _pair_ids(pairs)
    by_id = {pair.pair_id: pair for pair in pairs}
    for row, key in zip(rows, _row_ids(rows), strict=True):
        original = composite_to_epcr(by_id[keys[key]].composite_score)
        assert row.epcr <= original + 1e-9, (
            f"row EPCR {row.epcr} exceeds its own composite-derived {original}"
        )


# ---------------------------------------------------------------------------
# The scorer's own bounds, at the hardest input.
# ---------------------------------------------------------------------------


def test_ten_zero_scored_candidates_stay_inside_the_accepted_range() -> None:
    """Every composite 0.0: separation has to invent ten values above the floor.

    The scorer validates ``0 < epcr <= 1`` and rejects zero outright, so a naive
    "subtract 0.01 each time" pass would emit a negative EPCR by row 2 and lose
    the whole submission to a ValueError.
    """
    pairs = [
        _pair(_variant(40_200_000 + index * 1_000), None, composite=0.0) for index in range(10)
    ]
    rows = build_submission_rows(pairs, proband_id=ACCEPTED_PROBAND_ID)

    epcrs = [row.epcr for row in rows]
    assert len(rows) == 10
    assert len(set(epcrs)) == 10
    assert all(EPCR_FLOOR <= value <= EPCR_CEILING for value in epcrs), epcrs
    assert all(0.0 < value <= 1.0 for value in epcrs)
    assert epcrs == sorted(epcrs, reverse=True)
    assert epcrs[-1] == pytest.approx(EPCR_FLOOR)


def test_ten_perfect_candidates_stay_inside_the_accepted_range() -> None:
    pairs = [
        _pair(_variant(40_200_000 + index * 1_000), None, composite=1.0) for index in range(10)
    ]
    epcrs = [row.epcr for row in build_submission_rows(pairs, proband_id=ACCEPTED_PROBAND_ID)]
    assert epcrs[0] == pytest.approx(EPCR_CEILING)
    assert len(set(epcrs)) == 10
    assert all(0.0 < value <= 1.0 for value in epcrs)


# ---------------------------------------------------------------------------
# The rendered artifact, and the self-check that guards it.
# ---------------------------------------------------------------------------


def test_rendered_csv_carries_ten_distinct_epcr_strings() -> None:
    """Distinct floats are not enough — the scorer parses the rendered text."""
    rows = build_submission_rows(_colliding_pairs(), proband_id=ACCEPTED_PROBAND_ID)
    text = render_submission_csv(rows)
    rendered = [row["epcr"] for row in csv.DictReader(io.StringIO(text, newline=""))]
    assert len(set(rendered)) == len(rendered), rendered
    assert all(len(value.split(".")[1]) == EPCR_DECIMALS for value in rendered)

    ok, errors = validate_submission(text)
    assert ok, errors


def test_validate_submission_rejects_a_repeated_epcr() -> None:
    """The contract self-check must catch a tie even if it arrives from elsewhere."""
    rows = build_submission_rows(_colliding_pairs(), proband_id=ACCEPTED_PROBAND_ID)
    tied = (*rows[:5], replace(rows[5], epcr=rows[4].epcr), *rows[6:])
    ok, errors = validate_submission(render_submission_csv(tied))
    assert not ok
    assert any("repeats the value" in error for error in errors), errors


def test_separation_is_deterministic(  # GP-30
) -> None:
    """Byte-identical repeat runs. Integer EPCR arithmetic, not repeated float subtraction."""
    pairs = _colliding_pairs()
    first = render_submission_csv(build_submission_rows(pairs, proband_id=ACCEPTED_PROBAND_ID))
    second = render_submission_csv(
        build_submission_rows(list(reversed(pairs)), proband_id=ACCEPTED_PROBAND_ID)
    )
    assert first == second, "submission bytes depend on the order candidates were handed in"


def test_a_single_row_submission_is_left_alone() -> None:
    rows = build_submission_rows(
        [_pair(_variant(40_200_000), _variant(40_210_500, alt="A"), composite=0.5)],
        proband_id=ACCEPTED_PROBAND_ID,
    )
    assert len(rows) == 1
    assert rows[0].epcr == pytest.approx(composite_to_epcr(0.5))
