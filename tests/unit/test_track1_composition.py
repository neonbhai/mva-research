"""How the ten rows are laid out, once the ranking itself is settled.

Ranking is a scientific judgement made in `mva.prioritization`. *Composition* is
the separate question of how to lay that judgement out in a twelve-column CSV so
the challenge scorer reads it the way we mean it. Two composition rules exist,
and they are mirror images:

* a single-variant row that falls BELOW the pair it was carved out of is dropped
  (`_drop_subsumed`) — it adds no variant the F-max union does not already hold,
  and rank points go to the first full match, so it only burns a row;
* a single-variant row that outranks its own parent pair is NOT dropped. The pair
  is moved above it instead.

The second rule is a bet, and the bet is one-sided. Rank points come from the
best *full* match, matched by frozenset equality: if the answer is the pair,
pair-above-single scores 100 and single-above-pair scores 50, and no arrangement
gets both to rank 1. F-max is identical either way, because the pair re-emits the
single's variant. The challenge states the answer is a clinically validated
compound heterozygote and its scorer gates partial credit on
`len(true_variants) == 2`, so we bet on the pair — and keep the single directly
beneath it, because a row below the answer costs nothing.
"""

from __future__ import annotations

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
    MIN_EPCR_SEPARATION,
    SubmissionRow,
    build_submission_rows,
    render_submission_csv,
    validate_submission,
)

pytestmark = pytest.mark.unit

GENE = "SYNTHKIN1"


def _variant(position: int, *, alt: str = "T") -> VariantRecord:
    return VariantRecord(
        coordinate=GenomicCoordinate(
            build=GenomeBuild.GRCH38, contig="chr15", position=position, ref="C", alt=alt
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


def _candidate(
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


def _shape(row: SubmissionRow) -> tuple[int, int | None]:
    """A row as (first position, second position or None) — enough to identify it."""
    return (row.pos_1, int(row.pos_2) if row.pos_2 else None)


def _assert_strictly_separated(rows: tuple[SubmissionRow, ...]) -> None:
    epcrs = [row.epcr for row in rows]
    assert len(set(epcrs)) == len(epcrs), f"tied EPCRs after composition: {epcrs}"
    for index in range(1, len(epcrs)):
        gap = epcrs[index - 1] - epcrs[index]
        assert gap >= MIN_EPCR_SEPARATION - 1e-9, (
            f"rows {index} and {index + 1} are only {gap:.4f} apart"
        )


# ---------------------------------------------------------------------------
# The rule: the pair outranks the single, and the single stays.
# ---------------------------------------------------------------------------


def test_a_pair_outranks_a_single_variant_row_it_subsumes() -> None:
    """The hazard: the single scores higher, so it used to be emitted first.

    That put the full-match pair at rank 2 — 50 rank points instead of 100 — for
    no F-max gain whatsoever, since the pair re-emits the single's variant.
    """
    first = _variant(40_200_000)
    second = _variant(40_210_500, alt="A")
    single = _candidate(first, None, composite=0.90)
    pair = _candidate(first, second, composite=0.80)

    rows = build_submission_rows([single, pair], proband_id=ACCEPTED_PROBAND_ID)

    assert [_shape(row) for row in rows] == [(40_200_000, 40_210_500), (40_200_000, None)], (
        "the single-variant row still precedes the pair whose variant it carries; "
        "the pair can never reach rank 1 from there"
    )
    _assert_strictly_separated(rows)


def test_the_outranked_single_is_kept_not_dropped() -> None:
    """Rows below the answer are free upside; deleting one is a pure forfeit.

    Dropping the single would give up a full match in the world where the truth
    is a single variant, and gain nothing in the world where it is a pair.
    """
    first = _variant(40_200_000)
    single = _candidate(first, None, composite=0.90)
    pair = _candidate(first, _variant(40_210_500, alt="A"), composite=0.80)

    rows = build_submission_rows([single, pair], proband_id=ACCEPTED_PROBAND_ID)

    assert len(rows) == 2
    assert (40_200_000, None) in [_shape(row) for row in rows], (
        "the single-variant hypothesis was deleted rather than demoted"
    )


def test_two_pairs_carrying_the_same_variant_both_outrank_the_single() -> None:
    shared = _variant(40_200_000)
    single = _candidate(shared, None, composite=0.90)
    left = _candidate(shared, _variant(40_210_500, alt="A"), composite=0.80)
    right = _candidate(shared, _variant(40_211_000, alt="G"), composite=0.70)

    rows = build_submission_rows([single, left, right], proband_id=ACCEPTED_PROBAND_ID)

    shapes = [_shape(row) for row in rows]
    assert shapes.index((40_200_000, None)) == len(shapes) - 1, (
        f"the single must sit below every pair that carries its variant: {shapes}"
    )
    assert len(rows) == 3
    _assert_strictly_separated(rows)


# ---------------------------------------------------------------------------
# The mirror rule, and everything else, must be left alone.
# ---------------------------------------------------------------------------


def test_a_single_ranked_below_its_pair_is_still_dropped() -> None:
    """`_drop_subsumed` is unchanged: only the inverted case is a reorder."""
    first = _variant(40_200_000)
    pair = _candidate(first, _variant(40_210_500, alt="A"), composite=0.90)
    single = _candidate(first, None, composite=0.80)

    rows = build_submission_rows([pair, single], proband_id=ACCEPTED_PROBAND_ID)

    assert [_shape(row) for row in rows] == [(40_200_000, 40_210_500)], (
        "a single already covered by a higher-ranked pair adds no variant the "
        "F-max union lacks; it only burns one of the ten rows"
    )


def test_unrelated_rows_keep_their_epcr_order() -> None:
    """Promotion is surgical: a row with no subset relation must not move."""
    top = _candidate(_variant(40_100_000), _variant(40_100_500, alt="A"), composite=0.95)
    shared = _variant(40_200_000)
    single = _candidate(shared, None, composite=0.90)
    pair = _candidate(shared, _variant(40_210_500, alt="A"), composite=0.80)
    tail = _candidate(_variant(40_900_000), _variant(40_901_000, alt="A"), composite=0.40)

    rows = build_submission_rows([top, single, pair, tail], proband_id=ACCEPTED_PROBAND_ID)

    assert [_shape(row) for row in rows] == [
        (40_100_000, 40_100_500),
        (40_200_000, 40_210_500),
        (40_200_000, None),
        (40_900_000, 40_901_000),
    ]
    _assert_strictly_separated(rows)


def test_a_partially_overlapping_pair_is_not_promoted() -> None:
    """Only a STRICT subset triggers promotion; two pairs sharing one variant are
    genuinely different hypotheses and neither outranks the other by this rule."""
    shared = _variant(40_200_000)
    left = _candidate(shared, _variant(40_210_500, alt="A"), composite=0.90)
    right = _candidate(shared, _variant(40_211_000, alt="G"), composite=0.80)

    rows = build_submission_rows([left, right], proband_id=ACCEPTED_PROBAND_ID)
    assert [_shape(row) for row in rows] == [
        (40_200_000, 40_210_500),
        (40_200_000, 40_211_000),
    ]


# ---------------------------------------------------------------------------
# The composed file still satisfies the contract, deterministically.
# ---------------------------------------------------------------------------


def test_the_composed_submission_passes_its_own_contract_check() -> None:
    shared = _variant(40_200_000)
    candidates = [
        _candidate(shared, None, composite=0.90),
        _candidate(shared, _variant(40_210_500, alt="A"), composite=0.80),
        _candidate(_variant(40_300_000), _variant(40_301_000, alt="A"), composite=0.60),
    ]
    rows = build_submission_rows(candidates, proband_id=ACCEPTED_PROBAND_ID)
    ok, errors = validate_submission(render_submission_csv(rows))
    assert ok, errors


def test_composition_is_deterministic_under_input_order() -> None:  # GP-30
    shared = _variant(40_200_000)
    candidates = [
        _candidate(shared, None, composite=0.90),
        _candidate(shared, _variant(40_210_500, alt="A"), composite=0.80),
        _candidate(shared, _variant(40_211_000, alt="G"), composite=0.70),
        _candidate(_variant(40_300_000), _variant(40_301_000, alt="A"), composite=0.60),
    ]
    forward = render_submission_csv(
        build_submission_rows(candidates, proband_id=ACCEPTED_PROBAND_ID)
    )
    backward = render_submission_csv(
        build_submission_rows(list(reversed(candidates)), proband_id=ACCEPTED_PROBAND_ID)
    )
    assert forward == backward
