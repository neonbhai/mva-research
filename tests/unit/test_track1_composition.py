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
from hypothesis import given, settings
from hypothesis import strategies as st

from mva.models.genome import GenomeBuild, GenomicCoordinate
from mva.models.pair import (
    NO_SECOND_VARIANT,
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
    render_submission_csv_unvalidated,
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


def _assert_no_subset_precedes_its_superset(rows: tuple[SubmissionRow, ...]) -> None:
    """No emitted row may be a strict subset of a row ranked below it (ADR 0015/0023)."""
    keys = [row.variant_keys() for row in rows]
    for index, variants in enumerate(keys):
        below = [other for other in range(index + 1, len(keys)) if variants < keys[other]]
        assert not below, (
            f"row {index + 1} is a strict subset of row(s) {[i + 1 for i in below]}, which "
            "are ranked below it; only the superset can be a full match and it re-emits "
            "the subset's variants, so the subset above it is a rank demotion for nothing"
        )


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
# The exchange must leave the window in the state pass 1 put it in.
#
# Promotion runs over the WHOLE ranked list and then exchanges an out-of-window
# superset for an in-window subset (ADR 0023). The exchanged-in pair lands at the
# subset's slot — which can be BELOW a different subset of that same pair, one the
# exchange skipped because another superset already covered it inside the window.
# That re-creates exactly the inversion the first pass exists to remove, and the
# first pass was not run again.
# ---------------------------------------------------------------------------


def test_the_exchange_does_not_re_create_the_inversion_promotion_removes() -> None:
    """DEFECT 1. A pair swapped in from below the cut must not land under a subset.

    Five candidates into a four-row window:

    * ``{v1}`` — a single that outscores its own parent pair, so it is kept;
    * ``{v1, v2}`` — that parent pair, lifted above it by pass 1;
    * ``{v3}`` — a second single, with no superset inside the window;
    * ``{v4, v5}`` — an unrelated hypothesis;
    * ``{v1, v3}`` — ranked below the cut, and a superset of BOTH singles.

    Pass 2 exchanges ``{v3}`` out for ``{v1, v3}``, which is the whole point of
    ADR 0023. But ``{v1, v3}`` lands at slot 3, and ``{v1}`` is sitting at slot 2 —
    so the file now proposes a strict subset above its own superset, which is the
    arrangement pass 1 removed two passes earlier. If ``{v1, v3}`` is the true
    compound heterozygote, that costs it a rank for nothing: only the pair can
    satisfy ``row.variants == true_variants``, and it re-emits ``{v1}`` anyway, so
    the F-max union is identical either way.
    """
    v1 = _variant(40_200_000)
    v2 = _variant(40_210_500, alt="A")
    v3 = _variant(40_300_000, alt="G")
    v4 = _variant(40_400_000, alt="A")
    v5 = _variant(40_500_000, alt="G")

    candidates = [
        _candidate(v1, None, composite=0.95),
        _candidate(v1, v2, composite=0.90),
        _candidate(v3, None, composite=0.85),
        _candidate(v4, v5, composite=0.80),
        _candidate(v1, v3, composite=0.40),
    ]

    rows = build_submission_rows(candidates, proband_id=ACCEPTED_PROBAND_ID, max_rows=4)

    assert [_shape(row) for row in rows] == [
        (40_200_000, 40_210_500),
        (40_200_000, 40_300_000),
        (40_200_000, None),
        (40_400_000, 40_500_000),
    ], (
        "the pair exchanged in from below the cut is ranked beneath a single-variant "
        "row it subsumes; promotion was not re-run after the exchange moved it in"
    )
    _assert_no_subset_precedes_its_superset(rows)
    _assert_strictly_separated(rows)


def test_the_renderer_refuses_a_row_that_is_a_strict_subset_of_a_later_row() -> None:
    """DEFECT 1, closing assertion. The renderer must not trust the passes.

    The exact three-row file the exchange produced, hand-built so the check is on
    the renderer rather than on composition: row 2 proposes ``{chr1:100}`` and row
    3 proposes ``{chr1:100, chr2:200}``, so a strict subset is ranked above its own
    superset. Every per-row rule holds — right proband, chr-prefixed contigs,
    distinct in-range epcrs in descending order, no duplicated variant set, no
    compound-het hypothesis split across two rows — and `render_submission_csv`
    accepted it. Three composition passes vouching for an invariant is not the same
    as checking it, which is how the exchange re-introduced this silently.
    """
    rows = (
        SubmissionRow(
            "PROBAND01", "chr1", 100, "A", "T", "chr3", "300", "G", "A", 0.70, "primary", "GENA"
        ),
        SubmissionRow("PROBAND01", "chr1", 100, "A", "T", "", "", "", "", 0.69, "primary", "GENA"),
        SubmissionRow(
            "PROBAND01", "chr1", 100, "A", "T", "chr2", "200", "C", "G", 0.60, "primary", "GENB"
        ),
    )

    ok, errors = validate_submission(render_submission_csv_unvalidated(rows))
    assert not ok
    assert any("row 2 is a strict subset of row 3" in error for error in errors), errors

    with pytest.raises(ValueError, match="row 2 is a strict subset of row 3"):
        render_submission_csv(rows)


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


# ---------------------------------------------------------------------------
# The sort key is TOTAL, so input order decides nothing (GP-30).
#
# `CandidatePair.sort_key()` used to carry only (-composite, variant_a). Two
# pairs with the same score and the same FIRST variant therefore compared equal,
# every sort over them was merely *stable*, and the emitted order was whichever
# order the caller happened to supply. That is not a tidiness bug: whichever row
# landed first took the higher EPCR, and rank points go to the best FULL match
# ordered by EPCR descending — so if the demoted row was the true answer, the
# submission scored 50 rank points instead of 100 because of list order.
# ---------------------------------------------------------------------------


def test_composition_is_deterministic_for_tied_pairs_sharing_first_variant() -> None:
    """The exact shape the old two-component key could not separate.

    Same composite, same ``variant_a``, different ``variant_b``. Reversing the
    input must not change one byte of the CSV, and the tiebreak must fall to the
    second variant's coordinate — 40210500 before 40211000.
    """
    shared = _variant(40_200_000)
    left = _candidate(shared, _variant(40_210_500, alt="A"), composite=0.80)
    right = _candidate(shared, _variant(40_211_000, alt="G"), composite=0.80)

    assert left.sort_key() != right.sort_key(), (
        "two distinct candidates share a sort key; every sort over them is only "
        "as deterministic as the caller's input order"
    )

    forward = render_submission_csv(
        build_submission_rows([left, right], proband_id=ACCEPTED_PROBAND_ID)
    )
    reversed_ = render_submission_csv(
        build_submission_rows([right, left], proband_id=ACCEPTED_PROBAND_ID)
    )

    assert forward == reversed_, (
        "reversing two tied candidates changed the rendered submission; the "
        "second variant is not in the sort key"
    )
    rows = build_submission_rows([right, left], proband_id=ACCEPTED_PROBAND_ID)
    assert [_shape(row) for row in rows] == [
        (40_200_000, 40_210_500),
        (40_200_000, 40_211_000),
    ], "the tie must break on the second variant's coordinate, ascending"
    assert rows[0].epcr > rows[1].epcr


def test_a_single_and_a_pair_sharing_a_first_variant_order_deterministically() -> None:
    """The sentinel case: one candidate has no second variant at all.

    :data:`NO_SECOND_VARIANT` sorts below every real coordinate, so the single
    precedes the pair *by construction* rather than by whichever the caller
    listed first. (Composition then promotes the pair above it — ADR 0015 — which
    is a separate, deliberate rule; what matters here is that the input order
    cannot change the outcome.)
    """
    shared = _variant(40_200_000)
    single = _candidate(shared, None, composite=0.80)
    pair = _candidate(shared, _variant(40_210_500, alt="A"), composite=0.80)

    assert single.sort_key()[2] == NO_SECOND_VARIANT
    assert sorted([pair, single], key=lambda c: c.sort_key()) == [single, pair]
    assert sorted([single, pair], key=lambda c: c.sort_key()) == [single, pair]

    forward = render_submission_csv(
        build_submission_rows([single, pair], proband_id=ACCEPTED_PROBAND_ID)
    )
    reversed_ = render_submission_csv(
        build_submission_rows([pair, single], proband_id=ACCEPTED_PROBAND_ID)
    )
    assert forward == reversed_


def test_a_genuine_score_difference_still_decides_the_order() -> None:
    """Requirement: this only breaks ties. Composite remains the primary key.

    Two candidates that already differ on score must keep their existing relative
    order whatever the appended tiebreakers would have said on their own — here
    the lower-scoring candidate would win on both coordinate and pair id.
    """
    high = _candidate(_variant(40_900_000), _variant(40_901_000, alt="A"), composite=0.90)
    low = _candidate(_variant(40_100_000), _variant(40_100_500, alt="A"), composite=0.10)

    assert sorted([low, high], key=lambda c: c.sort_key()) == [high, low]
    assert high.sort_key()[0] < low.sort_key()[0], "composite must lead, descending"
    # The tiebreakers point the other way; they must not get a vote.
    assert high.sort_key()[1] > low.sort_key()[1]


def _tie_shaped_candidates() -> list[CandidatePair]:
    """A list containing every ordering hazard at once.

    Three candidates tied at 0.80 sharing a first variant (one of them a single),
    two more tied at 0.50 sharing a different first variant, and two untied rows
    on either side of them.
    """
    shared = _variant(40_200_000)
    other = _variant(40_300_000)
    return [
        _candidate(_variant(40_100_000), _variant(40_100_500, alt="A"), composite=0.95),
        _candidate(shared, None, composite=0.80),
        _candidate(shared, _variant(40_210_500, alt="A"), composite=0.80),
        _candidate(shared, _variant(40_211_000, alt="G"), composite=0.80),
        _candidate(other, _variant(40_310_000, alt="A"), composite=0.50),
        _candidate(other, _variant(40_311_000, alt="G"), composite=0.50),
        _candidate(_variant(40_900_000), _variant(40_901_000, alt="A"), composite=0.20),
    ]


@given(order=st.permutations(range(7)))
@settings(max_examples=100, deadline=None)
def test_any_permutation_of_the_candidate_list_renders_the_same_submission(
    order: list[int],
) -> None:
    """GP-30 as a property: the caller's list order carries no information.

    The candidate list is a *set* of hypotheses that happens to arrive in a list.
    If any permutation of it renders a different CSV, then some ordering decision
    is being made by the order of a Python list rather than by the science.
    """
    candidates = _tie_shaped_candidates()
    baseline = render_submission_csv(
        build_submission_rows(candidates, proband_id=ACCEPTED_PROBAND_ID)
    )
    permuted = render_submission_csv(
        build_submission_rows([candidates[i] for i in order], proband_id=ACCEPTED_PROBAND_ID)
    )
    assert permuted == baseline, f"permutation {order} rendered a different submission"


def test_every_candidate_in_a_tied_block_has_a_distinct_sort_key() -> None:
    """Totality, asserted directly rather than inferred from a rendered file."""
    keys = [candidate.sort_key() for candidate in _tie_shaped_candidates()]
    assert len(set(keys)) == len(keys), "sort_key() is not injective over these candidates"
