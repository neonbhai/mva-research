"""Unit tests for the variant prioritisation stage.

The behaviours locked here are the ones a plausible-looking refactor breaks
silently: that a common variant is *down-ranked* rather than deleted (GP-13),
that unknown phase survives instead of being optimistically upgraded (GP-15),
and that two runs over the same input produce the same list in the same order
(GP-30).

Fixtures are built inline by :func:`make_variant` rather than read from disk, so
each test states the exact genotype it depends on.
"""

from __future__ import annotations

from datetime import UTC, datetime

import pytest
from hypothesis import given
from hypothesis import strategies as st

from mva.clock import FixedClock
from mva.config import FrequencyThresholds, PhaseWeights, QualityThresholds, ScoringWeights
from mva.models.genome import GenomeBuild, GenomicCoordinate
from mva.models.pair import CandidatePair, ComponentScores, InheritanceModel, PhaseStatus
from mva.models.variant import (
    ConsequenceAnnotation,
    FilterStatus,
    Genotype,
    ImpactSeverity,
    PopulationFrequency,
    VariantRecord,
    Zygosity,
)
from mva.prioritization.filters import (
    FLAG_COMMON_VARIANT,
    FLAG_LOW_QUALITY_CALL,
    FLAG_PLAUSIBLE_CANDIDATE,
    REASON_NO_ALT_ALLELE,
    REASON_NON_CANONICAL_CONTIG,
    REASON_WRONG_GENOME_BUILD,
    apply_hard_filters,
    apply_soft_flags,
    select_candidate_variants,
)
from mva.prioritization.pairing import PairCandidate, generate_pairs, infer_phase
from mva.prioritization.ranking import rank_pairs
from mva.prioritization.scoring import (
    NEUTRAL_MECHANISM_SCORE,
    NEUTRAL_PHENOTYPE_SCORE,
    ScoredPair,
    composite_score,
    score_pair,
)

CLOCK = FixedClock(datetime(2026, 1, 1, tzinfo=UTC))
FREQUENCY = FrequencyThresholds()
QUALITY = QualityThresholds()
WEIGHTS = ScoringWeights()
PHASE_WEIGHTS = PhaseWeights()

#: The synthetic gene the fixtures below hang off. Fictional, per GP-20.
GENE = "SYNTHKIN1"


def make_variant(
    *,
    contig: str = "chr15",
    position: int = 40200000,
    ref: str = "C",
    alt: str = "T",
    build: GenomeBuild = GenomeBuild.GRCH38,
    zygosity: Zygosity = Zygosity.HET,
    genotype_string: str = "0/1",
    phased: bool = False,
    phase_set: int | None = None,
    depth: int | None = 45,
    ref_reads: int | None = 23,
    alt_reads: int | None = 22,
    genotype_quality: int | None = 99,
    filter_status: FilterStatus = FilterStatus.PASS,
    gene: str = GENE,
    impact: ImpactSeverity = ImpactSeverity.HIGH,
    consequence_terms: tuple[str, ...] = ("stop_gained",),
    splice_ai: float | None = 0.02,
    pathogenicity_scores: dict[str, float] | None = None,
    allele_frequency: float | None = 0.0,
    population: str = "global",
    homozygote_count: int = 0,
    annotate: bool = True,
    qc_flags: tuple[str, ...] = (),
) -> VariantRecord:
    """Build one fully-specified ``VariantRecord``.

    Defaults describe a clean, ultra-rare, HIGH-impact heterozygote — the shape
    a real candidate has — so each test only has to say how it differs.
    """
    frequencies: tuple[PopulationFrequency, ...] = ()
    if allele_frequency is not None:
        frequencies = (
            PopulationFrequency(
                source="SYNTH_gnomAD",
                version="v0.0-synthetic",
                population=population,
                allele_frequency=allele_frequency,
                allele_count=0,
                allele_number=152312,
                homozygote_count=homozygote_count,
                filter_status="PASS",
            ),
        )
    consequences: tuple[ConsequenceAnnotation, ...] = ()
    if annotate:
        consequences = (
            ConsequenceAnnotation(
                gene_symbol=gene,
                gene_id="SYNTHG0001",
                transcript_id="SYNTHT0001.1",
                is_canonical=True,
                is_mane_select=True,
                consequence_terms=consequence_terms,
                impact=impact,
                splice_ai_delta_max=splice_ai,
                pathogenicity_scores=pathogenicity_scores or {"CADD_phred": 35.0},
                source_tool="SYNTH_vep",
                source_tool_version="v0.0-synthetic",
            ),
        )
    return VariantRecord(
        coordinate=GenomicCoordinate(
            build=build, contig=contig, position=position, ref=ref, alt=alt
        ),
        genotype=Genotype(
            zygosity=zygosity,
            genotype_string=genotype_string,
            phased=phased,
            phase_set=phase_set,
            depth=depth,
            ref_reads=ref_reads,
            alt_reads=alt_reads,
            genotype_quality=genotype_quality,
        ),
        filter_status=filter_status,
        raw_filters=("PASS",) if filter_status is FilterStatus.PASS else ("LowQual",),
        quality=820.5,
        consequences=consequences,
        population_frequencies=frequencies,
        qc_flags=qc_flags,
        source_artifact="inline_test_fixture",
    )


def score_all(candidates: tuple[PairCandidate, ...]) -> list[ScoredPair]:
    return [
        score_pair(
            candidate,
            phenotype_score=NEUTRAL_PHENOTYPE_SCORE,
            mechanism_score=NEUTRAL_MECHANISM_SCORE,
            weights=WEIGHTS,
            phase_weights=PHASE_WEIGHTS,
            frequency=FREQUENCY,
            quality=QUALITY,
            clock=CLOCK,
        )
        for candidate in candidates
    ]


def run_pipeline(variants: list[VariantRecord]) -> tuple[CandidatePair, ...]:
    """Soft-flag, pair, score and rank — the stage as the composition root uses it."""
    flagged = apply_soft_flags(variants, frequency=FREQUENCY, quality=QUALITY)
    return rank_pairs(score_all(generate_pairs(flagged)), clock=CLOCK)


# ---------------------------------------------------------------------------
# Fixture variants shared by several tests
# ---------------------------------------------------------------------------


def rare_stop_gained() -> VariantRecord:
    return make_variant(position=40200000, ref="C", alt="T", allele_frequency=0.0)


def rare_splice_donor() -> VariantRecord:
    return make_variant(
        position=40210500,
        ref="G",
        alt="A",
        allele_frequency=5.3e-6,
        consequence_terms=("splice_donor_variant", "missense_variant"),
        splice_ai=0.91,
        pathogenicity_scores={"CADD_phred": 32.6, "REVEL": 0.86},
        depth=38,
        ref_reads=18,
        alt_reads=20,
        genotype_quality=98,
    )


def common_synonymous() -> VariantRecord:
    return make_variant(
        position=40205000,
        ref="A",
        alt="G",
        allele_frequency=0.1204,
        homozygote_count=1102,
        impact=ImpactSeverity.LOW,
        consequence_terms=("synonymous_variant",),
        splice_ai=0.0,
        pathogenicity_scores={"CADD_phred": 1.2},
        depth=52,
        ref_reads=26,
        alt_reads=26,
    )


# ---------------------------------------------------------------------------
# 1. Recovery of the rare compound-heterozygous pair
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_rare_compound_het_pair_is_recovered_from_a_mixed_variant_set() -> None:
    variants = [
        rare_stop_gained(),
        common_synonymous(),
        rare_splice_donor(),
        make_variant(
            contig="chr7",
            position=55000000,
            ref="C",
            alt="G",
            gene="SYNTHMET2",
            impact=ImpactSeverity.MODERATE,
            consequence_terms=("missense_variant",),
            allele_frequency=2.62e-5,
            pathogenicity_scores={"CADD_phred": 26.1, "REVEL": 0.62},
        ),
    ]
    ranked = run_pipeline(variants)

    top = ranked[0]
    assert top.rank == 1
    assert top.gene_symbol == GENE
    assert top.inheritance_model is InheritanceModel.COMPOUND_HETEROZYGOUS
    assert top.variant_ids == (
        "GRCh38:chr15:40200000:C:T",
        "GRCh38:chr15:40210500:G:A",
    )
    # The ranking is a research plan, not a league table.
    assert top.recommended_next_test
    assert top.discriminating_experiment


# ---------------------------------------------------------------------------
# 2. Common variants are down-ranked, never removed (GP-13)
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_common_variant_is_down_ranked_but_still_present() -> None:
    ranked = run_pipeline([rare_stop_gained(), rare_splice_donor(), common_synonymous()])

    by_ids = {pair.variant_ids: pair for pair in ranked}
    target = by_ids[("GRCh38:chr15:40200000:C:T", "GRCh38:chr15:40210500:G:A")]
    with_common = by_ids[("GRCh38:chr15:40200000:C:T", "GRCh38:chr15:40205000:A:G")]

    assert with_common.composite_score < target.composite_score
    assert FLAG_COMMON_VARIANT in with_common.flags
    assert with_common.rank is not None and with_common.rank > 1
    # Down-ranked, and the reason is recorded rather than merely implied.
    assert with_common.contradicting_evidence_ids


@pytest.mark.unit
def test_select_candidate_variants_marks_without_deleting() -> None:
    variants = apply_soft_flags(
        [rare_stop_gained(), common_synonymous()], frequency=FREQUENCY, quality=QUALITY
    )
    selected = select_candidate_variants(variants, frequency=FREQUENCY)

    assert [v.variant_id for v in selected] == ["GRCh38:chr15:40200000:C:T"]
    assert FLAG_PLAUSIBLE_CANDIDATE in selected[0].qc_flags
    # The caller still holds everything; selection removed nothing from it.
    assert len(variants) == 2


# ---------------------------------------------------------------------------
# 3. Low-quality calls are marked and retained
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_low_quality_call_is_flagged_and_retained() -> None:
    low_quality = make_variant(
        position=40211000,
        ref="G",
        alt="GAT",
        filter_status=FilterStatus.FILTERED,
        depth=4,
        ref_reads=3,
        alt_reads=1,
        genotype_quality=8,
        consequence_terms=("frameshift_variant",),
        allele_frequency=6.6e-6,
    )
    result = apply_hard_filters([low_quality], expected_build=GenomeBuild.GRCH38)
    assert result.removed == ()
    assert len(result.retained) == 1

    flagged = apply_soft_flags(result.flagged, frequency=FREQUENCY, quality=QUALITY)
    assert FLAG_LOW_QUALITY_CALL in flagged[0].qc_flags

    ranked = run_pipeline([rare_stop_gained(), rare_splice_donor(), low_quality])
    ids = {pair.variant_ids for pair in ranked}
    assert ("GRCh38:chr15:40200000:C:T", "GRCh38:chr15:40211000:G:GAT") in ids

    with_low_quality = next(
        pair for pair in ranked if "GRCh38:chr15:40211000:G:GAT" in pair.variant_ids
    )
    assert with_low_quality.scores.analytical_validity < 0.2
    assert FLAG_LOW_QUALITY_CALL in with_low_quality.flags
    assert any("real" in question.question for question in with_low_quality.blocking_questions)


# ---------------------------------------------------------------------------
# 4. Pairs are gene-scoped
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_pairs_form_within_a_gene_and_not_across_genes() -> None:
    same_gene_a = rare_stop_gained()
    same_gene_b = rare_splice_donor()
    other_gene = make_variant(
        contig="chr7",
        position=55000000,
        ref="C",
        alt="G",
        gene="SYNTHMET2",
        impact=ImpactSeverity.MODERATE,
        consequence_terms=("missense_variant",),
    )
    candidates = generate_pairs([same_gene_a, same_gene_b, other_gene])

    two_variant = [c for c in candidates if c.is_pair]
    assert len(two_variant) == 1
    assert two_variant[0].gene_symbol == GENE
    assert two_variant[0].variant_ids == (same_gene_a.variant_id, same_gene_b.variant_id)
    assert all(len({v.coordinate.contig for v in c.variants}) == 1 for c in candidates)
    # The MODERATE-impact lone het forms no single-variant hypothesis.
    assert not any(c.gene_symbol == "SYNTHMET2" for c in candidates)


@pytest.mark.unit
def test_a_variant_in_two_genes_pairs_in_each() -> None:
    shared = make_variant(position=40200000)
    overlapping = shared.with_annotations(
        consequences=(
            *shared.consequences,
            ConsequenceAnnotation(
                gene_symbol="SYNTHOTH3",
                transcript_id="SYNTHT0010.1",
                is_canonical=True,
                consequence_terms=("stop_gained",),
                impact=ImpactSeverity.HIGH,
            ),
        )
    )
    partner = make_variant(position=40201000, ref="C", alt="A", gene="SYNTHOTH3")
    candidates = generate_pairs([overlapping, partner])
    genes = {c.gene_symbol for c in candidates if c.is_pair}
    assert genes == {"SYNTHOTH3"}
    assert {c.gene_symbol for c in candidates} == {GENE, "SYNTHOTH3"}


@pytest.mark.unit
def test_gene_pair_cap_is_recorded_when_it_truncates() -> None:
    variants = [
        make_variant(position=40200000 + offset * 1000, ref="C", alt="T") for offset in range(5)
    ]
    candidates = generate_pairs(variants, max_pairs_per_gene=3)
    assert len(candidates) == 3
    assert all("gene_pair_cap_truncated" in c.flags for c in candidates)


# ---------------------------------------------------------------------------
# 5. In-cis pairs are detected and heavily down-ranked (GP-15)
# ---------------------------------------------------------------------------


def cis_pair_variants() -> tuple[VariantRecord, VariantRecord]:
    first = make_variant(
        position=40201000,
        ref="C",
        alt="A",
        gene="SYNTHOTH3",
        genotype_string="1|0",
        phased=True,
        phase_set=40201000,
        impact=ImpactSeverity.MODERATE,
        consequence_terms=("missense_variant",),
        allele_frequency=1.31e-5,
        pathogenicity_scores={"CADD_phred": 22.4, "REVEL": 0.41},
    )
    second = make_variant(
        position=40201050,
        ref="T",
        alt="G",
        gene="SYNTHOTH3",
        genotype_string="1|0",
        phased=True,
        phase_set=40201000,
        impact=ImpactSeverity.MODERATE,
        consequence_terms=("missense_variant",),
        allele_frequency=1.97e-5,
        pathogenicity_scores={"CADD_phred": 24.9, "REVEL": 0.55},
    )
    return first, second


@pytest.mark.unit
def test_same_phase_set_same_haplotype_resolves_to_cis_confirmed() -> None:
    first, second = cis_pair_variants()
    phase = infer_phase(first, second)
    assert phase.status is PhaseStatus.CIS_CONFIRMED
    assert phase.method == "phase_set"
    assert phase.distance_bp == 50


@pytest.mark.unit
def test_opposite_haplotypes_in_one_phase_set_resolve_to_trans() -> None:
    first, second = cis_pair_variants()
    in_trans = second.model_copy(
        update={"genotype": second.genotype.model_copy(update={"genotype_string": "0|1"})}
    )
    assert infer_phase(first, in_trans).status is PhaseStatus.TRANS_CONFIRMED


@pytest.mark.unit
def test_in_cis_pair_is_disqualifying_and_does_not_rank_first() -> None:
    first, second = cis_pair_variants()
    ranked = run_pipeline([rare_stop_gained(), rare_splice_donor(), first, second])

    cis = next(pair for pair in ranked if pair.gene_symbol == "SYNTHOTH3" and pair.is_pair)
    assert cis.phase.status is PhaseStatus.CIS_CONFIRMED
    assert cis.phase_is_disqualifying is True
    assert cis.rank is not None and cis.rank > 1
    assert cis.scores.inheritance_consistency <= PHASE_WEIGHTS.cis_confirmed
    assert cis.scores.contradiction_penalty > 0.0
    assert ranked[0].gene_symbol == GENE


# ---------------------------------------------------------------------------
# 6. Unknown phase is preserved, never upgraded (GP-15)
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_unknown_phase_is_preserved_and_raises_a_blocking_question() -> None:
    first, second = rare_stop_gained(), rare_splice_donor()
    phase = infer_phase(first, second)
    assert phase.status is PhaseStatus.UNKNOWN
    assert phase.method == "none"
    assert phase.distance_bp == 10500
    assert phase.notes is not None and "read" in phase.notes

    scored = score_all(generate_pairs([first, second]))
    pair = next(item for item in scored if item.candidate.is_pair)
    blocking = [question for question in pair.open_questions if question.blocking]
    trans = next(question for question in blocking if question.question == "Is the pair in trans?")
    assert trans.resolving_test == ("parental segregation testing or long-read/read-backed phasing")
    assert pair.scores.inheritance_consistency == pytest.approx(PHASE_WEIGHTS.unknown)


@pytest.mark.unit
def test_phase_is_not_upgraded_when_phase_sets_differ() -> None:
    first, second = cis_pair_variants()
    other_set = second.model_copy(
        update={"genotype": second.genotype.model_copy(update={"phase_set": 99999999})}
    )
    assert infer_phase(first, other_set).status is PhaseStatus.UNKNOWN


@pytest.mark.unit
def test_missing_frequency_data_raises_its_own_open_question() -> None:
    unknown_frequency = make_variant(position=40210500, ref="G", alt="A", allele_frequency=None)
    flagged = apply_soft_flags(
        [rare_stop_gained(), unknown_frequency], frequency=FREQUENCY, quality=QUALITY
    )
    scored = score_all(generate_pairs(flagged))
    pair = next(item for item in scored if item.candidate.is_pair)
    assert any("allele frequency" in q.question for q in pair.open_questions)
    # Absence of data scores neutrally; it is never treated as rarity (GP-14).
    assert pair.scores.rarity == pytest.approx(FREQUENCY.absent_frequency_score)


# ---------------------------------------------------------------------------
# 7. Hard filters remove only the invalid (GP-13)
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_hard_filters_remove_invalid_records_and_nothing_else() -> None:
    wrong_build = make_variant(position=40200000, build=GenomeBuild.GRCH37)
    non_canonical = make_variant(contig="chr99", position=1000)
    hom_ref = make_variant(position=40202000, zygosity=Zygosity.HOM_REF, genotype_string="0/0")
    common = common_synonymous()
    low_quality = make_variant(
        position=40211000,
        filter_status=FilterStatus.FILTERED,
        depth=4,
        ref_reads=3,
        alt_reads=1,
        genotype_quality=8,
    )
    benign = make_variant(
        position=40203000,
        impact=ImpactSeverity.MODIFIER,
        consequence_terms=("intron_variant",),
    )

    result = apply_hard_filters(
        [wrong_build, non_canonical, hom_ref, common, low_quality, benign],
        expected_build=GenomeBuild.GRCH38,
    )

    assert dict(result.removed) == {
        wrong_build.variant_id: REASON_WRONG_GENOME_BUILD,
        non_canonical.variant_id: REASON_NON_CANONICAL_CONTIG,
        hom_ref.variant_id: REASON_NO_ALT_ALLELE,
    }
    retained_ids = {variant.variant_id for variant in result.retained}
    assert retained_ids == {common.variant_id, low_quality.variant_id, benign.variant_id}
    assert result.counts["input"] == 6
    assert result.counts["retained"] == 3
    assert result.counts[f"removed_{REASON_WRONG_GENOME_BUILD}"] == 1


@pytest.mark.unit
def test_hard_filters_never_remove_for_frequency_quality_or_consequence() -> None:
    """The three reasons GP-13 forbids, asserted explicitly so a future
    'optimisation' that adds one of them fails loudly."""
    variants = [
        common_synonymous(),
        make_variant(
            position=40211000,
            filter_status=FilterStatus.FILTERED,
            depth=2,
            ref_reads=2,
            alt_reads=0,
            genotype_quality=3,
        ),
        make_variant(
            position=40203000,
            impact=ImpactSeverity.MODIFIER,
            consequence_terms=("intron_variant",),
        ),
    ]
    result = apply_hard_filters(variants, expected_build=GenomeBuild.GRCH38)
    assert result.removed == ()
    assert len(result.retained) == len(variants)


# ---------------------------------------------------------------------------
# 8. Determinism (GP-30)
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_ranking_is_byte_stable_across_repeat_runs() -> None:
    first, second = cis_pair_variants()
    variants = [
        rare_splice_donor(),
        common_synonymous(),
        second,
        rare_stop_gained(),
        first,
    ]
    run_a = run_pipeline(variants)
    run_b = run_pipeline(list(reversed(variants)))

    assert [pair.pair_id for pair in run_a] == [pair.pair_id for pair in run_b]
    assert [pair.rank for pair in run_a] == [pair.rank for pair in run_b]
    assert [pair.composite_score for pair in run_a] == [pair.composite_score for pair in run_b]
    assert [pair.scores.as_dict() for pair in run_a] == [pair.scores.as_dict() for pair in run_b]


@pytest.mark.unit
def test_max_rank_truncates_the_returned_list_only() -> None:
    variants = [rare_stop_gained(), rare_splice_donor(), common_synonymous()]
    full = run_pipeline(variants)
    flagged = apply_soft_flags(variants, frequency=FREQUENCY, quality=QUALITY)
    limited = rank_pairs(score_all(generate_pairs(flagged)), clock=CLOCK, max_rank=2)

    assert len(full) > 2
    assert len(limited) == 2
    assert [pair.pair_id for pair in limited] == [pair.pair_id for pair in full[:2]]


# ---------------------------------------------------------------------------
# 9. Composite is monotone in every component
# ---------------------------------------------------------------------------


COMPONENT_NAMES = (
    "analytical_validity",
    "rarity",
    "molecular_consequence",
    "inheritance_consistency",
    "phenotype_similarity",
    "mechanistic_relevance",
    "evidence_quality",
)

unit_interval = st.floats(min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False)


@given(
    base=st.fixed_dictionaries(
        dict.fromkeys((*COMPONENT_NAMES, "contradiction_penalty"), unit_interval)
    ),
    component=st.sampled_from(COMPONENT_NAMES),
    increase=st.floats(min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False),
)
@pytest.mark.unit
def test_composite_is_monotone_non_decreasing_in_each_positive_component(
    base: dict[str, float], component: str, increase: float
) -> None:
    low = ComponentScores(**base)
    high = ComponentScores(**{**base, component: min(1.0, base[component] + increase)})
    assert composite_score(high, WEIGHTS) >= composite_score(low, WEIGHTS)


@given(
    base=st.fixed_dictionaries(
        dict.fromkeys((*COMPONENT_NAMES, "contradiction_penalty"), unit_interval)
    ),
    increase=st.floats(min_value=0.0, max_value=1.0, allow_nan=False, allow_infinity=False),
)
@pytest.mark.unit
def test_composite_is_monotone_non_increasing_in_the_contradiction_penalty(
    base: dict[str, float], increase: float
) -> None:
    low = ComponentScores(**base)
    high = ComponentScores(
        **{**base, "contradiction_penalty": min(1.0, base["contradiction_penalty"] + increase)}
    )
    assert composite_score(high, WEIGHTS) <= composite_score(low, WEIGHTS)


@given(
    scores=st.fixed_dictionaries(
        dict.fromkeys((*COMPONENT_NAMES, "contradiction_penalty"), unit_interval)
    )
)
@pytest.mark.unit
def test_composite_stays_inside_the_unit_interval(scores: dict[str, float]) -> None:
    assert 0.0 <= composite_score(ComponentScores(**scores), WEIGHTS) <= 1.0


# ---------------------------------------------------------------------------
# 10. Contradictions reduce the composite (GP-19)
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_contradiction_penalty_reduces_the_composite() -> None:
    neutral = dict.fromkeys(COMPONENT_NAMES, 0.8)
    clean = ComponentScores(**neutral, contradiction_penalty=0.0)
    contradicted = ComponentScores(**neutral, contradiction_penalty=0.6)

    difference = composite_score(clean, WEIGHTS) - composite_score(contradicted, WEIGHTS)
    assert difference == pytest.approx(0.6 * WEIGHTS.contradiction_penalty_weight)


@pytest.mark.unit
def test_recorded_contradictions_are_persisted_as_evidence() -> None:
    flagged = apply_soft_flags(
        [rare_stop_gained(), common_synonymous()], frequency=FREQUENCY, quality=QUALITY
    )
    scored = score_all(generate_pairs(flagged))
    pair = next(item for item in scored if item.candidate.is_pair)

    codes = {
        item.payload.get("contradiction_code")
        for item in pair.contradicting_evidence
        if "contradiction_code" in item.payload
    }
    assert {"common_variant", "population_homozygotes", "benign_consequence"} <= codes
    assert pair.scores.contradiction_penalty > 0.0
    # GP-17: every emitted item states what it does not establish.
    for item in (*pair.supporting_evidence, *pair.contradicting_evidence):
        assert "Uncalibrated pipeline heuristic" in item.limitations


@pytest.mark.unit
def test_every_component_rationale_reaches_the_final_narrative() -> None:
    scored = score_all(generate_pairs([rare_stop_gained(), rare_splice_donor()]))
    pair = next(item for item in scored if item.candidate.is_pair)
    for label in (
        "Analytical validity",
        "Rarity",
        "Molecular consequence",
        "Inheritance consistency",
        "Phenotype similarity",
        "Mechanistic relevance",
        "Evidence quality",
        "Contradiction penalty",
    ):
        assert label in pair.rationale
    assert "neutral default" in pair.rationale
    assert len(pair.supporting_evidence) + len(pair.contradicting_evidence) >= len(COMPONENT_NAMES)
