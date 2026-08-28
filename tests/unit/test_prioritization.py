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

import itertools
from datetime import UTC, datetime

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from mva.clock import FixedClock
from mva.config import FrequencyThresholds, PhaseWeights, QualityThresholds, ScoringWeights
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
    FLAG_REF_ALLELE_MISMATCH,
    INGESTION_QC_FLAGS,
    OP_SPLIT_MULTIALLELIC,
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
    FLAG_POSSIBLE_MOSAIC,
    INGESTION_NEUTRAL_FLAGS,
    INGESTION_QUALITY_FLAGS,
    REASON_NO_ALT_ALLELE,
    REASON_NON_CANONICAL_CONTIG,
    REASON_WRONG_GENOME_BUILD,
    apply_hard_filters,
    apply_soft_flags,
    select_candidate_variants,
)
from mva.prioritization.pairing import (
    PRODUCED_INHERITANCE_MODELS,
    UNPRODUCED_INHERITANCE_MODELS,
    PairCandidate,
    generate_pairs,
    infer_phase,
)
from mva.prioritization.ranking import rank_pairs
from mva.prioritization.scoring import (
    NEUTRAL_MECHANISM_SCORE,
    NEUTRAL_PHENOTYPE_SCORE,
    ScoredPair,
    _variant_analytical_validity,
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
    allele_number: int = 152312,
    homozygote_count: int = 0,
    frequencies: tuple[PopulationFrequency, ...] | None = None,
    annotate: bool = True,
    extra_consequences: tuple[ConsequenceAnnotation, ...] = (),
    qc_flags: tuple[str, ...] = (),
    normalisation_ops: tuple[str, ...] = (),
    alt_allele_index: int | None = None,
) -> VariantRecord:
    """Build one fully-specified ``VariantRecord``.

    Defaults describe a clean, ultra-rare, HIGH-impact heterozygote — the shape
    a real candidate has — so each test only has to say how it differs.
    """
    if frequencies is None:
        frequencies = ()
    if frequencies == () and allele_frequency is not None:
        frequencies = (
            PopulationFrequency(
                source="SYNTH_gnomAD",
                version="v0.0-synthetic",
                population=population,
                allele_frequency=allele_frequency,
                allele_count=0,
                allele_number=allele_number,
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
            *extra_consequences,
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
            alt_allele_index=alt_allele_index,
        ),
        filter_status=filter_status,
        raw_filters=("PASS",) if filter_status is FilterStatus.PASS else ("LowQual",),
        quality=820.5,
        consequences=consequences,
        population_frequencies=frequencies,
        qc_flags=qc_flags,
        source_artifact="inline_test_fixture",
        normalisation_ops=normalisation_ops,
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


# ---------------------------------------------------------------------------
# GP-19 regression: contradiction semantics
#
# Added by the lead engineer during integration. A low component score is weak
# support, not evidence against. Conflating the two floods the dossier's
# contradiction section with noise, which is precisely how a reader learns to
# ignore it.
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_low_component_score_is_not_a_contradiction() -> None:
    """A weak component must not be recorded as evidence AGAINST the candidate.

    Regression guard: a candidate with nothing genuinely opposing it must report
    `has_contradictions is False`, even when several of its component scores are
    low. Only `collect_contradictions` findings — in-cis phase, common allele
    frequency, population homozygotes — count as contradictions, and those are the
    only things that feed the subtracted penalty.
    """
    from mva.models.evidence import EvidenceDirection

    # A rare, high-impact, good-quality pair with nothing opposing it, but with
    # deliberately low phenotype and mechanism inputs.
    variants = (
        make_variant(position=40_200_000, ref="C", alt="T"),
        make_variant(position=40_210_500, ref="G", alt="A"),
    )
    pair = generate_pairs(variants)[0]
    scored = score_pair(
        pair,
        phenotype_score=0.0,
        mechanism_score=0.0,
        weights=ScoringWeights(),
        phase_weights=PhaseWeights(),
        frequency=FrequencyThresholds(),
        quality=QualityThresholds(),
        clock=FixedClock(datetime(2026, 1, 1, tzinfo=UTC)),
    )

    assert scored.scores.contradiction_penalty == 0.0
    assert not any(
        item.direction is EvidenceDirection.CONTRADICTS for item in scored.contradicting_evidence
    ), (
        "a low component score was emitted as CONTRADICTS; reserve that direction "
        "for genuine opposing evidence (GP-19)"
    )


@pytest.mark.unit
def test_genuine_contradiction_is_still_recorded_and_penalised() -> None:
    """The counterpart: real opposing evidence must survive AND carry a penalty."""
    from mva.models.evidence import EvidenceDirection

    variants = tuple(
        apply_soft_flags(
            (
                make_variant(
                    position=40_205_000,
                    ref="A",
                    alt="G",
                    allele_frequency=0.12,
                    impact=ImpactSeverity.MODERATE,
                    consequence_terms=("missense_variant",),
                ),
                make_variant(
                    position=40_206_000,
                    ref="T",
                    alt="C",
                    allele_frequency=0.09,
                    impact=ImpactSeverity.MODERATE,
                    consequence_terms=("missense_variant",),
                ),
            ),
            frequency=FrequencyThresholds(),
            quality=QualityThresholds(),
        )
    )
    pair = generate_pairs(variants)[0]
    scored = score_pair(
        pair,
        phenotype_score=0.5,
        mechanism_score=0.5,
        weights=ScoringWeights(),
        phase_weights=PhaseWeights(),
        frequency=FrequencyThresholds(),
        quality=QualityThresholds(),
        clock=FixedClock(datetime(2026, 1, 1, tzinfo=UTC)),
    )

    assert scored.scores.contradiction_penalty > 0.0
    assert any(
        item.direction is EvidenceDirection.CONTRADICTS for item in scored.contradicting_evidence
    ), "a common-variant pair recorded no contradicting evidence"


# ---------------------------------------------------------------------------
# 10. Adversarial-review regressions
# ---------------------------------------------------------------------------


def _split_multiallelic_pair() -> tuple[VariantRecord, VariantRecord]:
    """The `chr11:5000000 A>G,GT  GT=1/2  DP=44  AD=2,21,21` site, decomposed.

    Each record holds the SITE's ref depth (2) and its own alt depth (21), so
    `alt/(ref+alt)` is 0.913 for both while the site allele fraction is 0.477.
    """

    def _allele(alt: str, index: int, impact: ImpactSeverity) -> VariantRecord:
        return make_variant(
            contig="chr11",
            position=5000000,
            ref="A",
            alt=alt,
            genotype_string="1/2",
            depth=44,
            ref_reads=2,
            alt_reads=21,
            genotype_quality=97,
            normalisation_ops=(OP_SPLIT_MULTIALLELIC,),
            alt_allele_index=index,
            gene="SYNTHMUL4",
            allele_frequency=None,
            impact=impact,
        )

    return (
        _allele("G", 1, ImpactSeverity.MODERATE),
        _allele("GT", 2, ImpactSeverity.HIGH),
    )


@pytest.mark.unit
def test_split_multiallelic_call_is_not_penalised_for_its_vcf_formatting() -> None:
    """G1: the prioritisation stage bands the SITE allele fraction, not alt/(ref+alt).

    Both stages read `Genotype.allele_balance` directly, undoing the site-aware
    fix that `mva.ingestion.qc` already applies. The result: a textbook compound
    heterozygote at a multiallelic site was flagged `low_quality_call`, given an
    analytical validity of 0.75 with a fabricated mosaicism note, and lost 0.178
    of composite to the shape of its VCF line.
    """
    variants = list(_split_multiallelic_pair())
    for variant in variants:
        assert variant.genotype.allele_balance == pytest.approx(0.913, abs=0.001)
        assert variant.allele_fraction == pytest.approx(0.477, abs=0.001)

    flagged = apply_soft_flags(variants, frequency=FREQUENCY, quality=QUALITY)
    for variant in flagged:
        assert FLAG_LOW_QUALITY_CALL not in variant.qc_flags
        assert FLAG_POSSIBLE_MOSAIC not in variant.qc_flags

        score, notes = _variant_analytical_validity(variant, QUALITY)
        assert score == pytest.approx(1.0)
        assert notes == []


@pytest.mark.unit
def test_split_multiallelic_pair_keeps_its_composite_and_asks_no_false_questions() -> None:
    """The same record, scored end to end: no call-validity question is fabricated."""
    flagged = apply_soft_flags(
        list(_split_multiallelic_pair()), frequency=FREQUENCY, quality=QUALITY
    )
    pair = next(c for c in generate_pairs(flagged) if c.is_pair)
    scored = score_all((pair,))[0]

    assert scored.scores.analytical_validity == pytest.approx(1.0)
    assert scored.scores.contradiction_penalty == pytest.approx(0.0)
    assert not any("call-validity" in q.question_id for q in scored.open_questions)
    assert "mosaic" not in scored.rationale.lower()


@pytest.mark.unit
@pytest.mark.parametrize(
    ("ref_reads", "alt_reads"),
    [
        (17, 3),  # fraction 0.15 — below the het band, above the mosaic floor
        (19, 1),  # fraction 0.05 — below the mosaic floor: artifact
        (20, 20),  # fraction 0.50 — inside the band: nothing to say
        (30, 10),  # fraction 0.25 — exactly at the band edge
        (2, 38),  # fraction 0.95 — ABOVE the band: never mosaicism
        (1, 49),  # fraction 0.98 — likewise
    ],
)
def test_only_a_low_allele_fraction_is_described_as_mosaic(ref_reads: int, alt_reads: int) -> None:
    """G2: mosaicism dilutes the alternate allele; it can never concentrate it.

    The branch was `if not (min <= b <= max): if b >= mosaic_floor: MOSAIC`, and
    since the floor (0.10) is below the ceiling (0.75) every high-balance het
    qualified. An allele balance of 0.95 was reported as "within the mosaic
    window, so treated as possible mosaicism rather than noise", and that
    sentence reached the dossier.
    """
    variant = make_variant(depth=ref_reads + alt_reads, ref_reads=ref_reads, alt_reads=alt_reads)
    fraction = variant.allele_fraction
    assert fraction is not None

    _, notes = _variant_analytical_validity(variant, QUALITY)
    joined = " ".join(notes).lower()

    # The word appears only where the alternate allele is DILUTED. Above the
    # band it must not appear at all, in any form a reader could take as a
    # mosaic finding.
    assert ("mosaic" in joined) is (fraction < QUALITY.min_allele_balance_het)
    if fraction > QUALITY.max_allele_balance_het:
        assert "above the" in joined
        assert "dropout" in joined


@pytest.mark.unit
def test_high_allele_balance_is_penalised_without_being_called_mosaic() -> None:
    """The above-band case keeps its own multiplier: narrowed, not disarmed."""
    high = make_variant(depth=40, ref_reads=2, alt_reads=38)
    clean = make_variant(depth=40, ref_reads=20, alt_reads=20)

    high_score, _ = _variant_analytical_validity(high, QUALITY)
    clean_score, _ = _variant_analytical_validity(clean, QUALITY)

    assert high_score < clean_score
    assert high_score == pytest.approx(0.60)


def _founder_allele() -> VariantRecord:
    """8e-6 across 125,000 global chromosomes; one allele in a 40-chromosome group."""
    return make_variant(
        frequencies=(
            PopulationFrequency(
                source="SYNTH_gnomAD",
                version="v0.0-synthetic",
                population="global",
                allele_frequency=8e-6,
                allele_count=1,
                allele_number=125000,
                homozygote_count=0,
                filter_status="PASS",
            ),
            PopulationFrequency(
                source="SYNTH_gnomAD",
                version="v0.0-synthetic",
                population="ami",
                allele_frequency=0.025,
                allele_count=1,
                allele_number=40,
                homozygote_count=0,
                filter_status="PASS",
            ),
        ),
    )


@pytest.mark.unit
def test_a_singleton_in_a_tiny_cohort_does_not_make_a_variant_common() -> None:
    """G3: popmax had no cohort-size guard, so AC=1/AN=40 read as AF 0.025.

    A genuine founder allele was flagged `common_variant`, its rarity component
    fell to 0.1196, and `select_candidate_variants` dropped it. See ADR 0010.
    """
    variant = _founder_allele()

    observed = variant.max_allele_frequency(min_allele_number=FREQUENCY.min_allele_number)
    assert observed is not None
    assert observed.population == "global"

    flagged = apply_soft_flags([variant], frequency=FREQUENCY, quality=QUALITY)[0]
    assert FLAG_COMMON_VARIANT not in flagged.qc_flags

    selected = select_candidate_variants([flagged], frequency=FREQUENCY)
    assert len(selected) == 1, "a founder allele was removed from the plausible set"

    scored = score_all(generate_pairs([flagged]))[0]
    assert scored.scores.rarity == pytest.approx(1.0)
    assert scored.scores.contradiction_penalty == pytest.approx(0.0)


@pytest.mark.unit
def test_an_excluded_population_is_named_rather_than_dropped_silently() -> None:
    """The guard is reported: a reader can see what it set aside and disagree."""
    selection = _founder_allele().select_max_allele_frequency(
        min_allele_number=FREQUENCY.min_allele_number
    )
    assert [f.population for f in selection.excluded] == ["ami"]

    scored = score_all(generate_pairs([_founder_allele()]))[0]
    assert "ami" in scored.rationale
    assert "AN 40" in scored.rationale


@pytest.mark.unit
def test_the_guard_still_calls_a_well_powered_common_allele_common() -> None:
    """G3 narrows the popmax; it must not blind it."""
    common = make_variant(allele_frequency=0.12, allele_number=152312, homozygote_count=1102)
    flagged = apply_soft_flags([common], frequency=FREQUENCY, quality=QUALITY)[0]
    assert FLAG_COMMON_VARIANT in flagged.qc_flags


@pytest.mark.unit
def test_ingestion_quality_flags_partition_the_shared_vocabulary() -> None:
    """G4: the set is derived from the emitted vocabulary, not re-typed by hand.

    Eight of the twelve hand-written names were emitted nowhere in `src/`, and
    three that ingestion does emit were missing.
    """
    recognised = INGESTION_QUALITY_FLAGS | INGESTION_NEUTRAL_FLAGS
    assert recognised >= INGESTION_QC_FLAGS, (
        "ingestion emits QC flags the prioritisation stage does not recognise: "
        + ", ".join(sorted(INGESTION_QC_FLAGS - recognised))
    )
    assert not (INGESTION_QUALITY_FLAGS & INGESTION_NEUTRAL_FLAGS)
    # Every name that means "the call may not be real" is a quality flag.
    assert {FLAG_REF_ALLELE_MISMATCH, "filtered_by_caller", "high_allele_balance"} <= (
        INGESTION_QUALITY_FLAGS
    )
    # And the ones that deliberately are not are exactly the documented three.
    documented_neutral = {"possible_mosaic", "no_quality_metrics", "no_caller_filter"}
    assert set(INGESTION_NEUTRAL_FLAGS) == documented_neutral


@pytest.mark.unit
def test_a_reference_allele_mismatch_is_penalised_rather_than_ignored() -> None:
    """A confidently wrong coordinate must not rank like a clean call.

    VCF says REF=C, the reference says A — a build or patch mismatch. QC flags it
    and raises a `contradicts` item; the prioritisation stage then scored it
    `analytical validity 1.0` with no flag, no contradiction and no open question.
    """
    mismatched = make_variant(qc_flags=(FLAG_REF_ALLELE_MISMATCH,))
    clean = make_variant()

    flagged = apply_soft_flags([mismatched], frequency=FREQUENCY, quality=QUALITY)[0]
    assert FLAG_LOW_QUALITY_CALL in flagged.qc_flags

    penalised = score_all(generate_pairs([flagged]))[0]
    baseline = score_all(generate_pairs([clean]))[0]

    assert penalised.scores.analytical_validity < baseline.scores.analytical_validity
    assert penalised.scores.contradiction_penalty > 0.0
    assert penalised.composite < baseline.composite
    assert any("call-validity" in q.question_id for q in penalised.open_questions)


@pytest.mark.unit
def test_the_rendered_consequence_term_belongs_to_the_worst_impact_transcript() -> None:
    """G6: the rationale used `min()` over strings, which sorts alphabetically.

    A `stop_gained` + `intron_variant` variant rendered as
    "intron_variant (HIGH impact), loss-of-function term stop_gained". The score
    was right; the sentence a clinician reads was not.
    """
    variant = make_variant(
        consequence_terms=("stop_gained",),
        impact=ImpactSeverity.HIGH,
        extra_consequences=(
            ConsequenceAnnotation(
                gene_symbol=GENE,
                gene_id="SYNTHG0001",
                transcript_id="SYNTHT0009.1",
                consequence_terms=("intron_variant",),
                impact=ImpactSeverity.MODIFIER,
                source_tool="SYNTH_vep",
                source_tool_version="v0.0-synthetic",
            ),
        ),
    )
    assert variant.worst_impact_for_gene(GENE) is ImpactSeverity.HIGH

    scored = score_all(generate_pairs([variant]))[0]
    assert "stop_gained (HIGH impact)" in scored.rationale
    assert "intron_variant (HIGH impact)" not in scored.rationale


@pytest.mark.unit
def test_a_phased_multiallelic_het_still_resolves_trans() -> None:
    """G7: `1|2` lost its own ALT index at the split, so resolved phase was discarded.

    Allele 1 is unambiguously on haplotype slot 0 and the partner allele on slot
    1. Without the recorded index `_alt_haplotype_indices('1|2')` returned
    `{0, 1}` and the pair bailed to UNKNOWN — failing safe, but throwing away a
    genuine *trans* call.
    """
    split = make_variant(
        contig="chr1",
        position=100,
        ref="A",
        alt="G",
        genotype_string="1|2",
        phased=True,
        phase_set=100,
        alt_allele_index=1,
        normalisation_ops=(OP_SPLIT_MULTIALLELIC,),
        depth=44,
        ref_reads=2,
        alt_reads=21,
    )
    partner = make_variant(
        contig="chr1",
        position=300,
        ref="C",
        alt="T",
        genotype_string="0|1",
        phased=True,
        phase_set=100,
        alt_allele_index=1,
    )

    phase = infer_phase(split, partner)
    assert phase.status is PhaseStatus.TRANS_CONFIRMED
    assert phase.method == "phase_set"


@pytest.mark.unit
def test_an_unrecorded_alt_index_still_refuses_to_guess() -> None:
    """The G7 fix must not become a licence to resolve phase it cannot see."""
    a = make_variant(contig="chr1", position=100, genotype_string="1|2", phased=True, phase_set=100)
    b = make_variant(
        contig="chr1",
        position=300,
        ref="C",
        alt="T",
        genotype_string="0|1",
        phased=True,
        phase_set=100,
    )
    assert infer_phase(a, b).status is PhaseStatus.UNKNOWN


@pytest.mark.unit
def test_every_inheritance_model_is_produced_or_documented_as_unreachable() -> None:
    """G8: five enum members were unreachable while the docs implied coverage.

    A `chrM` hemizygous call scored 0.90 for "a single call accounting for both
    gene copies" — mitochondrial DNA has no gene copies. Either the pipeline
    produces a model or it says, in code, why it cannot.
    """
    assert PRODUCED_INHERITANCE_MODELS.isdisjoint(UNPRODUCED_INHERITANCE_MODELS)
    assert PRODUCED_INHERITANCE_MODELS | set(UNPRODUCED_INHERITANCE_MODELS) == set(
        InheritanceModel
    ), "an InheritanceModel member is neither produced nor on the documented exclusion list"
    for model, reason in UNPRODUCED_INHERITANCE_MODELS.items():
        assert reason.strip(), f"{model.value} is excluded with no stated reason"

    corpus = [
        # compound heterozygous: two hets in one gene
        make_variant(position=40200000, ref="C", alt="T"),
        make_variant(position=40210500, ref="G", alt="A"),
        # homozygous recessive
        make_variant(
            position=40206000,
            ref="T",
            alt="C",
            zygosity=Zygosity.HOM_ALT,
            genotype_string="1/1",
            gene="SYNTHHOM1",
        ),
        # X-linked recessive
        make_variant(
            contig="chrX",
            position=1000,
            zygosity=Zygosity.HEMIZYGOUS,
            genotype_string="1",
            gene="SYNTHX1",
        ),
        # mitochondrial: chrM, whatever the zygosity label says
        make_variant(
            contig="chrM",
            position=3243,
            ref="A",
            alt="G",
            zygosity=Zygosity.HEMIZYGOUS,
            genotype_string="1",
            gene="SYNTHMT1",
        ),
        # mosaic: flagged by the QC stage on allele fraction
        make_variant(
            position=99000,
            gene="SYNTHMOS1",
            depth=40,
            ref_reads=34,
            alt_reads=6,
            qc_flags=(FLAG_POSSIBLE_MOSAIC,),
        ),
        # unknown: a lone HIGH-impact het with nothing to pair against
        make_variant(position=77000, gene="SYNTHLONE1"),
    ]
    produced = {candidate.inheritance_model for candidate in generate_pairs(corpus)}

    assert produced == PRODUCED_INHERITANCE_MODELS, (
        "PRODUCED_INHERITANCE_MODELS does not match what the pipeline produces; "
        f"missing {sorted(m.value for m in PRODUCED_INHERITANCE_MODELS - produced)}, "
        f"unexpected {sorted(m.value for m in produced - PRODUCED_INHERITANCE_MODELS)}"
    )


@pytest.mark.unit
def test_a_mitochondrial_call_is_not_described_as_accounting_for_both_gene_copies() -> None:
    """The rationale a reader sees must not claim mtDNA has two gene copies."""
    mt = make_variant(
        contig="chrM",
        position=3243,
        ref="A",
        alt="G",
        zygosity=Zygosity.HEMIZYGOUS,
        genotype_string="1",
        gene="SYNTHMT1",
    )
    scored = score_all(generate_pairs([mt]))[0]

    assert scored.candidate.inheritance_model is InheritanceModel.MITOCHONDRIAL
    assert "both gene copies" not in scored.rationale
    assert "heteroplasmy" in scored.rationale
    assert scored.scores.inheritance_consistency == pytest.approx(0.50)


# ---------------------------------------------------------------------------
# 12. The ranking key is TOTAL (GP-30)
#
# `rank_pairs` orders on (-composite, variant_a, variant_b-or-sentinel, pair_id)
# and always has. `CandidatePair.sort_key()` — the key every *downstream* sort
# uses — carried only the first two components, so two candidates with the same
# score and the same first variant compared equal and their emitted order was
# decided by the caller's list order. The tests below lock both keys against
# that shape, and lock the two to agree: a total key that only one of them uses
# leaves the bug alive at the other site.
# ---------------------------------------------------------------------------


def _tied_scored_pair(
    variant_a: VariantRecord, variant_b: VariantRecord | None, *, composite: float
) -> ScoredPair:
    """A ``ScoredPair`` with a composite chosen by the test, not by the scorer."""
    variant_ids = (
        (variant_a.variant_id,)
        if variant_b is None
        else (variant_a.variant_id, variant_b.variant_id)
    )
    candidate = PairCandidate(
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
        flags=(),
    )
    return ScoredPair(
        candidate=candidate,
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
        composite=composite,
        supporting_evidence=(),
        contradicting_evidence=(),
        open_questions=(),
        rationale="fixture",
    )


@pytest.mark.unit
def test_rank_pairs_is_deterministic_for_tied_pairs_sharing_first_variant() -> None:
    """Equal composite, equal first variant, different second variant.

    The tie must break on the second variant's coordinate, and reversing the
    input must not change the ranked list.
    """
    shared = rare_stop_gained()
    left = _tied_scored_pair(shared, rare_splice_donor(), composite=0.80)
    right = _tied_scored_pair(shared, common_synonymous(), composite=0.80)

    forward = rank_pairs([left, right], clock=CLOCK)
    backward = rank_pairs([right, left], clock=CLOCK)

    assert [pair.pair_id for pair in forward] == [pair.pair_id for pair in backward], (
        "reversing two tied candidates reordered the ranked list"
    )
    # common_synonymous() is at 40205000, rare_splice_donor() at 40210500, so the
    # coordinate tiebreak puts the former first.
    assert forward[0].variant_ids[1] == common_synonymous().variant_id
    assert [pair.rank for pair in forward] == [1, 2]


@pytest.mark.unit
def test_rank_pairs_orders_a_single_before_a_tied_pair_sharing_its_variant() -> None:
    """``NO_SECOND_VARIANT`` sorts below every real coordinate, by construction."""
    shared = rare_stop_gained()
    single = _tied_scored_pair(shared, None, composite=0.80)
    pair = _tied_scored_pair(shared, rare_splice_donor(), composite=0.80)

    forward = rank_pairs([single, pair], clock=CLOCK)
    backward = rank_pairs([pair, single], clock=CLOCK)

    assert [p.pair_id for p in forward] == [p.pair_id for p in backward]
    assert forward[0].variant_b is None, "the single-variant candidate must rank first"
    assert forward[0].sort_key()[2] == NO_SECOND_VARIANT


@pytest.mark.unit
def test_candidate_pair_sort_key_is_total_over_a_tied_block() -> None:
    """No two ranked candidates may share a sort key — and this block would.

    ``rank_pairs`` is not the only thing that sorts these objects: reporting and
    submission shaping sort them again via ``CandidatePair.sort_key()``. If that
    key admits ties, those later sorts are only as stable as the list handed to
    them and the rendered submission becomes a function of input order. All three
    candidates here share a composite *and* a first variant, which is exactly the
    shape the two-component key could not separate.
    """
    shared = rare_stop_gained()
    ranked = rank_pairs(
        [
            _tied_scored_pair(shared, rare_splice_donor(), composite=0.80),
            _tied_scored_pair(shared, common_synonymous(), composite=0.80),
            _tied_scored_pair(shared, None, composite=0.80),
        ],
        clock=CLOCK,
    )
    keys = [pair.sort_key() for pair in ranked]

    assert len(ranked) == 3
    assert len(set(keys)) == len(keys), "CandidatePair.sort_key() is not injective"
    assert keys == sorted(keys), (
        "rank_pairs' ordering and CandidatePair.sort_key() disagree; a downstream "
        "re-sort would reorder the ranked list"
    )


@pytest.mark.unit
def test_the_ranked_list_is_already_in_candidate_sort_key_order() -> None:
    """The two keys agree on real pipeline output, so re-sorting is a no-op."""
    first, second = cis_pair_variants()
    ranked = run_pipeline(
        [rare_splice_donor(), common_synonymous(), second, rare_stop_gained(), first]
    )
    keys = [pair.sort_key() for pair in ranked]

    assert len(ranked) > 1
    assert len(set(keys)) == len(keys), "CandidatePair.sort_key() is not injective"
    assert keys == sorted(keys), "rank order and sort_key() order disagree"


@pytest.mark.unit
def test_ranking_only_breaks_ties_and_never_reorders_on_score() -> None:
    """The appended tiebreakers must not outvote a genuine score difference."""
    first, second = cis_pair_variants()
    ranked = run_pipeline(
        [rare_splice_donor(), common_synonymous(), second, rare_stop_gained(), first]
    )
    for higher, lower in itertools.pairwise(ranked):
        assert higher.composite_score >= lower.composite_score, (
            "a tiebreak component reordered two candidates that differ on composite"
        )


@pytest.mark.unit
@given(order=st.permutations(range(5)))
@settings(max_examples=60, deadline=None)
def test_ranking_is_invariant_under_any_permutation_of_the_input(order: list[int]) -> None:
    """GP-30 as a property, over the whole score-and-rank stage."""
    first, second = cis_pair_variants()
    variants = [rare_splice_donor(), common_synonymous(), second, rare_stop_gained(), first]
    baseline = run_pipeline(variants)
    permuted = run_pipeline([variants[i] for i in order])

    assert [pair.pair_id for pair in permuted] == [pair.pair_id for pair in baseline]
    assert [pair.composite_score for pair in permuted] == [
        pair.composite_score for pair in baseline
    ]
