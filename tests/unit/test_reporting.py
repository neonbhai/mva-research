"""Reporting-layer behaviour that must never regress.

Each test here corresponds to a way the reports could quietly mislead a reader:
an unsourced sentence, an inference printed as a finding, a contradiction that
disappears, a rejected compound that is never mentioned, an unknown direction
filed under agreement, symptom management presented as mechanism correction.
"""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

import pytest
from jinja2 import UndefinedError

from mva.clock import FixedClock
from mva.errors import UnsourcedAssertionError
from mva.evidence.ledger import AssertionResolver, EvidenceLedger
from mva.models import (
    ApprovalStatus,
    AssertionTier,
    CandidatePair,
    Citation,
    ComponentScores,
    ConsequenceAnnotation,
    DrugHypothesis,
    EffectDirection,
    EvidenceCategory,
    EvidenceDirection,
    EvidenceItem,
    EvidenceStrength,
    EvidenceType,
    FilterStatus,
    GenomeBuild,
    GenomicCoordinate,
    Genotype,
    ImpactSeverity,
    InheritanceModel,
    InterventionClass,
    MechanismHypothesis,
    MechanismLink,
    MechanismNode,
    MechanismNodeKind,
    OpenQuestion,
    PediatricEvidence,
    PharmacokineticProfile,
    PhaseEvidence,
    PhaseStatus,
    PopulationFrequency,
    RejectionReason,
    VariantRecord,
    Zygosity,
)
from mva.reporting import (
    NOT_MEDICAL_ADVICE,
    TIER_MARKERS,
    Assertion,
    AssertionChecker,
    build_candidate_dossier,
    build_rejection_record,
    build_track2_report,
    render_template,
)

FIXED_INSTANT = datetime(2026, 1, 1, tzinfo=UTC)
RUN_ID = "RUN-TEST-0001"
GENE = "SYNTHKIN1"

CONTRADICTION_CLAIM = (
    "A published unaffected sibling carries the same biallelic SYNTHKIN1 genotype."
)


# ---------------------------------------------------------------------------
# Fixtures, built inline so a reader can see exactly what is being asserted.
# ---------------------------------------------------------------------------


def _clock() -> FixedClock:
    return FixedClock(FIXED_INSTANT)


def _evidence(
    evidence_id: str,
    subject_id: str,
    claim: str,
    *,
    category: EvidenceCategory,
    tier: AssertionTier = AssertionTier.LITERATURE_MECHANISM,
    evidence_type: EvidenceType = EvidenceType.EXPERT_REVIEW,
    direction: EvidenceDirection = EvidenceDirection.SUPPORTS,
    strength: EvidenceStrength = EvidenceStrength.MODERATE,
) -> EvidenceItem:
    return EvidenceItem(
        evidence_id=evidence_id,
        subject_id=subject_id,
        subject_kind="pair",
        claim=claim,
        category=category,
        direction=direction,
        strength=strength,
        evidence_type=evidence_type,
        tier=tier,
        citation=Citation(source="SyntheticRef", identifier="SR-1", version="2026.01"),
        method="synthetic fixture, hand-written for this test",
        tool="test-fixture",
        tool_version="0.0.0",
        limitations="Fabricated for testing; carries no biological validity whatsoever.",
        timestamp=FIXED_INSTANT,
        run_id=RUN_ID,
    )


def _ledger() -> EvidenceLedger:
    ledger = EvidenceLedger(run_id=RUN_ID)
    ledger.extend(
        [
            _evidence(
                "EV-ANAL-0001",
                "PAIR-1",
                "Both calls pass FILTER at depth above 30 with balanced allele fractions.",
                category=EvidenceCategory.ANALYTICAL,
                tier=AssertionTier.OBSERVED_DATA,
                evidence_type=EvidenceType.DIRECT_MEASUREMENT,
            ),
            _evidence(
                "EV-POPU-0001",
                "PAIR-1",
                "Neither allele is observed in the synthetic reference cohort.",
                category=EvidenceCategory.POPULATION,
                tier=AssertionTier.INFERENCE,
                evidence_type=EvidenceType.PIPELINE_INFERENCE,
            ),
            _evidence(
                "EV-MECH-0001",
                "MECH-1",
                "Loss of the SYNTHKIN1 product abolishes checkpoint arrest in cell models.",
                category=EvidenceCategory.MECHANISM,
                evidence_type=EvidenceType.CELL_LINE,
            ),
            _evidence(
                "EV-DRUG-0001",
                "DRUG-A",
                "Synthexostat raises checkpoint signalling in a relevant cell model.",
                category=EvidenceCategory.DRUG,
                evidence_type=EvidenceType.CELL_LINE,
            ),
            _evidence(
                "EV-DRUG-0002",
                "DRUG-C",
                "Synthicam reduces seizure frequency in paediatric epilepsy.",
                category=EvidenceCategory.DRUG,
                evidence_type=EvidenceType.HUMAN_TRIAL,
            ),
            _evidence(
                "EV-DRUG-0003",
                "DRUG-D",
                "Synthexadrug has a context-dependent effect on checkpoint signalling.",
                category=EvidenceCategory.DRUG,
                tier=AssertionTier.INFERENCE,
                evidence_type=EvidenceType.PIPELINE_INFERENCE,
            ),
            _evidence(
                "EV-CONT-0001",
                "PAIR-1",
                CONTRADICTION_CLAIM,
                category=EvidenceCategory.CONTRADICTION,
                direction=EvidenceDirection.CONTRADICTS,
                evidence_type=EvidenceType.HUMAN_CASE_REPORT,
                strength=EvidenceStrength.STRONG,
            ),
        ]
    )
    return ledger


def _resolver() -> AssertionResolver:
    return AssertionResolver(_ledger())


def _variant(position: int, ref: str, alt: str, *, contig: str = "chr15") -> VariantRecord:
    return VariantRecord(
        coordinate=GenomicCoordinate(
            build=GenomeBuild.GRCH38, contig=contig, position=position, ref=ref, alt=alt
        ),
        genotype=Genotype(
            zygosity=Zygosity.HET,
            genotype_string="0/1",
            depth=42,
            ref_reads=20,
            alt_reads=22,
            genotype_quality=99,
        ),
        filter_status=FilterStatus.PASS,
        consequences=(
            ConsequenceAnnotation(
                gene_symbol=GENE,
                transcript_id="ENST00000000001",
                consequence_terms=("missense_variant",),
                impact=ImpactSeverity.MODERATE,
                hgvs_p="p.Arg100Trp",
                source_tool="synthetic-consequence",
                source_tool_version="0.0.0",
            ),
        ),
        population_frequencies=(
            PopulationFrequency(
                source="SyntheticFrequencies",
                version="v0",
                population="global",
                allele_frequency=1e-6,
            ),
        ),
        source_artifact="synthetic-fixture",
    )


def _scores() -> ComponentScores:
    return ComponentScores(
        analytical_validity=0.92,
        rarity=0.88,
        molecular_consequence=0.75,
        inheritance_consistency=0.55,
        phenotype_similarity=0.70,
        mechanistic_relevance=0.60,
        evidence_quality=0.50,
        contradiction_penalty=0.20,
    )


def _pair() -> CandidatePair:
    return CandidatePair(
        pair_id="PAIR-1",
        gene_symbol=GENE,
        variant_a=_variant(40200000, "C", "T"),
        variant_b=_variant(40210500, "G", "A"),
        inheritance_model=InheritanceModel.COMPOUND_HETEROZYGOUS,
        phase=PhaseEvidence(status=PhaseStatus.UNKNOWN, method="none", distance_bp=10500),
        scores=_scores(),
        composite_score=0.71,
        supporting_evidence_ids=("EV-ANAL-0001", "EV-POPU-0001", "EV-MECH-0001"),
        contradicting_evidence_ids=("EV-CONT-0001",),
        missing_evidence=(
            OpenQuestion(
                question_id="OQ-PHASE-01",
                question="Are the two variants in trans?",
                why_it_matters="Without trans, one allele is intact and the model does not apply.",
                resolving_test="Parental segregation or long-read phasing",
                blocking=True,
            ),
        ),
        recommended_next_test="Trio phasing of both sites.",
        discriminating_experiment="Western blot of the SYNTHKIN1 product in patient fibroblasts.",
        rank_rationale="Highest composite in the synthetic cohort.",
        flags=("phase_unknown",),
    )


def _mechanism() -> MechanismHypothesis:
    return MechanismHypothesis(
        mechanism_id="MECH-1",
        gene_symbol=GENE,
        pair_id="PAIR-1",
        summary="Biallelic loss of SYNTHKIN1 removes checkpoint arrest, permitting aneuploidy.",
        nodes=(
            MechanismNode(
                node_id="N-PROTEIN",
                kind=MechanismNodeKind.PROTEIN,
                label="SYNTHKIN1 kinase",
                state_in_patient=EffectDirection.LOSS_OF_FUNCTION,
                deviation_is_pathological=True,
            ),
            MechanismNode(
                node_id="N-CHECKPOINT",
                kind=MechanismNodeKind.CELLULAR_PROCESS,
                label="Spindle assembly checkpoint",
                state_in_patient=EffectDirection.DECREASE,
                deviation_is_pathological=True,
            ),
        ),
        links=(
            MechanismLink(
                link_id="L1",
                source_node_id="N-PROTEIN",
                target_node_id="N-CHECKPOINT",
                relation="fails_to_sustain",
                direction=EffectDirection.DECREASE,
                tier=AssertionTier.LITERATURE_MECHANISM,
                strength=EvidenceStrength.MODERATE,
                evidence_ids=("EV-MECH-0001",),
                is_directly_demonstrated=False,
                uncertainty="Demonstrated in cell lines only; no patient-derived system.",
            ),
        ),
        disease_direction=EffectDirection.DECREASE,
        therapeutic_target_node_id="N-CHECKPOINT",
        required_correction=EffectDirection.RESTORE,
        supporting_evidence_ids=("EV-MECH-0001",),
        uncertainties=("Mechanistic heterogeneity across genes under this clinical label.",),
        developmental_window_caveat=(
            "Structural damage established in utero cannot be reversed by post-natal "
            "checkpoint modulation."
        ),
    )


def _drug(
    drug_id: str,
    name: str,
    *,
    intervention_class: InterventionClass,
    observed: EffectDirection,
    evidence_ids: tuple[str, ...],
    rejected: bool = False,
    rejection_reasons: tuple[RejectionReason, ...] = (),
    rejection_rationale: str = "",
) -> DrugHypothesis:
    return DrugHypothesis(
        drug_id=drug_id,
        name=name,
        approved_name=name.lower(),
        approval_status=ApprovalStatus.APPROVED_OTHER_INDICATION,
        intervention_class=intervention_class,
        target="checkpoint signalling",
        target_node_id="N-CHECKPOINT",
        mechanism_of_action="synthetic mechanism of action for testing",
        required_direction=EffectDirection.RESTORE,
        observed_direction=observed,
        is_direct_evidence=True,
        strongest_evidence_type=EvidenceType.CELL_LINE,
        pediatric_evidence=PediatricEvidence(
            has_pediatric_exposure=True, youngest_age_studied="4 years"
        ),
        pharmacokinetics=PharmacokineticProfile(
            route="oral",
            achievable_plasma_concentration_um=2.0,
            required_effective_concentration_um=1.0,
        ),
        worsens_chromosomal_instability=False,
        evidence_ids=evidence_ids,
        proposed_validation_experiment="Micronucleus assay in patient-derived fibroblasts.",
        score=0.60,
        rejected=rejected,
        rejection_reasons=rejection_reasons,
        rejection_rationale=rejection_rationale,
    )


def _accepted_agreeing() -> DrugHypothesis:
    return _drug(
        "DRUG-A",
        "Synthexostat",
        intervention_class=InterventionClass.DISEASE_MODIFYING,
        observed=EffectDirection.INCREASE,
        evidence_ids=("EV-DRUG-0001",),
    )


def _accepted_undetermined() -> DrugHypothesis:
    return _drug(
        "DRUG-D",
        "Synthexadrug",
        intervention_class=InterventionClass.DISEASE_MODIFYING,
        observed=EffectDirection.CONTEXT_DEPENDENT,
        evidence_ids=("EV-DRUG-0003",),
    )


def _symptomatic() -> DrugHypothesis:
    return _drug(
        "DRUG-C",
        "Synthicam",
        intervention_class=InterventionClass.SYMPTOMATIC,
        observed=EffectDirection.INCREASE,
        evidence_ids=("EV-DRUG-0002",),
    )


def _rejected_wrong_direction() -> DrugHypothesis:
    """A checkpoint inhibitor: right target, inverted sign. Constructible only as rejected."""
    return _drug(
        "DRUG-B",
        "Synthinib",
        intervention_class=InterventionClass.DISEASE_MODIFYING,
        observed=EffectDirection.DECREASE,
        evidence_ids=(),
        rejected=True,
        rejection_reasons=(RejectionReason.WRONG_DIRECTION,),
        rejection_rationale=(
            "Pushes checkpoint signalling further down, the same direction as the disease."
        ),
    )


def _section(report: str, heading: str) -> str:
    """The body of one '## ' section, up to the next one."""
    marker = f"\n## {heading}"
    start = report.index(marker) + len(marker)
    remainder = report[start:]
    end = remainder.find("\n## ")
    return remainder if end == -1 else remainder[:end]


def _full_track2_report() -> str:
    return build_track2_report(
        _mechanism(),
        [_accepted_agreeing(), _accepted_undetermined(), _symptomatic()],
        [_rejected_wrong_direction()],
        pair=_pair(),
        resolver=_resolver(),
        clock=_clock(),
    )


# ---------------------------------------------------------------------------
# 1. GP-10 — an unsupported assertion never reaches the report.
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_unsourced_assertion_is_refused() -> None:
    checker = AssertionChecker(_resolver(), strict=True)
    claim = Assertion(
        text="SYNTHKIN1 loss causes the proband's phenotype.",
        tier=AssertionTier.INFERENCE,
        evidence_ids=(),
    )
    with pytest.raises(UnsourcedAssertionError) as excinfo:
        checker.check(claim)
    assert "GP-10" in str(excinfo.value)


@pytest.mark.unit
def test_assertion_citing_an_unknown_id_is_refused() -> None:
    """A dangling citation is indistinguishable from a fabricated one."""
    checker = AssertionChecker(_resolver(), strict=True)
    with pytest.raises(UnsourcedAssertionError):
        checker.check(
            Assertion(
                text="A claim resting on evidence nobody recorded.",
                tier=AssertionTier.DATABASE_ASSERTION,
                evidence_ids=("EV-DOES-NOT-EXIST",),
            )
        )


@pytest.mark.unit
def test_non_strict_checker_does_not_raise() -> None:
    checker = AssertionChecker(_resolver(), strict=False)
    claim = Assertion(text="Draft sentence.", tier=AssertionTier.SPECULATION, evidence_ids=())
    assert checker.check(claim) is claim


# ---------------------------------------------------------------------------
# 2. Unproven tiers are visually marked wherever they render.
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize(
    ("tier", "marker"),
    [(AssertionTier.INFERENCE, "[inferred]"), (AssertionTier.SPECULATION, "[SPECULATIVE]")],
)
def test_unproven_tier_renders_with_a_leading_marker(tier: AssertionTier, marker: str) -> None:
    claim = Assertion(text="The chain is complete.", tier=tier, evidence_ids=("EV-X",))
    rendered = claim.rendered()
    assert marker in rendered
    assert rendered.startswith(marker), "an unproven tier must be marked before the sentence"


@pytest.mark.unit
def test_every_tier_has_a_marker() -> None:
    assert set(TIER_MARKERS) == set(AssertionTier)


@pytest.mark.unit
def test_inference_marker_survives_into_the_dossier() -> None:
    """The population evidence in the fixture is INFERENCE tier."""
    dossier = build_candidate_dossier([_pair()], resolver=_resolver(), clock=_clock(), top_n=1)
    assert "[inferred]" in dossier


# ---------------------------------------------------------------------------
# 3. Contradictions survive into the output (GP-19).
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_contradicting_evidence_appears_in_the_dossier() -> None:
    dossier = build_candidate_dossier([_pair()], resolver=_resolver(), clock=_clock(), top_n=1)
    assert CONTRADICTION_CLAIM in dossier
    assert "EV-CONT-0001" in dossier
    assert "GP-19" in dossier


@pytest.mark.unit
def test_dossier_shows_the_component_vector_not_only_the_composite() -> None:
    dossier = build_candidate_dossier([_pair()], resolver=_resolver(), clock=_clock(), top_n=1)
    for component in _scores().as_dict():
        assert component.replace("_", " ") in dossier
    assert "phase is UNKNOWN" in dossier or "Phase is UNKNOWN" in dossier


# ---------------------------------------------------------------------------
# 4-6. Drug triage: rejections, tri-state direction, intervention class.
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_wrong_direction_drug_appears_in_the_rejection_record_with_its_reason() -> None:
    record = build_rejection_record([_rejected_wrong_direction()], clock=_clock())
    assert "Synthinib" in record
    assert "wrong_direction" in record
    assert "DISAGREES" in record
    assert "same direction as the disease" in record


@pytest.mark.unit
def test_rejected_drug_is_not_silently_dropped_from_the_full_report() -> None:
    report = _full_track2_report()
    assert "Synthinib" in report
    assert "## Rejection record" in report


@pytest.mark.unit
def test_undeterminable_direction_is_not_presented_as_agreement() -> None:
    report = _full_track2_report()
    agreeing = _section(report, "Disease-modifying candidates — direction AGREES")
    undetermined = _section(report, "Direction CANNOT BE DETERMINED — not presented as agreement")
    assert "Synthexadrug" not in agreeing
    assert "Synthexadrug" in undetermined
    assert "CANNOT BE DETERMINED" in undetermined


@pytest.mark.unit
def test_symptomatic_drug_is_not_in_the_disease_modifying_section() -> None:
    report = _full_track2_report()
    agreeing = _section(report, "Disease-modifying candidates — direction AGREES")
    symptomatic = _section(report, "Symptomatic, supportive, surveillance and preventive measures")
    assert "Synthicam" not in agreeing
    assert "Synthicam" in symptomatic
    assert "Synthexostat" in agreeing


@pytest.mark.unit
def test_every_drug_answers_the_eight_mandatory_questions() -> None:
    report = _full_track2_report()
    for question in (
        "What node is modified?",
        "In what direction must it move?",
        "Does the drug act that way?",
        "What evidence tier supports it",
        "Was the concentration clinically achievable?",
        "Is there paediatric exposure evidence?",
        "Could it worsen chromosome instability or cancer susceptibility?",
        "What experiment should be performed before any clinical consideration?",
    ):
        assert question in report


# ---------------------------------------------------------------------------
# 7. The banner and the caveats a Track 2 report may never lose.
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_track2_report_carries_the_mandatory_header_and_caveats() -> None:
    report = _full_track2_report()
    assert NOT_MEDICAL_ADVICE in report
    assert "NOT MEDICAL ADVICE" in report.split("## Read this first")[0]
    assert "Structural damage established in utero" in report
    assert "uncalibrated heuristics" in report
    assert "inferred rather than directly demonstrated" in report


# ---------------------------------------------------------------------------
# 8-9. Determinism and loud failure.
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_rendering_twice_is_byte_identical() -> None:
    assert _full_track2_report() == _full_track2_report()
    first = build_candidate_dossier([_pair()], resolver=_resolver(), clock=_clock(), top_n=1)
    second = build_candidate_dossier([_pair()], resolver=_resolver(), clock=_clock(), top_n=1)
    assert first == second


@pytest.mark.unit
def test_missing_template_variable_is_a_loud_error(tmp_path: Path) -> None:
    """StrictUndefined: a typo must not render as an empty caveat."""
    (tmp_path / "broken.md.j2").write_text("Caveat: {{ missing_caveat }}\n", encoding="utf-8")
    with pytest.raises(UndefinedError):
        render_template("broken.md.j2", {"present": "value"}, templates_dir=tmp_path)


@pytest.mark.unit
def test_missing_template_is_a_loud_error(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        render_template("absent.md.j2", {}, templates_dir=tmp_path)
