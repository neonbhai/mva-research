"""Unit tests for the intervention stage.

These tests exist to hold one line: **a compound that acts on the right target in
the wrong direction must never be presented as a candidate.** Everything else here
protects the distinctions that keep that line meaningful — wrong direction vs
undeterminable direction, wrong direction vs unapproved, mechanism correction vs
symptom control — because each collapse produces a confident recommendation that
is wrong in a different way.

The expectations are locked against `tests/golden/expected_drug_outcomes.tsv`.
Changing them requires a decision record (GP-32).
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from mva.clock import demo_clock
from mva.determinism import stable_hash
from mva.errors import IngestionError
from mva.interventions import (
    DIRECTION_UNKNOWN_CREDIT,
    EVIDENCE_TYPE_WEIGHT,
    INDIRECT_EVIDENCE_MULTIPLIER,
    INTERVENTION_CLASS_WEIGHT,
    REASON_PRIORITY,
    CatalogEntry,
    DrugCatalog,
    DrugTriageResult,
    assess_evidence_quality,
    assess_safety,
    check_direction,
    check_target_in_mechanism,
    generate_drug_hypotheses,
    is_chromosomal_instability_context,
    is_neurological_context,
    required_direction_for_node,
)
from mva.mechanisms import MechanismLibrary, build_mechanism
from mva.models.drug import (
    ApprovalStatus,
    DrugHypothesis,
    InterventionClass,
    RejectionReason,
)
from mva.models.evidence import EvidenceDirection, EvidenceType
from mva.models.mechanism import (
    EffectDirection,
    MechanismHypothesis,
    MechanismNode,
    MechanismNodeKind,
    directions_agree,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
KNOWLEDGE = REPO_ROOT / "knowledge" / "public"
CATALOG_PATH = KNOWLEDGE / "drug_catalog.tsv"
VERSION = "synthetic-demo-2026.1"


def _mechanism() -> MechanismHypothesis:
    library = MechanismLibrary.from_tsv(
        KNOWLEDGE / "mechanisms.tsv", KNOWLEDGE / "mechanism_meta.tsv", version=VERSION
    )
    built = build_mechanism(
        "SYNTHKIN1", pair_id="PAIR-DEMO-01", library=library, clock=demo_clock()
    )
    assert built.hypothesis is not None
    return built.hypothesis


def _catalog() -> DrugCatalog:
    return DrugCatalog.from_tsv(CATALOG_PATH, version=VERSION)


def _triage() -> DrugTriageResult:
    return generate_drug_hypotheses(mechanism=_mechanism(), catalog=_catalog(), clock=demo_clock())


@pytest.fixture
def mechanism() -> MechanismHypothesis:
    return _mechanism()


@pytest.fixture
def catalog() -> DrugCatalog:
    return _catalog()


@pytest.fixture
def result() -> DrugTriageResult:
    return _triage()


def _by_id(result: DrugTriageResult, drug_id: str) -> DrugHypothesis:
    for hypothesis in result.all_hypotheses:
        if hypothesis.drug_id == drug_id:
            return hypothesis
    raise AssertionError(f"{drug_id} was dropped from the triage output entirely")


def _entry(catalog: DrugCatalog, drug_id: str) -> CatalogEntry:
    for entry in catalog.entries():
        if entry.drug_id == drug_id:
            return entry
    raise AssertionError(f"{drug_id} missing from the catalogue")


def _record(result: DrugTriageResult, drug_id: str) -> dict[str, object]:
    for row in result.rejection_record:
        if row["drug_id"] == drug_id:
            return row
    raise AssertionError(f"no audit row for {drug_id}")


# 7 — the headline test ------------------------------------------------------


@pytest.mark.unit
def test_wrong_direction_drug_is_rejected(result: DrugTriageResult) -> None:
    """SYNTH-DRUG-B binds the exact target named in the mechanism and pushes it the
    disease way. It is the compound a target-proximity search ranks first, and it
    must be rejected for WRONG_DIRECTION."""
    synthinib = _by_id(result, "SYNTH-DRUG-B")

    assert synthinib.rejected is True
    assert synthinib.rejection_reasons[0] is RejectionReason.WRONG_DIRECTION
    assert synthinib.directions_agree is False
    assert synthinib.required_direction is EffectDirection.RESTORE
    assert synthinib.observed_direction is EffectDirection.DECREASE
    assert synthinib.rank is None
    assert synthinib.drug_id not in {d.drug_id for d in result.accepted}
    assert "disease direction" in synthinib.rejection_rationale
    # The type system, not merely the pipeline, refuses to accept it: an accepted
    # wrong-direction hypothesis cannot even be constructed.
    with pytest.raises(ValueError, match="wrong-direction agent must be rejected"):
        DrugHypothesis(
            drug_id=synthinib.drug_id,
            name=synthinib.name,
            approval_status=synthinib.approval_status,
            intervention_class=synthinib.intervention_class,
            target=synthinib.target,
            target_node_id=synthinib.target_node_id,
            mechanism_of_action=synthinib.mechanism_of_action,
            required_direction=EffectDirection.RESTORE,
            observed_direction=EffectDirection.DECREASE,
            is_direct_evidence=synthinib.is_direct_evidence,
            strongest_evidence_type=synthinib.strongest_evidence_type,
            pediatric_evidence=synthinib.pediatric_evidence,
            pharmacokinetics=synthinib.pharmacokinetics,
            proposed_validation_experiment=synthinib.proposed_validation_experiment,
            rejected=False,
        )


@pytest.mark.unit
def test_wrong_direction_is_disqualifying_even_with_good_evidence(
    mechanism: MechanismHypothesis, catalog: DrugCatalog
) -> None:
    """Upgrading the evidence tier of a wrong-direction agent must not rescue it.

    A scoring layer alone could be re-weighted into recommending it; the rejection
    is structural instead.
    """
    upgraded = replace(
        _entry(catalog, "SYNTH-DRUG-B"),
        strongest_evidence_type=EvidenceType.HUMAN_TRIAL,
        has_pediatric_exposure=True,
        worsens_cin=False,
    )
    triaged = generate_drug_hypotheses(
        mechanism=mechanism,
        catalog=DrugCatalog([upgraded], version=VERSION),
        clock=demo_clock(),
    )
    assert triaged.accepted == ()
    assert triaged.rejected[0].rejection_reasons[0] is RejectionReason.WRONG_DIRECTION


# 8 -------------------------------------------------------------------------


@pytest.mark.unit
def test_tool_compound_is_rejected_for_approval_not_for_direction(
    result: DrugTriageResult,
) -> None:
    """SYNTH-DRUG-C points the right way. The reason it fails is regulatory, and the
    report must say so — the two findings have completely different follow-ups."""
    synthophore = _by_id(result, "SYNTH-DRUG-C")

    assert synthophore.rejected is True
    assert synthophore.approval_status is ApprovalStatus.TOOL_COMPOUND
    assert synthophore.is_repurposable is False
    assert synthophore.approved_name is None, "a tool compound has no approved name"
    assert synthophore.rejection_reasons[0] is RejectionReason.NOT_APPROVED
    assert RejectionReason.WRONG_DIRECTION not in synthophore.rejection_reasons
    assert synthophore.directions_agree is True
    assert _record(result, "SYNTH-DRUG-C")["primary_reason"] == "not_approved"


# 9 -------------------------------------------------------------------------


@pytest.mark.unit
def test_context_dependent_direction_is_undetermined_penalised_and_flagged(
    result: DrugTriageResult,
) -> None:
    """SYNTH-DRUG-E is bidirectional in the literature. That is neither agreement nor
    disagreement: it is accepted, penalised, and flagged DIRECTION_UNKNOWN."""
    synthaxel = _by_id(result, "SYNTH-DRUG-G")
    top = result.accepted[0]

    assert synthaxel.rejected is False, (
        "SYNTH-DRUG-G carries the ASSUMPTION-DRUG-02 lesson since ADR 0011: an "
        "undetermined direction is neither agreement nor a wrong-direction "
        "rejection. Its CIN risk IS assessed, so nothing else disqualifies it."
    )
    assert synthaxel.directions_agree is None
    assert synthaxel.observed_direction is EffectDirection.CONTEXT_DEPENDENT
    assert synthaxel.rejection_reasons == (), "an accepted drug carries no rejection reasons"

    row = _record(result, "SYNTH-DRUG-G")
    assert row["rejected"] is False
    reasons = str(row["reasons"])
    assert "direction_unknown" in reasons
    assert "wrong_direction" not in reasons

    # Penalised: the undetermined-direction credit is strictly between 0 and 1, and
    # the candidate ranks below one whose direction agrees.
    assert 0.0 < DIRECTION_UNKNOWN_CREDIT < 1.0
    assert synthaxel.score < top.score
    assert synthaxel.rank is not None and top.rank is not None and synthaxel.rank > top.rank

    # Flagged in words, and never as support.
    undetermined = [
        item
        for item in result.evidence
        if item.subject_id == "SYNTH-DRUG-G" and "UNDETERMINED" in item.claim
    ]
    assert len(undetermined) == 1
    assert undetermined[0].direction is EvidenceDirection.NEUTRAL


@pytest.mark.unit
def test_unsigned_direction_never_counts_as_agreement() -> None:
    """Every unsigned direction resolves to None, on both sides of the comparison.

    `NO_CHANGE` is deliberately NOT in this list. A demonstrated null is a signed,
    established finding, not an absence of one; see the measured-null test below.
    """
    for unsigned in (EffectDirection.UNKNOWN, EffectDirection.CONTEXT_DEPENDENT):
        assert check_direction(required=EffectDirection.RESTORE, observed=unsigned).agrees is None
        assert check_direction(required=unsigned, observed=EffectDirection.RESTORE).agrees is None
        verdict = check_direction(required=EffectDirection.RESTORE, observed=unsigned)
        assert verdict.rejection_reason is RejectionReason.DIRECTION_UNKNOWN


@pytest.mark.unit
def test_a_measured_null_disagrees_and_scores_below_an_unsigned_peer(
    mechanism: MechanismHypothesis, catalog: DrugCatalog
) -> None:
    """A demonstrated `no_change` is evidence AGAINST, not a gap in the record.

    Treating it as unsigned would hand it the undetermined-direction credit and a
    recommendation to go and measure the sign that has already been measured. It
    must therefore score strictly below a genuinely `context_dependent` peer, whose
    sign really is open.
    """
    verdict = check_direction(required=EffectDirection.RESTORE, observed=EffectDirection.NO_CHANGE)
    assert verdict.agrees is False
    assert "demonstrated null" in verdict.rationale
    assert "disease direction" not in verdict.rationale, (
        "a null result does not push the target anywhere; the rationale must not say it does"
    )

    base = _entry(catalog, "SYNTH-DRUG-A")
    null_result = replace(
        base, drug_id="SYNTH-DRUG-NULL", observed_direction=EffectDirection.NO_CHANGE
    )
    unsigned_peer = replace(
        base, drug_id="SYNTH-DRUG-CTX", observed_direction=EffectDirection.CONTEXT_DEPENDENT
    )
    triaged = generate_drug_hypotheses(
        mechanism=mechanism,
        catalog=DrugCatalog([null_result, unsigned_peer], version=VERSION),
        clock=demo_clock(),
    )
    scored = {d.drug_id: d for d in triaged.all_hypotheses}
    assert scored["SYNTH-DRUG-NULL"].directions_agree is False
    assert scored["SYNTH-DRUG-NULL"].rejected is True
    assert scored["SYNTH-DRUG-CTX"].directions_agree is None
    assert scored["SYNTH-DRUG-NULL"].score < scored["SYNTH-DRUG-CTX"].score


# 10 ------------------------------------------------------------------------


@pytest.mark.unit
def test_symptomatic_agent_is_never_presented_as_disease_modifying(
    result: DrugTriageResult, mechanism: MechanismHypothesis
) -> None:
    """SYNTH-DRUG-D controls seizures. It acts on the organismal phenotype, not on the
    checkpoint, and the output must make that impossible to misread."""
    synthazepam = _by_id(result, "SYNTH-DRUG-D")

    assert synthazepam.rejected is False
    assert synthazepam.intervention_class is InterventionClass.SYMPTOMATIC
    assert synthazepam.intervention_class is not InterventionClass.DISEASE_MODIFYING
    assert synthazepam.target_node_id != mechanism.therapeutic_target_node_id
    assert synthazepam.target_node_id == "N7_organism"

    # It is ranked below every disease-modifying candidate that was accepted.
    disease_modifying = [
        d for d in result.accepted if d.intervention_class is InterventionClass.DISEASE_MODIFYING
    ]
    assert disease_modifying
    assert synthazepam.score < max(d.score for d in disease_modifying)
    assert (
        INTERVENTION_CLASS_WEIGHT[InterventionClass.SYMPTOMATIC]
        < INTERVENTION_CLASS_WEIGHT[InterventionClass.DISEASE_MODIFYING]
    )

    # And it says so in an evidence row, in words a renderer will surface.
    statements = [
        item
        for item in result.evidence
        if item.subject_id == "SYNTH-DRUG-D" and "does NOT correct the mechanism" in item.claim
    ]
    assert len(statements) == 1
    assert "never be presented as disease-modifying" in statements[0].claim
    assert "mechanism_mismatch" in str(_record(result, "SYNTH-DRUG-D")["reasons"])


# 11 ------------------------------------------------------------------------


@pytest.mark.unit
def test_target_not_in_mechanism_is_rejected(
    result: DrugTriageResult, mechanism: MechanismHypothesis
) -> None:
    """SYNTH-DRUG-F targets N_absent. Off the chain, nothing can be said about it."""
    synthramide = _by_id(result, "SYNTH-DRUG-F")

    assert check_target_in_mechanism("N_absent", mechanism) is False
    assert check_target_in_mechanism("N5_process", mechanism) is True
    assert synthramide.rejected is True
    assert synthramide.rejection_reasons[0] is RejectionReason.TARGET_NOT_IN_MECHANISM
    assert synthramide.approval_status is ApprovalStatus.APPROVED, (
        "an approved agent is still rejected when its target is off-chain"
    )
    assert synthramide.required_direction is EffectDirection.UNKNOWN
    assert synthramide.directions_agree is None
    assert _record(result, "SYNTH-DRUG-F")["primary_reason"] == "target_not_in_mechanism"


@pytest.mark.unit
def test_required_direction_is_per_node(mechanism: MechanismHypothesis) -> None:
    """The therapeutic target uses the curated correction; other chain nodes use the
    inverse of their state in the patient."""
    assert required_direction_for_node(mechanism, "N5_process") is EffectDirection.RESTORE
    assert required_direction_for_node(mechanism, "N6_cellphen") is EffectDirection.DECREASE
    assert required_direction_for_node(mechanism, "N7_organism") is EffectDirection.DECREASE
    assert required_direction_for_node(mechanism, "N_absent") is EffectDirection.UNKNOWN


# 12 ------------------------------------------------------------------------


@pytest.mark.unit
def test_directions_agree_is_tri_state() -> None:
    """The contract the whole stage rests on (GP-16)."""
    assert directions_agree(EffectDirection.RESTORE, EffectDirection.STABILISE) is True
    assert directions_agree(EffectDirection.RESTORE, EffectDirection.DECREASE) is False
    assert directions_agree(EffectDirection.RESTORE, EffectDirection.UNKNOWN) is None
    # None is not falsy-equivalent to False anywhere it matters.
    assert directions_agree(EffectDirection.RESTORE, EffectDirection.UNKNOWN) is not False
    assert (
        check_direction(
            required=EffectDirection.RESTORE, observed=EffectDirection.DECREASE
        ).rejection_reason
        is RejectionReason.WRONG_DIRECTION
    )


# 13 ------------------------------------------------------------------------


@pytest.mark.unit
def test_rejected_drugs_are_preserved_with_their_reasons(
    result: DrugTriageResult, catalog: DrugCatalog
) -> None:
    """GP-19: nothing is dropped. The rejections are the most reviewable output."""
    assert {d.drug_id for d in result.rejected} == {
        "SYNTH-DRUG-B",
        "SYNTH-DRUG-C",
        "SYNTH-DRUG-E",  # unassessed CIN risk, disqualifying since ADR 0011
        "SYNTH-DRUG-F",
    }
    assert len(result.accepted) + len(result.rejected) == len(catalog)
    for rejected in result.rejected:
        assert rejected.rejection_reasons, "a silent rejection destroys the audit trail"
        assert rejected.rejection_rationale.strip()
        assert rejected.evidence_ids, "GP-10"
        assert _record(result, rejected.drug_id)["rejected"] is True
    # Every rejection is also a persisted contradiction row.
    contradictions = {item.subject_id for item in result.evidence if item.is_contradiction}
    assert {"SYNTH-DRUG-B", "SYNTH-DRUG-C", "SYNTH-DRUG-F"} <= contradictions


# 14 ------------------------------------------------------------------------


@pytest.mark.unit
def test_accepted_list_is_ranked_with_synthexostat_first(result: DrugTriageResult) -> None:
    """The golden outcome: the approved, right-direction, in-vivo-evidenced,
    paediatrically-exposed agent leads the table."""
    assert [d.rank for d in result.accepted] == list(range(1, len(result.accepted) + 1))
    top = result.accepted[0]
    assert top.drug_id == "SYNTH-DRUG-A"
    assert top.rank == 1
    assert top.directions_agree is True
    assert top.approval_status is ApprovalStatus.APPROVED_OTHER_INDICATION
    assert top.has_in_vivo_evidence is True
    assert top.pediatric_evidence.has_pediatric_exposure is True
    assert top.pharmacokinetics.concentration_achievable is True
    assert top.worsens_chromosomal_instability is False
    scores = [d.score for d in result.accepted]
    assert scores == sorted(scores, reverse=True)


@pytest.mark.unit
def test_golden_outcomes_match_the_locked_expectations(result: DrugTriageResult) -> None:
    """Lock against tests/golden/expected_drug_outcomes.tsv (GP-32).

    The expectations are READ FROM THE FILE, not restated as a literal here. A
    reproducibility review found the previous version kept a hardcoded dict whose
    docstring merely *claimed* to be locked against the TSV — so flipping two
    verdicts in the golden file left the whole suite green, and the artifact
    GP-32 protects was decorative.
    """
    expected = _read_golden_outcomes()
    assert expected, "no golden drug outcomes were loaded from the TSV"

    produced = {d.drug_id for d in result.all_hypotheses}
    assert produced == set(expected), (
        "the catalogue and the golden expectations have drifted apart: "
        f"only produced={sorted(produced - set(expected))} "
        f"only expected={sorted(set(expected) - produced)}"
    )

    for drug_id, (outcome, reason) in sorted(expected.items()):
        hypothesis = _by_id(result, drug_id)
        should_reject = outcome == "rejected"
        assert hypothesis.rejected is should_reject, (
            f"{drug_id} ({hypothesis.name}): expected {outcome!r}, "
            f"got rejected={hypothesis.rejected} "
            f"reasons={[r.value for r in hypothesis.rejection_reasons]}"
        )
        if reason:
            assert hypothesis.rejection_reasons, f"{drug_id} rejected with no reason"
            assert hypothesis.rejection_reasons[0].value == reason, (
                f"{drug_id}: leading rejection reason is "
                f"{hypothesis.rejection_reasons[0].value!r}, expected {reason!r} "
                f"(all: {[r.value for r in hypothesis.rejection_reasons]})"
            )
        if outcome == "accepted_direction_undetermined":
            assert hypothesis.directions_agree is None, (
                f"{drug_id} should carry an UNDETERMINED direction; got "
                f"{hypothesis.directions_agree!r}. 'Cannot determine' and 'agrees' "
                "are different findings (ASSUMPTION-DRUG-02)."
            )
        if outcome == "accepted_symptomatic_only":
            assert hypothesis.intervention_class is InterventionClass.SYMPTOMATIC, (
                f"{drug_id} must stay classified symptomatic, not presented as "
                "correcting the mechanism (ASSUMPTION-DRUG-04)"
            )


def _read_golden_outcomes() -> dict[str, tuple[str, str]]:
    """Parse tests/golden/expected_drug_outcomes.tsv -> {drug_id: (outcome, reason)}."""
    import csv

    path = Path(__file__).resolve().parents[1] / "golden" / "expected_drug_outcomes.tsv"
    rows = [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    return {
        row["drug_id"]: (row["expected_outcome"], row["expected_primary_reason"].strip())
        for row in csv.DictReader(rows, delimiter="\t")
    }


# 15 ------------------------------------------------------------------------


@pytest.mark.unit
def test_triage_is_deterministic() -> None:
    """GP-30: two runs are byte-identical, evidence IDs included."""
    first, second = _triage(), _triage()
    assert [d.drug_id for d in first.accepted] == [d.drug_id for d in second.accepted]
    assert [d.drug_id for d in first.rejected] == [d.drug_id for d in second.rejected]
    assert [e.evidence_id for e in first.evidence] == [e.evidence_id for e in second.evidence]
    assert stable_hash([d.model_dump(mode="json") for d in first.all_hypotheses]) == stable_hash(
        [d.model_dump(mode="json") for d in second.all_hypotheses]
    )
    assert stable_hash(first.rejection_record) == stable_hash(second.rejection_record)


# catalogue ------------------------------------------------------------------


@pytest.mark.unit
def test_catalogue_preserves_unknowns_as_unknown(catalog: DrugCatalog) -> None:
    """GP-14: a blank cell is None, never a convenient default."""
    synthaxel = _entry(catalog, "SYNTH-DRUG-E")
    assert synthaxel.worsens_cin is None, "unassessed CIN risk must not become False"
    synthophore = _entry(catalog, "SYNTH-DRUG-C")
    assert synthophore.achievable_plasma_um is None
    assert synthophore.concentration_achievable is None
    assert synthophore.has_administration_route is False
    assert _entry(catalog, "SYNTH-DRUG-B").approved_name is None
    assert len(catalog) == 7  # A-G; SYNTH-DRUG-G added by ADR 0011
    assert [e.drug_id for e in catalog.for_target_node("N5_process")] == [
        "SYNTH-DRUG-A",
        "SYNTH-DRUG-B",
        "SYNTH-DRUG-C",
    ]
    assert catalog.for_target_node("N_nothing") == ()


@pytest.mark.unit
def test_catalogue_rejects_an_unrecognised_enum(tmp_path: Path) -> None:
    path = tmp_path / "catalog.tsv"
    path.write_text(
        CATALOG_PATH.read_text(encoding="utf-8").replace("\ttool_compound\t", "\tprobably_fine\t"),
        encoding="utf-8",
    )
    with pytest.raises(IngestionError, match="approval_status"):
        DrugCatalog.from_tsv(path, version=VERSION)


# safety and evidence quality -------------------------------------------------


@pytest.mark.unit
def test_cin_worsening_agent_is_disqualified_in_a_cin_context(
    mechanism: MechanismHypothesis, catalog: DrugCatalog
) -> None:
    """The mechanism is a chromosomal-instability disorder; an aneuploidy-inducing
    agent is disqualified, not merely down-ranked."""
    assert is_chromosomal_instability_context(mechanism) is True
    verdict = assess_safety(
        _entry(catalog, "SYNTH-DRUG-B"), mechanism=mechanism, clock=demo_clock()
    )
    assert verdict.disqualifying is True
    assert RejectionReason.ONCOGENIC_RISK in verdict.reasons
    critical = [c for c in verdict.concerns if c.is_disqualifying]
    assert critical and critical[0].severity == "critical"


@pytest.mark.unit
def test_missing_paediatric_exposure_is_recorded_but_not_disqualifying(
    mechanism: MechanismHypothesis, catalog: DrugCatalog
) -> None:
    """It must appear, every time, and it must not silently kill the candidate."""
    entry = replace(_entry(catalog, "SYNTH-DRUG-A"), has_pediatric_exposure=False)
    verdict = assess_safety(entry, mechanism=mechanism, clock=demo_clock())
    assert RejectionReason.NO_PEDIATRIC_EVIDENCE in verdict.reasons
    assert verdict.disqualifying is False
    assert any(c.concern_id.endswith("-PAED") for c in verdict.concerns)


@pytest.mark.unit
def test_unreachable_concentration_is_a_concern(
    mechanism: MechanismHypothesis, catalog: DrugCatalog
) -> None:
    """Mandatory question 5: a mechanism that cannot be reached in vivo is not a therapy."""
    entry = replace(
        _entry(catalog, "SYNTH-DRUG-A"), achievable_plasma_um=0.1, required_effective_um=10.0
    )
    verdict = assess_safety(entry, mechanism=mechanism, clock=demo_clock())
    assert RejectionReason.CONCENTRATION_NOT_ACHIEVABLE in verdict.reasons
    assert verdict.disqualifying is True


@pytest.mark.unit
def test_unassessed_cin_risk_blocks_in_a_cin_context(
    mechanism: MechanismHypothesis, catalog: DrugCatalog
) -> None:
    """Mandatory question 7, unanswered, is a BLOCKING gap here — so it must block.

    The model field, the evidence schema, the rendered report and
    ASSUMPTION-DRUG-07 all call an unassessed oncogenic risk a blocking gap in a
    chromosomal-instability disorder. Recording it as a non-fatal concern let such
    an agent be presented as a ranked candidate in a report that told the reader, in
    the same paragraph, that the question "must be answered before any further
    consideration".
    """
    assert is_chromosomal_instability_context(mechanism) is True
    assert is_neurological_context(mechanism) is True
    entry = replace(_entry(catalog, "SYNTH-DRUG-A"), worsens_cin=None)
    verdict = assess_safety(entry, mechanism=mechanism, clock=demo_clock())

    assert any("not been assessed" in c.description for c in verdict.concerns)
    assert verdict.disqualifying is True
    assert RejectionReason.ONCOGENIC_RISK in verdict.reasons
    blocking = [c for c in verdict.concerns if c.is_disqualifying]
    assert blocking and blocking[0].concern_id.endswith("-CIN-UNASSESSED")

    triaged = generate_drug_hypotheses(
        mechanism=mechanism,
        catalog=DrugCatalog([entry], version=VERSION),
        clock=demo_clock(),
    )
    assert triaged.accepted == ()
    assert triaged.rejected[0].rejection_reasons[0] is RejectionReason.ONCOGENIC_RISK


@pytest.mark.unit
def test_unassessed_cin_risk_outside_a_cin_context_is_recorded_not_fatal(
    catalog: DrugCatalog,
) -> None:
    """The gate is context-sensitive, not a blanket rule: elsewhere it stays a concern."""
    benign = MechanismHypothesis(
        mechanism_id="MECH-BENIGN-01",
        gene_symbol="SYNTHOTHER1",
        summary="A synthetic transporter deficiency with no chromosome-segregation component.",
        nodes=(
            MechanismNode(
                node_id="N1_protein",
                kind=MechanismNodeKind.PROTEIN,
                label="Synthetic solute transporter",
                state_in_patient=EffectDirection.LOSS_OF_FUNCTION,
                deviation_is_pathological=True,
            ),
        ),
        links=(),
        disease_direction=EffectDirection.LOSS_OF_FUNCTION,
        therapeutic_target_node_id="N1_protein",
        required_correction=EffectDirection.RESTORE,
    )
    assert is_chromosomal_instability_context(benign) is False
    verdict = assess_safety(
        replace(_entry(catalog, "SYNTH-DRUG-A"), worsens_cin=None),
        mechanism=benign,
        clock=demo_clock(),
    )
    assert verdict.disqualifying is False
    assert RejectionReason.SAFETY_CONCERN in verdict.reasons


@pytest.mark.unit
def test_a_clean_agent_raises_no_concern_but_claims_no_clearance(
    mechanism: MechanismHypothesis, catalog: DrugCatalog
) -> None:
    verdict = assess_safety(
        _entry(catalog, "SYNTH-DRUG-A"), mechanism=mechanism, clock=demo_clock()
    )
    assert verdict.concerns == ()
    assert verdict.disqualifying is False
    assert "not a safety clearance" in verdict.rationale


@pytest.mark.unit
def test_evidence_tiers_are_ranked_honestly() -> None:
    """human trial > human observational > animal > primary cells > cell line >
    biochemical binding > in silico, and indirect is discounted against direct."""
    ladder = [
        EvidenceType.HUMAN_TRIAL,
        EvidenceType.HUMAN_OBSERVATIONAL,
        EvidenceType.HUMAN_CASE_REPORT,
        EvidenceType.ANIMAL_MODEL,
        EvidenceType.PRIMARY_PATIENT_CELLS,
        EvidenceType.CELL_LINE,
        EvidenceType.BIOCHEMICAL_BINDING,
        EvidenceType.IN_SILICO_PREDICTION,
    ]
    weights = [EVIDENCE_TYPE_WEIGHT[tier] for tier in ladder]
    assert weights == sorted(weights, reverse=True)
    assert len(set(weights)) == len(weights)


@pytest.mark.unit
def test_indirect_evidence_is_penalised_against_direct(catalog: DrugCatalog) -> None:
    entry = _entry(catalog, "SYNTH-DRUG-A")
    direct_score, direct_rationale = assess_evidence_quality(entry)
    indirect_score, indirect_rationale = assess_evidence_quality(
        replace(entry, is_direct_evidence=False)
    )
    assert indirect_score < direct_score
    assert indirect_score == pytest.approx(direct_score * INDIRECT_EVIDENCE_MULTIPLIER)
    assert "INDIRECT" in indirect_rationale
    assert "direct evidence" in direct_rationale


# cross-cutting ---------------------------------------------------------------


@pytest.mark.unit
def test_every_hypothesis_answers_the_eight_mandatory_questions(
    result: DrugTriageResult,
) -> None:
    """An incomplete hypothesis must be unconstructible, rejected ones included."""
    for hypothesis in result.all_hypotheses:
        assert hypothesis.target_node_id  # 1
        assert hypothesis.required_direction  # 2
        assert hypothesis.observed_direction  # 3
        assert hypothesis.strongest_evidence_type  # 4
        assert hypothesis.pharmacokinetics is not None  # 5
        assert hypothesis.pediatric_evidence is not None  # 6
        # 7 is tri-state: None means unassessed, which is itself an answer on record.
        assert hypothesis.worsens_chromosomal_instability in (True, False, None)
        assert len(hypothesis.proposed_validation_experiment) >= 10  # 8
        assert hypothesis.pediatric_evidence.caveat


@pytest.mark.unit
def test_every_evidence_item_states_its_limitations(result: DrugTriageResult) -> None:
    """GP-17 and GP-20: no unlimited claim, and the synthetic source is labelled."""
    assert result.evidence
    for item in result.evidence:
        assert len(item.limitations.strip()) > 3
        assert "SYNTHETIC" in item.limitations
        assert item.tool.startswith("mva.interventions")
        assert item.subject_kind == "drug"


@pytest.mark.unit
def test_safety_concerns_are_linked_to_their_evidence_rows(result: DrugTriageResult) -> None:
    """GP-10: a concern printed in a report resolves to a row in the evidence store."""
    known = {item.evidence_id for item in result.evidence}
    for hypothesis in result.all_hypotheses:
        for concern in hypothesis.safety_concerns:
            assert concern.evidence_ids
            assert set(concern.evidence_ids) <= known


@pytest.mark.unit
def test_reason_priority_covers_every_rejection_reason() -> None:
    """A new reason must be given a place in the ordering, or a report could lead with
    the wrong one."""
    assert set(REASON_PRIORITY) == set(RejectionReason)
    assert len(REASON_PRIORITY) == len(set(REASON_PRIORITY))
    assert REASON_PRIORITY[0] is RejectionReason.TARGET_NOT_IN_MECHANISM
    assert REASON_PRIORITY[1] is RejectionReason.WRONG_DIRECTION


@pytest.mark.unit
def test_audit_rows_are_renderable_scalars(result: DrugTriageResult) -> None:
    """The ledger must survive a CSV/Markdown renderer without special-casing."""
    assert result.rejection_record
    for row in result.rejection_record:
        assert set(row) >= {"drug_id", "rejected", "primary_reason", "reasons", "rationale"}
        for value in row.values():
            assert isinstance(value, str | bool | int | float | type(None))


# compensatory nodes ----------------------------------------------------------


def _with_compensatory_node(node_id: str) -> MechanismHypothesis:
    """The shipped chain, with one node re-marked as a compensatory response.

    Rebuilt through the constructor rather than `model_copy` so the hypothesis
    validator actually runs on the result.
    """
    base = _mechanism()
    nodes = tuple(
        node.model_copy(update={"deviation_is_pathological": False})
        if node.node_id == node_id
        else node
        for node in base.nodes
    )
    fields = {name: getattr(base, name) for name in MechanismHypothesis.model_fields}
    fields["nodes"] = nodes
    return MechanismHypothesis(**fields)


@pytest.mark.unit
def test_a_compensatory_node_yields_no_corrective_direction(
    mechanism: MechanismHypothesis, catalog: DrugCatalog
) -> None:
    """Not every deviation from wild type should be pushed back.

    Clearance of aneuploid progenitors deviates from wild type and is protective;
    "correcting" it suppresses the response. No corrective sign follows from a
    compensatory node's state alone, so the check must return UNKNOWN and send the
    candidate to the "cannot determine" path — never award it agreement.
    """
    assert required_direction_for_node(mechanism, "N6_cellphen") is EffectDirection.DECREASE
    compensatory = _with_compensatory_node("N6_cellphen")
    assert required_direction_for_node(compensatory, "N6_cellphen") is EffectDirection.UNKNOWN

    # The same agent flips from an apparent agreement to an explicit unknown.
    suppressor = replace(
        _entry(catalog, "SYNTH-DRUG-E"),
        observed_direction=EffectDirection.DECREASE,
        worsens_cin=False,
    )
    as_pathological = generate_drug_hypotheses(
        mechanism=mechanism,
        catalog=DrugCatalog([suppressor], version=VERSION),
        clock=demo_clock(),
    ).all_hypotheses[0]
    as_compensatory = generate_drug_hypotheses(
        mechanism=compensatory,
        catalog=DrugCatalog([suppressor], version=VERSION),
        clock=demo_clock(),
    ).all_hypotheses[0]

    assert as_pathological.directions_agree is True
    assert as_compensatory.directions_agree is None
    assert as_compensatory.score < as_pathological.score


# non-fatal reasons are carried, not discarded --------------------------------


@pytest.mark.unit
def test_an_accepted_hypothesis_carries_its_non_fatal_reasons(result: DrugTriageResult) -> None:
    """GP-19: a reason that did not kill the candidate is still a finding.

    Computing `MECHANISM_MISMATCH` and then dropping it on the floor is what left
    the symptomatic/disease-modifying separation resting on a single curated cell.
    """
    synthazepam = _by_id(result, "SYNTH-DRUG-D")
    assert synthazepam.rejected is False
    assert synthazepam.rejection_reasons == ()
    assert RejectionReason.MECHANISM_MISMATCH in synthazepam.concerns
    # A reason is fatal or a concern, never both.
    for hypothesis in result.all_hypotheses:
        assert not (set(hypothesis.concerns) & set(hypothesis.rejection_reasons))
        if hypothesis.rejected:
            assert hypothesis.concerns == ()


@pytest.mark.unit
def test_acting_off_the_therapeutic_target_is_stated_whatever_the_declared_class(
    mechanism: MechanismHypothesis, catalog: DrugCatalog
) -> None:
    """The structural claim must not rest on a label.

    Flipping one TSV cell — Synthazepam's `intervention_class` — moved a
    phenotype-level agent into the disease-modifying section at rank 1, score 1.000,
    with nothing in the report saying it acts on the organismal phenotype. The
    statement is therefore gated on the target node, not on the class.
    """
    relabelled = replace(
        _entry(catalog, "SYNTH-DRUG-D"), intervention_class=InterventionClass.DISEASE_MODIFYING
    )
    triaged = generate_drug_hypotheses(
        mechanism=mechanism,
        catalog=DrugCatalog([relabelled], version=VERSION),
        clock=demo_clock(),
    )
    hypothesis = triaged.all_hypotheses[0]
    assert hypothesis.intervention_class is InterventionClass.DISEASE_MODIFYING
    assert hypothesis.target_node_id != mechanism.therapeutic_target_node_id
    assert RejectionReason.MECHANISM_MISMATCH in hypothesis.concerns

    statements = [
        item for item in triaged.evidence if "does NOT correct the mechanism" in item.claim
    ]
    assert len(statements) == 1
    assert "never be presented as disease-modifying" in statements[0].claim
    assert mechanism.therapeutic_target_node_id in statements[0].claim


# copies must re-validate -----------------------------------------------------


@pytest.mark.unit
def test_a_rejected_wrong_direction_drug_cannot_be_copied_into_acceptance(
    result: DrugTriageResult,
) -> None:
    """`model_copy(update=...)` does not re-run validators; `revalidated_copy` does.

    Without this the contraindicated compound can be flipped to `rejected=False`
    with no error, after which it matches no report section at all and disappears
    from the output entirely rather than being flagged.
    """
    synthinib = _by_id(result, "SYNTH-DRUG-B")
    assert synthinib.directions_agree is False and synthinib.rejected is True

    with pytest.raises(ValueError, match="wrong-direction agent must be rejected"):
        synthinib.revalidated_copy(rejected=False, rejection_reasons=())

    # The bypass still exists in pydantic itself, which is exactly why the pipeline
    # never uses it and the renderer re-checks.
    smuggled = synthinib.model_copy(update={"rejected": False, "rejection_reasons": ()})
    assert smuggled.rejected is False, "model_copy is expected to skip validation"

    # And the rank stamp the pipeline applies goes through the validating path.
    top = result.accepted[0]
    assert top.revalidated_copy(rank=7).rank == 7
