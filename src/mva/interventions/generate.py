"""Drug triage: turn a mechanism plus a catalogue into ranked, audited hypotheses.

This is the orchestration layer of Track 2. It takes a `MechanismHypothesis` as an
**argument** — it never reaches into the mechanism stage (GP-03) — walks the
catalogue, and produces two lists that are equally important:

* `accepted`, ranked from 1, each answering all eight mandatory questions through
  its fields; and
* `rejected`, kept in full with the reasons that killed them (GP-19). A pipeline
  that silently drops candidates cannot be reviewed, and the single most valuable
  output of this stage is the *rejection* of a plausible-looking, contraindicated
  agent. Deleting it would delete the finding.

Three separations are enforced deliberately, because collapsing any of them is a
known way to produce a confident, wrong recommendation:

* **wrong direction vs unknown direction.** The first is disqualifying, the second
  is a penalty and a flag. See `mva.interventions.direction`.
* **wrong direction vs not approved.** A tool compound pointing the right way is
  rejected for regulatory status, not for biology; a phase-1 inhibitor pointing
  the wrong way is rejected for biology first. The primary reason is chosen by a
  fixed priority order, so the report never blames the wrong thing.
* **symptomatic vs disease-modifying.** A symptom-management agent is scored under
  a class multiplier and carries an explicit evidence row stating that it does not
  correct the mechanism. It can appear in the accepted list; it can never be
  presented as mechanism correction.

NOTHING here is medical advice. Every row is a research hypothesis requiring
pre-clinical validation, and every accepted row names the experiment to run first.
"""

from __future__ import annotations

from dataclasses import dataclass

from mva.clock import Clock
from mva.interventions.catalog import CatalogEntry, DrugCatalog
from mva.interventions.direction import (
    DirectionVerdict,
    check_direction,
    check_target_in_mechanism,
    required_direction_for_node,
)
from mva.interventions.safety import (
    SEVERITY_PENALTY,
    SafetyVerdict,
    assess_evidence_quality,
    assess_safety,
)
from mva.models.base import AssertionTier
from mva.models.drug import (
    ApprovalStatus,
    DrugHypothesis,
    InterventionClass,
    PediatricEvidence,
    PharmacokineticProfile,
    RejectionReason,
    SafetyConcern,
)
from mva.models.evidence import (
    Citation,
    EvidenceCategory,
    EvidenceDirection,
    EvidenceItem,
    EvidenceStrength,
    EvidenceType,
    make_evidence_id,
)
from mva.models.mechanism import EffectDirection, MechanismHypothesis

__all__ = [
    "APPROVAL_STATUS_WEIGHT",
    "APPROVAL_WEIGHT",
    "DIRECTION_UNKNOWN_CREDIT",
    "DIRECTION_WEIGHT",
    "EVIDENCE_QUALITY_FLOOR",
    "EVIDENCE_WEIGHT",
    "INTERVENTION_CLASS_WEIGHT",
    "REASON_PRIORITY",
    "SAFETY_WEIGHT",
    "DrugTriageResult",
    "generate_drug_hypotheses",
]

TOOL_NAME = "mva.interventions.generate"
TOOL_VERSION = "0.1.0"

#: GP-20. Appended to every limitation this module writes.
SYNTHETIC_CATALOGUE_LIMITATION = (
    "Derived from a SYNTHETIC drug catalogue describing fictional agents. It is not "
    "pharmacological fact, carries no clinical validity, and is not medical advice."
)

# --------------------------------------------------------------------- weights
# GP-32: every number below is a documented heuristic, not a calibrated parameter.
# Changing one requires a decision record, a test, and a before/after comparison
# against tests/golden/expected_drug_outcomes.tsv. They live here rather than in
# config/default.yaml only until the config layer grows a drug-triage section;
# they are module constants so the diff of a change is visible.

#: Does the agent push the target the way the mechanism requires? Largest single
#: weight: direction is the axis on which naive repurposing fails.
DIRECTION_WEIGHT = 0.40
#: How far is the evidence from a demonstrated human outcome?
EVIDENCE_WEIGHT = 0.30
#: What did the safety pass find, weighted by severity?
SAFETY_WEIGHT = 0.20
#: Is the agent actually available to redeploy?
APPROVAL_WEIGHT = 0.10

#: Fraction of the direction term retained when agreement is undeterminable. Well
#: below 1.0 (an unknown is never an agreement) and well above 0.0 (an unknown is
#: not a disagreement either).
DIRECTION_UNKNOWN_CREDIT = 0.35

#: Below this evidence-quality score the candidate is flagged INSUFFICIENT_EVIDENCE.
#: A flag, not a filter (GP-13): weak evidence down-ranks, it does not delete.
EVIDENCE_QUALITY_FLOOR = 0.30

APPROVAL_STATUS_WEIGHT: dict[ApprovalStatus, float] = {
    ApprovalStatus.APPROVED: 1.00,
    ApprovalStatus.APPROVED_OTHER_INDICATION: 0.90,
    ApprovalStatus.APPROVED_ADULT_ONLY: 0.70,
    ApprovalStatus.INVESTIGATIONAL_PHASE_3: 0.35,
    ApprovalStatus.INVESTIGATIONAL_PHASE_2: 0.25,
    ApprovalStatus.INVESTIGATIONAL_PHASE_1: 0.15,
    ApprovalStatus.PRECLINICAL: 0.05,
    ApprovalStatus.TOOL_COMPOUND: 0.05,
    ApprovalStatus.WITHDRAWN: 0.00,
    ApprovalStatus.UNKNOWN: 0.00,
}

#: Multiplier keeping symptom management below mechanism correction in the ranking.
#: A symptomatic agent is a legitimate, sometimes vital, recommendation; it is not
#: a disease-modifying one, and the ordering must show that at a glance.
INTERVENTION_CLASS_WEIGHT: dict[InterventionClass, float] = {
    InterventionClass.DISEASE_MODIFYING: 1.00,
    InterventionClass.PREVENTIVE: 0.80,
    InterventionClass.SYMPTOMATIC: 0.50,
    InterventionClass.SUPPORTIVE: 0.45,
    InterventionClass.SURVEILLANCE: 0.40,
}

#: Order in which reasons are reported. The first reason a candidate raises is the
#: one a report leads with, so this ordering decides whether the reader is told
#: "wrong direction" or "not approved" about the same compound — they are different
#: findings with different follow-ups and must not be interchangeable.
REASON_PRIORITY: tuple[RejectionReason, ...] = (
    RejectionReason.TARGET_NOT_IN_MECHANISM,
    RejectionReason.WRONG_DIRECTION,
    RejectionReason.ONCOGENIC_RISK,
    RejectionReason.CONCENTRATION_NOT_ACHIEVABLE,
    RejectionReason.NOT_APPROVED,
    RejectionReason.SAFETY_CONCERN,
    RejectionReason.PHARMACOKINETIC_BARRIER,
    RejectionReason.MECHANISM_MISMATCH,
    RejectionReason.INSUFFICIENT_EVIDENCE,
    RejectionReason.NO_PEDIATRIC_EVIDENCE,
    RejectionReason.DIRECTION_UNKNOWN,
)


@dataclass(frozen=True)
class DrugTriageResult:
    """The complete, auditable output of one triage pass.

    `rejection_record` is the renderable ledger: one row for every candidate that
    raised at least one `RejectionReason`, fatal or not. An accepted candidate with
    an undetermined direction therefore appears there too, flagged, because "we
    accepted it despite X" is exactly the sentence a reviewer needs to see.
    """

    accepted: tuple[DrugHypothesis, ...]
    rejected: tuple[DrugHypothesis, ...]
    evidence: tuple[EvidenceItem, ...]
    rejection_record: tuple[dict[str, object], ...]

    @property
    def all_hypotheses(self) -> tuple[DrugHypothesis, ...]:
        """Accepted then rejected, in a total order (GP-30)."""
        return tuple(sorted([*self.accepted, *self.rejected], key=lambda d: d.sort_key()))


def generate_drug_hypotheses(
    *, mechanism: MechanismHypothesis, catalog: DrugCatalog, clock: Clock
) -> DrugTriageResult:
    """Triage every catalogue entry against one mechanism.

    Every entry is evaluated; nothing is pre-filtered by target, because a
    candidate whose target is off-chain is a *finding* (`TARGET_NOT_IN_MECHANISM`)
    rather than a non-event, and because filtering before scoring would hide the
    near misses that make the ranking legible (GP-13).
    """
    citation = Citation(
        source="mva-knowledge/drug_catalog.tsv",
        identifier=mechanism.mechanism_id,
        version=catalog.version,
        title="Synthetic drug catalogue",
    )
    triaged = [
        _triage_entry(entry, mechanism=mechanism, clock=clock, citation=citation)
        for entry in catalog.entries()
    ]

    accepted = sorted(
        (t for t in triaged if not t.hypothesis.rejected), key=lambda t: t.hypothesis.sort_key()
    )
    rejected = sorted(
        (t for t in triaged if t.hypothesis.rejected), key=lambda t: t.hypothesis.sort_key()
    )
    ranked = tuple(
        t.hypothesis.model_copy(update={"rank": position})
        for position, t in enumerate(accepted, start=1)
    )
    evidence: list[EvidenceItem] = []
    for t in triaged:
        evidence.extend(t.evidence)
    record = tuple(t.record for t in triaged if t.record is not None)
    return DrugTriageResult(
        accepted=ranked,
        rejected=tuple(t.hypothesis for t in rejected),
        evidence=tuple(evidence),
        rejection_record=record,
    )


# -------------------------------------------------------------------- internals


@dataclass(frozen=True)
class _Triaged:
    hypothesis: DrugHypothesis
    evidence: tuple[EvidenceItem, ...]
    record: dict[str, object] | None


def _clamp01(value: float) -> float:
    return max(0.0, min(1.0, value))


def _ordered(reasons: list[RejectionReason]) -> tuple[RejectionReason, ...]:
    """Deduplicate and sort reasons by reporting priority (GP-30: total order)."""
    unique = set(reasons)
    known = [reason for reason in REASON_PRIORITY if reason in unique]
    unlisted = sorted((r for r in unique if r not in REASON_PRIORITY), key=lambda r: r.value)
    return (*known, *unlisted)


def _strength_for(score: float) -> EvidenceStrength:
    if score >= 0.80:
        return EvidenceStrength.STRONG
    if score >= 0.60:
        return EvidenceStrength.MODERATE
    if score >= 0.40:
        return EvidenceStrength.SUPPORTING
    if score >= 0.20:
        return EvidenceStrength.WEAK
    return EvidenceStrength.INSUFFICIENT


def _safety_term(concerns: tuple[SafetyConcern, ...]) -> float:
    penalty = sum(SEVERITY_PENALTY.get(concern.severity, 0.20) for concern in concerns)
    return _clamp01(1.0 - penalty)


def _score(
    entry: CatalogEntry,
    *,
    verdict: DirectionVerdict,
    safety: SafetyVerdict,
    evidence_quality: float,
) -> float:
    """Composite score in ``[0, 1]``.

    ``(direction, evidence, safety, approval)`` weighted, then multiplied by the
    intervention-class weight so symptom management cannot outrank mechanism
    correction on evidence strength alone — which it otherwise would, since a
    seizure drug has human-trial data and a mechanism-correcting agent rarely does.
    """
    if verdict.agrees is True:
        direction_term = 1.0
    elif verdict.agrees is None:
        direction_term = DIRECTION_UNKNOWN_CREDIT
    else:
        direction_term = 0.0
    raw = (
        DIRECTION_WEIGHT * direction_term
        + EVIDENCE_WEIGHT * evidence_quality
        + SAFETY_WEIGHT * _safety_term(safety.concerns)
        + APPROVAL_WEIGHT * APPROVAL_STATUS_WEIGHT.get(entry.approval_status, 0.0)
    )
    class_weight = INTERVENTION_CLASS_WEIGHT.get(entry.intervention_class, 0.5)
    return round(_clamp01(raw * class_weight), 6)


def _evidence_item(
    *,
    subject_id: str,
    claim: str,
    category: EvidenceCategory,
    direction: EvidenceDirection,
    strength: EvidenceStrength,
    evidence_type: EvidenceType,
    tier: AssertionTier,
    method: str,
    limitations: str,
    clock: Clock,
    citation: Citation | None,
    numeric_value: float | None = None,
    payload: dict[str, str | int | float | bool | None] | None = None,
) -> EvidenceItem:
    return EvidenceItem(
        evidence_id=make_evidence_id(
            subject_id=subject_id, category=category, claim=claim, tool=TOOL_NAME
        ),
        subject_id=subject_id,
        subject_kind="drug",
        claim=claim,
        category=category,
        direction=direction,
        strength=strength,
        evidence_type=evidence_type,
        tier=tier,
        citation=citation,
        method=method,
        tool=TOOL_NAME,
        tool_version=TOOL_VERSION,
        limitations=f"{limitations} {SYNTHETIC_CATALOGUE_LIMITATION}",
        timestamp=clock.now(),
        numeric_value=numeric_value,
        payload=payload or {},
    )


def _triage_entry(
    entry: CatalogEntry,
    *,
    mechanism: MechanismHypothesis,
    clock: Clock,
    citation: Citation,
) -> _Triaged:
    target_in_mechanism = check_target_in_mechanism(entry.target_node_id, mechanism)
    required = (
        required_direction_for_node(mechanism, entry.target_node_id)
        if target_in_mechanism
        else EffectDirection.UNKNOWN
    )
    verdict = check_direction(required=required, observed=entry.observed_direction)
    safety = assess_safety(entry, mechanism=mechanism, clock=clock)
    evidence_quality, evidence_rationale = assess_evidence_quality(entry)
    score = _score(entry, verdict=verdict, safety=safety, evidence_quality=evidence_quality)

    reasons: list[RejectionReason] = []
    fatal = False
    if not target_in_mechanism:
        reasons.append(RejectionReason.TARGET_NOT_IN_MECHANISM)
        fatal = True
    elif entry.target_node_id != mechanism.therapeutic_target_node_id:
        reasons.append(RejectionReason.MECHANISM_MISMATCH)
    if verdict.rejection_reason is not None:
        reasons.append(verdict.rejection_reason)
        fatal = fatal or verdict.agrees is False
    if not entry.is_repurposable:
        reasons.append(RejectionReason.NOT_APPROVED)
        fatal = True
    reasons.extend(safety.reasons)
    fatal = fatal or safety.disqualifying
    if evidence_quality < EVIDENCE_QUALITY_FLOOR:
        reasons.append(RejectionReason.INSUFFICIENT_EVIDENCE)

    ordered = _ordered(reasons)
    primary = ordered[0] if ordered else None
    rationale = _rationale(
        entry,
        rejected=fatal,
        ordered=ordered,
        verdict=verdict,
        safety=safety,
        target_in_mechanism=target_in_mechanism,
        mechanism=mechanism,
    )

    items = _drug_evidence(
        entry,
        mechanism=mechanism,
        verdict=verdict,
        safety=safety,
        evidence_quality=evidence_quality,
        evidence_rationale=evidence_rationale,
        score=score,
        rejected=fatal,
        ordered=ordered,
        rationale=rationale,
        clock=clock,
        citation=citation,
    )
    concern_ids = {
        concern.concern_id: item.evidence_id
        for concern, item in zip(safety.concerns, _concern_items(items, safety), strict=True)
    }
    linked_concerns = tuple(
        concern.model_copy(update={"evidence_ids": (concern_ids[concern.concern_id],)})
        for concern in safety.concerns
    )

    hypothesis = DrugHypothesis(
        drug_id=entry.drug_id,
        name=entry.name,
        approved_name=entry.approved_name if entry.is_repurposable else None,
        approval_status=entry.approval_status,
        intervention_class=entry.intervention_class,
        target=entry.target,
        target_node_id=entry.target_node_id,
        mechanism_of_action=entry.mechanism_of_action,
        required_direction=required,
        observed_direction=entry.observed_direction,
        is_direct_evidence=entry.is_direct_evidence,
        strongest_evidence_type=entry.strongest_evidence_type,
        pediatric_evidence=PediatricEvidence(
            has_pediatric_exposure=entry.has_pediatric_exposure,
            youngest_age_studied=entry.youngest_age_studied,
            indication=entry.pediatric_indication,
            tolerability_summary=None,
        ),
        pharmacokinetics=PharmacokineticProfile(
            route=entry.route,
            cns_penetrant=entry.cns_penetrant,
            achievable_plasma_concentration_um=entry.achievable_plasma_um,
            required_effective_concentration_um=entry.required_effective_um,
            half_life_hours=entry.half_life_hours,
            notes=_pk_notes(entry),
        ),
        safety_concerns=linked_concerns,
        worsens_chromosomal_instability=entry.worsens_cin,
        evidence_ids=tuple(item.evidence_id for item in items),
        contradicting_evidence_ids=tuple(
            item.evidence_id for item in items if item.direction is EvidenceDirection.CONTRADICTS
        ),
        proposed_validation_experiment=entry.validation_experiment,
        score=score,
        rank=None,
        rejected=fatal,
        rejection_reasons=ordered if fatal else (),
        rejection_rationale=rationale if fatal else "",
    )
    record = (
        _record(
            entry,
            hypothesis=hypothesis,
            ordered=ordered,
            primary=primary,
            rationale=rationale,
            evidence_quality=evidence_quality,
            target_in_mechanism=target_in_mechanism,
        )
        if ordered
        else None
    )
    return _Triaged(hypothesis=hypothesis, evidence=items, record=record)


def _concern_items(
    items: tuple[EvidenceItem, ...], safety: SafetyVerdict
) -> tuple[EvidenceItem, ...]:
    """The safety rows of `items`, in the same order as `safety.concerns`."""
    safety_items = tuple(item for item in items if item.category is EvidenceCategory.SAFETY)
    if len(safety_items) != len(safety.concerns):  # pragma: no cover - internal invariant
        msg = (
            f"Safety evidence rows ({len(safety_items)}) do not match concerns "
            f"({len(safety.concerns)}); the audit trail would be wrong."
        )
        raise ValueError(msg)
    return safety_items


def _pk_notes(entry: CatalogEntry) -> str:
    achievable = entry.concentration_achievable
    if achievable is True:
        return (
            f"Peak plasma {entry.achievable_plasma_um} uM against a required "
            f"{entry.required_effective_um} uM: reachable on the recorded figures."
        )
    if achievable is False:
        return (
            f"Peak plasma {entry.achievable_plasma_um} uM against a required "
            f"{entry.required_effective_um} uM: NOT reachable."
        )
    return (
        "Exposure relative to the effective concentration cannot be computed; a figure is missing."
    )


def _rationale(
    entry: CatalogEntry,
    *,
    rejected: bool,
    ordered: tuple[RejectionReason, ...],
    verdict: DirectionVerdict,
    safety: SafetyVerdict,
    target_in_mechanism: bool,
    mechanism: MechanismHypothesis,
) -> str:
    parts: list[str] = []
    if ordered:
        listed = ", ".join(reason.value for reason in ordered)
        verb = "REJECTED" if rejected else "ACCEPTED WITH CONCERNS"
        parts.append(
            f"{verb} ({entry.drug_id}, {entry.name}). Reasons, in priority order: {listed}."
        )
    if not target_in_mechanism:
        parts.append(
            f"Target node {entry.target_node_id} is not among the {len(mechanism.nodes)} nodes of "
            f"{mechanism.mechanism_id}, so no direction check is even possible."
        )
    elif entry.target_node_id != mechanism.therapeutic_target_node_id:
        parts.append(
            f"Acts on {entry.target_node_id}, not on the designated therapeutic target "
            f"{mechanism.therapeutic_target_node_id}; it does not correct the mechanism."
        )
    parts.append(verdict.rationale)
    if not entry.is_repurposable:
        parts.append(
            f"Approval status is {entry.approval_status.value}: not an approved agent, so it is "
            "not a repurposing candidate regardless of its direction of effect."
        )
    parts.append(safety.rationale)
    return " ".join(parts)


def _drug_evidence(
    entry: CatalogEntry,
    *,
    mechanism: MechanismHypothesis,
    verdict: DirectionVerdict,
    safety: SafetyVerdict,
    evidence_quality: float,
    evidence_rationale: str,
    score: float,
    rejected: bool,
    ordered: tuple[RejectionReason, ...],
    rationale: str,
    clock: Clock,
    citation: Citation,
) -> tuple[EvidenceItem, ...]:
    """One evidence row per finding. GP-10: nothing below is asserted uncited."""
    items: list[EvidenceItem] = []

    if verdict.agrees is True:
        direction_evidence = EvidenceDirection.SUPPORTS
        direction_strength = EvidenceStrength.MODERATE
        direction_limits = (
            "Sign agreement is not efficacy: it says the agent pushes the right way, not that "
            "pushing helps, nor by how much."
        )
    elif verdict.agrees is False:
        direction_evidence = EvidenceDirection.CONTRADICTS
        direction_strength = EvidenceStrength.STRONG
        direction_limits = (
            "Establishes that the agent is contraindicated on direction grounds; it does not "
            "establish the magnitude of harm."
        )
    else:
        direction_evidence = EvidenceDirection.NEUTRAL
        direction_strength = EvidenceStrength.INSUFFICIENT
        direction_limits = (
            "Records an absence of a determinable sign. It is neither support nor contradiction "
            "and must not be reported as either (GP-14, GP-16)."
        )
    items.append(
        _evidence_item(
            subject_id=entry.drug_id,
            claim=(
                f"Direction check for {entry.name} on {entry.target_node_id}: {verdict.rationale}"
            ),
            category=EvidenceCategory.DRUG,
            direction=direction_evidence,
            strength=direction_strength,
            evidence_type=EvidenceType.CURATED_DATABASE,
            tier=AssertionTier.DATABASE_ASSERTION,
            method=(
                "Compared the mechanism's required correction at the target node with the "
                "catalogue's signed observed direction (mva.interventions.direction)."
            ),
            limitations=direction_limits,
            clock=clock,
            citation=citation,
            payload={
                "required_direction": verdict.required.value,
                "observed_direction": verdict.observed.value,
                "directions_agree": verdict.agrees,
                "mechanism_id": mechanism.mechanism_id,
            },
        )
    )

    items.append(
        _evidence_item(
            subject_id=entry.drug_id,
            claim=evidence_rationale,
            category=EvidenceCategory.DRUG,
            direction=EvidenceDirection.NEUTRAL,
            strength=_strength_for(evidence_quality),
            evidence_type=entry.strongest_evidence_type,
            tier=AssertionTier.DATABASE_ASSERTION,
            method="Evidence-tier weighting with a discount for indirect evidence.",
            limitations=(
                "Grades the distance from a demonstrated human outcome. A high tier for a "
                "different indication says nothing about efficacy in this mechanism."
            ),
            clock=clock,
            citation=citation,
            numeric_value=evidence_quality,
        )
    )

    for concern in safety.concerns:
        items.append(
            _evidence_item(
                subject_id=entry.drug_id,
                claim=f"Safety concern {concern.concern_id}: {concern.description}",
                category=EvidenceCategory.SAFETY,
                direction=EvidenceDirection.CONTRADICTS,
                strength=(
                    EvidenceStrength.STRONG
                    if concern.severity == "critical"
                    else EvidenceStrength.MODERATE
                ),
                evidence_type=EvidenceType.CURATED_DATABASE,
                tier=AssertionTier.DATABASE_ASSERTION,
                method="Safety pass over the catalogue fields, in this disease context.",
                limitations=(
                    "Reflects only the fields the catalogue records. The absence of a concern "
                    "is not a safety clearance, and no interaction or long-term toxicity data "
                    "was consulted."
                ),
                clock=clock,
                citation=citation,
                payload={
                    "concern_id": concern.concern_id,
                    "severity": concern.severity,
                    "is_disqualifying": concern.is_disqualifying,
                },
            )
        )

    if entry.intervention_class is not InterventionClass.DISEASE_MODIFYING:
        items.append(
            _evidence_item(
                subject_id=entry.drug_id,
                claim=(
                    f"{entry.name} is classified {entry.intervention_class.value}: it acts on "
                    f"{entry.target_node_id} and does NOT correct the mechanism at "
                    f"{mechanism.therapeutic_target_node_id}. It must never be presented as "
                    "disease-modifying."
                ),
                category=EvidenceCategory.DRUG,
                direction=EvidenceDirection.NEUTRAL,
                strength=EvidenceStrength.MODERATE,
                evidence_type=EvidenceType.PIPELINE_INFERENCE,
                tier=AssertionTier.INFERENCE,
                method="Compared the catalogue's intervention class and target node with the "
                "mechanism's therapeutic target node.",
                limitations=(
                    "A classification, not an efficacy claim. Symptom control can be the right "
                    "clinical priority; it is simply a different claim from mechanism correction."
                ),
                clock=clock,
                citation=citation,
            )
        )

    listed = ", ".join(reason.value for reason in ordered) if ordered else "none"
    if rejected:
        triage_claim = (
            f"{entry.drug_id} ({entry.name}) was REJECTED with reasons [{listed}] and is retained "
            f"with its rationale at score {score:.3f}."
        )
    else:
        triage_claim = (
            f"{entry.drug_id} ({entry.name}) was ACCEPTED as a hypothesis at score {score:.3f}, "
            f"carrying [{listed}] as recorded concerns."
        )
    items.append(
        _evidence_item(
            subject_id=entry.drug_id,
            claim=triage_claim,
            category=(EvidenceCategory.CONTRADICTION if rejected else EvidenceCategory.DRUG),
            direction=(EvidenceDirection.CONTRADICTS if rejected else EvidenceDirection.NEUTRAL),
            strength=_strength_for(score),
            evidence_type=EvidenceType.PIPELINE_INFERENCE,
            tier=AssertionTier.INFERENCE,
            method=(
                f"score = clamp01({DIRECTION_WEIGHT}*direction + {EVIDENCE_WEIGHT}*evidence + "
                f"{SAFETY_WEIGHT}*safety + {APPROVAL_WEIGHT}*approval) * class_weight. "
                f"Rationale: {rationale}"
            ),
            limitations=(
                "A pipeline conclusion, never external support (GP-12). The score ranks "
                "hypotheses against each other; it is not a probability of benefit and cannot "
                "be compared across mechanisms."
            ),
            clock=clock,
            citation=None,
            numeric_value=score,
            payload={"rejected": rejected, "reasons": listed},
        )
    )
    return tuple(items)


def _record(
    entry: CatalogEntry,
    *,
    hypothesis: DrugHypothesis,
    ordered: tuple[RejectionReason, ...],
    primary: RejectionReason | None,
    rationale: str,
    evidence_quality: float,
    target_in_mechanism: bool,
) -> dict[str, object]:
    """One renderable audit row. Values are scalars so any renderer can print it."""
    return {
        "drug_id": entry.drug_id,
        "name": entry.name,
        "rejected": hypothesis.rejected,
        "primary_reason": primary.value if primary is not None else None,
        "reasons": ", ".join(reason.value for reason in ordered),
        "directions_agree": hypothesis.directions_agree,
        "required_direction": hypothesis.required_direction.value,
        "observed_direction": hypothesis.observed_direction.value,
        "target_node_id": entry.target_node_id,
        "target_in_mechanism": target_in_mechanism,
        "approval_status": entry.approval_status.value,
        "intervention_class": entry.intervention_class.value,
        "evidence_quality": evidence_quality,
        "n_safety_concerns": len(hypothesis.safety_concerns),
        "score": hypothesis.score,
        "rationale": rationale,
    }
