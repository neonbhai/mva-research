"""Drug-repurposing hypothesis models.

Everything in this module exists to make one class of error impossible to commit
silently: **recommending a compound that acts on the right target in the wrong
direction**. In a spindle-assembly-checkpoint deficiency, the drugs most likely
to surface from a naive pathway or target-proximity search are checkpoint
*inhibitors* — developed precisely to push chromosomally-unstable cancer cells
past their tolerance ceiling. For a patient whose non-tumour cells already sit
above that ceiling with no reserve, those compounds push in the disease direction.

Proximity is high and the sign is inverted. A pipeline that ranks by pathway
membership without a signed direction term will confidently recommend the most
contraindicated compound class available. Hence: direction is a required,
validated field, and disagreement is disqualifying.

NOTHING in this module is medical advice. Every output is a research hypothesis
requiring pre-clinical validation.
"""

from __future__ import annotations

from enum import StrEnum
from typing import Self

from pydantic import Field, computed_field, model_validator

from mva.models.base import FrozenModel
from mva.models.evidence import IN_VIVO_TYPES, EvidenceType
from mva.models.mechanism import EffectDirection, directions_agree


class ApprovalStatus(StrEnum):
    """Regulatory status. Drives eligibility for a *repurposing* claim.

    Repurposing means redeploying an already-approved agent. A tool compound that
    never entered humans is a valid research probe but is NOT a repurposing
    candidate, and must not be presented as one.
    """

    APPROVED = "approved"
    APPROVED_OTHER_INDICATION = "approved_other_indication"
    APPROVED_ADULT_ONLY = "approved_adult_only"
    INVESTIGATIONAL_PHASE_3 = "investigational_phase_3"
    INVESTIGATIONAL_PHASE_2 = "investigational_phase_2"
    INVESTIGATIONAL_PHASE_1 = "investigational_phase_1"
    PRECLINICAL = "preclinical"
    TOOL_COMPOUND = "tool_compound"
    WITHDRAWN = "withdrawn"
    UNKNOWN = "unknown"


#: Statuses that qualify as genuine drug *repurposing*.
REPURPOSABLE_STATUSES: frozenset[ApprovalStatus] = frozenset(
    {
        ApprovalStatus.APPROVED,
        ApprovalStatus.APPROVED_OTHER_INDICATION,
        ApprovalStatus.APPROVED_ADULT_ONLY,
    }
)


class InterventionClass(StrEnum):
    """What the intervention is trying to achieve.

    ``SYMPTOMATIC`` vs ``DISEASE_MODIFYING`` is tracked explicitly because
    conflating symptom management with correction of the disease mechanism is a
    recurring and consequential category error in repurposing write-ups.
    """

    DISEASE_MODIFYING = "disease_modifying"
    SYMPTOMATIC = "symptomatic"
    SURVEILLANCE = "surveillance"
    SUPPORTIVE = "supportive"
    PREVENTIVE = "preventive"


class RejectionReason(StrEnum):
    """Why a candidate was rejected. Rejections are persisted, never discarded."""

    WRONG_DIRECTION = "wrong_direction"
    DIRECTION_UNKNOWN = "direction_unknown"
    NOT_APPROVED = "not_approved"
    INSUFFICIENT_EVIDENCE = "insufficient_evidence"
    SAFETY_CONCERN = "safety_concern"
    ONCOGENIC_RISK = "oncogenic_risk"
    NO_PEDIATRIC_EVIDENCE = "no_pediatric_evidence"
    PHARMACOKINETIC_BARRIER = "pharmacokinetic_barrier"
    TARGET_NOT_IN_MECHANISM = "target_not_in_mechanism"
    CONCENTRATION_NOT_ACHIEVABLE = "concentration_not_achievable"
    MECHANISM_MISMATCH = "mechanism_mismatch"


class SafetyConcern(FrozenModel):
    """A specific, sourced safety signal."""

    concern_id: str
    description: str
    severity: str = Field(description="'critical' | 'major' | 'moderate' | 'minor'")
    population: str = Field(description="Population in which it was observed.")
    evidence_ids: tuple[str, ...] = ()
    is_disqualifying: bool = False


class PharmacokineticProfile(FrozenModel):
    """PK facts that determine whether a mechanism can be reached in vivo."""

    route: str | None = None
    cns_penetrant: bool | None = Field(
        default=None,
        description="None means unknown. Critical when the phenotype is neurological.",
    )
    achievable_plasma_concentration_um: float | None = Field(default=None, ge=0.0)
    required_effective_concentration_um: float | None = Field(default=None, ge=0.0)
    half_life_hours: float | None = Field(default=None, ge=0.0)
    notes: str = ""

    @computed_field  # type: ignore[prop-decorator]
    @property
    def concentration_achievable(self) -> bool | None:
        """Is the effective concentration clinically reachable?

        ``None`` when either figure is unknown. A compound that works at 10 uM in
        culture but peaks at 0.1 uM in plasma is not a therapy, however elegant the
        mechanism.
        """
        if (
            self.achievable_plasma_concentration_um is None
            or self.required_effective_concentration_um is None
        ):
            return None
        return self.achievable_plasma_concentration_um >= self.required_effective_concentration_um


class PediatricEvidence(FrozenModel):
    """Exposure evidence in children."""

    has_pediatric_exposure: bool
    youngest_age_studied: str | None = None
    indication: str | None = None
    tolerability_summary: str | None = None
    evidence_ids: tuple[str, ...] = ()
    caveat: str = Field(
        default=(
            "Tolerability in a paediatric oncology population is not evidence of "
            "safety in a germline chromosomal-instability population."
        )
    )


class DrugHypothesis(FrozenModel):
    """A single, fully-audited drug-repurposing hypothesis.

    The eight mandatory questions from the project's scientific safeguards are each
    answered by a required field, so an incomplete hypothesis cannot be constructed:

    1. What node is modified?              -> ``target_node_id``
    2. In what direction must it move?     -> ``required_direction``
    3. Does the drug act that way?         -> ``observed_direction`` / ``directions_agree``
    4. What evidence tier?                 -> ``strongest_evidence_type``
    5. Concentration achievable?           -> ``pharmacokinetics.concentration_achievable``
    6. Paediatric exposure?                -> ``pediatric_evidence``
    7. Could it worsen CIN/cancer risk?    -> ``worsens_chromosomal_instability``
    8. What experiment comes first?        -> ``proposed_validation_experiment``
    """

    drug_id: str
    name: str
    approved_name: str | None = Field(default=None, description="INN/generic name, if approved.")
    approval_status: ApprovalStatus
    intervention_class: InterventionClass

    target: str = Field(description="Molecular target, e.g. 'TTK/MPS1 kinase'.")
    target_node_id: str = Field(description="Mechanism node this acts on. Must exist.")
    mechanism_of_action: str

    required_direction: EffectDirection = Field(
        description="Direction the mechanism says the target must be pushed."
    )
    observed_direction: EffectDirection = Field(
        description="Direction the drug actually pushes it, per the literature."
    )

    is_direct_evidence: bool = Field(
        description="True if evidence is for THIS drug on THIS target in a relevant system."
    )
    strongest_evidence_type: EvidenceType
    pediatric_evidence: PediatricEvidence
    pharmacokinetics: PharmacokineticProfile
    safety_concerns: tuple[SafetyConcern, ...] = ()
    worsens_chromosomal_instability: bool | None = Field(
        default=None,
        description=(
            "Would this agent increase aneuploidy or cancer susceptibility? None means "
            "unassessed, which is itself a blocking gap for this disease context."
        ),
    )

    evidence_ids: tuple[str, ...] = ()
    contradicting_evidence_ids: tuple[str, ...] = ()
    proposed_validation_experiment: str = Field(
        min_length=10,
        description="The concrete pre-clinical experiment to run BEFORE any clinical thought.",
    )

    score: float = Field(default=0.0, ge=0.0, le=1.0)
    rank: int | None = Field(default=None, ge=1)
    rejected: bool = False
    rejection_reasons: tuple[RejectionReason, ...] = ()
    rejection_rationale: str = ""
    concerns: tuple[RejectionReason, ...] = Field(
        default=(),
        description=(
            "Reasons raised against this hypothesis that were NOT fatal. Held in a "
            "field of their own rather than merged into `rejection_reasons`, which "
            "an accepted hypothesis must leave empty: without somewhere to put them "
            "the triage stage computes the near misses and then discards them, and "
            "a candidate is presented as clean when it is merely not disqualified. "
            "A reason is either fatal or a concern, never both."
        ),
    )

    # ---------------------------------------------------------------- derived

    @computed_field  # type: ignore[prop-decorator]
    @property
    def directions_agree(self) -> bool | None:
        """Tri-state. ``None`` means undeterminable, which is NOT agreement."""
        return directions_agree(self.required_direction, self.observed_direction)

    @property
    def is_repurposable(self) -> bool:
        return self.approval_status in REPURPOSABLE_STATUSES

    @property
    def has_in_vivo_evidence(self) -> bool:
        return self.strongest_evidence_type in IN_VIVO_TYPES

    @property
    def disqualifying_safety(self) -> tuple[SafetyConcern, ...]:
        return tuple(c for c in self.safety_concerns if c.is_disqualifying)

    @property
    def evidence_is_indirect_only(self) -> bool:
        return not self.is_direct_evidence

    # ---------------------------------------------------------------- validators

    @model_validator(mode="after")
    def _rejection_bookkeeping(self) -> Self:
        if self.rejected and not self.rejection_reasons:
            msg = (
                f"Drug {self.drug_id!r} is marked rejected but lists no rejection reason. "
                "Silent rejection destroys the audit trail; record the reason."
            )
            raise ValueError(msg)
        if not self.rejected and self.rejection_reasons:
            msg = (
                f"Drug {self.drug_id!r} carries rejection reasons "
                f"{[r.value for r in self.rejection_reasons]} but rejected=False. "
                "A non-fatal reason belongs in `concerns`."
            )
            raise ValueError(msg)
        overlap = sorted(r.value for r in set(self.concerns) & set(self.rejection_reasons))
        if overlap:
            msg = (
                f"Drug {self.drug_id!r} lists {overlap} as both a rejection reason and a "
                "non-fatal concern. A reason is one or the other; recording it as both "
                "makes the rejection record contradict the candidate block."
            )
            raise ValueError(msg)
        return self

    @model_validator(mode="after")
    def _wrong_direction_must_be_rejected(self) -> Self:
        """A hard invariant, enforced at the type boundary rather than downstream.

        If the drug demonstrably pushes the target the wrong way, it cannot exist as
        an accepted hypothesis. This is deliberately not left to the scoring layer:
        a weight change must never be able to resurrect a contraindicated compound.
        """
        if self.directions_agree is False and not self.rejected:
            msg = (
                f"Drug {self.drug_id!r} ({self.name}) requires "
                f"{self.required_direction.value!r} but acts "
                f"{self.observed_direction.value!r} on {self.target!r}, yet is not "
                "rejected. A wrong-direction agent must be rejected with "
                "RejectionReason.WRONG_DIRECTION."
            )
            raise ValueError(msg)
        return self

    @model_validator(mode="after")
    def _non_approved_cannot_be_called_repurposing(self) -> Self:
        if not self.is_repurposable and self.approved_name is not None:
            msg = (
                f"Drug {self.drug_id!r} has approval_status "
                f"{self.approval_status.value!r} but carries an approved_name. A tool "
                "compound or investigational agent has no approved name."
            )
            raise ValueError(msg)
        return self

    def sort_key(self) -> tuple[int, float, str]:
        """Rejected candidates always sort after accepted ones, then by score desc."""
        return (1 if self.rejected else 0, -self.score, self.drug_id)
