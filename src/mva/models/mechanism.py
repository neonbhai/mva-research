"""Mechanism hypothesis models.

A mechanism is represented as an explicit **chain of links**, each independently
evidenced and each carrying its own direction of effect:

    variant -> transcript effect -> protein effect -> cellular process -> phenotype

Storing the chain rather than a prose paragraph is what makes Track 2 auditable.
It lets the pipeline ask mechanical questions a paragraph cannot answer: which
link is only inferred? which node would a drug act on? in which direction must
that node be pushed?
"""

from __future__ import annotations

from enum import StrEnum
from typing import Self

from pydantic import Field, model_validator

from mva.models.base import AssertionTier, FrozenModel
from mva.models.evidence import EvidenceStrength


class MechanismNodeKind(StrEnum):
    """The rung of biological organisation a node sits on."""

    VARIANT = "variant"
    TRANSCRIPT = "transcript"
    PROTEIN = "protein"
    COMPLEX = "complex"
    CELLULAR_PROCESS = "cellular_process"
    CELLULAR_PHENOTYPE = "cellular_phenotype"
    TISSUE = "tissue"
    ORGANISMAL_PHENOTYPE = "organismal_phenotype"


class EffectDirection(StrEnum):
    """Signed direction of an effect.

    This enum is the single most load-bearing type in the Track 2 pipeline. The
    canonical failure mode of naive drug repurposing is to find a compound acting
    on the right node in the WRONG direction — e.g. proposing a checkpoint
    inhibitor for a checkpoint-deficiency disorder, because both are "about the
    spindle assembly checkpoint" and a target-proximity search cannot see signs.
    """

    INCREASE = "increase"
    DECREASE = "decrease"
    LOSS_OF_FUNCTION = "loss_of_function"
    GAIN_OF_FUNCTION = "gain_of_function"
    STABILISE = "stabilise"
    DESTABILISE = "destabilise"
    RESTORE = "restore"
    NO_CHANGE = "no_change"
    UNKNOWN = "unknown"
    CONTEXT_DEPENDENT = "context_dependent"
    """Genuinely bidirectional in the literature (e.g. taxane effects on
    missegregation depend on dose and baseline CIN mechanism). Must never be
    silently resolved to a single sign."""


#: Directions that reduce the activity/abundance of their node.
DOWNWARD_DIRECTIONS: frozenset[EffectDirection] = frozenset(
    {EffectDirection.DECREASE, EffectDirection.LOSS_OF_FUNCTION, EffectDirection.DESTABILISE}
)

#: Directions that raise the activity/abundance of their node.
UPWARD_DIRECTIONS: frozenset[EffectDirection] = frozenset(
    {
        EffectDirection.INCREASE,
        EffectDirection.GAIN_OF_FUNCTION,
        EffectDirection.STABILISE,
        EffectDirection.RESTORE,
    }
)

#: Directions from which no sign can be derived. Never treat these as agreement.
#:
#: ``NO_CHANGE`` is deliberately NOT here. A demonstrated null result is a signed,
#: established finding — usually the strongest evidence *against* the hypothesis
#: that this agent corrects this node. Filing it as "unknown" would award it the
#: undetermined-direction credit and a recommendation to go and measure the sign
#: that has already been measured.
UNSIGNED_DIRECTIONS: frozenset[EffectDirection] = frozenset(
    {EffectDirection.UNKNOWN, EffectDirection.CONTEXT_DEPENDENT}
)

#: A direction that was measured and found to be null.
NULL_DIRECTIONS: frozenset[EffectDirection] = frozenset({EffectDirection.NO_CHANGE})


def directions_agree(required: EffectDirection, observed: EffectDirection) -> bool | None:
    """Do a required correction and an observed drug effect point the same way?

    Returns ``None`` — explicitly not ``False`` — when either direction is unsigned.
    "We cannot tell" and "they disagree" are different findings and are scored
    differently; conflating them would let an unknown masquerade as a rejection or,
    worse, as agreement.
    """
    if required in UNSIGNED_DIRECTIONS or observed in UNSIGNED_DIRECTIONS:
        return None
    if observed in NULL_DIRECTIONS:
        # A measured null does not push the node in the required direction. This
        # is a disagreement with evidence behind it, not an open question.
        return False
    if required in NULL_DIRECTIONS:
        return None
    both_up = required in UPWARD_DIRECTIONS and observed in UPWARD_DIRECTIONS
    both_down = required in DOWNWARD_DIRECTIONS and observed in DOWNWARD_DIRECTIONS
    return both_up or both_down


#: The direction a therapy must push a node, given how that node deviates in the
#: patient. Defined here rather than in the interventions layer so that the
#: `MechanismHypothesis` validator and the drug direction check cannot drift apart
#: -- two copies of a sign table is one copy too many.
CORRECTIVE_DIRECTION: dict[EffectDirection, EffectDirection] = {
    EffectDirection.DECREASE: EffectDirection.INCREASE,
    EffectDirection.INCREASE: EffectDirection.DECREASE,
    EffectDirection.LOSS_OF_FUNCTION: EffectDirection.RESTORE,
    EffectDirection.GAIN_OF_FUNCTION: EffectDirection.DECREASE,
    EffectDirection.DESTABILISE: EffectDirection.STABILISE,
    EffectDirection.STABILISE: EffectDirection.DESTABILISE,
    EffectDirection.RESTORE: EffectDirection.UNKNOWN,
    EffectDirection.NO_CHANGE: EffectDirection.UNKNOWN,
    EffectDirection.UNKNOWN: EffectDirection.UNKNOWN,
    EffectDirection.CONTEXT_DEPENDENT: EffectDirection.UNKNOWN,
}


def corrective_direction(state: EffectDirection) -> EffectDirection:
    """The direction a therapy must push a node deviating by ``state``."""
    return CORRECTIVE_DIRECTION.get(state, EffectDirection.UNKNOWN)


def same_sign(a: EffectDirection, b: EffectDirection) -> bool:
    """Whether two signed directions point the same way."""
    if a in UNSIGNED_DIRECTIONS or b in UNSIGNED_DIRECTIONS:
        return a is b
    return (a in UPWARD_DIRECTIONS and b in UPWARD_DIRECTIONS) or (
        a in DOWNWARD_DIRECTIONS and b in DOWNWARD_DIRECTIONS
    )


class MechanismNode(FrozenModel):
    """One entity in the mechanism chain."""

    node_id: str
    kind: MechanismNodeKind
    label: str
    identifier: str | None = Field(default=None, description="External ID: HGNC, UniProt, GO, HPO.")
    state_in_patient: EffectDirection = Field(
        description="How this node deviates from wild type in the affected individual."
    )
    deviation_is_pathological: bool = Field(
        description=(
            "Whether this node's deviation is part of the DISEASE or part of a "
            "COMPENSATORY response. Mandatory and un-defaulted, because the naive "
            "assumption -- that every deviation from wild type should be pushed back "
            "-- inverts the sign on compensatory nodes. Suppressing a protective "
            "response, such as clearance of aneuploid progenitors, 'corrects' a "
            "deviation and is exactly the wrong thing to do. Where a node is "
            "compensatory, no corrective direction can be derived from its state "
            "alone and the direction check must return UNKNOWN rather than guess."
        )
    )
    description: str = ""


class MechanismLink(FrozenModel):
    """A directed, evidenced edge between two mechanism nodes."""

    link_id: str
    source_node_id: str
    target_node_id: str
    relation: str = Field(description="e.g. 'reduces_abundance_of', 'fails_to_inhibit'.")
    direction: EffectDirection
    tier: AssertionTier
    strength: EvidenceStrength
    evidence_ids: tuple[str, ...] = ()
    contradicting_evidence_ids: tuple[str, ...] = ()
    is_directly_demonstrated: bool = Field(
        description=(
            "True only if experimental work demonstrated THIS link, in a relevant "
            "system. False means inferred by analogy or pathway membership."
        )
    )
    uncertainty: str = Field(description="What is unresolved about this specific link.")


class DiscriminatingExperiment(FrozenModel):
    """An experiment that would separate this hypothesis from a stated alternative."""

    experiment_id: str
    description: str
    measures: str = Field(description="The concrete readout.")
    distinguishes_from: str = Field(description="The alternative hypothesis it rules out.")
    expected_if_true: str
    expected_if_false: str
    feasibility: str = Field(description="'routine' | 'specialised' | 'research_only'")


class MechanismHypothesis(FrozenModel):
    """A full, auditable mechanistic account linking a candidate pair to a phenotype."""

    mechanism_id: str
    gene_symbol: str
    pair_id: str | None = None
    summary: str = Field(description="One-paragraph plain statement of the mechanism.")
    nodes: tuple[MechanismNode, ...]
    links: tuple[MechanismLink, ...]

    disease_direction: EffectDirection = Field(
        description=(
            "Net direction of the pathological change in the key node. The therapeutic "
            "requirement is its inverse; this field is what the drug direction check "
            "compares against."
        )
    )
    therapeutic_target_node_id: str = Field(description="Which node an intervention should act on.")
    required_correction: EffectDirection = Field(
        description="The direction a therapy must push the target node."
    )

    supporting_evidence_ids: tuple[str, ...] = ()
    contradicting_evidence_ids: tuple[str, ...] = ()
    uncertainties: tuple[str, ...] = ()
    discriminating_experiments: tuple[DiscriminatingExperiment, ...] = ()
    developmental_window_caveat: str = Field(
        default="",
        description=(
            "Whether post-natal correction could plausibly help. Structural damage "
            "established in utero is not reversible by later target modulation, and a "
            "mechanistically correct drug can still be therapeutically irrelevant."
        ),
    )

    @property
    def node_ids(self) -> frozenset[str]:
        return frozenset(n.node_id for n in self.nodes)

    @property
    def inferred_links(self) -> tuple[MechanismLink, ...]:
        """Links that are NOT directly demonstrated — the weak points of the chain."""
        return tuple(link for link in self.links if not link.is_directly_demonstrated)

    @property
    def is_fully_demonstrated(self) -> bool:
        return bool(self.links) and not self.inferred_links

    def target_node(self) -> MechanismNode:
        for node in self.nodes:
            if node.node_id == self.therapeutic_target_node_id:
                return node
        msg = (
            f"Mechanism {self.mechanism_id!r} names therapeutic target "
            f"{self.therapeutic_target_node_id!r}, which is not among its nodes."
        )
        raise ValueError(msg)

    @model_validator(mode="after")
    def _direction_triple_is_consistent(self) -> Self:
        """The three curated direction cells must agree.

        ``disease_direction``, ``required_correction`` and the target node's
        ``state_in_patient`` are authored independently, and the drug direction
        check compares an agent's observed direction against
        ``required_correction`` alone. Without this validator a single mistyped
        cell silently inverts the entire Track 2 gate: the contraindicated
        checkpoint inhibitor is then reported as "direction AGREES" and the
        correct stabiliser is rejected as wrong-direction, with the report
        printing its own refutation and nothing reading it.

        So the invariant is enforced at construction, where it cannot be skipped:
        the required correction must be the corrective direction of both the
        disease direction and the target node's patient state.
        """
        if not self.nodes:
            return self

        target = self.target_node()

        if not target.deviation_is_pathological:
            msg = (
                f"Mechanism {self.mechanism_id!r} designates {target.node_id!r} as the "
                "therapeutic target, but that node is marked compensatory "
                "(deviation_is_pathological=False). Pushing a compensatory response "
                "back towards wild type suppresses a protective mechanism. Choose a "
                "pathological node as the target, or correct the node's flag."
            )
            raise ValueError(msg)

        expected_from_state = corrective_direction(target.state_in_patient)
        if expected_from_state is not EffectDirection.UNKNOWN and not same_sign(
            self.required_correction, expected_from_state
        ):
            msg = (
                f"Mechanism {self.mechanism_id!r} requires "
                f"{self.required_correction.value!r} at target {target.node_id!r}, but "
                f"that node's patient state is {target.state_in_patient.value!r}, whose "
                f"correction is {expected_from_state.value!r}. These point opposite ways. "
                "A sign error here inverts the entire drug direction check, so it is "
                "refused at construction rather than reported."
            )
            raise ValueError(msg)

        expected_from_disease = corrective_direction(self.disease_direction)
        if expected_from_disease is not EffectDirection.UNKNOWN and not same_sign(
            self.required_correction, expected_from_disease
        ):
            msg = (
                f"Mechanism {self.mechanism_id!r} declares disease_direction "
                f"{self.disease_direction.value!r}, whose correction is "
                f"{expected_from_disease.value!r}, but required_correction is "
                f"{self.required_correction.value!r}. These point opposite ways."
            )
            raise ValueError(msg)
        return self
