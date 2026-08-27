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

from pydantic import Field

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
UNSIGNED_DIRECTIONS: frozenset[EffectDirection] = frozenset(
    {EffectDirection.UNKNOWN, EffectDirection.CONTEXT_DEPENDENT, EffectDirection.NO_CHANGE}
)


def directions_agree(required: EffectDirection, observed: EffectDirection) -> bool | None:
    """Do a required correction and an observed drug effect point the same way?

    Returns ``None`` — explicitly not ``False`` — when either direction is unsigned.
    "We cannot tell" and "they disagree" are different findings and are scored
    differently; conflating them would let an unknown masquerade as a rejection or,
    worse, as agreement.
    """
    if required in UNSIGNED_DIRECTIONS or observed in UNSIGNED_DIRECTIONS:
        return None
    both_up = required in UPWARD_DIRECTIONS and observed in UPWARD_DIRECTIONS
    both_down = required in DOWNWARD_DIRECTIONS and observed in DOWNWARD_DIRECTIONS
    return both_up or both_down


class MechanismNode(FrozenModel):
    """One entity in the mechanism chain."""

    node_id: str
    kind: MechanismNodeKind
    label: str
    identifier: str | None = Field(default=None, description="External ID: HGNC, UniProt, GO, HPO.")
    state_in_patient: EffectDirection = Field(
        description="How this node deviates from wild type in the affected individual."
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
