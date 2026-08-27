"""Direction-of-effect checking. The load-bearing logic of Track 2.

The failure this module exists to prevent is specific and well documented: in a
checkpoint-*deficiency* disease, the compounds a target-proximity search ranks
highest are checkpoint *inhibitors*, because a deficiency disorder and an
oncology programme share a pathway vocabulary. They bind the right node and push
it further in the disease direction. Proximity is maximal; the sign is inverted.

So the check is not "is the target related?" but three separate questions, kept
separate on purpose:

1. **Is the target even in the mechanism?** If not, nothing can be said about
   direction and the candidate is out (`TARGET_NOT_IN_MECHANISM`). A compound
   scored against a node that is not on the chain is scored against nothing.
2. **Which way must this node move?** For the designated therapeutic target that
   is the mechanism's `required_correction`. For any other node on the chain it is
   the inverse of that node's state in the patient — push it back toward wild type.
3. **Does the agent move it that way?** Tri-state, via `directions_agree`.

The tri-state is the whole point. ``False`` (demonstrably opposite) is
disqualifying. ``None`` (either side unsigned) is *not* a disagreement and is
emphatically *not* an agreement — it is an admission of ignorance that costs the
candidate score and is recorded as `DIRECTION_UNKNOWN`. Collapsing the two would
either resurrect a contraindicated compound or reject an untested one, and both
errors are expensive in opposite directions.
"""

from __future__ import annotations

from dataclasses import dataclass

from mva.models.drug import RejectionReason
from mva.models.mechanism import (
    DOWNWARD_DIRECTIONS,
    UNSIGNED_DIRECTIONS,
    UPWARD_DIRECTIONS,
    EffectDirection,
    MechanismHypothesis,
    directions_agree,
)

__all__ = [
    "DirectionVerdict",
    "check_direction",
    "check_target_in_mechanism",
    "corrective_direction",
    "required_direction_for_node",
]

#: The therapeutic requirement implied by a node's state in the patient: what a
#: therapy would have to do to move that node back toward wild type. Deliberately
#: partial — every unsigned state maps to `UNKNOWN` rather than to a guess.
_CORRECTIVE_DIRECTION: dict[EffectDirection, EffectDirection] = {
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


@dataclass(frozen=True)
class DirectionVerdict:
    """The signed comparison between what is needed and what the agent does.

    `agrees` is tri-state and mirrors `mva.models.mechanism.directions_agree`:
    ``True`` same sign, ``False`` opposite sign, ``None`` undeterminable because
    at least one side is unsigned. `rejection_reason` is populated for both
    ``False`` (`WRONG_DIRECTION`, disqualifying) and ``None`` (`DIRECTION_UNKNOWN`,
    a recorded concern and a score penalty) — it is the *caller* that decides
    which reasons are fatal, and it must never treat the two alike.
    """

    agrees: bool | None
    required: EffectDirection
    observed: EffectDirection
    rationale: str
    rejection_reason: RejectionReason | None


def corrective_direction(state_in_patient: EffectDirection) -> EffectDirection:
    """The direction a therapy must push a node given its state in the patient.

    Returns `EffectDirection.UNKNOWN` for any unsigned state. Inventing a sign for
    a node nobody has characterised is exactly how an unknown becomes an
    apparent agreement.
    """
    return _CORRECTIVE_DIRECTION.get(state_in_patient, EffectDirection.UNKNOWN)


def check_target_in_mechanism(target_node_id: str, mechanism: MechanismHypothesis) -> bool:
    """Is the drug's target node actually on this mechanism chain?

    Exact node-ID membership. No pathway-adjacency fallback: "related to the
    pathway" is not a mechanism, and a candidate whose target is off-chain cannot
    be checked for direction at all.
    """
    return target_node_id in mechanism.node_ids


def required_direction_for_node(
    mechanism: MechanismHypothesis, target_node_id: str
) -> EffectDirection:
    """The correction required at `target_node_id`.

    For the mechanism's designated therapeutic target this is the curated
    `required_correction`. For any other node on the chain it is derived from that
    node's `state_in_patient`: an agent may legitimately act somewhere else on the
    chain — including on the organismal phenotype, which is what a symptomatic
    drug does — and it still has to push in the corrective direction for the node
    it actually acts on.

    A node that is not on the chain yields `EffectDirection.UNKNOWN`; the caller
    is expected to have already rejected it via `check_target_in_mechanism`.
    """
    if target_node_id == mechanism.therapeutic_target_node_id:
        return mechanism.required_correction
    for node in mechanism.nodes:
        if node.node_id == target_node_id:
            return corrective_direction(node.state_in_patient)
    return EffectDirection.UNKNOWN


def _sign_word(direction: EffectDirection) -> str:
    if direction in UPWARD_DIRECTIONS:
        return "upward"
    if direction in DOWNWARD_DIRECTIONS:
        return "downward"
    return "unsigned"


def check_direction(*, required: EffectDirection, observed: EffectDirection) -> DirectionVerdict:
    """Compare the required correction with the agent's observed effect.

    The rationale is written to be quoted verbatim into a report: it names both
    directions and their signs, so a reader can check the conclusion rather than
    trust it.
    """
    agrees = directions_agree(required, observed)
    required_sign = _sign_word(required)
    observed_sign = _sign_word(observed)

    if agrees is True:
        return DirectionVerdict(
            agrees=True,
            required=required,
            observed=observed,
            rationale=(
                f"The mechanism requires {required.value} ({required_sign}) at the target and "
                f"the agent acts {observed.value} ({observed_sign}); the signs agree."
            ),
            rejection_reason=None,
        )

    if agrees is False:
        return DirectionVerdict(
            agrees=False,
            required=required,
            observed=observed,
            rationale=(
                f"DISQUALIFYING: the mechanism requires {required.value} ({required_sign}) at "
                f"the target but the agent acts {observed.value} ({observed_sign}). It pushes "
                "the target in the disease direction."
            ),
            rejection_reason=RejectionReason.WRONG_DIRECTION,
        )

    unsigned: list[str] = []
    if required in UNSIGNED_DIRECTIONS:
        unsigned.append(f"the required correction is {required.value}")
    if observed in UNSIGNED_DIRECTIONS:
        unsigned.append(f"the agent's observed effect is {observed.value}")
    return DirectionVerdict(
        agrees=None,
        required=required,
        observed=observed,
        rationale=(
            f"UNDETERMINED: {' and '.join(unsigned)}, which carries no sign. This is NOT "
            "agreement and NOT disagreement; the sign must be established experimentally "
            "before the candidate can be assessed, and until then the candidate is penalised."
        ),
        rejection_reason=RejectionReason.DIRECTION_UNKNOWN,
    )
