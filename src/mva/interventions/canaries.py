"""A live positive control for the drug-triage pipeline.

`mva.privacy.patterns.detector_canaries` gives every privacy rule a specimen it
MUST match, and `mva.privacy.audit.check_detector_liveness` pushes those specimens
through the real scanner on every run. The argument for it is that a detector which
has silently stopped detecting emits exactly the report a clean repository does.
From the outside the two are indistinguishable.

**The science side had no equivalent, and the asymmetry is dangerous in exactly the
same shape.** Every existing Track 2 test asserts that something is *rejected*: the
MPS1 inhibitor comes last, the highest-scoring agent does not win, "cannot
determine" is never rounded to agreement. Those are the right tests and they share
one blind spot — a triage stage that rejected literally everything would pass all of
them. A pipeline whose only demonstrated behaviour is refusal is not a careful
pipeline; it is an unfalsifiable one, and it scores identically to a working one.

So this module injects agents whose verdict is not in question and demands the
pipeline reach it, through the same `generate_drug_hypotheses` every real run calls:

``ACCEPT``
    A synthetic agent that pushes the therapeutic target the way the mechanism
    requires, on direct human-trial evidence, approved, with paediatric exposure,
    an achievable exposure, and — the part that usually blocks — a chromosomal
    instability question that was **asked and answered negative**. Nothing about it
    is borderline. If this is rejected, the pipeline cannot say yes to anything.

``INVERT``
    The same agent with its direction flipped and nothing else changed. It must be
    rejected, and rejected for `WRONG_DIRECTION` specifically. Rejecting it for
    some other reason would pass a naive test while proving the sign check dead.

``UNKNOWN_DIRECTION``
    The same agent with an unsigned direction, everything else unchanged. It is
    **accepted**, carrying `DIRECTION_UNKNOWN` as a visible concern and a score
    penalty -- and that is the point. Wrong direction is disqualifying; unknown
    direction is not (GP-16, ASSUMPTION-DRUG-02). A pipeline that rejected this
    would be treating "nobody measured it" as "it points the wrong way", which is
    the GP-14 error in its other, quieter direction. The canary pins both halves:
    accepted, and accepted *with the flag attached*, scoring below the agent whose
    sign is known.

``UNASSESSED_RISK``
    The accepted agent with one cell blanked: nobody asked the chromosomal
    instability question. Rejected for `ONCOGENIC_RISK` while its direction still
    AGREES — which is the GP-14 boundary made executable. The rejection must come
    from the unasked question, not from a direction the pipeline quietly downgraded.

``OFF_CHAIN``
    A negative control acting on a node that is not on the mechanism at all.
    `TARGET_NOT_IN_MECHANISM`, so that "accepted" cannot be an artefact of the
    triage stage ignoring the mechanism it was handed.

``SIBLING_GENE``
    The accepted agent, unchanged, against a **different synthetic gene** whose
    disease direction is the mirror image. The same drug must now be rejected for
    `WRONG_DIRECTION`. This is the control that proves the verdict is a property of
    the mechanism-plus-drug pair rather than of the drug: a pipeline that had
    memorised "this compound is good" passes every other case here and fails this
    one.

GP-20. Everything below is **fabricated**. `CANARYGENE1`, `CANARYGENE2`,
`Canarystatin` and `Canaryase` are not genes, proteins or compounds; the node
identifiers are not ontology terms; the "evidence" describes experiments nobody
ran. They are shaped to be unambiguous to the *pipeline*, not plausible as biology,
and they are never written into a knowledge table, an evidence store or a report.
Nothing here is biologically valid and nothing here is medical advice.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

from mva.clock import demo_clock
from mva.interventions.catalog import CatalogEntry, DrugCatalog
from mva.interventions.generate import generate_drug_hypotheses
from mva.models.base import AssertionTier
from mva.models.drug import ApprovalStatus, InterventionClass, RejectionReason
from mva.models.evidence import EvidenceStrength, EvidenceType
from mva.models.mechanism import (
    EffectDirection,
    MechanismHypothesis,
    MechanismLink,
    MechanismNode,
    MechanismNodeKind,
)

__all__ = [
    "CANARY_VERSION",
    "CanaryFailure",
    "CanaryOutcome",
    "PipelineCanary",
    "TriageLiveness",
    "canary_mechanism",
    "check_triage_liveness",
    "pipeline_canaries",
]

#: Stamped on the synthetic catalogue so no artifact can confuse a canary run with
#: `bub1b-literature-v1` or with the synthetic demo tables (GP-20, GP-31).
CANARY_VERSION: Final[str] = "pipeline-canary-SYNTHETIC-v1"

_SYNTHETIC = (
    "SYNTHETIC PIPELINE CANARY. Fabricated for the sole purpose of proving the "
    "triage stage still reaches a verdict. Not biology, not pharmacology, not "
    "medical advice."
)


def canary_mechanism(
    *,
    mechanism_id: str,
    gene_symbol: str,
    disease_direction: EffectDirection,
) -> MechanismHypothesis:
    """A two-node fabricated chain in a chromosomal-instability context.

    The CIN context is deliberate rather than decorative: the unassessed-risk gate
    (ADR 0011) only bites inside one, and that gate is the pipeline behaviour most
    worth proving is still alive in both directions -- it must block the agent whose
    risk nobody measured, and it must NOT block the agent whose risk was measured
    and came back negative.

    ``required_correction`` is derived from ``disease_direction`` rather than passed
    in, so a caller cannot construct a canary whose own expected answer is wrong.
    """
    corrective = {
        EffectDirection.DECREASE: EffectDirection.INCREASE,
        EffectDirection.INCREASE: EffectDirection.DECREASE,
    }[disease_direction]
    prefix = mechanism_id.rsplit("-", 1)[-1]
    target_id = f"{prefix}_protein"
    downstream_id = f"{prefix}_aneuploidy"
    return MechanismHypothesis(
        mechanism_id=mechanism_id,
        gene_symbol=gene_symbol,
        summary=(
            f"{_SYNTHETIC} A fabricated chain in which {gene_symbol} protein abundance "
            f"shows a {disease_direction.value} and downstream aneuploidy and chromosome "
            f"missegregation rise. It describes no real gene and no real disease."
        ),
        nodes=(
            MechanismNode(
                node_id=target_id,
                kind=MechanismNodeKind.PROTEIN,
                label=f"{gene_symbol} protein abundance (SYNTHETIC)",
                identifier=None,
                state_in_patient=disease_direction,
                deviation_is_pathological=True,
                description=_SYNTHETIC,
            ),
            MechanismNode(
                node_id=downstream_id,
                kind=MechanismNodeKind.CELLULAR_PHENOTYPE,
                label="Chromosome missegregation and aneuploidy (SYNTHETIC)",
                identifier=None,
                state_in_patient=EffectDirection.INCREASE,
                deviation_is_pathological=True,
                description=_SYNTHETIC,
            ),
        ),
        links=(
            MechanismLink(
                link_id=f"{prefix}_L1",
                source_node_id=target_id,
                target_node_id=downstream_id,
                relation="fabricated_link_for_a_pipeline_canary",
                direction=EffectDirection.INCREASE,
                tier=AssertionTier.INFERENCE,
                strength=EvidenceStrength.MODERATE,
                is_directly_demonstrated=False,
                uncertainty=_SYNTHETIC,
            ),
        ),
        disease_direction=disease_direction,
        therapeutic_target_node_id=target_id,
        required_correction=corrective,
        uncertainties=(_SYNTHETIC,),
        developmental_window_caveat=_SYNTHETIC,
    )


def _canary_entry(
    *,
    drug_id: str,
    name: str,
    target_node_id: str,
    observed_direction: EffectDirection,
    worsens_cin: bool | None,
) -> CatalogEntry:
    """An agent with every non-direction question already answered.

    Everything except ``observed_direction`` and ``worsens_cin`` is held at the
    unambiguously favourable setting, so that a rejection can only come from the
    axis the canary is varying. A canary that is rejected for three reasons at once
    proves nothing about any of them.
    """
    return CatalogEntry(
        drug_id=drug_id,
        name=name,
        approved_name=name.lower(),
        approval_status=ApprovalStatus.APPROVED,
        intervention_class=InterventionClass.DISEASE_MODIFYING,
        target=f"Fabricated target of {name} (SYNTHETIC)",
        target_node_id=target_node_id,
        mechanism_of_action=_SYNTHETIC,
        observed_direction=observed_direction,
        is_direct_evidence=True,
        strongest_evidence_type=EvidenceType.HUMAN_TRIAL,
        has_pediatric_exposure=True,
        youngest_age_studied="fabricated",
        pediatric_indication="fabricated indication (SYNTHETIC)",
        route="oral",
        cns_penetrant=True,
        achievable_plasma_um=100.0,
        required_effective_um=1.0,
        half_life_hours=12.0,
        worsens_cin=worsens_cin,
        validation_experiment=_SYNTHETIC,
        notes=_SYNTHETIC,
    )


@dataclass(frozen=True, slots=True)
class PipelineCanary:
    """One fabricated pair and the verdict the pipeline must reach on it."""

    canary_id: str
    what_it_proves: str
    mechanism: MechanismHypothesis
    entry: CatalogEntry
    expect_accepted: bool
    expect_primary_reason: RejectionReason | None
    expect_directions_agree: bool | None
    #: Reasons that must be recorded against the candidate WITHOUT killing it.
    #: An accepted candidate that silently drops its concerns is the failure this
    #: field exists to catch: "we accepted it despite X" is the sentence a reviewer
    #: needs, and it disappears the moment a non-fatal reason stops being carried.
    expect_concerns: tuple[RejectionReason, ...] = ()


@dataclass(frozen=True, slots=True)
class CanaryOutcome:
    """What the pipeline actually did with one canary."""

    canary_id: str
    accepted: bool
    primary_reason: RejectionReason | None
    directions_agree: bool | None
    concerns: tuple[RejectionReason, ...]
    score: float


@dataclass(frozen=True, slots=True)
class CanaryFailure:
    """One canary whose verdict was not the one that was demanded."""

    canary_id: str
    what_it_proves: str
    detail: str


@dataclass(frozen=True, slots=True)
class TriageLiveness:
    """The result of re-proving the triage stage against its own controls."""

    outcomes: tuple[CanaryOutcome, ...]
    failures: tuple[CanaryFailure, ...]

    @property
    def ok(self) -> bool:
        return not self.failures

    @property
    def accepted_count(self) -> int:
        return sum(1 for outcome in self.outcomes if outcome.accepted)

    def summary(self) -> str:
        if self.ok:
            return (
                f"{len(self.outcomes)} pipeline canary/canaries reached the demanded "
                f"verdict ({self.accepted_count} accepted, "
                f"{len(self.outcomes) - self.accepted_count} rejected)."
            )
        return "\n".join(
            [f"{len(self.failures)} pipeline canary/canaries reached the wrong verdict:"]
            + [
                f"  [{failure.canary_id}] {failure.detail}\n      proves: {failure.what_it_proves}"
                for failure in self.failures
            ]
        )


def pipeline_canaries() -> tuple[PipelineCanary, ...]:
    """The controls. Fabricated throughout; see this module's docstring (GP-20)."""
    primary = canary_mechanism(
        mechanism_id="MECH-CANARY-01",
        gene_symbol="CANARYGENE1",
        disease_direction=EffectDirection.DECREASE,
    )
    sibling = canary_mechanism(
        mechanism_id="MECH-CANARY-02",
        gene_symbol="CANARYGENE2",
        disease_direction=EffectDirection.INCREASE,
    )
    target = primary.therapeutic_target_node_id

    def entry(
        drug_id: str,
        name: str,
        direction: EffectDirection,
        *,
        node: str = target,
        worsens_cin: bool | None = False,
    ) -> CatalogEntry:
        return _canary_entry(
            drug_id=drug_id,
            name=name,
            target_node_id=node,
            observed_direction=direction,
            worsens_cin=worsens_cin,
        )

    return (
        PipelineCanary(
            canary_id="ACCEPT",
            what_it_proves=(
                "The pipeline can still say yes. Every other Track 2 test asserts a "
                "rejection, and a stage that rejected everything would pass all of them."
            ),
            mechanism=primary,
            entry=entry("DRUG-CANARY-POS", "Canarystatin", EffectDirection.INCREASE),
            expect_accepted=True,
            expect_primary_reason=None,
            expect_directions_agree=True,
        ),
        PipelineCanary(
            canary_id="INVERT",
            what_it_proves=(
                "The sign check is alive. The same agent with only its direction flipped "
                "must be rejected FOR THAT REASON, not for some other reason that happens "
                "to fire."
            ),
            mechanism=primary,
            entry=entry("DRUG-CANARY-NEG", "Canarystatin (inverted)", EffectDirection.DECREASE),
            expect_accepted=False,
            expect_primary_reason=RejectionReason.WRONG_DIRECTION,
            expect_directions_agree=False,
        ),
        PipelineCanary(
            canary_id="UNKNOWN_DIRECTION",
            what_it_proves=(
                "GP-16 / ASSUMPTION-DRUG-02: unknown is tri-state, and the two unknowns "
                "are not alike. Wrong direction kills; unsigned direction is a flag and a "
                "penalty. Accepted, but accepted WITH the flag still attached."
            ),
            mechanism=primary,
            entry=entry("DRUG-CANARY-UNK", "Canaryase", EffectDirection.UNKNOWN),
            expect_accepted=True,
            expect_primary_reason=None,
            expect_directions_agree=None,
            expect_concerns=(RejectionReason.DIRECTION_UNKNOWN,),
        ),
        PipelineCanary(
            canary_id="UNASSESSED_RISK",
            what_it_proves=(
                "GP-14 and ADR 0011: the accepted agent with one cell blanked is rejected "
                "on the UNASKED question while its direction still agrees. The rejection "
                "must not be smuggled in as a direction problem."
            ),
            mechanism=primary,
            entry=entry(
                "DRUG-CANARY-UNASSESSED",
                "Canarystatin (risk unassessed)",
                EffectDirection.INCREASE,
                worsens_cin=None,
            ),
            expect_accepted=False,
            expect_primary_reason=RejectionReason.ONCOGENIC_RISK,
            expect_directions_agree=True,
        ),
        PipelineCanary(
            canary_id="OFF_CHAIN",
            what_it_proves=(
                "The triage stage reads the mechanism it was handed. An agent acting on a "
                "node that is not on the chain is a finding, not a non-event (GP-13)."
            ),
            mechanism=primary,
            entry=entry(
                "DRUG-CANARY-OFFCHAIN",
                "Canarystatin (off-chain target)",
                EffectDirection.INCREASE,
                node="NOT_A_NODE_ON_THIS_CHAIN",
            ),
            expect_accepted=False,
            expect_primary_reason=RejectionReason.TARGET_NOT_IN_MECHANISM,
            # None, not False: a node that is not on the chain has no required
            # correction, so there is nothing for the observed direction to disagree
            # WITH. Inventing a sign here would be the mistake, not catching one.
            expect_directions_agree=None,
        ),
        PipelineCanary(
            canary_id="SIBLING_GENE",
            what_it_proves=(
                "The verdict is a property of the mechanism-plus-drug pair, not of the "
                "drug. The agent accepted against CANARYGENE1 is rejected unchanged "
                "against CANARYGENE2, whose disease direction is the mirror image."
            ),
            mechanism=sibling,
            entry=entry(
                "DRUG-CANARY-POS",
                "Canarystatin",
                EffectDirection.INCREASE,
                node=sibling.therapeutic_target_node_id,
            ),
            expect_accepted=False,
            expect_primary_reason=RejectionReason.WRONG_DIRECTION,
            expect_directions_agree=False,
        ),
    )


def run_canary(canary: PipelineCanary) -> CanaryOutcome:
    """Push one canary through the real triage function.

    Deliberately `generate_drug_hypotheses` and not a reimplementation of the
    direction check, for the same reason `check_detector_liveness` calls `scan_bytes`
    rather than `rule.pattern.search`: the gates, the fatal/non-fatal split and the
    reason-priority ordering all live between the two, and proving the inner
    comparison in isolation would leave that whole path unexercised.
    """
    result = generate_drug_hypotheses(
        mechanism=canary.mechanism,
        catalog=DrugCatalog([canary.entry], version=CANARY_VERSION),
        clock=demo_clock(),
    )
    (hypothesis,) = result.all_hypotheses
    return CanaryOutcome(
        canary_id=canary.canary_id,
        accepted=not hypothesis.rejected,
        primary_reason=(hypothesis.rejection_reasons[0] if hypothesis.rejection_reasons else None),
        directions_agree=hypothesis.directions_agree,
        concerns=tuple(hypothesis.concerns),
        score=hypothesis.score,
    )


def check_triage_liveness() -> TriageLiveness:
    """Re-prove the drug-triage stage against its own controls, live.

    A live probe rather than a static assertion, arriving at the science from the
    same direction `mva.privacy.audit.check_detector_liveness` arrives at the
    scanner. Every static case for "the triage works" is an argument about source
    code, and a triage stage that has stopped discriminating produces exactly the
    output a careful one does when everything it is shown happens to be bad.
    """
    outcomes: list[CanaryOutcome] = []
    failures: list[CanaryFailure] = []
    for canary in pipeline_canaries():
        outcome = run_canary(canary)
        outcomes.append(outcome)
        problems: list[str] = []
        if outcome.accepted is not canary.expect_accepted:
            problems.append(
                f"expected {'ACCEPTED' if canary.expect_accepted else 'REJECTED'}, "
                f"got {'ACCEPTED' if outcome.accepted else 'REJECTED'}"
            )
        if outcome.primary_reason != canary.expect_primary_reason:
            problems.append(
                f"expected primary reason "
                f"{getattr(canary.expect_primary_reason, 'value', None)!r}, got "
                f"{getattr(outcome.primary_reason, 'value', None)!r}"
            )
        if outcome.directions_agree is not canary.expect_directions_agree:
            problems.append(
                f"expected directions_agree {canary.expect_directions_agree!r}, "
                f"got {outcome.directions_agree!r}"
            )
        missing = [r for r in canary.expect_concerns if r not in outcome.concerns]
        if missing:
            problems.append(
                f"expected concern(s) {[r.value for r in missing]} to be carried on the "
                f"candidate, got {[r.value for r in outcome.concerns]}"
            )
        if problems:
            failures.append(
                CanaryFailure(
                    canary_id=canary.canary_id,
                    what_it_proves=canary.what_it_proves,
                    detail="; ".join(problems),
                )
            )
    return TriageLiveness(outcomes=tuple(outcomes), failures=tuple(failures))
