"""Safety and evidence-quality assessment for a candidate agent.

Two functions, two different jobs:

* `assess_safety` asks *could this hurt this particular patient?* It is
  context-sensitive on purpose. An agent that increases chromosomal instability is
  a routine oncology tool; in a germline chromosomal-instability disorder it is
  the disease mechanism in a bottle, and it is disqualified outright
  (`ONCOGENIC_RISK`).
* `assess_evidence_quality` asks *how far is this claim from a human outcome?* It
  ranks evidence tiers honestly, so "a compound binds the target in a plate" can
  never be scored like "a randomised trial in children showed benefit".

The rule that governs both: **an unknown is a concern, never a pass.** Missing
paediatric exposure, unassessed oncogenic risk and unmeasured exposure are all
recorded explicitly, and none of them is silently treated as fine — the commonest
way a gap disappears from a write-up is that nobody had a field to put it in.

Which unknowns are *fatal* is context-dependent, and deliberately so. An
unassessed oncogenic risk is a recorded concern in general and a **disqualifying**
one in a chromosomal-instability disorder, because there the unmeasured quantity
is the disease mechanism itself.

Nothing here is medical advice; every output is a research hypothesis.
"""

from __future__ import annotations

from dataclasses import dataclass

from mva.clock import Clock
from mva.interventions.catalog import CatalogEntry
from mva.models.drug import RejectionReason, SafetyConcern
from mva.models.evidence import EvidenceType
from mva.models.mechanism import MechanismHypothesis

__all__ = [
    "EVIDENCE_TYPE_WEIGHT",
    "INDIRECT_EVIDENCE_MULTIPLIER",
    "SEVERITY_PENALTY",
    "SafetyVerdict",
    "assess_evidence_quality",
    "assess_safety",
    "is_chromosomal_instability_context",
    "is_neurological_context",
]

#: Lower-cased substrings that mark a mechanism as a chromosomal-instability
#: context. A keyword scan over curated labels is a deliberately crude proxy; it
#: is documented as such and it errs toward flagging (GP-20). A negative result
#: means "not detected", never "this agent is safe here".
CIN_MARKERS: tuple[str, ...] = (
    "aneuploid",
    "missegregation",
    "chromosomal instability",
    "chromosome instability",
    "genome instability",
    "genomic instability",
    "micronucle",
    "spindle assembly checkpoint",
    "mitotic checkpoint",
)

#: Substrings that mark a phenotype as neurological, for CNS-penetration checks.
NEURO_MARKERS: tuple[str, ...] = (
    "microcephaly",
    "brain",
    "cerebral",
    "cerebellar",
    "neuro",
    "cognitive",
    "seizure",
    "epilep",
    "developmental delay",
    "intellectual",
)

#: Evidence tiers ranked by distance from a demonstrated human outcome.
#: human trial > human observational > case report > animal > primary patient
#: cells > cell line > biochemical binding > in silico. GP-32: heuristic weights.
EVIDENCE_TYPE_WEIGHT: dict[EvidenceType, float] = {
    EvidenceType.HUMAN_TRIAL: 1.00,
    EvidenceType.HUMAN_OBSERVATIONAL: 0.85,
    EvidenceType.HUMAN_CASE_REPORT: 0.70,
    EvidenceType.ANIMAL_MODEL: 0.60,
    EvidenceType.DIRECT_MEASUREMENT: 0.55,
    EvidenceType.PRIMARY_PATIENT_CELLS: 0.50,
    EvidenceType.CELL_LINE: 0.35,
    EvidenceType.EXPERT_REVIEW: 0.30,
    EvidenceType.CURATED_DATABASE: 0.30,
    EvidenceType.BIOCHEMICAL_BINDING: 0.25,
    EvidenceType.IN_SILICO_PREDICTION: 0.10,
    EvidenceType.PIPELINE_INFERENCE: 0.05,
}

#: Applied when the evidence is for a related agent, target or system rather than
#: for THIS agent on THIS target. Indirect evidence is real evidence, worth
#: roughly half of the direct article — not zero, and not the same.
INDIRECT_EVIDENCE_MULTIPLIER = 0.60

#: How much each recorded concern costs the safety term of a drug's score.
SEVERITY_PENALTY: dict[str, float] = {
    "critical": 1.00,
    "major": 0.35,
    "moderate": 0.20,
    "minor": 0.10,
}


@dataclass(frozen=True)
class SafetyVerdict:
    """Everything the safety pass concluded about one agent.

    `reasons` lists **every** reason raised, fatal or not; `disqualifying` says
    whether any of them is fatal. Keeping the two apart is what lets a candidate be
    accepted while still carrying, visibly, the concerns that were found — the
    alternative is a binary verdict that quietly discards the near misses.
    """

    concerns: tuple[SafetyConcern, ...]
    disqualifying: bool
    reasons: tuple[RejectionReason, ...]
    rationale: str


def _matches(mechanism: MechanismHypothesis, markers: tuple[str, ...]) -> bool:
    haystack = " ".join(
        [mechanism.summary, *(f"{node.label} {node.description}" for node in mechanism.nodes)]
    ).lower()
    return any(marker in haystack for marker in markers)


def is_chromosomal_instability_context(mechanism: MechanismHypothesis) -> bool:
    """Does this mechanism describe a chromosomal-instability disease?

    Keyword scan over the curated summary and node labels. Crude by design: the
    cost of a false positive is an extra flag on a report, the cost of a false
    negative is recommending an aneuploidy-inducing agent to a patient whose cells
    already sit at their aneuploidy-tolerance ceiling.
    """
    return _matches(mechanism, CIN_MARKERS)


def is_neurological_context(mechanism: MechanismHypothesis) -> bool:
    """Does the phenotype involve the CNS, making brain penetration decisive?"""
    return _matches(mechanism, NEURO_MARKERS)


def assess_safety(
    entry: CatalogEntry, *, mechanism: MechanismHypothesis, clock: Clock
) -> SafetyVerdict:
    """Assess one agent against this specific disease context.

    Concerns are emitted in a fixed order (oncogenic risk, paediatric exposure,
    exposure/PK, route, CNS penetration) so repeat runs are byte-identical (GP-30).
    The returned `SafetyConcern` values carry no evidence IDs; the triage stage
    emits one `EvidenceItem` per concern and links them back.
    """
    concerns: list[SafetyConcern] = []
    reasons: list[RejectionReason] = []
    notes: list[str] = []
    cin_context = is_chromosomal_instability_context(mechanism)

    # 1. Oncogenic / chromosomal-instability risk. Mandatory question 7.
    if entry.worsens_cin is True:
        concerns.append(
            SafetyConcern(
                concern_id=f"SC-{entry.drug_id}-CIN",
                description=(
                    f"{entry.name} is recorded as increasing chromosomal instability. In a "
                    "germline chromosomal-instability disorder this pushes the disease "
                    "mechanism itself, in non-tumour cells that have no tolerance reserve."
                ),
                severity="critical",
                population=(
                    "Children with a constitutional chromosomal-instability disorder"
                    if cin_context
                    else "Patients with elevated baseline aneuploidy"
                ),
                is_disqualifying=cin_context,
            )
        )
        reasons.append(RejectionReason.ONCOGENIC_RISK)
        notes.append(
            "Increases chromosomal instability; disqualifying in this disease context."
            if cin_context
            else "Increases chromosomal instability; no CIN context detected, flagged not fatal."
        )
    elif entry.worsens_cin is None:
        # An unassessed answer to mandatory question 7 is BLOCKING in a
        # chromosomal-instability context, and the whole pipeline says so: the model
        # field, the evidence schema, the report text and ASSUMPTION-DRUG-07 all call
        # it a blocking gap. It therefore has to block. Recording it as a non-fatal
        # concern let an agent whose oncogenic risk nobody had measured be presented
        # as a ranked candidate in a report that told the reader, in the same
        # paragraph, that the question "must be answered before any further
        # consideration". Outside a CIN context it stays a recorded concern.
        concerns.append(
            SafetyConcern(
                concern_id=f"SC-{entry.drug_id}-CIN-UNASSESSED",
                description=(
                    f"Whether {entry.name} increases aneuploidy or cancer susceptibility has "
                    "not been assessed. Unassessed is not negative: in a germline "
                    "chromosomal-instability disorder it is a blocking gap that must be "
                    "closed before any use."
                    if cin_context
                    else f"Whether {entry.name} increases aneuploidy or cancer susceptibility "
                    "has not been assessed. No chromosomal-instability context was detected "
                    "for this mechanism, so the gap is recorded rather than fatal."
                ),
                severity="critical" if cin_context else "major",
                population=(
                    "Children with a constitutional chromosomal-instability disorder"
                    if cin_context
                    else "Not assessed in any population"
                ),
                is_disqualifying=cin_context,
            )
        )
        reasons.append(
            RejectionReason.ONCOGENIC_RISK if cin_context else RejectionReason.SAFETY_CONCERN
        )
        notes.append(
            "Oncogenic/CIN risk unassessed; blocking in this disease context."
            if cin_context
            else "Oncogenic/CIN risk unassessed; no CIN context detected, flagged not fatal."
        )

    # 2. Paediatric exposure. Recorded always; never fatal on its own.
    if not entry.has_pediatric_exposure:
        concerns.append(
            SafetyConcern(
                concern_id=f"SC-{entry.drug_id}-PAED",
                description=(
                    f"No paediatric exposure is recorded for {entry.name}. Dosing, "
                    "tolerability and developmental toxicity in children are therefore "
                    "unknown rather than reassuring."
                ),
                severity="major",
                population="Children",
                is_disqualifying=False,
            )
        )
        reasons.append(RejectionReason.NO_PEDIATRIC_EVIDENCE)
        notes.append("No paediatric exposure on record.")

    # 3. Exposure. Mandatory question 5: is the effective concentration reachable?
    achievable = entry.concentration_achievable
    if achievable is False:
        concerns.append(
            SafetyConcern(
                concern_id=f"SC-{entry.drug_id}-CONC",
                description=(
                    f"{entry.name} requires {entry.required_effective_um} uM to act but peaks "
                    f"at {entry.achievable_plasma_um} uM in plasma. The mechanism cannot be "
                    "reached at a tolerated dose, however correct it is."
                ),
                severity="critical",
                population="Any",
                is_disqualifying=True,
            )
        )
        reasons.append(RejectionReason.CONCENTRATION_NOT_ACHIEVABLE)
        notes.append("Effective concentration is not clinically achievable.")
    elif achievable is None:
        concerns.append(
            SafetyConcern(
                concern_id=f"SC-{entry.drug_id}-CONC-UNKNOWN",
                description=(
                    f"Achievable and/or required concentration for {entry.name} is not "
                    "recorded, so it is unknown whether the mechanism can be reached in vivo."
                ),
                severity="major",
                population="Any",
                is_disqualifying=False,
            )
        )
        reasons.append(RejectionReason.PHARMACOKINETIC_BARRIER)
        notes.append("Exposure relative to the effective concentration is unknown.")

    # 4. Route of administration.
    if not entry.has_administration_route:
        concerns.append(
            SafetyConcern(
                concern_id=f"SC-{entry.drug_id}-ROUTE",
                description=(
                    f"No viable route of administration is recorded for {entry.name}; it has "
                    "no in vivo formulation."
                ),
                severity="major",
                population="Any",
                is_disqualifying=False,
            )
        )
        reasons.append(RejectionReason.PHARMACOKINETIC_BARRIER)
        notes.append("No viable administration route.")

    # 5. CNS penetration, only where the phenotype makes it decisive.
    if is_neurological_context(mechanism) and entry.cns_penetrant is not True:
        unknown = entry.cns_penetrant is None
        concerns.append(
            SafetyConcern(
                concern_id=f"SC-{entry.drug_id}-CNS",
                description=(
                    f"The phenotype is neurological and {entry.name} is "
                    + ("of unknown CNS penetration" if unknown else "not CNS-penetrant")
                    + ". Target engagement in the affected tissue is "
                    + ("unestablished." if unknown else "unlikely.")
                ),
                severity="moderate" if unknown else "major",
                population="Patients with a CNS phenotype",
                is_disqualifying=False,
            )
        )
        reasons.append(RejectionReason.PHARMACOKINETIC_BARRIER)
        notes.append("CNS penetration is inadequate or unknown for a neurological phenotype.")

    disqualifying = any(concern.is_disqualifying for concern in concerns)
    if not notes:
        notes.append("No safety concern was raised by the recorded fields.")
    rationale = (
        f"{len(concerns)} concern(s) recorded for {entry.drug_id} against mechanism "
        f"{mechanism.mechanism_id} as of {clock.now().date().isoformat()}: "
        + " ".join(notes)
        + " Absence of a concern here reflects the fields the catalogue records, not a "
        "safety clearance."
    )
    return SafetyVerdict(
        concerns=tuple(concerns),
        disqualifying=disqualifying,
        reasons=tuple(dict.fromkeys(reasons)),
        rationale=rationale,
    )


def assess_evidence_quality(entry: CatalogEntry) -> tuple[float, str]:
    """Score the evidence behind an agent in ``[0, 1]``, with a rationale.

    Formula (GP-32: the weights are documented heuristics, not measurements)::

        score = EVIDENCE_TYPE_WEIGHT[strongest_evidence_type]
              * (1.0 if is_direct_evidence else INDIRECT_EVIDENCE_MULTIPLIER)

    The tier ladder runs human trial (1.00) > human observational (0.85) > human
    case report (0.70) > animal model (0.60) > primary patient cells (0.50) >
    cell line (0.35) > biochemical binding (0.25) > in silico (0.10). The ordering
    is the point: it is what keeps a binding assay from being read as a result.
    """
    weight = EVIDENCE_TYPE_WEIGHT.get(entry.strongest_evidence_type, 0.05)
    multiplier = 1.0 if entry.is_direct_evidence else INDIRECT_EVIDENCE_MULTIPLIER
    score = round(max(0.0, min(1.0, weight * multiplier)), 6)
    directness = (
        "direct evidence for this agent on this target"
        if entry.is_direct_evidence
        else (
            "INDIRECT evidence only (a related agent, target or system), discounted by "
            f"{INDIRECT_EVIDENCE_MULTIPLIER:.2f}"
        )
    )
    rationale = (
        f"Strongest evidence for {entry.drug_id} is "
        f"{entry.strongest_evidence_type.value} (tier weight {weight:.2f}), and it is "
        f"{directness}; evidence quality scores {score:.3f}. This grades the distance from "
        "a demonstrated human outcome, not whether the agent would work."
    )
    return score, rationale
