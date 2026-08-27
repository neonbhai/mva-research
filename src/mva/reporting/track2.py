"""Track 2 reporting: mechanism, drug hypotheses, and the rejections.

Everything here exists to make one failure impossible to commit quietly.

In a loss-of-function checkpoint disorder, the compounds a naive target-proximity
search ranks first are checkpoint **inhibitors** — they bind exactly the proteins
the mechanism names. They were developed to push chromosomally unstable tumour
cells past their aneuploidy-tolerance ceiling, and the patient's non-tumour cells
already sit above that ceiling with no reserve. The proximity is high and the sign
is inverted. So the reports are built so that direction is never a footnote: a
wrong-direction agent cannot be constructed as accepted at all (the model
validator refuses it), and it still appears here, by name, in the rejection
record with its reason (GP-19).

Three separations are structural, not stylistic:

* **"Cannot determine" is its own section.** A drug whose ``directions_agree`` is
  ``None`` is never listed beside one where it is ``True``. An unknown that gets
  filed under agreement is an unknown that has been laundered into support
  (GP-16, ASSUMPTION-DRUG-02).
* **Symptomatic agents have their own section.** Presenting an anticonvulsant as
  addressing a chromosome-segregation defect is a category error, and the report
  cannot commit it because the sections are built from
  ``InterventionClass`` (ASSUMPTION-DRUG-04).
* **Rejections are printed, not dropped.** The rejection record is embedded in
  every drug report. A pipeline that silently discards candidates cannot be
  audited, and the discarded set is where the contraindicated compounds are.

**GP-10 in this module.** A headline claim must be sourced or the report is
refused: the mechanism summary, and the rationale for any *accepted* drug. An
individual chain link that cites nothing is a different case — it is rendered as
a gap rather than as a claim, and named in the "inferred links" section, because
an uncited link is a real state of the record and suppressing the whole report
would suppress the gap along with it.

Nothing in this module is medical advice.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

from mva.clock import Clock
from mva.errors import ReportCompletenessError
from mva.models import (
    AssertionTier,
    CandidatePair,
    DrugHypothesis,
    EvidenceItem,
    InterventionClass,
    MechanismHypothesis,
    MechanismLink,
    MechanismNode,
    RejectionReason,
)
from mva.reporting.assertions import TIER_MARKERS, Assertion, AssertionChecker, weakest_tier
from mva.reporting.render import NOT_RECORDED, format_cell, markdown_table, render_template

if TYPE_CHECKING:  # pragma: no cover - typing only
    from mva.evidence.ledger import AssertionResolver

#: Printed at the top of every Track 2 artifact. Not a footer, not a footnote.
NOT_MEDICAL_ADVICE = (
    "**RESEARCH HYPOTHESIS — NOT MEDICAL ADVICE.** Everything below is a "
    "computationally generated hypothesis requiring pre-clinical validation. It is "
    "not a treatment recommendation, not a prescription, and not a clinical "
    "opinion. No compound named here should be administered to any patient on the "
    "basis of this document."
)

#: Stated wherever a number appears, because a score that looks calibrated will be
#: read as calibrated.
UNCALIBRATED_WEIGHTS = (
    "**Scoring weights are uncalibrated heuristics.** They were chosen by reasoning "
    "about this disease context and are fitted to no labelled dataset "
    "(ASSUMPTION-SCORING-01). No score in this document is a probability, and no "
    "ranking here carries a claim of clinical validity."
)

#: The eight mandatory questions every drug must answer, in order.
DRUG_QUESTIONS: tuple[str, ...] = (
    "What node is modified?",
    "In what direction must it move?",
    "Does the drug act that way?",
    "What evidence tier supports it (binding / cells / animals / humans)?",
    "Was the concentration clinically achievable?",
    "Is there paediatric exposure evidence?",
    "Could it worsen chromosome instability or cancer susceptibility?",
    "What experiment should be performed before any clinical consideration?",
)

#: What would have to change for a rejection to be revisited. Keyed by reason so
#: the rejection record is actionable rather than merely final.
REVERSAL_CONDITIONS: dict[RejectionReason, str] = {
    RejectionReason.WRONG_DIRECTION: (
        "Nothing short of the mechanism itself being wrong. If the required correction "
        "is right, this agent pushes the target the wrong way and is contraindicated."
    ),
    RejectionReason.DIRECTION_UNKNOWN: (
        "A signed direction of effect on this target, measured in a relevant system."
    ),
    RejectionReason.NOT_APPROVED: (
        "Regulatory approval. Until then it is a research probe, not a repurposing "
        "candidate (ASSUMPTION-DRUG-03)."
    ),
    RejectionReason.INSUFFICIENT_EVIDENCE: (
        "Direct experimental evidence for this agent on this target in a relevant system."
    ),
    RejectionReason.SAFETY_CONCERN: (
        "Evidence that the concern does not apply to a germline chromosomal-instability "
        "population, or a monitorable mitigation."
    ),
    RejectionReason.ONCOGENIC_RISK: (
        "Evidence that the agent does not increase aneuploidy or cancer susceptibility "
        "in non-tumour cells."
    ),
    RejectionReason.NO_PEDIATRIC_EVIDENCE: (
        "Paediatric exposure data — noting that tolerability in a paediatric oncology "
        "population does not transfer to this population (ASSUMPTION-DRUG-06)."
    ),
    RejectionReason.PHARMACOKINETIC_BARRIER: (
        "A formulation or route that reaches the affected tissue."
    ),
    RejectionReason.TARGET_NOT_IN_MECHANISM: (
        "A mechanism revision that places this target in the chain, with evidence."
    ),
    RejectionReason.CONCENTRATION_NOT_ACHIEVABLE: (
        "Measured plasma or tissue exposure at or above the effective concentration "
        "(ASSUMPTION-DRUG-05)."
    ),
    RejectionReason.MECHANISM_MISMATCH: (
        "Evidence linking this target to the mechanism actually operating in this gene "
        "(ASSUMPTION-MECHANISM-03)."
    ),
}

_LINK_HEADERS = ("Link", "Relation", "Direction", "Tier", "Strength", "Demonstrated?", "Evidence")
_NODE_HEADERS = ("Node", "Kind", "State in patient", "Identifier")


# ---------------------------------------------------------------------------
# Public builders
# ---------------------------------------------------------------------------


def build_mechanism_report(
    mechanism: MechanismHypothesis,
    *,
    resolver: AssertionResolver,
    clock: Clock,
) -> str:
    """Render the mechanism chain, link by link, with its weak points named."""
    checker = AssertionChecker(resolver, strict=True)
    return render_template(
        "mechanism_report.md.j2", _mechanism_context(mechanism, resolver, checker, clock)
    )


def build_drug_report(
    accepted: Sequence[DrugHypothesis],
    rejected: Sequence[DrugHypothesis],
    *,
    mechanism: MechanismHypothesis,
    resolver: AssertionResolver,
    clock: Clock,
) -> str:
    """Render the drug hypotheses alone, with the rejection record embedded."""
    return _render_track2(
        mechanism=mechanism,
        accepted=accepted,
        rejected=rejected,
        pair=None,
        resolver=resolver,
        clock=clock,
        mechanism_block="",
        title=f"Drug repurposing hypotheses — {mechanism.gene_symbol}",
    )


def build_rejection_record(rejected: Sequence[DrugHypothesis], *, clock: Clock) -> str:
    """Render every rejected candidate with its reasons (GP-19).

    Takes no resolver: a rejection record states no scientific claim about
    efficacy. It records a decision, the reason for it, and what would have to
    change — which is exactly what makes the pipeline auditable rather than
    merely confident.
    """
    context: dict[str, object] = {
        "generated_at": clock.now().isoformat(),
        # The record is also written as a standalone file (Snakefile: rejection_record.md),
        # where it names compounds and lists what would have to change for each. Read in
        # isolation, with no banner, that reads as a shortlist of things to try next.
        "not_medical_advice": NOT_MEDICAL_ADVICE,
        "rejected": [_rejection_context(drug) for drug in sorted(rejected, key=_drug_sort_key)],
        "count": len(rejected),
    }
    return render_template("rejection_record.md.j2", context)


def build_track2_report(
    mechanism: MechanismHypothesis,
    accepted: Sequence[DrugHypothesis],
    rejected: Sequence[DrugHypothesis],
    *,
    pair: CandidatePair | None,
    resolver: AssertionResolver,
    clock: Clock,
) -> str:
    """The full Track 2 deliverable: variant pair, mechanism, drugs, rejections."""
    return _render_track2(
        mechanism=mechanism,
        accepted=accepted,
        rejected=rejected,
        pair=pair,
        resolver=resolver,
        clock=clock,
        mechanism_block=build_mechanism_report(mechanism, resolver=resolver, clock=clock),
        title=f"Track 2 — mechanism-grounded repurposing hypotheses for {mechanism.gene_symbol}",
    )


# ---------------------------------------------------------------------------
# Context assembly
# ---------------------------------------------------------------------------


def _render_track2(
    *,
    mechanism: MechanismHypothesis,
    accepted: Sequence[DrugHypothesis],
    rejected: Sequence[DrugHypothesis],
    pair: CandidatePair | None,
    resolver: AssertionResolver,
    clock: Clock,
    mechanism_block: str,
    title: str,
) -> str:
    checker = AssertionChecker(resolver, strict=True)
    presentable, all_rejected = _partition_rejected(accepted, rejected)

    agreeing = _select(presentable, InterventionClass.DISEASE_MODIFYING, agree=True)
    undetermined = _select(presentable, InterventionClass.DISEASE_MODIFYING, agree=None)
    other_class = _select_non_disease_modifying(presentable)
    _check_every_hypothesis_is_rendered(
        accepted=accepted,
        rejected=rejected,
        sections=(agreeing, undetermined, other_class, all_rejected),
    )

    context: dict[str, object] = {
        "title": title,
        "generated_at": clock.now().isoformat(),
        "not_medical_advice": NOT_MEDICAL_ADVICE,
        "uncalibrated_weights": UNCALIBRATED_WEIGHTS,
        "mechanism_block": mechanism_block,
        "mechanism": _mechanism_glance(mechanism, resolver, checker),
        "pair": _pair_context(pair),
        "disease_modifying": [
            _drug_context(drug, mechanism, resolver, checker, require_rationale=True)
            for drug in agreeing
        ],
        "direction_undetermined": [
            _drug_context(drug, mechanism, resolver, checker, require_rationale=True)
            for drug in undetermined
        ],
        "symptomatic": [
            _drug_context(drug, mechanism, resolver, checker, require_rationale=False)
            for drug in other_class
        ],
        "rejection_record": build_rejection_record(all_rejected, clock=clock),
        "rejected_count": len(all_rejected),
        "caveats": _caveats(mechanism),
    }
    return render_template("track2_report.md.j2", context)


def _check_every_hypothesis_is_rendered(
    *,
    accepted: Sequence[DrugHypothesis],
    rejected: Sequence[DrugHypothesis],
    sections: tuple[tuple[DrugHypothesis, ...], ...],
) -> None:
    """Refuse to render unless every input hypothesis lands in exactly one section.

    The sections are chosen by predicate — ``intervention_class``, then
    ``directions_agree`` matched against ``is True`` / ``is None``, then
    ``rejected``. Those predicates are not a partition, and the case they miss is
    the worst one available: a disease-modifying agent whose ``directions_agree``
    is ``False`` while ``rejected`` is ``False`` matches neither AGREES nor CANNOT
    BE DETERMINED, and `_partition_rejected` only rescues ``rejected=True``. Such
    an object cannot be constructed — the `DrugHypothesis` validator forbids it —
    but it can be *copied* into existence through ``model_copy(update=...)``, which
    does not re-run validators.

    The observed consequence was a report reading "Candidates presented: 0 ...
    Rejected: 0" for a run that had triaged a contraindicated compound: it was
    neither presented nor rejected, it was gone. Erasing a contraindicated
    compound is strictly worse than listing it in the wrong section, so a gap here
    raises instead of shortening the report.
    """
    counts: dict[str, int] = {}
    for section in sections:
        for drug in section:
            counts[drug.drug_id] = counts.get(drug.drug_id, 0) + 1
    expected = {drug.drug_id for drug in (*accepted, *rejected)}
    missing = sorted(expected - set(counts))
    duplicated = sorted(drug_id for drug_id, count in counts.items() if count > 1)
    unexpected = sorted(set(counts) - expected)
    if not (missing or duplicated or unexpected):
        return
    msg = (
        f"Track 2 rendering would not account for every hypothesis it was given: "
        f"{len(missing)} would appear in no section {missing}, {len(duplicated)} in more "
        f"than one {duplicated}, {len(unexpected)} in a section but not in the input "
        f"{unexpected}. A hypothesis in no section is deleted from the report rather than "
        "reported, so this is refused. The usual cause is a DrugHypothesis whose "
        "directions_agree is False while rejected is False, which model_copy can produce "
        "and the constructor cannot: build it with revalidated_copy instead."
    )
    raise ReportCompletenessError(msg)


def _drug_sort_key(drug: DrugHypothesis) -> tuple[int, float, str]:
    return drug.sort_key()


def _partition_rejected(
    accepted: Sequence[DrugHypothesis], rejected: Sequence[DrugHypothesis]
) -> tuple[tuple[DrugHypothesis, ...], tuple[DrugHypothesis, ...]]:
    """Move anything flagged rejected out of the accepted list.

    Defensive on purpose: a drug carrying ``rejected=True`` that reaches the
    accepted list is a caller bug, and the safe direction to resolve it is
    towards the rejection record, never towards a recommendation.
    """
    misfiled = tuple(drug for drug in accepted if drug.rejected)
    presentable = tuple(drug for drug in accepted if not drug.rejected)
    combined = {drug.drug_id: drug for drug in (*rejected, *misfiled)}
    return presentable, tuple(sorted(combined.values(), key=_drug_sort_key))


def _select(
    drugs: Sequence[DrugHypothesis], intervention: InterventionClass, *, agree: bool | None
) -> tuple[DrugHypothesis, ...]:
    return tuple(
        sorted(
            (
                drug
                for drug in drugs
                if drug.intervention_class is intervention and drug.directions_agree is agree
            ),
            key=_drug_sort_key,
        )
    )


def _select_non_disease_modifying(drugs: Sequence[DrugHypothesis]) -> tuple[DrugHypothesis, ...]:
    return tuple(
        sorted(
            (
                drug
                for drug in drugs
                if drug.intervention_class is not InterventionClass.DISEASE_MODIFYING
            ),
            key=_drug_sort_key,
        )
    )


def _mechanism_glance(
    mechanism: MechanismHypothesis,
    resolver: AssertionResolver,
    checker: AssertionChecker,
) -> dict[str, object]:
    """The parts of the mechanism that must appear in *every* Track 2 artifact."""
    summary = _summary_assertion(mechanism, resolver, checker)
    target = mechanism.target_node()
    return {
        "mechanism_id": mechanism.mechanism_id,
        "gene_symbol": mechanism.gene_symbol,
        "summary": summary,
        "target_node_id": target.node_id,
        "target_label": target.label,
        "target_state": target.state_in_patient.value,
        "disease_direction": mechanism.disease_direction.value,
        "required_correction": mechanism.required_correction.value,
        "developmental_window_caveat": (
            mechanism.developmental_window_caveat
            or (
                "No developmental-window caveat was recorded. This is a gap, not a "
                "clearance: aneuploidy and its structural consequences are established "
                "during early development, and post-natal modulation cannot reverse "
                "malformation already present at birth (ASSUMPTION-MECHANISM-02)."
            )
        ),
        "inferred_links": [_link_summary(link) for link in mechanism.inferred_links],
        "unsourced_links": [
            _link_summary(link) for link in mechanism.links if not link.evidence_ids
        ],
        "is_fully_demonstrated": mechanism.is_fully_demonstrated,
        "uncertainties": list(mechanism.uncertainties),
    }


def _summary_assertion(
    mechanism: MechanismHypothesis,
    resolver: AssertionResolver,
    checker: AssertionChecker,
) -> str:
    """The mechanism's headline claim, gated by GP-10.

    An unsourced mechanism summary is refused outright. This is the sentence the
    entire Track 2 report hangs from; if it cites nothing, there is nothing below
    it worth rendering.
    """
    items = resolver.resolve(mechanism.supporting_evidence_ids)
    return checker.check(
        Assertion(
            text=mechanism.summary,
            # SPECULATION is the placeholder for the no-evidence case; the checker
            # raises before it can ever be rendered.
            tier=weakest_tier(items) if items else AssertionTier.SPECULATION,
            evidence_ids=tuple(sorted(mechanism.supporting_evidence_ids)),
        )
    ).rendered()


def _mechanism_context(
    mechanism: MechanismHypothesis,
    resolver: AssertionResolver,
    checker: AssertionChecker,
    clock: Clock,
) -> dict[str, object]:
    nodes = sorted(mechanism.nodes, key=lambda node: node.node_id)
    links = sorted(mechanism.links, key=lambda link: link.link_id)
    return {
        "generated_at": clock.now().isoformat(),
        "not_medical_advice": NOT_MEDICAL_ADVICE,
        "mechanism": _mechanism_glance(mechanism, resolver, checker),
        "pair_id": mechanism.pair_id or NOT_RECORDED,
        "nodes_table": markdown_table(
            _NODE_HEADERS,
            [
                [node.label, node.kind.value, node.state_in_patient.value, node.identifier]
                for node in nodes
            ],
        ),
        "links_table": markdown_table(
            _LINK_HEADERS,
            [
                [
                    f"{link.source_node_id} -> {link.target_node_id}",
                    link.relation,
                    link.direction.value,
                    f"{link.tier.value} {TIER_MARKERS[link.tier]}",
                    link.strength.value,
                    link.is_directly_demonstrated,
                    ", ".join(sorted(link.evidence_ids)) or "UNSOURCED",
                ]
                for link in links
            ],
        ),
        "links": [_link_context(link, checker) for link in links],
        "contradictions": [
            _evidence_context(item)
            for item in resolver.resolve(mechanism.contradicting_evidence_ids)
        ],
        "experiments": [
            {
                "experiment_id": experiment.experiment_id,
                "description": experiment.description,
                "measures": experiment.measures,
                "distinguishes_from": experiment.distinguishes_from,
                "expected_if_true": experiment.expected_if_true,
                "expected_if_false": experiment.expected_if_false,
                "feasibility": experiment.feasibility,
            }
            for experiment in sorted(
                mechanism.discriminating_experiments, key=lambda e: e.experiment_id
            )
        ],
    }


def _link_summary(link: MechanismLink) -> str:
    return (
        f"{link.link_id}: {link.source_node_id} -{link.relation}-> {link.target_node_id} "
        f"({link.direction.value}, {link.tier.value}, {link.strength.value}) — "
        f"unresolved: {link.uncertainty}"
    )


def _link_context(link: MechanismLink, checker: AssertionChecker) -> dict[str, object]:
    claim = (
        f"{link.source_node_id} {link.relation} {link.target_node_id} "
        f"(direction: {link.direction.value})."
    )
    if link.evidence_ids:
        statement = checker.check(
            Assertion(text=claim, tier=link.tier, evidence_ids=tuple(sorted(link.evidence_ids)))
        ).rendered()
        basis = ", ".join(sorted(link.evidence_ids))
    else:
        # Rendered as a gap, not as a claim: see the module docstring.
        statement = f"{claim} [UNSOURCED — recorded as a gap, not asserted]"
        basis = "none"
    return {
        "link_id": link.link_id,
        "statement": statement,
        "basis": basis,
        "is_directly_demonstrated": link.is_directly_demonstrated,
        "uncertainty": link.uncertainty,
        "contradictions": sorted(link.contradicting_evidence_ids),
    }


def _pair_context(pair: CandidatePair | None) -> dict[str, object] | None:
    if pair is None:
        return None
    return {
        "pair_id": pair.pair_id,
        "gene_symbol": pair.gene_symbol,
        "inheritance_model": pair.inheritance_model.value,
        "phase_status": pair.phase.status.value,
        "composite_score": f"{pair.composite_score:.3f}",
        "variant_ids": list(pair.variant_ids),
        "flags": sorted(pair.flags),
    }


def _drug_context(
    drug: DrugHypothesis,
    mechanism: MechanismHypothesis,
    resolver: AssertionResolver,
    checker: AssertionChecker,
    *,
    require_rationale: bool,
) -> dict[str, object]:
    evidence = resolver.resolve(drug.evidence_ids)
    contradictions = resolver.resolve(drug.contradicting_evidence_ids)
    rationale = _drug_rationale(drug, checker, evidence, require_rationale=require_rationale)
    return {
        "drug_id": drug.drug_id,
        "name": drug.name,
        "approved_name": drug.approved_name or NOT_RECORDED,
        "approval_status": drug.approval_status.value,
        "is_repurposable": drug.is_repurposable,
        "intervention_class": drug.intervention_class.value,
        "score": f"{drug.score:.3f}",
        "rank": drug.rank,
        "rationale": rationale,
        "concerns": [reason.value for reason in drug.concerns],
        "direction_verdict": _direction_verdict(drug),
        "target_in_mechanism": drug.target_node_id in mechanism.node_ids,
        "questions": _drug_questions(drug, mechanism),
        "safety_concerns": [
            {
                "concern_id": concern.concern_id,
                "description": concern.description,
                "severity": concern.severity,
                "population": concern.population,
                "is_disqualifying": concern.is_disqualifying,
            }
            for concern in sorted(drug.safety_concerns, key=lambda c: c.concern_id)
        ],
        "evidence": [_evidence_context(item) for item in evidence],
        "contradictions": [_evidence_context(item) for item in contradictions],
    }


def _drug_rationale(
    drug: DrugHypothesis,
    checker: AssertionChecker,
    evidence: Sequence[EvidenceItem],
    *,
    require_rationale: bool,
) -> str:
    """The one sentence that says why this compound is here at all.

    For a presented candidate this is a claim about biology and must be sourced
    (GP-10). For a symptomatic or supportive agent the statement is about clinical
    management rather than about the disease mechanism, so it is rendered with its
    citations where they exist and as an explicit gap where they do not.
    """
    claim = (
        f"{drug.name} acts on {drug.target} by {drug.mechanism_of_action}, which bears on "
        f"mechanism node {drug.target_node_id}."
    )
    if not drug.evidence_ids and not require_rationale:
        return f"{claim} [UNSOURCED — recorded as a gap, not asserted]"
    return checker.check(
        Assertion(
            text=claim,
            tier=weakest_tier(evidence) if evidence else AssertionTier.SPECULATION,
            evidence_ids=tuple(sorted(drug.evidence_ids)),
        )
    ).rendered()


def _direction_verdict(drug: DrugHypothesis) -> str:
    agreement = drug.directions_agree
    if agreement is True:
        return (
            f"AGREES — required {drug.required_direction.value}, observed "
            f"{drug.observed_direction.value}."
        )
    if agreement is None:
        return (
            f"CANNOT BE DETERMINED — required {drug.required_direction.value}, observed "
            f"{drug.observed_direction.value}. An unsigned direction is not agreement and "
            "is never counted as support (GP-16, ASSUMPTION-DRUG-02)."
        )
    return (
        f"DISAGREES — required {drug.required_direction.value}, observed "
        f"{drug.observed_direction.value}. Disqualifying (GP-16)."
    )


def _drug_questions(
    drug: DrugHypothesis, mechanism: MechanismHypothesis
) -> list[dict[str, object]]:
    """The eight mandatory questions, as labelled fields."""
    paediatric = drug.pediatric_evidence
    node_label = next(
        (node.label for node in mechanism.nodes if node.node_id == drug.target_node_id),
        "node not present in this mechanism",
    )
    answers = (
        f"{drug.target} — mechanism node {drug.target_node_id} ({node_label}).",
        _required_direction_answer(drug, mechanism),
        _direction_verdict(drug),
        (
            f"{drug.strongest_evidence_type.value}; "
            f"{'direct' if drug.is_direct_evidence else 'INDIRECT'} evidence for this drug on "
            f"this target; in vivo: {format_cell(drug.has_in_vivo_evidence)}."
        ),
        _concentration_answer(drug),
        (
            f"{'Yes' if paediatric.has_pediatric_exposure else 'NO PAEDIATRIC EXPOSURE RECORDED'}"
            f" — youngest age studied: {paediatric.youngest_age_studied or NOT_RECORDED}; "
            f"indication: {paediatric.indication or NOT_RECORDED}; tolerability: "
            f"{paediatric.tolerability_summary or NOT_RECORDED}. {paediatric.caveat}"
        ),
        _instability_answer(drug),
        drug.proposed_validation_experiment,
    )
    return [
        {"number": index, "question": question, "answer": answer}
        for index, (question, answer) in enumerate(zip(DRUG_QUESTIONS, answers, strict=True), 1)
    ]


def _required_direction_answer(drug: DrugHypothesis, mechanism: MechanismHypothesis) -> str:
    """Mandatory question 2, stating the derivation that was actually used.

    The requirement is only the inverse of ``disease_direction`` for the mechanism's
    designated therapeutic target. For any other node it comes from that node's own
    ``state_in_patient``, and for a node off the chain it does not exist at all.
    Printing "the required correction is its inverse" for all three cases stated a
    derivation the pipeline had not performed — and did so about the compounds whose
    grounds are weakest, since those are exactly the ones acting away from the target.
    """
    target_node: MechanismNode | None = next(
        (node for node in mechanism.nodes if node.node_id == drug.target_node_id), None
    )
    if drug.target_node_id == mechanism.therapeutic_target_node_id:
        return (
            f"{drug.required_direction.value}. `{drug.target_node_id}` is this mechanism's "
            f"designated therapeutic target, whose disease direction is "
            f"{mechanism.disease_direction.value}; the required correction is its inverse."
        )
    if target_node is None:
        return (
            f"{drug.required_direction.value}. `{drug.target_node_id}` is not a node of "
            f"{mechanism.mechanism_id}, so no required direction can be derived for it at "
            f"all. The mechanism's own required correction "
            f"({mechanism.required_correction.value}) applies at "
            f"`{mechanism.therapeutic_target_node_id}` and NOT to this agent."
        )
    if not target_node.deviation_is_pathological:
        return (
            f"{drug.required_direction.value}. `{drug.target_node_id}` is not the therapeutic "
            f"target (`{mechanism.therapeutic_target_node_id}`), and its deviation is recorded "
            "as COMPENSATORY rather than pathological. No corrective direction follows from a "
            "compensatory node's state — pushing it back towards wild type would suppress a "
            "protective response — so the requirement is undetermined rather than derived."
        )
    return (
        f"{drug.required_direction.value}. `{drug.target_node_id}` is not the therapeutic "
        f"target (`{mechanism.therapeutic_target_node_id}`), so the requirement is derived "
        f"from that node's own state in the patient ({target_node.state_in_patient.value}) "
        f"and NOT from the mechanism's disease direction "
        f"({mechanism.disease_direction.value})."
    )


def _concentration_answer(drug: DrugHypothesis) -> str:
    profile = drug.pharmacokinetics
    achievable = profile.concentration_achievable
    detail = (
        f"achievable plasma {format_cell(profile.achievable_plasma_concentration_um)} uM vs "
        f"required {format_cell(profile.required_effective_concentration_um)} uM; route "
        f"{profile.route or NOT_RECORDED}; CNS penetrant: {format_cell(profile.cns_penetrant)}."
    )
    if achievable is None:
        return (
            f"CANNOT BE DETERMINED — {detail} An unknown exposure is recorded as a gap, not "
            "resolved optimistically (ASSUMPTION-DRUG-05)."
        )
    if achievable:
        return f"Yes — {detail}"
    return (
        f"NO — {detail} A compound effective in culture but not reachable in vivo is not a "
        "therapy, however elegant the mechanism."
    )


def _instability_answer(drug: DrugHypothesis) -> str:
    worsens = drug.worsens_chromosomal_instability
    if worsens is None:
        return (
            "UNASSESSED. In a chromosomal-instability disorder this is itself a blocking "
            "gap, not a pass (ASSUMPTION-DRUG-07). It must be answered before any further "
            "consideration."
        )
    if worsens:
        return (
            "YES — this agent could increase aneuploidy or cancer susceptibility. "
            "Disqualifying in this disease context."
        )
    return "No increase in aneuploidy or cancer susceptibility is expected on current evidence."


def _rejection_context(drug: DrugHypothesis) -> dict[str, object]:
    reasons = sorted(reason.value for reason in drug.rejection_reasons)
    return {
        "drug_id": drug.drug_id,
        "name": drug.name,
        "intervention_class": drug.intervention_class.value,
        "approval_status": drug.approval_status.value,
        "target": drug.target,
        "target_node_id": drug.target_node_id,
        "required_direction": drug.required_direction.value,
        "observed_direction": drug.observed_direction.value,
        "direction_verdict": _direction_verdict(drug),
        "reasons": reasons,
        "rationale": drug.rejection_rationale or NOT_RECORDED,
        "reversal_conditions": [
            REVERSAL_CONDITIONS[reason]
            for reason in sorted(set(drug.rejection_reasons), key=lambda r: r.value)
        ],
        "safety_concerns": [
            f"{concern.severity}: {concern.description} (observed in {concern.population})"
            for concern in sorted(drug.safety_concerns, key=lambda c: c.concern_id)
        ],
    }


def _evidence_context(item: EvidenceItem) -> dict[str, object]:
    return {
        "evidence_id": item.evidence_id,
        "claim": item.claim,
        "tier": item.tier.value,
        "marker": TIER_MARKERS[item.tier],
        "strength": item.strength.value,
        "evidence_type": item.evidence_type.value,
        "direction": item.direction.value,
        "limitations": item.limitations,
        "source": (
            item.citation.key if item.citation is not None else f"{item.tool} {item.tool_version}"
        ),
    }


def _caveats(mechanism: MechanismHypothesis) -> list[str]:
    caveats = [
        UNCALIBRATED_WEIGHTS,
        "Direction of effect is the primary gate. A compound acting on the right node in "
        "the wrong direction is contraindicated, not merely unhelpful (ASSUMPTION-DRUG-01).",
        "Approval status and direction are separate rejections and are never conflated. An "
        "investigational agent may be a valid research probe and still not be a repurposing "
        "candidate (ASSUMPTION-DRUG-03).",
        "Paediatric tolerability does not transfer across populations: exposure in a "
        "paediatric oncology cohort is not evidence of safety in a germline "
        "chromosomal-instability population (ASSUMPTION-DRUG-06).",
        "Drug, literature and mechanism sources in this pipeline may be local synthetic "
        "substitutes graded in docs/maturity-ledger.md. A synthetic substitute is never "
        "biologically valid evidence (GP-20).",
    ]
    if not mechanism.is_fully_demonstrated:
        caveats.append(
            "Mechanism links inferred rather than directly demonstrated: "
            f"{len(mechanism.inferred_links)} of {len(mechanism.links)}. Every drug "
            "hypothesis below inherits that uncertainty (ASSUMPTION-MECHANISM-01)."
        )
    return caveats


__all__ = [
    "DRUG_QUESTIONS",
    "NOT_MEDICAL_ADVICE",
    "REVERSAL_CONDITIONS",
    "UNCALIBRATED_WEIGHTS",
    "build_drug_report",
    "build_mechanism_report",
    "build_rejection_record",
    "build_track2_report",
]
