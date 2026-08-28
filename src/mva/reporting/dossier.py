"""The Track 1 candidate dossier: the answer to "why is this ranked first?".

The submission CSV is twelve columns of coordinates and one number. It cannot say
whether the phase is unknown, whether anything contradicts the hypothesis, or what
experiment would settle it. The dossier is where that lives, and it is structured
as nine fixed questions rather than prose because prose lets a weak answer hide in
a strong paragraph. A question with no answer renders as a question with no
answer.

The nine questions come from the project's scientific safeguards:

1. Is each variant analytically credible?
2. Is each sufficiently rare?
3. Are both predicted to damage the same gene?
4. Is phase confirmed, likely or unknown?
5. Does gene dysfunction explain the phenotype?
6. Is the mechanism direct or inferred?
7. What evidence contradicts the hypothesis?
8. What observation would falsify it?
9. What experiment would discriminate it from the next candidate?

Answers are readouts of recorded state — component scores, genotype quality,
phase evidence — and are labelled as such. Where evidence IDs are cited, the
answer is additionally run through the GP-10 gate and rendered with the tier
marker of its **weakest** citation, so a sentence resting partly on inference is
never printed as though it were measured.

**Sensitivity.** A dossier contains patient genotypes and quality metrics. It is
a SENSITIVE artifact: it stays in the external workspace and is never a candidate
for public export (GP-40, GP-43). The rendered header says so on the first line.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING

from mva.clock import Clock
from mva.models import (
    CIS_STATES,
    AssertionTier,
    CandidatePair,
    EvidenceCategory,
    EvidenceItem,
    InheritanceModel,
    PhaseStatus,
    VariantRecord,
)
from mva.reporting.assertions import TIER_MARKERS, Assertion, AssertionChecker, weakest_tier
from mva.reporting.render import NOT_RECORDED, format_cell, markdown_table, render_template
from mva.reporting.track1 import SubmissionComposition, composite_to_epcr, truncation_notice

if TYPE_CHECKING:  # pragma: no cover - typing only
    from mva.evidence.ledger import AssertionResolver

#: The nine mandatory questions, with the evidence category whose items speak to
#: each. ``None`` means the question is answered from recorded pipeline state and
#: from the contradiction set rather than from a category of supporting evidence.
DOSSIER_QUESTIONS: tuple[tuple[str, str, EvidenceCategory | None], ...] = (
    ("Q1", "Is each variant analytically credible?", EvidenceCategory.ANALYTICAL),
    ("Q2", "Is each variant sufficiently rare?", EvidenceCategory.POPULATION),
    ("Q3", "Are both variants predicted to damage the same gene?", EvidenceCategory.CONSEQUENCE),
    ("Q4", "Is phase confirmed, likely or unknown?", EvidenceCategory.INHERITANCE),
    ("Q5", "Does gene dysfunction explain the phenotype?", EvidenceCategory.PHENOTYPE),
    ("Q6", "Is the mechanism direct or inferred?", EvidenceCategory.MECHANISM),
    ("Q7", "What evidence contradicts the hypothesis?", None),
    ("Q8", "What observation would falsify this hypothesis?", None),
    ("Q9", "What experiment would discriminate it from the next candidate?", None),
)

#: Rendered when a question has no supporting evidence behind it. Phrased as a
#: gap, never as reassurance: absence of evidence is not evidence of absence
#: (GP-14).
NO_EVIDENCE_BASIS = (
    "no evidence cited — recorded as a gap, not as support (GP-10/GP-14); the "
    "statement above is a readout of pipeline state, not a sourced claim"
)

#: The question whose evidence is the contradiction set rather than a supporting
#: category. Named rather than hard-coded at the use site so the coupling is visible.
CONTRADICTION_QUESTION_ID = "Q7"

_SCORE_HEADERS = ("Component", "Score")
_VARIANT_HEADERS = (
    "Variant",
    "Coordinate",
    "Zygosity",
    "GT",
    "DP",
    "GQ",
    "Allele balance",
    "FILTER",
    "Worst impact in gene",
    "Max population AF",
)


def build_candidate_dossier(
    pairs: Sequence[CandidatePair],
    *,
    submission: SubmissionComposition,
    resolver: AssertionResolver,
    clock: Clock,
    top_n: int = 5,
) -> str:
    """Render the dossier for the top ``top_n`` candidates.

    Raises :class:`~mva.errors.UnsourcedAssertionError` if any pair cites evidence
    the ledger cannot resolve. That is the intended failure: a dossier is the
    document a reviewer trusts, and a dangling citation in it is indistinguishable
    from a fabricated one.

    ``submission`` is the composition that produced the submission CSV, and it is
    required rather than optional. The dossier's job here is to tell a reader why
    the submission is shorter than the ranked list above it, and the only honest
    source for that is the submission itself. This function used to derive the
    notice from ``len(pairs)`` and the format cap, which reported a submitted count
    it had not seen and blamed the format for rows that composition had dropped on
    their contents. A parameter that can be omitted is a parameter that will be.

    Raises:
        ValueError: if ``submission`` was composed from a different number of
            candidates than ``pairs`` holds — the notice would then describe a
            submission this dossier is not the ranked list for.
    """
    if top_n < 1:
        msg = f"top_n={top_n} must be at least 1."
        raise ValueError(msg)

    checker = AssertionChecker(resolver, strict=True)
    ordered = sorted(pairs, key=lambda pair: pair.sort_key())
    selected = ordered[:top_n]

    if submission.total_candidates != len(ordered):
        msg = (
            f"the submission was composed from {submission.total_candidates} candidates "
            f"but this dossier ranks {len(ordered)}. The truncation notice would then "
            "describe a submission that is not the one this ranking produced."
        )
        raise ValueError(msg)

    context: dict[str, object] = {
        "generated_at": clock.now().isoformat(),
        "total_candidates": len(ordered),
        "shown": len(selected),
        "top_n": top_n,
        "truncation_notice": truncation_notice(submission),
        "candidates": [
            _candidate_context(pair, rank=index, resolver=resolver, checker=checker)
            for index, pair in enumerate(selected, start=1)
        ],
        "tier_legend": [
            {"marker": TIER_MARKERS[tier], "tier": tier.value} for tier in AssertionTier
        ],
        "caveats": _caveats(),
    }
    return render_template("dossier.md.j2", context)


def _caveats() -> list[str]:
    """Caveats that apply to every dossier, stated where they cannot be skipped."""
    return [
        "Component scores and the composite are **uncalibrated heuristics**, not "
        "probabilities. The weights were chosen by reasoning about a severe paediatric "
        "recessive disorder and are fitted to no labelled dataset (ASSUMPTION-SCORING-01).",
        "EPCR in the submission is a **rank-ordering confidence**, not a calibrated "
        "probability. It is an affine, strictly increasing map of the composite score onto "
        "the range the challenge scorer accepts.",
        "Annotation, frequency and phenotype sources may be local synthetic substitutes. "
        "Their grade is recorded in docs/maturity-ledger.md; a synthetic substitute is "
        "never biologically valid evidence (GP-20).",
        "Ranking is not diagnosis. Nothing here is medical advice, and no candidate should "
        "be acted on without orthogonal confirmation and clinical review.",
        "Absence of a finding in this document means it was not assessed or not recorded — "
        "it is never evidence that the finding is absent (GP-14).",
    ]


def _candidate_context(
    pair: CandidatePair,
    *,
    rank: int,
    resolver: AssertionResolver,
    checker: AssertionChecker,
) -> dict[str, object]:
    supporting = (
        resolver.require(
            f"supporting evidence cited by candidate {pair.pair_id}",
            pair.supporting_evidence_ids,
        )
        if pair.supporting_evidence_ids
        else ()
    )
    contradicting = (
        resolver.require(
            f"contradicting evidence cited by candidate {pair.pair_id}",
            pair.contradicting_evidence_ids,
        )
        if pair.contradicting_evidence_ids
        else ()
    )

    by_category: dict[EvidenceCategory, list[EvidenceItem]] = {}
    for item in supporting:
        by_category.setdefault(item.category, []).append(item)

    answers = _answers(pair, contradicting)
    questions: list[dict[str, object]] = []
    for question_id, question, category in DOSSIER_QUESTIONS:
        if question_id == CONTRADICTION_QUESTION_ID:
            # Q7 is the one question whose evidence is the contradicting set. Citing
            # it here is what makes the answer carry the [contradicted] marker.
            items: tuple[EvidenceItem, ...] = contradicting
        elif category is not None:
            items = tuple(by_category.get(category, ()))
        else:
            items = ()
        questions.append(
            _question_context(question_id, question, answers[question_id], items, checker)
        )

    return {
        "rank": rank,
        "pair_id": pair.pair_id,
        "gene_symbol": pair.gene_symbol,
        "inheritance_model": pair.inheritance_model.value,
        "is_pair": pair.is_pair,
        "phase_status": pair.phase.status.value,
        "phase_disqualifying": pair.phase_is_disqualifying,
        "composite_score": f"{pair.composite_score:.3f}",
        "epcr": f"{composite_to_epcr(pair.composite_score):.4f}",
        "flags": sorted(pair.flags),
        "rank_rationale": pair.rank_rationale or NOT_RECORDED,
        "score_table": markdown_table(_SCORE_HEADERS, _score_rows(pair)),
        "variants_table": markdown_table(_VARIANT_HEADERS, _variant_rows(pair)),
        "questions": questions,
        "open_questions": [
            {
                "question_id": question.question_id,
                "question": question.question,
                "why_it_matters": question.why_it_matters,
                "resolving_test": question.resolving_test,
                "blocking": question.blocking,
            }
            for question in sorted(pair.missing_evidence, key=lambda q: q.question_id)
        ],
        "blocking_count": len(pair.blocking_questions),
        "contradictions": [_evidence_context(item) for item in contradicting],
        "supporting_count": len(supporting),
    }


def _question_context(
    question_id: str,
    question: str,
    answer: str,
    items: Sequence[EvidenceItem],
    checker: AssertionChecker,
) -> dict[str, object]:
    """Attach evidence to an answer and push it through the GP-10 gate."""
    if not items:
        return {
            "id": question_id,
            "question": question,
            "answer": answer,
            "basis": NO_EVIDENCE_BASIS,
            "citations": "none",
            "evidence": [],
        }
    checked = checker.check(
        Assertion(
            text=answer,
            tier=weakest_tier(items),
            evidence_ids=tuple(sorted(item.evidence_id for item in items)),
        )
    )
    return {
        "id": question_id,
        "question": question,
        "answer": checked.rendered(),
        "basis": f"tier of weakest citation: {weakest_tier(items).value}",
        "citations": checked.citation(),
        "evidence": [_evidence_context(item) for item in items],
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
        "source": item.citation.key
        if item.citation is not None
        else f"{item.tool} {item.tool_version}",
    }


def _score_rows(pair: CandidatePair) -> list[list[object]]:
    """The score vector, in declaration order — never collapsed to the composite."""
    rows: list[list[object]] = [
        [name.replace("_", " "), value] for name, value in pair.scores.as_dict().items()
    ]
    rows.append(["**composite (weighted)**", pair.composite_score])
    return rows


def _variant_rows(pair: CandidatePair) -> list[list[object]]:
    return [_variant_row(label, variant, pair.gene_symbol) for label, variant in _labelled(pair)]


def _labelled(pair: CandidatePair) -> tuple[tuple[str, VariantRecord], ...]:
    if pair.variant_b is None:
        return (("A", pair.variant_a),)
    return (("A", pair.variant_a), ("B", pair.variant_b))


def _variant_row(label: str, variant: VariantRecord, gene: str) -> list[object]:
    genotype = variant.genotype
    frequency = variant.max_allele_frequency()
    impact = variant.worst_impact_for_gene(gene)
    return [
        label,
        variant.variant_id,
        genotype.zygosity.value,
        genotype.genotype_string,
        genotype.depth,
        genotype.genotype_quality,
        genotype.allele_balance,
        variant.filter_status.value,
        impact.value if impact is not None else NOT_RECORDED,
        (
            f"{frequency.allele_frequency:.2e} ({frequency.provenance_key})"
            if frequency is not None
            else "no frequency data"
        ),
    ]


# ---------------------------------------------------------------------------
# Answers. Each is a readout of what the pipeline actually recorded.
# ---------------------------------------------------------------------------


def _answers(pair: CandidatePair, contradicting: Sequence[EvidenceItem]) -> Mapping[str, str]:
    return {
        "Q1": _answer_analytical(pair),
        "Q2": _answer_rarity(pair),
        "Q3": _answer_consequence(pair),
        "Q4": _answer_phase(pair),
        "Q5": _answer_phenotype(pair),
        "Q6": _answer_mechanism(pair),
        "Q7": _answer_contradictions(pair, contradicting),
        "Q8": _answer_falsification(pair),
        "Q9": _answer_discriminating(pair),
    }


def _answer_analytical(pair: CandidatePair) -> str:
    parts: list[str] = [f"Analytical validity scored {pair.scores.analytical_validity:.3f}."]
    for label, variant in _labelled(pair):
        genotype = variant.genotype
        flags = ", ".join(variant.qc_flags) if variant.qc_flags else "none"
        parts.append(
            f"Variant {label} ({variant.variant_id}): FILTER={variant.filter_status.value}, "
            f"DP={format_cell(genotype.depth)}, GQ={format_cell(genotype.genotype_quality)}, "
            f"allele balance={format_cell(genotype.allele_balance)}, QC flags: {flags}."
        )
    parts.append(
        "Quality flags down-rank; they never delete a candidate (GP-13). A skewed allele "
        "balance may be somatic mosaicism rather than an artifact in this disease context "
        "(ASSUMPTION-MOSAIC-01)."
    )
    return " ".join(parts)


def _answer_rarity(pair: CandidatePair) -> str:
    parts: list[str] = [f"Rarity scored {pair.scores.rarity:.3f}."]
    for label, variant in _labelled(pair):
        frequency = variant.max_allele_frequency()
        if frequency is None:
            parts.append(
                f"Variant {label}: no population frequency recorded. Scored at the "
                "mid-range default — absence of frequency data is not evidence of rarity "
                "(GP-14, ASSUMPTION-FREQUENCY-01)."
            )
        else:
            parts.append(
                f"Variant {label}: maximum population AF {frequency.allele_frequency:.3e} "
                f"in {frequency.provenance_key}"
                + (
                    f", {frequency.homozygote_count} homozygote(s)."
                    if frequency.homozygote_count is not None
                    else "."
                )
            )
    parts.append("Maximum across populations is used, not global AF (ASSUMPTION-FREQUENCY-02).")
    return " ".join(parts)


def _answer_consequence(pair: CandidatePair) -> str:
    parts: list[str] = [
        f"Molecular consequence scored {pair.scores.molecular_consequence:.3f} against "
        f"gene {pair.gene_symbol}."
    ]
    for label, variant in _labelled(pair):
        consequences = variant.consequences_for_gene(pair.gene_symbol)
        if not consequences:
            parts.append(
                f"Variant {label} has no recorded consequence on {pair.gene_symbol}: the "
                "two variants are not demonstrably hitting the same gene product."
            )
            continue
        worst = variant.worst_impact_for_gene(pair.gene_symbol)
        terms = ", ".join(sorted({term for csq in consequences for term in csq.consequence_terms}))
        hgvs = ", ".join(
            sorted({csq.hgvs_p or csq.hgvs_c or csq.transcript_id for csq in consequences})
        )
        parts.append(
            f"Variant {label}: worst impact {worst.value if worst else NOT_RECORDED} across "
            f"{len(consequences)} transcript(s); consequence terms: {terms}; {hgvs}."
        )
    parts.append(
        "Severity is taken across all transcripts, not the canonical one alone "
        "(ASSUMPTION-TRANSCRIPT-01). In-silico prediction is not proof of causality "
        "(ASSUMPTION-PREDICTION-01)."
    )
    return " ".join(parts)


def _answer_phase(pair: CandidatePair) -> str:
    phase = pair.phase
    detail = (
        f"Phase is {phase.status.value} (method: {phase.method}"
        + (
            f", {phase.supporting_reads} supporting read(s)"
            if phase.supporting_reads is not None
            else ""
        )
        + (f", sites {phase.distance_bp} bp apart" if phase.distance_bp is not None else "")
        + ")."
    )
    if phase.status in CIS_STATES:
        verdict = (
            "The variants are called in cis. One allele remains intact, so a "
            "compound-heterozygous recessive mechanism does not apply. The candidate is "
            "down-weighted to near-zero rather than deleted, because short-read phasing "
            "can be wrong (GP-13, ASSUMPTION-PHASE-01)."
        )
    elif phase.status is PhaseStatus.UNKNOWN:
        verdict = (
            "Phase is UNKNOWN and is not upgraded. Two heterozygous variants in one gene "
            "are a compound heterozygote only if in trans; assuming trans without evidence "
            "is the error this pipeline refuses to make (GP-15, ASSUMPTION-PHASE-01). "
            "Resolving test: parental segregation, read-backed phasing, or long reads."
        )
    elif phase.status is PhaseStatus.TRANS_LIKELY:
        verdict = (
            "Trans is statistically supported but not proven; the recessive model is "
            "plausible, not established."
        )
    else:
        verdict = "Trans is confirmed; the biallelic requirement is met."
    if pair.inheritance_model is not InheritanceModel.COMPOUND_HETEROZYGOUS:
        verdict += (
            f" The proposed model is {pair.inheritance_model.value}, for which phase is "
            "informative but not decisive."
        )
    return (
        f"{detail} {verdict} Inheritance consistency scored "
        f"{pair.scores.inheritance_consistency:.3f}, which already carries the phase penalty."
    )


def _answer_phenotype(pair: CandidatePair) -> str:
    return (
        f"Phenotype similarity scored {pair.scores.phenotype_similarity:.3f} between the "
        f"phenotype associated with {pair.gene_symbol} and the proband's recorded HPO "
        "profile. Only EXCLUDED terms contribute negative evidence; NOT_ASSESSED terms are "
        "excluded from the denominator entirely, so a gene is never penalised for features "
        "nobody checked (GP-14, ASSUMPTION-PHENOTYPE-01)."
    )


def _answer_mechanism(pair: CandidatePair) -> str:
    return (
        f"Mechanistic relevance scored {pair.scores.mechanistic_relevance:.3f} and evidence "
        f"quality {pair.scores.evidence_quality:.3f}. The tier marker on this answer states "
        "whether the mechanistic link is directly demonstrated (literature/observed) or "
        "inferred by this pipeline from pathway membership. Mechanism is resolved per gene, "
        "never per disease label (ASSUMPTION-MECHANISM-03); a chain is only as strong as its "
        "weakest link (ASSUMPTION-MECHANISM-01)."
    )


def _answer_contradictions(pair: CandidatePair, contradicting: Sequence[EvidenceItem]) -> str:
    if not contradicting:
        return (
            "No contradicting evidence is recorded for this candidate. That is an absence of "
            "recorded contradiction, not a demonstration of consistency (GP-14). The "
            f"contradiction penalty applied was {pair.scores.contradiction_penalty:.3f}."
        )
    claims = " ".join(f"({item.evidence_id}) {item.claim}" for item in contradicting)
    return (
        f"{len(contradicting)} contradicting item(s) are retained and subtracted from the "
        f"score (penalty {pair.scores.contradiction_penalty:.3f}); they are never discarded "
        f"(GP-19). {claims}"
    )


def _answer_falsification(pair: CandidatePair) -> str:
    observations: list[str] = []
    if pair.inheritance_model is InheritanceModel.COMPOUND_HETEROZYGOUS:
        observations.append(
            "phasing (trio, read-backed or long-read) placing both variants on the same "
            "haplotype, which leaves one intact allele"
        )
        observations.append("an unaffected parent or sibling carrying the same biallelic genotype")
    else:
        observations.append(
            "segregation showing the genotype in an unaffected first-degree relative"
        )
    observations.append(
        f"a functional assay showing normal {pair.gene_symbol} product level and activity in "
        "patient-derived cells"
    )
    observations.append(
        "population data placing either allele above the frequency compatible with a severe "
        "recessive disorder"
    )
    for question in sorted(pair.missing_evidence, key=lambda q: q.question_id):
        if question.blocking:
            observations.append(f"{question.resolving_test} returning a negative result")
    return "This hypothesis is falsified by any of: " + "; ".join(observations) + "."


def _answer_discriminating(pair: CandidatePair) -> str:
    experiment = pair.discriminating_experiment
    next_test = pair.recommended_next_test
    if experiment is None:
        return (
            f"No discriminating experiment is recorded for this candidate. The recommended "
            f"next test is: {next_test}. Until a discriminating experiment is specified, this "
            "candidate cannot be separated from the next-ranked one on evidence, only on the "
            "heuristic score."
        )
    return (
        f"Discriminating experiment: {experiment} Highest-value next test for this candidate: "
        f"{next_test}"
    )


__all__ = ["DOSSIER_QUESTIONS", "NO_EVIDENCE_BASIS", "build_candidate_dossier"]
