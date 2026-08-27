"""Component scoring for candidate pairs.

The score is a **vector, not a scalar**. Seven positive components and one
subtracted contradiction penalty are computed independently, each returns the
sentence that explains it, and every one of them is emitted as an
``EvidenceItem`` so that "why is this ranked first?" is answerable from the
evidence store alone.

Honesty constraints this module is built around:

* Every weight and cut-point is an **uncalibrated heuristic** (GP-32). Nothing
  here was fitted to labelled outcomes and no component is a probability.
* Absence of population-frequency data is scored at a documented mid-range
  value, never as rarity (GP-14).
* Phase is taken as given by :mod:`mva.prioritization.pairing`; unknown phase
  is penalised, never upgraded (GP-15).
* Contradictions are recorded as evidence, not quietly folded into the number
  (GP-19), and every emitted item states its limitations (GP-17).
* ``phenotype_score`` and ``mechanism_score`` are computed by other stages and
  passed in (GP-03). When a caller has nothing to supply it passes the
  documented neutral value, and the rationale says so.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import pairwise

from mva import __version__
from mva.clock import Clock
from mva.config import FrequencyThresholds, PhaseWeights, QualityThresholds, ScoringWeights
from mva.models.base import AssertionTier
from mva.models.evidence import (
    EvidenceCategory,
    EvidenceDirection,
    EvidenceItem,
    EvidenceStrength,
    EvidenceType,
    make_evidence_id,
)
from mva.models.pair import CIS_STATES, ComponentScores, InheritanceModel, OpenQuestion, PhaseStatus
from mva.models.variant import (
    FilterStatus,
    ImpactSeverity,
    PopulationFrequency,
    VariantRecord,
    Zygosity,
)
from mva.prioritization.filters import (
    FLAG_BENIGN_CONSEQUENCE,
    FLAG_COMMON_VARIANT,
    FLAG_LOW_QUALITY_CALL,
    FLAG_NO_FREQUENCY_DATA,
)
from mva.prioritization.pairing import PairCandidate

TOOL_NAME = "mva.prioritization.scoring"
TOOL_VERSION = __version__

#: Value a caller supplies when the phenotype stage produced nothing for this
#: gene. Mid-range on purpose: no phenotype comparison is neither support nor
#: contradiction (GP-14).
NEUTRAL_PHENOTYPE_SCORE = 0.5

#: Same contract for the mechanism stage.
NEUTRAL_MECHANISM_SCORE = 0.5

#: Shared caveat on every item this module emits (GP-17).
HEURISTIC_LIMITATION = (
    "Uncalibrated pipeline heuristic: the cut-points and weights were chosen by "
    "reasoning about severe paediatric recessive disease, not fitted to any labelled "
    "outcome dataset (GP-32). This number is not a probability of pathogenicity and "
    "carries no claim of clinical validity."
)

# ---------------------------------------------------------------------------
# Analytical validity
# ---------------------------------------------------------------------------

#: Multiplier when the caller filtered the site out.
_FILTERED_MULTIPLIER = 0.25
#: Multiplier when FILTER was absent entirely (no caller opinion recorded).
_FILTER_MISSING_MULTIPLIER = 0.70
#: Multiplier when depth or GQ is simply not reported.
_UNREPORTED_METRIC_MULTIPLIER = 0.70
#: Floor for a proportional depth/GQ shortfall, so a terrible call is never zero.
_METRIC_FLOOR = 0.15
#: Allele balance below the mosaic floor: likely artifact.
_ALLELE_BALANCE_ARTIFACT_MULTIPLIER = 0.40
#: Allele balance in the mosaic window: down-weighted, not dismissed.
_MOSAIC_MULTIPLIER = 0.75
#: Allele balance unknowable.
_ALLELE_BALANCE_UNKNOWN_MULTIPLIER = 0.85
#: Applied when an upstream stage flagged the call for a reason we cannot re-derive.
_UPSTREAM_FLAG_MULTIPLIER = 0.60

# ---------------------------------------------------------------------------
# Rarity: piecewise-linear in log10(AF) between the configured cut-points.
# ---------------------------------------------------------------------------

_RARITY_AT_ULTRA_RARE = 1.00
_RARITY_AT_RARE = 0.80
_RARITY_AT_LOW_FREQUENCY = 0.50
_RARITY_AT_MAX_PLAUSIBLE = 0.15
_RARITY_FLOOR = 0.02
#: Frequency at which the rarity score reaches its floor.
_RARITY_FLOOR_AF = 0.5

# ---------------------------------------------------------------------------
# Molecular consequence
# ---------------------------------------------------------------------------

_IMPACT_BASE: dict[ImpactSeverity, float] = {
    ImpactSeverity.HIGH: 0.90,
    ImpactSeverity.MODERATE: 0.50,
    ImpactSeverity.LOW: 0.15,
    ImpactSeverity.MODIFIER: 0.05,
}
#: No consequence annotation at all for this gene: unknown, not benign.
_IMPACT_UNANNOTATED = 0.10

#: Sequence Ontology terms that remove or truncate the gene product outright.
LOF_TERMS: frozenset[str] = frozenset(
    {
        "frameshift_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "start_lost",
        "stop_gained",
        "transcript_ablation",
    }
)
_LOF_UPLIFT = 0.04
_SPLICEAI_UPLIFT = 0.05
_SPLICEAI_FLOOR = 0.5
_CADD_UPLIFT = 0.04
_CADD_FLOOR = 20.0
_CADD_CEILING = 35.0
_REVEL_UPLIFT = 0.03
_REVEL_FLOOR = 0.5

# ---------------------------------------------------------------------------
# Inheritance
# ---------------------------------------------------------------------------

_INHERITANCE_COMPOUND_HET = 1.00
_INHERITANCE_HOMOZYGOUS = 0.90
_INHERITANCE_MIXED_PAIR = 0.60
#: A lone heterozygote with no de-novo evidence and no second hit is a weak
#: inheritance hypothesis whichever model you have in mind: recessive needs a
#: second allele, dominant needs segregation or de-novo status. Kept, not
#: deleted (GP-13), but scored for what it is.
_INHERITANCE_SINGLE_UNKNOWN = 0.20
_INHERITANCE_OTHER = 0.50

# ---------------------------------------------------------------------------
# Contradictions (GP-19). Combined by noisy-OR, so they accumulate without ever
# exceeding 1.0 and no single reason can be double-counted into dominance.
# ---------------------------------------------------------------------------

CONTRADICTION_COMMON_VARIANT = 0.45
CONTRADICTION_POPULATION_HOMOZYGOTES = 0.35
CONTRADICTION_CIS_PHASE = 0.60
CONTRADICTION_LOW_QUALITY = 0.30
CONTRADICTION_BENIGN_CONSEQUENCE = 0.25

#: Component score at or above which an emitted item is recorded as SUPPORTS.
#: A reporting convention, not a significance threshold.
SUPPORT_THRESHOLD = 0.5


@dataclass(frozen=True)
class Contradiction:
    """One recorded reason the candidate argues against itself."""

    code: str
    magnitude: float
    subject_id: str
    claim: str
    category: EvidenceCategory


@dataclass(frozen=True)
class ScoredPair:
    """A candidate with its score vector, its evidence and its unanswered questions."""

    candidate: PairCandidate
    scores: ComponentScores
    composite: float
    supporting_evidence: tuple[EvidenceItem, ...]
    contradicting_evidence: tuple[EvidenceItem, ...]
    open_questions: tuple[OpenQuestion, ...]
    rationale: str

    @property
    def pair_id(self) -> str:
        return self.candidate.pair_id

    @property
    def evidence(self) -> tuple[EvidenceItem, ...]:
        return (*self.supporting_evidence, *self.contradicting_evidence)


def _clamp(value: float) -> float:
    return min(1.0, max(0.0, value))


# ---------------------------------------------------------------------------
# Analytical validity
# ---------------------------------------------------------------------------


def _variant_analytical_validity(
    variant: VariantRecord, quality: QualityThresholds
) -> tuple[float, list[str]]:
    score = 1.0
    notes: list[str] = []
    genotype = variant.genotype

    if variant.filter_status is FilterStatus.FILTERED:
        score *= _FILTERED_MULTIPLIER
        notes.append(f"caller FILTER {'/'.join(variant.raw_filters) or 'filtered'}")
    elif variant.filter_status is FilterStatus.MISSING:
        score *= _FILTER_MISSING_MULTIPLIER
        notes.append("no caller FILTER recorded")

    if genotype.depth is None:
        score *= _UNREPORTED_METRIC_MULTIPLIER
        notes.append("depth not reported")
    elif genotype.depth < quality.min_depth:
        score *= max(_METRIC_FLOOR, genotype.depth / quality.min_depth)
        notes.append(f"depth {genotype.depth} below the {quality.min_depth}x minimum")

    gq = genotype.genotype_quality
    if gq is None:
        score *= _UNREPORTED_METRIC_MULTIPLIER
        notes.append("genotype quality not reported")
    elif gq < quality.min_genotype_quality:
        score *= max(_METRIC_FLOOR, gq / quality.min_genotype_quality)
        notes.append(f"GQ {gq} below the required {quality.min_genotype_quality}")

    balance = genotype.allele_balance
    if genotype.is_heterozygous:
        if balance is None:
            score *= _ALLELE_BALANCE_UNKNOWN_MULTIPLIER
            notes.append("allele balance unknowable from the recorded depths")
        elif not (quality.min_allele_balance_het <= balance <= quality.max_allele_balance_het):
            if balance >= quality.mosaic_allele_balance_floor:
                score *= _MOSAIC_MULTIPLIER
                notes.append(
                    f"allele balance {balance:.2f} outside the het band but within the "
                    "mosaic window, so treated as possible mosaicism rather than noise"
                )
            else:
                score *= _ALLELE_BALANCE_ARTIFACT_MULTIPLIER
                notes.append(f"allele balance {balance:.2f} below the mosaic floor")

    if FLAG_LOW_QUALITY_CALL in variant.qc_flags and score >= 0.9:
        score *= _UPSTREAM_FLAG_MULTIPLIER
        notes.append("flagged low quality upstream for a reason not re-derivable here")

    return _clamp(score), notes


def score_analytical_validity(pair: PairCandidate, quality: QualityThresholds) -> tuple[float, str]:
    """Do we believe the calls are real? Scored on the WEAKEST member call.

    A hypothesis that needs two variants is only as real as the less credible of
    the two; averaging would let an immaculate call launder a four-read one.
    """
    scored = [
        (variant, *_variant_analytical_validity(variant, quality)) for variant in pair.variants
    ]
    weakest = min(scored, key=lambda item: (item[1], item[0].variant_id))
    variant, score, notes = weakest
    if notes:
        detail = f"{variant.variant_id} has {', '.join(notes)}"
    else:
        detail = (
            f"{variant.variant_id} passes FILTER at depth {variant.genotype.depth} and "
            f"GQ {variant.genotype.genotype_quality} with balanced allele support"
        )
    scope = "weakest of the two calls" if pair.is_pair else "the single call"
    return score, (f"Analytical validity {score:.2f}: scored on the {scope}; {detail}.")


# ---------------------------------------------------------------------------
# Rarity
# ---------------------------------------------------------------------------


def _rarity_from_frequency(allele_frequency: float, frequency: FrequencyThresholds) -> float:
    """Piecewise-linear interpolation of the configured cut-points in log10 space."""
    if allele_frequency <= 0.0:
        return _RARITY_AT_ULTRA_RARE
    knots: list[tuple[float, float]] = [
        (math.log10(max(frequency.ultra_rare, 1e-12)), _RARITY_AT_ULTRA_RARE),
        (math.log10(max(frequency.rare, 1e-12)), _RARITY_AT_RARE),
        (math.log10(max(frequency.low_frequency, 1e-12)), _RARITY_AT_LOW_FREQUENCY),
        (math.log10(max(frequency.max_plausible_recessive, 1e-12)), _RARITY_AT_MAX_PLAUSIBLE),
        (math.log10(_RARITY_FLOOR_AF), _RARITY_FLOOR),
    ]
    value = math.log10(allele_frequency)
    if value <= knots[0][0]:
        return _RARITY_AT_ULTRA_RARE
    for (lo_x, lo_y), (hi_x, hi_y) in pairwise(knots):
        if value <= hi_x:
            span = hi_x - lo_x
            if span <= 0:
                return hi_y
            fraction = (value - lo_x) / span
            return _clamp(lo_y + fraction * (hi_y - lo_y))
    return _RARITY_FLOOR


def _variant_rarity(
    variant: VariantRecord, frequency: FrequencyThresholds
) -> tuple[float, PopulationFrequency | None]:
    observed = variant.max_allele_frequency()
    if observed is None:
        return frequency.absent_frequency_score, None
    return _rarity_from_frequency(observed.allele_frequency, frequency), observed


def score_rarity(pair: PairCandidate, frequency: FrequencyThresholds) -> tuple[float, str]:
    """How consistent the observed frequencies are with a severe recessive disorder.

    Scored on the most common member allele, using the MAXIMUM frequency across
    populations rather than the global figure. Absence of frequency data scores
    the configured mid-range value, because a site missing from a reference
    cohort is usually poorly covered there, not rare (GP-14).
    """
    scored = [(variant, *_variant_rarity(variant, frequency)) for variant in pair.variants]
    variant, score, observed = min(scored, key=lambda item: (item[1], item[0].variant_id))
    if observed is None:
        detail = (
            f"{variant.variant_id} has no population-frequency record and is scored at the "
            f"configured neutral {frequency.absent_frequency_score:.2f}; absence of "
            "frequency data is not evidence of rarity (GP-14)"
        )
    else:
        detail = (
            f"the most common member allele is {variant.variant_id} at "
            f"AF {observed.allele_frequency:.3g} ({observed.provenance_key}), against "
            f"cut-points ultra-rare {frequency.ultra_rare:.3g} / rare "
            f"{frequency.rare:.3g} / max-plausible-recessive "
            f"{frequency.max_plausible_recessive:.3g}"
        )
    return score, f"Rarity {score:.2f}: {detail}."


# ---------------------------------------------------------------------------
# Molecular consequence
# ---------------------------------------------------------------------------


def _best_numeric(variant: VariantRecord, gene: str, key: str) -> float | None:
    values = [
        csq.pathogenicity_scores[key]
        for csq in variant.consequences_for_gene(gene)
        if key in csq.pathogenicity_scores
    ]
    return max(values) if values else None


def _best_spliceai(variant: VariantRecord, gene: str) -> float | None:
    values = [
        csq.splice_ai_delta_max
        for csq in variant.consequences_for_gene(gene)
        if csq.splice_ai_delta_max is not None
    ]
    return max(values) if values else None


def _variant_consequence(variant: VariantRecord, gene: str) -> tuple[float, str]:
    impact = variant.worst_impact_for_gene(gene)
    if impact is None:
        return _IMPACT_UNANNOTATED, "no consequence annotation for this gene (unknown, not benign)"

    base = _IMPACT_BASE[impact]
    uplift = 0.0
    parts: list[str] = []

    terms = {term for csq in variant.consequences_for_gene(gene) for term in csq.consequence_terms}
    lof = sorted(terms & LOF_TERMS)
    most_severe = min(
        (csq.most_severe_term for csq in variant.consequences_for_gene(gene)), default="unknown"
    )
    if lof:
        uplift += _LOF_UPLIFT
        parts.append(f"loss-of-function term {lof[0]}")

    spliceai = _best_spliceai(variant, gene)
    if spliceai is not None and spliceai >= _SPLICEAI_FLOOR:
        uplift += _SPLICEAI_UPLIFT * _clamp((spliceai - _SPLICEAI_FLOOR) / (1.0 - _SPLICEAI_FLOOR))
        parts.append(f"SpliceAI {spliceai:.2f}")

    cadd = _best_numeric(variant, gene, "CADD_phred")
    if cadd is not None and cadd >= _CADD_FLOOR:
        uplift += _CADD_UPLIFT * _clamp((cadd - _CADD_FLOOR) / (_CADD_CEILING - _CADD_FLOOR))
        parts.append(f"CADD {cadd:.1f}")

    revel = _best_numeric(variant, gene, "REVEL")
    if revel is not None and revel >= _REVEL_FLOOR:
        uplift += _REVEL_UPLIFT * _clamp((revel - _REVEL_FLOOR) / (1.0 - _REVEL_FLOOR))
        parts.append(f"REVEL {revel:.2f}")

    detail = f"{most_severe} ({impact.value.upper()} impact)"
    if parts:
        detail += ", " + ", ".join(parts)
    return _clamp(base + uplift), detail


def score_consequence(pair: PairCandidate) -> tuple[float, str]:
    """Predicted damage to the gene product, scored on the LESS damaged allele.

    A biallelic hypothesis requires both copies to be affected, so the weaker
    prediction sets the ceiling. In-silico uplifts (SpliceAI, CADD, REVEL) are
    capped contributions on top of the ordinal impact class, never a substitute
    for it: a prediction is not an observation (GP-12).
    """
    scored = [
        (variant, *_variant_consequence(variant, pair.gene_symbol)) for variant in pair.variants
    ]
    variant, score, detail = min(scored, key=lambda item: (item[1], item[0].variant_id))
    scope = "less damaging of the two alleles" if pair.is_pair else "the single allele"
    return score, (
        f"Molecular consequence {score:.2f}: the {scope}, {variant.variant_id}, is {detail}."
    )


# ---------------------------------------------------------------------------
# Inheritance consistency
# ---------------------------------------------------------------------------


def _phase_multiplier(status: PhaseStatus, phase_weights: PhaseWeights) -> float:
    match status:
        case PhaseStatus.TRANS_CONFIRMED:
            return phase_weights.trans_confirmed
        case PhaseStatus.TRANS_LIKELY:
            return phase_weights.trans_likely
        case PhaseStatus.UNKNOWN:
            return phase_weights.unknown
        case PhaseStatus.CIS_LIKELY:
            return phase_weights.cis_likely
        case PhaseStatus.CIS_CONFIRMED:
            return phase_weights.cis_confirmed


def _inheritance_base(pair: PairCandidate) -> tuple[float, str]:
    match pair.inheritance_model:
        case InheritanceModel.COMPOUND_HETEROZYGOUS:
            return _INHERITANCE_COMPOUND_HET, "two heterozygous alleles in one gene"
        case InheritanceModel.HOMOZYGOUS_RECESSIVE | InheritanceModel.X_LINKED_RECESSIVE:
            return _INHERITANCE_HOMOZYGOUS, "a single call accounting for both gene copies"
        case InheritanceModel.UNKNOWN if pair.is_pair:
            return _INHERITANCE_MIXED_PAIR, "two alleles of mixed zygosity, model undetermined"
        case InheritanceModel.UNKNOWN:
            return (
                _INHERITANCE_SINGLE_UNKNOWN,
                "a lone heterozygote with neither a second hit nor de-novo evidence",
            )
        case _:
            return _INHERITANCE_OTHER, f"model {pair.inheritance_model.value}"


def score_inheritance(pair: PairCandidate, phase_weights: PhaseWeights) -> tuple[float, str]:
    """Fit to the inheritance model, including the phase penalty (GP-15).

    The phase multiplier applies only to two-variant candidates, where the
    haplotype relationship is what makes or breaks the hypothesis. It is applied
    exactly as configured — unknown phase is penalised, never quietly promoted.
    """
    base, base_detail = _inheritance_base(pair)
    if not pair.is_pair:
        return _clamp(base), (
            f"Inheritance consistency {_clamp(base):.2f}: {pair.inheritance_model.value} from "
            f"{base_detail}; no phase multiplier applies to a single-variant hypothesis."
        )
    multiplier = _phase_multiplier(pair.phase.status, phase_weights)
    score = _clamp(base * multiplier)
    return score, (
        f"Inheritance consistency {score:.2f}: {pair.inheritance_model.value} scored "
        f"{base:.2f} for {base_detail}, multiplied by the configured phase weight "
        f"{multiplier:.2f} for {pair.phase.status.value} (method "
        f"{pair.phase.method}); phase is never assumed (GP-15)."
    )


# ---------------------------------------------------------------------------
# Evidence quality
# ---------------------------------------------------------------------------


def _has_orthogonal_predictors(variant: VariantRecord, gene: str) -> bool:
    predictors = 0
    if _best_spliceai(variant, gene) is not None:
        predictors += 1
    predictors += len(
        {key for csq in variant.consequences_for_gene(gene) for key in csq.pathogenicity_scores}
    )
    return predictors >= 2


def score_evidence_quality(pair: PairCandidate) -> tuple[float, str]:
    """Strength and INDEPENDENCE of the evidence, not its quantity.

    Six equally weighted facets: versioned frequency provenance, annotation
    depth on a named transcript, orthogonal in-silico predictors, curated
    clinical assertions, whether any phasing information exists at all, and
    whether the hypothesis rests on more than one independently called site.
    """
    variants = pair.variants
    gene = pair.gene_symbol
    facets: dict[str, float] = {
        "frequency_provenance": sum(v.has_frequency_data for v in variants) / len(variants),
        "annotation_depth": sum(
            1.0
            if any(csq.is_mane_select or csq.is_canonical for csq in v.consequences_for_gene(gene))
            else (0.5 if v.consequences_for_gene(gene) else 0.0)
            for v in variants
        )
        / len(variants),
        "orthogonal_predictors": sum(_has_orthogonal_predictors(v, gene) for v in variants)
        / len(variants),
        "clinical_curation": sum(bool(v.clinical_assertions) for v in variants) / len(variants),
        "phase_information": 1.0 if pair.phase.method != "none" else 0.0,
        "independent_observations": 1.0 if pair.is_pair else 0.5,
    }
    score = _clamp(sum(facets.values()) / len(facets))
    weak = sorted(name for name, value in facets.items() if value == 0.0)
    detail = f"absent: {', '.join(weak)}" if weak else "all facets present"
    return score, f"Evidence quality {score:.2f}: mean of six independence facets ({detail})."


# ---------------------------------------------------------------------------
# Contradictions (GP-19)
# ---------------------------------------------------------------------------


def collect_contradictions(
    pair: PairCandidate, frequency: FrequencyThresholds, quality: QualityThresholds
) -> tuple[Contradiction, ...]:
    """Enumerate everything that argues AGAINST this candidate.

    Nothing here is ever discarded: each entry becomes a persisted
    contradicting ``EvidenceItem`` as well as feeding the penalty.
    """
    found: list[Contradiction] = []

    if pair.inheritance_model is InheritanceModel.COMPOUND_HETEROZYGOUS and (
        pair.phase.status in CIS_STATES
    ):
        found.append(
            Contradiction(
                code="in_cis_phase",
                magnitude=CONTRADICTION_CIS_PHASE,
                subject_id=pair.pair_id,
                claim=(
                    f"Phase evidence places both alleles of {pair.gene_symbol} on the same "
                    f"haplotype ({pair.phase.status.value}), leaving one gene copy intact; "
                    "a compound heterozygote is not possible under this call."
                ),
                category=EvidenceCategory.INHERITANCE,
            )
        )

    for variant in sorted(pair.variants, key=lambda v: v.variant_id):
        observed = variant.max_allele_frequency()
        if FLAG_COMMON_VARIANT in variant.qc_flags and observed is not None:
            found.append(
                Contradiction(
                    code="common_variant",
                    magnitude=CONTRADICTION_COMMON_VARIANT,
                    subject_id=variant.variant_id,
                    claim=(
                        f"{variant.variant_id} reaches AF {observed.allele_frequency:.3g} in "
                        f"{observed.provenance_key}, above the maximum plausible "
                        f"{frequency.max_plausible_recessive:.3g} for a severe recessive allele."
                    ),
                    category=EvidenceCategory.POPULATION,
                )
            )
        if (
            observed is not None
            and observed.homozygote_count is not None
            and observed.homozygote_count > 0
            and observed.allele_frequency > frequency.low_frequency
        ):
            found.append(
                Contradiction(
                    code="population_homozygotes",
                    magnitude=CONTRADICTION_POPULATION_HOMOZYGOTES,
                    subject_id=variant.variant_id,
                    claim=(
                        f"{observed.homozygote_count} homozygotes for {variant.variant_id} are "
                        f"recorded in {observed.provenance_key}; an unselected cohort carrying "
                        "homozygotes argues against a severe early-onset recessive effect."
                    ),
                    category=EvidenceCategory.POPULATION,
                )
            )
        if FLAG_LOW_QUALITY_CALL in variant.qc_flags:
            found.append(
                Contradiction(
                    code="low_quality_call",
                    magnitude=CONTRADICTION_LOW_QUALITY,
                    subject_id=variant.variant_id,
                    claim=(
                        f"{variant.variant_id} breached a call-quality threshold "
                        f"(FILTER {variant.filter_status.value}, depth "
                        f"{variant.genotype.depth}, GQ {variant.genotype.genotype_quality}; "
                        f"minimums {quality.min_depth}x / {quality.min_genotype_quality}), so "
                        "the genotype itself may not be real."
                    ),
                    category=EvidenceCategory.ANALYTICAL,
                )
            )
        if FLAG_BENIGN_CONSEQUENCE in variant.qc_flags:
            found.append(
                Contradiction(
                    code="benign_consequence",
                    magnitude=CONTRADICTION_BENIGN_CONSEQUENCE,
                    subject_id=variant.variant_id,
                    claim=(
                        f"Every predicted consequence of {variant.variant_id} is MODIFIER or "
                        "LOW impact, so no mechanism of protein disruption is predicted."
                    ),
                    category=EvidenceCategory.CONSEQUENCE,
                )
            )
    return tuple(found)


def contradiction_penalty(contradictions: tuple[Contradiction, ...]) -> float:
    """Combine contradiction magnitudes by noisy-OR.

    Independent reasons accumulate but can never sum past 1.0, and a single
    reason can never be counted twice into dominance.
    """
    remaining = 1.0
    for item in contradictions:
        remaining *= 1.0 - _clamp(item.magnitude)
    return _clamp(1.0 - remaining)


# ---------------------------------------------------------------------------
# Composite
# ---------------------------------------------------------------------------


def composite_score(scores: ComponentScores, weights: ScoringWeights) -> float:
    """Weighted sum of the seven positive components minus the penalty, clamped.

    Reported alongside the component vector, never instead of it: the number is
    a ranking device, and the components are what a reviewer actually argues
    with.
    """
    positive = (
        weights.analytical_validity * scores.analytical_validity
        + weights.rarity * scores.rarity
        + weights.molecular_consequence * scores.molecular_consequence
        + weights.inheritance_consistency * scores.inheritance_consistency
        + weights.phenotype_similarity * scores.phenotype_similarity
        + weights.mechanistic_relevance * scores.mechanistic_relevance
        + weights.evidence_quality * scores.evidence_quality
    )
    penalty = weights.contradiction_penalty_weight * scores.contradiction_penalty
    return _clamp(positive - penalty)


# ---------------------------------------------------------------------------
# Evidence emission
# ---------------------------------------------------------------------------


def _strength_for(score: float) -> EvidenceStrength:
    """Distance from neutral, capped at MODERATE.

    A pipeline inference is never STRONG or DEFINITIVE however extreme the
    number: the ceiling encodes that this module predicts, it does not measure.
    """
    distance = abs(score - SUPPORT_THRESHOLD)
    if distance >= 0.40:
        return EvidenceStrength.MODERATE
    if distance >= 0.25:
        return EvidenceStrength.SUPPORTING
    if distance >= 0.10:
        return EvidenceStrength.WEAK
    return EvidenceStrength.INSUFFICIENT


def _component_evidence(
    *,
    subject_id: str,
    category: EvidenceCategory,
    claim: str,
    score: float,
    method: str,
    limitation: str,
    clock: Clock,
    payload: dict[str, str | int | float | bool | None],
) -> EvidenceItem:
    # A low component score is WEAK SUPPORT, not opposition. Emitting CONTRADICTS
    # here would conflate "this component scored poorly" with "there is evidence
    # against this hypothesis", which floods the dossier's contradiction section
    # with noise and devalues GP-19. Genuine contradictions come from
    # `collect_contradictions` (in-cis phase, common allele frequency, population
    # homozygotes) and are the only things that feed the subtracted penalty.
    direction = (
        EvidenceDirection.SUPPORTS if score >= SUPPORT_THRESHOLD else EvidenceDirection.NEUTRAL
    )
    return EvidenceItem(
        evidence_id=make_evidence_id(
            subject_id=subject_id, category=category, claim=claim, tool=TOOL_NAME
        ),
        subject_id=subject_id,
        subject_kind="pair",
        claim=claim,
        category=category,
        direction=direction,
        strength=_strength_for(score),
        evidence_type=EvidenceType.PIPELINE_INFERENCE,
        tier=AssertionTier.INFERENCE,
        method=method,
        tool=TOOL_NAME,
        tool_version=TOOL_VERSION,
        limitations=f"{HEURISTIC_LIMITATION} {limitation}",
        timestamp=clock.now(),
        numeric_value=score,
        payload=payload,
    )


def _contradiction_evidence(
    contradiction: Contradiction, pair_id: str, clock: Clock
) -> EvidenceItem:
    """Build a contradicting evidence item.

    The payload carries ``pair_id`` **only when the contradiction is genuinely
    about the pair** (e.g. in-cis phase). A variant-level fact — "this allele is
    common", "the cohort contains homozygotes" — is true independently of which
    pair we happened to be considering when we noticed it, and the same variant
    participates in several pairs. Stamping the pair into a variant-level claim
    produced the same content-derived evidence ID with differing payloads, which
    the ledger correctly rejected as a collision. The pair -> evidence link is
    already carried by ``CandidatePair.contradicting_evidence_ids``; it does not
    need to be duplicated into the evidence itself.
    """
    is_pair_level = contradiction.subject_id == pair_id
    payload: dict[str, str | int | float | bool | None] = {"contradiction_code": contradiction.code}
    if is_pair_level:
        payload["pair_id"] = pair_id

    return EvidenceItem(
        evidence_id=make_evidence_id(
            subject_id=contradiction.subject_id,
            category=contradiction.category,
            claim=contradiction.claim,
            tool=TOOL_NAME,
        ),
        subject_id=contradiction.subject_id,
        subject_kind="pair" if is_pair_level else "variant",
        claim=contradiction.claim,
        category=contradiction.category,
        direction=EvidenceDirection.CONTRADICTS,
        strength=_strength_for(SUPPORT_THRESHOLD + contradiction.magnitude),
        evidence_type=EvidenceType.PIPELINE_INFERENCE,
        tier=AssertionTier.INFERENCE,
        method=(
            "Rule-based contradiction detection over the annotated variant record; "
            "magnitudes combine by noisy-OR into the subtracted contradiction penalty."
        ),
        tool=TOOL_NAME,
        tool_version=TOOL_VERSION,
        limitations=(
            f"{HEURISTIC_LIMITATION} A recorded contradiction lowers the rank; it does not "
            "remove the candidate, because the underlying annotation may itself be wrong "
            "(GP-13, GP-19)."
        ),
        timestamp=clock.now(),
        numeric_value=contradiction.magnitude,
        payload=payload,
    )


# ---------------------------------------------------------------------------
# Open questions
# ---------------------------------------------------------------------------


def _open_questions(
    pair: PairCandidate, contradictions: tuple[Contradiction, ...]
) -> tuple[OpenQuestion, ...]:
    """Enumerate what is still unknown. Absence of a question is itself a claim."""
    questions: list[OpenQuestion] = []
    codes = {item.code for item in contradictions}

    if pair.is_pair and pair.phase.status is PhaseStatus.UNKNOWN:
        questions.append(
            OpenQuestion(
                question_id=f"{pair.pair_id}::phase-trans",
                question="Is the pair in trans?",
                why_it_matters=(
                    "Two heterozygous variants in one gene only constitute a compound "
                    "heterozygote on opposite haplotypes. In cis, one gene copy is intact "
                    "and the recessive hypothesis fails outright."
                ),
                resolving_test="parental segregation testing or long-read/read-backed phasing",
                blocking=True,
            )
        )
    if "in_cis_phase" in codes:
        questions.append(
            OpenQuestion(
                question_id=f"{pair.pair_id}::phase-cis-confirmation",
                question="Is the in-cis phase call correct?",
                why_it_matters=(
                    "Short-read phase sets can be wrong. The candidate is heavily "
                    "down-ranked rather than deleted, so the cis call itself needs "
                    "orthogonal confirmation before the hypothesis is discarded."
                ),
                resolving_test="long-read sequencing or parental genotypes across both sites",
                blocking=False,
            )
        )
    for variant in sorted(pair.variants, key=lambda v: v.variant_id):
        if FLAG_NO_FREQUENCY_DATA in variant.qc_flags:
            questions.append(
                OpenQuestion(
                    question_id=f"{pair.pair_id}::frequency::{variant.variant_id}",
                    question=f"What is the population allele frequency of {variant.variant_id}?",
                    why_it_matters=(
                        "No frequency record exists, so the variant was scored at a neutral "
                        "mid-range value. Absence of data is not rarity (GP-14) and a common "
                        "allele cannot yet be excluded."
                    ),
                    resolving_test=(
                        "look the site up in a versioned population reference with adequate "
                        "coverage at this locus, ancestry-matched where possible"
                    ),
                    blocking=True,
                )
            )
        if FLAG_LOW_QUALITY_CALL in variant.qc_flags:
            questions.append(
                OpenQuestion(
                    question_id=f"{pair.pair_id}::call-validity::{variant.variant_id}",
                    question=f"Is the genotype call at {variant.variant_id} real?",
                    why_it_matters=(
                        "The call breached a quality threshold. No downstream conclusion is "
                        "safe while the genotype itself is in doubt."
                    ),
                    resolving_test="orthogonal confirmation by Sanger or targeted resequencing",
                    blocking=True,
                )
            )
    if not pair.is_pair and pair.inheritance_model is InheritanceModel.UNKNOWN:
        questions.append(
            OpenQuestion(
                question_id=f"{pair.pair_id}::second-allele",
                question="Is there a second hit, or is this variant de novo?",
                why_it_matters=(
                    "A lone heterozygote supports neither a recessive nor a dominant "
                    "conclusion without either a second allele or de-novo status."
                ),
                resolving_test=(
                    "trio sequencing for de-novo status, plus CNV/structural analysis of the "
                    "second allele in case a deletion is masked"
                ),
                blocking=True,
            )
        )
    if any(v.genotype.zygosity is Zygosity.HOM_ALT for v in pair.variants):
        questions.append(
            OpenQuestion(
                question_id=f"{pair.pair_id}::true-homozygosity",
                question="Is the homozygous call true homozygosity rather than hemizygosity?",
                why_it_matters=(
                    "A deletion on one allele makes a heterozygote look homozygous, which "
                    "changes both the inheritance model and the recurrence risk."
                ),
                resolving_test="copy-number analysis and runs-of-homozygosity assessment",
                blocking=False,
            )
        )
    if not any(variant.clinical_assertions for variant in pair.variants):
        questions.append(
            OpenQuestion(
                question_id=f"{pair.pair_id}::clinical-curation",
                question="Has either variant been curated by a clinical database?",
                why_it_matters=(
                    "Every score here is a pipeline inference; no independent curated "
                    "assertion is currently supporting or opposing the candidate."
                ),
                resolving_test="ClinVar/literature review of both loci at a pinned release",
                blocking=False,
            )
        )
    return tuple(questions)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def score_pair(
    pair: PairCandidate,
    *,
    phenotype_score: float,
    mechanism_score: float,
    weights: ScoringWeights,
    phase_weights: PhaseWeights,
    frequency: FrequencyThresholds,
    quality: QualityThresholds,
    clock: Clock,
) -> ScoredPair:
    """Produce the full score vector, its evidence and its unanswered questions.

    ``phenotype_score`` and ``mechanism_score`` are computed by the phenotype and
    mechanism stages and passed in (GP-03). Callers with nothing to supply pass
    :data:`NEUTRAL_PHENOTYPE_SCORE` / :data:`NEUTRAL_MECHANISM_SCORE`, and the
    rationale then states plainly that no comparison was made.
    """
    validity, validity_note = score_analytical_validity(pair, quality)
    rarity, rarity_note = score_rarity(pair, frequency)
    consequence, consequence_note = score_consequence(pair)
    inheritance, inheritance_note = score_inheritance(pair, phase_weights)
    evidence_quality, evidence_note = score_evidence_quality(pair)

    phenotype = _clamp(phenotype_score)
    mechanism = _clamp(mechanism_score)
    phenotype_note = f"Phenotype similarity {phenotype:.2f}: supplied by the phenotype stage" + (
        " as the documented neutral default, meaning no gene-phenotype comparison was "
        "available for this candidate rather than that the comparison failed."
        if phenotype == NEUTRAL_PHENOTYPE_SCORE
        else "."
    )
    mechanism_note = f"Mechanistic relevance {mechanism:.2f}: supplied by the mechanism stage" + (
        " as the documented neutral default, meaning no mechanism chain was available "
        "for this gene rather than that none exists."
        if mechanism == NEUTRAL_MECHANISM_SCORE
        else "."
    )

    contradictions = collect_contradictions(pair, frequency, quality)
    penalty = contradiction_penalty(contradictions)

    scores = ComponentScores(
        analytical_validity=validity,
        rarity=rarity,
        molecular_consequence=consequence,
        inheritance_consistency=inheritance,
        phenotype_similarity=phenotype,
        mechanistic_relevance=mechanism,
        evidence_quality=evidence_quality,
        contradiction_penalty=penalty,
    )
    composite = composite_score(scores, weights)

    items = _build_component_evidence(
        pair,
        clock=clock,
        notes={
            "analytical_validity": (validity, validity_note),
            "rarity": (rarity, rarity_note),
            "molecular_consequence": (consequence, consequence_note),
            "inheritance_consistency": (inheritance, inheritance_note),
            "phenotype_similarity": (phenotype, phenotype_note),
            "mechanistic_relevance": (mechanism, mechanism_note),
            "evidence_quality": (evidence_quality, evidence_note),
        },
    )
    contradiction_items = tuple(
        _contradiction_evidence(item, pair.pair_id, clock) for item in contradictions
    )
    supporting = tuple(item for item in items if item.direction is EvidenceDirection.SUPPORTS)
    contradicting = (
        *(item for item in items if item.direction is EvidenceDirection.CONTRADICTS),
        *contradiction_items,
    )

    penalty_note = (
        f"Contradiction penalty {penalty:.2f} from "
        + "; ".join(item.claim for item in contradictions)
        if contradictions
        else "Contradiction penalty 0.00: nothing recorded argues against this candidate."
    )
    questions = _open_questions(pair, contradictions)
    blocking = sum(1 for question in questions if question.blocking)
    rationale = " ".join(
        [
            validity_note,
            rarity_note,
            consequence_note,
            inheritance_note,
            phenotype_note,
            mechanism_note,
            evidence_note,
            penalty_note,
            f"Composite {composite:.3f} from the configured weights; "
            f"{len(questions)} open question(s), {blocking} blocking.",
        ]
    )

    return ScoredPair(
        candidate=pair,
        scores=scores,
        composite=composite,
        supporting_evidence=supporting,
        contradicting_evidence=contradicting,
        open_questions=questions,
        rationale=rationale,
    )


_COMPONENT_CATEGORIES: dict[str, EvidenceCategory] = {
    "analytical_validity": EvidenceCategory.ANALYTICAL,
    "rarity": EvidenceCategory.POPULATION,
    "molecular_consequence": EvidenceCategory.CONSEQUENCE,
    "inheritance_consistency": EvidenceCategory.INHERITANCE,
    "phenotype_similarity": EvidenceCategory.PHENOTYPE,
    "mechanistic_relevance": EvidenceCategory.MECHANISM,
    "evidence_quality": EvidenceCategory.PROVENANCE,
}

_COMPONENT_METHODS: dict[str, str] = {
    "analytical_validity": (
        "Multiplicative penalties on caller FILTER, depth, genotype quality and allele "
        "balance, evaluated on the weakest member call."
    ),
    "rarity": (
        "Piecewise-linear interpolation of log10(maximum population allele frequency) "
        "between the configured ultra-rare/rare/low-frequency/max-plausible cut-points."
    ),
    "molecular_consequence": (
        "Ordinal VEP-style impact class as the base, with capped uplifts for "
        "loss-of-function terms, SpliceAI, CADD and REVEL, taken on the less damaged allele."
    ),
    "inheritance_consistency": (
        "Model-fit base multiplied by the configured phase weight for the inferred phase status."
    ),
    "phenotype_similarity": "Supplied by the phenotype stage; not computed here (GP-03).",
    "mechanistic_relevance": "Supplied by the mechanism stage; not computed here (GP-03).",
    "evidence_quality": (
        "Unweighted mean of six independence facets over the candidate's member variants."
    ),
}

_COMPONENT_LIMITATIONS: dict[str, str] = {
    "analytical_validity": (
        "Depth and GQ are proxies for call correctness, not confirmation; only orthogonal "
        "resequencing establishes that a genotype is real."
    ),
    "rarity": (
        "Reference cohorts under-represent many ancestries, so a low observed frequency may "
        "reflect sampling rather than true rarity, and absence of a record establishes nothing."
    ),
    "molecular_consequence": (
        "In-silico impact and pathogenicity predictors are correlated with one another and "
        "with training-set ascertainment; agreement between them is not independent support."
    ),
    "inheritance_consistency": (
        "Derived from a caller phase set and zygosity alone. It establishes no segregation, "
        "no de-novo status and no parental genotype."
    ),
    "phenotype_similarity": (
        "This module did not compute the value and cannot vouch for how the phenotype stage "
        "handled unassessed terms."
    ),
    "mechanistic_relevance": (
        "This module did not compute the value; mechanism links are literature-derived and "
        "may not hold in the relevant tissue or developmental window."
    ),
    "evidence_quality": (
        "Counts the presence of evidence facets, not their correctness. A well-provenanced "
        "wrong annotation scores exactly as highly as a right one."
    ),
}


def _build_component_evidence(
    pair: PairCandidate, *, clock: Clock, notes: dict[str, tuple[float, str]]
) -> tuple[EvidenceItem, ...]:
    return tuple(
        _component_evidence(
            subject_id=pair.pair_id,
            category=_COMPONENT_CATEGORIES[component],
            claim=claim,
            score=score,
            method=_COMPONENT_METHODS[component],
            limitation=_COMPONENT_LIMITATIONS[component],
            clock=clock,
            payload={
                "component": component,
                "gene_symbol": pair.gene_symbol,
                "inheritance_model": pair.inheritance_model.value,
                "phase_status": pair.phase.status.value,
            },
        )
        for component, (score, claim) in notes.items()
    )
