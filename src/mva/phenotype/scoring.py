"""Phenotype-similarity scoring under four-valued observation logic.

This module exists to prevent one specific, common and very expensive bug: a
phenotype scorer that treats "we never assessed this" as "the patient does not
have this". That collapse turns a four-valued clinical logic into a boolean, and
it fails in the direction that matters — it manufactures evidence against exactly
the candidate genes whose associated features are the ones nobody thought to
check. The gene with the richest curated phenotype profile is the gene most
likely to have unassessed terms, so a naive denominator systematically penalises
the best candidates.

The four statuses are kept distinct end to end (GP-14):

===============  ==============================  ==========================
Patient status   Contribution to the score       Recorded as
===============  ==============================  ==========================
``OBSERVED``     positive, weighted by strength  ``matched_terms``
``EXCLUDED``     negative, weighted by strength  ``contradicted_terms``
``NOT_ASSESSED`` **exactly zero, both ways**     ``unassessed_terms``
``UNCERTAIN``    **exactly zero, both ways**     rationale only
===============  ==============================  ==========================

``EXCLUDED`` is the only status permitted to carry negative weight — the same
rule :data:`mva.models.phenotype.NEGATIVE_EVIDENCE_STATUSES` states as a type.

**Maturity (GP-20): synthetic-substitute.** The scoring logic is real and tested,
but it is driven by ``knowledge/public/gene_phenotype.tsv``, whose gene symbols
are fictional. It also computes term-level overlap only: there is no ontology
traversal, so "Microcephaly" and its parent "Abnormality of head size" are
unrelated strings to this module. A real deployment needs an HPO release with
ancestor closure and an information-content measure (Resnik/Phenomizer). Until
then no output here may be described as biologically valid.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Final

from mva import __version__
from mva.clock import Clock
from mva.models.base import AssertionTier
from mva.models.evidence import (
    Citation,
    EvidenceCategory,
    EvidenceDirection,
    EvidenceItem,
    EvidenceStrength,
    EvidenceType,
    make_evidence_id,
)
from mva.models.phenotype import ObservationStatus, PhenotypeProfile
from mva.phenotype.hpo import GeneAssociation, GenePhenotypeIndex

TOOL_NAME: Final[str] = "mva.phenotype.scoring"
TOOL_VERSION: Final[str] = __version__

#: Score returned when the gene has **no informative overlap** with the profile.
#:
#: Deliberately mid-range, and deliberately not 0.0. A gene about which this
#: pipeline holds no assessed phenotype information — because it has no curated
#: terms, or because every one of its terms is ``NOT_ASSESSED``/``UNCERTAIN`` in
#: this patient — has not been argued against. Returning 0.0 would let ignorance
#: masquerade as a negative finding and would push such genes below genes that
#: were actually contradicted by an assessment, which inverts the evidential
#: ordering.
#:
#: This mirrors :attr:`mva.config.FrequencyThresholds.absent_frequency_score`,
#: which takes the same position for missing population frequency, and is the
#: numeric expression of GP-14.
NEUTRAL_SCORE: Final[float] = 0.5

#: Weight of the *specificity* component (do the gene's assessed terms agree?).
#: Heuristic, not calibrated (GP-32); changing it needs a decision record.
SPECIFICITY_WEIGHT: Final[float] = 0.5

#: Weight of the *coverage* component (how much of the patient does the gene explain?).
COVERAGE_WEIGHT: Final[float] = 0.5

#: Decimal places the final score is rounded to, so repeat runs and cross-platform
#: float formatting cannot produce differing artifact bytes (GP-30).
SCORE_PRECISION: Final[int] = 6

_STRENGTH_TO_EVIDENCE: Final[dict[str, EvidenceStrength]] = {
    "definitive": EvidenceStrength.DEFINITIVE,
    "strong": EvidenceStrength.STRONG,
    "moderate": EvidenceStrength.MODERATE,
    "supporting": EvidenceStrength.SUPPORTING,
}

_CURATION_LIMITATION: Final[str] = (
    "HPO gene-phenotype association strength is a human curation judgement, not a "
    "measurement: it reflects how much has been published about a gene, which "
    "correlates with study effort as much as with biology. Term overlap here is "
    "exact-match only, with no ontology ancestor closure, so a parent/child term "
    "pair scores as unrelated."
)

_NOT_ASSESSED_LIMITATION: Final[str] = (
    "Absence of assessment is not absence of the feature. This item records an "
    "information gap and carries zero weight in both directions; it must never be "
    "read as evidence that the patient lacks the term, and it does not enter the "
    "score's denominator."
)


@dataclass(frozen=True, slots=True)
class PhenotypeScoreBreakdown:
    """The arithmetic behind a :class:`PhenotypeMatch`, exposed for auditing."""

    score: float
    specificity: float | None
    """``matched_weight / informative_weight``; ``None`` when nothing was informative."""
    coverage: float | None
    """Fraction of the patient's OBSERVED terms this gene accounts for."""
    matched_weight: float
    contradicted_weight: float
    informative_weight: float
    """The denominator of ``specificity``. Built from informative terms only."""
    observed_term_count: int


@dataclass(frozen=True, slots=True)
class PhenotypeMatch:
    """Result of scoring one gene against one patient's phenotype profile."""

    gene_symbol: str
    score: float
    """In ``[0, 1]``. See :func:`score_gene_phenotype` for the denominator rationale."""
    matched_terms: tuple[str, ...]
    """OBSERVED in the patient **and** associated with the gene."""
    contradicted_terms: tuple[str, ...]
    """EXCLUDED in the patient but associated with the gene — real negative evidence."""
    unassessed_terms: tuple[str, ...]
    """Associated with the gene, NOT_ASSESSED in the patient. Zero weight, both ways."""
    unexplained_terms: tuple[str, ...]
    """OBSERVED in the patient but not associated with the gene."""
    rationale: str
    evidence: tuple[EvidenceItem, ...]

    @property
    def has_contradiction(self) -> bool:
        return bool(self.contradicted_terms)


@dataclass(frozen=True, slots=True)
class _Partition:
    """Gene-associated terms split by the patient's four-valued status."""

    matched: tuple[GeneAssociation, ...]
    contradicted: tuple[GeneAssociation, ...]
    unassessed: tuple[GeneAssociation, ...]
    uncertain: tuple[GeneAssociation, ...]
    unexplained: tuple[str, ...]


def score_gene_phenotype(
    gene_symbol: str,
    *,
    profile: PhenotypeProfile,
    index: GenePhenotypeIndex,
    clock: Clock,
) -> PhenotypeMatch:
    """Score one gene against one phenotype profile.

    **The denominator, and why it is what it is.**

    The score is a weighted mean of two ratios, each with an explicitly chosen
    denominator::

        specificity = matched_weight / (matched_weight + contradicted_weight)
        coverage    = |matched terms| / |patient OBSERVED terms|
        score       = SPECIFICITY_WEIGHT * specificity + COVERAGE_WEIGHT * coverage

    *Specificity* asks: of the gene's curated terms that this patient was actually
    assessed for, how many came back positive? Its denominator — named
    ``informative_weight`` — is the summed association weight of the gene's
    ``OBSERVED`` and ``EXCLUDED`` terms and **nothing else**. ``NOT_ASSESSED`` and
    ``UNCERTAIN`` terms are absent from both numerator and denominator.

    The tempting alternative is to divide by the gene's *total* associated weight.
    That is the bug this module exists to prevent. Under it, a gene curated with
    twenty features would be capped at 4/20 = 0.2 when the clinic assessed four of
    them and all four matched — while a gene curated with exactly those four
    features scores 1.0 on identical patient data. The penalty tracks how
    thoroughly the gene has been *studied* and how much of the workup happened to
    be *ordered*, not how well it fits the patient. Excluding uninformative terms
    from the denominator makes the guarantee exact and mechanically testable: a
    gene associated with a ``NOT_ASSESSED`` term scores identically to the same
    gene without that association at all. ``tests/unit/test_phenotype.py`` asserts
    that equality directly.

    *Coverage* supplies the other half, which specificity alone cannot see: how
    much of the patient's presentation the gene accounts for. Without it, a gene
    matching one of four observed features scores the same as a gene matching all
    four. Its denominator is the patient's ``OBSERVED`` terms — again informative
    only, since ``NOT_ASSESSED`` and ``UNCERTAIN`` patient terms are excluded, so
    an unassessed term cannot shrink or inflate this ratio either. Coverage is
    count-based rather than weight-based because a patient's observations carry no
    curation strength; only the gene side does.

    **The neutral case.** If the gene has no informative overlap at all — no
    ``OBSERVED`` and no ``EXCLUDED`` associated term — both ratios are undefined
    and the function returns :data:`NEUTRAL_SCORE` (0.5), not 0.0. No information
    is not evidence of a poor match. A contradicted gene can therefore score
    *below* a gene nobody has phenotype data for, which is the correct ordering:
    one has been argued against, the other has not been argued about.

    Every score contribution is accompanied by an :class:`EvidenceItem` (GP-10)
    with a non-empty ``limitations`` field (GP-17). All outputs are sorted (GP-30).
    """
    associations = index.terms_for_gene(gene_symbol)
    partition = _partition(associations, profile=profile)
    breakdown = _compute_breakdown(partition, profile=profile)
    rationale = _build_rationale(gene_symbol, partition=partition, breakdown=breakdown)
    evidence = _build_evidence(
        gene_symbol,
        partition=partition,
        breakdown=breakdown,
        profile=profile,
        index=index,
        rationale=rationale,
        clock=clock,
    )
    return PhenotypeMatch(
        gene_symbol=gene_symbol,
        score=breakdown.score,
        matched_terms=tuple(assoc.hpo_id for assoc in partition.matched),
        contradicted_terms=tuple(assoc.hpo_id for assoc in partition.contradicted),
        unassessed_terms=tuple(assoc.hpo_id for assoc in partition.unassessed),
        unexplained_terms=partition.unexplained,
        rationale=rationale,
        evidence=evidence,
    )


def score_all_genes(
    gene_symbols: Sequence[str],
    *,
    profile: PhenotypeProfile,
    index: GenePhenotypeIndex,
    clock: Clock,
) -> dict[str, PhenotypeMatch]:
    """Score many genes, returning a dict in sorted key order.

    Duplicate symbols collapse to one entry. Insertion order is sorted rather than
    input order so that anything serialising this mapping — a JSON artifact, a
    Parquet table, a report table — is byte-stable across runs (GP-30).
    """
    ordered = sorted(set(gene_symbols), key=lambda name: (name.upper(), name))
    return {
        symbol: score_gene_phenotype(symbol, profile=profile, index=index, clock=clock)
        for symbol in ordered
    }


# ---------------------------------------------------------------------------
# Partitioning and arithmetic
# ---------------------------------------------------------------------------


def _partition(associations: Sequence[GeneAssociation], *, profile: PhenotypeProfile) -> _Partition:
    """Split the gene's associated terms by the patient's status for each.

    ``profile.status_of`` returns ``NOT_ASSESSED`` for a term the profile never
    mentions, which is the safe default: a term nobody recorded is a term nobody
    assessed. That is what keeps "the gene is associated with a feature this
    workup did not cover" out of the negative bucket.
    """
    buckets: dict[ObservationStatus, list[GeneAssociation]] = {
        status: [] for status in ObservationStatus
    }
    for assoc in associations:
        buckets[profile.status_of(assoc.hpo_id)].append(assoc)

    associated = {assoc.hpo_id for assoc in associations}
    unexplained = tuple(sorted(term for term in profile.observed_terms if term not in associated))

    return _Partition(
        matched=tuple(sorted(buckets[ObservationStatus.OBSERVED], key=lambda a: a.hpo_id)),
        contradicted=tuple(sorted(buckets[ObservationStatus.EXCLUDED], key=lambda a: a.hpo_id)),
        unassessed=tuple(sorted(buckets[ObservationStatus.NOT_ASSESSED], key=lambda a: a.hpo_id)),
        uncertain=tuple(sorted(buckets[ObservationStatus.UNCERTAIN], key=lambda a: a.hpo_id)),
        unexplained=unexplained,
    )


def _compute_breakdown(
    partition: _Partition, *, profile: PhenotypeProfile
) -> PhenotypeScoreBreakdown:
    """Apply the documented formula. See :func:`score_gene_phenotype`."""
    matched_weight = sum(assoc.weight for assoc in partition.matched)
    contradicted_weight = sum(assoc.weight for assoc in partition.contradicted)
    informative_weight = matched_weight + contradicted_weight
    observed_count = len(profile.observed_terms)

    if informative_weight <= 0.0:
        # No assessed, gene-associated term exists. Nothing has been established
        # for or against this gene; return the documented neutral value.
        return PhenotypeScoreBreakdown(
            score=NEUTRAL_SCORE,
            specificity=None,
            coverage=None,
            matched_weight=0.0,
            contradicted_weight=0.0,
            informative_weight=0.0,
            observed_term_count=observed_count,
        )

    specificity = matched_weight / informative_weight
    coverage = (len(partition.matched) / observed_count) if observed_count else NEUTRAL_SCORE
    raw = SPECIFICITY_WEIGHT * specificity + COVERAGE_WEIGHT * coverage
    score = round(min(1.0, max(0.0, raw)), SCORE_PRECISION)

    return PhenotypeScoreBreakdown(
        score=score,
        specificity=round(specificity, SCORE_PRECISION),
        coverage=round(coverage, SCORE_PRECISION),
        matched_weight=round(matched_weight, SCORE_PRECISION),
        contradicted_weight=round(contradicted_weight, SCORE_PRECISION),
        informative_weight=round(informative_weight, SCORE_PRECISION),
        observed_term_count=observed_count,
    )


def _aggregate_strength(breakdown: PhenotypeScoreBreakdown) -> EvidenceStrength:
    """Band the computed score into an ordinal strength.

    Never ``DEFINITIVE``: this is the pipeline's own inference over a curated
    substitute knowledge base, and no arrangement of term overlap makes that a
    definitive statement about a patient.
    """
    if breakdown.specificity is None:
        return EvidenceStrength.INSUFFICIENT
    if breakdown.score >= 0.85:
        return EvidenceStrength.STRONG
    if breakdown.score >= 0.65:
        return EvidenceStrength.MODERATE
    if breakdown.score >= 0.40:
        return EvidenceStrength.SUPPORTING
    return EvidenceStrength.WEAK


# ---------------------------------------------------------------------------
# Rationale
# ---------------------------------------------------------------------------


def _terms(associations: Sequence[GeneAssociation]) -> str:
    return ", ".join(assoc.hpo_id for assoc in associations) if associations else "none"


def _build_rationale(
    gene_symbol: str, *, partition: _Partition, breakdown: PhenotypeScoreBreakdown
) -> str:
    """One deterministic paragraph explaining the number, including what was ignored."""
    parts = [
        f"{gene_symbol}: matched (observed and gene-associated) {_terms(partition.matched)}; "
        f"contradicted (explicitly EXCLUDED after assessment) {_terms(partition.contradicted)}; "
        f"not assessed (no information, excluded from the denominator per GP-14) "
        f"{_terms(partition.unassessed)}; uncertain (assessed but equivocal, also "
        f"contributing zero) {_terms(partition.uncertain)}; observed but not associated "
        f"with this gene "
        f"{', '.join(partition.unexplained) if partition.unexplained else 'none'}."
    ]
    if breakdown.specificity is None:
        parts.append(
            f"No gene-associated term was informative in this profile, so the score is the "
            f"neutral {NEUTRAL_SCORE:.3f} rather than 0: absence of information is not "
            "evidence of a poor match."
        )
    else:
        parts.append(
            f"Specificity {breakdown.specificity:.3f} "
            f"(matched weight {breakdown.matched_weight:.2f} of informative weight "
            f"{breakdown.informative_weight:.2f}); coverage "
            f"{breakdown.coverage if breakdown.coverage is not None else NEUTRAL_SCORE:.3f} "
            f"({len(partition.matched)} of {breakdown.observed_term_count} observed terms "
            f"explained). Score {breakdown.score:.3f} = "
            f"{SPECIFICITY_WEIGHT:g}*specificity + {COVERAGE_WEIGHT:g}*coverage."
        )
    return " ".join(parts)


# ---------------------------------------------------------------------------
# Evidence (GP-10 / GP-17)
# ---------------------------------------------------------------------------


def _make_item(
    *,
    gene_symbol: str,
    claim: str,
    direction: EvidenceDirection,
    strength: EvidenceStrength,
    evidence_type: EvidenceType,
    tier: AssertionTier,
    method: str,
    limitations: str,
    citation: Citation,
    clock: Clock,
    numeric_value: float | None = None,
    payload: dict[str, str | int | float | bool | None] | None = None,
) -> EvidenceItem:
    """Build one PHENOTYPE evidence item with a content-derived, stable id."""
    return EvidenceItem(
        evidence_id=make_evidence_id(
            subject_id=gene_symbol,
            category=EvidenceCategory.PHENOTYPE,
            claim=claim,
            tool=TOOL_NAME,
        ),
        subject_id=gene_symbol,
        subject_kind="gene",
        claim=claim,
        category=EvidenceCategory.PHENOTYPE,
        direction=direction,
        strength=strength,
        evidence_type=evidence_type,
        tier=tier,
        citation=citation,
        method=method,
        tool=TOOL_NAME,
        tool_version=TOOL_VERSION,
        limitations=limitations,
        timestamp=clock.now(),
        run_id=None,
        numeric_value=numeric_value,
        payload=payload or {},
    )


def _build_evidence(
    gene_symbol: str,
    *,
    partition: _Partition,
    breakdown: PhenotypeScoreBreakdown,
    profile: PhenotypeProfile,
    index: GenePhenotypeIndex,
    rationale: str,
    clock: Clock,
) -> tuple[EvidenceItem, ...]:
    """Assemble the evidence ledger for one gene, in a fixed deterministic order.

    Order is: the computed match, then the per-term observations grouped by status
    (matched, contradicted, unassessed, uncertain), then the unexplained-features
    summary. Within each group terms are iterated in sorted HPO-id order, so the
    sequence is fully determined by the inputs (GP-30).
    """
    profile_citation = Citation(
        source="phenotype_profile",
        identifier=profile.source_artifact,
        version=profile.hpo_version,
        url=None,
        title=None,
    )
    index_citation = Citation(
        source="gene_phenotype_index",
        identifier=gene_symbol,
        version=index.version,
        url=None,
        title=None,
    )

    items: list[EvidenceItem] = [
        _computed_match_item(
            gene_symbol,
            breakdown=breakdown,
            partition=partition,
            rationale=rationale,
            source_artifact=profile.source_artifact,
            citation=index_citation,
            clock=clock,
        )
    ]
    items.extend(
        _observation_items(gene_symbol, partition=partition, citation=profile_citation, clock=clock)
    )
    items.extend(
        _information_gap_items(
            gene_symbol, partition=partition, citation=profile_citation, clock=clock
        )
    )
    unexplained_item = _unexplained_item(
        gene_symbol, partition=partition, citation=profile_citation, clock=clock
    )
    if unexplained_item is not None:
        items.append(unexplained_item)
    return tuple(items)


def _computed_match_item(
    gene_symbol: str,
    *,
    breakdown: PhenotypeScoreBreakdown,
    partition: _Partition,
    rationale: str,
    source_artifact: str,
    citation: Citation,
    clock: Clock,
) -> EvidenceItem:
    """The INFERENCE-tier item carrying the score itself."""
    if breakdown.score > NEUTRAL_SCORE:
        direction = EvidenceDirection.SUPPORTS
    elif breakdown.score < NEUTRAL_SCORE:
        direction = EvidenceDirection.CONTRADICTS
    else:
        direction = EvidenceDirection.NEUTRAL

    claim = (
        f"Phenotype similarity between {gene_symbol} and the profile in {source_artifact} "
        f"scores {breakdown.score:.3f} over {len(partition.matched)} matched, "
        f"{len(partition.contradicted)} contradicted and {len(partition.unassessed)} "
        f"unassessed gene-associated HPO terms."
    )
    return _make_item(
        gene_symbol=gene_symbol,
        claim=claim,
        direction=direction,
        strength=_aggregate_strength(breakdown),
        evidence_type=EvidenceType.PIPELINE_INFERENCE,
        tier=AssertionTier.INFERENCE,
        method=(
            f"score = {SPECIFICITY_WEIGHT:g}*specificity + {COVERAGE_WEIGHT:g}*coverage, where "
            "specificity = matched association weight / informative (OBSERVED + EXCLUDED) "
            "association weight, and coverage = matched terms / observed patient terms. "
            "NOT_ASSESSED and UNCERTAIN terms are excluded from both numerator and "
            f"denominator. Zero informative overlap returns the neutral {NEUTRAL_SCORE}. "
            f"Rationale: {rationale}"
        ),
        limitations=(
            f"{_CURATION_LIMITATION} The component weights ({SPECIFICITY_WEIGHT:g}/"
            f"{COVERAGE_WEIGHT:g}) are heuristics chosen for the demo case, not values "
            "calibrated against a labelled cohort, so the number ranks candidates and does "
            "not estimate a probability. The gene-phenotype knowledge base backing this "
            "score is a synthetic substitute (GP-20)."
        ),
        citation=citation,
        clock=clock,
        numeric_value=breakdown.score,
        payload={
            "specificity": breakdown.specificity,
            "coverage": breakdown.coverage,
            "matched_weight": breakdown.matched_weight,
            "contradicted_weight": breakdown.contradicted_weight,
            "informative_weight": breakdown.informative_weight,
            "matched_count": len(partition.matched),
            "contradicted_count": len(partition.contradicted),
            "unassessed_count": len(partition.unassessed),
            "uncertain_count": len(partition.uncertain),
            "unexplained_count": len(partition.unexplained),
            "observed_term_count": breakdown.observed_term_count,
        },
    )


def _observation_items(
    gene_symbol: str,
    *,
    partition: _Partition,
    citation: Citation,
    clock: Clock,
) -> list[EvidenceItem]:
    """OBSERVED_DATA items for terms the clinician actually assessed."""
    items: list[EvidenceItem] = []
    for assoc in partition.matched:
        items.append(
            _make_item(
                gene_symbol=gene_symbol,
                claim=(
                    f"{assoc.hpo_id} ({assoc.label}) is OBSERVED in this subject and is a "
                    f"{assoc.association_strength} phenotype association of {gene_symbol}."
                ),
                direction=EvidenceDirection.SUPPORTS,
                strength=_STRENGTH_TO_EVIDENCE[assoc.association_strength],
                evidence_type=EvidenceType.DIRECT_MEASUREMENT,
                tier=AssertionTier.OBSERVED_DATA,
                method=(
                    "Exact HPO term intersection between the subject's OBSERVED terms and "
                    f"the gene's curated associations ({assoc.source} {assoc.version})."
                ),
                limitations=(
                    f"{_CURATION_LIMITATION} A shared term shows the gene can produce this "
                    "feature, not that it did so in this subject; common features such as "
                    "microcephaly are associated with hundreds of genes."
                ),
                citation=citation,
                clock=clock,
                numeric_value=assoc.weight,
                payload={
                    "hpo_id": assoc.hpo_id,
                    "status": ObservationStatus.OBSERVED.value,
                    "association_strength": assoc.association_strength,
                },
            )
        )
    for assoc in partition.contradicted:
        items.append(
            _make_item(
                gene_symbol=gene_symbol,
                claim=(
                    f"{assoc.hpo_id} ({assoc.label}) was assessed and EXCLUDED in this "
                    f"subject, yet is a {assoc.association_strength} phenotype association "
                    f"of {gene_symbol}; this argues against {gene_symbol}."
                ),
                direction=EvidenceDirection.CONTRADICTS,
                strength=_STRENGTH_TO_EVIDENCE[assoc.association_strength],
                evidence_type=EvidenceType.DIRECT_MEASUREMENT,
                tier=AssertionTier.OBSERVED_DATA,
                method=(
                    "Exact HPO term intersection between the subject's EXCLUDED terms and "
                    f"the gene's curated associations ({assoc.source} {assoc.version}). "
                    "EXCLUDED is the only status permitted to carry negative weight."
                ),
                limitations=(
                    f"{_CURATION_LIMITATION} A curated association is not fully penetrant "
                    "and may be age-dependent: a feature genuinely absent today can appear "
                    "later, so this contradiction down-ranks the gene and must not be used "
                    "to delete it (GP-13, GP-19)."
                ),
                citation=citation,
                clock=clock,
                numeric_value=assoc.weight,
                payload={
                    "hpo_id": assoc.hpo_id,
                    "status": ObservationStatus.EXCLUDED.value,
                    "association_strength": assoc.association_strength,
                },
            )
        )
    return items


def _information_gap_items(
    gene_symbol: str,
    *,
    partition: _Partition,
    citation: Citation,
    clock: Clock,
) -> list[EvidenceItem]:
    """NEUTRAL items recording what was never assessed, or assessed inconclusively.

    These carry zero weight by construction — they are logged so that a reader can
    see *which* questions were never asked, which is usually the most actionable
    output of a phenotype match: it names the next clinical test to order.
    """
    items: list[EvidenceItem] = []
    for assoc in partition.unassessed:
        items.append(
            _make_item(
                gene_symbol=gene_symbol,
                claim=(
                    f"{assoc.hpo_id} ({assoc.label}) is a {assoc.association_strength} "
                    f"phenotype association of {gene_symbol} but was NOT ASSESSED in this "
                    "subject; it contributes zero to the score in either direction."
                ),
                direction=EvidenceDirection.NEUTRAL,
                strength=EvidenceStrength.INSUFFICIENT,
                evidence_type=EvidenceType.PIPELINE_INFERENCE,
                tier=AssertionTier.INFERENCE,
                method=(
                    "Term present in the gene's curated associations and absent from the "
                    "subject's assessed terms, or explicitly recorded as not_assessed. "
                    "Excluded from both the numerator and the denominator of the score."
                ),
                limitations=(
                    f"{_NOT_ASSESSED_LIMITATION} It is an actionable gap rather than a "
                    "finding: assessing this feature could move the score in either "
                    "direction."
                ),
                citation=citation,
                clock=clock,
                numeric_value=0.0,
                payload={
                    "hpo_id": assoc.hpo_id,
                    "status": ObservationStatus.NOT_ASSESSED.value,
                    "association_strength": assoc.association_strength,
                },
            )
        )
    for assoc in partition.uncertain:
        items.append(
            _make_item(
                gene_symbol=gene_symbol,
                claim=(
                    f"{assoc.hpo_id} ({assoc.label}) is a {assoc.association_strength} "
                    f"phenotype association of {gene_symbol} but was recorded as UNCERTAIN "
                    "in this subject; it contributes zero to the score in either direction."
                ),
                direction=EvidenceDirection.NEUTRAL,
                strength=EvidenceStrength.INSUFFICIENT,
                evidence_type=EvidenceType.PIPELINE_INFERENCE,
                tier=AssertionTier.INFERENCE,
                method=(
                    "Term present in the gene's curated associations and recorded with "
                    "status 'uncertain'. Excluded from both the numerator and the "
                    "denominator of the score."
                ),
                limitations=(
                    "An equivocal assessment is not a negative one. Treating it as either "
                    "presence or absence would fabricate certainty the record does not "
                    "contain; re-assessment, not re-weighting, is the remedy."
                ),
                citation=citation,
                clock=clock,
                numeric_value=0.0,
                payload={
                    "hpo_id": assoc.hpo_id,
                    "status": ObservationStatus.UNCERTAIN.value,
                    "association_strength": assoc.association_strength,
                },
            )
        )
    return items


def _unexplained_item(
    gene_symbol: str,
    *,
    partition: _Partition,
    citation: Citation,
    clock: Clock,
) -> EvidenceItem | None:
    """Summary item for observed features the gene does not account for."""
    if not partition.unexplained:
        return None
    return _make_item(
        gene_symbol=gene_symbol,
        claim=(
            f"{len(partition.unexplained)} OBSERVED feature(s) are not curated associations "
            f"of {gene_symbol}: {', '.join(partition.unexplained)}."
        ),
        direction=EvidenceDirection.NEUTRAL,
        strength=EvidenceStrength.INSUFFICIENT,
        evidence_type=EvidenceType.PIPELINE_INFERENCE,
        tier=AssertionTier.INFERENCE,
        method=(
            "Set difference between the subject's OBSERVED terms and the gene's curated "
            "associations. Lowers the coverage component of the score."
        ),
        limitations=(
            "Recorded as NEUTRAL, not as a contradiction: an unexplained feature may come "
            "from a second condition, from prematurity or from a treatment effect, and the "
            "curated association list for any gene is incomplete. It lowers coverage, and "
            "must not be reported as evidence that the gene is wrong."
        ),
        citation=citation,
        clock=clock,
        numeric_value=float(len(partition.unexplained)),
        payload={"unexplained_count": len(partition.unexplained)},
    )
