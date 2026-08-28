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

Two modes, and which one ran is always recorded
-----------------------------------------------

:class:`ScoringMode.ONTOLOGY`
    The real measure. A :class:`~mva.phenotype.semantics.PhenotypeSemantics` —
    one parsed HPO release plus an information-content model built from the
    published annotation corpus — is supplied by the composition root. Term
    comparison is graph-aware in both directions (see
    :mod:`mva.phenotype.propagation`) and graded by published semantic similarity
    (see :mod:`mva.phenotype.similarity`).

:class:`ScoringMode.EXACT`
    The degraded fallback used when no release is wired in. Terms match only if
    their identifiers are equal, so a parent/child pair scores as unrelated. Its
    arithmetic is unchanged from before ontology support existed, deliberately:
    the fallback must not quietly move a score that a golden expectation is
    locked to. Every evidence item it emits says, in its ``limitations``, that no
    ontology was loaded.

Both modes share the same shape::

    score = SPECIFICITY_WEIGHT * specificity + COVERAGE_WEIGHT * coverage

*Specificity* asks: of the gene's curated terms that this patient was actually
assessed for, how many came back positive? *Coverage* asks: how much of the
patient's presentation does the gene account for? What changes between modes is
how "assessed for" and "accounts for" are decided — by string equality, or by the
ontology.

**Maturity (GP-20).** In ``ONTOLOGY`` mode the ontology, the annotation corpus and
the similarity measure are all real published resources, and the release string is
recorded in provenance. The gene→term associations in
``knowledge/public/gene_phenotype.tsv`` remain a **synthetic substitute** with
fictional gene symbols, so no output of this module is biologically valid until a
real gene-phenotype knowledge base replaces it. The two are separate maturity
claims and must not be conflated.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
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
from mva.phenotype.hpo import (
    UNCURATED_ASSOCIATION_WEIGHT,
    GeneAssociation,
    GenePhenotypeIndex,
)
from mva.phenotype.propagation import InferenceBasis, PropagatedProfile
from mva.phenotype.semantics import PhenotypeSemantics
from mva.phenotype.similarity import AGGREGATION_CITATION, BestMatch

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


class ScoringMode(StrEnum):
    """How term comparison was performed. Always recorded on the result.

    Never inferred by a reader from the presence of an ontology argument: the two
    modes produce different numbers from the same inputs, so which one ran is part
    of the result, not part of the call site.
    """

    EXACT = "exact_term_match"
    """Identifier equality only. No ontology loaded; parent/child pairs do not match."""

    ONTOLOGY = "ontology_semantic"
    """DAG-aware entailment plus information-content semantic similarity."""


class GeneAnnotationStatus(StrEnum):
    """Whether this pipeline holds any phenotype knowledge about the gene at all.

    This distinction is GP-14 applied to the knowledge base rather than to the
    patient. "We have no curated features for this gene" and "we have curated
    features and none of them fit" are different facts that a bare number cannot
    separate, and a ranking that treats the first as the second buries every gene
    that is simply under-studied.
    """

    ANNOTATED = "annotated"
    """The knowledge base holds at least one curated term for this gene."""

    NO_ANNOTATIONS = "no_annotations"
    """The knowledge base holds nothing. The score is the neutral value and carries
    no evidential content in either direction."""


_STRENGTH_TO_EVIDENCE: Final[dict[str, EvidenceStrength]] = {
    "definitive": EvidenceStrength.DEFINITIVE,
    "strong": EvidenceStrength.STRONG,
    "moderate": EvidenceStrength.MODERATE,
    "supporting": EvidenceStrength.SUPPORTING,
    "limited": EvidenceStrength.WEAK,
    "disputed": EvidenceStrength.INSUFFICIENT,
    "refuted": EvidenceStrength.INSUFFICIENT,
    "no known disease relationship": EvidenceStrength.INSUFFICIENT,
}

#: What an association carries when NO source classifies the gene's validity.
#:
#: ``INSUFFICIENT`` is the honest reading of "nobody has curated this": the evidence
#: for the gene-disease relationship is not weak, it is unstated. It is deliberately
#: not ``SUPPORTING``, which would assert a curation nobody made (ADR 0021).
_UNCURATED_EVIDENCE_STRENGTH: Final[EvidenceStrength] = EvidenceStrength.INSUFFICIENT


def _evidence_strength(assoc: GeneAssociation) -> EvidenceStrength:
    """Evidence strength for one association, absence included (ADR 0021)."""
    if assoc.association_strength is None:
        return _UNCURATED_EVIDENCE_STRENGTH
    return _STRENGTH_TO_EVIDENCE[assoc.association_strength]


def _strength_phrase(assoc: GeneAssociation) -> str:
    """``"a definitive"`` / ``"an uncurated"`` — the article travels with the word."""
    if assoc.association_strength is None:
        return "an uncurated (no source classifies this gene's disease validity)"
    return f"a {assoc.association_strength}"


_CURATION_LIMITATION: Final[str] = (
    "Gene-phenotype association strength is a human curation judgement, not a "
    "measurement: it reflects how much has been published about a gene, which "
    "correlates with study effort as much as with biology. It is also a GENE-level "
    "clinical-validity classification, so it does not vary between the terms of one "
    "gene and cannot discriminate within it (ADR 0021)."
)

_EXACT_MODE_LIMITATION: Final[str] = (
    "No HPO release was supplied to this run, so term overlap is exact-match only: "
    "a parent/child term pair scores as unrelated and a true match through the "
    "ontology is invisible (tech debt TD-04). Supply a PhenotypeSemantics to the "
    "scorer to enable graph-aware comparison."
)

_ONTOLOGY_MODE_LIMITATION: Final[str] = (
    "Semantic similarity is computed over the ontology graph and an annotation "
    "corpus, so it measures how close two curated vocabularies are, not whether "
    "this gene caused this patient's presentation. Information content reflects "
    "how often a term is annotated in the published corpus, which tracks curation "
    "effort as well as true rarity, and the corpus is incomplete for recently "
    "described conditions."
)

_NOT_ASSESSED_LIMITATION: Final[str] = (
    "Absence of assessment is not absence of the feature. This item records an "
    "information gap and carries zero weight in both directions; it must never be "
    "read as evidence that the patient lacks the term, and it does not enter the "
    "score's denominator."
)

_NO_ANNOTATION_LIMITATION: Final[str] = (
    "This is a statement about the knowledge base, not about the gene. A gene with "
    "no curated phenotype associations has not been argued against; it scores the "
    "neutral value and must not be ranked below a gene that was actively "
    "contradicted by an assessment (GP-14)."
)


@dataclass(frozen=True, slots=True)
class PhenotypeScoreBreakdown:
    """The arithmetic behind a :class:`PhenotypeMatch`, exposed for auditing."""

    score: float
    specificity: float | None
    """``matched_weight / informative_weight``; ``None`` when nothing was informative."""
    coverage: float | None
    """How much of the patient's OBSERVED presentation the gene accounts for.

    In ``EXACT`` mode this is the fraction of observed terms the gene lists
    literally. In ``ONTOLOGY`` mode it is the best-match average of Köhler et al.
    (2009) over the patient's observed terms — a graded value, so a gene annotated
    with the parent of an observed term scores well rather than zero."""
    matched_weight: float
    contradicted_weight: float
    informative_weight: float
    """The denominator of ``specificity``. Built from informative terms only."""
    observed_term_count: int
    mode: ScoringMode
    annotation_status: GeneAnnotationStatus
    symmetric_similarity: float | None
    """The Phenomizer symmetric best-match average, for audit only. It is **not**
    the score: its patient-ward direction averages over every curated term of the
    gene, including ones nobody assessed, which is the GP-14 penalty this module
    refuses to apply. Reported so a reader can see both numbers. ``None`` outside
    ``ONTOLOGY`` mode."""
    ic_missing_term_count: int
    """Gene terms with no information content in the corpus. They keep their full
    curated weight rather than being dropped; recorded so the reader knows the
    IC-weighting was not applied to them."""
    uncomparable_observed_count: int
    """Observed patient terms the corpus could not place. Excluded from the
    coverage mean rather than averaged in as zero (GP-14)."""


@dataclass(frozen=True, slots=True)
class PhenotypeMatch:
    """Result of scoring one gene against one patient's phenotype profile."""

    gene_symbol: str
    score: float
    """In ``[0, 1]``. See :func:`score_gene_phenotype` for the denominator rationale."""
    matched_terms: tuple[str, ...]
    """OBSERVED in the patient (directly or by upward entailment) **and** associated
    with the gene."""
    contradicted_terms: tuple[str, ...]
    """EXCLUDED in the patient (directly or by downward entailment) but associated
    with the gene — real negative evidence."""
    unassessed_terms: tuple[str, ...]
    """Associated with the gene, NOT_ASSESSED in the patient. Zero weight, both ways."""
    unexplained_terms: tuple[str, ...]
    """OBSERVED in the patient and not accounted for by the gene."""
    rationale: str
    evidence: tuple[EvidenceItem, ...]
    breakdown: PhenotypeScoreBreakdown
    mode: ScoringMode
    annotation_status: GeneAnnotationStatus
    implied_matched_terms: tuple[str, ...]
    """Gene terms that matched **only** through the ontology — a descendant of the
    term was observed. Separated from direct matches because an entailed match is
    weaker evidence than a recorded one and a report must be able to say so."""
    implied_contradicted_terms: tuple[str, ...]
    """Gene terms contradicted **only** through the ontology — an ancestor of the
    term was assessed and excluded."""
    conflicting_terms: tuple[str, ...]
    """Terms the patient record entails as both present and absent. Non-empty means
    the source record is internally inconsistent; reported, never silently repaired."""
    unresolved_profile_terms: tuple[str, ...]
    """Patient terms absent from the loaded ontology release. An information gap."""
    best_matches: tuple[BestMatch, ...]
    """Per observed term, the gene term that matched it best and their most
    informative common ancestor. Empty outside ``ONTOLOGY`` mode."""

    @property
    def has_contradiction(self) -> bool:
        return bool(self.contradicted_terms)

    @property
    def is_ontology_aware(self) -> bool:
        return self.mode is ScoringMode.ONTOLOGY


@dataclass(frozen=True, slots=True)
class _Partition:
    """Gene-associated terms split by the patient's four-valued status."""

    matched: tuple[GeneAssociation, ...]
    contradicted: tuple[GeneAssociation, ...]
    unassessed: tuple[GeneAssociation, ...]
    uncertain: tuple[GeneAssociation, ...]
    unexplained: tuple[str, ...]
    implied_matched: tuple[str, ...]
    implied_contradicted: tuple[str, ...]


def score_gene_phenotype(
    gene_symbol: str,
    *,
    profile: PhenotypeProfile,
    index: GenePhenotypeIndex,
    clock: Clock,
    semantics: PhenotypeSemantics | None = None,
    propagated: PropagatedProfile | None = None,
) -> PhenotypeMatch:
    """Score one gene against one phenotype profile.

    ``semantics`` is the switch between the two modes documented on this module.
    It is an argument rather than a module-level lookup because the ontology
    release is a resource the composition root owns: a hidden global would make
    "which HPO release produced this number" invisible at the call site and would
    let one case's ontology survive into the next.

    ``propagated`` is an optional pre-built entailment closure for ``profile``,
    passed by :func:`score_all_genes` so the closures are built once per subject
    rather than once per gene. It must have been built from ``semantics``'
    ontology.

    **The denominator, and why it is what it is.**

    The score is a weighted mean of two ratios, each with an explicitly chosen
    denominator::

        specificity = matched_weight / (matched_weight + contradicted_weight)
        coverage    = how much of the patient's presentation the gene accounts for
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
    that equality directly, in both modes.

    In ``ONTOLOGY`` mode each term's contribution is additionally scaled by its
    normalised information content, so a match on "Phenotypic abnormality" counts
    for almost nothing and a match on "Premature chromatid separation" counts for
    almost everything. That is the point of an IC-based measure: specificity is a
    property of the corpus, never of graph depth.

    *Coverage* supplies the other half, which specificity alone cannot see: how
    much of the patient's presentation the gene accounts for. Without it, a gene
    matching one of four observed features scores the same as a gene matching all
    four. Its denominator is the patient's ``OBSERVED`` terms — again informative
    only, since ``NOT_ASSESSED`` and ``UNCERTAIN`` patient terms are excluded, so
    an unassessed term cannot shrink or inflate this ratio either.

    In ``ONTOLOGY`` mode coverage is the best-match average of Köhler et al.
    (2009) in the patient→gene direction. The reverse direction of the published
    symmetric measure is computed and reported for audit but is deliberately not
    scored on: averaging over *every* curated term of the gene re-introduces the
    "penalise the well-studied gene for questions nobody asked" bug that the
    specificity denominator above exists to remove. The
    assessment-restricted specificity component carries the gene→patient
    information instead.

    **The neutral case.** A component that is undefined contributes
    :data:`NEUTRAL_SCORE`, not zero. If both are undefined the score is exactly
    ``NEUTRAL_SCORE``. No information is not evidence of a poor match, and a
    contradicted gene can therefore score *below* a gene nobody has phenotype data
    for — which is the correct ordering: one has been argued against, the other
    has not been argued about.

    **A gene with no curated terms is not a gene that failed to match.** It is
    reported as :attr:`GeneAnnotationStatus.NO_ANNOTATIONS`, scores exactly
    ``NEUTRAL_SCORE``, and carries an evidence item saying so. A gene *with*
    curated terms that fit nothing scores below neutral on the coverage component,
    because that is a real finding about a real annotation set.

    Every score contribution is accompanied by an :class:`EvidenceItem` (GP-10)
    with a non-empty ``limitations`` field (GP-17). All outputs are sorted (GP-30).
    """
    associations = index.terms_for_gene(gene_symbol)
    annotation_status = (
        GeneAnnotationStatus.ANNOTATED if associations else GeneAnnotationStatus.NO_ANNOTATIONS
    )

    if semantics is None:
        partition = _partition_exact(associations, profile=profile)
        breakdown = _compute_breakdown_exact(
            partition, profile=profile, annotation_status=annotation_status
        )
        best_matches: tuple[BestMatch, ...] = ()
        conflicts: tuple[str, ...] = ()
        unresolved: tuple[str, ...] = ()
    else:
        closure = propagated if propagated is not None else semantics.propagate(profile)
        partition = _partition_propagated(associations, propagated=closure)
        breakdown, best_matches = _compute_breakdown_ontology(
            partition,
            associations=associations,
            propagated=closure,
            semantics=semantics,
            annotation_status=annotation_status,
        )
        conflicts = closure.conflicts
        unresolved = closure.unresolved_terms

    rationale = _build_rationale(gene_symbol, partition=partition, breakdown=breakdown)
    evidence = _build_evidence(
        gene_symbol,
        partition=partition,
        breakdown=breakdown,
        profile=profile,
        index=index,
        rationale=rationale,
        clock=clock,
        semantics=semantics,
        conflicts=conflicts,
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
        breakdown=breakdown,
        mode=breakdown.mode,
        annotation_status=annotation_status,
        implied_matched_terms=partition.implied_matched,
        implied_contradicted_terms=partition.implied_contradicted,
        conflicting_terms=conflicts,
        unresolved_profile_terms=unresolved,
        best_matches=best_matches,
    )


def score_all_genes(
    gene_symbols: Sequence[str],
    *,
    profile: PhenotypeProfile,
    index: GenePhenotypeIndex,
    clock: Clock,
    semantics: PhenotypeSemantics | None = None,
) -> dict[str, PhenotypeMatch]:
    """Score many genes, returning a dict in sorted key order.

    Duplicate symbols collapse to one entry. Insertion order is sorted rather than
    input order so that anything serialising this mapping — a JSON artifact, a
    Parquet table, a report table — is byte-stable across runs (GP-30).

    The subject's entailment closure is built once here and reused for every gene.
    Building it per gene would repeat the ontology traversal for every candidate
    and would be the difference between milliseconds and minutes on a real panel.
    """
    ordered = sorted(set(gene_symbols), key=lambda name: (name.upper(), name))
    closure = semantics.propagate(profile) if semantics is not None else None
    return {
        symbol: score_gene_phenotype(
            symbol,
            profile=profile,
            index=index,
            clock=clock,
            semantics=semantics,
            propagated=closure,
        )
        for symbol in ordered
    }


# ---------------------------------------------------------------------------
# Partitioning
# ---------------------------------------------------------------------------


def _partition_exact(
    associations: Sequence[GeneAssociation], *, profile: PhenotypeProfile
) -> _Partition:
    """Split the gene's terms by the patient's status, by identifier equality only.

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
        implied_matched=(),
        implied_contradicted=(),
    )


def _partition_propagated(
    associations: Sequence[GeneAssociation], *, propagated: PropagatedProfile
) -> _Partition:
    """Split the gene's terms by the patient's **entailed** four-valued status.

    The two entailments run in opposite directions and are computed by
    :class:`~mva.phenotype.propagation.PropagatedProfile`: a gene term is matched
    when it subsumes an observed feature (upward), and contradicted when it is
    subsumed by an excluded one (downward). Reversing either is the bug that makes
    every score in a run wrong, so the direction lives in one tested place rather
    than being re-derived here.

    An observed patient term counts as *accounted for* when some gene term stands
    in a subsumption relation to it in either direction — the gene lists the
    general form of what was seen, or a more specific form of it. That is a graph
    fact rather than a similarity threshold, so ``unexplained_terms`` stays
    reviewable and does not move when a weight is retuned.
    """
    ontology = propagated.ontology
    buckets: dict[ObservationStatus, list[GeneAssociation]] = {
        status: [] for status in ObservationStatus
    }
    implied_matched: list[str] = []
    implied_contradicted: list[str] = []

    for assoc in associations:
        inferred = propagated.status_of(assoc.hpo_id)
        buckets[inferred.status].append(assoc)
        if inferred.basis is InferenceBasis.ANCESTOR_OF_OBSERVED:
            implied_matched.append(assoc.hpo_id)
        elif inferred.basis is InferenceBasis.DESCENDANT_OF_EXCLUDED:
            implied_contradicted.append(assoc.hpo_id)

    gene_terms = sorted({assoc.hpo_id for assoc in associations})
    unexplained: list[str] = []
    for observed in propagated.recorded_observed:
        observed_ancestors = frozenset(ontology.ancestor_closure(observed))
        accounted = any(
            term in observed_ancestors or observed in ontology.ancestor_closure(term)
            for term in gene_terms
        )
        if not accounted:
            unexplained.append(observed)

    return _Partition(
        matched=tuple(sorted(buckets[ObservationStatus.OBSERVED], key=lambda a: a.hpo_id)),
        contradicted=tuple(sorted(buckets[ObservationStatus.EXCLUDED], key=lambda a: a.hpo_id)),
        unassessed=tuple(sorted(buckets[ObservationStatus.NOT_ASSESSED], key=lambda a: a.hpo_id)),
        uncertain=tuple(sorted(buckets[ObservationStatus.UNCERTAIN], key=lambda a: a.hpo_id)),
        unexplained=tuple(sorted(unexplained)),
        implied_matched=tuple(sorted(implied_matched)),
        implied_contradicted=tuple(sorted(implied_contradicted)),
    )


# ---------------------------------------------------------------------------
# Arithmetic
# ---------------------------------------------------------------------------


def _compute_breakdown_exact(
    partition: _Partition,
    *,
    profile: PhenotypeProfile,
    annotation_status: GeneAnnotationStatus,
) -> PhenotypeScoreBreakdown:
    """The pre-ontology formula, preserved bit-for-bit.

    Not refactored into the ontology path. Its numbers are what the golden
    expectation in ``tests/golden/`` is locked to (GP-32), and a shared code path
    that "should" produce the same values is exactly how a silent re-baseline
    happens.
    """
    # UNCURATED_ASSOCIATION_WEIGHT is named here, not defaulted inside the property:
    # what an unclassified gene-disease relationship contributes is a scientific
    # choice and stays visible at the call site (ADR 0021).
    matched_weight = sum(
        assoc.weight_or(UNCURATED_ASSOCIATION_WEIGHT) for assoc in partition.matched
    )
    contradicted_weight = sum(
        assoc.weight_or(UNCURATED_ASSOCIATION_WEIGHT) for assoc in partition.contradicted
    )
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
            mode=ScoringMode.EXACT,
            annotation_status=annotation_status,
            symmetric_similarity=None,
            ic_missing_term_count=0,
            uncomparable_observed_count=0,
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
        mode=ScoringMode.EXACT,
        annotation_status=annotation_status,
        symmetric_similarity=None,
        ic_missing_term_count=0,
        uncomparable_observed_count=0,
    )


def _compute_breakdown_ontology(
    partition: _Partition,
    *,
    associations: Sequence[GeneAssociation],
    propagated: PropagatedProfile,
    semantics: PhenotypeSemantics,
    annotation_status: GeneAnnotationStatus,
) -> tuple[PhenotypeScoreBreakdown, tuple[BestMatch, ...]]:
    """The graph-aware formula. See :func:`score_gene_phenotype` for the reasoning."""
    observed = propagated.recorded_observed
    gene_terms = sorted({assoc.hpo_id for assoc in associations})

    if annotation_status is GeneAnnotationStatus.NO_ANNOTATIONS:
        # Not "nothing matched" — nothing was known. Both components are undefined
        # and the score is exactly neutral (GP-14).
        return (
            PhenotypeScoreBreakdown(
                score=NEUTRAL_SCORE,
                specificity=None,
                coverage=None,
                matched_weight=0.0,
                contradicted_weight=0.0,
                informative_weight=0.0,
                observed_term_count=len(observed),
                mode=ScoringMode.ONTOLOGY,
                annotation_status=annotation_status,
                symmetric_similarity=None,
                ic_missing_term_count=0,
                uncomparable_observed_count=0,
            ),
            (),
        )

    information_content = semantics.information_content
    ic_missing = 0

    def scaled(assoc: GeneAssociation) -> float:
        """Curated weight scaled by the term's normalised information content.

        A term with no corpus support keeps its full curated weight rather than
        being scaled to zero: no annotation data is not evidence of generality, and
        silently zeroing it would erase a real curated association (GP-14). The
        count is reported so the reader knows the scaling did not apply.
        """
        nonlocal ic_missing
        # See _compute_breakdown_exact: the uncurated weight is named, never defaulted.
        curated = assoc.weight_or(UNCURATED_ASSOCIATION_WEIGHT)
        scale = information_content.normalised_ic(assoc.hpo_id)
        if scale is None:
            ic_missing += 1
            return curated
        return curated * scale

    matched_weight = sum(scaled(assoc) for assoc in partition.matched)
    contradicted_weight = sum(scaled(assoc) for assoc in partition.contradicted)
    informative_weight = matched_weight + contradicted_weight
    specificity = matched_weight / informative_weight if informative_weight > 0.0 else None

    bma = semantics.similarity.best_match_average(observed, gene_terms)
    coverage = bma.value
    symmetric = semantics.similarity.symmetric_best_match_average(observed, gene_terms)

    raw = SPECIFICITY_WEIGHT * _or_neutral(specificity) + COVERAGE_WEIGHT * _or_neutral(coverage)
    score = round(min(1.0, max(0.0, raw)), SCORE_PRECISION)

    return (
        PhenotypeScoreBreakdown(
            score=score,
            specificity=_round_optional(specificity),
            coverage=_round_optional(coverage),
            matched_weight=round(matched_weight, SCORE_PRECISION),
            contradicted_weight=round(contradicted_weight, SCORE_PRECISION),
            informative_weight=round(informative_weight, SCORE_PRECISION),
            observed_term_count=len(observed),
            mode=ScoringMode.ONTOLOGY,
            annotation_status=annotation_status,
            symmetric_similarity=_round_optional(symmetric),
            ic_missing_term_count=ic_missing,
            uncomparable_observed_count=len(bma.uncomparable),
        ),
        bma.matches,
    )


def _or_neutral(value: float | None) -> float:
    """An undefined component contributes the neutral value, never zero (GP-14)."""
    return NEUTRAL_SCORE if value is None else value


def _round_optional(value: float | None) -> float | None:
    return None if value is None else round(value, SCORE_PRECISION)


def _aggregate_strength(breakdown: PhenotypeScoreBreakdown) -> EvidenceStrength:
    """Band the computed score into an ordinal strength.

    Never ``DEFINITIVE``: this is the pipeline's own inference over a curated
    substitute knowledge base, and no arrangement of term overlap makes that a
    definitive statement about a patient.
    """
    if breakdown.annotation_status is GeneAnnotationStatus.NO_ANNOTATIONS:
        return EvidenceStrength.INSUFFICIENT
    if breakdown.specificity is None and breakdown.coverage is None:
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


def _listing(terms: Sequence[str]) -> str:
    return ", ".join(terms) if terms else "none"


def _build_rationale(
    gene_symbol: str, *, partition: _Partition, breakdown: PhenotypeScoreBreakdown
) -> str:
    """One deterministic paragraph explaining the number, including what was ignored."""
    parts = [
        f"{gene_symbol}: matched (observed and gene-associated) {_terms(partition.matched)}; "
        f"contradicted (explicitly EXCLUDED after assessment) {_terms(partition.contradicted)}; "
        f"not assessed (no information, excluded from the denominator per GP-14) "
        f"{_terms(partition.unassessed)}; uncertain (assessed but equivocal, also "
        f"contributing zero) {_terms(partition.uncertain)}; observed but not accounted "
        f"for by this gene {_listing(partition.unexplained)}."
    ]
    if breakdown.mode is ScoringMode.ONTOLOGY:
        parts.append(
            f"Comparison was ontology-aware ({breakdown.mode.value}): matched via upward "
            f"entailment only {_listing(partition.implied_matched)}; contradicted via "
            f"downward entailment only {_listing(partition.implied_contradicted)}."
        )
    if breakdown.annotation_status is GeneAnnotationStatus.NO_ANNOTATIONS:
        parts.append(
            f"The knowledge base holds no curated phenotype terms for {gene_symbol}, so "
            f"both components are undefined and the score is the neutral "
            f"{NEUTRAL_SCORE:.3f} rather than 0: absence of information is not evidence "
            "of a poor match. This is ignorance about the gene, and it is not the same "
            "state as a gene whose curated terms fit nothing (GP-14)."
        )
    elif breakdown.specificity is None and breakdown.coverage is None:
        parts.append(
            f"No gene-associated term was informative in this profile, so the score is the "
            f"neutral {NEUTRAL_SCORE:.3f} rather than 0: absence of information is not "
            "evidence of a poor match."
        )
    else:
        specificity_text = (
            f"{breakdown.specificity:.3f}"
            if breakdown.specificity is not None
            else f"undefined (contributes the neutral {NEUTRAL_SCORE:.3f})"
        )
        coverage_text = (
            f"{breakdown.coverage:.3f}"
            if breakdown.coverage is not None
            else f"undefined (contributes the neutral {NEUTRAL_SCORE:.3f})"
        )
        parts.append(
            f"Specificity {specificity_text} "
            f"(matched weight {breakdown.matched_weight:.4f} of informative weight "
            f"{breakdown.informative_weight:.4f}); coverage {coverage_text}. "
            f"Score {breakdown.score:.3f} = "
            f"{SPECIFICITY_WEIGHT:g}*specificity + {COVERAGE_WEIGHT:g}*coverage."
        )
    if breakdown.symmetric_similarity is not None:
        parts.append(
            f"Phenomizer symmetric best-match average {breakdown.symmetric_similarity:.3f} "
            "(reported for audit; not scored on, because its gene-ward direction "
            "penalises a gene for curated features nobody assessed)."
        )
    if breakdown.uncomparable_observed_count:
        parts.append(
            f"{breakdown.uncomparable_observed_count} observed term(s) could not be placed "
            "in the annotation corpus and were excluded from the coverage mean rather than "
            "scored as zero."
        )
    if breakdown.ic_missing_term_count:
        parts.append(
            f"{breakdown.ic_missing_term_count} gene term(s) have no information content in "
            "the corpus and kept their full curated weight unscaled."
        )
    return " ".join(parts)


# ---------------------------------------------------------------------------
# Evidence (GP-10 / GP-17)
# ---------------------------------------------------------------------------


def _mode_limitation(semantics: PhenotypeSemantics | None) -> str:
    """The mode-specific caveat every evidence item from this module carries."""
    if semantics is None:
        return _EXACT_MODE_LIMITATION
    return (
        f"{_ONTOLOGY_MODE_LIMITATION} Scored against HPO release "
        f"{semantics.data_version} with information content from "
        f"{semantics.corpus.source_name} ({semantics.information_content.entity_count} "
        f"annotated {semantics.corpus.kind.value}s)."
    )


def _mode_method(semantics: PhenotypeSemantics | None) -> str:
    """The mode-specific method prefix, naming the published measure in use."""
    if semantics is None:
        return (
            "Exact HPO identifier equality; no ontology release loaded, so no ancestor "
            "closure and no information content were used."
        )
    return (
        f"Ontology-aware comparison over HPO {semantics.data_version}. "
        f"OBSERVED terms entail their ancestors; EXCLUDED terms entail their descendants. "
        f"Pairwise similarity: {semantics.measure.value} over corpus information content. "
        f"{semantics.method_citation()}"
    )


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
    semantics: PhenotypeSemantics | None,
    conflicts: Sequence[str],
) -> tuple[EvidenceItem, ...]:
    """Assemble the evidence ledger for one gene, in a fixed deterministic order.

    Order is: the computed match, the no-annotation notice where it applies, then
    the per-term observations grouped by status (matched, contradicted, unassessed,
    uncertain), then the unexplained-features summary, then any record conflict.
    Within each group terms are iterated in sorted HPO-id order, so the sequence is
    fully determined by the inputs (GP-30).
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
    ontology_citation = (
        Citation(
            source="human_phenotype_ontology",
            identifier=semantics.corpus.source_name,
            version=semantics.data_version,
            url=None,
            title="Human Phenotype Ontology",
        )
        if semantics is not None
        else profile_citation
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
            semantics=semantics,
        )
    ]
    no_annotation_item = _no_annotation_item(
        gene_symbol, breakdown=breakdown, citation=index_citation, clock=clock, semantics=semantics
    )
    if no_annotation_item is not None:
        items.append(no_annotation_item)
    items.extend(
        _observation_items(
            gene_symbol,
            partition=partition,
            citation=profile_citation,
            clock=clock,
            semantics=semantics,
        )
    )
    items.extend(
        _information_gap_items(
            gene_symbol,
            partition=partition,
            citation=profile_citation,
            clock=clock,
            semantics=semantics,
        )
    )
    unexplained_item = _unexplained_item(
        gene_symbol,
        partition=partition,
        citation=profile_citation,
        clock=clock,
        semantics=semantics,
    )
    if unexplained_item is not None:
        items.append(unexplained_item)
    conflict_item = _conflict_item(
        gene_symbol,
        conflicts=conflicts,
        citation=ontology_citation,
        clock=clock,
        semantics=semantics,
    )
    if conflict_item is not None:
        items.append(conflict_item)
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
    semantics: PhenotypeSemantics | None,
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
    payload: dict[str, str | int | float | bool | None] = {
        "mode": breakdown.mode.value,
        "annotation_status": breakdown.annotation_status.value,
        "specificity": breakdown.specificity,
        "coverage": breakdown.coverage,
        "symmetric_similarity": breakdown.symmetric_similarity,
        "matched_weight": breakdown.matched_weight,
        "contradicted_weight": breakdown.contradicted_weight,
        "informative_weight": breakdown.informative_weight,
        "matched_count": len(partition.matched),
        "contradicted_count": len(partition.contradicted),
        "unassessed_count": len(partition.unassessed),
        "uncertain_count": len(partition.uncertain),
        "unexplained_count": len(partition.unexplained),
        "implied_matched_count": len(partition.implied_matched),
        "implied_contradicted_count": len(partition.implied_contradicted),
        "ic_missing_term_count": breakdown.ic_missing_term_count,
        "uncomparable_observed_count": breakdown.uncomparable_observed_count,
        "observed_term_count": breakdown.observed_term_count,
    }
    if semantics is not None:
        payload["hpo_data_version"] = semantics.data_version
        payload["similarity_measure"] = semantics.measure.value
        payload["ic_corpus"] = semantics.corpus.source_name

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
            "association weight, and coverage measures how much of the patient's observed "
            "presentation the gene accounts for. NOT_ASSESSED and UNCERTAIN terms are "
            "excluded from both numerator and denominator, and an undefined component "
            f"contributes the neutral {NEUTRAL_SCORE}. {_mode_method(semantics)} "
            f"Rationale: {rationale}"
        ),
        limitations=(
            f"{_CURATION_LIMITATION} {_mode_limitation(semantics)} The component weights "
            f"({SPECIFICITY_WEIGHT:g}/{COVERAGE_WEIGHT:g}) are heuristics chosen for the "
            "demo case, not values calibrated against a labelled cohort, so the number "
            "ranks candidates and does not estimate a probability. The gene-phenotype "
            "knowledge base backing this score is a synthetic substitute (GP-20)."
        ),
        citation=citation,
        clock=clock,
        numeric_value=breakdown.score,
        payload=payload,
    )


def _no_annotation_item(
    gene_symbol: str,
    *,
    breakdown: PhenotypeScoreBreakdown,
    citation: Citation,
    clock: Clock,
    semantics: PhenotypeSemantics | None,
) -> EvidenceItem | None:
    """Explicit notice that the knowledge base holds nothing for this gene.

    Emitted so that "no curated terms" appears in the ledger as its own fact. A
    reader who sees only the neutral score cannot otherwise tell it apart from a
    gene whose curated terms happened to be unassessed in this patient, and those
    are different reasons to be uncertain.
    """
    if breakdown.annotation_status is not GeneAnnotationStatus.NO_ANNOTATIONS:
        return None
    return _make_item(
        gene_symbol=gene_symbol,
        claim=(
            f"The gene-phenotype knowledge base holds no curated HPO terms for "
            f"{gene_symbol}; phenotype similarity is undefined and the neutral "
            f"{NEUTRAL_SCORE} was used."
        ),
        direction=EvidenceDirection.NEUTRAL,
        strength=EvidenceStrength.INSUFFICIENT,
        evidence_type=EvidenceType.PIPELINE_INFERENCE,
        tier=AssertionTier.INFERENCE,
        method=(
            "Lookup returned an empty association set for this gene symbol. No "
            f"comparison was attempted. {_mode_method(semantics)}"
        ),
        limitations=f"{_NO_ANNOTATION_LIMITATION} {_mode_limitation(semantics)}",
        citation=citation,
        clock=clock,
        numeric_value=0.0,
        payload={
            "annotation_status": breakdown.annotation_status.value,
            "curated_term_count": 0,
        },
    )


def _observation_items(
    gene_symbol: str,
    *,
    partition: _Partition,
    citation: Citation,
    clock: Clock,
    semantics: PhenotypeSemantics | None,
) -> list[EvidenceItem]:
    """OBSERVED_DATA items for terms the clinician actually assessed."""
    implied_matched = frozenset(partition.implied_matched)
    implied_contradicted = frozenset(partition.implied_contradicted)
    items: list[EvidenceItem] = []
    for assoc in partition.matched:
        via_graph = assoc.hpo_id in implied_matched
        items.append(
            _make_item(
                gene_symbol=gene_symbol,
                claim=(
                    f"{assoc.hpo_id} ({assoc.label}) is "
                    + (
                        "entailed OBSERVED in this subject (a more specific descendant term "
                        "was recorded as present)"
                        if via_graph
                        else "OBSERVED in this subject"
                    )
                    + f" and is {_strength_phrase(assoc)} phenotype association of "
                    f"{gene_symbol}."
                ),
                direction=EvidenceDirection.SUPPORTS,
                strength=_evidence_strength(assoc),
                evidence_type=(
                    EvidenceType.PIPELINE_INFERENCE
                    if via_graph
                    else EvidenceType.DIRECT_MEASUREMENT
                ),
                tier=AssertionTier.INFERENCE if via_graph else AssertionTier.OBSERVED_DATA,
                method=(
                    (
                        "Upward is_a entailment: an observed descendant of this term implies "
                        "it. Ancestors only — never siblings and never other descendants. "
                        if via_graph
                        else "HPO term intersection with the subject's OBSERVED terms. "
                    )
                    + f"Gene association source: {assoc.source} {assoc.version}. "
                    + _mode_method(semantics)
                ),
                limitations=(
                    f"{_CURATION_LIMITATION} {_mode_limitation(semantics)} A shared term "
                    "shows the gene can produce this feature, not that it did so in this "
                    "subject; common features such as microcephaly are associated with "
                    "hundreds of genes."
                    + (
                        " This match is entailed by the ontology rather than recorded "
                        "directly, so it is weaker than an exact term match: the clinician "
                        "asserted the specific finding, not this general one."
                        if via_graph
                        else ""
                    )
                ),
                citation=citation,
                clock=clock,
                numeric_value=assoc.weight_or(UNCURATED_ASSOCIATION_WEIGHT),
                payload={
                    "hpo_id": assoc.hpo_id,
                    "status": ObservationStatus.OBSERVED.value,
                    "association_strength": assoc.association_strength or "",
                    "association_strength_source": assoc.association_strength_source or "",
                    "hpo_frequency": (
                        assoc.hpo_frequency.raw if assoc.hpo_frequency is not None else ""
                    ),
                    "entailed": via_graph,
                },
            )
        )
    for assoc in partition.contradicted:
        via_graph = assoc.hpo_id in implied_contradicted
        items.append(
            _make_item(
                gene_symbol=gene_symbol,
                claim=(
                    f"{assoc.hpo_id} ({assoc.label}) was "
                    + (
                        "entailed absent in this subject (a more general ancestor term was "
                        "assessed and EXCLUDED)"
                        if via_graph
                        else "assessed and EXCLUDED in this subject"
                    )
                    + f", yet is {_strength_phrase(assoc)} phenotype association of "
                    f"{gene_symbol}; this argues against {gene_symbol}."
                ),
                direction=EvidenceDirection.CONTRADICTS,
                strength=_evidence_strength(assoc),
                evidence_type=(
                    EvidenceType.PIPELINE_INFERENCE
                    if via_graph
                    else EvidenceType.DIRECT_MEASUREMENT
                ),
                tier=AssertionTier.INFERENCE if via_graph else AssertionTier.OBSERVED_DATA,
                method=(
                    (
                        "Downward is_a entailment: excluding a term excludes every descendant "
                        "of it. Descendants only — never ancestors, since a specific finding "
                        "being absent says nothing about the general one. "
                        if via_graph
                        else "HPO term intersection with the subject's EXCLUDED terms. "
                    )
                    + "EXCLUDED is the only status permitted to carry negative weight. "
                    + f"Gene association source: {assoc.source} {assoc.version}. "
                    + _mode_method(semantics)
                ),
                limitations=(
                    f"{_CURATION_LIMITATION} {_mode_limitation(semantics)} A curated "
                    "association is not fully penetrant and may be age-dependent: a feature "
                    "genuinely absent today can appear later, so this contradiction "
                    "down-ranks the gene and must not be used to delete it (GP-13, GP-19)."
                    + (
                        " The exclusion is entailed from a broader negative finding rather "
                        "than recorded for this term, so it rests on the assessment having "
                        "been as thorough as the broad term implies."
                        if via_graph
                        else ""
                    )
                ),
                citation=citation,
                clock=clock,
                numeric_value=assoc.weight_or(UNCURATED_ASSOCIATION_WEIGHT),
                payload={
                    "hpo_id": assoc.hpo_id,
                    "status": ObservationStatus.EXCLUDED.value,
                    "association_strength": assoc.association_strength or "",
                    "association_strength_source": assoc.association_strength_source or "",
                    "hpo_frequency": (
                        assoc.hpo_frequency.raw if assoc.hpo_frequency is not None else ""
                    ),
                    "entailed": via_graph,
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
    semantics: PhenotypeSemantics | None,
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
                    f"{assoc.hpo_id} ({assoc.label}) is {_strength_phrase(assoc)} "
                    f"phenotype association of {gene_symbol} but was NOT ASSESSED in this "
                    "subject; it contributes zero to the score in either direction."
                ),
                direction=EvidenceDirection.NEUTRAL,
                strength=EvidenceStrength.INSUFFICIENT,
                evidence_type=EvidenceType.PIPELINE_INFERENCE,
                tier=AssertionTier.INFERENCE,
                method=(
                    "Term present in the gene's curated associations and neither recorded "
                    "nor entailed as OBSERVED or EXCLUDED for this subject. Excluded from "
                    f"both the numerator and the denominator of the score. "
                    f"{_mode_method(semantics)}"
                ),
                limitations=(
                    f"{_NOT_ASSESSED_LIMITATION} It is an actionable gap rather than a "
                    "finding: assessing this feature could move the score in either "
                    f"direction. {_mode_limitation(semantics)}"
                ),
                citation=citation,
                clock=clock,
                numeric_value=0.0,
                payload={
                    "hpo_id": assoc.hpo_id,
                    "status": ObservationStatus.NOT_ASSESSED.value,
                    "association_strength": assoc.association_strength or "",
                    "association_strength_source": assoc.association_strength_source or "",
                    "hpo_frequency": (
                        assoc.hpo_frequency.raw if assoc.hpo_frequency is not None else ""
                    ),
                },
            )
        )
    for assoc in partition.uncertain:
        items.append(
            _make_item(
                gene_symbol=gene_symbol,
                claim=(
                    f"{assoc.hpo_id} ({assoc.label}) is {_strength_phrase(assoc)} "
                    f"phenotype association of {gene_symbol} but was recorded as UNCERTAIN "
                    "in this subject; it contributes zero to the score in either direction."
                ),
                direction=EvidenceDirection.NEUTRAL,
                strength=EvidenceStrength.INSUFFICIENT,
                evidence_type=EvidenceType.PIPELINE_INFERENCE,
                tier=AssertionTier.INFERENCE,
                method=(
                    "Term present in the gene's curated associations and recorded with "
                    "status 'uncertain'. An equivocal finding entails nothing about its "
                    "ancestors or its descendants, so it is not propagated in either "
                    "direction. Excluded from both the numerator and the denominator of "
                    f"the score. {_mode_method(semantics)}"
                ),
                limitations=(
                    "An equivocal assessment is not a negative one. Treating it as either "
                    "presence or absence would fabricate certainty the record does not "
                    "contain; re-assessment, not re-weighting, is the remedy. "
                    f"{_mode_limitation(semantics)}"
                ),
                citation=citation,
                clock=clock,
                numeric_value=0.0,
                payload={
                    "hpo_id": assoc.hpo_id,
                    "status": ObservationStatus.UNCERTAIN.value,
                    "association_strength": assoc.association_strength or "",
                    "association_strength_source": assoc.association_strength_source or "",
                    "hpo_frequency": (
                        assoc.hpo_frequency.raw if assoc.hpo_frequency is not None else ""
                    ),
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
    semantics: PhenotypeSemantics | None,
) -> EvidenceItem | None:
    """Summary item for observed features the gene does not account for."""
    if not partition.unexplained:
        return None
    return _make_item(
        gene_symbol=gene_symbol,
        claim=(
            f"{len(partition.unexplained)} OBSERVED feature(s) are not accounted for by "
            f"{gene_symbol}: {', '.join(partition.unexplained)}."
        ),
        direction=EvidenceDirection.NEUTRAL,
        strength=EvidenceStrength.INSUFFICIENT,
        evidence_type=EvidenceType.PIPELINE_INFERENCE,
        tier=AssertionTier.INFERENCE,
        method=(
            (
                "Observed terms with no is_a subsumption relation, in either direction, to "
                "any curated term of this gene. "
                if semantics is not None
                else "Set difference between the subject's OBSERVED terms and the gene's "
                "curated associations. "
            )
            + f"Lowers the coverage component of the score. {_mode_method(semantics)}"
        ),
        limitations=(
            "Recorded as NEUTRAL, not as a contradiction: an unexplained feature may come "
            "from a second condition, from prematurity or from a treatment effect, and the "
            "curated association list for any gene is incomplete. It lowers coverage, and "
            "must not be reported as evidence that the gene is wrong. "
            f"{_mode_limitation(semantics)}"
        ),
        citation=citation,
        clock=clock,
        numeric_value=float(len(partition.unexplained)),
        payload={"unexplained_count": len(partition.unexplained)},
    )


def _conflict_item(
    gene_symbol: str,
    *,
    conflicts: Sequence[str],
    citation: Citation,
    clock: Clock,
    semantics: PhenotypeSemantics | None,
) -> EvidenceItem | None:
    """Notice that the patient record entails a term both present and absent.

    Emitted per gene because the ledger is keyed by subject, but the finding is
    about the record, not the gene. It is never used to adjust a score: an
    inconsistent record is a data-quality problem for a human to resolve, and a
    scorer that quietly picks a side hides it.
    """
    if not conflicts or semantics is None:
        return None
    listed = ", ".join(conflicts[:10])
    suffix = "" if len(conflicts) <= 10 else f" (and {len(conflicts) - 10} more)"
    return _make_item(
        gene_symbol=gene_symbol,
        claim=(
            f"The phenotype record entails {len(conflicts)} HPO term(s) as both present and "
            f"absent: {listed}{suffix}. The record is internally inconsistent."
        ),
        direction=EvidenceDirection.NEUTRAL,
        strength=EvidenceStrength.INSUFFICIENT,
        evidence_type=EvidenceType.PIPELINE_INFERENCE,
        tier=AssertionTier.INFERENCE,
        method=(
            "Two checks, unioned. (1) Rows grouped by resolved primary HPO id: a term "
            "asserted both OBSERVED and EXCLUDED under two spellings of one identifier "
            "is demoted to UNCERTAIN and contributes nothing in either direction. "
            "(2) Intersection of the upward closure of OBSERVED terms with the downward "
            "closure of EXCLUDED terms: a non-empty intersection means an observed "
            "finding sits under an ancestor that was assessed and excluded. "
            f"{_mode_method(semantics)}"
        ),
        limitations=(
            "Reported, not resolved. Choosing between two contradictory clinical "
            "assertions is not a decision a scorer can make, and doing it silently "
            "would hide a data-entry error a human needs to see. Where the "
            "contradiction is across two different terms the directly recorded status "
            "is kept for each and no score is adjusted; where it is about a single "
            "term under two spellings the term is withdrawn from scoring entirely "
            "rather than resolved in either direction. Both cases leave the affected "
            "term unable to support or argue against any gene until the record is "
            f"corrected. {_mode_limitation(semantics)}"
        ),
        citation=citation,
        clock=clock,
        numeric_value=float(len(conflicts)),
        payload={
            "conflicting_term_count": len(conflicts),
            "aggregation_citation": AGGREGATION_CITATION,
        },
    )
