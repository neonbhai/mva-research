"""Deterministic ranking and next-experiment assignment.

Two jobs. First, impose a **total** order on scored candidates so that a repeat
run is byte-identical (GP-30): composite descending, then genomic position
ascending, then the deterministic pair id as a final tiebreak. Sorting on the
composite alone is not a total order — ties are common once components saturate,
and dict/set iteration order would silently decide them.

Second, say what to do next. A ranked list with no experiment attached is a
league table; a ranked list where every row names the assay that would confirm
or kill it is a research plan. ``recommended_next_test`` answers "what would
resolve THIS candidate", ``discriminating_experiment`` answers "what would tell
rank N apart from rank N+1" — a different and usually cheaper question.
"""

from __future__ import annotations

from collections.abc import Sequence

from mva.clock import Clock
from mva.models.pair import (
    CIS_STATES,
    NO_SECOND_VARIANT,
    CandidatePair,
    InheritanceModel,
    PhaseStatus,
)
from mva.models.variant import Zygosity
from mva.prioritization.filters import (
    FLAG_COMMON_VARIANT,
    FLAG_LOW_QUALITY_CALL,
    FLAG_NO_FREQUENCY_DATA,
)
from mva.prioritization.pairing import PairCandidate
from mva.prioritization.scoring import ScoredPair

FLAG_HAS_CONTRADICTIONS = "has_contradictions"
FLAG_BLOCKING_QUESTIONS = "blocking_questions_open"

#: What a difference in each component would actually take to resolve.
DISCRIMINATING_EXPERIMENTS: dict[str, str] = {
    "analytical_validity": (
        "orthogonal confirmation of both candidates' genotypes by Sanger or targeted "
        "resequencing at adequate depth"
    ),
    "rarity": (
        "genotyping both candidate alleles in an ancestry-matched population reference "
        "with adequate coverage at these loci"
    ),
    "molecular_consequence": (
        "a functional readout of the predicted lesion — RT-PCR or RNA-seq of patient cells "
        "for splice predictions, protein assay for missense"
    ),
    "inheritance_consistency": (
        "parental segregation testing, or long-read sequencing to phase the alleles directly"
    ),
    "phenotype_similarity": (
        "deeper phenotyping against both genes' reported HPO profiles, recording explicit "
        "exclusions rather than leaving features unassessed"
    ),
    "mechanistic_relevance": (
        "a pathway-level functional readout in patient-derived cells for both genes"
    ),
    "evidence_quality": (
        "curated review of both loci against a pinned ClinVar/literature release to add "
        "evidence independent of this pipeline"
    ),
    "contradiction_penalty": (
        "resolving the recorded contradictions on the lower-ranked candidate before "
        "spending on either"
    ),
}

_DEFAULT_NEXT_TEST = (
    "Deep phenotyping and curated review of this gene, then a functional assay of the "
    "predicted lesion in patient-derived cells."
)


def _sort_key(
    scored: ScoredPair,
) -> tuple[float, tuple[int, int, str, str], tuple[int, int, str, str], str]:
    """The same total order as :meth:`mva.models.pair.CandidatePair.sort_key`.

    It is spelled separately because ranking runs *before* the ``CandidatePair``
    exists: the score lives on the :class:`~mva.prioritization.scoring.ScoredPair`
    wrapper and the candidate is still a
    :class:`~mva.prioritization.pairing.PairCandidate`. The component order and
    the :data:`~mva.models.pair.NO_SECOND_VARIANT` sentinel are shared with the
    model so the two keys cannot disagree about which candidate ranks first.
    """
    candidate = scored.candidate
    second = NO_SECOND_VARIANT if candidate.variant_b is None else candidate.variant_b.sort_key()
    return (-scored.composite, candidate.variant_a.sort_key(), second, candidate.pair_id)


def _recommended_next_test(candidate: PairCandidate) -> str:
    """The single highest-value next experiment for this candidate.

    Ordered by what would invalidate the candidate soonest. Call validity comes
    first: if the genotype is not real, nothing downstream of it is worth doing.
    """
    if any(FLAG_LOW_QUALITY_CALL in variant.qc_flags for variant in candidate.variants):
        return (
            "Orthogonal confirmation of the flagged genotype by Sanger or targeted "
            "resequencing before any further interpretation — the call itself is in doubt."
        )
    if (
        candidate.inheritance_model is InheritanceModel.COMPOUND_HETEROZYGOUS
        and candidate.phase.status in CIS_STATES
    ):
        return (
            "Long-read sequencing or parental genotypes across both sites, to confirm or "
            "overturn the in-cis phase call that currently disqualifies this hypothesis."
        )
    if candidate.variant_b is not None and candidate.phase.status is PhaseStatus.UNKNOWN:
        return (
            "Parental segregation testing (trio genotyping), or long-read/read-backed "
            "phasing, to establish whether the two variants lie in trans."
        )
    if any(FLAG_NO_FREQUENCY_DATA in variant.qc_flags for variant in candidate.variants):
        return (
            "Look the unrecorded site up in a versioned, ancestry-matched population "
            "reference with adequate coverage; it is currently scored neutrally, not rare."
        )
    if any(FLAG_COMMON_VARIANT in variant.qc_flags for variant in candidate.variants):
        return (
            "Check homozygote counts and phenotype in population reference cohorts; a "
            "healthy homozygote would retire this hypothesis outright."
        )
    if candidate.variant_b is None and candidate.inheritance_model is InheritanceModel.UNKNOWN:
        return (
            "Trio sequencing for de-novo status, plus copy-number analysis of the second "
            "allele in case a deletion is masking a compound heterozygote."
        )
    if any(v.genotype.zygosity is Zygosity.HOM_ALT for v in candidate.variants):
        return (
            "Copy-number and runs-of-homozygosity analysis, to confirm true homozygosity "
            "rather than hemizygosity created by a deletion on the other allele."
        )
    return _DEFAULT_NEXT_TEST


def _to_candidate_pair(scored: ScoredPair, *, rank: int, total: int, stamp: str) -> CandidatePair:
    candidate = scored.candidate
    flags = list(candidate.flags)
    if scored.contradicting_evidence:
        flags.append(FLAG_HAS_CONTRADICTIONS)
    if any(question.blocking for question in scored.open_questions):
        flags.append(FLAG_BLOCKING_QUESTIONS)

    return CandidatePair(
        pair_id=candidate.pair_id,
        gene_symbol=candidate.gene_symbol,
        variant_a=candidate.variant_a,
        variant_b=candidate.variant_b,
        inheritance_model=candidate.inheritance_model,
        phase=candidate.phase,
        scores=scored.scores,
        composite_score=scored.composite,
        rank=rank,
        supporting_evidence_ids=tuple(item.evidence_id for item in scored.supporting_evidence),
        contradicting_evidence_ids=tuple(
            item.evidence_id for item in scored.contradicting_evidence
        ),
        missing_evidence=scored.open_questions,
        recommended_next_test=_recommended_next_test(candidate),
        discriminating_experiment=None,
        rank_rationale=f"{scored.rationale} Ranked {rank} of {total} at {stamp}.",
        flags=tuple(flags),
    )


def rank_pairs(
    scored: Sequence[ScoredPair], *, clock: Clock, max_rank: int | None = None
) -> tuple[CandidatePair, ...]:
    """Order candidates and materialise them as ranked ``CandidatePair`` records.

    Ordering is composite descending, then genomic position ascending, then pair
    id — a total order, so two runs over the same input produce the same list in
    the same order with the same scores. ``max_rank`` truncates the returned
    list only; nothing is deleted from the evidence the caller already holds
    (GP-13, GP-19).
    """
    ordered = sorted(scored, key=_sort_key)
    stamp = clock.now().isoformat()
    total = len(ordered)
    ranked = tuple(
        _to_candidate_pair(item, rank=index, total=total, stamp=stamp)
        for index, item in enumerate(ordered, start=1)
    )
    if max_rank is not None:
        ranked = ranked[:max_rank]
    return assign_discriminating_experiments(ranked)


def _largest_component_gap(higher: CandidatePair, lower: CandidatePair) -> tuple[str, float, float]:
    """The component that most separates two adjacent ranks."""
    upper = higher.scores.as_dict()
    under = lower.scores.as_dict()
    component = max(sorted(upper), key=lambda name: (abs(upper[name] - under[name]), name))
    return component, upper[component], under[component]


def _weakest_component(pair: CandidatePair) -> str:
    scores = pair.scores.as_dict()
    positive = {name: value for name, value in scores.items() if name != "contradiction_penalty"}
    return min(sorted(positive), key=lambda name: (positive[name], name))


def assign_discriminating_experiments(
    pairs: Sequence[CandidatePair],
) -> tuple[CandidatePair, ...]:
    """Attach, to each rank, the experiment that separates it from the next one.

    For the last row there is no next candidate, so it instead names the assay
    that would raise its own weakest component — the cheapest way for it to move
    up the list rather than a comparison against nothing.
    """
    result: list[CandidatePair] = []
    for index, pair in enumerate(pairs):
        if index + 1 < len(pairs):
            following = pairs[index + 1]
            component, upper, under = _largest_component_gap(pair, following)
            experiment = (
                f"Ranks {pair.rank} and {following.rank} differ most in {component} "
                f"({upper:.2f} vs {under:.2f}): "
                f"{DISCRIMINATING_EXPERIMENTS[component]} would confirm or collapse that "
                "difference."
            )
        else:
            component = _weakest_component(pair)
            experiment = (
                f"Lowest-ranked candidate carried forward; its weakest component is "
                f"{component} ({pair.scores.as_dict()[component]:.2f}), so "
                f"{DISCRIMINATING_EXPERIMENTS[component]} is what would move it up the list."
            )
        result.append(pair.model_copy(update={"discriminating_experiment": experiment}))
    return tuple(result)
