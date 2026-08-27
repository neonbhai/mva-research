"""Variant prioritisation: filter, pair, score, rank.

The stage that turns a flagged variant set into a ranked list of candidate
hypotheses, each carrying its component score vector, its evidence (supporting
*and* contradicting), and the experiments that would resolve it.

Composition happens in ``mva.pipeline``: this package imports no peer stage, and
takes the phenotype and mechanism scores as arguments (GP-03).
"""

from __future__ import annotations

from mva.prioritization.filters import (
    FLAG_BENIGN_CONSEQUENCE,
    FLAG_COMMON_VARIANT,
    FLAG_HOMOZYGOUS_CALL,
    FLAG_LOW_FREQUENCY_VARIANT,
    FLAG_LOW_QUALITY_CALL,
    FLAG_NO_FREQUENCY_DATA,
    FLAG_PLAUSIBLE_CANDIDATE,
    FLAG_POSSIBLE_MOSAIC,
    HARD_FILTER_REASONS,
    FilterResult,
    apply_hard_filters,
    apply_soft_flags,
    select_candidate_variants,
)
from mva.prioritization.pairing import (
    FLAG_PHASE_CIS,
    FLAG_PHASE_TRANS,
    FLAG_PHASE_UNKNOWN,
    FLAG_SINGLE_VARIANT,
    PairCandidate,
    generate_pairs,
    infer_phase,
)
from mva.prioritization.ranking import (
    assign_discriminating_experiments,
    rank_pairs,
)
from mva.prioritization.scoring import (
    NEUTRAL_MECHANISM_SCORE,
    NEUTRAL_PHENOTYPE_SCORE,
    Contradiction,
    ScoredPair,
    collect_contradictions,
    composite_score,
    contradiction_penalty,
    score_analytical_validity,
    score_consequence,
    score_evidence_quality,
    score_inheritance,
    score_pair,
    score_rarity,
)

__all__ = [
    "FLAG_BENIGN_CONSEQUENCE",
    "FLAG_COMMON_VARIANT",
    "FLAG_HOMOZYGOUS_CALL",
    "FLAG_LOW_FREQUENCY_VARIANT",
    "FLAG_LOW_QUALITY_CALL",
    "FLAG_NO_FREQUENCY_DATA",
    "FLAG_PHASE_CIS",
    "FLAG_PHASE_TRANS",
    "FLAG_PHASE_UNKNOWN",
    "FLAG_PLAUSIBLE_CANDIDATE",
    "FLAG_POSSIBLE_MOSAIC",
    "FLAG_SINGLE_VARIANT",
    "HARD_FILTER_REASONS",
    "NEUTRAL_MECHANISM_SCORE",
    "NEUTRAL_PHENOTYPE_SCORE",
    "Contradiction",
    "FilterResult",
    "PairCandidate",
    "ScoredPair",
    "apply_hard_filters",
    "apply_soft_flags",
    "assign_discriminating_experiments",
    "collect_contradictions",
    "composite_score",
    "contradiction_penalty",
    "generate_pairs",
    "infer_phase",
    "rank_pairs",
    "score_analytical_validity",
    "score_consequence",
    "score_evidence_quality",
    "score_inheritance",
    "score_pair",
    "score_rarity",
    "select_candidate_variants",
]
