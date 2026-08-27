"""Mechanism construction: from a candidate gene to an auditable causal chain.

The stage turns curated knowledge into a `MechanismHypothesis` — an explicit
chain of evidenced links, each with its own signed direction — and reports what
the chain does *not* establish. It imports no peer stage (GP-03); the pipeline
passes its output into the intervention stage as an argument.

Maturity (GP-20): **synthetic-substitute**. The chain tables shipped here are
fabricated for a fictional demo gene and are not biologically valid.
"""

from __future__ import annotations

from mva.mechanisms.builder import (
    DEMONSTRATED_FRACTION_WEIGHT,
    INFERRED_LINK_PENALTY,
    LINK_STRENGTH_WEIGHT,
    MEAN_STRENGTH_WEIGHT,
    SYNTHETIC_CHAIN_LIMITATION,
    MechanismResult,
    build_mechanism,
    mechanism_relevance_score,
)
from mva.mechanisms.library import MechanismLibrary

__all__ = [
    "DEMONSTRATED_FRACTION_WEIGHT",
    "INFERRED_LINK_PENALTY",
    "LINK_STRENGTH_WEIGHT",
    "MEAN_STRENGTH_WEIGHT",
    "SYNTHETIC_CHAIN_LIMITATION",
    "MechanismLibrary",
    "MechanismResult",
    "build_mechanism",
    "mechanism_relevance_score",
]
