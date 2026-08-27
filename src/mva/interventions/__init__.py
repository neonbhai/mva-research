"""Intervention stage: mechanism-grounded, direction-checked drug hypotheses.

The stage receives a `MechanismHypothesis` as an argument and never imports the
mechanism stage (GP-03). Its job is to answer, for every catalogue entry, the
eight mandatory questions — which node, which direction, does the agent act that
way, on what evidence tier, at a reachable concentration, with what paediatric
exposure, at what oncogenic risk, and which experiment comes first — and to keep
the rejections (GP-19).

Maturity (GP-20): **synthetic-substitute**. The catalogue shipped here describes
fictional agents; nothing produced by this package is pharmacological fact or
medical advice.
"""

from __future__ import annotations

from mva.interventions.catalog import CATALOG_COLUMNS, CatalogEntry, DrugCatalog
from mva.interventions.direction import (
    DirectionVerdict,
    check_direction,
    check_target_in_mechanism,
    corrective_direction,
    required_direction_for_node,
)
from mva.interventions.generate import (
    APPROVAL_STATUS_WEIGHT,
    APPROVAL_WEIGHT,
    DIRECTION_UNKNOWN_CREDIT,
    DIRECTION_WEIGHT,
    EVIDENCE_QUALITY_FLOOR,
    EVIDENCE_WEIGHT,
    INTERVENTION_CLASS_WEIGHT,
    REASON_PRIORITY,
    SAFETY_WEIGHT,
    DrugTriageResult,
    generate_drug_hypotheses,
)
from mva.interventions.safety import (
    EVIDENCE_TYPE_WEIGHT,
    INDIRECT_EVIDENCE_MULTIPLIER,
    SEVERITY_PENALTY,
    SafetyVerdict,
    assess_evidence_quality,
    assess_safety,
    is_chromosomal_instability_context,
    is_neurological_context,
)

__all__ = [
    "APPROVAL_STATUS_WEIGHT",
    "APPROVAL_WEIGHT",
    "CATALOG_COLUMNS",
    "DIRECTION_UNKNOWN_CREDIT",
    "DIRECTION_WEIGHT",
    "EVIDENCE_QUALITY_FLOOR",
    "EVIDENCE_TYPE_WEIGHT",
    "EVIDENCE_WEIGHT",
    "INDIRECT_EVIDENCE_MULTIPLIER",
    "INTERVENTION_CLASS_WEIGHT",
    "REASON_PRIORITY",
    "SAFETY_WEIGHT",
    "SEVERITY_PENALTY",
    "CatalogEntry",
    "DirectionVerdict",
    "DrugCatalog",
    "DrugTriageResult",
    "SafetyVerdict",
    "assess_evidence_quality",
    "assess_safety",
    "check_direction",
    "check_target_in_mechanism",
    "corrective_direction",
    "generate_drug_hypotheses",
    "is_chromosomal_instability_context",
    "is_neurological_context",
    "required_direction_for_node",
]
