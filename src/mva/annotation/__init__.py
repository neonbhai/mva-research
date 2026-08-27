"""Annotation stage: attach consequences, frequencies and clinical assertions.

The stage is adapter-shaped on purpose. `base` defines the Protocols an
annotation source must satisfy, `local_tables` implements them over hash-pinned
synthetic TSVs under `knowledge/public/`, and `service` orchestrates them into
annotated records plus an evidence trail.

What ships here is a **synthetic substitute**, not a real annotation pipeline:
fictional genes, invented allele frequencies, no ClinVar. That grade is stated by
the adapters themselves (`name`/`version`), by every EvidenceItem's `limitations`
and by the run warnings, so it cannot quietly stop being true (GP-20).

See `knowledge/adapters/README.md` for how a real VEP / gnomAD / ClinVar /
SpliceAI adapter is slotted in, and why it may never receive proband coordinates
over a network (PRIV-05).
"""

from __future__ import annotations

from mva.annotation.base import (
    SYNTHETIC_STANDIN_LIMITATION,
    AdapterDescriptor,
    AdapterRole,
    AdapterSet,
    ClinicalAdapter,
    ConsequenceAdapter,
    FrequencyAdapter,
    MaturityAware,
    is_synthetic,
)
from mva.annotation.local_tables import (
    CONSEQUENCE_TABLE,
    FREQUENCY_TABLE,
    KnowledgeManifest,
    KnowledgeTable,
    LocalConsequenceAdapter,
    LocalFrequencyAdapter,
    NullClinicalAdapter,
    compute_manifest,
    load_default_adapters,
    load_manifest,
    render_manifest_yaml,
    resolve_table_path,
    verify_manifest,
)
from mva.annotation.service import AnnotationResult, annotate_variants

__all__ = [
    "CONSEQUENCE_TABLE",
    "FREQUENCY_TABLE",
    "SYNTHETIC_STANDIN_LIMITATION",
    "AdapterDescriptor",
    "AdapterRole",
    "AdapterSet",
    "AnnotationResult",
    "ClinicalAdapter",
    "ConsequenceAdapter",
    "FrequencyAdapter",
    "KnowledgeManifest",
    "KnowledgeTable",
    "LocalConsequenceAdapter",
    "LocalFrequencyAdapter",
    "MaturityAware",
    "NullClinicalAdapter",
    "annotate_variants",
    "compute_manifest",
    "is_synthetic",
    "load_default_adapters",
    "load_manifest",
    "render_manifest_yaml",
    "resolve_table_path",
    "verify_manifest",
]
