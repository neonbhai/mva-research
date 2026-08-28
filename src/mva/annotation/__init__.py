"""Annotation stage: attach consequences, frequencies and clinical assertions.

The stage is adapter-shaped on purpose. `base` defines the Protocols an
annotation source must satisfy, `local_tables` implements them over hash-pinned
synthetic TSVs under `knowledge/public/`, and `service` orchestrates them into
annotated records plus an evidence trail.

Two adapter sets live here and the case's `synthetic` flag chooses between them
(ADR 0027):

* `local_tables` — hash-pinned synthetic TSVs under `knowledge/public/`. Fictional
  genes, invented allele frequencies, no ClinVar. What every demo and every test
  runs on.
* `binding` — the real releases: ClinVar, gnomAD v4.1 exomes, SnpEff and MANE,
  resolved from `config.resources` under `$MVA_RESOURCES`. What a non-synthetic
  case runs on, or the run refuses to start.

Which one is in force is stated by the adapters themselves (`name`/`version`), by
every EvidenceItem's `limitations` and by the run warnings, so a synthetic
substitute cannot quietly stop announcing that it is one (GP-20).

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
from mva.annotation.binding import (
    AdapterBindingError,
    BoundAdapters,
    PartialCoverageFrequencyAdapter,
    ResolvedResources,
    build_clinical_adapter,
    build_consequence_adapter,
    build_frequency_adapter,
    build_real_adapter_set,
    representation_warnings,
    resolve_real_resources,
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
from mva.annotation.service import (
    DEFAULT_ANNOTATION_BATCH_SIZE,
    AnnotatedVariant,
    AnnotationResult,
    AnnotationStream,
    annotate_variants,
    iter_annotated,
)

__all__ = [
    "CONSEQUENCE_TABLE",
    "DEFAULT_ANNOTATION_BATCH_SIZE",
    "FREQUENCY_TABLE",
    "SYNTHETIC_STANDIN_LIMITATION",
    "AdapterBindingError",
    "AdapterDescriptor",
    "AdapterRole",
    "AdapterSet",
    "AnnotatedVariant",
    "AnnotationResult",
    "AnnotationStream",
    "BoundAdapters",
    "ClinicalAdapter",
    "ConsequenceAdapter",
    "FrequencyAdapter",
    "KnowledgeManifest",
    "KnowledgeTable",
    "LocalConsequenceAdapter",
    "LocalFrequencyAdapter",
    "MaturityAware",
    "NullClinicalAdapter",
    "PartialCoverageFrequencyAdapter",
    "ResolvedResources",
    "annotate_variants",
    "build_clinical_adapter",
    "build_consequence_adapter",
    "build_frequency_adapter",
    "build_real_adapter_set",
    "compute_manifest",
    "is_synthetic",
    "iter_annotated",
    "load_default_adapters",
    "load_manifest",
    "render_manifest_yaml",
    "representation_warnings",
    "resolve_real_resources",
    "resolve_table_path",
    "verify_manifest",
]
