"""VCF ingestion, normalisation and analytical QC.

The stage's contract, in order:

1. :func:`read_vcf` turns a VCF into sorted, build-anchored :class:`VariantRecord`s,
   decomposing multiallelic sites and hard-filtering only the un-analysable.
2. :func:`normalise_variants` trims to a parsimonious representation, left-aligns
   when a :class:`ReferenceLookup` is supplied, and flags REF disagreements.
3. :func:`assess_quality` adds non-destructive QC flags and one analytical
   :class:`~mva.models.evidence.EvidenceItem` per variant.

Nothing here reaches the network, and nothing here echoes a record's contents into
a warning, a metric or an exception message.
"""

from mva.ingestion.normalise import (
    OP_LEFT_ALIGN,
    OP_SPLIT_MULTIALLELIC,
    OP_TRIM,
    REF_ALLELE_MISMATCH_FLAG,
    NormalisationResult,
    ReferenceLookup,
    normalise_variants,
    split_multiallelic,
    trim_and_left_align,
)
from mva.ingestion.qc import (
    FLAG_FILTERED_BY_CALLER,
    FLAG_HIGH_ALLELE_BALANCE,
    FLAG_LOW_ALLELE_BALANCE,
    FLAG_LOW_DEPTH,
    FLAG_LOW_GQ,
    FLAG_NO_CALLER_FILTER,
    FLAG_NO_QUALITY_METRICS,
    FLAG_ORDER,
    FLAG_POSSIBLE_MOSAIC,
    NEUTRAL_FLAGS,
    QcResult,
    allele_fraction,
    assess_quality,
)
from mva.ingestion.reader import (
    BACKEND_AUTO,
    BACKEND_CYVCF2,
    BACKEND_TEXT,
    SUPPORTED_BACKENDS,
    IngestionResult,
    detect_backend,
    read_vcf,
)

__all__ = [
    "BACKEND_AUTO",
    "BACKEND_CYVCF2",
    "BACKEND_TEXT",
    "FLAG_FILTERED_BY_CALLER",
    "FLAG_HIGH_ALLELE_BALANCE",
    "FLAG_LOW_ALLELE_BALANCE",
    "FLAG_LOW_DEPTH",
    "FLAG_LOW_GQ",
    "FLAG_NO_CALLER_FILTER",
    "FLAG_NO_QUALITY_METRICS",
    "FLAG_ORDER",
    "FLAG_POSSIBLE_MOSAIC",
    "NEUTRAL_FLAGS",
    "OP_LEFT_ALIGN",
    "OP_SPLIT_MULTIALLELIC",
    "OP_TRIM",
    "REF_ALLELE_MISMATCH_FLAG",
    "SUPPORTED_BACKENDS",
    "IngestionResult",
    "NormalisationResult",
    "QcResult",
    "ReferenceLookup",
    "allele_fraction",
    "assess_quality",
    "detect_backend",
    "normalise_variants",
    "read_vcf",
    "split_multiallelic",
    "trim_and_left_align",
]
