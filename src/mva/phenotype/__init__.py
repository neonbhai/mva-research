"""Phenotype representation and gene-phenotype scoring.

Layer 4 (data in). Imports models, errors and clock only; it never imports a peer
stage (GP-03), and the composition root in ``mva.pipeline`` is what joins its
output to variant annotation.

The organising rule of the whole package is GP-14: **absence of information is
not negative information**. ``ObservationStatus`` is four-valued and stays that
way from the TSV parser through to the score's denominator. The two places that
would normally destroy the distinction — a status parser with a default branch,
and a similarity score that divides by the gene's total term count — are both
handled explicitly and both covered by regression tests.

Maturity (GP-20): ``synthetic-substitute``. The parsing, four-valued logic,
evidence emission and determinism are real; the gene-phenotype knowledge base
behind them is fabricated, and there is no ontology ancestor closure yet.

Typical use::

    profile = load_phenotype_profile(
        path, subject_id="PROBAND01", hpo_version="2025-05-06",
        source_artifact="synthetic_phenotype.tsv",
    )
    index = GenePhenotypeIndex.from_tsv(knowledge_path, version="v0.0-synthetic")
    matches = score_all_genes(["SYNTHKIN1"], profile=profile, index=index, clock=clock)
"""

from __future__ import annotations

from mva.phenotype.hpo import (
    GENE_PHENOTYPE_COLUMNS,
    HPO_ID_PATTERN,
    STRENGTH_WEIGHTS,
    GeneAssociation,
    GenePhenotypeIndex,
    normalise_hpo_id,
    read_tsv_rows,
)
from mva.phenotype.loader import (
    DEFAULT_EXTRACTION_CONFIDENCE,
    OPTIONAL_COLUMNS,
    REQUIRED_COLUMNS,
    STATUS_ALIASES,
    allowed_status_summary,
    load_phenotype_profile,
    parse_onset,
    parse_status,
)
from mva.phenotype.scoring import (
    COVERAGE_WEIGHT,
    NEUTRAL_SCORE,
    SPECIFICITY_WEIGHT,
    TOOL_NAME,
    TOOL_VERSION,
    PhenotypeMatch,
    PhenotypeScoreBreakdown,
    score_all_genes,
    score_gene_phenotype,
)

__all__ = [
    "COVERAGE_WEIGHT",
    "DEFAULT_EXTRACTION_CONFIDENCE",
    "GENE_PHENOTYPE_COLUMNS",
    "HPO_ID_PATTERN",
    "NEUTRAL_SCORE",
    "OPTIONAL_COLUMNS",
    "REQUIRED_COLUMNS",
    "SPECIFICITY_WEIGHT",
    "STATUS_ALIASES",
    "STRENGTH_WEIGHTS",
    "TOOL_NAME",
    "TOOL_VERSION",
    "GeneAssociation",
    "GenePhenotypeIndex",
    "PhenotypeMatch",
    "PhenotypeScoreBreakdown",
    "allowed_status_summary",
    "load_phenotype_profile",
    "normalise_hpo_id",
    "parse_onset",
    "parse_status",
    "read_tsv_rows",
    "score_all_genes",
    "score_gene_phenotype",
]
