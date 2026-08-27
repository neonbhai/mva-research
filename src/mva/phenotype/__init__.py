"""Phenotype representation and ontology-aware gene-phenotype scoring.

Layer 4 (data in). Imports models, errors, determinism and clock only; it never
imports a peer stage (GP-03), and the composition root in ``mva.pipeline`` is what
joins its output to variant annotation.

The organising rule of the whole package is GP-14: **absence of information is
not negative information**. ``ObservationStatus`` is four-valued and stays that
way from the TSV parser through the ontology closure to the score's denominator.
The places that would normally destroy the distinction — a status parser with a
default branch, a similarity score that divides by the gene's total term count,
and an ontology closure that propagates a negative finding the wrong way — are
each handled explicitly and each covered by regression tests.

Module map
----------

``hpo``
    Identifier normalisation, TSV reading, and the curated gene→term index.
``loader``
    Phenotype TSV → :class:`~mva.models.phenotype.PhenotypeProfile`.
``ontology``
    ``hp.obo`` → :class:`HpoOntology`, a DAG with multi-parent-safe, cycle-safe
    ancestor and descendant closures.
``corpus``
    ``phenotype.hpoa`` / ``genes_to_phenotype.txt`` → an annotation corpus and the
    corpus-derived :class:`InformationContent`. IC comes from annotation
    frequency, never from graph depth.
``similarity``
    Resnik / Lin / Jiang-Conrath pairwise measures and the Phenomizer best-match
    average, with citations.
``propagation``
    Four-valued status entailment over the DAG. OBSERVED propagates **up**,
    EXCLUDED propagates **down**, UNCERTAIN and NOT_ASSESSED propagate neither way.
``semantics``
    :class:`HpoResourceSet` (paths and pinned digests, supplied by the composition
    root) and the assembled :class:`PhenotypeSemantics` the scorer consumes.
``scoring``
    The score itself, in either :class:`ScoringMode`.

Maturity (GP-20). The ontology, the annotation corpus and the similarity measures
are **real** published resources and the release string is carried into
provenance. The gene→term knowledge base they are pointed at
(``knowledge/public/gene_phenotype.tsv``) is still a ``synthetic-substitute`` with
fictional gene symbols, so no output is biologically valid yet. Those are two
separate claims and must not be merged into one.

Typical use::

    profile = load_phenotype_profile(
        path, subject_id="PROBAND01", hpo_version="2026-06-23",
        source_artifact="synthetic_phenotype.tsv",
    )
    index = GenePhenotypeIndex.from_tsv(knowledge_path, version="v0.0-synthetic")

    # Graph-aware scoring; the composition root owns these paths.
    semantics = HpoResourceSet(
        ontology_path=resources / "hp.obo",
        annotation_path=resources / "phenotype.hpoa",
    ).load()
    matches = score_all_genes(
        ["SYNTHKIN1"], profile=profile, index=index, clock=clock, semantics=semantics
    )

Omitting ``semantics`` falls back to :attr:`ScoringMode.EXACT`, which compares
identifiers only and says so in every evidence item it emits.
"""

from __future__ import annotations

from mva.phenotype.corpus import (
    GENES_TO_PHENOTYPE_COLUMNS,
    HPOA_COLUMNS,
    NEGATION_QUALIFIER,
    PHENOTYPIC_ABNORMALITY_ASPECT,
    AnnotationCorpus,
    CorpusKind,
    CorpusStats,
    InformationContent,
)
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
from mva.phenotype.ontology import (
    DATA_VERSION_KEY,
    HpoOntology,
    HpoTerm,
    ontology_provenance,
)
from mva.phenotype.propagation import (
    InferenceBasis,
    InferredObservation,
    PropagatedProfile,
    PropagationSummary,
)
from mva.phenotype.scoring import (
    COVERAGE_WEIGHT,
    NEUTRAL_SCORE,
    SPECIFICITY_WEIGHT,
    TOOL_NAME,
    TOOL_VERSION,
    GeneAnnotationStatus,
    PhenotypeMatch,
    PhenotypeScoreBreakdown,
    ScoringMode,
    score_all_genes,
    score_gene_phenotype,
)
from mva.phenotype.semantics import HpoResourceSet, PhenotypeSemantics
from mva.phenotype.similarity import (
    AGGREGATION_CITATION,
    MEASURE_CITATIONS,
    BestMatch,
    BestMatchAverage,
    SimilarityMeasure,
    TermSimilarity,
)

__all__ = [
    "AGGREGATION_CITATION",
    "COVERAGE_WEIGHT",
    "DATA_VERSION_KEY",
    "DEFAULT_EXTRACTION_CONFIDENCE",
    "GENES_TO_PHENOTYPE_COLUMNS",
    "GENE_PHENOTYPE_COLUMNS",
    "HPOA_COLUMNS",
    "HPO_ID_PATTERN",
    "MEASURE_CITATIONS",
    "NEGATION_QUALIFIER",
    "NEUTRAL_SCORE",
    "OPTIONAL_COLUMNS",
    "PHENOTYPIC_ABNORMALITY_ASPECT",
    "REQUIRED_COLUMNS",
    "SPECIFICITY_WEIGHT",
    "STATUS_ALIASES",
    "STRENGTH_WEIGHTS",
    "TOOL_NAME",
    "TOOL_VERSION",
    "AnnotationCorpus",
    "BestMatch",
    "BestMatchAverage",
    "CorpusKind",
    "CorpusStats",
    "GeneAnnotationStatus",
    "GeneAssociation",
    "GenePhenotypeIndex",
    "HpoOntology",
    "HpoResourceSet",
    "HpoTerm",
    "InferenceBasis",
    "InferredObservation",
    "InformationContent",
    "PhenotypeMatch",
    "PhenotypeScoreBreakdown",
    "PhenotypeSemantics",
    "PropagatedProfile",
    "PropagationSummary",
    "ScoringMode",
    "SimilarityMeasure",
    "TermSimilarity",
    "allowed_status_summary",
    "load_phenotype_profile",
    "normalise_hpo_id",
    "ontology_provenance",
    "parse_onset",
    "parse_status",
    "read_tsv_rows",
    "score_all_genes",
    "score_gene_phenotype",
]
