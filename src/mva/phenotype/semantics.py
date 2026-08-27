"""The object the scorer is handed: one ontology release plus its similarity model.

Everything ontology-aware in this package is reachable from a single explicitly
constructed :class:`PhenotypeSemantics`. That is a deliberate shape:

* **The composition root owns resource paths.** :class:`HpoResourceSet` takes
  paths and optional pinned digests as plain constructor arguments. No absolute
  path and no digest is hardcoded in this package, so the same code serves the
  full 20k-term release, a committed test subgraph, and a future release without
  an edit.
* **There is no global.** No module-level cache, no ``lru_cache`` on a loader.
  A hidden ontology singleton would make the release invisible at every call
  site, would survive between cases in a long-lived process, and is exactly the
  layering violation this repository's architecture tests exist to catch.
* **Loading is expensive and happens once.** Parsing the release, the annotation
  corpus and the information-content table costs a few seconds and tens of
  megabytes. :meth:`HpoResourceSet.load` does it once and hands back an object the
  caller keeps for the run.

**No network (PRIV-05).** Local files only, supplied by path. Acquiring the
release is a separate offline step outside this package.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from mva.models.phenotype import PhenotypeProfile
from mva.phenotype.corpus import AnnotationCorpus, CorpusKind, InformationContent
from mva.phenotype.ontology import HpoOntology, ontology_provenance
from mva.phenotype.propagation import PropagatedProfile
from mva.phenotype.similarity import (
    AGGREGATION_CITATION,
    MEASURE_CITATIONS,
    SimilarityMeasure,
    TermSimilarity,
)


@dataclass(frozen=True, slots=True)
class HpoResourceSet:
    """Where the HPO release lives and which bytes are expected.

    ``*_sha256`` are optional pins. When supplied, a mismatch is a hard failure
    rather than a warning: a different ontology or annotation release changes every
    information content value and therefore every phenotype score, so an unnoticed
    swap would silently invalidate a golden expectation and a submission alike.
    When omitted, the digest is still computed and recorded, so provenance always
    names the exact bytes.
    """

    ontology_path: Path
    """Path to ``hp.obo``."""

    annotation_path: Path
    """Path to ``phenotype.hpoa`` (disease corpus) or ``genes_to_phenotype.txt``."""

    annotation_kind: CorpusKind = CorpusKind.DISEASE
    """Which file ``annotation_path`` is. Explicit rather than sniffed from the
    filename: the two produce different information content and therefore
    different scores, and provenance must not depend on a naming convention."""

    measure: SimilarityMeasure = SimilarityMeasure.LIN
    """Pairwise measure. Lin (1998) by default — natively bounded to ``[0, 1]`` and
    normalised by the specificity of both terms compared."""

    ontology_sha256: str | None = None
    annotation_sha256: str | None = None

    def load(self) -> PhenotypeSemantics:
        """Parse everything once and return the assembled semantics object."""
        ontology = HpoOntology.from_obo(self.ontology_path, expected_sha256=self.ontology_sha256)
        if self.annotation_kind is CorpusKind.DISEASE:
            corpus = AnnotationCorpus.from_hpoa(
                self.annotation_path,
                ontology=ontology,
                expected_sha256=self.annotation_sha256,
            )
        else:
            corpus = AnnotationCorpus.from_genes_to_phenotype(
                self.annotation_path,
                ontology=ontology,
                expected_sha256=self.annotation_sha256,
            )
        information_content = InformationContent.from_corpus(corpus, ontology=ontology)
        similarity = TermSimilarity(
            ontology=ontology,
            information_content=information_content,
            measure=self.measure,
        )
        return PhenotypeSemantics(
            ontology=ontology,
            corpus=corpus,
            information_content=information_content,
            similarity=similarity,
        )


@dataclass(frozen=True, slots=True)
class PhenotypeSemantics:
    """One loaded HPO release, its annotation corpus, and the similarity model.

    Held by the caller for the duration of a run and passed into
    :func:`mva.phenotype.scoring.score_gene_phenotype`. Cheap to pass around,
    expensive to build; build it once.
    """

    ontology: HpoOntology
    corpus: AnnotationCorpus
    information_content: InformationContent
    similarity: TermSimilarity

    @property
    def measure(self) -> SimilarityMeasure:
        return self.similarity.measure

    @property
    def data_version(self) -> str:
        """The HPO release string, verbatim from the OBO header."""
        return self.ontology.data_version

    def propagate(self, profile: PhenotypeProfile) -> PropagatedProfile:
        """Materialise one subject's entailments over this release.

        Convenience only: the scorer calls this once per profile so that the
        closures are built once rather than once per candidate gene.
        """
        return PropagatedProfile(profile, ontology=self.ontology)

    def method_citation(self) -> str:
        """One sentence naming both published choices behind a score."""
        return f"Pairwise: {MEASURE_CITATIONS[self.measure]}. Aggregation: {AGGREGATION_CITATION}."

    def provenance(self) -> Mapping[str, str]:
        """Everything a reader needs to reproduce these numbers, sorted and stringly typed."""
        fields: dict[str, str] = {
            **ontology_provenance(self.ontology),
            **self.information_content.provenance(),
            "similarity_measure": self.measure.value,
            "similarity_citation": MEASURE_CITATIONS[self.measure],
            "similarity_aggregation_citation": AGGREGATION_CITATION,
            "annotation_source": self.corpus.source_name,
            "annotation_rows_kept": str(self.corpus.stats.rows_kept),
            "annotation_rows_read": str(self.corpus.stats.rows_read),
        }
        if self.corpus.source_sha256 is not None:
            fields["annotation_sha256"] = self.corpus.source_sha256
        return dict(sorted(fields.items()))
