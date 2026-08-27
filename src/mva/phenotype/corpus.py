"""Annotation corpora and corpus-derived information content.

Semantic similarity over an ontology needs a notion of how *specific* a term is.
There are two ways to get one, and only one of them is defensible:

* **Depth in the graph.** Cheap, and wrong. HPO's branches are annotated to wildly
  different granularity — the nervous-system subtree is many levels deeper than
  the immune subtree — so depth measures how much curation effort a branch has
  received, not how surprising a finding is. A term four levels down in a deep
  branch and a term four levels down in a shallow one are not comparable, and
  treating them as such is a documented methodological error in the ontology
  literature.
* **Corpus frequency.** ``IC(t) = -ln p(t)``, where ``p(t)`` is the fraction of
  annotated entities in a real annotation release that carry ``t`` *or any of its
  descendants*. This is Resnik's definition, and it is what this module computes.

The descendant part is the **true-path rule**: an entity annotated with
"Microcephaly" is, by the semantics of ``is_a``, also an instance of "Abnormality
of skull size". Counting only the literal annotation would make every internal
term look rarer than every one of its children, which inverts the specificity
ordering the measure exists to provide.

References
----------
Resnik P. (1995) *Using information content to evaluate semantic similarity in a
taxonomy.* IJCAI-95, 448-453.

Köhler S., Schulz M.H., Krawitz P., et al. (2009) *Clinical diagnostics in human
genetics with semantic similarity searches in ontologies.* Am J Hum Genet
85(4):457-464. — the source of the disease-annotation corpus convention used by
:meth:`AnnotationCorpus.from_hpoa`.

**Maturity (GP-20): real.** Both corpora parsed here are the published HPO
annotation release. The gene symbols in ``knowledge/public/gene_phenotype.tsv``
that the scorer compares against remain a synthetic substitute; that is a separate
row in the maturity ledger.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Final

from mva.determinism import hash_file
from mva.errors import IngestionError
from mva.phenotype.ontology import HpoOntology

#: ``phenotype.hpoa`` column order, as declared by the release header. Read by
#: name from the header row rather than by position, so a future column insertion
#: shifts nothing.
HPOA_COLUMNS: Final[tuple[str, ...]] = (
    "database_id",
    "qualifier",
    "hpo_id",
    "aspect",
)

#: ``genes_to_phenotype.txt`` columns this module requires.
GENES_TO_PHENOTYPE_COLUMNS: Final[tuple[str, ...]] = ("gene_symbol", "hpo_id")

#: The ``aspect`` values kept when building an information-content corpus.
#:
#: ``P`` is the phenotypic-abnormality subontology. ``I`` (inheritance), ``C``
#: (clinical course), ``M`` (modifier) and ``H`` (past medical history) are
#: separate subontologies describing the *annotation*, not the patient's
#: phenotype. Mixing them in dilutes the frequency denominator with terms no
#: clinical profile in this pipeline will ever contain, which deflates the
#: information content of every real phenotype term.
PHENOTYPIC_ABNORMALITY_ASPECT: Final[str] = "P"

#: ``qualifier`` value marking an annotation as an explicit NOT.
#:
#: These rows say "this disease does *not* present this feature". They are real
#: information, but they are not evidence that the term occurs, so counting them
#: as occurrences would inflate the frequency — and therefore deflate the
#: information content — of exactly the terms curators bothered to negate.
NEGATION_QUALIFIER: Final[str] = "NOT"


class CorpusKind(StrEnum):
    """What kind of entity the corpus counts.

    Recorded rather than inferred, because the choice changes every information
    content value and therefore every score, and provenance has to name it.
    """

    DISEASE = "disease"
    """One entity per disease (``phenotype.hpoa``). The Phenomizer convention."""

    GENE = "gene"
    """One entity per gene (``genes_to_phenotype.txt``)."""


@dataclass(frozen=True, slots=True)
class CorpusStats:
    """Counts describing what a corpus load kept and what it dropped.

    Every one of these is a number a reader needs in order to judge the
    information content that comes out. A loader that silently drops a third of
    its rows produces a plausible-looking IC table and a wrong score.

    **The invariant, which is the point of this class:**
    ``rows_kept + rows_negated + rows_wrong_aspect + rows_unresolvable == rows_read``.
    It holds exactly, with no residual "dropped for other reasons" term, because
    every remaining way a row could fail now raises instead (see
    :func:`_iter_tsv` and :func:`_require_field`). If the identity ever stops
    holding, a row is being discarded somewhere without being accounted for, and
    the information content this corpus produces cannot be trusted.
    :meth:`is_balanced` states it in code so a test can assert it directly.
    """

    rows_read: int
    rows_kept: int
    rows_negated: int
    """Rows skipped because ``qualifier`` was ``NOT``."""
    rows_wrong_aspect: int
    """Rows skipped because they were not in the phenotypic-abnormality subontology."""
    rows_unresolvable: int
    """Rows whose HPO term is absent from the loaded release (obsoleted without a
    successor, or from a newer release than the ontology)."""
    entity_count: int

    @property
    def is_balanced(self) -> bool:
        """True when every row read was accounted for by exactly one outcome."""
        return (
            self.rows_kept + self.rows_negated + self.rows_wrong_aspect + self.rows_unresolvable
            == self.rows_read
        )


class AnnotationCorpus:
    """Entity -> HPO terms, resolved against one ontology release.

    Immutable after construction and deterministically ordered throughout. Terms
    are stored already resolved through :meth:`HpoOntology.resolve`, so an
    ``alt_id`` from a historical merge and its primary identifier count as one
    term rather than two.
    """

    __slots__ = ("_entities", "_kind", "_source_name", "_source_sha256", "_stats", "_unresolved")

    def __init__(
        self,
        entities: Mapping[str, Sequence[str]],
        *,
        kind: CorpusKind,
        stats: CorpusStats,
        source_name: str,
        source_sha256: str | None = None,
        unresolved_terms: Sequence[str] = (),
    ) -> None:
        self._entities: dict[str, tuple[str, ...]] = {
            entity: tuple(sorted(set(terms))) for entity, terms in sorted(entities.items())
        }
        self._kind: CorpusKind = kind
        self._stats: CorpusStats = stats
        self._source_name: str = source_name
        self._source_sha256: str | None = source_sha256
        self._unresolved: tuple[str, ...] = tuple(sorted(set(unresolved_terms)))

    # -- construction -------------------------------------------------------

    @classmethod
    def from_hpoa(
        cls,
        path: Path,
        *,
        ontology: HpoOntology,
        expected_sha256: str | None = None,
    ) -> AnnotationCorpus:
        """Load ``phenotype.hpoa`` as a per-disease corpus.

        This is the corpus Köhler et al. (2009) use, and the default for this
        package: information content answers "how surprising is this finding across
        the space of rare diseases?", which is the question a differential
        diagnosis asks.

        ``NOT``-qualified rows and non-``P`` aspects are excluded; both exclusions
        are counted in :class:`CorpusStats` rather than being silent. Every other
        row is either kept or raises — see :func:`_iter_tsv` for why nothing here
        may be dropped quietly.
        """
        digest = _verify(path, expected_sha256=expected_sha256, label="HPO disease annotations")
        entities: dict[str, list[str]] = {}
        unresolved: set[str] = set()
        rows = kept = negated = wrong_aspect = unresolvable = 0

        for lineno, fields in _iter_tsv(path, required=HPOA_COLUMNS):
            rows += 1
            # Validated before the qualifier/aspect filters so that a malformed row
            # cannot hide behind a bucket that would have discarded it anyway.
            entity = _require_field(
                fields["database_id"], column="database_id", path=path, lineno=lineno
            )
            if fields["qualifier"].strip().upper() == NEGATION_QUALIFIER:
                negated += 1
                continue
            if fields["aspect"].strip().upper() != PHENOTYPIC_ABNORMALITY_ASPECT:
                wrong_aspect += 1
                continue
            raw_term = _require_field(fields["hpo_id"], column="hpo_id", path=path, lineno=lineno)
            resolved = ontology.resolve(raw_term)
            if resolved is None:
                unresolvable += 1
                unresolved.add(raw_term.upper().replace("HP_", "HP:"))
                continue
            entities.setdefault(entity, []).append(resolved)
            kept += 1

        stats = CorpusStats(
            rows_read=rows,
            rows_kept=kept,
            rows_negated=negated,
            rows_wrong_aspect=wrong_aspect,
            rows_unresolvable=unresolvable,
            entity_count=len(entities),
        )
        return cls(
            entities,
            kind=CorpusKind.DISEASE,
            stats=stats,
            source_name=path.name,
            source_sha256=digest,
            unresolved_terms=sorted(unresolved),
        )

    @classmethod
    def from_genes_to_phenotype(
        cls,
        path: Path,
        *,
        ontology: HpoOntology,
        expected_sha256: str | None = None,
    ) -> AnnotationCorpus:
        """Load ``genes_to_phenotype.txt`` as a per-gene corpus.

        Offered alongside :meth:`from_hpoa` because this pipeline scores *genes*,
        so a gene-frequency prior is arguably the better-matched denominator. It is
        not the default: the file is derived from the disease annotations, so it
        weights a term by how many genes reach it, which over-counts terms in
        locus-heterogeneous conditions. Which corpus was used is recorded in
        provenance either way — the two give different information content and
        therefore different scores, and that must never be ambiguous.
        """
        digest = _verify(path, expected_sha256=expected_sha256, label="HPO gene annotations")
        entities: dict[str, list[str]] = {}
        unresolved: set[str] = set()
        rows = kept = unresolvable = 0

        for lineno, fields in _iter_tsv(path, required=GENES_TO_PHENOTYPE_COLUMNS):
            rows += 1
            entity = _require_field(
                fields["gene_symbol"], column="gene_symbol", path=path, lineno=lineno
            )
            raw_term = _require_field(fields["hpo_id"], column="hpo_id", path=path, lineno=lineno)
            resolved = ontology.resolve(raw_term)
            if resolved is None:
                unresolvable += 1
                unresolved.add(raw_term.upper().replace("HP_", "HP:"))
                continue
            entities.setdefault(entity, []).append(resolved)
            kept += 1

        stats = CorpusStats(
            rows_read=rows,
            rows_kept=kept,
            rows_negated=0,
            rows_wrong_aspect=0,
            rows_unresolvable=unresolvable,
            entity_count=len(entities),
        )
        return cls(
            entities,
            kind=CorpusKind.GENE,
            stats=stats,
            source_name=path.name,
            source_sha256=digest,
            unresolved_terms=sorted(unresolved),
        )

    # -- accessors ----------------------------------------------------------

    @property
    def kind(self) -> CorpusKind:
        return self._kind

    @property
    def stats(self) -> CorpusStats:
        return self._stats

    @property
    def source_name(self) -> str:
        return self._source_name

    @property
    def source_sha256(self) -> str | None:
        return self._source_sha256

    @property
    def entity_ids(self) -> tuple[str, ...]:
        """Every entity identifier, sorted (GP-30)."""
        return tuple(self._entities)

    @property
    def unresolved_terms(self) -> tuple[str, ...]:
        """Terms the ontology release could not resolve, sorted.

        Surfaced rather than swallowed: a large value here means the annotation
        release and the ontology release are mismatched, which quietly biases every
        information content value.
        """
        return self._unresolved

    def terms_for(self, entity_id: str) -> tuple[str, ...]:
        return self._entities.get(entity_id, ())

    def __len__(self) -> int:
        return len(self._entities)

    def __repr__(self) -> str:
        return (
            f"AnnotationCorpus(kind={self._kind.value!r}, source={self._source_name!r}, "
            f"entities={len(self._entities)})"
        )


class InformationContent:
    """Corpus-derived specificity for every term in an ontology release.

    ``IC(t) = -ln( n(t) / N )`` where ``n(t)`` is the number of annotated entities
    whose annotation set, closed **upward** under ``is_a``, contains ``t``, and
    ``N`` is the number of entities that contributed at least one resolvable term.

    Two properties follow from that definition and both matter downstream:

    * ``IC(root) == 0``. Every entity reaches the root, so "the patient has a
      phenotypic abnormality" carries no information. A similarity measure whose
      only common ancestor is the root therefore scores zero — which is a *computed*
      zero meaning "these terms share nothing", not a missing value.
    * ``IC`` is undefined, not zero, for a term no entity reaches. :meth:`ic`
      returns ``None`` there. This is GP-14 in numeric form: "no annotation data
      for this term" is a different fact from "this term is maximally general", and
      collapsing the two would let a curation gap read as a real finding.
    """

    __slots__ = ("_corpus_kind", "_corpus_source", "_counts", "_max_ic", "_ontology", "_total")

    def __init__(
        self,
        counts: Mapping[str, int],
        *,
        total: int,
        ontology: HpoOntology,
        corpus_kind: CorpusKind,
        corpus_source: str,
    ) -> None:
        if total <= 0:
            msg = (
                "Cannot build information content from an empty annotation corpus: "
                "every term would be equally and infinitely specific. Check that the "
                "annotation file and the ontology release refer to the same HPO version."
            )
            raise IngestionError(msg)
        self._counts: dict[str, int] = dict(sorted(counts.items()))
        self._total: int = total
        self._ontology: HpoOntology = ontology
        self._corpus_kind: CorpusKind = corpus_kind
        self._corpus_source: str = corpus_source
        self._max_ic: float = max(
            (-math.log(count / total) for count in self._counts.values()), default=0.0
        )

    @classmethod
    def from_corpus(cls, corpus: AnnotationCorpus, *, ontology: HpoOntology) -> InformationContent:
        """Count annotations under the true-path rule and derive IC.

        One pass over the corpus. For each entity the union of the reflexive
        ancestor closures of its terms is computed once and each member counted
        once, so an entity annotated with both a term and its parent contributes
        one to the parent, not two. Double counting there would make heavily
        annotated diseases dominate the frequency of every internal node.
        """
        counts: dict[str, int] = {}
        contributing = 0
        for entity_id in corpus.entity_ids:
            terms = corpus.terms_for(entity_id)
            closure: set[str] = set()
            for term in terms:
                closure.update(ontology.ancestor_closure(term))
            if not closure:
                continue
            contributing += 1
            for term in closure:
                counts[term] = counts.get(term, 0) + 1
        return cls(
            counts,
            total=contributing,
            ontology=ontology,
            corpus_kind=corpus.kind,
            corpus_source=corpus.source_name,
        )

    # -- accessors ----------------------------------------------------------

    def ic(self, hpo_id: str) -> float | None:
        """Information content in nats, or ``None`` when the corpus never reaches the term.

        ``None`` is not an error and must not be replaced by ``0.0`` at the call
        site: zero is the value of the *root*, the least informative term there is,
        and imputing it for a term nobody annotated would make an unstudied finding
        look like a universal one.
        """
        resolved = self._ontology.resolve(hpo_id)
        if resolved is None:
            return None
        count = self._counts.get(resolved)
        if count is None or count <= 0:
            return None
        # `-math.log(1.0)` is `-0.0`, which serialises as "-0.0" and would make an
        # artifact differ from one built by a path that produced `0.0` (GP-30).
        value = -math.log(count / self._total)
        return value if value != 0.0 else 0.0

    def normalised_ic(self, hpo_id: str) -> float | None:
        """:meth:`ic` rescaled to ``[0, 1]`` by the corpus maximum, or ``None``.

        Used where a bounded weight is needed. The rescaling is a presentation
        choice, not a modelling one: it preserves the ordering exactly and the
        divisor is recorded in provenance.
        """
        value = self.ic(hpo_id)
        if value is None:
            return None
        if self._max_ic <= 0.0:
            return 0.0
        return value / self._max_ic

    def annotation_count(self, hpo_id: str) -> int:
        """Entities reaching this term under the true-path rule. ``0`` if none."""
        resolved = self._ontology.resolve(hpo_id)
        if resolved is None:
            return 0
        return self._counts.get(resolved, 0)

    @property
    def max_ic(self) -> float:
        return self._max_ic

    @property
    def entity_count(self) -> int:
        """``N``: entities that contributed at least one resolvable term."""
        return self._total

    @property
    def covered_term_count(self) -> int:
        """Terms with at least one annotation. Terms outside this have no IC."""
        return len(self._counts)

    @property
    def corpus_kind(self) -> CorpusKind:
        return self._corpus_kind

    @property
    def corpus_source(self) -> str:
        return self._corpus_source

    @property
    def ontology(self) -> HpoOntology:
        return self._ontology

    def provenance(self) -> Mapping[str, str]:
        """Sorted, string-valued description of how this IC table was derived."""
        return dict(
            sorted(
                {
                    "ic_corpus_kind": self._corpus_kind.value,
                    "ic_corpus_source": self._corpus_source,
                    "ic_entity_count": str(self._total),
                    "ic_covered_terms": str(len(self._counts)),
                    "ic_max": f"{self._max_ic:.6f}",
                }.items()
            )
        )

    def __repr__(self) -> str:
        return (
            f"InformationContent(corpus={self._corpus_source!r}, entities={self._total}, "
            f"terms={len(self._counts)}, max_ic={self._max_ic:.3f})"
        )


# ---------------------------------------------------------------------------
# Tabular reading
# ---------------------------------------------------------------------------


def _verify(path: Path, *, expected_sha256: str | None, label: str) -> str:
    """Existence and digest check, returning the computed digest for provenance."""
    if not path.is_file():
        msg = (
            f"{label} file not found: {path.as_posix()}. This stage reads local files "
            "only (PRIV-05); the composition root supplies the path."
        )
        raise IngestionError(msg)
    digest = hash_file(path)
    if expected_sha256 is not None and digest != expected_sha256.strip().lower():
        msg = (
            f"{label} {path.name} has sha256 {digest}, but {expected_sha256.strip()} was "
            "pinned. A different annotation release changes every information content "
            "value and therefore every phenotype score; re-pin deliberately."
        )
        raise IngestionError(msg)
    return digest


def _require_field(value: str, *, column: str, path: Path, lineno: int) -> str:
    """Return a stripped mandatory field, or raise naming the file and line only.

    An empty entity identifier used to be skipped with a bare ``continue``, which
    put the row in no bucket at all: it was counted in ``rows_read`` and in nothing
    else, so the statistics silently stopped adding up. Failing closed is the same
    trade :func:`_iter_tsv` makes, for the same reason.

    The message carries the column name, the file name and the line number and
    **never the field's value** (PRIV-09/GP-41): these files are public reference
    data today, but this reader is the one a real annotation export would go
    through, and an error string is a disclosure vector.
    """
    token = value.strip()
    if not token:
        msg = (
            f"{path.name} line {lineno}: required column {column!r} is empty. Refusing to "
            "skip the row: a row that belongs to no entity is counted as read and as "
            "nothing else, which makes the parse statistics stop adding up while the "
            "information-content table quietly changes."
        )
        raise IngestionError(msg)
    return token


def _iter_tsv(path: Path, *, required: Sequence[str]) -> Iterable[tuple[int, dict[str, str]]]:
    """Stream a ``#``-commented TSV, yielding ``(line_number, row)`` pairs.

    A streaming generator rather than :func:`mva.phenotype.hpo.read_tsv_rows`
    because these files are 20-35 MB and 300k rows: materialising every row as a
    dict first costs several hundred megabytes for no benefit. Only the columns in
    ``required`` are retained, which keeps the per-row dict small.

    **A row whose field count differs from the header is rejected, not skipped.**
    This reader used to drop short rows with a bare ``continue``, before the
    caller's ``rows_read`` counter had even incremented. That is the worst shape a
    data bug can take here: information content is the *denominator* of every
    similarity score, so a lost annotation shifts the IC of the affected term and
    of every ancestor of it, moving every phenotype score in the run — while
    :class:`CorpusStats` still reported a clean parse. An unfalsifiable change to
    the science is worse than a loud failure, so this fails closed.

    Padding or truncating instead would be worse again: a shifted ``aspect`` or
    ``qualifier`` column silently reclassifies rows rather than losing them.

    The error names the file and the line number and nothing else — no field
    values, no row content (PRIV-09/GP-41).
    """
    with path.open("r", encoding="utf-8-sig") as handle:
        header: list[str] | None = None
        wanted: list[tuple[str, int]] = []
        for lineno, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if header is None:
                header = [name.strip() for name in fields]
                missing = [name for name in required if name not in header]
                if missing:
                    msg = (
                        f"{path.name} is missing required column(s): {', '.join(missing)}. "
                        f"Header found: {', '.join(header)}."
                    )
                    raise IngestionError(msg)
                wanted = [(name, header.index(name)) for name in required]
                continue
            if len(fields) != len(header):
                msg = (
                    f"{path.name} line {lineno}: expected {len(header)} tab-separated fields "
                    f"(from the header), found {len(fields)}. Refusing to skip or pad the row: "
                    "information content is the denominator of every phenotype similarity "
                    "score, so a silently dropped annotation moves every score in the run "
                    "while the parse statistics still report a clean load. Repair the file or "
                    "re-download the annotation release."
                )
                raise IngestionError(msg)
            yield lineno, {name: fields[position] for name, position in wanted}
        if header is None:
            msg = f"{path.name} contains no header row (only comments or blank lines)."
            raise IngestionError(msg)
