"""Published information-content semantic similarity over the HPO DAG.

Nothing here is invented. Three pairwise measures and one aggregation are
implemented, all from the literature, and the one actually used is recorded in
provenance so a reader never has to guess.

Pairwise measures
-----------------
All three are built on the **most informative common ancestor** (MICA): the term
of highest information content that subsumes both arguments. On a DAG the MICA is
not "the lowest common ancestor" — a pair can have several incomparable common
ancestors, and the informative one is not always the deepest. It is selected by
information content, never by depth.

``resnik``
    ``sim(a, b) = IC(MICA(a, b))``. Resnik P. (1995) *Using information content to
    evaluate semantic similarity in a taxonomy*, IJCAI-95, 448-453. Unbounded
    above; rescaled here by the corpus maximum so it can be combined with other
    ``[0, 1]`` components.

``lin``
    ``sim(a, b) = 2*IC(MICA) / (IC(a) + IC(b))``. Lin D. (1998) *An
    information-theoretic definition of similarity*, ICML-98, 296-304. Natively in
    ``[0, 1]`` and normalised by the specificity of the terms being compared, so a
    shared ancestor counts for more between two rare findings than between two
    common ones. **This is the default.**

``jiang_conrath``
    ``d(a, b) = IC(a) + IC(b) - 2*IC(MICA)``, converted to a similarity by
    ``1 / (1 + d)``. Jiang J.J., Conrath D.W. (1997) *Semantic similarity based on
    corpus statistics and lexical taxonomy*, ROCLING X, 19-33; the reciprocal
    conversion follows Lord P.W. et al. (2003), Bioinformatics 19(10):1275-1283.

Aggregation
-----------
:meth:`TermSimilarity.best_match_average` and
:meth:`TermSimilarity.symmetric_best_match_average` implement the best-match
average of Köhler S. et al. (2009) *Clinical diagnostics in human genetics with
semantic similarity searches in ontologies*, Am J Hum Genet 85(4):457-464 — the
Phenomizer measure. Each query term is scored against its best-matching target
term and those best matches are averaged; the symmetric form averages the two
directions.

Determinism (GP-30)
-------------------
Every ``argmax`` in this module breaks ties by the **lowest HPO identifier**, and
every iteration is over a sorted sequence. Two terms with equal information
content are a routine occurrence — any two terms annotated to the same number of
diseases tie exactly — so an unstable argmax here would produce a different MICA,
a different ``method`` string and different artifact bytes on the next run.

``None`` is not zero
--------------------
Every function returns ``None`` where the comparison could not be made — an
unresolvable term, or a term the annotation corpus never reaches — and a real
number where it could. Zero is a *result*: it means "these terms share only the
root". Collapsing the two would make a curation gap indistinguishable from a
mismatch, which is the GP-14 failure this package exists to prevent.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from enum import StrEnum

from mva.phenotype.corpus import InformationContent
from mva.phenotype.ontology import HpoOntology


class SimilarityMeasure(StrEnum):
    """Which published pairwise measure a :class:`TermSimilarity` applies."""

    RESNIK = "resnik"
    """IC of the most informative common ancestor (Resnik 1995), rescaled to [0,1]."""

    LIN = "lin"
    """Resnik normalised by the specificity of both terms (Lin 1998). Default."""

    JIANG_CONRATH = "jiang_conrath"
    """Reciprocal of the Jiang-Conrath distance (Jiang & Conrath 1997)."""


#: Human-readable citation per measure, embedded in evidence ``method`` strings so
#: a report states which published measure produced a number (GP-10, GP-17).
MEASURE_CITATIONS: dict[SimilarityMeasure, str] = {
    SimilarityMeasure.RESNIK: (
        "Resnik (1995), Using information content to evaluate semantic similarity "
        "in a taxonomy, IJCAI-95:448-453"
    ),
    SimilarityMeasure.LIN: (
        "Lin (1998), An information-theoretic definition of similarity, ICML-98:296-304"
    ),
    SimilarityMeasure.JIANG_CONRATH: (
        "Jiang & Conrath (1997), Semantic similarity based on corpus statistics and "
        "lexical taxonomy, ROCLING X:19-33; reciprocal conversion after Lord et al. "
        "(2003), Bioinformatics 19(10):1275-1283"
    ),
}

#: Citation for the aggregation, which is a separate published choice from the
#: pairwise measure and is reported separately.
AGGREGATION_CITATION: str = (
    "Köhler et al. (2009), Clinical diagnostics in human genetics with semantic "
    "similarity searches in ontologies, Am J Hum Genet 85(4):457-464 "
    "(symmetric best-match average, 'Phenomizer')"
)


@dataclass(frozen=True, slots=True)
class BestMatch:
    """One query term and the target term that matched it best.

    Carries the MICA as well as the number, because "these two terms are 0.62
    similar" is not reviewable and "they share *Abnormality of skull size*" is.
    """

    query: str
    target: str
    similarity: float
    mica: str | None
    """The most informative common ancestor. ``None`` when the measure did not need
    one (identical terms) or none carried information content."""


@dataclass(frozen=True, slots=True)
class BestMatchAverage:
    """Result of a one-directional best-match average.

    ``value`` is ``None`` when no query term could be compared at all. That is a
    different state from ``0.0`` and the two must not be merged: ``0.0`` means the
    gene's features and the patient's features share nothing beyond the root, and
    ``None`` means the comparison could not be attempted.
    """

    value: float | None
    matches: tuple[BestMatch, ...]
    """Best match per comparable query term, in query order."""
    uncomparable: tuple[str, ...]
    """Query terms with no computable match, sorted. An information gap, not a zero."""


class TermSimilarity:
    """Pairwise and aggregate semantic similarity for one ontology + corpus pair.

    Explicitly constructed and passed around; there is no module-level instance and
    no cached factory function. The ontology, the annotation corpus and the measure
    together determine every number this object produces, so making it a global
    would make those three invisible at the call site — and would let one case's
    ontology leak into the next.

    The pair cache lives on the instance for the same reason. It is bounded by the
    number of distinct term pairs a run actually compares, which for a clinical
    profile against a gene panel is small.
    """

    __slots__ = ("_ic", "_measure", "_ontology", "_pair_cache")

    def __init__(
        self,
        *,
        ontology: HpoOntology,
        information_content: InformationContent,
        measure: SimilarityMeasure = SimilarityMeasure.LIN,
    ) -> None:
        self._ontology: HpoOntology = ontology
        self._ic: InformationContent = information_content
        self._measure: SimilarityMeasure = measure
        self._pair_cache: dict[tuple[str, str], tuple[float | None, str | None]] = {}

    @property
    def measure(self) -> SimilarityMeasure:
        return self._measure

    @property
    def ontology(self) -> HpoOntology:
        return self._ontology

    @property
    def information_content(self) -> InformationContent:
        return self._ic

    @property
    def citation(self) -> str:
        """Citation for the pairwise measure in use."""
        return MEASURE_CITATIONS[self._measure]

    # -- MICA ---------------------------------------------------------------

    def most_informative_common_ancestor(self, first: str, second: str) -> str | None:
        """The common ancestor of highest information content, or ``None``.

        Intersects the two **reflexive upward** closures, so a term is its own
        ancestor and ``mica(t, t) == t``.

        Selection is by information content. Exact ties are routine — any two terms
        annotated to the same number of diseases have identical IC, which happens
        constantly between a term and a parent that nothing else sits under — so the
        tie-break is spelled out and deterministic (GP-30):

        1. higher information content wins;
        2. on an exact tie, the **more specific** candidate wins, meaning the one
           that has the other among its ancestors. Without this, ``mica(t, t)``
           could return an ancestor of ``t`` rather than ``t``;
        3. otherwise the candidate encountered first in the sorted closure wins,
           which is the lowest HPO identifier.

        Returns ``None`` when the terms share no ancestor the corpus ever reaches.
        On a single-rooted ontology that is unusual, and it means the corpus and
        the ontology release do not match.
        """
        left = self._ontology.ancestor_closure(first)
        if not left:
            return None
        right = frozenset(self._ontology.ancestor_closure(second))
        if not right:
            return None

        best_id: str | None = None
        best_ic = float("-inf")
        for candidate in left:  # already sorted -> first-wins is lowest-id-wins
            if candidate not in right:
                continue
            value = self._ic.ic(candidate)
            if value is None:
                continue
            if best_id is None or value > best_ic:
                best_ic = value
                best_id = candidate
            elif value == best_ic and best_id in self._ontology.ancestor_closure(candidate):
                best_id = candidate
        return best_id

    # -- pairwise -----------------------------------------------------------

    def pairwise(self, first: str, second: str) -> float | None:
        """Similarity in ``[0, 1]`` under the configured measure, or ``None``.

        ``None`` means "not computable from this ontology and corpus": one of the
        terms is not in the release, or the measure needs an information content
        that the corpus does not provide for it. It never means "no similarity".

        Which inputs each measure needs differs, and that difference is
        deliberate rather than smoothed over:

        * ``resnik`` needs only ``IC(MICA)``, so it still yields a value for a term
          the corpus never annotates directly, as long as an annotated ancestor
          exists.
        * ``lin`` and ``jiang_conrath`` divide by ``IC(a) + IC(b)`` and therefore
          need both. Substituting a placeholder IC for a missing one would invent
          a specificity the corpus does not support.
        """
        cached = self._pair_cache.get((first, second))
        if cached is not None:
            return cached[0]
        value, mica = self._compute_pair(first, second)
        self._pair_cache[(first, second)] = (value, mica)
        return value

    def explain(self, first: str, second: str) -> tuple[float | None, str | None]:
        """:meth:`pairwise` together with the MICA that produced it."""
        cached = self._pair_cache.get((first, second))
        if cached is not None:
            return cached
        result = self._compute_pair(first, second)
        self._pair_cache[(first, second)] = result
        return result

    def _compute_pair(self, first: str, second: str) -> tuple[float | None, str | None]:
        left = self._ontology.resolve(first)
        right = self._ontology.resolve(second)
        if left is None or right is None:
            return None, None
        if left == right:
            # A term is maximally similar to itself under every measure here. Stated
            # explicitly because Lin is 0/0 for two identical root-level terms.
            return 1.0, left

        mica = self.most_informative_common_ancestor(left, right)
        if mica is None:
            return None, None
        mica_ic = self._ic.ic(mica)
        if mica_ic is None:
            return None, None

        if self._measure is SimilarityMeasure.RESNIK:
            max_ic = self._ic.max_ic
            if max_ic <= 0.0:
                return 0.0, mica
            return _clamp(mica_ic / max_ic), mica

        left_ic = self._ic.ic(left)
        right_ic = self._ic.ic(right)
        if left_ic is None or right_ic is None:
            return None, None

        if self._measure is SimilarityMeasure.LIN:
            denominator = left_ic + right_ic
            if denominator <= 0.0:
                # Both terms are the root: no information on either side.
                return 0.0, mica
            return _clamp((2.0 * mica_ic) / denominator), mica

        distance = left_ic + right_ic - (2.0 * mica_ic)
        return _clamp(1.0 / (1.0 + max(0.0, distance))), mica

    # -- aggregation --------------------------------------------------------

    def best_match(self, query: str, targets: Sequence[str]) -> BestMatch | None:
        """Highest-scoring target for one query term, or ``None`` if none is comparable.

        Targets are iterated in sorted order and the comparison is strict ``>``, so
        an exact tie resolves to the lowest HPO identifier. Deterministic
        tie-breaking is load-bearing here rather than cosmetic: a gene annotated
        with two sibling terms will tie for many queries.
        """
        best: BestMatch | None = None
        for target in sorted(set(targets)):
            value, mica = self.explain(query, target)
            if value is None:
                continue
            if best is None or value > best.similarity:
                best = BestMatch(query=query, target=target, similarity=value, mica=mica)
        return best

    def best_match_average(
        self, queries: Iterable[str], targets: Sequence[str]
    ) -> BestMatchAverage:
        """Mean over query terms of each one's best match against ``targets``.

        This is ``sim(Q → D)`` in Köhler et al. (2009). Query terms with no
        computable match are **excluded from the mean and reported**, not scored as
        zero: averaging in a zero for a term the corpus cannot place would let a
        gap in the annotation release depress a real match (GP-14).
        """
        ordered = sorted(set(queries))
        target_list = sorted(set(targets))
        matches: list[BestMatch] = []
        uncomparable: list[str] = []
        for query in ordered:
            best = self.best_match(query, target_list) if target_list else None
            if best is None:
                uncomparable.append(query)
            else:
                matches.append(best)
        if not matches:
            return BestMatchAverage(
                value=None, matches=(), uncomparable=tuple(sorted(uncomparable))
            )
        total = sum(match.similarity for match in matches)
        return BestMatchAverage(
            value=total / len(matches),
            matches=tuple(matches),
            uncomparable=tuple(sorted(uncomparable)),
        )

    def symmetric_best_match_average(
        self, first: Sequence[str], second: Sequence[str]
    ) -> float | None:
        """``0.5 * [sim(A -> B) + sim(B -> A)]`` — the Phenomizer symmetric score.

        Reported for audit alongside the pipeline score rather than used as the
        score itself. The ``B → A`` direction averages over *every* term of the
        second set, which when the second set is a gene's curated profile means a
        well-studied gene is penalised for features nobody in this clinic assessed
        — the exact effect GP-14 forbids. :mod:`mva.phenotype.scoring` therefore
        uses the ``A → B`` direction plus a separate, assessment-restricted
        specificity component, and states so.

        ``None`` if either direction is uncomputable.
        """
        forward = self.best_match_average(first, second).value
        backward = self.best_match_average(second, first).value
        if forward is None or backward is None:
            return None
        return (forward + backward) / 2.0


def _clamp(value: float) -> float:
    """Bound a similarity to ``[0, 1]``, absorbing float error at the endpoints."""
    if value < 0.0:
        return 0.0
    if value > 1.0:
        return 1.0
    return value
