"""Information content and the published similarity measures.

Two properties dominate here.

**IC comes from the annotation corpus, not from graph depth.** The tests below
demonstrate the difference rather than asserting it in a comment: the fixture
contains terms that sit at the same depth with very different information content,
and terms whose depth is not even well defined because they reach the root by
paths of different lengths.

**``None`` is not zero.** A term the corpus never reaches has no information
content, and a pair that cannot be compared has no similarity. Both return
``None``. Zero is a real answer meaning "these terms share only the uninformative
root", and a scorer that cannot tell the two apart turns a curation gap into
evidence.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from mva.errors import IngestionError
from mva.phenotype.corpus import AnnotationCorpus, CorpusKind, CorpusStats, InformationContent
from mva.phenotype.ontology import HpoOntology
from mva.phenotype.similarity import (
    MEASURE_CITATIONS,
    SimilarityMeasure,
    TermSimilarity,
)

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
HPO_FIXTURES = REPO_ROOT / "tests" / "fixtures" / "hpo"
SUBSET_OBO = HPO_FIXTURES / "hp_subset.obo"
SUBSET_HPOA = HPO_FIXTURES / "phenotype_subset.hpoa"
SUBSET_G2P = HPO_FIXTURES / "genes_to_phenotype_subset.txt"

ROOT = "HP:0000001"
ABNORMALITY = "HP:0000118"
EYE = "HP:0000478"
LENS = "HP:0000517"
CATARACT = "HP:0000518"
MICROCEPHALY = "HP:0000252"
IUGR = "HP:0001511"
SEIZURE = "HP:0001250"
DANDY_WALKER = "HP:0001305"
PREMATURE_CHROMATID_SEPARATION = "HP:0200024"


@pytest.fixture(scope="module")
def ontology() -> HpoOntology:
    return HpoOntology.from_obo(SUBSET_OBO)


@pytest.fixture(scope="module")
def corpus(ontology: HpoOntology) -> AnnotationCorpus:
    return AnnotationCorpus.from_hpoa(SUBSET_HPOA, ontology=ontology)


@pytest.fixture(scope="module")
def information_content(corpus: AnnotationCorpus, ontology: HpoOntology) -> InformationContent:
    return InformationContent.from_corpus(corpus, ontology=ontology)


@pytest.fixture(scope="module")
def lin(ontology: HpoOntology, information_content: InformationContent) -> TermSimilarity:
    return TermSimilarity(ontology=ontology, information_content=information_content)


# ---------------------------------------------------------------------------
# Corpus loading
# ---------------------------------------------------------------------------


def test_corpus_reports_what_it_dropped(corpus: AnnotationCorpus) -> None:
    """A loader that silently drops rows produces a plausible, wrong IC table."""
    stats = corpus.stats
    assert stats.rows_read > 0
    assert stats.rows_kept > 0
    assert stats.is_balanced, "every row read must be accounted for by exactly one outcome"
    assert stats.entity_count == len(corpus)


def test_negated_annotations_are_not_counted_as_occurrences(
    ontology: HpoOntology, tmp_path: Path
) -> None:
    """A ``NOT`` row says the disease lacks the feature; counting it inflates frequency."""
    path = tmp_path / "tiny.hpoa"
    path.write_text(
        "#version: test\n"
        "database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\t"
        "frequency\tsex\tmodifier\taspect\tbiocuration\n"
        f"OMIM:1\tA\t\t{CATARACT}\t\t\t\t\t\t\tP\t\n"
        f"OMIM:2\tB\tNOT\t{CATARACT}\t\t\t\t\t\t\tP\t\n"
        f"OMIM:3\tC\t\t{MICROCEPHALY}\t\t\t\t\t\t\tI\t\n",
        encoding="utf-8",
    )
    built = AnnotationCorpus.from_hpoa(path, ontology=ontology)
    assert built.stats.rows_negated == 1
    assert built.stats.rows_wrong_aspect == 1
    assert built.entity_ids == ("OMIM:1",)


def test_missing_column_is_rejected_by_name(ontology: HpoOntology, tmp_path: Path) -> None:
    path = tmp_path / "broken.hpoa"
    path.write_text("database_id\thpo_id\n", encoding="utf-8")
    with pytest.raises(IngestionError, match="qualifier"):
        AnnotationCorpus.from_hpoa(path, ontology=ontology)


#: A well-formed 12-column ``phenotype.hpoa`` data row, for the malformed-row tests.
_HPOA_HEADER = (
    "database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\t"
    "frequency\tsex\tmodifier\taspect\tbiocuration\n"
)
_HPOA_GOOD_ROW = f"OMIM:1\tA\t\t{CATARACT}\t\t\t\t\t\t\tP\t\n"


def _hpoa_with(tmp_path: Path, row: str, *, name: str) -> Path:
    path = tmp_path / name
    path.write_text("#version: unit-test\n" + _HPOA_HEADER + _HPOA_GOOD_ROW + row, encoding="utf-8")
    return path


@pytest.mark.parametrize(
    ("case", "row"),
    [
        # One field short: the field count no longer matches the header, so every
        # column after the gap means something different from what it says.
        ("short", f"OMIM:2\tB\t\t{MICROCEPHALY}\t\t\t\t\t\tP\t\n"),
        # One field long: the same problem with the opposite sign.
        ("overlong", f"OMIM:3\tC\t\t{SEIZURE}\t\t\t\t\t\t\tP\t\tEXTRA\n"),
    ],
)
def test_malformed_data_row_is_never_silently_dropped(
    ontology: HpoOntology, tmp_path: Path, case: str, row: str
) -> None:
    """A row the reader cannot trust must fail the load, not vanish from it.

    This reader used to skip a short row with a bare ``continue``, *before* the
    caller's ``rows_read`` counter incremented — so the annotation disappeared and
    :class:`CorpusStats` still reported a clean parse. That is the worst shape a
    data bug can take here: information content is the denominator of every
    similarity score, so a lost annotation shifts the IC of that term and of every
    ancestor of it, moving every phenotype score in the run with nothing anywhere
    recording that it happened.

    The pinned release parses cleanly today, which is exactly why this needs a test
    rather than an inspection: nothing in the current data would reveal a
    regression here.
    """
    path = _hpoa_with(tmp_path, row, name=f"{case}.hpoa")
    with pytest.raises(IngestionError) as raised:
        AnnotationCorpus.from_hpoa(path, ontology=ontology)
    message = str(raised.value)
    assert "line 4" in message, "the error must name the offending line"
    assert path.name in message, "the error must name the file"


def test_row_with_no_entity_identifier_is_rejected(ontology: HpoOntology, tmp_path: Path) -> None:
    """A row belonging to no disease was counted as read and as nothing else.

    It keeps the header's field count, so the shape check above cannot catch it; it
    simply hit a bare ``continue`` and broke the stats identity.
    """
    path = _hpoa_with(tmp_path, f"\t\t\t{MICROCEPHALY}\t\t\t\t\t\t\tP\t\n", name="no_entity.hpoa")
    with pytest.raises(IngestionError, match="database_id"):
        AnnotationCorpus.from_hpoa(path, ontology=ontology)


def test_malformed_row_errors_carry_no_field_content(ontology: HpoOntology, tmp_path: Path) -> None:
    """GP-41/PRIV-09: the message names the file and the line, never the row.

    These files are public reference data today, but this is the reader a real
    annotation export would pass through, and an exception string is a disclosure
    vector.
    """
    path = _hpoa_with(
        tmp_path,
        "OMIM:SECRET\tPATIENT NARRATIVE\t\tHP:0000252\t\t\t\t\t\tP\t\n",
        name="leak.hpoa",
    )
    with pytest.raises(IngestionError) as raised:
        AnnotationCorpus.from_hpoa(path, ontology=ontology)
    message = str(raised.value)
    for secret in ("OMIM:SECRET", "PATIENT NARRATIVE", "HP:0000252"):
        assert secret not in message, f"error message leaked {secret!r}"


def test_gene_corpus_rejects_malformed_rows_too(ontology: HpoOntology, tmp_path: Path) -> None:
    """The second corpus loader shares the reader and therefore the guarantee."""
    path = tmp_path / "genes.txt"
    path.write_text(
        "ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id\n"
        f"1\tGENEA\t{CATARACT}\tx\t-\tOMIM:1\n"
        f"2\tGENEB\t{MICROCEPHALY}\tx\t-\n",
        encoding="utf-8",
    )
    with pytest.raises(IngestionError, match="line 3"):
        AnnotationCorpus.from_genes_to_phenotype(path, ontology=ontology)


@pytest.mark.parametrize("fixture_name", ["phenotype_subset.hpoa", "genes_to_phenotype_subset.txt"])
def test_committed_fixtures_account_for_every_row(ontology: HpoOntology, fixture_name: str) -> None:
    """The stats identity holds exactly, with no unexplained residue.

    ``rows_kept + negated + wrong_aspect + unresolvable == rows_read`` is the whole
    point of :class:`CorpusStats`. If it ever fails, a row is being discarded
    somewhere without being accounted for, and the information content this corpus
    produces cannot be trusted.
    """
    path = HPO_FIXTURES / fixture_name
    corpus = (
        AnnotationCorpus.from_hpoa(path, ontology=ontology)
        if fixture_name.endswith(".hpoa")
        else AnnotationCorpus.from_genes_to_phenotype(path, ontology=ontology)
    )
    assert corpus.stats.is_balanced
    assert corpus.stats.rows_read > 0


def test_gene_corpus_is_a_different_corpus_and_says_so(ontology: HpoOntology) -> None:
    """The choice of corpus changes every score, so it is recorded, not inferred."""
    genes = AnnotationCorpus.from_genes_to_phenotype(SUBSET_G2P, ontology=ontology)
    assert genes.kind is CorpusKind.GENE
    assert len(genes) > 0
    gene_ic = InformationContent.from_corpus(genes, ontology=ontology)
    assert gene_ic.corpus_kind is CorpusKind.GENE
    assert gene_ic.provenance()["ic_corpus_source"] == SUBSET_G2P.name


def test_empty_corpus_is_refused(ontology: HpoOntology) -> None:
    """Zero entities would make every term infinitely and equally specific."""
    empty = AnnotationCorpus(
        {},
        kind=CorpusKind.DISEASE,
        stats=CorpusStats(
            rows_read=0,
            rows_kept=0,
            rows_negated=0,
            rows_wrong_aspect=0,
            rows_unresolvable=0,
            entity_count=0,
        ),
        source_name="empty",
    )
    with pytest.raises(IngestionError, match="empty annotation corpus"):
        InformationContent.from_corpus(empty, ontology=ontology)


# ---------------------------------------------------------------------------
# Information content
# ---------------------------------------------------------------------------


def test_true_path_rule_gives_an_ancestor_at_least_its_descendant_count(
    information_content: InformationContent, ontology: HpoOntology
) -> None:
    """An entity annotated with a term is an instance of every ancestor of that term.

    Without this, every internal node would look rarer than its own children and the
    whole specificity ordering would invert.
    """
    for term in (CATARACT, MICROCEPHALY, DANDY_WALKER):
        child_count = information_content.annotation_count(term)
        assert child_count > 0
        for ancestor in ontology.ancestor_closure(term):
            assert information_content.annotation_count(ancestor) >= child_count


def test_information_content_decreases_toward_the_root(
    information_content: InformationContent, ontology: HpoOntology
) -> None:
    """IC is monotone non-increasing along any upward path, by construction."""
    for ancestor in ontology.ancestor_closure(CATARACT):
        parent_ic = information_content.ic(ancestor)
        child_ic = information_content.ic(CATARACT)
        assert parent_ic is not None and child_ic is not None
        assert parent_ic <= child_ic + 1e-12


def test_root_has_zero_information_content(information_content: InformationContent) -> None:
    """ "The patient has a phenotypic abnormality" tells a differential nothing."""
    assert information_content.ic(ROOT) == 0.0
    assert information_content.ic(ABNORMALITY) == 0.0
    assert math.copysign(1.0, information_content.ic(ROOT) or 0.0) > 0, "no negative zero (GP-30)"


def test_information_content_is_not_graph_depth(
    information_content: InformationContent, ontology: HpoOntology
) -> None:
    """The claim the module makes, demonstrated on real terms rather than asserted.

    Two terms one step below the same parent have the same depth. If depth were a
    usable proxy for specificity they would have to be equally specific; the corpus
    says otherwise.
    """
    siblings = [
        term for term in ontology.children(LENS) if information_content.ic(term) is not None
    ]
    values = sorted({round(information_content.ic(term) or 0.0, 6) for term in siblings})
    if len(values) < 2:  # pragma: no cover - fixture guard
        pytest.skip("fixture slice has too few annotated siblings to demonstrate this")
    assert values[0] != values[-1], "same depth, same IC: depth would have been enough"


def test_unannotated_term_has_no_information_content_rather_than_zero(
    information_content: InformationContent, ontology: HpoOntology
) -> None:
    """``None`` is the honest answer. Zero is the value of the root."""
    unannotated = [
        term for term in ontology.term_ids if information_content.annotation_count(term) == 0
    ]
    assert unannotated, "fixture slice has no unannotated term to exercise this path"
    for term in unannotated[:20]:
        assert information_content.ic(term) is None
        assert information_content.normalised_ic(term) is None


def test_normalised_ic_is_bounded(
    information_content: InformationContent, ontology: HpoOntology
) -> None:
    for term in ontology.term_ids:
        value = information_content.normalised_ic(term)
        if value is not None:
            assert 0.0 <= value <= 1.0


def test_ic_matches_the_resnik_definition(information_content: InformationContent) -> None:
    """-ln(n(t)/N), computed independently from the reported counts."""
    count = information_content.annotation_count(CATARACT)
    expected = -math.log(count / information_content.entity_count)
    actual = information_content.ic(CATARACT)
    assert actual is not None
    assert actual == pytest.approx(expected)


# ---------------------------------------------------------------------------
# MICA
# ---------------------------------------------------------------------------


def test_mica_of_a_term_with_itself_is_the_term(lin: TermSimilarity) -> None:
    assert lin.most_informative_common_ancestor(CATARACT, CATARACT) == CATARACT


def test_mica_of_parent_and_child_is_the_parent(lin: TermSimilarity) -> None:
    assert lin.most_informative_common_ancestor(CATARACT, LENS) == LENS


def test_mica_is_chosen_by_information_content_not_by_depth(
    lin: TermSimilarity, information_content: InformationContent, ontology: HpoOntology
) -> None:
    """On a DAG the shared ancestors are not totally ordered; IC picks among them."""
    mica = lin.most_informative_common_ancestor(CATARACT, MICROCEPHALY)
    assert mica is not None
    shared = set(ontology.ancestor_closure(CATARACT)) & set(ontology.ancestor_closure(MICROCEPHALY))
    best = max(
        (information_content.ic(term) or 0.0 for term in shared),
        default=0.0,
    )
    assert (information_content.ic(mica) or 0.0) == pytest.approx(best)


def test_mica_tie_breaks_on_the_lowest_identifier(
    ontology: HpoOntology, information_content: InformationContent
) -> None:
    """Equal-IC ties are routine; without a stable rule the MICA changes between runs."""
    similarity = TermSimilarity(ontology=ontology, information_content=information_content)
    repeats = {
        similarity.most_informative_common_ancestor(CATARACT, MICROCEPHALY) for _ in range(25)
    }
    assert len(repeats) == 1


# ---------------------------------------------------------------------------
# Pairwise measures
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("measure", list(SimilarityMeasure))
def test_every_measure_is_bounded_and_self_maximal(
    measure: SimilarityMeasure,
    ontology: HpoOntology,
    information_content: InformationContent,
) -> None:
    similarity = TermSimilarity(
        ontology=ontology, information_content=information_content, measure=measure
    )
    assert similarity.pairwise(CATARACT, CATARACT) == 1.0
    for other in (MICROCEPHALY, SEIZURE, DANDY_WALKER, IUGR):
        value = similarity.pairwise(CATARACT, other)
        assert value is not None
        assert 0.0 <= value <= 1.0
        assert value <= 1.0


@pytest.mark.parametrize("measure", list(SimilarityMeasure))
def test_measures_are_symmetric(
    measure: SimilarityMeasure,
    ontology: HpoOntology,
    information_content: InformationContent,
) -> None:
    similarity = TermSimilarity(
        ontology=ontology, information_content=information_content, measure=measure
    )
    for a, b in ((CATARACT, MICROCEPHALY), (SEIZURE, DANDY_WALKER), (IUGR, CATARACT)):
        assert similarity.pairwise(a, b) == similarity.pairwise(b, a)


def test_lin_matches_its_published_formula(
    lin: TermSimilarity, information_content: InformationContent
) -> None:
    """2*IC(MICA) / (IC(a) + IC(b)), Lin (1998). Recomputed here from the parts."""
    mica = lin.most_informative_common_ancestor(CATARACT, MICROCEPHALY)
    assert mica is not None
    mica_ic = information_content.ic(mica)
    a_ic = information_content.ic(CATARACT)
    b_ic = information_content.ic(MICROCEPHALY)
    assert mica_ic is not None and a_ic is not None and b_ic is not None
    expected = (2.0 * mica_ic) / (a_ic + b_ic)
    assert lin.pairwise(CATARACT, MICROCEPHALY) == pytest.approx(expected)


def test_a_closer_term_scores_higher(lin: TermSimilarity) -> None:
    """Parent/child must beat cousins-across-the-ontology. This is TD-04's whole point."""
    parent_child = lin.pairwise(CATARACT, LENS)
    distant = lin.pairwise(CATARACT, MICROCEPHALY)
    assert parent_child is not None and distant is not None
    assert parent_child > distant


def test_sharing_only_the_root_scores_a_computed_zero(lin: TermSimilarity) -> None:
    """Zero here is a result, not a missing value: the terms genuinely share nothing."""
    value = lin.pairwise(ABNORMALITY, ROOT)
    assert value == 0.0


def test_a_term_outside_the_release_is_uncomparable_not_dissimilar(
    lin: TermSimilarity,
) -> None:
    assert lin.pairwise("HP:9999999", CATARACT) is None
    assert lin.pairwise(CATARACT, "HP:9999999") is None


def test_lin_needs_both_information_contents_and_says_so_by_returning_none(
    ontology: HpoOntology, information_content: InformationContent
) -> None:
    """Substituting a placeholder IC would invent specificity the corpus lacks."""
    unannotated = next(
        term for term in ontology.term_ids if information_content.annotation_count(term) == 0
    )
    lin_measure = TermSimilarity(ontology=ontology, information_content=information_content)
    assert lin_measure.pairwise(unannotated, CATARACT) is None
    resnik = TermSimilarity(
        ontology=ontology,
        information_content=information_content,
        measure=SimilarityMeasure.RESNIK,
    )
    # Resnik needs only the MICA, so it can still place the term.
    assert resnik.pairwise(unannotated, CATARACT) is not None


def test_explain_returns_the_mica_behind_the_number(lin: TermSimilarity) -> None:
    """A bare similarity is not reviewable; the shared ancestor is."""
    value, mica = lin.explain(CATARACT, MICROCEPHALY)
    assert value is not None
    assert mica is not None
    assert mica in lin.ontology.ancestor_closure(CATARACT)
    assert mica in lin.ontology.ancestor_closure(MICROCEPHALY)


def test_each_measure_carries_a_citation() -> None:
    """GP-17: a number in a report must name the published method that produced it."""
    for measure in SimilarityMeasure:
        assert MEASURE_CITATIONS[measure]
        assert any(char.isdigit() for char in MEASURE_CITATIONS[measure]), "needs a year"


# ---------------------------------------------------------------------------
# Best-match aggregation
# ---------------------------------------------------------------------------


def test_best_match_average_of_identical_sets_is_one(lin: TermSimilarity) -> None:
    terms = [CATARACT, MICROCEPHALY, SEIZURE]
    result = lin.best_match_average(terms, terms)
    assert result.value == pytest.approx(1.0)
    assert result.uncomparable == ()


def test_best_match_average_credits_a_parent_target(lin: TermSimilarity) -> None:
    """The TD-04 fix at the aggregation layer: the parent is not a zero."""
    against_parent = lin.best_match_average([CATARACT], [LENS])
    against_stranger = lin.best_match_average([CATARACT], [SEIZURE])
    assert against_parent.value is not None and against_stranger.value is not None
    assert against_parent.value > against_stranger.value
    assert against_parent.value > 0.0


def test_adding_a_target_never_lowers_a_best_match(lin: TermSimilarity) -> None:
    """``max`` over targets: more curated terms can only help.

    This is the monotonicity the scorer relies on to guarantee that a NOT_ASSESSED
    association cannot penalise a gene (GP-14).
    """
    base = lin.best_match([CATARACT][0], [SEIZURE])
    widened = lin.best_match([CATARACT][0], [SEIZURE, LENS])
    assert base is not None and widened is not None
    assert widened.similarity >= base.similarity


def test_uncomparable_query_terms_are_reported_not_averaged_as_zero(
    lin: TermSimilarity,
) -> None:
    """Averaging in a zero for a term the corpus cannot place would fabricate a mismatch."""
    result = lin.best_match_average([CATARACT, "HP:9999999"], [CATARACT])
    assert result.uncomparable == ("HP:9999999",)
    assert result.value == pytest.approx(1.0)


def test_an_entirely_uncomparable_query_returns_none_not_zero(lin: TermSimilarity) -> None:
    result = lin.best_match_average(["HP:9999999"], [CATARACT])
    assert result.value is None
    assert result.matches == ()


def test_empty_target_set_is_uncomparable_not_zero(lin: TermSimilarity) -> None:
    """A gene with no curated terms has not been argued against (GP-14)."""
    result = lin.best_match_average([CATARACT], [])
    assert result.value is None
    assert result.uncomparable == (CATARACT,)


def test_symmetric_best_match_average_is_symmetric(lin: TermSimilarity) -> None:
    forward = lin.symmetric_best_match_average([CATARACT, MICROCEPHALY], [LENS, SEIZURE])
    backward = lin.symmetric_best_match_average([LENS, SEIZURE], [CATARACT, MICROCEPHALY])
    assert forward is not None
    assert forward == pytest.approx(backward)


def test_symmetric_average_lies_between_its_two_directions(lin: TermSimilarity) -> None:
    query = [CATARACT, MICROCEPHALY]
    target = [LENS, SEIZURE, DANDY_WALKER]
    forward = lin.best_match_average(query, target).value
    backward = lin.best_match_average(target, query).value
    symmetric = lin.symmetric_best_match_average(query, target)
    assert forward is not None and backward is not None and symmetric is not None
    assert min(forward, backward) <= symmetric <= max(forward, backward)


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------


def test_results_do_not_depend_on_input_order(lin: TermSimilarity) -> None:
    """Query and target sequences are sorted internally, so caller order is irrelevant."""
    forward = lin.best_match_average([CATARACT, MICROCEPHALY], [SEIZURE, LENS])
    reverse = lin.best_match_average([MICROCEPHALY, CATARACT], [LENS, SEIZURE])
    assert forward.value == reverse.value
    assert forward.matches == reverse.matches


def test_duplicate_terms_do_not_change_an_average(lin: TermSimilarity) -> None:
    once = lin.best_match_average([CATARACT, MICROCEPHALY], [LENS])
    twice = lin.best_match_average([CATARACT, CATARACT, MICROCEPHALY], [LENS, LENS])
    assert once.value == twice.value


def test_caching_does_not_change_an_answer(
    ontology: HpoOntology, information_content: InformationContent
) -> None:
    warm = TermSimilarity(ontology=ontology, information_content=information_content)
    cold = TermSimilarity(ontology=ontology, information_content=information_content)
    first = warm.pairwise(CATARACT, MICROCEPHALY)
    again = warm.pairwise(CATARACT, MICROCEPHALY)
    assert first == again == cold.pairwise(CATARACT, MICROCEPHALY)


def test_premature_chromatid_separation_is_more_specific_than_seizure(
    information_content: InformationContent,
) -> None:
    """A sanity check against the real corpus, not an arbitrary constant.

    If this ever inverts, the corpus and the ontology release have drifted apart or
    the true-path propagation has broken.
    """
    rare = information_content.ic(PREMATURE_CHROMATID_SEPARATION)
    common = information_content.ic(SEIZURE)
    assert rare is not None and common is not None
    assert rare > common
