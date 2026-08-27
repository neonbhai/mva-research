"""Structural tests for the HPO DAG loader.

These are not "does it parse" tests. Each one pins a property that a plausible
wrong implementation would violate:

* a term has several parents, and a loader that reads one loses most of the graph;
* the same ancestor is reachable by paths of different lengths, so depth is not a
  well-defined property of a term and must never stand in for specificity;
* a cyclic file must terminate rather than recurse to death;
* an obsoleted or aliased identifier must resolve or report itself unresolvable,
  never silently become a term that matches nothing.

The fixture under ``tests/fixtures/hpo/`` is a real slice of the published
release, not a hand-drawn graph, so these properties are the ontology's, not the
author's. See ``tests/fixtures/hpo/build_fixture.py``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mva.errors import IngestionError
from mva.phenotype.ontology import HpoOntology, HpoTerm, ontology_provenance

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
HPO_FIXTURES = REPO_ROOT / "tests" / "fixtures" / "hpo"
SUBSET_OBO = HPO_FIXTURES / "hp_subset.obo"

ROOT = "HP:0000001"
MICROCEPHALY = "HP:0000252"
DANDY_WALKER = "HP:0001305"
CATARACT = "HP:0000518"
EYE_ABNORMALITY = "HP:0000478"
LENS_ABNORMALITY = "HP:0000517"
#: Reachable from the root by paths of markedly different length; see the fixture
#: generator's docstring.
UNEVEN_DEPTH_TERM = "HP:0007370"


@pytest.fixture(scope="module")
def ontology() -> HpoOntology:
    return HpoOntology.from_obo(SUBSET_OBO)


def _tiny(*terms: HpoTerm) -> HpoOntology:
    """A hand-built graph, for cases the real release cannot contain (e.g. a cycle)."""
    return HpoOntology(terms, data_version="hp/releases/0000-00-00")


def _term(
    hpo_id: str,
    *,
    parents: tuple[str, ...] = (),
    name: str = "term",
    alt_ids: tuple[str, ...] = (),
    obsolete: bool = False,
    replaced_by: str | None = None,
) -> HpoTerm:
    return HpoTerm(
        hpo_id=hpo_id,
        name=name,
        parents=parents,
        alt_ids=alt_ids,
        is_obsolete=obsolete,
        replaced_by=replaced_by,
    )


# ---------------------------------------------------------------------------
# Release identity and provenance
# ---------------------------------------------------------------------------


def test_release_version_is_read_from_the_file_not_assumed(ontology: HpoOntology) -> None:
    """The release string comes from the OBO header and flows into provenance.

    A phenotype score is only interpretable against a named ontology release, so
    the version is parsed, not configured next to the file where the two can drift.
    """
    assert ontology.data_version.startswith("hp/releases/")
    assert ontology.release == ontology.data_version.rsplit("/", 1)[-1]
    provenance = ontology_provenance(ontology)
    assert provenance["hpo_data_version"] == ontology.data_version
    assert provenance["hpo_sha256"], "the parsed file's digest must reach provenance"
    assert list(provenance) == sorted(provenance), "provenance keys must be sorted (GP-30)"


def test_missing_data_version_is_rejected(tmp_path: Path) -> None:
    """An OBO with no release header is refused rather than labelled 'unknown'."""
    path = tmp_path / "no_version.obo"
    path.write_text("format-version: 1.2\n\n[Term]\nid: HP:0000001\nname: All\n", encoding="utf-8")
    with pytest.raises(IngestionError, match="data-version"):
        HpoOntology.from_obo(path)


def test_pinned_digest_mismatch_is_a_hard_failure(tmp_path: Path) -> None:
    """A swapped ontology release silently changes every score, so the pin is fatal."""
    path = tmp_path / "hp.obo"
    path.write_text(
        "data-version: hp/releases/2020-01-01\n\n[Term]\nid: HP:0000001\nname: All\n",
        encoding="utf-8",
    )
    with pytest.raises(IngestionError, match="sha256"):
        HpoOntology.from_obo(path, expected_sha256="0" * 64)


def test_missing_file_message_names_the_offline_requirement(tmp_path: Path) -> None:
    with pytest.raises(IngestionError, match="PRIV-05"):
        HpoOntology.from_obo(tmp_path / "absent.obo")


# ---------------------------------------------------------------------------
# The DAG is a DAG
# ---------------------------------------------------------------------------


def test_a_real_term_has_several_is_a_parents(ontology: HpoOntology) -> None:
    """HPO is not a tree. Reading a single parent loses most of the graph."""
    parents = ontology.parents(DANDY_WALKER)
    assert len(parents) >= 3, f"expected a genuinely multi-parent term, got {parents}"
    assert parents == tuple(sorted(parents)), "parents must be returned in a total order"
    assert len(ontology.parents(MICROCEPHALY)) >= 2


def test_many_terms_in_the_slice_are_multi_parent(ontology: HpoOntology) -> None:
    """Multi-parenthood is the norm, not a curiosity of one term."""
    multi = [term for term in ontology.term_ids if len(ontology.parents(term)) > 1]
    assert len(multi) >= 20


def test_ancestor_closure_follows_every_parent(ontology: HpoOntology) -> None:
    """The closure of a multi-parent term contains all of its parents' lineages."""
    closure = frozenset(ontology.ancestor_closure(DANDY_WALKER))
    for parent in ontology.parents(DANDY_WALKER):
        assert frozenset(ontology.ancestor_closure(parent)) <= closure


def test_a_term_reaches_the_root_by_paths_of_different_lengths(ontology: HpoOntology) -> None:
    """Depth is not a function of a term, which is why IC is not derived from it.

    Computed directly here rather than trusted: if the fixture ever loses this
    property, the comment in :mod:`mva.phenotype.corpus` explaining why depth is
    rejected stops being demonstrated by anything.
    """
    shortest = _path_length(ontology, UNEVEN_DEPTH_TERM, minimise=True)
    longest = _path_length(ontology, UNEVEN_DEPTH_TERM, minimise=False)
    assert shortest < longest, "fixture no longer contains an uneven-depth term"
    assert longest - shortest >= 3


def _path_length(ontology: HpoOntology, term: str, *, minimise: bool) -> int:
    """Shortest or longest ``is_a`` path from ``term`` to a root, memoised."""
    memo: dict[str, int] = {}

    def walk(node: str, seen: frozenset[str]) -> int:
        if node in memo:
            return memo[node]
        parents = tuple(p for p in ontology.parents(node) if p not in seen)
        if not parents:
            result = 0
        else:
            lengths = [1 + walk(p, seen | {node}) for p in parents]
            result = min(lengths) if minimise else max(lengths)
        memo[node] = result
        return result

    return walk(term, frozenset())


def test_closures_are_reflexive_and_sorted(ontology: HpoOntology) -> None:
    """A term is its own ancestor and its own descendant.

    Reflexivity is what makes ``sim(t, t)`` fall out of the general case instead of
    needing a special branch that someone will forget.
    """
    ancestors = ontology.ancestor_closure(CATARACT)
    descendants = ontology.descendant_closure(EYE_ABNORMALITY)
    assert CATARACT in ancestors
    assert EYE_ABNORMALITY in descendants
    assert ancestors == tuple(sorted(ancestors))
    assert descendants == tuple(sorted(descendants))


def test_ancestor_and_descendant_are_inverse_relations(ontology: HpoOntology) -> None:
    """``a in ancestors(b)`` if and only if ``b in descendants(a)``.

    Directly guards the confusion that inverts every negative finding.
    """
    ancestors = ontology.ancestor_closure(CATARACT)
    for candidate in ancestors:
        assert CATARACT in ontology.descendant_closure(candidate)
    for descendant in ontology.descendant_closure(LENS_ABNORMALITY):
        assert LENS_ABNORMALITY in ontology.ancestor_closure(descendant)


def test_ancestor_closure_does_not_leak_into_siblings(ontology: HpoOntology) -> None:
    """Going up must not go back down. Cataract's ancestors are not eye disorders at large."""
    ancestors = frozenset(ontology.ancestor_closure(CATARACT))
    siblings = [
        term
        for term in ontology.children(LENS_ABNORMALITY)
        if term != CATARACT and CATARACT not in ontology.descendant_closure(term)
    ]
    assert siblings, "fixture no longer provides a sibling of Cataract"
    for sibling in siblings:
        assert sibling not in ancestors


def test_root_is_the_ancestor_of_everything(ontology: HpoOntology) -> None:
    assert ontology.roots == (ROOT,)
    for term in ontology.term_ids:
        assert ROOT in ontology.ancestor_closure(term)


def test_is_ancestor_of_matches_the_closure(ontology: HpoOntology) -> None:
    assert ontology.is_ancestor_of(EYE_ABNORMALITY, CATARACT)
    assert not ontology.is_ancestor_of(CATARACT, EYE_ABNORMALITY)


# ---------------------------------------------------------------------------
# Hostile graphs
# ---------------------------------------------------------------------------


def test_a_cycle_terminates_instead_of_recursing() -> None:
    """A malformed file must not take the interpreter down with it.

    The release is acyclic, but a hand-edited or truncated one need not be, and a
    recursive closure would exhaust the stack rather than report a bad file.
    """
    ontology = _tiny(
        _term("HP:0000001"),
        _term("HP:0000002", parents=("HP:0000003",)),
        _term("HP:0000003", parents=("HP:0000002",)),
    )
    assert ontology.ancestor_closure("HP:0000002") == ("HP:0000002", "HP:0000003")
    assert ontology.descendant_closure("HP:0000002") == ("HP:0000002", "HP:0000003")


def test_diamond_is_visited_once_not_twice() -> None:
    """Re-convergent paths must not multiply work or duplicate output."""
    ontology = _tiny(
        _term("HP:0000001"),
        _term("HP:0000002", parents=("HP:0000001",)),
        _term("HP:0000003", parents=("HP:0000001",)),
        _term("HP:0000004", parents=("HP:0000002", "HP:0000003")),
    )
    closure = ontology.ancestor_closure("HP:0000004")
    assert closure == ("HP:0000001", "HP:0000002", "HP:0000003", "HP:0000004")
    assert len(closure) == len(set(closure))


def test_duplicate_term_definition_is_rejected() -> None:
    with pytest.raises(IngestionError, match="Duplicate HPO term"):
        _tiny(_term("HP:0000002"), _term("HP:0000002"))


def test_malformed_identifier_in_the_file_names_the_file(tmp_path: Path) -> None:
    path = tmp_path / "bad.obo"
    path.write_text(
        "data-version: hp/releases/2026-06-23\n\n[Term]\nid: HP:123\nname: Broken\n",
        encoding="utf-8",
    )
    with pytest.raises(IngestionError, match=r"bad\.obo"):
        HpoOntology.from_obo(path)


# ---------------------------------------------------------------------------
# Identifier resolution
# ---------------------------------------------------------------------------


def test_alt_id_resolves_to_the_primary_term() -> None:
    """A profile written against an older release must not silently lose terms."""
    ontology = _tiny(
        _term("HP:0000001"),
        _term("HP:0000002", parents=("HP:0000001",), alt_ids=("HP:0009999",)),
    )
    assert ontology.resolve("HP:0009999") == "HP:0000002"
    assert ontology.ancestor_closure("HP:0009999") == ("HP:0000001", "HP:0000002")


def test_obsolete_term_forwards_to_its_replacement() -> None:
    ontology = _tiny(
        _term("HP:0000001"),
        _term("HP:0000002", parents=("HP:0000001",)),
        _term("HP:0000003", obsolete=True, replaced_by="HP:0000002"),
    )
    assert ontology.resolve("HP:0000003") == "HP:0000002"
    assert ontology.obsolete_ids == ("HP:0000003",)


def test_obsolete_term_without_replacement_is_unresolvable_not_a_zero_match() -> None:
    """Unresolvable must be reportable, so the caller can call it an information gap.

    Returning a placeholder node instead would make the term score as "compared and
    found unrelated", which is a claim the release does not support.
    """
    ontology = _tiny(_term("HP:0000001"), _term("HP:0000003", obsolete=True))
    assert ontology.resolve("HP:0000003") is None
    assert ontology.ancestor_closure("HP:0000003") == ()
    assert not ontology.contains("HP:0000003")


def test_resolve_all_separates_known_from_unknown(ontology: HpoOntology) -> None:
    resolved, missing = ontology.resolve_all([CATARACT, "HP:9999999", "hp_0000252"])
    assert resolved == (MICROCEPHALY, CATARACT)
    assert missing == ("HP:9999999",)


def test_identifier_spellings_are_normalised(ontology: HpoOntology) -> None:
    assert ontology.resolve("hp_0000252") == MICROCEPHALY
    assert ontology.resolve("  hp:0000252  ") == MICROCEPHALY


def test_unknown_term_yields_empty_closures_not_a_self_singleton(ontology: HpoOntology) -> None:
    assert ontology.ancestor_closure("HP:9999999") == ()
    assert ontology.descendant_closure("HP:9999999") == ()
    assert ontology.parents("HP:9999999") == ()
    assert ontology.children("HP:9999999") == ()


# ---------------------------------------------------------------------------
# Determinism and shape
# ---------------------------------------------------------------------------


def test_repeat_loads_produce_identical_structure() -> None:
    """No parse-order or set-iteration dependence in the loaded graph (GP-30)."""
    first = HpoOntology.from_obo(SUBSET_OBO)
    second = HpoOntology.from_obo(SUBSET_OBO)
    assert first.term_ids == second.term_ids
    assert first.source_sha256 == second.source_sha256
    for term in first.term_ids:
        assert first.parents(term) == second.parents(term)
        assert first.ancestor_closure(term) == second.ancestor_closure(term)


def test_cached_closure_matches_the_uncached_one(ontology: HpoOntology) -> None:
    """Memoisation must not change an answer, only its cost."""
    fresh = HpoOntology.from_obo(SUBSET_OBO)
    for term in (CATARACT, DANDY_WALKER, MICROCEPHALY, UNEVEN_DEPTH_TERM):
        assert ontology.ancestor_closure(term) == fresh.ancestor_closure(term)
        assert ontology.ancestor_closure(term) == ontology.ancestor_closure(term)


def test_fixture_has_the_shape_the_tests_assume(ontology: HpoOntology) -> None:
    """Guards the fixture itself: a regenerated slice must stay usable."""
    assert 300 <= len(ontology) <= 2000
    for term in (MICROCEPHALY, DANDY_WALKER, CATARACT, EYE_ABNORMALITY, UNEVEN_DEPTH_TERM):
        assert ontology.contains(term), f"{term} missing from the fixture slice"
