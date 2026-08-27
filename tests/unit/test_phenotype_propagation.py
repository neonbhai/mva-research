"""The four-valued logic, propagated over the DAG, tested in both directions.

Every test here has a mirror image. It is not enough to check that observing a
child credits its parent; the same test must check that it does **not** credit the
parent's other children, because an implementation that propagates in both
directions passes the first half and inflates every score in the run.

The four statuses and their propagation:

===============  =====================  =========================================
Status           Propagates             Because
===============  =====================  =========================================
``OBSERVED``     upward (to ancestors)  the specific implies the general
``EXCLUDED``     downward (descendants) the general absent implies the specific
``UNCERTAIN``    neither                an equivocal finding entails nothing
``NOT_ASSESSED`` neither                nobody looked; nothing follows
===============  =====================  =========================================
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mva.models.phenotype import ObservationStatus, Onset, PhenotypeObservation, PhenotypeProfile
from mva.phenotype.ontology import HpoOntology
from mva.phenotype.propagation import InferenceBasis, PropagatedProfile

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
SUBSET_OBO = REPO_ROOT / "tests" / "fixtures" / "hpo" / "hp_subset.obo"

ROOT = "HP:0000001"
ABNORMALITY = "HP:0000118"
EYE = "HP:0000478"
LENS = "HP:0000517"
CATARACT = "HP:0000518"
MICROCEPHALY = "HP:0000252"
SEIZURE = "HP:0001250"
NEPHROBLASTOMA = "HP:0002667"
DANDY_WALKER = "HP:0001305"

#: A real term in the slice that carries a retired ``alt_id``. Both spellings name
#: one term, which is what makes a disagreement between them irreconcilable rather
#: than a matter of precedence.
URETER = "HP:0000069"
URETER_ALT_ID = "HP:0006001"


@pytest.fixture(scope="module")
def ontology() -> HpoOntology:
    return HpoOntology.from_obo(SUBSET_OBO)


def _observation(hpo_id: str, status: ObservationStatus) -> PhenotypeObservation:
    return PhenotypeObservation(
        hpo_id=hpo_id,
        label="fixture term",
        status=status,
        onset=Onset.UNKNOWN,
        provenance="unit_test",
        extraction_confidence=1.0,
        source_excerpt_hash=None,
        notes=None,
    )


def _profile(*pairs: tuple[str, ObservationStatus]) -> PhenotypeProfile:
    return PhenotypeProfile(
        subject_id="PROBAND_TEST",
        observations=tuple(_observation(term, status) for term, status in sorted(pairs)),
        source_artifact="unit_test_profile",
        hpo_version="hp/releases/fixture",
    )


def _propagate(ontology: HpoOntology, *pairs: tuple[str, ObservationStatus]) -> PropagatedProfile:
    return PropagatedProfile(_profile(*pairs), ontology=ontology)


# ---------------------------------------------------------------------------
# OBSERVED propagates UP, and only up
# ---------------------------------------------------------------------------


def test_observed_child_makes_every_ancestor_observed(ontology: HpoOntology) -> None:
    """The TD-04 case: an observed child must credit its associated parent."""
    propagated = _propagate(ontology, (CATARACT, ObservationStatus.OBSERVED))
    for ancestor in ontology.ancestor_closure(CATARACT):
        inferred = propagated.status_of(ancestor)
        assert inferred.status is ObservationStatus.OBSERVED, ancestor
    assert propagated.status_of(LENS).basis is InferenceBasis.ANCESTOR_OF_OBSERVED
    assert propagated.status_of(LENS).source_term == CATARACT
    assert propagated.status_of(CATARACT).basis is InferenceBasis.RECORDED


def test_observed_parent_does_not_make_its_children_observed(ontology: HpoOntology) -> None:
    """The mirror image. Observing the general does NOT imply the specific.

    An implementation that propagates OBSERVED downward credits every gene under a
    broad finding, which inflates the whole ranking.
    """
    propagated = _propagate(ontology, (EYE, ObservationStatus.OBSERVED))
    assert propagated.status_of(EYE).status is ObservationStatus.OBSERVED
    for descendant in ontology.descendant_closure(EYE):
        if descendant == EYE:
            continue
        inferred = propagated.status_of(descendant)
        assert inferred.status is ObservationStatus.NOT_ASSESSED, descendant
        assert inferred.basis is InferenceBasis.UNRECORDED


def test_observed_term_says_nothing_about_its_siblings(ontology: HpoOntology) -> None:
    """Up, then not back down: a sibling shares an ancestor, not a status."""
    propagated = _propagate(ontology, (CATARACT, ObservationStatus.OBSERVED))
    siblings = [
        term
        for term in ontology.children(LENS)
        if term != CATARACT and CATARACT not in ontology.descendant_closure(term)
    ]
    assert siblings, "fixture no longer provides a sibling of Cataract"
    for sibling in siblings:
        assert propagated.status_of(sibling).status is ObservationStatus.NOT_ASSESSED


def test_observed_closure_lists_ancestors_and_nothing_else(ontology: HpoOntology) -> None:
    propagated = _propagate(ontology, (MICROCEPHALY, ObservationStatus.OBSERVED))
    assert propagated.observed_closure == ontology.ancestor_closure(MICROCEPHALY)
    assert propagated.recorded_observed == (MICROCEPHALY,)
    assert propagated.summary.implied_observed == len(propagated.observed_closure) - 1


# ---------------------------------------------------------------------------
# EXCLUDED propagates DOWN, and only down
# ---------------------------------------------------------------------------


def test_excluded_parent_excludes_every_descendant(ontology: HpoOntology) -> None:
    """ "No abnormality of the eye" entails "no cataract"."""
    propagated = _propagate(ontology, (EYE, ObservationStatus.EXCLUDED))
    for descendant in ontology.descendant_closure(EYE):
        assert propagated.status_of(descendant).status is ObservationStatus.EXCLUDED, descendant
    assert propagated.status_of(CATARACT).basis is InferenceBasis.DESCENDANT_OF_EXCLUDED
    assert propagated.status_of(CATARACT).source_term == EYE


def test_excluded_child_does_not_exclude_its_ancestors(ontology: HpoOntology) -> None:
    """The mirror image, and the more dangerous direction.

    "No cataract" plainly does not mean "no abnormality of the eye". Propagating a
    negative upward would manufacture a broad negative finding from a narrow one
    and penalise every gene under that branch.
    """
    propagated = _propagate(ontology, (CATARACT, ObservationStatus.EXCLUDED))
    assert propagated.status_of(CATARACT).status is ObservationStatus.EXCLUDED
    for ancestor in ontology.ancestor_closure(CATARACT):
        if ancestor == CATARACT:
            continue
        inferred = propagated.status_of(ancestor)
        assert inferred.status is ObservationStatus.NOT_ASSESSED, ancestor
        assert inferred.basis is InferenceBasis.UNRECORDED


def test_excluded_closure_lists_descendants_and_nothing_else(ontology: HpoOntology) -> None:
    propagated = _propagate(ontology, (LENS, ObservationStatus.EXCLUDED))
    assert propagated.excluded_closure == ontology.descendant_closure(LENS)
    assert propagated.recorded_excluded == (LENS,)


def test_the_two_directions_are_not_interchangeable(ontology: HpoOntology) -> None:
    """Swapping the two closures would be caught here even if every other test passed.

    Observing X and excluding X are recorded identically apart from the status, so a
    scorer that used one closure for both would produce the same term set for each.
    """
    observed = _propagate(ontology, (CATARACT, ObservationStatus.OBSERVED))
    excluded = _propagate(ontology, (CATARACT, ObservationStatus.EXCLUDED))
    assert observed.observed_closure == ontology.ancestor_closure(CATARACT)
    assert excluded.excluded_closure == ontology.descendant_closure(CATARACT)
    assert observed.observed_closure != excluded.excluded_closure


# ---------------------------------------------------------------------------
# UNCERTAIN and NOT_ASSESSED propagate in neither direction
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("status", [ObservationStatus.UNCERTAIN, ObservationStatus.NOT_ASSESSED])
def test_uninformative_status_does_not_propagate_anywhere(
    ontology: HpoOntology, status: ObservationStatus
) -> None:
    """An equivocal or absent assessment entails nothing, in either direction."""
    propagated = _propagate(ontology, (LENS, status))
    assert propagated.status_of(LENS).status is status
    assert propagated.observed_closure == ()
    assert propagated.excluded_closure == ()
    for other in (*ontology.ancestor_closure(LENS), *ontology.descendant_closure(LENS)):
        if other == LENS:
            continue
        assert propagated.status_of(other).status is ObservationStatus.NOT_ASSESSED


def test_uncertain_is_not_observed_and_not_excluded(ontology: HpoOntology) -> None:
    """The distinction the whole package exists to preserve, at the closure layer."""
    propagated = _propagate(ontology, (SEIZURE, ObservationStatus.UNCERTAIN))
    inferred = propagated.status_of(SEIZURE)
    assert inferred.status is ObservationStatus.UNCERTAIN
    assert inferred.status is not ObservationStatus.OBSERVED
    assert inferred.status is not ObservationStatus.EXCLUDED
    assert SEIZURE not in propagated.observed_closure
    assert SEIZURE not in propagated.excluded_closure


def test_unrecorded_term_is_not_assessed_not_excluded(ontology: HpoOntology) -> None:
    """A term nobody wrote down is a term nobody assessed (GP-14)."""
    propagated = _propagate(ontology, (MICROCEPHALY, ObservationStatus.OBSERVED))
    inferred = propagated.status_of(NEPHROBLASTOMA)
    assert inferred.status is ObservationStatus.NOT_ASSESSED
    assert inferred.basis is InferenceBasis.UNRECORDED


def test_all_four_statuses_survive_one_profile(ontology: HpoOntology) -> None:
    """End-to-end: four statuses in, four distinct statuses out."""
    propagated = _propagate(
        ontology,
        (MICROCEPHALY, ObservationStatus.OBSERVED),
        (CATARACT, ObservationStatus.EXCLUDED),
        (SEIZURE, ObservationStatus.UNCERTAIN),
        (NEPHROBLASTOMA, ObservationStatus.NOT_ASSESSED),
    )
    assert propagated.status_of(MICROCEPHALY).status is ObservationStatus.OBSERVED
    assert propagated.status_of(CATARACT).status is ObservationStatus.EXCLUDED
    assert propagated.status_of(SEIZURE).status is ObservationStatus.UNCERTAIN
    assert propagated.status_of(NEPHROBLASTOMA).status is ObservationStatus.NOT_ASSESSED
    assert propagated.summary.recorded_observed == 1
    assert propagated.summary.recorded_excluded == 1
    assert propagated.summary.recorded_uncertain == 1
    assert propagated.summary.recorded_not_assessed == 1


# ---------------------------------------------------------------------------
# Precedence and conflict
# ---------------------------------------------------------------------------


def test_entailed_observed_overrides_a_recorded_not_assessed(ontology: HpoOntology) -> None:
    """If the specific finding was seen, the general one is not 'unassessed'.

    The form said nobody looked at the parent term; the record of the child says
    otherwise. Entailment wins over an absence of information, never the reverse.
    """
    propagated = _propagate(
        ontology,
        (CATARACT, ObservationStatus.OBSERVED),
        (LENS, ObservationStatus.NOT_ASSESSED),
    )
    inferred = propagated.status_of(LENS)
    assert inferred.status is ObservationStatus.OBSERVED
    assert inferred.basis is InferenceBasis.ANCESTOR_OF_OBSERVED


def test_a_recorded_status_outranks_an_entailed_one(ontology: HpoOntology) -> None:
    """A clinician's direct assertion about a term beats an inference about it."""
    propagated = _propagate(
        ontology,
        (EYE, ObservationStatus.EXCLUDED),
        (CATARACT, ObservationStatus.OBSERVED),
    )
    assert propagated.status_of(CATARACT).status is ObservationStatus.OBSERVED
    assert propagated.status_of(CATARACT).basis is InferenceBasis.RECORDED
    assert propagated.status_of(EYE).status is ObservationStatus.EXCLUDED
    assert propagated.status_of(EYE).basis is InferenceBasis.RECORDED


def test_contradictory_record_is_reported_not_silently_resolved(
    ontology: HpoOntology,
) -> None:
    """An inconsistent record is a data-quality problem for a human, not a tie-break."""
    propagated = _propagate(
        ontology,
        (EYE, ObservationStatus.EXCLUDED),
        (CATARACT, ObservationStatus.OBSERVED),
    )
    assert CATARACT in propagated.conflicts
    assert propagated.summary.conflicts == len(propagated.conflicts)
    assert propagated.conflicts == tuple(sorted(propagated.conflicts))


def test_a_consistent_record_reports_no_conflicts(ontology: HpoOntology) -> None:
    propagated = _propagate(
        ontology,
        (MICROCEPHALY, ObservationStatus.OBSERVED),
        (CATARACT, ObservationStatus.EXCLUDED),
    )
    assert propagated.conflicts == ()


# ---------------------------------------------------------------------------
# Terms the release does not contain
# ---------------------------------------------------------------------------


def test_unresolvable_term_is_an_information_gap_not_a_mismatch(
    ontology: HpoOntology,
) -> None:
    """A term outside the release contributes nothing and is named, not dropped."""
    propagated = _propagate(
        ontology,
        (MICROCEPHALY, ObservationStatus.OBSERVED),
        ("HP:9999999", ObservationStatus.OBSERVED),
    )
    assert propagated.unresolved_terms == ("HP:9999999",)
    assert "HP:9999999" not in propagated.observed_closure
    inferred = propagated.status_of("HP:9999999")
    assert inferred.status is ObservationStatus.NOT_ASSESSED
    assert inferred.basis is InferenceBasis.UNRESOLVED


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------


def test_closures_are_sorted_and_repeatable(ontology: HpoOntology) -> None:
    """No set iteration reaches an output position (GP-30)."""
    pairs = (
        (DANDY_WALKER, ObservationStatus.OBSERVED),
        (MICROCEPHALY, ObservationStatus.OBSERVED),
        (EYE, ObservationStatus.EXCLUDED),
    )
    first = _propagate(ontology, *pairs)
    second = _propagate(ontology, *pairs)
    assert first.observed_closure == second.observed_closure
    assert first.excluded_closure == second.excluded_closure
    assert first.observed_closure == tuple(sorted(first.observed_closure))
    assert first.excluded_closure == tuple(sorted(first.excluded_closure))
    assert first.recorded_observed == tuple(sorted(first.recorded_observed))


def test_input_order_does_not_change_the_result(ontology: HpoOntology) -> None:
    forward = _propagate(
        ontology,
        (MICROCEPHALY, ObservationStatus.OBSERVED),
        (CATARACT, ObservationStatus.EXCLUDED),
    )
    reverse = _propagate(
        ontology,
        (CATARACT, ObservationStatus.EXCLUDED),
        (MICROCEPHALY, ObservationStatus.OBSERVED),
    )
    assert forward.observed_closure == reverse.observed_closure
    assert forward.excluded_closure == reverse.excluded_closure
    assert forward.conflicts == reverse.conflicts


# ---------------------------------------------------------------------------
# One term recorded twice under two spellings
# ---------------------------------------------------------------------------


def test_the_fixture_really_does_alias_the_two_identifiers(ontology: HpoOntology) -> None:
    """Guard for the tests below: they are meaningless unless both ids are one term."""
    assert ontology.resolve(URETER_ALT_ID) == URETER
    term = ontology.term(URETER)
    assert term is not None
    assert URETER_ALT_ID in term.alt_ids


def test_primary_and_alias_with_opposite_statuses_are_conflicting(
    ontology: HpoOntology,
) -> None:
    """EXCLUDED under the retired id and OBSERVED under the current one is a conflict.

    This module used to rank the two statuses and keep OBSERVED. That did two
    things at once, both wrong: it converted a genuine contradiction into positive
    evidence, and it removed the term from ``conflicts`` so nothing downstream could
    see that anything had happened. ``alt_id`` exists because real annotations
    accumulate under retired identifiers, so this input is realistic, not exotic.
    """
    propagated = _propagate(
        ontology,
        (URETER, ObservationStatus.EXCLUDED),
        (URETER_ALT_ID, ObservationStatus.OBSERVED),
    )

    assert URETER in propagated.conflicts
    assert URETER in propagated.direct_conflicts
    assert propagated.summary.direct_conflicts == 1

    # Neither direction silently wins.
    inferred = propagated.status_of(URETER)
    assert inferred.status is not ObservationStatus.OBSERVED
    assert inferred.status is not ObservationStatus.EXCLUDED
    assert inferred.status is ObservationStatus.UNCERTAIN

    # And the term is withdrawn from scoring in both directions: it seeds neither
    # closure, so it can neither support nor argue against any gene.
    assert propagated.observed_closure == ()
    assert propagated.excluded_closure == ()
    assert propagated.recorded_observed == ()
    assert propagated.recorded_excluded == ()


def test_the_alias_spelling_reports_the_same_conflicted_status(
    ontology: HpoOntology,
) -> None:
    """Asking under either spelling must give the same answer; they are one term."""
    propagated = _propagate(
        ontology,
        (URETER, ObservationStatus.EXCLUDED),
        (URETER_ALT_ID, ObservationStatus.OBSERVED),
    )
    assert propagated.status_of(URETER_ALT_ID).status is ObservationStatus.UNCERTAIN
    assert propagated.status_of(URETER).status is ObservationStatus.UNCERTAIN


def test_the_conflict_does_not_depend_on_which_spelling_came_first(
    ontology: HpoOntology,
) -> None:
    """A precedence rule would be order-sensitive in spirit; a conflict is not."""
    forward = _propagate(
        ontology,
        (URETER, ObservationStatus.EXCLUDED),
        (URETER_ALT_ID, ObservationStatus.OBSERVED),
    )
    reverse = _propagate(
        ontology,
        (URETER_ALT_ID, ObservationStatus.OBSERVED),
        (URETER, ObservationStatus.EXCLUDED),
    )
    assert forward.conflicts == reverse.conflicts
    assert forward.status_of(URETER).status is reverse.status_of(URETER).status


@pytest.mark.parametrize(
    ("primary_status", "alias_status", "expected"),
    [
        # Unanimous: nothing is being decided, so collapse.
        (ObservationStatus.OBSERVED, ObservationStatus.OBSERVED, ObservationStatus.OBSERVED),
        (ObservationStatus.EXCLUDED, ObservationStatus.EXCLUDED, ObservationStatus.EXCLUDED),
        # One informative status plus an uninformative duplicate is a redundant row,
        # not a disagreement: discarding the finding would lose real evidence.
        (ObservationStatus.OBSERVED, ObservationStatus.NOT_ASSESSED, ObservationStatus.OBSERVED),
        (ObservationStatus.EXCLUDED, ObservationStatus.NOT_ASSESSED, ObservationStatus.EXCLUDED),
        (ObservationStatus.OBSERVED, ObservationStatus.UNCERTAIN, ObservationStatus.OBSERVED),
        (ObservationStatus.EXCLUDED, ObservationStatus.UNCERTAIN, ObservationStatus.EXCLUDED),
        # Both uninformative: neither carries weight, so no evidence can be created
        # or destroyed; keep the more specific description of the gap.
        (
            ObservationStatus.UNCERTAIN,
            ObservationStatus.NOT_ASSESSED,
            ObservationStatus.UNCERTAIN,
        ),
    ],
)
def test_only_identical_or_uninformative_duplicates_collapse(
    ontology: HpoOntology,
    primary_status: ObservationStatus,
    alias_status: ObservationStatus,
    expected: ObservationStatus,
) -> None:
    """Every duplicate pairing that is *not* a contradiction, and what it collapses to."""
    propagated = _propagate(ontology, (URETER, primary_status), (URETER_ALT_ID, alias_status))
    assert propagated.status_of(URETER).status is expected
    assert propagated.direct_conflicts == ()


def test_a_single_recorded_term_is_never_reported_as_conflicted(
    ontology: HpoOntology,
) -> None:
    """The common case must not acquire a conflict from the reconciliation path."""
    for status in ObservationStatus:
        propagated = _propagate(ontology, (URETER, status))
        assert propagated.direct_conflicts == ()
        assert propagated.status_of(URETER).status is status
