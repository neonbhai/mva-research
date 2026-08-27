"""Ontology-aware gene-phenotype scoring: the properties TD-04 was opened for.

The debt this pays down is one sentence: *an observed child term does not credit
its associated parent term, so true matches are under-scored*. That is the first
test below. Everything after it exists because the obvious fix breaks something
else:

* crediting a parent must not also credit siblings or children (direction);
* crediting through the graph must not resurrect ``NOT_ASSESSED`` as evidence
  (GP-14);
* ``EXCLUDED`` must still be able to argue against a gene, now through the
  *downward* closure;
* a gene the knowledge base says nothing about must stay distinguishable from a
  gene it says plenty about, none of which fits (GP-14 again);
* and the whole thing must produce identical bytes twice (GP-30).
"""

from __future__ import annotations

import os
import subprocess
import sys
import textwrap
from collections.abc import Sequence
from pathlib import Path

import pytest

from mva.clock import Clock, demo_clock
from mva.models.evidence import EvidenceDirection
from mva.models.phenotype import (
    ObservationStatus,
    Onset,
    PhenotypeObservation,
    PhenotypeProfile,
)
from mva.phenotype.corpus import CorpusKind
from mva.phenotype.hpo import GeneAssociation, GenePhenotypeIndex
from mva.phenotype.scoring import (
    NEUTRAL_SCORE,
    GeneAnnotationStatus,
    PhenotypeMatch,
    ScoringMode,
    score_all_genes,
    score_gene_phenotype,
)
from mva.phenotype.semantics import HpoResourceSet, PhenotypeSemantics
from mva.phenotype.similarity import SimilarityMeasure

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
HPO_FIXTURES = REPO_ROOT / "tests" / "fixtures" / "hpo"
SUBSET_OBO = HPO_FIXTURES / "hp_subset.obo"
SUBSET_HPOA = HPO_FIXTURES / "phenotype_subset.hpoa"
SUBSET_G2P = HPO_FIXTURES / "genes_to_phenotype_subset.txt"

EYE = "HP:0000478"
LENS = "HP:0000517"
CATARACT = "HP:0000518"
MICROCEPHALY = "HP:0000252"
SEIZURE = "HP:0001250"
NEPHROBLASTOMA = "HP:0002667"
DANDY_WALKER = "HP:0001305"
IUGR = "HP:0001511"
#: A real term in the slice carrying a retired ``alt_id``; both name one term.
URETER = "HP:0000069"
URETER_ALT_ID = "HP:0006001"

KNOWLEDGE_VERSION = "unit-test-v0"


@pytest.fixture(scope="module")
def semantics() -> PhenotypeSemantics:
    return HpoResourceSet(ontology_path=SUBSET_OBO, annotation_path=SUBSET_HPOA).load()


@pytest.fixture
def clock() -> Clock:
    return demo_clock()


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


def _index(*rows: tuple[str, str, str]) -> GenePhenotypeIndex:
    """Build an index from ``(gene, hpo_id, strength)`` triples."""
    return GenePhenotypeIndex(
        [
            GeneAssociation(
                gene_symbol=gene,
                hpo_id=hpo_id,
                label="fixture term",
                association_strength=strength,
                source="unit_test_curation",
                version="v0",
            )
            for gene, hpo_id, strength in rows
        ],
        version=KNOWLEDGE_VERSION,
    )


def _score(
    gene: str,
    *,
    profile: PhenotypeProfile,
    index: GenePhenotypeIndex,
    clock: Clock,
    semantics: PhenotypeSemantics | None,
) -> PhenotypeMatch:
    return score_gene_phenotype(
        gene, profile=profile, index=index, clock=clock, semantics=semantics
    )


# ---------------------------------------------------------------------------
# TD-04: the parent term is credited
# ---------------------------------------------------------------------------


def test_observed_child_credits_a_gene_annotated_with_the_parent(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """The debt itself. Under exact matching this gene scored as a non-match."""
    profile = _profile((CATARACT, ObservationStatus.OBSERVED))
    index = _index(("PARENTGENE", LENS, "strong"))

    exact = _score("PARENTGENE", profile=profile, index=index, clock=clock, semantics=None)
    ontology_aware = _score(
        "PARENTGENE", profile=profile, index=index, clock=clock, semantics=semantics
    )

    assert exact.matched_terms == (), "exact matching cannot see the parent (that is TD-04)"
    assert ontology_aware.matched_terms == (LENS,)
    assert ontology_aware.implied_matched_terms == (LENS,)
    assert ontology_aware.score > exact.score
    assert ontology_aware.mode is ScoringMode.ONTOLOGY
    assert exact.mode is ScoringMode.EXACT


def test_an_entailed_match_is_labelled_as_entailed_not_as_observed_data(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """An inferred match is weaker evidence than a recorded one and must say so."""
    profile = _profile((CATARACT, ObservationStatus.OBSERVED))
    index = _index(("PARENTGENE", LENS, "strong"), ("EXACTGENE", CATARACT, "strong"))

    entailed = _score("PARENTGENE", profile=profile, index=index, clock=clock, semantics=semantics)
    direct = _score("EXACTGENE", profile=profile, index=index, clock=clock, semantics=semantics)

    assert entailed.implied_matched_terms == (LENS,)
    assert direct.implied_matched_terms == ()
    entailed_claims = " ".join(item.claim for item in entailed.evidence)
    assert "entailed OBSERVED" in entailed_claims
    assert "entailed OBSERVED" not in " ".join(item.claim for item in direct.evidence)


def test_the_exact_match_still_outscores_a_more_general_one(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """Crediting an ancestor must not make it as good as the recorded finding.

    Information content is what keeps them apart: the ancestor is a more general
    statement, more diseases in the corpus reach it, and it therefore earns less
    credit. Nothing about graph distance is involved.
    """
    profile = _profile((CATARACT, ObservationStatus.OBSERVED))
    index = _index(("BROADGENE", EYE, "strong"), ("EXACTGENE", CATARACT, "strong"))
    entailed = _score("BROADGENE", profile=profile, index=index, clock=clock, semantics=semantics)
    direct = _score("EXACTGENE", profile=profile, index=index, clock=clock, semantics=semantics)
    assert entailed.implied_matched_terms == (EYE,)
    assert direct.score > entailed.score


def test_terms_the_corpus_cannot_distinguish_score_the_same(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """A deliberate, documented consequence of corpus-derived information content.

    In this slice every disease annotated with *Abnormality of the lens* is
    annotated with *Cataract*, so the two terms carry identical information and Lin
    rates them identical. That is the measure being honest about the corpus, not a
    bug — and it is exactly the behaviour a depth-based proxy would get wrong, since
    the two sit at different depths.
    """
    ic = semantics.information_content
    assert ic.ic(LENS) == ic.ic(CATARACT), "fixture no longer has an indistinguishable pair"
    profile = _profile((CATARACT, ObservationStatus.OBSERVED))
    index = _index(("PARENTGENE", LENS, "strong"), ("EXACTGENE", CATARACT, "strong"))
    entailed = _score("PARENTGENE", profile=profile, index=index, clock=clock, semantics=semantics)
    direct = _score("EXACTGENE", profile=profile, index=index, clock=clock, semantics=semantics)
    assert entailed.score == direct.score


def test_a_distant_gene_scores_below_a_parent_gene(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    profile = _profile((CATARACT, ObservationStatus.OBSERVED))
    index = _index(("PARENTGENE", LENS, "strong"), ("FARGENE", SEIZURE, "strong"))
    near = _score("PARENTGENE", profile=profile, index=index, clock=clock, semantics=semantics)
    far = _score("FARGENE", profile=profile, index=index, clock=clock, semantics=semantics)
    assert near.score > far.score


# ---------------------------------------------------------------------------
# Direction: up for OBSERVED, and only up
# ---------------------------------------------------------------------------


def test_observing_a_parent_does_not_credit_a_gene_annotated_with_the_child(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """The inflation bug. Observing the general does not establish the specific.

    A gene curated for *Cataract* must not be recorded as matched because the
    record says only *Abnormality of the lens*.
    """
    profile = _profile((LENS, ObservationStatus.OBSERVED))
    index = _index(("CHILDGENE", CATARACT, "strong"))
    match = _score("CHILDGENE", profile=profile, index=index, clock=clock, semantics=semantics)
    assert match.matched_terms == ()
    assert match.unassessed_terms == (CATARACT,)
    assert match.implied_matched_terms == ()


def test_observing_one_child_does_not_credit_a_sibling(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """Up, never back down: siblings share an ancestor, not a status."""
    profile = _profile((CATARACT, ObservationStatus.OBSERVED))
    ontology = semantics.ontology
    sibling = next(
        term
        for term in ontology.children(LENS)
        if term != CATARACT and CATARACT not in ontology.descendant_closure(term)
    )
    index = _index(("SIBLINGGENE", sibling, "strong"))
    match = _score("SIBLINGGENE", profile=profile, index=index, clock=clock, semantics=semantics)
    assert match.matched_terms == ()
    assert match.unassessed_terms == (sibling,)


# ---------------------------------------------------------------------------
# Direction: down for EXCLUDED, and only down
# ---------------------------------------------------------------------------


def test_excluding_a_parent_contradicts_a_gene_annotated_with_the_child(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """EXCLUDED is real negative evidence and it propagates downward."""
    profile = _profile(
        (MICROCEPHALY, ObservationStatus.OBSERVED),
        (EYE, ObservationStatus.EXCLUDED),
    )
    index = _index(("EYEGENE", CATARACT, "strong"))
    match = _score("EYEGENE", profile=profile, index=index, clock=clock, semantics=semantics)
    assert match.contradicted_terms == (CATARACT,)
    assert match.implied_contradicted_terms == (CATARACT,)
    assert match.has_contradiction
    assert match.score < NEUTRAL_SCORE


def test_excluding_a_child_does_not_contradict_a_gene_annotated_with_the_parent(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """The mirror image, and the dangerous one.

    "No cataract" does not mean "no abnormality of the lens". Propagating the
    negative upward would penalise every gene under that branch on the strength of
    one narrow finding.
    """
    profile = _profile((CATARACT, ObservationStatus.EXCLUDED))
    index = _index(("LENSGENE", LENS, "strong"))
    match = _score("LENSGENE", profile=profile, index=index, clock=clock, semantics=semantics)
    assert match.contradicted_terms == ()
    assert match.unassessed_terms == (LENS,)


def test_a_contradicted_gene_ranks_below_a_gene_nobody_has_data_for(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """Argued against must rank below not argued about. GP-14 as an ordering."""
    profile = _profile(
        (MICROCEPHALY, ObservationStatus.OBSERVED),
        (EYE, ObservationStatus.EXCLUDED),
    )
    index = _index(("EYEGENE", CATARACT, "strong"))
    contradicted = _score("EYEGENE", profile=profile, index=index, clock=clock, semantics=semantics)
    unknown = _score("NOSUCHGENE", profile=profile, index=index, clock=clock, semantics=semantics)
    assert unknown.score == NEUTRAL_SCORE
    assert contradicted.score < unknown.score


# ---------------------------------------------------------------------------
# All four statuses survive
# ---------------------------------------------------------------------------


def test_the_four_statuses_land_in_four_different_buckets(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    profile = _profile(
        (MICROCEPHALY, ObservationStatus.OBSERVED),
        (CATARACT, ObservationStatus.EXCLUDED),
        (SEIZURE, ObservationStatus.UNCERTAIN),
        (NEPHROBLASTOMA, ObservationStatus.NOT_ASSESSED),
    )
    index = _index(
        ("FOURGENE", MICROCEPHALY, "strong"),
        ("FOURGENE", CATARACT, "moderate"),
        ("FOURGENE", SEIZURE, "moderate"),
        ("FOURGENE", NEPHROBLASTOMA, "supporting"),
    )
    match = _score("FOURGENE", profile=profile, index=index, clock=clock, semantics=semantics)
    assert match.matched_terms == (MICROCEPHALY,)
    assert match.contradicted_terms == (CATARACT,)
    assert match.unassessed_terms == (NEPHROBLASTOMA,)
    # UNCERTAIN is in none of the three scoring buckets and is reported separately.
    assert SEIZURE not in match.matched_terms
    assert SEIZURE not in match.contradicted_terms
    assert SEIZURE not in match.unassessed_terms
    assert "UNCERTAIN" in " ".join(item.claim for item in match.evidence)


def test_uncertain_is_not_observed(semantics: PhenotypeSemantics, clock: Clock) -> None:
    """An equivocal report must not be worth as much as a recorded finding."""
    profile_uncertain = _profile((CATARACT, ObservationStatus.UNCERTAIN))
    profile_observed = _profile((CATARACT, ObservationStatus.OBSERVED))
    index = _index(("GENE", CATARACT, "strong"))
    uncertain = _score(
        "GENE", profile=profile_uncertain, index=index, clock=clock, semantics=semantics
    )
    observed = _score(
        "GENE", profile=profile_observed, index=index, clock=clock, semantics=semantics
    )
    assert uncertain.matched_terms == ()
    assert observed.matched_terms == (CATARACT,)
    assert uncertain.score < observed.score


def test_uncertain_is_not_excluded(semantics: PhenotypeSemantics, clock: Clock) -> None:
    profile = _profile((CATARACT, ObservationStatus.UNCERTAIN))
    index = _index(("GENE", CATARACT, "strong"))
    match = _score("GENE", profile=profile, index=index, clock=clock, semantics=semantics)
    assert match.contradicted_terms == ()
    assert not match.has_contradiction


def test_uncertain_does_not_propagate_to_ancestors(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    profile = _profile((CATARACT, ObservationStatus.UNCERTAIN))
    index = _index(("PARENTGENE", LENS, "strong"))
    match = _score("PARENTGENE", profile=profile, index=index, clock=clock, semantics=semantics)
    assert match.matched_terms == ()
    assert match.unassessed_terms == (LENS,)


# ---------------------------------------------------------------------------
# GP-14: NOT_ASSESSED carries no weight, in either direction
# ---------------------------------------------------------------------------


def test_not_assessed_association_does_not_change_the_specificity_component(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """The regression guard, stated exactly.

    A gene associated with a term nobody assessed must have the same specificity —
    the same numerator and the same denominator — as the same gene without that
    association. This is the numeric form of "absence of information is not
    negative information".
    """
    profile = _profile(
        (MICROCEPHALY, ObservationStatus.OBSERVED),
        (NEPHROBLASTOMA, ObservationStatus.NOT_ASSESSED),
    )
    lean = _index(("GENE", MICROCEPHALY, "strong"))
    padded = _index(("GENE", MICROCEPHALY, "strong"), ("GENE", NEPHROBLASTOMA, "definitive"))

    thin = _score("GENE", profile=profile, index=lean, clock=clock, semantics=semantics)
    thick = _score("GENE", profile=profile, index=padded, clock=clock, semantics=semantics)

    assert thick.unassessed_terms == (NEPHROBLASTOMA,)
    assert thin.breakdown.specificity == thick.breakdown.specificity
    assert thin.breakdown.matched_weight == thick.breakdown.matched_weight
    assert thin.breakdown.informative_weight == thick.breakdown.informative_weight


def test_not_assessed_association_never_lowers_the_score(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """The directional guarantee GP-14 actually asserts: no penalty, ever.

    Coverage is a maximum over the gene's curated terms, so an extra association can
    only raise a best match. It is therefore impossible for a gene to be punished
    for being curated with a feature nobody checked — which is the failure mode this
    package exists to prevent.
    """
    profile = _profile(
        (MICROCEPHALY, ObservationStatus.OBSERVED),
        (IUGR, ObservationStatus.OBSERVED),
    )
    lean = _index(("GENE", MICROCEPHALY, "strong"))
    for extra in (NEPHROBLASTOMA, SEIZURE, CATARACT, DANDY_WALKER):
        padded = _index(("GENE", MICROCEPHALY, "strong"), ("GENE", extra, "definitive"))
        thin = _score("GENE", profile=profile, index=lean, clock=clock, semantics=semantics)
        thick = _score("GENE", profile=profile, index=padded, clock=clock, semantics=semantics)
        assert thick.score >= thin.score, f"adding an unassessed {extra} lowered the score"


def test_an_unrecorded_term_is_treated_as_unassessed_not_absent(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """A term the workup never mentioned must not become negative evidence."""
    profile = _profile((MICROCEPHALY, ObservationStatus.OBSERVED))
    index = _index(("GENE", MICROCEPHALY, "strong"), ("GENE", NEPHROBLASTOMA, "strong"))
    match = _score("GENE", profile=profile, index=index, clock=clock, semantics=semantics)
    assert match.unassessed_terms == (NEPHROBLASTOMA,)
    assert match.contradicted_terms == ()


# ---------------------------------------------------------------------------
# GP-14: no knowledge is not the same as knowledge that does not fit
# ---------------------------------------------------------------------------


def test_no_annotations_is_distinguishable_from_annotations_that_do_not_fit(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """Requirement in one assertion: the two states differ in both status and score.

    "We hold nothing about this gene" is neutral. "We hold five features and none of
    them is anywhere near this patient" is a real, negative finding about a real
    annotation set, and burying it at the same neutral value throws away the only
    thing the knowledge base actually told us.
    """
    profile = _profile(
        (MICROCEPHALY, ObservationStatus.OBSERVED),
        (IUGR, ObservationStatus.OBSERVED),
    )
    index = _index(("MISFITGENE", CATARACT, "strong"), ("MISFITGENE", SEIZURE, "moderate"))

    silent = _score("UNKNOWNGENE", profile=profile, index=index, clock=clock, semantics=semantics)
    misfit = _score("MISFITGENE", profile=profile, index=index, clock=clock, semantics=semantics)

    assert silent.annotation_status is GeneAnnotationStatus.NO_ANNOTATIONS
    assert misfit.annotation_status is GeneAnnotationStatus.ANNOTATED
    assert silent.score == NEUTRAL_SCORE
    assert misfit.score < silent.score
    assert silent.breakdown.coverage is None, "no annotations means no computable coverage"
    assert misfit.breakdown.coverage is not None


def test_no_annotations_emits_its_own_evidence_item(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """The neutral score alone cannot tell a reader which kind of ignorance it is."""
    profile = _profile((MICROCEPHALY, ObservationStatus.OBSERVED))
    index = _index(("OTHER", CATARACT, "strong"))
    match = _score("UNKNOWNGENE", profile=profile, index=index, clock=clock, semantics=semantics)
    notices = [item for item in match.evidence if "no curated HPO terms" in item.claim]
    assert len(notices) == 1
    assert notices[0].direction is EvidenceDirection.NEUTRAL
    assert notices[0].limitations


def test_similarity_is_never_imputed_as_zero_for_missing_annotation_data(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """A term the corpus cannot place is excluded from the mean, not scored as 0.

    The gene corpus in the fixture has no annotation reaching
    ``HP:0200024``; scoring against it must report an uncomparable term rather than
    quietly averaging in a zero and depressing a real match.
    """
    gene_semantics = HpoResourceSet(
        ontology_path=SUBSET_OBO,
        annotation_path=SUBSET_G2P,
        annotation_kind=CorpusKind.GENE,
    ).load()
    unplaceable = "HP:0200024"
    assert gene_semantics.information_content.ic(unplaceable) is None

    profile = _profile(
        (MICROCEPHALY, ObservationStatus.OBSERVED),
        (unplaceable, ObservationStatus.OBSERVED),
    )
    index = _index(("GENE", MICROCEPHALY, "strong"))

    with_gap = _score("GENE", profile=profile, index=index, clock=clock, semantics=gene_semantics)
    without_gap = _score(
        "GENE",
        profile=_profile((MICROCEPHALY, ObservationStatus.OBSERVED)),
        index=index,
        clock=clock,
        semantics=gene_semantics,
    )
    assert with_gap.breakdown.uncomparable_observed_count == 1
    assert with_gap.breakdown.coverage == without_gap.breakdown.coverage
    assert "could not be placed" in with_gap.rationale


def test_a_gene_term_without_information_content_keeps_its_curated_weight(
    clock: Clock,
) -> None:
    """Scaling an unmeasurable term to zero would erase a real curated association."""
    gene_semantics = HpoResourceSet(
        ontology_path=SUBSET_OBO,
        annotation_path=SUBSET_G2P,
        annotation_kind=CorpusKind.GENE,
    ).load()
    unplaceable = "HP:0200024"
    profile = _profile((unplaceable, ObservationStatus.OBSERVED))
    index = _index(("GENE", unplaceable, "strong"))
    match = _score("GENE", profile=profile, index=index, clock=clock, semantics=gene_semantics)
    assert match.matched_terms == (unplaceable,)
    assert match.breakdown.ic_missing_term_count == 1
    assert match.breakdown.matched_weight > 0.0
    assert "no information content" in match.rationale


# ---------------------------------------------------------------------------
# Provenance and evidence (GP-10 / GP-17 / GP-20)
# ---------------------------------------------------------------------------


def test_the_real_release_version_reaches_provenance_and_evidence(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """A score is only interpretable against a named release."""
    assert semantics.data_version.startswith("hp/releases/")
    provenance = semantics.provenance()
    assert provenance["hpo_data_version"] == semantics.data_version
    assert provenance["similarity_measure"] == SimilarityMeasure.LIN.value
    assert provenance["hpo_sha256"]
    assert provenance["annotation_sha256"]
    assert list(provenance) == sorted(provenance)

    profile = _profile((MICROCEPHALY, ObservationStatus.OBSERVED))
    index = _index(("GENE", MICROCEPHALY, "strong"))
    match = _score("GENE", profile=profile, index=index, clock=clock, semantics=semantics)
    payloads = [item.payload for item in match.evidence if "hpo_data_version" in item.payload]
    assert payloads and payloads[0]["hpo_data_version"] == semantics.data_version


def test_every_evidence_item_states_a_limitation_and_cites_the_measure(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """GP-17, plus the requirement that a published measure is named, not implied."""
    profile = _profile(
        (MICROCEPHALY, ObservationStatus.OBSERVED),
        (CATARACT, ObservationStatus.EXCLUDED),
        (SEIZURE, ObservationStatus.UNCERTAIN),
    )
    index = _index(
        ("GENE", MICROCEPHALY, "strong"),
        ("GENE", CATARACT, "moderate"),
        ("GENE", SEIZURE, "moderate"),
        ("GENE", NEPHROBLASTOMA, "supporting"),
    )
    match = _score("GENE", profile=profile, index=index, clock=clock, semantics=semantics)
    assert match.evidence
    for item in match.evidence:
        assert item.limitations.strip(), item.claim
        assert "hp/releases/" in item.method, item.claim
    joined = " ".join(item.method for item in match.evidence)
    assert "Lin (1998)" in joined
    assert "Köhler et al. (2009)" in joined


def test_exact_mode_evidence_declares_that_no_ontology_was_loaded(clock: Clock) -> None:
    """The fallback must not read as though it did the real comparison."""
    profile = _profile((CATARACT, ObservationStatus.OBSERVED))
    index = _index(("GENE", LENS, "strong"))
    match = _score("GENE", profile=profile, index=index, clock=clock, semantics=None)
    assert match.mode is ScoringMode.EXACT
    for item in match.evidence:
        assert "TD-04" in item.limitations


def test_a_contradictory_record_is_surfaced_in_the_ledger(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    profile = _profile(
        (EYE, ObservationStatus.EXCLUDED),
        (CATARACT, ObservationStatus.OBSERVED),
    )
    index = _index(("GENE", CATARACT, "strong"))
    match = _score("GENE", profile=profile, index=index, clock=clock, semantics=semantics)
    assert CATARACT in match.conflicting_terms
    assert any("internally inconsistent" in item.claim for item in match.evidence)
    # The directly recorded status still wins for the term it was recorded about.
    assert match.matched_terms == (CATARACT,)


def test_primary_and_alias_with_opposite_statuses_are_conflicting(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """A term asserted present under one spelling and absent under another.

    The scorer must not silently take a side. The term appears in
    ``conflicting_terms``, the ledger carries an item naming the record as
    internally inconsistent, and the term contributes to neither the matched nor
    the contradicted bucket — so it can neither support nor argue against the gene
    until a human resolves the record.

    Resolving this by precedence (the previous behaviour: OBSERVED outranks
    EXCLUDED) turned a contradiction into positive evidence *and* removed it from
    the conflict report, which is the inversion the four-valued logic exists to
    prevent.
    """
    profile = _profile(
        (URETER, ObservationStatus.EXCLUDED),
        (URETER_ALT_ID, ObservationStatus.OBSERVED),
        (MICROCEPHALY, ObservationStatus.OBSERVED),
    )
    index = _index(("GENE", URETER, "definitive"), ("GENE", MICROCEPHALY, "strong"))
    match = _score("GENE", profile=profile, index=index, clock=clock, semantics=semantics)

    assert URETER in match.conflicting_terms
    assert URETER not in match.matched_terms
    assert URETER not in match.contradicted_terms

    contradiction_items = [item for item in match.evidence if "inconsistent" in item.claim]
    assert len(contradiction_items) == 1
    assert URETER in contradiction_items[0].claim
    assert contradiction_items[0].limitations.strip()


def test_a_conflicted_term_moves_the_score_in_neither_direction(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """ "Contributes nothing" stated numerically, against both alternatives.

    The conflicted association must score exactly as if it had never been curated —
    strictly between the score it would get if OBSERVED won and the one it would
    get if EXCLUDED won, and equal to neither.
    """
    conflicted = _profile(
        (URETER, ObservationStatus.EXCLUDED),
        (URETER_ALT_ID, ObservationStatus.OBSERVED),
        (MICROCEPHALY, ObservationStatus.OBSERVED),
    )
    if_observed = _profile(
        (URETER, ObservationStatus.OBSERVED),
        (MICROCEPHALY, ObservationStatus.OBSERVED),
    )
    if_excluded = _profile(
        (URETER, ObservationStatus.EXCLUDED),
        (MICROCEPHALY, ObservationStatus.OBSERVED),
    )
    index = _index(("GENE", URETER, "definitive"), ("GENE", MICROCEPHALY, "strong"))
    without = _index(("GENE", MICROCEPHALY, "strong"))

    scores = {
        name: _score("GENE", profile=prof, index=idx, clock=clock, semantics=semantics).breakdown
        for name, prof, idx in (
            ("conflicted", conflicted, index),
            ("observed", if_observed, index),
            ("excluded", if_excluded, index),
        )
    }
    baseline = _score(
        "GENE", profile=conflicted, index=without, clock=clock, semantics=semantics
    ).breakdown

    # Neither side won: the conflicted term is out of the specificity ratio entirely.
    assert scores["conflicted"].specificity == baseline.specificity
    assert scores["conflicted"].matched_weight == baseline.matched_weight
    assert scores["conflicted"].contradicted_weight == baseline.contradicted_weight
    assert scores["conflicted"].specificity != scores["excluded"].specificity
    assert scores["conflicted"].matched_weight != scores["observed"].matched_weight


def test_an_unresolvable_profile_term_is_reported(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    profile = _profile(
        (MICROCEPHALY, ObservationStatus.OBSERVED),
        ("HP:9999999", ObservationStatus.OBSERVED),
    )
    index = _index(("GENE", MICROCEPHALY, "strong"))
    match = _score("GENE", profile=profile, index=index, clock=clock, semantics=semantics)
    assert match.unresolved_profile_terms == ("HP:9999999",)


def test_symmetric_similarity_is_reported_but_not_scored_on(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    """Both numbers are visible; the one that penalises unassessed features is not used.

    The gene here is curated with a feature nobody assessed. The symmetric
    best-match average drops because of it; the pipeline score does not.
    """
    profile = _profile((MICROCEPHALY, ObservationStatus.OBSERVED))
    lean = _index(("GENE", MICROCEPHALY, "strong"))
    padded = _index(("GENE", MICROCEPHALY, "strong"), ("GENE", SEIZURE, "strong"))
    thin = _score("GENE", profile=profile, index=lean, clock=clock, semantics=semantics)
    thick = _score("GENE", profile=profile, index=padded, clock=clock, semantics=semantics)

    assert thin.breakdown.symmetric_similarity is not None
    assert thick.breakdown.symmetric_similarity is not None
    assert thick.breakdown.symmetric_similarity < thin.breakdown.symmetric_similarity
    assert thick.score >= thin.score


# ---------------------------------------------------------------------------
# Ordering and determinism (GP-30)
# ---------------------------------------------------------------------------


def test_all_outputs_are_sorted(semantics: PhenotypeSemantics, clock: Clock) -> None:
    profile = _profile(
        (MICROCEPHALY, ObservationStatus.OBSERVED),
        (IUGR, ObservationStatus.OBSERVED),
        (DANDY_WALKER, ObservationStatus.OBSERVED),
        (EYE, ObservationStatus.EXCLUDED),
    )
    index = _index(
        ("GENE", NEPHROBLASTOMA, "supporting"),
        ("GENE", CATARACT, "moderate"),
        ("GENE", MICROCEPHALY, "strong"),
        ("GENE", IUGR, "strong"),
    )
    matches = score_all_genes(
        ["ZGENE", "GENE", "agene"],
        profile=profile,
        index=index,
        clock=clock,
        semantics=semantics,
    )
    assert list(matches) == ["agene", "GENE", "ZGENE"]
    match = matches["GENE"]
    for terms in (
        match.matched_terms,
        match.contradicted_terms,
        match.unassessed_terms,
        match.unexplained_terms,
        match.implied_matched_terms,
        match.implied_contradicted_terms,
    ):
        assert list(terms) == sorted(terms)
    assert [best.query for best in match.best_matches] == sorted(
        best.query for best in match.best_matches
    )


def test_repeat_scoring_in_one_process_is_identical(
    semantics: PhenotypeSemantics, clock: Clock
) -> None:
    profile = _profile(
        (MICROCEPHALY, ObservationStatus.OBSERVED),
        (EYE, ObservationStatus.EXCLUDED),
    )
    index = _index(("GENE", MICROCEPHALY, "strong"), ("GENE", CATARACT, "moderate"))
    first = _score("GENE", profile=profile, index=index, clock=clock, semantics=semantics)
    second = _score("GENE", profile=profile, index=index, clock=clock, semantics=semantics)
    assert first == second


#: Scored in a subprocess so the hash seed is genuinely different, not merely set.
_DETERMINISM_SCRIPT = textwrap.dedent(
    """
    from pathlib import Path

    from mva.clock import demo_clock
    from mva.determinism import stable_hash
    from mva.models.phenotype import (
        ObservationStatus,
        Onset,
        PhenotypeObservation,
        PhenotypeProfile,
    )
    from mva.phenotype.hpo import GeneAssociation, GenePhenotypeIndex
    from mva.phenotype.scoring import score_all_genes
    from mva.phenotype.semantics import HpoResourceSet

    fixtures = Path({fixtures!r})
    semantics = HpoResourceSet(
        ontology_path=fixtures / "hp_subset.obo",
        annotation_path=fixtures / "phenotype_subset.hpoa",
    ).load()

    def observation(term, status):
        return PhenotypeObservation(
            hpo_id=term,
            label="fixture term",
            status=status,
            onset=Onset.UNKNOWN,
            provenance="determinism_probe",
            extraction_confidence=1.0,
            source_excerpt_hash=None,
            notes=None,
        )

    profile = PhenotypeProfile(
        subject_id="PROBAND_TEST",
        observations=(
            observation("HP:0000252", ObservationStatus.OBSERVED),
            observation("HP:0001305", ObservationStatus.OBSERVED),
            observation("HP:0001511", ObservationStatus.OBSERVED),
            observation("HP:0000478", ObservationStatus.EXCLUDED),
            observation("HP:0001250", ObservationStatus.UNCERTAIN),
            observation("HP:0002667", ObservationStatus.NOT_ASSESSED),
        ),
        source_artifact="determinism_probe",
        hpo_version="hp/releases/fixture",
    )

    rows = [
        ("GENE_A", "HP:0000517", "strong"),
        ("GENE_A", "HP:0000252", "definitive"),
        ("GENE_A", "HP:0002667", "supporting"),
        ("GENE_B", "HP:0000518", "moderate"),
        ("GENE_B", "HP:0001250", "supporting"),
        ("GENE_C", "HP:0001305", "strong"),
        ("gene_c", "HP:0001511", "strong"),
    ]
    index = GenePhenotypeIndex(
        [
            GeneAssociation(
                gene_symbol=gene,
                hpo_id=term,
                label="fixture term",
                association_strength=strength,
                source="determinism_probe",
                version="v0",
            )
            for gene, term, strength in rows
        ],
        version="determinism-v0",
    )

    matches = score_all_genes(
        ["GENE_C", "GENE_A", "GENE_B", "GENE_D"],
        profile=profile,
        index=index,
        clock=demo_clock(),
        semantics=semantics,
    )
    payload = [
        {{
            "gene": gene,
            "score": f"{{match.score:.6f}}",
            "mode": match.mode.value,
            "annotation_status": match.annotation_status.value,
            "matched": list(match.matched_terms),
            "contradicted": list(match.contradicted_terms),
            "unassessed": list(match.unassessed_terms),
            "unexplained": list(match.unexplained_terms),
            "implied_matched": list(match.implied_matched_terms),
            "implied_contradicted": list(match.implied_contradicted_terms),
            "conflicts": list(match.conflicting_terms),
            "best_matches": [
                [best.query, best.target, f"{{best.similarity:.6f}}", best.mica]
                for best in match.best_matches
            ],
            "rationale": match.rationale,
            "evidence": [item.model_dump(mode="json") for item in match.evidence],
        }}
        for gene, match in matches.items()
    ]
    print(stable_hash({{"matches": payload, "provenance": dict(semantics.provenance())}}))
    """
)


def _run_scorer(seed: str) -> str:
    """Score in a fresh interpreter under an explicit ``PYTHONHASHSEED``."""
    env = {**os.environ, "PYTHONHASHSEED": seed}
    result = subprocess.run(  # noqa: S603
        [sys.executable, "-c", _DETERMINISM_SCRIPT.format(fixtures=str(HPO_FIXTURES))],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
        timeout=300,
    )
    assert result.returncode == 0, (
        f"determinism probe failed under PYTHONHASHSEED={seed}:\n{result.stderr[-2000:]}"
    )
    return result.stdout.strip()


@pytest.mark.slow
def test_scores_are_identical_under_different_hash_seeds() -> None:
    """GP-30, proved rather than argued.

    ``PYTHONHASHSEED`` changes the iteration order of every ``set`` and of every
    ``dict`` keyed on strings whose insertion order is set-derived. Running the
    scorer twice inside one interpreter cannot observe that; two interpreters with
    different seeds can. The comparison covers the scores, every term list, the
    best-match selection, the rationale text, the full evidence ledger and the
    provenance block — so an unstable argmax or an unsorted closure anywhere in the
    stack shows up as a different digest.
    """
    seeds: Sequence[str] = ("0", "1", "524287")
    digests = [_run_scorer(seed) for seed in seeds]
    assert all(digest for digest in digests), "probe produced no output"
    assert len(set(digests)) == 1, (
        "phenotype scoring is not deterministic across hash seeds: "
        + ", ".join(f"seed {seed} -> {digest}" for seed, digest in zip(seeds, digests, strict=True))
        + "\n\nGP-30 remediation: find the set/dict iteration, unstable argmax or "
        "float formatting that reached an output position."
    )
