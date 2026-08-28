"""Track 2: the literature-derived BUB1B mechanism and its drug triage.

Every other test in this repository runs against the *synthetic* demo tables, where
the answers were authored to make the machinery legible. This module runs the same
machinery against the **real, published biology of BUB1B** (`knowledge/literature/bub1b/`,
citations in that directory's `SOURCES.md`) and asserts the scientific conclusions
that `docs/track2-hypothesis.md` reports. It is the executable half of that document:
if the reasoning changes, this fails, and if this changes silently, the report is
lying.

Three claims are load-bearing enough to be worth stating in prose before the code:

1. **The naive top hit is rejected.** An MPS1/TTK inhibitor binds the exact kinase
   that builds the checkpoint named in the mechanism — maximal target proximity —
   and pushes it further in the disease direction. It must come out last, with
   `WRONG_DIRECTION`, not first (ASSUMPTION-DRUG-01).
2. **The highest-scoring agent is not the recommended one.** NMN has the only
   direction that *agrees* at the therapeutic target, and it is still rejected,
   because approval status and an unassessed oncogenic-risk question are gates that
   a score cannot buy its way past (ASSUMPTION-DRUG-03, ADR 0011).
3. **"Cannot determine" survives to the top of the accepted list.** The rank-1
   candidate has `directions_agree is None`. It must never be reported as agreement
   (ASSUMPTION-DRUG-02).

These tests contain no patient data and require none: the whole Track 2 argument is
built from public literature about a gene (GP-40, GP-44, ADR 0008).
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from mva.clock import demo_clock
from mva.determinism import stable_hash
from mva.interventions import (
    DrugCatalog,
    check_target_in_mechanism,
    generate_drug_hypotheses,
    is_chromosomal_instability_context,
    is_neurological_context,
    required_direction_for_node,
)
from mva.mechanisms import MechanismLibrary, build_mechanism, mechanism_relevance_score
from mva.models.drug import ApprovalStatus, DrugHypothesis, RejectionReason
from mva.models.mechanism import (
    UNSIGNED_DIRECTIONS,
    EffectDirection,
    MechanismHypothesis,
)

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
LITERATURE = REPO_ROOT / "knowledge" / "literature" / "bub1b"
CHAIN_PATH = LITERATURE / "mechanisms.tsv"
META_PATH = LITERATURE / "mechanism_meta.tsv"
CATALOG_PATH = LITERATURE / "drug_catalog.tsv"
SOURCES_PATH = LITERATURE / "SOURCES.md"

#: Version stamped onto evidence produced from these tables. Distinct from
#: `synthetic-v0.0` on purpose: a reader must be able to tell at a glance whether a
#: claim came from the fictional demo tables or from the literature (GP-20).
LIBRARY_VERSION = "bub1b-literature-v1"

#: The node a therapy has to move: BubR1 protein abundance. NOT the checkpoint node.
#: The checkpoint is where the disease is *named*; the protein is where the patient's
#: deviation is a quantity a drug could plausibly change, and it is upstream of BOTH
#: routes to missegregation (checkpoint and kinetochore attachment). See
#: ASSUMPTION-MECHANISM-06.
TARGET_NODE = "N3_bubr1"

#: The compensatory node: p53/p16 arrest, senescence and immune clearance of
#: missegregating cells. It deviates from wild type and is NOT a therapeutic target
#: (ASSUMPTION-MECHANISM-04).
COMPENSATORY_NODE = "N10_clearance"


# --------------------------------------------------------------------- fixtures


@pytest.fixture(scope="module")
def mechanism() -> MechanismHypothesis:
    library = MechanismLibrary.from_tsv(CHAIN_PATH, META_PATH, version=LIBRARY_VERSION)
    found = library.for_gene("BUB1B")
    assert found is not None, "the literature tables must carry a chain for BUB1B"
    return found


@pytest.fixture(scope="module")
def catalog() -> DrugCatalog:
    return DrugCatalog.from_tsv(CATALOG_PATH, version=LIBRARY_VERSION)


@pytest.fixture(scope="module")
def triage(mechanism: MechanismHypothesis, catalog: DrugCatalog) -> dict[str, DrugHypothesis]:
    result = generate_drug_hypotheses(mechanism=mechanism, catalog=catalog, clock=demo_clock())
    return {hypothesis.drug_id: hypothesis for hypothesis in result.all_hypotheses}


# ------------------------------------------------------------------- the chain


def test_chain_is_graded_link_by_link(mechanism: MechanismHypothesis) -> None:
    """ASSUMPTION-MECHANISM-01: a chain presented at uniform confidence is rhetoric.

    Two of the thirteen links are inferred rather than demonstrated, and the report
    names them. `L12` — that progenitor attrition is what produces growth restriction
    and microcephaly — is the weakest joint, and it is the one carrying the entire
    therapeutic rationale for the growth phenotype.
    """
    assert len(mechanism.nodes) == 12
    assert len(mechanism.links) == 13
    inferred = {link.link_id for link in mechanism.inferred_links}
    assert inferred == {"L11", "L12"}, (
        "The inferred links are the weak points this chain advertises. If the set "
        "changed, docs/track2-hypothesis.md names them explicitly and must change in "
        "the same commit."
    )
    assert not mechanism.is_fully_demonstrated
    assert all(link.uncertainty for link in mechanism.links), (
        "GP-17: every link states what is unresolved about it. An empty uncertainty "
        "cell is a gap in the record, not a statement of confidence."
    )


def test_therapeutic_target_is_the_protein_not_the_checkpoint(
    mechanism: MechanismHypothesis,
) -> None:
    """The target is chosen where the deviation is a movable quantity.

    Picking the checkpoint node instead would be picking the node that best *names*
    the disease. Every approved agent that acts there pushes it down, and BubR1
    reaches missegregation by a second, checkpoint-independent route through
    PP2A-B56 recruitment that a checkpoint-level intervention would not touch.
    """
    assert mechanism.therapeutic_target_node_id == TARGET_NODE
    target = mechanism.target_node()
    assert target.state_in_patient is EffectDirection.DECREASE
    assert target.deviation_is_pathological
    assert mechanism.required_correction is EffectDirection.INCREASE
    assert mechanism.disease_direction is EffectDirection.DECREASE


def test_clearance_node_is_compensatory_and_yields_no_corrective_direction(
    mechanism: MechanismHypothesis,
) -> None:
    """ASSUMPTION-MECHANISM-04: not every deviation from wild type is the disease.

    Arrest, senescence and immune clearance of missegregating cells deviate from
    wild type, and pushing them back towards it suppresses a protective response —
    in a child whose cancer risk is the thing that kills. The literature is
    genuinely split (the same p16 response drives the progeroid phenotype and
    suppresses tumours), which is precisely when the honest answer is to derive no
    sign at all.
    """
    node = next(n for n in mechanism.nodes if n.node_id == COMPENSATORY_NODE)
    assert not node.deviation_is_pathological
    assert required_direction_for_node(mechanism, COMPENSATORY_NODE) is EffectDirection.UNKNOWN


def test_disease_context_detectors_fire(mechanism: MechanismHypothesis) -> None:
    """Both context switches must trip, because both change which gaps are fatal.

    The chromosomal-instability context is what makes an unassessed oncogenic-risk
    answer blocking rather than merely noted (ADR 0011). The neurological context is
    what makes CNS penetration a recorded concern.
    """
    assert is_chromosomal_instability_context(mechanism)
    assert is_neurological_context(mechanism)


def test_mechanism_relevance_is_reported_not_flattered(mechanism: MechanismHypothesis) -> None:
    """The score is a claim about the chain's evidence, not about its truth."""
    score = mechanism_relevance_score(mechanism)
    assert 0.6 < score < 0.8, (
        "A well-sourced chain with two inferred links should land in the upper-middle "
        "of the range. A score at either extreme means either the grading became "
        "uniform (which ASSUMPTION-MECHANISM-01 forbids) or the chain lost its links."
    )


def test_building_the_mechanism_emits_falsification_experiments(
    mechanism: MechanismHypothesis,
) -> None:
    """One discriminating experiment per inferred link, each naming its alternative.

    An experiment that cannot fail against a stated alternative is a demonstration,
    not a test.
    """
    result = build_mechanism(
        "BUB1B",
        pair_id=None,
        library=MechanismLibrary([mechanism], version=LIBRARY_VERSION),
        clock=demo_clock(),
    )
    assert result.hypothesis is not None
    experiments = result.hypothesis.discriminating_experiments
    assert len(experiments) == len(mechanism.inferred_links)
    assert all(experiment.distinguishes_from for experiment in experiments)
    assert any("INFERRED" in warning for warning in result.warnings)
    assert any("developmental-window" in warning for warning in result.warnings)


# ------------------------------------------------------------ the direction gate


def test_the_naive_top_hit_is_rejected_for_direction(
    triage: dict[str, DrugHypothesis],
) -> None:
    """ASSUMPTION-DRUG-01, worked on real compounds.

    An MPS1 inhibitor binds the kinase that builds the checkpoint named in this
    mechanism. No search that ranks by target proximity can score it lower than
    first. It must come out rejected, with `WRONG_DIRECTION` as the reason a reader
    is shown first, because the follow-up for "wrong direction" is "never" and the
    follow-up for "not approved" is "wait".
    """
    mps1 = triage["DRUG-MPS1"]
    assert mps1.rejected
    assert mps1.directions_agree is False
    assert mps1.rejection_reasons[0] is RejectionReason.WRONG_DIRECTION
    assert mps1.worsens_chromosomal_instability is True
    assert mps1.score == min(hypothesis.score for hypothesis in triage.values()), (
        "The compound with the highest target proximity in the catalogue must also "
        "be the lowest-scoring. That inversion is the entire point of the stage."
    )


@pytest.mark.parametrize(
    "drug_id",
    ["DRUG-MPS1", "DRUG-HDAC", "DRUG-CQ", "DRUG-HSP90", "DRUG-METF"],
    ids=["mps1-inhibitor", "hdac-inhibitor", "autophagy-inhibitor", "hsp90-inhibitor", "ampk"],
)
def test_agents_pushing_the_disease_direction_are_disqualified(
    drug_id: str, triage: dict[str, DrugHypothesis]
) -> None:
    """Five real, mostly approved agents, all wrong-signed for the soma.

    Three of them (vorinostat, hydroxychloroquine, metformin) are approved and have
    paediatric exposure, which is exactly what makes them dangerous to a repurposing
    search that scores availability and tolerability but cannot see a sign. Two of
    them are *right* for the tumour on the same evidence that rejects them here.
    """
    hypothesis = triage[drug_id]
    assert hypothesis.directions_agree is False
    assert hypothesis.rejected
    assert RejectionReason.WRONG_DIRECTION in hypothesis.rejection_reasons


def test_direction_agreement_does_not_override_the_gates(
    triage: dict[str, DrugHypothesis],
) -> None:
    """The highest-scoring agent is rejected, and that is correct.

    NMN is the only entry whose observed direction *agrees* at the therapeutic
    target, on the strongest direct evidence in the catalogue. It is still rejected:
    it is not an approved medicine (ASSUMPTION-DRUG-03) and nobody has asked whether
    it changes aneuploidy or cancer susceptibility, which is blocking in this disease
    context (ASSUMPTION-DRUG-07, ADR 0011). A hypothesis can be right about the
    biology and still not be a repurposing candidate.
    """
    nmn = triage["DRUG-NMN"]
    assert nmn.directions_agree is True
    assert nmn.required_direction is EffectDirection.INCREASE
    assert nmn.observed_direction is EffectDirection.INCREASE
    assert nmn.rejected
    assert nmn.approval_status is not ApprovalStatus.APPROVED
    assert nmn.worsens_chromosomal_instability is None
    assert set(nmn.rejection_reasons) >= {
        RejectionReason.NOT_APPROVED,
        RejectionReason.ONCOGENIC_RISK,
    }
    assert nmn.score == max(hypothesis.score for hypothesis in triage.values()), (
        "NMN outscores every other entry and is still rejected. If a future change "
        "makes the top-scoring agent also the accepted one, check that a gate was "
        "not quietly demoted into a score term."
    )


def test_cannot_determine_reaches_the_top_of_the_accepted_list(
    triage: dict[str, DrugHypothesis],
) -> None:
    """ASSUMPTION-DRUG-02: `None` is not agreement, and it is not rejection either.

    Nicotinamide is simultaneously an NAD+ precursor — which should support the
    SIRT2-dependent deacetylation that slows BubR1 turnover — and a product
    inhibitor of the same sirtuin. The two pull opposite ways on one node, so no
    single sign can honestly be assigned. It ranks first anyway, carrying the flag,
    because it is the only entry whose cancer-susceptibility question has any answer
    at all.
    """
    accepted = [h for h in triage.values() if not h.rejected]
    assert [h.drug_id for h in accepted] == ["DRUG-NAM"]
    top = accepted[0]
    assert top.rank == 1
    assert top.directions_agree is None
    assert top.observed_direction is EffectDirection.CONTEXT_DEPENDENT
    assert top.observed_direction in UNSIGNED_DIRECTIONS
    assert RejectionReason.DIRECTION_UNKNOWN in top.concerns, (
        "An accepted candidate with an undetermined direction must carry the concern "
        "beside it. 'We accepted it despite X' is the sentence a reviewer needs."
    )
    assert top.score < 0.5, (
        "The top-ranked candidate scores below half. The ranking is telling the "
        "reader there is no good answer yet; a report that presents it as a "
        "recommendation is over-reading its own output."
    )


def test_no_accepted_hypothesis_is_wrong_signed(triage: dict[str, DrugHypothesis]) -> None:
    """The invariant the whole stage exists for, asserted over real compounds."""
    for hypothesis in triage.values():
        if not hypothesis.rejected:
            assert hypothesis.directions_agree is not False


def test_every_candidate_names_the_experiment_that_comes_first(
    triage: dict[str, DrugHypothesis],
) -> None:
    """A hypothesis with no disconfirming experiment is not a scientific claim."""
    for hypothesis in triage.values():
        assert len(hypothesis.proposed_validation_experiment) > 40, hypothesis.drug_id


def test_exposure_gap_is_recorded_for_every_candidate(
    triage: dict[str, DrugHypothesis],
) -> None:
    """ASSUMPTION-DRUG-05, and an honest statement about the field.

    No publication reports an effective concentration at the BubR1-abundance node
    for any of these agents, so every concentration cell is blank and every
    candidate carries the gap. Blank is unknown, not achievable.
    """
    for hypothesis in triage.values():
        assert hypothesis.pharmacokinetics.concentration_achievable is None, hypothesis.drug_id


def test_off_target_node_agents_are_flagged_as_mechanism_mismatch(
    triage: dict[str, DrugHypothesis], mechanism: MechanismHypothesis
) -> None:
    """ASSUMPTION-DRUG-04: acting somewhere on the chain is not correcting the mechanism.

    The statement is emitted from the node comparison, not from the curated class
    cell, so it cannot be edited away in the catalogue.
    """
    for hypothesis in triage.values():
        assert check_target_in_mechanism(hypothesis.target_node_id, mechanism), hypothesis.drug_id
        off_target = hypothesis.target_node_id != mechanism.therapeutic_target_node_id
        raised = set(hypothesis.concerns) | set(hypothesis.rejection_reasons)
        if off_target and not hypothesis.rejected:
            assert RejectionReason.MECHANISM_MISMATCH in raised, hypothesis.drug_id


def test_triage_is_deterministic(mechanism: MechanismHypothesis, catalog: DrugCatalog) -> None:
    """GP-30. Two runs of the same tables produce byte-identical conclusions."""
    first = generate_drug_hypotheses(mechanism=mechanism, catalog=catalog, clock=demo_clock())
    second = generate_drug_hypotheses(mechanism=mechanism, catalog=catalog, clock=demo_clock())
    assert stable_hash([h.model_dump(mode="json") for h in first.all_hypotheses]) == stable_hash(
        [h.model_dump(mode="json") for h in second.all_hypotheses]
    )
    assert stable_hash([e.model_dump(mode="json") for e in first.evidence]) == stable_hash(
        [e.model_dump(mode="json") for e in second.evidence]
    )


def test_rejections_are_kept_not_discarded(
    mechanism: MechanismHypothesis, catalog: DrugCatalog
) -> None:
    """GP-19. The rejection of a plausible, contraindicated agent is the finding.

    Twelve of thirteen entries are rejected. A stage that dropped them would report
    one candidate and no reasoning, which is worse than reporting nothing.
    """
    result = generate_drug_hypotheses(mechanism=mechanism, catalog=catalog, clock=demo_clock())
    assert len(result.rejected) == 12
    assert len(result.accepted) == 1
    assert len(result.rejection_record) == len(catalog), (
        "Every entry raised at least one reason, fatal or not, so every entry appears "
        "in the renderable ledger."
    )
    assert all(hypothesis.rejection_rationale for hypothesis in result.rejected)


# ------------------------------------------------------- citations resolve (GP-10)


def test_every_curated_row_is_cited_in_sources(catalog: DrugCatalog) -> None:
    """GP-10 applied to the knowledge tables themselves.

    A literature-derived table whose rows are not traceable to a named publication
    is indistinguishable from a fabricated one, and the whole reason these tables sit
    outside `knowledge/public/` is that they claim to be the former (ADR 0022).
    """
    sources = SOURCES_PATH.read_text(encoding="utf-8")
    chain = CHAIN_PATH.read_text(encoding="utf-8")

    link_ids = sorted(set(re.findall(r"\bL\d{2}\b", chain)))
    assert len(link_ids) == 13
    uncited_links = [link_id for link_id in link_ids if link_id not in sources]
    assert not uncited_links, (
        f"Mechanism links with no row in SOURCES.md: {uncited_links}.\n\n"
        "Remediation: add the link to the 'Mechanism chain — link to reference' table "
        "with the reference IDs behind it, or mark it INFERRED there. A link with no "
        "entry either way is an unattributed causal claim."
    )

    uncited_drugs = [entry.drug_id for entry in catalog.entries() if entry.drug_id not in sources]
    assert not uncited_drugs, (
        f"Catalogue rows with no entry in SOURCES.md: {uncited_drugs}.\n\n"
        "Remediation: add the row to the 'Drug catalogue — row to reference' table. "
        "Every direction cell is a pharmacological claim and needs a citation."
    )


def test_every_reference_id_used_is_defined(catalog: DrugCatalog) -> None:
    """No dangling `Rnn` pointers, in either direction.

    A citation that resolves to nothing reads as authority while pointing at
    nothing, which is the failure mode `tests/unit/test_docs_integrity.py` exists to
    prevent for `GP-nn` and `ADR nnnn`. The reference list is checked the same way.
    """
    assert catalog.version == LIBRARY_VERSION
    text = SOURCES_PATH.read_text(encoding="utf-8")
    marker = "## Mechanism chain"
    assert marker in text, "SOURCES.md must keep the reference list above the mapping tables"
    reference_list, usage = text.split(marker, 1)

    defined = {
        match.group(1)
        for line in reference_list.splitlines()
        if (match := re.match(r"\|\s*(R\d{2})\s*\|", line))
    }
    assert len(defined) >= 50, "the BUB1B reference list should not shrink silently"
    used = set(re.findall(r"\bR\d{2}\b", usage))
    dangling = sorted(used - defined)
    assert not dangling, (
        f"Reference IDs cited but not defined in the reference list: {dangling}.\n\n"
        "Remediation: add the publication to the reference list with its PMID, or fix "
        "the identifier. Never invent one."
    )

    unused = sorted(defined - used)
    assert not unused, (
        f"References defined but cited by no row: {unused}.\n\n"
        "Remediation: either cite it from a link or catalogue row, or remove it. A "
        "reference list padded with unused citations makes the sourcing look denser "
        "than it is."
    )
