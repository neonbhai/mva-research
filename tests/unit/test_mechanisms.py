"""Unit tests for the mechanism stage.

The chain is the artifact the rest of Track 2 reasons over, so these tests pin
the three things a downstream stage relies on: the shape of the chain, the signed
therapeutic requirement, and the honest advertisement of which links are only
inferred. A chain that quietly upgraded an inference to a demonstration would
propagate a false confidence into every drug hypothesis built on it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mva.clock import demo_clock
from mva.determinism import stable_hash
from mva.errors import IngestionError
from mva.mechanisms import (
    DEMONSTRATED_FRACTION_WEIGHT,
    INFERRED_LINK_PENALTY,
    LINK_STRENGTH_WEIGHT,
    MEAN_STRENGTH_WEIGHT,
    MechanismLibrary,
    MechanismResult,
    build_mechanism,
    mechanism_relevance_score,
)
from mva.models.base import AssertionTier
from mva.models.evidence import EvidenceStrength
from mva.models.mechanism import (
    EffectDirection,
    MechanismHypothesis,
    MechanismLink,
    MechanismNode,
    MechanismNodeKind,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
CHAIN_PATH = REPO_ROOT / "knowledge" / "public" / "mechanisms.tsv"
META_PATH = REPO_ROOT / "knowledge" / "public" / "mechanism_meta.tsv"
LIBRARY_VERSION = "synthetic-demo-2026.1"


def _library() -> MechanismLibrary:
    return MechanismLibrary.from_tsv(CHAIN_PATH, META_PATH, version=LIBRARY_VERSION)


@pytest.fixture
def library() -> MechanismLibrary:
    return _library()


@pytest.fixture
def mechanism(library: MechanismLibrary) -> MechanismHypothesis:
    found = library.for_gene("SYNTHKIN1")
    assert found is not None, "the demo library must carry the SYNTHKIN1 chain"
    return found


@pytest.fixture
def result(library: MechanismLibrary) -> MechanismResult:
    return build_mechanism("SYNTHKIN1", pair_id="PAIR-DEMO-01", library=library, clock=demo_clock())


# 1 --------------------------------------------------------------------------


@pytest.mark.unit
def test_library_loads_the_synthkin1_chain(mechanism: MechanismHypothesis) -> None:
    """The chain table is ragged by hand; the loader must still yield 7 nodes and 6 links."""
    assert mechanism.mechanism_id == "MECH-SYNTHKIN1-01"
    assert mechanism.gene_symbol == "SYNTHKIN1"
    assert len(mechanism.nodes) == 7
    assert len(mechanism.links) == 6
    assert [node.node_id for node in mechanism.nodes] == [
        "N1_variant",
        "N2_transcript",
        "N3_protein",
        "N4_complex",
        "N5_process",
        "N6_cellphen",
        "N7_organism",
    ]
    assert [link.link_id for link in mechanism.links] == ["L1", "L2", "L3", "L4", "L5", "L6"]
    # Every link must join two declared nodes, or the chain is not a chain.
    for link in mechanism.links:
        assert link.source_node_id in mechanism.node_ids
        assert link.target_node_id in mechanism.node_ids


@pytest.mark.unit
def test_nodes_and_links_are_fully_typed(mechanism: MechanismHypothesis) -> None:
    """GP-02/GP-16: no stringly-typed directions survive the boundary."""
    assert mechanism.nodes[0].kind is MechanismNodeKind.VARIANT
    assert mechanism.nodes[0].state_in_patient is EffectDirection.LOSS_OF_FUNCTION
    assert mechanism.nodes[-1].kind is MechanismNodeKind.ORGANISMAL_PHENOTYPE
    l6 = mechanism.links[-1]
    assert l6.tier is AssertionTier.INFERENCE
    assert l6.strength is EvidenceStrength.MODERATE
    assert all(link.uncertainty for link in mechanism.links), "GP-17: no link without a caveat"


# 2 --------------------------------------------------------------------------


@pytest.mark.unit
def test_therapeutic_requirement_is_signed_and_targeted(mechanism: MechanismHypothesis) -> None:
    """The two fields the drug stage checks against (GP-16)."""
    assert mechanism.required_correction is EffectDirection.RESTORE
    assert mechanism.therapeutic_target_node_id == "N5_process"
    assert mechanism.disease_direction is EffectDirection.LOSS_OF_FUNCTION
    target = mechanism.target_node()
    assert target.node_id == "N5_process"
    assert target.kind is MechanismNodeKind.CELLULAR_PROCESS
    assert mechanism.developmental_window_caveat, "the post-natal-relevance caveat must survive"


# 3 --------------------------------------------------------------------------


@pytest.mark.unit
def test_inferred_links_are_exactly_l6(mechanism: MechanismHypothesis) -> None:
    """L6 is the only link with is_directly_demonstrated=0 and must be named as such."""
    assert [link.link_id for link in mechanism.inferred_links] == ["L6"]
    assert mechanism.is_fully_demonstrated is False
    assert all(link.is_directly_demonstrated for link in mechanism.links[:-1])


# 4 --------------------------------------------------------------------------


@pytest.mark.unit
def test_relevance_score_penalises_inferred_links(mechanism: MechanismHypothesis) -> None:
    """More inferred links must strictly lower the score, never merely not raise it."""
    baseline = mechanism_relevance_score(mechanism)

    more_inferred = mechanism.model_copy(
        update={
            "links": tuple(
                link.model_copy(update={"is_directly_demonstrated": False})
                if link.link_id in {"L5", "L6"}
                else link
                for link in mechanism.links
            )
        }
    )
    fully_demonstrated = mechanism.model_copy(
        update={
            "links": tuple(
                link.model_copy(update={"is_directly_demonstrated": True})
                for link in mechanism.links
            )
        }
    )

    assert len(more_inferred.inferred_links) == 2
    assert mechanism_relevance_score(more_inferred) < baseline
    assert baseline < mechanism_relevance_score(fully_demonstrated)
    assert 0.0 <= mechanism_relevance_score(more_inferred) <= 1.0


@pytest.mark.unit
def test_relevance_score_of_a_chainless_hypothesis_is_zero(
    mechanism: MechanismHypothesis,
) -> None:
    """An unlinked node set asserts no mechanism and must not score as if it did."""
    assert mechanism_relevance_score(None) == 0.0
    assert mechanism_relevance_score(mechanism.model_copy(update={"links": ()})) == 0.0


# 5 --------------------------------------------------------------------------


@pytest.mark.unit
def test_target_node_raises_when_the_target_is_not_a_node() -> None:
    """A mechanism pointing at a node it does not contain is incoherent, not merely empty."""
    orphan = MechanismHypothesis(
        mechanism_id="MECH-ORPHAN-01",
        gene_symbol="SYNTHKIN1",
        summary="A chain whose therapeutic target was never declared as a node.",
        nodes=(
            MechanismNode(
                node_id="N1_variant",
                kind=MechanismNodeKind.VARIANT,
                label="Variant",
                state_in_patient=EffectDirection.LOSS_OF_FUNCTION,
            ),
        ),
        links=(),
        disease_direction=EffectDirection.LOSS_OF_FUNCTION,
        therapeutic_target_node_id="N_does_not_exist",
        required_correction=EffectDirection.RESTORE,
    )
    with pytest.raises(ValueError, match="not among its nodes"):
        orphan.target_node()


@pytest.mark.unit
def test_loader_rejects_a_target_that_is_not_a_node(tmp_path: Path) -> None:
    """The same incoherence is refused at load time rather than shipped downstream."""
    meta = tmp_path / "meta.tsv"
    original = META_PATH.read_text(encoding="utf-8")
    meta.write_text(original.replace("\tN5_process\t", "\tN_absent\t"), encoding="utf-8")
    with pytest.raises(IngestionError, match="therapeutic target"):
        MechanismLibrary.from_tsv(CHAIN_PATH, meta, version=LIBRARY_VERSION)


@pytest.mark.unit
def test_loader_rejects_a_chain_with_no_metadata(tmp_path: Path) -> None:
    """A chain with no required_correction cannot be direction-checked (GP-16)."""
    meta = tmp_path / "meta.tsv"
    lines = META_PATH.read_text(encoding="utf-8").splitlines()
    meta.write_text(
        "\n".join(line for line in lines if not line.startswith("MECH-")) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(IngestionError, match="no metadata row"):
        MechanismLibrary.from_tsv(CHAIN_PATH, meta, version=LIBRARY_VERSION)


@pytest.mark.unit
def test_loader_rejects_a_link_to_an_undeclared_node(tmp_path: Path) -> None:
    chain = tmp_path / "chain.tsv"
    chain.write_text(
        CHAIN_PATH.read_text(encoding="utf-8").replace(
            "\tN6_cellphen\tN7_organism\t", "\tN6_cellphen\tN9_ghost\t"
        ),
        encoding="utf-8",
    )
    with pytest.raises(IngestionError, match="N9_ghost"):
        MechanismLibrary.from_tsv(chain, META_PATH, version=LIBRARY_VERSION)


# 6 --------------------------------------------------------------------------


@pytest.mark.unit
def test_loading_is_deterministic() -> None:
    """GP-30: two loads of the same tables are byte-identical, including evidence IDs."""
    first, second = _library(), _library()
    assert stable_hash([m.model_dump(mode="json") for m in first.all_mechanisms()]) == stable_hash(
        [m.model_dump(mode="json") for m in second.all_mechanisms()]
    )

    run_one = build_mechanism("SYNTHKIN1", pair_id="PAIR-1", library=first, clock=demo_clock())
    run_two = build_mechanism("SYNTHKIN1", pair_id="PAIR-1", library=second, clock=demo_clock())
    assert run_one.warnings == run_two.warnings
    assert [e.evidence_id for e in run_one.evidence] == [e.evidence_id for e in run_two.evidence]
    assert stable_hash([e.model_dump(mode="json") for e in run_one.evidence]) == stable_hash(
        [e.model_dump(mode="json") for e in run_two.evidence]
    )


# builder ---------------------------------------------------------------------


@pytest.mark.unit
def test_build_mechanism_binds_the_pair_and_cites_everything(result: MechanismResult) -> None:
    """GP-10/GP-17: every emitted claim is an EvidenceItem and states its limitations."""
    assert result.hypothesis is not None
    assert result.hypothesis.pair_id == "PAIR-DEMO-01"
    assert result.evidence, "a built mechanism with no evidence is an unsourced claim"
    assert set(result.hypothesis.supporting_evidence_ids) == {
        item.evidence_id for item in result.evidence
    }
    for item in result.evidence:
        assert item.limitations.strip(), "GP-17"
        assert item.subject_id in {"MECH-SYNTHKIN1-01", "SYNTHKIN1"}
        assert item.tool_version


@pytest.mark.unit
def test_inferred_link_is_flagged_in_warnings_and_evidence(result: MechanismResult) -> None:
    """The weak joint is advertised in words, not only in a boolean field."""
    assert any("INFERRED" in warning and "L6" in warning for warning in result.warnings)
    assert any("developmental-window" in warning for warning in result.warnings)
    inferred_rows = [item for item in result.evidence if "INFERRED link" in item.limitations]
    assert len(inferred_rows) == 1
    assert "L6" in inferred_rows[0].claim


@pytest.mark.unit
def test_inferred_link_gets_a_discriminating_experiment(result: MechanismResult) -> None:
    """Each unproven link is paired with the experiment that would falsify it."""
    assert result.hypothesis is not None
    experiments = result.hypothesis.discriminating_experiments
    assert [e.experiment_id for e in experiments] == ["DX-MECH-SYNTHKIN1-01-L6"]
    assert experiments[0].expected_if_true != experiments[0].expected_if_false


@pytest.mark.unit
def test_unknown_gene_yields_absence_not_an_empty_mechanism(library: MechanismLibrary) -> None:
    """GP-14: no curated chain is a reportable absence, not a chain of length zero."""
    missing = build_mechanism("NOT_A_GENE", pair_id=None, library=library, clock=demo_clock())
    assert missing.hypothesis is None
    assert len(missing.evidence) == 1
    assert "Absence from this library is absence of curation" in missing.evidence[0].limitations
    assert missing.warnings


@pytest.mark.unit
def test_for_gene_is_case_insensitive_and_all_mechanisms_is_ordered(
    library: MechanismLibrary,
) -> None:
    assert library.for_gene("synthkin1") is not None
    assert library.for_gene("  SYNTHKIN1 ") is not None
    ids = [m.mechanism_id for m in library.all_mechanisms()]
    assert ids == sorted(ids)
    assert library.version == LIBRARY_VERSION


@pytest.mark.unit
def test_orphan_link_row_reports_the_missing_field(tmp_path: Path) -> None:
    """A half-filled row is an error: guessing its block would fabricate an edge."""
    chain = tmp_path / "chain.tsv"
    text = CHAIN_PATH.read_text(encoding="utf-8")
    broken = text.replace(
        "L4\tN4_complex\tN5_process\tfails_to_inhibit_anaphase_promoting_activity",
        "L4\tN4_complex\tN5_process\t",
    )
    chain.write_text(broken, encoding="utf-8")
    with pytest.raises(IngestionError):
        MechanismLibrary.from_tsv(chain, META_PATH, version=LIBRARY_VERSION)


@pytest.mark.unit
def test_relevance_score_formula_matches_its_documentation(
    mechanism: MechanismHypothesis,
) -> None:
    """GP-32: the published formula is the implemented one, so a weight change is visible."""
    links = mechanism.links
    mean_strength = sum(LINK_STRENGTH_WEIGHT[link.strength] for link in links) / len(links)
    n_inferred = len(mechanism.inferred_links)
    expected = (
        MEAN_STRENGTH_WEIGHT * mean_strength
        + DEMONSTRATED_FRACTION_WEIGHT * ((len(links) - n_inferred) / len(links))
        - INFERRED_LINK_PENALTY * n_inferred
    )
    assert mechanism_relevance_score(mechanism) == pytest.approx(round(expected, 6))


@pytest.mark.unit
def test_link_models_are_frozen(mechanism: MechanismHypothesis) -> None:
    """An upstream stage cannot mutate a chain another stage already scored."""
    link: MechanismLink = mechanism.links[0]
    with pytest.raises(ValueError, match="frozen"):
        link.is_directly_demonstrated = False  # type: ignore[misc]
