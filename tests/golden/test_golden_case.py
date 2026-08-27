"""Golden expectation locks for the synthetic case.

These are the project's acceptance criteria expressed as tests. The expectations
live in `tests/golden/*.tsv` and are **never re-baselined to make a test pass**
(GP-32): changing one requires a decision record.

The point of a golden test here is not that the numbers are right — they are
uncalibrated heuristics — but that the *ordering the design intends* falls out of
general scoring rather than a special case.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pytest

from mva.config import CaseConfig, Workspace
from mva.models.pair import InheritanceModel
from mva.orchestrator import execute_pipeline

pytestmark = [pytest.mark.golden, pytest.mark.integration]


@pytest.fixture
def expected_ranking(
    golden: Callable[[str], list[dict[str, str]]],
) -> list[dict[str, str]]:
    return golden("expected_ranking.tsv")


@pytest.fixture
def expected_drugs(
    golden: Callable[[str], list[dict[str, str]]],
) -> list[dict[str, str]]:
    return golden("expected_drug_outcomes.tsv")


def test_synthetic_causal_pair_ranks_first(
    synthetic_config: CaseConfig,
    synthetic_workspace: Workspace,
    expected_ranking: list[dict[str, str]],
) -> None:
    """The known compound-heterozygous answer must rank #1.

    Nothing in the pipeline special-cases this pair. If this test fails, the
    scoring model has changed in a way that matters, and the correct response is a
    decision record explaining why — not an edit to expected_ranking.tsv.
    """
    result = execute_pipeline(synthetic_config, synthetic_workspace)
    assert result.ranked_pairs, "pipeline produced no ranked candidates"

    expected = expected_ranking[0]
    top = result.ranked_pairs[0]

    assert top.rank == 1
    assert top.gene_symbol == expected["gene_symbol"]
    assert set(top.variant_ids) == {expected["variant_a"], expected["variant_b"]}
    assert top.inheritance_model is InheritanceModel(expected["inheritance_model"])

    # Phase must remain UNKNOWN and be flagged: this is a proband-only VCF, and
    # claiming trans without evidence is the error GP-15 exists to prevent.
    must_flag = expected["must_be_flagged"]
    assert must_flag in top.flags, f"expected flag {must_flag!r} on the top candidate"


def test_top_candidate_carries_blocking_phase_question(
    synthetic_config: CaseConfig, synthetic_workspace: Workspace
) -> None:
    """An unresolved phase must surface as a blocking open question, not a footnote."""
    result = execute_pipeline(synthetic_config, synthetic_workspace)
    top = result.ranked_pairs[0]
    blocking = top.blocking_questions
    assert blocking, "top candidate has unknown phase but no blocking open question"
    assert any("trans" in q.question.lower() or "phase" in q.question.lower() for q in blocking)


def test_common_variant_pair_is_downranked_not_deleted(
    synthetic_config: CaseConfig, synthetic_workspace: Workspace
) -> None:
    """GP-13: the common-variant pair survives into the output, ranked below the answer."""
    result = execute_pipeline(synthetic_config, synthetic_workspace)
    common_ids = {
        "GRCh38:chr15:40205000:A:G",
        "GRCh38:chr15:40206000:T:C",
    }
    survivors = [p for p in result.ranked_pairs if common_ids & set(p.variant_ids)]
    assert survivors, "common variants were destroyed rather than down-ranked (GP-13)"
    assert all(p.rank is not None and p.rank > 1 for p in survivors)


def test_in_cis_pair_does_not_rank_first(
    synthetic_config: CaseConfig, synthetic_workspace: Workspace
) -> None:
    """An in-cis pair cannot be a compound heterozygote (GP-15)."""
    result = execute_pipeline(synthetic_config, synthetic_workspace)
    cis_ids = {"GRCh38:chr15:40201000:C:A", "GRCh38:chr15:40201050:T:G"}
    cis_pairs = [p for p in result.ranked_pairs if set(p.variant_ids) == cis_ids]
    assert cis_pairs, "the in-cis pair was destroyed rather than down-ranked"
    cis = cis_pairs[0]
    assert cis.phase_is_disqualifying
    assert cis.rank is not None and cis.rank > 1


def test_wrong_direction_drug_is_rejected(
    synthetic_config: CaseConfig,
    synthetic_workspace: Workspace,
    expected_drugs: list[dict[str, str]],
) -> None:
    """The headline Track 2 criterion.

    SYNTH-DRUG-B binds the correct target and pushes it the wrong way. It is the
    compound a naive target-proximity search ranks first, and it must be rejected.
    """
    result = execute_pipeline(synthetic_config, synthetic_workspace)
    assert result.mechanism is not None, "no mechanism was built for the top gene"

    # Assert against what the PIPELINE produced, not against the expectation file.
    # An earlier version of this test read both sides from expected_drugs and was
    # therefore a tautology: it stayed green with the direction check inverted.
    for row in expected_drugs:
        drug_id = row["drug_id"]
        hypothesis = result.drug(drug_id)
        assert hypothesis is not None, (
            f"{drug_id} was in the catalogue but the pipeline produced no hypothesis "
            "for it — a candidate must never be silently dropped (GP-19)"
        )

        outcome = row["expected_outcome"]
        if outcome == "rejected":
            assert hypothesis.rejected, (
                f"{drug_id} ({hypothesis.name}) was expected to be REJECTED but was accepted"
            )
            reasons = [r.value for r in hypothesis.rejection_reasons]
            expected_reason = row["expected_primary_reason"]
            assert reasons, f"{drug_id} is rejected but records no reason"
            assert reasons[0] == expected_reason, (
                f"{drug_id} rejected for {reasons[0]!r}, expected {expected_reason!r} "
                f"as the leading reason (all reasons: {reasons})"
            )
        else:
            assert not hypothesis.rejected, (
                f"{drug_id} was expected to be accepted ({outcome}) but was rejected "
                f"for {[r.value for r in hypothesis.rejection_reasons]}"
            )

    # The headline criterion, stated directly against the result.
    synthinib = result.drug("SYNTH-DRUG-B")
    assert synthinib is not None
    assert synthinib.rejected
    assert synthinib.directions_agree is False, (
        "the wrong-direction agent must be recorded as DISAGREEING, not as "
        "undetermined — 'cannot tell' and 'disagrees' are different findings"
    )
    assert "wrong_direction" in [r.value for r in synthinib.rejection_reasons]


def test_demo_artifacts_are_all_produced(
    synthetic_config: CaseConfig, synthetic_workspace: Workspace
) -> None:
    """Every artifact `just demo` promises must actually exist."""
    result = execute_pipeline(synthetic_config, synthetic_workspace)
    produced = {a.kind.value for a in result.manifest.artifacts}
    required = {
        "run_manifest",
        "normalised_variants",
        "annotated_variants",
        "candidate_pairs",
        "evidence_db",
        "submission",
        "dossier",
        "mechanism_report",
        "drug_hypotheses",
        "rejection_record",
        "track2_report",
        "provenance_manifest",
    }
    missing = required - produced
    # run_manifest is embedded in provenance_manifest rather than written twice.
    missing.discard("run_manifest")
    assert not missing, f"demo did not produce: {sorted(missing)}"


def test_repeat_run_is_byte_identical(
    synthetic_config: CaseConfig, synthetic_workspace: Workspace, tmp_path: Path
) -> None:
    """GP-30: determinism is an acceptance criterion, not an aspiration."""
    from mva.pipeline import verify_determinism

    first = execute_pipeline(synthetic_config, synthetic_workspace)
    second = execute_pipeline(synthetic_config, synthetic_workspace)
    identical, differences = verify_determinism(first.digest, second.digest)
    assert identical, "repeat run produced different artifacts:\n" + "\n".join(differences)


def test_direction_check_is_load_bearing(
    synthetic_config: CaseConfig, synthetic_workspace: Workspace
) -> None:
    """Guard against the direction gate becoming a no-op.

    A previous version of the Track 2 acceptance test passed with
    `agrees is False` inverted to `True`, because it only ever compared the
    expectation file against itself. This asserts the distribution the gate
    actually produces: at least one drug disagrees, at least one agrees, and at
    least one is undeterminable — so a gate that collapsed to a constant would
    fail here regardless of which constant it collapsed to.
    """
    result = execute_pipeline(synthetic_config, synthetic_workspace)
    verdicts = [d.directions_agree for d in (*result.accepted_drugs, *result.rejected_drugs)]
    assert verdicts, "no drug hypotheses were produced at all"
    assert False in verdicts, "no drug DISAGREES; the direction check may be a no-op"
    assert True in verdicts, "no drug AGREES; the direction check may be inverted"
    assert None in verdicts, (
        "no drug has an undeterminable direction; the tri-state may have collapsed "
        "to a boolean, which is how 'cannot tell' becomes 'agrees'"
    )


def test_no_accepted_drug_disagrees_with_the_required_direction(
    synthetic_config: CaseConfig, synthetic_workspace: Workspace
) -> None:
    """The central Track 2 claim, asserted over the whole accepted set."""
    result = execute_pipeline(synthetic_config, synthetic_workspace)
    for drug in result.accepted_drugs:
        assert drug.directions_agree is not False, (
            f"{drug.drug_id} ({drug.name}) acts {drug.observed_direction.value!r} "
            f"where {drug.required_direction.value!r} is required, yet was accepted. "
            "A wrong-direction agent must never reach the accepted set (GP-16)."
        )
