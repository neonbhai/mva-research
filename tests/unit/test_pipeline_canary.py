"""The positive control for Track 2, and the tests that keep the control honest.

Every other test in `tests/unit/test_track2_bub1b.py` asserts a *rejection*: the
MPS1 inhibitor comes last, the top-scoring agent does not win, "cannot determine"
is never rounded up. Those are the right assertions and they share one blind spot.
A triage stage that rejected every input would pass all of them, and would score
identically to a careful one.

`mva.interventions.canaries` closes that. It is the science-side counterpart of
`mva.privacy.patterns.detector_canaries`: fabricated cases whose verdict is not in
question, pushed through the real `generate_drug_hypotheses`, with the answer
demanded rather than observed.

The fixtures are deliberately not biology. `CANARYGENE1`, `Canarystatin` and the
rest are shaped to be unambiguous to the pipeline, and `test_the_fixtures_are_never
_presented_as_biology` is what keeps them from leaking anywhere a reader could
mistake them for a finding (GP-20).
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest

from mva.config import find_repo_root
from mva.interventions.canaries import (
    CANARY_VERSION,
    CanaryOutcome,
    PipelineCanary,
    check_triage_liveness,
    pipeline_canaries,
    run_canary,
)
from mva.models.drug import RejectionReason

pytestmark = pytest.mark.unit

REPO_ROOT = find_repo_root(Path(__file__))


@pytest.fixture(scope="module")
def canaries() -> dict[str, PipelineCanary]:
    return {canary.canary_id: canary for canary in pipeline_canaries()}


@pytest.fixture(scope="module")
def outcomes() -> dict[str, CanaryOutcome]:
    return {outcome.canary_id: outcome for outcome in check_triage_liveness().outcomes}


def test_the_triage_stage_still_reaches_every_demanded_verdict() -> None:
    """THE GATE. Runs all six controls through the real triage function."""
    liveness = check_triage_liveness()
    assert liveness.ok, liveness.summary()


def test_the_pipeline_can_still_say_yes(outcomes: dict[str, CanaryOutcome]) -> None:
    """The single most important assertion here.

    Without it, deleting the entire direction check and rejecting unconditionally
    would leave this repository's Track 2 suite green.
    """
    accept = outcomes["ACCEPT"]
    assert accept.accepted is True, (
        "The triage stage rejected an agent that is approved, has paediatric "
        "exposure, reaches its target at a fraction of its achievable exposure, has "
        "direct human-trial evidence, was measured to push the therapeutic target in "
        "exactly the required direction, and whose chromosomal-instability risk was "
        "assessed and came back negative. If this is rejected, nothing can be "
        "accepted, and every 'correctly rejected' test in this suite is vacuous."
    )
    assert accept.directions_agree is True
    assert accept.primary_reason is None
    assert accept.concerns == ()
    assert accept.score == pytest.approx(1.0), (
        "The unambiguous case should also score unambiguously. A ceiling that the "
        "best possible candidate cannot reach means the scale is not calibrated to "
        "anything."
    )


def test_at_least_one_canary_demands_acceptance(canaries: dict[str, PipelineCanary]) -> None:
    """Guards the control set itself.

    A canary suite that drifted into all-rejections would be indistinguishable from
    no canary suite at all, and would pass `test_the_triage_stage_still_reaches_
    every_demanded_verdict` while proving nothing.
    """
    accepting = [c for c in canaries.values() if c.expect_accepted]
    assert accepting, (
        "Every canary now expects a rejection. That is the exact failure this module "
        "exists to prevent: a pipeline that says no to everything passes a suite that "
        "only ever asks it to say no."
    )


def test_inverting_the_direction_rejects_for_the_direction(
    outcomes: dict[str, CanaryOutcome],
) -> None:
    """Rejected is not enough; rejected *for the right reason* is the assertion.

    An agent rejected for approval status while its sign check was broken would
    satisfy a naive "it was rejected" test and hide a dead gate.
    """
    invert = outcomes["INVERT"]
    assert invert.accepted is False
    assert invert.primary_reason is RejectionReason.WRONG_DIRECTION
    assert invert.directions_agree is False


def test_an_unsigned_direction_is_neither_accepted_blindly_nor_rejected(
    outcomes: dict[str, CanaryOutcome],
) -> None:
    """GP-16 and GP-14, on the same candidate.

    Wrong direction disqualifies. Unknown direction does not — it is a flag and a
    penalty. Rejecting this candidate would mean treating "nobody measured it" as
    "it points the wrong way", which is the GP-14 error running in its quieter
    direction, and it is just as wrong as the loud one.
    """
    unknown = outcomes["UNKNOWN_DIRECTION"]
    assert unknown.accepted is True
    assert unknown.directions_agree is None
    assert RejectionReason.DIRECTION_UNKNOWN in unknown.concerns, (
        "Accepted with the flag dropped is worse than rejected: the report then "
        "presents an unsigned direction as an unremarkable one."
    )
    assert unknown.score < outcomes["ACCEPT"].score, (
        "An unsigned direction must cost something. If it is free, the ranking "
        "cannot distinguish a measured agent from an unmeasured one."
    )


def test_the_unasked_safety_question_rejects_an_agent_whose_direction_agrees(
    outcomes: dict[str, CanaryOutcome],
) -> None:
    """GP-14 / ADR 0011, made executable.

    This is the whole Track 2 finding in one control: the agent's biology is right
    and it is still rejected, because in a chromosomal-instability disorder the
    question nobody asked is disqualifying. Crucially the direction verdict must
    stay AGREES — a pipeline that quietly downgraded the direction to justify the
    rejection would be inventing a biological objection it does not have.
    """
    unassessed = outcomes["UNASSESSED_RISK"]
    assert unassessed.accepted is False
    assert unassessed.primary_reason is RejectionReason.ONCOGENIC_RISK
    assert unassessed.directions_agree is True


def test_an_off_chain_target_is_a_finding_not_a_non_event(
    outcomes: dict[str, CanaryOutcome],
) -> None:
    """The negative control. Rules out "accepted" being an artefact of ignoring the
    mechanism, and pins that an off-chain node yields no direction verdict at all
    rather than a fabricated disagreement."""
    off_chain = outcomes["OFF_CHAIN"]
    assert off_chain.accepted is False
    assert off_chain.primary_reason is RejectionReason.TARGET_NOT_IN_MECHANISM
    assert off_chain.directions_agree is None


def test_the_same_agent_flips_verdict_against_a_sibling_gene(
    canaries: dict[str, PipelineCanary], outcomes: dict[str, CanaryOutcome]
) -> None:
    """The sibling-case control, and the sharpest one here.

    The verdict has to be a property of the (mechanism, drug) pair. A stage that had
    somehow memorised "this compound is good" passes every other canary in this file
    and fails only this one: the identical agent, unchanged, is rejected against a
    second synthetic gene whose disease direction is the mirror image.
    """
    assert canaries["SIBLING_GENE"].entry.drug_id == canaries["ACCEPT"].entry.drug_id
    assert (
        canaries["SIBLING_GENE"].mechanism.gene_symbol != canaries["ACCEPT"].mechanism.gene_symbol
    )
    assert (
        canaries["SIBLING_GENE"].mechanism.required_correction
        != canaries["ACCEPT"].mechanism.required_correction
    )
    sibling = outcomes["SIBLING_GENE"]
    assert sibling.accepted is False
    assert sibling.primary_reason is RejectionReason.WRONG_DIRECTION
    assert outcomes["ACCEPT"].accepted is True


def test_each_control_varies_exactly_one_thing(canaries: dict[str, PipelineCanary]) -> None:
    """Control hygiene. A canary rejected for three reasons at once proves none of them.

    Every catalogue-level control is the ACCEPT baseline with a single field moved,
    ignoring the identifying strings. If two fields differ, the canary no longer
    isolates the axis it claims to test.
    """
    baseline = canaries["ACCEPT"].entry
    identity = {"drug_id", "name", "approved_name", "target"}
    for canary_id in ("INVERT", "UNKNOWN_DIRECTION", "UNASSESSED_RISK", "OFF_CHAIN"):
        entry = canaries[canary_id].entry
        differing = {
            field.name
            for field in dataclasses.fields(baseline)
            if field.name not in identity
            and getattr(entry, field.name) != getattr(baseline, field.name)
        }
        assert len(differing) == 1, (
            f"{canary_id} differs from the ACCEPT baseline in {sorted(differing)}. A "
            f"control that moves more than one variable cannot attribute its verdict "
            f"to any of them."
        )


def test_running_a_canary_twice_gives_the_same_verdict(canaries: dict[str, PipelineCanary]) -> None:
    """GP-30. The canaries use a fixed clock, so they are a determinism probe too."""
    for canary in canaries.values():
        assert run_canary(canary) == run_canary(canary)


def test_the_fixtures_are_never_presented_as_biology(canaries: dict[str, PipelineCanary]) -> None:
    """GP-20. Obviously fabricated, and labelled so at every level."""
    assert "SYNTHETIC" in CANARY_VERSION
    for canary in canaries.values():
        mechanism = canary.mechanism
        assert mechanism.gene_symbol.startswith("CANARYGENE")
        assert "SYNTHETIC" in mechanism.summary
        assert all("SYNTHETIC" in node.label for node in mechanism.nodes)
        assert all("SYNTHETIC" in link.uncertainty for link in mechanism.links)
        assert canary.entry.drug_id.startswith("DRUG-CANARY")
        assert canary.entry.name.startswith("Canary")
        assert "SYNTHETIC" in canary.entry.notes
        assert "SYNTHETIC" in canary.entry.validation_experiment


def test_the_fixtures_never_leak_into_a_knowledge_table_or_a_report(
    canaries: dict[str, PipelineCanary],
) -> None:
    """A synthetic substitute is never described as biologically valid (GP-20).

    The strongest form of that rule is that these names cannot appear anywhere a
    reader would take them for a result. They live in one module and in this test.
    """
    tokens = {"CANARYGENE1", "CANARYGENE2", "Canarystatin", "Canaryase", "MECH-CANARY"}
    allowed = {
        REPO_ROOT / "src" / "mva" / "interventions" / "canaries.py",
        Path(__file__),
    }
    offenders: list[str] = []
    for root in (REPO_ROOT / "knowledge", REPO_ROOT / "docs", REPO_ROOT / "config"):
        for path in root.rglob("*"):
            if not path.is_file() or path in allowed:
                continue
            try:
                text = path.read_text(encoding="utf-8")
            except (OSError, UnicodeDecodeError):
                continue
            hits = sorted(token for token in tokens if token in text)
            if hits:
                offenders.append(f"  {path.relative_to(REPO_ROOT)}: {hits}")
    assert not offenders, (
        "Pipeline-canary fixtures appear outside the module that defines them:\n"
        + "\n".join(offenders)
        + "\n\nThese are fabricated agents and fabricated genes. Anywhere they can be "
        "read as a finding, they are misinformation (GP-20)."
    )
