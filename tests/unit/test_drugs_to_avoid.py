"""``docs/track2-drugs-to-avoid.md`` must say what the pipeline says.

A rejection list written by hand agrees with the pipeline on the day it is typed
and never again. The page is therefore rendered from a live triage run, and this
file is the thing that makes "generated" mean something: if the committed bytes and
the generator disagree, the build is red.

Beyond sync, three properties are asserted because they are the reason the page
exists rather than incidental to it:

* the GP-14 distinction between OPPOSES (measured wrong) and never-measured is
  present, with both failure modes of confusing them stated;
* the counts in the prose are the pipeline's counts, not a number somebody typed;
* every agent got a *specific* per-disease warning, not a generic fallback. A page
  that degrades to "this agent was rejected" has stopped warning anybody and looks
  identical to one that is working.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from mva.config import find_repo_root

REPO_ROOT = find_repo_root(Path(__file__))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools.track2.drugs_to_avoid import (  # noqa: E402
    DANGER_BY_NODE_AND_DIRECTION,
    OUTPUT_PATH,
    AvoidEntry,
    DangerRuleMissingError,
    build_entries,
    render,
)

pytestmark = pytest.mark.unit


@pytest.fixture(scope="module")
def entries() -> tuple[AvoidEntry, ...]:
    return build_entries(REPO_ROOT)[1]


@pytest.fixture(scope="module")
def committed() -> str:
    return (REPO_ROOT / OUTPUT_PATH).read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def flat(committed: str) -> str:
    """The page with line breaks collapsed.

    The generator hard-wraps its prose, so a sentence this file asserts on is
    usually split across two lines. Asserting against the wrapped text would make
    every test here a test of the wrap width.
    """
    return " ".join(committed.split())


def test_the_committed_page_is_what_the_generator_produces(committed: str) -> None:
    """THE SYNC GATE. Every generated artifact needs one of these or it is not generated."""
    assert committed == render(REPO_ROOT), (
        f"{OUTPUT_PATH} is out of sync with its generator.\n\n"
        "Remediation: `just drugs-to-avoid`, then commit the result. Do not edit the "
        "page by hand -- it is rendered from the pipeline's own rejection output, and a "
        "hand-edit is a claim about a pipeline run that never happened."
    )


def test_rendering_twice_produces_identical_bytes() -> None:
    """GP-30. A generator with a wall-clock read in it has an unreviewable diff."""
    assert render(REPO_ROOT) == render(REPO_ROOT)


def test_every_rejected_agent_is_on_the_page(entries, committed: str) -> None:
    assert len(entries) == 12, (
        "docs/track2-hypothesis.md §5 reports 12 rejections. If that number changed, "
        "both documents change in the same commit."
    )
    for entry in entries:
        assert entry.name in committed
        assert entry.drug_id in committed
        assert entry.node_id in committed


def test_the_page_separates_measured_wrong_from_never_measured(entries, flat: str) -> None:
    """The one thing this page exists to get right (GP-14).

    ``OPPOSES`` is a measurement that came back the wrong way. ``oncogenic_risk``
    on a blank cell is nobody having asked. The pipeline writes the *same*
    rejection label for a blank cell as for a measured positive, so if this page
    does not pull them apart, nothing does.
    """
    opposes = [e for e in entries if e.opposes]
    never_measured = [e for e in entries if e.cin_state == "never_measured"]
    measured_positive = [e for e in entries if e.cin_state == "measured_positive"]

    assert opposes, "a page about wrong directions with no wrong directions on it"
    assert never_measured, "a page about unasked questions with no unasked questions"
    assert measured_positive, (
        "the measured-positive column is what makes the never-measured column mean "
        "something; with nothing in it the distinction is theoretical"
    )

    assert f"**OPPOSES — measured, and measured the wrong way.** {len(opposes)} of" in flat
    assert f"**Never measured.** {len(never_measured)} of" in flat
    assert f"**Measured, and the answer was yes.** {len(measured_positive)} agents" in flat
    assert "absence of information is not negative information" in flat
    assert "GP-14" in flat

    # Both failure modes of confusing them, not just the alarming one.
    assert "invents evidence that does not exist" in flat
    assert "Reading it as reassurance is worse" in flat


def test_measured_and_never_measured_are_not_rendered_with_the_same_cell(entries) -> None:
    """If the two states print identically, separating them in prose is decoration."""
    cells = {e.cin_state: e.cin_cell for e in entries}
    assert len(set(cells.values())) == len(cells)
    assert cells["never_measured"] != cells.get("measured_positive")


def test_the_direction_verdict_never_rounds_unknown_to_agreement(entries) -> None:
    """GP-16 / ASSUMPTION-DRUG-02, read off the rendered text rather than the field."""
    for entry in entries:
        if entry.directions_agree is None:
            assert entry.direction_verdict == "cannot determine (unsigned)"
            assert not entry.opposes
        elif entry.directions_agree is False:
            assert entry.opposes
            assert "OPPOSES" in entry.direction_verdict


def test_the_highest_scoring_agent_is_on_this_page(entries, flat: str) -> None:
    """A gate is a gate. A pipeline where the top score always wins cannot express
    "right about the biology, not yet a candidate"."""
    best = max(entries, key=lambda e: e.score)
    assert best.directions_agree is True
    assert best.cin_state == "never_measured"
    assert f"({best.name}, {best.score:.3f})" in flat


def test_every_agent_has_a_disease_specific_warning(entries) -> None:
    """No generic fallback sentence is allowed to reach the page."""
    for entry in entries:
        assert (entry.node_id, entry.observed) in DANGER_BY_NODE_AND_DIRECTION
        assert len(entry.danger) > 120, (
            f"{entry.name}'s warning is too short to say anything specific about this "
            f"disease: {entry.danger!r}"
        )


def test_an_uncovered_node_and_direction_fails_loudly() -> None:
    """The canary for the warning table.

    A new catalogue row landing on an uncovered node must stop the build, not
    quietly render a page with a blank where the reason should be.
    """
    from tools.track2.drugs_to_avoid import _danger

    with pytest.raises(DangerRuleMissingError, match="No danger rule for"):
        _danger("N_NOT_A_REAL_NODE", "increase", "Fabricatinib")


def test_the_page_says_it_is_generated_and_says_by_what(flat: str) -> None:
    assert "GENERATED FILE — do not edit" in flat
    assert "tools/track2/drugs_to_avoid.py" in flat
    assert "just drugs-to-avoid" in flat


def test_the_page_refuses_to_be_read_as_medical_advice(flat: str) -> None:
    """ "Avoid" is a statement about this pipeline, not about a patient."""
    assert "Nothing here is medical advice" in flat
    assert "never means" in flat and "useless in this disease" in flat
    assert "**tumour**-directed candidates" in flat
