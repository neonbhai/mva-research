"""``docs/track2-ranking.md`` must say what the pipeline says.

The rejections already had a page that could not drift (`test_drugs_to_avoid.py`). The
acceptance did not: "1 accepted, nicotinic acid, rank 1, 0.500" lived in
`docs/track2-hypothesis.md` §5 as hand-typed prose, which agrees with the pipeline on
the day it is typed and never again. That was the single most quotable number in the
repository and the one number no command could reproduce. The page is therefore
rendered from a live triage run, and this file is the thing that makes "generated" mean
something: if the committed bytes and the generator disagree, the build is red.

Beyond sync, five properties are asserted because they are the reason the page exists
rather than incidental to it:

* the accepted score is never printed without the unmeasured direction beside it,
  checked block by block over the rendered page rather than taken on trust;
* the hand-written §5 table and the live run agree on every one of the thirteen scores,
  and a disagreement is reported rather than reconciled;
* `unknown` is kept apart from both agreement and disagreement (GP-16), and a blank
  chromosomal-instability cell is kept apart from a measured one (GP-14);
* the catalogue's provenance is read off the loader's own flag, so a synthetic
  substitute can never be described as literature-curated (GP-20);
* the published scoring weights are the pipeline's weights, so the page cannot keep
  explaining a score the code no longer computes.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

from mva.config import find_repo_root
from mva.interventions import (
    APPROVAL_WEIGHT,
    DIRECTION_UNKNOWN_CREDIT,
    DIRECTION_WEIGHT,
    EVIDENCE_WEIGHT,
    SAFETY_WEIGHT,
)
from mva.models.drug import RejectionReason

REPO_ROOT = find_repo_root(Path(__file__))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools.track2.ranking import (  # noqa: E402
    CONCERN_CLAUSE,
    OUTPUT_PATH,
    PROBE_RATIONALE_BY_NODE_AND_DIRECTION,
    UNMEASURED_DIRECTION_MARKERS,
    AcceptedProbe,
    ConcernClauseMissingError,
    ProbeRationaleMissingError,
    Ranking,
    build_ranking,
    render,
)

pytestmark = pytest.mark.unit

#: The hand-written scientific case. Read here, never written: if it and a live run
#: disagree about a score, that is a finding to report, not a file to edit.
HYPOTHESIS_DOC = "docs/track2-hypothesis.md"


@pytest.fixture(scope="module")
def ranking() -> Ranking:
    return build_ranking(REPO_ROOT)


@pytest.fixture(scope="module")
def probe(ranking: Ranking) -> AcceptedProbe:
    return ranking.accepted[0]


@pytest.fixture(scope="module")
def committed() -> str:
    return (REPO_ROOT / OUTPUT_PATH).read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def flat(committed: str) -> str:
    """The page with line breaks collapsed.

    The generator hard-wraps its prose, so a sentence this file asserts on is usually
    split across two lines. Asserting against the wrapped text would make every test
    here a test of the wrap width.
    """
    return " ".join(committed.split())


@pytest.fixture(scope="module")
def blocks(committed: str) -> tuple[str, ...]:
    """The page split into the smallest units a claim can sit in, each unwrapped.

    A paragraph, a single bullet with its continuation lines, one table row, one
    heading. Finer than paragraphs on purpose: a whole bullet list treated as one
    block would let a score in one bullet be excused by a caveat in a different one,
    which is exactly the adjacency the score rule is supposed to guarantee.
    """
    out: list[str] = []
    buffer: list[str] = []

    def flush() -> None:
        if buffer:
            out.append(" ".join(" ".join(buffer).split()))
            buffer.clear()

    for line in committed.splitlines():
        stripped = line.strip()
        if not stripped:
            flush()
        elif line.startswith("|") or stripped.startswith("#"):
            flush()
            out.append(" ".join(stripped.split()))
        elif stripped.startswith("- "):
            flush()
            buffer.append(stripped)
        else:
            buffer.append(stripped)
    flush()
    return tuple(out)


def test_the_committed_page_is_what_the_generator_produces(committed: str) -> None:
    """THE SYNC GATE. Every generated artifact needs one of these or it is not generated."""
    assert committed == render(REPO_ROOT), (
        f"{OUTPUT_PATH} is out of sync with its generator.\n\n"
        "Remediation: `just ranking`, then commit the result. Do not edit the page by "
        "hand -- it is rendered from the pipeline's own ranked output, and a hand-edit "
        "is a claim about a pipeline run that never happened."
    )


def test_rendering_twice_produces_identical_bytes() -> None:
    """GP-30. A generator with a wall-clock read in it has an unreviewable diff."""
    assert render(REPO_ROOT) == render(REPO_ROOT)


def test_the_pipeline_still_accepts_something(ranking: Ranking) -> None:
    """The positive control for the page's whole subject.

    A triage stage that rejected everything would satisfy every rejection test in the
    suite, and this page would have nothing to render. `build_ranking` refuses to
    produce an empty ranking rather than emitting a page that asserts a result which
    does not exist; this is that refusal, asserted from the outside.
    """
    assert ranking.accepted, "no accepted candidate: see `just science-canary`"
    assert len(ranking.rows) == 13, (
        "docs/track2-hypothesis.md §5 reports 13 agents screened. If that number "
        "changed, both documents change in the same commit."
    )
    assert [probe.row.rank for probe in ranking.accepted] == [1]


def test_the_accepted_score_is_never_printed_without_its_unmeasured_direction(
    ranking: Ranking, blocks: tuple[str, ...]
) -> None:
    """The one thing this page exists to get right.

    "Nicotinic acid, rank 1, 0.500" reads as a finding about a drug and is a finding
    about a pipeline. The direction of effect on BubR1 abundance has never been
    measured in any genotype, so a score quoted on its own silently answers the
    question the hypothesis is asking. Every block carrying the score must therefore
    carry the sign's status too -- in the table row, in the bullet and in the prose.
    """
    score = f"{ranking.accepted[0].row.score:.3f}"
    carrying = [block for block in blocks if score in block]

    assert carrying, (
        f"the accepted score {score} does not appear on the page at all, which means "
        f"the headline result is still not reproducible from this document"
    )
    naked = [
        block
        for block in carrying
        if not any(marker in block for marker in UNMEASURED_DIRECTION_MARKERS)
    ]
    assert not naked, (
        f"the accepted score {score} is printed without the unmeasured direction "
        f"beside it, in:\n\n"
        + "\n\n".join(f"  {block}" for block in naked)
        + "\n\nRemediation: say in the SAME block that the direction is unmeasured, "
        f"using one of {list(UNMEASURED_DIRECTION_MARKERS)}. The score is the number a "
        "reader quotes; detached from the sign it reports a result nobody has produced."
    )


def test_the_live_run_and_the_hand_written_case_agree_on_every_score(ranking: Ranking) -> None:
    """Cross-document check against the one document that is NOT generated.

    `docs/track2-hypothesis.md` §5 prints the same thirteen agents and scores as prose
    a person typed. This asserts the live run still produces each of them.

    If this fails, the two documents disagree about what the pipeline produced. That is
    a real finding about the repository and the remediation is NOT to edit either number
    into agreement -- the generated page would then carry a score the pipeline did not
    compute, or the scientific case would carry one it did not either.
    """
    text = (REPO_ROOT / HYPOTHESIS_DOC).read_text(encoding="utf-8")
    lines = text.splitlines()
    missing: list[str] = []
    for row in ranking.rows:
        score = f"{row.score:.3f}"
        if not any(row.name in line and score in line for line in lines):
            missing.append(f"  {row.name} — live run scores {score}")

    assert not missing, (
        f"{HYPOTHESIS_DOC} and a live triage run disagree about these agents:\n"
        + "\n".join(missing)
        + f"\n\nSTOP AND REPORT THIS. Do not change either file to match the other. "
        f"{HYPOTHESIS_DOC} is the hand-written scientific case and {OUTPUT_PATH} is "
        "rendered from the pipeline; a disagreement means the catalogue, a weight or a "
        "gate moved underneath a published claim, and which of the two is now wrong is "
        "a review decision rather than a text edit."
    )


def test_the_three_direction_states_render_distinctly(ranking: Ranking) -> None:
    """GP-16 / ASSUMPTION-DRUG-02, read off the rendered text rather than the field.

    `unknown` is neither agreement nor disagreement. If the three states print the
    same words, the tri-state is preserved in the code and lost in the document, which
    is the only place a reader looks.
    """
    verdicts = {row.directions_agree: row.direction_verdict for row in ranking.rows}
    assert set(verdicts) == {True, False, None}, (
        "this run no longer exercises all three direction states, so the page can no "
        "longer demonstrate that it keeps them apart"
    )
    assert len(set(verdicts.values())) == 3
    assert "OPPOSES" in verdicts[False]
    assert verdicts[None] == "cannot determine (unsigned)"
    assert "agrees" in verdicts[True] and "OPPOSES" not in verdicts[True]

    for row in ranking.rows:
        if row.directions_agree is None:
            assert row.direction_verdict == "cannot determine (unsigned)"


def test_a_blank_chromosomal_instability_cell_never_renders_as_measured_harm(
    ranking: Ranking, flat: str
) -> None:
    """GP-14. Nobody having asked is not a finding, and must not print like one."""
    cells = {row.cin_state: row.cin_cell for row in ranking.rows}
    assert set(cells) == {"measured_positive", "measured_negative", "never_measured"}
    assert len(set(cells.values())) == 3
    assert cells["never_measured"] != cells["measured_positive"]

    assert "GP-14" in flat
    assert "absence of a finding, not a finding of harm" in flat
    assert "Do not read the rejected rows as a list of harmful agents" in flat


def test_the_accepted_candidate_is_not_the_top_scoring_one(ranking: Ranking, flat: str) -> None:
    """A gate is a gate.

    The counterpart of `test_the_highest_scoring_agent_is_on_this_page` in
    `test_drugs_to_avoid.py`, asserted from the accepting side: a pipeline where the
    top score always wins cannot express "right about the biology, not yet a
    candidate", and the page has to show that rather than claim it.
    """
    accepted = ranking.accepted[0].row
    top = ranking.top_scoring
    assert top.rejected
    assert top.score > accepted.score
    assert top.directions_agree is True
    assert accepted.directions_agree is None

    assert f"The highest score of the run is {top.score:.3f}" in flat
    assert f"{top.name} — is **rejected**" in flat


def test_the_accepted_candidate_carries_its_concerns_onto_the_page(
    probe: AcceptedProbe, flat: str
) -> None:
    """Accepted means not disqualified, which is a weaker claim than it looks.

    An accepted hypothesis must leave `rejection_reasons` empty, so its near misses
    live in `concerns`. A page that prints the rank and drops the concerns presents the
    candidate as cleaner than the pipeline found it.
    """
    assert probe.concerns, "an accepted candidate with no concerns needs a different test"
    for concern in probe.concerns:
        assert f"`{concern.value}`" in flat
        assert CONCERN_CLAUSE[concern].split(".")[0] in flat
    assert "Accepted means not disqualified" in flat


def test_every_accepted_agent_has_a_specific_probe_rationale(ranking: Ranking) -> None:
    """No generic fallback sentence is allowed to reach the page."""
    for probe in ranking.accepted:
        key = (probe.row.node_id, probe.row.observed)
        assert key in PROBE_RATIONALE_BY_NODE_AND_DIRECTION
        assert len(probe.rationale) > 120, (
            f"{probe.row.name}'s rationale is too short to say anything specific about "
            f"this node and this direction: {probe.rationale!r}"
        )


def test_an_uncovered_node_and_direction_fails_loudly() -> None:
    """The canary for the rationale table.

    A new catalogue row that survives the gates on an uncovered node must stop the
    build, not quietly render a page with a blank where the experiment's purpose
    should be.
    """
    from tools.track2.ranking import _rationale

    with pytest.raises(ProbeRationaleMissingError, match="No probe rationale for"):
        _rationale("N_NOT_A_REAL_NODE", "increase", "Fabricatinib")


def test_an_unexplained_concern_fails_loudly() -> None:
    """The canary for the concern table.

    A concern the pipeline learns to raise must not reach the page as a bare enum
    value, and must not be dropped either. Both of those make the candidate look
    better than it is, so the generator refuses instead.
    """
    from tools.track2.ranking import _concern_clause

    assert RejectionReason.SAFETY_CONCERN not in CONCERN_CLAUSE
    with pytest.raises(ConcernClauseMissingError, match="has no clause in CONCERN_CLAUSE"):
        _concern_clause(RejectionReason.SAFETY_CONCERN, "Fabricatinib")


def test_a_synthetic_catalogue_can_never_be_described_as_literature_curated(
    ranking: Ranking, flat: str
) -> None:
    """GP-20, tested on both branches rather than trusted on the live one.

    The live run reads `knowledge/literature/bub1b/`, so the happy path only ever
    exercises the literature wording. A hardcoded provenance line would keep asserting
    biological validity if this generator were pointed at the demo tables, which is
    precisely the sentence GP-20 forbids -- so the wording is read off the loader's own
    flag, and this pins both outcomes.
    """
    from tools.track2.ranking import _provenance_sentence

    assert ranking.catalogue_is_synthetic is False
    assert "literature-curated" in flat
    assert "SYNTHETIC" not in flat

    pretend_synthetic = Ranking(
        mechanism_id=ranking.mechanism_id,
        gene_symbol=ranking.gene_symbol,
        catalogue_version="synthetic-v0.0",
        catalogue_is_synthetic=True,
        n_evidence=ranking.n_evidence,
        rows=(),
        accepted=(),
    )
    sentence = _provenance_sentence(pretend_synthetic)
    assert "SYNTHETIC" in sentence
    assert "fictional agents" in sentence
    assert "literature-curated" not in sentence
    assert "GP-20" in sentence


def test_the_published_weights_are_the_pipelines_weights(flat: str) -> None:
    """The score frame cannot drift away from the code that computes the score.

    The page explains where 0.500 comes from. If a weight moved in
    `src/mva/interventions/generate.py` and the explanation did not, the page would be
    describing a scoring function that no longer exists -- a quieter failure than a
    wrong score, because the number still looks right.
    """
    for weight in (DIRECTION_WEIGHT, EVIDENCE_WEIGHT, SAFETY_WEIGHT, APPROVAL_WEIGHT):
        assert f"| {weight:.2f} |" in flat
    assert f"collects {DIRECTION_UNKNOWN_CREDIT:.2f} of that term" in flat
    assert "DIRECTION_UNKNOWN_CREDIT" in flat


def test_the_page_says_it_is_generated_and_says_by_what(flat: str) -> None:
    assert "GENERATED FILE — do not edit" in flat
    assert "tools/track2/ranking.py" in flat
    assert "just ranking" in flat
    assert "tests/unit/test_ranking.py" in flat


#: A number immediately followed by a dose or concentration unit. A range such as
#: "750-1,000 mg" is caught on its second number, so no dash needs to be matched.
_DOSE = re.compile(r"\b\d[\d.,]*\s*(mg|mcg|[\u00b5\u03bc]g|ug|ml|mL|IU|[\u00b5\u03bc]M|uM|nM|mM)\b")


def test_the_page_refuses_to_be_read_as_medical_advice(flat: str) -> None:
    """ "Accepted" is a statement about this pipeline, not about a patient.

    The dose scan is not decoration. This page renders the catalogue's own
    `proposed_validation_experiment` text verbatim, so a dose added to a knowledge
    table upstream would arrive on a document that promises there is no dose on it,
    and nobody would have made a decision about that.
    """
    assert "Nothing here is medical advice" in flat
    assert "there is no dose on this page" in flat
    assert "experimental probe" in flat
    assert "has never been measured in any genotype" in flat
    assert "research hypothesis requiring pre-clinical validation" in flat
    assert "Not a recommendation and not a dose" in flat

    found = _DOSE.findall(flat)
    assert not found, (
        f"a dose-like quantity reached {OUTPUT_PATH}: {found}. The page states that no "
        "dose appears on it, so either the claim or the content is now wrong. "
        "Remediation: find which upstream field introduced it (most likely a "
        "`knowledge/literature/bub1b/` cell rendered verbatim) and decide deliberately "
        "whether this document should carry it. Do not delete the promise to make this "
        "pass."
    )
