"""Generate ``docs/track2-ranking.md`` from the pipeline's own accepted ranking.

Why this is generated and not written
-------------------------------------
``tools/track2/drugs_to_avoid.py`` already renders the twelve rejections from a live
triage run, so a reader can reproduce *why an agent was rejected* by command. The
other half of the same headline could only be taken on trust: the accepted agent and
its score of 0.500 appeared in `docs/track2-hypothesis.md` §5 and nowhere else, typed
by hand, with nothing in the repository able to notice when the catalogue, the
mechanism or a scoring weight moved underneath them. The ranking is the claim a
reader is most likely to quote and was the one claim no command could reproduce.

So this page is rendered from a live
:func:`~mva.interventions.generate.generate_drug_hypotheses` run and
``tests/unit/test_ranking.py`` fails if the committed file and the generator disagree.
It reads the same tables, in the same order, under the same fixed clock as its
sibling page, and it imports that page's source constants rather than restating them:
two generated pages describing two different runs while both claim to describe the
pipeline would be worse than one hand-written page.

The one thing this page exists to get right
-------------------------------------------
A score is not a recommendation, and rank 1 is the sentence most easily misread in
this repository. Three separations hold the page up:

* **accepted is not recommended.** It means the agent cleared every gate this
  pipeline applies, on public literature about a gene, with no patient data. The
  accepted agent is an experimental probe. There is no dose on this page.
* **accepted is not clean.** An accepted hypothesis must leave ``rejection_reasons``
  empty, so its near misses live in ``concerns``, and this page prints them. Without
  somewhere to put them the triage stage computes the near misses and then discards
  them, and a candidate reads as unblemished when it is merely not disqualified.
* **the top score does not win.** Approval status and an unassessed
  chromosomal-instability question are gates rather than score terms. The
  highest-scoring agent of the run is on the sibling page, rejected, and it is the
  only one of the thirteen whose direction at its node was measured and agreed.

And one rule the renderer enforces mechanically: the accepted score is never printed
without the unmeasured direction beside it. The direction of effect on BubR1
abundance has never been measured in any genotype, which is the whole question the
first experiment asks — a score quoted on its own silently answers it.
"""

from __future__ import annotations

import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from mva.clock import demo_clock
from mva.interventions import (
    APPROVAL_WEIGHT,
    DIRECTION_UNKNOWN_CREDIT,
    DIRECTION_WEIGHT,
    EVIDENCE_WEIGHT,
    SAFETY_WEIGHT,
    DrugCatalog,
    DrugTriageResult,
    generate_drug_hypotheses,
)
from mva.mechanisms import MechanismLibrary
from mva.models.drug import DrugHypothesis, RejectionReason
from mva.models.mechanism import MechanismHypothesis

# The tables, their version label and the sibling page's name come from the sibling
# generator rather than being retyped here. Both pages must describe ONE run: if these
# two constants ever drifted apart, the repository would present an acceptance from run
# A beside rejections from run B, with both pages claiming to be rendered from the
# pipeline and nothing in the test suite able to tell.
from tools.track2.drugs_to_avoid import LIBRARY_VERSION, LITERATURE_DIR
from tools.track2.drugs_to_avoid import OUTPUT_PATH as AVOID_PATH

__all__ = [
    "CONCERN_CLAUSE",
    "LIBRARY_VERSION",
    "LITERATURE_DIR",
    "OUTPUT_PATH",
    "PROBE_RATIONALE_BY_NODE_AND_DIRECTION",
    "UNMEASURED_DIRECTION_MARKERS",
    "AcceptedProbe",
    "ConcernClauseMissingError",
    "ProbeRationaleMissingError",
    "RankRow",
    "Ranking",
    "build_ranking",
    "render",
    "run_triage",
    "write",
]

#: Relative to the repository root.
OUTPUT_PATH: Final[str] = "docs/track2-ranking.md"

GENERATOR: Final[str] = "tools/track2/ranking.py"
DRIFT_TEST: Final[str] = "tests/unit/test_ranking.py"
RECIPE: Final[str] = "just ranking"

#: Phrases that count as stating the accepted agent's direction is unmeasured.
#:
#: The renderer may not print the accepted score in any block that contains none of
#: them, and ``tests/unit/test_ranking.py`` checks the rendered page block by block.
#: Every phrase here is about the *direction*, deliberately: "never measured" on its
#: own would also match the chromosomal-instability cell, and a block could then pass
#: the check while saying nothing at all about the sign. The rule is cheap and the
#: failure it prevents is the one that matters — "nicotinic acid, 0.500, rank 1" is a
#: sentence that reads as a finding about a drug and is a finding about a pipeline.
UNMEASURED_DIRECTION_MARKERS: Final[tuple[str, ...]] = (
    "cannot determine (unsigned)",
    "unsigned direction",
    "direction unmeasured",
)

#: What a first experiment at this node, pushed this way, can and cannot establish.
#:
#: Keyed on the mechanism rather than on the drug, and on the same ``(node,
#: direction)`` key space as ``DANGER_BY_NODE_AND_DIRECTION`` in the sibling
#: generator — the two tables sit on opposite sides of the same gate and a pair may
#: legitimately appear in both, because what makes an agent interesting at a node and
#: what makes it dangerous there are frequently the same fact with the sign flipped.
#: Only pairs a run has actually *accepted* are listed. A new catalogue row that
#: survives the gates on an uncovered pair is a loud failure rather than a bland
#: sentence (see :func:`_rationale`).
PROBE_RATIONALE_BY_NODE_AND_DIRECTION: Final[dict[tuple[str, str], str]] = {
    ("N3_bubr1", "unknown"): (
        "The agent acts on the therapeutic target and nobody has measured which way it "
        "moves it. That is what makes it a probe rather than a candidate, and it fixes "
        "what the first experiment is for: measuring the sign, not trying a treatment. "
        "The readout that decides the hypothesis is whether BubR1 abundance moves at "
        "all, and in which direction, in cells carrying the patient's genotype — an "
        "experiment built that way is informative whichever way it comes out, and one "
        "that assumed the sign would be informative only if it guessed right. The "
        "agent is here because the question is answerable, not because the answer is "
        "known."
    ),
}

#: Why a recorded, non-fatal concern is still worth a reader's attention.
#:
#: ``concerns`` is the field an accepted hypothesis puts its near misses in, so these
#: clauses are the difference between "accepted" and "accepted, and here is what is
#: still wrong with it". Keyed on the pipeline's own enum, so a concern the pipeline
#: learns to raise cannot reach the page without someone writing down what it means
#: (see :func:`_concern_clause`).
CONCERN_CLAUSE: Final[dict[RejectionReason, str]] = {
    RejectionReason.DIRECTION_UNKNOWN: (
        "The direction this agent pushes the target has never been measured. Recorded "
        "as a flag and a score penalty, never as agreement and never as a "
        "wrong-direction rejection: `directions_agree` is tri-state and `None` is "
        "neither end of it (GP-16, ASSUMPTION-DRUG-02). This is the concern the first "
        "experiment exists to retire."
    ),
    RejectionReason.PHARMACOKINETIC_BARRIER: (
        "Exposure relative to the effective concentration cannot be computed, because "
        "one of the two figures has never been recorded. `concentration_achievable` is "
        "therefore `None` and stays `None`: the gap is kept as a concern rather than "
        "resolved optimistically (ASSUMPTION-DRUG-05). An agent that works in culture "
        "at a concentration the body never reaches is not a therapy, and nothing on "
        "this page can tell you which case this is."
    ),
}

#: Read against an *accepted* agent, the three chromosomal-instability states carry a
#: different emphasis from the same three on the rejection page, which is why the
#: wording is written out again rather than shared. There the question is what a gap
#: does not license; here it is what a closed gate does not promise.
_CIN_CLAUSE: Final[dict[str, str]] = {
    "measured_positive": (
        "The chromosomal-instability question **was asked of this agent and the answer "
        "was yes**. This is a measurement, not a gap."
    ),
    "never_measured": (
        "The chromosomal-instability question has **never been asked of this agent**. "
        "GP-14: that is the absence of a finding and not a finding, and in a "
        "chromosomal-instability disorder it blocks (ADR 0011, ASSUMPTION-DRUG-07)."
    ),
    "measured_negative": (
        "The chromosomal-instability question was asked and came back negative — in "
        "adults, a population that is not this one. That closes the gate and is not a "
        "clean bill of health: a negative measured elsewhere is not a measurement here "
        "(ADR 0011, ASSUMPTION-DRUG-06). Micronucleus frequency for this agent on a "
        "BUB1B-deficient background remains unmeasured, and the first experiment reads "
        "it out rather than deferring it to a later study."
    ),
}


@dataclass(frozen=True, slots=True)
class RankRow:
    """One screened agent as the ranked table prints it — accepted or rejected."""

    rank: int | None
    drug_id: str
    name: str
    node_id: str
    required: str
    observed: str
    directions_agree: bool | None
    worsens_cin: bool | None
    score: float
    rejected: bool
    primary_reason: str
    n_concerns: int

    @property
    def direction_verdict(self) -> str:
        if self.directions_agree is False:
            return "OPPOSES (measured)"
        if self.directions_agree is True:
            return "agrees (measured)"
        return "cannot determine (unsigned)"

    @property
    def cin_state(self) -> str:
        if self.worsens_cin is True:
            return "measured_positive"
        if self.worsens_cin is False:
            return "measured_negative"
        return "never_measured"

    @property
    def cin_cell(self) -> str:
        return {
            "measured_positive": "**measured — yes**",
            "measured_negative": "measured — no (other population)",
            "never_measured": "**never measured**",
        }[self.cin_state]

    @property
    def rank_cell(self) -> str:
        return "—" if self.rank is None else f"**{self.rank}**"

    @property
    def outcome(self) -> str:
        """What the pipeline did with it, and at what cost in the accepted case."""
        if self.rejected:
            return f"rejected — `{self.primary_reason}`"
        if self.n_concerns:
            return f"**accepted**, {self.n_concerns} recorded concerns"
        return "**accepted**"


@dataclass(frozen=True, slots=True)
class SafetySignal:
    """A non-disqualifying safety signal the accepted agent still carries."""

    concern_id: str
    severity: str
    population: str
    description: str


@dataclass(frozen=True, slots=True)
class AcceptedProbe:
    """An accepted agent, reduced to what this page has to say about it.

    Holds its :class:`RankRow` rather than copying its fields, so the table and the
    detail block below it can never disagree about a score or a verdict.
    """

    row: RankRow
    node_label: str
    approved_name: str | None
    approval_status: str
    intervention_class: str
    target: str
    evidence_type: str
    is_direct_evidence: bool
    concerns: tuple[RejectionReason, ...]
    safety_signals: tuple[SafetySignal, ...]
    validation_experiment: str
    rationale: str

    @property
    def evidence_sentence(self) -> str:
        """Tier and directness, both printed, with the right gloss for each case.

        The two are separate fields and the gloss must follow the second one, not the
        first. An earlier draft hardcoded the indirect wording onto whatever tier came
        back, which on a direct-evidence agent would have printed "no measurement of
        this agent on this target" directly underneath the measurement.
        """
        if self.is_direct_evidence:
            return (
                f"strongest tier `{self.evidence_type}`, and it is **direct** — the "
                f"evidence is for this agent on this target in a relevant system. "
                f"Directness is not strength: it says the measurement was made here, "
                f"not that it was made well or in this population."
            )
        return (
            f"strongest tier `{self.evidence_type}`, and it is **indirect** — no "
            f"measurement of this agent on this target in a relevant system. The tier "
            f"and the directness are separate fields and both are printed, because the "
            f"tier alone reads as strength it does not have here: a human trial is the "
            f"highest tier this catalogue records, and this one was not run on this "
            f"target, in this disease, or in this population."
        )


@dataclass(frozen=True, slots=True)
class Ranking:
    """One triage pass, reduced to what the page needs. Nothing is recomputed here."""

    mechanism_id: str
    gene_symbol: str
    catalogue_version: str
    catalogue_is_synthetic: bool
    n_evidence: int
    rows: tuple[RankRow, ...]
    accepted: tuple[AcceptedProbe, ...]

    @property
    def rejected_rows(self) -> tuple[RankRow, ...]:
        return tuple(row for row in self.rows if row.rejected)

    @property
    def top_scoring(self) -> RankRow:
        """Highest score of the whole run, accepted or not. Frequently not accepted."""
        return max(self.rows, key=lambda row: row.score)


class ProbeRationaleMissingError(RuntimeError):
    """No rationale covers an accepted agent's (node, direction) pair.

    Raised rather than papered over with a generic sentence. A page that quietly
    degrades to "this agent was accepted" has stopped saying what the experiment is
    for, which on the accepting side is the whole content, and it looks exactly like
    a page that is working.
    """


class ConcernClauseMissingError(RuntimeError):
    """An accepted agent carries a recorded concern this page cannot explain.

    The concern list is the only thing separating "accepted" from "clean". Printing a
    bare enum value a reader cannot interpret spends the reader's trust without
    informing them, and dropping it silently would be worse.
    """


def _rationale(node_id: str, observed: str, name: str) -> str:
    try:
        return PROBE_RATIONALE_BY_NODE_AND_DIRECTION[(node_id, observed)]
    except KeyError:
        msg = (
            f"No probe rationale for {name} acting on {node_id!r} in the {observed!r} "
            f"direction. Add the pair to PROBE_RATIONALE_BY_NODE_AND_DIRECTION in "
            f"{GENERATOR}, stating what a first experiment on THAT node in THAT "
            f"direction can and cannot establish in this disease. Do not write a "
            f"generic sentence, and do not describe an unmeasured direction as a "
            f"therapeutic rationale (GP-16): the experiment's purpose is the only "
            f"thing on this page a reader cannot get from the catalogue."
        )
        raise ProbeRationaleMissingError(msg) from None


def _concern_clause(reason: RejectionReason, name: str) -> str:
    try:
        return CONCERN_CLAUSE[reason]
    except KeyError:
        msg = (
            f"{name} was accepted carrying the concern {reason.value!r}, which has no "
            f"clause in CONCERN_CLAUSE in {GENERATOR}. Add one saying what the concern "
            f"means for a reader and what would retire it. A concern is what stands "
            f"between 'accepted' and 'clean'; printing the bare enum value, or leaving "
            f"it off the page, both present the candidate as better than the pipeline "
            f"found it."
        )
        raise ConcernClauseMissingError(msg) from None


def run_triage(repo_root: Path) -> tuple[MechanismHypothesis, DrugCatalog, DrugTriageResult]:
    """One deterministic triage pass over the literature tables (GP-30).

    Returns the whole result, not a slice of it: this page's subject is the ranking,
    so it needs the accepted side, the rejected side and the catalogue's own
    provenance flag in one pass. A second triage run to fetch the rest would be a
    second chance to disagree.
    """
    literature = repo_root / LITERATURE_DIR
    mechanism = MechanismLibrary.from_tsv(
        literature / "mechanisms.tsv",
        literature / "mechanism_meta.tsv",
        version=LIBRARY_VERSION,
    ).for_gene("BUB1B")
    if mechanism is None:
        msg = f"{literature}/mechanisms.tsv carries no chain for BUB1B."
        raise RuntimeError(msg)
    # synthetic=False is mandatory here: these rows are literature-curated (ADR 0022)
    # and the default is True, so omitting it labels every evidence row this triage
    # produces as describing fictional agents.
    catalog = DrugCatalog.from_tsv(
        literature / "drug_catalog.tsv", version=LIBRARY_VERSION, synthetic=False
    )
    result = generate_drug_hypotheses(mechanism=mechanism, catalog=catalog, clock=demo_clock())
    return mechanism, catalog, result


def _row(hypothesis: DrugHypothesis) -> RankRow:
    return RankRow(
        rank=hypothesis.rank,
        drug_id=hypothesis.drug_id,
        name=hypothesis.name,
        node_id=hypothesis.target_node_id,
        required=hypothesis.required_direction.value,
        observed=hypothesis.observed_direction.value,
        directions_agree=hypothesis.directions_agree,
        worsens_cin=hypothesis.worsens_chromosomal_instability,
        score=hypothesis.score,
        rejected=hypothesis.rejected,
        primary_reason=(
            hypothesis.rejection_reasons[0].value if hypothesis.rejection_reasons else ""
        ),
        n_concerns=len(hypothesis.concerns),
    )


def _probe(hypothesis: DrugHypothesis, *, node_label: str) -> AcceptedProbe:
    return AcceptedProbe(
        row=_row(hypothesis),
        node_label=node_label,
        approved_name=hypothesis.approved_name,
        approval_status=hypothesis.approval_status.value,
        intervention_class=hypothesis.intervention_class.value,
        target=hypothesis.target,
        evidence_type=hypothesis.strongest_evidence_type.value,
        is_direct_evidence=hypothesis.is_direct_evidence,
        concerns=hypothesis.concerns,
        safety_signals=tuple(
            SafetySignal(
                concern_id=concern.concern_id,
                severity=concern.severity,
                population=concern.population,
                description=concern.description,
            )
            for concern in hypothesis.safety_concerns
        ),
        validation_experiment=hypothesis.proposed_validation_experiment,
        rationale=_rationale(
            hypothesis.target_node_id,
            hypothesis.observed_direction.value,
            hypothesis.name,
        ),
    )


def build_ranking(repo_root: Path) -> Ranking:
    """Every screened agent in the pipeline's own total order, accepted first."""
    mechanism, catalog, result = run_triage(repo_root)
    labels = {node.node_id: node.label for node in mechanism.nodes}
    if not result.accepted:
        msg = (
            "The triage run accepted nothing, so there is no ranking to render and "
            f"{OUTPUT_PATH} would be a page asserting a result that does not exist. "
            "This is the failure `just science-canary` exists to catch: a stage that "
            "rejects everything passes every rejection test in the suite. Fix the "
            "stage or the tables; do not render an empty page."
        )
        raise RuntimeError(msg)
    accepted = tuple(
        _probe(
            hypothesis,
            node_label=labels.get(hypothesis.target_node_id, "(node not on this chain)"),
        )
        for hypothesis in result.accepted
    )
    return Ranking(
        mechanism_id=mechanism.mechanism_id,
        gene_symbol=mechanism.gene_symbol,
        catalogue_version=catalog.version,
        catalogue_is_synthetic=catalog.synthetic,
        n_evidence=len(result.evidence),
        rows=tuple(_row(hypothesis) for hypothesis in result.all_hypotheses),
        accepted=accepted,
    )


#: Prose wrap width. Interpolated counts change line lengths, so paragraphs are
#: wrapped at render time rather than laid out by hand -- otherwise a page that is
#: regenerated after a catalogue change comes back with ragged, unreviewable diffs.
#: Matches the sibling page, so the two read as one document set.
_WRAP: Final[int] = 88


def _para(text: str) -> str:
    return _fill(text)


def _bullet(text: str) -> str:
    return _fill(text, initial_indent="- ", subsequent_indent="  ")


def _fill(text: str, **kwargs: str) -> str:
    """Wrap prose without ever splitting on a hyphen.

    ``break_on_hyphens`` is off because the identifiers this page cites are
    hyphenated -- ``ASSUMPTION-DRUG-05``, ``GP-16``, ``N3_bubr1`` -- and a wrap inside
    one of them silently breaks the docs-integrity lint that checks every cited
    identifier resolves. A citation split across a newline reads as authority and
    points at nothing.
    """
    return textwrap.fill(" ".join(text.split()), width=_WRAP, break_on_hyphens=False, **kwargs)


def _provenance_sentence(ranking: Ranking) -> str:
    """GP-20, read off the catalogue rather than asserted.

    A hardcoded "literature-curated" line would keep claiming biological validity if
    this generator were ever pointed at the demo tables, which is the exact sentence
    GP-20 exists to forbid. So the flag the loader set decides the wording.
    """
    if ranking.catalogue_is_synthetic:
        return (
            "**The catalogue loaded for this render is SYNTHETIC and describes "
            "fictional agents.** Nothing below is pharmacological fact, it carries no "
            "clinical validity, and the ranking is a demonstration that the machinery "
            "runs — not a result about any real drug (GP-20)."
        )
    return (
        "The catalogue is **literature-curated**: every row cites a specific "
        "publication, listed in `knowledge/literature/bub1b/SOURCES.md`. That is real "
        "pharmacology and it is not a pharmacological database — it is one author's "
        "reading of the literature on one date, with no completeness guarantee, so the "
        "absence of an agent is not evidence against it (ADR 0022, GP-20)."
    )


def _probe_block(probe: AcceptedProbe) -> str:
    row = probe.row
    approved = f"`{probe.approved_name}`" if probe.approved_name else "none (not approved)"
    lines = [
        f"### Rank {row.rank} — {row.name} — `{row.drug_id}`",
        "",
        _bullet(
            f"**Regulatory standing:** approval status `{probe.approval_status}`, approved "
            f"name {approved}, intervention class `{probe.intervention_class}`. Repurposing "
            f"means redeploying an agent that already exists; this field is what makes the "
            f"word applicable, and it says nothing about whether the agent works here."
        ),
        _bullet(f"**Node acted on:** `{row.node_id}` — {probe.node_label}"),
        _bullet(f"**Molecular target:** {probe.target}"),
        _bullet(
            f"**Required correction:** `{row.required}` — **observed:** `{row.observed}` — "
            f"**verdict:** {row.direction_verdict}"
        ),
        _bullet(
            f"**Score:** {row.score:.3f} of 1, carrying an unsigned direction at the node "
            f"it acts on. The score ranks hypotheses against each other inside one run. It "
            f"is not a probability, not a strength of evidence and not a recommendation."
        ),
        _bullet(
            f"**Chromosomal-instability question:** {row.cin_cell}. "
            f"{_CIN_CLAUSE[row.cin_state]}"
        ),
        _bullet(f"**Evidence:** {probe.evidence_sentence}"),
    ]
    if probe.concerns:
        lines.append(
            _bullet(
                f"**Recorded concerns ({len(probe.concerns)}, none fatal):** "
                f"{', '.join(f'`{c.value}`' for c in probe.concerns)}. An accepted "
                f"hypothesis must leave `rejection_reasons` empty, so these live in a "
                f"field of their own instead of being discarded — \"accepted despite "
                f"this\" is the sentence a reviewer needs and the one a bare rank hides."
            )
        )
        for concern in probe.concerns:
            lines.append(
                _fill(
                    f"`{concern.value}` — {_concern_clause(concern, probe.row.name)}",
                    initial_indent="  - ",
                    subsequent_indent="    ",
                )
            )
    else:
        lines.append(
            _bullet(
                "**Recorded concerns:** none. That means no gate raised a non-fatal flag, "
                "not that no concern exists: the pipeline records only the questions it "
                "knows how to ask (GP-14)."
            )
        )
    lines.extend(
        [
            _bullet(
                "**Safety signals carried (none disqualifying):** "
                + (
                    "; ".join(
                        f"`{s.concern_id}` ({s.severity}, {s.population})"
                        for s in probe.safety_signals
                    )
                    if probe.safety_signals
                    else "none recorded"
                )
            ),
        ]
    )
    lines.extend(
        _fill(
            f"`{signal.concern_id}` — {signal.description}",
            initial_indent="  - ",
            subsequent_indent="    ",
        )
        for signal in probe.safety_signals
    )
    lines.extend(
        [
            _bullet(
                "**What a first experiment here can and cannot establish:** "
                f"{probe.rationale}"
            ),
            _bullet(
                "**The experiment that comes first, in the pipeline's own words:** "
                f"{probe.validation_experiment}"
            ),
            "",
        ]
    )
    return "\n".join(lines)


def render(repo_root: Path) -> str:
    """The whole page. Pure: same inputs, same bytes (GP-30)."""
    ranking = build_ranking(repo_root)
    top = ranking.top_scoring
    accepted_rows = tuple(row for row in ranking.rows if not row.rejected)
    agreeing = [row for row in ranking.rows if row.directions_agree is True]
    unsigned = [row for row in ranking.rows if row.directions_agree is None]
    opposing = [row for row in ranking.rows if row.directions_agree is False]
    never_measured_cin = [
        row for row in ranking.rejected_rows if row.cin_state == "never_measured"
    ]

    lines: list[str] = [
        "# Track 2 — the accepted drug ranking in BUB1B-related mosaic variegated aneuploidy",
        "",
        f"<!-- GENERATED by {GENERATOR}. Do not edit this file. -->",
        "",
        _para(
            f"**GENERATED FILE — do not edit.** Rendered from a live triage run by "
            f"`{GENERATOR}`; regenerate with `{RECIPE}`. `{DRIFT_TEST}` fails if this file "
            f"and the generator disagree, so the page cannot drift away from the pipeline "
            f"it describes."
        ),
        "",
        _para(
            "**Nothing here is medical advice, and there is no dose on this page.** "
            '"Accepted" means *this agent cleared every gate this pipeline applies*, on '
            "public literature about a gene, with no patient data of any kind. It does not "
            "mean recommended, it does not mean safe and it does not mean a treatment "
            "exists. The accepted agent is an **experimental probe**: its direction of "
            "effect on BubR1 abundance has never been measured in any genotype, so the "
            "first experiment measures that sign instead of testing a treatment. Every "
            "row below is a research hypothesis requiring pre-clinical validation."
        ),
        "",
        _para(
            f"Source: `{LITERATURE_DIR}/` — mechanism `{ranking.mechanism_id}`, gene "
            f"{ranking.gene_symbol}, catalogue version `{ranking.catalogue_version}`."
        ),
        "",
        _para(_provenance_sentence(ranking)),
        "",
        "---",
        "",
        "## What the run produced",
        "",
        _para(
            f"{len(ranking.rows)} agents screened, {len(accepted_rows)} accepted, "
            f"{len(ranking.rejected_rows)} rejected, {ranking.n_evidence} evidence rows, "
            f"every rejection kept with its reasons (GP-19). The rejections are rendered "
            f"in full on `{AVOID_PATH}` from this same pass over these same tables under "
            f"the same fixed clock; this page is the accepting side of it."
        ),
        "",
        "## The accepted agent is not the highest-scoring one",
        "",
        _para(
            f"Rank 1 scores {accepted_rows[0].score:.3f} with an unsigned direction at the "
            f"therapeutic target. The highest score of the run is {top.score:.3f} and that "
            f"agent — {top.name} — is **rejected**. It is also the only one of the "
            f"{len(ranking.rows)} whose direction at its node was measured and found to "
            f"agree."
        ),
        "",
        _para(
            "Because approval status and an unassessed chromosomal-instability question "
            "are **gates** rather than score terms, an agent can be right about the "
            "biology and still not be a candidate (ADR 0011, ASSUMPTION-DRUG-07). A "
            "pipeline in which the top score always wins has no way to say that sentence, "
            "and this run is the case that shows the gate firing: the agent it stops is "
            "the one the ranking likes most."
        ),
        "",
        _para(
            f"The direction column is tri-state throughout and is never rounded (GP-16, "
            f"ASSUMPTION-DRUG-02). Of {len(ranking.rows)} agents, {len(agreeing)} was "
            f"measured to move its node the right way, {len(opposing)} were measured to "
            f"move it the wrong way, and {len(unsigned)} carry a direction nobody has "
            f"measured. The third group is neither of the first two."
        ),
        "",
        _para(
            f"**Do not read the rejected rows as a list of harmful agents.** "
            f"{len(never_measured_cin)} of the {len(ranking.rejected_rows)} rejections have "
            f"a blank chromosomal-instability cell: nobody ever asked whether they increase "
            f"aneuploidy. GP-14 — that is the absence of a finding, not a finding of harm. "
            f"In this disease an unasked question is disqualifying by design, which is a "
            f"fact about the gate and not about the agent. `{AVOID_PATH}` separates the two "
            f"cases agent by agent."
        ),
        "",
        "---",
        "",
        "## Ranked order",
        "",
        "| Rank | Agent | Node | Required | Observed | Direction verdict | Worsens CIN? | "
        "Score | Outcome |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    lines.extend(
        f"| {row.rank_cell} | {row.name} | `{row.node_id}` | `{row.required}` | "
        f"`{row.observed}` | "
        f"{'**OPPOSES**' if row.directions_agree is False else row.direction_verdict} | "
        f"{row.cin_cell} | {row.score:.3f} | {row.outcome} |"
        for row in ranking.rows
    )
    lines.extend(
        [
            "",
            _para(
                "Rejected rows are ordered by score, which is the order the pipeline "
                "produced and not a ranking of how unsuitable they are: a gate does not "
                "score, and two agents stopped by the same gate can sit at opposite ends "
                "of the table."
            ),
            "",
            "---",
            "",
            f"## The accepted candidate ({len(ranking.accepted)})",
            "",
            _para(
                "Every mandatory question this pipeline asks is answered by a field below, "
                'including the ones whose answer is "nobody has measured it". An '
                "incomplete hypothesis cannot be constructed; a hypothesis that is "
                "complete and still unsigned at its target node can, and this is one."
            ),
            "",
        ]
    )
    lines.extend(_probe_block(probe) for probe in ranking.accepted)
    lines.extend(
        [
            "---",
            "",
            "## Where the score comes from",
            "",
            _para(
                "The four weights below are module constants in "
                "`src/mva/interventions/generate.py`, read by this page rather than retyped "
                "into it. They are documented heuristics, not calibrated parameters: "
                "changing one requires a decision record, a test and a before/after "
                "comparison (GP-32)."
            ),
            "",
            "| Term | Weight | What it asks |",
            "|---|---|---|",
            f"| Direction | {DIRECTION_WEIGHT:.2f} | Does the agent push the node the way "
            f"the mechanism requires? |",
            f"| Evidence | {EVIDENCE_WEIGHT:.2f} | How far is the evidence from a "
            f"demonstrated human outcome? |",
            f"| Safety | {SAFETY_WEIGHT:.2f} | What did the safety pass find, weighted by "
            f"severity? |",
            f"| Approval | {APPROVAL_WEIGHT:.2f} | Is the agent actually available to "
            f"redeploy? |",
            "",
            _para(
                f"Direction carries the largest single weight, and an agent whose direction "
                f"is unmeasured collects {DIRECTION_UNKNOWN_CREDIT:.2f} of that term — "
                f"`DIRECTION_UNKNOWN_CREDIT`, the fraction retained when agreement is "
                f"undeterminable. Well below 1.0 because an unknown is never an agreement, "
                f"and well above 0.0 because an unknown is not a disagreement either "
                f"(GP-14, GP-16). That is the term the first experiment would move."
            ),
            "",
            "---",
            "",
            '## What "accepted" does not mean',
            "",
            _bullet(
                "**Not a recommendation and not a dose.** This pipeline models no exposure "
                "a patient would receive. It reads public literature about a gene and ranks "
                "research hypotheses; prescribing is outside everything it does."
            ),
            _bullet(
                "**Not a measured direction.** `unknown` is not agreement. It is a flag and "
                "a score penalty, and the hypothesis survives carrying it (GP-16, "
                "ASSUMPTION-DRUG-02). Reading rank 1 as \"this raises BubR1\" states a "
                "result nobody has produced."
            ),
            _bullet(
                "**Not a clean bill of health.** The chromosomal-instability gate is closed "
                "on a negative measured in another population, which is enough to clear a "
                "gate and not enough to support a claim (ASSUMPTION-DRUG-06)."
            ),
            _bullet(
                "**Not unblemished.** The recorded concerns above are the pipeline's own "
                "near misses, kept rather than discarded. Accepted means not disqualified."
            ),
            _bullet(
                "**One context only.** The ranking scores the soma. An agent's direction can "
                "be desirable in the tumour and harmful in every other tissue of the same "
                "child, so a rejection here is not a verdict on tumour-directed use and this "
                "acceptance is not a verdict on it either (ASSUMPTION-DRUG-08)."
            ),
            "",
            "---",
            "",
            "## Provenance",
            "",
            _bullet(f"Generator: `{GENERATOR}`"),
            _bullet(f"Drift test: `{DRIFT_TEST}` — regenerate with `{RECIPE}`"),
            _bullet(
                f"Rendered from `{LITERATURE_DIR}/mechanisms.tsv`, "
                f"`{LITERATURE_DIR}/mechanism_meta.tsv` and "
                f"`{LITERATURE_DIR}/drug_catalog.tsv`"
            ),
            _bullet(
                f"Citations for every row: `{LITERATURE_DIR}/SOURCES.md`, machine-verified "
                f"against PubMed by `tests/unit/test_citations.py`"
            ),
            _bullet(
                "Clock: `mva.clock.demo_clock()` — fixed, so repeat renders are "
                "byte-identical (GP-30)"
            ),
            _bullet(
                f"Scores, ranks, concerns and rejection reasons are "
                f"`{DrugHypothesis.__module__}.DrugHypothesis` fields written by the "
                f"triage stage, never by this generator"
            ),
            _bullet(f"The twelve rejections, in full: `{AVOID_PATH}`"),
            "",
        ]
    )
    return "\n".join(lines).rstrip("\n") + "\n"


def write(repo_root: Path) -> Path:
    """Render and write the page. Returns the path written."""
    path = repo_root / OUTPUT_PATH
    path.write_text(render(repo_root), encoding="utf-8")
    return path
