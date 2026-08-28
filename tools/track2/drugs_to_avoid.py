"""Generate ``docs/track2-drugs-to-avoid.md`` from the pipeline's own rejections.

Why this is generated and not written
-------------------------------------
``docs/track2-hypothesis.md`` §5 lists the rejected agents inline. That list was
correct when it was typed. The moment the catalogue, the mechanism or a scoring
weight changes, it is a claim about a pipeline that no longer exists, and nothing
in the repository would notice. So this page is rendered from a live
:func:`~mva.interventions.generate.generate_drug_hypotheses` run and
``tests/unit/test_drugs_to_avoid.py`` fails if the committed file and the generator
disagree.

The one thing this page exists to get right
-------------------------------------------
Two sentences read almost identically in a rejection list and mean opposite things:

* **OPPOSES** -- somebody measured what the agent does at this node, and it does the
  reverse of what the disease needs. That is a *finding*.
* **never measured** -- nobody has asked. That is the *absence* of a finding.

GP-14 says absence of information is not negative information, and the direction of
the mistake matters in both directions. Reading "never measured" as a safety
finding invents evidence that does not exist. Reading it as reassurance is worse:
in a chromosomal-instability syndrome the unasked question is disqualifying by
design (ADR 0011), and an agent sitting in that column may be harmless or may be
the most dangerous on the page. The record cannot tell you which. That is the whole
point of keeping the two columns apart.

There is a third state the pipeline's single ``oncogenic_risk`` label collapses,
and this page pulls it back out: for two agents the chromosomal-instability
question *was* asked and the answer was yes. Same rejection reason, opposite
epistemic status.

Nothing here is medical advice. "Avoid" means "rejected by this pipeline for
systemic use in the affected child, on public literature about a gene" -- never
"useless in this disease", and several of these agents are retained in the record
precisely because they are plausible tumour-directed candidates under a separate
risk assessment (ASSUMPTION-DRUG-08).
"""

from __future__ import annotations

import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from mva.clock import demo_clock
from mva.interventions import DrugCatalog, generate_drug_hypotheses
from mva.mechanisms import MechanismLibrary
from mva.models.drug import DrugHypothesis, RejectionReason
from mva.models.mechanism import MechanismHypothesis

__all__ = [
    "DANGER_BY_NODE_AND_DIRECTION",
    "LITERATURE_DIR",
    "OUTPUT_PATH",
    "AvoidEntry",
    "build_entries",
    "render",
    "run_triage",
]

#: Relative to the repository root.
OUTPUT_PATH: Final[str] = "docs/track2-drugs-to-avoid.md"
LITERATURE_DIR: Final[str] = "knowledge/literature/bub1b"

#: Matches ``tests/unit/test_track2_bub1b.py``. Distinct from ``synthetic-v0.0`` so a
#: reader can tell literature-derived evidence from the fictional demo tables (GP-20).
LIBRARY_VERSION: Final[str] = "bub1b-literature-v1"

GENERATOR: Final[str] = "tools/track2/drugs_to_avoid.py"

#: Why an agent is dangerous *in this disease*, keyed by the node it acts on and the
#: direction it was observed to push that node.
#:
#: Keyed on the mechanism rather than on the drug, deliberately. A per-agent table
#: would be the hand-written list this file exists to replace, and would say nothing
#: about the next agent that lands on the same node. Keyed this way, a new catalogue
#: row acting on ``N5_sac`` in the ``decrease`` direction inherits the correct
#: warning the day it is added, and an agent that lands on a node with no rule is a
#: loud failure rather than a bland sentence (see :func:`_danger`).
DANGER_BY_NODE_AND_DIRECTION: Final[dict[tuple[str, str], str]] = {
    ("N5_sac", "decrease"): (
        "It lowers spindle assembly checkpoint signalling — the exact signalling this "
        "patient's BubR1 deficiency already cannot sustain. This is not a side effect "
        "that resembles the disease; it is the disease mechanism applied "
        "pharmacologically. Maximal target proximity, inverted sign: any ranking by "
        "pathway membership or target closeness that cannot read signs puts an agent "
        "like this first."
    ),
    ("N9_stress", "increase"): (
        "It raises aneuploidy-associated proteotoxic and metabolic stress. That is "
        "selectively lethal to aneuploid cells, which is exactly why the class is "
        "attractive against the rhabdomyosarcoma and exactly why it is dangerous "
        "everywhere else: in this syndrome every tissue is mosaically aneuploid and "
        "the child is already growth-restricted. Rejected for the soma on the same "
        "evidence that recommends it for the tumour."
    ),
    ("N9_stress", "unknown"): (
        "Its net effect on the stress node in an aneuploid background has never been "
        "measured, and the two plausible signs point opposite ways: lowering "
        "translational load could relieve the burden, or removing the survival "
        "programme could kill the cells carrying it. Until that sign exists there is "
        "nothing here to weigh."
    ),
    ("N10_clearance", "decrease"): (
        "It removes p53- and p16-dependent arrest, senescence and immune clearance of "
        "missegregating cells. In BubR1-hypomorphic mice that same response drives the "
        "progeroid phenotype *and* suppresses aneuploidy-driven tumorigenesis, so "
        "removing it in a cancer-predisposed child may remove a barrier rather than a "
        "burden. The node is flagged compensatory: no corrective direction is derived "
        "from its state, and the direction check returns cannot-determine rather than "
        "guessing (ASSUMPTION-MECHANISM-04)."
    ),
    ("N3_bubr1", "increase"): (
        "The direction is right — this agent pushes the therapeutic target the way the "
        "mechanism needs. It is on this page for the questions nobody has answered "
        "about it, not for its biology, and it is the clearest demonstration that the "
        "highest-scoring candidate is not automatically the recommended one."
    ),
    ("N3_bubr1", "unknown"): (
        "It acts on the therapeutic target and nobody has measured which way it moves "
        "it. An unmeasured sign at the target node is not a detail outstanding: it is "
        "the entire question the hypothesis asks."
    ),
    ("N3_bubr1", "context_dependent"): (
        "It moves the therapeutic target in both directions at once. On this axis that "
        "is the NAD+ precursor that is also a sirtuin product inhibitor: it feeds the "
        "pool SIRT2 consumes and supplies the molecule that blocks SIRT2, in the same "
        "cell at the same time, and which term dominates at a given intracellular "
        "concentration has never been measured in any genotype. A context-dependent "
        "sign is worse than a missing one for a first experiment, because it destroys "
        "interpretability in both outcomes: a rise in BubR1 cannot distinguish NAD+ "
        "supply from an off-axis effect, and a fall cannot distinguish a wrong "
        "mechanism from having inhibited the enzyme the mechanism runs on. The only "
        "direct measurement of this agent class on this endpoint runs the wrong way — "
        "the sirtuin-inhibiting concentration lowered BubR1 protein — at an exposure "
        "well above anything a patient reaches, which fixes the sign of the inhibitory "
        "arm without fixing the dose at which it wins."
    ),
}

_CIN_CLAUSE: Final[dict[str | None, str]] = {
    "measured_positive": (
        "The chromosomal-instability question **was asked of this agent and the answer "
        "was yes**. This is a measurement, not a gap."
    ),
    "never_measured": (
        "The chromosomal-instability question has **never been asked of this agent**. "
        "GP-14: that is not a clean bill of health and must never be filled in "
        "optimistically. In a chromosomal-instability disorder it blocks (ADR 0011, "
        "ASSUMPTION-DRUG-07)."
    ),
    "measured_negative": (
        "The chromosomal-instability question was asked and came back negative — in a "
        "population that is not this one. Enough to clear a gate, not enough to "
        "support a claim (ASSUMPTION-DRUG-06)."
    ),
}


@dataclass(frozen=True, slots=True)
class AvoidEntry:
    """One rejected agent, reduced to what this page has to say about it."""

    drug_id: str
    name: str
    node_id: str
    node_label: str
    required: str
    observed: str
    directions_agree: bool | None
    worsens_cin: bool | None
    primary_reason: str
    all_reasons: tuple[str, ...]
    score: float
    danger: str

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
    def opposes(self) -> bool:
        """Measured, and measured to point the wrong way."""
        return self.directions_agree is False


class DangerRuleMissingError(RuntimeError):
    """No warning rule covers a rejected agent's (node, direction) pair.

    Raised rather than papered over with a generic sentence. A page that quietly
    degrades to "this agent was rejected" for a new agent is a page that stops
    warning anybody, and it looks exactly like a page that is working.
    """


def _danger(node_id: str, observed: str, name: str) -> str:
    try:
        return DANGER_BY_NODE_AND_DIRECTION[(node_id, observed)]
    except KeyError:
        msg = (
            f"No danger rule for {name} acting on {node_id!r} in the {observed!r} "
            f"direction. Add the pair to DANGER_BY_NODE_AND_DIRECTION in "
            f"{GENERATOR}, stating why an agent that moves THAT node THAT way is "
            f"dangerous in this disease specifically. Do not write a generic "
            f"sentence: the per-disease reasoning is the only thing on this page a "
            f"reader cannot get from the catalogue."
        )
        raise DangerRuleMissingError(msg) from None


def run_triage(repo_root: Path) -> tuple[MechanismHypothesis, tuple[DrugHypothesis, ...]]:
    """One deterministic triage pass over the literature tables (GP-30)."""
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
    return mechanism, result.rejected


def build_entries(repo_root: Path) -> tuple[MechanismHypothesis, tuple[AvoidEntry, ...]]:
    """Every rejected agent, in the pipeline's own order (worst-scoring last)."""
    mechanism, rejected = run_triage(repo_root)
    labels = {node.node_id: node.label for node in mechanism.nodes}
    entries = tuple(
        AvoidEntry(
            drug_id=hypothesis.drug_id,
            name=hypothesis.name,
            node_id=hypothesis.target_node_id,
            node_label=labels.get(hypothesis.target_node_id, "(node not on this chain)"),
            required=hypothesis.required_direction.value,
            observed=hypothesis.observed_direction.value,
            directions_agree=hypothesis.directions_agree,
            worsens_cin=hypothesis.worsens_chromosomal_instability,
            primary_reason=(
                hypothesis.rejection_reasons[0].value if hypothesis.rejection_reasons else ""
            ),
            all_reasons=tuple(reason.value for reason in hypothesis.rejection_reasons),
            score=hypothesis.score,
            danger=_danger(
                hypothesis.target_node_id,
                hypothesis.observed_direction.value,
                hypothesis.name,
            ),
        )
        for hypothesis in rejected
    )
    return mechanism, entries


#: Prose wrap width. Interpolated counts change line lengths, so paragraphs are
#: wrapped at render time rather than laid out by hand -- otherwise a page that is
#: regenerated after a catalogue change comes back with ragged, unreviewable diffs.
_WRAP: Final[int] = 88


def _para(text: str) -> str:
    return _fill(text)


def _bullet(text: str) -> str:
    return _fill(text, initial_indent="- ", subsequent_indent="  ")


def _fill(text: str, **kwargs: str) -> str:
    """Wrap prose without ever splitting on a hyphen.

    ``break_on_hyphens`` is off because the identifiers this page cites are
    hyphenated -- ``ASSUMPTION-DRUG-08``, ``GP-14``, ``N3_bubr1`` -- and a wrap
    inside one of them silently breaks the docs-integrity lint that checks every
    cited identifier resolves. A citation split across a newline reads as authority
    and points at nothing.
    """
    return textwrap.fill(" ".join(text.split()), width=_WRAP, break_on_hyphens=False, **kwargs)


def _names(entries: list[AvoidEntry]) -> str:
    """Name the agents rather than hardcoding one, so the prose follows the run."""
    names = [entry.name for entry in entries]
    if not names:
        return "No agent"
    if len(names) == 1:
        return names[0]
    return ", ".join(names[:-1]) + f" and {names[-1]}"


def _detail_block(entry: AvoidEntry) -> str:
    return "\n".join(
        [
            f"#### {entry.name} — `{entry.drug_id}`",
            "",
            _bullet(f"**Node acted on:** `{entry.node_id}` — {entry.node_label}"),
            _bullet(
                f"**Required correction:** `{entry.required}` — **observed:** "
                f"`{entry.observed}` — **verdict:** {entry.direction_verdict}"
            ),
            _bullet(
                f"**Chromosomal-instability question:** {entry.cin_cell}. "
                f"{_CIN_CLAUSE[entry.cin_state]}"
            ),
            _bullet(
                f"**Rejected for:** `{entry.primary_reason}` (all reasons: "
                f"{', '.join(f'`{r}`' for r in entry.all_reasons)})"
            ),
            _bullet(f"**Why it is dangerous here specifically:** {entry.danger}"),
            "",
        ]
    )


def render(repo_root: Path) -> str:
    """The whole page. Pure: same inputs, same bytes (GP-30)."""
    mechanism, entries = build_entries(repo_root)
    opposes = [e for e in entries if e.opposes]
    unsigned = [e for e in entries if not e.opposes]
    measured_cin = [e for e in entries if e.cin_state == "measured_positive"]
    never_measured_cin = [e for e in entries if e.cin_state == "never_measured"]
    agreeing = [e for e in unsigned if e.directions_agree is True]

    lines: list[str] = [
        "# Track 2 — agents to avoid in BUB1B-related mosaic variegated aneuploidy",
        "",
        f"<!-- GENERATED by {GENERATOR}. Do not edit this file. -->",
        "",
        _para(
            f"**GENERATED FILE — do not edit.** Rendered from a live triage run by "
            f"`{GENERATOR}`; regenerate with `just drugs-to-avoid`. "
            f"`tests/unit/test_drugs_to_avoid.py` fails if this file and the generator "
            f"disagree, so the page cannot drift away from the pipeline it describes."
        ),
        "",
        _para(
            '**Nothing here is medical advice.** "Avoid" means *rejected by this '
            "pipeline for systemic use in the affected child*, derived from public "
            "literature about a gene, with no patient data of any kind. It never means "
            '"useless in this disease": several agents below are retained in the record '
            "precisely because they are plausible **tumour**-directed candidates under "
            "an entirely separate risk assessment (ASSUMPTION-DRUG-08). Every row is a "
            "research hypothesis requiring pre-clinical validation."
        ),
        "",
        _para(
            f"Source: `{LITERATURE_DIR}/` — mechanism `{mechanism.mechanism_id}`, gene "
            f"{mechanism.gene_symbol}, catalogue version `{LIBRARY_VERSION}`."
        ),
        "",
        "---",
        "",
        "## Read this first: two rejections that look the same and are not",
        "",
        _para(
            "The pipeline rejected every agent below. It did **not** reject them for "
            "the same kind of reason, and collapsing the difference is the fastest way "
            "to turn this page into misinformation."
        ),
        "",
        _para(
            f"**OPPOSES — measured, and measured the wrong way.** {len(opposes)} of "
            f"{len(entries)} agents. Somebody ran the experiment at this node and the "
            f"agent moves it away from what the disease needs. That is a *finding*, and "
            f"it is disqualifying on the biology."
        ),
        "",
        _para(
            f"**Never measured.** {len(never_measured_cin)} of {len(entries)} agents "
            f"have a blank chromosomal-instability cell, and several also carry an "
            f"unsigned direction. Nobody has asked. **GP-14: absence of information is "
            f"not negative information.** Two mistakes are available here and both are "
            f"bad:"
        ),
        "",
        _bullet(
            "Reading it as a safety finding invents evidence that does not exist. "
            "Nobody has shown these agents do anything harmful at this node."
        ),
        _bullet(
            "Reading it as reassurance is worse. In a chromosomal-instability syndrome "
            "the unasked question is disqualifying **by design** (ADR 0011, "
            "ASSUMPTION-DRUG-07). An agent in this column may be harmless or may be the "
            "most dangerous on the page, and the record cannot tell you which."
        ),
        "",
        _para(
            f"**Measured, and the answer was yes.** {len(measured_cin)} agents. The "
            f"pipeline raises the same `oncogenic_risk` label for a measured positive as "
            f"it does for a blank cell, so this page separates them: the label is "
            f"identical and the epistemic status is opposite."
        ),
        "",
        _para(
            "The practical consequence is §5 of `docs/track2-hypothesis.md`: the "
            "bottleneck in repurposing for chromosomal-instability disorders is not "
            "finding targets. It is that the safety question that matters most has not "
            "been asked of the candidate drugs."
        ),
        "",
        "---",
        "",
        "## Summary",
        "",
        "| Agent | Node | Required | Observed | Direction verdict | Worsens CIN? | "
        "Primary reason | Score |",
        "|---|---|---|---|---|---|---|---|",
    ]
    lines.extend(
        f"| {e.name} | `{e.node_id}` | `{e.required}` | `{e.observed}` | "
        f"{'**OPPOSES**' if e.opposes else e.direction_verdict} | {e.cin_cell} | "
        f"`{e.primary_reason}` | {e.score:.3f} |"
        for e in entries
    )
    lines.extend(
        [
            "",
            "---",
            "",
            f"## A. Measured to point the wrong way ({len(opposes)})",
            "",
            _para(
                "For each of these the experiment was done. The direction is a "
                "measurement, and it disagrees with the correction this mechanism "
                "requires."
            ),
            "",
        ]
    )
    lines.extend(_detail_block(entry) for entry in opposes)
    lines.extend(
        [
            "---",
            "",
            f"## B. Direction not disqualifying — rejected on an unanswered "
            f"question ({len(unsigned)})",
            "",
            _para(
                f"Not one of these was measured to move its node the wrong way. "
                f"{_names(agreeing)} was measured to move it the **right** way, at the "
                f"therapeutic target, and is rejected anyway. The remaining "
                f"{len(unsigned) - len(agreeing)} carry an unsigned direction, which is "
                f"tri-state and is **never** rounded to agreement (GP-16): it is a "
                f"recorded concern and a score penalty, not a verdict."
            ),
            "",
            _para(
                f"What ends all {len(unsigned)} is the same thing — a question with no "
                f"answer rather than an answer nobody likes. That is a gate, not a score "
                f"term, which is why the highest-scoring agent of the whole run "
                f"({_names(agreeing)}, {max(e.score for e in unsigned):.3f}) is on this "
                f"page and the accepted candidate scores lower."
            ),
            "",
        ]
    )
    lines.extend(_detail_block(entry) for entry in unsigned)
    lines.extend(
        [
            "---",
            "",
            "## Provenance",
            "",
            _bullet(f"Generator: `{GENERATOR}`"),
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
                f"Rejection reasons are `{RejectionReason.__module__}.RejectionReason` "
                f"values written by the pipeline, not by this generator"
            ),
            "",
        ]
    )
    return "\n".join(lines).rstrip("\n") + "\n"


def write(repo_root: Path) -> Path:
    """Render and write the page. Returns the path written."""
    path = repo_root / OUTPUT_PATH
    path.write_text(render(repo_root), encoding="utf-8")
    return path
