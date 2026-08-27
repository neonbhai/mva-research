"""Construct the mechanism hypothesis for a candidate pair, with its evidence.

This stage answers one question — *by what chain of steps would this gene produce
this phenotype?* — and answers it as a structure rather than a paragraph, so that
the intervention stage can interrogate it mechanically: which node do we push,
which way, and which links are we only guessing at.

Two commitments shape the code:

* **Nothing is asserted without an `EvidenceItem` (GP-10)**, including the
  relevance score itself. A number that influences a ranking is a claim.
* **The weak links are advertised, not buried (GP-17).** `mechanism_relevance_score`
  penalises inferred links explicitly, and every inferred link gets its own
  evidence row saying, in words, that it was not demonstrated.

The builder imports no peer stage (GP-03): it receives a `MechanismLibrary` and a
`Clock`, and its output is passed *into* the intervention stage by the pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass

from mva.clock import Clock
from mva.mechanisms.library import MechanismLibrary
from mva.models.base import AssertionTier
from mva.models.evidence import (
    Citation,
    EvidenceCategory,
    EvidenceDirection,
    EvidenceItem,
    EvidenceStrength,
    EvidenceType,
    make_evidence_id,
)
from mva.models.mechanism import (
    DiscriminatingExperiment,
    MechanismHypothesis,
    MechanismLink,
)

__all__ = [
    "DEMONSTRATED_FRACTION_WEIGHT",
    "INFERRED_LINK_PENALTY",
    "LINK_STRENGTH_WEIGHT",
    "MEAN_STRENGTH_WEIGHT",
    "MechanismResult",
    "build_mechanism",
    "mechanism_relevance_score",
]

TOOL_NAME = "mva.mechanisms.builder"
TOOL_VERSION = "0.1.0"

#: GP-20. Appended to the limitations of every evidence item this module emits.
SYNTHETIC_CHAIN_LIMITATION = (
    "The mechanism chain behind this claim is a SYNTHETIC curated table shipped for "
    "a fictional demo case. It is not derived from the literature, not biologically "
    "valid, and must never be read as a real mechanistic finding."
)

#: Ordinal strength mapped to a weight for scoring only. GP-32: these are tunable
#: heuristics, not measurements; changing one requires a decision record and a
#: before/after comparison against the golden expectations.
LINK_STRENGTH_WEIGHT: dict[EvidenceStrength, float] = {
    EvidenceStrength.DEFINITIVE: 1.00,
    EvidenceStrength.STRONG: 0.85,
    EvidenceStrength.MODERATE: 0.60,
    EvidenceStrength.SUPPORTING: 0.45,
    EvidenceStrength.WEAK: 0.25,
    EvidenceStrength.INSUFFICIENT: 0.05,
}

#: Weight on the mean per-link strength of the chain.
MEAN_STRENGTH_WEIGHT = 0.55
#: Weight on the fraction of links that were directly demonstrated.
DEMONSTRATED_FRACTION_WEIGHT = 0.45
#: Flat deduction per inferred link, on top of the fraction term.
INFERRED_LINK_PENALTY = 0.05


@dataclass(frozen=True)
class MechanismResult:
    """What the mechanism stage hands to the composition root.

    `hypothesis` is ``None`` when no chain is curated for the gene. That is a
    reportable absence, not a failure: it comes with an evidence row saying so
    (GP-14), so a downstream reader can distinguish "no mechanism known" from
    "mechanism not looked for".
    """

    hypothesis: MechanismHypothesis | None
    evidence: tuple[EvidenceItem, ...]
    warnings: tuple[str, ...]


def _clamp01(value: float) -> float:
    return max(0.0, min(1.0, value))


def mechanism_relevance_score(hypothesis: MechanismHypothesis | None) -> float:
    """How much weight the chain can bear, in ``[0, 1]``.

    Formula (all three terms are module constants; GP-32)::

        mean_strength         = mean(LINK_STRENGTH_WEIGHT[link.strength] for link in links)
        demonstrated_fraction = (n_links - n_inferred) / n_links
        score = clamp01(
              MEAN_STRENGTH_WEIGHT         * mean_strength          # 0.55
            + DEMONSTRATED_FRACTION_WEIGHT * demonstrated_fraction  # 0.45
            - INFERRED_LINK_PENALTY        * n_inferred             # 0.05 each
        )

    Inferred links are penalised **twice**, deliberately. The fraction term asks
    "what proportion of this chain is demonstrated?", which a long chain can
    dilute; the flat per-link deduction asks "how many places could this chain
    break?", which it cannot. A chain is only as good as its weakest joint, and
    `MechanismHypothesis.inferred_links` names those joints exactly.

    Returns ``0.0`` for a missing hypothesis and for a chain with no links: an
    unlinked node set makes no mechanistic claim and must not score as if it did.
    """
    if hypothesis is None or not hypothesis.links:
        return 0.0
    links = hypothesis.links
    mean_strength = sum(LINK_STRENGTH_WEIGHT[link.strength] for link in links) / len(links)
    n_inferred = len(hypothesis.inferred_links)
    demonstrated_fraction = (len(links) - n_inferred) / len(links)
    raw = (
        MEAN_STRENGTH_WEIGHT * mean_strength
        + DEMONSTRATED_FRACTION_WEIGHT * demonstrated_fraction
        - INFERRED_LINK_PENALTY * n_inferred
    )
    return round(_clamp01(raw), 6)


def build_mechanism(
    gene_symbol: str,
    *,
    pair_id: str | None,
    library: MechanismLibrary,
    clock: Clock,
) -> MechanismResult:
    """Bind the curated chain for `gene_symbol` to a candidate pair.

    The chain itself is knowledge, not a finding about this patient; binding it to
    a `pair_id` is what turns it into a hypothesis *about this case*. The evidence
    rows emitted here therefore describe the chain and its gaps — never the
    patient's genotype, which lives in the ingestion and prioritisation stages.
    """
    hypothesis = library.for_gene(gene_symbol)
    citation_version = library.version
    if hypothesis is None:
        absence = _evidence(
            subject_id=gene_symbol,
            subject_kind="gene",
            claim=(
                f"No curated mechanism chain is available for {gene_symbol} in mechanism "
                f"library {citation_version}."
            ),
            category=EvidenceCategory.MECHANISM,
            direction=EvidenceDirection.NEUTRAL,
            strength=EvidenceStrength.INSUFFICIENT,
            evidence_type=EvidenceType.PIPELINE_INFERENCE,
            tier=AssertionTier.INFERENCE,
            method="Lookup of the gene symbol in the curated mechanism library.",
            limitations=(
                "Absence from this library is absence of curation, NOT evidence that no "
                "mechanism exists (GP-14). The library covers a single synthetic demo gene."
            ),
            clock=clock,
            citation=None,
        )
        warning = (
            f"No mechanism chain for gene {gene_symbol}; downstream drug triage cannot run "
            "for this candidate without one."
        )
        return MechanismResult(hypothesis=None, evidence=(absence,), warnings=(warning,))

    citation = Citation(
        source="mva-knowledge/mechanisms.tsv",
        identifier=hypothesis.mechanism_id,
        version=citation_version,
        title="Synthetic curated mechanism chain",
    )
    items: list[EvidenceItem] = [_chain_evidence(hypothesis, clock=clock, citation=citation)]
    items.extend(
        _link_evidence(hypothesis, link, clock=clock, citation=citation)
        for link in hypothesis.links
    )
    score = mechanism_relevance_score(hypothesis)
    items.append(_score_evidence(hypothesis, score=score, clock=clock))

    bound = hypothesis.model_copy(
        update={
            "pair_id": pair_id,
            "supporting_evidence_ids": tuple(item.evidence_id for item in items),
            "discriminating_experiments": _discriminating_experiments(hypothesis),
        }
    )
    return MechanismResult(
        hypothesis=bound,
        evidence=tuple(items),
        warnings=_warnings(bound, score=score),
    )


# --------------------------------------------------------------------- internals


def _evidence(
    *,
    subject_id: str,
    subject_kind: str,
    claim: str,
    category: EvidenceCategory,
    direction: EvidenceDirection,
    strength: EvidenceStrength,
    evidence_type: EvidenceType,
    tier: AssertionTier,
    method: str,
    limitations: str,
    clock: Clock,
    citation: Citation | None,
    numeric_value: float | None = None,
    payload: dict[str, str | int | float | bool | None] | None = None,
) -> EvidenceItem:
    """Build one evidence row with a content-derived, reproducible ID (GP-30)."""
    return EvidenceItem(
        evidence_id=make_evidence_id(
            subject_id=subject_id, category=category, claim=claim, tool=TOOL_NAME
        ),
        subject_id=subject_id,
        subject_kind=subject_kind,
        claim=claim,
        category=category,
        direction=direction,
        strength=strength,
        evidence_type=evidence_type,
        tier=tier,
        citation=citation,
        method=method,
        tool=TOOL_NAME,
        tool_version=TOOL_VERSION,
        limitations=limitations,
        timestamp=clock.now(),
        numeric_value=numeric_value,
        payload=payload or {},
    )


def _weakest(links: tuple[MechanismLink, ...]) -> EvidenceStrength:
    """The chain's weakest link strength — the honest strength of the whole chain."""
    return min(links, key=lambda link: LINK_STRENGTH_WEIGHT[link.strength]).strength


def _chain_evidence(
    hypothesis: MechanismHypothesis, *, clock: Clock, citation: Citation
) -> EvidenceItem:
    terminal = hypothesis.nodes[-1]
    n_inferred = len(hypothesis.inferred_links)
    claim = (
        f"A {len(hypothesis.nodes)}-node, {len(hypothesis.links)}-link curated chain connects "
        f"{hypothesis.gene_symbol} to {terminal.label}, with the therapeutic target at "
        f"{hypothesis.therapeutic_target_node_id} requiring "
        f"{hypothesis.required_correction.value}."
    )
    limitations = (
        f"{n_inferred} of {len(hypothesis.links)} links are inferred rather than directly "
        "demonstrated, so the chain is not a proven causal path. The chain is curated "
        "knowledge about the gene, not a measurement in this individual. "
        f"{SYNTHETIC_CHAIN_LIMITATION}"
    )
    if hypothesis.developmental_window_caveat:
        limitations += " Developmental-window caveat: " + hypothesis.developmental_window_caveat
    return _evidence(
        subject_id=hypothesis.mechanism_id,
        subject_kind="mechanism",
        claim=claim,
        category=EvidenceCategory.MECHANISM,
        direction=EvidenceDirection.SUPPORTS,
        strength=_weakest(hypothesis.links) if hypothesis.links else EvidenceStrength.INSUFFICIENT,
        evidence_type=EvidenceType.CURATED_DATABASE,
        tier=AssertionTier.DATABASE_ASSERTION,
        method="Assembled from the curated mechanism chain and metadata tables.",
        limitations=limitations,
        clock=clock,
        citation=citation,
        payload={
            "gene_symbol": hypothesis.gene_symbol,
            "n_nodes": len(hypothesis.nodes),
            "n_links": len(hypothesis.links),
            "n_inferred_links": n_inferred,
            "therapeutic_target_node_id": hypothesis.therapeutic_target_node_id,
            "required_correction": hypothesis.required_correction.value,
        },
    )


def _link_evidence(
    hypothesis: MechanismHypothesis,
    link: MechanismLink,
    *,
    clock: Clock,
    citation: Citation,
) -> EvidenceItem:
    """One row per link, so a reader can audit the chain joint by joint."""
    claim = (
        f"{link.link_id}: {link.source_node_id} {link.relation} {link.target_node_id} "
        f"({link.direction.value})."
    )
    if link.is_directly_demonstrated:
        limitations = f"{link.uncertainty} {SYNTHETIC_CHAIN_LIMITATION}"
        evidence_type = EvidenceType.CURATED_DATABASE
        direction = EvidenceDirection.SUPPORTS
    else:
        limitations = (
            "INFERRED link: no experiment demonstrated this step in a relevant system; it "
            f"is carried by analogy or pathway membership. {link.uncertainty} "
            f"{SYNTHETIC_CHAIN_LIMITATION}"
        )
        evidence_type = EvidenceType.PIPELINE_INFERENCE
        direction = EvidenceDirection.NEUTRAL
    return _evidence(
        subject_id=hypothesis.mechanism_id,
        subject_kind="mechanism",
        claim=claim,
        category=EvidenceCategory.MECHANISM,
        direction=direction,
        strength=link.strength,
        evidence_type=evidence_type,
        tier=link.tier,
        method="Read from the curated mechanism chain table; one row per link.",
        limitations=limitations,
        clock=clock,
        citation=citation,
        payload={
            "link_id": link.link_id,
            "source_node_id": link.source_node_id,
            "target_node_id": link.target_node_id,
            "direction": link.direction.value,
            "is_directly_demonstrated": link.is_directly_demonstrated,
        },
    )


def _score_evidence(hypothesis: MechanismHypothesis, *, score: float, clock: Clock) -> EvidenceItem:
    """GP-10: the relevance score influences a ranking, so it is itself a cited claim."""
    n_inferred = len(hypothesis.inferred_links)
    claim = (
        f"Mechanism relevance for {hypothesis.mechanism_id} scores {score:.3f} from "
        f"{len(hypothesis.links)} links, {n_inferred} of them inferred."
    )
    return _evidence(
        subject_id=hypothesis.mechanism_id,
        subject_kind="mechanism",
        claim=claim,
        category=EvidenceCategory.MECHANISM,
        direction=EvidenceDirection.NEUTRAL,
        strength=EvidenceStrength.SUPPORTING,
        evidence_type=EvidenceType.PIPELINE_INFERENCE,
        tier=AssertionTier.INFERENCE,
        method=(
            f"score = clamp01({MEAN_STRENGTH_WEIGHT} * mean_link_strength + "
            f"{DEMONSTRATED_FRACTION_WEIGHT} * demonstrated_fraction - "
            f"{INFERRED_LINK_PENALTY} * n_inferred_links)"
        ),
        limitations=(
            "A heuristic weighting of curated strength labels, not a probability. It "
            "measures how well-evidenced the chain is, NOT whether the chain is true, and "
            "not whether correcting the target would help the patient."
        ),
        clock=clock,
        citation=None,
        numeric_value=score,
        payload={"n_links": len(hypothesis.links), "n_inferred_links": n_inferred},
    )


def _discriminating_experiments(
    hypothesis: MechanismHypothesis,
) -> tuple[DiscriminatingExperiment, ...]:
    """Propose one falsification test per inferred link.

    The inferred links are precisely where the chain could be wrong, so they are
    where an experiment buys the most. Each proposal names the alternative it would
    rule out; an experiment that cannot fail against a stated alternative is a
    demonstration, not a test.
    """
    experiments: list[DiscriminatingExperiment] = []
    labels = {node.node_id: node.label for node in hypothesis.nodes}
    for link in hypothesis.inferred_links:
        source = labels.get(link.source_node_id, link.source_node_id)
        target = labels.get(link.target_node_id, link.target_node_id)
        experiments.append(
            DiscriminatingExperiment(
                experiment_id=f"DX-{hypothesis.mechanism_id}-{link.link_id}",
                description=(
                    f"Quantify {target} as a function of {source} in a relevant model system, "
                    f"testing the inferred link {link.link_id}."
                ),
                measures=f"Magnitude of {target} across a graded range of {source}.",
                distinguishes_from=(
                    f"The alternative that {target} arises independently of {source}, and that "
                    f"{link.link_id} is pathway co-membership rather than causation."
                ),
                expected_if_true=f"{target} tracks {source} monotonically.",
                expected_if_false=f"{target} is unchanged across the range of {source}.",
                feasibility="specialised",
            )
        )
    return tuple(experiments)


def _warnings(hypothesis: MechanismHypothesis, *, score: float) -> tuple[str, ...]:
    warnings: list[str] = []
    inferred = hypothesis.inferred_links
    if inferred:
        named = ", ".join(link.link_id for link in inferred)
        warnings.append(
            f"{len(inferred)} of {len(hypothesis.links)} links in {hypothesis.mechanism_id} are "
            f"INFERRED, not directly demonstrated ({named}); the chain is not a proven causal path."
        )
    if hypothesis.developmental_window_caveat:
        warnings.append(
            f"{hypothesis.mechanism_id} carries a developmental-window caveat: a mechanistically "
            "correct agent may still be therapeutically irrelevant post-natally."
        )
    if score < 0.5:
        warnings.append(
            f"Mechanism relevance for {hypothesis.mechanism_id} is low ({score:.3f}); treat any "
            "drug hypothesis built on it as exploratory."
        )
    return tuple(warnings)
