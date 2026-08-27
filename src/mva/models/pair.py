"""Candidate variant-pair models.

A ``CandidatePair`` is the Track 1 unit of prediction. Its defining design choice
is that the score is a **vector, not a scalar**. Component scores stay visible all
the way into the report, because "why is this ranked first?" must be answerable
without re-running the pipeline.
"""

from __future__ import annotations

import hashlib
from enum import StrEnum
from typing import Self

from pydantic import Field, computed_field, model_validator

from mva.models.base import FrozenModel
from mva.models.variant import VariantRecord


class InheritanceModel(StrEnum):
    """Candidate mode of inheritance under consideration."""

    COMPOUND_HETEROZYGOUS = "compound_heterozygous"
    HOMOZYGOUS_RECESSIVE = "homozygous_recessive"
    DE_NOVO_DOMINANT = "de_novo_dominant"
    AUTOSOMAL_DOMINANT = "autosomal_dominant"
    X_LINKED_RECESSIVE = "x_linked_recessive"
    X_LINKED_DOMINANT = "x_linked_dominant"
    MITOCHONDRIAL = "mitochondrial"
    MOSAIC = "mosaic"
    UNKNOWN = "unknown"


class PhaseStatus(StrEnum):
    """Whether the two variants are on opposite haplotypes.

    ``UNKNOWN`` is the honest default for a proband-only short-read VCF and must
    be preserved, not optimistically upgraded. Two heterozygous variants in the
    same gene are only a compound heterozygote **if they are in trans**; assuming
    trans without evidence is a headline scientific error this codebase refuses to
    make. See docs/scientific-assumptions.md (ASSUMPTION-PHASE-01).
    """

    TRANS_CONFIRMED = "trans_confirmed"
    """Proven opposite haplotypes: parental data, read-backed phasing, or long reads."""

    TRANS_LIKELY = "trans_likely"
    """Statistically supported but not proven (e.g. population phasing)."""

    UNKNOWN = "unknown"
    """No phasing information. The expected state for proband-only exome/genome."""

    CIS_LIKELY = "cis_likely"
    CIS_CONFIRMED = "cis_confirmed"
    """Same haplotype. Cannot constitute a compound heterozygote — one allele is
    still intact — so this is a near-disqualifying finding for a recessive model."""


#: Phase states incompatible with a compound-heterozygous recessive mechanism.
CIS_STATES: frozenset[PhaseStatus] = frozenset({PhaseStatus.CIS_LIKELY, PhaseStatus.CIS_CONFIRMED})


class PhaseEvidence(FrozenModel):
    """How a phase determination was reached (or why it could not be)."""

    status: PhaseStatus
    method: str = Field(
        description="'parental_trio' | 'read_backed' | 'long_read' | 'phase_set' | 'none'"
    )
    supporting_reads: int | None = Field(default=None, ge=0)
    distance_bp: int | None = Field(
        default=None,
        ge=0,
        description=(
            "Distance between the two sites. Read-backed phasing is only possible "
            "within roughly a fragment length; a large distance explains why phase "
            "is UNKNOWN rather than implying anything about it."
        ),
    )
    notes: str | None = None


class ComponentScores(FrozenModel):
    """The score vector. Every component is on [0, 1] and independently inspectable.

    These are NOT calibrated probabilities. They are transparent heuristics whose
    weights live in configuration; see docs/decisions/0005-separate-filter-and-rank.md.
    """

    analytical_validity: float = Field(ge=0.0, le=1.0)
    """Do we believe the calls are real? Depth, allele balance, FILTER, mappability."""

    rarity: float = Field(ge=0.0, le=1.0)
    """How consistent the allele frequencies are with a severe recessive disorder."""

    molecular_consequence: float = Field(ge=0.0, le=1.0)
    """Predicted damage to the gene product."""

    inheritance_consistency: float = Field(ge=0.0, le=1.0)
    """Fit to the inheritance model, including the phase penalty."""

    phenotype_similarity: float = Field(ge=0.0, le=1.0)
    """Overlap between gene-associated phenotype and the proband's HPO profile."""

    mechanistic_relevance: float = Field(ge=0.0, le=1.0)
    """Whether the gene's known biology plausibly produces this phenotype."""

    evidence_quality: float = Field(ge=0.0, le=1.0)
    """Strength and independence of the supporting evidence, not its quantity."""

    contradiction_penalty: float = Field(
        ge=0.0,
        le=1.0,
        description="Subtracted, not multiplied. 0 means nothing contradicts the hypothesis.",
    )

    def as_dict(self) -> dict[str, float]:
        return {
            "analytical_validity": self.analytical_validity,
            "rarity": self.rarity,
            "molecular_consequence": self.molecular_consequence,
            "inheritance_consistency": self.inheritance_consistency,
            "phenotype_similarity": self.phenotype_similarity,
            "mechanistic_relevance": self.mechanistic_relevance,
            "evidence_quality": self.evidence_quality,
            "contradiction_penalty": self.contradiction_penalty,
        }


class OpenQuestion(FrozenModel):
    """A named gap in the evidence for this candidate.

    Enumerating what is *missing* is treated as a first-class output. A candidate
    that ranks first with four unanswered questions is a different object from one
    that ranks first with none, and the report must be able to say so.
    """

    question_id: str
    question: str
    why_it_matters: str
    resolving_test: str = Field(description="The concrete assay/data that would answer it.")
    blocking: bool = Field(
        default=False, description="True if a clinical conclusion cannot be drawn without it."
    )


class CandidatePair(FrozenModel):
    """Two variants proposed as jointly causal, plus everything known about them.

    ``variant_b`` is optional so that single-variant hypotheses (dominant, de novo,
    homozygous) share one ranked list with compound-het pairs. The challenge
    submission format encodes exactly this shape: a row with blank ``*_2`` columns
    is a single-variant proposal.
    """

    pair_id: str
    gene_symbol: str
    variant_a: VariantRecord
    variant_b: VariantRecord | None = None
    inheritance_model: InheritanceModel
    phase: PhaseEvidence

    scores: ComponentScores
    composite_score: float = Field(
        ge=0.0,
        le=1.0,
        description="Weighted combination. Reported alongside components, never instead of them.",
    )
    rank: int | None = Field(default=None, ge=1)

    supporting_evidence_ids: tuple[str, ...] = ()
    contradicting_evidence_ids: tuple[str, ...] = ()
    missing_evidence: tuple[OpenQuestion, ...] = ()
    recommended_next_test: str = Field(
        description="The single highest-value next experiment for THIS candidate."
    )
    discriminating_experiment: str | None = Field(
        default=None,
        description="What would distinguish this candidate from the next-ranked one.",
    )
    rank_rationale: str = Field(default="", description="Human-readable justification.")
    flags: tuple[str, ...] = Field(
        default=(), description="Soft markers, e.g. 'low_quality_call', 'phase_unknown'."
    )

    # ---------------------------------------------------------------- accessors

    @computed_field  # type: ignore[prop-decorator]
    @property
    def is_pair(self) -> bool:
        return self.variant_b is not None

    @property
    def variants(self) -> tuple[VariantRecord, ...]:
        return (self.variant_a,) if self.variant_b is None else (self.variant_a, self.variant_b)

    @property
    def variant_ids(self) -> tuple[str, ...]:
        return tuple(v.variant_id for v in self.variants)

    @property
    def has_contradictions(self) -> bool:
        return bool(self.contradicting_evidence_ids)

    @property
    def phase_is_disqualifying(self) -> bool:
        """In-cis pairs cannot be compound heterozygotes."""
        return (
            self.inheritance_model is InheritanceModel.COMPOUND_HETEROZYGOUS
            and self.phase.status in CIS_STATES
        )

    @property
    def blocking_questions(self) -> tuple[OpenQuestion, ...]:
        return tuple(q for q in self.missing_evidence if q.blocking)

    # ---------------------------------------------------------------- validators

    @model_validator(mode="after")
    def _both_variants_same_build(self) -> Self:
        if self.variant_b is not None:
            self.variant_a.coordinate.assert_same_build(self.variant_b.coordinate)
        return self

    @model_validator(mode="after")
    def _pair_variants_distinct(self) -> Self:
        if self.variant_b is not None and self.variant_a.variant_id == self.variant_b.variant_id:
            msg = (
                f"Pair {self.pair_id!r} lists the same variant twice "
                f"({self.variant_a.variant_id}). A variant cannot pair with itself; a "
                "homozygous call must use InheritanceModel.HOMOZYGOUS_RECESSIVE."
            )
            raise ValueError(msg)
        return self

    @model_validator(mode="after")
    def _compound_het_requires_two_hets(self) -> Self:
        if self.inheritance_model is InheritanceModel.COMPOUND_HETEROZYGOUS and (
            self.variant_b is None
        ):
            msg = f"Pair {self.pair_id!r} is COMPOUND_HETEROZYGOUS but has only one variant."
            raise ValueError(msg)
        return self

    def sort_key(self) -> tuple[float, tuple[int, int, str, str]]:
        """Deterministic ordering: score desc, then genomic position asc.

        The positional tiebreak is what makes repeat runs byte-identical; relying on
        dict/set iteration order here would break the determinism guarantee.
        """
        return (-self.composite_score, self.variant_a.coordinate.sort_key())


def make_pair_id(gene_symbol: str, variant_ids: tuple[str, ...]) -> str:
    """Deterministic, order-insensitive pair identifier."""
    body = "|".join(sorted(variant_ids))
    digest = hashlib.blake2b(f"{gene_symbol}::{body}".encode(), digest_size=6).hexdigest()
    return f"PAIR-{gene_symbol}-{digest}"
