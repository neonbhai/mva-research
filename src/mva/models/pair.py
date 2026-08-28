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

from mva.models.base import FrozenModel, error_token
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

#: Stand-in coordinate key for the absent second variant of a single-variant
#: hypothesis, so that :meth:`CandidatePair.sort_key` is defined on both shapes.
#:
#: It sorts before every real coordinate and cannot be produced by one:
#: :func:`~mva.models.genome.contig_sort_key` is non-negative and
#: ``GenomicCoordinate.position`` is constrained ``> 0``. A single-variant
#: candidate therefore orders **ahead of** a pair sharing its first variant at an
#: equal composite score — a documented rule rather than an accident of whichever
#: one the caller listed first. :class:`mva.prioritization.pairing.PairCandidate`
#: uses the same sentinel for the same reason.
NO_SECOND_VARIANT: tuple[int, int, str, str] = (-1, -1, "", "")


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
        """A pairing defect must not pay for itself with a disclosure (PRIV-09).

        ``pair_id`` stays in the message and the variant ID does not, which is the
        whole distinction: :func:`make_pair_id` is a blake2b digest that is already
        written into every artifact this run produces, so naming it costs nothing
        and is the only way to tie the failure to the candidate. A ``variant_id``
        is a coordinate and an allele pair — under this project's threat model a
        handful of those identifies an individual and their parents — and an
        ordinary ``ValueError`` carries it into the terminal, the log file, a crash
        report and an agent's context.
        """
        if self.variant_b is not None and self.variant_a.variant_id == self.variant_b.variant_id:
            msg = (
                f"Pair {self.pair_id!r} lists the same variant as both of its members "
                f"(variant handle <variant:{error_token(self.variant_a.variant_id)}>; the "
                "coordinate is tokenised rather than echoed, PRIV-09). A variant cannot "
                "pair with itself: two copies of one allele is a homozygous call, and it "
                "must be built as InheritanceModel.HOMOZYGOUS_RECESSIVE with variant_b "
                "left None. This normally means the pair generator emitted a self-pair "
                "for a gene holding a single qualifying variant."
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

    def sort_key(self) -> tuple[float, tuple[int, int, str, str], tuple[int, int, str, str], str]:
        """A **total** order over candidates. Four components, in this order:

        1. ``-composite_score`` — the scientific judgement, descending. This is and
           stays the primary component; everything below it is a tiebreak, never a
           re-ranking. Two candidates with different scores are ordered by score
           alone, exactly as before.
        2. ``variant_a.coordinate.sort_key()`` — first variant, by karyotype
           position ascending, so a tied block reads down the chromosome rather
           than in whatever order it was assembled.
        3. the second variant's coordinate, or :data:`NO_SECOND_VARIANT` when there
           is none. **This component is why the key is total.** Without it, two
           pairs sharing a first variant and a composite score compare *equal*, so
           every sort over them is merely stable and their order is decided by the
           caller's input order — which changes the rendered submission CSV, swaps
           which row receives the higher EPCR, and can cost 50 rank points if the
           demoted row is the true answer. The sentinel sorts below every real
           coordinate, so a single-variant candidate deterministically precedes a
           pair sharing its first variant.
        4. ``pair_id`` — the terminator. :func:`make_pair_id` derives it from the
           gene symbol and the sorted variant ids via blake2b, so it is unique per
           hypothesis and identical across processes; two distinct candidates can
           never compare equal. It is a content digest, never :func:`hash`, which
           varies with ``PYTHONHASHSEED`` (GP-30).

        Totality is the requirement, not tidiness. A key that admits ties promotes
        input order to a scientific decision and breaks the byte-identical
        repeat-run guarantee (GP-30); relying on dict or set iteration order here
        would do the same.
        """
        second = NO_SECOND_VARIANT if self.variant_b is None else self.variant_b.sort_key()
        return (
            -self.composite_score,
            self.variant_a.coordinate.sort_key(),
            second,
            self.pair_id,
        )


def make_pair_id(gene_symbol: str, variant_ids: tuple[str, ...]) -> str:
    """Deterministic, order-insensitive pair identifier."""
    body = "|".join(sorted(variant_ids))
    digest = hashlib.blake2b(f"{gene_symbol}::{body}".encode(), digest_size=6).hexdigest()
    return f"PAIR-{gene_symbol}-{digest}"
