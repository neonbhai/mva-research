"""The evidence ledger's atomic unit.

Central rule of this codebase: **no claim without an EvidenceItem**. Scores,
report sentences, drug rankings and mechanism links all cite evidence IDs. A
renderer that encounters an uncited assertion refuses to emit it.

Evidence is append-only and explicitly signed: it records not just *what* was
concluded but *by which tool at which version from which source*, and — unusually
but deliberately — what its **limitations** are. An evidence item with no stated
limitation is treated as a smell, not a strength.
"""

from __future__ import annotations

import hashlib
from datetime import datetime
from enum import StrEnum
from typing import Self

from pydantic import Field, model_validator

from mva.models.base import AssertionTier, FrozenModel


class EvidenceCategory(StrEnum):
    """What aspect of the hypothesis this evidence speaks to."""

    ANALYTICAL = "analytical"
    """Is the call technically real? Depth, allele balance, mappability, FILTER."""

    POPULATION = "population"
    """Is it rare enough? Always paired with a named population and dataset version."""

    INHERITANCE = "inheritance"
    """Segregation, phase, de-novo status, parental genotypes."""

    CONSEQUENCE = "consequence"
    """Predicted molecular effect on transcript/protein."""

    PHENOTYPE = "phenotype"
    """Match between gene-associated phenotype and observed HPO profile."""

    MECHANISM = "mechanism"
    """Biological pathway/process linkage."""

    DRUG = "drug"
    """Pharmacological intervention evidence."""

    SAFETY = "safety"
    """Toxicity, contraindication, pediatric exposure, oncogenic risk."""

    CONTRADICTION = "contradiction"
    """Evidence that actively opposes a hypothesis. Never discarded."""

    PROVENANCE = "provenance"
    """Bookkeeping about how an artifact was produced."""


class EvidenceDirection(StrEnum):
    """Whether the evidence argues for or against its subject."""

    SUPPORTS = "supports"
    CONTRADICTS = "contradicts"
    NEUTRAL = "neutral"


class EvidenceStrength(StrEnum):
    """Weight of the evidence.

    Ordinal, not interval: ``STRONG`` is not "3x MODERATE". Numeric mapping happens
    only in the scoring layer, where the mapping is configuration and is documented
    as a heuristic.
    """

    DEFINITIVE = "definitive"
    STRONG = "strong"
    MODERATE = "moderate"
    SUPPORTING = "supporting"
    WEAK = "weak"
    INSUFFICIENT = "insufficient"


class EvidenceType(StrEnum):
    """The empirical modality that produced the evidence.

    Ordered roughly by how far it is from a human clinical outcome. This is what
    stops "a compound binds the target in vitro" from being scored like "a
    randomised trial in children showed benefit".
    """

    DIRECT_MEASUREMENT = "direct_measurement"
    """Measured in this patient's own data."""

    CURATED_DATABASE = "curated_database"
    IN_SILICO_PREDICTION = "in_silico_prediction"
    BIOCHEMICAL_BINDING = "biochemical_binding"
    CELL_LINE = "cell_line"
    PRIMARY_PATIENT_CELLS = "primary_patient_cells"
    ANIMAL_MODEL = "animal_model"
    HUMAN_CASE_REPORT = "human_case_report"
    HUMAN_OBSERVATIONAL = "human_observational"
    HUMAN_TRIAL = "human_trial"
    EXPERT_REVIEW = "expert_review"
    PIPELINE_INFERENCE = "pipeline_inference"
    """Produced by this pipeline's own logic. Never treated as external support."""


#: Evidence types that are predictions or bookkeeping, not observations of biology.
NON_EMPIRICAL_TYPES: frozenset[EvidenceType] = frozenset(
    {EvidenceType.IN_SILICO_PREDICTION, EvidenceType.PIPELINE_INFERENCE}
)

#: Evidence types demonstrating an effect in a whole living organism.
IN_VIVO_TYPES: frozenset[EvidenceType] = frozenset(
    {
        EvidenceType.ANIMAL_MODEL,
        EvidenceType.HUMAN_CASE_REPORT,
        EvidenceType.HUMAN_OBSERVATIONAL,
        EvidenceType.HUMAN_TRIAL,
    }
)


class Citation(FrozenModel):
    """A resolvable pointer to the source of a claim."""

    source: str = Field(description="Database or publication venue, e.g. 'gnomAD', 'PubMed'")
    identifier: str = Field(description="Accession/PMID/DOI/URL fragment")
    version: str | None = Field(
        default=None,
        description="Dataset release or access date. Required for database assertions.",
    )
    url: str | None = None
    title: str | None = None

    @property
    def key(self) -> str:
        return f"{self.source}:{self.identifier}" + (f"@{self.version}" if self.version else "")


class EvidenceItem(FrozenModel):
    """One structured, attributable, falsifiable claim.

    ``subject_id`` is a free-form pointer to whatever the claim is *about*: a
    canonical variant ID, a pair ID, a gene symbol, a drug ID, a mechanism ID. The
    evidence store indexes on it, which is what allows a report renderer to ask
    "show me everything known about this pair, including what argues against it".
    """

    evidence_id: str = Field(description="Deterministic, content-derived. See make_evidence_id.")
    subject_id: str = Field(min_length=1)
    subject_kind: str = Field(
        description="variant | pair | gene | phenotype | mechanism | drug | run | artifact"
    )
    claim: str = Field(min_length=3, description="One falsifiable sentence.")
    category: EvidenceCategory
    direction: EvidenceDirection
    strength: EvidenceStrength
    evidence_type: EvidenceType
    tier: AssertionTier
    citation: Citation | None = None
    method: str = Field(description="How the claim was produced, concretely.")
    tool: str = Field(description="Software or database that produced it.")
    tool_version: str
    limitations: str = Field(
        min_length=3,
        description=(
            "What this evidence does NOT establish. Mandatory. An in-silico score with "
            "no stated limitation is how prediction gets mistaken for proof."
        ),
    )
    timestamp: datetime
    run_id: str | None = None
    numeric_value: float | None = Field(
        default=None, description="Optional machine-readable magnitude (AF, score, depth)."
    )
    payload: dict[str, str | int | float | bool | None] = Field(
        default_factory=dict,
        description="Structured detail. Must never carry raw patient identifiers.",
    )

    @model_validator(mode="after")
    def _database_claims_need_versioned_citation(self) -> Self:
        """A database assertion without a version is not reproducible.

        Population frequency is the canonical trap: 'AF = 0.0001' is meaningless
        without knowing gnomAD v2.1.1 vs v4.1 and which population.
        """
        if self.tier is AssertionTier.DATABASE_ASSERTION and (
            self.citation is None or self.citation.version is None
        ):
            msg = (
                f"Evidence {self.evidence_id!r} is a DATABASE_ASSERTION but lacks a "
                "versioned citation. Record source + version, or downgrade the tier."
            )
            raise ValueError(msg)
        return self

    @model_validator(mode="after")
    def _observed_data_must_be_measured(self) -> Self:
        """Guard against laundering a prediction as an observation."""
        if self.tier is AssertionTier.OBSERVED_DATA and self.evidence_type in NON_EMPIRICAL_TYPES:
            msg = (
                f"Evidence {self.evidence_id!r} claims tier OBSERVED_DATA but has "
                f"evidence_type {self.evidence_type.value!r}. A prediction is not an observation."
            )
            raise ValueError(msg)
        return self

    @property
    def is_contradiction(self) -> bool:
        return self.direction is EvidenceDirection.CONTRADICTS


def make_evidence_id(
    *,
    subject_id: str,
    category: EvidenceCategory,
    claim: str,
    tool: str,
) -> str:
    """Deterministic evidence ID.

    Content-derived rather than random so that a re-run of the same pipeline on the
    same input yields byte-identical evidence tables (see the determinism tests).
    Timestamp is deliberately excluded from the hash: the same conclusion drawn
    twice is the same evidence, not two.
    """
    digest = hashlib.blake2b(
        "\x1f".join([subject_id, category.value, claim, tool]).encode("utf-8"),
        digest_size=8,
    ).hexdigest()
    return f"EV-{category.value[:4].upper()}-{digest}"
