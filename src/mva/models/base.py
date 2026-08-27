"""Shared base types for every mva domain model.

Design rules enforced here:

* Domain models are **frozen**. Pipeline stages produce new objects rather than
  mutating upstream ones, so an artifact written to disk can never drift from the
  object that was hashed into the provenance record.
* ``extra="forbid"``. A typo in a config key or an adapter payload is a loud error,
  not a silently dropped field. In a clinical-adjacent context, silent field loss is
  the single most dangerous default a schema can have.
* Enums are ``StrEnum`` so Parquet/DuckDB/CSV round-trips are lossless and readable.
"""

from __future__ import annotations

from enum import StrEnum

from pydantic import BaseModel, ConfigDict


class FrozenModel(BaseModel):
    """Immutable, strictly-validated base for all domain records."""

    model_config = ConfigDict(
        frozen=True,
        extra="forbid",
        str_strip_whitespace=True,
        validate_default=True,
        use_enum_values=False,
    )


class MutableModel(BaseModel):
    """Base for accumulator/builder objects that are never persisted directly."""

    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)


class Sensitivity(StrEnum):
    """Data-handling classification attached to every artifact.

    ``SENSITIVE`` artifacts are patient-derived and may never leave the external
    workspace. ``PUBLIC`` artifacts are safe to commit or submit. ``DERIVED_SAFE``
    is the deliberately narrow middle tier: computed from sensitive input but
    provably free of patient genotypes (e.g. aggregate QC counts). Promotion to
    ``PUBLIC`` always requires an explicit export gate, never an implicit default.
    """

    SENSITIVE = "sensitive"
    DERIVED_SAFE = "derived_safe"
    PUBLIC = "public"


class AssertionTier(StrEnum):
    """How a statement in a report came to be believed.

    Reports must label every claim with one of these. The ordering is deliberate:
    it is the epistemic ladder from "we measured it" down to "we are guessing".
    A report renderer refuses to emit an unlabelled assertion.
    """

    OBSERVED_DATA = "observed_data"
    """Directly present in the patient's own data (e.g. this genotype was called)."""

    DATABASE_ASSERTION = "database_assertion"
    """A curated third-party database says so (ClinVar, gnomAD, OMIM)."""

    COMPUTATIONAL_PREDICTION = "computational_prediction"
    """An in-silico tool predicted it (VEP consequence, SpliceAI, CADD)."""

    LITERATURE_MECHANISM = "literature_mechanism"
    """Published experimental work supports the mechanistic link."""

    INFERENCE = "inference"
    """Derived by this pipeline's logic from the tiers above."""

    SPECULATION = "speculation"
    """Plausible, unsupported. Must be visually marked in every rendered report."""


#: Tiers that may never be presented as established fact without a hedge.
UNPROVEN_TIERS: frozenset[AssertionTier] = frozenset(
    {AssertionTier.INFERENCE, AssertionTier.SPECULATION}
)
