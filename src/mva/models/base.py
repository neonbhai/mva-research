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

import hmac
import secrets
from enum import StrEnum

from pydantic import BaseModel, ConfigDict

# ---------------------------------------------------------------------------
# Non-disclosing identifiers for error messages (PRIV-09)
# ---------------------------------------------------------------------------

#: Per-process salt. Never persisted, never logged. A plain truncated hash of a
#: low-entropy value -- a coordinate, an HPO term, a short identifier -- is
#: brute-forceable in seconds, so the token is an HMAC under a random key that
#: dies with the process. Two messages in one run about the same record share a
#: token; nothing outside the run can reverse it.
_ERROR_TOKEN_SALT: bytes = secrets.token_bytes(16)


def error_token(value: object) -> str:
    """A stable-within-run, non-reversible handle for use in error messages.

    Exception messages travel further than anyone intends: to the terminal, into
    logs, into a crash report, and into an AI agent's context. A message that
    names `chr15:40200000 C>T 0/1` has disclosed a genotype to every one of those.
    This gives a debuggable handle instead -- `variant<a3f19c2b>` correlates across
    messages in a single run without carrying the record.
    """
    digest = hmac.new(_ERROR_TOKEN_SALT, str(value).encode("utf-8"), "sha256")
    return digest.hexdigest()[:8]


class FrozenModel(BaseModel):
    """Immutable, strictly-validated base for all domain records."""

    model_config = ConfigDict(
        frozen=True,
        extra="forbid",
        str_strip_whitespace=True,
        validate_default=True,
        use_enum_values=False,
        # PRIV-09. By default pydantic appends `input_value=...` to every
        # validation error, which for these models is the whole patient record --
        # coordinates, genotype, depths. A ValidationError is an ordinary
        # ValueError: nothing catches it, so the traceback reaches the terminal,
        # the logs, and an agent's context. The field name and the constraint that
        # failed are enough to debug with; the value is not ours to echo.
        hide_input_in_errors=True,
    )


class MutableModel(BaseModel):
    """Base for accumulator/builder objects that are never persisted directly."""

    model_config = ConfigDict(
        extra="forbid",
        str_strip_whitespace=True,
        hide_input_in_errors=True,  # PRIV-09, as above
    )


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
