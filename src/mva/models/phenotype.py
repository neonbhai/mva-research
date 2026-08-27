"""Phenotype representation.

The single most important distinction in this module: **"not assessed" is not
"absent"**. Treating a missing phenotype as a negative finding silently
manufactures evidence against candidate genes, and is one of the classic ways a
rare-disease pipeline arrives at a confident wrong answer.
"""

from __future__ import annotations

import re
from enum import StrEnum

from pydantic import Field, field_validator

from mva.models.base import FrozenModel

_HPO_RE = re.compile(r"^HP:\d{7}$")


class ObservationStatus(StrEnum):
    """Four-valued phenotype logic. Never collapse to a boolean."""

    OBSERVED = "observed"
    """Clinician recorded the feature as present."""

    EXCLUDED = "excluded"
    """Clinician explicitly recorded the feature as ABSENT after assessment.
    This is real, usable negative evidence."""

    UNCERTAIN = "uncertain"
    """Assessed, but the finding was equivocal."""

    NOT_ASSESSED = "not_assessed"
    """No information. Carries NO evidential weight in either direction."""


#: Statuses that may contribute negative evidence to phenotype scoring.
NEGATIVE_EVIDENCE_STATUSES: frozenset[ObservationStatus] = frozenset({ObservationStatus.EXCLUDED})

#: Statuses that must be treated as absence of information, not information.
UNINFORMATIVE_STATUSES: frozenset[ObservationStatus] = frozenset(
    {ObservationStatus.NOT_ASSESSED, ObservationStatus.UNCERTAIN}
)


class Onset(StrEnum):
    """Coarse onset bucketing (HPO onset subontology, simplified)."""

    ANTENATAL = "antenatal"
    CONGENITAL = "congenital"
    NEONATAL = "neonatal"
    INFANTILE = "infantile"
    CHILDHOOD = "childhood"
    JUVENILE = "juvenile"
    ADULT = "adult"
    UNKNOWN = "unknown"


class PhenotypeObservation(FrozenModel):
    """One HPO term as recorded for one individual."""

    hpo_id: str
    label: str = Field(min_length=1)
    status: ObservationStatus
    onset: Onset = Onset.UNKNOWN
    provenance: str = Field(
        description="Where this came from: 'clinical_summary', 'phenopacket', 'synthetic_fixture'."
    )
    extraction_confidence: float = Field(
        ge=0.0,
        le=1.0,
        description=(
            "Confidence that the term was correctly extracted from the source. This is "
            "an EXTRACTION confidence, not a clinical certainty; the two are different "
            "and must not be multiplied together without saying so."
        ),
    )
    source_excerpt_hash: str | None = Field(
        default=None,
        description=(
            "Hash of the source text, never the text itself. Free-text clinical notes "
            "are directly identifying and must not enter derived artifacts."
        ),
    )
    notes: str | None = None

    @field_validator("hpo_id")
    @classmethod
    def _valid_hpo(cls, value: str) -> str:
        token = value.strip().upper().replace("HP_", "HP:")
        if not _HPO_RE.match(token):
            msg = f"Invalid HPO identifier {value!r}; expected the form 'HP:0001250'."
            raise ValueError(msg)
        return token

    @property
    def is_present(self) -> bool:
        return self.status is ObservationStatus.OBSERVED

    @property
    def is_explicitly_absent(self) -> bool:
        """True ONLY for EXCLUDED. NOT_ASSESSED deliberately returns False."""
        return self.status is ObservationStatus.EXCLUDED

    @property
    def is_informative(self) -> bool:
        return self.status not in UNINFORMATIVE_STATUSES


class PhenotypeProfile(FrozenModel):
    """The full phenotype picture for one individual."""

    subject_id: str
    observations: tuple[PhenotypeObservation, ...]
    source_artifact: str
    hpo_version: str = Field(description="HPO release used for term validation/labels.")

    @property
    def observed_terms(self) -> tuple[str, ...]:
        return tuple(o.hpo_id for o in self.observations if o.is_present)

    @property
    def excluded_terms(self) -> tuple[str, ...]:
        return tuple(o.hpo_id for o in self.observations if o.is_explicitly_absent)

    @property
    def not_assessed_terms(self) -> tuple[str, ...]:
        return tuple(
            o.hpo_id for o in self.observations if o.status is ObservationStatus.NOT_ASSESSED
        )

    def status_of(self, hpo_id: str) -> ObservationStatus:
        """Status lookup that returns NOT_ASSESSED for unknown terms.

        This is the safe default: a term nobody recorded is a term nobody assessed.
        """
        for obs in self.observations:
            if obs.hpo_id == hpo_id:
                return obs.status
        return ObservationStatus.NOT_ASSESSED
