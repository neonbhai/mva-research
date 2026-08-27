"""Parse a phenotype TSV into a typed :class:`PhenotypeProfile` (GP-02).

The whole risk surface of this file is one column: ``status``. Real clinical
extracts spell negation a dozen ways (``absent``, ``negated``, ``no``,
``not_present``, ``ruled out``), and a reader that does not recognise a spelling
has exactly two options. It can guess a default — which either invents a finding
the clinician never recorded, or erases a negative finding that is genuine,
usable evidence — or it can refuse. This module refuses (GP-14).

The second trap is the distinction the whole package exists to preserve:
``not_assessed`` is **not** ``excluded``. "Nobody looked" and "somebody looked and
it was not there" are different facts with opposite evidential consequences, and
they are kept apart from the first line of parsing through to the final score.

**Privacy.** The ``notes`` column is carried through verbatim because the
synthetic fixture uses it, but it is free text: anything derived from a real
clinical narrative belongs in ``source_excerpt_hash``, never in ``notes``
(PRIV-09). Error messages here are built from file names, line numbers and
controlled-vocabulary tokens only.
"""

from __future__ import annotations

from pathlib import Path
from typing import Final

from mva.errors import IngestionError
from mva.models.phenotype import (
    ObservationStatus,
    Onset,
    PhenotypeObservation,
    PhenotypeProfile,
)
from mva.phenotype.hpo import normalise_hpo_id, read_tsv_rows

#: Columns a phenotype TSV must provide. Everything else is optional.
REQUIRED_COLUMNS: Final[tuple[str, ...]] = ("hpo_id", "label", "status")

#: Optional columns and the documented behaviour when they are absent.
OPTIONAL_COLUMNS: Final[tuple[str, ...]] = (
    "onset",
    "provenance",
    "extraction_confidence",
    "notes",
)

#: Accepted spellings of each status, normalised to lower-case with ``-``, ``/``,
#: ``.`` and whitespace folded to ``_``.
#:
#: This table is the negation-handling contract. It is deliberately explicit and
#: deliberately finite: an unlisted spelling raises rather than falling back to a
#: default, because both plausible defaults are wrong.
#:
#: * ``observed`` — the clinician recorded the feature as present.
#: * ``excluded`` — the clinician assessed the feature and recorded it ABSENT.
#:   This is the negation bucket, and the only status that may contribute
#:   negative evidence downstream.
#: * ``uncertain`` — assessed, equivocal. Contributes nothing in either direction.
#: * ``not_assessed`` — no information at all. Contributes nothing in either
#:   direction, and must never be read as absence.
#:
#: Note where the ambiguous tokens land. ``unknown``, ``na`` and ``missing`` map to
#: ``not_assessed`` (no information), never to ``excluded``; ``suspected`` and
#: ``probable`` map to ``uncertain``, never to ``observed``. In both cases the
#: mapping is chosen so that a mis-parse loses evidence rather than manufacturing
#: it.
STATUS_ALIASES: Final[dict[str, ObservationStatus]] = {
    # --- present ----------------------------------------------------------
    "observed": ObservationStatus.OBSERVED,
    "present": ObservationStatus.OBSERVED,
    "positive": ObservationStatus.OBSERVED,
    "affected": ObservationStatus.OBSERVED,
    "yes": ObservationStatus.OBSERVED,
    # --- explicitly absent (negation) -------------------------------------
    "excluded": ObservationStatus.EXCLUDED,
    "absent": ObservationStatus.EXCLUDED,
    "negated": ObservationStatus.EXCLUDED,
    "negative": ObservationStatus.EXCLUDED,
    "no": ObservationStatus.EXCLUDED,
    "not_present": ObservationStatus.EXCLUDED,
    "not_observed": ObservationStatus.EXCLUDED,
    "ruled_out": ObservationStatus.EXCLUDED,
    "denied": ObservationStatus.EXCLUDED,
    # --- assessed but equivocal -------------------------------------------
    "uncertain": ObservationStatus.UNCERTAIN,
    "equivocal": ObservationStatus.UNCERTAIN,
    "unclear": ObservationStatus.UNCERTAIN,
    "questionable": ObservationStatus.UNCERTAIN,
    "borderline": ObservationStatus.UNCERTAIN,
    "possible": ObservationStatus.UNCERTAIN,
    "probable": ObservationStatus.UNCERTAIN,
    "suspected": ObservationStatus.UNCERTAIN,
    # --- no information ---------------------------------------------------
    "not_assessed": ObservationStatus.NOT_ASSESSED,
    "unassessed": ObservationStatus.NOT_ASSESSED,
    "not_evaluated": ObservationStatus.NOT_ASSESSED,
    "not_examined": ObservationStatus.NOT_ASSESSED,
    "no_information": ObservationStatus.NOT_ASSESSED,
    "unknown": ObservationStatus.NOT_ASSESSED,
    "missing": ObservationStatus.NOT_ASSESSED,
    "na": ObservationStatus.NOT_ASSESSED,
    "n_a": ObservationStatus.NOT_ASSESSED,
    "nd": ObservationStatus.NOT_ASSESSED,
}

#: Extraction confidence recorded when the source file omits the column.
#:
#: 1.0 is correct *here specifically* and would be dishonest elsewhere: this is
#: confidence that the term was transcribed correctly out of a structured tabular
#: file, where no natural-language extraction step exists to be wrong about. It
#: says nothing about clinical certainty — the two are different quantities and
#: :class:`PhenotypeObservation` documents why they must not be multiplied.
DEFAULT_EXTRACTION_CONFIDENCE: Final[float] = 1.0


def parse_status(raw: str) -> ObservationStatus:
    """Map a raw status token to an :class:`ObservationStatus`, or raise.

    Case-insensitive; ``-``, ``/``, ``.`` and whitespace are folded to ``_``, so
    ``"Not Assessed"``, ``"NOT-ASSESSED"`` and ``"not_assessed"`` are one token, and
    ``"n/a"`` normalises to ``"n_a"``.

    Raises :class:`~mva.errors.IngestionError` on any unrecognised value, including
    the empty string. There is no default. A four-valued logic with a fallback
    branch is a three-valued logic wearing a disguise, and the missing value is
    always the one that matters.
    """
    token = _normalise_status_token(raw)
    status = STATUS_ALIASES.get(token)
    if status is None:
        msg = (
            f"Unrecognised phenotype status {raw.strip()!r}. Refusing to default: "
            "guessing 'observed' invents a finding, and guessing 'excluded' turns a "
            "missing assessment into negative evidence against a candidate gene "
            "(GP-14). Canonical values are "
            f"{', '.join(sorted(member.value for member in ObservationStatus))}. "
            f"Accepted spellings: {allowed_status_summary()}."
        )
        raise IngestionError(msg)
    return status


def allowed_status_summary() -> str:
    """Human-readable listing of every accepted status spelling, grouped."""
    groups: dict[str, list[str]] = {}
    for alias, status in STATUS_ALIASES.items():
        groups.setdefault(status.value, []).append(alias)
    return "; ".join(
        f"{status} <- {', '.join(sorted(aliases))}" for status, aliases in sorted(groups.items())
    )


def parse_onset(raw: str) -> Onset:
    """Map a raw onset token to an :class:`Onset`, or raise.

    An empty cell is :attr:`Onset.UNKNOWN` — the enum has an explicit member for
    "not stated", so no information is representable without inventing a bucket.
    An unrecognised non-empty token still raises: ``"neonatal?"`` silently becoming
    ``unknown`` would lose a real clinical detail without telling anyone.
    """
    token = _normalise_status_token(raw)
    if not token:
        return Onset.UNKNOWN
    try:
        return Onset(token)
    except ValueError as exc:
        allowed = ", ".join(sorted(onset.value for onset in Onset))
        msg = f"Unrecognised onset {raw.strip()!r}; allowed values are: {allowed}."
        raise IngestionError(msg) from exc


def load_phenotype_profile(
    path: Path,
    *,
    subject_id: str,
    hpo_version: str,
    source_artifact: str,
) -> PhenotypeProfile:
    """Load a phenotype TSV into a validated, deterministically ordered profile.

    Lines beginning with ``#`` are comments; the first remaining line is the
    header. ``hpo_id``, ``label`` and ``status`` are required; ``onset``,
    ``provenance``, ``extraction_confidence`` and ``notes`` are optional.

    Observations are sorted by HPO identifier so that two runs over the same file
    produce the same object and the same artifact bytes (GP-30). Duplicate
    identifiers are rejected: two rows for one term either agree (and one is
    redundant) or disagree (and no reader can pick the right one), and a duplicate
    would additionally be double-counted in the scoring denominator.
    """
    rows = read_tsv_rows(path, required_columns=REQUIRED_COLUMNS)

    observations: list[PhenotypeObservation] = []
    seen: dict[str, int] = {}

    for lineno, row in rows:
        location = f"{path.name} line {lineno}"
        hpo_id = normalise_hpo_id(row["hpo_id"], context=location)
        if hpo_id in seen:
            msg = (
                f"{location}: duplicate HPO term {hpo_id} (first seen on line {seen[hpo_id]}). "
                "Two rows for one term cannot be reconciled by a parser; resolve the "
                "conflict in the source file."
            )
            raise IngestionError(msg)
        seen[hpo_id] = lineno
        observations.append(
            _build_observation(
                row,
                hpo_id=hpo_id,
                location=location,
                source_artifact=source_artifact,
            )
        )

    observations.sort(key=lambda obs: obs.hpo_id)
    return PhenotypeProfile(
        subject_id=subject_id,
        observations=tuple(observations),
        source_artifact=source_artifact,
        hpo_version=hpo_version,
    )


def _build_observation(
    row: dict[str, str],
    *,
    hpo_id: str,
    location: str,
    source_artifact: str,
) -> PhenotypeObservation:
    """Turn one validated row into a typed observation."""
    label = row["label"].strip()
    if not label:
        msg = f"{location}: HPO term {hpo_id} has an empty label."
        raise IngestionError(msg)

    notes = row.get("notes", "").strip()
    provenance = row.get("provenance", "").strip() or source_artifact

    return PhenotypeObservation(
        hpo_id=hpo_id,
        label=label,
        status=parse_status(row["status"]),
        onset=parse_onset(row.get("onset", "")),
        provenance=provenance,
        extraction_confidence=_parse_confidence(row.get("extraction_confidence", ""), location),
        source_excerpt_hash=None,
        notes=notes or None,
    )


def _parse_confidence(raw: str, location: str) -> float:
    """Parse the extraction-confidence column, defaulting only when it is absent."""
    token = raw.strip()
    if not token:
        return DEFAULT_EXTRACTION_CONFIDENCE
    try:
        value = float(token)
    except ValueError as exc:
        msg = f"{location}: extraction_confidence {token!r} is not a number."
        raise IngestionError(msg) from exc
    if not 0.0 <= value <= 1.0:
        msg = f"{location}: extraction_confidence {value} is outside [0, 1]."
        raise IngestionError(msg)
    return value


def _normalise_status_token(raw: str) -> str:
    """Fold case and separator punctuation so one concept is one token."""
    token = raw.strip().lower()
    for char in (" ", "\t", "-", "/", "."):
        token = token.replace(char, "_")
    while "__" in token:
        token = token.replace("__", "_")
    return token.strip("_")
