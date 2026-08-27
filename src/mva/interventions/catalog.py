"""Load the local drug catalogue.

The catalogue is the intervention stage's only source of pharmacological fact,
and it is deliberately **signed**: every row carries the direction the agent
pushes its target, not merely the target it binds. That single column is what
lets the stage reject a compound acting on the right node the wrong way (GP-16),
which a target-proximity search structurally cannot do.

Parsing happens once, here, at the boundary, into a typed row (GP-02). Nothing
downstream ever sees a dict of strings. Three parsing rules matter:

* A blank cell becomes ``None``, never a default. "Not recorded" and "recorded as
  zero" are different facts and stay different (GP-14) — `worsens_cin` blank means
  *unassessed*, which is a blocking gap, not a clean bill of health.
* Enum cells are parsed strictly. An unrecognised approval status is an error, not
  an `UNKNOWN` fallback, because silently downgrading it would hide a data fault.
* Numeric cells are parsed to float and rejected if negative.

GP-20: the shipped catalogue is **synthetic**; every agent in it is fictional.
"""

from __future__ import annotations

import csv
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

from mva.errors import IngestionError
from mva.models.drug import ApprovalStatus, InterventionClass
from mva.models.evidence import EvidenceType
from mva.models.mechanism import EffectDirection

__all__ = ["CATALOG_COLUMNS", "CatalogEntry", "DrugCatalog"]

#: The catalogue schema, in file order. The column names are the contract.
CATALOG_COLUMNS: tuple[str, ...] = (
    "drug_id",
    "name",
    "approved_name",
    "approval_status",
    "intervention_class",
    "target",
    "target_node_id",
    "mechanism_of_action",
    "observed_direction",
    "is_direct_evidence",
    "strongest_evidence_type",
    "has_pediatric_exposure",
    "youngest_age_studied",
    "pediatric_indication",
    "route",
    "cns_penetrant",
    "achievable_plasma_um",
    "required_effective_um",
    "half_life_hours",
    "worsens_cin",
    "validation_experiment",
    "notes",
)

_TRUE_TOKENS: frozenset[str] = frozenset({"1", "true", "yes", "y", "t"})
_FALSE_TOKENS: frozenset[str] = frozenset({"0", "false", "no", "n", "f"})

#: Route values that mean "this agent has never been given to a living subject".
_NO_ROUTE_TOKENS: frozenset[str] = frozenset({"", "none", "n/a", "not applicable"})


@dataclass(frozen=True)
class CatalogEntry:
    """One typed catalogue row.

    A faithful mirror of the TSV: this type adds no judgement and drops no field.
    Interpretation — direction agreement, safety, evidence quality — happens in
    the sibling modules, so that a reader can always separate *what the catalogue
    says* from *what the pipeline concluded*.
    """

    drug_id: str
    name: str
    approved_name: str | None
    approval_status: ApprovalStatus
    intervention_class: InterventionClass
    target: str
    target_node_id: str
    mechanism_of_action: str
    observed_direction: EffectDirection
    is_direct_evidence: bool
    strongest_evidence_type: EvidenceType
    has_pediatric_exposure: bool
    youngest_age_studied: str | None
    pediatric_indication: str | None
    route: str | None
    cns_penetrant: bool | None
    achievable_plasma_um: float | None
    required_effective_um: float | None
    half_life_hours: float | None
    worsens_cin: bool | None
    validation_experiment: str
    notes: str

    @property
    def is_repurposable(self) -> bool:
        """Whether the agent is approved somewhere, i.e. eligible for a repurposing claim."""
        return self.approval_status in {
            ApprovalStatus.APPROVED,
            ApprovalStatus.APPROVED_OTHER_INDICATION,
            ApprovalStatus.APPROVED_ADULT_ONLY,
        }

    @property
    def has_administration_route(self) -> bool:
        """False when the catalogue records no viable route of administration."""
        return (self.route or "").strip().lower() not in _NO_ROUTE_TOKENS

    @property
    def concentration_achievable(self) -> bool | None:
        """Tri-state. ``None`` when either concentration figure is missing."""
        if self.achievable_plasma_um is None or self.required_effective_um is None:
            return None
        return self.achievable_plasma_um >= self.required_effective_um


class DrugCatalog:
    """An immutable, versioned set of catalogue rows."""

    def __init__(self, entries: Sequence[CatalogEntry], *, version: str) -> None:
        if not version.strip():
            msg = "DrugCatalog requires a non-empty version string (GP-20, GP-31)."
            raise IngestionError(msg)
        self._version = version
        self._entries: tuple[CatalogEntry, ...] = tuple(
            sorted(entries, key=lambda entry: entry.drug_id)
        )
        seen: set[str] = set()
        for entry in self._entries:
            if entry.drug_id in seen:
                msg = f"Drug catalogue contains a duplicate drug_id: {entry.drug_id}."
                raise IngestionError(msg)
            seen.add(entry.drug_id)

    @classmethod
    def from_tsv(cls, path: Path, *, version: str) -> DrugCatalog:
        """Parse a `#`-commented TSV catalogue."""
        if not path.is_file():
            msg = f"Drug catalogue not found: {path.name}"
            raise IngestionError(msg)
        lines = [
            line
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        if not lines:
            msg = f"Drug catalogue {path.name} contains no header row."
            raise IngestionError(msg)
        rows = [tuple(cell.strip() for cell in row) for row in csv.reader(lines, delimiter="\t")]
        header = rows[0]
        missing = [name for name in CATALOG_COLUMNS if name not in header]
        if missing:
            msg = (
                f"Drug catalogue {path.name} is missing {len(missing)} required "
                f"column(s): {', '.join(missing)}."
            )
            raise IngestionError(msg)
        index = {name: position for position, name in enumerate(header)}
        entries = [
            _build_entry(row, index=index, ctx=f"{path.name} row {number}")
            for number, row in enumerate(rows[1:], start=1)
        ]
        return cls(entries, version=version)

    def entries(self) -> tuple[CatalogEntry, ...]:
        """Every row, ordered by `drug_id` (GP-30)."""
        return self._entries

    def for_target_node(self, node_id: str) -> tuple[CatalogEntry, ...]:
        """Rows whose target node is exactly `node_id`.

        Exact match only. There is no fuzzy "nearby node" fallback: pathway
        proximity is precisely the heuristic that surfaces contraindicated agents,
        and it is not going to be smuggled in through a lookup helper.
        """
        return tuple(entry for entry in self._entries if entry.target_node_id == node_id)

    @property
    def version(self) -> str:
        """Version of the catalogue table behind this object."""
        return self._version

    def __len__(self) -> int:
        return len(self._entries)


# --------------------------------------------------------------------- parsing


def _cell(row: Sequence[str], *, index: dict[str, int], name: str) -> str:
    position = index[name]
    return row[position].strip() if position < len(row) else ""


def _require(value: str, *, field: str, ctx: str) -> str:
    if not value:
        msg = f"{ctx}: required field {field!r} is empty."
        raise IngestionError(msg)
    return value


def _parse_enum[E: StrEnum](enum_cls: type[E], raw: str, *, field: str, ctx: str) -> E:
    """Strict enum parse. The offending value is not echoed (PRIV-09)."""
    try:
        return enum_cls(raw)
    except ValueError as exc:
        allowed = ", ".join(sorted(str(member.value) for member in enum_cls))
        msg = f"{ctx}: field {field!r} holds an unrecognised value. Allowed: {allowed}."
        raise IngestionError(msg) from exc


def _parse_bool(raw: str, *, field: str, ctx: str) -> bool:
    lowered = raw.lower()
    if lowered in _TRUE_TOKENS:
        return True
    if lowered in _FALSE_TOKENS:
        return False
    msg = f"{ctx}: field {field!r} is not a recognised boolean. Use 1/0."
    raise IngestionError(msg)


def _parse_optional_bool(raw: str, *, field: str, ctx: str) -> bool | None:
    """Blank means *unassessed*, which is a fact in its own right (GP-14)."""
    return None if not raw else _parse_bool(raw, field=field, ctx=ctx)


def _parse_optional_float(raw: str, *, field: str, ctx: str) -> float | None:
    if not raw:
        return None
    try:
        value = float(raw)
    except ValueError as exc:
        msg = f"{ctx}: field {field!r} is not a number."
        raise IngestionError(msg) from exc
    if value < 0:
        msg = f"{ctx}: field {field!r} is negative."
        raise IngestionError(msg)
    return value


def _build_entry(row: Sequence[str], *, index: dict[str, int], ctx: str) -> CatalogEntry:
    if len(row) > len(index):
        msg = f"{ctx}: {len(row)} fields for a {len(index)}-column header."
        raise IngestionError(msg)

    def get(name: str) -> str:
        return _cell(row, index=index, name=name)

    drug_id = _require(get("drug_id"), field="drug_id", ctx=ctx)
    ctx = f"{ctx} ({drug_id})"
    return CatalogEntry(
        drug_id=drug_id,
        name=_require(get("name"), field="name", ctx=ctx),
        approved_name=get("approved_name") or None,
        approval_status=_parse_enum(
            ApprovalStatus,
            _require(get("approval_status"), field="approval_status", ctx=ctx),
            field="approval_status",
            ctx=ctx,
        ),
        intervention_class=_parse_enum(
            InterventionClass,
            _require(get("intervention_class"), field="intervention_class", ctx=ctx),
            field="intervention_class",
            ctx=ctx,
        ),
        target=_require(get("target"), field="target", ctx=ctx),
        target_node_id=_require(get("target_node_id"), field="target_node_id", ctx=ctx),
        mechanism_of_action=_require(
            get("mechanism_of_action"), field="mechanism_of_action", ctx=ctx
        ),
        observed_direction=_parse_enum(
            EffectDirection,
            _require(get("observed_direction"), field="observed_direction", ctx=ctx),
            field="observed_direction",
            ctx=ctx,
        ),
        is_direct_evidence=_parse_bool(
            _require(get("is_direct_evidence"), field="is_direct_evidence", ctx=ctx),
            field="is_direct_evidence",
            ctx=ctx,
        ),
        strongest_evidence_type=_parse_enum(
            EvidenceType,
            _require(get("strongest_evidence_type"), field="strongest_evidence_type", ctx=ctx),
            field="strongest_evidence_type",
            ctx=ctx,
        ),
        has_pediatric_exposure=_parse_bool(
            _require(get("has_pediatric_exposure"), field="has_pediatric_exposure", ctx=ctx),
            field="has_pediatric_exposure",
            ctx=ctx,
        ),
        youngest_age_studied=get("youngest_age_studied") or None,
        pediatric_indication=get("pediatric_indication") or None,
        route=get("route") or None,
        cns_penetrant=_parse_optional_bool(get("cns_penetrant"), field="cns_penetrant", ctx=ctx),
        achievable_plasma_um=_parse_optional_float(
            get("achievable_plasma_um"), field="achievable_plasma_um", ctx=ctx
        ),
        required_effective_um=_parse_optional_float(
            get("required_effective_um"), field="required_effective_um", ctx=ctx
        ),
        half_life_hours=_parse_optional_float(
            get("half_life_hours"), field="half_life_hours", ctx=ctx
        ),
        worsens_cin=_parse_optional_bool(get("worsens_cin"), field="worsens_cin", ctx=ctx),
        validation_experiment=_require(
            get("validation_experiment"), field="validation_experiment", ctx=ctx
        ),
        notes=get("notes"),
    )
