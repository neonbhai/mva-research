"""Load curated mechanism chains from the local knowledge tables.

The chain table (`knowledge/public/mechanisms.tsv`) is a **two-record-type table
sharing one header**: a row either populates the *node block*
(``node_id .. deviation_is_pathological``) or the *link block*
(``link_id .. uncertainty``), never both. Rows are therefore matched by block
rather than by absolute column index, which keeps the loader robust against the
commonest authoring artefact in a hand-maintained TSV — one tab too many or too
few in the padding between the two blocks. A row that populates both blocks, or
neither, is an error rather than a guess.

Everything here is a *loader*, not an inference engine: it turns a curated table
into `MechanismHypothesis` values (GP-02) and refuses malformed input loudly. No
scoring, no evidence emission, no drug knowledge — the chain is passed into the
intervention stage as an argument (GP-03).

GP-20: the tables shipped in this repository are **synthetic**. The library
carries a `version` string that is stamped onto every downstream evidence item so
a reader can always tell which table produced a claim.
"""

from __future__ import annotations

import csv
from collections.abc import Sequence
from enum import StrEnum
from pathlib import Path

from mva.errors import IngestionError
from mva.models.base import AssertionTier
from mva.models.evidence import EvidenceStrength
from mva.models.mechanism import (
    EffectDirection,
    MechanismHypothesis,
    MechanismLink,
    MechanismNode,
    MechanismNodeKind,
)

__all__ = ["MechanismLibrary"]

#: Key columns present on every row of the chain table.
_KEY_FIELDS: tuple[str, ...] = ("mechanism_id", "gene_symbol")

#: The node block of the chain table, in file order.
_NODE_FIELDS: tuple[str, ...] = (
    "node_id",
    "node_kind",
    "node_label",
    "node_identifier",
    "state_in_patient",
    "deviation_is_pathological",
)

#: Node columns that must carry an explicit value on every node row.
#: ``node_identifier`` is the only optional one: an external ID is genuinely
#: absent for some nodes. ``deviation_is_pathological`` is required precisely
#: because its convenient default -- "everything that deviates is disease" --
#: inverts the corrective sign on a compensatory node (see MechanismNode).
_REQUIRED_NODE_FIELDS: tuple[str, ...] = tuple(
    field for field in _NODE_FIELDS if field != "node_identifier"
)

#: The link block of the chain table, in file order.
_LINK_FIELDS: tuple[str, ...] = (
    "link_id",
    "source_node_id",
    "target_node_id",
    "relation",
    "direction",
    "tier",
    "strength",
    "is_directly_demonstrated",
    "uncertainty",
)

#: Columns of the mechanism metadata table.
_META_FIELDS: tuple[str, ...] = (
    "mechanism_id",
    "gene_symbol",
    "summary",
    "disease_direction",
    "therapeutic_target_node_id",
    "required_correction",
    "developmental_window_caveat",
)

_TRUE_TOKENS: frozenset[str] = frozenset({"1", "true", "yes", "y", "t"})
_FALSE_TOKENS: frozenset[str] = frozenset({"0", "false", "no", "n", "f"})

#: Stand-in text used when a curator left a link's `uncertainty` cell empty.
#: An empty uncertainty is a gap in the record, not a statement of confidence
#: (GP-17), so it is rendered as such rather than as an empty string.
_UNSTATED_UNCERTAINTY = "Not stated by the curator of this link; treat as unquantified."


def _read_table(path: Path, *, expected_header: Sequence[str]) -> tuple[tuple[str, ...], ...]:
    """Read a `#`-commented TSV into stripped rows, header excluded.

    Cells are stripped, so a cell containing only whitespace is indistinguishable
    from an empty one — which is what the source tables mean by it.
    """
    if not path.is_file():
        msg = f"Mechanism table not found: {path.name}"
        raise IngestionError(msg)
    lines = [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not lines:
        msg = f"Mechanism table {path.name} contains no header row."
        raise IngestionError(msg)
    rows = [tuple(cell.strip() for cell in row) for row in csv.reader(lines, delimiter="\t")]
    header = rows[0]
    missing = [name for name in expected_header if name not in header]
    if missing:
        msg = (
            f"Mechanism table {path.name} is missing {len(missing)} required "
            f"column(s): {', '.join(missing)}. Add the column and give every row an "
            "explicit value. A missing column is never filled with a default here: the "
            "fields this loader requires are exactly the ones whose wrong-by-default "
            "value would invert a downstream sign check."
        )
        raise IngestionError(msg)
    return tuple(rows[1:])


def _pad(values: Sequence[str], width: int) -> tuple[str, ...]:
    return tuple(values) + ("",) * max(0, width - len(values))


def _resolve_row(row: Sequence[str], *, source: str, line_no: int) -> tuple[str, dict[str, str]]:
    """Classify one chain row as a node row or a link row and map its block.

    Returns ``("node" | "link", mapping)``. Raises `IngestionError` when the row
    is ambiguous, which is the only safe response: silently guessing which block
    a half-filled row belongs to would fabricate a mechanism edge.
    """
    if len(row) <= len(_KEY_FIELDS):
        msg = f"{source} line {line_no}: row has {len(row)} fields; expected at least 3."
        raise IngestionError(msg)
    keys = dict(zip(_KEY_FIELDS, row[: len(_KEY_FIELDS)], strict=False))
    if not keys["mechanism_id"] or not keys["gene_symbol"]:
        msg = f"{source} line {line_no}: row is missing mechanism_id or gene_symbol."
        raise IngestionError(msg)

    tail = list(row[len(_KEY_FIELDS) :])
    offset = 0
    while offset < len(tail) and tail[offset] == "":
        offset += 1
    if offset >= len(tail):
        msg = f"{source} line {line_no}: row populates neither the node block nor the link block."
        raise IngestionError(msg)

    if offset == 0:
        block = _pad(tail[: len(_NODE_FIELDS)], len(_NODE_FIELDS))
        if any(cell for cell in tail[len(_NODE_FIELDS) :]):
            msg = f"{source} line {line_no}: row populates both the node and the link block."
            raise IngestionError(msg)
        mapping = {**keys, **dict(zip(_NODE_FIELDS, block, strict=True))}
        for field in _REQUIRED_NODE_FIELDS:
            if not mapping[field]:
                msg = (
                    f"{source} line {line_no}: node row is missing required field {field!r}. "
                    "It has no default; supply the value."
                )
                raise IngestionError(msg)
        return "node", mapping

    block = _pad(tail[offset : offset + len(_LINK_FIELDS)], len(_LINK_FIELDS))
    if any(cell for cell in tail[offset + len(_LINK_FIELDS) :]):
        msg = f"{source} line {line_no}: link row carries {len(row)} fields; trailing data ignored."
        raise IngestionError(msg)
    mapping = {**keys, **dict(zip(_LINK_FIELDS, block, strict=True))}
    for field in _LINK_FIELDS[:-1]:
        if not mapping[field]:
            msg = f"{source} line {line_no}: link row is missing required field {field!r}."
            raise IngestionError(msg)
    return "link", mapping


def _parse_enum[E: StrEnum](enum_cls: type[E], raw: str, *, field: str, ctx: str) -> E:
    """Parse an enum cell.

    The offending value is deliberately NOT echoed into the exception message
    (PRIV-09): messages carry identifiers, field names and the allowed set only.
    """
    try:
        return enum_cls(raw)
    except ValueError as exc:
        allowed = ", ".join(sorted(str(member.value) for member in enum_cls))
        msg = f"{ctx}: field {field!r} holds an unrecognised value. Allowed: {allowed}."
        raise IngestionError(msg) from exc


def _parse_bool(raw: str, *, field: str, ctx: str) -> bool:
    lowered = raw.strip().lower()
    if lowered in _TRUE_TOKENS:
        return True
    if lowered in _FALSE_TOKENS:
        return False
    msg = f"{ctx}: field {field!r} is not a recognised boolean. Use 1/0."
    raise IngestionError(msg)


class MechanismLibrary:
    """An immutable, versioned collection of curated mechanism chains.

    Lookup is by gene symbol because that is the key the upstream candidate-pair
    stage carries. A gene with no curated chain returns ``None`` rather than an
    empty mechanism: "we have no chain" and "the chain is empty" are different
    facts, and only the first is honest (GP-14).
    """

    def __init__(self, mechanisms: Sequence[MechanismHypothesis], *, version: str) -> None:
        if not version.strip():
            msg = "MechanismLibrary requires a non-empty version string (GP-20, GP-31)."
            raise IngestionError(msg)
        self._version = version
        self._mechanisms: tuple[MechanismHypothesis, ...] = tuple(
            sorted(mechanisms, key=lambda m: m.mechanism_id)
        )
        by_gene: dict[str, list[MechanismHypothesis]] = {}
        for mechanism in self._mechanisms:
            by_gene.setdefault(mechanism.gene_symbol.upper(), []).append(mechanism)
        self._by_gene: dict[str, tuple[MechanismHypothesis, ...]] = {
            gene: tuple(items) for gene, items in by_gene.items()
        }

    # ------------------------------------------------------------------ loading

    @classmethod
    def from_tsv(cls, chain_path: Path, meta_path: Path, *, version: str) -> MechanismLibrary:
        """Build a library from the chain table and the mechanism metadata table.

        Both tables are required: a chain without metadata has no therapeutic
        target and no required correction, and a mechanism the intervention stage
        cannot ask "which way must this move?" is not usable (GP-16).
        """
        chain_rows = _read_table(chain_path, expected_header=(*_KEY_FIELDS, *_NODE_FIELDS))
        meta_rows = _read_table(meta_path, expected_header=_META_FIELDS)

        nodes: dict[str, list[MechanismNode]] = {}
        links: dict[str, list[MechanismLink]] = {}
        genes: dict[str, str] = {}
        for index, row in enumerate(chain_rows, start=1):
            kind, mapping = _resolve_row(row, source=chain_path.name, line_no=index)
            mechanism_id = mapping["mechanism_id"]
            ctx = f"{chain_path.name} row {index} (mechanism {mechanism_id})"
            genes.setdefault(mechanism_id, mapping["gene_symbol"])
            if kind == "node":
                nodes.setdefault(mechanism_id, []).append(_build_node(mapping, ctx=ctx))
            else:
                links.setdefault(mechanism_id, []).append(_build_link(mapping, ctx=ctx))

        meta_by_id = _parse_meta(meta_rows, source=meta_path.name)
        built: list[MechanismHypothesis] = []
        for mechanism_id in sorted(nodes | links):
            built.append(
                _assemble(
                    mechanism_id,
                    gene_symbol=genes[mechanism_id],
                    nodes=tuple(nodes.get(mechanism_id, ())),
                    links=tuple(links.get(mechanism_id, ())),
                    meta=meta_by_id.get(mechanism_id),
                    source=meta_path.name,
                )
            )
        return cls(built, version=version)

    # ------------------------------------------------------------------ queries

    def for_gene(self, gene_symbol: str) -> MechanismHypothesis | None:
        """The curated chain for a gene, or ``None`` when none is curated.

        Matching is case-insensitive. If a gene somehow carries several chains the
        lowest `mechanism_id` wins, so the answer never depends on file order.
        """
        found = self._by_gene.get(gene_symbol.strip().upper())
        return found[0] if found else None

    def all_mechanisms(self) -> tuple[MechanismHypothesis, ...]:
        """Every chain, ordered by `mechanism_id` (GP-30)."""
        return self._mechanisms

    @property
    def version(self) -> str:
        """Version of the knowledge tables behind this library."""
        return self._version

    def __len__(self) -> int:
        return len(self._mechanisms)


def _build_node(mapping: dict[str, str], *, ctx: str) -> MechanismNode:
    return MechanismNode(
        node_id=mapping["node_id"],
        kind=_parse_enum(MechanismNodeKind, mapping["node_kind"], field="node_kind", ctx=ctx),
        label=mapping["node_label"],
        identifier=mapping["node_identifier"] or None,
        state_in_patient=_parse_enum(
            EffectDirection, mapping["state_in_patient"], field="state_in_patient", ctx=ctx
        ),
        deviation_is_pathological=_parse_bool(
            mapping["deviation_is_pathological"], field="deviation_is_pathological", ctx=ctx
        ),
    )


def _build_link(mapping: dict[str, str], *, ctx: str) -> MechanismLink:
    return MechanismLink(
        link_id=mapping["link_id"],
        source_node_id=mapping["source_node_id"],
        target_node_id=mapping["target_node_id"],
        relation=mapping["relation"],
        direction=_parse_enum(EffectDirection, mapping["direction"], field="direction", ctx=ctx),
        tier=_parse_enum(AssertionTier, mapping["tier"], field="tier", ctx=ctx),
        strength=_parse_enum(EvidenceStrength, mapping["strength"], field="strength", ctx=ctx),
        is_directly_demonstrated=_parse_bool(
            mapping["is_directly_demonstrated"], field="is_directly_demonstrated", ctx=ctx
        ),
        uncertainty=mapping["uncertainty"] or _UNSTATED_UNCERTAINTY,
    )


def _parse_meta(rows: Sequence[Sequence[str]], *, source: str) -> dict[str, dict[str, str]]:
    parsed: dict[str, dict[str, str]] = {}
    for index, row in enumerate(rows, start=1):
        values = _pad(row, len(_META_FIELDS))
        if len(values) > len(_META_FIELDS):
            msg = f"{source} row {index}: {len(values)} fields; expected {len(_META_FIELDS)}."
            raise IngestionError(msg)
        mapping = dict(zip(_META_FIELDS, values, strict=True))
        mechanism_id = mapping["mechanism_id"]
        if not mechanism_id:
            msg = f"{source} row {index}: missing mechanism_id."
            raise IngestionError(msg)
        if mechanism_id in parsed:
            msg = f"{source} row {index}: duplicate metadata for mechanism {mechanism_id}."
            raise IngestionError(msg)
        parsed[mechanism_id] = mapping
    return parsed


def _assemble(
    mechanism_id: str,
    *,
    gene_symbol: str,
    nodes: tuple[MechanismNode, ...],
    links: tuple[MechanismLink, ...],
    meta: dict[str, str] | None,
    source: str,
) -> MechanismHypothesis:
    """Combine one mechanism's nodes, links and metadata, validating referential integrity."""
    ctx = f"mechanism {mechanism_id}"
    if meta is None:
        msg = f"{ctx}: chain rows exist but {source} carries no metadata row for it."
        raise IngestionError(msg)
    if not nodes:
        msg = f"{ctx}: no node rows."
        raise IngestionError(msg)
    node_ids = {node.node_id for node in nodes}
    if len(node_ids) != len(nodes):
        msg = f"{ctx}: duplicate node_id among {len(nodes)} node rows."
        raise IngestionError(msg)
    link_ids = {link.link_id for link in links}
    if len(link_ids) != len(links):
        msg = f"{ctx}: duplicate link_id among {len(links)} link rows."
        raise IngestionError(msg)
    for link in links:
        for endpoint in (link.source_node_id, link.target_node_id):
            if endpoint not in node_ids:
                msg = (
                    f"{ctx}: link {link.link_id} references node {endpoint}, which is not declared."
                )
                raise IngestionError(msg)

    target = meta["therapeutic_target_node_id"]
    if target not in node_ids:
        msg = f"{ctx}: therapeutic target {target} is not among the {len(nodes)} declared nodes."
        raise IngestionError(msg)

    uncertainties = tuple(
        f"{link.link_id}: {link.uncertainty}" for link in links if link.uncertainty
    )
    return MechanismHypothesis(
        mechanism_id=mechanism_id,
        gene_symbol=gene_symbol,
        pair_id=None,
        summary=meta["summary"],
        nodes=nodes,
        links=links,
        disease_direction=_parse_enum(
            EffectDirection, meta["disease_direction"], field="disease_direction", ctx=ctx
        ),
        therapeutic_target_node_id=target,
        required_correction=_parse_enum(
            EffectDirection, meta["required_correction"], field="required_correction", ctx=ctx
        ),
        uncertainties=uncertainties,
        developmental_window_caveat=meta["developmental_window_caveat"],
    )
