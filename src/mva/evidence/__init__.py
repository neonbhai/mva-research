"""Persistence layer: the evidence ledger and the DuckDB/Parquet evidence store.

Layer 6 in the GP-01 stack. Everything upstream produces typed models; this
package is where they become durable, queryable and exportable, and where GP-10
("no claim without evidence") stops being a convention and becomes a gate.

Two halves:

* :class:`EvidenceLedger` / :class:`AssertionResolver` — in-memory, no I/O. Stages
  accumulate evidence as they work; reporting resolves citations through the
  resolver, which refuses unsourced claims.
* :class:`EvidenceStore` — the DuckDB file, its schema, idempotent writers for
  every domain model, and a byte-reproducible Parquet export (GP-30).

The graph is stored relationally in ``graph_edges``; this project deliberately
does not adopt a graph database. See ``schema.sql`` for the reasoning.
"""

from __future__ import annotations

from mva.evidence.ledger import AssertionResolver, EvidenceLedger
from mva.evidence.store import (
    PARQUET_COMPRESSION,
    PARQUET_COMPRESSION_LEVEL,
    PARQUET_ROW_GROUP_SIZE,
    PARQUET_VERSION,
    PRIMARY_KEYS,
    TABLES,
    EvidenceStore,
    GraphEdge,
)

__all__ = [
    "PARQUET_COMPRESSION",
    "PARQUET_COMPRESSION_LEVEL",
    "PARQUET_ROW_GROUP_SIZE",
    "PARQUET_VERSION",
    "PRIMARY_KEYS",
    "TABLES",
    "AssertionResolver",
    "EvidenceLedger",
    "EvidenceStore",
    "GraphEdge",
]
