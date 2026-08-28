"""Persistence layer: the evidence ledger and the DuckDB/Parquet evidence store.

Layer 6 in the GP-01 stack. Everything upstream produces typed models; this
package is where they become durable, queryable and exportable, and where GP-10
("no claim without evidence") stops being a convention and becomes a gate.

Two halves:

* :class:`EvidenceLedger` / :class:`AssertionResolver` — the run's write buffer
  and the GP-10 gate. In memory by default; a ledger handed a ``spill_dir``
  overflows into :class:`SqliteEvidenceSpill` so that a whole-genome run's tens
  of millions of items do not have to be resident at once.
* :class:`EvidenceStore` — the DuckDB file, its schema, idempotent writers for
  every domain model, and a byte-reproducible Parquet export (GP-30).

The graph is stored relationally in ``graph_edges``; this project deliberately
does not adopt a graph database. See ``schema.sql`` for the reasoning.
"""

from __future__ import annotations

from mva.evidence.ledger import (
    DEFAULT_SPILL_THRESHOLD,
    AssertionResolver,
    EvidenceLedger,
)
from mva.evidence.spill import DEFAULT_FLUSH_BATCH, SqliteEvidenceSpill
from mva.evidence.store import (
    EVIDENCE_WRITE_BATCH,
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
    "DEFAULT_FLUSH_BATCH",
    "DEFAULT_SPILL_THRESHOLD",
    "EVIDENCE_WRITE_BATCH",
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
    "SqliteEvidenceSpill",
]
