# ADR 0002 — DuckDB + Parquet instead of a conventional database

**Status:** accepted · **Date:** 2026-08-27

## Context
We need to store variants, annotations, pairs, phenotypes, evidence, mechanisms,
drugs, citations, runs, provenance and graph edges, and to query across them
analytically. Options: PostgreSQL/SQLite, DuckDB+Parquet, a document store, or a
graph database.

## Decision
**DuckDB** as the analytical query engine, **Parquet** for large typed artifacts.
No server. Graph structure is stored **relationally** in a `graph_edges` table.

## Rationale
- Zero-configuration and embedded. A server process holding patient genotypes,
  listening on a port and writing its own WAL somewhere outside our control is a
  privacy liability. DuckDB is a file we can place inside the encrypted
  workspace and delete atomically.
- Columnar and vectorised: the natural shape of variant tables. Queries that
  would need hand-written indices in SQLite are fast by default.
- Parquet is typed, compressed, portable and *hashable* — which is what lets the
  determinism test compare artifacts byte for byte.
- Postgres would add an install, a daemon, a network surface and a backup path we
  would then have to reason about deleting.

## Why not a graph database
The domain is genuinely a graph (variant → gene → mechanism → drug), and Neo4j
would model it directly. We are not adopting one because:
- it reintroduces a server and a second storage system to secure and delete;
- the queries we actually run are shallow (1–3 hops), which SQL joins handle;
- the graph is small — thousands of edges, not billions.

We keep the option open by storing subject–predicate–object edges relationally
with an export path. If traversal depth grows, the edges are already in the right
shape to load into a graph engine.

## Consequences
- Byte-identical Parquet requires care: sort by primary key, pin compression and
  row-group size, and let the *writer* embed no timestamp of its own. That much is
  asserted in a test. It is not a claim about the data: `evidence_items.parquet`
  carries `timestamp`/`timestamp_iso` columns and three more tables carry a
  rendered timestamp inside a rationale string, so those four differ between real
  runs — in recorded time only, never in content (`docs/handoff-integrity.md` §4).
- Analysts get SQL over the evidence store without a migration or a server.
