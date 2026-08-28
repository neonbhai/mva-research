# Architecture

## The one pipeline

```
VCF ─▶ ingest ─▶ normalise ─▶ QC ─▶ annotate ─▶ select ─▶ pair ─▶ score ─▶ rank
                                                                            │
 phenotype ─▶ profile ─▶ gene-phenotype match ──────────────────────────────┤
                                                                            ▼
                                                          ┌─── Track 1: submission + dossier
                                                          │
                                            ranked pair ──┤
                                                          │
                                                          └──▶ mechanism ─▶ intervention
                                                                    │            │
                                                                    │       direction check
                                                                    │            │
                                                                    │       safety filter
                                                                    │            │
                                                                    └────▶ Track 2 report
```

One pipeline, two exits. Track 1 leaves after ranking; Track 2 continues from a
selected or highly-ranked pair. There is no separate Track 2 codebase, because a
drug hypothesis that is not anchored to the ranked variant pair is not grounded
in anything.

Every arrow writes a hashed artifact with provenance. Every claim along the way
becomes an `EvidenceItem`.

## Layering (GP-01)

```
8  cli, pipeline              composition root — the only layer that knows everything
7  reporting                  rendering, assertion checking
6  evidence                   DuckDB + Parquet persistence
5  prioritization, mechanisms, interventions      analysis
4  ingestion, annotation, phenotype               data in
3  privacy                    cross-cutting; imports layers 0-2 only
2  config                     typed settings + the workspace boundary
1  errors, determinism, clock   foundation utilities
0  models                     the shared vocabulary; imports nothing
```

Imports go down, never up. Peer stages never import each other — they compose
only in `src/mva/pipeline.py` (GP-03). Both rules are enforced by AST-walking
tests in `tests/unit/test_architecture.py`, whose failure messages carry their own
remediation, because the reader is often an agent whose only view of the rule is
that string.

This layering was not free: the first run of the layering test caught `config`
importing `determinism`, which is why foundation utilities sit *below* config
rather than above it.

## The three contracts everything speaks

**`VariantRecord`** — a build-anchored coordinate, a genotype with its analytical
evidence, and *lists* of transcript consequences and population frequencies.
Immutable and additive: annotation returns a new record, so ingested and
annotated are separately hashable artifacts.

**`EvidenceItem`** — one structured, attributable, falsifiable claim, with
mandatory `limitations`. Content-derived IDs mean the same conclusion drawn twice
is one piece of evidence, so the ID is stable across runs. The *record* is not
byte-identical across runs under a real clock: `timestamp` is excluded from the ID
but is still written to `evidence_items.parquet` (`docs/handoff-integrity.md` §4).

**`CandidatePair`** — the Track 1 unit of prediction. Its defining choice is that
the score is a **vector, not a scalar**: eight component scores stay visible into
the report, so "why is this ranked first?" is answerable without re-running
anything.

## Adapter boundaries

Every external tool or data source sits behind a `Protocol`, with a local,
hash-pinned implementation today and a documented slot for the real thing:

| Boundary | Today | Future |
|---|---|---|
| VCF reading | `cyvcf2` backend + pure-Python fallback | — |
| Normalisation | in-process trim/left-align | `bcftools norm` |
| Consequence | local TSV | VEP / SnpEff |
| Frequency | local TSV | gnomAD sites-only |
| Clinical | null adapter (returns nothing, honestly) | ClinVar release VCF |
| Splice | column on the consequence table | SpliceAI |
| Phenotype ranking | local gene→HPO index | Exomiser / LIRICAL |
| SV / CNV | not implemented | dedicated callers |
| Drugs / literature | local TSV | DrugBank / ChEMBL / PubMed |

`docs/maturity-ledger.md` grades each of these `real`, `synthetic-substitute` or
`stub`. A synthetic substitute is never described as biologically valid (GP-20).

**A real remote adapter may never live on the sensitive path.** A structural test
forbids importing a network client anywhere under `ingestion`, `annotation`,
`phenotype` or `prioritization`. Reference data is acquired in a separate,
public-only step; proband coordinates are never transmitted.

## Storage

DuckDB for querying, Parquet for large typed artifacts, no server (ADR 0002).
Graph structure — variant → gene → mechanism → drug — is stored **relationally**
in a `graph_edges` subject/predicate/object table. The domain is a graph, but the
queries are shallow and the data is small; keeping edges relational avoids a
second system to secure and delete, while leaving the export path open.

## Determinism

Byte-identical repeat runs are an acceptance criterion, not an aspiration
(GP-30) — **under a fixed clock**, which is what `config.synthetic` selects and
what every determinism check in this repo runs under. For a real case the clock
is `SystemClock` (`pipeline.py:450`), so the scientific content is identical and
eleven artifacts differ in recorded time only; the artifact-by-artifact map is in
`docs/handoff-integrity.md` §4 and the fix is TD-21. The mechanisms:

- all timestamps come from an injected `Clock`; a lint forbids `datetime.now()`;
- evidence and pair IDs are content-derived;
- every sort key is total (score, then genomic position);
- Parquet writes are sorted by primary key with pinned compression and row-group
  size, and the *writer* embeds no timestamp of its own — though
  `evidence_items` carries the evidence timestamp as data;
- canonical JSON for every hash.

`just demo-determinism` runs the **demo** twice and compares artifact hashes,
skipping `provenance.json` and the DuckDB container by filename. Nothing has
ever run it against a non-synthetic config, and pointed at one it would report
a failure on the eleven time-bearing artifacts.

## Privacy in the architecture

The workspace is external and validated before any file is read (ADR 0006). The
`privacy` package sits at layer 3 so any stage can call it while it stays unable
to import upward. Public export is gated twice — an allowlist plus a post-render
content re-scan — because classification is a claim and the scan is the
verification (GP-43).

The deliberate asymmetry: this repo is built to be maximally legible to an agent,
*except* for the patient-data path, which is built to be illegible to one
(ADR 0008).

## Workflow

Snakemake expresses the DAG (ADR 0001), but the rules are thin — they shell out
to `mva` subcommands so the logic stays in tested Python. The same DAG runs
without Snakemake via `mva run all`, so the workflow engine is a convenience
rather than a dependency.
