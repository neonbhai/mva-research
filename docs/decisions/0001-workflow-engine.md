# ADR 0001 — Snakemake for the reproducible workflow

**Status:** accepted · **Date:** 2026-08-27

## Context
The pipeline is a DAG of file-producing stages over patient data, and the
challenge states organizers may rerun submissions. We need a workflow layer that
makes the DAG explicit, resumable and reproducible. Candidates: Snakemake,
Nextflow, a plain Makefile, Prefect/Airflow, or a hand-rolled Python driver.

## Decision
Use **Snakemake** as the workflow engine, with the DAG also executable through
`mva run` so that Snakemake is a convenience, not a hard dependency.

## Rationale
- Snakemake is the dominant convention in academic genomics; a reviewer reading
  this repo already knows how to read a `Snakefile`.
- It is a Python library, so it installs cleanly on macOS/arm64 via `uv`
  (verified) and does not need a JVM as Nextflow does.
- File-based DAG semantics match our provenance model: every stage's output is a
  hashed artifact.
- Airflow/Prefect are orchestration servers. We want a local, offline, single-
  machine run with no daemon and no network listener — a scheduler service would
  be a liability under the privacy model, not an asset.
- A plain Makefile cannot express the wildcard/config-driven structure without
  becoming unreadable.

## Consequences
- Every stage must be a pure function of (inputs, config) to files. That is a
  constraint we wanted anyway for GP-30.
- The rules stay thin: they shell out to `mva` subcommands, so the logic lives in
  tested Python rather than in the workflow file.
- Because the same stages are callable via `mva run all`, a machine without
  Snakemake can still run the pipeline. Snakemake is an optional extra.
