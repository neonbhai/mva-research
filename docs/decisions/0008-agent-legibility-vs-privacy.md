# ADR 0008 — Agent legibility stops at the patient-data boundary

**Status:** accepted · **Date:** 2026-08-27

## Context
Current practice for agent-built codebases optimises for **legibility to the
agent**: push all knowledge into the repo, because anything the agent cannot
access in context effectively does not exist. Structured docs, checked-in plans,
custom lints whose error messages carry their own remediation. We have adopted
all of that.

It is in direct tension with our privacy model, which requires part of the system
to be **deliberately illegible** to the agent. The general advice assumes the
agent should be able to see everything it works on. Here, one category of thing
must never enter its context at all.

## Decision
Split the two.

**Legible to agents — maximise:** architecture, contracts, golden principles,
decision records, the synthetic case, test fixtures, knowledge tables, lint
remediation messages, maturity grades. All in-repo, all mechanically verified.

**Illegible to agents — enforce:** anything under `MVA_WORKSPACE`. Patient VCFs,
phenotype extracts, reads, derived variant tables, internal reports, logs.

Concretely:
- Introspection commands (`mva inspect`, QC summaries, audit output) operate on
  **synthetic runs and aggregate counts only** — never records.
- Error messages, warnings and skip reasons carry IDs, counts and reason codes,
  never record text.
- The privacy scanner reports paths, line numbers and span lengths, never matched
  content (GP-41).
- Debugging is done on the synthetic case. Always.

## Consequences
- A genuine cost: an agent cannot debug a real-data failure by looking at it. The
  mitigation is that the synthetic fixture is built to reproduce the structural
  shapes that fail, and that a human can inspect the real data directly.
- Some conveniences are deliberately unavailable, e.g. printing a dataframe head.
  `safe_repr()` gives shape and column names instead.
- This asymmetry is the design, not an oversight. Record it when someone asks why
  a debugging affordance is missing.
