# ADR 0009 — Determinism and privacy checks are blocking gates

**Status:** accepted · **Date:** 2026-08-27

## Context
High-throughput agent-built repos often adopt minimal blocking merge gates:
short-lived PRs, flaky tests retried rather than blocking, auto-merge on green.
The reasoning is sound where throughput is the constraint and a flake is noise.

## Decision
**Reject that trade-off here.** `just verify` is a single blocking gate covering
lint, typecheck, unit/integration/golden tests, the architecture lints, the docs
integrity check and the privacy audit. Nothing merges past a failure.

## Rationale
- **A flaky test in this repo is a determinism bug, not noise.** GP-30 makes
  byte-identical repeat runs an acceptance criterion. A test that passes on retry
  has just told us the pipeline is non-deterministic — precisely the signal a
  retry policy would discard.
- **Privacy failures are not recoverable after the fact.** A disclosed genotype
  cannot be un-disclosed, and git history preserves what `rm` removes. There is
  no throughput argument that outweighs this.
- Throughput is not our constraint. This is a small codebase over a short
  horizon, not a million-line product at 1,500 PRs.

## Consequences
- Slower iteration when a check fails. Accepted.
- The gate must stay fast enough to run constantly (currently seconds), or it
  will be bypassed — that is the real risk, and it is a reason to keep the test
  suite lean rather than to weaken the gate.
- Review findings are promoted into tests or lints so the gate keeps enforcing
  them after the reviewer has gone.
