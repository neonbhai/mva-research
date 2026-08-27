# ADR 0007 — Synthetic-first development

**Status:** accepted · **Date:** 2026-08-27

## Context
The pipeline must be built, debugged and demonstrated before — and arguably
without ever — touching the real case. Debugging genomics code involves printing
records, and printing a real record puts it in a terminal, a log, a test failure
diff, and an agent's context.

## Decision
Develop entirely against a **fully synthetic case** with a known correct answer.
Every test, the demo, and all golden expectations use it. The real case is
connected only at the end, by a human, with the network profile enforced.

## The fixture is designed adversarially
It is not a happy path. It contains the specific traps this pipeline must not
fall into:
- two rare damaging heterozygous variants in one fictional gene (the answer);
- common distractors in the same gene, which a naive ranker pairs enthusiastically;
- a low-quality call that must be flagged, not deleted;
- an **in-cis** pair with phase-set evidence, which must not be called a compound
  heterozygote;
- a multiallelic site requiring correct splitting and per-allele depth assignment;
- a phenotype profile with an explicitly excluded term, a not-assessed term and an
  uncertain term, so the four-valued logic is exercised;
- a drug catalogue containing a right-target/wrong-direction agent, a tool
  compound, a symptomatic agent, a context-dependent agent, and an off-mechanism
  agent.

If the pipeline gets the synthetic case right, it is right for structural
reasons, not because it was tuned to a leaderboard.

## Consequences
- Fictional gene symbols (`SYNTH*`) make it impossible to mistake demo output for
  a real finding.
- The synthetic answer must never be special-cased. A test asserting the correct
  pair ranks first is only meaningful if the general scoring produces it.
- Golden files are locked by GP-32 and are not re-baselined to make a test pass.
