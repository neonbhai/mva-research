# ADR 0003 — Deterministic local pipeline; no LLM in the patient-data path

**Status:** accepted · **Date:** 2026-08-27

## Context
An LLM could plausibly do parts of this work at runtime: reading the VCF,
judging variant plausibility, drafting the mechanism, ranking drugs. The tooling
to do so is readily available and would have been faster to build.

## Decision
The patient-data path contains **no model inference**. Every stage from VCF
ingestion through candidate ranking is deterministic Python. Claude is the
engineer that *built* the machinery, never a component that *runs* inside it.

## Rationale
1. **Privacy.** Sending a child's genotypes to a hosted model is a disclosure to
   a third party. Even with a local model, prompt and output logging becomes a
   new copy of patient data in a new place. The cleanest guarantee is that the
   data never enters a model at all.
2. **Reproducibility.** The challenge may rerun our submission. Sampling
   temperature, model version drift and prompt sensitivity make that
   irreproducible in a way that no manifest can capture.
3. **Auditability.** "Why is this pair ranked first?" must be answerable by
   reading code and component scores. A model's ranking is a judgement we cannot
   decompose, and GP-10 requires every claim to resolve to structured evidence.
4. **Eloquence is not evidence.** A fluent mechanistic narrative is exactly the
   failure mode this project is built to resist. Model fluency correlates with
   persuasiveness, not with correctness.

## Consequences
- Heuristics must be written explicitly, which makes them criticisable. That is
  the point.
- We forgo the flexibility of natural-language reasoning over literature. The
  mechanism and drug knowledge instead lives in versioned, citable tables.
- No hosted-LLM dependency is required to run anything in this repo.

## What an LLM may still do
Build-time engineering, documentation, adversarial review of *code and
assumptions*, and — with explicit human approval and public data only —
literature triage that produces citations for the knowledge tables. Never
inference over patient records.
