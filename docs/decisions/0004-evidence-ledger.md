# ADR 0004 — The evidence ledger

**Status:** accepted · **Date:** 2026-08-27

## Context
Rare-disease analysis produces confident-sounding conclusions from stacked weak
inferences. The usual failure is not a wrong number but an unattributable one:
a report states "this variant is rare and damaging" and nothing records which
database at which version said so, or what that claim does not establish.

## Decision
Every claim is an **`EvidenceItem`**: a structured, append-only record with
subject, claim, category, direction (supports/contradicts/neutral), strength,
evidence type, assertion tier, citation, method, tool, tool version,
**limitations**, and timestamp. Nothing is stated in a report unless it resolves
to evidence IDs (`AssertionResolver.require`).

## Notable design choices
- **`limitations` is mandatory.** An in-silico score with no stated limitation is
  how prediction becomes mistaken for proof.
- **`direction` includes `contradicts`.** Contradicting evidence is stored
  alongside supporting evidence and surfaced in reports (GP-19). Discarding it
  would be the most consequential possible omission.
- **`tier` separates epistemics from modality.** `AssertionTier` says how we came
  to believe it (observed / database / prediction / literature / inference /
  speculation); `EvidenceType` says what kind of experiment it was. A validator
  forbids labelling a prediction as an observation.
- **IDs are content-derived**, excluding the timestamp: the same conclusion drawn
  twice is one piece of evidence, not two. So the *ID* is stable across runs, which
  is what makes the ledger idempotent. It does not make repeat runs byte-identical:
  `EvidenceItem.timestamp` still varies under a real clock and still lands in
  `evidence_items.parquet`. Nothing else does — see `docs/handoff-integrity.md` §4.
- **Database assertions require a versioned citation.** "AF = 0.0001" is
  meaningless without knowing which gnomAD release and which population.

## Consequences
- Writing a stage is more work: you must say what you know, how, and what it
  fails to establish.
- Reports become auditable end to end, and "what argues against this?" is a query
  rather than an interview.
