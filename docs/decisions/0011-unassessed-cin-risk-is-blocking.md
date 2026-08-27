# ADR 0011 — An unassessed chromosomal-instability risk is disqualifying

**Status:** accepted · **Date:** 2026-08-28 · **Supersedes:** the SYNTH-DRUG-E row
of `tests/golden/expected_drug_outcomes.tsv`

## Context

`docs/scientific-assumptions.md` (ASSUMPTION-DRUG-07) states that
`worsens_chromosomal_instability = None` — meaning nobody has assessed whether an
agent increases aneuploidy or cancer susceptibility — "is itself a blocking gap in
this disease context, not a pass". The evidence schema, the drug model docstring
and the rendered Track 2 report all repeat that word: **blocking**.

An adversarial pharmacology review found it blocked nothing.
`interventions/safety.py` recorded the gap as a non-fatal concern, and the
reviewer verified a candidate with a blank `worsens_cin` cell was **accepted at
rank 2** and rendered in the "direction AGREES" section, directly beneath a
sentence telling the reader this was a blocking gap that must be answered before
any further consideration.

Prose was doing work the code refused to do. That is the worst configuration: a
reader trusts the sentence, and the sentence is not enforced.

## Decision

Make it true rather than delete the claim. In a chromosomal-instability disease
context, `worsens_chromosomal_instability is None` on a **disease-modifying**
candidate is **disqualifying**, recorded as `RejectionReason.ONCOGENIC_RISK` at
severity `critical`. Outside a CIN context it remains a non-fatal concern.

## Why blocking rather than a penalty

The patient's non-tumour cells already sit above the aneuploidy-tolerance ceiling
with no reserve. The specific harm this pipeline exists to prevent is advancing a
compound that pushes chromosome instability further. "We never checked whether
this agent does the thing that would hurt this child most" is not a rounding
error on a score — it is the one question that must be answered first.

Rejection is not deletion. The candidate is preserved in the rejection record
with its reasons (GP-19), and its "what would have to change" entry names the
missing measurement. The output is therefore *more* useful, not less: it says
exactly which experiment unblocks the candidate.

## Before / after

| Drug | Before | After |
|---|---|---|
| SYNTH-DRUG-A Synthexostat | accepted, rank 1 | accepted, rank 1 (unchanged) |
| SYNTH-DRUG-B Synthinib | rejected, `wrong_direction` | unchanged |
| SYNTH-DRUG-C Synthophore | rejected, `not_approved` | unchanged |
| SYNTH-DRUG-D Synthazepam | accepted, symptomatic | unchanged (`worsens_cin = 0`) |
| **SYNTH-DRUG-E Synthaxel** | **accepted, rank 3, direction undetermined** | **rejected, `oncogenic_risk`** |
| SYNTH-DRUG-F Synthramide | rejected, `target_not_in_mechanism` | unchanged |

## The consequence that had to be handled

E was the fixture's only carrier of a distinct and important lesson —
ASSUMPTION-DRUG-02, that a `CONTEXT_DEPENDENT` direction resolves to "cannot
determine", which is neither agreement nor a wrong-direction rejection. Rejecting
E for an unrelated reason would have removed that coverage silently, which is
exactly the kind of quiet erosion GP-32 exists to prevent.

So a seventh catalogue entry, **SYNTH-DRUG-G (Synthemod)**, was added: a
`CONTEXT_DEPENDENT` agent whose CIN effect **has** been assessed
(`worsens_cin = 0`). It is accepted with `directions_agree is None`, preserving
the lesson while E carries the new one. The fixture now distinguishes three
states that were previously conflated in one row: direction disagrees, direction
undetermined, and CIN risk unassessed.

## Alternatives rejected

- **Set E's `worsens_cin` to 0.** Rejected: E's own mechanism-of-action text says
  its effects on missegregation are bidirectional. Declaring it CIN-safe to keep
  a test green would be fabricating data to fit an expectation.
- **Soften the doc instead of the code.** Rejected: the claim is correct for this
  disease. The code was wrong, not the assumption.
- **Down-rank rather than reject.** Rejected here specifically because the gap is
  about the disease mechanism itself. GP-13's flag-don't-delete rule governs
  *variant* filtering, where the cost of a false negative is a missed diagnosis.
  Here the cost of a false positive is advancing a compound that may worsen the
  child's chromosome instability.

## Consequences

- `tests/golden/expected_drug_outcomes.tsv` is re-baselined for E, and the hash
  lock in `tests/golden/test_locked_files.py` updated in the same commit — which
  is the GP-32 procedure working as intended: the change forced this record.
- Real catalogues rarely characterise CIN effects, so on real data this rule will
  reject many candidates. That is the honest answer for this disease, and each
  rejection names the assay that would unblock it.
