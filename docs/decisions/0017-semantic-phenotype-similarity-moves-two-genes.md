# ADR 0017 — Real ontology similarity moves two genes; the movement is the point

**Status:** accepted
**Date:** 2026-08-28
**Touches:** GP-14, GP-32, TD-04, ASSUMPTION-PHENOTYPE-03/04/05

## Context

TD-04 said an observed child term did not credit its associated parent, so true
matches were under-scored. Replacing exact-term matching with real HPO graph
semantics — information content from the actual annotation corpus, Lin
similarity, Köhler symmetric best-match — necessarily changes scores. Under GP-32
that requires arbitration and a record rather than a re-baseline.

Measured on the golden case by running `execute_pipeline` with the real release:

| Gene | phenotype before → after | Δ | composite Δ | rank |
|---|---|---|---|---|
| SYNTHKIN1 (all 14 pairs) | 1.000000 → 1.000000 | 0 | 0 | unchanged |
| SYNTHMET2 | 0.625000 → **0.682407** | +0.057407 | +0.008037 | 4 → 4 |
| SYNTHMUL4 (2 pairs) | 0.500000 → **0.300333** | −0.199667 | −0.027953 | 7 → 7, 15 → 15 |
| SYNTHOTH3 | 0.000000 → 0.000000 | 0 | 0 | 14 → 14 |
| SYNTHSOL5 (unranked) | 0.500000 → 0.250000 | −0.250000 | — | — |

**Rank 1 is byte-identical** and the complete ordering 1..18 is unchanged.
`expected_ranking.tsv` and `expected_drug_outcomes.tsv` remain correct as written
and were not touched.

## Decision

**Both movements are accepted as correct. The scoring rule is not adjusted to
suppress them.**

SYNTHMUL4's −0.1997 is the whole purpose of the change. Its only curated term
(Seizure) is `UNCERTAIN` in the **synthetic** case's profile
(`tests/fixtures/synthetic/synthetic_phenotype.tsv`), so specificity stays undefined and
contributes the neutral 0.5 exactly as before. The entire drop is the new
coverage component discovering that Seizure is semantically **distant** from all
four observed features (Lin best-match average 0.1007).

That is *"we know this gene's feature and it does not fit"* — which the old exact
matcher could not distinguish from *"we know nothing about this gene"*. Both
scored 0.5. Making them different is the fix, not a regression.

The obvious lever — make an undefined component contribute `NEUTRAL_SCORE` so
SYNTHMUL4 stays at 0.5 — was rejected. It would restore precisely the conflation
TD-04 exists to remove, and would do so in the direction that flatters weak
candidates.

## Why this is not a GP-14 violation

The distinction GP-14 protects is preserved, and can be read off the numbers:

- `GeneAnnotationStatus.NO_ANNOTATIONS` — nothing curated — scores **0.500** with
  `coverage=None` and emits its own evidence item saying so.
- A gene *with* annotations that fit nothing scores **0.250**.

Demonstrated on the golden case: `GENE_WITH_NO_TERMS` 0.500 vs SYNTHSOL5 0.250.
Absence still costs nothing; only knowledge that fails to fit costs something.
A `NOT_ASSESSED` association is tested to leave `specificity`, `matched_weight`
and `informative_weight` *exactly* unchanged.

## Two properties worth recording

**Contradictions are withdrawn from scoring, not resolved.** When the same term
is recorded both OBSERVED and EXCLUDED — including under an `alt_id` and its
primary, which resolve to one term — the result is `UNCERTAIN` and conflicted:
the only status that contributes zero in both directions and propagates nowhere.
Pinned three ways: identical to the gene *without* that association, and
different from both the if-OBSERVED-won and if-EXCLUDED-won scores.

**The corpus reader now fails closed.** A ragged annotation row raises at load
rather than silently vanishing. Information content is the denominator of every
similarity score, so a dropped row perturbs *all* phenotype scores while the
stats claim a clean parse. The trade is a new way the pipeline can refuse to
start, on a future release shipping a malformed row. That is accepted
deliberately: refusing to start is recoverable, and silently wrong science is not.
The pinned release parses bit-identically (`rows_read=285598, rows_kept=267062`,
max IC 9.467692).

## Consequences

- Scored similarity is deliberately asymmetric (ASSUMPTION-PHENOTYPE-05); the
  reverse direction is computed and reported but not scored on.
- IC comes from the corpus, never from graph depth (ASSUMPTION-PHENOTYPE-04),
  and is undefined — not zero — for terms the corpus never reaches.
- **These numbers are not calibrated.** No known-answer validation set exists
  (TD-02), so "more correct" here means "represents the distinction it claims to
  represent", not "measurably more sensitive".
