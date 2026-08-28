# ADR 0014 — EPCR must be injective; magnitude is not a scientific claim

**Status:** accepted
**Date:** 2026-08-28
**Touches:** GP-32, TD-10, ASSUMPTION-SCORING-03

## Context

The Track 1 scorer was read directly from the challenge Space
(`evaluation.py`) and its logic re-implemented and executed. Two facts follow,
and they are verified rather than inferred:

1. F-max thresholds are **the unique EPCR values present in our own submission**,
   swept descending, compared `row.epcr >= t`. Not a fixed grid.
2. Rank is `epcr` descending with **ties broken by file order**.

Therefore both metrics depend only on the **order and tie structure** of our
EPCRs, never on their magnitude. Executed results, true pair at 0.90:

| submission shape | rank points | F-max |
|---|---|---|
| true pair, no ties | 100.0 | 1.0000 |
| true pair tied with a wrong pair | 100.0 | 0.6667 |
| same tie, wrong row first in file order | **50.0** | 0.6667 |

The rendered demo submission was emitting **four tied EPCRs** (`0.4420` twice,
`0.3900` twice). A tie at the answer's threshold pulls a wrong row into the same
prediction set, halving precision; if file order goes against us it also takes
rank 1. The exposure is up to **50 rank points and 0.333 F-max simultaneously**.

## Decision

**Emitted EPCRs are strictly decreasing with a minimum separation of 0.01,
enforced after truncation to `max_rows`, order-preserving by construction.**

- Arithmetic is integer (`_EPCR_UNIT_SCALE`); repeated float subtraction of 0.01
  is not byte-reproducible and GP-30 requires it to be.
- A reserved floor guarantees `epcr ∈ (0, 1]` even with ten rows at composite
  0.0. A naive "subtract 0.01 per row" emits a negative EPCR by row 2, which the
  scorer rejects with `ValueError` — losing the entire submission.
- **No row is ever dropped to break a tie.** Rows below the answer cannot lower
  either metric, so dropping one is pure forfeited upside.
- `validate_submission` rejects any repeated EPCR on the rendered bytes.

**This is a rendering decision, not a scoring one.** The pass cannot reorder: it
walks rows in the order given and assigns each a value strictly below its
predecessor. It turns a collision into a separation and does nothing else. No
weight or threshold in `config/` changed, so ASSUMPTION-SCORING-03 is untouched.

**TD-10 is re-scoped.** It claimed EPCR miscalibration costs points. The
invariance above shows calibration is worth *exactly zero*: any strictly
increasing reparameterisation scores identically. What is worth a great deal is
**injectivity**. TD-10 becomes a reporting-honesty debt — an EPCR presented to a
clinician as a probability had better be one — not a scoring debt.

## Before / after on the golden case

| row | before | after |
|---|---|---|
| 4 `chr15:40210500 / chr15:40211000` | 0.6148 | **0.6110** |
| 8 `chr15:40206000 / chr15:40210500` | 0.4420 | **0.4320** |
| 10 `chr15:40205000 / chr15:40210500` | 0.3900 | **0.3800** |

Rows 1, 2, 3, 5, 6, 7, 9 byte-identical. No row moves, no coordinate changes.
Ten distinct values, minimum gap exactly 0.0100.

**The invariant this had to demonstrate is that the separation pass never
reorders**, and in particular never displaces the top-ranked candidate. That is
asserted on the synthetic golden case rather than on any real result:
`tests/golden/test_golden_case.py::test_synthetic_causal_pair_ranks_first`, plus
`tests/unit/test_track1_composition.py::test_unrelated_rows_keep_their_epcr_order`
and `::test_a_genuine_score_difference_still_decides_the_order`. Real-case EPCRs
and ranks are patient-derived and are not reproduced in this repository.

## The enforcement gap this exposed — and it is the important part

These values changed across runs and **`tests/golden` still passed.**
`expected_ranking.tsv` locks `rank, gene_symbol, variant_a, variant_b,
inheritance_model, must_be_flagged` — it does **not** lock EPCR. So the one
quantity the challenge actually scores on was free to move without tripping the
gate that exists to force a decision record.

GP-32 says golden expectations are never silently re-baselined. That guarantee
was only as good as the columns under lock, and the scored quantity was not one
of them. A rule enforced over the wrong projection of the output is not enforced.

**Required follow-up:** add a golden expectation covering the emitted submission
EPCRs, hash-locked in `tests/golden/test_locked_files.py`. Deferred only until
the submission-composition ordering settles (the pair-before-subset-single rule),
because locking a value that is about to change on purpose teaches the team to
re-baseline — the precise habit GP-32 exists to prevent. Tracked in
`docs/next-actions.md`.

## Consequences

- Enforced by `tests/unit/test_epcr_separation.py` (10 tests), including that
  emitted order equals an independently recomputed pre-fix order — the test does
  not compare the pass against itself.
- Ten rows consume 0.09 of the 0.99 available range.
