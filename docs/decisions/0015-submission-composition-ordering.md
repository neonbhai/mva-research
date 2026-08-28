# ADR 0015 — A pair outranks the single it subsumes

**Status:** accepted; amended by ADR 0023
**Date:** 2026-08-28
**Touches:** ADR 0014 (EPCR injectivity), ADR 0023, GP-32

> **Amended by ADR 0023 in two places.** (1) "Applied after truncation" below was
> a defect worth 50 rank points: it deleted any pair ranked below the cut before
> promotion could look for it. Promotion now runs on the full ranked list, before
> truncation. (2) "The single is kept, not dropped" holds only while keeping it is
> free. When a pair has to cross the ten-row cut, it is **exchanged** with its own
> subset rather than displacing an unrelated row.

## Context

The Track 1 scorer awards rank points only for the **best full match**, and a
full match is frozenset equality: `row.variants == true_variants`. Verified
against the real `evaluation.py`.

`_drop_subsumed` dropped a single-variant candidate that ranked *below* a pair
containing it. It did nothing in the opposite case. So when a single outscored
its own parent pair, both were emitted with the single first:

```
chr15:40200000 / -                  epcr=0.9010   ← single, row 1
chr15:40200000 / chr15:40210500     epcr=0.8020   ← the pair, row 2
```

If that pair is the answer, this costs **100 → 50 rank points** and buys nothing:
F-max is 1.0 either way, because the pair re-emits the single's variant so the
prediction set at the top threshold is identical.

## Decision

**When a single-variant candidate is a strict subset of a pair candidate that is
also being emitted, the pair ranks above the single. Both rows are kept.**

Applied **before** EPCR separation, so separation is handed the order we actually
mean to submit and remains order-preserving. (This ADR said "after truncation";
ADR 0023 corrected it to *before* truncation, and priced the difference.)

Three reasons for this shape rather than the alternatives:

1. **It is a bet, and the bet is one-sided.** No arrangement puts both at rank 1.
   The challenge states the answer is a clinically validated
   **compound-heterozygous pair**; `evaluation.py` computes
   `is_compound_het = len(true_variants) == 2` and offers partial credit only in
   that branch; `groundtruth.py`'s own local fallback holds two variants.
   Promoting the pair bets with the challenge's explicit statement.
2. **The single is kept, not dropped** — *while keeping it is free*. Rows below
   the answer cannot lower either metric (verified), so dropping the single would
   forfeit a full match in the low-probability world where the truth is
   single-variant and gain nothing in the likely one. ADR 0023 bounds this: once a
   pair has to cross the ten-row cut, keeping the single costs an unrelated pair's
   slot, and the pair takes the single's slot instead.
3. **It is composition, not scoring.** No composite score is touched. Adjusting
   scores to achieve this ordering would make a scoring change masquerade as a
   rendering fix, and the ranking is a scientific judgement that must stay
   auditable as one.

## Before / after on the golden case

*(As of this ADR. ADR 0023 moved promotion before truncation and records the
one further row that changed.)*

Promotion fires where predicted — the single `chr15:40206000` sat above two pairs
carrying it:

| row | before | after |
|---|---|---|
| 6 | `40206000` single, 0.5273 | `40200000 / 40206000` pair, 0.4420 |
| 7 | `40200000 / 40206000` pair, 0.4420 | `40206000 / 40210500` pair, 0.4320 |
| 8 | `40206000 / 40210500` pair, 0.4420 | `40206000` single, **0.4220** |

Same ten rows, same coordinates, nothing added or removed. Rows 1–5 and 9–10
unchanged. **The causal pair remains row 1 at REDACTED-EPCR.** Ten distinct EPCRs,
minimum gap 0.0100.

Cumulative delta from the pre-ADR-0014 baseline: three EPCR values moved
(0.6148→0.6110, 0.4420→0.4320, 0.3900→0.3800) and one row repositioned (the
`40206000` single, row 6 → row 8). No coordinate changed.

## Termination

A promoted pair can have no superset of its own — two variants is the maximum a
row carries — and each promotion removes one subset-above-superset inversion
without creating another. The pass strictly decreases inversions over a fixed
scan order, so it terminates and is deterministic (GP-30).

## Consequences

- Enforced by `tests/unit/test_track1_composition.py` (8 tests), including that
  the mirror rule still drops a single *below* its pair, that a
  partially-overlapping pair is **not** promoted (strict subset only), and that
  unrelated rows keep their order — promotion is surgical, not a re-sort.
- If the ground truth ever turns out to be single-variant, this decision costs 50
  rank points and should be revisited. That is the accepted risk, stated here so
  the reversal is cheap.
