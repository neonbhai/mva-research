# ADR 0005 — Filtering and ranking are separate concerns

**Status:** accepted · **Date:** 2026-08-27

## Context
The conventional rare-disease pipeline filters hard and early: drop everything
above an allele-frequency cut-off, drop non-PASS calls, drop low-impact
consequences, then rank what survives. It is fast and it usually works.

It also fails silently. A filtered variant leaves no trace, so a causal variant
removed at step one is indistinguishable from one that never existed. Published
pathogenic variants exceed naive frequency expectations in under-represented
populations; real causal variants get non-PASS flags in repetitive regions.

## Decision
**Hard filters remove only invalid or genuinely impossible records.** Everything
else is flagged and down-ranked.

Hard filter (removal) is limited to: wrong genome build, non-canonical contig,
`*`/`.` alt allele, hom-ref or missing genotype. These are not weak candidates —
they are records the pipeline cannot form a coherent hypothesis about.

Soft flags (retention + down-ranking): common, low-frequency, no frequency data,
benign consequence, low depth, low GQ, caller-filtered, allele-balance outliers,
possible mosaic.

## Consequences
- The output carries weak candidates. That is intended: a ranked list with a
  visible bad candidate at position 40 is more informative than a short list
  whose omissions are invisible.
- Every removal is counted and reported by reason code, so the number of dropped
  records is always visible.
- Scoring must be robust to junk rather than assuming a pre-cleaned input.
- Reviewers can ask "where did variant X go?" and get an answer.
