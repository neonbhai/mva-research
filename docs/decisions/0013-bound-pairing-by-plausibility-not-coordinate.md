# ADR 0013 — Bound the pairing hypothesis space by plausibility, not by coordinate

**Status:** accepted
**Date:** 2026-08-28
**Touches:** ADR 0005 (filtering ≠ ranking), ADR 0010 (filtering allele frequency), GP-13, GP-32

## Context

`generate_pairs` built every within-gene pair, sorted them by
`(gene, first coordinate, second coordinate, pair_id)`, and truncated at
`max_pairs_per_gene = 20`. Truncation was therefore **by genomic position**: the
kept set was "every pair involving the leftmost variant, then the second, …", so
once the cap bit, variants at the right-hand end of a gene became unreachable in
any pairing. At ten variants in a gene 44% of pairings survived; at twenty,
10.5%; at fifty, 1.6%.

Nothing bounded the input. `apply_hard_filters` removes only invalid records by
design, and `select_candidate_variants` — which does bound the set — was
exported, tested, and never called. So every common SNP and synonymous change in
a gene counted toward n.

Reproduced by execution: twenty calls in one gene, the two ultra-rare
HIGH-impact calls placed at the right-hand end. All twenty kept candidates were
anchored on the two leftmost variants; **the plausible pair was deleted before it
could be scored.** Truncated candidates already carried
`gene_pair_cap_truncated`, but nothing read it, so the run's output was
indistinguishable from a complete search.

For a challenge whose answer *is* a compound-heterozygous pair, silently deleting
the correct pair is the worst failure this pipeline can have.

## Decision

1. Truncation is ordered by **plausibility**: recessive-plausible allele
   frequency (membership of `select_candidate_variants`, ADR 0010's
   `min_allele_number` guard included), then call quality, then predicted impact
   for the gene, taking a candidate's **weakest** member. `sort_key()` remains
   the total, stable tiebreak (GP-30). Position is no longer a criterion.
2. A second bound, `max_pairing_variants` (default 24), caps how many of a gene's
   variants enter the quadratic pair enumeration, by the same ordering. A variant
   excluded from *pairing* still forms its own single-variant hypothesis.
3. Both bounds report themselves. `generate_pair_candidates` returns a
   `PairingResult` carrying `GeneCapEvent`s and a warning naming the flag, the
   genes and the counts — **gene symbols and counts only, never coordinates**
   (GP-41). The composition root surfaces it through `PipelineResult.warnings`.
4. The frequency criterion is used as an **ordering, never as a membership
   gate.** Nothing is removed from the hypothesis space on frequency grounds.

## Rejected alternative: feed `select_candidate_variants` into `generate_pairs`

Measured on the synthetic case: 12 flagged variants → 9 selected; 18 candidate
hypotheses → **9**. Half the hypothesis space deleted. The nine deleted include
every candidate `test_common_variant_pair_is_downranked_not_deleted` requires to
survive, so a hard gate fails a hash-locked golden test.

It would also shrink the submission from ten rows to roughly five — pure
forfeited upside, since the verified scorer confirms rows below the answer
**cannot lower either metric** (see the executed results in
`docs/references/track1-submission-contract.md`).

Making a frequency selection the boundary of the hypothesis space is the exact
asymmetry ADR 0005 exists to prevent: a deleted candidate is invisible and
unfalsifiable, and published pathogenic alleles above naive frequency
expectations in under-represented populations are precisely what a frequency
filter loses silently. **ADR 0005 stands unamended.**

## Consequences

- A cap can still delete a hypothesis. It now deletes the *least plausible* one
  and says so, rather than the rightmost one in silence.
- No weight or threshold in `config/` changed (ASSUMPTION-SCORING-03 untouched).
- Golden output unchanged: SYNTHKIN1 carries 5 variants (max 15 candidates), so
  neither bound fires on the synthetic case. The 18-candidate ranked list is
  identical in content, order and composites; `tests/golden` passes unmodified.
- Enforced by `tests/unit/test_pairing_cap.py` (11 tests), including that the
  kept set is not a coordinate prefix and that the warning carries no coordinates.
