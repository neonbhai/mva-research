# ADR 0010 — A population must be large enough to set the maximum allele frequency

**Status:** accepted · **Date:** 2026-08-27

## Context
Rarity is scored on the **maximum** allele frequency across recorded
populations, not the global figure (ASSUMPTION-FREQUENCY-02). The reasoning is
sound: a variant common in any single ancestry is not a plausible ultra-rare
cause, and the global AF systematically under-estimates frequency for alleles
enriched in cohorts that reference panels under-sample.

Taken raw, though, "maximum" hands the decision to whichever population happened
to be sampled *least*. `VariantRecord.max_allele_frequency` maximised
`allele_frequency` and ignored `allele_number`, a field the model already
stores. A review reproduced the failure on a plausible founder allele:

| population | AC | AN | AF |
|---|---|---|---|
| `global` | 1 | 125,000 | 8×10⁻⁶ |
| `ami` | 1 | 40 | 0.025 |

Both rows describe **one observed allele**. The second was allowed to set the
maximum, so the variant was flagged `common_variant`, its rarity component fell
from ~1.00 to 0.1196 (≈ −0.27 composite), it acquired a
`CONTRADICTS` evidence item asserting the allele exceeds the maximum plausible
frequency for a severe recessive disorder, and `select_candidate_variants`
dropped it from the plausible set entirely. A genuine founder allele in an
under-represented population — precisely the case the maximum-AF rule exists to
protect — is the case it destroyed.

This is a known problem with a known name. gnomAD's `grpmax` requires a genetic
ancestry group to be adequately sampled before it may set the group maximum, and
ACMG's BA1/BS1 rules are applied to a **filtering allele frequency** — the lower
bound of the 95% confidence interval on the observed AF — for exactly this
reason. One allele in forty chromosomes has a point estimate of 0.025 and a 95%
CI lower bound around 1×10⁻³; the point estimate is not the quantity anyone
should be thresholding.

## Decision
Add a configuration key `frequency.min_allele_number`, **default 2000**.

A population whose reported `allele_number` is below it may not set the maximum.
It is not discarded: `VariantRecord.select_max_allele_frequency` returns a
`FrequencySelection` carrying both the chosen row and every excluded row, and
the rarity rationale names each excluded cohort with its AC/AN. A reader who
disagrees with the threshold can see exactly what it cost them.

Three details are deliberate:

* **A population with no `allele_number` at all stays eligible.** An unreported
  cohort size is unknown, not small, and excluding it would silently discard the
  only record a source supplied (GP-14).
* **When no population is large enough, `observed` is `None`.** Downstream this
  reads as absence of frequency data — scored at the configured mid-range
  `absent_frequency_score`, flagged `no_frequency_data`, and raised as a blocking
  open question. That is the honest answer: a frequency measured on forty
  chromosomes establishes neither commonness nor rarity.
* **The default on the model method is `0`**, i.e. no guard. Only a caller that
  passes the configured threshold gets the new behaviour, so no code path
  acquires it by accident.

### Why 2000

`min_allele_number` is a heuristic like every other cut-point here (GP-32); it
is not fitted to anything. The reasoning:

* At AN = 2000 a **singleton** yields AF = 5×10⁻⁴. That sits inside the
  low-frequency band and is 20× below `max_plausible_recessive` (0.01), so no
  *single* observation in an eligible population can push a variant across the
  commonness threshold on its own. It takes AC = 20 to reach 0.01, which is a
  claim about a real recurrence rather than about sampling.
* The rule of three gives a 95% upper bound of ~1.5×10⁻³ for zero observations
  at AN = 2000, so an eligible population that reports AF = 0 is genuinely
  informative rather than merely under-sampled.
* It is far below the AN of any headline gnomAD ancestry group (tens of
  thousands), so on real reference data the guard binds only on the small
  sub-populations where it is meant to.

Anything in the 1,000–5,000 range would behave similarly. The value is not
load-bearing to the demo case and is expected to be revisited against real
gnomAD releases; see `docs/tech-debt.md`.

## Consequences
- Before/after on the synthetic case: **no change**. Every population in
  `knowledge/public/frequencies.tsv` reports AN ≥ 67,200, so no row is excluded
  and the golden ranking is byte-identical (SYNTHKIN1
  `chr15:40200000:C:T` + `chr15:40210500:G:A` still ranks #1 at composite
  0.8779). `tests/golden/expected_ranking.tsv` is untouched.
- `config/default.yaml` changes, so its hash lock in
  `tests/golden/test_locked_files.py` is updated in the same commit — which is
  what that lock is for.
- The guard can now suppress the only frequency record a variant has. That is a
  visible, testable outcome (`no_frequency_data` + a blocking question), not a
  silent one, and it never removes a candidate (GP-13).
- `max_allele_frequency()` gained a keyword argument. Callers that do not pass
  it keep the unguarded behaviour, so a future call site can forget the guard.
  The four call sites that matter (`_frequency_flag`,
  `select_candidate_variants`, `_variant_rarity`, `collect_contradictions`) all
  pass it, and the rarity test asserts the AC=1/AN=40 row is excluded end to end.
