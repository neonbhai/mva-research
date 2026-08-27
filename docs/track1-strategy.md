# Track 1 scoring strategy

How the challenge's two metrics actually respond to the shape of our submission,
where our ranking is most likely to bury the answer, and what to change.

**Provenance of every claim below.** The scorer's mechanics are taken from
`docs/references/track1-submission-contract.md`, which was derived by reading the
Space's source once. **I did not execute the real scorer.** Everything labelled
*derived* is arithmetic over a re-implementation of the contract as written
(rank tiers, EPCR-descending rank, frozenset matching, per-variant F-max union);
everything labelled *inferred* is a reading of contract prose where the prose
admits more than one implementation. Inferred items are collected in §6 as a
re-verification list. Everything about our own pipeline is read from code in
this repo and is checkable.

**No number here is a sensitivity estimate.** We have no known-answer validation
set (TD-02), so no statement of the form "this raises our hit rate to X%" appears
anywhere in this document, and none should be added. What *is* computable is the
conditional arithmetic: *given* the answer lands at rank r, or *given* a tie
occurs, here is exactly what the score does.

---

## 0. The scorer, reduced to two facts

Write `G` for the ground-truth variant set. The contract says it is one
clinically validated compound-heterozygous pair, so `|G| = 2`.

**Fact 1 — rank points depend only on the EPCR *order* of rows.**
Rank is `sorted(enumerate(rows), key=lambda x: (-epcr, file_index))`; tiers are
`[(1,100),(3,50),(5,25),(10,10)]`; a partial match scores half.

**Fact 2 — F-max depends only on the EPCR order *and the tie structure*.**
Thresholds are the unique emitted EPCR values, and the predicted set is
`predicted_variants |= row.variants` over rows at or above the threshold. So the
family of sets the sweep can ever see is exactly the family of *prefixes of the
EPCR-descending row order that end on a tie boundary*. Rows sharing an EPCR value
cannot be separated by any threshold.

### The invariance theorem, and what it does to TD-10

Both facts are properties of an ordering. Therefore:

> Let `φ` be any strictly increasing map from `(0, 1]` to `(0, 1]`. Replacing
> every emitted `epcr` with `φ(epcr)` changes neither the rank points nor the
> F-max.

Verified by construction over identity, `e²`, `√e`, an affine squeeze of all ten
rows into `[0.50, 0.52]`, and a pure rank index `1.00, 0.91, …, 0.19`: all five
produce identical rank points and identical F-max to six decimal places.

**Consequence: EPCR *calibration* is worth zero points under this scorer.** TD-10
says "miscalibration costs points even when the ranking is right". Under the
contract as written that is false. `epcr = 0.60` meaning "we rank this sixth"
rather than "60% likely" costs nothing at all. TD-10 is a **reporting** debt, not
a scoring debt, and should be re-scoped to say so — with one exception, which is
the whole of §3: the map must be **injective**. Ties are the only property of the
EPCR vector that the scorer can see, and they are pure loss.

---

## 1. Row budget: submit ten, under one condition

### 1.1 Rank points are monotone

Adding a row whose EPCR is **strictly below** every EPCR already present appends
it to the sorted order and shifts nothing. No existing row's rank changes, so the
best matching row's rank cannot increase. Rank points are non-decreasing.

The only way an added row costs rank points is a **tie**: a wrong row emitted at
the same EPCR as the true row, ordered ahead of it by file order, demotes the
true row by one. Across a tier boundary that is expensive — rank 1 → rank 2 is
−50, rank 3 → rank 4 is −25, rank 5 → rank 6 is −15.

### 1.2 F-max is monotone under the same condition

Adding a row at a strictly lower EPCR than all existing rows leaves every
previously available threshold with an unchanged predicted set (the new row is
below all of them) and adds one new, lower threshold. The sweep therefore
maximises over a strict superset of the values it saw before. F-max is
non-decreasing. Derived, and confirmed empirically: holding the true pair at
rank 3 and growing the submission from 3 to 10 rows leaves rank points at 50.0
and F-max at 0.5000 throughout, unchanged at every step.

### 1.3 The exact cost of a wrong row *above* the answer

Let `k` be the number of non-truth variants in the predicted set at the true
row's own threshold — i.e. the distinct wrong variants carried by rows ranked
above it, plus any wrong partner in the true row itself. Then at that threshold:

| what is present | precision | recall | F1 |
|---|---|---|---|
| both true variants, plus `k` wrong | `2/(2+k)` | `1` | **`4/(4+k)`** |
| one true variant, plus `k` wrong | `1/(1+k)` | `0.5` | **`2/(3+k)`** |

Derived and machine-checked. So, with every row above the answer being a wrong
**pair** (`k = 2(r−1)`):

| true pair at rank | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 |
|---|---|---|---|---|---|---|---|---|
| rank points | 100 | 50 | 50 | 25 | 25 | 10 | 10 | 10 |
| F-max | 1.000 | 0.667 | 0.500 | 0.400 | 0.333 | 0.286 | 0.222 | 0.182 |

**"What happens to F-max precision when a row's two variants are both wrong?"** —
each such row above the answer adds exactly 2 to `k`. The marginal cost is
`4/(4+k) − 4/(4+k+2)`, which is 0.333 for the first such row, 0.167 for the
second, 0.100 for the third, and under 0.02 by the eighth. **A wrong pair
ranked *below* the answer costs exactly nothing**, because the threshold at the
answer's own EPCR is still in the sweep and the maximum is taken over it.

This asymmetry is the whole row-budget argument: cost is incurred only *above*
the answer, and rows 6–10 are by construction below rows 1–5.

### 1.4 Recommendation

**Submit ten rows, always.** There is no scenario in the contract where a row
placed strictly below every existing row reduces either metric. The only failure
mode is a tie, and a tie is a rendering defect we control, not a property of
having ten rows.

**The tie hazard is live today.** `.demo-workspace/runs/synthetic-demo-*/submission/track1_submission.csv`
— the golden output — contains four tied EPCRs:

```
... chr15:40200000 / chr15:40206000 , 0.4420
... chr15:40206000 / chr15:40210500 , 0.4420
... chr15:40200000 / chr15:40205000 , 0.3900
... chr15:40205000 / chr15:40210500 , 0.3900
```

They are harmless there because the answer is at rank 1. Move the answer to
either of those pairs and the tied partner's variants are unseparable from it at
every threshold. A tie at the answer's own threshold costs, in the worst case:

```
true pair alone at 0.9000            rank_pts=100.0  Fmax=1.0000
+ wrong pair tied at 0.9000          rank_pts=100.0  Fmax=0.6667   (−0.333 F-max)
+ wrong pair separated at 0.8999     rank_pts=100.0  Fmax=1.0000   (no loss)
wrong pair first, tied at 0.9000     rank_pts= 50.0  Fmax=0.6667   (−50 and −0.333)
```

Ties arise because `composite_to_epcr` rounds to `EPCR_DECIMALS = 4` and
`build_submission_rows` then falls back to file order — the module docstring says
so explicitly. Two composites within ~1.01×10⁻⁴ collide. See §5, item 2.

---

## 2. Pair rows versus single-variant rows

### 2.1 The reading that decides it

Contract: *"Matching is exact frozenset equality over the row's variants; the
first (best) full match wins. Partial credit applies to compound-het rows on set
intersection."* Two readings (both **inferred**, §6 item 4):

- **Strict** — partial credit requires the *row* to be a compound-het row (both
  `_2` columns populated). A single-variant row can then never score: `{g₁}` is
  not frozenset-equal to `{g₁, g₂}` and is not eligible for intersection credit.
- **Lenient** — partial credit is computed from the intersection with a
  compound-het *ground truth*, whatever the row's shape. A single-variant row
  carrying one true variant then scores `0.5 × rank_points`.

### 2.2 Arithmetic for the three realistic cases

One pair row at rank 1, versus the same two variants split across two single rows
at ranks 1 and 2:

| case | pair row: rank pts | pair row: F-max | split: rank pts (strict / lenient) | split: F-max |
|---|---|---|---|---|
| both right | **100** | 1.000 | 0 / 50 | 1.000 |
| one right | **50** | 0.500 | 0 / 50 or 25 | 0.667 |
| both wrong | 0 | 0.000 | 0 / 0 | 0.000 |

Derived. Splitting is **weakly worse under the lenient reading and catastrophic
under the strict one**, and it burns a second slot. In the "one right" case
splitting buys `2/3 − 1/2 = +0.167` F-max — the wrong partner is one fewer false
positive at the good threshold — which is not remotely worth 50–100 rank points.

### 2.3 The hedge case, which has a clean dominance answer

Suppose we are confident in variant `a` and unsure of its partner:

```
single {a}         rank_pts =   0 (strict) / 50 (lenient)   Fmax = 0.667
pair  {a, wrong}   rank_pts =  50 (strict) /  50 (lenient)  Fmax = 0.500
```

**Pairing weakly dominates under both readings** — equal under lenient, +50 under
strict — for a cost of 0.167 F-max. Since the two metrics' combination weight is
unknown (§6 item 5), taking a guaranteed 50 rank points over a certain 0.167
F-max is the right side of that bet under any plausible weighting.

**Rule: never emit a single-variant row as a hedge on a pair hypothesis.** Emit
single-variant rows only where the *hypothesis itself* is single-variant —
homozygous, hemizygous, mitochondrial, mosaic, or a lone HIGH-impact het under a
dominant model. `_wants_single_candidate` in `src/mva/prioritization/pairing.py`
already draws exactly that line, and it is the right one. `_drop_subsumed` in
`src/mva/reporting/track1.py` already removes a single row subsumed by a kept
pair, which is also correct: under the strict reading such a row can never score,
and under either reading it adds no variant to the F-max union.

### 2.4 The corollary worth exploiting: recombination rows are free

A row whose variants are **all already present in higher-ranked rows** adds
nothing to the predicted set at any threshold at or above its own. Its `k`
contribution is zero everywhere it could matter. It is therefore an F-max-free
extra ticket at an exact-frozenset match.

This is not hypothetical — the golden demo already does it by accident. Its ten
rows carry only **eight distinct variants**: rows 4, 7, 8 and 10 introduce no new
coordinate at all, and rows 3, 9 introduce one each.

The actionable version: **when confidence is concentrated on a gene rather than
on a pairing, spend slots 2–10 on that gene's other within-gene pairings before
spending them on a second gene's pair.** A second gene's pair costs 2 false
positives at every threshold below it; a same-gene repairing of already-emitted
variants costs 0 or 1. Both give one extra shot at the exact match.

The one caveat is §6 item 3: if the scorer short-circuits on the *first* match
rather than taking the best, a fan-out row that partially overlaps `G` and sits
*above* the full match could cap us at half credit. Keeping fan-out rows inside
one rank tier bounds that loss at half the tier's value.

---

## 3. EPCR: the optimisation is degenerate, and that is the finding

The question was: *what shape of EPCR distribution across our ten rows maximises
expected F-max, given we do not know which row is right?*

By the invariance theorem in §0, **the shape is irrelevant.** Expected F-max is

```
E[F-max] = Σ_r  P(answer at rank r) · 4/(4 + k(r))
```

and neither term contains an EPCR magnitude. `P(answer at rank r)` is a property
of our *ranking*; `k(r)` is a property of which *variants* the rows above carry.
The only free parameter the EPCR vector actually exposes is its **tie structure**,
and the optimum there is trivially "no ties": ten distinct values give the sweep
ten thresholds, i.e. maximum resolution, and every collision merges two prefixes
that the sweep could otherwise have separated.

So the whole optimisation collapses to two rules:

1. **Injective.** Ten rows, ten distinct EPCR values, strictly decreasing in
   file order. Non-negotiable; see §1.4 for the cost of violating it.
2. **Well-separated, as insurance.** Under the contract the magnitudes cannot
   matter. Under one plausible alternative implementation they can: if the sweep
   runs over a **fixed grid** (`np.linspace`, or two-decimal steps) rather than
   over our emitted values, then rows closer together than the grid step become
   indistinguishable exactly as if they were tied. This is §6 item 2.

The cost of that second failure mode is measurable. Take ten realistically
bunched composites (`0.7201, 0.7199, 0.7195, 0.6903, …`) with the answer at
rank 1 and re-round the affine map at coarser resolution:

| rounding | distinct EPCRs of 10 | rank points | F-max |
|---|---|---|---|
| 4 dp (current) | 10 | 100 | 1.000 |
| 3 dp | 6 | 100 | 0.667 |
| 2 dp | 4 | 100 | 0.500 |
| 1 dp | 3 | 100 | 0.286 |

Derived. Note rank points never move — only F-max. Since our composites *do*
bunch (the demo's rows 7–10 span 0.052), buying separation costs nothing and
insures against a grid-based sweep.

**Recommended map.** Keep `composite_to_epcr` affine and monotone — it is right
and its docstring's reasoning is right — then apply a deterministic post-pass in
`build_submission_rows` that walks the EPCR-descending rows and enforces
`epcr[i] ≤ epcr[i−1] − δ` with `δ = 0.01`, clamped at `EPCR_FLOOR`. Ten rows need
0.09 of the 0.99 available range, so this never hits the floor. The pass is
order-preserving by construction, so it cannot change a scientific conclusion; it
can only turn a collision into a separation.

**And re-scope TD-10.** Its stated cost — "the F-max metric sweeps EPCR
thresholds, so miscalibration costs points even when the ranking is right" — does
not hold under the contract. The real debt is that EPCR is presented in a public
artifact under a name that reads as a probability. That is a reporting-honesty
issue, and the docstring already handles it in words.

---

## 4. Where our ranking puts the answer below rank 1, ranked by expected cost

Config values from `config/default.yaml`; composite deltas computed from
`composite_score` with `rarity`, `phenotype_similarity` and `mechanistic_relevance`
held at their neutral 0.5 (which, per item 1, is where they actually sit).

### 4.1 Most of the weight is dead on a real case — 0.40 of 1.00

`src/mva/orchestrator.py` builds `GenePhenotypeIndex` from
`knowledge/public/gene_phenotype.tsv` and `MechanismLibrary` from
`knowledge/public/mechanisms.tsv`. Both are synthetic: their gene symbols are
`SYNTHKIN1`, `SYNTHMET2`, `SYNTHOTH3`… On a real VCF, **no candidate gene matches
either table**, so `score_gene_phenotype` returns its `NEUTRAL_SCORE` of 0.5 for
every gene and `mechanism_relevance_score` its neutral for every gene. That pins
`phenotype_similarity` (0.14) and `mechanistic_relevance` (0.08) to a constant
0.11 contribution shared by every candidate: **0.22 of weight with exactly zero
discriminating power.**

With TD-01 unfixed, `knowledge/public/frequencies.tsv` is likewise synthetic, so
`select_max_allele_frequency` finds nothing, `_variant_rarity` returns
`absent_frequency_score = 0.5`, and `rarity` (0.18) is a second constant. Two of
the six `evidence_quality` facets — `frequency_provenance` and
`clinical_curation` — go to zero for every candidate for the same reason.

*Status note, checked against the working tree.* Real-source modules exist —
`src/mva/annotation/clinvar_vcf.py`, `src/mva/annotation/gnomad_sites.py`,
`src/mva/annotation/snpeff_local.py`, and real HPO ontology/propagation/
similarity modules under `src/mva/phenotype/` — but **none of them is imported by
`src/mva/orchestrator.py`**. Until they are wired into the runtime path, the
analysis in this section describes what actually runs. Wiring them is what
converts §5 item 7 from the most expensive recommendation into the cheapest.

**Live discriminating weight on a real case: 0.60** (analytical validity 0.20,
molecular consequence 0.18, inheritance 0.16, evidence quality 0.06), and the
first three take values from small discrete sets. Enumerating the reachable
composites over the documented multiplier products gives 728 distinct values from
3,136 component combinations — **ties are the norm**, and `_sort_key` in
`src/mva/prioritization/ranking.py` breaks them by ascending genomic position.
Ascending position over `CANONICAL_CONTIGS` means **chr1 systematically wins
ties**. On a real exome this is close to ranking by chromosome number.

This is the largest single gap and it is not a weight problem — no reweighting
fixes a component that is constant.

### 4.2 The within-gene candidate cap can delete the answer before it is scored

`generate_pairs(flagged, max_pairs_per_gene=20)` builds **all** `C(n,2)`
within-gene pairs, sorts them by `PairCandidate.sort_key()` — which is
`(gene, first coordinate, second coordinate, pair_id)` — and truncates to the
first 20. **Truncation is by genomic coordinate, not by plausibility.** The kept
20 are "every pair involving the leftmost variant, then every pair involving the
second, …", so once the cap bites, whole variants at the right-hand end of the
gene become unreachable in any pairing.

| variants in gene `n` | `C(n,2)` | kept | fraction of pairings reachable |
|---|---|---|---|
| 6 | 15 | 15 | 100% |
| 7 | 21 | 20 | 95% |
| 8 | 28 | 20 | 71% |
| 10 | 45 | 20 | 44% |
| 15 | 105 | 20 | 19% |
| 20 | 190 | 20 | 10.5% |
| 50 | 1225 | 20 | 1.6% |

Single-variant candidates consume cap slots too and sort *ahead* of the pairs, so
the real figures are lower. Whether `n` reaches 7 depends on input density, and
here is the thing: **nothing in the pipeline reduces `n` before pairing.**
`apply_hard_filters` removes only invalid records by design (ADR 0005), and
`select_candidate_variants` — which *does* bound the set to alleles below
`max_plausible_recessive` — is exported, tested, and **never called by
`src/mva/orchestrator.py`**. So every common SNP, synonymous change and (on WGS)
intronic call in a gene counts toward `n`.

The failure is silent in the score and loud in the artifacts: truncated
candidates carry `FLAG_CAP_TRUNCATED` (`gene_pair_cap_truncated`). Nothing
currently gates on it. That is the cheapest possible fix (§5 item 1).

### 4.3 Homozygous singles systematically outrank compound-het pairs — +0.051

This is the TD-07 consequence, but the framing in the tech-debt table understates
what actually happens. Phase is `UNKNOWN` for every pair in a proband-only VCF,
so `score_inheritance` returns `1.00 × phase_weights.unknown = 0.55`. The
absolute ceiling loss is `0.16 × (1.00 − 0.55) = 0.072`, but that is uniform
across compound-het candidates and does not reorder them.

What reorders is the comparison against candidates that pay **no phase
multiplier at all**, because `score_inheritance` applies it only when
`pair.is_pair`:

| inheritance model | component | weighted |
|---|---|---|
| homozygous / X-linked recessive (single) | 0.90 | 0.1440 |
| compound het × TRANS_CONFIRMED | 1.00 | 0.1600 |
| **compound het × UNKNOWN (our only reachable state)** | **0.55** | **0.0880** |
| mitochondrial / mosaic (single) | 0.50 | 0.0800 |
| mixed-zygosity pair × UNKNOWN | 0.33 | 0.0528 |
| lone het (single) | 0.20 | 0.0320 |

Holding everything else equal and accounting for the `evidence_quality`
`independent_observations` facet (pair 1.0, single 0.5, worth +0.005 to the
pair):

```
compound-het pair, phase UNKNOWN   composite = 0.6900
homozygous single                  composite = 0.7410
→ every rare homozygous call outranks the answer class by +0.0510
```

Derived. In a real genome — especially with any run of homozygosity — rare
homozygous calls are not scarce. Each one is a row above the answer.

**The principled fix is not to raise `phase_weights.unknown`.** Raising it
because we know the answer is a compound het is precisely the reasoning
ASSUMPTION-SCORING-03 forbids. The defensible observation is an *internal
inconsistency*: a compound-het pair is discounted 0.55× for an unverified phase
assumption, while a homozygous single keeps 0.90 despite an equally unverified
zygosity assumption that the pipeline **itself raises as an open question** —
`_open_questions` emits `::true-homozygosity` ("Is the homozygous call true
homozygosity rather than hemizygosity?") for exactly these candidates and then
does not price it. Pricing that open question symmetrically is a change justified
by our own evidence model, not by the answer key.

### 4.4 Mosaic allele fraction read as artifact — up to −0.195

`_allele_balance_penalty` bands `VariantRecord.allele_fraction`:

| band | validity multiplier | contradiction | composite delta |
|---|---|---|---|
| 0.25–0.75 | ×1.00 | — | 0.000 |
| 0.10–0.25 (mosaic window) | ×0.75 | — | **−0.050** |
| > 0.75 (allelic imbalance) | ×0.60 | `low_quality_call` 0.30 | **−0.155** |
| < 0.10 (below mosaic floor) | ×0.40 | `low_quality_call` 0.30 | **−0.195** |

Derived, on a compound-het pair, scored on the **weakest** of the two calls.

Now the disease-specific part. In a proband with mosaic chromosome loss, let `f`
be the fraction of cells missing one homolog. For a heterozygous site on that
chromosome, the expected ALT read fraction is:

- homolog carrying **REF** lost: `1 / (2 − f)` → 0.667 at `f = 0.5`, **0.833 at
  `f = 0.8`** — above the 0.75 ceiling, −0.155.
- homolog carrying **ALT** lost: `(1 − f) / (2 − f)` → 0.333 at `f = 0.5`, 0.167
  at `f = 0.8` (mosaic window, −0.050), **0.091 at `f = 0.9`** — below the 0.10
  floor, −0.195.

Both tails are reachable in exactly the disease we are analysing, and they are
reachable *for the causal variants themselves*, because the causal gene's locus
is on a chromosome as liable to missegregate as any other.

`ASSUMPTION-MOSAIC-01` is right that mosaicism *dilutes* the alternate allele and
cannot concentrate it, and the code is right to give the above-band case its own
non-mosaic wording. But "not mosaicism" is not the same as "not signal": in a
mosaic-aneuploidy proband, allelic imbalance from **chromosome loss** is a
disease mechanism, and the pipeline currently reads it as a call-quality defect
and charges for it twice — once through `analytical_validity` and again through
the `low_quality_call` contradiction. `ASSUMPTION-SCORING-02` documents that
double count as deliberate and names its cost. Here it is at its most expensive
and least justified: the same observation is doing both jobs.

### 4.5 A causal allele we cannot represent at all — bounded, and smaller than TD-03 implies

TD-03 (no SV/CNV calling) and TD-13 (no repeat expansions) are real, but for
**Track 1 specifically** the submission format bounds the exposure. The answer
key stores `(chrom, pos, ref, alt)` with `ref`/`alt` compared as uppercased
strings. A large deletion, an inversion, or an STR expansion has no canonical
`ref`/`alt` we could ever guess even if we called it. The answer must therefore
be representable as two point coordinates — which means it is almost certainly
SNV/indel on both alleles, the class we can already see.

If one allele is nonetheless invisible, the ceiling is:

```
single {visible} at rank 1          rank_pts = 0 (strict) / 50 (lenient)   Fmax = 0.667
pair {visible, guess} at rank 1     rank_pts = 50 (both readings)          Fmax = 0.500
```

so roughly half of everything, and unrecoverable inside our scope. Given the
format argument, I rank this **below** items 4.1–4.4 on expected cost, which is a
deliberate disagreement with the ordering in `docs/tech-debt.md`.

### 4.6 The causal gene missing from a panel — not the risk it sounds like

`knowledge/public/gene_panel.tsv` is loaded nowhere in `src/mva/` outside
`src/mva/config.py`. Pairing is scoped by `variant.gene_symbols`, which is
derived purely from `VariantRecord.consequences`. So there is no panel gate to
fall outside of.

The real version of this risk is one line in the `pairing` module docstring:
*"candidates are gene-scoped, so a variant with no gene annotation forms no
candidate."* A variant the annotation source does not map to a gene is invisible
to pairing entirely — not down-ranked, absent — and with TD-01 the annotation
source is a local TSV. Annotation *coverage*, not panel membership, is the
binding constraint, and it collapses into item 4.1.

### 4.7 Ranked by expected cost

1. **4.1** — 0.40 of the weight is constant on a real case; ties broken by chromosome number.
2. **4.2** — the answer can be deleted before scoring, invisibly, with no bound on how often.
3. **4.3** — a systematic +0.051 handed to every rare homozygous call over the answer class.
4. **4.4** — up to −0.195 charged to the answer for the disease's own signature.
5. **1.4 / §3** — EPCR ties: −0.33 F-max and possibly −50 rank points, already occurring in the golden output.
6. **4.5** — bounded at roughly half credit, and argued down by the submission format.

---

## 5. Recommended changes, ranked by expected score per unit of risk

GP-32 requires a decision record, a test and an explicit before/after comparison
for any weight or threshold change. Each item states whether it trips that.

### 1. Gate on `gene_pair_cap_truncated` before submitting — *do this first*

**Change.** Surface the count of genes where `generate_pairs` truncated, and make
a non-zero count on the real case a blocking condition for a submitted row. The
flag already exists and is already propagated to every surviving candidate for
that gene; nothing computes with it.

**Expected effect.** No direct score change. It converts §4.2 from an unknown
into a known — a single boolean that tells us whether the answer could have been
deleted before scoring. If it fires, items 2 and 3 below become mandatory; if it
never fires, §4.2 is closed and we stop worrying about it.

**Cost.** Very low: a counter through the run report plus one assertion.
**Risk of making things worse.** None; it computes nothing new.
**GP-32.** No — no weight or threshold moves.

### 2. Guarantee strictly-decreasing, ≥0.01-separated EPCR

**Change.** A post-pass in `build_submission_rows` enforcing
`epcr[i] ≤ epcr[i−1] − 0.01`, floored at `EPCR_FLOOR`, plus a rule in
`validate_submission` rejecting any duplicate EPCR. `composite_to_epcr` stays
affine and untouched.

**Expected effect.** Removes a live defect (four tied EPCRs in the golden output
today) worth up to −0.333 F-max, and up to a further −50 rank points if the tie
lands adjacent to the answer and file order goes against us. Also insures against
a grid-based threshold sweep (§6 item 2), where our bunched composites would lose
0.33–0.71 F-max at 1–3 dp grid resolution.

**Cost.** Roughly ten lines and one test.
**Risk of making things worse.** Effectively zero — the pass is order-preserving
by construction, so no candidate can overtake another and no rank can change.
**GP-32.** Yes, formally: it changes an emitted numeric value and will move the
golden submission bytes. The before/after is mechanical — same order, ten
distinct values instead of eight.

### 3. Bound the pairing input, and stop truncating by coordinate

**Change.** Two parts. (a) Feed `generate_pairs` the output of
`select_candidate_variants` — already written, documented and tested, and already
carrying ADR 0010's `min_allele_number` guard — so that a biallelic recessive
hypothesis is built only from alleles below `max_plausible_recessive`. (b) Order
the within-gene candidate list by a plausibility proxy (worst impact severity
descending, then maximum population AF ascending) before applying
`max_pairs_per_gene`, and raise the cap now that the input is bounded.

**Expected effect.** Directly addresses the highest-severity fixable failure in
§4.2: turns a potential structural zero into a scored candidate. Part (b) alone
removes the "chr-position decides which hypotheses exist" behaviour.

**Cost.** Moderate. Part (a) is wiring; part (b) is a new ordering function and
its test.
**Risk of making things worse.** Real and specific: part (a) makes a frequency
*selection* into the boundary of the hypothesis space, which is the exact
asymmetry ADR 0005 exists to prevent. A founder allele above 1% in the proband's
ancestry would be excluded from pairing. Mitigations: `min_allele_number` is
already applied inside `select_candidate_variants`; the full flagged set is still
carried forward for the dossier (GP-13); and the count of variants excluded from
pairing must be reported. If that risk is judged unacceptable, take part (b) only
— it is pure improvement, since it changes *which* candidates a fixed cap keeps
rather than which exist.
**GP-32.** Yes for both parts, and it also needs an ADR: it touches the ADR 0005
filter/rank boundary directly and that boundary is load-bearing.

### 4. Stop double-charging allele-fraction deviation in a mosaic-aneuploidy case

**Change.** Narrowest defensible version: keep every `analytical_validity`
multiplier in `_allele_balance_penalty` exactly as it is, but stop
`_quality_flags` raising `low_quality_call` on **allele-fraction grounds alone**,
so the 0.30 `CONTRADICTION_LOW_QUALITY` term is not added on top. Depth, GQ,
caller FILTER and the ingestion `UNTRUSTED_CALL_FLAGS` continue to raise it.

**Expected effect.** Recovers `0.25 × 0.30 = 0.075` composite for candidates
whose only defect is an off-band allele fraction — i.e. exactly the signature of
mosaic chromosome loss (§4.4). The above-band case goes from −0.155 to −0.080 and
the below-floor case from −0.195 to −0.120. The doubt about the call is still
priced, once, through the 0.20-weighted validity component.

**Expected effect on ranking, honestly.** Unknown in sign for the submission as a
whole: it lifts every off-band candidate, not only the right ones. What makes it
defensible is that `ASSUMPTION-SCORING-02` already names this double count as a
deliberate cost, and this is the one context where the two counts are not
answering different questions — both are reading the same allele fraction.

**Cost.** Low.
**Risk of making things worse.** Moderate. Reference-allele dropout and
mis-genotyped homozygotes are genuine artifacts and this makes them cheaper. Do
not extend it to the below-floor band without CNV evidence (TD-03) that the site
sits in an aneuploid region.
**GP-32.** Yes — plus an amendment to `ASSUMPTION-MOSAIC-01` and
`ASSUMPTION-SCORING-02`, since both currently describe the behaviour being
changed.

### 5. Price the unverified-zygosity open question symmetrically

**Change.** Apply a discount to `_INHERITANCE_HOMOZYGOUS` (0.90) for a
homozygous single-variant candidate where true homozygosity is unconfirmed —
which, without CNV calling, is *every* such candidate. The pipeline already emits
the open question `::true-homozygosity` and already treats an unverified
haplotype assumption as worth a 0.55× multiplier on the pair side.

**Expected effect.** A discount to 0.75 closes `0.16 × 0.15 = 0.024` of the
0.051 gap in §4.3. It does not eliminate it, and it should not — a homozygous
call really does account for both gene copies with one observation.

**Cost.** Low: one constant and a conditional.
**Risk of making things worse.** This is the recommendation closest to the
ASSUMPTION-SCORING-03 line, and the argument has to stay on the right side of it.
The justification must be *"our evidence model prices unverified phase but not
unverified zygosity, and both are unverified for the same reason"* — never
*"the answer is a compound het"*. If that argument does not survive review,
**drop this item**; it is the smallest effect in the list.
**GP-32.** Yes, unambiguously.

### 6. Spend spare slots on same-gene fan-out, not on a tenth gene

**Change.** In `build_submission_rows`, prefer — among candidates whose
composites are within rounding distance of one another — those introducing fewer
variants not already emitted by a higher-ranked row.

**Expected effect.** Strictly non-negative: such rows contribute zero to `k` at
every threshold at or above their own (§2.4) while adding a shot at the exact
frozenset match.

**Cost.** Low.
**Risk of making things worse.** Only via §6 item 3 — if the scorer
short-circuits on a partial match, a partially-overlapping fan-out row above the
full match halves the credit. Bound it by keeping fan-out rows within one rank
tier. **Verify item 3 before implementing this.**
**GP-32.** No weight moves, but it changes the submission's composition; a short
decision record is cheap insurance.

### 7. Replace the synthetic knowledge and annotation layers (TD-01, TD-04)

**Change.** Real VEP/gnomAD/ClinVar annotation; a real HPO release with ancestor
closure; a real gene-phenotype table covering chromosomal-instability genes.

**Expected effect.** The largest of anything in this list — it reactivates 0.40
of the scoring weight and stops chromosome number from breaking ties (§4.1). It
is also the only change that makes `rarity`, `phenotype_similarity` and
`mechanistic_relevance` mean anything on a real case.

**Cost.** Was the highest in this list; may now be among the lowest. The adapter
modules named in §4.1 already exist in the tree — the remaining work is wiring
them through `src/mva/orchestrator.py` and replacing the synthetic knowledge
tables, not building them. **Re-price this item against the working tree before
scheduling anything else in this list.**
**Risk of making things worse.** Low in kind, high in schedule: a half-wired
annotation source that silently drops gene symbols is worse than an honestly
neutral one, because it converts item 4.6's "invisible to pairing" from a
theoretical to an actual loss.
**GP-32.** No weight change, but GP-20 maturity grades and the TD-01/TD-04 rows
must move with it.

### What I recommend *against*

- **Raising `phase_weights.unknown`.** The only reason to do it is that we know
  the answer is a compound het. That is fitting the key (ASSUMPTION-SCORING-03).
  Item 5 reaches part of the same place through an argument that survives without
  knowing the answer.
- **Using the six attempts to search configurations.** "Best of six counts" makes
  this tempting and it is leaderboard fitting with extra steps. The legitimate
  use is **mechanical verification**: spend attempt 1 confirming the submission
  parses, the `chr` prefix assumption holds, and both metrics come back
  non-degenerate. That verifies the vendored contract rather than fitting the key,
  and it is the highest-value non-scientific use of an attempt.

---

## 6. Re-verify against the live Space before submitting

The contract file says to. These are the specific items where the arithmetic
above would change if the implementation differs from the prose. Ordered by how
much of §1–§3 they would invalidate.

1. **F-max comparison operator and threshold set.** Is a row included at
   `epcr >= t` or `epcr > t`, and does the sweep include a threshold at or below
   the minimum emitted EPCR? Under `>` with thresholds drawn from emitted values,
   the lowest-EPCR row can never enter the predicted set and full recall through
   the last row is unreachable. §1.2's monotonicity result depends on this.
2. **Are thresholds our unique EPCR values, or a fixed grid?** A fixed grid makes
   EPCR *magnitude* score-relevant and turns the invariance theorem in §0 false.
   This is the single assumption most load-bearing for §3, and it is exactly what
   §5 item 2 insures against either way.
3. **Does match scanning short-circuit on a partial match?** "The first (best)
   full match wins" admits both readings. If it stops at the first row with any
   intersection, a partial-overlap row ranked above a full match caps us at
   `0.5 × rank_points`. This gates §5 item 6 and bounds the fan-out strategy.
4. **Does partial credit require the *row* to be a compound-het row, or only the
   ground truth?** Decides whether a single-variant row can score at all. §2's
   conclusion (always pair) holds under both readings, so this is informational
   rather than blocking — but it changes the magnitude of what §2.3 buys.
5. **How are rank points and F-max combined?** Two columns, a weighted sum, or a
   primary metric with a tiebreak? Every trade-off in §2 and §5 assumes rank
   points dominate; if F-max carries most of the weight, the single-variant
   hedging arithmetic in §2.3 flips.
6. **Does the ground truth for `PROBAND01` contain exactly two variants?** Every
   recall denominator here assumes `|G| = 2`. Secondary or incidental findings in
   the key would change all of §1.3.
7. **Is `finding_type` used in scoring?** If `secondary` rows are excluded from
   rank matching or from the F-max union, that would matter — though not to us
   today: `INCIDENTAL_FLAG` is defined in `src/mva/reporting/track1.py` and set
   nowhere, so every row we emit is `primary`.
8. **The `chr` prefix in the gold key.** The contract already flags this as the
   single most dangerous detail. Note that it is inferred from the CSV template
   and the repo's local fallback ground truth, **not** from the gated key itself.
   `_validate_contig` protects us against emitting a bare contig; nothing
   protects us if the key holds bare contigs and we emit prefixed ones.
9. **`PROBAND01` as the only accepted id, and the 10-row limit.** Cheap to
   re-check, and both hard-fail the whole submission.
