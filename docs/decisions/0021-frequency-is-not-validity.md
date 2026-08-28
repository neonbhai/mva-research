# ADR 0021 — Phenotype frequency is not gene-disease validity; they get separate columns

**Status:** accepted
**Date:** 2026-08-28
**Touches:** GP-13, GP-14, GP-19, GP-20, GP-30, GP-31, GP-32, TD-17, TD-20, ASSUMPTION-PHENOTYPE-07

## Context

`knowledge/real/gene_phenotype.tsv` (221,789 associations; 221,813 lines including
the header block) could not be loaded at all.
`GenePhenotypeIndex.from_tsv` raised on the first data row:

```
IngestionError: Unknown association_strength '' for A2ML1; allowed values are:
definitive, moderate, strong, supporting.
```

TD-17 recorded this as "31,169 empty rows". That understated it. The column named
`association_strength` did not hold association strength *at all*. It held HPO's
phenotype-**frequency** vocabulary, verbatim:

| value in `association_strength` | rows | what it actually meant |
|---|---:|---|
| `HP:0040283` | 47,097 | Occasional — present in 5–29% of cases |
| `HP:0040282` | 35,110 | Frequent — 30–79% of cases |
| `HP:0040281` | 19,201 | Very frequent — 80–99% of cases |
| `HP:0040284` | 6,996 | Very rare — 1–4% of cases |
| `HP:0040280` | 337 | Obligate — 100% of cases |
| `n/m` fractions | 81,778 | n of m annotated cases had the feature |
| percentages | 102 | recorded directly in the source |
| *(empty)* | 31,169 | HPO recorded `-`: frequency not stated |

Neither module was confused about its own half. `tools/build_knowledge/gene_phenotype.py`
said in its docstring that it was writing frequency and that frequency and strength
"are different scientific concepts". `src/mva/phenotype/hpo.py` validated the same
column as curated gene-disease evidence strength and weighted it 1.0/0.8/0.6/0.4.
**Two modules agreed on a column name and disagreed about what it meant.** The
column name was the entire interface, and it was ambiguous.

The reader failing closed was correct behaviour and is preserved. A lenient reader
would have weighted "occurs in 5–29% of cases" as "definitively curated", which is
the exact class of fabricated precision this repository exists to prevent.

Two quantities were in play, and only one of them was in the file:

* **Phenotype frequency** — how often a feature occurs *among cases of a disease*.
  Real, useful, and what HPO's `genes_to_phenotype.txt` actually carries.
* **Gene-disease clinical validity** — how confident an expert panel is that
  variation in a gene causes disease *at all*. Not in `genes_to_phenotype.txt`
  anywhere. It lives in ClinGen Gene-Disease Validity and EBI Gene2Phenotype
  (DDG2P), both already on disk and already the source of `evidence_tier` in
  `knowledge/disease/mva_panel.tsv`.

## Decision

### 1. Each quantity gets its own column, its own vocabulary and its own parser

`knowledge/*/gene_phenotype.tsv` columns are now:

```
gene_symbol  hpo_id  label  association_strength  association_strength_source
             hpo_frequency  source  version
```

* **`hpo_frequency`** carries HPO's own token verbatim — one of the five
  `HP:00402xx` terms, an `n/m` fraction, or a percentage — parsed by
  `mva.phenotype.hpo.parse_hpo_frequency` into a typed `HpoFrequency` with the
  proportion bounds `hp.obo` states. Empty means HPO recorded `-`: **unmeasured,
  not "never occurs"** (GP-14). `HP:0040285` (Excluded, 0% of cases) is a *negated*
  annotation, not an association: the generator drops such rows and prints the
  count in the table header (0 in this release, measured not assumed), and the
  reader refuses one by name.
* **`association_strength`** carries curated clinical validity from ClinGen and
  DDG2P, in those sources' own vocabulary, case-folded (ClinGen's `Definitive` and
  DDG2P's `definitive` are one concept; folding case is not a semantic remap).
* **`association_strength_source`** names which panel made the call. `source`
  stays `HPO`, which is where the gene→term *annotation* came from. One provenance
  field covering both would attribute an expert-panel judgement to HPO (GP-31).

Both new columns are **optional** to the reader, so one reader loads the real table
and the synthetic demo table (`knowledge/public/gene_phenotype.tsv`, which has
neither). A file without `hpo_frequency` is not a file whose phenotypes never
occur; it is a file that does not state how often they do.

### 2. `association_strength` is gene-level, and the table says so

ClinGen and DDG2P curate the **(gene, disease)** pair; this column reduces to the
**gene**, taking the strongest classification either panel made. A gene curated
`definitive` for one disease and `limited` for another reads `definitive` on every
row, including phenotypes annotated via the weaker disease. That over-statement is
real and is stated in the table header, in the generator docstring and here.

The alternative — joining each HPO annotation to the curation for *its own*
disease — was measured, not assumed: it reaches **75,379 of 275,046** candidate
source rows (27%). HPO annotates against OMIM and ORPHA; ClinGen curates against
MONDO; no crosswalk exists in the downloaded resources. A column that is precise
for a quarter of rows and absent for the rest is not more honest than a gene-level
one whose granularity is documented — it is the same approximation with the
documentation removed. Recovering (gene, disease) granularity needs a MONDO↔OMIM
crosswalk and is recorded as tech debt.

### 3. The weight vocabulary becomes the real one (GP-32)

`STRENGTH_WEIGHTS` was a four-value scale the synthetic demo table invented. It is
now the union of the two real curation vocabularies plus the demo's own:

| classification | weight | before | source |
|---|---:|---:|---|
| `definitive` | 1.0 | 1.0 | ClinGen, DDG2P |
| `strong` | 0.8 | 0.8 | ClinGen, DDG2P |
| `moderate` | 0.6 | 0.6 | ClinGen, DDG2P |
| `limited` | 0.4 | — | ClinGen, DDG2P |
| `supporting` | 0.4 | 0.4 | synthetic demo table only |
| `disputed` | 0.2 | — | ClinGen, DDG2P |
| `refuted` | 0.1 | — | ClinGen, DDG2P |
| `no known disease relationship` | 0.1 | — | ClinGen |

**Before/after (GP-32):** every value that existed before is unchanged, so the
synthetic case is unmoved. Measured, not asserted: `tests/golden/expected_ranking.tsv`
and `tests/golden/expected_drug_outcomes.tsv` are byte-identical and their hash
locks in `tests/golden/test_locked_files.py` were not touched. The new entries are
reachable only from real data, which no golden expectation covers yet.

The floor is 0.1, not 0.0. A refuted gene-disease claim is still a curated
statement and the HPO annotation hanging off it is still a real annotation; zero
would delete a candidate that GP-13 and GP-19 say must be down-ranked and kept.
An unlisted classification still raises: a new ClinGen tier is a change in the
source's curation model and must be placed in the ladder deliberately.

### 4. Absence is representable, and what it is worth is said out loud

`GeneAssociation.association_strength` is now `str | None`. `None` means **no
source classifies this gene** — which is not `limited`, and not `refuted` (GP-14).

The two tempting defaults are both wrong and both silent:

* **0.0 erases** a real HPO annotation. The association drops out of the score's
  denominator entirely and the gene becomes indistinguishable from one this
  pipeline knows nothing about.
* **1.0 invents** a definitive expert classification that nobody made.

So the choice is made visibly and cannot be made by omission:

* `GeneAssociation.weight` returns `float | None`. Under pyright strict, `None` is
  not a number a caller can accidentally sum — every consumer must handle it.
* `GeneAssociation.weight_or(uncurated: float) -> float` has **no default
  argument**. Every call site names the value:
  `assoc.weight_or(UNCURATED_ASSOCIATION_WEIGHT)`, which is greppable.
* `UNCURATED_ASSOCIATION_WEIGHT = 0.5`, deliberately **off the curated ladder**.
  Reading 0.5 back out of an evidence item's `numeric_value` is itself the signal
  that nobody classified the gene. `test_uncurated_weight_is_not_mistakable_for_a_curated_one`
  enforces that it never collides with a curated weight.
* Uncurated associations report `EvidenceStrength.INSUFFICIENT` and say
  "an uncurated (no source classifies this gene's disease validity)" in their claim
  text, rather than borrowing a tier.

Note that the score's specificity component is a *ratio* of association weights, so
a gene-level constant strength cancels within a gene: the choice above changes what
a reader sees in the evidence store, not the ranking. That is a property of this
design, not an accident, and it is why the within-gene discrimination the old
per-term strengths pretended to provide has to come from `hpo_frequency` and
information content instead.

## Consequences

* `knowledge/real/gene_phenotype.tsv` loads: **221,789 rows across 3,634 genes**
  (`test_real_gene_phenotype_table_loads_through_the_reader`). That is the same row
  count as before — the "221,813" in the bug report is the file's line count, which
  includes 23 comment lines and the header. **No association was dropped by this
  change**: the 31,169 rows with no stated frequency are all still there, now with an
  empty `hpo_frequency` and a real `association_strength`.
* Every row in this build carries a curated strength, because the table is
  restricted to genes curated by ClinGen or DDG2P. The *absence* path is
  nonetheless live, tested and reachable — `restrict_to_genes` and
  `curated_strengths` are deliberately separate arguments to
  `parse_gene_to_phenotype`, so widening the restriction cannot silently invent
  curation for the genes it adds.
* Usefulness, not just loadability. Scoring `BUB1B` (MVA1) against the 8-term
  public HPO profile HP:0002859, HP:0000121, HP:0004322, HP:0001508, HP:0003202,
  HP:0001622, HP:0001518, HP:0200067, with real `PhenotypeSemantics` (hp.obo +
  phenotype.hpoa, release `hp/releases/2026-06-23`):

  ```
  score 0.871096   mode ontology_semantic   annotation_status annotated
  specificity 1.0  coverage 0.742192        symmetric best-match average 0.457203
  matched 5 / contradicted 0 / unassessed 107 / unexplained 5
  ```

  | matched term | label | strength (source) | hpo_frequency |
  |---|---|---|---|
  | HP:0001510 | Growth delay | definitive (ClinGen) | HP:0040283 |
  | HP:0001518 | Small for gestational age | definitive (ClinGen) | 9/9 |
  | HP:0002664 | Neoplasm | definitive (ClinGen) | HP:0040283 |
  | HP:0002859 | Rhabdomyosarcoma | definitive (ClinGen) | HP:0040283 |
  | HP:0004322 | Short stature | definitive (ClinGen) | HP:0040281 |

  Two of the five (HP:0001510, HP:0002664) are reached only by upward `is_a`
  entailment, so this exercises the real ontology path, not string equality.
* `gene_panel.tsv` and `gene_disease.tsv` regenerate **byte-identical** from the
  same generator run (GP-30); only `gene_phenotype.tsv` changed.
* The mutation was executed, not asserted. Reintroducing the conflation in the
  generator (`association_strength=frequency`) and rebuilding the table fails four
  named tests: `test_gene_phenotype_frequency_never_lands_in_the_strength_column`,
  `test_real_gene_phenotype_association_strength_is_curated_validity`,
  `test_real_gene_phenotype_table_loads_through_the_reader` and
  `test_gene_phenotype_strength_is_absent_when_no_source_curates_the_gene` — with
  the reader's message naming the confusion rather than reporting an unknown value.
* The conflation cannot come back quietly. A frequency token offered as an
  `association_strength` is refused with a message naming the confusion
  (`test_hpo_frequency_token_is_refused_as_an_association_strength`), the committed
  table is asserted to hold curated validity and no frequency-shaped value
  (`test_real_gene_phenotype_association_strength_is_curated_validity`), and the
  generator is asserted to write each quantity to its own column
  (`test_gene_phenotype_frequency_never_lands_in_the_strength_column`).

## Alternatives rejected

**Loosen the reader to accept anything.** This is what made the bug expensive in
the first place. The reader failing closed is the only reason the conflation was
found rather than shipped as a number.

**Drop the rows with no frequency.** 31,169 real gene-phenotype annotations
deleted because a *different* column was missing. GP-13: hard filters remove
invalid records, not inconvenient ones.

**Map HPO frequency onto the strength scale** (e.g. Obligate→definitive). This is
the conflation with extra steps. "Occurs in 100% of patients who have this disease"
and "we are certain this gene causes the disease" are independent claims; a gene
can be refuted and still have obligate annotations hanging off the disease it was
refuted for.

**Leave `association_strength` empty for every HPO-derived row.** Honest, but it
throws away curation that is on disk and already trusted elsewhere in the repo
(`mva_panel.tsv`), and it would make `UNCURATED_ASSOCIATION_WEIGHT` the weight of
every real association — a constant standing in for data we have.

**Keep one column and rename it `hpo_frequency` only.** Then the reader has no
strength at all and `STRENGTH_WEIGHTS` becomes dead code, while ClinGen and DDG2P
sit unused on disk.
