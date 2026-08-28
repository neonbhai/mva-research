# ADR 0019 — A named selection stage in front of pairing, and unknown-frequency variants are retained

**Status:** accepted · **Date:** 2026-08-28

## Context

`apply_hard_filters` removes only invalid or impossible records: wrong build,
non-canonical contig, placeholder allele, hom-ref, no-call. That is GP-13 and
ADR 0005 working exactly as intended — commonness, poor quality and
predicted-benign consequence are soft flags that down-rank, never deletions.

The consequence at whole-genome scale was measured, not assumed
(`docs/scale-report.md` §5). On a 4.5 M-record callset:

- **2,020,500 variants reach pairing** — 44.9% of the callset falls inside a MANE
  gene locus, measured with `ManeGeneIndex.genes_at` over the phantom's
  positions. Mean load: **103.6 variants per gene** across 19,500 genes.
- `generate_pair_candidates` would build **7,402,500 candidate objects**, keep
  390,000 after `max_pairs_per_gene=20`, and **discard 94.7% of them**.
- **Every one of the 19,500 genes exceeds the cap.** So every candidate the
  pipeline emits carries `gene_pair_cap_truncated`.

A flag raised on 100% of output distinguishes nothing. Its stated purpose — "a
silently shortened hypothesis list is indistinguishable from a complete one"
(ADR 0013) — is defeated when the list is always shortened. And the ordering
doing the shortening is a truncation heuristic that was designed as a backstop
for hyper-variable genes, not as the primary filter for two million variants.

The cap was already performing selection. It was doing it silently, on the wrong
input, by the wrong rule, and announcing it on every result.

## Decision

**1. Add an explicit selection stage between annotation and pairing.**
`mva.prioritization.selection` keeps only variants that are plausibly rare and
plausibly coding or splice-relevant. It is a named stage with recorded
provenance, not a hidden filter:

- Every dropped variant is counted under one of five reason codes and every
  retained variant under one of five keep codes. The two sets partition their
  side of the decision, so the counts sum to the input and a reader can check
  the arithmetic.
- The counts, the thresholds that produced them and the stage's warnings are
  written as `selection/selection_report.json` (`ArtifactKind.SELECTION_REPORT`,
  classified DERIVED_SAFE — counts only, never a variant ID) and registered with
  `ArtifactProvenance` like any other artifact.
- Aggregate `EvidenceItem`s carry the same counts into the ledger, so a report
  can cite what selection did (GP-10).
- The stage never mutates a record. A selected variant is the same object that
  arrived, so nothing downstream changes shape because selection ran.
- Thresholds live in configuration, not in constants inside the filter (GP-32).

**2. Variants with NO frequency data are RETAINED.**
`retain_unknown_frequency` defaults to True, and so does the equivalent for an
observation that exists but is under-powered (every reporting population below
`min_allele_number`, ADR 0010).

**3. Variants whose impact was never assessed are RETAINED.**
`retain_unassessed_impact` defaults to True. `impact is None` means NOT ASSESSED
(ADR 0016), not MODIFIER.

**4. The selection cut-point is separate from, and looser than, the ranking
cut-point.** `SelectionThresholds.max_population_frequency` defaults to 0.02
against `FrequencyThresholds.max_plausible_recessive` of 0.01. A configuration
that inverts that relationship is reported as a warning.

**5. A curated pathogenic assertion overrides both gates.** A known pathogenic
allele that is common in some cohort is still the answer.

## Why unknown-frequency variants are retained

This is the decision that decides whether the pipeline can find the answer at
all, so the reasoning is written out rather than asserted.

**Absence of a gnomAD record is not AF = 0 (GP-14).** The pipeline already
refuses that conflation twice — `annotate_variants` emits an explicit NEUTRAL
"frequency data is UNAVAILABLE" evidence item, and `FrequencyThresholds`
`absent_frequency_score = 0.5` scores the gap mid-range. A selection stage that
dropped unknown-frequency variants would reintroduce the same error with the
opposite sign: treating "we have never seen it" as "it is common enough to
delete".

**A genuinely novel pathogenic variant is precisely the target.** A severe
paediatric recessive disorder is not caused by an allele gnomAD has catalogued
at 3%. The most likely shape of the answer is a variant seen in no reference
cohort — which, under a drop-the-unknown rule, is the one variant guaranteed to
be deleted.

**Absence correlates with the wrong things.** A variant missing from gnomAD is
most often missing because the site is poorly covered, because the ancestry is
under-represented in the panel, because an indel's representation did not join,
or because a liftover failed. `docs/next-actions.md` records the measurement:
on the real gnomAD chr21 shard, **0 of 1,029 right-shifted indel spellings
joined without a reference against 989 with one**, and thirty of those are
variants gnomAD calls *common*. Under a drop-the-unknown rule, an indel
normalisation bug becomes indistinguishable from a scientific finding, and the
loss is concentrated in the patients whose ancestries the references serve
worst.

**The cost of the opposite error is bounded and visible.** Retaining unknowns
inflates the selected set — measured on the phantom, unknown-frequency variants
are 38% of the input and the selected set is still 1.6% of the callset, a few
hundred candidates. They are then *ranked*, where the gap costs them: rarity
scores 0.5, not 1.0. A false positive is a candidate at rank 40. A false
negative is a deleted candidate, which is invisible and unfalsifiable.

The switch exists (`retain_unknown_frequency = False`) because the policy should
be arguable, not because it should be used. Setting it False produces a warning
naming the failure mode in as many words.

## Consequences

**GP-13 is bent here, deliberately, and that is stated.** This is the only stage
entitled to delete a valid record. Three things bound the damage:

1. The selection cut-point is looser than the ranking cut-point, so nothing
   deleted here was still rankable above zero on rarity.
2. Every unknown — no frequency data, no assessed impact — is retained, so the
   stage only ever deletes on a *positive* finding (a measured common frequency,
   or an assessed non-coding consequence), never on an absence.
3. The stage can be switched off (`enabled: false`) and still runs, still counts
   and still reports, so "what did selection cost me" is answerable without
   rewiring anything.

**GP-19 is partially satisfied, and the shortfall is named.** GP-19 wants failed
candidates persisted with their reasons. Four million individual rejection
records is not something this pipeline can hold — the ledger spill exists
because ten million evidence items already do not fit. The persisted form is the
per-reason counts, in the report artifact and in aggregate evidence items. Per-
variant rejection records for the dropped set remain available by re-running with
`enabled: false`, which is not the same thing and is not claimed to be.

**`synonymous_variant` is not in the coding term set.** It is coding but does not
alter the protein, and including it would roughly double the surviving set. A
synonymous call that *also* carries `splice_region_variant` is retained, because
matching is over the whole term list across every transcript — which is exactly
the case `prioritization.filters` cites when it declines to use predicted
consequence at all.

**Variants with no gene assignment are dropped, and counted.** Roughly 55% of a
whole-genome callset lands outside any MANE gene locus. Gene-scoped pairing
already ignores them, invisibly. Dropping them under
`dropped_no_gene_assignment` changes nothing downstream and makes the largest
single loss in the pipeline a number in an artifact. It is not a claim that they
are benign.

**Selection depends on annotation having run.** With no consequence adapter
bound, every variant has zero consequences and `drop_without_gene_assignment`
would delete the entire callset. The stage is therefore only correct downstream
of a working gene-assignment path (`docs/next-actions.md` item 1), and the
report makes the failure loud rather than silent: `dropped_no_gene_assignment`
equal to the input count is unmistakable.

## Alternatives rejected

**Leave the cap to do it.** Measured: 94.7% of hypotheses discarded, every gene
flagged, ranking decided by a truncation heuristic. The status quo is not a
neutral option; it is an unnamed filter with worse behaviour.

**Filter inside `apply_hard_filters`.** That module's entire purpose is that it
removes only invalid records, and its docstring says so at length. Adding a
frequency rule there would make GP-13 unenforceable by reading the code.

**Drop unknown-frequency variants to shrink the set further.** Rejected above.
It optimises the number this stage reports at the cost of the result the
pipeline exists to produce.

**Emit one `EvidenceItem` per dropped variant.** Four million items of the exact
kind the ledger spill exists to survive, describing candidates nobody will read.
Counts by reason carry the same information at four orders of magnitude less
cost.
