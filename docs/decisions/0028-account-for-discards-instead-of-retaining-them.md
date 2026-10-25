# ADR 0028 — Annotation is staged, and a dropped variant is accounted for rather than retained

**Status:** accepted · **Date:** 2026-08-29 · **Amends:** GP-13 · **Extends:** ADR 0005, ADR 0019

## Context

A real whole-genome run was killed after 7 h 30 m. It had annotated 38.5% of
4,962,060 variants and was projecting **~19.5 hours and ~95 GB of intermediates
onto a 60 GB volume**. It could not have finished. Nothing in the pipeline said
so; the only way to know the stage was alive was to sample its stack.

An ad-hoc script did the same callset end to end in ~14 minutes. Reading it and
measuring against it produced the number this decision turns on:

| | |
|---|---|
| Records in the ad-hoc slimmed callset | 5,012,204 |
| Records our pipeline annotates | 4,962,060 |
| SnpEff over all 5.01 M | **~14 minutes** |
| Survivors of a HIGH/MODERATE impact filter | **14,228 — 0.284%** |

**The difference was not filtering before annotation.** The ad-hoc prefilter drops
hom-ref calls, missing genotypes and no-alt records — all three already in
`HARD_FILTER_REASONS`, all annotation-independent — and the callset was already
non-ref, so it removed almost nothing.

The difference is that we did per-variant gnomAD lookups, ClinVar lookups, ~6 KB
of JSON and ~3 evidence rows for **all 4,962,060** variants, when 14,228 of them
can possibly matter. SnpEff is cheap and 350x selective. We ran it and threw its
selectivity away.

Two things confirmed the diagnosis while this was being built, and both narrow it
rather than widen it:

* **SnpEff is not the slow part.** Run genome-wide over the same 5,012,204
  records it is a tens-of-minutes job — the ad-hoc script's ~14 minutes — against
  our pipeline's 19.5-hour projection on the same tool and the same input.
  Independent consequence-filtering tools land in the same place on this callset,
  **14–18 k survivors**, and our measured 14,228 sits inside that. The bottleneck
  is **population frequency, not consequence annotation**.
* An audit of `src/mva/annotation/snpeff_local.py` confirms the SnpEff path is
  **not** the defect. It launches the JVM 20 times for a 4.96 M callset, not 993;
  it streams SnpEff's stdout rather than buffering it; it passes `-noStats`,
  `-noLog` and `-nodownload`. What remains there is small and is recorded as
  TD-26.

That matches our own stall exactly: 8% CPU during 4.96 M random tabix seeks
across 22 gnomAD shards.

### Why GP-13 was in the way

GP-13 said hard filters remove only invalid records and **everything else is
flagged and down-ranked, never deleted**. ADR 0019 had already bent it once, for
this same reason at a smaller scale, and said so plainly: the selection stage is
"the one stage entitled to delete a valid record", and its persisted form is
per-reason counts because "four million individual rejection records is not
something this pipeline can hold". It also named the shortfall — GP-19 wants
failed candidates persisted with their reasons, and counts are not that.

The principle protects against a candidate becoming **invisible and
unfalsifiable**. It does not require retention to do that. A one-line record
naming the variant and the rule that dropped it satisfies the protection at
roughly 75 bytes instead of roughly 13 KB.

## Decision

**1. Annotation is staged. The consequence pass is whole-callset; the expensive
half is partitioned.**

`mva.annotation.triage` decides, from a variant's consequence annotations alone,
whether it gets a gnomAD lookup, a full `EvidenceItem` set and a row in
`variants/annotated.json`. The verdict has to be reachable before those calls,
because avoiding them is the entire saving.

**2. The bar is generous by default, and the default is stated.**

A variant is kept when any transcript annotation says any of: impact at least
`min_impact` (default **LOW**); any coding-altering or splice-relevant SO term,
whatever the impact; a SpliceAI delta at or above 0.2; impact never assessed on
any transcript (ADR 0016); or a curated clinical assertion. A variant with **no
consequence annotation at all** is dropped, and that one is provable rather than
heuristic — pairing is gene-scoped, so a variant with no gene symbol forms no
candidate pair under any configuration.

At 14,228 survivors on HIGH/MODERATE there is room to be wrong, so the default is
looser than the measurement requires on three axes at once. Impact class is an
ordinal bucket produced by a predictor, and the class this pipeline can least
afford to lose — a synonymous variant that disrupts splicing — is exactly the one
a tight impact bar loses. `prioritization.filters` declines to use predicted
consequence at all, for the same reason, in as many words.

**Measured survivor counts.** On the scale phantom at 200,000 records
(`tools/scale/stage_harness.py triage`), the default bar retained **9,078 of
200,000 — 4.54%, a 22.0x record reduction**. The phantom's own rates are the
reason it is not 0.284%: 44.9% of its records fall in a gene and 8% of those carry
an impact of `None`, which the bar retains by policy. On the real callset's
measured impact distribution the same bar sits between the 14,228 HIGH/MODERATE
survivors and the ~2.02 M genic records — nearer the former, because the LOW band
is small and MODIFIER is most of what is genic.

**3. Every dropped variant gets one row in a discard ledger, at every gate.**

`discards/discard_ledger.tsv`: `variant_id`, `stage`, `rule`, `value`. Measured at
**64.3 bytes per row**, so 5 M rows cost ~322 MB rather than the ~65 GB the same
information as `EvidenceItem`s would.

Four stages write to it — `ingest`, `triage`, `hard_filter`, `selection` — and
`ingest` is there deliberately. Discard auditing fails in a predictable place: it
is complete *downstream* of the first region gate and absent *at* it. Ours had
exactly that hole — the reader's `skipped_count` removed 116,485 decoy and
unplaced records on this callset, more than eight times the number of variants
that survive to the end, as a bare integer in `qc_report.json`. The distinction
the ledger now keeps is the one that closes it: **"never extracted" and
"extracted then dropped" stay distinguishable.**

`mva discards explain <ledger> <variant_id>` answers "why isn't variant X in the
output?" by scanning the file. An empty result is reported as **ambiguous**, not
as good news: the ledger records removals, not a census, so absence from it means
either the variant survived or it was never in the input (GP-14).

**4. GP-13 is amended: the requirement is that a dropped record be ACCOUNTED FOR,
not that it be retained.** The discard ledger is the mechanism that discharges it,
and it is now the thing GP-13 points at.

**5. The bar can never be narrower than the selection stage behind it.**

`triage.narrow_to_selection` widens the configured bar until it cannot drop
anything `SelectionThresholds` would keep, and every narrowing is a run warning
rather than a silent adjustment. When `selection.enabled` is false — which is the
synthetic demo's configuration, ADR 0027 §7 — the bar collapses to the
gene-assignment test, which is the only drop that stays provable without a
selection stage behind it.

**6. The clinical lookup runs in front of the bar; the frequency lookup does not.**

This was wrong in the first implementation and the scale phantom caught it.
`selection.retain_pathogenic_clinical` lets a curated pathogenic assertion rescue
a variant both other gates would drop — and that rescue needs a ClinVar lookup.
With the lookup behind the bar, the selected set went **3,616 → 3,600**: sixteen
variants lost, exactly the ones the rescue existed for.

So ClinVar is asked about every record and its answer is an input to the verdict.
ClinVar is one ~100 MB tabix file; gnomAD is 22 multi-gigabyte shards, and gnomAD
stays behind the bar. The cheap half of the deferral was the one that was wrong.
`triage.clinical_before_bar: false` restores the hole and warns about it in the
run report, because someone may measure the cost as unaffordable and should be
able to say so deliberately.

**7. A stage that cannot possibly complete says so in its first minute.**

`mva.progress.StageMonitor` projects total time and bytes from the first batch
against free disk, prints the projection, and raises `InsufficientDiskError` when
the projected output does not fit. Fed the measured run's own shape it prints:

```
[annotate] pre-flight projection from 250,000 of 4,962,060 records (5.0% sampled)
[annotate]   measured : 18 rec/s, 19.0 kB/record, 4.8 GB written in 3h 55m
[annotate]   projected: 78h 00m total, 94.3 GB of output
[annotate]   disk     : 60.0 GB free, 54.0 GB budget after 10% headroom
```

and refuses. **This matters more than the speedup.** A pre-flight projection of
runtime and disk is an unusual thing for a pipeline to carry, because it only
earns its keep in a design that generates enough intermediate data for the
projection to matter — and this design does, which is exactly why it had to be
added here.

The same object emits throttled progress — records/s, count, ETA, bytes — every
30 seconds. It reads the clock once per 4,096 records so the meter is not a
measurable share of the stage it measures, and it writes to stderr so it cannot
contaminate a stdout that is a CSV.

## Consequences

**What is given up.** `variants/annotated.json` is no longer a complete record of
the callset. It holds what cleared the bar, and everything else is a line in the
discard ledger with no annotations attached. A reviewer who wants the full
annotated set re-runs with `triage.enabled: false`, which is the same recourse
ADR 0019 offers for selection and has the same honest limit: it is a re-run, not a
retrieval.

**The residual risk, named.** A MODIFIER-impact variant can still matter — through
splicing, through a regulatory element, through a mechanism the annotator does not
model. Decision 6 closes the ClinVar-shaped part of that. What remains is a
variant that is MODIFIER on every transcript, carries no coding or splice term, no
SpliceAI signal above 0.2, and no ClinVar record at all. Such a variant is dropped
here and would have been dropped by `selection` anyway under
`dropped_not_coding_or_splice` — so this decision does not create that exposure,
it inherits it from ADR 0019. What it adds is that the variant is now named in a
file instead of counted in a total.

The bar must therefore stay generous, and the discards must stay inspectable.
Both are enforced: `narrow_to_selection` bounds the first, and
`test_every_dropped_variant_is_named_in_the_ledger` bounds the second by diffing
the annotated artifact across two runs rather than by trusting a counter.

**Determinism (GP-30) is proven, not asserted.**
`tests/integration/test_staged_annotation_determinism.py` runs the synthetic case
through both paths and compares `submission/track1_submission.csv` and
`candidates/ranked_pairs.json` **byte for byte**, under `selection.enabled` both
false and true. A third test asserts the bar is not vacuous, so the first two
cannot pass for the wrong reason. On the scale phantom the selected set is
identical: 3,616 with the bar on and off.

**Measured, on the phantom at 200,000 records, bar off → on:**

| | off | on | |
|---|---:|---:|---|
| records reaching the frequency lookup | 200,000 | 9,078 | **22.0x** |
| `annotated.json` bytes | 176,263,861 | 10,608,383 | **16.6x** |
| evidence items | 304,578 | 19,715 | **15.4x** |
| wall seconds | 27.2 | 5.2 | **5.2x** |
| selected variants | 3,616 | 3,616 | **identical** |

The wall-clock ratio is a **lower bound and is labelled as one**. The phantom's
frequency and clinical adapters are hash functions with no index to open and no
bytes to decompress, so 5.2x is the cost of evidence construction and JSON
serialisation alone. The part the phantom cannot see is the part the real run died
of: a gnomAD region query merges the spans of the ids it is given, so asking about
fourteen thousand coordinates instead of five million collapses the decompressed
byte count by orders of magnitude rather than by the record ratio. Claiming 5.2x
as the whole speedup would be claiming a measurement that was not made, and
claiming 350x would be claiming one that was not made either (GP-20).

**Two smaller changes ride along, both answer-preserving.**
`variants/normalised.json` is streamed rather than built as one list of dicts —
byte-identical output, one fewer full materialisation of a 2.6 GB artifact. And
`_vcf_unescape` in the SnpEff parser gained a `if "%" not in value` guard: it ran
eight unconditional `str.replace` passes over five fields per ANN entry, on the
order of 4x10^9 scans across a whole-genome callset, essentially all finding
nothing.

## Alternatives rejected

**Filter before annotation.** Measured and rejected: the ad-hoc script's own
prefilter removed almost nothing, because the callset was already non-ref. The
selectivity that matters is only computable *after* consequence annotation, which
is why the partition is where it is.

**Emit one `EvidenceItem` per dropped variant.** ADR 0019 rejected this and was
right: ~65 GB of items nobody will cite, of exactly the kind the ledger spill
exists to survive. The discard row is the middle term it could not afford —
per-variant, but not per-variant prose.

**Keep per-reason counts and nothing else.** That is the status quo, and it is
what makes GP-19 only partly satisfied today. A count cannot answer "why isn't
variant X in the output?", which is the question a reviewer actually asks, and
`docs/repo-followup-backlog.md` §1.3 asked for the answer independently of any
performance concern.

**Defer the ClinVar lookup behind the bar as well.** Implemented, measured, and
reverted within the same change: 3,616 → 3,600. The saving was small and the cost
was the rescue that `retain_pathogenic_clinical` exists to provide.

**Switch to `bcftools csq`.** Not rejected, deferred, and recorded as TD-27 with
the evidence: it annotates a genome-wide callset in **~10–20 minutes**, and it is
fast for an architectural reason — it streams the VCF against a GFF3 and writes a
`BCSQ` tag in place, so no per-variant annotated intermediate is ever
materialised. It is also haplotype-aware, which is directly relevant to a
compound-heterozygote question. Two caveats stop it being a drop-in: it has no
MANE/canonical transcript selection, so the obvious "take the first
protein-coding consequence" loop mis-annotates whichever candidate happens to be
listed first, and does it silently; and `--local-csq`, which is easy to pass by
reflex, turns off the haplotype awareness the tool is being credited
for. Our SnpEff path is already at 20 JVM launches, so it is not the emergency the
19.5-hour number made it look like.

**Retain everything and buy a bigger disk.** The volume is 60 GB and the
projection was 95 GB, so this is not available today. It is also the wrong shape
of answer: the ad-hoc script did the same work with the same tool in ~14 minutes
against this run's 19.5-hour projection, and more disk does not fix that.
