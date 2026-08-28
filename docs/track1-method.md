# Track 1 method — what was actually run

**No patient record appears in this file** — no coordinate, no genotype, no
candidate gene symbol, no sample identifier. The finding itself lives on the
encrypted case volume, outside this repository (ADR 0006, GP-40).

That is a narrower claim than "no patient data", and the difference is
deliberate. Aggregate measurements derived from a patient's genome are still
derived from that patient: a per-chromosome dosage table discloses genetic sex,
and an exact record count fingerprints the callset. Both appeared in an earlier
draft of this file under a blanket "no patient data" banner, which was wrong.
Counts below are therefore reported at the precision the method needs and no
finer, and results that would disclose a genetic characteristic are described by
their bearing on the method rather than by their value.

Executed 2026-08-28 against the challenge proband WGS VCF (GRCh38, single sample,
GATK-called, ~5.0 million records).

---

## The route

```
VCF inspection -> curated MVA panel -> rare damaging same-gene pairs -> genome-wide compound-het
```

Every stage produces a submittable ranking on its own; the later stages were used
to *check* the earlier ones, not to unblock them.

### 0. Inspection, without reading records

`tools/inspect/vcf_schema.py` emits header declarations, contig style, sample
**count** and salted digests, INFO/FORMAT/FILTER **ids only**, and aggregate scan
counts. It never prints a record. Two facts from it shaped everything downstream:

* **Contigs are bare** (`15`, not `chr15`), Ensembl/NCBI style. The scorer compares
  chromosome strings raw, so emitting a correct answer with a bare contig scores
  exactly zero. Conversion happens once, in the renderer, and a test asserts it.
* **A large minority** of records are indel-bearing, and multi-allelic records are **not**
  split. Both make allele representation a first-order correctness problem rather
  than a formatting detail.

### 1. Panel derivation, from files rather than from memory

The MVA gene set was derived by querying the HPO release on disk for the genes
annotated to the MVA disease series (OMIM:257300, 614114, 617598, 620153, 620189
and ORPHA:1052). This reproduced, independently, the tier-1 set of the
repository's own 116-gene curated panel (`knowledge/disease/mva_panel.tsv`, built
from DDG2P, ClinGen, HPO and ClinVar). Gene coordinates came from MANE v1.5, read
per gene — not typed from memory. That precaution is not theoretical: an earlier
query in this project used a GRCh37 coordinate against a GRCh38 resource and
returned **zero records**, silently, which is indistinguishable from a real
negative.

### 2. Allele representation, once, for every join

Every allele — from the proband, from gnomAD, and from ClinVar — is reduced to
its minimal left-aligned representation against the GRCh38 FASTA before any
comparison, through the single rule in `src/mva/alleles.py` (ADR 0018).

This is the difference between a working join and a plausible-looking failure.
Measured earlier on the real gnomAD chr21 shard, over a 520 kb exonic window:
1,382 indel records, of which **1,029 (74.5%) sit in repeat tracts**. Joining
their right-shifted spellings without a reference matched **0 of 1,029**; with a
reference, **989**. Thirty of the recovered records are variants gnomAD calls
*common*, 22 of them above 5% allele frequency. Without canonicalisation those
thirty would have been scored as novel and ultra-rare — the exact profile of a
causal variant.

### 3. Annotation

| Resource | Release | Used for |
|---|---|---|
| SnpEff | 5.4c, database GRCh38.115 | consequence and impact |
| gnomAD | v4.1 **exomes**, 24 chromosome shards | population frequency |
| ClinVar | VCF 2026-08-22 | clinical assertions |
| MANE | GRCh38 v1.5 | gene intervals, transcript identity |
| HPO | 2026-06-23 (`hp.obo`, `phenotype.hpoa`, `genes_to_phenotype.txt`) | phenotype semantics |
| Reference | GRCh38 no-alt analysis set | left-alignment |

Consequence was derived twice: once restricted to the canonical transcript, and
once across **all** annotated transcripts. Collapsing to MANE Select is a
data-loss bug — a variant can be benign on the canonical transcript and
splice-disrupting on a tissue-relevant isoform — so the all-transcript pass is
what the conclusion rests on.

### 4. Selection

> **What this section documents.** The analysis that produced the submitted
> ranking was executed as a scripted pipeline over the resources listed above.
> The library's own `iter_selected` stage (ADR 0019) is **not identical to it**:
> ADR 0019 retains `impact is None` by default and treats biallelic pairing as a
> separate downstream stage, whereas the run described here applied an
> impact filter and a same-gene biallelic rule in one step. Wiring the two into
> agreement is tracked as the composition-root item in `docs/next-actions.md`.
> Until that lands, this document describes **what was run**, not what
> `mva run all` currently does — and saying so is the point, because a method
> section that describes code the run did not use is worse than no method
> section.
>
> **Updated 2026-08-28 — what closed, and what has not.** The composition root is
> now wired (ADR 0027). `mva run all` on a case with `synthetic: false` binds the
> resources in the table above — SnpEff 5.4c/GRCh38.115 pinned from the reviewed
> installer manifest, gnomAD v4.1 exomes, ClinVar 2026-08-22, MANE v1.5 — and
> passes the GRCh38 no-alt FASTA to `normalise_variants` **and** to both joining
> adapters, so the left-alignment the scripted run relied on is now what the
> library does too. A case whose resources are absent refuses to start rather than
> falling back to the synthetic tables.
>
> Two divergences remain, and both are stated rather than closed:
>
> 1. **Selection.** `iter_selected` still applies ADR 0019's rules — `impact is
>    None` retained, biallelic pairing as a separate downstream stage — where the
>    scripted run applied an impact filter and a same-gene biallelic rule in one
>    step. The library is the *more* conservative of the two (it deletes less), so
>    the difference can add candidates below the answer and cannot remove the
>    answer, but the two are not step-for-step identical.
> 2. **Nothing here has been executed end to end by the library against the real
>    callset.** The proband VCF lives on the encrypted case volume, outside what
>    the agents building this pipeline may read (ADR 0008, GP-40). The wiring is
>    verified against the real public releases and the synthetic fixtures; the
>    real run is a human step in `docs/submission-runbook.md`.
>
> So a reader reproducing this should expect the same resources, the same
> representation rule and the same ranking inputs, and should verify the ranking
> itself against `submission/track1_submission.csv` rather than assuming it.


Retained: PASS-filtered, HIGH or MODERATE impact, rare, and biallelic within a
gene (>=2 heterozygous or >=1 homozygous damaging allele).

**A variant with no frequency record is retained, not dropped.** Absence of a
gnomAD row means the site was never assessed; it does not mean the allele is
absent from the population, and it certainly does not mean AF = 0. Treating
"unknown" as "common" would discard the novel private variant that a recessive
diagnosis most often turns on; treating it as "rare" would flood the ranking. It
is carried forward with its ignorance recorded (GP-14, ASSUMPTION-FREQUENCY-01).

Impact `None` means NOT ASSESSED and is never read as MODIFIER (ADR 0016).
MODIFIER is a positive prediction of negligible effect; `None` is the absence of
a prediction. Conflating them silently converts ignorance into evidence of
benignity.

### 5. Ranking

Candidate genes were ranked by HPO semantic similarity against the proband's
phenotype profile: Lin (1998) pairwise similarity, information content computed
from the `phenotype.hpoa` annotation corpus (267,062 rows) under the true-path
rule — **never from graph depth** — aggregated as best-match average in the
proband -> gene direction (ASSUMPTION-PHENOTYPE-04, -05).

The reverse direction is deliberately not used for the score. Averaging over
every term of a *gene's* curated profile penalises a well-studied gene for
features nobody in this clinic assessed, which is the same category error as
treating a missing annotation as a negative finding.

### 6. Stated exclusions

* **The extended MHC** (GRCh38 chr6:28,510,120-33,480,577) is excluded before
  ranking. It is the most polymorphic region of the genome and a well-known
  source of short-read alignment artefacts, and its genes carry hundreds of
  curated immune phenotype annotations — so annotation density inflates their
  similarity score for reasons unrelated to the patient. Applied as a region rule
  to every gene, never as a judgement about a particular candidate.
* **Structural and copy-number variation was not assessed.** The input is an
  SNV/indel callset; an exon-level deletion is invisible to it. This is a real
  gap in coverage, not a claim that none exists.
* **gnomAD genomes was not used** (524 GB, beyond this machine's storage). Deep
  intronic variants are therefore not frequency-assessable, and candidates that
  depend on intronic rarity are carried explicitly rather than dismissed.

---

## Aneuploidy was tested for directly

A report naming an aneuploidy syndrome should say whether aneuploidy is visible
in the data. Two orthogonal tests were run over the whole callset:

* **Allele balance** at every heterozygous biallelic PASS SNV in a bounded
  depth window. Median allele balance sits at the balanced expectation on all
  autosomes.
* **Relative depth** per chromosome against the autosomal median, flat across
  all autosomes. The sex chromosomes act as the positive control that makes the
  autosomal null meaningful: the method demonstrably resolves a real
  copy-number difference, so its silence on the autosomes is evidence rather
  than insensitivity. Both the counts and the per-chromosome values are
  genetic characteristics of the child and stay on the case volume; only the
  method and its verdict are reported here.

**Neither test detects whole-chromosome imbalance, and that is the expected
result.** MVA aneuploidy is *variegated*: different cells carry different wrong
chromosomes, so gains and losses cancel in bulk DNA, and the mosaic fraction in
peripheral blood is often low. The diagnostic assay is a karyotype scored cell by
cell, not short-read WGS. Reporting this as a negative, rather than omitting it,
bounds what this pipeline can and cannot see.

---

## Submission shaping

The scorer's mechanics were verified by re-implementing and **executing**
`evaluation.py` from the challenge Space (see
`docs/references/track1-submission-contract.md`). Four rules follow, and all four
are enforced by the renderer and asserted before any file is written:

1. **Fill all ten rows.** F-max is a maximum over thresholds, so a row below the
   answer's EPCR cannot lower either metric. Unused rows forfeit upside.
2. **Never emit two equal EPCRs.** A tie pulls a wrong row into the same
   prediction set: F-max falls 1.0 -> 0.667, and if file order puts the wrong row
   first it also takes rank 1, costing 50 further points.
3. **Never split a pair across two rows.** A compound-het proposal is one row
   using the `_2` columns; splitting it caps the result at partial credit,
   100 -> 50 rank points for nothing.
4. **Every chromosome `chr`-prefixed.** The single highest-consequence formatting
   detail in the project.

The rendered submission is re-parsed and scored against a replica of the real
scorer before it is accepted, under several hypotheses about the answer key —
including hypotheses in which our own answer is wrong — to confirm the row
ordering degrades gracefully rather than cliff-edging to zero.

`notes` is capped at 120 characters against a restricted character set. The
submission is published CC-BY 4.0 under the challenge terms, so it carries
structured mechanism-class values only. Free clinical narrative must not reach a
public artifact, and the renderer refuses to write one that does.
