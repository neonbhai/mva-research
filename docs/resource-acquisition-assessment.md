# Resource acquisition assessment — allele-frequency reference for a WGS proband

**Date:** 2026-08-28 · **Status:** recommendation, not yet executed
**Scope:** whether to continue, narrow or cancel the running gnomAD v4.1 exomes pull, and what to
substitute. All sizes below are **measured** by HTTP `HEAD`/range request on the dates shown, not
quoted from documentation.

> **Boundary observed.** Nothing in `SageBio/mva-hackathon-2026-data` was accessed — not the VCF,
> not the FASTQs, not the dataset README, with or without a token. That dataset is patient data
> under WCG IRB protocol #REDACTED-PROTOCOL. Every claim about the proband VCF below is inference from the
> **public** challenge Space source, and is labelled as such. The task brief also mentioned reading
> the dataset card "if publicly rendered"; that conflicts with the explicit prohibition on the
> README, so the prohibition was followed and the card was not fetched.

---

## Verdict

> **NARROW.** Let the gnomAD v4.1 **exomes** pull run to completion — it is **184.8 GB, not 250 GB**,
> and leaves **~360 GB free**, so it never threatens the 200–250 GB reserve. Narrow it by
> (a) **re-ordering** the queue so chr15 (BUB1B) / chr11 (CEP57) / chr5 (TRIP13) land next and
> fetching every `.tbi`; (b) **not acquiring gnomAD genomes (524.4 GB) or joint (~800 GB)** — they do
> not fit on this machine at any reserve, so genome-wide gnomAD is off the table regardless of what
> we decide; (c) **not acquiring dbSNP** — its `FREQ` field has no allele number and is disqualified
> by ADR-0010; (d) adding **1000G 30x (30.2 GB)** as the genome-wide non-coding common-variant
> filter. A slim-and-delete step can cut the 184.8 GB to ~5.4 GB and is documented below, but it
> collides with the existing `resources.yaml` / `gnomad_sites.py` contract — **hold it in reserve**
> rather than running it, since the disk headroom makes it unnecessary.

The three numbers that carry the verdict:

| | measured |
|---|---|
| gnomAD v4.1 exomes, all 24 shards | **184.8 GB** (the 250 GB premise was 65 GB high) |
| gnomAD v4.1 genomes, all 24 shards | **524.4 GB** — against 524 GB free, i.e. genome-wide gnomAD is impossible here at *any* reserve |
| Same exome content as a slim per-ancestry AF table | **5.4 GB** (31.9 B/variant gzip, measured) — a **34×** reduction |

And the finding that changes the shape of the problem:

> **Remote `tabix` range queries against the gnomAD bucket work today, with the pysam already in
> `.venv`.** BUB1B (chr15:40,161,000–40,320,000) returned **50,474 gnomAD exome variants with full
> per-ancestry AC/AN in 28.6 s**, storing zero bytes. Three MVA genes on three chromosomes:
> **64,756 variants in 75.7 s**. This removes the "second download round mid-analysis" risk that
> made deferral unattractive.

---

## Measured state of the running job

Sampled twice, 90 s and 100 s windows, while 6 `curl` streams were live:

| quantity | measured |
|---|---|
| free space | 524 GiB (was 531 GiB at brief time) |
| on disk in `v4.1_exomes/` | 13.6 GB |
| aggregate throughput | **10.49 MB/s**, then **9.76 MB/s** → call it 10 MB/s |
| remaining | 171.2 GB |
| ETA to completion | **~4.8–5.0 h** |
| **free space on completion** | **~360 GB** |

Two corrections to the brief:

1. **The total is 184.8 GB, not ~250 GB.** So `524 − 171 = 353` GB free, not `532 − 250 = 282`. The
   200–250 GB reserve survives with ~110 GB of slack. **Disk was never the binding constraint.**
2. **chr17–22 are not "already down".** All six are still in flight — they are the six running
   `curl`s. Their true combined size is **40.1 GB**, of which 13.6 GB has landed. The "~6 GB"
   figure was a mid-transfer snapshot.

The binding constraint is **time and bandwidth**, not gigabytes: for ~4.8 h the link is saturated,
starving every other acquisition and every remote probe on the critical path.

---

## Q1 — Is the challenge VCF already population-annotated?

**Public evidence says: almost certainly not, and it must be verified from the header the moment the
file is legitimately in hand.**

What the public Space source actually supports
(`https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/raw/main/`):

- `tabs/about.py`: *"You'll be given the child's genomic data (FASTQ + VCF) and a description of
  their symptoms."* — describes the deliverable as raw reads plus a called VCF. **No annotation is
  mentioned anywhere in the Space.**
- `tabs/rules.py`: *"File types: genomic data in VCF format; raw sequencing data optionally in
  BAM/CRAM. Phenotypic data as standardized HPO terms. Data size: ~85 GB (single-subject dataset)."*
  — "VCF format" as a bare file-type declaration, with no annotation claim.
- `tabs/faq.py`, `config.py`, `app.py`, `utils.py`, `tabs/submit_track1.py`: no mention of gnomAD,
  VEP, dbSNP, annotation, DRAGEN or GATK. The Space is a registration/leaderboard front end and says
  nothing about the upstream pipeline.

Inference from the naming convention `WGS_SPECIMEN_FLOWCELL.vcf.gz` — `FLOWCELL` is a NovaSeq
(S4-family) flowcell ID — is that this is standard secondary-analysis output named by flowcell,
i.e. DRAGEN or GATK germline. Both emit **population-annotation-free** VCFs by default. A 300 MB
bgzipped WGS VCF is also consistent with a plain single-sample SNV/indel callset (~4.5–5M records);
a gnomAD-annotated WGS VCF would typically be considerably larger.

**The trap to check for explicitly.** DRAGEN and GATK both write an `INFO/AF`, and it is the
**sample-level allele fraction**, not a population frequency. Reading it as population AF would
mark every het variant as "AF ≈ 0.5 → common" and destroy the analysis. The header
`##INFO=<ID=AF,...>` description distinguishes them unambiguously.

**Required verification, first 60 seconds of legitimate access** (reads only the header):

```bash
# from the repo root, using the venv that already has pysam
./.venv/bin/python -c "
import pysam,sys
h=pysam.VariantFile(sys.argv[1]).header
print('--- source/reference ---')
for r in h.records:
    if r.key in ('source','reference','DRAGENCommandLine','GATKCommandLine','fileDate'): print(r)
print('--- INFO keys ---'); print(sorted(h.info.keys()))
print('--- AF description (is it population or sample?) ---')
if 'AF' in h.info: print(h.info['AF'].description)
for k in ('gnomad_AF','gnomAD_AF','AF_popmax','AF_grpmax','CSQ','ANN'):
    if k in h.info: print('POPULATION ANNOTATION PRESENT:', k, '->', h.info[k].description[:200])
" /path/to/WGS_SPECIMEN_FLOWCELL.vcf.gz
```

**If that header shows `CSQ`/`ANN` with gnomAD fields, or `AF_grpmax`, stop and cancel the pull
entirely** — most of this document becomes moot and the remaining need is only the fields the
annotation omitted. Treat this as a live branch, not a formality.

---

## Q2 — Smaller, officially published, locally queryable AF resources

All sizes measured 2026-08-28 by `HEAD`. "Per-anc AC/AN" is the column that matters: ADR-0010's
`min_allele_number` guard needs an **allele number**, not just a frequency.

| resource | measured size | scope | samples | per-anc AC **and** AN? | licence | verdict |
|---|---|---|---|---|---|---|
| **gnomAD v4.1 exomes** (running) | **184.8 GB** (chr1 17.50, chr2 13.30 … chrY 0.11) | coding only (UKB ∩ Broad capture, each padded 150 bp) | 730,947 | **Yes** — `AC_/AN_/AF_` × {afr, amr, asj, eas, fin, mid, nfe, sas, remaining} + `grpmax` + **`faf95` / `fafmax_faf95_max`** | CC0 / no restriction | **Keep, slimmed** |
| gnomAD v4.1 genomes | **524.4 GB** (chr1 41.05, chr2 43.36 …) | **genome-wide** | 76,215 | Yes (+ `ami`) | CC0 | **Impossible** — 524.4 GB against 524 GB free leaves zero reserve |
| gnomAD v4.1 joint | chr1 **67.11**, chr21 **11.03**, chr22 **14.57** → **~800 GB** extrapolated | genome-wide | 807,162 | Yes, incl. joint FAF | CC0 | **Impossible** |
| **dbSNP b157 GRCh38** | **27.52 GB** (`GCF_000001405.40.gz`) + 3.1 MB `.tbi` | genome-wide | varies by study | **NO — AF only, no AC/AN** | US public domain | **Reject for rarity** (see below) |
| **1000G 30x GRCh38** (3,202-sample phased panel) | **30.2 GB** (chr1 2.23 … chrX 2.66) | genome-wide | 3,202 | Yes, but **AN too small** (below) | IGSR — open, no restriction | **Accept for stage A only** |
| VEP cache 114 GRCh38 | **24.91 GB** (merged: 27.27 GB) | genome-wide, but see Q3 | n/a | **NO — AF only** | Apache-2.0 (code); data unrestricted | Consequence yes, frequency no |
| Exomiser 2512 hg38 | **21.94 GB** + 12.23 GB phenotype = **34.17 GB** | genome-wide | n/a | **Yes** — but Java-only access | AGPL-3.0 (software); bundle terms unverified | Only real alternative; heavy integration |

### There is no official gnomAD slim / AF-only release

Enumerated the bucket: `release/4.1/` contains only `constraint/`, `exome_cnv/`, `genome_sv/`,
`ht/`, `local_ancestry/`, `lof_curation/`, `pext/`, `tsv/`, `vcf/`, and `vcf/` contains exactly
`exomes/`, `genomes/`, `joint/`. **No AF-only, sites-slim or minimal sites product exists.** The
5.4 GB slim table below therefore has to be derived locally — which is cheap, and is the recommendation.

### Why dbSNP's `FREQ` is disqualifying, not merely weaker

Decoded the real header and first data records from `GCF_000001405.40.gz`:

```
##INFO=<ID=FREQ,Number=.,Type=String,Description="An ordered list of allele frequencies as
  reported by various genomic studies, starting with the reference allele...">

NC_000001.11 10001 rs1570391677 T A,C . . ...;FREQ=KOREAN:0.9891,0.0109,.|SGDP_PRJ:0,1,.|dbGaP_PopFreq:1,.,0
```

Three fatal properties:

1. **No allele number, anywhere.** `FREQ` is a study→frequency string. ADR-0010's
   `min_allele_number` guard is *unimplementable* against it. Every dbSNP row would enter as
   "AN unreported", which ADR-0010 deliberately keeps **eligible** ("an unreported cohort size is
   unknown, not small") — so dbSNP rows would flow straight past the guard.
2. **It contains exactly the pathology ADR-0010 exists to stop, amplified.** `SGDP_PRJ:0,1` is
   AF = **1.0** from the Simons Genome Diversity Panel — a few hundred individuals, often a handful
   per population. Under ADR-0010's maximum-across-populations rule with no AN to guard on, a
   tiny-cohort AF of 1.0 sets the maximum and flags a genuinely ultra-rare variant
   `common_variant`. This is the AC=1/AN=40 founder-allele failure from the ADR, at larger scale
   and with the guard structurally disabled.
3. Contigs are RefSeq accessions (`NC_000001.11`), needing a rename map before any join.

**Use dbSNP for rsID assignment if wanted. Never as a frequency source in this pipeline.**

### Why 1000G is a stage-A filter and not a rarity scorer — with the guard arithmetic

The panel does carry proper per-superpopulation `AC_*`/`AN_*`. Measured from real records on chr21:

```
chr21:5030578 C>T   global AC=74 AN=6404 AF=0.01156
  AFR AC=74  AN=1786    AMR AC=0 AN=980    EAS AC=0 AN=1170
  EUR AC=0   AN=1266    SAS AC=0 AN=1202
```

**Every superpopulation AN (980 – 1,786) is below `frequency.min_allele_number = 2000`.** Under the
repo's own configured guard, *no 1000G ancestry row is eligible to set the maximum* — only the
global row (AN = 6,404) survives. So in this pipeline 1000G behaves as a **global-AF-only**
resource, and the per-ancestry columns are decorative. That is a direct, non-obvious consequence of
ADR-0010 that should be recorded rather than discovered at runtime.

Its resolution ceiling, quantified:

- Smallest non-zero AF representable: 1/6,404 = **1.56 × 10⁻⁴**.
- Rule of three, zero observations at AN = 6,404: 95% upper bound ≈ **4.7 × 10⁻⁴**.
- So 1000G can assert *"AF < ~5 × 10⁻⁴"* and nothing finer. **It cannot distinguish 10⁻⁴ from 10⁻⁶
  — which is exactly the band a causal ultra-rare allele lives in.**

Conversely it is near-perfect at the job being asked of it in stage A: a 1% allele has expected
count 64 in 6,404 chromosomes, so P(missing it entirely) ≈ e⁻⁶⁴ ≈ 10⁻²⁸. **Commonness is decidable
at 3,202 samples; rarity is not.**

---

## Q3 — Do VEP cache or Exomiser suffice for initial population filtering?

**Neither replaces gnomAD for the frequency contract this pipeline has written. VEP cannot support
ADR-0010 at all; Exomiser can, but only from Java.**

### VEP GRCh38 cache — 24.91 GB measured

Two independent disqualifiers, both from Ensembl's own cache documentation:

1. **Coverage is gated on dbSNP accessioning, not on gnomAD.** Ensembl states plainly: *"gnomAD is
   not a variant accessioning body... any gnomAD variant that are not accessioned will not be
   available in the cache."* Their worked example shows `5-32100960-ATAAG-A` returning **no
   frequency information at all** from the cache. The cache's frequency layer is a *subset* of
   gnomAD — silently, and with no flag distinguishing "absent from gnomAD" (informative: rare) from
   "present in gnomAD but unaccessioned in dbSNP" (uninformative). In a pipeline whose ADR-0010
   treats missing frequency as a **blocking open question**, a silently truncated frequency source
   is worse than none.
2. **Granularity is AF only.** VEP emits `gnomADe_AF`, `gnomADe_AFR_AF`, `gnomADg_NFE_AF` and so on
   — per-population **frequencies with no AC and no AN**. `min_allele_number` cannot be evaluated.

The cache does bundle gnomAD exomes **v4.1** and genomes **v4.1** (release 116 table), which is
current — the problem is the accessioning filter and the missing AC/AN, not staleness.

**Verdict: acquire the VEP cache for consequence/transcript annotation, which is what it is good at.
Do not source frequency from it.**

### Exomiser 2512 — 21.94 GB (hg38) + 12.23 GB (phenotype) = 34.17 GB measured

Better than expected. Its frequency model **does** carry allele counts —
`exomiser-core/.../frequency/Frequency.java`:

```java
public record Frequency(FrequencySource source, float frequency, int ac, int an, int homs) {
```

and `FrequencySource` enumerates `GNOMAD_E_{AFR,AMR,ASJ,EAS,FIN,NFE,MID,OTH,SAS}` and
`GNOMAD_G_{...,AMI,...}` — i.e. **per-ancestry AC/AN from both gnomAD exomes and genomes**, plus
ClinVar and HPO in the same bundle. On the science, this is the one alternative that could satisfy
ADR-0010, and it is genome-wide.

The cost is integration, and on this deadline it is decisive:

- The data ships as H2 / MVStore `.mv.db` files designed for in-process access from the Exomiser
  Java application. There is no supported way to query them from this Python/Snakemake/DuckDB
  pipeline without either embedding Exomiser or reverse-engineering MVStore.
- It duplicates gnomAD rather than replacing it — the same underlying numbers, behind a JVM.
- Bundle redistribution terms (it aggregates ClinVar, gnomAD, HPO, and historically HGMD-PUBLIC)
  were **not verified in this pass** and should be before any reliance. The Exomiser software is
  AGPL-3.0.
- `openjdk@21` was observed installing on this machine, so a JVM is arriving anyway — but "arriving"
  is not "integrated and tested".

**Verdict: not for the primary path now. It is the strongest fallback if the slim-gnomAD route
fails, and it is the right answer for a project with more than a day left.**

---

## Q4 — Deferring shard acquisition until candidate chromosomes are known

**The scheme is sound, and the risk the brief correctly identified has been largely eliminated by
measurement.**

### The two-stage scheme

**Stage A — genome-wide common-variant removal (no gnomAD needed).**
Source: 1000G 30x, 30.2 GB, or its ~2 GB slimmed sites-only derivative. Drop anything with global
AF > 1% (and, for a recessive hypothesis, > 0.5%). Of ~4.5–5M variants in a WGS callset, the large
majority are common polymorphisms present in 1000G; typical rare-disease practice retains
**~150,000–400,000** variants after a 1% population filter alone. Justified by the arithmetic
above: at AN = 6,404 a 1% allele is missed with probability ~10⁻²⁸, so this stage has essentially
no false-negative risk for the class it removes.

**Stage B — consequence and gene filtering (no frequency needed).**
Restrict to protein-altering / splice-region / ClinVar-flagged, and to genes with phenotype or
constraint support. This is where WGS collapses hardest: of ~4.5M variants only **~20,000–25,000**
are coding, and after a rare + protein-altering filter a single exome/genome typically leaves
**~150–500** candidate variants, spread over **~100–300 genes on 15–23 distinct chromosomes** for a
genome-wide scan, or **a handful of chromosomes** once a phenotype-driven gene panel is applied.

**Stage C — exact per-ancestry gnomAD frequency for survivors only.**

### The risk, and why it is now small

The brief's concern was: *deferring means a second download round mid-analysis, and if the causal
variant's chromosome is fetched late we lose time at the worst moment.* That was the right worry
about a **download**-based deferral. It does not apply to a **range-query**-based one.

**Measured, on this machine, against the live gnomAD bucket, with `.venv`'s pysam 0.24.0:**

```
BUB1B    chr15:40,161,000-40,320,000   50,474 variants   28.6 s
CEP57    chr11:95,792,000-95,840,000    9,353 variants   30.2 s
TRIP13   chr5:900,000-940,000           4,929 variants   17.0 s
TOTAL    64,756 variants, 3 chromosomes, 75.7 s, 0 bytes stored
```

Records come back with the full INFO payload — `AF`, `AC_*`/`AN_*` for all nine ancestry groups,
`grpmax`, `AN_grpmax`, `faf95`, `fafmax_faf95_max`. Example returned live:

```
chr15:40161009 A>AGG  AF=4.90e-06  AF_grpmax=1.23e-04  AN_grpmax=16256  grpmax=amr  faf95max=2.13e-05
```

So **stage C for a candidate gene set costs ~30 s per chromosome, right now, with nothing
downloaded**. There is no "second round" and no worst-moment stall. Note this was measured while six
downloads were saturating the link; on a quiet link it is faster.

**One caveat that matters to this repo:** TD-06 records that the pipeline runs under Python-level
network denial. Remote tabix is therefore appropriate for **acquisition, triage and exploration**,
but the **reproducible pipeline run must read from a local artifact**. That is precisely why the
recommendation keeps a local slim table rather than making the pipeline network-dependent —
determinism (and the byte-identity guarantees the golden tests rest on) requires a fixed local file.

---

## Q5 — Disk plan, and whether re-aligning 79 GB of FASTQ is worth it

### Allocation

| bucket | GB | note |
|---|---|---|
| Reserve for FASTQ/BAM/intermediates (user requirement) | **250** | untouched, and now comfortable |
| gnomAD v4.1 exomes **slim AF table** | **~6** | replaces 184.8 GB; measured 5.4 GB gzip |
| Transient: in-flight shards during slimming | **~45 peak** | reclaimed continuously; peak = largest concurrent set |
| 1000G 30x (or ~2 GB slimmed) | **30** | genome-wide stage-A filter |
| VEP GRCh38 cache (consequence only) | **25** | optional but recommended |
| Already held (ClinVar, HPO, constraint, DDG2P, ClinGen, MANE, GTF) | **~0.5** | |
| Headroom | **~165** | |

Peak transient stays under ~50 GB, so free space never drops below ~430 GB and the 300 GB disk
guard never fires. Compare with the naive plan: keeping all 184.8 GB of raw shards would leave
~360 GB — still above the reserve, but spending 180 GB to store 413 INFO fields per record when the
pipeline reads about 22 of them.

### Re-aligning the FASTQs: **no, for Track 1** — and the public evidence is unusually direct

Cost, on this machine:

- 79 GB FASTQ in, ~90–110 GB BAM out (~25–30 GB as CRAM), plus sort/merge intermediates at roughly
  1.5–2× — **peak ~250–260 GB**, i.e. the *entire* reserve.
- 30× WGS alignment (bwa-mem2 / DRAGMAP) + sort + dedup + recalibration on an Apple Silicon laptop:
  **~12–30 h wall clock**, single machine. That exceeds the time remaining and would run at the same
  time as everything else.

Benefit for Track 1: **essentially zero**, and the Space says why. `tabs/faq.py`:

> *"rank points (based on how high the true variant(s) land in your ranked list, **with partial
> credit if you recover only one of two compound-heterozygous variants**)"*

and `tabs/about.py`: *"If you correctly identify one of the two variants..."*. **The answer key is a
compound-heterozygous pair** — two SNV/indel alleles in one gene. That is exactly what the provided
called VCF contains, and it is consistent with MVA's established genetics (biallelic *BUB1B*, and
more rarely *CEP57* / *TRIP13*). Re-deriving a callset the NHS lab already produced and validated
would consume the entire reserve and most of the remaining time to recover variants you were handed.

Note in passing: the repo's own synthetic golden case is `chr15:40200000:C:T` + `chr15:40210500:G:A`
— a chr15 compound het sitting inside the BUB1B locus. The demo was already built against the right
shape of answer.

### The TD-03 exception — real, but it does **not** justify the disk here

TD-03 ("no structural/copy-number variant calling") is a genuine gap, and correctly noted as
uncloseable from an SNV/indel VCF: Manta/DELLY/GRIDSS and CNV callers all need alignments. Three
reasons it still does not buy realignment now:

1. **The scored answer is a compound het of two point variants**, per the FAQ text above. SV calling
   cannot improve the metric being scored.
2. **In MVA, the aneuploidy is the phenotype, not the genotype.** Mosaic variegated aneuploidy is
   *caused* by germline biallelic loss of mitotic-checkpoint function (BUB1B/CEP57/TRIP13); the
   whole-chromosome mosaic gains and losses are the downstream cellular consequence. Chasing
   chromosome-scale events in the proband's own WGS would be characterising the phenotype, and a
   single bulk WGS sample is in any case a poor instrument for low-fraction mosaicism.
3. **Check before you build.** `tabs/rules.py` says *"raw sequencing data optionally in BAM/CRAM."*
   If alignments are in the dataset, SV/CNV calling costs **no alignment at all** — only caller time
   and ~30 GB. Confirm what the dataset actually contains before spending 250 GB and 20 h to
   manufacture a BAM that may already exist.

**Residual risk, stated plainly:** if one of the two compound-het alleles were an exon-level
deletion, it would be invisible to both the provided VCF and this plan. That happens in real
diagnostics. It is inconsistent with the FAQ's framing of two rankable variants, so the risk is
judged low — but it is not zero, and it is the one scenario in which TD-03 would have cost the
answer. Revisit only if the dataset turns out to ship BAM/CRAM (cheap) — not by realigning (expensive).

---

## Exact commands

### 1. Do NOT kill the running pull — reprioritise it

The six live `curl`s are chr17–22 (40.1 GB combined, 13.6 GB landed). Let them finish; killing them
discards partial transfers, and `-C -` resume means restarting costs nothing but wastes what is done.
Queue the high-prior chromosomes next — **chr15 (BUB1B), chr11 (CEP57), chr5 (TRIP13)**:

```bash
DEST=/Users/someshwar-tripathi/Contri/bio-hackathon/mva-resources/gnomad/v4.1_exomes
BASE=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes
LOG=/Users/someshwar-tripathi/Contri/bio-hackathon/mva-resources/logs/gnomad.log

# wait for the current batch to drain, then take the MVA-gene chromosomes first
for c in 15 11 5; do
  ( curl -fsS --retry 8 --retry-delay 5 --retry-all-errors -C - --max-time 0 \
      -o "$DEST/gnomad.exomes.v4.1.sites.chr${c}.vcf.bgz" \
      "$BASE/gnomad.exomes.v4.1.sites.chr${c}.vcf.bgz"
    curl -fsS --retry 8 -C - -o "$DEST/gnomad.exomes.v4.1.sites.chr${c}.vcf.bgz.tbi" \
      "$BASE/gnomad.exomes.v4.1.sites.chr${c}.vcf.bgz.tbi" ) >>"$LOG" 2>&1 &
done
```

Also fetch the `.tbi` for every shard — they are ~14 KB each and make the slim step and any local
range query possible:

```bash
for c in $(seq 1 22) X Y; do
  curl -fsS --retry 8 -C - -o "$DEST/gnomad.exomes.v4.1.sites.chr${c}.vcf.bgz.tbi" \
    "$BASE/gnomad.exomes.v4.1.sites.chr${c}.vcf.bgz.tbi" &
done; wait
```

### 2. Slim each shard as it lands, then delete the `.bgz` — 184.8 GB → ~5.4 GB

Measured on the BUB1B region: **122.0 B/variant raw, 31.9 B/variant gzip**, projecting to
**20.8 GB raw / 5.4 GB gzip** for gnomAD v4.1 exomes. Uses `.venv`'s cyvcf2 — no bcftools needed
(bcftools/tabix/bgzip are **not installed** on this machine; `pysam` 0.24.0 and `cyvcf2` 0.34.0 are).

```bash
cat > /private/tmp/claude-502/-Users-someshwar-tripathi-Contri-bio-hackathon/617483ae-6c2e-4a4e-8cec-5c7909799af0/scratchpad/slim_gnomad.py <<'PY'
import gzip, os, sys
from cyvcf2 import VCF
POPS = ["afr","amr","asj","eas","fin","mid","nfe","sas","remaining"]
COLS = (["chrom","pos","ref","alt","AC","AN","AF","AF_grpmax","AN_grpmax","AC_grpmax",
         "grpmax","faf95","fafmax_faf95_max"]
        + [f"{p}_{f}" for p in POPS for f in ("AC","AN")])
src, dst = sys.argv[1], sys.argv[2]
def g(v, k):
    try:
        x = v.INFO.get(k)
    except KeyError:
        return ""
    if x is None: return ""
    return str(x[0]) if isinstance(x, tuple) else str(x)
n = 0
with gzip.open(dst, "wt", compresslevel=6) as out:
    out.write("\t".join(COLS) + "\n")
    for v in VCF(src):
        alt = v.ALT[0] if v.ALT else "."
        row = [v.CHROM, str(v.POS), v.REF, alt,
               g(v,"AC"), g(v,"AN"), g(v,"AF"), g(v,"AF_grpmax"), g(v,"AN_grpmax"),
               g(v,"AC_grpmax"), g(v,"grpmax"), g(v,"faf95"), g(v,"fafmax_faf95_max")]
        for p in POPS:
            row += [g(v, "AC_" + p), g(v, "AN_" + p)]
        out.write("\t".join(row) + "\n"); n += 1
print(f"{os.path.basename(src)}: {n} variants -> {os.path.getsize(dst)/2**20:.1f} MB", flush=True)
PY

# run per shard, delete the .bgz only after the slim file is written successfully
DEST=/Users/someshwar-tripathi/Contri/bio-hackathon/mva-resources/gnomad/v4.1_exomes
SLIM=/Users/someshwar-tripathi/Contri/bio-hackathon/mva-resources/gnomad/v4.1_exomes_slim
mkdir -p "$SLIM"
cd /Users/someshwar-tripathi/Contri/bio-hackathon/mva-research
for f in "$DEST"/gnomad.exomes.v4.1.sites.chr*.vcf.bgz; do
  c=$(basename "$f" .vcf.bgz); o="$SLIM/${c}.af.tsv.gz"
  [ -s "$o" ] && continue
  # only process a shard that is fully downloaded (no live curl writing to it)
  pgrep -fl "curl.*$(basename "$f")" >/dev/null && { echo "skip in-flight $c"; continue; }
  ./.venv/bin/python /private/tmp/claude-502/-Users-someshwar-tripathi-Contri-bio-hackathon/617483ae-6c2e-4a4e-8cec-5c7909799af0/scratchpad/slim_gnomad.py "$f" "$o.part" \
    && mv "$o.part" "$o" && rm -f "$f" && echo "reclaimed $(basename "$f")"
done
```

Re-run that loop after each batch completes (or wrap it in a `while` poll). **Verify a slim file
before deleting its shard** — the `&&` chain above already does, but check row counts on the first
one before trusting the loop.

### 3. Acquire 1000G 30x for the genome-wide stage-A filter — 30.2 GB

Exact base URL:
`https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/`

```bash
KG=/Users/someshwar-tripathi/Contri/bio-hackathon/mva-resources/1000g_30x
mkdir -p "$KG"
B=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV
for c in $(seq 1 22) X; do
  f=1kGP_high_coverage_Illumina.chr${c}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
  curl -fsS --retry 8 --retry-all-errors -C - -o "$KG/$f" "$B/$f" &
  while [ "$(jobs -r | wc -l)" -ge 4 ]; do wait -n; done
done; wait
```

Note chrX is published as `...chrX...v2.vcf.gz` in some mirrors — if chrX 404s, retry with the `.v2`
name. Then slim it the same way (keep `AC`, `AN`, `AF` global plus the superpopulation columns for
provenance, even though ADR-0010's guard will exclude them) — the slimmed panel is ~2 GB and the
30.2 GB of genotype payload can then be deleted.

### 4. Available immediately, no download — targeted per-ancestry gnomAD frequency

Use during triage, before any shard finishes:

```bash
cd /Users/someshwar-tripathi/Contri/bio-hackathon/mva-research
./.venv/bin/python - <<'PY'
import pysam
URL = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr%s.vcf.bgz"
POPS = ["afr","amr","asj","eas","fin","mid","nfe","sas","remaining"]
for chrom, start, end, name in [("15",40161000,40320000,"BUB1B"),
                                ("11",95792000,95840000,"CEP57"),
                                ("5",900000,940000,"TRIP13")]:
    vf = pysam.VariantFile(URL % chrom)
    for rec in vf.fetch("chr"+chrom, start, end):
        pass  # -> rec.info['AC_nfe'], rec.info['AN_nfe'], rec.info['fafmax_faf95_max'], ...
PY
```

### 5. Do NOT start these

```
# gnomAD v4.1 genomes  -- 524.4 GB measured, does not fit (524 GB free)
# gnomAD v4.1 joint    -- ~800 GB extrapolated, does not fit
# dbSNP GRCh38         -- 27.5 GB, FREQ has no AC/AN; disqualified by ADR-0010
```

---

## The scientific cost of this recommendation, stated honestly

This is the part not to skim. The recommendation **loses real allele-frequency coverage**, and a
wrong "rare" call is how false positives are manufactured.

**1. No gnomAD-grade frequency for non-coding variants — the largest loss.**
The proband is WGS; gnomAD exomes covers only the UKB ∩ Broad capture (each padded 150 bp),
i.e. roughly 1–2% of the genome. Of ~4.5–5M variants in a WGS callset only **~20,000–25,000 are
coding** — so gnomAD exomes supplies exact per-ancestry frequency for **well under 1% of the
proband's variants by count**. For the other ~99%, the best available number after this plan is
1000G's global AF at AN = 6,404. That resolves commonness and nothing finer. A deep-intronic,
promoter, UTR or regulatory variant that is genuinely ultra-rare and genuinely causal will arrive at
the ranking stage with either "common/not-common" from 1000G or, if absent from 1000G, **no
frequency data at all** — which ADR-0010 correctly scores at `absent_frequency_score` and flags
`no_frequency_data` with a blocking open question. That is the honest behaviour, but it is not the
same as knowing the variant is rare, and it means **the pipeline cannot distinguish a novel causal
regulatory variant from the tens of thousands of unremarkable rare non-coding variants every genome
carries.** If the MVA answer is regulatory, this plan does not find it — and neither does the
alternative, because gnomAD genomes does not fit on this disk.

**2. Ancestry resolution collapses outside coding regions.**
Inside the exome capture: nine ancestry groups with AC and AN, `grpmax`, and precomputed `faf95`.
Outside it: one global 1000G number, because every 1000G superpopulation AN (980–1,786) falls below
`min_allele_number = 2000`. A founder allele enriched in one ancestry — **the exact case ADR-0010
was written to protect** — is invisible as such in non-coding space. Concretely, a variant at global
AF 0.004 but AFR AF 0.014 passes a 1% global filter and is scored "rare" when it is common in the
relevant population. That failure mode is live for every non-coding variant under this plan.

**3. Reduced power even where covered, at the extreme tail.**
gnomAD v4.1 exomes has 730,947 samples (AN up to ~1.46M) versus 76,215 genomes — so for coding
variants the exomes are the *more* powerful resource, and choosing exomes over genomes is a gain,
not a loss, in the region it covers. The loss is purely of territory, not of depth.

**4. Variant classes still entirely invisible — unchanged by this plan, but they bound the claim.**
SVs and CNVs (TD-03), repeat expansions (TD-13) and mtDNA heteroplasmy (TD-14) are not called and
not scorable. Declining the realignment does not create these gaps — they exist regardless of disk —
but it forecloses closing TD-03 in this round. **A negative result from this pipeline cannot
distinguish "no structural variant" from "never looked."** Any write-up must say so.

**5. One upgrade this plan buys, which partly offsets the above.**
gnomAD ships `faf95` and `fafmax_faf95_max` — the Poisson 95% CI filtering allele frequency, computed
per ancestry group. That is precisely the instrument TD-15 says should replace the unfitted
2000-allele heuristic, and it arrives free in the slim table (both fields are in the extraction
above). Within coding regions the pipeline can move from a reasoned cut-point to the ACMG BA1/BS1
standard quantity. **TD-15 is closeable this round, for coding variants, at no additional download
cost.**

---

## Integration caveat — this collides with `resources.yaml` and `gnomad_sites.py`

Discovered while checking the working tree; **not resolved here**, and it must be settled before the
slim-and-delete step in command 2 is run.

- `knowledge/manifests/resources.yaml` declares one entry per shard (`gnomad_exomes_chr1`,
  `gnomad_exomes_chr10`, `..._tbi`, ...) each pinning a `path` such as
  `gnomad/v4.1_exomes/gnomad.exomes.v4.1.sites.chr1.vcf.bgz` plus `sha256`, `size_bytes` and a
  `status` of `fetched` / `not_fetched`. **Deleting the `.bgz` after slimming would flip every one
  of those entries back to "no local file at the declared path"**, and any provenance or
  acquisition gate keyed on `status: fetched` would fail.
- `src/mva/annotation/gnomad_sites.py` reads the ancestry groups **out of the VCF header** with
  `re.compile(r"^##INFO=<ID=AF_([a-z]+),")` and expects `AF_<grp>` / `AN_<grp>` INFO keys on the
  record. It consumes the **VCF**, not a TSV. A slim TSV is not a drop-in for it.

Three ways to reconcile, in increasing order of work:

1. **Keep the `.bgz` files and skip the slimming.** Costs 184.8 GB, still leaves ~360 GB, still
   inside the 200–250 GB reserve. This is the zero-risk option and the right one if the deadline is
   tight — the slimming is a disk optimisation, not a scientific one.
2. **Slim to a bgzip-compressed *VCF* rather than a TSV**, keeping the `AF_<grp>`/`AN_<grp>` INFO
   keys and header lines that `gnomad_sites.py` already parses, then re-index with `pysam.tabix_index`.
   The adapter keeps working unchanged; the manifest's `path` stays valid; only `sha256` and
   `size_bytes` need re-recording, and the entry should gain a note that it is a locally derived
   subset rather than the upstream artifact. Larger than 5.4 GB (VCF framing and the header cost
   more than TSV) but still roughly an order of magnitude below the raw shards.
3. Teach the adapter a TSV backend. Most work, least appropriate on this deadline.

**Recommendation: option 1 now, option 2 if disk pressure actually materialises.** The measured
headroom (~360 GB free on completion, reserve intact) means the 184.8 GB is affordable, and
preserving the manifest/provenance contract is worth more than 180 GB of a disk that is not full.
This downgrades command 2 from "do this" to "hold in reserve" — the verdict above is otherwise
unchanged, since its load-bearing claims are that genomes/joint are impossible, that dbSNP and VEP
cannot satisfy ADR-0010, and that realignment is not worth 250 GB.

---

## What would change this verdict

- **The VCF header shows gnomAD annotation already present** (Q1) → cancel the pull outright; keep
  only what the existing annotation omits.
- **The dataset ships BAM/CRAM** → TD-03 becomes cheap (~30 GB + caller time, no alignment) and
  should be done; realignment stays rejected.
- **The candidate list lands in non-coding space** → the honest answer is that this machine cannot
  support that analysis at gnomAD grade, and the finding must be reported with its frequency
  evidence marked absent rather than rare. Escalate to cloud compute rather than pretending
  1000G's AN = 6,404 resolves it.
