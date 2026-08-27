# Track 1 fast path — from VCF arrival to a first ranking

**Purpose.** The proband WGS VCF (~300 MB, GRCh38, already called) is not on this
machine yet. This document is the runbook for the first hour after it lands: what to
look at, in what order, with what thresholds, using only what is already on disk.

**Required sequence, in priority order.**
`VCF inspection → known MVA genes → rare damaging same-gene pairs → genome-wide compound-het.`
Every stage produces a submittable ranking on its own. Do not wait for gnomAD or VEP
before producing one.

Companion artifacts: `knowledge/disease/mva_panel.tsv` (116 genes, cited),
`knowledge/disease/mva_panel_sources.md`, `docs/references/track1-submission-contract.md`.

---

## 0. Preconditions and the two rules that bound everything below

- **The VCF never enters the repository.** It lives under `$MVA_WORKSPACE`, outside
  the working tree (ADR 0006, GP-40). Every path below is relative to that workspace.
- **These commands print counts and vocabulary, never records.** No `CHROM/POS/REF/ALT`,
  no genotype string, no sample name is ever pasted into the repo, a commit message,
  a log, or an agent transcript (GP-41, GP-42, GP-44). Where a step must show a sample
  ID, it is for the human operator's terminal only. Snippets below are written to obey
  this: they emit `Counter` summaries and set sizes.
- **`bcftools` and `tabix` are NOT installed on this machine** (`which bcftools` → not
  found). Everything is `cyvcf2` 0.34.0 / `pysam` 0.24.0, both installed and verified.
  `pysam` bundles htslib, so **`pysam.tabix_index(path, preset="vcf", force=True)` can
  build a missing index in-process** — you do not need the CLI tool to get random access.

Run every snippet as `uv run python - <<'PY' … PY` from the repo root.

---

## 1. VCF inspection — the exact checklist, in order

Set `VCF` once:

```python
VCF_PATH = os.environ["MVA_WORKSPACE"] + "/inputs/<proband>.vcf.gz"
```

### Step 1 — Header: build, contig naming, sample count

This is the first command on the real file. It answers the three questions that
change everything downstream.

```python
from cyvcf2 import VCF
import re, collections

v = VCF(VCF_PATH)

hdr = v.raw_header
print("n_samples          :", len(v.samples))
print(
    "sample_ids_hash    :", [hash(s) & 0xFFFF for s in v.samples]
)  # never print names into a transcript
print("n_contigs_declared :", len(v.seqnames))
print("contig_style       :", "chr-prefixed" if v.seqnames[0].startswith("chr") else "BARE")
print("first_5_contigs    :", v.seqnames[:5])
print(
    "reference_lines    :",
    [l for l in hdr.splitlines() if l.startswith(("##reference", "##assembly"))],
)
print(
    "source_lines       :",
    [
        l[:160]
        for l in hdr.splitlines()
        if l.startswith(
            (
                "##source",
                "##GATKCommandLine",
                "##DRAGEN",
                "##bcftools",
                "##deepvariant",
                "##fileformat",
                "##FILTER=<ID=Lo",
                "##pipeline",
            )
        )
    ],
)
```

> **`contig_style` is the single most dangerous value in this project.** The scorer
> compares `chrom` with `.strip()` only — no prefix or case normalisation — and the
> answer key is `chr`-prefixed (`chr15`, `chr1`, `chr7`; see the submission contract).
> A bare-contig VCF scores **zero while looking correct to a human reader**.
> The repo already defends this: `mva.models.genome.normalise_contig` forces UCSC
> style at parse time and a test asserts the renderer emits it. Record the answer and
> move on — do not "fix" the VCF.
>
> **Confirmed hazard, not theoretical:** the ClinVar VCF already on disk uses **bare**
> contigs (`v.seqnames[:5] == ['1','2','3','4','5']`). Any ClinVar join therefore
> crosses a naming boundary in *at least one* direction whatever the proband VCF uses.
> Normalise both sides through `normalise_contig` before comparing.

### Step 2 — Is it single-sample or a trio?

`len(v.samples)` from step 1. This is the highest-value single fact after contig style.

| `n_samples` | Consequence |
|---|---|
| 1 | Phase is UNKNOWN for every pair (TD-07). §5 applies in full. |
| 3 | **Segregation is available.** TD-07 is lifted; a trans-confirmed compound het becomes provable and the inheritance component starts discriminating instead of sitting flat at 0.55. Escalate immediately — this changes the ranking design, not just its inputs. |

### Step 3 — What is already in INFO? (this decides whether the gnomAD download matters)

```python
recs = [h.info() for h in v.header_iter()]
info_ids = sorted(h["ID"] for h in recs if h.get("HeaderType") == "INFO")
format_ids = sorted(h["ID"] for h in recs if h.get("HeaderType") == "FORMAT")
filter_ids = sorted(h["ID"] for h in recs if h.get("HeaderType") == "FILTER")
print("INFO  :", info_ids)
print("FORMAT:", format_ids)
print("FILTER:", filter_ids)

FREQ = {
    "AF",
    "AF_popmax",
    "AF_grpmax",
    "gnomAD_AF",
    "gnomADg_AF",
    "gnomADe_AF",
    "MAX_AF",
    "AC",
    "AN",
    "POPMAX",
}
ANNO = {"CSQ", "ANN", "BCSQ", "EFF", "LOF", "NMD"}
print("carries_population_frequency:", sorted(FREQ & set(info_ids)))
print("carries_consequence_annot   :", sorted(ANNO & set(info_ids)))
```

> **Why this matters more than it looks.** A `CSQ`/`ANN` field means the VCF is already
> VEP- or snpEff-annotated, which supplies both the gene assignment and the consequence
> severity — the two things this pipeline is otherwise missing entirely (§4). A
> `gnomAD_AF`/`AF_grpmax` field makes the multi-hundred-GB gnomAD sites download
> **redundant**, and that download is currently 5/24 chromosomes complete with 1/24
> indexed (§2). Check this before spending bandwidth.
>
> Beware the ambiguous case: a bare `AF` emitted by GATK is the **sample's** alt allele
> frequency, not a population frequency. Read its `Description` before trusting it:
> `[h for h in recs if h.get("ID") == "AF"]`.

### Step 4 — Does FORMAT carry phase? (this decides whether §5 is a hard problem)

From step 3's `format_ids`. Look for `PS`, `PGT`, `PID`.

| Present | Consequence |
|---|---|
| `PS` (and/or `PGT`/`PID`) | GATK-style physical phasing ran. `mva.prioritization.pairing.infer_phase` **already resolves these** — same phase set, opposite haplotype slots ⇒ `TRANS_CONFIRMED`. Any candidate pair within an active region is decidable *today*, with no new code. |
| absent | Every pair is `PhaseStatus.UNKNOWN`. §5 applies in full. |

Physical phasing only bridges variants the caller saw in one active region — tens to a
few hundred bp — so it will resolve a minority of pairs. That minority is worth a great
deal: a `TRANS_CONFIRMED` pair takes `phase_weights.trans_confirmed = 1.00` instead of
`unknown = 0.55` and separates from the field immediately.

### Step 5 — FILTER vocabulary and how much of the file is PASS

```python
from collections import Counter

filt, n = Counter(), 0
for rec in VCF(VCF_PATH):
    n += 1
    filt[rec.FILTER or "PASS"] += 1  # cyvcf2 returns None for PASS
print("total_records:", n)
for k, c in filt.most_common():
    print(f"  {k:30s} {c:>10,}  {100 * c / n:5.2f}%")
```

A single-sample GRCh38 WGS is normally 4.0–5.5 M records. Far fewer means an exome or a
pre-filtered file; far more means multi-allelics are unsplit or the file is unfiltered.
Note the vocabulary itself — `LowQual`, `RefCall`, VQSR tranches, DRAGEN's `lod_fstar`
— and feed it to `mva.ingestion.qc`, which maps a non-PASS FILTER to `low_quality_call`
**without deleting the record** (GP-13).

### Step 6 — Multi-allelic split, normalisation, and record shape

```python
from collections import Counter

shape = Counter()
n = 0
for rec in VCF(VCF_PATH):
    n += 1
    nalt = len(rec.ALT)
    shape["multiallelic" if nalt > 1 else "biallelic"] += 1
    if nalt == 1:
        r, a = rec.REF, rec.ALT[0]
        if a in ("*", "."):
            shape["star_or_missing_alt"] += 1
        elif len(r) == len(a) == 1:
            shape["snv"] += 1
        elif len(r) != len(a):
            shape["indel"] += 1
            # left-aligned + parsimonious ⇒ exactly one shared leading base, none trailing
            if not (r[0] == a[0] and len(r) > 1 and len(a) > 1 and r[-1] == a[-1]):
                pass
            else:
                shape["indel_not_parsimonious"] += 1
        else:
            shape["mnv"] += 1
    if n >= 500_000:
        break  # a 500k-record sample is plenty to characterise shape
print(shape)
```

- `multiallelic > 0` ⇒ the file is **not** split. The repo's reader decomposes them, and
  `VariantRecord.allele_fraction` (not `Genotype.allele_balance`) is the field that stays
  correct after the split — see TD-16 and `filters._quality_flags`. Nothing to do beyond
  confirming the reader is used.
- `indel_not_parsimonious > 0` ⇒ not `norm`-ed. Matters for the **ClinVar join and the
  submission key only**: the scorer compares `(chrom, pos, ref, alt)` literally, so the
  submitted representation must be the one in the proband VCF, **not** a re-normalised
  one. Never renormalise a variant you intend to submit.
- `star_or_missing_alt > 0` ⇒ spanning deletions present; `filters.apply_hard_filters`
  already drops these as un-analysable.

### Step 7 — Genotype and depth distribution (sanity + the mosaicism question)

```python
import numpy as np
from collections import Counter

zyg, dp, ab = Counter(), [], []
for i, rec in enumerate(VCF(VCF_PATH)):
    t = rec.gt_types[0]  # 0 HOM_REF 1 HET 2 UNKNOWN 3 HOM_ALT
    zyg[{0: "hom_ref", 1: "het", 2: "unknown", 3: "hom_alt"}[int(t)]] += 1
    if t in (1, 3):
        d = rec.gt_depths[0]
        if d and d > 0:
            dp.append(int(d))
            if t == 1:
                ad = rec.format("AD")
                if ad is not None and ad.shape[1] >= 2:
                    tot = int(ad[0][:2].sum())
                    if tot > 0:
                        ab.append(float(ad[0][1]) / tot)
    if i >= 500_000:
        break
print(zyg, "het/hom_alt =", round(zyg["het"] / max(zyg["hom_alt"], 1), 2))
print("depth  median/p10/p90:", np.percentile(dp, [50, 10, 90]).round(1) if dp else "n/a")
print(
    "het AB median:",
    round(float(np.median(ab)), 3) if ab else "n/a",
    " fraction 0.10<=AB<0.25:",
    round(float(np.mean([(0.10 <= x < 0.25) for x in ab])), 4) if ab else "n/a",
)
```

Expected het/hom-alt ≈ 1.5–2.0 for WGS. Median het AB should sit near 0.50; median
depth near 30× for a clinical WGS.

> **Read the low-AB fraction carefully in this case.** `config/default.yaml` flags
> `0.10 ≤ AB < 0.25` as `possible_mosaic` rather than as noise, deliberately, because
> this proband has a **mosaic** aneuploidy phenotype. An elevated low-AB fraction here
> may be biology, not artifact. It is also why the allele-balance phase proxy in §5 is
> specifically untrustworthy for *this* patient.

### Step 8 — Index the file if it is not already indexed

Region queries make the tier-1/tier-2 gene sweeps instant instead of a full pass.

```python
import os, pysam

if not (os.path.exists(VCF_PATH + ".tbi") or os.path.exists(VCF_PATH + ".csi")):
    pysam.tabix_index(VCF_PATH, preset="vcf", force=True)  # bundled htslib; no tabix binary needed
```

Requires the file to be bgzip-compressed. If it is plain `.vcf` or plain `.gz`,
recompress with `pysam.tabix_compress` first. **Write the index next to the VCF inside
`$MVA_WORKSPACE`, never into the repo.**

### Step 9 — The one-pass census that sizes every later stage

Run this once and keep the numbers; §5's combinatorics are read off them rather than
guessed.

```python
from collections import Counter
import pickle
from mva.models.genome import normalise_contig

per_contig = Counter()
het_sites = Counter()
n_pass = n_het = 0
for rec in VCF(VCF_PATH):
    if rec.FILTER is not None:  # non-PASS
        continue
    n_pass += 1
    c = normalise_contig(rec.CHROM) if not rec.CHROM.startswith(("chrUn", "HLA")) else None
    if c is None:
        continue
    per_contig[c] += 1
    if rec.gt_types[0] == 1:
        n_het += 1
        het_sites[c] += 1
print("PASS records:", n_pass, " PASS hets:", n_het)
print(dict(sorted(per_contig.items())))
```

### Step 10 — Write the inspection record

Produce `$MVA_WORKSPACE/runs/<run>/qc/vcf_inspection.json` with **only** the derived
facts: `n_samples`, `contig_style`, `declared_build`, `caller`, `info_ids`,
`format_ids`, `filter_vocabulary`, `n_records`, `pass_fraction`, `multiallelic_split`,
`het_hom_ratio`, `median_depth`, `carries_population_frequency`,
`carries_consequence_annotation`, `has_phase_set`. No coordinates, no sample names.
This file is safe to reason over and to hand to another agent.

### Decision table — what each answer changes

| Finding | Immediate consequence |
|---|---|
| contig style is bare | Renderer must add `chr`. Already enforced; verify the test still passes. |
| `n_samples == 3` | TD-07 lifted. Redesign inheritance scoring before ranking. Highest-impact possible finding. |
| `CSQ`/`ANN` present | Gene assignment and consequence severity are solved. Skip §4's MANE-interval fallback entirely. |
| `gnomAD_AF`/`AF_grpmax` present | **Cancel the gnomAD sites download.** Use the in-file AF. |
| `PS`/`PGT` present | Close-together pairs become phase-decidable with existing code. |
| multi-allelics unsplit | Use the repo reader, and read `allele_fraction`, not `allele_balance` (TD-16). |
| record count ≪ 4 M | It is an exome, or pre-filtered. Non-coding hypotheses are off the table; say so in the report. |

---

## 2. What is on disk, and what actually works today

| Resource | Path under `mva-resources/` | State | Usable now? |
|---|---|---|---|
| ClinVar VCF + index | `clinvar/clinvar.vcf.gz{,.tbi}` | complete, `fileDate=2026-08-22`, GRCh38, 4 467 990 records | **yes** — region queries verified. **Bare contigs.** Carries `CLNSIG`, `CLNREVSTAT`, `GENEINFO`, and `MC` (SO-term molecular consequence) |
| gnomAD gene constraint | `gnomad/gnomad.v4.1.constraint_metrics.tsv` | complete, 95 MB, 211 523 rows | **yes** — pLI, LOEUF (`lof.oe_ci.upper`), `mis.z_score` per transcript |
| gnomAD exome sites VCFs | `gnomad/v4.1_exomes/` | **5 of 24 chromosomes present, only chr21 complete and indexed** (download log records one `OK`) | **no**, except chr21. And they are *exomes*: no non-coding AF for a WGS proband |
| MANE summary | `mane/MANE.GRCh38.v1.5.summary.txt.gz` | complete, 19 363 MANE Select genes | **yes** — HGNC/Ensembl/NCBI ids, GRCh38 contig |
| MANE genomic GTF | `mane/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz` | complete, 19 363 `gene` features, **`chr`-prefixed** | **yes** — this is the gene-assignment fallback (§4) |
| HPO | `hpo/{hp.obo,phenotype.hpoa,genes_to_phenotype.txt,phenotype_to_genes.txt}` | complete, v2026-06-23 | **yes** |
| DDG2P / ClinGen | `genepanels/` | complete (`DDG2P.csv`; the `.gz` beside it is an HTML error page) | **yes** — already distilled into `knowledge/disease/mva_panel.tsv` |
| snpEff | `snpeff/snpEff/snpEff.jar`, `snpeff/jdk` (openjdk 21 symlink) | **core unpacked and a JDK is linked**; the `GRCh38.115` database was still downloading in 8 parts at 2026-08-28 04:21 | **check, do not assume** — see the probe below |
| VEP | — | not present, and its cache is ~26 GB | **no** |

**This table is a snapshot taken 2026-08-28 and it moves under you** — downloads were
in flight while it was written and other agents are wiring adapters
(`src/mva/annotation/{clinvar_vcf,gnomad_sites,snpeff_local}.py` now exist). **Probe the
state, never recall it:**

```python
import pathlib, subprocess

RES = pathlib.Path("../mva-resources")
jar = RES / "snpeff/snpEff/snpEff.jar"
jdk = RES / "snpeff/jdk/bin/java"
db = list((RES / "snpeff").glob("snpEff_v5_4_GRCh38*.zip"))
part = list((RES / "snpeff").glob("*.part*"))
print("snpeff_jar      :", jar.exists())
print("jdk             :", jdk.exists())
print("db_zip_complete :", bool(db), " still_downloading:", bool(part))
print("gnomad_indexed  :", sorted(p.name for p in (RES / "gnomad/v4.1_exomes").glob("*.tbi")))
if jdk.exists() and jar.exists():
    print(
        subprocess.run(
            [str(jdk), "-jar", str(jar), "databases"], capture_output=True, text=True, timeout=120
        ).stdout[:200]
    )
```

**Net position at the time of writing: no consequence annotator was runnable end to
end** — the snpEff core and a JDK were present but the GRCh38 database was incomplete,
and only 1 of 24 gnomAD sites files was indexed. Plan around that (§4). If the probe
above says otherwise when the VCF lands, so much the better: use the real annotator and
skip the fallback.

---

## 3. Tiered candidate generation

Thresholds and where each comes from. Values named `config.frequency.*` /
`config.quality.*` are in `config/default.yaml` and are documented heuristics, not
calibrated cut-points (GP-32); changing one requires an ADR.

### Tier 1 — biallelic hits in the established MVA genes

**Gene set:** `panel_tier == 1` in `knowledge/disease/mva_panel.tsv` — 7 genes:
`BUB1B`, `CEP57`, `TRIP13`, `CENATAC`, `MAD1L1`, `BUB1`, `BUB3`.
All are autosomal recessive or curated only as recessive, which matches a
compound-heterozygous answer.

**Selection:** every record whose position falls in one of the 7 MANE gene spans
(padded ±5 kb to catch promoter/UTR/splice-region calls) **and** carries an alt allele.
Then keep genes with either one HOM_ALT call or ≥ 2 HET calls.

**Thresholds — deliberately almost none.** Seven genes span roughly 400 kb; a WGS yields
on the order of a few hundred variants there. That is small enough to look at *all* of
it, so at tier 1:

| Filter | Tier 1 setting | Source |
|---|---|---|
| FILTER | non-PASS **flagged, not dropped** | GP-13, `filters.apply_hard_filters` |
| Population AF | none applied. If AF is available, record it; do not filter on it | the causal allele may be a founder allele, and gnomAD coverage is absent for 23/24 chromosomes anyway (§2) |
| Consequence | none applied | there is no annotator (§2); filtering on an absent predicate would delete the answer |
| Depth / GQ | `config.quality.min_depth = 10`, `min_genotype_quality = 20` — **flag only** | `config/default.yaml` |
| Allele balance | het band 0.25–0.75; 0.10–0.25 ⇒ `possible_mosaic`, not noise | `config/default.yaml`, chosen for this disease |

**Cost:** one indexed region query per gene, or a single pass. Seconds.

**Output:** if any tier-1 gene shows a plausible biallelic configuration, that is the
submission's rank-1 row, and everything below it is insurance.

### Tier 2 — the broader CIN / repair / centrosome / microcephaly panel

**Gene set:** `panel_tier == 2` — 109 genes. Filter to the 92 rows with
`inheritance ∈ {autosomal_recessive, x_linked}` for the compound-het hypothesis; the 12
`autosomal_dominant` rows (`NIPBL`, `SMC3`, `RAD21`, `WT1`, `KIF11`, `TUBB`, `GMNN`, …)
stay in scope for a single-hit dominant/de-novo hypothesis and must not be paired.

**Thresholds:**

| Filter | Tier 2 setting | Source |
|---|---|---|
| Population AF (max across populations, `min_allele_number = 2000`) | ≤ `config.frequency.max_plausible_recessive = 0.01`; **absent AF is not a filter** and scores `absent_frequency_score = 0.5` | `config/default.yaml`, ADR 0010, GP-14 |
| ClinVar | `CLNSIG ∈ {Pathogenic, Likely_pathogenic, Pathogenic/Likely_pathogenic}` is a **promotion**, never a requirement | `clinvar.vcf.gz`, GP-12 |
| Consequence | if `CSQ`/`ANN` present: `IMPACT ∈ {HIGH, MODERATE}` or splice-region. Otherwise: in a MANE CDS or within 8 bp of an exon boundary (§4) | VEP severity ladder / MANE GTF |
| Quality | as tier 1, flag-only | `config/default.yaml` |
| Gene prior | order genes by `panel_tier`, then `evidence_tier` (definitive > strong > moderate > limited > candidate) | `mva_panel.tsv` |

**Hard rule:** a hit in one of the 11 `evidence_tier == candidate` genes (`TTK`,
`AURKB`, `MAD2L1`, `ESPL1`, `ZW10`, `ZWILCH`, `KNTC1`, `NDC80`, `NUF2`, `SPC24`,
`SPC25`) is reportable as a **research observation only** — no curated gene–disease
validity exists for any of them (see `mva_panel_sources.md` §5.3). It may occupy a
low-EPCR row; it may not be the primary finding.

**Cost:** 109 region queries, or one pass. Under a minute either way.

### Tier 3 — genome-wide compound het

**Gene set:** all 19 363 MANE Select genes.

| Filter | Tier 3 setting | Source |
|---|---|---|
| FILTER | PASS only (this is the one tier where the volume justifies it) | pragmatic; record the count dropped |
| Population AF | ≤ `config.frequency.rare = 0.0001` where AF is available | `config/default.yaml` |
| Consequence | HIGH impact, or ClinVar P/LP, or MANE CDS + nonsynonymous | §4 |
| Gene constraint | prefer `lof.oe_ci.upper` (LOEUF) < 0.6 **or** `lof.pLI` > 0.9 for a dominant hypothesis; for recessive, constraint is *weak* evidence — recessive genes are frequently unconstrained. Use as a tie-break only | `gnomad.v4.1.constraint_metrics.tsv` |
| Zygosity | ≥ 2 HET, or 1 HOM_ALT | — |

**Cost:** one full pass plus an interval join. Minutes.

**Honest expectation:** tier 3 is currently the weakest stage, because two of its three
discriminators (population AF, consequence severity) are unavailable for most of the
genome (§2). If the VCF carries `CSQ` and `gnomAD_AF` (step 3), tier 3 becomes strong
immediately. If it does not, tier 3 is a completeness check, not a ranking.

---

## 4. Needed vs nice-to-have, per tier

| Capability | Tier 1 | Tier 2 | Tier 3 | Have it now? |
|---|---|---|---|---|
| The VCF | required | required | required | pending |
| `mva_panel.tsv` | **required** | **required** | not used | **yes** |
| MANE gene spans (gene assignment) | **required** | **required** | **required** | **yes** — `MANE…ensembl_genomic.gtf.gz`, chr-prefixed |
| ClinVar P/LP lookup | nice | **strongly wanted** | **strongly wanted** | **yes**, indexed |
| HPO phenotype match | nice | wanted | wanted | **yes** |
| gnomAD gene constraint | not used | tie-break | tie-break | **yes** |
| Population allele frequency | **not needed** (7 genes) | wanted | **required for ranking** | **no** — unless the VCF carries it (step 3) |
| Consequence severity (VEP/snpEff) | not needed | **strongly wanted** | **required for ranking** | **no** — unless the VCF carries `CSQ`/`ANN` |
| Trio genotypes | would settle it | would settle it | would settle it | unknown until step 2 |

### The two gaps, and the fallback for each

**Gap 1 — gene assignment.** `pairing.generate_pairs` is gene-scoped: *"a variant with
no gene annotation forms no candidate."* With no annotator, `VariantRecord.gene_symbols`
is empty for every record and **the entire prioritisation stage silently returns zero
candidates.** This is the first thing to fix and it needs no external download.

*Fallback (implement this):* parse `gene` features out of
`MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz` (19 363 rows, already `chr`-prefixed, drop
`*_fix`/`*_alt` contigs), build a per-contig sorted interval list, and assign
`gene_symbols` by `bisect` on POS with ±5 kb padding. ~40 lines, no dependency, runs in
under a second. It gives gene assignment and nothing else — no consequence, no
transcript position. That is sufficient for tiers 1–3 to *generate* candidates.

**Gap 2 — consequence severity.** Ranked fallbacks, best first:
0. A working local snpEff (probe in §2). `mva.annotation.snpeff_local` already wraps
   it, offline-pinned. If the GRCh38 database finished downloading, this gap closes
   entirely and tier 3 becomes a real ranking stage.
1. `CSQ`/`ANN` already in the VCF (step 3) — use it, it is the real thing.
2. ClinVar `MC` (SO terms) and `CLNSIG` for variants already in ClinVar — real evidence,
   but only covers known variants.
3. MANE CDS intervals from the same GTF: in-CDS vs intronic vs intergenic, plus
   "within 8 bp of an exon boundary" as a splice-region proxy. **This is a location
   class, not a consequence prediction** — it cannot tell missense from nonsense. Label
   it as such wherever it is used (GP-12: prediction is not observation).
4. Nothing. Then rank on gene prior, rarity and ClinVar alone, and say so in the report.

---

## 5. Compound-het logic without phase

### 5.1 What already exists — do not rebuild it

`src/mva/prioritization/pairing.py` already implements, and is tested:

- `generate_pairs(variants, *, max_pairs_per_gene=20)` — groups alt-carrying variants by
  gene, emits **all** unordered pairs `C(n,2)` per gene plus single-variant candidates
  for hom-alt/hemizygous calls, mtDNA calls, `possible_mosaic` calls and HIGH-impact
  lone hets. Deterministic ordering, truncation flagged with `gene_pair_cap_truncated`.
- `infer_phase(a, b)` — resolves `CIS_CONFIRMED` / `TRANS_CONFIRMED` **from caller phase
  sets only** (`PS` + phased GT + opposite haplotype slots), correctly handling
  multiallelic split records via `alt_allele_index`. Everything else is `UNKNOWN` with a
  recorded reason and inter-site distance. Nothing upgrades UNKNOWN to trans (GP-15).
- `PRODUCED_INHERITANCE_MODELS` / `UNPRODUCED_INHERITANCE_MODELS` — an explicit,
  test-asserted statement of which inheritance models a proband-only VCF can support.
- Flag promotion onto the candidate row (`common_variant`, `possible_mosaic`, …).

**The pairing and phase logic is done.** Do not write a second one.

### 5.2 The delta — four changes, in priority order

**D1. Populate `gene_symbols`.** Gap 1 above. Without it `generate_pairs` returns an
empty tuple and every tier produces nothing. *This is the single highest-value code
change in the project right now.* Nothing else on this list matters until it is done.

**D2. Make `max_pairs_per_gene` truncation damage-aware.** Today, when a gene exceeds
the cap, candidates are sorted by `sort_key()` — **genomic position** — and the first 20
are kept. A gene with 30 rare hets keeps the 20 leftmost, which is arbitrary with
respect to how damaging they are: the true pair can be silently dropped, flagged only by
`gene_pair_cap_truncated`. *Fix:* rank member variants by (impact severity, rarity,
ClinVar significance) **before** `combinations`, so truncation removes the least
plausible pairs rather than the rightmost ones. Keep the flag. This is a change inside
`generate_pairs`; the honest alternative is to raise the cap for panel genes only.

**D3. Add a `PhaseStatus.TRANS_LIKELY` producer, or explicitly decline to.**
`config/default.yaml` already defines `phase_weights.trans_likely = 0.85`, but §5.1's
`infer_phase` never emits it — it is currently unreachable configuration. Either
implement a legitimate producer (see 5.4) or delete the weight so the config does not
imply a capability that does not exist.

**D4. Add a gene-prior input to scoring.** `mva_panel.tsv`'s `panel_tier` and
`evidence_tier` are the strongest discriminator available in the no-phase, no-AF,
no-consequence situation, and nothing currently reads them (`paths.gene_panel` is
declared in `src/mva/config.py:299` but **no code reads the file** — grep finds no
consumer). Wire the panel into `mechanistic_relevance` (weight 0.08) or, better, add a
`gene_prior` component with its own weight and an ADR.

### 5.3 How many pairs this generates, and how to bound it

The count is `Σ_genes C(n_g, 2)`, so it is entirely governed by the filter applied
*before* pairing. Order-of-magnitude, for a single-sample GRCh38 WGS — **treat these as
estimates to be replaced by the step-9 census, not as facts**:

| Pre-pairing filter | het variants entering | pairs generated | Feasible? |
|---|---|---|---|
| none (all PASS hets, whole gene span ±5 kb) | ~2.0 M, ~100+ per gene | ~10⁸ | no — this is why the cap exists |
| PASS + coding/splice only | ~20 k, ~1–3 per gene | ~10⁴ | yes |
| PASS + coding/splice + AF ≤ 0.01 | ~2–5 k | ~10³ | yes, comfortably |
| PASS + coding/splice + AF ≤ 1e-4 (tier 3 setting) | ~150–500 | **~10–100** | trivially |
| tier 1 (7 genes, no filter) | tens | **< 100** | trivially |

Three bounds, applied in this order:

1. **Filter before pairing, not after.** `combinations` on an unfiltered gene is the
   only way to reach 10⁸. `select_candidate_variants` already exists for exactly this
   and applies the `max_plausible_recessive` ceiling without deleting anything the
   caller holds.
2. **Keep the per-gene cap**, but make it damage-aware (D2), and raise it for the 116
   panel genes where the search space is small and completeness matters more.
3. **Cap the reported set, not the computed set.** The submission takes at most 10 rows
   (hard limit in the scorer). Compute thousands, rank them, emit ≤ 10. The evidence
   store keeps the rest for the report.

### 5.4 Ranking without phase — what actually discriminates

With `n_samples == 1` and no `PS`, every two-het candidate receives
`phase_weights.unknown = 0.55`. **A constant discriminates nothing.** The phase component
contributes exactly zero ordering information across the compound-het candidate set, so
the ranking has to be carried elsewhere. In priority order:

1. **Gene prior** (`panel_tier`, then `evidence_tier`). A biallelic hit in *BUB1B*
   outranks one in an unpanelled gene by an enormous margin, and this is the only signal
   that is fully available today. See D4.
2. **Phenotype match.** `phenotype.hpoa` terms for the proband's HPO profile against the
   gene's disease. `HP:0002667` (nephroblastoma) is the sharpest single discriminator in
   this differential — annotated 6/6 to MVA3/*TRIP13* and to MVA1/*BUB1B*, and **not at
   all** to MVA2/*CEP57*.
3. **ClinVar significance** of either member allele. Real, curated, and on disk.
4. **Rarity**, where AF exists.
5. **Consequence severity**, where it exists. A HIGH+HIGH pair (two LoF alleles) outranks
   HIGH+MODERATE, which outranks MODERATE+MODERATE.
6. **Inter-site distance**, as a *reporting* aid only. `PhaseEvidence.distance_bp` is
   already recorded. Two sites ≤ 500 bp apart (`READ_BACKED_PHASING_SPAN_BP`) are
   resolvable by targeted re-analysis of the alignments — that is the
   `recommended_next_test`, not a score.

**Two proxies that must NOT be used to upgrade phase:**

- *Allele-balance symmetry.* Two true trans hets in a diploid region both sit near 0.5,
  so a 0.50/0.95 pair looks wrong. In **this** proband it is unusable: mosaic aneuploidy
  produces genuinely skewed allele fractions, which is precisely why
  `mosaic_allele_balance_floor` exists. Using AB as a phase proxy here would
  preferentially down-rank the mosaic signal the case is about.
- *"Two rare hets in a recessive disease gene are probably in trans."* True as a
  population statement, false as evidence about this patient. GP-15 forbids it, and the
  scoring already grants UNKNOWN a generous 0.55 for exactly this reason.

**What *would* legitimately produce `TRANS_LIKELY` (D3):** parental genotypes (step 2),
caller phase sets (step 4), or read-backed phasing from the alignments. The BAM/FASTQ
are in the gated dataset; if the team has access, a read-backed phase call on the top
five pairs is the single highest-yield follow-up. Nothing derivable from a proband-only
VCF alone justifies it.

### 5.5 Mapping a ranked pair to the submission

Per `docs/references/track1-submission-contract.md`:

- **One compound-het pair is ONE row**, using `chrom_2/pos_2/ref_2/alt_2`. Not two rows.
- `proband_id` must be exactly `PROBAND01`.
- `chrom` must be `chr`-prefixed and is compared **raw** — this is where step 1 pays off.
- `ref`/`alt` are upper-cased and compared literally: submit the representation that is
  in the proband VCF (see step 6 — do not renormalise).
- `epcr ∈ (0, 1]`, float. Rank is derived from `epcr` descending.
- Maximum 10 rows. Rank tiers are `1 → 100`, `2–3 → 50`, `4–5 → 25`, `6–10 → 10`.
- Partial credit exists for a compound-het row on set intersection, so a row with one
  correct member still scores half — **an uncertain second allele is worth submitting,
  not worth omitting.**

---

## 6. The first thirty minutes, in order

| # | Action | Blocking? |
|---|---|---|
| 1 | Steps 1–4 of §1 (header, samples, INFO, FORMAT) | yes — everything branches on these |
| 2 | Steps 5–7 (FILTER, shape, genotype census) | no, run while step 3 proceeds |
| 3 | Implement D1: MANE-interval gene assignment | **yes — nothing produces candidates without it** |
| 4 | Tier 1 sweep over the 7 MVA genes | — |
| 5 | Tier 2 sweep over the 109 panel genes + ClinVar join | — |
| 6 | Render a submission CSV from whatever tiers 1–2 produced | do this before starting tier 3 |
| 7 | Tier 3 genome-wide | — |
| 8 | Re-rank as gnomAD / a consequence annotator land | — |

Step 6 is not optional. Six submission attempts are allowed and the best counts; a
tier-1/tier-2 submission banked early costs nothing and protects against everything
after it going wrong.

---

## 7. The single biggest scientific risk in this fast path

**Almost all of the discriminating power sits in the 116-gene panel, for as long as the
two signals that would work genome-wide — population allele frequency and predicted
consequence — are unavailable.** At the time of writing, gnomAD sites data was 1 of 24
chromosomes indexed and exome-only, and the snpEff GRCh38 database was still
downloading. Both are being fixed by other work in flight; **run the §2 probe before
accepting this risk as live.** While it is live:

- Tiers 1 and 2 are strong: small gene sets, real curated priors, real ClinVar evidence,
  a real phenotype match.
- Tier 3 is, today, close to unranked. Without AF and without consequence severity, a
  genome-wide list of "genes with ≥ 2 rare-ish hets" has no defensible ordering beyond
  gene constraint, which is weak evidence for a recessive disorder.

**Therefore: if the causal gene is not on the panel, this fast path will very likely not
find it.** The panel is built to make that unlikely — 116 genes covering the MVA
phenotype series, the SAC and kinetochore machinery, centrosome and centriole biology,
cohesin and condensin, the Fanconi and DNA-damage-response differentials, primary
microcephaly, replication licensing and primordial dwarfism, and Wilms-tumour
predisposition — but "unlikely" is not "no", and MVA has an OMIM phenotype-series gap
at MVA5/MVA6 that we could not resolve from any file on disk. A novel gene is a live
possibility in a challenge built around a hard case.

**Mitigations, in order of value:**

1. **Check step 3 first.** If the VCF already carries `CSQ`/`ANN` and `gnomAD_AF`, this
   entire risk evaporates and tier 3 becomes a real ranking stage. Do this before
   anything else.
2. **Finish the snpEff GRCh38 database download and verify the JAR runs** — separating
   HIGH from MODIFIER is most of the value, and both the adapter
   (`mva.annotation.snpeff_local`) and a JDK are already in place.
3. Finish the gnomAD download (or index what is already there), or substitute a lighter
   source of AF. `mva.annotation.gnomad_sites` is already written against it.
4. Keep tier 3 in the run and in the report even when it is unranked, so that a
   negative genome-wide result is recorded as *"searched, could not rank"* rather than
   as *"nothing there"* (GP-14).

Second-order risks, all already tracked: contig-prefix mismatch (guarded by
`normalise_contig` plus a test, but re-verify against the real header), structural
variants and repeat expansions being invisible (TD-13, TD-03 — a compound het of one SNV
and one deletion cannot be assembled from this VCF at all), and mtDNA heteroplasmy being
unmeasured (TD-14).
