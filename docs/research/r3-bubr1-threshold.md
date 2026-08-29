# R3 — Where is the residual-BubR1 threshold?

Scope: public literature only. Every claim carries author, year, journal, PMID. Numbers I could
not verify in an abstract or full text are marked **unresolved**, not estimated.

Date of search: 2026-08-28.

---

## 1. Verdict

No single threshold is defensible, because the published data do not describe a step. They
describe a graded dose–response whose steepness depends on which readout you score: premature
chromatid separation appears at reductions as small as a quarter, whole-cell aneuploidy
accumulates somewhere below 58% residual, and per-division missegregation in transformed cells
only becomes gross below about 10%. The "10–60%" range in our document is not one threshold
measured imprecisely by two studies; it is two different quantities measured precisely by two
studies that never measured the same thing.

The most useful correction is not a better number. It is that **"the threshold" is the wrong
noun** — the document has to say which function, in which cell type, over how many divisions.

---

## 2. Every quantitative datapoint found

Residual level is "% of wild-type / control" unless stated. "Unresolved" means the paper is
paywalled, or the number lives only inside a figure image I could not read.

### 2.1 Human cultured-cell titrations

| Source | System | Residual | How quantified | Readout | Result |
|---|---|---|---|---|---|
| Bohers 2008, *Hum Genet* 124:473-8, PMID 18932004 | HeLa, stable shRNA-transduced lines | 77% | **unresolved** (full text paywalled) | PCS, metaphase spreads | PCS induced |
| " | " | 58% | unresolved | PCS; aneuploidy | PCS induced; aneuploidy not "clearly detected" |
| " | " | 14% | unresolved | PCS; aneuploidy | PCS induced; aneuploidy clearly detected |
| " | " | 10% | unresolved | PCS; aneuploidy | as above |
| " | " | 8.5% | unresolved | PCS; aneuploidy | as above |
| Suijkerbuijk 2010, *Cancer Res* 70:4891-900, PMID 20516114 (PMC2887387) | HeLa "HTR34", doxycycline-inducible shRNA + H2B-EYFP | ~13% (1 ng/ml dox) | quantitative immunoblot, BUBR1/tubulin, relative to untreated, mean of 3 experiments | % cells with chromosome missegregation, live imaging | "virtually no impact on the fidelity of chromosome segregation" |
| " | " | ~6% (1000 ng/ml dox) | same | same | missegregation "in the majority of cells" |
| " | HeLa/U2OS, transient shRNA, no rescue | near-zero | — | live imaging; nocodazole arrest | ">80% of the cells showed massive chromosome missegregations"; no mitotic delay |
| Ochiai 2014, *PNAS* 111:1461-6, PMID 24344301 (PMC3910577) | HCT116, biallelic single-base-pair edit 44 kb upstream of *BUB1B* (TALEN) | reduced; **densitometry values printed only inside Fig. 3A — unresolved** | western, β-tubulin normalised | PCS frequency; metaphase chromosome counts; colcemid mitotic index | increased PCS, mosaic aneuploidy, checkpoint defect |

### 2.2 Human patient and carrier cells

| Source | System | Residual | How quantified | Readout | Result |
|---|---|---|---|---|---|
| Suijkerbuijk 2010, PMID 20516114 | MVA patient lymphoblastoid + SV40-fibroblast lines (absent/I909T, 386X/R727C, 753X/R814H, 731X/Y155C) | full-length BubR1 "decreased 2-6 fold" vs control ⇒ **≈17–50% residual** | quantitative immunoblot, "band intensity of BUBR1/Tubulin relative to healthy control cells", mean of 4 experiments | nocodazole and STLC checkpoint (MPM-2 by flow); chromosome alignment under MG132 | checkpoint deficient; alignment defective |
| " | carrier parents (one mutant allele) | "consistently in between those of patients and unrelated healthy controls" — no number given | same | same | not reported as defective |
| " | patient fibroblasts | kinetochore-bound BubR1 quantified separately (BUBR1/ACA relative to control) | immunofluorescence, 25 kinetochores × 5 cells × 3 experiments | kinetochore loading | low, but "residual protein could be recruited to kinetochores, indicative of at least partial functionality" |
| Matsuura 2006, *Am J Med Genet A* 140:358-67, PMID 16411201 | seven Japanese PCS-syndrome families, monoallelic *BUB1B* | "modest decrease" of transcript and protein — **not quantified as a %** | RT-PCR + western | colcemid checkpoint; Cdc20/p55cdc kinetochore localisation | reduced BubR1, diminished kinetochore localisation, checkpoint failure; restored by chr-15 transfer |
| García-Castillo 2008, *Am J Med Genet A* 146A:1687-95, PMID 18548531 | human heterozygous *BUB1B* carrier parents | ~50% by genotype; **not measured** | — | PCS on metaphase spreads | "PCS was observed in 8/8 propositi and **7/7 carrier parents**"; carriers clinically unaffected |
| Ochiai 2014, PMID 24344301 | patient KH01, homozygous regulatory variant | "significant reduction" — not quantified | western | PCS; mosaic aneuploidy | PCS in 67% of lymphocytes; mosaic aneuploidy in 17% |

### 2.3 Mouse allelic series

| Source | System | Residual | How quantified | Readout | Result |
|---|---|---|---|---|---|
| Baker 2004, *Nat Genet* 36:744-9, PMID 15208629 | BubR1^H/H MEFs | ~11% (as cited by Simmons 2019); ~10% (as cited by Wijshake 2012) | **primary paper paywalled — unresolved**; both figures are secondary citations | aneuploidy; senescence | "Graded reduction of BubR1 expression in mouse embryonic fibroblasts causes increased aneuploidy and senescence" (abstract) |
| Wijshake 2012, *PLoS Genet* 8:e1003138, PMID 23300461 (PMC3531486) | BubR1^−/H mice | "about 4% of normal BubR1 protein levels" | secondary citation within paper | viability | die at birth, respiratory failure |
| " | BubR1^H/H mice | "around 10% of normal BubR1 protein levels" | secondary citation within paper | aneuploidy; PCS; progeria | systemic near-diploid aneuploidy, PCS, tumour susceptibility |
| " | BubR1^+/GTTA MEFs (heterozygous MVA knock-in, 2211insGTTA) | "reduced amounts of wildtype BubR1"; **no % given**, blot only | western (Fig. 1E) | aneuploidy, chromosome counts, passage 5 | WT MEFs ~12% aneuploid; +/GTTA "significantly increased by 6%" (i.e. ≈18%) |
| " | BubR1^+/− MEFs | "residual levels of wildtype BubR1 protein are similar in BubR1^+/GTTA and BubR1^+/− MEFs" | western | aneuploidy; PCS | aneuploidy "modestly increased" — same as +/GTTA; **PCS *not* elevated** |
| " | BubR1^+/GTTA MEFs | same level as +/− | western | PCS | PCS "significantly higher... than in wildtype MEFs" |
| Simmons 2019, *Hum Mol Genet* 28:1822-36, PMID 30668728 (PMC6522074) | mouse E15.5 **cortex**, BubR1^H/H | **"about 1.8% of that of WT"** | western of cortical lysate, GAPDH-normalised, N=3 per genotype | cortical size; progenitor division; microcephaly | normal — no microcephaly |
| " | mouse E15.5 cortex, BubR1 CKO (Emx1-Cre on hypomorph background) | **"about 0.5% of that of WT"** | same blot, same normalisation | microcephaly; metaphase fraction; apoptosis | severe microcephaly; reduced metaphase fraction; massive apoptosis |
| Sterling 2023, *Front Cell Dev Biol* 11, PMC10602889 | same BubR1 cKO cortex | near-null (as above) | — | anaphase defects, en face apical explant | segregation defects in "an average of 23% of total mitotic cells... compared to 0% of WT cells"; γH2AX in 10% of prophase cells |
| Yang 2017, *Aging Cell* 16:598-601, PMID 28383136 (PMC5418205) | BubR1 low-expression progeroid mice, adult hippocampus | "low amounts of BubR1" — not quantified in abstract | — | adult hippocampal neurogenesis | deficits in progenitor proliferation and maturation |

### 2.4 The one datapoint that constrains the shape of the curve

Wijshake 2012 reports two genotypes with **the same residual wild-type BubR1** that differ in
PCS: `+/GTTA` has elevated PCS, `+/−` does not, while both have the same modestly-raised
aneuploidy rate. So PCS is not a function of BubR1 level alone. That single comparison is enough
to break "the threshold" as a pure-dosage concept.

---

## 3. Is the Bohers/Suijkerbuijk disagreement real?

**Mostly an artefact of different readouts, plus one low-powered negative and one number that was
never measured. But it hides a sharper real disagreement, which is inside Suijkerbuijk's own
paper.**

### 3.1 The two studies share no measurement

| | Bohers 2008 | Suijkerbuijk 2010 |
|---|---|---|
| Depletion | stable shRNA transduction, selected clones | doxycycline-inducible shRNA in one line (HTR34) |
| Exposure | chronic, many divisions | acute; duration is in the supplement — **unresolved** |
| Readout 1 | PCS on metaphase spreads | % of cells with missegregation, live imaging |
| Readout 2 | aneuploidy (accumulated) | — |
| Confound | clonal selection: five separate clones, no isogenic dose control | none of that; but only ~350 cells across 3 experiments and 3 conditions |

Aneuploidy is the time-integral of missegregation. Missegregation is the per-division rate. A rate
too small to see in roughly 100–120 imaged divisions still accumulates into visible aneuploidy in a
clonal line carried through weeks of passage. The 4× "disagreement" is close to what you would
predict from comparing an integral to a derivative.

### 3.2 Neither paper measured 50%

Bohers' points are 8.5, 10, 14, 58, 77. "Aneuploidy was clearly detected in cells with residual
BUBR1 below 50%" (abstract) is an interpolation across a gap they did not sample. The defensible
restatement of their own data is narrower and less useful:

> aneuploidy present at ≤14% residual, not clearly detected at ≥58%.

That brackets the transition inside `(14%, 58%)` and says nothing about where in it.

### 3.3 Suijkerbuijk's negative at 13% is weak

Their Fig. 4D lower panel reports "average of 3 experiments, a total of at least 350 cells" across
three dox conditions. Roughly 100–120 divisions per condition. If baseline HeLa missegregation runs
a few percent, a rise to 6–8% would not reliably clear significance. "Virtually no impact" is
compatible with a real but modest elevation.

It is also the only readout they titrated. **They never titrated the checkpoint.** The nocodazole
and alignment assays in that paper are mutant-replacement experiments, not dose–response
experiments. So Suijkerbuijk 2010 contains no measurement of how much BubR1 the spindle assembly
checkpoint needs.

### 3.4 The real contradiction is elsewhere

Suijkerbuijk 2010 reports that MVA patient cells carry full-length BubR1 "decreased 2-6 fold"
against healthy controls — roughly 17–50% residual — and that those cells are frankly checkpoint
deficient and misalign chromosomes. Their own HeLa at ~13% is not. Patient cells at a *higher*
residual level are more defective than HeLa at a lower one.

Three explanations, none excluded by the data:

1. **HeLa is the wrong system.** Already aneuploid, p53-null, fast-cycling, and adapted. Its
   tolerance of BubR1 loss need not resemble a diploid cell's.
2. **The residual protein in patients is qualitatively defective**, not merely scarce. Several MVA
   alleles are missense.
3. **The readouts differ again** — patient cells were scored on nocodazole arrest and alignment,
   the HeLa titration on live-imaging missegregation.

Whichever holds, it undercuts using a HeLa titration to set a threshold for the disease.

---

## 4. Is the threshold readout-specific? Yes.

Ranked by how much BubR1 you can lose before the readout moves — most sensitive first.

| Readout | Where it breaks | Evidence |
|---|---|---|
| **PCS / centromeric cohesion** | already at ~50–77% residual | Bohers 2008: "Induction of PCS was detected in **all** transduced cells" — the least-depleted clone was 77%. García-Castillo 2008: PCS in 7/7 human heterozygous carriers. Wijshake 2012: PCS elevated in +/GTTA MEFs. |
| **Accumulated aneuploidy** | somewhere below 58%; clear at ≤14% | Bohers 2008. Wijshake 2012: +/GTTA MEFs 12% → ~18% aneuploid, "relatively small" against H/H. |
| **Per-division missegregation** | below ~13%, gross at ~6% | Suijkerbuijk 2010, Fig. 4D lower. |
| **Spindle assembly checkpoint (nocodazole arrest)** | **not titrated by anyone I found** | Suijkerbuijk 2010 tested SAC only in mutant-replacement, not against a dose series. |
| **Attachment stability / error correction** | **not titrated** | — |
| **PP2A-B56 recruitment via KARD** | **not titrated** | — |
| **Cortical development (mouse)** | tolerates ~1.8%; fails at ~0.5% | Simmons 2019. |

Two further reasons "the threshold" cannot be a single number:

**The denominator is wrong.** Whole-cell immunoblot percentage is not the functional variable.
Suijkerbuijk 2010 quantified kinetochore-bound BubR1 separately and found residual patient protein
still loads onto kinetochores. Simmons 2019 makes the same point in mouse: "the residual BubR1
protein was localized to the chromosomes in BubR1^H/H, but largely absent in BubR1 CKO", and they
attribute the hypomorph's tolerance to exactly that. A cell at 10% whole-cell BubR1 with normal
kinetochore loading is not in the same state as a cell at 10% with none.

**Level does not determine phenotype.** Wijshake 2012's +/GTTA versus +/− comparison, above.

---

## 5. Is it cell-type specific? Yes — and the disease tissues are the worst covered.

**The same allele gives different levels in different tissues.** The BubR1 hypomorph is cited at
~11% in MEFs (Baker 2004, via Simmons 2019) but measured at ~1.8% in E15.5 cortex (Simmons 2019).

**Tolerance also differs.** BubR1^H/H mice are dwarfed, progeroid, aneuploid and tumour-prone, yet
Simmons 2019 found cortical size and progenitor division normal, and concluded:

> "Taken together, this data supports the idea that substantial loss, approaching complete absence,
> of BubR1 in cortical progenitors is required to induce microcephaly in mice."

and

> "Of note, cortex size is only decreased when BubR1 level is further eliminated by Cre-mediated
> cortex specific deletion."

Caveats, which matter:

- Mouse BubR1^H/H does not get microcephaly; human MVA patients do. The mouse cortex result may
  simply not transfer.
- Distinguishing 1.8% from 0.5% on a western is near the noise floor of densitometry. Two labs,
  two antibodies, two blots.
- Simmons 2019 themselves invoke "functional compensation" in the chronic hypomorph background —
  adaptation over development, which an acute knockdown would not reproduce.
- The same H/H allele *does* impair adult hippocampal neurogenesis (Yang 2017), so "brain tolerates
  it" is too coarse. Embryonic cortical neurogenesis and adult hippocampal neurogenesis differ.

**Growth plate: I did not find any study.** Not one paper measured BubR1 level against a segregation
or proliferation readout in chondrocytes or growth plate. Queries in §7.

---

## 6. What our document should say instead of "~10–60%"

Replace the range with a statement about the *shape* of the relationship and the *identity* of the
gap. Proposed sentence for the boxed hypothesis in §4 and for L03:

> **"Residual BubR1 and chromosome missegregation are related by a graded dose–response, not a
> step, and no published experiment locates a threshold: Bohers 2008 (PMID 18932004) scored
> premature chromatid separation and accumulated aneuploidy in chronically depleted HeLa clones
> with no data point between 14% and 58% residual, while Suijkerbuijk 2010 (PMID 20516114) scored
> per-division missegregation by live imaging in the same cell type and found ~13% residual had
> 'virtually no impact', with gross missegregation only at ~6% — the two never measured the same
> quantity, so the apparent 4× disagreement is not evidence about where a threshold sits."**

If one sentence is needed for the hypothesis box, use the shorter form:

> **"…if it raises residual BubR1 far enough to move the dose–response — whose position nobody has
> measured in a disease-relevant cell type, and which differs by readout, with cohesion (PCS)
> affected at ~50–77% residual and per-division missegregation only below ~13%."**

Two supporting edits the ledger should carry:

1. Drop "~50% residual threshold" wherever it survives. It is an interpolation across an unsampled
   gap in one paper, restated by us as a measurement.
2. Add the number that actually bears on the therapeutic claim: **MVA patient cells sit at roughly
   17–50% of control full-length BubR1** (Suijkerbuijk 2010, Fig. 2A, immunoblot BUBR1/tubulin
   against healthy controls). That is the level our intervention must move, and it sits inside the
   region nobody titrated.

---

## 7. The experiment that would actually locate it

### System

Endogenous degron knock-in at *BUB1B* in **hTERT-RPE1** — diploid, p53-proficient,
karyotypically stable, baseline missegregation ~1%. Both alleles tagged with an auxin-inducible
degron (AID2: OsTIR1(F74G) + 5-Ph-IAA) or dTAG.

Why not repeat the HeLa experiments: HeLa is already aneuploid, so its baseline swamps the effect
being measured; stable shRNA clones confound dose with clonal drift; and a degron is continuous,
fast, reversible, and does not select.

### Dosing

At least eight levels spanning 5–100% residual, with spacing ≤10 percentage points **across the
14–58% region** — the region neither paper sampled and where MVA patient cells actually sit.

### Quantifying residual — two measures, not one

1. Endogenous tag fluorescence, per cell, by quantitative live imaging.
2. **Kinetochore-bound BubR1** by immunofluorescence normalised to ACA or CENP-C.

Report both. Suijkerbuijk 2010 and Simmons 2019 both indicate kinetochore-loaded protein, not
whole-cell abundance, is the functional quantity.

### Readouts, all in the same cells at every level

1. Per-division missegregation: live H2B imaging, lagging chromosomes and micronuclei, **≥500
   divisions per level** — powered to detect 1% → 3%.
2. SAC strength: mitotic arrest duration in nocodazole and in low-dose taxol.
3. Attachment stability: cold-stable K-fibres; interkinetochore distance.
4. PP2A-B56 recruitment: B56 kinetochore intensity; KARD phosphorylation.
5. PCS: metaphase spreads, scored blind.
6. Accumulated aneuploidy: single-cell sequencing after N passages held at a fixed level — reported
   **separately** from (1), because conflating these is what produced the current confusion.

Fit a dose–response per readout and report EC50 with confidence intervals. That converts "the
threshold" into six numbers, which is what it actually is.

### The arm that tests our claim

A degron only lowers BubR1. Our hypothesis is about **raising** it. So add: doxycycline-titratable,
shRNA-resistant *BUB1B* re-expression in MVA patient-derived cells against an isogenic corrected
control, stepping residual up from patient baseline, scored on the same six readouts. Measuring the
direction we propose to move is not optional; the descending and ascending curves need not
coincide, because Simmons 2019's "functional compensation" argument implies chronic adaptation.

### Disease-tissue arm

Repeat readouts 1, 2 and 6 in iPSC-derived neural progenitors from an MVA patient and an isogenic
corrected control. Nobody has done this in growth plate; chondrocyte differentiation from the same
iPSC line would be the first such measurement.

---

## 8. Searches run

PubMed E-utilities (esearch → esummary/efetch), Europe PMC, PMC BioC full text, Unpaywall,
Semantic Scholar, HAL, OMIM, WebSearch/WebFetch.

### Productive

| Query | Count | Kept |
|---|---|---|
| `BUBR1 dosage aneuploidy` | 4 | 18932004 |
| `BUB1B hypomorph` | 8 | 18932004, 24825348, 27313317 |
| `BUBR1 levels premature chromatid separation` | 7 | 16411201, 18932004, 21389084 |
| `BUB1B expression level chromosome missegregation` | 4 | 19117986 |
| `BubR1 neural progenitor` | 6 | 30668728, 37900274, 28383136, 27848932 |
| `Bub1b +/- mice aneuploidy tumor` | 13 | 23300461, 15767571 |
| `BubR1 haploinsufficient mice aneuploidy` | 6 | 24981203, 16355229 |
| `BUB1B carrier parents premature chromatid separation` | 1 | 18548531 |
| `BubR1 dwarfism mice` | 3 | 15208629 |
| `mitotic checkpoint complex BubR1 concentration model` | 6 | 26256776, 35767360 |
| `BUB1B[Title] AND CRISPR` | 1 | 34545188 (not a titration) |
| citations of PMID 18932004 via Europe PMC | 20 | none reproduce its internal numbers |

### Found nothing — exact queries

These returned **zero** PubMed records. I did not find the evidence; I am not claiming it does not
exist.

| Query | Count |
|---|---|
| `BUBR1 auxin degron` | 0 |
| `BUB1B degron` | 0 |
| `BubR1 chondrocyte` | 0 |
| `BubR1 chondrocyte growth plate` | 0 |
| `BUB1B osteochondro` | 0 |
| `BubR1 cell type specific sensitivity` | 0 |
| `aneuploidy tolerance tissue specific neural progenitor` | 0 |
| `BUB1B heterozygous knockout human cell line` | 0 |
| `BUB1B HCT116 knockout` | 0 |
| `BubR1 KARD PP2A B56 depletion level` | 0 |
| `BubR1 cohesion protection shugoshin PP2A premature chromatid` | 0 |
| `partial knockdown spindle checkpoint dose response RNAi` | 0 |
| `absolute quantification spindle assembly checkpoint proteins copy number` | 0 |
| `mitotic checkpoint complex stoichiometry quantitative proteomics` | 0 |
| `BubR1 excess over Cdc20` | 0 |
| `BubR1 growth plate` | 2 — neither about BubR1 dosage (RSK2 mitosis; Arabidopsis SAC) |

### Full text I could not obtain

| Item | Status |
|---|---|
| **Bohers 2008, PMID 18932004** | **Unresolved — paywalled.** Not in PMC (`inPMC: N`). Europe PMC: "Subscription required". Unpaywall: `is_oa: false`, zero OA locations. HAL (author is INSERM Rouen): 0 hits for both author and title queries. Springer direct: 303 redirect to `idp.springer.com`, then a 3 KB block page. OMIM 602860: HTTP 403. **Consequence: how Bohers quantified residual BubR1 (which blot, which loading control), how PCS and aneuploidy were scored, how many metaphases, and the per-clone numbers are all unverified.** Everything I state about Bohers comes from the abstract. |
| Baker 2004, PMID 15208629 | Paywalled. The ~10–11% figure for BubR1^H/H is taken from secondary citations in Wijshake 2012 and Simmons 2019, not from the primary blot. |
| Suijkerbuijk 2010 supplement | Not in the PMC author manuscript. Doxycycline exposure duration for the titration is **unresolved**. |
| Ochiai 2014 Fig. 3A | Densitometry values are printed inside the figure image. Residual % **unresolved**. |
| Wijshake 2012 Table S1 | BubR1^H/H aneuploidy rate not in the accessible text. |

---

## 9. What changed against the current ledger

- `CL-11` / `WK-11`: the weakening was right but the reason was incomplete. The problem is not only
  the 14–58% gap. It is that Bohers and Suijkerbuijk measure different quantities, so citing them as
  two estimates of one threshold is a category error.
- `L03` should stop describing a threshold and start describing a dose–response with readout-specific
  slopes.
- New, and load-bearing for the therapeutic claim: MVA patient cells sit at **≈17–50% of control**
  (Suijkerbuijk 2010 Fig. 2A). This is higher than the level at which HeLa tolerates depletion, which
  is an argument against reasoning about this disease from HeLa titrations at all.
- New: PCS is not a pure function of BubR1 level (Wijshake 2012, `+/GTTA` vs `+/−`).
- New: mouse cortex tolerates ~1.8% residual without microcephaly (Simmons 2019) — a tissue-specific
  result in the exact tissue our document cares about, pointing the opposite way from a 50% threshold.
