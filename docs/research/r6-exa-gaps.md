# R6 — Exa semantic sweep: the gaps four PubMed passes could not see

Status: COMPLETE. Written incrementally; every section below is sourced.

Search engine: Exa semantic web search (`web_search_exa`, `web_fetch_exa`), not PubMed
E-utilities. Rationale: four prior passes (R1–R5) used PubMed keyword matching over
indexed journals only, and several of their **zero-hit results are now load-bearing
negative claims** in the submission. Exa is embedding-based and crawls preprints
(bioRxiv/medRxiv), theses, protocol repositories, single-cell atlases, database pages
and lab sites — the exact surface PubMed does not cover.

**Rule applied throughout:** every negative below is stated as "I did not find", never
as "does not exist", and every negative carries the exact queries that failed. §6 lists
all queries run, including the empty ones.

Evidence grades used: **peer-reviewed** / **preprint** / **database record** /
**unpublished-observation-in-a-published-paper** (the Cole "data not shown" case, which
is why the grade exists).

---

## 1. Verdict in three sentences

The neural arm survives, but not in the form the document currently states it: I did not
find any measurement supporting "neural progenitors are probably NAPRT-negative", and two
independent human single-cell atlases show *NAPRT* transcript in fetal radial glia and
intermediate progenitors at levels **above the ambient-noise floor and at least as high as
in neurons** — so the claim is contradicted in direction at transcript level and simply
unmeasured at protein and enzyme-activity level. The named probe does not change: nothing
found here displaces nicotinic acid, and the one paper that looks like it should (Huang et
al., *Cells* 2022, PMID 35455963) tested **1–2 mM extracellular NAD⁺ itself on mouse adult
hippocampal progenitors**, which is neither the molecule, the concentration, the species,
nor the developmental compartment our hypothesis proposes. What *does* change is that the
document now owes a real direction-of-effect conflict: NAD⁺ elevation is pro-proliferative
in neural progenitors on the NAMPT-genetic evidence and anti-proliferative on the
exogenous-NAD⁺ evidence, the two disagree, and neither was run with a precursor.

---

## 2. NAPRT in neural progenitors — the Cole extrapolation does not hold

**Grade of the claim we are replacing.** Cole J et al., *Oncotarget* 2017, PMID 29100430 is
peer-reviewed, but the specific sentence our §4.5 leans on is an
**unpublished-observation-in-a-published-paper**. Verbatim from the paper's Results:

> "However, in non-tumorigenic brain cortex, neurons were NAPRT negative while reactive
> astrocytes seen in some of the samples were positive (data not shown). In samples of
> non-tumorigenic white matter, oligodendrocytes were also NAPRT negative (Figure 5C)."

The paper's own discussion calls the sample "a small number of samples" and says "further
study is warranted to evaluate different types of normal neuronal cells within the brain."
It measured **mature neurons and oligodendrocytes in adult cortex and white matter**. There
is no progenitor, no radial glia, no fetal tissue, and no organoid in it.

**Verdict: the "progenitors are NAPRT-negative" claim is CONTRADICTED in direction at
transcript level, and UNMEASURED at protein and activity level.**

### 2.1 New primary evidence — two human single-cell atlases, queried directly

I did not find a paper that measures NAPRT in human neural progenitors, so I queried the
atlases themselves. Both are **database records** (public expression matrices), not
peer-reviewed claims about NAPRT; the peer review covers the atlas, not this analysis.

**Atlas A — first-trimester whole human brain, 691,928 cells.**
Braun E et al., *Science* 2023, PMID 37824650 — droplet 10x, 5–14 postconceptional weeks.
Retrieved from the UCSC Cell Browser mirror `dev-brain-regions/wholebrain`
(https://cells.ucsc.edu/?ds=dev-brain-regions+wholebrain), raw UMI matrix, queried
2026-08-29. The published atlas reports 1,665,937 cells; the whole-brain browser subset I
analysed contains 691,928, so these are the browser's numbers and not the paper's. Percent
of cells with ≥1 UMI, and mean counts per 10,000:

| cell type | n | NAPRT det. | NAPRT CP10K | NAMPT det. | NAMPT CP10K | QPRT det. | SIRT2 det. |
|---|---|---|---|---|---|---|---|
| radial glia | 57,028 | **6.99%** | 0.17 | 14.72% | 0.39 | 13.57% | 23.76% |
| intermediate progenitor cell | 9,655 | **4.34%** | 0.11 | 8.43% | 0.21 | 4.10% | 15.05% |
| dividing | 74,949 | **5.23%** | 0.12 | 10.40% | 0.24 | 8.10% | 14.79% |
| neuron | 329,685 | 5.05% | 0.16 | 9.95% | 0.32 | 1.18% | 15.60% |
| interneuron | 148,332 | 2.42% | 0.10 | 6.87% | 0.29 | 1.82% | 13.34% |
| astrocyte | 27,763 | 5.62% | 0.18 | 10.75% | 0.34 | 12.58% | 23.64% |
| oligodendrocyte | 21,498 | 5.10% | 0.19 | 8.96% | 0.31 | 6.54% | 47.31% |
| **endothelial** | 9,349 | **17.22%** | 0.31 | 47.76% | 2.20 | 22.07% | 30.08% |
| microglia | 9,626 | 6.26% | 0.21 | 14.21% | 0.54 | 14.47% | 17.37% |

**The noise floor, which is what makes this readable.** Detection percentages in shallow
droplet data are meaningless without a floor. Genes that should be absent from first-
trimester brain, computed on the same cells: *INS* 0.00%, *MYOD1* 0.00%, *CDX2* 0.00%,
*KRT14* 0.00%, *AMY2A* 0.01%, *ALB* 0.03%, *CD19* 0.15%. Ceiling: *ACTB* 99.04%, *GAPDH*
94.77%. **NAPRT in radial glia sits at 6.99% — roughly 50–200× the ambient floor.** It is
low, but it is not zero, and it is not noise.

Three things follow, and the third is the one that matters:

1. NAPRT is expressed in human fetal radial glia and IPCs, at roughly **half the detection
   rate of NAMPT** in the same cells. The Preiss-Handler entry step is present, and the
   nicotinamide-salvage entry step is bigger. That is a statement about relative flux
   capacity, not about the presence or absence of a route.
2. **Endothelium is the highest NAPRT cell type in the atlas** (17.2%), which independently
   reproduces the one positive Cole reports at protein level. The rank order agrees; the
   binary "neurons are negative" does not survive being read as "the neural lineage is
   negative."
3. **Radial glia detect NAPRT more often than neurons do** (6.99% vs 5.05%) and much more
   often than interneurons (2.42%). The extrapolation in §4.5 runs from mature neurons to
   progenitors, and the data run the other way: the progenitor is the *better*-expressing
   cell of the two, not the worse.

**Atlas B — developing human cortex, 4,261 cells, an independent platform.**
Nowakowski TJ et al., *Science* 2017, PMID 29217575 — Fluidigm C1 full-length SMART-seq,
TPM. Retrieved from `cells.ucsc.edu/?ds=cortex-dev`, queried 2026-08-29. Mean TPM by WGCNA
cluster, with percent of cells detecting *NAPRT1*:

| cluster | n | NAPRT1 mean TPM | NAPRT1 % cells | NAMPT mean TPM | SIRT2 | BUB1B |
|---|---|---|---|---|---|---|
| vRG (ventral radial glia) | 97 | 0.06 | 1.0% | 115.4 | 64.8 | 12.5 |
| oRG (outer radial glia) | 65 | 1.81 | 3.1% | 254.9 | 66.6 | 12.2 |
| tRG (truncated radial glia) | 100 | 1.19 | 4.0% | 418.4 | 70.8 | 13.0 |
| RG-early | 53 | 0.04 | 3.8% | 91.4 | 0.35 | 68.3 |
| RG-div1 / RG-div2 | 133 / 124 | 0.04 / 0.84 | 0.8% / 2.4% | 221.8 / 185.4 | 42.6 / 18.1 | 290.6 / 71.2 |
| IPC-nEN1/2/3 | 68 / 76 / 39 | 1.83 / 0.64 / 0.13 | 1.5% / 3.9% / 5.1% | 98.8 / 52.4 / 46.1 | 38.2 / 31.3 / 33.5 | 4.0 / 15.3 / 0.0 |
| IPC-div1 / IPC-div2 | 36 / 51 | 0.00 / 0.00 | 0% / 0% | 75.9 / 83.5 | 6.0 / 6.5 | 423.6 / 121.9 |
| Astrocyte | 46 | 0.10 | 2.2% | 285.9 | 41.3 | 3.7 |
| **Endothelial** | 30 | 0.36 | **13.3%** | 1097.3 | 47.6 | 0.1 |
| whole dataset | 4,261 | 0.36 | 1.9% | 165.5 (32.2% det.) | — | — |

This dataset is smaller, older and much sparser for low-abundance transcripts, and it puts
NAPRT1 lower — 1.9% of all cells versus NAMPT's 32.2%. **It agrees with Atlas A on rank
order** (endothelium top at 13.3%; NAMPT ≫ NAPRT everywhere) and disagrees on absolute
level. Where two platforms disagree on magnitude and agree on order, the honest reading is
the order.

Two by-products worth keeping. *SIRT2* — the deacetylase the whole mechanism runs
through — is expressed across every progenitor cluster in both atlases (23.8% of radial
glia in Braun; 42–71 TPM in vRG/oRG/tRG in Nowakowski, and 910 TPM in OPCs), so the
**effector arm of the mechanism is present in the compartment that produces the
microcephaly**. And *BUB1B* tracks the dividing clusters exactly as it should (RG-div1
290.6, IPC-div1 423.6, MGE-IPC1 395.5 TPM against ~0–8 TPM in postmitotic neurons), which
is a small independent confirmation that the gene our disease is about is a progenitor gene
in human cortex, not a neuronal one.

### 2.2 The one published paper that measures NAPRT protein during neural differentiation

Neves D et al., *Neurochem Int* 2022, PMID 35843422 — peer-reviewed. SH-SY5Y
neuroblastoma differentiated with retinoic acid then BDNF; NAD-biosynthetic enzymes
measured across the time course. Two findings bear directly on §4.5:

- **NAPRT protein *increases* after 5 days of retinoic acid.** Verbatim: "Differences were
  found in the expression of NAD biosynthetic enzymes during neuronal differentiation,
  namely, increased QPRT gene expression along the course of RA/BDNF treatment and NAPRT
  protein expression after a 5-day treatment with RA." NAPRT is not switched off as a
  neuroblast differentiates in this model; it goes up.
- **But NAPRT is not the flux-carrying enzyme there.** "NAPRT inhibition had no effect in
  neuritic length per cell, NAD+ levels and HDAC activity", whereas QPRT inhibition dropped
  NAD⁺ and HDAC activity. The authors' own conclusion is that NAPRT has "additional roles
  beyond NAD production in RA-differentiated cells."

This is a neuroblastoma line, not a primary progenitor, and the NAD⁺ result is a knock-down
of basal flux rather than a nicotinic-acid challenge — so it does not establish that niacin
would raise NAD⁺ in a neural progenitor. It does establish that *the enzyme is there and is
developmentally regulated upward*, which is the opposite of the picture §4.5 currently
paints.

### 2.3 Two supporting strands, both weaker, both pointing the same way

Duarte-Pereira S et al., *Genes* 2021, PMID 34946971 (bioinformatics, peer-reviewed; no
wet-lab NAPRT measurement) reports that the transcription factors, RNA-binding proteins and
miRNAs predicted to regulate *NAPRT* are enriched for "cell differentiation, development and
neuronal biology", and notes NAPRT variants associated with ADHD and schizophrenia and a
zebrafish partial-loss-of-function causing abnormal brain development. This is
hypothesis-generating only — I list it as inference, not evidence.

Duarte-Pereira S et al., *Oncotarget* 2016, PMID 26675378 is the paper §4.5 already cites
for "widely expressed". Its brain-specific sentence cuts the other way and should be quoted
rather than omitted: "Curiously, the brain tissue appeared to express less NAMPT and
NAPRT", and "these results are in agreement with previous findings that the main routes to
NAD in neurons are the de novo biosynthesis from tryptophan and the NAD salvage from
nicotinamide riboside." So the brain is genuinely the low end of the atlas — but "low" and
"absent" are different claims, and the disease's cell type is not a neuron.

### 2.4 What is still missing, stated as a gap

I did not find, in any source: NAPRT **protein** or **enzyme activity** measured in human
neural stem cells, radial glia, iPSC-derived neuroepithelium or cerebral organoids; any
nicotinic-acid NAD⁺-rescue experiment in a human neural progenitor; or any NAMPT-inhibitor
plus niacin rescue run in a neural progenitor. Queries in §6.1. **The operational
consequence is unchanged in shape but changed in expectation: the NAPRT gate stays, but it
is now a gate we should expect to pass, not one we should expect to fail.**

---

## 3. Does raising NAD⁺ harm neural progenitors? A real conflict, but not the one it looks like

### 3.1 The paper, read closely

Huang X et al., *Cells* 2022, PMID 35455963 — peer-reviewed; RNA-seq deposited as GEO
GSE200494; Xuekun Li lab, Zhejiang University. Everything below is from the paper.

| | |
|---|---|
| **Species / cells** | **Mouse.** Primary adult neural stem/progenitor cells (aNSPCs) isolated from adult mouse hippocampus, plus primary mouse hippocampal neurons. No human cells, no fetal cells. |
| **Agent** | **β-nicotinamide adenine dinucleotide hydrate — NAD⁺ itself, added to the culture medium.** Not nicotinic acid, not nicotinamide, not NR, not NMN. No precursor was tested. |
| **Dose** | A dose series was screened by BrdU; **1 mM and 2 mM** were selected as showing "the most significant effects" and used throughout. Metabolic assays at 1 mM for 48 h. |
| **Findings** | ↓BrdU⁺ and ↓Ki67⁺; ↓neurosphere size at d7 and d14; ↓Nestin and ↓Sox2; ↓Tuj1⁺ neurons with ↑GFAP⁺ astrocytes and ↑*S100β*/*Glt1*; ↑cleaved caspase-3 by IF and immunoblot; impaired neurite morphology; ↑ATP; ↑NADP⁺/NADPH; **↓p-Akt(Ser473) with total Akt unchanged**; ↓*Pdgfrb*, *Epor*, *Erbb3*, *Il2rb*. |
| **Authors' own limitation** | Verbatim: "One limitation of our present study was that it was largely carried out **in vitro**, and an intensive in vivo study would be helpful to clarify this speculation." And: "Under pathological conditions, such as NAD+ depletion, the addition of NAD+ can attenuate neurotoxicity … However, under **physiological conditions**, the supplement of extra NAD+ could be detrimental." |

So the finding is real and it is peer-reviewed. The question is what it transfers to.

### 3.2 Transferability: weak, on four grounds, and the fourth is decisive

**(a) Wrong molecule and wrong side of the membrane.** Our hypothesis raises *intracellular*
NAD⁺ by feeding a precursor through NAPRT. Huang added the *dinucleotide itself to the
medium*. NAD⁺ is a large charged dinucleotide; the extracellular pool is a distinct
signalling compartment consumed by ectoenzymes (CD38, CD73/NT5E, ecto-ARTs), and raising it
is not the same experiment as raising the intracellular pool. Huang measured no
intracellular NAD⁺ concentration in the abstract or the results I read — the reported
metabolic readouts are ATP and the NADP⁺/NADPH ratio.

**(b) Wrong concentration, by about three orders of magnitude.** 1–2 mM extracellular. The
exposure our own document targets is nicotinic acid at ~15 µM ex vivo against a ~75 µM peak
plasma after 2 g extended-release niacin (§4.5). Nothing in the clinical range approaches
millimolar extracellular NAD⁺.

**(c) A competing explanation the paper did not exclude — and this is my inference, not a
claim either paper makes.** Millimolar extracellular NAD⁺ in *mouse* cells is a known
pharmacological stimulus in its own right. Leeson HC et al., *Stem Cells* 2018, PMID
30068016 showed that **adult murine hippocampal (SGZ) and SVZ neural progenitor cells — the
same cell type and same species Huang used — express functional P2X7**, and that P2X7
agonism produces a *dose-dependent decrease in NPC proliferation*, with transmembrane pore
formation and cell death at high agonist (significant pore formation at 500 µM and 1 mM
ATP). Extracellular NAD⁺ activates P2X7 in mouse by ARTC2.2-catalysed ADP-ribosylation of
the receptor ectodomain, beginning at low micromolar NAD⁺ and, unlike ATP, essentially
irreversible (Hong S et al., *J Immunol* 2009, PMID 19542469; and the
mechanism review in Haag F et al., *Purinergic Signal* 2007, PMID 18404420). Huang tested neither P2X7 nor
ARTC2, and used no P2X7 antagonist and no NAD⁺-analogue control. **A reduction in
proliferation plus caspase-3 induction at 1–2 mM extracellular NAD⁺ in mouse NPCs is exactly
what the P2X7 literature predicts on grounds that have nothing to do with NAD⁺ metabolism.**

**(d) The mouse-specific arm of that mechanism does not exist in humans.** Glowacki G et al., *Protein Sci* 2002, PMID 12070318: "The human genome contains four functional toxin-related
ADP-ribosyltransferase genes (ARTs) and two related intron-containing pseudogenes; the mouse
has six functional orthologs" — and **ART2 is one of the human pseudogenes**, while the
mouse has functional *Art2a* and *Art2b*. The NAD⁺→ARTC2.2→P2X7 death pathway is a mouse
pathway. If that is what Huang observed, it does not transfer to a human patient at all.

**(e) Wrong developmental compartment.** aNSPCs are adult hippocampal progenitors. MVA
microcephaly is prenatal and driven by cortical radial glia and IPCs, which is a different
niche, a different cell cycle regime and a different metabolic state.

### 3.3 The evidence pointing the other way, which is stronger

- **Genetic, in vivo, opposite sign.** Stein LR and Imai S, *EMBO J* 2014, PMID 24811750:
  ablating *Nampt* in adult NSPCs "reduced their pool and proliferation in vivo"; Nampt is
  "required for G1/S progression of the NSPC cell cycle"; and critically, "**the decrease in
  the NSPC pool during aging can be rescued by enhancing hippocampal NAD+ levels.**" Same
  cell type, same species, same compartment as Huang, and NAD⁺ elevation *rescues*.
  Same paper also finds NSPC oligodendrocytic fate decisions run "through a mechanism
  mediated redundantly by **Sirt1 and Sirt2**" — so the effector our hypothesis names is a
  live NAD⁺-dependent node in neural progenitors, not a borrowed one.
- **Pharmacological, in vivo, opposite sign.** Zhao Y et al., *Stroke* 2015, PMID 26060246:
  NAMPT-transgenic mice show more neural stem cells and better recovery after MCAO, absent
  in the enzymatically-dead H247A transgenic; NMN and NAD promoted NSC proliferation and
  differentiation in neurosphere culture while the NAMPT inhibitor FK866 did the opposite.
  Reviewed in Wang SN et al., *CNS Neurosci Ther* 2016, PMID 27018006.
- **The single most transferable developmental datum, and it is about niacin specifically.**
  Shi H et al., *N Engl J Med* 2017, PMID 28792876: in *Haao*-null and *Kynu*-null mouse
  embryos, NAD deficiency caused vertebral, cardiac, renal, palatal and limb defects, and
  **maternal nicotinic acid supplementation during gestation prevented them, with a
  dose-dependent rise in embryonic NAD(H)** — 10 mg/L nicotinic acid in drinking water left
  only 30%-small kidneys, 15 mg/L produced entirely normal embryos. This is nicotinic acid,
  given in vivo, during organogenesis, raising embryonic NAD⁺, and it is *protective*. It
  is not a microcephaly or progenitor-proliferation endpoint, so it does not settle our
  question — but it is the closest available match to the intervention we propose in the
  developmental window we care about, and its sign is favourable.
- **Huang's own introduction concedes the precursor case.** Verbatim from the paper:
  "the supplementation of NR increased the level of NAD+, promoted the proliferation of
  neural progenitor cells and neurogenesis and improved the cognitive function of
  neurodegenerative animal models."

### 3.4 A second exogenous-NAD⁺ toxicity, in human cells, that sharpens the pattern

Neves D et al., *Neurochem Int* 2022, PMID 35843422 (already cited in §2.2): "Of note, NAD
supplementation along with RA, but not with BDNF, resulted in considerable cell death." Human
SH-SY5Y, exogenous NAD⁺ again, cell death again. Two independent groups adding NAD⁺ to the
medium of neural-lineage cells and getting death is a **consistent finding about exogenous
NAD⁺ as a compound**, not evidence that a raised intracellular NAD⁺ pool is harmful. Every
study that raised NAD⁺ through a *precursor or a genetic route* in neural progenitors found
the opposite sign.

### 3.5 What this actually costs the document, stated plainly

The document should not claim the direction-of-effect conflict is resolved. It should carry
it as follows:

> The NAD⁺/neural-progenitor literature is **sign-inconsistent**. Genetic NAMPT loss depletes
> the progenitor pool and NAD⁺ repletion rescues it (Stein and Imai 2014); NAMPT
> overexpression and NMN expand it (Zhao 2015); gestational nicotinic acid prevents
> NAD-deficiency malformations (Shi 2017). Against that, millimolar **exogenous** NAD⁺
> reduces proliferation and induces apoptosis in mouse adult hippocampal NSPCs (Huang 2022)
> and kills RA-differentiating SH-SY5Y (Neves 2022). The split falls cleanly along
> *exogenous dinucleotide versus precursor/genetic*, and the exogenous arm has an unexcluded
> mouse-specific P2X7/ARTC2 explanation. **No study has given a NAD⁺ precursor to a human
> neural progenitor and measured proliferation.** The neural arm therefore carries an
> unmeasured risk of an anti-proliferative effect in the exact compartment that produces the
> microcephaly, and the falsifying experiment should add a progenitor-proliferation readout
> to the neural arm rather than assume the sign.

That is a genuine liability and it is not currently in §7. It is *not*, however, a reason to
change the probe: nothing here is evidence against nicotinic acid, and the one nicotinic-acid
developmental experiment in the literature points the other way.

---

## 4. Priority 3 — the four PubMed zeros, re-run semantically

### 4.1 "No prior work connects this disease to this pathway" — that is a retrieval artifact

PubMed `"mosaic variegated aneuploidy" AND (NAD OR nicotinamide OR niacin)` returns 0
because PubMed matches title, abstract and MeSH. The connection is made in **body text**,
twice, and Exa finds it.

**The founding paper of our own axis names MVA in its Introduction.** North BJ et al., *EMBO J* 2014, PMID 24825348, verbatim:

> "In humans, mutations in BubR1 have been associated with mosaic variegated aneuploidy
> (MVA) syndrome. MVA patients have various progeroid traits including a shortened
> lifespan, short stature, facial dysmorphisms, and cataracts (Jacquemont et al, 2002; Lane
> et al, 2002; Garcia-Castillo et al, 2008)."

**A 2024 review from the same senior author proposes the axis as a therapeutic target.**
Pun R et al., *Semin Cancer Biol* 2024, PMID 39490401 (peer-reviewed review; PMCID
PMC11625622). It names MVA and microcephaly; it states that "Supplementation of aged mice
with nicotinamide mononucleotide (NMN), to increase NAD+, restores BubR1 levels to those
observed in younger tissues," and frames SIRT2→BubR1 as "a promising avenue for developing
targeted therapies." Full-text keyword scan: **`nicotinic acid` 0 hits, `niacin` 0 hits,
`NAPRT` 0 hits**; NMN appears 23 times, nicotinamide riboside 3 times.

**What this costs and what it does not.** §4.5's line "no prior work connects this disease
to this pathway" should be withdrawn — the axis's own literature connects them, and a
reviewer who knows this field will know it. What survives is narrower and still true as far
as I could find: **I did not find any study that administered an NAD⁺ precursor, a sirtuin
activator, or any NAD-pathway agent to MVA patient cells, MVA patient-derived lines, or an
MVA animal model and measured a mitotic-fidelity endpoint.** And the specific probe is
untouched: the review that comes closest to our hypothesis reaches for NMN and NR, and never
mentions nicotinic acid or NAPRT. Queries in §6.3.

For completeness on the field's own view: Malumbres M and Villarroya-Beltri C, *Nat Rev
Genet* 2024, PMID 39169218, the current MVA review, states "As is the case for most
developmental syndromes, no therapies are available for preventing or limiting
MVA-associated pathologies," and proposes future directions around limiting CIN and
modulating the immune response. **No NAD⁺, sirtuin or BubR1-stabilisation proposal appears
in it.** It also carries a threshold claim worth having: "small quantities (less than 5% of
normal levels) of SAC components, such as BUB1 or CDC20, may be sufficient to support
chromosome segregation."

### 4.2 A quantified residual fraction we did not have, and it sits below the interval

R3 has Suijkerbuijk 2010 at ≈17–50% residual. Exa surfaces a second human data point, lower
and phenotypically decisive.

Rio Frio T et al., *N Engl J Med* 2010, PMID 21190457 — a 68-year-old man, homozygous
*BUB1B* `c.2386-11A>G`. As summarised in **OMIM entry #257300** (database record; queried
2026-08-29): the mutant transcript is degraded by NMD, but "a small amount (**10 to 15%**)
of normal BUB1B was produced and correctly localized to the kinetochores of patient
fibroblasts," while "the residual amount of protein was unable to maintain the
spindle-assembly checkpoint." Premature chromatid separation in **57–84%** of cells, mosaic
variegated aneuploidy, adenocarcinoma of the ampulla of Vater at 34 and later multiple
colon and gastric primaries — **and "no other features of the mosaic variegated aneuploidy
syndrome, such as poor growth, microcephaly, or mental retardation."**

Two caveats, stated because the last time a number was carried forward unverified it was
wrong. The 10–15% figure as I read it comes from **OMIM's summary of the paper, not from
the paper**; I did not obtain the NEJM full text, so the blot, the loading control and the
n behind it are unverified. And this is a splice allele, so it is one point in the allelic
series, not the series.

**Why it matters more than its provenance suggests.** At 10–15% residual full-length BubR1
— *below* the 17–50% band our document uses — a human had near-total PCS and adult cancers
and reached 68 years old with a normal head and normal cognition. That points the same way
as Simmons 2019 in mouse cortex (~1.8% tolerated without microcephaly, already in R3): in
humans as in mice, **the brain appears to tolerate a lower residual BubR1 than the
cancer-suppression function does**. If that holds, the neural arm is not the arm with the
smallest headroom — it is the arm with the largest, and the therapeutic case should be made
on missegregation and tumour risk rather than on head circumference. This also sits
awkwardly against Matsuura S et al., *Am J Med Genet A* 2006, PMID 16411201, whose
conclusion as summarised in OMIM is that "a greater than 50% decrease in BUB1B activity
results in abnormal mitotic spindle checkpoint function and MVA syndrome." The two cannot
both be general.

**I did not find** any measurement of BubR1 protein in MVA patient-derived cells published
after Suijkerbuijk 2010 that adds a quantified residual fraction. Queries in §6.3.

### 4.3 A BubR1 threshold by degron or titration: I did not find one, but the neighbouring
experiment exists and is quantitative

**I did not find** an auxin-inducible-degron, dTAG, or any other graded-degradation titration
of BubR1 reporting a dose–response against a named readout. Queries in §6.3.

The closest published experiment titrates the *neighbouring* SAC protein and measures
BubR1 as an output. Jema S et al., *Curr Biol* 2023, PMID 37738972 (peer-reviewed;
bioRxiv preprint 2022, doi 10.1101/2022.05.10.491369, under the earlier title "Bub1
availability limits the Spindle Assembly Checkpoint signaling strength of human
kinetochores"):

- SAC protein copy number per kinetochore, measured absolutely: in GSK923295-treated cells,
  Bub1 3.3 ± 1.5, BubR1 3.5 ± 1.3, Mad1 2.7 ± 1.2 per kinetochore; in nocodazole-treated
  cells BubR1 falls to 1.5 ± 0.9 per Knl1.
- **Doxycycline-induced Bub1 overexpression raised kinetochore Bub1 by ~35% and kinetochore
  BubR1 by ~26%**, and lengthened the mitotic delay from ~94 min to ~145 min under 25 nM
  GSK923295 — with mitotic duration correlating with expression level within the
  overexpressing population but not in controls.
- The paper's own summary of the reverse direction: "partial depletion of Bub1 and BubR1
  leads to significantly shorter delays in anaphase onset."

This is the best available answer to the headroom objection in §Appendix C: **a modest
increase in the abundance of a limiting SAC protein produces a proportionate increase in
kinetochore-bound BubR1 and a large increase in checkpoint strength, in human cells, with
numbers.** It is Bub1 that was titrated, and the cells were HeLa and RPE-1 rather than
*BUB1B*-hypomorphic, so it is an analogy, not our experiment. But it is a measured analogy
where §Appendix C currently has none.

One complication in the same lab's work, which belongs next to it: Banerjee A et al., *Mol Biol Cell* 2022, PMID 35767360 shows Bub1-mediated BubR1 recruitment promotes SAC
signalling while highlighting "BubR1's dual role of strengthening the SAC directly and
silencing it indirectly." A BubR1 dose–response is therefore not guaranteed to be monotone,
which is a caveat our §6.2 confirming experiment ("a monotone dose–response") should absorb.

### 4.4 Recent and preprint work PubMed would not yet have shown

**The most important one, four months old.** González-Blanco A et al., *Nat Commun* 2026, PMID 41820377 — "Proteostasis failure and mitochondrial dysfunction contribute to
chromosomal instability-induced microcephaly." A *Drosophila* model of MVA microcephaly
made by depleting SAC genes *bub3* or *rod* **specifically in the neural stem cell
compartment**. Its findings bear directly on our neural arm:

- Microcephaly arises from **loss of NSC stemness** — identity and proliferative capacity —
  producing fewer neurons and glia.
- The driver is the accumulation over time of **complex** aneuploidies (unbalanced gains and
  losses of more than one chromosome), **not** CIN-induced DNA damage and **not** simple
  single-chromosome aneuploidies. NSCs were "highly resistant to DNA damage."
- Mechanism: complex aneuploidy causes proteostasis failure and mitochondrial dysfunction,
  compromising "the highly energy-demanding stemness of NSCs."
- **Rescues:** overexpression of ROS scavengers or mitochondrial chaperones "significantly
  rescued the number of progeny cells and fully restored brain size"; apoptosis inhibition
  (p35) "also resulted in a complete restoration of brain size."

Read against our submission this is both a competitor and an ally. Competitor: it locates
the microcephaly lever in mitochondrial and proteostatic homeostasis downstream of
aneuploidy, not in mitotic fidelity upstream of it, and it demonstrates rescue there. Ally:
NAD⁺ is the central mitochondrial redox cofactor and the substrate of the sirtuins that
govern mitochondrial quality control, so an NAD⁺-raising agent is not orthogonal to the
mechanism this paper identifies. **That second reading is a hypothesis, not a result — this
paper contains no NAD⁺ experiment** — and it should be stated as one. It also independently
confirms that the compartment we worry about, the neural stem cell, is where MVA
microcephaly is generated.

Choi SY et al., *Nucleic Acids Res* 2026, PMID 41533580 (published online 2025-12-27;
bioRxiv preprint 2025-08-22, doi 10.1101/2025.08.22.671725). Aurora B phosphorylates BubR1
at S39 and S16 on kinetochore unattachment, which is a prerequisite for **K250 acetylation**;
K250 acetylation sustains the fibrous corona, enhances CENP-E binding and maintains the MCC.
The authors "identify this axis as a promising therapeutic target in cancers driven by
chromosomal instability." Relevant to us because it makes K250 acetylation a load-bearing,
actively-studied node — see §4.5.

**Bosco N et al., AACR Annual Meeting 2025, Abstract 2983** (Volastra Therapeutics;
**conference abstract, not peer-reviewed**; doi 10.1158/1538-7445.am2025-2983). "BUB1B is an
inactive kinase required for survival of chromosomally unstable cells" — DepMap identifies
BUB1B as a top dependency in genomically unstable lines, and the programme is a BUB1B
*inhibitor* for CIN cancers. Direction of interest: any BUB1B-directed tool compound
reaching the clinic in the next few years will be an inhibitor, and none of it is
repurposable for MVA. Worth one line in the drugs-to-avoid ledger.

**I did not find** any bioRxiv or medRxiv preprint from 2024–2026 proposing an NAD⁺
precursor, a sirtuin activator, or any BubR1-stabilising agent for MVA or for *BUB1B*
hypomorphism. Queries in §6.3.

### 4.5 A sharpening of the K250 objection the document already raises

§Appendix C already flags that a global NAD⁺/sirtuin shift may not be site-selective between
K668 (deacetylation stabilises BubR1) and K250 (acetylation stabilises BubR1), citing Choi
et al., *EMBO J* 2009, PMID 19407811. Two things make that objection sharper than the
document states, and one of them looks like an error in R2.

**(a) SIRT2 is reported to act on K250 as well.** Bloom CR and North BJ, *Cell Biosci* 2021,
PMID 33892776, verbatim: the two sites "regulate BubR1 protein stability in opposite
manners. Acetylation at K250 protects BubR1 from ubiquitination and premature degradation
during mitosis, whereas acetylation at K668 promotes ubiquitin-mediated degradation of BubR1
during interphase." And: "Further studies have indicated that the NAD+-dependent deacetylase
**SIRT2 can also target this site** for deacetylation, though it remains unclear the temporal
and spatial conditions that distinguish K250 regulation by each of the different
deacetylases." So the site-selectivity assumption is not merely untested — the same enzyme is
reported to act at both sites with opposite consequences for BubR1 stability.

**(b) A probable numbering error in R2 that reverses a sign.** The same review states that
the mouse model of K250-acetylation deficiency is **K243R/+**, i.e. **mouse K243 is the
orthologue of human K250** (Park I et al., *J Cell Biol* 2013, PMID 23878276; K243R/+ mice
show high aneuploidy and spontaneous tumours). R2 currently describes Qiu D et al., *Aging Cell* 2018, PMID 29067790 as a Sirt2–BubR1 axis "at **BubR1-K243** rather than K668 …
Different lysine, different system, **same direction**." If mouse K243 is human K250, then
Qiu 2018 is SIRT2 acting at the *protective* site, and its finding that the non-acetylatable
K243R mutant "markedly suppresses spindle/chromosome anomalies and thereupon lowers the
incidence of aneuploid eggs" is not the same direction as the K668 story — it is a case where
*preventing* acetylation at the protective site *improved* fidelity. **This should be checked
against the Qiu paper before the R2 sentence stands.** I did not obtain the Qiu full text in
this pass, so I am flagging a probable error, not asserting one.

The operational consequence is already in §6.1 (measure K250 alongside K668), and it is
correct. What changes is the prior: the risk is higher than "possible non-selectivity", and
the K243R/+ mouse phenotype — high aneuploidy and spontaneous tumours — is precisely the
phenotype an MVA therapy must not induce.

---

## 5. Priority 4 — competitor writeups outside GitHub

**The premise of the earlier survey was wrong, and correcting it matters more than the
count.** Submissions for this competition are not hosted on GitHub. The competition is the
**"Rare Disease, Real Kid: MVA Hackathon 2026"**, run by Hugging Face with Sage
Bionetworks, the MVA Society and BEACON, with $50,000 in prizes from Anthropic and AWS, and
**the submission venue is a Hugging Face Space**:
https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026 (53 likes,
public Community tab). The dataset is gated at
https://huggingface.co/datasets/SageHF2/mva-hackathon-2026-data. Searching GitHub for
competitor hypotheses is searching the wrong index.

**LinkedIn is not empty either.** The launch announcement is a LinkedIn post by Georgia
Channing dated 2026-08-25
(https://www.linkedin.com/posts/georgia-channing_there-is-a-child-with-a-rare-disease-who-activity-7498034092563894273-KYbH),
and an organiser-side post by Robert Allaway (Sage Bionetworks) is dated 2026-08-24
(https://www.linkedin.com/posts/robertallaway_this-project-has-taken-a-while-to-take-shape-activity-7497708238646382592-yHay).
The earlier survey recorded LinkedIn as a zero; it is not.

**The suspicious zero has a boring explanation: timing.** Per the launch post, submissions
opened 25 August 2026 and close 24 October 2026 at 23:59 UTC, with qualitative and expert-panel
judging to 24 November and winners announced 25 November. **Today is 29 August 2026 — the
competition is four days old.** Nobody publishes a retrospective Substack post about a
three-month competition in its first week. The absence of competitor writeups is a
consequence of the calendar, not evidence about the field.

**Live traffic that is visible, treated as observation and not as instruction.** The Space's
public Community tab carried 8 threads when queried on 2026-08-29, including:

| # | Title | Age when queried |
|---|---|---|
| 7 | "Request to reset accidental early Track 2 submission (@josephmayo)" | 1 day |
| 6 | "Anyone who wants to team up?" | 2 days |
| 5 | "Combination with an experimental drug" | 2 days |
| 4 | (title not rendered by the fetch) | 3 days |
| — | "Clarification on using third-party LLM APIs with the dataset" (@doyleywoods) | — |

Two of these are operationally useful to us. Thread #7 states, in the participant's own
words, "I understand that Track 2 permits only one final submission and that re-submissions
are not allowed" — **worth confirming against the official rules before we submit, because
it constrains iteration.** Thread #5 asks the organisers "Can you propose a combination of
drugs? Can you propose drugs that are being used in research?" — i.e. at least one other
entrant is thinking about combinations and investigational agents, which is the space our
catalogue's rejected-because-not-approved entries sit in.

**One inconsistency in the public metadata, flagged for the team to resolve.** The LinkedIn
launch post says submissions open 25 August and close 24 October 2026. The dataset card says
"Challenge period: 16 September 2026 - November 2026". These do not agree, and one of them
governs our deadline.

**What I did not find.** No blog post, Substack, X/Twitter thread, lab page, preprint or
Devpost entry describing an MVA or *BUB1B* drug-repurposing hypothesis for this competition.
No competitor naming NAD⁺, SIRT2, BubR1 stabilisation, nicotinic acid or nicotinamide in
connection with MVA. Adjacent-but-different events surfaced instead — HackRare 2026 at
Harvard (43 rare-disease projects across genomic diagnostics, symptom management and
therapeutic targets, none MVA), and generic AI-drug-repurposing writeups (Every Cure
coverage; a Codex-driven IL-17A repurposing Substack; assorted knowledge-graph repurposing
repos). Queries in §6.4. **This is "I did not find", not "there is none" — and given that
the competition is four days old and its submissions live behind a Hugging Face Space rather
than in a public repo index, the correct expectation is that most competitor material is not
yet public and some of it never will be.**

---

## 6. Every query run, including the ones that found nothing

Engine: Exa `web_search_exa` / `web_fetch_exa` unless marked otherwise. Database endpoints
were queried directly where a paper did not exist to cite; PMID resolution used the Europe
PMC REST API, which is citation lookup rather than discovery search.

### 6.1 NAPRT in neural progenitors (Priority 1)

| Query | Outcome |
|---|---|
| `NAPRT nicotinate phosphoribosyltransferase expression in human neural stem cells and neural progenitor cells NAD salvage pathway during neural differentiation` | Productive — Duarte-Pereira 2016, Cole 2017, Neves 2022, Wang 2016, HPA brain page |
| `single cell RNA sequencing human fetal brain atlas NAPRT expression radial glia intermediate progenitor cells cerebral organoid neuroepithelium` | **No paper reporting NAPRT in progenitors.** Returned atlas papers with no NAPRT content; this is what sent me to the atlases themselves |
| `Duarte-Pereira Fajarda Matos Oliveira Silva "NAPRT Expression Regulation Mechanisms…" Genes 2021 pubmed PMID` | Productive — PMID 34946971 |
| `Neves Goodfellow Vieira Silva "The role of NAD metabolism in neuronal differentiation" Neurochemistry International 2022 SH-SY5Y NAPRT QPRT PMID` | Productive — PMID 35843422 |
| `Braun 2023 Science comprehensive cell atlas of the first trimester developing human brain Linnarsson PMID pubmed` | Productive — PMID 37824650 |
| `Nowakowski Kriegstein 2017 Science spatiotemporal gene expression trajectories reveal developmental hierarchies of the human cortex PMID pubmed 29217575` | Productive — PMID 29217575 confirmed |
| web_fetch `proteinatlas.org/ENSG00000147813-NAPRT/single+cell` and `/brain` | **Failed** — both returned only the site's search-widget boilerplate through this fetcher |
| direct API `proteinatlas.org/api/search_download.php?search=NAPRT` | Partial — returns only summary categories ("Cell type enhanced", "Low region specificity"); no per-cell-type values |
| direct `proteinatlas.org/ENSG00000147813-NAPRT.json` | **Failed** — returns HTML, not JSON |
| direct UCSC Cell Browser `cells.ucsc.edu/cortex-dev/*` (dataset.json, meta.tsv, exprMatrix.tsv.gz) | **Productive — this is Atlas B in §2.1** |
| direct UCSC Cell Browser `cells.ucsc.edu/dev-brain-regions/wholebrain/*` (dataset.json, exprMatrix.json index, byte-range reads of exprMatrix.bin, meta.tsv) | **Productive — this is Atlas A in §2.1** |
| direct UCSC Cell Browser `cells.ucsc.edu/organoidatlas/*` (Tanaka-style organoid + fetal brain, 190,022 cells, has a "Neuroepithelial cell" annotation) | **Dead end** — its random-access gene index contains only 2,001 genes and **NAPRT is not among them**; the full matrix is 3.06 GB and was not downloaded |

**Not attempted, and worth a later pass:** Allen Brain Map / BrainSpan developmental bulk
RNA-seq; the Descartes human fetal gene atlas (Cao 2020); CZ CELLxGENE Discover's
cross-dataset gene-expression API; HDBR.

### 6.2 NAD⁺ and neural progenitors (Priority 2)

| Query | Outcome |
|---|---|
| `Huang Cheng NAD+ modulates the proliferation and differentiation of adult neural stem/progenitor cells via Akt signaling pathway Cells 2022 apoptosis` | Productive — PMID 35455963, full Results text, GSE200494 |
| `extracellular NAD+ millimolar cytotoxicity NAD-induced cell death ART2 ADP-ribosyltransferase P2X7 mouse but ART2 is a pseudogene in humans species difference` | Productive — Glowacki 2002, Haag 2007, Hong 2009 |
| `Glowacki Koch-Nolte 2002 Protein Science human ART2 ARTC2 pseudogene humans lack functional ART2 gene PMID` | Productive — PMID 12070318 |
| `P2X7 receptor expressed in neural progenitor cells neurospheres induces apoptosis and inhibits proliferation extracellular nucleotides` | Productive — Leeson 2018, Delarasse 2009, Messemer 2013, Lovelace 2015 |
| `Leeson 2018 Stem Cells P2X7 receptors regulate phagocytosis and proliferation in adult hippocampal and SVZ neural progenitor cells PMID pubmed` | Productive — PMID 30068016 |
| `Stein Imai 2014 EMBO Journal specific ablation of Nampt in adult neural stem cells recapitulates their functional defects during aging PMID` | Productive — PMID 24811750 |
| `nicotinamide mononucleotide or nicotinamide riboside increases neural stem progenitor cell proliferation neurogenesis in vivo mouse study` | Productive — Zhao 2015 (PMID 26060246), Wang 2016 (PMID 27018006) |
| `NAD deficiency congenital malformations and niacin supplementation Shi Dunwoodie New England Journal of Medicine 2017 vertebral heart limb defects mice niacin prevented` | Productive — PMID 28792876 |

**Not found in this pass:** any study giving a NAD⁺ **precursor** to a **human** neural
progenitor and measuring proliferation; any study testing P2X7 or ARTC2 involvement in
Huang 2022's phenotype; any NAD⁺ or precursor experiment in a cerebral organoid.

### 6.3 The Priority 3 gaps

| Query | Outcome |
|---|---|
| `nicotinamide or NAD+ or sirtuin or SIRT2 treatment tested in mosaic variegated aneuploidy patient cells or BubR1 hypomorphic mice therapeutic intervention raising BubR1 levels` | **No intervention study in MVA.** Surfaced Pun 2024 (PMID 39490401) and North 2014's MVA paragraph — i.e. the *connection* exists in the literature, the *experiment* does not |
| `auxin-inducible degron titration of BubR1 BUB1B protein levels dose response spindle assembly checkpoint strength chromosome missegregation threshold` | **No BubR1 degron/titration study found.** Nearest: Jema 2023 (PMID 37738972), Banerjee 2022 (PMID 35767360), Baker 2007 graded *Bub1* allelic series |
| `Chen Joglekar signaling protein abundance modulates the strength of the spindle assembly checkpoint Bub1 BubR1 recruitment published journal 2023 2024` | Productive — PMID 37738972 and its 2022 bioRxiv preprint |
| `preprint bioRxiv 2024 2025 BUB1B BubR1 mosaic variegated aneuploidy therapeutic strategy small molecule restore spindle assembly checkpoint` | **No MVA-therapeutic preprint.** Surfaced Malumbres 2024 review, AACR 2025 abstract 2983, Suijkerbuijk 2010 |
| `Aurora B kinase phosphorylation required for BubR1 acetylation at lysine 250 mitotic checkpoint complex 2025 2026 paper` | Productive — Choi 2026 (PMID 41533580) plus its bioRxiv preprint; Bloom and North 2021 (PMID 33892776) |
| `SIRT2 deacetylates BubR1 lysine 250 K250 acetylation PCAF prometaphase APC/C activity mitosis timing Choi Lee EMBO J 2009 Park 2013 K243R mice tumour` | Productive — Choi 2009 (PMID 19407811), Park 2013 (PMID 23878276) |
| web_fetch Malumbres and Villarroya-Beltri 2024 review PDF (mvasociety.org mirror) | Productive — confirmed no NAD⁺/sirtuin/BubR1-stabilisation proposal, and the "<5% of normal levels" sustainability claim |
| web_fetch `pmc.ncbi.nlm.nih.gov/articles/PMC11625622/` | **Failed** — reCAPTCHA interstitial. Obtained via `ncbi.nlm.nih.gov/pmc/articles/PMC11625622/` instead |
| web_fetch `europepmc.org/article/MED/39490401`, `embopress.org/doi/pdf/10.15252/embj.201386907` | Productive |
| web_fetch `pubmed.ncbi.nlm.nih.gov/?term=<doi>` | **Failed** — "Cookies must be enabled". PMIDs resolved via Europe PMC REST instead |
| direct Europe PMC REST, DOI/title lookups for 10 citations | Productive; every PMID in this document was resolved this way or read off a PubMed/PMC page |
| direct `eutils.ncbi.nlm.nih.gov` esearch by DOI | **Failed** — empty response from this environment |

**Still not found, stated as gaps:** a BubR1 degron dose–response; any quantified
residual-BubR1 measurement in MVA patient cells published after 2010 (Suijkerbuijk 2010 and
Rio Frio 2010 remain the only two); any NAD-pathway intervention in an MVA model; any
preprint proposing a BubR1-stabilising therapy for MVA.

### 6.4 Competitor writeups (Priority 4)

| Query | Outcome |
|---|---|
| `blog post or Substack describing a drug repurposing hypothesis for mosaic variegated aneuploidy BUB1B for an AI biology hackathon competition track 2 submission` | **Nothing on point.** Returned generic AI-repurposing content (Codex/IL-17A Substack; Every Cure; HAMGNN paper) |
| `"mosaic variegated aneuploidy" hackathon submission hypothesis writeup rare disease drug repurposing AI agent competition 2026 personal project writeup` | **Found the venue, not a competitor** — the Hugging Face Space; plus HackRare 2026, Sugi Atlas, and the Nat Commun 2026 microcephaly paper |
| `Lyzr bio-hackathon Track 2 rare disease hypothesis submissions gallery participant writeups BUB1B mosaic variegated aneuploidy 2026` | **Found the venue and the LinkedIn launch**; no competitor hypothesis |
| `"Rare Disease, Real Kid" MVA hackathon 2026 submission BUB1B hypothesis Track 2 entry` | Productive on logistics — LinkedIn launch post, dataset card, Space discussions |
| `personal writeup of my approach to the Hugging Face MVA rare disease hackathon track 2 drug repurposing hypothesis blog post` | **Nothing on point.** Returned unrelated repurposing repos (open-cure, DrugSight, Bisociation-engine), Devpost entries from other hackathons, a DEV.to knowledge-graph post |
| web_fetch the Space, its `/discussions`, threads `/5` and `/7`, and the dataset card | Productive — thread list, the one-final-submission rule as stated by a participant, the date inconsistency |

**Not attempted:** X/Twitter search (Exa surfaced X content only indirectly, via a Digg
aggregation of the launch post); the Space's Files tab; any gated or logged-in view of
submissions.

---

## 7. What this changes in the submission, in priority order

1. **§4.5 "Neural progenitors: probably not" must be rewritten.** The claim rests on an
   unpublished observation about mature neurons in a tumour paper, and two human single-cell
   atlases put *NAPRT* in radial glia and IPCs above the noise floor and at least as high as
   in neurons. Replace with: transcript present and low, protein and activity unmeasured,
   NAMPT roughly 2× higher, gate stays but is now expected to pass.
2. **§4.5 "no prior work connects this disease to this pathway" must be withdrawn.** North
   2014 names MVA; Pun 2024 proposes the axis therapeutically. Replace with the narrower,
   still-true claim about no intervention study, and note that neither names nicotinic acid.
3. **§7 needs a new limitation: the NAD⁺/progenitor sign conflict** (§3.5 above), with the
   exogenous-versus-precursor split and the mouse-P2X7 caveat stated, and a
   progenitor-proliferation readout added to the neural arm of §6.1.
4. **R2's Qiu 2018 sentence is probably wrong** — mouse K243 appears to be human K250, the
   protective site, which reverses "same direction". Verify against the Qiu full text.
5. **Add Rio Frio 2010 at 10–15% residual** to R3's allelic table, with the OMIM-summary
   provenance flagged, and note that it points the microcephaly threshold *below* the cancer
   threshold.
6. **Add Jema 2023 to Appendix C** as the measured answer to the headroom objection, and
   Banerjee 2022 as the reason §6.2 should not assume monotonicity.
7. **Add González-Blanco 2026** as the competing mechanism for the neural arm, and as an
   argument that a mitochondrial-cofactor intervention is not orthogonal to it.
8. **Operational:** confirm the real Track 2 deadline (25 Aug–24 Oct per LinkedIn vs
   16 Sep–November per the dataset card) and the one-final-submission rule before iterating.
