# Track 2 — A drug-repurposing hypothesis for BUB1B-related mosaic variegated aneuploidy

**Mosaic variegated aneuploidy syndrome 1** · OMIM 257300 · MONDO:0009759 ·
Orphanet:1052 · *BUB1B* (HGNC:1149, ENSG00000156970, UniProt O60566, 15q15.1)

---

## 0. Scope, and what the citations guarantee

A research hypothesis about a gene, assembled from public literature and run through this
repository's mechanism and intervention stages. Every factual claim carries an author,
year, journal and PMID; claims with no publication behind them are labelled **inference**.

**An earlier version of this paragraph claimed that every citation had been resolved
against PubMed before being written down. That claim was itself insufficiently checked,
and it is the first thing this document has to correct.** Resolving a citation establishes
that the paper exists and is the one named. It does not establish that the paper says what
the sentence attached to it says, and those are different questions with different failure
modes.

Every claim has since been re-checked against the source rather than the identifier.
`docs/track2-claims-ledger.md` records the result: 77 claims, of which 47 hold, 22 are
weakened, 3 are withdrawn and 5 could not be verified. Appendix I gives the four most
instructive catches. So the guarantee is two-tier and only the first tier is mechanical: a
test enforces that every row resolves to a cited publication, and reading enforces that
every claim matches what that publication says, which failed 25 times out of 77.

This is not medical advice, and it is not derived from patient data: no genotype,
coordinate or clinical record was read, requested or inferred (ADR 0008, GP-44). Because
the patient workspace is deliberately illegible to the agent that wrote this, the
hypothesis is stated genotype-conditionally (§7.3).

---

## 1. The answer

*BUB1B* pairs a null or truncating allele with a hypomorphic one. The deficient quantity
is BubR1 protein. The disease is named after the spindle assembly checkpoint, but the
checkpoint fails because the protein is scarce.

So the therapeutic target is protein abundance. The route with the most evidence behind it
is NAD⁺ → SIRT2 → deacetylation of BubR1 at K668, which slows the protein's degradation.

Abundance is a floor rather than a full explanation. Sieben et al., *J Clin Invest* 2020,
PMID 31738183 compared mouse genotypes with indistinguishable total BubR1 and
near-identical aneuploidy, and found median lifespans differing about 1.75-fold. Total
protein predicts whether cells missegregate. It does not predict how severe the disease is
(§4.2).

Direction of effect splits at two levels, and the two disagree.

| Level | Verdict | Basis |
|---|---|---|
| Node: does raising BubR1 correct the disease? | Agrees | Ectopic BubR1 restored checkpoint activity in MVA patient cells (Suijkerbuijk et al., *Cancer Res* 2010, PMID 20516114). The readout was a 2-fold rise in nocodazole response against the carrier parent; neither aneuploidy nor missegregation was measured. Sustained elevated BubR1 in mice reduced aneuploidy and tumorigenesis (Baker et al., *Nat Cell Biol* 2013, PMID 23242215), in animals that were not BubR1-deficient. |
| Molecule: does any approved drug do it? | Cannot determine, for all 13 | The one agent with a measured, agreeing direction is NMN (North et al., *EMBO J* 2014, PMID 24825348), which is not an approved medicine. Nicotinamide raises NAD⁺ and is also a measured competitive inhibitor of human SIRT2 (Marcotte et al., *Anal Biochem* 2004, PMID 15301953), so it feeds and blocks the same axis. Nobody has measured the net effect of doing both at once. |

Of 13 agents evaluated, 1 was accepted and 12 rejected. The accepted agent carries an
undetermined direction, and the highest-scoring agent in the catalogue is one of the
rejected ones.

The accepted agent is **nicotinic acid**, and it is the mechanistic probe for one reason
above the others: it moves a single arm of the mechanism, so a result from it is readable.
Nicotinamide moves both arms at once — it raises NAD⁺ and inhibits the sirtuin that spends
it — so neither a positive nor a negative result from nicotinamide would identify the
pathway. It stays in the experiment as a deliberate **high-dose negative control**, where
the mechanism predicts BubR1 should *fall* (Appendix G.3).

That ordering is the reverse of what this document printed at first, and the correction is
the most useful thing in it. Nicotinamide was ranked first and nicotinic acid rejected, on
two catalogue cells that were both wrong: one filled from a skin-cancer trial that measures
no chromosome number, the other left blank although the measurements existed. §5 shows the
cells, the sources and the ranking that moved once they were fixed.

The deliverable is an experiment (§6) and a statement of what would falsify the hypothesis.

---

## 2. The mechanism

The chain is machine-readable: `knowledge/literature/bub1b/mechanisms.tsv`, 12 nodes and 13
links, each carrying a tier, a strength, an `is_directly_demonstrated` flag and a stated
uncertainty, because a chain presented at uniform confidence is a rhetorical device
(ASSUMPTION-MECHANISM-01). Appendix A gives every link with its evidence.

```
BUB1B alleles → transcript ↓ → BubR1 protein ↓ ─┬→ MCC ↓ → SAC ↓ → PCS ↑ ─┐
                                  [TARGET]      │                          ├→ aneuploidy ↑
                                                └→ KT–MT attachment ↓ ─────┘
                                                                            │
   ┌────────────────────────────────────────────────────────────────────────┤
   ↓                                    ↓                                   ↓
proteotoxic/metabolic stress ↑    p53/p16 arrest, senescence,        cancer risk ↑
                    └────(inferred)───▶ immune clearance ↑           (rhabdomyosarcoma,
                                        [COMPENSATORY]                Wilms, leukaemia)
                                              │
                                        (inferred, weakest)
                                              ↓
                              growth restriction, microcephaly
```

Eleven links are directly demonstrated; two are not. L11 holds that aneuploidy-associated
stress *causes* the arrest and senescence response. L12 holds that progenitor attrition
produces growth restriction and microcephaly, against a published alternative through
ciliogenesis (Miyamoto et al., *Hum Mol Genet* 2011, PMID 21389084).

One node is compensatory and cannot be a target. `N10_clearance` — p53/p16 arrest,
senescence and immune clearance of missegregating cells — deviates from wild type, and
pushing it back would suppress a response holding a cancer-predisposed child's aneuploid
cells in check. The same response drives the progeroid phenotype in BubR1-hypomorphic mice
(Baker et al., *Nat Cell Biol* 2008, PMID 18516091) and suppresses aneuploidy-driven
tumorigenesis (Li et al., *PNAS* 2010, PMID 20663956), so the direction check returns
`UNKNOWN` and the model refuses the node as a target (ASSUMPTION-MECHANISM-04).

The target is the protein rather than the checkpoint because BubR1 reaches missegregation
twice: patient cells show an impaired checkpoint *and* alignment defects (Suijkerbuijk et
al., *Cancer Res* 2010, PMID 20516114), and alignment runs through the KARD–PP2A-B56 arm,
independent of MCC assembly (Suijkerbuijk et al., *Dev Cell* 2012, PMID 23079597).
Abundance is the last point at which one quantity governs both (ADR 0025), and activity is
no alternative, since human BUBR1's retained catalytic triad is dispensable for error-free
segregation (Suijkerbuijk et al., *Dev Cell* 2012, PMID 22698286).

---

## 3. Direction of effect

The naive move is to strengthen the checkpoint. No approved agent strengthens the spindle
assembly checkpoint; the class does not exist.

The naive counter-move is the oncology one: inhibit MPS1, weaken the checkpoint further,
drive lethal aneuploidy. Reversine does exactly that (Santaguida et al., *J Cell Biol*
2010, PMID 20624901), and the class is in active anticancer development (Wang et al., *Eur
J Med Chem* 2019, PMID 31121430). It inverts here because the affected soma is not a
tumour: oncology pushes chromosomally unstable cells past an aneuploidy-tolerance ceiling,
and in this disorder a substantial and variable fraction of non-tumour cells already sits
at that ceiling in every tissue. Target proximity is maximal and the sign is inverted;
reversine scores 0.110, last of 13.

The same inversion catches three approved drugs, which is the harder case. Tang et al.,
*Cell* 2011, PMID 21315436 found AICAR, 17-AAG and chloroquine selectively lethal to
aneuploid cells, because aneuploid cells' autophagic capacity is already limited
(Santaguida et al., *Genes Dev* 2015, PMID 26404941) and their HSP90-dependent folding
already impaired (Donnelly et al., *EMBO J* 2014, PMID 25205676). Against an MVA-associated
embryonal tumour that is attractive; in the soma it spares no compartment.
Hydroxychloroquine, metformin and tanespimycin are cheap, well tolerated and paediatrically
exposed, which is what makes them dangerous to a search that scores availability and cannot
see a sign, and all three are rejected `WRONG_DIRECTION` on the evidence that recommends
their analogues for the tumour arm.

Two consequences. This pipeline scores one context — the soma — so every rejection means
"not for systemic use in the affected child", never "useless in this disease"
(ASSUMPTION-DRUG-08), and what decides the tumour arm is a therapeutic index rather than a
mechanism chain. And the primary hypothesis escapes the tension, because sustained elevated
BubR1 reduced aneuploidy *and* tumorigenesis in the same animals, even against oncogenic Ras
(Baker 2013). Appendix B lists all ten classes considered.

---

## 4. The primary hypothesis

> In *BUB1B*-MVA alleles of the abundance-lowering class, pharmacological elevation of the
> NAD⁺ pool raises BubR1 protein abundance through SIRT2-dependent deacetylation of BubR1
> lysine-668, and — if it raises residual BubR1 far enough to move the dose–response, whose
> position nobody has measured in a disease-relevant cell type and which differs by readout
> — reduces the ongoing rate of chromosome missegregation.

### 4.1 There is no single threshold, because there is no single readout

"The threshold" is the wrong noun. Residual BubR1 and mitotic failure are related by a
graded dose–response whose steepness depends on which function is scored, in which cell
type, over how many divisions. Bohers et al., *Hum Genet* 2008, PMID 18932004 scored
premature chromatid separation and accumulated aneuploidy in depleted HeLa clones;
Suijkerbuijk 2010 scored per-division missegregation in the same cell type; neither locates
a threshold, and because they never measured the same quantity their apparent fourfold
disagreement is not evidence about where one sits. The number that bears on the therapy is
in neither series: MVA patient cells sit at roughly 17–50% residual full-length BubR1
(Suijkerbuijk 2010), inside the interval nobody titrated. Appendix A.3 ranks the readouts.

### 4.2 Abundance is the floor, not the severity predictor

Sieben 2020 compared `BubR1^H/L1002P` and `BubR1^H/H` mice, whose BUBR1 protein and
tissue-FISH aneuploidy rates are near-identical and whose median lifespans differ 1.75-fold,
343 days against 196; their conclusion is that "BUBR1 allelic effects beyond protein level
and aneuploidy contribute to disease heterogeneity". Level still sets the regime in the same
data: a 40–90% BUBR1 reduction gives tumour predisposition or mild premature ageing, more
than 90% gives progeroid phenotypes or early postnatal lethality.

That narrows the thesis without refuting it. Raising total BubR1 remains the only
intervention that has restored checkpoint function to the L1012P protein in human cells
(Suijkerbuijk 2010) and the only one that has corrected progeria in an MVA mouse model
(Baker 2013), so it stays the primary pharmacodynamic variable — as a target-engagement
readout rather than an outcome measure, and raising it past a regime boundary should move
the phenotype toward the mild end of a regime rather than to normal. Appendix A.4 has the
comparison, including the one reading of Sieben 2020 its own figures contradict.

### 4.3 Why the NAD⁺/SIRT2 axis, and where it runs out

North 2014 showed that BubR1 abundance is set by the acetylation state of lysine-668, which
SIRT2 holds deacetylated in an NAD⁺-dependent manner. The same paper hit a ceiling, and the
ceiling is the most important result here: SIRT2 overexpression raised median lifespan 123%
in male `BubR1^H/H` mice, but could not restore BubR1 to wild-type levels, because the
hypomorph's deficit is transcriptional and deacetylation "is only able to stabilize the
protein that is present". That is the lesion *BUB1B*-MVA has, so "enough residual BubR1
exists for a stabilisation-only mechanism to push it over a functional threshold" is a
separate, load-bearing premise — and the one experiment bearing on it came out against it.
Appendix G.1 has the doses and the tissue.

### 4.4 The sign paradox at the molecule, and which molecule the experiment uses

Four objections are answered in Appendix C; the strongest is that the NAD⁺ evidence comes
from an ageing context, so a young child has not lost NAD⁺ and the pharmacology may have
nothing to push against. I do not have an answer to it.

The molecule is where the direction stops being signable. The agent with the best paediatric
safety evidence on this axis is nicotinamide — ENDIT randomised 552 first-degree relatives,
children included, to 1.2 g/m²/day or placebo for five years and reported no effect on
growth in children (Gale et al., *Lancet* 2004, PMID 15043959) — and it is the one whose
direction cannot be signed, because it inhibits the sirtuin this mechanism runs on. That
inhibition is measured against human SIRT2 rather than inferred from the family's conserved
pocket (Marcotte 2004). Nicotinamide replenishes the NAD⁺ pool SIRT2 consumes and supplies
the product that inhibits SIRT2, in the same cell at the same time, and which term dominates
at a given intracellular concentration has never been measured in any genotype:
`context_dependent` at the node, `cannot determine` at the molecule.

**The probe is nicotinic acid.** Not because nicotinamide would fail — because nicotinamide
cannot produce an interpretable result, and interpretability is the whole job of a first
experiment. Nicotinamide moves both arms of the proposed mechanism in opposite directions
at once. If BubR1 rises, NAD⁺ supply won or something off-axis did, and the assay cannot
tell you which. If BubR1 falls or does not move, the mechanism is wrong or the enzyme it
runs on was inhibited, and the assay cannot tell you which either. Both outcomes leave the
hypothesis where it started. Nicotinic acid enters through the Preiss-Handler route, is not
a sirtuin reaction product, and moves one arm, so its result is readable in both
directions.

**Nicotinamide stays in the experiment as the high-dose negative control**, and it is a
good one, because the mechanism makes a signed prediction about it: above the measured
SIRT2 IC₅₀ range, nicotinamide should *lower* BubR1. That prediction has been met once
already, in the founding paper for this axis — 20 mM nicotinamide for 16 h reduced BubR1
protein in HeLa cells, with no change in cell-cycle distribution and none in *BUB1B* mRNA
(North et al., *EMBO J* 2014, PMID 24825348, Figure 1A and 1B). 20 mM is 14–25× the peak
plasma a patient reaches, so it fixes the sign of the inhibitory arm and not the dose at
which that arm wins. Appendix G.3 has the exposures, the IC₅₀ spread and the four reasons
the comparison is not settled; Appendix G.2 gives all four NAD⁺-axis agents.

### 4.5 Does niacin reach the cell? NAPRT, and a claim that turned out to be wrong

Nicotinic acid reaches NAD⁺ only through the Preiss-Handler pathway, whose entry enzyme is
NAPRT. If the cells this disease damages do not express NAPRT, niacin does nothing in them
and a null result would be uninterpretable — which would not rescue nicotinamide, but would
change what the first experiment has to confirm before it confirms anything else. So the
tissue distribution of NAPRT is load-bearing, and the version of it this project had been
carrying is wrong.

**Where the "gut, liver and kidney only" claim comes from, and what it actually says.** It
traces to Hara et al., *J Biol Chem* 2007, PMID 17604275, whose sentence is: NAPRT "was
abundantly expressed in the small intestine, liver, and kidney in mice". That is a
statement about **where expression is highest, in mice**. It is not a statement that human
expression is exclusive to those organs, and it has been restated as one — here included.
The same paper is the origin of the two facts that make niacin attractive: in cells
expressing endogenous NAPRT, adding nicotinic acid but *not* nicotinamide almost doubled
cellular NAD⁺ and the effect was reversed by NAPRT knockdown; and NAPRT, unlike NAMPT, is
insensitive to the physiological concentration of NAD⁺, so the Preiss-Handler entry step is
not feedback-inhibited.

**What the human data say.** NAPRT is broadly expressed in normal human tissue. An RT-PCR
panel across normal human tissues found the gene "widely expressed under normal conditions"
(Duarte-Pereira et al., *Oncotarget* 2016, PMID 26675378). At protein level, the only
cell-type-resolved survey done with a validated antibody rather than the poor commercial
ones states that "most normal tissues express NAPRT while a significant proportion of
malignant cells do not" (Cole et al., *Oncotarget* 2017, PMID 29100430). The reason NAPRT
looks restricted in the literature is that most of that literature is oncology: loss of
NAPRT is a **tumour** phenomenon, caused by tumour-specific promoter methylation and found
in 5–65% of lines depending on cancer type (Shames et al., *Clin Cancer Res* 2013, PMID
24097869). Reading a tumour-cell-line frequency as a normal-tissue distribution is the
error.

**Lymphocytes: positive, and this is the strongest evidence on the list.** Primary human
peripheral blood mononuclear cells raised total cellular NAD⁺ in a concentration-dependent
way on nicotinic acid, significantly above 5 µM and 2.1-fold on average at 15 µM, in 13 of
13 donors (Weidele et al., *Biochem Pharmacol* 2010, PMID 20599792). The micronucleus
result this document leans on in §5.2 is from the same cell type and the same 15 µM
exposure (Weidele et al., *DNA Repair (Amst)* 2017, PMID 28216063), and that paper states
that PBMCs express nicotinate phosphoribosyltransferase, which is why the model was chosen.
NAPRT enzyme activity had already been measured directly in normal human lymphocytes
(Sestini et al., *Arch Biochem Biophys* 2000, PMID 10898945). For scale: 15 µM is about a
fifth of the ~75 µM peak plasma nicotinic acid after 2 g of extended-release niacin (Menon
et al., *Int J Clin Pharmacol Ther* 2007, PMID 17725178), so the effective ex vivo exposure
is comfortably reachable in a patient.

**Fibroblasts: positive on transcript and on function, but not directly measured.** GTEx v8
puts *NAPRT* at a median 23.5 TPM in cultured primary fibroblasts and 32.1 TPM in
EBV-transformed lymphocytes, against a floor of 9.1 TPM in the lowest of all 54 tissues and
a whole-atlas maximum of 106.7 — skin, at 94.8 TPM, outranks liver at 92.9 (GTEx Portal, v8
median gene expression, *NAPRT* ENSG00000147813; queried 2026-08-29). Functionally,
nicotinic acid rescued FK866-induced NAD⁺ depletion and senescence in Hs68 human
fibroblasts (Song et al., *Biogerontology* 2015, PMID 26330291), and nicotinic acid rescue
of a NAMPT inhibitor is precisely the assay that reports NAPRT competence, because it works
only in NAPRT-expressing cells (Watson et al., *Mol Cell Biol* 2009, PMID 19703994).

Three caveats travel with that, and the first experiment has to absorb them. No published
measurement of NAPRT protein or enzyme activity in primary human **dermal** fibroblasts was
found, and no study measuring the NAD⁺ response to nicotinic acid in them; searches are
listed below. Immortalised fibroblast lines are frequently NAPRT-low, which is consistent
with methylation on culture adaptation, so a result in an hTERT line does not transfer to a
patient line or the reverse. And antibody quality in this field is poor enough that
protein-level negatives from the older commercial reagents should not be trusted. **The
operational consequence is a gate, not a caveat: confirm NAPRT in the specific patient
lines before running the arm, by NAPRT immunoblot or by the nicotinic-acid-rescue assay.**
The catalogue's `validation_experiment` cell for nicotinic acid now says so.

**Neural progenitors: probably not, and the neural arm cannot use niacin.** With the
validated antibody, neurons in non-tumorigenic human cortex and white-matter
oligodendrocytes were NAPRT-negative, while reactive astrocytes and endothelium were
positive (Cole et al., *Oncotarget* 2017, PMID 29100430); more than 70% of gliomas were
negative in the same study, and glioblastomas and neuroblastomas are among the tumour types
specifically noted as NAPRT-deficient and not rescuable by nicotinic acid (Watson et al.,
*Mol Cell Biol* 2009, PMID 19703994). Brain is the bottom of the GTEx atlas for this gene.
This matters because microcephaly is the phenotype the therapeutic argument is weakest on
already (§7.2): if the neural compartment is NAPRT-negative, niacin cannot address it even
if the mechanism is right everywhere else, and a neural arm would need a different
precursor.

**The searches that came back empty**, stated as gaps rather than as findings. PubMed
`NAPRT AND ("neural progenitor" OR "neural stem" OR iPSC OR "neuroepithelial")` returns **0
hits** — there is no dedicated study of NAPRT in human neural progenitors, and the
paragraph above is inference from adult neuronal protein staining, neural tumour lines and
the transcript atlas. PubMed `NAPRT AND fibroblast` returns **1 hit**, a murine colitis
paper, so the fibroblast conclusion rests on the transcript atlas plus one indirect rescue
experiment. And PubMed `"mosaic variegated aneuploidy" AND (NAD OR nicotinamide OR niacin)`
returns **0 hits**: no prior work connects this disease to this pathway, which is a reason
to state the hypothesis carefully and not a reason to state it confidently.

---

## 5. What the pipeline produced

Reproduce with `uv run pytest -q tests/unit/test_track2_bub1b.py` (26 tests). Tables:
`knowledge/literature/bub1b/`. Output: 1 accepted, 12 rejected, 89 evidence rows, every
rejection kept with its reasons (GP-19).

| Rank | Agent | Node | Required | Observed | Direction verdict | Score | Outcome |
|---|---|---|---|---|---|---|---|
| **1** | Nicotinic acid | `N3_bubr1` | increase | unknown | **cannot determine** | 0.500 | accepted, 2 concerns |
| — | Nicotinamide mononucleotide | `N3_bubr1` | increase | increase | **AGREES** | 0.605 | rejected — `oncogenic_risk` |
| — | Gentamicin | `N3_bubr1` | increase | unknown | cannot determine | 0.420 | rejected — `oncogenic_risk` |
| — | Dasatinib | `N10_clearance` | *unknown (compensatory)* | decrease | cannot determine | 0.348 | rejected — `oncogenic_risk` |
| — | Nicotinamide riboside | `N3_bubr1` | increase | unknown | cannot determine | 0.345 | rejected — `oncogenic_risk` |
| — | Sirolimus | `N9_stress` | decrease | unknown | cannot determine | 0.303 | rejected — `oncogenic_risk` |
| — | Nicotinamide | `N3_bubr1` | increase | context_dependent | cannot determine | 0.275 | rejected — `oncogenic_risk` |
| — | Ataluren | `N3_bubr1` | increase | unknown | cannot determine | 0.203 | rejected — `oncogenic_risk` |
| — | Hydroxychloroquine | `N9_stress` | decrease | increase | **OPPOSES** | 0.163 | rejected — `wrong_direction` |
| — | Metformin | `N9_stress` | decrease | increase | **OPPOSES** | 0.163 | rejected — `wrong_direction` |
| — | Vorinostat | `N5_sac` | restore | decrease | **OPPOSES** | 0.133 | rejected — `wrong_direction` |
| — | Tanespimycin | `N9_stress` | decrease | increase | **OPPOSES** | 0.130 | rejected — `wrong_direction` |
| — | Reversine (MPS1i) | `N5_sac` | restore | decrease | **OPPOSES** | 0.110 | rejected — `wrong_direction` |

**This is the second ranking this document has printed, and the first one was wrong.** An
earlier run put nicotinamide at rank 1 with 0.365 and rejected nicotinic acid at 0.410 —
the rejected agent outscoring the accepted one — on a single cell: `worsens_cin`. Both
cells were wrong, in opposite directions. Nicotinamide's had been filled from a
skin-cancer chemoprevention trial that measures no chromosome number (§5.1). Nicotinic
acid's had been left blank although two randomised malignancy tabulations and a human
micronucleus assay already existed (§5.2). Correcting the two cells was the whole
intervention: no weight, no gate, no rule and no line of `src/mva/interventions/` was
touched, and the ranking inverted on its own. The gate was right and the catalogue was
wrong, which is the failure mode this pipeline is least able to catch by itself — a
literature-curated table is exactly as good as the reading behind it.

The highest-scoring agent is still rejected. Approval status and the unassessed
oncogenic-risk question are gates rather than score terms, so NMN outscores everything and
does not survive; a pipeline where the top score always wins cannot express "right about
the biology, not yet a candidate". The accepted candidate scores 0.500 out of 1 and still
carries an undetermined direction: nobody has measured what nicotinic acid does to BubR1
abundance, and `unknown` is not agreement (ASSUMPTION-DRUG-02).

Seven of the twelve rejections are `oncogenic_risk`, and in all seven the risk was never
measured rather than found. In a chromosomal-instability syndrome "nobody has asked whether
this agent increases aneuploidy" is disqualifying (ADR 0011), and it disqualifies almost
everything — the most transferable finding here, because the bottleneck is an unasked safety
question rather than target identification.

### 5.1 The gate that should not have been open

Nicotinamide's `worsens_cin` cell had been recorded as assessed-and-negative on Chen et
al., *N Engl J Med* 2015, PMID 26488693 — a phase 3 trial, 386 adults, 12 months, which
found 23% fewer new non-melanoma skin cancers. It is a real cancer endpoint and it is more
than most agents in this catalogue have. It is also the wrong evidence for this gate, on
four counts, and none of them is a technicality.

It measures the wrong tumour biology: UV-mutagen-driven keratinocyte carcinoma in a field
of sun-damaged skin, through enhanced repair of nucleotide-level lesions. MVA1 tumours are
embryonal, on a germline background of whole-chromosome missegregation. Nothing in the
trial measures chromosome number, missegregation or aneuploidy. It is the wrong population
and the wrong exposure length: immunocompetent adults with two prior skin cancers, for 12
months, against a child exposed from early childhood for a decade or more. It did not
replicate: the same drug, dose and endpoint in 158 organ-transplant recipients gave 207
new keratinocyte cancers on nicotinamide against 210 on placebo, rate ratio 1.0, 95% CI
0.8 to 1.3 (Allen et al., *N Engl J Med* 2023, PMID 36856616) — underpowered, stopped
early for recruitment, and still a phase 3 null in a higher-risk group.

And the fourth count is the one that matters here. The gate cleared nicotinamide on
evidence orthogonal to the agent's *own* most plausible harm mechanism. Nicotinamide is a
sirtuin product inhibitor; SIRT2 is what holds BubR1 K668 deacetylated (North et al.,
*EMBO J* 2014, PMID 24825348); and Sirt2-null mice show centrosome amplification,
aneuploidy and tumours (Kim et al., *Cancer Cell* 2011, PMID 22014574). A skin-cancer
chemoprevention result says nothing about that axis. The cell is now blank, which is what
"nobody has measured this" is supposed to look like (GP-14) — and blank blocks.

No micronucleus, chromosome-aberration, sister-chromatid-exchange or aneuploidy
measurement of therapeutic nicotinamide was found in any species, in humans or in animals.
That is a not-found and not a demonstrated absence: the searches run are listed in
`docs/research/r4-niacin-safety.md` §9 (queries S21, S22, S23 against PubMed) and
`docs/research/r2-nicotinamide-sirt2.md` §8.

### 5.2 The gate that should not have been shut

An earlier version of this section said niacin's cancer question was "answerable from
existing data" without a citation, was corrected to say the opposite — that it was a data
request against unreported safety capture — and the correction was also wrong. The data
are reported, and public today.

Both large modern niacin trials tabulate malignancy by randomised arm. HPS2-THRIVE
(Landray et al., *N Engl J Med* 2014, PMID 25014686; 25,673 adults, median 3.9 years)
posts neoplasms benign, malignant and unspecified as a MedDRA System Organ Class serious
adverse event on ClinicalTrials.gov NCT00461630: **751/12,838 on extended-release
niacin-laropiprant against 749/12,835 on placebo** over a mean 3.6 years, a rate ratio of
1.00 (95% CI 0.91–1.11, computed here from the posted counts). AIM-HIGH publishes the same
category in a dedicated safety letter: **72/1,718 against 73/1,696, P = 0.87** (Anderson et
al., *N Engl J Med* 2014, PMID 25014706). Roughly 1,645 malignancy events across 29,087
adults, split almost exactly evenly.

On the question a chromosomal-instability syndrome actually cares about, niacin is the
better-characterised of the two agents, not the worse. Nicotinic acid supplementation
raised NAD⁺ and **lowered** micronucleus frequency in primary human peripheral blood
mononuclear cells after X-irradiation (Weidele et al., *DNA Repair (Amst)* 2017, PMID
28216063). It has been through micronucleus and sister-chromatid-exchange assays in rat
bone marrow in vivo (Spronck and Kirkland, *Mutat Res* 2002, PMID 12379464). For
nicotinamide there is no equivalent measurement at all.

Three limits travel with that `0` and the row carries all three. Spronck and Kirkland's own
conclusion is that "pharmacological intakes of niacin do not appear to be further
protective compared to adequate intakes" — deficiency destabilises chromosomes and
adequacy protects, but above adequacy more bought nothing, which is a direct caution
against the supraphysiological-dosing hope. **Both niacin megatrials enrolled zero
participants aged 18 or under**, mean age 64.9 and 63.7 years. And neither agent has ever
been measured in a chromosomal-instability background, so the cell reads
assessed-and-negative *in adults*, never cleared for a child with germline aneuploidy
(ASSUMPTION-DRUG-06).

State exactly what those tables are. A serious-adverse-event category by MedDRA System
Organ Class is a safety-reporting convention, not a pre-specified oncology endpoint: cancer
appears among neither trial's primary nor secondary outcomes on ClinicalTrials.gov, the
counts are not adjudicated for malignancy, and there is no incidence or mortality analysis.
What they establish is bounded and worth writing down without inflation — **malignancy
counts by randomised arm exist, are public, and are null in adults over three to four
years.** A registry-linked cancer-incidence analysis by arm would be stronger, and asking
for one is a real data request. It is no longer what stands between this agent and the
gate.

---

## 6. The experiments

A hypothesis with no stated disconfirming experiment is not a scientific claim.

### 6.1 The falsifying experiment

Patient-derived fibroblasts or LCLs (biobanked lines, not the index child), isogenic
*BUB1B*-hypomorphic hTERT-RPE1 lines representing each published mutation class, and an
allele-corrected control, treated across a clinically achievable NAD⁺-precursor exposure
range and scored blind. Nine readouts run in causal order, from intracellular NAD⁺ to the
fate of the aneuploid cells; Appendix D gives each with its reason, including why the
missegregation *rate* is the endpoint and aneuploidy burden is not.

The precursor is **nicotinic acid**, because it moves one arm of the mechanism and its
result is therefore readable (§4.4). **Nicotinamide runs as a high-dose arm and it is a
negative control, not a second candidate**: above the measured SIRT2 IC₅₀ range the
mechanism predicts BubR1 should *fall*, and a fall would support the axis while
disqualifying that molecule as the therapy. A rise on nicotinamide at high dose would be
evidence against the SIRT2 route rather than for the drug.

Two gates come before the readouts and both can invalidate the run rather than merely
qualify it. Confirm NAPRT expression in each line used, by immunoblot or by
nicotinic-acid rescue of a NAMPT inhibitor, because a NAPRT-negative line cannot convert
nicotinic acid to NAD⁺ and would produce a null that means nothing (§4.5) — this matters
most for any neural arm, where the enzyme is probably absent. And confirm the allele class
first (§6.3), because an allele that yields no transcript is a different experiment.

The hypothesis is falsified if NAD⁺ rises and K668 acetylation falls while total BubR1 does
not; if BubR1 rises and the missegregation rate does not fall; or if the rescue persists
under SIRT2 knockdown. Two results bound the claim without falsifying it, and are kept out
of that list on purpose: an effect confined to the abundance-lowering allele class, which
contracts the eligible population; and a fall in aneuploidy burden without a fall in the
measured missegregation rate, which is not a rescue. Controls: ectopic BubR1 re-expression
as the ceiling, SIRT2 knockdown for mechanism specificity, the allele-corrected isogenic
line as the floor.

### 6.2 The confirming experiment

A monotone dose–response in which precursor exposure raises BubR1 and produces a
proportional, blinded reduction in the directly measured missegregation rate, abolished by
SIRT2 knockdown and not exceeded by ectopic BubR1. Then, and only then, in vivo in
`Bub1b^H/H` mice at a clinically scalable exposure, with missegregation rate and tumour
incidence as co-primary endpoints and aneuploidy burden secondary, read jointly with
aneuploid-cell survival. Tumour incidence is co-primary because in a cancer-predisposition
syndrome a healthspan-only readout cannot answer the question that matters, and Baker 2011's
genetic senescent-cell clearance is the cautionary example of an intervention that looked
good on healthspan while the cancer question sat unasked.

### 6.3 The experiment that comes before either

Determine the allele class (ASSUMPTION-DRUG-09). Quantify *BUB1B* transcript and BubR1
protein in the patient's own cells against a control and, if abundance is reduced, run the
Suijkerbuijk 2010 functional replacement assay to establish whether the residual protein is
quantitatively or qualitatively defective. That single assay decides whether this hypothesis
is applicable, and it is run by the team holding the sample.

Two cheap measurements belong in the same session because they gate the interpretation of
everything after them, and neither exists anywhere in the literature for these cells:
**NAPRT** in the patient lines, which decides whether the probe can act at all (§4.5), and
a baseline **micronucleus frequency** under each agent in a *BUB1B*-deficient background,
which is the measurement whose absence blocks every NAD⁺ precursor in this catalogue
(§5.1).

---

## 7. Limitations

### 7.1 What the model could not express

Three gaps, worked out in Appendix E.1. Surveillance has no home in a drug catalogue: it
moves no node, it changes stage at diagnosis, and *BUB1B*-MVA clears the >5% Wilms-risk rule
that triggers renal ultrasonography every 3–4 months (Scott et al., *Arch Dis Child* 2006,
PMID 16857697). If one thing here has near-term value for a real child it is that
recommendation, and it is not a drug — a judgement labelled **inference**, undercut by the
same paper, which says there is no clear evidence screening reduces mortality or morbidity.
Agent-specific toxicity has no field, so gentamicin's ototoxicity and the bone growth events
in 4% of children on dasatinib (Gore et al., *J Clin Oncol* 2018, PMID 29498925) live in
`notes` cells. And single-node scoring cannot see cross-node harm: tanespimycin was rejected
at the stress node when HSP90 inhibition strips, selectively in the mutant, the exact protein
this hypothesis tries to raise.

### 7.2 The weakest link

L12. That progenitor attrition produces growth restriction and microcephaly is inferred, not
demonstrated in human *BUB1B*-mutant tissue, and there is a published alternative route
through ciliogenesis (Miyamoto 2011). Everything this document says about the growth
phenotype rests on it. It matters less than it looks, because the developmental-window
caveat already removes growth from the achievable endpoints (ASSUMPTION-MECHANISM-02):
growth restriction in MVA is established in utero (Bohers 2008; Hanks et al., *Nat Genet*
2004, PMID 15475955), though that the microcephaly is structural at birth and irreversible
afterwards is itself an **inference**, and Appendix E.2 shows the animal evidence for
irreversibility is weak in both directions. What post-natal correction could plausibly
change is the rate at which new aneuploid cells are generated, so a missegregation-rate or
tumour-incidence endpoint is arguable and a head-circumference endpoint is not. The link
that would hurt most if it broke is L03, that BubR1 abundance is the governing quantity, and
it is the best-evidenced link in the chain.

### 7.3 One child, and the grade of this work

A mechanism-grounded hypothesis is not evidence of efficacy in an individual. Everything
above is reasoning about a *gene* from case series, cell lines and mouse models; the largest
human dataset behind any of it is a handful of MVA patient cell lines; and there has never
been a clinical trial in MVA, nor will there be one powered conventionally. A mechanism
correct in aggregate can be irrelevant in an individual, particularly here, where the
published allele classes predict opposite answers to "would more protein help?" and the
class is unknown to me by design. The n-of-1 designs that could generate real evidence need
a validated, responsive biomarker, which §6 is partly an argument for building.

Under GP-20 the chain and catalogue are **literature-derived**, a third maturity category
defined in ADR 0022, distinct from this repository's synthetic demo tables and its
machine-acquired public reference data. Their guarantee is "every row resolves to a cited
publication", enforced by two tests rather than by a release hash. The second tier cannot be
asserted by a test at all, and `docs/track2-claims-ledger.md` — the record of someone reading
the papers instead — failed 25 times out of 77. A table can be fully cited and still be wrong
in a quarter of its sentences, and no amount of green CI will say so.

### 7.4 The misfolding half of the turnover problem

The mechanism in §4.3 stabilises protein that already exists: SIRT2 holds K668
deacetylated, and deacetylation "is only able to stabilize the protein that is present"
(North et al., *EMBO J* 2014, PMID 24825348). That is a claim about how fast a folded BubR1
molecule is degraded, and it is not the only reason the abundance-lowering alleles are
scarce. The kinase-domain substitutions behave as HSP90 clients — their levels are severely
decreased within five hours of HSP90 inhibition, and cycloheximide and geldanamycin
together remove virtually all mutant protein while wild-type shows only minor changes
(Suijkerbuijk et al., *Cancer Res* 2010, PMID 20516114) — which is the signature of
chaperone-dependent triage rather than of a folded protein being tagged too fast.

Nothing establishes that holding K668 deacetylated counteracts misfolding-driven turnover.
No study measuring BubR1 or BUB1B protein abundance under proteostasis, autophagy, mTOR or
mitochondrial manipulation was found in any species
(`docs/research/r8-proteostasis-axis.md` §8), and that is a not-found rather than a
demonstrated absence (GP-14). So the proposal may be acting on one route of a two-route
turnover problem in the exact allele class it selects for, and §6.1 already names the result
that would expose it: NAD⁺ rises, K668 acetylation falls, total BubR1 does not.

Stating this does not move the proposal onto the proteostasis axis, and it is not an
argument that the axis should be the target instead. Appendix E.1 records why: if the
substitution mutants are cleared because they misfold, the agents that would act there are
chaperone-supporting or proteasome-sparing, and no approved agent obviously fills that
space.

### 7.5 Nobody has measured NAD⁺ in an aneuploid cell

SIRT2 is an NAD⁺-dependent deacetylase, so the axis in §4.3 can only run as fast as the
NAD⁺ supply in the cells it has to act on. Those cells are aneuploid, and their NAD⁺
concentration has not been measured. A PubMed sweep of `aneuploidy AND "NAD+" AND cells`
returned 70 records, none of which, on title screening, reports an NAD⁺ measurement in an
aneuploid or chromosomally unstable background (`docs/research/r8-proteostasis-axis.md` §8,
Q14). That is a not-found and not a demonstrated absence (GP-14): it is not evidence that
NAD⁺ in those cells is normal, and it is not evidence that it is low.

What makes the gap worth stating is the mitochondrial finding beside it. In a *Drosophila*
model of MVA-driven microcephaly, neuroblasts made chromosomally unstable showed compromised
mitophagy and reduced mitochondrial-genome transcription
(González-Blanco et al., *Nat Commun* 2026, PMID 41820377). If mitochondrial dysfunction of
that kind lowers NAD⁺ in the cells this hypothesis most wants to treat, SIRT2 activity falls
with it and the axis is weakest in the compartment that needs it most. Nothing joins those
two statements in one experiment: that model induces instability by depleting *bub3* or
*rod*, and it measures neither NAD⁺ nor BubR1, so the connection is an **inference**.

The two NAD⁺ questions in this document point in opposite directions, and the first
experiment settles which one applies. Appendix C's strongest objection is that a young child
has not lost NAD⁺, so a precursor has nothing to push against; this one is that the
aneuploid cells may have lost it through mitochondrial dysfunction. Readout 1 of §6.1
measures intracellular NAD⁺ in the patient lines before anything else.

---

## 8. Scalability

A new disease is two tables. Nothing under `src/` was added or changed for BUB1B: the gene
contributed three TSVs (12 nodes and 13 links, 13 catalogue rows, one metadata row), a
`SOURCES.md`, a deterministic generator, a test module and two decision records. That module
loads the tables through the same public API and runs the same `build_mechanism` and
`generate_drug_hypotheses` functions the synthetic demo runs, so the direction gate, the
compensatory-node rule, the tri-state agreement and the safety context detectors were all
written for a fictional gene and worked unmodified on a real one. Appendix F lists the
reusable components and three routes to curating signs at volume.

The obstacle to scaling is not the adapter boundary, since swapping the TSV catalogue for
ChEMBL, DrugBank or Open Targets is a boundary change rather than a logic change (GP-02). It
is that `observed_direction` is not a column in any public drug database: ChEMBL records
targets and potencies, Open Targets records associations, and none records which way a
compound pushes a node, which is what the whole gate depends on.

Mechanism is resolved per gene and never per disease name (ASSUMPTION-MECHANISM-03), and MVA
is the example: it is also caused by *CEP57* (Snape et al., *Nat Genet* 2011, PMID 21552266)
and *TRIP13* (Yost et al., *Nat Genet* 2017, PMID 28553959), so everything in §4 is
irrelevant to two of the three and the same catalogue would produce different verdicts
against a *CEP57* chain. The class this generalises to is any monogenic disorder that is a
dosage problem rather than a signalling problem, where the same three questions apply: is the
relationship graded, where does it break, and is there an approved agent that moves the
quantity the right way? For each of them the oncology literature will offer compounds that
bind the right protein and push the wrong way, and the gate in §3 is what stops them.

---

## 9. Methods, as the submission form asks for them

The Track 2 form ships a ten-question methods questionnaire with no upload slot for the
answers, so they live here. The questions are quoted verbatim, and where the honest answer
is "none" it says so.

**1. Team name.** mva-research. The repository carries no separate team identity: the
project name in `pyproject.toml` is `mva-research`, and all 34 commits are authored `MVA
Research`.

**2. Please describe your approach in detail (how you went from variant/mechanism to
candidate medication(s)).** Not from a variant. Track 2 begins where Track 1's ranked pair
ends, at a gene symbol, and the coordinates behind it were never read, requested or
inferred (ADR 0008). Six steps: *BUB1B* as a gene, everything downstream public literature;
a 12-node, 13-link chain in `mechanisms.tsv`, each link graded and cited to a PMID; choice
of target node — `N3_bubr1`, not `N5_sac` — with the direction triple authored in three
separate cells and validated for agreement at construction, so a mistyped cell fails loudly
instead of inverting the gate (ADR 0025); a 13-row catalogue carrying a named
`target_node_id`, a signed `observed_direction`, an evidence type, paediatric exposure, a
`worsens_cin` cell left blank when nobody has asked, and the experiment that would move it;
deterministic triage, asking per row whether the target is on the chain, which way that node
must move and whether the agent moves it that way, then applying gates in fixed priority
order and scoring survivors (direction 0.40, evidence 0.30, safety 0.20, approval 0.10); and
two lists, both kept. What arrived is a ranked list whose top-scoring agent is rejected and
whose accepted one carries `directions_agree = None`, not a medication. The chain is
deliberately not wired into the composition root, because a pipeline whose behaviour depends
on which gene it is looking at is not a pipeline.

**3. Was candidate identification automated (e.g., via a drug-target database query or
LLM-assisted search), or primarily manual literature review and expert curation?** Both,
and the automation is in the checking rather than in the finding. No drug-target database
was queried — there is no ChEMBL, DrugBank, Open Targets or ClinicalTrials.gov client
anywhere in `src/` or `tools/`, because `observed_direction` is not a column any of them
carries (§8). Which agents enter the catalogue, every sign with its citation, and every
mechanism link with its grade came from LLM-assisted literature curation; the target node,
the rule that unassessed CIN risk blocks, and the weights are human-reviewable decision
records (ADR 0025, ADR 0011, ASSUMPTION-SCORING-01). Everything after that is plain Python
with no model inference (ADR 0003), byte-identical under a fixed clock (GP-30); Appendix H
tabulates the boundary. What the pipeline structurally could not do: produce a sign nobody
curated, see harm at a node other than the one an agent targets, record agent-specific
toxicity, distinguish soma from tumour (ASSUMPTION-DRUG-08), or score a non-drug. No model
read patient data at any point, and Track 2 had none to read.

*Which models, on what terms.* The organisers asked entrants to name "the provider, the plan
or tier, and the relevant setting" for any hosted model, in this methods description
(C. Suver, Sage Bionetworks, Community tab, 2026-08-26). Literature curation and drafting
used Anthropic's Claude Code on Claude Opus 5. The adversarial reviews checked in under
`prompts/` were run through OpenAI's Codex CLI on gpt-5.6, and literature sweeps used Exa's
search API. **TODO before submission — account holder to complete:** the plan or tier for
each, and whether training on submitted content is off. Suver flags two specific traps worth
checking against the account rather than assumed: content rated or flagged for safety review
may still be used, and credits or grant programmes can carry terms that override account
defaults. Nothing sent to any of the three was patient-derived. The case volume is outside
what these tools may read (ADR 0008, GP-40), and Track 2 starts from a gene symbol.

**4. Please describe any manual literature review or expert curation in detail.** 61
references in `knowledge/literature/bub1b/SOURCES.md`, each resolved against PubMed on
title, authors, journal, volume, pages and PMID, none recalled from memory — an *identity*
check, with the separate content check in `docs/track2-claims-ledger.md`. Ontology
identifiers were resolved live against HGNC, UniProt, Ensembl, GO, HPO and MONDO, and two
tests make "every row is cited" a machine-checked property. Candidates were found by class
rather than by query: for each class a mechanistic reader would reach for, we looked for the
publication that would *kill* it, and two catalogue rows exist because that read succeeded —
vorinostat, where the acetylation argument reverses once someone measures the mitotic
readout (Magnaghi-Jaulin et al., *Cancer Res* 2007, PMID 17616695), and reversine, where the
oncology logic inverts in a non-tumour soma. This is not a systematic review: no
pre-registered search string, no PRISMA flow, no coverage statement, no second reviewer, no
inter-rater agreement, and no credentialed pharmacologist or clinical geneticist reviewed
these rows. The adversarial reviews are checked in as re-runnable briefs under `prompts/`;
the pharmacology brief ran against the *synthetic* catalogue, so it reviewed the logic rather
than the biology, and it still produced the most consequential change in Track 2 — it found
the prose claiming an unassessed CIN risk was blocking while the code let such a candidate
through at rank 2, and ADR 0011 made the code true rather than deleting the sentence.

**5. Please state whether you relied solely on publicly available drug/target databases and
literature (preferred), or whether proprietary data sources were also used.** Solely public:
peer-reviewed literature indexed in PubMed, public ontology and nomenclature authorities,
and one public regulatory record. No licensed database, no commercial dataset, no internal
corpus, and no patient data.

**6. Please describe the public data sources used in detail (e.g., DrugBank, ChEMBL,
ClinicalTrials.gov, OMIM, PubMed literature search).** PubMed supplied the entire evidence
base — 60 of the 61 references, each with a PMID. HGNC supplied HGNC:1149 (*BUB1B*,
15q15.1), UniProt O60566 and Ensembl ENSG00000156970; Gene Ontology supplied GO:0033597,
GO:0007094 and GO:0008608 as the `node_identifier` cells; HPO supplied HP:0001511 and
HP:0002664; OMIM, MONDO and Orphanet supplied the disease label (OMIM 257300,
MONDO:0009759, Orphanet:1052) for identification only. The one non-literature record is the
European Medicines Agency's public entry for Translarna (ataluren), whose authorisation was
not renewed. DrugBank, ChEMBL, Open Targets and ClinicalTrials.gov were not queried: the
first three record targets, potencies and associations but not which way a compound pushes a
node, and clinical-trial evidence entered through the published reports rather than registry
entries, because the report carries the result and the registry entry does not.

**7. Please describe any proprietary data sources used in detail.** None were used.

**8. How did you characterize the variant's mechanism (loss-of-function vs
gain-of-function, pathway disrupted, downstream biological consequence)?** Loss of function,
and the qualifier is load-bearing: a *dosage* deficiency rather than a null. Biallelic MVA
alleles pair a null or truncating allele with a hypomorphic missense or regulatory one, so
cells carry a reduced but non-zero amount of BubR1 (Suijkerbuijk 2010). Three results
establish that the deficient quantity is the protein rather than any particular activity:
graded knockdown in HeLa places aneuploidy below 50% residual protein with PCS across the
whole range (Bohers 2008), a bound rather than a located threshold for the reasons in §4.1;
an intergenic substitution 44 kb upstream of a transcription start site reproduced the
syndrome when edited into cultured human cells (Ochiai et al., *Proc Natl Acad Sci U S A*
2014, PMID 24344301); and ectopic BubR1 restored checkpoint activity in patient cells. Gain
of function is excluded by that same rescue and by BubR1 being a pseudokinase whose retained
catalytic triad is dispensable for error-free segregation (Suijkerbuijk 2012). Two pathways
are disrupted, not one — an MCC arm restraining APC/C-Cdc20 and a KARD–PP2A-B56 arm
stabilising kinetochore attachments — with only the abundance node upstream of each (ADR
0025). The downstream consequence chain is §2's, ending in growth restriction and
microcephaly on one branch and non-monotonic embryonal tumour risk on the other. This
characterises the *gene and its published allele classes*, not this child's variants, which
were never seen (ADR 0008), so the hypothesis is genotype-conditional (ASSUMPTION-DRUG-09)
and §6.3 runs first.

**9. Please provide an estimate of time or effort spent on this analysis (if possible).**
Elapsed wall clock for the whole repository is 18 h 21 min — first commit 2026-08-27 22:52
+0530, last commit 2026-08-28 17:13 +0530, 34 commits across two calendar days. Track 2 is a
subset: 8 of those commits touch `knowledge/literature/`, `src/mva/mechanisms/`,
`src/mva/interventions/`, `tools/build_knowledge/`, this document or its test module, and it
added 25 chain rows, 13 catalogue rows, a `SOURCES.md`, a deterministic generator, a 22-test
module, two decision records, this document and its two companions, and no line under
`src/`. Elapsed wall clock is not person-hours, and the literature reading is not separable
from the commit record because it happened in the same window; this is an agent-assisted
build, so a headcount-hours figure would be invented rather than derived.

**10. Provide a method abstract (up to 500 words), including strengths and limitations of
your approach as well as methodological detail.** Track 2 started from a gene, not a variant.
No genotype, coordinate or clinical record was read, requested or inferred (ADR 0008), and
the hypothesis is stated genotype-conditionally because the allele class is unknown to us by
design. We curated a 12-node, 13-link mechanism chain into a TSV, grading each link
independently for tier, strength, direct demonstration and stated uncertainty, and citing
each to a PMID. Every claim was afterwards re-checked against what the cited paper says
rather than against the identifier; 25 of 77 claims were weakened, withdrawn or left
unverified. Two links are labelled inference. One node — p53/p16 arrest, senescence and
immune clearance — is flagged compensatory, and the model refuses to designate it a
therapeutic target, because pushing it back would suppress a response holding a
cancer-predisposed child's aneuploid cells in check. The therapeutic target is BubR1 protein
abundance, not the checkpoint the disease is named after: BubR1 reaches missegregation by
checkpoint restraint of APC/C-Cdc20 and by PP2A-B56 recruitment that stabilises kinetochore
attachments, only the protein is upstream of both, and the checkpoint node's entire
accessible pharmacology is inhibitors, which push the disease direction. Abundance is a floor
rather than a severity predictor: mouse genotypes with indistinguishable total BubR1 and
near-identical aneuploidy differ 1.75-fold in median lifespan (Sieben 2020), so total protein
is the target-engagement readout and not the outcome measure. Thirteen real agents were
hand-curated into a catalogue, each with a signed observed direction at a named node, cited;
no drug-target database was queried, because the sign is not a column any of them carries. A
deterministic Python stage then asked three separate questions per row — is the target on the
chain, which way must that node move, does the agent move it that way. Agreement is
tri-state, and "cannot determine" is never read as agreement. Gates follow in priority order:
a wrong direction disqualifies, and in a chromosomal-instability disease so does an
unmeasured effect on chromosomal instability. Output: 1 accepted, 12 rejected, 89 evidence
rows, every rejection kept with its reasons and the experiment that would unblock it;
byte-identical under a fixed clock, with 26 tests asserting the conclusions. **Strengths.**
The sign is the primary gate, so the compounds a target-proximity search ranks highest — MPS1
inhibitors — come out last. Gates are not score terms, so the highest-scoring agent, the only
one with a measured agreement, is still rejected, and the accepted candidate carries an
undetermined direction rather than a manufactured one. Nothing under `src/` was written for
this gene, so a new disease is two tables. **Limitations.** The catalogue is n=13 and
hand-curated, with no pre-registered search, no second reviewer and no credentialed
pharmacologist. Scoring is single-node, so cross-node harm is invisible to it, and
agent-specific toxicity has no field. One context is scored — the soma — so every rejection
means "not for systemic use in the affected child". The weakest link is inferred and carries
the entire growth rationale. And no approved agent has a measured effect on BubR1 abundance,
which is why the deliverable is a falsifying experiment rather than a medication.

---

## 10. Reproducing this

```bash
uv run pytest -q tests/unit/test_track2_bub1b.py   # 26 tests: the chain, the gate, the citations
just verify                                        # the full acceptance gate
```

| Artefact | Path |
|---|---|
| Mechanism chain (12 nodes, 13 links, individually graded) | `knowledge/literature/bub1b/mechanisms.tsv` |
| Target node and required correction | `knowledge/literature/bub1b/mechanism_meta.tsv` |
| Drug catalogue (13 agents, signed directions) | `knowledge/literature/bub1b/drug_catalog.tsv` |
| Every citation, mapped row by row | `knowledge/literature/bub1b/SOURCES.md` |
| Every claim re-checked against its source, with verdicts and corrections | `docs/track2-claims-ledger.md` |
| The class base rate and the five-condition conjunction, without a fabricated posterior | `docs/track2-quantified-uncertainty.md` |
| The allelic-series evidence behind §4.2 and Appendix A.4 | `docs/research/r1-allelic-series.md` |
| The dose–response evidence behind §4.1 and Appendix A.3 | `docs/research/r3-bubr1-threshold.md` |
| Catalogue generator | `tools/build_knowledge/bub1b_catalog.py` |
| Executable assertions | `tests/unit/test_track2_bub1b.py` |
| Why literature tables sit outside `knowledge/public/` | `docs/decisions/0022-literature-tables-are-a-third-knowledge-category.md` |
| Why the target is abundance, not the checkpoint | `docs/decisions/0025-the-target-node-is-bubr1-abundance-not-the-checkpoint.md` |
| Why an unassessed CIN risk blocks | `docs/decisions/0011-unassessed-cin-risk-is-blocking.md` |
| The assumptions this rests on | `docs/scientific-assumptions.md` |

---

# Appendices

Everything below supports a claim in the body. It is here rather than there because the
body is the argument and this is the working.

## Appendix A — The mechanism chain

### A.1 Demonstrated links

| Link | Claim | Evidence |
|---|---|---|
| L01 | Truncating and upstream-regulatory *BUB1B* alleles reduce full-length transcript | Suijkerbuijk et al., *Cancer Res* 2010;70:4891-900 (PMID 20516114) reports absence of transcript from truncating mutants. Ochiai et al., *PNAS* 2014;111:1461-6 (PMID 24344301) identified an intergenic substitution 44 kb upstream of a *BUB1B* transcription start site and reproduced the syndrome by editing that single base pair into cultured human cells, which produced reduced transcript, increased PCS and MVA. |
| L02 | Reduced transcript plus accelerated turnover of the missense protein lowers BubR1 abundance | Suijkerbuijk 2010: "Low protein abundance is the direct result of the absence of transcripts from truncating mutants combined with high protein turnover of missense mutants." |
| L03 | Insufficient BubR1 impairs the mitotic checkpoint complex, and the effect is dosage-graded | Suijkerbuijk 2010 (MVA patient cell lines show impaired checkpoint, alignment defects and low BubR1). Bohers et al., *Hum Genet* 2008;124:473-8 (PMID 18932004) produced a residual-BubR1 gradient (8.5 / 10 / 14 / 58 / 77%) in shRNA-transduced HeLa — a transformed, already-aneuploid line — and found PCS across the whole range, with aneuploidy clearly detected below 50% residual protein. The relationship is a graded dose–response rather than a step, and its position differs by readout; §4.1 and A.3 give the reconciliation. |
| L04 | The checkpoint complex restrains APC/C-Cdc20; failure permits premature anaphase | Lara-Gonzalez, Pines & Desai, *Semin Cell Dev Biol* 2021;117:86-98 (PMID 34210579); Kops, Weaver & Cleveland, *Nat Rev Cancer* 2005;5:773-85 (PMID 16195750). |
| L05 | BubR1 stabilises kinetochore–microtubule attachments by recruiting PP2A-B56 — a route independent of the checkpoint | Suijkerbuijk et al., *Dev Cell* 2012;23:745-55 (PMID 23079597); Kruse et al., *J Cell Sci* 2013;126:1086-92 (PMID 23345399); Lampson & Kapoor, *Nat Cell Biol* 2005;7:93-8 (PMID 15592459). |
| L06–L08 | Checkpoint failure → premature chromatid separation → whole-chromosome missegregation; destabilised attachments do the same by a second route | Bohers 2008; Suijkerbuijk 2010; Matsuura et al., *Am J Med Genet A* 2006;140:358-67 (PMID 16411201). |
| L09 | Aneuploidy imposes proteotoxic and metabolic stress | Oromendia, Dodgson & Amon, *Genes Dev* 2012;26:2696-708 (PMID 23222101); Donnelly et al., *EMBO J* 2014;33:2374-87 (PMID 25205676); Santaguida et al., *Genes Dev* 2015;29:2010-21 (PMID 26404941); reviewed in Santaguida & Amon, *Nat Rev Mol Cell Biol* 2015;16:473-85 (PMID 26204159). |
| L10 | Missegregation triggers p53/p16 arrest, senescence and immune clearance | Thompson & Compton, *J Cell Biol* 2010;188:369-81 (PMID 20123995); Santaguida et al., *Dev Cell* 2017;41:638-651 (PMID 28633018); Li et al., *PNAS* 2010;107:14188-93 (PMID 20663956); Baker et al., *Nat Cell Biol* 2008;10:825-36 (PMID 18516091). The clearance arm is demonstrated as pro-inflammatory signalling plus NK-cell killing **in co-culture**; the extension to in-vivo human immunosurveillance, which §7.1 and Appendix E lean on, is an **inference**. |
| L13 | Aneuploidy raises cancer risk, non-monotonically | Hanks et al., *Nat Genet* 2004;36:1159-61 (PMID 15475955); Jacquemont et al., *Am J Med Genet* 2002;109:17-21 (PMID 11932988) — 3 of 14 reported MVA cases (21%) developed a malignancy: rhabdomyosarcoma, acute lymphoblastic leukaemia and nephroblastoma. The non-monotonicity is from Weaver et al., *Cancer Cell* 2007;11:25-36 (PMID 17189716) and Silk et al., *PNAS* 2013;110:E4134-41 (PMID 24133140): low-rate missegregation promotes tumours, very high-rate missegregation suppresses them. |

### A.2 The two inferred links, and the compensatory node

**L11 — aneuploidy-associated stress *causes* the arrest/senescence response.** Both
phenomena are documented; that one causes the other, rather than both being parallel
consequences of missegregation, is not established.

**L12 — progenitor attrition produces growth restriction and microcephaly.** The weakest
joint in the chain, carrying the entire therapeutic rationale for the growth phenotype.
Baker et al., *Nat Genet* 2004;36:744-9 (PMID 15208629) establishes that BubR1
insufficiency produces dwarfism and progeroid phenotypes, but by what route is not shown
in human *BUB1B*-mutant tissue, and there is a published alternative: BubR1 insufficiency
also impairs ciliogenesis (Miyamoto et al., *Hum Mol Genet* 2011;20:2058-70, PMID
21389084), a developmental mechanism with nothing to do with progenitor depletion. A
cell-autonomous growth defect is a third possibility. Which tissues and which
developmental windows dominate is unresolved.

**The compensatory node, in full.** In BubR1-hypomorphic mice, p16-driven senescence
*drives* the progeroid phenotype and deleting p16 attenuated it (Baker 2008, PMID
18516091), while clearing p16-positive senescent cells delayed age-related deterioration
(Baker et al., *Nature* 2011;479:232-6, PMID 22048312, performed in that same BubR1
progeroid model). The same p53 arm *suppresses* aneuploidy-driven tumorigenesis (Li 2010,
PMID 20663956), and immune clearance removes missegregating cells (Santaguida 2017, PMID
28633018 — in-vitro NK co-culture; the in-vivo human reading is an inference). One
response, pathological for growth and ageing, protective against cancer. This is exactly
the case where a signed state field cannot say which, which is why the direction check
returns `UNKNOWN` and the model refuses the node as a target (ASSUMPTION-MECHANISM-04).

### A.3 The dose–response, readout by readout

Ranked by how much BubR1 can be lost before the readout moves, most sensitive first.
Compiled in `docs/research/r3-bubr1-threshold.md`, which also records the searches that
found nothing.

| Readout | Where it moves | Evidence |
|---|---|---|
| PCS / centromeric cohesion | already at ~50–77% residual | Bohers 2008: PCS "detected in all transduced cells", the least-depleted clone being 77%. García-Castillo et al., *Am J Med Genet A* 2008 (PMID 18548531) found PCS in 7 of 7 human heterozygous carrier parents, who are clinically unaffected. |
| Accumulated aneuploidy | somewhere below 58%; clear at ≤14% | Bohers 2008. The "below 50%" figure is an interpolation across a gap the series does not sample. |
| Per-division missegregation | below ~13%, gross at ~6% | Suijkerbuijk 2010's doxycycline titration in HeLa. |
| Spindle assembly checkpoint, attachment stability, PP2A-B56 recruitment | **never titrated by anyone** | Suijkerbuijk 2010 tested the checkpoint only in mutant-replacement experiments, not against a dose series. |
| Cortical development (mouse) | tolerates ~1.8%; fails at ~0.5% | Simmons et al., *Hum Mol Genet* 2019;28:1822-1836 (PMID 30668728): `BubR1^H/H` E15.5 cortex sits at about 1.8% of wild-type BubR1 with no microcephaly, and only conditional deletion to about 0.5% produces it. |

Two further reasons "the threshold" cannot be one number. The denominator is wrong:
whole-cell immunoblot percentage is not the functional variable, and both Suijkerbuijk
2010 and Simmons 2019 report that residual protein still loads onto kinetochores or
chromosomes, so a cell at 10% with normal loading is not in the state of a cell at 10%
with none. And level does not determine phenotype: Wijshake et al., *PLoS Genet*
2012;8:e1003138 (PMID 23300461) reports two mouse genotypes with the same residual
wild-type BubR1 that differ in PCS — `+/GTTA` elevated, `+/−` not — at the same modestly
raised aneuploidy rate. One comparison is enough to break the pure-dosage reading of PCS.

### A.4 The allelic series, in full

Sieben 2020 combined four *BubR1* alleles pairwise: `L1002P` (a new knock-in modelling
the human L1012P missense allele), `X753` (a frameshift knock-in producing truncated,
unstable protein), `H` (the pre-existing hypomorph at about 10% of normal) and a null.
`X753/L1002P`, the exact human patient genotype, is embryonic lethal in mice, though the
human patient lived 3.6 months. `H/L1002P` mice are viable with a median lifespan of 343
days against 691 for wild-type; `H/H` mice reach 196 days.

The comparison that matters is `H/L1002P` against `H/H`. Protein reduction "appeared to
be very similar between the mutants" in MEFs and multiple tissues, with a significant
difference only in spleen. Tissue FISH for chromosomes 4 and 7 showed karyotypic
abnormality rates "prominently increased" in both, "with minimal differences observed
between the 2 genotypes". Median lifespan differs 1.75-fold, kyphosis onset 238 against
175 days, and cardiac stress tolerance differs sharply. Three qualifications the headline
elides: several phenotypes are equally severe in both, including lipodystrophy,
cataractogenesis at 161 against 168 days, and carcinogen-induced tumour susceptibility;
in mitotic splenocytes the *more* aneuploid genotype lived longer, which damages an
aneuploidy-centric reading more than a dosage-centric one; and the paper itself identifies
critical thresholds in vivo, with a 40–90% BUBR1 reduction giving tumour predisposition or
mild premature ageing and greater than 90% giving progeroid phenotypes or early postnatal
lethality. So level sets the regime and allele identity modulates severity within it.

**One reading of the paper does not survive its own figures.** The claim that L1002P
protein accumulates normally at unattached kinetochores despite low total abundance —
and therefore that abundance is the wrong variable — rests on a sentence about
`BubR1^+/L1002P` heterozygotes, which carry an intact wild-type allele, measured with an
antibody that does not distinguish wild-type from mutant protein. Nothing in the paper
attributes that kinetochore signal to the mutant. In the biallelic genotype the same
assay reports BUBR1 accumulation at unattached kinetochores severely compromised in both
`H/L1002P` and `H/H`, so kinetochore BubR1 tracks total abundance downward. Human patient
cells agree: kinetochore BUBR1 is low in proportion to total, and the residual protein is
still competent to localise (Suijkerbuijk 2010).

**What follows for the marker stack.** Total BubR1 is the target-engagement readout, not
the severity endpoint. Appendix D lists the readouts that separate the two.

## Appendix B — The therapeutic classes considered

Direction is carried across by class in two rows, and `SOURCES.md` records the
substitution: hydroxychloroquine stands for the chloroquine Tang 2011 screened, and
metformin for AICAR through AMPK activation. Only the tanespimycin row is exact, since
17-AAG is the compound Tang screened.

| Class | Representative | Node | Direction | Desired in **soma**? | Desired in **tumour**? | Outcome |
|---|---|---|---|---|---|---|
| **Restore BubR1 abundance** (NAD⁺/SIRT2 axis) | NMN, nicotinamide riboside, nicotinamide, nicotinic acid | BubR1 abundance | ↑ (NMN measured); unmeasured or bidirectional for the rest | **Yes** | **Yes** | **Primary hypothesis.** Best agent is unapproved; approved agents have undetermined sign. |
| **Readthrough / restore full-length protein** | Ataluren, gentamicin | BubR1 abundance | ↑ *if* a nonsense allele produces a transcript | Yes, conditionally | Yes | Rejected. Genotype-conditional and mechanistically undercut: truncating *BUB1B* alleles yield no transcript (Suijkerbuijk 2010), so NMD removes the substrate before the ribosome reaches the stop codon. Ataluren's EU authorisation was also not renewed (EMA; EC decision effective 2025). |
| **Reduce aneuploidy consequences — proteostasis** | Sirolimus | Proteotoxic stress | Unsigned at the node, but no longer untested in an aneuploid background: neuroblast number rescued, brain size not | Unclear | Unclear | Rejected, and now on a measured negative rather than a missing measurement. Rapamycin at 40 µM, fed to larvae whose neuroblasts had been made chromosomally unstable, rescued the number of aneuploid neuroblasts and rescued neither progeny cell number nor larval brain size, the microcephaly endpoint the paper is named for; the authors read TOR inhibition as a possible harm to the proliferative capacity of highly aneuploid stem cells (González-Blanco et al., *Nat Commun* 2026, PMID 41820377). That model depletes *bub3* or *rod* in fly neuroblasts rather than the BubR1 orthologue and measures no BubR1, so it bears on MVA-class instability and not on *BUB1B* dosage. Also anti-growth in a growth-restricted child and immunosuppressive in a child whose aneuploid cells are partly held in check by immune clearance — that last step is an **inference**: Santaguida 2017 demonstrates pro-inflammatory signalling and NK-cell killing of missegregating cells in co-culture, and the extension to in-vivo human immunosurveillance is ours, not the paper's. Two objections at *other* nodes, which single-node scoring cannot express. |
| **Aneuploidy-selective antiproliferatives** | Hydroxychloroquine, tanespimycin (17-AAG), metformin | Proteotoxic / metabolic stress | ↑ stress (kills aneuploid cells) | **No — catastrophic** | **Yes** | Rejected for the soma with `WRONG_DIRECTION`. Retained in the record as tumour-arm candidates requiring a therapeutic-index study. |
| **Cell-cycle / checkpoint modulation** | Reversine (MPS1i) | SAC | ↓ checkpoint | **No — pushes the disease mechanism** | Yes | Rejected, lowest score of all. The naive top hit. |
| **Mitotic timing — restrain APC/C directly** | proTAME, apcin | SAC / APC/C | ↑ restraint on APC/C-Cdc20 | Yes, in principle | — | Not modelled, and the one class whose direction is right and whose absence is purely regulatory. Partial APC/C inhibition substitutes for the restraint the failing checkpoint cannot supply. Both are tool compounds with no record of human administration — an **inference** from the absence of any registered trial or published clinical report, not a sourced negative — so under ASSUMPTION-DRUG-03 they are research probes; they would add a third row saying only `NOT_APPROVED`, which reversine and tanespimycin already demonstrate. Recorded in `SOURCES.md` under "considered and left out". The obvious objection is a therapeutic window: too much APC/C restraint is mitotic arrest and slippage, which generates aneuploidy of its own. |
| **Acetylation / degradation modulation** | Vorinostat (HDACi class) | SAC | ↓ checkpoint | **No** | — | Rejected. Argued for by analogy — BubR1 is stabilised by K250 acetylation (Choi et al., *EMBO J* 2009;28:2077-89, PMID 19407811) and HDAC inhibitors raise acetylation — but measured at the mitotic readout, the class induces premature sister-chromatid separation and overrides the SAC (Magnaghi-Jaulin et al., *Cancer Res* 2007;67:6360-7, PMID 17616695). A mechanism-shaped argument that reverses when someone measures the actual readout. |
| **Senolytics** | Dasatinib (with quercetin; Zhu et al., *Aging Cell* 2015;14:644-58, PMID 25754370) | Clearance / senescence (**compensatory**) | ↓ senescent-cell burden | **Cannot determine** | Possibly harmful | Rejected. Motivated by the best animal evidence in this document — Baker 2011's genetic senescent-cell clearance in the exact BubR1 model (INK-ATTAC, a dimerizer-triggered transgene, not a senolytic drug; the pharmacophore line is Zhu 2015) — and blocked by the compensatory-node rule in §2: the response it removes may be a barrier rather than a burden. |
| **Antioxidants** | N-acetylcysteine | — | ↓ ROS | Plausible | **No** | Not modelled. The oxidative-stress node is not separately represented, and the strongest relevant finding points the wrong way for a cancer-predisposed child: antioxidants accelerated tumour progression in mice (Sayin et al., *Sci Transl Med* 2014;6:221ra15, PMID 24477002). |
| **Surveillance** (non-pharmacological) | Renal ultrasound, clinical review | — | Does not move any node | Yes — highest near-term value | — | Deliberately kept out of the drug catalogue (§7.1, Appendix E). |

## Appendix C — The four objections, in full

**"The NAD⁺ evidence comes from an ageing context. A young child has not lost NAD⁺, so
where is the headroom?"** The strongest objection, and I do not have a good answer. North
2014 frames the axis as an explanation for the age-related decline in BubR1; the
intervention restores something that was lost. In a young child NAD⁺ is presumably not
the limiting quantity, and the pharmacology may have nothing to push against. This is the
first thing the experiment in §6.1 measures, and a plausible route to a clean negative.

**"Suijkerbuijk 2010 describes two mutation classes. Half of them are qualitative defects
that more protein will not fix."** Partly, and the split is more favourable than "half".
Of the substitutions that paper characterises, five kinase-domain alleles (I909T, R727C,
R814H, L844F, L1012P) are abundance-lowering, showing a 5–10 fold decrease in protein
level; R550Q "rescued both functions of BUBR1"; Q921H was "indistinguishable from
wild-type LAP-BUBR1, both on protein level and functionality"; and only Y155C is a
genuine separation-of-function allele, restoring chromosome alignment while failing to
reconstitute checkpoint activity. The paper's own conclusion runs the other way from the
objection: chromosomal instability in *BUB1B*-MVA "is a result of low BUBR1 protein
abundance". The hypothesis remains explicitly genotype-conditional (ASSUMPTION-DRUG-09,
§7.3) — the class still has to be determined per patient, which is what §6.3 is for — but
the eligible fraction is larger than a coin flip. The sharper version of the objection is
better evidenced: even in the abundance-lowering class, the kinase-domain mutants are
destabilised by **misfolding**, their levels collapsing under HSP90 inhibition and being
rescued by proteasome inhibition (Suijkerbuijk 2010, Fig. 3C). Whether the K668
deacetylation route can counteract misfolding-driven turnover has never been tested.

**"Is more BubR1 automatically safe?"** No. Mad2 overexpression is tumorigenic in mice
(Sotillo et al., *Cancer Cell* 2007;11:9-23, PMID 17189715), so "more checkpoint protein"
is not monotonically good as a class statement. BubR1 overexpression specifically was
protective (Baker 2013), but any dose-finding has to look for a ceiling as well as a
floor.

**"You are pushing a whole cell toward more deacetylation. Which lysine does that hit?"**
Not answered anywhere in the literature this document cites, and it sits between two of
its own sections. BubR1 is **stabilised** by acetylation at K250 — PCAF acetylates it at
prometaphase, the acetyl-mimetic K250Q resists ubiquitin-dependent degradation and the
acetylation-deficient K250R is unstable (Choi 2009) — which is the entire argument
Appendix B records for HDAC inhibitors. BubR1 is **destabilised** by acetylation at K668,
which is why SIRT2-driven deacetylation raises abundance (North 2014). Two lysines,
opposite signs on the same protein's stability, both correctly cited. Nothing establishes
that a global NAD⁺/sirtuin shift is site-selective for K668 and spares K250; if it is not,
the intervention buys stability at one residue and spends it at another. Readout 2 of
§6.1 therefore measures K250 acetylation alongside K668.

## Appendix D — The falsifying experiment's readouts

1. Intracellular NAD⁺ (LC-MS/MS).
2. BubR1 K668 acetylation **and** K250 acetylation (site-specific immunoprecipitation or
   targeted MS) — K250 because acetylation there stabilises BubR1 and a global sirtuin
   shift may not be site-selective (Appendix C, fourth objection).
3. Total BubR1 protein (quantitative immunoblot plus targeted MS), expressed as a
   percentage of the isogenic corrected control. Target engagement, not outcome
   (Appendix A.4).
4. **Chromosome missegregation rate, measured directly as a rate** — live-cell imaging of
   H2B-labelled cells scoring lagging chromosomes and anaphase bridges *per anaphase*, or
   paired daughter-cell single-cell sequencing with lineage tracking so that each error is
   attributed to a division that was watched. Blinded. Micronucleus rate is scored
   alongside and is not a substitute: micronuclei arise from several processes besides
   whole-chromosome missegregation.
5. Premature chromatid separation frequency, scored separately from aneuploidy rather than
   pooled with it. Sieben 2020 found PCS diverging dramatically between `BubR1^H/L1002P`
   and `BubR1^H/H` at identical aneuploidy, which makes PCS the readout sensitive to a
   mutant protein's interference rather than to its scarcity, and Wijshake 2012 shows PCS
   is not a pure function of BubR1 level.
6. Checkpoint strength (mitotic index under nocodazole challenge). This is the assay in
   which raising L1012P to wild-type levels produced full rescue (Suijkerbuijk 2010), and
   the closest existing thing to a validated response endpoint for that allele.
7. Kinetochore-bound BubR1 normalised to centromere signal, as a check that added protein
   goes where it is needed. Suijkerbuijk 2010, Simmons 2019 and Sieben 2020 all indicate
   that kinetochore- or chromosome-loaded protein, not whole-cell abundance, is the
   functional quantity.
8. Aneuploidy burden (single-cell karyotyping or low-pass single-cell WGS) — **secondary
   and slow**, for the flux-versus-stock reason in §6.1, and kept as a safety and mechanism
   check rather than a severity surrogate.
9. **Survival and proliferation of the aneuploid cells themselves**, analysed separately
   from burden: time-lapse fate tracking of aneuploid daughters, plus burden expressed
   against a division-count denominator rather than against elapsed time.

Two further readouts are worth building, and neither has been done for a *BUB1B* MVA
allele. Allele-resolved BubR1 abundance — targeted MS on the tryptic peptide spanning the
substitution — because total protein cannot distinguish "more wild-type" from "more
stabilised mutant", and the stabilisation mechanism proposed here acts on the mutant pool.
And PP2A-B56 recruitment via KARD phosphorylation, because the pseudokinase domain that
L1012P sits in controls that arm, and a nocodazole checkpoint assay would not surface a
defect there. That is the single highest-value experiment for de-risking the thesis.

## Appendix E — What the model could not express, in full

### E.1 The three gaps

**Non-pharmacological interventions have no home in a drug catalogue.** Tumour
surveillance has to be argued, not delegated to a guideline. Scott et al., *Arch Dis
Child* 2006;91:995-9 (PMID 16857697) recommend surveillance for any child at >5% risk of
Wilms tumour, by renal ultrasonography every 3–4 months, while stating plainly that the
risks and benefits are finely balanced, that there is no clear evidence that screening
reduces mortality or morbidity, and that the recommendations are pragmatic rather than
evidence-based. *BUB1B*-MVA clears the 5% rule on any published series: Suijkerbuijk
2010's introduction puts it at 37% of patients developing cancers "mostly within the first
3 years of life"; García-Castillo et al., *Am J Med Genet A* 2008 (PMID 18548531) reports
Wilms tumour in 7 of 8 patients with monoallelic *BUB1B* mutations; and Jacquemont et al.,
*Am J Med Genet* 2002 (PMID 11932988) reports malignancy in 3 of 14 clinically ascertained
cases (21%) — the weakest of the three as a rate for this child, because the series was
ascertained clinically before *BUB1B* was implicated and the index patient's leukaemia is
inside its own numerator. Surveillance does not move any node on this chain; it changes
*stage at diagnosis*. Forcing it into the catalogue would have produced either a category
error or a spurious `WRONG_DIRECTION` rejection, so it is discussed and excluded. If one
thing in this document has near-term value for a real child, it is this paragraph, and it
is not a drug — a judgement, labelled **inference**, and one the surveillance literature
does not underwrite, since the paper the recommendation rests on is the same paper saying
the mortality benefit is unproven.

**Agent-specific toxicity has no field.** Gentamicin's cumulative ototoxicity and
nephrotoxicity make chronic dosing implausible for a lifelong condition — standard
aminoglycoside pharmacology, carried here as an **inference** rather than a sourced claim,
because no citation is attached to it. Sirolimus is anti-growth in a growth-restricted
child and immunosuppressive in a child whose aneuploid cells are partly held in check by
immunity — that last is also an **inference**, from the in-vitro NK-cell co-culture of
Santaguida 2017 to in-vivo human immunosurveillance. And the paediatric dasatinib trial
reported bone growth and development events in 4% of patients (Gore et al., *J Clin Oncol*
2018;36:1330-1338, PMID 29498925), which is a specific concern in exactly this phenotype.
None of these is expressible in the schema, so they live in `notes` cells and here. The
safety pass says so itself: "Absence of a concern here reflects the fields the catalogue
records, not a safety clearance."

**Single-node scoring cannot see cross-node harm.** An agent is scored against the one
node it targets, and sirolimus's objections all live at other nodes. The sharper worked
example is tanespimycin, because there the missed evidence is not hypothetical. The
catalogue rejects 17-AAG at `N9_stress`, as an aneuploidy-selective antiproliferative. The
stronger objection was sitting at the *target* node the whole time, in a paper the
catalogue already cites: MVA kinase-domain substitution mutants are HSP90 clients, and in
Suijkerbuijk 2010 Fig. 3C "levels of both mutant proteins were severely decreased after 5h
treatments with the HSP90 inhibitor", with cycloheximide and geldanamycin together
removing "virtually all mutant BUBR1 protein" while wild-type showed "only minor changes".
HSP90 inhibition destroys the exact protein the primary hypothesis is trying to raise,
selectively in the mutant — a `WRONG_DIRECTION` at `N3_bubr1`, measured in the disease's
own alleles rather than carried across from a generic aneuploidy model. The rejection
outcome does not change; the reason does, and single-node scoring routed around it. The
same finding has a corollary this document does not otherwise pursue: if the substitution
mutants are cleared because they misfold, then chaperone-supporting or proteasome-sparing
approaches act at the target node by a route with no NAD⁺ in it. No approved agent
obviously fills that space, which is why it is a note rather than a catalogue row — but
"the acetylation axis is the only handle on abundance" is not established, and §4 should
not be read as though it were.

### E.2 Why the animal evidence on growth is weaker than it looks

SIRT2 overexpression extended male lifespan in `BubR1^H/H` mice and left the animals
**small** (North 2014), which reads as evidence that growth is developmentally fixed. Two
things undercut it. BubR1 was never restored to wild-type levels in those animals (§4.3),
so their staying small is as consistent with too little of the mechanism as with
irreversibility, and no published experiment separates the two. And the mouse cortex is a
poor witness in the other direction: `BubR1^H/H` E15.5 cortex sits at roughly 1.8% of
wild-type BubR1 with normal cortical size and progenitor division, and microcephaly
appears only when conditional deletion drops it to about 0.5% (Simmons et al., *Hum Mol
Genet* 2019;28:1822-1836, PMID 30668728). Mouse `BubR1^H/H` does not get microcephaly;
human MVA patients do. So the conclusion that no post-natal increase in BubR1 will reverse
growth restriction or microcephaly rests on the developmental inference alone, not on the
animal data, which is compatible with either story.

## Appendix F — The reusable machinery

| Component | Why it transfers |
|---|---|
| `EffectDirection` + `directions_agree` (tri-state) | A signed vocabulary with an explicit "cannot determine". Applies to any mechanism where a therapy has a sign. |
| The direction-triple validator on `MechanismHypothesis` | `disease_direction`, `required_correction` and the target node's `state_in_patient` are authored independently and must agree at construction. One mistyped cell would otherwise invert the entire gate silently. Gene-agnostic. |
| `deviation_is_pathological` (mandatory, un-defaulted) | The compensatory-node concept. Every disease with a damage response has this problem: DNA-repair disorders (p53 attrition), lysosomal disease (autophagy upregulation), mitochondrial disease (mitophagy). |
| Rejection reasons as a priority-ordered vocabulary | "Wrong direction" and "not approved" are different findings with different follow-ups ("never" vs "wait"). The ordering makes reports say the right one. |
| Context detectors (`is_chromosomal_instability_context`) | Which gaps are *fatal* is disease-dependent. The mechanism itself selects the gate. |
| Evidence-tier ladder + `is_direct_evidence` discount | Keeps a binding assay from being read as a result. Universal. |
| The two-TSV loader contract | A new disease is two tables, not a code change. |

Three tractable routes to signs at volume, since `observed_direction` is in no public
database: harvest MoA verbs from ChEMBL/DrugBank text, where inhibitor/agonist maps
cleanly onto ↓/↑ for a *direct* target and not at all for downstream nodes; use
perturbation signatures (LINCS L1000) to read the sign empirically at a transcriptional
node; or restrict automated ingestion to direct-target rows and require curation for
everything else. The tri-state means the last option degrades gracefully, because unknown
signs become flagged candidates rather than fabricated ones.

## Appendix G — The NAD⁺ axis at the molecule

### G.1 The axis, and the ceiling North 2014 hit

North et al., *EMBO J* 2014;33:1438-53 (PMID 24825348) established the axis: SIRT2 holds
BubR1 K668 deacetylated in an NAD⁺-dependent manner, CBP acetylates it, and the balance
sets turnover. Treating mice with NMN increased BubR1 abundance in vivo — 500 mg/kg/day
intraperitoneally for 7 days, in wild-type mice aged 3 and 30 months, with the effect
framed as restoration of an age-related decline, BubR1 in treated aged animals coming back
to the level of the 3-month-old ones. SIRT2 overexpression in `BubR1^H/H` hypomorphic
animals increased median lifespan by 123% in males, with no change in females, and the
animals remained small.

The ceiling is the part that bears on this hypothesis. SIRT2 overexpression raised BubR1
in the tissues where BubR1 could be detected at all, but the authors report that they
"were not able to see a restoration of BubR1 levels observed in the BubR1^H/H mice back to
wild-type levels", and they say why: the hypomorph's deficit is **transcriptional**, and
SIRT2-mediated deacetylation "is only able to stabilize the protein that is present"
rather than drive more mRNA from the *BubR1* promoter. They then offer that shortfall as
the explanation for their own mixed result — lifespan moved while lordokyphosis, cataract
formation and small size did not — reasoning that SIRT2 overexpression "may not have
achieved a level of BubR1 able to reverse the senescence effects of BubR1 depletion".

Read that against what this document proposes. The closest intervention anyone has run in
a BubR1-deficient animal was genetic, constitutive, present from development, and applied
to the SIRT2 arm directly rather than to NAD⁺ upstream of it — which I take to be the
stronger of the two perturbations, though that is an **inference** and not something
either paper measures. It still could not close the gap, because the lesion it worked
against was deficient *synthesis*, which is the lesion *BUB1B*-MVA has.

Where the effect was seen also matters. The NMN result that carries the strongest
direction claim was read out in **testes**, chosen because they express high levels of
BubR1; in the same experiment BubR1 protein levels "were difficult to determine in the
heart tissue of aged animals". The tissue that supplies the measurement is not one of the
tissues in which MVA does its damage, and no reading of the paper closes that gap.

### G.2 Direction at the molecule, agent by agent

| Agent | Approval | Direction on BubR1 abundance | Verdict | Fate |
|---|---|---|---|---|
| **NMN** | Investigational | **increase** — measured in vivo (North 2014) | **AGREES** | Rejected: not an approved medicine (ASSUMPTION-DRUG-03), and its effect on aneuploidy and cancer risk has never been assessed (blocking, ADR 0011). Scores 0.605, the highest in the catalogue, and is still rejected. |
| **Nicotinamide riboside** | Food ingredient | **unknown** — raises human NAD⁺ (Trammell et al., *Nat Commun* 2016;7:12948, PMID 27721479; Martens et al., *Nat Commun* 2018;9:1286, PMID 29599478) but BubR1 never measured | cannot determine | Rejected: not an approved medicine. |
| **Nicotinamide** | Approved (pellagra) | **context_dependent** — NAD⁺ precursor *and* a product inhibitor of sirtuins, measured against human SIRT2 itself alongside SIRT1 (Marcotte et al., *Anal Biochem* 2004;332:90-9, PMID 15301953), and separately against yeast Sir2 and human SIRT1 (Bitterman et al., *J Biol Chem* 2002;277:45099-107, PMID 12297502); the pocket conserved across the family is why (Avalos, Bever & Wolberger, *Mol Cell* 2005;17:855-68, PMID 15780941). Potency is assay- and substrate-dependent, so there is no single number — 12–600 µM across nine substrate conditions on one enzyme preparation (Zessin et al., *Int J Mol Sci* 2023, PMID 37108579) and 138 µM independently (Chiang and Lin, *Org Biomol Chem* 2016, PMID 26796034). It feeds and inhibits the same axis, and the *net* effect is what is unmeasured | **cannot determine** | **Rejected**, score 0.275, on `oncogenic_risk`: no micronucleus or aneuploidy measurement of this agent exists in any species, and the skin-cancer trial that had filled that cell does not bear on it (§5.1). Kept in the experiment as the **high-dose negative control** (G.3). |
| **Nicotinic acid** | Approved (dyslipidaemia) | **unknown** — raises human blood NAD⁺ up to 8-fold at 750–1,000 mg/day (Pirinen et al., *Cell Metab* 2020;31:1078-1090, PMID 32386566). The accompanying functional benefit in mitochondrial myopathy will not carry weight: the trial is open-label with no placebo arm, strength and mitochondrial biogenesis rose "in all subjects" — the matched controls included — so no disease-specific benefit can be read out of it, and some patients showed an anaemia tendency. "No sirtuin-inhibiting liability" is an **inference** from nicotinic acid not being the sirtuin reaction product; no direct dose-response measurement was located, and Avalos 2005 engineered a Sir2 variant that *is* inhibited by nicotinic acid, so the claim is about wild-type enzymes. BubR1 never measured | cannot determine | **Accepted, rank 1, score 0.500**, with the unsigned-direction flag attached. The oncogenic-risk cell is assessed-and-negative in adults (§5.2), which is a closed gate and not a clean bill of health. |

### G.3 Which of the two is the probe, and why the other one stays in the experiment

The two rows above changed places, and the reason is worth separating from the reason the
*pipeline* changed places. The pipeline moved because two `worsens_cin` cells were
corrected (§5). The scientific reading points the same way for a different reason, and it
would hold even if the scores had not moved.

**Nicotinic acid is the mechanistic probe, because it moves one arm.** It enters the NAD⁺
pool through the Preiss-Handler route, it is not a sirtuin reaction product, and no
measurement of it inhibiting a wild-type sirtuin was located. A result from it is
therefore readable in both directions: BubR1 rises and the missegregation rate falls, and
the NAD⁺ → SIRT2 → K668 axis is supported; nothing moves, and the axis is in trouble.

**Nicotinamide cannot produce a readable result at all**, and that — not its safety file —
is the disqualifying argument against it as the first experiment. It raises NAD⁺ *and*
inhibits the sirtuin that spends it, in the same cell at the same time. If BubR1 goes up,
supply won or something off-axis did, and the experiment cannot say which. If BubR1 goes
down or does not move, the mechanism is wrong or the enzyme it runs on was inhibited, and
the experiment cannot say which. Both outcomes leave the hypothesis exactly where it
started. A probe with that property is not a weak probe; it is not a probe.

**So it becomes the negative control, and it is a good one.** The mechanism makes a signed,
falsifiable prediction about it: at an exposure above the measured SIRT2 IC₅₀ range,
nicotinamide should *lower* BubR1. That prediction has already been met once, in the
founding paper for this axis — 20 mM nicotinamide for 16 h in HeLa cells reduced BubR1
protein, with no change in cell-cycle distribution and no change in *BUB1B* mRNA (North et
al., *EMBO J* 2014, PMID 24825348, Figure 1A and 1B and Methods; sirtinol at 25 µM did the
same, and nicotinamide lowered BubR1 in MEFs in the same paper's supplement). Read that
carefully in both directions. It is the sharpest single fact against nicotinamide as a
therapy and it is a *tool-compound* concentration: 20 mM is roughly 14–25× the 790–1,400 µM
peak plasma reached at 60–80 mg/kg oral in patients (Hoskin et al., *Int J Radiat Oncol
Biol Phys* 1995, PMID 7607932; Bussink et al., *Radiother Oncol* 2002, PMID 12142092). It
fixes the *sign* of the inhibitory arm. It does not fix the dose at which that arm wins,
and this document does not claim it does.

The dose comparison is the part that keeps nicotinamide dangerous rather than merely
uninterpretable. Peak plasma at 60–80 mg/kg oral sits **above every measured SIRT2 IC₅₀**,
including the 600 µM top of Zessin's nine-condition range. Taking nicotinamide by the
niacin route instead inverts the exposure: after 2 g of extended-release niacin the parent
acid peaks at about 76 µM and the nicotinamide metabolite below that (Menon et al., *Int J
Clin Pharmacol Ther* 2007, PMID 17725178) — roughly an order of magnitude under the
inhibitory range, not zero. Niacin is the cleaner arm on this axis, not the clean one.

Four things stop this from being a settled argument, and all four are real. Plasma is not
the compartment the enzyme sees, and chronic oral nicotinamide did not raise steady-state
hepatic nicotinamide in mice at all (Mitchell et al., *Cell Metab* 2018, PMID 29514072) —
it went into methylated and oxidised catabolites; no human study of intracellular free
nicotinamide after dosing was located. Peaks are brief. The IC₅₀ is a property of the assay
more than of the enzyme — 12 µM to 600 µM on one enzyme preparation, and BubR1-K668 has
never been the substrate in such an assay, so the number that matters has not been
measured. And the one in vivo readout of a SIRT2 substrate under chronic high-dose oral
nicotinamide went the *favourable* way: tubulin acetylation fell, meaning SIRT2 output rose
(Mitchell et al., *Cell Metab* 2018, PMID 29514072). The honest verdict on the inhibitory
arm is *unresolved at therapeutic doses and unresolvable by nicotinamide itself*, which is
precisely why the first experiment should be run with the other molecule.

## Appendix H — Where the automation boundary falls

| Stage | Who did it | Deterministic? |
|---|---|---|
| Which agents and classes enter the catalogue at all | LLM-assisted literature curation (Claude) | No |
| Every `observed_direction` sign, and the citation behind it | LLM-assisted literature curation | No |
| Every mechanism link, its grade and its uncertainty cell | LLM-assisted literature curation | No |
| Choice of target node; the rule that unassessed CIN risk blocks; the weights | Human-reviewable decision records (ADR 0025, ADR 0011, ASSUMPTION-SCORING-01) | No — but written down and testable |
| Target-on-chain check, required direction per node, tri-state agreement, gates, score, rank, rejection-reason priority, the 89 evidence rows | `src/mva/interventions/` — plain Python, no model inference (ADR 0003) | **Yes.** Byte-identical under a fixed clock (GP-30); `test_triage_is_deterministic` hashes the output |

## Appendix I — Three corrections, and one correction that was itself wrong

The claims ledger's most instructive catches, referred to from §0.

**A sirtuin claim cited to the wrong enzyme.** The nicotinamide/sirtuin inhibition claim
was cited to Bitterman et al., *J Biol Chem* 2002;277:45099-107 (PMID 12297502), which
assayed yeast Sir2 and human SIRT1 — not the SIRT2 this document's mechanism runs on.

**The sequel to that correction was also wrong.** Having caught the mis-citation, the
audit went on to assert that nobody had measured nicotinamide against SIRT2, on the
strength of not having found such a measurement. That is GP-14 violated in the act of
correcting a citation, and a later external review found the measurement in a single
search: Marcotte et al., *Anal Biochem* 2004;332:90-9 (PMID 15301953) assayed nicotinamide
as a competitive inhibitor of human SIRT1 and SIRT2. The verdict on nicotinamide did not
change; the reason for it did, and the new reason is better, because the inhibition is
measured rather than inferred from the conserved pocket (Avalos, Bever & Wolberger, *Mol
Cell* 2005;17:855-68, PMID 15780941).

**A source cited for the opposite of what its title says.** A claim that aneuploid cells
*depend on* autophagic degradation was cited to a paper titled "Aneuploidy-induced
cellular stresses limit autophagic degradation" (Santaguida et al., *Genes Dev*
2015;29:2010-21, PMID 26404941). The direction reverses, and it happens to reverse in
favour of the argument that inhibiting autophagy is selectively lethal to aneuploid cells
(§3), which is why it survived unnoticed.

**A guideline that never mentions this syndrome.** A tumour-surveillance guideline was
cited as the operative protocol for *BUB1B*-MVA; its text covers Beckwith-Wiedemann,
WAGR, *WT1*-related disorders and nine other predisposition syndromes, and mentions MVA,
*BUB1B* and PCS syndrome nowhere. It was withdrawn and replaced by Scott 2006's general
>5% Wilms-risk rule, argued from this syndrome's own published tumour rates
(Appendix E.1).

---

*Nothing in this document is medical advice. Every output is a research hypothesis
requiring pre-clinical validation. No patient data was accessed in producing it.*
