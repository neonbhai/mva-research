# R2: Nicotinamide vs SIRT2 — does product inhibition sink it as a mechanistic probe?

Scope: public literature only. Every factual claim carries author, year, journal, PMID. Anything I could not resolve is marked "unresolved", not estimated. Section 8 lists every search run, including the ones that found nothing.

Question under test. Track 2 proposes raising BubR1 protein by NAD+ -> SIRT2 -> deacetylation of BubR1-K668. The two approved NAD+ precursors are nicotinamide and nicotinic acid. Nicotinamide is also the reaction product of every sirtuin and a feedback inhibitor of them. A reviewer says ranking nicotinamide first is indefensible. This document collects the measured numbers.

---

## 1. Verdict

Yes at the doses actually used in humans: peak plasma nicotinamide after 60-80 mg/kg oral is 790-1400 uM (Hoskin PJ et al., *Int J Radiat Oncol Biol Phys* 1995, PMID 7607932; Bussink J et al., *Radiother Oncol* 2002, PMID 12142092), which sits above every measured SIRT2 nicotinamide IC50 (12-600 uM across nine substrate conditions in Zessin M et al., *Int J Mol Sci* 2023, PMID 37108579; 138 uM in Chiang YL et al., *Org Biomol Chem* 2016, PMID 26796034), and the founding paper for our own mechanism used 20 mM nicotinamide as a sirtuin inhibitor and saw BubR1 protein fall (North BJ et al., *EMBO J* 2014, PMID 24825348). But the inference is not secure, because the one in vivo experiment with a SIRT2-specific readout under chronic high-dose oral nicotinamide found tubulin acetylation went **down** — SIRT2 output up, not down — while steady-state tissue nicotinamide did not rise at all (Mitchell SJ et al., *Cell Metab* 2018, PMID 29514072). So: **yes** on plasma-versus-IC50 and yes as a risk that cannot be dismissed; **unknown** at the intracellular level, unknown at doses below ~25 mg/kg, and unknown for BubR1-K668 specifically, which has never been used as the substrate in a nicotinamide inhibition assay.

**Consequence for the ranking.** The reviewer's operational point stands even where his premise is shaky. Nicotinamide moves both arms of the mechanism at once, so no result it produces identifies the pathway. It should not be first. See section 7.

---

## 2. Measured potency of nicotinamide against SIRT2

### 2a. The anchor paper: Marcotte 2004

Marcotte PA et al., *Anal Biochem* 2004, PMID 15301953. This is the paper the pipeline cited. The abstract confirms scope but not numbers:

> "Using this method, the kinetics of the reaction of the cosubstrate nicotinamide adenine dinucleotide and the competitive inhibitor nicotinamide with SIRT1 and SIRT2 has been analyzed."

Two flags. First, the abstract calls nicotinamide a **competitive** inhibitor. The mechanistic literature calls it **non-competitive** with respect to NAD+ — it acts by base exchange onto a covalent ADP-ribosyl-enzyme-acetyl-peptide intermediate, regenerating NAD+ (Jackson MD et al., *J Biol Chem* 2003, PMID 14522996; Sauve AA et al., *Biochemistry* 2003, PMID 12899610; Avalos JL et al., *Mol Cell* 2005, PMID 15780941). One later group re-measured and argued the kinetics look *apparently competitive* with NAD+ after all (Guan X et al., *PLoS One* 2014, PMID 25221980). So the inhibition mode itself is contested, which matters: a non-competitive inhibitor is not outrun by raising NAD+.

Second, Marcotte's own assay uses "nicotinamide at millimolar concentrations" as the reaction stop reagent. The paper's working assumption is therefore that *complete* block needs mM.

**The numeric Ki/IC50 for SIRT2 in Marcotte 2004 is UNRESOLVED.** The paper is paywalled (Elsevier, *Anal Biochem* 332:90-99); no PMC deposit exists (PubMed-to-PMC link returns zero). I checked 26 papers that cite it and none quote its SIRT2 Ki. Searches listed in section 8. Do not quote a number for Marcotte until someone reads the PDF.

### 2b. Primary measurements

| Enzyme | Substrate / acyl group | Assay format | NAM potency | Source |
|---|---|---|---|---|
| human SIRT2 | C1, acetyl-Lys fluorogenic peptide, 20 µM, 100 nM enzyme | discontinuous, trypsin developer | IC50 = **12 ± 1 µM** | Zessin M et al., *Int J Mol Sci* 2023, PMID 37108579 (Table 2) |
| human SIRT2 | Mcm1, myristoyl-Lys TNFα peptide, 1 µM | continuous fluorescence | IC50 = **43 ± 6 µM** | Zessin M et al., *Int J Mol Sci* 2023, PMID 37108579 |
| human SIRT2 | Mcm1 + 0.1% Tween20, 1 µM | continuous fluorescence | IC50 = **45 ± 3 µM** | Zessin M et al., *Int J Mol Sci* 2023, PMID 37108579 |
| human SIRT2 | C2, myristoyl-Lys, 1 µM | discontinuous, trypsin developer | IC50 = **42 ± 6 µM** | Zessin M et al., *Int J Mol Sci* 2023, PMID 37108579 |
| human SIRT2 | F4, myristoyl-Lys fluorescein peptide, 1 µM | continuous fluorescence | IC50 = **100 ± 5 µM** | Zessin M et al., *Int J Mol Sci* 2023, PMID 37108579 |
| human SIRT2 | S1, 1 µM | continuous | IC50 = **120 ± 10 µM** | Zessin M et al., *Int J Mol Sci* 2023, PMID 37108579 |
| human SIRT2 | F4, 40 nM (near Km) | continuous fluorescence | IC50 = **211 ± 7 µM** | Zessin M et al., *Int J Mol Sci* 2023, PMID 37108579 |
| human SIRT2 | S2, thioacetyl-aminoundecanoyl-Lys, 1 µM | continuous | IC50 = **340 ± 7 µM** | Zessin M et al., *Int J Mol Sci* 2023, PMID 37108579 |
| human SIRT2 | S2, 10 nM | continuous | IC50 = **600 ± 70 µM** | Zessin M et al., *Int J Mol Sci* 2023, PMID 37108579 |
| human SIRT2 | myristoyl-Lys aminocoumarin peptide (substrate 4), 5 µM, 1 mM NAD+, 1 µM enzyme, 37 °C | two-step fluorogenic, trypsin developer | IC50 = **138 µM** | Chiang YL et al., *Org Biomol Chem* 2016, PMID 26796034 |
| human SIRT1 | same | same | IC50 = **79 µM** | Chiang YL et al., *Org Biomol Chem* 2016, PMID 26796034 |
| human SIRT3 | same | same | IC50 = **84 µM** | Chiang YL et al., *Org Biomol Chem* 2016, PMID 26796034 |
| human SIRT2 | acetyl / hexanoyl / decanoyl / myristoyl H3K9 peptides, saturating NAD+, 37 °C, pH 7.5 | steady-state, quench, HPLC | IC50 varies with acyl chain; **highest for acetyl, falls as chain lengthens**. Numeric values are in Figure 4 only — UNRESOLVED (see 2c) | Feldman JL et al., *Biochemistry* 2015, PMID 25897714 |
| human SIRT1 | ACS2-derived acetyl peptide | MS-based deacetylation assay | IC50 = **62 ± 4 µM** | Fischer F et al., *PLoS One* 2012, PMID 23028781 |
| human SIRT3 | ACS2-Lys642 acetyl peptide | MS-based | IC50 = **43 ± 3 µM** | Fischer F et al., *PLoS One* 2012, PMID 23028781 |
| human SIRT5 | CPS1-Lys527 acetyl peptide | MS-based | IC50 = **1.6 ± 0.3 mM** | Fischer F et al., *PLoS One* 2012, PMID 23028781 |
| human SIRT5 | CPS1-Lys527 succinyl peptide | MS-based | IC50 = **21 ± 4 µM** | Fischer F et al., *PLoS One* 2012, PMID 23028781 |
| yeast Hst2 | acetyl peptide | MS-based | IC50 = **49 ± 2 µM** | Fischer F et al., *PLoS One* 2012, PMID 23028781 |
| Sir2Tm (bacterial) | acetyl peptide | MS-based | IC50 = **39 ± 5 µM** | Fischer F et al., *PLoS One* 2012, PMID 23028781 |
| human SIRT3 | acetyl peptide, 1 mM NAD+, 90 min, 37 °C | fluorogenic | IC50 = **36.7 ± 1.3 µM** | Guan X et al., *PLoS One* 2014, PMID 25221980 |
| human SIRT1 | same | fluorogenic | IC50 = **68.1 ± 1.8 µM** | Guan X et al., *PLoS One* 2014, PMID 25221980 |
| yeast Sir2 / human SIRT1 | acetyl-p53 peptide | in vitro deacetylation | IC50 **< 50 µM**, non-competitive | Bitterman KJ et al., *J Biol Chem* 2002, PMID 12297502 |

### 2c. How much does it vary?

A lot. This is the single most important number in the document.

Within **one paper, one lab, one recombinant human SIRT2 preparation**, the nicotinamide IC50 spans **12 µM to 600 µM — a 50-fold range** — driven only by which peptide substrate is used and at what concentration (Zessin M et al., *Int J Mol Sci* 2023, PMID 37108579, Table 2). The authors state the direction explicitly: inhibitor potency for lysine-channel binders rose as substrate concentration fell, whereas

> "In contrast, NAM inhibition showed the opposite effect, specifically for substrates F4 and S2".

That is, nicotinamide gets **weaker**, not stronger, as you approach physiological (sub-Km) substrate occupancy on those substrates. The most potent value in the whole table, 12 µM, comes from the acetylated substrate C1 read out in a **discontinuous trypsin-developer format** — which is the same class of assay Marcotte 2004 introduced. That is a warning: the assay format that produced the low number is the one most likely to have been used for the number our pipeline inherited.

Feldman JL et al., *Biochemistry* 2015, PMID 25897714 independently show the acyl-group dependence and give the mechanism:

> "We found that the IC 50 of nicotinamide varied with the Sirtuin and acyl-substrate analyzed... For SIRT1 and SIRT2, the IC 50 for nicotinamide is the highest for the acetylated substrate and decreases as the acyl chain length increases"

and

> "The slower conversion of the long-chain O-alkylimidate intermediates (k9) relative to deacetylation likely increases the sensitivity of SIRT2 to inhibition by nicotinamide"

Note this contradicts the Zessin ordering for acetyl vs myristoyl (Feldman: acetyl least inhibited; Zessin: acetyl-C1 most inhibited at 12 µM). The two used different peptides, different readouts and different substrate concentrations. Taken together they say the same operational thing: **the nicotinamide IC50 for SIRT2 is a property of the assay, not of the enzyme.** A single quoted IC50 is not a usable input to a ranking algorithm.

Fischer F et al., *PLoS One* 2012, PMID 23028781 make the substrate point even harder within one enzyme: SIRT5 deacetylation IC50 = 1.6 mM but SIRT5 desuccinylation IC50 = 21 µM — a **76-fold shift from changing only the acyl group on the same lysine of the same peptide**.

A second independent lab, using a myristoyl-aminocoumarin substrate at 5 µM with 1 mM NAD+, reports SIRT2 IC50 = 138 µM and notes it was "comparable to those previously reported" (Chiang YL et al., *Org Biomol Chem* 2016, PMID 26796034). So the two independent, clearly documented SIRT2 measurements I could reach are 138 µM (Chiang 2016) and a nine-condition spread of 12-600 µM (Zessin 2023).

The widely repeated round number "SIRT2 IC50 = 100 µM" traces to a review sentence in Guan X et al., *PLoS One* 2014, PMID 25221980:

> "The IC 50 values for nicotinamide inhibition of bacterial Sir2, yeast Sir2, mouse Sir2, SIRT1, SIRT2, SIRT3, and SIRT5 were measured to be 26, 120, 160, 50, 100, 36.7 µM, and 1.6 mM respectively"

whose citation for the SIRT2 entry is Tervo AJ et al., *J Med Chem* 2004, PMID 15566299 — a virtual-screening paper in which nicotinamide is a reference compound, not the subject. I could not open that paper (paywalled, ACS). **Treat "SIRT2 = 100 µM" as a secondary, uncorroborated figure**, not as the primary measurement. Coincidentally Zessin's F4-at-1-µM condition also gives 100 ± 5 µM, but that is one condition out of nine.

### 2d. What is NOT measured

I did not find any measurement of nicotinamide inhibition of SIRT2 **on BubR1-K668 or on any BubR1-derived peptide**. Every SIRT2 nicotinamide number above uses histone H3, p53, TNFα or synthetic model peptides. Given the 50-fold substrate dependence documented in 2c, the SIRT2/BubR1-K668 IC50 is unknown. Searches that failed are in section 8.

---

## 3. The crossover: at what nicotinamide concentration does inhibition beat supply?

### 3a. Short answer

There is no published crossover concentration for SIRT2. I found no study that titrates nicotinamide against a SIRT2-dependent output and reports the concentration at which net deacetylation turns from up to down. What exists is a set of point observations that go in **both** directions, and a review that says the question is unsettled. Searches that failed are in section 8.

### 3b. The best in vivo data point, and it cuts against the reviewer

Mitchell SJ et al., *Cell Metab* 2018, PMID 29514072. One-year-old male C57BL/6J mice fed standard diet or high-fat diet supplemented with nicotinamide at 0.5 or 1.0 g/kg of diet for life. This is the most directly relevant experiment I found: chronic oral nicotinamide, mammalian tissue, sirtuin-substrate readouts including a SIRT2 substrate.

Findings, in the authors' words:

> "Though neither hepatic NAD + nor NADP + were boosted by NAM, acetylation of some SIRT1 targets was enhanced by NAM supplementation in a diet and NAM dose-dependent manner."

> "Immunoblotting of total liver lysates with acetyl-lysine antibody revealed a clear increase in global protein acetylation in response to NAM treatment, especially in HFD-fed mice ... suggesting overall sirtuin inhibition."

So far this supports the reviewer. But the SIRT2-specific readout goes the other way:

> "Unlike SIRT1 targets, acetylation of tubulin (SIRT2 target) was significantly reduced in HFDH liver, whereas SOD2 acetylation (SIRT3 target) was unresponsive to NAM"

Acetyl-α-tubulin K40 is the canonical SIRT2 substrate (North BJ et al., *Mol Cell* 2003, PMID 12620231). Reduced tubulin acetylation at the **high** nicotinamide dose means SIRT2-directed deacetylation went **up**, not down, under chronic high-dose oral nicotinamide in vivo. Three sirtuins, three different directions, in the same livers. This is the closest thing in the literature to the experiment our reviewer says cannot work, and for SIRT2 specifically it worked in the direction we want.

Caveats I will not paper over: this is mouse liver, not human tissue; it is one figure panel (S3F/S3G) in a supplement; the effect reached significance only in the high-fat-diet high-NAM group; and tubulin acetylation is also controlled by HDAC6 and by αTAT1, so it is not a clean SIRT2 reporter.

### 3c. Why a fixed crossover may not exist: nicotinamide does not accumulate

Also from Mitchell SJ et al., *Cell Metab* 2018, PMID 29514072:

> "Though NAM supplementation did not elevate the steady-state concentration of hepatic NAM, it increased in a dose-dependent manner the levels of methylated and oxidized NAM metabolites."

Lifelong oral nicotinamide did **not** raise steady-state tissue nicotinamide. It was disposed of — into NAD+ salvage, into N1-methylnicotinamide, into 2-PY/4-PY oxidation products. The same paper reports that NAMPT protein fell in NAM-treated liver, i.e. the salvage step feeding NAD+ is itself downregulated by the substrate. So the intracellular nicotinamide concentration is a *buffered* variable, not one you set by dosing. This is the strongest single reason to think that the naive "plasma Cmax vs in vitro IC50" comparison in section 4 overstates the risk.

### 3d. The review that examined this question directly

Hwang ES et al., *Cell Mol Life Sci* 2017, PMID 28417163. Title: "Nicotinamide is an inhibitor of SIRT1 in vitro, but can be a stimulator in cells." From the abstract:

> "once administered to cells, NAM is rapidly converted to NAD+ and, therefore, the cellular concentration of NAM decreases rapidly while that of NAD+ increases. The result would be an inhibition of SIRT1 for a limited duration, followed by an increase in the activity. This possibility raises a concern on the validity of the interpretation of the results in the studies that use NAM as a SIRT1 inhibitor."

> "we reviewed published literature in which treatment with NAM was used to inhibit SIRT1 and found that the expected inhibitory effect of NAM was either unreliable or muted in many cases. In addition, studies demonstrated NAM administration stimulates SIRT1 activity and improves the functions of cells and organs."

> "We conclude that NAM treatment can hypothetically be stimulatory to SIRT1."

Note the hedge in the last sentence: "hypothetically". This review argues the crossover is time-dependent (early inhibition, later stimulation) rather than a fixed concentration, and its recommendation for anyone using NAM as an inhibitor is to assay within 1-2 h of treatment. The full text is behind a Springer paywall; I read the abstract from PubMed and the body via the PMC HTML view. Numbers I take from the body — cellular NAM 0.7-13.7 µM in human fibroblasts, cellular NAD+ 360-540 µM, NAM at 5-20 mM raising NAD+ by 50-100% over 8-12 h — are **second-hand within a review** and I did not resolve their primary sources. Marked accordingly.

### 3e. What "crossover" even means here is contested

The two mechanistic pictures give different answers:

- Non-competitive / base-exchange model (Jackson MD et al., *J Biol Chem* 2003, PMID 14522996; Sauve AA et al., *Biochemistry* 2003, PMID 12899610; Avalos JL et al., *Mol Cell* 2005, PMID 15780941). Nicotinamide intercepts a covalent intermediate. Raising NAD+ does **not** relieve it. Under this model there is a real crossover and it is unfavourable.
- Apparent-competitive model (Guan X et al., *PLoS One* 2014, PMID 25221980): "we demonstrate that NAM inhibition of SIRT3 involves apparent competition between the inhibitor and the enzyme cofactor NAD+, contrary to the traditional characterization of base exchange as noncompetitive inhibition." Under this model the NAD+/NAM ratio is what matters, and a precursor that raises NAD+ more than it raises free NAM comes out ahead.

Sauve AA et al., *Biochemistry* 2003, PMID 12899610 also quantify how much of the intermediate nicotinamide can claw back: "Saturation of the nicotinamide site for mouse, yeast, and bacterial Sir2s causes 95, 65, and 21% of the intermediate, respectively, to return to acetylated protein." The species spread is large and **no value for human SIRT2 is given there**. I did not find this "nicotinamide switching" fraction measured for human SIRT2.

### 3f. A second in vitro result that widens the range further

Rymarchyk S et al., *Biomolecules* 2021, PMID 33670751 (Cen and Sinclair labs) measured Ki(NAM) for human SIRT1 against four physiological peptide substrates by HPLC, at 500 µM NAD+, 37 °C:

| SIRT1 substrate | Ki (NAM) |
|---|---|
| p53-K382Ac | 294.4 ± 44.6 µM |
| H3K9Ac | 247.9 ± 35.1 µM |
| H4K16Ac | 356.6 ± 38.2 µM |
| p300-K1024Ac | 336.1 ± 68.7 µM, and inhibition is **incomplete** |

These are 4-7x higher than the 50-68 µM commonly quoted for SIRT1. Two consequences. First, the popular SIRT1 number is soft too. Second, and more useful: for p300-K1024Ac, "the deacetylation of p300K1024Ac cannot be completely inhibited at physiological NAM concentrations, in sharp contrast to the behaviors of other peptides." Nicotinamide can be a **partial** inhibitor, and which substrate you pick decides that.

Rymarchyk S et al., *Biomolecules* 2021, PMID 33670751 also state the position for SIRT2 explicitly, citing earlier work: "For mouse Sir2α and human SIRT2, complete inhibition can be achieved with increasing concentrations of NAM". So SIRT2 is not one of the escape cases — at high enough nicotinamide it does shut off completely. The open question is only where "high enough" sits for the BubR1-K668 substrate, and that has not been measured.

---

## 4. Achievable human concentrations vs inhibitory concentrations

### 4a. Nicotinamide, measured in humans

| Dose / regimen | Population | Measured plasma nicotinamide | Source |
|---|---|---|---|
| 80 mg/kg/day oral, during radiotherapy | adults, head/neck and rectal cancer, n=11, HPLC | peak **0.4-1.4 µmol/mL (400-1400 µM)**, mean 0.9 ± 0.3 µmol/mL; levels >= 0.7 µmol/mL held for a mean 2.8 h | Hoskin PJ et al., *Int J Radiat Oncol Biol Phys* 1995, PMID 7607932 |
| 60 mg/kg oral (ARCON) | adults, n=22 and n=10 | mean Cmax **793 nmol/mL (793 µM)** without domperidone, 776 nmol/mL with; tmax ~0.5 h | Bussink J et al., *Radiother Oncol* 2002, PMID 12142092 |
| 6 g single oral | adults with cancer, n=12, escalation 3-10 g | mean peak **156.4 ± 33.6 µg/mL = 1281 ± 275 µM** (my conversion, MW 122.12); peak by 30 min; dose-peak relationship linear to 6 g (r=0.91) | Dragovic J et al., *Radiother Oncol* 1995, PMID 8532910 |
| 2.5 vs 25 mg/kg, immediate vs sustained release | healthy young adult men | **kinetics are non-linear**: a 10-fold dose increase produced a **62-fold AUC increase**. Cmax values are in the paper's figures/tables — I could not open the full text, so specific Cmax at 25 mg/kg is UNRESOLVED | Petley A et al., *Diabetes* 1995, PMID 7859933 |

Conversion used throughout: nicotinamide MW 122.12 g/mol, so 1 µg/mL = 8.19 µM.

### 4b. Paediatric

The relevant paediatric precedent is the type 1 diabetes prevention programme: oral nicotinamide **1.2 g/m2/day**, given for 5 years to relatives from age 3 upward, modified release in ENDIT (Gale EA et al., *Lancet* 2004, PMID 15043959) and standard release in a parallel Cuban trial (Cabrera-Rode E et al., *Autoimmunity* 2006, PMID 16891222). Safety at this dose in children is documented: growth was unaffected, serious adverse events did not differ from placebo (Gale EA et al., *Lancet* 2004, PMID 15043959), and the safety review found the therapeutic index wide with reversible hepatotoxicity only at very high doses, flagging adult doses above 3 g/day as the caution threshold (Knip M et al., *Diabetologia* 2000, PMID 11126400).

**I did not find any measured plasma nicotinamide concentration in children at 1.2 g/m2/day.** Searches listed in section 8. What can be said without measurement: 1.2 g/m2 for a 5-year-old of roughly 0.75 m2 and 18 kg is about 50 mg/kg/day, which sits between the 60 mg/kg and the sub-radiosensitiser range, so plasma in the high-hundreds of µM is a reasonable expectation. That is an inference from adult mg/kg data, not a measurement, and it is exactly the kind of number our document should not assert.

### 4c. Nicotinic acid, measured in humans

| Dose / regimen | Population | Measured plasma | Source |
|---|---|---|---|
| 2000 mg extended-release niacin, single dose after low-fat snack | 12 healthy men | niacin Cmax **9.3 µg/mL = ~75.5 µM** (MW 123.11), AUC(0-t) 26.2 µg·h/mL; tmax 4.6 h; t1/2 0.9 h. Niacin was "the highest of all analytes measured", so the **nicotinamide metabolite peaked below ~76 µM**; its tmax was 8.6 h and t1/2 4.3 h | Menon RM et al., *Int J Clin Pharmacol Ther* 2007, PMID 17725178 |
| 500 / 1000 / 1500 mg/day extended-release niacin | 10 dialysis and 8 CKD patients | nicotinic acid, nicotinamide and nicotinuric acid all measured over 24 h; the abstract reports only that CKD levels were higher than dialysis and "not higher than reported in patients without renal impairment". Numeric Cmax values are in the full text — UNRESOLVED | Reiche I et al., *Nephrol Dial Transplant* 2011, PMID 20562093 |

### 4d. The comparison

Put the two columns next to each other.

- SIRT2 nicotinamide IC50, measured: **12-600 µM** across nine substrate conditions (Zessin M et al., *Int J Mol Sci* 2023, PMID 37108579) and **138 µM** from an independent lab (Chiang YL et al., *Org Biomol Chem* 2016, PMID 26796034), plus an uncorroborated secondary value of 100 µM traced to Tervo AJ et al., *J Med Chem* 2004, PMID 15566299.
- Plasma nicotinamide after therapeutic nicotinamide dosing: **~790-1400 µM at peak** (Hoskin 1995; Bussink 2002; Dragovic 1995).
- Plasma nicotinamide after therapeutic nicotinic acid dosing: **below ~76 µM** (Menon 2007).

At nicotinamide doses in the 60-80 mg/kg range, peak plasma nicotinamide exceeds the top of the measured SIRT2 IC50 range by roughly 1.5-3x and the bottom by up to 100x. That is not a marginal overlap. On the plasma-vs-IC50 comparison alone, the reviewer is right.

Four things weaken that conclusion, and all four are real:

1. **Plasma is not the compartment.** The enzyme sees intracellular free nicotinamide. Chronic oral nicotinamide in mice **did not raise steady-state hepatic nicotinamide at all** (Mitchell SJ et al., *Cell Metab* 2018, PMID 29514072); the dose went into methylated and oxidised catabolites instead. I did not find any human study measuring intracellular or tissue free nicotinamide after dosing. That is the single biggest gap in this analysis.
2. **Peaks are brief.** tmax is 0.5-2 h and levels above 700 µM were held for a mean of 2.8 h (Hoskin PJ et al., *Int J Radiat Oncol Biol Phys* 1995, PMID 7607932). Between peaks, plasma falls. BubR1 protein stability is a slow variable; a few hours per day of SIRT2 suppression is not obviously equivalent to continuous suppression.
3. **The IC50 is not a fixed property.** Section 2c: 50-fold spread on substrate alone, and the BubR1-K668 value is unmeasured.
4. **The one in vivo SIRT2 readout under chronic high-dose nicotinamide went the favourable way** (tubulin acetylation fell; Mitchell SJ et al., *Cell Metab* 2018, PMID 29514072).

Dose is the lever. At 500 mg twice daily — the dose used in the skin-cancer chemoprevention trial and the common supplement dose — the exposure is roughly a tenth of the radiosensitiser regimens, and I found no measured plasma level at that dose. **UNRESOLVED, and worth resolving**, because a low-dose nicotinamide arm may sit below the inhibitory range while a high-dose arm certainly does not.

---

## 5. Nicotinic acid: liability or clean?

Our document asserted, without a citation, that nicotinic acid is free of the product-inhibition problem. The assertion is **substantially correct but needs three qualifications**, and here are the citations it was missing.

### 5a. Evidence that nicotinic acid does not inhibit sirtuins

**Direct, and the best single citation.** Gallo CM et al., *Mol Cell Biol* 2004, PMID 14729974. Pnc1 is a nicotinamidase that converts nicotinamide to nicotinic acid. From the abstract:

> "Pnc1 is a nicotinamidase that converts NAM to nicotinic acid. Here we show that recombinant Pnc1 stimulates Sir2 HDAC activity in vitro by preventing the accumulation of NAM produced by Sir2."

Converting nicotinamide into nicotinic acid **relieves** sirtuin inhibition. If nicotinic acid were itself an inhibitor at the concentrations generated, the coupled reaction would not have stimulated Sir2. This is the cleanest available demonstration that the acid is not a substitute inhibitor for the amide. Caveat: yeast Sir2, in vitro, and the paper shows relief of inhibition rather than a direct nicotinic-acid dose-response against a sirtuin.

**Structural / mechanistic.** Avalos JL et al., *Mol Cell* 2005, PMID 15780941:

> "Based on our structures, we engineered a mutant that deacetylates peptides by using nicotinic acid adenine dinucleotide (NAAD) as a cosubstrate and is inhibited by nicotinic acid."

Making a sirtuin nicotinic-acid-inhibited required protein engineering. That is only a sensible experiment if the wild-type enzyme is not inhibited by nicotinic acid. This is a strong inference, not an explicit statement — I did not open the full text, and the paper may state it outright somewhere I could not see. Marked as inference.

**Chemical selectivity of the C-pocket.** Schmidt MT et al., *J Biol Chem* 2004, PMID 15269219:

> "Sir2 enzymes showed an exquisite selectivity for the nicotinamide base coenzyme, with the most dramatic losses in binding affinity/reactivity resulting from relatively minor changes in the nicotinamide ring, either by reduction, as in NADH, or by converting the amide to its acid analogue."

The nicotinamide C-pocket reads the carboxamide. Swapping it for a carboxylate loses most of the binding. This paper also tested the salvage-pathway metabolites broadly and concluded "Nicotinamide was the most potent inhibitor examined." Numeric values per compound are in the full text, which I could not open (JBC returned HTTP 403). Whether free nicotinic acid was among the compounds tested and what its IC50 was: **UNRESOLVED**.

**Indirect, assay-design evidence.** Smith BC et al., *Anal Biochem* 2009, PMID 19615966. Their continuous sirtuin assay works by hydrolysing the nicotinamide product to nicotinic acid and ammonia with nicotinamidase, then reading the ammonia. The assay accumulates nicotinic acid in the cuvette while measuring undistorted steady-state sirtuin kinetics that "agreed very favorably" with HPLC and charcoal-binding gold standards. A working assay of this design is only possible if the nicotinic acid it generates does not inhibit the sirtuin. Again inference, not a stated result.

### 5b. The three qualifications

1. **Nicotinic acid is metabolised to nicotinamide in vivo.** After 2000 mg extended-release niacin, plasma nicotinamide is a measured metabolite with tmax 8.6 h and t1/2 4.3 h (Menon RM et al., *Int J Clin Pharmacol Ther* 2007, PMID 17725178). Niacin therefore does deliver nicotinamide to tissue — just far less of it. Since parent niacin was "the highest of all analytes measured" at 9.3 µg/mL, the nicotinamide metabolite peaked somewhere below ~76 µM. Compare against nicotinamide dosing at 790-1400 µM. Niacin is roughly an order of magnitude cleaner on this axis, not absolutely clean.

2. **Nicotinic acid needs NAPRT, and NAPRT is tissue-restricted.** Hara N et al., *J Biol Chem* 2007, PMID 17604275: NAPRT "was abundantly expressed in the small intestine, liver, and kidney in mice"; "In cells expressing endogenous NAPRT, the addition of NA but not Nam almost doubled cellular NAD contents"; and both effects "were reversed by knockdown of NAPRT expression." Many cell types are NAPRT-deficient — the phenomenon is well enough established in oncology that NAPRT1 status is used to stratify tumours for NAMPT-inhibitor plus niacin rescue (O'Brien T et al., *Neoplasia* 2013, PMID 24403854). **Whether the cell type that matters for MVA aneuploidy expresses NAPRT is a question I did not resolve and it is the decisive one for niacin.** If the target cells lack NAPRT, niacin will not raise their NAD+ at all, and its freedom from product inhibition is worth nothing.

3. **Nicotinic acid has the better feedback biology.** Also from Hara N et al., *J Biol Chem* 2007, PMID 17604275: "NAPRT, but not Nam phosphoribosyltransferase (NamPRT ...), is insensitive to the physiological concentration of NAD" — the Preiss-Handler entry step is not feedback-inhibited by NAD+, whereas the nicotinamide salvage step is. Independently, Mitchell SJ et al., *Cell Metab* 2018, PMID 29514072 found that chronic nicotinamide feeding **lowered NAMPT protein**, i.e. nicotinamide suppresses the enzyme that would convert it into NAD+. The paper's own headline is that "neither hepatic NAD+ nor NADP+ were boosted by NAM". Hara also notes the older literature: "classical studies in mammals showed that nicotinic acid (NA) is a better precursor than nicotinamide (Nam) in elevating tissue NAD levels."

**Verdict on 5.** Nicotinic acid is clean on direct sirtuin product inhibition (three converging lines, none of them a direct nicotinic-acid dose-response, so state it as "no evidence of inhibition" rather than "does not inhibit"). Its real liabilities are different ones: NAPRT tissue restriction, flushing, and the fact that it still generates some nicotinamide. Our document should now cite Gallo 2004 and Hara 2007 rather than asserting it bare.

---

## 6. NMN and NR, for comparison

### 6a. NMN is the only agent with a demonstrated BubR1 effect

North BJ et al., *EMBO J* 2014, PMID 24825348. The founding paper for our mechanism. What it actually did:

- **Sirtuin specificity.** All seven sirtuins were tested. SIRT2 and SIRT3 were the only ones that deacetylated BubR1 in vitro; only SIRT2 did so in cells, because SIRT3 is mitochondrial. So the mechanism is genuinely SIRT2-specific in the cytoplasm.
- **NMN in vivo.** NMN 500 mg/kg/day by intraperitoneal injection for 7 days raised NAD+ in heart and testes in young and aged mice, and raised BubR1 protein; in aged animals BubR1 was "brought back to the levels of the 3-month-old animals."
- **And this, which the reviewer will cite at us.** HeLa cells were treated with 20 mM nicotinamide or 25 µM sirtinol for 16 h. "Treatment with either of these inhibitors resulted in a reduction in BubR1 protein levels (Fig 1A and B)."

That last point is the strongest single argument against nicotinamide as our probe: **in the founding paper, nicotinamide was used as the sirtuin inhibitor and it lowered BubR1** — the exact opposite of our intended effect, on our exact endpoint. It has to be stated plainly in our submission.

Its force is limited by dose. 20 mM is roughly 14-25x the peak plasma concentration reached by 60-80 mg/kg oral nicotinamide in patients (790-1400 µM; Hoskin 1995, Bussink 2002) and about 16x the peak after a 6 g oral dose (Dragovic 1995). It establishes the direction of the effect at a supraphysiological concentration. It does not establish where the sign flips.

### 6b. Human exposure data for NMN and NR

| Agent | Regimen | Result | Source |
|---|---|---|---|
| NMN | single oral 100, 250, 500 mg, 10 healthy Japanese men | safe; no significant clinical or lab changes; plasma N-methyl-2-pyridone-5-carboxamide and N-methyl-4-pyridone-5-carboxamide rose dose-dependently. Note the readouts are **downstream nicotinamide catabolites**, i.e. NMN does convert onward to nicotinamide | Irie J et al., *Endocr J* 2020, PMID 31685720 |
| NR | single oral 100, 300, 1000 mg; plus a one-person pilot | dose-dependent rise in the blood NAD+ metabolome; blood NAD+ rose "as much as 2.7-fold with a single oral dose of NR"; NR has "distinct and superior pharmacokinetics to those of nicotinic acid and nicotinamide" for elevating mouse hepatic NAD+ | Trammell SA et al., *Nat Commun* 2016, PMID 27721479 |
| NR | 250 mg titrated to 1000 mg twice daily over 9 days, 8 healthy volunteers | whole-blood NAD+ increased by 100% at steady state (p=0.001); no adverse events | Airhart SE et al., *PLoS One* 2017, PMID 29211728 |

### 6c. Two things that constrain NMN and NR here

**Regulatory.** The brief asks for **approved** agents. Nicotinamide and nicotinic acid qualify: nicotinic acid is a marketed prescription drug in extended-release form, and nicotinamide has decades of trial use including a 5-year paediatric-inclusive randomised trial at 1.2 g/m2/day (Gale EA et al., *Lancet* 2004, PMID 15043959). Nicotinamide's regulatory status is worth quoting precisely, because it cuts both ways: "Nicotinamide is for regulatory purposes classed as a food additive rather than a drug and has not therefore required the formal safety evaluation normally expected of a new therapy" (Knip M et al., *Diabetologia* 2000, PMID 11126400). NR and NMN are neither approved drugs nor supported by paediatric safety data of that scale. I did not check current regulatory classifications of NR or NMN in any jurisdiction and am not asserting one.

**Route into the cell.** Nikiforov A et al., *J Biol Chem* 2011, PMID 21504897: "besides nicotinamide and nicotinic acid, only the corresponding nucleosides readily enter the cells. Nucleotides (e.g. NAD and NMN) undergo extracellular degradation resulting in the formation of permeable precursors." In that system NMN is dephosphorylated outside the cell and enters as NR. Downstream, both NR and NMN feed NAD+ and are eventually degraded back to nicotinamide — Irie's pyridone metabolites are the proof of that in humans. So **no NAD+ precursor is free of nicotinamide generation**; they differ in how much and how fast.

### 6d. Independent support for the SIRT2-BubR1 axis

Qiu D et al., *Aging Cell* 2018, PMID 29067790 report a second Sirt2-BubR1 acetylation axis in mouse oocytes, at **BubR1-K243** rather than K668, where Sirt2 depletion disrupts kinetochore-microtubule attachment and a non-acetylatable K243R mutant "markedly suppresses spindle/chromosome anomalies and thereupon lowers the incidence of aneuploid eggs." Different lysine, different system, same direction. This paper did not use nicotinamide (no occurrences in the full text), so it adds nothing to the inhibition question, but it does mean the SIRT2-BubR1-aneuploidy link is not resting on one paper.

---

## 7. If I had to pick ONE approved agent as the first mechanistic probe

**Nicotinic acid.**

The reason is not that nicotinamide would fail. It is that **nicotinamide cannot produce an interpretable result**, and interpretability is the entire job of a mechanistic probe. Nicotinamide moves both arms of the proposed mechanism in opposite directions at the same time. If BubR1 goes up, you cannot say whether NAD+ supply won or whether some SIRT2-independent effect did. If BubR1 goes down or does not move, you cannot say whether the mechanism is wrong or whether you inhibited the enzyme you were trying to drive. Either outcome leaves the hypothesis exactly where it was. On this narrow point the reviewer is right, and the right response is not to defend nicotinamide's ranking but to change it.

Nicotinic acid moves one arm only:

- It enters NAD+ synthesis by Preiss-Handler via NAPRT, which is **not feedback-inhibited by physiological NAD+**, unlike the NAMPT step that nicotinamide must use (Hara N et al., *J Biol Chem* 2007, PMID 17604275).
- Converting nicotinamide to nicotinic acid **relieves** sirtuin inhibition rather than causing it (Gallo CM et al., *Mol Cell Biol* 2004, PMID 14729974), and the sirtuin C-pocket discriminates sharply against the acid form of the pyridine ring (Schmidt MT et al., *J Biol Chem* 2004, PMID 15269219).
- Nicotinamide exposure from a niacin dose is roughly an order of magnitude lower than from dosing nicotinamide directly: below ~76 µM after 2000 mg extended-release niacin (Menon RM et al., *Int J Clin Pharmacol Ther* 2007, PMID 17725178) versus 790-1400 µM after 60-80 mg/kg nicotinamide (Hoskin PJ et al., *Int J Radiat Oncol Biol Phys* 1995, PMID 7607932; Bussink J et al., *Radiother Oncol* 2002, PMID 12142092).
- It is an actual approved prescription drug with characterised human PK, which nicotinamide (regulated as a food additive per Knip M et al., *Diabetologia* 2000, PMID 11126400) and NR/NMN are not.

**The condition that makes or breaks this choice.** NAPRT is "abundantly expressed in the small intestine, liver, and kidney" and NA raises NAD+ only "In cells expressing endogenous NAPRT" (Hara N et al., *J Biol Chem* 2007, PMID 17604275). Before running anything, confirm NAPRT expression in the cell type you intend to probe. In a NAPRT-negative cell, nicotinic acid does nothing at all and the experiment is wasted. I did not resolve NAPRT expression in the cell types relevant to MVA (lymphocytes, fibroblasts, bone marrow). **That is the single most valuable next lookup.**

**Two honest caveats.** Flushing is a real tolerability problem, especially in children, and it is the reason nicotinic acid is rarely used paediatrically. And no NAD+ precursor is nicotinamide-free downstream: niacin's own metabolism generates nicotinamide (Menon 2007), just much less of it.

**If the question is a cell-based first experiment rather than a patient-facing ranking**, the answer is different and simpler: use NMN or NR, because North BJ et al., *EMBO J* 2014, PMID 24825348 already established the readout with NMN, and "approved" is not a constraint in a dish. The approved-agent constraint is what forces the nicotinic acid answer.

**What our submission should change.** Demote nicotinamide from first. State the North 2014 20 mM result explicitly rather than letting a reviewer find it. Promote nicotinic acid with the Gallo 2004 and Hara 2007 citations that were previously missing, and attach the NAPRT precondition to it. If a nicotinamide arm is kept at all, keep it as a **deliberate negative control** at a high dose — the mechanism predicts BubR1 should fall, and North 2014 says it does.

---

## 8. Searches run, including those that found nothing

All PubMed queries via E-utilities esearch (`db=pubmed`), then efetch/esummary. Full text via efetch `db=pmc`, WebFetch, or publisher site as noted.

### 8a. Queries that found what I needed

| Query | Outcome |
|---|---|
| `SIRT2[TIAB] AND nicotinamide[TIAB] AND (IC50[TIAB] OR Ki[TIAB] OR kinetic*[TIAB])` | 18 hits; led to Zessin 2023, Chiang 2016, Feldman 2015 |
| `sirtuin[TIAB] AND nicotinamide inhibition[TIAB]` | 13 hits; led to Fischer 2012, Rymarchyk 2021 |
| `nicotinamide AND SIRT2 AND (inhibition OR inhibitor OR Ki OR IC50)` | 60 hits, relevance-sorted |
| `Bitterman AND nicotinamide AND SIRT1 AND 2002` | PMID 12297502 |
| `Jackson AND Denu AND nicotinamide inhibition transglycosidation Sir2` | PMID 14522996 |
| `Avalos JL[Author] AND nicotinamide` | PMID 15780941 |
| `Schmidt MT AND Denu AND coenzyme specificity Sir2` | PMID 15269219 |
| `"nicotinic acid" AND deacetylase AND (inhibition OR Ki OR IC50) AND NAD` | 11 hits; found Gallo 2004, PMID 14729974 |
| `Hara N AND "nicotinic acid phosphoribosyltransferase" AND NAD` | PMID 17604275 |
| `Menon RM AND niacin AND pharmacokinetics` | PMID 17725178 |
| `nicotinamide AND (radiosensitiz* OR ARCON) AND plasma` | 11 hits; Hoskin 1995, Bussink 2002, Dragovic 1995 |
| `North BJ AND SIRT2 AND tubulin deacetylase` | PMID 12620231 |
| `(NMN OR "nicotinamide mononucleotide") AND BubR1` | 1 hit, PMID 24825348 |
| `nicotinamide mononucleotide AND (pharmacokinetics OR plasma) AND human*` | PMID 31685720 |
| elink `pubmed_pmc` for the mechanism papers; efetch `db=pmc` for PMC4470489, PMC3446968, PMC2886189, PMC4164625, PMC10138348, PMC7922766, PMC5854409, PMC4860349 | full texts obtained |

### 8b. Queries that failed — what I looked for and did not find

**Marcotte 2004 numeric SIRT2 Ki/IC50. I did not find it.** Searches run:
- `Marcotte[Author] AND SIRT` — 1 hit, the paper itself, abstract only.
- esearch `db=pmc` with `15301953[PMID]` — Count = 0. The paper is not in PMC.
- elink `pubmed_pubmed_citedin` for PMID 15301953 — 26 citing papers, screened; none quotes its SIRT2 Ki.
- WebSearch "Marcotte 2004 Analytical Biochemistry SIRT2 nicotinamide Ki value fluorescence trypsin assay" — nothing numeric.
- WebSearch `"Marcotte" 2004 SIRT2 nicotinamide "Ki" OR "IC50" value cited sirtuin assay Abbott` — nothing numeric.
- The ScienceDirect page (S0003269704004580) is abstract-only without a subscription.
Conclusion: **UNRESOLVED**, not absent. Someone with Elsevier access should read pages 90-99 and Figure/Table with the nicotinamide kinetics.

**Tervo 2004 nicotinamide IC50 for SIRT2 (the source of the "100 µM" figure). I did not find it.** ACS *J Med Chem* paywall; abstract does not mention nicotinamide. Searched `Tervo 2004 J Med Chem in silico inhibitors human sirtuin type 2 nicotinamide IC50 100 uM reference compound`. **UNRESOLVED.**

**Feldman 2015 numeric nicotinamide IC50 values for SIRT2 by acyl group. I did not find them in text.** The values exist only in Figure 4 of PMID 25897714; the PMC XML carries the caption but not the plotted numbers. **UNRESOLVED.**

**A nicotinamide IC50 or Ki against SIRT2 using BubR1 or a BubR1-derived peptide. I did not find one.** Searches:
- `SIRT2 AND BubR1 AND nicotinamide` — **exactly 1 hit**, North 2014, which uses 20 mM NAM in cells and reports no enzyme kinetics.
- `BubR1 AND (K668 OR lysine-668 OR acetylation)` — 35 hits, none reporting nicotinamide kinetics on that site.

**A titration of nicotinamide against a SIRT2-dependent output identifying a crossover concentration. I did not find one.** Searches:
- `nicotinamide[Title] AND (biphasic OR dose-dependent) AND (SIRT1 OR sirtuin)` — 9 hits.
- `nicotinamide AND (alpha-tubulin OR tubulin) AND acetylation AND SIRT2` — 14 hits.
- `nicotinamide AND dose[TIAB] AND SIRT2[TIAB]` — 8 hits.
- `nicotinamide treatment AND intracellular NAD+ AND (acetylation OR deacetylation) AND cells[TIAB]` — 24 hits.
The closest is Mitchell 2018 (two doses, in vivo, mouse liver, one figure panel), which is a two-point comparison, not a titration.

**A human measurement of intracellular or tissue free nicotinamide after nicotinamide dosing. I did not find one.** Searches:
- `intracellular nicotinamide concentration AND (measurement OR quantification OR tissue)` — 25 shown of a larger set.
- `nicotinamide AND (intracellular OR tissue) AND concentration AND (supplementation OR administration) AND human*` — 92 hits, top results screened, not exhaustive.
This is the largest remaining hole in section 4.

**Measured plasma nicotinamide in children at 1.2 g/m2/day. I did not find one.** Searches:
- `nicotinamide AND (children OR paediatric OR pediatric) AND plasma concentration` — 20 shown.
- `nicotinamide AND children AND (plasma OR serum) AND (concentration OR level*) AND diabetes` — 15 shown.
- `nicotinamide AND (ENDIT OR "diabetes prevention") AND plasma concentration AND relatives` — 1 hit, no PK.
- `nicotinamide AND "1.2 g" AND (plasma OR serum)` — 4 hits, no levels.

**Measured plasma nicotinamide at 500 mg twice daily. I did not find one.** `nicotinamide AND "500 mg" AND (plasma OR serum)` — 67 hits, top results screened only. Not exhaustive; worth one more pass.

**An immediate-release nicotinic acid Cmax in humans. I did not resolve one.** `niacin AND immediate-release AND plasma AND (concentration OR Cmax)` — 6 hits, none giving a clean Cmax. Reiche I et al., *Nephrol Dial Transplant* 2011, PMID 20562093 has extended-release numbers in its full text which I did not open. **UNRESOLVED.**

**A direct dose-response of free nicotinic acid against any purified sirtuin. I did not find one.** Searches:
- `(isonicotinamide OR "nicotinic acid") AND (Sir2 OR sirtuin) AND (base exchange OR inhibitor*)` — 15 hits.
- `"nicotinic acid"[TIAB] AND (SIRT1 OR SIRT2 OR SIRT3)[TIAB]` — 15 hits.
- `"nicotinic acid" AND "base exchange"` — 12 hits.
- `pyridine* AND Sir2 AND (transglycosidation OR base exchange)` — 1 hit (Jackson 2003), whose Bronsted-plot panel of 3-substituted pyridines may or may not include nicotinic acid; full text not opened. **UNRESOLVED.**
The nicotinic-acid-is-clean conclusion therefore rests on three indirect lines (Gallo 2004, Avalos 2005, Schmidt 2004) rather than a direct measurement. Say "no evidence of inhibition", not "does not inhibit".

**Any NAD+ precursor study in BUB1B-related mosaic variegated aneuploidy. I did not find one.** `(BUB1B OR "mosaic variegated aneuploidy") AND (NAD OR nicotinamide OR sirtuin)` — 6 hits total: an intervertebral disc GWAS, a BubR1 hypertension paper, Qiu 2018 (oocytes), North 2014, the North 2014 commentary, and a 2010 endothelial polyploidy/NAMPT paper. None is a therapeutic NAD+ study in MVA.

**NR or NMN acting as sirtuin inhibitors. I did not find any report.** `("nicotinamide riboside" OR "nicotinamide mononucleotide") AND (sirtuin OR SIRT1 OR SIRT2) AND inhibit*[TIAB]` — 15 hits, none reporting direct inhibition.

### 8c. Access failures (not the same as absence)

| Target | Attempt | Result |
|---|---|---|
| Marcotte 2004 full text | ScienceDirect S0003269704004580 | abstract only, paywalled |
| Schmidt 2004 full text | WebFetch jbc.org | HTTP 403 |
| Petley 1995 full text (Cmax at 25 mg/kg) | WebFetch diabetesjournals.org | HTTP 403 |
| Hwang & Song 2017 full text | efetch db=pmc PMC11107671 | abstract-only deposit; body read via PMC HTML through WebFetch, numbers therefore second-hand |
| North 2014 full text | embopress.org, then link.springer.com | 301 then 303 to an auth endpoint; obtained instead via PMC4194088 HTML |
| Europe PMC `fullTextXML` for PMC4470489 / PMC11107671 / PMC4194088 | REST | returned 0 bytes; NCBI efetch `db=pmc` used instead |
| Tervo 2004 full text | ACS J Med Chem | paywalled |
| Reiche 2011 full text (niacin Cmax in CKD) | not attempted beyond abstract | unresolved |

