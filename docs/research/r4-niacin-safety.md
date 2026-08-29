# R4: Niacin, nicotinamide, and the oncogenic-risk gate

Status: COMPLETE. Written incrementally during research; every section below is sourced.

Scope: public literature only. Every factual claim carries author, year, journal, PMID or NCT.

Question under test: our pipeline cleared **nicotinamide** through the `oncogenic_risk` gate on the strength of one phase 3 trial with a cancer endpoint, and rejected **nicotinic acid (niacin)** as unassessed. Our own document asserts, without citation, that niacin's cancer question is "answerable from existing data" because decades of large cardiovascular outcome trials exist. Nobody checked whether those trials captured cancer.

---

## 1. Verdict

**Yes — niacin's cancer question is answerable from existing data, and in large part it is already answered.** Malignancy counts by randomised arm are public for both large modern niacin trials: HPS2-THRIVE posted 751/12,838 versus 749/12,835 neoplasm serious adverse events on ClinicalTrials.gov over a mean 3.6-year follow-up (RR 1.00, 95% CI 0.91–1.11), and AIM-HIGH published 72/1,718 versus 73/1,696 (P = 0.87) in the New England Journal of Medicine, with HPS2-THRIVE's cancers actively ascertained through national cancer registries and blind-adjudicated per its published design. What does *not* exist, for niacin or for nicotinamide, is any evidence bearing on germline chromosomal instability or embryonal tumour risk in a child — no study of either agent measures chromosome missegregation or aneuploidy in a chromosomal-instability background.

Our internal ledger entry CL-63 is therefore **half right and half wrong**: right that the claim needed a citation, wrong that "malignancies may sit in unreported SAE capture; that is a data-request, not existing published data." The data are reported, and no request is needed.

---

## 2. The large niacin outcome trials: was cancer captured?

Short answer: **yes, in the two largest modern trials, and the per-arm malignancy counts are already public.** They were not published as a cancer endpoint in the primary journal article; they were posted as MedDRA serious-adverse-event data on ClinicalTrials.gov, and in AIM-HIGH's case also printed in a follow-up NEJM safety letter.

| Trial | n | Duration | Cancer endpoint captured? | Result |
|---|---|---|---|---|
| **HPS2-THRIVE** (HPS2-THRIVE Collaborative Group, *N Engl J Med* 2014;371:203-12, PMID 25014686; NCT00461630) | 25,673 | median 3.9 y on treatment; AE window mean 3.6 y | **Not pre-specified** — primary and all four secondary outcomes are vascular or all-cause mortality. But **actively ascertained**: the design paper states the coordinating centre sought information "from national registries (where available) about cancers", blind-adjudicated against pre-specified criteria. Published as a MedDRA System Organ Class **AE category** in the posted ClinicalTrials.gov results. | Neoplasms benign, malignant and unspecified (incl cysts and polyps), serious AEs: **751/12,838 (5.85%) niacin-laropiprant vs 749/12,835 (5.84%) placebo. RR 1.00, 95% CI 0.91–1.11.** No signal. |
| **AIM-HIGH** (AIM-HIGH Investigators, *N Engl J Med* 2011;365:2255-67, PMID 22085343; NCT00120289) | 3,414 | mean 3 y, stopped early for futility | **Not pre-specified.** Published as a MedDRA SOC **AE category** in a dedicated NEJM safety letter (Anderson TJ et al., *N Engl J Med* 2014;371:288-90, PMID 25014706) and, broken out into 77 preferred terms, in the posted ClinicalTrials.gov results. | Neoplasms SOC serious AEs: **72/1,718 (4.19%) niacin vs 73/1,696 (4.30%), P = 0.87. RR 0.97, 95% CI 0.71–1.34.** The only sub-term at P<0.10 favoured placebo (non-small-cell lung cancer 0 vs 4, P=0.06). No signal. |
| **Coronary Drug Project** (Canner PL et al., *J Am Coll Cardiol* 1986;8:1245-55, PMID 3782631) | 8,341 men (niacin arm 1,119) | ~6.2 y on drug; mortality follow-up to mean 15 y | **Not found.** The 1986 paper's abstract reports all-cause mortality only. I did not find a published niacin-vs-placebo cancer-incidence or cancer-mortality analysis from CDP; the trial predates ClinicalTrials.gov. See §7 and §9. | All-cause mortality at 15 y: 52.0% niacin vs 58.2% placebo, 11% relative reduction, P = 0.0004. Mortality moved the right way, which is weak indirect evidence against a large cancer excess. |
| **Stockholm Ischaemic Heart Disease** (Carlson LA and Rosenhamer G, *Acta Med Scand* 1988;223:405-18, PMID 3287837) | 555 | 5 y | **Not captured** in the report I read. Also open-label and niacin was given **combined with clofibrate**, so niacin cannot be isolated. | Total mortality 61 vs 82, 26% reduction, P<0.05. |
| **Cochrane meta-analysis** (Schandelmaier S et al., *Cochrane Database Syst Rev* 2017;6:CD009744, PMID 28616955) | 39,195 across 23 RCTs | median 11.5 months | **Cancer is not a review outcome.** Closest proxy is non-cardiovascular mortality. | Non-cardiovascular mortality **RR 1.12, 95% CI 0.98–1.28** (32,966 participants, 5 studies, I²=0%, high-quality evidence). All-cause mortality RR 1.05 (0.97–1.12). Discontinuation for side effects **RR 2.17 (1.70–2.77)**. |
| **FDA Mini-Sentinel** (Gagne JJ et al., *Pharmacoepidemiol Drug Saf* 2018;27:30-37, PMID 29108128) | 234,242 ER-niacin initiators; 210,389 propensity-matched | 2007–2013 | **Not captured.** Only major GI bleeding and intracranial haemorrhage were assessed. | GI bleeding HR 0.98 (0.82–1.18); ICH HR 1.21 (0.66–2.22). A cancer analysis in this cohort is feasible and has not been done. |

**Every one of these participants was an adult.** Posted baseline data give a mean age of 64.9 years in HPS2-THRIVE and 63.7 in AIM-HIGH, with **zero participants aged 18 or under in either trial**. Whatever the malignancy nulls establish, they do not extend to a child.

Reading the two big numbers together: 751 vs 749 malignancy serious AEs in 25,673 patients over ~3.6 years, plus 72 vs 73 in 3,414 over 3 years. Roughly 1,645 malignancy events, split almost exactly evenly. That excludes anything worse than about a 10% relative excess in malignancy over 3–4 years in older adults.

**What this is not.** None of it is a pre-specified cancer endpoint with an incidence analysis. "Serious adverse events by System Organ Class" is a safety-reporting convention, not an oncology endpoint: it counts events serious enough to be reported, does not systematically screen, and is not adjudicated for cancer specifically — HPS2-THRIVE's registry linkage is the partial exception, and its linked results have not been published. The correct description is: **malignancy counts by randomised arm exist, are public, and are null.**

---

## 3. Verified source notes

### 3.1 HPS2-THRIVE — Landray et al. 2014

HPS2-THRIVE Collaborative Group (Landray MJ et al.), *N Engl J Med* 2014;371:203-12, PMID 25014686 — NCT00461630.

- n = 25,673 adults with established vascular disease.
- Randomised to extended-release niacin 2 g + laropiprant 40 mg daily vs matching placebo.
- Median follow-up 3.9 years.
- Primary outcome: first major vascular event (nonfatal MI, coronary death, stroke, arterial revascularisation). 13.2% vs 13.7%, rate ratio 0.96 (95% CI 0.90–1.03), P=0.29.
- Excess serious adverse events reported in the abstract, by body system: disturbances of diabetes control (+3.7 percentage points), new diabetes diagnoses (+1.3), gastrointestinal (+1.0), musculoskeletal (+0.7), skin (+0.3), infection (+1.4), bleeding (+0.7).
- **The abstract lists no cancer endpoint and no neoplasm SAE category.** Our internal ledger stopped here and concluded that malignancy data "may sit in unreported SAE capture; that is a data-request, not existing published data." **That conclusion is wrong.**

Two things the abstract does not tell you:

**(a) Cancers were actively ascertained by registry linkage.** From the trial design paper, HPS2-THRIVE Collaborative Group (Haynes R et al.), *Eur Heart J* 2013;34:1279-91, PMID 23444397 (PMC3640201), Methods, verbatim:

> "The coordinating centre sought further details from the participant's medical records about all reports that might relate to MVEs or safety outcomes, and from national registries (where available) about cancers and the certified causes of any deaths. All such information was reviewed by coordinating centre clinicians (blind to treatment allocation) and events adjudicated according to pre-specified criteria."

Cancer ascertainment was designed into the trial. The same paper also tabulates cancer as a reason for stopping study treatment during randomised follow-up: **65 (niacin-laropiprant) vs 62 (placebo)**.

**(b) The per-arm malignancy count is posted publicly.** ClinicalTrials.gov NCT00461630, results section, adverse events module, analysed at System Organ Class level over a mean 3.6-year follow-up:

> Neoplasms benign, malignant and unspecified (incl cysts and polyps) — **751/12,838 ER niacin/laropiprant vs 749/12,835 placebo**

RR 1.002, 95% CI 0.91–1.11 (computed here from the posted counts; the trialists publish no test). Total serious AEs, for context: 7,137/12,838 vs 6,762/12,835.

The dedicated HPS2-THRIVE safety paper Haynes R et al., *Clin Ther* 2019;41:1767-1777, PMID 31447131 examines new-onset diabetes, disturbances of diabetes control, bleeding, infection and gastrointestinal upset in depth. It does not analyse cancer. So the deep safety analysis exists but stops short of malignancy; the malignancy numbers are on the registry instead.

### 3.2 Colletti et al. 1993 — paediatric niacin — VERIFIED

Colletti RB et al., *Pediatrics* 1993;92:78-82, PMID 8516088.

Retrospective review, two university referral lipid clinics, all children given single-drug niacin for severe hypercholesterolaemia 1980–1991. Verified against the abstract, verbatim:

- 21 children, aged 4 to 14 years, niacin 500–2250 mg daily.
- Adverse effects in **16 of 21 (76%)**.
- Reversible dose-related aminotransferase elevation in **6 of 21 (29%)**.
- Niacin **discontinued in 8 of 21 (38%)** for flushing, abdominal pain, vomiting, headache, or raised aminotransferases.
- Authors' own conclusion: "Until further study demonstrates long-term safety, niacin treatment should be reserved for the closely-supervised treatment of severe hypercholesterolemia by a lipid specialist."

Our internal figure of 76% / 38% in 21 children is **correct as cited**.

### 3.3 AIM-HIGH — cancer WAS captured and IS published — VERIFIED

Main trial: AIM-HIGH Investigators (Boden WE et al.), *N Engl J Med* 2011;365:2255-67, PMID 22085343 — NCT00120289.
Safety report: Anderson TJ et al., *N Engl J Med* 2014;371:288-90, PMID 25014706 (PMC4156937, open access).

- n = 3,414 (1,718 extended-release niacin 1500–2000 mg/day + simvastatin; 1,696 placebo + simvastatin).
- Stopped early for futility at a mean follow-up of 3 years.
- Primary endpoint: composite of coronary death, nonfatal MI, ischaemic stroke, hospitalisation for ACS, symptom-driven revascularisation. 16.4% vs 16.2%, HR 1.02 (95% CI 0.87–1.21), P=0.79.

The 2014 safety letter publishes a full MedDRA System Organ Class serious-adverse-event table. It contains a neoplasm row, quoted verbatim from Table 1:

> Neoplasms benign, malignant, and unspecified, including cysts and polyps — **72 (4.2%) niacin vs 73 (4.3%) placebo, P = 0.87**
> Non-small-cell lung cancer — 0 niacin vs 4 (0.2%) placebo, P = 0.06

This is the single most important fact in this document. **Cancer data from a randomised, placebo-controlled niacin trial are already published and publicly readable.** They show no excess of neoplasms on niacin over 3 years in 3,414 adults.

What this is, precisely: a MedDRA SOC serious-adverse-event category, coded from standard case-report forms, MedDRA v15.0. It is **not** a pre-specified cancer endpoint, there is no adjudication of malignancy, and no incidence/mortality analysis. It is one line in a safety table. But it is a published count of malignant events by randomised arm, which is exactly what the pipeline's gate was asking for and did not find.

Caveats that limit its weight:
- 3,414 patients, 3 years. Underpowered for cancer as an outcome; the 95% CI around a 4.2% vs 4.3% difference admits a substantial relative excess.
- Adults with established atherosclerotic cardiovascular disease. Per the posted baseline data (NCT00120289), mean age 63.7 years and **zero participants aged 18 or under**. Not children, not a chromosomal-instability background.
- Anderson et al. urge caution in interpreting the serious-AE analysis generally, and explicitly label the serious-infection findings "provisional and exploratory".

### 3.4 Coronary Drug Project — Canner et al. 1986

Canner PL et al., *J Am Coll Cardiol* 1986;8:1245-55, PMID 3782631.

- n = 8,341 men aged 30–64 with ECG-documented prior MI; conducted 1966–1975 across five lipid-influencing drugs plus placebo. (Per-arm sizes are not in the 1986 abstract and I did not verify them from the CDP monographs.)
- Mortality follow-up to a mean of 15 years, roughly 9 years after the trial ended and the drug was stopped.
- Result: all-cause mortality 52.0% niacin vs 58.2% placebo, an 11% relative reduction, P = 0.0004.
- The abstract reports all-cause mortality only. Whether cause-specific (cancer) mortality was tabulated in the paper is checked below.

### 3.5 Cochrane systematic review — Schandelmaier et al. 2017

Schandelmaier S et al., *Cochrane Database Syst Rev* 2017;6:CD009744, PMID 28616955.

- 23 randomised trials, 39,195 participants, published 1968–2015. Median treatment duration 11.5 months, median dose 2 g/day.
- Pre-specified review outcomes are mortality, cardiovascular events, and side effects. **Cancer is not among the listed review outcomes in the abstract.**
- Relevant to the cancer question by proxy: **non-cardiovascular mortality RR 1.12 (95% CI 0.98–1.28), 32,966 participants, 5 studies, I² = 0%, graded high-quality evidence.** Cancer deaths sit inside that composite. The point estimate is above 1 and the interval touches it; it does not exclude a modest excess of non-cardiovascular death.
- Overall mortality RR 1.05 (0.97–1.12), 35,543 participants, 12 studies.
- Discontinuation for side effects: **RR 2.17 (95% CI 1.70–2.77)**, 33,539 participants, 17 studies. Niacin is more than twice as likely as control to be stopped for intolerance.

---

## 4. What the nicotinamide skin-cancer trial does and does not establish

Chen AC et al., *N Engl J Med* 2015;373:1618-26, PMID 26488693 (ONTRAC).

What it is:
- Phase 3, double-blind, 1:1, **386 participants** who had had at least two non-melanoma skin cancers in the previous 5 years. Immunocompetent adults, high-risk, sun-damaged.
- Nicotinamide 500 mg twice daily or placebo for **12 months**; dermatologist review every 3 months to 18 months.
- Primary endpoint: number of new non-melanoma skin cancers (BCC + SCC) during the 12-month intervention.
- Result: rate of new non-melanoma skin cancers **23% lower** (95% CI 4 to 38), P = 0.02. SCC 30% lower (95% CI 0 to 51, P = 0.05); BCC 20% lower (95% CI -6 to 39, P = 0.12). Actinic keratoses 13–20% lower at each timepoint.
- No noteworthy between-group difference in adverse events over 12 months.
- **No benefit persisted after nicotinamide was stopped.**

What it establishes: that 1 g/day oral nicotinamide for 12 months reduces the rate of new UV-driven keratinocyte carcinomas in immunocompetent adults with heavy prior sun damage, and is well tolerated in adults over that period.

What it does **not** establish, and each of these matters for MVA1:

1. **Wrong tumour biology.** Keratinocyte cancers are UV-mutagen-driven epithelial cancers arising in a field of somatically damaged skin. MVA1 tumours are embryonal — rhabdomyosarcoma, Wilms tumour, leukaemia — arising on a germline background of whole-chromosome missegregation. Nothing in Chen 2015 measures chromosome number, missegregation, or aneuploidy. The proposed protective mechanism in the paper is enhancement of repair of UV-induced DNA damage, i.e. nucleotide-level lesions, not mitotic fidelity.
2. **Wrong population and wrong exposure length.** ONTRAC enrolled adults with at least two prior non-melanoma skin cancers — a heavily sun-damaged, older population (I did not verify the exact median age from the paper) — and the intervention ran 12 months. An MVA1 child would be exposed from early childhood, potentially for a decade or more, on a germline instability background. There is no cancer-endpoint data at that duration for either agent.
3. **The effect does not generalise even within skin cancer.** Allen NC et al., *N Engl J Med* 2023;388:804-812, PMID 36856616 (ONTRANS) ran the same drug, dose, and endpoint in 158 organ-transplant recipients with at least two prior keratinocyte cancers: **207 new keratinocyte cancers on nicotinamide vs 210 on placebo, rate ratio 1.0 (95% CI 0.8 to 1.3), P = 0.96.** No difference in SCC, BCC, or actinic keratoses. The trial stopped early for poor recruitment, so it is underpowered, but it is a phase 3 null in a higher-risk, immunosuppressed group.
4. **The pooled picture is heterogeneous.** Mainville L et al., *J Cutan Med Surg* 2022;26:297-308, PMID 35134311 — systematic review, 29 trials, 3,039 patients. Skin cancer rate ratio 0.50 (95% CI 0.29–0.85) but I² = 64%, from only 552 patients across 5 trials, and the authors state heterogeneity was explained by risk of bias. GRADE: weak recommendation, moderate-quality evidence.

**Bottom line for the gate.** Chen 2015 is a real cancer endpoint, and it is genuinely more than niacin has in the form of a designed oncology outcome. But it measures a different cancer, driven by a different mutational process, in a different population, over a different timescale, and its effect did not replicate in the one other phase 3 population tested. It is not evidence about embryonal tumour risk in a child with constitutional aneuploidy.

---

## 5. Paediatric tolerability

### 5.1 Nicotinic acid (niacin) in children

Colletti RB et al., *Pediatrics* 1993;92:78-82, PMID 8516088 remains the primary paediatric tolerability source: 21 children aged 4–14, 500–2250 mg/day, adverse effects in 76%, transaminase elevation in 29%, discontinuation in 38%. Retrospective, uncontrolled, single-drug niacin, two referral clinics, 1980–1991.

I searched for further paediatric niacin experience and found no larger or more recent controlled series (searches listed in §9). What exists is context, not a second cohort:

- The adult evidence points the same way. Schandelmaier S et al., *Cochrane Database Syst Rev* 2017;6:CD009744, PMID 28616955 — discontinuation for side effects **RR 2.17 (95% CI 1.70–2.77)** across 17 trials and 33,539 adults. Haynes R et al., *Clin Ther* 2019;41:1767-1777, PMID 31447131 — in HPS2-THRIVE, hazard ratios of 1.32 for new-onset diabetes, 1.56 for serious disturbance of diabetes control among those already diabetic, 1.38 for serious bleeding, 1.22 for serious infection. Note that these excesses were seen **after** a run-in that had already excluded people who could not tolerate niacin.
- The magnitude of the hepatotoxicity concern in children (29% transaminase rise in Colletti) is the single most relevant number for a chronic paediatric indication, and nothing has been published since to revise it.

**One paediatric case series of 21 children, from 1993, with a 38% discontinuation rate, is the entire paediatric safety base for chronic niacin.** That is thinner than the paediatric base for nicotinamide by a wide margin.

### 5.2 Nicotinamide in children

Nicotinamide has been given to children at high dose for years, in randomised trials, with growth measured.

- Gale EA et al., *Lancet* 2004;363:925-31, PMID 15043959 (ENDIT). 552 relatives of people with type 1 diabetes, with confirmed islet-cell antibodies, recruited across 18 European countries plus Canada and the USA. The cohort mixed children and adults; the abstract does not give the age range or the number of children, and I did not verify those from the full text — but the growth finding quoted below establishes that children were enrolled. Participants were randomised to modified-release nicotinamide **1.2 g/m²/day or placebo for 5 years**. Primary endpoint was diabetes onset; nicotinamide had no effect (adjusted HR 1.01, 95% CI 0.73–1.38). Relevant here for safety, not efficacy: "The number of serious adverse events did not differ between treatment groups. **Nicotinamide treatment did not affect growth in children** or first-phase insulin secretion."
- Lampeter EF et al., *Diabetes* 1998;47:980-4, PMID 9604880 (DENIS). 55 **siblings aged 3–12 years** with high-titre islet-cell antibodies, randomised to slow-release nicotinamide 1.2 g/m²/day or placebo, followed to a maximum of 3.8 years (median 2.1). No effect on diabetes onset. The nicotinamide group showed **decreased first-phase insulin secretion** in response to intravenous glucose (P = 0.03); "no other side effects were observed."
- Further paediatric randomised experience exists in the IMDIAB series, e.g. Crinò A et al., *Eur J Endocrinol* 2004;150:719-24, PMID 15132730 and Visalli N et al., *Diabetes Metab Res Rev* 1999;15:181-5, PMID 10441040.

**Asymmetry on tolerability is real and runs against niacin.** Nicotinamide has multi-year, gram-scale, randomised paediatric exposure with growth data and a comparable serious-AE rate. Niacin has 21 retrospectively reviewed children and a 38% stop rate. This part of the pipeline's preference for nicotinamide is defensible — but it is a *tolerability* argument, and the pipeline recorded it under an *oncogenic-risk* gate. Those are not the same gate.

---

## 6. Does either agent affect chromosomal instability, aneuploidy or micronuclei?

This is the question the `oncogenic_risk` gate should be asking in a chromosomal-instability syndrome, and it is not the question either drug's headline trial answers. Here is what exists.

### 6.1 Nicotinic acid — measured, in micronucleus assays, and the direction is protective

**Human cells, direct micronucleus measurement.**
Weidele K et al., *DNA Repair (Amst)* 2017;52:12-23, PMID 28216063. Nicotinic acid supplementation of resting and proliferating human peripheral blood mononuclear cells raised NAD⁺, increased DNA repair efficiency, and **enhanced genomic stability measured by the micronucleus test after X-irradiation**. Reduced p53 acetylation at higher NAD⁺ suggested increased SIRT1 activity. This is the cleanest single result: the agent our pipeline rejected as unassessed has been assayed for micronucleus formation in primary human cells, and it lowered it.

**Human, observational, dietary.**
Fenech M et al., *Carcinogenesis* 2005;26:991-9, PMID 15705599. 190 healthy adults, food-frequency questionnaire against lymphocyte micronucleus index. Highest tertile of **preformed nicotinic acid intake was associated with a 46% lower micronucleus frequency** than the lowest tertile (P < 0.005) — the second-largest protective association in the study after calcium. Observational, confounded, dietary-range intakes, not therapeutic doses.

**Rodent, in vivo, bone marrow.**
Spronck JC and Kirkland JB, *Mutat Res* 2002;508:83-97, PMID 12379464. Long-Evans rats on niacin-deficient, niacin-replete pair-fed, or nicotinic acid-supplemented (4 g/kg diet) diets for 3 weeks; comet, micronucleus and sister-chromatid-exchange assays in bone marrow, with and without etoposide. Niacin deficiency alone caused **6.2-fold higher spontaneous micronucleus formation and 2.8-fold higher SCE frequency**. Two findings matter for us, and the second cuts against the mechanistic hope:
- Deficiency destabilises chromosomes. Adequacy protects.
- The authors' own conclusion: "**Pharmacological intakes of niacin do not appear to be further protective compared to adequate intakes.**" In a niacin-replete animal, supraphysiological nicotinic acid bought no additional chromosomal stability.

Kostecki LM et al., *Mutat Res* 2007;625:50-61, PMID 17618655 — niacin deficiency delays DNA excision repair and increases spontaneous and nitrosourea-induced chromosomal instability in rat bone marrow. Same direction.

Boyonoski AC et al., *J Nutr* 2002;132:115-20, PMID 11773517 — pharmacological niacin intake increases bone marrow poly(ADP-ribose) and **lengthens the latency of ethylnitrosourea-induced carcinogenesis** in rats. Companion paper Boyonoski AC et al., *J Nutr* 2002;132:108-14, PMID 11773516 shows deficiency shortens that latency.

Bartleman AP et al., *Nutr Cancer* 2008;60:251-8, PMID 18444158 — niacin supplementation **decreases the incidence of alkylation-induced non-lymphocytic leukaemia** in Long-Evans rats.

Reviews synthesising this body of work: Kirkland JB, *Mutat Res* 2012;733:14-20, PMID 22138132 ("Niacin requirements for genomic stability") and Kirkland JB, *Subcell Biochem* 2012;56:21-36, PMID 22116692; earlier, Hageman GJ and Stierum RH, *Mutat Res* 2001;475:45-56, PMID 11295153.

**Summary for nicotinic acid.** Every measurement located points the same way: niacin status supports chromosomal stability through NAD⁺-dependent poly(ADP-ribose) metabolism, deficiency raises micronuclei and shortens carcinogenesis latency, and supplementation is protective or neutral. Nothing found suggests niacin destabilises chromosomes. The important caveat is Spronck and Kirkland's own: above adequacy, more did not help.

### 6.2 Nicotinamide — the direct measurements are thinner, and its pharmacology touches the MVA1 axis

I found **no** study measuring chromosome missegregation, aneuploidy, or micronucleus frequency under therapeutic nicotinamide in humans, and no such study in a chromosomal-instability background. Searches listed in §9.

What does exist is a mechanistic concern that runs in the *wrong* direction for MVA1, and that nobody has closed:

1. **Nicotinamide is a sirtuin product inhibitor.** Bitterman KJ et al., *J Biol Chem* 2002;277:45099-107, PMID 12297502 — "Inhibition of silencing and accelerated aging by nicotinamide, a putative negative regulator of yeast Sir2 and human SIRT1."
2. **SIRT2 is the enzyme that keeps BubR1 up.** North BJ et al., *EMBO J* 2014;33:1438-53, PMID 24825348 — the age-related fall in BubR1 is due to declining NAD⁺ and loss of SIRT2's ability to keep BubR1 lysine-668 deacetylated; SIRT2 overexpression or NMN raises BubR1 in vivo and extends median lifespan in BubR1-hypomorphic mice.
3. **Losing SIRT2 causes exactly the phenotype MVA1 already has.** Kim HS et al., *Cancer Cell* 2011;20:487-99, PMID 22014574 — Sirt2-null mice show elevated Aurora-A and -B, centrosome amplification, **aneuploidy**, mitotic cell death, and develop tumours; SIRT2 acts through deacetylating the APC/C coactivators CDH1 and CDC20.
4. Recent synthesis of the same axis: Pun R et al., *Semin Cancer Biol* 2024;106-107:201-216, PMID 39490401 — "BubR1 and SIRT2: Insights into aneuploidy, aging, and cancer."

Put together: the drug our pipeline cleared through the oncogenic-risk gate is a product inhibitor of the sirtuin that maintains BubR1, in a disease defined by insufficient BubR1. Whether plasma nicotinamide at 500 mg twice daily reaches SIRT2-inhibitory concentrations in the relevant tissues, and what the net effect on BubR1 abundance is, **I did not establish** — I found no study that measures it. That is the point. Nicotinamide's cancer file (Chen 2015) is about UV-induced DNA lesions in skin. It says nothing about the mitotic axis, and the mitotic axis is where nicotinamide has a plausible mechanism for making things worse.

There is one directly relevant piece of positive evidence for nicotinamide on the DNA-lesion side, from the same Sydney group that ran ONTRAC: Surjana D et al., *Carcinogenesis* 2013;34:1144-9, PMID 23349012 — nicotinamide enhances repair of ultraviolet-induced DNA damage in human keratinocytes and ex vivo skin. That is nucleotide-excision-repair biology on point lesions. It is not chromosome-segregation biology, and it does not speak to whole-chromosome missegregation.

---

## 7. Other niacin trials, and what I could not find

- **Coronary Drug Project.** Canner PL et al., *J Am Coll Cardiol* 1986;8:1245-55, PMID 3782631. I did **not** find a published niacin-versus-placebo cancer-incidence or cancer-mortality analysis from CDP. The 1986 paper's abstract reports all-cause mortality only (52.0% vs 58.2% at 15 years, P = 0.0004). Searching PubMed for `"Coronary Drug Project" AND (cancer OR neoplasm OR "cause of death" OR "cause-specific mortality")` returned 2 hits, neither of which is a CDP cause-specific-mortality report. CDP predates ClinicalTrials.gov, so there is no posted results record. I did not obtain the CDP monographs, so I write "I did not find", not "does not exist". The 11% all-cause mortality reduction is weak indirect evidence against a large cancer excess, since cancer deaths sit inside all-cause mortality.
- **Stockholm Ischaemic Heart Disease Secondary Prevention Study.** Carlson LA and Rosenhamer G, *Acta Med Scand* 1988;223:405-18, PMID 3287837. n = 555 post-MI patients, 5 years, **open-label**, and nicotinic acid was given **in combination with clofibrate**, so niacin cannot be isolated. Total mortality 61 (treatment) vs 82 (control), a 26% reduction, P < 0.05. No cancer analysis in the abstract; I did not find one.
- **ARBITER 6-HALTS** (Villines TC et al., *J Am Coll Cardiol* 2010;55:2721-6, PMID 20399059) and the surrogate-endpoint niacin trials (HATS, CLAS, FATS) are imaging or biomarker studies of a few hundred patients over 1–3 years. They are too small and too short for a cancer signal, and I did not find cancer analyses in them.
- **FDA Mini-Sentinel.** Gagne JJ et al., *Pharmacoepidemiol Drug Saf* 2018;27:30-37, PMID 29108128. 234,242 extended-release niacin initiators, 210,389 propensity-matched to fenofibrate initiators. This is the largest real-world niacin safety cohort I located. It assessed **only** major gastrointestinal bleeding (HR 0.98, 95% CI 0.82–1.18) and intracranial haemorrhage (HR 1.21, 95% CI 0.66–2.22). **Cancer was not an outcome.** Note the size: a cancer analysis in this cohort is feasible and has not been done — that would be a genuine data request, and a cheap one.
- **A caution that applies to both agents.** Nakazzi F et al., *Cancer Lett* 2026;645:218334, PMID 41724424 reports that vitamin B3 derivatives support pancreatic cancer cell survival and chemotherapy resistance. The wider preclinical NAD⁺-and-cancer literature is genuinely mixed — NAD⁺ supplementation suppresses metastasis in some models and supports tumour cell survival in others. None of it is specific to nicotinic acid versus nicotinamide, and none of it is in a CIN background. It is a reason to keep the gate open for both agents, not a reason to discriminate between them.

---

## 8. Direct answer: is our pipeline's asymmetry justified?

**No. As an `oncogenic_risk` decision, it is an artefact of which trial happened to exist.**

The gate records whether a publication with a cancer endpoint exists. It does not record whether that evidence transfers to this patient. On the substance:

1. **The factual premise for rejecting niacin is false.** Niacin does have per-arm malignancy data from randomised, placebo-controlled trials — 25,673 patients in HPS2-THRIVE and 3,414 in AIM-HIGH, roughly 1,645 malignancy events, split evenly. Nicotinamide's phase 3 cancer trial enrolled 386 people for 12 months. On sample size and exposure, niacin's cancer safety data are the larger body.
2. **On the specific question a CIN syndrome cares about, niacin is the better-characterised agent, not the worse one.** Nicotinic acid has been put through the micronucleus assay in primary human cells (Weidele K et al., *DNA Repair (Amst)* 2017;52:12-23, PMID 28216063), in rat bone marrow in vivo with micronucleus and sister-chromatid-exchange readouts (Spronck JC and Kirkland JB, *Mutat Res* 2002;508:83-97, PMID 12379464), and in carcinogenesis-latency and leukaemia-incidence models (Boyonoski AC et al., *J Nutr* 2002;132:115-20, PMID 11773517; Bartleman AP et al., *Nutr Cancer* 2008;60:251-8, PMID 18444158). Every result points protective or neutral. For nicotinamide I found no equivalent measurement at all.
3. **Nicotinamide carries an unclosed mechanistic liability on exactly the MVA1 axis.** It is a sirtuin product inhibitor (Bitterman KJ et al., *J Biol Chem* 2002;277:45099-107, PMID 12297502); SIRT2 is what maintains BubR1 (North BJ et al., *EMBO J* 2014;33:1438-53, PMID 24825348); and losing SIRT2 produces aneuploidy and tumours in mice (Kim HS et al., *Cancer Cell* 2011;20:487-99, PMID 22014574). Chen 2015 does not test this. A skin-cancer chemoprevention result in sun-damaged older adults is not evidence that the drug is safe for the spindle checkpoint in a BUB1B-deficient child. The gate cleared nicotinamide on evidence that is orthogonal to its own most plausible harm mechanism.
4. **The one asymmetry that is real runs on tolerability, not oncogenicity.** Nicotinamide has 552 subjects (children and adults) on 1.2 g/m²/day for 5 years with growth measured in the children and no serious-AE excess (Gale EA et al., *Lancet* 2004;363:925-31, PMID 15043959) plus 55 children aged 3–12 for up to 3.8 years (Lampeter EF et al., *Diabetes* 1998;47:980-4, PMID 9604880). Niacin has 21 retrospectively reviewed children with 76% adverse effects and 38% discontinuation (Colletti RB et al., *Pediatrics* 1993;92:78-82, PMID 8516088), plus an adult discontinuation risk ratio of 2.17 (Schandelmaier S et al., *Cochrane Database Syst Rev* 2017;6:CD009744, PMID 28616955) and HPS2-THRIVE's excesses of diabetes, infection and bleeding. **Preferring nicotinamide on chronic paediatric tolerability is well supported.** It just is not what the `oncogenic_risk` gate is for.

**Recommended correction to the pipeline.** Either (a) niacin should pass the `oncogenic_risk` gate on the HPS2-THRIVE and AIM-HIGH neoplasm data and then be ranked down on a separate paediatric-tolerability gate, or (b) both agents should fail the gate, because neither has any evidence about chromosomal instability or embryonal tumours in a child and the gate's true question is unanswered for both. What is not defensible is the current state, where one agent passes on evidence that does not transfer and the other fails on evidence that was never looked for.

Concretely, three cheap actions would close the gate properly:
1. Cite the two posted neoplasm SAE tables (NCT00461630, NCT00120289) — no data request needed, they are public today.
2. Request cancer-registry-linked incidence by arm from the HPS2-THRIVE coordinating centre. Per HPS2-THRIVE Collaborative Group, *Eur Heart J* 2013;34:1279-91, PMID 23444397, that linkage was performed and blind-adjudicated. This is the only genuine data request on the list.
3. Run a micronucleus / aneuploidy assay on both agents in BUB1B-deficient patient cells. Neither has been tested in a CIN background, and this is the only experiment that answers the question actually being asked.

---

## 9. Searches run

Databases: PubMed E-utilities (esearch → esummary → efetch), PubMed Central efetch, ClinicalTrials.gov API v2, Europe PMC REST.

| # | Database | Query | Result |
|---|---|---|---|
| S1 | PubMed esummary | `26488693,25014686,8516088` | All three PMIDs resolve to the cited papers. Citations correct as written in our documents. |
| S2 | PubMed | `HPS2-THRIVE` (44 hits) | Found the dedicated safety paper, Haynes 2019, PMID 31447131, and the design paper, PMID 23444397. |
| S3 | PubMed | `niacin AND cancer AND randomized controlled trial[pt] AND (incidence OR mortality)` (16 hits) | **Found no trial with niacin and a pre-specified cancer endpoint.** Most hits are the Linxian nutrition-intervention trials, where niacin appears only as one component of a multi-vitamin/mineral combination and cannot be isolated. |
| S4 | PubMed | `AIM-HIGH niacin trial` (56 hits) | Located the AIM-HIGH publication family. |
| S5 | PubMed | `"AIM-HIGH Investigators"[cn]` (7 hits) | Found Anderson 2014, PMID 25014706 — the safety letter containing the neoplasm SOC row. |
| S6 | PubMed | `Coronary Drug Project niacin` (31 hits) | Canner metabolic-syndrome and glycaemic-status sub-analyses; **no cancer analysis**. |
| S7 | PubMed | `Canner PL Coronary Drug Project fifteen year mortality` (1 hit) | PMID 3782631, all-cause mortality only. |
| S8 | PubMed | `"Coronary Drug Project" AND (cancer OR neoplasm OR "cause of death" OR "cause-specific mortality")` (2 hits) | **Found nothing.** Neither hit is a CDP cause-specific-mortality report. I did not find CDP cancer data; I did not obtain the CDP monographs, so this is "not found", not "does not exist". |
| S9 | PubMed | `Clofibrate and niacin in coronary heart disease Coronary Drug Project Research Group` | **0 hits.** The 1975 CDP primary report is not retrievable by this phrasing. |
| S10 | PubMed | `Schandelmaier niacin Cochrane` (1 hit) | PMID 28616955. Cancer is not a review outcome; non-cardiovascular mortality is the closest proxy. |
| S11 | PMC efetch | PMC4156937 (Anderson 2014) | Retrieved Table 1 in full including the neoplasm SOC row, 72 vs 73, P = 0.87. |
| S12 | PMC efetch | PMC3640201 (HPS2-THRIVE design, Eur Heart J 2013) | Found the registry-linkage sentence and the cancer-as-reason-for-stopping counts, 65 vs 62. |
| S13 | Europe PMC REST | `PMC6481694/fullTextXML` (Cochrane review) | **Returned 0 bytes.** Full text of the Cochrane review could not be retrieved this way. Conclusions about the Cochrane review rest on its structured abstract only. |
| S14 | ClinicalTrials.gov API v2 | `studies/NCT00461630` | HPS2-THRIVE. Outcomes module: primary and all four secondary outcomes are vascular or all-cause mortality — **cancer is not a pre-specified outcome**. Results module: neoplasm SOC serious AEs 751/12,838 vs 749/12,835. |
| S15 | ClinicalTrials.gov API v2 | `studies/NCT00120289` | AIM-HIGH. Cancer not a pre-specified outcome; 77 individual neoplasm preferred terms posted per arm. |
| S15b | ClinicalTrials.gov API v2 | baseline characteristics modules of NCT00461630 and NCT00120289 | Mean age 64.9 y (HPS2-THRIVE) and 63.7 y (AIM-HIGH); **zero participants aged 18 or under in either trial**. |
| S16 | PubMed | `(niacin OR "nicotinic acid") AND (child OR children OR pediatric OR paediatric OR adolescent) AND (hypercholesterolemia OR dyslipidemia OR tolerability OR safety)` (196 hits, top 25 reviewed) | **Found no paediatric niacin cohort larger or more recent than Colletti 1993.** Hits are adult dyslipidaemia, dietary-intake surveys, and unrelated topics. One adult glioblastoma phase I/II niacin trial (PMID 41313494) and one pregnancy case report (PMID 39949380) are the only new clinical niacin exposures found. |
| S17 | PubMed | `nicotinamide AND (ENDIT OR "type 1 diabetes") AND children AND randomized` (14 hits) | Found ENDIT (PMID 15043959), DENIS (PMID 9604880), IMDIAB IX (PMID 15132730), IMDIAB VI (PMID 10441040). Substantial paediatric nicotinamide exposure confirmed. |
| S18 | PubMed | `(niacin OR "nicotinic acid" OR nicotinamide) AND (micronucleus OR micronuclei OR aneuploidy OR "chromosomal instability" OR "chromosome missegregation" OR "genomic instability")` (109 hits, top 30 reviewed) | Too noisy; dominated by NAD⁺-metabolism reviews. Refined below. |
| S19 | PubMed | `"niacin deficiency" AND (chromosom* OR micronucle* OR "genomic instability" OR leukemia)` (15 hits) | **The productive query.** Yielded Spronck 2002, Kostecki 2007, Boyonoski 2000/2002 ×2, Bartleman 2008, Kirkland 2012 ×2, Hageman 2001. |
| S20 | PubMed | `Fenech micronucleus micronutrient intake genome instability lymphocytes` (1 hit) | PMID 15705599 — nicotinic acid intake and human lymphocyte micronucleus index. |
| S21 | PubMed | `nicotinamide AND micronucleus AND (lymphocyte OR keratinocyte OR human)` (10 hits) | **Found no study measuring micronucleus frequency under therapeutic nicotinamide in humans.** The nicotinamide hits are tumour radiosensitisation studies where nicotinamide is used for its blood-flow effect (PMID 8248558, PMID 7842826), not genotoxicity assessments. The one nicotinic-acid hit, Weidele 2017, is the human-cell micronucleus result. |
| S22 | PubMed | `nicotinamide AND ("sister chromatid exchange" OR "chromosome aberration" OR clastogenic OR genotoxic)` (104 hits, top 15 reviewed) | **Found no chromosome-aberration or SCE study of therapeutic nicotinamide.** Hits are NMN toxicology packages and unrelated compounds. |
| S23 | PubMed | `nicotinamide AND ("mitotic" OR "spindle" OR "chromosome segregation") AND (oocyte OR embryo OR cell)` (92 hits, top 15 reviewed) | **Nothing on nicotinamide itself.** All relevant hits concern NMN or nicotinamide riboside in oocytes, not nicotinamide, and none in a CIN background. |
| S24 | PubMed | `SIRT2 AND (aneuploidy OR "genome integrity" OR "chromosomal instability" OR mitotic)` (55 hits) | Kim 2011 (PMID 22014574), Pun 2024 (PMID 39490401). |
| S25 | PubMed | `Kim HS SIRT2 genome integrity tumorigenesis APC/C`; `North BJ SIRT2 BubR1 deacetylation lifespan` | PMID 22014574, PMID 24825348. |
| S26 | PubMed | `nicotinamide inhibition of silencing accelerated aging Sir2 SIRT1` (1 hit) | PMID 12297502 (Bitterman 2002). An earlier phrasing including "Bitterman nicotinamide inhibition Sir2" returned **0 hits** — the author-plus-concept form failed where the title-phrase form worked. |
| S27 | PubMed | `nicotinamide skin cancer chemoprevention transplant recipients randomized` (5 hits) | ONTRANS (PMID 36856616), the phase II transplant trial (PMID 27061568), and the Mainville meta-analysis (PMID 35134311). |
| S28 | PubMed | `(Stockholm ... Carlson) OR (HATS ...) OR (ARBITER 6-HALTS ...)` (7 hits) | Stockholm study PMID 3287837 located. **No cancer analysis found in any of these trials.** |
| S29 | PubMed | `(niacin OR "nicotinic acid") AND cancer AND (cohort OR "case-control" OR observational OR risk)` (193 hits, top 20 reviewed) | **Found no observational study of therapeutic niacin exposure and cancer incidence.** Hits are dietary-niacin-intake epidemiology and NAMPT-inhibitor oncology, which are different questions. |
| S30 | PubMed | `NAD precursor supplementation AND (tumor promotion OR tumour growth OR carcinogenesis OR metastasis)` (24 hits) | Mixed preclinical literature; nothing agent-specific to nicotinic acid vs nicotinamide and nothing in a CIN background. |
| S31 | PubMed esummary | 22-PMID bulk verification | Volume/issue/page details confirmed for every citation in this document. One correction made: Pun 2024 is *Semin Cancer Biol* 106-107:201-216, not 107-108:113-125. |

### Negatives, stated precisely

Each of these is "I did not find", not "does not exist":

- I did not find any randomised trial in which cancer incidence or cancer mortality was a **pre-specified** endpoint for nicotinic acid.
- I did not find a published cancer-incidence or cancer-mortality analysis by treatment arm from the Coronary Drug Project.
- I did not find a paediatric niacin cohort larger or more recent than the 21 children in Colletti 1993.
- I did not find any study measuring micronucleus frequency, chromosome aberrations, sister-chromatid exchange, aneuploidy, or chromosome missegregation under therapeutic nicotinamide, in humans or in animals.
- I did not find any study of either agent in a chromosomal-instability background — no BUB1B-mutant, MVA, or aneuploidy-syndrome model, and no measurement of BubR1 abundance under either drug.
- I did not find the cancer-registry-linked incidence results that the HPS2-THRIVE design paper says were collected. The posted SOC-level SAE count is the closest public substitute.
- I did not retrieve the full text of the HPS2-THRIVE primary NEJM article or the Cochrane review, so statements about what they contain rest on the structured abstract plus the ClinicalTrials.gov posting in the first case, and the structured abstract alone in the second.
