# Track 2 — Claims ledger: every cited claim, re-checked against its source

Companion to `docs/track2-hypothesis.md`. Nothing here is new science. It is the
result of going back to **77 substantive factual claims** in that document and asking a
harder question than "does the PMID resolve":

> **Does the cited paper actually support the specific claim attached to it?**

---

## 0. Method, and what the grades mean

Every PMID in `docs/track2-hypothesis.md` and `knowledge/literature/bub1b/SOURCES.md`
(54 unique) was re-fetched from PubMed. Where the claim turned on something an abstract
cannot settle — a dose, a model organism, a sex effect, a count of mutations, whether a
guideline covers a population — the PMC or Europe PMC full text was read. Regulatory and
registry claims were checked against the EMA product record and the ClinicalTrials.gov
API. Dates of check: 2026-08-27/28.

| Grade | Meaning |
|---|---|
| **HOLDS** | The cited source supports the claim as written. |
| **WEAKENED** | The source supports a narrower, more qualified, or differently-scoped claim. The sentence must be rewritten; the conclusion survives. |
| **WITHDRAWN** | The cited source does not support the claim. It must be removed or re-sourced before it can stand. |
| **UNVERIFIED** | Could not be checked from the source available. Reason stated. Not a criticism of the claim — a statement about what I could see. |

**Headline:**

| Verdict | Count |
|---|---|
| HOLDS | **47** |
| WEAKENED | **22** |
| WITHDRAWN | **3** |
| UNVERIFIED | **5** |
| **Total re-checked** | **77** |

Six of the seventy-seven carry **no citation and no inference label**, which the rule in §0
(Scope, and what the citations guarantee) forbids; they are graded UNVERIFIED or flagged in
§7 of this ledger.

> **Status note, added 2026-08-29.** This ledger is a dated audit of
> `docs/track2-hypothesis.md` as it stood on 2026-08-27/28, and it is kept as a record
> rather than rewritten. Two things it says have since been overtaken by follow-up
> research and by a rerun of the pipeline, and a reader should not take either at face
> value. **CL-63 is now resolved and its verdict was wrong in the other direction:** the
> per-arm malignancy counts for niacin *are* published — ClinicalTrials.gov NCT00461630 and
> Anderson et al., *N Engl J Med* 2014, PMID 25014706 — so no data request is needed, and
> `worsens_cin` for nicotinic acid is now recorded as assessed-and-negative in adults.
> **CL-62 still holds as a description of its trial but no longer supports the cell it was
> filling:** nicotinamide's oncogenic-risk cell has been blanked, because a UV-driven
> skin-cancer endpoint does not bear on chromosomal instability and did not replicate in
> transplant recipients. The consequence is that the ranking inverted — nicotinic acid is
> now the accepted rank-1 candidate and nicotinamide is rejected — so every cross-reference
> below to "§5.1 makes it the reason nicotinamide is the accepted candidate" describes a
> superseded version of that section. The current reasoning is in §4.4 (The sign paradox at
> the molecule, and which molecule the experiment uses), §4.5 (Does niacin reach the cell?
> NAPRT, and a claim that turned out to be wrong), §5.1 (The gate that should not have been
> open) and §5.2 (The gate that should not have been shut) of the hypothesis document.
>
> **Status of the corrections queue, re-checked 2026-08-29.** Every correction C-1 to C-16
> in §7 was read back against `docs/track2-hypothesis.md` as it now stands, and the queue is
> almost closed. **Thirteen are applied** — C-1 (in the corrected form the W-1 CORRECTION
> describes, not the wording quoted in the entry), C-2, C-3, C-5, C-6, C-7, C-9, C-10, C-12,
> C-13, C-14, C-15 and C-16. **Two are overtaken** by later rewrites of the sections they
> targeted, C-4 and C-11, and are marked as such. **One is partial, C-8**, and it is the only
> item a reader should still chase: §3 (Direction of effect) now rejects hydroxychloroquine
> and metformin on the evidence that recommends their *analogues*, which was the half that
> mattered, but it still does not name the model system, and the primary demonstration in
> **trisomic mouse embryonic fibroblasts** is a point in the hypothesis's favour, since
> constitutional trisomy is the better model for MVA than a transformed line. The
> "Label, don't rewrite" list at the end of §7 is also closed: all six claims that carried
> neither a citation nor an inference label now carry one. Each entry in §7 carries its own
> marker and the original text underneath it, because what was asked for is part of the
> record.

Three withdrawals and twenty-two weakenings out of seventy-seven is not a good result and
is not meant to look like one. The document's *conclusions* survive every one of them —
in two cases the correction makes the hypothesis look **better**, not worse, which is
itself worth noticing about how the errors were distributed. What does not survive is the
document's own claim in §0 (Scope, and what the citations guarantee) that "every factual
claim below carries an author, year, journal and PMID". Several claims carry a PMID that
does not carry them.

---

## 1. The three withdrawals, in full

### W-1 — Nicotinamide's inhibition of **SIRT2** is not what Bitterman 2002 measured

**Claim as written** (§1, The answer; Appendix B, The therapeutic classes considered;
§4.4, The sign paradox at the molecule): nicotinamide is "simultaneously an inhibitor of
the sirtuin the mechanism depends on", cited to Bitterman et al., *J Biol Chem*
2002;277:45099-107 (PMID 12297502).

**What the source says.** Bitterman measured **yeast Sir2 and human SIRT1**. Title:
"…a putative negative regulator of yeast sir2 and human SIRT1". Abstract: "physiological
concentrations of nicotinamide noncompetitively inhibit both Sir2 and SIRT1 in vitro
(IC50 < 50 µM)". **SIRT2 is not assayed and is not mentioned.**

The mechanism in §4.3 (Why the NAD⁺/SIRT2 axis, and where it runs out) depends on **SIRT2**,
not SIRT1 (North 2014, PMID 24825348).

**Why this matters more than a footnote.** §4.4 (The sign paradox at the molecule, and which
molecule the experiment uses) calls this "the nicotinamide sign paradox" and §5.1 (The gate
that should not have been open) makes it the reason nicotinamide is the accepted candidate
with an undeterminable direction. The paradox is the document's stated non-obvious finding,
and its only citation is about a different enzyme.

**What can be said instead.** Nicotinamide inhibition is a **class property of the
conserved sirtuin nicotinamide-binding pocket**, established structurally and
mechanistically: Avalos, Bever & Wolberger, *Mol Cell* 2005;17:855-68 (PMID 15780941)
show free nicotinamide binds "a conserved pocket that participates in NAD⁺ binding and
catalysis" and drives base exchange at the expense of deacetylation; Jackson, Schmidt,
Oppenheimer & Denu, *J Biol Chem* 2003;278:50985-98 (PMID 14522996) establish the
transglycosidation mechanism across Sir2 enzymes; Schmidt, Smith, Jackson & Denu,
*J Biol Chem* 2004;279:40122-9 (PMID 15269219) find nicotinamide "the most potent
inhibitor examined" of Sir2 enzymes. **I did not find a primary paper reporting a
nicotinamide IC50 against human SIRT2 specifically.**

So the honest sentence is *weaker and still sufficient*: nicotinamide inhibits sirtuins
through a pocket conserved across the family, so SIRT2 inhibition is a well-founded
**class inference** — not a measurement, and therefore exactly the kind of unsigned
direction the document's own tri-state is built to represent. The `context_dependent`
verdict and the `cannot determine` outcome are unchanged. The word that has to go is the
implied certainty.

#### W-1 CORRECTION — the paragraph above is wrong, and an external review found it

The original text is kept as written, because how the error was made is the point.

Marcotte PA, Richardson PL, Guo J, Barrett LW, Xu N, Gunasekera A, Glaser KB (2004),
*Anal Biochem* 332(1):90-9, Marcotte et al., *Anal Biochem* 2004, PMID 15301953, states in its abstract that "the kinetics of
the reaction of the cosubstrate nicotinamide adenine dinucleotide and the competitive
inhibitor nicotinamide with **SIRT1 and SIRT2** has been analyzed". Nicotinamide
inhibition of human SIRT2 is **measured**, not a class inference. Added as R58.

**The mechanism of the error.** The sentence above — "I did not find a primary paper
reporting a nicotinamide IC50 against human SIRT2 specifically" — is a true report of a
search. §1 (The answer) and §4.4 (The sign paradox at the molecule) of the hypothesis
document then rendered it as *"SIRT2 specifically has never been measured against it"* and
*"I could not find a primary paper … at all.
Which is the point."* That conversion is **GP-14 violated inside a citation correction**:
an unsuccessful search became a negative finding, and it survived because the audit that
produced it was checking the citation it already had rather than looking for a better one.
A correction is a claim and needs the same evidentiary discipline as the claim it replaces.

**What changes, and what does not.** The potency of nicotinamide against SIRT2 is
substrate- and assay-dependent, so no single IC50 transfers. The verdict is unchanged —
`context_dependent` at the node, `cannot determine` at the molecule — but the reason is
now the right one: nicotinamide simultaneously replenishes the NAD⁺ pool SIRT2 consumes
and supplies the product that inhibits it, and the **net intracellular effect on BubR1
abundance** has never been measured in any genotype. Applied in `docs/track2-hypothesis.md`
§0 (Scope, and what the citations guarantee), §1 (The answer) and §4.4 (The sign paradox at
the molecule), and in `knowledge/literature/bub1b/SOURCES.md` (DRUG-NAM row).

---

### W-2 — "Half" of the Suijkerbuijk mutation classes are qualitative defects. They are not.

**Claim as written** (Appendix C, The four objections, in full — second objection):
*"Suijkerbuijk 2010 describes two mutation classes. Half of them are qualitative defects
that more protein will not fix."* — answered **"Correct"**.

**What the source says** (PMC2887387, full text):

- **Five** kinase-domain-localised substitutions (I909T, R727C, R814H, L844F, L1012P)
  show a 5–10 fold decrease in protein level and are **abundance-lowering**.
- Truncations are abundance-lowering through NMD.
- Of the two stable substitutions located away from the kinase domain: **"expression of
  LAP-R550Q rescued both functions of BUBR1"** — R550Q is functionally normal.
  **Y155C alone** is a genuine separation-of-function allele: it "restored chromosome
  alignment … but failed to reconstitute mitotic checkpoint activity".
- Q921H "was indistinguishable from wild-type LAP-BUBR1, both on protein level and
  functionality" — a passenger on the L844F allele.

So of the characterised substitutions, **one** (Y155C) is a qualitative defect that more
protein will not fix, not half. The paper's own conclusion is the opposite of the
objection's implication: *"we propose that chromosomal instability in MVA patients
carrying BUB1B mutations is a result of low BUBR1 protein abundance."* And even the
mixed group is described as "a reduction in protein **combined with** direct functional
impairment" — the abundance node is deviated in every patient.

**Effect on the hypothesis: this correction is in our favour and we got it wrong against
ourselves.** The eligible allele population is substantially larger than the objection
conceded (Appendix C of the hypothesis document, The four objections, in full). The
genotype-conditionality of ASSUMPTION-DRUG-09 and the first-experiment requirement of §6.3
(The experiment that comes before either) both stand — the class still has to be determined
per patient — but the prior that a given patient falls in the treatable class is much better
than "half".

One caveat that the correction does **not** remove, and which Appendix C (The four
objections, in full) should say instead:
even in the abundance class, the kinase-domain mutants are destabilised by **misfolding**
— their levels collapse under HSP90 inhibition and are rescued by proteasome inhibition
(Suijkerbuijk 2010, Fig. 3C). Whether the K668 deacetylation route can counteract
misfolding-driven turnover is untested. That is the sharper version of the objection
Appendix C is trying to raise, and it is better evidenced than the one it does raise.

---

### W-3 — Kalish 2024 does not cover MVA or *BUB1B*

**Claim as written** (§7.1, What the model could not express; Appendix E.1, The three
gaps): tumour surveillance is supported by "established protocols for
Wilms tumour surveillance in at-risk children (Scott et al., *Arch Dis Child* 2006;91:995-9,
PMID 16857697; **updated in Kalish et al., *Clin Cancer Res* 2024;30:5260-5269, PMID
39320341**)".

**What the source says.** Full text of PMC11611621 checked programmatically for
`BUB1B`, `variegated`, `chromatid`, `MVA`: **zero occurrences**, in body, tables and
reference list alike. The guideline covers Beckwith-Wiedemann, WAGR, *WT1*-related
disorders, Simpson-Golabi-Behmel, Perlman, Bohring-Opitz, trisomy 18, osteopathia striata
with cranial sclerosis, Mulibrey nanism, *PIK3CA*-related overgrowth, Sotos and Weaver.
MVA and *BUB1B* are not among them.

The Kalish citation therefore cannot be used to say that surveillance protocols have been
"updated" for this child's condition. It establishes a protocol for other syndromes.

**What survives, and how it must be argued instead.** Scott 2006's operative rule is
population-agnostic: *"Surveillance should be offered to children at >5% risk of Wilms
tumour"*, by renal ultrasonography every 3–4 months. *BUB1B*-MVA clears that threshold on
any published series — see WK-65 below. The argument is available; it has to be **made**,
by naming the risk estimate and the >5% rule, not delegated to a guideline that does not
list the syndrome.

---

## 2. The material weakenings

The full table is §3 of this ledger. These are the ones that change what the document
should say.

### WK-11 / WK-46 — The "~50% residual threshold" is not a measured threshold

Appears in Appendix A.1 (Demonstrated links, L03), §4.1 (There is no single threshold,
because there is no single readout), §6.1 (The falsifying experiment) readout 3, and — most
consequentially — inside the boxed primary hypothesis in §4 (The primary hypothesis): *"if
it crosses the **~50% residual threshold**"*.

Bohers 2008 (PMID 18932004) is a **shRNA knockdown series in HeLa cells** — a transformed,
already-aneuploid line; the document does not say so. The gradient is
**8.5 / 10 / 14 / 58 / 77 %**. There is **no data point between 14% and 58%.** The paper's
sentence is "Aneuploidy was clearly detected in cells with residual BUBR1 below 50%" — but
the experiment can only locate the breakpoint somewhere in the open interval (14%, 58%).
"~50%" reads the upper edge of an unsampled gap as a measurement.

A second study in the same cell type points lower. Suijkerbuijk 2010 titrated doxycycline
against an inducible *BUB1B* shRNA in HeLa and found that **~13% residual BubR1 "had
virtually no impact on the fidelity of chromosome segregation"**, while ~6% caused
missegregation in the majority of cells. Different readouts (cumulative karyotype versus
per-division time-lapse) plausibly reconcile the two, but the honest summary is:

> The relationship is dosage-graded and there is a threshold. **Where it sits is not
> established.** Published knockdown data bound it somewhere between roughly 10% and 60%
> residual protein, and the two available HeLa series disagree about where in that range.
> Matsuura 2006 (PMID 16411201) independently concludes from human PCS-syndrome cells that
> "a >50% decrease in expression (or activity) of BubR1 is involved", which is convergent
> but is an association across patients, not a dose–response.

This is the single most quantitatively load-bearing number in the document and it should be
stated as a range with a named uncertainty, not a figure. Locating it is the
second-highest-value measurement in §6 (The experiments; see
`docs/track2-quantified-uncertainty.md` §5, What would raise the number, in order).

### WK-52 — SIRT2 overexpression extended median lifespan **in male mice only**

§4.3 (Why the NAD⁺/SIRT2 axis, and where it runs out) and Appendix G.1: "SIRT2
overexpression in **BubR1^H/H** hypomorphic animals increased median lifespan."

North 2014 full text: *"Lifespan extension was observed preferentially in male mice, with
a 123% increase in median lifespan… **There was no change in median lifespan for
females.**"* The pooled figure is +58% median, +21% maximal.

A sex-restricted effect stated as an unqualified effect is exactly the failure mode this
document criticises elsewhere. It also matters mechanistically: a null in half the animals
is evidence that the axis is not robustly load-bearing in the disease model.

The same passage carries a second fact the document should want: **"SIRT2tg/BubR1H/H mice
remained small."** Boosting the SIRT2 arm in the closest available animal model did **not**
correct the growth phenotype. That is direct support for the position in §7.2 (The weakest
link) that growth is not an achievable endpoint — better support than the citations §7.2
currently uses (see WK-70).

### WK-51 — The NMN result is a restoration in **aged wild-type** mice, at a high dose

§4.3 (Why the NAD⁺/SIRT2 axis, and where it runs out) and Appendix G.1: "Treating mice with
the NAD⁺ precursor **NMN increased BubR1 abundance in vivo**."

True, and the strongest single fact in the document. The full text adds three things it
does not say:

- The animals were **3-month-old and 30-month-old wild-type C57BL/6**, not BubR1-deficient.
- Dose was **500 mg/kg/day intraperitoneally for 7 days**.
- The headline effect is restorative: BubR1 in aged treated animals was "brought back to
  the levels of the 3-month-old animals".

Appendix C (The four objections, in full) already names the headroom objection and calls it
the strongest one. These three
details are what make it strong, and they belong next to the claim rather than one section
later. The dose in particular bears on condition 3 of the conjunction in
`docs/track2-quantified-uncertainty.md`.

### WK-03 — "Demonstrated" at the node overstates what the node evidence is

The direction table in §1 (The answer) answers *"does raising BubR1 correct the disease?"*
with **"AGREES. Signed, demonstrated."** What is demonstrated is narrower:

- Suijkerbuijk 2010's rescue readout is **a 2-fold increase in the nocodazole response**
  in patient cells, with near-full restoration shown against the **carrier parent** control
  (731X/WT), not a healthy control, in **one** patient line. Aneuploidy and missegregation
  were not readouts of the rescue.
- Baker 2013 is transgenic overexpression from development onward in mice that are not
  BubR1-deficient.

So: *the direction at the node is signed and demonstrated*; *correction of the disease* is
not. The second falsification criterion of §6.1 (The falsifying experiment) — "BubR1
increases but missegregation does not fall" — is precisely the gap that the foundational
human experiment leaves open. The document should say so, because a reviewer who reads
Suijkerbuijk 2010 will notice.

### WK-34 — Santaguida 2015 shows aneuploid cells' autophagy is **limited**, not that they depend on it

§3 (Direction of effect): "Aneuploid cells depend on autophagic degradation to manage their
protein burden (Santaguida 2015)."

The paper is titled *"Aneuploidy-induced cellular stresses **limit** autophagic
degradation"*: autophagosomal cargo accumulates within lysosomes, and aneuploid cells
respond by activating TFEB to raise autophagy/lysosomal gene expression. The reliance is
visible in the compensatory upregulation; the demonstrated finding is that capacity is
already compromised. The therapeutic direction is unaffected — a cell whose degradative
capacity is already saturated is *more* vulnerable to autophagy inhibition, not less — but
the sentence inverts the paper's result.

### WK-36 — Tang 2011 did not test hydroxychloroquine or metformin

§3 (Direction of effect): "Hydroxychloroquine and metformin are approved, cheap, extremely
well tolerated… Both are rejected here with `WRONG_DIRECTION` at the proteotoxic-stress
node, **on the same evidence that recommends them for the tumour.**"

Tang 2011 identified **AICAR, 17-AAG and chloroquine**. Its primary in-vivo-adjacent
demonstration is AICAR-induced p53-mediated apoptosis in **primary mouse embryonic
fibroblasts trisomic for chromosome 1, 13, 16 or 19**, plus AICAR and 17-AAG against
aneuploid human cancer lines. Neither hydroxychloroquine nor metformin appears.
17-AAG *is* tanespimycin, so that row is exact. The other two are class substitutions:
HCQ for chloroquine, metformin for AICAR across AMPK activation.

`SOURCES.md` discloses this ("direction carried across the AMPK class"). The hypothesis
document does not, and "the same evidence" is the phrase that makes it a claim rather than
an inference.

### WK-65 — The 21% malignancy figure is the wrong number, and better ones are already in our own reference list

§7.1 (What the model could not express) uses "a 21% malignancy rate in early case series
(Jacquemont 2002)". The paper says
exactly that — 3 of 14 (Jacquemont et al., *Am J Med Genet* 2002, PMID 11932988) — so the *reporting* holds. As a **rate for this
child** it does not:

- The 14 cases are **clinically ascertained MVA from 1992–2002, before *BUB1B* was
  implicated** (Hanks 2004). It is not a *BUB1B*-genotyped denominator.
- The malignancy that prompted the report — the index patient's ALL — is inside the
  numerator.
- Two papers **already cited in this document** give higher, genotype-anchored figures:
  Suijkerbuijk 2010's introduction states **"37% of patients develop cancers including
  rhabdomyosarcoma, Wilms tumor and leukemia, mostly within the first 3 years of life"**;
  García-Castillo 2008 (PMID 18548531) reports **Wilms tumour in 7 of 8** patients with
  monoallelic *BUB1B* mutations. Yost 2017 (PMID 28553959) adds that biallelic *BUB1B* and
  *TRIP13* carry high embryonal-tumour risk whereas *CEP57*-MVA does not.

Using 21% understates the case for the one intervention §7.1 (What the model could not
express) says has near-term value.
Fixing this and fixing W-3 are the same edit.

### WK-66 — Scott 2006 says there is no clear evidence that Wilms surveillance reduces mortality

§7.1 (What the model could not express) calls surveillance "the highest expected-value
intervention available today" and "if
one thing in this document has near-term value for a real child, it is this paragraph".

Scott 2006's own results section: *"The potential risks and benefits of Wilms tumour
surveillance are finely balanced and **there is no clear evidence that screening reduces
mortality or morbidity**. Prospective evidence-based data … would be difficult and costly
to generate and are unlikely to become available in the foreseeable future."* The
recommendations are explicitly **pragmatic**, not evidence-based.

The hedge in §7.1 (What the model could not express) — "It changes *stage at diagnosis*" —
is the right hedge and is what the
evidence supports. The superlative sitting next to it is not sourced, and the source it is
sitting next to argues against it. Keep the recommendation; surface the caveat.

### WK-70 — "Structural at birth" is not established by any of the three citations

§7.2 (The weakest link): "Growth restriction and microcephaly in MVA are established in
utero and are structural at birth (Hanks 2004; Bohers 2008; García-Castillo 2008)."

- **Bohers 2008** supports the growth half: the phenotype includes "severe **intrauterine**
  growth retardation".
- **Hanks 2004** lists "growth retardation, microcephaly" without onset or trajectory.
- **García-Castillo 2008** is the weakest fit of the three. Its two index patients are
  ***BUB1B*-mutation-negative** and "exhibited growth retardation … **but not
  microcephaly**". It is a clinical-subtype delineation, not a natal-anthropometry study.

None of the three reports head circumference at birth or its post-natal trajectory. That
microcephaly is structural and irreversible after birth is an **inference** — a reasonable
one, and the rule in §0 (Scope, and what the citations guarantee) is that inferences are
labelled.

This one matters because §7.2 (The weakest link) uses it to *reduce* the damage from L12
being weak: "L12
being weak weakens a phenotype the therapy was never going to reach." An unlabelled
inference doing load-bearing work in a self-criticism section is the worst place for it.
**Better evidence for the same conclusion exists and is already in the reference list:**
SIRT2 overexpression in BubR1^H/H mice extended male lifespan and left the animals small
(North 2014). That is a direct experimental result showing the intervention arm does not
rescue growth.

### WK-29 — "There is no catalytic activity to agonise" is stronger than the pseudokinase paper

§2 (The mechanism), citing Suijkerbuijk 2012 (PMID 22698286). The abstract:
*"Although BUBR1 evolved to a typical pseudokinase in some vertebrates, **it retained the
catalytic triad in humans**. However, we show that putative catalysis by human BUBR1 is
**dispensable** for error-free chromosome segregation."* The paper's own closing phrase is "an unusual, **triad-containing**
pseudokinase".

Dispensable is not absent. The therapeutic conclusion — do not look for an activator of
BubR1 enzyme function — survives, for a slightly different reason: not that there is
nothing to activate, but that activating it would not help.

There is a second consequence the document misses. Suijkerbuijk 2010 proposes that "the
kinase domain is involved in regulating BUBR1 protein stability in a non-catalytic
fashion", and shows the MVA kinase-domain mutants are HSP90 clients cleared by the
proteasome. That opens a **stability/chaperone-directed** pharmacology at the target node
which the document does not consider — and it flags a second, more direct objection to
tanespimycin than the one recorded (see §4 below).

### WK-45 — Baker 2011 is genetic ablation, not a senolytic drug

Appendix B (The therapeutic classes considered) and §6.2 (The confirming experiment) refer
to "Baker 2011's senolytic result". Baker 2011 (PMID 22048312) uses
**INK-ATTAC**, a transgene for inducible elimination of p16^Ink4a-positive cells on
administration of a dimerizer. It is drug-*triggered* genetic ablation, not pharmacological
senolysis, and the senolytic pharmacophore literature (Zhu 2015, PMID 25754370) is a
separate line. The use of it in §6.2 (The confirming experiment) as "the cautionary example
of an intervention that looked good on healthspan while the cancer question sat unasked" is
correct and unaffected — the
paper reports no tumour endpoint — but "senolytic result" should read "genetic
senescent-cell clearance".

### WK-59 — Pirinen 2020 is open-label and the controls improved too

Appendix G.2 (Direction at the molecule, agent by agent), nicotinic acid row: "raises human
NAD⁺ several-fold with functional benefit in mitochondrial myopathy (Pirinen 2020)".

NAD⁺ rise: **holds** — "Blood NAD+ increased in all subjects, up to 8-fold", at
750–1,000 mg/day. "Functional benefit" is where it slips. The design is **open-label with
no placebo arm**; patients and their matched controls both received niacin; and the result
is *"muscle strength and mitochondrial biogenesis increased **in all subjects**"* — controls
included. There is no controlled comparison from which disease-specific benefit can be read.
The abstract also reports **"Some patients showed anemia tendency"**, which the row's
"no sirtuin-inhibiting liability" framing does not balance.

### WK-06 / WK-38 — "Truncating alleles yield no transcript" is established for two of five

Appendix A.1 (Demonstrated links, L01) and Appendix B (The therapeutic classes considered,
ataluren row, where it carries the mechanistic case against readthrough agents).

Suijkerbuijk 2010 full text: mRNA from **386X** could not be detected; absence of **731X**
transcript was **inferred** from a 2-fold reduction in the 731X/WT carrier; and — the
document's phrasing does not admit this — *"**It is, however, unclear what the status of
194X, 483X or 753X is in patient cells.**"* Matsuura 2006 adds that transcripts from
1833delT and the splice-site allele were "significantly reduced, **probably** due to
nonsense-mediated mRNA decay".

The NMD objection to ataluren is sound for the alleles where it was measured and is an
extrapolation for the rest. Since ataluren is rejected on approval status independently,
nothing turns on it — but the sentence claims more than two Northern blots can carry.

### The remaining weakenings, briefly

- **WK-16** — L09 says aneuploidy imposes "proteotoxic **and metabolic** stress". R21–R24
  demonstrate proteotoxic, protein-folding and lysosomal/autophagic stress. None of the four
  demonstrates metabolic or energy stress; that comes from AICAR's selectivity in Tang 2011,
  which is a different paper and a different argument. Also: **Oromendia 2012 is yeast**, and
  is the first citation on the row.
- **WK-18 / WK-27 / WK-44** — Santaguida 2017's immune clearance is demonstrated as
  pro-inflammatory signalling plus **in-vitro NK-cell co-culture**. Appendix B (The
  therapeutic classes considered) and §7.1 (What the model could not express) extend it
  to "a child whose aneuploid cells are partly held in check by immunity" — an in-vivo human
  claim, used as an objection to sirolimus. It is a reasonable inference and is not labelled
  as one.
- **WK-30** — §3 (Direction of effect) says reversine "abrogates the checkpoint and causes
  missegregation".
  Santaguida 2010 reports **dose-dependent** SAC inhibition and failure of attachment error
  correction; it does not report a missegregation rate. The direction is not in doubt; the
  verb is one step ahead of the measurement.
- **WK-33** — the sentence in §3 (Direction of effect), "Each kills cells by exploiting a
  vulnerability that arises *from being aneuploid*", is Tang 2011's proposal ("Our results
  **suggest**…"), demonstrated for AICAR in trisomic MEFs and for AICAR+17-AAG in aneuploid
  human cancer lines. The model system — mouse trisomic MEFs — is not stated, and it is a
  good one for MVA, so stating it helps.
- **WK-73** — §7.3 (One child, and the grade of this work): "the largest human dataset
  behind any of this is a handful of MVA patient cell lines". True for the *mechanistic*
  claims. Not true for the phenotype and cancer-risk claims, which rest on Hanks 2004
  (5 families), Matsuura 2006 (7 families) and García-Castillo 2008 (35 patients reviewed,
  21 *BUB1B*-screened).

---

## 3. The full ledger

`R##` refers to the reference list in `knowledge/literature/bub1b/SOURCES.md`.

### §1 — The answer

| ID | Claim as written | Source | Verdict | What the source says instead |
|---|---|---|---|---|
| CL-01 | Ectopic BubR1 restored checkpoint activity in MVA patient cells | R03 | HOLDS | The paper's own words. Readout is a 2-fold rise in nocodazole response; near-full restoration shown in one line (731X/Y155C) against the carrier parent. |
| CL-02 | Sustained elevated BubR1 in mice reduced aneuploidy and tumorigenesis and extended healthy lifespan | R16 | HOLDS | — |
| CL-03 | Node verdict: raising BubR1 corrects the disease — "AGREES. Signed, demonstrated" | R03, R16 | **WEAKENED** | See WK-03. Direction demonstrated; disease correction not. |
| CL-04 | NMN is the one agent with a measured, agreeing direction, and is not an approved medicine | R19 | HOLDS | — |
| CL-05 | Nicotinamide is an inhibitor of the sirtuin the mechanism depends on | R36 | **WITHDRAWN** | See W-1. Bitterman measured yeast Sir2 and human SIRT1. |

### Appendix A.1 — Demonstrated links

| ID | Claim as written | Source | Verdict | What the source says instead |
|---|---|---|---|---|
| CL-06 | L01: absence of transcript from truncating mutants | R03 | **WEAKENED** | See WK-06. Measured for 386X, inferred for 731X, "unclear" for 194X/483X/753X. |
| CL-07 | L01: intergenic substitution 44 kb upstream; single-bp edit into cultured human cells produced reduced transcript, increased PCS and MVA | R05 | HOLDS | Verbatim ("hallmarks of the syndrome"). |
| CL-08 | L02: "Low protein abundance is the direct result of the absence of transcripts … combined with high protein turnover of missense mutants" | R03 | HOLDS | Direct quotation. |
| CL-09 | L03: MVA patient cell lines show impaired checkpoint, alignment defects and low BubR1 | R03 | HOLDS | — |
| CL-10 | L03: residual-BubR1 gradient 8.5 / 10 / 14 / 58 / 77 %, PCS across the range | R04 | HOLDS | Exact. Model system is **HeLa shRNA**; not stated in the document. |
| CL-11 | L03: "aneuploidy appearing below ~50% residual protein" | R04 | **WEAKENED** | See WK-11. No data point between 14% and 58%; Suijkerbuijk 2010's HeLa titration puts ~13% at "virtually no impact". |
| CL-12 | L04: the checkpoint complex restrains APC/C-Cdc20; failure permits premature anaphase | R13, R14 | HOLDS | — |
| CL-13 | L05: BubR1 stabilises KT–MT attachments by recruiting PP2A-B56, a route independent of the checkpoint | R09, R11, R12 | HOLDS | Suijkerbuijk 2012 KARD–PP2A-B56α; Kruse 2013 direct B56 binding; Lampson & Kapoor via Aurora regulation. |
| CL-14 | L06–L08: checkpoint failure → PCS → whole-chromosome missegregation; destabilised attachments do the same | R03, R04, R06 | HOLDS | — |
| CL-15 | L09: aneuploidy imposes proteotoxic stress | R21–R24 | HOLDS | R21 is **yeast**; R22, R23 human. |
| CL-16 | L09: "…and metabolic stress" | R21–R24 | **WEAKENED** | None of the four demonstrates metabolic/energy stress. |
| CL-17 | L10: missegregation triggers p53/p16 arrest and senescence | R26, R17 | HOLDS | Note R17 links **BubR1 insufficiency** to Cdkn2a activation in skeletal muscle and fat, not missegregation directly. |
| CL-18 | L10: "…and immune clearance" | R25 | **WEAKENED** | Demonstrated as pro-inflammatory signalling plus in-vitro NK co-culture. |
| CL-19 | L13: 3 of 14 reported MVA cases (21%) developed a malignancy — RMS, ALL, nephroblastoma | R08 | HOLDS | Faithful to the paper. See WK-65 for why it is the wrong denominator for this child. |
| CL-20 | L13: low-rate missegregation promotes tumours, very high-rate suppresses them | R28, R29 | HOLDS | Silk 2013 states it in those terms. |

### Appendix A.2 — The two inferred links

| ID | Claim as written | Source | Verdict | What the source says instead |
|---|---|---|---|---|
| CL-21 | L11 is not directly demonstrated | — | HOLDS | Correctly labelled inference. Thompson & Compton's own abstract flags tissue-context exceptions, which supports the flag. |
| CL-22 | L12: Baker 2004 establishes BubR1 insufficiency produces dwarfism and progeroid phenotypes | R15 | HOLDS | Mouse; the document says so. |
| CL-23 | L12: BubR1 insufficiency also impairs ciliogenesis | R20 | HOLDS | Human cells plus **medaka fish** morpholino; organism not stated. |

### Appendix A.2 — The compensatory node

| ID | Claim as written | Source | Verdict | What the source says instead |
|---|---|---|---|---|
| CL-24 | In BubR1-hypomorphic mice, deleting p16 attenuated the progeroid phenotype | R17 | HOLDS | Restricted to skeletal muscle and fat. The document omits the opposing arm: **p19Arf inactivation exacerbated** senescence and ageing — which strengthens the compensatory-node reading in Appendix A.2 (The two inferred links, and the compensatory node). |
| CL-25 | Clearing p16-positive senescent cells delayed age-related deterioration, in the BubR1 progeroid model | R18 | HOLDS | Exact. No lifespan extension claimed, and the document does not claim one. |
| CL-26 | The same p53 arm suppresses aneuploidy-driven tumorigenesis | R27 | HOLDS | Mouse; two independent aneuploidy models. |
| CL-27 | Immune clearance removes missegregating cells | R25 | **WEAKENED** | As CL-18. |

### §2 — The mechanism: non-checkpoint roles

| ID | Claim as written | Source | Verdict | What the source says instead |
|---|---|---|---|---|
| CL-28 | MVA patient cells show both an impaired checkpoint and chromosome alignment defects | R03 | HOLDS | 2–3 fold increase in misalignment under MG132. |
| CL-29 | BubR1 is a pseudokinase; there is no catalytic activity to agonise | R10 | **WEAKENED** | See WK-29. Human BUBR1 "retained the catalytic triad"; catalysis is *dispensable*. |

### §3 — Direction of effect

| ID | Claim as written | Source | Verdict | What the source says instead |
|---|---|---|---|---|
| CL-30 | Reversine abrogates the checkpoint and causes missegregation | R32 | **WEAKENED** | Dose-dependent SAC inhibition and loss of attachment error correction; missegregation not measured. |
| CL-31 | The MPS1-inhibitor class is in active anticancer development | R33 | HOLDS | "developed and approved to begin clinical trials for advanced nonhaematologic malignancies". |
| CL-32 | Tang 2011 found three aneuploidy-selective compounds: AICAR, 17-AAG, chloroquine | R31 | HOLDS | Exact. |
| CL-33 | Each kills cells by exploiting a vulnerability that arises from being aneuploid | R31 | **WEAKENED** | The paper's proposal; demonstrated for AICAR in trisomic MEFs and AICAR+17-AAG in aneuploid human cancer lines. |
| CL-34 | Aneuploid cells depend on autophagic degradation to manage their protein burden | R23 | **WEAKENED** | See WK-34. The paper shows autophagic degradation is *limited*. |
| CL-35 | Aneuploid cells have impaired HSP90-dependent folding | R22 | HOLDS | Title-level. |
| CL-36 | HCQ and metformin rejected "on the same evidence that recommends them for the tumour" | R31 | **WEAKENED** | See WK-36. Tang tested chloroquine and AICAR. |
| CL-37 | Baker 2013: raising BubR1 reduced aneuploidy **and** tumorigenesis, even against oncogenic Ras | R16 | HOLDS | Exact. |
| CL-38 | Truncating *BUB1B* alleles yield no transcript, so NMD removes the substrate | R03 | **WEAKENED** | As CL-06. |
| CL-39 | Ataluren's EU authorisation was not renewed (EC decision effective 2025) | R53 | HOLDS | EMA record: EC decision **28 March 2025**, status "Expired"; CHMP recommended non-renewal because "the effectiveness of the medicine has not been confirmed". |
| CL-40 | BubR1 is stabilised by K250 acetylation | R35 | HOLDS | K250Q degradation inhibited; K250R unstable. |
| CL-41 | HDAC inhibitors induce premature sister-chromatid separation and override the SAC | R34 | HOLDS | Title and abstract, at class level. |
| CL-42 | Dasatinib with quercetin as senolytic | R46 | HOLDS | Dasatinib alone eliminated senescent human fat-cell progenitors; the combination was needed for MEFs. |
| CL-43 | Antioxidants accelerated tumour progression in mice | R52 | HOLDS | NAC and vitamin E, B-RAF and K-RAS lung models; mechanism is ROS–p53 disruption. |
| CL-44 | Sirolimus is immunosuppressive "in a child whose aneuploid cells are partly held in check by immunity" | R25 | **WEAKENED** | As CL-18; unlabelled inference from in-vitro NK data to in-vivo human immunosurveillance. |
| CL-45 | Baker 2011 is the motivating "senolytic" evidence, in the exact BubR1 model | R18 | **WEAKENED** | See WK-45. INK-ATTAC genetic ablation, not pharmacological senolysis. |
| CL-77 | Vorinostat specifically was among the HDAC inhibitors tested by Magnaghi-Jaulin | R34 | **UNVERIFIED** | *Cancer Res* 2007 is not open access; the abstract says "HDACI treatment" without naming agents. The class-level claim (CL-41) stands; the agent-level one is unchecked. |

### §4 — The primary hypothesis

| ID | Claim as written | Source | Verdict | What the source says instead |
|---|---|---|---|---|
| CL-46 | The boxed hypothesis: "if it crosses the ~50% residual threshold" | R04 | **WEAKENED** | As CL-11. This is the most load-bearing instance. |
| CL-47 | The rescue is "a rescue experiment in the actual human disease" | R03 | HOLDS | Patient-derived SV40-transformed fibroblasts and EBV-transformed lymphoblastoids. |
| CL-48 | Baker 2013 corrected both checkpoint impairment and microtubule–kinetochore attachment defects | R16 | HOLDS | Verbatim. |
| CL-49 | Lowering BubR1 in vivo reproduces the phenotype | R15 | HOLDS | — |
| CL-50 | SIRT2 holds K668 deacetylated NAD⁺-dependently; CBP acetylates it; the balance determines BubR1 turnover | R19 | HOLDS | Full text: K668Q raises ubiquitination and degradation rate; K668R lowers both. |
| CL-51 | Treating mice with NMN increased BubR1 abundance in vivo | R19 | **WEAKENED** | See WK-51. 500 mg/kg/day i.p., 7 days, wild-type C57BL/6 aged 3 and 30 months; effect framed as restoration in aged animals. |
| CL-52 | SIRT2 overexpression in BubR1^H/H animals increased median lifespan | R19 | **WEAKENED** | See WK-52. +123% median in males; **no change in females**. |
| CL-53 | North 2014 frames the axis as an explanation for the age-related decline in BubR1 | R19 | HOLDS | The paper's framing exactly. Supports the headroom concession in Appendix C (The four objections, in full). |
| CL-54 | "Half" of the Suijkerbuijk classes are qualitative defects more protein will not fix | R03 | **WITHDRAWN** | See W-2. One characterised allele (Y155C) is a genuine separation-of-function defect; R550Q rescued both functions; Q921H is wild-type-like. |
| CL-55 | Mad2 overexpression is tumorigenic in mice | R30 | HOLDS | — |
| CL-56 | NMN's effect on aneuploidy and cancer risk has never been assessed | R19 | HOLDS | Confirmed against the full text: no aneuploidy or tumour endpoint appears. |
| CL-57 | Nicotinamide riboside raises human NAD⁺ | R37, R38 | HOLDS | Trammell's 2.7-fold is an **n=1 pilot**; the dose-ranging trial (100/300/1,000 mg) is the substantive result. Martens is middle-aged and older adults. |
| CL-58 | Nicotinic acid raises human NAD⁺ several-fold | R39 | HOLDS | Up to 8-fold in blood at 750–1,000 mg/day. |
| CL-59 | "…with functional benefit in mitochondrial myopathy" | R39 | **WEAKENED** | See WK-59. Open-label, no placebo; strength rose "in all subjects", controls included. |
| CL-60 | Nicotinic acid has "no sirtuin-inhibiting liability" | — | **UNVERIFIED** | No citation attached. Mechanistically plausible (nicotinic acid is not the sirtuin reaction product) but not sourced, and Avalos 2005 engineered a Sir2 variant that *is* inhibited by nicotinic acid — so the claim is about wild-type enzymes and should say so. |
| CL-61 | ENDIT randomised 552 first-degree relatives, children included, to 1.2 g/m²/day or placebo for five years, and reported no effect on growth in children | R41 | HOLDS | "Nicotinamide treatment did not affect growth in children." 552 ICA-positive relatives; 18 European countries, Canada and the USA. |

### §5 — What the pipeline produced

| ID | Claim as written | Source | Verdict | What the source says instead |
|---|---|---|---|---|
| CL-62 | A phase 3 randomised trial with a cancer-incidence endpoint found fewer cancers, in adults with UV-driven skin cancer | R40 | HOLDS | n=386, 500 mg b.d. 12 months, 23% lower rate of new non-melanoma skin cancers (95% CI 4–38, P=0.02). Also: "no evidence of benefit after nicotinamide was discontinued". |
| CL-63 | For nicotinic acid the oncogenic-risk gate is answerable **from existing data** — large cardiovascular outcome trials with cancer-incidence data already collected | — | **UNVERIFIED** | No citation attached, and the largest such trial does not report it. HPS2-THRIVE (Landray et al., *N Engl J Med* 2014;371:203-12, PMID 25014686; n=25,673, median 3.9 y) reports serious adverse events by body system — diabetes, GI, musculoskeletal, skin, infection, bleeding — with **no cancer-incidence endpoint in the primary publication**. Malignancies may sit in unreported SAE capture; that is a data-request, not existing published data. See §5 of `docs/track2-quantified-uncertainty.md` (What would raise the number, in order). |

### §6 — The experiments

| ID | Claim as written | Source | Verdict | What the source says instead |
|---|---|---|---|---|
| CL-64 | The Suijkerbuijk 2010 functional replacement assay can establish whether residual protein is quantitatively or qualitatively defective | R03 | HOLDS | shRNA depletion plus shRNA-insensitive LAP-BUBR1 mutant re-expression, scored on nocodazole response and chromosome alignment. Exactly the assay described. |

### §7 — Limitations

| ID | Claim as written | Source | Verdict | What the source says instead |
|---|---|---|---|---|
| CL-65 | "a 21% malignancy rate in early case series" | R08 | **WEAKENED** | See WK-65. Pre-molecular clinical ascertainment; index case in the numerator; higher genotype-anchored figures exist in our own reference list. |
| CL-66 | "established protocols for Wilms tumour surveillance in at-risk children" (Scott 2006) | R49 | **WEAKENED** | See WK-66. The paper states there is no clear evidence screening reduces mortality or morbidity, and calls its recommendations pragmatic. Its operative rule (>5% risk → 3–4 monthly ultrasound) is usable and should be quoted. |
| CL-67 | "updated in Kalish et al. 2024" | R54 | **WITHDRAWN** | See W-3. MVA, *BUB1B* and PCS syndrome appear nowhere in that guideline. |
| CL-68 | The paediatric dasatinib trial reported bone growth and development events in 4% of patients | R45 | HOLDS | Verbatim. n=113 CML-CP. |
| CL-69 | Gentamicin's cumulative ototoxicity and nephrotoxicity make chronic dosing implausible | — | **UNVERIFIED** | No citation attached. Uncontroversial pharmacology, but the rule in §0 (Scope, and what the citations guarantee) of the hypothesis document requires either a citation or an inference label. |
| CL-70 | Growth restriction and microcephaly are established in utero and structural at birth | R01, R04, R07 | **WEAKENED** | See WK-70. IUGR supported by R04. Nothing in the three establishes microcephaly at birth or its irreversibility; R07's index patients are *BUB1B*-negative and **lacked** microcephaly. |
| CL-71 | L03 is the best-evidenced link in the chain | R03, R04 | HOLDS | Consistent with the sources, with the CL-11 caveat about the threshold's precision. |
| CL-72 | "There has never been a clinical trial in MVA" | — | HOLDS | Verified 2026-08-28 against the ClinicalTrials.gov v2 API: `mosaic variegated aneuploidy` → 0 studies; `variegated aneuploidy` → 0; `BUB1B` → 0 (control queries `aneuploidy` → 313, `Fanconi anemia` → 169). ClinicalTrials.gov is not a complete global registry, so this is "no registered US/international trial", not a world-historical negative. |
| CL-73 | "The largest human dataset behind any of this is a handful of MVA patient cell lines" | — | **WEAKENED** | True of the mechanistic claims; not of the phenotype and cancer-risk claims, which rest on multi-family case series (R01, R06, R07). |

### Appendix B — The therapeutic classes considered, and §8 — Scalability: remaining rows

| ID | Claim as written | Source | Verdict | What the source says instead |
|---|---|---|---|---|
| CL-74 | *CEP57* causes MVA2 through a centrosomal mechanism | R50 | HOLDS | "CEP57 is a centrosomal protein and is involved in nucleating and stabilizing microtubules." |
| CL-75 | *TRIP13* causes a third MVA form through checkpoint silencing | R51 | HOLDS | Yost 2017 reports substantial SAC impairment in TRIP13-null patient cells; "checkpoint silencing" describes TRIP13's normal MAD2-disassembly function, which is the correct gloss. |
| CL-76 | proTAME and apcin are "tool compounds never given to humans" | — | **UNVERIFIED** | No citation attached. Consistent with the absence of any registered trial, but asserted rather than sourced. |

---

## 4. Two things I found that are not citation errors, and matter more than most of them

### 4.1 The document holds two opposite acetylation arguments and never reconciles them

- **Appendix B (The therapeutic classes considered), vorinostat row**: BubR1 is
  **stabilised by K250 acetylation** (Choi 2009,
  PMID 19407811) — the argument for HDAC inhibitors, which the document then defeats on the
  mitotic readout.
- **§4.3 (Why the NAD⁺/SIRT2 axis, and where it runs out) and Appendix G.1, the primary
  hypothesis**: BubR1 is **destabilised by K668 acetylation**; SIRT2
  deacetylation of K668 raises abundance (North 2014, full text: K668Q raises ubiquitination
  and degradation rate).

Two lysines, opposite signs on the same protein's stability. The primary hypothesis proposes
to shift the cell toward **more** sirtuin-driven deacetylation. Nothing in either paper
establishes that a global NAD⁺/sirtuin shift is site-selective for K668 and spares K250
(a PCAF site). If it is not, the intervention buys stability at one lysine and spends it at
another.

This is not a mis-citation — both claims are correctly sourced. It is an unexamined
mechanistic vulnerability sitting between two sections, and it is exactly the kind of thing
a reviewer who reads both citations will find. It belongs in Appendix C (The four
objections, in full) as a fourth objection.

### 4.2 Tanespimycin has a second, more direct wrong-direction objection — at the target node

Appendix B (The therapeutic classes considered) rejects tanespimycin (17-AAG) at the
**proteotoxic-stress node** as an
aneuploidy-selective antiproliferative. Suijkerbuijk 2010 (already cited, R03) supplies a
stronger objection at the **target node**: MVA kinase-domain substitution mutants are HSP90
clients. Fig. 3C — *"levels of both mutant proteins were severely decreased after 5h
treatments with the HSP90 inhibitor"*, and CHX + geldanamycin together "removed virtually
all mutant BUBR1 protein", while wild-type BubR1 showed "only minor changes".

So HSP90 inhibition **directly destroys the exact protein the primary hypothesis is trying
to raise**, selectively in the mutant. That is a `WRONG_DIRECTION` at `N3_bubr1`, not only
at `N9_stress`, and it is measured in the disease's own mutants rather than carried across
from a generic aneuploidy model.

The rejection outcome does not change. What this shows is a genuine limitation of the
architecture that §7.1 (What the model could not express) half-anticipates ("single-node
scoring cannot see cross-node harm"):
here the *stronger* evidence sat at a different node in a paper the catalogue already cites,
and single-node scoring routed around it. That is a better worked example of the limitation
than the sirolimus one currently in §7.1 (What the model could not express), because the
missed evidence is not hypothetical.

**A corollary worth stating:** the same finding opens a pharmacology the document does not
consider. If MVA substitution mutants are destabilised by misfolding and cleared by the
proteasome, then chaperone-supporting or proteasome-sparing approaches act at the target
node by a route entirely independent of NAD⁺/SIRT2. No approved agent obviously fills it,
which is why it is a note and not a row — but "the only handle on abundance is the
acetylation axis" is not established, and §4 (The primary hypothesis) currently reads as
though it were.

---

## 5. What I could not check

| Item | Why |
|---|---|
| Which specific HDAC inhibitors Magnaghi-Jaulin 2007 tested (CL-77) | *Cancer Res* 2007;67:6360-7 is not open access; abstract says "HDACI treatment" without naming agents. |
| Whether niacin cardiovascular outcome trials collected malignancy data (CL-63) | HPS2-THRIVE's primary publication reports SAEs by body system, not cancer incidence. Confirming would need trial-level SAE tables or a data request. |
| Nicotinic acid's sirtuin liability (CL-60) | No primary measurement located. |
| Gentamicin chronic-toxicity claim (CL-69), proTAME/apcin human exposure (CL-76) | No citation attached in the source document; nothing to check against. |

None of these is a reason to doubt the claim. They are reasons not to write it down as
verified.

---

## 6. What this exercise did not do

- It did not check **citation identity** — that each PMID resolves to the author, year,
  journal and volume written down. That is a separate pass with a separate owner.
- It did not re-derive the pipeline's scores, which are arithmetic over the curated tables
  and are covered by `tests/unit/test_track2_bub1b.py`.
- It read abstracts for every claim and full text where the claim turned on something an
  abstract cannot settle. A claim graded HOLDS on an abstract could still be weakened by a
  figure I did not see. Where that risk seemed material — the threshold, the mutation
  classes, the NMN experiment, the lifespan result, the surveillance guideline — I went to
  full text and it changed the verdict five times out of five.

That last observation is the honest summary of this whole document: **the errors were
where the full text was.**

---

## 7. Corrections queue for `docs/track2-hypothesis.md`

Exact edits, ordered by how badly they are needed. Nothing here has been applied — this
ledger does not modify the hypothesis document.

> **Applied since, noted 2026-08-29.** The sentence above was true on the day it was
> written and is not true now. Most of this queue has been worked through in
> `docs/track2-hypothesis.md`: thirteen applied, two overtaken, one partial. Each entry
> below now opens with a status marker; the correction text underneath every marker is the
> original and has not been rewritten. The counts, and the one item still worth acting on,
> are in the status note in §0.

### Must fix — the three withdrawals

**C-1 (§1, The answer; §4.4, The sign paradox at the molecule — W-1). SUPERSEDED — see the
W-1 CORRECTION above; APPLIED in that corrected form, in §1 (The answer), §4.4 (The sign
paradox at the molecule) and Appendix G.2 (Direction at the molecule, agent by agent) of
the hypothesis document.** The remedy below
was applied, and its second half was wrong: SIRT2 *has* been measured against nicotinamide
(Marcotte 2004, PMID 15301953, R58). The wording actually now in the document is the
corrected wording, not the wording quoted here. Kept for the record.

Original remedy: anywhere nicotinamide is described as inhibiting "the sirtuin
the mechanism depends on" on the authority of Bitterman et al., *J Biol Chem* 2002, PMID 12297502, change the evidence status
from measurement to class inference and add the supporting reference.

- §1 (The answer) table, molecule row: *"one of them — nicotinamide — is simultaneously an
  inhibitor of
  the sirtuin the mechanism depends on"* →
  *"one of them — nicotinamide — is a product inhibitor of sirtuins as a class, so it feeds
  and plausibly inhibits the same axis; SIRT2 specifically has never been measured against
  it."*
- Appendix G.2 of the hypothesis document (Direction at the molecule, agent by agent) —
  the agent-by-agent table that stood at §4.4 when this ledger was written —
  nicotinamide row: after "(Bitterman et al., *J Biol Chem* 2002;277:45099-107,
  Bitterman et al., *J Biol Chem* 2002, PMID 12297502)", add: "— measured against yeast Sir2 and human SIRT1, not SIRT2;
  extension to SIRT2 rests on the conserved nicotinamide-binding pocket 
  (Avalos et al., *Mol Cell* 2005, 17:855-68, PMID 15780941)."
- Add R56 to `SOURCES.md`: Avalos JL, Bever KM, Wolberger C (2005) Mechanism of sirtuin
  inhibition by nicotinamide: altering the NAD(+) cosubstrate specificity of a Sir2 enzyme.
  Avalos et al., *Mol Cell* 2005, 17(6):855-68, PMID 15780941.

**C-2 (Appendix C, The four objections, in full — second objection; W-2). APPLIED** — the
second objection in Appendix C now carries the five-allele split, R550Q, Q921H, Y155C, the
paper's own conclusion and the misfolding version of the objection. Replace the
objection and the answer:

> **"Suijkerbuijk 2010 describes two mutation classes. Half of them are qualitative defects
> that more protein will not fix."** Partly. The split is more favourable than "half": of
> the substitutions characterised, five kinase-domain alleles are abundance-lowering, R550Q
> rescued both BubR1 functions and Q921H is wild-type-like, and only **Y155C** is a genuine
> separation-of-function allele. The paper's own conclusion is that CIN in *BUB1B*-MVA "is a
> result of low BUBR1 protein abundance". The hypothesis is still explicitly
> **genotype-conditional** (ASSUMPTION-DRUG-09, §7.3) — the class has to be determined per
> patient — but the eligible fraction is larger than a coin flip. The sharper version of the
> objection is different: the abundance-class mutants are destabilised by **misfolding**
> (HSP90-dependent, proteasome-cleared), and whether K668 deacetylation can counteract
> misfolding-driven turnover has never been tested.

**C-3 (§7.1, What the model could not express; W-3). APPLIED** — the Kalish citation is
gone from the hypothesis document; §7.1 and Appendix E.1 carry the >5% rule with all three
tumour rates; and `SOURCES.md` records R54 under "considered and left out".
In the tumour-surveillance bullet,
delete "updated in Kalish et al.,
*Clin Cancer Res* 2024;30:5260-5269, PMID 39320341" — that guideline does not mention MVA,
*BUB1B* or PCS syndrome. Replace the surveillance sentence with the operative rule and the
correct risk figure:

> Tumour surveillance — Scott et al. recommend surveillance for any child at **>5% risk of
> Wilms tumour**, by renal ultrasonography every 3–4 months (*Arch Dis Child* 2006;91:995-9,
> Scott et al., *Arch Dis Child* 2006, PMID 16857697), while stating plainly that the risks and benefits are finely balanced and
> that there is no clear evidence screening reduces mortality or morbidity. *BUB1B*-MVA
> clears the 5% rule on any published series: Jacquemont 2002 reports malignancy in 3 of 14
> clinically ascertained cases (Jacquemont et al., *Am J Med Genet* 2002, PMID 11932988); Suijkerbuijk 2010 puts it at 37%; and
> García-Castillo 2008 reports Wilms tumour in 7 of 8 patients with monoallelic *BUB1B*
> mutations (García-Castillo et al., *Am J Med Genet A* 2008, PMID 18548531). It does not move any node on this chain. It changes *stage at
> diagnosis*.

Update `SOURCES.md` R54 to "considered and left out" or delete it.

### Should fix — weakenings that change a stated number or a stated certainty

| # | Where | Change |
|---|---|---|
| C-4 — **OVERTAKEN** | §4 (The primary hypothesis) boxed hypothesis, Appendix A.1 (Demonstrated links) L03, §4.1 (There is no single threshold, because there is no single readout), §6.1 (The falsifying experiment) readout 3 | Replace "~50% residual threshold" with a range and a named uncertainty. Suggested: *"if it crosses the residual-BubR1 threshold — which published knockdown series place somewhere between roughly 10% and 60%, and do not agree on"*. Add to Appendix A.1 L03 that Bohers is a **HeLa shRNA series with no data point between 14% and 58%**, and that Suijkerbuijk 2010's own HeLa titration found ~13% residual had "virtually no impact on the fidelity of chromosome segregation". **Overtaken, 2026-08-29:** the phrase "~50% residual threshold" no longer appears anywhere in the hypothesis document. §4.1 (There is no single threshold, because there is no single readout) replaced it with a dose–response whose position is stated as unmeasured; Appendix A.1 L03 now names the shRNA HeLa series and its transformed, already-aneuploid background; and Appendix A.3 (The dose–response, readout by readout) carries both the unsampled gap and the ~13% titration. |
| C-5 — **APPLIED** | §4.3 (Why the NAD⁺/SIRT2 axis, and where it runs out) and Appendix G.1 | "SIRT2 overexpression in **BubR1^H/H** hypomorphic animals increased median lifespan" → *"…increased median lifespan by 123% in males, with no change in females; the animals remained small."* |
| C-6 — **APPLIED** | §4.3 (Why the NAD⁺/SIRT2 axis, and where it runs out) and Appendix G.1 | After "NMN increased BubR1 abundance in vivo", add *"(500 mg/kg/day i.p. for 7 days, in wild-type mice aged 3 and 30 months; the effect is framed as restoration of an age-related decline)"* — which is the headroom objection of Appendix C (The four objections, in full) made concrete at the point of the claim. |
| C-7 — **APPLIED** | §7.2 (The weakest link) | "established in utero and are structural at birth (Hanks 2004; Bohers 2008; García-Castillo 2008)" → keep Bohers for intrauterine growth restriction; drop García-Castillo (its index patients are *BUB1B*-negative and lacked microcephaly); label the irreversibility of microcephaly as **inference**. Add the better evidence: SIRT2 overexpression in BubR1^H/H mice left the animals small (North 2014). |
| C-8 — **PARTIAL** | §3 (Direction of effect) | "Hydroxychloroquine and metformin … on the same evidence that recommends them for the tumour" → say that Tang 2011 tested **chloroquine and AICAR** and that the approved analogues carry the direction across by class. Add that the primary demonstration is in **trisomic mouse embryonic fibroblasts** — a point in the document's favour, since constitutional trisomy is the better model. **Partial, 2026-08-29:** the first half is applied — §3 rejects all three agents "on the evidence that recommends their analogues for the tumour arm", and Appendix B records hydroxychloroquine standing for chloroquine and metformin for AICAR. The second half is not: the model system is still unstated, and "trisomic" appears nowhere in the hypothesis document. |
| C-9 — **APPLIED** | §3 (Direction of effect) | "Aneuploid cells depend on autophagic degradation to manage their protein burden (Santaguida 2015)" → *"Aneuploid cells accumulate autophagosomal cargo in lysosomes and compensate by activating TFEB; their autophagic degradation capacity is already limited (Santaguida 2015)"* — which is what the paper shows and a stronger argument for the same conclusion. |
| C-10 — **APPLIED** | §2 (The mechanism) | "There is no catalytic activity to agonise" → *"Human BubR1 retains the catalytic triad but its catalysis is dispensable for error-free segregation (Suijkerbuijk 2012), so there is no activity worth agonising."* |
| C-11 — **OVERTAKEN** | §5.2 (The gate that should not have been shut) and §7.1 (What the model could not express) | The claim that nicotinic acid's cancer question is answerable "from existing data" needs a citation and a caveat: HPS2-THRIVE (Landray et al., *N Engl J Med* 2014, PMID 25014686; n=25,673, median 3.9 y) reports serious adverse events by body system and **no cancer endpoint**. It is a data request against existing safety capture, not a re-read of published tables. **Overtaken, 2026-08-29:** §5.2 (The gate that should not have been shut) now cites HPS2-THRIVE and states the serious-adverse-event caveat this asks for, but reaches the opposite conclusion — the per-arm malignancy counts are posted and public, so it is no longer a data request. See the note on CL-63 in §0. |
| C-12 — **APPLIED** | Appendix B (The therapeutic classes considered), §6.2 (The confirming experiment) | "senolytic result" → "genetic senescent-cell clearance" for Baker 2011 (INK-ATTAC, not a senolytic drug). |
| C-13 — **APPLIED** | Appendix G.2 (Direction at the molecule, agent by agent) | Nicotinic acid row: qualify "functional benefit in mitochondrial myopathy" as **open-label, no placebo arm, with strength rising in treated controls as well as patients**, and note the reported anaemia tendency. |
| C-14 — **APPLIED** | §1 (The answer) direction table | "AGREES. Signed, demonstrated" → note that what is demonstrated is the direction at the node; the human rescue readout is checkpoint response, not missegregation. |
| C-15 — **APPLIED** | Appendix C (The four objections, in full) | Add a fourth objection: BubR1 is **stabilised** by K250 acetylation (Choi 2009) and **destabilised** by K668 acetylation (North 2014). Nothing establishes that a global NAD⁺/sirtuin shift is site-selective for K668. See ledger §4.1. |
| C-16 — **APPLIED** | §7.1 (What the model could not express) | Replace or supplement the sirolimus example of cross-node blindness with the tanespimycin one: HSP90 inhibition directly destroys mutant BubR1 in MVA patient-derived cells (Suijkerbuijk 2010 Fig. 3C), a `WRONG_DIRECTION` at `N3_bubr1` that single-node scoring routed around. See ledger §4.2. |

### Label, don't rewrite

Six claims carry no citation and no inference label, which the rule in §0 (Scope, and what
the citations guarantee) forbids: CL-60
(nicotinic acid sirtuin liability), CL-63 (existing CVOT cancer data), CL-69 (gentamicin
chronic toxicity), CL-76 (proTAME/apcin never given to humans), plus the in-vivo human
immunosurveillance inference at CL-18/27/44 and the "highest expected-value intervention"
judgement at CL-66. Each needs either a source or the word **inference**.

> **APPLIED, all six, noted 2026-08-29.** In the hypothesis document CL-60, CL-69, CL-76 and
> CL-66 are now labelled **inference** where they are stated; the immunosurveillance
> extension is labelled at Appendix A.1 L10, in the sirolimus row of Appendix B and in
> Appendix E.1; and CL-63 is sourced in §5.2, in the other direction, as the §0 status note
> records.

---

*No patient data was accessed in producing this ledger. Public literature only.*
