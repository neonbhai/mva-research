# Track 2 — A drug-repurposing hypothesis for BUB1B-related mosaic variegated aneuploidy

**Mosaic variegated aneuploidy syndrome 1** · OMIM 257300 · MONDO:0009759 ·
Orphanet:1052 · *BUB1B* (HGNC:1149, ENSG00000156970, UniProt O60566, 15q15.1)

---

## 0. What this is, and four things it is not

This is a **research hypothesis about a gene**, assembled from public literature and
run through this repository's mechanism and intervention stages. Every factual claim
below carries an author, year, journal and PMID; every citation was resolved against
PubMed before being written down. Where a claim has no publication behind it, it is
labelled **inference** and the label is not decoration.

It is **not medical advice.** Nothing here should reach a prescription pad. A
mechanism-grounded repurposing hypothesis is a proposal for pre-clinical work.

It is **not derived from patient data.** No genotype, no coordinates, no clinical
record were read, requested or inferred. Under ADR 0008 and GP-44 the patient
workspace is deliberately illegible to the agent that wrote this. That constraint
turned out to be productive rather than merely restrictive, and §7.3 explains why:
the honest form of this hypothesis is *genotype-conditional*, and stating the
condition explicitly is better science than quietly assuming it away.

It is **not evidence of efficacy in any individual.** See §8 on the n-of-1 problem.

---

## 1. The answer, up front

**The mechanism is a dosage disease.** *BUB1B* pairs a null or truncating allele with
a hypomorphic missense or regulatory allele, and the quantity that is deficient is
**BubR1 protein**. The disease is named after the spindle assembly checkpoint, but the
checkpoint is a consequence node; the protein is where the deviation is a number.

**The primary hypothesis** is therefore: *raise residual BubR1 protein abundance above
the threshold at which chromosome missegregation emerges, via the NAD⁺–SIRT2–BubR1
K668 deacetylation axis.*

**Direction-of-effect verdict, stated at two levels, because they differ:**

| Level | Verdict | Basis |
|---|---|---|
| At the **node** (does raising BubR1 correct the disease?) | **AGREES.** Signed, demonstrated, and — unusually — desirable in *both* soma and tumour. | Ectopic BubR1 restored checkpoint activity in MVA patient cells (Suijkerbuijk 2010); sustained elevated BubR1 in mice reduced aneuploidy **and** tumorigenesis and extended healthy lifespan (Baker 2013). |
| At the **molecule** (does any approved drug do it?) | **CANNOT DETERMINE, for every approved agent.** Not "no", not "yes". | The one agent with a measured, agreeing direction (NMN, North 2014) is not an approved medicine. Every approved NAD⁺-raising agent has an unmeasured effect on BubR1, and one of them — nicotinamide — is simultaneously an inhibitor of the sirtuin the mechanism depends on. |

The pipeline's output reflects that split exactly: of 13 real agents evaluated,
**1 was accepted and 12 rejected**; the accepted candidate carries
`directions_agree = None`; and the highest-scoring agent in the whole catalogue is one
of the rejected ones.

**The deliverable is not a prescription. It is a first experiment** (§6), and a
statement of what would falsify the hypothesis.

---

## 2. The mechanism, graded link by link

A chain presented at uniform confidence is a rhetorical device, not science
(ASSUMPTION-MECHANISM-01). The machine-readable chain is
`knowledge/literature/bub1b/mechanisms.tsv`: 12 nodes, 13 links, each with its own
tier, strength, `is_directly_demonstrated` flag and stated uncertainty. Citations map
to links in `knowledge/literature/bub1b/SOURCES.md`.

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

### 2.1 Established links

| Link | Claim | Evidence |
|---|---|---|
| L01 | Truncating and upstream-regulatory *BUB1B* alleles reduce full-length transcript | Suijkerbuijk et al., *Cancer Res* 2010;70:4891-900 (PMID 20516114) reports absence of transcript from truncating mutants. Ochiai et al., *PNAS* 2014;111:1461-6 (PMID 24344301) identified an intergenic substitution 44 kb upstream of a *BUB1B* transcription start site and **reproduced the syndrome by editing that single base pair into cultured human cells**, which produced reduced transcript, increased PCS and MVA. |
| L02 | Reduced transcript plus accelerated turnover of the missense protein lowers BubR1 abundance | Suijkerbuijk 2010: "Low protein abundance is the direct result of the absence of transcripts from truncating mutants combined with high protein turnover of missense mutants." |
| L03 | Insufficient BubR1 impairs the mitotic checkpoint complex, and the effect is dosage-graded | Suijkerbuijk 2010 (MVA patient cell lines show impaired checkpoint, alignment defects and low BubR1). Bohers et al., *Hum Genet* 2008;124:473-8 (PMID 18932004) produced a residual-BubR1 gradient (8.5 / 10 / 14 / 58 / 77%) and found PCS across the range with **aneuploidy appearing below ~50% residual protein**. |
| L04 | The checkpoint complex restrains APC/C-Cdc20; failure permits premature anaphase | Lara-Gonzalez, Pines & Desai, *Semin Cell Dev Biol* 2021;117:86-98 (PMID 34210579); Kops, Weaver & Cleveland, *Nat Rev Cancer* 2005;5:773-85 (PMID 16195750). |
| L05 | BubR1 stabilises kinetochore–microtubule attachments by recruiting PP2A-B56 — **a route independent of the checkpoint** | Suijkerbuijk et al., *Dev Cell* 2012;23:745-55 (PMID 23079597); Kruse et al., *J Cell Sci* 2013;126:1086-92 (PMID 23345399); Lampson & Kapoor, *Nat Cell Biol* 2005;7:93-8 (PMID 15592459). |
| L06–L08 | Checkpoint failure → premature chromatid separation → whole-chromosome missegregation; destabilised attachments do the same by a second route | Bohers 2008; Suijkerbuijk 2010; Matsuura et al., *Am J Med Genet A* 2006;140:358-67 (PMID 16411201). |
| L09 | Aneuploidy imposes proteotoxic and metabolic stress | Oromendia, Dodgson & Amon, *Genes Dev* 2012;26:2696-708 (PMID 23222101); Donnelly et al., *EMBO J* 2014;33:2374-87 (PMID 25205676); Santaguida et al., *Genes Dev* 2015;29:2010-21 (PMID 26404941); reviewed in Santaguida & Amon, *Nat Rev Mol Cell Biol* 2015;16:473-85 (PMID 26204159). |
| L10 | Missegregation triggers p53/p16 arrest, senescence and immune clearance | Thompson & Compton, *J Cell Biol* 2010;188:369-81 (PMID 20123995); Santaguida et al., *Dev Cell* 2017;41:638-651 (PMID 28633018); Li et al., *PNAS* 2010;107:14188-93 (PMID 20663956); Baker et al., *Nat Cell Biol* 2008;10:825-36 (PMID 18516091). |
| L13 | Aneuploidy raises cancer risk, **non-monotonically** | Hanks et al., *Nat Genet* 2004;36:1159-61 (PMID 15475955); Jacquemont et al., *Am J Med Genet* 2002;109:17-21 (PMID 11932988) — 3 of 14 reported MVA cases (21%) developed a malignancy: rhabdomyosarcoma, acute lymphoblastic leukaemia and nephroblastoma. The non-monotonicity is from Weaver et al., *Cancer Cell* 2007;11:25-36 (PMID 17189716) and Silk et al., *PNAS* 2013;110:E4134-41 (PMID 24133140): low-rate missegregation promotes tumours, very high-rate missegregation suppresses them. |

### 2.2 Inferred links — the weak points, named

Two of thirteen links are **not** directly demonstrated, and the pipeline flags both.

- **L11 — aneuploidy-associated stress *causes* the arrest/senescence response.** Both
  phenomena are documented; that one causes the other, rather than both being parallel
  consequences of missegregation, is not established.
- **L12 — progenitor attrition produces growth restriction and microcephaly. This is
  the weakest joint in the chain, and it carries the entire therapeutic rationale for
  the growth phenotype.** The mouse work (Baker et al., *Nat Genet* 2004;36:744-9,
  PMID 15208629) establishes that BubR1 insufficiency produces dwarfism and progeroid
  phenotypes, but *by what route* is not shown in human *BUB1B*-mutant tissue — and
  there is at least one published alternative: BubR1 insufficiency also impairs
  **ciliogenesis** (Miyamoto et al., *Hum Mol Genet* 2011;20:2058-70, PMID 21389084),
  which is a developmental mechanism that has nothing to do with progenitor depletion.
  A cell-autonomous growth defect is a third possibility. Which tissues and which
  developmental windows dominate is unresolved.

The mechanism relevance score is **0.706** — a heuristic weighting of curated strength
labels, not a probability, and explicitly not a claim that the chain is *true*.

### 2.3 One node is compensatory, and it is not a target

`N10_clearance` — p53/p16-dependent arrest, senescence and immune clearance of
missegregating cells — deviates from wild type in the patient, and the naive move is to
push every deviation back. **That would suppress a protective response in a child whose
cancer risk is the thing most likely to kill.** The literature is genuinely split:

- In BubR1-hypomorphic mice, p16-driven senescence **drives** the progeroid phenotype;
  deleting p16 attenuated it (Baker 2008, PMID 18516091), and clearing p16-positive
  senescent cells delayed age-related deterioration (Baker et al., *Nature*
  2011;479:232-6, PMID 22048312 — performed in the BubR1 progeroid model).
- The same p53 arm **suppresses** aneuploidy-driven tumorigenesis (Li 2010,
  PMID 20663956), and immune clearance removes missegregating cells (Santaguida 2017,
  PMID 28633018).

Same response; pathological for growth and ageing, protective against cancer. This is
exactly the case where a signed state field cannot tell you which. The node is flagged
`deviation_is_pathological = 0`, the direction check returns `UNKNOWN` rather than
guessing, and the model **refuses to let it be designated the therapeutic target at
all** (ASSUMPTION-MECHANISM-04).

### 2.4 Do BubR1's non-checkpoint roles matter here? Yes, decisively

They are the reason the therapeutic target is the protein and not the checkpoint. MVA
patient cells show *both* an impaired checkpoint *and* chromosome alignment defects
(Suijkerbuijk 2010). Alignment depends on the KARD–PP2A-B56 arm (L05), which is
independent of MCC assembly. Any intervention aimed at the checkpoint addresses one
arm and leaves the other running. The protein-abundance node is the last point at which
one quantity governs both. Written up as ADR 0025 and generalised as
ASSUMPTION-MECHANISM-06.

A second non-checkpoint fact narrows the pharmacology sharply: **BubR1 is a
pseudokinase** (Suijkerbuijk et al., *Dev Cell* 2012;22:1321-9, PMID 22698286). There
is no catalytic activity to agonise. Whatever the therapy is, it cannot be an activator
of BubR1 enzyme function, because there isn't any.

---

## 3. Direction of effect: the trap, and the counter-trap

### 3.1 The naive move and the naive counter-move

The naive move is *"the checkpoint is weak, so strengthen it."* There is no approved
agent that strengthens the SAC; the class does not exist.

The naive counter-move is the oncology logic: *"inhibit MPS1 or Aurora B, weaken the
checkpoint further, and drive lethal aneuploidy."* MPS1 inhibitors do exactly what they
say — reversine abrogates the checkpoint and causes missegregation (Santaguida et al.,
*J Cell Biol* 2010;190:73-87, PMID 20624901) — and the class is in active anticancer
development (Wang et al., *Eur J Med Chem* 2019;175:247-268, PMID 31121430).

**This child is not a tumour.** The oncology logic works by pushing chromosomally
unstable cells past an aneuploidy-tolerance ceiling. This child's non-tumour cells are
already at that ceiling in every tissue, with no reserve. Target proximity is maximal
and the sign is inverted. In the pipeline, reversine scores **0.110, the lowest of all
13 agents**, and is rejected with `WRONG_DIRECTION` — the inversion is mechanical, not
editorial, and no weight change can resurrect it (a model validator forbids constructing
a wrong-direction hypothesis as accepted).

### 3.2 The tension the challenge asks about, confronted

The therapeutic window genuinely differs between tumour and soma, and this is not a
detail — it is the central structural fact about treating this child.

Tang et al., *Cell* 2011;144:499-512 (PMID 21315436) screened for **aneuploidy-selective
antiproliferation compounds** and found three: AICAR (energy stress), 17-AAG (HSP90 /
protein-folding inhibition), and chloroquine (autophagy inhibition). Each kills cells by
exploiting a vulnerability that arises *from being aneuploid*. Aneuploid cells depend on
autophagic degradation to manage their protein burden (Santaguida 2015) and have
impaired HSP90-dependent folding (Donnelly 2014).

For the rhabdomyosarcoma, that is an attractive property. For the child, **the same
property is a description of every cell in their body.** Hydroxychloroquine and
metformin are approved, cheap, extremely well tolerated, and have extensive paediatric
exposure — which is precisely what makes them dangerous to a repurposing search that
scores availability and tolerability but cannot see a sign. Both are rejected here with
`WRONG_DIRECTION` at the proteotoxic-stress node, **on the same evidence that
recommends them for the tumour.**

What this implies, stated plainly:

1. **A tumour-directed and a soma-directed hypothesis are two different analyses with
   opposite signs**, and merging them produces a recommendation that is wrong for
   whichever context the reader is actually in. This pipeline scores **one** context —
   the soma — and says so (ASSUMPTION-DRUG-08). Every rejection above means "rejected
   for systemic use in the affected child", never "useless in this disease".
2. The question that decides the tumour arm is not answerable from a mechanism chain.
   It is a **therapeutic index**: how much more sensitive is the tumour than the
   patient's own non-tumour cells, which share the aneuploidy that makes the tumour
   sensitive? Every tumour-plausible agent's `validation_experiment` cell names that
   comparison, because it is the whole question.
3. **The primary hypothesis was chosen partly because it escapes the tension.** Raising
   BubR1 abundance is the rare intervention whose direction is desirable in *both*
   contexts: Baker 2013 showed sustained elevated BubR1 reduced aneuploidy **and
   reduced tumorigenesis**, even against oncogenic Ras. That is the opposite of the
   checkpoint-inhibitor case, and it is the strongest argument for preferring an
   upstream, abundance-directed hypothesis over any consequence-directed one.

### 3.3 The classes considered, with directions

| Class | Representative | Node | Direction | Desired in **soma**? | Desired in **tumour**? | Outcome |
|---|---|---|---|---|---|---|
| **Restore BubR1 abundance** (NAD⁺/SIRT2 axis) | NMN, nicotinamide riboside, nicotinamide, nicotinic acid | BubR1 abundance | ↑ (NMN measured); unmeasured or bidirectional for the rest | **Yes** | **Yes** | **Primary hypothesis.** Best agent is unapproved; approved agents have undetermined sign. |
| **Readthrough / restore full-length protein** | Ataluren, gentamicin | BubR1 abundance | ↑ *if* a nonsense allele produces a transcript | Yes, conditionally | Yes | Rejected. Genotype-conditional **and** mechanistically undercut: truncating *BUB1B* alleles yield no transcript (Suijkerbuijk 2010), so NMD removes the substrate before the ribosome reaches the stop codon. Ataluren's EU authorisation was also not renewed (EMA; EC decision effective 2025). |
| **Reduce aneuploidy consequences — proteostasis** | Sirolimus | Proteotoxic stress | Unknown in an aneuploid background | Unclear | Unclear | Rejected. Also anti-growth in a growth-restricted child and immunosuppressive in a child whose aneuploid cells are partly held in check by immune clearance — two objections at *other* nodes, which single-node scoring cannot express. |
| **Aneuploidy-selective antiproliferatives** | Hydroxychloroquine, tanespimycin (17-AAG), metformin | Proteotoxic / metabolic stress | ↑ stress (kills aneuploid cells) | **No — catastrophic** | **Yes** | Rejected for the soma with `WRONG_DIRECTION`. Retained in the record as tumour-arm candidates requiring a therapeutic-index study. |
| **Cell-cycle / checkpoint modulation** | Reversine (MPS1i) | SAC | ↓ checkpoint | **No — pushes the disease mechanism** | Yes | Rejected, lowest score of all. The naive top hit. |
| **Mitotic timing — restrain APC/C directly** | proTAME, apcin | SAC / APC/C | ↑ restraint on APC/C-Cdc20 | Yes, in principle | — | **Not modelled, and this is the one class whose direction is right and whose absence is purely regulatory.** Partial APC/C inhibition substitutes for the restraint the failing checkpoint cannot supply. Both are tool compounds never given to humans, so under ASSUMPTION-DRUG-03 they are research probes, not repurposing candidates; they would add a third row saying only `NOT_APPROVED`, which reversine and tanespimycin already demonstrate. Recorded in `SOURCES.md` under "considered and left out". The obvious objection is a therapeutic window: too much APC/C restraint is mitotic arrest and slippage, which generates aneuploidy of its own. |
| **Acetylation / degradation modulation** | Vorinostat (HDACi class) | SAC | ↓ checkpoint | **No** | — | Rejected. Argued for by analogy — BubR1 is stabilised by K250 acetylation (Choi et al., *EMBO J* 2009;28:2077-89, PMID 19407811) and HDAC inhibitors raise acetylation — but measured at the mitotic readout, the class *induces premature sister-chromatid separation and overrides the SAC* (Magnaghi-Jaulin et al., *Cancer Res* 2007;67:6360-7, PMID 17616695). A mechanism-shaped argument that reverses when someone measures the actual readout. |
| **Senolytics** | Dasatinib (with quercetin; Zhu et al., *Aging Cell* 2015;14:644-58, PMID 25754370) | Clearance / senescence (**compensatory**) | ↓ senescent-cell burden | **Cannot determine** | Possibly harmful | Rejected. Motivated by the best animal evidence in this whole document (Baker 2011 in the exact BubR1 model) and blocked by §2.3: the response it removes may be a barrier rather than a burden. |
| **Antioxidants** | N-acetylcysteine | — | ↓ ROS | Plausible | **No** | Not modelled. The oxidative-stress node is not separately represented, and the strongest relevant finding points the wrong way for a cancer-predisposed child: antioxidants accelerated tumour progression in mice (Sayin et al., *Sci Transl Med* 2014;6:221ra15, PMID 24477002). |
| **Surveillance** (non-pharmacological) | Renal ultrasound, clinical review | — | Does not move any node | Yes — highest near-term value | — | Deliberately kept out of the drug catalogue (§7.1). |

---

## 4. The primary hypothesis, stated so it can be attacked

> **In BUB1B-MVA alleles of the *abundance-lowering* class, pharmacological elevation
> of the NAD⁺ pool raises BubR1 protein abundance through SIRT2-dependent deacetylation
> of BubR1 lysine-668, and the resulting increase — if it crosses the ~50% residual
> threshold — reduces the ongoing rate of chromosome missegregation.**

### 4.1 Why this direction is the right one

Four independent results, all pointing the same way:

1. **Ectopic BubR1 expression restored mitotic checkpoint activity in MVA patient
   cells.** Suijkerbuijk 2010 states it as proof that BubR1 dysfunction *causes* the
   segregation errors in these patients. This is a rescue experiment in the actual
   human disease, which is more than most rare-disease repurposing hypotheses have.
2. **Raising BubR1 in vivo is beneficial and not oncogenic.** Baker et al., *Nat Cell
   Biol* 2013;15:96-102 (PMID 23242215): sustained high-level BubR1 preserved genomic
   integrity, reduced
   tumorigenesis even against oncogenic Ras, corrected both checkpoint impairment
   *and* microtubule–kinetochore attachment defects, and extended lifespan.
3. **Lowering BubR1 in vivo reproduces the phenotype.** Baker 2004 (PMID 15208629).
4. **The relationship is graded with a threshold to climb back over.** Bohers 2008
   (PMID 18932004): aneuploidy emerges below ~50% residual BubR1. Partial restoration
   is therefore a coherent goal rather than an all-or-nothing one.

### 4.2 Why the NAD⁺/SIRT2 axis specifically

North et al., *EMBO J* 2014;33:1438-53 (PMID 24825348) showed that BubR1 abundance is
set by the acetylation state of **lysine-668**: SIRT2 holds K668 deacetylated in an
NAD⁺-dependent manner, CBP acetylates it, and the balance determines BubR1 turnover.
Treating mice with the NAD⁺ precursor **NMN increased BubR1 abundance in vivo**, and
SIRT2 overexpression in **BubR1^H/H** hypomorphic animals increased median lifespan.

That is a druggable, upstream handle on the exact quantity the disease is short of.

### 4.3 The three objections a sceptical expert raises first — and my answers

**"The NAD⁺ evidence comes from an ageing context. A young child has not lost NAD⁺, so
where is the headroom?"** This is the strongest objection and I do not have a good
answer. North 2014 frames the axis as an explanation for the *age-related decline* in
BubR1; the intervention restores something that was lost. In a young child, NAD⁺ is
presumably not the limiting quantity, and the pharmacology may simply have nothing to
push against. This is the first thing the experiment in §6 measures, and it is a
plausible route to a clean negative.

**"Suijkerbuijk 2010 describes two mutation classes. Half of them are qualitative
defects that more protein will not fix."** Correct, and the hypothesis is therefore
explicitly **genotype-conditional** (ASSUMPTION-DRUG-09). See §7.3.

**"Is more BubR1 automatically safe?"** No, and I will not claim it is. Mad2
overexpression is tumorigenic in mice (Sotillo et al., *Cancer Cell* 2007;11:9-23,
PMID 17189715), so "more checkpoint protein" is not monotonically good as a class
statement. BubR1 overexpression specifically was protective (Baker 2013), but any
dose-finding has to look for a ceiling as well as a floor.

### 4.4 The direction verdict, at the molecule

| Agent | Approval | Direction on BubR1 abundance | Verdict | Fate |
|---|---|---|---|---|
| **NMN** | Investigational | **increase** — measured in vivo (North 2014) | **AGREES** | Rejected: not an approved medicine (ASSUMPTION-DRUG-03), and its effect on aneuploidy/cancer risk has never been assessed (blocking, ADR 0011). Scores **0.605 — the highest in the catalogue — and is still rejected.** |
| **Nicotinamide riboside** | Food ingredient | **unknown** — raises human NAD⁺ (Trammell, *Nat Commun* 2016;7:12948, PMID 27721479; Martens, *Nat Commun* 2018;9:1286, PMID 29599478) but BubR1 never measured | cannot determine | Rejected: not an approved medicine. |
| **Nicotinamide** | Approved (pellagra) | **context_dependent** — NAD⁺ precursor *and* a product inhibitor of sirtuins (Bitterman et al., *J Biol Chem* 2002;277:45099-107, PMID 12297502), so it feeds and inhibits the same axis | **cannot determine** | **Accepted, rank 1, score 0.365**, with the flag attached. |
| **Nicotinic acid** | Approved (dyslipidaemia) | **unknown** — raises human NAD⁺ several-fold with functional benefit in mitochondrial myopathy (Pirinen et al., *Cell Metab* 2020;31:1078-1090, PMID 32386566); no sirtuin-inhibiting liability; BubR1 never measured | cannot determine | Rejected on the unassessed oncogenic-risk gate — which is unusually cheap to close (§5.2). |

**The nicotinamide sign paradox is the non-obvious finding here.** The cheapest,
safest, best-paediatric-evidenced NAD⁺-related agent — ENDIT randomised 552
first-degree relatives, children included, to 1.2 g/m²/day or placebo for five years,
and reported no effect on growth in children (Gale et al., *Lancet* 2004;363:925-31,
PMID 15043959) — is the one whose direction cannot be signed, because it inhibits the
very sirtuin the mechanism depends on. A repurposing search that reads "NAD⁺ precursor
→ good" gets this exactly backwards, and might get the *sign* backwards.

---

## 5. What the pipeline actually produced

Reproduce with `uv run pytest -q tests/unit/test_track2_bub1b.py` (22 tests). Tables:
`knowledge/literature/bub1b/`. Output: **1 accepted, 12 rejected, 89 evidence rows**,
every rejection kept with its reasons (GP-19).

| Rank | Agent | Node | Required | Observed | Direction verdict | Score | Outcome |
|---|---|---|---|---|---|---|---|
| **1** | Nicotinamide | `N3_bubr1` | increase | context_dependent | **cannot determine** | 0.365 | accepted, 3 concerns |
| — | Nicotinamide mononucleotide | `N3_bubr1` | increase | increase | **AGREES** | 0.605 | rejected — `oncogenic_risk` |
| — | Gentamicin | `N3_bubr1` | increase | unknown | cannot determine | 0.420 | rejected — `oncogenic_risk` |
| — | Nicotinic acid | `N3_bubr1` | increase | unknown | cannot determine | 0.410 | rejected — `oncogenic_risk` |
| — | Dasatinib | `N10_clearance` | *unknown (compensatory)* | decrease | cannot determine | 0.348 | rejected — `oncogenic_risk` |
| — | Nicotinamide riboside | `N3_bubr1` | increase | unknown | cannot determine | 0.345 | rejected — `oncogenic_risk` |
| — | Sirolimus | `N9_stress` | decrease | unknown | cannot determine | 0.303 | rejected — `oncogenic_risk` |
| — | Ataluren | `N3_bubr1` | increase | unknown | cannot determine | 0.203 | rejected — `oncogenic_risk` |
| — | Hydroxychloroquine | `N9_stress` | decrease | increase | **OPPOSES** | 0.163 | rejected — `wrong_direction` |
| — | Metformin | `N9_stress` | decrease | increase | **OPPOSES** | 0.163 | rejected — `wrong_direction` |
| — | Vorinostat | `N5_sac` | restore | decrease | **OPPOSES** | 0.133 | rejected — `wrong_direction` |
| — | Tanespimycin | `N9_stress` | decrease | increase | **OPPOSES** | 0.130 | rejected — `wrong_direction` |
| — | Reversine (MPS1i) | `N5_sac` | restore | decrease | **OPPOSES** | 0.110 | rejected — `wrong_direction` |

### 5.1 Three things this table says that a prose write-up would hide

**The highest-scoring agent is rejected.** NMN outscores everything and does not
survive, because approval status and the unassessed oncogenic-risk question are *gates*
rather than score terms. A pipeline where the top score always wins cannot express
"right about the biology, not yet a candidate".

**The accepted candidate scores 0.365 out of 1 and carries an undetermined direction.**
The ranking is telling the reader there is no good answer yet. Presenting rank 1 as a
recommendation would be over-reading the pipeline's own output, and this document does
not do that.

**Nicotinamide ranks first largely because its cancer-susceptibility question has an
answer at all** — a phase 3 randomised trial with a cancer-incidence endpoint that found
*fewer* cancers (Chen et al., *N Engl J Med* 2015;373:1618-26, PMID 26488693). That trial was in
adults with UV-driven skin cancer and does not transfer to a germline
chromosomal-instability population (ASSUMPTION-DRUG-06). It is enough to clear a gate,
not enough to support a claim.

### 5.2 The binding constraint is not target identification

Seven of twelve rejections are `oncogenic_risk`, and in every one of those seven cases
the risk was not *found* — it was **never measured**. In a chromosomal-instability
syndrome, "nobody has asked whether this agent increases aneuploidy" is disqualifying
(ADR 0011), and it disqualifies almost everything.

This is the most transferable finding in the document. **The bottleneck in repurposing
for chromosomal-instability disorders is not finding targets. It is that the one safety
question that matters most has not been asked of the candidate drugs.** For nicotinic
acid it is answerable *from existing data* — decades of large cardiovascular outcome
trials with cancer-incidence data already collected — without running a new study. That
is a concrete, cheap, high-value next step that fell out of the analysis rather than
being assumed into it.

---

## 6. The experiments: one to falsify, one to confirm

A hypothesis with no stated disconfirming experiment is not a scientific claim.

### 6.1 The falsifying experiment

**Design.** Primary fibroblasts or LCLs from BUB1B-MVA patients (biobanked lines; this
does not require the index child), plus isogenic *BUB1B*-hypomorphic hTERT-RPE1 lines
engineered to represent each of the two published mutation classes, and an
allele-corrected control. Treat across a clinically achievable NAD⁺-precursor exposure
range. Blinded scoring throughout.

**Readouts, in causal order:**

1. Intracellular NAD⁺ (LC-MS/MS)
2. BubR1 K668 acetylation (site-specific immunoprecipitation / targeted MS)
3. Total BubR1 protein (quantitative immunoblot + targeted MS), expressed as a
   percentage of the isogenic corrected control — the Bohers 2008 scale
4. Premature chromatid separation frequency and micronucleus rate, blinded
5. Checkpoint strength (mitotic index under nocodazole challenge)
6. Aneuploidy burden (single-cell karyotyping or low-pass single-cell WGS)

**The hypothesis is falsified if any of these holds:**

- NAD⁺ rises and K668 acetylation falls, **but total BubR1 does not increase.** The
  axis exists but has no headroom in these cells — the objection in §4.3, made concrete.
- **BubR1 increases but missegregation does not fall.** Right biomarker, no phenotype.
  This is the failure mode a surrogate-endpoint-only study would miss, which is why
  readouts 4–6 are mandatory rather than exploratory.
- Rescue **persists under SIRT2 knockdown.** Whatever is happening is not the proposed
  mechanism, and the hypothesis as stated is wrong even if the effect is real.
- The effect is present only in the abundance-lowering allele class **and** the
  qualitative class shows no benefit — which does not falsify the mechanism but does
  bound the eligible population, and must be reported as such rather than averaged away.

**Controls that make it a test rather than a demonstration:** ectopic BubR1
re-expression as the ceiling control (how much of the deficit is recoverable at all?);
SIRT2 knockdown as the mechanism-specificity control; the allele-corrected isogenic line
as the floor.

### 6.2 The confirming experiment

A monotone dose–response in which NAD⁺ precursor exposure raises BubR1 across the
Bohers threshold and produces a **proportional, blinded** reduction in PCS and
micronucleus frequency, abolished by SIRT2 knockdown and not exceeded by ectopic BubR1
re-expression. Then, and only then, in vivo in **Bub1b^H/H** mice at a clinically
scalable exposure, with **aneuploidy burden and tumour incidence as co-primary
endpoints.** Tumour incidence is co-primary, not secondary: in a cancer-predisposition
syndrome, a healthspan-only readout cannot answer the question that matters, and Baker
2011's senolytic result is the cautionary example of an intervention that looked good on
healthspan while the cancer question sat unasked.

### 6.3 The experiment that must come first, before either

**Determine the allele class** (ASSUMPTION-DRUG-09). Quantify *BUB1B* transcript and
BubR1 protein in the patient's own cells against a control, and — if abundance is
reduced — perform the Suijkerbuijk 2010 functional replacement assay to establish
whether the residual protein is quantitatively or qualitatively defective. That single
assay decides whether this entire hypothesis is applicable, and it is run by the team
holding the sample, not by anyone reading this.

---

## 7. Honest limitations

### 7.1 What the model could not express

- **Non-pharmacological interventions have no home in a drug catalogue.** Tumour
  surveillance — the highest expected-value intervention available today, given a 21%
  malignancy rate in early case series (Jacquemont 2002) and established protocols for
  Wilms tumour surveillance in at-risk children (Scott et al., *Arch Dis Child*
  2006;91:995-9, PMID 16857697; updated in Kalish et al., *Clin Cancer Res*
  2024;30:5260-5269, PMID 39320341) — does not move any node on this chain. It changes
  *stage at diagnosis*. Forcing it into the catalogue would have produced either a
  category error or a spurious `WRONG_DIRECTION` rejection, so it is discussed here and
  excluded there. **If one thing in this document has near-term value for a real child,
  it is this paragraph, and it is not a drug.**
- **Agent-specific toxicity has no field.** Gentamicin's cumulative ototoxicity and
  nephrotoxicity make chronic dosing implausible for a lifelong condition; sirolimus is
  anti-growth in a growth-restricted child and immunosuppressive in a child whose
  aneuploid cells are partly held in check by immunity; and the paediatric dasatinib
  trial reported bone growth and development events in 4% of patients (Gore et al.,
  *J Clin Oncol* 2018;36:1330-1338, PMID 29498925) — which is a specific concern in
  exactly this phenotype. None of these is expressible in the schema, so they live in
  the `notes` cells and in this section. The safety pass says so itself: "Absence of a
  concern here reflects the fields the catalogue records, not a safety clearance."
- **Single-node scoring cannot see cross-node harm.** An agent is scored against the one
  node it targets. Sirolimus's objections all live at *other* nodes.

### 7.2 The weakest link, named

**L12.** That progenitor attrition is what produces growth restriction and microcephaly
is inferred, not demonstrated in human *BUB1B*-mutant tissue, and there is a published
alternative route through ciliogenesis (Miyamoto 2011). Everything this document says
about the *growth* phenotype rests on it.

This matters less than it looks, for a reason worth stating: **the developmental-window
caveat already removes growth from the achievable endpoints** (ASSUMPTION-MECHANISM-02).
Growth restriction and microcephaly in MVA are established in utero and are structural
at birth (Hanks 2004; Bohers 2008; García-Castillo et al., *Am J Med Genet A*
2008;146A:1687-95, PMID 18548531). No post-natal increase in BubR1 can reverse them.
What post-natal correction could plausibly change is the *rate at which new aneuploid
cells are generated from that point on* — and therefore the trajectory of cancer risk
and progeroid tissue attrition. **A trial endpoint has to match that claim: a
cytogenetic or tumour-incidence endpoint is arguable; a head-circumference endpoint is
not.** L12 being weak weakens a phenotype the therapy was never going to reach.

The link that would hurt most if it broke is **L03** — that BubR1 abundance is the
governing quantity — and it is the best-evidenced link in the chain.

### 7.3 The n-of-1 problem, stated plainly

**This is one child, and a mechanism-grounded hypothesis is not evidence of efficacy in
that child.** Everything above is reasoning about a *gene* from a literature of case
series, cell lines and mouse models. Specifically:

- The largest human dataset behind any of this is a handful of MVA patient cell lines.
  There has never been a clinical trial in MVA, and there will not be one powered
  conventionally: the condition is too rare.
- A mechanism that is correct in aggregate can be irrelevant in an individual —
  particularly here, where the published allele classes predict *opposite* answers to
  "would more protein help?" and the class is unknown to me by design.
- The n-of-1 designs that could generate real evidence here (N-of-1 crossover on a
  cytogenetic biomarker, or an aggregated-N-of-1 series across the international MVA
  cohort) require a validated, responsive biomarker — which §6 is, in part, an argument
  for building.
- Absence of a signal in one child is not absence of an effect, and presence of one is
  not proof of it. This is the reason the deliverable is an experiment and not a
  recommendation.

### 7.4 Grading of this work

Under GP-20: the mechanism chain and drug catalogue are **literature-derived** — a
third maturity category, defined in ADR 0022, distinct from the repository's synthetic
demo tables and from its machine-acquired public reference data. They are real biology,
curated by hand, with the integrity guarantee "every row resolves to a cited
publication" (enforced by test) rather than "these bytes match a release hash".

---

## 8. Scalability: what generalises and what does not

### 8.1 Reusable machinery (nothing below is BUB1B-specific)

| Component | Why it transfers |
|---|---|
| `EffectDirection` + `directions_agree` (tri-state) | A signed vocabulary with an explicit "cannot determine". Applies to any mechanism where a therapy has a sign. |
| The **direction-triple validator** on `MechanismHypothesis` | `disease_direction`, `required_correction` and the target node's `state_in_patient` are authored independently and must agree at construction. One mistyped cell would otherwise invert the entire gate silently. Gene-agnostic. |
| `deviation_is_pathological` (mandatory, un-defaulted) | The compensatory-node concept. Every disease with a damage response has this problem: DNA-repair disorders (p53 attrition), lysosomal disease (autophagy upregulation), mitochondrial disease (mitophagy). |
| Rejection reasons as a **priority-ordered** vocabulary | "Wrong direction" and "not approved" are different findings with different follow-ups ("never" vs "wait"). The ordering makes reports say the right one. |
| Context detectors (`is_chromosomal_instability_context`) | Which gaps are *fatal* is disease-dependent. The mechanism itself selects the gate. |
| Evidence-tier ladder + `is_direct_evidence` discount | Keeps a binding assay from being read as a result. Universal. |
| The two-TSV loader contract | A new disease is **two tables**, not a code change. |

### 8.2 What is BUB1B-specific

Exactly three files: 12 nodes and 13 links, 13 catalogue rows, one metadata row —
about 40 lines of curation and a `SOURCES.md`. **No Python was written for BUB1B.** The
Track 2 test module loads the literature tables through the same public API the pipeline
uses (`MechanismLibrary.from_tsv`, `DrugCatalog.from_tsv`) and runs the same
`build_mechanism` / `generate_drug_hypotheses` functions the synthetic demo runs. That is
the scalability claim, and it is checkable by inspection: the files added for BUB1B are
three TSVs, a `SOURCES.md`, a deterministic generator under `tools/`, a test module and
two decision records. **Not one line under `src/` was added or changed for this gene** —
the direction gate, the compensatory-node rule, the tri-state agreement and the safety
context detectors were all already there, written for a fictional gene, and they worked
unmodified on the real one.

### 8.3 The adapter boundary, and the honest bottleneck

The catalogue is a `DrugCatalog` behind a TSV loader; swapping it for ChEMBL, DrugBank
or Open Targets is a boundary change, not a logic change (GP-02: parsing happens at the
boundary, once, into typed models).

But scaling has a real obstacle, and gesturing past it would be dishonest:
**`observed_direction` is not a column in any public drug database.** ChEMBL records
targets and potencies; Open Targets records associations. None records *which way the
compound pushes the node*, and the whole Track 2 gate depends on that sign. Scaling this
approach means curating or inferring signs at volume, which is where the work is. Three
tractable routes: harvest MoA verbs from ChEMBL/DrugBank text (inhibitor/agonist maps
cleanly onto ↓/↑ for a *direct* target and not at all for downstream nodes); use
perturbation signatures (LINCS L1000) to read the sign empirically at a transcriptional
node; or restrict automated ingestion to direct-target rows and require curation for
everything else. The tri-state means the last option degrades gracefully — unknown
signs become flagged candidates rather than fabricated ones.

### 8.4 Same disease label, different mechanisms

MVA is caused by *BUB1B* (MVA1), *CEP57* (MVA2; Snape et al., *Nat Genet* 2011;43:527-9,
PMID 21552266) and *TRIP13* (Yost et al., *Nat Genet* 2017;49:1148-1151, PMID 28553959).
Same clinical label; a checkpoint-signalling mechanism, a centrosomal mechanism and a
checkpoint-silencing mechanism. **Everything in §4 is irrelevant to two of the three** —
raising BubR1 abundance does nothing for a CEP57 patient. Because mechanism is resolved
per gene and never per disease name (ASSUMPTION-MECHANISM-03), the same catalogue run
against a *CEP57* chain would produce different verdicts for the same compounds, and
that is a demonstrable property rather than a claim.

The broader class this generalises to: chromosomal-instability and mitotic-checkpoint
disorders (MVA, Roberts syndrome, Warsaw breakage, PCS syndrome), DNA-repair
instability syndromes (Fanconi anaemia, Bloom, Nijmegen, ataxia-telangiectasia) and, in
general, **any monogenic disorder that is a dosage problem rather than a signalling
problem** — where the same question applies: is there a threshold, is it graded, and is
there an approved agent that moves the quantity the right way? For every one of them, the
oncology literature will offer compounds that bind the right protein and push the wrong
way, and the gate in §3.1 is what stops them.

---

## 9. Reproducing this

```bash
uv run pytest -q tests/unit/test_track2_bub1b.py   # 22 tests: the chain, the gate, the citations
just verify                                        # the full acceptance gate
```

| Artefact | Path |
|---|---|
| Mechanism chain (12 nodes, 13 links, individually graded) | `knowledge/literature/bub1b/mechanisms.tsv` |
| Target node and required correction | `knowledge/literature/bub1b/mechanism_meta.tsv` |
| Drug catalogue (13 agents, signed directions) | `knowledge/literature/bub1b/drug_catalog.tsv` |
| Every citation, mapped row by row | `knowledge/literature/bub1b/SOURCES.md` |
| Catalogue generator | `tools/build_knowledge/bub1b_catalog.py` |
| Executable assertions | `tests/unit/test_track2_bub1b.py` |
| Why literature tables sit outside `knowledge/public/` | `docs/decisions/0022-literature-tables-are-a-third-knowledge-category.md` |
| Why the target is abundance, not the checkpoint | `docs/decisions/0025-the-target-node-is-bubr1-abundance-not-the-checkpoint.md` |
| Why an unassessed CIN risk blocks | `docs/decisions/0011-unassessed-cin-risk-is-blocking.md` |
| The assumptions this rests on | `docs/scientific-assumptions.md` |

---

*Nothing in this document is medical advice. Every output is a research hypothesis
requiring pre-clinical validation. No patient data was accessed in producing it.*
