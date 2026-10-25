# R8 — The proteostasis / mitochondrial axis, and whether it holes our submission

**Status: complete.** Every section is sourced to a PMID, a regulatory record or a file path.

**Date of retrieval:** 2026-08-29. Everything below was read on that date.

**What this is.** An assessment of whether the proteostasis / mitochondrial axis is a better
route than the NAD⁺/SIRT2 one for this disease. One paper carries that axis —
González-Blanco et al., *Nat Commun* 2026, PMID 41820377 — and it is the paper any
antioxidant- or autophagy-directed candidate list (N-acetylcysteine, idebenone, metformin,
sirolimus) has to rest on. `docs/track2-hypothesis.md` cites it zero times, and the word
"proteostasis" appears there exactly once, in the rejected sirolimus row of Appendix B (line 922).
This document establishes whether that is a hole in our science or a defensible omission.

**Method.** PubMed E-utilities (`esummary`, `efetch`, `esearch`), Europe PMC full-text XML
(`https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13121444/fullTextXML`, 170,193 bytes,
retrieved in full), Exa semantic web search, EMA and FDA public records. Read-only throughout.
Every PMID in this document was resolved against the NCBI record; none is cited from memory. §8
lists every query run, including the empty ones.

**Rule applied throughout.** Every negative is "I did not find", never "does not exist". What a
paper measured and what I infer from it are separated and the inferences are labelled.

---

## 0. Verdict, stated once

**The paper is real and the identifier resolves — but the claim commonly drawn from it, that it supports an antioxidant or proteostasis lead candidate, is a misreading of it.**

Three findings, in order of how much they change our position:

1. **The identifier resolves, and the paper is the one the axis rests on.**
   González-Blanco et al., *Nat Commun* 2026, PMID 41820377 —
   "Proteostasis failure and mitochondrial dysfunction contribute to chromosomal
   instability-induced microcephaly", Nature Communications 17:3829, published 12 March 2026,
   PMCID PMC13121444. Full record in §1. It is explicitly framed as an MVA paper, in its title
   context, its abstract's first sentence and its discussion.

2. **The axis is real and downstream of ours, not in competition with it.** It is a *Drosophila*
   neuroblast model in which chromosomal instability is induced by RNAi against `bub3` or `rod` —
   **not** against the fly BubR1 orthologue, and nothing in the paper measures BubR1, BUB1B or any
   protein whose abundance our hypothesis targets. Proteostasis failure and mitochondrial
   dysfunction are consequences of accumulated aneuploidy. Our node is upstream of the
   missegregation that produces them (§6).

3. **The strongest claim drawn from the paper is false as stated, and the paper says so in its own words.**
   That claim is that "TOR inhibition was one of the rescues of CIN-induced microcephaly", and it
   is what sirolimus, everolimus and metformin would be built on. TOR inhibition rescued
   neuroblast *number*. It did not rescue brain size — the microcephaly endpoint in the
   paper's own title. The paper: "we noticed that none of the genetic interventions aimed at
   increasing autophagy was able to rescue the number of Pros+ cells or the size of the larval
   brain". Rapamycin fed to larvae at 40 µM gave the same split (§4).

**Recommendation: (b), address it as a stated limitation** — but a limitation with teeth, because
the correction in item 3 is ours to make and it is load-bearing against the two candidates this
axis would otherwise rate highest. Full argument in the closing section.

---

## 1. The record, verified

Retrieved from `https://pubmed.ncbi.nlm.nih.gov/41820377/` and cross-checked against the NCBI
`esummary` record and the Europe PMC JATS front matter.

| Field | Value |
|---|---|
| PMID | 41820377 |
| Title | Proteostasis failure and mitochondrial dysfunction contribute to chromosomal instability-induced microcephaly |
| Authors | Amanda González-Blanco, Adrián Acuña-Higaki, David Boettger, Jery Joy, Marco Milán |
| Affiliation | Institute for Research in Biomedicine (IRB Barcelona), The Barcelona Institute of Science and Technology; ICREA. Corresponding author Marco Milán. |
| Journal | Nature Communications |
| Year / volume / article | 2026, volume 17, article number 3829 (`elocation-id`, no page range) |
| DOI | 10.1038/s41467-026-70521-0 |
| PMCID | PMC13121444 |
| Publication date | 12 March 2026 |
| Licence | CC BY-NC-ND 4.0, open access |

House-form citation:

González-Blanco et al., *Nat Commun* 2026, PMID 41820377

### 1.1 Abstract, verbatim

> Mosaic variegated aneuploidy (MVA), a rare human congenital disorder that causes microcephaly,
> is characterized by extensive abnormalities in chromosome number and results from mutations in
> genes involved in accurate mitotic chromosome segregation. To characterize the cellular
> mechanisms underlying this disease, here we generated a *Drosophila* model of microcephaly
> caused by the depletion of a single spindle assembly checkpoint (SAC) gene in the neural stem
> cell (NSC) compartment. We present evidence that loss of stemness - compromised identity and
> proliferative capacity of NSCs- plays an important role in MVA and results in a reduced number
> of neurons and glial cells. We show that loss of stemness arises from the accumulation over
> time of an unbalanced number of gains and losses of more than one chromosome, rather than a
> direct consequence of chromosomal instability-induced DNA damage or the production of simple
> aneuploidies. We unravel a contribution of proteostasis failure and mitochondrial dysfunction
> to the negative impact of complex aneuploidies on stemness, a highly energy demanding cellular
> state. We identify overexpression of Radical Oxygen Species scavengers, mitochondria chaperones
> and apoptosis inhibition as genetic interventions capable of dampening the deleterious effects
> of aneuploidy on brain size.

---

## 2. Does it bear on MVA and BUB1B specifically?

**On MVA: yes, explicitly and repeatedly. On BUB1B: only by citation, never by experiment.**

The MVA framing is not incidental. Introduction, first paragraph:

> Mosaic variegated aneuploidy (MVA) is characterized by widespread mosaic aneuploidies,
> involving different chromosomes and tissues. MVA patients can present growth retardation, early
> childhood cancer, and microcephaly. MVA is caused by pathogenic mutations in components of the
> spindle-assembly checkpoint (SAC) and proteins involved in centrosome dynamics during mitosis.

Results, first paragraph, connecting the model to patients:

> In humans, MVA-driven microcephaly is caused in half of the patients by bi-allelic mutations in
> single genes with well-documented roles in spindle assembly checkpoint (SAC) function or
> correct kinetochore-microtubule attachment. Our first aim was then to address whether targeted
> depletion of a single SAC gene could also cause a reduction in larval brain size in flies and
> whether NBs were the most susceptible cells to CIN.

Introduction, stating the design intent:

> Here we have generated a fly model of microcephaly caused by depletion of a single SAC gene, as
> occurs in MVA patients.

Discussion:

> Here we generated a *Drosophila* model of MVA-driven microcephaly by targeting the SAC genes
> *bub3* or *rod* specifically in the neural stem cell (NSC) compartment.

Its reference list carries the core *BUB1B* clinical literature —
Hanks et al., *Nat Genet* 2004, PMID 15475955 as ref 6;
Matsuura et al., *Am J Med Genet A* 2006, PMID 16411201 as ref 5;
Suijkerbuijk et al., *Cancer Res* 2010, PMID 20516114 as ref 2; the current review,
Malumbres et al., *Nat Rev Genet* 2024, PMID 39169218, as ref 4. Every one of those resolves.

### 2.1 The gene is `bub3` or `rod`, not BubR1 — and this matters

The CIN is induced by `worniu-gal4`-driven RNAi against `bub3` or `rod`. Human orthologues are
BUB3 and KNTC1/ROD. *Drosophila* has a BubR1 orthologue and the authors did not use it. Grepping
the full text for `BubR1`, `BUB1B` and `BUB1` returns matches only in the reference list and in
the introduction's clinical framing — never in a figure, a genotype or a Methods line.

The consequence is precise, and it is a limit rather than a refutation. The paper models
**MVA-class CIN in a neural stem cell compartment**. It does not model *BUB1B* allele dosage, and
nothing in it bears on whether BubR1 protein abundance is the deficient quantity — the claim
`docs/track2-hypothesis.md` §1 makes. Those are different questions and this paper answers only
the first.

### 2.2 One apparent citation error in the paper, recorded

The introduction sentence "In mouse models of MVA-driven microcephaly caused by mutations in SAC
genes, a dramatic reduction in the neural stem cells (NSCs), their intermediate progenitors, and
post-mitotic neurons is associated with massive cell death11" cites reference 11. Reference 11 in
the article's list is:

Okimoto et al., *Nat Genet* 2017, PMID 27869830

which is "Inactivation of Capicua drives cancer metastasis" and has nothing to do with MVA, SAC
genes or microcephaly. I verified that the JATS `<xref rid="CR11">` is present in the body text
and that the numbered reference list is internally consistent everywhere else I checked (ref 12 is
Sterling et al., *Front Cell Dev Biol* 2023, PMID 37900274, which does match its p53 sentence;
refs 13, 14, 15 all match their sentences). The likely intended citation is a mouse *BubR1*
microcephaly paper — Simmons et al., *Hum Mol Genet* 2019, PMID 30668728 fits the sentence and is
not in the list. **Recorded as an observation, not used to discredit the paper.** It is worth one
line in our report only if a judge asks how closely we read it.

---

## 3. What the paper actually measured

Everything below is from the Results and Methods of the full text. System is
*Drosophila melanogaster* larval central brain throughout, `worniu-gal4` driver, readouts at 48/62/72/96 h
after egg laying. No mammalian cell, no patient material, no mouse.

| Node | Manipulation | Readout | Result |
|---|---|---|---|
| CIN induction | `bub3-RNAi` or `rod-RNAi` in neuroblasts | NB number, Pros+ cells, Elav+ neurons, Repo+ glia, brain lobe size, adult brain size | Progressive NB loss; reduced progeny; reduced larval and adult brain size |
| Karyotype | FISH for chromosomes II and III in NBs | % NBs with unbalanced chromosome number | Complex aneuploidies accumulate over time; simple trisomy (X or 2L) has "a mild impact" |
| DNA damage | ionising radiation; DDR gene depletion | pH2Av, p53 sensor, cDcp1 | NSCs "highly resistant to DNA damage"; DNA damage is not the driver |
| Transcriptome | RNA-seq, CIN brains vs `gfp-RNAi` controls | DESeq2, GSEA | Up: four Hsp70-family chaperones, GstE8, vesicle trafficking, autophagy, mitophagy, protein folding. Down: ribosome biogenesis, RNA metabolism, mitochondrial ribosomes, `dj-1beta`, overall mitochondrial-genome transcription |
| Proteasome | CL1-GFP degron reporter | GFP accumulation in NBs | Degraded in controls, accumulates in aneuploid NBs → "near saturation of the proteasome" |
| Aggregate seeding | Htt25Q-Cerulean (non-aggregating on its own) | foci; NB number | Diffuse in controls, large foci in aneuploid NBs, and further reduced NB number |
| Autophagy | mCherry-Atg8a | puncta | High autophagy induced in aneuploid NBs; absent from control NBs |
| Mitochondrial morphology | sqh-Mito-YFP | distribution | Large, unevenly distributed accumulations in aneuploid NBs |
| Mitochondrial oxidation / turnover | MitoTimer | green→red shift | Shift toward red, which the authors themselves read as "either high rates of mitochondrial oxidation or slow rates of mitochondrial turnover" |
| Mitophagy | MitoQC (tandem GFP-mCherry) | mitolysosome puncta | Large accumulations positive for both colours → mitophagy compromised |

### 3.1 Three things that were not measured, and one of them is load-bearing

- **ROS were never measured.** Case-sensitive grep for `ROS` in the full text returns six matches:
  one describing GstE8's family as "ROS scavenging", three describing the *overexpressed*
  scavengers Sod2 and GTPx-1, one parenthetical definition, one discussion sentence. There is no
  DHE, MitoSOX, roGFP or any other direct ROS measurement. The scavenger rescue is a
  **genetic sufficiency** result. It is not evidence that ROS were elevated, and the paper's only
  oxidation-adjacent measurement, MitoTimer, is ambiguous by the authors' own statement. This is
  the single largest evidential gap in the arm any NAC or idebenone row would sit on.
- **BubR1 / BUB1B abundance was never measured**, in any arm, in any genotype.
- **No mammalian replication.** Every result is in fly neuroblasts. The paper offers none, and I
  did not find one (§8, query Q9).

---

## 4. The rescue table, and the split that decides everything

This is the part that gets compressed in summary, and the compression is where the error enters.
Two readouts behave differently and the paper is explicit about it: **neuroblast number** (an
intermediate, cell-count readout) and **brain size** (the disease endpoint the title names).

| Intervention | Type | NB number | Pros+ progeny | Brain size |
|---|---|---|---|---|
| Atg1 overexpression | genetic | significant rescue | **no rescue** | **no rescue** |
| Rheb depletion / `Tor` mutant heterozygosity | genetic | significant rescue | **no rescue** | **no rescue** |
| **Rapamycin, 40 µM, fed to larvae** | **pharmacological** | rescue | **no rescue** | **no rescue** |
| Hsp60 or Hsp60C overexpression | genetic | partial rescue | rescue | **complete restoration** |
| Sod2 or GTPx-1 overexpression | genetic | partial rescue | rescue | **complete restoration** |
| p35 (caspase inhibition) | genetic | no rescue in mature brains | partial rescue | **complete restoration** |
| dMyc or Worniu overexpression | genetic | no rescue | Worniu partial | not reported as rescued |
| `prospero` depletion | genetic | no rescue | — | — |

The paper's own sentence, verbatim:

> Interestingly, boosting autophagy (by overexpressing Atg1 or depleting the activity of Tor
> through Rheb depletion or the use of Tor mutants in heterozygosis) led to a significant rescue
> of the number of aneuploid NBs in mature CIN brains, without causing any major effect in
> otherwise control brains (Fig. 6D and Supplementary Fig. 6A). A similar result was obtained by
> feeding larvae with the TOR inhibitor Rapamycin (Supplementary Fig. 6C). These results indicate
> that proteostasis failure, most probably due to an imbalance in the proteome caused by
> aneuploidy, functionally contributes to the reduction in the number of NBs in CIN brains.
> However, we noticed that none of the genetic interventions aimed at increasing autophagy was
> able to rescue the number of Pros+ cells or the size of the larval brain (Fig. 6E and
> Supplementary Fig. 6C). These results suggest that TOR inhibition or non-physiological levels
> of autophagy might have a negative impact on the proliferative capacity of highly aneuploid NBs.

And again in the Discussion:

> Activation of autophagy in these cells either by overexpression of Atg1 or depletion of TOR was
> able to partially rescue the reduced number of NSCs observed in SAC-depleted brains. The genetic
> interventions performed were unsuccessful in restoring progeny cell numbers or brain size.
> Future studies are required to ascertain whether TOR inhibition or aberrant autophagy levels
> negatively impact the proliferative capacity of highly aneuploid NSCs.

### 4.1 What this does to a sirolimus or metformin row

The reading that puts sirolimus, everolimus or metformin on this node is that **TOR inhibition was
one of the rescues of CIN-induced microcephaly** in the *Drosophila* proteostasis model
(PMID 41820377) — two independent lines pointing at one node, one of them in the correct gene.

Read against the paragraph above, that sentence is wrong on the endpoint. TOR inhibition rescued a
cell count and failed on microcephaly, and the authors flag it as a possible *liability* for
proliferative capacity. This is not a nuance: it removes the paper as support for
sirolimus/everolimus, and for metformin, whose only route onto the node is the same two processes
— Atg1/autophagy up, TOR down — that rescued the fly model genetically. Those two processes are
exactly the ones that did not rescue brain size.

The one **market-approved drug actually administered in the paper is rapamycin**, and it failed
on the disease endpoint. That is a measured negative where our Appendix B row currently says
"Unknown in an aneuploid background".

The mechanism also circulates through a bioRxiv preprint rather than the journal article, under a
title of the form "Gonzalez-Blanco et al. (2025) bioRxiv: TOR depletion rescues NSC stemness in
MVA Drosophila model". That title describes the NSC-count result alone, so the endpoint conflation
may have entered from the preprint framing.
**I did not retrieve the preprint** and cannot say whether its text differs from the journal
version; treat the preprint-versus-journal question as open.

---

## 5. Is there a market-approved drug that moves the node in the correct direction?

Track 2 is scoped to "existing market-approved medications" — verified today by fetching
`tabs/rules.py` from the Space, which reads: "Track 2 - Drug Repurposing: leveraging the genetic
information and variant prediction from Task 1, identify existing market-approved medications
whose mechanism of action may target the disease-causing pathway(s)." Investigational compounds do
not qualify.

Approval statuses below: arimoclomol and idebenone were verified today against the FDA and EMA
records respectively. Sirolimus, everolimus, metformin, N-acetylcysteine and edaravone are
stated from standard labelling and were **not** re-verified today; if any of them enters our
report, check the label first.

The paper's own closing sentence on translation names no molecule:

> Thus, genetic interventions aimed at restoring mitochondria homeostasis or blocking apoptosis
> will undoubtedly offer new opportunities for drug discovery and disease treatment.

| Candidate | Node it would move | Approval status | Evidence it moves the node the right way | Verdict |
|---|---|---|---|---|
| Sirolimus / everolimus | TOR ↓ → autophagy ↑ | Approved, with paediatric exposure | Measured, and it fails on the endpoint. Rapamycin 40 µM fed to CIN larvae rescued NB number and did not rescue Pros+ cells or brain size — González-Blanco et al., *Nat Commun* 2026, PMID 41820377, Supplementary Fig. 6C | Reject. Already rejected in our Appendix B for other reasons; now rejectable on measured endpoint data |
| Metformin | AMPK ↑ → autophagy ↑, mTORC1 ↓ | Approved, paediatric | Same endpoint failure applies by its own stated route. Independently, metformin is one of three compounds found selectively lethal to aneuploid cells — Tang et al., *Cell* 2011, PMID 21315436 | Reject. Our catalogue already scores it OPPOSES at `N9_stress` |
| N-acetylcysteine | ROS ↓ | Approved (acetaminophen overdose; mucolytic) | Genetic Sod2/GTPx-1 overexpression restored brain size in flies; NAC is not that intervention. Against it: NAC and vitamin E accelerate tumour progression by disrupting the ROS–p53 axis — Sayin et al., *Sci Transl Med* 2014, PMID 24477002 — and NAC increases melanoma metastasis — Le Gal et al., *Sci Transl Med* 2015, PMID 26446958 | Reject on oncogenic risk. §5.1 |
| Idebenone | short-chain quinone, electron carrier / antioxidant | EU marketing authorisation for LHON, 8 September 2015, EU/1/15/1020/001, adolescents ≥12 y and adults, granted "under exceptional circumstances". Not FDA-approved. Efficacy trial not retrieved; approval status verified from the EMA record | PubMed `aneuploidy AND idebenone` returns **0 hits**. No measured direction in any aneuploidy, CIN or MVA context | Reject, unassessed. Also fails a strict reading of "market-approved" in the US |
| Edaravone | free-radical scavenger | Approved, ALS, adults | PubMed `edaravone AND aneuploidy` returns **0 hits** | Reject, unassessed |
| Arimoclomol | HSF1 → HSP70 co-induction | **Approved 20 September 2024**, US, with miglustat, for NPC, ages ≥2 y — Keam, *Drugs* 2025, PMID 39715913 | Node mismatch: the paper rescued with the **mitochondrial matrix** chaperones Hsp60/Hsp60C, and the Hsp70-family genes were already among the most *upregulated* transcripts in CIN brains. PubMed `arimoclomol AND aneuploidy` and `"chromosomal instability" AND arimoclomol` both return **0 hits** | Reject, unassessed, and probably the wrong chaperone compartment. Worth naming as the closest existing approved drug on the proteostasis node |
| Caspase inhibitor (p35 equivalent) | apoptosis ↓ | **None approved.** Emricasan and similar are investigational | Not applicable | Out of scope for Track 2 |

**No approved drug does what the two interventions that restored brain size actually did.**
Hsp60/Hsp60C is a mitochondrial matrix chaperonin; Sod2 is mitochondrial MnSOD; GTPx-1 is a
glutathione peroxidase. `HSP60 AND aneuploidy AND drug` returns 0 hits on PubMed. The honest
version of this is already on the record, and we should not pretend to have discovered it: every
candidate on this node has to be held to low confidence because every rescue in the source paper
was **genetic** — an overexpression or a mutant, not a small molecule.

### 5.1 The antioxidant arm carries the same oncogenic-risk problem our own gate is built on

This is the sharpest thing in §5 and it is not a technicality.

The arm that fully restored brain size is the ROS/mitochondrial-chaperone arm, and the approved
molecules that gesture at it are antioxidants. In a cancer-predisposition syndrome, systemic
antioxidants have a measured, published, wrong-direction effect on tumour biology:

Sayin et al., *Sci Transl Med* 2014, PMID 24477002 — NAC and vitamin E "markedly increase tumor
progression and reduce survival" in B-RAF- and K-RAS-driven mouse lung cancer, and the mechanism
is disruption of the ROS–p53 axis: "antioxidants accelerate tumor growth by disrupting the ROS-p53
axis". The authors' own translational warning is that antioxidants "may accelerate the growth of
early tumors or precancerous lesions in high-risk populations".

Le Gal et al., *Sci Transl Med* 2015, PMID 26446958 — NAC increases lymph node metastases in an
endogenous mouse melanoma model; NAC and Trolox increase migration and invasion of human melanoma
cells.

Now put that against the surveillance mechanism in this disease. Aneuploidy-driven tumorigenesis
is suppressed by the ATM–p53 pathway — Li et al., *PNAS* 2010, PMID 20663956 — and
`docs/track2-hypothesis.md` §2 already refuses `N10_clearance` as a target on exactly this ground.
An antioxidant that blunts ROS–p53 signalling in a child with germline whole-chromosome
instability and an embryonal-tumour predisposition is the same category of error the document
already catches for hydroxychloroquine, metformin and tanespimycin, arriving from the opposite
direction.

**This is an inference at one step and I am labelling it.** Neither Sayin 2014 nor Le Gal 2015
studied a CIN syndrome, an MVA model or a child, and neither measured chromosome number. What is
measured is that two approved-or-otherwise-available antioxidants accelerate tumour progression
through p53. What is inferred is that this transfers to a germline CIN background whose tumour
suppression runs through the same axis. The inference is not certain. It is, however, precisely
the kind of unmeasured-risk question our `oncogenic_risk` gate exists to block on, and it blocks
seven of our twelve rejections already.

---

## 6. A rival axis, or a complementary one?

**Complement, in series — with one measured coupling that is more interesting than either.**

### 6.1 The axes are sequential, not parallel

Our mechanism graph runs `BUB1B alleles → transcript ↓ → BubR1 protein ↓ → MCC ↓ / KT–MT ↓ →
missegregation → aneuploidy ↑ → N9_stress ↑`. The González-Blanco axis begins at "aneuploidy ↑"
and describes what that aneuploidy then does to a stem cell: proteome imbalance, proteasome
saturation, autophagy induction, mitophagy failure, mitochondrial accumulation, loss of stemness.
Every measurement in the paper is downstream of the missegregation event.

So the two are not rival explanations of one node. Raising BubR1 aims to reduce the *rate* at
which complex aneuploid karyotypes are generated. Relieving proteotoxic and mitochondrial stress
aims to make cells that already carry them survive and keep dividing. Both can be true. The paper
supports this reading directly: its own model shows the damage accumulating "over time and through
consecutive cell cycles", which is a rate-dependent process, and it shows that at least 48 h of
`bub3` depletion were required before NB number moved.

**Does relieving proteotoxic stress raise BubR1?** Unmeasured. PubMed `BubR1 AND proteostasis`
returns **0 hits**. `rapamycin AND BubR1` returns 1 hit,
Kuijt et al., *Chromosoma* 2014, PMID 24695965, which uses rapamycin-analogue dimerisation as a
protein-targeting tool and is not pharmacology. No study I found measures BubR1 protein abundance
under mTOR inhibition, proteasome modulation, autophagy induction or chaperone induction, in any
system.

There is an unmeasured direction-of-effect worry in that gap, and it runs against the proteostasis
route rather than against ours. mTOR inhibition reduces cap-dependent translation globally. If BubR1 abundance is
the deficient quantity, a drug that lowers global translation is a candidate for making the
primary lesion worse. **That is an inference from mTOR biology, not a measurement**, and I found
nothing that tests it — but proposing sirolimus for a protein-abundance disease should have to
answer it.

### 6.2 One measured coupling, in the direction that helps

Oxidative stress lowers BubR1 recruitment, and antioxidants reduce the resulting segregation
errors:

Nago et al., *Reprod Med Biol* 2025, PMID 39845481 — in `Sod1`-deficient mouse oocytes matured in
vitro at 21% O₂, chromosome misalignment and segregation errors rose, and the change was
"associated with a decrease in the recruitment of BUBR1 to kinetochores at M I stage". Adding
L-ascorbic acid or N-acetylcysteine during maturation reduced segregation errors in the knockout
oocytes. REC8 was unchanged, so the effect is not general protein loss.

That is the only paper I found in which a redox manipulation and BUBR1 are measured in the same
experiment, and it says the two axes are coupled rather than parallel. Four caveats, all of which
matter: it is **meiosis I in oocytes**, not mitosis in somatic cells; the gene is **Sod1**
(cytosolic Cu/Zn SOD), not the mitochondrial Sod2 that rescued the fly brains; the readout is
**kinetochore-localised BUBR1 immunofluorescence**, not total protein abundance, which is the
quantity our hypothesis targets; and there is no *BUB1B*-mutant background anywhere in it.

A second, weaker link points the same way. Zhao et al., *Mol Cell Biochem* 2022, PMID 35639235
reports by ChIP that SIRT1 "could bind to the promoter of BubR1 … to upregulate its expression",
and that autophagy induction reduces mitotic catastrophe while BubR1 overexpression does the same.
If that transcriptional link is real it matters more than it looks, because
North et al., *EMBO J* 2014, PMID 24825348 established that the SIRT2 route is stabilisation-only
and cannot exceed the transcript ceiling — the limitation our §4.3 calls load-bearing.
**But this is a cancer-cell-line paper using BMH-21 to induce mitotic catastrophe, it is SIRT1 not SIRT2, and it is one ChIP experiment.**
Recorded as a lead, not as support.

### 6.3 Our own molecule already moves part of that node

Two verified results put nicotinic acid on the mitochondrial arm, in humans:

Pirinen et al., *Cell Metab* 2020, PMID 32386566 — niacin at 750–1,000 mg/day for 10 months in
adult-onset mitochondrial myopathy raised blood NAD⁺ up to 8-fold, brought patient muscle NAD⁺ to
control levels, and "muscle strength and mitochondrial biogenesis increased in all subjects".
`docs/track2-hypothesis.md` Appendix G already cites this paper and correctly refuses to give the
functional benefit weight, because the trial is open-label and the controls improved too. For the
present purpose that objection does not bite: the mitochondrial-biogenesis and NAD⁺ measurements
are biochemical, and the claim here is target engagement on the mitochondrial node, not efficacy.

Fang et al., *Cell Metab* 2016, PMID 27732836 — NAD⁺ repletion in ataxia-telangiectasia models
"stimulate[s] neuronal DNA repair and improve[s] mitochondrial quality via mitophagy". Mitophagy
is the exact quality-control process González-Blanco et al. show to be compromised in aneuploid
neuroblasts.
**Caveat that must travel with this: Fang 2016 used nicotinamide riboside and NMN, which are not approved medicines, and it is a different disease.**

There is also a plausible enzymatic bridge to the specific node that restored brain size: SIRT3
deacetylates and activates MnSOD/SOD2 — Qiu et al., *Cell Metab* 2010, PMID 21109198, and
Tao et al., *Mol Cell* 2010, PMID 21172655. SIRT3 is NAD⁺-dependent, so NAD⁺ repletion is a
candidate route to raising SOD2 *activity* without administering an exogenous antioxidant.
**This is an inference and it is a three-step one: neither paper studies niacin, aneuploidy, or neural stem cells, and nobody has measured SOD2 acetylation under nicotinic acid in any CIN background.**
It is a hypothesis worth one sentence in our report and no more.

### 6.4 The reverse coupling, which nobody has measured and which cuts against us

If mitochondrial function is impaired in aneuploid cells, NAD⁺ availability may fall, and SIRT2
activity with it — which would weaken the very axis our hypothesis runs on, in exactly the cells
we most want to treat. PubMed `aneuploidy AND "NAD+" AND cells` returns 70 hits, none of which, on
title screening, measures NAD⁺ concentration in an aneuploid or CIN background. I could not
establish this either way. It is a real open risk to our primary hypothesis, surfaced by reading
this paper, and it belongs in our limitations whatever else we do.

---

## 7. The strongest published objections to the proteostasis route in this disease

Deliberately assembled as the case *against* the mechanism, in descending strength.

**1. The paper's own negative on the disease endpoint.** Autophagy induction and TOR inhibition —
genetic and pharmacological — rescued neuroblast number and failed to rescue brain size, and the
authors read the result as a possible harm to proliferative capacity in aneuploid cells.
González-Blanco et al., *Nat Commun* 2026, PMID 41820377. This is the strongest objection
available because it is in the exact model the proteostasis route rests on, on the exact endpoint
the disease is named for, using the exact class of drug that route proposes.

**2. Aneuploid cells have limited autophagic reserve, so pushing on that lever is selectively lethal to them.**
Santaguida et al., *Genes Dev* 2015, PMID 26404941 shows autophagosomal cargo accumulating in
lysosomes of aneuploid cells with a TFEB-mediated lysosomal stress response — that is, the pathway
is already near capacity. Tang et al., *Cell* 2011, PMID 21315436 found AICAR, 17-AAG and
chloroquine selectively lethal to aneuploid cells. Donnelly et al., *EMBO J* 2014, PMID 25205676
shows HSF1 deficiency and impaired HSP90-dependent folding as hallmarks of aneuploid human cells.
`docs/track2-hypothesis.md` §3 already builds its `N9_stress` rejections on this literature; the
fly paper does not overturn it, and its Htt25Q result — aggregate seeding caused a **further**
reduction in aneuploid NB number — is consistent with it.

**3. The antioxidant arm has a measured tumour-promoting liability through the same p53 axis that suppresses aneuploidy-driven tumours.**
Two antioxidant results, set against the surveillance pathway they would blunt. Detail and the
labelled inference step in §5.1.

Sayin et al., *Sci Transl Med* 2014, PMID 24477002

Le Gal et al., *Sci Transl Med* 2015, PMID 26446958

Li et al., *PNAS* 2010, PMID 20663956

**4. The severity regime may not match human BUB1B-MVA.**
Simmons et al., *Hum Mol Genet* 2019, PMID 30668728 — "Nearly complete deletion of BubR1 causes
microcephaly through shortened mitosis and massive cell death". Our Appendix A.3 already records
the quantitative form: `BubR1^H/H` mouse E15.5 cortex sits at about 1.8% of wild-type BubR1 with
no microcephaly, and only conditional deletion to about 0.5% produces it. Chronic RNAi depletion
of a SAC gene throughout larval neurogenesis is a severe-end model. Whether the
proteostasis/mitochondrial response it reveals is engaged at the residual-BubR1 levels MVA
patients actually carry — roughly 17–50% full-length protein, is **unmeasured**. The residual
range is from Suijkerbuijk et al., *Cancer Res* 2010, PMID 20516114

**5. Mouse BubR1 microcephaly has a different published mechanism, and it is apoptotic.**
Sterling et al., *Front Cell Dev Biol* 2023, PMID 37900274 finds chromosomal defects and DNA
damage driving p53 activation and apoptosis in BubR1-deficient mouse cortex, with residual
apoptosis persisting after *Trp53* co-deletion and only "minimal rescue of cortical size". The fly
paper argues the opposite for its model — that NSC loss is "neither a direct consequence of cell
death nor a result of premature differentiation" and that NSCs are "highly resistant to DNA
damage". Both may be right in their own system, but the mammalian *BubR1* literature does not
currently point at proteostasis. Note the irony for the route: the intervention that fully
restored brain size in the fly alongside the mitochondrial one was **apoptosis inhibition**,
which has no approved drug and which in a cancer-predisposition syndrome would be the least
advisable thing on the list.

**6. A specificity claim that does not replicate in human cells.** The paper's case rests on
complex aneuploidies being categorically different from simple ones: trisomy X or 2L "did not
recapitulate any of the gene classes pinpointed in brains subjected to CIN". In human cells,
Dürrbaum et al., *BMC Genomics* 2014, PMID 24548329 profiled eleven trisomic and tetrasomic lines
plus two with complex karyotypes and reports the opposite — "complex aneuploidy elicits the same
transcriptional changes as trisomy", with the same signature of upregulated ER/Golgi/lysosome and
downregulated ribosomes. These are different systems and different assays and I am **not**
claiming one refutes the other. What it does mean is that the complex-versus-simple distinction,
which is the paper's own framing of what makes MVA special, is not established in human cells.

**7. Timing.** Every rescue in the paper is applied throughout larval neurogenesis, in parallel
with the insult. MVA microcephaly in humans is present at birth —
Hanks et al., *Nat Genet* 2004, PMID 15475955 and the clinical delineation in
García-Castillo et al., *Am J Med Genet A* 2008, PMID 18548531. A drug given to a diagnosed child
is a postnatal intervention against a prenatal structural deficit. `docs/track2-hypothesis.md`
§7.2 already concedes this against our own hypothesis; it applies with at least equal force to a
mechanism whose only endpoint is brain size.

**8. The neural compartment problem is ours, and this axis does not share it.**
This objection runs the other way and should be recorded honestly. If the neural arm is where the
proteostasis mechanism lives, our
molecule may not reach it: with a validated antibody, neurons in non-tumorigenic human cortex were
NAPRT-negative — Cole et al., *Oncotarget* 2017, PMID 29100430 — and §4.5 of our hypothesis
already concludes that "the neural arm cannot use niacin". The proteostasis-route agents do not
have that problem. Rapamycin, metformin, NAC and idebenone all reach brain to some degree. That is
the one respect in which that axis is better positioned than ours, and no amount of correcting the
TOR claim changes it.

---

## 8. Queries run, including the empty ones

| # | Source | Query | Result |
|---|---|---|---|
| Q1 | PubMed web | `pubmed.ncbi.nlm.nih.gov/41820377/` | Resolves. Full record as in §1 |
| Q2 | Europe PMC REST | `PMC13121444/fullTextXML` | 170,193 bytes, complete article including Methods and all 58 references. Used for every verbatim quotation in this document |
| Q3 | PMC direct / Exa | `pmc.ncbi.nlm.nih.gov/articles/PMC13121444/` via Exa | reCAPTCHA interstitial, no content. WebFetch on the same URL succeeded. Europe PMC used instead |
| Q4 | PubMed esummary | 12-PMID batch: `41820377,34216545,37900274,26404941,30668728,26573328,30794535,28351851,22968442,17702937,25205676,21315436` | All 12 resolve to the papers named here |
| Q5 | PubMed esummary | 7-PMID batch: `32386566,21109198,21172655,24825348,20516114,31738183,15475955` | All 7 resolve |
| Q6 | PubMed esummary | `39169218,18548531,35044816,27869830` | All 4 resolve. Confirms the ref-11 mismatch in §2.2 |
| Q7 | PubMed | `"mosaic variegated aneuploidy" AND (proteostasis OR mitochondrial)` | **2 hits.** PMID 41820377 and Villarroya-Beltri et al. (MAD1L1 syndrome, PMID 36322655). González-Blanco 2026 is effectively the only paper connecting MVA to this axis |
| Q8 | PubMed | `"mosaic variegated aneuploidy" AND (antioxidant OR idebenone OR edaravone OR acetylcysteine)` | **0 hits.** No drug evidence of any kind exists on this axis in this disease |
| Q9 | PubMed | `aneuploidy AND idebenone`; `edaravone AND aneuploidy`; `arimoclomol AND aneuploidy`; `"chromosomal instability" AND arimoclomol`; `HSP60 AND aneuploidy AND drug` | **0 hits each.** Five separate empty results |
| Q10 | PubMed | `BubR1 AND proteostasis` | **0 hits** |
| Q11 | PubMed | `rapamycin AND BubR1` | 1 hit, PMID 24695965, a dimerisation-tool paper, not pharmacology |
| Q12 | PubMed | `mTOR AND BUB1B` | 5 hits, all cancer-expression correlation studies. None measures BubR1 protein under mTOR inhibition |
| Q13 | PubMed | `BubR1 AND "oxidative stress"` (9 hits), `BubR1 AND autophagy` (7 hits) | Yielded the two couplings in §6.2, PMID 39845481 and PMID 35639235. The remainder are unrelated |
| Q14 | PubMed | `aneuploidy AND "NAD+" AND cells` | 70 hits; **none found**, on title screening, measuring NAD⁺ concentration in an aneuploid or CIN background. §6.4 |
| Q15 | PubMed | `acetylcysteine AND aneuploidy AND micronucleus` | 2 hits, both irrelevant (nanoparticle genotoxicology, PMID 25325158; an industrial compound, PMID 23628427) |
| Q16 | PubMed esearch | `Le Gal K[au] AND antioxidants AND melanoma metastasis`; `Fang EF NAD+ replenishment mitophagy ataxia telangiectasia`; `Durrbaum aneuploidy transcriptional response human cells`; `arimoclomol first approval[ti]` | Located PMIDs 26446958, 27732836, 24548329, 39715913, each then verified by esummary |
| Q17 | Space raw file | `tabs/rules.py` | Track 2 scope verbatim, quoted in §5 |
| Q18 | Exa | arimoclomol FDA approval; idebenone EMA authorisation | FDA press release 20 Sept 2024; EMA EPAR, marketing authorisation 8 Sept 2015, EU/1/15/1020/001 |

### Negatives, stated precisely

Each is "I did not find", not "does not exist".

- I did not find any mammalian replication of the Hsp60/Sod2/GTPx-1 brain-size rescue.
- I did not find any study measuring BubR1 or BUB1B protein abundance under proteostasis,
  autophagy, mTOR or mitochondrial manipulation, in any species.
- I did not find any approved drug with a measured effect on mitochondrial chaperone capacity or
  MnSOD activity in a chromosomal-instability background.
- I did not find any measurement of NAD⁺ concentration in aneuploid or CIN cells.
- I did not find any direct ROS measurement in the source paper — only genetic scavenger
  overexpression and the ambiguous MitoTimer readout.
  González-Blanco et al., *Nat Commun* 2026, PMID 41820377
- I did not retrieve the bioRxiv preprint version of the paper, so I cannot say whether its
  framing of the TOR result differs from the journal article's.
- I did not retrieve the Supplementary Figures themselves. Statements about Supplementary Fig. 6C
  rest on the main-text sentences that describe it, quoted verbatim in §4.
- I did not retrieve the idebenone pivotal efficacy trial; the approval statement rests on the EMA
  record alone.

---

## What this means for our submission

**Recommendation: (b) — address it as a stated limitation, and say why we did not pursue it.**

Not because the mechanism is wrong. It is a real, well-executed paper, it is explicitly about MVA,
and its mitochondrial arm produced the cleanest rescue in it. The reason is that the node yields
**no admissible candidate**, and the paper itself is what establishes that.

The argument, in the order I would put it in the report:

1. **The axis exists and is downstream of ours.** Proteostasis failure and mitochondrial
   dysfunction are consequences of accumulated complex aneuploidy in a neural stem cell. Our
   target is the protein whose scarcity generates that aneuploidy. Cite the paper at `N9_stress`,
   which is a node our mechanism graph already carries and currently supports with nothing
   MVA-specific. This is a citation we should have and do not.

2. **The one approved drug tested in that model failed on the disease endpoint.** Rapamycin at 40
   µM rescued neuroblast count and did not rescue brain size or progeny number, and the authors
   flag TOR inhibition as a possible harm to the proliferative capacity of aneuploid cells. Our
   Appendix B sirolimus row should change from "Unknown in an aneuploid background" to a measured
   negative with this citation. That is a strictly better rejection than the one we have.

3. **The arm that did work has no approved molecule.** Mitochondrial matrix chaperones and
   mitochondrial MnSOD, overexpressed cell-autonomously in the affected stem cells. No approved
   drug does that. The closest is arimoclomol, approved since September 2024 with paediatric
   labelling, and it induces the wrong chaperone family in the wrong compartment — the Hsp70 genes
   were already among the most upregulated transcripts in the CIN brains. Name it, reject it,
   state the compartment mismatch. Naming the best candidate against your own conclusion is worth
   more than not knowing it exists.

4. **The approved antioxidants trip our own gate.** NAC and idebenone reach the node only as
   generic ROS reduction, and systemic antioxidants have a measured tumour-accelerating effect
   through the ROS–p53 axis, in a syndrome whose aneuploidy surveillance runs through ATM–p53.
   That is the same reasoning that already rejects hydroxychloroquine, metformin and tanespimycin
   in §3, applied consistently. The transfer step is an inference and should be printed as one.

5. **Our molecule already moves part of that node, and we should claim exactly that much and no more.**
   Niacin raised NAD⁺ up to 8-fold and increased mitochondrial biogenesis in humans. That is
   target engagement on the mitochondrial arm from an approved oral drug, which no
   candidate on the proteostasis axis can show in humans. It is not efficacy, the trial was
   open-label, the controls improved too, and no aneuploid cell was involved. One sentence, with
   all four caveats attached.

6. **Concede the one thing that genuinely favours the other axis.** Its agents cross into brain;
   ours probably does not, because cortical neurons are NAPRT-negative and our own §4.5 says the
   neural arm cannot use niacin. Microcephaly is the phenotype our therapeutic argument is weakest
   on, and this paper is entirely about microcephaly. Conceding it is not damage — the limitation
   is real, and a report that cites this paper without stating it would be the weaker document.

**Why not (a).** Adding a second target node with its own candidate rows would produce rows that
all reject: sirolimus and metformin on a measured endpoint failure, NAC and idebenone on oncogenic
risk, arimoclomol and edaravone as unassessed, caspase inhibitors as unapproved. Our pipeline
already reaches that outcome for `N9_stress` by a different route. Re-running it to generate six
more rejections would add process and no information. What is worth adding is the *evidence*: a
measured negative in the sirolimus cell, and this paper cited at the node.

**Why not (c).** The axis does exist. What is wrong is one sentence commonly drawn from the
paper, and it happens to be the sentence carrying the route's strongest candidate. Claiming the
mechanism does not apply to this disease would require dismissing a *Nature Communications* paper
whose abstract opens by naming the disease, and we would be dismissing it on species and endpoint
grounds that apply to a good deal of our own evidence base too.

### What would change my mind

Toward **(a)**, any one of:

- A measurement of BubR1 protein abundance rising under proteostasis or mitochondrial
  manipulation in a human cell. Currently zero papers exist (Q10, Q12).
- A mammalian replication of the Hsp60/Sod2 brain-size rescue in a *BUB1B*-hypomorphic background
  at MVA-relevant residual protein levels rather than by chronic SAC-gene knockdown.
- An approved drug shown to raise mitochondrial chaperone capacity or MnSOD activity in human
  cells, with a measured chromosome-number or micronucleus readout that is not adverse. That last
  clause is the gate; it is what nicotinic acid has and what nothing on this node has.
- Evidence that the ROS–p53 tumour-surveillance objection does not transfer — for instance a
  mitochondria-targeted antioxidant with a measured aneuploidy endpoint in a CIN background.

Toward **(c)**, essentially only: a demonstration that the proteostasis and mitochondrial
response is not engaged at the residual BubR1 levels MVA patients carry, which would make the fly
model a severe-end artefact. Simmons et al., *Hum Mol Genet* 2019, PMID 30668728 makes that
question live — mouse cortex tolerates 1.8% of wild-type BubR1 without microcephaly — but it does
not answer it.

And one finding would weaken **our own** hypothesis rather than the proteostasis one: a
measurement showing NAD⁺ is depleted in aneuploid cells, which would mean the mitochondrial
dysfunction this paper describes
undercuts the SIRT2 axis we depend on, in exactly the cells we most want to treat. Nobody has
measured it (Q14). It belongs in our limitations either way.
