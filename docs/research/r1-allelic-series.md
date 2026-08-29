# R1 — The allelic-series challenge to the dosage thesis

Scope: public literature only. Prepared 2026-08-28.

Note on attribution: the paper the reviewer calls "Jeganathan et al. 2019" is
**Sieben CJ, Jeganathan KB, Nelson GG, Sturmlechner I, Zhang C, van Deursen WH, Bakker B,
Foijer F, Li H, Baker DJ, van Deursen JM. "BubR1 allelic effects drive phenotypic
heterogeneity in mosaic-variegated aneuploidy progeria syndrome." *J Clin Invest*
2020;130(1):171–188. PMID 31738183. PMC6934189. doi:10.1172/JCI126863.**
Jeganathan is second author, not first. Published online 18 Nov 2019, print issue Jan 2020.
A corrigendum exists (*J Clin Invest* 2020;130(11):6188, PMID 33136097, jci.org/articles/view/144781);
it corrects a duplicated image panel in Figure 2E and the Methods text for lymphoma mitotic
index. It does not touch any claim discussed here.

---

## 1. Verdict

The reviewer is right that total BubR1 abundance does not fully predict phenotype severity,
and Sieben 2020 shows this cleanly: two mouse genotypes with indistinguishable total BUBR1
and near-identical tissue aneuploidy differ ~1.75-fold in median lifespan. But the reviewer's
third claim is wrong as stated — the "normal kinetochore accumulation" result is from a
*heterozygote* carrying a wild-type allele, measured with a total-BUBR1 antibody, and in the
disease-relevant biallelic genotype the same paper reports kinetochore BUBR1 accumulation is
*severely compromised*. The dosage thesis survives in a narrowed form: raising total BubR1
is the only intervention that has ever corrected an MVA-model phenotype in vivo
(Baker 2013, Nat Cell Biol, PMID 23242215) and the only one that has ever restored checkpoint
function to the L1012P protein specifically (Suijkerbuijk 2010, Cancer Res, PMID 20516114) —
so abundance is a defensible primary pharmacodynamic variable, but it is a *floor*, not a
sufficient predictor of severity, and it should not be sold as one.

---

## 2. What Sieben 2020 actually is

**Model system.** Mice, not human cells. Four *BubR1* alleles combined pairwise:

| Allele | What it is |
| --- | --- |
| `L1002P` | New knock-in. Mouse L1002 = human L1012. Models the human `BUBR1 3035T>C` (L1012P) missense allele from Hanks 2004 (Nat Genet, PMID 15475955). |
| `X753` | Pre-existing knock-in of human `2211insGTTA` frameshift; truncated, unstable protein. From Wijshake 2012 (PLoS Genet, PMID 23300461). |
| `H` | Pre-existing hypomorph, ~10% of normal wild-type BUBR1. From Baker 2004 (Nat Genet, PMID 15208629). Models the "no discernible mutation, low wild-type protein" allele seen in Japanese/Asian MVA patients. |
| `–` | Null. |

Background: mixed 129/Sv × C57BL/6.

**What was measured.** BUBR1 protein by Western blot (BD Transduction 612503, normalised to
α-tubulin or Ponceau S) in passage-5 MEFs and in testes, thymus, spleen, gastrocnemius, fat.
Transcript by qRT-PCR. Protein turnover by cycloheximide ± MG132. Aneuploidy three ways:
metaphase chromosome counts (MEFs, splenocytes), interphase FISH for chromosomes 4 and 7
(liver, lung, spleen, skeletal muscle, and others), and single-cell DNA sequencing of MEFs.
Mitotic function: colcemid-challenge SAC assay, monastrol-washout error correction, live H2B-mRFP
imaging, centrosome movement, spindle geometry, kinetochore BUBR1 immunofluorescence normalised
to centromere signal. Phenotype: lifespan, weight, cataract onset, kyphosis onset, grip strength,
muscle fibre cross-sectional area, treadmill work, echo-MRI body fat, adipocyte size,
isoproterenol cardiac stress test, DMBA-induced lung tumours, spontaneous tumour spectrum.
RNA-seq of gastrocnemius and inguinal adipose.

**Headline results.**

- `X753/L1002P` (the exact human patient genotype) is embryonic lethal in mice, dying before
  E13.5, though the human patient lived 3.6 months. Mice tolerate MVA alleles worse than humans.
- `H/X753` mice die within 18 hours of birth. `H/L1002P` mice are viable, median lifespan 343 d
  (vs 691 d wild-type), and show cancer predisposition, dwarfism, lipodystrophy, sarcopenia,
  cataracts, kyphosis and low cardiac stress tolerance.
- `H/H` mice: median lifespan 196 d.
- Heterozygotes `+/L1002P` and `+/–` get reduced tumour latency and mild lifespan reduction but
  little progeria. `+/X753` gets overt progeria (from Wijshake 2012).

---

## 3. The reviewer's three claims, checked

### Claim 1 — "similar total BubR1 reduction and similar aneuploidy, different phenotype severity" — **SUPPORTED**

The comparison is `BubR1^H/L1002P` vs `BubR1^H/H`.

Protein:
> "despite the significant differences in body size between BubR1^H/L1002P and BubR1^H/H mice,
> the level of BUBR1 protein reduction in MEFs and multiple tissues appeared to be very similar
> between the mutants, with significant differences observed only in the spleen (Figure 6E)."

Aneuploidy in MEFs:
> "Aneuploidy rates were also relatively high in BubR1^H/L1002P MEFs, reaching virtually the same
> levels as observed in BubR1^H/H MEFs"

Aneuploidy in tissues (the measure the authors themselves call most reliable):
> "FISH analysis for chromosomes 4 and 7 showed that rates of karyotypic abnormalities were
> prominently increased in a broad spectrum of tissues from both BubR1^H/L1002P and BubR1^H/H
> mice, with minimal differences observed between the 2 genotypes (Figure 9B…)"

And the authors' own summary sentence:
> "aneuploidy rate or degree alone does not clearly correlate with the severity of the progeroid
> phenotypes observed in BubR1^H/L1002P and BubR1^H/H MVA models, implying that other pathological
> events are key contributors or drivers."

Phenotype divergence: median lifespan 343 d vs 196 d; kyphosis onset 238 d vs 175 d;
isoproterenol lethality 50% after 18 injections vs 50% within 5 injections; growth retardation
milder.

**Three qualifications the reviewer omits.**

(a) The divergence is *partial*. Several phenotypes were the same severity:
> "Although several phenotypes of BubR1^H/L1002P are similar in severity to those of BubR1^H/H
> mice, including lipodystrophy, cataractogenesis, and carcinogen-induced tumor susceptibility,
> others are definitely milder, including growth retardation, cardiac stress sensitivity,
> lifespan shortening, and muscle wasting."

Cataract median onset was 161 d (`H/L1002P`) vs 168 d (`H/H`) — effectively identical.

(b) The divergence runs in the direction that *hurts* an aneuploidy-centric model more than a
dosage-centric one. In mitotic splenocytes, aneuploidy was 38% in `H/L1002P` vs 15% in `H/H` —
the mouse with *more* aneuploidy lived *longer*.

(c) The paper still endorses a dosage threshold explicitly:
> "although BUBR1 levels do not closely track with phenotypic severity, this study identifies
> several critical BUBR1 thresholds in vivo. A reduction of approximately 40%–90% in BUBR1
> protein levels results in tumor predisposition or mild premature-aging phenotypes, whereas a
> reduction of greater than 90% results in progeroid phenotypes or early postnatal lethality."

So the paper's position is: level sets the *regime*; allele identity modulates *severity within*
the regime.

There is a second, independent same-protein/different-outcome pair in the paper, which the
reviewer did not cite and which is arguably stronger for their case, at a matched *monoallelic*
level:
> "Skeletal muscle from BubR1^+/L1002P and BubR1^+/X753 mice had similarly reduced BUBR1 protein
> levels compared with BubR1^+/+ mice, indicating that differences in BUBR1 protein levels are
> unlikely to account for the differential impact of these MVA mutations on mTORC1 signaling."

`+/X753` gets mTORC1 hyperactivation and sarcopenia; `+/L1002P`, at the same BUBR1 level, has a
skeletal-muscle transcriptome with **3** differentially expressed genes vs wild-type.

### Claim 2 — "authors concluded allelic effects beyond protein level and aneuploidy contribute to heterogeneity" — **SUPPORTED, verbatim**

Abstract, final sentence:
> "Together, these data demonstrate that BUBR1 allelic effects beyond protein level and aneuploidy
> contribute to disease heterogeneity in both MVA patients and heterozygous carriers of MVA
> mutations."

Results, monoallelic section:
> "Overall, the above experiments indicate that allelic effects beyond BUBR1 levels alone
> contribute to the mitotic phenotypes of monoallelic MVA mutants."

Discussion:
> "The mere difference between total BUBR1 protein levels in BubR1^H/H and BubR1^H/L1002P mice is
> that in the latter strain the total pool consists of a mixture of BUBR1^WT and BUBR1^L1002P
> protein rather than BUBR1^WT alone."

Supplemental Table 5 is titled: "BubR1 MVA mono- and biallelic mice exhibit phenotypic diversity
independent of BUBR1 levels and aneuploidy rates."

This claim is not overstated. It is the paper's thesis.

### Claim 3 — "L1002P mutant protein accumulated normally at unattached kinetochores despite reduced total abundance" — **OVERSTATED / misattributed**

The sentence exists:
> "Nevertheless, mitotic BubR1^+/L1002P MEFs were still able to accumulate normal amounts of
> BUBR1 at unattached kinetochores, as demonstrated by immunofluorescence (IF) labeling for
> BUBR1 (Supplemental Figure 2F)."

Three reasons the reviewer's version does not follow.

**(i) Wrong genotype.** It is `BubR1^+/L1002P` — a heterozygote with one intact wild-type allele.
The paper says "BUBR1", not "BUBR1^L1002P". The antibody (BD Transduction 612503) does not
distinguish wild-type from mutant protein, and Supplemental Figure 2F's legend describes only
"Quantification of BUBR1 levels at the kinetochore in indicated MEFs, normalized to centromeric
signal." Nothing in the paper attributes the kinetochore signal to the mutant protein. A wild-type
allele supplying a normal kinetochore pool is a fully consistent — and simpler — reading.

**(ii) The disease-relevant genotype shows the opposite.** In `BubR1^H/L1002P`, where the
wild-type contribution is reduced to hypomorph levels:
> "SAC activity was severely compromised in both BubR1^H/L1002P and BubR1^H/H MEFs (Figure 8H),
> as was BUBR1 accumulation at unattached kinetochores at the onset of mitosis
> (Supplemental Figure 10, A and B)."

Supplemental Figure 10 legend confirms this is the same assay: "Quantification of BubR1 levels at
the kinetochore in BubR1^+/+, BubR1^H/L1002P and BubR1^H/H MEFs, normalized to centromeric signal."
So in the genotype that actually models the disease, kinetochore BUBR1 *tracks total abundance
downward*, in both mutants, equally.

**(iii) Human patient cells agree with (ii), not with the reviewer.** Suijkerbuijk 2010
(Cancer Res, PMID 20516114):
> "Low overall protein abundance translated to low levels of BUBR1 on unattached kinetochores
> (Fig. 2B). Enhancing contrast settings, however, showed that the residual protein could be
> recruited to kinetochores, indicative of at least partial functionality of mutant BUBR1."

That is: *less* kinetochore BUBR1 in proportion to total, with the residual protein still
competent to localise. Which supports the dosage reading — the mutant protein is not
qualitatively broken at localisation, there is just too little of it.

**Verdict on Claim 3:** the quoted observation is real but is about a heterozygote and about total
BUBR1. The inference "kinetochore-localised BubR1 is normal while total abundance is low,
therefore abundance is the wrong variable" is not supported by Sieben 2020 and is contradicted by
its own Supplemental Figure 10.

---

## 4. Other evidence bearing on "is MVA1 reducible to total BubR1 abundance?"

### 4a. Evidence FOR abundance as the dominant causal variable

**Suijkerbuijk SJ, van Osch MH, Bos FL, Hanks S, Rahman N, Kops GJ. Cancer Res 2010;70:4891–4900.
PMID 20516114. PMC2887387.** The single most important paper for our thesis, and it tests the
exact allele the reviewer cites.

> "forced overexpression of the poorly expressed substitution mutants I909T and L1012P to levels
> comparable to wild-type BUBR1 fully restored the response to nocodazole (Fig. 4D, upper graph).
> This showed that these mutations do not impose significant constraints on BUBR1 function other
> than affecting overall BUBR1 protein abundance."

Also from the same paper:
> "Ectopic expression of BUBR1 restored mitotic checkpoint activity, proving that BUBR1
> dysfunction causes chromosome segregation errors in the patients."

> "all MVA mutations fall in two distinct classes: those that impose specific defects in checkpoint
> activity or microtubule attachment and those that lower BUBR1 protein abundance."

And a threshold, by doxycycline titration of a BUBR1 shRNA in HeLa:
> "Whereas saturating amounts of doxycyline reduced BUBR1 protein to ~6% and caused chromosome
> mis-segregation in the majority of cells, low amounts of doxycyline, while still significantly
> reducing BUBR1 levels to ~13%, had virtually no impact on the fidelity of chromosome segregation."

Mechanism of the low abundance for kinase-domain missense alleles: HSP90-dependent folding,
proteasomal clearance. Geldanamycin removed the mutant protein; MG132 rescued it. This is a
tractable stabilisation target.

**Baker DJ, Dawlaty MM, Wijshake T, et al. Nat Cell Biol 2013;15(1):96–102. PMID 23242215.
PMC3707109.** This is the in vivo dosage-rescue
result, and it is easy to miss because it is reported as a transgene validation control:

> "Expression of Flag-BubR1 corrected all premature aging phenotypes of BubR1^H/H mice
> (Fig. 1c, and data not shown), confirming that Flag-BubR1 protein was functional and adequately
> expressed."

Figure 1c legend: "Both Flag-BubR1 transgenes correct the growth retardation and aging-associated
phenotypes of BubR1^H/H mice."

So: raising total BubR1 in the reference MVA progeria mouse corrects the phenotype. That is the
strongest existing support for "raise it above a threshold and the phenotype corrects." Caveat:
`BubR1^H/H` is a pure hypomorph — the pool is all wild-type protein. It does not test whether
raising a mixed WT + mutant pool works.

The same paper also shows the up-direction of the dose-response in wild-type animals: sustained
high BubR1 reduces aneuploidy, reduces tumours, extends median and maximum lifespan. Effect was
tissue-selective — muscle, heart, kidney, retina improved; osteoporosis and cataractogenesis
were not delayed.

**Bohers E, Sarafan-Vasseur N, Drouet A, et al. Hum Genet 2008;124(5):473–478. PMID 18932004.** A dose-response series in HeLa via
shRNA: residual BUBR1 at 8.5, 10, 14, 58 and 77%.
> "Induction of PCS was detected in all transduced cells and its level was correlated to the
> decrease of BUBR1. Aneuploidy was clearly detected in cells with residual BUBR1 below 50%.
> Our data demonstrate that the function of the human BUBR1 protein in the spindle checkpoint is
> remarkably dosage-dependent"

**Matsuura S, et al. Am J Med Genet A 2006;140:358–367. PMID 16411201.** Human patient cells.
Reduced BubR1 protein *and* diminished kinetochore localisation; Cdc20 kinetochore association
abolished.
> "Microcell-mediated transfer of chromosome 15 (containing BUB1B) into the cells restored normal
> BubR1 levels, kinetochore localization of p55cdc, and the normal responses to colcemid treatment."
> "…suggest that >50% decrease in expression (or activity) of BubR1 is involved in the PCS syndrome."

A gene-dosage rescue in human patient cells.

**Ochiai H, et al. PNAS 2014;111:1461–1466. PMID 24344301.** The second allele in the Japanese
PCS/MVA families is a purely regulatory intergenic G>A 44 kb upstream of a *BUB1B* transcription
start site. Editing it biallelically into cultured human cells reproduced reduced *BUB1B*
transcript, increased PCS and MVA. A pure expression-level lesion, no coding change, sufficient
to cause the syndrome. Same for the Alu/SVA insertion cases: Kato M, et al. Hum Genome Var
2017;4:17021, PMID 28611924; and Kawakami R, et al. J Hum Genet 2024;69(2):91–99, PMID 38102195,
which reports "immunoblotting showed severely diminished BUB1B protein level in patient derived cells."

**Baker DJ, et al. Nat Genet 2004;36:744–749. PMID 15208629** (as summarised in Sieben 2020's
introduction): a graded protein-level allelic series with sharp viability thresholds —
0% (`–/–`) dies by E3.5; ~5% (`–/H`) dies by P1; ~10% (`H/H`) viable but progeroid.

### 4b. Evidence AGAINST reducing MVA1 to total abundance

**Sieben 2020 itself** (section 3 above), plus one more finding: the mutant protein does
something the absence of protein does not. Comparing `H/L1002P` and `H/H` at matched aneuploidy:
> "The only notable difference between the 2 genotypes was the dramatic difference in PCS
> (Figure 8A), indicating that the presence of BUBR1^L1002P protein interferes with the cell's
> ability to sustain strong bonds between duplicated chromosomes before anaphase onset."

That is an interference effect of the mutant protein, not a shortfall of the wild-type protein.
It argues that in a mixed pool, mutant protein is not merely inert filler.

**Not all MVA alleles are dosage alleles.** Suijkerbuijk 2010 identified a separation-of-function
class: Y155C (TPR domain) has *normal* protein level, localises normally to unattached
kinetochores, restores chromosome alignment, but fails to restore checkpoint activity. R550Q has
normal level and normal function. For a Y155C-carrying patient, raising abundance would not help.
The dosage argument applies to the kinase/pseudokinase-domain missense class (which includes
L1012P, I909T, R727C, R814H, L844F) and to the truncating and regulatory alleles — not to all.

**The pseudokinase domain has a real, non-structural function, discovered after Suijkerbuijk 2010.**
Gama Braga L, et al. Cell Rep 2020;33(7):108397. PMID 33207204: the C-terminal pseudokinase domain of
human BUBR1 is required to promote KARD phosphorylation and PP2A-B56 recruitment; mutation or
removal causes attenuated checkpoint silencing and chromosome alignment errors via Aurora B
imbalance. L1012 sits in this domain. Suijkerbuijk's rescue assay was a nocodazole checkpoint
response, so it would not necessarily surface a PP2A-B56-recruitment defect. **This is an
unresolved gap**: whether L1012P protein raised to normal abundance is fully functional in the
KARD/PP2A arm has not, in what I found, been tested.

**A second hypomorph gives a different phenotype.** Sieben 2020 notes that an independently
generated *BubR1* hypomorph, `BubR1^L/L` (Kyuragi R, et al. Arterioscler Thromb Vasc Biol
2015;35(2):341–347, PMID 25524773), "did not show any of the overtly detectable progeroid features of
BubR1^H/H mice, such as kyphosis, lipodystrophy, and cataract formation." I have this only as
cited in Sieben 2020; I did not read the Kyuragi paper's own quantification of BubR1 level, so the
implied "same level, different phenotype" is unverified. Mark unresolved.

**Non-mitotic BubR1 functions are real and could contribute to the phenotype.**

- *Insulin receptor endocytosis.* Choi E, Zhang X, Xing C, Yu H. Cell 2016;166:567–581.
  PMID 27374329. "MAD2 directly binds to IR and facilitates BUBR1-dependent recruitment of the
  clathrin adaptor AP2 to IR… BUBR1 deficiency enhances insulin sensitivity in mice." Follow-up:
  Park J, et al. Diabetes 2023;72(12):1781–1794. PMID 37725942. This is an interphase function.
  Sieben 2020 flags it as a candidate link to the mTORC1 hyperactivation in `+/X753` muscle.
- *Separable in vivo, in Drosophila.* Tang R, et al. Cell Rep 2020;31(12):107794. PMID 32579921.
  BubR1 kinase-dead knock-ins are "homozygously viable and fertile and show no defects in mitotic
  progression" but have attenuated insulin signalling and increased lifespan. A clean
  demonstration that a BubR1 function can be lost with mitosis intact.
- *Lipolysis.* Liu M, et al. Aging (Albany NY) 2024;16(4):3257–3279. PMID 38334966. Fat-body-specific
  BubR1 depletion in Drosophila increases lipid degradation via IMD/Relish. Relevant given
  lipodystrophy is a core MVA1 phenotype.
- *Ciliogenesis.* Miyamoto T, et al. Hum Mol Genet 2011;20(10):2058–2070. PMID 21389084.
  "BUBR1 is essential for the primary cilium formation, and… the PCS (MVA) syndrome is thus a
  novel ciliopathy." Proposed as the explanation for Dandy-Walker complex and polycystic kidneys.
  See also Miyamoto T, Matsuura S. Oncotarget 2015, PMID 26309087 (no abstract in
  PubMed; I did not read the full text).
- *Cardiac.* Pun R, et al. Aging Cell 2025;24(9):e70160. PMID 40607964 — BubR1 insufficiency drives
  cardiac hypertrophy, fibrosis and senescence, and BubR1 protein declines in human heart failure.
  Pun R, et al. J Am Heart Assoc 2025;14(6):e038286, PMID 40055864 — "BubR1 Controls Heart Development by
  Promoting Expression of Cardiogenesis Regulators."
- *Microcephaly.* Simmons AJ, et al. Hum Mol Genet 2019;28(11):1822–1836. PMID 30668728 — conditional
  cortical BubR1 loss causes microcephaly through shortened mitosis and massive apoptosis.
  Sterling NA, et al. Front Cell Dev Biol 2023;11:1282182. PMID 37900274 — a substantial fraction
  of that cell death is p53-independent, and co-deleting *Trp53* gives only minimal rescue of
  cortical size.

**Downstream mediator, not BubR1 itself.** Sieben 2020's own mechanistic finding is that the
severity difference between `H/L1002P` and `H/H` skeletal muscle tracks SASP complexity, not
senescence burden:
> "rates of senescence are similar between BubR1^H/H and BubR1^H/L1002P mice based on levels of
> Cdkn2a expression… Because BubR1^H/H muscles seem to have a more complex secretome than their
> BubR1^H/L1002P counterparts, one plausible explanation would be that senescent cells that
> accumulate in BubR1^H/H muscles have more profound degenerative properties."

That points at senolytic/SASP-directed intervention as an orthogonal lever, independent of BubR1 level.

### 4c. Human genotype–phenotype correlation in BUB1B-MVA — what exists

Thin, and the field says so. Sieben 2020:
> "Although MVA cases with BUBR1 mutations constitute the largest cohort, the overall size of this
> group is still too small for meaningful genotype-phenotype correlation analyses."

What I found:

- **García-Castillo H, Vásquez-Velásquez AI, Rivera H, Barros-Núñez P. Am J Med Genet A
  2008;146A(13):1687–1695. PMID 18548531.** The only clinical-subtype analysis I located. Reviewed
  21/35 MVA patients screened for *BUB1B*. Reported three pictures: patients with **monoallelic**
  *BUB1B* mutations "were severely affected with Dandy-Walker complex (7/8), cataracts (6/6), and
  Wilms' tumor (7/8)"; patients with **biallelic** mutations "showed a moderate phenotype";
  patients **without** *BUB1B* mutations "were mildly affected." Important caveat: the
  "monoallelic" group is now understood to carry a second, non-coding hypomorphic allele
  (Ochiai 2014, PMID 24344301), so the mono/bi contrast is not a clean protein-level contrast. It
  does not settle the dosage question either way, but it is the one human severity ranking on record.
- **Pavone P, et al. Neurol Sci 2022;43:6763–6775. PMID 35804254.** Reviews 31 MVA1 cases plus 4
  prenatal, claims "a new MVA1 genotype-phenotype correlation" around a novel `p.Arg893Ser`
  variant. I read the abstract only; the abstract does not state a quantitative correlation, and
  I did not obtain the full text. **Unresolved.**
- **Hanks S, et al. Nat Genet 2004;36:1159–1161. PMID 15475955** — the founding paper; five
  families. **Hanks S, et al. Cancer Lett 2006;239(2):234–238. PMID 16182441** — CGH and *BUB1B*
  analysis in MVA-associated childhood cancers. **Rio Frio T, et al. N Engl J Med
  2010;363:2628–2637. PMID 21190457** — a homozygous intronic `c.2386-11A>G` splice variant
  presenting as adult-onset GI neoplasia rather than childhood MVA, which widens the phenotypic
  range attributable to *BUB1B* alone.
- **Malumbres M, Villarroya-Beltri C. Nat Rev Genet 2024;25(12):864–878. PMID 39169218** — current
  review of MVA in development, ageing and cancer. Paywalled; I read the abstract only and did not
  extract any dosage or genotype–phenotype claim from it. **Unresolved.**

---

## 5. What this implies for the pharmacodynamic marker

**Do not use total BubR1 protein as the sole PD endpoint.** Sieben 2020 makes it indefensible as
a severity predictor: `H/L1002P` and `H/H` are indistinguishable on it and differ 1.75-fold in
lifespan; `+/L1002P` and `+/X753` are indistinguishable on it and differ completely in mTORC1
signalling and muscle transcriptome.

**Do keep total BubR1 as the primary target-engagement readout.** It is the variable that
Suijkerbuijk 2010 raised to rescue L1012P function, the variable Baker 2013 raised to correct
`BubR1^H/H` progeria in vivo, the variable Matsuura 2006 restored by chromosome transfer in
patient cells, and the variable that Bohers 2008 and Baker 2004 show has hard functional
thresholds (~50% for aneuploidy onset in HeLa; ~10% for mouse viability). Nothing displaces it as
the thing an intervention must move.

A defensible marker stack, in order:

1. **Mutant-allele-specific BubR1 abundance** (primary target engagement). Total protein cannot
   distinguish "more wild-type" from "more stabilised mutant," and for the stabilisation
   mechanism we are proposing (HSP90/proteasome axis, per Suijkerbuijk 2010 Fig. 3C) the mutant
   pool is the pool being acted on. Requires an allele-discriminating assay — targeted MS/PRM on
   the tryptic peptide spanning the substitution is the obvious route. I did not find any paper
   that has done this for a BUB1B MVA allele.
2. **Kinetochore-localised BubR1, normalised to centromere signal** (mechanistic proximal
   readout). This is the assay Sieben 2020 uses (Supp. Figs 2F, 10A–B), Matsuura 2006 uses, and
   Suijkerbuijk 2010 uses. In the biallelic genotypes it does move with total abundance, so it is
   a check on whether the added protein is going where it needs to go rather than an independent
   variable.
3. **SAC function** — colcemid/nocodazole challenge with a mitotic-index readout. This is the
   assay in which raising L1012P produced full rescue (Suijkerbuijk 2010 Fig. 4D). It is the
   closest existing thing to a validated response endpoint for this allele.
4. **PP2A-B56 recruitment / KARD phosphorylation** — the arm the pseudokinase domain controls
   (Gama Braga 2020, PMID 33207204). Add this because it is the specific way an
   L1012P-in-the-pseudokinase-domain protein could be defective at normal abundance and pass
   assay 3 anyway. This is the single highest-value experiment to de-risk the thesis.
5. **PCS rate**, separately from aneuploidy rate. Sieben 2020 shows PCS diverges dramatically
   between `H/L1002P` and `H/H` at identical aneuploidy, which makes it the readout that is
   sensitive to the mutant protein's interference effect.
6. **Aneuploidy rate**: keep as a safety/mechanism check, drop as a severity surrogate. The
   authors' own data anticorrelate it with severity across the `H/L1002P` vs `H/H` pair, and they
   warn that cultured-MEF aneuploidy is unreliable in vivo. Tissue FISH is what they trust.
7. **Downstream:** *Cdkn2a* expression and SASP composition, and for muscle, mTORC1 substrate
   phosphorylation (p70 S6K T389, 4EBP1 T37/46), per Sieben 2020 Figs 5C and 10.

**How the thesis should be worded.** Not "MVA1 is a dosage disease." Rather:

> For the kinase/pseudokinase-domain missense class of BUB1B alleles — the largest class, and the
> one containing L1012P — low BubR1 abundance is the established proximate cause of mitotic
> failure, is the only variable that has been shown to rescue the defect in patient cells and to
> correct progeria in an MVA mouse model, and is therefore the right thing to intervene on.
> Abundance sets the disease regime; allele identity modulates severity within it. Raising
> abundance above threshold is expected to move the phenotype toward the mild end of the regime,
> not to normal.

That claim survives everything in this document. "Total BubR1 predicts severity" does not.

---

## 6. Searches run

E-utilities `esearch` + `esummary`/`efetch` against PubMed, plus PMC full-text retrieval.
Full-text read in full: Sieben 2020 (PMC6934189) and its supplemental PDF
(jci.org/articles/view/126863/sd/1); Suijkerbuijk 2010 (PMC2887387); Baker 2013 (PMC3707109).

**Searches that returned relevant material**

| Query | Count | Notes |
| --- | --- | --- |
| `PMC6934189` (db=pmc) | 1 | Resolved to PMID 31738183. |
| `Jeganathan AND BubR1 AND allelic` | 3 | Confirmed the target paper is Sieben 2020, not "Jeganathan 2019". |
| `BUB1B AND mosaic variegated aneuploidy` | 37 | Core set. |
| `BUB1B AND (genotype-phenotype OR genotype phenotype)` | 9 | Yielded García-Castillo 2008, Pavone 2022. |
| `BUBR1 AND (residual OR hypomorph*) AND (severity OR phenotype)` | 6 | Yielded Bohers 2008. |
| `mosaic variegated aneuploidy AND (clinical OR patients OR cohort OR review)` | 101 | Yielded Malumbres 2024 review, clinical case series. |
| `BubR1 AND (non-mitotic OR interphase OR "insulin receptor" OR megakaryocyte OR meiosis)` | 107 | Yielded Park 2023. Mostly oocyte-meiosis noise. |
| `(premature chromatid separation OR PCS syndrome) AND BUB1B` | 20 | Yielded Ochiai 2014, Kato 2017, Miyamoto 2011, Matsuura 2006, Rio Frio 2010. |
| `BubR1 AND (overexpress* OR "increased expression") AND (lifespan OR aging OR healthspan OR aneuploidy)` | 32 | Yielded Baker 2013, Weaver 2016. |
| `BubR1 AND (acetylation OR PCAF OR SIRT2 OR ubiquitin OR degradation)` | 179 | Yielded FBXW7/BubR1, SIRT2/BubR1 reviews. |
| `(mitotic checkpoint complex) AND BUBR1 AND (stoichiometry OR quantitative OR concentration OR titration)` | 21 | Yielded Banerjee 2022, Ibrahim 2015. |
| `Choi E AND Yu H AND (insulin OR "metabolic homeostasis")` | 7 | Yielded Choi 2016 Cell. |
| `BUBR1 AND (pseudokinase OR "kinase activity" OR "kinase-dead")` | 55 | Yielded Gama Braga 2020, Tang 2020. |
| `Hanks S AND BUB1B` | 4 | Yielded Hanks 2004, Hanks 2006. |
| `Rahman N AND ("mosaic variegated aneuploidy" OR BUB1B)` | 8 | Confirmed the UK cohort papers. |
| `BUB1B AND ("mosaic variegated aneuploidy" OR "premature chromatid separation") AND 2020:2026[dp]` | 11 | Confirms nothing new on L1012P since Sieben 2020. |
| PubMed `citedin` link for PMID 31738183 | 12 | No follow-up study on the L1002P allele; no rescue attempt. |

**Searches that found nothing relevant — these are the basis for the "I did not find" statements below**

| Query | Count | Result |
| --- | --- | --- |
| `("mosaic variegated aneuploidy" OR "PCS syndrome" OR MVA1) AND ("clinical severity" OR "phenotypic severity" OR "disease severity")` | **0** | Nothing. |
| `("mosaic variegated aneuploidy" OR "PCS syndrome") AND ("BUBR1 protein" OR "BubR1 level*" OR expression) AND (correlat* OR sever*)` | 6 | All already known; none correlates measured patient BUBR1 protein with clinical severity. |
| `mosaic variegated aneuploidy AND (treatment OR therapy OR therapeutic)` | 23 | Only supportive/oncologic care — chemo, nephrectomy, HSCT. No molecular intervention. |
| `(BUB1B OR BubR1) AND ("gene therapy" OR AAV OR "adeno-associated" OR antisense OR ASO OR "mRNA therapy")` | 18 | All cancer-biomarker papers. No MVA1 therapeutic. |
| `(BUBR1 OR BUB1B) AND (chaperone OR "proteasome inhibitor" OR stabiliz* OR "small molecule" OR therapeutic OR rescue)` | 340 | Overwhelmingly cancer bioinformatics. No BubR1-stabilisation therapeutic programme. |
| `BUBR1 AND (HSP90 OR geldanamycin OR chaperone OR proteostasis)` | 5 | None follows up Suijkerbuijk 2010's HSP90 finding for therapeutic stabilisation. |
| `("mitotic checkpoint complex" OR "APC/C") AND (MVA OR "mosaic variegated aneuploidy" OR "patient cells") AND BUBR1` | 1 | A CDC20 progeroid paper, unrelated. |
| `("mosaic variegated aneuploidy" OR "premature chromatid separation") AND ("natural history" OR registry OR series OR "case series")` | 3 | No MVA1 natural-history registry. |
| `("aneuploidy rate" OR "aneuploidy burden" OR "percentage of aneuploid cells") AND (phenotype OR severity) AND (mosaic OR MVA)` | 8 | Nothing on MVA1 severity. |
| `(BUB1B OR BUBR1) AND (iPSC OR "induced pluripotent" OR organoid)` | 2 | Neither is an MVA1 patient-derived model. |
| `BubR1 AND mouse AND aneuploidy AND (uncoupl* OR dissociat* OR "does not correlate" OR independent)` | 7 | Nothing beyond Sieben 2020 on level/severity uncoupling. |

**Statements of the form "I did not find" — with the searches above as the record:**

- I did not find a study quantifying residual BUBR1 protein across a series of human MVA1 patients
  and correlating it with clinical severity.
- I did not find any small-molecule, chaperone, proteostasis, gene-therapy or ASO programme aimed
  at raising BubR1 in MVA1.
- I did not find any experiment raising BubR1 in a *BubR1* MVA-missense mouse (`H/L1002P`,
  `H/X753`, `+/X753`). The only in vivo rescue on record is the `BubR1^H/H` pure hypomorph
  (Baker 2013, PMID 23242215).
- I did not find any test of whether L1012P/L1002P protein restored to normal abundance is
  competent in PP2A-B56 recruitment / checkpoint silencing (the pseudokinase-domain function
  described by Gama Braga 2020, PMID 33207204).
- I did not find an allele-resolved measurement of mutant vs wild-type BubR1 abundance in any
  MVA1 sample.
- I did not find a follow-up study of the *BubR1^L1002P* allele since Sieben 2020.

None of these is a claim that the work does not exist. They are the limits of the searches listed.

**Marked unresolved (paywalled or not read in full):**

- Pavone P, et al. Neurol Sci 2022, PMID 35804254 — abstract asserts a genotype-phenotype
  correlation; abstract does not state what it is. Full text not obtained.
- Malumbres M, Villarroya-Beltri C. Nat Rev Genet 2024, PMID 39169218 — paywalled; abstract only.
- Kyuragi R, et al. ATVB 2015, PMID 25524773 — read only as cited in Sieben 2020; the
  `BubR1^L/L` protein level was not independently verified.
- Miyamoto T, Matsuura S. Oncotarget 2015, PMID 26309087 — no abstract in PubMed; full text not read.
