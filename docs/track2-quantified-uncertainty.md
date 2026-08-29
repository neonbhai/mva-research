# Track 2 — Quantified uncertainty, without a fabricated probability

Companion to `docs/track2-hypothesis.md`. It answers the question a reader is entitled
to ask after §4 — *how likely is any of this to help a child?* — without answering it
with a number that would be invented.

---

## 1. The objection this document has to answer first

`docs/track2-hypothesis.md` §2.2 says:

> The mechanism relevance score is **0.706** — a heuristic weighting of curated strength
> labels, not a probability, and explicitly not a claim that the chain is *true*.

That is a good line, and it forbids something. If a heuristic weighting of our own
confidence labels cannot be read as a probability that the chain is true, then **a
"P(benefit) ≈ 0.15" attached to this hypothesis cannot be either.** It would be the same
object — our own beliefs, multiplied, and re-presented as evidence about the world — with
a more persuasive decimal point on it. A reviewer who reads §2.2 and then finds a posterior
one section later has found a document that criticises a move and then makes it.

So this document does three things and refuses a fourth:

| | |
|---|---|
| ✓ | A **class comparison** for this hypothesis, against outcome data that exists independently of us — a ranking against measured proxies, not a rate of its own |
| ✓ | An explicit **conjunction** of the conditions that must all hold, priced by the reader, not by us |
| ✓ | A statement of **which single measurement moves the number most** |
| ✗ | A posterior probability that *this* hypothesis benefits *this* child |

### Why a class base rate is legitimate where a hypothesis-specific posterior is not

The difference is where the inputs come from.

A **hypothesis-specific posterior** takes as input our own graded confidence in each link —
the same labels that produced 0.706. Its inputs are the beliefs under test. Computing it
would be scoring our own conviction and reporting the score as a finding. There is no
external observation anywhere in the calculation, so there is nothing that could make it
come out differently if we were wrong.

A **class base rate** takes as input somebody else's completed outcomes: how often did
drug-development programmes of a stated shape end in an approved, benefit-demonstrating
therapy? Those outcomes were recorded before this hypothesis existed and cannot be
influenced by it. The base rate says nothing about whether *our* chain is true. It says
what usually happens to hypotheses of this shape — which is a fact about the world, and a
fact this hypothesis has to argue its way out of rather than assume its way past.

Two honest limits on that legitimacy, stated now rather than in a footnote:

1. **Reference-class choice is a judgement.** Different defensible classes give different
   numbers, by a factor of several. §3 lists the classes considered and what each gives,
   so the choice is visible rather than buried.
2. **A base rate is not a forecast for one case.** It is the right thing to hold *before*
   evidence specific to this hypothesis exists — which is exactly the situation §6 of the
   hypothesis document is designed to end.
3. **A modelled figure is not an outcome record, and only outcome records earn the
   argument above.** The whole legitimacy of a class base rate rests on somebody having
   *counted* something that happened. A number produced by a model from assumed error
   rates, assumed power and an assumed prior does not qualify, however well published,
   because its inputs are beliefs and it moves when they do. §3 marks the one row where
   that applies, and §3 does not use it as an anchor.

---

## 2. There is no published base rate for exactly this class, and that is not a gap I can fill

The question I would most like to answer is: *of mechanism-grounded drug-repurposing
hypotheses published for rare monogenic disease, what fraction reach demonstrated clinical
benefit?*

**No such number exists, and it is not obtainable.** The denominator is unobservable.
Hypotheses are not registered; the ones nobody pursues are never counted; and the ones that
fail in a cell-culture experiment leave no record. Any figure of that form in circulation
is either a proxy from a later pipeline stage or made up.

Everything in §3 is therefore a **proxy conditioned on a later stage than we are at**, and
every one of them is consequently an over-estimate for a hypothesis that has not yet had
its first experiment. That is stated at each row rather than corrected for, because
correcting for it would mean inventing the correction.

---

## 3. The proxies, and what each one actually measures

All figures below were read from the cited sources, not recalled.

| Reference class | Figure | Source | Is it the right class for us? |
|---|---|---|---|
| Any drug programme, Phase 1 → approval, all indications | **13.8%** | Wong, Siah & Lo, *Biostatistics* 2019;20:273-286 (PMID 29394327); 406,038 trial entries, 21,143 compounds, 15,102 programmes, 2000–2015 | No — conditioned on reaching Phase 1. We are years short of that. |
| **Orphan drugs, all indications**, Phase 1 → approval | **6.2%** | Wong 2019, Table 4 | Closer, but dominated by oncology. |
| **Orphan drugs excluding oncology**, Phase 1 → approval | **13.6%** | Wong 2019, §4.3 | The closest *pipeline* analogue: a rare, non-cancer indication that has already reached first-in-human. |
| Candidate new molecular entities in early-stage trials → approval | **~10%** | King, Davis & Degner, *PLoS Genet* 2019;15:e1008489 (PMID 31830040) | Same conditioning problem. |
| **Uplift for a genetically supported target** | **>2×**, and for **Mendelian** associations the effect **holds prospectively** | King 2019 | **Yes — and this is our one strong claim to be above the class median.** *BUB1B*-MVA is a Mendelian association with an unambiguous causal gene and a coding mechanism. |
| Preclinical target–disease relationships that are true | **~7%** — and it is a **model output, not an observed rate** | Hingorani et al., *Sci Rep* 2019;9:18911 (PMID 31827124), a modelling study. The 7% (and the 92.6% false-discovery figure beside it) fall out of a probabilistic model whose inputs are *assumed*: a prior γ ≈ 0.005 that any given protein–disease pairing is causal, an assumed type-I error rate, assumed statistical power, an assumed number of causal genes per disease, and reported programme success rates. Move any of those and the 7% moves with it. | **Right stage, wrong kind of number.** Nobody followed a cohort of preclinical target–disease relationships to an outcome and counted. Usable as an order-of-magnitude sanity check on a prior; **not** usable as a base rate, and not used as one below. |
| Rare diseases with an FDA-approved drug after 40 years of the Orphan Drug Act | **~5%**; up to **15%** have at least one drug "developed and shown promise" | Fermaglich & Miller, *Orphanet J Rare Dis* 2023;18:163 (PMID 37353796); 6,340 designations across 1,079 diseases → 882 designations approved, covering 392 diseases | **Yes, as a population reality check.** MVA is not among the 392. |
| Animal-to-human concordance | Discordant in both directions | Perel et al., *BMJ* 2007;334:197 (PMID 17175568) | Directly relevant: tirilazad reduced infarct volume by 29% in animal models and **worsened** human outcome; corticosteroids benefited animal head-injury models and did not benefit patients. |

Wong 2019 adds the sentence that governs how all of this should be read:

> "The overall POS presented in this study … [is] much higher than the 1% to 3% that is
> colloquially seen **as it is conditioned on the drug development program entering
> Phase 1.**"

### What the proxies support — a ranking, not a rate

§2 said the denominator for this class does not exist. That forbids something this section
originally did anyway, which was to read a numeric range off the table as though the table
had produced it. It had not. So the conclusion is stated in the form the evidence actually
takes:

> **For the class this hypothesis belongs to — a mechanism-grounded, genetically supported,
> approved-agent repurposing proposal in an ultra-rare monogenic disease, stated before its
> first validating experiment — the probability of reaching demonstrated clinical benefit
> sits *below every measured proxy in the table above*, and none of those proxies measures
> it.**

That is a ranking, and the ranking is what the rows support. Every measured proxy is
conditioned on a stage this hypothesis has not reached: Wong's 13.6% starts at
first-in-human; King's ~10% starts at an early-stage trial; Fermaglich's 5% counts diseases
that had a statutory incentive and forty years of attempts behind them. Each is an upper
bound. The direction of every correction is downward, and nothing in the table says how far.

**If a number is wanted, it has to be labelled for what it is.** My judgement is that the
right order of magnitude is a few percent. That is a **judgement, not a derivation**: I got
there by taking Wong's 13.6% orphan non-oncology figure as a ceiling, moving down for the
several preclinical stages between here and Phase 1, and moving up for Mendelian genetic
support (King 2019). There is no arithmetic behind it, the three inputs are weighed by
feel, and I could not defend 3% against 1% or against 8%. It should be read as an
order-of-magnitude judgement and nothing more.

**Hingorani 2019's ~7% is specifically not available as the anchor for that**, and an
earlier version of this section used it as one. It is a model's output under an assumed
prior and assumed error rates, not a count of outcomes; anchoring on it would borrow the
credibility of an observation for the content of an assumption. Nor is it legitimate to
multiply the proxies together — Wong's 13.8% already contains programmes whose mechanism
turned out to be false, so multiplying by a mechanism-failure figure double-counts the same
failure.

**This is a statement about a reference class. It is not a probability that this hypothesis
is true, and no arithmetic in this document produces one.**

---

## 4. The conjunction: five conditions, all of which must hold

The proxies in §3 tell you where to start. The conjunction tells you *why* rates of this
kind come out where they do, and it is the part a reader can weigh against our evidence.

For a benefit to occur in a child with *BUB1B*-MVA, **every one** of these must be true.
Not most. All.

| # | Condition | What would have to be true | Current evidence state |
|---|---|---|---|
| **1** | **The mechanism is correct in this patient's cells** | BubR1 *abundance*, not a qualitative defect, is the quantity governing missegregation in this child's allele class | **Strong at the gene, unknown for the individual.** Suijkerbuijk 2010's own conclusion is that CIN in *BUB1B*-MVA "is a result of low BUBR1 protein abundance", and the eligible allele class is larger than §4.3 concedes (claims ledger W-2). But the class is unknown here by design — §6.3 is the assay that decides it. |
| **2** | **The direction at the molecule is correct** | An **approved** NAD⁺-raising agent raises BubR1 protein in a human *BUB1B*-mutant cell | **`cannot determine`, for every approved agent.** This is not "probably yes". Measured only for NMN — unapproved, in wild-type mice (claims ledger WK-51). Never measured in a human cell of any genotype. Nicotinamide's sign is ambiguous by construction (claims ledger W-1). |
| **3** | **Exposure is achievable at a tolerable chronic paediatric dose** | A dose that moves intracellular NAD⁺ enough to shift K668 acetylation, given for years, in a growing child | **Uncertain and under-examined.** The only in-vivo direction evidence is 500 mg/kg/day i.p. NMN for 7 days in mice. The only paediatric niacin series reports adverse effects in **76%** of 21 children and discontinuation in **38%** (Colletti et al., *Pediatrics* 1993;92:78-82, PMID 8516088). The largest niacin outcome trial found excess serious adverse events across five body systems (the HPS2-THRIVE trial; Landray et al., *N Engl J Med* 2014, PMID 25014686). Nicotinamide is the exception: ENDIT gave 1.2 g/m²/day to children for five years with no effect on growth (Gale et al., *Lancet* 2004, PMID 15043959) — and nicotinamide is the agent whose sign we cannot determine. |
| **4** | **The effect size is sufficient** | The achievable rise crosses whatever the real residual-BubR1 threshold is, by enough to change missegregation rate | **Unknown, and the threshold's location is unknown too.** The "~50%" figure is the upper edge of an unsampled gap in a HeLa knockdown series (claims ledger WK-11). And in the closest genetic experiment available — SIRT2 overexpression in BubR1^H/H mice — BubR1 was **not restored to wild-type levels** and the animals **remained small**; North 2014 attributes the shortfall to the hypomorph's deficit being transcriptional, which SIRT2 cannot correct because it can only stabilise protein that already exists (hypothesis document §4.2). That is the closest thing to a direct test of this condition, and it came out against it. |
| **5** | **The safety question is answerable** | Whether the agent raises aneuploidy or cancer risk can be established before chronic dosing in a germline chromosomal-instability syndrome | **Currently unmet for six of seven candidate agents.** Seven of twelve pipeline rejections are `oncogenic_risk`, and in every case the risk was never measured (ADR 0011). The one agent with an answer has it from a 12-month trial in adults with UV-driven skin cancer (Chen et al., *N Engl J Med* 2015, PMID 26488693), which §5.1 already says does not transfer. |

### Pricing it

I am not assigning numbers to these, because the honest number for conditions 2 and 4 is
*undetermined*, not *probable*, and putting a decimal on "undetermined" is the exact error
§1 forbids.

**What follows is a fact about multiplication, not a fact about this hypothesis**, and two
things have to stay in view while reading it. Both are fatal to treating any output of it
as a probability, and calling the exercise "illustrative" does not neutralise either — a
number that is quoted works like a number that was derived, whatever label sits beside it.

- **The inputs are arbitrary.** 0.9, 0.8, 0.7 and 0.5 are round numbers chosen to span a
  range. They are not estimates of conditions 1–5, nobody has grounds for such estimates,
  and that is precisely why the left-hand column is a hypothetical rather than a column of
  values. Two of the five conditions are `cannot determine` — not "likely", not "0.7".
- **The conditions are dependent, so the products are wrong even if the inputs were right.**
  A product of five terms assumes independence, and these five do not have it. The
  dependence runs against the hypothesis: if condition 1 fails (qualitative allele class),
  conditions 2 and 4 are moot; if condition 3 forces a low dose, condition 4 gets harder.
  Correcting for that is not something I can do either, because the correlations are as
  unmeasured as the terms.

| If you priced each of the five at… | The product of five such numbers would be |
|---|---|
| 0.9 | 0.59 |
| 0.8 | 0.33 |
| 0.7 | 0.17 |
| 0.5 | 0.031 |

**The single thing that table shows is that a conjunction of five conditions is small even
when no individual one looks unreasonable.** That is its whole content. None of its outputs
is a posterior, none should be quoted on its own, and none should be set beside §3 for
agreement — §3 is a ranking against proxies that measure something else, arrived at a
different way, and two numbers landing near each other by construction is not corroboration.
The move this section is built to support is not picking a row. It is forcing a reader who
wants a larger number to say **which condition they are raising, and on what evidence.**

---

## 5. What would raise the number, in order

The value of writing the conjunction down is that it identifies which measurement is worth
the most. These are not equal.

**1 — Does an achievable NAD⁺ exposure raise total BubR1 protein in a *BUB1B*-mutant human
cell?** (§6.1 readouts 1→3.)

This is the highest-value single measurement, for three reasons. It attacks **condition 2**,
which is the only one currently `cannot determine` rather than merely uncertain — the
pipeline's own accepted candidate carries `directions_agree = None` because of it. It is
cheap: patient-derived and isogenic lines, LC-MS/MS and a quantitative immunoblot. And it
is **symmetric** — a null result ends the hypothesis rather than deferring it, which is the
§4.3 headroom objection made concrete. Every other measurement is worth less because it is
conditional on this one.

**2 — Where is the residual-BubR1 threshold, actually?**

Attacks **condition 4**, which is currently priced on a number the source does not support
(claims ledger WK-11). Published knockdown data bound it somewhere between roughly 10% and
60% residual protein, and the two available HeLa series disagree. §6.1's isogenic
*BUB1B*-hypomorphic RPE1 allelic series measures this as a by-product of running the
experiment at all. If the threshold is near 15% rather than 50%, condition 4 becomes much
easier and the whole hypothesis improves; if it is near 50%, partial restoration is a much
harder ask.

**3 — Determine the allele class in the index patient.** (§6.3.)

Attacks **condition 1**. It does not change the *hypothesis*, but it decides whether the
hypothesis is about this child at all, and it is run by the team holding the sample. It is
first in time and second in information, which is why §6.3 has it first and this list has
it third.

**4 — Ask the aneuploidy question of nicotinic acid from existing trial data.**

Attacks **condition 5**, and §5.2 is right that it is unusually cheap. One correction from
the claims ledger (CL-63): the primary publication of the largest such trial —
HPS2-THRIVE, 25,673 adults, median 3.9 years — reports serious adverse events by body
system and **does not report cancer incidence**. So this is a **data request** against
existing trial safety capture, not a re-read of published tables. Still cheap; not free.
And a 3.9-year adult cardiovascular cohort is a weak instrument for embryonal-tumour risk
in a child with germline chromosomal instability, so a positive answer clears the gate
without answering the question.

**What would *not* move the number, and is worth naming because it is tempting:** more
mechanism. Another paper strengthening the NAD⁺–SIRT2–BubR1 chain raises confidence in
condition 1, which is already the best-supported of the five. The binding constraints are
2 and 5, and neither is a literature problem.

---

## 6. What this document does not license

- It does not license "≈15%", or any other single figure, as the probability that this
  hypothesis helps a child. §3 is a ranking against proxies that measure a later stage,
  plus one explicitly labelled order-of-magnitude judgement; §4 is arithmetic on arbitrary
  inputs over dependent conditions. Neither is a posterior about this chain, and combining
  them does not make one.
- It does not license optimism from the Mendelian genetic-support uplift. King 2019's
  >2× applies to *approval*, over a portfolio, for targets with human genetic evidence. It
  is a reason to place this hypothesis above the class median. It is not a reason to place
  it near 1.
- It does not license treating the class ranking as unchangeable. That is what §5 is for — and
  it is the one legitimate route from a class comparison to a hypothesis-specific number:
  **make a measurement.** A posterior computed from our own confidence labels is not
  evidence; a measurement of whether BubR1 rises in a mutant human cell is.

The reason the deliverable in §6 of the hypothesis document is an experiment rather than a
recommendation is precisely this. Right now the honest answer to "how likely is this to
help?" is the class ranking in §3 and nothing sharper, because nothing hypothesis-specific
has been measured yet. The
purpose of the experiment is to make that answer obsolete.

---

*Nothing in this document is medical advice. No patient data was accessed in producing it.*
