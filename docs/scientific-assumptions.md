# Scientific assumptions

Every assumption this pipeline makes, stated so it can be argued with. Each has
an ID cited from code and tests.

**Nothing in this repository is medical advice.** Track 2 output is a research
hypothesis requiring pre-clinical validation.

---

## Inheritance and phase

### ASSUMPTION-PHASE-01 — Two heterozygous variants in one gene are NOT assumed to be in trans
A compound heterozygote requires the two variants on **opposite haplotypes**. In
cis, one allele remains intact and the recessive mechanism does not apply.
Proband-only short-read data usually cannot resolve this.

We therefore preserve `PhaseStatus.UNKNOWN` end to end and apply a configured
penalty (`phase_weights.unknown`, default 0.55). We never upgrade to trans. Every
pair with unknown phase carries a **blocking** open question naming the resolving
test: parental segregation, read-backed phasing, or long reads.

In-cis pairs are down-weighted to near-zero (0.10 / 0.02) rather than deleted,
because short-read phasing calls can be wrong (GP-13).

### ASSUMPTION-INHERITANCE-01 — Compound-het is a prior, not a constraint
The challenge's ground truth is a clinically validated compound-heterozygous
pair, which makes biallelic recessive the right prior. Single-variant candidates
share the same ranked list rather than being filtered out. A prior that becomes
a filter stops being a prior.

`InheritanceModel` is the **vocabulary of the domain, not a list of things this
pipeline infers**, and the difference is worth stating plainly. From a
proband-only VCF, `mva.prioritization.pairing` actually produces six of its nine
members: `COMPOUND_HETEROZYGOUS`, `HOMOZYGOUS_RECESSIVE`, `X_LINKED_RECESSIVE`,
`MITOCHONDRIAL` (any `chrM` call — mtDNA dosage is per-cell heteroplasmy, not
two gene copies), `MOSAIC` (a call carrying `possible_mosaic`) and `UNKNOWN`.

The other three are **unreachable by construction**, each because it needs data
this pipeline never receives: `DE_NOVO_DOMINANT` requires trio genotypes;
`AUTOSOMAL_DOMINANT` and `X_LINKED_DOMINANT` require segregation in affected
relatives or de-novo status, without which a lone heterozygote is scored
`UNKNOWN` rather than promoted. The exclusion list lives in code as
`UNPRODUCED_INHERITANCE_MODELS`, with the reason attached to each member, and a
test asserts that every enum member is either produced by the pipeline or on
that list — so this paragraph cannot quietly become false.

### ASSUMPTION-MOSAIC-01 — Skewed allele balance may be signal
In a mosaic aneuploidy disorder, a heterozygous call with allele balance below
the usual band may reflect genuine somatic mosaicism rather than an artifact. We
flag `possible_mosaic` in preference to `low_allele_balance` above a configured
floor, so mosaic candidates are not down-ranked as noise. This is a deliberate
sensitivity/specificity trade specific to this disease context.

**The mosaic window is one-sided, and the band order is load-bearing.**
Mosaicism dilutes the alternate allele across a fraction of cells; it can never
concentrate it. An allele fraction *above* `max_allele_balance_het` is therefore
allelic imbalance, reference-allele dropout, or a homozygote mis-called as a
heterozygote — never mosaicism. Both the QC stage and the scoring stage test the
above-band case **first** and give it its own multiplier and its own note. A
scoring copy that tested "outside the band" before "above the ceiling" reported
an allele fraction of 0.95 as "within the mosaic window, so treated as possible
mosaicism rather than noise", and that sentence reached the dossier.

### ASSUMPTION-MOSAIC-02 — The het band is applied to the site allele fraction
The quantity banded is `VariantRecord.allele_fraction`, not
`Genotype.allele_balance`. For a biallelic call the two are identical. For a
record decomposed from a multiallelic site they are not: after splitting
`A>G,GT` with `GT=1/2` and `AD=2,21,21`, each split record holds the *site's*
REF depth (2) and its own ALT depth (21), so `alt/(ref+alt)` is 0.91 and a
textbook compound heterozygote looks homozygous. `alt/DP` is 0.48, which is what
the band is for. Where the site depth is unknown the fraction is `None` and no
allele-balance finding may be raised at all, because a confident wrong flag is
worse than no flag (GP-14).

This is enforced structurally: a lint in `tests/unit/test_architecture.py`
forbids reading `genotype.allele_balance` outside `src/mva/models/` and
`src/mva/ingestion/qc.py`, because the original fix was written once in the QC
stage and then undone by two prioritisation stages that re-derived the wrong
quantity from the raw field.

---

## Variant interpretation

### ASSUMPTION-TRANSCRIPT-01 — Canonical transcripts are not sufficient
A variant can be benign on MANE-Select and splice-disrupting on a
tissue-relevant isoform. We store **all** transcript annotations on the variant
and take the maximum severity across transcripts for a gene, rather than reading
the canonical row. Canonical/MANE status is retained as metadata, not used as a
filter.

### ASSUMPTION-FREQUENCY-01 — Absence of frequency data is not rarity
A variant missing from the frequency reference is scored at a mid-range default
(0.5), not 1.0. Absence usually means poor coverage in reference cohorts, not
that the variant is vanishingly rare. Recorded as a `no_frequency_data` flag and
an open question.

### ASSUMPTION-FREQUENCY-02 — Maximum population AF, not global AF
Rarity uses the highest AF across recorded populations. Global AF systematically
under-estimates frequency for variants enriched in ancestries
under-represented in reference cohorts, which would make common variants in those
populations look like plausible ultra-rare causes. Every frequency carries
source, version and population (GP-18).

**A population must be large enough to set that maximum.** Taken raw, "maximum"
hands the decision to whichever cohort was sampled least: one allele observed
once in 40 chromosomes is an AF of 0.025, above any plausible recessive
cut-point, while the same site sits at 8×10⁻⁶ across 125,000 global chromosomes.
The maximum then reports a genuine founder allele as `common_variant`, sinks its
rarity component, and drops it from the plausible-candidate set — destroying
exactly the case the maximum-AF rule exists to protect. Populations reporting
fewer than `frequency.min_allele_number` alleles (default 2000) are excluded
from the maximum and **recorded** in the rarity rationale rather than silently
dropped; a population that reports no allele number at all stays eligible,
because an unrecorded cohort size is unknown, not small (GP-14). gnomAD's
`grpmax` and the ACMG BA1/BS1 filtering allele frequency address the same
problem. Rationale and the choice of 2000:
`docs/decisions/0010-filtering-allele-frequency.md`.

### ASSUMPTION-PREDICTION-01 — In-silico prediction is never proof of causality
CADD, REVEL and SpliceAI scores are `COMPUTATIONAL_PREDICTION` tier and
`IN_SILICO_PREDICTION` type. A type-level validator forbids labelling them as
observed data. Reports state the tier next to the claim.

---

## Phenotype

### ASSUMPTION-PHENOTYPE-01 — Four-valued logic, never a boolean
`OBSERVED`, `EXCLUDED`, `UNCERTAIN`, `NOT_ASSESSED` are distinct.

- Only `EXCLUDED` contributes negative evidence — a clinician looked and did not
  find it.
- `NOT_ASSESSED` contributes **nothing in either direction** and is excluded from
  the scoring denominator. Including it would penalise a gene simply for being
  associated with features nobody checked, which manufactures evidence against
  the correct gene.
- `UNCERTAIN` is likewise uninformative for scoring, but reported separately.

A term absent from the record is `NOT_ASSESSED`, never `EXCLUDED`.

### ASSUMPTION-PHENOTYPE-02 — Extraction confidence ≠ clinical certainty
`extraction_confidence` describes how sure we are the term was correctly pulled
from the source document. It is not the clinician's confidence in the finding.
The two are not multiplied.

### ASSUMPTION-PHENOTYPE-03 — Entailment is directional
`OBSERVED` propagates **up** the ontology: observing a child term entails every
ancestor. `EXCLUDED` propagates **down**: excluding a parent excludes its
descendants. `UNCERTAIN` and `NOT_ASSESSED` propagate in neither direction.

The two directions are opposite and are not interchangeable. Propagating
`OBSERVED` downward would credit a gene for features nobody observed; propagating
`EXCLUDED` upward would exclude an entire organ system on the strength of one
absent sign. A record that entails a term both present and absent is **reported
as a conflict, never resolved by rule** — including when the two arrive under an
`alt_id` and its primary, which resolve to the same term.

### ASSUMPTION-PHENOTYPE-04 — Specificity comes from a corpus, not from depth
Information content is `IC(t) = -ln(n(t)/N)` over the real HPO annotation corpus
under the true-path rule, where `n(t)` counts diseases reaching `t` or any
descendant. It is **not** derived from graph depth, which is a known
methodological error: depth measures how finely an area of the ontology has been
subdivided, not how rare a finding is.

Two honest limitations. IC tracks **curation effort** as well as true rarity, so
a well-studied phenotype looks less specific than an equally rare neglected one.
And IC is **undefined — not zero** — for a term the corpus never reaches; such
terms are excluded from the mean and counted, never scored as uninformative.

### ASSUMPTION-PHENOTYPE-05 — The scored similarity is asymmetric on purpose
The scored coverage component is Köhler et al.'s `sim(Q→D)`: patient terms
matched against the gene's. The reverse direction is computed and reported in the
evidence payload but **is not scored on**, because averaging over a gene's
curated terms penalises well-studied genes for features nobody assessed in this
patient — which is GP-14 inverted. Both numbers are visible so a reader can audit
the choice.

---

## Mechanism

### ASSUMPTION-PHENOTYPE-06 — A family-history term is scored as the proband's own, deliberately
The referral document for this case lists **HP:0200067 Recurrent spontaneous
abortion** among the proband's clinical features. The proband did not have
miscarriages; the proband's parents did. Recording a parental finding in a
subject's phenotype profile is, strictly, a category error, and the loader has no
column that distinguishes "observed in the subject" from "observed in a
first-degree relative".

We keep the term as `observed` anyway, and say so here rather than hiding it,
because in a chromosome-instability disorder parental reproductive loss is
**evidence about the proband's genotype**: recurrent aneuploid conceptions are
the expected correlate of a segregation defect segregating in that family. The
referral document also directs participants to treat it as phenotypic input
rather than background.

The cost is real and bounded. HP:0200067 contributes information content to the
similarity score for genes annotated with it, so a gene causing isolated
recurrent pregnancy loss and nothing else would score higher than it should. The
mitigation is that it is one of eight terms, and the other seven are all proband
findings.

The honest fix is a `subject` column on the phenotype profile distinguishing
proband from relative, and a scorer that treats relative-observed terms as
evidence about the family's segregating genotype rather than the proband's
presentation. That is not implemented, and this assumption is the record of the
gap. See ASSUMPTION-PHENOTYPE-01 for why the status vocabulary is
four-valued, which is the same class of problem solved properly.

### ASSUMPTION-PHENOTYPE-07 — Association strength is gene-level validity, not a per-term claim
`association_strength` on a gene-phenotype association answers "how confident is
an expert panel that variation in this gene causes disease **at all**" — ClinGen
Gene-Disease Validity and DDG2P confidence, carried verbatim. It does **not**
answer "how strongly is this gene tied to *this particular feature*", and no
source in this pipeline answers that question.

Two consequences follow and are stated rather than hidden.

First, the value is constant across every term of a gene, so it cannot
discriminate *within* a gene. Because the score's specificity component is a
ratio of association weights, a per-gene constant cancels: the strength changes
what a reader sees in the evidence store, not the ranking. Within-gene
discrimination comes from the ontology's information content, and could come from
`hpo_frequency` — an obligate feature matching is stronger evidence than a very
rare one — but frequency is **not** an input to the score today. That is a
deliberate omission, not an oversight; wiring it in is a scoring change requiring
its own decision record and before/after.

Second, ClinGen and DDG2P curate the (gene, **disease**) pair while this column
reduces to the gene, taking the strongest classification across all of a gene's
diseases. A gene curated `definitive` for one disease and `limited` for another
reads `definitive` on the phenotypes annotated via the weaker one. The
disease-level join was measured, not assumed: it covers 27% of candidate rows
(HPO annotates against OMIM/ORPHA, ClinGen against MONDO, no crosswalk on disk).
See ADR 0021 and TD-20.

The related claim this replaces is the one the pipeline used to make by accident:
that HPO's phenotype **frequency** (`HP:0040283`, `12/45`, `50%`) was a curation
confidence. It was not, and the two now occupy separate columns.

### ASSUMPTION-MECHANISM-01 — Chain links are individually graded
A mechanism is a chain of typed links, each carrying its own tier, strength and
`is_directly_demonstrated` flag. `MechanismHypothesis.inferred_links` exposes the
weak points, and the relevance score is penalised for them. A chain is only as
strong as its weakest link, and a report must show which link that is.

### ASSUMPTION-MECHANISM-02 — Developmental window
Aneuploidy and its structural consequences are established during early
development. Post-natal modulation of a checkpoint cannot reverse malformation
already present at birth. A mechanistically correct agent can therefore still be
therapeutically irrelevant, and the mechanism record carries an explicit
`developmental_window_caveat` that reports must surface.

### ASSUMPTION-MECHANISM-03 — Mechanistic heterogeneity within a disease label
Genes grouped under one clinical label may act through different mechanisms
(checkpoint signalling, centrosome function, splicing, replication stress). A
drug rationale valid for one may be irrelevant or inverted for another. Mechanism
is resolved per gene, never per disease name.

### ASSUMPTION-MECHANISM-04 — Not every deviation from wild type is the disease
A node can deviate because it is broken or because it is *compensating*. Clearance
of aneuploid progenitors is the worked example: it deviates from wild type, and
pushing it back "corrects" a protective response. Nothing in a signed state field
distinguishes the two cases, so every `MechanismNode` carries a mandatory,
un-defaulted `deviation_is_pathological` flag. Where a node is compensatory, no
corrective direction is derived from its state — the direction check returns
`UNKNOWN` rather than guessing — and such a node may not be designated the
therapeutic target at all.

### ASSUMPTION-MECHANISM-05 — The direction triple is checked, not trusted
`disease_direction`, `required_correction` and the therapeutic target node's
`state_in_patient` are authored independently in two curated tables, while the
drug check compares an agent against `required_correction` alone. A single
mistyped cell therefore inverts the whole Track 2 gate silently. The three are
required to agree at construction time (`MechanismHypothesis`), so an inconsistent
chain cannot be built, let alone reported.

### ASSUMPTION-MECHANISM-06 — The target node is where the deviation is a movable quantity
The therapeutic target is chosen at the node whose deviation is a **quantity a drug
could plausibly change in the required direction**, not at the node that best *names*
the disease. The two are usually different, and picking the second is how a mechanism
write-up ends up measuring every candidate against a target whose entire accessible
pharmacology is contraindicated.

BUB1B-MVA is the worked example (ADR 0025). The disease is named after the spindle
assembly checkpoint, and every compound that binds the checkpoint is an *inhibitor*,
because the oncology programme wants the checkpoint weakened. BubR1 itself is a
pseudokinase, so there is no activity to agonise. Two rungs upstream, the deviation is
"there is not enough BubR1 protein" — graded, with a published threshold, demonstrated
to be correctable in patient cells, and upstream of **both** routes by which BubR1 loss
produces missegregation. That is the target.

Three tests a candidate node must pass: the deviation is quantitative rather than
categorical; the node is upstream of every branch the phenotype travels through; and
at least one real agent moves it in the corrective direction. A node failing the third
test is not disqualified as *biology*, but designating it the target converts the
direction gate from a filter into a guarantee of rejection, and the report then reads
as if nothing could ever help.

---

## Drug hypotheses

### ASSUMPTION-DRUG-01 — Direction of effect is the primary gate
This is the assumption the whole Track 2 design exists to enforce.

In a loss-of-function checkpoint disorder, a naive target-proximity or pathway
search surfaces **checkpoint inhibitors** at the top, because they bind exactly
the proteins named in the mechanism. Those compounds are developed to push
chromosomally unstable *tumour* cells past their aneuploidy-tolerance ceiling.
The patient's non-tumour cells already sit above that ceiling with no reserve. The
proximity is high and the sign is inverted.

Every drug therefore carries a signed `observed_direction` checked against the
mechanism's `required_correction`. Disagreement is disqualifying, enforced by a
model validator so no weight change can resurrect a contraindicated compound.

### ASSUMPTION-DRUG-02 — "Cannot determine" is not "agrees"
`directions_agree` is tri-state. `None` — when either direction is `UNKNOWN` or
`CONTEXT_DEPENDENT` — is a penalty and a flag, never agreement and never a
wrong-direction rejection. Some real agents genuinely have bidirectional,
dose- and context-dependent effects; forcing them to a single sign would be
fabrication.

A measured `NO_CHANGE` is **not** in that category. A demonstrated null is a
signed, established finding — usually the strongest evidence available *against*
the claim that this agent corrects this node — so it resolves to `False`, not to
`None`. Filing it as unsigned would award it the undetermined-direction credit and
a recommendation to go and measure a sign that has already been measured.

### ASSUMPTION-DRUG-03 — Repurposing requires approval
An investigational agent or tool compound may be a valid research probe but is
not a repurposing candidate. Approval status is tracked separately from
direction, and the two rejection reasons are never conflated.

### ASSUMPTION-DRUG-04 — Symptom management is not mechanism correction
`InterventionClass` separates `DISEASE_MODIFYING` from `SYMPTOMATIC`,
`SURVEILLANCE`, `SUPPORTIVE` and `PREVENTIVE`. Presenting an anticonvulsant as
addressing a chromosome-segregation defect is a category error the type system
prevents.

The class is a **curated label**, so the separation does not rest on it. The
statement that an agent "does not correct the mechanism" is emitted whenever its
`target_node_id` differs from the mechanism's `therapeutic_target_node_id`,
whatever the class cell says, and the non-fatal `MECHANISM_MISMATCH` reason is
carried on the hypothesis as a `concern` and printed beside the candidate rather
than discarded once the candidate survives triage.

### ASSUMPTION-DRUG-05 — Achievable concentration is required
A compound effective at 10 µM in culture that peaks at 0.1 µM in plasma is not a
therapy. Where either figure is unknown, `concentration_achievable` is `None` and
that gap is recorded as a concern, not resolved optimistically.

### ASSUMPTION-DRUG-06 — Paediatric tolerability does not transfer across populations
Tolerability in a paediatric *oncology* population is not evidence of safety in a
germline chromosomal-instability population. This caveat is a default field on
`PediatricEvidence`.

### ASSUMPTION-DRUG-07 — Could it worsen the disease?
Every candidate must answer whether it could increase aneuploidy or cancer
susceptibility. `None` (unassessed) is itself a blocking gap in this disease
context, not a pass — and it **blocks**: in a chromosomal-instability mechanism an
unassessed answer is disqualifying with `ONCOGENIC_RISK`, not merely recorded.
Outside such a context it stays a recorded, non-fatal concern. A report that calls
a gap blocking while the pipeline ranks the candidate anyway is telling the reader
two different things at once.

### ASSUMPTION-DRUG-08 — The tumour and the soma are two therapeutic contexts, not one
A child with a chromosomal-instability cancer-predisposition syndrome carries two
populations of cells with **opposite therapeutic requirements**, and an agent's
direction can be desirable in one and catastrophic in the other. This is not a
hypothetical: the compounds identified as *aneuploidy-selectively lethal* — autophagy
inhibition, HSP90 inhibition, AMPK-driven energy stress (Tang et al., *Cell* 2011,
PMID 21315436) — kill cells by exploiting exactly the property the patient's every
tissue already has. They are attractive against the rhabdomyosarcoma on precisely the
evidence that disqualifies them for the soma.

The pipeline scores a drug against **one** context: the soma. It is not a tumour-board
tool and must not be read as one. A rejection here therefore means "rejected for
systemic use in the affected child" and never "this compound is useless in this
disease". Where an agent is plausible for the tumour arm, the catalogue's
`validation_experiment` cell says so and names the study that would be required — a
therapeutic-index comparison against the patient's own non-tumour cells, which share
the aneuploidy that makes the tumour sensitive. That comparison is the whole question,
and it is not one a mechanism chain can answer.

Encoding the tumour arm as a second mechanism with its own target and its own signs is
the honest structural fix. It is not implemented, and this assumption is the record of
the gap.

### ASSUMPTION-DRUG-09 — A genotype-conditional hypothesis must declare its condition
Some repurposing hypotheses are valid only for a subset of alleles, and the subset is
usually not knowable from the disease label. Presenting such a hypothesis without its
condition is how a proposal that applies to a minority of patients gets read as
applying to the one in front of the reader.

BUB1B-MVA alleles fall into two published classes: those that lower BubR1 *abundance*,
and those that impose a *qualitative* checkpoint or attachment defect at normal
abundance (Suijkerbuijk et al., *Cancer Res* 2010, PMID 20516114). Every
abundance-raising hypothesis in this repository is conditional on the first class and
would be expected to fail against the second. Readthrough agents are conditional on a
nonsense allele *that still produces a transcript*, which the same paper reports
truncating BUB1B alleles largely do not.

The rule: state the condition in the hypothesis, name the assay that resolves it, and
put that assay **before** the efficacy experiment in the proposed order. This is also
what makes a mechanism-grounded hypothesis writable without patient data at all — the
condition is declared as a condition rather than silently assumed away, and the
resolving assay is run by whoever holds the genotype.

---

## Scoring

### ASSUMPTION-SCORING-01 — The weights are heuristics
Every weight in `config/default.yaml` was chosen by reasoning about a severe
paediatric recessive disorder. **None is calibrated against a labelled dataset.**
No composite score here is a probability, and none carries a claim of clinical
validity. Component scores stay visible so a reader can disagree with the
weighting rather than with an opaque number.

### ASSUMPTION-SCORING-02 — Contradictions subtract
The contradiction penalty is subtracted from the weighted sum, not multiplied
into it, and is weighted above the smallest positive components so that
contradicting evidence can actually move a candidate rather than being averaged
into irrelevance.

**`common_variant`, `low_quality_call` and `benign_consequence` are each counted
twice on purpose:** they lower a positive component (rarity, analytical validity,
molecular consequence respectively) *and* add a term to the subtracted penalty.
This is deliberate double counting, and it is stated here because an undocumented
double count is indistinguishable from a bug. The rationale is that the two
counts answer different questions. The component asks "how strong is this line of
support?" and degrades smoothly; the penalty asks "is there a recorded reason
this hypothesis is wrong?" and is emitted as a persisted `CONTRADICTS`
`EvidenceItem` a reviewer can argue with (GP-19). A common allele that scored
0.12 for rarity would otherwise still be outranked by an unannotated variant
scoring the neutral 0.5, and nothing in the output would say why the first one
should not be trusted. The cost is that these three findings move a candidate
further than a single-count model would; the penalty magnitudes were chosen with
that overlap in mind and are, like every number here, uncalibrated (GP-32).

### ASSUMPTION-SCORING-03 — No leaderboard tuning
Weights are not tuned against the live leaderboard during scaffolding. Doing so
would fit the answer key rather than the biology, and would invalidate the golden
tests as a measure of anything.
