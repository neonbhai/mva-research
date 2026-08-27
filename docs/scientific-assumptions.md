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

---

## Mechanism

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
