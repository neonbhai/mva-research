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
pair, which makes biallelic recessive the right prior. The architecture still
represents homozygous-recessive, de-novo dominant, X-linked, mitochondrial and
mosaic models, and single-variant candidates share the same ranked list. A prior
that becomes a filter stops being a prior.

### ASSUMPTION-MOSAIC-01 — Skewed allele balance may be signal
In a mosaic aneuploidy disorder, a heterozygous call with allele balance below
the usual band may reflect genuine somatic mosaicism rather than an artifact. We
flag `possible_mosaic` in preference to `low_allele_balance` above a configured
floor, so mosaic candidates are not down-ranked as noise. This is a deliberate
sensitivity/specificity trade specific to this disease context.

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

### ASSUMPTION-DRUG-03 — Repurposing requires approval
An investigational agent or tool compound may be a valid research probe but is
not a repurposing candidate. Approval status is tracked separately from
direction, and the two rejection reasons are never conflated.

### ASSUMPTION-DRUG-04 — Symptom management is not mechanism correction
`InterventionClass` separates `DISEASE_MODIFYING` from `SYMPTOMATIC`,
`SURVEILLANCE`, `SUPPORTIVE` and `PREVENTIVE`. Presenting an anticonvulsant as
addressing a chromosome-segregation defect is a category error the type system
prevents.

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
context, not a pass.

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

### ASSUMPTION-SCORING-03 — No leaderboard tuning
Weights are not tuned against the live leaderboard during scaffolding. Doing so
would fit the answer key rather than the biology, and would invalidate the golden
tests as a measure of anything.
