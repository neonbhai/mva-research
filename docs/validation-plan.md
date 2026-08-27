# Validation plan

How we would establish that this pipeline is right — and what we have and have
not yet done. The distinction between "the code works" and "the science is
correct" is maintained throughout; the first is testable today, the second is
not.

## Level 1 — Software correctness (DONE)

Unit, integration and golden tests over the synthetic case. The list of required
behaviours is in `docs/current-status.md`. Key locks:

- the synthetic causal pair ranks first, produced by general scoring rather than
  a special case;
- the wrong-direction drug is rejected;
- repeat runs are byte-identical;
- the submission CSV round-trips against the verified challenge contract.

**What this establishes:** the implementation does what it says. **What it does
not establish:** that what it says is biologically right.

## Level 2 — Determinism and reproducibility (DONE)

`just demo-determinism` runs the demo twice into separate workspaces and compares
artifact hashes. The run manifest records config hash, input hashes, git commit,
dirty flag, tool versions and the network profile, so a third party can tell
whether a rerun is expected to match.

**Gap:** reproducibility is verified on one machine and one OS. Cross-platform
byte-identity (linux/amd64 vs macOS/arm64) is untested; the container exists for
this and has not been exercised.

## Level 3 — Known-answer validation on public cases (NOT DONE)

The honest next step, and the one that would most change our confidence.

Take published, solved rare-disease cases with public variant data and known
causal genes, run the pipeline blind, and measure where the true answer ranks.
Candidate sources: ClinVar pathogenic compound-het pairs embedded in synthetic
backgrounds; the CAGI Rare Genomes Project public materials (the challenge's own
scoring methodology derives from CAGI6, Stenton et al. 2024); DDD/DECIPHER
published cases.

Metrics: top-1 and top-10 recall, rank distribution of the true pair, and — more
informative than either — a per-component breakdown of *why* the true answer
ranked where it did.

**Until this is done, no claim about real-world sensitivity is supportable.** The
synthetic case proves the machinery, not the calibration.

## Level 4 — Weight calibration (NOT DONE, deliberately)

Current weights are heuristics (ASSUMPTION-SCORING-01). Calibrating them requires
a labelled set from Level 3. Constraints we are imposing on ourselves:

- no tuning against the live leaderboard (ASSUMPTION-SCORING-03) — that fits the
  answer key, not the biology, and invalidates the golden tests as a measure of
  anything;
- any weight change needs a decision record, a test, and a before/after
  comparison (GP-32);
- golden expectations are never re-baselined to accommodate a new weighting.

## Level 5 — Scientific review (PARTIALLY DONE)

Adversarial review of genomic assumptions and pharmacological reasoning by
reviewers independent of the implementation. Findings are promoted into tests or
lints rather than left as prose, so the gate keeps enforcing them.

**Gap:** no human domain expert has reviewed this. Model-based adversarial review
is a filter, not a substitute — and treating model fluency as evidence quality is
one of the failure modes this project exists to resist.

## Level 6 — Wet-lab validation (OUT OF SCOPE, specified)

For the real case, every top candidate and every drug hypothesis carries the
experiment that would test it. The pipeline's job is to produce a falsifiable,
prioritised list of experiments — not conclusions.

Track 1 candidates: parental segregation testing (resolves the blocking phase
question), RNA studies for splice predictions, functional assay of the gene
product.

Track 2 hypotheses: the `proposed_validation_experiment` field is mandatory on
every `DrugHypothesis`. For the lead candidate class, the shape is a
dose-response missegregation/micronucleus assay in patient-derived cells across a
clinically achievable concentration range, with an isogenic corrected control and
a blinded scorer.

## What would falsify our top candidate

Stated up front, because a hypothesis that cannot be falsified is not a
hypothesis:

- parental testing showing the two variants in **cis**;
- a functional assay showing preserved gene-product activity;
- the same genotype in an unaffected family member;
- a better-fitting candidate in a gene whose phenotype match explains features
  ours does not.

## Standing limitations

- Annotation is a synthetic substitute, not VEP/ClinVar/gnomAD. Consequence and
  frequency values in the demo are fabricated (`docs/maturity-ledger.md`).
- No structural or copy-number variant calling. In a chromosomal-instability
  disorder this is a material gap, not a rounding error.
- Proband-only: no segregation, no de-novo detection, no trio-based phasing.
- Phenotype matching is a flat term-overlap score with no ontology-graph
  propagation — an observed child term does not currently credit an associated
  parent term.
