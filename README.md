# mva-research

A drug search for a child whose cell-division checkpoint is already failing will rank
checkpoint *inhibitors* at the top — oncology compounds built to kill chromosomally unstable
cells, which is the opposite of what this child needs. This repo is the pipeline that catches
that, by making every drug carry a signed direction of effect and refusing to build a
hypothesis when the signs disagree.

Built for **Rare Disease, Real Kid: The MVA Hackathon 2026** (Sage Bionetworks), Track 2.

- The scientific case, with every citation: **[`docs/track2-hypothesis.md`](docs/track2-hypothesis.md)**
- Pitch video, 3:45: **<https://youtu.be/lG6TIhnSnPU>**

*No patient data is in this repository. The config layer refuses a workspace inside the repo
tree or under a cloud-synced directory, and everything demoable runs on a synthetic case with
fictional genes. Nothing here is medical advice.* ([how that's enforced](docs/privacy-model.md))

## What the screen found

**13 agents in. One kept. Twelve rejected.**

Seven of those twelve were not rejected for measured harm. They were rejected because nobody
has ever asked whether they make chromosomal instability worse, and in this disease an unasked
question is not a clean answer.

The survivor, nicotinic acid, is carried as an experimental probe. Which way it moves BubR1
abundance has never been measured, the report says so out loud, and no dose is proposed
anywhere.

## How a drug gets rejected

**Direction of effect is signed and mandatory.** Every drug records which way the target has
to move and which way the published evidence says that drug moves it. When those disagree the
hypothesis object cannot be constructed at all — a validator refuses it, so no amount of
re-weighting brings a contraindicated compound back (`GP-16`).

**Unmeasured is its own state.** If either direction was never measured the value is
`unknown`, which is neither agreement nor disagreement and costs the candidate (`GP-14`).

**A protective response can't be a target.** If the body's own reaction to the disease turns
out to be protective, it cannot be nominated as something to reverse. That rule rejected a
senolytic backed by the best animal evidence in the entire document.

## My first ranking was wrong

A safety check had been marked clear on the strength of a study that never counted
chromosomes. Review caught it. Correcting **two cells in an evidence table** reversed the
ranking, with no pipeline code changed.

I'd rather you heard that from me, because it is the design working. Weights live in
configuration and evidence lives in tables, so you can see which row moved the decision and
move it yourself.

Afterwards I re-read all 77 claims in the report against the actual papers instead of against
the identifiers: **47 held, 22 needed weakening, 3 I withdrew, 5 I could not verify either
way.**

## Run it

```bash
just bootstrap   # Python 3.12, deps, pre-commit privacy hook
just verify      # lint, types, 1,782 tests, architecture, docs, privacy audit
just demo        # end-to-end synthetic run
```

Needs macOS or Linux, [`uv`](https://docs.astral.sh/uv/) and `just`. No genomics tooling
required to bootstrap.

The demo fixture is built adversarially rather than as a happy path: an **in-cis** pair with
phase-set evidence that must not be reported as a compound heterozygote, a multiallelic site
needing correct per-allele depth, a low-quality call that has to be flagged instead of
deleted, and a drug catalogue seeded with a right-target/wrong-direction agent. If the
pipeline gets this fixture right it's right for structural reasons, not because it was tuned
to a leaderboard. `just demo-determinism` runs it twice and compares artifact hashes.

## What this does and does not establish

Reproducible computation is not biological evidence, so, plainly:

- The BUB1B mechanism chain and 13-agent catalogue in `knowledge/literature/bub1b/` are real,
  every row cited to a PMID and checked in both directions by test. Real biology, still a
  hypothesis.
- Every scoring weight is heuristic, calibrated against no labelled dataset. No composite
  score is a probability.
- Annotation and the catalogue the *demo* runs on are synthetic stand-ins for
  VEP/gnomAD/ClinVar/DrugBank. `docs/maturity-ledger.md` grades every module `real`,
  `synthetic-substitute` or `stub`.
- There is no copy-number or structural variant calling, which in a chromosomal-instability
  disorder is a material gap rather than a rounding error (TD-03).
- No human domain expert has reviewed any of this. Adversarial model review is a filter, not
  a substitute, and treating model fluency as evidence quality is the exact failure this
  project exists to resist (TD-09).

## Further reading

| | |
|---|---|
| The rules, with IDs | `docs/golden-principles.md` |
| Why a choice was made | `docs/decisions/` (28 records) |
| Architecture and layering | `docs/architecture.md` |
| Threat model, and deletion on macOS | `docs/privacy-model.md` |
| What's real and what's mocked | `docs/maturity-ledger.md` |
| Deferred work, with costs | `docs/tech-debt.md` |

Track 1 scoring methodology is adapted from the CAGI6 Rare Genomes Project assessment
(Stenton et al., 2024). HPO identifiers in the synthetic fixture are real public ontology
terms; the phenotype profile and all gene symbols are fabricated.
