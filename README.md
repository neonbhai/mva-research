# mva-research

A provenance-first drug-repurposing pipeline built for **Rare Disease, Real Kid: The
MVA Hackathon 2026** (Sage Bionetworks), Track 2.

## Run it

```bash
just all      # bootstrap, then the gate, then the synthetic demo
```

macOS or Linux, [`uv`](https://docs.astral.sh/uv/) and
[`just`](https://just.systems). That is the whole install list: no patient data, no
credentials, no encrypted volume, no genomics toolchain. `bootstrap` fetches
dependencies, so it needs the network once; no run touches it.

**Budget about six minutes.** Five of them are the test suite: 1,817 tests, a
progress bar, and nothing else to read. Kill it at four minutes and you file a
failure that never happened. Success is the last two lines — `verify: all gates
passed`, then a path to 33 artifacts under `$TMPDIR/mva-research-demo/runs/`, which
`just demo-artifacts` will list.

Shorter doors: `just demo` alone takes about 7 seconds and writes those same
artifacts; `just verify` alone is the gate.

## The failure mode this exists to catch

A drug search for a child whose cell-division checkpoint is already failing will rank
checkpoint *inhibitors* at the top — oncology compounds built to kill chromosomally unstable
cells, which is the opposite of what this child needs. This repo is the pipeline that catches
that, by making every drug carry a signed direction of effect and refusing to build a
hypothesis when the signs disagree.

- The scientific case, with every citation: **[`docs/track2-hypothesis.md`](docs/track2-hypothesis.md)**
- Pitch video, 3:45: **<https://youtu.be/lG6TIhnSnPU>**
- Architecture as a page: **<https://neonbhai.github.io/mva-research/architecture.html>**

*Everything demoable runs on a synthetic case with fictional genes, and no patient data is in
this repository. Nothing here is medical advice.*
([how that's enforced](docs/privacy-model.md))

## What the screen found

**13 agents in. One kept. Twelve rejected.**

Seven of those twelve were not rejected for measured harm. They were rejected because nobody
has ever asked whether they make chromosomal instability worse, and in this disease an unasked
question is not a clean answer.

The survivor, nicotinic acid, is carried as an experimental probe. Which way it moves BubR1
abundance has never been measured, the report says so out loud, and no dose is proposed
anywhere.

## How a drug gets rejected

**Direction of effect is signed and mandatory.** Every drug records which way the target must
move and which way the published evidence says the drug moves it. When those disagree the
hypothesis object cannot be constructed at all: a validator refuses it, so no amount of
re-weighting brings a contraindicated compound back (`GP-16`).

**Unmeasured is its own state.** If either direction was never measured the value is
`unknown`, which is neither agreement nor disagreement and costs the candidate (`GP-14`).

**A protective response can't be a target.** If the body's own reaction to the disease is
protective, it cannot be nominated as something to reverse. That rule rejected a senolytic
backed by the best animal evidence in the document.

## My first ranking was wrong

A safety check had been marked clear on the strength of a study that never counted
chromosomes. Review caught it. Correcting **two cells in an evidence table** reversed the
ranking, with no pipeline code changed.

That is the design working: weights live in configuration and evidence lives in tables, so you
can see which row moved the decision and move it yourself.

## Reproducing what's claimed

Four claims, four different answers to who can check them.

**1. The synthetic pipeline behaves as described.** Anyone: `just demo`, then
`just demo-determinism`, which runs it twice and compares artifact hashes. The fixture is
adversarial rather than a happy path — an **in-cis** pair with phase-set evidence that must
not be reported as a compound heterozygote, a multiallelic site needing correct per-allele
depth, a low-quality call that has to be flagged instead of deleted, and a drug catalogue
seeded with a right-target/wrong-direction agent. Passing it is structural, not leaderboard
tuning. The linux/amd64 container in
[`workflow/containers/README.md`](workflow/containers/README.md) reproduces the same 33
artifacts under the same run id.

**2. The 13-agent screen: twelve rejected, nicotinic acid kept at 0.500.** Anyone, with no
patient data: `just ranking`, `just drugs-to-avoid`, `just science-canary`. Both pages are
generated from a live run and locked against drift: `tests/unit/test_ranking.py` and
`tests/unit/test_drugs_to_avoid.py` fail when a committed page stops matching the generator,
and one of them cross-checks the live run against all thirteen scores typed by hand into
`docs/track2-hypothesis.md` §5. `just science-canary` is the inverse check: 12 controls
proving the direction gate still reaches "accept", since a stage that rejected everything
would pass every rejection test I have.

**3. The Track 1 variant ranking.** Nobody, without their own gated access from Sage; no copy
of that data is here and there will not be one. Attempt 1 banked 100.0 rank points and F-max
1.000 (`docs/repo-followup-backlog.md` §0) — the ceiling for both. Once a score cannot move,
it has stopped distinguishing anything, so I read it as a format check.

**4. The 77-claim audit: 47 held, 22 weakened, 3 withdrawn, 5 I could not verify.** Anyone
with library access, and there is no command. I re-read every claim against the papers rather
than against the identifiers, and `docs/track2-claims-ledger.md` records each one with its
source and grade. `just check-citations` proves only that the identifiers resolve against a
pinned ledger, a far weaker claim than "this paper supports this sentence".

## Running your own case

`mva case` is the front door; you do not need to read the schema first.

```bash
mva case init                   # guided; --no-input for scripts, --run to run it
mva case check -c mycase.yaml   # says whether it could run; runs nothing
mva case phenotype-template     # writes the phenotype TSV header to fill in
```

- **Phenotype is a TSV.** `hpo_id`, `label`, `status` required; `onset`, `provenance`,
  `extraction_confidence`, `notes` optional — the terms come out of a clinical document, so
  their provenance travels with them.
- **The VCF goes in as supplied**: `.vcf`, bgzipped or BCF, with `.tbi`/`.csi` used for
  per-contig streaming.
- **Raw reads are never read.** No BAM, CRAM or FASTQ path exists, so no alignment and no SV
  or CNV calling — in a chromosomal-instability disorder, a material gap rather than a
  rounding error (TD-03).
- **The workspace is any directory outside the repo and outside a cloud-synced root.**
  `export MVA_WORKSPACE=~/mva-workspace` is enough. An encrypted volume is the stronger
  choice, and nothing enforces it ([why](docs/privacy-model.md)).
- Gated challenge data needs your own accepted terms (`docs/submission-runbook.md` §1.3).

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
- No human domain expert has reviewed any of this. Adversarial model review is a filter, not
  a substitute, and treating model fluency as evidence quality is the exact failure this
  project exists to resist (TD-09).

## Further reading

| | |
|---|---|
| Commands, invariants, gate timings | `AGENTS.md` |
| The rules, with IDs | `docs/golden-principles.md` |
| Why a choice was made | `docs/decisions/` (29 records) |
| Architecture and layering | `docs/architecture.md` |
| Threat model, and deletion on macOS | `docs/privacy-model.md` |
| What's real and what's mocked | `docs/maturity-ledger.md` |
| Deferred work, with costs | `docs/tech-debt.md` |

Track 1 scoring methodology is adapted from the CAGI6 Rare Genomes Project assessment
(Stenton et al., 2024). HPO identifiers in the synthetic fixture are real public ontology
terms; the phenotype profile and all gene symbols are fabricated.
