# mva-research

A provenance-first pipeline for rare-disease variant prioritisation and
mechanism-grounded drug-repurposing hypotheses, built for
**Rare Disease, Real Kid: The MVA Hackathon 2026**.

> **This repository contains no patient data and never will.** Everything
> demonstrable here runs on a fully synthetic case with fictional genes. Real
> patient data lives in an external, encrypted workspace that the code refuses to
> locate inside this repo.
>
> **Nothing here is medical advice.** Track 2 output is a research hypothesis
> requiring pre-clinical validation.

---

## The problem

Mosaic Variegated Aneuploidy is a rare disorder in which a large fraction of a
child's cells carry the wrong number of chromosomes. The mechanism, where it is
understood, runs through the machinery that normally stops a cell dividing before
its chromosomes are correctly attached — so when that machinery is impaired,
cells divide early, chromosomes missegregate, and aneuploidy accumulates.

The challenge provides one real child's clinical phenotype, VCF and raw reads,
and asks two things:

- **Track 1** — a ranked prediction of the causal variant, or compound-heterozygous
  variant pair.
- **Track 2** — a scientifically rigorous drug-repurposing hypothesis grounded in
  the disrupted mechanism.

The public evaluation code shows the Track 1 answer key is a clinically validated
compound-heterozygous pair. We treat that as an inheritance-model *prior* while
keeping the architecture able to represent other models and other variant classes
— a prior that becomes a filter has stopped being a prior.

## The specific way this problem goes wrong

Both tracks have a characteristic failure mode, and the architecture is built
to resist them rather than to produce an answer quickly.

**Track 1 fails by filtering.** The conventional pipeline drops everything above a
frequency cut-off, drops non-PASS calls, drops low-impact consequences, then ranks
the survivors. A variant removed at step one leaves no trace, so a causal variant
that was filtered is indistinguishable from one that never existed. We separate
filtering from ranking: hard filters remove only records that are *invalid or
impossible*; everything else is flagged and down-ranked (ADR 0005).

**Track 2 fails by sign error.** In a checkpoint-*deficiency* disorder, the drugs
that surface at the top of any target-proximity or pathway search are checkpoint
*inhibitors* — compounds developed to push chromosomally unstable tumour cells
past their aneuploidy-tolerance ceiling. They bind exactly the proteins named in
the mechanism. But the patient's non-tumour cells already sit above that ceiling
with no reserve. The proximity is high and the sign is inverted.

So direction of effect is a mandatory, signed, type-checked field, and a drug
whose observed direction disagrees with the required correction **cannot be
constructed as an accepted hypothesis** — the validator refuses. No weight change
can resurrect a contraindicated compound.

Run against the real BUB1B literature, that gate does what it was built for: of
13 real agents, an MPS1 inhibitor comes out **last** with `WRONG_DIRECTION`, five
agents are rejected for pushing the disease direction, and the highest-scoring
agent in the catalogue is rejected anyway because it is not an approved medicine
and nobody has measured whether it worsens chromosomal instability. The full case,
with every citation, is `docs/track2-hypothesis.md`.

## Architecture

One pipeline, two exits:

```
variant → variant pair → gene → mechanism → drug hypothesis
                  │                              │
            Track 1 exits here            Track 2 continues
```

Track 1 leaves after candidate-pair ranking. Track 2 continues from a ranked pair
through mechanism characterisation, intervention discovery, direction checking,
safety filtering and experimental design. There is no separate Track 2 codebase,
because a drug hypothesis not anchored to the ranked pair is not grounded in
anything.

Three ideas carry the design:

**A vector score, not a number.** A `CandidatePair` carries eight visible
component scores — analytical validity, rarity, consequence, inheritance,
phenotype, mechanism, evidence quality, contradiction penalty. "Why is this ranked
first?" is answerable without re-running anything.

**No claim without evidence.** Every score contribution and every report sentence
resolves to an `EvidenceItem` carrying its source, tool, version, and its
mandatory **limitations**. A renderer refuses to emit an unsourced assertion.
Evidence that *contradicts* a hypothesis is stored alongside evidence that
supports it and is surfaced, never discarded.

**Deterministic and local.** The patient-data path contains no model inference.
Repeat runs are byte-identical under the frozen demo clock, which is checked, not
asserted; for a real case the scientific content is identical and only recorded
timestamps move (scope: `docs/handoff-integrity.md` §4).

Full detail: `docs/architecture.md`. Decisions and their rationale:
`docs/decisions/`.

## Setup

Requires macOS or Linux, [`uv`](https://docs.astral.sh/uv/), and `just`.

```bash
git clone <this-repo> && cd mva-research
just bootstrap        # installs Python 3.12, syncs deps, installs the pre-commit privacy hook
just verify           # the gate: lint + typecheck + tests + architecture + docs + privacy audit
```

`just bootstrap` needs no genomics tooling. `cyvcf2`, `pysam` and Snakemake all
install cleanly on Apple Silicon (verified), and the pipeline falls back to
pure-Python backends where they are absent.

## The synthetic demo

```bash
just demo               # end-to-end synthetic run
just demo-artifacts     # list what it produced
just demo-determinism   # run twice, compare artifact hashes
just clean-demo
```

The fixture is built adversarially rather than as a happy path. It contains a
known compound-heterozygous answer in a fictional gene, plus: common distractors
in the same gene, a low-quality call that must be flagged rather than deleted, an
**in-cis** pair with phase-set evidence that must not be called a compound
heterozygote, a multiallelic site requiring correct per-allele depth assignment, a
phenotype profile exercising all four observation states, and a drug catalogue
containing a right-target/wrong-direction agent, a tool compound, a symptomatic
agent, a context-dependent agent and an off-mechanism agent.

If the pipeline gets this right, it is right for structural reasons — not because
it was tuned to a leaderboard.

`just demo` produces: run manifest, normalised variants, annotated variants,
ranked candidate pairs, evidence database, Track 1 submission CSV, candidate
dossier, mechanism report, ranked drug hypotheses, the rejection record for the
wrong-direction drug, a Track 2 report, a provenance manifest and a privacy-audit
report.

## Privacy boundaries

The asset is a child's genome. Disclosure is not revocable, and the adversary is
mostly accident — git, the shell, the network stack, and macOS itself.

- Patient data lives in an **external workspace** (`MVA_WORKSPACE`). The config
  layer refuses a path resolving inside the repo (symlinks resolved) or under a
  cloud-synced root. On macOS `~/Desktop` and `~/Documents` are iCloud-synced by
  default — the most likely silent disclosure path on this platform.
- `.gitignore` denies genomic and tabular formats **repo-wide**, with narrow
  negations for synthetic fixtures that the audit independently re-validates.
- A pre-commit hook blocks any commit failing the privacy audit.
- The privacy scanner reports **paths, line numbers and span lengths — never
  matched content.** Its own output is a leak vector, because an agent runs it.
- Logs cannot carry genomic records: a redaction filter is installed on every
  handler, plus a record factory.
- Public export is gated twice — an allowlist *and* a post-render content re-scan.
  Classification is a claim; the scan is the verification.
- Annotation is local. A structural test forbids importing a network client
  anywhere on the patient-data path.

Threat model, controls and macOS-specific deletion guidance:
`docs/privacy-model.md`.

## What is deterministic, and what is hypothetical

**Deterministic** — ingestion, normalisation, QC, annotation lookup, filtering,
pairing, phase inference, scoring, ranking, submission rendering, evidence
persistence, provenance. Same inputs and config produce identical scientific
content; artifact bytes are identical too under a fixed clock (TD-21).

**Heuristic** — every scoring weight. They were chosen by reasoning about a severe
paediatric recessive disorder and are **not calibrated against any labelled
dataset**. No composite score is a probability.

**Hypothetical** — the entire Track 2 output. A mechanism chain with an inferred
link is labelled as inferred; a drug hypothesis is a proposal for an experiment,
not a treatment.

**Mocked** — annotation, and the drug catalogue and mechanism library **that the
demo runs on**, are synthetic stand-ins for VEP/gnomAD/ClinVar/DrugBank.
`docs/maturity-ledger.md` grades every module `real`, `synthetic-substitute` or
`stub`. A synthetic substitute is never described as biologically valid.

**Literature-derived** — neither machine-acquired nor synthetic, so ADR 0022
makes it a third kind of knowledge source alongside those two. `knowledge/literature/bub1b/`
holds a real BUB1B/BubR1 mechanism chain and a 13-agent drug catalogue in which
every row is cited to a PMID, checked in both directions by test. The scientific
case built on it is `docs/track2-hypothesis.md`. It is real biology and it is still
a hypothesis: nothing in it is medical advice, and no patient data was used to
produce it.

## Connecting real data later

Real data never enters git. The sequence:

1. Create an encrypted APFS sparse bundle and mount it (`docs/privacy-model.md`).
2. `export MVA_WORKSPACE=/Volumes/MVACASE/case01` — outside the repo, outside
   cloud-synced directories. The config layer verifies this and refuses otherwise.
3. Place the VCF and phenotype under `$MVA_WORKSPACE/inputs/`.
4. Write a case config with `synthetic: false` and `network_profile:
   offline_enforced`. A non-synthetic case with the network open is rejected by a
   validator.
5. Replace the synthetic annotation adapters with real, hash-pinned local
   resources (TD-01). **Until this is done, no real-data claim is supportable.**
6. Run under an OS-level network control as well as the Python guard.
7. `just privacy-audit` before any export.

## Limitations

Stated plainly, because the alternative is implying otherwise:

- **No known-answer validation.** The synthetic case proves the machinery, not the
  calibration. No real-world sensitivity claim is supportable yet (TD-02).
- **No structural or copy-number variant calling.** In a chromosomal-instability
  disorder this is a material gap, not a rounding error (TD-03).
- **Proband-only.** No segregation, no de-novo detection, no trio phasing — so
  phase stays UNKNOWN for every pair, which caps the inheritance component (TD-07).
- **Phenotype matching is flat term overlap** with no ontology-graph propagation,
  so true matches are under-scored (TD-04).
- **Network denial is a tripwire, not a boundary.** C extensions and subprocesses
  bypass the Python guard (TD-06).
- **No human domain-expert review.** Adversarial model review is a filter, not a
  substitute — and treating model fluency as evidence quality is precisely the
  failure this project exists to resist (TD-09).

## Deletion obligations

The challenge data-use terms require deletion **within 30 days of close
(2026-10-24)** from all environments including derived datasets, confirmed by
email to the organizers. `docs/privacy-model.md` documents the procedure and is
honest about what macOS can and cannot guarantee: on APFS with copy-on-write, SSD
wear-levelling and Time Machine local snapshots, overwrite-based erasure is false
assurance. The only defensible guarantee is cryptographic erasure — work inside an
encrypted volume from day one, then destroy the key.

Nothing in this repository is subject to that obligation; it is all synthetic.

## Acknowledgements

Built for the MVA Hackathon 2026 (Sage Bionetworks). Track 1 scoring methodology
is adapted from the CAGI6 Rare Genomes Project assessment (Stenton et al., 2024).
HPO identifiers used in the synthetic fixture are real public ontology terms; the
phenotype profile and all gene symbols are fabricated.
