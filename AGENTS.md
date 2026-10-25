# AGENTS.md

Router, not an encyclopedia. The `docs/` tree is the system of record; this file
tells you which door to open.

One file, on purpose. `CLAUDE.md` is a single `@AGENTS.md` import and holds no
content of its own, because two files saying the same thing stay synchronised
right up until they don't. TD-17 is this repo's own worked example: two modules
agreed on a column name while disagreeing about what it meant, and the table was
dead — no code path could use it. If you are adding guidance for a tool, add it
here.

## Run it

```bash
just all      # bootstrap, then the gate, then the synthetic demo
```

That is the whole answer for a fresh clone. Needs macOS or Linux, `uv` and
`just`. No patient data, no credentials, no encrypted volume. The network is used
once, by `bootstrap`, to fetch dependencies; no run uses it.

## What this is

A provenance-first pipeline for the MVA Hackathon 2026:
`variant → variant pair → gene → mechanism → drug hypothesis`.
Track 1 consumes it through candidate-pair ranking. Track 2 continues from a
ranked pair through mechanism, intervention, direction and safety checks.

The patient-data path is **deterministic and local**. An agent is the engineer
building the machinery, never a runtime component reading the genome.

## The gate

`just verify` is the acceptance gate (ADR 0009). **It is green right now and the
job is to keep it that way.** Do not merge past a failure; determinism and
privacy checks are blocking by design.

Seven steps, in this order. Measured on macOS/Apple Silicon with a warm `uv`
cache, 2026-10-25 — your absolute numbers will differ, the proportions will not:

| Step | What it runs | Measured |
|---|---|---|
| `lint` | `ruff check` + `ruff format --check` | under 1s |
| `typecheck` | `pyright src tests`, pinned 1.1.411 | 8s |
| `arch` | layering + docs-integrity lints, 23 tests | 2s |
| `test` | the full suite: 1817 passed, 13 skipped | **5m07s** |
| `workflow-check` | Snakemake dry run; asserts the DAG resolves | 1s |
| `demo-determinism` | two runs, 33 artifacts hash-compared | 11s |
| `privacy-audit` | `mva privacy audit --repo .` | 5s |
| | **`just verify` end to end** | **5m34s** |

`just demo` on its own is 7s and writes the same 33 artifacts. `just all` is
`verify` plus that, plus a `bootstrap` whose cost is entirely whether the
dependencies are already downloaded.

Two things an agent gets wrong here:

- **Patience.** The test step runs for five minutes behind a progress bar.
  Killing it yields a failure report about a run that had not failed.
- **Making the gate agree with the code instead of the reverse.**
  `PYRIGHT_PYTHON_FORCE_VERSION` is pinned to `1.1.411` in the justfile because
  the `pyright` PyPI package is only a downloader and would otherwise fetch
  whatever is newest — a gate whose rules arrive from the network is not a gate.
  Do not bump that pin, and do not loosen `pyproject.toml`, to make
  `just typecheck` pass. Fix the code, or make the version bump a deliberate,
  committed change with the fixes alongside it.

## Non-negotiable invariants

Full list with IDs: `docs/golden-principles.md`. The ones you will hit first:

1. **GP-11** A genomic coordinate without a genome build is invalid.
2. **GP-13** Filtering ≠ ranking. Hard filters remove only invalid/impossible
   records. Weak candidates are flagged and down-ranked. Where a stage does drop a
   valid record — selection and the annotation relevance bar, and only those — it
   writes one row to the discard ledger naming the variant, the rule and the
   value. Accounted for, not retained (ADR 0028).
3. **GP-14** Absence of information is not negative information. `NOT_ASSESSED`
   ≠ `EXCLUDED`. No frequency data ≠ rare.
4. **GP-15** Phase is never assumed. `UNKNOWN` survives end to end.
5. **GP-16** Direction of effect is mandatory and signed. Tri-state: `None`
   means undeterminable and is *not* agreement.
6. **GP-10 / GP-17** No claim without an `EvidenceItem`; every item states its
   limitations.
7. **GP-30** Repeat runs are byte-identical *under a fixed clock*. No RNG, no
   set/dict-order dependence, no time-derived id or hash. Timestamps come from
   `mva.clock`, which is `SystemClock` for a real case — so a real run's
   scientific content is identical and eleven artifacts differ in recorded time
   only. Scope: `docs/handoff-integrity.md` §4, TD-21.
8. **GP-40** Patient data never enters the repo. The workspace is external.
9. **GP-41** Privacy tooling reports paths and counts, never matched content.

## Where things live

| Need | Path |
|---|---|
| Architecture and layering | `docs/architecture.md` |
| Enforceable rules (GP-nn) | `docs/golden-principles.md` |
| Why a choice was made | `docs/decisions/` |
| Scientific claims and their caveats | `docs/scientific-assumptions.md` |
| Threat model and controls | `docs/privacy-model.md` |
| What is real vs mocked | `docs/maturity-ledger.md` |
| Challenge submission contract | `docs/references/track1-submission-contract.md` |
| Track 2 scientific case (BUB1B) | `docs/track2-hypothesis.md` |
| Every cited claim, re-graded | `docs/track2-claims-ledger.md` |
| Literature-derived knowledge tables | `knowledge/literature/bub1b/` |
| Current state, blockers, commands | `docs/current-status.md` |
| What to do next | `docs/next-actions.md` |
| Deferred work with costs | `docs/tech-debt.md` |
| Integrity fixes that cross owners | `docs/handoff-integrity.md` |
| Container build and parity assertion | `workflow/containers/README.md` |
| Reviewer briefs (re-runnable) | `prompts/` |
| Domain contracts | `src/mva/models/` |

## Commands

```bash
just all            # bootstrap + verify + demo, the whole thing
just bootstrap      # install toolchain + sync deps
just verify         # THE GATE: lint + typecheck + arch + tests + workflow + determinism + privacy
just demo           # end-to-end synthetic run
just demo-determinism  # run the demo twice, compare artifact hashes
just privacy-audit  # standalone privacy scan
just science-canary # 12 controls proving the direction gate can still say yes
just ranking        # regenerate docs/track2-ranking.md from a live run
just drugs-to-avoid # regenerate docs/track2-drugs-to-avoid.md from a live run
just clean-demo     # remove demo workspace
```

```bash
mva verify reproduces <submitted.csv> <produced.csv>
```

Answers whether a run reproduces a banked Track 1 submission. Both arguments are
patient data; it emits counts, flags and a correlation only, never a row, so it is
safe for an agent to run under GP-44 and ADR 0008. Exits non-zero on anything
short of an exact reproduction. See `docs/submission-runbook.md` §4a.

## Setting up a case

`mva case` is the front door. It does not require reading the schema first.

```bash
mva case init                      # guided; --no-input for scripts, --run to run it
mva case check -c <case>.yaml      # says whether it could run; runs nothing
mva case phenotype-template        # writes a phenotype TSV header to fill in
```

What a case needs, if you are writing the YAML by hand
(`config/schemas/case.schema.json` is the authority):

- Required: `case_id`, `proband_id`, `genome_build`, `synthetic`, `inputs`.
  `inputs` requires `vcf` and `phenotype`, both workspace-relative.
- `genome_build` is `GRCh38` or `GRCh37`. GRCh37 exists so mismatches can be
  *rejected*, not silently accepted (GP-11).
- `synthetic` is explicit and never inferred. A real case mislabelled synthetic
  would skip privacy enforcement.

Inputs, as the code actually consumes them:

- **Phenotype is a TSV.** `hpo_id`, `label`, `status` required; `onset`,
  `provenance`, `extraction_confidence`, `notes` optional
  (`src/mva/phenotype/loader.py`). The provenance and confidence columns are
  there because terms are extracted from a clinical document — that is the
  intended workflow, not a workaround for a missing structured feed.
- **The VCF goes in as supplied**: plain `.vcf`, `.gz`/`.bgz`/`.bgzf`, or BCF.
  A `.tbi`/`.csi` index next to it is used for per-contig streaming.
- **Raw reads are not consumed at all.** No BAM, CRAM or FASTQ path exists:
  no alignment, no variant calling, no SV/CNV calling (TD-03). Nothing will
  produce a structural-variant result, so do not wait for one.
- Gated challenge data, if you have accepted the terms:
  `hf download SageBio/mva-hackathon-2026-data --repo-type dataset --local-dir <outside-the-repo>`
  (`docs/submission-runbook.md` §1.3).

## Rules for changing things

- **Scoring weights** are configuration (`config/default.yaml`), not code.
  Changing one needs a decision record, a test, and a before/after comparison.
  Golden expectations in `tests/golden/` are **never silently re-baselined**
  (GP-32). If a golden file has to move, the diff and the reason are the change;
  re-baselining to get a green run destroys the only record of what the old
  answer was.
- **Two documents are generated, not written**: `docs/track2-ranking.md` and
  `docs/track2-drugs-to-avoid.md`. Hand-editing either is a wasted change —
  `tests/unit/test_ranking.py` and `tests/unit/test_drugs_to_avoid.py` fail when a
  page stops matching a live run. Change the generator or the inputs, then
  `just ranking` / `just drugs-to-avoid`.
- **`docs/architecture.md` is canonical.** `docs/architecture.html` is a view of
  it. Edit the Markdown; never let a fact exist only in the HTML.
- **New stage?** It goes in its own package and does not import peer stages.
  Composition happens only in `src/mva/pipeline.py` (GP-03).
- **New external tool or data source?** It gets an adapter behind a Protocol in
  the relevant package, plus a row in `docs/maturity-ledger.md` grading it
  `real` / `synthetic-substitute` / `stub`. Never describe a synthetic
  substitute as biologically valid (GP-20).
- **Custom lints** must put the remediation in the failure message. The reader
  may be an agent whose only view of the rule is that string.
- Review findings are promoted into a test or a lint, never left as prose.

## Working with the synthetic case

Everything demoable runs on `config/synthetic-case.yaml` and
`tests/fixtures/synthetic/`. Genes prefixed `SYNTH*` are fictional. Debugging is
always done on the synthetic case — never paste real case data anywhere.

## Privacy boundary

`MVA_WORKSPACE` points outside this repo. The config layer refuses a workspace
inside the repo tree or under a cloud-synced directory (`~/Desktop` and
`~/Documents` are iCloud-synced by default on macOS). It checks location only;
encryption is a recommendation it cannot verify, and it is not required to run.
Do not read, cat, head or grep files under the workspace. Introspection commands
operate on synthetic runs and aggregate counts only (GP-44, ADR 0008).
