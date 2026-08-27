# CLAUDE.md

Router, not an encyclopedia. The `docs/` tree is the system of record; this file
tells you which door to open. Keep it under ~100 lines.

## What this is

A provenance-first pipeline for the MVA Hackathon 2026:
`variant → variant pair → gene → mechanism → drug hypothesis`.
Track 1 consumes it through candidate-pair ranking. Track 2 continues from a
ranked pair through mechanism, intervention, direction and safety checks.

The patient-data path is **deterministic and local**. Claude is the engineer
building the machinery, never a runtime component reading the genome.

## Non-negotiable invariants

Full list with IDs: `docs/golden-principles.md`. The ones you will hit first:

1. **GP-11** A genomic coordinate without a genome build is invalid.
2. **GP-13** Filtering ≠ ranking. Hard filters remove only invalid/impossible
   records. Weak candidates are flagged and down-ranked, never deleted.
3. **GP-14** Absence of information is not negative information. `NOT_ASSESSED`
   ≠ `EXCLUDED`. No frequency data ≠ rare.
4. **GP-15** Phase is never assumed. `UNKNOWN` survives end to end.
5. **GP-16** Direction of effect is mandatory and signed. Tri-state: `None`
   means undeterminable and is *not* agreement.
6. **GP-10 / GP-17** No claim without an `EvidenceItem`; every item states its
   limitations.
7. **GP-30** Repeat runs are byte-identical. No wall clock, no RNG, no
   set/dict-order dependence. Timestamps come from `mva.clock`.
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
| Current state, blockers, commands | `docs/current-status.md` |
| What to do next | `docs/next-actions.md` |
| Deferred work with costs | `docs/tech-debt.md` |
| Reviewer briefs (re-runnable) | `prompts/` |
| Domain contracts | `src/mva/models/` |

## Commands

```bash
just bootstrap      # install toolchain + sync deps
just verify         # THE GATE: lint + typecheck + tests + arch + docs + privacy
just demo           # end-to-end synthetic run
just privacy-audit  # standalone privacy scan
just clean-demo     # remove demo workspace
```

`just verify` is the acceptance gate. Do not merge past a failure; determinism
and privacy checks are blocking by design (see `docs/decisions/0009-*`).

## Rules for changing things

- **Scoring weights** are configuration (`config/default.yaml`), not code.
  Changing one needs a decision record, a test, and a before/after comparison.
  Golden expectations in `tests/golden/` are never silently re-baselined (GP-32).
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
`~/Documents` are iCloud-synced by default on macOS). Do not read, cat, head or
grep files under the workspace. Introspection commands operate on synthetic runs
and aggregate counts only (GP-44, ADR 0008).
