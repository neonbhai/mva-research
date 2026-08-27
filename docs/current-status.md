# Current status

**Updated:** 2026-08-27 · **Phase:** 3 (synthetic vertical slice) — integration in
progress.

## What exists

- **Domain contracts** — `src/mva/models/`: variant, genome/coordinate, evidence,
  candidate pair, phenotype, mechanism, drug, provenance. Frozen, strictly
  validated, `extra="forbid"`.
- **Foundation** — injectable clock, canonical hashing, typed errors, typed config
  with the workspace privacy boundary.
- **Structural lints** — `tests/unit/test_architecture.py` enforces the import
  layer map, the no-peer-stage-imports rule, models-as-leaf, no-network-clients on
  the sensitive path, and no direct wall-clock reads. Failure messages carry their
  own remediation.
- **Privacy protections** — deny-by-default `.gitignore` with narrow, audited
  negations; privacy package (scanner, redaction, netguard, export gate).
- **Synthetic fixture** — adversarial by design: known compound-het answer, common
  distractors, low-quality call, in-cis pair, multiallelic site, four-valued
  phenotype profile, and a drug catalogue containing a wrong-direction agent.
- **Knowledge tables** — `knowledge/public/`: consequences, frequencies, gene
  panel, gene–phenotype, mechanism chain, mechanism metadata, drug catalogue. All
  synthetic, all versioned.
- **Golden expectations** — `tests/golden/`: the expected #1 pair and the expected
  drug triage outcomes.
- **Documentation** — architecture, 9 decision records, scientific assumptions,
  privacy model, validation plan, maturity ledger, tech debt, golden principles,
  the vendored Track 1 submission contract, and five re-runnable reviewer briefs.

## What works

Verified by running:

- `uv run pytest tests/unit/test_architecture.py` — 5 passed. The layering test
  caught a real ordering error on its first run (foundation utilities had to move
  below `config`), which is the intended behaviour.
- `uv run ruff check src tests` — clean.
- `uv run pyright src tests` — 0 errors, strict mode.
- `.gitignore` verified with `git check-ignore -v --no-index`: sensitive paths
  ignored, synthetic fixtures correctly re-included through negations.

## What is incomplete

- **Composition root** — `src/mva/pipeline.py` and `src/mva/cli.py` not yet
  written. They are the last integration step and depend on the stage APIs now
  landing.
- **Stage integration** — ingestion, annotation, phenotype, prioritization,
  mechanisms, interventions, evidence and reporting are implemented but not yet
  wired end to end or run together.
- **`just demo`** — not yet executable (needs the CLI).
- **Phase 4 adversarial review** — briefs written, reviews not yet run.
- **Test count** — the 21 required behaviours are specified and mostly
  implemented per-stage; the end-to-end golden and determinism tests are pending
  the composition root.

## Commands

```bash
just bootstrap        # install toolchain, sync deps, install pre-commit hook
just verify           # THE GATE: lint + typecheck + arch + tests + privacy audit
just lint             # ruff check + format check
just typecheck        # pyright strict
just test             # full pytest suite
just arch             # structural + docs-integrity lints only
just demo             # end-to-end synthetic run
just demo-determinism # run twice, compare artifact hashes
just privacy-audit    # standalone privacy scan
just clean-demo       # remove demo workspace
```

## Blockers

None external. No decision is waiting on the user.

## Notes for the next session

Start by reading this file and `docs/next-actions.md`. The repo is the system of
record; nothing important lives only in conversation history.

The one design tension worth knowing about: this repo is built to be maximally
legible to an agent **except** for the patient-data path, which is deliberately
illegible to one. If a debugging affordance seems missing, that is usually why —
see ADR 0008.
