# Current status

**Updated:** 2026-08-27 · **Phase:** 5 (verification) — Phase 4 adversarial
review in progress.

## What exists

The complete vertical slice, running end to end on the synthetic case.

- **Domain contracts** — `src/mva/models/`: variant, genome/coordinate, evidence,
  candidate pair, phenotype, mechanism, drug, provenance. Frozen,
  `extra="forbid"`, with invariants enforced at the type boundary.
- **Stages** — ingestion (cyvcf2 + pure-Python backends, multiallelic splitting,
  trimming, QC flagging), annotation adapters over hash-pinned local tables,
  four-valued phenotype scoring, filter/pair/score/rank prioritisation, mechanism
  library, drug triage with signed direction checking, DuckDB+Parquet evidence
  store, reporting (Track 1 CSV, dossier, mechanism and Track 2 reports).
- **Composition root** — `pipeline.py` (artifacts, provenance, determinism),
  `orchestrator.py` (stage graph), `cli.py` (Typer). Stages never import each
  other.
- **Workflow** — an 11-job Snakemake DAG whose rules are one `mva` subcommand
  each; the same DAG runs without Snakemake via `mva run all`.
- **Privacy** — scanner (11 checks), redaction filter, netguard, export gate,
  deny-by-default `.gitignore`, pre-commit hook.
- **Structural lints** — import layering, no peer-stage imports, models-as-leaf,
  no network clients on the sensitive path, no wall-clock reads, docs integrity.
- **Documentation** — architecture, 9 ADRs, scientific assumptions, privacy
  model, validation plan, maturity ledger, tech debt, golden principles, the
  vendored Track 1 contract, 5 re-runnable reviewer briefs.

## What works

Verified by running, not asserted:

| Command | Result |
|---|---|
| `just verify` | **all gates passed** |
| `uv run pytest -q` | **296 passed** |
| `uv run ruff check src tests` | clean |
| `uv run pyright src tests` (strict) | 0 errors |
| `just demo` | 30 artifacts, all promised outputs produced |
| `just demo-determinism` | **30 artifacts byte-identical across two runs** |
| `just privacy-audit` | **11 checks pass** |
| `just snakemake -n` | 11-job DAG builds; full run completes |

**Both acceptance criteria confirmed in the artifacts, not just in tests:**

1. The synthetic causal pair `chr15:40200000 C>T` + `chr15:40210500 G>A`
   (SYNTHKIN1) is **row 1** of `track1_submission.csv` at epcr REDACTED-EPCR, produced
   by general scoring with nothing special-cased.
2. The wrong-direction agent `SYNTH-DRUG-B` is **rejected**, leading reason
   `wrong_direction`, with the reversal condition recorded as "nothing short of
   the mechanism itself being wrong".

## What is incomplete

- **Phase 4 adversarial review** — four reviewers running (genomic validity,
  pharmacology, privacy, reproducibility). Findings will be promoted into tests
  or lints, never left as prose.
- **Annotation is a synthetic substitute.** Consequence and frequency values are
  fabricated. This is the blocking prerequisite for any real-data claim (TD-01).
- **No known-answer validation.** The synthetic case proves the machinery, not
  the calibration (TD-02).
- **No structural or copy-number variant calling** (TD-03).
- Full list with costs and triggers: `docs/tech-debt.md`.

## Commands

```bash
just bootstrap        # install toolchain, sync deps, install pre-commit hook
just verify           # THE GATE: lint + typecheck + arch + tests + privacy audit
just demo             # end-to-end synthetic run
just demo-artifacts   # list what it produced
just demo-audit       # show the run's privacy audit
just demo-determinism # run twice, compare artifact hashes
just privacy-audit    # standalone privacy scan
just snakemake track1 # run the DAG to the Track 1 target
just clean-demo
```

The demo workspace is `$TMPDIR/mva-research-demo` — deliberately outside the
repo, because that is what a real run requires and the privacy audit fails an
in-repo workspace.

## Blockers

None. No decision is waiting on the user.

## Notes for the next session

Read this file and `docs/next-actions.md` first. The repo is the system of
record; nothing important lives only in conversation history.

Four integration bugs were found by wiring the stages together rather than by
unit tests, which is worth knowing when judging coverage: an evidence-ID
collision from stamping a pair ID onto a variant-level fact; a low component
score being emitted as `CONTRADICTS`; a `.gitignore` negation that re-opened the
`.env.*` deny rule; and the orchestrator bypassing the public-export gate. All
four are now locked by regression tests.

The one design tension worth knowing: this repo is built to be maximally legible
to an agent **except** for the patient-data path, which is deliberately illegible
to one. If a debugging affordance seems missing, that is usually why — ADR 0008.
