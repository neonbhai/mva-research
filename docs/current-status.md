# Current status

**Updated:** 2026-08-28 · **Phase:** 6 (handoff). Phases 1-5 complete.

## What exists

A complete, provenance-first vertical slice running end to end on a synthetic
case, with Phase 4 adversarial-review findings integrated.

- **Domain contracts** — `src/mva/models/`: variant, genome/coordinate, evidence,
  candidate pair, phenotype, mechanism, drug, provenance. Frozen,
  `extra="forbid"`, invariants enforced at construction rather than downstream.
- **Stages** — ingestion (cyvcf2 + pure-Python backends, multiallelic splitting,
  parsimony trimming, QC flagging), annotation adapters over hash-pinned local
  tables, four-valued phenotype scoring, filter/pair/score/rank prioritisation,
  mechanism library, drug triage with signed direction checking, DuckDB+Parquet
  evidence store, reporting (Track 1 CSV, dossier, mechanism and Track 2 reports).
- **Composition root** — `pipeline.py`, `orchestrator.py`, `cli.py`. Stages never
  import each other; an AST lint enforces it.
- **Workflow** — an 11-job Snakemake DAG whose rules are one `mva` subcommand
  each. The same DAG runs without Snakemake via `mva run all`.
- **Privacy** — 11-check scanner, redaction filter, netguard, export gate,
  deny-by-default `.gitignore`, pre-commit hook.
- **Documentation** — architecture, 11 ADRs, scientific assumptions, privacy
  model, validation plan, maturity ledger, tech debt, golden principles, the
  vendored Track 1 contract, 5 re-runnable reviewer briefs.

## What works

Verified by running, not asserted:

| Command | Result |
|---|---|
| `just verify` | **all gates passed** |
| `uv run pytest -q` | **411 passed** |
| `uv run ruff check src tests` | clean |
| `uv run pyright src tests` (strict) | 0 errors |
| `just demo` | 30 artifacts, every promised output produced |
| `just demo-determinism` | **30 artifacts byte-identical** |
| cross-process determinism | **identical under differing PYTHONHASHSEED and TZ** |
| `just privacy-audit` | **11 checks pass, 0 fail-severity findings** |
| `just workflow-check` | 11-job DAG resolves |

**Both acceptance criteria confirmed in the artifacts, not just in tests:**

1. The synthetic causal pair `chr15:40200000 C>T` + `chr15:40210500 G>A`
   (SYNTHKIN1) is **row 1** of `track1_submission.csv` at epcr REDACTED-EPCR, produced
   by general scoring with nothing special-cased.
2. The wrong-direction agent `SYNTH-DRUG-B` is **rejected**, leading reason
   `wrong_direction`, reversal condition recorded as "nothing short of the
   mechanism itself being wrong".

## What Phase 4 changed

Four adversarial reviewers, each reproducing claims by running code. They found
real defects, including one in a test written by the lead engineer.

- **The central Track 2 claim was breakable.** Nothing validated that
  `disease_direction`, `required_correction` and the target node's state agreed,
  so a one-cell typo inverted the whole direction gate — the contraindicated
  checkpoint inhibitor rendered as "direction AGREES". Now a model validator.
- **A wrong-direction drug could be silently erased.** `model_copy` bypasses
  pydantic validators; the renderer then matched such a hypothesis to no section
  at all. Now re-validated, and the renderer raises on an unrendered hypothesis.
- **The privacy content allowlist was an unconditional bypass.** A real VCF
  committed under `knowledge/public/` or `tests/golden/` passed every check.
- **Case-variant APFS paths defeated workspace containment.** Now compares
  `(st_dev, st_ino)`.
- **The Track 1 export gate discarded its own verdict.**
- **The multiallelic allele-balance correction was being thrown away** by the two
  stages that act on it, costing 0.178 composite on VCF formatting alone.
- **The headline Track 2 golden test was a tautology** — it compared the
  expectation file against itself and stayed green with the direction check
  inverted. Now asserts on the pipeline result.
- **`config/default.yaml` and both golden TSVs were unpinned**, so a weight edit
  or a re-baseline passed `just verify` unnoticed. Now sha256-locked.

One finding forced a deliberate behaviour change with a decision record:
ADR 0011 makes an unassessed chromosomal-instability risk disqualifying, which
re-baselined SYNTH-DRUG-E and required a new catalogue entry (SYNTH-DRUG-G) to
preserve the "direction undetermined is not disagreement" lesson.

## What is incomplete

- **Annotation is a synthetic substitute.** Consequence and frequency values are
  fabricated. This is the blocking prerequisite for any real-data claim (TD-01).
- **No known-answer validation.** The synthetic case proves the machinery, not
  the calibration. No real-world sensitivity claim is supportable (TD-02).
- **No structural or copy-number variant calling** (TD-03); no repeat expansions
  (TD-13); no mtDNA heteroplasmy (TD-14).
- **Proband-only**: phase is UNKNOWN for every pair, capping the inheritance
  component (TD-07).
- **No human domain-expert review** (TD-09).
- Full list with costs and triggers: `docs/tech-debt.md`.

## Commands

```bash
just bootstrap        # install toolchain, sync deps, install pre-commit hook
just verify           # THE GATE: lint, typecheck, arch, tests, workflow,
                      #   determinism, privacy audit
just verify-fast      # inner loop; skips the two end-to-end runs
just demo             # end-to-end synthetic run
just demo-artifacts   # list what it produced
just demo-audit       # show the run's privacy audit
just demo-determinism # run twice, compare artifact hashes
just privacy-audit    # standalone privacy scan
just snakemake track1 # run the DAG to the Track 1 target
just clean-demo
```

The demo workspace is `$TMPDIR/mva-research-demo` — outside the repo, because
that is what a real run requires and the privacy audit fails an in-repo
workspace.

## Blockers

None. No decision is waiting on the user.

## Notes for the next session

Read this file and `docs/next-actions.md` first. The repo is the system of
record; nothing important lives only in conversation history.

Worth knowing when judging test coverage: **most of the serious defects were
found by wiring stages together and by adversarial review, not by unit tests.**
The unit suite was green while the direction gate was invertible, the export gate
was inert, and a real VCF could be committed. Every one of those is now locked by
a regression test — but the lesson is that a green suite is evidence about the
tests, not about the system.

The one design tension worth knowing: this repo is built to be maximally legible
to an agent **except** for the patient-data path, which is deliberately illegible
to one. If a debugging affordance seems missing, that is usually why — ADR 0008.
