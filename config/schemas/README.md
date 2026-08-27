# Config schemas

| File | Mirrors | Validated by |
|---|---|---|
| `case.schema.json` | `src/mva/config.py::CaseConfig` | `snakemake.utils.validate(config, ...)` at the top of the `Snakefile` |

## Why a JSON Schema exists when Pydantic already validates

Pydantic is the authority. It runs inside `mva run validate` and it enforces
things a JSON Schema cannot — that the seven positive scoring weights sum to
1.0, that a non-synthetic case may not request `network_profile: online`.

The JSON Schema buys one thing Pydantic cannot: **the error arrives before
Snakemake builds the DAG**, which is before any directory is created and before
any patient file is opened. A typo in a case file should not cost you a
half-written run directory.

So the schema is the early, cheap check, not the only one. It is deliberately
allowed to be *weaker* than the model. It must never be *wider* in the places
that matter — the two below.

## The two rules the schema encodes exactly

**Input paths are workspace-relative.** `$defs/workspaceRelativePath` rejects a
leading `/`, a leading `~`, any `..` segment and any backslash. This is the same
rule as `InputPaths._must_be_relative`, and it exists so that no committed config
and no provenance manifest ever records where patient data lives on disk
(ADR 0006).

**`max_submission_rows` is capped at 10.** The challenge scorer raises
`ValueError` above 10 rows, so an 11-row config is not a degraded submission, it
is a zero (`docs/references/track1-submission-contract.md`).

`synthetic` is `required` rather than defaulted, for the same reason it has no
default in `CaseConfig`: a real case that was silently treated as synthetic would
skip privacy enforcement, and the failure mode of forgetting the key must be a
loud error rather than the unsafe branch.

## The three workflow-only keys

`workspace`, `run_id` and `defaults_config` are accepted by this schema and
rejected by `CaseConfig`, which sets `extra="forbid"`. That asymmetry is
intentional: they are Snakemake overrides passed as `--config key=value`, never
keys in a committed case file. Snakemake merges them into the same dict it
validates, so the schema has to know about them; `mva` never sees them.

## Keeping the two in sync

There is no generator. When `CaseConfig` gains a field, add it here in the same
commit. `tests/unit/test_workflow.py` asserts that `config/synthetic-case.yaml`
still validates, which catches a field being *removed* from the schema but not a
field being *added* to the model — so the sync is a review responsibility, and
`additionalProperties: false` is what makes a forgotten field fail loudly on the
next config that uses it.
