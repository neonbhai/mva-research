# Reviewer brief — reproducibility and tests

You are a research-software engineer reviewing `mva-research` adversarially. Do
NOT modify files.

## What to attack
1. **Determinism (GP-30).** Find every source of nondeterminism: `datetime.now`,
   `time.time`, `random`, `uuid`, `set`/`dict` iteration feeding output order,
   unstable sorts, `hash()` on strings, float repr, parallel iteration, filesystem
   ordering from `glob`/`iterdir`, Parquet metadata, dict ordering in JSON.
   The lint catches direct clock reads — what does it miss?
2. **Test quality.** Which tests would still pass if the implementation were
   subtly wrong? Look for tests asserting a function ran rather than that it was
   correct; tests whose fixture is built by the same code path they test;
   over-mocking that removes the thing under test.
3. **Golden tests.** Is the golden expectation genuinely independent of the
   implementation, or was it generated from it? Would it catch a regression?
4. **The `just verify` gate.** Does it actually cover everything claimed? Can it
   pass while something important is broken?
5. **Provenance.** Does the run manifest capture enough to reproduce a run? Input
   hashes, config hash, git commit AND dirty flag, tool versions, reference
   versions, network profile. What is missing?
6. **Environment pinning.** Is `uv.lock` committed? Are tool versions recorded as
   observed at runtime rather than as declared?
7. **Coverage gaps.** Which of the 21 required test behaviours (see
   `docs/current-status.md`) are covered by a test that does not actually exercise
   the behaviour?

## Output
Severity-ranked findings with file:line and the concrete failure scenario. For
each, say whether it is a test gap or an implementation bug. Every valid finding
must map to a new or strengthened test.
