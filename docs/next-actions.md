# Next actions

## Right now

**First command to run:**
```bash
cd mva-research && just verify
```

If it passes, the repo is in a good state and the next step is below. If it
fails, fix the failure before anything else — `just verify` is a blocking gate by
design (ADR 0009).

## Immediate

Phases 1-5 are complete: the pipeline runs end to end, `just verify` is green,
and both acceptance criteria are confirmed in the artifacts.

The highest-value next step is **TD-02, known-answer validation** — see
`docs/validation-plan.md` Level 3. Until it is done, the synthetic case proves
the machinery and nothing about calibration, and no real-world sensitivity claim
is supportable.

## Phase 4 review: one brief not yet run

Four of the five reviewer briefs have been run and their findings integrated
(genomic validity, pharmacology, privacy, reproducibility). **The architecture
brief has not been run:**

```
prompts/review-architecture.md
```

Run it the same way, and promote every accepted finding into a test or a lint —
never leave one as prose. A finding that cannot be mechanically enforced has not
really been fixed. The other four briefs are re-runnable and worth re-running
after any significant change.

## Before connecting real data — the exact sequence

Do these in order. Steps 1–4 are privacy prerequisites; step 5 is a scientific
prerequisite.

1. Create and mount the encrypted APFS sparse bundle
   (`docs/privacy-model.md` has the commands). Exclude it from Time Machine and
   Spotlight.
2. `export MVA_WORKSPACE=/Volumes/MVACASE/case01` and
   `export TMPDIR=/Volumes/MVACASE/tmp`.
3. Write the real case config with `synthetic: false` and `network_profile:
   offline_enforced`. Run under an OS-level network control as well as the Python
   guard — the Python guard is a tripwire, not a boundary (TD-06).
4. `just privacy-audit` — must pass clean.
5. **Replace the synthetic annotation adapters with real, hash-pinned local
   resources** (VEP or SnpEff, a gnomAD sites-only release, a ClinVar release
   VCF). Until this is done the pipeline's consequence and frequency values are
   fabricated and **no real-data claim is supportable** (TD-01).

Then, and only then, run the real case.

## Before submitting

- Re-verify the Track 1 contract against the live Space. Ours is vendored in
  `docs/references/track1-submission-contract.md` and was correct on 2026-08-27,
  but the Space is mutable.
- Confirm `proband_id` is exactly `PROBAND01` and every chromosome carries the
  `chr` prefix — the scorer compares chrom strings raw, so a missing prefix scores
  zero while looking correct.
- Confirm ≤10 rows and every `epcr` in `(0, 1]`.
- Get human domain-expert review (TD-09).

## Deliberately not doing

- Tuning weights against the live leaderboard (ASSUMPTION-SCORING-03).
- Adding a graph database (ADR 0002) or an LLM to the runtime path (ADR 0003).
- Relaxing the blocking gate for throughput (ADR 0009).
