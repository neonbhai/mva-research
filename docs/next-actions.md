# Next actions

## Right now

**First command to run:**
```bash
cd mva-research && just verify
```

If it passes, the repo is in a good state and the next step is below. If it
fails, fix the failure before anything else — `just verify` is a blocking gate by
design (ADR 0009).

## Immediate (Phase 3 completion)

1. Write `src/mva/pipeline.py` (composition root) and `src/mva/cli.py` (Typer).
2. Run `just demo` and fix whatever breaks.
3. Confirm the two acceptance criteria by inspecting artifacts, not by trusting
   tests: the synthetic causal pair ranks **first**, and the wrong-direction drug
   is **rejected with `WRONG_DIRECTION`**.
4. Run `just demo-determinism` and confirm byte-identical artifacts.

## Phase 4 (adversarial review)

Run each brief in `prompts/` as a separate reviewer with no stake in the
implementation:

```
prompts/review-architecture.md
prompts/review-genomic-validity.md
prompts/review-pharmacology.md
prompts/review-privacy.md
prompts/review-reproducibility.md
```

**Every accepted finding is promoted into a test or a lint — never left as prose.**
A finding that cannot be mechanically enforced has not really been fixed.

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
