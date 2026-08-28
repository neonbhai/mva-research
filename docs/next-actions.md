# Next actions

Updated 2026-08-28. The proband data is here; the pipeline is not yet wired to it.

## First command

```bash
just verify
```

Blocking by design (ADR 0009). Fix a failure before anything else.

## The critical path to a submission

In order. Nothing downstream is meaningful until (1) holds.

1. **Gene assignment must reach `VariantRecord.consequences`.**
   `gene_symbols` is a *derived* property reading `csq.gene_symbol`, and
   `pairing._variants_by_gene` groups on it. With no annotator the pipeline emits
   **zero candidate pairs** and scores zero. Two independent implementations are
   in flight: `annotation/gene_intervals.py` (MANE interval join — the robust
   fallback) and `annotation/snpeff_local.py` (full consequence prediction).
   Either unblocks the path; both is better.
2. **Wire the real adapters at the composition root.** They exist and are tested
   but nothing binds them, so the pipeline still runs on synthetic tables. See
   the integration backlog below.
3. **Normalise indels against the reference FASTA** before any frequency or
   clinical join. a substantial fraction of the proband's records are indel-bearing and a failed
   join is indistinguishable from "novel and ultra-rare".
4. **Run the real case, then read the run's warnings**, not just its output.

## Integration backlog — owned by the orchestrator, not by any agent

| # | Change | File | Blocked by |
|---|---|---|---|
| 1 | Bind ClinVar / gnomAD / SnpEff / MANE adapters into `AdapterSet` — **and pass `reference=` to ClinVar and gnomAD** (see below) | `orchestrator.py`, `annotation/__init__.py` | adapter agents in flight |
| 2 | Bind `HpoResourceSet` and pass `semantics=` to `score_all_genes`; thread the real `hpo_version` into `load_phenotype_profile` | `orchestrator.py` | HPO agent in flight |
| 3 | Load `config.inputs.reference_fasta` and pass the lookup into `normalise_variants` | `orchestrator.py`, `config.py` | reference FASTA download; join agent in flight |
| 4 | A resource-root setting for out-of-repo releases, with digests | `config.py`, `config/` | — |
| 5 | `.gitignore` negations for real public-reference fixtures (ADR 0012) | `.gitignore` | fixture agents in flight |
| 6 | Golden expectation locking the **emitted submission EPCRs** | `tests/golden/` | submission ordering settling — see ADR 0014 |
| 7 | Maturity-ledger rows: split HPO (`real`) from `gene_phenotype.tsv` (`synthetic-substitute`); regrade `annotation` per adapter | `docs/maturity-ledger.md` | adapters landing |
| 8 | ASSUMPTION-PHENOTYPE-03/04/05 (directional entailment; IC from corpus not depth; deliberate asymmetry) | `docs/scientific-assumptions.md` | — |
| 9 | Decision record for the SYNTHMUL4 / SYNTHMET2 phenotype movement | `docs/decisions/` | HPO agent's final delta |

Item 6 is the one most likely to be skipped and most costly to skip — see the
enforcement-gap section of ADR 0014.

**Item 1 has a trap.** Both the ClinVar and gnomAD adapters take an optional
`reference=`. Constructed without one they still work, still pass their tests,
and silently run in their degraded state — which is the *default*. Measured on
the real gnomAD chr21 shard: of 1,382 indel records in a 520 kb exonic window,
1,029 sit in repeat tracts, and **0 of their right-shifted spellings join without
a reference against 989 with one**. Thirty of those are variants gnomAD calls
common. Wiring the adapters without `reference=` therefore leaves the exact
defect ADR 0018 was written to close, while looking fully wired.

Surface `representation_limitation` into the run warnings at the same time, so a
run that *is* degraded says so.

## Open review findings not yet closed

From the Codex adversarial review (2026-08-28), routed to file owners:

- gnomAD: incomplete shards return `{}` instead of failing closed — **live while
  the 184.8 GB download is in flight**.
- gnomAD: backend exceptions disclose the queried patient coordinate (verified).
- ClinVar: allele representation not canonicalised before joining (verified).
- MANE: placeholder `impact` can be reported with `synthetic=False` (verified).
- SnpEff: malformed subprocess output becomes silent absence (verified); install
  is not content-pinned.
- HPO: alias-vs-primary status conflicts resolve silently (verified); malformed
  corpus rows are dropped without appearing in degradation stats.

From the earlier architecture review, still open: C1 (three divergent `run_id`
derivations; the Snakemake DAG fails end to end), M1, M3, M4, M5, M6, M8.
**C1 matters before anyone reruns our submission** — reproducibility is a stated
challenge requirement.

## Before submitting

- Re-verify the contract against the live Space. Ours is vendored and was
  re-verified against the real `evaluation.py` on 2026-08-28.
- `proband_id` exactly `PROBAND01`; every chromosome `chr`-prefixed. **The proband
  VCF uses bare contigs** — the conversion is tested, but any new emission path
  must be re-checked. A correct answer emitted as `15` scores exactly zero.
- ≤10 rows, every `epcr` in `(0, 1]`, all distinct.
- Fill all ten rows: rows below the answer cannot lower either metric.
- Never split a candidate pair across two rows (100 → 50 rank points).
- Human domain-expert review (TD-09).

## Deliberately not doing

- Tuning weights against the live leaderboard (ASSUMPTION-SCORING-03).
- Fetching the 79 GB of FASTQs or realigning (`docs/resource-acquisition-assessment.md`).
- gnomAD genomes (524 GB — impossible here) or dbSNP (no allele number, so
  ADR 0010's guard is unimplementable against it).
- A graph database (ADR 0002) or an LLM in the runtime path (ADR 0003).
