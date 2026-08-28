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
2. ~~**Wire the real adapters at the composition root.**~~ **DONE 2026-08-28**
   (ADR 0027). `synthetic: false` binds ClinVar, gnomAD, SnpEff and MANE from
   `$MVA_RESOURCES`; `synthetic: true` keeps the demo tables. A real case with
   missing resources now REFUSES rather than falling back.
3. ~~**Normalise indels against the reference FASTA.**~~ **DONE 2026-08-28.** One
   `FastaReference` handle is passed to `normalise_variants` and to **both**
   joining adapters, and `qc/qc_report.json` now carries the typed
   `left_alignment` state.
4. **Run the real case, then read the run's warnings**, not just its output.
   Still open, and now the only remaining step: the proband VCF is outside what
   agents may read (ADR 0008), so this is a human step —
   `docs/submission-runbook.md`.

## Integration backlog — owned by the orchestrator, not by any agent

| # | Change | File | Blocked by |
|---|---|---|---|
| ~~1~~ | ~~Bind ClinVar / gnomAD / SnpEff / MANE adapters into `AdapterSet` — and pass `reference=` to ClinVar and gnomAD~~ **CLOSED 2026-08-28**, ADR 0027. `annotation/binding.py` + `orchestrator._resolve_case_resources`. `reference=` is a keyword with **no default** on `build_real_adapter_set`, and `tests/unit/test_adapter_binding.py` asserts both adapters actually consulted it, not merely that they were constructed. | `orchestrator.py`, `annotation/binding.py` | — |
| 2 | Bind `HpoResourceSet` and pass `semantics=` to `score_all_genes`; thread the real `hpo_version` into `load_phenotype_profile` | `orchestrator.py` | HPO agent in flight |
| ~~3~~ | ~~Load the reference FASTA and pass the lookup into `normalise_variants`~~ **CLOSED 2026-08-28.** Resolved by `mva.resources.reference_fasta_path` (per-case override wins, else the shared release); one handle serves normalisation and both adapters. | `orchestrator.py` | — |
| ~~4~~ | ~~A resource-root setting for out-of-repo releases, with digests~~ **CLOSED** — `ResourceSettings` + ADR 0020; this change added `resources.snpeff_pins` and consumes the whole block. | `config.py`, `config/` | — |
| 5 | `.gitignore` negations for real public-reference fixtures (ADR 0012) | `.gitignore` | fixture agents in flight |
| 6 | Golden expectation locking the **emitted submission EPCRs** | `tests/golden/` | submission ordering settling — see ADR 0014 |
| 7 | Maturity-ledger rows: split HPO (`real`) from `gene_phenotype.tsv` (`synthetic-substitute`); regrade `annotation` per adapter | `docs/maturity-ledger.md` | adapters landing |
| 8 | ASSUMPTION-PHENOTYPE-03/04/05 (directional entailment; IC from corpus not depth; deliberate asymmetry) | `docs/scientific-assumptions.md` | — |
| 9 | Decision record for the SYNTHMUL4 / SYNTHMET2 phenotype movement | `docs/decisions/` | HPO agent's final delta |

Item 6 is the one most likely to be skipped and most costly to skip — see the
enforcement-gap section of ADR 0014.

**Item 1 had a trap, and it is now closed by construction.** Both the ClinVar and
gnomAD adapters take an optional `reference=`. Constructed without one they still
work, still pass their tests, and silently run in their degraded state — which is
the *default*. Measured on the real gnomAD chr21 shard: of 1,382 indel records in
a 520 kb exonic window, 1,029 sit in repeat tracts, and **0 of their right-shifted
spellings join without a reference against 989 with one**. Thirty of those are
variants gnomAD calls common. On ClinVar the same omission costs 1,761
Pathogenic/Likely-pathogenic assertions over one 520 kb window.

`build_real_adapter_set(resolved, *, reference)` therefore takes the keyword with
**no default**, so it cannot be omitted by accident; when it is `None`, both
adapters' own `representation_limitation` is surfaced into the run warnings with
the measured cost attached. Adding a third adapter that joins on coordinates? Give
it the same treatment.

## Facts established 2026-08-28 that change the wiring

Measured, not assumed. Each was found by reproducing a failure first.

- **gnomAD v4.1 exomes has no chrM shard.** `frequencies()` therefore raises on
  any mitochondrial variant *even against a fully complete release*. This is not a
  bug to fix; it is a property of the release. **Handled 2026-08-28** by
  `PartialCoverageFrequencyAdapter`, which drives `lookup_partial()`, leaves the
  unanswerable variants absent, and names the contigs and counts in the run
  warnings (never a coordinate — PRIV-09).
- **SnpEff must be pinned with `SnpEffArtifactPins.from_manifest(...)`, never
  `measure()`.** `measure()` hashes whatever happens to be installed, which pins
  the install to itself and verifies nothing.
- **A complete data file with a short or stale index returned `{}`** while every
  completeness check reported "fine" — a second door past the EOF check. Now
  closed by `tabix_index_reach()`, which compares the index's furthest reachable
  BGZF offset against file size. Note mtime is unusable as a gate here: all 24
  real shards trip htslib's "index older than data" warning and all 24 are sound.
- **The cost of an unreferenced ClinVar join, measured** on the real 2026-08-22
  release over chr17:43,000,000-43,520,000: 3,595 indel ALT alleles, 2,215 (61.6%)
  in repeat tracts, **0** right-shifted spellings joining without a reference
  versus **2,211** with one — of which **1,761 are assertions ClinVar calls
  Pathogenic/Likely pathogenic**. This is what item 1's `reference=` trap costs on
  the clinical slot, not the frequency slot.
- **The evidence spill file carries proband coordinates** (evidence subjects are
  variant IDs) and is therefore SENSITIVE under GP-40. It must be created inside
  the workspace and unlinked on close — never in the repo, never in a temp dir
  that outlives the run.

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
