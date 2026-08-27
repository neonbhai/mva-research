# Maturity ledger

What is real, what is a synthetic stand-in, and what is a stub. Machine-readable
intent: `tests/unit/test_docs_integrity.py` asserts every package under
`src/mva/` has a row here.

**Grades**
- `real` — genuine implementation; behaves correctly on real inputs.
- `synthetic-substitute` — a functioning implementation over fabricated data. The
  *machinery* is real; the *biology* is not. Never describe output from one of
  these as biologically valid (GP-20).
- `stub` — placeholder; does not do the job yet.

| Module | Grade | What is real | What is not | Debt |
|---|---|---|---|---|
| `models` | real | All domain contracts, validators, invariants. | — | — |
| `config` | real | Typed config, workspace boundary, cloud-sync rejection. | — | — |
| `clock` / `determinism` / `errors` | real | Injectable clock, canonical hashing, typed errors. | — | — |
| `privacy` | real | Scanner (11 checks), redaction on every handler, export gate acting on its verdict, containment by (st_dev, st_ino). | `netguard` is best-effort and must be labelled so: C extensions and subprocesses bypass it. The declared-synthetic downgrade is defeatable by a deliberate lie, not by accident. | TD-06 |
| `ingestion` | real | VCF parsing (cyvcf2 + pure-Python), multiallelic splitting, trimming, QC flagging. | Left-alignment requires a reference FASTA; without one only trimming is claimed. | — |
| `annotation` | **synthetic-substitute** | Adapter protocols, manifest hash-pinning, multi-transcript preservation, "no data ≠ rare" handling. | **Consequence and frequency values are fabricated.** Not VEP, not gnomAD, not ClinVar, not SpliceAI. | TD-01 |
| `phenotype` | **synthetic-substitute** | Four-valued observation logic, gene–phenotype scoring, contradiction handling. | Gene–phenotype associations are fabricated. No HPO ontology graph, so no ancestor propagation. | TD-04 |
| `prioritization` | real | Hard/soft filter separation, pairing, phase inference, component scoring, ranking. | Weights are uncalibrated heuristics. Phase is UNKNOWN for all pairs without trio data. | TD-02, TD-07 |
| `mechanisms` | **synthetic-substitute** | Chain loading, link grading, inferred-link penalty, target resolution. | The SYNTHKIN1 mechanism is fictional. Shape mirrors real checkpoint biology; content does not. | TD-08 |
| `interventions` | **synthetic-substitute** | Direction checking, tri-state agreement, safety filtering, evidence tiering, rejection recording. | **The drug catalogue is fabricated.** No real compound is evaluated. | TD-08 |
| `evidence` | real | DuckDB schema, Parquet export, idempotent writes, assertion resolution. | Cross-platform byte-identity untested. | TD-05 |
| `reporting` | real | Assertion checking, tier marking, Track 1 CSV against the verified contract, dossier and Track 2 rendering. | Reports are only as good as the synthetic inputs they render. | — |
| `cli` / `pipeline` | real | Composition root, stage orchestration, determinism verification. | — | — |
| Structural / CNV calling | **absent** | — | Not implemented at all. A material gap for a chromosomal-instability disorder. | TD-03 |
| Repeat-expansion calling | **absent** | — | Not implemented at all. Short-tandem-repeat expansions are neither called nor representable, so a negative result at a repeat locus means "not looked at", not "not present". | TD-13 |
| mtDNA heteroplasmy | **absent** | `chrM` calls are recognised and modelled as `InheritanceModel.MITOCHONDRIAL` instead of being described in nuclear two-copy language. | The heteroplasmy fraction is never measured or represented, and neither is the sampled tissue. Every mitochondrial candidate scores the same flat inheritance value regardless of load. | TD-14 |

## The rule this table enforces

A `synthetic-substitute` grade must be visible wherever its output is consumed.
Reports state the grade of anything they depend on, and the demo output is not
presentable as a scientific finding. Replacing a substitute with the real thing
is tracked in `docs/tech-debt.md`, not left implicit.
