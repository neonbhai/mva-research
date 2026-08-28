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
| `annotation` | **mixed, by case**: a `synthetic: false` case is **real**; the demo path is a **synthetic-substitute** | Adapter protocols, manifest hash-pinning, multi-transcript preservation, "no data ≠ rare" handling. Since ADR 0027 a non-synthetic case binds real SnpEff 5.4c/GRCh38.115 (pinned from the reviewed installer manifest, never `measure()`), gnomAD v4.1 exomes, ClinVar 2026-08-22 and MANE v1.5, with the GRCh38 no-alt FASTA passed to both joining adapters so repeat-tract indels actually join (ADR 0018). Missing resources are a refusal, never a silent fallback. | **The demo and every test still run on fabricated consequence and frequency tables** with fictional `SYNTH*` genes, deliberately (ADR 0027): a suite that depended on 200 GB of reference data would pass in one place and fail in another. No demo output is biologically valid. SpliceAI is bound in neither branch. The real branch has never been executed end to end against the real callset here — the proband VCF is outside what agents may read (ADR 0008). | TD-01 |
| `phenotype` | **mixed**: ontology and `knowledge/real/gene_phenotype.tsv` are real; the demo path is a **synthetic-substitute** | Four-valued observation logic, gene–phenotype scoring, contradiction handling. Real HPO release, real annotation corpus and information-content similarity (ADR 0017). `knowledge/real/gene_phenotype.tsv` is a real HPO table (221,789 rows, 3,634 genes) whose `association_strength` is real ClinGen/DDG2P clinical validity and whose `hpo_frequency` is HPO's own frequency vocabulary (ADR 0021). | The gene–phenotype table the **demo** is wired to (`knowledge/public/gene_phenotype.tsv`) is still fabricated — fictional `SYNTH*` gene symbols — so no demo output is biologically valid. The real table is loadable and scored against, but is not yet wired into the composition root. `association_strength` is gene-level, not (gene, disease)-level. | TD-20 |
| `prioritization` | real | Hard/soft filter separation, the named candidate-selection stage in front of pairing (ADR 0019), pairing, phase inference, component scoring, ranking. | Weights are uncalibrated heuristics. Phase is UNKNOWN for all pairs without trio data. Selection's coding/splice rule rests on predicted consequence, the weakest link in the annotation chain. | TD-02, TD-07 |
| `mechanisms` | **mixed**: the demo path is a **synthetic-substitute**; the BUB1B chain is **literature-derived** | Chain loading, link grading, inferred-link penalty, target resolution. `knowledge/literature/bub1b/mechanisms.tsv` is a real 12-node, 13-link BUB1B/BubR1 chain, each link graded and cited to a PMID in that directory's `SOURCES.md` (ADR 0022, ADR 0025). | The SYNTHKIN1 mechanism the **demo** is wired to is fictional — shape mirrors real checkpoint biology, content does not. The literature chain is loadable and tested but is deliberately **not** wired into the composition root: a pipeline whose behaviour depends on which gene it is looking at is not a pipeline. | TD-08 |
| `interventions` | **mixed**: the demo path is a **synthetic-substitute**; the BUB1B catalogue is **literature-derived** | Direction checking, tri-state agreement, safety filtering, evidence tiering, rejection recording. `knowledge/literature/bub1b/drug_catalog.tsv` holds 13 real agents with signed, cited directions; the triage produces 1 accepted and 12 rejected (`tests/unit/test_track2_bub1b.py`). | The catalogue the **demo** renders is fabricated. The literature catalogue evaluates real compounds, but it is hand-curated at n=13 and is **not** a systematic screen of any drug database — `observed_direction` is not a column any public drug database carries (`docs/track2-hypothesis.md` §8.3). It scores **one** context, the soma; it is not a tumour-board tool (ASSUMPTION-DRUG-08). | TD-08, TD-25 |
| `evidence` | real | DuckDB schema, Parquet export, idempotent writes, assertion resolution, and a compressed SQLite spill so the ledger survives a whole-genome run's ~7 M items in bounded memory (ADR 0019 handoff, `docs/scale-report.md` §9). | Cross-platform byte-identity untested. The spill is scratch, not an artifact: it is never hashed and is unlinked on close. | TD-05, TD-18, TD-19 |
| `reporting` | real | Assertion checking, tier marking, Track 1 CSV against the verified contract, dossier and Track 2 rendering. | Reports are only as good as the synthetic inputs they render. | — |
| `cli` / `pipeline` | real | Composition root, stage orchestration, determinism verification. | Determinism verification is synthetic-only: `verify determinism` inherits the clock `config.synthetic` picks, so it has only ever run frozen. A real run's scientific content is byte-identical; eleven artifacts differ in recorded time (`docs/handoff-integrity.md` §4). | TD-21 |
| Structural / CNV calling | **absent** | — | Not implemented at all. A material gap for a chromosomal-instability disorder. | TD-03 |
| Repeat-expansion calling | **absent** | — | Not implemented at all. Short-tandem-repeat expansions are neither called nor representable, so a negative result at a repeat locus means "not looked at", not "not present". | TD-13 |
| mtDNA heteroplasmy | **absent** | `chrM` calls are recognised and modelled as `InheritanceModel.MITOCHONDRIAL` instead of being described in nuclear two-copy language. | The heteroplasmy fraction is never measured or represented, and neither is the sampled tissue. Every mitochondrial candidate scores the same flat inheritance value regardless of load. | TD-14 |

## The rule this table enforces

A `synthetic-substitute` grade must be visible wherever its output is consumed.
Reports state the grade of anything they depend on, and the demo output is not
presentable as a scientific finding. Replacing a substitute with the real thing
is tracked in `docs/tech-debt.md`, not left implicit.
