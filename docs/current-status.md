# Current status

Updated 2026-08-28. Phase 7: real data connected.

## What exists

**The proband data is here.** Access to `SageBio/mva-hackathon-2026-data` was
granted and the case files were fetched into the encrypted volume:

| File | Size | sha256 (first 16) |
|---|---|---|
| `WGS_SPECIMEN_FLOWCELL.vcf.gz` | 305 MB | `REDACTED-HASH` |
| `WGS_SPECIMEN_FLOWCELL.vcf.gz.tbi` | 2.2 MB | `REDACTED-HASH` |
| `Challenge_Clinical_Phenotype_1.docx` | 17 KB | — |

They live only in `/Volumes/MVACASE/raw`, mode 700. The 79 GB of FASTQs was
deliberately not fetched (`docs/resource-acquisition-assessment.md`).

`/Volumes/MVACASE` is an AES-256 APFS sparse bundle, Spotlight indexing off,
`nobrowse`, key at `~/private/REDACTED-KEY-PATH` (mode 600). FileVault is on; no Time
Machine destination is configured. The bundle exists for the **30-day deletion
obligation** — on APFS only cryptographic erasure proves deletion.

**Real reference data** at `~/Contri/bio-hackathon/mva-resources` (outside the
repo; `knowledge/manifests/resources.yaml` is the committed record): ClinVar
GRCh38 `2026-08-22` (185 MB), HPO `2026-06-23` (66 MB), gnomAD v4.1 constraint
(92 MB), DDG2P, ClinGen, MANE v1.5, SnpEff GRCh38.115, gnomAD v4.1 exome sites
(184.8 GB, downloading), GRCh38 reference FASTA (3.3 GB, downloading).

## What works

- **ClinVar adapter, verified against the real release.** Correct significance,
  star ratings and review status; **absent variants are omitted, not called
  benign** (GP-14). 2,989 P/LP variants across the CIN gene panel.
- **Real HPO semantics** — ontology DAG, IC from the actual annotation corpus
  (not graph depth), Lin/Resnik/Jiang-Conrath, Köhler symmetric best-match.
  Directional propagation: OBSERVED up, EXCLUDED down. Determinism proven across
  three `PYTHONHASHSEED` values.
- **Contig conversion, the highest-stakes detail in the project.** The proband
  VCF uses **bare contigs**; the scorer compares chrom strings raw. Verified:
  `15` → `chr15`, `MT` → `chrM`, decoys and scaffolds rejected. The renderer
  forces UCSC and `validate_submission` re-checks the rendered bytes.
- **EPCR ties eliminated** — ten strictly-decreasing values, minimum gap 0.0100,
  order-preserving by construction. Causal pair still row 1 at REDACTED-EPCR.
- **Pairing cap no longer truncates by coordinate** — it orders by plausibility
  and reports when it fires (ADR 0013).
- **Track 2 runs on real BUB1B literature, not the synthetic tables.**
  `knowledge/literature/bub1b/` holds a 12-node, 13-link cited mechanism chain and
  13 real agents with signed directions (ADR 0022, ADR 0025). Triage output:
  1 accepted, 12 rejected, 89 evidence rows; the MPS1 inhibitor ranks **last**
  with `WRONG_DIRECTION`, and the highest-scoring agent is rejected anyway. Case:
  `docs/track2-hypothesis.md`. Asserted by
  `tests/unit/test_track2_bub1b.py` (22 tests), including a both-directions
  citation check over `SOURCES.md`. **Deliberately not wired into the composition
  root** — see ADR 0022.
- 756 tests passing.

## What is incomplete

- **Composition-root wiring.** The real adapters (ClinVar, gnomAD, SnpEff, MANE
  gene intervals) and HPO semantics exist but are **not yet bound in
  `orchestrator.py`**. Until they are, the pipeline still runs on synthetic
  tables. This is the single largest remaining gap.
- **`generate_pairs` yields nothing without gene assignment.**
  `VariantRecord.gene_symbols` is derived purely from `consequences`, so with no
  annotator the pipeline emits **zero candidate pairs**. Two fixes in flight.
- **Indel left-alignment.** a substantial fraction of scanned records are indel-bearing and no
  reference FASTA is configured, so minimal representation is not reached. Indels
  may silently fail to join gnomAD/ClinVar — which looks exactly like "rare".
- **gnomAD constraint has no chrX/chrY** (verified byte-complete download). An
  X-linked candidate gets no pLI/LOEUF, and absence must not read as unconstrained.
- Exome-only frequencies for a WGS proband: exact per-ancestry AF is available
  for well under 1% of the callset. A regulatory causal variant is out of reach.
- TD-01 partially paid: clinical slot done, frequency and consequence in flight.

## Commands

```bash
just verify                       # the blocking gate
just demo                         # synthetic end-to-end
just privacy-audit                # standalone privacy scan
uv run python tools/inspect/vcf_schema.py <vcf>   # schema metadata only, never records
```

Real-case environment:
```bash
export MVA_WORKSPACE=/Volumes/MVACASE/case01
export TMPDIR=/Volumes/MVACASE/tmp
```

## Blockers

1. **Composition-root wiring** — nothing real reaches the pipeline until
   `orchestrator.py` binds the new adapters. Blocked only on concurrent edits.
2. **Gene assignment must land** or the submission is empty.
3. **Reference FASTA download** must finish before indels can be normalised to
   minimal representation.

## Known-dangerous details

- The proband VCF's reference is `..._no_chr.fasta` — **bare contigs**. The
  scorer does no normalisation. Emitting `15` scores exactly 0.0 while looking
  correct. Verified handled, but any new emission path must re-check.
- `INFO/AF` in this VCF is GATK's **sample** allele frequency, not a population
  frequency. Reading it as population AF would mark every het as common.
- Single sample, so phase is UNKNOWN for every pair (TD-07 is binding, not
  theoretical).
- Multi-allelics are **not** split in the source (1.46%); normalisation must.
