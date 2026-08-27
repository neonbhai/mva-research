# Golden principles

Mechanically-checkable rules for this repository. Each has an ID. Lint failures,
test failures and review findings **cite the ID**, so a fix instruction is always
traceable to the rule it enforces.

New principles are added only with a decision record. Removing one requires the
same.

---

## Layering

**GP-01 — Dependencies point one way.**
Import layers, low to high:

```
0. models         (pure types; imports nothing from mva except models)
1. errors, determinism, clock, logging_setup   (foundation utilities)
2. config         (may import models + utilities)
3. privacy        (cross-cutting; may import layers 0-2 only)
4. ingestion, annotation, phenotype     (data in)
5. prioritization, mechanisms, interventions  (analysis)
6. evidence       (persistence)
7. reporting      (rendering)
8. cli, pipeline  (composition root)
```

A module may import from the same layer or any lower layer. Never higher.
Enforced by `tests/unit/test_architecture.py`.

**GP-02 — Adapters return typed models, never dicts.**
Any function crossing a boundary into an external tool or data source returns a
Pydantic model from `mva.models`. Parsing happens at the boundary, once.

**GP-03 — The composition root is the only place that knows about everything.**
`mva.pipeline` and `mva.cli` wire stages together. A stage never imports another
stage.

---

## Scientific integrity

**GP-10 — No claim without evidence.**
Any statement rendered into a report, and any score contribution, is backed by an
`EvidenceItem` ID that resolves in the evidence store. `mva.reporting.assertions`
refuses to emit unsourced claims.

**GP-11 — A coordinate without a genome build is invalid.**
Enforced in `GenomicCoordinate`. Cross-build comparison raises; it never coerces.

**GP-12 — Prediction is not observation.**
`AssertionTier.OBSERVED_DATA` may not be paired with an in-silico or inference
evidence type. Enforced in `EvidenceItem`.

**GP-13 — Filtering and ranking are separate.**
Hard filters remove only *invalid or impossible* records (malformed, wrong build,
non-canonical contig). Everything else is flagged and down-ranked, never deleted.
See ADR 0005.

**GP-14 — Absence of information is not negative information.**
`ObservationStatus.NOT_ASSESSED` carries no evidential weight. Only `EXCLUDED`
contributes negative phenotype evidence. Likewise, absence of population-frequency
data is not evidence of rarity.

**GP-15 — Phase is never assumed.**
Two heterozygous variants in one gene are a compound heterozygote only if in
*trans*. `PhaseStatus.UNKNOWN` is preserved end-to-end and penalised, never
silently upgraded.

**GP-16 — Direction of effect is mandatory and signed.**
Every mechanism link and every drug hypothesis carries a signed
`EffectDirection`. A drug whose observed direction disagrees with the required
correction cannot be constructed as an accepted hypothesis. "Unknown" is
tri-state and never counts as agreement.

**GP-17 — Every EvidenceItem states its limitations.**
A non-empty `limitations` field is required. Evidence with no stated limitation is
a smell, not a strength.

**GP-18 — Population frequency carries source, version and population.**
A bare float is rejected by the type system.

**GP-19 — Contradictions are persisted, never discarded.**
Rejected drugs, contradicting evidence and failed candidates remain in the
evidence store with their reasons.

**GP-20 — Mocked adapters are labelled as mocked.**
Every module has a row in `docs/maturity-ledger.md` grading it `real`,
`synthetic-substitute` or `stub`. Reports surface the grade of anything they
depend on. A synthetic substitute is never described as biologically valid.

---

## Determinism & provenance

**GP-30 — Repeat runs are byte-identical.**
No wall-clock timestamps, RNG, dict-order or set-order dependence in artifact
content. Sorting is explicit and total. Timestamps come from the injected clock.

**GP-31 — Every artifact has provenance.**
Written artifacts register an `ArtifactProvenance` with content hash, upstream
IDs, producing stage and sensitivity classification.

**GP-32 — Weights are configuration, not code.**
Changing a scoring weight requires: a decision record, a test, and an explicit
before/after comparison. Golden expectations are never silently re-baselined.

---

## Privacy

**GP-40 — Patient data never enters the repository.**
The workspace is an external path. `config` rejects a workspace resolving inside
the repo tree, and rejects cloud-synced locations.

**GP-41 — The scanner never echoes what it finds.**
Privacy tooling reports path, line, span length and rule ID. Never the matched
bytes.

**GP-42 — Logs cannot carry genomic records.**
A redaction filter is installed on every handler, plus a record factory.

**GP-43 — Public export is allowlist-gated and re-scanned.**
Classification is a claim; the post-render scan is the verification. Both must
pass.

**GP-44 — The sensitive path is deliberately illegible to agents.**
Introspection commands operate on synthetic runs and aggregate counts only. See
ADR 0008.
