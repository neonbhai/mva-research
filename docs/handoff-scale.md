# Handoff: making the pipeline survive a whole-genome callset

> **APPLIED 2026-08-28.** §1.1-§1.7 are in `src/mva/orchestrator.py` and §2's
> `SelectionThresholds` is in `src/mva/config.py` (`CaseConfig.selection`), with
> `config/default.yaml` and its hash lock updated in the same commit.
> `as_payload()` moved onto the Pydantic model rather than being replaced by
> `model_dump`, so `SelectionStream.report()` was untouched and its
> `dict[str, bool | float]` typing survived. §3's "set `selection.enabled: false`
> in `config/synthetic-case.yaml`" was taken: the twelve-record fixture is the
> GP-13 demonstration and enabling selection deletes the two common variants
> `tests/golden/test_golden_case.py` asserts survive (ADR 0027 §7). The demo
> therefore still emits 18 candidates and ten submission rows. The rest of this
> document is kept as the record of what was measured and why.

Everything below is implemented and tested **except** the `orchestrator.py` and
`config.py` edits, which are owned by other agents and are specified here rather
than applied. Nothing in this document is a proposal: the code it calls exists,
the tests pass, and the numbers are measured.

Three blockers were fixed, and a fourth was found on the way.

| # | Blocker | Fix | Where |
|---|---|---|---|
| 1 | `annotate_variants` materialises the callset twice | `iter_annotated()` streams; `annotate_variants` is now a drain of it | `src/mva/annotation/service.py` |
| 2 | `EvidenceLedger` cannot hold ~7 M items | spills to a compressed SQLite file past a threshold | `src/mva/evidence/ledger.py`, `spill.py` |
| 3 | no rarity/coding selection before pairing | a named, counted selection stage | `src/mva/prioritization/selection.py`, ADR 0019 |
| 4 | *(found)* `apply_hard_filters` / `apply_soft_flags` also materialise | `iter_hard_filtered()`, `iter_soft_flagged()` | `src/mva/prioritization/filters.py` |
| 4b | *(found)* an artifact written mid-chain cannot be pulled | `RunContext.open_json_rows_artifact()` push sink | `src/mva/pipeline.py` |

---

## 1. The call-site diff for `orchestrator.py`

`orchestrator.py` was not touched. This is the change, against the file as of
2026-08-28. It assumes the streaming ingest handed over in `docs/scale-report.md`
§8.7 is applied first; if it is not, substitute `qc.variants` for `qc_stream`
below and the annotate half still works (it just keeps the ingest copy alive).

### 1.1 Imports

```diff
+from collections.abc import Iterator
@@
-from mva.annotation.service import annotate_variants
+from mva.annotation.service import iter_annotated
@@
-from mva.prioritization.filters import apply_hard_filters, apply_soft_flags
+from mva.prioritization.filters import iter_hard_filtered, iter_soft_flagged
+from mva.prioritization.selection import iter_selected
@@
+from mva.models.variant import VariantRecord
```

`apply_hard_filters` and `apply_soft_flags` remain exported and unchanged; only
this call site stops using them.

### 1.2 The ledger gets somewhere to spill

```diff
-    ledger = EvidenceLedger(run_id=context.run_id)
+    # The spill directory is inside the workspace: evidence subjects are variant
+    # IDs, so the file carries proband coordinates and is SENSITIVE (GP-40).
+    # `close()` in `_finish` unlinks it. Without a spill_dir the ledger keeps its
+    # pre-existing all-in-memory behaviour, which is what every test gets.
+    ledger = EvidenceLedger(run_id=context.run_id, spill_dir=workspace.tmp_dir)
```

### 1.3 The annotate stage, and selection with it

Replace the whole `# ---- annotate` block (currently lines 268-281) with:

```python
    # ---------------------------------------------------------------- annotate
    adapters = load_default_adapters(knowledge_root, manifest_path)
    annotated = iter_annotated(qc.variants, adapters=adapters, clock=context.clock)

    # ONE pass drives everything: annotation feeds the ledger and the artifact,
    # then the hard filter, the soft flags and selection. Nothing but the selected
    # variants accumulates, and that is a few hundred records rather than 4.5 M.
    #
    # The artifact is written from the point BEFORE hard filtering, so
    # `variants/annotated.json` still holds every annotated record exactly as it
    # does today. That is why this uses the PUSH sink: a puller cannot both write
    # every record and hand a filtered subset onward in the same pass.
    with context.open_json_rows_artifact(
        "variants/annotated.json",
        kind=ArtifactKind.ANNOTATED_VARIANTS,
        stage="annotate",
        upstream=[normalised_art.artifact_id],
    ) as sink:

        def _recorded() -> Iterator[VariantRecord]:
            for item in annotated:
                ledger.extend(item.evidence)
                sink.write(item.variant.model_dump(mode="json"))
                yield item.variant

        filtered = iter_hard_filtered(_recorded(), expected_build=config.genome_build)
        flagged = iter_soft_flagged(
            filtered, frequency=config.frequency, quality=config.quality
        )
        selection = iter_selected(
            flagged,
            frequency=config.frequency,
            thresholds=config.selection,
            clock=context.clock,
        )
        selected = list(selection)

    annotated_art = sink.provenance
    assert annotated_art is not None  # set by the context manager on exit
    context.warnings.extend(annotated.warnings())

    # ---------------------------------------------------------------- selection
    selection_report = selection.report()
    ledger.extend(selection.evidence())
    context.warnings.extend(selection_report.warnings)
    context.write_json_artifact(
        "selection/selection_report.json",
        {
            **selection_report.as_payload(),
            "hard_filter": filtered.counts(),
        },
        kind=ArtifactKind.SELECTION_REPORT,
        stage="select",
        upstream=[annotated_art.artifact_id],
        row_count=selection_report.input_count,
    )

    if not _should_run("prioritise", stop_after):
        return _finish(context, repo_root=repo_root, ledger=ledger, ranked=[])
```

Until `CaseConfig.selection` exists (§2), pass
`thresholds=SelectionThresholds()` and add
`from mva.prioritization.selection import SelectionThresholds, iter_selected` —
the defaults are the same values §2 specifies, so the swap is one line later.

Note `annotated.warnings()` is a **method** on the stream, where
`AnnotationResult.warnings` was an attribute. `selection.report()` and
`selection.evidence()` are only meaningful after `list(selection)` has drained the
chain, which the block above guarantees.

### 1.4 The prioritise stage loses its own filtering

```diff
     # ---------------------------------------------------------------- prioritise
-    filtered = apply_hard_filters(annotation.variants, expected_build=config.genome_build)
-    flagged = apply_soft_flags(
-        filtered.retained, frequency=config.frequency, quality=config.quality
-    )
     pairing = generate_pair_candidates(
-        flagged,
+        selected,
         max_pairs_per_gene=config.max_pairs_per_gene,
         frequency=config.frequency,
     )
```

Hard filtering and soft flagging now happen inside the annotate pass, above; the
records in `selected` have already been through both, carrying identical flags.

### 1.5 `STAGES`

```diff
 STAGES: tuple[str, ...] = (
     "validate",
     "ingest",
     "annotate",
+    "select",
     "prioritise",
     "mechanism",
     "drugs",
     "report",
 )
```

`--stop-after select` then stops after the selection report is written. The
`_should_run("prioritise", ...)` guard shown in §1.3 replaces the one currently
sitting between annotate and phenotype.

`mva.cli` needs no change: `run <stage>` takes a free-form string and
`_should_run` validates it against `STAGES`. The Snakemake DAG does: add a
`SELECTION_REPORT = f"{RUN_DIR}/selection/selection_report.json"` path and a
`workflow/rules/select.smk` shaped exactly like `annotate.smk` — input
`ANNOTATED_VARIANTS`, output `SELECTION_REPORT`, `shell: mva_run("select")` — and
make `rule prioritise` depend on it. `just workflow-check` is part of the gate, so
a `STAGES` entry with no rule behind it will be noticed there rather than later.

### 1.6 Persisting a spilled ledger

```diff
     with EvidenceStore(db_path) as store:
         store.initialise()
-        store.write_evidence(ledger.items())
+        # `items()` materialises; a spilled ledger holds more items than the
+        # process has memory for, which is why it spilled.
+        store.write_evidence_stream(ledger.iter_items())
```

and, at the end of `_finish`, after `_persist_evidence` has run:

```diff
     _write_privacy_audit(context, repo_root=repo_root)
+    # The spill file carries proband coordinates and the run is over.
+    ledger.close()
```

`close()` is a no-op on a ledger that never spilled, so it is safe
unconditionally, and `len(ledger)` keeps working afterwards — `PipelineResult`
reads `evidence_count=len(ledger)` several lines later and does not need
reordering. Reading the ledger's *contents* after close raises `EvidenceError`
rather than returning an empty ledger, which would read as a run that produced no
evidence.

### 1.7 The one test that must change with it

`tests/unit/test_privacy.py::test_the_offline_profile_is_armed_while_the_annotate_stage_executes`
monkeypatches `orchestrator.annotate_variants`. That name will no longer be bound.
Patch `orchestrator.iter_annotated` instead — the probe must return the stream
object unchanged, so:

```python
    real = orchestrator.iter_annotated

    def probe(*args: object, **kwargs: object) -> object:
        observed.append(netguard.is_armed())
        return real(*args, **kwargs)  # type: ignore[arg-type]

    orchestrator.iter_annotated = probe  # type: ignore[assignment]
```

The assertion `observed == [True]` still holds: `iter_annotated` is called once,
inside the armed block. (It is called before the stream is drained, which is the
same instant the old call happened.)

---

## 2. The config fields, for `config.py`

Add one model and one field on `CaseConfig`. Until they exist,
`mva.prioritization.selection.SelectionThresholds` is a frozen dataclass owned by
that module with the identical field set and defaults, and `iter_selected`
defaults to it — so the stage works today and the move is a type swap, not a
rewrite.

```python
class SelectionThresholds(StrictModel):
    """Cut-points for the candidate-selection stage (ADR 0019).

    Separate from `FrequencyThresholds`, which holds *ranking* cut-points. This
    stage is the only one entitled to delete a valid record, so its thresholds are
    stated apart from the ones that merely down-rank (GP-13).
    """

    enabled: bool = Field(
        default=True,
        description=(
            "Off retains every variant with reason 'retained_selection_disabled'. "
            "The stage still counts and still reports, so 'what did selection cost "
            "me' is answerable without rewiring."
        ),
    )
    max_population_frequency: float = Field(
        default=0.02,
        gt=0.0,
        le=1.0,
        description=(
            "Selection cut-point on the maximum AF across adequately powered "
            "populations. Deliberately LOOSER than "
            "frequency.max_plausible_recessive (0.01), so nothing this stage "
            "deletes was still rankable above zero on rarity. Inverting that "
            "relationship is reported as a run warning, not silently honoured."
        ),
    )
    retain_unknown_frequency: bool = Field(
        default=True,
        description=(
            "GP-14. A variant with no usable frequency observation is UNKNOWN, "
            "not rare and not common, and is kept. False is the setting that "
            "loses a novel pathogenic variant, preferentially for ancestries the "
            "reference cohorts under-sample. See ADR 0019."
        ),
    )
    retain_unassessed_impact: bool = Field(
        default=True,
        description=(
            "ADR 0016. `impact is None` means NOT ASSESSED, not MODIFIER. A MANE "
            "interval join produces exactly this shape for every variant it places."
        ),
    )
    retain_pathogenic_clinical: bool = Field(
        default=True,
        description="A curated pathogenic assertion overrides both other gates.",
    )
    drop_without_gene_assignment: bool = Field(
        default=True,
        description=(
            "Drop variants with no consequence annotation at all. Gene-scoped "
            "pairing already ignores them; dropping them here makes the count "
            "visible. ~55% of a whole-genome callset by volume."
        ),
    )
    min_splice_ai_delta: float = Field(
        default=0.2,
        ge=0.0,
        le=1.0,
        description=(
            "A SpliceAI delta at or above this retains the variant whatever its "
            "terms say. SpliceAI's own high-recall cut-point."
        ),
    )


class CaseConfig(StrictModel):
    ...
    selection: SelectionThresholds = Field(default_factory=SelectionThresholds)
```

and in `config/default.yaml`:

```yaml
selection:
  enabled: true
  max_population_frequency: 0.02
  retain_unknown_frequency: true
  retain_unassessed_impact: true
  retain_pathogenic_clinical: true
  drop_without_gene_assignment: true
  min_splice_ai_delta: 0.2
```

**`config/default.yaml` is hash-locked** (`tests/golden/test_locked_files.py`).
Adding this block changes its sha256, so the lock must be updated in the same
commit — that is the gate working, and ADR 0019 is the decision record it wants.

Two call sites need the type swap once the field lands:
`orchestrator.py` passes `thresholds=config.selection`, and
`selection.iter_selected`'s default parameter becomes `config.SelectionThresholds`.
`SelectionThresholds.as_payload()` has no equivalent on a Pydantic model — replace
its use in `SelectionStream.report()` with `model_dump(mode="json")`, which emits
the same keys.

---

## 3. What this does to the synthetic demo

Measured, not projected: `apply_hard_filters` → `apply_soft_flags` →
`select_variants` → `generate_pair_candidates` → `score_pair` → `rank_pairs`
over `tests/fixtures/synthetic/synthetic_case.vcf` with the real local adapters.

```
12 records in, 12 survive hard filtering
selection: 8 retained, 4 dropped
    dropped_common_in_population    3
    dropped_not_coding_or_splice    1
    retained_coding_or_splice       8
    (note) frequency_unknown        1   <- retained, GP-14
candidates: 18 -> 9
```

**The ranking is unchanged.** Before and after selection, ranks 1-5 are the same
pairs with the same composite scores to six decimal places:

```
1  SYNTHKIN1  0.877899  chr15:40200000:C:T + chr15:40210500:G:A
2  SYNTHKIN1  0.823133  chr15:40200000:C:T                       (single-variant)
3  SYNTHKIN1  0.820499  chr15:40210500:G:A                       (single-variant)
4  SYNTHMET2  0.667470  chr7:55000000:C:G + chr7:55001000:A:T
5  SYNTHKIN1  0.617133  chr15:40200000:C:T + chr15:40211000:G:GAT
```

Rank 1 is exactly `tests/golden/expected_ranking.tsv`, so **the golden
expectation still passes** and needs no re-baselining.

One consequence that is not a test failure but is a behaviour change: the demo
now has 9 candidates rather than 18, so `submission/track1_submission.csv` carries
**9 data rows instead of 10**. No test asserts a row count from the demo (the
row-count tests in `tests/integration/test_track1_format.py` build their own
pairs), but `docs/next-actions.md` says "fill all ten rows: rows below the answer
cannot lower either metric". On a real 4.5 M callset there will be far more than
ten candidates, so this is an artifact of a 12-record fixture — but if you want
the demo to keep emitting ten rows, set `selection.enabled: false` in
`config/synthetic-case.yaml` and the stage will report what it *would* have done
without acting on it.

---

## 4. The APIs, in one place

### `mva.annotation.service`

```python
iter_annotated(
    variants: Iterable[VariantRecord], *, adapters: AdapterSet, clock: Clock,
    batch_size: int | None = DEFAULT_ANNOTATION_BATCH_SIZE,   # 5_000
) -> AnnotationStream
```
`AnnotationStream` iterates `AnnotatedVariant(variant, evidence)` and offers
`.coverage() -> dict[str, float]`, `.warnings() -> tuple[str, ...]`, `.drained`.
Single-pass; re-iterating raises `AnnotationError`. `batch_size=None` means one
whole-input batch, which is exactly the pre-streaming behaviour and is what
`annotate_variants` now uses.

`annotate_variants(...) -> AnnotationResult` is unchanged in signature and
output. It is a drain of the stream, so the two cannot drift.

### `mva.evidence.ledger`

```python
EvidenceLedger(
    *, run_id: str,
    spill_dir: Path | None = None,          # None = today's behaviour, in memory
    spill_threshold: int = DEFAULT_SPILL_THRESHOLD,   # 50_000
    flush_batch: int = DEFAULT_FLUSH_BATCH,           # 10_000
)
```
New: `.iter_items()`, `.spilled`, `.close()`, context-manager support.
`items()`, `for_subject()`, `contradictions()`, `get()`, `__len__`,
`__contains__`, `__iter__` are unchanged in behaviour — `for_subject` and
`contradictions` additionally stopped sorting the whole ledger to answer about
one subject.

### `mva.evidence.store`

```python
EvidenceStore.write_evidence_stream(
    items: Iterable[EvidenceItem], *, batch_size: int = EVIDENCE_WRITE_BATCH  # 50_000
) -> int
```
Byte-identical Parquet export to `write_evidence` over the same items; a test
asserts the exported files match byte for byte.

### `mva.prioritization.selection`

```python
iter_selected(
    variants: Iterable[VariantRecord], *, frequency: FrequencyThresholds,
    thresholds: SelectionThresholds | None = None, clock: Clock,
) -> SelectionStream
select_variants(...) -> SelectionResult        # the fixture-scale API
```
`SelectionStream` iterates the **retained** records; `.decisions()` iterates every
verdict instead (one or the other, not both). `.report() -> SelectionReport` and
`.evidence() -> tuple[EvidenceItem, ...]` after draining.

### `mva.prioritization.filters`

```python
iter_hard_filtered(variants, *, expected_build) -> HardFilterStream   # .counts()
iter_soft_flagged(variants, *, frequency, quality) -> Iterator[VariantRecord]
```
Same records, same flags, same counts as the batch functions. The one thing the
stream does **not** keep is `FilterResult.removed`, the `(variant_id, reason)`
list of everything discarded: at WGS scale that is millions of proband
coordinates held in order to be counted, and the counts are what a report can
state anyway (GP-41).

### `mva.pipeline`

```python
with context.open_json_rows_artifact(rel, kind=..., stage=..., upstream=...) as sink:
    sink.write(row)          # bytes identical to write_json_rows_artifact
artifact = context.artifacts[-1]
```

---

## 5. Things to know that I did not fix

- **`_persist_evidence` writes every ledger item to DuckDB and then to Parquet,
  and `verify_determinism` hashes the Parquet.** With the ledger holding ~6.9 M
  items, that is a large DuckDB file and a large hash on every run, twice for
  `just verify`'s determinism check. `write_evidence_stream` makes it *possible*;
  it does not make it cheap. If the run time matters more than completeness, the
  composition to consider is pushing QC and annotation evidence into the ledger
  only for variants that survive selection — but that weakens GP-19 further than
  ADR 0019 already concedes, so it is a decision, not an optimisation.
- **Two concurrent runs of the same case in the same workspace would collide on
  the spill filename** (`evidence-ledger-<run_id>.sqlite`). They already collide
  on the run directory, since `run_id` is derived from the config, so this is a
  pre-existing property rather than a new one — but the spill makes it a data
  hazard rather than an overwrite. A pid or a lock file would close it.
- **`generate_pair_candidates` still takes a `Sequence`.** That is correct: after
  selection it receives hundreds of records, not millions. It is only safe
  *because* selection runs in front of it, which is the point of ADR 0019.
- **`selection` depends on annotation having run.** With no consequence adapter
  bound, every variant has zero consequences and `drop_without_gene_assignment`
  deletes the entire callset. The report makes that unmistakable
  (`dropped_no_gene_assignment` equal to the input count) rather than silent, but
  it is a real coupling and `docs/next-actions.md` item 1 has to land first.
- **The scale phantom is not biologically valid.** `tools/scale/stage_harness.py`
  fabricates every coordinate, allele, gene, consequence and frequency from a
  seeded hash. It measures bytes and objects per second and says nothing about
  variants, genes or disease (GP-20).
