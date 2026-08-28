# Handoff: integrity fixes that land outside one owner's files

Four defects were found by an independent adversarial review, each of the same
shape: **the system states something it has not established.** Three are fixed in
place. The parts that belong in files owned by another agent are written out here
as exact diffs, and the fourth — GP-30's real scope — is a documentation defect
whose whole remedy is to write down what is true.

Nothing in this file is a plan. Every diff below is complete and applies to the
code as it stands; every claim is traced to a file and a line.

| # | Defect | Fixed here | Handed off |
|---|---|---|---|
| 1 | A failing FASTA degraded to trim-only while adapters reported `APPLIED` | `src/mva/alleles.py`, `src/mva/errors.py` | `annotation/gnomad_sites.py`, `annotation/clinvar_vcf.py`, `ingestion/normalise.py` (§1) |
| 2 | Patient coordinates and genotypes in exception messages | `models/variant.py`, `models/pair.py` | `phenotype/loader.py` (§2) |
| 3 | Annotation is not structurally network-isolated | `tests/unit/test_architecture.py`, `docs/privacy-model.md`, TD-06 | Nothing code-side; the residual hole is an OS control (§3) |
| 4 | Real runs are not byte-reproducible | Documentation only | `pipeline.py`, `prioritization/ranking.py`, `interventions/safety.py` (§4) |

---

## 1. The reference-failure lie

### What was wrong

`alleles.py` caught **every** exception from `ReferenceLookup.fetch` and returned
`None`, which the shift loop read as "no base here, stop". Four different
conditions arrived at that one `None`:

| Condition | What it means | What it used to do |
|---|---|---|
| `position < 1` | before the start of a contig | stop — **correct** |
| read past the end of a contig | the chromosome ends here | stop — **correct** |
| accessor raised (`OSError`, missing contig) | the reference is broken | stop — **a lie** |
| empty or non-nucleotide sequence mid-contig | the reference is broken | stop — **a lie** |

Both adapters then reported left-alignment as `APPLIED` on the strength of the
reference object merely being non-`None` (`gnomad_sites.py:1247`,
`clinvar_vcf.py:628`), so a FASTA that raised on every read produced trim-only
join keys underneath a claim of reference-backed canonicalisation.

The cost is not cosmetic. On the real gnomAD v4.1 exomes chr21 shard, 1,029 of
1,382 indel records in a 520 kb exonic window (74.5%) sit in a repeat tract and
therefore have another legal spelling; without a working reference **0 of those
1,029** right-shifted spellings join, against 989 with one, and 30 of the
recovered records are alleles gnomAD calls common (ADR 0018). A failed join
presents as "no frequency data", which the ranker scores as novel and ultra-rare
— the strongest promoting signal it has. So the lie manufactures false positives
and deletes true ClinVar assertions at the same time, silently, in both
directions.

### What was fixed in `src/mva/alleles.py`

* `_read_one_base` replaces `_fetch_base` and **classifies** instead of
  collapsing: `None` only for a position the reference genuinely does not cover;
  `mva.errors.ReferenceUnusableError` for a raise, a missing contig, an empty
  read where a base must exist, or a non-nucleotide.
* `_required_base` is the leftward reader. Every base the left shift asks for is
  at `position - 1` with `position > 1`, i.e. inside a contig the record claims to
  sit on, so absence there is breakage.
* `ReferenceStatus` (`NOT_SUPPLIED` / `USABLE` / `UNUSABLE`) is carried on
  `CanonicalAllele`, with `CanonicalAllele.left_alignment_proven` as the property
  an adapter must branch on. This is `CanonicalAllele.operations`' discipline
  applied to the failure case: `operations` answers "did the position move", which
  is ambiguous on its own, because a position that did not move may be already
  left-most *or* may have hit an unreadable reference.
* `canonicalise_allele` catches `ReferenceUnusableError` at exactly one
  documented boundary and keeps the **trimmed** form, not the partially shifted
  one — an allele abandoned half way through a repeat tract joins against neither
  the source nor another run of this pipeline.
* `rightmost_equivalent_bound` returns a `QueryBound` (position + status).
  `rightmost_equivalent_position` remains as an `int`-returning shim for the two
  adapters that call it; its docstring says plainly that a caller of it cannot
  tell a proven bound from where a broken reference stopped.
* The raised message is PRIV-09-safe: `error_token` handle, no contig, no
  position, and `from None` so a backend's own coordinate-bearing message is not
  chained through into the traceback.

Tests: `tests/unit/test_allele_reference_integrity.py` (20 cases). ADR 0026 records
the decision and the limits that remain.

### 1a. `src/mva/annotation/gnomad_sites.py` — REQUIRED

The adapter must stop deriving `APPLIED` from `self._reference is not None`.

```diff
@@ class GnomadSitesFrequencyAdapter.__init__
         self._reference = reference
+        # Per-run count of alleles whose canonicalisation could not read the
+        # reference. Not a boolean: the report states how many records are
+        # affected, and "one unreadable base" and "the FASTA is gone" are
+        # different operator problems.
+        self._unusable_reference_alleles = 0

@@ def representation_status(self) -> LeftAlignmentStatus:
-        ``APPLIED`` when a reference was supplied, ``UNAVAILABLE_NO_REFERENCE``
-        otherwise. Trimming is unconditional either way, so the non-minimal class
+        ``APPLIED`` only when a reference was supplied **and every read the rule
+        needed from it succeeded**. A reference object that is merely non-``None``
+        proves nothing: a FASTA that raises on every read yields trim-only keys,
+        and reporting those as reference-backed is a provenance lie no downstream
+        consumer can detect. Trimming is unconditional either way, so the non-minimal class
         of mismatch joins in both states; what the degraded state costs is the
         repeat-tract class, where gnomAD and the caller place one insertion at
         different positions. Typed rather than logged so a report can state the
         limitation instead of a reader having to infer it from a missing key —
         which is precisely the inference GP-14 forbids.
         """
         if self._reference is None:
             return LeftAlignmentStatus.UNAVAILABLE_NO_REFERENCE
+        if self._unusable_reference_alleles:
+            return LeftAlignmentStatus.INCOMPLETE_REFERENCE_UNUSABLE
         return LeftAlignmentStatus.APPLIED

@@ def representation_limitation(self) -> str | None:
     def representation_limitation(self) -> str | None:
         """One sentence for a report footer, or ``None`` when nothing is degraded."""
-        if self._reference is not None:
+        if self._reference is not None and not self._unusable_reference_alleles:
             return None
+        if self._reference is not None:
+            return (
+                f"A reference FASTA was supplied but could not be read for "
+                f"{self._unusable_reference_alleles} gnomAD lookup(s), which were "
+                "canonicalised by trimming only. Those variants may fail to join "
+                "against gnomAD's left-aligned alleles and would then be reported as "
+                "having no frequency data — absence of information, not evidence of "
+                "rarity. Check that the FASTA is the same assembly and patch release "
+                "as the callset, and that its .fai index matches the file."
+            )
         return (

@@ def canonicalise(self, contig, position, ref, alt) -> CanonicalAllele:
-        return canonicalise_allele(
+        canonical = canonicalise_allele(
             contig=contig,
             position=position,
             ref=ref,
             alt=alt,
             reference=self._reference,
         )
+        if canonical.reference_status is ReferenceStatus.UNUSABLE:
+            self._unusable_reference_alleles += 1
+        return canonical
```

and the import at the top of the module:

```diff
 from mva.alleles import (
     CanonicalAllele,
     LeftAlignmentStatus,
     ReferenceLookup,
+    ReferenceStatus,
     canonicalise_allele,
-    rightmost_equivalent_position,
+    rightmost_equivalent_bound,
 )
```

`_parse_variant_id` (`gnomad_sites.py:725-745`) should take the bound rather than
the bare int, so an under-reaching query window is observable rather than
inferred:

```diff
     search_end = canonical.position
     if reference is not None:
-        search_end = max(
-            search_end,
-            rightmost_equivalent_position(
-                contig=canonical_contig,
-                position=canonical.position,
-                ref=canonical.ref,
-                alt=canonical.alt,
-                reference=reference,
-            ),
+        # `.proven` is False when the reference could not be read to the right of
+        # the record. The window is then the furthest point the search reached,
+        # which may be short — but the key is trim-only for the same reason, so
+        # the limitation is the one already reported by representation_status
+        # rather than a second, undeclared one.
+        bound = rightmost_equivalent_bound(
+            contig=canonical_contig,
+            position=canonical.position,
+            ref=canonical.ref,
+            alt=canonical.alt,
+            reference=reference,
         )
+        search_end = max(search_end, bound.position)
```

**Correction, found while applying this diff: the counter above is not enough on
its own.** `gnomad_sites.py` reaches `canonicalise_allele` by two paths, and the
diff only counts one. The release side goes through `self.canonicalise`
(`_record_frequencies`), but gnomAD's release is already left-aligned, so
`_left_shift` exits on its first comparison and *no reference base is ever read*
for the overwhelming majority of records — the counter stays 0. The reads that
actually fail are on the **query** side, in the module-level `_parse_variant_id`,
which does not touch the adapter. Applied as published, the test below fails with
`APPLIED`. The fix is to carry `reference_status` on `_Query` (set from the
canonical allele, and from `QueryBound.proven` for the fetch window) and to add the
count in `lookup_partial`, which is the single point both `frequencies` and
`lookup_partial` pass through. `clinvar_vcf.py` does not have this problem: its
`_parse_query` is a method and already calls `self.canonicalise`, so there the
published diff is sufficient — only the `QueryBound` half needed a guard against
double-counting one allele that fails on both the key and the bound.

**Do not make `_parse_variant_id` raise.** `tests/unit/test_gnomad_adapter.py::test_a_reference_that_raises_cannot_leak_a_coordinate_or_break_the_lookup`
pins the opposite behaviour deliberately, and it is right to: a lookup that dies
on an unreadable base is worse than one that degrades and says so. That is why
`rightmost_equivalent_bound` reports rather than raises.

Test to add alongside (`tests/unit/test_gnomad_adapter.py`):

```python
def test_a_broken_reference_is_not_reported_as_applied() -> None:
    """The adapter must not claim reference-backed canonicalisation it did not get.

    Reproduced before the fix: `representation_status` was derived from
    `self._reference is not None`, so a FASTA raising on every read reported
    APPLIED over trim-only join keys.
    """

    class _Broken:
        def fetch(self, contig: str, start: int, end: int) -> str:
            raise OSError("handle closed")

    with open_adapter(FIXTURE, reference=_Broken()) as instance:
        instance.frequencies([REPEAT_INSERTION_SHIFTED_ONE])
        assert instance.representation_status is (
            LeftAlignmentStatus.INCOMPLETE_REFERENCE_UNUSABLE
        )
        assert instance.representation_limitation is not None
```

### 1b. `src/mva/annotation/clinvar_vcf.py` — REQUIRED

Structurally identical. `representation_status` is at `clinvar_vcf.py:617`,
`representation_limitation` immediately below it, `canonicalise` at
`clinvar_vcf.py:766`, and the `rightmost_equivalent_position` call at
`clinvar_vcf.py:825-834`.

```diff
@@ __init__
         self._reference = reference
+        self._unusable_reference_alleles = 0

@@ def representation_status(self) -> LeftAlignmentStatus:
         if self._reference is None:
             return LeftAlignmentStatus.UNAVAILABLE_NO_REFERENCE
+        if self._unusable_reference_alleles:
+            return LeftAlignmentStatus.INCOMPLETE_REFERENCE_UNUSABLE
         return LeftAlignmentStatus.APPLIED

@@ def representation_limitation(self) -> str | None:
-        if self._reference is not None:
+        if self._reference is not None and not self._unusable_reference_alleles:
             return None
+        if self._reference is not None:
+            return (
+                f"A reference FASTA was supplied but could not be read for "
+                f"{self._unusable_reference_alleles} ClinVar lookup(s), which were "
+                "canonicalised by trimming only. Those variants may fail to join "
+                "against ClinVar's left-aligned alleles and would then be reported "
+                "as having no ClinVar record — absence of information, not evidence "
+                "of benignity."
+            )
         return (

@@ def canonicalise(self, contig, position, ref, alt) -> CanonicalAllele:
-        return canonicalise_allele(
+        canonical = canonicalise_allele(
             contig=contig,
             position=position,
             ref=ref,
             alt=alt,
             reference=self._reference,
         )
+        if canonical.reference_status is ReferenceStatus.UNUSABLE:
+            self._unusable_reference_alleles += 1
+        return canonical
```

plus the same import change and the same `rightmost_equivalent_bound` migration
in the `_Query` builder.

Note for whoever applies this: the counter belongs on the adapter, not on the
call, because `representation_status` is a *run-level* claim. A per-lookup return
value would let one degraded record be reported and the next one not.

### 1c. `src/mva/ingestion/normalise.py` — RECOMMENDED

Ingestion is mostly protected already: `_reference_matches` (`normalise.py:389`)
reads the record's own REF span first and returns `None` on any exception, which
sets `aligner = None`, raises `WARN_REFERENCE_LOOKUP_FAILED` and marks the record
`reference_consulted=False`. So a FASTA that is broken *everywhere* is caught.

The residual gap is narrow but real: if the REF span reads fine and a base the
*shift* needs does not, `reference_consulted` stays `True`, `unaligned` stays 0,
and `summarise_left_alignment` returns `APPLIED` for the batch.

```diff
@@ def _normalise_one(...)
-    return _Outcome(
-        record=trim_and_left_align(current, aligner),
-        reference_consulted=aligner is not None,
-    )
+    normalised, status = trim_and_left_align_with_status(current, aligner)
+    if status is ReferenceStatus.UNUSABLE:
+        # The REF span read cleanly but a base the shift needed did not. Without
+        # this the batch reports APPLIED for a record that was trimmed only.
+        warnings[WARN_REFERENCE_LOOKUP_FAILED] += 1
+    return _Outcome(
+        record=normalised,
+        reference_consulted=aligner is not None and status is not ReferenceStatus.UNUSABLE,
+    )

@@ after trim_and_left_align
+def trim_and_left_align_with_status(
+    record: VariantRecord, reference: ReferenceLookup | None
+) -> tuple[VariantRecord, ReferenceStatus]:
+    """:func:`trim_and_left_align`, plus how far the reference could be trusted.
+
+    Separate from the record-only form because ``VariantRecord`` has nowhere to
+    put the status: ``normalisation_ops`` records what happened, and a record that
+    did not move is ambiguous between "already left-most" and "the reference could
+    not be read". The batch report needs the second, so it has to be returned.
+    """
+    coordinate = record.coordinate
+    canonical = canonicalise_allele(
+        contig=coordinate.contig,
+        position=coordinate.position,
+        ref=coordinate.ref,
+        alt=coordinate.alt,
+        reference=reference,
+    )
+    return _apply_canonical(record, canonical), canonical.reference_status
```

with `trim_and_left_align` refactored to share `_apply_canonical` so there is one
implementation (ADR 0018 forbids a second).

---

## 2. Patient data in exception messages

Fixed in `models/variant.py:_phase_consistency` and
`models/pair.py:_pair_variants_distinct`. `hide_input_in_errors=True` suppresses
only pydantic's own `input_value=` suffix; it does not touch a string a validator
built by hand, and a `model_validator(mode="after")` runs with the whole record in
scope. Both now name the fields and the constraint, carry an `error_token` handle,
and echo nothing. A new AST lint over `src/mva/models/` (`test_error_message_privacy.py::test_no_model_validator_interpolates_a_record_field_into_a_message`)
replaces a per-line regex that could not match either defect, because both were
written as multi-line parenthesised assignments.

### 2a. `src/mva/phenotype/loader.py` — an additional finding, not in the review

The same AST sweep over the whole of `src/mva` found two more sites in the
loader for the **patient's own phenotype file**:

* `loader.py:200` — `f"{location}: duplicate HPO term {hpo_id} (first seen on line {seen[hpo_id]})."`
* `loader.py:235` — `f"{location}: HPO term {hpo_id} has an empty label."`

An HPO term drawn from the proband's phenotype TSV is patient data; a handful of
specific terms is identifying in the same way a handful of rare coordinates is.
`PhenotypeObservation`'s own validator already tokenises (pinned by an existing
test), so these two are the remaining raw interpolations on that path.

```diff
-            msg = (
-                f"{location}: duplicate HPO term {hpo_id} (first seen on line {seen[hpo_id]}). "
-                "Two rows for one term cannot be reconciled by a parser; resolve the "
-                "conflict in the source file."
-            )
+            msg = (
+                f"{location}: duplicate HPO term <hpo:{error_token(hpo_id)}>, first seen "
+                f"on line {seen[hpo_id]}. The term is tokenised rather than echoed "
+                "(PRIV-09); the two line numbers locate it in the source. Two rows for "
+                "one term cannot be reconciled by a parser; resolve the conflict in the "
+                "source file."
+            )
```

```diff
-        msg = f"{location}: HPO term {hpo_id} has an empty label."
+        msg = (
+            f"{location}: HPO term <hpo:{error_token(hpo_id)}> has an empty `label` "
+            "column. The term is tokenised rather than echoed (PRIV-09); the line "
+            "number locates the row."
+        )
```

`loader.py` does not currently import `error_token`; add
`from mva.models.base import error_token` alongside the existing
`from mva.models.phenotype import ...`.

**Not** a defect, and deliberately left alone: `phenotype/ontology.py:117` and
`phenotype/hpo.py:535` interpolate HPO IDs too, but those read the public HPO
ontology release and the public gene-phenotype association file. Those IDs are
reference data, not observations about a person, and redacting them would make a
knowledge-file error undebuggable for no privacy gain. GP-14's converse applies:
overstating a limitation is its own dishonesty.

The model-layer lint is scoped to `src/mva/models/` for exactly this reason —
extending it repo-wide needs the public/patient distinction encoded first, since
a blanket rule would either flag the ontology loaders or be silently disabled.

---

## 3. Network isolation: what is now proven, and what cannot be

Fixed in place — see `tests/unit/test_architecture.py`, `docs/privacy-model.md`
("Network denial: honest limits") and TD-06. Summarised because the *scope* is
the deliverable:

**Now enforced.** No module in `ingestion`, `annotation`, `phenotype` or
`prioritization` imports a network client, `socket`/`socketserver`/`ssl`/
`asyncio`, or `ctypes`/`cffi`/`subprocess`/`importlib`, except two entries in
`PATIENT_PATH_IMPORT_EXEMPTIONS` that each carry their reason and residual risk.
Three new tests: the matcher is exercised against sources written to defeat it
(the old lint was only ever run against a tree it already passed on, which is how
`socket` and `ctypes` went unlisted while the docstring claimed a structural
guarantee); every exemption is checked to still correspond to a real import; and
the set of patient-path modules permitted to spawn a process is pinned to exactly
`{annotation/snpeff_local.py}`.

**Still not true, and now said so in three places.** The exempted module spawns a
JVM (`snpeff_local.py:1600`) and writes a VCF of proband coordinates to its
stdin. Audit hooks are per-interpreter, so `netguard` cannot see it;
`OfflineProfile(strict=True)` would block the spawn but the pipeline runs
non-strict (`orchestrator.py:181`) because strict also blocks the `git` calls the
provenance manifest needs. `-nodownload/-noStats/-noLog` are flags to a
cooperative program, not a boundary around an uncooperative one. Separately,
`pysam`/`cyvcf2` are C extensions on the same path and emit no audit event at
all; the compensating control there is the adapters' refusal of non-`file://`
URIs, which is a real check but a different one.

**The feasible OS control on macOS, stated but deliberately not implemented.**
Wrap the *whole* invocation, so the JVM inherits it — a child cannot escape its
parent's Seatbelt profile, which is exactly why it has to be outside the Python
process:

```bash
sandbox-exec -p '(version 1)(allow default)(deny network-outbound)' \
  uv run mva run all --config <case>.yaml --workspace "$MVA_WORKSPACE"
```

`sandbox-exec` is deprecated by Apple but functional on the target machine; a
`pf` deny rule scoped to the run user is the non-deprecated equivalent. It is not
wired into `justfile` or the CLI on purpose: a wrapper the program applies to
itself is a wrapper the program can fail to apply, and `NetworkProfile.OFFLINE_ENFORCED`
would then read as verified when nothing verified it. Until something in the run
can *observe* the control, it stays an operator assertion and the CLI prints it as
one.

---

## 4. GP-30: what is byte-reproducible, and what is not

GP-30 says "repeat runs are byte-identical". For a **real** case that is false as
written, and the docs that repeat it are overclaiming. This section is the honest
scope; TD-21 tracks the fix.

### Why it holds for the demo and not for a real case

`pipeline.py:450` is the whole of the clock decision:

```python
active_clock = clock or (demo_clock() if config.synthetic else SystemClock())
```

`config.synthetic` is the only switch. The `clock` parameter (`pipeline.py:441`)
exists but **no caller ever passes it** and there is no `--clock` flag. The only
case config in the repo is `config/synthetic-case.yaml` (`synthetic: true`), and
`just demo-determinism`, `tests/integration/test_determinism_cross_process.py`
and the scale harness all hard-code it. So every determinism check ever run has
been frozen at `2026-01-01T00:00:00Z`.

`verify_determinism` (`pipeline.py:573`) compares `{relative_path: sha256}` and
skips **whole files** by suffix — `provenance.json`, `*.duckdb`, `*.duckdb.wal`.
It does not freeze the clock and does not mask fields. Pointed at a real case it
would report a failure, on ten artifacts.

### Artifact-by-artifact

**(a) Byte-identical even under `SystemClock`** — 17 artifacts:
`variants/normalised.json`, `variants/annotated.json`, `qc/qc_report.json`,
`submission/track1_submission.csv`, `privacy/privacy_audit.md`, and
`evidence/parquet/{variants, consequences, frequencies, clinical_assertions,
phenotype_observations, genes, mechanism_nodes, mechanism_links, mechanisms,
citations, graph_edges, pipeline_runs, artifact_provenance}.parquet`.

**(b) Differ ONLY in recorded time** — 11 artifacts:

| Artifact | Time-bearing bytes |
|---|---|
| `candidates/ranked_pairs.json` | `rank_rationale`, one per row (`ranking.py:169`) |
| `reports/candidate_dossier.md` | `generated_at` + the same `rank_rationale` |
| `reports/mechanism_report.md` | `generated_at` (`track2.py:438`) |
| `reports/drug_hypotheses.md` | `generated_at` ×2, plus a **date** stamp |
| `reports/rejection_record.md` | `generated_at`, plus the date stamp |
| `reports/track2_report.md` | `generated_at` ×3, plus the date stamp |
| `evidence/parquet/evidence_items.parquet` | `timestamp`, `timestamp_iso` columns |
| `evidence/parquet/candidate_pairs.parquet` | `rank_rationale` column |
| `evidence/parquet/drugs.parquet` | `rejection_rationale` (date stamp) |
| `evidence/parquet/drug_rejections.parquet` | `rationale` (date stamp) |
| `provenance.json` | `started_at`, `completed_at`, `artifacts[].created_at` ×~30 |

**(c) Differ in scientific content — none.** This is the strongest positive claim
available and it should be made, because it is the one that matters for a rerun:
no RNG anywhere on the data path; no timestamp feeds any identifier or hash
(`run_id` from `config_hash` + reference fingerprint, `artifact_id` from the
relative path, `evidence_id` a blake2b of subject/category/claim/tool with the
timestamp deliberately excluded, `edge_id` content-derived); every ordering is
total and time-free; verified across `PYTHONHASHSEED` and `TZ` variation by
`tests/integration/test_determinism_cross_process.py`.

`ArtifactProvenance.content_hash` differs for the (b) artifacts, but *because*
the recorded time differs — no hash is computed over a timestamp as an input, and
no ordering, id or score depends on one.

### The midnight trap

`interventions/safety.py:323` bakes a **date** rather than a full timestamp into
free text:

```python
f"{mechanism.mechanism_id} as of {clock.now().date().isoformat()}: "
```

It reaches five artifacts. A future real-clock determinism check would therefore
be green at 14:00 and red at 00:00 — the worst possible failure signature,
because it looks like flakiness and ADR 0009 correctly says a flaky test here is
a determinism bug. Whoever fixes this should fix `safety.py` first.

### Remediation

Two options; the first is cheaper and more honest.

**Preferred — a run-scoped fixed clock, chosen explicitly.** Thread the existing,
never-used `clock` parameter to the CLI:

```diff
@@ src/mva/cli.py — `run all` and `verify determinism`
+    reproducible: bool = typer.Option(
+        False,
+        "--reproducible",
+        help=(
+            "Freeze the run clock at the case's config hash epoch so repeat runs are "
+            "byte-identical, including reports and provenance. Recorded in the run "
+            "manifest so a reader knows the timestamps are nominal, not observed."
+        ),
+    ),
```

passing `clock=FixedClock(...)` into `execute_pipeline`. The important part is
the manifest field: a frozen timestamp that is *not* labelled as frozen is a
second provenance lie, and this whole exercise is about not adding those.
`verify determinism` should set it unconditionally — it exists to compare bytes,
and comparing bytes under a moving clock is not a check, it is a coin toss.

**Alternative — a timestamp-excluded digest.** Replace the whole-file suffix skip
in `verify_determinism` with field-level normalisation of the known time-bearing
keys (`generated_at`, `timestamp`, `timestamp_iso`, `started_at`, `completed_at`,
`created_at`, and the `Ranked N of M at ...` / `as of YYYY-MM-DD` suffixes).
Weaker: it needs a maintained list of keys, and a new timestamp added anywhere is
silently excluded rather than caught. Prefer it only as a supplement.

Either way, **GP-30 must be quoted with its scope until one exists.**

### Documents that overclaim

Corrected where the file is owned; listed here where it is not, so nothing is
lost:

| File:line | Quote | Status |
|---|---|---|
| `docs/golden-principles.md:97` | "No wall-clock timestamps ... in artifact content" | corrected |
| `CLAUDE.md:30` | "Repeat runs are byte-identical. No wall clock" | corrected |
| `README.md:92` | "Repeat runs are byte-identical, which is checked, not asserted." | corrected |
| `README.md:166` | "Same inputs and config produce byte-identical artifacts." | corrected |
| `docs/architecture.md:106` | "Byte-identical repeat runs are an acceptance criterion, not an aspiration" | corrected |
| `docs/validation-plan.md:16, 22` | "repeat runs are byte-identical"; "Level 2 — Determinism and reproducibility (**DONE**)" | corrected |
| `src/mva/determinism.py:4-6` | "every artifact this pipeline writes is byte-identical. That is checked directly by the repeat-run test" | corrected |
| `src/mva/models/provenance.py:5-7` | "must produce byte-identical outputs, and the golden tests assert precisely that" | corrected |
| `src/mva/pipeline.py:194` | `"""Write canonical JSON so repeat runs are byte-identical (GP-30)."""` | corrected — the docstring now says it stabilises *serialisation*, not *content* |
| `docs/decisions/0004-evidence-ledger.md:29` and `docs/architecture.md:62` | "Content-derived IDs ... is also what makes repeat runs byte-identical" | corrected in both — the ID half kept, the trailing clause replaced by what `EvidenceItem.timestamp` actually does |
| `docs/decisions/0002-duckdb-parquet-evidence-store.md:39` | "embed no timestamps. This is asserted in a test." | corrected — scoped to the writer, with the four tables whose *data* carries a time |
| `docs/maturity-ledger.md:28` | `cli`/`pipeline` row has an empty "not real" cell | already corrected by the ledger's owner; the cell now states the synthetic-only scope |
| `docs/scale-report.md:585` | "Determinism (GP-30), proved not asserted" | corrected — heading and lead now scope the proof to the fixed clock it ran under |
