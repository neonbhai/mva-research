# ADR 0027 — A non-synthetic case binds the real adapters or refuses to run

**Status:** accepted · **Date:** 2026-08-28

## Context

Until this change the composition root called `load_default_adapters(...)` —
the **synthetic** TSV adapters under `knowledge/public/` — and
`normalise_variants` without a reference FASTA, for every case, unconditionally.
The real adapters existed, were tested, and were bound by nothing.

Three consequences, all invisible from the outside:

1. **The executable pipeline did not implement the documented method.** An
   organiser following `docs/track1-method.md` and the CLI in `README.md` would
   not reproduce our ranking, because the code they ran annotated against
   fictional genes (`SYNTH*`) and invented allele frequencies. This was flagged in
   review as the single most serious defect in the repository.
2. **No indel could join.** With no reference, `normalise_variants` cannot
   left-align, so a right-shifted proband indel keeps a spelling ClinVar and
   gnomAD do not hold. A failed join is indistinguishable from "no record", which
   this pipeline scores as novel and ultra-rare — its strongest promoting signal
   (ADR 0018, GP-14).
3. **Both joining adapters take `reference=` as an OPTIONAL keyword.** Wiring them
   without it produces a fully "wired" pipeline that is still in the degraded
   state ADR 0018 was written to close, and every test still passes.

What the third one costs, measured on the real releases rather than estimated:

| Source | Window | Indels in repeat tracts | Join **without** a reference | **With** one |
|---|---|---|---|---|
| gnomAD v4.1 exomes chr21 | 520 kb exonic | 1,029 of 1,382 | **0** | 989 — 30 of them variants gnomAD calls *common*, 22 above 5% AF |
| ClinVar 2026-08-22 | chr17:43,000,000-43,520,000 | 2,215 of 3,595 ALT alleles | **0** | 2,211 — **1,761 Pathogenic / Likely pathogenic** |

## Decision

**1. The case's `synthetic` flag selects the adapter set, and nothing else does.**

* `synthetic: true` → the hash-pinned demo tables, exactly as before. The demo,
  `just demo-determinism` and the whole suite keep running with no resource root,
  so a green suite means the same thing on every machine. Fictional `SYNTH*` genes
  have no place in a real gene model anyway.
* `synthetic: false` → ClinVar, gnomAD, SnpEff and MANE, resolved from
  `config.resources` under `$MVA_RESOURCES`, **required**.

**2. `reference=` is passed to both joining adapters, and to `normalise_variants`.**
`build_real_adapter_set` takes `reference` as a keyword with **no default**, so
the trap cannot be re-entered by omission. One `FastaReference` handle serves
normalisation and both adapters — it is one immutable file, and two handles would
be two caches over it.

**3. A real case with missing resources FAILS. It never falls back.**
No `$MVA_RESOURCES`, an unregistered file, an unfetched manifest entry or a failed
integrity pin each raise before the VCF is opened. Every message names what was
refused *and what it refused to substitute*.

**4. A degraded representation is still representable, and says so.**
`reference=None` remains a legal argument, because refusing to model a degraded
state only relocates the degradation. When it is `None`, each adapter's own
`representation_limitation` is surfaced into `context.warnings` with the measured
cost above. Under decision 3 a real run always has a reference, so a
representation warning on a real run is a tripwire, not an expected state.

**5. A coverage hole in a release is reported, not fatal.**
gnomAD v4.1 exomes ships **no chrM shard** — a property of the release, not a
failed download — so `frequencies()` raises on any mitochondrial variant even
against a complete install. `PartialCoverageFrequencyAdapter` wraps
`lookup_partial()`: answerable variants are answered, unanswerable ones are
absent **and counted**, and the run warning names the contigs and counts (never a
coordinate — PRIV-09). Absence of a resource is not evidence of rarity (GP-14).

**6. SnpEff is pinned from `snpeff_pins.json`, never from `measure()`.**
`measure()` hashes whatever is installed, pinning the install to itself and
verifying nothing. `SnpEffArtifactPins.from_manifest` reads digests the installer
wrote only after checking each artifact against the reviewed constants in the
script.

**7. The synthetic demo runs with `selection.enabled: false`.**
The candidate-selection stage (ADR 0019) exists to make a 4.5 M-record callset
pairable. The demo fixture is twelve records: there is no scale problem to solve,
and enabling it deletes the two common variants at chr15:40205000 (AF 0.120) and
chr15:40206000 (AF 0.082) — which are not noise in that fixture, they **are** the
GP-13 demonstration that `tests/golden/test_golden_case.py` locks. It would also
drop the submission from ten rows to nine, and rows below the answer cannot lower
either challenge metric. Disabled does not mean absent: the stage still runs,
counts, writes `selection/selection_report.json` and warns that it is off.

## Alternatives rejected

**Fall back to the synthetic tables with a warning.** The failure mode is a real
proband ranked against fabricated evidence, producing a submission, a dossier and
a provenance manifest that all look entirely healthy. A warning in a scrollback is
not a match for that. A loud failure costs minutes; a quiet one is not recoverable
after the answer has been submitted (GP-20).

**Bind the real adapters for synthetic cases too when resources are present.**
Rejected: it makes `just verify` depend on whether a machine holds 200 GB of
reference data, so the suite would pass in one place and fail in another — and the
`SYNTH*` genes exist in no real gene model, so the demo would rank nothing.

**Make the reference FASTA optional on a real run.** Rejected on the numbers in
the table above. 1,761 Pathogenic/Likely-pathogenic assertions unreachable is not
a degraded mode worth offering as a convenience.

**Require a fixed list of 24 gnomAD contigs.** Rejected: which contigs a release
ships is a property of the release. Requiring chrM would fail every complete
install. The required set is derived from what the manifest registers.

## Consequences

* `mva run all --config <a real case>` now needs `$MVA_RESOURCES` and ~50 s of
  binding (JVM probe, MANE index, 24 shard opens, spot verification of 33
  registered resources) before the first variant is annotated. Paid once per run.
* The `annotation` package is no longer honestly graded `synthetic-substitute`
  outright; `docs/maturity-ledger.md` now grades it per branch.
* **Not checked:** that the registered releases are the ones the *submission* was
  built from. The manifest pins bytes; matching a submission to a manifest is
  `docs/submission-runbook.md`'s job.
* A real case has still never been executed end to end in this repository: the
  proband data is outside the boundary agents may read (GP-40, ADR 0008), so this
  wiring is verified against the public releases and the synthetic fixtures only.
  That divergence is recorded in `docs/track1-method.md`.
