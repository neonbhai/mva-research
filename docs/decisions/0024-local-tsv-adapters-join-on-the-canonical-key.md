# ADR 0024 — The local TSV adapters join on the canonical key

**Status:** accepted
**Date:** 2026-08-28
**Touches:** ADR 0018 (one canonicalisation rule), GP-11, GP-13, GP-14, GP-20

## Context

ADR 0018 established that minimal representation and left-alignment live in
exactly one module, `src/mva/alleles.py`, and that **every join key goes through
it**. Ingestion, the ClinVar adapter and the gnomAD adapter were converted. The
local TSV adapters in `mva.annotation.local_tables` were not.

They stored the table's `variant_id` column verbatim and looked it up by direct
string membership:

```python
grouped.setdefault(variant_id, []).append(annotation)   # index
...
if variant_id in self._index                            # lookup
```

So a table holding `GRCh38:chr1:100:AT:AG` did not answer a query for the
equivalent minimal `GRCh38:chr1:101:T:G`, and neither side ever noticed. **A failed
join does not raise. It returns absence** — which this pipeline reads as "no
frequency data" (GP-14) and, for a missing consequence, as grounds for dropping the
allele in selection. A failed join is indistinguishable from "novel and
ultra-rare", the exact profile of a causal variant.

The tables shipped here are synthetic (GP-20), which is why this was missed and
why it still matters: `load_default_adapters` is the **default executable adapter
path**. It is what `just demo` runs and what an organizer reruns, and the challenge
requires reproducibility.

## Decision

**`_join_key` canonicalises both sides of the join, and delegates.** It is applied
to the table's `variant_id` when the index is built and to the caller's ID when it
is looked up, so the two agree by construction rather than by whoever spelled it
first. Three normalisations, each a comparison that would otherwise fail on a
spelling difference:

* build via `GenomeBuild.parse`, so `hg38` and `GRCh38` name one assembly while
  `GRCh37` still never joins `GRCh38` (GP-11);
* contig via `normalise_contig`, so `15` and `chr15` are one chromosome;
* alleles via `mva.alleles.canonicalise_allele`.

There is deliberately **no allele arithmetic in this module**. The structural lint
`test_no_annotation_adapter_defines_its_own_canonicalisation` fails the build over
any that appears, which is what makes ADR 0018 enforced rather than asserted.

An ID that is not `build:contig:pos:ref:alt`, names an unrecognised assembly, or
carries a non-canonical contig is returned **verbatim**. Refusing to guess is the
point: a key reinterpreted on a guess joins the wrong record, which is worse than
joining none.

**Results stay keyed by the caller's own ID string.** `mva.annotation.service`
resolves the mapping with `record.variant_id`; re-keying to the canonical form
would move the silent miss one layer up instead of removing it.

## Trim-only is honest, and says so

A TSV has no reference genome. Trimming needs none and is always correct;
left-alignment means rolling an indel leftwards through a repeat tract, which
requires reading the bases to its left. So `_join_key` passes no `ReferenceLookup`,
and `canonicalise_allele` therefore does not record `left_align` in `operations` —
by design, because "recording an operation that was never performed is a provenance
lie no downstream consumer can detect".

That leaves a real residual limitation: `chr1:100 A>AT` and `chr1:104 T>TT` are the
same insertion into a homopolymer and will still not join. **Declaring it is part
of the fix.** Both adapters publish a `left_alignment` property returning the same
`LeftAlignmentReport` the ingestion stage and the real adapters use, derived from
the index through `summarise_left_alignment` so no adapter can label its own batch
by hand:

* a table containing an indel reports `UNAVAILABLE_NO_REFERENCE` and
  `is_degraded=True`;
* a table of SNVs reports `NOT_REQUIRED`.

"Could not left-align" and "had nothing to left-align" are opposite claims about
how far to trust the rarity of every indel in a run (GP-14), and they must not
share a value. `test_an_indel_written_at_another_position_in_a_repeat_tract_still_misses`
locks the residual miss in as a known limitation rather than leaving a future
reader to infer from the passing join tests that these adapters are
representation-complete.

## Consequences

- `tests/unit/test_local_table_join.py` (14 tests) asserts the agreement
  **directly against `canonicalise_allele`**, not by observing that a join
  happened to succeed. Agreement inferred from a passing join is exactly the
  evidence that was available while the two implementations disagreed.
- Two table rows spelling one variant two ways now **merge** into that variant's
  annotations instead of becoming two entries only one of which a query can find.
- No shipped table changes. Every `variant_id` in `knowledge/public/` is already
  minimal and UCSC-spelled, so the committed manifest hashes and the golden case
  are untouched — the fix removes a latent failure, it does not re-baseline
  anything (GP-32).
- Making these adapters take an optional `ReferenceLookup`, as the gnomAD and
  ClinVar adapters do, would close the repeat-tract gap. It is not done here
  because these tables are a synthetic substitute and the report makes the
  limitation visible; see `docs/tech-debt.md`.
