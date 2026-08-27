# ADR 0016 — "Impact not assessed" must be representable

**Status:** accepted
**Date:** 2026-08-28
**Touches:** GP-14, GP-20, ADR 0005

## Context

`VariantRecord.gene_symbols` is derived purely from `consequences`, and
`pairing._variants_by_gene` groups on it. With no consequence annotator the
pipeline emits **zero candidate pairs** — verified by execution:

```
gene_symbols per record: [(), ()]
generate_pairs -> ()          # two het variants inside BUB1B
```

For a challenge whose answer is a compound-heterozygous pair, that is a total
failure. Gene assignment is therefore load-bearing.

A MANE interval join answers *"which gene is this variant in?"* — a locational
fact — but says nothing about molecular consequence. `ConsequenceAnnotation.impact`
was a required `ImpactSeverity`, so such an adapter could not construct one
without **inventing a severity it had not computed**. An adversarial review caught
the first attempt doing exactly that behind a constructor argument, while
declaring `synthetic = False`.

Tracing every enum value showed no honest choice existed:

| value | `scoring._IMPACT_BASE` | in `filters.BENIGN_IMPACTS`? | effect |
|---|---|---|---|
| MODIFIER | 0.05 | **yes** | fabricates a `benign_consequence` flag |
| LOW | 0.15 | **yes** | same |
| MODERATE | 0.50 | no | invents severity |
| HIGH | 0.90 | no | invents severity **and** promotes every lone heterozygote |

Two of the four actively assert benignity the source never claimed; the other two
assert severity it never claimed.

## Decision

**`ConsequenceAnnotation.impact` becomes `ImpactSeverity | None`, defaulting to
`None`, which means NOT ASSESSED.**

`None` is emphatically **not** `MODIFIER`. MODIFIER is a positive prediction of
negligible effect and the filter layer treats it as one. Collapsing "nobody
predicted" into "predicted harmless" is precisely the GP-14 failure this field
now exists to keep representable — and it is the failure mode that would quietly
discard a real candidate.

`worst_impact_for_gene` drops unassessed transcripts rather than counting them as
least-severe, so a gene-assignment-only annotation cannot dilute a real
prediction from another adapter.

## Why the consumers needed no new logic

The pipeline already handled the honest answer correctly:

- `filters._consequence_flag` returns early on empty consequences, and
  `None in BENIGN_IMPACTS` is `False` — so an unassessed variant is never flagged
  benign. Verified.
- `scoring._variant_consequence` sees `worst_impact_for_gene() -> None` and
  returns `_IMPACT_UNANNOTATED`, with the rationale *"no consequence annotation
  for this gene (unknown, not benign)"* — text that was already written for this
  case before anything could produce it.
- `pairing._IMPACT_ORDER_UNKNOWN` already sorts unannotated alongside LOW, on
  exactly these grounds.

Only three sites needed changing, all rendering or persistence:
`service._impact_text` (renders "not assessed"), the evidence strength for an
unassessed annotation (`INSUFFICIENT`), and the Parquet writer (nullable column).

## Consequences

- A gene-assignment adapter can now be honest and still be useful: it supplies
  `gene_symbol`, which is what `generate_pairs` needs, and withholds what it did
  not compute.
- Evidence text renders "predicted impact: not assessed", so a reader is never
  shown a severity nobody computed.
- `tests/unit/test_gene_intervals.py` goes from 102 passed / 19 skipped to
  **121 passed / 0 skipped** — the skipped tests were gated on this model change.
- **A future adapter must not reintroduce a default.** The value of this field is
  that it distinguishes three states — assessed-and-severe, assessed-and-benign,
  and not-assessed — and any default collapses the third into one of the others.
