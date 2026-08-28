# ADR 0018 — One allele-canonicalisation rule, at layer 1

**Status:** accepted
**Date:** 2026-08-28
**Touches:** GP-01, GP-03, GP-14, ADR 0010

## Context

An independent review found that variant *representation* did not agree between
ingestion and the annotation adapters, so equivalent variants failed to join.
Reproduced by execution:

```
clinvar record key : GRCh38:chr1:100:AT:AG
proband trimmed key: GRCh38:chr1:101:T:G
join result        : {}          <- the pathogenic assertion silently gone
```

and, separately, with no reference FASTA in the production path a proband
insertion stays right-shifted at `chr1:101 A>AA` while gnomAD stores the
left-aligned `chr1:100 A>AA`.

**A failed join is indistinguishable from "novel, ultra-rare variant"** — the
strongest promoting signal in our ranking. It simultaneously manufactures false
positives and discards true pathogenic assertions, and it does both silently.
**A large minority of the records in a real short-read callset are
indel-bearing** — this one included — so a rule that silently mishandles indels
is a main-path defect, not a corner case.

The root cause was two implementations of one rule.

## Decision

**Minimal representation and left-alignment live in exactly one module,
`src/mva/alleles.py`, at layer 1**, implementing Tan, Abecasis & Kang 2015
(*Bioinformatics* 31(13):2202-4, Algorithm 1) and VCF 4.2 §1.4.1/§5. Ingestion
and every annotation adapter call it; none keeps a copy. That was stated here
before it was true: the gnomAD adapter retained a private `minimal_representation`
for one release cycle after this ADR was accepted. It trimmed correctly and could
not left-align, so a proband indel that ingestion had left-aligned was looked up
under a key gnomAD does not hold — in the highest-weight signal the ranker has.
The adapter now takes the same optional `ReferenceLookup` as ingestion and
ClinVar, and a structural lint over `src/mva/annotation` fails on allele string
surgery outside `mva.alleles`, so the claim is enforced rather than asserted.

Layer 1 is forced, not chosen. `ingestion` and `annotation` are peer stages at
layer 4 and GP-03 forbids them importing each other, so a rule both must share
has to sit below both. It is a pure function of (coordinate, alleles, optional
reference lookup) and depends only on `models`, which is what layer 1 is for.

**The layer map now enforces this rather than merely permitting it.** `alleles`
is a mapped owner, and a new test makes an *unmapped* module a failure instead of
a skip. That test also removed the phantom `knowledge` entry and added
`orchestrator` — which, because both layer tests skipped unmapped owners, meant
the one module importing every stage was exempt from GP-01/GP-03 entirely.

## Absence of a reference is a declared state, not a silent skip

Left-alignment is impossible without the reference genome. That now produces a
typed `LeftAlignmentReport` with an explicit status
(`APPLIED` / `NOT_REQUIRED` / `UNAVAILABLE_NO_REFERENCE` /
`INCOMPLETE_REFERENCE_UNUSABLE`), reaching a reader through the run warnings,
`RunManifest.warnings`, and `qc/qc_report.json`.

`NOT_REQUIRED` — a callset with no indels — is deliberately **not** degraded.
GP-14 cuts both ways: overstating a limitation is its own dishonesty.

## An honest correction to the review's framing

The pinned ClinVar release was scanned in full: **0 non-minimal ALT entries in
4,467,990 records.** So the ClinVar-side asymmetry is a genuine code-level defect
— reproduced above with a constructed record — that costs nothing on *this*
release. The half that bites today is the production path, and it bit twice over: a
right-shifted proband indel missed ClinVar *and* gnomAD, and until the gnomAD
adapter was given a reference it missed gnomAD **even when ingestion had
left-aligned correctly**, because the adapter re-derived its own key.

Measured on the real v4.1 exomes chr21 shard against GRCh38: of 1,382 indel
records in a 520 kb exonic window, **1,029 (74.5%) sit in a repeat tract** and so
have at least one other legal spelling. Under the private rule **0 of those
1,029** right-shifted spellings joined; under the shared rule with a reference
**989** do — the other 40 being the AN=0 records GP-14 requires to stay omitted.
**30 of the recovered records are ones gnomAD calls common** (max population
AF > 1%), 22 of them above 5%: variants that were being scored as novel and
ultra-rare.

Both facts are pinned as tests, so a future release that stops being minimally
represented fails a test instead of quietly losing assertions.

## Known limits, stated rather than hidden

- **Symbolic and missing alleles** (`*`, `.`, `<DEL>`) are returned verbatim.
  There is no defined shift for them, and mangling one moves a coordinate on a
  guess.
- **Complex substitutions** (both alleles >1 base after trimming) are not
  shiftable indels; leaving them is the standard, not a gap.
- **`MAX_SHIFT_BP = 1000`.** A tandem repeat longer than 1 kb stops
  un-proven-leftmost.
- The right-shift bound may over-reach by one base. Over-reaching costs one index
  query; under-reaching would re-open the silent miss. The asymmetry is deliberate.
- **A region query is part of the rule.** Fixing a join key does not help if the
  fetch window never reaches the record. Both adapters query from the left-most
  position out to `rightmost_equivalent_position`; a window stopping at the query
  position loses a right-shifted record to the *fetch* rather than to the key, and
  from outside the two failures are indistinguishable.

## Consequences

- A `hypothesis` property test (400 examples × 2 reference states) asserts that
  ingestion, the ClinVar adapter **and the gnomAD adapter** can never disagree.
  The previous version compared only two of the three, which is precisely how the
  gnomAD copy survived this ADR. A structural lint in the same module fails any
  adapter that reintroduces a local trim, matching on the *shape* of allele string
  surgery rather than on variable names — because Tan Algorithm 1 *is*
  `ref[-1] == alt[-1]` and a name list only catches the copy spelled the way the
  last one was.
- A moved POS is recorded in the record's normalisation operations. Our submission
  is scored on exact coordinates, so an untraceable coordinate change is a
  correctness risk as well as an auditability one.
