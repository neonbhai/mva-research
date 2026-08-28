# ADR 0026 — A broken reference is a declared state, not a silent skip

**Status:** accepted · **Date:** 2026-08-28
**Touches:** GP-14, GP-31, ADR 0018

## Context

ADR 0018 established that the *absence* of a reference FASTA is a declared state:
without one, left-alignment cannot happen, `left_align` is not recorded in the
record's operations, and a typed `LeftAlignmentReport` carries the limitation to a
reader. That was the right rule and it was incompletely applied. It covered
`reference is None`. It did not cover the reference that is *present and broken*.

`mva.alleles._fetch_base` caught every exception from `ReferenceLookup.fetch` and
returned `None`, and the shift loop read `None` as "no base here, stop". Four
different conditions arrived at that one value:

| Condition | Truth |
|---|---|
| `position < 1` | before the start of a contig — genuine absence |
| an empty read past the end of a contig | the chromosome ends here — genuine absence |
| the accessor raised, or the contig is missing from the index | **the reference is broken** |
| an empty or non-nucleotide read inside a contig | **the reference is broken** |

The last two were converted into the first two, which is to say an I/O error was
converted into a quality downgrade that no caller could observe. Both annotation
adapters then reported left-alignment as `APPLIED` on the strength of their
reference object being non-`None` — `gnomad_sites.py:1247`, `clinvar_vcf.py:628`
— so a FASTA that raised on every read produced trim-only join keys underneath an
explicit claim of reference-backed canonicalisation.

The cost is the one ADR 0018 measured. On the real gnomAD v4.1 exomes chr21
shard, **1,029 of 1,382 indel records (74.5%) in a 520 kb exonic window sit in a
repeat tract**, and without a working reference **0 of those 1,029** right-shifted
spellings join, against 989 with one. Thirty of the recovered records are alleles
gnomAD calls common. A failed join presents as "no frequency data", which the
ranker scores as novel and ultra-rare.

**The provenance lie is worse than the failure it hides.** A run that cannot
left-align and says so is a run whose indel results a reader knows to distrust. A
run that cannot left-align and reports `APPLIED` has removed the reader's ability
to tell.

## Decision

**Two things are separated, and both are carried in the return type.**

1. **What happened** — `CanonicalAllele.operations`, unchanged. It names only
   operations actually performed.
2. **How far the reference could be trusted** — `CanonicalAllele.reference_status`,
   new, one of `NOT_SUPPLIED` / `USABLE` / `UNUSABLE`.

The second exists because the first is ambiguous on its own. A position that did
not move may be already left-most, or may have hit an unreadable reference, and
those have opposite implications for whether an empty join result means the source
holds no record. `CanonicalAllele.left_alignment_proven` is the single property an
adapter branches on before claiming `LeftAlignmentStatus.APPLIED`.

`USABLE` is vacuously true when the rule needed no read at all — a SNV, or an
indel whose last bases already differ. The claim being made is "nothing the rule
asked for was denied", not "a read happened", and for representation purposes
those have the same consequence: the result *is* the left-most spelling. Reporting
such a record as degraded would overstate the limitation, which ADR 0018 is
explicit is its own kind of dishonesty.

**A function whose return type cannot carry the degraded state raises it.**
`_left_shift` and `_required_base` raise `mva.errors.ReferenceUnusableError`.
`canonicalise_allele` catches it at exactly one documented boundary and converts
it into `UNUSABLE`, because its return type *can* carry the state and because a
per-record exception in a whole-callset loop would turn one bad base into a dead
run — GP-13 forbids that for a record-level problem.

**The degraded result is the trimmed form, never a partial shift.** An allele
abandoned half way through a repeat tract sits at neither the input position nor
the left-most one, so it would join against neither the source nor another run of
this pipeline over the same input. Trim-only is at least defined and reproducible.

**The direction of the read decides what an empty answer means.** The left shift
reads at `position - 1` with `position > 1`, inside a contig the record claims to
sit on, so nothing there is breakage. The right shift reads *beyond* the record,
where nothing there is the end of the chromosome. This is not a heuristic; it
follows from which side of the record is being read, and it is the reason
`rightmost_equivalent_bound` reports rather than raises on an empty read.

## The query window, and why it reports instead of raising

`rightmost_equivalent_position` decides how far right an index query reaches. A
bound that stops early makes the fetch window too small, and a source record lost
to the *fetch* is indistinguishable from one lost to the key, which is in turn
indistinguishable from the source genuinely having nothing.

It returns a bare `int`, so it cannot say "I could not compute this". The obvious
fix — raise — is wrong here, and a test says why:
`test_a_reference_that_raises_cannot_leak_a_coordinate_or_break_the_lookup` pins
that an unreadable reference must not kill the adapter's lookup, and that is the
right call. A lookup that dies on one bad base is worse than one that degrades and
declares it.

So `rightmost_equivalent_bound` returns a typed `QueryBound` carrying both the
position and its `reference_status`, and `rightmost_equivalent_position` remains
as an `int`-returning shim over it for the two adapters that call it, with a
docstring stating exactly what it cannot tell a caller. Migrating those two
adapters is one line each.

An under-reaching bound is not, on its own, a new failure: an unreadable
reference has already left the *query key* trimmed-only, so a right-shifted
source record would not have joined even if it had been fetched. What matters is
that the caller can tell.

## Honest limits

- **The adapters still report `APPLIED` today.** The primitives are honest; the
  two callers that consume them are owned elsewhere. The exact diffs are in
  `docs/handoff-integrity.md` §1a and §1b, with the tests to add. Until they are
  applied, the run-level claim remains wrong even though the per-allele value is
  now right.
- **A rightward empty read cannot be classified.** A FASTA returning nothing
  everywhere is indistinguishable, on a rightward read alone, from a variant on
  the last base of a contig. The *leftward* read catches the same FASTA, so the
  record is still reported `UNUSABLE` — but the bound is reported proven. Pinned
  as a test rather than left implicit.
- **`MAX_SHIFT_BP` exhaustion is still silent.** A tandem repeat longer than 1 kb
  stops un-proven-leftmost with a perfectly healthy reference. ADR 0018 states
  this as a known limit; it is *not* folded into `UNUSABLE`, because attributing
  it to the reference would misname the cause.
- **Ingestion has a narrow residual gap.** `_reference_matches` reads the record's
  REF span first and disables alignment on any failure, so a wholly-broken FASTA
  is already caught. If the REF span reads and a shift base does not,
  `reference_consulted` stays `True` and the batch reports `APPLIED`. Diff in
  `docs/handoff-integrity.md` §1c.

## Consequences

- `mva.errors.ReferenceUnusableError` is new. Its message is PRIV-09-safe: an
  `error_token` locus handle, the failing accessor's exception *type*, and no
  contig or position. The original exception is suppressed with `from None`
  rather than chained, because a genomics backend routinely puts the region it was
  asked for into its own message and a chained traceback carries that everywhere
  the traceback goes.
- Twenty cases in `tests/unit/test_allele_reference_integrity.py` pin the
  distinction, including that the degraded result is byte-identical to the
  no-reference result and that the four failure modes read differently from each
  other — a message that only says "reference error" moves the debugging cost
  rather than removing it.
- No behaviour changes for a healthy reference. Every existing allele test passes
  unmodified, which is the point: this makes a silent failure loud without making
  a working path noisier.
