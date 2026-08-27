# ADR 0012 — Public reference fixtures are their own category, not "synthetic"

**Status:** accepted
**Date:** 2026-08-28
**Supersedes:** nothing. **Amends:** the privacy audit's exemption model.

## Context

Real annotation adapters have landed (TD-01). Each needs a test fixture cut from
the release it reads: a slice of the NCBI ClinVar VCF, of gnomAD sites, of the
MANE GTF. These are **real data** — that is the entire point, since a fixture of
invented records would not test the parser against what the source actually
emits.

The repository's `.gitignore` is deny-by-default for genomic formats, with narrow
negations. Section 12 says, verbatim:

> NEVER add a negation outside tests/fixtures/synthetic/ or knowledge/public/.

and the privacy audit backs it with `synthetic_fixtures_marked`, requiring every
negated file to carry `##mva_synthetic=true` or a `SYNTH-` sample ID. A real
ClinVar slice **cannot honestly carry that marker.** It is not synthetic.

So there were three options, and only one of them is honest:

1. **Stretch the synthetic exemption** to cover these files — i.e. mark a real
   ClinVar slice `mva_synthetic=true`. Rejected outright. The marker is a factual
   assertion, and the audit check exists precisely so that a real VCF dropped into
   the fixtures directory **fails** instead of sliding in. Teaching the codebase
   to lie in the one place designed to catch a lie destroys the check's value
   everywhere, and it would do so at the exact moment we are about to handle a
   real proband's genome.
2. **Do not commit the fixtures** — vendor nothing, regenerate on demand. Rejected:
   the tests then require a 185 MB download to run, `just verify` stops being
   hermetic, and CI cannot reproduce a failure.
3. **Add a distinct category with its own, narrower predicate.** Accepted.

## Decision

Introduce **`public-reference-fixture`** as a category separate from
`declared-synthetic`. A file qualifies only if **all four** hold:

1. It lives under an allowlisted prefix — `tests/fixtures/<source>/` where
   `<source>` names a resource registered in `knowledge/manifests/resources.yaml`.
2. It is **derived from a hash-pinned public release**: the manifest entry it
   claims descent from exists and is `status: fetched` with a real `sha256`.
3. It contains **no sample columns and no genotypes**. For VCF this is mechanical:
   the `#CHROM` line has exactly 8 fields (through `INFO`), no `FORMAT`, no sample
   names. A sites-only file cannot carry a genotype by construction.
4. A committed, re-runnable generator script sits beside it recording exactly
   which regions were cut from which release.

Condition 3 is the one doing the real work. The privacy risk being managed is not
"real data exists in the repo" — ClinVar is public domain and already on every
clinical genomics machine on earth. It is **"a file with per-sample genotypes
gets committed"**, because that is what discloses an individual. A sites-only
slice of a public reference release has no individual to disclose.

Condition 2 is what stops the category from becoming a general-purpose escape
hatch: you cannot claim it for a file that is not traceable to a public release
we recorded fetching.

## Consequences

- `src/mva/privacy/audit.py` gains this category with its own predicate. The
  synthetic marker requirement is **unchanged** for `tests/fixtures/synthetic/`.
  The two exemptions do not share code paths, so widening one cannot widen the
  other by accident.
- The `.gitignore` negations for reference fixtures are per-source and per-extension.
  A negation for `tests/fixtures/clinvar/**/*.vcf.gz` does not admit anything else.
- **The audit must fail loudly if a file under a reference-fixture prefix carries
  sample columns.** That is the whole guarantee; it is a test, not a comment.
- The proband's VCF can never qualify: it has sample columns, and it is not in
  `resources.yaml` — patient data is never a registered public resource.

## Rejected alternative worth naming

Allowing the category based on *provenance alone* ("it came from ClinVar, so it
is fine") without condition 3. Rejected because provenance is a claim in a
comment and genotype absence is a property of the bytes. When the two disagree,
only one of them is checkable — so the check is written against the bytes.
