# ADR 0029 — The pre-commit hook blocks on `fail`, not on every warning

**Status:** accepted · **Date:** 2026-10-25 · **Amends:** PRIV-01 · **Extends:** ADR 0009

## Context

`just install-hooks` wrote a pre-commit hook that ran
`mva privacy audit --staged --strict`. Under `--strict`, `_result` in
`mva/privacy/audit.py` treats any `warn`-severity finding as disqualifying:

```python
passed = not has_fail and not (strict and has_warn)
```

The hook therefore could not pass on a clean tree. Three documentation files and
one source file, all free of patient data, were refused on five findings:

| Staged file | Rule | Why it matched |
|---|---|---|
| `AGENTS.md` | `iso_date_bare` | the gate-timing table's measurement date |
| `docs/privacy-model.md` | `iso_date_bare` | a dated paragraph |
| `workflow/containers/Dockerfile` | `iso_date_bare` | the date the base-image digest was verified against docker.io |
| `src/mva/cli.py` (×2) | `hpo_term` | `HP:0000252` in a commented-out phenotype-template example row |

Neither rule was asking to be fatal. Both say so where they are defined:

* `iso_date_bare.false_positive_risk` — *"Very high by construction: provenance
  timestamps, lockfiles and changelogs are made of these. Kept as `warn` for
  visibility only."* Its pattern is a bare ISO-8601 date with no keying, so it
  matches the `**Date:**` line of all 29 records in this directory, including this
  one.
* `hpo_term` is one of `CONTEXT_GATED_RULES`. The whole-file logic in
  `_refine_severity` already decides the question — `fail` at
  `HPO_DISTINCT_FAIL_THRESHOLD` distinct terms, `warn` below it — and had answered
  `warn` at 2. `--strict` overrode that answer with a blanket promotion, which
  discards the refinement the rule exists to perform.

`run_audit`'s own docstring had already scoped the flag correctly: *"a strict run
is for CI on a release branch, not for the inner loop."* The hook **is** the inner
loop. The defect was the hook's invocation, not the detectors and not the
threshold.

A gate that cannot pass on a clean tree does not make a repository safer. It
teaches everyone who works in it to reach for `--no-verify` by reflex, and a
bypass that is habitual is a bypass that is also used on the one commit that
mattered. That is the real exposure here, and it is larger than the one this
configuration was buying.

## Decision

The generated hook runs `mva privacy audit --staged`. It blocks on `fail`
severity. Warnings are still counted and printed; they are no longer fatal.

`just privacy-audit-staged` keeps `--strict` and is unchanged, so the strict read
is still one command away for CI or a release branch — the use the docstring
reserves it for.

**Nothing was weakened to reach this.** No detector, pattern, severity floor,
threshold or allowlist was touched. `CONTENT_ALLOWLIST_PREFIXES` is untouched and
still scoped to `tests/fixtures/synthetic/`. The audit's own remediation text —
*do not weaken a check to make it pass* — is the reason the fix is at the call
site.

## What still blocks a commit

Every rule that marks an actual disclosure carries `fail` severity and so fails
the hook without `--strict`: `vcf_header`, `vcf_chrom_line`, `genotype_field`,
`person_name_keyed`, `mrn_keyed`, `dob_keyed`, and the `magic:*` container rules
that catch BAM/CRAM/BCF by signature rather than by extension. So does a
`warn`-floor rule the whole-file logic promotes: a `vcf_data_line` in a file that
also matches a VCF header or a genotype field, or `hpo_term` at or above the
distinct-term threshold without a reviewed-path exemption.

Verified rather than asserted. A fabricated four-line VCF — `##fileformat`
header, `#CHROM` row, one `0/1:10,9:19` genotype record, no real data — was staged
and the non-strict staged audit **failed it, exit 2**, on ten findings. The same
audit exits 0 on the four files above. The probe was removed and is not in the
tree; `git_staged_sensitive` reports both runs under
`tests/unit/test_privacy.py`'s existing coverage of the staged path.

## Consequences

* The hook passes on a clean tree, so `--no-verify` returns to meaning what it
  should: a deliberate, explained exception.
* A privacy finding that is genuinely advisory now reads as advisory. Five
  findings printed under a `pass` is information; the same five under a `FAIL`
  that everyone bypasses is noise.
* `iso_date_bare` and an under-threshold `hpo_term` no longer stop a commit. This
  is the accepted cost, and it is the cost of a rule the detector set never
  intended to be blocking. The strict recipe remains for when the stronger read is
  wanted.
* The asymmetry that hid this is gone: `just verify`'s `privacy-audit` step was
  never strict and ran with an empty index, so the gate stayed green while the
  hook could not pass. Both paths now agree about what blocks.
