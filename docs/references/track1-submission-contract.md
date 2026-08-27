# Track 1 submission contract (verified)

Vendored into the repo deliberately. This was derived by reading the challenge
Space's source once; the Space is gated and mutable, so a future run — human or
agent — cannot reliably re-derive it. Treat this file as the contract of record
and re-verify before final submission.

**Source:** raw files under
`https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/raw/main/`
(`evaluation.py`, `config.py`, `groundtruth.py`, `tabs/rules.py`,
`tabs/submit_track1.py`, `static/templates/track1_submission_template.csv`).
Verified 2026-08-27. The rendered Gradio page returns only a shell to a fetcher;
the source files are authoritative.

## CSV format

Header required. Exact column order:

```
proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr,finding_type,notes
```

Template rows, verbatim:

```
PROBAND01,chr1,100000,A,T,chr1,200000,G,C,0.90,primary,
PROBAND01,chr7,300000,C,T,,,,,0.30,secondary,Incidental finding unrelated to primary phenotype
```

| Field | Required | Rule |
|---|---|---|
| `proband_id` | yes | **Only `PROBAND01` is accepted.** Anything else hard-fails with "Unknown proband_id". |
| `chrom_1`,`pos_1`,`ref_1`,`alt_1` | yes | `KeyError` if absent. `pos` cast with `int()`. |
| `chrom_2`…`alt_2` | no | Read with `.get()`. Blank = single-variant proposal. |
| `epcr` | yes | `float`, validated `0 < epcr <= 1`, else `ValueError`. |
| `finding_type` | no | Must be `primary` or `secondary` if present. |
| `notes` | no | Never read by the scorer. |

**Row limit: 10.** `if len(rows) > 10: raise ValueError`. No minimum.

**A compound-het pair is ONE ROW**, using the `_2` columns — not two rows.

**No rank column.** Rank is derived from `epcr` descending, ties broken by file
order:
```python
rows_sorted = sorted(enumerate(rows), key=lambda x: (-x[1][1], x[0]))
```
Pre-sorting is not required, but we do it anyway for determinism.

## Scoring

Variant key: `(chrom, pos, ref, alt)`, built as
```python
(chrom.strip(), int(pos), ref.strip().upper(), alt.strip().upper())
```

> **The single most dangerous detail.** `chrom` is only `.strip()`-ed. There is
> **no chr-prefix or case normalisation**. Emitting `15` where the answer key
> holds `chr15` scores zero while looking correct to a human reader. Ensembl-style
> GRCh38 VCFs use bare contig names, so this is a live hazard, not a theoretical
> one. Our renderer forces UCSC style and a test asserts it.

Matching is exact frozenset equality over the row's variants; the first (best)
full match wins. Partial credit applies to compound-het rows on set intersection.

Two metrics:
- **Rank points** — tiers `[(1,100),(3,50),(5,25),(10,10)]`: rank 1 → 100,
  ranks 2–3 → 50, 4–5 → 25, 6–10 → 10, else 0. A partial match scores
  `0.5 × rank_points`.
- **F-max** — swept over every unique EPCR threshold, computed **per individual
  variant, not per row** (`predicted_variants |= row.variants`).

Ground truth is a private gated dataset (`gold_standard_track1.json`); the answer
is a clinically validated compound-heterozygous pair. Methodology adapted from
the CAGI6 Rare Genomes Project assessment (Stenton et al., 2024).

## Genome build

**GRCh38**, stated explicitly ("Coordinates must use GRCh38"). Chromosome naming
is **`chr`-prefixed**: the docs example is `chr15`, the template uses
`chr1`/`chr7`, and the repo's local fallback ground truth uses `("chr2", …)`.

## Other verified constraints

- Track 1 also requires a report (PDF or Markdown only) and a public GitHub URL.
- Track 2: report (PDF/MD), public GitHub URL, 3-minute video URL. Rubric —
  Scientific Rigor 35%, Potential Impact 25%, Innovation 25%, Scalability 15%.
- Submission limits: Track 1 six attempts (best counts), Track 2 one.
- Close: 2026-10-24 23:59 UTC.
- Data is gated under WCG IRB protocol #REDACTED-PROTOCOL; no resharing; **deletion within
  30 days of close from all environments including derived datasets**, confirmed
  by email. See `docs/privacy-model.md` for how we meet that obligation.
- No stated rule restricting LLM use or external data. Reproducibility IS
  required: organizers may rerun submissions.

## Two gotchas

1. The scorer evaluates only `next(iter(submissions))` — the first proband
   encountered. A stray `proband_id` hard-fails the whole submission.
2. Chromosome strings are compared raw. See the warning above.

---

# Re-verified against the real scorer source — 2026-08-28

The whole of `evaluation.py` was read directly from the Space
(`raw/main/evaluation.py`, public, no authentication) and its logic
re-implemented and executed. Everything above is confirmed. The results below
are **verified by execution**, not inferred, and they determine how we submit.

## Confirmed mechanics

- F-max thresholds are **the unique EPCR values present in our own submission**,
  swept descending — not a fixed grid. The comparison is `row.epcr >= t`.
- Consequence: both metrics depend only on the **order and tie structure** of our
  EPCRs, never on their magnitude. Any strictly increasing reparameterisation
  scores identically. Calibrating EPCR toward a "true probability" is worth
  exactly zero points; **injectivity is worth a great deal** (see below).
- Full match is `row.variants == true_variants` (frozenset equality), `break` on
  the first. Partial match is searched **only if no full match exists**, and only
  for a compound-het answer.
- `_parse_variant` applies `.strip()` to chrom and `.upper()` to ref/alt. There is
  **no chr-prefix normalisation** anywhere.

## The four results that decide our submission

| Submission shape | Rank points | F-max |
|---|---|---|
| True pair alone, one row | **100.0** | **1.0000** |
| True pair + 4 wrong rows at *lower* EPCR | **100.0** | **1.0000** |
| True pair tied with a wrong pair | 100.0 | 0.6667 |
| Same tie, wrong row first in file order | **50.0** | 0.6667 |
| True pair split into two single-variant rows | **50.0** | 1.0000 |
| One true + one wrong variant, rank 1 | 50.0 | 0.5000 |
| True pair, bare contig `15` instead of `chr15` | **0.0** | **0.0000** |

**1. Fill all ten rows.** Rows with EPCR strictly below the answer's cannot lower
either metric — F-max is a maximum over thresholds, and a lower threshold never
displaces the better one. Unused rows are forfeited upside, not saved risk.

**2. Never emit two equal EPCRs.** A tie at the answer's threshold pulls a wrong
row into the same prediction set: precision halves and F-max falls 1.0 → 0.667.
If file order puts the wrong row first it also takes rank 1, costing a further
50 rank points. Strict separation is order-preserving and therefore free.

**3. Never split a candidate pair across two rows.** Two single-variant rows can
still reach F-max 1.0, but `row.variants == true_variants` is then never true, so
the best available is partial credit: **100 → 50 rank points for nothing**.
A compound-het proposal is one row, using the `_2` columns.

**4. The chr prefix is worth the entire submission.** A correct answer emitted as
`15` rather than `chr15` scores 0.0 on both metrics while looking right to a human
reader. This is the single highest-consequence formatting detail in the project.

## Still not verifiable from here

The answer key (`gold_standard_track1.json`) lives in the private
`SageBio/mva-hackathon-2026-gt` dataset. Nothing above depends on its contents —
only on its *shape*, which `score_proband` fixes as a frozenset of 1 or 2
`(chrom, pos, ref, alt)` tuples.
