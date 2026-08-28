# ADR 0023 — Promotion runs before truncation, and the renderer is closed

**Status:** accepted
**Date:** 2026-08-28
**Touches:** ADR 0014 (EPCR injectivity), ADR 0015 (a pair outranks its single),
GP-15, GP-30, GP-32

## Context

An independent adversarial review of `mva.reporting.track1` found two defects that
cost points directly. Both were priced against the scorer, not argued about: the
replica in `tests/unit/test_track1_scoring.py` is transcribed from
`docs/references/track1-submission-contract.md` and is pinned by
`test_the_replica_reproduces_every_verified_scorer_result` to the seven submission
shapes the contract records as *executed* results.

### 1. Truncation ran before promotion

`build_submission_rows` sliced to ten rows and only then called
`_promote_pairs_above_subsets`:

```python
submitted, promoted = _promote_pairs_above_subsets(rows[:max_rows])
```

ADR 0015 exists because a single-variant candidate routinely outscores the pair it
was carved out of — the compound-het hypothesis carries a phase penalty its own
halves do not (GP-15) — and rank points go only to the best *full* match. But the
pair that outscores nothing is also the pair most likely to sit **below** the
ten-row cut, and slicing first deleted it before promotion could look for it.

Reproduced by execution (`test_the_promoted_pair_is_worth_100_rank_points_not_50`):
rank 1 is one half of the answer, ranks 2–10 are unrelated hypotheses, rank 11 is
the pair carrying both causal alleles.

| | rank points | F-max |
|---|---|---|
| truncate, then promote | **50.0** | **0.6667** |
| promote, then truncate | **100.0** | **1.0000** |

Both metrics, from the order of two statements. The rank half is the partial-credit
branch: a one-variant row can never satisfy `row.variants == true_variants`. The
F-max half is worse in kind — the second causal allele appeared in no submitted
row at all, so no threshold in the sweep could ever predict it.

### 2. `render_submission_csv` was an exported validation bypass

`SubmissionRow` is a plain dataclass with no invariants. The module's public
renderer serialised any sequence of them with no checks, so eleven rows, a bare
contig, a duplicate EPCR, an out-of-range EPCR and an arbitrary `proband_id` all
had a one-call public path to submission bytes. The reviewer confirmed
`write_submission` correctly rejects every one of those — which is the point: the
checks existed and the renderer stood beside them.

Separately, cross-row validation checked duplicate variant sets and EPCR ties but
never detected **one compound-het hypothesis split across two single-variant
rows**. Each row is independently valid, so `write_submission` accepted it. The
contract's executed result for that shape is 50 rank points where 100 was
available, at F-max 1.0000 — the F-max is what makes it silent.

## Decision

### Promotion runs on the full ranked list, before truncation

Composition order is now: `_drop_subsumed` → `_promote_pairs_above_subsets` →
truncate → `_enforce_epcr_separation`. Separation stays last and stays after
truncation: it is order-preserving by construction and only the submitted rows are
scored.

Reaching past the cut changes what promotion *costs*, so the pass is in three
parts and they are justified differently.

**1. Reordering inside the window is free.** No row enters or leaves. Every
superset already submitted is lifted above every subset it covers — ADR 0015
unchanged, including several pairs carrying one variant.

*Amended 2026-08-28.* "Free" is true of the lift in isolation and was read as
true of the pass as a whole, which it is not: step 2 brings a new superset into
the window, and that superset can land below a subset the lift had already
placed. An adversarial review reproduced exactly that — a five-row list into a
four-row window emitted a single above its own superset, costing the rank tier
this ADR exists to protect. The lift is therefore **re-run after the exchange**
(step 3). Once suffices: the lift only permutes the window, so it creates no new
exchange opportunity and there is no fixpoint to iterate to.

**2. Crossing the cut costs a row, so it is rationed.** A subset row with no
superset inside the window is **exchanged** with the highest-ranked superset
outside it: the pair takes the subset's slot, the subset takes the pair's. Exactly
one row enters, exactly one leaves, the row that leaves is the subset itself, and
**no unrelated row moves by a single position**.

**3. Re-lift, because step 2 can undo step 1.** The exchanged-in superset is
placed at the departing subset's rank, which may sit below a subset the lift
already ordered. Step 1 is therefore repeated. `validate_submission` then makes
the invariant unfalsifiable rather than trusted: a submission in which any row is
a strict subset of a lower-ranked row is **refused**, on every path to bytes —
`render_submission_csv`, `write_submission`, and the orchestrator's re-check —
not in the one renderer that happened to be reviewed.

### The exchange inverts one sentence of ADR 0015, on purpose

ADR 0015 keeps the outranked single *because keeping it is free* — a row below the
answer cannot lower either metric, verified by execution. Crossing the cut, it is
not free: the only way to keep it is to evict something else. And in the world the
challenge tells us we are in, a single-variant row is worth nothing —
`row.variants == true_variants` cannot hold for a one-variant row against a
two-variant answer — while the row it would displace is another *pair*, which can.
Its variant is not lost either: the promoted pair re-emits it, so the F-max union
at that threshold is unchanged.

**Keep the single when it is free; drop it when it costs a slot.**

Two alternatives were rejected, and the second is the one worth recording:

* *Slide the pair up and let the last row fall off the end.* Spends a submitted
  pair to keep a single that can never full-match.
* *Promote every out-of-window superset.* Measured on the golden case, this pulled
  four pairs above one single and evicted two **higher-scoring** unrelated pairs.
  That is a re-ranking on the strength of a subset relation, and re-ranking is a
  scientific judgement that does not belong in a rendering pass (the same
  reasoning ADR 0015 used to refuse implementing promotion as a score adjustment).

### The renderer validates; the raw path is named for what it skips

`render_submission_csv` now runs `validate_submission` and raises `ValueError`
listing every failing check. `render_submission_csv_unvalidated` keeps the raw
serialiser for the two callers that need it — tests that must feed the validator
bytes it should reject, and tests that price a hypothetical shape against the
scorer replica. It is in `mva.reporting.track1.__all__` and deliberately **not** in
`mva.reporting`'s package exports, so the shortest spelling of "serialise these
rows" is the one that checks them.

Three shapes were considered. *Stop exporting it* was rejected because a validator
has to be testable against invalid bytes and a test reaching for a private name is
a bypass with extra steps. *Rename it and leave it public* was rejected on its own,
because a frightening name is still the shortest path from a dataclass to a file.
*Make the safe spelling the obvious one* was adopted, which is both of the others
done together.

### Splitting a pair across two rows is a validation error

`_validate_split_pairs` rejects two single-variant rows the pipeline attributed to
**one gene** when no row proposes the two variants together.

The gene, not the coordinates: two variants are a compound heterozygote when they
sit in the same gene. Requiring every unordered pair of single-variant rows to be
joined would demand 45 pair rows for ten singles, which is neither possible nor a
claim we would want to make. Gene identity is read from `notes` — the twelve
columns the challenge specifies have no gene field, so a validator working on CSV
text has nothing else to key on — and a row with no attribution is skipped rather
than guessed at. `finding_type=secondary` is exempt: an incidental finding is
stated to be unrelated to the primary phenotype, and joining it would assert a
compound heterozygote the pipeline explicitly did not claim.

## Before / after on the golden case

The exchange rule was chosen partly *because* its diff is one row. Same synthetic
case, same ranked list, same fifteen candidates after `_drop_subsumed`:

| row | before | after |
|---|---|---|
| 5 | `chr11:5000000:A:GT` single, 0.5941 | `chr11:5000000:A:G / chr11:5000000:A:GT` pair, 0.3469 |

Rows 1–4 and 6–10 carry identical coordinates. The pair had ranked below the cut
and was not submitted at all; the single it subsumes is the row it replaced. EPCRs
below row 5 shift downward because `_enforce_epcr_separation` cascades from the
promoted row's own value — magnitude carries no information to the scorer, only
order and injectivity do (ADR 0014), and the order is unchanged.

**The top-ranked pair keeps rank 1 and all ten EPCRs stay distinct**, which is
the invariant this change had to preserve. It is proved on the synthetic case by
`tests/golden/test_golden_case.py::test_synthetic_causal_pair_ranks_first` and
`tests/unit/test_track1_composition.py::test_unrelated_rows_keep_their_epcr_order`;
the real case's own score and rank are patient-derived and stay on the case
volume.

## Consequences

- `tests/unit/test_track1_scoring.py` (17 tests) is new and asserts **points**,
  not layout, through a replica pinned to the contract's executed results.
- `_promote_pairs_above_subsets` returns a third value, the count of exchanges,
  logged as a count and never as a coordinate (GP-41).
- If the ground truth ever turns out to be single-variant, the exchange costs the
  full match that the retained single would have provided. That is ADR 0015's
  accepted risk, now taken one step further, and stated here so the reversal is
  cheap: delete the exchange step of `_promote_pairs_above_subsets`, keeping the
  lift.
- The split-pair check is best-effort on an optional column. A submission whose
  rows carry no `notes` is not checked for splits — recorded in
  `docs/tech-debt.md` rather than papered over.
