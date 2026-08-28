# Submission ledger

Six Track 1 attempts exist, one Track 2 attempt exists, and only the best Track 1
score counts. All seven slots are written down before they are used, because an
attempt budget that is only tracked in someone's memory gets spent reactively —
and Track 2 has no second chance at all.

**What is deliberately not here.** No variant coordinates, no allele
representations, no EPCR values, no gene-level result. Those are the child's
data and live on the case volume (GP-40). Scores that the public leaderboard
already displays are recorded; the confidence numbers behind them are not.

---

## Track 1 — 6 attempts, best score counts

| # | Date | Approach | Rank pts | F-max | Repo revision | Notes |
|---|---|---|---|---|---|---|
| 1 | 2026-08-28 | `provenance-first-acmg-comphet` | **100.0** | **1.000** | `576a8e2` | Full match at rank 1. Ranking produced by the scripted pipeline, not the library — stated in `track1-method.md`. |
| 2 | — | *unspent* | | | | Reserved for the library's own ranking, once it has completed a real-callset run end to end. |
| 3 | — | *unspent* | | | | |
| 4 | — | *unspent* | | | | |
| 5 | — | *unspent* | | | | |
| 6 | — | *unspent* | | | | |

**Attempt 1 detail.** Submitted under HF account `lordneo7`. Ten rows, the
candidate pair carried as a single row using the `_2` columns rather than split
across two, every contig `chr`-prefixed, EPCRs distinct and descending within
`(0, 1]`. The scorer returned a full match at rank 1 and reported an F-max
threshold equal to the top row's own EPCR, meaning the winning prediction set
was that row alone.

**Why attempt 2 is being held rather than spent.** The automated metric is at
its ceiling and cannot be improved on — see below — so a resubmission can only
change what a human reviewer reads. The one change worth spending a slot on is a
ranking the library produced itself, which is also the blind run described in
`docs/track1-method.md`. Until that run completes there is nothing a new
submission would say that this one does not.

## Track 2 — 1 attempt, no resubmission

| # | Date | Status | Notes |
|---|---|---|---|
| 1 | — | *unspent* | Requires report, public GitHub URL and a 3-minute video URL. One shot: everything is verified before it is sent, not after. |

---

## The automated metric is not the contest

As of 2026-08-27 the public Track 1 leaderboard showed **every scored entry at
100.0 rank points and F-max 1.000** — the mathematical ceiling — with the
earliest such entry roughly 22 hours after the data was released, and most on
the submitter's first attempt. That was read directly off the public leaderboard
tab, not inferred.

Two things follow, and both are reasons to stop optimising the CSV:

1. A perfect automated score is a **participation gate**, not a differentiator.
   Attempt 1 clearing it confirms the submission format is correct and says
   nothing else about the work.
2. Every remaining distinction is decided by the qualitative review window and by
   Track 2. Effort spent moving a number that is already at its maximum is
   effort not spent on the parts that are still open.

This is recorded here rather than in a chat log because the temptation to burn
the remaining five attempts on a score that cannot move is exactly what a ledger
is for.

## Rules for spending a slot

- The predictions file passes the structural pre-flight in
  `docs/submission-runbook.md` §3 before upload, every time, including on
  attempt 6.
- The filename carries the account name and a short approach name — that string
  is what the public leaderboard displays as the model.
- The repo revision that produced the attempt is recorded in the table above,
  so a score can always be traced to the code that earned it.
- The submitted CSV and report stay on the case volume. Copies made to a plain
  disk for a file-picker's sake are deleted afterwards: they sit outside the
  encrypted volume and therefore outside the cryptographic-erasure guarantee
  that `docs/privacy-model.md` relies on.
