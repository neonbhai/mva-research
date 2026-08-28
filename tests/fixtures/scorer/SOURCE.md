# Vendored Track 1 scorer

`evaluation_vendored.py` is a byte-for-byte copy of the challenge's own scoring
code, kept here as a **differential oracle**. It is not run in production and
nothing in `src/` imports it.

| | |
|---|---|
| Source | `https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/raw/main/evaluation.py` |
| Space revision at fetch | `1112710080520b3a2848d11b6ce1327dafbe79cf` |
| File sha256 | `6d18b581e65a45e1ccc120071d588e740c2e42e983ff50704c60a40232b19180` |
| Fetched | 2026-08-28 |

## Why a copy and not a hash

Pinning a hash tells you the upstream file changed. It does not tell you whether
your own understanding of the scoring rules is still correct, and that is the
thing that can silently rot: `tests/unit/test_track1_scoring.py` **reimplements**
the scorer, and 13 tests assert points against that reimplementation. Those tests
verify our copy against our copy. They would stay green if the real scorer's
behaviour moved.

So the copy is used as an oracle. `test_scorer_differential.py` runs both
implementations over randomised submissions and asserts they agree on both
metrics. A divergence fails the suite and names the case.

The Space is public and mutable, and the challenge runs until 2026-10-24. Re-fetch
and re-run before a late submission; `test_vendored_scorer_matches_its_recorded_digest`
fails loudly if this file is edited, but it is offline and cannot know that
upstream has moved. That check is deliberately a human step, not a network call
in the test suite (PRIV-05: no network on the patient-data path, and a test that
reaches the internet is a test that fails on a plane).
