"""The reimplemented scorer is checked against the challenge's own code.

``tests/unit/test_track1_scoring.py`` reimplements the challenge scorer so that
the submission renderer can be reasoned about offline. That reimplementation is
the thing this repository's ranking strategy is built on: EPCR injectivity
(ADR 0014), composition ordering (ADR 0015) and promotion-before-truncation
(ADR 0023) all exist to satisfy rules that live in someone else's file.

Thirteen tests assert points against that reimplementation. Every one of them
verifies our copy against our copy, so all thirteen would stay green if the real
scorer's behaviour moved and ours did not. That is the same failure the privacy
audit had -- a check that cannot tell you it has stopped checking -- and it gets
the same treatment: an oracle, and a differential run against it.

The oracle is ``tests/fixtures/scorer/evaluation_vendored.py``, a byte-for-byte
copy of the published scorer. Its provenance and the reason it is a copy rather
than a hash are in ``tests/fixtures/scorer/SOURCE.md``.
"""

from __future__ import annotations

import hashlib
import importlib.util
import random
import sys
from pathlib import Path
from types import ModuleType

import pytest

# No package __init__ under tests/, so the sibling module is loaded by path
# rather than imported. Keeping it path-based means this file does not depend on
# how pytest happens to have configured sys.path.
_SIBLING = Path(__file__).resolve().parent / "test_track1_scoring.py"
_spec = importlib.util.spec_from_file_location("mva_track1_scoring_ref", _SIBLING)
assert _spec is not None and _spec.loader is not None
_ref = importlib.util.module_from_spec(_spec)
sys.modules["mva_track1_scoring_ref"] = _ref
_spec.loader.exec_module(_ref)
RANK_TIERS = _ref.RANK_TIERS
f_max = _ref.f_max
rank_points = _ref.rank_points

FIXTURE_DIR = Path(__file__).resolve().parents[1] / "fixtures" / "scorer"
VENDORED = FIXTURE_DIR / "evaluation_vendored.py"

#: sha256 of the vendored scorer, recorded in SOURCE.md at fetch time. This does
#: not prove the upstream file is unchanged -- nothing offline can -- it proves
#: THIS file is the one whose provenance SOURCE.md describes.
VENDORED_SHA256 = "6d18b581e65a45e1ccc120071d588e740c2e42e983ff50704c60a40232b19180"

#: Fabricated loci. Three contigs x three positions is enough to exercise full
#: match, partial match, no match, and rows that overlap the truth on one allele.
_VARIANTS: list[tuple[str, int, str, str]] = [
    (f"chr{contig}", pos, "A", "G") for contig in (1, 2, 15) for pos in (100, 200, 300)
]


def _load_oracle() -> ModuleType:
    """Import the vendored scorer under a real module name.

    The name must be registered in ``sys.modules`` before ``exec_module`` runs:
    the file uses ``@dataclass``, and dataclasses resolves annotations through
    ``sys.modules[cls.__module__]``, which raises for a module that is not there.
    """
    spec = importlib.util.spec_from_file_location("mva_vendored_scorer", VENDORED)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["mva_vendored_scorer"] = module
    spec.loader.exec_module(module)
    return module


@pytest.mark.unit
def test_vendored_scorer_matches_its_recorded_digest() -> None:
    """The oracle is the file SOURCE.md documents, not whatever landed here."""
    digest = hashlib.sha256(VENDORED.read_bytes()).hexdigest()
    assert digest == VENDORED_SHA256, (
        "the vendored scorer does not match the digest recorded in SOURCE.md. "
        "If it was re-fetched deliberately, update both and re-read the diff: "
        "a changed scorer can silently invalidate ADR 0014, ADR 0015 and ADR 0023."
    )


@pytest.mark.unit
def test_the_rank_tiers_are_the_published_ones() -> None:
    oracle = _load_oracle()
    assert [tuple(tier) for tier in oracle.RANK_POINT_TIERS] == [tuple(t) for t in RANK_TIERS]


def _random_case(rng: random.Random) -> tuple[list[tuple[frozenset, float]], frozenset]:
    """One submission and one ground truth.

    Truth is one or two variants, because those are the two shapes the scorer
    branches on: partial credit exists only for a compound-het answer.
    """
    truth = frozenset(rng.sample(_VARIANTS, rng.choice([1, 2])))
    rows: list[tuple[frozenset, float]] = []
    for i in range(rng.randint(1, 7)):
        variants = frozenset(rng.sample(_VARIANTS, rng.choice([1, 1, 2])))
        epcr = round(1.0 - i * 0.09 - rng.random() * 0.02, 4)
        if epcr > 0:
            rows.append((variants, epcr))
    return rows, truth


@pytest.mark.unit
def test_both_metrics_agree_with_the_published_scorer_over_random_submissions() -> None:
    """Differential test. A divergence names the case rather than just failing.

    The seed is fixed so a failure is reproducible; the point of randomising at
    all is to reach the branch combinations a hand-written table does not --
    partial credit under a compound-het truth, a full match below a partial one,
    and F-max thresholds that coincide.
    """
    oracle = _load_oracle()
    # S311: not cryptography. A fixed-seed PRNG is the point -- the same 2,000
    # cases every run, so a divergence is reproducible rather than flaky.
    rng = random.Random(20260828)  # noqa: S311

    for trial in range(2000):
        rows, truth = _random_case(rng)
        if not rows:
            continue
        # The scorer ranks by EPCR descending, ties broken by file order, and it
        # does that itself -- rank is assigned here the same way so the oracle
        # sees the ordering it would have computed.
        order = sorted(range(len(rows)), key=lambda i: (-rows[i][1], i))
        oracle_rows = [
            oracle.SubmissionRow(variants=rows[i][0], epcr=rows[i][1], rank=rank)
            for rank, i in enumerate(order, start=1)
        ]
        expected = oracle.score_proband("PROBAND01", oracle_rows, truth)

        assert rank_points(rows, truth) == pytest.approx(expected.rank_points), (
            f"rank_points diverged on trial {trial}: {len(truth)}-variant truth, {len(rows)} row(s)"
        )
        assert f_max(rows, truth) == pytest.approx(expected.f_max), (
            f"f_max diverged on trial {trial}: {len(truth)}-variant truth, {len(rows)} row(s)"
        )
