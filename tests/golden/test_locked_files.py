"""Hash locks on the files that silently change results.

GP-32 says a scoring-weight change requires a decision record, a test and a
before/after comparison, and that golden expectations are never silently
re-baselined. A reproducibility review demonstrated that neither was enforced:
setting `rarity: 0.18 -> 0.00` and flipping two verdicts in
`expected_drug_outcomes.tsv` both left `just verify` fully green.

A promise nothing checks is not a promise. These locks make the two most
result-altering edits in the repository impossible to make *quietly*: the diff
must include the new hash, which is the moment a reviewer asks for the ADR.

**Updating a hash here is not a chore to work around.** It is the gate. If you
are changing one, you should already have written the decision record.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from mva.config import find_repo_root

pytestmark = [pytest.mark.golden, pytest.mark.unit]

REPO = find_repo_root(Path(__file__))

#: path -> sha256. Regenerate with:
#:     shasum -a 256 config/default.yaml tests/golden/*.tsv
LOCKED: dict[str, str] = {
    # Updated by ADR 0010 (frequency.min_allele_number). Golden expectations
    # below are unchanged: the synthetic case ranks identically before and after.
    "config/default.yaml": "bba9451138e1b3eb38d3622a0ab1911b3b599d10e3afbd314e28e343d8ff8058",
    "tests/golden/expected_ranking.tsv": "ec082f5843ed0f86b37aaa370d48b7857497639cc0a253b5eaf87be153cd802b",  # noqa: E501 - a sha256 does not wrap
    "tests/golden/expected_drug_outcomes.tsv": "3a8184e2379c69b9c3b04c95a4ba0bba1866b7d68b5afe6e11787acd024d6d8a",  # noqa: E501 - a sha256 does not wrap
}

REMEDIATION = """

GP-32 remediation: this file changes pipeline results, so it is hash-locked.
If you meant to change it:
  1. write a decision record under docs/decisions/ saying what changed and why;
  2. record the before/after comparison (ranking and drug triage) in that ADR;
  3. update the hash below in the SAME commit.
If you did NOT mean to change it, revert the edit. Do not update the hash to
make this test pass — that is precisely the silent re-baseline the lock exists
to prevent.
"""


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.mark.parametrize("relative", sorted(LOCKED))
def test_result_altering_file_is_unchanged(relative: str) -> None:
    path = REPO / relative
    assert path.is_file(), f"locked file {relative} is missing"
    actual = _sha256(path)
    expected = LOCKED[relative]
    assert expected, f"No hash recorded for {relative}. Set it to {actual!r}." + REMEDIATION
    assert actual == expected, (
        f"{relative} changed.\n  expected sha256 {expected}\n  actual   sha256 {actual}"
        + REMEDIATION
    )


def test_scoring_weights_still_sum_to_one() -> None:
    """A weight edit that keeps the sum valid is the dangerous kind."""
    import yaml

    raw = yaml.safe_load((REPO / "config" / "default.yaml").read_text(encoding="utf-8"))
    weights = raw["weights"]
    positive = sum(value for key, value in weights.items() if key != "contradiction_penalty_weight")
    assert abs(positive - 1.0) < 1e-9, (
        f"positive scoring weights sum to {positive}, not 1.0" + REMEDIATION
    )
