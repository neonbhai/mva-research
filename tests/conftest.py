"""Shared fixtures.

The synthetic case is copied into a `tmp_path` workspace for every test that runs
the pipeline, so no test ever writes into the repository and no test depends on a
previous run's artifacts.
"""

from __future__ import annotations

import shutil
from collections.abc import Callable, Iterator
from pathlib import Path

import pytest

from mva.cli import load_case_config_with_defaults
from mva.config import CaseConfig, Workspace, find_repo_root

REPO_ROOT = find_repo_root(Path(__file__))
FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "synthetic"


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return REPO_ROOT


@pytest.fixture
def synthetic_config() -> CaseConfig:
    """The demo case config, with shared defaults merged beneath it."""
    return load_case_config_with_defaults(
        REPO_ROOT / "config" / "synthetic-case.yaml",
        REPO_ROOT / "config" / "default.yaml",
    )


@pytest.fixture
def synthetic_workspace(tmp_path: Path) -> Iterator[Workspace]:
    """A throwaway workspace seeded with the synthetic inputs.

    `allow_inside_repo` is not needed here: `tmp_path` is outside the repo, which
    is exactly the arrangement the privacy model requires of a real run.
    """
    inputs = tmp_path / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    shutil.copy(FIXTURE_DIR / "synthetic_case.vcf", inputs / "synthetic_case.vcf")
    shutil.copy(FIXTURE_DIR / "synthetic_phenotype.tsv", inputs / "synthetic_phenotype.tsv")
    yield Workspace(root=tmp_path.resolve(), repo_root=REPO_ROOT)


def _read_golden_tsv(name: str) -> list[dict[str, str]]:
    """Read a golden expectation file, skipping comment lines."""
    import csv

    path = REPO_ROOT / "tests" / "golden" / name
    lines = [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    return list(csv.DictReader(lines, delimiter="\t"))


@pytest.fixture(scope="session")
def golden() -> Callable[[str], list[dict[str, str]]]:
    """Reader for `tests/golden/*.tsv`.

    Exposed as a fixture because `tests/` is deliberately not an importable
    package -- test modules should not import from each other.
    """
    return _read_golden_tsv


# --------------------------------------------------------------------------- network


#: Setting this to ``1`` opts the suite into the tests that reach the internet.
NETWORK_OPT_IN_ENV = "MVA_VERIFY_CITATIONS"


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Skip `network`-marked tests unless someone explicitly asked for them.

    PRIV-05 makes "offline by default" a property of this repository rather than a
    habit, so a test that opens a socket has to be opted into twice over: it must
    be marked, and the run must ask for it. Either ask works --
    ``pytest -m network`` (the marker expression) or ``MVA_VERIFY_CITATIONS=1``
    (the environment) -- because the first is what a human types and the second is
    what a scheduled refresh job sets.

    ``-m "not network"`` is honoured as the opt-*out* it obviously is, rather than
    being read as the substring "network" and treated as consent.
    """
    import os

    markexpr = str(config.getoption("markexpr", default="") or "")
    asked_by_marker = "network" in markexpr and "not network" not in markexpr
    asked_by_env = os.environ.get(NETWORK_OPT_IN_ENV) == "1"
    if asked_by_marker or asked_by_env:
        return

    skip = pytest.mark.skip(
        reason=(
            "network test: opt in with `-m network` or "
            f"`{NETWORK_OPT_IN_ENV}=1`, or run `just verify-citations`"
        )
    )
    for item in items:
        if "network" in item.keywords:
            item.add_marker(skip)
