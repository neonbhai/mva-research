"""Shared fixtures.

The synthetic case is copied into a `tmp_path` workspace for every test that runs
the pipeline, so no test ever writes into the repository and no test depends on a
previous run's artifacts.
"""

from __future__ import annotations

import shutil
from collections.abc import Iterator
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


def read_golden_tsv(name: str) -> list[dict[str, str]]:
    """Read a golden expectation file, skipping comment lines."""
    import csv

    path = REPO_ROOT / "tests" / "golden" / name
    lines = [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    return list(csv.DictReader(lines, delimiter="\t"))
