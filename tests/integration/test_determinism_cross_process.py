"""Determinism verified ACROSS processes, not within one.

`mva verify determinism` calls `execute_pipeline` twice in a single interpreter:
same hash seed, same timezone, same already-imported modules. It therefore cannot
observe the class of bug it exists to catch. A reproducibility review made the
point concretely, and separately confirmed by experiment that the pipeline does
survive `PYTHONHASHSEED` and `TZ` variation — this test is what keeps that true.

Marked `slow`: it spawns two full pipeline runs.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from mva.config import find_repo_root

pytestmark = [pytest.mark.integration, pytest.mark.slow]

REPO = find_repo_root(Path(__file__))
FIXTURES = REPO / "tests" / "fixtures" / "synthetic"

#: Artifacts excluded from byte comparison, with the reason each is excluded.
#: `provenance.json` is compared FIELD-WISE below instead of being skipped, so
#: the exclusion cannot hide a real difference.
_BINARY_CONTAINERS = (".duckdb", ".duckdb.wal")


def _seed(workspace: Path) -> None:
    (workspace / "inputs").mkdir(parents=True, exist_ok=True)
    for name in ("synthetic_case.vcf", "synthetic_phenotype.tsv"):
        shutil.copy(FIXTURES / name, workspace / "inputs" / name)


def _run(workspace: Path, env_overrides: dict[str, str]) -> None:
    env = {**os.environ, **env_overrides}
    result = subprocess.run(  # noqa: S603
        [
            sys.executable,
            "-m",
            "mva.cli",
            "run",
            "all",
            "--config",
            str(REPO / "config" / "synthetic-case.yaml"),
            "--defaults",
            str(REPO / "config" / "default.yaml"),
            "--workspace",
            str(workspace),
        ],
        cwd=REPO,
        env=env,
        capture_output=True,
        text=True,
        check=False,
        timeout=600,
    )
    assert result.returncode == 0, (
        f"pipeline failed (exit {result.returncode}). stderr tail:\n"
        + "\n".join(result.stderr.splitlines()[-15:])
    )


def _digest(workspace: Path) -> dict[str, str]:
    digests: dict[str, str] = {}
    runs = workspace / "runs"
    for path in sorted(runs.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(runs).as_posix()
        # Strip the run-id directory so two workspaces are comparable.
        relative = relative.split("/", 1)[1] if "/" in relative else relative
        if relative.endswith(_BINARY_CONTAINERS):
            continue
        digests[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
    return digests


def test_pipeline_is_deterministic_across_processes(tmp_path: Path) -> None:
    """Two separate interpreters, different hash seeds and timezones, same bytes."""
    first = tmp_path / "run-a"
    second = tmp_path / "run-b"
    _seed(first)
    _seed(second)

    _run(first, {"PYTHONHASHSEED": "1", "TZ": "UTC", "LC_ALL": "C"})
    _run(second, {"PYTHONHASHSEED": "2", "TZ": "Asia/Kolkata", "LC_ALL": "en_US.UTF-8"})

    a, b = _digest(first), _digest(second)
    assert a, "first run produced no artifacts"
    assert set(a) == set(b), (
        "the two runs produced different artifact sets: "
        f"only in A={sorted(set(a) - set(b))} only in B={sorted(set(b) - set(a))}"
    )

    differing = [name for name in sorted(a) if a[name] != b[name]]
    # provenance.json legitimately differs on completed_at and git_dirty; it is
    # compared field-wise below rather than skipped outright.
    differing = [name for name in differing if not name.endswith("provenance.json")]
    assert not differing, (
        "artifacts differ between processes with different PYTHONHASHSEED/TZ:\n  "
        + "\n  ".join(differing)
        + "\n\nGP-30 remediation: find the unordered iteration, unstable sort, "
        "locale-dependent format or environment read that leaked into output."
    )


def test_provenance_differs_only_in_declared_volatile_fields(tmp_path: Path) -> None:
    """Compare the manifest field-wise instead of excluding the whole file."""
    first = tmp_path / "p-a"
    second = tmp_path / "p-b"
    _seed(first)
    _seed(second)
    _run(first, {"PYTHONHASHSEED": "1"})
    _run(second, {"PYTHONHASHSEED": "3"})

    def _load(workspace: Path) -> dict[str, object]:
        path = next((workspace / "runs").rglob("provenance.json"))
        parsed = json.loads(path.read_text(encoding="utf-8"))
        assert isinstance(parsed, dict)
        return parsed

    a, b = _load(first), _load(second)

    #: Fields allowed to differ between two runs of the same inputs.
    volatile = {"completed_at", "started_at", "git_dirty", "artifacts"}
    for key in sorted(set(a) | set(b)):
        if key in volatile:
            continue
        assert a.get(key) == b.get(key), (
            f"provenance field {key!r} differs between two runs of identical inputs; "
            "it is not in the declared volatile set, so this is a determinism bug"
        )

    assert a["config_hash"] == b["config_hash"]
    assert a["run_id"] == b["run_id"]
    assert a["reference_versions"] == b["reference_versions"]
