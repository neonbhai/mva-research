"""Provenance must distinguish runs that produced different answers.

Regression guard for a reproducibility-review finding: two runs over *different
reference data* were indistinguishable. Same case config produced the same
`config_hash`, the same `run_id`, the same input hashes and the same run
directory — while emitting a different submission. They silently overwrote each
other. A manifest that cannot tell those apart is not a reproduction record.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from mva.config import CaseConfig, Workspace, find_repo_root
from mva.orchestrator import execute_pipeline
from mva.pipeline import make_run_id, reference_versions_from_manifest

pytestmark = pytest.mark.integration

REPO = find_repo_root(Path(__file__))
MANIFEST = REPO / "knowledge" / "manifests" / "knowledge.yaml"


def test_reference_versions_are_recorded(
    synthetic_config: CaseConfig, synthetic_workspace: Workspace
) -> None:
    """Every knowledge table the run consumed appears in the manifest."""
    result = execute_pipeline(synthetic_config, synthetic_workspace)
    refs = result.manifest.reference_versions
    assert refs, (
        "RunManifest.reference_versions is empty. The field is declared, persisted "
        "to the evidence store and rendered — but if nothing populates it, a run "
        "cannot be told apart from one using different reference data."
    )
    names = {r.name for r in refs}
    assert {"consequences", "frequencies"} <= names, (
        f"annotation tables missing from reference_versions: {sorted(names)}"
    )
    for ref in refs:
        assert ref.checksum, f"reference {ref.name!r} recorded without a content hash"


def test_run_id_changes_when_reference_data_changes(tmp_path: Path) -> None:
    """Different reference content must not collide in the same run directory."""
    from mva.cli import load_case_config_with_defaults

    config = load_case_config_with_defaults(
        REPO / "config" / "synthetic-case.yaml", REPO / "config" / "default.yaml"
    )

    real_versions = reference_versions_from_manifest(MANIFEST)
    assert real_versions, "no reference versions parsed from the manifest"

    # A manifest whose content hashes differ must yield a different run id, even
    # though the case config is byte-identical.
    altered = tmp_path / "knowledge.yaml"
    text = MANIFEST.read_text(encoding="utf-8")
    marker = next(r.checksum for r in real_versions if r.checksum)
    assert marker
    altered.write_text(text.replace(marker, "0" * len(marker), 1), encoding="utf-8")
    altered_versions = reference_versions_from_manifest(altered)

    from mva.pipeline import reference_fingerprint

    baseline = make_run_id(config, reference_fingerprint_value=reference_fingerprint(real_versions))
    changed = make_run_id(
        config, reference_fingerprint_value=reference_fingerprint(altered_versions)
    )
    assert baseline != changed, (
        "the run id is identical for two different sets of reference data; the two "
        "runs would share a directory and silently overwrite each other"
    )


def test_manifest_records_git_state_and_tool_versions(
    synthetic_config: CaseConfig, synthetic_workspace: Workspace
) -> None:
    """A run must be traceable to the code that produced it."""
    result = execute_pipeline(synthetic_config, synthetic_workspace)
    manifest = result.manifest
    assert manifest.tool_versions, "no tool versions recorded"
    assert any(t.name == "python" for t in manifest.tool_versions)
    assert manifest.config_hash
    assert manifest.inputs, "no input hashes recorded"
    for record in manifest.inputs:
        assert record.content_hash
        assert not Path(record.relative_path).is_absolute(), (
            "an input path is absolute; provenance must not record where patient data lives on disk"
        )


def test_synthetic_workspace_seeded(tmp_path: Path) -> None:
    """Sanity: the fixture copy helper behaves (guards the tests above)."""
    src = REPO / "tests" / "fixtures" / "synthetic" / "synthetic_case.vcf"
    dst = tmp_path / "synthetic_case.vcf"
    shutil.copy(src, dst)
    assert "##mva_synthetic=true" in dst.read_text(encoding="utf-8")
