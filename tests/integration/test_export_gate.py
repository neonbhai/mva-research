"""The public-export gate must actually run on the artifacts that leave.

Regression guard for an integration hole found in Phase 4 review: the orchestrator
originally wrote the submission with `render_submission_csv` + a raw file write,
which bypassed `write_submission` and therefore bypassed both the contract
self-check and the GP-43 export gate. The submission is the one artifact that
necessarily leaves the workspace, so it is exactly the one that must not be
trusted on its declared classification alone.
"""

from __future__ import annotations

import pytest

from mva.config import CaseConfig, Workspace
from mva.models.base import Sensitivity
from mva.orchestrator import PUBLIC_EXPORT_ALLOWLIST, execute_pipeline

pytestmark = pytest.mark.integration


def test_every_public_artifact_is_on_the_export_allowlist(
    synthetic_config: CaseConfig, synthetic_workspace: Workspace
) -> None:
    """Deny by default: a PUBLIC artifact not on the allowlist is a gap."""
    result = execute_pipeline(synthetic_config, synthetic_workspace)
    public = [a for a in result.manifest.artifacts if a.sensitivity is Sensitivity.PUBLIC]
    assert public, "the run produced no public artifacts at all"
    for artifact in public:
        name = artifact.relative_path.rsplit("/", 1)[-1]
        assert name in PUBLIC_EXPORT_ALLOWLIST, (
            f"{name!r} is classified PUBLIC but is not on PUBLIC_EXPORT_ALLOWLIST. "
            "Add it deliberately or reclassify the artifact; classification alone "
            "must never be enough to publish (GP-43)."
        )


def test_sensitive_artifacts_are_never_classified_public(
    synthetic_config: CaseConfig, synthetic_workspace: Workspace
) -> None:
    """Variant-level artifacts carry genotypes and must stay SENSITIVE."""
    result = execute_pipeline(synthetic_config, synthetic_workspace)
    must_be_sensitive = {
        "normalised_variants",
        "annotated_variants",
        "candidate_pairs",
        "evidence_db",
        "dossier",
    }
    for artifact in result.manifest.artifacts:
        if artifact.kind.value in must_be_sensitive:
            assert artifact.sensitivity is Sensitivity.SENSITIVE, (
                f"{artifact.kind.value} is classified {artifact.sensitivity.value}; "
                "it carries patient genotypes and must be SENSITIVE"
            )


def test_submission_passes_its_own_contract_check(
    synthetic_config: CaseConfig, synthetic_workspace: Workspace
) -> None:
    """The written submission re-validates against the Track 1 contract."""
    from mva.reporting.track1 import validate_submission

    result = execute_pipeline(synthetic_config, synthetic_workspace)
    submission = next(a for a in result.manifest.artifacts if a.kind.value == "submission")
    text = (synthetic_workspace.root / submission.relative_path).read_text(encoding="utf-8")
    ok, errors = validate_submission(text)
    assert ok, f"written submission violates the Track 1 contract: {errors}"


def test_gate_refuses_a_sensitive_artifact(tmp_path) -> None:  # type: ignore[no-untyped-def]
    """The gate refuses on declared sensitivity even for an allowlisted name."""
    from mva.privacy.export import gate_public_export

    path = tmp_path / "track1_submission.csv"
    path.write_text("proband_id\nPROBAND01\n", encoding="utf-8")

    decision = gate_public_export(
        path, declared=Sensitivity.SENSITIVE, allowlist=PUBLIC_EXPORT_ALLOWLIST
    )
    assert not decision.allowed
    assert any("sensitiv" in reason.lower() for reason in decision.reasons)
