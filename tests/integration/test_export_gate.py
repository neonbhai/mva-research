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


# ---------------------------------------------------------------------------
# The gate must RUN, and PUBLIC must mean what the label says (M7).
#
# `_gate_public_artifact` was called for the submission and for nothing else,
# while MECHANISM_REPORT, DRUG_HYPOTHESES, REJECTION_RECORD and TRACK2_REPORT
# were all classified PUBLIC and all on the export allowlist. Four of the five
# artifacts a reader would publish were gated zero times, and one of them --
# `reports/track2_report.md` -- rendered the proband's coordinates:
#
#     reports/track2_report.md:49:  - Variants: GRCh38:chr15:40200000:C:T, ...
#
# `docs/architecture.md` says public export is "gated twice -- an allowlist plus a
# post-render content re-scan". That described the submission alone.
# ---------------------------------------------------------------------------


def test_every_public_artifact_passes_the_gate_during_the_run(
    synthetic_config: CaseConfig, synthetic_workspace: Workspace
) -> None:
    """Every PUBLIC artifact was gated as it was registered, and passed.

    `execute_pipeline` wires `_gate_public_artifact` to `RunContext.on_register`,
    which raises `ExportBlockedError` on a refusal -- so a completed run is itself
    the assertion that every artifact was gated. This test proves the hook is
    actually attached, and that the gate would have refused had it not been: it
    re-runs the identical verdict afterwards and requires an allowed decision for
    each one.
    """
    from mva.pipeline import RunContext
    from mva.privacy.export import gate_public_export

    gated: list[str] = []
    real_register = RunContext.register_artifact

    def spy(self: RunContext, **kwargs: object) -> object:
        artifact = real_register(self, **kwargs)  # type: ignore[arg-type]
        if artifact.sensitivity is Sensitivity.PUBLIC:
            gated.append(artifact.relative_path)
        return artifact

    RunContext.register_artifact = spy  # type: ignore[assignment,method-assign]
    try:
        result = execute_pipeline(synthetic_config, synthetic_workspace)
    finally:
        RunContext.register_artifact = real_register  # type: ignore[method-assign]

    public = [a for a in result.manifest.artifacts if a.sensitivity is Sensitivity.PUBLIC]
    assert public, "the run produced no public artifacts at all"
    assert {a.relative_path for a in public} == set(gated), (
        "a PUBLIC artifact was registered without passing through the gate hook"
    )

    for artifact in public:
        decision = gate_public_export(
            synthetic_workspace.root / artifact.relative_path,
            declared=artifact.sensitivity,
            allowlist=PUBLIC_EXPORT_ALLOWLIST,
        )
        assert decision.allowed, (
            f"{artifact.relative_path} is classified PUBLIC but the export gate "
            f"refuses it: {decision.reasons}. The run should not have completed."
        )


def test_no_public_artifact_other_than_the_submission_contains_a_variant_id(
    synthetic_config: CaseConfig, synthetic_workspace: Workspace
) -> None:
    """PUBLIC must mean 'carries no proband coordinate' -- one documented exception.

    The submission is that exception and cannot be otherwise: chrom/pos/ref/alt IS
    the Track 1 format. Every other PUBLIC artifact is a narrative about public
    gene/mechanism/drug knowledge, and a coordinate in one buys nothing the pair
    identifier does not already give. Two exceptions would be one too many: a label
    that sometimes means 'no patient data' stops being readable at a glance, and
    the reader glancing is the person about to attach the file to an email.
    """
    result = execute_pipeline(synthetic_config, synthetic_workspace)
    assert result.ranked_pairs, "no ranked pairs, so this test would assert nothing"

    for artifact in result.manifest.artifacts:
        if artifact.sensitivity is not Sensitivity.PUBLIC:
            continue
        if artifact.kind.value == "submission":
            continue
        text = (synthetic_workspace.root / artifact.relative_path).read_text(encoding="utf-8")
        for variant_id in result.ranked_pairs[0].variant_ids:
            assert variant_id not in text, (
                f"{artifact.relative_path} is classified PUBLIC and on the export "
                "allowlist, and it contains a proband variant coordinate. Either "
                "stop rendering it (preferred: name the pair by `pair_id`, which "
                "resolves in the SENSITIVE candidates artifact) or reclassify the "
                "artifact -- but do not leave a PUBLIC file carrying a genotype."
            )


def test_track2_report_still_names_the_pair_it_is_anchored_to(
    synthetic_config: CaseConfig, synthetic_workspace: Workspace
) -> None:
    """Removing the coordinates must not remove the anchoring.

    A drug hypothesis not tied to a ranked pair is not grounded in anything, so the
    fix had to keep the anchor while dropping the part of it that identifies the
    patient. `pair_id` is the join key into `candidates/ranked_pairs.json`.
    """
    result = execute_pipeline(synthetic_config, synthetic_workspace)
    track2 = next(a for a in result.manifest.artifacts if a.kind.value == "track2_report")
    text = (synthetic_workspace.root / track2.relative_path).read_text(encoding="utf-8")

    pair = result.ranked_pairs[0]
    assert pair.pair_id in text, "the report no longer names the pair it is anchored to"
    assert pair.gene_symbol in text
    assert "ranked_pairs.json" in text, (
        "the report drops the coordinates without telling the reader where the "
        "full record lives, which reads as an omission rather than a boundary"
    )


def test_the_gate_would_refuse_a_public_artifact_carrying_a_variant_id(
    tmp_path,  # type: ignore[no-untyped-def]
) -> None:
    """The second half of 'gated twice' has to be able to see this class of leak.

    When `track2_report.md` was rendering coordinates, running the gate over it
    returned `allowed=True`: no rule in the shared audit battery matches a
    colon-delimited `GRCh38:chr...:C:T`. `variant_id` is REDACTION_ONLY there for
    good reason -- the repo's own public knowledge tables are full of
    build-qualified coordinates -- but that reasoning does not reach a file about
    to leave the workspace, so the gate adds it back (`EXPORT_ONLY_RULES`).
    """
    from mva.privacy.export import gate_public_export

    path = tmp_path / "track2_report.md"
    path.write_text("# Track 2\n\n- Variants: GRCh38:chr15:40200000:C:T\n", encoding="utf-8")

    decision = gate_public_export(
        path, declared=Sensitivity.PUBLIC, allowlist=PUBLIC_EXPORT_ALLOWLIST
    )
    assert not decision.allowed
    assert "variant_id" in decision.scanned_rules
    # The refusal names the failing CHECK, never the coordinate (GP-41).
    assert "40200000" not in " ".join(decision.reasons)


def test_the_track1_submission_is_not_caught_by_the_export_only_rule(
    synthetic_config: CaseConfig, synthetic_workspace: Workspace
) -> None:
    """The one artifact that must carry coordinates still passes.

    The submission writes chrom/pos/ref/alt as four separate CSV columns, never as
    this project's `build:chrom:pos:ref:alt` rendering, so `EXPORT_ONLY_RULES`
    tightens the gate without needing an exception carved out for it. If that ever
    changes, this test fails before the run does.
    """
    from mva.privacy.export import EXPORT_ONLY_RULES

    result = execute_pipeline(synthetic_config, synthetic_workspace)
    submission = next(a for a in result.manifest.artifacts if a.kind.value == "submission")
    data = (synthetic_workspace.root / submission.relative_path).read_bytes()
    for rule in EXPORT_ONLY_RULES:
        assert rule.pattern.search(data) is None, (
            f"the submission now matches export-only rule {rule.rule_id!r}. The gate "
            "will refuse it, blocking every run. Either the submission format "
            "changed or the rule is too broad -- decide which, deliberately."
        )
