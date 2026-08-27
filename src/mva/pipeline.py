"""Composition root.

This is the only module that knows the whole graph (GP-03). Stages never import
each other; they are wired together here and communicate through the typed models
in `mva.models`.

Every stage writes an artifact and registers `ArtifactProvenance` for it, so the
run manifest is a complete account of what was produced from what. Artifacts are
written under the external workspace, never inside the repo (ADR 0006).
"""

from __future__ import annotations

import json
import platform
import subprocess
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

import yaml

from mva.clock import Clock, SystemClock, demo_clock
from mva.config import CaseConfig, NetworkProfile, Workspace
from mva.determinism import canonical_json, hash_file, hash_text, short_hash
from mva.errors import ConfigError
from mva.models.base import Sensitivity
from mva.models.provenance import (
    ArtifactKind,
    ArtifactProvenance,
    CommandRecord,
    InputRecord,
    ReferenceVersion,
    RunManifest,
    ToolVersion,
)

# ---------------------------------------------------------------------------
# Artifact registry
# ---------------------------------------------------------------------------

#: Which artifacts are permitted to leave the workspace. Everything absent from
#: this mapping defaults to SENSITIVE — fail closed (GP-43).
ARTIFACT_SENSITIVITY: dict[ArtifactKind, Sensitivity] = {
    ArtifactKind.INPUT_VCF: Sensitivity.SENSITIVE,
    ArtifactKind.INPUT_PHENOTYPE: Sensitivity.SENSITIVE,
    ArtifactKind.NORMALISED_VARIANTS: Sensitivity.SENSITIVE,
    ArtifactKind.ANNOTATED_VARIANTS: Sensitivity.SENSITIVE,
    ArtifactKind.CANDIDATE_PAIRS: Sensitivity.SENSITIVE,
    ArtifactKind.EVIDENCE_DB: Sensitivity.SENSITIVE,
    ArtifactKind.DOSSIER: Sensitivity.SENSITIVE,
    # Aggregate counts only; no genotypes.
    ArtifactKind.QC_REPORT: Sensitivity.DERIVED_SAFE,
    ArtifactKind.RUN_MANIFEST: Sensitivity.DERIVED_SAFE,
    ArtifactKind.PROVENANCE_MANIFEST: Sensitivity.DERIVED_SAFE,
    ArtifactKind.PRIVACY_AUDIT: Sensitivity.DERIVED_SAFE,
    # The submission is the ONE artifact permitted to carry proband coordinates
    # while classified PUBLIC: that IS the challenge format, and there is no way
    # to submit without them. It writes chrom/pos/ref/alt as separate CSV columns
    # and passes the export gate, which re-scans the bytes.
    ArtifactKind.SUBMISSION: Sensitivity.PUBLIC,
    # Track 2 outputs concern public gene/mechanism/drug knowledge, not the
    # proband's genotypes.
    #
    # That sentence used to be false for TRACK2_REPORT, which rendered
    # `pair.variant_ids` into its "Anchoring variant pair" section — a
    # build-qualified coordinate and both alleles, for a real proband, in a file
    # classified PUBLIC and on the export allowlist. The report now names the pair
    # by `pair_id` only; the coordinates stay in `candidate_pairs.json`, which is
    # SENSITIVE and never leaves the workspace. See
    # `mva.reporting.track2._pair_context`.
    ArtifactKind.MECHANISM_REPORT: Sensitivity.PUBLIC,
    ArtifactKind.DRUG_HYPOTHESES: Sensitivity.PUBLIC,
    ArtifactKind.REJECTION_RECORD: Sensitivity.PUBLIC,
    ArtifactKind.TRACK2_REPORT: Sensitivity.PUBLIC,
}


@dataclass
class RunContext:
    """Mutable state threaded through a single pipeline run.

    Deliberately the only mutable object in the data path: stages consume and
    return frozen models, and this accumulates the bookkeeping around them.
    """

    config: CaseConfig
    workspace: Workspace
    clock: Clock
    run_id: str
    run_dir: Path
    started_at: datetime
    inputs: list[InputRecord] = field(default_factory=list)
    artifacts: list[ArtifactProvenance] = field(default_factory=list)
    commands: list[CommandRecord] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    allow_workspace_in_repo: bool = False
    reference_versions: list[ReferenceVersion] = field(default_factory=list)

    #: Called with every artifact the instant it is registered, before the caller
    #: gets it back. `mva.orchestrator` sets this to the public-export gate.
    #:
    #: It is a hook on the funnel rather than a call at each site because the
    #: failure it replaces was exactly a missed call site: the gate ran on the
    #: submission and on nothing else, while four more artifacts were classified
    #: PUBLIC. A per-site call is a rule a future stage author must remember; a
    #: hook on `register_artifact` is a rule they cannot get wrong, because there
    #: is no other way to produce an artifact.
    on_register: Callable[[RunContext, ArtifactProvenance], None] | None = None

    def artifact_path(self, relative: str) -> Path:
        path = self.run_dir / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def register_artifact(
        self,
        *,
        kind: ArtifactKind,
        path: Path,
        stage: str,
        upstream: Sequence[str] = (),
        row_count: int | None = None,
        tool_versions: Sequence[ToolVersion] = (),
        notes: str = "",
    ) -> ArtifactProvenance:
        """Record how an artifact was produced and how sensitive it is."""
        relative = path.relative_to(self.workspace.root).as_posix()
        provenance = ArtifactProvenance(
            artifact_id=f"ART-{short_hash(relative, 10)}",
            kind=kind,
            relative_path=relative,
            sensitivity=ARTIFACT_SENSITIVITY.get(kind, Sensitivity.SENSITIVE),
            content_hash=hash_file(path),
            size_bytes=path.stat().st_size,
            produced_by_stage=stage,
            upstream_artifact_ids=tuple(upstream),
            tool_versions=tuple(tool_versions),
            created_at=self.clock.now(),
            row_count=row_count,
            notes=notes,
        )
        self.artifacts.append(provenance)
        if self.on_register is not None:
            self.on_register(self, provenance)
        return provenance

    def write_text_artifact(
        self,
        relative: str,
        content: str,
        *,
        kind: ArtifactKind,
        stage: str,
        upstream: Sequence[str] = (),
        row_count: int | None = None,
        notes: str = "",
    ) -> ArtifactProvenance:
        path = self.artifact_path(relative)
        path.write_text(content, encoding="utf-8")
        return self.register_artifact(
            kind=kind,
            path=path,
            stage=stage,
            upstream=upstream,
            row_count=row_count,
            notes=notes,
        )

    def write_json_artifact(
        self,
        relative: str,
        payload: object,
        *,
        kind: ArtifactKind,
        stage: str,
        upstream: Sequence[str] = (),
        row_count: int | None = None,
    ) -> ArtifactProvenance:
        """Write canonical JSON so repeat runs are byte-identical (GP-30)."""
        text = canonical_json(payload)
        return self.write_text_artifact(
            relative,
            text + "\n",
            kind=kind,
            stage=stage,
            upstream=upstream,
            row_count=row_count,
        )


# ---------------------------------------------------------------------------
# Run identity
# ---------------------------------------------------------------------------


def reference_versions_from_manifest(manifest_path: Path) -> tuple[ReferenceVersion, ...]:
    """Record every knowledge table this run depended on, with its content hash.

    Without this, two runs over *different reference data* were indistinguishable:
    the same case config produced the same `config_hash`, the same `run_id`, the
    same input hashes and the same run directory, while emitting a different
    submission. A manifest that cannot tell those apart is not a reproduction
    record, so the reference hashes are recorded here AND folded into the run
    identity by `make_run_id`.
    """
    try:
        raw = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return ()
    if not isinstance(raw, dict):
        return ()
    tables = raw.get("tables")
    if not isinstance(tables, dict):
        return ()

    versions: list[ReferenceVersion] = []
    for name in sorted(tables):
        entry = tables[name]
        if not isinstance(entry, dict):
            continue
        versions.append(
            ReferenceVersion(
                name=str(name),
                version=str(entry.get("version", "unknown")),
                checksum=str(entry.get("sha256", "")) or None,
                retrieved=str(entry.get("retrieved", "")) or None,
            )
        )
    return tuple(versions)


def reference_fingerprint(versions: Sequence[ReferenceVersion]) -> str:
    """A single stable digest over every reference the run consumed."""
    if not versions:
        return "no-references"
    return short_hash([(v.name, v.version, v.checksum) for v in versions], 12)


def make_run_id(config: CaseConfig, *, reference_fingerprint_value: str = "") -> str:
    """Deterministic run identifier.

    Derived from the case id and the config hash rather than from the clock or a
    UUID, so that the same case run twice with the same configuration lands in the
    same directory and produces comparable artifacts (GP-30). A caller wanting
    separate directories supplies a separate workspace.

    Note the deliberate absence of a `clock` parameter: taking one would invite a
    timestamped run id, which is exactly what breaks byte-identical repeat runs.
    """
    if reference_fingerprint_value:
        combined = short_hash([config.config_hash(), reference_fingerprint_value], 12)
        return f"{config.case_id}-{combined}"
    return f"{config.case_id}-{config.config_hash()[:12]}"


def git_state(repo_root: Path) -> tuple[str | None, bool]:
    """Current commit and whether the tree is dirty.

    A dirty tree means the run is not reproducible from the commit alone, which
    the manifest records rather than hides.
    """
    try:
        commit = subprocess.run(  # noqa: S603
            ["git", "-C", str(repo_root), "rev-parse", "HEAD"],  # noqa: S607
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
        )
        status = subprocess.run(  # noqa: S603
            ["git", "-C", str(repo_root), "status", "--porcelain"],  # noqa: S607
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
        )
    except (OSError, subprocess.SubprocessError):
        return None, False
    if commit.returncode != 0:
        return None, False
    return commit.stdout.strip() or None, bool(status.stdout.strip())


def collect_tool_versions() -> tuple[ToolVersion, ...]:
    """Record versions as OBSERVED at runtime, not as declared in pyproject."""
    versions: list[ToolVersion] = [
        ToolVersion(name="python", version=platform.python_version(), source="binary"),
    ]
    for package in ("pydantic", "duckdb", "pyarrow", "cyvcf2", "pysam", "typer", "jinja2"):
        try:
            module = __import__(package)
        except ImportError:
            continue
        version = getattr(module, "__version__", None)
        if isinstance(version, str):
            versions.append(ToolVersion(name=package, version=version))
    return tuple(sorted(versions, key=lambda t: t.name))


def hash_input(path: Path, *, role: str, workspace: Workspace, sensitive: bool) -> InputRecord:
    return InputRecord(
        role=role,
        relative_path=path.relative_to(workspace.root).as_posix(),
        content_hash=hash_file(path),
        size_bytes=path.stat().st_size,
        sensitivity=Sensitivity.SENSITIVE if sensitive else Sensitivity.PUBLIC,
    )


# ---------------------------------------------------------------------------
# Stage 1: case configuration validation
# ---------------------------------------------------------------------------


def validate_case(
    config: CaseConfig,
    workspace: Workspace,
    *,
    clock: Clock | None = None,
    allow_workspace_in_repo: bool = False,
    reference_versions: Sequence[ReferenceVersion] = (),
) -> RunContext:
    """Validate the case and open a run context.

    Runs before any patient file is read: by the time data is in the wrong place
    the mistake is already irreversible.
    """
    active_clock = clock or (demo_clock() if config.synthetic else SystemClock())
    run_id = make_run_id(
        config, reference_fingerprint_value=reference_fingerprint(reference_versions)
    )
    run_dir = workspace.runs_dir / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    context = RunContext(
        config=config,
        workspace=workspace,
        clock=active_clock,
        run_id=run_id,
        run_dir=run_dir,
        started_at=active_clock.now(),
        allow_workspace_in_repo=allow_workspace_in_repo,
        reference_versions=list(reference_versions),
    )

    vcf_path = workspace.path(config.inputs.vcf)
    phenotype_path = workspace.path(config.inputs.phenotype)
    for label, path in (("vcf", vcf_path), ("phenotype", phenotype_path)):
        if not path.is_file():
            msg = (
                f"Case {config.case_id!r} declares input {label!r} at a workspace-relative "
                f"path that does not exist. Place the file under the workspace and retry. "
                f"(The absolute path is deliberately not printed.)"
            )
            raise ConfigError(msg)

    context.inputs.append(
        hash_input(vcf_path, role="vcf", workspace=workspace, sensitive=not config.synthetic)
    )
    context.inputs.append(
        hash_input(
            phenotype_path, role="phenotype", workspace=workspace, sensitive=not config.synthetic
        )
    )

    if allow_workspace_in_repo:
        context.warnings.append(
            "Workspace is inside the repository (--allow-workspace-in-repo). This is "
            "permitted only for synthetic demos and tests; it must never be used for a "
            "real case."
        )
    if not config.synthetic and config.network_profile is not NetworkProfile.OFFLINE_ENFORCED:
        context.warnings.append(
            f"Real case running with network_profile={config.network_profile.value!r}. "
            "The Python guard is a tripwire, not a boundary; pair it with an OS-level "
            "control (see docs/privacy-model.md)."
        )
    return context


def build_run_manifest(context: RunContext, *, repo_root: Path) -> RunManifest:
    """Assemble the complete, self-describing record of the run."""
    commit, dirty = git_state(repo_root)
    return RunManifest(
        run_id=context.run_id,
        case_id=context.config.case_id,
        genome_build=context.config.genome_build.value,
        started_at=context.started_at,
        completed_at=context.clock.now(),
        config_hash=context.config.config_hash(),
        config_snapshot=_json_safe(context.config.model_dump(mode="json")),
        git_commit=commit,
        git_dirty=dirty,
        inputs=tuple(context.inputs),
        artifacts=tuple(context.artifacts),
        commands=tuple(context.commands),
        tool_versions=collect_tool_versions(),
        reference_versions=tuple(context.reference_versions),
        python_version=platform.python_version(),
        platform=f"{platform.system()}-{platform.machine()}",
        network_profile=context.config.network_profile.value,
        synthetic=context.config.synthetic,
        warnings=tuple(context.warnings),
    )


def _json_safe(value: Any) -> dict[str, object]:
    parsed: object = json.loads(canonical_json(value))
    return parsed if isinstance(parsed, dict) else {}


def artifact_digest(context: RunContext) -> dict[str, str]:
    """Path → content hash for every artifact, for the determinism check."""
    ordered = sorted(context.artifacts, key=lambda a: a.relative_path)
    return {a.relative_path: a.content_hash for a in ordered}


def stage_command(
    stage: str, argv: Sequence[str], *, clock: Clock, workspace: Workspace
) -> CommandRecord:
    """Record a command with workspace paths redacted (PRIV-09).

    A provenance manifest may travel with a submission; it must not disclose where
    patient files live or how they are named.
    """
    root = workspace.root.as_posix()
    redacted = tuple(arg.replace(root, "$MVA_WORKSPACE") for arg in argv)
    return CommandRecord(
        stage=stage,
        argv_redacted=redacted,
        exit_code=0,
        duration_seconds=0.0,
        started_at=clock.now(),
    )


def python_environment() -> str:
    return f"{sys.executable} ({platform.python_version()})"


def write_provenance_manifest(context: RunContext, manifest: RunManifest) -> ArtifactProvenance:
    """Write the provenance manifest LAST, so it describes everything before it."""
    return context.write_json_artifact(
        "provenance.json",
        manifest.model_dump(mode="json"),
        kind=ArtifactKind.PROVENANCE_MANIFEST,
        stage="provenance",
    )


def verify_determinism(first: dict[str, str], second: dict[str, str]) -> tuple[bool, list[str]]:
    """Compare two runs' artifact digests (GP-30).

    Two exclusions, both principled:

    * **The provenance manifest** contains the artifact hashes of everything else,
      so comparing it would be circular, and it legitimately differs on the
      git-dirty flag between invocations.
    * **The DuckDB file** is a storage-engine container. Its byte layout depends on
      page allocation and internal bookkeeping the engine does not promise to
      reproduce even for identical logical content. That would be a loophole if it
      were the only record of the evidence store -- so it is not: every table is
      additionally exported to Parquet with pinned compression and row-group size,
      and those exports ARE compared here. The determinism claim is made about the
      data, not about a container's internal layout.
    """
    skip = {"provenance.json", ".duckdb", ".duckdb.wal"}
    differences: list[str] = []
    keys = sorted(set(first) | set(second))
    for key in keys:
        if any(key.endswith(s) for s in skip):
            continue
        a, b = first.get(key), second.get(key)
        if a is None:
            differences.append(f"{key}: missing from first run")
        elif b is None:
            differences.append(f"{key}: missing from second run")
        elif a != b:
            differences.append(f"{key}: hash differs ({a[:12]} vs {b[:12]})")
    return (not differences), differences


def hash_of_text(text: str) -> str:
    return hash_text(text)
