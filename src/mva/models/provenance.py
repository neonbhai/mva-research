"""Provenance models.

Every artifact this pipeline writes is accompanied by a record of exactly how it
came to exist. This is what makes "deterministic" a checkable claim rather than an
aspiration: two runs whose manifests agree on input hashes, config hash, tool
versions and git commit must produce byte-identical outputs, and the golden tests
assert precisely that.

Provenance also carries the ``sensitivity`` classification, which is the machine
-readable basis for the public-export gate.
"""

from __future__ import annotations

from datetime import datetime
from enum import StrEnum

from pydantic import Field

from mva.models.base import FrozenModel, Sensitivity


class ArtifactKind(StrEnum):
    """What an artifact is, for routing and export decisions."""

    INPUT_VCF = "input_vcf"
    INPUT_PHENOTYPE = "input_phenotype"
    NORMALISED_VARIANTS = "normalised_variants"
    ANNOTATED_VARIANTS = "annotated_variants"
    QC_REPORT = "qc_report"
    CANDIDATE_PAIRS = "candidate_pairs"
    EVIDENCE_DB = "evidence_db"
    SUBMISSION = "submission"
    DOSSIER = "dossier"
    MECHANISM_REPORT = "mechanism_report"
    DRUG_HYPOTHESES = "drug_hypotheses"
    REJECTION_RECORD = "rejection_record"
    TRACK2_REPORT = "track2_report"
    RUN_MANIFEST = "run_manifest"
    PROVENANCE_MANIFEST = "provenance_manifest"
    PRIVACY_AUDIT = "privacy_audit"


class ToolVersion(FrozenModel):
    """A pinned tool or library version, as actually observed at runtime."""

    name: str
    version: str
    source: str = Field(default="python_package", description="python_package | binary | database")
    path: str | None = Field(default=None, description="Resolved path for binaries.")


class ReferenceVersion(FrozenModel):
    """A reference dataset version with an integrity hash."""

    name: str
    version: str
    build: str | None = None
    checksum: str | None = Field(default=None, description="sha256 of the manifest or file.")
    retrieved: str | None = Field(default=None, description="ISO date of retrieval.")
    url: str | None = None


class CommandRecord(FrozenModel):
    """One executed command.

    ``argv`` is stored redacted: workspace paths are replaced with the
    ``$MVA_WORKSPACE`` placeholder so that a provenance manifest — which may be
    published alongside a submission — cannot leak the location or naming of
    patient files.
    """

    stage: str
    argv_redacted: tuple[str, ...]
    exit_code: int
    duration_seconds: float = Field(ge=0.0)
    started_at: datetime


class ArtifactProvenance(FrozenModel):
    """How one artifact was produced, and whether it may leave the workspace."""

    artifact_id: str
    kind: ArtifactKind
    relative_path: str = Field(
        description="Path relative to the workspace root. Never an absolute path."
    )
    sensitivity: Sensitivity
    content_hash: str = Field(description="sha256 of the artifact bytes.")
    size_bytes: int = Field(ge=0)
    produced_by_stage: str
    upstream_artifact_ids: tuple[str, ...] = ()
    tool_versions: tuple[ToolVersion, ...] = ()
    created_at: datetime
    row_count: int | None = Field(default=None, ge=0)
    notes: str = ""

    @property
    def is_exportable(self) -> bool:
        """Only PUBLIC artifacts are candidates for export.

        Note this is necessary, not sufficient: the export gate additionally applies
        an allowlist and re-scans the bytes. Classification is a claim; the scan is
        the verification.
        """
        return self.sensitivity is Sensitivity.PUBLIC


class InputRecord(FrozenModel):
    """A hashed input to the run."""

    role: str = Field(description="e.g. 'vcf', 'phenotype', 'config'")
    relative_path: str
    content_hash: str
    size_bytes: int = Field(ge=0)
    sensitivity: Sensitivity


class RunManifest(FrozenModel):
    """The complete, self-describing record of one pipeline execution."""

    run_id: str
    case_id: str
    genome_build: str
    started_at: datetime
    completed_at: datetime | None = None

    config_hash: str = Field(description="sha256 over the canonicalised resolved config.")
    config_snapshot: dict[str, object] = Field(
        default_factory=dict, description="Fully resolved config, for exact reproduction."
    )
    git_commit: str | None = None
    git_dirty: bool = Field(
        default=False,
        description="True if the working tree had uncommitted changes. A dirty run is "
        "not reproducible from the commit alone and is marked as such.",
    )

    inputs: tuple[InputRecord, ...] = ()
    artifacts: tuple[ArtifactProvenance, ...] = ()
    commands: tuple[CommandRecord, ...] = ()
    tool_versions: tuple[ToolVersion, ...] = ()
    reference_versions: tuple[ReferenceVersion, ...] = ()

    python_version: str = ""
    platform: str = ""
    network_profile: str = Field(
        default="unknown",
        description="'offline_enforced' | 'offline_best_effort' | 'online'. Recorded so a "
        "reader can tell whether any external service could have seen patient coordinates.",
    )
    synthetic: bool = Field(
        default=False,
        description="True for the synthetic demo case. Gates several safety checks.",
    )
    warnings: tuple[str, ...] = ()

    @property
    def sensitive_artifacts(self) -> tuple[ArtifactProvenance, ...]:
        return tuple(a for a in self.artifacts if a.sensitivity is Sensitivity.SENSITIVE)

    @property
    def public_artifacts(self) -> tuple[ArtifactProvenance, ...]:
        return tuple(a for a in self.artifacts if a.sensitivity is Sensitivity.PUBLIC)

    def artifact_by_id(self, artifact_id: str) -> ArtifactProvenance | None:
        for artifact in self.artifacts:
            if artifact.artifact_id == artifact_id:
                return artifact
        return None

    @property
    def is_reproducible(self) -> bool:
        """Whether this run can be exactly reproduced from recorded state alone."""
        return self.git_commit is not None and not self.git_dirty
