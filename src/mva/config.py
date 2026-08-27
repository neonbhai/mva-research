"""Typed configuration and the workspace privacy boundary.

Two responsibilities:

1. Load and validate a case configuration into typed models (Pydantic), so a typo
   in a YAML key is a loud error rather than a silently ignored setting.
2. Enforce GP-40: the workspace holding patient data must resolve **outside** the
   repository and outside cloud-synced locations. This check runs before any file
   is read, because by the time data has been written to the wrong place the
   mistake is already irreversible.
"""

from __future__ import annotations

import os
from enum import StrEnum
from pathlib import Path
from typing import Any, Self

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from mva.determinism import stable_hash
from mva.errors import ConfigError, WorkspaceError
from mva.models.genome import GenomeBuild

# ---------------------------------------------------------------------------
# Privacy boundary constants (PRIV-07)
# ---------------------------------------------------------------------------

#: Directory names that are cloud-synced by default on macOS or common on any OS.
#: `~/Desktop` and `~/Documents` are iCloud-synced whenever "Desktop & Documents
#: Folders" is enabled, which is the default on a new Mac. A patient VCF placed
#: there is uploaded to Apple within seconds and is outside the researcher's
#: control from that moment.
CLOUD_SYNCED_MARKERS: tuple[str, ...] = (
    "Library/Mobile Documents",
    "Library/CloudStorage",
    "Dropbox",
    "Google Drive",
    "GoogleDrive",
    "OneDrive",
    "Box Sync",
    "pCloud",
    "Sync.com",
    "Nextcloud",
    "iCloud Drive",
)

#: Home-relative directories that are cloud-synced by default on macOS.
CLOUD_SYNCED_HOME_DIRS: tuple[str, ...] = ("Desktop", "Documents")


def path_is_within(child: Path, parent: Path) -> bool:
    """Is ``child`` the same directory as ``parent``, or inside it?

    GP-40's containment check cannot be a string comparison, and
    ``Path.resolve().is_relative_to()`` is a string comparison wearing a costume.
    ``resolve()`` expands symlinks and ``..`` but does **not** case-fold, and the
    default macOS filesystem (APFS) is case-INsensitive: ``/x/REPO/ws`` and
    ``/x/repo/ws`` are one directory with two spellings, and only the second was
    recognised as being inside ``/x/repo``. Spelling a directory in the wrong case
    is not an exotic attack — shell completion, a copied path from Finder and a
    ``$HOME`` written by a different tool all produce it — and the consequence was
    that a workspace physically inside the repository was accepted.

    ``os.path.normcase`` does not help either: on POSIX it is the identity
    function, so it is exactly as blind on macOS as ``resolve()`` is.

    So containment is decided on filesystem IDENTITY. ``(st_dev, st_ino)`` names
    the directory itself rather than a path to it, which is simultaneously robust
    to case, to a symlinked intermediate component, to a hardlinked directory and
    to ``..`` games. Every ancestor of the child is compared against the parent,
    because "inside" means "the parent is one of my ancestors".

    Two fallbacks remain for paths that do not exist yet (a workspace the user has
    not created): the ordinary resolved-prefix test, and a case-folded prefix test.
    The case-folded test can in principle produce a false positive on a
    genuinely case-sensitive filesystem holding two directories that differ only
    in case. That failure direction is the safe one — it refuses a workspace — and
    it is far less likely than the failure it replaces.
    """
    child_real = Path(os.path.realpath(child))
    parent_real = Path(os.path.realpath(parent))

    if child_real == parent_real or child_real.is_relative_to(parent_real):
        return True

    parent_id = _path_identity(parent_real)
    if parent_id is not None:
        for ancestor in (child_real, *child_real.parents):
            if _path_identity(ancestor) == parent_id:
                return True

    folded_child = str(child_real).casefold()
    folded_parent = str(parent_real).casefold()
    return folded_child == folded_parent or folded_child.startswith(folded_parent + os.sep)


def _path_identity(path: Path) -> tuple[int, int] | None:
    """``(st_dev, st_ino)`` — the filesystem's own name for a directory, or None."""
    try:
        info = path.stat()
    except OSError:
        return None
    return (info.st_dev, info.st_ino)


class NetworkProfile(StrEnum):
    """How strictly outbound network is denied during sensitive stages."""

    OFFLINE_ENFORCED = "offline_enforced"
    """Audit hook armed AND an OS-level control asserted."""

    OFFLINE_BEST_EFFORT = "offline_best_effort"
    """Python-level guard only. Honest default: C extensions and subprocesses can
    still reach the network, so this is a tripwire, not a boundary."""

    ONLINE = "online"
    """Permitted only for synthetic cases and public-reference acquisition."""


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True, str_strip_whitespace=True)


# ---------------------------------------------------------------------------
# Scoring weights (GP-32)
# ---------------------------------------------------------------------------


class ScoringWeights(StrictModel):
    """Weights for combining component scores into a composite.

    **Every value here is a heuristic starting point, not a calibrated
    parameter.** They were chosen by reasoning about what matters for a severe
    paediatric recessive disorder, not fitted to any labelled dataset, and they
    carry no claim of clinical validity.

    Changing any weight requires a decision record, a test, and a documented
    before/after comparison (GP-32). Golden expectations are never silently
    re-baselined to accommodate a new weighting.
    """

    analytical_validity: float = Field(default=0.20, ge=0.0, le=1.0)
    """Highest single weight. A variant that is not analytically real cannot be
    causal, however attractive its biology."""

    rarity: float = Field(default=0.18, ge=0.0, le=1.0)
    molecular_consequence: float = Field(default=0.18, ge=0.0, le=1.0)
    inheritance_consistency: float = Field(default=0.16, ge=0.0, le=1.0)
    phenotype_similarity: float = Field(default=0.14, ge=0.0, le=1.0)
    mechanistic_relevance: float = Field(default=0.08, ge=0.0, le=1.0)
    evidence_quality: float = Field(default=0.06, ge=0.0, le=1.0)

    contradiction_penalty_weight: float = Field(
        default=0.25,
        ge=0.0,
        le=1.0,
        description=(
            "Subtracted from the weighted sum. Deliberately larger than the smallest "
            "positive weights: contradicting evidence should be able to move a "
            "candidate meaningfully, not be averaged into irrelevance."
        ),
    )

    @model_validator(mode="after")
    def _positive_weights_sum_to_one(self) -> Self:
        total = (
            self.analytical_validity
            + self.rarity
            + self.molecular_consequence
            + self.inheritance_consistency
            + self.phenotype_similarity
            + self.mechanistic_relevance
            + self.evidence_quality
        )
        if abs(total - 1.0) > 1e-6:
            msg = (
                f"Scoring weights must sum to 1.0 (got {total:.6f}). Adjust the "
                "components explicitly rather than relying on implicit renormalisation, "
                "so that the effect of a change is visible in the diff."
            )
            raise ValueError(msg)
        return self


class PhaseWeights(StrictModel):
    """Multipliers applied to inheritance consistency by phase state (GP-15).

    In-cis is near-disqualifying rather than fully disqualifying: phasing calls
    from short reads can be wrong, and destroying the candidate outright would
    violate GP-13. It is down-ranked hard and flagged, so a human sees it.
    """

    trans_confirmed: float = Field(default=1.00, ge=0.0, le=1.0)
    trans_likely: float = Field(default=0.85, ge=0.0, le=1.0)
    unknown: float = Field(default=0.55, ge=0.0, le=1.0)
    cis_likely: float = Field(default=0.10, ge=0.0, le=1.0)
    cis_confirmed: float = Field(default=0.02, ge=0.0, le=1.0)


class FrequencyThresholds(StrictModel):
    """Allele-frequency cut-points for the rarity component.

    These are ranking cut-points, not filters (GP-13). A variant above
    `max_plausible_recessive` is scored ~0 for rarity but remains in the output,
    flagged, because published pathogenic variants occasionally exceed naive
    frequency expectations in under-represented populations.
    """

    ultra_rare: float = Field(default=1e-5, ge=0.0, le=1.0)
    rare: float = Field(default=1e-4, ge=0.0, le=1.0)
    low_frequency: float = Field(default=1e-3, ge=0.0, le=1.0)
    max_plausible_recessive: float = Field(default=1e-2, ge=0.0, le=1.0)
    min_allele_number: int = Field(
        default=2000,
        ge=0,
        description=(
            "Minimum allele number (sampled chromosomes) a population must report "
            "before it may set the maximum AF. Below it the population is recorded "
            "as excluded rather than used. Guards the popmax against tiny cohorts: "
            "AC=1 in AN=40 is an AF of 0.025 and reads as 'common', which is how a "
            "genuine founder allele gets discarded (ADR 0010, "
            "ASSUMPTION-FREQUENCY-02). A population that reports no allele number "
            "at all stays eligible — an unrecorded cohort size is unknown, not "
            "small (GP-14)."
        ),
    )
    absent_frequency_score: float = Field(
        default=0.5,
        ge=0.0,
        le=1.0,
        description=(
            "Score when NO frequency data exists. Deliberately mid-range, not 1.0: "
            "absence of frequency data is not evidence of rarity (GP-14). It often "
            "means the site is poorly covered in reference cohorts."
        ),
    )


class QualityThresholds(StrictModel):
    """Analytical-validity cut-points. Breaching one flags; it does not delete."""

    min_depth: int = Field(default=10, ge=0)
    min_genotype_quality: int = Field(default=20, ge=0)
    min_allele_balance_het: float = Field(default=0.25, ge=0.0, le=1.0)
    max_allele_balance_het: float = Field(default=0.75, ge=0.0, le=1.0)
    mosaic_allele_balance_floor: float = Field(
        default=0.10,
        ge=0.0,
        le=1.0,
        description=(
            "Below the het band but above this floor, a call is flagged "
            "'possible_mosaic' rather than 'low_allele_balance'. In a mosaic "
            "aneuploidy disorder, skewed allele balance may be signal, not noise."
        ),
    )


# ---------------------------------------------------------------------------
# Case configuration
# ---------------------------------------------------------------------------


class InputPaths(StrictModel):
    """Workspace-relative input locations. Absolute paths are rejected.

    Forcing relative paths keeps the absolute location of patient data out of every
    committed config and every provenance manifest.
    """

    vcf: str
    phenotype: str
    pedigree: str | None = None
    reference_fasta: str | None = None

    @field_validator("vcf", "phenotype", "pedigree", "reference_fasta")
    @classmethod
    def _must_be_relative(cls, value: str | None) -> str | None:
        if value is None:
            return None
        if Path(value).is_absolute():
            msg = (
                f"Input path {value!r} is absolute. Input paths must be relative to the "
                "workspace root so that no committed artifact records where patient "
                "data lives on disk."
            )
            raise ValueError(msg)
        if ".." in Path(value).parts:
            msg = f"Input path {value!r} escapes the workspace via '..'."
            raise ValueError(msg)
        return value


class KnowledgeSources(StrictModel):
    """Which local knowledge adapters to use, and at what version."""

    gene_panel: str = "knowledge/public/gene_panel.tsv"
    hpo_terms: str = "knowledge/public/hpo_terms.tsv"
    gene_phenotype: str = "knowledge/public/gene_phenotype.tsv"
    drug_catalog: str = "knowledge/public/drug_catalog.tsv"
    mechanism_library: str = "knowledge/public/mechanisms.tsv"
    manifest: str = "knowledge/manifests/knowledge.yaml"


class CaseConfig(StrictModel):
    """Everything needed to run one case."""

    case_id: str = Field(min_length=1, pattern=r"^[A-Za-z0-9_\-]+$")
    proband_id: str = Field(
        description=(
            "Identifier used in the submission. The challenge scorer accepts only "
            "'PROBAND01'; any other value hard-fails with 'Unknown proband_id'."
        )
    )
    genome_build: GenomeBuild
    synthetic: bool = Field(
        description=(
            "Explicit, never inferred. Gates safety checks that must not be relaxed by "
            "accident: a real case mislabelled synthetic would skip privacy enforcement."
        )
    )
    description: str = ""

    inputs: InputPaths
    knowledge: KnowledgeSources = Field(default_factory=KnowledgeSources)

    weights: ScoringWeights = Field(default_factory=ScoringWeights)
    phase_weights: PhaseWeights = Field(default_factory=PhaseWeights)
    frequency: FrequencyThresholds = Field(default_factory=FrequencyThresholds)
    quality: QualityThresholds = Field(default_factory=QualityThresholds)

    max_submission_rows: int = Field(
        default=10,
        ge=1,
        le=10,
        description="Challenge hard limit: >10 rows raises ValueError in the scorer.",
    )
    max_pairs_per_gene: int = Field(default=20, ge=1)
    network_profile: NetworkProfile = NetworkProfile.OFFLINE_BEST_EFFORT

    @model_validator(mode="after")
    def _real_cases_must_be_offline(self) -> Self:
        """A real patient case may never run with the network open (PRIV-05)."""
        if not self.synthetic and self.network_profile is NetworkProfile.ONLINE:
            msg = (
                f"Case {self.case_id!r} is not synthetic but requests "
                "network_profile=online. Patient coordinates must never be sent to an "
                "external service. Use offline_enforced, and pre-download any reference "
                "data as a separate public-only acquisition step."
            )
            raise ValueError(msg)
        return self

    def config_hash(self) -> str:
        """Stable hash of the resolved configuration, for the run manifest."""
        return stable_hash(self.model_dump(mode="json"))


# ---------------------------------------------------------------------------
# Workspace resolution (GP-40 / PRIV-02 / PRIV-07)
# ---------------------------------------------------------------------------


class Workspace(StrictModel):
    """A validated external workspace root."""

    root: Path
    repo_root: Path

    def path(self, relative: str) -> Path:
        """Resolve a workspace-relative path, refusing escapes."""
        candidate = (self.root / relative).resolve()
        if not candidate.is_relative_to(self.root.resolve()):
            msg = f"Path {relative!r} escapes the workspace root."
            raise WorkspaceError(msg)
        return candidate

    @property
    def runs_dir(self) -> Path:
        return self.root / "runs"

    @property
    def tmp_dir(self) -> Path:
        return self.root / "tmp"


def find_repo_root(start: Path | None = None) -> Path:
    """Locate the repository root by walking up to the directory holding pyproject.toml."""
    current = (start or Path(__file__)).resolve()
    for parent in [current, *current.parents]:
        if (parent / "pyproject.toml").is_file():
            return parent
    msg = "Could not locate repository root (no pyproject.toml found above this file)."
    raise ConfigError(msg)


def resolve_workspace(
    workspace: str | Path | None = None,
    *,
    repo_root: Path | None = None,
    allow_inside_repo: bool = False,
) -> Workspace:
    """Validate and return the external workspace root.

    Checks, in order:

    * the path is provided (via argument or ``MVA_WORKSPACE``);
    * it exists;
    * after full symlink resolution it is **not** inside the repository;
    * it is not in a cloud-synced location.

    ``allow_inside_repo`` exists solely so the test suite and the synthetic demo can
    use a temporary directory. It is refused for non-synthetic cases by the caller,
    and the privacy audit reports whenever it was used.
    """
    repo = (repo_root or find_repo_root()).resolve()

    raw = workspace if workspace is not None else os.environ.get("MVA_WORKSPACE")
    if raw is None:
        msg = (
            "No workspace configured. Set MVA_WORKSPACE to a directory OUTSIDE this "
            "repository (and outside ~/Desktop and ~/Documents, which are iCloud-synced "
            "by default on macOS), or pass --workspace. Patient data must never be "
            "written inside the repo tree."
        )
        raise WorkspaceError(msg)

    root = Path(raw).expanduser()
    if not root.exists():
        msg = f"Workspace {root.as_posix()!r} does not exist. Create it before running."
        raise WorkspaceError(msg)

    resolved = root.resolve()

    if path_is_within(resolved, repo) and not allow_inside_repo:
        msg = (
            f"Workspace resolves inside the repository ({resolved.as_posix()}). Patient "
            "data inside the repo tree is one `git add -A` away from being committed, "
            "and remains recoverable from git history afterwards. Move the workspace "
            "outside the repo. This check follows symlinks, so a symlinked shortcut "
            "into the repo is caught too."
        )
        raise WorkspaceError(msg)

    _assert_not_cloud_synced(resolved)

    return Workspace(root=resolved, repo_root=repo)


def _assert_not_cloud_synced(resolved: Path) -> None:
    """Refuse workspaces under a known cloud-sync root (PRIV-07)."""
    posix = resolved.as_posix()
    for marker in CLOUD_SYNCED_MARKERS:
        if marker in posix:
            msg = (
                f"Workspace {posix!r} is inside a cloud-synced location ({marker!r}). "
                "Patient data placed there is uploaded to a third party automatically "
                "and cannot be recalled. Choose a local, non-synced directory."
            )
            raise WorkspaceError(msg)

    home = Path.home().resolve()
    for name in CLOUD_SYNCED_HOME_DIRS:
        synced = home / name
        if resolved == synced or resolved.is_relative_to(synced):
            msg = (
                f"Workspace {posix!r} is under ~/{name}, which macOS syncs to iCloud "
                "Drive by default ('Desktop & Documents Folders'). Use a directory "
                "outside the synced set, e.g. ~/mva-workspace, or an encrypted disk "
                "image (see docs/privacy-model.md)."
            )
            raise WorkspaceError(msg)


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def load_case_config(path: Path) -> CaseConfig:
    """Load and validate a case configuration file."""
    if not path.is_file():
        msg = f"Case configuration not found: {path.as_posix()}"
        raise ConfigError(msg)
    try:
        raw: Any = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        msg = f"Case configuration {path.name} is not valid YAML: {exc.__class__.__name__}"
        raise ConfigError(msg) from exc
    if not isinstance(raw, dict):
        msg = f"Case configuration {path.name} must be a YAML mapping at the top level."
        raise ConfigError(msg)
    try:
        return CaseConfig.model_validate(raw)
    except Exception as exc:
        msg = f"Case configuration {path.name} is invalid: {exc}"
        raise ConfigError(msg) from exc
