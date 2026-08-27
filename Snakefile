# =============================================================================
# mva-research — Snakemake workflow
# =============================================================================
#
# WHY SNAKEMAKE (ADR 0001 — docs/decisions/0001-workflow-engine.md)
#   The pipeline is a DAG of file-producing stages and the organizers may rerun a
#   submission, so the DAG has to be explicit, resumable and reproducible.
#
# WHY THESE RULES ARE ONE LINE LONG
#   Every rule shells out to an `mva` subcommand. There is deliberately NO
#   analysis logic in this file: the logic lives in tested, type-checked Python
#   under src/mva/, and the identical DAG runs without Snakemake via
#   `mva run all`. If you find yourself writing a `run:` block here, the code
#   belongs in a stage package instead (ADR 0001, GP-03).
#
# WHERE OUTPUT GOES — ADR 0006, GP-40
#   EVERY output path in this workflow is rooted at {WORKSPACE}/runs/{RUN_ID}/.
#   The workspace is an EXTERNAL absolute path supplied by `MVA_WORKSPACE` or by
#   `--config workspace=/path`; it is NEVER a path inside this repository and is
#   never a committed config value. Patient data inside the repo tree is one
#   `git add -A` away from being committed and stays recoverable from git history
#   afterwards. `tests/unit/test_workflow.py` enforces this mechanically by
#   scanning every `output:` block in this file and in workflow/rules/*.smk.
#
# WHAT RUNS THIS
#   just snakemake [ARGS]    -> uv run snakemake --cores 1 --configfile ...
#   just dag                 -> uv run snakemake --dag --configfile ...
#
# =============================================================================

import os
import shlex
from pathlib import Path

from snakemake.exceptions import WorkflowError
from snakemake.logging import logger


# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
# The synthetic case is the default so that `snakemake --list` and `just dag`
# work with no arguments (ADR 0007: everything demoable runs on the synthetic
# case). A `--configfile` on the command line takes precedence.
configfile: "config/synthetic-case.yaml"


REPO_ROOT = Path(workflow.basedir).resolve()

#: JSON Schema mirroring src/mva/config.py::CaseConfig. Validating here means a
#: typo in a case file is a loud error before any patient byte is opened, rather
#: than a silently ignored setting.
CASE_SCHEMA = "config/schemas/case.schema.json"

try:
    from snakemake.utils import validate as _validate_config

    _validate_config(config, CASE_SCHEMA)
except ImportError:  # pragma: no cover - jsonschema is an optional transitive dep
    # Not fatal: `mva run validate` re-validates with Pydantic and is the
    # authority. This is the early, cheap check, not the only one.
    pass


def _config_path(key: str, fallback: str) -> str:
    """A repo-relative config path, overridable via `--config key=...`."""
    value = config.get(key) or fallback
    return str(value)


def _case_config_path() -> str:
    """The case file actually in force, so `mva` is handed exactly what we read."""
    loaded = list(getattr(workflow, "configfiles", []) or [])
    if loaded:
        return os.path.relpath(str(loaded[-1]), start=str(REPO_ROOT))
    return "config/synthetic-case.yaml"


CASE_CONFIG = _case_config_path()
DEFAULTS_CONFIG = _config_path("defaults_config", "config/default.yaml")

CASE_ID = str(config["case_id"])
GENOME_BUILD = str(config["genome_build"])
SYNTHETIC = bool(config["synthetic"])
NETWORK_PROFILE = str(config.get("network_profile", "offline_best_effort"))


# --------------------------------------------------------------------------
# The workspace boundary (ADR 0006 / GP-40 / PRIV-02 / PRIV-07)
# --------------------------------------------------------------------------
#: Used when no workspace is configured. It is an absolute path that cannot
#: exist, chosen so that `--list`, `--dag` and schema validation still work
#: without a workspace while any attempt to actually RUN fails immediately with a
#: message that names the missing variable. It is deliberately NOT a temp dir and
#: deliberately NOT inside the repo: a silent default is how patient data ends up
#: somewhere nobody chose.
WORKSPACE_UNSET = "/MVA_WORKSPACE-IS-UNSET"

#: Escape hatch used only by `just demo`, which runs the fully synthetic case in
#: a throwaway in-repo directory. Never set this for a real case.
ALLOW_WORKSPACE_IN_REPO = os.environ.get("MVA_ALLOW_WORKSPACE_IN_REPO") == "1"


def _resolve_workspace() -> str:
    """Resolve the external workspace root, refusing a path inside the repo.

    This mirrors `mva.config.resolve_workspace`, which remains the authority and
    additionally rejects cloud-synced roots. Duplicating the containment check
    here buys a failure *before* Snakemake creates a single directory.
    """
    raw = config.get("workspace") or os.environ.get("MVA_WORKSPACE")
    if not raw:
        return WORKSPACE_UNSET

    root = Path(str(raw)).expanduser()
    if not root.is_absolute():
        msg = (
            f"Workspace {str(raw)!r} is relative. The workspace must be an absolute "
            "path outside this repository (ADR 0006). Set MVA_WORKSPACE=/abs/path or "
            "pass --config workspace=/abs/path."
        )
        raise WorkflowError(msg)

    resolved = root.resolve()
    if resolved.is_relative_to(REPO_ROOT) and not ALLOW_WORKSPACE_IN_REPO:
        msg = (
            f"Workspace resolves inside the repository ({resolved.as_posix()}). "
            "Patient data inside the repo tree is one `git add -A` away from being "
            "committed and remains recoverable from git history afterwards (ADR 0006, "
            "GP-40). Move the workspace outside the repo, e.g. "
            "MVA_WORKSPACE=~/mva-workspace. The synthetic demo may opt out with "
            "MVA_ALLOW_WORKSPACE_IN_REPO=1; never do that for a real case."
        )
        raise WorkflowError(msg)
    return resolved.as_posix()


WORKSPACE = _resolve_workspace()

if WORKSPACE == WORKSPACE_UNSET:
    # Loud on stderr so `just dag` (stdout) stays pipeable. The DAG still builds,
    # which is what makes `--list` and `--dag` usable on a fresh clone; only an
    # actual run is refused, by the `onstart` handler at the foot of this file.
    logger.warning(
        "MVA_WORKSPACE is not set. The DAG below is rooted at the placeholder "
        f"{WORKSPACE_UNSET!r} and cannot be executed. Set MVA_WORKSPACE to a "
        "directory OUTSIDE this repository, or pass --config workspace=/abs/path "
        "(ADR 0006)."
    )


def _run_id_from_pipeline() -> str:
    """Ask the pipeline what it will call this run, rather than guessing.

    Snakemake and `mva` must agree on the run directory or the DAG would track
    files that nothing writes. `mva.pipeline.make_run_id` derives it from the case
    id and the hash of the RESOLVED config — case merged over defaults — and
    deliberately takes no clock, because a timestamped run id would break
    byte-identical repeat runs (GP-30). Calling it here means the two cannot
    drift; re-implementing the rule in this file is how they would.

    This is config resolution, not analysis logic, so it does not violate the
    thin-rules constraint of ADR 0001: no rule body gains a line from it.
    """
    try:
        from mva.cli import load_case_config_with_defaults
        from mva.pipeline import make_run_id
    except ImportError:
        # A checkout where `mva` is not installed can still be inspected with
        # --list and --dag; it just cannot be run.
        logger.warning(
            "`mva` is not importable, so the run id cannot be derived from the "
            f"resolved config. Falling back to {CASE_ID!r}, which will NOT match the "
            "directory the pipeline writes to. Install the project (`uv sync`) before "
            "running this workflow for real."
        )
        return CASE_ID

    try:
        case = load_case_config_with_defaults(
            REPO_ROOT / CASE_CONFIG, REPO_ROOT / DEFAULTS_CONFIG
        )
    except Exception as exc:
        msg = (
            f"Could not resolve {CASE_CONFIG} over {DEFAULTS_CONFIG}: {exc}\n"
            "The workflow needs a valid case configuration before it can name the run "
            "directory. Fix the config; `mva run validate --config <cfg> --defaults "
            "<def> --workspace <ws>` reports the same error in full."
        )
        raise WorkflowError(msg) from exc
    return make_run_id(case)


def _resolve_run_id() -> str:
    """The run directory name — the only outside value that reaches a path.

    There are no wildcards anywhere in this workflow, by design. This value is the
    single exception, so it is checked here: nothing supplied through
    `--config run_id=...` or `MVA_RUN_ID` may steer an output out of the
    workspace (ADR 0006).
    """
    override = config.get("run_id") or os.environ.get("MVA_RUN_ID")
    run_id = str(override) if override else _run_id_from_pipeline()
    if run_id != Path(run_id).name or run_id in {"", ".", ".."} or run_id.startswith("-"):
        msg = (
            f"run_id {run_id!r} is not a single safe path segment. It must contain no "
            "'/', no '..' and no leading '-', because it is interpolated directly "
            "into every output path under the workspace (ADR 0006)."
        )
        raise WorkflowError(msg)
    return run_id


RUN_ID = _resolve_run_id()


# --------------------------------------------------------------------------
# Artifact layout — every path below is rooted at WORKSPACE (ADR 0006)
# --------------------------------------------------------------------------
# These paths MIRROR what src/mva/orchestrator.py writes. They are not a naming
# scheme this file is free to choose: Snakemake's DAG is only meaningful if the
# files it tracks are the files the pipeline actually produces. If a stage's
# output moves, it moves here in the same commit.
RUN_DIR = f"{WORKSPACE}/runs/{RUN_ID}"

VARIANTS_DIR = f"{RUN_DIR}/variants"
QC_DIR = f"{RUN_DIR}/qc"
CANDIDATES_DIR = f"{RUN_DIR}/candidates"
SUBMISSION_DIR = f"{RUN_DIR}/submission"
REPORTS_DIR = f"{RUN_DIR}/reports"
EVIDENCE_DIR = f"{RUN_DIR}/evidence"

# Workflow-only, not pipeline artifacts: gate markers and stage logs. They live
# beside the artifacts (inside the workspace, never in the repo) and are excluded
# from the provenance digest, which enumerates registered artifacts rather than
# scanning the directory.
STATUS_DIR = f"{RUN_DIR}/status"
LOG_DIR = f"{RUN_DIR}/logs"

# Inputs are workspace-relative in the case config, so the absolute location of
# patient data appears in no committed file and in no provenance manifest.
INPUT_VCF = f"{WORKSPACE}/{config['inputs']['vcf']}"
INPUT_PHENOTYPE = f"{WORKSPACE}/{config['inputs']['phenotype']}"

# Rules declare the case inputs through these lists rather than directly. With no
# workspace configured the files cannot exist, and Snakemake would refuse to
# build the DAG at all — so `just dag` and `snakemake --list` would fail on a
# fresh clone, where drawing the DAG is exactly what someone wants to do. In that
# state the workflow is documentation, not a plan, so the inputs are omitted.
# This cannot cause a run without its inputs: the `onstart` handler at the foot
# of this file refuses to execute whenever the workspace is unset.
VCF_INPUT = [] if WORKSPACE == WORKSPACE_UNSET else [INPUT_VCF]
PHENOTYPE_INPUT = [] if WORKSPACE == WORKSPACE_UNSET else [INPUT_PHENOTYPE]

# Gate markers. `mva run validate` and `mva privacy audit` are checks, not
# producers; a touched marker is how a check becomes a DAG node.
VALIDATED_OK = f"{STATUS_DIR}/validated.ok"
PROVENANCE_OK = f"{STATUS_DIR}/provenance.ok"
PRIVACY_AUDIT_OK = f"{STATUS_DIR}/privacy-audit.ok"

PROVENANCE_MANIFEST = f"{RUN_DIR}/provenance.json"
EVIDENCE_DB = f"{EVIDENCE_DIR}/evidence.duckdb"

NORMALISED_VARIANTS = f"{VARIANTS_DIR}/normalised.json"
QC_REPORT = f"{QC_DIR}/qc_report.json"
ANNOTATED_VARIANTS = f"{VARIANTS_DIR}/annotated.json"
RANKED_PAIRS = f"{CANDIDATES_DIR}/ranked_pairs.json"
TRACK1_SUBMISSION = f"{SUBMISSION_DIR}/track1_submission.csv"
TRACK1_DOSSIER = f"{REPORTS_DIR}/candidate_dossier.md"
MECHANISM_REPORT = f"{REPORTS_DIR}/mechanism_report.md"
DRUG_HYPOTHESES = f"{REPORTS_DIR}/drug_hypotheses.md"
REJECTION_RECORD = f"{REPORTS_DIR}/rejection_record.md"
TRACK2_REPORT = f"{REPORTS_DIR}/track2_report.md"

TRACK1_ARTIFACTS = [TRACK1_SUBMISSION, TRACK1_DOSSIER]
TRACK2_ARTIFACTS = [TRACK2_REPORT]


# --------------------------------------------------------------------------
# The single seam between this workflow and the pipeline
# --------------------------------------------------------------------------
#: `uv run snakemake` puts the project venv on PATH, so the console script is
#: directly callable. Override with MVA_CMD when driving a container.
MVA = os.environ.get("MVA_CMD", "mva")

#: Extra flags appended to every `mva run` call, e.g.
#: MVA_ALLOW_WORKSPACE_IN_REPO=1 MVA_EXTRA_ARGS=--allow-workspace-in-repo.
MVA_EXTRA_ARGS = os.environ.get("MVA_EXTRA_ARGS", "")


def mva_run(stage: str) -> str:
    """The shell command for one pipeline stage.

    This is the ONLY place a command line is constructed. Keeping it here is what
    keeps the rules one line long (ADR 0001) and guarantees every stage sees the
    same config, defaults and workspace.

    KNOWN COST, stated rather than hidden: `mva run <stage>` executes the whole
    PREFIX up to and including that stage (see `_run_partial` in src/mva/cli.py).
    Stages share a run directory and each depends on its predecessors' artifacts,
    so the CLI does not pretend they are independently invocable — and it is right
    not to. The consequence here is that a full `snakemake` run recomputes the
    early stages once per rule: measured on the synthetic case, roughly 20s
    against 3s for a single `mva run all`.

    Correctness is unaffected — every re-derivation is byte-identical (GP-30), a
    completed run reports "nothing to be done" on the next invocation, and a later
    rule rewriting an earlier rule's output writes the same bytes. So what the DAG
    buys today is explicitness, a rendered `just dag`, and resumability at stage
    granularity after a failure; it does not buy incremental computation.

    Getting that would need a resume-from-artifacts mode in the CLI, not a change
    to these rules. It is deliberately NOT faked here with `ancient()`, which
    would suppress re-runs by declaring inputs timeless and would therefore also
    suppress the legitimate re-run when the VCF or the config actually changes.
    """
    parts = [
        MVA,
        "run",
        stage,
        "--config",
        shlex.quote(CASE_CONFIG),
        "--defaults",
        shlex.quote(DEFAULTS_CONFIG),
        "--workspace",
        shlex.quote(WORKSPACE),
    ]
    command = " ".join(parts)
    if MVA_EXTRA_ARGS:
        command = f"{command} {MVA_EXTRA_ARGS}"
    # Stage output goes to a log INSIDE the workspace, never to the terminal.
    # A pasted traceback is a disclosure path (PRIV-03); Snakemake prints the log
    # path on failure, and the path is all an agent ever needs to see (GP-44).
    return f"{command} > {{log}} 2>&1"


# --------------------------------------------------------------------------
# Targets
# --------------------------------------------------------------------------
# `rule all` is first, so it is the default target.
ALL_ARTIFACTS = [
    VALIDATED_OK,
    NORMALISED_VARIANTS,
    QC_REPORT,
    ANNOTATED_VARIANTS,
    RANKED_PAIRS,
    TRACK1_SUBMISSION,
    TRACK1_DOSSIER,
    MECHANISM_REPORT,
    DRUG_HYPOTHESES,
    REJECTION_RECORD,
    TRACK2_REPORT,
    PROVENANCE_OK,
    PRIVACY_AUDIT_OK,
]


rule all:
    """Everything: both tracks, the provenance join point and the privacy gate."""
    input:
        ALL_ARTIFACTS,


rule track1:
    """Stop after the Track 1 submission and dossier."""
    input:
        RANKED_PAIRS,
        TRACK1_ARTIFACTS,


rule track2:
    """Stop after the Track 2 mechanism-and-intervention report."""
    input:
        DRUG_HYPOTHESES,
        REJECTION_RECORD,
        TRACK2_ARTIFACTS,


# Cheap rules: gates, join points and target aliases. None of these is worth a
# cluster submission even when this workflow is run with a cluster executor.
localrules:
    all,
    track1,
    track2,
    validate,
    provenance,
    privacy_audit,


include: "workflow/rules/ingest.smk"
include: "workflow/rules/annotate.smk"
include: "workflow/rules/prioritise.smk"
include: "workflow/rules/track2.smk"
include: "workflow/rules/report.smk"


onstart:
    if WORKSPACE == WORKSPACE_UNSET:
        raise WorkflowError(
            "No workspace configured. Set MVA_WORKSPACE to a directory OUTSIDE this "
            "repository (and outside ~/Desktop and ~/Documents, which are iCloud-synced "
            "by default on macOS), or pass --config workspace=/abs/path. Patient data "
            "must never be written inside the repo tree (ADR 0006, GP-40)."
        )
