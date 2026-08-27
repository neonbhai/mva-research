"""Command-line interface.

Together with `mva.pipeline` this is the composition root: the only place that
knows the whole stage graph (GP-03).

Two properties of this module are load-bearing rather than stylistic:

* **Output is counts and paths, never records.** An agent runs these commands and
  the output enters its context (PRIV-03, ADR 0008). Nothing here prints a
  genotype, a coordinate from patient data, or a workspace absolute path.
* **The redaction filter is installed before anything else runs**, so a stray
  library log cannot escape ahead of it (PRIV-04).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

from mva.config import (
    CaseConfig,
    NetworkProfile,
    Workspace,
    find_repo_root,
    load_case_config,
    resolve_workspace,
)
from mva.errors import MvaError, PrivacyViolationError

app = typer.Typer(
    name="mva",
    help="Provenance-first rare-disease variant and mechanism analysis.",
    no_args_is_help=True,
    add_completion=False,
)
run_app = typer.Typer(name="run", help="Execute pipeline stages.", no_args_is_help=True)
privacy_app = typer.Typer(name="privacy", help="Privacy boundary tooling.", no_args_is_help=True)
verify_app = typer.Typer(name="verify", help="Verification checks.", no_args_is_help=True)
app.add_typer(run_app)
app.add_typer(privacy_app)
app.add_typer(verify_app)

console = Console(stderr=False)
err_console = Console(stderr=True)


# ---------------------------------------------------------------------------
# Shared options
# ---------------------------------------------------------------------------

ConfigOpt = Annotated[
    Path, typer.Option("--config", "-c", help="Case configuration YAML.", exists=True)
]
DefaultsOpt = Annotated[
    Path | None, typer.Option("--defaults", help="Shared defaults YAML merged under the case.")
]
WorkspaceOpt = Annotated[
    Path | None,
    typer.Option("--workspace", "-w", help="External workspace root. Defaults to $MVA_WORKSPACE."),
]
AllowInRepoOpt = Annotated[
    bool,
    typer.Option(
        "--allow-workspace-in-repo",
        help="Permit a workspace inside the repo. SYNTHETIC DEMOS AND TESTS ONLY.",
    ),
]


def _load(
    config_path: Path, defaults_path: Path | None, workspace: Path | None, allow_in_repo: bool
) -> tuple[CaseConfig, Workspace]:
    """Load config and resolve the workspace, failing closed on privacy checks."""
    config = load_case_config_with_defaults(config_path, defaults_path)
    if allow_in_repo and not config.synthetic:
        msg = (
            f"--allow-workspace-in-repo was passed for case {config.case_id!r}, which is "
            "not marked synthetic. A real case may never write patient data inside the "
            "repository: it is one `git add -A` from being committed, and git history "
            "preserves what `rm` removes."
        )
        raise PrivacyViolationError(msg)
    ws = resolve_workspace(workspace, allow_inside_repo=allow_in_repo)
    return config, ws


def load_case_config_with_defaults(config_path: Path, defaults_path: Path | None) -> CaseConfig:
    """Merge shared defaults beneath the case config.

    Case values win. Merge is shallow-per-section, which is enough for the current
    schema and keeps the precedence obvious to a reader.
    """
    if defaults_path is None:
        return load_case_config(config_path)

    import yaml  # local import: keeps the module import graph flat

    defaults_raw = yaml.safe_load(defaults_path.read_text(encoding="utf-8")) or {}
    case_raw = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    if not isinstance(defaults_raw, dict) or not isinstance(case_raw, dict):
        return load_case_config(config_path)

    merged: dict[str, object] = dict(defaults_raw)
    for key, value in case_raw.items():
        base = merged.get(key)
        if isinstance(base, dict) and isinstance(value, dict):
            combined = dict(base)
            combined.update(value)
            merged[key] = combined
        else:
            merged[key] = value
    return CaseConfig.model_validate(merged)


def _install_privacy_guards(config: CaseConfig, workspace: Workspace) -> str:
    """Install log redaction and, where required, the offline profile."""
    from mva.privacy.redact import install_redaction

    install_redaction()

    if config.network_profile is NetworkProfile.ONLINE:
        return "online"
    from mva.privacy.netguard import configure_reference_cache

    configure_reference_cache(workspace.root)
    return config.network_profile.value


def _fail(message: str) -> None:
    err_console.print(f"[bold red]error[/bold red] {message}")
    raise typer.Exit(code=1)


# ---------------------------------------------------------------------------
# Top-level commands
# ---------------------------------------------------------------------------


@app.command()
def info() -> None:
    """Print the environment the pipeline sees. Counts and versions only."""
    from mva import __version__
    from mva.pipeline import collect_tool_versions, git_state

    repo = find_repo_root()
    commit, dirty = git_state(repo)

    table = Table(title="mva-research environment", show_header=True)
    table.add_column("component")
    table.add_column("value")
    table.add_row("mva version", __version__)
    table.add_row("git commit", (commit[:12] if commit else "unknown"))
    table.add_row("git dirty", "yes" if dirty else "no")
    for tool in collect_tool_versions():
        table.add_row(tool.name, tool.version)
    workspace_set = "set" if _workspace_env_set() else "NOT SET"
    table.add_row("MVA_WORKSPACE", workspace_set)
    console.print(table)
    if not _workspace_env_set():
        console.print(
            "\n[yellow]MVA_WORKSPACE is not set.[/yellow] Point it at a directory outside "
            "this repository and outside ~/Desktop and ~/Documents (iCloud-synced by "
            "default on macOS). See docs/privacy-model.md."
        )


def _workspace_env_set() -> bool:
    import os

    return bool(os.environ.get("MVA_WORKSPACE"))


# ---------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------


@run_app.command("validate")
def run_validate(
    config: ConfigOpt,
    defaults: DefaultsOpt = None,
    workspace: WorkspaceOpt = None,
    allow_workspace_in_repo: AllowInRepoOpt = False,
) -> None:
    """Validate the case configuration and the workspace boundary."""
    from mva.pipeline import validate_case

    try:
        case, ws = _load(config, defaults, workspace, allow_workspace_in_repo)
        context = validate_case(case, ws, allow_workspace_in_repo=allow_workspace_in_repo)
    except MvaError as exc:
        _fail(str(exc))
        return

    console.print(f"[green]ok[/green] case {case.case_id!r} validated")
    console.print(f"  build       {case.genome_build.value}")
    console.print(f"  synthetic   {case.synthetic}")
    console.print(f"  run id      {context.run_id}")
    console.print(f"  inputs      {len(context.inputs)} hashed")
    for warning in context.warnings:
        console.print(f"  [yellow]warning[/yellow] {warning}")


@run_app.command("all")
def run_all(
    config: ConfigOpt,
    defaults: DefaultsOpt = None,
    workspace: WorkspaceOpt = None,
    allow_workspace_in_repo: AllowInRepoOpt = False,
) -> None:
    """Run every stage end to end and write all artifacts."""
    from mva.orchestrator import execute_pipeline

    try:
        case, ws = _load(config, defaults, workspace, allow_workspace_in_repo)
        profile = _install_privacy_guards(case, ws)
        result = execute_pipeline(case, ws, allow_workspace_in_repo=allow_workspace_in_repo)
    except MvaError as exc:
        _fail(str(exc))
        return

    console.print(f"[green]ok[/green] run {result.run_id} complete  (network: {profile})")
    table = Table(show_header=True, title="artifacts")
    table.add_column("kind")
    table.add_column("path (workspace-relative)")
    table.add_column("sensitivity")
    table.add_column("rows", justify="right")
    for artifact in result.manifest.artifacts:
        table.add_row(
            artifact.kind.value,
            artifact.relative_path,
            artifact.sensitivity.value,
            str(artifact.row_count if artifact.row_count is not None else "-"),
        )
    console.print(table)
    console.print(f"\n  top candidate gene   {result.top_gene or 'none'}")
    console.print(f"  candidate pairs      {result.pair_count}")
    console.print(f"  evidence items       {result.evidence_count}")
    console.print(f"  drugs accepted       {result.drugs_accepted}")
    console.print(f"  drugs rejected       {result.drugs_rejected}")
    for warning in result.warnings:
        console.print(f"  [yellow]warning[/yellow] {warning}")


@run_app.command("ingest")
def run_ingest(
    config: ConfigOpt,
    defaults: DefaultsOpt = None,
    workspace: WorkspaceOpt = None,
    allow_workspace_in_repo: AllowInRepoOpt = False,
) -> None:
    """Ingest, normalise and QC the VCF."""
    _run_partial("ingest", config, defaults, workspace, allow_workspace_in_repo)


@run_app.command("annotate")
def run_annotate(
    config: ConfigOpt,
    defaults: DefaultsOpt = None,
    workspace: WorkspaceOpt = None,
    allow_workspace_in_repo: AllowInRepoOpt = False,
) -> None:
    """Annotate variants from local, hash-pinned knowledge tables."""
    _run_partial("annotate", config, defaults, workspace, allow_workspace_in_repo)


@run_app.command("prioritise")
def run_prioritise(
    config: ConfigOpt,
    defaults: DefaultsOpt = None,
    workspace: WorkspaceOpt = None,
    allow_workspace_in_repo: AllowInRepoOpt = False,
) -> None:
    """Filter, pair, score and rank candidates."""
    _run_partial("prioritise", config, defaults, workspace, allow_workspace_in_repo)


@run_app.command("mechanism")
def run_mechanism(
    config: ConfigOpt,
    defaults: DefaultsOpt = None,
    workspace: WorkspaceOpt = None,
    allow_workspace_in_repo: AllowInRepoOpt = False,
) -> None:
    """Construct the mechanism hypothesis for the top candidate."""
    _run_partial("mechanism", config, defaults, workspace, allow_workspace_in_repo)


@run_app.command("drugs")
def run_drugs(
    config: ConfigOpt,
    defaults: DefaultsOpt = None,
    workspace: WorkspaceOpt = None,
    allow_workspace_in_repo: AllowInRepoOpt = False,
) -> None:
    """Generate, direction-check and safety-filter drug hypotheses."""
    _run_partial("drugs", config, defaults, workspace, allow_workspace_in_repo)


@run_app.command("report")
def run_report(
    config: ConfigOpt,
    defaults: DefaultsOpt = None,
    workspace: WorkspaceOpt = None,
    allow_workspace_in_repo: AllowInRepoOpt = False,
) -> None:
    """Render submission, dossier and Track 2 reports."""
    _run_partial("report", config, defaults, workspace, allow_workspace_in_repo)


def _run_partial(
    stage: str,
    config: Path,
    defaults: Path | None,
    workspace: Path | None,
    allow_in_repo: bool,
) -> None:
    """Run the pipeline up to and including `stage`.

    Stages share a run directory and each depends on its predecessors' artifacts,
    so a partial run executes the prefix rather than pretending stages are
    independently invocable.
    """
    from mva.orchestrator import execute_pipeline

    try:
        case, ws = _load(config, defaults, workspace, allow_in_repo)
        _install_privacy_guards(case, ws)
        result = execute_pipeline(case, ws, allow_workspace_in_repo=allow_in_repo, stop_after=stage)
    except MvaError as exc:
        _fail(str(exc))
        return
    console.print(
        f"[green]ok[/green] ran through stage {stage!r}: "
        f"{len(result.manifest.artifacts)} artifacts in {result.run_id}"
    )


# ---------------------------------------------------------------------------
# privacy
# ---------------------------------------------------------------------------


@privacy_app.command("audit")
def privacy_audit(
    repo: Annotated[Path, typer.Option("--repo", help="Repository root.")] = Path(),
    workspace: WorkspaceOpt = None,
    staged: Annotated[bool, typer.Option("--staged", help="Audit staged changes only.")] = False,
    strict: Annotated[bool, typer.Option("--strict", help="Promote warnings to failures.")] = False,
    output: Annotated[
        Path | None, typer.Option("--output", help="Write the markdown report here.")
    ] = None,
) -> None:
    """Scan for privacy boundary violations.

    Reports paths, line numbers and counts. Never the matched content (GP-41) —
    this command's own output is a disclosure vector, because an agent runs it.
    """
    from mva.privacy.audit import run_audit
    from mva.privacy.redact import install_redaction

    # Armed here, by the composition root, and NOT by the audit. The
    # log_redaction_probe check now asserts the state it finds before installing
    # anything, so a command that logs without arming GP-42 is reported as the
    # defect it is. That only works if every entry point arms it — including this
    # one, which logs while it shells out to git.
    install_redaction()

    try:
        report = run_audit(repo.resolve(), workspace=workspace, staged_only=staged, strict=strict)
    except MvaError as exc:
        _fail(str(exc))
        return

    table = Table(show_header=True, title="privacy audit")
    table.add_column("check")
    table.add_column("result")
    table.add_column("findings", justify="right")
    table.add_column("summary")
    for result in report.results:
        status = "[green]pass[/green]" if result.passed else "[bold red]FAIL[/bold red]"
        table.add_row(result.name, status, str(len(result.findings)), result.summary)
    console.print(table)

    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(report.to_markdown(), encoding="utf-8")
        console.print(f"  report written to {output}")

    if not report.passed:
        err_console.print(
            f"\n[bold red]privacy audit FAILED[/bold red]: {', '.join(report.failed_checks)}"
        )
        raise typer.Exit(code=2)
    console.print("\n[green]privacy audit passed[/green]")


@privacy_app.command("classify")
def privacy_classify(
    path: Annotated[Path, typer.Argument(help="Path to classify.")],
) -> None:
    """Show how a path would be classified. Does not read the file."""
    from mva.privacy.classify import classify_path

    console.print(f"{path} -> [bold]{classify_path(path).value}[/bold]")


# ---------------------------------------------------------------------------
# verify
# ---------------------------------------------------------------------------


@verify_app.command("determinism")
def verify_determinism_cmd(
    config: ConfigOpt,
    defaults: DefaultsOpt = None,
    workspace: WorkspaceOpt = None,
    allow_workspace_in_repo: AllowInRepoOpt = False,
) -> None:
    """Run the pipeline twice and compare artifact hashes (GP-30)."""
    from mva.orchestrator import execute_pipeline
    from mva.pipeline import verify_determinism

    try:
        case, ws = _load(config, defaults, workspace, allow_workspace_in_repo)
        _install_privacy_guards(case, ws)
        first = execute_pipeline(case, ws, allow_workspace_in_repo=allow_workspace_in_repo)
        second = execute_pipeline(case, ws, allow_workspace_in_repo=allow_workspace_in_repo)
    except MvaError as exc:
        _fail(str(exc))
        return

    identical, differences = verify_determinism(first.digest, second.digest)
    if identical:
        console.print(
            f"[green]ok[/green] determinism verified: {len(first.digest)} artifacts "
            "byte-identical across two runs"
        )
        return
    err_console.print("[bold red]DETERMINISM FAILURE[/bold red]")
    for difference in differences:
        err_console.print(f"  {difference}")
    raise typer.Exit(code=1)


@verify_app.command("submission")
def verify_submission(
    path: Annotated[Path, typer.Argument(help="Submission CSV to check.", exists=True)],
) -> None:
    """Validate a submission CSV against the verified Track 1 contract."""
    from mva.reporting.track1 import validate_submission

    ok, errors = validate_submission(path.read_text(encoding="utf-8"))
    if ok:
        console.print("[green]ok[/green] submission conforms to the Track 1 contract")
        return
    err_console.print("[bold red]submission INVALID[/bold red]")
    for error in errors:
        err_console.print(f"  {error}")
    raise typer.Exit(code=1)


def main() -> None:
    try:
        app()
    except PrivacyViolationError as exc:
        err_console.print(f"[bold red]privacy violation[/bold red] {exc}")
        sys.exit(2)


if __name__ == "__main__":
    main()
