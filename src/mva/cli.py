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

import json
import sys
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

import typer
from pydantic import BaseModel, TypeAdapter, ValidationError
from rich.console import Console
from rich.table import Table

from mva.config import (
    CaseConfig,
    InputPaths,
    NetworkProfile,
    Workspace,
    find_repo_root,
    load_case_config,
    path_is_within,
    resolve_workspace,
)
from mva.errors import ConfigError, MvaError, PrivacyViolationError
from mva.models.genome import GenomeBuild

app = typer.Typer(
    name="mva",
    help="Provenance-first rare-disease variant and mechanism analysis.",
    no_args_is_help=True,
    add_completion=False,
)
case_app = typer.Typer(
    name="case",
    help="Set up and run one case. The front door: guided or fully scripted.",
    no_args_is_help=True,
)
run_app = typer.Typer(name="run", help="Execute pipeline stages.", no_args_is_help=True)
privacy_app = typer.Typer(name="privacy", help="Privacy boundary tooling.", no_args_is_help=True)
verify_app = typer.Typer(name="verify", help="Verification checks.", no_args_is_help=True)
discards_app = typer.Typer(
    name="discards",
    help="Ask a run's discard ledger why a variant is not in the output (ADR 0028).",
    no_args_is_help=True,
)
app.add_typer(case_app)
app.add_typer(run_app)
app.add_typer(privacy_app)
app.add_typer(verify_app)
app.add_typer(discards_app)

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
    """Install log redaction and, for any non-ONLINE profile, ARM the offline guard.

    Until this actually armed, the function installed a permanently inert audit
    hook, set two htslib variables, and returned the string ``"offline_enforced"``
    for the CLI to print. `netguard.is_armed()` was `False` for the whole run and a
    DNS lookup succeeded: the run reported a guarantee that no line of code
    provided. `OfflineProfile` existed and was exercised only by its own unit
    tests.

    Arming here is process-lifetime and one-way (see
    :func:`~mva.privacy.netguard.arm_for_process`), because a CLI process on the
    patient-data path has no "after". `mva.orchestrator.execute_pipeline`
    additionally arms a nested, unwinding profile around the stage graph, so a
    caller that reaches the pipeline without going through this CLI — Snakemake, a
    test, a notebook — is covered too. `OfflineProfile` is nesting-safe, so the two
    compose rather than fight.

    Returns the profile label for display only. The label is now a report of state
    the caller can verify with `netguard.is_armed()`, not a claim.
    """
    from mva.privacy.redact import install_redaction

    install_redaction()

    if config.network_profile is NetworkProfile.ONLINE:
        return NetworkProfile.ONLINE.value
    from mva.privacy.netguard import arm_for_process

    arm_for_process(workspace.root)
    return config.network_profile.value


def _report_network_guard(config: CaseConfig) -> None:
    """Print what the guard VERIFIABLY did, and name the half nobody verified.

    The profile label alone reads as a guarantee. Two thirds of it now is one:
    `is_armed()` is observed state, not a claim. The remaining third is not, and
    saying so here is the difference between a control and a reassurance —
    `NetworkProfile.OFFLINE_ENFORCED` is documented as "audit hook armed AND an
    OS-level control asserted", and nothing in this process can confirm the second
    half. A `sandbox-exec` profile, a `pf` rule or a network namespace is asserted
    by the operator outside this program; if it was not, the run had a tripwire and
    the operator believed it had a boundary (TD-06).
    """
    from mva.privacy.netguard import is_armed

    if config.network_profile is NetworkProfile.ONLINE:
        return
    state = "ARMED" if is_armed() else "[bold red]NOT ARMED[/bold red]"
    console.print(f"  network guard        python-level audit hook: {state}")
    if config.network_profile is NetworkProfile.OFFLINE_ENFORCED:
        console.print(
            "  [yellow]note[/yellow] offline_enforced also claims an OS-level control "
            "(sandbox-exec / pf / netns). That half is asserted by the operator and is "
            "NOT verified here. C extensions (pysam/htslib) and subprocesses bypass the "
            "hook entirely — see docs/privacy-model.md and TD-06."
        )


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
    _report_network_guard(case)
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


@verify_app.command("reproduces")
def verify_reproduces(
    submitted: Annotated[
        Path, typer.Argument(help="The submission CSV that was uploaded.", exists=True)
    ],
    produced: Annotated[
        Path, typer.Argument(help="The submission CSV a later run produced.", exists=True)
    ],
) -> None:
    """Compare two Track 1 submissions. Agreement statistics only, never rows.

    `docs/track1-method.md` §4 says the submitted ranking came from a scripted run
    and that the library had never been executed end to end against the real
    callset. This answers the one question that leaves: does the library's own
    output rank the same candidates in the same order?

    Both files are lists of one child's variant coordinates, so nothing printed
    here is drawn from either of them. Every number below is a count, a flag or a
    correlation; files are named by role and basename, and a malformed file is
    reported as a line NUMBER (GP-41, GP-44, PRIV-09, ADR 0008). Exits 1 on
    anything other than an exact reproduction, so it can gate CI.
    """
    from mva.reporting.reproduction import VERDICT_MEANINGS, compare_submission_files

    try:
        report = compare_submission_files(submitted, produced)
    except MvaError as exc:
        _fail(str(exc))
        return

    table = Table(show_header=True, title="track 1 reproduction check")
    table.add_column("statistic")
    table.add_column("value", justify="right")
    table.add_row("submitted rows", str(report.submitted_rows))
    table.add_row("produced rows", str(report.produced_rows))
    table.add_row("rows in both", str(report.intersection))
    table.add_row("submitted only", str(report.submitted_only))
    table.add_row("produced only", str(report.produced_only))
    table.add_row("top-1 identical", "yes" if report.top1_identical else "no")
    for agreement in report.top_k:
        # Both windows are shown because "0 of 3" and "0 of 0" are different
        # facts, and only the second one is not a disagreement (GP-14).
        table.add_row(
            f"top-{agreement.k} agreement",
            f"{agreement.agreeing} of {agreement.produced_window} "
            f"(submitted window {agreement.submitted_window})",
        )
    spearman = "undefined (<2 shared rows)" if report.spearman is None else f"{report.spearman:.4f}"
    table.add_row("spearman rho", spearman)
    console.print(table)

    verdict = report.verdict
    console.print(f"\n  verdict  [bold]{verdict.value}[/bold] — {VERDICT_MEANINGS[verdict]}")
    if report.reproduced:
        console.print("\n[green]ok[/green] the produced ranking reproduces the submitted one")
        return
    err_console.print(
        "\n[bold red]NOT REPRODUCED[/bold red] the two rankings are not the same; "
        "the counts above say how far apart they are. Row content is withheld (GP-41)."
    )
    raise typer.Exit(code=1)


# ---------------------------------------------------------------------------
# discards
# ---------------------------------------------------------------------------


@discards_app.command("explain")
def discards_explain(
    ledger: Annotated[
        Path, typer.Argument(help="Path to a run's discards/discard_ledger.tsv.", exists=True)
    ],
    variant_id: Annotated[
        str, typer.Argument(help="Canonical variant id, e.g. GRCh38:chr15:40200000:C:T.")
    ],
) -> None:
    """Say why one variant is not in a run's output, or that it was never dropped.

    The reviewer-facing half of GP-13 as amended (ADR 0028). "It was filtered" is
    not an answer; "the hard filter removed it because its genotype is hom_ref" is.

    An empty result is deliberately reported as AMBIGUOUS rather than as good news.
    The ledger records removals, not a census: a variant absent from it either
    survived to the output or was never in the input at all, and those are
    different facts that this file cannot tell apart (GP-14).
    """
    from mva.evidence.discards import explain_variant

    rows = explain_variant(ledger, variant_id)
    if not rows:
        console.print(
            f"[yellow]no discard record[/yellow] for {variant_id}\n"
            "  This means one of two things and the ledger cannot distinguish them:\n"
            "  (a) the variant survived every gate and is in variants/annotated.json, or\n"
            "  (b) it was never in the input callset at all.\n"
            "  Check the annotated-variants artifact to tell them apart (GP-14)."
        )
        return
    console.print(f"[bold]{variant_id}[/bold] was dropped:")
    for row in rows:
        value = f" (value: {row.value})" if row.value else ""
        console.print(f"  stage [bold]{row.stage}[/bold] · rule [bold]{row.rule}[/bold]{value}")


@discards_app.command("summary")
def discards_summary(
    ledger: Annotated[
        Path, typer.Argument(help="Path to a run's discards/discard_ledger.tsv.", exists=True)
    ],
) -> None:
    """Counts by stage and rule, recomputed from the file's bytes.

    Recomputed rather than read from `annotation/triage_report.json`, so the two
    disagreeing is detectable: a mismatch means the ledger was truncated, which is
    exactly what happens when the disk fills, and is the one condition under which
    the report's counts describe a file that no longer says what they say.

    Counts only. No coordinate is printed (GP-41, GP-44).
    """
    from mva.evidence.discards import count_discards

    counts = count_discards(ledger)
    if not counts:
        console.print("no discards recorded in this run")
        return
    total = sum(counts.values())
    console.print(f"[bold]{total}[/bold] discarded variant(s)")
    for key, count in counts.items():
        stage, rule = key.split("/", 1)
        console.print(f"  {stage:<12} {rule:<32} {count}")


# ---------------------------------------------------------------------------
# case — the guided front door
# ---------------------------------------------------------------------------
#
# Assembling a valid case used to mean reading four unrelated places first:
# `docs/submission-runbook.md` §1.3 for where the gated data comes from,
# `config/schemas/case.schema.json` for which keys are required,
# `src/mva/phenotype/loader.py` for the phenotype columns, and
# `src/mva/config.py` for the workspace rules. This group is those four answers
# in one command.
#
# Three properties are load-bearing rather than stylistic:
#
# * **One code path, two modes.** `resolve_case` does every check. The only
#   difference between a guided run and a scripted one is the `Prompter` it is
#   handed; `None` means "there is nobody to ask", so a missing value is an
#   immediate error with a non-zero exit and never a hang.
# * **Nothing here opens an input file.** Paths and extensions are checked by
#   name. The contents belong to the pipeline (GP-44, ADR 0008).
# * **No value is validated twice.** The workspace rules come from
#   `resolve_workspace`, the build from `GenomeBuild.parse`, the identifiers from
#   the constraints declared on `CaseConfig`, and the relative-path rule from
#   `InputPaths`. A second copy of a rule here would be a second rule.

#: Name suffixes a variant-call input may plausibly carry.
#:
#: Checked on the NAME only. A path and an extension are all that can be
#: verified without reading patient data, and reading it is not this command's
#: job: it exists to get a valid case assembled, not to inspect a genome.
VCF_SUFFIXES: tuple[str, ...] = (
    ".vcf",
    ".vcf.gz",
    ".vcf.bgz",
    ".gvcf",
    ".gvcf.gz",
    ".g.vcf",
    ".g.vcf.gz",
    ".bcf",
)

#: Name suffixes the phenotype loader can actually read. `read_tsv_rows` calls
#: `Path.read_text`, so a gzipped TSV is refused here rather than three stages
#: later, after the workspace and the VCF have already been accepted.
PHENOTYPE_SUFFIXES: tuple[str, ...] = (".tsv", ".txt", ".tab")

#: Raw reads, refused by name WITH the reason rather than quietly accepted.
#: This pipeline does no alignment and no SV/CNV calling (TD-03), so a BAM, CRAM
#: or FASTQ is not an input it consumes — and a command that accepted one would
#: be implying otherwise.
RAW_READ_SUFFIXES: tuple[str, ...] = (
    ".bam",
    ".cram",
    ".sam",
    ".ubam",
    ".fastq",
    ".fastq.gz",
    ".fq",
    ".fq.gz",
    ".sra",
)

#: How many invalid answers one question takes before the command gives up. An
#: unbounded re-prompt loop is a hang when stdin is a file of the wrong thing.
MAX_PROMPT_ATTEMPTS: int = 5

#: Default workspace-relative destination for the phenotype template.
DEFAULT_PHENOTYPE_TEMPLATE: str = "inputs/phenotype.tsv"

#: The shape of one phenotype row, shipped COMMENTED OUT in the template. The
#: term is a real public ontology identifier and the row is fabricated; it is
#: here to be copied, not to be run. See `phenotype_template_text`.
EXAMPLE_PHENOTYPE_ROW: tuple[str, ...] = (
    "HP:0000252",
    "Microcephaly",
    "observed",
    "congenital",
    "clinical_summary",
    "1.0",
    "",
)

#: The one `proband_id` the challenge scorer accepts. Not enforced — the schema
#: permits any non-empty string and this command does not narrow it — but
#: reported, because the failure mode otherwise is a scored submission of zero.
EXPECTED_PROBAND_ID: str = "PROBAND01"

#: Asks one question and returns the raw answer. Injected so that the guided and
#: the scripted paths are one path with different sources of answers.
Prompter = Callable[[str, str | None], str]


@dataclass(frozen=True)
class CaseRequest:
    """Raw, unvalidated values as supplied on the command line.

    `None` means "not supplied", which is the only thing that triggers a prompt.
    A value that was supplied is validated exactly as a typed answer is.
    """

    workspace: str | None
    vcf: str | None
    phenotype: str | None
    genome_build: str | None
    case_id: str | None
    proband_id: str | None
    synthetic: bool
    description: str


@dataclass(frozen=True)
class ResolvedCase:
    """A case that passed every existing validator, ready to render."""

    workspace: Workspace
    case_id: str
    proband_id: str
    genome_build: GenomeBuild
    synthetic: bool
    description: str
    vcf: str
    """Workspace-relative, as `InputPaths` requires (ADR 0006)."""
    phenotype: str


@dataclass(frozen=True)
class CheckRow:
    """One line of `mva case check`. `state` is 'ok', 'FAIL' or 'note'."""

    name: str
    state: str
    detail: str


# --------------------------------------------------------------- field parsers


def _declared_field_parser(model: type[BaseModel], field: str) -> Callable[[str], str]:
    """A parser for one field, built from the model's OWN declared constraints.

    `CaseConfig.model_fields["case_id"]` carries the pattern and the length bound
    written in `mva.config`. Validating an answer through it means this command
    cannot drift from the schema it is writing for: there is no second copy of
    the pattern here to forget to update.
    """
    adapter: TypeAdapter[str] = TypeAdapter(Annotated[str, model.model_fields[field]])
    return adapter.validate_python


def _parse_workspace(raw: str, repo_root: Path) -> Workspace:
    """Validate a workspace with the existing boundary check.

    Every rule — exists, outside the repo after symlink resolution, not under a
    cloud-synced root — belongs to `resolve_workspace` (GP-40, PRIV-07). There is
    deliberately no bypass flag: if that function refuses a directory, this
    command reports why and asks again.
    """
    if not raw.strip():
        msg = "No workspace given. It is a directory outside this repo, e.g. ~/mva-workspace."
        raise ConfigError(msg)
    return resolve_workspace(raw.strip(), repo_root=repo_root)


def _parse_proband_id(raw: str) -> str:
    """Validate `proband_id` against `CaseConfig`, plus the emptiness the schema adds.

    `config/schemas/case.schema.json` gives this field `minLength: 1`; the
    pydantic model does not, and pydantic is the authority (so the empty string
    would reach Snakemake's schema check and fail there instead). Refusing it
    here is the same rule reported earlier, not a new one.
    """
    value = _declared_field_parser(CaseConfig, "proband_id")(raw).strip()
    if not value:
        msg = (
            "proband_id is empty. case.schema.json requires at least one character, "
            f"and the challenge scorer accepts only {EXPECTED_PROBAND_ID!r}."
        )
        raise ConfigError(msg)
    return value


def _matches_suffix(name: str, suffixes: tuple[str, ...]) -> bool:
    lowered = name.lower()
    return any(lowered.endswith(suffix) for suffix in suffixes)


def _parse_workspace_input(
    raw: str,
    *,
    workspace: Workspace,
    label: str,
    suffixes: tuple[str, ...],
) -> str:
    """Resolve one input path to the workspace-relative form a case config stores.

    Accepts either an absolute path or a path already relative to the workspace,
    because both are what a person has to hand. Returns the relative form, which
    is the only form `InputPaths` permits: a committed config must not record
    where patient data lives (ADR 0006).

    **No path is echoed back.** `mva.pipeline.validate_case` sets the precedent
    for the same reason — this output goes to a terminal, a log and an agent's
    context (GP-44, ADR 0008), and the workspace layout is itself sensitive. The
    label and the reason are enough to act on.
    """
    token = raw.strip()
    if not token:
        msg = f"No path given for {label}."
        raise ConfigError(msg)

    if _matches_suffix(token, RAW_READ_SUFFIXES):
        msg = (
            f"The path given for {label} names a raw-read file. This pipeline consumes "
            "called variants only: it does no alignment and no SV/CNV calling (TD-03), "
            f"so a BAM, CRAM, SAM or FASTQ is not an input it can use. Supply the VCF "
            "produced from those reads."
        )
        raise ConfigError(msg)
    if not _matches_suffix(token, suffixes):
        msg = (
            f"The path given for {label} does not end in one of: {', '.join(suffixes)}. "
            "The extension is checked by name only — nothing here opens the file — so a "
            "file with the wrong name must be renamed or a different one supplied."
        )
        raise ConfigError(msg)

    candidate = Path(token).expanduser()
    if not candidate.is_absolute():
        candidate = workspace.root / candidate
    resolved = candidate.resolve()

    if not path_is_within(resolved, workspace.root):
        msg = (
            f"The path given for {label} is outside the workspace. A case config records "
            "workspace-relative paths only (ADR 0006), so the file has to live under the "
            "workspace root. Move or copy it there and retry — this command will not copy "
            "patient data on your behalf."
        )
        raise ConfigError(msg)

    try:
        relative = resolved.relative_to(workspace.root)
    except ValueError as exc:
        msg = (
            f"The path given for {label} is inside the workspace but cannot be expressed "
            "relative to it, which happens when the two are spelled in different letter "
            "cases on a case-insensitive filesystem. Spell both the same way and retry."
        )
        raise ConfigError(msg) from exc

    # Re-resolve through the workspace itself, which refuses an escape, and then
    # ask whether the file is there. `is_file` is a stat, not a read.
    if not workspace.path(relative.as_posix()).is_file():
        msg = (
            f"No file at the path given for {label}. The path is deliberately not echoed "
            "(GP-44); check it with `ls` in your own shell. If the file is not downloaded "
            "yet, see docs/submission-runbook.md §1.3 for the gated dataset."
        )
        raise ConfigError(msg)
    return relative.as_posix()


# --------------------------------------------------------------- asking


def _console_prompter(question: str, default: str | None) -> str:
    """Read one answer from stdin. Works on a terminal and on a pipe alike."""
    if default is None:
        return str(typer.prompt(question))
    return str(typer.prompt(question, default=default))


def _select_prompter(*, no_input: bool, interactive: bool) -> Prompter | None:
    """Decide whether anybody can be asked a question.

    `--no-input` wins outright: a caller that passes it has guaranteed itself a
    non-zero exit instead of a prompt, and that guarantee is the point.
    `--interactive` forces prompting when stdin is a pipe, which is how this
    command is driven by a test or an expect-style script. Otherwise prompting
    follows the TTY, so a script, a CI job or an agent never gets asked.
    """
    if no_input:
        return None
    if interactive:
        return _console_prompter
    return _console_prompter if sys.stdin.isatty() else None


def _explain(exc: Exception) -> str:
    """One readable line for a rejected answer, carrying no echoed input value.

    Pydantic appends `input_value=...` to a `ValidationError` by default, and for
    a path under the workspace that is exactly the thing not to print, so the
    message is rebuilt from the field location and the constraint (PRIV-09).
    """
    if isinstance(exc, ValidationError):
        parts: list[str] = []
        for error in exc.errors():
            location = ".".join(str(piece) for piece in error["loc"]) or "value"
            parts.append(f"{location}: {error['msg']}")
        return "; ".join(parts)
    return _without_doc_links(" ".join(str(exc).split()))


#: Pydantic appends this plus a URL to every error in a multi-error message. In a
#: three-column table it crowds out the part that says which field is wrong.
_DOC_LINK_MARKER: str = "For further information visit"


def _without_doc_links(message: str) -> str:
    """Drop pydantic's per-error documentation URLs, keeping the field and reason."""
    while _DOC_LINK_MARKER in message:
        head, _, tail = message.partition(_DOC_LINK_MARKER)
        _, _, rest = tail.strip().partition(" ")
        message = f"{head.strip()} {rest.strip()}".strip()
    return message


def _read_answer(prompter: Prompter, question: str, default: str | None) -> str:
    try:
        return prompter(question, default)
    except (EOFError, KeyboardInterrupt, typer.Abort) as exc:
        msg = (
            "Input ended before every required value was supplied, so nothing was "
            "written. Pass the remaining values as flags, or use --no-input to get a "
            "clear error instead of a question."
        )
        raise ConfigError(msg) from exc


def _ask[T](
    *,
    label: str,
    flag: str,
    supplied: str | None,
    question: str,
    default: str | None,
    parse: Callable[[str], T],
    prompter: Prompter | None,
) -> T:
    """Obtain one validated value, from a flag or from whoever is answering.

    The validation is identical either way, which is the whole reason this
    function exists: a flag value is not trusted more than a typed one, and a
    typed one gets the same validator rather than a friendlier approximation of
    it. With no prompter, a missing or invalid value raises.
    """
    raw = supplied
    for _ in range(MAX_PROMPT_ATTEMPTS):
        if raw is None:
            if prompter is None:
                raise ConfigError(_missing_value_message(label, flag))
            raw = _read_answer(prompter, question, default)
        try:
            return parse(raw)
        except (MvaError, ValueError) as exc:
            reason = _explain(exc)
            if prompter is None:
                msg = f"--{flag} ({label}): {reason}"
                raise ConfigError(msg) from exc
            err_console.print(f"  [bold red]invalid[/bold red] {label}: {reason}")
            raw = None
    msg = (
        f"Five invalid answers for {label!r}; giving up without writing anything. "
        f"The accepted form is described by --{flag} in `mva case init --help`."
    )
    raise ConfigError(msg)


def _missing_value_message(label: str, flag: str) -> str:
    return (
        f"No value for {label}. Supply --{flag}, or run from a terminal without "
        "--no-input to be asked for it."
    )


def _missing_flags(request: CaseRequest) -> list[str]:
    """Every required flag the caller left out, reported together.

    Together rather than one at a time because the caller is usually a script or
    an agent, which cannot iterate through six round trips as cheaply as a person
    can answer six questions.
    """
    supplied: tuple[tuple[str, str | None], ...] = (
        ("workspace", request.workspace),
        ("vcf", request.vcf),
        ("phenotype", request.phenotype),
        ("genome-build", request.genome_build),
        ("case-id", request.case_id),
        ("proband-id", request.proband_id),
    )
    return [f"--{flag}" for flag, value in supplied if value is None]


def resolve_case(
    request: CaseRequest,
    *,
    prompter: Prompter | None,
    repo_root: Path,
) -> ResolvedCase:
    """Validate a case request into a `ResolvedCase`, asking for what is missing.

    The workspace is resolved first because the two input paths are meaningless
    without it: they are stored relative to it and checked against it.

    Two prompts carry a default (`GRCh38`, and `PROBAND01` for the proband id).
    A default shown to a person who can see and overtype it is not the same thing
    as a silent flag default, and neither is available in `--no-input` mode —
    GP-11 makes a guessed genome build a defect, not a convenience.
    """
    if prompter is None:
        missing = _missing_flags(request)
        if missing:
            msg = (
                f"Missing required value(s): {', '.join(missing)}. Supply them as flags, "
                "or drop --no-input and run from a terminal to be asked. Nothing was "
                "written."
            )
            raise ConfigError(msg)

    workspace = _ask(
        label="workspace",
        flag="workspace",
        supplied=request.workspace,
        question="Workspace root (outside this repo, e.g. ~/mva-workspace)",
        default=None,
        parse=lambda raw: _parse_workspace(raw, repo_root),
        prompter=prompter,
    )
    vcf = _ask(
        label="vcf",
        flag="vcf",
        supplied=request.vcf,
        question="VCF path (absolute, or relative to the workspace)",
        default=None,
        parse=lambda raw: _parse_workspace_input(
            raw, workspace=workspace, label="vcf", suffixes=VCF_SUFFIXES
        ),
        prompter=prompter,
    )
    phenotype = _ask(
        label="phenotype",
        flag="phenotype",
        supplied=request.phenotype,
        question="Phenotype TSV path (see `mva case phenotype-template`)",
        default=None,
        parse=lambda raw: _parse_workspace_input(
            raw, workspace=workspace, label="phenotype", suffixes=PHENOTYPE_SUFFIXES
        ),
        prompter=prompter,
    )
    genome_build = _ask(
        label="genome build",
        flag="genome-build",
        supplied=request.genome_build,
        question="Genome build (GRCh38 or GRCh37)",
        default=GenomeBuild.GRCH38.value,
        parse=GenomeBuild.parse,
        prompter=prompter,
    )
    case_id = _ask(
        label="case id",
        flag="case-id",
        supplied=request.case_id,
        question="Case id (letters, digits, '-' and '_'; used as a directory name)",
        default=None,
        parse=_declared_field_parser(CaseConfig, "case_id"),
        prompter=prompter,
    )
    proband_id = _ask(
        label="proband id",
        flag="proband-id",
        supplied=request.proband_id,
        question="Proband id",
        default=EXPECTED_PROBAND_ID,
        parse=_parse_proband_id,
        prompter=prompter,
    )
    return ResolvedCase(
        workspace=workspace,
        case_id=case_id,
        proband_id=proband_id,
        genome_build=genome_build,
        synthetic=request.synthetic,
        description=request.description,
        vcf=vcf,
        phenotype=phenotype,
    )


# --------------------------------------------------------------- rendering


def case_config_model(resolved: ResolvedCase) -> CaseConfig:
    """Build the typed config, so the file is written only if the model accepts it.

    Pydantic is the authority on a case configuration (`config/schemas/case.schema.json`
    says so itself), and this is the call that consults it. `InputPaths` is what
    refuses an absolute or `..`-bearing input path.
    """
    return CaseConfig(
        case_id=resolved.case_id,
        proband_id=resolved.proband_id,
        genome_build=resolved.genome_build,
        synthetic=resolved.synthetic,
        description=resolved.description,
        inputs=InputPaths(vcf=resolved.vcf, phenotype=resolved.phenotype),
    )


def render_case_config(resolved: ResolvedCase) -> str:
    """Render the case YAML.

    The single writer, which is what makes the guided and the scripted mode
    byte-identical for identical answers: neither of them formats anything
    itself. Hand-rendered rather than dumped so the file carries the two
    sentences a reader of it needs, and so key order is a property of this
    function rather than of a serialiser's defaults (GP-30).
    """
    lines = [
        "# Case configuration written by `mva case init`.",
        "#",
        "# Input paths are WORKSPACE-RELATIVE, and the workspace root is deliberately",
        "# absent: it is per-machine, and a committed config must not record where",
        "# patient data lives (ADR 0006, GP-40). Supply it at run time with --workspace",
        "# or MVA_WORKSPACE.",
        f"case_id: {_yaml_scalar(resolved.case_id)}",
        f"proband_id: {_yaml_scalar(resolved.proband_id)}",
        f"genome_build: {resolved.genome_build.value}",
        f"synthetic: {'true' if resolved.synthetic else 'false'}",
    ]
    if resolved.description:
        lines.append(f"description: {_yaml_scalar(resolved.description)}")
    lines.extend(
        [
            "",
            "inputs:",
            f"  vcf: {_yaml_scalar(resolved.vcf)}",
            f"  phenotype: {_yaml_scalar(resolved.phenotype)}",
        ]
    )
    return "\n".join(lines) + "\n"


def _yaml_scalar(value: str) -> str:
    """Quote a string as a YAML double-quoted scalar.

    `json.dumps` is exactly that syntax for this subset, and it is total: a value
    holding a quote, a colon, a `#` or a leading digit survives it, which a bare
    scalar would not.
    """
    return json.dumps(value)


def phenotype_template_text() -> str:
    """The phenotype TSV header, with the vocabulary taken from the loader itself.

    The columns and the accepted `status` spellings are read from
    `mva.phenotype.loader`, so this template cannot describe a schema the parser
    does not implement. Guessing that vocabulary is the thing this file exists to
    stop: the loader refuses an unrecognised status rather than defaulting it
    (GP-14), and "status" is the one column where a wrong guess silently inverts
    the evidence.
    """
    from mva.models.phenotype import Onset
    from mva.phenotype.loader import (
        OPTIONAL_COLUMNS,
        REQUIRED_COLUMNS,
        allowed_status_summary,
    )

    columns = (*REQUIRED_COLUMNS, *OPTIONAL_COLUMNS)
    lines = [
        "# Phenotype profile. One row per HPO term, tab-separated.",
        "# '#' comment lines and blank lines are ignored; the first other line is the header.",
        "#",
        f"# Required: {', '.join(REQUIRED_COLUMNS)}",
        f"# Optional: {', '.join(OPTIONAL_COLUMNS)}",
        "#",
        "# hpo_id   'HP:' then exactly seven digits ('HP_0000252' is accepted too). A",
        "#          truncated id is an error, never zero-padded.",
        "# status   An unlisted spelling is refused, not defaulted: guessing 'observed'",
        "#          invents a finding and guessing 'excluded' turns a missing assessment",
        "#          into evidence against a gene (GP-14). 'not_assessed' is NOT",
        "#          'excluded'. Accepted spellings, by canonical value:",
    ]
    lines.extend(f"#            {group}" for group in allowed_status_summary().split("; "))
    lines.extend(
        [
            f"# onset    One of: {', '.join(onset.value for onset in Onset)}. Empty means unknown.",
            "# extraction_confidence  0.0-1.0. Empty means 1.0, which claims only that the",
            "#          term was transcribed correctly from a table, not clinical certainty.",
            "# notes    Free text, and NOT the place for a clinical narrative: anything",
            "#          derived from one belongs in a source_excerpt_hash (PRIV-09).",
            "#",
            "# The example row below is commented out deliberately. An uncommented one is",
            "# a finding no clinician recorded, and it becomes data the moment somebody",
            "# forgets to delete it — which is the mistake this whole file exists to",
            "# prevent (GP-14). Copy its shape; do not ship it.",
            "\t".join(columns),
            "#" + "\t".join(EXAMPLE_PHENOTYPE_ROW),
        ]
    )
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------- reporting


def _display_path(path: Path, *, repo_root: Path, workspace: Workspace | None) -> str:
    """A copy-pasteable spelling of a path that never reveals the workspace root.

    A path under the workspace is printed through `$MVA_WORKSPACE` rather than
    expanded, because the workspace layout is sensitive and this output reaches a
    log and an agent's context (GP-44, ADR 0008).
    """
    resolved = path.expanduser().resolve()
    if workspace is not None and path_is_within(resolved, workspace.root):
        relative = resolved.relative_to(workspace.root).as_posix()
        return f'"$MVA_WORKSPACE/{relative}"'
    if path_is_within(resolved, repo_root):
        return resolved.relative_to(repo_root).as_posix()
    return resolved.as_posix()


def _next_command(config_arg: str, defaults_arg: str | None) -> str:
    parts = ["mva run all", f"--config {config_arg}"]
    if defaults_arg is not None:
        parts.append(f"--defaults {defaults_arg}")
    parts.append('--workspace "$MVA_WORKSPACE"')
    return " ".join(parts)


def _print_case_summary(resolved: ResolvedCase, out_path: Path, repo_root: Path) -> None:
    """Report what was written. Workspace-relative input paths, nothing absolute."""
    console.print(
        f"[green]ok[/green] wrote "
        f"{_display_path(out_path, repo_root=repo_root, workspace=resolved.workspace)}"
    )
    console.print(f"  case id        {resolved.case_id}")
    console.print(f"  proband id     {resolved.proband_id}")
    console.print(f"  genome build   {resolved.genome_build.value}")
    console.print(f"  synthetic      {'true' if resolved.synthetic else 'false'}")
    console.print(f"  vcf            {resolved.vcf}   (workspace-relative)")
    console.print(f"  phenotype      {resolved.phenotype}   (workspace-relative)")
    console.print("  workspace      validated; its path is not printed (GP-44)")
    if resolved.proband_id != EXPECTED_PROBAND_ID:
        console.print(
            f"  [yellow]note[/yellow] the challenge scorer accepts only "
            f"proband_id={EXPECTED_PROBAND_ID!r} and hard-fails anything else with "
            "'Unknown proband_id' (docs/references/track1-submission-contract.md)."
        )
    if not resolved.synthetic:
        console.print(
            "  [yellow]note[/yellow] network_profile is left at its default, so the run "
            "will warn that the Python guard is a tripwire rather than a boundary. Set "
            "`network_profile: offline_enforced` only once an OS-level control is "
            "actually asserted (docs/privacy-model.md, TD-06)."
        )


# --------------------------------------------------------------- commands


@case_app.command("init")
def case_init(
    *,
    workspace: Annotated[
        str | None,
        typer.Option("--workspace", "-w", help="External workspace root. Outside this repo."),
    ] = None,
    vcf: Annotated[
        str | None,
        typer.Option("--vcf", help="VCF path: absolute, or relative to the workspace."),
    ] = None,
    phenotype: Annotated[
        str | None,
        typer.Option("--phenotype", help="Phenotype TSV path, same two forms."),
    ] = None,
    genome_build: Annotated[
        str | None,
        typer.Option("--genome-build", help="GRCh38 or GRCh37. Never guessed (GP-11)."),
    ] = None,
    case_id: Annotated[
        str | None,
        typer.Option("--case-id", help="Case id; used as a directory name under runs/."),
    ] = None,
    proband_id: Annotated[
        str | None,
        typer.Option("--proband-id", help=f"Proband id. The scorer accepts {EXPECTED_PROBAND_ID}."),
    ] = None,
    out: Annotated[
        Path | None,
        typer.Option("--out", "-o", help="Where to write. Default: config/<case-id>.yaml."),
    ] = None,
    synthetic: Annotated[
        bool,
        typer.Option("--synthetic", help="Mark the case synthetic. Default: NOT synthetic."),
    ] = False,
    description: Annotated[str, typer.Option("--description", help="Free-text note.")] = "",
    defaults: Annotated[
        Path | None,
        typer.Option("--defaults", help="Shared defaults YAML to validate and run against."),
    ] = None,
    no_input: Annotated[
        bool,
        typer.Option("--no-input", help="Never prompt. A missing value is an error, not a wait."),
    ] = False,
    interactive: Annotated[
        bool,
        typer.Option("--interactive", help="Prompt even when stdin is a pipe."),
    ] = False,
    yes: Annotated[
        bool, typer.Option("--yes", "-y", help="Assume yes; skip confirmation.")
    ] = False,
    force: Annotated[
        bool, typer.Option("--force", help="Overwrite an existing --out file.")
    ] = False,
    run: Annotated[
        bool, typer.Option("--run", help="Run the pipeline once the case is written.")
    ] = False,
) -> None:
    """Assemble a valid case configuration, and optionally run it.

    Guided when stdin is a terminal and something is missing; silent when every
    value is supplied. The same validators run either way — the workspace
    boundary from `mva.config.resolve_workspace`, the build from
    `GenomeBuild.parse`, the identifiers from the constraints on `CaseConfig`,
    the relative-path rule from `InputPaths` — so the two modes cannot disagree
    about what is valid, and they write identical bytes for identical answers.

    No input file is opened. Paths and extensions are checked by name, and no
    path from the workspace is printed back (GP-44, ADR 0008).

    Scripted:

        mva case init --no-input --workspace ~/mva-workspace \\
            --vcf inputs/proband.vcf --phenotype inputs/phenotype.tsv \\
            --genome-build GRCh38 --case-id real-case-01 --proband-id PROBAND01

    Guided: run it with no flags at all.
    """
    repo_root = find_repo_root()
    request = CaseRequest(
        workspace=workspace,
        vcf=vcf,
        phenotype=phenotype,
        genome_build=genome_build,
        case_id=case_id,
        proband_id=proband_id,
        synthetic=synthetic,
        description=description,
    )
    prompter = _select_prompter(no_input=no_input, interactive=interactive)
    confirmed = yes or no_input

    try:
        resolved = resolve_case(request, prompter=prompter, repo_root=repo_root)
        case_config_model(resolved)
        text = render_case_config(resolved)
        out_path = (out or repo_root / "config" / f"{resolved.case_id}.yaml").expanduser()
        _confirm_write(out_path, resolved, prompter=prompter, confirmed=confirmed, force=force)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(text, encoding="utf-8")
        # Read it back through the ordinary loader, so what is reported as valid
        # is the file on disk and its merge with --defaults, not the values in
        # memory a moment ago.
        load_case_config_with_defaults(out_path, defaults)
    except (MvaError, ValidationError) as exc:
        _fail(_explain(exc))
        return

    _print_case_summary(resolved, out_path, repo_root)
    config_arg = _display_path(out_path, repo_root=repo_root, workspace=resolved.workspace)
    defaults_arg = (
        None
        if defaults is None
        else _display_path(defaults, repo_root=repo_root, workspace=resolved.workspace)
    )
    if not run:
        console.print("\n  next:")
        if not _workspace_env_set():
            console.print(
                "    [yellow]note[/yellow] MVA_WORKSPACE is not set in this shell. Export "
                "it to the workspace you just gave (it is not echoed here), or pass "
                "--workspace again."
            )
        console.print(f"    {_next_command(config_arg, defaults_arg)}")
        return
    console.print("")
    run_all(
        config=out_path,
        defaults=defaults,
        workspace=resolved.workspace.root,
        allow_workspace_in_repo=False,
    )


def _confirm_write(
    out_path: Path,
    resolved: ResolvedCase,
    *,
    prompter: Prompter | None,
    confirmed: bool,
    force: bool,
) -> None:
    """Ask before writing, and before overwriting. Never ask when nobody can answer."""
    if out_path.exists() and not force:
        if prompter is None or confirmed:
            msg = (
                "The output file already exists. Pass --force to overwrite it, or --out to "
                "write somewhere else. Nothing was changed."
            )
            raise ConfigError(msg)
        answer = _read_answer(prompter, "Output file exists. Overwrite? [y/N]", "N")
        if answer.strip().lower() not in {"y", "yes"}:
            msg = "Declined to overwrite the existing file. Nothing was changed."
            raise ConfigError(msg)
    if prompter is None or confirmed:
        return
    console.print("")
    console.print(f"  case {resolved.case_id} / proband {resolved.proband_id}")
    console.print(f"  build {resolved.genome_build.value}, synthetic={resolved.synthetic}")
    console.print(f"  vcf {resolved.vcf}, phenotype {resolved.phenotype} (workspace-relative)")
    answer = _read_answer(prompter, "Write this case configuration? [Y/n]", "Y")
    if answer.strip().lower() not in {"y", "yes", ""}:
        msg = "Declined. Nothing was written."
        raise ConfigError(msg)


@case_app.command("check")
def case_check(
    config: ConfigOpt,
    defaults: DefaultsOpt = None,
    workspace: WorkspaceOpt = None,
) -> None:
    """Say whether a case config could run, without running anything.

    A dry run in the strict sense: it parses the config, resolves the workspace
    and asks whether each declared input is where the config says it is. It does
    not open an input, does not hash one, and does not create a run directory —
    `mva run validate` is the next step up, and it does all three.

    Exits 1 if anything is wrong. 'note' rows are reports, not failures.
    """
    rows: list[CheckRow] = []
    case: CaseConfig | None = None
    try:
        case = load_case_config_with_defaults(config, defaults)
    except (MvaError, ValidationError) as exc:
        rows.append(CheckRow("config parses", "FAIL", _explain(exc)))

    if case is not None:
        rows.extend(_check_config_rows(case))
        resolved_ws, workspace_row = _check_workspace_row(workspace)
        rows.append(workspace_row)
        rows.extend(_check_input_rows(case, resolved_ws))

    _render_check_rows(rows)
    if any(row.state == "FAIL" for row in rows):
        err_console.print(
            "\n[bold red]case NOT runnable[/bold red] — fix the FAIL rows above. "
            "`mva case init` writes a config that passes this check."
        )
        raise typer.Exit(code=1)
    console.print("\n[green]ok[/green] the case config is complete and its inputs are present")


def _check_config_rows(case: CaseConfig) -> list[CheckRow]:
    """Rows for the parts of the case the config file alone decides."""
    rows = [
        CheckRow(
            "required fields",
            "ok",
            "case_id, proband_id, genome_build, synthetic, inputs.vcf, inputs.phenotype",
        ),
        CheckRow("case id", "ok", case.case_id),
    ]
    if case.proband_id == EXPECTED_PROBAND_ID:
        rows.append(CheckRow("proband id", "ok", case.proband_id))
    else:
        rows.append(
            CheckRow(
                "proband id",
                "note",
                f"{case.proband_id} — the scorer accepts only {EXPECTED_PROBAND_ID}",
            )
        )
    if case.genome_build is GenomeBuild.GRCH38:
        rows.append(CheckRow("genome build", "ok", case.genome_build.value))
    else:
        rows.append(
            CheckRow("genome build", "note", f"{case.genome_build.value} — the challenge is GRCh38")
        )
    rows.append(
        CheckRow(
            "synthetic",
            "ok",
            "true — privacy checks are relaxed" if case.synthetic else "false — real-case rules",
        )
    )
    return rows


def _check_workspace_row(workspace: Path | None) -> tuple[Workspace | None, CheckRow]:
    try:
        resolved = resolve_workspace(workspace)
    except MvaError as exc:
        return None, CheckRow("workspace", "FAIL", _explain(exc))
    return resolved, CheckRow("workspace", "ok", "resolved; outside the repo, not cloud-synced")


def _check_input_rows(case: CaseConfig, workspace: Workspace | None) -> list[CheckRow]:
    """One row per declared input: is it there, and is its name plausible?"""
    declared: tuple[tuple[str, str, tuple[str, ...]], ...] = (
        ("inputs.vcf", case.inputs.vcf, VCF_SUFFIXES),
        ("inputs.phenotype", case.inputs.phenotype, PHENOTYPE_SUFFIXES),
    )
    rows: list[CheckRow] = []
    for name, relative, suffixes in declared:
        if workspace is None:
            rows.append(CheckRow(name, "FAIL", f"{relative} — not checkable without a workspace"))
            continue
        try:
            present = workspace.path(relative).is_file()
        except MvaError as exc:
            rows.append(CheckRow(name, "FAIL", _explain(exc)))
            continue
        if not present:
            rows.append(CheckRow(name, "FAIL", f"{relative} — no file at that workspace path"))
        elif not _matches_suffix(relative, suffixes):
            rows.append(
                CheckRow(name, "note", f"{relative} — expected one of {', '.join(suffixes)}")
            )
        else:
            rows.append(CheckRow(name, "ok", relative))
    return rows


def _render_check_rows(rows: list[CheckRow]) -> None:
    table = Table(show_header=True, title="case check")
    table.add_column("check")
    table.add_column("result")
    table.add_column("detail")
    states = {
        "ok": "[green]ok[/green]",
        "FAIL": "[bold red]FAIL[/bold red]",
        "note": "[yellow]note[/yellow]",
    }
    for row in rows:
        table.add_row(row.name, states.get(row.state, row.state), row.detail)
    console.print(table)


@case_app.command("phenotype-template")
def case_phenotype_template(
    relative_path: Annotated[
        str,
        typer.Argument(help="Workspace-relative destination."),
    ] = DEFAULT_PHENOTYPE_TEMPLATE,
    workspace: WorkspaceOpt = None,
    force: Annotated[bool, typer.Option("--force", help="Overwrite an existing file.")] = False,
) -> None:
    """Write a phenotype TSV with the correct header, ready to fill in.

    The columns and the accepted `status` spellings come from
    `mva.phenotype.loader`, so the template cannot describe a schema the parser
    does not implement — which is the point, since the loader refuses an
    unrecognised status rather than defaulting it (GP-14).

    It is written INTO the workspace, through `resolve_workspace`, because you
    are about to fill it with patient phenotype data and that never enters the
    repo tree (GP-40).
    """
    try:
        resolved = resolve_workspace(workspace)
        if Path(relative_path).is_absolute():
            msg = (
                "The destination is workspace-relative, e.g. 'inputs/phenotype.tsv'. Point "
                "--workspace at the workspace root instead of giving an absolute path."
            )
            raise ConfigError(msg)
        destination = resolved.path(relative_path)
        if destination.exists() and not force:
            msg = f"{relative_path} already exists in the workspace. Pass --force to overwrite."
            raise ConfigError(msg)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(phenotype_template_text(), encoding="utf-8")
    except MvaError as exc:
        _fail(str(exc))
        return

    console.print(f"[green]ok[/green] wrote {relative_path} in the workspace (path not echoed)")
    console.print("  fill one row per HPO term, then pass it to `mva case init --phenotype`")


def main() -> None:
    try:
        app()
    except PrivacyViolationError as exc:
        err_console.print(f"[bold red]privacy violation[/bold red] {exc}")
        sys.exit(2)


if __name__ == "__main__":
    main()
