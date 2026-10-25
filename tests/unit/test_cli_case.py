"""Tests for `mva case` — the guided front door (`src/mva/cli.py`).

Three properties are worth more than the rest of this file put together, because
each of them is a thing a future edit could break silently:

* **The two modes cannot disagree.** A guided run and a scripted run share one
  resolver, so for identical answers they must write identical bytes. If that
  stops being true, the interactive path has grown a rule of its own.
* **A missing value is an error, never a wait.** `--no-input` is a promise to a
  script, a CI job and an agent: a hang there is not a slow failure, it is a
  stuck pipeline nobody is watching.
* **Nothing opens an input file.** Asserted by making the file unreadable and
  requiring the command to succeed anyway (GP-44, ADR 0008). A test that merely
  checks the output has no genotype in it would pass for a command that read the
  whole VCF and said nothing about it.
"""

from __future__ import annotations

import stat
from pathlib import Path

import pytest
from typer.testing import CliRunner

from mva.cli import (
    EXPECTED_PROBAND_ID,
    ResolvedCase,
    phenotype_template_text,
    render_case_config,
)
from mva.cli import app as cli_app
from mva.config import CaseConfig, Workspace, load_case_config
from mva.models.genome import GenomeBuild
from mva.models.phenotype import ObservationStatus, Onset

pytestmark = [pytest.mark.unit]

#: Rich wraps to the terminal width, which would break substring assertions on a
#: narrow one. Also keeps colour escapes out of the captured text.
WIDE_TERMINAL: dict[str, str] = {"COLUMNS": "400", "TERM": "dumb", "NO_COLOR": "1"}

VCF_RELATIVE = "inputs/synthetic_case.vcf"
PHENOTYPE_RELATIVE = "inputs/synthetic_phenotype.tsv"


def _flags(workspace: Workspace, out: Path, *, case_id: str = "unit-case") -> list[str]:
    """Every required value as a flag: the shape a script or an agent uses."""
    return [
        "case",
        "init",
        "--no-input",
        "--workspace",
        str(workspace.root),
        "--vcf",
        VCF_RELATIVE,
        "--phenotype",
        PHENOTYPE_RELATIVE,
        "--genome-build",
        "GRCh38",
        "--case-id",
        case_id,
        "--proband-id",
        EXPECTED_PROBAND_ID,
        "--out",
        str(out),
    ]


def _flat(text: str) -> str:
    """Collapse whitespace so a wrapped message still matches a substring."""
    return " ".join(text.split())


def _invoke(args: list[str], *, stdin: str | None = None) -> tuple[int, str]:
    """Run the CLI and return the exit code with stdout and stderr together.

    The two streams are joined because what matters in these tests is whether a
    reason reached the operator at all, not which stream carried it.
    """
    result = CliRunner().invoke(cli_app, args, input=stdin, env=WIDE_TERMINAL)
    return result.exit_code, _flat(result.stdout + result.stderr)


# --------------------------------------------------------------------- shapes


def test_non_interactive_writes_a_config_the_loader_accepts(
    synthetic_workspace: Workspace, tmp_path: Path
) -> None:
    out = tmp_path / "written" / "case.yaml"

    exit_code, output = _invoke(_flags(synthetic_workspace, out))

    assert exit_code == 0, output
    case = load_case_config(out)
    assert case.case_id == "unit-case"
    assert case.proband_id == EXPECTED_PROBAND_ID
    assert case.genome_build is GenomeBuild.GRCH38
    assert case.synthetic is False, "a case is real unless --synthetic says otherwise"
    assert case.inputs.vcf == VCF_RELATIVE
    assert case.inputs.phenotype == PHENOTYPE_RELATIVE


def test_the_config_records_no_absolute_path_and_no_workspace_key(
    synthetic_workspace: Workspace, tmp_path: Path
) -> None:
    """ADR 0006: the workspace is passed at run time, never committed in a case file.

    `CaseConfig` forbids a `workspace` key outright, so writing one would make
    the file unloadable — but the absolute paths are the real hazard, because
    they would load fine and record where the patient data lives.
    """
    out = tmp_path / "case.yaml"

    exit_code, output = _invoke(_flags(synthetic_workspace, out))

    assert exit_code == 0, output
    text = out.read_text(encoding="utf-8")
    assert str(synthetic_workspace.root) not in text
    assert "workspace:" not in text
    for line in text.splitlines():
        assert not line.strip().startswith("/"), line


def test_rendered_config_is_exactly_this(synthetic_workspace: Workspace) -> None:
    """Lock the bytes. Key order and quoting are properties of one function (GP-30)."""
    resolved = ResolvedCase(
        workspace=synthetic_workspace,
        case_id="locked-case",
        proband_id=EXPECTED_PROBAND_ID,
        genome_build=GenomeBuild.GRCH38,
        synthetic=False,
        description="",
        vcf=VCF_RELATIVE,
        phenotype=PHENOTYPE_RELATIVE,
    )

    rendered = render_case_config(resolved)

    assert rendered.endswith("\n")
    body = [line for line in rendered.splitlines() if not line.startswith("#")]
    assert body == [
        'case_id: "locked-case"',
        'proband_id: "PROBAND01"',
        "genome_build: GRCh38",
        "synthetic: false",
        "",
        "inputs:",
        '  vcf: "inputs/synthetic_case.vcf"',
        '  phenotype: "inputs/synthetic_phenotype.tsv"',
    ]


def test_both_modes_write_identical_bytes(synthetic_workspace: Workspace, tmp_path: Path) -> None:
    """The property that keeps the guided path from growing rules of its own."""
    scripted = tmp_path / "scripted.yaml"
    guided = tmp_path / "guided.yaml"
    answers = "\n".join(
        [
            str(synthetic_workspace.root),
            VCF_RELATIVE,
            PHENOTYPE_RELATIVE,
            "",  # accept the offered GRCh38
            "unit-case",
            "",  # accept the offered PROBAND01
            "y",  # confirm the write
        ]
    )

    scripted_code, scripted_output = _invoke(_flags(synthetic_workspace, scripted))
    guided_code, guided_output = _invoke(
        ["case", "init", "--interactive", "--out", str(guided)], stdin=answers
    )

    assert scripted_code == 0, scripted_output
    assert guided_code == 0, guided_output
    assert guided.read_bytes() == scripted.read_bytes()


# --------------------------------------------------------------------- refusals


def test_no_input_with_a_missing_value_fails_and_names_every_missing_flag(
    synthetic_workspace: Workspace, tmp_path: Path
) -> None:
    """The promise `--no-input` makes: an error, not a question, and not a wait.

    All the missing flags are named at once because the caller is a script, which
    cannot iterate through six round trips as cheaply as a person answers six
    questions.
    """
    out = tmp_path / "never.yaml"

    exit_code, output = _invoke(
        [
            "case",
            "init",
            "--no-input",
            "--workspace",
            str(synthetic_workspace.root),
            "--vcf",
            VCF_RELATIVE,
            "--out",
            str(out),
        ]
    )

    assert exit_code == 1
    assert "--phenotype" in output
    assert "--genome-build" in output
    assert "--case-id" in output
    assert "--proband-id" in output
    assert not out.exists(), "a failed init must not leave a partial file"


def test_a_pipe_is_not_a_terminal_so_nothing_is_consumed_from_it(
    synthetic_workspace: Workspace, tmp_path: Path
) -> None:
    """Without a TTY and without --interactive, answers on stdin are not read.

    This is what stops a CI job from silently feeding a log file into six
    prompts. `--interactive` is the explicit opt in.
    """
    out = tmp_path / "never.yaml"

    exit_code, output = _invoke(
        ["case", "init", "--out", str(out)],
        stdin=f"{synthetic_workspace.root}\n{VCF_RELATIVE}\n{PHENOTYPE_RELATIVE}\n\nc\n\ny\n",
    )

    assert exit_code == 1
    assert "--workspace" in output
    assert not out.exists()


def test_a_workspace_inside_the_repo_is_refused(repo_root: Path, tmp_path: Path) -> None:
    """GP-40, through the existing validator. There is deliberately no bypass flag."""
    out = tmp_path / "never.yaml"

    exit_code, output = _invoke(
        [
            "case",
            "init",
            "--no-input",
            "--workspace",
            str(repo_root),
            "--vcf",
            VCF_RELATIVE,
            "--phenotype",
            PHENOTYPE_RELATIVE,
            "--genome-build",
            "GRCh38",
            "--case-id",
            "in-repo",
            "--proband-id",
            EXPECTED_PROBAND_ID,
            "--out",
            str(out),
        ]
    )

    assert exit_code == 1
    assert "inside the repository" in output
    assert not out.exists()
    assert "--allow" not in output, "the CLI must not advertise a way around GP-40"


@pytest.mark.parametrize("raw_read", ["inputs/proband.bam", "inputs/x.cram", "inputs/r1.fastq.gz"])
def test_raw_reads_are_refused_by_name_with_the_reason(
    synthetic_workspace: Workspace, tmp_path: Path, raw_read: str
) -> None:
    """TD-03: no alignment and no SV/CNV calling, so a BAM is not an input.

    Refused on the name alone, before anything touches the filesystem — these
    paths do not exist, and the refusal must not depend on their existing.
    """
    out = tmp_path / "never.yaml"
    flags = _flags(synthetic_workspace, out)
    flags[flags.index("--vcf") + 1] = raw_read

    exit_code, output = _invoke(flags)

    assert exit_code == 1
    assert "TD-03" in output
    assert not out.exists()


def test_an_input_outside_the_workspace_is_refused(
    synthetic_workspace: Workspace, repo_root: Path, tmp_path: Path
) -> None:
    out = tmp_path / "never.yaml"
    flags = _flags(synthetic_workspace, out)
    flags[flags.index("--vcf") + 1] = str(
        repo_root / "tests" / "fixtures" / "synthetic" / "synthetic_case.vcf"
    )

    exit_code, output = _invoke(flags)

    assert exit_code == 1
    assert "outside the workspace" in output
    assert not out.exists()


def test_an_unknown_genome_build_is_refused_rather_than_assumed(
    synthetic_workspace: Workspace, tmp_path: Path
) -> None:
    """GP-11. `GenomeBuild.parse` owns this rule; the CLI only reports it."""
    out = tmp_path / "never.yaml"
    flags = _flags(synthetic_workspace, out)
    flags[flags.index("--genome-build") + 1] = "hg99"

    exit_code, output = _invoke(flags)

    assert exit_code == 1
    assert "Unrecognised genome build" in output
    assert not out.exists()


def test_a_build_alias_resolves_to_its_canonical_spelling(
    synthetic_workspace: Workspace, tmp_path: Path
) -> None:
    """`hg38` is a real spelling in real VCF headers; the config stores GRCh38."""
    out = tmp_path / "case.yaml"
    flags = _flags(synthetic_workspace, out)
    flags[flags.index("--genome-build") + 1] = "hg38"

    exit_code, output = _invoke(flags)

    assert exit_code == 0, output
    assert load_case_config(out).genome_build is GenomeBuild.GRCH38


def test_an_existing_file_is_not_overwritten_without_force(
    synthetic_workspace: Workspace, tmp_path: Path
) -> None:
    out = tmp_path / "case.yaml"
    out.write_text("# already here\n", encoding="utf-8")

    exit_code, output = _invoke(_flags(synthetic_workspace, out))

    assert exit_code == 1
    assert "--force" in output
    assert out.read_text(encoding="utf-8") == "# already here\n"


def test_force_overwrites(synthetic_workspace: Workspace, tmp_path: Path) -> None:
    out = tmp_path / "case.yaml"
    out.write_text("# already here\n", encoding="utf-8")

    exit_code, output = _invoke([*_flags(synthetic_workspace, out), "--force"])

    assert exit_code == 0, output
    assert load_case_config(out).case_id == "unit-case"


# --------------------------------------------------------------------- privacy


def test_no_input_file_is_ever_opened(synthetic_workspace: Workspace, tmp_path: Path) -> None:
    """The input is made UNREADABLE and the command must still succeed.

    `mva case init` validates a path and an extension; the bytes behind it belong
    to the pipeline (GP-44, ADR 0008). Permissions are the only way to assert
    that from outside: a command that opened the file would fail here, and one
    that merely declined to print what it read would not.
    """
    locked = synthetic_workspace.root / "inputs" / "locked.vcf"
    locked.write_text("", encoding="utf-8")
    locked.chmod(0o000)
    out = tmp_path / "case.yaml"
    flags = _flags(synthetic_workspace, out)
    flags[flags.index("--vcf") + 1] = "inputs/locked.vcf"

    try:
        # Positive control: the file really is unreadable to this process. Without
        # it the assertions below pass trivially wherever chmod does not bite
        # (running as root, or a filesystem mounted without permissions).
        with pytest.raises(PermissionError):
            locked.read_text(encoding="utf-8")

        init_code, init_output = _invoke(flags)
        check_code, check_output = _invoke(
            [
                "case",
                "check",
                "--config",
                str(out),
                "--workspace",
                str(synthetic_workspace.root),
            ]
        )
    finally:
        locked.chmod(stat.S_IRUSR | stat.S_IWUSR)

    assert init_code == 0, init_output
    assert check_code == 0, check_output


def test_a_rejected_path_is_not_echoed_back(synthetic_workspace: Workspace, tmp_path: Path) -> None:
    """A path under the workspace is itself sensitive, so a refusal withholds it.

    `mva.pipeline.validate_case` sets the precedent for the same reason. The
    positive control matters: without it this test would pass for a command that
    printed nothing at all.
    """
    sentinel = "CASE-VOLUME-7f3a"
    out = tmp_path / "never.yaml"
    flags = _flags(synthetic_workspace, out)
    flags[flags.index("--vcf") + 1] = f"inputs/{sentinel}/proband.vcf"

    exit_code, output = _invoke(flags)

    assert exit_code == 1
    assert "No file at the path given for vcf" in output, output
    assert sentinel not in output


def test_a_rejected_identifier_is_reported_without_pydantic_echoing_it(
    synthetic_workspace: Workspace, tmp_path: Path
) -> None:
    """Pydantic appends `input_value=...` by default; the CLI rebuilds the message.

    The identifier here is operator-typed rather than patient-derived, so this is
    belt and braces — but the same code path reports a rejected `InputPaths`
    value, and that one is a path into the workspace (PRIV-09).
    """
    sentinel = "PATIENT-ID-9c2b!"
    out = tmp_path / "never.yaml"
    flags = _flags(synthetic_workspace, out)
    flags[flags.index("--case-id") + 1] = sentinel

    exit_code, output = _invoke(flags)

    assert exit_code == 1
    assert "String should match pattern" in output
    assert sentinel not in output


# --------------------------------------------------------------------- guided


def test_an_invalid_answer_is_re_prompted_with_its_reason(
    synthetic_workspace: Workspace, repo_root: Path, tmp_path: Path
) -> None:
    """Re-prompting, not aborting: the reason is the whole value of the prompt."""
    out = tmp_path / "case.yaml"
    answers = "\n".join(
        [
            str(repo_root),  # refused: inside the repo
            str(synthetic_workspace.root),
            "inputs/not-there.vcf",  # refused: no such file
            VCF_RELATIVE,
            PHENOTYPE_RELATIVE,
            "",
            "guided-case",
            "",
            "y",
        ]
    )

    exit_code, output = _invoke(["case", "init", "--interactive", "--out", str(out)], stdin=answers)

    assert exit_code == 0, output
    assert "inside the repository" in output
    assert "No file at the path given for vcf" in output
    assert load_case_config(out).case_id == "guided-case"


def test_the_prompt_loop_gives_up_rather_than_looping_forever(
    synthetic_workspace: Workspace, tmp_path: Path
) -> None:
    """A bounded loop. Unbounded, a file of the wrong thing on stdin is a hang."""
    out = tmp_path / "never.yaml"

    exit_code, output = _invoke(
        ["case", "init", "--interactive", "--out", str(out)],
        stdin="/nowhere/at/all\n" * 40,
    )

    assert exit_code == 1
    assert "giving up" in output
    assert not out.exists()


def test_declining_the_confirmation_writes_nothing(
    synthetic_workspace: Workspace, tmp_path: Path
) -> None:
    out = tmp_path / "never.yaml"
    answers = "\n".join(
        [
            str(synthetic_workspace.root),
            VCF_RELATIVE,
            PHENOTYPE_RELATIVE,
            "",
            "declined-case",
            "",
            "n",
        ]
    )

    exit_code, output = _invoke(["case", "init", "--interactive", "--out", str(out)], stdin=answers)

    assert exit_code == 1
    assert "Nothing was written" in output
    assert not out.exists()


# --------------------------------------------------------------------- check


def test_check_passes_on_a_config_init_just_wrote(
    synthetic_workspace: Workspace, tmp_path: Path
) -> None:
    """The two halves of the front door must agree about what is valid."""
    out = tmp_path / "case.yaml"
    assert _invoke(_flags(synthetic_workspace, out))[0] == 0

    exit_code, output = _invoke(
        ["case", "check", "--config", str(out), "--workspace", str(synthetic_workspace.root)]
    )

    assert exit_code == 0, output
    assert "FAIL" not in output


def test_check_names_a_missing_input_and_exits_non_zero(
    synthetic_workspace: Workspace, tmp_path: Path
) -> None:
    config = tmp_path / "case.yaml"
    config.write_text(
        "case_id: gap-case\n"
        f"proband_id: {EXPECTED_PROBAND_ID}\n"
        "genome_build: GRCh38\n"
        "synthetic: false\n"
        "inputs:\n"
        '  vcf: "inputs/absent.vcf"\n'
        f'  phenotype: "{PHENOTYPE_RELATIVE}"\n',
        encoding="utf-8",
    )

    exit_code, output = _invoke(
        ["case", "check", "--config", str(config), "--workspace", str(synthetic_workspace.root)]
    )

    assert exit_code == 1
    assert "inputs/absent.vcf" in output
    assert "FAIL" in output


def test_check_reports_a_missing_required_field(tmp_path: Path) -> None:
    config = tmp_path / "case.yaml"
    config.write_text(
        "case_id: no-proband\ngenome_build: GRCh38\nsynthetic: false\n"
        f'inputs:\n  vcf: "{VCF_RELATIVE}"\n  phenotype: "{PHENOTYPE_RELATIVE}"\n',
        encoding="utf-8",
    )

    exit_code, output = _invoke(["case", "check", "--config", str(config)])

    assert exit_code == 1
    assert "proband_id" in output
    assert "config parses" in output


def test_check_notes_a_proband_id_the_scorer_will_reject(
    synthetic_workspace: Workspace, tmp_path: Path
) -> None:
    """A note, not a failure: the schema permits it and the scorer does not.

    Reported rather than refused, because narrowing `proband_id` is a change to
    `mva.config`, not to its front door.
    """
    out = tmp_path / "case.yaml"
    flags = _flags(synthetic_workspace, out)
    flags[flags.index("--proband-id") + 1] = "CHILD-7"
    assert _invoke(flags)[0] == 0

    exit_code, output = _invoke(
        ["case", "check", "--config", str(out), "--workspace", str(synthetic_workspace.root)]
    )

    assert exit_code == 0, output
    assert EXPECTED_PROBAND_ID in output
    assert "note" in output


def test_check_fails_when_no_workspace_is_configured(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, synthetic_workspace: Workspace
) -> None:
    monkeypatch.delenv("MVA_WORKSPACE", raising=False)
    out = tmp_path / "case.yaml"
    assert _invoke(_flags(synthetic_workspace, out))[0] == 0

    exit_code, output = _invoke(["case", "check", "--config", str(out)])

    assert exit_code == 1
    assert "No workspace configured" in output


# --------------------------------------------------------------------- template


def test_the_phenotype_template_parses_with_the_real_loader(
    synthetic_workspace: Workspace,
) -> None:
    """The template's whole purpose is to stop someone guessing the schema.

    So it is checked against the parser rather than against a copy of the header:
    the loader refuses an unrecognised `status` instead of defaulting it (GP-14),
    which makes a wrong template an error the user has to debug.
    """
    from mva.phenotype.loader import OPTIONAL_COLUMNS, REQUIRED_COLUMNS, load_phenotype_profile

    exit_code, output = _invoke(
        [
            "case",
            "phenotype-template",
            "inputs/template.tsv",
            "--workspace",
            str(synthetic_workspace.root),
        ]
    )

    assert exit_code == 0, output
    destination = synthetic_workspace.path("inputs/template.tsv")
    header = next(
        line
        for line in destination.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    )
    assert header.split("\t") == [*REQUIRED_COLUMNS, *OPTIONAL_COLUMNS]

    # The template carries NO observations: its example row is commented out, so
    # a forgotten template contributes nothing rather than a fabricated finding.
    blank = load_phenotype_profile(
        destination,
        subject_id="template",
        hpo_version="test",
        source_artifact="template",
    )
    assert blank.observations == ()

    # And the commented example is a CORRECT row, not decoration: uncomment it
    # and the loader accepts it. Without this the example could rot unnoticed,
    # which is exactly as bad as having no example.
    filled = destination.parent / "filled.tsv"
    filled.write_text(
        destination.read_text(encoding="utf-8").replace("\n#HP:", "\nHP:"), encoding="utf-8"
    )
    profile = load_phenotype_profile(
        filled, subject_id="template", hpo_version="test", source_artifact="template"
    )
    assert len(profile.observations) == 1
    assert profile.observations[0].status is ObservationStatus.OBSERVED
    assert profile.observations[0].onset is Onset.CONGENITAL


def test_the_template_text_names_every_accepted_status() -> None:
    """Read from the loader's own table, so the two cannot drift apart."""
    from mva.phenotype.loader import STATUS_ALIASES

    text = phenotype_template_text()

    for alias in STATUS_ALIASES:
        assert alias in text, alias
    assert "'not_assessed' is NOT" in text, "the template must state the GP-14 distinction"


def test_the_template_will_not_be_written_into_the_repo(repo_root: Path) -> None:
    """GP-40: it is about to be filled with patient phenotype data.

    Enforced by routing the destination through `resolve_workspace`, which is the
    same refusal the pipeline itself gets — not a second check written here.
    """
    exit_code, output = _invoke(
        ["case", "phenotype-template", "inputs/template.tsv", "--workspace", str(repo_root)]
    )

    assert exit_code == 1
    assert "inside the repository" in output
    assert not (repo_root / "inputs" / "template.tsv").exists()


def test_the_template_refuses_an_absolute_destination(
    synthetic_workspace: Workspace, tmp_path: Path
) -> None:
    exit_code, output = _invoke(
        [
            "case",
            "phenotype-template",
            str(tmp_path / "anywhere.tsv"),
            "--workspace",
            str(synthetic_workspace.root),
        ]
    )

    assert exit_code == 1
    assert "workspace-relative" in output
    assert not (tmp_path / "anywhere.tsv").exists()


# --------------------------------------------------------------------- wiring


def test_the_case_group_is_reachable_from_the_root_help() -> None:
    """A front door nobody can find is not a front door."""
    root_code, root_output = _invoke(["--help"])
    group_code, group_output = _invoke(["case", "--help"])

    assert root_code == 0
    assert "case Set up and run one case" in root_output, root_output
    assert group_code == 0
    for subcommand in ("init", "check", "phenotype-template"):
        assert subcommand in group_output, subcommand


def test_init_writes_a_config_the_pipeline_contract_accepts(
    synthetic_workspace: Workspace, tmp_path: Path
) -> None:
    """Round-trip through the typed model, including the defaults merge.

    `CaseConfig` is the authority on a case configuration, so the test of a
    writer is whether the model accepts what it wrote — not whether the YAML
    looks right.
    """
    out = tmp_path / "case.yaml"

    exit_code, output = _invoke(
        [*_flags(synthetic_workspace, out), "--synthetic", "--description", "unit test case"]
    )

    assert exit_code == 0, output
    case: CaseConfig = load_case_config(out)
    assert case.synthetic is True
    assert case.description == "unit test case"
    assert case.config_hash()
