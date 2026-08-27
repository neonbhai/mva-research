"""Structural tests for the Snakemake workflow layer (ADR 0001, ADR 0006).

The workflow file is the one place in this repository that is neither type-checked
nor covered by the import-layering lints: `Snakefile` and `workflow/rules/*.smk`
are Python-with-a-preprocessor, so pyright and ruff never see them. These tests
are the substitute — most importantly the one that proves no rule can write an
artifact into the repository.

Same style as `tests/unit/test_architecture.py`: custom lints whose failure
messages carry their own remediation.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO / "Snakefile"
RULES_DIR = REPO / "workflow" / "rules"
CASE_CONFIG = REPO / "config" / "synthetic-case.yaml"
CASE_SCHEMA = REPO / "config" / "schemas" / "case.schema.json"

try:  # jsonschema arrives as a Snakemake dependency; it is not required directly.
    from jsonschema import ValidationError
    from jsonschema import validate as jsonschema_validate
except ImportError:  # pragma: no cover - exercised only where snakemake is absent
    jsonschema_validate = None
    ValidationError = None


def _snakemake_available() -> bool:
    return (
        shutil.which("snakemake") is not None
        or (Path(sys.executable).parent / "snakemake").exists()
    )


# ---------------------------------------------------------------------------
# 1. The workflow parses
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_snakefile_parses() -> None:
    """`snakemake --list` exercises the whole file: every rule, every include.

    `--list` is used rather than `--dryrun` on purpose. A dry run additionally
    requires the input files to exist, which means it requires a workspace — and
    the workspace is external, absolute and not present in CI (ADR 0006). Parsing
    is the part that must hold everywhere.
    """
    if not _snakemake_available():
        pytest.skip(
            "snakemake is not installed in this environment. It is an optional extra "
            "(`uv sync --extra workflow`) because the identical DAG runs without it "
            "via `mva run all` (ADR 0001)."
        )

    result = subprocess.run(  # noqa: S603 - fixed argv, no shell, no user input
        [sys.executable, "-m", "snakemake", "--configfile", str(CASE_CONFIG), "--list"],
        cwd=REPO,
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )

    assert result.returncode == 0, (
        "The Snakefile does not parse.\n"
        f"  exit code: {result.returncode}\n"
        f"  stderr:\n{result.stderr.strip()}\n"
        "\nRemediation: reproduce with\n"
        "    uv run snakemake --configfile config/synthetic-case.yaml --list\n"
        "A Snakefile is Python with a preprocessor, so this is the only thing that "
        "type-checks it at all — ruff and pyright are configured over src/ and tests/ "
        "and never see it. Common causes: a rule referencing a constant defined below "
        "it (the file is executed top to bottom), an `include:` path that does not "
        "exist, or a config key read at parse time that the case file does not define."
    )


# ---------------------------------------------------------------------------
# 2. Every rules file is actually wired in
# ---------------------------------------------------------------------------

_INCLUDE = re.compile(r'^\s*include:\s*["\']([^"\']+)["\']', re.MULTILINE)


@pytest.mark.unit
def test_every_rules_file_is_included() -> None:
    """A .smk file nobody includes is dead code that still looks authoritative."""
    included = {
        (REPO / target).resolve() for target in _INCLUDE.findall(SNAKEFILE.read_text("utf-8"))
    }
    present = {path.resolve() for path in sorted(RULES_DIR.glob("*.smk"))}

    orphans = sorted(path.relative_to(REPO).as_posix() for path in present - included)
    assert not orphans, (
        "Rules files that the Snakefile never includes:\n"
        + "\n".join(f"  {name}" for name in orphans)
        + '\n\nRemediation: add `include: "workflow/rules/<name>.smk"` to the '
        "Snakefile, or delete the file. An orphaned .smk is worse than a missing one: "
        "it reads as part of the pipeline, is found by grep, and is cited in review, "
        "while contributing nothing to the DAG. Note that includes are executed in "
        "order, so a file must come after the constants its rules reference."
    )

    dangling = sorted(
        path.relative_to(REPO).as_posix() for path in included - present if not path.is_file()
    )
    assert not dangling, (
        "The Snakefile includes files that do not exist:\n"
        + "\n".join(f"  {name}" for name in dangling)
        + "\n\nRemediation: create the file or remove the `include:` line."
    )


# ---------------------------------------------------------------------------
# 3. No rule writes inside the repository (ADR 0006 / GP-40)
# ---------------------------------------------------------------------------

_ASSIGNMENT = re.compile(r'^([A-Z][A-Z0-9_]*)\s*=\s*f?"([^"]*)"', re.MULTILINE)
_INTERPOLATION = re.compile(r"\{([A-Za-z_][A-Za-z0-9_]*)")
_IDENTIFIER = re.compile(r"\b([A-Z][A-Z0-9_]{2,})\b")

ADR_0006_REMEDIATION = (
    "\n\nADR 0006 remediation: every output path must be rooted at the external "
    "workspace, i.e. built from WORKSPACE or from a constant derived from it "
    "(RUN_DIR, INGEST_DIR, REPORT_DIR, ...). Define a new constant next to the "
    "others in the Snakefile's 'Artifact layout' block and reference it by name "
    "from the rule.\n"
    "Why this is a lint and not a convention: the convenient layout puts data next "
    "to code (./data/, ./runs/), and that is exactly how patient data ends up in "
    "git history, where `rm` does not remove it and a later push discloses it "
    "permanently. The workspace is an external absolute path from MVA_WORKSPACE, "
    "validated by mva.config.resolve_workspace(), which also rejects cloud-synced "
    "roots. A single hardcoded relative output in a rule defeats all of that, and "
    "would do so silently on the first real case rather than on the synthetic one."
)


def _workspace_rooted_names(text: str) -> set[str]:
    """Constants whose value is transitively built from WORKSPACE."""
    names = {"WORKSPACE"}
    assignments = _ASSIGNMENT.findall(text)
    changed = True
    while changed:
        changed = False
        for name, value in assignments:
            if name in names:
                continue
            if set(_INTERPOLATION.findall(value)) & names:
                names.add(name)
                changed = True
    return names


def _output_entries(text: str) -> list[tuple[int, str]]:
    """Every line inside an `output:` block, with its line number.

    Block extent is taken from indentation, which is how Snakemake itself reads a
    rule body.
    """
    entries: list[tuple[int, str]] = []
    lines = text.splitlines()
    index = 0
    while index < len(lines):
        line = lines[index]
        stripped = line.strip()
        if stripped.startswith("output:"):
            indent = len(line) - len(line.lstrip())
            inline = stripped[len("output:") :].strip()
            if inline:
                entries.append((index + 1, inline))
            index += 1
            while index < len(lines):
                body = lines[index]
                if not body.strip():
                    index += 1
                    continue
                body_indent = len(body) - len(body.lstrip())
                if body_indent <= indent:
                    break
                entries.append((index + 1, body.strip()))
                index += 1
            continue
        index += 1
    return entries


@pytest.mark.unit
def test_no_rule_writes_inside_the_repository() -> None:
    """Every declared output resolves under the external workspace."""
    snakefile_text = SNAKEFILE.read_text("utf-8")
    rooted = _workspace_rooted_names(snakefile_text)

    assert len(rooted) > 1, (
        "No workspace-rooted constants found in the Snakefile.\n\n"
        "Remediation: output paths are expected to be named constants built from "
        'WORKSPACE, e.g. RUN_DIR = f"{WORKSPACE}/runs/{RUN_ID}". If that layout '
        "changed, update _workspace_rooted_names() here in the same commit — but do "
        "not weaken the rule." + ADR_0006_REMEDIATION
    )

    violations: list[str] = []
    for path in [SNAKEFILE, *sorted(RULES_DIR.glob("*.smk"))]:
        rel = path.relative_to(REPO).as_posix()
        for lineno, entry in _output_entries(path.read_text("utf-8")):
            names = set(_IDENTIFIER.findall(entry))
            if names & rooted:
                continue
            violations.append(f"  {rel}:{lineno} — output: {entry}")

    assert not violations, (
        "Rule outputs that are not rooted at the external workspace.\n"
        + "\n".join(violations)
        + ADR_0006_REMEDIATION
    )


@pytest.mark.unit
def test_no_output_is_a_literal_relative_path() -> None:
    """A quoted output path must interpolate, never be a bare literal."""
    offenders: list[str] = []
    for path in [SNAKEFILE, *sorted(RULES_DIR.glob("*.smk"))]:
        rel = path.relative_to(REPO).as_posix()
        for lineno, entry in _output_entries(path.read_text("utf-8")):
            for literal in re.findall(r'"([^"]*)"', entry):
                if "{" not in literal and "/" in literal:
                    offenders.append(f'  {rel}:{lineno} — output: "{literal}"')

    assert not offenders, (
        "Hardcoded literal output paths found.\n" + "\n".join(offenders) + ADR_0006_REMEDIATION
    )


# ---------------------------------------------------------------------------
# 4. The case schema is valid and the synthetic case satisfies it
# ---------------------------------------------------------------------------

REQUIRED_CASE_KEYS = ("case_id", "proband_id", "genome_build", "synthetic", "inputs")


def _load_schema() -> dict[str, Any]:
    raw: Any = json.loads(CASE_SCHEMA.read_text("utf-8"))
    assert isinstance(raw, dict), (
        f"{CASE_SCHEMA.relative_to(REPO)} must contain a JSON object at the top level.\n\n"
        "Remediation: a JSON Schema document is an object. Restore the outer braces."
    )
    return raw


def _load_case() -> dict[str, Any]:
    raw: Any = yaml.safe_load(CASE_CONFIG.read_text("utf-8"))
    assert isinstance(raw, dict), (
        f"{CASE_CONFIG.relative_to(REPO)} must be a YAML mapping at the top level.\n\n"
        "Remediation: the case config is a mapping of settings; see "
        "src/mva/config.py::load_case_config, which raises the same error."
    )
    return raw


@pytest.mark.unit
def test_case_schema_is_valid_json() -> None:
    """The schema parses and declares the fields CaseConfig makes mandatory."""
    schema = _load_schema()
    required = set(schema.get("required", []))
    missing = [key for key in REQUIRED_CASE_KEYS if key not in required]

    assert not missing, (
        "config/schemas/case.schema.json does not require: " + ", ".join(missing) + "\n\n"
        "Remediation: these five keys have no default in "
        "src/mva/config.py::CaseConfig, so a case file that omits one is invalid, not "
        "merely under-specified. `synthetic` matters most: it is explicit and never "
        "inferred, because a real case mislabelled synthetic would skip privacy "
        "enforcement, and the failure mode of a forgotten key must be a loud error "
        "rather than the unsafe branch."
    )


@pytest.mark.unit
def test_synthetic_case_validates_against_the_schema() -> None:
    """config/synthetic-case.yaml is the case every demo and golden test runs."""
    schema = _load_schema()
    case = _load_case()

    if jsonschema_validate is not None and ValidationError is not None:
        try:
            jsonschema_validate(instance=case, schema=schema)
        except ValidationError as exc:
            path = "/".join(str(part) for part in exc.absolute_path) or "<root>"
            pytest.fail(
                "config/synthetic-case.yaml does not satisfy "
                "config/schemas/case.schema.json.\n"
                f"  at: {path}\n"
                f"  problem: {exc.message}\n"
                "\nRemediation: fix whichever of the two is wrong — and decide which "
                "by reading src/mva/config.py::CaseConfig, which is the authority. The "
                "schema exists so a bad case file fails BEFORE Snakemake builds the "
                "DAG and before any patient file is opened; it is allowed to be weaker "
                "than the Pydantic model but never wider on the workspace-relative "
                "path rule or the 10-row submission cap."
            )
        return

    # jsonschema absent: fall back to a structural check rather than adding a
    # dependency for one test. Weaker, but it still catches a renamed key.
    missing = [key for key in REQUIRED_CASE_KEYS if key not in case]
    assert not missing, (
        "config/synthetic-case.yaml is missing required keys: "
        + ", ".join(missing)
        + "\n\nRemediation: add them. (jsonschema is not installed here, so only the "
        "presence of required keys was checked, not their types. Install the workflow "
        "extra — `uv sync --extra workflow` — for the full schema check.)"
    )


@pytest.mark.unit
def test_schema_forbids_absolute_and_escaping_input_paths() -> None:
    """The workspace-relative rule is encoded in the schema, not just in Pydantic."""
    if jsonschema_validate is None or ValidationError is None:
        pytest.skip("jsonschema is not installed; this check needs a real validator.")

    schema = _load_schema()
    case = _load_case()
    hostile = ["/etc/passwd", "../../outside.vcf", "~/Documents/case.vcf"]

    accepted: list[str] = []
    for candidate in hostile:
        probe = dict(case)
        inputs = dict(case["inputs"])
        inputs["vcf"] = candidate
        probe["inputs"] = inputs
        try:
            jsonschema_validate(instance=probe, schema=schema)
        except ValidationError:
            continue
        accepted.append(candidate)

    assert not accepted, (
        "config/schemas/case.schema.json accepts input paths it must reject:\n"
        + "\n".join(f"  {candidate}" for candidate in accepted)
        + "\n\nRemediation: restore the pattern on $defs/workspaceRelativePath so that "
        "a leading '/', a leading '~' and any '..' segment are refused. Input paths are "
        "workspace-relative so that no committed config and no provenance manifest "
        "records where patient data lives on disk (ADR 0006). An absolute path in a "
        "case file is not just a validation miss: it is a disclosure the moment the "
        "file is committed or the manifest is published alongside a submission."
    )
