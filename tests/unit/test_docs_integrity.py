"""Mechanical documentation checks — the anti-doc-rot lint.

Documentation in this repository is load-bearing: `CLAUDE.md` routes a reader to
`docs/`, lint failures cite `GP-nn`, threat controls cite `PRIV-nn`, and scientific
caveats cite `ASSUMPTION-AREA-nn`. A citation that resolves to nothing is worse
than no citation, because it reads as authority while pointing at nothing. These
tests make every cross-reference in the repository *mechanically* checkable, so
prose rots loudly instead of quietly.

Same style as `tests/unit/test_architecture.py`: text- and AST-walking custom
lints whose failure messages carry their own remediation, because the reader is
often an agent whose only view of the rule is that string.

What is deliberately NOT checked: whether the prose is *true*. That is a review
responsibility. These tests only prove that everything named exists.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Iterator
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
DOCS = REPO / "docs"
SRC = REPO / "src"
TESTS = REPO / "tests"
DECISIONS = DOCS / "decisions"

CLAUDE_MD = REPO / "CLAUDE.md"
JUSTFILE = REPO / "justfile"
GOLDEN_PRINCIPLES = DOCS / "golden-principles.md"
PRIVACY_MODEL = DOCS / "privacy-model.md"
SCIENTIFIC_ASSUMPTIONS = DOCS / "scientific-assumptions.md"
CURRENT_STATUS = DOCS / "current-status.md"
MATURITY_LEDGER = DOCS / "maturity-ledger.md"

#: Roots swept for citations. Broader than any single rule requires, on purpose:
#: a stale `GP-nn` is equally wrong wherever it is written.
CITATION_ROOTS: tuple[Path, ...] = (
    DOCS,
    SRC,
    TESTS,
    REPO / "workflow",
    REPO / "config",
)

#: Individual files swept alongside those roots.
CITATION_FILES: tuple[Path, ...] = (
    CLAUDE_MD,
    JUSTFILE,
    REPO / "Snakefile",
    REPO / "README.md",
)

#: Only text formats. Binary fixtures and Parquet artifacts have no citations.
TEXT_SUFFIXES = frozenset(
    {".py", ".md", ".smk", ".sql", ".yaml", ".yml", ".json", ".toml", ".txt", ".cfg", ""}
)

_INLINE_CODE = re.compile(r"`([^`\n]+)`")


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _iter_text_files() -> Iterator[Path]:
    """Every text file that could plausibly carry a citation."""
    for root in CITATION_ROOTS:
        if not root.is_dir():
            continue
        for path in sorted(root.rglob("*")):
            if not path.is_file() or "__pycache__" in path.parts:
                continue
            if path.suffix in TEXT_SUFFIXES:
                yield path
    for path in CITATION_FILES:
        if path.is_file():
            yield path


def _citations(pattern: re.Pattern[str]) -> dict[str, list[str]]:
    """Map every distinct identifier matched by ``pattern`` to where it is cited."""
    found: dict[str, list[str]] = {}
    for path in _iter_text_files():
        rel = path.relative_to(REPO).as_posix()
        for lineno, line in enumerate(_read(path).splitlines(), start=1):
            for match in pattern.finditer(line):
                found.setdefault(match.group(0), []).append(f"{rel}:{lineno}")
    return found


def _undefined(
    cited: dict[str, list[str]], defined: Iterable[str], *, limit_sites: int = 3
) -> list[str]:
    """Render one line per citation that resolves to nothing."""
    known = set(defined)
    lines: list[str] = []
    for identifier in sorted(cited):
        if identifier in known:
            continue
        sites = cited[identifier]
        shown = ", ".join(sites[:limit_sites])
        if len(sites) > limit_sites:
            shown = f"{shown} (+{len(sites) - limit_sites} more)"
        lines.append(f"  {identifier} — cited at {shown}")
    return lines


# ---------------------------------------------------------------------------
# 1. Every repo-relative path named in CLAUDE.md exists
# ---------------------------------------------------------------------------


def _looks_like_repo_path(token: str) -> bool:
    """Whether a backticked token is a claim about a file in this repository.

    Conservative by design. `mva.clock` is a module, `SYNTH*` is a gene prefix and
    `~/Documents` is not ours; none of those is a claim this repo can falsify.
    """
    if not token or token[0] in {"~", "$", "/", "."}:
        return False
    if token.startswith(("http://", "https://")):
        return False
    if any(char in token for char in " \t→<>|="):
        return False
    if "/" in token:
        return True
    return Path(token).suffix in {
        ".md",
        ".py",
        ".yaml",
        ".yml",
        ".toml",
        ".json",
        ".sql",
        ".tsv",
        ".csv",
    }


def _path_is_present(token: str) -> bool:
    cleaned = token.rstrip("/")
    if any(char in cleaned for char in "*?["):
        return any(REPO.glob(cleaned))
    return (REPO / cleaned).exists()


@pytest.mark.unit
def test_claude_md_paths_exist() -> None:
    """Every repo-relative path CLAUDE.md points at is really there.

    CLAUDE.md is a router. A router pointing at a file that was renamed sends
    every subsequent reader — human or agent — down a dead end, and the reader
    usually assumes the mistake is theirs.
    """
    text = _read(CLAUDE_MD)
    missing: list[str] = []
    for lineno, line in enumerate(text.splitlines(), start=1):
        for match in _INLINE_CODE.finditer(line):
            token = match.group(1).strip()
            if _looks_like_repo_path(token) and not _path_is_present(token):
                missing.append(f"  CLAUDE.md:{lineno} — `{token}`")

    assert not missing, (
        "CLAUDE.md references paths that do not exist.\n"
        + "\n".join(missing)
        + "\n\nRemediation: either create the file/directory, or update CLAUDE.md to "
        "name the path that actually exists. Do NOT delete the reference to make this "
        "pass if the document is genuinely still owed — CLAUDE.md is the router into "
        "docs/, and a router that quietly forgets a destination is how a whole "
        "document stops being read. If the path is not a repo path at all (a module "
        "name, a glob for something outside the tree), un-backtick it or spell it as a "
        "module path so this lint stops treating it as a file claim."
    )


# ---------------------------------------------------------------------------
# 2. Decision records are well-formed
# ---------------------------------------------------------------------------

_ADR_FILENAME = re.compile(r"^(\d{4})-[a-z0-9]+(?:-[a-z0-9]+)*\.md$")
_ADR_TITLE = re.compile(r"^# ADR (\d{4}) — \S")
_ADR_STATUS = re.compile(r"^\*\*Status:\*\*\s*\S")


@pytest.mark.unit
def test_decision_records_are_well_formed() -> None:
    """Every ADR has a conforming filename, a numbered title and a status."""
    problems: list[str] = []
    for path in sorted(DECISIONS.glob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(REPO).as_posix()
        name_match = _ADR_FILENAME.match(path.name)
        if name_match is None:
            problems.append(f"  {rel} — filename is not NNNN-lowercase-slug.md")
            continue

        lines = _read(path).splitlines()
        title_match = next((_ADR_TITLE.match(line) for line in lines[:5] if line), None)
        if title_match is None:
            problems.append(
                f"  {rel} — no title line of the form '# ADR {name_match.group(1)} "
                "— <what was decided>' in the first five lines"
            )
        elif title_match.group(1) != name_match.group(1):
            problems.append(
                f"  {rel} — title says ADR {title_match.group(1)} but the filename "
                f"says {name_match.group(1)}"
            )

        if not any(_ADR_STATUS.match(line) for line in lines[:10]):
            problems.append(f"  {rel} — no '**Status:**' line in the first ten lines")

    assert not problems, (
        "Malformed decision records.\n"
        + "\n".join(problems)
        + "\n\nRemediation: an ADR is addressed by its number from code comments, test "
        "messages and CLAUDE.md, so the number has to be in the filename AND in the "
        "title, and the two must agree. Copy the shape of "
        "docs/decisions/0001-workflow-engine.md:\n"
        "    filename: NNNN-short-lowercase-slug.md\n"
        "    line 1:   # ADR NNNN — What was decided\n"
        "    line 3:   **Status:** accepted · **Date:** YYYY-MM-DD\n"
        "The status line is what tells a reader whether they are looking at a live "
        "constraint or a superseded one; an ADR without it cannot be relied on or "
        "safely ignored."
    )


# ---------------------------------------------------------------------------
# 3. Every ADR number cited resolves to a file
# ---------------------------------------------------------------------------

_ADR_REFERENCE = re.compile(r"ADR (\d{4})|docs/decisions/(\d{4})")


@pytest.mark.unit
def test_adr_references_resolve() -> None:
    """`ADR 0006` and `docs/decisions/0006-*` both have to point at something."""
    numbers: dict[str, list[str]] = {}
    for path in _iter_text_files():
        rel = path.relative_to(REPO).as_posix()
        for lineno, line in enumerate(_read(path).splitlines(), start=1):
            for match in _ADR_REFERENCE.finditer(line):
                number = match.group(1) or match.group(2)
                numbers.setdefault(number, []).append(f"{rel}:{lineno}")

    existing = {
        match.group(1)
        for path in DECISIONS.glob("*.md")
        if (match := _ADR_FILENAME.match(path.name)) is not None
    }
    dangling = _undefined(numbers, existing)

    assert not dangling, (
        "Citations of decision records that do not exist.\n"
        + "\n".join(dangling)
        + "\n\nRemediation: write the ADR, or correct the number. A dangling ADR "
        "citation is the most expensive kind of doc rot here, because the whole point "
        "of the number is that a reader can look up WHY a constraint exists before "
        "deciding whether to change it. If the decision was made but never written "
        "down, that is the bug — write it down; the format is in "
        "docs/decisions/0001-workflow-engine.md."
    )


# ---------------------------------------------------------------------------
# 4-6. Every GP / PRIV / ASSUMPTION identifier cited is defined
# ---------------------------------------------------------------------------

_GP_REFERENCE = re.compile(r"\bGP-\d+\b")
_GP_DEFINITION = re.compile(r"^\*\*(GP-\d+) —")

_PRIV_REFERENCE = re.compile(r"\bPRIV-\d+\b")
_PRIV_DEFINITION = re.compile(r"^\|\s*(PRIV-\d+)\s*\|")

_ASSUMPTION_REFERENCE = re.compile(r"\bASSUMPTION-[A-Z]+-\d+\b")
_ASSUMPTION_DEFINITION = re.compile(r"^#{1,6}\s+(ASSUMPTION-[A-Z]+-\d+)\b")


def _defined_ids(path: Path, pattern: re.Pattern[str]) -> set[str]:
    return {
        match.group(1)
        for line in _read(path).splitlines()
        if (match := pattern.match(line)) is not None
    }


@pytest.mark.unit
def test_golden_principle_citations_resolve() -> None:
    """Every `GP-nn` cited anywhere is defined in docs/golden-principles.md."""
    defined = _defined_ids(GOLDEN_PRINCIPLES, _GP_DEFINITION)
    assert defined, (
        "No golden principles parsed out of docs/golden-principles.md.\n\n"
        "Remediation: definitions must be bold headings of the exact form "
        "'**GP-nn — One-line statement.**' at the start of a line. That shape is "
        "what makes the IDs machine-checkable; if the document has been restyled, "
        "restyle it back or update _GP_DEFINITION in this file in the same commit."
    )

    dangling = _undefined(_citations(_GP_REFERENCE), defined)
    assert not dangling, (
        "Undefined golden-principle citations.\n"
        + "\n".join(dangling)
        + "\n\nRemediation: a lint, a test message or a docstring that cites GP-nn is "
        "telling its reader 'there is a rule, and it is written down'. If the rule is "
        "real, add it to docs/golden-principles.md as "
        "'**GP-nn — statement.**' — which per that document also requires a "
        "decision record. If the citation was a typo, fix the number. Do not remove "
        "the citation to silence this: an unattributed rule is one nobody can argue "
        "with or repeal."
    )


@pytest.mark.unit
def test_privacy_control_citations_resolve() -> None:
    """Every `PRIV-nn` cited anywhere is a row in the docs/privacy-model.md threat table."""
    defined = _defined_ids(PRIVACY_MODEL, _PRIV_DEFINITION)
    assert defined, (
        "No privacy controls parsed out of docs/privacy-model.md.\n\n"
        "Remediation: controls are the leftmost column of the threat table, so a row "
        "must begin '| PRIV-nn |'. Restore that shape, or update _PRIV_DEFINITION here "
        "in the same commit that restyles the table."
    )

    dangling = _undefined(_citations(_PRIV_REFERENCE), defined)
    assert not dangling, (
        "Undefined privacy-control citations.\n"
        + "\n".join(dangling)
        + "\n\nRemediation: add the threat to the table in docs/privacy-model.md with "
        "its vector, likelihood, impact and control, or fix the number. The table is "
        "the threat model; a control that exists in code but not in the table has "
        "never been reasoned about against a threat, and a citation pointing at a "
        "missing row implies an analysis that was never done."
    )


@pytest.mark.unit
def test_scientific_assumption_citations_resolve() -> None:
    """Every `ASSUMPTION-AREA-nn` cited anywhere is defined in docs/scientific-assumptions.md."""
    defined = _defined_ids(SCIENTIFIC_ASSUMPTIONS, _ASSUMPTION_DEFINITION)
    assert defined, (
        "No assumptions parsed out of docs/scientific-assumptions.md.\n\n"
        "Remediation: each assumption is a markdown heading of the form "
        "'### ASSUMPTION-AREA-nn — What is assumed'. Restore that shape or update "
        "_ASSUMPTION_DEFINITION here in the same commit."
    )

    dangling = _undefined(_citations(_ASSUMPTION_REFERENCE), defined)
    assert not dangling, (
        "Undefined scientific-assumption citations.\n"
        + "\n".join(dangling)
        + "\n\nRemediation: document the assumption in docs/scientific-assumptions.md "
        "with its rationale and its consequence if wrong, or fix the identifier. Code "
        "that cites an assumption is claiming a scientific choice was made "
        "deliberately and is written down where a domain expert can disagree with it. "
        "An unresolvable citation turns that claim into a bluff."
    )


# ---------------------------------------------------------------------------
# 7. The status document exists and keeps its sections
# ---------------------------------------------------------------------------

REQUIRED_STATUS_SECTIONS: tuple[str, ...] = (
    "## What exists",
    "## What works",
    "## What is incomplete",
    "## Commands",
    "## Blockers",
)


@pytest.mark.unit
def test_current_status_has_required_sections() -> None:
    """docs/current-status.md is the hand-off document; its shape is the contract."""
    assert CURRENT_STATUS.is_file(), (
        "docs/current-status.md is missing.\n\n"
        "Remediation: create it with these headings, in this order:\n"
        + "\n".join(f"    {heading}" for heading in REQUIRED_STATUS_SECTIONS)
        + "\nCLAUDE.md routes 'current state, blockers, commands' here, and it is the "
        "first file a new session (human or agent) reads. Its value is entirely in "
        "being honest about the gap between what exists and what works, so 'What is "
        "incomplete' and 'Blockers' are required sections rather than optional ones."
    )

    lines = {line.rstrip() for line in _read(CURRENT_STATUS).splitlines()}
    missing = [heading for heading in REQUIRED_STATUS_SECTIONS if heading not in lines]
    assert not missing, (
        "docs/current-status.md is missing required sections:\n"
        + "\n".join(f"  {heading}" for heading in missing)
        + "\n\nRemediation: add the heading verbatim (exact case and wording). These "
        "five sections are separated on purpose: 'What exists' is a file listing and "
        "'What works' is a claim about behaviour, and collapsing them is precisely how "
        "a repo comes to look more finished than it is. An empty section with an "
        "honest 'nothing yet' underneath passes; a deleted section does not."
    )


# ---------------------------------------------------------------------------
# 8. Every package is graded in the maturity ledger
# ---------------------------------------------------------------------------


def _source_packages() -> list[str]:
    """Top-level packages under src/mva/ (directories with an __init__.py)."""
    root = SRC / "mva"
    return sorted(
        path.name
        for path in root.iterdir()
        if path.is_dir() and (path / "__init__.py").is_file() and path.name != "__pycache__"
    )


@pytest.mark.unit
def test_every_package_is_graded_in_the_maturity_ledger() -> None:
    """GP-20: every package is graded real / synthetic-substitute / stub."""
    packages = _source_packages()
    assert packages, (
        "No packages found under src/mva/.\n\n"
        "Remediation: this lint locates packages by looking for directories "
        "containing __init__.py. If the layout changed, update _source_packages() "
        "here in the same commit."
    )

    assert MATURITY_LEDGER.is_file(), (
        "docs/maturity-ledger.md is missing.\n\n"
        "Remediation: create it with one table row per package under src/mva/, "
        "grading each `real`, `synthetic-substitute` or `stub` (GP-20). Packages "
        "currently needing a row:\n"
        + "\n".join(f"    | `{name}` | stub | | | |" for name in packages)
        + "\nGrade honestly and downward when unsure: a synthetic substitute described "
        "as real is the one documentation error in this repo that can produce a false "
        "scientific claim."
    )

    ledger = _read(MATURITY_LEDGER)
    rows = [line for line in ledger.splitlines() if line.lstrip().startswith("|")]
    ungraded = [name for name in packages if not any(f"`{name}`" in row for row in rows)]

    assert not ungraded, (
        "Packages with no row in docs/maturity-ledger.md:\n"
        + "\n".join(f"  src/mva/{name}/" for name in ungraded)
        + "\n\nRemediation: add a table row naming the package in backticks and "
        "grading it `real`, `synthetic-substitute` or `stub`, with what IS real and "
        "what is NOT (GP-20). One row may cover several packages, e.g. "
        "'| `clock` / `determinism` | real | ... |'. This is not bookkeeping: reports "
        "surface the grade of everything they depend on, so an ungraded package is a "
        "claim rendered into a report with no maturity attached to it. If you cannot "
        "say what is real about it yet, `stub` is the correct and safe answer."
    )


# ---------------------------------------------------------------------------
# 9. Every `just` recipe named in CLAUDE.md exists
# ---------------------------------------------------------------------------

_JUST_INVOCATION = re.compile(r"\bjust\s+([a-z][a-z0-9-]*)")
#: A recipe header at column zero: `name:`, `name DEP:` or `name *ARGS:`. The
#: negative lookahead keeps just's own settings (`set shell := ...`) out of the
#: recipe set, since `just set` is not a runnable command.
_JUST_RECIPE = re.compile(
    r"^(?!set\s|alias\s|import\s|export\s)([a-z][a-z0-9-]*)\s*(?:[A-Za-z_*+$][^:]*)?:"
)


@pytest.mark.unit
def test_just_recipes_named_in_claude_md_exist() -> None:
    """CLAUDE.md's command list is copy-pasteable or it is worse than nothing."""
    recipes = {
        match.group(1)
        for line in _read(JUSTFILE).splitlines()
        if (match := _JUST_RECIPE.match(line)) is not None
    }
    assert recipes, (
        "No recipes parsed out of the justfile.\n\n"
        "Remediation: recipes are lines starting at column zero as 'name:' or "
        "'name ARGS:'. If the justfile layout changed, update _JUST_RECIPE here."
    )

    cited: dict[str, list[str]] = {}
    for lineno, line in enumerate(_read(CLAUDE_MD).splitlines(), start=1):
        for match in _JUST_INVOCATION.finditer(line):
            cited.setdefault(match.group(1), []).append(f"CLAUDE.md:{lineno}")

    missing = _undefined(cited, recipes)
    assert not missing, (
        "CLAUDE.md documents `just` recipes that the justfile does not define.\n"
        + "\n".join(missing)
        + "\n\nRemediation: add the recipe to the justfile, or correct the name in "
        "CLAUDE.md. The Commands block in CLAUDE.md is the first thing anyone runs; a "
        "command that errors with 'Justfile does not contain recipe' on line one costs "
        "the reader their trust in the rest of the document. Renaming a recipe is a "
        "two-file change, always."
    )
