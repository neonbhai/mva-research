"""Deterministic markdown rendering (GP-30).

Two rules govern everything in this module.

**Loud, not blank.** The Jinja environment uses ``StrictUndefined``. A report that
silently renders an empty string where a caveat should have been is worse than no
report at all: the reader cannot tell the difference between "no safety concern"
and "the safety field was misspelled in the template". A missing variable is
therefore an exception, not whitespace.

**No hidden state.** Output depends only on the template text and the context
passed in. There is no wall-clock read here — every timestamp arrives
pre-formatted from an injected :class:`~mva.clock.Clock` — and every iteration a
caller performs is over an already-ordered sequence. Two runs with the same
context produce byte-identical bytes, which is what makes
``tests/unit/test_reporting.py`` able to assert determinism at all.

Autoescaping is off on purpose: the output is markdown, not HTML, and escaping
``&`` into ``&amp;`` inside a gene name or a HGVS string would corrupt the very
identifiers a reader needs to copy. Nothing rendered here is ever served as HTML;
the artifacts are ``.md`` and ``.csv`` files.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from enum import Enum
from functools import lru_cache
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, StrictUndefined, TemplateNotFound

#: Repository ``templates/`` directory. ``render.py`` lives at
#: ``<root>/src/mva/reporting/render.py``, so the root is three parents up.
DEFAULT_TEMPLATES_DIR: Path = Path(__file__).resolve().parents[3] / "templates"

#: Rendered stand-in for a value that was never recorded. Deliberately wordy:
#: an em-dash reads as "not applicable", and GP-14 says absence of information is
#: not negative information, so the report says which one it means.
NOT_RECORDED = "not recorded"


@lru_cache(maxsize=8)
def _environment(templates_dir: Path) -> Environment:
    """Build (and cache) the Jinja environment for one template directory.

    Cached because template compilation is the expensive part and the environment
    is immutable in practice; keyed by directory so a test can point at a
    temporary directory without disturbing the production one.
    """
    return Environment(
        loader=FileSystemLoader(str(templates_dir), encoding="utf-8"),
        autoescape=False,  # noqa: S701 - markdown/CSV output; HTML escaping would corrupt HGVS
        trim_blocks=True,
        lstrip_blocks=True,
        keep_trailing_newline=True,
        undefined=StrictUndefined,
    )


def render_template(
    name: str,
    context: Mapping[str, object],
    *,
    templates_dir: Path | None = None,
) -> str:
    """Render one template with a strict, fully-supplied context.

    Raises ``jinja2.UndefinedError`` if the template references anything the
    context does not provide, and ``FileNotFoundError`` if the template directory
    or the template itself is missing. Both are deliberate: a report is a
    deliverable, and a half-rendered deliverable is a defect that must surface at
    build time rather than in a reviewer's hands.
    """
    directory = (templates_dir or DEFAULT_TEMPLATES_DIR).resolve()
    if not directory.is_dir():
        msg = (
            f"Template directory {directory.as_posix()!r} does not exist. Reporting "
            "templates live in the repository's templates/ directory; pass "
            "templates_dir= explicitly when rendering from an installed package."
        )
        raise FileNotFoundError(msg)
    try:
        template = _environment(directory).get_template(name)
    except TemplateNotFound as exc:
        available = sorted(p.name for p in directory.glob("*.j2"))
        msg = f"Template {name!r} not found in {directory.as_posix()!r}. Available: {available}."
        raise FileNotFoundError(msg) from exc
    return template.render(dict(context))


def format_cell(value: object) -> str:
    """Render one value for a markdown table, deterministically.

    Floats are fixed at three decimals rather than left to ``repr``: component
    scores are heuristics, and printing ``0.7000000000000001`` implies a precision
    the number does not have. Tri-state booleans keep their third state visible —
    ``None`` never renders as ``False`` (GP-14, GP-16).
    """
    if value is None:
        return NOT_RECORDED
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, Enum):
        return str(value.value)
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def _escape_cell(text: str) -> str:
    """Make a cell safe for a pipe table without dropping information."""
    return text.replace("\\", "\\\\").replace("|", "\\|").replace("\n", " ").replace("\r", " ")


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    """A GitHub-flavoured pipe table. Row order is the caller's; nothing is sorted here.

    Sorting is left to the caller on purpose: the ordering of a table is a
    scientific statement (rank order, karyotype order) and hiding it inside a
    formatting helper would make it invisible in review.
    """
    if not headers:
        msg = "markdown_table requires at least one header column."
        raise ValueError(msg)
    width = len(headers)
    header_line = "| " + " | ".join(_escape_cell(h) for h in headers) + " |"
    divider = "| " + " | ".join("---" for _ in headers) + " |"
    lines = [header_line, divider]
    for index, row in enumerate(rows):
        if len(row) != width:
            msg = (
                f"markdown_table row {index} has {len(row)} cells but the table has "
                f"{width} columns. A ragged table silently mis-aligns values under "
                "the wrong headings, which in a scores table means reporting the "
                "wrong number for the wrong component."
            )
            raise ValueError(msg)
        lines.append("| " + " | ".join(_escape_cell(format_cell(c)) for c in row) + " |")
    return "\n".join(lines)


__all__ = [
    "DEFAULT_TEMPLATES_DIR",
    "NOT_RECORDED",
    "format_cell",
    "markdown_table",
    "render_template",
]
