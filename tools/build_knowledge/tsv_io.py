"""Deterministic TSV rendering shared by every generator in this package.

GP-30: repeat runs must be byte-identical. Callers pass rows that are already sorted by
an explicit key and already formatted to strings (via :func:`format_float` /
:func:`format_cell`, so the missing-vs-zero decision stays visible at the call site,
not buried in a shared formatter); this module only concatenates and writes.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Final

#: The acquisition date recorded in every real table's header comment. A literal
#: constant, not a wall-clock read (mirrors ``mva.annotation.local_tables``'s own
#: ``SYNTHETIC_RETRIEVED_DATE`` pattern): a regenerated table must render byte-identical
#: text to the committed one, which a ``datetime.now()`` call would break.
RETRIEVED_DATE: Final = "2026-08-28"

#: Six significant figures, fixed notation choice (Python's ``g`` format switches to
#: scientific notation only when that is shorter/more precise than fixed) -- applied
#: identically on every run of every table, so two runs over the same bytes format the
#: same float the same way.
_FLOAT_FORMAT: Final = "{:.6g}"


def format_float(value: float | None) -> str:
    """Render a numeric cell. ``None`` (unmeasured/unknown) becomes an EMPTY string.

    Never ``"0"`` for a missing value (GP-14): a gene with no pLI in the source must
    stay silent, not read as "pLI = 0" (which would mean something specific and false).
    """
    if value is None:
        return ""
    return _FLOAT_FORMAT.format(value)


def format_cell(value: str | None) -> str:
    """Render a text cell. ``None`` becomes an EMPTY string, never a placeholder."""
    if value is None:
        return ""
    return value


def render_comment_block(lines: Sequence[str]) -> tuple[str, ...]:
    """Prefix each line with ``# ``; a blank line becomes a bare ``#`` (no trailing space)."""
    return tuple(f"# {line}" if line else "#" for line in lines)


def _check_cell_is_safe(value: str, *, column: str, row_index: int) -> None:
    """These TSVs are unquoted by construction (see ``local_tables.py``'s ``_read_rows``
    docstring): a tab or newline inside a cell would silently corrupt column alignment
    for every row after it, so it is refused here instead.
    """
    if "\t" in value or "\n" in value or "\r" in value:
        msg = (
            f"Value for column {column!r} in row {row_index} contains a tab or newline "
            f"character ({value!r}); these TSVs are unquoted and cannot carry one."
        )
        raise ValueError(msg)


def render_tsv(
    *,
    comment_lines: Sequence[str],
    columns: Sequence[str],
    rows: Sequence[Mapping[str, str]],
) -> str:
    """Render full TSV text: ``#``-prefixed provenance comments, header, then data rows.

    ``rows`` must already be in final, deterministic order (GP-30) and fully formatted
    to strings. Every row must supply every column; a missing key is a bug in the
    caller, not something to silently pad.
    """
    lines: list[str] = list(render_comment_block(comment_lines))
    lines.append("\t".join(columns))
    for row_index, row in enumerate(rows):
        cells: list[str] = []
        for column in columns:
            if column not in row:
                msg = f"Row {row_index} is missing column {column!r}."
                raise ValueError(msg)
            value = row[column]
            _check_cell_is_safe(value, column=column, row_index=row_index)
            cells.append(value)
        lines.append("\t".join(cells))
    return "\n".join(lines) + "\n"


def write_tsv(
    path: Path,
    *,
    comment_lines: Sequence[str],
    columns: Sequence[str],
    rows: Sequence[Mapping[str, str]],
) -> None:
    """Render and write a TSV with explicit ``\\n`` line endings (GP-30: no CRLF drift)."""
    content = render_tsv(comment_lines=comment_lines, columns=columns, rows=rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(content)
