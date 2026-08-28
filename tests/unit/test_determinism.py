"""Byte-identity guarantees for the foundation serialisation helpers (GP-30).

Repeat runs being byte-identical is an acceptance criterion for this project, and
every artifact hash in the run manifest rests on `canonical_json` producing the
same bytes for the same value. These tests hold that floor.
"""

from __future__ import annotations

import json
from collections.abc import Iterator
from pathlib import Path

import pytest

from mva.determinism import canonical_json, write_canonical_json_rows


@pytest.mark.unit
def test_streamed_rows_are_byte_identical_to_canonical_json(tmp_path: Path) -> None:
    """The streaming writer must not become a second, subtly different encoder.

    Its whole justification is avoiding a 2.6 GB in-memory string while producing
    the same bytes. If the two ever diverge, every artifact hash and the GP-30
    byte-identity claim quietly stop meaning what they say -- so the equivalence
    is asserted here rather than reasoned about in a docstring.
    """
    rows: list[object] = [
        {"b": 2, "a": 1},
        {"z": [3, 1, 2], "y": None},
        {"nested": {"d": 4, "c": 3}},
        {"unicode": "é", "float": 1.5, "bool": True},
    ]
    path = tmp_path / "rows.json"
    written = write_canonical_json_rows(path, iter(rows))

    assert written == len(rows)
    assert path.read_text(encoding="utf-8") == canonical_json(rows) + "\n"


@pytest.mark.unit
def test_streamed_rows_handle_the_empty_case(tmp_path: Path) -> None:
    """An empty callset must still produce a valid, parseable array."""
    path = tmp_path / "empty.json"

    assert write_canonical_json_rows(path, iter(())) == 0
    assert path.read_text(encoding="utf-8") == canonical_json([]) + "\n"
    assert json.loads(path.read_text(encoding="utf-8")) == []


@pytest.mark.unit
def test_streamed_rows_consume_the_iterable_exactly_once(tmp_path: Path) -> None:
    """A generator must not be walked twice; at WGS scale a rewind is not possible."""
    consumed: list[int] = []

    def rows() -> Iterator[dict[str, int]]:
        for i in range(3):
            consumed.append(i)
            yield {"i": i}

    write_canonical_json_rows(tmp_path / "once.json", rows())
    assert consumed == [0, 1, 2]
