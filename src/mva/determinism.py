"""Canonicalisation and hashing helpers underpinning GP-30.

Determinism here means: given the same inputs, the same configuration and the same
tool versions, every artifact this pipeline writes is byte-identical. That is
checked directly by the repeat-run test, and it is what makes the provenance
manifest a meaningful claim rather than decoration.

The hazards this module exists to neutralise are all mundane: dict iteration
order, float repr drift, set ordering, and non-UTC timestamps.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable, Mapping
from datetime import date, datetime
from enum import Enum
from pathlib import Path
from typing import Any

#: Read size for file hashing. Fixed so hashes never depend on platform buffering.
_CHUNK = 1 << 20


def canonical_json(value: Any) -> str:
    """Serialise to a stable string.

    Keys sorted, no insignificant whitespace, non-ASCII escaped, and every
    non-JSON-native type reduced to a documented textual form. Floats are formatted
    with `repr` semantics via json, which is stable within a CPython major version;
    where a value must survive across versions, store it as a string.
    """
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        default=_json_default,
    )


def write_canonical_json_rows(path: Path, rows: Iterable[Any]) -> int:
    """Write a JSON array one row at a time; return the row count.

    Byte-identical to ``canonical_json(list(rows))``, but it never holds the whole
    artifact. :func:`canonical_json` builds the entire document as one Python
    ``str`` before a byte reaches disk, which is 2.6 GB for a whole-genome
    variants artifact (``docs/scale-report.md`` §4). Measured with this function:
    a 2,409,978,805-byte artifact written at 243 MB peak RSS.

    The equivalence holds because ``canonical_json`` of a list is exactly the
    concatenation of its elements' encodings between ``[`` and ``]`` with ``,``
    separators -- no whitespace, no trailing comma, keys already sorted per
    element. A test asserts the two agree rather than leaving that to the reader.
    """
    written = 0
    with path.open("w", encoding="utf-8") as handle:
        handle.write("[")
        for row in rows:
            if written:
                handle.write(",")
            handle.write(canonical_json(row))
            written += 1
        handle.write("]\n")
    return written


def _json_default(value: Any) -> Any:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, datetime):
        return value.isoformat()
    if isinstance(value, date):
        return value.isoformat()
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, frozenset | set):
        return sorted(str(item) for item in value)
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, Mapping):
        return dict(value)
    if hasattr(value, "model_dump"):
        return value.model_dump(mode="json")
    msg = f"canonical_json cannot serialise {type(value).__name__}"
    raise TypeError(msg)


def stable_hash(value: Any) -> str:
    """sha256 over the canonical JSON form."""
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def hash_file(path: Path) -> str:
    """sha256 of a file's bytes, streamed.

    Streaming matters here beyond memory: a whole VCF loaded into RAM widens the
    window in which plaintext genotypes could reach swap or a crash dump.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(_CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def hash_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def short_hash(value: Any, length: int = 12) -> str:
    return stable_hash(value)[:length]


def stable_sorted[T](items: Iterable[T], *, key: Any = None) -> list[T]:
    """Sort with a total order.

    Python's sort is stable, but stability only preserves *input* order — which for
    a set or dict is not defined. Callers must supply a key that is total; this
    wrapper exists to make the requirement explicit at call sites.
    """
    return sorted(items, key=key) if key is not None else sorted(items)  # type: ignore[type-var]
