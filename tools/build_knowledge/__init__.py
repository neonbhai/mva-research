"""Deterministic, offline generators that turn downloaded public resources into the
hash-pinnable local knowledge tables under ``knowledge/real/``.

This package reads only local files (never the network) under a resource root supplied
by the caller (CLI flag or the ``MVA_RESOURCES`` environment variable — see
``tools/build_knowledge/build.py``), and writes plain, ``#``-commented TSVs with real,
verifiable provenance in their header comments: source name, true release/version
string (read from the data or the download, never invented), retrieval date and
licence.

Three rules every module here follows, mirroring ``knowledge/adapters/README.md``:

1. **Never invent or impute (GP-14).** A missing source value is an empty TSV cell,
   never a fabricated ``0``. ``mva.annotation.local_tables`` already treats an empty
   cell as ``None``; nothing here may compromise that by defaulting.
2. **Preserve the source's own vocabulary.** Confidence/classification strings
   (ClinGen's ``Definitive``/``Disputed``/``Refuted``, HPO's own frequency terms) are
   carried through verbatim, never remapped onto an invented scale.
3. **Deterministic output (GP-30).** Every generator sorts by an explicit key before
   writing; see ``tools/build_knowledge/tsv_io.py``.

Not part of the installed ``mva`` package (see ``pyproject.toml``'s
``[tool.hatch.build.targets.wheel]``) and imports nothing from ``src/mva``: this is a
standalone, one-shot acquisition-to-table tool, run with::

    uv run python -m tools.build_knowledge.build
"""

from __future__ import annotations
