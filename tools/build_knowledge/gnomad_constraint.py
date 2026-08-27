"""Parse gnomAD gene constraint metrics into one row per gene.

Source shape: ``gnomad.vX.Y.constraint_metrics.tsv`` has one row per **transcript**
(a gene can have several -- an Ensembl-annotated transcript and, separately, a
RefSeq/Entrez-annotated representation of the same biology), so this module selects a
single representative transcript per gene deterministically: MANE Select first, then
canonical, then the lexicographically smallest transcript ID -- the same three-level
tie-break ``mva.annotation.local_tables._consequence_sort_key`` uses for consequence
annotations, applied here for consistency across the repo's knowledge tables.

Only rows whose ``gene_id`` starts with ``ENSG`` are considered: gnomAD v4.1 also
carries a RefSeq/Entrez-ID mirror of many transcripts (``gene_id`` holding a bare
Entrez integer, e.g. ``"1"`` for *A1BG*), and mixing two ID systems in one column
without a discriminator is exactly the join hazard this project's rules warn about. A
gene with no ``ENSG``-tagged row at all is *excluded* and reported, never silently
folded in with a blank ID.

Metric choice (verified against gnomAD's own public documentation and browser
behaviour, not guessed): ``pli``/``loeuf`` come from the ``lof.`` column group, not
``lof_hc_lc.`` -- ``lof.`` is gnomAD's standard, browser-displayed metric computed from
high-confidence LOFTEE pLoF variants only; ``lof_hc_lc.`` merges in low-confidence
calls and is not what "pLI"/"LOEUF" ordinarily refers to. ``mis_z`` comes from
``mis.z_score`` (the sample-size-calibrated score), not ``mis.z_raw``: ``z_score`` is
the direct v4 successor of the single ``mis_z`` column gnomAD v2 shipped.

A source cell of ``"NA"`` becomes ``None``, never ``0.0`` (GP-14): roughly a fifth of
genes have no computable LOEUF/pLI (typically because they are too short for a
meaningful expected-variant count -- flagged ``no_exp_lof`` in ``constraint_flags``),
and reporting that as ``0`` would silently claim "maximally constrained" for a gene
that was never actually measured.
"""

from __future__ import annotations

import csv
import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Final

GENE_PANEL_COLUMNS: Final[tuple[str, ...]] = (
    "gene_symbol",
    "ensembl_gene_id",
    "transcript_id",
    "canonical",
    "mane_select",
    "pli",
    "loeuf",
    "mis_z",
    "constraint_flags",
    "source",
    "version",
)

#: gnomAD's own missing-value sentinel in this file.
_NA: Final = "NA"

_REQUIRED_SOURCE_COLUMNS: Final[tuple[str, ...]] = (
    "gene",
    "gene_id",
    "transcript",
    "canonical",
    "mane_select",
    "lof.pLI",
    "lof.oe_ci.upper",
    "mis.z_score",
    "constraint_flags",
)

_FILENAME_VERSION_RE: Final[re.Pattern[str]] = re.compile(
    r"^gnomad\.v(?P<version>[0-9]+(?:\.[0-9]+)*)\.constraint_metrics\.tsv$"
)


@dataclass(frozen=True, slots=True)
class GenePanelRow:
    """One gene's selected representative-transcript constraint metrics."""

    gene_symbol: str
    ensembl_gene_id: str
    transcript_id: str
    canonical: bool
    mane_select: bool
    pli: float | None
    loeuf: float | None
    mis_z: float | None
    constraint_flags: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class GnomadConstraintResult:
    """Parsed gene-level constraint table plus the honest accounting of what didn't join."""

    rows: tuple[GenePanelRow, ...]
    version: str
    total_transcript_rows: int
    na_gene_rows: int
    genes_seen: int
    genes_without_ensembl_row: tuple[str, ...]

    @property
    def genes_included(self) -> int:
        return len(self.rows)


def extract_gnomad_version(path: Path) -> str:
    """The gnomAD release, read from the filename (the download names its own release;
    the TSV body carries no version column). Refuses to guess if the name is unexpected.
    """
    match = _FILENAME_VERSION_RE.match(path.name)
    if match is None:
        msg = (
            f"Cannot determine the gnomAD release from filename {path.name!r}; expected "
            "the form 'gnomad.vX.Y.constraint_metrics.tsv'. Refusing to guess a version "
            "(GP-18: every table needs a real, traceable version)."
        )
        raise ValueError(msg)
    return f"v{match.group('version')}"


_TRUE_LITERAL: Final = "true"
_FALSE_LITERAL: Final = "false"


def _parse_float(raw_value: str) -> float | None:
    if raw_value in (_NA, ""):
        return None
    return float(raw_value)


def _parse_bool(raw_value: str, *, context: str) -> bool:
    if raw_value == _TRUE_LITERAL:
        return True
    if raw_value == _FALSE_LITERAL:
        return False
    msg = f"{context}: expected {_TRUE_LITERAL!r} or {_FALSE_LITERAL!r}, found {raw_value!r}."
    raise ValueError(msg)


def _parse_constraint_flags(raw_value: str) -> tuple[str, ...]:
    """Turn ``'["no_exp_lof","outlier_mis"]'`` into ``("no_exp_lof", "outlier_mis")``;
    ``'[]'`` becomes ``()``.
    """
    stripped = raw_value.strip("[]")
    if not stripped:
        return ()
    return tuple(part.strip().strip('"').strip("'") for part in stripped.split(","))


def _select(
    best_by_gene: dict[str, tuple[tuple[int, int, str], GenePanelRow]],
    row: GenePanelRow,
) -> None:
    # MANE Select first, then canonical, then transcript ID: a total, explicit order
    # so the choice is the same regardless of file/iteration order (GP-30).
    sort_key = (0 if row.mane_select else 1, 0 if row.canonical else 1, row.transcript_id)
    existing = best_by_gene.get(row.gene_symbol)
    if existing is None or sort_key < existing[0]:
        best_by_gene[row.gene_symbol] = (sort_key, row)


def parse_gnomad_constraint(path: Path) -> GnomadConstraintResult:
    if not path.is_file():
        msg = f"gnomAD constraint file not found: {path}"
        raise FileNotFoundError(msg)
    version = extract_gnomad_version(path)

    best_by_gene: dict[str, tuple[tuple[int, int, str], GenePanelRow]] = {}
    total_rows = 0
    na_gene_rows = 0
    genes_seen: set[str] = set()
    genes_with_ensembl: set[str] = set()

    with path.open(newline="", encoding="utf-8") as handle:
        reader: Iterable[list[str]] = csv.reader(handle, delimiter="\t")
        header = next(reader)
        idx = {name: position for position, name in enumerate(header)}
        missing_columns = [c for c in _REQUIRED_SOURCE_COLUMNS if c not in idx]
        if missing_columns:
            msg = f"{path.name} is missing column(s) {missing_columns}."
            raise ValueError(msg)

        for fields in reader:
            total_rows += 1
            gene = fields[idx["gene"]]
            if gene == _NA:
                na_gene_rows += 1
                continue
            genes_seen.add(gene)
            gene_id = fields[idx["gene_id"]]
            if not gene_id.startswith("ENSG"):
                continue
            genes_with_ensembl.add(gene)
            transcript_id = fields[idx["transcript"]]
            row = GenePanelRow(
                gene_symbol=gene,
                ensembl_gene_id=gene_id,
                transcript_id=transcript_id,
                canonical=_parse_bool(fields[idx["canonical"]], context=f"{path.name} gene {gene}"),
                mane_select=_parse_bool(
                    fields[idx["mane_select"]], context=f"{path.name} gene {gene}"
                ),
                pli=_parse_float(fields[idx["lof.pLI"]]),
                loeuf=_parse_float(fields[idx["lof.oe_ci.upper"]]),
                mis_z=_parse_float(fields[idx["mis.z_score"]]),
                constraint_flags=_parse_constraint_flags(fields[idx["constraint_flags"]]),
            )
            _select(best_by_gene, row)

    genes_without_ensembl = tuple(sorted(genes_seen - genes_with_ensembl))
    rows = tuple(sorted((entry[1] for entry in best_by_gene.values()), key=lambda r: r.gene_symbol))
    return GnomadConstraintResult(
        rows=rows,
        version=version,
        total_transcript_rows=total_rows,
        na_gene_rows=na_gene_rows,
        genes_seen=len(genes_seen),
        genes_without_ensembl_row=genes_without_ensembl,
    )
