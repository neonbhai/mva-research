"""Parse gene-disease clinical validity curations into one row per (gene, disease, source).

Two sources feed this table, kept distinguishable by their own ``source`` value:

* **ClinGen** Gene-Disease Validity curations (``clingen_gene_validity.csv``): a clean
  CSV with a three-line preamble and ``+``-rule separator lines around the header (see
  :data:`_CLINGEN_HEADER_MARKER`). Every classification ClinGen recorded is kept
  verbatim in ``confidence`` -- including ``Disputed``, ``Refuted`` and
  ``No Known Disease Relationship`` -- because a curated-and-refuted claim is real
  evidence, and dropping it would make absence indistinguishable from refutation
  (GP-14, generalised from ``knowledge/adapters/README.md``'s frequency-table
  discussion). ``inheritance`` carries ClinGen's own MOI code unmodified; see
  :data:`CLINGEN_MOI_GLOSSARY` for what each code means (documentation only -- the
  value written to the table is always the source's own code, never the expansion).

* **EBI Gene2Phenotype (DDG2P)** (``DDG2P.csv``): parsed defensively by
  :func:`try_parse_ddg2p`, which content-sniffs the file (gzip magic bytes, or plain
  CSV text) rather than trusting its extension, and returns an explicit "unavailable"
  result -- never a fabricated row -- for anything else. This matters concretely: the
  first attempt to acquire this file (``DDG2P.csv.gz``, from
  ``https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz``) silently saved the
  Gene2Phenotype single-page-app's HTML shell (HTTP 200, ~1KB) instead of the dataset;
  a status-code-only check ("curl -f didn't fail") could not have caught that. The
  working source is the per-panel bulk export,
  ``https://www.ebi.ac.uk/gene2phenotype/api/panel/DD/download/``, which returns plain
  (not gzip) CSV with a different, richer column set than G2P's older bulk-download
  format documented publicly -- ``_DDG2P_REQUIRED_COLUMNS`` reflects what was actually
  found in the real file, not a guess.
"""

from __future__ import annotations

import csv
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Final

GENE_DISEASE_COLUMNS: Final[tuple[str, ...]] = (
    "gene_symbol",
    "hgnc_id",
    "disease_id",
    "disease_name",
    "inheritance",
    "mutation_consequence",
    "confidence",
    "classification_date",
    "gcep",
    "panel",
    "source",
    "version",
)

#: ClinGen's own MOI abbreviations (Gene-Disease Validity SOP), for documentation in
#: the rendered table's header comment. Never used to transform a value -- the table
#: always carries the source's raw code.
CLINGEN_MOI_GLOSSARY: Final[dict[str, str]] = {
    "AD": "autosomal dominant",
    "AR": "autosomal recessive",
    "XL": "X-linked",
    "SD": "semidominant",
    "MT": "mitochondrial",
    "UD": "undetermined",
}

_CLINGEN_HEADER_MARKER: Final = "GENE SYMBOL"
_CLINGEN_FILE_CREATED_PREFIX: Final = "FILE CREATED:"
_CLINGEN_REQUIRED_COLUMNS: Final[tuple[str, ...]] = (
    "GENE SYMBOL",
    "GENE ID (HGNC)",
    "DISEASE LABEL",
    "DISEASE ID (MONDO)",
    "MOI",
    "CLASSIFICATION",
    "CLASSIFICATION DATE",
    "GCEP",
)

#: DDG2P's actual column names (EBI Gene2Phenotype per-panel CSV export, verified
#: against the real downloaded file -- NOT the older bulk-download schema documented
#: in some published G2P papers, which this export does not match). Matched
#: case-insensitively against the header actually present; a mismatch fails closed
#: (see :func:`try_parse_ddg2p`) rather than guessing a positional mapping.
_DDG2P_REQUIRED_COLUMNS: Final[tuple[str, ...]] = (
    "gene symbol",
    "hgnc id",
    "disease name",
    "disease mim",
    "disease mondo",
    "allelic requirement",
    "confidence",
    "variant consequence",
    "date of last review",
    "panel",
)
_GZIP_MAGIC: Final = b"\x1f\x8b"


@dataclass(frozen=True, slots=True)
class AssociationRow:
    """One gene-disease association, from either ClinGen or DDG2P."""

    gene_symbol: str
    hgnc_id: str | None
    disease_id: str | None
    disease_name: str
    inheritance: str | None
    mutation_consequence: str | None
    confidence: str
    classification_date: str | None
    gcep: str | None
    panel: str | None
    source: str
    version: str


@dataclass(frozen=True, slots=True)
class ClinGenResult:
    rows: tuple[AssociationRow, ...]
    version: str
    total_rows: int
    unique_genes: tuple[str, ...]
    moi_counts: dict[str, int]
    classification_counts: dict[str, int]


@dataclass(frozen=True, slots=True)
class DDG2PResult:
    """Either a parsed DDG2P table, or an honest explanation of why there isn't one."""

    rows: tuple[AssociationRow, ...]
    available: bool
    reason: str | None
    total_rows: int = 0
    unique_genes: tuple[str, ...] = ()
    confidence_counts: dict[str, int] | None = None


def _sort_key(row: AssociationRow) -> tuple[str, str, str, str, str]:
    return (
        row.gene_symbol,
        row.disease_id or "",
        row.source,
        row.classification_date or "",
        row.gcep or "",
    )


def parse_clingen(path: Path) -> ClinGenResult:
    if not path.is_file():
        msg = f"ClinGen gene-disease validity file not found: {path}"
        raise FileNotFoundError(msg)
    with path.open(newline="", encoding="utf-8") as handle:
        raw_rows = list(csv.reader(handle))

    created_date: str | None = None
    for raw in raw_rows[:5]:
        if raw and raw[0].startswith(_CLINGEN_FILE_CREATED_PREFIX):
            created_date = raw[0].removeprefix(_CLINGEN_FILE_CREATED_PREFIX).strip()
            break
    if created_date is None:
        msg = (
            f"{path.name}: no {_CLINGEN_FILE_CREATED_PREFIX!r} line found in the first 5 "
            "rows; refusing to invent a version for this table (GP-18)."
        )
        raise ValueError(msg)

    header_index = next(
        (i for i, raw in enumerate(raw_rows) if raw and raw[0] == _CLINGEN_HEADER_MARKER), None
    )
    if header_index is None:
        msg = f"{path.name}: could not find the {_CLINGEN_HEADER_MARKER!r} header row."
        raise ValueError(msg)
    header = raw_rows[header_index]
    idx = {name: position for position, name in enumerate(header)}
    missing_columns = [c for c in _CLINGEN_REQUIRED_COLUMNS if c not in idx]
    if missing_columns:
        msg = f"{path.name} is missing column(s) {missing_columns}."
        raise ValueError(msg)

    def _is_rule_line(raw: list[str]) -> bool:
        return all(cell.strip("+") == "" for cell in raw)

    data_rows = [raw for raw in raw_rows[header_index + 1 :] if raw and not _is_rule_line(raw)]

    rows: list[AssociationRow] = []
    moi_counts: dict[str, int] = {}
    classification_counts: dict[str, int] = {}
    genes: set[str] = set()
    for raw in data_rows:
        gene = raw[idx["GENE SYMBOL"]]
        moi = raw[idx["MOI"]]
        classification = raw[idx["CLASSIFICATION"]]
        genes.add(gene)
        moi_counts[moi] = moi_counts.get(moi, 0) + 1
        classification_counts[classification] = classification_counts.get(classification, 0) + 1
        rows.append(
            AssociationRow(
                gene_symbol=gene,
                hgnc_id=raw[idx["GENE ID (HGNC)"]] or None,
                disease_id=raw[idx["DISEASE ID (MONDO)"]] or None,
                disease_name=raw[idx["DISEASE LABEL"]],
                inheritance=moi or None,
                mutation_consequence=None,  # ClinGen's curation model has no such field.
                confidence=classification,
                classification_date=raw[idx["CLASSIFICATION DATE"]] or None,
                gcep=raw[idx["GCEP"]] or None,
                panel=None,  # DDG2P-specific concept; ClinGen has GCEP instead.
                source="ClinGen",
                version=created_date,
            )
        )

    rows.sort(key=_sort_key)
    return ClinGenResult(
        rows=tuple(rows),
        version=created_date,
        total_rows=len(rows),
        unique_genes=tuple(sorted(genes)),
        moi_counts=moi_counts,
        classification_counts=classification_counts,
    )


def _sniff_and_decode(path: Path) -> str:
    """Return the file's text content, transparently gunzipping if it is gzip.

    Content is sniffed by magic bytes, never trusted from the filename: the DDG2P
    resource has been seen both as ``.csv.gz`` (when the download actually returns
    gzip) and as plain ``.csv`` (the working per-panel export) -- and, worse, as an
    HTML page saved under a ``.csv.gz`` name. Extension-based dispatch would get the
    first two right and the third catastrophically wrong (gzip.open on non-gzip bytes
    raises a low-level BadGzipFile with no context; plain-text-decoding an HTML shell
    "succeeds" and silently produces a garbage table).
    """
    head = path.read_bytes()[:2]
    if head == _GZIP_MAGIC:
        with gzip.open(path, mode="rt", encoding="utf-8", newline="") as handle:
            return handle.read()
    return path.read_text(encoding="utf-8")


def try_parse_ddg2p(path: Path, *, version: str) -> DDG2PResult:
    """Parse DDG2P if, and only if, the file is genuinely CSV (gzip-compressed or
    plain) with the expected columns. Any other shape is reported as unavailable with
    a specific, actionable reason -- never as an empty-but-silent success (GP-14: a
    source that could not be read is not the same claim as a source that was read and
    had nothing).
    """
    if not path.is_file():
        return DDG2PResult(rows=(), available=False, reason=f"{path} does not exist.")

    try:
        text = _sniff_and_decode(path)
    except (OSError, gzip.BadGzipFile, UnicodeDecodeError) as exc:
        return DDG2PResult(
            rows=(), available=False, reason=f"{path.name}: {exc.__class__.__name__}: {exc}"
        )

    try:
        reader = csv.reader(text.splitlines())
        header = next(reader)
        data_rows = list(reader)
    except (csv.Error, StopIteration) as exc:
        return DDG2PResult(
            rows=(), available=False, reason=f"{path.name}: {exc.__class__.__name__}: {exc}"
        )

    if not any(cell.strip() for cell in header):
        return DDG2PResult(
            rows=(),
            available=False,
            reason=(
                f"{path.name} decoded as text but its first line has no recognisable CSV "
                "header (all cells blank); this looks like an HTML page or other non-CSV "
                "content saved under this filename, not the DDG2P export."
            ),
        )

    lowered = {cell.strip().lower(): position for position, cell in enumerate(header)}
    missing_columns = [c for c in _DDG2P_REQUIRED_COLUMNS if c not in lowered]
    if missing_columns:
        return DDG2PResult(
            rows=(),
            available=False,
            reason=(
                f"{path.name} is readable CSV but is missing expected column(s) "
                f"{missing_columns} (found {list(header)!r}); the DDG2P schema may have "
                "changed. Refusing to guess a positional mapping."
            ),
        )

    rows: list[AssociationRow] = []
    genes: set[str] = set()
    confidence_counts: dict[str, int] = {}
    for raw in data_rows:
        if not raw or not raw[lowered["gene symbol"]].strip():
            continue
        gene = raw[lowered["gene symbol"]].strip()
        confidence = raw[lowered["confidence"]].strip()
        hgnc_id = raw[lowered["hgnc id"]].strip()
        mondo = raw[lowered["disease mondo"]].strip()
        disease_mim = raw[lowered["disease mim"]].strip()
        # MONDO preferred (matches ClinGen's own choice of disease vocabulary); an OMIM
        # MIM is the fallback; a small fraction of rows (see the rendered table's
        # header comment for the real count) have neither and are left with no
        # disease_id at all rather than a fabricated one.
        disease_id: str | None
        if mondo:
            disease_id = mondo
        elif disease_mim:
            disease_id = f"OMIM:{disease_mim}"
        else:
            disease_id = None

        genes.add(gene)
        confidence_counts[confidence] = confidence_counts.get(confidence, 0) + 1
        rows.append(
            AssociationRow(
                gene_symbol=gene,
                hgnc_id=f"HGNC:{hgnc_id}" if hgnc_id else None,
                disease_id=disease_id,
                disease_name=raw[lowered["disease name"]].strip(),
                inheritance=raw[lowered["allelic requirement"]].strip() or None,
                mutation_consequence=raw[lowered["variant consequence"]].strip() or None,
                confidence=confidence,
                classification_date=raw[lowered["date of last review"]].strip() or None,
                gcep=None,  # ClinGen-specific concept; DDG2P has "panel" instead.
                panel=raw[lowered["panel"]].strip() or None,
                source="DDG2P",
                version=version,
            )
        )
    rows.sort(key=_sort_key)
    return DDG2PResult(
        rows=tuple(rows),
        available=True,
        reason=None,
        total_rows=len(rows),
        unique_genes=tuple(sorted(genes)),
        confidence_counts=confidence_counts,
    )
