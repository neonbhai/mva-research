"""TSV-backed annotation adapters — the synthetic substitute shipped with this repo.

These adapters read hash-pinned tables under ``knowledge/public/``. They are a
**synthetic stand-in** for VEP/SnpEff (consequences), gnomAD (frequencies) and
ClinVar (clinical assertions): the genes are fictional, the numbers are invented,
and nothing here is biologically valid. Every adapter's ``name`` and ``version``
says so out loud, because those two strings are what end up on the face of every
EvidenceItem and in the run manifest (GP-20).

Why local tables rather than a web service, permanently:

* **PRIV-05.** A remote annotation call is a re-identification vector — the
  proband's exact coordinates would be handed to a third party. This module (and
  the whole ``mva.annotation`` package) is structurally forbidden from importing a
  network client; ``tests/unit/test_architecture.py`` enforces it.
* **GP-30.** A pinned local file gives byte-identical repeat runs. A live API does
  not, and cannot be made to.

Design notes that matter more than they look:

* A variant absent from a table is **omitted from the result mapping**, never
  defaulted. AF = 0 for an unqueried site is a fabricated negative (GP-14).
* Empty cells parse to ``None``, never to ``0`` or ``""``.
* Ordering within a variant is an explicit total order, not file order.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Final

import yaml

from mva.annotation.base import AdapterSet
from mva.determinism import hash_file
from mva.errors import AdapterUnavailableError
from mva.models.base import FrozenModel
from mva.models.variant import (
    ClinicalAssertion,
    ConsequenceAnnotation,
    ImpactSeverity,
    PopulationFrequency,
)

# --------------------------------------------------------------------------- names
# Deliberately unglamorous and self-incriminating: an adapter called
# "local-tsv-consequence" at version "synthetic-v0.0" cannot be mistaken in a report
# for a real VEP run (GP-20).

CONSEQUENCE_ADAPTER_NAME: Final = "local-tsv-consequence"
FREQUENCY_ADAPTER_NAME: Final = "local-tsv-frequency"
NULL_CLINICAL_ADAPTER_NAME: Final = "null-clinical"
NULL_CLINICAL_ADAPTER_VERSION: Final = "unavailable-v0.0"

#: Manifest entry names the default adapter set requires.
CONSEQUENCE_TABLE: Final = "consequences"
FREQUENCY_TABLE: Final = "frequencies"

#: Version stamped on every table by :func:`compute_manifest`. Every table shipped
#: here is fabricated; a real acquisition step must write the real release string.
SYNTHETIC_TABLE_VERSION: Final = "synthetic-v0.0"

#: Fixed "retrieval" date for the fabricated tables. Constant, not a wall-clock read:
#: a regenerated manifest must be byte-identical to the committed one (GP-30).
SYNTHETIC_RETRIEVED_DATE: Final = "2026-08-27"

MANIFEST_VERSION: Final = 1
MANIFEST_GENERATOR: Final = "mva.annotation.local_tables.compute_manifest"
KNOWLEDGE_PUBLIC_DIR: Final = "public"

_CONSEQUENCE_COLUMNS: Final[tuple[str, ...]] = (
    "variant_id",
    "gene_symbol",
    "gene_id",
    "transcript_id",
    "transcript_biotype",
    "is_canonical",
    "is_mane_select",
    "consequence_terms",
    "impact",
    "hgvs_c",
    "hgvs_p",
    "exon",
    "protein_position",
    "amino_acids",
    "splice_ai_delta_max",
    "cadd_phred",
    "revel",
)

_FREQUENCY_COLUMNS: Final[tuple[str, ...]] = (
    "variant_id",
    "source",
    "version",
    "population",
    "allele_frequency",
    "allele_count",
    "allele_number",
    "homozygote_count",
    "filter_status",
)

_TRUE_TOKENS: Final[frozenset[str]] = frozenset({"1", "true", "yes", "y"})
_FALSE_TOKENS: Final[frozenset[str]] = frozenset({"0", "false", "no", "n"})


# --------------------------------------------------------------------------- manifest


class KnowledgeTable(FrozenModel):
    """One hash-pinned local table, as declared in ``knowledge/manifests/knowledge.yaml``.

    ``synthetic`` is a written-down field rather than something inferred from a
    filename: GP-20 requires a fabricated resource to *declare* itself fabricated,
    and a value that must be typed out cannot be forgotten by accident.
    """

    name: str
    path: str
    """Relative to the knowledge root (e.g. ``public/consequences.tsv``)."""

    version: str
    sha256: str
    retrieved: str
    description: str
    synthetic: bool


class KnowledgeManifest(FrozenModel):
    """The versioned index of local knowledge tables.

    Strict (``extra="forbid"``): an unknown key is a loud error rather than a
    silently ignored setting, because a manifest is an integrity claim.
    """

    manifest_version: int
    paths_relative_to: str
    generated_by: str
    tables: dict[str, KnowledgeTable]


def compute_manifest(
    knowledge_root: Path,
    *,
    retrieved: str = SYNTHETIC_RETRIEVED_DATE,
    table_version: str = SYNTHETIC_TABLE_VERSION,
) -> dict[str, Any]:
    """Build the manifest structure for every table under ``<knowledge_root>/public``.

    This is the *only* supported way to produce the committed
    ``knowledge/manifests/knowledge.yaml``: hand-edited hashes drift from the files
    they claim to pin, which turns an integrity check into decoration. Regenerate
    after any table change and commit the diff.

    Deterministic by construction: tables sorted by name, no wall-clock read, and
    the description taken from the table's own first comment line.
    """
    public = knowledge_root / KNOWLEDGE_PUBLIC_DIR
    if not public.is_dir():
        msg = (
            f"Knowledge directory {public.as_posix()!r} does not exist; cannot compute a "
            "knowledge manifest."
        )
        raise AdapterUnavailableError(msg)

    tables: dict[str, Any] = {}
    for path in sorted(public.glob("*.tsv"), key=lambda p: p.name):
        tables[path.stem] = {
            "name": path.stem,
            "path": f"{KNOWLEDGE_PUBLIC_DIR}/{path.name}",
            "version": table_version,
            "sha256": hash_file(path),
            "retrieved": retrieved,
            "description": _first_comment_line(path) or f"Local knowledge table {path.name}.",
            # Every table shipped in this repository is fabricated. A real resource
            # must be added with synthetic: false, a real release version and a URL.
            "synthetic": True,
        }
    return {
        "manifest_version": MANIFEST_VERSION,
        "paths_relative_to": "knowledge",
        "generated_by": MANIFEST_GENERATOR,
        "tables": tables,
    }


#: Prepended to the rendered manifest. The regeneration command lives in the file it
#: generates, so nobody has to guess how the hashes were produced.
MANIFEST_HEADER: Final = """\
# Versioned index of the local knowledge tables this pipeline is allowed to read.
#
# GENERATED FILE - do not hand-edit. Regenerate after ANY table change with:
#
#   uv run python -c "from pathlib import Path; \\
#     from mva.annotation import render_manifest_yaml; \\
#     Path('knowledge/manifests/knowledge.yaml').write_text( \\
#       render_manifest_yaml(Path('knowledge')), encoding='utf-8')"
#
# `path` is relative to the knowledge root (this file lives at
# <knowledge_root>/manifests/knowledge.yaml). Every sha256 is verified before any
# table is read: `load_default_adapters` refuses to annotate against bytes that no
# longer match what the manifest pins, and names the offending file.
#
# `synthetic: true` on every entry is a fact, not a placeholder. Each table here is
# fabricated for the demo case - fictional genes, invented allele frequencies - and
# is NOT biologically valid (GP-20). A real resource is added with synthetic: false,
# its true release version and a retrieval date, via the offline, public-only
# acquisition step described in knowledge/adapters/README.md.
"""


def render_manifest_yaml(
    knowledge_root: Path,
    *,
    retrieved: str = SYNTHETIC_RETRIEVED_DATE,
    table_version: str = SYNTHETIC_TABLE_VERSION,
) -> str:
    """Render the committed manifest text. Deterministic: sorted keys, no wrapping."""
    body = yaml.safe_dump(
        compute_manifest(knowledge_root, retrieved=retrieved, table_version=table_version),
        sort_keys=True,
        default_flow_style=False,
        allow_unicode=False,
        width=10_000,
    )
    return f"{MANIFEST_HEADER}{body}"


def load_manifest(manifest_path: Path) -> KnowledgeManifest:
    """Parse and validate the knowledge manifest."""
    if not manifest_path.is_file():
        msg = (
            f"Knowledge manifest {manifest_path.as_posix()!r} not found. Regenerate it with "
            "mva.annotation.local_tables.compute_manifest."
        )
        raise AdapterUnavailableError(msg)
    try:
        raw: Any = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        msg = f"Knowledge manifest {manifest_path.name} is not valid YAML: {exc.__class__.__name__}"
        raise AdapterUnavailableError(msg) from exc
    if not isinstance(raw, dict):
        msg = f"Knowledge manifest {manifest_path.name} must be a YAML mapping at the top level."
        raise AdapterUnavailableError(msg)
    try:
        return KnowledgeManifest.model_validate(raw)
    except ValueError as exc:
        msg = f"Knowledge manifest {manifest_path.name} does not match the expected schema: {exc}"
        raise AdapterUnavailableError(msg) from exc


def resolve_table_path(knowledge_root: Path, table: KnowledgeTable) -> Path:
    """Resolve a manifest-declared path against the knowledge root.

    Paths are knowledge-root-relative (``public/x.tsv``); a leading ``knowledge/``
    is tolerated so that repo-root-relative spellings — the form used in
    ``config.KnowledgeSources`` — also resolve. Escapes via ``..`` are refused.
    """
    relative = Path(table.path)
    if relative.is_absolute():
        msg = (
            f"Knowledge table {table.name!r} declares an absolute path. Manifest paths must "
            "be relative to the knowledge root so the manifest stays portable."
        )
        raise AdapterUnavailableError(msg)
    parts = relative.parts
    if parts and parts[0] == knowledge_root.name:
        relative = Path(*parts[1:])
    if ".." in relative.parts:
        msg = f"Knowledge table {table.name!r} declares a path escaping the knowledge root."
        raise AdapterUnavailableError(msg)
    return knowledge_root / relative


def verify_manifest(knowledge_root: Path, manifest: KnowledgeManifest) -> None:
    """Check every declared table's sha256 against the bytes on disk.

    Raises :class:`AdapterUnavailableError` naming the offending *file* — never its
    contents (PRIV-09). Verification covers the whole manifest, not just the tables
    this stage happens to read: a partially-verified manifest is not an integrity
    guarantee.
    """
    for name in sorted(manifest.tables):
        table = manifest.tables[name]
        path = resolve_table_path(knowledge_root, table)
        if not path.is_file():
            msg = (
                f"Knowledge table {name!r} declared in the manifest is missing from disk "
                f"(expected {table.path!r} under {knowledge_root.as_posix()!r})."
            )
            raise AdapterUnavailableError(msg)
        actual = hash_file(path)
        if actual != table.sha256:
            msg = (
                f"Knowledge table {name!r} ({table.path!r}) failed its manifest integrity "
                f"check: expected sha256 {table.sha256}, found {actual}. The file changed "
                "without the manifest being regenerated; refusing to annotate against an "
                "unpinned resource. Regenerate with compute_manifest and review the diff."
            )
            raise AdapterUnavailableError(msg)


# --------------------------------------------------------------------------- parsing


def _first_comment_line(path: Path) -> str | None:
    """The table's own one-line self-description, used as the manifest description."""
    for raw in path.read_text(encoding="utf-8").splitlines():
        stripped = raw.strip()
        if stripped.startswith("#"):
            return stripped.lstrip("#").strip() or None
        if stripped:
            return None
    return None


def _read_rows(path: Path, columns: Sequence[str]) -> tuple[tuple[int, dict[str, str]], ...]:
    """Read a ``#``-commented TSV into (line number, row) pairs.

    Split on tabs directly rather than via :mod:`csv`: these tables are unquoted by
    construction, and a dialect sniffer is one more thing that can silently reinterpret
    the data. Row width is checked so a truncated line fails loudly instead of
    shifting every downstream column by one.
    """
    if not path.is_file():
        msg = (
            f"Annotation table {path.as_posix()!r} not found. The local adapters require a "
            "pre-downloaded, hash-pinned table under knowledge/; they never fetch anything."
        )
        raise AdapterUnavailableError(msg)

    header: list[str] | None = None
    rows: list[tuple[int, dict[str, str]]] = []
    for lineno, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        fields = raw.split("\t")
        if header is None:
            header = [field.strip() for field in fields]
            missing = [column for column in columns if column not in header]
            if missing:
                msg = (
                    f"Annotation table {path.name} is missing required column(s) "
                    f"{sorted(missing)}. The header is the schema; see the table's comment block."
                )
                raise AdapterUnavailableError(msg)
            continue
        if len(fields) != len(header):
            msg = (
                f"Annotation table {path.name} line {lineno}: expected {len(header)} "
                f"tab-separated columns, found {len(fields)}."
            )
            raise AdapterUnavailableError(msg)
        rows.append((lineno, dict(zip(header, fields, strict=True))))

    if header is None:
        msg = f"Annotation table {path.name} contains no header row."
        raise AdapterUnavailableError(msg)
    return tuple(rows)


def _cell(row: Mapping[str, str], column: str) -> str | None:
    """A trimmed cell, with the empty string mapped to ``None``.

    The distinction is the whole point: an empty REVEL cell means "not scored", and
    must not become ``0.0``, which reads as "scored, and benign".
    """
    value = row.get(column, "").strip()
    return value or None


def _required_cell(row: Mapping[str, str], column: str, *, path: Path, lineno: int) -> str:
    value = _cell(row, column)
    if value is None:
        msg = f"Annotation table {path.name} line {lineno}: required column {column!r} is empty."
        raise AdapterUnavailableError(msg)
    return value


def _parse_float(row: Mapping[str, str], column: str, *, path: Path, lineno: int) -> float | None:
    value = _cell(row, column)
    if value is None:
        return None
    try:
        return float(value)
    except ValueError as exc:
        msg = f"Annotation table {path.name} line {lineno}: column {column!r} is not a float."
        raise AdapterUnavailableError(msg) from exc


def _parse_int(row: Mapping[str, str], column: str, *, path: Path, lineno: int) -> int | None:
    value = _cell(row, column)
    if value is None:
        return None
    try:
        return int(value)
    except ValueError as exc:
        msg = f"Annotation table {path.name} line {lineno}: column {column!r} is not an integer."
        raise AdapterUnavailableError(msg) from exc


def _parse_bool(row: Mapping[str, str], column: str, *, path: Path, lineno: int) -> bool:
    value = _cell(row, column)
    if value is None:
        return False
    token = value.lower()
    if token in _TRUE_TOKENS:
        return True
    if token in _FALSE_TOKENS:
        return False
    msg = (
        f"Annotation table {path.name} line {lineno}: column {column!r} is not a boolean "
        f"(expected one of {sorted(_TRUE_TOKENS | _FALSE_TOKENS)})."
    )
    raise AdapterUnavailableError(msg)


def _parse_impact(row: Mapping[str, str], *, path: Path, lineno: int) -> ImpactSeverity:
    value = _required_cell(row, "impact", path=path, lineno=lineno).lower()
    try:
        return ImpactSeverity(value)
    except ValueError as exc:
        msg = (
            f"Annotation table {path.name} line {lineno}: column 'impact' is not one of "
            f"{[i.value for i in ImpactSeverity]}."
        )
        raise AdapterUnavailableError(msg) from exc


def _pathogenicity_scores(row: Mapping[str, str], *, path: Path, lineno: int) -> dict[str, float]:
    """Only scores that were actually present. A missing score is absent, not zero."""
    scores: dict[str, float] = {}
    for column, key in (("cadd_phred", "CADD_phred"), ("revel", "REVEL")):
        value = _parse_float(row, column, path=path, lineno=lineno)
        if value is not None:
            scores[key] = value
    return scores


def _build_consequence(
    row: Mapping[str, str], *, path: Path, lineno: int, tool: str, tool_version: str
) -> ConsequenceAnnotation:
    terms = tuple(
        term.strip()
        for term in _required_cell(row, "consequence_terms", path=path, lineno=lineno).split(",")
        if term.strip()
    )
    if not terms:
        msg = f"Annotation table {path.name} line {lineno}: 'consequence_terms' is empty."
        raise AdapterUnavailableError(msg)
    biotype = _cell(row, "transcript_biotype")
    return ConsequenceAnnotation(
        gene_symbol=_required_cell(row, "gene_symbol", path=path, lineno=lineno),
        gene_id=_cell(row, "gene_id"),
        transcript_id=_required_cell(row, "transcript_id", path=path, lineno=lineno),
        transcript_biotype=biotype if biotype is not None else "protein_coding",
        is_canonical=_parse_bool(row, "is_canonical", path=path, lineno=lineno),
        is_mane_select=_parse_bool(row, "is_mane_select", path=path, lineno=lineno),
        consequence_terms=terms,
        impact=_parse_impact(row, path=path, lineno=lineno),
        hgvs_c=_cell(row, "hgvs_c"),
        hgvs_p=_cell(row, "hgvs_p"),
        exon=_cell(row, "exon"),
        intron=_cell(row, "intron"),
        protein_position=_parse_int(row, "protein_position", path=path, lineno=lineno),
        amino_acids=_cell(row, "amino_acids"),
        splice_ai_delta_max=_parse_float(row, "splice_ai_delta_max", path=path, lineno=lineno),
        pathogenicity_scores=_pathogenicity_scores(row, path=path, lineno=lineno),
        source_tool=tool,
        source_tool_version=tool_version,
    )


def _consequence_sort_key(annotation: ConsequenceAnnotation) -> tuple[int, int, str, str]:
    """MANE first, then canonical, then transcript ID.

    A total order, so the annotation list is reproducible regardless of file order —
    and *all* transcripts are kept; ordering is presentation, not selection.
    """
    return (
        0 if annotation.is_mane_select else 1,
        0 if annotation.is_canonical else 1,
        annotation.transcript_id,
        annotation.gene_symbol,
    )


def _build_frequency(row: Mapping[str, str], *, path: Path, lineno: int) -> PopulationFrequency:
    allele_frequency = _parse_float(row, "allele_frequency", path=path, lineno=lineno)
    if allele_frequency is None:
        msg = (
            f"Annotation table {path.name} line {lineno}: 'allele_frequency' is empty. A row "
            "with no frequency is not a frequency of zero; omit the row instead (GP-14)."
        )
        raise AdapterUnavailableError(msg)
    return PopulationFrequency(
        # GP-18: source, version and population always come from the table itself.
        # The adapter never supplies them, because a frequency attributed to the
        # wrong cohort or release is worse than no frequency at all.
        source=_required_cell(row, "source", path=path, lineno=lineno),
        version=_required_cell(row, "version", path=path, lineno=lineno),
        population=_required_cell(row, "population", path=path, lineno=lineno),
        allele_frequency=allele_frequency,
        allele_count=_parse_int(row, "allele_count", path=path, lineno=lineno),
        allele_number=_parse_int(row, "allele_number", path=path, lineno=lineno),
        homozygote_count=_parse_int(row, "homozygote_count", path=path, lineno=lineno),
        filter_status=_cell(row, "filter_status"),
    )


def _frequency_sort_key(frequency: PopulationFrequency) -> tuple[str, str, str]:
    return (frequency.source, frequency.version, frequency.population)


def _unique_ids(variant_ids: Sequence[str]) -> tuple[str, ...]:
    """Deduplicate while preserving caller order."""
    seen: dict[str, None] = {}
    for variant_id in variant_ids:
        seen.setdefault(variant_id, None)
    return tuple(seen)


# --------------------------------------------------------------------------- adapters


class LocalConsequenceAdapter:
    """Consequence annotations from a local TSV. A synthetic stand-in for VEP.

    The whole table is parsed once at construction: it is small, and eager parsing
    means a malformed row fails at wiring time rather than halfway through a run.
    """

    def __init__(self, table_path: Path, *, version: str) -> None:
        self._table_path = table_path
        self._version = version
        self._index: dict[str, tuple[ConsequenceAnnotation, ...]] = _index_consequences(
            table_path, tool=CONSEQUENCE_ADAPTER_NAME, tool_version=version
        )

    @property
    def name(self) -> str:
        return CONSEQUENCE_ADAPTER_NAME

    @property
    def version(self) -> str:
        return self._version

    @property
    def synthetic(self) -> bool:
        """Always True. The table is fabricated; see GP-20 and the maturity ledger."""
        return True

    @property
    def table_path(self) -> Path:
        return self._table_path

    def annotate(
        self, variant_ids: Sequence[str]
    ) -> Mapping[str, tuple[ConsequenceAnnotation, ...]]:
        """Return annotations for the variants this table knows about.

        Unknown variants are **omitted**, not returned with an empty tuple: this
        adapter cannot distinguish "intergenic" from "not in my table", and it is
        not entitled to imply the former.
        """
        return {
            variant_id: self._index[variant_id]
            for variant_id in _unique_ids(variant_ids)
            if variant_id in self._index
        }


class LocalFrequencyAdapter:
    """Population frequencies from a local TSV. A synthetic stand-in for gnomAD."""

    def __init__(self, table_path: Path, *, version: str) -> None:
        self._table_path = table_path
        self._version = version
        self._index: dict[str, tuple[PopulationFrequency, ...]] = _index_frequencies(table_path)

    @property
    def name(self) -> str:
        return FREQUENCY_ADAPTER_NAME

    @property
    def version(self) -> str:
        return self._version

    @property
    def synthetic(self) -> bool:
        """Always True. These are invented allele frequencies, not gnomAD."""
        return True

    @property
    def table_path(self) -> Path:
        return self._table_path

    def frequencies(
        self, variant_ids: Sequence[str]
    ) -> Mapping[str, tuple[PopulationFrequency, ...]]:
        """Return frequency observations for variants present in the table.

        A variant absent from the table is absent from this mapping. It has no
        frequency data — which is emphatically not an allele frequency of zero, and
        not evidence of rarity (GP-14). The caller must handle the missing key.
        """
        return {
            variant_id: self._index[variant_id]
            for variant_id in _unique_ids(variant_ids)
            if variant_id in self._index
        }


class NullClinicalAdapter:
    """A clinical adapter with nothing to say, saying so honestly.

    This repository ships no ClinVar substitute: inventing clinical significance
    for fictional variants would be the single most misleading thing it could do.
    Rather than omit the slot, the slot is filled by an adapter whose name and
    version state that no source is available, so the absence is visible in the
    coverage table and in the run's warnings instead of being invisible.

    It is *not* evidence of benignity, and the service says so.
    """

    @property
    def name(self) -> str:
        return NULL_CLINICAL_ADAPTER_NAME

    @property
    def version(self) -> str:
        return NULL_CLINICAL_ADAPTER_VERSION

    @property
    def synthetic(self) -> bool:
        """True: this is a stand-in for a real ClinVar adapter, not one."""
        return True

    def assertions(self, variant_ids: Sequence[str]) -> Mapping[str, tuple[ClinicalAssertion, ...]]:
        """Always empty. No source is configured, so nothing is on record."""
        _ = variant_ids
        return {}


def _index_consequences(
    table_path: Path, *, tool: str, tool_version: str
) -> dict[str, tuple[ConsequenceAnnotation, ...]]:
    grouped: dict[str, list[ConsequenceAnnotation]] = {}
    for lineno, row in _read_rows(table_path, _CONSEQUENCE_COLUMNS):
        variant_id = _required_cell(row, "variant_id", path=table_path, lineno=lineno)
        annotation = _build_consequence(
            row, path=table_path, lineno=lineno, tool=tool, tool_version=tool_version
        )
        grouped.setdefault(variant_id, []).append(annotation)
    return {
        variant_id: tuple(sorted(annotations, key=_consequence_sort_key))
        for variant_id, annotations in grouped.items()
    }


def _index_frequencies(table_path: Path) -> dict[str, tuple[PopulationFrequency, ...]]:
    grouped: dict[str, list[PopulationFrequency]] = {}
    for lineno, row in _read_rows(table_path, _FREQUENCY_COLUMNS):
        variant_id = _required_cell(row, "variant_id", path=table_path, lineno=lineno)
        grouped.setdefault(variant_id, []).append(
            _build_frequency(row, path=table_path, lineno=lineno)
        )
    return {
        variant_id: tuple(sorted(frequencies, key=_frequency_sort_key))
        for variant_id, frequencies in grouped.items()
    }


def load_default_adapters(knowledge_root: Path, manifest_path: Path) -> AdapterSet:
    """Wire the default (synthetic) adapter set from hash-pinned local tables.

    Order of operations is deliberate: the manifest is parsed, **every** declared
    table's sha256 is verified, and only then is anything read as data. An adapter
    built over an unverified file would be annotating against unknown bytes.

    Raises :class:`AdapterUnavailableError` if the manifest is missing/invalid, a
    required entry is absent, or any table's hash disagrees with the manifest.
    """
    manifest = load_manifest(manifest_path)
    verify_manifest(knowledge_root, manifest)

    consequences = _require_table(manifest, CONSEQUENCE_TABLE, manifest_path=manifest_path)
    frequencies = _require_table(manifest, FREQUENCY_TABLE, manifest_path=manifest_path)

    return AdapterSet(
        consequence=LocalConsequenceAdapter(
            resolve_table_path(knowledge_root, consequences), version=consequences.version
        ),
        frequency=LocalFrequencyAdapter(
            resolve_table_path(knowledge_root, frequencies), version=frequencies.version
        ),
        clinical=NullClinicalAdapter(),
    )


def _require_table(
    manifest: KnowledgeManifest, name: str, *, manifest_path: Path
) -> KnowledgeTable:
    table = manifest.tables.get(name)
    if table is None:
        msg = (
            f"Knowledge manifest {manifest_path.name} declares no table {name!r}; the default "
            f"annotation adapters require {CONSEQUENCE_TABLE!r} and {FREQUENCY_TABLE!r}."
        )
        raise AdapterUnavailableError(msg)
    return table
