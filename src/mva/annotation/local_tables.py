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
* **Both sides of the join go through :func:`mva.alleles.canonicalise_allele`**
  (ADR 0018). See :func:`_join_key`.

Representation, and what these adapters may not claim
-----------------------------------------------------

``chr1:100 AT>AG`` and ``chr1:101 T>G`` are one substitution written two ways.
These adapters used to store the table's ``variant_id`` column verbatim and look
it up by direct string membership, so the two spellings did not join — and a
failed join does not raise. It returns *absence*, which this pipeline reads as
"no frequency data" (GP-14) and, for a missing consequence, as grounds for
dropping the allele in selection. A failed join is indistinguishable from "novel
and ultra-rare", which is the profile of a causal variant. ADR 0018 established
one canonicalisation rule at layer 1 and required every join key to go through
it; the real gnomAD and ClinVar adapters were converted then and these were not,
which mattered because these are the **default executable adapter path**.

A TSV has no reference genome, so what is available here is **trimming, not
left-alignment**. Trimming needs no reference and is always correct; rolling an
indel leftwards through a repeat tract requires reading the bases to its left.
The honest consequence is that an indel spelled at a different legal position
within a repeat still will not join, and neither adapter is entitled to imply
otherwise — so each publishes a
:class:`~mva.alleles.LeftAlignmentReport` on its ``left_alignment`` property
stating exactly that (``UNAVAILABLE_NO_REFERENCE`` whenever the table holds an
indel). Trim-only is strictly better than raw string comparison; it is not
complete, and the difference is declared rather than hoped over.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Final

import yaml

from mva.alleles import (
    LeftAlignmentReport,
    canonicalise_allele,
    is_sequence_allele,
    summarise_left_alignment,
)
from mva.annotation.base import AdapterSet
from mva.determinism import hash_file
from mva.errors import AdapterUnavailableError
from mva.models.base import FrozenModel
from mva.models.genome import ContigStyle, GenomeBuild, normalise_contig
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
    if not rows:
        # A header-only table used to build cleanly. The manifest hash verified, the
        # schema check passed, the AdapterSet was constructed, and every lookup then
        # missed — silently, because a miss in these adapters is the *correct*
        # answer for a variant the table does not hold (GP-14). Every variant lost
        # its `gene_symbols`, and compound-heterozygous pairing groups by gene, so
        # the whole run produced zero candidate pairs while reporting success.
        # Refused here for the same reason `gene_intervals._read_gtf` refuses an
        # empty gene model and `GnomadSitesFrequencyAdapter` refuses a release with
        # no complete shard: an index over nothing answers every question with
        # absence, and absence is what this pipeline is least able to detect.
        msg = (
            f"Annotation table {path.name} has a header but no data rows. An adapter over "
            "an empty index answers every lookup with 'not in my table', which is "
            "indistinguishable from a variant genuinely having no annotation and would "
            "silently strip the gene assignment off every variant in the run — collapsing "
            "compound-heterozygous pairing to zero pairs while the run reported success. "
            "Refusing to build it. Re-run the knowledge-table build, or check that the "
            "file was not truncated to its header."
        )
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


# ------------------------------------------------------------------- the join key


def _join_key(variant_id: str) -> str:
    """The one representation both sides of the join are reduced to (ADR 0018).

    Applied to the table's ``variant_id`` column when the index is built **and** to
    the caller's ID when it is looked up. Both sides therefore agree by
    construction rather than by whoever happened to spell it first — which is the
    whole content of ADR 0018, and the reason this function delegates instead of
    trimming anything itself. There is deliberately no allele arithmetic in this
    module; a second implementation of the rule is the defect, not the cure, and
    ``tests/unit/test_normalise_representation.py`` fails the build over any that
    appears here.

    Three normalisations, in order, each one a comparison that would otherwise fail
    on a spelling difference:

    * **Build** through :meth:`~mva.models.genome.GenomeBuild.parse`, so ``hg38``
      and ``GRCh38`` name one assembly — while ``GRCh37`` still never joins
      ``GRCh38`` (GP-11: the same locus differs by megabases).
    * **Contig** through :func:`~mva.models.genome.normalise_contig`, so ``15`` and
      ``chr15`` are one chromosome. (Only the *challenge scorer* compares contig
      strings raw; that is a property of the submission CSV, not of this join.)
    * **Alleles** through :func:`mva.alleles.canonicalise_allele`.

    **No reference is passed, and none is available.** These adapters read a TSV;
    there is no FASTA to roll an indel leftwards against. The result is trimmed and
    honestly not left-aligned — ``operations`` will not claim otherwise, and
    :func:`_left_alignment_report` turns that into a statement the caller receives.

    An ID that is not ``build:contig:pos:ref:alt``, names an unrecognised assembly,
    or carries a non-canonical contig is returned **verbatim**. Refusing to guess is
    the point: a key reinterpreted on a guess joins the wrong record, which is
    worse than joining none.
    """
    parts = variant_id.split(":")
    if len(parts) != 5:
        return variant_id
    build_token, contig_token, position_token, ref, alt = parts
    try:
        build = GenomeBuild.parse(build_token)
        contig = normalise_contig(contig_token, ContigStyle.UCSC)
        position = int(position_token)
    except ValueError:
        return variant_id
    canonical = canonicalise_allele(
        contig=contig,
        position=position,
        ref=ref.strip().upper(),
        alt=alt.strip().upper(),
    )
    return f"{build.value}:{contig}:{canonical.position}:{canonical.ref}:{canonical.alt}"


def _key_is_indel(key: str) -> bool:
    """True for a length-changing key. Symbolic alleles are not indels."""
    parts = key.split(":")
    if len(parts) != 5:
        return False
    ref, alt = parts[3], parts[4]
    if not is_sequence_allele(ref) or not is_sequence_allele(alt):
        return False
    return len(ref) != len(alt)


def _left_alignment_report(indel_queries: int) -> LeftAlignmentReport:
    """What this adapter's keys may and may not be trusted to join (GP-14).

    Derived from the counts rather than asserted, and derived through
    :func:`~mva.alleles.summarise_left_alignment` so no adapter can label its own
    batch by hand. ``reference_available=False`` is not a placeholder: a TSV
    adapter has no reference and never will, so every indel it is asked about is
    trimmed-only and the status is ``UNAVAILABLE_NO_REFERENCE``. A run that asks
    only about SNVs reports ``NOT_REQUIRED`` instead — "could not left-align" and
    "had nothing to left-align" are opposite claims about how far to trust the
    rarity of every indel in the run, and they must not share a value.

    ``indel_queries`` counts the **variants this run looked up**, not the rows in
    the table. This function used to take the index and count
    ``sum(1 for key in index if _key_is_indel(key))``, which is the wrong side of
    the join in the direction that matters: an SNV-only table asked about a run
    full of indels reported ``NOT_REQUIRED`` and asserted "this run contains no
    indel records" — over precisely the indels that were, at that moment, silently
    returning absence because a trim-only key cannot match a left-aligned one. The
    table's own composition is not a fact about the run; the query set is.
    """
    return summarise_left_alignment(
        indel_count=indel_queries,
        shifted_count=0,
        unaligned_indel_count=indel_queries,
        reference_available=False,
    )


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
        #: Indel variants this adapter has been ASKED about, accumulated across
        #: calls. Not a property of the table: see :func:`_left_alignment_report`.
        self._indel_queries = 0

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

    @property
    def left_alignment(self) -> LeftAlignmentReport:
        """Whether this run's indel keys can be trusted to join (GP-14).

        Always ``reference_available=False``: a TSV adapter has no FASTA, so its
        keys are trimmed and not left-aligned. Counted over the variants this
        adapter has been asked about rather than over the table's own rows, so it
        accumulates as the run proceeds. See :func:`_left_alignment_report`.
        """
        return _left_alignment_report(self._indel_queries)

    def annotate(
        self, variant_ids: Sequence[str]
    ) -> Mapping[str, tuple[ConsequenceAnnotation, ...]]:
        """Return annotations for the variants this table knows about.

        Unknown variants are **omitted**, not returned with an empty tuple: this
        adapter cannot distinguish "intergenic" from "not in my table", and it is
        not entitled to imply the former.

        The lookup is by :func:`_join_key`, but the result is keyed by the
        **caller's own ID string**. ``mva.annotation.service`` resolves the mapping
        with ``record.variant_id``; re-keying to the canonical form would move the
        silent miss one layer up instead of removing it.
        """
        found: dict[str, tuple[ConsequenceAnnotation, ...]] = {}
        for variant_id in _unique_ids(variant_ids):
            key = _join_key(variant_id)
            # Counted on every query, hit or miss. A miss is exactly the case the
            # left-alignment report exists to qualify: it is the shape a
            # representation mismatch takes.
            if _key_is_indel(key):
                self._indel_queries += 1
            entry = self._index.get(key)
            if entry is not None:
                found[variant_id] = entry
        return found


class LocalFrequencyAdapter:
    """Population frequencies from a local TSV. A synthetic stand-in for gnomAD."""

    def __init__(self, table_path: Path, *, version: str) -> None:
        self._table_path = table_path
        self._version = version
        self._index: dict[str, tuple[PopulationFrequency, ...]] = _index_frequencies(table_path)
        #: Indel variants this adapter has been ASKED about, accumulated across
        #: calls. Not a property of the table: see :func:`_left_alignment_report`.
        self._indel_queries = 0

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

    @property
    def left_alignment(self) -> LeftAlignmentReport:
        """Whether this run's indel keys can be trusted to join (GP-14).

        The expensive direction: a frequency that fails to join is scored as no
        frequency data at all, and absence of frequency is not evidence of rarity —
        but it is the input to the rarity signal, which is the strongest promoting
        term the ranker has. Counted over the variants this adapter has been asked
        about rather than over the table's own rows. See
        :func:`_left_alignment_report`.
        """
        return _left_alignment_report(self._indel_queries)

    def frequencies(
        self, variant_ids: Sequence[str]
    ) -> Mapping[str, tuple[PopulationFrequency, ...]]:
        """Return frequency observations for variants present in the table.

        A variant absent from the table is absent from this mapping. It has no
        frequency data — which is emphatically not an allele frequency of zero, and
        not evidence of rarity (GP-14). The caller must handle the missing key.

        Looked up by :func:`_join_key`, returned under the caller's own ID: see
        :meth:`LocalConsequenceAdapter.annotate`.
        """
        found: dict[str, tuple[PopulationFrequency, ...]] = {}
        for variant_id in _unique_ids(variant_ids):
            key = _join_key(variant_id)
            # Counted on every query, hit or miss: a miss is exactly the shape a
            # representation mismatch takes, and it is the case the left-alignment
            # report exists to qualify.
            if _key_is_indel(key):
                self._indel_queries += 1
            entry = self._index.get(key)
            if entry is not None:
                found[variant_id] = entry
        return found


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
        # ADR 0018: the table's spelling is canonicalised, not trusted. Two rows
        # spelling one variant two ways merge into that variant's annotations
        # rather than becoming two entries only one of which a query can find.
        grouped.setdefault(_join_key(variant_id), []).append(annotation)
    return {
        variant_id: tuple(sorted(annotations, key=_consequence_sort_key))
        for variant_id, annotations in grouped.items()
    }


def _index_frequencies(table_path: Path) -> dict[str, tuple[PopulationFrequency, ...]]:
    grouped: dict[str, list[PopulationFrequency]] = {}
    for lineno, row in _read_rows(table_path, _FREQUENCY_COLUMNS):
        variant_id = _required_cell(row, "variant_id", path=table_path, lineno=lineno)
        grouped.setdefault(_join_key(variant_id), []).append(  # ADR 0018
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
