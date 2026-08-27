"""VCF ingestion: turn a case VCF into typed, build-anchored ``VariantRecord``s.

Two interchangeable backends are provided:

* ``cyvcf2`` — htslib under the hood, handles bgzip/BCF and is what production uses;
* ``text``  — a dependency-free pure-Python parser, so the pipeline still runs (and
  still produces *identical* records) on a machine where the native wheel will not
  build. ``tests/unit/test_ingestion.py::test_backends_agree_on_fixture`` locks the
  equivalence in place.

Both backends emit the same private ``_RawRecord`` intermediate and share one
conversion routine, which is the only reason the equivalence is maintainable: the
interesting logic (contig normalisation, allele validation, multiallelic
decomposition, zygosity derivation, AD assignment) exists exactly once.

Two rules shape everything here.

**GP-13 — filtering and ranking are separate.** The only records dropped are the
ones that are *un-analysable*: a non-canonical contig, a symbolic/structural or
spanning-deletion allele, a malformed line. Anything merely *suspicious* — low
depth, a caller FILTER, skewed allele balance — survives and is flagged downstream
by :mod:`mva.ingestion.qc`.

**PRIV-09 — the reader never echoes what it read.** ``warnings``, ``skipped_reasons``
and every exception message carry reason codes and counts only. No VCF line, no
sample name, no genotype string, no coordinate ever reaches them, because a
traceback or a log line travels far further than the artifact it describes.
"""

from __future__ import annotations

import contextlib
import gzip
import importlib.util
from collections import Counter
from collections.abc import Generator, Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final, TextIO, cast

from pydantic import ValidationError

from mva.errors import AdapterUnavailableError, GenomeBuildMismatchError, IngestionError
from mva.models.genome import GenomeBuild, GenomicCoordinate, normalise_contig
from mva.models.variant import FilterStatus, Genotype, VariantRecord, Zygosity

# ---------------------------------------------------------------------------
# Backend identifiers
# ---------------------------------------------------------------------------

BACKEND_AUTO: Final = "auto"
BACKEND_CYVCF2: Final = "cyvcf2"
BACKEND_TEXT: Final = "text"

#: Every backend token ``read_vcf`` accepts.
SUPPORTED_BACKENDS: Final[tuple[str, ...]] = (BACKEND_AUTO, BACKEND_CYVCF2, BACKEND_TEXT)

# ---------------------------------------------------------------------------
# Reason codes. Deliberately coordinate-free and record-free (PRIV-09).
# ---------------------------------------------------------------------------

SKIP_NON_CANONICAL_CONTIG: Final = "non_canonical_contig"
SKIP_SYMBOLIC_ALLELE: Final = "symbolic_or_structural_allele"
SKIP_SPANNING_DELETION: Final = "spanning_deletion_allele"
SKIP_MISSING_ALT: Final = "missing_alt_allele"
SKIP_MALFORMED_RECORD: Final = "malformed_record"
SKIP_INVALID_COORDINATE: Final = "invalid_coordinate"

WARN_BUILD_NOT_DECLARED: Final = "genome_build_not_declared_in_header"
WARN_BUILD_UNRECOGNISED: Final = "genome_build_header_unrecognised"
WARN_NO_SAMPLE_COLUMN: Final = "no_sample_column"
WARN_MULTIPLE_SAMPLES: Final = "multiple_sample_columns_first_used"
WARN_MISSING_GENOTYPE: Final = "missing_format_gt"
WARN_MISSING_ALLELIC_DEPTH: Final = "missing_format_ad"
WARN_ALLELIC_DEPTH_LENGTH: Final = "allelic_depth_length_mismatch"
WARN_NO_CALL_GENOTYPE: Final = "no_call_genotype"

_GZIP_SUFFIXES: Final[frozenset[str]] = frozenset({".gz", ".bgz", ".bgzf"})
_MISSING_TOKENS: Final[frozenset[str]] = frozenset({".", ""})
_NUCLEOTIDES: Final[frozenset[str]] = frozenset("ACGTN")

#: VCF 4.2 spanning-deletion ALT: real syntax, but not an independently rankable
#: candidate, so it is a hard filter rather than a flag.
_SPANNING_DELETION_ALLELE: Final = "*"

#: VCF has 8 fixed columns; genotype data needs FORMAT + at least one sample.
_FIXED_COLUMNS: Final = 8
_FORMAT_COLUMN: Final = 8
_FIRST_SAMPLE_COLUMN: Final = 9

#: QUAL is stored as float32 by htslib but as decimal text in a plain VCF.
#: Rounding both to six significant digits makes the two backends agree exactly
#: (540.2 vs 540.2000122070312) without inventing precision that was never there.
_QUAL_SIGNIFICANT_DIGITS: Final = 6


# ---------------------------------------------------------------------------
# Public result type
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class IngestionResult:
    """Everything one VCF yielded, plus an honest account of what it did not.

    ``skipped_count`` and ``skipped_reasons`` exist so that "11 lines in, 12
    variants out" is never a silent claim: a dropped record is reported as an
    aggregate reason code and a count, never as the text that was dropped.
    """

    variants: tuple[VariantRecord, ...]
    declared_build: GenomeBuild | None
    sample_id: str
    warnings: tuple[str, ...]
    skipped_count: int
    skipped_reasons: tuple[str, ...]


# ---------------------------------------------------------------------------
# Private intermediates shared by both backends
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class _RawRecord:
    """One VCF data line, decoded but not yet interpreted or validated."""

    contig: str
    position: int
    ref: str
    alts: tuple[str, ...]
    quality: float | None
    filters: tuple[str, ...]
    genotype_string: str | None
    depth: int | None
    allelic_depths: tuple[int, ...] | None
    genotype_quality: int | None
    phase_set: int | None
    index: int


@dataclass(frozen=True, slots=True)
class _VcfHeader:
    """The only parts of the header this stage is entitled to act on."""

    declared_build: GenomeBuild | None
    build_header_present: bool
    samples: tuple[str, ...]


@dataclass(slots=True)
class _Tally:
    """Accumulator for privacy-safe aggregate reporting."""

    warnings: Counter[str]
    skipped: Counter[str]

    @classmethod
    def empty(cls) -> _Tally:
        return cls(warnings=Counter(), skipped=Counter())


# ---------------------------------------------------------------------------
# Backend selection
# ---------------------------------------------------------------------------


def detect_backend() -> str:
    """Return ``"cyvcf2"`` when the native backend is importable, else ``"text"``.

    Uses :func:`importlib.util.find_spec` rather than a trial import so that
    selecting a backend never has import side effects.
    """
    try:
        spec = importlib.util.find_spec(BACKEND_CYVCF2)
    except (ImportError, ValueError):  # pragma: no cover - defensive
        return BACKEND_TEXT
    return BACKEND_CYVCF2 if spec is not None else BACKEND_TEXT


def _resolve_backend(backend: str) -> str:
    if backend == BACKEND_AUTO:
        return detect_backend()
    if backend not in SUPPORTED_BACKENDS:
        msg = f"Unknown VCF backend {backend!r}; expected one of {list(SUPPORTED_BACKENDS)}."
        raise IngestionError(msg)
    if backend == BACKEND_CYVCF2 and detect_backend() != BACKEND_CYVCF2:
        msg = (
            "Backend 'cyvcf2' was requested explicitly but the package is not "
            "importable. Install the 'genomics' extra, or pass backend='text' to use "
            "the pure-Python parser."
        )
        raise AdapterUnavailableError(msg)
    return backend


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def read_vcf(
    path: Path,
    *,
    expected_build: GenomeBuild,
    source_artifact: str,
    backend: str = BACKEND_AUTO,
) -> IngestionResult:
    """Read a VCF into sorted, build-anchored variant records.

    A multiallelic site is decomposed here rather than downstream, because a
    ``VariantRecord`` is biallelic by construction (``GenomicCoordinate.alt`` admits
    one allele) and this is the last point at which the full ALT list still exists.
    Each product carries ``normalisation_ops=("split_multiallelic",)`` and the
    per-allele AD slice; the original GT string is preserved verbatim on the
    genotype so that phase and the original call remain auditable.

    Raises:
        GenomeBuildMismatchError: the header declares a build other than
            ``expected_build``. Never coerced — see GP-11.
        IngestionError: the file is absent, unreadable, or the backend token is
            unknown.
        AdapterUnavailableError: ``backend="cyvcf2"`` was demanded but is missing.
    """
    resolved = _resolve_backend(backend)
    if not path.is_file():
        msg = "VCF input does not exist or is not a regular file."
        raise IngestionError(msg)

    header = _read_header(path)
    tally = _Tally.empty()
    _check_declared_build(header, expected_build, tally)
    sample_id = _resolve_sample(header, tally)

    raws = _parse_records(path, resolved, tally)

    records: list[VariantRecord] = []
    for raw in raws:
        records.extend(
            _records_from_raw(
                raw, build=expected_build, source_artifact=source_artifact, tally=tally
            )
        )
    records.sort(key=lambda record: record.sort_key())

    return IngestionResult(
        variants=tuple(records),
        declared_build=header.declared_build,
        sample_id=sample_id,
        warnings=_format_counts(tally.warnings),
        skipped_count=sum(tally.skipped.values()),
        skipped_reasons=_format_counts(tally.skipped),
    )


def _parse_records(path: Path, backend: str, tally: _Tally) -> tuple[_RawRecord, ...]:
    if backend == BACKEND_CYVCF2:
        return _parse_with_cyvcf2(path)
    return _parse_with_text(path, tally)


def _check_declared_build(header: _VcfHeader, expected_build: GenomeBuild, tally: _Tally) -> None:
    """GP-11: a declared build that disagrees is fatal; an absent one is a warning."""
    if header.declared_build is None:
        code = WARN_BUILD_UNRECOGNISED if header.build_header_present else WARN_BUILD_NOT_DECLARED
        tally.warnings[code] += 1
        return
    if header.declared_build is not expected_build:
        msg = (
            f"VCF header declares genome build {header.declared_build.value} but this "
            f"run expects {expected_build.value}. Positions differ by megabases between "
            "assemblies, so the build is never coerced: re-run against the matching "
            "build, or lift over as an explicit, provenance-tracked stage."
        )
        raise GenomeBuildMismatchError(msg)


def _resolve_sample(header: _VcfHeader, tally: _Tally) -> str:
    """Single-sample stage: the first sample column is the proband."""
    samples = header.samples
    if not samples:
        tally.warnings[WARN_NO_SAMPLE_COLUMN] += 1
    elif len(samples) > 1:
        tally.warnings[WARN_MULTIPLE_SAMPLES] += len(samples)
    return next(iter(samples), "")


def _format_counts(counts: Counter[str]) -> tuple[str, ...]:
    """Render a tally as sorted ``"code (n=N)"`` strings — codes and counts only."""
    return tuple(f"{code} (n={counts[code]})" for code in sorted(counts))


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------


@contextlib.contextmanager
def _open_text(path: Path) -> Generator[TextIO]:
    """Open a plain or gzip/bgzip VCF as UTF-8 text."""
    handle: TextIO
    if path.suffix.lower() in _GZIP_SUFFIXES:
        handle = cast("TextIO", gzip.open(path, mode="rt", encoding="utf-8"))  # noqa: SIM115
    else:
        handle = cast("TextIO", path.open("r", encoding="utf-8"))
    try:
        yield handle
    finally:
        handle.close()


def _read_header(path: Path) -> _VcfHeader:
    """Parse the header from the file itself, whichever backend reads the body.

    Sharing this step is what keeps ``declared_build`` and ``sample_id`` identical
    across backends; htslib rewrites and reorders header lines when it re-emits
    them, so ``raw_header`` is deliberately not used.
    """
    lines: list[str] = []
    with _open_text(path) as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            lines.append(line.rstrip("\r\n"))
    return _parse_header_lines(lines)


def _parse_header_lines(lines: Iterable[str]) -> _VcfHeader:
    build_value: str | None = None
    samples: tuple[str, ...] = ()
    for line in lines:
        if line.startswith("#CHROM"):
            fields = line.lstrip("#").split("\t")
            samples = tuple(f.strip() for f in fields[_FIRST_SAMPLE_COLUMN:] if f.strip())
        elif build_value is None:
            for key in ("##reference=", "##assembly="):
                if line.startswith(key):
                    build_value = line[len(key) :].strip()
                    break
    declared = _build_from_header_value(build_value) if build_value else None
    return _VcfHeader(
        declared_build=declared,
        build_header_present=build_value is not None and build_value not in _MISSING_TOKENS,
        samples=samples,
    )


def _build_from_header_value(value: str) -> GenomeBuild | None:
    """Resolve ``##reference=`` to a build, tolerating a path-shaped value.

    ``##reference=GRCh38`` resolves directly. ``##reference=file:///refs/hg38.fa``
    resolves by alias substring — but only for aliases long enough to be
    unambiguous, because matching a bare ``38`` inside a path would be a coin flip,
    and guessing the assembly is exactly the failure GP-11 exists to prevent.
    """
    with contextlib.suppress(ValueError):
        return GenomeBuild.parse(value)
    lowered = value.lower()
    for build in (GenomeBuild.GRCH38, GenomeBuild.GRCH37):
        for alias in sorted(build.aliases):
            if len(alias) > 2 and alias.lower() in lowered:
                return build
    return None


# ---------------------------------------------------------------------------
# Text backend
# ---------------------------------------------------------------------------


def _parse_with_text(path: Path, tally: _Tally) -> tuple[_RawRecord, ...]:
    """Pure-Python VCF body parser. No third-party dependency, no network."""
    records: list[_RawRecord] = []
    index = 0
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            raw = _raw_from_text_line(line, index)
            if raw is None:
                tally.skipped[SKIP_MALFORMED_RECORD] += 1
            else:
                records.append(raw)
            index += 1
    return tuple(records)


def _raw_from_text_line(line: str, index: int) -> _RawRecord | None:
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < _FIXED_COLUMNS:
        return None
    position = _positive_int(fields[1])
    if position is None:
        return None
    sample = _sample_fields(fields)
    return _RawRecord(
        contig=fields[0].strip(),
        position=position,
        ref=fields[3].strip(),
        alts=_split_alts(fields[4]),
        quality=_round_quality(_float_or_none(fields[5])),
        filters=_split_filters(fields[6]),
        genotype_string=sample.get("GT"),
        depth=_int_or_none(sample.get("DP")),
        allelic_depths=_int_tuple_or_none(sample.get("AD")),
        genotype_quality=_int_or_none(sample.get("GQ")),
        phase_set=_int_or_none(sample.get("PS")),
        index=index,
    )


def _sample_fields(fields: Sequence[str]) -> dict[str, str]:
    """Zip FORMAT keys against the first sample column, tolerating truncation."""
    if len(fields) <= _FIRST_SAMPLE_COLUMN:
        return {}
    keys = fields[_FORMAT_COLUMN].strip().split(":")
    values = fields[_FIRST_SAMPLE_COLUMN].strip().split(":")
    return dict(zip(keys, values, strict=False))


def _split_alts(field: str) -> tuple[str, ...]:
    token = field.strip()
    if token in _MISSING_TOKENS:
        return ()
    return tuple(part.strip() for part in token.split(","))


def _split_filters(field: str) -> tuple[str, ...]:
    token = field.strip()
    if token in _MISSING_TOKENS:
        return ()
    return tuple(part.strip() for part in token.split(";") if part.strip())


# ---------------------------------------------------------------------------
# cyvcf2 backend
# ---------------------------------------------------------------------------


def _parse_with_cyvcf2(path: Path) -> tuple[_RawRecord, ...]:
    try:
        from cyvcf2 import VCF  # noqa: PLC0415 - optional native backend, imported on demand
    except ImportError as exc:  # pragma: no cover - guarded by _resolve_backend
        msg = "The cyvcf2 backend is not installed; install the 'genomics' extra."
        raise AdapterUnavailableError(msg) from exc

    try:
        reader = VCF(str(path))
    except (OSError, ValueError) as exc:
        msg = "htslib could not open the VCF input; it may be truncated or malformed."
        raise IngestionError(msg) from exc

    records: list[_RawRecord] = []
    try:
        for index, variant in enumerate(reader):
            records.append(_raw_from_cyvcf2(variant, index))
    except (OSError, ValueError) as exc:
        msg = "htslib failed while iterating VCF records; the input may be malformed."
        raise IngestionError(msg) from exc
    finally:
        reader.close()
    return tuple(records)


def _raw_from_cyvcf2(variant: object, index: int) -> _RawRecord:
    filters = getattr(variant, "FILTERS", None) or ()
    return _RawRecord(
        contig=str(getattr(variant, "CHROM", "")).strip(),
        position=int(getattr(variant, "POS", 0)),
        ref=str(getattr(variant, "REF", "")).strip(),
        alts=tuple(str(alt).strip() for alt in getattr(variant, "ALT", ())),
        quality=_round_quality(_float_or_none(getattr(variant, "QUAL", None))),
        filters=tuple(str(item).strip() for item in filters),
        genotype_string=_gt_string_from_cyvcf2(variant),
        depth=_cyvcf2_scalar(variant, "DP"),
        allelic_depths=_cyvcf2_vector(variant, "AD"),
        genotype_quality=_cyvcf2_scalar(variant, "GQ"),
        phase_set=_cyvcf2_scalar(variant, "PS"),
        index=index,
    )


def _cyvcf2_format(variant: object, key: str) -> object | None:
    """``variant.format(key)`` returns ``None`` when absent and raises when unknown."""
    fmt = getattr(variant, "format", None)
    if fmt is None:  # pragma: no cover - defensive
        return None
    try:
        return fmt(key)
    except (KeyError, ValueError):
        return None


def _cyvcf2_scalar(variant: object, key: str) -> int | None:
    values = _cyvcf2_vector(variant, key)
    return values[0] if values else None


def _cyvcf2_vector(variant: object, key: str) -> tuple[int, ...] | None:
    """First sample's values for an integer FORMAT tag, or ``None`` if unusable.

    htslib encodes missing integers as large negative sentinels, so any negative
    value means "not reported" rather than a real count and the whole tag is
    discarded — a fabricated 0 would read as a measured zero.
    """
    array = _cyvcf2_format(variant, key)
    if array is None:
        return None
    try:
        row = array[0]  # type: ignore[index]
        values = tuple(int(value) for value in row)
    except (IndexError, TypeError, ValueError):  # pragma: no cover - defensive
        return None
    if not values or any(value < 0 for value in values):
        return None
    return values


def _gt_string_from_cyvcf2(variant: object) -> str | None:
    """Rebuild the VCF GT text from cyvcf2's ``[allele, ..., phased]`` encoding."""
    genotypes = getattr(variant, "genotypes", None)
    if not genotypes:
        return None
    try:
        call = list(genotypes[0])
    except (IndexError, TypeError):  # pragma: no cover - defensive
        return None
    if len(call) < 2:  # pragma: no cover - defensive
        return None
    alleles = call[:-1]
    separator = "|" if bool(call[-1]) else "/"
    return separator.join("." if int(a) < 0 else str(int(a)) for a in alleles)


# ---------------------------------------------------------------------------
# Raw record -> VariantRecord
# ---------------------------------------------------------------------------


def _records_from_raw(
    raw: _RawRecord, *, build: GenomeBuild, source_artifact: str, tally: _Tally
) -> list[VariantRecord]:
    """Decompose one VCF line into zero or more biallelic ``VariantRecord``s."""
    try:
        contig = normalise_contig(raw.contig)
    except ValueError:
        tally.skipped[SKIP_NON_CANONICAL_CONTIG] += len(raw.alts) or 1
        return []

    if not raw.alts:
        tally.skipped[SKIP_MISSING_ALT] += 1
        return []
    if raw.genotype_string is None:
        tally.warnings[WARN_MISSING_GENOTYPE] += 1
    if raw.allelic_depths is None:
        tally.warnings[WARN_MISSING_ALLELIC_DEPTH] += 1
    elif len(raw.allelic_depths) != len(raw.alts) + 1:
        tally.warnings[WARN_ALLELIC_DEPTH_LENGTH] += 1

    ops: tuple[str, ...] = ("split_multiallelic",) if len(raw.alts) > 1 else ()
    records: list[VariantRecord] = []
    for offset, alt in enumerate(raw.alts):
        allele_index = offset + 1
        reason = _unusable_allele_reason(alt)
        if reason is not None:
            tally.skipped[reason] += 1
            continue
        record = _build_record(
            raw,
            contig=contig,
            alt=alt.strip().upper(),
            allele_index=allele_index,
            build=build,
            source_artifact=source_artifact,
            ops=ops,
            tally=tally,
        )
        if record is None:
            tally.skipped[SKIP_INVALID_COORDINATE] += 1
            continue
        records.append(record)
    return records


def _unusable_allele_reason(alt: str) -> str | None:
    """Hard-filter reasons only: alleles this pipeline cannot reason about at all."""
    allele = alt.strip().upper()
    if allele in _MISSING_TOKENS:
        return SKIP_MISSING_ALT
    if allele == _SPANNING_DELETION_ALLELE:
        return SKIP_SPANNING_DELETION
    if allele.startswith("<") or "[" in allele or "]" in allele:
        return SKIP_SYMBOLIC_ALLELE
    if not set(allele) <= _NUCLEOTIDES:
        return SKIP_MALFORMED_RECORD
    return None


def _build_record(
    raw: _RawRecord,
    *,
    contig: str,
    alt: str,
    allele_index: int,
    build: GenomeBuild,
    source_artifact: str,
    ops: tuple[str, ...],
    tally: _Tally,
) -> VariantRecord | None:
    try:
        coordinate = GenomicCoordinate(
            build=build, contig=contig, position=raw.position, ref=raw.ref, alt=alt
        )
    except (ValidationError, ValueError):
        return None
    genotype = _genotype_for_allele(raw, allele_index, tally)
    try:
        return VariantRecord(
            coordinate=coordinate,
            genotype=genotype,
            filter_status=_filter_status(raw.filters),
            raw_filters=raw.filters,
            quality=raw.quality,
            source_artifact=source_artifact,
            source_line_index=raw.index,
            normalisation_ops=ops,
        )
    except (ValidationError, ValueError):  # pragma: no cover - defensive
        return None


def _filter_status(filters: tuple[str, ...]) -> FilterStatus:
    if not filters:
        return FilterStatus.MISSING
    if len(filters) == 1 and filters[0].upper() == "PASS":
        return FilterStatus.PASS
    return FilterStatus.FILTERED


def _genotype_for_allele(raw: _RawRecord, allele_index: int, tally: _Tally) -> Genotype:
    """Derive this allele's call, keeping the site's original GT text verbatim."""
    gt_string = raw.genotype_string if raw.genotype_string is not None else "./."
    alleles = _parse_gt(gt_string)
    zygosity = _zygosity_for_allele(alleles, allele_index)
    if zygosity is Zygosity.UNKNOWN:
        tally.warnings[WARN_NO_CALL_GENOTYPE] += 1

    ref_reads, alt_reads = _allele_depths(raw, allele_index)
    depth = raw.depth
    if depth is None and raw.allelic_depths is not None:
        depth = sum(raw.allelic_depths)

    return Genotype(
        zygosity=zygosity,
        genotype_string=gt_string,
        phased="|" in gt_string,
        phase_set=raw.phase_set if raw.phase_set and raw.phase_set > 0 else None,
        depth=depth,
        ref_reads=ref_reads,
        alt_reads=alt_reads,
        genotype_quality=raw.genotype_quality,
    )


def _allele_depths(raw: _RawRecord, allele_index: int) -> tuple[int | None, int | None]:
    """AD is ``[ref, alt1, alt2, ...]``; allele *i* takes ``AD[0]`` and ``AD[i]``."""
    depths = raw.allelic_depths
    if depths is None or len(depths) != len(raw.alts) + 1:
        return (None, None)
    return (depths[0], depths[allele_index])


def _parse_gt(gt_string: str) -> tuple[int | None, ...]:
    tokens = gt_string.strip().replace("|", "/").split("/")
    alleles: list[int | None] = []
    for token in tokens:
        if token in _MISSING_TOKENS:
            alleles.append(None)
            continue
        try:
            alleles.append(int(token))
        except ValueError:
            alleles.append(None)
    return tuple(alleles)


def _zygosity_for_allele(alleles: tuple[int | None, ...], allele_index: int) -> Zygosity:
    """Zygosity of one ALT allele.

    A no-call stays :attr:`Zygosity.UNKNOWN` rather than collapsing to hom-ref:
    "we did not observe this" and "we observed no variant" are different claims and
    only one of them is true.
    """
    if not alleles or any(allele is None for allele in alleles):
        return Zygosity.UNKNOWN
    called = [allele for allele in alleles if allele is not None]
    copies = sum(1 for allele in called if allele == allele_index)
    if len(called) == 1:
        return Zygosity.HEMIZYGOUS if copies == 1 else Zygosity.HOM_REF
    if copies == 0:
        return Zygosity.HOM_REF
    if copies == len(called):
        return Zygosity.HOM_ALT
    return Zygosity.HET


# ---------------------------------------------------------------------------
# Scalar coercions
# ---------------------------------------------------------------------------


def _positive_int(token: str) -> int | None:
    value = _int_or_none(token)
    return value if value is not None and value > 0 else None


def _int_or_none(token: str | None) -> int | None:
    """Parse an integer FORMAT value; negatives and sentinels become ``None``."""
    if token is None:
        return None
    stripped = token.strip()
    if stripped in _MISSING_TOKENS:
        return None
    try:
        value = int(stripped)
    except ValueError:
        try:
            value = round(float(stripped))
        except ValueError:
            return None
    return value if value >= 0 else None


def _int_tuple_or_none(token: str | None) -> tuple[int, ...] | None:
    if token is None:
        return None
    stripped = token.strip()
    if stripped in _MISSING_TOKENS:
        return None
    values: list[int] = []
    for part in stripped.split(","):
        value = _int_or_none(part)
        if value is None:
            return None
        values.append(value)
    return tuple(values) or None


def _float_or_none(value: object) -> float | None:
    if value is None:
        return None
    if isinstance(value, str):
        stripped = value.strip()
        if stripped in _MISSING_TOKENS:
            return None
        try:
            return float(stripped)
        except ValueError:
            return None
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):  # pragma: no cover - defensive
        return None


def _round_quality(value: float | None) -> float | None:
    """Collapse float32/decimal-text disagreement so both backends agree exactly."""
    if value is None:
        return None
    return float(f"{value:.{_QUAL_SIGNIFICANT_DIGITS}g}")
