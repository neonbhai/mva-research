"""VCF ingestion: turn a case VCF into typed, build-anchored ``VariantRecord``s.

Two interchangeable backends are provided:

* ``text``  — a dependency-free pure-Python parser. **This is the default**, on
  measurement rather than taste: see :func:`detect_backend`.
* ``cyvcf2`` — htslib under the hood, handles bgzip/BCF, and is what most
  genomics code reaches for first. It produces *identical* records; it is simply
  slower here, because this stage needs the text representation htslib works hard
  to decode away. ``tests/unit/test_ingestion.py::
  test_backends_produce_identical_records_for_the_fixture`` locks the equivalence
  in place.

Both backends emit the same private ``_RawRecord`` intermediate and share one
conversion routine, which is the only reason the equivalence is maintainable: the
interesting logic (contig normalisation, allele validation, multiallelic
decomposition, zygosity derivation, AD assignment) exists exactly once.

Two entry points, one parser
----------------------------

:func:`read_vcf` materialises the whole callset and globally sorts it. That is
correct for any input and is what the fixture-scale callers and every existing
test use — but it costs ~4.2 KB of resident memory per record, which is ~19 GB
for a 4.5 M-variant WGS callset (``docs/scale-report.md`` §2).

:func:`iter_vcf` streams. It reads a bgzip+tabix-indexed VCF **one contig at a
time, in explicit karyotype order**, so the global sort that forced full
materialisation is never needed: records within a contig are already
position-sorted in the file, and the contig visit order is fixed here rather than
inherited from the file, the index or a dict. Peak memory becomes a function of
what the *caller* retains, not of the callset size.

Three rules shape everything here.

**GP-13 — filtering and ranking are separate.** The only records dropped are the
ones that are *un-analysable*: a non-canonical contig, a symbolic/structural or
spanning-deletion allele, a malformed line. Anything merely *suspicious* — low
depth, a caller FILTER, skewed allele balance — survives and is flagged downstream
by :mod:`mva.ingestion.qc`.

**GP-30 — repeat runs are byte-identical.** Streaming is where determinism usually
dies, so the ordering contract is stated rather than assumed. See
:data:`STRATEGY_INDEXED` and :meth:`VcfStream.summary`.

**PRIV-09 — the reader never echoes what it read.** ``warnings``,
``skipped_reasons`` and every exception message carry reason codes and counts
only. No VCF line, no sample name, no genotype string, no coordinate ever reaches
them, because a traceback or a log line travels far further than the artifact it
describes.
"""

from __future__ import annotations

import contextlib
import gzip
import importlib.util
from collections import Counter
from collections.abc import Generator, Iterable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final, TextIO, cast

from pydantic import ValidationError

from mva.errors import AdapterUnavailableError, GenomeBuildMismatchError, IngestionError
from mva.models.base import error_token
from mva.models.genome import (
    CANONICAL_CONTIGS,
    GenomeBuild,
    GenomicCoordinate,
    contig_sort_key,
    normalise_contig,
)
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
# Ordering strategies
# ---------------------------------------------------------------------------

#: Per-contig reads through a tabix/CSI index, contigs visited in the order of
#: :data:`~mva.models.genome.CANONICAL_CONTIGS`. Memory is bounded by one
#: coordinate's worth of records; nothing else is retained.
STRATEGY_INDEXED: Final = "indexed_per_contig"

#: Read everything, then sort globally. Correct for any input — including an
#: unsorted or contig-shuffled one — and is what :func:`read_vcf` always uses.
#: Resident memory is proportional to the callset.
STRATEGY_BUFFERED: Final = "buffered_sort"

#: Records per batch from :meth:`VcfStream.chunks`. Large enough to amortise a
#: per-batch call into a downstream stage, small enough that a batch is a few MB.
DEFAULT_CHUNK_SIZE: Final = 4096

#: Minimum coordinate gap at which :meth:`VcfStream.chunks` is allowed to cut.
#: A downstream per-chunk transform may move a record's position; concatenating
#: independently-sorted chunks stays globally sorted only if no chunk's output can
#: overtake the next chunk's. ``mva.alleles.MAX_SHIFT_BP`` (1000) bounds
#: left-alignment's leftward movement and trimming only moves a record right by
#: less than its REF length, so 4096 clears both with room to spare.
#: ``tests/unit/test_streaming.py`` asserts that relationship against the real
#: constant rather than trusting this comment.
DEFAULT_CHUNK_BOUNDARY_GAP: Final = 4096

#: Hard ceiling on a chunk. Reached only if records stay denser than one per
#: ``boundary_gap`` for this many rows, at which point the batch is no longer
#: bounded memory. Fail closed rather than cut at an unsafe boundary and emit a
#: silently mis-ordered stream.
DEFAULT_CHUNK_MAX_SIZE: Final = 250_000

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

#: Index suffixes appended to the VCF path, in the order they are looked for.
_INDEX_SUFFIXES: Final[tuple[str, ...]] = (".tbi", ".csi")

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

#: FORMAT tags this stage reads, in the order ``_RawRecord`` wants them.
_FORMAT_KEYS: Final[tuple[str, ...]] = ("DP", "AD", "GQ", "PS")


# ---------------------------------------------------------------------------
# Public result types
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class IngestionResult:
    """Everything one VCF yielded, plus an honest account of what it did not.

    ``skipped_count`` and ``skipped_reasons`` exist so that "11 lines in, 12
    variants out" is never a silent claim: a dropped record is reported as an
    aggregate reason code and a count, never as the text that was dropped.

    ``backend`` records which parser actually ran. The two backends agree on every
    well-formed record but not on what they *do* with a malformed one — htslib
    refuses the file, the text parser counts the line under
    ``malformed_record`` — so "which one ran" is a provenance fact, not a
    configuration detail.
    """

    variants: tuple[VariantRecord, ...]
    declared_build: GenomeBuild | None
    sample_id: str
    warnings: tuple[str, ...]
    skipped_count: int
    skipped_reasons: tuple[str, ...]
    backend: str = BACKEND_TEXT


@dataclass(frozen=True, slots=True)
class IngestionSummary:
    """What a :class:`VcfStream` observed. Everything except the records.

    This is the streaming counterpart of :class:`IngestionResult`: identical
    warning, skip and provenance reporting, with ``record_count`` standing in for
    the records themselves. ``complete`` is ``False`` when the stream was
    abandoned before its end, because a partial skip tally understates how many
    variants were discarded and an uncounted discard is unrecoverable.
    """

    declared_build: GenomeBuild | None
    sample_id: str
    backend: str
    strategy: str
    record_count: int
    warnings: tuple[str, ...]
    skipped_count: int
    skipped_reasons: tuple[str, ...]
    complete: bool


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
    """Return the default backend: ``"text"``.

    **Measured, not assumed.** On the 4,492,805-line WGS phantom
    (``docs/scale-report.md`` §1) the pure-Python text parser read the whole file
    in 142.4 s (32,005 rec/s) against cyvcf2's 453.6 s (10,044 rec/s) — 3.2x
    slower — and at the 1 M point that did not page, 55,298 rec/s against 42,545.
    Isolating the body parse on the 100 k file makes the reason plain: htslib's
    own iteration is fast (0.25 s), but this stage needs the *text*
    representation, and converting htslib's decoded arrays back into it costs more
    than parsing the text directly (0.89 s for the cyvcf2 adapter against 0.34 s
    for the text parser; 0.51 s after the adapter was fixed to stop allocating a
    NumPy array per FORMAT tag per record).

    cyvcf2 remains fully supported via ``backend="cyvcf2"`` and is the right choice
    for BCF input or for a file the text parser cannot open; it is simply not the
    default. :func:`cyvcf2_available` is the availability probe.

    Uses :func:`importlib.util.find_spec` rather than a trial import so that
    selecting a backend never has import side effects.
    """
    return BACKEND_TEXT


def cyvcf2_available() -> bool:
    """Whether the native backend can be imported, without importing it."""
    try:
        return importlib.util.find_spec(BACKEND_CYVCF2) is not None
    except (ImportError, ValueError):  # pragma: no cover - defensive
        return False


def _resolve_backend(backend: str) -> str:
    if backend == BACKEND_AUTO:
        return detect_backend()
    if backend not in SUPPORTED_BACKENDS:
        msg = f"Unknown VCF backend {backend!r}; expected one of {list(SUPPORTED_BACKENDS)}."
        raise IngestionError(msg)
    if backend == BACKEND_CYVCF2 and not cyvcf2_available():
        msg = (
            "Backend 'cyvcf2' was requested explicitly but the package is not "
            "importable. Install the 'genomics' extra, or pass backend='text' to use "
            "the pure-Python parser."
        )
        raise AdapterUnavailableError(msg)
    return backend


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------


def read_vcf(
    path: Path,
    *,
    expected_build: GenomeBuild,
    source_artifact: str,
    backend: str = BACKEND_AUTO,
) -> IngestionResult:
    """Read a VCF into sorted, build-anchored variant records, all at once.

    Always uses :data:`STRATEGY_BUFFERED`: every record is materialised and then
    globally sorted, which is correct for *any* input, including one whose contigs
    are shuffled or whose positions are unsorted. It is also ~4.2 KB of resident
    memory per record. For a WGS callset use :func:`iter_vcf`.

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
    stream = iter_vcf(
        path,
        expected_build=expected_build,
        source_artifact=source_artifact,
        backend=backend,
        streaming=False,
    )
    variants = tuple(stream)
    summary = stream.summary()
    return IngestionResult(
        variants=variants,
        declared_build=summary.declared_build,
        sample_id=summary.sample_id,
        warnings=summary.warnings,
        skipped_count=summary.skipped_count,
        skipped_reasons=summary.skipped_reasons,
        backend=summary.backend,
    )


def iter_vcf(
    path: Path,
    *,
    expected_build: GenomeBuild,
    source_artifact: str,
    backend: str = BACKEND_AUTO,
    streaming: bool = True,
) -> VcfStream:
    """Open a VCF for streaming reads. Header now, records on demand.

    The returned :class:`VcfStream` yields exactly the records :func:`read_vcf`
    would return, in exactly the same order, with exactly the same warnings and
    skip reasons — see ``tests/unit/test_streaming.py``, which asserts that
    element by element on a fixture rather than asserting it in prose. What
    differs is that the records are never all alive at once.

    Ordering, and why it is not inherited from the file (GP-30):

    * Contigs are visited in the order of
      :data:`~mva.models.genome.CANONICAL_CONTIGS` — chr1..chr22, chrX, chrY, chrM
      — written out explicitly here, followed by any non-canonical contigs the
      index lists (they yield no records; they are visited so their
      ``non_canonical_contig`` skips are still counted). The index's own contig
      order, `dict` iteration order and `PYTHONHASHSEED` cannot reach this.
    * Within a contig the file is already position-sorted — tabix cannot index it
      otherwise — and that is *verified* as the stream runs, not assumed: a
      position that goes backwards raises rather than emitting a mis-ordered
      stream.
    * Records sharing one coordinate (a multiallelic site's products, or two lines
      at the same POS) are buffered and ordered by ``(ref, alt)`` exactly as the
      global sort would, so a `A>G,C` site still emits `C` before `G`.

    Streaming needs a tabix/CSI index and ``pysam``. Without either — or when the
    file's own contig blocks are not in karyotype order, which would make
    ``source_line_index`` disagree with a sequential read — the stream falls back
    to :data:`STRATEGY_BUFFERED` and says so in :attr:`IngestionSummary.strategy`.
    The records and tallies are identical either way; only the memory profile
    changes.

    Args:
        streaming: set ``False`` to force :data:`STRATEGY_BUFFERED`. Used by
            :func:`read_vcf` and by the tests that compare the two paths.

    Raises:
        GenomeBuildMismatchError: as :func:`read_vcf`.
        IngestionError: as :func:`read_vcf`, and additionally if the input turns
            out not to be coordinate-sorted while streaming.
        AdapterUnavailableError: as :func:`read_vcf`.
    """
    resolved = _resolve_backend(backend)
    if not path.is_file():
        msg = "VCF input does not exist or is not a regular file."
        raise IngestionError(msg)

    header = _read_header(path)
    tally = _Tally.empty()
    _check_declared_build(header, expected_build, tally)
    sample_id = _resolve_sample(header, tally)

    contigs = _contig_visit_plan(path) if streaming else None
    return VcfStream(
        path=path,
        header=header,
        sample_id=sample_id,
        backend=resolved,
        contigs=contigs,
        build=expected_build,
        source_artifact=source_artifact,
        tally=tally,
    )


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
# The stream
# ---------------------------------------------------------------------------


class VcfStream:
    """A VCF opened for reading: header available now, records pulled on demand.

    Iterate it once. The warning and skip tallies are only complete once iteration
    has run to the end, which is why :meth:`summary` refuses to report them before
    then — a truncated skip count would understate how many variants were
    discarded, and a discarded variant nobody counted is invisible (GP-14).

    Construct via :func:`iter_vcf`.
    """

    __slots__ = (
        "_build",
        "_contigs",
        "_count",
        "_exhausted",
        "_header",
        "_path",
        "_ranks",
        "_resolved_backend",
        "_sample_id",
        "_source_artifact",
        "_started",
        "_tally",
    )

    def __init__(
        self,
        *,
        path: Path,
        header: _VcfHeader,
        sample_id: str,
        backend: str,
        contigs: tuple[str, ...] | None,
        build: GenomeBuild,
        source_artifact: str,
        tally: _Tally,
    ) -> None:
        self._path = path
        self._header = header
        self._sample_id = sample_id
        self._resolved_backend = backend
        self._contigs = contigs
        self._build = build
        self._source_artifact = source_artifact
        self._tally = tally
        self._ranks: dict[str, int] = {}
        self._started = False
        self._exhausted = False
        self._count = 0

    # ---------------------------------------------------------------- header

    @property
    def path(self) -> Path:
        return self._path

    @property
    def declared_build(self) -> GenomeBuild | None:
        return self._header.declared_build

    @property
    def sample_id(self) -> str:
        return self._sample_id

    @property
    def backend(self) -> str:
        """The backend that will actually parse the body."""
        return self._resolved_backend

    @property
    def strategy(self) -> str:
        """:data:`STRATEGY_INDEXED` or :data:`STRATEGY_BUFFERED`, decided on open."""
        return STRATEGY_INDEXED if self._contigs is not None else STRATEGY_BUFFERED

    @property
    def exhausted(self) -> bool:
        return self._exhausted

    # ---------------------------------------------------------------- records

    def __iter__(self) -> Iterator[VariantRecord]:
        if self._started:
            msg = (
                "This VCF stream has already been iterated. A stream is single-pass, "
                "and re-iterating would double-count every warning and skip. Call "
                "iter_vcf() again for a second pass, or materialise with read_vcf()."
            )
            raise IngestionError(msg)
        self._started = True
        return self._drive()

    def _drive(self) -> Iterator[VariantRecord]:
        count = 0
        for record in self._ordered():
            count += 1
            yield record
        self._count = count
        self._exhausted = True

    def chunks(
        self,
        *,
        target_size: int = DEFAULT_CHUNK_SIZE,
        boundary_gap: int = DEFAULT_CHUNK_BOUNDARY_GAP,
        max_size: int = DEFAULT_CHUNK_MAX_SIZE,
    ) -> Iterator[tuple[VariantRecord, ...]]:
        """Yield the same records, batched at coordinate gaps wide enough to be safe.

        A downstream stage that transforms a *batch* and sorts its own output —
        :func:`mva.ingestion.normalise.normalise_variants` is the motivating one —
        produces a globally sorted stream when the batches are concatenated only if
        no batch's output can overtake the next batch's. Trimming moves a record
        right by less than its REF length; left-alignment moves it left by at most
        ``mva.alleles.MAX_SHIFT_BP``. So a cut is only made where the next record is
        more than ``boundary_gap`` bases away, or on a contig change, which is what
        makes "normalise each chunk, concatenate" equal to "normalise everything,
        sort once".

        Args:
            target_size: cut at the first safe boundary at or after this many rows.
            boundary_gap: the gap, in bases, that makes a boundary safe. Must exceed
                the largest coordinate movement any per-chunk consumer can apply.
            max_size: refuse rather than cut unsafely. Reached only where records
                stay denser than one per ``boundary_gap`` for this many rows.

        Raises:
            IngestionError: ``max_size`` was reached with no safe boundary in sight.
        """
        batch: list[VariantRecord] = []
        previous: VariantRecord | None = None
        for record in self:
            if (
                previous is not None
                and len(batch) >= target_size
                and _is_safe_cut(previous, record, boundary_gap)
            ):
                yield tuple(batch)
                batch = []
            if len(batch) >= max_size:
                msg = (
                    f"No safe chunk boundary was found within {max_size} records at a "
                    f"gap of {boundary_gap} bases. Cutting anyway would let one chunk's "
                    "normalised output overtake the next chunk's and silently "
                    "mis-order the stream, so the read stops instead. Raise "
                    "boundary_gap, raise max_size, or read this file with read_vcf(). "
                    f"Region handle <region:{error_token(record.coordinate.contig)}>."
                )
                raise IngestionError(msg)
            batch.append(record)
            previous = record
        if batch:
            yield tuple(batch)

    # ---------------------------------------------------------------- tallies

    def summary(self) -> IngestionSummary:
        """Header facts, counts and reason codes — valid only once fully consumed.

        Raises:
            IngestionError: the stream was not read to the end. Use
                :meth:`partial_summary` if a partial count is genuinely what you
                want, and label it as one.
        """
        if not self._exhausted:
            msg = (
                "The VCF stream was not read to its end, so its warning and skip "
                "tallies are incomplete. Reporting them as final would understate how "
                "many records were discarded, and a discarded variant that nobody "
                "counted cannot be recovered or even noticed. Iterate the stream to "
                "completion, or call partial_summary() and label the counts partial."
            )
            raise IngestionError(msg)
        return self._summary()

    def partial_summary(self) -> IngestionSummary:
        """The tallies so far, with ``complete=False`` when iteration stopped early."""
        return self._summary()

    def _summary(self) -> IngestionSummary:
        return IngestionSummary(
            declared_build=self._header.declared_build,
            sample_id=self._sample_id,
            backend=self._resolved_backend,
            strategy=self.strategy,
            record_count=self._count,
            warnings=_format_counts(self._tally.warnings),
            skipped_count=sum(self._tally.skipped.values()),
            skipped_reasons=_format_counts(self._tally.skipped),
            complete=self._exhausted,
        )

    # ---------------------------------------------------------------- ordering

    def _ordered(self) -> Iterator[VariantRecord]:
        if self._contigs is None:
            records = list(self._variant_records())
            records.sort(key=lambda record: record.sort_key())
            yield from records
            return
        yield from self._streamed()

    def _streamed(self) -> Iterator[VariantRecord]:
        """Emit in sort_key order using one coordinate's worth of buffer.

        The whole memory argument lives in this method: the only thing held is the
        group of records sharing the current ``(contig, position)``, because that is
        the only ambiguity a coordinate-sorted file leaves. Everything coarser is
        settled by the contig visit plan; everything finer is settled by sorting the
        group on ``(ref, alt)``, which is what the global sort's key does after
        position, with Python's stable sort preserving file order on a full tie.
        """
        group: list[VariantRecord] = []
        group_key: tuple[int, int] = (-1, -1)
        for record in self._variant_records():
            coordinate = record.coordinate
            key = (self._rank(coordinate.contig), coordinate.position)
            if key != group_key:
                if key < group_key:
                    raise self._unsorted_error(coordinate.contig)
                yield from _ordered_group(group)
                group = []
                group_key = key
            group.append(record)
        yield from _ordered_group(group)

    def _rank(self, contig: str) -> int:
        rank = self._ranks.get(contig)
        if rank is None:
            rank = contig_sort_key(contig)
            self._ranks[contig] = rank
        return rank

    def _unsorted_error(self, contig: str) -> IngestionError:
        return IngestionError(
            "Records are not in coordinate order within a contig, which a "
            "tabix-indexed VCF cannot be. Streaming relies on that order to emit "
            "sorted output without holding the callset, so it stops here rather than "
            "emit a mis-ordered stream. Re-sort and re-index the input, or read it "
            "with read_vcf(), which sorts globally. Contig handle "
            f"<contig:{error_token(contig)}> — the name and position are tokenised "
            "rather than echoed (PRIV-09)."
        )

    # ---------------------------------------------------------------- parsing

    def _variant_records(self) -> Iterator[VariantRecord]:
        build = self._build
        source_artifact = self._source_artifact
        tally = self._tally
        for raw in self._raw_records():
            yield from _records_from_raw(
                raw, build=build, source_artifact=source_artifact, tally=tally
            )

    def _raw_records(self) -> Iterator[_RawRecord]:
        if self._resolved_backend == BACKEND_CYVCF2:
            return _iter_raw_cyvcf2(self._path, self._contigs)
        return _iter_raw_text(self._path, self._contigs, self._tally)


def _ordered_group(group: list[VariantRecord]) -> Iterable[VariantRecord]:
    """Order records sharing one coordinate exactly as the global sort would."""
    if len(group) < 2:
        return group
    return sorted(group, key=lambda record: (record.coordinate.ref, record.coordinate.alt))


def _is_safe_cut(previous: VariantRecord, following: VariantRecord, boundary_gap: int) -> bool:
    if previous.coordinate.contig != following.coordinate.contig:
        return True
    return following.coordinate.position - previous.coordinate.position > boundary_gap


# ---------------------------------------------------------------------------
# Contig visit plan
# ---------------------------------------------------------------------------


def _contig_visit_plan(path: Path) -> tuple[str, ...] | None:
    """The karyotype-ordered contig plan for a streaming read, or ``None``.

    ``None`` means "this file cannot be streamed" and the caller falls back to
    :data:`STRATEGY_BUFFERED`. That happens when there is no tabix/CSI index, when
    ``pysam`` is absent, or when the file's own contig blocks are not already in
    karyotype order.

    That last condition is about ``source_line_index``, not about output order.
    The plan is built from :data:`~mva.models.genome.CANONICAL_CONTIGS`, so the
    *records* would come out correctly ordered whatever the file does; but
    ``source_line_index`` is the record's ordinal among the file's data lines, and
    a per-contig read can only reproduce it by counting lines in the order the file
    stores them. Rather than emit a provenance field that quietly means something
    different from the one :func:`read_vcf` emits, streaming declines.
    """
    if not any(Path(f"{path}{suffix}").is_file() for suffix in _INDEX_SUFFIXES):
        return None
    try:
        if importlib.util.find_spec("pysam") is None:
            return None
    except (ImportError, ValueError):  # pragma: no cover - defensive
        return None

    indexed = _indexed_contigs(path)
    if not indexed:
        return None

    by_rank: dict[int, list[str]] = {}
    extras: list[str] = []
    for name in indexed:
        canonical = _canonical_contig(name)
        if canonical is None:
            extras.append(name)
        else:
            by_rank.setdefault(contig_sort_key(canonical), []).append(name)

    # Written out from CANONICAL_CONTIGS, in that literal order. Not from the
    # index, not from a dict's iteration order, not from the file (GP-30).
    plan: list[str] = []
    for canonical in CANONICAL_CONTIGS:
        plan.extend(by_rank.get(contig_sort_key(canonical), ()))
    plan.extend(extras)

    if plan != list(indexed):
        return None
    return tuple(plan)


def _indexed_contigs(path: Path) -> tuple[str, ...]:
    """Contig names the index lists, in the order the index stores them (file order)."""
    import pysam  # noqa: PLC0415 - optional native index reader, imported on demand

    try:
        handle = pysam.TabixFile(str(path))
    except (OSError, ValueError):  # pragma: no cover - unreadable index
        return ()
    try:
        return tuple(str(name) for name in handle.contigs)
    finally:
        handle.close()


_CONTIG_CANONICAL: dict[str, str | None] = {}


def _canonical_contig(raw: str) -> str | None:
    """:func:`normalise_contig` memoised; ``None`` for a non-canonical contig.

    ``normalise_contig`` is a regex match plus string surgery and the profile in
    ``docs/scale-report.md`` §3 caught it running three times per record. A VCF
    holds a few dozen distinct contig strings, so a dict lookup answers all but the
    first of those calls. The mapping is a pure function of its input, so memoising
    it cannot make a run depend on anything but the input (GP-30).
    """
    try:
        return _CONTIG_CANONICAL[raw]
    except KeyError:
        pass
    try:
        canonical: str | None = normalise_contig(raw)
    except ValueError:
        canonical = None
    _CONTIG_CANONICAL[raw] = canonical
    return canonical


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


def _text_lines(path: Path, contigs: tuple[str, ...] | None) -> Iterator[str]:
    """Body lines, either sequentially or one contig at a time through tabix.

    Both forms yield exactly the file's data lines, in the file's own order, which
    is what makes the running ordinal in :func:`_iter_raw_text` the same
    ``source_line_index`` either way.
    """
    if contigs is None:
        with _open_text(path) as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                yield line
        return

    import pysam  # noqa: PLC0415 - optional native index reader, imported on demand

    try:
        handle = pysam.TabixFile(str(path))
    except (OSError, ValueError) as exc:
        msg = "The tabix index could not be opened; it may be absent or truncated."
        raise IngestionError(msg) from exc
    try:
        for contig in contigs:
            for line in handle.fetch(contig):
                yield str(line)
    except (OSError, ValueError) as exc:
        msg = "The tabix index disagreed with the VCF body; the pair may be stale."
        raise IngestionError(msg) from exc
    finally:
        handle.close()


def _iter_raw_text(
    path: Path, contigs: tuple[str, ...] | None, tally: _Tally
) -> Iterator[_RawRecord]:
    """Pure-Python VCF body parser. No third-party dependency, no network."""
    for index, line in enumerate(_text_lines(path, contigs)):
        raw = _raw_from_text_line(line, index)
        if raw is None:
            tally.skipped[SKIP_MALFORMED_RECORD] += 1
        else:
            yield raw


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


def _iter_raw_cyvcf2(path: Path, contigs: tuple[str, ...] | None) -> Iterator[_RawRecord]:
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

    try:
        index = 0
        if contigs is None:
            for variant in reader:
                yield _raw_from_cyvcf2(variant, index)
                index += 1
        else:
            for contig in contigs:
                for variant in reader(contig):
                    yield _raw_from_cyvcf2(variant, index)
                    index += 1
    except (OSError, ValueError) as exc:
        msg = "htslib failed while iterating VCF records; the input may be malformed."
        raise IngestionError(msg) from exc
    finally:
        reader.close()


def _raw_from_cyvcf2(variant: object, index: int) -> _RawRecord:
    """Adapt one htslib record.

    The shape of this function is the whole cyvcf2 performance story
    (``docs/scale-report.md`` §1). It used to call ``variant.format(key)`` four
    times per record — each allocating a NumPy array that was unpacked element by
    element with ``int()`` and thrown away — including for tags the record does not
    carry, where the call raises and the exception is caught. Now the record's own
    FORMAT list is consulted once, absent tags are never asked for, and each array
    is converted with a single ``tolist()``. Measured on the 100 k phantom, the
    body parse fell from 0.89 s to 0.51 s. It is still slower than the text parser
    at 0.34 s, which is why :func:`detect_backend` returns ``"text"``.
    """
    filters = getattr(variant, "FILTERS", None) or ()
    depth, allelic_depths, genotype_quality, phase_set = _cyvcf2_sample_fields(variant)
    return _RawRecord(
        contig=str(getattr(variant, "CHROM", "")).strip(),
        position=int(getattr(variant, "POS", 0)),
        ref=str(getattr(variant, "REF", "")).strip(),
        alts=tuple(str(alt).strip() for alt in getattr(variant, "ALT", ())),
        quality=_round_quality(_float_or_none(getattr(variant, "QUAL", None))),
        filters=tuple(str(item).strip() for item in filters),
        genotype_string=_gt_string_from_cyvcf2(variant),
        depth=_first_or_none(depth),
        allelic_depths=allelic_depths,
        genotype_quality=_first_or_none(genotype_quality),
        phase_set=_first_or_none(phase_set),
        index=index,
    )


def _cyvcf2_sample_fields(
    variant: object,
) -> tuple[
    tuple[int, ...] | None, tuple[int, ...] | None, tuple[int, ...] | None, tuple[int, ...] | None
]:
    """First sample's values for DP, AD, GQ and PS, in that order.

    htslib encodes missing integers as large negative sentinels, so any negative
    value means "not reported" rather than a real count and the whole tag is
    discarded — a fabricated 0 would read as a measured zero.
    """
    present = getattr(variant, "FORMAT", ())
    fmt = getattr(variant, "format", None)
    values: list[tuple[int, ...] | None] = []
    for key in _FORMAT_KEYS:
        values.append(None if fmt is None else _cyvcf2_vector(fmt, key, present))
    return (values[0], values[1], values[2], values[3])


def _cyvcf2_vector(fmt: object, key: str, present: object) -> tuple[int, ...] | None:
    if key not in present:  # type: ignore[operator]
        return None
    try:
        array = fmt(key)  # type: ignore[operator]
    except (KeyError, ValueError):  # pragma: no cover - FORMAT list already checked
        return None
    if array is None:
        return None
    try:
        row = cast("list[int]", array[0].tolist())
    except (AttributeError, IndexError, TypeError, ValueError):  # pragma: no cover - defensive
        return None
    if not row or any(value < 0 for value in row):
        return None
    return tuple(row)


def _first_or_none(values: tuple[int, ...] | None) -> int | None:
    return values[0] if values else None


#: GT strings repeat endlessly across a callset — a WGS file holds a handful of
#: distinct ones — so the join is done once per distinct call rather than per
#: record. Keyed on htslib's own ``[allele, ..., phased]`` encoding, which is a
#: pure function of the record, so this cannot make a run non-deterministic.
_GT_STRINGS: dict[tuple[int, ...], str] = {}


def _gt_string_from_cyvcf2(variant: object) -> str | None:
    """Rebuild the VCF GT text from cyvcf2's ``[allele, ..., phased]`` encoding."""
    genotypes = getattr(variant, "genotypes", None)
    if not genotypes:
        return None
    try:
        call = tuple(int(value) for value in genotypes[0])
    except (IndexError, TypeError, ValueError):  # pragma: no cover - defensive
        return None
    if len(call) < 2:  # pragma: no cover - defensive
        return None
    cached = _GT_STRINGS.get(call)
    if cached is not None:
        return cached
    separator = "|" if bool(call[-1]) else "/"
    rendered = separator.join("." if allele < 0 else str(allele) for allele in call[:-1])
    _GT_STRINGS[call] = rendered
    return rendered


# ---------------------------------------------------------------------------
# Raw record -> VariantRecord
# ---------------------------------------------------------------------------


def _records_from_raw(
    raw: _RawRecord, *, build: GenomeBuild, source_artifact: str, tally: _Tally
) -> list[VariantRecord]:
    """Decompose one VCF line into zero or more biallelic ``VariantRecord``s."""
    contig = _canonical_contig(raw.contig)
    if contig is None:
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
    """Derive this allele's call, keeping the site's original GT text verbatim.

    ``allele_index`` is 1-based over the site's ALT list and is stored on the
    genotype, because the verbatim GT alone cannot say which ALT this record is.
    """
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
        # The GT text is kept verbatim, so after a multiallelic split it no longer
        # says WHICH of its allele numbers is this record's. Recording the index
        # is what lets a phased '1|2' still resolve a haplotype slot instead of
        # reading as "both haplotypes carry an alternate allele".
        alt_allele_index=allele_index,
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
