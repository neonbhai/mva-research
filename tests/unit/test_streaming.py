"""The streaming ingestion path must be indistinguishable from the batch one.

``iter_vcf`` exists for one reason: ``read_vcf`` holds the whole callset, which is
~19 GB for a 4.5 M-variant WGS file (``docs/scale-report.md`` §2). Replacing a
global sort with per-contig iteration is exactly the kind of change that trades a
memory problem for a silent ordering problem, so the contract asserted here is the
strongest one available: **the same records, in the same order, with the same
warnings and the same skip counts**, element by element, on a fixture built to
contain every case the two paths could disagree on.

The rest of the file defends the three ways a streaming rewrite usually goes
wrong: order inherited from the file rather than fixed (GP-30), discarded records
that stop being counted (GP-14), and an exception that names the coordinate it
choked on (PRIV-09).
"""

from __future__ import annotations

import subprocess
import sys
from dataclasses import replace
from datetime import UTC, datetime
from itertools import pairwise
from pathlib import Path

import pytest

from mva.alleles import MAX_SHIFT_BP
from mva.clock import FixedClock
from mva.config import QualityThresholds
from mva.determinism import canonical_json
from mva.errors import IngestionError
from mva.ingestion.normalise import normalise_variants
from mva.ingestion.qc import assess_quality, iter_assessed
from mva.ingestion.reader import (
    BACKEND_CYVCF2,
    BACKEND_TEXT,
    DEFAULT_CHUNK_BOUNDARY_GAP,
    STRATEGY_BUFFERED,
    STRATEGY_INDEXED,
    IngestionSummary,
    iter_vcf,
    read_vcf,
)
from mva.models.genome import CANONICAL_CONTIGS, GenomeBuild
from mva.models.variant import VariantRecord

pysam = pytest.importorskip("pysam")

BUILD = GenomeBuild.GRCH38
ARTIFACT = "streaming-fixture"
CLOCK = FixedClock(datetime(2026, 1, 1, tzinfo=UTC))
THRESHOLDS = QualityThresholds()

_HEADER = (
    "##fileformat=VCFv4.2",
    "##reference=GRCh38",
    '##FILTER=<ID=PASS,Description="All filters passed">',
    '##FILTER=<ID=LowQual,Description="Low quality variant call">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH_PROBAND01",
)

#: One line per way the two read paths could disagree, in the coordinate order a
#: tabix-indexable file must be written in. Contig blocks are in karyotype order
#: (chr1, chr2, chr10, chrX) with the non-canonical block last, which is the shape
#: `_contig_visit_plan` requires before it will stream.
_ROWS: tuple[str, ...] = (
    # plain het
    "chr1\t1000\t.\tA\tG\t820.5\tPASS\t.\tGT:DP:AD:GQ\t0/1:45:23,22:99",
    # multiallelic whose ALT order (T,A) is NOT its sorted order (A,T): the split
    # products share one coordinate and must come out A before T.
    "chr1\t2000\t.\tC\tT,A\t455.0\tPASS\t.\tGT:DP:AD:GQ\t1/2:44:2,21,21:97",
    # a second line at the SAME position with a longer REF: sorts after the two
    # above on `ref`, which only a per-position buffer can get right.
    "chr1\t2000\t.\tCAT\tC\t310.6\tLowQual\t.\tGT:DP:AD:GQ\t0/1:12:7,5:44",
    # spanning deletion: skipped, counted
    "chr1\t3000\t.\tG\tA,*\t500.0\tPASS\t.\tGT:DP:AD:GQ\t0/1:30:15,15,0:80",
    # symbolic allele: skipped, counted
    "chr1\t4000\t.\tG\t<DEL>\t500.0\tPASS\t.\tGT:DP\t0/1:30",
    # malformed line (fewer than the 8 fixed columns): skipped, counted, and it
    # still consumes a source_line_index in both paths
    "chr1\t5000\t.\tC\tT",
    # no FORMAT/sample columns at all
    "chr1\t6000\t.\tT\tC\t100.0\t.\t.",
    # explicit no-call
    "chr1\t7000\t.\tA\tT\t90.0\tPASS\t.\tGT:DP:AD:GQ\t./.:9:5,4:20",
    # missing ALT
    "chr1\t8000\t.\tA\t.\t90.0\tPASS\t.\tGT:DP\t0/0:31",
    # hom-alt, no FILTER opinion
    "chr1\t9000\t.\tT\tC\t880.1\t.\t.\tGT:DP:AD:GQ\t1/1:48:1,47:99",
    # phased, with a phase set
    "chr2\t5000\t.\tC\tA\t540.2\tPASS\t.\tGT:DP:AD:GQ:PS\t1|0:40:21,19:95:5000",
    "chr2\t5050\t.\tT\tG\t515.8\tPASS\t.\tGT:DP:AD:GQ:PS\t1|0:39:20,19:93:5000",
    # chr10 sorts AFTER chr2 by karyotype rank and BEFORE it lexically
    "chr10\t100\t.\tG\tA\t760.4\tPASS\t.\tGT:DP:AD:GQ\t0/1:38:18,20:98",
    "chr10\t900000\t.\tGA\tG\t410.0\tPASS\t.\tGT:DP:AD:GQ\t0/1:22:12,10:70",
    # hemizygous
    "chrX\t200\t.\tC\tT\t640.0\tPASS\t.\tGT:DP:AD:GQ\t1:29:0,29:99",
    # non-canonical contig: skipped, counted, and its block must still be visited
    "chrUn_KI270302v1\t300\t.\tC\tT\t500.0\tPASS\t.\tGT:DP\t0/1:30",
    "chrUn_KI270302v1\t400\t.\tC\tG,T\t500.0\tPASS\t.\tGT:DP\t1/2:30",
)


#: The subset every VCF parser must agree on: at least the 8 fixed columns plus
#: FORMAT and one sample. The two rows this drops are the two the backends handle
#: differently, and they are asserted separately.
_WELL_FORMED_ROWS: tuple[str, ...] = tuple(row for row in _ROWS if row.count("\t") >= 9)


# ---------------------------------------------------------------------------
# Fixture construction
# ---------------------------------------------------------------------------


def _write_plain(path: Path, rows: tuple[str, ...] = _ROWS) -> Path:
    path.write_text("\n".join((*_HEADER, *rows)) + "\n", encoding="utf-8")
    return path


def _write_indexed(directory: Path, rows: tuple[str, ...] = _ROWS, *, name: str = "case") -> Path:
    """bgzip + tabix a VCF, which is what makes the streaming strategy available."""
    plain = _write_plain(directory / f"{name}.vcf", rows)
    compressed = directory / f"{name}.vcf.gz"
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)
    return compressed


@pytest.fixture
def indexed_vcf(tmp_path: Path) -> Path:
    return _write_indexed(tmp_path)


def _stream(
    path: Path, backend: str = BACKEND_TEXT
) -> tuple[list[VariantRecord], IngestionSummary]:
    handle = iter_vcf(path, expected_build=BUILD, source_artifact=ARTIFACT, backend=backend)
    records = list(handle)
    return records, handle.summary()


# ---------------------------------------------------------------------------
# 1. The headline: equal records, equal order
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_streaming_and_batch_produce_equal_records_in_equal_order(indexed_vcf: Path) -> None:
    """The whole point. Every field of every record, in the same position."""
    batch = read_vcf(indexed_vcf, expected_build=BUILD, source_artifact=ARTIFACT)
    streamed, summary = _stream(indexed_vcf)

    assert summary.strategy == STRATEGY_INDEXED, "the fixture must exercise streaming"
    assert len(streamed) == len(batch.variants)
    for index, (left, right) in enumerate(zip(streamed, batch.variants, strict=True)):
        assert left == right, f"record {index} differs between the streaming and batch paths"
    assert streamed == list(batch.variants)


@pytest.mark.unit
def test_streaming_reports_the_same_warnings_and_skips(indexed_vcf: Path) -> None:
    """A discarded variant that stops being counted is invisible and unrecoverable."""
    batch = read_vcf(indexed_vcf, expected_build=BUILD, source_artifact=ARTIFACT)
    _, summary = _stream(indexed_vcf)

    assert summary.warnings == batch.warnings
    assert summary.skipped_reasons == batch.skipped_reasons
    assert summary.skipped_count == batch.skipped_count
    assert summary.declared_build == batch.declared_build
    assert summary.sample_id == batch.sample_id
    assert summary.complete is True
    # The fixture is only worth anything if it actually discards things.
    assert summary.skipped_count > 0
    assert len(summary.skipped_reasons) >= 4


@pytest.mark.unit
def test_non_canonical_contig_block_is_still_visited_and_counted(indexed_vcf: Path) -> None:
    """Streaming skips a contig's records, not the contig itself."""
    _, summary = _stream(indexed_vcf)
    # Three ALT alleles across the two chrUn lines, all counted.
    assert "non_canonical_contig (n=3)" in summary.skipped_reasons


@pytest.mark.unit
def test_source_line_index_matches_the_sequential_read(indexed_vcf: Path) -> None:
    """The per-contig read must reproduce the file's own data-line ordinal.

    ``source_line_index`` is a provenance claim about which line a record came
    from. A streaming reader that renumbered records in visit order would emit a
    field that quietly means something different.
    """
    batch = read_vcf(indexed_vcf, expected_build=BUILD, source_artifact=ARTIFACT)
    streamed, _ = _stream(indexed_vcf)
    assert [r.source_line_index for r in streamed] == [r.source_line_index for r in batch.variants]


@pytest.mark.unit
def test_both_backends_stream_identically(tmp_path: Path) -> None:
    """Backend equivalence has to survive the streaming path too, not only read_vcf.

    Run against ``_WELL_FORMED_ROWS``: htslib and the text parser agree on every
    record a VCF is allowed to contain, and disagree about what to do with one it
    is not (see ``test_the_backends_diverge_on_a_line_htslib_refuses``). That
    disagreement is exactly why the resolved backend is recorded as provenance.
    """
    pytest.importorskip("cyvcf2")
    path = _write_indexed(tmp_path, _WELL_FORMED_ROWS, name="wellformed")
    text, text_summary = _stream(path, BACKEND_TEXT)
    native, native_summary = _stream(path, BACKEND_CYVCF2)

    assert native_summary.strategy == STRATEGY_INDEXED
    assert native == text
    assert replace(native_summary, backend=BACKEND_TEXT) == text_summary


@pytest.mark.unit
def test_the_backends_diverge_on_a_line_htslib_refuses(indexed_vcf: Path) -> None:
    """Pre-existing, unchanged, and the reason `backend` is on the result.

    The fixture holds a line with no FORMAT/sample columns. The text parser reads
    it as a record with no genotype and warns; htslib raises out of
    ``variant.genotypes`` and the file is unreadable through that backend. Nothing
    here decides which is right — this test exists so the difference is a recorded
    fact rather than a surprise in a WGS run.
    """
    pytest.importorskip("cyvcf2")
    streamed, summary = _stream(indexed_vcf, BACKEND_TEXT)
    assert summary.backend == BACKEND_TEXT
    assert any("missing_format_gt" in warning for warning in summary.warnings)
    assert any(r.coordinate.contig == "chr1" and r.coordinate.position == 6000 for r in streamed)
    with pytest.raises(Exception, match="genotypes"):
        _stream(indexed_vcf, BACKEND_CYVCF2)


# ---------------------------------------------------------------------------
# 2. Ordering is fixed here, not inherited
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_a_multiallelic_site_is_emitted_in_sorted_allele_order(indexed_vcf: Path) -> None:
    """`C>T,A` at chr1:2000 must emit A before T, as the global sort would."""
    streamed, _ = _stream(indexed_vcf)
    at_2000 = [
        r for r in streamed if r.coordinate.contig == "chr1" and r.coordinate.position == 2000
    ]
    assert [(r.coordinate.ref, r.coordinate.alt) for r in at_2000] == [
        ("C", "A"),
        ("C", "T"),
        ("CAT", "C"),
    ]


@pytest.mark.unit
def test_contigs_come_out_in_karyotype_order_not_lexical_order(indexed_vcf: Path) -> None:
    streamed, _ = _stream(indexed_vcf)
    seen: list[str] = []
    for record in streamed:
        if not seen or seen[-1] != record.coordinate.contig:
            seen.append(record.coordinate.contig)
    assert seen == ["chr1", "chr2", "chr10", "chrX"]
    # ...and that order is the one CANONICAL_CONTIGS declares, not the file's and
    # not a lexical sort of it.
    assert seen == [c for c in CANONICAL_CONTIGS if c in set(seen)]
    assert seen != sorted(seen)


@pytest.mark.unit
def test_a_contig_shuffled_file_falls_back_rather_than_renumbering(tmp_path: Path) -> None:
    """chr2 before chr1 in the file: correct records, no streaming, no silent drift."""
    shuffled = (
        "chr2\t5000\t.\tC\tA\t540.2\tPASS\t.\tGT:DP:AD:GQ\t0/1:40:21,19:95",
        "chr1\t1000\t.\tA\tG\t820.5\tPASS\t.\tGT:DP:AD:GQ\t0/1:45:23,22:99",
    )
    path = _write_indexed(tmp_path, shuffled, name="shuffled")
    streamed, summary = _stream(path)

    assert summary.strategy == STRATEGY_BUFFERED
    batch = read_vcf(path, expected_build=BUILD, source_artifact=ARTIFACT)
    assert streamed == list(batch.variants)
    assert [r.coordinate.contig for r in streamed] == ["chr1", "chr2"]


@pytest.mark.unit
def test_an_unindexed_file_falls_back_to_the_buffered_path(tmp_path: Path) -> None:
    path = _write_plain(tmp_path / "plain.vcf")
    streamed, summary = _stream(path)
    assert summary.strategy == STRATEGY_BUFFERED
    batch = read_vcf(path, expected_build=BUILD, source_artifact=ARTIFACT)
    assert streamed == list(batch.variants)


@pytest.mark.unit
def test_repeat_streams_are_byte_identical_under_different_hash_seeds(indexed_vcf: Path) -> None:
    """GP-30, proved across processes rather than within one.

    Set ordering, dict ordering and `hash()` are the classic ways a streaming
    rewrite becomes seed-dependent, and none of them show up inside a single
    interpreter that always uses the same seed.
    """
    script = (
        "import sys;"
        "from dataclasses import asdict;"
        "from pathlib import Path;"
        "from mva.determinism import canonical_json;"
        "from mva.ingestion.reader import iter_vcf;"
        "from mva.models.genome import GenomeBuild;"
        "s = iter_vcf(Path(sys.argv[1]), expected_build=GenomeBuild.GRCH38,"
        " source_artifact='streaming-fixture');"
        "rows = [v.model_dump(mode='json') for v in s];"
        "print(canonical_json({'rows': rows, 'summary': asdict(s.summary())}))"
    )
    outputs = [
        subprocess.run(  # noqa: S603 - fixed argv, no shell
            [sys.executable, "-c", script, str(indexed_vcf)],
            check=True,
            capture_output=True,
            text=True,
            env={"PYTHONHASHSEED": seed, "PATH": "/usr/bin:/bin"},
            cwd=Path(__file__).resolve().parents[2],
        ).stdout
        for seed in ("0", "1", "12345")
    ]
    assert outputs[0] == outputs[1] == outputs[2]
    assert outputs[0].strip()


# ---------------------------------------------------------------------------
# 3. Chunking, and the composition it exists for
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_the_chunk_boundary_gap_clears_the_normalisation_shift_window() -> None:
    """The reason chunked normalisation is safe, asserted against the real constant.

    Left-alignment can move a record up to ``MAX_SHIFT_BP`` bases left. If a chunk
    boundary were narrower than that, one chunk's normalised output could overtake
    the next chunk's and the concatenation would not be sorted.
    """
    assert DEFAULT_CHUNK_BOUNDARY_GAP > MAX_SHIFT_BP


@pytest.mark.unit
def test_chunks_concatenate_back_into_the_same_stream(indexed_vcf: Path) -> None:
    handle = iter_vcf(indexed_vcf, expected_build=BUILD, source_artifact=ARTIFACT)
    batched = [
        record for chunk in handle.chunks(target_size=2, boundary_gap=100) for record in chunk
    ]
    batch = read_vcf(indexed_vcf, expected_build=BUILD, source_artifact=ARTIFACT)
    assert batched == list(batch.variants)
    assert handle.summary().skipped_count == batch.skipped_count


@pytest.mark.unit
def test_chunks_only_cut_at_a_safe_boundary(indexed_vcf: Path) -> None:
    handle = iter_vcf(indexed_vcf, expected_build=BUILD, source_artifact=ARTIFACT)
    chunks = list(handle.chunks(target_size=1, boundary_gap=500))
    assert len(chunks) > 1
    for previous, following in pairwise(chunks):
        last = previous[-1].coordinate
        first = following[0].coordinate
        assert last.contig != first.contig or first.position - last.position > 500


@pytest.mark.unit
def test_chunks_refuse_to_cut_unsafely_rather_than_mis_order(indexed_vcf: Path) -> None:
    handle = iter_vcf(indexed_vcf, expected_build=BUILD, source_artifact=ARTIFACT)
    with pytest.raises(IngestionError, match="safe chunk boundary"):
        list(handle.chunks(target_size=1, boundary_gap=10_000_000, max_size=2))


@pytest.mark.unit
def test_chunked_normalisation_equals_whole_callset_normalisation(indexed_vcf: Path) -> None:
    """The composition the orchestrator is asked to adopt, checked against the old one."""
    batch = read_vcf(indexed_vcf, expected_build=BUILD, source_artifact=ARTIFACT)
    expected = normalise_variants(batch.variants)

    handle = iter_vcf(indexed_vcf, expected_build=BUILD, source_artifact=ARTIFACT)
    streamed = [
        record
        for chunk in handle.chunks(target_size=2)
        for record in normalise_variants(chunk).variants
    ]
    assert streamed == list(expected.variants)


# ---------------------------------------------------------------------------
# 4. Streaming QC
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_streaming_qc_equals_batch_qc(indexed_vcf: Path) -> None:
    batch = read_vcf(indexed_vcf, expected_build=BUILD, source_artifact=ARTIFACT)
    normalised = normalise_variants(batch.variants)
    expected = assess_quality(normalised.variants, thresholds=THRESHOLDS, clock=CLOCK)

    stream = iter_assessed(normalised.variants, thresholds=THRESHOLDS, clock=CLOCK)
    variants = []
    evidence = []
    for item in stream:
        variants.append(item.variant)
        evidence.append(item.evidence)

    assert variants == list(expected.variants)
    assert evidence == list(expected.evidence)
    assert stream.metrics() == expected.metrics
    # The metrics block is written into an artifact, so equality of the dicts is
    # not enough: the serialised bytes must match too.
    assert canonical_json(stream.metrics()) == canonical_json(expected.metrics)


@pytest.mark.unit
def test_end_to_end_streaming_matches_end_to_end_batch(indexed_vcf: Path) -> None:
    """Read, normalise and QC — streamed and chunked against materialised and sorted."""
    batch = read_vcf(indexed_vcf, expected_build=BUILD, source_artifact=ARTIFACT)
    expected = assess_quality(
        normalise_variants(batch.variants).variants, thresholds=THRESHOLDS, clock=CLOCK
    )

    handle = iter_vcf(indexed_vcf, expected_build=BUILD, source_artifact=ARTIFACT)
    normalised = (
        record
        for chunk in handle.chunks(target_size=3)
        for record in normalise_variants(chunk).variants
    )
    stream = iter_assessed(normalised, thresholds=THRESHOLDS, clock=CLOCK)
    assessed = list(stream)

    assert [a.variant for a in assessed] == list(expected.variants)
    assert [a.evidence for a in assessed] == list(expected.evidence)
    assert stream.metrics() == expected.metrics
    assert handle.summary().skipped_count == batch.skipped_count


@pytest.mark.unit
def test_streaming_qc_refuses_unsorted_input() -> None:
    """It cannot sort, so it must not pretend the order was fine."""
    path = Path(__file__).resolve().parents[1] / "fixtures" / "synthetic" / "synthetic_case.vcf"
    batch = read_vcf(path, expected_build=BUILD, source_artifact=ARTIFACT)
    reversed_order = list(reversed(batch.variants))
    stream = iter_assessed(reversed_order, thresholds=THRESHOLDS, clock=CLOCK)
    with pytest.raises(IngestionError, match="out of coordinate order"):
        list(stream)


# ---------------------------------------------------------------------------
# 5. Fail-closed accounting and PRIV-09
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_summary_refuses_to_report_before_the_stream_is_drained(indexed_vcf: Path) -> None:
    handle = iter_vcf(indexed_vcf, expected_build=BUILD, source_artifact=ARTIFACT)
    iterator = iter(handle)
    next(iterator)
    with pytest.raises(IngestionError, match="not read to its end"):
        handle.summary()
    partial = handle.partial_summary()
    assert partial.complete is False


@pytest.mark.unit
def test_a_stream_refuses_a_second_iteration(indexed_vcf: Path) -> None:
    handle = iter_vcf(indexed_vcf, expected_build=BUILD, source_artifact=ARTIFACT)
    list(handle)
    with pytest.raises(IngestionError, match="already been iterated"):
        list(handle)


@pytest.mark.unit
def test_streaming_errors_never_echo_a_coordinate_or_a_genotype(indexed_vcf: Path) -> None:
    """PRIV-09. An exception message reaches the terminal, the log and an agent."""
    forbidden = (
        "chr1",
        "chr2",
        "chr10",
        "chrX",
        "SYNTH_PROBAND01",
        "0/1",
        "1|0",
        "./.",
        "CAT",
        "2000",
        "5000",
    )
    messages: list[str] = []

    handle = iter_vcf(indexed_vcf, expected_build=BUILD, source_artifact=ARTIFACT)
    with pytest.raises(IngestionError) as unsafe_cut:
        list(handle.chunks(target_size=1, boundary_gap=10_000_000, max_size=2))
    messages.append(str(unsafe_cut.value))

    batch = read_vcf(indexed_vcf, expected_build=BUILD, source_artifact=ARTIFACT)
    with pytest.raises(IngestionError) as unsorted:
        list(iter_assessed(list(reversed(batch.variants)), thresholds=THRESHOLDS, clock=CLOCK))
    messages.append(str(unsorted.value))

    drained = iter_vcf(indexed_vcf, expected_build=BUILD, source_artifact=ARTIFACT)
    with pytest.raises(IngestionError) as premature:
        drained.summary()
    messages.append(str(premature.value))

    for message in messages:
        for token in forbidden:
            assert token not in message, f"{token!r} leaked into: {message}"


@pytest.mark.unit
def test_a_stream_that_is_abandoned_says_its_counts_are_partial(indexed_vcf: Path) -> None:
    handle = iter_vcf(indexed_vcf, expected_build=BUILD, source_artifact=ARTIFACT)
    for _ in zip(handle, range(2), strict=False):
        pass
    summary = handle.partial_summary()
    assert summary.complete is False
    assert summary.strategy == STRATEGY_INDEXED
