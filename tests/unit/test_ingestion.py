"""Unit tests for the ingestion stage (reader, normalisation, QC).

The invariants under test are the ones that would fail *silently* in production:
a coerced genome build, a multiallelic site whose allelic depths were assigned to
the wrong allele, a low-quality call quietly dropped instead of flagged, a
left-alignment claimed but never performed, and a warning string that leaks the
record it was complaining about.
"""

from __future__ import annotations

from collections.abc import Sequence
from datetime import UTC, datetime
from pathlib import Path

import pytest
from hypothesis import assume, given, settings
from hypothesis import strategies as st

import mva.ingestion.normalise as normalise_module
import mva.ingestion.qc as qc_module
import mva.ingestion.reader as reader_module
from mva.clock import FixedClock
from mva.config import QualityThresholds
from mva.errors import GenomeBuildMismatchError, ReferenceMismatchError
from mva.ingestion import (
    BACKEND_CYVCF2,
    BACKEND_TEXT,
    FLAG_FILTERED_BY_CALLER,
    FLAG_HIGH_ALLELE_BALANCE,
    FLAG_LOW_ALLELE_BALANCE,
    FLAG_LOW_DEPTH,
    FLAG_LOW_GQ,
    FLAG_POSSIBLE_MOSAIC,
    OP_LEFT_ALIGN,
    OP_SPLIT_MULTIALLELIC,
    OP_TRIM,
    REF_ALLELE_MISMATCH_FLAG,
    SUPPORTED_BACKENDS,
    IngestionResult,
    allele_fraction,
    assess_quality,
    detect_backend,
    normalise_variants,
    read_vcf,
    split_multiallelic,
    trim_and_left_align,
)
from mva.models.evidence import (
    AssertionTier,
    EvidenceCategory,
    EvidenceDirection,
    EvidenceType,
)
from mva.models.genome import GenomeBuild, GenomicCoordinate
from mva.models.variant import FilterStatus, Genotype, VariantRecord, Zygosity

FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "synthetic" / "synthetic_case.vcf"
ARTIFACT = "tests/fixtures/synthetic/synthetic_case.vcf"
CLOCK = FixedClock(datetime(2026, 1, 1, tzinfo=UTC))
THRESHOLDS = QualityThresholds()

_FORMAT_HEADERS = (
    '##FILTER=<ID=PASS,Description="All filters passed">',
    '##FILTER=<ID=LowQual,Description="Low quality variant call">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
)
_COLUMNS = (
    "#CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
    "FORMAT",
    "SYNTH_PROBAND01",
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _row(*fields: str) -> str:
    return "\t".join(fields)


def _write_vcf(path: Path, rows: Sequence[str], *, reference: str | None = "GRCh38") -> Path:
    lines = ["##fileformat=VCFv4.2"]
    if reference is not None:
        lines.append(f"##reference={reference}")
    lines.extend(_FORMAT_HEADERS)
    lines.append(_row(*_COLUMNS))
    lines.extend(rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _variant(
    *,
    contig: str = "chr1",
    position: int = 100,
    ref: str = "C",
    alt: str = "T",
    zygosity: Zygosity = Zygosity.HET,
    genotype_string: str = "0/1",
    phased: bool = False,
    phase_set: int | None = None,
    depth: int | None = 40,
    ref_reads: int | None = 20,
    alt_reads: int | None = 20,
    genotype_quality: int | None = 99,
    filter_status: FilterStatus = FilterStatus.PASS,
    normalisation_ops: tuple[str, ...] = (),
) -> VariantRecord:
    return VariantRecord(
        coordinate=GenomicCoordinate(
            build=GenomeBuild.GRCH38, contig=contig, position=position, ref=ref, alt=alt
        ),
        genotype=Genotype(
            zygosity=zygosity,
            genotype_string=genotype_string,
            phased=phased,
            phase_set=phase_set,
            depth=depth,
            ref_reads=ref_reads,
            alt_reads=alt_reads,
            genotype_quality=genotype_quality,
        ),
        filter_status=filter_status,
        raw_filters=("PASS",) if filter_status is FilterStatus.PASS else ("LowQual",),
        source_artifact="unit-test",
        source_line_index=0,
        normalisation_ops=normalisation_ops,
    )


class _FakeReference:
    """A tiny in-memory ``ReferenceLookup`` over one contig.

    ``fetch`` is 1-based inclusive, matching the Protocol; out-of-range requests
    raise, so the "lookup failed" path is exercised as well as the happy path.
    """

    def __init__(self, contig: str, sequence: str, *, origin: int = 1) -> None:
        self._contig = contig
        self._sequence = sequence
        self._origin = origin

    def fetch(self, contig: str, start: int, end: int) -> str:
        if contig != self._contig:
            msg = f"unknown contig {contig}"
            raise KeyError(msg)
        lo = start - self._origin
        hi = end - self._origin + 1
        if lo < 0 or hi > len(self._sequence):
            msg = "requested range is outside the loaded sequence"
            raise IndexError(msg)
        return self._sequence[lo:hi]


def _read_fixture(backend: str = BACKEND_TEXT) -> IngestionResult:
    return read_vcf(
        FIXTURE,
        expected_build=GenomeBuild.GRCH38,
        source_artifact=ARTIFACT,
        backend=backend,
    )


def _find(variants: Sequence[VariantRecord], contig: str, position: int) -> VariantRecord:
    matches = [
        v for v in variants if v.coordinate.contig == contig and v.coordinate.position == position
    ]
    assert matches, f"no record at {contig}:{position}"
    return matches[0]


# ---------------------------------------------------------------------------
# 1. Multiallelic decomposition
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_multiallelic_site_splits_into_two_het_records() -> None:
    """chr11:5000000 A>G,GT with GT=1/2 becomes two het records, AD split per allele."""
    result = _read_fixture()
    chr11 = [v for v in result.variants if v.coordinate.contig == "chr11"]

    assert len(chr11) == 2
    by_alt = {v.coordinate.alt: v for v in chr11}
    assert set(by_alt) == {"G", "GT"}

    for record in chr11:
        assert record.genotype.zygosity is Zygosity.HET
        assert record.normalisation_ops == (OP_SPLIT_MULTIALLELIC,)
        # The original site call is preserved verbatim for auditability.
        assert record.genotype.genotype_string == "1/2"
        assert record.genotype.ref_reads == 2
        assert record.genotype.depth == 44

    # AD is [ref, alt1, alt2] -> allele 1 takes AD[1], allele 2 takes AD[2].
    assert by_alt["G"].genotype.alt_reads == 21
    assert by_alt["GT"].genotype.alt_reads == 21


@pytest.mark.unit
def test_multiallelic_allelic_depths_are_assigned_per_allele(tmp_path: Path) -> None:
    """With distinct per-allele depths, each record must take its own AD slot."""
    path = _write_vcf(
        tmp_path / "multi.vcf",
        [
            _row(
                "chr11",
                "5000000",
                ".",
                "A",
                "G,GT",
                "455.0",
                "PASS",
                ".",
                "GT:DP:AD:GQ",
                "1/2:44:2,30,7",
            )
        ],
    )
    result = read_vcf(
        path, expected_build=GenomeBuild.GRCH38, source_artifact="tmp", backend=BACKEND_TEXT
    )
    by_alt = {v.coordinate.alt: v for v in result.variants}

    assert by_alt["G"].genotype.ref_reads == 2
    assert by_alt["G"].genotype.alt_reads == 30
    assert by_alt["GT"].genotype.ref_reads == 2
    assert by_alt["GT"].genotype.alt_reads == 7
    assert by_alt["G"].genotype.allele_balance == pytest.approx(30 / 32)

    normalised = normalise_variants(result.variants)
    assert normalised.operations_applied[OP_SPLIT_MULTIALLELIC] == 2


@pytest.mark.unit
def test_split_multiallelic_call_is_not_flagged_for_allele_balance() -> None:
    """After splitting, ref_reads excludes the other ALT's reads.

    Banding ``alt/(ref+alt)`` would read 21/23 = 0.91 and flag a textbook compound
    het as ``high_allele_balance``. The site allele fraction ``alt/DP`` is banded
    instead, so the candidate is not down-ranked for an artifact of representation.
    """
    ingested = _read_fixture()
    chr11 = [v for v in ingested.variants if v.coordinate.contig == "chr11"]

    for record in chr11:
        assert record.genotype.allele_balance == pytest.approx(21 / 23)
        assert allele_fraction(record) == pytest.approx(21 / 44)

    qc = assess_quality(ingested.variants, thresholds=THRESHOLDS, clock=CLOCK)
    for record in [v for v in qc.variants if v.coordinate.contig == "chr11"]:
        assert FLAG_HIGH_ALLELE_BALANCE not in record.qc_flags
        assert record.qc_flags == ()


@pytest.mark.unit
def test_split_multiallelic_still_flags_a_genuinely_skewed_allele(tmp_path: Path) -> None:
    path = _write_vcf(
        tmp_path / "skewed.vcf",
        [
            _row(
                "chr11",
                "5000000",
                ".",
                "A",
                "G,GT",
                "455.0",
                "PASS",
                ".",
                "GT:DP:AD:GQ",
                "1/2:40:2,36,2",
            )
        ],
    )
    ingested = read_vcf(
        path, expected_build=GenomeBuild.GRCH38, source_artifact="tmp", backend=BACKEND_TEXT
    )
    qc = assess_quality(ingested.variants, thresholds=THRESHOLDS, clock=CLOCK)
    by_alt = {v.coordinate.alt: v for v in qc.variants}

    assert FLAG_HIGH_ALLELE_BALANCE in by_alt["G"].qc_flags  # 36/40 = 0.90
    assert FLAG_LOW_ALLELE_BALANCE in by_alt["GT"].qc_flags  # 2/40 = 0.05


@pytest.mark.unit
def test_split_multiallelic_is_idempotent() -> None:
    """Re-splitting an already-decomposed record must not change or duplicate it."""
    record = _variant(
        genotype_string="1/2",
        zygosity=Zygosity.HET,
        normalisation_ops=(OP_SPLIT_MULTIALLELIC,),
    )
    assert split_multiallelic(record) == (record,)

    unresolved = _variant(genotype_string="1/2", zygosity=Zygosity.UNKNOWN)
    (resolved,) = split_multiallelic(unresolved)
    assert resolved.genotype.zygosity is Zygosity.HET
    assert resolved.normalisation_ops == (OP_SPLIT_MULTIALLELIC,)
    assert split_multiallelic(resolved) == (resolved,)


# ---------------------------------------------------------------------------
# 2. Genome build
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_declared_build_conflict_raises_and_names_both_builds(tmp_path: Path) -> None:
    path = _write_vcf(
        tmp_path / "b37.vcf",
        [_row("chr15", "40200000", ".", "C", "T", "820.5", "PASS", ".", "GT:DP", "0/1:45")],
        reference="GRCh37",
    )
    with pytest.raises(GenomeBuildMismatchError) as excinfo:
        read_vcf(
            path,
            expected_build=GenomeBuild.GRCH38,
            source_artifact="tmp",
            backend=BACKEND_TEXT,
        )
    message = str(excinfo.value)
    assert "GRCh37" in message
    assert "GRCh38" in message


@pytest.mark.unit
def test_absent_build_header_warns_but_proceeds(tmp_path: Path) -> None:
    path = _write_vcf(
        tmp_path / "nobuild.vcf",
        [_row("chr15", "40200000", ".", "C", "T", "820.5", "PASS", ".", "GT:DP", "0/1:45")],
        reference=None,
    )
    result = read_vcf(
        path, expected_build=GenomeBuild.GRCH38, source_artifact="tmp", backend=BACKEND_TEXT
    )
    assert result.declared_build is None
    assert len(result.variants) == 1
    assert result.variants[0].build is GenomeBuild.GRCH38
    assert any("genome_build_not_declared_in_header" in w for w in result.warnings)


@pytest.mark.unit
def test_fixture_declares_grch38() -> None:
    assert _read_fixture().declared_build is GenomeBuild.GRCH38


# ---------------------------------------------------------------------------
# 3. Reference-allele validation
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_reference_mismatch_flags_and_contradicts_by_default() -> None:
    record = _variant(contig="chr1", position=5, ref="C", alt="T")
    reference = _FakeReference("chr1", "TCGAAAAGCTGATCGATCGA")  # position 5 is 'A', not 'C'

    result = normalise_variants([record], reference=reference)

    assert len(result.variants) == 1, "a REF mismatch is flagged, never dropped (GP-13)"
    assert REF_ALLELE_MISMATCH_FLAG in result.variants[0].qc_flags
    assert any(REF_ALLELE_MISMATCH_FLAG in warning for warning in result.warnings)

    qc = assess_quality(result.variants, thresholds=THRESHOLDS, clock=CLOCK)
    item = qc.evidence[0]
    assert item.direction is EvidenceDirection.CONTRADICTS
    assert item.is_contradiction
    assert item.category is EvidenceCategory.ANALYTICAL
    assert "reference" in item.limitations.lower()


@pytest.mark.unit
def test_reference_mismatch_raises_only_in_strict_mode() -> None:
    record = _variant(contig="chr1", position=5, ref="C", alt="T")
    reference = _FakeReference("chr1", "TCGAAAAGCTGATCGATCGA")

    with pytest.raises(ReferenceMismatchError):
        normalise_variants([record], reference=reference, strict_reference=True)


@pytest.mark.unit
def test_matching_reference_raises_no_flag() -> None:
    record = _variant(contig="chr1", position=3, ref="G", alt="T")
    reference = _FakeReference("chr1", "TCGAAAAGCTGATCGATCGA")  # position 3 is 'G'

    result = normalise_variants([record], reference=reference, strict_reference=True)
    assert result.variants[0].qc_flags == ()


@pytest.mark.unit
def test_failed_reference_lookup_is_not_treated_as_a_mismatch() -> None:
    """GP-14: a lookup that could not answer is absence of information."""
    record = _variant(contig="chr2", position=5, ref="C", alt="T")
    reference = _FakeReference("chr1", "TCGAAAAGCTGATCGATCGA")

    result = normalise_variants([record], reference=reference, strict_reference=True)
    assert result.variants[0].qc_flags == ()
    assert any("reference_lookup_failed" in warning for warning in result.warnings)


# ---------------------------------------------------------------------------
# 4. QC flags a low-quality call without deleting it
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_low_quality_call_is_flagged_and_retained() -> None:
    ingested = _read_fixture()
    qc = assess_quality(ingested.variants, thresholds=THRESHOLDS, clock=CLOCK)

    assert len(qc.variants) == len(ingested.variants), "QC never removes a record (GP-13)"
    record = _find(qc.variants, "chr15", 40211000)

    assert record.genotype.depth == 4
    assert record.genotype.genotype_quality == 8
    assert record.filter_status is FilterStatus.FILTERED
    assert FLAG_LOW_DEPTH in record.qc_flags
    assert FLAG_LOW_GQ in record.qc_flags
    assert FLAG_FILTERED_BY_CALLER in record.qc_flags

    item = next(e for e in qc.evidence if e.subject_id == record.variant_id)
    assert item.direction is EvidenceDirection.CONTRADICTS
    assert item.tier is AssertionTier.OBSERVED_DATA
    assert item.evidence_type is EvidenceType.DIRECT_MEASUREMENT
    assert len(item.limitations) > 30
    assert qc.metrics["flag_low_depth"] >= 1


@pytest.mark.unit
def test_every_variant_gets_exactly_one_analytical_evidence_item() -> None:
    ingested = _read_fixture()
    qc = assess_quality(ingested.variants, thresholds=THRESHOLDS, clock=CLOCK)

    assert len(qc.evidence) == len(qc.variants)
    assert [e.subject_id for e in qc.evidence] == [v.variant_id for v in qc.variants]
    for item in qc.evidence:
        assert item.category is EvidenceCategory.ANALYTICAL
        assert item.tier is AssertionTier.OBSERVED_DATA
        assert item.evidence_type is EvidenceType.DIRECT_MEASUREMENT
        assert item.timestamp == CLOCK.now()
        assert item.limitations.strip()


# ---------------------------------------------------------------------------
# 5. Mosaic precedence
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_possible_mosaic_beats_low_allele_balance() -> None:
    """AB below the het band but at/above the mosaic floor is signal, not noise."""
    mosaic = _variant(position=100, ref_reads=17, alt_reads=3)  # AB = 0.15
    very_low = _variant(position=200, ref_reads=19, alt_reads=1)  # AB = 0.05
    high = _variant(position=300, ref_reads=2, alt_reads=38)  # AB = 0.95

    qc = assess_quality([mosaic, very_low, high], thresholds=THRESHOLDS, clock=CLOCK)
    flagged = {v.coordinate.position: v.qc_flags for v in qc.variants}

    assert FLAG_POSSIBLE_MOSAIC in flagged[100]
    assert FLAG_LOW_ALLELE_BALANCE not in flagged[100]
    assert FLAG_LOW_ALLELE_BALANCE in flagged[200]
    assert FLAG_POSSIBLE_MOSAIC not in flagged[200]
    assert FLAG_HIGH_ALLELE_BALANCE in flagged[300]

    mosaic_evidence = next(e for e in qc.evidence if ":100:" in e.subject_id)
    assert mosaic_evidence.direction is EvidenceDirection.NEUTRAL
    assert "mosaic" in mosaic_evidence.limitations.lower()


@pytest.mark.unit
def test_allele_balance_flags_apply_only_to_heterozygous_calls() -> None:
    hom_alt = _variant(zygosity=Zygosity.HOM_ALT, genotype_string="1/1", ref_reads=1, alt_reads=47)
    qc = assess_quality([hom_alt], thresholds=THRESHOLDS, clock=CLOCK)
    assert FLAG_HIGH_ALLELE_BALANCE not in qc.variants[0].qc_flags


# ---------------------------------------------------------------------------
# 6. Trimming and left-alignment
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_anchored_insertion_is_already_minimal() -> None:
    record = _variant(position=40211000, ref="G", alt="GAT")
    result = trim_and_left_align(record, None)

    assert result.coordinate.position == 40211000
    assert result.coordinate.ref == "G"
    assert result.coordinate.alt == "GAT"
    assert OP_TRIM not in result.normalisation_ops


@pytest.mark.unit
def test_shared_prefix_is_trimmed() -> None:
    record = _variant(position=100, ref="AAT", alt="AAG")
    result = trim_and_left_align(record, None)

    assert (result.coordinate.position, result.coordinate.ref, result.coordinate.alt) == (
        102,
        "T",
        "G",
    )
    assert OP_TRIM in result.normalisation_ops
    assert OP_LEFT_ALIGN not in result.normalisation_ops


@pytest.mark.unit
def test_shared_suffix_and_prefix_are_both_trimmed() -> None:
    record = _variant(position=100, ref="GATTACA", alt="GATTA")
    result = trim_and_left_align(record, None)

    assert (result.coordinate.position, result.coordinate.ref, result.coordinate.alt) == (
        103,
        "TAC",
        "T",
    )


@pytest.mark.unit
def test_left_alignment_requires_a_reference() -> None:
    """Without a reference the record is trimmed but never *claimed* to be aligned."""
    record = _variant(position=6, ref="AA", alt="A")

    without = normalise_variants([record])
    assert OP_LEFT_ALIGN not in without.variants[0].normalisation_ops
    assert without.variants[0].coordinate.position == 6
    assert any("left_alignment_skipped" in warning for warning in without.warnings)
    assert OP_LEFT_ALIGN not in without.operations_applied

    reference = _FakeReference("chr1", "TCGAAAAGCTGATCGATCGA")
    with_reference = normalise_variants([record], reference=reference)
    shifted = with_reference.variants[0]
    assert (shifted.coordinate.position, shifted.coordinate.ref, shifted.coordinate.alt) == (
        3,
        "GA",
        "G",
    )
    assert OP_LEFT_ALIGN in shifted.normalisation_ops
    assert with_reference.operations_applied[OP_LEFT_ALIGN] == 1


# ---------------------------------------------------------------------------
# 7. Backend equivalence
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_detect_backend_returns_a_supported_token() -> None:
    backend = detect_backend()
    assert backend in SUPPORTED_BACKENDS
    assert backend != "auto"


@pytest.mark.unit
def test_backends_produce_identical_records_for_the_fixture() -> None:
    pytest.importorskip("cyvcf2")
    text = _read_fixture(BACKEND_TEXT)
    native = _read_fixture(BACKEND_CYVCF2)

    assert native.variants == text.variants
    assert native == text


@pytest.mark.unit
def test_auto_backend_matches_the_text_backend() -> None:
    assert read_vcf(
        FIXTURE, expected_build=GenomeBuild.GRCH38, source_artifact=ARTIFACT
    ) == _read_fixture(BACKEND_TEXT)


# ---------------------------------------------------------------------------
# 8. Phase
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_phased_record_preserves_phase_and_phase_set() -> None:
    result = _read_fixture()
    record = _find(result.variants, "chr15", 40201000)

    assert record.genotype.phased is True
    assert record.genotype.genotype_string == "1|0"
    assert record.genotype.phase_set == 40201000
    assert record.genotype.zygosity is Zygosity.HET

    unphased = _find(result.variants, "chr15", 40200000)
    assert unphased.genotype.phased is False
    assert unphased.genotype.phase_set is None


# ---------------------------------------------------------------------------
# 9. Determinism, ordering and hard filters
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_reading_twice_yields_equal_tuples() -> None:
    first = _read_fixture()
    second = _read_fixture()
    assert first.variants == second.variants
    assert first == second


@pytest.mark.unit
def test_full_stage_is_deterministic_across_repeat_runs() -> None:
    def run() -> tuple[tuple[str, ...], tuple[str, ...]]:
        ingested = _read_fixture()
        normalised = normalise_variants(ingested.variants)
        qc = assess_quality(normalised.variants, thresholds=THRESHOLDS, clock=CLOCK)
        return (
            tuple(v.variant_id for v in qc.variants),
            tuple(e.evidence_id for e in qc.evidence),
        )

    assert run() == run()


@pytest.mark.unit
def test_output_is_sorted_by_coordinate_sort_key() -> None:
    variants = _read_fixture().variants
    assert list(variants) == sorted(variants, key=lambda v: v.sort_key())
    assert [v.coordinate.contig for v in variants[:4]] == ["chr3", "chr7", "chr7", "chr11"]


@pytest.mark.unit
def test_non_canonical_contig_is_a_hard_filter(tmp_path: Path) -> None:
    path = _write_vcf(
        tmp_path / "alt.vcf",
        [
            _row("chr15", "40200000", ".", "C", "T", "820.5", "PASS", ".", "GT:DP", "0/1:45"),
            _row(
                "chr1_KI270766v1_alt",
                "1000",
                ".",
                "C",
                "T",
                "500.0",
                "PASS",
                ".",
                "GT:DP",
                "0/1:30",
            ),
        ],
    )
    result = read_vcf(
        path, expected_build=GenomeBuild.GRCH38, source_artifact="tmp", backend=BACKEND_TEXT
    )
    assert len(result.variants) == 1
    assert result.skipped_count == 1
    assert any("non_canonical_contig" in reason for reason in result.skipped_reasons)


@pytest.mark.unit
def test_symbolic_and_spanning_alleles_are_hard_filtered(tmp_path: Path) -> None:
    path = _write_vcf(
        tmp_path / "symbolic.vcf",
        [
            _row("chr15", "1000", ".", "C", "<DEL>", "500.0", "PASS", ".", "GT:DP", "0/1:30"),
            _row("chr15", "2000", ".", "C", "*", "500.0", "PASS", ".", "GT:DP", "0/1:30"),
            _row("chr15", "3000", ".", "C", ".", "500.0", "PASS", ".", "GT:DP", "0/0:30"),
        ],
    )
    result = read_vcf(
        path, expected_build=GenomeBuild.GRCH38, source_artifact="tmp", backend=BACKEND_TEXT
    )
    assert result.variants == ()
    assert result.skipped_count == 3
    codes = " ".join(result.skipped_reasons)
    assert "symbolic_or_structural_allele" in codes
    assert "spanning_deletion_allele" in codes
    assert "missing_alt_allele" in codes


@pytest.mark.unit
def test_warnings_metrics_and_skip_reasons_never_echo_record_text() -> None:
    """PRIV-09: reason codes and counts only — no line, sample, genotype or locus."""
    ingested = _read_fixture()
    normalised = normalise_variants(ingested.variants)
    qc = assess_quality(normalised.variants, thresholds=THRESHOLDS, clock=CLOCK)

    leaks = ("SYNTH_PROBAND01", "0/1", "1|0", "1/2", "chr15", "40200000", "\t")
    strings = [*ingested.warnings, *ingested.skipped_reasons, *normalised.warnings]
    for text in strings:
        for leak in leaks:
            assert leak not in text, f"privacy leak {leak!r} in {text!r}"

    for key, value in qc.metrics.items():
        assert isinstance(value, int | float)
        for leak in leaks:
            assert leak not in key


@pytest.mark.unit
def test_ingestion_module_imports_no_network_clients() -> None:
    """PRIV-05, checked at the module level as well as structurally."""
    for module in (reader_module, normalise_module, qc_module):
        source = Path(str(module.__file__)).read_text(encoding="utf-8")
        for forbidden in ("import requests", "import httpx", "import urllib", "import aiohttp"):
            assert forbidden not in source


# ---------------------------------------------------------------------------
# 10. Property-based normalisation invariants
# ---------------------------------------------------------------------------


@pytest.mark.unit
@given(
    ref=st.text(alphabet="ACGT", min_size=1, max_size=8),
    alt=st.text(alphabet="ACGT", min_size=1, max_size=8),
    position=st.integers(min_value=1, max_value=10_000),
)
@settings(max_examples=300, deadline=None)
def test_trimming_is_idempotent_and_never_empties_an_allele(
    ref: str, alt: str, position: int
) -> None:
    assume(ref != alt)
    record = _variant(position=position, ref=ref, alt=alt)

    once = trim_and_left_align(record, None)
    twice = trim_and_left_align(once, None)

    assert once.coordinate == twice.coordinate
    assert once.normalisation_ops == twice.normalisation_ops

    assert len(once.coordinate.ref) >= 1
    assert len(once.coordinate.alt) >= 1
    assert once.coordinate.ref != once.coordinate.alt
    assert once.coordinate.position >= position
    # Trimming only removes shared bases; it can never lengthen an allele.
    assert len(once.coordinate.ref) <= len(ref)
    assert len(once.coordinate.alt) <= len(alt)
    # The reference span may only shrink from the right by what was trimmed.
    assert once.coordinate.end <= record.coordinate.end
    assert OP_LEFT_ALIGN not in once.normalisation_ops
