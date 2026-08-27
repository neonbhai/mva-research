"""One representation rule, two callers, and the state where it cannot be applied.

These tests exist because of a specific, verified defect. Variant representation
did not agree between ingestion and the annotation adapters, so equivalent
variants failed to join — and a failed join is not an error. It is silence, which
this pipeline reads as "no ClinVar record" and "no population frequency", which is
scored as novel and ultra-rare: the strongest promoting signal the ranker has. The
same bug therefore manufactures false positives and deletes true pathogenic
assertions at once.

Three things are pinned here:

1. **The rule is shared.** :mod:`mva.alleles` holds it; ingestion and the ClinVar
   adapter both call that one function. The property test below generates allele
   pairs and asserts the two callers can never disagree, because two
   implementations that agree today is exactly the state the repository was in
   when the defect was introduced.
2. **Left-alignment against a real indexed FASTA works**, including the boundary
   translations that are easy to get silently wrong: 1-based inclusive versus
   pysam's 0-based half-open, and a ``chr``-prefixed FASTA against a bare-contig
   VCF.
3. **The absence of a reference is a typed, surfaced state** (GP-14), not a skipped
   step. A run that could not left-align its indels has to say so.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import pysam
import pytest
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

from mva.alleles import (
    LeftAlignmentStatus,
    canonicalise_allele,
    rightmost_equivalent_position,
    trim_parsimoniously,
)
from mva.annotation.clinvar_vcf import ClinvarVcfAdapter
from mva.determinism import hash_file, stable_hash
from mva.errors import AdapterUnavailableError
from mva.ingestion.normalise import (
    FastaReference,
    normalise_variants,
    open_reference_fasta,
    trim_and_left_align,
)
from mva.models.genome import GenomeBuild, GenomicCoordinate
from mva.models.variant import (
    OP_LEFT_ALIGN,
    OP_TRIM,
    FilterStatus,
    Genotype,
    VariantRecord,
    Zygosity,
)

pytestmark = pytest.mark.unit

# ---------------------------------------------------------------------------
# A tiny reference with a repeat in it
# ---------------------------------------------------------------------------
#
# Positions are 1-based. The tract at 22-26 is five C's preceded by an A, which is
# the shape that makes a single-C insertion ambiguous: it can legally be written
# anywhere from POS 21 to POS 26 and every spelling is the same event. gnomAD and
# ClinVar store the left-most one. A caller that stores any other one does not
# join, and is told it has found a novel variant.
#
#          1234567890123456789 0 12345 678901234567890
CHR21_SEQ = "GATCGATCGATCGATCGATCACCCCCGTTAACCGGTTAACGGATCC"
REPEAT_TRACT_START = 21
REPEAT_TRACT_END = 26

#: The same insertion, spelled at each end of the tract. These two strings are the
#: whole problem in miniature.
LEFTMOST_INSERTION = (21, "A", "AC")
RIGHTMOST_INSERTION = (26, "C", "CC")


def write_fasta(directory: Path, sequences: dict[str, str], *, name: str = "ref.fa") -> Path:
    """A plain FASTA on disk, wrapped at 60 columns like every real one."""
    path = directory / name
    lines: list[str] = []
    for contig, sequence in sequences.items():
        lines.append(f">{contig}")
        lines.extend(sequence[index : index + 60] for index in range(0, len(sequence), 60))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


class InMemoryReference:
    """A ``ReferenceLookup`` with no file behind it, for the property test.

    1-based inclusive, like the Protocol requires. Deliberately not the FASTA
    reader: this test is about the rule, and mixing in htslib would let a reader
    wonder whether a disagreement came from the reference or from the algorithm.
    """

    def __init__(self, contig: str, sequence: str, *, start: int = 1) -> None:
        self._contig = contig
        self._sequence = sequence
        self._start = start

    def fetch(self, contig: str, start: int, end: int) -> str:
        if contig != self._contig:
            raise KeyError(contig)
        offset = start - self._start
        if offset < 0 or end - self._start >= len(self._sequence):
            raise IndexError(start)
        return self._sequence[offset : end - self._start + 1]


def make_record(
    position: int,
    ref: str,
    alt: str,
    *,
    contig: str = "chr21",
    ops: tuple[str, ...] = (),
) -> VariantRecord:
    return VariantRecord(
        coordinate=GenomicCoordinate(
            build=GenomeBuild.GRCH38, contig=contig, position=position, ref=ref, alt=alt
        ),
        genotype=Genotype(
            zygosity=Zygosity.HET,
            genotype_string="0/1",
            depth=40,
            ref_reads=20,
            alt_reads=20,
            genotype_quality=99,
        ),
        filter_status=FilterStatus.PASS,
        source_artifact="tests/unit/test_normalise_representation.py",
        normalisation_ops=ops,
    )


@pytest.fixture(scope="module")
def fasta_reference(tmp_path_factory: pytest.TempPathFactory) -> Iterator[FastaReference]:
    """An indexed FASTA whose contig is ``chr``-prefixed, as the real one is."""
    directory = tmp_path_factory.mktemp("reference")
    reference = open_reference_fasta(write_fasta(directory, {"chr21": CHR21_SEQ}))
    yield reference
    reference.close()


# ---------------------------------------------------------------------------
# 1. The FASTA boundary: coordinates and contig names
# ---------------------------------------------------------------------------


def test_open_reference_fasta_builds_its_own_index(tmp_path: Path) -> None:
    """No samtools binary exists on the target machine, so pysam must do it."""
    path = write_fasta(tmp_path, {"chr21": CHR21_SEQ})
    assert not path.with_name(path.name + ".fai").is_file()

    reference = open_reference_fasta(path)
    try:
        assert path.with_name(path.name + ".fai").is_file()
    finally:
        reference.close()


def test_fetch_is_one_based_inclusive(fasta_reference: FastaReference) -> None:
    """The single off-by-one that would shift every left-aligned indel by a base.

    pysam's ``fetch`` is 0-based half-open. If that translation is wrong, nothing
    raises: every indel simply lands one base away from where ClinVar and gnomAD
    hold it, and the whole run silently stops joining.
    """
    assert fasta_reference.fetch("chr21", 1, 1) == "G"
    assert fasta_reference.fetch("chr21", 1, 4) == "GATC"
    assert fasta_reference.fetch("chr21", REPEAT_TRACT_START, REPEAT_TRACT_START) == "A"
    assert fasta_reference.fetch("chr21", 22, 26) == "CCCCC"
    assert fasta_reference.fetch("chr21", len(CHR21_SEQ), len(CHR21_SEQ)) == CHR21_SEQ[-1]


def test_a_chr_prefixed_fasta_serves_a_bare_contig_vcf(tmp_path: Path) -> None:
    """The naming mismatch the real inputs have: FASTA ``chr21``, VCF ``21``.

    Ingestion canonicalises every contig to the UCSC spelling before a coordinate
    exists, so the lookup is always asked for ``chr21``. What varies is the FASTA.
    Both spellings must resolve, and a contig the FASTA does not hold must raise
    rather than silently answer with the wrong sequence.
    """
    bare = open_reference_fasta(write_fasta(tmp_path, {"21": CHR21_SEQ}, name="bare.fa"))
    try:
        assert bare.contig_map["chr21"] == "21"
        assert bare.fetch("chr21", 22, 26) == "CCCCC"
        assert "chr7" not in bare.contig_map
        with pytest.raises(KeyError):
            bare.fetch("chr7", 1, 1)
    finally:
        bare.close()


def test_the_block_cache_never_changes_a_base(tmp_path: Path) -> None:
    """The read cache is an optimisation, so it must be invisible in the output.

    Left-alignment walks one base at a time; a per-base htslib call over a
    whole-genome VCF is the difference between minutes and hours. A caching bug,
    though, produces confidently wrong reference bases and therefore confidently
    wrong coordinates, so every fetch is compared against pysam directly —
    including ranges that straddle and exceed the block boundary.
    """
    sequence = "".join("ACGT"[(index * 7 + index // 13) % 4] for index in range(20_000))
    path = write_fasta(tmp_path, {"chr21": sequence}, name="long.fa")
    reference = open_reference_fasta(path)
    raw = pysam.FastaFile(str(path))
    try:
        for start in (1, 4095, 4096, 4097, 8191, 8192, 8193, 12_000, 19_990):
            for length in (1, 2, 100, 4096, 5000):
                end = min(start + length - 1, len(sequence))
                assert reference.fetch("chr21", start, end) == raw.fetch("chr21", start - 1, end)
    finally:
        raw.close()
        reference.close()


def test_missing_reference_fasta_fails_loudly(tmp_path: Path) -> None:
    """Never degrade by accident: choosing to run without a reference is explicit."""
    with pytest.raises(AdapterUnavailableError):
        open_reference_fasta(tmp_path / "not-there.fa")


# ---------------------------------------------------------------------------
# 2. Left-alignment against a real reference
# ---------------------------------------------------------------------------


def test_reference_fasta_left_aligns_repeat_indel_before_frequency_join(
    fasta_reference: FastaReference,
) -> None:
    """A repeat-expansion insertion reaches the representation gnomAD holds.

    **This is the join that was failing.** The caller emits the insertion at the
    right-hand end of a five-base C tract; gnomAD and ClinVar store the left-most
    spelling. Trimming cannot reconcile those — undoing a right shift needs the
    bases to the left — so without a reference the variant is looked up under a key
    no population database has, comes back empty, and is scored as a novel,
    ultra-rare allele.
    """
    right_shifted = make_record(*RIGHTMOST_INSERTION)
    assert right_shifted.variant_id == "GRCh38:chr21:26:C:CC"

    result = normalise_variants([right_shifted], reference=fasta_reference)
    aligned = result.variants[0]

    position, ref, alt = LEFTMOST_INSERTION
    assert aligned.variant_id == f"GRCh38:chr21:{position}:{ref}:{alt}"
    assert result.left_alignment.status is LeftAlignmentStatus.APPLIED
    assert result.left_alignment.shifted_count == 1
    assert not result.left_alignment.is_degraded

    # Same event without a reference: still trimmed, still un-joinable, and the
    # result says which of those two it is.
    unaligned = normalise_variants([right_shifted])
    assert unaligned.variants[0].variant_id == "GRCh38:chr21:26:C:CC"
    assert unaligned.left_alignment.status is LeftAlignmentStatus.UNAVAILABLE_NO_REFERENCE


def test_every_spelling_in_the_repeat_tract_normalises_to_one_key(
    fasta_reference: FastaReference,
) -> None:
    """The insertion is legal at six positions; all six must produce one join key.

    Written as a sweep rather than a single case because the defect is not "one
    coordinate is wrong". It is that a family of equally valid spellings maps to a
    family of keys, and only one member of that family finds the database record.
    """
    spellings = [make_record(REPEAT_TRACT_START, "A", "AC")]
    spellings.extend(
        make_record(position, "C", "CC")
        for position in range(REPEAT_TRACT_START + 1, REPEAT_TRACT_END + 1)
    )
    keys = {
        normalise_variants([record], reference=fasta_reference).variants[0].variant_id
        for record in spellings
    }
    assert keys == {"GRCh38:chr21:21:A:AC"}


def test_a_moved_coordinate_always_leaves_a_trace(fasta_reference: FastaReference) -> None:
    """GP-31: left-alignment changes POS, and the submission is scored on POS.

    A coordinate that moved with no record of the move is unauditable. Both
    operations are recorded, and ``left_align`` is recorded *only* when a shift
    actually happened — an operation claimed but not performed is worse than none,
    because nothing downstream can detect it.
    """
    shifted = trim_and_left_align(make_record(*RIGHTMOST_INSERTION), fasta_reference)
    assert shifted.coordinate.position == REPEAT_TRACT_START
    assert OP_LEFT_ALIGN in shifted.normalisation_ops

    already_leftmost = trim_and_left_align(make_record(*LEFTMOST_INSERTION), fasta_reference)
    assert already_leftmost.coordinate.position == REPEAT_TRACT_START
    assert OP_LEFT_ALIGN not in already_leftmost.normalisation_ops

    trimmed_only = trim_and_left_align(make_record(21, "AC", "AG"), None)
    assert OP_TRIM in trimmed_only.normalisation_ops
    assert OP_LEFT_ALIGN not in trimmed_only.normalisation_ops


def test_left_alignment_is_idempotent(fasta_reference: FastaReference) -> None:
    """Re-normalising an aligned record must not move it or re-record the op."""
    once = trim_and_left_align(make_record(*RIGHTMOST_INSERTION), fasta_reference)
    twice = trim_and_left_align(once, fasta_reference)
    assert once.coordinate == twice.coordinate
    assert once.normalisation_ops == twice.normalisation_ops == (OP_LEFT_ALIGN,)


def test_the_shift_stops_at_the_start_of_a_contig(tmp_path: Path) -> None:
    """A tract running off the left edge must terminate, not walk into position 0."""
    reference = open_reference_fasta(write_fasta(tmp_path, {"chr21": "AAAAAAAAAACGT"}))
    try:
        aligned = trim_and_left_align(make_record(9, "A", "AA"), reference)
        assert aligned.coordinate.position == 1
        assert aligned.coordinate.ref == "A"
    finally:
        reference.close()


def test_rightmost_equivalent_position_bounds_the_repeat_tract(
    fasta_reference: FastaReference,
) -> None:
    """The mirror shift, which is what makes an index query complete.

    A source record spelling this insertion at the right-hand end of the tract
    occupies a span disjoint from a query at the left-hand end, so it is never
    fetched and the join fails exactly as if the source held nothing. Querying out
    to this bound closes that hole with a number read from the reference rather
    than a guessed padding constant.
    """
    position, ref, alt = LEFTMOST_INSERTION
    # One past the tract, not at its end: ``27 G>CG`` describes the same haplotype
    # with the padding base on the right. VCF convention puts the pad on the left,
    # so no source writes it that way — but a search bound that under-reaches
    # re-opens the miss, and one surplus base of index query costs nothing.
    assert (
        rightmost_equivalent_position(
            contig="chr21", position=position, ref=ref, alt=alt, reference=fasta_reference
        )
        == REPEAT_TRACT_END + 1
    )
    # A substitution is not in a tract and cannot move.
    assert (
        rightmost_equivalent_position(
            contig="chr21", position=5, ref="G", alt="T", reference=fasta_reference
        )
        == 5
    )


# ---------------------------------------------------------------------------
# 3. The degraded state (GP-14)
# ---------------------------------------------------------------------------


def test_a_missing_reference_produces_an_explicit_degraded_state() -> None:
    """Absence of a reference is a state a caller receives, not a step that is skipped.

    The previous behaviour emitted a warning code and moved on. A code in a list
    nobody branches on is indistinguishable from a clean run to every consumer of
    the result, which is how a whole-run indel-join failure stayed invisible.
    """
    result = normalise_variants([make_record(*RIGHTMOST_INSERTION), make_record(5, "G", "T")])

    report = result.left_alignment
    assert report.status is LeftAlignmentStatus.UNAVAILABLE_NO_REFERENCE
    assert report.is_degraded is True
    assert report.reference_available is False
    assert report.indel_count == 1  # the SNV is unaffected by left-alignment
    assert report.unaligned_indel_count == 1

    # It is legible without reading this code, and it names the consequence.
    described = report.describe()
    assert "DEGRADED" in described
    assert "reference_fasta" in described
    assert any(described == warning for warning in result.warnings)

    payload = report.as_dict()
    assert payload["status"] == "unavailable_no_reference"
    assert payload["degraded"] is True


def test_a_run_with_no_indels_is_not_reported_as_degraded() -> None:
    """GP-14 cuts both ways: 'nothing to align' must not read as 'could not align'."""
    result = normalise_variants([make_record(5, "G", "T")])
    assert result.left_alignment.status is LeftAlignmentStatus.NOT_REQUIRED
    assert result.left_alignment.is_degraded is False
    assert result.warnings == ()


def test_an_unreadable_reference_degrades_partially_rather_than_silently() -> None:
    """A reference that cannot answer is absence of information, not a match.

    The shift must stop and the record must be counted as un-aligned. Treating an
    unreadable base as agreement would move a coordinate on no evidence at all.
    """

    class BrokenReference:
        def fetch(self, contig: str, start: int, end: int) -> str:
            raise OSError("index unavailable")

    result = normalise_variants([make_record(*RIGHTMOST_INSERTION)], reference=BrokenReference())
    assert result.left_alignment.status is LeftAlignmentStatus.INCOMPLETE_REFERENCE_UNUSABLE
    assert result.left_alignment.is_degraded is True
    assert result.variants[0].coordinate.position == REPEAT_TRACT_END


def test_no_degraded_message_echoes_a_coordinate_or_an_allele() -> None:
    """PRIV-09: warnings reach terminals, logs, crash reports and agent context."""
    result = normalise_variants(
        [make_record(*RIGHTMOST_INSERTION), make_record(REPEAT_TRACT_START, "A", "AC")]
    )
    text = " ".join((*result.warnings, result.left_alignment.describe()))
    for forbidden in ("chr21", "26", "GRCh38:chr21", "C:CC", "0/1"):
        assert forbidden not in text, f"{forbidden!r} leaked into a surfaced message"


def test_repeat_runs_are_byte_identical(fasta_reference: FastaReference) -> None:
    """GP-30, over the whole result object including the degraded report."""
    records = [
        make_record(*RIGHTMOST_INSERTION),
        make_record(5, "G", "T"),
        make_record(REPEAT_TRACT_START, "A", "AC"),
    ]

    def digest(reference: FastaReference | None) -> str:
        result = normalise_variants(records, reference=reference)
        return stable_hash(
            {
                "variants": [record.model_dump(mode="json") for record in result.variants],
                "operations": result.operations_applied,
                "warnings": list(result.warnings),
                "left_alignment": result.left_alignment.as_dict(),
            }
        )

    assert digest(fasta_reference) == digest(fasta_reference)
    assert digest(None) == digest(None)
    assert digest(fasta_reference) != digest(None), "the two states must be distinguishable"


# ---------------------------------------------------------------------------
# 4. One rule, two callers
# ---------------------------------------------------------------------------

MINIMAL_CLINVAR_VCF = (
    "##fileformat=VCFv4.1\n"
    "##fileDate=2026-08-22\n"
    "##source=ClinVar\n"
    "##reference=GRCh38\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "21\t21\t1\tA\tAC\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel\n"
)


def _indexed_vcf(directory: Path, text: str, *, name: str = "mini.vcf") -> Path:
    plain = directory / name
    plain.write_text(text, encoding="utf-8")
    compressed = directory / f"{name}.gz"
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)
    return compressed


@pytest.fixture(scope="module")
def paired_adapters(
    tmp_path_factory: pytest.TempPathFactory,
) -> Iterator[tuple[ClinvarVcfAdapter, ClinvarVcfAdapter, InMemoryReference]]:
    """The ClinVar adapter with and without a reference, plus that reference."""
    directory = tmp_path_factory.mktemp("paired")
    path = _indexed_vcf(directory, MINIMAL_CLINVAR_VCF)
    digest = hash_file(path)
    reference = InMemoryReference("chr21", CHR21_SEQ)
    without = ClinvarVcfAdapter(path, expected_sha256=digest)
    with_reference = ClinvarVcfAdapter(path, expected_sha256=digest, reference=reference)
    yield without, with_reference, reference
    without.close()
    with_reference.close()


@given(
    ref=st.text(alphabet="ACGT", min_size=1, max_size=8),
    alt=st.text(alphabet="ACGT", min_size=1, max_size=8),
    position=st.integers(min_value=1, max_value=len(CHR21_SEQ) - 8),
)
@settings(
    max_examples=400, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
def test_adapter_and_ingestion_canonicalisation_can_never_disagree(
    paired_adapters: tuple[ClinvarVcfAdapter, ClinvarVcfAdapter, InMemoryReference],
    ref: str,
    alt: str,
    position: int,
) -> None:
    """The property the whole refactor exists to guarantee.

    Before this change the two sides were separate implementations: ingestion
    trimmed, the adapter did not, and the disagreement surfaced only as a missing
    clinical assertion. A single shared function makes agreement structural, and
    this test is what stops a future edit from quietly re-introducing a second
    copy — including in the harder direction, where a reference is available and
    both sides must also left-align identically.
    """
    assume(ref != alt)
    without_reference, with_reference, reference = paired_adapters
    record = make_record(position, ref, alt)

    for adapter, lookup in ((without_reference, None), (with_reference, reference)):
        from_adapter = adapter._canonicalise("chr21", position, ref, alt)
        from_ingestion = trim_and_left_align(record, lookup).coordinate
        assert (from_adapter.position, from_adapter.ref, from_adapter.alt) == (
            from_ingestion.position,
            from_ingestion.ref,
            from_ingestion.alt,
        )


@given(
    ref=st.text(alphabet="ACGT", min_size=1, max_size=8),
    alt=st.text(alphabet="ACGT", min_size=1, max_size=8),
    position=st.integers(min_value=1, max_value=len(CHR21_SEQ) - 8),
)
@settings(max_examples=400, deadline=None)
def test_canonicalisation_is_minimal_stable_and_never_empties_an_allele(
    ref: str, alt: str, position: int
) -> None:
    """The invariants of Tan et al. (2015), stated as properties rather than cases."""
    assume(ref != alt)
    reference = InMemoryReference("chr21", CHR21_SEQ)
    once = canonicalise_allele(
        contig="chr21", position=position, ref=ref, alt=alt, reference=reference
    )
    twice = canonicalise_allele(
        contig="chr21",
        position=once.position,
        ref=once.ref,
        alt=once.alt,
        reference=reference,
    )

    assert (once.position, once.ref, once.alt) == (twice.position, twice.ref, twice.alt)
    assert len(once.ref) >= 1
    assert len(once.alt) >= 1
    assert once.ref != once.alt
    assert once.position >= 1
    # Minimal: the canonical form shares neither a trimmable prefix nor suffix.
    assert trim_parsimoniously(once.position, once.ref, once.alt) == (
        once.position,
        once.ref,
        once.alt,
    )
    # A left shift may only move the coordinate leftwards; a trim only rightwards.
    if once.left_aligned:
        assert once.position < position or once.trimmed


def test_symbolic_and_missing_alleles_are_returned_untouched() -> None:
    """``*`` and ``.`` are legal VCF and have no defined trimming or shifting.

    Mangling them would move a coordinate on a guess. They are filtered earlier in
    ingestion; the rule still has to be safe if one reaches it.
    """
    for allele in ("*", ".", "<DEL>", ""):
        result = canonicalise_allele(contig="chr21", position=100, ref="A", alt=allele)
        assert (result.position, result.ref, result.alt) == (100, "A", allele)
        assert result.operations == ()
