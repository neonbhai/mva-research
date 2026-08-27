"""Unit tests for the MANE gene-interval adapter.

Everything here runs against ``tests/fixtures/mane/MANE.GRCh38.v1.5.slice.*`` — a
25-gene slice cut out of the genuine MANE GRCh38 v1.5 release by
``tests/fixtures/mane/make_fixture.py``, which documents every region and why it
was chosen. The gene spans, the Ensembl and HGNC identifiers, the RefSeq
accessions and the symbol drift are all the release's own, not values written to
match the parser.

The interesting failures in a gene-assignment adapter are not "wrong symbol
returned". They are:

* a variant in two overlapping genes being reported in one, which deletes the
  compound-heterozygote hypothesis for the other gene entirely;
* an intergenic variant coming back as an empty-but-present entry, or as a
  silently invented nearest gene;
* ``chr15`` being compared raw against ``NC_000015.10``, which fails by finding
  nothing and is therefore indistinguishable from "this variant is intergenic";
* a one-base error at a gene boundary, which loses variants at exactly the exon
  edges where damaging variants concentrate;
* a GRCh37 coordinate joined against a GRCh38 gene model, which returns the wrong
  gene or no gene while looking entirely successful.

Each of those has a test below, and the last section proves the point of the
whole module: ``generate_pairs``, which returns nothing at all today, returns a
BUB1B compound-heterozygous candidate once these annotations are attached.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import pytest
from pydantic import ValidationError

from mva.annotation.base import AdapterRole, AdapterSet, ConsequenceAdapter, is_synthetic
from mva.annotation.gene_intervals import (
    GENE_LOCUS_TERM,
    IMPACT_NOT_ASSESSED_REMEDIATION,
    INTERGENIC_TERM,
    MANE_ADAPTER_NAME,
    MANE_PLUS_CLINICAL_STATUS,
    MANE_SELECT_STATUS,
    REFSEQ_ACCESSION_TO_CONTIG,
    UNREPRESENTABLE_MANE_FIELDS,
    GeneBackfillConsequenceAdapter,
    GeneInterval,
    GeneRelation,
    ManeGeneAdapter,
    ManeGeneIndex,
    contig_for_refseq_accession,
    gzip_stored_filename,
    model_accepts_unassessed_impact,
)
from mva.determinism import hash_file, stable_hash
from mva.errors import AdapterUnavailableError, GenomeBuildMismatchError
from mva.models.genome import GenomeBuild, GenomicCoordinate
from mva.models.variant import (
    ConsequenceAnnotation,
    FilterStatus,
    Genotype,
    ImpactSeverity,
    VariantRecord,
    Zygosity,
)
from mva.prioritization.pairing import generate_pairs

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "mane"
FIXTURE_GTF = FIXTURE_DIR / "MANE.GRCh38.v1.5.slice.ensembl_genomic.gtf.gz"
FIXTURE_SUMMARY = FIXTURE_DIR / "MANE.GRCh38.v1.5.slice.summary.txt.gz"

#: The release the fixture was cut from, as the adapter reads it back.
EXPECTED_VERSION = "GRCh38-MANE-v1.5"

# --------------------------------------------------------------------------- anchors
#
# Real MANE v1.5 gene spans, 1-based and inclusive at both ends, exactly as the
# GTF states them. These are public reference coordinates, not patient data.


@dataclass(frozen=True)
class Locus:
    """One real gene from the fixture, named so a test reads as biology."""

    symbol: str
    contig: str
    start: int
    end: int
    strand: str


#: The MVA gene this project exists for. Note the *gene* span: the MANE Select
#: transcript is 40,161,069-40,221,123, 85 bases narrower at the 5' end.
BUB1B = Locus("BUB1B", "chr15", 40_160_984, 40_221_137, "+")
BUB1B_TRANSCRIPT_START = 40_161_069

#: Overlaps BUB1B's last 3,710 bases, on the same strand.
PAK6 = Locus("PAK6", "chr15", 40_217_428, 40_277_489, "+")

#: Ends 52,056 bases before BUB1B begins; the gap between them is intergenic.
BMF = Locus("BMF", "chr15", 40_087_889, 40_108_928, "-")

#: ANKRD63 sits entirely inside PLCB2 — a nested locus, which is what breaks a
#: sorted-by-start scan that stops at the first interval that does not overlap.
PLCB2 = Locus("PLCB2", "chr15", 40_278_176, 40_307_965, "-")
ANKRD63 = Locus("ANKRD63", "chr15", 40_278_372, 40_283_064, "-")

#: An MVA gene that overlaps a different gene on each side, both on the opposite
#: strand: FAM76B below it and MTMR2 above it.
CEP57 = Locus("CEP57", "chr11", 95_789_965, 95_837_070, "+")
FAM76B = Locus("FAM76B", "chr11", 95_768_953, 95_790_409, "-")
MTMR2 = Locus("MTMR2", "chr11", 95_821_766, 95_925_315, "-")

#: An MVA gene whose GRCh37 coordinate is around chr5:600,000 — see
#: ``GRCH37_TRIP13_POSITION``, which is 293 kb away and in no gene at all.
TRIP13 = Locus("TRIP13", "chr5", 892_849, 919_357, "+")
BRD9 = Locus("BRD9", "chr5", 850_291, 892_838, "-")
ZDHHC11 = Locus("ZDHHC11", "chr5", 795_605, 858_973, "-")
CEP72 = Locus("CEP72", "chr5", 612_305, 667_703, "+")

#: The gene HGNC renamed from CENPJ. ``CENPJ`` is absent from MANE v1.5 entirely.
CPAP = Locus("CPAP", "chr13", 24_882_276, 24_922_889, "-")

#: Two MANE transcripts for one gene: a Select and a Plus Clinical.
MECP2 = Locus("MECP2", "chrX", 154_021_573, 154_137_103, "-")

#: 51 bases, the smallest interval in the fixture, and the one gene whose two
#: files disagree about its name.
TINY = Locus("RP11-148K1.15", "chr7", 151_061_928, 151_061_978, "+")
TINY_NCBI_SYMBOL = "LOC128092247"

#: The 10-base intergenic gap between BRD9's last base and TRIP13's first.
INTERGENIC_GAP_START = BRD9.end + 1
INTERGENIC_GAP_END = TRIP13.start - 1

#: The coordinate a GRCh37 TRIP13 lookup lands on. In GRCh38 it is in no gene,
#: 12,305 bases short of CEP72 and 292,849 short of TRIP13 itself.
GRCH37_TRIP13_POSITION = 600_000

#: Inside the 52 kb BMF-to-BUB1B gap.
INTERGENIC_CHR15 = 40_130_000


# --------------------------------------------------------------------------- helpers


@pytest.fixture(scope="module")
def index() -> ManeGeneIndex:
    """The fixture release, loaded once. Construction verifies both hashes."""
    return open_index()


def open_index(
    gtf: Path = FIXTURE_GTF, summary: Path = FIXTURE_SUMMARY, **kwargs: object
) -> ManeGeneIndex:
    """Open a release, pinning it to whatever bytes it currently holds."""
    return ManeGeneIndex(
        gtf,
        summary,
        expected_gtf_sha256=hash_file(gtf),
        expected_summary_sha256=hash_file(summary),
        **kwargs,  # type: ignore[arg-type]
    )


#: Whether ``ConsequenceAnnotation`` can yet hold "impact not assessed".
#:
#: ``ManeGeneAdapter`` refuses to construct while this is False, and offers no way
#: to override that refusal, so every test of the Protocol shim below is gated on
#: it. They are skips rather than deletions on purpose: the day the remediation in
#: :data:`IMPACT_NOT_ASSESSED_REMEDIATION` lands, this flips and the whole shim is
#: covered with no edit here. The interval join itself — where all the logic lives —
#: is tested unconditionally through ``ManeGeneIndex``.
MODEL_CAN_REPRESENT_UNASSESSED_IMPACT = model_accepts_unassessed_impact()

requires_nullable_impact = pytest.mark.skipif(
    not MODEL_CAN_REPRESENT_UNASSESSED_IMPACT,
    reason=(
        "ConsequenceAnnotation.impact is still a required ImpactSeverity, so "
        "ManeGeneAdapter refuses to construct rather than fabricate one. Apply the "
        "remediation in IMPACT_NOT_ASSESSED_REMEDIATION to enable these."
    ),
)


@pytest.fixture(scope="module")
def adapter(index: ManeGeneIndex) -> ManeGeneAdapter:
    """The Protocol shim, once the model can represent an unassessed impact."""
    if not MODEL_CAN_REPRESENT_UNASSESSED_IMPACT:
        pytest.skip("ConsequenceAnnotation.impact cannot yet represent 'not assessed'.")
    return ManeGeneAdapter(index)


def symbols_at(index: ManeGeneIndex, contig: str, start: int, end: int | None = None) -> list[str]:
    return [gene.gene_symbol for gene in index.genes_at(contig, start, end)]


def gene(index: ManeGeneIndex, symbol: str) -> GeneInterval:
    found = index.genes_for_symbol(symbol)
    assert len(found) == 1, f"{symbol} should resolve to exactly one gene, got {len(found)}"
    return found[0]


def variant_id(contig: str, position: int, ref: str = "C", alt: str = "T") -> str:
    return f"GRCh38:{contig}:{position}:{ref}:{alt}"


def make_record(
    contig: str,
    position: int,
    consequences: tuple[ConsequenceAnnotation, ...] = (),
    *,
    ref: str = "C",
    alt: str = "T",
    zygosity: Zygosity = Zygosity.HET,
) -> VariantRecord:
    return VariantRecord(
        coordinate=GenomicCoordinate(
            build=GenomeBuild.GRCH38, contig=contig, position=position, ref=ref, alt=alt
        ),
        genotype=Genotype(zygosity=zygosity, genotype_string="0/1"),
        filter_status=FilterStatus.PASS,
        source_artifact="test",
        consequences=consequences,
    )


GTF_TEMPLATE = (
    "{contig}\tensembl_havana\tgene\t{start}\t{end}\t.\t{strand}\t.\t"
    'gene_id "{gene_id}"; gene_type "protein_coding"; gene_name "{symbol}";\n'
)

SUMMARY_HEADER = (
    "#NCBI_GeneID\tEnsembl_Gene\tHGNC_ID\tsymbol\tname\tRefSeq_nuc\tRefSeq_prot\t"
    "Ensembl_nuc\tEnsembl_prot\tMANE_status\tGRCh38_chr\tchr_start\tchr_end\tchr_strand\n"
)

SUMMARY_TEMPLATE = (
    "GeneID:1\t{gene_id}\tHGNC:1\t{symbol}\ta gene\tNM_1.1\tNP_1.1\tENST0000000{n}.1\t"
    "ENSP0000000{n}.1\t{status}\t{accession}\t{start}\t{end}\t{strand}\n"
)


def write_release(
    directory: Path,
    *,
    gtf_rows: Iterable[str],
    summary_rows: Iterable[str],
    prefix: str = "MANE.GRCh38.v1.5.",
    gtf_prefix: str | None = None,
    stored_prefix: str | None = None,
) -> tuple[Path, Path]:
    """Write a small, hand-built MANE release pair for a failure-path test.

    ``stored_prefix`` overrides the release name gzip records *inside* the
    compressed bytes, which is how the "file was renamed after distribution" case
    is reproduced.
    """
    gtf_path = directory / f"{gtf_prefix or prefix}ensembl_genomic.gtf.gz"
    summary_path = directory / f"{prefix}summary.txt.gz"
    _write_gz(gtf_path, "".join(gtf_rows), stored_prefix, "ensembl_genomic.gtf")
    _write_gz(
        summary_path,
        SUMMARY_HEADER + "".join(summary_rows),
        stored_prefix,
        "summary.txt",
    )
    return gtf_path, summary_path


def _write_gz(path: Path, text: str, stored_prefix: str | None, suffix: str) -> None:
    stored_name = None if stored_prefix is None else f"{stored_prefix}{suffix}"
    with (
        path.open("wb") as raw,
        gzip.GzipFile(
            filename=stored_name or path.name.removesuffix(".gz"), mode="wb", fileobj=raw, mtime=0
        ) as handle,
    ):
        handle.write(text.encode("utf-8"))


def minimal_release(directory: Path, **kwargs: object) -> tuple[Path, Path]:
    """One gene, valid in both files. The baseline the failure paths deviate from."""
    return write_release(
        directory,
        gtf_rows=[
            GTF_TEMPLATE.format(
                contig="chr15", start=100, end=200, strand="+", gene_id="ENSG1.1", symbol="AAA"
            )
        ],
        summary_rows=[
            SUMMARY_TEMPLATE.format(
                gene_id="ENSG1.1",
                symbol="AAA",
                n=1,
                status=MANE_SELECT_STATUS,
                accession="NC_000015.10",
                start=100,
                end=200,
                strand="+",
            )
        ],
        **kwargs,  # type: ignore[arg-type]
    )


# --------------------------------------------------------------------------- fixture identity


def test_fixture_is_a_real_mane_release_slice(index: ManeGeneIndex) -> None:
    """The committed slice really is MANE v1.5, not rows written to fit the parser."""
    assert index.gene_count == 25
    bub1b = gene(index, "BUB1B")
    assert (bub1b.contig, bub1b.start, bub1b.end, bub1b.strand) == (
        BUB1B.contig,
        BUB1B.start,
        BUB1B.end,
        BUB1B.strand,
    )
    assert bub1b.gene_id == "ENSG00000156970.15"
    assert bub1b.hgnc_id == "HGNC:1149"
    assert bub1b.ncbi_gene_id == "GeneID:701"
    assert bub1b.mane_select.transcript_id == "ENST00000287598.11"
    assert bub1b.mane_select.refseq_nuc.startswith("NM_")


def test_every_mva_gene_in_the_fixture_is_at_its_grch38_span(index: ManeGeneIndex) -> None:
    for locus in (BUB1B, CEP57, TRIP13, CPAP):
        found = gene(index, locus.symbol)
        assert (found.contig, found.start, found.end) == (locus.contig, locus.start, locus.end)


def test_version_is_read_from_the_release_not_invented(index: ManeGeneIndex) -> None:
    assert index.version == EXPECTED_VERSION


def test_the_gtf_carries_no_in_band_version_which_is_why_the_name_is_used() -> None:
    """The premise of :func:`_release_version`, asserted rather than assumed.

    If a future MANE GTF gains a ``#!`` header line, this fails and the version
    should start being read from it instead of from the filename.
    """
    with gzip.open(FIXTURE_GTF, "rt", encoding="utf-8") as handle:
        assert not any(line.startswith("#") for line in handle)


def test_the_summary_carries_its_release_name_inside_the_compressed_bytes() -> None:
    """One copy of the release string lives under the sha256 pin, not beside it."""
    stored = gzip_stored_filename(FIXTURE_SUMMARY)
    assert stored is not None
    assert stored.startswith("MANE.GRCh38.v1.5.")


@requires_nullable_impact
def test_adapter_declares_itself_real(adapter: ManeGeneAdapter) -> None:
    assert adapter.synthetic is False
    assert is_synthetic(adapter) is False
    assert adapter.name == MANE_ADAPTER_NAME


@requires_nullable_impact
def test_adapter_name_cannot_be_mistaken_for_an_effect_predictor(
    adapter: ManeGeneAdapter,
) -> None:
    lowered = adapter.name.lower()
    assert "join" in lowered
    for predictor in ("vep", "snpeff", "nirvana", "annovar"):
        assert predictor not in lowered


@requires_nullable_impact
def test_satisfies_the_consequence_adapter_protocol(adapter: ManeGeneAdapter) -> None:
    assert isinstance(adapter, ConsequenceAdapter)


@requires_nullable_impact
def test_adapter_set_reports_the_consequence_slot_as_real(adapter: ManeGeneAdapter) -> None:
    class _Frequency:
        name = "stub-frequency"
        version = "v0"

        def frequencies(self, variant_ids: Sequence[str]) -> Mapping[str, tuple[()]]:
            return {}

    descriptors = AdapterSet(consequence=adapter, frequency=_Frequency()).descriptors()
    consequence = next(d for d in descriptors if d.role is AdapterRole.CONSEQUENCE)
    assert consequence.synthetic is False
    assert consequence.label.startswith(f"{MANE_ADAPTER_NAME}@{EXPECTED_VERSION}")


def test_unrepresentable_fields_are_named_not_silently_dropped() -> None:
    assert "HGNC_ID" in UNREPRESENTABLE_MANE_FIELDS
    assert "strand" in UNREPRESENTABLE_MANE_FIELDS
    assert "distance_bp" in UNREPRESENTABLE_MANE_FIELDS


# --------------------------------------------------------------------------- contig naming


def test_refseq_map_covers_every_primary_grch38_sequence() -> None:
    assert contig_for_refseq_accession("NC_000001.11") == "chr1"
    assert contig_for_refseq_accession("NC_000022.11") == "chr22"
    assert contig_for_refseq_accession("NC_000023.11") == "chrX"
    assert contig_for_refseq_accession("NC_000024.10") == "chrY"
    assert contig_for_refseq_accession("NC_012920.1") == "chrM"
    assert len(REFSEQ_ACCESSION_TO_CONTIG) == 25


def test_the_two_spellings_of_a_chromosome_share_no_substring() -> None:
    """Why the mapping must be explicit: neither name can be derived from the other."""
    accession, contig = "NC_000023.11", "chrX"
    assert "X" not in accession
    assert "23" not in contig
    assert contig_for_refseq_accession(accession) == contig


def test_accession_version_is_ignored_so_a_patch_bump_does_not_lose_a_chromosome() -> None:
    """A version-pinned map would fail by finding no genes, which reads as intergenic."""
    assert contig_for_refseq_accession("NC_000015.10") == "chr15"
    assert contig_for_refseq_accession("NC_000015.11") == "chr15"
    assert contig_for_refseq_accession("NC_000015") == "chr15"


def test_scaffold_accessions_resolve_to_no_contig() -> None:
    for accession in ("NT_187633.1", "NW_015148966.2", "NC_000025.1", "banana"):
        assert contig_for_refseq_accession(accession) is None


def test_chrx_genes_join_across_both_spellings(index: ManeGeneIndex) -> None:
    """MECP2 says ``chrX`` in the GTF and ``NC_000023.11`` in the summary."""
    mecp2 = gene(index, "MECP2")
    assert mecp2.contig == "chrX"
    assert mecp2.hgnc_id == "HGNC:6990"
    assert symbols_at(index, "chrX", MECP2.start) == ["MECP2"]


def test_both_files_skip_exactly_the_same_off_assembly_genes(index: ManeGeneIndex) -> None:
    """MUC2 on a fix patch and GSTT1 on an alt haplotype, dropped by both readers."""
    counts = index.skipped_gene_counts
    assert counts["gtf"] == counts["summary"] == 2
    assert index.genes_for_symbol("MUC2") == ()
    assert index.genes_for_symbol("GSTT1") == ()


def test_only_canonical_contigs_are_indexed(index: ManeGeneIndex) -> None:
    assert set(index.contigs) == {"chr5", "chr7", "chr11", "chr13", "chr15", "chrX"}
    assert all("_" not in contig for contig in index.contigs)


def test_a_contig_disagreement_between_the_files_is_refused(tmp_path: Path) -> None:
    """The single most dangerous detail: the two spellings must describe one genome."""
    gtf_path, summary_path = write_release(
        tmp_path,
        gtf_rows=[
            GTF_TEMPLATE.format(
                contig="chr15", start=100, end=200, strand="+", gene_id="ENSG1.1", symbol="AAA"
            )
        ],
        summary_rows=[
            SUMMARY_TEMPLATE.format(
                gene_id="ENSG1.1",
                symbol="AAA",
                n=1,
                status=MANE_SELECT_STATUS,
                accession="NC_000011.10",  # chr11, not chr15
                start=100,
                end=200,
                strand="+",
            )
        ],
    )
    with pytest.raises(AdapterUnavailableError, match="different contigs"):
        open_index(gtf_path, summary_path)


def test_a_grch37_variant_id_is_refused_not_silently_missed(index: ManeGeneIndex) -> None:
    """A silent miss would read as 'every variant in this batch is intergenic'."""
    with pytest.raises(GenomeBuildMismatchError, match="GRCh37"):
        index.assign_variants([f"GRCh37:chr15:{BUB1B.start}:C:T"])


def test_the_grch37_trip13_coordinate_finds_nothing_in_grch38(index: ManeGeneIndex) -> None:
    """The live hazard, reproduced. It fails by finding nothing, not by raising."""
    assert index.genes_at("chr5", GRCH37_TRIP13_POSITION) == ()
    assert index.assign_variants([variant_id("chr5", GRCH37_TRIP13_POSITION)]) == {}
    # ... while the same lookup at the true GRCh38 span does find the gene.
    assert symbols_at(index, "chr5", TRIP13.start) == ["TRIP13"]


def test_ensembl_style_contigs_are_normalised_on_the_way_in(index: ManeGeneIndex) -> None:
    """A caller passing ``15`` gets the same answer as one passing ``chr15``."""
    assert symbols_at(index, "15", BUB1B.start) == symbols_at(index, "chr15", BUB1B.start)


def test_a_non_canonical_contig_is_refused_rather_than_missed(index: ManeGeneIndex) -> None:
    with pytest.raises(ValueError, match="Non-canonical contig"):
        index.genes_at("chr15_KI270905v1_alt", 1000)


def test_a_contig_with_no_indexed_genes_returns_nothing(index: ManeGeneIndex) -> None:
    assert index.genes_at("chr1", 1_000_000) == ()


# --------------------------------------------------------------------------- boundaries
#
# GTF start/end are 1-based and inclusive at both ends; VCF POS is 1-based. The
# module does no half-open conversion anywhere, so these four assertions per gene
# are the whole of the boundary contract.


@pytest.mark.parametrize(
    "locus", [BUB1B, CPAP, TRIP13, CEP72, TINY], ids=lambda locus: locus.symbol
)
def test_first_and_last_base_of_a_gene_are_inside_it(index: ManeGeneIndex, locus: Locus) -> None:
    assert locus.symbol in symbols_at(index, locus.contig, locus.start)
    assert locus.symbol in symbols_at(index, locus.contig, locus.end)


@pytest.mark.parametrize(
    "locus", [BUB1B, CPAP, TRIP13, CEP72, TINY], ids=lambda locus: locus.symbol
)
def test_one_base_either_side_of_a_gene_is_outside_it(index: ManeGeneIndex, locus: Locus) -> None:
    assert locus.symbol not in symbols_at(index, locus.contig, locus.start - 1)
    assert locus.symbol not in symbols_at(index, locus.contig, locus.end + 1)


def test_the_boundary_holds_for_a_minus_strand_gene(index: ManeGeneIndex) -> None:
    """Strand does not flip the interval: MANE states start < end regardless."""
    cpap = gene(index, "CPAP")
    assert cpap.strand == "-"
    assert cpap.start < cpap.end
    assert cpap.contains(CPAP.start) and cpap.contains(CPAP.end)
    assert not cpap.contains(CPAP.start - 1)
    assert not cpap.contains(CPAP.end + 1)


def test_a_single_base_gene_lookup_defaults_to_a_one_base_span(index: ManeGeneIndex) -> None:
    assert index.genes_at("chr15", BUB1B.start) == index.genes_at("chr15", BUB1B.start, BUB1B.start)


def test_the_gene_span_is_wider_than_the_mane_transcript_span(index: ManeGeneIndex) -> None:
    """Using the summary's transcript span would drop 85 bases of BUB1B's 5' end.

    That window is small, and it is exactly at the UTR boundary, which is why the
    index is built from the GTF ``gene`` rows rather than the summary's
    ``chr_start``/``chr_end``.
    """
    assert BUB1B.start < BUB1B_TRANSCRIPT_START
    for position in (BUB1B.start, BUB1B_TRANSCRIPT_START - 1):
        assert "BUB1B" in symbols_at(index, "chr15", position)


def test_a_deletion_reaching_into_a_gene_from_outside_overlaps_it(index: ManeGeneIndex) -> None:
    """A REF span, not a point: ``[POS, POS + len(REF) - 1]``."""
    outside = BUB1B.start - 3
    assert index.genes_at("chr15", outside) == ()
    assigned = index.assign_variants([variant_id("chr15", outside, ref="CGTA", alt="C")])
    assert [a.gene.gene_symbol for a in assigned[variant_id("chr15", outside, "CGTA", "C")]] == [
        "BUB1B"
    ]


def test_a_deletion_ending_one_base_short_of_a_gene_does_not_overlap_it(
    index: ManeGeneIndex,
) -> None:
    outside = BUB1B.start - 4
    assert index.assign_variants([variant_id("chr15", outside, ref="CGT", alt="C")]) == {}


def test_a_span_end_before_its_start_is_refused(index: ManeGeneIndex) -> None:
    with pytest.raises(ValueError, match="precedes span start"):
        index.genes_at("chr15", 500, 400)


def test_gap_to_is_zero_when_the_span_touches_the_gene(index: ManeGeneIndex) -> None:
    bub1b = gene(index, "BUB1B")
    assert bub1b.gap_to(BUB1B.start, BUB1B.start) == 0
    assert bub1b.gap_to(BUB1B.end, BUB1B.end) == 0
    assert bub1b.gap_to(BUB1B.start - 1, BUB1B.start - 1) == 1
    assert bub1b.gap_to(BUB1B.end + 1, BUB1B.end + 1) == 1
    assert bub1b.gap_to(BUB1B.end + 10, BUB1B.end + 20) == 10


def test_gene_length_counts_both_endpoints(index: ManeGeneIndex) -> None:
    assert gene(index, TINY.symbol).length == TINY.end - TINY.start + 1 == 51


# --------------------------------------------------------------------------- multi-gene


def test_a_variant_in_two_overlapping_genes_reports_both(index: ManeGeneIndex) -> None:
    """CEP57 and FAM76B overlap, on opposite strands. Both are real answers."""
    assert symbols_at(index, "chr11", CEP57.start) == ["FAM76B", "CEP57"]


def test_the_same_gene_overlaps_a_different_neighbour_on_its_other_side(
    index: ManeGeneIndex,
) -> None:
    assert symbols_at(index, "chr11", MTMR2.start) == ["CEP57", "MTMR2"]


def test_an_mva_gene_shares_its_tail_with_the_next_gene(index: ManeGeneIndex) -> None:
    """A BUB1B variant in the shared window is also a PAK6 variant, and stays one."""
    assert symbols_at(index, "chr15", PAK6.start) == ["BUB1B", "PAK6"]
    assert symbols_at(index, "chr15", BUB1B.end) == ["BUB1B", "PAK6"]


def test_a_nested_gene_does_not_hide_its_container(index: ManeGeneIndex) -> None:
    """ANKRD63 is entirely inside PLCB2; a naive sorted scan returns only one."""
    assert symbols_at(index, "chr15", ANKRD63.start) == ["PLCB2", "ANKRD63"]
    assert symbols_at(index, "chr15", ANKRD63.end) == ["PLCB2", "ANKRD63"]
    # ... and the container alone, just outside the nested gene.
    assert symbols_at(index, "chr15", ANKRD63.end + 1) == ["PLCB2"]


def test_two_minus_strand_genes_that_overlap_are_both_reported(index: ManeGeneIndex) -> None:
    """ZDHHC11 and BRD9 share 8,683 bases, both on the minus strand."""
    assert symbols_at(index, "chr5", BRD9.start) == ["ZDHHC11", "BRD9"]
    assert symbols_at(index, "chr5", ZDHHC11.end) == ["ZDHHC11", "BRD9"]
    assert symbols_at(index, "chr5", ZDHHC11.end + 1) == ["BRD9"]


def test_a_span_covering_several_genes_returns_every_one(index: ManeGeneIndex) -> None:
    assert symbols_at(index, "chr5", ZDHHC11.start, TRIP13.end) == [
        "ZDHHC11",
        "BRD9",
        "TRIP13",
    ]


def test_overlapping_genes_come_back_in_the_documented_total_order(
    index: ManeGeneIndex,
) -> None:
    """``(distance, contig, start, end, gene_id)`` — never set or dict order."""
    assignments = index.assign("chr11", CEP57.start)
    assert [a.gene.gene_symbol for a in assignments] == ["FAM76B", "CEP57"]
    keys = [a.sort_key() for a in assignments]
    assert keys == sorted(keys)
    assert all(key[0] == 0 for key in keys)  # every overlap has distance 0


@requires_nullable_impact
def test_no_overlapping_gene_is_dropped_by_the_adapter(adapter: ManeGeneAdapter) -> None:
    identifier = variant_id("chr11", CEP57.start)
    annotations = adapter.annotate([identifier])[identifier]
    assert [a.gene_symbol for a in annotations] == ["FAM76B", "CEP57"]


@requires_nullable_impact
def test_a_multi_gene_variant_reaches_gene_symbols_as_two_genes(
    adapter: ManeGeneAdapter,
) -> None:
    identifier = variant_id("chr11", CEP57.start)
    record = make_record("chr11", CEP57.start, adapter.annotate([identifier])[identifier])
    assert record.gene_symbols == ("FAM76B", "CEP57")


# --------------------------------------------------------------------------- absence


def test_an_intergenic_position_is_in_no_gene(index: ManeGeneIndex) -> None:
    """The 10-base gap between BRD9 and TRIP13 belongs to neither."""
    assert INTERGENIC_GAP_END - INTERGENIC_GAP_START + 1 == 10
    for position in (INTERGENIC_GAP_START, INTERGENIC_GAP_END):
        assert index.genes_at("chr5", position) == ()
    # One base either side of the gap is inside a gene, so the gap is exact.
    assert symbols_at(index, "chr5", INTERGENIC_GAP_START - 1) == ["BRD9"]
    assert symbols_at(index, "chr5", INTERGENIC_GAP_END + 1) == ["TRIP13"]


def test_an_intergenic_variant_is_omitted_from_the_mapping_not_empty_in_it(
    index: ManeGeneIndex,
) -> None:
    """GP-14: 'no gene here' and 'we did not look' must stay different facts."""
    identifier = variant_id("chr15", INTERGENIC_CHR15)
    assigned = index.assign_variants([identifier, variant_id("chr15", BUB1B.start)])
    assert identifier not in assigned
    assert assigned.get(identifier) is None
    assert variant_id("chr15", BUB1B.start) in assigned


def test_an_entirely_intergenic_batch_returns_an_empty_mapping(index: ManeGeneIndex) -> None:
    assigned = index.assign_variants(
        [variant_id("chr5", INTERGENIC_GAP_START), variant_id("chr15", INTERGENIC_CHR15)]
    )
    assert assigned == {}


@requires_nullable_impact
def test_the_adapter_omits_intergenic_variants_too(adapter: ManeGeneAdapter) -> None:
    identifier = variant_id("chr15", INTERGENIC_CHR15)
    assert adapter.annotate([identifier]) == {}


def test_nearest_gene_is_off_by_default(index: ManeGeneIndex) -> None:
    assert index.assign("chr15", INTERGENIC_CHR15) == ()
    assert index.assign_variants([variant_id("chr15", INTERGENIC_CHR15)]) == {}


def test_nearest_gene_when_asked_for_is_labelled_as_an_inference(
    index: ManeGeneIndex,
) -> None:
    assignments = index.assign("chr5", INTERGENIC_GAP_START, nearest_within_bp=100)
    assert [a.gene.gene_symbol for a in assignments] == ["BRD9"]
    assert assignments[0].relation is GeneRelation.NEAREST
    assert assignments[0].is_overlap is False
    assert assignments[0].distance_bp == 1


def test_nearest_gene_never_displaces_an_overlap(index: ManeGeneIndex) -> None:
    """Inside a gene, asking for nearest changes nothing at all."""
    with_nearest = index.assign("chr15", BUB1B.start, nearest_within_bp=1_000_000)
    without = index.assign("chr15", BUB1B.start)
    assert with_nearest == without
    assert all(a.relation is GeneRelation.OVERLAP for a in with_nearest)


def test_nearest_gene_at_the_grch37_coordinate_names_the_wrong_gene(
    index: ManeGeneIndex,
) -> None:
    """Exactly why the label matters: the closest gene here is not TRIP13."""
    assignments = index.assign("chr5", GRCH37_TRIP13_POSITION, nearest_within_bp=50_000)
    assert [a.gene.gene_symbol for a in assignments] == ["CEP72"]
    assert assignments[0].relation is GeneRelation.NEAREST
    assert assignments[0].distance_bp == CEP72.start - GRCH37_TRIP13_POSITION
    assert assignments[0].gene.gene_symbol != "TRIP13"


def test_nearest_gene_respects_its_distance_limit(index: ManeGeneIndex) -> None:
    just_short = CEP72.start - GRCH37_TRIP13_POSITION - 1
    assert index.assign("chr5", GRCH37_TRIP13_POSITION, nearest_within_bp=just_short) == ()
    assert index.assign("chr5", GRCH37_TRIP13_POSITION, nearest_within_bp=just_short + 1) != ()


def test_every_gene_tied_at_the_minimum_distance_is_returned(index: ManeGeneIndex) -> None:
    """IRAK1 ends 922 bases before MECP2 begins; a two-base span sits centred.

    The gap is even, so no single base is equidistant — a two-base REF span (a
    one-base deletion, in VCF terms) is, at 461 bases from each. Both genes are
    returned rather than whichever the scan reached first.
    """
    irak1_end, mecp2_start = 154_020_650, MECP2.start
    span_start = (irak1_end + mecp2_start - 1) // 2
    assignments = index.assign("chrX", span_start, span_start + 1, nearest_within_bp=1000)
    assert [a.gene.gene_symbol for a in assignments] == ["IRAK1", "MECP2"]
    assert {a.distance_bp for a in assignments} == {461}
    assert all(a.relation is GeneRelation.NEAREST for a in assignments)


@requires_nullable_impact
def test_a_nearest_gene_annotation_carries_the_intergenic_term(
    index: ManeGeneIndex,
) -> None:
    """In the lossier annotation shape, the SO term is what marks the inference."""
    near = ManeGeneAdapter(index, nearest_within_bp=50_000)
    identifier = variant_id("chr5", GRCH37_TRIP13_POSITION)
    annotations = near.annotate([identifier])[identifier]
    assert [a.consequence_terms for a in annotations] == [(INTERGENIC_TERM,)]
    assert annotations[0].gene_symbol == "CEP72"


@requires_nullable_impact
def test_an_overlap_annotation_carries_the_gene_locus_term(adapter: ManeGeneAdapter) -> None:
    identifier = variant_id("chr15", BUB1B.start)
    annotations = adapter.annotate([identifier])[identifier]
    assert annotations[0].consequence_terms == (GENE_LOCUS_TERM,)


def test_a_negative_nearest_distance_is_refused(index: ManeGeneIndex) -> None:
    with pytest.raises(ValueError, match="must not be negative"):
        index.assign("chr15", BUB1B.start, nearest_within_bp=-1)


# --------------------------------------------------------------------------- transcripts


def test_both_mane_transcripts_of_a_gene_are_kept(index: ManeGeneIndex) -> None:
    """Collapsing Plus Clinical into Select is the loss rule 3 of the adapter README forbids."""
    mecp2 = gene(index, "MECP2")
    assert len(mecp2.transcripts) == 2
    assert [t.mane_status for t in mecp2.transcripts] == [
        MANE_SELECT_STATUS,
        MANE_PLUS_CLINICAL_STATUS,
    ]
    assert mecp2.mane_select.transcript_id == "ENST00000453960.7"
    assert mecp2.transcripts[1].transcript_id == "ENST00000303391.11"


def test_the_select_transcript_always_sorts_first(index: ManeGeneIndex) -> None:
    for interval in index.genes:
        assert interval.transcripts[0].is_select
        assert interval.mane_select is interval.transcripts[0]


@requires_nullable_impact
def test_the_adapter_emits_one_annotation_per_transcript(adapter: ManeGeneAdapter) -> None:
    identifier = variant_id("chrX", MECP2.start)
    annotations = adapter.annotate([identifier])[identifier]
    assert len(annotations) == 2
    assert [a.transcript_id for a in annotations] == [
        "ENST00000453960.7",
        "ENST00000303391.11",
    ]
    assert [a.is_mane_select for a in annotations] == [True, False]
    assert {a.gene_symbol for a in annotations} == {"MECP2"}


@requires_nullable_impact
def test_canonicality_is_not_inferred_from_mane_select(adapter: ManeGeneAdapter) -> None:
    """MANE asserts a matched transcript pair, not an Ensembl canonical flag."""
    identifier = variant_id("chr15", BUB1B.start)
    annotation = adapter.annotate([identifier])[identifier][0]
    assert annotation.is_mane_select is True
    assert annotation.is_canonical is False
    assert "is_canonical" in UNREPRESENTABLE_MANE_FIELDS


# --------------------------------------------------------------------------- identifiers


def test_cenpj_is_absent_and_cpap_is_the_same_gene(index: ManeGeneIndex) -> None:
    """The confirmed symbol-drift trap, on real MANE rows."""
    assert index.genes_for_symbol("CENPJ") == ()
    cpap = gene(index, "CPAP")
    assert cpap.gene_id == "ENSG00000151849.18"
    assert cpap.hgnc_id == "HGNC:17272"
    assert cpap.ncbi_gene_id == "GeneID:55835"


def test_a_panel_gene_that_cannot_be_located_is_reported_not_dropped(
    index: ManeGeneIndex,
) -> None:
    resolution = index.resolve_panel(["BUB1B", "CENPJ", "CEP57", "TRIP13", "CASC5"])
    assert resolution.missing == ("CENPJ", "CASC5")
    assert set(resolution.found) == {"BUB1B", "CEP57", "TRIP13"}
    assert resolution.is_complete is False
    assert "CENPJ" in resolution.describe_missing()
    assert "CPAP" in resolution.describe_missing()


def test_a_fully_resolvable_panel_reports_nothing_missing(index: ManeGeneIndex) -> None:
    resolution = index.resolve_panel(["BUB1B", "CEP57", "TRIP13", "CPAP"])
    assert resolution.missing == ()
    assert resolution.is_complete is True
    assert resolution.describe_missing() == ""


def test_symbol_lookup_is_exact_and_invents_no_aliases(index: ManeGeneIndex) -> None:
    """Guessing that CENPJ means CPAP is a curation claim with no source here."""
    assert index.genes_for_symbol("bub1b") == ()
    assert index.genes_for_symbol("BUB1B ") == ()
    assert index.genes_for_symbol("BUB1") == ()


def test_both_spellings_a_release_gives_a_gene_resolve(index: ManeGeneIndex) -> None:
    """47 genes in v1.5 are named differently by the GTF and the summary."""
    by_ensembl = index.genes_for_symbol(TINY.symbol)
    by_ncbi = index.genes_for_symbol(TINY_NCBI_SYMBOL)
    assert by_ensembl == by_ncbi
    assert by_ensembl[0].gene_id == "ENSG00000288608.1"
    assert by_ensembl[0].symbols == (TINY.symbol, TINY_NCBI_SYMBOL)


def test_an_absent_hgnc_id_is_none_not_an_empty_string(index: ManeGeneIndex) -> None:
    assert gene(index, TINY.symbol).hgnc_id is None


def test_a_gene_resolves_by_its_stable_ensembl_id(index: ManeGeneIndex) -> None:
    found = index.gene_for_id("ENSG00000156970.15")
    assert found is not None
    assert found.gene_symbol == "BUB1B"
    assert index.gene_for_id("ENSG00000000000.1") is None


@requires_nullable_impact
def test_the_annotation_carries_the_ensembl_gene_id_not_only_the_symbol(
    adapter: ManeGeneAdapter,
) -> None:
    identifier = variant_id("chr15", BUB1B.start)
    annotation = adapter.annotate([identifier])[identifier][0]
    assert annotation.gene_id == "ENSG00000156970.15"
    assert annotation.gene_symbol == "BUB1B"
    assert annotation.source_tool == MANE_ADAPTER_NAME
    assert annotation.source_tool_version.startswith(EXPECTED_VERSION)


# --------------------------------------------------------------------------- determinism


def test_repeat_lookups_are_byte_identical(index: ManeGeneIndex) -> None:
    ids = [
        variant_id("chr15", BUB1B.start),
        variant_id("chr11", CEP57.start),
        variant_id("chr5", TRIP13.end),
    ]
    first = stable_hash(_dump(index.assign_variants(ids)))
    for _ in range(3):
        assert stable_hash(_dump(index.assign_variants(ids))) == first


def test_two_independently_built_indexes_agree_byte_for_byte() -> None:
    """No module-level cache: two indexes are two objects and must still agree."""
    left, right = open_index(), open_index()
    assert left is not right
    ids = [variant_id("chr11", CEP57.start), variant_id("chr15", PAK6.start)]
    assert stable_hash(_dump(left.assign_variants(ids))) == stable_hash(
        _dump(right.assign_variants(ids))
    )
    assert [g.gene_id for g in left.genes] == [g.gene_id for g in right.genes]


def test_input_order_does_not_change_per_variant_results(index: ManeGeneIndex) -> None:
    a, b = variant_id("chr15", BUB1B.start), variant_id("chr11", CEP57.start)
    forward = index.assign_variants([a, b])
    reverse = index.assign_variants([b, a])
    assert forward[a] == reverse[a]
    assert forward[b] == reverse[b]


def test_result_keys_follow_caller_order(index: ManeGeneIndex) -> None:
    a, b = variant_id("chr11", CEP57.start), variant_id("chr15", BUB1B.start)
    assert list(index.assign_variants([a, b])) == [a, b]
    assert list(index.assign_variants([b, a])) == [b, a]


def test_duplicate_ids_are_collapsed(index: ManeGeneIndex) -> None:
    identifier = variant_id("chr15", BUB1B.start)
    assert list(index.assign_variants([identifier, identifier, identifier])) == [identifier]


def test_the_gene_list_is_in_a_stable_total_order(index: ManeGeneIndex) -> None:
    keys = [g.sort_key() for g in index.genes]
    assert keys == sorted(keys)
    assert len(set(keys)) == len(keys)


@requires_nullable_impact
def test_annotations_for_one_variant_are_stably_ordered(adapter: ManeGeneAdapter) -> None:
    identifier = variant_id("chr11", CEP57.start)
    first = adapter.annotate([identifier])[identifier]
    assert first == adapter.annotate([identifier])[identifier]
    assert [a.gene_symbol for a in first] == ["FAM76B", "CEP57"]


def _dump(assigned: Mapping[str, object]) -> list[list[object]]:
    """A canonical-JSON-safe projection of an assignment mapping."""
    rendered: list[list[object]] = []
    for key in assigned:
        assignments = assigned[key]
        assert isinstance(assignments, tuple)
        rendered.append(
            [
                key,
                [
                    [a.gene.gene_id, a.gene.gene_symbol, a.relation.value, a.distance_bp]
                    for a in assignments
                ],
            ]
        )
    return rendered


# --------------------------------------------------------------------------- integrity


def test_construction_without_an_integrity_pin_fails_closed() -> None:
    with pytest.raises(AdapterUnavailableError, match="without an integrity pin"):
        ManeGeneIndex(FIXTURE_GTF, FIXTURE_SUMMARY)


def test_each_file_needs_its_own_pin() -> None:
    with pytest.raises(AdapterUnavailableError, match="summary"):
        ManeGeneIndex(FIXTURE_GTF, FIXTURE_SUMMARY, expected_gtf_sha256=hash_file(FIXTURE_GTF))


def test_a_sha256_mismatch_refuses_to_open_the_release() -> None:
    with pytest.raises(AdapterUnavailableError, match="failed its sha256"):
        ManeGeneIndex(
            FIXTURE_GTF,
            FIXTURE_SUMMARY,
            expected_gtf_sha256="0" * 64,
            expected_summary_sha256=hash_file(FIXTURE_SUMMARY),
        )


def test_an_integrity_failure_never_echoes_the_file_contents() -> None:
    with pytest.raises(AdapterUnavailableError) as excinfo:
        ManeGeneIndex(
            FIXTURE_GTF,
            FIXTURE_SUMMARY,
            expected_gtf_sha256="0" * 64,
            expected_summary_sha256=hash_file(FIXTURE_SUMMARY),
        )
    message = str(excinfo.value)
    assert "BUB1B" not in message
    assert "ENSG" not in message
    assert FIXTURE_GTF.name in message


def test_a_missing_file_is_a_distinct_failure(tmp_path: Path) -> None:
    """Absent is reported before unpinned: the file is named, its bytes are not."""
    absent_gtf = tmp_path / "MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz"
    absent_summary = tmp_path / "MANE.GRCh38.v1.5.summary.txt.gz"
    with pytest.raises(AdapterUnavailableError, match="not found"):
        ManeGeneIndex(
            absent_gtf,
            FIXTURE_SUMMARY,
            expected_gtf_sha256="0" * 64,
            expected_summary_sha256="0" * 64,
        )
    with pytest.raises(AdapterUnavailableError, match="not found"):
        ManeGeneIndex(
            FIXTURE_GTF,
            absent_summary,
            expected_gtf_sha256="0" * 64,
            expected_summary_sha256="0" * 64,
        )


def test_a_release_without_a_mane_prefix_has_no_identifiable_version(tmp_path: Path) -> None:
    gtf_path, summary_path = minimal_release(tmp_path, prefix="mane-latest.")
    with pytest.raises(AdapterUnavailableError, match="release prefix"):
        open_index(gtf_path, summary_path)


def test_two_files_from_different_releases_are_refused(tmp_path: Path) -> None:
    """A v1.3 summary against a v1.5 GTF drops genes rather than failing."""
    gtf_path, summary_path = minimal_release(tmp_path, gtf_prefix="MANE.GRCh38.v1.3.")
    with pytest.raises(AdapterUnavailableError, match="release mismatch"):
        open_index(gtf_path, summary_path)


def test_a_renamed_release_is_caught_by_the_name_inside_the_bytes(tmp_path: Path) -> None:
    gtf_path, summary_path = minimal_release(
        tmp_path, prefix="MANE.GRCh38.v9.9.", stored_prefix="MANE.GRCh38.v1.5."
    )
    with pytest.raises(AdapterUnavailableError, match="was renamed"):
        open_index(gtf_path, summary_path)


def test_a_grch37_release_is_refused_for_a_grch38_index(tmp_path: Path) -> None:
    gtf_path, summary_path = minimal_release(tmp_path, prefix="MANE.GRCh37.v1.5.")
    with pytest.raises(AdapterUnavailableError, match="GRCh37"):
        open_index(gtf_path, summary_path)


def test_a_summary_missing_a_required_column_is_refused(tmp_path: Path) -> None:
    gtf_path = tmp_path / "MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz"
    summary_path = tmp_path / "MANE.GRCh38.v1.5.summary.txt.gz"
    _write_gz(
        gtf_path,
        GTF_TEMPLATE.format(
            contig="chr15", start=100, end=200, strand="+", gene_id="ENSG1.1", symbol="AAA"
        ),
        None,
        "ensembl_genomic.gtf",
    )
    _write_gz(summary_path, "#Ensembl_Gene\tsymbol\n", None, "summary.txt")
    with pytest.raises(AdapterUnavailableError, match="missing MANE summary column"):
        open_index(gtf_path, summary_path)


def test_a_gtf_row_with_the_wrong_column_count_is_refused(tmp_path: Path) -> None:
    gtf_path = tmp_path / "MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz"
    summary_path = tmp_path / "MANE.GRCh38.v1.5.summary.txt.gz"
    _write_gz(gtf_path, "chr15\tsource\tgene\t100\t200\n", None, "ensembl_genomic.gtf")
    _write_gz(summary_path, SUMMARY_HEADER, None, "summary.txt")
    with pytest.raises(AdapterUnavailableError, match="tab-separated"):
        open_index(gtf_path, summary_path)


def test_a_gene_with_no_summary_row_is_refused(tmp_path: Path) -> None:
    """The gene would otherwise have no stable identity to join on downstream."""
    gtf_path, summary_path = write_release(
        tmp_path,
        gtf_rows=[
            GTF_TEMPLATE.format(
                contig="chr15", start=100, end=200, strand="+", gene_id="ENSG1.1", symbol="AAA"
            ),
            GTF_TEMPLATE.format(
                contig="chr15", start=300, end=400, strand="+", gene_id="ENSG2.1", symbol="BBB"
            ),
        ],
        summary_rows=[
            SUMMARY_TEMPLATE.format(
                gene_id="ENSG1.1",
                symbol="AAA",
                n=1,
                status=MANE_SELECT_STATUS,
                accession="NC_000015.10",
                start=100,
                end=200,
                strand="+",
            )
        ],
    )
    with pytest.raises(AdapterUnavailableError, match="no row in"):
        open_index(gtf_path, summary_path)


def test_a_release_with_no_primary_contig_genes_is_refused(tmp_path: Path) -> None:
    """An empty gene model reads exactly like a genome of intergenic variants."""
    gtf_path, summary_path = write_release(
        tmp_path,
        gtf_rows=[
            GTF_TEMPLATE.format(
                contig="chr15_KI270905v1_alt",
                start=100,
                end=200,
                strand="+",
                gene_id="ENSG1.1",
                symbol="AAA",
            )
        ],
        summary_rows=[
            SUMMARY_TEMPLATE.format(
                gene_id="ENSG1.1",
                symbol="AAA",
                n=1,
                status=MANE_SELECT_STATUS,
                accession="NT_187633.1",
                start=100,
                end=200,
                strand="+",
            )
        ],
    )
    with pytest.raises(AdapterUnavailableError, match="no gene rows"):
        open_index(gtf_path, summary_path)


def test_one_gene_on_two_primary_contigs_is_refused(tmp_path: Path) -> None:
    gtf_path, summary_path = write_release(
        tmp_path,
        gtf_rows=[
            GTF_TEMPLATE.format(
                contig="chr15", start=100, end=200, strand="+", gene_id="ENSG1.1", symbol="AAA"
            ),
            GTF_TEMPLATE.format(
                contig="chr11", start=100, end=200, strand="+", gene_id="ENSG1.1", symbol="AAA"
            ),
        ],
        summary_rows=[],
    )
    with pytest.raises(AdapterUnavailableError, match="more than one primary contig"):
        open_index(gtf_path, summary_path)


def test_a_malformed_variant_id_is_refused_without_echoing_it(index: ManeGeneIndex) -> None:
    for bad in ("chr15:40160984:C:T", "GRCh38:chr15:40160984:C:T:extra"):
        with pytest.raises(ValueError, match="colon-separated fields") as excinfo:
            index.assign_variants([bad])
        assert "40160984" not in str(excinfo.value)


def test_a_variant_id_with_an_unparseable_build_is_refused(index: ManeGeneIndex) -> None:
    with pytest.raises(ValueError, match="recognised genome build"):
        index.assign_variants(["T2T:chr15:40160984:C:T"])


def test_a_variant_id_with_a_non_numeric_position_is_refused(index: ManeGeneIndex) -> None:
    with pytest.raises(ValueError, match="positive integer"):
        index.assign_variants(["GRCh38:chr15:forty:C:T"])


# --------------------------------------------------------------------------- the blocker


def test_the_remediation_message_names_every_file_it_touches() -> None:
    """A blocker whose message is not actionable is just a crash.

    Three files consume ``ConsequenceAnnotation.impact`` without a None branch:
    ``models/variant.py`` (the field and ``worst_impact_for_gene``),
    ``annotation/service.py::_consequence_evidence`` (the claim, the payload, the
    direction and the strength) and ``evidence/store.py`` (``csq.impact.value``).
    All three are named, so the fix does not have to be rediscovered.
    """
    for path in (
        "src/mva/models/variant.py",
        "src/mva/annotation/service.py",
        "src/mva/evidence/store.py",
    ):
        assert path in IMPACT_NOT_ASSESSED_REMEDIATION
    assert "ImpactSeverity | None" in IMPACT_NOT_ASSESSED_REMEDIATION
    assert "worst_impact_for_gene" in IMPACT_NOT_ASSESSED_REMEDIATION
    # It must also say why no enum value would have done instead.
    assert "benign_consequence" in IMPACT_NOT_ASSESSED_REMEDIATION
    assert "GP-14" in IMPACT_NOT_ASSESSED_REMEDIATION


def test_the_adapter_refuses_rather_than_fabricating_an_impact(index: ManeGeneIndex) -> None:
    """An interval join computes no impact, and every enum value would be a claim.

    MODIFIER and LOW are both members of ``prioritization.filters.BENIGN_IMPACTS``,
    so either would attach a ``benign_consequence`` flag to a variant nobody
    assessed — absence of information rendered as negative information (GP-14).
    MODERATE and HIGH invent severity instead, and HIGH additionally promotes every
    lone heterozygote in ``pairing._wants_single_candidate``.

    Passes in both worlds on purpose: it refuses while the model cannot hold the
    truth, and emits the truth once it can.
    """
    if MODEL_CAN_REPRESENT_UNASSESSED_IMPACT:
        identifier = variant_id("chr15", BUB1B.start)
        assert ManeGeneAdapter(index).annotate([identifier])[identifier][0].impact is None
        return
    with pytest.raises(AdapterUnavailableError) as excinfo:
        ManeGeneAdapter(index)
    assert str(excinfo.value) == IMPACT_NOT_ASSESSED_REMEDIATION


def test_the_model_probe_agrees_with_what_the_model_actually_does(
    index: ManeGeneIndex,
) -> None:
    """The probe gates 19 tests and the adapter's refusal, so it is checked directly.

    Construction, not pydantic introspection: the moment
    ``ConsequenceAnnotation.impact`` is widened, the probe must notice on its own
    and the adapter must start working with no edit to this module.
    """
    payload = {
        "gene_symbol": "BUB1B",
        "transcript_id": "ENST00000287598.11",
        "consequence_terms": [GENE_LOCUS_TERM],
        "impact": None,
    }
    if MODEL_CAN_REPRESENT_UNASSESSED_IMPACT:
        assert ConsequenceAnnotation.model_validate(payload).impact is None
        assert ManeGeneAdapter(index).version == EXPECTED_VERSION
    else:
        with pytest.raises(ValidationError):
            ConsequenceAnnotation.model_validate(payload)


def test_placeholder_impact_cannot_be_reported_as_real(index: ManeGeneIndex) -> None:
    """There is no configuration that makes this adapter emit an invented impact.

    An earlier draft offered ``unassessed_impact=`` as a last-resort override,
    stamped into the version string. That was the fabrication with a switch in
    front of it: ``ManeGeneAdapter(index, unassessed_impact=HIGH)`` would have
    turned every gene overlap into a HIGH-impact consequence — promoting lone
    heterozygotes in ``pairing._wants_single_candidate`` and scoring 0.90 in
    ``scoring._IMPACT_BASE`` — while ``synthetic`` and ``is_synthetic()`` both
    still returned False. ``is_synthetic`` fails closed precisely so a mock cannot
    present itself as real, and a hand-rolled opt-out defeats that mechanism.

    The knob is gone. This test holds it gone: no keyword on this constructor may
    supply an impact severity.
    """
    import inspect

    parameters = inspect.signature(ManeGeneAdapter.__init__).parameters
    assert "unassessed_impact" not in parameters
    assert not any("impact" in name for name in parameters)

    if MODEL_CAN_REPRESENT_UNASSESSED_IMPACT:
        # The model can hold the truth, so the adapter emits it and stays real.
        built = ManeGeneAdapter(index)
        assert built.synthetic is False
        assert is_synthetic(built) is False
        assert built.version == EXPECTED_VERSION
        identifier = variant_id("chr15", BUB1B.start)
        assert built.annotate([identifier])[identifier][0].impact is None
    else:
        # The model cannot, so the adapter does not exist rather than lying.
        with pytest.raises(AdapterUnavailableError):
            ManeGeneAdapter(index)


@requires_nullable_impact
def test_no_annotation_this_adapter_emits_carries_an_impact(index: ManeGeneIndex) -> None:
    """Every field it fills is a MANE fact; the one it cannot compute stays empty."""
    built = ManeGeneAdapter(index)
    ids = [
        variant_id("chr15", BUB1B.start),
        variant_id("chr11", CEP57.start),
        variant_id("chrX", MECP2.start),
    ]
    annotations = [a for group in built.annotate(ids).values() for a in group]
    assert annotations
    assert all(a.impact is None for a in annotations)
    assert all(a.hgvs_c is None and a.hgvs_p is None for a in annotations)
    assert all(a.splice_ai_delta_max is None for a in annotations)
    assert all(a.pathogenicity_scores == {} for a in annotations)


@requires_nullable_impact
def test_an_unassessed_impact_never_becomes_a_benign_call(index: ManeGeneIndex) -> None:
    """The reason no enum value would do: LOW and MODIFIER both flag benign."""
    from mva.config import FrequencyThresholds, QualityThresholds
    from mva.prioritization.filters import BENIGN_IMPACTS, FLAG_BENIGN_CONSEQUENCE, apply_soft_flags

    assert ImpactSeverity.MODIFIER in BENIGN_IMPACTS
    assert ImpactSeverity.LOW in BENIGN_IMPACTS

    identifier = variant_id("chr15", BUB1B.start)
    annotations = ManeGeneAdapter(index).annotate([identifier])[identifier]
    record = make_record("chr15", BUB1B.start, annotations)
    flagged = apply_soft_flags(
        [record], frequency=FrequencyThresholds(), quality=QualityThresholds()
    )
    assert FLAG_BENIGN_CONSEQUENCE not in flagged[0].qc_flags
    # ... and the gene is known while its consequence honestly is not.
    assert flagged[0].gene_symbols == ("BUB1B",)
    assert flagged[0].worst_impact_for_gene("BUB1B") is None


# --------------------------------------------------------------------------- composition


def test_backfill_merge_logic_runs_without_either_real_adapter() -> None:
    """The composition itself, independent of what fills the two slots."""
    covered = variant_id("chr15", BUB1B.start)
    uncovered = variant_id("chr11", CEP57.start)
    combined = GeneBackfillConsequenceAdapter(
        _StubConsequenceAdapter({covered: "PRIMARY_GENE"}),
        _StubConsequenceAdapter({covered: "NEVER_USED", uncovered: "FALLBACK_GENE"}),
    )
    result = combined.annotate([covered, uncovered, variant_id("chr5", TRIP13.start)])
    assert [a.gene_symbol for a in result[covered]] == ["PRIMARY_GENE"]
    assert [a.gene_symbol for a in result[uncovered]] == ["FALLBACK_GENE"]
    assert len(result) == 2  # the third variant is known to neither, so it is absent
    assert combined.name == "stub-predictor+stub-predictor"


@requires_nullable_impact
def test_backfill_uses_the_primary_answer_where_it_has_one(
    adapter: ManeGeneAdapter,
) -> None:
    primary = _StubConsequenceAdapter({variant_id("chr15", BUB1B.start): "SNPEFF_GENE"})
    combined = GeneBackfillConsequenceAdapter(primary, adapter)
    identifier = variant_id("chr15", BUB1B.start)
    assert [a.gene_symbol for a in combined.annotate([identifier])[identifier]] == ["SNPEFF_GENE"]


@requires_nullable_impact
def test_backfill_fills_only_the_variants_the_primary_omitted(
    adapter: ManeGeneAdapter,
) -> None:
    covered = variant_id("chr15", BUB1B.start)
    uncovered = variant_id("chr11", CEP57.start)
    combined = GeneBackfillConsequenceAdapter(
        _StubConsequenceAdapter({covered: "SNPEFF_GENE"}), adapter
    )
    result = combined.annotate([covered, uncovered])
    assert [a.gene_symbol for a in result[covered]] == ["SNPEFF_GENE"]
    assert [a.gene_symbol for a in result[uncovered]] == ["FAM76B", "CEP57"]


@requires_nullable_impact
def test_backfill_leaves_a_variant_neither_adapter_knows_absent(
    adapter: ManeGeneAdapter,
) -> None:
    intergenic = variant_id("chr15", INTERGENIC_CHR15)
    combined = GeneBackfillConsequenceAdapter(_StubConsequenceAdapter({}), adapter)
    assert combined.annotate([intergenic]) == {}


@requires_nullable_impact
def test_backfill_maturity_fails_closed(adapter: ManeGeneAdapter) -> None:
    """A real predictor plus an undeclared partner is still labelled synthetic."""
    combined = GeneBackfillConsequenceAdapter(_StubConsequenceAdapter({}), adapter)
    assert combined.synthetic is True  # the stub declares nothing
    assert is_synthetic(combined) is True
    assert combined.name.endswith(MANE_ADAPTER_NAME)
    assert EXPECTED_VERSION in combined.version


class _StubConsequenceAdapter:
    """A stand-in primary predictor. Deliberately declares no ``synthetic``."""

    def __init__(self, genes_by_variant: Mapping[str, str]) -> None:
        self._genes = genes_by_variant

    @property
    def name(self) -> str:
        return "stub-predictor"

    @property
    def version(self) -> str:
        return "v0"

    def annotate(
        self, variant_ids: Sequence[str]
    ) -> Mapping[str, tuple[ConsequenceAnnotation, ...]]:
        return {
            identifier: (
                ConsequenceAnnotation(
                    gene_symbol=self._genes[identifier],
                    transcript_id="ENST_STUB",
                    consequence_terms=("missense_variant",),
                    impact=ImpactSeverity.MODERATE,
                ),
            )
            for identifier in variant_ids
            if identifier in self._genes
        }


# --------------------------------------------------------------------------- the payoff


def test_generate_pairs_returns_nothing_without_gene_assignment() -> None:
    """The defect this module exists to fix, reproduced first."""
    starved = [make_record("chr15", 40_200_000), make_record("chr15", 40_210_000)]
    assert all(record.gene_symbols == () for record in starved)
    assert generate_pairs(starved) == ()


def test_the_index_assigns_both_bub1b_variants_that_generate_pairs_needs(
    index: ManeGeneIndex,
) -> None:
    """The join that unstarves pairing, asserted without the Protocol shim."""
    positions = (40_200_000, 40_210_000)
    ids = [variant_id("chr15", position) for position in positions]
    assigned = index.assign_variants(ids)
    assert [[a.gene.gene_symbol for a in assigned[i]] for i in ids] == [["BUB1B"], ["BUB1B"]]


def test_a_bub1b_compound_het_appears_as_soon_as_the_gene_is_attached(
    index: ManeGeneIndex,
) -> None:
    """gene -> gene_symbols -> _variants_by_gene -> generate_pairs, end to end.

    The consequence attached here comes from a stand-in *effect predictor*, not
    from MANE: the impact is that predictor's claim, and the gene symbol is the
    one this module's interval join resolved. That is exactly the composition
    ``GeneBackfillConsequenceAdapter`` produces, and it keeps this test from
    depending on an impact severity MANE never computed.
    """
    positions = (40_200_000, 40_210_000)
    ids = [variant_id("chr15", position) for position in positions]
    assigned = index.assign_variants(ids)
    records = [
        make_record(
            "chr15",
            position,
            tuple(
                ConsequenceAnnotation(
                    gene_symbol=assignment.gene.gene_symbol,
                    gene_id=assignment.gene.gene_id,
                    transcript_id=assignment.gene.mane_select.transcript_id,
                    consequence_terms=("missense_variant",),
                    impact=ImpactSeverity.MODERATE,
                )
                for assignment in assigned[identifier]
            ),
        )
        for position, identifier in zip(positions, ids, strict=True)
    ]
    assert [record.gene_symbols for record in records] == [("BUB1B",), ("BUB1B",)]

    candidates = generate_pairs(records)
    assert candidates != ()
    pairs = [c for c in candidates if c.is_pair]
    assert [c.gene_symbol for c in pairs] == ["BUB1B"]
    assert pairs[0].variant_ids == tuple(ids)


def test_a_variant_in_two_genes_would_pair_under_each(index: ManeGeneIndex) -> None:
    """The multi-gene guarantee, carried into the hypothesis space itself."""
    sites = (CEP57.start, MTMR2.start - 1, FAM76B.start)
    ids = [variant_id("chr11", position) for position in sites]
    assigned = index.assign_variants(ids)
    records = [
        make_record(
            "chr11",
            position,
            tuple(
                ConsequenceAnnotation(
                    gene_symbol=assignment.gene.gene_symbol,
                    transcript_id=assignment.gene.mane_select.transcript_id,
                    consequence_terms=("missense_variant",),
                    impact=ImpactSeverity.MODERATE,
                )
                for assignment in assigned[identifier]
            ),
        )
        for position, identifier in zip(sites, ids, strict=True)
    ]
    assert records[0].gene_symbols == ("FAM76B", "CEP57")
    genes = {c.gene_symbol for c in generate_pairs(records) if c.is_pair}
    assert genes == {"CEP57", "FAM76B"}


# ----------------------------------------------------------------- real proband VCF
#
# The proband callset is GATK output with BARE contigs (`15`, not `chr15`), a
# full decoy-and-HLA contig set in its header, and an indel-heavy record mix — a
# large minority of records are indel-bearing.
# These four tests are that shape, not the tidy one.


def test_a_bare_contig_variant_id_resolves_exactly_as_a_prefixed_one(
    index: ManeGeneIndex,
) -> None:
    """The proband VCF says ``15``; the canonical key says ``chr15``."""
    bare = f"GRCh38:15:{BUB1B.start}:C:T"
    prefixed = variant_id("chr15", BUB1B.start)
    bare_result = index.assign_variants([bare])
    prefixed_result = index.assign_variants([prefixed])
    assert list(bare_result) == [bare]
    assert list(prefixed_result) == [prefixed]
    assert [a.gene.gene_id for a in bare_result[bare]] == [
        a.gene.gene_id for a in prefixed_result[prefixed]
    ]


def test_the_result_is_keyed_by_the_callers_own_contig_spelling(
    index: ManeGeneIndex,
) -> None:
    """Rekeying to a normalised form would make every caller lookup miss."""
    bare = f"GRCh38:X:{MECP2.start}:C:T"
    assigned = index.assign_variants([bare])
    assert list(assigned) == [bare]
    assert [a.gene.gene_symbol for a in assigned[bare]] == ["MECP2"]


def test_the_mitochondrial_contig_resolves_from_both_spellings() -> None:
    """``chrM`` and ``MT`` are the same sequence; MANE v1.5 places no gene on it."""
    assert contig_for_refseq_accession("NC_012920.1") == "chrM"


def test_a_decoy_or_hla_contig_is_refused_rather_than_silently_missed(
    index: ManeGeneIndex,
) -> None:
    """Thousands of contigs in the callset header; this pipeline reasons about 25.

    Refusing is the right failure: a decoy that silently returned no gene would be
    indistinguishable from a real intergenic call. Ingestion already cannot build a
    ``GenomicCoordinate`` on such a contig, so this path is a backstop.
    """
    for contig in ("chrUn_KI270302v1", "chr15_KI270905v1_alt", "chrEBV"):
        with pytest.raises(ValueError, match="Non-canonical contig"):
            index.assign_variants([f"GRCh38:{contig}:1000:C:T"])


def test_an_indel_bearing_callset_joins_on_its_ref_span(index: ManeGeneIndex) -> None:
    """Nearly one record in five is an indel, so the span path is the common path."""
    deletion = f"GRCh38:15:{BUB1B.end - 2}:CATG:C"
    insertion = f"GRCh38:15:{BUB1B.end}:C:CATG"
    past_the_end = f"GRCh38:15:{BUB1B.end + 1}:CATG:C"
    assigned = index.assign_variants([deletion, insertion, past_the_end])
    assert [a.gene.gene_symbol for a in assigned[deletion]] == ["BUB1B", "PAK6"]
    assert [a.gene.gene_symbol for a in assigned[insertion]] == ["BUB1B", "PAK6"]
    assert [a.gene.gene_symbol for a in assigned[past_the_end]] == ["PAK6"]


# --------------------------------------------------------------------------- performance


def test_lookups_do_not_scale_with_the_number_of_genes(index: ManeGeneIndex) -> None:
    """A guard on the algorithm, not a benchmark.

    A linear scan would touch every gene on the contig for every query. The
    binary search plus bounded backward walk touches a handful, so 20,000 lookups
    spread across the fixture's contigs stay well inside a second even on a slow
    machine. The measured figure on the full 19,299-gene release is ~868,000
    lookups per second, reported alongside this module.
    """
    import time

    sites = [("chr15", BUB1B.start + (step % 60_000)) for step in range(0, 20_000 * 7, 7)]
    started = time.perf_counter()
    hits = sum(1 for contig, position in sites if index.genes_at(contig, position))
    elapsed = time.perf_counter() - started
    assert hits > 0
    assert elapsed < 5.0, f"20,000 lookups took {elapsed:.2f}s; the index is not being used"


def test_the_index_holds_no_module_level_state() -> None:
    """Two indexes over the same release are independent objects, not one cache."""
    first, second = open_index(), open_index()
    assert first is not second
    assert first.genes is not second.genes
    assert [g.gene_id for g in first.genes] == [g.gene_id for g in second.genes]
