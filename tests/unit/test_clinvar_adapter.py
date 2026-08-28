"""Unit tests for the real ClinVar VCF adapter.

Everything here runs against ``tests/fixtures/clinvar/clinvar_slice.vcf.gz`` — a
1425-record slice cut out of the genuine NCBI ClinVar GRCh38 release of
2026-08-22 by ``tests/fixtures/clinvar/make_fixture.py``, which documents the
exact regions and why each was chosen. The records are real, so the CLNSIG
spellings, the review statuses, the underscore encoding and the presence of
oncogenicity-only records are all the source's own, not something written to
match the parser.

The interesting failures in a clinical adapter are not "wrong string returned".
They are:

* a variant ClinVar has never seen coming back as an empty-but-present entry, or
  worse as a benign call (GP-14);
* a 0-star single unreviewed submission and a 3-star expert-panel review being
  rendered as the same "Pathogenic";
* a conflicting classification being quietly resolved to its most severe member;
* ``chr15`` being compared raw against ClinVar's ``15``, which fails by finding
  nothing and is therefore indistinguishable from "no record".

Each of those has a test below.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import re
import traceback
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import pysam
import pytest

from mva.alleles import (
    LeftAlignmentStatus,
    ReferenceLookup,
    canonicalise_allele,
    is_sequence_allele,
    trim_parsimoniously,
)
from mva.annotation import (
    SYNTHETIC_STANDIN_LIMITATION,
    AdapterRole,
    AdapterSet,
    ClinicalAdapter,
    LocalConsequenceAdapter,
    LocalFrequencyAdapter,
    annotate_variants,
    is_synthetic,
    load_default_adapters,
)
from mva.annotation.clinvar_vcf import (
    CLINVAR_ADAPTER_NAME,
    CLINVAR_STAR_RATINGS,
    ClinvarVcfAdapter,
    merge_query_regions,
    merge_query_spans,
    read_shipped_md5,
)
from mva.clock import FixedClock
from mva.determinism import hash_file, stable_hash
from mva.errors import AdapterUnavailableError, GenomeBuildMismatchError
from mva.ingestion.normalise import open_reference_fasta
from mva.models.base import AssertionTier
from mva.models.evidence import EvidenceDirection, EvidenceStrength
from mva.models.genome import GenomeBuild, GenomicCoordinate
from mva.models.variant import (
    ClinicalAssertion,
    FilterStatus,
    Genotype,
    VariantRecord,
    Zygosity,
)

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE = REPO_ROOT / "tests" / "fixtures" / "clinvar" / "clinvar_slice.vcf.gz"
KNOWLEDGE_ROOT = REPO_ROOT / "knowledge"

#: The release the fixture was cut from, read out of its own ``##`` headers.
EXPECTED_VERSION = "GRCh38-2026-08-22"

# --------------------------------------------------------------------------- anchors
#
# Real records from the fixture, each chosen for the trap it exercises. Coordinates
# are public ClinVar reference data, not patient data.

#: CFTR, ``practice_guideline`` — the top of the review-depth scale (4 stars).
PRACTICE_GUIDELINE = "GRCh38:chr7:117509123:G:A"

#: CFTR, ``reviewed_by_expert_panel`` Pathogenic — 3 stars.
EXPERT_PANEL_PATHOGENIC = "GRCh38:chr7:117509033:G:A"

#: BUB1B, ``no_assertion_criteria_provided`` Pathogenic — 0 stars. Same leading
#: significance word as the record above and a completely different evidential
#: weight; this pair is the whole of "review status is not significance".
UNREVIEWED_PATHOGENIC = "GRCh38:chr15:40200239:A:G"

#: CFTR, ``Conflicting_classifications_of_pathogenicity`` with a CLNSIGCONF
#: breakdown of ``Likely_pathogenic_(1)|Uncertain_significance_(3)``.
CONFLICTING = "GRCh38:chr7:117509053:G:C"

#: BUB1B, Likely benign, a single clean condition — the decoding anchor.
BUB1B_LIKELY_BENIGN = "GRCh38:chr15:40199604:A:C"

#: CFTR one-base deletion, Pathogenic.
INDEL_PATHOGENIC = "GRCh38:chr7:117509031:CA:C"

#: BRCA1 3'UTR. The record exists in ClinVar but carries only an *oncogenicity*
#: classification (ONC) and no germline CLNSIG at all.
ONCOGENICITY_ONLY = "GRCh38:chr17:43045670:G:A"

#: Inside a region the fixture covers, but no ClinVar record at this allele.
NOT_IN_CLINVAR = "GRCh38:chr15:40199604:A:T"

#: Same locus as ``NOT_IN_CLINVAR`` — proves the miss is the allele, not the region.
ANOTHER_ABSENT = "GRCh38:chr15:40201234:G:A"

MINIMAL_HEADER = (
    "##fileformat=VCFv4.1\n"
    "##fileDate=2026-08-22\n"
    "##source=ClinVar\n"
    "##reference=GRCh38\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


# --------------------------------------------------------------------------- helpers


def write_indexed_vcf(directory: Path, text: str, *, name: str = "mini.vcf") -> Path:
    """bgzip + tabix a VCF built in-test, returning the compressed path.

    Used only for shapes the real release does not contain (multi-ALT, percent
    escapes, two records at one allele) or must not contain (a GRCh37 header).
    Everything the release *does* contain is tested against the release itself.
    """
    plain = directory / name
    plain.write_text(text, encoding="utf-8")
    compressed = directory / f"{name}.gz"
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)
    return compressed


def open_adapter(path: Path, **kwargs: object) -> ClinvarVcfAdapter:
    """Construct an adapter pinned to the file's actual sha256."""
    return ClinvarVcfAdapter(path, expected_sha256=hash_file(path), **kwargs)  # type: ignore[arg-type]


@pytest.fixture(scope="module")
def adapter() -> Iterator[ClinvarVcfAdapter]:
    instance = open_adapter(FIXTURE)
    yield instance
    instance.close()


def fixture_records() -> list[list[str]]:
    """Every record in the committed slice, as raw column lists."""
    rows: list[list[str]] = []
    with gzip.open(FIXTURE, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                rows.append(line.rstrip("\n").split("\t"))
    return rows


def info_of(columns: list[str]) -> dict[str, str]:
    fields: dict[str, str] = {}
    for part in columns[7].split(";"):
        key, separator, value = part.partition("=")
        fields[key] = value if separator else ""
    return fields


def make_record(variant_id: str) -> VariantRecord:
    """A minimal proband-shaped record at a canonical variant ID."""
    build, contig, position, ref, alt = variant_id.split(":")
    return VariantRecord(
        coordinate=GenomicCoordinate(
            build=GenomeBuild(build), contig=contig, position=int(position), ref=ref, alt=alt
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
        source_artifact="tests/fixtures/clinvar/clinvar_slice.vcf.gz",
    )


def dump(result: Mapping[str, tuple[ClinicalAssertion, ...]]) -> str:
    """Canonical hash of an assertions mapping, for byte-identity comparison."""
    return stable_hash(
        {key: [item.model_dump(mode="json") for item in value] for key, value in result.items()}
    )


# --------------------------------------------------------------------------- identity


def test_fixture_is_the_real_release_slice() -> None:
    """Guard the premise of every other test: this is genuine ClinVar."""
    with gzip.open(FIXTURE, "rt", encoding="utf-8") as handle:
        header = [line for line in handle if line.startswith("##")]
    assert "##source=ClinVar\n" in header
    assert "##reference=GRCh38\n" in header
    assert "##fileDate=2026-08-22\n" in header
    assert len(fixture_records()) == 1425


def test_version_is_read_from_the_release_header(adapter: ClinvarVcfAdapter) -> None:
    """The release string is taken from the file, never supplied by the caller."""
    assert adapter.version == EXPECTED_VERSION
    assert adapter.name == CLINVAR_ADAPTER_NAME


def test_adapter_declares_itself_real(adapter: ClinvarVcfAdapter) -> None:
    """GP-20: ``is_synthetic`` fails closed, so ``synthetic = False`` is deliberate."""
    assert adapter.synthetic is False
    assert is_synthetic(adapter) is False


def test_satisfies_the_clinical_adapter_protocol(adapter: ClinvarVcfAdapter) -> None:
    assert isinstance(adapter, ClinicalAdapter)


def test_adapter_set_reports_the_clinical_slot_as_real(adapter: ClinvarVcfAdapter) -> None:
    """The descriptor the run manifest and report footers are built from."""
    adapters = AdapterSet(
        consequence=LocalConsequenceAdapter(
            KNOWLEDGE_ROOT / "public" / "consequences.tsv", version="synthetic-v0.0"
        ),
        frequency=LocalFrequencyAdapter(
            KNOWLEDGE_ROOT / "public" / "frequencies.tsv", version="synthetic-v0.0"
        ),
        clinical=adapter,
    )
    clinical = next(d for d in adapters.descriptors() if d.role is AdapterRole.CLINICAL)
    assert clinical.synthetic is False
    assert clinical.label == f"{CLINVAR_ADAPTER_NAME}@{EXPECTED_VERSION}"


# --------------------------------------------------------------------- contig mapping


def test_contig_map_is_resolved_from_the_index_not_assumed(adapter: ClinvarVcfAdapter) -> None:
    """ClinVar is Ensembl-named; the pipeline's join key is UCSC-named.

    Asserted directly rather than inferred from a successful lookup: a wrong map
    fails by finding nothing, which reads exactly like "ClinVar has no record".
    """
    mapping = adapter.contig_map
    assert mapping["chr15"] == "15"
    assert mapping["chr7"] == "7"
    assert mapping["chr17"] == "17"
    # The full release also carries X/Y/MT; the slice's index only holds 7/15/17,
    # so absence here proves the map is read from the index rather than fabricated.
    assert set(mapping) == {"chr7", "chr15", "chr17"}


def test_full_release_maps_sex_and_mitochondrial_contigs(tmp_path: Path) -> None:
    """``chrX -> X``, ``chrY -> Y`` and — the easy one to get wrong — ``chrM -> MT``."""
    path = write_indexed_vcf(
        tmp_path,
        MINIMAL_HEADER
        + "X\t100\t1\tA\tG\t.\t.\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter\n"
        + "Y\t100\t2\tA\tG\t.\t.\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter\n"
        + "MT\t100\t3\tA\tG\t.\t.\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter\n",
    )
    instance = open_adapter(path)
    try:
        assert instance.contig_map["chrX"] == "X"
        assert instance.contig_map["chrY"] == "Y"
        assert instance.contig_map["chrM"] == "MT"
        assert set(instance.assertions(["GRCh38:chrM:100:A:G"])) == {"GRCh38:chrM:100:A:G"}
    finally:
        instance.close()


def test_ucsc_variant_id_joins_to_a_bare_clinvar_contig(adapter: ClinvarVcfAdapter) -> None:
    """The headline case: a ``chr15`` key finds a record stored under ``15``."""
    result = adapter.assertions([BUB1B_LIKELY_BENIGN])
    assert BUB1B_LIKELY_BENIGN in result
    assert result[BUB1B_LIKELY_BENIGN][0].significance == "Likely benign"


def test_result_is_keyed_by_the_caller_s_own_spelling(adapter: ClinvarVcfAdapter) -> None:
    """A caller who writes the contig bare gets its key back, not a rewritten one."""
    bare = "GRCh38:15:40199604:A:C"
    result = adapter.assertions([bare])
    assert set(result) == {bare}
    assert result[bare] == adapter.assertions([BUB1B_LIKELY_BENIGN])[BUB1B_LIKELY_BENIGN]


# ----------------------------------------------------------- GP-14: absence is absence


def test_variants_clinvar_has_never_seen_are_omitted_entirely(
    adapter: ClinvarVcfAdapter,
) -> None:
    """Not an empty tuple, not a benign call — no key at all (GP-14).

    An empty-but-present entry says "ClinVar looked and found nothing to say",
    which downstream is one short step from "nothing pathogenic is known here".
    Absence of a curated assertion is absence of evidence in both directions.
    """
    result = adapter.assertions([NOT_IN_CLINVAR, ANOTHER_ABSENT, BUB1B_LIKELY_BENIGN])
    assert NOT_IN_CLINVAR not in result
    assert ANOTHER_ABSENT not in result
    assert set(result) == {BUB1B_LIKELY_BENIGN}
    assert all(value for value in result.values()), "no key may map to an empty tuple"


def test_an_entirely_unknown_batch_returns_an_empty_mapping(adapter: ClinvarVcfAdapter) -> None:
    assert adapter.assertions([NOT_IN_CLINVAR, ANOTHER_ABSENT]) == {}


def test_no_absent_variant_is_given_a_fabricated_benign_call(
    adapter: ClinvarVcfAdapter,
) -> None:
    """The failure this test exists for would look like a helpful default."""
    for variant_id, assertions in adapter.assertions([NOT_IN_CLINVAR]).items():
        pytest.fail(f"absent variant {variant_id} was given {len(assertions)} assertion(s)")


def test_oncogenicity_only_records_yield_no_germline_assertion(
    adapter: ClinvarVcfAdapter,
) -> None:
    """A record with ONC but no CLNSIG has no *germline* classification.

    ClinVar classifies germline pathogenicity, oncogenicity and somatic clinical
    impact on three separate axes. ``ClinicalAssertion`` models the first.
    Reading an oncogenicity call as a germline one would attribute to ClinVar a
    statement it did not make; defaulting it to benign would be worse.
    """
    columns = next(row for row in fixture_records() if row[1] == "43045670")
    assert "CLNSIG" not in info_of(columns), "fixture drifted: this record gained a CLNSIG"
    assert "ONC" in info_of(columns)
    assert adapter.assertions([ONCOGENICITY_ONLY]) == {}


def test_the_slice_really_does_contain_such_records() -> None:
    """Guard the test above against a fixture that silently stopped covering it."""
    germline_absent = sum(1 for row in fixture_records() if "CLNSIG" not in info_of(row))
    assert germline_absent > 100


# ------------------------------------------------- review status is not significance


def test_star_rating_survives_into_the_model(adapter: ClinvarVcfAdapter) -> None:
    """Two real Pathogenic calls, 0 stars and 3 stars, must not read alike.

    ``annotation.service._significance_strength`` grades evidence on the star
    rating alone, so flattening it here would silently promote one unreviewed
    submission to expert-panel weight.
    """
    result = adapter.assertions([UNREVIEWED_PATHOGENIC, EXPERT_PANEL_PATHOGENIC])
    unreviewed = result[UNREVIEWED_PATHOGENIC][0]
    expert = result[EXPERT_PANEL_PATHOGENIC][0]

    assert unreviewed.significance.startswith("Pathogenic")
    assert expert.significance.startswith("Pathogenic")

    assert unreviewed.star_rating == 0
    assert unreviewed.review_status == "no assertion criteria provided"
    assert expert.star_rating == 3
    assert expert.review_status == "reviewed by expert panel"


def test_practice_guideline_is_four_stars(adapter: ClinvarVcfAdapter) -> None:
    assertion = adapter.assertions([PRACTICE_GUIDELINE])[PRACTICE_GUIDELINE][0]
    assert assertion.review_status == "practice guideline"
    assert assertion.star_rating == 4
    assert assertion.significance == "Pathogenic"


def test_every_review_status_in_the_release_has_a_star_rating() -> None:
    """The lookup table must cover what the real release actually contains.

    Keyed on the raw VCF token, so this fails loudly if the table is ever written
    against the prettified display text instead.
    """
    seen = {info_of(row)["CLNREVSTAT"] for row in fixture_records() if "CLNREVSTAT" in info_of(row)}
    unmapped = sorted(token for token in seen if token not in CLINVAR_STAR_RATINGS)
    assert not unmapped, f"review statuses with no star rating: {unmapped}"
    assert seen >= {
        "practice_guideline",
        "reviewed_by_expert_panel",
        "no_assertion_criteria_provided",
    }


def test_unknown_review_status_scores_none_not_zero(tmp_path: Path) -> None:
    """A status ClinVar invents next year is ungraded, not zero-star (GP-14).

    ``None`` grades WEAK downstream, which is conservative. A fabricated ``0``
    would be indistinguishable from a genuinely unreviewed submission.
    """
    path = write_indexed_vcf(
        tmp_path,
        MINIMAL_HEADER + "15\t100\t9\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_llm\n",
    )
    instance = open_adapter(path)
    try:
        assertion = instance.assertions(["GRCh38:chr15:100:A:G"])["GRCh38:chr15:100:A:G"][0]
        assert assertion.star_rating is None
        assert assertion.review_status == "reviewed by llm"
    finally:
        instance.close()


# ------------------------------------------------------------------------- conflicts


def test_conflicting_classifications_are_reported_as_conflict(
    adapter: ClinvarVcfAdapter,
) -> None:
    """Never resolved to the most severe or the most common member call.

    The fixture record's CLNSIGCONF is ``Likely_pathogenic_(1)|Uncertain_significance_(3)``.
    A "resolve to most severe" bug returns "Likely pathogenic"; a "resolve to
    majority" bug returns "Uncertain significance". Both are wrong, and both look
    like a cleaner answer than the truth.
    """
    assertion = adapter.assertions([CONFLICTING])[CONFLICTING][0]
    assert assertion.significance == "Conflicting classifications of pathogenicity"
    assert assertion.review_status == "criteria provided, conflicting classifications"
    assert assertion.star_rating == 1
    assert assertion.significance not in {"Likely pathogenic", "Uncertain significance"}


def test_the_annotation_service_files_this_adapter_s_output_correctly(
    adapter: ClinvarVcfAdapter,
) -> None:
    """End-to-end through ``annotate_variants``, with nothing in the service changed.

    Three things are locked at once, and all three are properties of the *strings
    and numbers this adapter emits* rather than of the service:

    * a conflicting call reads NEUTRAL, because the service looks for the
      substring "conflict" before it looks for "pathogenic". Emitting a resolved
      "Likely pathogenic" would silently flip the direction to SUPPORTS.
    * the star rating decides evidence strength, so 3 stars must reach STRONG and
      0 stars must not.
    * ``EvidenceItem`` refuses a DATABASE_ASSERTION whose citation has no version,
      so this also proves the release string survives the whole way.
    """
    base = load_default_adapters(KNOWLEDGE_ROOT, KNOWLEDGE_ROOT / "manifests" / "knowledge.yaml")
    adapters = AdapterSet(consequence=base.consequence, frequency=base.frequency, clinical=adapter)
    records = [
        make_record(CONFLICTING),
        make_record(EXPERT_PANEL_PATHOGENIC),
        make_record(UNREVIEWED_PATHOGENIC),
        make_record(NOT_IN_CLINVAR),
    ]
    result = annotate_variants(
        records, adapters=adapters, clock=FixedClock(datetime(2026, 1, 1, tzinfo=UTC))
    )
    clinical = {
        item.subject_id: item
        for item in result.evidence
        if item.tier is AssertionTier.DATABASE_ASSERTION and item.tool == CLINVAR_ADAPTER_NAME
    }

    assert clinical[CONFLICTING].direction is EvidenceDirection.NEUTRAL
    assert clinical[EXPERT_PANEL_PATHOGENIC].direction is EvidenceDirection.SUPPORTS
    assert clinical[EXPERT_PANEL_PATHOGENIC].strength is EvidenceStrength.STRONG
    assert clinical[UNREVIEWED_PATHOGENIC].strength is not EvidenceStrength.STRONG
    assert NOT_IN_CLINVAR not in clinical

    citation = clinical[EXPERT_PANEL_PATHOGENIC].citation
    assert citation is not None
    assert citation.source == "ClinVar"
    assert citation.version == EXPECTED_VERSION
    assert citation.identifier.startswith("VariationID:")


def test_a_real_adapter_carries_no_synthetic_disclosure(adapter: ClinvarVcfAdapter) -> None:
    """GP-20 cuts both ways: a real source must not be labelled a mock either."""
    base = load_default_adapters(KNOWLEDGE_ROOT, KNOWLEDGE_ROOT / "manifests" / "knowledge.yaml")
    adapters = AdapterSet(consequence=base.consequence, frequency=base.frequency, clinical=adapter)
    result = annotate_variants(
        [make_record(EXPERT_PANEL_PATHOGENIC)],
        adapters=adapters,
        clock=FixedClock(datetime(2026, 1, 1, tzinfo=UTC)),
    )
    item = next(i for i in result.evidence if i.tool == CLINVAR_ADAPTER_NAME)
    assert SYNTHETIC_STANDIN_LIMITATION not in item.limitations
    assert item.limitations, "a real source still has to state its limitations (GP-17)"


def test_multi_value_significance_keeps_clinvar_s_own_separator(
    adapter: ClinvarVcfAdapter,
) -> None:
    """``Pathogenic|Affects`` stays two visible calls rather than collapsing to one."""
    assertion = adapter.assertions([UNREVIEWED_PATHOGENIC])[UNREVIEWED_PATHOGENIC][0]
    assert assertion.significance == "Pathogenic|Affects"


# --------------------------------------------------------------------------- parsing


def test_conditions_are_decoded_from_clinvar_s_vcf_escaping(
    adapter: ClinvarVcfAdapter,
) -> None:
    """``_`` is ClinVar's encoding of a space, so decoding restores its own text."""
    assertion = adapter.assertions([BUB1B_LIKELY_BENIGN])[BUB1B_LIKELY_BENIGN][0]
    assert assertion.conditions == ("Mosaic variegated aneuploidy syndrome 1",)


def test_multiple_conditions_are_all_retained(adapter: ClinvarVcfAdapter) -> None:
    assertion = adapter.assertions([PRACTICE_GUIDELINE])[PRACTICE_GUIDELINE][0]
    assert "Cystic fibrosis" in assertion.conditions
    assert len(assertion.conditions) > 1


def test_percent_escapes_are_decoded(tmp_path: Path) -> None:
    """VCF percent escapes appear in the real release (``%3D``, ``%3B``).

    Decoded with a single-pass scanner rather than ``urllib`` — the whole
    ``urllib`` package is a forbidden import on the patient-data path (PRIV-05).
    """
    path = write_indexed_vcf(
        tmp_path,
        MINIMAL_HEADER + "15\t100\t9\tA\tG\t.\t.\tCLNSIG=Benign;CLNDN=Weight%3D40kg%3Bobese;"
        "CLNREVSTAT=criteria_provided,_single_submitter\n",
    )
    instance = open_adapter(path)
    try:
        assertion = instance.assertions(["GRCh38:chr15:100:A:G"])["GRCh38:chr15:100:A:G"][0]
        assert assertion.conditions == ("Weight=40kg;obese",)
    finally:
        instance.close()


def test_accession_names_the_kind_of_identifier_it_is(adapter: ClinvarVcfAdapter) -> None:
    """The VCF ID column is ClinVar's Variation ID; a bare integer is unresolvable."""
    assertion = adapter.assertions([BUB1B_LIKELY_BENIGN])[BUB1B_LIKELY_BENIGN][0]
    assert assertion.accession == "VariationID:2742393"
    assert assertion.source == "ClinVar"
    assert assertion.version == EXPECTED_VERSION


def test_indels_are_matched_on_the_same_key_as_snvs(adapter: ClinvarVcfAdapter) -> None:
    """ClinVar uses the same anchored-base convention the proband VCF does."""
    assertion = adapter.assertions([INDEL_PATHOGENIC])[INDEL_PATHOGENIC][0]
    assert assertion.significance == "Pathogenic"
    assert assertion.star_rating == 1


def test_null_placeholder_conditions_are_dropped(tmp_path: Path) -> None:
    """ClinVar writes ``.`` where a submitter named no condition."""
    path = write_indexed_vcf(
        tmp_path,
        MINIMAL_HEADER + "15\t100\t9\tA\tG\t.\t.\tCLNSIG=Benign;CLNDN=.|Cystic_fibrosis|.;"
        "CLNREVSTAT=criteria_provided,_single_submitter\n",
    )
    instance = open_adapter(path)
    try:
        assertion = instance.assertions(["GRCh38:chr15:100:A:G"])["GRCh38:chr15:100:A:G"][0]
        assert assertion.conditions == ("Cystic fibrosis",)
    finally:
        instance.close()


# ------------------------------------------------------------------ multi-allelic


def test_multi_alt_records_are_split_per_allele(tmp_path: Path) -> None:
    """VCF permits several ALTs on one line; ClinVar's 2026-08-22 release uses none.

    A full scan of that release found 0 multi-ALT records in 4,467,990, so this
    shape cannot be cut from the real file — it is built here instead. A reader
    that assumed one ALT per line would attach the assertion to ``G`` and drop
    ``T`` entirely, which downstream is indistinguishable from ClinVar having no
    record for the second allele.
    """
    path = write_indexed_vcf(
        tmp_path,
        MINIMAL_HEADER + "15\t100\t9\tA\tG,T\t.\t.\tCLNSIG=Pathogenic;"
        "CLNREVSTAT=reviewed_by_expert_panel;CLNDN=Mosaic_variegated_aneuploidy_syndrome_1\n",
    )
    instance = open_adapter(path)
    try:
        result = instance.assertions(["GRCh38:chr15:100:A:G", "GRCh38:chr15:100:A:T"])
        assert set(result) == {"GRCh38:chr15:100:A:G", "GRCh38:chr15:100:A:T"}
        assert result["GRCh38:chr15:100:A:G"][0] == result["GRCh38:chr15:100:A:T"][0]
        assert result["GRCh38:chr15:100:A:G"][0].star_rating == 3
    finally:
        instance.close()


def test_an_unlisted_alt_at_a_listed_position_is_still_absent(tmp_path: Path) -> None:
    """Position match is not allele match. The join key is all five fields."""
    path = write_indexed_vcf(
        tmp_path,
        MINIMAL_HEADER + "15\t100\t9\tA\tG,T\t.\t.\tCLNSIG=Pathogenic;"
        "CLNREVSTAT=reviewed_by_expert_panel\n",
    )
    instance = open_adapter(path)
    try:
        assert instance.assertions(["GRCh38:chr15:100:A:C"]) == {}
    finally:
        instance.close()


# ------------------------------------------------------------------- determinism


def test_repeat_lookups_are_byte_identical(adapter: ClinvarVcfAdapter) -> None:
    """GP-30. Same inputs, same bytes — including across a reopened adapter."""
    ids = [BUB1B_LIKELY_BENIGN, CONFLICTING, PRACTICE_GUIDELINE, NOT_IN_CLINVAR]
    first = adapter.assertions(ids)
    second = adapter.assertions(ids)
    other = open_adapter(FIXTURE)
    try:
        third = other.assertions(ids)
    finally:
        other.close()
    assert dump(first) == dump(second) == dump(third)


def test_input_order_does_not_change_per_variant_results(adapter: ClinvarVcfAdapter) -> None:
    """Query order is a caller detail; the assertions attached must not depend on it."""
    ids = [BUB1B_LIKELY_BENIGN, CONFLICTING, PRACTICE_GUIDELINE, INDEL_PATHOGENIC]
    forward = adapter.assertions(ids)
    backward = adapter.assertions(list(reversed(ids)))
    assert forward == backward


def test_result_keys_follow_caller_order(adapter: ClinvarVcfAdapter) -> None:
    """Deterministic iteration order, taken from the caller rather than from a set."""
    ids = [PRACTICE_GUIDELINE, BUB1B_LIKELY_BENIGN, NOT_IN_CLINVAR, CONFLICTING]
    assert list(adapter.assertions(ids)) == [PRACTICE_GUIDELINE, BUB1B_LIKELY_BENIGN, CONFLICTING]


def test_duplicate_ids_are_collapsed(adapter: ClinvarVcfAdapter) -> None:
    single = adapter.assertions([BUB1B_LIKELY_BENIGN])
    repeated = adapter.assertions([BUB1B_LIKELY_BENIGN] * 3)
    assert single == repeated


def test_multiple_assertions_at_one_allele_are_ordered_by_accession(tmp_path: Path) -> None:
    """A total order over ties, independent of file order (GP-30).

    Written by hand because the committed slice has exactly one record per
    allele; the ordering path would otherwise never run.
    """
    body = (
        "15\t100\t900\tA\tG\t.\t.\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter\n"
        "15\t100\t100\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel\n"
    )
    path = write_indexed_vcf(tmp_path, MINIMAL_HEADER + body)
    instance = open_adapter(path)
    try:
        assertions = instance.assertions(["GRCh38:chr15:100:A:G"])["GRCh38:chr15:100:A:G"]
        assert [item.accession for item in assertions] == [
            "VariationID:100",
            "VariationID:900",
        ]
    finally:
        instance.close()


# ---------------------------------------------------------------------- integrity


def test_sha256_mismatch_refuses_to_open_the_release(tmp_path: Path) -> None:
    wrong = "0" * 64
    with pytest.raises(AdapterUnavailableError) as excinfo:
        ClinvarVcfAdapter(FIXTURE, expected_sha256=wrong)
    message = str(excinfo.value)
    assert FIXTURE.name in message
    assert wrong in message


def test_integrity_failure_never_echoes_the_file_s_contents() -> None:
    """PRIV-09: the message names the file and the digests, and nothing else."""
    with pytest.raises(AdapterUnavailableError) as excinfo:
        ClinvarVcfAdapter(FIXTURE, expected_sha256="0" * 64)
    message = str(excinfo.value)
    for leaked in ("CLNSIG", "Pathogenic", "40199604", "BUB1B", "NC_000015"):
        assert leaked not in message


def test_construction_without_an_integrity_pin_fails_closed() -> None:
    """An unpinned 185 MB download is exactly what quietly becomes a different file."""
    with pytest.raises(AdapterUnavailableError, match="without an integrity pin"):
        ClinvarVcfAdapter(FIXTURE)


def test_the_shipped_md5_sidecar_is_an_accepted_pin(tmp_path: Path) -> None:
    """NCBI ships ``clinvar.vcf.gz.md5``; the composition root may wire either source."""
    digest = hashlib.new("md5", FIXTURE.read_bytes(), usedforsecurity=False).hexdigest()
    sidecar = tmp_path / "clinvar.vcf.gz.md5"
    sidecar.write_text(f"{digest}  /netmnt/vast01/vcf_GRCh38/clinvar_20260822.vcf.gz\n", "utf-8")
    assert read_shipped_md5(sidecar) == digest

    instance = ClinvarVcfAdapter(FIXTURE, expected_md5=digest)
    try:
        assert instance.version == EXPECTED_VERSION
    finally:
        instance.close()

    with pytest.raises(AdapterUnavailableError, match="md5 check"):
        ClinvarVcfAdapter(FIXTURE, expected_md5="f" * 32)


def test_malformed_md5_sidecar_is_rejected(tmp_path: Path) -> None:
    sidecar = tmp_path / "clinvar.vcf.gz.md5"
    sidecar.write_text("not-a-digest  clinvar.vcf.gz\n", encoding="utf-8")
    with pytest.raises(AdapterUnavailableError, match="hex digest"):
        read_shipped_md5(sidecar)


def test_missing_release_and_missing_index_are_distinct_failures(tmp_path: Path) -> None:
    with pytest.raises(AdapterUnavailableError, match="never fetches anything"):
        ClinvarVcfAdapter(tmp_path / "absent.vcf.gz", expected_sha256="0" * 64)

    orphan = tmp_path / "orphan.vcf.gz"
    orphan.write_bytes(FIXTURE.read_bytes())
    with pytest.raises(AdapterUnavailableError, match="tabix index"):
        ClinvarVcfAdapter(orphan, expected_sha256=hash_file(orphan))


# ------------------------------------------------------------------ header guards


def test_a_non_clinvar_vcf_is_refused(tmp_path: Path) -> None:
    """Parsing another VCF's INFO with ClinVar's vocabulary invents provenance."""
    text = MINIMAL_HEADER.replace("##source=ClinVar", "##source=SomeOtherTool")
    path = write_indexed_vcf(tmp_path, text + "15\t100\t9\tA\tG\t.\t.\t.\n")
    with pytest.raises(AdapterUnavailableError, match="##source=ClinVar"):
        open_adapter(path)


def test_a_grch37_release_is_refused_for_a_grch38_adapter(tmp_path: Path) -> None:
    """GP-11: the same locus differs by megabases between assemblies."""
    text = MINIMAL_HEADER.replace("##reference=GRCh38", "##reference=GRCh37")
    path = write_indexed_vcf(tmp_path, text + "15\t100\t9\tA\tG\t.\t.\t.\n")
    with pytest.raises(AdapterUnavailableError, match="GRCh37 release"):
        open_adapter(path)


def test_a_release_without_a_filedate_is_refused(tmp_path: Path) -> None:
    """No release string means no citable version, and EvidenceItem refuses one."""
    text = MINIMAL_HEADER.replace("##fileDate=2026-08-22\n", "")
    path = write_indexed_vcf(tmp_path, text + "15\t100\t9\tA\tG\t.\t.\t.\n")
    with pytest.raises(AdapterUnavailableError, match="fileDate"):
        open_adapter(path)


def test_a_grch37_adapter_joins_grch37_ids(tmp_path: Path) -> None:
    """The build in the join key comes from the release, never a hardcoded GRCh38.

    Regression guard for a bug that fails silently in the worst way: a hardcoded
    build makes the query key and the record key disagree, so every lookup
    returns nothing and the run reports "ClinVar has no record" for the entire
    genome while raising no error at all.
    """
    text = MINIMAL_HEADER.replace("##reference=GRCh38", "##reference=GRCh37")
    path = write_indexed_vcf(
        tmp_path,
        text + "15\t100\t9\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel\n",
    )
    instance = ClinvarVcfAdapter(path, expected_sha256=hash_file(path), build=GenomeBuild.GRCH37)
    try:
        assert instance.version == "GRCh37-2026-08-22"
        assert set(instance.assertions(["GRCh37:chr15:100:A:G"])) == {"GRCh37:chr15:100:A:G"}
        with pytest.raises(GenomeBuildMismatchError):
            instance.assertions(["GRCh38:chr15:100:A:G"])
    finally:
        instance.close()


def test_cross_build_variant_ids_are_refused_not_silently_missed(
    adapter: ClinvarVcfAdapter,
) -> None:
    """A GRCh37 key against a GRCh38 release must raise, not return "no record"."""
    with pytest.raises(GenomeBuildMismatchError):
        adapter.assertions(["GRCh37:chr15:40199604:A:C"])


def test_malformed_variant_ids_are_refused_without_echoing_them(
    adapter: ClinvarVcfAdapter,
) -> None:
    with pytest.raises(ValueError, match="colon-separated fields") as excinfo:
        adapter.assertions(["chr15:40199604:A:C"])
    assert "40199604" not in str(excinfo.value)

    with pytest.raises(ValueError, match="not an integer") as excinfo:
        adapter.assertions(["GRCh38:chr15:not-a-position:A:C"])
    assert "not-a-position" not in str(excinfo.value)


def test_build_property_reports_the_release_assembly(adapter: ClinvarVcfAdapter) -> None:
    assert adapter.build is GenomeBuild.GRCH38


# ------------------------------------------------------------------- batch lookup


def test_region_merging_coalesces_only_within_the_window() -> None:
    """Query-count optimisation; correctness still rests on the exact key match."""
    assert merge_query_regions([], 100) == ()
    assert merge_query_regions([5, 5, 5], 100) == ((5, 5),)
    assert merge_query_regions([300, 100, 150], 100) == ((100, 150), (300, 300))
    assert merge_query_regions([300, 100, 150], 200) == ((100, 300),)
    assert merge_query_regions([100, 150, 400], 100) == ((100, 150), (400, 400))
    assert merge_query_regions([100, 150, 400], 0) == ((100, 100), (150, 150), (400, 400))


def test_a_record_spanning_two_query_regions_is_counted_once(tmp_path: Path) -> None:
    """Regression: tabix returns everything that *overlaps* a region.

    A deletion anchored at 100 with a 400-base REF overlaps both the (100, 100)
    and the (400, 400) query spans, so without a start-position guard its
    assertion is appended twice and the variant appears to carry duplicate
    ClinVar calls — inflating the apparent weight of a single submission.
    """
    long_ref = "A" * 400
    path = write_indexed_vcf(
        tmp_path,
        MINIMAL_HEADER + f"15\t100\t9\t{long_ref}\tA\t.\t.\tCLNSIG=Pathogenic;"
        "CLNREVSTAT=reviewed_by_expert_panel\n",
    )
    instance = ClinvarVcfAdapter(path, expected_sha256=hash_file(path), merge_window_bp=0)
    try:
        deletion = f"GRCh38:chr15:100:{long_ref}:A"
        result = instance.assertions([deletion, "GRCh38:chr15:400:A:T"])
        assert len(result[deletion]) == 1
    finally:
        instance.close()


def test_every_record_in_the_slice_is_retrievable_in_one_batch(
    adapter: ClinvarVcfAdapter,
) -> None:
    """The batch path, exercised over the whole fixture at once.

    Also pins the germline/non-germline split: 1425 records in, and exactly the
    ones carrying a CLNSIG come back.
    """
    rows = fixture_records()
    contig_to_ucsc = {"7": "chr7", "15": "chr15", "17": "chr17"}
    ids = [
        f"GRCh38:{contig_to_ucsc[row[0]]}:{row[1]}:{row[3]}:{row[4]}"
        for row in rows
        if "," not in row[4]
    ]
    expected = sum(1 for row in rows if info_of(row).get("CLNSIG"))

    result = adapter.assertions(ids)
    assert len(result) == expected
    assert len(result) < len(ids), "the non-germline records must not have been invented into calls"
    assert all(len(value) == 1 for value in result.values())


def test_a_large_batch_of_misses_costs_no_more_than_a_batch_of_hits(
    adapter: ClinvarVcfAdapter,
) -> None:
    """A 5,000-variant batch must not be 5,000 scans. Correctness proxy for that.

    Timing is not asserted (it would be flaky in CI); what is asserted is that a
    batch far larger than the fixture still resolves exactly the records that
    exist, which is only tractable because lookups are tabix region queries.
    """
    present = {(row[0], row[1], row[3], row[4]) for row in fixture_records()}
    misses = [
        f"GRCh38:chr15:{position}:A:T"
        for position in range(40_199_000, 40_206_000)
        if ("15", str(position), "A", "T") not in present
    ][:5_000]
    assert len(misses) == 5_000
    result = adapter.assertions([*misses, BUB1B_LIKELY_BENIGN])
    assert set(result) == {BUB1B_LIKELY_BENIGN}


# --------------------------------------------------------------- representation
#
# The adapter used to build its join key from the raw VCF columns while ingestion
# built the proband's key from a trimmed, minimal representation. Both were
# internally consistent and together they could not join. The failure has no error
# and no log line: the assertion is simply absent, which this pipeline reads as
# "ClinVar has no record", which promotes the variant as novel. Every test below
# is a spelling of that one failure.


class _RepeatReference:
    """A 1-based inclusive reference over one short contig, for shift tests.

    ``chr15`` 40200230..40200245 is ``TTTTTGCCCCCGATCG``: a five-base C tract at
    40200236-40200240 in which a single-C insertion has six legal spellings.
    """

    START = 40_200_230
    SEQUENCE = "TTTTTGCCCCCGATCG"

    def fetch(self, contig: str, start: int, end: int) -> str:
        if contig != "chr15":
            raise KeyError(contig)
        offset = start - self.START
        if offset < 0 or end - self.START >= len(self.SEQUENCE):
            raise IndexError(start)
        return self.SEQUENCE[offset : end - self.START + 1]


def test_equivalent_nonminimal_alleles_join(tmp_path: Path) -> None:
    """A non-minimal ClinVar record reaches the trimmed proband key.

    ClinVar's ``40200239 AT>AG`` and the normalised proband's ``40200240 T>G`` are
    the same substitution at the same base. Keyed on the raw columns they are two
    different strings and the pathogenic assertion disappears; canonicalised
    through the shared rule they are one key.

    Both directions are asserted. The caller's ID is canonicalised too, so the
    non-minimal spelling of the *query* also resolves — an adapter that only joins
    when its caller happened to normalise first has put its correctness somewhere
    else.
    """
    path = write_indexed_vcf(
        tmp_path,
        MINIMAL_HEADER + "15\t40200239\t7777\tAT\tAG\t.\t.\tCLNSIG=Pathogenic;"
        "CLNREVSTAT=reviewed_by_expert_panel;CLNDN=Mosaic_variegated_aneuploidy\n",
    )
    instance = open_adapter(path)
    try:
        trimmed = "GRCh38:chr15:40200240:T:G"
        result = instance.assertions([trimmed])
        assert trimmed in result, "the trimmed proband key must receive the assertion"
        assertion = result[trimmed][0]
        assert assertion.significance == "Pathogenic"
        assert assertion.star_rating == 3
        assert assertion.accession == "VariationID:7777"

        non_minimal = "GRCh38:chr15:40200239:AT:AG"
        assert instance.assertions([non_minimal])[non_minimal] == result[trimmed]
    finally:
        instance.close()


def test_a_non_minimal_record_is_not_also_reachable_under_its_raw_key(tmp_path: Path) -> None:
    """Canonicalisation moves the key; it does not add a second one.

    Two keys for one record would double-count a single submission wherever the
    caller passed both spellings, which inflates the apparent weight of one
    unreviewed assertion. The record answers to the canonical form only, and the
    caller's own spelling is what the result is keyed by.
    """
    path = write_indexed_vcf(
        tmp_path,
        MINIMAL_HEADER + "15\t40200239\t7777\tAT\tAG\t.\t.\tCLNSIG=Pathogenic\n",
    )
    instance = open_adapter(path)
    try:
        both = ["GRCh38:chr15:40200240:T:G", "GRCh38:chr15:40200239:AT:AG"]
        result = instance.assertions(both)
        # Both spellings answer, each once. Keying the canonical form to only the
        # last caller ID that produced it would answer one and silently drop the
        # other — the same class of silent miss, one level up.
        assert set(result) == set(both)
        assert all(len(value) == 1 for value in result.values())
        assert result[both[0]] == result[both[1]]
    finally:
        instance.close()


def test_query_spans_merge_as_intervals_not_points() -> None:
    """A query is a range once a reference is configured, so merging must be too.

    ``merge_query_regions`` is the point-wise special case and keeps its behaviour;
    the interval form is what lets a repeat-tract query cover every position an
    equivalent source record could be spelled at.
    """
    assert merge_query_spans([], 100) == ()
    assert merge_query_spans([(100, 106)], 0) == ((100, 106),)
    assert merge_query_spans([(100, 106), (104, 110)], 0) == ((100, 110),)
    assert merge_query_spans([(100, 106), (300, 300)], 0) == ((100, 106), (300, 300))
    assert merge_query_spans([(300, 300), (100, 106)], 200) == ((100, 300),)
    # A contained span must not shrink the one that holds it.
    assert merge_query_spans([(100, 400), (150, 160)], 0) == ((100, 400),)


def test_the_pinned_release_is_already_minimally_represented() -> None:
    """Why the fix above is a guard rather than a repair, stated rather than assumed.

    Every one of the 1425 records in the committed slice is already in minimal
    form, and a full scan of the 2026-08-22 release found 0 non-minimal ALT
    entries in 4,467,990 records. So on *this* release the trimming asymmetry
    costs nothing today — the exposure is a future release, a differently-produced
    VCF, or the left-alignment class below, which the release does not protect
    against at all. Recorded here so nobody reads the fix as dead code, and so a
    release that stops being minimal fails a test instead of losing assertions.
    """
    non_minimal = [
        row
        for row in fixture_records()
        for alt in row[4].split(",")
        if len(row[3]) > 1 and len(alt) > 1 and (row[3][0] == alt[0] or row[3][-1] == alt[-1])
    ]
    assert non_minimal == []


def test_a_reference_reconciles_a_shifted_indel_spelling(tmp_path: Path) -> None:
    """The class the release's own minimality does NOT protect against.

    ClinVar stores the left-most spelling of an insertion in a repeat. A caller
    whose VCF spells the identical event at the right-hand end of the same tract
    has a different key, and trimming cannot bridge it — undoing a right shift
    needs the reference bases to the left. Without a reference the lookup misses;
    with one it joins.
    """
    path = write_indexed_vcf(
        tmp_path,
        MINIMAL_HEADER + "15\t40200235\t8888\tG\tGC\t.\t.\tCLNSIG=Likely_pathogenic\n",
    )
    right_shifted = "GRCh38:chr15:40200240:C:CC"

    without = open_adapter(path)
    try:
        assert without.assertions([right_shifted]) == {}
        assert without.representation_status is LeftAlignmentStatus.UNAVAILABLE_NO_REFERENCE
        limitation = without.representation_limitation
        assert limitation is not None
        assert "absence of information" in limitation
    finally:
        without.close()

    with_reference = open_adapter(path, reference=_RepeatReference())
    try:
        joined = with_reference.assertions([right_shifted])
        assert right_shifted in joined, "the shifted spelling must reach the ClinVar record"
        assert joined[right_shifted][0].significance == "Likely pathogenic"
        assert with_reference.representation_status is LeftAlignmentStatus.APPLIED
        assert with_reference.representation_limitation is None
    finally:
        with_reference.close()


def test_a_broken_reference_is_not_reported_as_applied(tmp_path: Path) -> None:
    """The adapter must not claim reference-backed canonicalisation it did not get.

    Reproduced before the fix: `representation_status` was derived from
    `self._reference is not None`, so a FASTA raising on every read reported
    APPLIED over trim-only join keys. A ClinVar miss then reads as "no assertion
    on record", which is absence of information and not evidence of benignity
    (GP-14).
    """

    class _Broken:
        def fetch(self, contig: str, start: int, end: int) -> str:
            raise OSError("handle closed")

    path = write_indexed_vcf(
        tmp_path,
        MINIMAL_HEADER + "15\t40200235\t8888\tG\tGC\t.\t.\tCLNSIG=Likely_pathogenic\n",
    )
    instance = open_adapter(path, reference=_Broken())
    try:
        # Degraded, not dead: an unreadable base must not propagate as an exception.
        assert instance.assertions(["GRCh38:chr15:40200240:C:CC"]) == {}
        assert instance.representation_status is (LeftAlignmentStatus.INCOMPLETE_REFERENCE_UNUSABLE)
        limitation = instance.representation_limitation
        assert limitation is not None
        assert "trimming only" in limitation
    finally:
        instance.close()


def test_a_shifted_release_record_is_still_found_from_the_left_most_query(tmp_path: Path) -> None:
    """The mirror case: the *release* holds the right-shifted spelling.

    An index query finds records by the span they occupy, so a record spelled at
    the right of a repeat sits outside a span queried at the left of it and is
    never fetched at all — canonicalising it afterwards would be too late. The
    query therefore reaches out to the right-most equivalent position, a bound read
    from the reference rather than a guessed padding constant.
    """
    path = write_indexed_vcf(
        tmp_path,
        MINIMAL_HEADER + "15\t40200240\t9999\tC\tCC\t.\t.\tCLNSIG=Pathogenic\n",
    )
    instance = ClinvarVcfAdapter(
        path,
        expected_sha256=hash_file(path),
        merge_window_bp=0,
        reference=_RepeatReference(),
    )
    try:
        left_most = "GRCh38:chr15:40200235:G:GC"
        result = instance.assertions([left_most])
        assert left_most in result
        assert result[left_most][0].significance == "Pathogenic"
    finally:
        instance.close()


def test_canonicalisation_does_not_duplicate_an_assertion_across_query_spans(
    tmp_path: Path,
) -> None:
    """The de-duplication that replaced the raw-position guard.

    Records are no longer dropped by comparing their raw POS to the queried span —
    that guard is exactly what would discard a ``100 AT>AG`` record answering a
    query for 101. Identity is now (join key, ClinVar record), so a single
    submission returned by two overlapping fetches is still counted once while two
    genuinely distinct submissions at one allele are both kept.
    """
    long_ref = "A" * 400
    path = write_indexed_vcf(
        tmp_path,
        MINIMAL_HEADER
        + f"15\t100\t9\t{long_ref}\tA\t.\t.\t"
        + "CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel\n"
        + "15\t500\t10\tG\tA\t.\t.\tCLNSIG=Benign\n"
        + "15\t500\t11\tG\tA\t.\t.\tCLNSIG=Likely_benign\n",
    )
    instance = ClinvarVcfAdapter(path, expected_sha256=hash_file(path), merge_window_bp=0)
    try:
        deletion = f"GRCh38:chr15:100:{long_ref}:A"
        result = instance.assertions([deletion, "GRCh38:chr15:400:A:T", "GRCh38:chr15:500:G:A"])
        assert len(result[deletion]) == 1
        assert len(result["GRCh38:chr15:500:G:A"]) == 2, "two submissions are two assertions"
    finally:
        instance.close()


def test_repeat_lookups_are_byte_identical_with_and_without_a_reference(tmp_path: Path) -> None:
    """GP-30 over the new code path, in both representation states."""
    path = write_indexed_vcf(
        tmp_path,
        MINIMAL_HEADER + "15\t40200235\t8888\tG\tGC\t.\t.\tCLNSIG=Likely_pathogenic\n",
    )
    ids = ["GRCh38:chr15:40200240:C:CC", "GRCh38:chr15:40200235:G:GC"]
    for reference in (None, _RepeatReference()):
        instance = ClinvarVcfAdapter(path, expected_sha256=hash_file(path), reference=reference)
        try:
            assert dump(instance.assertions(ids)) == dump(instance.assertions(ids))
        finally:
            instance.close()


# ------------------------------------------- adversarial review, findings 2 and 3
#
# Finding 2 was raised against the gnomAD adapter, where it is closed: htslib puts
# the queried region — a proband coordinate — into the message of anything it
# raises, and an unwrapped backend failure therefore prints it to the terminal,
# the log, a crash report and an agent's context. The *same* backend is reached
# the same way here, and this adapter's fetch was not wrapped, so the leak was
# live in the clinical slot as well. Verified against the real pysam:
#
#     >>> pysam.TabixFile(release).fetch("nope", 1, 2)
#     ValueError: could not create iterator for region 'nope:2-2'
#
# Finding 3 is the ClinVar half of ADR 0018, measured here the way the gnomAD half
# was measured, against the real 2026-08-22 release and GRCh38.


class _ExplodingTabix:
    """A tabix handle whose fetch fails the way htslib's does: with the region."""

    def __init__(self, inner: object, *, on_iteration: bool) -> None:
        self._inner = inner
        self._on_iteration = on_iteration

    @property
    def contigs(self) -> list[str]:
        return list(getattr(self._inner, "contigs"))  # noqa: B009 - untyped backend

    @property
    def header(self) -> Iterator[str]:
        return iter(getattr(self._inner, "header"))  # noqa: B009 - untyped backend

    def fetch(self, reference: str, start: int, end: int) -> Iterator[str]:
        region = f"{reference}:{start + 1}-{end}"
        if not self._on_iteration:
            msg = f"could not create iterator for region '{region}'"
            raise ValueError(msg)

        def _iterate() -> Iterator[str]:
            msg = f"htslib failed reading {region}"
            raise OSError(msg)
            yield ""  # pragma: no cover - unreachable, makes this a generator

        return _iterate()

    def close(self) -> None:
        closer = getattr(self._inner, "close", None)
        if closer is not None:
            closer()


@pytest.mark.parametrize("on_iteration", [False, True])
def test_a_backend_failure_never_echoes_the_queried_region(
    monkeypatch: pytest.MonkeyPatch, on_iteration: bool
) -> None:
    """PRIV-09 on the error path, both at the call and at the first ``next()``.

    htslib defers work to the iterator, so guarding only the ``fetch`` call would
    move the leak rather than remove it. The replacement is raised with the
    context suppressed: chaining would print the original message again.
    """
    instance = open_adapter(FIXTURE)
    # Replaced after construction so the header parse and the contig map are the
    # real ones: what is under test is the lookup path, not the open.
    monkeypatch.setattr(
        instance, "_tabix", _ExplodingTabix(instance._tabix, on_iteration=on_iteration)
    )
    with pytest.raises(AdapterUnavailableError) as excinfo:
        instance.assertions([UNREVIEWED_PATHOGENIC])
    instance.close()

    rendered = "".join(
        traceback.format_exception(type(excinfo.value), excinfo.value, excinfo.value.__traceback__)
    )
    assert "40200239" not in rendered, "the position reached the traceback"
    assert "could not create iterator for region" not in rendered, "htslib's text was chained"
    assert "htslib failed reading" not in rendered, "htslib's text was chained"
    assert "ClinVar" in rendered
    # The backend's exception *class* survives — a diagnostic carrying no patient
    # data — while its message and its frame do not.
    assert ("ValueError" in rendered) is not on_iteration
    assert ("OSError" in rendered) is on_iteration
    assert re.search(r"<region:[0-9a-f]{8}>", rendered), "no correlation handle was given"


def test_a_lookup_after_close_is_refused_rather_than_answered_as_absence() -> None:
    """A closed handle must not answer 'ClinVar has no record' for a whole batch."""
    instance = open_adapter(FIXTURE)
    assert instance.assertions([UNREVIEWED_PATHOGENIC])
    instance.close()
    with pytest.raises(AdapterUnavailableError) as excinfo:
        instance.assertions([UNREVIEWED_PATHOGENIC])
    message = str(excinfo.value)
    assert "closed" in message
    assert "40200239" not in message
    # Idempotent, as before.
    instance.close()


# --------------------------------------- finding 3: the ClinVar half of ADR 0018
#
# ADR 0018 measured the gnomAD half against the real v4.1 exomes chr21 shard. This
# is the same measurement on the same terms against the real ClinVar release, and
# it is a test rather than a paragraph so that a release which stops behaving this
# way fails instead of quietly losing assertions.
#
# The window is chr17:43,000,000-43,520,000 — 520 kb, the same width as the gnomAD
# measurement, over the BRCA1 neighbourhood, which is the densest curated exonic
# region ClinVar has. Public reference data throughout.
#
# Measured 2026-08-28 against ClinVar 2026-08-22 and GRCh38_no_alt:
#
#   ClinVar records in the window          15,862
#   indel ALT alleles                       3,595
#     in a repeat tract                     2,215  (61.6%)
#     of those, no germline CLNSIG              4
#   distinct right-shifted spellings        2,215
#   join WITHOUT a reference                    0
#   join WITH a reference                   2,211  (= 2,215 - the 4 with no CLNSIG)
#   of the recovered, Pathogenic/LP         1,761
#
# The accounting closes exactly: every right-shifted spelling that ClinVar holds a
# germline classification for is recovered by the reference, and the only ones
# that are not are the four records that carry no germline classification at all,
# which GP-14 requires to stay omitted.

FULL_CLINVAR = Path(
    os.environ.get(
        "MVA_CLINVAR_VCF",
        str(REPO_ROOT.parent / "mva-resources" / "clinvar" / "clinvar.vcf.gz"),
    )
)
FULL_REFERENCE = Path(
    os.environ.get(
        "MVA_REFERENCE_FASTA",
        str(REPO_ROOT.parent / "mva-resources" / "reference" / "GRCh38_no_alt.fa"),
    )
)

requires_full_clinvar = pytest.mark.skipif(
    not (FULL_CLINVAR.is_file() and FULL_REFERENCE.is_file()),
    reason=f"full ClinVar release or GRCh38 FASTA absent ({FULL_CLINVAR}, {FULL_REFERENCE})",
)

MEASURED_CONTIG = "chr17"
MEASURED_START = 43_000_000
MEASURED_END = 43_520_000
MEASURED = {
    "records": 15_862,
    "indel_alleles": 3_595,
    "in_repeat": 2_215,
    "no_germline_classification": 4,
    "joins_without_reference": 0,
    "joins_with_reference": 2_211,
    "recovered_pathogenic": 1_761,
}


def roll_right_max(
    contig: str, position: int, ref: str, alt: str, reference: ReferenceLookup
) -> tuple[int, str, str, int]:
    """The right-most VCF-conventional spelling of one anchored indel.

    The mirror of :func:`mva.alleles._left_shift`, written out here rather than
    imported because it is the *wrong* representation on purpose: it manufactures
    the spelling a VCF that was never left-aligned would carry, which is the input
    the adapter has to reconcile. ``steps == 0`` means the event has exactly one
    legal spelling and so is not in a repeat tract.

    ``rightmost_equivalent_position`` is deliberately not used for the repeat-tract
    count: it over-reaches by one base by design (its docstring says so), because
    it bounds a *fetch window* and under-reaching there would re-open the silent
    miss. Counting repeats with it would call every indel a repeat.
    """
    steps = 0
    for _ in range(1000):
        if len(ref) > 1 and len(alt) > 1:
            break  # a complex substitution: not a shiftable indel
        if ref[0] != alt[0]:
            break  # a delins, not an anchored pure indel: it has no other spelling
        insertion = len(alt) > len(ref)
        sequence = alt[1:] if insertion else ref[1:]
        probe = position + 1 if insertion else position + len(ref)
        try:
            following = reference.fetch(contig, probe, probe).upper()
        except (KeyError, ValueError):
            break
        if following != sequence[0]:
            break  # the tract ends here; this is the right-most spelling
        rotated = sequence[1:] + sequence[0]
        position += 1
        anchor = sequence[0]
        ref, alt = (anchor, anchor + rotated) if insertion else (anchor + rotated, anchor)
        steps += 1
    return position, ref, alt, steps


@dataclass(frozen=True)
class _WindowScan:
    """What the measured window holds, before any adapter is asked anything."""

    records: int
    indel_alleles: int
    in_repeat: int
    no_germline: int
    right_shifted_ids: tuple[str, ...]


def _scan_measured_window(reference: ReferenceLookup) -> _WindowScan:
    """Read the window and manufacture one right-shifted spelling per repeat indel."""
    records = 0
    indel_alleles = 0
    in_repeat = 0
    no_germline = 0
    queries: list[str] = []
    seen: set[str] = set()

    handle = pysam.TabixFile(str(FULL_CLINVAR))
    try:
        for line in handle.fetch(
            MEASURED_CONTIG.removeprefix("chr"), MEASURED_START - 1, MEASURED_END
        ):
            columns = line.split("\t", 8)
            records += 1
            position = int(columns[1])
            ref = columns[3].strip().upper()
            classified = "CLNSIG=" in columns[7]
            for raw_alt in columns[4].split(","):
                alt = raw_alt.strip().upper()
                if not is_sequence_allele(ref) or not is_sequence_allele(alt):
                    continue
                trimmed = trim_parsimoniously(position, ref, alt)
                if len(trimmed[1]) == len(trimmed[2]):
                    continue  # a substitution, not an indel
                indel_alleles += 1
                canonical = canonicalise_allele(
                    contig=MEASURED_CONTIG,
                    position=position,
                    ref=ref,
                    alt=alt,
                    reference=reference,
                )
                shifted = roll_right_max(
                    MEASURED_CONTIG,
                    canonical.position,
                    canonical.ref,
                    canonical.alt,
                    reference,
                )
                if shifted[3] == 0:
                    continue  # one legal spelling only: no representation risk
                in_repeat += 1
                if not classified:
                    no_germline += 1
                variant_id = f"GRCh38:{MEASURED_CONTIG}:{shifted[0]}:{shifted[1]}:{shifted[2]}"
                if variant_id not in seen:
                    seen.add(variant_id)
                    queries.append(variant_id)
    finally:
        handle.close()
    return _WindowScan(
        records=records,
        indel_alleles=indel_alleles,
        in_repeat=in_repeat,
        no_germline=no_germline,
        right_shifted_ids=tuple(queries),
    )


@pytest.mark.slow
@requires_full_clinvar
def test_the_measured_cost_of_an_unreferenced_clinvar_join() -> None:
    """ADR 0018's gnomAD measurement, repeated on the clinical slot.

    Every number in ``MEASURED`` is asserted, so the claim in the report is the
    claim the code produces. A ClinVar release that stops being reconcilable this
    way fails here rather than quietly returning fewer pathogenic assertions.
    """
    reference = open_reference_fasta(FULL_REFERENCE)
    scan = _scan_measured_window(reference)
    queries = list(scan.right_shifted_ids)

    md5 = read_shipped_md5(FULL_CLINVAR.with_name(FULL_CLINVAR.name + ".md5"))
    degraded = ClinvarVcfAdapter(FULL_CLINVAR, expected_md5=md5)
    try:
        without_reference = degraded.assertions(queries)
        assert degraded.representation_status is LeftAlignmentStatus.UNAVAILABLE_NO_REFERENCE
        assert degraded.representation_limitation is not None
    finally:
        degraded.close()

    referenced = ClinvarVcfAdapter(FULL_CLINVAR, expected_md5=md5, reference=reference)
    try:
        with_reference = referenced.assertions(queries)
        assert referenced.representation_status is LeftAlignmentStatus.APPLIED
        assert referenced.representation_limitation is None
    finally:
        referenced.close()

    recovered = set(with_reference) - set(without_reference)
    strongly_pathogenic = {
        variant_id
        for variant_id in recovered
        if any(
            assertion.significance.startswith(("Pathogenic", "Likely_pathogenic"))
            and "Conflicting" not in assertion.significance
            for assertion in with_reference[variant_id]
        )
    }

    assert scan.records == MEASURED["records"]
    assert scan.indel_alleles == MEASURED["indel_alleles"]
    assert scan.in_repeat == MEASURED["in_repeat"]
    assert scan.no_germline == MEASURED["no_germline_classification"]
    assert len(queries) == MEASURED["in_repeat"], "every repeat indel has a distinct spelling here"
    assert len(without_reference) == MEASURED["joins_without_reference"]
    assert len(with_reference) == MEASURED["joins_with_reference"]
    # A reference only ever adds joins; it never moves or removes one.
    assert not set(without_reference) - set(with_reference)
    # And the gap closes exactly: the only right-shifted spellings that stay
    # unjoined are the records carrying no germline classification at all, which
    # GP-14 requires to stay omitted rather than be reported as benign.
    assert len(queries) - len(with_reference) == MEASURED["no_germline_classification"]
    assert len(strongly_pathogenic) == MEASURED["recovered_pathogenic"]


# ------------------------------------- the same short-index door, on the clinical slot
#
# Found while closing finding 1 on the gnomAD adapter, and reproduced here before
# being fixed: the identical failure exists on this adapter and costs more. The
# release's *data* is sha256-pinned; its *index* is only checked to exist. An
# index built from a shorter file region-queries the complete release perfectly
# happily and answers "no record" for everything past its reach — which is not
# evidence of benignity (GP-14), and which the measurement above prices at 1,761
# Pathogenic/Likely_pathogenic assertions in a single 520 kb window.


def stale_indexed_release(tmp_path: Path) -> tuple[Path, str]:
    """A COMPLETE ClinVar slice beside an index built from only its first half."""
    with gzip.open(FIXTURE, "rt", encoding="utf-8") as handle:
        lines = handle.readlines()
    header = [line for line in lines if line.startswith("#")]
    records = [line for line in lines if not line.startswith("#")]

    half = tmp_path / "half.vcf"
    half.write_text("".join(header + records[: len(records) // 2]), encoding="utf-8")
    pysam.tabix_compress(str(half), str(tmp_path / "half.vcf.gz"), force=True)
    pysam.tabix_index(str(tmp_path / "half.vcf.gz"), preset="vcf", force=True)

    release = tmp_path / "clinvar.vcf.gz"
    release.write_bytes(FIXTURE.read_bytes())
    release.with_name(release.name + ".tbi").write_bytes(
        (tmp_path / "half.vcf.gz.tbi").read_bytes()
    )
    columns = records[-1].split("\t")
    return release, f"GRCh38:chr{columns[0]}:{columns[1]}:{columns[3]}:{columns[4]}"


def test_an_index_that_does_not_reach_the_end_of_the_release_is_refused(
    tmp_path: Path,
) -> None:
    """Reproduced by execution: the assertion was silently absent, not reported missing.

    Before this check the adapter constructed cleanly, the sha256 pin passed
    (the data really is intact), and the query returned ``{}`` — indistinguishable
    from "ClinVar has nothing on record here".
    """
    release, past_the_index = stale_indexed_release(tmp_path)
    with pytest.raises(AdapterUnavailableError) as excinfo:
        open_adapter(release)
    message = str(excinfo.value)
    assert "index" in message
    assert "does not reach" in message
    # PRIV-09: the file and the shape of the failure, never a coordinate.
    assert past_the_index.split(":")[2] not in message


def test_a_healthy_release_and_the_real_one_are_not_rejected(tmp_path: Path) -> None:
    """The control, on both the committed fixture and every in-test release.

    A check that rejected healthy pairs would be worse than no check: it would push
    the next operator to delete it.
    """
    instance = open_adapter(FIXTURE)
    try:
        assert instance.assertions([UNREVIEWED_PATHOGENIC])
    finally:
        instance.close()

    built = write_indexed_vcf(tmp_path, MINIMAL_HEADER + "15\t100\t1\tA\tG\t.\t.\tCLNSIG=Benign\n")
    instance = open_adapter(built)
    try:
        assert instance.assertions(["GRCh38:chr15:100:A:G"])
    finally:
        instance.close()


@pytest.mark.slow
@requires_full_clinvar
def test_the_real_release_index_reaches_the_end_of_the_release() -> None:
    """The 193 MB release on disk, measured rather than assumed."""
    md5 = read_shipped_md5(FULL_CLINVAR.with_name(FULL_CLINVAR.name + ".md5"))
    instance = ClinvarVcfAdapter(FULL_CLINVAR, expected_md5=md5)
    try:
        assert instance.assertions(["GRCh38:chr15:40200239:A:G"])
    finally:
        instance.close()
