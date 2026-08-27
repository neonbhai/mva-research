"""Unit tests for the annotation adapter layer.

These are mostly *integrity* tests rather than behaviour tests: the interesting
failures in an annotation stage are not "wrong number returned", they are
"absence quietly became zero", "two transcripts quietly became one" and "a
synthetic table quietly stopped saying it was synthetic". Each test below locks
one of those doors.
"""

from __future__ import annotations

import shutil
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from pathlib import Path

import pytest
from pydantic import ValidationError

from mva.annotation import (
    SYNTHETIC_STANDIN_LIMITATION,
    AdapterSet,
    AnnotationResult,
    LocalConsequenceAdapter,
    NullClinicalAdapter,
    annotate_variants,
    compute_manifest,
    load_default_adapters,
    load_manifest,
    resolve_table_path,
)
from mva.clock import Clock, FixedClock
from mva.determinism import stable_hash
from mva.errors import AdapterUnavailableError
from mva.models.base import AssertionTier
from mva.models.evidence import (
    Citation,
    EvidenceCategory,
    EvidenceDirection,
    EvidenceItem,
    EvidenceType,
)
from mva.models.genome import GenomeBuild, GenomicCoordinate
from mva.models.variant import (
    ClinicalAssertion,
    FilterStatus,
    Genotype,
    ImpactSeverity,
    VariantRecord,
    Zygosity,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
KNOWLEDGE_ROOT = REPO_ROOT / "knowledge"
MANIFEST_PATH = KNOWLEDGE_ROOT / "manifests" / "knowledge.yaml"

Coordinate = tuple[str, int, str, str]

#: Two transcript rows in consequences.tsv; a stop_gained on both.
MULTI_TRANSCRIPT: Coordinate = ("chr15", 40_200_000, "C", "T")

#: Present in consequences.tsv, deliberately ABSENT from frequencies.tsv. This is
#: the GP-14 fixture: no frequency data, which is not the same as AF = 0.
NO_FREQUENCY_DATA: Coordinate = ("chr11", 5_000_000, "A", "GT")

#: Two population rows (global + nfe) — checks that populations are not merged.
TWO_POPULATIONS: Coordinate = ("chr15", 40_210_500, "G", "A")

#: An intronic row whose hgvs_p / exon / protein_position cells are blank.
SPARSE_CELLS: Coordinate = ("chr3", 10_000_000, "G", "A")


# --------------------------------------------------------------------------- fixtures


@pytest.fixture
def adapters() -> AdapterSet:
    return load_default_adapters(KNOWLEDGE_ROOT, MANIFEST_PATH)


@pytest.fixture
def clock() -> Clock:
    return FixedClock(datetime(2026, 1, 1, tzinfo=UTC))


def make_record(coordinate: Coordinate, *, line_index: int = 0) -> VariantRecord:
    contig, position, ref, alt = coordinate
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
        source_artifact="tests/fixtures/synthetic/synthetic_case.vcf",
        source_line_index=line_index,
    )


def run(adapters: AdapterSet, clock: Clock, *coordinates: Coordinate) -> AnnotationResult:
    records = [make_record(c, line_index=i) for i, c in enumerate(coordinates)]
    return annotate_variants(records, adapters=adapters, clock=clock)


def evidence_for(result: AnnotationResult, variant_id: str) -> tuple[EvidenceItem, ...]:
    return tuple(item for item in result.evidence if item.subject_id == variant_id)


# --------------------------------------------------------------------- GP-14: absence


@pytest.mark.unit
def test_missing_frequency_is_recorded_as_absent_never_as_zero(
    adapters: AdapterSet, clock: Clock
) -> None:
    """GP-14. The single most dangerous default in this stage is AF = 0.

    A variant the reference cohort has never seen must come out with *no*
    frequency data and an explicit, citable statement that the data is missing.
    Defaulting to zero would turn a coverage gap into the strongest possible
    rarity signal — and would do it most often for the ancestries reference
    cohorts represent worst.
    """
    result = run(adapters, clock, NO_FREQUENCY_DATA)
    (variant,) = result.variants

    assert variant.population_frequencies == ()
    assert variant.has_frequency_data is False
    assert variant.max_allele_frequency() is None

    population_items = [
        item
        for item in evidence_for(result, variant.variant_id)
        if item.category is EvidenceCategory.POPULATION
    ]
    assert len(population_items) == 1, "the gap must be recorded exactly once, not zero times"
    no_data = population_items[0]

    assert no_data.direction is EvidenceDirection.NEUTRAL
    claim = no_data.claim.lower()
    assert "unavailable" in claim
    assert "not evidence of rarity" in claim
    # The claim must not smuggle a number in: there is no number to report.
    assert no_data.numeric_value is None
    assert no_data.payload["frequency_data_available"] is False


@pytest.mark.unit
def test_missing_frequency_is_reported_in_coverage_and_warnings(
    adapters: AdapterSet, clock: Clock
) -> None:
    """The gap is visible at run level too, not only per variant."""
    result = run(adapters, clock, MULTI_TRANSCRIPT, NO_FREQUENCY_DATA)

    assert result.coverage[adapters.frequency.name] == pytest.approx(0.5)
    assert result.coverage[adapters.consequence.name] == pytest.approx(1.0)
    assert any(
        "no population-frequency data" in warning and "GP-14" in warning
        for warning in result.warnings
    )


@pytest.mark.unit
def test_observed_zero_frequency_is_distinct_from_missing(
    adapters: AdapterSet, clock: Clock
) -> None:
    """AF = 0.0 over 152,312 alleles is data; an absent row is not.

    Both variants would look identical under a naive "default to zero" scheme.
    They must not look identical here.
    """
    result = run(adapters, clock, MULTI_TRANSCRIPT, NO_FREQUENCY_DATA)
    observed, missing = result.variants

    assert observed.has_frequency_data is True
    observed_af = observed.max_allele_frequency()
    assert observed_af is not None
    assert observed_af.allele_frequency == 0.0
    assert observed_af.allele_number == 152_312

    assert missing.has_frequency_data is False


# ------------------------------------------------------------- transcript preservation


@pytest.mark.unit
def test_multi_transcript_variant_keeps_every_transcript(
    adapters: AdapterSet, clock: Clock
) -> None:
    """GP-20 / ASSUMPTION-TRANSCRIPT-01: never collapse to the canonical transcript.

    A variant can be benign on MANE-Select and splice-disrupting on the
    tissue-relevant isoform, so both rows survive to the record and both produce
    their own evidence item.
    """
    result = run(adapters, clock, MULTI_TRANSCRIPT)
    (variant,) = result.variants

    transcripts = [c.transcript_id for c in variant.consequences]
    assert transcripts == ["SYNTHT0001.1", "SYNTHT0002.1"]
    assert [c.is_canonical for c in variant.consequences] == [True, False]

    consequence_items = [
        item
        for item in evidence_for(result, variant.variant_id)
        if item.tool == adapters.consequence.name
    ]
    assert len(consequence_items) == 2
    assert {str(item.payload["transcript_id"]) for item in consequence_items} == set(transcripts)


@pytest.mark.unit
def test_worst_impact_for_gene_is_high_for_the_stop_gained_variant(
    adapters: AdapterSet, clock: Clock
) -> None:
    result = run(adapters, clock, MULTI_TRANSCRIPT)
    (variant,) = result.variants

    assert variant.worst_impact_for_gene("SYNTHKIN1") is ImpactSeverity.HIGH
    assert variant.gene_symbols == ("SYNTHKIN1",)


# ------------------------------------------------------------------- GP-18: provenance


@pytest.mark.unit
def test_population_frequency_carries_source_version_and_population(
    adapters: AdapterSet, clock: Clock
) -> None:
    """GP-18. A bare float is unusable: 0.001 is common in Finns and vanishing globally."""
    result = run(adapters, clock, TWO_POPULATIONS)
    (variant,) = result.variants

    assert len(variant.population_frequencies) == 2
    for frequency in variant.population_frequencies:
        assert frequency.source == "SYNTH_gnomAD"
        assert frequency.version == "v0.0-synthetic"
        assert frequency.population
    populations = [f.population for f in variant.population_frequencies]
    assert populations == ["global", "nfe"], "populations are kept separate, not averaged"
    assert [f.provenance_key for f in variant.population_frequencies] == [
        "SYNTH_gnomAD/v0.0-synthetic/global",
        "SYNTH_gnomAD/v0.0-synthetic/nfe",
    ]

    # max_allele_frequency is the conservative rarity signal: the highest across
    # populations, not the global value.
    maximum = variant.max_allele_frequency()
    assert maximum is not None
    assert maximum.population == "nfe"


@pytest.mark.unit
def test_empty_cells_parse_to_none_not_to_zero(adapters: AdapterSet, clock: Clock) -> None:
    """An unscored REVEL is absent, not 0.0 — which would read as 'scored, benign'."""
    result = run(adapters, clock, SPARSE_CELLS)
    (variant,) = result.variants
    (annotation,) = variant.consequences

    assert annotation.hgvs_p is None
    assert annotation.exon is None
    assert annotation.protein_position is None
    assert annotation.amino_acids is None
    assert "REVEL" not in annotation.pathogenicity_scores
    assert annotation.pathogenicity_scores["CADD_phred"] == pytest.approx(2.1)


# ------------------------------------------------------------------ GP-12: tiering


@pytest.mark.unit
def test_consequence_evidence_is_filed_as_computational_prediction(
    adapters: AdapterSet, clock: Clock
) -> None:
    """GP-12: a tool's opinion about a transcript is a prediction, not an observation."""
    result = run(adapters, clock, MULTI_TRANSCRIPT, TWO_POPULATIONS, SPARSE_CELLS)
    consequence_items = [item for item in result.evidence if item.tool == adapters.consequence.name]

    assert consequence_items
    for item in consequence_items:
        assert item.tier is AssertionTier.COMPUTATIONAL_PREDICTION
        assert item.evidence_type is EvidenceType.IN_SILICO_PREDICTION
        assert item.category is EvidenceCategory.CONSEQUENCE
        # The model itself refuses the laundering; assert we never even try.
        assert item.tier is not AssertionTier.OBSERVED_DATA


@pytest.mark.unit
def test_frequency_evidence_carries_a_versioned_citation(
    adapters: AdapterSet, clock: Clock
) -> None:
    """A DATABASE_ASSERTION without a release version is not reproducible."""
    result = run(adapters, clock, TWO_POPULATIONS)
    database_items = [
        item for item in result.evidence if item.tier is AssertionTier.DATABASE_ASSERTION
    ]

    assert database_items
    for item in database_items:
        assert item.evidence_type is EvidenceType.CURATED_DATABASE
        assert item.citation is not None
        assert item.citation.version == "v0.0-synthetic"
        assert item.citation.source == "SYNTH_gnomAD"


@pytest.mark.unit
def test_database_assertion_without_versioned_citation_is_rejected(
    adapters: AdapterSet, clock: Clock
) -> None:
    """The guarantee above is enforced by the type, not by our good manners."""
    result = run(adapters, clock, TWO_POPULATIONS)
    item = next(i for i in result.evidence if i.tier is AssertionTier.DATABASE_ASSERTION)

    without_citation = item.model_dump()
    without_citation["citation"] = None
    with pytest.raises(ValidationError, match="versioned citation"):
        EvidenceItem.model_validate(without_citation)

    unversioned = item.model_dump()
    unversioned["citation"] = Citation(source="SYNTH_gnomAD", identifier="x").model_dump()
    with pytest.raises(ValidationError, match="versioned citation"):
        EvidenceItem.model_validate(unversioned)


# ------------------------------------------------------------ GP-17/GP-20: honest mocks


@pytest.mark.unit
def test_every_evidence_item_declares_its_synthetic_origin(
    adapters: AdapterSet, clock: Clock
) -> None:
    """GP-17 + GP-20: limitations are mandatory, and a mock says it is a mock."""
    result = run(adapters, clock, MULTI_TRANSCRIPT, NO_FREQUENCY_DATA, TWO_POPULATIONS)

    assert result.evidence
    for item in result.evidence:
        assert SYNTHETIC_STANDIN_LIMITATION in item.limitations
        assert "NOT biologically valid" in item.limitations
        assert item.tool and item.tool_version

    assert any("SYNTHETIC stand-in" in warning for warning in result.warnings)


@pytest.mark.unit
def test_adapter_identities_advertise_that_they_are_synthetic(adapters: AdapterSet) -> None:
    """The names alone must be unmistakable in a report footer."""
    assert adapters.consequence.name == "local-tsv-consequence"
    assert adapters.consequence.version == "synthetic-v0.0"
    assert adapters.frequency.name == "local-tsv-frequency"
    assert adapters.frequency.version == "synthetic-v0.0"
    assert all(descriptor.synthetic for descriptor in adapters.descriptors())


@pytest.mark.unit
def test_null_clinical_adapter_returns_nothing_and_says_so(
    adapters: AdapterSet, clock: Clock
) -> None:
    """No ClinVar substitute is shipped; inventing one would be the worst option."""
    assert isinstance(adapters.clinical, NullClinicalAdapter)
    assert adapters.clinical.assertions(["GRCh38:chr15:40200000:C:T"]) == {}

    result = run(adapters, clock, MULTI_TRANSCRIPT)
    assert result.coverage["null-clinical"] == pytest.approx(0.0)
    assert any("not evidence of benignity" in warning for warning in result.warnings)


class _StubClinicalAdapter:
    """A clinical source with one assertion, for exercising the populated path.

    Deliberately does NOT declare a ``synthetic`` property: `is_synthetic` fails
    closed, so this stub still gets the mock disclosure.
    """

    @property
    def name(self) -> str:
        return "stub-clinvar"

    @property
    def version(self) -> str:
        return "test-v0.0"

    def assertions(self, variant_ids: Sequence[str]) -> Mapping[str, tuple[ClinicalAssertion, ...]]:
        return {
            variant_ids[0]: (
                ClinicalAssertion(
                    source="ClinVar",
                    version="2026-01-01",
                    accession="VCV000000001",
                    significance="Pathogenic",
                    review_status="criteria provided, multiple submitters",
                    star_rating=2,
                ),
            )
        }


@pytest.mark.unit
def test_clinical_assertions_are_attached_with_a_versioned_citation(clock: Clock) -> None:
    base = load_default_adapters(KNOWLEDGE_ROOT, MANIFEST_PATH)
    adapters = AdapterSet(
        consequence=base.consequence, frequency=base.frequency, clinical=_StubClinicalAdapter()
    )

    result = run(adapters, clock, MULTI_TRANSCRIPT)
    (variant,) = result.variants
    (assertion,) = variant.clinical_assertions
    assert assertion.significance == "Pathogenic"

    item = next(i for i in result.evidence if i.tool == "stub-clinvar")
    assert item.citation is not None
    assert item.citation.version == "2026-01-01"
    assert item.direction is EvidenceDirection.SUPPORTS
    # Undeclared maturity is treated as a mock (fail closed).
    assert SYNTHETIC_STANDIN_LIMITATION in item.limitations


# ------------------------------------------------------------- manifest integrity


@pytest.mark.unit
def test_committed_manifest_hashes_match_the_files_on_disk() -> None:
    """A hand-edited hash turns an integrity check into decoration."""
    committed = load_manifest(MANIFEST_PATH)
    computed = compute_manifest(KNOWLEDGE_ROOT)
    computed_tables: dict[str, dict[str, object]] = computed["tables"]

    assert set(committed.tables) == set(computed_tables)
    for name, table in sorted(committed.tables.items()):
        assert table.sha256 == computed_tables[name]["sha256"], f"{name} hash is stale"
        assert table.synthetic is True
        assert table.version == "synthetic-v0.0"
        assert resolve_table_path(KNOWLEDGE_ROOT, table).is_file()


@pytest.mark.unit
def test_manifest_sha256_mismatch_refuses_to_load_adapters(tmp_path: Path) -> None:
    """Changed bytes with an unchanged manifest is a stop, not a warning."""
    knowledge_root = tmp_path / "knowledge"
    shutil.copytree(KNOWLEDGE_ROOT, knowledge_root)
    tampered = knowledge_root / "public" / "consequences.tsv"
    added_row = "GRCh38:chr9:99999999:A:T\tTAMPERED_GENE"
    with tampered.open("a", encoding="utf-8") as handle:
        handle.write(f"{added_row}\n")

    with pytest.raises(AdapterUnavailableError) as excinfo:
        load_default_adapters(knowledge_root, knowledge_root / "manifests" / "knowledge.yaml")

    message = str(excinfo.value)
    assert "consequences" in message
    assert "public/consequences.tsv" in message
    # PRIV-09: name the file, never echo what is in it.
    assert "TAMPERED_GENE" not in message


@pytest.mark.unit
def test_missing_table_file_is_named_not_guessed(tmp_path: Path) -> None:
    knowledge_root = tmp_path / "knowledge"
    shutil.copytree(KNOWLEDGE_ROOT, knowledge_root)
    (knowledge_root / "public" / "frequencies.tsv").unlink()

    with pytest.raises(AdapterUnavailableError, match="frequencies"):
        load_default_adapters(knowledge_root, knowledge_root / "manifests" / "knowledge.yaml")


@pytest.mark.unit
def test_missing_manifest_is_an_adapter_error(tmp_path: Path) -> None:
    with pytest.raises(AdapterUnavailableError, match="not found"):
        load_default_adapters(KNOWLEDGE_ROOT, tmp_path / "absent.yaml")


@pytest.mark.unit
def test_malformed_table_fails_loudly_at_construction(tmp_path: Path) -> None:
    """A short row must not shift every downstream column by one, silently."""
    table = tmp_path / "broken.tsv"
    table.write_text(
        "# broken fixture\nvariant_id\tgene_symbol\nGRCh38:chr1:1:A:T\n", encoding="utf-8"
    )
    with pytest.raises(AdapterUnavailableError, match="missing required column"):
        LocalConsequenceAdapter(table, version="synthetic-v0.0")


# ------------------------------------------------------------------- GP-30: determinism


@pytest.mark.unit
def test_repeat_runs_are_identical(clock: Clock) -> None:
    """GP-30. Two independently-wired runs over the same inputs agree byte for byte."""
    coordinates = (MULTI_TRANSCRIPT, NO_FREQUENCY_DATA, TWO_POPULATIONS, SPARSE_CELLS)

    first = run(load_default_adapters(KNOWLEDGE_ROOT, MANIFEST_PATH), clock, *coordinates)
    second = run(load_default_adapters(KNOWLEDGE_ROOT, MANIFEST_PATH), clock, *coordinates)

    assert first == second
    assert stable_hash([v.model_dump(mode="json") for v in first.variants]) == stable_hash(
        [v.model_dump(mode="json") for v in second.variants]
    )
    assert stable_hash([e.model_dump(mode="json") for e in first.evidence]) == stable_hash(
        [e.model_dump(mode="json") for e in second.evidence]
    )
    assert [e.evidence_id for e in first.evidence] == [e.evidence_id for e in second.evidence]
    assert first.warnings == second.warnings


@pytest.mark.unit
def test_input_order_and_membership_are_preserved(adapters: AdapterSet, clock: Clock) -> None:
    """GP-13: annotation annotates. It never drops a record for being uninformative."""
    coordinates = (SPARSE_CELLS, NO_FREQUENCY_DATA, MULTI_TRANSCRIPT)
    result = run(adapters, clock, *coordinates)

    assert len(result.variants) == len(coordinates)
    assert [v.coordinate.contig for v in result.variants] == [c[0] for c in coordinates]


@pytest.mark.unit
def test_unknown_variant_is_annotated_with_nothing_and_still_returned(
    adapters: AdapterSet, clock: Clock
) -> None:
    """Absent from both tables: no consequences, no frequencies, still present."""
    result = run(adapters, clock, ("chr1", 12_345, "A", "G"))
    (variant,) = result.variants

    assert variant.consequences == ()
    assert variant.population_frequencies == ()
    assert variant.has_frequency_data is False
    assert any("no consequence annotation" in warning for warning in result.warnings)
