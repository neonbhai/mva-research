"""Behavioural tests for the evidence store, the ledger and the GP-10 gate.

These are the guarantees the rest of the pipeline is allowed to rely on:

* the schema initialises to an empty, fully-enumerated database;
* every writer is idempotent on its primary key;
* an EvidenceItem survives a write/read round trip with every enum, tuple and
  payload value intact;
* contradictions are stored and retrievable (GP-19);
* a Parquet export is byte-identical on repeat (GP-30), checked with sha256;
* an unsourced or dangling citation refuses to render (GP-10).
"""

from __future__ import annotations

import hashlib
from collections.abc import Iterator
from datetime import UTC, datetime
from pathlib import Path

import pytest

from mva.errors import UnsourcedAssertionError
from mva.evidence import AssertionResolver, EvidenceLedger, EvidenceStore, GraphEdge
from mva.evidence.store import TABLES
from mva.models import (
    ApprovalStatus,
    AssertionTier,
    CandidatePair,
    Citation,
    ComponentScores,
    ConsequenceAnnotation,
    DrugHypothesis,
    EffectDirection,
    EvidenceCategory,
    EvidenceDirection,
    EvidenceItem,
    EvidenceStrength,
    EvidenceType,
    FilterStatus,
    GenomeBuild,
    GenomicCoordinate,
    Genotype,
    ImpactSeverity,
    InheritanceModel,
    InterventionClass,
    ObservationStatus,
    PediatricEvidence,
    PharmacokineticProfile,
    PhaseEvidence,
    PhaseStatus,
    PhenotypeObservation,
    PhenotypeProfile,
    PopulationFrequency,
    RejectionReason,
    VariantRecord,
    Zygosity,
    make_evidence_id,
    make_pair_id,
)

#: Fixed instant. A wall-clock read here would make the determinism test a
#: coin flip (GP-30).
INSTANT = datetime(2026, 1, 1, 12, 30, 0, tzinfo=UTC)
RUN_ID = "run-test-0001"


# --------------------------------------------------------------------- helpers


def make_variant(
    *,
    contig: str = "chr15",
    position: int = 100_000,
    ref: str = "A",
    alt: str = "T",
    zygosity: Zygosity = Zygosity.HET,
    gene: str = "SYNTHA",
) -> VariantRecord:
    """A minimal but fully valid annotated variant."""
    return VariantRecord(
        coordinate=GenomicCoordinate(
            build=GenomeBuild.GRCH38, contig=contig, position=position, ref=ref, alt=alt
        ),
        genotype=Genotype(
            zygosity=zygosity,
            genotype_string="0/1",
            depth=40,
            ref_reads=21,
            alt_reads=19,
            genotype_quality=99,
        ),
        filter_status=FilterStatus.PASS,
        raw_filters=("PASS",),
        quality=880.0,
        consequences=(
            ConsequenceAnnotation(
                gene_symbol=gene,
                gene_id="ENSG00000000001",
                transcript_id="ENST00000000001",
                is_canonical=True,
                is_mane_select=True,
                consequence_terms=("stop_gained", "coding_sequence_variant"),
                impact=ImpactSeverity.HIGH,
                hgvs_c="c.100A>T",
                hgvs_p="p.Lys34Ter",
                pathogenicity_scores={"CADD_phred": 34.0, "REVEL": 0.91},
                source_tool="synthetic_vep",
                source_tool_version="0.0.1",
            ),
        ),
        population_frequencies=(
            PopulationFrequency(
                source="gnomAD_genomes",
                version="v4.1.0",
                population="global",
                allele_frequency=0.000_01,
                allele_count=2,
                allele_number=152_000,
            ),
        ),
        qc_flags=("synthetic_fixture",),
        source_artifact="synthetic_case.vcf",
        source_line_index=7,
        normalisation_ops=("split_multiallelic", "left_align"),
    )


def make_evidence(
    *,
    subject_id: str,
    claim: str,
    category: EvidenceCategory = EvidenceCategory.CONSEQUENCE,
    direction: EvidenceDirection = EvidenceDirection.SUPPORTS,
    tier: AssertionTier = AssertionTier.COMPUTATIONAL_PREDICTION,
    evidence_type: EvidenceType = EvidenceType.IN_SILICO_PREDICTION,
    citation: Citation | None = None,
    run_id: str | None = RUN_ID,
) -> EvidenceItem:
    return EvidenceItem(
        evidence_id=make_evidence_id(
            subject_id=subject_id, category=category, claim=claim, tool="synthetic_vep"
        ),
        subject_id=subject_id,
        subject_kind="variant",
        claim=claim,
        category=category,
        direction=direction,
        strength=EvidenceStrength.MODERATE,
        evidence_type=evidence_type,
        tier=tier,
        citation=citation,
        method="Deterministic synthetic annotation over a fixed table.",
        tool="synthetic_vep",
        tool_version="0.0.1",
        limitations="Synthetic substitute; not biologically validated.",
        timestamp=INSTANT,
        run_id=run_id,
        numeric_value=34.0,
        payload={"transcript": "ENST00000000001", "rank": 1, "flagged": True, "note": None},
    )


def make_pair(
    *,
    variant_a: VariantRecord,
    variant_b: VariantRecord | None,
    composite_score: float,
    gene: str = "SYNTHA",
) -> CandidatePair:
    variant_ids = (
        (variant_a.variant_id,)
        if variant_b is None
        else (variant_a.variant_id, variant_b.variant_id)
    )
    return CandidatePair(
        pair_id=make_pair_id(gene, variant_ids),
        gene_symbol=gene,
        variant_a=variant_a,
        variant_b=variant_b,
        inheritance_model=(
            InheritanceModel.COMPOUND_HETEROZYGOUS
            if variant_b is not None
            else InheritanceModel.DE_NOVO_DOMINANT
        ),
        phase=PhaseEvidence(status=PhaseStatus.UNKNOWN, method="none"),
        scores=ComponentScores(
            analytical_validity=0.9,
            rarity=0.95,
            molecular_consequence=0.8,
            inheritance_consistency=0.5,
            phenotype_similarity=0.7,
            mechanistic_relevance=0.6,
            evidence_quality=0.5,
            contradiction_penalty=0.1,
        ),
        composite_score=composite_score,
        supporting_evidence_ids=("EV-CONS-aaaa", "EV-CONS-bbbb"),
        contradicting_evidence_ids=("EV-CONT-cccc",),
        recommended_next_test="Trio sequencing to establish phase.",
        rank_rationale="Synthetic fixture.",
        flags=("phase_unknown",),
    )


def make_drug(*, rejected: bool) -> DrugHypothesis:
    """A wrong-direction agent (which MUST be rejected) or a compliant one."""
    return DrugHypothesis(
        drug_id="DRUG-REJECT-001" if rejected else "DRUG-OK-001",
        name="Synthetic MPS1 inhibitor" if rejected else "Synthetic antioxidant",
        approved_name=None if rejected else "synthoxidant",
        approval_status=(
            ApprovalStatus.TOOL_COMPOUND if rejected else ApprovalStatus.APPROVED_OTHER_INDICATION
        ),
        intervention_class=InterventionClass.DISEASE_MODIFYING,
        target="TTK/MPS1 kinase" if rejected else "Mitochondrial ROS",
        target_node_id="node-checkpoint",
        mechanism_of_action="Synthetic fixture mechanism of action.",
        required_direction=EffectDirection.RESTORE,
        # A checkpoint inhibitor pushes the target the same way the disease does.
        observed_direction=(EffectDirection.DECREASE if rejected else EffectDirection.STABILISE),
        is_direct_evidence=False,
        strongest_evidence_type=EvidenceType.CELL_LINE,
        pediatric_evidence=PediatricEvidence(has_pediatric_exposure=False),
        pharmacokinetics=PharmacokineticProfile(route="oral", cns_penetrant=None),
        worsens_chromosomal_instability=True if rejected else None,
        proposed_validation_experiment="Micronucleus assay in patient-derived fibroblasts.",
        score=0.0 if rejected else 0.4,
        rejected=rejected,
        rejection_reasons=(
            (RejectionReason.WRONG_DIRECTION, RejectionReason.ONCOGENIC_RISK) if rejected else ()
        ),
        rejection_rationale=(
            "Acts in the disease direction on an already-deficient checkpoint." if rejected else ""
        ),
    )


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def store(tmp_path: Path) -> Iterator[EvidenceStore]:
    with EvidenceStore(tmp_path / "evidence.duckdb") as opened:
        opened.initialise()
        yield opened


# ----------------------------------------------------------------------- tests


@pytest.mark.unit
def test_schema_initialises_empty(store: EvidenceStore) -> None:
    """Every declared table exists and starts at zero rows."""
    counts = store.counts()
    assert set(counts) == set(TABLES)
    assert all(count == 0 for count in counts.values()), counts
    # initialise() is idempotent: re-applying must not raise or repopulate.
    store.initialise()
    assert store.counts() == counts


@pytest.mark.unit
def test_evidence_round_trip_preserves_every_field(store: EvidenceStore) -> None:
    """A read reconstructs the identical model, enums, payload and timestamp."""
    original = make_evidence(
        subject_id="GRCh38:chr15:100000:A:T",
        claim="Stop-gained on the MANE-Select transcript.",
        tier=AssertionTier.DATABASE_ASSERTION,
        evidence_type=EvidenceType.CURATED_DATABASE,
        citation=Citation(
            source="ClinVar",
            identifier="VCV000000001",
            version="2026-01",
            url="https://example.invalid/vcv1",
            title="Synthetic record",
        ),
    )
    assert store.write_evidence([original]) == 1

    restored = store.evidence_for("GRCh38:chr15:100000:A:T")
    assert len(restored) == 1
    item = restored[0]

    assert item == original
    # Spelled out, because "==" on a pydantic model is easy to satisfy accidentally
    # if a field silently became a string.
    assert item.category is EvidenceCategory.CONSEQUENCE
    assert item.direction is EvidenceDirection.SUPPORTS
    assert item.strength is EvidenceStrength.MODERATE
    assert item.evidence_type is EvidenceType.CURATED_DATABASE
    assert item.tier is AssertionTier.DATABASE_ASSERTION
    assert item.citation is not None
    assert item.citation.version == "2026-01"
    assert item.timestamp == INSTANT
    assert item.payload == {
        "transcript": "ENST00000000001",
        "rank": 1,
        "flagged": True,
        "note": None,
    }
    assert item.limitations == "Synthetic substitute; not biologically validated."
    # The citation was also normalised into the bibliography table.
    assert store.counts()["citations"] == 1


@pytest.mark.unit
def test_writers_are_idempotent(store: EvidenceStore) -> None:
    """Writing the same input twice converges instead of accumulating."""
    variants = [make_variant(), make_variant(position=200_000, ref="G", alt="C")]

    first = store.write_variants(variants)
    store.write_consequences(variants)
    store.write_frequencies(variants)
    counts_after_first = store.counts()

    second = store.write_variants(variants)
    store.write_consequences(variants)
    store.write_frequencies(variants)
    counts_after_second = store.counts()

    assert first == second == 2
    assert counts_after_first == counts_after_second
    assert counts_after_first["variants"] == 2
    assert counts_after_first["consequences"] == 2
    assert counts_after_first["frequencies"] == 2
    assert counts_after_first["genes"] == 1

    rows = store.query("SELECT variant_id, qc_flags FROM variants ORDER BY variant_id")
    assert len(rows) == 2
    assert rows[0]["qc_flags"] == ["synthetic_fixture"]


@pytest.mark.unit
def test_contradicting_evidence_is_preserved(store: EvidenceStore) -> None:
    """GP-19: what argues against a hypothesis survives the round trip."""
    subject = "PAIR-SYNTHA-abc123"
    supporting = make_evidence(
        subject_id=subject,
        claim="Both variants are rare in every recorded population.",
        category=EvidenceCategory.POPULATION,
    )
    contradicting = make_evidence(
        subject_id=subject,
        claim="Read-backed phasing places both variants on the same haplotype.",
        category=EvidenceCategory.CONTRADICTION,
        direction=EvidenceDirection.CONTRADICTS,
        tier=AssertionTier.OBSERVED_DATA,
        evidence_type=EvidenceType.DIRECT_MEASUREMENT,
    )
    store.write_evidence([supporting, contradicting])

    assert len(store.evidence_for(subject)) == 2

    contradictions = store.contradictions_for(subject)
    assert len(contradictions) == 1
    assert contradictions[0] == contradicting
    assert contradictions[0].direction is EvidenceDirection.CONTRADICTS
    assert contradictions[0].is_contradiction

    # And the ledger agrees with the store.
    ledger = EvidenceLedger(run_id=RUN_ID)
    ledger.extend([supporting, contradicting])
    assert ledger.contradictions() == (contradicting,)
    assert len(ledger) == 2


@pytest.mark.unit
def test_assertion_resolver_enforces_gp10() -> None:
    """An unsourced or dangling claim never reaches a report."""
    ledger = EvidenceLedger(run_id=RUN_ID)
    known = ledger.add(
        make_evidence(
            subject_id="GRCh38:chr15:100000:A:T",
            claim="Predicted loss of function.",
            run_id=None,
        )
    )
    resolver = AssertionResolver(ledger)

    assert known.run_id == RUN_ID  # stamped by the ledger

    with pytest.raises(UnsourcedAssertionError, match="cites no evidence"):
        resolver.require("This gene causes the phenotype.", [])

    with pytest.raises(UnsourcedAssertionError, match="do not resolve"):
        resolver.require("This gene causes the phenotype.", ["EV-CONS-deadbeef"])

    resolved = resolver.require("Predicted loss of function.", [known.evidence_id])
    assert resolved == (known,)
    # resolve() is best-effort and does not raise on an unknown ID.
    assert resolver.resolve(["EV-CONS-deadbeef"]) == ()


@pytest.mark.unit
def test_export_parquet_is_byte_identical(store: EvidenceStore, tmp_path: Path) -> None:
    """GP-30: two exports of the same data hash identically, file by file."""
    variants = [make_variant(), make_variant(position=200_000, ref="G", alt="C")]
    store.write_variants(variants)
    store.write_consequences(variants)
    store.write_frequencies(variants)
    store.write_evidence(
        [
            make_evidence(subject_id=variants[0].variant_id, claim="Stop gained."),
            make_evidence(
                subject_id=variants[1].variant_id,
                claim="Splice region variant of uncertain effect.",
                direction=EvidenceDirection.CONTRADICTS,
            ),
        ]
    )
    store.write_pairs(
        [make_pair(variant_a=variants[0], variant_b=variants[1], composite_score=0.8)]
    )
    store.write_edges(
        [
            GraphEdge(
                subject_id=variants[0].variant_id,
                subject_kind="variant",
                predicate="disrupts",
                object_id="SYNTHA",
                object_kind="gene",
                evidence_ids=("EV-CONS-0001",),
                confidence=0.9,
            )
        ]
    )
    store.write_phenotypes(
        PhenotypeProfile(
            subject_id="CASE-001",
            observations=(
                PhenotypeObservation(
                    hpo_id="HP:0001250",
                    label="Seizure",
                    status=ObservationStatus.OBSERVED,
                    provenance="synthetic_fixture",
                    extraction_confidence=0.9,
                ),
            ),
            source_artifact="synthetic_phenotype.tsv",
            hpo_version="2026-01-01",
        )
    )

    first = store.export_parquet(tmp_path / "export-a")
    second = store.export_parquet(tmp_path / "export-b")

    assert set(first) == set(TABLES)
    assert set(second) == set(TABLES)
    differing = [table for table in TABLES if sha256_of(first[table]) != sha256_of(second[table])]
    assert not differing, f"non-deterministic parquet for: {differing}"

    # Re-exporting into the same directory must also be stable, which is the case
    # a re-run actually hits.
    again = store.export_parquet(tmp_path / "export-a")
    assert sha256_of(again["variants"]) == sha256_of(first["variants"])

    # A selective export writes only what was asked for.
    selective = store.export_parquet(tmp_path / "export-c", tables=["evidence_items"])
    assert set(selective) == {"evidence_items"}


@pytest.mark.unit
def test_graph_edges_round_trip(store: EvidenceStore) -> None:
    """The relational graph keeps its list column and both traversal directions."""
    edges = [
        GraphEdge(
            subject_id="GRCh38:chr15:100000:A:T",
            subject_kind="variant",
            predicate="disrupts",
            object_id="SYNTHA",
            object_kind="gene",
            evidence_ids=("EV-CONS-0001", "EV-CONS-0002"),
            confidence=0.85,
        ),
        GraphEdge(
            subject_id="SYNTHA",
            subject_kind="gene",
            predicate="associated_with",
            object_id="HP:0001250",
            object_kind="phenotype",
        ),
    ]
    assert store.write_edges(edges, run_id=RUN_ID) == 2
    # Idempotent on the derived edge_id.
    assert store.write_edges(edges, run_id=RUN_ID) == 2
    assert store.counts()["graph_edges"] == 2

    rows = store.query(
        "SELECT * FROM graph_edges WHERE subject_id = ? ORDER BY predicate",
        ("GRCh38:chr15:100000:A:T",),
    )
    assert len(rows) == 1
    row = rows[0]
    assert row["edge_id"] == edges[0].edge_id
    assert row["predicate"] == "disrupts"
    assert row["object_kind"] == "gene"
    assert row["evidence_ids"] == ["EV-CONS-0001", "EV-CONS-0002"]
    assert row["confidence"] == pytest.approx(0.85)
    assert row["run_id"] == RUN_ID

    # An edge with no confidence stores NULL, not an invented 1.0.
    reverse = store.query(
        "SELECT confidence, evidence_ids FROM graph_edges WHERE object_id = ?", ("HP:0001250",)
    )
    assert reverse == [{"confidence": None, "evidence_ids": []}]


@pytest.mark.unit
def test_rejected_drugs_persist_with_reasons(store: EvidenceStore) -> None:
    """GP-19: a wrong-direction agent is kept, with every reason it was dropped."""
    rejected = make_drug(rejected=True)
    accepted = make_drug(rejected=False)
    assert rejected.directions_agree is False
    assert store.write_drugs([rejected, accepted], run_id=RUN_ID) == 2

    counts = store.counts()
    assert counts["drugs"] == 2
    assert counts["drug_rejections"] == 2  # one row per reason

    rows = store.query(
        "SELECT reason, drug_name, rationale, required_direction, observed_direction "
        "FROM drug_rejections WHERE drug_id = ? ORDER BY reason",
        (rejected.drug_id,),
    )
    assert [row["reason"] for row in rows] == ["oncogenic_risk", "wrong_direction"]
    assert all(row["rationale"] for row in rows)
    assert rows[0]["required_direction"] == "restore"
    assert rows[0]["observed_direction"] == "decrease"

    stored = store.query(
        "SELECT drug_id, rejected, directions_agree, worsens_chromosomal_instability "
        "FROM drugs ORDER BY drug_id"
    )
    by_id = {row["drug_id"]: row for row in stored}
    assert by_id[rejected.drug_id]["rejected"] is True
    assert by_id[rejected.drug_id]["directions_agree"] is False
    # Tri-state survives: the accepted drug is unassessed, not "safe".
    assert by_id[accepted.drug_id]["worsens_chromosomal_instability"] is None

    assert store.write_drugs([rejected, accepted], run_id=RUN_ID) == 2
    assert store.counts() == counts


@pytest.mark.unit
def test_ranked_pairs_orders_and_limits(store: EvidenceStore) -> None:
    """Ranking is by composite score descending, and the limit is respected."""
    variant_a = make_variant()
    variant_b = make_variant(position=200_000, ref="G", alt="C")
    variant_c = make_variant(contig="chr7", position=50_000, gene="SYNTHB")

    pairs = [
        make_pair(variant_a=variant_a, variant_b=variant_b, composite_score=0.42),
        make_pair(variant_a=variant_c, variant_b=None, composite_score=0.91, gene="SYNTHB"),
        make_pair(variant_a=variant_b, variant_b=None, composite_score=0.67, gene="SYNTHC"),
    ]
    assert store.write_pairs(pairs) == 3

    ranked = store.ranked_pairs()
    assert [row["composite_score"] for row in ranked] == [0.91, 0.67, 0.42]

    top = store.ranked_pairs(limit=2)
    assert len(top) == 2
    assert [row["composite_score"] for row in top] == [0.91, 0.67]
    assert top[0]["gene_symbol"] == "SYNTHB"
    assert top[0]["contradicting_evidence_ids"] == ["EV-CONT-cccc"]
    assert store.ranked_pairs(limit=0) == []
