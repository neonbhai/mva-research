"""Submission shaping.

Covers finding F-01, raised by the prioritisation stage during integration: the
ranked list legitimately contains single-variant hypotheses carved out of a
higher-ranked pair, and submitting both wastes one of the ten available rows.
"""

from __future__ import annotations

import pytest

from mva.models.genome import GenomeBuild, GenomicCoordinate
from mva.models.pair import (
    CandidatePair,
    ComponentScores,
    InheritanceModel,
    PhaseEvidence,
    PhaseStatus,
    make_pair_id,
)
from mva.models.variant import (
    ConsequenceAnnotation,
    FilterStatus,
    Genotype,
    ImpactSeverity,
    VariantRecord,
    Zygosity,
)
from mva.reporting.track1 import build_submission_rows

pytestmark = pytest.mark.integration


def _variant(contig: str, position: int, ref: str, alt: str) -> VariantRecord:
    return VariantRecord(
        coordinate=GenomicCoordinate(
            build=GenomeBuild.GRCH38, contig=contig, position=position, ref=ref, alt=alt
        ),
        genotype=Genotype(
            zygosity=Zygosity.HET,
            genotype_string="0/1",
            depth=45,
            ref_reads=23,
            alt_reads=22,
            genotype_quality=99,
        ),
        filter_status=FilterStatus.PASS,
        consequences=(
            ConsequenceAnnotation(
                gene_symbol="SYNTHKIN1",
                transcript_id="SYNTHT0001.1",
                consequence_terms=("stop_gained",),
                impact=ImpactSeverity.HIGH,
            ),
        ),
        source_artifact="test",
    )


def _pair(
    variant_a: VariantRecord, variant_b: VariantRecord | None, composite: float, rank: int
) -> CandidatePair:
    variants = (
        (variant_a.variant_id,)
        if variant_b is None
        else (variant_a.variant_id, variant_b.variant_id)
    )
    return CandidatePair(
        pair_id=make_pair_id("SYNTHKIN1", variants),
        gene_symbol="SYNTHKIN1",
        variant_a=variant_a,
        variant_b=variant_b,
        inheritance_model=(
            InheritanceModel.COMPOUND_HETEROZYGOUS
            if variant_b is not None
            else InheritanceModel.UNKNOWN
        ),
        phase=PhaseEvidence(status=PhaseStatus.UNKNOWN, method="none"),
        scores=ComponentScores(
            analytical_validity=0.9,
            rarity=0.9,
            molecular_consequence=0.9,
            inheritance_consistency=0.5,
            phenotype_similarity=0.8,
            mechanistic_relevance=0.7,
            evidence_quality=0.6,
            contradiction_penalty=0.0,
        ),
        composite_score=composite,
        rank=rank,
        recommended_next_test="Parental segregation testing.",
    )


def test_no_submission_row_is_subsumed_by_an_earlier_row() -> None:
    """A row already covered by an earlier row adds nothing the scorer can use.

    F-max unions variants across rows (`predicted_variants |= row.variants`), so a
    subset row contributes no new variant, and rank points go to the first full
    match. It only consumes one of the ten rows.
    """
    a = _variant("chr15", 40_200_000, "C", "T")
    b = _variant("chr15", 40_210_500, "G", "A")
    top = _pair(a, b, composite=0.90, rank=1)
    subset = _pair(a, None, composite=0.80, rank=2)
    other = _pair(
        _variant("chr7", 55_000_000, "C", "G"),
        _variant("chr7", 55_001_000, "A", "T"),
        composite=0.70,
        rank=3,
    )

    rows = build_submission_rows([top, subset, other], proband_id="PROBAND01")

    seen: list[frozenset[str]] = []
    for row in rows:
        variants = {f"{row.chrom_1}:{row.pos_1}:{row.ref_1}:{row.alt_1}"}
        if row.chrom_2:
            variants.add(f"{row.chrom_2}:{row.pos_2}:{row.ref_2}:{row.alt_2}")
        frozen = frozenset(variants)
        assert not any(frozen <= earlier for earlier in seen), (
            "a submission row is subsumed by an earlier row; it consumes one of the "
            "ten available rows without adding a variant the scorer can use (F-01)"
        )
        seen.append(frozen)

    assert len(rows) == 2, "the subsumed single-variant row should have been dropped"
    assert any(row.chrom_1 == "chr7" for row in rows), (
        "the distinct second-gene hypothesis must survive the dedup"
    )


def test_dedup_does_not_drop_a_partially_overlapping_pair() -> None:
    """Only strict subsets are dropped.

    Two pairs sharing one variant but differing in the other are genuinely
    different hypotheses, and both must reach the submission.
    """
    shared = _variant("chr15", 40_200_000, "C", "T")
    first = _pair(shared, _variant("chr15", 40_210_500, "G", "A"), composite=0.90, rank=1)
    second = _pair(shared, _variant("chr15", 40_211_000, "G", "T"), composite=0.85, rank=2)

    rows = build_submission_rows([first, second], proband_id="PROBAND01")
    assert len(rows) == 2, "a partially overlapping pair is a distinct hypothesis"
