"""Compatibility tests for the Track 1 submission format.

Every coordinate in this file is **fabricated**. Nothing here is derived from
patient data, and the file is safe to commit (GP-40).

The contract implemented is `docs/references/track1-submission-contract.md`. These
tests are a self-check against our reading of the challenge scorer, and they exist
because the failure mode they guard against is silent: a submission with a bare
contig name or a zero EPCR uploads cleanly, looks correct to a human, scores
nothing, and burns one of six attempts.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path

import pytest

from mva.models import (
    CandidatePair,
    ComponentScores,
    FilterStatus,
    GenomeBuild,
    GenomicCoordinate,
    Genotype,
    InheritanceModel,
    PhaseEvidence,
    PhaseStatus,
    VariantRecord,
    Zygosity,
)
from mva.reporting.track1 import (
    ACCEPTED_PROBAND_ID,
    EPCR_FLOOR,
    MAX_SUBMISSION_ROWS,
    TRACK1_COLUMNS,
    SubmissionRow,
    build_submission_rows,
    composite_to_epcr,
    render_submission_csv,
    render_submission_csv_unvalidated,
    truncation_notice,
    validate_submission,
    write_submission,
)

#: The header line the contract specifies, verbatim.
EXPECTED_HEADER = (
    "proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr,finding_type,notes"
)

SECOND_VARIANT_COLUMNS = ("chrom_2", "pos_2", "ref_2", "alt_2")


def _variant(contig: str, position: int, ref: str, alt: str) -> VariantRecord:
    return VariantRecord(
        coordinate=GenomicCoordinate(
            build=GenomeBuild.GRCH38, contig=contig, position=position, ref=ref, alt=alt
        ),
        genotype=Genotype(zygosity=Zygosity.HET, genotype_string="0/1", depth=30),
        filter_status=FilterStatus.PASS,
        source_artifact="fabricated-fixture",
    )


def _scores() -> ComponentScores:
    return ComponentScores(
        analytical_validity=0.5,
        rarity=0.5,
        molecular_consequence=0.5,
        inheritance_consistency=0.5,
        phenotype_similarity=0.5,
        mechanistic_relevance=0.5,
        evidence_quality=0.5,
        contradiction_penalty=0.0,
    )


def _pair(
    *,
    pair_id: str,
    gene: str,
    composite: float,
    contig: str = "chr7",
    position: int = 100_000,
    second: bool = True,
) -> CandidatePair:
    """A fabricated candidate. ``second=False`` gives a single-variant proposal."""
    variant_b = _variant(contig, position + 5_000, "G", "C") if second else None
    return CandidatePair(
        pair_id=pair_id,
        gene_symbol=gene,
        variant_a=_variant(contig, position, "A", "T"),
        variant_b=variant_b,
        inheritance_model=(
            InheritanceModel.COMPOUND_HETEROZYGOUS if second else InheritanceModel.DE_NOVO_DOMINANT
        ),
        phase=PhaseEvidence(status=PhaseStatus.UNKNOWN, method="none"),
        scores=_scores(),
        composite_score=composite,
        recommended_next_test="Fabricated next test for a fabricated candidate.",
    )


def _parse(csv_text: str) -> list[dict[str, str]]:
    reader = csv.DictReader(io.StringIO(csv_text, newline=""))
    assert reader.fieldnames is not None
    return [{key: (value or "") for key, value in row.items()} for row in reader]


# ---------------------------------------------------------------------------
# 10. Round-trip against the public format.
# ---------------------------------------------------------------------------


@pytest.mark.integration
def test_rendered_csv_round_trips_through_dictreader() -> None:
    rows = build_submission_rows(
        [
            _pair(pair_id="P1", gene="FAKEGENE1", composite=0.90, contig="chr1", position=100_000),
            _pair(
                pair_id="P2",
                gene="FAKEGENE2",
                composite=0.30,
                contig="chr7",
                position=300_000,
                second=False,
            ),
        ],
        proband_id=ACCEPTED_PROBAND_ID,
    )
    text = render_submission_csv(rows)

    assert text.splitlines()[0] == EXPECTED_HEADER
    parsed = _parse(text)
    assert len(parsed) == 2
    assert len(parsed) <= MAX_SUBMISSION_ROWS

    for row in parsed:
        assert tuple(row) == TRACK1_COLUMNS
        assert row["proband_id"] == "PROBAND01"
        assert int(row["pos_1"]) > 0
        assert 0 < float(row["epcr"]) <= 1
        assert row["finding_type"] in {"primary", "secondary"}

    compound, single = parsed
    assert all(compound[column] for column in SECOND_VARIANT_COLUMNS), (
        "a compound-het pair is ONE row using the _2 columns, not two rows"
    )
    assert int(compound["pos_2"]) > 0
    assert all(single[column] == "" for column in SECOND_VARIANT_COLUMNS)

    ok, errors = validate_submission(text)
    assert ok, errors


@pytest.mark.integration
def test_header_is_present_and_exact() -> None:
    text = render_submission_csv([])
    assert text == EXPECTED_HEADER + "\n"
    assert tuple(EXPECTED_HEADER.split(",")) == TRACK1_COLUMNS
    assert len(TRACK1_COLUMNS) == 12


# ---------------------------------------------------------------------------
# 11. The chr prefix — the single most dangerous detail in the contract.
# ---------------------------------------------------------------------------


@pytest.mark.integration
def test_bare_contig_input_still_renders_the_chr_prefix() -> None:
    """Ensembl-style '15' must never reach the file; the scorer compares raw strings."""
    rows = build_submission_rows(
        [_pair(pair_id="P1", gene="FAKEGENE1", composite=0.5, contig="15", position=40_200_000)],
        proband_id=ACCEPTED_PROBAND_ID,
    )
    assert rows[0].chrom_1 == "chr15"
    assert rows[0].chrom_2 == "chr15"

    parsed = _parse(render_submission_csv(rows))
    for row in parsed:
        for column in ("chrom_1", "chrom_2"):
            assert row[column] == "" or row[column].startswith("chr")


@pytest.mark.integration
def test_validation_rejects_a_chromosome_without_the_prefix() -> None:
    """A hand-built row bypassing the builder must still be caught."""
    bare = SubmissionRow(
        proband_id=ACCEPTED_PROBAND_ID,
        chrom_1="15",
        pos_1=40_200_000,
        ref_1="C",
        alt_1="T",
        chrom_2="",
        pos_2="",
        ref_2="",
        alt_2="",
        epcr=0.5,
    )
    # Rendered through the unvalidated path on purpose: `render_submission_csv`
    # now refuses invalid rows outright (ADR 0023), and a validator can only be
    # shown to work on bytes that fail it.
    ok, errors = validate_submission(render_submission_csv_unvalidated([bare]))
    assert not ok
    assert any("chr" in error and "prefix" in error for error in errors)


# ---------------------------------------------------------------------------
# 12. EPCR bounds.
# ---------------------------------------------------------------------------


@pytest.mark.integration
def test_composite_to_epcr_stays_inside_the_accepted_range() -> None:
    assert composite_to_epcr(0.0) > 0, "the scorer rejects epcr == 0 outright"
    assert composite_to_epcr(0.0) == EPCR_FLOOR
    assert composite_to_epcr(1.0) <= 1.0
    assert composite_to_epcr(1.0) == 1.0
    for composite in (0.0, 0.01, 0.25, 0.5, 0.75, 0.999, 1.0):
        assert 0 < composite_to_epcr(composite) <= 1


@pytest.mark.integration
def test_composite_to_epcr_is_strictly_increasing() -> None:
    """Rank order is the whole prediction; the map must not reorder it."""
    values = [composite_to_epcr(step / 100) for step in range(101)]
    assert values == sorted(values)
    assert values[0] < values[-1]


@pytest.mark.integration
def test_a_zero_scored_candidate_still_produces_a_valid_row() -> None:
    rows = build_submission_rows(
        [_pair(pair_id="P0", gene="FAKEGENE0", composite=0.0)], proband_id=ACCEPTED_PROBAND_ID
    )
    ok, errors = validate_submission(render_submission_csv(rows))
    assert ok, errors
    assert float(_parse(render_submission_csv(rows))[0]["epcr"]) > 0


# ---------------------------------------------------------------------------
# 13-14. Ordering and the ten-row limit.
# ---------------------------------------------------------------------------


@pytest.mark.integration
def test_rows_are_sorted_by_epcr_descending() -> None:
    rows = build_submission_rows(
        [
            _pair(pair_id="P1", gene="G1", composite=0.10, contig="chr1", position=1_000),
            _pair(pair_id="P2", gene="G2", composite=0.90, contig="chr2", position=2_000),
            _pair(pair_id="P3", gene="G3", composite=0.50, contig="chr3", position=3_000),
        ],
        proband_id=ACCEPTED_PROBAND_ID,
    )
    epcrs = [float(row["epcr"]) for row in _parse(render_submission_csv(rows))]
    assert epcrs == sorted(epcrs, reverse=True)
    assert [row.notes.split(";")[0] for row in rows] == ["G2", "G3", "G1"]


@pytest.mark.integration
def test_more_than_ten_candidates_truncate_to_exactly_ten() -> None:
    pairs = [
        _pair(
            pair_id=f"P{index}",
            gene=f"GENE{index}",
            composite=1.0 - index / 100,
            contig="chr1",
            position=1_000 + index * 1_000,
        )
        for index in range(13)
    ]
    rows = build_submission_rows(pairs, proband_id=ACCEPTED_PROBAND_ID)
    assert len(rows) == MAX_SUBMISSION_ROWS

    parsed = _parse(render_submission_csv(rows))
    assert len(parsed) == 10
    # The ten kept are the ten highest-scoring, in order.
    assert [row["notes"].split(";")[0] for row in parsed] == [f"GENE{i}" for i in range(10)]

    notice = truncation_notice(len(pairs))
    assert notice is not None
    assert "13" in notice and "10" in notice
    assert truncation_notice(MAX_SUBMISSION_ROWS) is None


@pytest.mark.integration
def test_validation_rejects_an_eleven_row_file() -> None:
    rows = build_submission_rows(
        [
            _pair(
                pair_id=f"P{index}",
                gene=f"GENE{index}",
                composite=0.5,
                contig="chr1",
                position=1_000 + index * 1_000,
            )
            for index in range(11)
        ],
        proband_id=ACCEPTED_PROBAND_ID,
        max_rows=MAX_SUBMISSION_ROWS,
    )
    assert len(rows) == 10
    eleven = render_submission_csv_unvalidated([*rows, rows[0]])
    ok, errors = validate_submission(eleven)
    assert not ok
    assert any("11 data rows" in error for error in errors)


# ---------------------------------------------------------------------------
# Privacy and the write path.
# ---------------------------------------------------------------------------


@pytest.mark.integration
def test_notes_carry_only_gene_and_mechanism_level_statements() -> None:
    pair = _pair(pair_id="P1", gene="FAKEGENE1", composite=0.7).model_copy(
        update={"rank_rationale": "Child presented at 14 months with seizures; MRN 123456."}
    )
    rows = build_submission_rows([pair], proband_id=ACCEPTED_PROBAND_ID)
    note = rows[0].notes
    assert "seizures" not in note
    assert "123456" not in note
    assert note.startswith("FAKEGENE1")
    assert "phase=unknown" in note


@pytest.mark.integration
def test_write_submission_round_trips_to_disk(tmp_path: Path) -> None:
    rows = build_submission_rows(
        [_pair(pair_id="P1", gene="FAKEGENE1", composite=0.8)], proband_id=ACCEPTED_PROBAND_ID
    )
    # The filename is load-bearing, not incidental: write_submission now acts on
    # the public-export gate's verdict, and the gate's allowlist is deny-by-default
    # (mva.privacy.export.PUBLIC_EXPORT_ALLOWLIST). A submission written under any
    # other name is refused and deleted.
    path = write_submission(rows, tmp_path / "nested" / "track1_submission.csv")
    text = path.read_text(encoding="utf-8")
    assert text == render_submission_csv(rows)
    ok, errors = validate_submission(text)
    assert ok, errors


@pytest.mark.integration
def test_write_submission_refuses_an_invalid_file(tmp_path: Path) -> None:
    """A malformed submission that exists is worse than one that does not."""
    bad = SubmissionRow(
        proband_id="NOT-THE-PROBAND",
        chrom_1="chr1",
        pos_1=1_000,
        ref_1="A",
        alt_1="T",
        chrom_2="",
        pos_2="",
        ref_2="",
        alt_2="",
        epcr=0.5,
    )
    target = tmp_path / "submission.csv"
    with pytest.raises(ValueError, match="contract check"):
        write_submission([bad], target)
    assert not target.exists()
