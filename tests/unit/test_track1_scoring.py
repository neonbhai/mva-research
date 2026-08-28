"""What the submission is actually WORTH, computed with the challenge's own arithmetic.

Every other Track 1 test asserts a *shape*: ten rows, a `chr` prefix, distinct
EPCRs, a pair above the single it subsumes. Shape tests are necessary and they are
not sufficient, because the thing that decides the competition is a number produced
by `evaluation.py` — and a composition change can preserve every shape assertion in
the suite while moving that number from 100 to 50.

So this module carries a **replica of the scorer** transcribed from
`docs/references/track1-submission-contract.md` (re-verified against the real
`evaluation.py` on 2026-08-28), and the defect tests below assert POINTS, not
layout. :func:`test_the_replica_reproduces_every_verified_scorer_result` pins the
replica against the seven submission shapes the contract records as executed
results, so a wrong replica fails loudly here instead of quietly blessing a
regression somewhere else.

The replica is deliberately a transcription rather than an import: there is nothing
to import (the Space is gated), and a paraphrase that "looks equivalent" is exactly
how a scoring assumption drifts away from the scorer.
"""

from __future__ import annotations

import csv
import io
from collections.abc import Sequence
from dataclasses import replace
from pathlib import Path

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
from mva.reporting.track1 import (
    ACCEPTED_PROBAND_ID,
    MAX_SUBMISSION_ROWS,
    SubmissionRow,
    build_submission_rows,
    render_submission_csv,
    render_submission_csv_unvalidated,
    validate_submission,
    write_submission,
)

pytestmark = pytest.mark.unit

GENE = "SYNTHKIN1"

#: The scorer's variant key: ``(chrom, pos, ref, alt)``. ``chrom`` is only
#: ``.strip()``-ed — there is no chr-prefix normalisation on either side.
VariantKey = tuple[str, int, str, str]

#: ``[(rank_ceiling, points)]`` exactly as `evaluation.py` holds them.
RANK_TIERS: tuple[tuple[int, int], ...] = ((1, 100), (3, 50), (5, 25), (10, 10))


# ---------------------------------------------------------------------------
# The replica
# ---------------------------------------------------------------------------


def _tier_points(rank: int) -> int:
    for ceiling, points in RANK_TIERS:
        if rank <= ceiling:
            return points
    return 0


def _parse_variant(chrom: str, pos: str, ref: str, alt: str) -> VariantKey:
    """`_parse_variant` from the scorer: strip the contig, upper-case the alleles."""
    return (chrom.strip(), int(pos), ref.strip().upper(), alt.strip().upper())


def parse_submission(csv_text: str) -> tuple[tuple[frozenset[VariantKey], float], ...]:
    """Read rendered CSV the way the scorer does: (variant set, epcr) in file order."""
    reader = csv.DictReader(io.StringIO(csv_text, newline=""))
    rows: list[tuple[frozenset[VariantKey], float]] = []
    for row in reader:
        keys = {
            _parse_variant(
                row["chrom_1"] or "", row["pos_1"] or "", row["ref_1"] or "", row["alt_1"] or ""
            )
        }
        if (row.get("chrom_2") or "").strip():
            keys.add(
                _parse_variant(
                    row["chrom_2"] or "", row["pos_2"] or "", row["ref_2"] or "", row["alt_2"] or ""
                )
            )
        rows.append((frozenset(keys), float(row["epcr"] or "")))
    return tuple(rows)


def rank_points(
    rows: Sequence[tuple[frozenset[VariantKey], float]], truth: frozenset[VariantKey]
) -> float:
    """`score_proband`'s rank component.

    Rank is ``sorted(enumerate(rows), key=lambda x: (-x[1][1], x[0]))`` — EPCR
    descending, ties broken by file order. The first FULL match (frozenset
    equality) wins and the loop breaks. Partial credit is searched **only if no
    full match exists at all**, and only when the answer is a compound het; it
    scores ``0.5 x`` the tier points of the first intersecting row.
    """
    ordered = [row for _, row in sorted(enumerate(rows), key=lambda x: (-x[1][1], x[0]))]
    for rank, (variants, _epcr) in enumerate(ordered, start=1):
        if variants == truth:
            return float(_tier_points(rank))
    if len(truth) == 2:
        for rank, (variants, _epcr) in enumerate(ordered, start=1):
            if variants & truth:
                return 0.5 * _tier_points(rank)
    return 0.0


def f_max(
    rows: Sequence[tuple[frozenset[VariantKey], float]], truth: frozenset[VariantKey]
) -> float:
    """`score_proband`'s F-max component, computed per VARIANT and not per row.

    Thresholds are the unique EPCR values present in the submission itself, swept
    with ``row.epcr >= t`` and unioning ``predicted |= row.variants``. Both metrics
    therefore depend only on the ORDER and TIE STRUCTURE of the EPCR column, never
    on its magnitude.
    """
    best = 0.0
    for threshold in sorted({epcr for _, epcr in rows}, reverse=True):
        predicted: set[VariantKey] = set()
        for variants, epcr in rows:
            if epcr >= threshold:
                predicted |= set(variants)
        hits = len(predicted & truth)
        if not hits:
            continue
        precision = hits / len(predicted)
        recall = hits / len(truth)
        best = max(best, 2 * precision * recall / (precision + recall))
    return best


def score(csv_text: str, truth: frozenset[VariantKey]) -> tuple[float, float]:
    """The pair the challenge reports: (rank points, F-max)."""
    rows = parse_submission(csv_text)
    return rank_points(rows, truth), f_max(rows, truth)


# ---------------------------------------------------------------------------
# The replica is pinned to the contract's executed results
# ---------------------------------------------------------------------------

A: VariantKey = ("chr15", 40_200_000, "C", "T")
B: VariantKey = ("chr15", 40_210_500, "G", "A")
WRONG_1: VariantKey = ("chr7", 55_000_000, "C", "G")
WRONG_2: VariantKey = ("chr7", 55_001_000, "A", "T")
TRUTH: frozenset[VariantKey] = frozenset({A, B})


def test_the_replica_reproduces_every_verified_scorer_result() -> None:
    """The seven rows of the contract's "four results that decide our submission".

    Those numbers were produced by executing the real `evaluation.py`. If this
    replica disagrees with any of them it is the replica that is wrong, and every
    points assertion below it is worthless.
    """
    cases: tuple[tuple[str, list[tuple[frozenset[VariantKey], float]], float, float], ...] = (
        ("true pair alone", [(TRUTH, 0.90)], 100.0, 1.0),
        (
            "true pair + 4 wrong rows at lower epcr",
            [
                (TRUTH, 0.90),
                (frozenset({WRONG_1}), 0.80),
                (frozenset({WRONG_2}), 0.70),
                (frozenset({WRONG_1, WRONG_2}), 0.60),
                (frozenset({("chr3", 10_000_000, "G", "A")}), 0.50),
            ],
            100.0,
            1.0,
        ),
        (
            "true pair tied with a wrong pair, true first",
            [(TRUTH, 0.90), (frozenset({WRONG_1, WRONG_2}), 0.90)],
            100.0,
            2 / 3,
        ),
        (
            "same tie, wrong row first in file order",
            [(frozenset({WRONG_1, WRONG_2}), 0.90), (TRUTH, 0.90)],
            50.0,
            2 / 3,
        ),
        (
            "true pair split into two single-variant rows",
            [(frozenset({A}), 0.90), (frozenset({B}), 0.80)],
            50.0,
            1.0,
        ),
        ("one true + one wrong variant at rank 1", [(frozenset({A, WRONG_1}), 0.90)], 50.0, 0.5),
        (
            "true pair with a bare contig",
            [(frozenset({("15", 40_200_000, "C", "T"), ("15", 40_210_500, "G", "A")}), 0.90)],
            0.0,
            0.0,
        ),
    )
    for label, rows, expected_rank, expected_fmax in cases:
        assert rank_points(rows, TRUTH) == pytest.approx(expected_rank), label
        assert f_max(rows, TRUTH) == pytest.approx(expected_fmax, abs=1e-4), label


# ---------------------------------------------------------------------------
# Candidate construction
# ---------------------------------------------------------------------------


def _variant(
    position: int, *, contig: str = "chr15", ref: str = "C", alt: str = "T"
) -> VariantRecord:
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
                gene_symbol=GENE,
                transcript_id="SYNTHT0001.1",
                consequence_terms=("stop_gained",),
                impact=ImpactSeverity.HIGH,
            ),
        ),
        source_artifact="test",
    )


def _candidate(
    variant_a: VariantRecord,
    variant_b: VariantRecord | None,
    *,
    composite: float,
    gene: str = GENE,
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
        recommended_next_test="Parental segregation testing.",
    )


# ---------------------------------------------------------------------------
# DEFECT 1 — truncation ran before promotion, so the full-match row was discarded
# ---------------------------------------------------------------------------


def _eleven_candidates_with_the_pair_last() -> list[CandidatePair]:
    """Rank 1 is one half of the answer; rank 11 is the answer itself.

    Nine unrelated hypotheses sit between them, so the pair falls one row past the
    ten-row limit. Nothing here is contrived: a compound-het pair is scored partly
    on phase, which is `UNKNOWN` in this pipeline by default (GP-15) and therefore
    penalised, while each constituent single is scored as a dominant/de-novo
    hypothesis with no such penalty. A single outscoring its own parent pair is the
    normal case, not the exotic one — that is the whole premise of ADR 0015.
    """
    half_a = _variant(40_200_000)
    half_b = _variant(40_210_500, ref="G", alt="A")
    filler = [
        _candidate(
            _variant(41_000_000 + 1_000 * index, contig="chr7"),
            _variant(41_500_000 + 1_000 * index, contig="chr7", ref="G", alt="A"),
            composite=0.90 - 0.05 * index,
            gene=f"SYNTHFIL{index}",
        )
        for index in range(9)
    ]
    return [
        _candidate(half_a, None, composite=0.99),
        *filler,
        _candidate(half_a, half_b, composite=0.10),
    ]


def test_a_pair_below_the_row_limit_is_promoted_into_the_submission() -> None:
    """The reproduction: eleven candidates, the answer eleventh.

    Truncating to ten rows *before* promotion deleted the only row that can ever
    equal the answer as a frozenset. What survived was the single at rank 1 — a
    partial match — and the second causal allele appeared nowhere in the file, so
    F-max was capped too. Both metrics, from one line of ordering.
    """
    rows = build_submission_rows(
        _eleven_candidates_with_the_pair_last(), proband_id=ACCEPTED_PROBAND_ID
    )
    assert len(rows) == MAX_SUBMISSION_ROWS

    shapes = [(row.pos_1, int(row.pos_2) if row.pos_2 else None) for row in rows]
    assert (40_200_000, 40_210_500) in shapes, (
        "the pair carrying both causal alleles was truncated away before it could be "
        f"promoted above its own subset: {shapes}"
    )
    assert shapes[0] == (40_200_000, 40_210_500), (
        "the promoted pair must sit immediately above the single it subsumes, which "
        f"was rank 1: {shapes}"
    )


def test_the_promoted_pair_is_worth_100_rank_points_not_50() -> None:
    """The same defect, priced. This is the assertion that matters.

    Before the fix: 50.0 rank points and F-max 0.6667, because the submission held
    one half of the answer and never the whole of it. After: 100.0 and 1.0000.
    """
    text = render_submission_csv_unvalidated(
        build_submission_rows(
            _eleven_candidates_with_the_pair_last(), proband_id=ACCEPTED_PROBAND_ID
        )
    )
    points, fmax = score(text, TRUTH)

    assert points == pytest.approx(100.0), (
        f"the answer scored {points} rank points; a full match at rank 1 is 100.0 and "
        "a partial match at rank 1 is 50.0"
    )
    assert fmax == pytest.approx(1.0), (
        f"F-max was {fmax:.4f}; it is capped below 1.0 whenever the second causal "
        "allele appears in no submitted row"
    )


def test_the_promoted_pair_takes_its_own_subset_slot_and_no_other_row_moves() -> None:
    """ADR 0023's exchange rule: the pair pays for its slot out of the single's.

    Crossing the ten-row cut costs a row, so the row it costs is the subset itself.
    The alternative — slide the pair up and let the last row fall off the end —
    spends a submitted *pair* to keep a single that can never equal a two-variant
    answer as a frozenset. All nine unrelated hypotheses must keep their exact
    positions: a rendering pass does not get to re-rank the science.
    """
    rows = build_submission_rows(
        _eleven_candidates_with_the_pair_last(), proband_id=ACCEPTED_PROBAND_ID
    )
    shapes = [(row.pos_1, int(row.pos_2) if row.pos_2 else None) for row in rows]

    assert (40_200_000, None) not in shapes, (
        "the single was kept at the price of a filler pair. Its variant is already "
        f"re-emitted by the promoted pair, so it bought nothing: {shapes}"
    )
    assert shapes[1:] == [(41_000_000 + 1_000 * i, 41_500_000 + 1_000 * i) for i in range(9)], (
        f"an unrelated row moved; the exchange must touch exactly two rows: {shapes}"
    )


def test_promotion_across_the_cut_is_deterministic_under_input_order() -> None:  # GP-30
    candidates = _eleven_candidates_with_the_pair_last()
    forward = render_submission_csv_unvalidated(
        build_submission_rows(candidates, proband_id=ACCEPTED_PROBAND_ID)
    )
    backward = render_submission_csv_unvalidated(
        build_submission_rows(list(reversed(candidates)), proband_id=ACCEPTED_PROBAND_ID)
    )
    assert forward == backward


# ---------------------------------------------------------------------------
# DEFECT 2 — a pair split across two single-variant rows
# ---------------------------------------------------------------------------


def _split_pair_rows() -> tuple[SubmissionRow, ...]:
    """The two halves of one compound-het hypothesis, as two independent rows.

    Each row is individually legal: right proband, `chr`-prefixed contig, in-range
    EPCR, distinct variant sets, no tie. Nothing in the per-row contract is
    violated, which is exactly why cross-row validation has to catch it.
    """
    return (
        SubmissionRow(
            proband_id=ACCEPTED_PROBAND_ID,
            chrom_1="chr15",
            pos_1=40_200_000,
            ref_1="C",
            alt_1="T",
            chrom_2="",
            pos_2="",
            ref_2="",
            alt_2="",
            epcr=0.9000,
            notes=f"{GENE}; unknown; phase=unknown",
        ),
        SubmissionRow(
            proband_id=ACCEPTED_PROBAND_ID,
            chrom_1="chr15",
            pos_1=40_210_500,
            ref_1="G",
            alt_1="A",
            chrom_2="",
            pos_2="",
            ref_2="",
            alt_2="",
            epcr=0.8000,
            notes=f"{GENE}; unknown; phase=unknown",
        ),
    )


def test_a_split_pair_costs_exactly_half_the_rank_points() -> None:
    """Why the validator has to reject the split, priced with the scorer's arithmetic.

    F-max is a perfect 1.0 — both causal alleles are predicted, and F-max is
    computed per variant, not per row. It is the rank component that collapses:
    ``row.variants == true_variants`` is never true for a single-variant row, so
    the best available outcome is partial credit. **100 -> 50 for nothing.**
    """
    text = render_submission_csv_unvalidated(_split_pair_rows())
    points, fmax = score(text, TRUTH)

    assert points == pytest.approx(50.0)
    assert fmax == pytest.approx(1.0)

    joined = render_submission_csv_unvalidated(
        (
            SubmissionRow(
                proband_id=ACCEPTED_PROBAND_ID,
                chrom_1="chr15",
                pos_1=40_200_000,
                ref_1="C",
                alt_1="T",
                chrom_2="chr15",
                pos_2="40210500",
                ref_2="G",
                alt_2="A",
                epcr=0.9000,
                notes=f"{GENE}; compound_heterozygous; phase=unknown",
            ),
        )
    )
    assert score(joined, TRUTH) == (pytest.approx(100.0), pytest.approx(1.0)), (
        "the identical prediction, emitted as ONE row with the _2 columns, is worth twice as much"
    )


def test_validate_submission_rejects_a_pair_split_across_two_rows() -> None:
    """The reproduction: every per-row rule passes and the submission is still wrong.

    Right proband, `chr`-prefixed contigs, in-range and untied EPCRs, distinct
    variant sets — the duplicate-set check and the tie check both see nothing,
    because neither is looking across rows for a hypothesis that was *taken apart*.
    """
    ok, errors = validate_submission(render_submission_csv_unvalidated(_split_pair_rows()))

    assert not ok, "a pair split across two rows passed the contract self-check"
    assert any("split" in error.lower() for error in errors), errors


def test_write_submission_refuses_a_split_pair(tmp_path: Path) -> None:
    """The gate that matters: the split must not be able to reach a file."""
    target = tmp_path / "track1_submission.csv"
    with pytest.raises(ValueError, match="split"):
        write_submission(_split_pair_rows(), target)
    assert not target.exists(), "the rejected submission was written anyway"


def test_the_same_two_variants_joined_into_one_row_validate_cleanly() -> None:
    """The check must accept the fix it is asking for, or it is just noise."""
    joined = (
        SubmissionRow(
            proband_id=ACCEPTED_PROBAND_ID,
            chrom_1="chr15",
            pos_1=40_200_000,
            ref_1="C",
            alt_1="T",
            chrom_2="chr15",
            pos_2="40210500",
            ref_2="G",
            alt_2="A",
            epcr=0.9000,
            notes=f"{GENE}; compound_heterozygous; phase=unknown",
        ),
    )
    ok, errors = validate_submission(render_submission_csv_unvalidated(joined))
    assert ok, errors


def test_a_pair_row_plus_its_own_halves_is_not_a_split() -> None:
    """ADR 0015's own output must not trip the new check.

    A pair row with both of its single-variant rows beneath it proposes the joined
    hypothesis explicitly. Nothing has been taken apart, so there is nothing to
    report — and a check that fired here would fail the golden case.
    """
    rows = (
        SubmissionRow(
            proband_id=ACCEPTED_PROBAND_ID,
            chrom_1="chr15",
            pos_1=40_200_000,
            ref_1="C",
            alt_1="T",
            chrom_2="chr15",
            pos_2="40210500",
            ref_2="G",
            alt_2="A",
            epcr=0.9000,
            notes=f"{GENE}; compound_heterozygous; phase=unknown",
        ),
        *_split_pair_rows(),
    )
    rows = (rows[0], replace(rows[1], epcr=0.8000), replace(rows[2], epcr=0.7000))
    ok, errors = validate_submission(render_submission_csv_unvalidated(rows))
    assert ok, errors


def test_two_singles_in_different_genes_are_not_a_split() -> None:
    """Two unrelated dominant hypotheses are two hypotheses, not a dismembered pair.

    The check keys on the gene the pipeline attributed the row to. Requiring every
    pair of single-variant rows to be joined would demand 45 pair rows for ten
    singles, which is neither possible nor meaningful.
    """
    left, right = _split_pair_rows()
    rows = (left, replace(right, notes="SYNTHMET2; unknown; phase=unknown"))
    ok, errors = validate_submission(render_submission_csv_unvalidated(rows))
    assert ok, errors


# ---------------------------------------------------------------------------
# DEFECT 2 (the other half) — the exported renderer was a validation bypass
# ---------------------------------------------------------------------------


def _valid_row(position: int, *, epcr: float, proband: str = ACCEPTED_PROBAND_ID) -> SubmissionRow:
    return SubmissionRow(
        proband_id=proband,
        chrom_1="chr15",
        pos_1=position,
        ref_1="C",
        alt_1="T",
        chrom_2="chr7",
        pos_2=str(55_000_000 + position),
        ref_2="A",
        alt_2="G",
        epcr=epcr,
        notes=f"SYNTHFIL{position}; compound_heterozygous; phase=unknown",
    )


@pytest.mark.parametrize(
    ("label", "rows", "fragment"),
    [
        (
            "eleven rows",
            tuple(_valid_row(40_200_000 + i, epcr=0.99 - 0.02 * i) for i in range(11)),
            "row limit",
        ),
        (
            "bare contig",
            (replace(_valid_row(40_200_000, epcr=0.9), chrom_1="15"),),
            "UCSC 'chr' prefix",
        ),
        (
            "duplicate epcr",
            (_valid_row(40_200_000, epcr=0.9), _valid_row(40_300_000, epcr=0.9)),
            "repeats the value",
        ),
        (
            "arbitrary proband id",
            (_valid_row(40_200_000, epcr=0.9, proband="CHILD-7"),),
            "Unknown proband_id",
        ),
        (
            "epcr out of range",
            (replace(_valid_row(40_200_000, epcr=0.9), epcr=0.0),),
            "0 < epcr <= 1",
        ),
    ],
    ids=lambda value: value if isinstance(value, str) else "",
)
def test_the_renderer_refuses_every_state_the_validator_rejects(
    label: str, rows: tuple[SubmissionRow, ...], fragment: str
) -> None:
    """`render_submission_csv` used to serialise anything it was handed.

    `SubmissionRow` carries no invariants of its own, so the raw serialiser was a
    one-call public path from an arbitrary dataclass to submission bytes, around
    every rule in the contract. Each state below is silent at upload time and fatal
    at scoring time; the reviewer confirmed `write_submission` already caught them
    all, which is exactly the point — the checks existed and the renderer stood
    beside them.
    """
    with pytest.raises(ValueError, match=fragment):
        render_submission_csv(rows)

    # The escape hatch still renders — it is named for what it skips, not removed —
    # and `validate_submission` still names the same failure on those bytes.
    ok, errors = validate_submission(render_submission_csv_unvalidated(rows))
    assert not ok, label
    assert any(fragment in error for error in errors), (label, errors)


def test_the_unvalidated_renderer_is_not_reachable_from_the_package_namespace() -> None:
    """The escape hatch is module-level only, and named for what it does not do.

    `mva.reporting` is the import surface the rest of the pipeline uses. Keeping
    the raw path out of it means the shortest spelling of "serialise these rows" is
    the one that checks them, while a test that genuinely needs invalid bytes has
    to say `unvalidated` out loud.
    """
    from mva import reporting
    from mva.reporting import track1

    assert not hasattr(reporting, "render_submission_csv_unvalidated")
    assert "render_submission_csv_unvalidated" not in reporting.__all__
    assert "render_submission_csv_unvalidated" in track1.__all__
