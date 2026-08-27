"""The per-gene cap must never delete the answer for being on the right.

`generate_pairs` used to build every within-gene pair, sort by
``(gene, first coordinate, second coordinate, pair_id)`` and cut at
``max_pairs_per_gene``. The kept set was therefore "every pair involving the
leftmost variant, then the second, ..." — so once the cap bit, whole variants at
the right-hand end of a gene became unreachable in *any* pairing. At ten variants
in a gene 44% of pairings survived; at twenty, 10.5%. Which survived was decided
by chromosomal position, which is not evidence of anything.

Two properties are locked here:

1. a gene with many variants does not lose its most plausible pairing to a
   variant's coordinate, and
2. when a cap does fire it is **surfaced** — flagged on every surviving candidate
   *and* reported as an event a run report can print. A cap that silently deletes
   the answer is the worst failure mode in this project; a cap that says so is a
   known limitation.
"""

from __future__ import annotations

import pytest

from mva.config import FrequencyThresholds
from mva.models.genome import GenomeBuild, GenomicCoordinate
from mva.models.variant import (
    ConsequenceAnnotation,
    FilterStatus,
    Genotype,
    ImpactSeverity,
    PopulationFrequency,
    VariantRecord,
    Zygosity,
)
from mva.prioritization.filters import FLAG_LOW_QUALITY_CALL
from mva.prioritization.pairing import (
    DEFAULT_MAX_PAIRING_VARIANTS,
    FLAG_CAP_TRUNCATED,
    generate_pair_candidates,
    generate_pairs,
)

pytestmark = pytest.mark.unit

GENE = "SYNTHKIN1"
FREQUENCY = FrequencyThresholds()


def make_variant(
    position: int,
    *,
    impact: ImpactSeverity = ImpactSeverity.HIGH,
    allele_frequency: float = 0.0,
    zygosity: Zygosity = Zygosity.HET,
    qc_flags: tuple[str, ...] = (),
    gene: str = GENE,
) -> VariantRecord:
    """One alt-carrying call in ``gene``. Defaults describe a clean, rare, HIGH one."""
    terms = ("stop_gained",) if impact is ImpactSeverity.HIGH else ("synonymous_variant",)
    return VariantRecord(
        coordinate=GenomicCoordinate(
            build=GenomeBuild.GRCH38, contig="chr15", position=position, ref="C", alt="T"
        ),
        genotype=Genotype(
            zygosity=zygosity,
            genotype_string="1/1" if zygosity is Zygosity.HOM_ALT else "0/1",
            depth=45,
            ref_reads=0 if zygosity is Zygosity.HOM_ALT else 23,
            alt_reads=45 if zygosity is Zygosity.HOM_ALT else 22,
            genotype_quality=99,
        ),
        filter_status=FilterStatus.PASS,
        consequences=(
            ConsequenceAnnotation(
                gene_symbol=gene,
                transcript_id="SYNTHT0001.1",
                is_canonical=True,
                consequence_terms=terms,
                impact=impact,
            ),
        ),
        population_frequencies=(
            PopulationFrequency(
                source="SYNTH_gnomAD",
                version="v0.0-synthetic",
                population="global",
                allele_frequency=allele_frequency,
                allele_number=152_312,
            ),
        ),
        qc_flags=qc_flags,
        source_artifact="test",
    )


def buried_answer_corpus() -> tuple[list[VariantRecord], tuple[str, str]]:
    """Twenty calls in one gene; the only plausible pair sits at the right-hand end.

    Eighteen common, LOW-impact calls occupy the left of the gene. The two
    ultra-rare HIGH-impact calls that a clinician would actually pair are the
    last two coordinates. Under coordinate truncation at a cap of 20 they are
    unreachable: every kept candidate is anchored on one of the two leftmost
    variants.
    """
    common = [
        make_variant(40_200_000 + index * 1_000, impact=ImpactSeverity.LOW, allele_frequency=0.12)
        for index in range(18)
    ]
    first = make_variant(40_300_000)
    second = make_variant(40_301_000)
    return [*common, first, second], (first.variant_id, second.variant_id)


# ---------------------------------------------------------------------------
# 1. The defect: a pairing lost to genomic position.
# ---------------------------------------------------------------------------


def test_the_plausible_pairing_survives_a_cap_that_used_to_delete_it() -> None:
    variants, answer = buried_answer_corpus()
    result = generate_pair_candidates(variants, max_pairs_per_gene=20, frequency=FREQUENCY)

    surviving = {candidate.variant_ids for candidate in result.candidates}
    assert answer in surviving, (
        "the two ultra-rare HIGH-impact calls were dropped by the per-gene cap. "
        "They are the last two coordinates in the gene, which is the only thing "
        "wrong with them — truncation must order by plausibility, not position."
    )


def test_truncation_keeps_plausible_candidates_not_leftmost_ones() -> None:
    """The kept set must not be a prefix of the gene's coordinate order."""
    variants, _ = buried_answer_corpus()
    result = generate_pair_candidates(variants, max_pairs_per_gene=20, frequency=FREQUENCY)

    anchors = {candidate.variant_a.coordinate.position for candidate in result.candidates}
    assert 40_300_000 in anchors and 40_301_000 in anchors
    assert len(anchors) > 2, (
        f"every kept candidate is anchored on one of {sorted(anchors)} — the cap is "
        "still slicing the coordinate-sorted list"
    )


def test_a_low_quality_call_is_capped_before_a_clean_one() -> None:
    """Call quality outranks position too: analytical validity comes first."""
    clean = [make_variant(40_400_000 + index * 1_000) for index in range(3)]
    suspect = [
        make_variant(40_200_000 + index * 1_000, qc_flags=(FLAG_LOW_QUALITY_CALL,))
        for index in range(3)
    ]
    result = generate_pair_candidates([*suspect, *clean], max_pairs_per_gene=3, frequency=FREQUENCY)
    kept_positions = {
        variant.coordinate.position
        for candidate in result.candidates
        for variant in candidate.variants
    }
    assert kept_positions <= {40_400_000, 40_401_000, 40_402_000}, (
        f"a low-quality call was kept over a clean one: {sorted(kept_positions)}"
    )


# ---------------------------------------------------------------------------
# 2. When a cap fires, it says so.
# ---------------------------------------------------------------------------


def test_the_cap_flag_fires_and_reaches_a_reportable_event() -> None:
    variants, _ = buried_answer_corpus()
    result = generate_pair_candidates(variants, max_pairs_per_gene=20, frequency=FREQUENCY)

    assert result.truncated
    assert result.truncated_genes == (GENE,)
    assert all(FLAG_CAP_TRUNCATED in c.flags for c in result.candidates), (
        "a candidate from a truncated gene must carry the flag: the hypothesis "
        "list it belongs to is known to be incomplete"
    )

    (event,) = result.cap_events
    assert event.gene_symbol == GENE
    assert event.variants == 20
    assert event.kept == 20
    assert event.generated > event.kept
    assert event.dropped == event.generated - event.kept
    assert result.dropped_candidates == event.dropped


def test_the_cap_surfaces_as_a_warning_naming_the_flag_and_the_counts() -> None:
    """The signal has to reach a human, not just a field on a record.

    The pre-fix code set the flag on every surviving candidate and stopped there,
    so a run could delete the answer and produce output indistinguishable from a
    complete one.
    """
    variants, _ = buried_answer_corpus()
    result = generate_pair_candidates(variants, max_pairs_per_gene=20, frequency=FREQUENCY)

    (warning,) = result.warnings
    assert FLAG_CAP_TRUNCATED in warning
    assert GENE in warning
    assert str(result.dropped_candidates) in warning
    # GP-41: counts and gene symbols, never coordinates or variant identifiers.
    assert "chr15" not in warning
    assert "40200000" not in warning


def test_no_cap_means_no_flag_no_event_and_no_warning() -> None:
    variants = [make_variant(40_200_000 + index * 1_000) for index in range(4)]
    result = generate_pair_candidates(variants, max_pairs_per_gene=20, frequency=FREQUENCY)

    assert not result.truncated
    assert result.cap_events == ()
    assert result.warnings == ()
    assert not any(FLAG_CAP_TRUNCATED in c.flags for c in result.candidates)


# ---------------------------------------------------------------------------
# 3. The input bound removes pairings, never variants.
# ---------------------------------------------------------------------------


def test_a_variant_excluded_from_pairing_still_forms_its_own_hypothesis() -> None:
    """Bounding the pairing input must not make a variant wholly unreachable.

    A homozygous call explains both gene copies by itself; it needs no partner. If
    the input bound could delete it outright, the bound would be a hard filter on
    the hypothesis space, which is exactly what ADR 0005 and GP-13 forbid.
    """
    plausible = [make_variant(40_200_000 + index * 1_000) for index in range(3)]
    excluded = make_variant(
        40_900_000,
        impact=ImpactSeverity.LOW,
        allele_frequency=0.30,
        zygosity=Zygosity.HOM_ALT,
    )
    result = generate_pair_candidates(
        [*plausible, excluded], max_pairs_per_gene=20, max_pairing_variants=3, frequency=FREQUENCY
    )

    singles = {c.variant_ids for c in result.candidates if not c.is_pair}
    pairs = {c.variant_ids for c in result.candidates if c.is_pair}
    assert (excluded.variant_id,) in singles, (
        "the variant bounded out of PAIR enumeration also lost its single-variant "
        "hypothesis; the bound deleted a candidate rather than a combination"
    )
    assert not any(excluded.variant_id in ids for ids in pairs)
    assert result.truncated, "an input bound that fires must be reported like a cap"


def test_the_input_bound_prefers_plausible_variants() -> None:
    common = [
        make_variant(40_200_000 + index * 1_000, impact=ImpactSeverity.LOW, allele_frequency=0.20)
        for index in range(5)
    ]
    rare = [make_variant(40_800_000 + index * 1_000) for index in range(3)]
    result = generate_pair_candidates(
        [*common, *rare], max_pairs_per_gene=50, max_pairing_variants=3, frequency=FREQUENCY
    )
    paired = {
        variant.variant_id
        for candidate in result.candidates
        if candidate.is_pair
        for variant in candidate.variants
    }
    assert paired == {variant.variant_id for variant in rare}


# ---------------------------------------------------------------------------
# 4. Nothing moves when no cap fires, and repeat runs are identical (GP-30).
# ---------------------------------------------------------------------------


def test_output_is_unchanged_when_no_bound_or_cap_applies() -> None:
    """The uncapped path must be byte-for-byte what it was, with or without thresholds."""
    variants = [make_variant(40_200_000 + index * 1_000) for index in range(5)]
    assert len(variants) < DEFAULT_MAX_PAIRING_VARIANTS

    without = generate_pairs(variants)
    with_thresholds = generate_pairs(variants, frequency=FREQUENCY)
    assert [c.pair_id for c in without] == [c.pair_id for c in with_thresholds]
    assert [c.flags for c in without] == [c.flags for c in with_thresholds]
    # Canonical order: gene, first coordinate, second coordinate, pair id.
    assert list(without) == sorted(without, key=lambda c: c.sort_key())


def test_repeat_calls_are_identical_regardless_of_input_order() -> None:
    variants, _ = buried_answer_corpus()
    forward = generate_pair_candidates(variants, max_pairs_per_gene=20, frequency=FREQUENCY)
    backward = generate_pair_candidates(
        list(reversed(variants)), max_pairs_per_gene=20, frequency=FREQUENCY
    )
    assert [c.pair_id for c in forward.candidates] == [c.pair_id for c in backward.candidates]
    assert forward.cap_events == backward.cap_events


def test_a_cap_below_one_is_refused() -> None:
    with pytest.raises(ValueError, match="max_pairs_per_gene"):
        generate_pair_candidates([make_variant(40_200_000)], max_pairs_per_gene=0)
    with pytest.raises(ValueError, match="max_pairing_variants"):
        generate_pair_candidates([make_variant(40_200_000)], max_pairing_variants=1)
