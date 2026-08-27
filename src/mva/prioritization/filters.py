"""Hard filters and soft flags for candidate variants (GP-13).

The single rule this module exists to enforce: **filtering and ranking are
different operations**. A hard filter may only remove a record that is *invalid
or genuinely impossible* — the wrong assembly, a contig this pipeline is not
entitled to reason about, a placeholder allele, or a genotype that carries no
alternate allele at all. Nothing is ever removed for being common, low quality
or predicted-benign; those are soft flags that travel with the record and
down-rank it, so a human reviewer can still see the candidate and disagree.

The asymmetry is deliberate. A deleted candidate is invisible and therefore
unfalsifiable; a flagged candidate ranked 47th is merely unlikely. Published
pathogenic alleles that exceed naive frequency expectations in
under-represented populations are exactly the case a frequency *filter* loses
silently, which is why this module has none.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from mva.config import FrequencyThresholds, QualityThresholds
from mva.models.genome import CANONICAL_CONTIGS, GenomeBuild
from mva.models.variant import (
    INGESTION_QC_FLAGS,
    UNTRUSTED_CALL_FLAGS,
    FilterStatus,
    ImpactSeverity,
    VariantRecord,
    Zygosity,
)

# ---------------------------------------------------------------------------
# Reason codes for HARD removal. Every one of these describes a record that is
# invalid or carries no alternate allele — never one that is merely unlikely.
# ---------------------------------------------------------------------------

REASON_WRONG_GENOME_BUILD = "wrong_genome_build"
"""Coordinates from another assembly. Positions differ by megabases; comparing
them would be confidently wrong rather than merely uncertain (GP-11)."""

REASON_NON_CANONICAL_CONTIG = "non_canonical_contig"
"""Alt/decoy/unplaced scaffold, or a contig outside chr1-22,X,Y,M."""

REASON_UNINFORMATIVE_ALT = "uninformative_alt_allele"
"""'*' (spanning deletion) or '.' (missing) — not an independently rankable allele."""

REASON_NO_ALT_ALLELE = "no_alt_allele_called"
"""Homozygous reference: the sample carries no copy of the alternate allele."""

REASON_MISSING_GENOTYPE = "missing_genotype"
"""No call at this site. Retained nowhere, but also never coerced to hom-ref."""

#: Ordered so that ``counts`` keys and the removal log are stable across runs.
HARD_FILTER_REASONS: tuple[str, ...] = (
    REASON_WRONG_GENOME_BUILD,
    REASON_NON_CANONICAL_CONTIG,
    REASON_UNINFORMATIVE_ALT,
    REASON_NO_ALT_ALLELE,
    REASON_MISSING_GENOTYPE,
)

# ---------------------------------------------------------------------------
# Soft flags. These down-rank. They never delete.
# ---------------------------------------------------------------------------

FLAG_COMMON_VARIANT = "common_variant"
FLAG_LOW_FREQUENCY_VARIANT = "low_frequency_variant"
FLAG_NO_FREQUENCY_DATA = "no_frequency_data"
FLAG_BENIGN_CONSEQUENCE = "benign_consequence"
FLAG_LOW_QUALITY_CALL = "low_quality_call"
FLAG_HOMOZYGOUS_CALL = "homozygous_call"
FLAG_POSSIBLE_MOSAIC = "possible_mosaic"
FLAG_PLAUSIBLE_CANDIDATE = "plausible_candidate"

#: Applied in this order so that a record's ``qc_flags`` tuple is deterministic.
SOFT_FLAGS: tuple[str, ...] = (
    FLAG_COMMON_VARIANT,
    FLAG_LOW_FREQUENCY_VARIANT,
    FLAG_NO_FREQUENCY_DATA,
    FLAG_BENIGN_CONSEQUENCE,
    FLAG_LOW_QUALITY_CALL,
    FLAG_POSSIBLE_MOSAIC,
    FLAG_HOMOZYGOUS_CALL,
)

#: QC flags an upstream ingestion stage may attach that mean "do not trust this
#: call". Recognised by name so that `low_quality_call` genuinely *inherits*
#: upstream findings rather than silently re-deriving a subset of them.
#:
#: DERIVED, not re-typed. The hand-written version of this set named eight flags
#: nothing in `src/` ever emitted (`strand_bias`, `low_mappability`, ...) while
#: missing three that ingestion really does raise — so a `ref_allele_mismatch`
#: variant, whose REF disagrees with the reference sequence, reached scoring with
#: an analytical validity of 1.0 and no contradiction at all: a confidently wrong
#: coordinate ranking like a clean call. Two packages that may not import each
#: other (GP-03) share the vocabulary through :mod:`mva.models.variant` instead.
INGESTION_QUALITY_FLAGS: frozenset[str] = UNTRUSTED_CALL_FLAGS | {FLAG_LOW_QUALITY_CALL}

#: Ingestion findings that are recognised but are NOT quality concerns. Named
#: rather than omitted, so "considered and deliberately not treated as a defect"
#: stays distinguishable from "forgotten". Together with
#: :data:`INGESTION_QUALITY_FLAGS` this partitions :data:`INGESTION_QC_FLAGS`,
#: and a test holds the two sets to that.
INGESTION_NEUTRAL_FLAGS: frozenset[str] = INGESTION_QC_FLAGS - INGESTION_QUALITY_FLAGS

#: Impacts that, on their own, do not support a severe loss-of-function hypothesis.
BENIGN_IMPACTS: frozenset[ImpactSeverity] = frozenset({ImpactSeverity.LOW, ImpactSeverity.MODIFIER})

#: VCF placeholder alleles.
UNINFORMATIVE_ALT_ALLELES: frozenset[str] = frozenset({"*", "."})


@dataclass(frozen=True)
class FilterResult:
    """Outcome of the hard-filter pass.

    ``removed`` lists ONLY invalid records, with the reason code that removed
    them, so the run report can state exactly what was discarded and why.
    ``flagged`` is ``retained`` carrying the soft flags that need no
    configuration (zygosity, absence of frequency data, inherited upstream
    quality findings); the threshold-dependent flags are added by
    :func:`apply_soft_flags`, which is idempotent with respect to these.
    """

    retained: tuple[VariantRecord, ...]
    removed: tuple[tuple[str, str], ...]
    flagged: tuple[VariantRecord, ...]
    counts: dict[str, int]

    @property
    def removed_ids(self) -> tuple[str, ...]:
        return tuple(variant_id for variant_id, _ in self.removed)


def _hard_filter_reason(variant: VariantRecord, expected_build: GenomeBuild) -> str | None:
    """Return the reason this record is invalid, or ``None`` if it is usable.

    Checks run in a fixed order so the recorded reason is deterministic when a
    record fails more than one.
    """
    if variant.build is not expected_build:
        return REASON_WRONG_GENOME_BUILD
    if variant.coordinate.contig not in CANONICAL_CONTIGS:
        return REASON_NON_CANONICAL_CONTIG
    if variant.coordinate.alt in UNINFORMATIVE_ALT_ALLELES:
        return REASON_UNINFORMATIVE_ALT
    if variant.genotype.zygosity is Zygosity.HOM_REF:
        return REASON_NO_ALT_ALLELE
    if variant.genotype.zygosity is Zygosity.UNKNOWN:
        return REASON_MISSING_GENOTYPE
    return None


def _configuration_free_flags(variant: VariantRecord) -> tuple[str, ...]:
    """Soft flags derivable without any threshold configuration."""
    flags: list[str] = []
    if not variant.has_frequency_data:
        flags.append(FLAG_NO_FREQUENCY_DATA)
    if variant.filter_status is FilterStatus.FILTERED or (
        set(variant.qc_flags) & INGESTION_QUALITY_FLAGS
    ):
        flags.append(FLAG_LOW_QUALITY_CALL)
    if variant.genotype.zygosity is Zygosity.HOM_ALT:
        flags.append(FLAG_HOMOZYGOUS_CALL)
    return tuple(flags)


def apply_hard_filters(
    variants: Sequence[VariantRecord], *, expected_build: GenomeBuild
) -> FilterResult:
    """Remove only invalid or impossible records (GP-13).

    Removal is restricted to: wrong genome build, non-canonical contig, a ``*``
    or ``.`` alternate allele, a homozygous-reference genotype (no alternate
    allele is present in the sample) and a missing genotype. Commonness, poor
    call quality and benign predicted consequence are NOT removal reasons here
    and never will be — see :func:`apply_soft_flags`.
    """
    retained: list[VariantRecord] = []
    removed: list[tuple[str, str]] = []
    reason_counts: dict[str, int] = dict.fromkeys(HARD_FILTER_REASONS, 0)

    for variant in variants:
        reason = _hard_filter_reason(variant, expected_build)
        if reason is None:
            retained.append(variant)
        else:
            removed.append((variant.variant_id, reason))
            reason_counts[reason] += 1

    flagged = tuple(
        variant.with_qc_flags(*_configuration_free_flags(variant)) for variant in retained
    )

    counts: dict[str, int] = {
        "input": len(variants),
        "retained": len(retained),
        "removed": len(removed),
    }
    for reason in HARD_FILTER_REASONS:
        counts[f"removed_{reason}"] = reason_counts[reason]

    return FilterResult(
        retained=tuple(retained), removed=tuple(removed), flagged=flagged, counts=counts
    )


def _frequency_flag(variant: VariantRecord, frequency: FrequencyThresholds) -> str | None:
    """Frequency band flag, using the MAXIMUM AF across populations.

    Maximum rather than global: a variant common in any single ancestry is not a
    plausible ultra-rare cause, and the global figure systematically
    under-estimates frequency for alleles enriched in under-represented cohorts.

    Populations reporting fewer than ``frequency.min_allele_number`` alleles are
    skipped, because a single observation in a 40-chromosome subpopulation is an
    AF of 0.025 and would flag a genuine founder allele ``common_variant``
    (ADR 0010). Where no population is large enough, the result is treated as an
    absence of frequency data rather than as either rarity or commonness.
    """
    observed = variant.max_allele_frequency(min_allele_number=frequency.min_allele_number)
    if observed is None:
        return FLAG_NO_FREQUENCY_DATA
    if observed.allele_frequency > frequency.max_plausible_recessive:
        return FLAG_COMMON_VARIANT
    if observed.allele_frequency > frequency.low_frequency:
        return FLAG_LOW_FREQUENCY_VARIANT
    return None


def _quality_flags(variant: VariantRecord, quality: QualityThresholds) -> tuple[str, ...]:
    """Call-quality flags. Breaching a threshold flags; it never deletes."""
    genotype = variant.genotype
    low_quality = variant.filter_status is FilterStatus.FILTERED or bool(
        set(variant.qc_flags) & INGESTION_QUALITY_FLAGS
    )
    if genotype.depth is not None and genotype.depth < quality.min_depth:
        low_quality = True
    if (
        genotype.genotype_quality is not None
        and genotype.genotype_quality < quality.min_genotype_quality
    ):
        low_quality = True

    mosaic = False
    # `VariantRecord.allele_fraction`, never `Genotype.allele_balance`: on a
    # record split out of a multiallelic site the latter is alt/(ref+alt) over
    # the SITE's ref depth, which reads 0.91 for a textbook compound het and
    # flags it low quality for the shape of its VCF line (ADR 0010 sibling fix).
    fraction = variant.allele_fraction
    if genotype.is_heterozygous and fraction is not None:
        if fraction > quality.max_allele_balance_het:
            low_quality = True
        elif fraction < quality.min_allele_balance_het:
            # Below the het band but above the mosaic floor is flagged as
            # possible mosaicism rather than as noise: in a mosaic aneuploidy
            # disorder a skewed allele balance may be the signal itself.
            if fraction >= quality.mosaic_allele_balance_floor:
                mosaic = True
            else:
                low_quality = True

    flags: list[str] = []
    if low_quality:
        flags.append(FLAG_LOW_QUALITY_CALL)
    if mosaic:
        flags.append(FLAG_POSSIBLE_MOSAIC)
    return tuple(flags)


def _consequence_flag(variant: VariantRecord) -> str | None:
    """Flag records whose *only* predicted impacts are MODIFIER or LOW."""
    if not variant.consequences:
        return None
    if all(csq.impact in BENIGN_IMPACTS for csq in variant.consequences):
        return FLAG_BENIGN_CONSEQUENCE
    return None


def apply_soft_flags(
    variants: Sequence[VariantRecord],
    *,
    frequency: FrequencyThresholds,
    quality: QualityThresholds,
) -> tuple[VariantRecord, ...]:
    """Annotate records with down-ranking markers. Nothing is removed (GP-13).

    Flags applied: ``common_variant``, ``low_frequency_variant``,
    ``no_frequency_data``, ``benign_consequence``, ``low_quality_call``,
    ``possible_mosaic`` and ``homozygous_call``. Flags are additive and
    order-stable, so re-running this pass is a no-op.
    """
    flagged: list[VariantRecord] = []
    for variant in variants:
        found: set[str] = set()
        frequency_flag = _frequency_flag(variant, frequency)
        if frequency_flag is not None:
            found.add(frequency_flag)
        found.update(_quality_flags(variant, quality))
        consequence_flag = _consequence_flag(variant)
        if consequence_flag is not None:
            found.add(consequence_flag)
        if variant.genotype.zygosity is Zygosity.HOM_ALT:
            found.add(FLAG_HOMOZYGOUS_CALL)
        ordered = tuple(flag for flag in SOFT_FLAGS if flag in found)
        flagged.append(variant.with_qc_flags(*ordered))
    return tuple(flagged)


def select_candidate_variants(
    variants: Sequence[VariantRecord], *, frequency: FrequencyThresholds
) -> tuple[VariantRecord, ...]:
    """Mark the subset that is *plausible* for a severe recessive hypothesis.

    This is a selection, not a filter: the caller keeps the full flagged set and
    is expected to carry it forward (GP-13). The returned records are the same
    objects with an added ``plausible_candidate`` flag, so a downstream stage
    can prefer them without losing the ability to see everything else.

    Plausibility is deliberately narrow — it uses only the two things a variant
    must have to be a candidate at all: an alternate allele in this sample, and
    a frequency that does not already exceed the maximum plausible for a severe
    recessive disorder. Predicted consequence is intentionally NOT used: impact
    prediction is the weakest link in the chain, and excluding LOW-impact calls
    here would silently discard splice-region synonymous variants.
    """
    selected: list[VariantRecord] = []
    for variant in variants:
        if not variant.genotype.carries_alt:
            continue
        observed = variant.max_allele_frequency(min_allele_number=frequency.min_allele_number)
        if observed is not None and observed.allele_frequency > frequency.max_plausible_recessive:
            continue
        selected.append(variant.with_qc_flags(FLAG_PLAUSIBLE_CANDIDATE))
    return tuple(selected)
