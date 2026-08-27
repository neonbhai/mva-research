"""Analytical quality control: flag calls, never delete them.

GP-13 in its most literal form. Every threshold breach here produces a *flag* on
the record and an ``EvidenceItem`` explaining it; nothing is removed. Ranking
happens later and is free to push a flagged call to the bottom, but a human
reviewing the run can always see what was down-ranked and why.

The one judgement call worth stating explicitly is ``possible_mosaic``. A
heterozygous call whose allele balance sits below the expected band is, in a
typical germline pipeline, an artifact. In **Mosaic Variegated Aneuploidy** it may
be the finding itself: a somatic variant present in a fraction of cells produces
exactly that signature. So an allele balance below the het band but at or above
``mosaic_allele_balance_floor`` is flagged ``possible_mosaic`` *instead of*
``low_allele_balance``, and its evidence item is NEUTRAL rather than
contradicting. Calling it noise would be assuming the answer.

Every item carries a non-empty ``limitations`` string (GP-17), a
``DIRECT_MEASUREMENT`` type at ``OBSERVED_DATA`` tier (GP-12 — these really are
measurements, not predictions), and a timestamp from the injected clock (GP-30).
"""

from __future__ import annotations

import statistics
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime
from typing import Final

from mva.clock import Clock
from mva.config import QualityThresholds
from mva.ingestion.normalise import OP_SPLIT_MULTIALLELIC, REF_ALLELE_MISMATCH_FLAG
from mva.models.base import AssertionTier
from mva.models.evidence import (
    EvidenceCategory,
    EvidenceDirection,
    EvidenceItem,
    EvidenceStrength,
    EvidenceType,
    make_evidence_id,
)
from mva.models.variant import FilterStatus, VariantRecord, Zygosity

# ---------------------------------------------------------------------------
# Flag vocabulary
# ---------------------------------------------------------------------------

FLAG_FILTERED_BY_CALLER: Final = "filtered_by_caller"
FLAG_LOW_DEPTH: Final = "low_depth"
FLAG_LOW_GQ: Final = "low_gq"
FLAG_POSSIBLE_MOSAIC: Final = "possible_mosaic"
FLAG_LOW_ALLELE_BALANCE: Final = "low_allele_balance"
FLAG_HIGH_ALLELE_BALANCE: Final = "high_allele_balance"
FLAG_NO_QUALITY_METRICS: Final = "no_quality_metrics"

#: Canonical emission order, so ``qc_flags`` is deterministic (GP-30).
FLAG_ORDER: Final[tuple[str, ...]] = (
    FLAG_FILTERED_BY_CALLER,
    FLAG_LOW_DEPTH,
    FLAG_LOW_GQ,
    FLAG_POSSIBLE_MOSAIC,
    FLAG_LOW_ALLELE_BALANCE,
    FLAG_HIGH_ALLELE_BALANCE,
    FLAG_NO_QUALITY_METRICS,
)

#: Flags that argue *against* the call being analytically real.
CONTRADICTING_FLAGS: Final[frozenset[str]] = frozenset(
    {
        REF_ALLELE_MISMATCH_FLAG,
        FLAG_FILTERED_BY_CALLER,
        FLAG_LOW_DEPTH,
        FLAG_LOW_GQ,
        FLAG_LOW_ALLELE_BALANCE,
        FLAG_HIGH_ALLELE_BALANCE,
    }
)

TOOL_NAME: Final = "mva.ingestion.qc"

_BASE_LIMITATIONS: Final = (
    "Derived only from the caller's own DP, GQ, AD and FILTER for a single sample. "
    "It does not establish that the variant is real: there is no orthogonal assay, no "
    "re-alignment, and no assessment of mappability, repeat context or caller-specific "
    "systematic error at this locus. The thresholds are configuration defaults, not "
    "clinically validated cut-points."
)
_MOSAIC_LIMITATION: Final = (
    " Allele balance below the germline het band is consistent with somatic mosaicism "
    "but is equally consistent with allelic dropout or a mapping artifact; this "
    "evidence cannot distinguish them and asserts no mosaic fraction."
)
_REF_MISMATCH_LIMITATION: Final = (
    " The REF allele disagrees with the supplied reference sequence, which usually "
    "indicts the assembly or the reference FASTA rather than the genotype call; this "
    "evidence cannot tell those cases apart."
)
_NO_METRICS_LIMITATION: Final = (
    " No depth, genotype-quality or allelic-depth values were reported, so the absence "
    "of a quality flag here is absence of information, not evidence of quality."
)


# ---------------------------------------------------------------------------
# Public result type
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class QcResult:
    """The same variants, now flagged, with one evidence item each.

    ``metrics`` holds aggregate counts and summary statistics only. It is
    deliberately safe to write to a ``DERIVED_SAFE`` artifact: no coordinate, no
    genotype and no sample identifier appears in it.
    """

    variants: tuple[VariantRecord, ...]
    evidence: tuple[EvidenceItem, ...]
    metrics: dict[str, int | float]


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def assess_quality(
    variants: Sequence[VariantRecord],
    *,
    thresholds: QualityThresholds,
    clock: Clock,
    tool_version: str = "mva-qc/0.1.0",
) -> QcResult:
    """Flag analytical-quality concerns and emit one evidence item per variant.

    No variant is ever dropped. Records are returned sorted by
    ``coordinate.sort_key()`` with their evidence in the same order.
    """
    timestamp = clock.now()
    assessed: list[VariantRecord] = []
    evidence: list[EvidenceItem] = []

    for record in sorted(variants, key=lambda item: item.sort_key()):
        flags = _flags_for(record, thresholds)
        flagged = record.with_qc_flags(*flags) if flags else record
        assessed.append(flagged)
        evidence.append(
            _evidence_for(flagged, thresholds=thresholds, timestamp=timestamp, version=tool_version)
        )

    return QcResult(
        variants=tuple(assessed),
        evidence=tuple(evidence),
        metrics=_aggregate_metrics(assessed, evidence),
    )


# ---------------------------------------------------------------------------
# Flagging
# ---------------------------------------------------------------------------


def _flags_for(record: VariantRecord, thresholds: QualityThresholds) -> tuple[str, ...]:
    genotype = record.genotype
    raised: set[str] = set()

    if record.filter_status is not FilterStatus.PASS:
        raised.add(FLAG_FILTERED_BY_CALLER)
    if genotype.depth is not None and genotype.depth < thresholds.min_depth:
        raised.add(FLAG_LOW_DEPTH)
    if (
        genotype.genotype_quality is not None
        and genotype.genotype_quality < thresholds.min_genotype_quality
    ):
        raised.add(FLAG_LOW_GQ)

    fraction = allele_fraction(record)
    if fraction is not None and genotype.zygosity is Zygosity.HET:
        balance_flag = _allele_balance_flag(fraction, thresholds)
        if balance_flag is not None:
            raised.add(balance_flag)

    if (
        genotype.depth is None
        and genotype.genotype_quality is None
        and genotype.allele_balance is None
    ):
        raised.add(FLAG_NO_QUALITY_METRICS)

    return tuple(flag for flag in FLAG_ORDER if flag in raised)


def allele_fraction(record: VariantRecord) -> float | None:
    """The ALT read fraction in a form the het band can actually be applied to.

    For a biallelic call this is exactly ``Genotype.allele_balance``.

    For a record decomposed from a multiallelic site it is not, and the difference
    matters. After splitting ``A>G,GT`` with ``GT=1/2`` and ``AD=2,21,21``, allele 1
    holds ``ref_reads=2`` and ``alt_reads=21`` — but those two reads are the *site's*
    REF depth and exclude the 21 reads carrying the other ALT. ``alt/(ref+alt)`` is
    then 0.91 and looks homozygous, so a textbook compound-heterozygous call would
    be flagged ``high_allele_balance`` and down-ranked for an artifact of the
    representation. Where the site depth is known, the site allele fraction
    ``alt/DP`` (0.48 here) is banded instead; where it is not, ``None`` is returned
    and no allele-balance flag is raised, because a wrong flag is worse than none
    (GP-14).
    """
    genotype = record.genotype
    balance = genotype.allele_balance
    if OP_SPLIT_MULTIALLELIC not in record.normalisation_ops:
        return balance
    depth, alt_reads, ref_reads = genotype.depth, genotype.alt_reads, genotype.ref_reads
    if depth is None or alt_reads is None:
        return None
    if depth <= 0 or depth < (ref_reads or 0) + alt_reads:
        return balance
    return alt_reads / depth


def _allele_balance_flag(balance: float, thresholds: QualityThresholds) -> str | None:
    """Classify a heterozygous allele balance.

    ``possible_mosaic`` deliberately takes precedence over ``low_allele_balance``:
    in a mosaic aneuploidy disorder a skewed het is a candidate finding, and
    labelling it "low quality" would bury the signal this pipeline exists to find.
    """
    if balance > thresholds.max_allele_balance_het:
        return FLAG_HIGH_ALLELE_BALANCE
    if balance >= thresholds.min_allele_balance_het:
        return None
    if balance >= thresholds.mosaic_allele_balance_floor:
        return FLAG_POSSIBLE_MOSAIC
    return FLAG_LOW_ALLELE_BALANCE


# ---------------------------------------------------------------------------
# Evidence
# ---------------------------------------------------------------------------


def _evidence_for(
    record: VariantRecord,
    *,
    thresholds: QualityThresholds,
    timestamp: datetime,
    version: str,
) -> EvidenceItem:
    concerns = tuple(flag for flag in record.qc_flags if flag in CONTRADICTING_FLAGS)
    claim = _claim_for(record)
    return EvidenceItem(
        evidence_id=make_evidence_id(
            subject_id=record.variant_id,
            category=EvidenceCategory.ANALYTICAL,
            claim=claim,
            tool=TOOL_NAME,
        ),
        subject_id=record.variant_id,
        subject_kind="variant",
        claim=claim,
        category=EvidenceCategory.ANALYTICAL,
        direction=_direction_for(record, concerns),
        strength=_strength_for(record, concerns, thresholds),
        evidence_type=EvidenceType.DIRECT_MEASUREMENT,
        tier=AssertionTier.OBSERVED_DATA,
        citation=None,
        method=(
            "Comparison of the caller-reported depth, genotype quality, allelic depths "
            "and FILTER for this call against the configured analytical thresholds. No "
            "re-alignment, re-genotyping or orthogonal confirmation was performed."
        ),
        tool=TOOL_NAME,
        tool_version=version,
        limitations=_limitations_for(record),
        timestamp=timestamp,
        numeric_value=_numeric_value_for(record),
        payload=_payload_for(record, thresholds),
    )


def _claim_for(record: VariantRecord) -> str:
    genotype = record.genotype
    verdict = (
        "no analytical concerns were raised"
        if not record.qc_flags
        else "analytical concerns raised: " + ", ".join(record.qc_flags)
    )
    return (
        f"{record.variant_id}: caller-reported depth={_render(genotype.depth)}, "
        f"GQ={_render(genotype.genotype_quality)}, "
        f"allele_fraction={_render_float(allele_fraction(record))}, "
        f"FILTER={record.filter_status.value} — {verdict}."
    )


def _direction_for(record: VariantRecord, concerns: tuple[str, ...]) -> EvidenceDirection:
    if concerns:
        return EvidenceDirection.CONTRADICTS
    if FLAG_NO_QUALITY_METRICS in record.qc_flags or FLAG_POSSIBLE_MOSAIC in record.qc_flags:
        return EvidenceDirection.NEUTRAL
    return EvidenceDirection.SUPPORTS


def _strength_for(
    record: VariantRecord, concerns: tuple[str, ...], thresholds: QualityThresholds
) -> EvidenceStrength:
    if REF_ALLELE_MISMATCH_FLAG in record.qc_flags:
        return EvidenceStrength.STRONG
    if concerns:
        hard = {FLAG_LOW_DEPTH, FLAG_LOW_GQ, FLAG_FILTERED_BY_CALLER}
        return (
            EvidenceStrength.MODERATE
            if hard.intersection(concerns)
            else EvidenceStrength.SUPPORTING
        )
    if FLAG_NO_QUALITY_METRICS in record.qc_flags:
        return EvidenceStrength.INSUFFICIENT
    if FLAG_POSSIBLE_MOSAIC in record.qc_flags:
        return EvidenceStrength.SUPPORTING
    genotype = record.genotype
    comfortable = (
        genotype.depth is not None
        and genotype.depth >= 2 * thresholds.min_depth
        and genotype.genotype_quality is not None
        and genotype.genotype_quality >= thresholds.min_genotype_quality
    )
    return EvidenceStrength.STRONG if comfortable else EvidenceStrength.MODERATE


def _limitations_for(record: VariantRecord) -> str:
    limitations = _BASE_LIMITATIONS
    if FLAG_POSSIBLE_MOSAIC in record.qc_flags:
        limitations += _MOSAIC_LIMITATION
    if REF_ALLELE_MISMATCH_FLAG in record.qc_flags:
        limitations += _REF_MISMATCH_LIMITATION
    if FLAG_NO_QUALITY_METRICS in record.qc_flags:
        limitations += _NO_METRICS_LIMITATION
    return limitations


def _numeric_value_for(record: VariantRecord) -> float | None:
    fraction = allele_fraction(record)
    if fraction is not None:
        return round(fraction, 6)
    return float(record.genotype.depth) if record.genotype.depth is not None else None


def _payload_for(
    record: VariantRecord, thresholds: QualityThresholds
) -> dict[str, str | int | float | bool | None]:
    """Structured detail. Aggregate-safe fields only — never the raw GT string."""
    balance = record.genotype.allele_balance
    fraction = allele_fraction(record)
    return {
        "depth": record.genotype.depth,
        "genotype_quality": record.genotype.genotype_quality,
        "allele_balance": None if balance is None else round(balance, 6),
        "allele_fraction": None if fraction is None else round(fraction, 6),
        "zygosity": record.genotype.zygosity.value,
        "phased": record.genotype.phased,
        "filter_status": record.filter_status.value,
        "qc_flags": ";".join(record.qc_flags),
        "qc_flag_count": len(record.qc_flags),
        "min_depth_threshold": thresholds.min_depth,
        "min_genotype_quality_threshold": thresholds.min_genotype_quality,
    }


def _render(value: int | None) -> str:
    return "unreported" if value is None else str(value)


def _render_float(value: float | None) -> str:
    return "unreported" if value is None else f"{value:.3f}"


# ---------------------------------------------------------------------------
# Aggregate metrics
# ---------------------------------------------------------------------------


def _aggregate_metrics(
    variants: Sequence[VariantRecord], evidence: Sequence[EvidenceItem]
) -> dict[str, int | float]:
    depths = [v.genotype.depth for v in variants if v.genotype.depth is not None]
    het_fractions = [
        fraction
        for fraction in (allele_fraction(v) for v in variants if v.genotype.is_heterozygous)
        if fraction is not None
    ]
    flag_names = (*FLAG_ORDER, REF_ALLELE_MISMATCH_FLAG)

    metrics: dict[str, int | float] = {
        "total_variants": len(variants),
        "flagged_variants": sum(1 for v in variants if v.qc_flags),
        "unflagged_variants": sum(1 for v in variants if not v.qc_flags),
    }
    for name in flag_names:
        metrics[f"flag_{name}"] = sum(1 for v in variants if name in v.qc_flags)

    metrics["evidence_items"] = len(evidence)
    for direction in (
        EvidenceDirection.SUPPORTS,
        EvidenceDirection.CONTRADICTS,
        EvidenceDirection.NEUTRAL,
    ):
        metrics[f"evidence_{direction.value}"] = sum(
            1 for item in evidence if item.direction is direction
        )

    metrics["variants_with_depth"] = len(depths)
    metrics["variants_with_het_allele_fraction"] = len(het_fractions)
    if depths:
        metrics["mean_depth"] = round(statistics.fmean(depths), 3)
        metrics["median_depth"] = round(float(statistics.median(depths)), 3)
        metrics["min_observed_depth"] = min(depths)
        metrics["max_observed_depth"] = max(depths)
    if het_fractions:
        metrics["mean_het_allele_fraction"] = round(statistics.fmean(het_fractions), 6)
    return metrics
