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
from array import array
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from datetime import datetime
from typing import Final

from mva.clock import Clock
from mva.config import QualityThresholds
from mva.errors import IngestionError
from mva.ingestion.normalise import REF_ALLELE_MISMATCH_FLAG
from mva.models.base import AssertionTier, error_token
from mva.models.evidence import (
    EvidenceCategory,
    EvidenceDirection,
    EvidenceItem,
    EvidenceStrength,
    EvidenceType,
    make_evidence_id,
)
from mva.models.genome import contig_sort_key
from mva.models.variant import (
    FLAG_FILTERED_BY_CALLER,
    FLAG_HIGH_ALLELE_BALANCE,
    FLAG_LOW_ALLELE_BALANCE,
    FLAG_LOW_DEPTH,
    FLAG_LOW_GQ,
    FLAG_NO_CALLER_FILTER,
    FLAG_NO_QUALITY_METRICS,
    FLAG_POSSIBLE_MOSAIC,
    FilterStatus,
    VariantRecord,
    Zygosity,
)

# ---------------------------------------------------------------------------
# Flag vocabulary
# ---------------------------------------------------------------------------

# The flag names themselves are imported from :mod:`mva.models.variant`, which
# owns the vocabulary shared with the prioritisation stage; only the emission
# ORDER and the contradiction classification are this module's business.

#: Canonical emission order, so ``qc_flags`` is deterministic (GP-30).
FLAG_ORDER: Final[tuple[str, ...]] = (
    FLAG_FILTERED_BY_CALLER,
    FLAG_NO_CALLER_FILTER,
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
_NO_CALLER_FILTER_LIMITATION: Final = (
    " The FILTER column was absent, so the caller recorded no opinion either way. That "
    "is not a passing call and not a rejected one; nothing about site-level filtering "
    "can be concluded from this record."
)

#: Flags that are noted but argue neither for nor against the call. Reported so a
#: reader can see they were considered, and kept out of the concern list so the
#: claim sentence does not describe them as problems.
NEUTRAL_FLAGS: Final[frozenset[str]] = frozenset(
    {FLAG_NO_CALLER_FILTER, FLAG_POSSIBLE_MOSAIC, FLAG_NO_QUALITY_METRICS}
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


@dataclass(frozen=True, slots=True)
class AssessedVariant:
    """One record and the analytical evidence item that describes it."""

    variant: VariantRecord
    evidence: EvidenceItem


# ---------------------------------------------------------------------------
# Entry points
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

    This holds the whole callset, its flagged copy and one evidence item each, all
    at once — measured at 7.13 KB per record, which is ~32 GB for a WGS callset
    (``docs/scale-report.md`` §2). It is the right shape at fixture scale and the
    wrong one at genome scale; :func:`iter_assessed` is the same work, streamed.
    """
    stream = iter_assessed(
        sorted(variants, key=lambda item: item.sort_key()),
        thresholds=thresholds,
        clock=clock,
        tool_version=tool_version,
    )
    assessed: list[VariantRecord] = []
    evidence: list[EvidenceItem] = []
    for item in stream:
        assessed.append(item.variant)
        evidence.append(item.evidence)
    return QcResult(variants=tuple(assessed), evidence=tuple(evidence), metrics=stream.metrics())


def iter_assessed(
    variants: Iterable[VariantRecord],
    *,
    thresholds: QualityThresholds,
    clock: Clock,
    tool_version: str = "mva-qc/0.1.0",
) -> QcStream:
    """Assess a stream of records without holding them.

    The caller supplies records **already in ``sort_key()`` order** — which is what
    :func:`mva.ingestion.reader.iter_vcf` emits and what
    :func:`~mva.ingestion.normalise.normalise_variants` returns for each chunk of
    such a stream. This function cannot sort them, so it *verifies* the order
    instead of assuming it: a record that arrives out of order raises rather than
    producing a run whose artifact is silently in a different order from the one
    :func:`assess_quality` would have written (GP-30).

    Metrics that need every value — the median depth, and the exact ``fsum``-based
    means — are accumulated into two fixed-width arrays rather than lists of Python
    objects, so the whole summary costs 16 bytes per record and comes out
    bit-identical to the batch path.
    """
    return QcStream(variants, thresholds=thresholds, clock=clock, tool_version=tool_version)


class QcStream:
    """Flagged records and their evidence, one at a time; metrics at the end.

    Construct via :func:`iter_assessed`.
    """

    __slots__ = ("_metrics", "_source", "_started", "_thresholds", "_timestamp", "_version")

    def __init__(
        self,
        variants: Iterable[VariantRecord],
        *,
        thresholds: QualityThresholds,
        clock: Clock,
        tool_version: str,
    ) -> None:
        self._source = variants
        self._thresholds = thresholds
        # Sampled once, before any record is seen, so every evidence item in a run
        # carries the same timestamp however long the stream takes (GP-30).
        self._timestamp = clock.now()
        self._version = tool_version
        self._metrics = _MetricsAccumulator()
        self._started = False

    def __iter__(self) -> Iterator[AssessedVariant]:
        if self._started:
            msg = (
                "This QC stream has already been iterated. Re-iterating would "
                "double-count every metric. Call iter_assessed() again for a second "
                "pass over a re-read source."
            )
            raise IngestionError(msg)
        self._started = True
        return self._drive()

    def _drive(self) -> Iterator[AssessedVariant]:
        thresholds = self._thresholds
        metrics = self._metrics
        previous: tuple[int, int, str, str] | None = None
        ranks: dict[str, int] = {}
        for record in self._source:
            coordinate = record.coordinate
            rank = ranks.get(coordinate.contig)
            if rank is None:
                rank = contig_sort_key(coordinate.contig)
                ranks[coordinate.contig] = rank
            key = (rank, coordinate.position, coordinate.ref, coordinate.alt)
            if previous is not None and key < previous:
                raise _unsorted_error(coordinate.contig)
            previous = key

            flags = _flags_for(record, thresholds)
            flagged = record.with_qc_flags(*flags) if flags else record
            evidence = _evidence_for(
                flagged,
                thresholds=thresholds,
                timestamp=self._timestamp,
                version=self._version,
            )
            metrics.observe(flagged, evidence)
            yield AssessedVariant(variant=flagged, evidence=evidence)

    def metrics(self) -> dict[str, int | float]:
        """The same aggregate metrics :class:`QcResult` carries, for what was seen.

        Meaningful only once the stream has been consumed; a caller that stops early
        gets the metrics for the prefix it read, and the record counts inside say so.
        """
        return self._metrics.finalise()


def _unsorted_error(contig: str) -> IngestionError:
    return IngestionError(
        "Records reached analytical QC out of coordinate order. Streaming QC cannot "
        "reorder them — it never holds more than one — so it stops rather than write "
        "an artifact ordered differently from the batch path, which would break the "
        "byte-identical repeat-run guarantee (GP-30). Sort the input, or use "
        "assess_quality(), which sorts. Contig handle "
        f"<contig:{error_token(contig)}> — the coordinate is tokenised rather than "
        "echoed (PRIV-09)."
    )


# ---------------------------------------------------------------------------
# Flagging
# ---------------------------------------------------------------------------


def _flags_for(record: VariantRecord, thresholds: QualityThresholds) -> tuple[str, ...]:
    genotype = record.genotype
    raised: set[str] = set()

    # Only a real FILTER entry means the caller objected. A missing FILTER
    # column means the caller expressed no opinion at all, and reporting that as
    # "filtered by caller" invents a rejection nobody made — it produced a
    # CONTRADICTS item on immaculate DP=40/GQ=99 calls.
    if record.filter_status is FilterStatus.FILTERED:
        raised.add(FLAG_FILTERED_BY_CALLER)
    elif record.filter_status is FilterStatus.MISSING:
        raised.add(FLAG_NO_CALLER_FILTER)
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
    """The ALT read fraction the heterozygous band may be applied to.

    Delegates to :attr:`VariantRecord.allele_fraction`, which is where the
    site-aware formula lives so that the prioritisation stage can use the same
    one without importing this package (GP-01/GP-03). Kept as a function because
    it is this module's published name for the quantity.
    """
    return record.allele_fraction


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
    concerns = [flag for flag in record.qc_flags if flag in CONTRADICTING_FLAGS]
    noted = [flag for flag in record.qc_flags if flag not in CONTRADICTING_FLAGS]
    if concerns:
        verdict = "analytical concerns raised: " + ", ".join(concerns)
        if noted:
            verdict += " (also noted: " + ", ".join(noted) + ")"
    elif noted:
        # A flag that is not a concern must not be rendered as one. Saying
        # "analytical concerns raised: no_caller_filter" reads as a rejection.
        verdict = "no analytical concerns were raised (noted: " + ", ".join(noted) + ")"
    else:
        verdict = "no analytical concerns were raised"
    return (
        f"{record.variant_id}: caller-reported depth={_render(genotype.depth)}, "
        f"GQ={_render(genotype.genotype_quality)}, "
        f"allele_fraction={_render_float(allele_fraction(record))}, "
        f"FILTER={record.filter_status.value} — {verdict}."
    )


def _direction_for(record: VariantRecord, concerns: tuple[str, ...]) -> EvidenceDirection:
    if concerns:
        return EvidenceDirection.CONTRADICTS
    if NEUTRAL_FLAGS.intersection(record.qc_flags):
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
    if FLAG_NO_CALLER_FILTER in record.qc_flags:
        limitations += _NO_CALLER_FILTER_LIMITATION
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
    fraction = record.allele_fraction
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


#: Every flag the metrics block counts, in emission order.
_FLAG_NAMES: Final[tuple[str, ...]] = (*FLAG_ORDER, REF_ALLELE_MISMATCH_FLAG)

#: Evidence directions the metrics block counts, in reporting order.
_DIRECTIONS: Final[tuple[EvidenceDirection, ...]] = (
    EvidenceDirection.SUPPORTS,
    EvidenceDirection.CONTRADICTS,
    EvidenceDirection.NEUTRAL,
)


class _MetricsAccumulator:
    """Aggregate metrics built one record at a time.

    Every count is a running integer. The two quantities that genuinely need all
    the data — the median depth, and the ``math.fsum``-exact means behind
    ``statistics.fmean`` — are kept in ``array`` buffers of machine ints and
    doubles rather than Python lists, which is 8 bytes per value instead of ~36 and
    keeps a whole-genome summary at tens of megabytes. Feeding the same values to
    the same ``statistics`` functions is what makes the streamed metrics identical
    to the batch ones rather than merely close.
    """

    __slots__ = ("_depths", "_directions", "_evidence", "_flagged", "_flags", "_het", "_total")

    def __init__(self) -> None:
        self._total = 0
        self._flagged = 0
        self._evidence = 0
        self._flags: dict[str, int] = dict.fromkeys(_FLAG_NAMES, 0)
        self._directions: dict[str, int] = {d.value: 0 for d in _DIRECTIONS}
        self._depths = array("q")
        self._het = array("d")

    def observe(self, variant: VariantRecord, evidence: EvidenceItem) -> None:
        self._total += 1
        flags = variant.qc_flags
        if flags:
            self._flagged += 1
            for name in flags:
                if name in self._flags:
                    self._flags[name] += 1
        self._evidence += 1
        direction = evidence.direction.value
        if direction in self._directions:
            self._directions[direction] += 1

        depth = variant.genotype.depth
        if depth is not None:
            self._depths.append(depth)
        if variant.genotype.is_heterozygous:
            fraction = allele_fraction(variant)
            if fraction is not None:
                self._het.append(fraction)

    def finalise(self) -> dict[str, int | float]:
        metrics: dict[str, int | float] = {
            "total_variants": self._total,
            "flagged_variants": self._flagged,
            "unflagged_variants": self._total - self._flagged,
        }
        for name in _FLAG_NAMES:
            metrics[f"flag_{name}"] = self._flags[name]

        metrics["evidence_items"] = self._evidence
        for direction in _DIRECTIONS:
            metrics[f"evidence_{direction.value}"] = self._directions[direction.value]

        depths = self._depths
        het_fractions = self._het
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
