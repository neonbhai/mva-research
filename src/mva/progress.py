"""Say how long a stage will take and how much disk it will need — in its first minute.

A whole-genome run was killed after 7 h 30 m. It had annotated 38.5% of 4,962,060
variants and was projecting ~19.5 hours and ~95 GB of intermediates onto a 60 GB
volume. It could not have finished, and every fact needed to know that was
available inside the first two minutes: the input count came from ingestion, and
the per-record time and bytes came from the first batch.

Nothing measured them. The only way to know the stage was alive at all was to
sample its stack.

This module is the instrument. It is deliberately one object rather than two,
because the projection and the progress line are the same measurement asked at
different intervals: records seen, bytes written, seconds elapsed, against a total
that ingestion already counted.

**A stage that cannot possibly complete should say so in its first minute, not its
twentieth hour.** :meth:`StageMonitor.observe` raises
:class:`~mva.errors.InsufficientDiskError` the moment the projection exceeds free
space, and the message carries both numbers and the two ways out.

Determinism (GP-30)
-------------------
Everything here is wall-clock and machine-dependent, and **none of it reaches an
artifact**. Output goes to a caller-supplied stream, stderr by default, so it
cannot contaminate stdout for a command whose stdout is a CSV. The monitor never
returns a value that a stage's result depends on, and the only control flow it
influences is the refusal — which stops a run rather than changing one. Two runs
of the same case therefore print different numbers and produce identical bytes.

The clock is :func:`time.monotonic` by default and injectable, so tests drive it
without sleeping and ``mva.clock`` stays what it is: the source of *recorded*
time. A duration is not a timestamp and must not come from a clock that a demo
run freezes.
"""

from __future__ import annotations

import os
import shutil
import sys
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Final, TextIO

from mva.errors import InsufficientDiskError

#: Set to ``0`` to silence progress output. The refusal is NOT silenced by it: a
#: run that cannot complete must fail whether or not anyone wanted the chatter.
PROGRESS_ENV: Final = "MVA_PROGRESS"

#: Seconds between progress lines. Long enough that a stage streaming 100k
#: records/s does not spend measurable time formatting, short enough that a human
#: watching a terminal sees movement.
DEFAULT_INTERVAL_SECONDS: Final = 30.0

#: Records that must be seen before a projection is trusted enough to refuse on.
#: A twelve-record fixture cannot project anything, and refusing a run on a sample
#: of one is worse than not refusing at all.
DEFAULT_MIN_SAMPLE_RECORDS: Final = 1_000

#: Fraction of free space left unclaimed. The projection is a projection, the
#: volume holds more than this stage's output, and a run that fills a disk to the
#: last byte takes the machine down with it rather than just itself.
DEFAULT_HEADROOM_FRACTION: Final = 0.10

#: Records between clock reads. `observe` is called once per record — 4.96 M times
#: on a whole-genome run — and `time.monotonic()` is a syscall on some platforms,
#: so consulting it per record would make the progress meter a measurable share of
#: the stage it is measuring. A stride of 4,096 puts the check at roughly 25 Hz for
#: a stage streaming 100k records/s, which is far more often than the 30-second
#: report interval needs and still ~4,000x cheaper than per-record.
_REPORT_STRIDE: Final = 4_096


def free_bytes_at(path: Path) -> int:
    """Free bytes on the volume holding ``path``, or its nearest existing parent.

    Returns 0 when the volume cannot be interrogated. Zero disables the refusal
    rather than triggering it: a stage that cannot see the disk must not conclude
    the disk is full, and "we could not look" is not "we looked and found
    nothing" (GP-14).
    """
    for candidate in [path, *path.parents]:
        if candidate.exists():
            try:
                return shutil.disk_usage(candidate).free
            except OSError:
                return 0
    return 0


def _human_bytes(value: float) -> str:
    """Bytes as a short decimal string. GB, not GiB — disks are sold in GB."""
    for unit, scale in (("TB", 1e12), ("GB", 1e9), ("MB", 1e6), ("kB", 1e3)):
        if abs(value) >= scale:
            return f"{value / scale:.1f} {unit}"
    return f"{value:.0f} B"


def _human_duration(seconds: float) -> str:
    """A duration a human reads at a glance: ``19h 32m``, ``4m 07s``, ``12s``."""
    if seconds < 60:
        return f"{seconds:.0f}s"
    minutes, secs = divmod(int(seconds), 60)
    if minutes < 60:
        return f"{minutes}m {secs:02d}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h {minutes:02d}m"


@dataclass(frozen=True, slots=True)
class Projection:
    """What the stage will cost, extrapolated linearly from what it has done.

    Linear because the work is linear: one adapter lookup and one artifact row per
    record, batched at a fixed size. The projection is stated as a projection
    everywhere it is printed, and :attr:`is_trustworthy` says when the sample is
    too small to act on.
    """

    stage: str
    records_seen: int
    total_records: int
    elapsed_seconds: float
    bytes_written: int
    free_bytes: int
    min_sample_records: int = DEFAULT_MIN_SAMPLE_RECORDS
    headroom_fraction: float = DEFAULT_HEADROOM_FRACTION

    @property
    def records_per_second(self) -> float:
        return self.records_seen / self.elapsed_seconds if self.elapsed_seconds > 0 else 0.0

    @property
    def bytes_per_record(self) -> float:
        return self.bytes_written / self.records_seen if self.records_seen else 0.0

    @property
    def fraction_done(self) -> float:
        return self.records_seen / self.total_records if self.total_records else 0.0

    @property
    def projected_seconds(self) -> float:
        rate = self.records_per_second
        return self.total_records / rate if rate > 0 else 0.0

    @property
    def remaining_seconds(self) -> float:
        return max(0.0, self.projected_seconds - self.elapsed_seconds)

    @property
    def projected_bytes(self) -> int:
        return int(self.bytes_per_record * self.total_records)

    @property
    def budget_bytes(self) -> int:
        """Free space this stage is allowed to claim, after headroom."""
        return int(self.free_bytes * (1.0 - self.headroom_fraction))

    @property
    def is_trustworthy(self) -> bool:
        """Whether the sample is large enough to refuse a run on."""
        return self.records_seen >= self.min_sample_records and self.elapsed_seconds > 0

    @property
    def fits(self) -> bool:
        """Whether the projected output fits in the budget.

        True when the disk could not be interrogated (``free_bytes == 0``) — see
        :func:`free_bytes_at`. True when nothing has been written, because a stage
        that produces no artifact cannot exhaust a volume.
        """
        if self.free_bytes <= 0 or self.projected_bytes <= 0:
            return True
        return self.projected_bytes <= self.budget_bytes

    def as_lines(self) -> tuple[str, ...]:
        """The projection as printed. Every number is labelled as measured or projected."""
        return (
            f"[{self.stage}] pre-flight projection from {self.records_seen:,} of "
            f"{self.total_records:,} records ({self.fraction_done:.1%} sampled)",
            f"[{self.stage}]   measured : {self.records_per_second:,.0f} rec/s, "
            f"{_human_bytes(self.bytes_per_record)}/record, "
            f"{_human_bytes(self.bytes_written)} written in "
            f"{_human_duration(self.elapsed_seconds)}",
            f"[{self.stage}]   projected: {_human_duration(self.projected_seconds)} total, "
            f"{_human_bytes(self.projected_bytes)} of output",
            f"[{self.stage}]   disk     : {_human_bytes(self.free_bytes)} free, "
            f"{_human_bytes(self.budget_bytes)} budget after "
            f"{self.headroom_fraction:.0%} headroom",
        )

    def refusal_message(self) -> str:
        """Why the run is being stopped, and the two ways to make it possible."""
        return (
            f"The {self.stage} stage cannot complete on this volume, and stopping now "
            f"costs minutes rather than hours.\n"
            f"  measured   {self.records_seen:,} of {self.total_records:,} records in "
            f"{_human_duration(self.elapsed_seconds)} "
            f"({self.records_per_second:,.0f} rec/s, "
            f"{_human_bytes(self.bytes_per_record)} per record)\n"
            f"  projected  {_human_bytes(self.projected_bytes)} of output over "
            f"{_human_duration(self.projected_seconds)}\n"
            f"  available  {_human_bytes(self.free_bytes)} free, "
            f"{_human_bytes(self.budget_bytes)} after "
            f"{self.headroom_fraction:.0%} headroom\n"
            "Either give the workspace a larger volume, or tighten the annotation "
            "relevance bar so fewer records are materialised: `triage.min_impact: "
            "moderate` in the case config keeps HIGH and MODERATE only, and "
            "`triage.retain_without_gene_assignment` is already false by default "
            "(ADR 0028). The discard ledger records every variant either setting "
            "drops, by id and by rule, so nothing becomes invisible.\n"
            "No coordinate is named in this message (PRIV-09)."
        )


class StageMonitor:
    """Throughput, ETA and a disk refusal for one long stage.

    Fed one call per record — :meth:`observe` — and does almost nothing on most of
    them: an integer add, and a clock read only after the interval could plausibly
    have elapsed. The refusal fires once, on the first trustworthy projection that
    does not fit, which for a whole-genome annotate pass is at the end of the first
    batch.
    """

    __slots__ = (
        "_bytes",
        "_checked",
        "_enabled",
        "_free_bytes",
        "_headroom",
        "_interval",
        "_last_report",
        "_min_sample",
        "_next_check_at",
        "_now",
        "_out",
        "_project_after",
        "_records",
        "_stage",
        "_started",
        "_total",
    )

    def __init__(
        self,
        stage: str,
        *,
        total_records: int,
        free_bytes: int,
        out: TextIO | None = None,
        now: Callable[[], float] | None = None,
        interval_seconds: float = DEFAULT_INTERVAL_SECONDS,
        project_after_records: int = DEFAULT_MIN_SAMPLE_RECORDS,
        min_sample_records: int = DEFAULT_MIN_SAMPLE_RECORDS,
        headroom_fraction: float = DEFAULT_HEADROOM_FRACTION,
        enabled: bool | None = None,
    ) -> None:
        import time  # noqa: PLC0415 - injected in tests; not a module-level dependency

        self._stage = stage
        self._total = max(0, total_records)
        self._free_bytes = max(0, free_bytes)
        self._out = out if out is not None else sys.stderr
        self._now = now or time.monotonic
        self._interval = interval_seconds
        self._project_after = max(1, project_after_records)
        self._min_sample = max(1, min_sample_records)
        self._headroom = headroom_fraction
        self._enabled = _progress_enabled() if enabled is None else enabled
        self._records = 0
        self._bytes = 0
        self._started = self._now()
        self._last_report = self._started
        self._next_check_at = self._project_after
        self._checked = False

    # ------------------------------------------------------------------ driving

    def observe(self, *, records: int = 1, bytes_written: int = 0) -> None:
        """Account for work done, and refuse the run if the projection does not fit.

        Called once per record on the hot path, so the common case is two integer
        adds, a comparison and a modulo. The clock is read at most once every
        :data:`_REPORT_STRIDE` records; the projection — and with it the refusal —
        is computed once every ``project_after_records``.
        """
        self._records += records
        self._bytes += bytes_written
        if not records:
            # Bytes-only: the record these bytes belong to was counted when it was
            # examined, upstream of any bar. Accounting only -- no clock read and no
            # projection, so a caller may report bytes as often as it likes.
            return
        if self._records >= self._next_check_at:
            self._next_check_at = self._records + self._project_after
            self._enforce()
        elif self._enabled and self._records % _REPORT_STRIDE == 0:
            self._maybe_report()

    def _enforce(self) -> None:
        projection = self.projection()
        if not self._checked:
            self._checked = True
            self._emit(projection.as_lines())
            if projection.is_trustworthy and not projection.fits:
                raise InsufficientDiskError(projection.refusal_message())
        elif self._enabled:
            self._report(projection)

    def _maybe_report(self) -> None:
        now = self._now()
        if now - self._last_report < self._interval:
            return
        self._last_report = now
        self._report(self.projection())

    def _report(self, projection: Projection) -> None:
        self._last_report = self._now()
        self._emit(
            (
                f"[{self._stage}] {projection.records_seen:,}/{projection.total_records:,} "
                f"({projection.fraction_done:.1%}) "
                f"{projection.records_per_second:,.0f} rec/s "
                f"eta {_human_duration(projection.remaining_seconds)} "
                f"{_human_bytes(projection.bytes_written)} written",
            )
        )

    def _emit(self, lines: tuple[str, ...]) -> None:
        if not self._enabled:
            return
        for line in lines:
            self._out.write(line + "\n")
        self._out.flush()

    # ------------------------------------------------------------------ reporting

    def projection(self) -> Projection:
        """The current extrapolation. Cheap; allocates one frozen dataclass."""
        return Projection(
            stage=self._stage,
            records_seen=self._records,
            total_records=max(self._total, self._records),
            elapsed_seconds=max(0.0, self._now() - self._started),
            bytes_written=self._bytes,
            free_bytes=self._free_bytes,
            min_sample_records=self._min_sample,
            headroom_fraction=self._headroom,
        )

    def finish(self) -> Projection:
        """Final line, and the measurement the caller may want to record.

        The returned projection is the stage's actual cost, not an extrapolation:
        ``records_seen == total_records`` once the stream is drained.
        """
        projection = self.projection()
        self._emit(
            (
                f"[{self._stage}] done: {projection.records_seen:,} records in "
                f"{_human_duration(projection.elapsed_seconds)} "
                f"({projection.records_per_second:,.0f} rec/s), "
                f"{_human_bytes(projection.bytes_written)} written",
            )
        )
        return projection

    @property
    def records(self) -> int:
        return self._records

    @property
    def bytes_written(self) -> int:
        return self._bytes


def _progress_enabled() -> bool:
    """Progress is on unless ``MVA_PROGRESS=0``, and off when stderr is not a stream."""
    return os.environ.get(PROGRESS_ENV, "1") not in {"0", "false", "no"}
