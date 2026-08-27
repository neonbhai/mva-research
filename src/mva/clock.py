"""Injectable clock.

GP-30 requires byte-identical repeat runs. A wall-clock read embedded in an
artifact makes that impossible, so every timestamp in this codebase comes from a
`Clock` passed down from the composition root. Production uses `SystemClock`;
tests and the deterministic demo use `FixedClock`.

`tests/unit/test_architecture.py::test_no_naive_datetime_now` enforces that
nothing bypasses this module.
"""

from __future__ import annotations

from datetime import UTC, datetime
from typing import Protocol


class Clock(Protocol):
    """Source of the current time."""

    def now(self) -> datetime:
        """Timezone-aware current time."""
        ...


class SystemClock:
    """Real wall-clock time, always UTC-aware."""

    def now(self) -> datetime:
        return datetime.now(UTC)


class FixedClock:
    """A clock frozen at a chosen instant.

    Used by the synthetic demo so that `just demo` run twice produces identical
    manifests, and by tests that assert determinism.
    """

    def __init__(self, instant: datetime) -> None:
        if instant.tzinfo is None:
            msg = "FixedClock requires a timezone-aware datetime."
            raise ValueError(msg)
        self._instant = instant

    def now(self) -> datetime:
        return self._instant


#: The canonical instant used for deterministic synthetic runs.
DEMO_INSTANT = datetime(2026, 1, 1, 0, 0, 0, tzinfo=UTC)


def demo_clock() -> FixedClock:
    return FixedClock(DEMO_INSTANT)
