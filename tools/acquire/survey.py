"""Cross-referencing the declarative registry against what is actually on disk.

``tools.acquire.catalog.KNOWN_RESOURCES`` says what this tool knows how to fetch.
This module says what is *actually true right now* for each of those entries: a
complete, hash-verified local file (``FETCHED``); or nothing trustworthy yet
(``NOT_FETCHED``, with ``notes`` saying why -- never started, still growing, or
present but wrong).

This is the only place in the tool that decides a download is "done". Getting that
decision wrong in either direction is bad: calling a growing file done bakes a wrong
hash into a file this project commits (see ``knowledge/manifests/resources.yaml``'s
own header); calling a genuinely finished file "still growing" just means a
re-inspection is needed later, which costs nothing.
"""

from __future__ import annotations

import time
from collections.abc import Callable, Iterable
from datetime import UTC, datetime
from pathlib import Path
from typing import Final

from mva.determinism import hash_file
from tools.acquire.fetch import is_download_stable, sniff_content_mismatch
from tools.acquire.models import ResourceEntry, ResourceStatus

#: How long to watch a present-but-unverified file for growth before trusting it.
#: Short by design: a real in-progress download (multi-hundred-MB gnomAD VCFs
#: observed growing at tens of MB/s during this tool's own development) moves far
#: enough in this window to be caught; a genuinely finished file costs one stat-sleep-stat
#: round trip per resource.
DEFAULT_STABILITY_INTERVAL_SECONDS: Final = 2.0


def _mtime_date(path: Path) -> str:
    """The file's own last-modified date, as an ISO date string.

    Used as ``retrieved`` instead of "whenever this function happened to run":
    re-running the survey without re-downloading anything must not change the
    committed manifest's retrieval dates.
    """
    return datetime.fromtimestamp(path.stat().st_mtime, tz=UTC).date().isoformat()


def survey_resource(
    resource_root: Path,
    entry: ResourceEntry,
    *,
    stability_interval: float = DEFAULT_STABILITY_INTERVAL_SECONDS,
    sleep: Callable[[float], None] = time.sleep,
) -> ResourceEntry:
    """Return ``entry`` updated to reflect the real state of its local file, if any."""
    dest = resource_root / entry.path

    if not dest.is_file():
        return entry.model_copy(
            update={
                "status": ResourceStatus.NOT_FETCHED,
                "sha256": None,
                "size_bytes": None,
                "retrieved": None,
                "notes": "not yet fetched: no local file at the declared path",
            }
        )

    if not is_download_stable(dest, interval_seconds=stability_interval, sleep=sleep):
        observed = dest.stat().st_size
        return entry.model_copy(
            update={
                "status": ResourceStatus.NOT_FETCHED,
                "sha256": None,
                "size_bytes": None,
                "retrieved": None,
                "notes": (
                    f"download in progress: size grew during a {stability_interval:g}s "
                    f"observation window (observed {observed} bytes); refusing to hash a "
                    "moving target"
                ),
            }
        )

    mismatch = sniff_content_mismatch(dest)
    if mismatch is not None:
        return entry.model_copy(
            update={
                "status": ResourceStatus.NOT_FETCHED,
                "sha256": None,
                "size_bytes": None,
                "retrieved": None,
                "notes": f"local file present but failed a content sanity check: {mismatch}",
            }
        )

    return entry.model_copy(
        update={
            "status": ResourceStatus.FETCHED,
            "sha256": hash_file(dest),
            "size_bytes": dest.stat().st_size,
            "retrieved": _mtime_date(dest),
            "notes": "",
        }
    )


def survey_all(
    resource_root: Path,
    entries: Iterable[ResourceEntry],
    *,
    stability_interval: float = DEFAULT_STABILITY_INTERVAL_SECONDS,
    sleep: Callable[[float], None] = time.sleep,
) -> tuple[ResourceEntry, ...]:
    """:func:`survey_resource` over every declared entry, in declaration order."""
    return tuple(
        survey_resource(resource_root, entry, stability_interval=stability_interval, sleep=sleep)
        for entry in entries
    )
