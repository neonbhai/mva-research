"""Cross-referencing the declarative registry against what is actually on disk.

``tools.acquire.catalog.KNOWN_RESOURCES`` says what this tool knows how to fetch.
This module says what is *actually true right now* for each of those entries: a
complete, hash-pinned, format-verified local artifact (``FETCHED``); or nothing
trustworthy yet (``NOT_FETCHED``, with ``notes`` saying why — never started, still
growing, or present but wrong).

This is the only place in the tool that decides a download is "done". Getting that
decision wrong in either direction is bad: calling a growing file done bakes a wrong
hash into a file this project commits (see ``knowledge/manifests/resources.yaml``'s
own header); calling a genuinely finished file "still growing" just means a
re-inspection is needed later, which costs nothing.

Registering a resource performs, in order:

1. a stability window, so a file still being written is never hashed;
2. a shallow magic-number sniff (``sniff_content_mismatch``);
3. a **deep format probe** — opening the file as the format its name claims
   (``tools.acquire.formats``);
4. a **full sha256** over every byte, plus the sampled ``spot_sha256`` the run-time
   check will compare against (``mva.resources``).

Step 3 is not redundant with step 4. A hash pins whatever arrived; only step 3
notices that what arrived was an HTML error page. Both outcomes are written into
the manifest, so a reader can see what was proven rather than inferring it.
"""

from __future__ import annotations

import time
from collections.abc import Callable, Iterable
from datetime import UTC, datetime
from pathlib import Path
from typing import Final

from mva.resources import SPOT_PLAN, IntegrityRecord, ResourceKind, ResourceStatus
from tools.acquire.digest import DigestCache
from tools.acquire.fetch import is_download_stable, sniff_content_mismatch
from tools.acquire.formats import probe_format
from tools.acquire.models import ResourceEntry

#: How long to watch a present-but-unverified file for growth before trusting it.
#: Short by design: a real in-progress download (multi-hundred-MB gnomAD VCFs
#: observed growing at tens of MB/s during this tool's own development) moves far
#: enough in this window to be caught; a genuinely finished file costs one stat-sleep-stat
#: round trip per resource.
DEFAULT_STABILITY_INTERVAL_SECONDS: Final = 2.0

_UNFETCHED_FIELDS: Final[dict[str, None]] = {
    "sha256": None,
    "size_bytes": None,
    "retrieved": None,
    "integrity": None,
}


def _mtime_date(path: Path) -> str:
    """The artifact's own last-modified date, as an ISO date string.

    Used as ``retrieved`` instead of "whenever this function happened to run":
    re-running the survey without re-downloading anything must not change the
    committed manifest's retrieval dates.
    """
    if path.is_dir():
        members = [p for p in path.rglob("*") if p.is_file()]
        stamp = max((p.stat().st_mtime for p in members), default=path.stat().st_mtime)
    else:
        stamp = path.stat().st_mtime
    return datetime.fromtimestamp(stamp, tz=UTC).date().isoformat()


def _directory_size(path: Path) -> int:
    return sum(p.stat().st_size for p in path.rglob("*") if p.is_file() and not p.is_symlink())


def _not_fetched(entry: ResourceEntry, note: str) -> ResourceEntry:
    return entry.model_copy(
        update={"status": ResourceStatus.NOT_FETCHED, **_UNFETCHED_FIELDS, "notes": note}
    )


def survey_resource(
    resource_root: Path,
    entry: ResourceEntry,
    *,
    stability_interval: float = DEFAULT_STABILITY_INTERVAL_SECONDS,
    sleep: Callable[[float], None] = time.sleep,
    cache: DigestCache | None = None,
    verified_at: str | None = None,
) -> ResourceEntry:
    """Return ``entry`` updated to reflect the real state of its local artifact, if any.

    ``cache`` memoises full digests across a registration pass so an interrupted
    202.8 GB run resumes rather than restarts; see ``tools.acquire.digest`` for why
    that cache is emphatically not an integrity check. ``verified_at`` is injected
    rather than read from the wall clock so a re-survey of unchanged bytes produces
    an unchanged manifest (GP-30).
    """
    dest = resource_root / entry.path
    exists = dest.is_dir() if entry.kind is ResourceKind.DIRECTORY else dest.is_file()

    if not exists:
        noun = "directory" if entry.kind is ResourceKind.DIRECTORY else "file"
        return _not_fetched(entry, f"not yet fetched: no local {noun} at the declared path")

    if entry.kind is ResourceKind.FILE and not is_download_stable(
        dest, interval_seconds=stability_interval, sleep=sleep
    ):
        observed = dest.stat().st_size
        return _not_fetched(
            entry,
            f"download in progress: size grew during a {stability_interval:g}s observation "
            f"window (observed {observed} bytes); refusing to hash a moving target",
        )

    if entry.kind is ResourceKind.FILE:
        mismatch = sniff_content_mismatch(dest)
        if mismatch is not None:
            return _not_fetched(
                entry, f"local file present but failed a content sanity check: {mismatch}"
            )

    probe = probe_format(dest, entry.kind)
    if not probe.ok:
        return _not_fetched(
            entry, f"local artifact present but failed its format check: {probe.problem}"
        )

    digests = cache if cache is not None else DigestCache(Path("/nonexistent"), enabled=False)
    full_sha256, sampled = digests.digests(dest, entry.kind)
    size = _directory_size(dest) if entry.kind is ResourceKind.DIRECTORY else dest.stat().st_size

    return entry.model_copy(
        update={
            "status": ResourceStatus.FETCHED,
            "sha256": full_sha256,
            "size_bytes": size,
            "retrieved": _mtime_date(dest),
            "integrity": IntegrityRecord(
                verified_at=verified_at or _mtime_date(dest),
                spot_plan=SPOT_PLAN,
                spot_sha256=sampled,
                format_check=probe.check,
                format_detail=probe.detail,
                index_check=probe.index_check,
                index_detail=probe.index_detail,
            ),
            "notes": "",
        }
    )


def survey_all(
    resource_root: Path,
    entries: Iterable[ResourceEntry],
    *,
    stability_interval: float = DEFAULT_STABILITY_INTERVAL_SECONDS,
    sleep: Callable[[float], None] = time.sleep,
    cache: DigestCache | None = None,
    verified_at: str | None = None,
) -> tuple[ResourceEntry, ...]:
    """:func:`survey_resource` over every declared entry, in declaration order."""
    return tuple(
        survey_resource(
            resource_root,
            entry,
            stability_interval=stability_interval,
            sleep=sleep,
            cache=cache,
            verified_at=verified_at,
        )
        for entry in entries
    )
