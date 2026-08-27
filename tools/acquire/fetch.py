"""Downloading a registered resource to the external resource root.

Two safety properties matter more here than anywhere else in this tool:

1. **Nothing is ever written inside this repository.** A ClinVar VCF is ~190MB and
   a single gnomAD exome sites file is hundreds of MB; committing one by accident
   is exactly the kind of mistake ``git`` makes irrecoverably easy. The containment
   check reuses ``mva.config.path_is_within`` -- the same ``(st_dev, st_ino)``-based
   check that keeps patient workspaces out of the repo (GP-40) -- rather than
   reimplementing a prefix comparison that a symlink or a case-folding quirk could
   quietly defeat.
2. **Only allowlisted hosts are ever contacted.** Enforced again here even though
   ``ResourceEntry`` already checks it at construction time (``tools.acquire.hosts``)
   -- defence in depth against a URL that reaches this function some other way.
"""

from __future__ import annotations

import os
import time
import urllib.request
from collections.abc import Callable
from datetime import UTC, datetime
from pathlib import Path
from typing import Final

from mva.config import find_repo_root, path_is_within
from mva.determinism import hash_file
from tools.acquire.errors import ResourceFetchError, ResourceRootError
from tools.acquire.hosts import assert_allowed_host
from tools.acquire.models import ResourceEntry, ResourceStatus

#: Streamed in fixed-size chunks, matching mva.determinism.hash_file's own chunk size,
#: so neither hashing nor downloading ever holds a whole multi-GB file in memory.
_CHUNK: Final = 1 << 20

#: $MVA_RESOURCES falls back to this. Chosen to match where this project's own
#: reference downloads actually live for the hackathon; override with $MVA_RESOURCES
#: for any other layout.
DEFAULT_RESOURCE_ROOT: Final[Path] = Path("~/Contri/bio-hackathon/mva-resources").expanduser()

#: Magic bytes for a gzip/BGZF member (BGZF is gzip with a required extra field --
#: every BGZF stream is also a valid gzip stream and starts with the same two bytes).
_GZIP_MAGIC: Final = b"\x1f\x8b"

#: The first bytes of an HTML error/landing page, lowercased. A resource is declared
#: by its true content type (.vcf.gz, .csv, .tsv, .obo, ...); a server returning one
#: of these instead of the promised bytes is the exact failure this check exists to
#: catch -- it happened for real during this tool's own registration pass (a stale
#: Gene2Phenotype bulk-download URL returned a JS app shell with HTTP 200).
_HTML_SNIFFS: Final[tuple[bytes, ...]] = (b"<!doctype html", b"<html")

_GZIP_SUFFIXES: Final[frozenset[str]] = frozenset({".gz", ".bgz", ".tbi"})


def resolve_resource_root(
    explicit: Path | str | None = None,
    *,
    env: dict[str, str] | None = None,
    repo_root: Path | None = None,
) -> Path:
    """Resolve the external resource root, refusing one that resolves inside the repo.

    Order: ``explicit`` argument, then ``$MVA_RESOURCES``, then
    :data:`DEFAULT_RESOURCE_ROOT`. ``env`` is injectable (defaults to
    ``os.environ``-like access via the caller) purely so tests never depend on the
    real process environment.
    """
    variables = env if env is not None else os.environ
    raw = explicit if explicit is not None else variables.get("MVA_RESOURCES")
    root = Path(raw).expanduser() if raw else DEFAULT_RESOURCE_ROOT
    repo = (repo_root or find_repo_root()).resolve()
    _assert_outside_repo(root, repo)
    return root


def _assert_outside_repo(path: Path, repo_root: Path) -> None:
    if path_is_within(path, repo_root):
        msg = (
            f"Resource root {path.as_posix()!r} resolves inside the repository "
            f"({repo_root.as_posix()!r}). Public reference downloads can be hundreds of "
            "megabytes to several gigabytes each; writing them inside the repo tree "
            "risks a 'git add -A' committing them permanently. Point $MVA_RESOURCES "
            "somewhere outside the repo."
        )
        raise ResourceRootError(msg)


def is_download_stable(
    path: Path,
    *,
    interval_seconds: float = 2.0,
    sleep: Callable[[float], None] = time.sleep,
) -> bool:
    """True if ``path``'s size is unchanged across a short interval.

    A file being written by a resumable ``curl -C -`` process grows monotonically
    while the download is in flight; comparing size before and after a short pause is
    the cheapest signal available that nothing is actively writing to it -- no PID
    file, no lock, and it works even when the writer is a process this tool never
    started (exactly the situation this tool was built against: gnomAD exome VCFs
    downloading in the background while the manifest was generated).

    ``sleep`` is injectable so tests can simulate concurrent growth deterministically
    instead of racing a real timer.
    """
    before = path.stat().st_size
    sleep(interval_seconds)
    after = path.stat().st_size
    return before == after


def sniff_content_mismatch(path: Path) -> str | None:
    """A human-readable reason ``path``'s bytes don't match what its name promises, or None.

    Deliberately shallow -- this is a tripwire against a server silently substituting
    an HTML error page for the promised file (observed for real: a stale
    Gene2Phenotype URL returned `text/html` with HTTP 200), not a format validator.
    """
    try:
        with path.open("rb") as handle:
            head = handle.read(64)
    except OSError as exc:
        return f"could not be read ({exc.__class__.__name__})"

    lowered = head.lower()
    if any(lowered.startswith(sniff) for sniff in _HTML_SNIFFS):
        return (
            "starts with an HTML document, not the expected data format (likely an "
            "error/landing page)"
        )

    if path.suffix in _GZIP_SUFFIXES and not head.startswith(_GZIP_MAGIC):
        return f"named {path.suffix!r} but does not start with the gzip/BGZF magic bytes"

    return None


def _today_iso() -> str:
    return datetime.now(UTC).date().isoformat()


def _http_download(url: str, dest: Path, *, resume_from: int, timeout: float) -> None:
    """The actual network transport. Isolated so tests can replace it without a socket.

    Sends nothing but the URL itself and a standard ``Range`` header -- no proband
    coordinate, sample ID or workspace path ever has a code path into this function's
    arguments (see ``ResourceEntry``, whose fields are exhaustively public-reference
    metadata).
    """
    request = urllib.request.Request(url, method="GET")  # noqa: S310 - host is allowlisted by the caller
    if resume_from:
        request.add_header("Range", f"bytes={resume_from}-")
    with urllib.request.urlopen(request, timeout=timeout) as response:  # noqa: S310 - see above
        resumed = resume_from and response.status == 206
        mode = "ab" if resumed else "wb"
        dest.parent.mkdir(parents=True, exist_ok=True)
        with dest.open(mode) as handle:
            while chunk := response.read(_CHUNK):
                handle.write(chunk)


def fetch_resource(
    entry: ResourceEntry,
    resource_root: Path,
    *,
    timeout: float = 120.0,
    max_attempts: int = 3,
    retry_delay: Callable[[float], None] = time.sleep,
    repo_root: Path | None = None,
) -> ResourceEntry:
    """Download ``entry`` into ``resource_root``, resuming a partial file if present.

    Returns a new :class:`ResourceEntry` (this type is frozen) with ``status``,
    ``sha256``, ``size_bytes`` and ``retrieved`` filled in from the bytes actually
    written to disk -- never from the caller's claim about what the resource is.

    Raises :class:`~tools.acquire.errors.DisallowedHostError` before opening any
    connection if ``entry.url``'s host is not on the allowlist, and
    :class:`~tools.acquire.errors.ResourceRootError` if ``resource_root`` resolves
    inside this repository. ``repo_root`` is injectable (defaults to
    :func:`mva.config.find_repo_root`) purely for testability.
    """
    assert_allowed_host(entry.url)
    repo = (repo_root or find_repo_root()).resolve()
    _assert_outside_repo(resource_root, repo)

    dest = resource_root / entry.path
    _assert_outside_repo(dest, repo)
    dest.parent.mkdir(parents=True, exist_ok=True)

    last_error: OSError | None = None
    for attempt in range(1, max_attempts + 1):
        resume_from = dest.stat().st_size if dest.is_file() else 0
        try:
            _http_download(entry.url, dest, resume_from=resume_from, timeout=timeout)
            last_error = None
            break
        except OSError as exc:  # network error, timeout, etc.
            last_error = exc
            if attempt < max_attempts:
                retry_delay(min(2**attempt, 30))
    if last_error is not None:
        msg = (
            f"Failed to fetch resource {entry.name!r} after {max_attempts} attempt(s): "
            f"{last_error.__class__.__name__}"
        )
        raise ResourceFetchError(msg) from last_error

    mismatch = sniff_content_mismatch(dest)
    if mismatch is not None:
        return entry.model_copy(
            update={
                "status": ResourceStatus.NOT_FETCHED,
                "sha256": None,
                "size_bytes": None,
                "retrieved": None,
                "notes": (
                    f"fetch completed but the local file failed a content sanity check: {mismatch}"
                ),
            }
        )

    return entry.model_copy(
        update={
            "status": ResourceStatus.FETCHED,
            "sha256": hash_file(dest),
            "size_bytes": dest.stat().st_size,
            "retrieved": _today_iso(),
            "notes": "",
        }
    )
