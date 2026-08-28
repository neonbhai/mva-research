"""A resumable digest cache for the registration pass.

Registering the reference set means computing a full sha256 over 202.8 GB. That
takes about 160 s here and cannot be made to take less, because the
whole point is that every byte was read. What it *can* be made to survive is
interruption: a crashed process, a machine that slept, a second resource added a
day later. Without a cache, each of those restarts the whole pass.

So the survey optionally consults a cache keyed on ``(path, size, mtime_ns)`` and
holding the digest that was computed for exactly those bytes.

**This cache is not an integrity check and is never allowed to become one.**
``size + mtime`` is a weak identity: any tool that rewrites a file and restores
its timestamp defeats it, and in-place corruption preserves both. That is fine
*here* and nowhere else, because the only thing the cache decides is whether to
re-read bytes during registration — the value it returns was itself produced by a
full read, and the manifest it lands in is the thing later checks compare
against. It is a build cache, not a verifier. ``--rehash`` ignores it entirely.

Kept out of ``src/mva`` deliberately: nothing on the run-time path may treat a
timestamp as evidence about content (see the module docstring of
``mva.resources``, and ADR 0020).
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any, Final

from mva.determinism import hash_file
from mva.resources import ResourceKind, spot_digest, spot_tree_digest, tree_digest
from tools.acquire.errors import AcquisitionError

#: Lives beside the data it describes, under the resource root, so it travels with
#: the resources and is never mistaken for a committed artifact. Dot-prefixed and
#: outside the repo, so nothing can commit it.
CACHE_FILENAME: Final = ".mva-digest-cache.json"

CACHE_VERSION: Final = 1


class DigestCache:
    """Full and sampled digests, memoised on ``(size, mtime_ns)`` within a registration.

    Not thread-safe and not intended to be: registration is a single pass, and the
    cost being avoided is I/O, not CPU contention.
    """

    def __init__(self, path: Path, *, enabled: bool = True) -> None:
        self._path = path
        self._enabled = enabled
        self._entries: dict[str, dict[str, Any]] = {}
        self._hits = 0
        self._misses = 0
        if enabled and path.is_file():
            self._entries = _load(path)

    @property
    def hits(self) -> int:
        return self._hits

    @property
    def misses(self) -> int:
        return self._misses

    def _key(self, target: Path) -> str:
        return target.as_posix()

    def _fingerprint(self, target: Path, kind: ResourceKind) -> tuple[int, int]:
        """``(size, mtime_ns)`` for a file; the summed/maximum pair for a directory."""
        if kind is ResourceKind.DIRECTORY:
            members = [p for p in target.rglob("*") if p.is_file() and not p.is_symlink()]
            stats = [p.stat() for p in members]
            return (
                sum(s.st_size for s in stats),
                max((s.st_mtime_ns for s in stats), default=0),
            )
        info = target.stat()
        return (info.st_size, info.st_mtime_ns)

    def digests(self, target: Path, kind: ResourceKind = ResourceKind.FILE) -> tuple[str, str]:
        """``(full_sha256, spot_sha256)`` for ``target``, reading it only if it has to."""
        size, mtime_ns = self._fingerprint(target, kind)
        key = self._key(target)
        cached = self._entries.get(key) if self._enabled else None
        if (
            cached is not None
            and cached.get("size") == size
            and cached.get("mtime_ns") == mtime_ns
            and isinstance(cached.get("sha256"), str)
            and isinstance(cached.get("spot_sha256"), str)
        ):
            self._hits += 1
            return (str(cached["sha256"]), str(cached["spot_sha256"]))

        self._misses += 1
        full_fn: Callable[[Path], str]
        spot_fn: Callable[[Path], str]
        full_fn, spot_fn = (
            (tree_digest, spot_tree_digest)
            if kind is ResourceKind.DIRECTORY
            else (hash_file, spot_digest)
        )
        # Spot first: it is a strict subset of the reads the full digest performs, so
        # doing it while the file is warm in the page cache costs approximately nothing.
        sampled = spot_fn(target)
        full = full_fn(target)
        self._entries[key] = {
            "size": size,
            "mtime_ns": mtime_ns,
            "sha256": full,
            "spot_sha256": sampled,
        }
        return (full, sampled)

    def save(self) -> None:
        """Persist the cache. Written atomically: a half-written cache is a wrong one."""
        if not self._enabled:
            return
        payload = {"cache_version": CACHE_VERSION, "entries": self._entries}
        self._path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self._path.with_suffix(".tmp")
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        temporary.replace(self._path)


def _load(path: Path) -> dict[str, dict[str, Any]]:
    """Read the cache, treating anything unreadable or stale as an empty cache.

    A corrupt cache costs one re-hashing pass; a corrupt cache that raises costs a
    failed registration. Discarding it silently is the right trade only
    because the cache is never the source of an integrity claim.
    """
    try:
        raw: Any = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    if not isinstance(raw, dict) or raw.get("cache_version") != CACHE_VERSION:
        return {}
    entries = raw.get("entries")
    if not isinstance(entries, dict):
        return {}
    return {str(key): value for key, value in entries.items() if isinstance(value, dict)}


def cache_path_for(resource_root: Path) -> Path:
    """Where the digest cache for ``resource_root`` lives."""
    if not resource_root.is_dir():
        msg = f"Resource root {resource_root.as_posix()!r} is not a directory."
        raise AcquisitionError(msg)
    return resource_root / CACHE_FILENAME
