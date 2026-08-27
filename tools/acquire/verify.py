"""Verifying a committed ``resources.yaml`` against what is actually on disk.

This is the acquisition-side counterpart of
``mva.annotation.local_tables.verify_manifest``: given a manifest, confirm every
resource it claims to have pinned still matches. It differs from that function in
one way that matters -- ``knowledge/manifests/knowledge.yaml`` never contains an
unfetched table, so ``verify_manifest`` only ever needs MISSING/MISMATCH/OK. Here a
resource can legitimately be declared and not yet fetched (``sha256`` is ``None`` by
design; see ``ResourceEntry``), so verification needs a fourth outcome that means
"nothing to check, and that's expected" rather than "the pin failed".
"""

from __future__ import annotations

from collections.abc import Iterable
from enum import StrEnum
from pathlib import Path

from mva.determinism import hash_file
from mva.models.base import FrozenModel
from tools.acquire.errors import ResourceVerificationError
from tools.acquire.models import ResourceEntry, ResourceStatus


class VerificationStatus(StrEnum):
    OK = "ok"
    """The local file exists and its sha256 matches the manifest."""

    MISSING = "missing"
    """The manifest declares this resource FETCHED, but no local file exists."""

    MISMATCH = "mismatch"
    """The manifest declares this resource FETCHED, but the local sha256 disagrees."""

    PENDING = "pending"
    """The manifest declares this resource NOT_FETCHED. There is nothing to verify --
    this is the expected, honest state for a resource that has not been downloaded
    (or completed downloading) yet."""


class VerificationResult(FrozenModel):
    """The outcome of checking one resource. Never carries file bytes (PRIV-09)."""

    name: str
    path: str
    status: VerificationStatus
    message: str


def verify_resource(resource_root: Path, entry: ResourceEntry) -> VerificationResult:
    """Check one declared resource's local bytes against its manifest pin.

    Names the resource and its declared path in every message; never reads the
    mismatching file's content into the message (PRIV-09) -- a sha256 hex digest is
    a one-way summary, not "contents", and is safe to include in full.
    """
    if entry.status is ResourceStatus.NOT_FETCHED or entry.sha256 is None:
        return VerificationResult(
            name=entry.name,
            path=entry.path,
            status=VerificationStatus.PENDING,
            message=f"resource {entry.name!r} is not yet fetched; nothing to verify",
        )

    local_path = resource_root / entry.path
    if not local_path.is_file():
        return VerificationResult(
            name=entry.name,
            path=entry.path,
            status=VerificationStatus.MISSING,
            message=(
                f"resource {entry.name!r} is pinned in the manifest but missing from disk "
                f"(expected {entry.path!r} under {resource_root.as_posix()!r})"
            ),
        )

    actual = hash_file(local_path)
    if actual != entry.sha256:
        return VerificationResult(
            name=entry.name,
            path=entry.path,
            status=VerificationStatus.MISMATCH,
            message=(
                f"resource {entry.name!r} ({entry.path!r}) failed its manifest integrity "
                f"check: expected sha256 {entry.sha256}, found {actual}. The file changed "
                "without the manifest being regenerated."
            ),
        )

    return VerificationResult(
        name=entry.name,
        path=entry.path,
        status=VerificationStatus.OK,
        message=f"resource {entry.name!r} matches its manifest pin",
    )


def verify_all(
    resource_root: Path, entries: Iterable[ResourceEntry]
) -> tuple[VerificationResult, ...]:
    """:func:`verify_resource` over every declared entry, in declaration order."""
    return tuple(verify_resource(resource_root, entry) for entry in entries)


def assert_verified(
    resource_root: Path, entries: Iterable[ResourceEntry]
) -> tuple[VerificationResult, ...]:
    """Verify every entry, raising if any is MISSING or MISMATCH.

    ``PENDING`` entries never fail this: a resource that is honestly declared as not
    yet fetched is not an integrity violation. Only a claimed-and-broken pin is.
    """
    results = verify_all(resource_root, entries)
    broken = (VerificationStatus.MISSING, VerificationStatus.MISMATCH)
    failures = [r for r in results if r.status in broken]
    if failures:
        names = ", ".join(f"{r.name!r} ({r.status.value})" for r in failures)
        msg = f"Resource manifest verification failed for: {names}."
        raise ResourceVerificationError(msg)
    return results
