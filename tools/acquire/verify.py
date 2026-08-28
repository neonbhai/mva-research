"""Verifying a committed ``resources.yaml`` against what is actually on disk.

The checking itself lives in :mod:`mva.resources`, not here, and that placement is
deliberate. The pipeline has to verify these pins at run time — it is the party
with something to lose from annotating against the wrong bytes — and it cannot
import this package, which opens sockets (PRIV-05). Two verifiers would then have
to agree forever about what a mismatch is, and would not. So there is one, in the
layer both sides can see, and this module is a thin adapter that maps its outcomes
onto the vocabulary the acquisition CLI already speaks.

The vocabularies differ in exactly one place, which is why the mapping is not the
identity: ``knowledge/manifests/knowledge.yaml`` never contains an unfetched
table, so ``mva.annotation.local_tables.verify_manifest`` only ever needs
MISSING/MISMATCH/OK. Here a resource can legitimately be declared and not yet
fetched (``sha256`` is ``None`` by design; see :class:`ResourceEntry`), so
verification needs a fourth outcome meaning "nothing to check, and that's
expected" rather than "the pin failed".

Defaults to :attr:`~mva.resources.IntegrityMode.FULL`. This is the registration
side: it re-reads every byte, because the thing it is deciding is whether to keep
believing a hash that is about to be committed.
"""

from __future__ import annotations

from collections.abc import Iterable
from enum import StrEnum
from pathlib import Path

from mva.models.base import FrozenModel
from mva.resources import IntegrityMode, ResourceCheck
from mva.resources import verify_resource as _verify
from tools.acquire.errors import ResourceVerificationError
from tools.acquire.models import ResourceEntry


class VerificationStatus(StrEnum):
    OK = "ok"
    """The local artifact exists and its digest matches the manifest."""

    MISSING = "missing"
    """The manifest declares this resource FETCHED, but nothing is on disk."""

    MISMATCH = "mismatch"
    """The manifest declares this resource FETCHED, but the local bytes disagree —
    wrong size, wrong digest, or a spot digest recorded under a sampling plan this
    build cannot evaluate."""

    PENDING = "pending"
    """The manifest declares this resource NOT_FETCHED. There is nothing to verify --
    this is the expected, honest state for a resource that has not been downloaded
    (or completed downloading) yet."""


#: Every outcome the shared checker can produce maps to exactly one status here.
#: Exhaustive by construction rather than by ``else``: a new ``ResourceCheck``
#: member added upstream raises a KeyError at the mapping site, which is a loud
#: failure at the place that needs updating — instead of silently falling into
#: whichever bucket the fallback happened to name.
_STATUS_BY_CHECK: dict[ResourceCheck, VerificationStatus] = {
    ResourceCheck.OK: VerificationStatus.OK,
    ResourceCheck.UNPINNED: VerificationStatus.PENDING,
    ResourceCheck.MISSING: VerificationStatus.MISSING,
    ResourceCheck.SIZE_MISMATCH: VerificationStatus.MISMATCH,
    ResourceCheck.DIGEST_MISMATCH: VerificationStatus.MISMATCH,
    ResourceCheck.PLAN_UNKNOWN: VerificationStatus.MISMATCH,
}


class VerificationResult(FrozenModel):
    """The outcome of checking one resource. Never carries file bytes (PRIV-09)."""

    name: str
    path: str
    status: VerificationStatus
    message: str


def verify_resource(
    resource_root: Path,
    entry: ResourceEntry,
    *,
    mode: IntegrityMode = IntegrityMode.FULL,
) -> VerificationResult:
    """Check one declared resource's local bytes against its manifest pin.

    Names the resource and its declared path in every message; never reads the
    mismatching file's content into the message (PRIV-09) — a sha256 hex digest is
    a one-way summary, not "contents", and is safe to include in full.
    """
    outcome = _verify(resource_root, entry, mode=mode)
    return VerificationResult(
        name=outcome.name,
        path=outcome.path,
        status=_STATUS_BY_CHECK[outcome.check],
        message=outcome.message,
    )


def verify_all(
    resource_root: Path,
    entries: Iterable[ResourceEntry],
    *,
    mode: IntegrityMode = IntegrityMode.FULL,
) -> tuple[VerificationResult, ...]:
    """:func:`verify_resource` over every declared entry, in declaration order."""
    return tuple(verify_resource(resource_root, entry, mode=mode) for entry in entries)


def assert_verified(
    resource_root: Path,
    entries: Iterable[ResourceEntry],
    *,
    mode: IntegrityMode = IntegrityMode.FULL,
) -> tuple[VerificationResult, ...]:
    """Verify every entry, raising if any is MISSING or MISMATCH.

    ``PENDING`` entries never fail this: a resource that is honestly declared as not
    yet fetched is not an integrity violation. Only a claimed-and-broken pin is.
    """
    results = verify_all(resource_root, entries, mode=mode)
    broken = (VerificationStatus.MISSING, VerificationStatus.MISMATCH)
    failures = [r for r in results if r.status in broken]
    if failures:
        names = ", ".join(f"{r.name!r} ({r.status.value})" for r in failures)
        detail = "\n".join(f"  {r.message}" for r in failures)
        msg = f"Resource manifest verification failed for: {names}.\n{detail}"
        raise ResourceVerificationError(msg)
    return results
