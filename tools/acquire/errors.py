"""Typed exceptions for the public-reference acquisition tool.

Kept separate from ``mva.errors`` deliberately: this tool is not part of the
``mva`` package and must stay usable (and testable) even while ``src/mva`` is
being edited concurrently by other work. See ``knowledge/adapters/README.md``,
"Why no adapter here may touch the network (PRIV-05)", for why this tool exists
as a standalone step rather than a module under ``mva.annotation``.
"""

from __future__ import annotations


class AcquisitionError(Exception):
    """Base for every error this tool raises."""


class DisallowedHostError(AcquisitionError):
    """A URL's host is not on the public-reference allowlist.

    Raised before any connection is opened -- the whole point of the allowlist
    is that a bad host never gets as far as a socket.
    """


class ResourceRootError(AcquisitionError):
    """The resolved external resource root violates a safety boundary.

    Raised when the configured root resolves inside this repository. Large
    reference binaries (a 190MB ClinVar VCF, multi-GB gnomAD sites files) must
    never land somewhere a stray ``git add -A`` could commit them.
    """


class ResourceVerificationError(AcquisitionError):
    """One or more declared resources failed manifest verification.

    The message names the offending resource and path -- never file contents
    (PRIV-09). A sha256 digest is not "contents": it is a one-way summary, and
    is safe to print in full.
    """


class ResourceFetchError(AcquisitionError):
    """A download failed after exhausting retries.

    Wraps the underlying transport error rather than re-raising its class directly:
    some ``OSError`` subclasses (``urllib.error.HTTPError`` in particular) require
    constructor arguments this tool has no reason to fabricate.
    """
