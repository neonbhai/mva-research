"""The acquisition half of the resource registry model.

The *shape* of a registered resource lives in :mod:`mva.resources`, because the
pipeline has to read it at run time and ``src/mva`` may not import this package
(it opens sockets; the patient-data path is structurally forbidden a network
client — PRIV-05). :class:`ResourceEntry` here **subclasses**
:class:`mva.resources.ReferenceResource` and adds the two fields only the
acquisition side has any business knowing: the upstream ``url`` and the
``license``.

Subclassing rather than restating the fields is the point. Two hand-maintained
copies of a manifest schema drift, and the drift shows up as a resource that
writes fine and reads as ``extra="ignore"``-shaped silence. One class, one set of
validators, and ``extra="forbid"`` restored here so that a typo in a hand-edited
manifest is caught where it can be fixed — at the writer — rather than at the
reader, where it can only be reported.

This is the counterpart of ``mva.annotation.local_tables.KnowledgeTable``, and
differs from it in one load-bearing way: ``sha256``, ``size_bytes`` and
``retrieved`` are *nullable*, because ``knowledge/manifests/resources.yaml``
legitimately records resources that are declared but not yet fetched, or
mid-download. A null hash means exactly what it says: nothing is pinned yet. It
must never be filled in with a guess.
"""

from __future__ import annotations

from typing import Self

from pydantic import ConfigDict, Field, model_validator

from mva.resources import (
    FormatCheck,
    IntegrityRecord,
    ReferenceResource,
    ResourceKind,
    ResourceStatus,
)
from tools.acquire.hosts import assert_allowed_host

__all__ = [
    "FormatCheck",
    "IntegrityRecord",
    "ResourceEntry",
    "ResourceKind",
    "ResourceStatus",
]


class ResourceEntry(ReferenceResource):
    """One declared public reference resource, with its provenance.

    Strict (``extra="forbid"``): an unrecognised key in a hand-edited manifest is a
    loud error, not a silently ignored field — this file makes an integrity claim,
    same as ``knowledge/manifests/knowledge.yaml``.
    """

    model_config = ConfigDict(
        frozen=True,
        extra="forbid",
        str_strip_whitespace=True,
        validate_default=True,
        hide_input_in_errors=True,
    )

    url: str = Field(description="The single https source URL this resource is fetched from.")
    license: str = Field(
        min_length=1, description="The resource's actual license/attribution terms."
    )

    @model_validator(mode="after")
    def _url_is_allowlisted(self) -> Self:
        # Runs at construction time so a bad host can never even enter the registry --
        # not only at fetch time. See tools/acquire/hosts.py.
        assert_allowed_host(self.url)
        return self
