"""The resource registry model: one hash-pinned public reference file.

This is the acquisition-side counterpart of ``mva.annotation.local_tables.KnowledgeTable``,
and deliberately mirrors it: same immutability, same ``extra="forbid"`` strictness, same
"synthetic must be written out, never inferred" philosophy (GP-20). It differs in one
load-bearing way -- ``sha256``, ``size_bytes`` and ``retrieved`` are *nullable*, because
unlike ``knowledge/manifests/knowledge.yaml`` (where every table is always fully present),
``knowledge/manifests/resources.yaml`` legitimately records resources that are declared but
not yet fetched, or mid-download. A null hash here means exactly what it says: nothing is
pinned yet. It must never be filled in with a guess.
"""

from __future__ import annotations

import re
from enum import StrEnum
from pathlib import Path
from typing import Self

from pydantic import Field, model_validator

from mva.models.base import FrozenModel
from tools.acquire.hosts import assert_allowed_host

_SHA256_RE: re.Pattern[str] = re.compile(r"^[0-9a-f]{64}$")


class ResourceStatus(StrEnum):
    """What this tool actually knows about a declared resource's local bytes."""

    FETCHED = "fetched"
    """A complete, format-sane, hash-verified local file backs this entry."""

    NOT_FETCHED = "not_fetched"
    """No trustworthy local file: never started, still growing, or failed a sanity
    check. ``sha256``/``size_bytes``/``retrieved`` are all ``None`` -- see ``notes``
    for why."""


class ResourceEntry(FrozenModel):
    """One declared public reference resource.

    Strict (``extra="forbid"``): an unrecognised key in a hand-edited manifest is a
    loud error, not a silently ignored field -- this file makes an integrity claim,
    same as ``knowledge/manifests/knowledge.yaml``.
    """

    name: str = Field(min_length=1)
    url: str = Field(description="The single https source URL this resource is fetched from.")
    version: str = Field(min_length=1, description="The real upstream release identifier.")
    path: str = Field(description="Relative to the resource root (e.g. 'clinvar/clinvar.vcf.gz').")
    license: str = Field(
        min_length=1, description="The resource's actual license/attribution terms."
    )
    description: str = Field(min_length=1)
    synthetic: bool = Field(
        default=False,
        description="Always False here: this registry exists to hold REAL public reference "
        "data. A synthetic table belongs in knowledge/manifests/knowledge.yaml, not here.",
    )
    status: ResourceStatus = ResourceStatus.NOT_FETCHED
    sha256: str | None = Field(
        default=None, description="sha256 of the local bytes. Null until fetched."
    )
    size_bytes: int | None = Field(default=None, ge=0)
    retrieved: str | None = Field(
        default=None, description="ISO date the local bytes were written."
    )
    notes: str = Field(
        default="",
        description="Why a not-yet-fetched resource has no hash (not started / still growing / "
        "failed a content check), or empty for a clean fetch.",
    )

    @model_validator(mode="after")
    def _url_is_allowlisted(self) -> Self:
        # Runs at construction time so a bad host can never even enter the registry --
        # not only at fetch time. See tools/acquire/hosts.py.
        assert_allowed_host(self.url)
        return self

    @model_validator(mode="after")
    def _path_is_relative_and_contained(self) -> Self:
        candidate = Path(self.path)
        if candidate.is_absolute():
            msg = (
                f"Resource {self.name!r} declares an absolute path {self.path!r}; must be "
                "relative to the resource root."
            )
            raise ValueError(msg)
        if ".." in candidate.parts:
            msg = (
                f"Resource {self.name!r} declares a path {self.path!r} escaping the resource "
                "root via '..'."
            )
            raise ValueError(msg)
        return self

    @model_validator(mode="after")
    def _sha256_shape(self) -> Self:
        if self.sha256 is not None and not _SHA256_RE.match(self.sha256):
            msg = (
                f"Resource {self.name!r} has a malformed sha256 {self.sha256!r} (expected 64 "
                "lowercase hex characters)."
            )
            raise ValueError(msg)
        return self

    @model_validator(mode="after")
    def _status_matches_fields(self) -> Self:
        """The whole point of ``sha256`` being nullable: it may ONLY be non-null when
        the status says this resource was actually, verifiably fetched. A manifest
        claiming a hash for a not-fetched resource -- or claiming FETCHED with no
        hash -- is not a data-entry slip, it is an integrity claim that cannot be
        checked, which is worse than an honest gap.
        """
        fetched_fields = (self.sha256, self.size_bytes, self.retrieved)
        if self.status is ResourceStatus.FETCHED:
            if any(field is None for field in fetched_fields):
                msg = (
                    f"Resource {self.name!r} is marked {ResourceStatus.FETCHED.value!r} but is "
                    "missing sha256, size_bytes or retrieved. A fetched resource must carry "
                    "all three."
                )
                raise ValueError(msg)
        elif any(field is not None for field in fetched_fields):
            msg = (
                f"Resource {self.name!r} is marked {self.status.value!r} but declares sha256, "
                "size_bytes or retrieved. Those fields may only be set once the resource is "
                "actually fetched and verified -- never write a hash for bytes that are absent, "
                "still growing, or failed a content sanity check."
            )
            raise ValueError(msg)
        return self
