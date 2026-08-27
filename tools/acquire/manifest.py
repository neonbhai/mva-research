"""Rendering and loading ``knowledge/manifests/resources.yaml``.

This is the public record of what this tool has fetched: every entry in
``tools.acquire.catalog.KNOWN_RESOURCES``, in its current, disk-verified state
(``tools.acquire.survey.survey_all``), rendered deterministically so the file is
meaningful to commit and diff. It sits in ``knowledge/`` -- public, committed --
for the same reason ``knowledge/manifests/knowledge.yaml`` does: what was fetched,
from where, under what license, must be reviewable by anyone reading the repo.

Deliberately parallel to ``mva.annotation.local_tables``'s manifest handling
(``compute_manifest`` / ``render_manifest_yaml`` / ``load_manifest``), because this
is the same idea applied one step earlier in the pipeline: a generated,
hash-pinned index that is never hand-edited.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Final

import yaml

from mva.models.base import FrozenModel
from tools.acquire.errors import AcquisitionError
from tools.acquire.models import ResourceEntry

MANIFEST_VERSION: Final = 1
MANIFEST_GENERATOR: Final = "tools.acquire.manifest.render_resources_yaml"

#: What `path` on each entry is relative to. Deliberately NOT an absolute machine
#: path: the external resource root is per-machine ($MVA_RESOURCES or the default),
#: so committing one contributor's absolute path would make the manifest wrong on
#: every other machine. The manifest names the convention; `tools.acquire.fetch
#: .resolve_resource_root` resolves it at run time.
PATHS_RELATIVE_TO: Final = "resource_root"

MANIFEST_HEADER: Final = """\
# Versioned index of the public reference resources this project's ACQUISITION TOOL
# has fetched (or attempted to fetch) into the external resource root.
#
# GENERATED FILE - do not hand-edit. Regenerate with:
#
#   uv run python -m tools.acquire write-manifest
#
# `path` on each entry is relative to the resource root: $MVA_RESOURCES if set,
# otherwise ~/Contri/bio-hackathon/mva-resources (tools.acquire.fetch
# .resolve_resource_root). This file intentionally does NOT store that root as an
# absolute path -- it differs per machine, and the manifest must stay correct on
# everyone's.
#
# `sha256`, `size_bytes` and `retrieved` are null for any resource whose status is
# `not_fetched`: not yet started, still downloading, or present but failed a content
# sanity check (see `notes` on that entry). NEVER hand-fill these for a resource that
# is not actually, verifiably complete -- see tools/acquire/survey.py and
# knowledge/adapters/README.md, "Why no adapter here may touch the network (PRIV-05)".
#
# Every resource here is real public reference data: `synthetic: false` throughout.
# The synthetic demo tables live in knowledge/manifests/knowledge.yaml, not this file.
"""


class ResourceManifest(FrozenModel):
    """The versioned index of registered public reference resources.

    Strict (``extra="forbid"``): matches ``mva.annotation.local_tables.KnowledgeManifest``
    -- an unknown top-level key is a loud error, because this file makes an
    integrity claim.
    """

    manifest_version: int
    paths_relative_to: str
    generated_by: str
    resources: dict[str, ResourceEntry]


def render_resources_yaml(entries: tuple[ResourceEntry, ...]) -> str:
    """Render the committed manifest text. Deterministic: sorted keys, no wrapping.

    Duplicate names would silently collapse into one dict entry, which is exactly
    the kind of thing a generated file must never do quietly -- checked here even
    though ``tools.acquire.catalog`` already asserts uniqueness at import time,
    because this function must stay correct for any caller, not only that catalog.
    """
    names = [entry.name for entry in entries]
    if len(names) != len(set(names)):
        duplicates = sorted({name for name in names if names.count(name) > 1})
        msg = f"Cannot render a resources manifest with duplicate resource name(s): {duplicates}"
        raise AcquisitionError(msg)

    resources: dict[str, Any] = {
        entry.name: entry.model_dump(mode="json", exclude={"name"}) for entry in entries
    }
    body = yaml.safe_dump(
        {
            "manifest_version": MANIFEST_VERSION,
            "paths_relative_to": PATHS_RELATIVE_TO,
            "generated_by": MANIFEST_GENERATOR,
            "resources": resources,
        },
        sort_keys=True,
        default_flow_style=False,
        allow_unicode=False,
        width=10_000,
    )
    return f"{MANIFEST_HEADER}{body}"


def write_resources_manifest(manifest_path: Path, entries: tuple[ResourceEntry, ...]) -> None:
    """Render and write the manifest, guaranteeing a trailing newline."""
    text = render_resources_yaml(entries)
    if not text.endswith("\n"):
        text += "\n"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(text, encoding="utf-8")


def load_resources_manifest(manifest_path: Path) -> ResourceManifest:
    """Parse and validate the resources manifest."""
    if not manifest_path.is_file():
        msg = (
            f"Resources manifest {manifest_path.as_posix()!r} not found. Generate it with "
            "write_resources_manifest."
        )
        raise AcquisitionError(msg)
    try:
        raw: Any = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        msg = f"Resources manifest {manifest_path.name} is not valid YAML: {exc.__class__.__name__}"
        raise AcquisitionError(msg) from exc
    if not isinstance(raw, dict):
        msg = f"Resources manifest {manifest_path.name} must be a YAML mapping at the top level."
        raise AcquisitionError(msg)

    resources_raw = raw.get("resources")
    if isinstance(resources_raw, dict):
        for key, value in resources_raw.items():
            if isinstance(value, dict) and "name" not in value:
                value["name"] = key

    try:
        return ResourceManifest.model_validate(raw)
    except ValueError as exc:
        msg = f"Resources manifest {manifest_path.name} does not match the expected schema: {exc}"
        raise AcquisitionError(msg) from exc
