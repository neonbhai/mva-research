"""Command-line interface for the public-reference acquisition tool.

Run it as ``uv run python -m tools.acquire <command>``. This is the separate,
public-only acquisition step ``knowledge/adapters/README.md`` describes: it runs
before and independently of any patient data, and every command here takes only
public-reference identifiers -- a resource name, a resource root path -- never a
variant coordinate, sample ID, or workspace path. See ``tools/acquire/models.py``
for the full field set a resource carries, and ``tools/acquire/hosts.py`` for the
allowlist that gates every fetch.
"""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Final

import typer

from mva.config import find_repo_root
from tools.acquire.catalog import KNOWN_RESOURCES
from tools.acquire.errors import AcquisitionError
from tools.acquire.fetch import fetch_resource, resolve_resource_root
from tools.acquire.manifest import load_resources_manifest, write_resources_manifest
from tools.acquire.models import ResourceStatus
from tools.acquire.survey import survey_all
from tools.acquire.verify import VerificationStatus, assert_verified, verify_all

app = typer.Typer(
    name="acquire",
    help="Fetch, hash-pin and record public reference resources. Never patient data.",
    no_args_is_help=True,
    add_completion=False,
)

#: Relative to the repo root.
DEFAULT_MANIFEST_RELATIVE: Final[Path] = Path("knowledge/manifests/resources.yaml")

ResourceRootOpt = Annotated[
    Path | None,
    typer.Option(
        "--resource-root",
        help=(
            "External resource root. Defaults to $MVA_RESOURCES, or "
            "~/Contri/bio-hackathon/mva-resources."
        ),
    ),
]
ManifestOpt = Annotated[
    Path | None,
    typer.Option("--manifest", help="Path to resources.yaml. Defaults to the repo's own copy."),
]


def _default_manifest_path() -> Path:
    return find_repo_root() / DEFAULT_MANIFEST_RELATIVE


@app.command()
def status(resource_root: ResourceRootOpt = None) -> None:
    """Show the fetched/not-fetched state of every registered resource, from disk."""
    root = resolve_resource_root(resource_root)
    surveyed = survey_all(root, KNOWN_RESOURCES)
    fetched = sum(1 for entry in surveyed if entry.status is ResourceStatus.FETCHED)
    typer.echo(f"resource root: {root.as_posix()}")
    typer.echo(f"{fetched}/{len(surveyed)} resources fetched\n")
    for entry in surveyed:
        marker = "FETCHED    " if entry.status is ResourceStatus.FETCHED else "NOT_FETCHED"
        detail = entry.sha256[:12] if entry.sha256 else (entry.notes or "-")
        typer.echo(f"  {marker}  {entry.name:<32} {detail}")


@app.command()
def fetch(
    name: Annotated[
        str | None, typer.Option("--name", help="Fetch only this resource; default: all of them.")
    ] = None,
    resource_root: ResourceRootOpt = None,
) -> None:
    """Download one (or every) registered resource, resuming a partial file if present."""
    root = resolve_resource_root(resource_root)
    targets = [e for e in KNOWN_RESOURCES if name is None or e.name == name]
    if name is not None and not targets:
        typer.echo(f"no such resource: {name!r}", err=True)
        raise typer.Exit(code=1)

    failures = 0
    for entry in targets:
        typer.echo(f"fetching {entry.name} ...")
        try:
            result = fetch_resource(entry, root)
        except AcquisitionError as exc:
            typer.echo(f"  FAILED: {exc}", err=True)
            failures += 1
            continue
        if result.status is ResourceStatus.FETCHED:
            typer.echo(f"  OK  sha256={result.sha256}")
        else:
            typer.echo(f"  NOT COMPLETE: {result.notes}", err=True)
            failures += 1

    if failures:
        raise typer.Exit(code=1)


@app.command(name="write-manifest")
def write_manifest_command(
    manifest: ManifestOpt = None,
    resource_root: ResourceRootOpt = None,
) -> None:
    """Survey the resource root and (re)write the committed resources.yaml."""
    root = resolve_resource_root(resource_root)
    manifest_path = manifest or _default_manifest_path()
    surveyed = survey_all(root, KNOWN_RESOURCES)
    write_resources_manifest(manifest_path, surveyed)
    fetched = sum(1 for entry in surveyed if entry.status is ResourceStatus.FETCHED)
    typer.echo(f"wrote {manifest_path.as_posix()} ({fetched}/{len(surveyed)} resources fetched)")


@app.command()
def verify(
    manifest: ManifestOpt = None,
    resource_root: ResourceRootOpt = None,
    strict: Annotated[
        bool, typer.Option("--strict", help="Exit non-zero if any resource is MISSING or MISMATCH.")
    ] = True,
) -> None:
    """Check the committed manifest's hashes against the resource root on disk."""
    root = resolve_resource_root(resource_root)
    manifest_path = manifest or _default_manifest_path()
    parsed = load_resources_manifest(manifest_path)
    entries = tuple(parsed.resources.values())

    if strict:
        try:
            results = assert_verified(root, entries)
        except AcquisitionError as exc:
            typer.echo(str(exc), err=True)
            raise typer.Exit(code=1) from exc
    else:
        results = verify_all(root, entries)

    for result in results:
        typer.echo(f"  {result.status.value.upper():<9} {result.name:<32} {result.message}")

    broken = (VerificationStatus.MISSING, VerificationStatus.MISMATCH)
    failed = sum(1 for r in results if r.status in broken)
    if failed and not strict:
        raise typer.Exit(code=1)
