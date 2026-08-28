"""The out-of-repo reference releases: where they live, and whether they are intact.

Three responsibilities, in the order a run needs them:

1. **Resolve the resource root.** gnomAD v4.1 exomes alone is 184.8 GB. It cannot
   live in the repository (ADR 0006's reasoning about size, not privacy: the
   patient workspace is kept out because it is sensitive; these are kept out
   because a stray ``git add -A`` on 185 GB of public reference data is
   unrecoverable in a different way). The root is configured **explicitly** and is
   never guessed, exactly like the workspace — see :func:`resolve_resource_root`.

2. **Read the pinned manifest.** ``knowledge/manifests/resources.yaml`` records
   every release by ``sha256``, with its true upstream version string. It is
   written by ``tools/acquire`` — the public-only, network-facing acquisition step
   that is structurally separate from anything that sees patient data (PRIV-05) —
   and only *read* here.

3. **Verify, affordably, and fail closed.** This is the interesting part; see
   below and ADR 0020.

Why verification could not simply be "re-hash everything"
--------------------------------------------------------

Full sha256 over the registered set is 202.8 GB of reads. Measured here on
2026-08-28, end to end over all 65 registered resources: **159 s**. Doing that at
the start of every pipeline run would cost more than the run itself. Doing it
*never* means the pins are decoration.

So there are two tiers, and the manifest records which one was applied and when:

* :attr:`IntegrityMode.FULL` — every byte. What registration does, and what a
  release check should do. Correct; 159 s for the registered set.
* :attr:`IntegrityMode.SPOT` — the default at run time. Exact byte size, plus a
  sha256 over a **fixed, published, reproducible sample** of the file (head, tail
  and eight evenly-spaced interior windows: :data:`SPOT_PLAN`). At most 24 MiB per
  file regardless of size: 1.27 GB and **1.7 s** for all 65 registered resources,
  which is 0.63% of the bytes and 94x cheaper.

The honest statement of what SPOT proves — written here because a reader of the
manifest deserves it in the same place as the value:

* It **does** catch a truncated or still-downloading file, a file replaced by an
  HTML error page, a file swapped for a different release, and corruption that
  lands in a sampled window. Those are the failures that actually happen.
* It **does not** prove whole-file identity. For the 18.79 GB chr1 shard it reads
  24 MiB — 0.134% of the file. An adversary who can write to the resource root can
  defeat it deliberately,
  and so can bit-rot that misses every window. It is a checksum over a sample,
  not a proof, and it must never be described as one.
* ``size + mtime`` alone — the obvious cheap check — is **not** honest and is not
  offered as a mode. mtime is trivially preserved by any tool that rewrites a
  file, and size is preserved by any in-place corruption. It appears here only as
  a cache key inside the registration tool, where a full digest is computed
  anyway.

Every failure message names the resource and its path and never its contents
(PRIV-09). A sha256 digest is a one-way summary, not content, and is printed in
full so a mismatch is diagnosable.
"""

from __future__ import annotations

import hashlib
import os
from collections.abc import Iterable, Sequence
from enum import StrEnum
from pathlib import Path
from typing import Any, Final, Self

import yaml
from pydantic import ConfigDict, Field, field_validator, model_validator

from mva.config import CaseConfig, Workspace, find_repo_root, path_is_within
from mva.determinism import hash_file
from mva.errors import MvaError
from mva.models.base import FrozenModel

# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------


class ResourceError(MvaError):
    """A configured reference resource is unusable."""


class ResourceRootError(ResourceError):
    """The resource root is unset, missing, or violates a containment boundary."""


class ResourceIntegrityError(ResourceError):
    """A registered resource failed its manifest pin.

    Names the resource and its declared path, and quotes both digests. Never the
    file's contents (PRIV-09).
    """


# ---------------------------------------------------------------------------
# The spot-check sampling plan (GP-30: deterministic, total, no RNG)
# ---------------------------------------------------------------------------

#: Bytes read from the start of a file.
SPOT_HEAD_BYTES: Final = 8 * 1024 * 1024

#: Bytes read from the end of a file. The tail is not optional: a resumable
#: download that stalled leaves a file whose head is perfect, and for BGZF the
#: 28-byte end-of-file block lives here.
SPOT_TAIL_BYTES: Final = 8 * 1024 * 1024

#: Interior sample windows, evenly spaced. Eight rather than one because a single
#: midpoint window is a fixed target; evenly-spaced windows spread the sampled
#: region across the whole file for the same read cost.
SPOT_WINDOW_COUNT: Final = 8
SPOT_WINDOW_BYTES: Final = 1024 * 1024

#: The plan string recorded in the manifest next to every ``spot_sha256``.
#: It is part of the value: a digest computed under a different plan is a
#: different number, so changing any constant above REQUIRES bumping this and
#: re-registering. Verification refuses a digest whose plan it does not recognise
#: rather than comparing two incomparable numbers.
SPOT_PLAN: Final = "spot-v1:head=8MiB,tail=8MiB,windows=8x1MiB,size-bound"

#: Below this size the "sample" would cover the whole file anyway, so the plan
#: degenerates to a single window over all of it — and for such a file SPOT and
#: FULL are the same check. Stated rather than left implicit, because "the cheap
#: check happens to be exhaustive here" is a fact worth being able to rely on.
SPOT_EXHAUSTIVE_BELOW: Final = (
    SPOT_HEAD_BYTES + SPOT_TAIL_BYTES + (SPOT_WINDOW_COUNT * SPOT_WINDOW_BYTES)
)

_HASH_BLOCK: Final = 1 << 20


def spot_windows(size: int) -> tuple[tuple[int, int], ...]:
    """The ``(offset, length)`` byte ranges :func:`spot_digest` reads, for a file of ``size``.

    Deterministic and total: the same size always yields the same ranges, in the
    same order, on every platform. Exposed (rather than kept private inside
    :func:`spot_digest`) so the sampling claim in the manifest is auditable — a
    reviewer can print the ranges and check the coverage arithmetic themselves.
    """
    if size <= 0:
        return ()
    if size <= SPOT_EXHAUSTIVE_BELOW:
        return ((0, size),)

    interior_start = SPOT_HEAD_BYTES
    interior_end = size - SPOT_TAIL_BYTES - SPOT_WINDOW_BYTES
    span = interior_end - interior_start
    windows: list[tuple[int, int]] = [(0, SPOT_HEAD_BYTES)]
    for index in range(SPOT_WINDOW_COUNT):
        # Integer arithmetic only. Float division here would make the plan
        # depend on rounding mode, which is precisely the kind of thing GP-30
        # exists to keep out of a value that gets committed to a manifest.
        offset = interior_start + (span * index) // (SPOT_WINDOW_COUNT - 1)
        windows.append((offset, SPOT_WINDOW_BYTES))
    windows.append((size - SPOT_TAIL_BYTES, SPOT_TAIL_BYTES))
    return tuple(windows)


def spot_digest(path: Path) -> str:
    """sha256 over the sampled bytes of ``path``, bound to its exact size.

    The file's size is mixed into the digest before any content, so a file
    truncated or extended between the sampled windows cannot produce a matching
    value even if every sampled byte is unchanged. The offset and length of each
    window are mixed in too, which makes the digest meaningless to compare across
    plans — deliberately, since :data:`SPOT_PLAN` is what says which plan produced
    it.

    Reads at most :data:`SPOT_EXHAUSTIVE_BELOW` bytes (24 MiB) however large the
    file is. This is a sample, not a proof — see the module docstring.
    """
    size = path.stat().st_size
    digest = hashlib.sha256()
    digest.update(SPOT_PLAN.encode("ascii"))
    digest.update(b"\nsize=")
    digest.update(str(size).encode("ascii"))
    digest.update(b"\n")
    with path.open("rb") as handle:
        for offset, length in spot_windows(size):
            digest.update(f"@{offset}+{length}\n".encode("ascii"))
            handle.seek(offset)
            remaining = length
            while remaining > 0:
                chunk = handle.read(min(_HASH_BLOCK, remaining))
                if not chunk:
                    break
                digest.update(chunk)
                remaining -= len(chunk)
    return digest.hexdigest()


def tree_digest(root: Path) -> str:
    """sha256 over a directory's *shape and contents*, for a multi-file resource.

    A SnpEff genome database is 28 files in one directory; pinning only
    ``snpEffectPredictor.bin`` would leave the 26 ``sequence.*.bin`` files that
    supply every codon and HGVS string unpinned. This digests the sorted list of
    ``(relative path, size, sha256)`` triples, so adding, removing, renaming or
    editing any member changes the value.

    Sorted by POSIX relative path so the result does not depend on directory
    iteration order (GP-30). Symlinks and non-regular files are skipped and their
    absence is part of the digest only insofar as they were never included; a
    resource that needs them is not a candidate for this function.
    """
    digest = hashlib.sha256()
    digest.update(b"tree-v1\n")
    members = sorted(
        (p for p in root.rglob("*") if p.is_file() and not p.is_symlink()),
        key=lambda p: p.relative_to(root).as_posix(),
    )
    for member in members:
        relative = member.relative_to(root).as_posix()
        digest.update(f"{relative}\0{member.stat().st_size}\0{hash_file(member)}\n".encode())
    return digest.hexdigest()


def spot_tree_digest(root: Path) -> str:
    """:func:`tree_digest`'s cheap counterpart: sizes exhaustively, contents sampled."""
    digest = hashlib.sha256()
    digest.update(b"tree-spot-v1\n")
    digest.update(SPOT_PLAN.encode("ascii"))
    digest.update(b"\n")
    members = sorted(
        (p for p in root.rglob("*") if p.is_file() and not p.is_symlink()),
        key=lambda p: p.relative_to(root).as_posix(),
    )
    for member in members:
        relative = member.relative_to(root).as_posix()
        digest.update(f"{relative}\0{member.stat().st_size}\0{spot_digest(member)}\n".encode())
    return digest.hexdigest()


# ---------------------------------------------------------------------------
# The manifest read model
# ---------------------------------------------------------------------------


class ResourceKind(StrEnum):
    """Whether a registered resource is one file or a whole directory."""

    FILE = "file"
    DIRECTORY = "directory"
    """Digested as a tree (:func:`tree_digest`). Used for the SnpEff genome
    database, which is 28 interdependent files that are only meaningful together."""


class ResourceStatus(StrEnum):
    """What is actually known about a declared resource's local bytes."""

    FETCHED = "fetched"
    """A complete, format-checked, hash-pinned local artifact backs this entry."""

    NOT_FETCHED = "not_fetched"
    """Nothing trustworthy locally: never started, still growing, or failed a
    check. ``sha256`` / ``size_bytes`` / ``retrieved`` are all ``None`` by
    construction — see ``notes``."""


class FormatCheck(StrEnum):
    """What was proven about a resource's *format*, beyond its bytes hashing right.

    A sha256 pins whatever was downloaded. It does not notice that what was
    downloaded is an HTML error page — which happened here for real: a stale
    Gene2Phenotype bulk-download URL returned a JavaScript app shell with HTTP
    200, and ``file -b`` said "HTML document text". A hash of that page is a
    perfectly valid hash of the wrong thing. These values record that something
    actually opened the file as the format its name claims.
    """

    NOT_CHECKED = "not_checked"
    """No format check has been run. The honest default; never a passing state."""

    BGZF_TABIX_INDEXED = "bgzf_tabix_indexed"
    """A BGZF-compressed VCF whose tabix index opens and lists contigs, and whose
    first data record parses."""

    TABIX_INDEX = "tabix_index"
    """A tabix index: gzip member carrying the ``TBI\\1`` magic, loadable."""

    FASTA_FAIDX_CONSISTENT = "fasta_faidx_consistent"
    """A FASTA whose ``.fai`` was cross-checked against the file: every sampled
    contig's declared offset really lands on that contig's sequence."""

    FASTA_INDEX = "fasta_index"
    """A ``.fai``: every line parses as the five-column faidx record."""

    GZIP_TEXT = "gzip_text"
    """A gzip member that decompresses and whose leading lines match the expected
    shape (GTF columns, a TSV header, an OBO stanza)."""

    PLAIN_TEXT = "plain_text"
    """An uncompressed text table whose header row parses."""

    JAVA_ARCHIVE = "java_archive"
    """A zip/jar whose central directory opens and lists entries."""

    SNPEFF_DATABASE = "snpeff_database"
    """A SnpEff genome data directory: ``snpEffectPredictor.bin`` present and
    non-empty, and the per-contig ``sequence.*.bin`` files present."""

    OPAQUE_BINARY = "opaque_binary"
    """Non-empty binary with no independently checkable structure. Recorded
    honestly rather than dressed up as a pass."""


class IndexCheck(StrEnum):
    """Whether a companion index was proven to describe *this* data file.

    A tabix index and its VCF are two separate downloads. They can end up
    mismatched in ways that nothing downstream notices: the index resumes from an
    older attempt, the data file is re-fetched while the index is not, or an
    interrupted transfer leaves a complete index beside a shorter file. htslib
    warns about this when the index's **mtime precedes the data file's** — which
    is precisely the case for all 25 gnomAD shards here, because each ``.tbi`` was
    fetched before its multi-gigabyte ``.bgz`` finished.

    That warning is not evidence. mtime says nothing about content (ADR 0020),
    and in this case it is a false alarm produced by ordinary download ordering.
    But the *failure it gestures at* is real, so it is checked properly instead of
    either trusting the timestamps or ignoring the warning: the index is used to
    seek into the data file, and the records it returns must land where it said
    they would.
    """

    NOT_APPLICABLE = "not_applicable"
    """This resource has no companion index — or it IS one."""

    NOT_CHECKED = "not_checked"
    """An index exists but nothing opened it (pysam absent). Never a pass."""

    CONSISTENT = "consistent"
    """Random access through the index returned records at the coordinates the
    index claims. The index describes this data file."""

    STALE = "stale"
    """The index does not describe this data file. Fails closed: a resource in this
    state is recorded ``not_fetched`` rather than pinned."""


class IntegrityRecord(FrozenModel):
    """Exactly what was verified about one resource, and when.

    Written by the registration tool, read by the run-time check. Its existence is
    the point of ADR 0020: a manifest that pins a hash but does not say *when the
    bytes were last actually looked at, and how hard* invites the reader to assume
    more than was done.
    """

    verified_at: str = Field(
        min_length=1,
        description="ISO date the full digest and the format check below were performed.",
    )
    spot_plan: str = Field(
        min_length=1,
        description=(
            "The sampling plan `spot_sha256` was computed under (mva.resources.SPOT_PLAN). "
            "A digest computed under a different plan is a different number and is refused, "
            "not compared."
        ),
    )
    spot_sha256: str = Field(description="Digest over the sampled bytes. A checksum, not a proof.")
    format_check: FormatCheck = FormatCheck.NOT_CHECKED
    format_detail: str = Field(
        default="",
        description=(
            "What the format check actually observed — contig counts, index entry counts, "
            "the first parsed record's coordinate. Public reference data only; no patient "
            "coordinate can reach this field."
        ),
    )
    index_check: IndexCheck = Field(
        default=IndexCheck.NOT_APPLICABLE,
        description=(
            "Whether a companion index (.tbi/.csi/.fai) was proven to describe THIS data "
            "file, by seeking through it rather than by comparing timestamps."
        ),
    )
    index_detail: str = Field(
        default="",
        description="What the index cross-check observed, including any mtime skew and why it "
        "was not treated as evidence.",
    )

    @field_validator("spot_sha256")
    @classmethod
    def _hex64(cls, value: str) -> str:
        if len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
            msg = "spot_sha256 must be 64 lowercase hex characters."
            raise ValueError(msg)
        return value


class ReferenceResource(FrozenModel):
    """One registered public reference release, as the pipeline reads it.

    ``extra="ignore"`` here, uniquely among this project's models, and the reason
    is a real one rather than laziness. This is the *consumer* half of a
    producer/consumer manifest contract: ``tools/acquire`` writes acquisition
    metadata the run time deliberately does not model (the upstream URL, the
    license text) and its own :class:`~tools.acquire.models.ResourceEntry`
    **subclasses this model** with ``extra="forbid"`` restored. So a typo in a
    manifest key is still a loud error — it is caught at the place that writes it,
    where it can be fixed, rather than at the place that reads it, where it can
    only be reported. The subclass relationship is asserted by a test, so the two
    halves cannot drift.
    """

    model_config = ConfigDict(
        frozen=True,
        extra="ignore",
        str_strip_whitespace=True,
        validate_default=True,
        hide_input_in_errors=True,
    )

    name: str = Field(min_length=1)
    version: str = Field(
        min_length=1,
        description=(
            "The TRUE upstream release identifier — 'v4.1', a ClinVar weekly date, "
            "'GRCh38.115'. Stamped onto every EvidenceItem derived from this resource, so "
            "a package version or an invented string here makes a claim unreproducible."
        ),
    )
    path: str = Field(description="Relative to the resource root. Never absolute.")
    description: str = Field(min_length=1)
    kind: ResourceKind = ResourceKind.FILE
    synthetic: bool = Field(
        default=False,
        description=(
            "False for every entry in this registry: it holds real public reference data. "
            "Synthetic demo tables live in knowledge/manifests/knowledge.yaml (GP-20)."
        ),
    )
    status: ResourceStatus = ResourceStatus.NOT_FETCHED
    sha256: str | None = Field(default=None, description="Full digest. Null until fetched.")
    size_bytes: int | None = Field(default=None, ge=0)
    retrieved: str | None = Field(
        default=None, description="ISO date the local bytes were written."
    )
    integrity: IntegrityRecord | None = Field(
        default=None,
        description="What was verified and when. Null for a resource that is not fetched.",
    )
    derived_from: str | None = Field(
        default=None,
        description=(
            "For an artifact produced locally from another registered resource — the "
            "decompressed FASTA, its faidx index, an unzipped database. Names the parent "
            "resource and the transformation, so the URL on this entry is understood as "
            "the origin of the bytes rather than a claim that this exact file was served."
        ),
    )
    notes: str = Field(
        default="",
        description="Why an unfetched resource has no hash, or empty for a clean fetch.",
    )

    @field_validator("path")
    @classmethod
    def _path_is_relative_and_contained(cls, value: str) -> str:
        candidate = Path(value)
        if candidate.is_absolute():
            msg = (
                f"Resource declares an absolute path {value!r}; resource paths must be "
                "relative to the resource root."
            )
            raise ValueError(msg)
        if ".." in candidate.parts:
            msg = f"Resource declares a path {value!r} escaping the resource root via '..'."
            raise ValueError(msg)
        return value

    @field_validator("sha256")
    @classmethod
    def _sha256_shape(cls, value: str | None) -> str | None:
        if value is None:
            return None
        if len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
            msg = (
                f"Resource has a malformed sha256 {value!r} (expected 64 lowercase hex characters)."
            )
            raise ValueError(msg)
        return value

    @model_validator(mode="after")
    def _status_matches_fields(self) -> Self:
        """``sha256`` may only be non-null when the status says it was really fetched.

        A manifest claiming a hash for a resource that is absent, still growing or
        failed a sanity check is not a data-entry slip: it is an integrity claim
        that cannot be checked, which is worse than an honest gap.
        """
        fetched_fields = (self.sha256, self.size_bytes, self.retrieved)
        if self.status is ResourceStatus.FETCHED:
            if any(field is None for field in fetched_fields):
                msg = (
                    f"Resource {self.name!r} is marked 'fetched' but is missing sha256, "
                    "size_bytes or retrieved. A fetched resource must carry all three."
                )
                raise ValueError(msg)
        elif any(field is not None for field in fetched_fields) or self.integrity is not None:
            msg = (
                f"Resource {self.name!r} is marked {self.status.value!r} but declares "
                "sha256, size_bytes, retrieved or integrity. Those may only be set once "
                "the resource is actually fetched and verified — never write a hash for "
                "bytes that are absent, still growing, or failed a content check."
            )
            raise ValueError(msg)
        return self


class ResourceManifest(FrozenModel):
    """``knowledge/manifests/resources.yaml``, parsed."""

    model_config = ConfigDict(
        frozen=True, extra="ignore", str_strip_whitespace=True, hide_input_in_errors=True
    )

    manifest_version: int
    paths_relative_to: str
    generated_by: str
    resources: dict[str, ReferenceResource]

    def entries(self) -> tuple[ReferenceResource, ...]:
        """Every entry, in a total order (GP-30)."""
        return tuple(self.resources[name] for name in sorted(self.resources))

    def require(self, name: str) -> ReferenceResource:
        """The named entry, or a message that lists what IS registered."""
        try:
            return self.resources[name]
        except KeyError as exc:
            known = ", ".join(sorted(self.resources)[:12])
            msg = (
                f"No resource named {name!r} in the manifest ({len(self.resources)} registered; "
                f"first few: {known}). Register it in tools/acquire/catalog.py and re-run "
                "`uv run python -m tools.acquire write-manifest`."
            )
            raise ResourceError(msg) from exc


def load_resource_manifest(path: Path) -> ResourceManifest:
    """Parse and validate the resources manifest. Read-only; never writes."""
    if not path.is_file():
        msg = (
            f"Resource manifest not found: {path.as_posix()}. Generate it with "
            "`uv run python -m tools.acquire write-manifest`."
        )
        raise ResourceError(msg)
    try:
        raw: Any = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        msg = f"Resource manifest {path.name} is not valid YAML: {exc.__class__.__name__}"
        raise ResourceError(msg) from exc
    if not isinstance(raw, dict):
        msg = f"Resource manifest {path.name} must be a YAML mapping at the top level."
        raise ResourceError(msg)

    resources_raw = raw.get("resources")
    if isinstance(resources_raw, dict):
        for key, value in resources_raw.items():
            if isinstance(value, dict) and "name" not in value:
                value["name"] = key

    try:
        return ResourceManifest.model_validate(raw)
    except ValueError as exc:
        msg = f"Resource manifest {path.name} does not match the expected schema: {exc}"
        raise ResourceError(msg) from exc


# ---------------------------------------------------------------------------
# The resource root (mirrors mva.config.resolve_workspace)
# ---------------------------------------------------------------------------


class ResourceRoot(FrozenModel):
    """A validated external root holding the reference releases."""

    root: Path
    repo_root: Path

    def path(self, relative: str) -> Path:
        """Resolve a root-relative path, refusing an escape."""
        candidate = (self.root / relative).resolve()
        if not candidate.is_relative_to(self.root.resolve()):
            msg = f"Resource path {relative!r} escapes the resource root."
            raise ResourceRootError(msg)
        return candidate

    def locate(self, resource: ReferenceResource) -> Path:
        """The absolute path of a registered resource under this root."""
        return self.path(resource.path)


def resolve_resource_root(
    resource_root: str | Path | None = None,
    *,
    repo_root: Path | None = None,
    env: dict[str, str] | None = None,
    must_exist: bool = True,
) -> ResourceRoot:
    """Validate and return the external reference-data root.

    Deliberately shaped like :func:`mva.config.resolve_workspace`, and deliberately
    **without a default**. An earlier version of the acquisition tool fell back to a
    hard-coded ``~/Contri/bio-hackathon/mva-resources``, which is one contributor's
    layout: on any other machine it silently resolved to a directory that did not
    exist, and "resource missing" is indistinguishable from "resource not
    registered" once it reaches an adapter. A guessed path that is wrong is worse
    than no path at all, because it produces a plausible error somewhere else.

    Checks, in order: the path is provided (argument or ``$MVA_RESOURCES``); it
    exists; after full symlink resolution it is **not** inside the repository.

    Unlike the workspace, there is no cloud-sync check: everything under this root
    is public reference data, so the PRIV-07 threat (patient data uploaded to a
    third party) does not apply. Putting 185 GB in a synced folder is a bad idea
    for other reasons, and is the operator's to make.
    """
    repo = (repo_root or find_repo_root()).resolve()
    variables = env if env is not None else os.environ

    raw = resource_root if resource_root is not None else variables.get("MVA_RESOURCES")
    if raw is None:
        msg = (
            "No resource root configured. Set MVA_RESOURCES to the directory holding the "
            "public reference releases (gnomAD, ClinVar, the GRCh38 FASTA, MANE, SnpEff), "
            "or pass --resource-root. It must be OUTSIDE this repository: the registered "
            "set is ~203 GB and a stray 'git add -A' over it is not recoverable. There is "
            "deliberately no default — a guessed path that is wrong fails later, "
            "somewhere less obvious."
        )
        raise ResourceRootError(msg)

    root = Path(raw).expanduser()
    if must_exist and not root.is_dir():
        msg = (
            f"Resource root {root.as_posix()!r} does not exist (or is not a directory). "
            "Fetch the releases first: `uv run python -m tools.acquire fetch`."
        )
        raise ResourceRootError(msg)

    resolved = root.resolve()
    if path_is_within(resolved, repo):
        msg = (
            f"Resource root {resolved.as_posix()!r} resolves inside the repository "
            f"({repo.as_posix()}). Reference releases are hundreds of megabytes to tens of "
            "gigabytes each; inside the repo tree they are one 'git add -A' from being "
            "committed permanently. This check follows symlinks. Point MVA_RESOURCES "
            "outside the repo."
        )
        raise ResourceRootError(msg)

    return ResourceRoot(root=resolved, repo_root=repo)


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------


class IntegrityMode(StrEnum):
    """How hard to check a registered resource before trusting it.

    There is no ``off``. A mode that skips verification would be the one everybody
    selects under deadline pressure, and the pins would become decoration
    (ADR 0009: gates are blocking by design).
    """

    SPOT = "spot"
    """Size + sampled digest (:data:`SPOT_PLAN`). ~24 MiB per file. The default."""

    FULL = "full"
    """Every byte. 202.8 GB and 159 s for the registered set as it stands."""


class ResourceCheck(StrEnum):
    OK = "ok"
    MISSING = "missing"
    """Declared fetched, but nothing is on disk at the declared path."""

    SIZE_MISMATCH = "size_mismatch"
    """Present, but not the size the manifest pinned. Almost always a truncated or
    still-running download."""

    DIGEST_MISMATCH = "digest_mismatch"
    """Present and the right size, but the bytes are not the pinned bytes."""

    UNPINNED = "unpinned"
    """The manifest itself says this resource was never fetched. Nothing to verify;
    an honest gap, not an integrity failure."""

    PLAN_UNKNOWN = "plan_unknown"
    """The entry's spot digest was computed under a sampling plan this build does
    not implement. Refused rather than compared — two digests under different
    plans are not comparable, and treating a mismatch between them as corruption
    would be a false alarm while treating it as a pass would be a lie."""


class ResourceCheckResult(FrozenModel):
    """The outcome of checking one resource. Never carries file bytes (PRIV-09)."""

    name: str
    path: str
    check: ResourceCheck
    mode: IntegrityMode
    message: str

    @property
    def ok(self) -> bool:
        return self.check in {ResourceCheck.OK, ResourceCheck.UNPINNED}


def _as_root(root: ResourceRoot | Path) -> ResourceRoot:
    """Accept either a validated :class:`ResourceRoot` or a bare directory.

    Verification only needs a base directory to resolve paths against; the
    containment and existence rules are :func:`resolve_resource_root`'s job and are
    not re-litigated here. Accepting a bare ``Path`` keeps the acquisition tool
    (which resolves the root itself, before the directory may even exist) able to
    reuse this exact code rather than growing a second, divergent verifier.
    """
    if isinstance(root, ResourceRoot):
        return root
    resolved = Path(root).resolve()
    return ResourceRoot(root=resolved, repo_root=resolved)


def _measured_size(target: Path, kind: ResourceKind) -> int:
    if kind is ResourceKind.DIRECTORY:
        return sum(
            p.stat().st_size for p in target.rglob("*") if p.is_file() and not p.is_symlink()
        )
    return target.stat().st_size


def verify_resource(
    root: ResourceRoot | Path,
    resource: ReferenceResource,
    *,
    mode: IntegrityMode = IntegrityMode.SPOT,
) -> ResourceCheckResult:
    """Check one registered resource's local bytes against its manifest pin."""
    root = _as_root(root)
    target = root.locate(resource)
    unpinned = ResourceCheckResult(
        name=resource.name,
        path=resource.path,
        check=ResourceCheck.UNPINNED,
        mode=mode,
        message=(
            f"resource {resource.name!r} is not pinned in the manifest "
            f"({resource.notes or 'never fetched'}); nothing to verify"
        ),
    )
    if resource.status is not ResourceStatus.FETCHED or resource.sha256 is None:
        return unpinned

    exists = target.is_dir() if resource.kind is ResourceKind.DIRECTORY else target.is_file()
    if not exists:
        return ResourceCheckResult(
            name=resource.name,
            path=resource.path,
            check=ResourceCheck.MISSING,
            mode=mode,
            message=(
                f"resource {resource.name!r} is pinned in the manifest but absent from disk "
                f"(expected {resource.path!r} under {root.root.as_posix()})"
            ),
        )

    # Size first: it is two stat calls, it is the failure that actually happens
    # (a truncated or in-flight download), and it gives a far more useful message
    # than a digest mismatch would.
    actual_size = _measured_size(target, resource.kind)
    if resource.size_bytes is not None and actual_size != resource.size_bytes:
        short = resource.size_bytes - actual_size
        hint = (
            f" — {short} bytes short, which is what an interrupted download looks like"
            if short > 0
            else ""
        )
        return ResourceCheckResult(
            name=resource.name,
            path=resource.path,
            check=ResourceCheck.SIZE_MISMATCH,
            mode=mode,
            message=(
                f"resource {resource.name!r} ({resource.path!r}) is {actual_size} bytes; the "
                f"manifest pins {resource.size_bytes}{hint}. Re-fetch it, then regenerate the "
                "manifest — do not edit the pinned size."
            ),
        )

    if mode is IntegrityMode.FULL:
        expected, actual, label = (
            resource.sha256,
            tree_digest(target) if resource.kind is ResourceKind.DIRECTORY else hash_file(target),
            "sha256",
        )
    else:
        if resource.integrity is None:
            return unpinned.model_copy(
                update={
                    "check": ResourceCheck.PLAN_UNKNOWN,
                    "message": (
                        f"resource {resource.name!r} carries no integrity record, so there is "
                        "no spot digest to compare against. Regenerate the manifest with "
                        "`uv run python -m tools.acquire write-manifest`, or verify with "
                        "--mode full."
                    ),
                }
            )
        if resource.integrity.spot_plan != SPOT_PLAN:
            return ResourceCheckResult(
                name=resource.name,
                path=resource.path,
                check=ResourceCheck.PLAN_UNKNOWN,
                mode=mode,
                message=(
                    f"resource {resource.name!r} pins a spot digest under plan "
                    f"{resource.integrity.spot_plan!r}, but this build implements "
                    f"{SPOT_PLAN!r}. Digests under different sampling plans are not "
                    "comparable. Re-register the manifest, or verify with --mode full."
                ),
            )
        expected, actual, label = (
            resource.integrity.spot_sha256,
            spot_tree_digest(target)
            if resource.kind is ResourceKind.DIRECTORY
            else spot_digest(target),
            "spot digest",
        )

    if actual != expected:
        coverage = (
            "the full file" if mode is IntegrityMode.FULL else f"a {SPOT_PLAN} sample of the file"
        )
        return ResourceCheckResult(
            name=resource.name,
            path=resource.path,
            check=ResourceCheck.DIGEST_MISMATCH,
            mode=mode,
            message=(
                f"resource {resource.name!r} ({resource.path!r}) failed its manifest "
                f"integrity check over {coverage}: expected {label} {expected}, found "
                f"{actual}. The bytes changed without the manifest being regenerated."
            ),
        )

    return ResourceCheckResult(
        name=resource.name,
        path=resource.path,
        check=ResourceCheck.OK,
        mode=mode,
        message=f"resource {resource.name!r} matches its pinned {label}",
    )


def verify_resources(
    root: ResourceRoot | Path,
    resources: Iterable[ReferenceResource],
    *,
    mode: IntegrityMode = IntegrityMode.SPOT,
) -> tuple[ResourceCheckResult, ...]:
    """:func:`verify_resource` over every entry, in the given order."""
    return tuple(verify_resource(root, resource, mode=mode) for resource in resources)


def assert_resources_verified(
    root: ResourceRoot | Path,
    resources: Iterable[ReferenceResource],
    *,
    mode: IntegrityMode = IntegrityMode.SPOT,
) -> tuple[ResourceCheckResult, ...]:
    """Verify every entry, raising :class:`ResourceIntegrityError` on any failure.

    Fails closed: ``UNPINNED`` passes (an honestly-declared gap is not corruption),
    everything else raises. The message names each failing resource and its
    declared path, and never its contents (PRIV-09).
    """
    results = verify_resources(root, resources, mode=mode)
    failures = [result for result in results if not result.ok]
    if failures:
        detail = "\n".join(
            f"  {result.check.value.upper()}: {result.message}" for result in failures
        )
        msg = (
            f"{len(failures)} of {len(results)} reference resources failed integrity "
            f"verification (mode={mode.value}):\n{detail}"
        )
        raise ResourceIntegrityError(msg)
    return results


def required_resources(
    manifest: ResourceManifest, names: Sequence[str]
) -> tuple[ReferenceResource, ...]:
    """The named entries, raising if any is unregistered or unfetched.

    The point of separating this from :func:`verify_resources` is GP-14: a run that
    needs ClinVar must fail loudly when ClinVar is unfetched, but the *manifest* is
    allowed to carry unfetched entries without that being an error. Which absences
    matter is a property of the run, not of the manifest.
    """
    selected: list[ReferenceResource] = []
    unfetched: list[str] = []
    for name in names:
        resource = manifest.require(name)
        if resource.status is not ResourceStatus.FETCHED:
            unfetched.append(f"{name} ({resource.notes or 'not fetched'})")
        selected.append(resource)
    if unfetched:
        msg = (
            "Reference resources required by this run are not fetched: "
            + "; ".join(unfetched)
            + ". Fetch them with `uv run python -m tools.acquire fetch`, then regenerate the "
            "manifest. Running without them is not a degraded mode that was chosen — it is "
            "one that would happen silently."
        )
        raise ResourceError(msg)
    return tuple(selected)


# ---------------------------------------------------------------------------
# The one resolution the composition root actually asks for
# ---------------------------------------------------------------------------


def reference_fasta_path(
    config: CaseConfig,
    *,
    workspace: Workspace | None = None,
    resource_root: ResourceRoot | None = None,
) -> Path | None:
    """Where the GRCh38 FASTA for this case is, or ``None`` if there isn't one.

    Two sources, in precedence order, because they mean different things:

    1. ``config.inputs.reference_fasta`` — a **workspace-relative** override, for a
       case that ships its own reference (a different build, a subsetted FASTA).
       Explicit per case, so it wins.
    2. ``config.resources.reference_fasta`` — the shared public GRCh38 release under
       the resource root. This is the normal path.

    Returning ``None`` rather than raising is deliberate and is GP-14 applied to
    configuration: normalisation without a reference is a *degraded* mode, not an
    impossible one, and the caller must be the one to decide whether to proceed and
    to surface ``representation_limitation`` in the run warnings. A function that
    raised here would push that decision into an import-time accident.
    """
    override = config.inputs.reference_fasta
    if override is not None:
        if workspace is None:
            msg = (
                f"Case {config.case_id!r} sets inputs.reference_fasta={override!r}, which is "
                "workspace-relative, but no workspace was supplied to resolve it against."
            )
            raise ResourceError(msg)
        return workspace.path(override)

    if resource_root is None:
        return None
    return resource_root.path(config.resources.reference_fasta)
