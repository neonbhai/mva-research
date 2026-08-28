"""Building the REAL adapter set out of the out-of-repo reference releases.

``local_tables`` wires the synthetic stand-ins from hash-pinned TSVs under
``knowledge/public/``. This module is its counterpart for a real case: it turns
:class:`~mva.config.ResourceSettings` plus a validated
:class:`~mva.resources.ResourceRoot` into a :class:`~mva.annotation.base.AdapterSet`
backed by ClinVar, gnomAD, SnpEff and MANE. The *policy* — which of the two a
given case gets — stays at the composition root in :mod:`mva.orchestrator`, where
a reviewer looks for it (GP-03).

Three things here are load-bearing, and each exists because getting it wrong is
invisible rather than loud.

**1. The reference FASTA reaches both joining adapters.**
:class:`~mva.annotation.clinvar_vcf.ClinvarVcfAdapter` and
:class:`~mva.annotation.gnomad_sites.GnomadSitesFrequencyAdapter` both take
``reference=`` as an *optional* keyword. Constructed without one they still work,
still pass their own tests, and silently run in a trim-only mode that cannot
reconcile a right-shifted indel — which is the default. Measured on the real
releases:

* gnomAD, chr21, a 520 kb exonic window: 1,382 indel records, 1,029 (74.5%) in
  repeat tracts. Right-shifted spellings joining **without** a reference: 0 of
  1,029. **With** one: 989. Thirty of the recovered records are variants gnomAD
  itself calls common, 22 above 5% AF — scored as novel and ultra-rare without
  the reference.
* ClinVar, chr17:43,000,000-43,520,000 against the 2026-08-22 release: 3,595
  indel ALT alleles, 2,215 (61.6%) in repeat tracts, **0** joining without a
  reference against **2,211** with one — of which **1,761 are assertions ClinVar
  calls Pathogenic or Likely pathogenic**.

So ``reference`` is a *required, non-defaulted* keyword on
:func:`build_real_adapter_set` and on the two per-slot builders. Passing ``None``
is still permitted, because a degraded run is a real state (GP-14) and refusing
to represent it would only push the degradation somewhere less visible — but it
has to be typed out at the call site, and it is reported in
:attr:`BoundAdapters.warnings` rather than left for a reader to infer from an
absent assertion. See ADR 0018 and ADR 0027.

**2. gnomAD v4.1 exomes has no chrM shard.** Not a download that failed: the
release does not contain one. :meth:`GnomadSitesFrequencyAdapter.frequencies`
therefore *raises* on any mitochondrial variant even against a complete release,
which would abort a whole run over one variant. Refusing the run is wrong and
silently omitting the variant is worse, so the frequency slot is wrapped in
:class:`PartialCoverageFrequencyAdapter`, which answers what can be answered and
**names the hole** in the run warnings (GP-14).

**3. SnpEff is pinned from the installer's manifest, never from ``measure()``.**
:meth:`SnpEffArtifactPins.measure` hashes whatever happens to be installed, which
pins the installation to itself and verifies nothing.
:meth:`SnpEffArtifactPins.from_manifest` reads the digests
``tools/setup/install_snpeff.sh`` wrote *after* checking each artifact against the
reviewed constants in the script, so it is an assertion about which bytes were
reviewed rather than a recording of which bytes are present.

Nothing in this module fetches anything, and no adapter it builds may receive a
proband coordinate over a network: every source is a local file the acquisition
step already downloaded (PRIV-05).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final, Protocol

from mva.alleles import ReferenceLookup
from mva.annotation.base import AdapterSet, ConsequenceAdapter
from mva.annotation.clinvar_vcf import ClinvarVcfAdapter
from mva.annotation.gene_intervals import (
    GeneBackfillConsequenceAdapter,
    ManeGeneAdapter,
    ManeGeneIndex,
)
from mva.annotation.gnomad_sites import GnomadSitesFrequencyAdapter
from mva.annotation.snpeff_local import SnpEffArtifactPins, SnpEffConsequenceAdapter
from mva.config import CaseConfig
from mva.models.genome import GenomeBuild, contig_sort_key, normalise_contig
from mva.models.variant import PopulationFrequency
from mva.resources import (
    IntegrityMode,
    ReferenceResource,
    ResourceError,
    ResourceManifest,
    ResourceRoot,
    assert_resources_verified,
    required_resources,
)

# ---------------------------------------------------------------------------
# Which registered resources a real run cannot proceed without
# ---------------------------------------------------------------------------

#: Config attributes naming the single-file resources every real run needs, in
#: the order they are reported. Named as *config attributes* rather than as
#: manifest keys because the config is what a case actually reads: a manifest key
#: renamed upstream must not silently unbind an adapter, and looking the entry up
#: by the path the config declares makes the two impossible to disagree about.
_REQUIRED_PATH_SETTINGS: Final[tuple[str, ...]] = (
    "reference_fasta",
    "reference_fasta_index",
    "clinvar_vcf",
    "clinvar_vcf_index",
    "mane_gtf",
    "mane_summary",
    "snpeff_jar",
    "snpeff_config",
)


class AdapterBindingError(ResourceError):
    """A real adapter set cannot be built from the configured resources.

    A subclass of :class:`~mva.resources.ResourceError` because that is what it
    is: every cause is a reference release that is missing, unregistered, or not
    what the manifest pinned. Messages name resources and paths, never file
    contents and never a patient coordinate (PRIV-09).
    """


# ---------------------------------------------------------------------------
# The frequency slot's coverage wrapper
# ---------------------------------------------------------------------------


class PartialCoverageFrequencyAdapter:
    """A frequency adapter that survives a release's coverage holes, and names them.

    gnomAD v4.1 exomes ships shards for chr1-22, X and Y and **none for chrM**.
    That is a property of the release, not a failed download, so
    :meth:`GnomadSitesFrequencyAdapter.frequencies` — which fails closed on any
    contig it cannot query — raises on a mitochondrial variant against an
    otherwise complete install. A proband callset routinely holds a few dozen
    chrM calls, so the composition root has exactly three options: abort the whole
    run over them, drop them silently, or answer what is answerable and say what
    was not asked. This is the third.

    It delegates to :meth:`GnomadSitesFrequencyAdapter.lookup_partial`, which
    returns the gap as data. Variants on an unqueryable contig are **absent from
    the result**, which the annotation stage already renders as "no frequency
    data" — a neutral no-data evidence item and a mid-range rarity score, never an
    allele frequency of zero (GP-14). What this class adds is that the absence is
    *counted and reported*: :meth:`warnings` states which contigs could not be
    queried and how many variants that affected, so the run manifest carries the
    hole instead of a reader having to deduce it.

    Identity is delegated, not invented: ``name`` and ``version`` are the wrapped
    adapter's, because this **is** the gnomAD adapter with its coverage gap made
    explicit, and every ``EvidenceItem`` must name the release it came from
    (GP-18).
    """

    def __init__(self, inner: GnomadSitesFrequencyAdapter) -> None:
        self._inner = inner
        #: contig -> variants that could not be looked up on it. A dict rather
        #: than a set so the count survives, and reported in karyotype order so
        #: repeat runs emit byte-identical warnings (GP-30).
        self._gaps: dict[str, int] = {}
        self._batches_with_gaps = 0

    # ------------------------------------------------------------------ identity

    @property
    def name(self) -> str:
        return self._inner.name

    @property
    def version(self) -> str:
        return self._inner.version

    @property
    def synthetic(self) -> bool:
        """False: every number this returns was read out of a gnomAD sites VCF."""
        return self._inner.synthetic

    @property
    def inner(self) -> GnomadSitesFrequencyAdapter:
        """The wrapped adapter, for a caller that needs its release-level facts."""
        return self._inner

    # -------------------------------------------------------------------- lookup

    def frequencies(
        self, variant_ids: Sequence[str]
    ) -> Mapping[str, tuple[PopulationFrequency, ...]]:
        """Frequencies for every variant on a contig this release covers.

        A variant gnomAD holds no record for is absent, exactly as
        :meth:`GnomadSitesFrequencyAdapter.frequencies` leaves it. A variant on a
        contig the release does not cover is *also* absent — a different fact,
        which is why it is counted here and stated by :meth:`warnings` rather than
        being allowed to look like the first.
        """
        lookup = self._inner.lookup_partial(variant_ids)
        if not lookup.is_complete:
            self._batches_with_gaps += 1
            for contig in lookup.unqueryable_contigs:
                self._gaps[contig] = self._gaps.get(contig, 0) + sum(
                    1
                    for variant_id in lookup.unqueryable_variant_ids
                    if _contig_of(variant_id) == contig
                )
        return lookup.frequencies

    def close(self) -> None:
        """Release the wrapped adapter's shard handles. Idempotent."""
        self._inner.close()

    # ----------------------------------------------------------------- reporting

    @property
    def unqueryable_contigs(self) -> tuple[str, ...]:
        """Contigs asked about that this release does not cover, karyotype-ordered."""
        return tuple(sorted(self._gaps, key=contig_sort_key))

    @property
    def unqueryable_count(self) -> int:
        """How many variants could not be looked up at all. A count, never an ID."""
        return sum(self._gaps.values())

    def warnings(self) -> tuple[str, ...]:
        """The coverage hole, in words, for the run manifest. Counts only (GP-41)."""
        if not self._gaps:
            return ()
        breakdown = ", ".join(
            f"{contig}={self._gaps[contig]}" for contig in self.unqueryable_contigs
        )
        return (
            f"{self.unqueryable_count} variant(s) could not be looked up in "
            f"{self._inner.name}@{self._inner.version}: it has no shard covering "
            f"{breakdown}. gnomAD v4.1 exomes ships no chrM shard at all, so a "
            "mitochondrial call reaches this path even against a complete release. "
            "These variants are recorded as having NO frequency data, which is absence "
            "of a resource and not evidence of rarity, and they must never be scored as "
            "allele frequency 0 (GP-14). No coordinate is echoed (PRIV-09).",
        )


def _contig_of(variant_id: str) -> str | None:
    """The contig of a canonical ``build:contig:pos:ref:alt`` ID, or ``None``.

    Tolerant on purpose: this is used only to attribute a coverage gap to a
    contig for a *count*, and a malformed ID must not be able to turn a reporting
    line into a crash on the patient path.
    """
    parts = variant_id.split(":")
    if len(parts) < 2:
        return None
    try:
        return normalise_contig(parts[1])
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# Resolution: paths, versions and pins, all verified before anything is opened
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ResolvedResources:
    """Absolute paths and release versions for one real binding, already verified.

    Deliberately separate from :func:`build_real_adapter_set`. Resolution is where
    "is this release present, registered and intact?" is answered and where a run
    fails closed; construction is where handles are opened. Keeping them apart is
    what makes the verification non-optional: there is no ``verify=False``
    parameter to reach for, because verification is a different function that the
    production path always calls (ADR 0020, ADR 0027).
    """

    reference_fasta: Path
    clinvar_vcf: Path
    clinvar_sha256: str
    gnomad_dir: Path
    gnomad_release: str
    gnomad_contigs: tuple[str, ...]
    mane_gtf: Path
    mane_gtf_sha256: str
    mane_summary: Path
    mane_summary_sha256: str
    snpeff_jar: Path
    snpeff_config: Path
    snpeff_data_dir: Path
    snpeff_genome: str
    snpeff_java_binary: Path
    snpeff_pins: Path
    build: GenomeBuild = GenomeBuild.GRCH38
    verified_count: int = 0
    """How many registered resources passed their manifest pin. A count, so the
    provenance can say 'n resources verified' without carrying 60 paths."""

    warnings: tuple[str, ...] = ()
    """Anything true about this resolution that is not a refusal. Folded into
    :attr:`BoundAdapters.warnings`, so it reaches the run manifest."""


def resolve_real_resources(
    config: CaseConfig,
    *,
    resource_root: ResourceRoot,
    manifest: ResourceManifest,
) -> ResolvedResources:
    """Resolve, require and verify every release the real adapter set needs.

    Order matters and is the point: entries are located, *then* required to be
    fetched, *then* verified against their manifest pins, and only then is a path
    handed to a constructor. An adapter opened over an unverified file is
    annotating against unknown bytes while reporting a pinned release
    (ADR 0020).

    Verification honours ``config.resources.integrity_mode``: ``spot`` (the
    default) reads at most 24 MiB per file, ``full`` re-reads every byte. There is
    deliberately no third mode.

    Raises:
        AdapterBindingError: a resource this run needs is unregistered, not
            fetched, or fails its pin. Every message names the resource and its
            declared path, and never its contents (PRIV-09).
    """
    entries: list[ReferenceResource] = []
    paths: dict[str, Path] = {}
    for setting in _REQUIRED_PATH_SETTINGS:
        relative = getattr(config.resources, setting)
        if not isinstance(relative, str):  # pragma: no cover - config is typed
            msg = f"resources.{setting} is not a path string."
            raise AdapterBindingError(msg)
        entries.append(_entry_for_path(manifest, relative, setting=setting))
        paths[setting] = resource_root.path(relative)

    database_relative = f"{config.resources.snpeff_data_dir}/{config.resources.snpeff_genome}"
    entries.append(_entry_for_path(manifest, database_relative, setting="snpeff genome database"))

    shards = _gnomad_shards(config, manifest)
    entries.extend(shards)

    mode = IntegrityMode(config.resources.integrity_mode)
    try:
        required = required_resources(manifest, [entry.name for entry in entries])
        results = assert_resources_verified(resource_root, required, mode=mode)
    except ResourceError as exc:
        msg = (
            f"Case {config.case_id!r} is not synthetic, so its annotation must come from the "
            f"real reference releases under {resource_root.root.as_posix()} — and they are not "
            f"usable: {exc}\n"
            "Refusing to fall back to the synthetic tables in knowledge/public/: they hold "
            "fictional genes and invented allele frequencies, so the run would rank a real "
            "proband against fabricated evidence while every artifact looked healthy "
            "(GP-20, ADR 0027)."
        )
        raise AdapterBindingError(msg) from exc

    clinvar_entry = _entry_for_path(manifest, config.resources.clinvar_vcf, setting="clinvar_vcf")
    gtf_entry = _entry_for_path(manifest, config.resources.mane_gtf, setting="mane_gtf")
    summary_entry = _entry_for_path(manifest, config.resources.mane_summary, setting="mane_summary")
    release = _gnomad_release(shards)
    warnings = _unregistered_shard_warning(config, resource_root, shards)

    return ResolvedResources(
        reference_fasta=paths["reference_fasta"],
        clinvar_vcf=paths["clinvar_vcf"],
        clinvar_sha256=_pinned_sha256(clinvar_entry),
        gnomad_dir=resource_root.path(config.resources.gnomad_exomes_dir),
        gnomad_release=release,
        gnomad_contigs=_shard_contigs(config, shards),
        mane_gtf=paths["mane_gtf"],
        mane_gtf_sha256=_pinned_sha256(gtf_entry),
        mane_summary=paths["mane_summary"],
        mane_summary_sha256=_pinned_sha256(summary_entry),
        snpeff_jar=paths["snpeff_jar"],
        snpeff_config=paths["snpeff_config"],
        snpeff_data_dir=resource_root.path(config.resources.snpeff_data_dir),
        snpeff_genome=config.resources.snpeff_genome,
        snpeff_java_binary=_jre_path(resource_root, config.resources.snpeff_java_binary),
        snpeff_pins=resource_root.path(config.resources.snpeff_pins),
        build=config.genome_build,
        verified_count=len(results),
        warnings=warnings,
    )


def _jre_path(resource_root: ResourceRoot, relative: str) -> Path:
    """The configured JRE under the resource root, WITHOUT following symlinks.

    The one place this module declines to use :meth:`ResourceRoot.path`, and the
    reason is specific rather than convenient. That method resolves symlinks
    before testing containment, which is right for data — a manifest entry must
    not be able to point outside the root at something unpinned. The JRE is not
    data. ``tools/setup/install_snpeff.sh`` deliberately symlinks a real system
    JDK to ``<root>/snpeff/jdk`` precisely so the 300 MB runtime is not copied
    into the reference tree, so resolving it lands on ``/opt/...`` and the
    containment check refuses a correct, intended install.

    Containment is still enforced, lexically and earlier:
    :class:`~mva.config.ResourceSettings` refuses an absolute path and any ``..``
    segment for this setting, so the joined path cannot escape by spelling. What
    it may do is point *through* an operator-created symlink, which is the whole
    arrangement. :func:`~mva.annotation.snpeff_local.resolve_java_binary` then
    requires it to be a real file and the adapter runs it before trusting it.
    """
    return resource_root.root / relative


def _entry_for_path(
    manifest: ResourceManifest, relative: str, *, setting: str
) -> ReferenceResource:
    """The registered entry at ``relative``, or a message naming what to register."""
    for entry in manifest.entries():
        if entry.path == relative:
            return entry
    msg = (
        f"No resource registered at {relative!r} (configured as {setting}), so its bytes are "
        "unpinned and nothing can state which release they are. Register it in "
        "tools/acquire/catalog.py and re-run `uv run python -m tools.acquire write-manifest`. "
        "An unregistered reference file is not a smaller problem than a missing one: it makes "
        "every claim derived from it unreproducible (GP-18)."
    )
    raise AdapterBindingError(msg)


def _pinned_sha256(entry: ReferenceResource) -> str:
    """The entry's full digest, which :func:`required_resources` has already proven."""
    if entry.sha256 is None:  # pragma: no cover - required_resources raises first
        msg = f"Resource {entry.name!r} carries no sha256 despite being marked fetched."
        raise AdapterBindingError(msg)
    return entry.sha256


def _gnomad_shards(config: CaseConfig, manifest: ResourceManifest) -> tuple[ReferenceResource, ...]:
    """Every registered sites shard under the configured gnomAD directory.

    Derived from the manifest rather than from a hard-coded contig list, because
    which contigs a release ships is a property of the release: v4.1 exomes has
    chr1-22, X and Y and no chrM. Requiring a fixed 24 would fail every complete
    install; requiring none would let a run start against an empty directory and
    report the whole genome as novel.
    """
    prefix = f"{config.resources.gnomad_exomes_dir}/"
    shards = tuple(
        entry
        for entry in manifest.entries()
        if entry.path.startswith(prefix) and not entry.path.endswith(".tbi")
    )
    if not shards:
        msg = (
            f"No gnomAD sites shard is registered under {config.resources.gnomad_exomes_dir!r}. "
            "With no frequency source every variant would report as having no frequency data, "
            "and the rarity component would score the entire callset on absence rather than on "
            "evidence (GP-14). Fetch the release and regenerate the manifest."
        )
        raise AdapterBindingError(msg)
    return shards


def _unregistered_shard_warning(
    config: CaseConfig,
    resource_root: ResourceRoot,
    shards: Sequence[ReferenceResource],
) -> tuple[str, ...]:
    """Warn when the gnomAD directory holds a sites file the manifest does not pin.

    The adapter opens **every** sites file it finds under the directory, while
    verification can only check what is registered. Those two sets are normally
    identical and there is nothing to say. When they are not, the run would be
    answering frequency questions out of bytes nothing has pinned, while the
    provenance manifest reported n resources verified — a verified-something,
    used-something-else failure, which is the shape this repository treats as
    worse than an outright gap (GP-18, ADR 0020). Reported rather than refused:
    an extra shard is not evidence of tampering, and deleting a caller's data on
    that suspicion would be its own failure.
    """
    directory = resource_root.path(config.resources.gnomad_exomes_dir)
    if not directory.is_dir():
        return ()
    registered = {shard.path.rsplit("/", 1)[-1] for shard in shards}
    present = sorted(
        path.name
        for path in directory.iterdir()
        if path.is_file() and path.name.endswith((".vcf.bgz", ".vcf.gz"))
    )
    unregistered = [name for name in present if name not in registered]
    if not unregistered:
        return ()
    return (
        f"{len(unregistered)} gnomAD sites file(s) under "
        f"{config.resources.gnomad_exomes_dir!r} are not registered in the resource "
        f"manifest and were therefore NOT integrity-checked: {', '.join(unregistered)}. "
        "The adapter opens every sites file it finds, so frequencies may come from "
        "unpinned bytes while the provenance manifest reports a verified release. "
        "Register them with `uv run python -m tools.acquire write-manifest`, or remove "
        "them from the directory.",
    )


def _gnomad_release(shards: Sequence[ReferenceResource]) -> str:
    """The release token every shard agrees on, in the ``v4.1`` form the adapter checks.

    The gnomAD sites VCF header carries no release string, so the version comes
    from the acquisition step and is then cross-checked against every filename by
    the adapter itself. A directory mixing releases is refused here, before that
    check, because the message can name both versions.
    """
    versions = sorted({shard.version for shard in shards})
    if len(versions) != 1:
        msg = (
            f"Registered gnomAD shards declare more than one release ({', '.join(versions)}). "
            "Frequencies from different releases carry different cohort sizes and different "
            "ancestry coverage; labelling them as one misstates every allele number the "
            "ranking guards on (ADR 0010, GP-18)."
        )
        raise AdapterBindingError(msg)
    version = versions[0]
    return version if version.startswith("v") else f"v{version}"


def _shard_contigs(config: CaseConfig, shards: Sequence[ReferenceResource]) -> tuple[str, ...]:
    """Contigs recovered from the shard filenames, karyotype-ordered.

    Passed to the adapter as ``require_contigs`` so that a registered shard that
    is truncated, still arriving or unreadable fails at construction rather than
    at the first lookup that touches it — before the run has done any work.
    """
    template = config.resources.gnomad_exomes_filename
    head, _, tail = template.partition("{contig}")
    contigs: set[str] = set()
    for shard in shards:
        name = shard.path.rsplit("/", 1)[-1]
        if not name.startswith(head) or not name.endswith(tail):
            continue
        token = name[len(head) : len(name) - len(tail)] if tail else name[len(head) :]
        try:
            contigs.add(normalise_contig(token))
        except ValueError:
            continue
    return tuple(sorted(contigs, key=contig_sort_key))


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------


class _Closeable(Protocol):
    """Anything holding an OS handle this module is responsible for releasing."""

    def close(self) -> None: ...


@dataclass(slots=True)
class BoundAdapters:
    """An adapter set, what it cannot do, and how to let go of its handles.

    Returned by both branches of the composition root — the real one here and the
    synthetic one in :mod:`mva.orchestrator` — so the call site does not grow a
    conditional around closing handles or collecting warnings.
    """

    adapters: AdapterSet
    warnings: tuple[str, ...] = ()
    """Stated at BIND time: what a reader needs to know about this set's
    limitations before it has annotated anything. Holes discovered *during* the
    run are added by :meth:`run_warnings`."""

    coverage: PartialCoverageFrequencyAdapter | None = None
    """The frequency slot's coverage reporter, when the frequency source is one
    that can have holes. ``None`` for the synthetic tables, which cover exactly
    what they contain."""

    clinical: ClinvarVcfAdapter | None = None
    frequency: GnomadSitesFrequencyAdapter | None = None
    """The two joining adapters themselves, kept so their representation status can
    be read **after** the run rather than copied at construction. ``None`` for the
    synthetic branch, whose tables have no reference to break.

    Not a tuple of pre-rendered strings, deliberately.
    ``representation_limitation`` is derived from ``_unusable_reference_alleles``,
    a counter that is zero until a lookup has failed — and every lookup happens
    after this object exists. Anything captured here as a *value* answers "was the
    reference broken before we started", to which the answer is always no."""

    closers: tuple[_Closeable, ...] = ()
    """Everything holding an OS handle, closed in construction order."""

    def run_warnings(self) -> tuple[str, ...]:
        """Bind-time warnings plus anything the run itself discovered.

        Call it after annotation has drained. Two kinds of hole can only be known
        by then, and both used to be reported as if they could not happen:

        * :class:`PartialCoverageFrequencyAdapter` cannot know which contigs were
          asked about until they have been.
        * A joining adapter cannot know its reference was unreadable until it has
          tried to read it. ``representation_warnings`` was previously evaluated
          inside the ``BoundAdapters(...)`` call in
          :func:`build_real_adapter_set`, which is to say before the first lookup,
          so a run whose FASTA raised on **every** read still reported
          ``APPLIED`` with no limitation. Indels it could not left-align came back
          as "no gnomAD record" and were scored novel and ultra-rare — GP-14
          inverted into a manufactured false positive.
        """
        gaps = self.coverage.warnings() if self.coverage is not None else ()
        representation = (
            representation_warnings(clinical=self.clinical, frequency=self.frequency)
            if self.clinical is not None and self.frequency is not None
            else ()
        )
        return (*self.warnings, *representation, *gaps)

    def close(self) -> None:
        """Release every open file handle. Idempotent; safe in a ``finally``."""
        for closer in self.closers:
            closer.close()


def build_real_adapter_set(
    resolved: ResolvedResources, *, reference: ReferenceLookup | None
) -> BoundAdapters:
    """Bind ClinVar, gnomAD, SnpEff and MANE into one :class:`AdapterSet`.

    ``reference`` has no default. That is the whole point of this function: both
    joining adapters accept ``reference=None`` and degrade silently, so the
    keyword is made impossible to *forget* even though it remains possible to pass
    ``None`` deliberately. When it is ``None``, both adapters' own
    ``representation_limitation`` sentences are surfaced in
    :attr:`BoundAdapters.warnings`, so a degraded run says so instead of looking
    healthy (ADR 0018).

    Handles are opened in cheapest-to-most-expensive order and every one opened so
    far is closed if a later one fails, so a misconfigured SnpEff does not leak a
    ClinVar tabix handle for the life of the process.
    """
    closers: list[_Closeable] = []
    try:
        clinical = build_clinical_adapter(resolved, reference=reference)
        closers.append(clinical)
        frequency = build_frequency_adapter(resolved, reference=reference)
        closers.append(frequency)
        consequence = build_consequence_adapter(resolved)
    except BaseException:
        for closer in closers:
            closer.close()
        raise

    # ``warnings`` carries only what is knowable NOW. The representation status is
    # deliberately NOT among it: it is derived from a per-lookup failure counter
    # that no lookup has yet had the chance to increment, so evaluating it here
    # freezes a healthy answer over a run that has not started. Both adapters are
    # handed to BoundAdapters instead, and `run_warnings()` asks them at the end.
    return BoundAdapters(
        adapters=AdapterSet(consequence=consequence, frequency=frequency, clinical=clinical),
        warnings=tuple(resolved.warnings),
        coverage=frequency,
        clinical=clinical,
        frequency=frequency.inner,
        closers=tuple(closers),
    )


def representation_warnings(
    *, clinical: ClinvarVcfAdapter, frequency: GnomadSitesFrequencyAdapter
) -> tuple[str, ...]:
    """Each adapter's own account of what it cannot reconcile, or ``()``.

    Reads ``representation_limitation`` off the constructed adapters rather than
    re-deriving the condition from whether a reference was passed. The adapters
    are the authority on their own degraded state, and asking them means this
    warning cannot drift away from what they actually do.

    Under the composition root's current policy a real run always has a reference
    (ADR 0027), so an empty tuple is the expected result and a non-empty one is a
    tripwire: something bound a real adapter without the FASTA and the run is
    about to score right-shifted indels as novel.
    """
    warnings: list[str] = []
    for label, limitation, status in (
        ("clinical", clinical.representation_limitation, clinical.representation_status),
        ("frequency", frequency.representation_limitation, frequency.representation_status),
    ):
        if limitation is None:
            continue
        warnings.append(
            f"DEGRADED {label} join ({status.value}): {limitation} Measured cost on the real "
            "releases: 2,211 of 2,215 repeat-tract ClinVar indel assertions unreachable, 1,761 "
            "of them Pathogenic or Likely pathogenic; 989 of 1,029 repeat-tract gnomAD indel "
            "records unreachable, 30 of them variants gnomAD calls common. See ADR 0018."
        )
    return tuple(warnings)


def build_clinical_adapter(
    resolved: ResolvedResources, *, reference: ReferenceLookup | None
) -> ClinvarVcfAdapter:
    """The ClinVar slot, pinned to the manifest digest and given the reference."""
    return ClinvarVcfAdapter(
        resolved.clinvar_vcf,
        expected_sha256=resolved.clinvar_sha256,
        build=resolved.build,
        reference=reference,
    )


def build_frequency_adapter(
    resolved: ResolvedResources, *, reference: ReferenceLookup | None
) -> PartialCoverageFrequencyAdapter:
    """The gnomAD slot, given the reference and wrapped for the chrM-shaped hole."""
    inner = GnomadSitesFrequencyAdapter(
        resolved.gnomad_dir,
        release=resolved.gnomad_release,
        subset="exomes",
        require_contigs=resolved.gnomad_contigs,
        reference=reference,
    )
    return PartialCoverageFrequencyAdapter(inner)


def build_consequence_adapter(resolved: ResolvedResources) -> GeneBackfillConsequenceAdapter:
    """SnpEff in front, the MANE interval join behind it.

    :class:`GeneBackfillConsequenceAdapter` exists for exactly this pairing.
    SnpEff predicts the molecular effect and is the primary; the MANE join places
    a variant in a gene without claiming an effect, and backfills the variants
    SnpEff leaves unannotated. Neither alone is sufficient — SnpEff can return
    nothing for a contig its database does not know, and MANE can never say
    whether a variant is damaging — and a gene assignment is what
    ``VariantRecord.gene_symbols`` needs for pairing to produce any candidate at
    all.
    """
    pins = SnpEffArtifactPins.from_manifest(
        resolved.snpeff_pins, genome_database=resolved.snpeff_genome
    )
    primary: ConsequenceAdapter = SnpEffConsequenceAdapter(
        jar_path=resolved.snpeff_jar,
        data_dir=resolved.snpeff_data_dir,
        genome_database=resolved.snpeff_genome,
        pins=pins,
        java_binary=resolved.snpeff_java_binary,
        config_path=resolved.snpeff_config,
        build=resolved.build,
        mane_summary=resolved.mane_summary,
    )
    fallback: ConsequenceAdapter = ManeGeneAdapter(
        ManeGeneIndex(
            resolved.mane_gtf,
            resolved.mane_summary,
            expected_gtf_sha256=resolved.mane_gtf_sha256,
            expected_summary_sha256=resolved.mane_summary_sha256,
            build=resolved.build,
        )
    )
    return GeneBackfillConsequenceAdapter(primary, fallback)


__all__ = [
    "AdapterBindingError",
    "BoundAdapters",
    "PartialCoverageFrequencyAdapter",
    "ResolvedResources",
    "build_clinical_adapter",
    "build_consequence_adapter",
    "build_frequency_adapter",
    "build_real_adapter_set",
    "representation_warnings",
    "resolve_real_resources",
]
