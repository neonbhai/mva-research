"""Unit tests for the out-of-repo reference-release layer (`mva.resources`).

Nothing here touches the real 202.8 GB resource root, the network, or
``$MVA_RESOURCES``. Every test builds its own tiny root under ``tmp_path``, which
is the same arrangement the privacy model requires of a real run and is what keeps
these runnable on a machine that has never downloaded gnomAD.

What is being held to account, in order of how expensive the mistake would be:

* the resource root is never guessed (ADR 0020);
* the sampling plan is deterministic, total and size-bound (GP-30), so a committed
  ``spot_sha256`` means the same thing on every machine and in every process;
* verification fails CLOSED and names the file without echoing it (PRIV-09);
* the manifest's write model and read model cannot drift apart.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
from pydantic import ValidationError

from mva.config import CaseConfig, ResourceSettings, Workspace
from mva.resources import (
    SPOT_EXHAUSTIVE_BELOW,
    SPOT_PLAN,
    FormatCheck,
    IndexCheck,
    IntegrityMode,
    IntegrityRecord,
    ReferenceResource,
    ResourceCheck,
    ResourceError,
    ResourceIntegrityError,
    ResourceKind,
    ResourceRoot,
    ResourceRootError,
    ResourceStatus,
    assert_resources_verified,
    load_resource_manifest,
    reference_fasta_path,
    required_resources,
    resolve_resource_root,
    spot_digest,
    spot_windows,
    tree_digest,
    verify_resource,
)

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
COMMITTED_MANIFEST = REPO_ROOT / "knowledge" / "manifests" / "resources.yaml"


def _root(tmp_path: Path) -> ResourceRoot:
    resources = tmp_path / "resources"
    resources.mkdir(exist_ok=True)
    return ResourceRoot(root=resources.resolve(), repo_root=(tmp_path / "repo").resolve())


def _fetched(target: Path, root: ResourceRoot, **overrides: object) -> ReferenceResource:
    """A registered entry pinned to whatever ``target`` currently contains."""
    from mva.determinism import hash_file

    fields: dict[str, object] = {
        "name": "example",
        "version": "1.0",
        "path": target.relative_to(root.root).as_posix(),
        "description": "A registered example resource.",
        "status": ResourceStatus.FETCHED,
        "sha256": hash_file(target),
        "size_bytes": target.stat().st_size,
        "retrieved": "2026-08-28",
        "integrity": IntegrityRecord(
            verified_at="2026-08-28",
            spot_plan=SPOT_PLAN,
            spot_sha256=spot_digest(target),
            format_check=FormatCheck.PLAIN_TEXT,
            format_detail="2-column table",
        ),
    }
    fields.update(overrides)
    return ReferenceResource(**fields)  # type: ignore[arg-type]


# --------------------------------------------------------------------------- root


class TestResolveResourceRoot:
    def test_absent_configuration_raises_rather_than_guessing(self, tmp_path: Path) -> None:
        """ADR 0020: the hard-coded fallback to one contributor's home directory is
        gone. A guessed root that is wrong fails later, somewhere less obvious --
        'resource missing' and 'resource not registered' are indistinguishable by
        the time either reaches an adapter."""
        with pytest.raises(ResourceRootError, match="No resource root configured"):
            resolve_resource_root(env={}, repo_root=tmp_path)

    def test_env_var_is_used_when_no_explicit_path(self, tmp_path: Path) -> None:
        outside = tmp_path / "outside" / "resources"
        outside.mkdir(parents=True)
        repo = tmp_path / "repo"
        repo.mkdir()
        root = resolve_resource_root(env={"MVA_RESOURCES": str(outside)}, repo_root=repo)
        assert root.root == outside.resolve()

    def test_explicit_argument_beats_the_env_var(self, tmp_path: Path) -> None:
        explicit = tmp_path / "explicit"
        explicit.mkdir()
        env_dir = tmp_path / "from_env"
        env_dir.mkdir()
        repo = tmp_path / "repo"
        repo.mkdir()
        root = resolve_resource_root(explicit, env={"MVA_RESOURCES": str(env_dir)}, repo_root=repo)
        assert root.root == explicit.resolve()

    def test_root_inside_the_repo_is_refused(self, tmp_path: Path) -> None:
        repo = tmp_path / "repo"
        (repo / "resources").mkdir(parents=True)
        with pytest.raises(ResourceRootError, match="inside the repository"):
            resolve_resource_root(repo / "resources", env={}, repo_root=repo)

    def test_missing_root_is_refused(self, tmp_path: Path) -> None:
        repo = tmp_path / "repo"
        repo.mkdir()
        with pytest.raises(ResourceRootError, match="does not exist"):
            resolve_resource_root(tmp_path / "nope", env={}, repo_root=repo)

    def test_path_escape_is_refused(self, tmp_path: Path) -> None:
        root = _root(tmp_path)
        with pytest.raises(ResourceRootError, match="escapes"):
            root.path("../elsewhere/secret")


# --------------------------------------------------------------------------- sampling plan


class TestSpotPlan:
    def test_small_files_are_covered_completely(self) -> None:
        """At or below the sample size, SPOT and FULL read the same bytes. Asserted
        rather than left implicit: 'the cheap check happens to be exhaustive here'
        is a property downstream reasoning is allowed to rely on."""
        for size in (1, 1024, SPOT_EXHAUSTIVE_BELOW):
            assert spot_windows(size) == ((0, size),)

    def test_large_file_windows_are_ordered_in_range_and_bounded(self) -> None:
        size = 18_789_775_480  # the real gnomAD chr1 shard
        windows = spot_windows(size)
        offsets = [offset for offset, _ in windows]
        assert offsets == sorted(offsets), "windows must be read in file order"
        assert all(offset + length <= size for offset, length in windows)
        assert sum(length for _, length in windows) <= SPOT_EXHAUSTIVE_BELOW
        assert windows[0][0] == 0, "the head must be sampled"
        assert windows[-1][0] + windows[-1][1] == size, "the tail must reach the last byte"

    def test_plan_is_deterministic(self) -> None:
        assert spot_windows(10**10) == spot_windows(10**10)

    def test_digest_is_bound_to_size_not_only_content(self, tmp_path: Path) -> None:
        """The failure this guards is the real one: an interrupted download whose
        every surviving byte is correct. Without the size bound, a file truncated
        between two sample windows would digest identically."""
        target = tmp_path / "payload.bin"
        target.write_bytes(b"A" * 4096)
        before = spot_digest(target)
        target.write_bytes(b"A" * 4095)
        assert spot_digest(target) != before

    def test_digest_changes_when_a_sampled_byte_changes(self, tmp_path: Path) -> None:
        target = tmp_path / "payload.bin"
        target.write_bytes(b"A" * 4096)
        before = spot_digest(target)
        target.write_bytes(b"B" + b"A" * 4095)
        assert spot_digest(target) != before

    def test_digest_is_stable_across_calls(self, tmp_path: Path) -> None:
        target = tmp_path / "payload.bin"
        target.write_bytes(bytes(range(256)) * 64)
        assert spot_digest(target) == spot_digest(target)


class TestTreeDigest:
    def _tree(self, base: Path) -> Path:
        base.mkdir(parents=True)
        (base / "a.bin").write_bytes(b"alpha")
        (base / "nested").mkdir()
        (base / "nested" / "b.bin").write_bytes(b"beta")
        return base

    def test_same_contents_digest_the_same(self, tmp_path: Path) -> None:
        one = self._tree(tmp_path / "one")
        two = self._tree(tmp_path / "two")
        assert tree_digest(one) == tree_digest(two)

    def test_editing_any_member_changes_the_digest(self, tmp_path: Path) -> None:
        tree = self._tree(tmp_path / "db")
        before = tree_digest(tree)
        (tree / "nested" / "b.bin").write_bytes(b"gamma")
        assert tree_digest(tree) != before

    def test_adding_a_member_changes_the_digest(self, tmp_path: Path) -> None:
        """The whole reason the SnpEff database is pinned as a tree: a per-file pin
        on snpEffectPredictor.bin alone would not notice the sequence bins."""
        tree = self._tree(tmp_path / "db")
        before = tree_digest(tree)
        (tree / "c.bin").write_bytes(b"delta")
        assert tree_digest(tree) != before


# --------------------------------------------------------------------------- the read model


class TestReferenceResource:
    def test_absolute_path_is_refused(self) -> None:
        with pytest.raises(ValidationError):
            ReferenceResource(name="x", version="1", path="/etc/passwd", description="absolute")

    def test_parent_escape_is_refused(self) -> None:
        with pytest.raises(ValidationError):
            ReferenceResource(name="x", version="1", path="../outside.txt", description="escape")

    def test_malformed_sha256_is_refused(self) -> None:
        with pytest.raises(ValidationError):
            ReferenceResource(
                name="x",
                version="1",
                path="a.txt",
                description="bad hash",
                status=ResourceStatus.FETCHED,
                sha256="NOTAHASH",
                size_bytes=1,
                retrieved="2026-08-28",
            )

    def test_fetched_without_a_hash_is_refused(self) -> None:
        with pytest.raises(ValidationError, match="missing sha256"):
            ReferenceResource(
                name="x",
                version="1",
                path="a.txt",
                description="claimed fetched",
                status=ResourceStatus.FETCHED,
            )

    def test_unfetched_with_a_hash_is_refused(self) -> None:
        """An integrity claim that cannot be checked is worse than an honest gap."""
        with pytest.raises(ValidationError, match="may only be set once"):
            ReferenceResource(
                name="x",
                version="1",
                path="a.txt",
                description="claimed hash while unfetched",
                sha256="a" * 64,
                size_bytes=1,
                retrieved="2026-08-28",
            )


@pytest.mark.unit
def test_acquisition_entry_is_the_same_schema_as_the_runtime_read_model() -> None:
    """The write model must remain a strict extension of the read model.

    ``tools/acquire`` writes ``resources.yaml`` and ``mva.resources`` reads it. Two
    hand-maintained copies of that schema would drift, and the drift would surface
    as a field silently ignored on the read side (``extra="ignore"`` there is
    deliberate; see ``ReferenceResource``). Subclassing makes drift impossible, and
    this test is what stops someone from un-subclassing it later.
    """
    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))
    from tools.acquire.models import ResourceEntry

    assert issubclass(ResourceEntry, ReferenceResource)
    missing = set(ReferenceResource.model_fields) - set(ResourceEntry.model_fields)
    assert not missing, (
        f"ResourceEntry has lost read-model field(s): {sorted(missing)}.\n\n"
        "Remediation: ResourceEntry must stay a subclass of ReferenceResource. The "
        "runtime reads what this tool writes; a field the writer drops is a field the "
        "reader silently defaults."
    )
    assert ResourceEntry.model_config.get("extra") == "forbid", (
        "ResourceEntry must keep extra='forbid'. The read model ignores unknown keys "
        "on purpose, so the WRITER is the only place a typo in resources.yaml can be "
        "caught -- which is also the only place it can be fixed."
    )


# --------------------------------------------------------------------------- verification


class TestVerifyResource:
    def test_matching_bytes_pass_in_both_modes(self, tmp_path: Path) -> None:
        root = _root(tmp_path)
        target = root.root / "table.tsv"
        target.write_text("a\tb\n1\t2\n", encoding="utf-8")
        entry = _fetched(target, root)
        for mode in (IntegrityMode.SPOT, IntegrityMode.FULL):
            assert verify_resource(root, entry, mode=mode).check is ResourceCheck.OK

    def test_unfetched_entry_is_unpinned_not_a_failure(self, tmp_path: Path) -> None:
        root = _root(tmp_path)
        entry = ReferenceResource(
            name="pending", version="1", path="absent.txt", description="not yet fetched"
        )
        result = verify_resource(root, entry)
        assert result.check is ResourceCheck.UNPINNED
        assert result.ok

    def test_pinned_but_absent_is_missing(self, tmp_path: Path) -> None:
        root = _root(tmp_path)
        target = root.root / "table.tsv"
        target.write_text("a\tb\n", encoding="utf-8")
        entry = _fetched(target, root)
        target.unlink()
        assert verify_resource(root, entry).check is ResourceCheck.MISSING

    def test_truncation_is_reported_as_a_size_mismatch(self, tmp_path: Path) -> None:
        """Size is checked before any digest because it is two stat calls and it is
        the failure that actually happens -- and because 'N bytes short' is a far
        more actionable message than 'digest differs'."""
        root = _root(tmp_path)
        target = root.root / "big.bin"
        target.write_bytes(b"x" * 5000)
        entry = _fetched(target, root)
        target.write_bytes(b"x" * 4000)
        result = verify_resource(root, entry)
        assert result.check is ResourceCheck.SIZE_MISMATCH
        assert "1000 bytes short" in result.message

    def test_same_size_different_bytes_is_a_digest_mismatch(self, tmp_path: Path) -> None:
        root = _root(tmp_path)
        target = root.root / "same_size.bin"
        target.write_bytes(b"original content")
        entry = _fetched(target, root)
        target.write_bytes(b"tampered content")
        for mode in (IntegrityMode.SPOT, IntegrityMode.FULL):
            assert verify_resource(root, entry, mode=mode).check is ResourceCheck.DIGEST_MISMATCH

    def test_unknown_sampling_plan_is_refused_not_compared(self, tmp_path: Path) -> None:
        """Digests computed under different plans are not comparable. Reporting a
        difference between them as corruption would be a false alarm; reporting it
        as a pass would be a lie."""
        root = _root(tmp_path)
        target = root.root / "table.tsv"
        target.write_text("a\tb\n", encoding="utf-8")
        entry = _fetched(target, root)
        stale = entry.model_copy(
            update={
                "integrity": entry.integrity.model_copy(  # type: ignore[union-attr]
                    update={"spot_plan": "spot-v0:head=1MiB"}
                )
            }
        )
        result = verify_resource(root, stale, mode=IntegrityMode.SPOT)
        assert result.check is ResourceCheck.PLAN_UNKNOWN
        assert "not comparable" in result.message
        # FULL does not depend on the plan at all, so it still decides the question.
        assert verify_resource(root, stale, mode=IntegrityMode.FULL).check is ResourceCheck.OK

    def test_entry_without_an_integrity_record_cannot_be_spot_checked(self, tmp_path: Path) -> None:
        root = _root(tmp_path)
        target = root.root / "table.tsv"
        target.write_text("a\tb\n", encoding="utf-8")
        entry = _fetched(target, root, integrity=None)
        assert (
            verify_resource(root, entry, mode=IntegrityMode.SPOT).check
            is ResourceCheck.PLAN_UNKNOWN
        )

    def test_directory_resource_is_verified_as_a_tree(self, tmp_path: Path) -> None:
        root = _root(tmp_path)
        database = root.root / "db"
        database.mkdir()
        (database / "predictor.bin").write_bytes(b"model")
        (database / "sequence.1.bin").write_bytes(b"acgt")

        entry = ReferenceResource(
            name="db",
            version="GRCh38.115",
            path="db",
            description="A multi-file database.",
            kind=ResourceKind.DIRECTORY,
            status=ResourceStatus.FETCHED,
            sha256=tree_digest(database),
            size_bytes=len(b"model") + len(b"acgt"),
            retrieved="2026-08-28",
        )
        assert verify_resource(root, entry, mode=IntegrityMode.FULL).check is ResourceCheck.OK

        (database / "sequence.2.bin").write_bytes(b"tgca")
        assert (
            verify_resource(root, entry, mode=IntegrityMode.FULL).check
            is ResourceCheck.SIZE_MISMATCH
        )


class TestAssertResourcesVerified:
    def test_failure_raises_and_names_the_resource(self, tmp_path: Path) -> None:
        root = _root(tmp_path)
        target = root.root / "table.tsv"
        target.write_text("a\tb\n", encoding="utf-8")
        entry = _fetched(target, root, name="clinvar_vcf")
        target.unlink()
        with pytest.raises(ResourceIntegrityError, match="clinvar_vcf"):
            assert_resources_verified(root, [entry])

    def test_failure_message_never_echoes_file_contents(self, tmp_path: Path) -> None:
        """PRIV-09. A digest is a one-way summary and is safe to print in full; the
        bytes that produced it are not."""
        record = "PROBAND01\tchr15\t40200000\tC\tT"
        root = _root(tmp_path)
        target = root.root / "table.tsv"
        target.write_text(record, encoding="utf-8")
        entry = _fetched(target, root)
        target.write_text(record.replace("40200000", "40200001"), encoding="utf-8")
        with pytest.raises(ResourceIntegrityError) as excinfo:
            assert_resources_verified(root, [entry])
        message = str(excinfo.value)
        assert "40200000" not in message
        assert "PROBAND01" not in message
        assert entry.path in message

    def test_unpinned_entries_do_not_fail_the_gate(self, tmp_path: Path) -> None:
        root = _root(tmp_path)
        pending = ReferenceResource(
            name="pending", version="1", path="absent.txt", description="not yet fetched"
        )
        assert assert_resources_verified(root, [pending])[0].check is ResourceCheck.UNPINNED


# --------------------------------------------------------------------------- the committed manifest


class TestCommittedManifest:
    def test_it_parses_with_the_runtime_reader(self) -> None:
        manifest = load_resource_manifest(COMMITTED_MANIFEST)
        assert manifest.resources, "the committed manifest registers no resources"
        assert manifest.paths_relative_to == "resource_root"

    def test_every_pinned_entry_records_what_was_verified(self) -> None:
        """ADR 0020's central claim. A hash without a record of what was checked
        invites the reader to assume more was done than was."""
        manifest = load_resource_manifest(COMMITTED_MANIFEST)
        naked = [
            entry.name
            for entry in manifest.entries()
            if entry.status is ResourceStatus.FETCHED and entry.integrity is None
        ]
        assert not naked, (
            "Pinned resources with no integrity record: "
            + ", ".join(naked)
            + "\n\nRemediation: regenerate with `uv run python -m tools.acquire "
            "write-manifest`. A sha256 alone does not say whether anything ever opened "
            "the file as the format it claims to be (ADR 0020)."
        )

    def test_no_pinned_entry_admits_an_unchecked_format(self) -> None:
        manifest = load_resource_manifest(COMMITTED_MANIFEST)
        unchecked = [
            entry.name
            for entry in manifest.entries()
            if entry.integrity is not None
            and entry.integrity.format_check is FormatCheck.NOT_CHECKED
        ]
        assert not unchecked, (
            "Pinned resources whose format was never verified: "
            + ", ".join(unchecked)
            + "\n\nNOT_CHECKED is the honest default, never a passing state. Either add "
            "a probe for that format in tools/acquire/formats.py, or accept the resource "
            "as OPAQUE_BINARY deliberately."
        )

    def test_every_entry_is_real_reference_data(self) -> None:
        """GP-20: this registry holds real public releases. Synthetic demo tables
        belong in knowledge/manifests/knowledge.yaml, and mixing the two is how a
        fabricated frequency ends up described as gnomAD."""
        manifest = load_resource_manifest(COMMITTED_MANIFEST)
        assert all(not entry.synthetic for entry in manifest.entries())

    def test_settings_defaults_point_at_registered_resources(self) -> None:
        """The config defaults and the manifest are two descriptions of one layout.
        When they disagree, the pipeline opens an unpinned path and nothing says so."""
        manifest = load_resource_manifest(COMMITTED_MANIFEST)
        registered = {entry.path for entry in manifest.entries()}
        settings = ResourceSettings()
        expected = {
            settings.reference_fasta,
            settings.reference_fasta_index,
            settings.clinvar_vcf,
            settings.clinvar_vcf_index,
            settings.gnomad_constraint_metrics,
            settings.mane_summary,
            settings.mane_gtf,
            settings.hpo_obo,
            settings.hpo_annotations,
            settings.hpo_genes_to_phenotype,
            settings.hpo_phenotype_to_genes,
            settings.snpeff_jar,
            settings.snpeff_config,
            settings.gnomad_shard("chr21"),
            f"{settings.snpeff_data_dir}/{settings.snpeff_genome}",
        }
        unregistered = sorted(expected - registered)
        assert not unregistered, (
            "ResourceSettings names paths that no manifest entry pins:\n"
            + "\n".join(f"  {path}" for path in unregistered)
            + "\n\nRemediation: register them in tools/acquire/catalog.py and re-run "
            "`uv run python -m tools.acquire write-manifest`, or correct the default in "
            "ResourceSettings. A configured path with no pin behind it is read as data "
            "with no integrity check at all."
        )

    def test_required_resources_refuses_an_unfetched_dependency(self) -> None:
        manifest = load_resource_manifest(COMMITTED_MANIFEST)
        with pytest.raises(ResourceError, match="No resource named"):
            required_resources(manifest, ["not_a_registered_resource"])


# --------------------------------------------------------------------------- reference FASTA


class TestReferenceFastaPath:
    def _config(self, **inputs: object) -> CaseConfig:
        return CaseConfig.model_validate(
            {
                "case_id": "test-case",
                "proband_id": "PROBAND01",
                "genome_build": "GRCh38",
                "synthetic": True,
                "inputs": {"vcf": "in.vcf", "phenotype": "p.tsv", **inputs},
            }
        )

    def test_resource_root_supplies_the_shared_release(self, tmp_path: Path) -> None:
        root = _root(tmp_path)
        resolved = reference_fasta_path(self._config(), resource_root=root)
        assert resolved == root.root / "reference" / "GRCh38_no_alt.fa"

    def test_case_level_override_wins(self, tmp_path: Path) -> None:
        root = _root(tmp_path)
        workspace = Workspace(root=tmp_path.resolve(), repo_root=(tmp_path / "repo").resolve())
        resolved = reference_fasta_path(
            self._config(reference_fasta="refs/custom.fa"),
            workspace=workspace,
            resource_root=root,
        )
        assert resolved == (tmp_path / "refs" / "custom.fa").resolve()

    def test_no_reference_configured_returns_none_rather_than_raising(self) -> None:
        """GP-14 applied to configuration: normalising without a reference is a
        DEGRADED mode, not an impossible one. The caller decides whether to proceed
        and must surface `representation_limitation` if it does -- a function that
        raised here would make that an import-time accident instead of a decision."""
        assert reference_fasta_path(self._config()) is None

    def test_override_without_a_workspace_is_a_loud_error(self) -> None:
        with pytest.raises(ResourceError, match="workspace-relative"):
            reference_fasta_path(self._config(reference_fasta="refs/custom.fa"))


# --------------------------------------------------------------------------- config surface


class TestResourceSettings:
    def test_absolute_resource_path_is_refused(self) -> None:
        with pytest.raises(ValidationError):
            ResourceSettings(reference_fasta="/opt/genomes/GRCh38.fa")

    def test_parent_escape_is_refused(self) -> None:
        with pytest.raises(ValidationError):
            ResourceSettings(clinvar_vcf="../../etc/passwd")

    def test_integrity_mode_is_a_closed_set(self) -> None:
        """There is deliberately no 'off' (ADR 0020): a mode that skips verification
        is the one selected under deadline pressure."""
        assert ResourceSettings().integrity_mode == IntegrityMode.SPOT.value
        assert ResourceSettings(integrity_mode="full").integrity_mode == IntegrityMode.FULL.value
        with pytest.raises(ValidationError):
            ResourceSettings(integrity_mode="off")

    def test_gnomad_shard_expands_the_template(self) -> None:
        settings = ResourceSettings()
        assert settings.gnomad_shard("chr21") == (
            "gnomad/v4.1_exomes/gnomad.exomes.v4.1.sites.chr21.vcf.bgz"
        )

    def test_case_config_carries_resource_settings_by_default(self) -> None:
        config = CaseConfig.model_validate(
            {
                "case_id": "test-case",
                "proband_id": "PROBAND01",
                "genome_build": "GRCh38",
                "synthetic": True,
                "inputs": {"vcf": "in.vcf", "phenotype": "p.tsv"},
            }
        )
        assert config.resources.snpeff_genome == "GRCh38.115"
        assert config.inputs.reference_fasta is None

    def test_resource_settings_change_the_config_hash(self) -> None:
        """GP-30/GP-31: which reference release a run used is part of what makes the
        run reproducible, so it has to be inside the hash stamped on the manifest."""
        base = {
            "case_id": "test-case",
            "proband_id": "PROBAND01",
            "genome_build": "GRCh38",
            "synthetic": True,
            "inputs": {"vcf": "in.vcf", "phenotype": "p.tsv"},
        }
        default = CaseConfig.model_validate(base)
        altered = CaseConfig.model_validate({**base, "resources": {"snpeff_genome": "GRCh38.105"}})
        assert default.config_hash() != altered.config_hash()


# --------------------------------------------------------------------------- stale indexes


class TestIndexStaleness:
    """A tabix index and its VCF are two separate downloads and can disagree.

    htslib warns whenever the index's mtime precedes the data file's, which is true
    for all 25 gnomAD shards here purely because each small `.tbi` finished
    downloading before its multi-gigabyte `.bgz` did. That warning is noise; the
    failure it gestures at is not. These tests hold the real check to account.
    """

    FIXTURE = (
        REPO_ROOT / "tests" / "fixtures" / "gnomad" / "gnomad.exomes.v4.1.sites.chr21.slice.vcf.bgz"
    )

    def _pair(self, tmp_path: Path) -> Path:
        pytest.importorskip("pysam", reason="index cross-checks need the 'genomics' extra")
        data = tmp_path / self.FIXTURE.name
        data.write_bytes(self.FIXTURE.read_bytes())
        index = tmp_path / f"{self.FIXTURE.name}.tbi"
        index.write_bytes(Path(f"{self.FIXTURE}.tbi").read_bytes())
        return data

    def test_a_matching_pair_is_proven_consistent(self, tmp_path: Path) -> None:
        from tools.acquire.formats import probe_format

        result = probe_format(self._pair(tmp_path))
        assert result.ok
        assert result.index_check is IndexCheck.CONSISTENT
        assert "deep fetch" in result.index_detail

    def test_mtime_skew_alone_does_not_change_the_verdict(self, tmp_path: Path) -> None:
        """The exact situation on disk for every gnomAD shard. An index that is
        older than its data is not thereby wrong, and this project will not record
        a timestamp comparison as if it were evidence (ADR 0020)."""
        import os

        from tools.acquire.formats import probe_format

        data = self._pair(tmp_path)
        index = Path(f"{data}.tbi")
        os.utime(index, (1, 1))  # index far older than the data file

        result = probe_format(data)
        assert result.index_check is IndexCheck.CONSISTENT, (
            "an index older than its data must still be judged by access, not by mtime"
        )
        assert "mtime" in result.index_detail, (
            "the skew should be REPORTED -- silently ignoring it hides the fact that a "
            "human reading htslib's warning was right to ask"
        )

    def test_plain_truncation_is_caught_before_anything_can_be_pinned(self, tmp_path: Path) -> None:
        """A half-written BGZF file has no end-of-file block. That is diagnosed
        first, and more precisely, than any index question: the stream was never
        finished, so there is nothing to have an opinion about the index of."""
        from tools.acquire.formats import probe_format

        data = self._pair(tmp_path)
        payload = data.read_bytes()
        data.write_bytes(payload[: len(payload) // 2])

        result = probe_format(data)
        assert not result.ok
        assert result.problem is not None
        assert "never finished" in result.problem

    def test_index_referencing_bytes_past_eof_is_stale_and_fails_closed(
        self, tmp_path: Path
    ) -> None:
        """The structural line of evidence, isolated.

        Content is removed from the middle and the mandatory BGZF end-of-file block
        is restored, so the file passes every check that looks only at the stream:
        the magic bytes are right, the EOF block is present, the size is
        self-consistent. Only comparing the index's deepest reference against the
        file's actual length reveals that this index was built from something
        larger. That is the case a hash cannot catch either -- the hash of a wrong
        pair is a perfectly valid hash.
        """
        from tools.acquire.formats import _BGZF_EOF, probe_format

        data = self._pair(tmp_path)
        payload = data.read_bytes()
        data.write_bytes(payload[: len(payload) // 2] + _BGZF_EOF)

        result = probe_format(data)
        assert not result.ok
        assert result.index_check is IndexCheck.STALE
        assert result.problem is not None
        assert "beyond the end" in result.problem

    def test_out_of_window_records_are_reported_as_a_mismatch(self) -> None:
        """The behavioural line of evidence, exercised directly.

        A stale index of the SAME size passes the structural check, so the probe
        also seeks deep and insists the records that come back are inside the
        window it asked for. Driven with a stub here rather than a crafted file
        pair: the property under test is the decision, not htslib's seeking.
        """
        from tools.acquire.formats import _deep_fetch_probe, _TabixIndex

        class _Record:
            def __init__(self, pos: int) -> None:
                self.pos = pos

        class _Wrong:
            def fetch(self, _contig: str, _start: int, _end: int) -> list[_Record]:
                return [_Record(10)]  # nowhere near the requested deep window

        parsed = _TabixIndex(["chr21"], 0, {"chr21": 100})
        probed = _deep_fetch_probe(parsed, ["chr21"], _Wrong())
        assert probed is not None
        ok, detail = probed
        assert ok is False
        assert "out-of-window" in detail

    def test_in_window_records_are_accepted(self) -> None:
        from tools.acquire.formats import _deep_fetch_probe, _TabixIndex

        class _Record:
            def __init__(self, pos: int) -> None:
                self.pos = pos

        start = (100 - 1) << 14

        class _Right:
            def fetch(self, _contig: str, _start: int, _end: int) -> list[_Record]:
                return [_Record(start + 1)]

        parsed = _TabixIndex(["chr21"], 0, {"chr21": 100})
        probed = _deep_fetch_probe(parsed, ["chr21"], _Right())
        assert probed is not None
        assert probed[0] is True


@pytest.mark.unit
def test_every_indexed_release_records_its_index_verdict() -> None:
    """ADR 0020: the manifest must say whether each index was checked against its
    data file, not merely that the data file hashed correctly.

    Both halves matter. A pinned VCF whose index is stale returns the wrong records
    for exactly the rare coordinates this pipeline exists to find, and every hash in
    the manifest still verifies while it does.
    """
    manifest = load_resource_manifest(COMMITTED_MANIFEST)
    indexed = [
        entry
        for entry in manifest.entries()
        if entry.integrity is not None
        and entry.integrity.format_check
        in {FormatCheck.BGZF_TABIX_INDEXED, FormatCheck.FASTA_FAIDX_CONSISTENT}
    ]
    assert indexed, "no indexed releases found in the manifest"

    unproven = [
        entry.name
        for entry in indexed
        if entry.integrity is None or entry.integrity.index_check is not IndexCheck.CONSISTENT
    ]
    assert not unproven, (
        "Indexed releases whose index was never proven to describe their data: "
        + ", ".join(unproven)
        + "\n\nRemediation: regenerate with `uv run python -m tools.acquire write-manifest`. "
        "If one reports STALE, re-fetch the pair -- do not pin them together."
    )
