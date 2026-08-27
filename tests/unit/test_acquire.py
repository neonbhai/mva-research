"""Unit tests for the public-reference acquisition tool (``tools/acquire``).

``tools/acquire`` lives outside ``src/mva`` on purpose: it is the one place in this
project allowed to touch the network, in exchange for a structural promise that it
will only ever talk to a small allowlist of public reference hosts and will never
see, log or transmit anything patient-derived. See
``knowledge/adapters/README.md``, "Why no adapter here may touch the network
(PRIV-05)", for the contract these tests hold the tool to.

No test here makes a real network call. ``fetch_resource``'s transport
(``tools.acquire.fetch._http_download``) is monkeypatched everywhere it would
otherwise open a socket, so what is actually exercised is the safety machinery
around it: the host allowlist, the outside-the-repo assertion, the growth/format
sanity checks, and the manifest round trip -- exactly the properties that matter
if this tool is ever pointed at something it shouldn't be.
"""

from __future__ import annotations

import gzip
import inspect
import sys
from collections.abc import Callable
from pathlib import Path

# tools/ is not an installed package (unlike `mva`, which is editable-installed via
# `uv sync`), so it is only importable once the repo root is on sys.path. Bootstrap
# that here rather than relying on pytest's own import-mode path insertion, which
# only guarantees `tests/unit` (this file's own directory) is on the path.
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import pytest  # noqa: E402 - must follow the sys.path bootstrap above
import tools.acquire.fetch as fetch_module  # noqa: E402
from pydantic import ValidationError  # noqa: E402
from tools.acquire import (  # noqa: E402
    ALLOWED_HOSTS,
    KNOWN_RESOURCES,
    AcquisitionError,
    DisallowedHostError,
    ResourceEntry,
    ResourceFetchError,
    ResourceRootError,
    ResourceStatus,
    ResourceVerificationError,
    VerificationStatus,
    assert_allowed_host,
    assert_verified,
    fetch_resource,
    is_download_stable,
    load_resources_manifest,
    render_resources_yaml,
    resolve_resource_root,
    sniff_content_mismatch,
    survey_all,
    survey_resource,
    verify_all,
    verify_resource,
    write_resources_manifest,
)

pytestmark = pytest.mark.unit

_HEX64 = "a" * 64


def _make_entry(**overrides: object) -> ResourceEntry:
    """A minimal, valid, allowlisted, not-yet-fetched resource entry."""
    fields: dict[str, object] = {
        "name": "test_resource",
        "url": "https://ftp.ncbi.nlm.nih.gov/pub/some/file.txt",
        "version": "1.0",
        "path": "somewhere/file.txt",
        "license": "Public domain",
        "description": "A test resource.",
    }
    fields.update(overrides)
    return ResourceEntry(**fields)  # type: ignore[arg-type]


# --------------------------------------------------------------------------- hosts


class TestHostAllowlist:
    def test_every_allowlisted_host_is_accepted(self) -> None:
        for host in ALLOWED_HOSTS:
            assert_allowed_host(f"https://{host}/some/path")  # must not raise

    def test_unknown_host_is_rejected(self) -> None:
        with pytest.raises(DisallowedHostError, match="not on the public-reference allowlist"):
            assert_allowed_host("https://evil.example.com/clinvar.vcf.gz")

    def test_plain_http_is_rejected_even_for_an_allowlisted_host(self) -> None:
        with pytest.raises(DisallowedHostError, match="https"):
            assert_allowed_host("http://ftp.ncbi.nlm.nih.gov/pub/x")

    def test_non_http_scheme_is_rejected(self) -> None:
        with pytest.raises(DisallowedHostError):
            assert_allowed_host("ftp://ftp.ncbi.nlm.nih.gov/pub/x")

    def test_host_matching_is_case_insensitive(self) -> None:
        assert_allowed_host("https://FTP.NCBI.NLM.NIH.GOV/pub/x")

    def test_allowlist_message_never_echoes_query_string(self) -> None:
        """The rejection message names the host, not the full URL -- a query string
        is exactly where a careless caller could put patient-shaped data."""
        with pytest.raises(DisallowedHostError) as excinfo:
            assert_allowed_host("https://evil.example.com/?patient=secret-value")
        assert "secret-value" not in str(excinfo.value)


# --------------------------------------------------------------------------- ResourceEntry


class TestResourceEntry:
    def test_valid_entry_constructs(self) -> None:
        entry = _make_entry()
        assert entry.status is ResourceStatus.NOT_FETCHED
        assert entry.synthetic is False
        assert entry.sha256 is None

    def test_construction_rejects_disallowed_host(self) -> None:
        """Defence in depth: a bad host can never even enter the registry, not only
        at fetch time. Raised directly (not wrapped in ValidationError) because
        pydantic only wraps ValueError/AssertionError from a validator."""
        with pytest.raises(DisallowedHostError):
            _make_entry(url="https://evil.example.com/x")

    def test_absolute_path_is_rejected(self) -> None:
        with pytest.raises(ValidationError, match="absolute path"):
            _make_entry(path="/etc/passwd")

    def test_path_escape_is_rejected(self) -> None:
        with pytest.raises(ValidationError, match="escaping the resource root"):
            _make_entry(path="../../etc/passwd")

    def test_malformed_sha256_is_rejected(self) -> None:
        with pytest.raises(ValidationError, match="malformed sha256"):
            _make_entry(sha256="not-a-hash")

    @pytest.mark.parametrize("missing_field", ["sha256", "size_bytes", "retrieved"])
    def test_fetched_status_requires_all_three_fields(self, missing_field: str) -> None:
        fields: dict[str, object] = {
            "status": ResourceStatus.FETCHED,
            "sha256": _HEX64,
            "size_bytes": 10,
            "retrieved": "2026-08-27",
        }
        fields[missing_field] = None
        with pytest.raises(ValidationError, match="missing sha256, size_bytes or retrieved"):
            _make_entry(**fields)

    def test_fetched_status_with_all_fields_is_valid(self) -> None:
        entry = _make_entry(
            status=ResourceStatus.FETCHED, sha256=_HEX64, size_bytes=10, retrieved="2026-08-27"
        )
        assert entry.status is ResourceStatus.FETCHED

    @pytest.mark.parametrize(
        "extra", [{"sha256": _HEX64}, {"size_bytes": 10}, {"retrieved": "2026-08-27"}]
    )
    def test_not_fetched_status_rejects_any_fetched_field(self, extra: dict[str, object]) -> None:
        with pytest.raises(ValidationError, match="may only be set once the resource is"):
            _make_entry(status=ResourceStatus.NOT_FETCHED, **extra)

    def test_extra_field_is_rejected(self) -> None:
        with pytest.raises(ValidationError):
            _make_entry(made_up_field="surprise")

    def test_entry_is_frozen(self) -> None:
        entry = _make_entry()
        with pytest.raises(ValidationError):
            entry.name = "changed"  # type: ignore[misc]


# --------------------------------------------------------------------------- resolve_resource_root


class TestResolveResourceRoot:
    def test_default_falls_back_to_documented_path(self, tmp_path: Path) -> None:
        # repo_root injected so this never depends on where the real repo lives.
        root = resolve_resource_root(env={}, repo_root=tmp_path)
        assert root == Path("~/Contri/bio-hackathon/mva-resources").expanduser()

    def test_env_var_overrides_default(self, tmp_path: Path) -> None:
        outside = tmp_path / "outside" / "resources"
        outside.mkdir(parents=True)
        repo = tmp_path / "repo"
        repo.mkdir()
        root = resolve_resource_root(env={"MVA_RESOURCES": str(outside)}, repo_root=repo)
        assert root == outside

    def test_explicit_argument_overrides_env_var(self, tmp_path: Path) -> None:
        explicit = tmp_path / "explicit"
        env_dir = tmp_path / "from_env"
        repo = tmp_path / "repo"
        repo.mkdir()
        root = resolve_resource_root(explicit, env={"MVA_RESOURCES": str(env_dir)}, repo_root=repo)
        assert root == explicit

    def test_root_inside_repo_is_refused(self, tmp_path: Path) -> None:
        repo = tmp_path / "repo"
        repo.mkdir()
        inside = repo / "resources"
        with pytest.raises(ResourceRootError, match="inside the repository"):
            resolve_resource_root(inside, env={}, repo_root=repo)

    def test_root_equal_to_repo_is_refused(self, tmp_path: Path) -> None:
        repo = tmp_path / "repo"
        repo.mkdir()
        with pytest.raises(ResourceRootError):
            resolve_resource_root(repo, env={}, repo_root=repo)


# --------------------------------------------------------------------------- is_download_stable


class TestIsDownloadStable:
    def test_unchanged_file_is_stable(self, tmp_path: Path) -> None:
        target = tmp_path / "steady.bin"
        target.write_bytes(b"x" * 100)
        assert is_download_stable(target, interval_seconds=0.01, sleep=lambda _s: None) is True

    def test_growing_file_is_not_stable(self, tmp_path: Path) -> None:
        target = tmp_path / "growing.bin"
        target.write_bytes(b"x" * 100)

        def grow(_seconds: float) -> None:
            target.write_bytes(target.read_bytes() + b"y" * 100)

        assert is_download_stable(target, interval_seconds=0.01, sleep=grow) is False


# --------------------------------------------------------------------------- sniff_content_mismatch


class TestSniffContentMismatch:
    def test_real_gzip_passes(self, tmp_path: Path) -> None:
        target = tmp_path / "real.vcf.gz"
        with gzip.open(target, "wb") as handle:
            handle.write(b"##fileformat=VCFv4.2\n")
        assert sniff_content_mismatch(target) is None

    def test_html_error_page_named_gz_is_caught(self, tmp_path: Path) -> None:
        """This is the exact failure this tool hit for real while being built: a
        Gene2Phenotype download URL returned an HTML app shell with HTTP 200 for a
        `.csv.gz` request. This check is what makes that non-silent."""
        target = tmp_path / "DDG2P.csv.gz"
        target.write_text("<!doctype html>\n<html><body>nope</body></html>", encoding="utf-8")
        reason = sniff_content_mismatch(target)
        assert reason is not None
        assert "HTML" in reason

    def test_plain_text_named_gz_without_magic_bytes_is_caught(self, tmp_path: Path) -> None:
        target = tmp_path / "not_really.tbi"
        target.write_bytes(b"this is not bgzf")
        reason = sniff_content_mismatch(target)
        assert reason is not None
        assert "magic bytes" in reason

    def test_plain_text_file_with_expected_extension_passes(self, tmp_path: Path) -> None:
        target = tmp_path / "phenotype.hpoa"
        target.write_text("#description: test\ndatabase_id\tdisease_name\n", encoding="utf-8")
        assert sniff_content_mismatch(target) is None


# --------------------------------------------------------------------------- survey_resource


class TestSurveyResource:
    def test_missing_file_is_not_fetched(self, tmp_path: Path) -> None:
        entry = _make_entry(path="does/not/exist.txt")
        surveyed = survey_resource(tmp_path, entry)
        assert surveyed.status is ResourceStatus.NOT_FETCHED
        assert surveyed.sha256 is None
        assert "no local file" in surveyed.notes

    def test_growing_file_is_not_fetched_and_not_hashed(self, tmp_path: Path) -> None:
        """The core guarantee this tool was built for: never hash a moving target."""
        (tmp_path / "sub").mkdir()
        target = tmp_path / "sub" / "growing.vcf.bgz"
        target.write_bytes(b"\x1f\x8b" + b"0" * 100)
        entry = _make_entry(path="sub/growing.vcf.bgz")

        def grow(_seconds: float) -> None:
            target.write_bytes(target.read_bytes() + b"0" * 100)

        surveyed = survey_resource(tmp_path, entry, stability_interval=0.01, sleep=grow)
        assert surveyed.status is ResourceStatus.NOT_FETCHED
        assert surveyed.sha256 is None
        assert surveyed.size_bytes is None
        assert "in progress" in surveyed.notes

    def test_content_mismatch_is_not_fetched(self, tmp_path: Path) -> None:
        (tmp_path / "sub").mkdir()
        target = tmp_path / "sub" / "DDG2P.csv.gz"
        target.write_text("<!doctype html>", encoding="utf-8")
        entry = _make_entry(path="sub/DDG2P.csv.gz")
        surveyed = survey_resource(tmp_path, entry, stability_interval=0.01, sleep=lambda _s: None)
        assert surveyed.status is ResourceStatus.NOT_FETCHED
        assert "content sanity check" in surveyed.notes

    def test_stable_valid_file_is_fetched_and_hashed(self, tmp_path: Path) -> None:
        (tmp_path / "sub").mkdir()
        target = tmp_path / "sub" / "clean.tsv"
        target.write_text("a\tb\n1\t2\n", encoding="utf-8")
        entry = _make_entry(path="sub/clean.tsv")
        surveyed = survey_resource(tmp_path, entry, stability_interval=0.01, sleep=lambda _s: None)
        assert surveyed.status is ResourceStatus.FETCHED
        assert surveyed.sha256 is not None
        assert len(surveyed.sha256) == 64
        assert surveyed.size_bytes == target.stat().st_size
        assert surveyed.retrieved is not None
        assert surveyed.notes == ""

    def test_survey_all_preserves_order_and_count(self, tmp_path: Path) -> None:
        entries = tuple(_make_entry(name=f"r{i}", path=f"r{i}.txt") for i in range(3))
        surveyed = survey_all(tmp_path, entries, stability_interval=0.01, sleep=lambda _s: None)
        assert [e.name for e in surveyed] == [e.name for e in entries]
        assert all(e.status is ResourceStatus.NOT_FETCHED for e in surveyed)


# --------------------------------------------------------------------------- fetch_resource


class TestFetchResource:
    def test_disallowed_host_never_reaches_the_network(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Defence in depth: even a hand-built entry that bypasses ResourceEntry's own
        construction-time check (via model_construct) must be refused before any
        download is attempted."""
        called = False

        def fail_if_called(*_args: object, **_kwargs: object) -> None:
            nonlocal called
            called = True

        monkeypatch.setattr(fetch_module, "_http_download", fail_if_called)
        bad = ResourceEntry.model_construct(
            name="bad",
            url="https://evil.example.com/x",
            version="1",
            path="x.txt",
            license="L",
            description="D",
            synthetic=False,
            status=ResourceStatus.NOT_FETCHED,
            sha256=None,
            size_bytes=None,
            retrieved=None,
            notes="",
        )
        with pytest.raises(DisallowedHostError):
            fetch_resource(bad, tmp_path, repo_root=tmp_path.parent)
        assert called is False

    def test_resource_root_inside_repo_is_refused(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        def unreachable(_url: str, _dest: Path, *, resume_from: int, timeout: float) -> None:
            del resume_from, timeout
            pytest.fail("the containment check must reject this before any download is attempted")

        monkeypatch.setattr(fetch_module, "_http_download", unreachable)
        repo = tmp_path / "repo"
        repo.mkdir()
        entry = _make_entry()
        with pytest.raises(ResourceRootError):
            fetch_resource(entry, repo / "resources", repo_root=repo)

    def test_successful_fetch_records_real_hash_and_size(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        payload = b"hello world, this is fake ClinVar-shaped data\n"

        def fake_download(_url: str, dest: Path, *, resume_from: int, timeout: float) -> None:
            del resume_from, timeout
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(payload)

        monkeypatch.setattr(fetch_module, "_http_download", fake_download)
        repo = tmp_path / "repo"
        repo.mkdir()
        root = tmp_path / "resources"
        entry = _make_entry(path="sub/file.txt")

        result = fetch_resource(entry, root, repo_root=repo)

        assert result.status is ResourceStatus.FETCHED
        assert result.size_bytes == len(payload)
        assert result.retrieved is not None
        assert (root / "sub" / "file.txt").read_bytes() == payload

        import hashlib

        assert result.sha256 == hashlib.sha256(payload).hexdigest()

    def test_fetch_that_returns_html_is_not_marked_fetched(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Reproduces the real Gene2Phenotype failure: the transport succeeds (no
        exception, HTTP 200) but the bytes are wrong. This must never be recorded as
        a successful, hash-pinned fetch."""

        def fake_download(_url: str, dest: Path, *, resume_from: int, timeout: float) -> None:
            del resume_from, timeout
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_text("<!doctype html><html>nope</html>", encoding="utf-8")

        monkeypatch.setattr(fetch_module, "_http_download", fake_download)
        repo = tmp_path / "repo"
        repo.mkdir()
        root = tmp_path / "resources"
        entry = _make_entry(path="sub/DDG2P.csv.gz")

        result = fetch_resource(entry, root, repo_root=repo)

        assert result.status is ResourceStatus.NOT_FETCHED
        assert result.sha256 is None
        assert "content sanity check" in result.notes

    def test_transport_failure_raises_resource_fetch_error_after_retries(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        attempts = 0

        def always_fails(_url: str, _dest: Path, *, resume_from: int, timeout: float) -> None:
            del resume_from, timeout
            nonlocal attempts
            attempts += 1
            raise OSError("connection reset")

        monkeypatch.setattr(fetch_module, "_http_download", always_fails)
        repo = tmp_path / "repo"
        repo.mkdir()
        entry = _make_entry()

        with pytest.raises(ResourceFetchError):
            fetch_resource(
                entry,
                tmp_path / "resources",
                repo_root=repo,
                max_attempts=2,
                retry_delay=lambda _s: None,
            )
        assert attempts == 2

    def test_resume_uses_existing_partial_size(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        repo = tmp_path / "repo"
        repo.mkdir()
        root = tmp_path / "resources"
        dest = root / "sub" / "file.txt"
        dest.parent.mkdir(parents=True)
        dest.write_bytes(b"partial-bytes")

        seen_resume_from: list[int] = []

        def fake_download(_url: str, d: Path, *, resume_from: int, timeout: float) -> None:
            del timeout
            seen_resume_from.append(resume_from)
            d.write_bytes(b"partial-bytes-completed")

        monkeypatch.setattr(fetch_module, "_http_download", fake_download)
        entry = _make_entry(path="sub/file.txt")
        fetch_resource(entry, root, repo_root=repo)
        assert seen_resume_from == [len(b"partial-bytes")]


# --------------------------------------------------------------------------- verify


class TestVerify:
    def test_not_fetched_entry_is_pending(self, tmp_path: Path) -> None:
        entry = _make_entry()
        result = verify_resource(tmp_path, entry)
        assert result.status is VerificationStatus.PENDING

    def test_fetched_entry_missing_from_disk_is_missing(self, tmp_path: Path) -> None:
        entry = _make_entry(
            status=ResourceStatus.FETCHED, sha256=_HEX64, size_bytes=1, retrieved="2026-08-27"
        )
        result = verify_resource(tmp_path, entry)
        assert result.status is VerificationStatus.MISSING

    def test_fetched_entry_with_wrong_hash_is_mismatch(self, tmp_path: Path) -> None:
        target = tmp_path / "somewhere" / "file.txt"
        target.parent.mkdir(parents=True)
        target.write_text("real content", encoding="utf-8")
        entry = _make_entry(
            status=ResourceStatus.FETCHED, sha256=_HEX64, size_bytes=1, retrieved="2026-08-27"
        )
        result = verify_resource(tmp_path, entry)
        assert result.status is VerificationStatus.MISMATCH
        # PRIV-09: the message names the file and carries hash digests, never content.
        assert "real content" not in result.message

    def test_fetched_entry_with_correct_hash_is_ok(self, tmp_path: Path) -> None:
        target = tmp_path / "somewhere" / "file.txt"
        target.parent.mkdir(parents=True)
        target.write_text("real content", encoding="utf-8")

        import hashlib

        digest = hashlib.sha256(b"real content").hexdigest()
        entry = _make_entry(
            status=ResourceStatus.FETCHED, sha256=digest, size_bytes=12, retrieved="2026-08-27"
        )
        result = verify_resource(tmp_path, entry)
        assert result.status is VerificationStatus.OK

    def test_assert_verified_raises_on_missing(self, tmp_path: Path) -> None:
        entry = _make_entry(
            status=ResourceStatus.FETCHED, sha256=_HEX64, size_bytes=1, retrieved="2026-08-27"
        )
        with pytest.raises(ResourceVerificationError):
            assert_verified(tmp_path, [entry])

    def test_assert_verified_passes_for_pending_and_ok(self, tmp_path: Path) -> None:
        target = tmp_path / "somewhere" / "file.txt"
        target.parent.mkdir(parents=True)
        target.write_text("real content", encoding="utf-8")

        import hashlib

        digest = hashlib.sha256(b"real content").hexdigest()
        ok_entry = _make_entry(
            name="ok",
            path="somewhere/file.txt",
            status=ResourceStatus.FETCHED,
            sha256=digest,
            size_bytes=12,
            retrieved="2026-08-27",
        )
        pending_entry = _make_entry(name="pending", path="elsewhere/other.txt")
        results = assert_verified(tmp_path, [ok_entry, pending_entry])
        assert {r.status for r in results} == {VerificationStatus.OK, VerificationStatus.PENDING}

    def test_verify_all_matches_length_and_order(self, tmp_path: Path) -> None:
        entries = [_make_entry(name=f"r{i}", path=f"r{i}.txt") for i in range(4)]
        results = verify_all(tmp_path, entries)
        assert [r.name for r in results] == [f"r{i}" for i in range(4)]


# --------------------------------------------------------------------------- manifest


class TestManifest:
    def test_render_is_deterministic(self) -> None:
        entries = (_make_entry(name="a"), _make_entry(name="b", path="other.txt"))
        first = render_resources_yaml(entries)
        second = render_resources_yaml(entries)
        assert first == second

    def test_render_sorts_by_name_regardless_of_input_order(self) -> None:
        forward = render_resources_yaml(
            (_make_entry(name="zzz"), _make_entry(name="aaa", path="a.txt"))
        )
        backward = render_resources_yaml(
            (_make_entry(name="aaa", path="a.txt"), _make_entry(name="zzz"))
        )
        assert forward == backward

    def test_render_rejects_duplicate_names(self) -> None:
        entries = (_make_entry(name="dup"), _make_entry(name="dup", path="other.txt"))
        with pytest.raises(AcquisitionError, match="duplicate"):
            render_resources_yaml(entries)

    def test_write_manifest_ends_with_trailing_newline(self, tmp_path: Path) -> None:
        manifest_path = tmp_path / "resources.yaml"
        write_resources_manifest(manifest_path, (_make_entry(),))
        text = manifest_path.read_text(encoding="utf-8")
        assert text.endswith("\n")
        assert not text.endswith("\n\n")

    def test_write_then_load_round_trips(self, tmp_path: Path) -> None:
        manifest_path = tmp_path / "resources.yaml"
        entries = (
            _make_entry(
                name="a",
                status=ResourceStatus.FETCHED,
                sha256=_HEX64,
                size_bytes=5,
                retrieved="2026-08-27",
            ),
            _make_entry(name="b", path="other.txt"),
        )
        write_resources_manifest(manifest_path, entries)
        loaded = load_resources_manifest(manifest_path)
        assert loaded.manifest_version == 1
        assert set(loaded.resources) == {"a", "b"}
        assert loaded.resources["a"].status is ResourceStatus.FETCHED
        assert loaded.resources["a"].sha256 == _HEX64
        assert loaded.resources["b"].status is ResourceStatus.NOT_FETCHED

    def test_load_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(AcquisitionError, match="not found"):
            load_resources_manifest(tmp_path / "does_not_exist.yaml")

    def test_load_invalid_yaml_raises(self, tmp_path: Path) -> None:
        manifest_path = tmp_path / "resources.yaml"
        manifest_path.write_text("resources: [this, is, not, a, mapping", encoding="utf-8")
        with pytest.raises(AcquisitionError):
            load_resources_manifest(manifest_path)

    def test_load_rejects_unknown_top_level_key(self, tmp_path: Path) -> None:
        manifest_path = tmp_path / "resources.yaml"
        write_resources_manifest(manifest_path, (_make_entry(),))
        mutated = manifest_path.read_text(encoding="utf-8") + "\nsome_unexpected_key: true\n"
        manifest_path.write_text(mutated, encoding="utf-8")
        with pytest.raises(AcquisitionError, match="does not match the expected schema"):
            load_resources_manifest(manifest_path)

    def test_committed_manifest_parses_and_matches_the_catalog(self) -> None:
        """The actual committed knowledge/manifests/resources.yaml, if present, must
        stay in sync with the declarative catalog: same set of resource names."""
        manifest_path = _REPO_ROOT / "knowledge" / "manifests" / "resources.yaml"
        if not manifest_path.is_file():
            pytest.skip("knowledge/manifests/resources.yaml has not been generated yet")
        manifest = load_resources_manifest(manifest_path)
        assert set(manifest.resources) == {entry.name for entry in KNOWN_RESOURCES}
        for entry in manifest.resources.values():
            assert entry.synthetic is False
            if entry.status is ResourceStatus.FETCHED:
                assert entry.sha256 is not None
                assert len(entry.sha256) == 64


# --------------------------------------------------------------------------- catalog


class TestCatalog:
    def test_all_names_are_unique(self) -> None:
        names = [entry.name for entry in KNOWN_RESOURCES]
        assert len(names) == len(set(names))

    def test_every_entry_is_declared_not_fetched_and_real(self) -> None:
        """The catalog is purely declarative: it never claims a hash itself. Every
        entry is also real public data (synthetic=False) -- this registry exists
        only for that; the synthetic demo tables live in knowledge.yaml instead."""
        for entry in KNOWN_RESOURCES:
            assert entry.status is ResourceStatus.NOT_FETCHED
            assert entry.sha256 is None
            assert entry.synthetic is False

    def test_every_url_is_allowlisted(self) -> None:
        # Redundant with ResourceEntry's own validator (which would have raised at
        # import time otherwise), but kept as an explicit, independent assertion.
        for entry in KNOWN_RESOURCES:
            assert_allowed_host(entry.url)

    def test_includes_the_full_gnomad_exome_chromosome_pattern(self) -> None:
        names = {entry.name for entry in KNOWN_RESOURCES}
        for chrom in (*[str(n) for n in range(1, 23)], "X", "Y"):
            assert f"gnomad_exomes_chr{chrom}" in names
            assert f"gnomad_exomes_chr{chrom}_tbi" in names


# --------------------------------------------------------------------------- the hard safety guard


class TestHardSafetyGuard:
    """PRIV-05: this tool must be structurally incapable of accepting, logging or
    transmitting anything that could be patient-derived. It takes no variant
    coordinates, no sample/proband/pedigree identifiers, and no path under a
    patient workspace -- checked here by walking every public callable's actual
    signature, so a future parameter addition that violates this fails a test
    instead of a code review."""

    _FORBIDDEN_PARAM_NAMES = frozenset(
        {
            "variant_id",
            "variant_ids",
            "sample_id",
            "sample_ids",
            "genotype",
            "genotypes",
            "pedigree",
            "proband_id",
            "patient_id",
            "workspace",
            "workspace_root",
        }
    )

    def _public_callables(self) -> list[Callable[..., object]]:
        import tools.acquire as pkg

        return [member for name in pkg.__all__ if callable(member := getattr(pkg, name))]

    def test_no_public_function_accepts_patient_shaped_parameters(self) -> None:
        violations: list[str] = []
        for obj in self._public_callables():
            try:
                signature = inspect.signature(obj)
            except (TypeError, ValueError):
                continue
            for param_name in signature.parameters:
                if param_name.lower() in self._FORBIDDEN_PARAM_NAMES:
                    violations.append(f"{obj!r} accepts parameter {param_name!r}")
        assert not violations, violations

    def test_resource_entry_has_no_patient_shaped_field(self) -> None:
        violations = [
            field
            for field in ResourceEntry.model_fields
            if field.lower() in self._FORBIDDEN_PARAM_NAMES
        ]
        assert not violations, violations

    def test_disallowed_host_is_refused_before_any_write(self, tmp_path: Path) -> None:
        """The end-to-end version of the allowlist guarantee: given an arbitrary,
        non-allowlisted URL, this tool raises and writes nothing -- it does not fetch."""
        bad = ResourceEntry.model_construct(
            name="bad",
            url="https://not-a-real-reference-host.example.com/data.vcf.gz",
            version="1",
            path="somewhere/data.vcf.gz",
            license="L",
            description="D",
            synthetic=False,
            status=ResourceStatus.NOT_FETCHED,
            sha256=None,
            size_bytes=None,
            retrieved=None,
            notes="",
        )
        with pytest.raises(DisallowedHostError):
            fetch_resource(bad, tmp_path, repo_root=tmp_path.parent)
        assert list(tmp_path.rglob("*")) == []
