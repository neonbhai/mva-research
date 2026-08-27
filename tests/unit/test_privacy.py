"""Privacy boundary enforcement tests (GP-40..GP-44).

Every payload used here is fabricated and lives only under ``tmp_path``. Where a
canary has to be recognisable to the scanner it is ASSEMBLED AT RUNTIME from inert
fragments — a bare coordinate, a bare digit run, a nucleotide unit — so that this
file does not itself become a scanner hit in ``content_scan``. That is not a
trick: the fragments are inert precisely because the rules are keyword- and
structure-anchored, which is the property the tests below verify.
"""

from __future__ import annotations

import logging
import shutil
import socket
import subprocess
from pathlib import Path

import pytest

from mva.config import find_repo_root, path_is_within, resolve_workspace
from mva.errors import (
    ExportBlockedError,
    NetworkDeniedError,
    PrivacyViolationError,
    WorkspaceError,
)
from mva.models.base import Sensitivity
from mva.privacy.audit import (
    CONTENT_ALLOWLIST_PREFIXES,
    PATH_DOWNGRADABLE_RULES,
    check_gitignore_effectiveness,
    check_log_redaction_probe,
    check_notebook_output_purity,
    check_symlink_escape,
    check_synthetic_fixtures_marked,
    check_workspace_containment,
    declares_synthetic,
    redact_path,
    run_audit,
    scan_bytes,
)
from mva.privacy.classify import classify_path, is_sensitive_extension
from mva.privacy.export import (
    PUBLIC_EXPORT_ALLOWLIST,
    export_public_artifact,
    gate_public_export,
)
from mva.privacy.netguard import BLOCKED_EVENTS, OfflineProfile, is_armed
from mva.privacy.patterns import RULES, decode_scrubbed, rule_by_id, sniff_binary
from mva.privacy.redact import install_redaction, redact_text

pytestmark = [pytest.mark.unit, pytest.mark.privacy]

# --- inert fragments -------------------------------------------------------
POS = "40200000"
CHROM = "chr15"
RECORD_ID = "4457821"
HPO_DIGITS = "0001250"
HPO_DIGITS_2 = "0000252"
HPO_DIGITS_3 = "0001511"
SEQ = "ACGT" * 12


def vcf_data_line() -> str:
    """A fabricated VCF record line, tab-separated at runtime."""
    return "\t".join((CHROM, POS, ".", "C", "T", "820.5", "PASS", "SYNTH", "GT:DP", "0/1:45"))


def coordinate_line() -> str:
    """A public-style coordinate row: VCF-shaped, but with no FORMAT/GT columns."""
    return "\t".join((CHROM, POS, ".", "C", "T", "820.5", "PASS", "AF=0.01"))


def vcf_document() -> str:
    return "\n".join(
        (
            "##fileformat=VCFv4.2",
            "##mva_synthetic=true",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH_PROBAND01",
            vcf_data_line(),
            "",
        )
    )


def hpo(digits: str) -> str:
    return f"HP:{digits}"


GIT = shutil.which("git") or "git"


def git(root: Path, *args: str, check: bool = False) -> subprocess.CompletedProcess[bytes]:
    """Run git against a throwaway repository under ``tmp_path``."""
    return subprocess.run(  # noqa: S603
        [GIT, "-C", str(root), *args], check=check, capture_output=True
    )


def matched_gitignore_pattern(root: Path, probe: str) -> tuple[int, str]:
    """``git check-ignore -v`` exit code and the pattern text it reported."""
    proc = git(root, "check-ignore", "-v", "--no-index", "--", probe)
    if proc.returncode != 0 or not proc.stdout:
        return proc.returncode, ""
    left = proc.stdout.decode().splitlines()[0].split("\t")[0]
    return proc.returncode, left.split(":", 2)[2]


def _init_repo(root: Path) -> None:
    """A throwaway git repo carrying this project's real .gitignore."""
    root.mkdir(parents=True, exist_ok=True)
    for args in (
        ("init", "-q"),
        ("config", "user.email", "t@example.invalid"),
        ("config", "user.name", "t"),
    ):
        git(root, *args, check=True)
    real = find_repo_root() / ".gitignore"
    (root / ".gitignore").write_bytes(real.read_bytes())


# ---------------------------------------------------------------------------
# 1. Sensitive artifacts cannot be exported (GP-43)
# ---------------------------------------------------------------------------


def test_sensitive_artifact_rejected_from_public_export(tmp_path: Path) -> None:
    artifact = tmp_path / "proband.vcf"
    artifact.write_text(vcf_document(), encoding="utf-8")

    decision = gate_public_export(artifact, declared=Sensitivity.PUBLIC, allowlist=["*.vcf"])
    assert not decision.allowed
    # Even with the allowlist and the declaration satisfied, the re-scan refuses.
    assert any(reason.startswith("content_scan") for reason in decision.reasons)
    assert "vcf_header" in decision.scanned_rules

    with pytest.raises(ExportBlockedError) as exc:
        export_public_artifact(
            artifact,
            tmp_path / "out" / "proband.vcf",
            declared=Sensitivity.PUBLIC,
            allowlist=["*.vcf"],
        )
    # The refusal names checks, never content (GP-41).
    assert POS not in str(exc.value)
    assert not (tmp_path / "out" / "proband.vcf").exists()


def test_declared_sensitive_is_refused_even_when_bytes_are_clean(tmp_path: Path) -> None:
    """Classification and scan are independent gates; failing either is enough."""
    artifact = tmp_path / "summary.csv"
    artifact.write_text("gene,count\nSYNTH1,3\n", encoding="utf-8")

    assert (
        gate_public_export(artifact, declared=Sensitivity.SENSITIVE, allowlist=["*.csv"]).allowed
        is False
    )
    assert (
        gate_public_export(artifact, declared=Sensitivity.PUBLIC, allowlist=["*.tsv"]).allowed
        is False
    )
    assert gate_public_export(artifact, declared=Sensitivity.PUBLIC, allowlist=["*.csv"]).allowed


def test_export_writes_only_when_allowed(tmp_path: Path) -> None:
    src = tmp_path / "submission.csv"
    src.write_text("rank,gene\n1,SYNTH1\n", encoding="utf-8")
    dest = tmp_path / "public" / "submission.csv"
    assert (
        export_public_artifact(src, dest, declared=Sensitivity.PUBLIC, allowlist=["submission.csv"])
        == dest
    )
    assert dest.read_text(encoding="utf-8") == src.read_text(encoding="utf-8")

    with pytest.raises(ExportBlockedError):
        export_public_artifact(
            src,
            tmp_path / "public" / "submission.bam",
            declared=Sensitivity.PUBLIC,
            allowlist=["submission.csv"],
        )


# ---------------------------------------------------------------------------
# 2 & 3. .gitignore actually ignores, and negations still work where intended
# ---------------------------------------------------------------------------


def test_sensitive_paths_are_git_ignored(tmp_path: Path) -> None:
    """Runs the real check, which requires a match that is NOT a negation rule."""
    repo = tmp_path / "repo"
    _init_repo(repo)
    result = check_gitignore_effectiveness(repo)
    assert result.passed, [(f.path, f.detail) for f in result.findings]
    assert not result.findings


@pytest.mark.parametrize(
    "probe",
    [
        "workspace/proband.vcf",
        "runs/case01/out.bam",
        "a/b/sample_R1.fastq.gz",
        "evidence.duckdb",
        "variants.parquet",
    ],
)
def test_individual_sensitive_probe_is_ignored_by_a_positive_rule(
    tmp_path: Path, probe: str
) -> None:
    repo = tmp_path / "repo"
    _init_repo(repo)
    code, pattern = matched_gitignore_pattern(repo, probe)
    assert code == 0, f"{probe} matched no rule at all"
    assert not pattern.startswith("!"), (
        f"{probe} matched NEGATION rule {pattern!r}; exit code 0 from check-ignore "
        "means 'a rule matched', which includes re-inclusions."
    )


def test_synthetic_fixtures_are_trackable(tmp_path: Path) -> None:
    """The narrow negations must still let synthetic fixtures be committed."""
    repo = tmp_path / "repo"
    _init_repo(repo)
    fixture = repo / "tests" / "fixtures" / "synthetic" / "case.vcf"
    fixture.parent.mkdir(parents=True)
    fixture.write_text(vcf_document(), encoding="utf-8")

    code, pattern = matched_gitignore_pattern(repo, "tests/fixtures/synthetic/case.vcf")
    assert code == 0
    assert pattern.startswith("!"), "the synthetic-fixture negation is no longer in force"

    status = git(repo, "status", "--porcelain", "--untracked-files=all", check=True)
    assert "tests/fixtures/synthetic/case.vcf" in status.stdout.decode()


# ---------------------------------------------------------------------------
# 4 & 5. Logs and text redaction (GP-42)
# ---------------------------------------------------------------------------


def test_logs_do_not_dump_sensitive_records(caplog: pytest.LogCaptureFixture) -> None:
    payload = vcf_data_line()
    buffer = logging.StreamHandler()
    records: list[str] = []

    class _Sink(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            records.append(self.format(record))

    sink = _Sink(level=logging.DEBUG)
    root = logging.getLogger()
    root.addHandler(sink)
    try:
        install_redaction()
        with caplog.at_level(logging.DEBUG):
            logging.getLogger("mva.ingestion.vcf").debug("row: %s", payload)
            logging.getLogger("mva.ingestion.vcf").debug(payload)
    finally:
        root.removeHandler(sink)
        buffer.close()

    captured = "\n".join(records) + "\n" + caplog.text
    assert POS not in captured
    assert f"{CHROM}\t{POS}" not in captured
    assert "REDACTED" in captured


def test_exception_tracebacks_are_stripped_from_records() -> None:
    """Frame locals in a traceback routinely hold the record that caused the error."""
    records: list[str] = []

    class _Sink(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            records.append(self.format(record))

    sink = _Sink(level=logging.DEBUG)
    root = logging.getLogger()
    root.addHandler(sink)
    try:
        install_redaction()
        try:
            raise ValueError(vcf_data_line())
        except ValueError:
            logging.getLogger("mva.ingestion.vcf").exception("parse failed")
    finally:
        root.removeHandler(sink)

    joined = "\n".join(records)
    assert POS not in joined
    assert "REDACTED:traceback" in joined


def test_redact_text_never_returns_the_original_substring() -> None:
    for payload, forbidden in (
        (vcf_data_line(), POS),
        (vcf_document(), "##fileformat=VCFv4.2"),
        ("\n".join(("@read_1", SEQ, "+", "I" * len(SEQ))), SEQ),
        (hpo(HPO_DIGITS), hpo(HPO_DIGITS)),
        (f"MRN: {RECORD_ID}", RECORD_ID),
        (f"patient_name: {'Alice'}", "Alice"),
        (f"DOB: {'2015-03-02'}", "2015-03-02"),
        ("@RG\tID:1\tSM:{}".format("HOSP-991"), "HOSP-991"),
    ):
        redacted = redact_text(payload)
        assert forbidden not in redacted, f"{forbidden!r} survived redaction"
        assert "<REDACTED:" in redacted


def test_redact_text_is_idempotent_and_leaves_ordinary_text_alone() -> None:
    plain = "ranked 4 candidate pairs for case SYNTH01"
    assert redact_text(plain) == plain
    once = redact_text(vcf_data_line())
    assert redact_text(once) == once


# ---------------------------------------------------------------------------
# 6. GP-41: findings never carry matched content
# ---------------------------------------------------------------------------


def test_finding_detail_never_contains_matched_content(tmp_path: Path) -> None:
    target = tmp_path / "leaky.txt"
    target.write_text(vcf_document(), encoding="utf-8")

    findings = scan_bytes(target.read_bytes(), check="content_scan", path_label="leaky.txt")
    assert findings
    assert any(f.severity == "fail" for f in findings)
    for finding in findings:
        assert POS not in finding.detail
        assert CHROM not in finding.detail
        assert finding.span_len is None or finding.span_len > 0

    repo = tmp_path / "repo"
    _init_repo(repo)
    (repo / "leaky.txt").write_text(vcf_document(), encoding="utf-8")
    report = run_audit(repo)
    markdown = report.to_markdown()
    assert "leaky.txt" in markdown
    assert POS not in markdown
    assert f"{CHROM}\t{POS}" not in markdown
    assert "##fileformat" not in markdown
    assert "content_scan" in report.failed_checks

    serialised = repr(report.to_dict())
    assert POS not in serialised


# ---------------------------------------------------------------------------
# 7. Offline profile (netguard)
# ---------------------------------------------------------------------------


def test_offline_profile_blocks_dns_and_only_while_armed() -> None:
    assert not is_armed()
    # Disarmed: a loopback lookup must still work, or the guard is unusable.
    socket.getaddrinfo("127.0.0.1", 80, proto=socket.IPPROTO_TCP)

    with pytest.raises(NetworkDeniedError) as exc, OfflineProfile():
        assert is_armed()
        socket.getaddrinfo("example.invalid", 443)
    assert "socket.getaddrinfo" in str(exc.value)
    # Event args (hostnames, ports) are never included (GP-41).
    assert "example.invalid" not in str(exc.value)

    assert not is_armed()
    socket.getaddrinfo("127.0.0.1", 80, proto=socket.IPPROTO_TCP)


def test_offline_profile_restores_nested_state() -> None:
    with OfflineProfile():
        with OfflineProfile(strict=True):
            assert is_armed()
        assert is_armed(), "inner profile must not disarm the outer one"
    assert not is_armed()


# ---------------------------------------------------------------------------
# 8 & 9. False-positive discipline
# ---------------------------------------------------------------------------


def test_mrn_rule_does_not_match_a_bare_genomic_position() -> None:
    rule = rule_by_id("mrn")
    for benign in (
        f"{CHROM}\t{POS}".encode(),
        b"POS\t40200000\t.\tC\tT",
        b"allele_count 1234567",
        b"the patient id field is populated downstream",
        b"offset 987654321 bytes",
    ):
        assert rule.pattern.search(benign) is None, f"mrn false-positived on {benign[:20]!r}"
    # It must still fire when the keyword anchor is present.
    assert rule.pattern.search(f"MRN: {RECORD_ID}".encode()) is not None


def test_hpo_rule_does_not_fail_on_a_single_term(tmp_path: Path) -> None:
    single = f"phenotype of interest: {hpo(HPO_DIGITS)}\n".encode()
    findings = scan_bytes(single, check="content_scan", path_label="notes.txt")
    hpo_findings = [f for f in findings if f.rule_id == "hpo_term"]
    assert hpo_findings
    assert all(f.severity == "warn" for f in hpo_findings)

    three = "\n".join(hpo(d) for d in (HPO_DIGITS, HPO_DIGITS_2, HPO_DIGITS_3)).encode()
    escalated = [
        f
        for f in scan_bytes(three, check="content_scan", path_label="profile.txt")
        if f.rule_id == "hpo_term"
    ]
    assert escalated and all(f.severity == "fail" for f in escalated)


def test_vcf_data_line_is_warn_alone_and_fail_beside_a_header() -> None:
    coordinates = "\n".join(coordinate_line() for _ in range(3)).encode()
    alone = [
        f
        for f in scan_bytes(coordinates, check="content_scan", path_label="public_table.tsv")
        if f.rule_id == "vcf_data_line"
    ]
    assert alone and all(f.severity == "warn" for f in alone)

    with_header = [
        f
        for f in scan_bytes(vcf_document().encode(), check="content_scan", path_label="x.vcf")
        if f.rule_id == "vcf_data_line"
    ]
    assert with_header and all(f.severity == "fail" for f in with_header)


def test_iso_date_is_warn_but_keyed_dob_fails() -> None:
    bare = scan_bytes(b"created_at: 2026-01-01T00:00:00Z", check="c", path_label="m.json")
    assert bare and all(f.severity == "warn" for f in bare)
    keyed = scan_bytes(f"DOB: {'2015-03-02'}".encode(), check="c", path_label="m.json")
    assert any(f.rule_id == "dob" and f.severity == "fail" for f in keyed)


def test_magic_byte_sniffers_are_not_regexes() -> None:
    assert sniff_binary(b"CRAM\x03\x00rest") == "cram"
    assert sniff_binary(b"not a container") is None
    bgzf_head = b"\x1f\x8b\x08\x04" + b"\x00" * 8 + b"BC\x02\x00" + b"\x1b\x00"
    assert sniff_binary(bgzf_head) == "bgzf"


def test_every_rule_declares_its_false_positive_risk() -> None:
    for rule in RULES:
        assert rule.false_positive_risk.strip(), rule.rule_id
        assert rule.severity in {"fail", "warn"}
    with pytest.raises(KeyError):
        rule_by_id("no_such_rule")


# ---------------------------------------------------------------------------
# 10 & 11. Workspace boundary (GP-40)
# ---------------------------------------------------------------------------


def test_workspace_inside_repo_is_rejected(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    (repo / "workspace").mkdir(parents=True)
    (repo / "pyproject.toml").write_text("[project]\nname='x'\n", encoding="utf-8")

    with pytest.raises(WorkspaceError) as exc:
        resolve_workspace(repo / "workspace", repo_root=repo)
    assert "inside the repository" in str(exc.value)

    # A symlink into the repo must not launder the containment check.
    outside = tmp_path / "link-to-inside"
    outside.symlink_to(repo / "workspace", target_is_directory=True)
    with pytest.raises(WorkspaceError):
        resolve_workspace(outside, repo_root=repo)

    external = tmp_path / "real-workspace"
    external.mkdir()
    assert resolve_workspace(external, repo_root=repo).root == external.resolve()


def test_symlink_into_sensitive_path_is_detected(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    _init_repo(repo)
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    patient = workspace / "proband.vcf"
    patient.write_text(vcf_document(), encoding="utf-8")

    (repo / "shortcut.vcf").symlink_to(patient)
    (repo / "ws").symlink_to(workspace, target_is_directory=True)

    result = check_symlink_escape(repo, workspace=workspace)
    assert not result.passed
    offenders = {f.path for f in result.findings if f.severity == "fail"}
    assert {"shortcut.vcf", "ws"} <= offenders
    for finding in result.findings:
        # The absolute workspace location is itself disclosive.
        assert str(workspace) not in finding.detail
        assert POS not in finding.detail

    clean = check_symlink_escape(repo, workspace=None)
    assert any(f.rule_id is None and f.severity == "fail" for f in clean.findings)


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("path", "sensitive"),
    [
        ("proband.vcf", True),
        ("proband.vcf.gz", True),
        ("s_R1.fastq.gz", True),
        ("aln.cram", True),
        ("evidence.duckdb", True),
        ("variants.parquet", True),
        ("family.ped", True),
        ("notes.md", False),
        ("knowledge/public/gene_panel.tsv", False),
    ],
)
def test_extension_classification(path: str, sensitive: bool) -> None:
    assert is_sensitive_extension(Path(path)) is sensitive


def test_classify_path_fails_closed() -> None:
    assert classify_path(Path("knowledge/public/gene_panel.tsv")) is Sensitivity.PUBLIC
    assert classify_path(Path("docs/architecture.md")) is Sensitivity.PUBLIC
    assert classify_path(Path("runs/case01/summary.tsv")) is Sensitivity.SENSITIVE
    assert classify_path(Path("workspace/notes.txt")) is Sensitivity.SENSITIVE
    # Unrecognised is SENSITIVE, never "probably fine".
    assert classify_path(Path("scratch/whatever.bin")) is Sensitivity.SENSITIVE


# ---------------------------------------------------------------------------
# End-to-end audit shape
# ---------------------------------------------------------------------------


def test_audit_detects_a_tracked_sensitive_file(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    _init_repo(repo)
    leaked = repo / "proband.vcf"
    leaked.write_text(vcf_document(), encoding="utf-8")
    git(repo, "add", "-f", "proband.vcf", check=True)

    staged = run_audit(repo, staged_only=True)
    assert "git_staged_sensitive" in staged.failed_checks

    git(repo, "commit", "-qm", "leak", check=True)
    full = run_audit(repo)
    assert "git_tracked_sensitive" in full.failed_checks
    assert POS not in full.to_markdown()


def test_audit_report_is_serialisable_and_content_free(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    _init_repo(repo)
    report = run_audit(repo)
    payload = report.to_dict()
    assert set(payload) == {"passed", "failed_checks", "results"}
    assert isinstance(payload["results"], list)
    assert report.to_markdown().startswith("# Privacy audit")


def test_log_redaction_probe_passes_on_the_configured_stack() -> None:
    # Armed explicitly, the way a composition root must: the probe no longer
    # installs redaction before measuring it, so this test would otherwise be
    # asserting on whichever earlier test happened to arm the process.
    install_redaction()
    report = run_audit(find_repo_root(), staged_only=True)
    probe = next(r for r in report.results if r.name == "log_redaction_probe")
    assert probe.passed, [f.detail for f in probe.findings]


# ---------------------------------------------------------------------------
# 12. The content allowlist is not a bypass (V1)
#
# Every payload below is fabricated. As everywhere in this file, the recognisable
# strings are ASSEMBLED AT RUNTIME so that this module does not itself become a
# fail-severity hit in the audit's own content_scan.
# ---------------------------------------------------------------------------


def real_shaped_vcf(sample: str = "PROBAND") -> str:
    """A structurally complete VCF: header, column row, and a called genotype.

    This is what a real call set looks like to the scanner, and what a reviewer
    committed under every allowlisted prefix to prove the allowlist was an
    unconditional bypass. It carries NO synthetic declaration, which is the whole
    point: an undeclared file must not inherit a directory's reputation.
    """
    return "\n".join(
        (
            "##fileformat=VCFv4.2",
            "##reference=GRCh38",
            "\t".join(
                (
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "INFO",
                    "FORMAT",
                    sample,
                )
            ),
            "\t".join((CHROM, POS, ".", "C", "T", "50", "PASS", ".", "GT:DP", "0/1:30")),
            "",
        )
    )


@pytest.mark.parametrize("prefix", CONTENT_ALLOWLIST_PREFIXES)
def test_a_real_shaped_vcf_fails_under_every_allowlisted_prefix(
    tmp_path: Path, prefix: str
) -> None:
    """The allowlist may soften documented false positives, never a whole file.

    Before this, `_resolve_context` downgraded EVERY fail to warn for any path
    under an allowlisted prefix, so a complete VCF committed to
    `knowledge/public/`, `knowledge/manifests/` or `tests/golden/` produced a
    dozen rule hits and `passed=True`.
    """
    repo = tmp_path / "repo"
    _init_repo(repo)
    for name in ("variants.tsv", "expected.csv", "patient.vcf"):
        target = repo / prefix / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(real_shaped_vcf(), encoding="utf-8")

    report = run_audit(repo)
    assert report.passed is False
    assert "content_scan" in report.failed_checks
    failed_rules = {f.rule_id for f in report.findings if f.severity == "fail"}
    assert {"vcf_header", "vcf_chrom_line", "genotype_field"} <= failed_rules


def test_allowlist_downgrades_only_the_documented_false_positive_rules() -> None:
    """`vcf_data_line` and `hpo_term` soften on the path; nothing else does."""
    assert {"vcf_data_line", "hpo_term"} == PATH_DOWNGRADABLE_RULES
    findings = scan_bytes(
        real_shaped_vcf().encode(),
        check="content_scan",
        path_label="knowledge/public/variants.tsv",
        allowlisted=True,
    )
    by_rule = {f.rule_id: f.severity for f in findings}
    assert by_rule["vcf_header"] == "fail"
    assert by_rule["vcf_chrom_line"] == "fail"
    assert by_rule["genotype_field"] == "fail"
    # The two documented FP-prone rules are still softened.
    assert by_rule["vcf_data_line"] == "warn"


def test_a_declared_synthetic_fixture_is_still_downgraded() -> None:
    """The escape hatch that keeps the repo's own fixture legal, and only that.

    A declaration is believed only when the marker is in the HEAD of the file and
    every named subject carries a synthetic prefix.
    """
    declared = "\n".join(
        (
            "##fileformat=VCFv4.2",
            "##mva_synthetic=true",
            "\t".join(
                (
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "INFO",
                    "FORMAT",
                    "SYNTH_PROBAND01",
                )
            ),
            "\t".join((CHROM, POS, ".", "C", "T", "50", "PASS", ".", "GT:DP", "0/1:30")),
            "",
        )
    )
    assert declares_synthetic(declared.encode())
    softened = scan_bytes(
        declared.encode(),
        check="content_scan",
        path_label="tests/fixtures/synthetic/case.vcf",
        allowlisted=True,
    )
    assert all(f.severity == "warn" for f in softened)

    # Same marker, but the sample column names an undeclared subject.
    swapped = declared.replace("SYNTH_PROBAND01", "PROBAND01")
    assert not declares_synthetic(swapped.encode())
    assert any(
        f.severity == "fail"
        for f in scan_bytes(
            swapped.encode(),
            check="content_scan",
            path_label="tests/fixtures/synthetic/case.vcf",
            allowlisted=True,
        )
    )


def test_container_magic_is_never_softened_by_the_allowlist() -> None:
    """A BAM inside the synthetic fixture directory is a BAM."""
    findings = scan_bytes(
        b"CRAM\x03\x00" + b"\x00" * 64,
        check="content_scan",
        path_label="tests/fixtures/synthetic/aln.cram",
        allowlisted=True,
    )
    assert any(f.rule_id == "magic:cram" and f.severity == "fail" for f in findings)


# ---------------------------------------------------------------------------
# 13. Case-insensitive paths do not launder workspace containment (V2, GP-40)
# ---------------------------------------------------------------------------


def _case_variant_is_the_same_dir(repo: Path) -> bool:
    """True on a case-insensitive filesystem (APFS, NTFS), where this bug lives."""
    upper = repo.parent / repo.name.upper()
    try:
        return upper.is_dir() and upper.stat().st_ino == repo.stat().st_ino
    except OSError:
        return False


def test_case_variant_spelling_does_not_escape_the_repo(tmp_path: Path) -> None:
    """`Path.resolve()` does not case-fold, and on APFS that is a containment hole.

    `.../ci/REPO/WS_INSIDE` and `.../ci/repo/ws_inside` are ONE directory on the
    default macOS filesystem. `is_relative_to` compared strings and reported them
    as unrelated, so a workspace physically inside the repository was accepted.
    """
    repo = tmp_path / "repo"
    inside = repo / "ws_inside"
    inside.mkdir(parents=True)
    (repo / "pyproject.toml").write_text("[project]\nname='x'\n", encoding="utf-8")

    if not _case_variant_is_the_same_dir(repo):
        pytest.skip("case-sensitive filesystem; the case-variant hole cannot occur here")

    variant = tmp_path / "REPO" / "WS_INSIDE"
    assert not variant.resolve().is_relative_to(repo.resolve()), (
        "precondition: this is exactly the comparison that used to be trusted"
    )
    assert path_is_within(variant, repo)

    with pytest.raises(WorkspaceError, match="inside the repository"):
        resolve_workspace(variant, repo_root=repo)

    result = check_workspace_containment(repo, workspace=variant)
    assert not result.passed
    assert any(f.severity == "fail" for f in result.findings)
    # The location of the patient directory is itself disclosive.
    assert all(str(variant) not in f.detail for f in result.findings)


def test_containment_survives_a_symlinked_intermediate_component(tmp_path: Path) -> None:
    """Identity comparison, not string comparison: a linked parent is still the parent."""
    repo = tmp_path / "repo"
    (repo / "sub").mkdir(parents=True)
    link = tmp_path / "repo-link"
    link.symlink_to(repo, target_is_directory=True)
    assert path_is_within(link / "sub", repo)
    assert not path_is_within(tmp_path / "elsewhere", repo)


# ---------------------------------------------------------------------------
# 14. The report never echoes a filename that IS the identifier (V4, GP-41)
# ---------------------------------------------------------------------------

NHS_DIGITS = "9999999999"
FAKE_SURNAME = "Faketon"


def phi_filename() -> str:
    """A fabricated sequencing filename of the shape real ones have."""
    return f"NHS{NHS_DIGITS}_{FAKE_SURNAME}.vcf.txt"


def test_a_filename_carrying_an_identifier_is_redacted_in_the_report(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    _init_repo(repo)
    leaky = repo / "data" / phi_filename()
    leaky.parent.mkdir(parents=True)
    leaky.write_text(real_shaped_vcf(), encoding="utf-8")

    report = run_audit(repo)
    assert report.passed is False

    markdown = report.to_markdown()
    serialised = repr(report.to_dict())
    for rendered in (markdown, serialised):
        assert NHS_DIGITS not in rendered
        assert FAKE_SURNAME not in rendered
        assert phi_filename() not in rendered
    # The directory chain and the extension survive: enough to act on, not enough
    # to reconstruct.
    assert "data/<REDACTED:path:" in markdown
    assert ".vcf.txt" in markdown


def test_the_path_itself_is_scanned_not_only_the_bytes(tmp_path: Path) -> None:
    """A file whose NAME is the PHI carries no genomic content at all."""
    repo = tmp_path / "repo"
    _init_repo(repo)
    target = repo / "data" / f"NHS{NHS_DIGITS}_{FAKE_SURNAME}.txt"
    target.parent.mkdir(parents=True)
    target.write_text("aggregate counts only\n", encoding="utf-8")

    report = run_audit(repo)
    assert "content_scan" in report.failed_checks
    assert any(f.rule_id == "path_identifier" and f.severity == "fail" for f in report.findings)
    assert NHS_DIGITS not in report.to_markdown()


def test_redact_path_leaves_ordinary_paths_alone() -> None:
    for benign in (
        "src/mva/privacy/audit.py",
        "knowledge/public/gene_panel.tsv",
        "tests/golden/expected_ranking.tsv",
        ".",
        "<logging>",
        "$MVA_WORKSPACE",
    ):
        assert redact_path(benign) == benign


# ---------------------------------------------------------------------------
# 15. The Track 1 export gate acts on its verdict (V5, GP-43)
# ---------------------------------------------------------------------------


def test_track1_submission_is_deleted_when_the_gate_refuses(tmp_path: Path) -> None:
    """`gate_public_export` RETURNS a verdict and never raises for a refusal.

    The caller wrapped it in `try/except` and ignored `decision.allowed`, so a
    verdict of `allowed=False` deleted nothing and raised nothing, while the
    docstring claimed "the file is deleted and the error propagates: fail closed".
    """
    from mva.reporting.track1 import ACCEPTED_PROBAND_ID, SubmissionRow, write_submission

    row = SubmissionRow(
        proband_id=ACCEPTED_PROBAND_ID,
        chrom_1="chr1",
        pos_1=1_000,
        ref_1="A",
        alt_1="T",
        chrom_2="",
        pos_2="",
        ref_2="",
        alt_2="",
        epcr=0.5,
    )
    # A name that is NOT on the deny-by-default allowlist.
    target = tmp_path / "not_on_the_allowlist.csv"
    with pytest.raises(ExportBlockedError, match="allowlist"):
        write_submission([row], target)
    assert not target.exists(), "a refused submission must not survive on disk"

    allowed = tmp_path / "track1_submission.csv"
    assert write_submission([row], allowed) == allowed
    assert allowed.exists()


def test_the_export_allowlist_is_not_derived_from_the_file(tmp_path: Path) -> None:
    """`allowlist=(path.name,)` is self-satisfying: every file is on it."""
    artifact = tmp_path / "anything_at_all.csv"
    artifact.write_text("gene,count\nSYNTH1,3\n", encoding="utf-8")
    assert gate_public_export(
        artifact, declared=Sensitivity.PUBLIC, allowlist=(artifact.name,)
    ).allowed, "precondition: a self-derived allowlist cannot refuse anything"
    decision = gate_public_export(
        artifact, declared=Sensitivity.PUBLIC, allowlist=PUBLIC_EXPORT_ALLOWLIST
    )
    assert not decision.allowed
    assert any(reason.startswith("allowlist") for reason in decision.reasons)


def test_the_shared_allowlist_has_not_drifted_from_the_composition_root() -> None:
    """One policy, one list. `mva.privacy.export` owns it; the orchestrator reuses it."""
    from mva.orchestrator import PUBLIC_EXPORT_ALLOWLIST as ORCHESTRATOR_ALLOWLIST

    assert set(ORCHESTRATOR_ALLOWLIST) == set(PUBLIC_EXPORT_ALLOWLIST)


# ---------------------------------------------------------------------------
# 16. Regex false negatives (V6)
#
# Each case below was verified as producing NO fail finding before the widening.
# ---------------------------------------------------------------------------

REFSEQ_CONTIG = "NC_000012.12"


def fail_rules(payload: bytes, *, label: str = "leak.txt") -> set[str | None]:
    return {
        f.rule_id
        for f in scan_bytes(payload, check="content_scan", path_label=label)
        if f.severity == "fail"
    }


def test_refseq_contig_names_are_recognised_as_vcf_records() -> None:
    """CHROM is a contig NAME, not a chromosome number. `NC_000012.12` is a contig."""
    body = "\n".join(
        (
            "\t".join((REFSEQ_CONTIG, POS, ".", "C", "T", "50", "PASS", ".", "GT:DP", "0/1:30")),
            "",
        )
    ).encode()
    assert rule_by_id("vcf_data_line").pattern.search(body) is not None
    assert "vcf_data_line" in fail_rules(body)


def test_mutect2_format_keys_are_recognised() -> None:
    """`[A-Z]{1,4}` rejected `F1R2`, so a somatic call set matched nothing."""
    body = "\n".join(
        (
            "\t".join(
                (CHROM, POS, ".", "C", "T", "50", "PASS", ".", "GT:AD:AF:F1R2:F2R1", "0/1:5,5")
            ),
            "",
        )
    ).encode()
    assert "genotype_field" in fail_rules(body)


def test_a_haploid_genotype_is_recognised() -> None:
    """chrY, chrM and a hemizygous X are written as a single allele, no separator."""
    body = "\n".join(
        ("\t".join(("chrY", POS, ".", "C", "T", "50", "PASS", ".", "GT:DP", "1:30")), "")
    ).encode()
    assert "genotype_field" in fail_rules(body)


def test_a_plain_gzipped_vcf_is_decompressed_and_scanned() -> None:
    """Compression is not a privacy control. Plain gzip matched nothing at all."""
    import gzip

    payload = gzip.compress(real_shaped_vcf().encode())
    assert sniff_binary(payload[:64]) is None, "precondition: not BGZF, so magic sniffing is blind"
    rules = fail_rules(payload, label="calls.vcf.gz")
    assert {"vcf_header", "vcf_chrom_line", "genotype_field"} <= rules


def test_a_fasta_record_is_recognised() -> None:
    """No rule existed for FASTA at all."""
    body = "\n".join((">" + "synthetic_consensus", SEQ * 2, "")).encode()
    assert "fasta_record" in fail_rules(body, label="consensus.txt")


def test_a_comma_separated_phi_table_is_recognised() -> None:
    """The keyed rules required `[:=#]`, so a CSV — the likeliest carrier — was blind."""
    body = "\n".join(
        (
            f"mrn,{RECORD_ID}",
            f"dob,{'2015-03-02'}",
            f"patient_surname,{FAKE_SURNAME}",
            "",
        )
    ).encode()
    assert {"mrn", "dob", "person_name_keyed"} <= fail_rules(body, label="cohort.csv")

    tabbed = "\n".join((f"MRN\t{RECORD_ID}", "")).encode()
    assert "mrn" in fail_rules(tabbed, label="cohort.tsv")


def test_the_obo_underscore_form_of_an_hpo_term_is_recognised() -> None:
    """`HP_0000001` is how ontology dumps, RDF exports and file names write it."""
    profile = "\n".join(f"HP_{digits}" for digits in (HPO_DIGITS, HPO_DIGITS_2, HPO_DIGITS_3))
    findings = scan_bytes(profile.encode(), check="content_scan", path_label="profile.txt")
    hpo_findings = [f for f in findings if f.rule_id == "hpo_term"]
    assert hpo_findings and all(f.severity == "fail" for f in hpo_findings)


def test_a_plink_ped_body_is_recognised() -> None:
    """A .ped line is whitespace-delimited, so no VCF rule ever saw it."""
    body = "\n".join(
        (" ".join(("FAM01", "SYNTH_IND01", "0", "0", "1", "2", *"AGCTAATT")), "")
    ).encode()
    assert "plink_ped_line" in fail_rules(body, label="family.txt")


def test_a_space_delimited_read_group_is_recognised() -> None:
    """SAM headers survive round-trips through tools that collapse tabs to spaces."""
    body = "\n".join((" ".join(("@RG", "ID:1", "SM:HOSP-991")), "")).encode()
    assert "sam_rg_sample" in fail_rules(body, label="header.txt")


# ---------------------------------------------------------------------------
# 17. Log redaction covers extra= and handlers added later (V7, GP-42)
# ---------------------------------------------------------------------------


class _Recorder(logging.Handler):
    """Captures both the formatted line and the record's own attribute dict."""

    def __init__(self) -> None:
        super().__init__(level=logging.DEBUG)
        self.lines: list[str] = []
        self.attributes: list[dict[str, object]] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.lines.append(self.format(record))
        self.attributes.append(dict(record.__dict__))


def test_structured_extra_fields_are_scrubbed() -> None:
    """`Logger.makeRecord` writes `extra` into `record.__dict__` AFTER the factory.

    The factory therefore never saw it, `scrub_record` never walked `__dict__`,
    and a JSON/OTLP handler serialised exactly those keys verbatim.
    """
    install_redaction()
    sink = _Recorder()
    root = logging.getLogger()
    root.addHandler(sink)
    try:
        logging.getLogger("mva.ingestion.vcf").warning(
            "loaded", extra={"extra_sample": f"MRN: {RECORD_ID} {hpo(HPO_DIGITS)}"}
        )
    finally:
        root.removeHandler(sink)

    assert sink.attributes
    survived = str(sink.attributes[-1]["extra_sample"])
    assert RECORD_ID not in survived
    assert hpo(HPO_DIGITS) not in survived
    assert "REDACTED" in survived


def test_a_handler_added_after_installation_is_still_armed() -> None:
    """Libraries add handlers on import, on first use and in `basicConfig`."""
    install_redaction()
    sink = _Recorder()  # constructed AFTER install_redaction()
    record = logging.makeLogRecord(
        {
            "msg": f"MRN: {RECORD_ID}",
            "levelno": logging.WARNING,
            "levelname": "WARNING",
            "name": "third.party",
        }
    )
    sink.handle(record)
    assert RECORD_ID not in "\n".join(sink.lines)
    assert "REDACTED" in "\n".join(sink.lines)


# ---------------------------------------------------------------------------
# 18. redact_text knows this project's own canonical formats (V8)
# ---------------------------------------------------------------------------


def test_redact_text_covers_the_canonical_variant_id() -> None:
    """The exact string other modules interpolate into messages."""
    variant_id = ":".join(("GRCh38", CHROM, POS, "C", "T"))
    redacted = redact_text(f"no evidence for {variant_id}")
    assert POS not in redacted
    assert variant_id not in redacted
    assert "REDACTED" in redacted


def test_redact_text_covers_a_bare_genotype_token() -> None:
    for token in ("0/1", "1|0", "./."):
        redacted = redact_text(f"genotype {token} phased")
        assert token not in redacted
        assert "REDACTED" in redacted
    # The lookarounds keep it to one allele a side.
    assert redact_text("a 50/50 split on 1/2/2020") == "a 50/50 split on 1/2/2020"


def test_the_greedy_redaction_rules_are_not_audit_rules() -> None:
    """A `0/1` in a code comment must not turn the audit report into noise."""
    assert {"variant_id", "genotype_token"}.isdisjoint({rule.rule_id for rule in RULES})


# ---------------------------------------------------------------------------
# 19. The medium-severity controls that were theatre
# ---------------------------------------------------------------------------


def test_the_log_probe_reports_a_stack_that_was_never_armed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The probe used to call `install_redaction()` before measuring anything.

    Measuring a stack it had just repaired, it could only ever report success —
    so the one failure it exists to report, "the application forgot to arm
    GP-42", was structurally undetectable.
    """
    monkeypatch.setattr("mva.privacy.audit.redaction_installed", lambda: False)
    result = check_log_redaction_probe()
    assert not result.passed
    assert any("NOT armed" in f.detail and f.severity == "fail" for f in result.findings)


def test_the_log_probe_walks_the_whole_logger_registry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A handler on a library's own logger is never consulted by the root."""
    monkeypatch.setattr("mva.privacy.audit.unfiltered_handlers", lambda: ["ThirdPartyHandler"])
    result = check_log_redaction_probe()
    assert any("ThirdPartyHandler" in f.detail for f in result.findings)


@pytest.mark.parametrize("event", ["socket.sendto", "socket.bind", "socket.gethostbyaddr"])
def test_netguard_blocks_the_pure_python_datagram_routes(event: str) -> None:
    assert event in BLOCKED_EVENTS


def test_a_udp_datagram_cannot_leave_an_armed_profile() -> None:
    """`sendto` never emits `socket.connect`: a few lines of pure Python sufficed."""
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    try:
        with pytest.raises(NetworkDeniedError) as exc, OfflineProfile():
            sock.sendto(hpo(HPO_DIGITS).encode(), ("192.0.2.1", 53))
        assert "socket.sendto" in str(exc.value)
        # Event args (the destination) are never included (GP-41).
        assert "192.0.2.1" not in str(exc.value)
    finally:
        sock.close()


def test_decode_scrubbed_has_a_real_call_site(tmp_path: Path) -> None:
    """A control with no call site protects nothing.

    `UnicodeDecodeError` embeds the offending bytes in `str(exc)`, and a
    traceback travels to the terminal, the CI log and a model context. The
    notebook check now decodes through the guard instead of handing raw bytes to
    `json.loads`.
    """
    undecodable = b'{"cells": [' + bytes([0xFF, 0xFE, 0xFD]) + b"]}"
    with pytest.raises(PrivacyViolationError) as exc:
        decode_scrubbed(undecodable, path=Path("notes.ipynb"))
    assert "Content withheld" in str(exc.value)

    repo = tmp_path / "repo"
    _init_repo(repo)
    notebook = repo / "notes.ipynb"
    notebook.write_bytes(undecodable)
    git(repo, "add", "-f", "notes.ipynb", check=True)

    result = check_notebook_output_purity(repo)
    assert result.findings
    for finding in result.findings:
        assert "\\xff" not in finding.detail
        assert "0xff" not in finding.detail.lower()


def test_a_synthetic_marker_must_be_in_the_head_of_the_file(tmp_path: Path) -> None:
    """`SYNTH_` anywhere in 8 MiB was a declaration. A declaration is a header."""
    repo = tmp_path / "repo"
    _init_repo(repo)
    fixture = repo / "tests" / "fixtures" / "synthetic" / "buried.vcf"
    fixture.parent.mkdir(parents=True)
    fixture.write_text("# padding\n" * 2000 + "# SYNTH_marker\n", encoding="utf-8")

    result = check_synthetic_fixtures_marked(repo)
    assert not result.passed
    assert any("first" in f.detail for f in result.findings)


def test_a_declared_fixture_must_name_only_synthetic_subjects(tmp_path: Path) -> None:
    """The marker is a claim; the sample column is the field a real export fills in."""
    repo = tmp_path / "repo"
    _init_repo(repo)
    fixture = repo / "tests" / "fixtures" / "synthetic" / "case.vcf"
    fixture.parent.mkdir(parents=True)
    marked_but_real = real_shaped_vcf().replace(
        "##reference=GRCh38", "##mva_synthetic=true\n##reference=GRCh38"
    )
    fixture.write_text(marked_but_real, encoding="utf-8")

    result = check_synthetic_fixtures_marked(repo)
    assert not result.passed
    offender = next(f for f in result.findings)
    assert "do not carry a SYNTH_ prefix" in offender.detail
    # Sample names are counted, never emitted (GP-41).
    assert "PROBAND" not in offender.detail

    fixture.write_text(marked_but_real.replace("\tPROBAND", "\tSYNTH_PROBAND"), encoding="utf-8")
    assert check_synthetic_fixtures_marked(repo).passed


def test_the_pre_commit_subset_also_scans_the_path(tmp_path: Path) -> None:
    """The staged subset is what the git hook runs; the path rule must be in it."""
    repo = tmp_path / "repo"
    _init_repo(repo)
    staged = repo / "data" / f"NHS{NHS_DIGITS}_{FAKE_SURNAME}.txt"
    staged.parent.mkdir(parents=True)
    staged.write_text("aggregate counts only\n", encoding="utf-8")
    git(repo, "add", "-f", str(staged.relative_to(repo)), check=True)

    report = run_audit(repo, staged_only=True)
    assert "git_staged_sensitive" in report.failed_checks
    assert NHS_DIGITS not in report.to_markdown()
