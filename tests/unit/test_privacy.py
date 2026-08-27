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

from mva.config import find_repo_root, resolve_workspace
from mva.errors import ExportBlockedError, NetworkDeniedError, WorkspaceError
from mva.models.base import Sensitivity
from mva.privacy.audit import (
    check_gitignore_effectiveness,
    check_symlink_escape,
    run_audit,
    scan_bytes,
)
from mva.privacy.classify import classify_path, is_sensitive_extension
from mva.privacy.export import export_public_artifact, gate_public_export
from mva.privacy.netguard import OfflineProfile, is_armed
from mva.privacy.patterns import RULES, rule_by_id, sniff_binary
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
    report = run_audit(find_repo_root(), staged_only=True)
    probe = next(r for r in report.results if r.name == "log_redaction_probe")
    assert probe.passed, [f.detail for f in probe.findings]
