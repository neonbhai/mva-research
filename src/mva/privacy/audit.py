"""The privacy audit engine.

Every check answers one question about one concrete failure mode, and every check
reports **paths, line numbers, span lengths and rule IDs only** (GP-41). The audit
is run by an agent and by CI; its output lands in a model context window and a
build log, so a check that echoed what it found would be the largest leak in the
system.

The check set is chosen around how patient data actually escapes a research repo,
which is almost never "someone published a VCF on purpose":

* it is already **tracked**, so ``.gitignore`` is irrelevant to it
  (``git_tracked_sensitive``);
* it is **staged right now**, and the worktree copy has since been cleaned
  (``git_staged_sensitive`` reads the staged blob, not the file on disk);
* the ignore rule that was supposed to catch it is **cancelled by a negation**
  further down the file (``gitignore_effectiveness``, ``gitignore_negation_safety``);
* the workspace is a **symlink** into the repo, or a symlink in the repo points
  out into the workspace (``workspace_containment``, ``symlink_escape``);
* it is in a file whose extension says nothing (``content_scan``);
* a real file was dropped into the synthetic-fixture directory, where the
  ``.gitignore`` negations deliberately re-admit files
  (``synthetic_fixtures_marked``);
* it is in a notebook output cell (``notebook_output_purity``);
* it went to a log (``log_redaction_probe``, which is a live runtime probe, not a
  static check — the only way to know redaction works is to try to defeat it);
* the workspace is outside the repo but inside iCloud Drive
  (``cloud_sync_location``).
"""

from __future__ import annotations

import json
import logging
import os
import shutil
import subprocess
from bisect import bisect_right
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final, cast

from mva.config import (
    CLOUD_SYNCED_HOME_DIRS,
    CLOUD_SYNCED_MARKERS,
)
from mva.privacy.classify import is_sensitive_extension
from mva.privacy.patterns import (
    BINARY_SEVERITY,
    HEAD_BYTES,
    HPO_DISTINCT_FAIL_THRESHOLD,
    MAX_SCAN_BYTES,
    RULES,
    Severity,
    correlation_id,
    decode_lossy,
    read_capped,
    sniff_binary,
)
from mva.privacy.redact import GenomicRedactionFilter, install_redaction

# ---------------------------------------------------------------------------
# Policy constants
# ---------------------------------------------------------------------------

#: Tracked files here are asserted synthetic/public and audited by other checks.
TRACKED_EXEMPT_PREFIXES: Final[tuple[str, ...]] = (
    "tests/fixtures/synthetic/",
    "knowledge/public/",
)

#: Files whose *content* legitimately looks like genomic data. Findings inside
#: these are recorded at `warn` rather than suppressed, because visibility is the
#: point; the guarantee that they really are synthetic comes from
#: ``synthetic_fixtures_marked`` and from curation review of knowledge/.
CONTENT_ALLOWLIST_PREFIXES: Final[tuple[str, ...]] = (
    "tests/fixtures/synthetic/",
    "knowledge/public/",
    "knowledge/manifests/",
    "tests/golden/",
)

#: Prefixes where an HPO identifier is a CONSTANT rather than a phenotype profile:
#: reviewed source, its tests, docs and prompt briefs. This downgrade is specific
#: to ``hpo_term`` — every other fail-severity rule keeps full strength here.
#:
#: The reasoning: an HPO ID is a public ontology term. It becomes disclosive only
#: as part of a patient's *profile*, and by GP-40 a patient profile lives in the
#: external workspace and never in the repository at all — which is enforced by
#: ``workspace_containment`` and ``git_tracked_sensitive``, not by counting terms
#: in a ``.py`` file. Without this, a docstring showing the canonical ``HP:``
#: format and a test fixture listing three phenotypes both fail forever, and an
#: audit that is permanently red is an audit that gets switched off.
HPO_CODE_PREFIXES: Final[tuple[str, ...]] = ("src/", "tests/", "docs/", "prompts/")

#: The only places a ``!`` negation may re-admit a FILE.
NEGATION_ALLOWED_PREFIXES: Final[tuple[str, ...]] = (
    "tests/fixtures/synthetic/",
    "knowledge/",
    "config/",
    "templates/",
    "tests/golden/",
)

#: Directories never walked: build caches and the object store itself.
SKIP_DIR_NAMES: Final[frozenset[str]] = frozenset(
    {
        ".git",
        ".venv",
        "venv",
        "env",
        "__pycache__",
        ".pytest_cache",
        ".ruff_cache",
        ".mypy_cache",
        ".hypothesis",
        "node_modules",
        ".tox",
        ".nox",
        "htmlcov",
    }
)

#: Synthetic probe paths for ``gitignore_effectiveness``. These files are never
#: created: ``git check-ignore --no-index`` evaluates the patterns against a
#: hypothetical path, so the probe cannot itself write a sensitive-looking file
#: into the repo it is auditing.
GITIGNORE_PROBES: Final[tuple[str, ...]] = (
    "workspace/proband.vcf",
    "workspace/runs/r1/calls.vcf.gz",
    "runs/case01/out.bam",
    "runs/case01/out.cram",
    "a/b/sample_R1.fastq.gz",
    "a/b/sample_R2.fastq",
    "evidence.duckdb",
    "variants.parquet",
    "data/cohort.sqlite",
    "notes/family.ped",
    "raw/aligned.sam",
    "case.phenopacket.json",
    "secrets.yaml",
    "pipeline.log",
    "analysis.ipynb",
    ".env",
)

_MAX_FINDINGS_PER_RULE: Final[int] = 3
_MAX_MATCHES_PER_RULE: Final[int] = 500
_MAX_FINDINGS_PER_CHECK: Final[int] = 200

_EMPTY_TREE: Final[str] = "4b825dc642cb6eb9a060e54bf8d69288fbee4904"
_GIT: Final[str] = shutil.which("git") or "git"

#: Stand-in for the absolute workspace path in every emitted finding. The location
#: and naming of the patient directory is itself disclosive (``.../smith_family/``),
#: and this report is an artifact that may be committed.
WORKSPACE_LABEL: Final[str] = "$MVA_WORKSPACE"


# ---------------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Finding:
    """One observation. ``detail`` is prose about the RULE, never about the match."""

    check: str
    path: str
    line: int | None
    rule_id: str | None
    span_len: int | None
    severity: Severity
    detail: str


@dataclass(frozen=True, slots=True)
class CheckResult:
    name: str
    passed: bool
    severity: str
    findings: tuple[Finding, ...]
    summary: str


@dataclass(frozen=True, slots=True)
class AuditReport:
    results: tuple[CheckResult, ...]
    passed: bool
    failed_checks: tuple[str, ...]

    @property
    def findings(self) -> tuple[Finding, ...]:
        return tuple(f for result in self.results for f in result.findings)

    def to_markdown(self) -> str:
        """Render for a human or an agent. Contains no matched content by construction."""
        status = "PASS" if self.passed else "FAIL"
        lines: list[str] = [
            "# Privacy audit",
            "",
            f"**Result: {status}** — {len(self.results)} checks, "
            f"{len(self.failed_checks)} failed, {len(self.findings)} findings.",
            "",
            "> This report intentionally contains no matched content: only paths, "
            "line numbers, span lengths and rule IDs (GP-41).",
            "",
            "| Check | Status | Findings | Summary |",
            "| --- | --- | --- | --- |",
        ]
        for result in self.results:
            mark = "pass" if result.passed else "**FAIL**"
            summary = result.summary.replace("|", "\\|")
            lines.append(f"| `{result.name}` | {mark} | {len(result.findings)} | {summary} |")

        for result in self.results:
            if not result.findings:
                continue
            lines += ["", f"## {result.name}", ""]
            for finding in result.findings:
                location = finding.path
                if finding.line is not None:
                    location = f"{location}:{finding.line}"
                bits = [f"- `{finding.severity}` `{location}`"]
                if finding.rule_id is not None:
                    bits.append(f"rule=`{finding.rule_id}`")
                if finding.span_len is not None:
                    bits.append(f"span_len={finding.span_len}")
                bits.append(f"— {finding.detail}")
                lines.append(" ".join(bits))

        if not self.passed:
            lines += [
                "",
                "## Remediation",
                "",
                "Failing checks are listed above. Do NOT weaken a check to make it "
                "pass: fix the repository, or record a decision in `docs/decisions/` "
                "explaining why a specific finding is acceptable.",
            ]
        return "\n".join(lines) + "\n"

    def to_dict(self) -> dict[str, object]:
        return {
            "passed": self.passed,
            "failed_checks": list(self.failed_checks),
            "results": [
                {
                    "name": result.name,
                    "passed": result.passed,
                    "severity": result.severity,
                    "summary": result.summary,
                    "findings": [
                        {
                            "check": f.check,
                            "path": f.path,
                            "line": f.line,
                            "rule_id": f.rule_id,
                            "span_len": f.span_len,
                            "severity": f.severity,
                            "detail": f.detail,
                        }
                        for f in result.findings
                    ],
                }
                for result in self.results
            ],
        }


def _result(name: str, findings: Sequence[Finding], summary: str, *, strict: bool) -> CheckResult:
    """Derive pass/fail. Under ``strict`` a warning is also disqualifying."""
    trimmed = tuple(findings[:_MAX_FINDINGS_PER_CHECK])
    if len(findings) > _MAX_FINDINGS_PER_CHECK:
        summary = f"{summary} (showing {_MAX_FINDINGS_PER_CHECK} of {len(findings)} findings)"
    has_fail = any(f.severity == "fail" for f in findings)
    has_warn = any(f.severity == "warn" for f in findings)
    passed = not has_fail and not (strict and has_warn)
    severity = "fail" if has_fail else ("warn" if has_warn else "pass")
    return CheckResult(
        name=name, passed=passed, severity=severity, findings=trimmed, summary=summary
    )


# ---------------------------------------------------------------------------
# Git plumbing
# ---------------------------------------------------------------------------


def _git(repo_root: Path, args: Sequence[str], *, stdin: bytes | None = None) -> tuple[int, bytes]:
    """Run a git plumbing command. Never raises on non-zero; callers decide."""
    try:
        proc = subprocess.run(  # noqa: S603
            [_GIT, "-C", str(repo_root), *args],
            capture_output=True,
            input=stdin,
            check=False,
        )
    except OSError:
        return 127, b""
    return proc.returncode, proc.stdout


def _split_nul(data: bytes) -> list[str]:
    return [chunk.decode("utf-8", errors="replace") for chunk in data.split(b"\0") if chunk]


def _git_available(repo_root: Path) -> bool:
    code, _ = _git(repo_root, ["rev-parse", "--git-dir"])
    return code == 0


def _exempt(path: str, prefixes: Sequence[str]) -> bool:
    return any(path.startswith(prefix) for prefix in prefixes)


def _unavailable(name: str, reason: str, *, strict: bool) -> CheckResult:
    finding = Finding(
        check=name,
        path=".",
        line=None,
        rule_id=None,
        span_len=None,
        severity="warn",
        detail=reason,
    )
    return _result(name, [finding], f"Skipped: {reason}", strict=strict)


# ---------------------------------------------------------------------------
# Content scanning
# ---------------------------------------------------------------------------


def _line_offsets(data: bytes) -> list[int]:
    offsets: list[int] = []
    start = 0
    while True:
        index = data.find(b"\n", start)
        if index == -1:
            return offsets
        offsets.append(index)
        start = index + 1


def _line_of(offsets: Sequence[int], position: int) -> int:
    return bisect_right(offsets, position - 1) + 1


def _resolve_context(
    *,
    rule_id: str,
    base_severity: Severity,
    description: str,
    matched_ids: set[str],
    distinct_hpo: int,
    allowlisted: bool,
    hpo_is_constant: bool,
) -> tuple[Severity, str]:
    """Refine a rule's isolated severity using whole-file context.

    Kept separate from :func:`scan_bytes` because this is the part reviewers should
    argue about: it is where a `warn` becomes a `fail`, and where a documented
    false-positive class is held down. Every branch states its reason in the
    returned note, so the report explains its own severities.
    """
    severity: Severity = base_severity
    note = description
    if rule_id == "vcf_data_line":
        promoted = bool({"vcf_header", "genotype_field"} & matched_ids)
        severity = "fail" if promoted else "warn"
        note = (
            f"{note} Promoted to fail: the file also matches a VCF header or a genotype field."
            if promoted
            else f"{note} Held at warn: no VCF header or genotype field in this file."
        )
    elif rule_id == "hpo_term":
        over_threshold = distinct_hpo >= HPO_DISTINCT_FAIL_THRESHOLD
        severity = "fail" if over_threshold and not hpo_is_constant else "warn"
        note = f"{note} {distinct_hpo} distinct term(s) in this file."
        if over_threshold and hpo_is_constant:
            note = f"{note} Held at warn: reviewed source/docs, where HPO IDs are constants."
    if allowlisted and severity == "fail":
        severity = "warn"
        note = f"{note} Downgraded: path is on the audited public/synthetic allowlist."
    return severity, note


def scan_bytes(
    data: bytes,
    *,
    check: str,
    path_label: str,
    allowlisted: bool = False,
    hpo_is_constant: bool = False,
) -> list[Finding]:
    """Apply the rule battery to one buffer and return content-free findings.

    Two rules are resolved here rather than in :mod:`mva.privacy.patterns`, because
    only this function has whole-file context:

    * ``vcf_data_line`` is promoted from ``warn`` to ``fail`` when the same file
      also matches ``vcf_header`` or ``genotype_field``. A tab-separated
      coordinate table is a public resource; the same table next to a VCF header
      is a patient's call set.
    * ``hpo_term`` fails only at :data:`HPO_DISTINCT_FAIL_THRESHOLD` DISTINCT terms
      in a non-allowlisted file. Distinctness is counted over
      :func:`~mva.privacy.patterns.correlation_id` values, so the terms themselves
      are never accumulated in a data structure that could later be printed.
    """
    offsets = _line_offsets(data)
    spans_by_rule: dict[str, list[tuple[int, int]]] = {}
    distinct_hpo = 0

    for rule in RULES:
        spans: list[tuple[int, int]] = []
        seen: set[str] = set()
        for match in rule.pattern.finditer(data):
            spans.append(match.span())
            if rule.rule_id == "hpo_term":
                # The matched bytes exist only as an argument to the salted HMAC;
                # they are never bound to a name that outlives this expression.
                seen.add(correlation_id(match.group(0)))
            if len(spans) >= _MAX_MATCHES_PER_RULE:
                break
        if spans:
            spans_by_rule[rule.rule_id] = spans
            if rule.rule_id == "hpo_term":
                distinct_hpo = len(seen)

    matched_ids = set(spans_by_rule)
    findings: list[Finding] = []

    for rule in RULES:
        spans = spans_by_rule.get(rule.rule_id, [])
        if not spans:
            continue
        severity, note = _resolve_context(
            rule_id=rule.rule_id,
            base_severity=rule.severity,
            description=rule.description,
            matched_ids=matched_ids,
            distinct_hpo=distinct_hpo,
            allowlisted=allowlisted,
            hpo_is_constant=hpo_is_constant,
        )
        total = len(spans)
        for start, end in spans[:_MAX_FINDINGS_PER_RULE]:
            findings.append(
                Finding(
                    check=check,
                    path=path_label,
                    line=_line_of(offsets, start),
                    rule_id=rule.rule_id,
                    span_len=end - start,
                    severity=severity,
                    detail=f"{note} ({total} occurrence(s) in this file.)",
                )
            )

    kind = sniff_binary(data[:HEAD_BYTES])
    if kind is not None:
        severity = BINARY_SEVERITY[kind]
        if allowlisted and severity == "fail":
            severity = "warn"
        findings.append(
            Finding(
                check=check,
                path=path_label,
                line=None,
                rule_id=f"magic:{kind}",
                span_len=None,
                severity=severity,
                detail=(
                    f"Container identified from magic bytes as {kind!r}. "
                    "Extension checks are defeated by a rename; this is not."
                ),
            )
        )
    return findings


def scan_file(
    path: Path, *, check: str, path_label: str, allowlisted: bool, hpo_is_constant: bool = False
) -> list[Finding]:
    """Scan one file on disk, capped, in bytes, never decoding for detection."""
    try:
        size = path.stat().st_size
    except OSError:
        return []
    try:
        data = read_capped(path, MAX_SCAN_BYTES)
    except OSError as exc:
        return [
            Finding(
                check=check,
                path=path_label,
                line=None,
                rule_id=None,
                span_len=None,
                severity="warn",
                detail=f"Unreadable ({type(exc).__name__}); not scanned.",
            )
        ]
    findings = scan_bytes(
        data,
        check=check,
        path_label=path_label,
        allowlisted=allowlisted,
        hpo_is_constant=hpo_is_constant,
    )
    if size > MAX_SCAN_BYTES:
        findings.append(
            Finding(
                check=check,
                path=path_label,
                line=None,
                rule_id=None,
                span_len=None,
                severity="warn",
                detail=(
                    f"File is {size} bytes; only the first {MAX_SCAN_BYTES} were scanned. "
                    "Magic-byte identification still applied to the head."
                ),
            )
        )
    return findings


# ---------------------------------------------------------------------------
# Workspace resolution shared by two checks
# ---------------------------------------------------------------------------


def _configured_workspace(workspace: Path | None) -> Path | None:
    raw: str | Path | None = workspace if workspace is not None else os.environ.get("MVA_WORKSPACE")
    if raw is None:
        return None
    try:
        return Path(raw).expanduser().resolve()
    except OSError:
        return None


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------


def check_git_tracked_sensitive(repo_root: Path, *, strict: bool = False) -> CheckResult:
    """Tracked files with a sensitive extension.

    This is the check that matters most, and the one people assume ``.gitignore``
    already covers. It does not: ignore rules are consulted only for *untracked*
    files. Once a path is in the index, git tracks it forever regardless of any
    pattern, and `git rm --cached` leaves it recoverable in history.
    """
    name = "git_tracked_sensitive"
    if not _git_available(repo_root):
        return _unavailable(name, "not a git repository (or git unavailable)", strict=strict)
    code, out = _git(repo_root, ["ls-files", "-z"])
    if code != 0:
        return _unavailable(name, "git ls-files failed", strict=strict)
    paths = _split_nul(out)
    findings = [
        Finding(
            check=name,
            path=path,
            line=None,
            rule_id=None,
            span_len=None,
            severity="fail",
            detail=(
                "Tracked file has a patient-data extension. .gitignore does not "
                "apply to tracked paths. Remove it from the index and purge history."
            ),
        )
        for path in sorted(paths)
        if is_sensitive_extension(Path(path)) and not _exempt(path, TRACKED_EXEMPT_PREFIXES)
    ]
    return _result(
        name,
        findings,
        f"{len(paths)} tracked paths inspected; {len(findings)} with sensitive extensions.",
        strict=strict,
    )


def check_git_staged_sensitive(repo_root: Path, *, strict: bool = False) -> CheckResult:
    """Sensitive extensions or content in the STAGED blobs.

    Deliberately reads ``git cat-file blob`` rather than the worktree. The staged
    content is what a commit will contain; a file can be staged dirty and then
    cleaned on disk, and only the index still holds the payload. This is the
    pre-commit hook's check.
    """
    name = "git_staged_sensitive"
    if not _git_available(repo_root):
        return _unavailable(name, "not a git repository (or git unavailable)", strict=strict)

    args = ["diff", "--cached", "--name-only", "-z", "--diff-filter=ACMR"]
    code, out = _git(repo_root, args)
    if code != 0:
        # No HEAD yet: diff against the empty tree so the first commit is audited too.
        code, out = _git(repo_root, [*args, _EMPTY_TREE])
        if code != 0:
            return _unavailable(name, "git diff --cached failed", strict=strict)

    paths = sorted(_split_nul(out))
    findings: list[Finding] = []
    for path in paths:
        allowlisted = _exempt(path, CONTENT_ALLOWLIST_PREFIXES)
        if is_sensitive_extension(Path(path)) and not _exempt(path, TRACKED_EXEMPT_PREFIXES):
            findings.append(
                Finding(
                    check=name,
                    path=path,
                    line=None,
                    rule_id=None,
                    span_len=None,
                    severity="fail",
                    detail="Staged file has a patient-data extension. Unstage it.",
                )
            )
        sha_code, sha_out = _git(repo_root, ["rev-parse", f":{path}"])
        if sha_code != 0:
            continue
        sha = sha_out.decode("ascii", errors="replace").strip()
        blob_code, blob = _git(repo_root, ["cat-file", "blob", sha])
        if blob_code != 0:
            continue
        findings.extend(
            scan_bytes(
                blob[:MAX_SCAN_BYTES],
                check=name,
                path_label=path,
                allowlisted=allowlisted,
                hpo_is_constant=_exempt(path, HPO_CODE_PREFIXES),
            )
        )
    return _result(name, findings, f"{len(paths)} staged blob(s) inspected.", strict=strict)


def _parse_check_ignore(line: str) -> tuple[str, str] | None:
    """Parse one ``git check-ignore -v`` line into (source_location, pattern)."""
    if "\t" not in line:
        return None
    left, _, _ = line.partition("\t")
    parts = left.split(":", 2)
    if len(parts) != 3:
        return None
    source, lineno, pattern = parts
    return f"{source}:{lineno}", pattern


def check_gitignore_effectiveness(repo_root: Path, *, strict: bool = False) -> CheckResult:
    """Does ``.gitignore`` actually ignore the paths it is supposed to?

    The subtlety that makes this check worth writing: ``git check-ignore`` exits 0
    when *a rule matched*, and a NEGATION (``!pattern``) is a rule. So exit 0 alone
    means "some rule had an opinion", not "the file is ignored". A pass therefore
    requires exit 0 **and** a reported pattern that does not begin with ``!``.

    ``--no-index`` evaluates the patterns against a hypothetical path, so no probe
    file is ever written — which matters, because writing ``workspace/proband.vcf``
    to test whether it would be ignored is exactly the accident being guarded
    against.
    """
    name = "gitignore_effectiveness"
    if not _git_available(repo_root):
        return _unavailable(name, "not a git repository (or git unavailable)", strict=strict)

    findings: list[Finding] = []
    for probe in GITIGNORE_PROBES:
        code, out = _git(repo_root, ["check-ignore", "-v", "--no-index", "--", probe])
        text = out.decode("utf-8", errors="replace").strip()
        if code != 0 or not text:
            findings.append(
                Finding(
                    check=name,
                    path=probe,
                    line=None,
                    rule_id=None,
                    span_len=None,
                    severity="fail",
                    detail=(
                        "No .gitignore rule matches this probe path; a real file here "
                        "would be committable."
                    ),
                )
            )
            continue
        parsed = _parse_check_ignore(text.splitlines()[0])
        if parsed is None:
            findings.append(
                Finding(
                    check=name,
                    path=probe,
                    line=None,
                    rule_id=None,
                    span_len=None,
                    severity="warn",
                    detail="Unparseable `git check-ignore -v` output.",
                )
            )
            continue
        source, pattern = parsed
        if pattern.startswith("!"):
            findings.append(
                Finding(
                    check=name,
                    path=probe,
                    line=None,
                    rule_id=None,
                    span_len=None,
                    severity="fail",
                    detail=(
                        f"Matched by a NEGATION rule at {source} — the path is "
                        "RE-INCLUDED, not ignored. Exit code 0 from check-ignore does "
                        "not mean ignored."
                    ),
                )
            )
    return _result(
        name,
        findings,
        f"{len(GITIGNORE_PROBES)} synthetic probe paths evaluated with --no-index.",
        strict=strict,
    )


def check_gitignore_negation_safety(repo_root: Path, *, strict: bool = False) -> CheckResult:
    """No ``!`` rule may re-admit a file outside the audited public directories.

    A negation is the only way to punch a hole in a deny-by-default ignore file,
    so the blast radius of each one is reviewed here rather than trusted.

    One carve-out, stated explicitly because it is a real weakening: a *pure
    directory* re-inclusion (``!dir/`` or ``!dir/**/``) that is an ancestor of an
    allowed prefix is permitted. Git will not descend into an excluded directory,
    so those lines are structurally required for the narrow file negations beneath
    them to have any effect — and a directory negation cannot by itself re-admit a
    file, because the file patterns still apply inside.
    """
    name = "gitignore_negation_safety"
    gitignore = repo_root / ".gitignore"
    if not gitignore.is_file():
        return _unavailable(name, "no .gitignore at the repository root", strict=strict)

    findings: list[Finding] = []
    text = decode_lossy(read_capped(gitignore, MAX_SCAN_BYTES))
    negations = 0
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line.startswith("!"):
            continue
        negations += 1
        target = line[1:].strip().lstrip("/")
        if _exempt(target, NEGATION_ALLOWED_PREFIXES):
            continue
        directory = target.removesuffix("**/").removesuffix("**").rstrip("/")
        is_directory_rule = target.endswith("/")
        if is_directory_rule and any(
            allowed.startswith(f"{directory}/") for allowed in NEGATION_ALLOWED_PREFIXES
        ):
            continue
        findings.append(
            Finding(
                check=name,
                path=".gitignore",
                line=lineno,
                rule_id=None,
                span_len=None,
                severity="fail",
                detail=(
                    f"Negation `{line}` re-admits a path outside the audited set "
                    f"({', '.join(NEGATION_ALLOWED_PREFIXES)}). Every negation is a "
                    "hole in a deny-by-default policy and needs a decision record."
                ),
            )
        )
    return _result(name, findings, f"{negations} negation rule(s) reviewed.", strict=strict)


def check_workspace_containment(
    repo_root: Path, *, workspace: Path | None = None, strict: bool = False
) -> CheckResult:
    """GP-40: the patient workspace must resolve outside the repository.

    Symlinks are resolved first. A workspace that is a symlink into the repo tree
    passes a naive string comparison and fails this one, which is the entire point:
    data written through the link lands inside the repo and is one ``git add -A``
    from being committed and permanently recoverable.

    The absolute path is never echoed — the workspace directory name is itself
    disclosive — so findings report the relationship, not the location.
    """
    name = "workspace_containment"
    resolved = _configured_workspace(workspace)
    if resolved is None:
        return _unavailable(
            name,
            "no workspace configured (pass workspace= or set MVA_WORKSPACE)",
            strict=strict,
        )
    repo = repo_root.resolve()
    findings: list[Finding] = []
    if resolved == repo or resolved.is_relative_to(repo):
        findings.append(
            Finding(
                check=name,
                path=WORKSPACE_LABEL,
                line=None,
                rule_id=None,
                span_len=None,
                severity="fail",
                detail=(
                    "Workspace resolves INSIDE the repository (symlinks followed). "
                    "Move it out; data written there is one `git add -A` from the "
                    "index and stays recoverable from history afterwards."
                ),
            )
        )
    if not resolved.exists():
        findings.append(
            Finding(
                check=name,
                path=WORKSPACE_LABEL,
                line=None,
                rule_id=None,
                span_len=None,
                severity="warn",
                detail="Configured workspace does not exist yet.",
            )
        )
    return _result(name, findings, "Workspace containment evaluated.", strict=strict)


def check_symlink_escape(
    repo_root: Path, *, workspace: Path | None = None, strict: bool = False
) -> CheckResult:
    """Symlinks in the repo that reach into the workspace or at sensitive data.

    Walked with ``followlinks=False``: following them would both loop and,
    ironically, read the patient data this check exists to keep at arm's length.
    A link is a leak because git stores the *link*, but every tool that reads the
    repo — including an agent — reads through it.
    """
    name = "symlink_escape"
    repo = repo_root.resolve()
    resolved_workspace = _configured_workspace(workspace)
    findings: list[Finding] = []

    for dirpath, dirnames, filenames in os.walk(repo, followlinks=False):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIR_NAMES]
        base = Path(dirpath)
        for entry in sorted(dirnames) + sorted(filenames):
            link = base / entry
            if not link.is_symlink():
                continue
            label = link.relative_to(repo).as_posix()
            target = Path(os.path.realpath(link))
            if resolved_workspace is not None and (
                target == resolved_workspace or target.is_relative_to(resolved_workspace)
            ):
                findings.append(
                    Finding(
                        check=name,
                        path=label,
                        line=None,
                        rule_id=None,
                        span_len=None,
                        severity="fail",
                        detail=(
                            f"Symlink resolves into {WORKSPACE_LABEL}. The repo now "
                            "reads patient data through a tracked path."
                        ),
                    )
                )
            elif is_sensitive_extension(target):
                findings.append(
                    Finding(
                        check=name,
                        path=label,
                        line=None,
                        rule_id=None,
                        span_len=None,
                        severity="fail",
                        detail="Symlink resolves to a path with a patient-data extension.",
                    )
                )
            elif not target.is_relative_to(repo):
                findings.append(
                    Finding(
                        check=name,
                        path=label,
                        line=None,
                        rule_id=None,
                        span_len=None,
                        severity="warn",
                        detail="Symlink resolves outside the repository.",
                    )
                )
    return _result(name, findings, "Repository walked for symlinks.", strict=strict)


def _scannable_paths(repo_root: Path) -> list[str]:
    """Tracked plus untracked-but-not-ignored paths, deduplicated and sorted.

    Ignored files are excluded on purpose: an ignored patient VCF sitting in a
    developer's worktree is expected and is not this check's business. The checks
    that DO care about it are ``workspace_containment`` and ``gitignore_effectiveness``.
    """
    seen: set[str] = set()
    for args in (["ls-files", "-z"], ["ls-files", "-z", "--others", "--exclude-standard"]):
        code, out = _git(repo_root, args)
        if code == 0:
            seen.update(_split_nul(out))
    return sorted(seen)


def check_content_scan(repo_root: Path, *, strict: bool = False) -> CheckResult:
    """The rule battery over every file git would let you commit."""
    name = "content_scan"
    if not _git_available(repo_root):
        return _unavailable(name, "not a git repository (or git unavailable)", strict=strict)
    paths = _scannable_paths(repo_root)
    findings: list[Finding] = []
    scanned = 0
    for path in paths:
        candidate = repo_root / path
        if not candidate.is_file() or candidate.is_symlink():
            continue
        if any(part in SKIP_DIR_NAMES for part in Path(path).parts):
            continue
        scanned += 1
        findings.extend(
            scan_file(
                candidate,
                check=name,
                path_label=path,
                allowlisted=_exempt(path, CONTENT_ALLOWLIST_PREFIXES),
                hpo_is_constant=_exempt(path, HPO_CODE_PREFIXES),
            )
        )
    fails = sum(1 for f in findings if f.severity == "fail")
    return _result(
        name,
        findings,
        f"{scanned} file(s) scanned (cap {MAX_SCAN_BYTES} bytes); {fails} fail-severity hit(s).",
        strict=strict,
    )


def check_synthetic_fixtures_marked(repo_root: Path, *, strict: bool = False) -> CheckResult:
    """Everything under ``tests/fixtures/synthetic/`` must declare itself synthetic.

    This directory is the one place ``.gitignore`` deliberately re-admits ``.vcf``,
    ``.ped`` and ``.bed`` files. Without this check, the safest-looking way to leak
    a real patient VCF is to save it here, because every other control assumes the
    directory is synthetic. The marker requirement inverts that: an unmarked file
    FAILS rather than sliding in.
    """
    name = "synthetic_fixtures_marked"
    root = repo_root / "tests" / "fixtures" / "synthetic"
    if not root.is_dir():
        return _result(name, [], "No synthetic fixture directory present.", strict=strict)

    markers = (b"mva_synthetic=true", b"SYNTH_", b"SYNTH-")
    findings: list[Finding] = []
    checked = 0
    for candidate in sorted(root.rglob("*")):
        if not candidate.is_file() or candidate.name == ".DS_Store":
            continue
        checked += 1
        label = candidate.relative_to(repo_root).as_posix()
        try:
            head = read_capped(candidate, MAX_SCAN_BYTES)
        except OSError:
            head = b""
        if not any(marker in head for marker in markers):
            findings.append(
                Finding(
                    check=name,
                    path=label,
                    line=None,
                    rule_id=None,
                    span_len=None,
                    severity="fail",
                    detail=(
                        "No synthetic marker found. Add `mva_synthetic=true` (a VCF "
                        "header line, a leading comment) or use SYNTH_-prefixed sample "
                        "IDs. If the file is not synthetic it must not be here at all."
                    ),
                )
            )
    return _result(name, findings, f"{checked} fixture file(s) checked.", strict=strict)


def check_notebook_output_purity(repo_root: Path, *, strict: bool = False) -> CheckResult:
    """No tracked notebook may carry executed outputs.

    A notebook output cell is a verbatim dump of whatever the cell printed: a
    DataFrame head of variant rows, a rendered pedigree, a stack trace with frame
    locals. It is also invisible in a normal diff review.
    """
    name = "notebook_output_purity"
    if not _git_available(repo_root):
        return _unavailable(name, "not a git repository (or git unavailable)", strict=strict)
    code, out = _git(repo_root, ["ls-files", "-z", "--", "*.ipynb"])
    if code != 0:
        return _unavailable(name, "git ls-files failed", strict=strict)

    findings: list[Finding] = []
    paths = sorted(_split_nul(out))
    for path in paths:
        candidate = repo_root / path
        if not candidate.is_file():
            continue
        try:
            document = cast(dict[str, Any], json.loads(read_capped(candidate, MAX_SCAN_BYTES)))
        except (ValueError, UnicodeDecodeError) as exc:
            findings.append(
                Finding(
                    check=name,
                    path=path,
                    line=None,
                    rule_id=None,
                    span_len=None,
                    severity="warn",
                    detail=f"Not parseable as a notebook ({type(exc).__name__}).",
                )
            )
            continue
        cells = document.get("cells")
        if not isinstance(cells, list):
            continue
        for index, cell in enumerate(cast(list[Any], cells)):
            if not isinstance(cell, dict):
                continue
            outputs = cast(dict[str, Any], cell).get("outputs")
            if isinstance(outputs, list) and outputs:
                findings.append(
                    Finding(
                        check=name,
                        path=path,
                        line=None,
                        rule_id=None,
                        span_len=None,
                        severity="fail",
                        detail=(
                            f"Cell {index} carries {len(cast(list[Any], outputs))} stored "
                            "output(s). Strip outputs before committing."
                        ),
                    )
                )
    return _result(name, findings, f"{len(paths)} tracked notebook(s) checked.", strict=strict)


# ---------------------------------------------------------------------------
# Runtime probe (GP-42)
#
# The canary payloads below are ASSEMBLED AT RUNTIME from fragments. Written as
# plain literals they would make this very file trip `content_scan`, and lowering
# the rules to accommodate the audit's own source would be exactly the kind of
# self-serving weakening the audit exists to prevent. Assembly also documents the
# rules' shape: a bare position or a bare digit run is intentionally not
# detectable, so the fragments alone are inert.
# ---------------------------------------------------------------------------

_CANARY_POS: Final[str] = "40200000"
_CANARY_ID: Final[str] = "4457821"
_CANARY_HPO_DIGITS: Final[str] = "0001250"
_CANARY_SEQ: Final[str] = "ACGT" * 12


def canary_payloads() -> tuple[tuple[str, str, tuple[str, ...]], ...]:
    """(label, payload, substrings that must NOT survive redaction)."""
    vcf_line = "\t".join(
        ("chr15", _CANARY_POS, ".", "C", "T", "820.5", "PASS", "SYNTH", "GT:DP", "0/1:45")
    )
    fastq = "\n".join(("@synthetic_read_1", _CANARY_SEQ, "+", "I" * len(_CANARY_SEQ)))
    hpo = f"HP:{_CANARY_HPO_DIGITS}"
    mrn = f"MRN: {_CANARY_ID}"
    return (
        ("vcf_data_line", vcf_line, (f"chr15\t{_CANARY_POS}", _CANARY_POS)),
        ("fastq_record", fastq, (_CANARY_SEQ,)),
        ("hpo_term", hpo, (hpo,)),
        ("mrn", mrn, (_CANARY_ID,)),
    )


class _CaptureHandler(logging.Handler):
    """Collects what a handler would have written, without writing it."""

    def __init__(self, sink: list[str]) -> None:
        super().__init__(level=logging.DEBUG)
        self._sink = sink

    def emit(self, record: logging.LogRecord) -> None:
        try:
            self._sink.append(self.format(record))
        except Exception:  # a formatting failure must not abort the audit
            self._sink.append("<format-error>")


def check_log_redaction_probe(*, strict: bool = False) -> CheckResult:
    """Push real canaries through the CONFIGURED logging stack and look for survivors.

    This is a live probe rather than a static assertion, because every static
    argument for "our logs are clean" has failed in practice: the filter was
    attached to a logger instead of a handler, a library added its own handler
    after configuration, a formatter cached ``exc_text``. The only reliable test is
    to log something that must not appear and then read every handler's output.

    Each handler's ``emit`` is temporarily swapped for a capture that formats the
    record and discards it. Filters still run — they are applied in
    ``Handler.handle`` before ``emit`` — so this measures the real pipeline while
    keeping the canaries off the terminal, out of the log file, and out of any
    model context reading this run.
    """
    name = "log_redaction_probe"
    canaries = canary_payloads()
    captured: list[str] = []
    root = logging.getLogger()
    probe = _CaptureHandler(captured)

    original_level = root.level
    swapped: list[tuple[logging.Handler, bool, Any]] = []

    root.addHandler(probe)
    for handler in list(root.handlers):
        if handler is probe:
            continue
        had_own = "emit" in handler.__dict__
        swapped.append((handler, had_own, handler.__dict__.get("emit")))

        def _capture_for(bound: logging.Handler) -> Callable[[logging.LogRecord], None]:
            def _emit(record: logging.LogRecord) -> None:
                try:
                    captured.append(bound.format(record))
                except Exception:  # never let a formatter abort the probe
                    captured.append("<format-error>")

            return _emit

        cast(Any, handler).emit = _capture_for(handler)

    install_redaction()
    unfiltered = [
        type(handler).__name__
        for handler in root.handlers
        if not any(isinstance(f, GenomicRedactionFilter) for f in handler.filters)
    ]

    logger = logging.getLogger("mva.privacy.audit.canary")
    logger.setLevel(logging.NOTSET)
    logger.propagate = True
    root.setLevel(logging.DEBUG)
    try:
        for _, payload, _ in canaries:
            # Both paths matter: %-args are scrubbed separately from record.msg.
            logger.debug("canary via args: %s", payload)
            logger.debug(payload)
    finally:
        root.setLevel(original_level)
        root.removeHandler(probe)
        for handler, had_own, previous in swapped:
            if had_own:
                cast(Any, handler).emit = previous
            else:
                handler.__dict__.pop("emit", None)

    haystack = "\n".join(captured)
    findings: list[Finding] = [
        Finding(
            check=name,
            path="<logging>",
            line=None,
            rule_id=rule_id,
            span_len=len(needle),
            severity="fail",
            detail=(
                "A canary substring survived the logging pipeline. Redaction is not "
                "in force for at least one handler."
            ),
        )
        for rule_id, _, needles in canaries
        for needle in needles
        if needle in haystack
    ]
    if not captured:
        findings.append(
            Finding(
                check=name,
                path="<logging>",
                line=None,
                rule_id=None,
                span_len=None,
                severity="warn",
                detail="No handler output captured; the probe was inconclusive.",
            )
        )
    findings.extend(
        Finding(
            check=name,
            path="<logging>",
            line=None,
            rule_id=None,
            span_len=None,
            severity="warn",
            detail=f"Root handler {handler_name} carries no GenomicRedactionFilter.",
        )
        for handler_name in unfiltered
    )
    return _result(
        name,
        findings,
        f"{len(canaries)} canary payload(s) pushed at DEBUG; "
        f"{len(captured)} handler emission(s) captured.",
        strict=strict,
    )


def check_cloud_sync_location(
    repo_root: Path, *, workspace: Path | None = None, strict: bool = False
) -> CheckResult:
    """The workspace must not sit under a cloud-synced root.

    ``~/Desktop`` and ``~/Documents`` are included because macOS syncs both to
    iCloud Drive by default ("Desktop & Documents Folders"). A VCF written there is
    uploaded within seconds, is outside the researcher's control from that moment,
    and cannot be recalled. Findings name the matched marker, never the path.
    """
    name = "cloud_sync_location"
    findings: list[Finding] = []
    home = Path.home().resolve()

    resolved = _configured_workspace(workspace)
    if resolved is not None:
        posix = resolved.as_posix()
        findings.extend(
            Finding(
                check=name,
                path=WORKSPACE_LABEL,
                line=None,
                rule_id=None,
                span_len=None,
                severity="fail",
                detail=(
                    f"Workspace is under a cloud-synced root ({marker!r}). Patient "
                    "data placed there is uploaded to a third party automatically."
                ),
            )
            for marker in CLOUD_SYNCED_MARKERS
            if marker in posix
        )
        for dirname in CLOUD_SYNCED_HOME_DIRS:
            synced = home / dirname
            if resolved == synced or resolved.is_relative_to(synced):
                findings.append(
                    Finding(
                        check=name,
                        path=WORKSPACE_LABEL,
                        line=None,
                        rule_id=None,
                        span_len=None,
                        severity="fail",
                        detail=(
                            f"Workspace is under ~/{dirname}, which macOS syncs to "
                            "iCloud Drive by default. Use a non-synced directory."
                        ),
                    )
                )

    repo_posix = repo_root.resolve().as_posix()
    findings.extend(
        Finding(
            check=name,
            path=".",
            line=None,
            rule_id=None,
            span_len=None,
            severity="warn",
            detail=(
                f"The repository itself is under a cloud-synced root ({marker!r}). "
                "Audit reports and logs written here leave the machine."
            ),
        )
        for marker in CLOUD_SYNCED_MARKERS
        if marker in repo_posix
    )
    summary = (
        "Workspace and repository checked against cloud-sync markers."
        if resolved is not None
        else "No workspace configured; only the repository was checked."
    )
    return _result(name, findings, summary, strict=strict)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

#: Checks that inspect the index/worktree only — the pre-commit subset.
STAGED_CHECKS: Final[tuple[str, ...]] = (
    "git_staged_sensitive",
    "gitignore_effectiveness",
    "gitignore_negation_safety",
    "synthetic_fixtures_marked",
    "notebook_output_purity",
    "log_redaction_probe",
)


def run_audit(
    repo_root: Path,
    *,
    workspace: Path | None = None,
    staged_only: bool = False,
    strict: bool = False,
) -> AuditReport:
    """Run the privacy audit.

    ``staged_only`` runs the fast pre-commit subset (:data:`STAGED_CHECKS`): what is
    about to be committed, plus the ignore-file and logging invariants. It skips the
    whole-tree content scan and the symlink walk, which are the slow checks and
    which say nothing about *this* commit.

    ``strict`` promotes warnings to failures. Off by default because the
    highest-value warnings (``vcf_data_line`` in a public coordinate table, ISO
    dates everywhere) are dominated by legitimate matches; a strict run is for CI
    on a release branch, not for the inner loop.

    ``workspace`` defaults to ``$MVA_WORKSPACE``. When neither is set, the two
    workspace checks report a warning and are skipped rather than failing: an
    audit of the repository alone is a legitimate thing to want.
    """
    root = repo_root.resolve()
    results: list[CheckResult] = []
    if staged_only:
        results = [
            check_git_staged_sensitive(root, strict=strict),
            check_gitignore_effectiveness(root, strict=strict),
            check_gitignore_negation_safety(root, strict=strict),
            check_synthetic_fixtures_marked(root, strict=strict),
            check_notebook_output_purity(root, strict=strict),
            check_log_redaction_probe(strict=strict),
        ]
    else:
        results = [
            check_git_tracked_sensitive(root, strict=strict),
            check_git_staged_sensitive(root, strict=strict),
            check_gitignore_effectiveness(root, strict=strict),
            check_gitignore_negation_safety(root, strict=strict),
            check_workspace_containment(root, workspace=workspace, strict=strict),
            check_symlink_escape(root, workspace=workspace, strict=strict),
            check_content_scan(root, strict=strict),
            check_synthetic_fixtures_marked(root, strict=strict),
            check_notebook_output_purity(root, strict=strict),
            check_log_redaction_probe(strict=strict),
            check_cloud_sync_location(root, workspace=workspace, strict=strict),
        ]

    failed = tuple(result.name for result in results if not result.passed)
    return AuditReport(results=tuple(results), passed=not failed, failed_checks=failed)
