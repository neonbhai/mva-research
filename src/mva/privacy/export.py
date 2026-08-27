"""The public-export gate (GP-43).

Nothing leaves the workspace except through :func:`export_public_artifact`, and
nothing passes it without satisfying **all** of:

1. an explicit allowlist match — someone decided, in advance and in writing, that
   an artifact of this name is publishable;
2. ``declared is Sensitivity.PUBLIC`` — the producing stage claims it;
3. a fresh scan of the actual bytes finding no ``fail``-severity hit — the claim
   is verified against the file as it exists right now;
4. no patient-data extension on either endpoint.

The reason (2) and (3) are separate is the whole design. **Classification is a
claim; the scan is the verification.** A stage that computes aggregate QC counts
legitimately declares its output ``PUBLIC``; a later change that adds a
``worst_variant`` column to the same artifact does not update the declaration, and
nothing about the declaration would catch it. Equally, the scan alone is not
enough: it is a set of heuristics with known blind spots, and a file it happens
not to match is not thereby proven safe. Requiring both means a leak needs two
independent failures.

Refusals name the failing check and never the content: an ``ExportBlockedError``
travels to a terminal and a model context exactly like a scanner report does.
"""

from __future__ import annotations

import fnmatch
import shutil
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from mva.errors import ExportBlockedError
from mva.models.base import Sensitivity
from mva.privacy.audit import scan_bytes
from mva.privacy.classify import is_sensitive_extension
from mva.privacy.patterns import MAX_SCAN_BYTES, read_capped

#: Filenames permitted to leave the workspace. Deny by default: an artifact absent
#: from this list is refused even if its classification says PUBLIC (GP-43).
#:
#: It lives HERE, in the privacy layer, rather than only in the composition root,
#: because the allowlist is the policy and every caller of the gate must share
#: one. A caller that passes its own list — ``allowlist=(path.name,)`` was the
#: real example — has written a check that can only pass, since the file is
#: always on a list built from the file.
PUBLIC_EXPORT_ALLOWLIST: tuple[str, ...] = (
    "track1_submission.csv",
    "mechanism_report.md",
    "drug_hypotheses.md",
    "rejection_record.md",
    "track2_report.md",
)


@dataclass(frozen=True, slots=True)
class ExportDecision:
    """The gate's verdict. ``reasons`` names failing CHECKS, never content."""

    allowed: bool
    reasons: tuple[str, ...]
    scanned_rules: tuple[str, ...]


def _matches_allowlist(path: Path, allowlist: Sequence[str]) -> bool:
    """Glob match against the file name and the POSIX path.

    Both forms are tried so an allowlist can be written either as
    ``"submission.csv"`` or as ``"runs/*/submission.csv"``.
    """
    posix = path.as_posix()
    return any(
        fnmatch.fnmatch(path.name, pattern) or fnmatch.fnmatch(posix, pattern)
        for pattern in allowlist
    )


def gate_public_export(
    path: Path, *, declared: Sensitivity, allowlist: Sequence[str]
) -> ExportDecision:
    """Decide whether ``path`` may be published. Never raises for a refusal.

    Returns the full set of reasons rather than short-circuiting on the first one,
    because a caller fixing an export wants to see everything that is wrong with
    it, not to discover the problems one run at a time.
    """
    reasons: list[str] = []

    if not _matches_allowlist(path, allowlist):
        reasons.append("allowlist: file name does not match any entry in the export allowlist")
    if declared is not Sensitivity.PUBLIC:
        reasons.append(
            f"declared_sensitivity: artifact is declared {declared.value!r}, not 'public'"
        )
    if is_sensitive_extension(path):
        reasons.append("extension: path carries a patient-data extension")

    scanned_rules: tuple[str, ...] = ()
    if not path.is_file():
        reasons.append("readability: path is not a readable file")
    else:
        try:
            data = read_capped(path, MAX_SCAN_BYTES)
        except OSError as exc:
            reasons.append(f"readability: could not read the artifact ({type(exc).__name__})")
        else:
            findings = scan_bytes(data, check="export_scan", path_label=path.name)
            failing = tuple(
                sorted({f.rule_id for f in findings if f.severity == "fail" and f.rule_id})
            )
            scanned_rules = failing
            if failing:
                reasons.append(
                    "content_scan: fail-severity rule(s) matched the bytes: " + ", ".join(failing)
                )
            if path.stat().st_size > MAX_SCAN_BYTES:
                reasons.append(
                    "content_scan: artifact exceeds the scan cap, so it cannot be verified"
                )

    return ExportDecision(allowed=not reasons, reasons=tuple(reasons), scanned_rules=scanned_rules)


def export_public_artifact(
    src: Path, dest: Path, *, declared: Sensitivity, allowlist: Sequence[str]
) -> Path:
    """Copy ``src`` to ``dest`` only if the gate allows it.

    Fails closed and fails BEFORE writing anything: a refusal must not leave a
    partial copy outside the workspace, because a partial copy of a VCF is still a
    VCF. The destination is gated too — exporting to a ``.bam`` path is refused
    even when the source is clean, since the extension is what downstream tools
    and humans will trust.
    """
    decision = gate_public_export(src, declared=declared, allowlist=allowlist)
    if not decision.allowed:
        msg = (
            f"Export of {src.name!r} blocked by the public-export gate. Failing "
            f"checks: {'; '.join(decision.reasons)}. Content is withheld under GP-41."
        )
        raise ExportBlockedError(msg)
    if is_sensitive_extension(dest):
        msg = (
            f"Export destination {dest.name!r} carries a patient-data extension. "
            "Refused before writing."
        )
        raise ExportBlockedError(msg)

    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)
    return dest
