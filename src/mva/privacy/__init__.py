"""Privacy enforcement: the controls that make GP-40..GP-44 checkable.

Layer 3 in the GP-01 stack. It may import ``models``, the foundation utilities and
``config``, and nothing above.

The five concerns, and the one rule that binds them:

* :mod:`~mva.privacy.patterns` — the shared detection battery (bytes regexes plus
  container magic-byte sniffers), each rule carrying an explicit false-positive
  note.
* :mod:`~mva.privacy.redact` — text redaction and a logging pipeline that cannot
  carry a genomic record (GP-42).
* :mod:`~mva.privacy.netguard` — an offline tripwire for the patient-data path,
  with its limits stated honestly in the module docstring.
* :mod:`~mva.privacy.classify` — cheap, fail-closed path classification.
* :mod:`~mva.privacy.audit` — the check engine.
* :mod:`~mva.privacy.export` — the allowlist-and-rescan gate (GP-43).

**GP-41 binds all of them: nothing in this package ever emits what it matched.**
Findings carry a path, a line number, a span length and a rule ID. The scanner's
own output is a leak vector, because an agent runs the audit and the result enters
a model context, a CI log and a terminal scrollback.
"""

from __future__ import annotations

from mva.privacy.audit import (
    AuditReport,
    CheckResult,
    Finding,
    run_audit,
)
from mva.privacy.classify import (
    SENSITIVE_EXTENSIONS,
    classify_path,
    is_sensitive_extension,
)
from mva.privacy.export import (
    PUBLIC_EXPORT_ALLOWLIST,
    ExportDecision,
    export_public_artifact,
    gate_public_export,
)
from mva.privacy.netguard import (
    OfflineProfile,
    arm_audit_hook,
    configure_reference_cache,
    is_armed,
)
from mva.privacy.patterns import RULES, Rule, rule_by_id, sniff_binary
from mva.privacy.redact import (
    GenomicRedactionFilter,
    install_redaction,
    redact_text,
    safe_repr,
)

__all__ = [
    "PUBLIC_EXPORT_ALLOWLIST",
    "RULES",
    "SENSITIVE_EXTENSIONS",
    "AuditReport",
    "CheckResult",
    "ExportDecision",
    "Finding",
    "GenomicRedactionFilter",
    "OfflineProfile",
    "Rule",
    "arm_audit_hook",
    "classify_path",
    "configure_reference_cache",
    "export_public_artifact",
    "gate_public_export",
    "install_redaction",
    "is_armed",
    "is_sensitive_extension",
    "redact_text",
    "rule_by_id",
    "run_audit",
    "safe_repr",
    "sniff_binary",
]
