"""Structural tests that enforce GP-01..GP-03 (see docs/golden-principles.md).

These are custom lints. When one fails it prints the offending import *and the
remediation*, because the reader may be an agent whose only view of the rule is
this error message.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[2] / "src" / "mva"

#: GP-01 layering. Lower numbers may not import higher numbers.
LAYERS: dict[str, int] = {
    "models": 0,
    # Foundation utilities: pure, depend on models at most. They sit BELOW config
    # because config needs them (hashing for config_hash, typed errors to raise).
    "errors": 1,
    "determinism": 1,
    "clock": 1,
    # Allele canonicalisation (minimal representation + left-alignment). Layer 1
    # because BOTH ingestion and annotation must use the identical rule, and they
    # are peer stages forbidden to import each other (GP-03). Two copies of this
    # rule is what let equivalent variants fail to join — which looks exactly
    # like "novel and ultra-rare" (ADR 0018).
    "alleles": 1,
    "config": 2,
    # Out-of-repo reference releases: root resolution, the read side of
    # knowledge/manifests/resources.yaml, and the tiered integrity check (ADR 0020).
    # Layer 2 alongside config because it consumes CaseConfig and is consumed by the
    # composition root; it imports no stage and no network client.
    "resources": 2,
    "privacy": 3,
    "ingestion": 4,
    "annotation": 4,
    "phenotype": 4,
    "prioritization": 5,
    "mechanisms": 5,
    "interventions": 5,
    "evidence": 6,
    "reporting": 7,
    "pipeline": 8,
    # The composition root. It was absent, and because both layer tests skip an
    # unmapped owner it was silently exempt from GP-01/GP-03 entirely — the one
    # module that imports every stage was the one nothing checked.
    "orchestrator": 8,
    "cli": 8,
}

REMEDIATION = (
    "\n\nGP-01 remediation: dependencies point one way only "
    "(models -> config -> utils -> privacy -> data-in -> analysis -> evidence -> "
    "reporting -> cli). To fix, either (a) move the shared type down into "
    "src/mva/models/, (b) invert the dependency by passing the needed value in as "
    "a parameter from the composition root (src/mva/pipeline.py), or (c) define a "
    "Protocol in the lower layer that the higher layer implements. Do NOT add the "
    "import and do NOT widen the layer map without a decision record."
)


def _top_level_package(module_path: Path) -> str:
    """The layer name a source file belongs to."""
    rel = module_path.relative_to(SRC)
    return rel.parts[0] if len(rel.parts) > 1 else rel.stem


def _iter_source_files() -> list[Path]:
    return sorted(p for p in SRC.rglob("*.py") if "__pycache__" not in p.parts)


def _imported_mva_modules(tree: ast.AST) -> list[tuple[str, int]]:
    """Every `mva.<layer>` referenced by an import, with its line number."""
    found: list[tuple[str, int]] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.startswith("mva."):
                    found.append((alias.name.split(".")[1], node.lineno))
        elif isinstance(node, ast.ImportFrom):
            if node.level:  # relative import
                continue
            if node.module and node.module.startswith("mva."):
                found.append((node.module.split(".")[1], node.lineno))
    return found


@pytest.mark.unit
def test_no_backward_layer_imports() -> None:
    """GP-01: a module never imports from a higher layer."""
    violations: list[str] = []
    for path in _iter_source_files():
        owner = _top_level_package(path)
        if owner not in LAYERS:
            continue
        own_layer = LAYERS[owner]
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for imported, lineno in _imported_mva_modules(tree):
            if imported not in LAYERS or imported == owner:
                continue
            if LAYERS[imported] > own_layer:
                violations.append(
                    f"  {path.relative_to(SRC.parent.parent)}:{lineno} — "
                    f"'{owner}' (layer {own_layer}) imports "
                    f"'{imported}' (layer {LAYERS[imported]})"
                )
    assert not violations, (
        "GP-01 violated: backward layer imports found.\n" + "\n".join(violations) + REMEDIATION
    )


@pytest.mark.unit
def test_stages_do_not_import_each_other() -> None:
    """GP-03: analysis stages compose only at the composition root."""
    peers = {
        "ingestion",
        "annotation",
        "phenotype",
        "prioritization",
        "mechanisms",
        "interventions",
    }
    violations: list[str] = []
    for path in _iter_source_files():
        owner = _top_level_package(path)
        if owner not in peers:
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for imported, lineno in _imported_mva_modules(tree):
            if imported in peers and imported != owner:
                violations.append(
                    f"  {path.relative_to(SRC.parent.parent)}:{lineno} — "
                    f"stage '{owner}' imports peer stage '{imported}'"
                )
    assert not violations, (
        "GP-03 violated: stages must not import each other.\n"
        + "\n".join(violations)
        + "\n\nGP-03 remediation: wire the two stages together in src/mva/pipeline.py "
        "and pass the result of one into the other as an argument. If they share a "
        "helper, move the helper into src/mva/models/ or a layer-appropriate module."
    )


@pytest.mark.unit
def test_models_are_leaf_modules() -> None:
    """GP-01: `models` is the foundation and imports nothing else from mva."""
    violations: list[str] = []
    for path in sorted((SRC / "models").rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for imported, lineno in _imported_mva_modules(tree):
            if imported != "models":
                violations.append(f"  {path.name}:{lineno} — models imports '{imported}'")
    assert not violations, (
        "GP-01 violated: mva.models must be a leaf package.\n"
        + "\n".join(violations)
        + "\n\nGP-01 remediation: models is the shared vocabulary; it cannot depend on "
        "anything that depends on it. Move the needed logic into the caller."
    )


# ---------------------------------------------------------------------------
# PRIV-05: what may be imported on the patient-data path.
#
# Read the docstring of `test_no_network_clients_in_sensitive_stages` before
# relying on any of this. It is an import lint. It is not a network boundary.
# ---------------------------------------------------------------------------

#: Packages that speak a network protocol on the caller's behalf. The realistic
#: accident: someone adds `requests.get(...)` to an annotation step.
NETWORK_CLIENT_IMPORTS: frozenset[str] = frozenset(
    {
        "requests",
        "httpx",
        "urllib",
        "urllib3",
        "aiohttp",
        "http",
        "ftplib",
        "smtplib",
        "poplib",
        "imaplib",
        "nntplib",
        "telnetlib",
        "xmlrpc",
        "webbrowser",
        "paramiko",
        "pycurl",
        "websockets",
        "grpc",
        "boto3",
    }
)

#: The layer underneath those. `socket` is what every client above is built on,
#: and `socket.send` on an already-connected socket emits no audit event at all,
#: so `mva.privacy.netguard` cannot see it (see that module's honest-limits list).
#: `asyncio` carries `open_connection` and `loop.sock_connect`.
DIRECT_SOCKET_IMPORTS: frozenset[str] = frozenset({"socket", "socketserver", "ssl", "asyncio"})

#: Routes that reach the network *without* the runtime guard ever being consulted.
#: These are not network clients; they are the documented ways around the audit
#: hook, and `mva.privacy.netguard` names all three in its own docstring:
#:
#: * `ctypes`/`cffi` — `CDLL("libc").connect(...)` never touches the socket module.
#: * `subprocess` — audit hooks are per-interpreter. A spawned child has an
#:   unrestricted network, and no amount of care in this process changes that.
#: * `importlib` — `import_module("socket")` is invisible to an AST import lint,
#:   which is to say invisible to this test.
AUDIT_HOOK_BYPASS_IMPORTS: frozenset[str] = frozenset({"ctypes", "cffi", "subprocess", "importlib"})

FORBIDDEN_ON_PATIENT_PATH: frozenset[str] = (
    NETWORK_CLIENT_IMPORTS | DIRECT_SOCKET_IMPORTS | AUDIT_HOOK_BYPASS_IMPORTS
)

SENSITIVE_PACKAGES: frozenset[str] = frozenset(
    {"ingestion", "annotation", "phenotype", "prioritization"}
)

#: `(repo-relative path, import root) -> why this one is allowed`.
#:
#: Every entry is a decision that someone has to defend, which is the point: an
#: exemption list makes a hole visible, whereas a short forbidden list makes the
#: same hole invisible by simply never asking about it. Widening this map is a
#: reviewable event; it must never be widened to make a red test green.
PATIENT_PATH_IMPORT_EXEMPTIONS: dict[tuple[str, str], str] = {
    (
        "annotation/snpeff_local.py",
        "subprocess",
    ): (
        "SnpEff is a Java program; there is no in-process alternative. The child "
        "gets the proband's coordinates as a VCF on stdin (never a file, never an "
        "argv), a hand-built 6-key environment rather than os.environ, a scratch "
        "cwd/HOME/TMPDIR, and -nodownload/-noStats/-noLog. Those constrain SnpEff's "
        "EXPECTED behaviour; none of them is a network boundary for the child, and "
        "this test cannot see inside it at all. See TD-06."
    ),
    (
        "ingestion/reader.py",
        "importlib",
    ): (
        "`importlib.util.find_spec('pysam')` — a capability probe that imports "
        "nothing. It is exempted rather than ignored because `import_module` in the "
        "same package would defeat every import assertion in this file."
    ),
}


def forbidden_imports(source: str, *, filename: str = "<test>") -> set[tuple[str, int]]:
    """Import roots this lint objects to, as `(root, lineno)`.

    Extracted so the matcher can be exercised against sources written for the
    purpose. A lint asserted only against the tree it already passes on is a lint
    whose blind spots are structurally invisible — which is exactly how `socket`
    and `ctypes` went unlisted while the docstring claimed a structural guarantee.
    """
    tree = ast.parse(source, filename=filename)
    found: set[tuple[str, int]] = set()
    for node in ast.walk(tree):
        names: list[str] = []
        if isinstance(node, ast.Import):
            names = [alias.name.split(".")[0] for alias in node.names]
        elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
            names = [node.module.split(".")[0]]
        else:
            continue
        found.update((name, node.lineno) for name in names if name in FORBIDDEN_ON_PATIENT_PATH)
    return found


@pytest.mark.unit
def test_no_network_clients_in_sensitive_stages() -> None:
    """PRIV-05, and an honest statement of how far it reaches.

    **What this proves.** No module in `ingestion`, `annotation`, `phenotype` or
    `prioritization` contains a top-level or function-level `import` of a network
    client, of `socket`/`ssl`/`asyncio`, or of `ctypes`/`cffi`/`subprocess`/
    `importlib` — except for the entries in
    :data:`PATIENT_PATH_IMPORT_EXEMPTIONS`, each of which states its reason. That
    is a real property and it catches the realistic accident: someone adding
    `requests.get(...)` to an annotation step, or reaching for a raw socket.

    **What this does NOT prove, and must not be read as proving.**

    * **It says nothing about child processes.** `annotation/snpeff_local.py` is
      exempted above and spawns a JVM with the proband's coordinates on its stdin.
      Audit hooks are per-interpreter; that child has a completely unrestricted
      network. `-nodownload/-noStats/-noLog` constrain what SnpEff is *expected*
      to do — they are flags to a cooperative program, not a boundary around an
      uncooperative one. No test in this repository observes the child's sockets.
    * **It says nothing about C extensions.** `pysam`/htslib and `cyvcf2` are
      imported by five modules on this path and call `connect(2)` and libcurl from
      C, emitting no Python audit event. They are deliberately *not* forbidden
      here: they are the primary genomics I/O path and there is no alternative.
      htslib's CRAM reference auto-fetch is handled by an environment control
      (`REF_PATH=/dev/null`), not by this lint and not by `netguard`.
    * **It is an AST lint over import statements only.** `getattr(__builtins__,
      '__import__')('socket')`, an `exec` of a computed string, or a socket handed
      in as a constructor argument all pass it. Adding `importlib` above raises the
      cost of the obvious evasion; it does not close the class.
    * **It is a lint, not a runtime control.** It observes source text at test
      time. It cannot observe a run.

    The only control that actually bounds the child process is at the OS. On the
    macOS target that is `sandbox-exec` with a Seatbelt profile denying
    `network-outbound`, wrapped around the *whole* `mva run` invocation so the
    JVM inherits it — a child cannot escape its parent's sandbox. Nothing in this
    repository applies or verifies it; `NetworkProfile.OFFLINE_ENFORCED` is the
    operator asserting they did, which the CLI prints as an assertion rather than
    an observation. That gap is TD-06, and it is the reason this docstring is long.
    """
    violations: list[str] = []
    for path in _iter_source_files():
        if _top_level_package(path) not in SENSITIVE_PACKAGES:
            continue
        relative = path.relative_to(SRC).as_posix()
        for name, lineno in sorted(forbidden_imports(path.read_text(encoding="utf-8"))):
            if (relative, name) in PATIENT_PATH_IMPORT_EXEMPTIONS:
                continue
            violations.append(
                f"  {path.relative_to(SRC.parent.parent)}:{lineno} — imports '{name}' "
                f"on the patient-data path"
            )
    assert not violations, (
        "PRIV-05 violated: a network-reaching import on the patient-data path.\n"
        + "\n".join(sorted(violations))
        + "\n\nPRIV-05 remediation: annotation of patient coordinates must be LOCAL. "
        "Use a pre-downloaded, hash-pinned resource under knowledge/ via an adapter. "
        "If a remote source is genuinely needed, it belongs in a separate offline "
        "acquisition step that fetches PUBLIC reference data only and never sends "
        "proband coordinates. If the import is genuinely unavoidable — a local tool "
        "with no in-process equivalent — add it to PATIENT_PATH_IMPORT_EXEMPTIONS in "
        "tests/unit/test_architecture.py WITH the reason and the residual risk, and "
        "update docs/privacy-model.md. Do not widen the forbidden set's complement."
    )


@pytest.mark.unit
def test_the_patient_path_lint_catches_the_direct_and_bypass_routes() -> None:
    """The lint is exercised against sources written to defeat it.

    `socket`, `ctypes` and `importlib` were absent from the forbidden set while
    the test's own docstring claimed the rule was "structural, not conventional".
    Asserting a lint only against a tree it already passes on cannot detect that:
    the missing cases produce no output to be wrong about. So the matcher is fed
    the code it is supposed to object to.
    """
    must_catch = {
        "import socket": "socket",
        "import ctypes": "ctypes",
        "from ctypes import CDLL": "ctypes",
        "import importlib": "importlib",
        "import subprocess": "subprocess",
        "import ssl": "ssl",
        "import asyncio": "asyncio",
        "import requests": "requests",
        "import urllib3": "urllib3",
        "from http import client": "http",
        "def f():\n    import socket\n": "socket",
    }
    missed = [
        source
        for source, expected in must_catch.items()
        if expected not in {name for name, _ in forbidden_imports(source)}
    ]
    assert not missed, (
        "the PRIV-05 import lint does not object to:\n"
        + "\n".join(f"  {source!r}" for source in missed)
        + "\n\nAdd the missing root to NETWORK_CLIENT_IMPORTS, DIRECT_SOCKET_IMPORTS "
        "or AUDIT_HOOK_BYPASS_IMPORTS in this file."
    )

    # And must not object to the genomics backends, which are the point of the
    # stage. Over-forbidding would force the exemption list to swallow the rule.
    for allowed in ("import pysam", "import cyvcf2", "import gzip", "import struct"):
        assert not forbidden_imports(allowed), f"{allowed!r} is not a network route"


@pytest.mark.unit
def test_every_patient_path_import_exemption_is_real_and_reasoned() -> None:
    """An exemption for an import that no longer exists is a rule quietly relaxed.

    Both halves matter. A stale entry means the forbidden set is weaker than it
    reads. A reason-less entry means nobody can tell whether the hole was decided
    or inherited.
    """
    stale: list[str] = []
    for (relative, name), reason in sorted(PATIENT_PATH_IMPORT_EXEMPTIONS.items()):
        path = SRC / relative
        assert path.is_file(), f"exemption names a file that does not exist: {relative}"
        assert len(reason) > 80, (
            f"exemption ({relative}, {name}) has no substantive reason; an exemption "
            "without a stated residual risk is indistinguishable from an oversight"
        )
        if name not in {found for found, _ in forbidden_imports(path.read_text(encoding="utf-8"))}:
            stale.append(f"  ({relative}, {name}) — no longer imported")
    assert not stale, (
        "PATIENT_PATH_IMPORT_EXEMPTIONS has entries that no longer apply:\n"
        + "\n".join(stale)
        + "\n\nDelete them. An exemption that outlives its import silently re-opens "
        "the rule for the next person who adds that import back."
    )


@pytest.mark.unit
def test_the_only_subprocess_on_the_patient_data_path_is_the_declared_one() -> None:
    """The subprocess hole is bounded to one file, and that is all it is.

    `mva.privacy.netguard` states plainly that a spawned child has an unrestricted
    network, and the pipeline runs the offline profile non-strict by default
    because `strict=True` would block `subprocess.Popen` and so break the
    provenance manifest's `git` calls. So the guard cannot stop the SnpEff JVM, and
    is not configured to try.

    What is enforceable from here is that the exposure does not grow silently: a
    second stage shelling out on this path fails this test rather than inheriting
    the first one's justification. It is a containment assertion about the source
    tree, and explicitly NOT a claim that the one remaining child is contained.
    """
    spawning = {
        relative for (relative, name) in PATIENT_PATH_IMPORT_EXEMPTIONS if name == "subprocess"
    }
    assert spawning == {"annotation/snpeff_local.py"}, (
        "the set of patient-path modules permitted to spawn a process changed.\n"
        f"  now: {sorted(spawning)}\n"
        "\nEach one is an unpoliced network peer holding proband coordinates. "
        "Adding another needs a decision record and a docs/privacy-model.md update, "
        "not an exemption-list edit."
    )


@pytest.mark.unit
def test_no_naive_datetime_now() -> None:
    """GP-30: wall-clock reads break determinism and must go through the clock module."""
    violations: list[str] = []
    for path in _iter_source_files():
        if path.name == "clock.py":
            continue
        text = path.read_text(encoding="utf-8")
        for lineno, line in enumerate(text.splitlines(), start=1):
            stripped = line.strip()
            if stripped.startswith("#"):
                continue
            if "datetime.now(" in line or "datetime.utcnow(" in line or "time.time(" in line:
                violations.append(f"  {path.relative_to(SRC.parent.parent)}:{lineno}")
    assert not violations, (
        "GP-30 violated: direct wall-clock access.\n"
        + "\n".join(violations)
        + "\n\nGP-30 remediation: inject a Clock (src/mva/clock.py) and call "
        "`clock.now()`. Determinism tests compare two runs byte-for-byte; an "
        "embedded wall-clock timestamp makes that impossible."
    )


# ---------------------------------------------------------------------------
# The allele-balance lint. A review finding promoted into a rule (CLAUDE.md).
# ---------------------------------------------------------------------------

#: Files entitled to read the raw ``Genotype.allele_balance``.
ALLELE_BALANCE_OWNERS: frozenset[str] = frozenset(
    {
        "models/variant.py",  # defines it, and defines the site-aware fraction over it
        "ingestion/qc.py",  # reports it alongside the fraction in its evidence payload
        # The two below RENDER the raw measured field verbatim into an artifact
        # column; neither applies the heterozygous band to it, so neither can
        # produce the misjudgement this lint exists to prevent. They are still
        # showing a reader 0.913 where the site fraction is 0.477, which is
        # tracked as TD-16 — the fix is to render both numbers, and it belongs to
        # whoever owns those packages.
        "evidence/store.py",
        "reporting/dossier.py",
    }
)

_ALLELE_BALANCE_READ = "allele_balance"

ALLELE_BALANCE_REMEDIATION = (
    "\n\nRemediation: read `VariantRecord.allele_fraction`, not "
    "`genotype.allele_balance`. On a record decomposed from a multiallelic site "
    "the raw balance is alt/(ref+alt) over the SITE's reference depth and "
    "excludes the reads carrying the other ALT: a textbook compound heterozygote "
    "at AD=2,21,21 reads 0.913 instead of 0.477, gets flagged `low_quality_call` "
    "and `possible_mosaic`, and loses ~0.18 of composite to the shape of its VCF "
    "line. `mva.ingestion.qc` was written to fix exactly this, and two "
    "prioritisation stages then undid it by re-deriving the wrong quantity from "
    "the raw attribute. The lint exists because the fix is invisible at the call "
    "site — both spellings compile, both return a float, and only one is right. "
    "If a new module genuinely needs the raw per-allele balance, add it to "
    "ALLELE_BALANCE_OWNERS in this file with a comment saying why "
    "(ASSUMPTION-MOSAIC-02)."
)


@pytest.mark.unit
def test_allele_balance_is_read_only_where_it_is_defined() -> None:
    """The site-aware allele fraction is the only quantity the het band is applied to.

    Structural rather than conventional: an attribute read is easy to write by
    accident and impossible to spot in review, which is how this regressed once
    already. Every exemption in ``ALLELE_BALANCE_OWNERS`` carries the reason it
    is there, so widening the set is a visible decision rather than a diff.
    """
    violations: list[str] = []
    for path in _iter_source_files():
        relative = path.relative_to(SRC).as_posix()
        if relative in ALLELE_BALANCE_OWNERS:
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Attribute) and node.attr == _ALLELE_BALANCE_READ:
                violations.append(
                    f"  {path.relative_to(SRC.parent.parent)}:{node.lineno} — "
                    f"reads `.{_ALLELE_BALANCE_READ}`"
                )

    assert not violations, (
        "Raw allele balance read outside its owning modules.\n"
        + "\n".join(violations)
        + ALLELE_BALANCE_REMEDIATION
    )


@pytest.mark.unit
def test_every_source_module_is_in_the_layer_map() -> None:
    """An unmapped module must FAIL, not be silently skipped.

    Both layer tests do ``if owner not in LAYERS: continue``. That is a sensible
    guard for ``__init__``, but it meant the composition root — the one module
    that imports every stage — was exempt from GP-01 and GP-03 entirely, simply
    because nobody had added it. A rule that silently excuses whatever it does
    not recognise is not enforced; it is decorative.

    Remediation when this fails: put the new module in LAYERS at the layer it
    genuinely belongs to. Do NOT widen an existing layer to accommodate an import
    that should not exist — that inverts the check into a record of the code's
    current shape instead of a constraint on it.
    """
    owners = {_top_level_package(path) for path in _iter_source_files()}
    owners.discard("__init__")

    unmapped = sorted(owners - set(LAYERS))
    assert not unmapped, (
        "Source modules missing from the enforced layer map:\n"
        + "\n".join(f"    src/mva/{name}" for name in unmapped)
        + "\n\nThey are currently exempt from GP-01/GP-03. Add each to LAYERS."
    )

    phantom = sorted(set(LAYERS) - owners)
    assert not phantom, (
        "LAYERS declares owners with no matching source:\n"
        + "\n".join(f"    {name}" for name in phantom)
        + "\n\nA phantom entry makes the map look more complete than it is. "
        "Remove it, or add the module it was reserving a place for."
    )
