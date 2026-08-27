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
    "config": 2,
    "privacy": 3,
    "ingestion": 4,
    "annotation": 4,
    "phenotype": 4,
    "knowledge": 4,
    "prioritization": 5,
    "mechanisms": 5,
    "interventions": 5,
    "evidence": 6,
    "reporting": 7,
    "pipeline": 8,
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


@pytest.mark.unit
def test_no_network_clients_in_sensitive_stages() -> None:
    """PRIV-05: modules on the patient-data path must not import network clients.

    A remote annotation service that receives a proband's coordinates is a
    re-identification vector. The rule is structural, not conventional: the import
    simply may not exist in these packages.
    """
    forbidden = {"requests", "httpx", "urllib", "aiohttp", "http", "ftplib", "smtplib"}
    sensitive = {"ingestion", "annotation", "phenotype", "prioritization"}
    violations: list[str] = []
    for path in _iter_source_files():
        if _top_level_package(path) not in sensitive:
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            names: list[str] = []
            if isinstance(node, ast.Import):
                names = [a.name.split(".")[0] for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
                names = [node.module.split(".")[0]]
            else:
                continue
            for name in names:
                if name in forbidden:
                    violations.append(
                        f"  {path.relative_to(SRC.parent.parent)}:{node.lineno} — "
                        f"imports network client '{name}'"
                    )
    assert not violations, (
        "PRIV-05 violated: network client on the patient-data path.\n"
        + "\n".join(violations)
        + "\n\nPRIV-05 remediation: annotation of patient coordinates must be LOCAL. "
        "Use a pre-downloaded, hash-pinned resource under knowledge/ via an adapter. "
        "If a remote source is genuinely needed, it belongs in a separate offline "
        "acquisition step that fetches PUBLIC reference data only and never sends "
        "proband coordinates."
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
