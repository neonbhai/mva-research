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
