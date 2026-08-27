# mva-research command runner.
# `just verify` is the acceptance gate (ADR 0009). Everything else is a shortcut into it.

set shell := ["bash", "-uc"]

# Workspace for the synthetic demo.
#
# Deliberately OUTSIDE the repository (ADR 0006). The privacy audit's
# `workspace_containment` check fails an in-repo workspace, and it is right to:
# a workspace inside the tree is one `git add -A` from being committed. Running
# the demo the same way a real case must be run keeps the honest path the
# default path. Override with: just demo DEMO_ROOT=/some/other/dir
DEMO_ROOT := env_var_or_default("TMPDIR", "/tmp")
DEMO_WORKSPACE := DEMO_ROOT / "mva-research-demo"
DET_WORKSPACE := DEMO_ROOT / "mva-research-demo-det"

default:
    @just --list

# ---------------------------------------------------------------- setup

# Install the toolchain and sync dependencies.
bootstrap:
    @echo "==> checking uv"
    @command -v uv >/dev/null || { echo "uv not found. Install: https://docs.astral.sh/uv/"; exit 1; }
    @echo "==> installing Python 3.12"
    uv python install 3.12
    @echo "==> syncing dependencies (incl. genomics + workflow extras)"
    uv sync --extra genomics --extra workflow
    @echo "==> installing git hooks"
    @just install-hooks
    @echo "==> bootstrap complete. Next: just verify"

# Install the pre-commit privacy gate.
install-hooks:
    @mkdir -p .githooks
    @printf '#!/usr/bin/env bash\n# Blocking privacy gate (PRIV-01). Bypass only with --no-verify, and never for real data.\nexec uv run mva privacy audit --staged --strict\n' > .githooks/pre-commit
    @chmod +x .githooks/pre-commit
    @git config core.hooksPath .githooks
    @echo "git hooks installed (core.hooksPath=.githooks)"

# ---------------------------------------------------------------- quality

# Lint and format check.
lint:
    uv run ruff check src tests
    uv run ruff format --check src tests

# Apply formatting and safe fixes.
fmt:
    uv run ruff check --fix src tests
    uv run ruff format src tests

# Static type check (strict).
typecheck:
    uv run pyright src tests

# Full test suite.
test:
    uv run pytest -q

# Fast unit tests only.
test-unit:
    uv run pytest -q -m unit

# Golden expectation locks.
test-golden:
    uv run pytest -q -m golden

# Structural lints: layering, docs integrity, no-network-on-sensitive-path.
arch:
    uv run pytest -q tests/unit/test_architecture.py tests/unit/test_docs_integrity.py

# THE GATE. Everything must pass; nothing here is advisory (ADR 0009).
verify: lint typecheck arch test privacy-audit
    @echo ""
    @echo "  verify: all gates passed"

# ---------------------------------------------------------------- privacy

# Scan the repository for privacy boundary violations.
privacy-audit:
    uv run mva privacy audit --repo .

# Scan only what is staged for commit.
privacy-audit-staged:
    uv run mva privacy audit --repo . --staged --strict

# ---------------------------------------------------------------- demo

# Run the synthetic case end to end and write all artifacts.
demo: clean-demo
    @just _seed-workspace "{{DEMO_WORKSPACE}}"
    uv run mva run all \
        --config config/synthetic-case.yaml \
        --defaults config/default.yaml \
        --workspace "{{DEMO_WORKSPACE}}"
    @echo ""
    @echo "  demo artifacts: {{DEMO_WORKSPACE}}/runs/"

# Copy the synthetic inputs into a workspace. Synthetic fixtures only.
_seed-workspace WS:
    @mkdir -p "{{WS}}/inputs"
    @cp tests/fixtures/synthetic/synthetic_case.vcf "{{WS}}/inputs/"
    @cp tests/fixtures/synthetic/synthetic_phenotype.tsv "{{WS}}/inputs/"

# Show what the demo produced.
demo-artifacts:
    @find "{{DEMO_WORKSPACE}}/runs" -type f 2>/dev/null | sort || echo "no demo run found; try: just demo"

# Show the privacy audit from the most recent demo run.
demo-audit:
    @cat "{{DEMO_WORKSPACE}}"/runs/*/privacy/privacy_audit.md 2>/dev/null || echo "no demo run found; try: just demo"

# Prove determinism: run the demo twice and compare artifact hashes.
demo-determinism:
    @just _seed-workspace "{{DET_WORKSPACE}}"
    uv run mva verify determinism \
        --config config/synthetic-case.yaml \
        --defaults config/default.yaml \
        --workspace "{{DET_WORKSPACE}}"

# Remove the demo workspace.
clean-demo:
    @rm -rf "{{DEMO_WORKSPACE}}" "{{DET_WORKSPACE}}"
    @echo "demo workspace removed: {{DEMO_WORKSPACE}}"

# ---------------------------------------------------------------- workflow

# Run the Snakemake DAG for the synthetic case.
#
# `--configfile` is placed LAST on purpose: it takes nargs="+", so anything
# following it is greedily swallowed. With it last, a bare target parses
# correctly. Usage: `just snakemake track1`, `just snakemake -n`.
snakemake *ARGS:
    MVA_WORKSPACE="{{DEMO_WORKSPACE}}" uv run snakemake --cores 1 {{ARGS}} \
        --configfile config/synthetic-case.yaml

# Render the workflow DAG.
dag:
    @MVA_WORKSPACE="{{DEMO_WORKSPACE}}" uv run snakemake --dag \
        --configfile config/synthetic-case.yaml | head -50

# ---------------------------------------------------------------- housekeeping

clean: clean-demo
    @rm -rf .pytest_cache .ruff_cache .coverage htmlcov dist build
    @find . -name '__pycache__' -type d -prune -exec rm -rf {} + 2>/dev/null || true
    @echo "clean"

# Print the environment the pipeline sees.
info:
    @uv run mva info
