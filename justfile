# mva-research command runner.
# `just verify` is the acceptance gate (ADR 0009). Everything else is a shortcut into it.

set shell := ["bash", "-uc"]

# Workspace for the synthetic demo. Overridable: `just demo WORKSPACE=/path`.
DEMO_WORKSPACE := justfile_directory() + "/.demo-workspace"

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
    @mkdir -p "{{DEMO_WORKSPACE}}/inputs"
    @cp tests/fixtures/synthetic/synthetic_case.vcf "{{DEMO_WORKSPACE}}/inputs/"
    @cp tests/fixtures/synthetic/synthetic_phenotype.tsv "{{DEMO_WORKSPACE}}/inputs/"
    uv run mva run all \
        --config config/synthetic-case.yaml \
        --defaults config/default.yaml \
        --workspace "{{DEMO_WORKSPACE}}" \
        --allow-workspace-in-repo
    @echo ""
    @echo "  demo artifacts: {{DEMO_WORKSPACE}}/runs/"

# Show what the demo produced.
demo-artifacts:
    @find "{{DEMO_WORKSPACE}}/runs" -type f 2>/dev/null | sort || echo "no demo run found; try: just demo"

# Prove determinism: run the demo twice and compare artifact hashes.
demo-determinism:
    uv run mva verify determinism \
        --config config/synthetic-case.yaml \
        --defaults config/default.yaml \
        --workspace "{{DEMO_WORKSPACE}}-det" \
        --allow-workspace-in-repo

# Remove the demo workspace.
clean-demo:
    @rm -rf "{{DEMO_WORKSPACE}}" "{{DEMO_WORKSPACE}}-det"
    @echo "demo workspace removed"

# ---------------------------------------------------------------- workflow

# Run the Snakemake DAG for the synthetic case.
snakemake *ARGS:
    uv run snakemake --cores 1 --configfile config/synthetic-case.yaml {{ARGS}}

# Render the workflow DAG.
dag:
    uv run snakemake --dag --configfile config/synthetic-case.yaml | head -50

# ---------------------------------------------------------------- housekeeping

clean: clean-demo
    @rm -rf .pytest_cache .ruff_cache .coverage htmlcov dist build
    @find . -name '__pycache__' -type d -prune -exec rm -rf {} + 2>/dev/null || true
    @echo "clean"

# Print the environment the pipeline sees.
info:
    @uv run mva info
