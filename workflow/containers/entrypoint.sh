#!/bin/sh
# =============================================================================
# Container entrypoint for mva-research.
# =============================================================================
#
# Two modes, and the default one is unchanged from a bare `ENTRYPOINT ["mva"]`:
#
#   docker run --rm -v "$PWD/out:/workspace" IMAGE demo
#       Seed the mounted workspace with the synthetic fixtures and run the
#       synthetic case end to end. This is the container equivalent of
#       `just demo`.
#
#   docker run --rm ... IMAGE run all --config ... --workspace /workspace
#   docker run --rm ... IMAGE verify determinism ...
#   docker run --rm ... IMAGE                      (-> mva --help)
#       Anything that is not `demo` is passed straight through to `mva`, so
#       every invocation documented in README.md keeps working.
#
# DIVERGENCE IS A BUILD FAILURE. The `mva run all` arguments and the fixture
# list below are asserted against the justfile's `demo` and `_seed-workspace`
# recipes at image build time by assert-demo-parity.sh. Change one without the
# other and the image does not build.
#
# WHY THIS DOES NOT CLEAN FIRST. `just demo` depends on `clean-demo`, which
# `rm -rf`s the demo workspace. This script does not: /workspace is the caller's
# bind mount, and a container that deletes the directory it was handed is a
# container nobody should point at real data. Re-seeding is idempotent instead.
# =============================================================================

set -eu

: "${MVA_WORKSPACE:=/workspace}"

# PRIV-08: temp files belong inside the workspace, never the host's /tmp. The
# image sets TMPDIR, but Python's tempfile needs the directory to exist, and on
# a fresh bind mount it does not.
mkdir -p "${TMPDIR:-$MVA_WORKSPACE/tmp}" 2>/dev/null || true

# Copy the synthetic fixtures into the workspace, exactly as the justfile's
# `_seed-workspace` recipe does.
#
# `cp -f` rather than `cp`: /app is mode a-w, so the fixtures are 0444 and so is
# any copy of them. Without -f, a second `demo` run into the same workspace
# fails on the unwritable destination.
seed_workspace() {
    mkdir -p "$MVA_WORKSPACE/inputs"
    cp -f /app/tests/fixtures/synthetic/synthetic_case.vcf "$MVA_WORKSPACE/inputs/"
    cp -f /app/tests/fixtures/synthetic/synthetic_phenotype.tsv "$MVA_WORKSPACE/inputs/"
}

case "${1-}" in
    demo)
        shift
        seed_workspace
        exec mva run all \
            --config config/synthetic-case.yaml \
            --defaults config/default.yaml \
            --workspace "$MVA_WORKSPACE" \
            "$@"
        ;;
    *)
        exec mva "$@"
        ;;
esac
