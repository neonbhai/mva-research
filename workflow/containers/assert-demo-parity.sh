#!/usr/bin/env bash
# =============================================================================
# Build-time assertion: `docker run IMAGE demo` == `just demo`.
# =============================================================================
#
# The container cannot simply call `just demo`. That recipe runs `uv run`, which
# wants to sync a project environment, and it writes to $TMPDIR/mva-research-demo
# rather than to the mounted workspace. So the entrypoint restates the run.
#
# A restatement is a copy, and a copy drifts. This script is the thing that stops
# it: it extracts the `mva run all` argument list and the seeded fixture list
# from BOTH the justfile and the entrypoint, and fails the image build if they
# disagree. The two cannot diverge silently, because the divergence is a red
# build rather than a demo that quietly runs a different case.
#
# Usage: assert-demo-parity.sh <justfile> <entrypoint.sh>
# =============================================================================

set -euo pipefail

justfile="${1:?usage: assert-demo-parity.sh <justfile> <entrypoint.sh>}"
entrypoint="${2:?usage: assert-demo-parity.sh <justfile> <entrypoint.sh>}"

for f in "$justfile" "$entrypoint"; do
    [ -r "$f" ] || { echo "assert-demo-parity: cannot read $f" >&2; exit 1; }
done

fail() {
    echo "" >&2
    echo "=============================================================" >&2
    echo " DEMO PARITY ASSERTION FAILED" >&2
    echo "=============================================================" >&2
    printf '%s\n' "$@" >&2
    echo "" >&2
    echo " Remediation: make these two agree, then rebuild." >&2
    echo "   native    justfile, recipes 'demo' and '_seed-workspace'" >&2
    echo "   container workflow/containers/entrypoint.sh" >&2
    echo "" >&2
    echo " They are asserted equal here so that 'docker run IMAGE demo'" >&2
    echo " cannot quietly run a different case than 'just demo'." >&2
    echo "=============================================================" >&2
    exit 1
}

# Join backslash continuations into one logical line, then collapse whitespace.
normalise() {
    awk '{
        line = $0
        if (sub(/\\[[:space:]]*$/, "", line)) { printf "%s ", line; next }
        print line
    }' | tr '\t' ' ' | tr -s ' ' | sed -e 's/^ //' -e 's/ *$//'
}

# One named recipe's body, not the whole justfile.
#
# Scoped deliberately: other recipes invoke the pipeline too, and a new one that
# happens to call `mva run all` must not be mistaken for the demo's own
# invocation. A recipe header sits at column 0 and the name is followed by `:` or
# by its parameters, so `demo` matches `demo:` and not `demo-determinism:`.
recipe_block() {
    awk -v want="$2" '
        /^[^ \t#]/ { inblock = ($0 ~ "^" want "([ \t]|:)") }
        inblock    { print }
    ' "$1"
}

# The `mva run all` invocation, minus the parts that legitimately differ:
# --workspace (a justfile variable vs. the mounted path) and a passthrough "$@".
run_args() {
    normalise \
        | grep -o 'mva run all .*' \
        | sed -e 's/ --workspace.*//' -e 's/ *"\$@" *//' \
        | sed -e 's/ *$//'
}

# Every synthetic fixture path either side seeds. Compared as a sorted set, so
# the leading /app/ in the container paths and the recipe's ordering do not
# matter, but adding, dropping or renaming a fixture does.
fixtures() {
    normalise | grep -o 'tests/fixtures/synthetic/[A-Za-z0-9_.-]*' | sort -u
}

# --- the run invocation ------------------------------------------------------

just_run_count=$(recipe_block "$justfile" demo | run_args | grep -c . || true)
entry_run_count=$(run_args < "$entrypoint" | grep -c . || true)

[ "$just_run_count" = "1" ] || fail \
    " Expected exactly one 'mva run all' invocation in the justfile's 'demo'" \
    " recipe, found ${just_run_count}." \
    " Either the recipe was renamed, or it no longer runs the case itself."
[ "$entry_run_count" = "1" ] || fail \
    " Expected exactly one 'mva run all' invocation in the entrypoint, found ${entry_run_count}."

just_run=$(recipe_block "$justfile" demo | run_args)
entry_run=$(run_args < "$entrypoint")

[ "$just_run" = "$entry_run" ] || fail \
    " The demo run arguments differ." \
    "" \
    "   just demo:  ${just_run}" \
    "   container:  ${entry_run}"

# --- the seeded fixtures ----------------------------------------------------

just_fixtures=$(recipe_block "$justfile" _seed-workspace | fixtures)
entry_fixtures=$(fixtures < "$entrypoint")

[ -n "$just_fixtures" ] || fail \
    " Found no tests/fixtures/synthetic/* paths in the justfile's" \
    " '_seed-workspace' recipe. Has it been renamed or moved?"

[ "$just_fixtures" = "$entry_fixtures" ] || fail \
    " The seeded fixture lists differ." \
    "" \
    "   just _seed-workspace:  $(echo "$just_fixtures" | tr '\n' ' ')" \
    "   container:             $(echo "$entry_fixtures" | tr '\n' ' ')"

echo "demo parity: ok -- '${entry_run} --workspace \$MVA_WORKSPACE'"
echo "demo parity: fixtures $(echo "$entry_fixtures" | tr '\n' ' ')"
