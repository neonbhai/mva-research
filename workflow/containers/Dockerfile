# syntax=docker/dockerfile:1.7
# =============================================================================
# mva-research — reproducible linux/amd64 image
# =============================================================================
#
# WHAT THIS IS FOR: CI and third-party reproduction on linux. It is NOT the local
# development environment on macOS/Apple Silicon — see README.md in this
# directory for why, and for the bind-mount rules that apply at runtime.
#
# NO PATIENT DATA IS EVER COPIED IN. There is no `COPY` of anything outside the
# repository, no `ADD` from a URL, and no `VOLUME` pre-populated with data. An
# image containing patient data is a *distributable copy* of that data: it is
# pushed to registries, cached in layers that survive deletion of the file, and
# shared by digest. Inputs arrive at `docker run` time as a read-only bind mount
# of the external workspace (ADR 0006, GP-40).
#
# BASE IMAGE PIN. `python:3.12-slim` is a moving tag, so the FROM line below
# pins the linux/amd64 manifest by digest.
#
# A digest that was never verified against a registry is worse than none,
# because it looks pinned. This one was not pasted from memory: it was read back
# from Docker Hub on the date recorded below, with
#
#   docker buildx imagetools inspect python:3.12-slim \
#     --format '{{ range .Manifest.Manifests }}{{ if eq .Platform.Architecture "amd64" }}{{ .Digest }}{{ end }}{{ end }}'
#
# and then confirmed to be the amd64 manifest and not a sibling, by inspecting
# the digest itself and checking what it reports:
#
#   docker buildx imagetools inspect \
#     python:3.12-slim@sha256:e8a38829209aea82774882e5157e1f7c3e8c8393ff53a697b3580896e57f345c \
#     --format '{{ .Image.OS }}/{{ .Image.Architecture }}'      # -> linux/amd64
#
#   # pinned digest (linux/amd64, verified against docker.io on date: 2026-10-25):
#   # sha256:e8a38829209aea82774882e5157e1f7c3e8c8393ff53a697b3580896e57f345c
#   # upstream at that digest: python 3.12 / slim-trixie,
#   # docker-library/python commit 2a3b794c223ab067d122719541cdd54a068732a5
#
# NOTE ON THE COMMAND ABOVE: the Go template field is `.Platform.Architecture`,
# capitalised. An earlier revision of this comment used `.Platform.architecture`,
# which buildx v0.34 rejects outright --
# `can't evaluate field architecture in type *v1.Platform`. Left corrected rather
# than as prose, because the next person to re-verify will paste it verbatim.
#
# RE-VERIFYING. Upstream retags `3.12-slim` on every patch release and security
# rebuild, so this digest goes stale rather than wrong: the pin keeps building
# the same image, it just stops being the newest one. To move it, run the command
# again, update both the FROM line and the date, and rebuild.
# =============================================================================

FROM --platform=linux/amd64 python:3.12-slim@sha256:e8a38829209aea82774882e5157e1f7c3e8c8393ff53a697b3580896e57f345c AS base

# uv is pinned by version AND copied from its own published image, so the build
# never curls an install script over the network at an unpinned URL.
COPY --from=ghcr.io/astral-sh/uv:0.9.7 /uv /uvx /usr/local/bin/

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONHASHSEED=0 \
    # Byte-identical repeat runs are an acceptance criterion (GP-30). A salted
    # hash seed would make set and dict iteration order vary per container start.
    TZ=UTC \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    # PRIV-06: htslib fetches missing CRAM reference sequences from
    # www.ebi.ac.uk by default. Deny it at the image level, not per invocation.
    REF_PATH=/dev/null \
    UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PROJECT_ENVIRONMENT=/opt/venv \
    PATH=/opt/venv/bin:$PATH

WORKDIR /app

# ---------------------------------------------------------------------------
# Dependencies, from the lockfile only. This layer is cached on the lockfile, so
# a source edit does not re-resolve anything, and a resolution can never drift:
# --frozen refuses to update uv.lock.
# ---------------------------------------------------------------------------
COPY pyproject.toml uv.lock README.md ./
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --no-install-project --extra genomics --extra workflow

# ---------------------------------------------------------------------------
# The project itself. Only tracked source: no workspace, no runs/, no fixtures
# beyond the synthetic ones already in the repo.
# ---------------------------------------------------------------------------
COPY src/ ./src/
COPY config/ ./config/
COPY knowledge/ ./knowledge/
COPY workflow/ ./workflow/
COPY Snakefile ./Snakefile

# Jinja report templates. Four .j2 files, and the reporting stage cannot run
# without them: `render.py` resolves them as `parents[3]/"templates"` relative to
# its own location, which is /app/templates for the editable install below. The
# dossier stage raises FileNotFoundError without this line, after the whole
# pipeline has already run.
COPY templates/ ./templates/

# ---------------------------------------------------------------------------
# The two synthetic demo fixtures. THE ONLY GENOMIC FILES IN THE IMAGE.
#
# Both are fully fabricated -- the gene symbols are fictional `SYNTH*` names --
# which is what makes them safe to build in when no real input ever is. They are
# here so the image can run its own demo without the caller having to supply a
# case first.
#
# NAMED ONE FILE AT A TIME, never `COPY tests/`. The enumeration IS the control:
# a glob or a directory copy would silently absorb whatever a future commit puts
# next to them, and the sibling `tests/fixtures/synthetic/consequence/` directory
# is already there. `.dockerignore` enforces the same restriction a second time,
# from the other direction, by keeping everything else under `tests/` out of the
# build context entirely.
# ---------------------------------------------------------------------------
COPY tests/fixtures/synthetic/synthetic_case.vcf \
     tests/fixtures/synthetic/synthetic_phenotype.tsv \
     ./tests/fixtures/synthetic/

# An EDITABLE install, deliberately, and this is load-bearing rather than a
# convenience.
#
# `mva` locates the repository by walking up from its own `__file__` to the
# directory holding `pyproject.toml` (`find_repo_root` in src/mva/config.py), and
# the result is not cosmetic: it is where the pipeline reads `knowledge/`
# (orchestrator), where the resource manifest is resolved, and what the privacy
# audit scans. With `--no-editable` the package is installed as a built copy
# under /opt/venv, nothing above it holds a `pyproject.toml`, and every run fails
# at startup with "Could not locate repository root".
#
# Editable keeps exactly one copy of the code, at /app/src, where the repo root
# resolves to /app -- and /app is read-only to the running user two layers down,
# so this costs nothing in provenance. `--frozen` still refuses to update
# uv.lock, so the resolution cannot drift either way.
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --extra genomics --extra workflow

# ---------------------------------------------------------------------------
# `docker run IMAGE demo` must stay the same run as `just demo`.
#
# The entrypoint restates the demo's arguments because it cannot call `just`
# (that recipe shells out to `uv run` and targets $TMPDIR, not the mount), and a
# restatement drifts. This assertion is what stops it: it compares the `mva run
# all` arguments and the seeded fixture list between the justfile's recipes and
# the entrypoint, and FAILS THE BUILD if they disagree.
#
# The justfile is bind-mounted for the check rather than copied, so the image
# does not carry a second, staleable copy of it.
# ---------------------------------------------------------------------------
RUN --mount=type=bind,source=justfile,target=/tmp/justfile \
    bash /app/workflow/containers/assert-demo-parity.sh \
        /tmp/justfile /app/workflow/containers/entrypoint.sh

# ---------------------------------------------------------------------------
# Non-root. The container writes only to the bind-mounted workspace, so it must
# not be able to write to /app either — a process that cannot modify the code it
# is running is a process whose provenance record stays true.
# ---------------------------------------------------------------------------
RUN groupadd --system --gid 10001 mva \
 && useradd --system --uid 10001 --gid mva --home-dir /home/mva --create-home mva \
 && chown -R root:root /app \
 && chmod -R a-w /app \
 && chmod 0555 /app/workflow/containers/entrypoint.sh \
 && mkdir -p /workspace \
 && chown 10001:10001 /workspace

USER 10001:10001

# The workspace is mounted here at runtime, read-only for inputs. Nothing is
# baked in; the directory exists only so the mount target is owned correctly.
ENV MVA_WORKSPACE=/workspace
VOLUME ["/workspace"]

# TMPDIR inside the container, never the host: OS temp files are a disclosure
# path (PRIV-08).
ENV TMPDIR=/workspace/tmp

# The entrypoint is a thin shim over `mva`, not a replacement for it: every
# argument except the literal `demo` is passed straight through, so
# `docker run IMAGE run all --config ...` and a bare `docker run IMAGE` (which
# still reaches `mva --help` via CMD) behave exactly as they did when this was
# `ENTRYPOINT ["mva"]`. `demo` is the one added verb:
#
#   docker run --rm -v "$PWD/out:/workspace" mva-research:local demo
#
# See workflow/containers/entrypoint.sh for what it seeds, and README.md in this
# directory for what it produces.
ENTRYPOINT ["/app/workflow/containers/entrypoint.sh"]
CMD ["--help"]
