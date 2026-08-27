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
# BASE IMAGE PIN. `python:3.12-slim` is a moving tag. Pin the digest below and
# reference it in the FROM line before publishing anything reproducible:
#
#   docker buildx imagetools inspect python:3.12-slim \
#     --format '{{ range .Manifest.Manifests }}{{ if eq .Platform.architecture "amd64" }}{{ .Digest }}{{ end }}{{ end }}'
#
#   # pinned digest (linux/amd64, verified <YYYY-MM-DD>):
#   # sha256:<PASTE HERE>
#   # FROM --platform=linux/amd64 python:3.12-slim@sha256:<PASTE HERE> AS base
#
# The digest is left unfilled deliberately rather than pasted from memory: a
# digest that was never verified against a registry is worse than none, because
# it looks pinned. Fill it in the same commit that first publishes the image.
# =============================================================================

FROM --platform=linux/amd64 python:3.12-slim AS base

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

RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --no-editable --extra genomics --extra workflow

# ---------------------------------------------------------------------------
# Non-root. The container writes only to the bind-mounted workspace, so it must
# not be able to write to /app either — a process that cannot modify the code it
# is running is a process whose provenance record stays true.
# ---------------------------------------------------------------------------
RUN groupadd --system --gid 10001 mva \
 && useradd --system --uid 10001 --gid mva --home-dir /home/mva --create-home mva \
 && chown -R root:root /app \
 && chmod -R a-w /app \
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

ENTRYPOINT ["mva"]
CMD ["--help"]
