# Containers

Two environments, deliberately different, for two different jobs.

## The development profile: native on macOS / Apple Silicon

**Local development does not use this container.** It runs natively through `uv`
on arm64, because the whole stack installs cleanly there — verified, not assumed:

| Package | arm64 wheel on macOS | Notes |
|---|---|---|
| `cyvcf2` | yes | vendors htslib; no Homebrew prerequisite |
| `pysam` | yes | same |
| `snakemake` | yes | pure Python; it is a library, not a service (ADR 0001) |
| `duckdb`, `pyarrow`, `pydantic` | yes | — |

That is the reason ADR 0001 chose Snakemake over Nextflow: no JVM, no daemon, and
nothing that has to be emulated. `just bootstrap && just verify && just demo` is
the loop, and it is a native-speed loop.

Running the pipeline inside a linux/amd64 container on Apple Silicon means
qemu emulation for every stage. It is several times slower, it changes the
platform pyright is configured for (`pythonPlatform = "Darwin"`), and it puts a
filesystem translation layer between the code and the workspace — which is a
correctness hazard exactly where correctness matters, because case sensitivity
and `mtime` granularity differ across that boundary and Snakemake's
up-to-dateness checks are `mtime`-based.

So: **arm64 native for development, linux/amd64 container for reproduction.**

## The reproduction profile: linux/amd64 container

The image exists so that a third party — the challenge organizers, who state
they may rerun submissions — can reproduce a run without reproducing our laptop.
It pins the base image by digest, installs from `uv.lock` with `--frozen`, sets
`PYTHONHASHSEED=0` and `TZ=UTC`, and runs as a non-root user.

```bash
docker buildx build --platform linux/amd64 -f workflow/containers/Dockerfile -t mva-research:local .
```

CI builds it; developers on Apple Silicon generally should not. The base image is
pinned by digest to the linux/amd64 manifest of `python:3.12-slim`, read back
from Docker Hub rather than recalled; the Dockerfile header records the date and
the command to re-verify it.

## Running the synthetic demo in one command

```bash
mkdir -p ~/mva-demo
docker run --rm -v "$HOME/mva-demo:/workspace" mva-research:local demo
```

`demo` is the one verb the entrypoint adds. It copies the two synthetic fixtures
into `inputs/` and runs the synthetic case end to end, which takes about 13
seconds under emulation on Apple Silicon and writes 33 artifacts:

| What | Path under `runs/<run-id>/` | Scale |
|---|---|---|
| Track 1 submission | `submission/track1_submission.csv` | 10 rows |
| Ranked candidate pairs | `candidates/ranked_pairs.json` | 18 pairs |
| Dossier, mechanism, Track 2, drug-hypothesis and rejection reports | `reports/` | 5 files |
| Evidence database | `evidence/evidence.duckdb`, `evidence/parquet/` | 17 tables |
| Discard ledger (GP-13) | `discards/discard_ledger.tsv` | — |
| Privacy audit | `privacy/privacy_audit.md` | — |
| Run provenance and tool versions | `provenance.json` | — |

Two details that will otherwise cost you a few minutes. Docker Desktop rejects a
bind mount whose host directory is absent instead of creating it, so the `mkdir`
is load-bearing. And any host directory works, `-v "$PWD/out:/workspace"`
included, because the fixtures are fabricated; the example keeps the workspace
outside the repository anyway, for the same reason ADR 0006 keeps the real one
there. The habit is the control, and `out/` is not in `.gitignore`.

The run id is derived from the case rather than the clock, so it stays
`synthetic-demo-964738e1d69a` across runs, and a second `demo` into the same
directory overwrites that run in place. Unlike `just demo`, which depends on
`clean-demo`, the container does not clear the workspace first: /workspace
belongs to the caller, and a container that `rm -rf`s its own mount is one nobody
should later aim at a real case.

Only the literal `demo` is intercepted. Everything else is handed to `mva`
unchanged, so the read-only-inputs invocation below, a bare `docker run IMAGE`
(which still reaches `mva --help`), and `docker run IMAGE verify determinism ...`
all behave as they did before the verb existed.

### Why the container demo cannot drift from `just demo`

The image cannot simply call `just`. That recipe shells out to `uv run`, which
wants to sync a project environment the image has already frozen, and it writes
to `$TMPDIR/mva-research-demo` rather than to the mount. So
`entrypoint.sh` restates the run, and because a restatement is a copy that drifts,
`assert-demo-parity.sh` runs during the build: it extracts the `mva run all`
arguments and the seeded fixture list from the justfile's own `demo` and
`_seed-workspace` recipes, compares them against the entrypoint, and fails the
build when they disagree. Change the case config on one side only and you get a
red build with both versions printed, rather than a demo that quietly analyses a
different case. The extraction is scoped to those two recipes by name, so a new
recipe that also invokes the pipeline does not trip it.

Argument parity is checked mechanically on every build. Output parity was checked
by running both: linux/amd64 under qemu and `just demo` natively on macOS arm64
produce the same run id and byte-identical submission, reports, ranked pairs,
annotated variants, triage report and discard ledger. Since the pipeline draws on
no RNG, no set-iteration order and no time-derived id (GP-30), the architecture
underneath does not reach the output — which is the property the organizers are
being asked to take on trust, and the one this demo lets them check in one
command.

## Patient data and images: the rule

**Patient data is bind-mounted read-only at runtime. It is never baked into an
image.**

```bash
docker run --rm \
  --network none \
  --mount type=bind,source="$MVA_WORKSPACE/inputs",target=/workspace/inputs,readonly \
  --mount type=bind,source="$MVA_WORKSPACE/runs",target=/workspace/runs \
  mva-research:local run all \
    --config config/synthetic-case.yaml \
    --defaults config/default.yaml \
    --workspace /workspace
```

Why this is a hard rule and not a preference:

- **An image containing patient data is a distributable copy of that data.**
  Images are pushed to registries, pulled by digest, cached across machines and
  shared by URL. Every one of those is a disclosure event. `docker save` turns
  the whole thing into a single portable file.
- **Deleting the file does not delete the data.** Docker layers are immutable and
  content-addressed. A `COPY case.vcf` followed by a `RUN rm case.vcf` in a later
  layer leaves the bytes intact and readable in the earlier layer, in exactly the
  way `rm` after `git add` leaves them in git history (ADR 0006 makes the same
  argument about the repository).
- **The deletion obligation cannot be met.** The challenge requires deletion
  within 30 days of close, from all environments including derived datasets
  (`docs/references/track1-submission-contract.md`). A published image is not an
  environment we control.

Concretely, in this Dockerfile: no `COPY` reads anything outside the repository,
there is no `ADD` from a URL, and the only writable path at runtime is the
mounted workspace. `--network none` is worth adding on every real-data run;
`REF_PATH=/dev/null` is already set in the image because htslib will otherwise
fetch CRAM reference sequences from `www.ebi.ac.uk` when it opens a patient file
(PRIV-06).

The inputs mount is `readonly` on purpose: a stage that cannot write to the input
directory cannot accidentally write a derived artifact next to the source data,
where it would inherit none of the workspace's handling.
