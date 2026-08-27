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

CI builds it; developers on Apple Silicon generally should not.

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
