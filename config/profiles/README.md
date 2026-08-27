# Snakemake profiles

A profile is a file of default command-line flags. It exists so that "how we run
the workflow" is reviewable in the repository instead of living in one person's
shell history.

## What is here

| File | Purpose |
|---|---|
| `local.yaml` | Single-machine, single-core, fail-fast local run. The only profile the synthetic demo needs. |

## Using it

Snakemake resolves `--profile DIR` to a **directory** containing `config.yaml`
(or `config.v8+.yaml`), not to a file. `local.yaml` is kept as a plain file here
because it is primarily documentation of the defaults; to use it as a profile,
link it into the shape Snakemake expects:

```bash
mkdir -p config/profiles/local
ln -sf ../local.yaml config/profiles/local/config.yaml
uv run snakemake --profile config/profiles/local --configfile config/synthetic-case.yaml
```

Today `just snakemake` passes the equivalent flags explicitly:

```bash
just snakemake                 # uv run snakemake --cores 1 --configfile config/synthetic-case.yaml
just snakemake -n              # dry run; any extra ARGS are forwarded
just dag                       # render the DAG
```

Explicit flags are preferred for the demo because the command a reviewer copies
should say what it does without requiring them to open a second file.

## What a profile must never contain

**The workspace.** `MVA_WORKSPACE` is an external absolute path supplied by the
environment and is never a committed config value (ADR 0006, GP-40). A profile
that pinned a workspace path would record where patient data lives on someone's
disk, in git, permanently — which is the exact failure the external-workspace
decision exists to prevent.

**Credentials or endpoints.** Nothing on the patient-data path talks to a network
service (PRIV-05). A profile with a remote executor or a storage endpoint in it
would be the first place that stops being true.

## Adding a cluster profile

If this ever needs to run on a shared cluster, the new profile belongs here
beside `local.yaml`, and the review question is not "does it work" but "who else
can read the workspace". A scheduler that stages files through a shared spool
directory moves patient data outside the boundary `resolve_workspace()` enforces,
and the `localrules:` list in the `Snakefile` exists so the cheap rules — the
gates and join points — never get submitted anywhere in the first place.
