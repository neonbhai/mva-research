# ADR 0006 — Patient data lives in an external workspace

**Status:** accepted · **Date:** 2026-08-27

## Context
The convenient layout puts data next to code: `./data/`, `./workspace/`,
`./runs/`. It is also how patient data ends up in git history, where `rm` does
not remove it and a later push discloses it permanently.

## Decision
The workspace is an **external, absolute path** supplied via `MVA_WORKSPACE` or
`--workspace`, never a committed config value. `resolve_workspace()` refuses:
- a path that resolves inside the repository (symlinks resolved first);
- a path under a cloud-synced root.

All input paths in a case config are **workspace-relative**, so no committed file
and no provenance manifest records where patient data lives on disk.

## Why the cloud-sync check exists
On macOS, `~/Desktop` and `~/Documents` are synced to iCloud Drive by default
whenever "Desktop & Documents Folders" is enabled — the default on a new machine.
A VCF placed there is uploaded within seconds and is outside the researcher's
control from that moment. This is the single most likely accidental disclosure
path on this platform, and it is silent. We also reject Dropbox, Google Drive,
OneDrive, Box and `Library/CloudStorage`.

## Consequences
- Running the pipeline requires one environment variable. Worth it.
- The demo needs a workspace too; it creates a temporary one and the audit
  records that the in-repo escape hatch was used.
- Defence in depth: the `.gitignore` denies genomic extensions repo-wide anyway,
  and the privacy audit checks tracked files, which `.gitignore` does not protect.
