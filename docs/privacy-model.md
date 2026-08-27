# Privacy model

The asset is one child's clinical VCF, phenotype and raw reads. It is
re-identifiable and the disclosure is **not revocable**. The adversary model is
mostly *accident* — the agent, the shell, git, the network stack, and macOS
itself — rather than a targeted attacker.

Boundary: `$MVA_WORKSPACE`, an external absolute path (ADR 0006).

---

## Threat table

| ID | Threat | Path | L | I | Control |
|---|---|---|---|---|---|
| PRIV-01 | Patient file committed | `git add -A` from repo root | H | Critical | Deny-by-default `.gitignore`; pre-commit hook running `mva privacy audit --staged --strict`; `git config core.hooksPath .githooks` |
| PRIV-02 | Agent/subagent reads patient files | Agent globs `**/*.vcf`; workspace symlinked into repo | H | Critical | Workspace outside repo, path only via env; `workspace_containment` + `symlink_escape` checks; ADR 0008 (illegibility by design) |
| PRIV-03 | Terminal output enters model context | `head file.vcf`, dataframe print, pytest failure diff | H | Critical | No stage prints records; CLI is counts/paths only; `safe_repr()`; pytest `--tb=line`, never `--showlocals` |
| PRIV-04 | Logs dump VCF records | `log.debug("variant=%s", rec)` | H | High | `GenomicRedactionFilter` on **every handler** + `setLogRecordFactory`; `log_redaction_probe` check |
| PRIV-05 | Annotation APIs leak coordinates | VEP REST, ClinVar, gnomAD, MyVariant | M | Critical | Structural test forbids network-client imports on the sensitive path; annotation is local, hash-pinned tables only |
| PRIV-06 | CRAM reference auto-fetch | htslib `REF_PATH` defaults to `www.ebi.ac.uk/ena/cram/md5/%s` | M | Med | `REF_PATH=/dev/null`, `REF_CACHE` inside workspace, explicit `-T ref.fa` |
| PRIV-07 | Cloud sync | `~/Desktop`, `~/Documents` are iCloud-synced by default | H | Critical | `resolve_workspace` rejects synced roots; `cloud_sync_location` check |
| PRIV-08 | OS temp files | `tempfile` defaults to `/var/folders/…` | M | High | `TMPDIR` pinned inside the workspace at CLI entry |
| PRIV-09 | Crash dumps / tracebacks embed genotypes | Exception carrying record locals | M | Critical | Exceptions built from IDs and counts only; `exc_info` redacted by the log filter; no `rich` locals |
| PRIV-10 | Reports embed raw coordinates | Templates rendering variant rows | H | High | Two render targets: internal (workspace-only) and public (allowlisted fields); public output re-scanned after render |
| PRIV-11 | Notebook outputs | `.ipynb` cell outputs | M | Critical | `*.ipynb` gitignored; `notebook_output_purity` check |
| PRIV-12 | LLM prompt content | Pasting a VCF slice to debug | H | Critical | Policy + design: debugging is always on the synthetic case (ADR 0007) |
| PRIV-13 | The public export IS the leak | Submission CSV / report | H | Critical | Allowlist + per-column validators + post-render content re-scan, fail-closed |
| PRIV-14 | Time Machine backups | Local APFS snapshots pin deleted blocks | H | High | `tmutil addexclusion`; encrypted volume; snapshot purge at teardown |
| PRIV-15 | Spotlight indexing | `.Spotlight-V100` stores content excerpts | M | Med | Mount `-noindex`; `mdutil -i off`; `.metadata_never_index` sentinel |
| PRIV-16 | Deletion obligation unmet | Residue in snapshots, swap, agent transcripts | H | High | Cryptographic erasure (below); retention manifest with a delete-by date |
| PRIV-17 | Git history poisoning | File committed then `rm`'d — still in history and packfiles | M | Critical | Treat any hit as an incident: `git filter-repo`, expire reflog, `gc --prune=now`. If ever pushed, consider it disclosed. |

---

## The scanner must not become the leak

The privacy scanner is run by an agent, and its output enters model context. So
the scanner's own output contract is a control, not a nicety (GP-41):

- Never emit `match.group()`, matched bytes, or surrounding context.
- Permitted derivations of a match: `match.span()`, `len(match.group(0))`, the
  rule ID, the line number, the file path.
- Fixed placeholder form: `<REDACTED:{rule_id}:len={n}>`.
- Correlation IDs, if needed, use an HMAC with an ephemeral per-process salt.
  A plain truncated hash of a low-entropy value (an HPO term, a short MRN) is
  brute-forceable.
- Read as bytes with a size cap and `errors="replace"`. `UnicodeDecodeError`
  embeds the offending bytes in its own message, so every read is wrapped and
  re-raised scrubbed.

---

## Why the log filter, not a formatter

A `Formatter` is per-handler and optional. A handler added later with a default
formatter, `logging.lastResort`, or third-party code calling
`record.getMessage()` all bypass it. A `Filter` mutates the `LogRecord` itself,
so every downstream consumer — including handlers we never configured — sees
only redacted text.

The filter is attached to **every handler**, not to the logger: a filter on a
logger is consulted only for records logged *directly* to it. Records
propagating up from child loggers skip ancestor filters entirely. That is the
classic hole, and it is exactly how a library's debug log would escape.

---

## Network denial: honest limits

`sys.addaudithook` is the stronger guard because it fires inside CPython's
`_socket` module and cannot be defeated by a library holding a pre-bound
`from socket import socket` reference. Audit hooks cannot be removed, so the
profile gates on a module flag rather than installing and uninstalling.

**It is a tripwire and a developer guardrail, not a security boundary.** It is
blind to:
- C extensions — `pysam`/htslib and `cyvcf2` call `connect(2)` directly and raise
  no Python audit event;
- subprocesses, once spawned;
- `ctypes`/`cffi` calling libc directly.

The real boundary is OS-level. For a real-data run, use one of:
```bash
sandbox-exec -p '(version 1)(allow default)(deny network*)' <cmd>
networksetup -setairportpower en0 off
```
Whichever is used is recorded in the run manifest's `network_profile`.

---

## Secure deletion on macOS APFS

What does **not** work, stated plainly:

- `rm` unlinks a directory entry. The blocks keep their contents until reused.
- `rm -P` overwrites *logical* blocks. APFS is copy-on-write, so the overwrite
  lands on **new** blocks and the originals are simply released. `srm` was
  removed from macOS in 10.11. Overwrite-based erasure is false assurance here.
- `diskutil secureErase freespace` is unsupported on APFS and refused on SSDs.
  SSD wear-levelling and over-provisioning mean the flash controller can retain
  copies at physical addresses no LBA can reach. TRIM is best-effort.
- **Local Time Machine snapshots are the biggest real leak**: after `rm`, the
  file remains fully readable inside the snapshot. `.DocumentRevisions-V100` and
  `.Spotlight-V100` are further independent copies.

**The only defensible guarantee is cryptographic erasure.** Set this up on day
one, not at cleanup time.

```bash
# Create an encrypted sparse bundle. Record the passphrase OUT OF BAND -- not the
# login keychain, which is itself backed up and snapshotted.
hdiutil create -size 250g -type SPARSEBUNDLE -fs APFS \
  -encryption AES-256 -stdinpass -volname MVACASE ~/private/mvacase.sparsebundle

sudo tmutil addexclusion -p ~/private/mvacase.sparsebundle
tmutil isexcluded ~/private/mvacase.sparsebundle          # expect [Excluded]

hdiutil attach -stdinpass -nobrowse -noindex \
  -mountpoint /Volumes/MVACASE ~/private/mvacase.sparsebundle
sudo mdutil -i off -d /Volumes/MVACASE
touch /Volumes/MVACASE/.metadata_never_index

export MVA_WORKSPACE=/Volumes/MVACASE/case01
export TMPDIR=/Volumes/MVACASE/tmp
```

Teardown — the ciphertext is what remains, and the key is gone:

```bash
hdiutil detach /Volumes/MVACASE
rm -rf ~/private/mvacase.sparsebundle        # plain rm is correct: this is ciphertext

tmutil listlocalsnapshots /
tmutil listlocalsnapshots /System/Volumes/Data
sudo tmutil deletelocalsnapshots /
sudo tmutil thinlocalsnapshots / 999999999999 4

mdfind -onlyin ~ 'kMDItemFSName == "*.vcf"'   # verify nothing indexed
qlmanage -r cache                             # Quick Look thumbnails
```

Also purge derived copies outside the image: agent transcripts under
`~/.claude/projects/`, shell history entries naming workspace paths, and any
editor search index. Keep FileVault on — it covers swap, which can hold plaintext
from process memory. Stream VCF records rather than loading whole files to shrink
that window.

---

## Challenge deletion obligation

The data-use terms require deletion **within 30 days of the hackathon close
(2026-10-24)**, from all environments including derived datasets, confirmed by
email to the organizers. Concretely that means:

1. Destroy the encrypted volume and its passphrase (above).
2. Purge local snapshots and verify none postdate first ingest.
3. Delete derived artifacts: evidence DB, Parquet exports, internal reports — all
   of which live inside the volume and go with it.
4. Delete agent transcripts and shell history referencing the case.
5. Record the teardown, then send the confirmation email.

Everything in this repository is synthetic and is unaffected by that obligation.
