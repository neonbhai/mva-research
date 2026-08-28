# ADR 0020 — Reference releases are pinned in full and verified by sample

**Status:** accepted
**Date:** 2026-08-28
**Supersedes:** nothing. **Amends:** the integrity model in
`knowledge/adapters/README.md` ("The knowledge manifest and integrity"), which
assumed every pinned artifact is small enough to re-hash on demand.

## Context

The pipeline now depends on public reference releases that live outside the
repository, at `$MVA_RESOURCES`:

| Release | Size |
|---|---|
| gnomAD v4.1 exomes, 24 sites VCF shards + tabix indexes | 184.8 GB |
| GRCh38 no-alt FASTA (decompressed) + `.fai` | 3.14 GB |
| SnpEff GRCh38.115 genome database (28 files) | 0.60 GB |
| ClinVar weekly VCF + index, MANE v1.5, HPO, gene panels | 0.35 GB |
| gnomAD v4.1 constraint metrics | 0.10 GB |
| **Total registered** | **202.8 GB across 89 artifacts** |

`knowledge/manifests/knowledge.yaml` already establishes the rule these have to
obey: nothing is opened as data until its `sha256` has been checked against the
manifest, and a mismatch raises without echoing the file (PRIV-09). That rule was
written for TSVs of a few kilobytes. Applied literally to 202.8 GB it is
unaffordable, and an unaffordable check is one that gets a `--skip` flag within a
week and then gets passed by default.

Two things make this more than a performance footnote.

**First, silent-wrong-bytes is a failure that already happened here.** A
Gene2Phenotype bulk-download URL returned a JavaScript app shell with HTTP 200;
`file -b` reported "HTML document text". A sha256 of that page is a perfectly
valid sha256 — of the wrong thing. Hashing alone would have pinned the error page
forever and every subsequent check would have passed.

**Second, the failure that actually recurs is truncation, not tampering.** These
files arrive over hours, in parallel ranged parts, on a laptop that sleeps. The
manifest in this repository has itself carried entries reading `download in
progress: size grew during a 2s observation window` — that is the real threat
model, and it is one a cheap check can catch completely.

## Decision

**Pin in full at registration. Verify by published sample at run time. Record in
the manifest exactly what was checked, by which method, and when.**

Concretely:

1. **Registration** (`uv run python -m tools.acquire write-manifest`) computes a
   full `sha256` over every byte, runs a **deep format probe**, and stores both
   outcomes plus a sampled digest in each entry's `integrity` record.

2. **Run time** defaults to `IntegrityMode.SPOT`: exact byte size, plus a sha256
   over a fixed, published, reproducible sample of the file —

   ```
   spot-v1:head=8MiB,tail=8MiB,windows=8x1MiB,size-bound
   ```

   the first 8 MiB, the last 8 MiB, and eight 1 MiB windows at evenly-spaced
   interior offsets, with the file's exact size mixed into the digest before any
   content. Files at or below 24 MiB are read whole, so for them SPOT and FULL are
   the same check.

3. **`IntegrityMode.FULL`** re-reads everything and stays one flag away
   (`--mode full`), for a release check or after any suspicion.

4. **There is no `off`.** A mode that skips verification is the one that gets
   selected under deadline pressure (ADR 0009: gates are blocking by design).

5. **A mismatch raises `ResourceIntegrityError`** naming the resource and its
   declared path and quoting both digests, never file contents (PRIV-09). A
   digest is a one-way summary, not content.

6. **A spot digest recorded under an unrecognised `spot_plan` is refused, not
   compared.** Digests under different sampling plans are incomparable; treating a
   difference as corruption would be a false alarm and treating it as a pass would
   be a lie.

## What each mode costs, and what it proves

Measured on this machine, 2026-08-28, over the registered 202.8 GB:

| | bytes read | wall time | proves |
|---|---|---|---|
| `full` | 202.8 GB (100%) | **159 s** | whole-file identity |
| `spot` | 1.27 GB (0.63%) | **1.7 s** | size exactly; sampled content; completeness |

Both figures are end-to-end `uv run python -m tools.acquire verify --mode …` over
all 65 registered resources, both reporting 65/65 OK. `spot` is **94x cheaper**.
Registration itself (full sha256 + deep format probe + index cross-check) is the
same order as `full`.

For the largest single artifact — the 18.79 GB gnomAD chr1 shard — `spot` reads
24 MiB, **0.134%** of it.

`spot` **does** catch:

* a truncated or still-downloading file (the size is bound into the digest, and
  the tail window covers the BGZF end-of-file block that is written last);
* a file replaced by an HTML error page;
* a file swapped for a different release of the same dataset;
* corruption that lands anywhere in the ~24 MiB sampled.

`spot` **does not** prove whole-file identity. Corruption confined to the 99.4%
of bytes it does not read will pass, and an adversary with write access to the
resource root can defeat it deliberately. It is a checksum over a sample, and
this repository will not describe it as anything else.

## A pinned file can still have a lying index

Hashing every byte of a VCF says nothing about whether the `.tbi` beside it
describes *that* VCF. They are two separate downloads. The index can resume from
an earlier attempt, the data can be re-fetched while the index is not, or an
interrupted transfer can leave a complete index next to a shorter file. In every
case both files hash exactly as the manifest says, and random access returns the
wrong records — or none — for precisely the rare coordinates this pipeline
exists to find. A missing frequency join is indistinguishable from "novel and
ultra-rare" (GP-14), so this failure is silent all the way to the submission.

htslib gestures at it: it warns `index file is older than the data file`, which
it currently does for **all 25 gnomAD shards here**, because each small `.tbi`
finished downloading before its multi-gigabyte `.bgz` did. Under this ADR's own
reasoning that warning is not evidence — mtime says nothing about content — and
in this case it is a false alarm from ordinary download ordering. But "the
warning is noise" is not the same as "the failure is imaginary", so the index is
checked properly instead of being either trusted or ignored:

1. **Structural.** Parse the `.tbi` and take the deepest BGZF block offset it
   references. It must lie inside the data file. An index built from a larger file
   blows past the end, and this is arithmetic over ~100 kB of index — no seeking,
   no htslib. (Strictly greater than the file size: an offset *equal* to it is the
   canonical "the last record runs to the end" terminator, which every gnomAD shard
   uses. An earlier `>=` here flagged all 25 as stale and would have failed 184.8 GB
   of good data closed.)
2. **Behavioural.** Seek to the deepest region the index's own linear index claims
   to cover, and require the records that come back to lie inside the window that
   was asked for. Aimed deep on purpose: a fetch near the start is served by the
   first BGZF block, which is where a wrong index is most likely to be accidentally
   right. This is what catches a same-size index built from different data.

The FASTA gets the equivalent treatment: its `.fai` is cross-read against the file
itself, so a declared offset that lands on different bases is caught. A `.fai`
built from a different FASTA still opens and still returns sequence — silently
offset sequence, which produces confidently wrong left-alignment (ADR 0018).

Each entry records the verdict as `integrity.index_check`, one of
`consistent` / `stale` / `not_checked` / `not_applicable`, with
`integrity.index_detail` saying what was observed, **including the mtime skew and
why it was not treated as evidence**. `stale` fails closed: the resource is
recorded `not_fetched` rather than pinned, because a hash of a mismatched pair is
a perfectly valid hash of the wrong thing.

## What was rejected

**Size + mtime alone.** The obvious cheap check, and dishonest. `mtime` is
restored by any tool that rewrites a file (`rsync -a`, a restore from backup, a
re-extraction of an archive), and in-place corruption preserves both fields. It
would report OK for a file whose bytes had changed, which is worse than no check
because it produces a green result. It survives in exactly one place — as the key
of the registration-time digest cache in `tools/acquire/digest.py`, where the only
thing it decides is whether to re-read bytes whose full digest was computed in the
same session, and where `--rehash` bypasses it. That module says so at the top.

**Full hashing on every run.** 159 s before every pipeline invocation, for a
pipeline whose synthetic demo runs in seconds. It would be disabled and then
forgotten.

**Hashing only the small resources.** This inverts the risk. gnomAD and the FASTA
are both the largest artifacts and the ones whose corruption is hardest to notice
downstream: a bad frequency join looks exactly like "novel and ultra-rare"
(GP-14), and a subtly wrong reference produces confidently wrong left-alignment
(ADR 0018).

**A Merkle tree per file.** Strictly better than sampling and genuinely
affordable to verify incrementally, but it needs a full pass to build, a chunk
list to store per file, and a format nothing else in this project reads. The
sampling plan buys most of the detection for one line of manifest.

## Consequences

* Every entry in `knowledge/manifests/resources.yaml` gains an `integrity` block:
  `verified_at`, `spot_plan`, `spot_sha256`, `format_check`, `format_detail`. A
  reader can see what was proven instead of inferring it from a hash's presence.
* Indexed releases additionally carry `index_check` / `index_detail`. A test
  (`test_every_indexed_release_records_its_index_verdict`) fails the build if any
  indexed release is pinned without its index having been proven against its data.
* `format_check` is a closed vocabulary
  (`mva.resources.FormatCheck`) whose default is `not_checked` — never a passing
  state. A resource that fails its format probe is recorded `not_fetched` with the
  reason, so it can never be pinned while being the wrong thing.
* Changing any sampling constant **requires** bumping `SPOT_PLAN` and
  re-registering. Verification refuses an unknown plan rather than comparing
  incomparable numbers.
* The SnpEff genome database is pinned as a `directory` resource under a tree
  digest over sorted `(relative path, size, sha256)` triples. Pinning only
  `snpEffectPredictor.bin` would have left the 25 `sequence*.bin` files unpinned,
  and those supply every codon — without them SnpEff still annotates, but HGVS.c
  and HGVS.p are silently absent from every consequence.
* The resource root has **no default path**. It comes from `$MVA_RESOURCES` or an
  explicit flag, must exist, and must resolve outside the repository — the same
  shape as the workspace rule in ADR 0006, for a different reason (size, not
  sensitivity). The previous hard-coded fallback to one contributor's home
  directory is removed: a guessed root that is wrong fails later, somewhere less
  obvious, and "resource missing" is indistinguishable from "resource not
  registered" by the time it reaches an adapter.

## How to check this decision is still being honoured

```bash
# cheap: what a run does
uv run python -m tools.acquire verify --mode spot

# expensive: what a release should do
uv run python -m tools.acquire verify --mode full
```

Both exit non-zero on any `MISSING` or `MISMATCH`. `PENDING` (honestly declared
not-fetched) never fails either.
