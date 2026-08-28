# Scale report: ingesting a WGS-scale VCF

**Question.** The pipeline has only ever run on a 12-record fixture. A real
300 MB whole-genome VCF with ~4.5 M variants arrives shortly. Nobody has
measured what ingestion does at that scale.

**Answer in one line.** Every algorithm here is linear and nothing is
quadratic — but the design holds the entire callset in memory as Pydantic
models, and at 4.5 M records that needs about **20 GB for `read_vcf` alone and
32 GB to reach the end of QC**, on a machine with 24 GB. It does not fit, and
it fails before annotation is reached.

Everything below was observed on this machine unless explicitly labelled
`projected` or `assumption`. Where something could not be run, that is said
rather than filled in.

Reproduce with `tools/scale/generate_wgs_vcf.py` and `tools/scale/harness.py`.

---

## 0. What was measured against

`tools/scale/generate_wgs_vcf.py` builds a **synthetic throughput phantom**.

| property | value |
|---|---|
| records | 4,492,805 VCF data lines → 4,556,006 `VariantRecord`s after multi-allelic decomposition |
| contigs | 24 (`chr1`–`chr22`, `chrX`, `chrY`), real GRCh38 lengths, `chr`-prefixed |
| density | ~1,500 variants/Mb (1 per ~660 bp); chrX at 0.48x, chrY at 0.06x (male-sample model) |
| sample | one, `SYNTH_PROBAND01`, `FORMAT=GT:AD:DP:GQ` throughout |
| het:hom-alt | 1.49:1 |
| indels | 10.1% |
| multi-allelic | 1.31% of sites (2 or 3 ALTs) |
| non-PASS FILTER | 7.95% (LowQual, VQSR tranches, one two-filter combination) |
| depth | Gamma(9, 3.6): mean ~32x, sd ~11x, right-skewed |
| GQ | 78% at 99 for confident calls, spread below |
| INFO | GATK-shaped, 14 keys — roughly two-thirds of every line |
| line length | 217 bytes mean, uncompressed |
| on disk | 194.8 MB bgzf + 1.8 MB `.tbi`, written through pysam (no `bgzip`/`tabix` binary on this machine) |
| marker | `##mva_synthetic=true` on line 2, verified present in the first 4096 inflated bytes — which is where `mva.privacy.audit` looks, and it does inflate gzip members |
| seed | `20260828`; output is a pure function of `(seed, variants)` |

Files, outside the repo as required, under
`/Users/someshwar-tripathi/Contri/bio-hackathon/mva-resources/scale/`:
`wgs_synth_{100000,500000,1000000,4500000}.vcf.gz` (+ `.tbi`).

> **This file is not biologically valid and must never be described as such.**
> No reference genome was consulted: every REF base is a random nucleotide, not
> the base at that locus. Every coordinate, allele, depth and genotype is
> fabricated. It exists to measure bytes and objects per second. It says nothing
> about whether the science is right, and no result in this document is evidence
> about variants, genes or disease.

**Adjusting to the real file.** The phantom is 195 MB bgzf against the real
file's stated 300 MB at a similar variant count, so the real lines are ~1.5x
longer — almost certainly a longer INFO column. Both backends pay for INFO (the
text backend splits it, htslib parses it), so **read timings below should be
scaled up by roughly 1.5x**. The memory figures are driven by record count, not
line length, and should carry over directly.

---

## 1. Measured throughput and memory

Method: `tools/scale/harness.py sweep`. Each configuration runs in its own child
interpreter, because `ru_maxrss` is a high-water mark and a shared process would
report one backend's peak as the other's. Peak RSS is
`resource.getrusage(RUSAGE_SELF).ru_maxrss`; current RSS is `psutil`. Stages are
held the way `orchestrator.py` holds them — the `read_vcf`, `normalise_variants`
and `assess_quality` results all alive at once.

| file | backend | records out | read | rec/s | peak RSS, read | normalise | QC | peak RSS through QC |
|---|---|---|---|---|---|---|---|---|
| 100k | cyvcf2 | 100,806 | 2.15 s | 46,784 | 0.49 GB | 0.14 s | 1.34 s | 0.86 GB |
| 100k | text | 100,806 | 1.50 s | **67,127** | 0.47 GB | 0.14 s | 1.43 s | 0.85 GB |
| 500k | cyvcf2 | 505,068 | 10.82 s | 46,694 | 2.23 GB | 0.88 s | 7.01 s | 3.70 GB |
| 500k | text | 505,068 | 8.65 s | 58,421 | 2.20 GB | 0.62 s | 6.99 s | 3.68 GB |
| 1M | cyvcf2 | 1,011,274 | 23.94 s | 42,242 | 4.37 GB | 1.72 s | 16.91 s | 6.13 GB |
| 1M | text | 1,011,274 | 19.33 s | 52,307 | 4.37 GB | 1.44 s | 18.14 s | 6.13 GB |
| 1M re-run § | cyvcf2 | 1,011,274 | 23.77 s | 42,545 | 4.39 GB | 2.30 s | 18.54 s | **7.25 GB** |
| 1M re-run § | text | 1,011,274 | 18.29 s | 55,298 | 4.37 GB | 1.71 s | 10.97 s | **7.28 GB** |
| **4.5M** | **cyvcf2** | **4,556,006** | **453.61 s** | **10,044** | 6.80 GB † | not run ‡ | not run ‡ | — |
| **4.5M** | **text** | **4,556,006** | **142.35 s** | **32,005** | 6.55 GB † | not run ‡ | not run ‡ | — |

† These two peaks are **not** the memory the process wanted. See §2.

‡ Read-only at full scale, deliberately: the read alone drove the machine into
sustained paging, and running normalise+QC on top would have added another
~12 GB of demand to a machine already 10 GB oversubscribed. Their cost is
projected from the clean linear points instead, and labelled as a projection.

§ **Code moved under the measurement.** `src/mva/ingestion/normalise.py` and
`src/mva/models/variant.py` were rewritten by a concurrent agent at 04:57–04:58,
after the first sweeps, and a new `src/mva/alleles.py` appeared with them.
`reader.py`, `qc.py` and `models/genome.py` were **not** touched. The 1M point
was therefore re-measured against the current tree; that is the `re-run` row and
it is the one the projections use. Read throughput is unchanged (within run
noise), which is why the full-file read figures — taken before the change —
remain representative. What did change is memory: **peak through QC rose from
6.13 GB to 7.25 GB per million records, +18%.**

Both backends produce **identical** record counts, identical normalisation
operation counts (`split_multiallelic: 26,649`, `trim: 8,775` at 1M) and
identical warnings at every size, including at full scale
(`no_call_genotype (n=13553)`, 0 skipped, both backends). The equivalence the
reader claims holds at WGS scale.

### The text backend is faster than cyvcf2, by a lot

Not a rounding difference: **1.4x faster at 100k, 1.24x at 1M, and 3.2x faster
on the full file**. `detect_backend()` prefers cyvcf2 whenever it is importable,
so the default path is the slow one.

The profile says why. Per record the cyvcf2 path calls `variant.format(key)`
four times — DP, AD, GQ, PS — and each returns a fresh NumPy array that is
immediately unpacked into a Python tuple and thrown away, plus `variant.genotypes`
to rebuild the GT **string** that htslib had already parsed into integers:

```
397492  0.207  reader.py:483(_cyvcf2_format)     <- 4 calls per record
397492  0.202  reader.py:499(_cyvcf2_vector)     <- 0.641s cumulative
1096666 0.088  {built-in method builtins.getattr}
```

The adapter converts htslib's fast binary representation back into the text
representation the shared `_RawRecord` wants. htslib's advantage is spent
undoing htslib. `_parse_with_cyvcf2` costs 1.820 s cumulative against
`_parse_with_text`'s 1.275 s on the same 100k file.

---

## 2. The callset is fully materialised — and that is the failure

`read_vcf` returns `IngestionResult.variants: tuple[VariantRecord, ...]`. Every
record for the whole genome is alive at once. So is every record returned by
`normalise_variants`, and every record plus one `EvidenceItem` each from
`assess_quality`. This is a design property, not a tuning knob.

It gets worse inside `read_vcf` itself: `_parse_records` returns a **complete**
`tuple[_RawRecord, ...]` for the entire file, and that tuple stays alive for the
whole `_records_from_raw` loop that builds the `VariantRecord` list from it.
Two full representations of the callset coexist. Nothing streams and nothing is
released early.

Measured cost per record, from the 1M re-run against the current tree (38 MB
baseline subtracted):

| quantity | measured | per record |
|---|---|---|
| resident after `read_vcf` | 4.29 GB | **4.20 KB** |
| peak during `read_vcf` | 4.37 GB | 4.29 KB |
| peak through `assess_quality` | 7.21 GB | **7.13 KB** |

`harness.py modelbench` attributes it. Building 500,000 bare `VariantRecord`s
with nothing else alive, run twice:

```
construction : 4.29s / 7.33s = 116,464 / 68,251 records/s   (noisy)
resident cost:         1.84GB =  3,672 bytes/record          (identical both runs)
model_dump   :   197,422 / 175,958 records/s
sort_key()   : 2,583,220 / 1,827,842 calls/s
```

**3,672 bytes for one record** — three nested Pydantic models
(`VariantRecord` + `GenomicCoordinate` + `Genotype`), each with its own
`__dict__` and `__pydantic_fields_set__`. That figure was byte-identical across
both runs and across the `variant.py` rewrite. The reader's intermediates add
the remaining ~500 bytes. **The model is the cost; the parser is not.**

### Why the full-file peak RSS reads low

The 4.5M rows above show peaks of 6.55 and 6.80 GB, which is *less* than the
6.13 GB the 1M run needed through QC and far below what 4.5x the data implies.
That is not the process being frugal. `ru_maxrss` counts **resident** pages, and
under pressure macOS caps residency by compressing and paging out the rest. The
evidence that the demand was real:

- **Swap grew from 3.1 GB to 33.0 GB** across the two full-file reads, sampled
  every 10 s while they ran.
- **Throughput collapsed exactly where paging began.** Text fell from 52,307 to
  32,005 rec/s (−39%) and cyvcf2 from 42,242 to 10,044 rec/s (−76%) between the
  1M and 4.5M files. Both files are the same shape; only the working set
  changed. cyvcf2 degrades far worse because its per-record NumPy allocate/free
  churn is brutal against a swapped-out heap.
- The 100k/500k/1M points, none of which paged, are **linear to within 2%** in
  memory: 0.49 / 2.23 / 4.37 GB peak for 1x / 5x / 10x the records.

So the honest figures for the real file (4,556,006 records), projected from the
last clean point:

| stage | projected peak | this machine has |
|---|---|---|
| `read_vcf` alone | **~19.6 GB** | 24 GB |
| through `assess_quality` | **~32.5 GB** | 24 GB |

`read_vcf` alone would fit only if nothing else were running. The ingest stage
as the orchestrator composes it would not fit at all.

### Projected wall time for the real 300 MB file

Two answers, because they differ by 6x and the difference is entirely memory:

| condition | text backend | cyvcf2 backend |
|---|---|---|
| measured on this 24 GB machine, paging (read only) | **2 min 22 s** | **7 min 34 s** |
| ×1.5 for the real file's longer lines | ~3 min 30 s | ~11 min |
| if the working set fitted in RAM (at the clean 1M re-run rate) | ~82 s | ~107 s |
| + `normalise_variants`, projected from the 1M re-run | +8 s | +10 s |
| + `assess_quality`, projected from the 1M re-run | +49 s | +84 s |
| **total ingest, if it fitted** | **~2 min 20 s** | **~3 min 20 s** |

**Nobody should plan on six hours.** The pipeline is not slow. Read, normalise
and QC together are a 2–3 minute job at 1M-scale rates. Every minute beyond that
is paging, and paging is the memory design surfacing as a time cost.

(`assess_quality` timing is the noisiest measurement here — 10.97 s and 18.54 s
for byte-identical work on the same 1M input, in two child processes minutes
apart. That spread is GC and memory-pressure noise on an already-loaded machine,
not a property of the code; the QC projections use the slower figure where the
distinction matters.)

---

## 3. The bottleneck, from the profile rather than from a guess

`harness.py profile` over 100k records, text backend, sorted by `tottime`
(3.485 s under the profiler against 1.50 s without it):

```
ncalls   tottime  cumtime  function
302418     0.351    0.697  pydantic_core SchemaValidator.validate_python
 99373     0.278    1.129  reader.py:388(_raw_from_text_line)
302418     0.271    0.968  pydantic/main.py:253(BaseModel.__init__)
300985     0.153    0.273  genome.py:75(normalise_contig)
2500378    0.147    0.147  str.strip
 99373     0.133    2.018  reader.py:540(_records_from_raw)
100806     0.105    1.609  reader.py:599(_build_record)
502597     0.086    0.086  re.Pattern.match
```

**There is no hotspot.** The largest single family is Pydantic model
construction — `302,418 = 3 × 100,806`, one validation each for
`GenomicCoordinate`, `Genotype` and `VariantRecord` — at 0.968 s cumulative,
about 28% of the profiled read. `_records_from_raw` (the model-building half)
is 2.018 s cumulative against `_parse_with_text`'s 1.275 s: **building the
objects costs more than parsing the text.**

Two smaller, fixable items visible in the same profile:

- `normalise_contig` runs **3x per record** (300,985 calls for 100,806 records)
  and each call is a regex match plus string surgery — 502,597 `re.match` calls
  in total. It is called from `_records_from_raw`, again from
  `GenomicCoordinate`'s `_canonical_contig` validator on the value that was just
  normalised, and again from `contig_sort_key` during the sort. Two of the three
  are re-deriving a known answer. A 24-entry dict lookup would erase it.
- `str.strip` is called 2.5 M times — 25 per record.

**Verdict on complexity: linear, and comfortably so.** Every measured curve is
linear or sub-linear in record count; the only super-linear term anywhere is GC
pressure, which is a symptom of §2 and not a separate defect. Nothing here is
quadratic and nothing here is fatal on time. This is "slow but linear and
acceptable" — provided the working set fits, which it does not.

---

## 4. Intermediate artifacts

`json.dumps([v.model_dump(mode="json") for v in variants])` measures **577 bytes
per record**, stable at every file size. For 4,556,006 records:

- **2.63 GB for one variants artifact.**
- `orchestrator.py` writes `variants/normalised.json` and later
  `variants/annotated.json`. The annotated one is strictly larger — it adds
  consequences, population frequencies and clinical assertions per record. A
  realistic total for the pair is **6–10 GB of JSON**.

The disk cost is not the dangerous part. `PipelineContext.write_json_artifact`
calls `canonical_json(payload)`, which builds **the whole artifact as a single
Python `str`** before a byte is written. Serialising the normalised callset
therefore needs, simultaneously: 4.5 M `VariantRecord`s (~16.5 GB), 4.5 M plain
`dict`s from `model_dump`, and one 2.6 GB string. DuckDB and pyarrow are already
first-class dependencies and would write this incrementally.

---

## 5. Post-filter candidates, and whether the per-gene cap matters

`harness.py pairs`. The cascade rates are **assumptions**, marked as such, chosen
to reproduce the two figures that are stable across diagnostic WGS practice —
roughly 20–24 k coding+splice variants per genome, of which a few hundred
survive a 1% rarity cut. They are declared as rates in `CASCADE` so a reader who
disagrees can substitute and re-run.

```
4,500,000 records
  x 0.92    PASS                 ->  4,140,000
  x 0.0053  coding or splice     ->     21,942
  x 0.035   gnomAD popmax <= 1%  ->        768
  x 0.60    missense/LoF/splice  ->        461
```

**≈460 candidate variants survive a standard rare-disease cascade.** Spread over
19,500 genes that is 0.024 per gene, and `generate_pairs` would build ~466
objects. **The cap of 20 per gene never fires** — not once.

It fires in exactly one place: the heavy tail the uniform model cannot see. TTN,
MUC16, OBSCN, NEB, RYR1 and the HLA locus reliably carry 4–14 rare coding
variants each in any single genome. Across ten such genes the cap truncates 6 of
them and discards 146 of 316 hypotheses. That is where the cap earns its keep,
and it is a small, bounded, real effect.

### But that is not what the pipeline feeds it

`apply_hard_filters` deliberately removes only invalid / hom-ref / no-call
records — filtering is not ranking (GP-13, ADR 0005). There is no rarity or
coding hard filter before pairing. So `generate_pairs` receives **every
alt-carrying annotated variant that has a gene symbol**.

Measured, not assumed: `ManeGeneIndex.genes_at` over the phantom's positions
returns **0.449 genes per variant**, so ~44.9% of a WGS callset falls inside a
MANE gene locus. (The phantom's positions are uniform along each contig, so this
measures the fraction of the *assembly* that is genic; real variants are close
enough to uniform genome-wide for that to be a fair proxy, but it is a proxy.)
That is **2,020,500 variants**, mean 103.6 per gene:

```
objects built    :  7,402,500      <- before the cap
kept after cap=20:    390,000
genes truncated  :     19,500      <- every single one
discarded        :  7,012,500      (94.7%)
```

**Two consequences, one measured and one structural.**

*Time is fine.* `harness.py pairbench` measured `generate_pairs` directly:

```
one gene, growing            spread over 19,500 genes
  1,000 ->  0.008s             100,000 ->  1.91s
  4,000 ->  0.028s             200,000 ->  7.63s   (4.0x for 2x input)
 16,000 ->  0.119s             400,000 -> 27.28s   (3.6x for 2x input)
 64,000 ->  0.520s             800,000 -> 41.53s   (1.44x for 2x input)
128,000 ->  1.030s
```

The combinatorial risk **has already been fixed**, by a concurrent agent, in
`pairing.py` at 04:36 today (ADR 0013): `max_pairing_variants=24` bounds the
pairing *input* per gene, so pair enumeration is at most C(24,2)=276 per gene however
many variants it holds. The growth is quadratic only while the per-gene load is
below 24 (visible in the 4.0x and 3.6x ratios above, at means of 10 and 20 per
gene) and then saturates — the 1.44x at 800,000 is the cap taking hold. Before
that change the same input would have enumerated **~104 million** pairs. At the
full 2.02 M genic load the projection from the 800 k measurement is **~50
seconds**. Not a problem.

*The flag is the problem.* Because every one of the 19,500 genes exceeds the
cap, **every candidate the pipeline emits carries `gene_pair_cap_truncated`**.
A flag raised on 100% of output distinguishes nothing, and its whole purpose —
"a silently shortened hypothesis list is indistinguishable from a complete one"
— is defeated. 94.7% of hypotheses are discarded by a cap whose ordering was
never meant to be the primary filter.

**So: does the cap matter?** With a rarity/coding filter in front of pairing,
barely — it fires on a handful of hyper-variable genes and that is correct
behaviour. As the pipeline is wired today, it matters enormously, because it is
doing the job a filter should be doing, on 2 million variants, and announcing it
on every single result.

---

## 6. Gene-interval index (`mva.annotation.gene_intervals`)

The module landed during this work and imports cleanly. Measured against the
real MANE GRCh38 v1.5 release in `mva-resources/mane/`:

```
ManeGeneIndex built     : 0.43 s, 19,299 genes, 80 MB RSS
genes_at point queries  : 1,443,933 lookups/s
projected for 4.5M      : 3.1 s
hits                    : 0.449 genes per variant
```

**This one is not a problem and will not become one.** 3 seconds and 80 MB for
the whole genome. The binary search its docstring promises is real. (The
harness computes the files' SHA-256 to satisfy the constructor's fail-closed
integrity pin; that is enough for a throughput measurement and is explicitly not
an integrity check.)

---

## 7. Verdict

**Is this pipeline ready to ingest the real VCF? No — on memory alone.**

Nothing is quadratic. Nothing is unexpectedly slow. Read + normalise + QC is a
2–3 minute job at rates measured on data that fits. The single blocking defect is
that the whole callset is materialised as Pydantic models, three times over.

### What must change first, in order

1. **Stop materialising the callset. `read_vcf` must stream.**
   The blocker is `records.sort(key=...)` at the end — the only reason every
   record must exist at once. A caller's VCF is already coordinate-sorted, so
   that sort exists to impose karyotype contig order across the file, which a
   tabix-driven per-contig read in `CANONICAL_CONTIGS` order gives for free.
   Add `iter_vcf(...) -> Iterator[VariantRecord]` beside `read_vcf` and keep
   `read_vcf` as `tuple(iter_vcf(...))` for the fixture-scale callers and the
   tests. **This is the change that decides whether the real file runs at all.**
   Expected effect: ~20 GB → bounded by chunk size.

2. **Do not hold three copies.** `orchestrator.py` keeps `ingestion.variants`,
   `normalised.variants` and `qc.variants` alive simultaneously — measured at
   7.13 KB/record against read's 4.20 KB. Fuse normalise and QC into the
   streaming pass, or spill to Parquet/DuckDB, both already dependencies.

3. **Switch the default backend, or fix the adapter.** `detect_backend()`
   prefers cyvcf2, which measured **3.2x slower** on the full file. The
   immediate fix is one line at the call site (`backend="text"`); the real fix
   is to stop calling `variant.format()` four times per record and to stop
   rebuilding a GT string htslib had already decoded.

4. **Do not build a 2.6 GB string to write a file.** `canonical_json` must
   stream, or the variants artifacts must move to Parquet.

5. **Put a rarity/coding filter in front of pairing** — as a *selection* step
   that reports what it set aside, not as a hard filter (GP-13/GP-14 both
   survive that). Today the cap silently does this job on 2 M variants and
   stamps `gene_pair_cap_truncated` on 100% of output, which makes the flag
   meaningless. With ~460 candidates in front of it, the cap fires only on the
   handful of hyper-variable genes it was designed for.

Items 1 and 2 are architectural and block the real file. Items 3–5 are
contained changes worth making anyway.

### Cheap wins visible in the profile, not blocking

- `normalise_contig` runs 3x per record, each a regex match. Two of the three
  re-derive a value already known. A dict lookup removes ~500 k regex matches
  per 100 k records.
- 2.5 M `str.strip` calls per 100 k records — 25 per record, most on fields
  already known to be clean.

### What was not measured, and why

- **`normalise_variants` and `assess_quality` at full 4.5 M scale.** Read-only
  was run at full scale; adding ~13 GB of further demand to a machine already
  10 GB oversubscribed would have measured the swap subsystem, not the code.
  Their cost is projected from the clean 100k/500k/1M points, which are linear
  to within 2%, and is labelled as a projection throughout. The 1M point was
  re-measured against the current tree after `normalise.py` and `variant.py`
  changed mid-run.
- **End-to-end `mva run all`.** Out of scope by instruction — `annotation/`,
  `prioritization/`, `privacy/`, `cli.py`, `orchestrator.py` and `pipeline.py`
  were being edited concurrently while this ran.
- **Annotation throughput.** Same reason. Note that it is the stage immediately
  after the one that already exhausts memory.

---

## 8. After: the streaming ingest, measured

**Appended, not substituted.** Everything above §7 is the pre-change measurement
and stands as written. This section reports what the same harness measures
against the same phantom after `reader.py` and `qc.py` were rewritten to stream.

Reproduce the batch rows with `tools/scale/harness.py sweep` as before. The
streaming rows come from a driver that imports `peak_rss_bytes` from that same
harness — same `resource.getrusage(RUSAGE_SELF).ru_maxrss` method, same
one-child-per-configuration discipline — and runs the composition described in
§8.3. Its source is reproduced in §8.8, so every number here is re-runnable.

### 8.1 The headline

`iter_vcf` + chunked `normalise_variants` + `iter_assessed`, text backend, the
full 4,492,805-line phantom (4,556,006 records out):

| stage | before (§1, §2) | after | change |
|---|---|---|---|
| `read` peak RSS | **~19.6 GB** (projected; 6.55 GB measured under paging) | **59.0 MB** measured | **332x less** |
| through `assess_quality`, peak RSS | **~32.5 GB** (projected) | **261.9 MB** measured | **124x less** |
| through QC **+ the 2.41 GB variants artifact written** | not reachable | **242.9 MB** measured | — |
| `read` wall | 142.35 s / 32,005 rec/s | 37.35 s / 121,997 rec/s | 3.8x faster |
| through QC, wall | not run at full scale ‡ | 93.83 s / 48,557 rec/s | — |
| swap growth | 3.1 GB -> 33.0 GB | **0 swaps, 0 page faults** (`/usr/bin/time -l`) | — |

**The 4 GB target is met with 15x headroom: 261.9 MB through QC.** Record
count, warnings and skips are unchanged at full scale —
`no_call_genotype (n=13553)`, 0 skipped — matching §1 exactly.

The artifact row is the one that closes §4: `canonical_json` per *record* into an
open file handle produced a 2,409,978,805-byte artifact while resident memory
never exceeded 243 MB. The 2.6 GB single Python `str` is not needed; it was only
ever an artifact of building the whole payload before writing any of it.

> **Wall-clock caveat, stated because it matters.** The 37.35 s read was measured
> on a quiet machine. The later runs in §8.2 and §8.4 were taken while a
> concurrent agent was downloading reference data (load average 1.9, Spotlight
> indexing, sustained disk I/O), and the *identical* read took 73.3 s — 1.96x
> slower for the same work. Every throughput comparison below is therefore
> **load-matched**: before and after measured back to back, minutes apart, on the
> same machine state. The memory figures need no such caveat; they reproduced to
> within 4% across every run.

### 8.2 Load-matched before/after at 1 M, both paths

Both halves measured within minutes of each other under the same background load.
The "before" column runs the *original* `reader.py` and `qc.py`, restored into the
tree as `_reader_before` / `_qc_before` for the measurement and removed after.

| configuration | before | after | change |
|---|---|---|---|
| `read`, text | 24.64 s / 41,042 rec/s / **4.373 GB** | 16.04 s / 63,037 rec/s / **0.059 GB** | 1.54x faster, **74x less memory** |
| `read`, cyvcf2 | 32.52 s / 31,092 rec/s / 4.39 GB | (see §8.4) | — |
| read+normalise+QC, text | 57.06 s / **7.209 GB** | 38.92 s / **0.132 GB** | 1.47x faster, **55x less memory** |

The memory ratio grows with the file: it is bounded work against linear work.

### 8.3 What the streaming path actually is

Three changes, each of which had to be provably equivalent to what it replaced.

1. **`read_vcf` no longer holds the file twice.** `_parse_records` used to return a
   complete `tuple[_RawRecord, ...]` that stayed alive for the whole
   record-building loop. Both parsers are now generators. `read_vcf` itself still
   materialises and still sorts globally — its contract is unchanged and every
   existing test passes untouched — but its peak fell from **4.37 GB to 3.92 GB
   per million records (-10%)**, which is the second copy going away.

2. **`iter_vcf` removes the global sort rather than optimising it.** Contigs are
   visited in `CANONICAL_CONTIGS` order through the tabix index; within a contig
   the file is already position-sorted; the only thing buffered is the group of
   records sharing one coordinate, ordered on `(ref, alt)` exactly as the global
   sort's key would. Measured chunk sizes on the 1 M file: 247 chunks, min 2,963,
   max 4,112 records — a **15 MB** upper bound on what is ever alive.

3. **`iter_assessed` streams QC.** Counts are running integers; the two statistics
   that genuinely need every value (median depth, and the `math.fsum`-exact means
   behind `statistics.fmean`) go into an `array('q')` and an `array('d')` — 16
   bytes per record, ~73 MB at full scale — and are handed to the *same*
   `statistics` functions, so the metrics block is bit-identical rather than
   approximately equal.

`normalise_variants` was not modified. It is called on chunks that `VcfStream.chunks`
cuts only where the next record is more than 4,096 bases away or on a contig
change, which clears `MAX_SHIFT_BP = 1000` — so concatenating independently
sorted chunk outputs is provably the same sequence as sorting everything once.
`tests/unit/test_streaming.py::test_chunked_normalisation_equals_whole_callset_normalisation`
asserts it, and `::test_the_chunk_boundary_gap_clears_the_normalisation_shift_window`
asserts the constant relationship the argument rests on.

### 8.4 Backend: re-decided on evidence, after fixing the adapter

§1 found cyvcf2 3.19x slower and blamed the adapter, so the adapter was fixed
before the decision was re-taken: `variant.FORMAT` is consulted once per record
and absent tags are never requested (the old code called `variant.format("PS")`
on every record and caught the resulting exception), each array is converted with
one `tolist()` instead of an element-by-element `int()` loop, and the GT string is
rebuilt once per *distinct* call rather than once per record.

Body parse, 100 k phantom, old and new alternating in one process, best of four:

```
old cyvcf2 adapter   1.532s    64,864 rec/s
new cyvcf2 adapter   1.208s    82,272 rec/s      -21%
old text parser      0.612s   162,303 rec/s
new text parser      0.553s   179,732 rec/s
```

Full file, streaming, back to back, twice each:

| backend | wall | rec/s | peak RSS | records | warnings | skipped |
|---|---|---|---|---|---|---|
| text | 73.25 s / 73.48 s | 62,200 / 62,004 | 59.1 / 59.4 MB | 4,556,006 | `no_call_genotype (n=13553)` | 0 |
| cyvcf2 | 109.04 s / 108.80 s | 41,783 / 41,874 | 81.0 / 77.4 MB | 4,556,006 | `no_call_genotype (n=13553)` | 0 |

**The adapter fix was worth 21% and did not change the answer.** text stays
1.49x faster (down from 3.19x), so `detect_backend()` now returns `"text"` and
says why, with these numbers, in its docstring. `read_vcf(backend="cyvcf2")` is
unchanged and still supported; `IngestionResult.backend` and
`IngestionSummary.backend` now record which parser ran, because the two do *not*
agree on malformed input — a line with no FORMAT/sample column is a warned-about
record to the text parser and an unreadable file to htslib, which is a fact a run
must be able to state rather than a detail to discover later.

### 8.5 Determinism (GP-30), proved not asserted — under a fixed clock

The proof below is real and reproducible, and its scope is the ingestion stream
under the clock `config.synthetic` selects, which is `demo_clock()`. It is not
evidence that a real run's artifact tree is byte-identical: that run takes
`SystemClock`, and 11 of its artifacts differ in recorded time (none in scientific
content). See `docs/handoff-integrity.md` §4 and TD-21.

The full 4,556,006-record stream — read, normalise, QC — was digested end to end
(sha256 over `canonical_json` of every record *and* every evidence item, in
emission order) in two processes under different `PYTHONHASHSEED`:

```
PYTHONHASHSEED=0            stream 6ef1aae2dd7827bc4777561c7aa343a36b2affd16c2aad6081cd02b5fb410500
PYTHONHASHSEED=987654321    stream 6ef1aae2dd7827bc4777561c7aa343a36b2affd16c2aad6081cd02b5fb410500
                           metrics 00cd7accb73a1703b06b9decfcfd26202aab24220a93b4e22669eee05162033b  (both)
```

Byte-identical. The contig visit order is written out from `CANONICAL_CONTIGS` in
`_contig_visit_plan`, never read back from the index, a `set` or a `dict`; the
same proof runs on the fixture in
`test_repeat_streams_are_byte_identical_under_different_hash_seeds`.

Where the file's own contig blocks are *not* in karyotype order, streaming
declines and falls back to the buffered sort. That is not about output order —
the plan would order the records correctly regardless — it is about
`source_line_index`, which is an ordinal over the file's data lines and can only
be reproduced by counting them in the order the file stores them. Emitting a
provenance field that quietly means something different from the one `read_vcf`
emits would be worse than not streaming.

### 8.6 QC metrics at full scale, for the record

Identical under both hash seeds; the first whole-genome QC summary this pipeline
has produced:

```
total_variants 4,556,006   flagged 399,149   unflagged 4,156,857
mean_depth 31.9   median_depth 31.0   min 3   max 118
mean_het_allele_fraction 0.499792   (2,654,126 het calls with a fraction)
filtered_by_caller 364,984   possible_mosaic 13,632   high_allele_balance 13,690
low_depth 10,325   low_gq 872   low_allele_balance 388
evidence: supports 4,156,857   contradicts 386,980   neutral 12,169
```

### 8.7 What is still not fixed

- **`orchestrator.py` still holds three collections at once** and still calls
  `read_vcf`, `ledger.extend(qc.evidence)` and `write_json_artifact` with a
  materialised list. Every measurement in §8.1 was taken through the streaming
  API directly; the orchestrator has to adopt it before the real file will run.
  The exact diff is in the handover, not applied here — `orchestrator.py` is
  owned by another agent.
- **`canonical_json` still builds one string.** §8.1's artifact row was written by
  calling it per record into an open handle. `PipelineContext.write_json_artifact`
  needs a streaming sibling; the diff is in the handover.
- **4.5 M `EvidenceItem`s cannot go into `EvidenceLedger`.** `iter_assessed` emits
  them one at a time, but the ledger's `extend` still wants them all. At WGS scale
  the ledger needs to spill or to hold aggregates.
- **§5's finding is untouched.** No rarity/coding selection was added in front of
  pairing; `gene_pair_cap_truncated` would still fire on 100% of output.

### 8.8 The streaming driver

`tools/scale/harness.py` measures `read_vcf`, which is the batch path by
definition. The streaming rows above were produced by this driver, which reuses
the harness's own RSS functions rather than re-implementing them. It is the
composition the orchestrator is being asked to adopt, so it doubles as the
worked example.

```python
from mva.clock import demo_clock
from mva.config import QualityThresholds
from mva.determinism import canonical_json
from mva.ingestion.normalise import normalise_variants
from mva.ingestion.qc import iter_assessed
from mva.ingestion.reader import iter_vcf
from mva.models.genome import GenomeBuild
from tools.scale.harness import current_rss_bytes, peak_rss_bytes

stream = iter_vcf(vcf, expected_build=GenomeBuild.GRCH38,
                  source_artifact="scale-phantom", backend=backend)

normalised = (record
              for chunk in stream.chunks()
              for record in normalise_variants(chunk).variants)

assessed = iter_assessed(normalised, thresholds=QualityThresholds(), clock=demo_clock())

with out_path.open("w", encoding="utf-8") as handle:          # the §4 fix
    handle.write("[")
    for n, item in enumerate(assessed):
        if n:
            handle.write(",")
        handle.write(canonical_json(item.variant.model_dump(mode="json")))
    handle.write("]\n")

summary = stream.summary()      # refuses to report unless the stream was drained
metrics = assessed.metrics()
peak = peak_rss_bytes()
```

Run as one child process per configuration
(`--stages read` / `read,normalise,qc` / `read,normalise,qc,artifact`), because
`ru_maxrss` is a high-water mark and a shared process reports one run's peak as
the next one's.

---

## 9. After: annotation, the evidence ledger and selection, measured

**Appended, not substituted.** Sections 1-8 are the ingest measurement and stand
as written. This section covers three things the ingest measurement left open: annotation
materialising the callset (§7 lists its throughput as unmeasured; it is not one
of §8.7's four bullets), the ledger unable to hold the evidence it is handed
(§8.7), and no selection in front of pairing (§8.7).

Reproduce with `tools/scale/stage_harness.py` (`sweep`, `annotate`, `ledger`,
`select`, `pipeline`, `digest`). It measures against fabricated `VariantRecord`
objects built in memory rather than parsed from a VCF, because the question here
is what one annotated record plus its evidence costs and whether a stage holds
all of them at once; a parser in front of that adds noise and 40 seconds per run.
Peak RSS is `ru_maxrss`, taken in a child process per configuration, exactly as
in §1. A watchdog thread kills a child at 12 GiB, so a configuration that cannot
fit is reported as killed rather than allowed to take the machine into swap —
that is a measurement, not a failure to measure.

The phantom's rates are matched to the ones measured earlier on real data: 44.9%
of records genic (§6), 38% with no frequency data at all, 1.2% of genic records
coding or splice-relevant, 8% carrying an impact of `None` — NOT ASSESSED, the
shape a MANE interval join produces (ADR 0016).

> **Not biologically valid.** Every coordinate, allele, gene symbol, consequence
> term and allele frequency in the phantom is a seeded hash. Nothing in this
> section is evidence about variants, genes or disease (GP-20).

**Load caveat, stated because it matters.** These runs shared a 24 GB machine
with several concurrent agents; load average sat between 6 and 9 throughout, and
one full `pytest` run overlapped part of the sweep. **Wall times are therefore an
upper bound** and are not comparable with §1's, which were taken on a quieter
machine. Peak RSS is unaffected by CPU contention and is what this section is
about.

### 9.1 The sweep

Each row is one child process. `KILLED` means the 12 GiB watchdog fired.

| stage | mode | records | wall | peak RSS | |
|---|---|---:|---:|---:|---|
| annotate | batch | 100,000 | 7.4s | 1.046 GiB | ok |
| annotate | stream | 100,000 | 4.6s | 0.080 GiB | ok |
| ledger | memory | 100,000 | 7.3s | 0.652 GiB | ok |
| ledger | spill | 100,000 | 20.9s | 0.337 GiB | ok |
| select | enabled | 100,000 | 4.2s | 0.088 GiB | ok |
| pipeline | stream | 100,000 | 16.4s | 0.411 GiB | ok |
| annotate | batch | 500,000 | 28.6s | 5.592 GiB | ok |
| annotate | stream | 500,000 | 22.8s | 0.080 GiB | ok |
| ledger | memory | 500,000 | 43.1s | 3.124 GiB | ok |
| ledger | spill | 500,000 | 87.3s | 0.351 GiB | ok |
| select | enabled | 500,000 | 24.6s | 0.125 GiB | ok |
| pipeline | stream | 500,000 | 86.1s | 0.413 GiB | ok |
| annotate | batch | 1,000,000 | 94.0s | 7.787 GiB | ok |
| annotate | stream | 1,000,000 | 44.9s | 0.080 GiB | ok |
| ledger | memory | 1,000,000 | 75.3s | 5.807 GiB | ok |
| ledger | spill | 1,000,000 | 246.8s | 0.353 GiB | ok |
| select | enabled | 1,000,000 | 112.9s | 0.165 GiB | ok |
| pipeline | stream | 1,000,000 | 248.8s | 0.359 GiB | ok |
| annotate | batch | 4,500,000 | — | — | **KILLED** |

**What the sweep settles.** Batch annotation is the stage that cannot reach WGS
scale: peak memory grows with input across the measured range (1.05 → 5.59 →
7.79 GiB) and the 4.5 M configuration hit the watchdog. The growth is **not**
proportional — per record it is 11.0, then 11.7, then 8.2 KiB, so the 1 M run
cost less per record than the 500 K one. It does not flatten either: holding the
best of those three rates, 4.5 M is `projected` to need roughly 35 GiB, against a
12 GiB watchdog and 24 GB of machine. Streaming annotation is **flat at 0.080 GiB
across a 10x range of input** — it holds a batch, not a callset. The in-memory
ledger has the same defect as batch annotation (0.65 → 3.12 → 5.81 GiB); the
spilled ledger is likewise flat, at ~0.35 GiB, and buys that with wall time
(20.9s → 246.8s, roughly 3-4x the in-memory ledger).

That trade is the intended one: **the spilled ledger is slower and finishes;
the in-memory ledger is faster and dies.**

### 9.2 At 4.5 M

Two stages were then run at the full target scale.

| stage | mode | records | wall | throughput | peak RSS |
|---|---|---:|---:|---:|---:|
| annotate | stream | 4,500,000 | 316.1s | 14,238 rec/s | 0.080 GiB |
| select | enabled | 4,500,000 | 363.0s | 12,397 rec/s | 0.476 GiB |

Streaming annotation held 0.080 GiB at 4.5 M, the **same** figure as at 100,000:
peak memory is now a function of batch size, not of callset size.

Selection at full scale, with the thresholds in `config/default.yaml`:

| | count | share |
|---|---:|---:|
| in | 4,500,000 | |
| **retained** | **80,250** | **1.78%** |
| dropped, common in population | 2,733,944 | 60.8% |
| dropped, no gene assignment | 999,764 | 22.2% |
| dropped, not coding or splice | 686,042 | 15.2% |
| dropped, frequency unknown | 0 | 0% |
| dropped, impact NOT ASSESSED | 0 | 0% |

Retained breaks down as 63,867 NOT-ASSESSED-impact, 15,513 coding-or-splice, and
870 on a pathogenic clinical assertion.

**The two zeroes are the point.** `dropped_frequency_unknown = 0` and
`dropped_impact_not_assessed = 0` are GP-14 and ADR 0016 holding under load:
**1,709,554 variants reached selection with no usable frequency observation and
every one was retained.** Absence of a frequency record is not evidence of
rarity, so it cannot be scored as AF 0 and cannot be used to drop. Had those
been treated as common, selection would have discarded 38% of the callset on
missing data alone.

Selection emits both facts as warnings on its own report rather than leaving
them to be discovered — including that it deleted 4,419,750 valid records, which
it is the only stage entitled to do (GP-13, ADR 0019).

### 9.3 What was not measured

**The spilled ledger was never run at 4.5 M.** The run was started and stopped
before it completed; no number was recorded and none is estimated here. What
stands is the flat ~0.35 GiB across 100 K–1 M in §9.1, which is evidence that it
does not grow with input, but it is **not** a 4.5 M measurement and must not be
cited as one.

Also still unmeasured: the real (non-phantom) adapters at scale, and the pairing
stage downstream of selection.
