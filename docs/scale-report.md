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
