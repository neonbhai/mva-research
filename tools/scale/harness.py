"""Scale harness: measure ``mva.ingestion`` against a WGS-scale VCF.

Everything reported by this module is observed, not modelled. Where a number is
extrapolated it is labelled ``projected_*`` and the measurement it was projected
from is reported next to it.

Peak memory is ``resource.getrusage(RUSAGE_SELF).ru_maxrss`` (bytes on macOS,
kibibytes on Linux; normalised here). It is measured in a **child process per
configuration**, because a peak is a high-water mark: once one backend has
allocated 8 GB, every later measurement in the same process inherits that mark
and reports it as its own.

Stages run in the order the orchestrator runs them, and — deliberately — hold
their results the way the orchestrator holds them: ``read_vcf`` result, then
``normalise_variants`` result, then ``assess_quality`` result, all alive at
once. That is the shape of the real memory question, and measuring the stages
in isolation would hide it.

Subcommands
-----------
``bench``     one configuration, in this process; prints a JSON result.
``sweep``     a matrix of ``bench`` child processes; prints a table + JSON.
``profile``   ``cProfile`` over one read, to find the bottleneck rather than guess.
``pairs``     the post-filter candidate and same-gene pair arithmetic.
``modelbench`` build time and resident cost of one ``VariantRecord``, in isolation.
``pairbench`` measured scaling of ``generate_pairs`` against one growing gene.
``intervals`` gene-interval lookup throughput, if that module exists yet.
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import resource
import subprocess
import sys
import time
from collections.abc import Callable, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Final

REPO_ROOT: Final = Path(__file__).resolve().parents[2]

#: macOS reports ru_maxrss in bytes, Linux in kibibytes.
_MAXRSS_SCALE: Final = 1 if sys.platform == "darwin" else 1024


def peak_rss_bytes() -> int:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * _MAXRSS_SCALE


def current_rss_bytes() -> int:
    """Current RSS via psutil; falls back to the peak when psutil is absent."""
    try:
        import psutil  # noqa: PLC0415 - optional, and only used for reporting
    except ImportError:  # pragma: no cover - psutil is installed here
        return peak_rss_bytes()
    return int(psutil.Process(os.getpid()).memory_info().rss)


# ---------------------------------------------------------------------------
# Result record
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class StageTiming:
    stage: str
    seconds: float
    rows_out: int
    rss_after: int
    peak_rss: int


@dataclass(slots=True)
class BenchResult:
    vcf: str
    vcf_bytes: int
    backend: str
    baseline_rss: int
    stages: list[StageTiming] = field(default_factory=list)
    notes: dict[str, Any] = field(default_factory=dict)

    def as_json(self) -> str:
        payload = asdict(self)
        return json.dumps(payload, indent=2, sort_keys=True)


# ---------------------------------------------------------------------------
# bench: one configuration, measured in this process
# ---------------------------------------------------------------------------


def _timed[T](label: str, fn: Callable[[], T], rows: Callable[[T], int]) -> tuple[T, StageTiming]:
    gc.collect()
    start = time.perf_counter()
    value = fn()
    elapsed = time.perf_counter() - start
    return value, StageTiming(
        stage=label,
        seconds=elapsed,
        rows_out=rows(value),
        rss_after=current_rss_bytes(),
        peak_rss=peak_rss_bytes(),
    )


def run_bench(
    vcf: Path,
    *,
    backend: str,
    stages: Sequence[str],
    artifact_sample: int,
) -> BenchResult:
    """Read (and optionally normalise/QC/serialise) one VCF, measuring each step."""
    from mva.clock import demo_clock  # noqa: PLC0415 - keep import cost out of the baseline
    from mva.config import QualityThresholds  # noqa: PLC0415
    from mva.ingestion import assess_quality, normalise_variants, read_vcf  # noqa: PLC0415
    from mva.models.genome import GenomeBuild  # noqa: PLC0415

    gc.collect()
    result = BenchResult(
        vcf=str(vcf),
        vcf_bytes=vcf.stat().st_size,
        backend=backend,
        baseline_rss=current_rss_bytes(),
    )

    ingestion, timing = _timed(
        "read_vcf",
        lambda: read_vcf(
            vcf,
            expected_build=GenomeBuild.GRCH38,
            source_artifact="scale-phantom",
            backend=backend,
        ),
        lambda r: len(r.variants),
    )
    result.stages.append(timing)
    result.notes["skipped_count"] = ingestion.skipped_count
    result.notes["skipped_reasons"] = list(ingestion.skipped_reasons)
    result.notes["reader_warnings"] = list(ingestion.warnings)

    normalised = None
    if "normalise" in stages:
        normalised, timing = _timed(
            "normalise_variants",
            lambda: normalise_variants(ingestion.variants),
            lambda r: len(r.variants),
        )
        result.stages.append(timing)
        result.notes["normalisation_operations"] = normalised.operations_applied

    if "qc" in stages and normalised is not None:
        qc, timing = _timed(
            "assess_quality",
            lambda: assess_quality(
                normalised.variants, thresholds=QualityThresholds(), clock=demo_clock()
            ),
            lambda r: len(r.variants),
        )
        result.stages.append(timing)
        result.notes["qc_evidence_items"] = len(qc.evidence)

    if "artifact" in stages:
        source = (normalised.variants if normalised is not None else ingestion.variants)[
            :artifact_sample
        ]
        blob, timing = _timed(
            "artifact_json_sample",
            lambda: json.dumps([v.model_dump(mode="json") for v in source]),
            len,
        )
        timing.rows_out = len(source)
        result.stages.append(timing)
        result.notes["artifact_sample_records"] = len(source)
        result.notes["artifact_sample_bytes"] = len(blob)
        result.notes["artifact_bytes_per_record"] = len(blob) / max(len(source), 1)
        del blob

    result.notes["final_peak_rss"] = peak_rss_bytes()
    return result


# ---------------------------------------------------------------------------
# sweep: a matrix of bench children
# ---------------------------------------------------------------------------


def _bench_child(
    vcf: Path, backend: str, stages: Sequence[str], artifact_sample: int, timeout: float
) -> dict[str, Any]:
    """Run one ``bench`` in a fresh interpreter so its peak RSS is its own."""
    cmd = [
        sys.executable,
        "-m",
        "tools.scale.harness",
        "bench",
        "--vcf",
        str(vcf),
        "--backend",
        backend,
        "--stages",
        ",".join(stages),
        "--artifact-sample",
        str(artifact_sample),
    ]
    started = time.perf_counter()
    proc = subprocess.run(  # noqa: S603 - fixed argv, no shell, paths from the caller
        cmd,
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
    )
    wall = time.perf_counter() - started
    if proc.returncode != 0:
        return {
            "vcf": str(vcf),
            "backend": backend,
            "error": proc.stderr.strip().splitlines()[-1] if proc.stderr.strip() else "failed",
            "returncode": proc.returncode,
            "child_wall_seconds": wall,
        }
    payload: dict[str, Any] = json.loads(proc.stdout)
    payload["child_wall_seconds"] = wall
    return payload


def _row(payload: dict[str, Any]) -> str:
    if "error" in payload:
        return f"{Path(payload['vcf']).name:<28} {payload['backend']:<8} ERROR {payload['error']}"
    stages: list[dict[str, Any]] = payload["stages"]
    read = next(s for s in stages if s["stage"] == "read_vcf")
    rate = read["rows_out"] / read["seconds"] if read["seconds"] else 0.0
    parts = [
        f"{Path(payload['vcf']).name:<28}",
        f"{payload['backend']:<8}",
        f"{read['rows_out']:>9,}",
        f"{read['seconds']:>8.2f}s",
        f"{rate:>10,.0f}/s",
        f"{read['peak_rss'] / 1e9:>7.2f}GB",
    ]
    for name in ("normalise_variants", "assess_quality"):
        stage = next((s for s in stages if s["stage"] == name), None)
        parts.append(f"{stage['seconds']:>8.2f}s" if stage else f"{'-':>9}")
    parts.append(f"{payload['notes']['final_peak_rss'] / 1e9:>7.2f}GB")
    return " ".join(parts)


def run_sweep(
    vcfs: Sequence[Path],
    backends: Sequence[str],
    stages: Sequence[str],
    *,
    artifact_sample: int,
    timeout: float,
    out: Path | None,
) -> int:
    header = (
        f"{'file':<28} {'backend':<8} {'records':>9} {'read':>9} "
        f"{'rec/s':>11} {'read_peak':>9} {'normalise':>9} {'qc':>9} {'pipe_peak':>9}"
    )
    print(header)  # noqa: T201 -- this tool's CLI report
    print("-" * len(header))  # noqa: T201 -- this tool's CLI report
    payloads: list[dict[str, Any]] = []
    for vcf in vcfs:
        for backend in backends:
            payload = _bench_child(vcf, backend, stages, artifact_sample, timeout)
            payloads.append(payload)
            print(_row(payload))  # noqa: T201 -- this tool's CLI report
            sys.stdout.flush()
    if out is not None:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(payloads, indent=2, sort_keys=True), encoding="utf-8")
        print(f"\nwrote {out}")  # noqa: T201 -- this tool's CLI report
    return 0


# ---------------------------------------------------------------------------
# profile
# ---------------------------------------------------------------------------


def run_profile(vcf: Path, backend: str, *, top: int, sort: str) -> int:
    """cProfile one ``read_vcf``. Use a subsample; the profiler costs ~2x."""
    import cProfile  # noqa: PLC0415 - profiling-only import
    import pstats  # noqa: PLC0415

    from mva.ingestion import read_vcf  # noqa: PLC0415
    from mva.models.genome import GenomeBuild  # noqa: PLC0415

    profiler = cProfile.Profile()
    profiler.enable()
    result = read_vcf(
        vcf, expected_build=GenomeBuild.GRCH38, source_artifact="scale-phantom", backend=backend
    )
    profiler.disable()
    print(  # noqa: T201 -- this tool's CLI report
        f"# profile: {vcf.name} backend={backend} records={len(result.variants):,}"
    )
    stats = pstats.Stats(profiler, stream=sys.stdout)
    stats.sort_stats(sort).print_stats(top)
    return 0


# ---------------------------------------------------------------------------
# pairs: post-filter candidate and pair arithmetic
# ---------------------------------------------------------------------------

#: Attrition of a **standard clinical rare-disease cascade** on a single-sample
#: WGS callset. Every entry is a rate so a reader who disagrees can substitute
#: their own and re-run. These are the routine ballparks of diagnostic WGS
#: filtering practice, anchored on the two figures that are stable across
#: studies -- roughly 20-24k coding variants per genome, of which a few hundred
#: survive a 1% rarity cut -- NOT measurements of this pipeline or of the
#: phantom. They are labelled as assumptions throughout the report.
#:
#: NOTE: this cascade is NOT what ``mva.prioritization.filters.apply_hard_filters``
#: does. That function deliberately removes only invalid/hom-ref/no-call records
#: (GP-13: filtering is not ranking). The gap between the two is measured by the
#: ``as_implemented`` scenario below and is the number that actually matters.
CASCADE: Final[tuple[tuple[str, float, str], ...]] = (
    ("pass_filter", 0.92, "FILTER=PASS"),
    (
        "coding_or_splice",
        0.0053,
        "in CDS or within a splice region; anchors on ~22-24k coding+splice "
        "variants in a 4.5M-variant genome",
    ),
    (
        "not_common",
        0.035,
        "gnomAD popmax AF <= 1%; most coding variation in any one genome is common",
    ),
    (
        "functional_class",
        0.60,
        "missense / LoF / splice rather than synonymous or UTR",
    ),
)

#: Genes that reliably accumulate many rare coding variants in any single
#: genome -- large, repetitive, or in a hyper-polymorphic locus. The uniform
#: Poisson model says the per-gene cap never binds; these genes are why that
#: conclusion is wrong, and the count is what the cap is actually for.
HYPERVARIABLE_GENE_LOAD: Final[tuple[tuple[str, int], ...]] = (
    ("TTN", 14),
    ("MUC16", 10),
    ("OBSCN", 8),
    ("NEB", 7),
    ("RYR1", 6),
    ("HLA-DRB1", 6),
    ("HLA-B", 5),
    ("AHNAK2", 5),
    ("PLEC", 4),
    ("SYNE1", 4),
)


def _candidates_built(k: int, max_pairing_variants: int) -> int:
    """Objects ``generate_pair_candidates`` builds for a gene holding ``k`` variants.

    Mirrors the implementation as it stands: a single candidate for every member
    that warrants one (worst case, all of them), plus ``C(min(k, P), 2)`` pairs,
    where ``P`` is ``max_pairing_variants``. That second bound is what makes the
    function linear rather than quadratic in ``k`` -- before it existed the term
    was ``C(k, 2)``.
    """
    paired = min(k, max_pairing_variants)
    return k + paired * (paired - 1) // 2


def _poisson_pair_load(
    genes: int, mean_per_gene: float, max_pairs_per_gene: int, max_pairing_variants: int
) -> tuple[float, float, float]:
    """Exact (unsampled, reproducible) Poisson expectation of the enumeration."""
    import math  # noqa: PLC0415 - local to this calculation

    before = capped = over_cap = 0.0
    log_mean = math.log(mean_per_gene) if mean_per_gene > 0 else float("-inf")
    horizon = max(400, int(mean_per_gene * 6) + 400)
    for k in range(horizon):
        if mean_per_gene <= 0:
            prob = 1.0 if k == 0 else 0.0
        else:
            prob = math.exp(-mean_per_gene + k * log_mean - math.lgamma(k + 1))
        expected_genes = genes * prob
        total = _candidates_built(k, max_pairing_variants)
        before += expected_genes * total
        capped += expected_genes * min(total, max_pairs_per_gene)
        if total > max_pairs_per_gene or k > max_pairing_variants:
            over_cap += expected_genes
    return before, capped, over_cap


def _report_pair_load(
    label: str, genes: int, candidates: int, cap: int, max_pairing_variants: int
) -> None:
    mean = candidates / genes if genes else 0.0
    before, capped, over = _poisson_pair_load(genes, mean, cap, max_pairing_variants)
    print(  # noqa: T201 -- this tool's CLI report
        f"\n## {label}\n"
        f"    candidate variants        : {candidates:>15,}\n"
        f"    genes assumed             : {genes:>15,}\n"
        f"    mean candidates per gene  : {mean:>15.3f}\n"
        f"    objects built             : {before:>15,.0f}   <- built BEFORE the cap\n"
        f"    kept after cap={cap:<3d}       : {capped:>15,.0f}\n"
        f"    genes truncated           : {over:>15,.1f}\n"
        f"    hypotheses discarded      : {before - capped:>15,.0f}"
    )


def run_pairs(
    *,
    callset: int,
    genes: int,
    max_pairs_per_gene: int,
    max_pairing_variants: int,
    measured_pass_rate: float | None,
    genic_fraction: float,
) -> int:
    """Project post-cascade candidate counts and the load ``generate_pairs`` faces.

    Three scenarios, because the answer differs by four orders of magnitude
    between them and only one of them is what the pipeline does today.
    """
    surviving = float(callset)
    print(f"# callset: {callset:,} records")  # noqa: T201 -- this tool's CLI report
    for name, rate, why in CASCADE:
        effective = measured_pass_rate if (name == "pass_filter" and measured_pass_rate) else rate
        surviving *= effective
        print(f"  x {effective:<8.4f} {name:<20} -> {surviving:>12,.0f}   ({why})")  # noqa: T201

    candidates = round(surviving)
    _report_pair_load(
        "A. standard rare-disease cascade, uniform Poisson gene load",
        genes,
        candidates,
        max_pairs_per_gene,
        max_pairing_variants,
    )

    # The tail the Poisson model cannot see.
    tail_before = sum(
        _candidates_built(n, max_pairing_variants) for _, n in HYPERVARIABLE_GENE_LOAD
    )
    tail_capped = sum(
        min(_candidates_built(n, max_pairing_variants), max_pairs_per_gene)
        for _, n in HYPERVARIABLE_GENE_LOAD
    )
    tail_over = sum(
        1
        for _, n in HYPERVARIABLE_GENE_LOAD
        if _candidates_built(n, max_pairing_variants) > max_pairs_per_gene
    )
    print(  # noqa: T201 -- this tool's CLI report
        f"\n## B. heavy tail the Poisson model misses "
        f"({len(HYPERVARIABLE_GENE_LOAD)} hyper-variable genes)\n"
        f"    objects built             : {tail_before:>15,}\n"
        f"    kept after cap={max_pairs_per_gene:<3d}       : {tail_capped:>15,}\n"
        f"    genes truncated           : {tail_over:>15,}\n"
        f"    hypotheses discarded      : {tail_before - tail_capped:>15,}"
    )

    # What the pipeline actually feeds generate_pairs today: apply_hard_filters
    # removes only invalid/hom-ref/no-call (GP-13), so every alt-carrying
    # annotated variant with a gene symbol arrives.
    genic = round(callset * genic_fraction)
    _report_pair_load(
        f"C. AS IMPLEMENTED -- no rarity or coding hard filter "
        f"({genic_fraction:.1%} of the callset falls in a MANE gene, measured)",
        genes,
        genic,
        max_pairs_per_gene,
        max_pairing_variants,
    )
    return 0


# ---------------------------------------------------------------------------
# modelbench: what one VariantRecord costs to build and to hold
# ---------------------------------------------------------------------------


def run_modelbench(count: int) -> int:
    """Isolate ``VariantRecord`` construction from everything the reader does.

    ``read_vcf``'s memory could be its ``_RawRecord`` intermediates, its two
    container copies, or the models themselves. Building ``count`` bare records
    with nothing else alive attributes the cost.
    """
    from mva.models.genome import GenomeBuild, GenomicCoordinate  # noqa: PLC0415
    from mva.models.variant import FilterStatus, Genotype, VariantRecord, Zygosity  # noqa: PLC0415

    gc.collect()
    before = current_rss_bytes()
    start = time.perf_counter()
    records = [
        VariantRecord(
            coordinate=GenomicCoordinate(
                build=GenomeBuild.GRCH38,
                contig="chr1",
                position=1_000 + i,
                ref="A",
                alt="G",
            ),
            genotype=Genotype(
                zygosity=Zygosity.HET,
                genotype_string="0/1",
                depth=30,
                ref_reads=15,
                alt_reads=15,
                genotype_quality=99,
            ),
            filter_status=FilterStatus.PASS,
            source_artifact="scale-phantom",
            source_line_index=i,
        )
        for i in range(count)
    ]
    build_seconds = time.perf_counter() - start
    gc.collect()
    after = current_rss_bytes()

    start = time.perf_counter()
    dumps = [record.model_dump(mode="json") for record in records[: min(count, 200_000)]]
    dump_seconds = time.perf_counter() - start
    dumped = len(dumps)
    del dumps

    start = time.perf_counter()
    keys = [record.sort_key() for record in records]
    sort_key_seconds = time.perf_counter() - start
    del keys

    per_record = (after - before) / count
    print(  # noqa: T201 -- this tool's CLI report
        f"# VariantRecord x {count:,} (three nested pydantic models each)\n"
        f"    construction        : {build_seconds:>8.2f}s "
        f"= {count / build_seconds:>10,.0f} records/s\n"
        f"    resident cost       : {(after - before) / 1e9:>8.2f}GB "
        f"= {per_record:>10,.0f} bytes/record\n"
        f"    model_dump(json)    : {dump_seconds:>8.2f}s for {dumped:,} "
        f"= {dumped / dump_seconds:>10,.0f} records/s\n"
        f"    sort_key() x{count:<8,}: {sort_key_seconds:>8.2f}s "
        f"= {count / sort_key_seconds:>10,.0f} calls/s "
        f"(the reader, normalise and QC each call it once per record)\n"
        f"    projected for 4.5M  : {4_500_000 * per_record / 1e9:.1f}GB resident, "
        f"{4_500_000 / (count / build_seconds):.0f}s to construct"
    )
    del records
    return 0


# ---------------------------------------------------------------------------
# pairbench: is generate_pairs quadratic in practice, or only on paper?
# ---------------------------------------------------------------------------


def _synthetic_gene_variants(gene: str, count: int, genes: int = 1) -> list[Any]:
    """``count`` het, HIGH-impact records spread over ``genes`` gene symbols.

    Built in memory; nothing under ``src/`` is touched. HIGH impact is chosen
    because ``_wants_single_candidate`` keeps a lone het only at HIGH, so every
    member produces a single candidate as well as entering pair enumeration --
    the worst case for the enumeration cost.
    """
    from mva.models.genome import GenomeBuild, GenomicCoordinate  # noqa: PLC0415
    from mva.models.variant import (  # noqa: PLC0415
        ConsequenceAnnotation,
        FilterStatus,
        Genotype,
        ImpactSeverity,
        VariantRecord,
        Zygosity,
    )

    consequences = [
        (
            ConsequenceAnnotation(
                gene_symbol=gene if genes == 1 else f"{gene}{g:05d}",
                transcript_id="ENST00000000001",
                consequence_terms=("missense_variant",),
                impact=ImpactSeverity.HIGH,
            ),
        )
        for g in range(genes)
    ]
    genotype = Genotype(
        zygosity=Zygosity.HET,
        genotype_string="0/1",
        depth=30,
        ref_reads=15,
        alt_reads=15,
        genotype_quality=99,
    )
    records: list[Any] = []
    for i in range(count):
        records.append(
            VariantRecord(
                coordinate=GenomicCoordinate(
                    build=GenomeBuild.GRCH38,
                    contig="chr1",
                    position=1_000_000 + i * 7,
                    ref="A",
                    alt="G",
                ),
                genotype=genotype,
                filter_status=FilterStatus.PASS,
                source_artifact="scale-phantom",
                consequences=consequences[i % genes],
            )
        )
    return records


def run_pairbench(sizes: Sequence[int], max_pairs_per_gene: int, genes: int) -> int:
    """Time ``generate_pairs`` against a growing input.

    With ``genes=1`` the whole load lands in one gene, which is where a
    quadratic pair enumeration would show itself. With ``genes>1`` the load is
    spread, which is the shape a whole callset actually has.
    """
    try:
        from mva.prioritization.pairing import generate_pairs  # noqa: PLC0415
    except ImportError as exc:
        print(  # noqa: T201 -- this tool's CLI report
            f"mva.prioritization.pairing is not importable ({exc}); that package is "
            "being edited concurrently. generate_pairs was NOT measured."
        )
        return 0

    print(  # noqa: T201 -- this tool's CLI report
        f"{'variants_in':>13} {'genes':>8} {'kept':>10} {'seconds':>10} "
        f"{'us_per_object':>14} {'ratio_vs_prev':>14}"
    )
    previous: tuple[int, float] | None = None
    for n in sizes:
        variants = _synthetic_gene_variants("SYNTHGENE", n, genes=genes)
        gc.collect()
        start = time.perf_counter()
        kept = generate_pairs(variants, max_pairs_per_gene=max_pairs_per_gene)
        elapsed = time.perf_counter() - start
        enumerated = max(genes, 1) * _candidates_built(n // max(genes, 1), 24)
        ratio = ""
        if previous is not None and previous[1] > 0:
            observed = elapsed / previous[1]
            expected = (n / previous[0]) ** 2
            ratio = f"{observed:.2f}x (n^2={expected:.2f}x)"
        print(  # noqa: T201 -- this tool's CLI report
            f"{n:>13,} {genes:>8,} {len(kept):>10,} {elapsed:>10.3f} "
            f"{elapsed / max(enumerated, 1) * 1e6:>14.3f} {ratio:>14}"
        )
        previous = (n, elapsed)
        del variants
    return 0


# ---------------------------------------------------------------------------
# intervals: gene-interval index throughput, if the module exists
# ---------------------------------------------------------------------------


def _sha256(path: Path) -> str:
    """Digest of the file as it is on disk right now.

    ``ManeGeneIndex`` fails closed without an integrity pin. This harness has no
    manifest to pin against, so it digests the file it is about to read. That
    satisfies the constructor's contract for a *throughput* measurement and is
    explicitly NOT an integrity check -- it cannot detect a substituted file,
    because it would simply digest the substitute.
    """
    import hashlib  # noqa: PLC0415 - only needed for this measurement

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _query_positions(vcf: Path, limit: int) -> list[tuple[str, int]]:
    """(contig, position) pairs taken from the phantom, so locality is realistic."""
    import gzip  # noqa: PLC0415 - only needed for this measurement

    out: list[tuple[str, int]] = []
    opener = gzip.open if vcf.suffix in {".gz", ".bgz"} else open
    with opener(vcf, "rt", encoding="utf-8") as handle:  # type: ignore[operator]
        for line in handle:
            if line.startswith("#"):
                continue
            contig, position, _ = line.split("\t", 2)
            out.append((contig, int(position)))
            if len(out) >= limit:
                break
    return out


def run_intervals(lookups: int, mane_dir: Path, vcf: Path | None) -> int:
    """Measure ``ManeGeneIndex`` construction and ``genes_at`` throughput.

    Imported dynamically on purpose: the module is being written concurrently and
    may not exist or may not import. An absence is reported as an absence, never
    papered over with an invented number.
    """
    import importlib  # noqa: PLC0415 - existence of the target is what is being tested

    try:
        module = importlib.import_module("mva.annotation.gene_intervals")
    except Exception as exc:
        print(  # noqa: T201 -- this tool's CLI report
            f"mva.annotation.gene_intervals did NOT import ({type(exc).__name__}: {exc}). "
            "No lookup throughput was measured."
        )
        return 0

    index_cls = getattr(module, "ManeGeneIndex", None)
    if index_cls is None:
        names = sorted(n for n in dir(module) if not n.startswith("_"))
        print(  # noqa: T201 -- this tool's CLI report
            f"# imported, but no ManeGeneIndex; public names: {names}. Nothing measured."
        )
        return 0

    gtf = next(iter(sorted(mane_dir.glob("MANE.*ensembl_genomic.gtf.gz"))), None)
    summary = next(iter(sorted(mane_dir.glob("MANE.*summary.txt.gz"))), None)
    if gtf is None or summary is None:
        print(f"# no MANE release under {mane_dir}; nothing measured.")  # noqa: T201
        return 0

    start = time.perf_counter()
    index = index_cls(
        gtf,
        summary,
        expected_gtf_sha256=_sha256(gtf),
        expected_summary_sha256=_sha256(summary),
    )
    build_seconds = time.perf_counter() - start
    print(  # noqa: T201 -- this tool's CLI report
        f"# ManeGeneIndex({gtf.name}): built in {build_seconds:.2f}s, "
        f"{index.gene_count:,} genes, RSS now {current_rss_bytes() / 1e6:.0f} MB"
    )

    if vcf is None:
        print("# no --vcf given; genes_at throughput NOT measured.")  # noqa: T201
        return 0

    queries = _query_positions(vcf, lookups)
    gc.collect()
    hits = 0
    start = time.perf_counter()
    for contig, position in queries:
        hits += len(index.genes_at(contig, position))
    elapsed = time.perf_counter() - start
    rate = len(queries) / elapsed if elapsed else 0.0
    print(  # noqa: T201 -- this tool's CLI report
        f"# genes_at: {len(queries):,} point queries in {elapsed:.2f}s "
        f"= {rate:,.0f} lookups/s ({hits:,} gene hits, "
        f"{hits / max(len(queries), 1):.3f} genes/variant)\n"
        f"# projected for 4,500,000 variants: {4_500_000 / rate:.1f}s"
    )
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="harness", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    bench = sub.add_parser("bench", help="one configuration, in this process")
    bench.add_argument("--vcf", type=Path, required=True)
    bench.add_argument("--backend", default="cyvcf2")
    bench.add_argument("--stages", default="read,normalise,qc,artifact")
    bench.add_argument("--artifact-sample", type=int, default=50_000)

    sweep = sub.add_parser("sweep", help="a matrix of bench child processes")
    sweep.add_argument("--vcf", type=Path, nargs="+", required=True)
    sweep.add_argument("--backends", default="cyvcf2,text")
    sweep.add_argument("--stages", default="read,normalise,qc,artifact")
    sweep.add_argument("--artifact-sample", type=int, default=50_000)
    sweep.add_argument("--timeout", type=float, default=7200.0)
    sweep.add_argument("--out", type=Path, default=None)

    prof = sub.add_parser("profile", help="cProfile one read_vcf")
    prof.add_argument("--vcf", type=Path, required=True)
    prof.add_argument("--backend", default="cyvcf2")
    prof.add_argument("--top", type=int, default=30)
    prof.add_argument("--sort", default="tottime")

    pairs = sub.add_parser("pairs", help="post-filter candidate and pair arithmetic")
    pairs.add_argument("--callset", type=int, default=4_500_000)
    pairs.add_argument("--genes", type=int, default=19_500)
    pairs.add_argument("--max-pairs-per-gene", type=int, default=20)
    pairs.add_argument("--measured-pass-rate", type=float, default=None)
    pairs.add_argument("--genic-fraction", type=float, default=0.449)
    pairs.add_argument("--max-pairing-variants", type=int, default=24)

    modelbench = sub.add_parser("modelbench", help="VariantRecord build + hold cost")
    modelbench.add_argument("--count", type=int, default=500_000)

    pairbench = sub.add_parser("pairbench", help="measure generate_pairs scaling")
    pairbench.add_argument("--sizes", default="50,100,200,400,800,1600")
    pairbench.add_argument("--max-pairs-per-gene", type=int, default=20)
    pairbench.add_argument("--genes", type=int, default=1)

    intervals = sub.add_parser("intervals", help="gene-interval lookup throughput")
    intervals.add_argument("--lookups", type=int, default=200_000)
    intervals.add_argument(
        "--mane-dir",
        type=Path,
        default=Path("/Users/someshwar-tripathi/Contri/bio-hackathon/mva-resources/mane"),
    )
    intervals.add_argument("--vcf", type=Path, default=None)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    if args.command == "bench":
        result = run_bench(
            args.vcf,
            backend=args.backend,
            stages=tuple(args.stages.split(",")),
            artifact_sample=args.artifact_sample,
        )
        print(result.as_json())  # noqa: T201 -- machine-read by `sweep`
        return 0
    if args.command == "sweep":
        return run_sweep(
            args.vcf,
            tuple(args.backends.split(",")),
            tuple(args.stages.split(",")),
            artifact_sample=args.artifact_sample,
            timeout=args.timeout,
            out=args.out,
        )
    if args.command == "modelbench":
        return run_modelbench(args.count)
    if args.command == "pairbench":
        return run_pairbench(
            tuple(int(v) for v in args.sizes.split(",")),
            args.max_pairs_per_gene,
            args.genes,
        )
    if args.command == "profile":
        return run_profile(args.vcf, args.backend, top=args.top, sort=args.sort)
    if args.command == "pairs":
        return run_pairs(
            callset=args.callset,
            genes=args.genes,
            max_pairs_per_gene=args.max_pairs_per_gene,
            max_pairing_variants=args.max_pairing_variants,
            measured_pass_rate=args.measured_pass_rate,
            genic_fraction=args.genic_fraction,
        )
    return run_intervals(args.lookups, args.mane_dir, args.vcf)


if __name__ == "__main__":
    raise SystemExit(main())
