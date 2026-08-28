"""Scale harness for the stages AFTER ingestion: annotate, ledger, select.

``tools/scale/harness.py`` measures ``mva.ingestion`` against a real bgzf file.
This one measures the three stages that come next, against synthetic
``VariantRecord`` objects built in memory — no VCF, no reference, no patient data
of any kind. That is deliberate: the question here is "how many bytes does one
annotated record plus its evidence cost, and does the stage hold all of them at
once", and a parser in front of it only adds noise and 40 seconds per run.

Every number this prints is observed. Nothing is modelled. Where a configuration
could not complete, it says so and reports where it died.

Method
------
Peak memory is ``resource.getrusage(RUSAGE_SELF).ru_maxrss``, measured **in a
child process per configuration**, because a peak is a high-water mark: once one
mode has allocated 12 GB every later measurement in the same process inherits
that mark and reports it as its own. The RSS helpers are imported from
``tools.scale.harness`` rather than re-implemented, so the two harnesses'
numbers are comparable.

A watchdog thread samples RSS and kills the child at ``--rss-cap`` GB. The
machine this runs on has 24 GB and is shared with other work; a run that pages
30 GB of Pydantic models measures the swap subsystem and takes the machine down
with it. A killed configuration is reported as killed, with the last RSS
observed — which is a measurement, not a failure to measure.

The phantom
-----------
``synthetic_records`` is a **throughput phantom**, exactly as
``tools/scale/generate_wgs_vcf.py`` is. Every coordinate, allele, genotype,
gene symbol, consequence term and allele frequency is fabricated from a seeded
hash. It is not biologically valid, no result from it says anything about
variants, genes or disease, and it must never be described otherwise. Its only
job is to have the right *shape* and the right *rates*:

* 44.9% of records fall inside a gene, matching the measured MANE hit rate
  (``docs/scale-report.md`` §6);
* 62% of those carry frequency data, so 38% exercise the GP-14 unknown path;
* 1.2% of genic records are coding or splice-relevant;
* 8% of genic records carry an impact of ``None`` — NOT ASSESSED, the shape a
  MANE interval join produces (ADR 0016);
* a handful carry a curated pathogenic assertion.

Subcommands
-----------
``annotate``  batch (``annotate_variants``) vs streaming (``iter_annotated``).
``ledger``    in-memory vs spilled ``EvidenceLedger``.
``select``    the selection stage's throughput and its drop counts by reason.
``pipeline``  annotate -> select -> artifact, the composition the orchestrator
              is being asked to adopt.
``digest``    sha256 over the whole streamed output, for the GP-30 hash-seed proof.
``sweep``     runs a matrix of the above, one child each, and prints a table.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import subprocess
import sys
import threading
import time
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Final

from mva.annotation.base import AdapterSet
from mva.models.variant import (
    ClinicalAssertion,
    ConsequenceAnnotation,
    FilterStatus,
    Genotype,
    ImpactSeverity,
    PopulationFrequency,
    VariantRecord,
    Zygosity,
)
from tools.scale.harness import REPO_ROOT, current_rss_bytes, peak_rss_bytes

GIB: Final = 1 << 30

#: Real GRCh38 lengths, so the phantom's coordinates land in a plausible range.
_CONTIGS: Final[tuple[tuple[str, int], ...]] = (
    ("chr1", 248_956_422),
    ("chr2", 242_193_529),
    ("chr3", 198_295_559),
    ("chr4", 190_214_555),
    ("chr5", 181_538_259),
    ("chr6", 170_805_979),
    ("chr7", 159_345_973),
    ("chr8", 145_138_636),
    ("chr9", 138_394_717),
    ("chr10", 133_797_422),
    ("chr11", 135_086_622),
    ("chr12", 133_275_309),
    ("chr13", 114_364_328),
    ("chr14", 107_043_718),
    ("chr15", 101_991_189),
    ("chr16", 90_338_345),
    ("chr17", 83_257_441),
    ("chr18", 80_373_285),
    ("chr19", 58_617_616),
    ("chr20", 64_444_167),
    ("chr21", 46_709_983),
    ("chr22", 50_818_468),
    ("chrX", 156_040_895),
    ("chrY", 57_227_415),
)

_BASES: Final = ("A", "C", "G", "T")

_CODING_TERMS: Final = ("missense_variant", "stop_gained", "frameshift_variant")
_SPLICE_TERMS: Final = ("splice_region_variant&synonymous_variant", "splice_acceptor_variant")
_QUIET_TERMS: Final = ("intron_variant", "synonymous_variant", "downstream_gene_variant")


def _draw(seed: int, index: int, salt: int) -> int:
    """A deterministic pseudo-random 32-bit draw. No RNG object, no global state."""
    digest = hashlib.blake2b(
        f"{seed}:{index}:{salt}".encode(), digest_size=4, usedforsecurity=False
    )
    return int.from_bytes(digest.digest(), "big")


def synthetic_records(count: int, *, seed: int = 20260828) -> Iterator[VariantRecord]:
    """Yield ``count`` fabricated records in coordinate order. NOT biological data."""
    per_contig = max(count // len(_CONTIGS), 1)
    emitted = 0
    for contig_index, (contig, length) in enumerate(_CONTIGS):
        target = per_contig if contig_index < len(_CONTIGS) - 1 else count - emitted
        step = max(length // max(target, 1), 1)
        position = 1
        for i in range(target):
            if emitted >= count:
                return
            draw = _draw(seed, emitted, 0)
            ref = _BASES[draw & 3]
            alt = _BASES[(draw >> 2) & 3]
            if alt == ref:
                alt = _BASES[(draw >> 4) & 3]
                if alt == ref:
                    alt = "A" if ref != "A" else "C"
            depth = 10 + ((draw >> 6) % 60)
            alt_reads = max(1, depth // 2 + ((draw >> 12) % 5) - 2)
            hom = (draw >> 17) % 5 == 0
            yield VariantRecord(
                coordinate={  # type: ignore[arg-type] - pydantic coerces the mapping
                    "build": "GRCh38",
                    "contig": contig,
                    "position": position + (i * step) % max(length - 1, 1),
                    "ref": ref,
                    "alt": alt,
                },
                genotype=Genotype(
                    zygosity=Zygosity.HOM_ALT if hom else Zygosity.HET,
                    genotype_string="1/1" if hom else "0/1",
                    depth=depth,
                    ref_reads=depth - alt_reads,
                    alt_reads=alt_reads,
                    genotype_quality=99,
                ),
                filter_status=FilterStatus.PASS,
                source_artifact="scale-stage-phantom",
                source_line_index=emitted,
            )
            emitted += 1
    return


# ---------------------------------------------------------------------------
# Phantom adapters
# ---------------------------------------------------------------------------


class PhantomConsequenceAdapter:
    """A fabricated consequence source with realistic hit and impact rates."""

    name = "phantom-consequence"
    version = "0.0.0-synthetic"
    synthetic = True

    def annotate(
        self, variant_ids: Sequence[str]
    ) -> Mapping[str, tuple[ConsequenceAnnotation, ...]]:
        out: dict[str, tuple[ConsequenceAnnotation, ...]] = {}
        for variant_id in variant_ids:
            draw = int.from_bytes(
                hashlib.blake2b(variant_id.encode(), digest_size=4, usedforsecurity=False).digest(),
                "big",
            )
            if draw % 1000 >= 449:  # 44.9% genic, matching the measured MANE rate
                continue
            gene = f"SYNTHG{draw % 19500:05d}"
            bucket = (draw >> 10) % 1000
            if bucket < 12:
                terms, impact = (_CODING_TERMS[bucket % 3],), ImpactSeverity.HIGH
            elif bucket < 20:
                terms, impact = (_SPLICE_TERMS[bucket % 2],), ImpactSeverity.LOW
            elif bucket < 100:
                # NOT ASSESSED: gene located, no consequence computed (ADR 0016).
                terms, impact = ("transcript_variant",), None
            else:
                terms, impact = (_QUIET_TERMS[bucket % 3],), ImpactSeverity.MODIFIER
            out[variant_id] = tuple(
                ConsequenceAnnotation(
                    gene_symbol=gene,
                    transcript_id=f"ENST{(draw >> 4) % 10**11:011d}.{n + 1}",
                    consequence_terms=terms,
                    impact=impact,
                    is_mane_select=n == 0,
                    source_tool=self.name,
                    source_tool_version=self.version,
                )
                # Two transcripts for one variant in six: every transcript is
                # retained (GP-14 / adapter rule 3), so this exercises the fan-out.
                for n in range(2 if (draw >> 20) % 6 == 0 else 1)
            )
        return out


class PhantomFrequencyAdapter:
    """Fabricated frequencies. 38% of variants get nothing — the GP-14 path."""

    name = "phantom-frequency"
    version = "0.0.0-synthetic"
    synthetic = True

    def frequencies(
        self, variant_ids: Sequence[str]
    ) -> Mapping[str, tuple[PopulationFrequency, ...]]:
        out: dict[str, tuple[PopulationFrequency, ...]] = {}
        for variant_id in variant_ids:
            draw = int.from_bytes(
                hashlib.blake2b(variant_id.encode(), digest_size=4, usedforsecurity=False).digest(),
                "big",
            )
            if (draw >> 3) % 100 >= 62:
                continue
            frequency = ((draw >> 7) % 100_000) / 100_000.0
            out[variant_id] = (
                PopulationFrequency(
                    source="phantom_gnomad",
                    version="v0.0-synthetic",
                    population="global",
                    allele_frequency=frequency,
                    allele_count=int(frequency * 152_312),
                    allele_number=152_312,
                ),
            )
        return out


class PhantomClinicalAdapter:
    """A curated-assertion source that knows about roughly one variant in 5,000."""

    name = "phantom-clinical"
    version = "0.0.0-synthetic"
    synthetic = True

    def assertions(self, variant_ids: Sequence[str]) -> Mapping[str, tuple[ClinicalAssertion, ...]]:
        out: dict[str, tuple[ClinicalAssertion, ...]] = {}
        for variant_id in variant_ids:
            draw = int.from_bytes(
                hashlib.blake2b(variant_id.encode(), digest_size=4, usedforsecurity=False).digest(),
                "big",
            )
            if (draw >> 11) % 5000:
                continue
            out[variant_id] = (
                ClinicalAssertion(
                    source="phantom_clinvar",
                    version="v0.0-synthetic",
                    accession=f"SCV{draw % 10**9:09d}",
                    significance="Pathogenic",
                    review_status="criteria_provided,_single_submitter",
                    star_rating=1,
                ),
            )
        return out


def phantom_adapters(*, clinical: bool = True) -> AdapterSet:
    return AdapterSet(
        consequence=PhantomConsequenceAdapter(),
        frequency=PhantomFrequencyAdapter(),
        clinical=PhantomClinicalAdapter() if clinical else None,
    )


# ---------------------------------------------------------------------------
# Watchdog + result record
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class StageResult:
    subcommand: str
    mode: str
    records: int
    seconds: float = 0.0
    records_per_second: float = 0.0
    peak_rss: int = 0
    rss_after: int = 0
    baseline_rss: int = 0
    killed: bool = False
    notes: dict[str, Any] = field(default_factory=dict)

    def as_json(self) -> str:
        return json.dumps(asdict(self), indent=2, sort_keys=True)


def _install_watchdog(cap_bytes: int, result: StageResult) -> None:
    """Kill this process before it drives the machine into swap."""
    if cap_bytes <= 0:
        return

    def watch() -> None:
        while True:
            rss = current_rss_bytes()
            if rss > cap_bytes:
                result.killed = True
                result.peak_rss = max(result.peak_rss, peak_rss_bytes())
                result.rss_after = rss
                result.notes["killed_at_gib"] = rss / GIB
                sys.stdout.write(result.as_json() + "\n")
                sys.stdout.flush()
                os._exit(9)
            time.sleep(0.2)

    threading.Thread(target=watch, daemon=True).start()


# ---------------------------------------------------------------------------
# annotate
# ---------------------------------------------------------------------------


def run_annotate(records: int, *, mode: str, batch_size: int, cap: int) -> StageResult:
    """Batch ``annotate_variants`` against streaming ``iter_annotated``."""
    from mva.annotation.service import annotate_variants, iter_annotated  # noqa: PLC0415
    from mva.clock import demo_clock  # noqa: PLC0415

    result = StageResult(subcommand="annotate", mode=mode, records=records)
    _install_watchdog(cap, result)
    adapters = phantom_adapters()
    clock = demo_clock()
    gc.collect()
    result.baseline_rss = current_rss_bytes()
    start = time.perf_counter()

    if mode == "batch":
        # Exactly what orchestrator.py does today: a materialised callset in, a
        # second materialised callset plus every evidence item out, both alive.
        source = list(synthetic_records(records))
        annotation = annotate_variants(source, adapters=adapters, clock=clock)
        result.notes["variants_out"] = len(annotation.variants)
        result.notes["evidence_out"] = len(annotation.evidence)
        result.notes["coverage"] = annotation.coverage
        result.notes["warnings"] = list(annotation.warnings)
    else:
        stream = iter_annotated(
            synthetic_records(records),
            adapters=adapters,
            clock=clock,
            batch_size=batch_size,
        )
        variants = 0
        evidence = 0
        for annotated in stream:
            variants += 1
            evidence += len(annotated.evidence)
        result.notes["variants_out"] = variants
        result.notes["evidence_out"] = evidence
        result.notes["coverage"] = stream.coverage()
        result.notes["warnings"] = list(stream.warnings())
        result.notes["batch_size"] = batch_size

    result.seconds = time.perf_counter() - start
    result.records_per_second = records / result.seconds if result.seconds else 0.0
    result.rss_after = current_rss_bytes()
    result.peak_rss = peak_rss_bytes()
    return result


# ---------------------------------------------------------------------------
# ledger
# ---------------------------------------------------------------------------


def run_ledger(records: int, *, mode: str, spill_dir: Path, cap: int) -> StageResult:
    """An EvidenceLedger fed the annotation evidence for ``records`` variants."""
    from mva.annotation.service import iter_annotated  # noqa: PLC0415
    from mva.clock import demo_clock  # noqa: PLC0415
    from mva.evidence.ledger import EvidenceLedger  # noqa: PLC0415

    result = StageResult(subcommand="ledger", mode=mode, records=records)
    _install_watchdog(cap, result)
    # No spill_threshold override: measure the shipped default.
    ledger = EvidenceLedger(
        run_id="SCALE-PHANTOM", spill_dir=spill_dir if mode == "spill" else None
    )
    gc.collect()
    result.baseline_rss = current_rss_bytes()
    start = time.perf_counter()

    stream = iter_annotated(
        synthetic_records(records), adapters=phantom_adapters(), clock=demo_clock()
    )
    for annotated in stream:
        ledger.extend(annotated.evidence)
    result.notes["ledger_items"] = len(ledger)
    result.notes["spilled"] = ledger.spilled
    write_seconds = time.perf_counter() - start

    scan_start = time.perf_counter()
    digest = hashlib.sha256()
    scanned = 0
    for item in ledger.iter_items():
        digest.update(item.evidence_id.encode())
        scanned += 1
    result.notes["ordered_scan_items"] = scanned
    result.notes["ordered_scan_sha256"] = digest.hexdigest()
    result.notes["write_seconds"] = write_seconds
    result.notes["ordered_scan_seconds"] = time.perf_counter() - scan_start
    if ledger.spilled:
        # Report the spill file's size before close() removes it.
        spill_file = spill_dir / "evidence-ledger-SCALE-PHANTOM.sqlite"
        result.notes["spill_bytes"] = spill_file.stat().st_size if spill_file.is_file() else 0

    result.seconds = time.perf_counter() - start
    result.records_per_second = records / result.seconds if result.seconds else 0.0
    result.rss_after = current_rss_bytes()
    result.peak_rss = peak_rss_bytes()
    ledger.close()
    return result


# ---------------------------------------------------------------------------
# select
# ---------------------------------------------------------------------------


def run_select(records: int, *, mode: str, cap: int) -> StageResult:
    """Annotation straight into selection, counting what reaches pairing."""
    from mva.annotation.service import iter_annotated  # noqa: PLC0415
    from mva.clock import demo_clock  # noqa: PLC0415
    from mva.config import FrequencyThresholds  # noqa: PLC0415
    from mva.prioritization.selection import SelectionThresholds, iter_selected  # noqa: PLC0415

    result = StageResult(subcommand="select", mode=mode, records=records)
    _install_watchdog(cap, result)
    clock = demo_clock()
    gc.collect()
    result.baseline_rss = current_rss_bytes()
    start = time.perf_counter()

    annotated = iter_annotated(synthetic_records(records), adapters=phantom_adapters(), clock=clock)
    selection = iter_selected(
        (item.variant for item in annotated),
        frequency=FrequencyThresholds(),
        thresholds=SelectionThresholds(enabled=mode != "disabled"),
        clock=clock,
    )
    kept = list(selection)
    report = selection.report()
    result.notes["selected"] = len(kept)
    result.notes["report"] = report.as_payload()
    result.notes["evidence_items"] = len(selection.evidence())
    result.seconds = time.perf_counter() - start
    result.records_per_second = records / result.seconds if result.seconds else 0.0
    result.rss_after = current_rss_bytes()
    result.peak_rss = peak_rss_bytes()
    return result


# ---------------------------------------------------------------------------
# pipeline: annotate -> artifact + ledger -> select
# ---------------------------------------------------------------------------


def run_pipeline(records: int, *, mode: str, out_dir: Path, cap: int) -> StageResult:
    """The composition ``docs/handoff-scale.md`` asks the orchestrator to adopt.

    One pass, end to end:

        records -> annotate -> ledger -> annotated.json -> hard filter
                -> soft flags -> selection -> the variants pairing will see

    The artifact is written from the point BEFORE hard filtering, through the push
    sink (``mva.pipeline.JsonRowsSink``), because it must still hold every
    annotated record while the same pass feeds the stages that drop most of them.
    """
    from mva.annotation.service import iter_annotated  # noqa: PLC0415
    from mva.clock import demo_clock  # noqa: PLC0415
    from mva.config import FrequencyThresholds, QualityThresholds  # noqa: PLC0415
    from mva.evidence.ledger import EvidenceLedger  # noqa: PLC0415
    from mva.models.genome import GenomeBuild  # noqa: PLC0415
    from mva.pipeline import JsonRowsSink  # noqa: PLC0415
    from mva.prioritization.filters import iter_hard_filtered, iter_soft_flagged  # noqa: PLC0415
    from mva.prioritization.selection import SelectionThresholds, iter_selected  # noqa: PLC0415

    result = StageResult(subcommand="pipeline", mode=mode, records=records)
    _install_watchdog(cap, result)
    out_dir.mkdir(parents=True, exist_ok=True)
    artifact = out_dir / "annotated.json"
    clock = demo_clock()
    ledger = EvidenceLedger(run_id="SCALE-PHANTOM", spill_dir=out_dir)
    gc.collect()
    result.baseline_rss = current_rss_bytes()
    start = time.perf_counter()

    annotated = iter_annotated(synthetic_records(records), adapters=phantom_adapters(), clock=clock)
    with artifact.open("w", encoding="utf-8") as handle:
        sink = JsonRowsSink(handle)
        recorded = _record_and_emit(annotated, ledger, sink)
        filtered = iter_hard_filtered(recorded, expected_build=GenomeBuild.GRCH38)
        flagged = iter_soft_flagged(
            filtered, frequency=FrequencyThresholds(), quality=QualityThresholds()
        )
        selection = iter_selected(
            flagged,
            frequency=FrequencyThresholds(),
            thresholds=SelectionThresholds(),
            clock=clock,
        )
        selected = list(selection)
        sink.close()

    result.notes["hard_filter_counts"] = filtered.counts()
    result.notes["artifact_rows"] = sink.row_count
    result.notes["artifact_bytes"] = artifact.stat().st_size
    result.notes["selected"] = len(selected)
    result.notes["ledger_items"] = len(ledger)
    result.notes["spilled"] = ledger.spilled
    spill_file = out_dir / "evidence-ledger-SCALE-PHANTOM.sqlite"
    result.notes["spill_bytes"] = spill_file.stat().st_size if spill_file.is_file() else 0
    result.notes["report"] = selection.report().as_payload()
    result.notes["annotation_warnings"] = list(annotated.warnings())
    result.seconds = time.perf_counter() - start
    result.records_per_second = records / result.seconds if result.seconds else 0.0
    result.rss_after = current_rss_bytes()
    result.peak_rss = peak_rss_bytes()
    ledger.close()
    artifact.unlink(missing_ok=True)
    return result


def _record_and_emit(annotated: Any, ledger: Any, sink: Any) -> Iterator[VariantRecord]:
    """Ledger the evidence, write the artifact row, pass the record on."""
    for item in annotated:
        ledger.extend(item.evidence)
        sink.write(item.variant.model_dump(mode="json"))
        yield item.variant


def _tee(annotated: Any, ledger: Any) -> Iterator[VariantRecord]:
    """Push evidence into the ledger and pass the record on. One pass, no copies."""
    for item in annotated:
        ledger.extend(item.evidence)
        yield item.variant


# ---------------------------------------------------------------------------
# digest: the GP-30 proof
# ---------------------------------------------------------------------------


def run_digest(records: int, *, cap: int) -> StageResult:
    """sha256 over every annotated record, every evidence item and every verdict.

    Run twice under different ``PYTHONHASHSEED`` values. Equal digests are the
    byte-identity claim; anything else is a determinism defect (GP-30).
    """
    from mva.annotation.service import iter_annotated  # noqa: PLC0415
    from mva.clock import demo_clock  # noqa: PLC0415
    from mva.config import FrequencyThresholds  # noqa: PLC0415
    from mva.determinism import canonical_json  # noqa: PLC0415
    from mva.evidence.ledger import EvidenceLedger  # noqa: PLC0415
    from mva.prioritization.selection import SelectionThresholds, iter_selected  # noqa: PLC0415

    result = StageResult(subcommand="digest", mode="stream", records=records)
    _install_watchdog(cap, result)
    clock = demo_clock()
    ledger = EvidenceLedger(run_id="SCALE-PHANTOM", spill_dir=Path(os.environ["MVA_SCALE_TMP"]))
    variants_digest = hashlib.sha256()
    verdicts_digest = hashlib.sha256()
    gc.collect()
    result.baseline_rss = current_rss_bytes()
    start = time.perf_counter()

    annotated = iter_annotated(synthetic_records(records), adapters=phantom_adapters(), clock=clock)
    selection = iter_selected(
        _tee(annotated, ledger),
        frequency=FrequencyThresholds(),
        thresholds=SelectionThresholds(),
        clock=clock,
    )
    for decision in selection.decisions():
        variants_digest.update(canonical_json(decision.variant.model_dump(mode="json")).encode())
        verdicts_digest.update(f"{decision.reason}|{','.join(decision.notes)}\n".encode())

    ledger_digest = hashlib.sha256()
    for item in ledger.iter_items():
        ledger_digest.update(canonical_json(item.model_dump(mode="json")).encode())

    result.notes["variants_sha256"] = variants_digest.hexdigest()
    result.notes["verdicts_sha256"] = verdicts_digest.hexdigest()
    result.notes["ledger_sha256"] = ledger_digest.hexdigest()
    result.notes["report_sha256"] = hashlib.sha256(
        canonical_json(selection.report().as_payload()).encode()
    ).hexdigest()
    result.notes["annotation_warnings_sha256"] = hashlib.sha256(
        canonical_json(list(annotated.warnings())).encode()
    ).hexdigest()
    result.notes["ledger_items"] = len(ledger)
    result.notes["spilled"] = ledger.spilled
    result.notes["hash_seed"] = os.environ.get("PYTHONHASHSEED", "<unset>")
    result.seconds = time.perf_counter() - start
    result.rss_after = current_rss_bytes()
    result.peak_rss = peak_rss_bytes()
    ledger.close()
    return result


# ---------------------------------------------------------------------------
# sweep
# ---------------------------------------------------------------------------


def _child(args: Sequence[str], *, timeout: float, env: Mapping[str, str] | None = None) -> Any:
    cmd = [sys.executable, "-m", "tools.scale.stage_harness", *args]
    child_env = dict(os.environ)
    if env:
        child_env.update(env)
    started = time.perf_counter()
    proc = subprocess.run(  # noqa: S603 - fixed argv, no shell
        cmd,
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
        timeout=timeout,
        env=child_env,
    )
    elapsed = time.perf_counter() - started
    # The child prints exactly one INDENTED JSON object, so parse the whole of
    # stdout rather than its last line -- which is "}" and taught this harness to
    # report every configuration as zero seconds and no peak.
    text = proc.stdout.strip()
    payload: dict[str, Any]
    try:
        payload = json.loads(text) if text else {}
    except json.JSONDecodeError:
        start = text.find("{")
        try:
            payload = json.loads(text[start:]) if start >= 0 else {}
        except json.JSONDecodeError:
            payload = {"parse_failed": True, "stdout_tail": text.splitlines()[-3:]}
    payload.setdefault("killed", proc.returncode != 0)
    payload["exit_code"] = proc.returncode
    payload["wall_seconds"] = elapsed
    if proc.returncode not in (0, 9):
        payload["stderr_tail"] = proc.stderr.strip().splitlines()[-6:]
    return payload


def _gib(value: object) -> str:
    return f"{float(value or 0) / GIB:8.3f}" if isinstance(value, int | float) else "       -"


def run_sweep(sizes: Sequence[int], *, cap: int, tmp: Path, timeout: float) -> list[Any]:
    rows: list[Any] = []
    for size in sizes:
        for subcommand, mode in (
            ("annotate", "batch"),
            ("annotate", "stream"),
            ("ledger", "memory"),
            ("ledger", "spill"),
            ("select", "enabled"),
            ("pipeline", "stream"),
        ):
            row = _child(
                [
                    subcommand,
                    "--records",
                    str(size),
                    "--mode",
                    mode,
                    "--rss-cap",
                    str(cap),
                    "--tmp",
                    str(tmp),
                ],
                timeout=timeout,
            )
            row.setdefault("subcommand", subcommand)
            row.setdefault("mode", mode)
            row.setdefault("records", size)
            rows.append(row)
            status = "KILLED" if row.get("killed") else "ok"
            print(  # noqa: T201 - this is a CLI
                f"{subcommand:9s} {mode:8s} {size:>9,}  "
                f"{row.get('seconds', 0.0):8.2f}s  "
                f"peak {_gib(row.get('peak_rss'))} GiB  {status}",
                flush=True,
            )
    return rows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "subcommand", choices=("annotate", "ledger", "select", "pipeline", "digest", "sweep")
    )
    parser.add_argument("--records", type=int, default=100_000)
    parser.add_argument("--mode", default="stream")
    parser.add_argument("--batch-size", type=int, default=5_000)
    parser.add_argument(
        "--rss-cap",
        type=float,
        default=14.0,
        help="Gibibytes. The child kills itself above this. 0 disables the watchdog.",
    )
    parser.add_argument("--tmp", type=Path, default=Path("/tmp/mva-scale"))  # noqa: S108
    parser.add_argument("--sizes", default="100000,500000,1000000,4500000")
    parser.add_argument("--timeout", type=float, default=3600.0)
    args = parser.parse_args(argv)

    cap = int(args.rss_cap * GIB)
    tmp = Path(args.tmp)
    tmp.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MVA_SCALE_TMP", str(tmp))

    if args.subcommand == "sweep":
        sizes = [int(value) for value in args.sizes.split(",") if value]
        rows = run_sweep(sizes, cap=cap, tmp=tmp, timeout=args.timeout)
        print(json.dumps(rows, indent=2, sort_keys=True))  # noqa: T201
        return 0

    result: StageResult
    if args.subcommand == "annotate":
        result = run_annotate(args.records, mode=args.mode, batch_size=args.batch_size, cap=cap)
    elif args.subcommand == "ledger":
        result = run_ledger(args.records, mode=args.mode, spill_dir=tmp, cap=cap)
    elif args.subcommand == "select":
        result = run_select(args.records, mode=args.mode, cap=cap)
    elif args.subcommand == "pipeline":
        result = run_pipeline(args.records, mode=args.mode, out_dir=tmp, cap=cap)
    else:
        result = run_digest(args.records, cap=cap)

    print(result.as_json())  # noqa: T201 - this is a CLI
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
