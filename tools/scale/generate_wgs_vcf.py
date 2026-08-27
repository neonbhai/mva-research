"""Deterministic generator for a WGS-scale synthetic VCF.

The output is a **throughput phantom**. It is shaped like the output of a
single-sample GRCh38 germline WGS pipeline -- 24 ``chr``-prefixed contigs,
Poisson-spaced sites at a realistic per-chromosome density, a GATK-style INFO
column, ``GT:AD:DP:GQ`` for one sample, ~10% indels, a small multi-allelic
fraction, non-PASS FILTER values, right-skewed depth and a 99-heavy GQ pile --
so that a parser doing realistic per-record work is measured doing realistic
per-record work.

Every REF base, every ALT base, every coordinate, every depth and every
genotype is fabricated from a seeded PRNG. **Nothing here is biologically
valid.** There is no reference genome behind it: a REF base is a random
nucleotide, not the base at that locus. The file measures bytes and objects per
second and nothing else, and no statement about variants, genes, inheritance or
disease may be derived from it.

The file carries ``##mva_synthetic=true`` because :mod:`mva.privacy.audit`
requires that marker on any genomic artifact, and a file that lacks it is
correctly refused.

Determinism: output is a pure function of ``(seed, variants)``. The same pair
produces a byte-identical body.
"""

from __future__ import annotations

import argparse
import math
import random
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Final

# ---------------------------------------------------------------------------
# Contig model
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Contig:
    """One GRCh38 primary contig and its relative variant density.

    ``density`` is a multiplier on the genome-average per-base variant rate.
    The values encode three real effects in the crudest defensible form: gene-
    and GC-rich chromosomes (17, 19, 22) call somewhat fewer variants per base;
    chrX in a male sample is haploid outside the PARs and yields roughly half
    the site count; chrY is haploid, largely heterochromatic and poorly
    mappable, so its callable yield is a small fraction of its length. They are
    approximations chosen to make the per-contig record distribution look like a
    real callset, not measurements.
    """

    name: str
    length: int
    density: float
    ploidy: int


#: GRCh38 primary assembly, chr1-22 + X + Y. Lengths are the real assembly
#: lengths (they go in the header and are what tabix indexes against); the
#: density multipliers are modelling choices, documented on :class:`Contig`.
CONTIGS: Final[tuple[Contig, ...]] = (
    Contig("chr1", 248_956_422, 1.00, 2),
    Contig("chr2", 242_193_529, 1.05, 2),
    Contig("chr3", 198_295_559, 1.05, 2),
    Contig("chr4", 190_214_555, 1.10, 2),
    Contig("chr5", 181_538_259, 1.05, 2),
    Contig("chr6", 170_805_979, 1.05, 2),
    Contig("chr7", 159_345_973, 1.05, 2),
    Contig("chr8", 145_138_636, 1.08, 2),
    Contig("chr9", 138_394_717, 0.95, 2),
    Contig("chr10", 133_797_422, 1.02, 2),
    Contig("chr11", 135_086_622, 1.00, 2),
    Contig("chr12", 133_275_309, 0.98, 2),
    Contig("chr13", 114_364_328, 1.05, 2),
    Contig("chr14", 107_043_718, 0.98, 2),
    Contig("chr15", 101_991_189, 0.95, 2),
    Contig("chr16", 90_338_345, 0.92, 2),
    Contig("chr17", 83_257_441, 0.85, 2),
    Contig("chr18", 80_373_285, 1.02, 2),
    Contig("chr19", 58_617_616, 0.80, 2),
    Contig("chr20", 64_444_167, 0.95, 2),
    Contig("chr21", 46_709_983, 0.98, 2),
    Contig("chr22", 50_818_468, 0.85, 2),
    Contig("chrX", 156_040_895, 0.48, 1),
    Contig("chrY", 57_227_415, 0.06, 1),
)

SAMPLE_ID: Final = "SYNTH_PROBAND01"
DEFAULT_SEED: Final = 20_260_828
DEFAULT_VARIANTS: Final = 4_500_000

_BASES: Final[tuple[str, ...]] = ("A", "C", "G", "T")

#: Site-class probabilities. ~10% indels is the usual single-sample WGS figure.
_P_INDEL: Final = 0.10
#: Fraction of sites carrying more than one ALT allele.
_P_MULTIALLELIC: Final = 0.013
#: Of the multi-allelic sites, the fraction carrying three ALTs rather than two.
_P_TRIALLELIC: Final = 0.08
#: Fraction of indels emitted in a non-parsimonious (shared-suffix) form, which
#: is what makes ``normalise.trim_and_left_align`` do real work rather than
#: no-op on every record.
_P_UNTRIMMED_INDEL: Final = 0.08

#: het:hom-alt of ~1.5:1 on the diploid contigs.
_P_HET: Final = 0.60
#: Phased calls, emitted with '|'. FORMAT stays exactly GT:AD:DP:GQ.
_P_PHASED: Final = 0.02
#: No-calls. Rare in a single-sample callset but non-zero, and they exercise the
#: reader's Zygosity.UNKNOWN path.
_P_NO_CALL: Final = 0.003

_P_PASS: Final = 0.92
_NON_PASS_FILTERS: Final[tuple[str, ...]] = (
    "LowQual",
    "VQSRTrancheSNP99.00to99.90",
    "VQSRTrancheSNP99.90to100.00",
    "VQSRTrancheINDEL99.00to99.90",
    "LowQual;VQSRTrancheSNP99.90to100.00",
)

#: Depth ~ Gamma(9, 3.6): mean ~32x, sd ~11x, right-skewed, as a 30x WGS run is.
_DP_SHAPE: Final = 9.0
_DP_SCALE: Final = 3.6

_MAX_INDEL_LEN: Final = 24


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------


def build_header(seed: int, variants: int) -> str:
    """The VCF header, including the synthetic marker the privacy audit requires."""
    lines = [
        "##fileformat=VCFv4.2",
        "##mva_synthetic=true",
        "##mva_synthetic_note=FULLY SYNTHETIC THROUGHPUT PHANTOM. Fabricated coordinates, "
        "fabricated alleles, fabricated depths and genotypes. Generated by "
        "tools/scale/generate_wgs_vcf.py to measure parser and model throughput. "
        "NOT biologically valid: no reference genome was consulted and no REF base "
        "corresponds to the assembly. No scientific conclusion may be drawn from it.",
        f"##mva_synthetic_seed={seed}",
        f"##mva_synthetic_target_records={variants}",
        "##reference=GRCh38",
        "##source=mva-scale-phantom",
    ]
    lines += [f"##contig=<ID={c.name},length={c.length}>" for c in CONTIGS]
    lines += [
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FILTER=<ID=LowQual,Description="Low quality variant call">',
        '##FILTER=<ID=VQSRTrancheSNP99.00to99.90,Description="Synthetic VQSR tranche">',
        '##FILTER=<ID=VQSRTrancheSNP99.90to100.00,Description="Synthetic VQSR tranche">',
        '##FILTER=<ID=VQSRTrancheINDEL99.00to99.90,Description="Synthetic VQSR tranche">',
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">',
        '##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Base quality rank sum">',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">',
        '##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Excess heterozygosity">',
        '##INFO=<ID=FS,Number=1,Type=Float,Description="Strand bias Fisher score">',
        '##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Max-likelihood allele count">',
        '##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Max-likelihood allele frequency">',
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS mapping quality">',
        '##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Mapping quality rank sum">',
        '##INFO=<ID=QD,Number=1,Type=Float,Description="Quality by depth">',
        '##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Read position rank sum">',
        '##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric odds ratio">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE_ID}",
    ]
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Allele construction
# ---------------------------------------------------------------------------


def _bases(rng: random.Random, length: int) -> str:
    return "".join(rng.choices(_BASES, k=length))


def _snv_alt(rng: random.Random, ref: str) -> str:
    """A single-base ALT that is not the REF base."""
    alt = _bases(rng, 1)
    while alt == ref:
        alt = _bases(rng, 1)
    return alt


def _indel_alleles(rng: random.Random, anchor: str) -> tuple[str, str]:
    """An anchored indel, occasionally emitted non-parsimoniously.

    The shared-suffix variant is deliberate: it is what a real caller emits at a
    decomposed multi-allelic site, and it is the only reason
    ``normalise.trim_and_left_align`` has anything to trim.
    """
    size = 1 + int(rng.expovariate(0.45))
    size = min(size, _MAX_INDEL_LEN)
    payload = _bases(rng, size)
    if rng.random() < 0.5:
        ref, alt = anchor, anchor + payload  # insertion
    else:
        ref, alt = anchor + payload, anchor  # deletion
    if rng.random() < _P_UNTRIMMED_INDEL:
        suffix = _bases(rng, 1)
        ref, alt = ref + suffix, alt + suffix
    return ref, alt


def _site_alleles(rng: random.Random) -> tuple[str, tuple[str, ...]]:
    """REF and the site's ALT list. Multi-allelic sites get 2 or 3 ALTs."""
    ref = _bases(rng, 1)
    n_alt = 1
    if rng.random() < _P_MULTIALLELIC:
        n_alt = 3 if rng.random() < _P_TRIALLELIC else 2

    if rng.random() >= _P_INDEL:
        alts: list[str] = []
        while len(alts) < n_alt:
            candidate = _snv_alt(rng, ref)
            if candidate not in alts:
                alts.append(candidate)
        return ref, tuple(alts)

    ref, first = _indel_alleles(rng, ref)
    alts = [first]
    while len(alts) < n_alt:
        # A second ALT at an indel site shares the site REF, so it is built by
        # extending or truncating that REF rather than re-anchoring.
        extra = ref + _bases(rng, 1 + int(rng.expovariate(0.6)))
        if extra not in alts:
            alts.append(extra)
    return ref, tuple(alts)


# ---------------------------------------------------------------------------
# Genotype construction
# ---------------------------------------------------------------------------


def _depth(rng: random.Random) -> int:
    return max(1, min(400, int(rng.gammavariate(_DP_SHAPE, _DP_SCALE))))


def _binomial(rng: random.Random, n: int, p: float) -> int:
    """Normal approximation to Binomial(n, p), clamped. Fast enough for 4.5M draws."""
    mean = n * p
    sd = math.sqrt(max(n * p * (1.0 - p), 1e-9))
    return max(0, min(n, round(rng.gauss(mean, sd))))


def _genotype_quality(rng: random.Random, depth: int, *, passing: bool) -> int:
    """GQ piles hard at 99 for confident calls and spreads out for the rest."""
    if passing and depth >= 20 and rng.random() < 0.78:
        return 99
    if depth < 8:
        return rng.randint(3, 45)
    return rng.randint(20, 98)


@dataclass(frozen=True, slots=True)
class _Call:
    gt: str
    ad: tuple[int, ...]
    depth: int
    allele_count: int


def _call(rng: random.Random, n_alt: int, ploidy: int, depth: int) -> _Call:
    """The sample column's genotype and per-allele depths."""
    if rng.random() < _P_NO_CALL:
        gt = "." if ploidy == 1 else "./."
        return _Call(gt=gt, ad=tuple([0] * (n_alt + 1)), depth=depth, allele_count=0)

    if ploidy == 1:
        chosen = rng.randint(1, n_alt)
        alt_reads = _binomial(rng, depth, 0.97)
        ad = [depth - alt_reads] + [0] * n_alt
        ad[chosen] = alt_reads
        return _Call(gt=str(chosen), ad=tuple(ad), depth=depth, allele_count=1)

    if n_alt > 1 and rng.random() < 0.45:
        # A genuine 1/2-style call: two different ALTs and few REF reads. This is
        # the case that makes the reader's multi-allelic AD slicing non-trivial.
        second = 2 if n_alt == 2 else rng.randint(2, n_alt)
        ref_reads = _binomial(rng, depth, 0.03)
        share = _binomial(rng, depth - ref_reads, 0.5)
        ad = [0] * (n_alt + 1)
        ad[0] = ref_reads
        ad[1] = share
        ad[second] = depth - ref_reads - share
        return _Call(gt=f"1/{second}", ad=tuple(ad), depth=depth, allele_count=2)

    het = rng.random() < _P_HET
    chosen = rng.randint(1, n_alt)
    alt_reads = _binomial(rng, depth, 0.5 if het else 0.98)
    ad = [depth - alt_reads] + [0] * n_alt
    ad[chosen] = alt_reads
    if not het:
        return _Call(gt=f"{chosen}/{chosen}", ad=tuple(ad), depth=depth, allele_count=2)
    if rng.random() < _P_PHASED:
        gt = f"0|{chosen}" if rng.random() < 0.5 else f"{chosen}|0"
    else:
        gt = f"0/{chosen}"
    return _Call(gt=gt, ad=tuple(ad), depth=depth, allele_count=1)


# ---------------------------------------------------------------------------
# Record rendering
# ---------------------------------------------------------------------------


def _info(rng: random.Random, call: _Call, n_alt: int, depth: int, qual: float) -> str:
    """A GATK-shaped INFO column.

    Present because it is roughly two thirds of a real VCF line, and both
    backends pay for it: the text backend splits it and htslib parses it. A
    phantom without it would flatter the parser by a wide margin.
    """
    per_allele = [call.allele_count if call.ad[i + 1] > 0 else 0 for i in range(n_alt)]
    ac = ",".join(str(value) for value in per_allele)
    af = ",".join(f"{value / 2.0:.3f}" for value in per_allele)
    return (
        f"AC={ac};AF={af};AN=2;"
        f"BaseQRankSum={rng.gauss(0, 0.9):.3f};"
        f"DP={depth};ExcessHet=0.0000;"
        f"FS={abs(rng.gauss(0, 2.0)):.3f};"
        f"MLEAC={ac};MLEAF={af};"
        f"MQ={min(60.0, rng.gauss(58.5, 3.0)):.2f};"
        f"MQRankSum={rng.gauss(0, 0.6):.3f};"
        f"QD={qual / max(depth, 1):.2f};"
        f"ReadPosRankSum={rng.gauss(0, 0.8):.3f};"
        f"SOR={abs(rng.gauss(0.7, 0.4)):.3f}"
    )


def _record(rng: random.Random, contig: Contig, position: int) -> str:
    ref, alts = _site_alleles(rng)
    n_alt = len(alts)
    depth = _depth(rng)
    passing = rng.random() < _P_PASS
    filter_field = "PASS" if passing else rng.choice(_NON_PASS_FILTERS)
    call = _call(rng, n_alt, contig.ploidy, depth)
    span = (9.0, 26.0) if passing else (0.4, 5.0)
    qual = round(depth * rng.uniform(*span), 1)
    gq = _genotype_quality(rng, depth, passing=passing)
    sample = f"{call.gt}:{','.join(str(v) for v in call.ad)}:{depth}:{gq}"
    return (
        f"{contig.name}\t{position}\t.\t{ref}\t{','.join(alts)}\t{qual}\t"
        f"{filter_field}\t{_info(rng, call, n_alt, depth, qual)}\tGT:AD:DP:GQ\t{sample}\n"
    )


# ---------------------------------------------------------------------------
# Generation
# ---------------------------------------------------------------------------


def plan_counts(variants: int) -> tuple[tuple[Contig, int], ...]:
    """Split a target record count across the contigs by ``length * density``."""
    weights = [c.length * c.density for c in CONTIGS]
    total = sum(weights)
    counts = [int(variants * w / total) for w in weights]
    # Push the rounding remainder onto chr1 so the total is exact and reproducible.
    counts[0] += variants - sum(counts)
    return tuple(zip(CONTIGS, counts, strict=True))


def _positions(rng: random.Random, contig: Contig, count: int) -> list[int]:
    """Poisson-spaced, strictly increasing positions along one contig.

    Exponential gaps rather than sort-a-random-sample: it is the right process
    for a variant landscape, it is O(n) with no sort, and it comes out already
    coordinate-sorted, which tabix requires.
    """
    if count <= 0:
        return []
    mean_gap = contig.length / (count + 1)
    positions: list[int] = []
    pos = 1
    for _ in range(count):
        pos += 1 + int(-mean_gap * math.log(rng.random()))
        if pos >= contig.length:
            break
        positions.append(pos)
    return positions


def write_vcf(path: Path, *, seed: int, variants: int, chunk_lines: int = 20_000) -> dict[str, int]:
    """Write a bgzipped, tabix-indexed synthetic WGS VCF. Returns per-contig counts.

    bgzf is written through pysam because this machine has no ``bgzip`` or
    ``tabix`` binary; the index is built by ``pysam.tabix_index``.
    """
    import pysam  # noqa: PLC0415 - optional native dependency, imported on demand

    rng = random.Random(seed)  # noqa: S311 - reproducible phantom data, not a secret
    written: dict[str, int] = {}

    handle = pysam.BGZFile(str(path), "wb", index=None)
    try:
        handle.write(build_header(seed, variants).encode("ascii"))
        for contig, count in plan_counts(variants):
            buffer: list[str] = []
            emitted = 0
            for position in _positions(rng, contig, count):
                buffer.append(_record(rng, contig, position))
                emitted += 1
                if len(buffer) >= chunk_lines:
                    handle.write("".join(buffer).encode("ascii"))
                    buffer.clear()
            if buffer:
                handle.write("".join(buffer).encode("ascii"))
            written[contig.name] = emitted
    finally:
        handle.close()

    pysam.tabix_index(str(path), preset="vcf", force=True)
    return written


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="generate_wgs_vcf",
        description="Generate a deterministic WGS-scale synthetic VCF (bgzip + tabix).",
    )
    parser.add_argument("--out", type=Path, required=True, help="Destination .vcf.gz path.")
    parser.add_argument("--variants", type=int, default=DEFAULT_VARIANTS)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    out: Path = args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    counts = write_vcf(out, seed=args.seed, variants=args.variants)
    size = out.stat().st_size
    print(  # noqa: T201 -- this tool's CLI report
        f"wrote {sum(counts.values())} records to {out} "
        f"({size / 1e6:.1f} MB bgzf, seed={args.seed})",
        file=sys.stdout,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
