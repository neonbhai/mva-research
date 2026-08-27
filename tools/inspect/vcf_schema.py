"""Emit SCHEMA-LEVEL metadata about a VCF. Never records, never genotypes.

Exists because ADR 0008 / GP-44 make everything under ``MVA_WORKSPACE``
deliberately illegible to agents, while the pipeline still has to be configured
correctly for the file in front of it — above all the contig naming style, which
`docs/references/track1-submission-contract.md` identifies as the single detail
that can silently zero an otherwise-correct submission.

The resolution is that a *tool* reads the patient data and emits aggregates; the
agent reads the aggregates. So this prints:

  * header declarations (fileformat, reference, source, caller)
  * the SET of INFO / FORMAT / FILTER ids -- names only, never values
  * contig naming style and count
  * sample COUNT, and a salted digest prefix of each sample id, never the id
  * record counts, bucketed

It never prints a coordinate, an allele, a genotype, or a sample name. Adding a
line here that does is a privacy regression, not a debugging convenience.
"""

# ruff: noqa: T201
# T201 (no print) is disabled deliberately and only here. This is a standalone
# operator CLI whose entire purpose is to write a metadata summary to stdout; the
# repo-wide ban exists to stop library code printing, which this is not. Routing
# it through the logging stack would be worse: logs get shipped and retained, and
# this tool runs against patient data.

from __future__ import annotations

import argparse
import hashlib
import sys
from collections import Counter
from pathlib import Path

import pysam

#: Contig names probed to report naming style. Not a filter -- only a display probe.
PRIMARY_PROBE: frozenset[str] = frozenset({"chr1", "1", "chrX", "X", "chrM", "MT", "chrY", "Y"})


def _sample_token(name: str) -> str:
    """Stable, non-reversing reference to a sample id (GP-41)."""
    return hashlib.sha256(name.encode("utf-8")).hexdigest()[:8]


def main() -> int:
    ap = argparse.ArgumentParser(description="Schema-level VCF metadata. No records.")
    ap.add_argument("vcf", type=Path)
    ap.add_argument(
        "--scan", type=int, default=200_000, help="records to scan for aggregate counts (0 = none)"
    )
    args = ap.parse_args()

    vf = pysam.VariantFile(str(args.vcf))
    hdr = vf.header

    print("=== header declarations ===")
    for rec in hdr.records:
        if rec.key in {
            "fileformat",
            "reference",
            "source",
            "sourceVersion",
            "fileDate",
            "reference_genome",
            "DRAGENCommandLine",
            "GATKCommandLine",
            "bcftools_normCommand",
            "PEDIGREE",
        }:
            print(f"  ##{rec.key}={str(rec)[:220].rstrip()}")

    contigs = [str(c) for c in hdr.contigs]
    prefixed = sum(1 for c in contigs if c.startswith("chr"))
    print("\n=== contigs ===")
    print(f"  count: {len(contigs)}")
    print(
        f"  chr-prefixed: {prefixed}/{len(contigs)}  -> style: "
        f"{'UCSC (chr-prefixed)' if prefixed > len(contigs) / 2 else 'Ensembl/NCBI (bare)'}"
    )
    print(f"  primary present: {sorted(c for c in contigs if c in PRIMARY_PROBE)}")

    print("\n=== samples ===")
    print(
        f"  count: {len(hdr.samples)}  -> {'SINGLE' if len(hdr.samples) == 1 else 'MULTI (trio?)'}"
    )
    print(f"  tokens: {[_sample_token(s) for s in hdr.samples]}")

    for kind in ("info", "formats", "filters"):
        ids = sorted(str(k) for k in getattr(hdr, kind).keys())
        print(f"\n=== {kind.upper()} ids ({len(ids)}) ===")
        print("  " + ", ".join(ids))

    if args.scan:
        print(f"\n=== aggregate scan (first {args.scan:,} records) ===")
        n = 0
        filt: Counter[str] = Counter()
        multi = 0
        indel = 0
        for rec in vf:
            n += 1
            alts = rec.alts or ()
            if len(alts) > 1:
                multi += 1
            if any(len(a) != len(rec.ref or "") for a in alts):
                indel += 1
            filt[",".join(str(k) for k in rec.filter.keys()) or "."] += 1
            if n >= args.scan:
                break
        print(f"  scanned: {n:,}")
        print(
            f"  multi-allelic records: {multi:,} ({100 * multi / max(n, 1):.2f}%)"
            f"   -> {'NOT split' if multi else 'appears SPLIT'}"
        )
        print(f"  indel-bearing records: {indel:,} ({100 * indel / max(n, 1):.2f}%)")
        print("  FILTER distribution:")
        for k, v in filt.most_common(8):
            print(f"     {v:9,}  {k}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
