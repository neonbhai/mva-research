"""Regenerate ``clinvar_slice.vcf.gz`` — a small, real slice of the ClinVar release.

The 185 MB NCBI ClinVar GRCh38 VCF is NOT committed. This script cuts a fixed,
documented set of regions out of it so the adapter's tests run against genuine
ClinVar records — real CLNSIG spellings, real review statuses, real percent- and
underscore-escaping — instead of records someone hand-wrote to match the parser.

Reproduce (the input is the file the acquisition step downloads):

    uv run python tests/fixtures/clinvar/make_fixture.py \
        --source ~/Contri/bio-hackathon/mva-resources/clinvar/clinvar.vcf.gz

Equivalent with htslib on the PATH (this repo has no bgzip/tabix binary, which is
why the script uses pysam instead):

    { tabix -H "$SRC"; tabix "$SRC" 7:117509000-117510000 15:40199000-40206000 \
        17:43045000-43048000; } > clinvar_slice.vcf
    bgzip -f clinvar_slice.vcf && tabix -p vcf clinvar_slice.vcf.gz

Why these three regions, and nothing else:

* ``7:117509000-117510000`` — CFTR. Carries the only ``practice_guideline``
  record in reach, i.e. the 4-star case. There are 51 in the whole release.
* ``15:40199000-40206000`` — BUB1B, the gene this project exists for. Supplies
  pathogenic, benign, VUS, conflicting and indel records at a locus whose name a
  reviewer will recognise.
* ``17:43045000-43048000`` — BRCA1 3' end. The densest available source of
  ``reviewed_by_expert_panel`` (3-star) records, of ``CLNSIGCONF`` conflict
  breakdowns, and — importantly — of records that carry an *oncogenicity*
  classification (``ONC``) and no germline ``CLNSIG`` at all, which the adapter
  must treat as "no germline classification on record", not as benign.

Region order is 7, 15, 17 so each contig's records stay contiguous and sorted,
which is what tabix requires of an input it is asked to index.

What the slice deliberately CANNOT cover: multi-ALT records. A full scan of the
2026-08-22 release found 0 of them in 4,467,990 records — ClinVar writes one ALT
per line. The adapter still splits per-ALT because VCF permits it, and that path
is exercised by a hand-built VCF inside the test module rather than faked here.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pysam

#: (contig, 1-based inclusive start, 1-based inclusive end). Bare ClinVar contigs.
REGIONS: tuple[tuple[str, int, int], ...] = (
    ("7", 117_509_000, 117_510_000),
    ("15", 40_199_000, 40_206_000),
    ("17", 43_045_000, 43_048_000),
)

FIXTURE_NAME = "clinvar_slice.vcf.gz"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path, help="Path to clinvar.vcf.gz")
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()

    source: Path = args.source
    out_dir: Path = args.out_dir
    plain = out_dir / FIXTURE_NAME.removesuffix(".gz")
    compressed = out_dir / FIXTURE_NAME

    tabix = pysam.TabixFile(str(source))
    try:
        lines: list[str] = list(tabix.header)
        for contig, start, end in REGIONS:
            # pysam.fetch takes 0-based half-open coordinates.
            lines.extend(tabix.fetch(contig, start - 1, end))
    finally:
        tabix.close()

    plain.write_text("\n".join(lines) + "\n", encoding="utf-8")
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)
    plain.unlink()

    with gzip.open(compressed, "rt", encoding="utf-8") as handle:
        records = sum(1 for line in handle if not line.startswith("#"))
    print(f"{compressed}: {records} records, {compressed.stat().st_size} bytes")


if __name__ == "__main__":
    main()
