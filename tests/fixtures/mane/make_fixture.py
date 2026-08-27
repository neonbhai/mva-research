"""Regenerate the committed MANE slice used by ``tests/unit/test_gene_intervals.py``.

The full MANE GRCh38 v1.5 release (an 8.6 MB GTF over 486,796 lines, plus a
1.1 MB summary over 19,437 rows) is NOT committed. This script cuts a fixed,
documented set of regions out of both files so the gene-interval adapter's tests
run against *genuine* MANE records — real gene spans, real Ensembl and HGNC
identifiers, the real RefSeq-accession spelling of the chromosome, and the real
symbol drift — instead of rows hand-written to match the parser.

Reproduce (the inputs are the files the acquisition step downloads):

    uv run python tests/fixtures/mane/make_fixture.py \
        --gtf ~/Contri/bio-hackathon/mva-resources/mane/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz \
        --summary ~/Contri/bio-hackathon/mva-resources/mane/MANE.GRCh38.v1.5.summary.txt.gz

Both outputs keep the ``MANE.<build>.v<release>.`` filename prefix, because that
prefix is the only place a MANE distribution states its own release: the GTF
carries no comment or header line at all (verified: zero lines beginning ``#``),
and the summary's single header line names the columns, not the version. The
adapter parses the release from that prefix and cross-checks it against the
release name gzip stores *inside* the compressed bytes, so a fixture whose name
disagreed with its content would fail construction rather than mislabel itself.

Why each region, and what it exists to test
-------------------------------------------

Every region is a real, contiguous window of GRCh38. Genes are selected by
*overlap* with the window and then emitted in full, so no gene span is truncated.

* ``chr5:560,000-925,000`` — CEP72(+), TPPP(-), ZDHHC11B(-), ZDHHC11(-),
  BRD9(-), TRIP13(+).

  - **TRIP13**, an MVA gene, at its true GRCh38 span 892,849-919,357. Its GRCh37
    coordinate is around chr5:600,000, which this window also covers and which
    is *intergenic* in GRCh38 — 12 kb short of CEP72 and 293 kb short of TRIP13.
    A build mix-up there returns nothing, silently, which is exactly the failure
    the submission contract calls the most dangerous detail in the project.
  - **Overlapping genes**: ZDHHC11 (795,605-858,973) and BRD9 (850,291-892,838)
    share 8,683 bases.
  - **A 10 bp intergenic gap**: BRD9 ends at 892,838 and TRIP13 begins at
    892,849, so 892,839-892,848 belongs to no gene at all. Small enough that an
    off-by-one in either direction lands inside a gene and is caught.
  - **Minus-strand genes**: TPPP, ZDHHC11B, ZDHHC11 and BRD9.

* ``chr11:95,650,000-95,950,000`` — FAM76B(-), CEP57(+), MTMR2(-).

  **CEP57**, an MVA gene, overlaps *both* of its neighbours: FAM76B on its left
  (95,768,953-95,790,409 against CEP57's 95,789,965) and MTMR2 on its right
  (95,821,766 against CEP57's 95,837,070). A single variant at 95,790,000 sits
  in two genes at once, on opposite strands. Collapsing that to one gene is the
  data-loss bug this fixture exists to make visible.

* ``chr13:24,760,000-24,930,000`` — RNF17(+), CPAP(-).

  The **symbol-drift trap**. The gene a clinician still calls *CENPJ* is present
  in MANE v1.5 only under its current HGNC-approved symbol **CPAP**
  (ENSG00000151849.18, HGNC:17272, GeneID:55835). A panel that names ``CENPJ``
  resolves to nothing, and the adapter must say so rather than drop it.

* ``chr15:40,080,000-40,345,000`` — BMF(-), BUB1B(+), PAK6(+), PLCB2(-),
  ANKRD63(-), INAFM2(+), CCDC9B(-).

  - **BUB1B**, the MVA gene this project exists for, at 40,160,984-40,221,137 —
    note the *gene* span, which is wider than the 40,161,069-40,221,123 MANE
    Select *transcript* span, and is what an interval join must use if variants
    in the UTR flanks are not to be lost.
  - **Overlapping genes on the same strand**: BUB1B and PAK6 (40,217,428) share
    3,710 bases, so a variant in BUB1B's last exon is also in PAK6.
  - **A nested locus**: ANKRD63 (40,278,372-40,283,064) sits entirely inside
    PLCB2 (40,278,176-40,307,965). Nesting breaks a naive sorted-by-start scan
    that stops at the first non-overlapping interval.
  - **A 52 kb intergenic gap** between BMF (ends 40,108,928) and BUB1B.

* ``chrX:153,940,000-154,145,000`` — HCFC1(-), TMEM187(+), IRAK1(-), MECP2(-).

  - The **sex-chromosome contig mapping**: these rows say ``chrX`` in the GTF and
    ``NC_000023.11`` in the summary. Nothing in either spelling contains the
    letter X, so the mapping has to be stated explicitly and can be asserted.
  - **MECP2 carries two summary rows** — one ``MANE Select`` and one ``MANE Plus
    Clinical`` — for one Ensembl gene. Both are real transcripts and both are
    kept; a parser that keyed a dict on the gene would silently drop one.
  - Another overlapping pair (HCFC1/TMEM187) and a 922 bp intergenic gap between
    IRAK1 (ends 154,020,650) and MECP2 (begins 154,021,573).

Two extra genes are included by identifier rather than by region, because both
files must agree about which loci are *out of scope*:

* ``ENSG00000293532.1`` (MUC2) — GTF contig ``chr11_KQ759759v2_fix``, summary
  accession ``NW_015148966.2``.
* ``ENSG00000277656.3`` (GSTT1) — GTF contig ``chr22_KI270879v1_alt``, summary
  accession ``NT_187633.1``.

Neither is a canonical chromosome, so both must be skipped — and skipped
*consistently*, by both the chr-prefixed and the RefSeq-accession code path. A
build that dropped one file's copy but not the other's would leave a gene whose
two halves disagree about where it is, which is the shape of the contig bug this
whole module is defending against.
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

#: (contig as spelled in the GTF, 1-based inclusive start, 1-based inclusive end).
REGIONS: tuple[tuple[str, int, int], ...] = (
    ("chr5", 560_000, 925_000),
    ("chr11", 95_650_000, 95_950_000),
    ("chr13", 24_760_000, 24_930_000),
    ("chr15", 40_080_000, 40_345_000),
    ("chrX", 153_940_000, 154_145_000),
)

#: Genes pulled in by identifier rather than by region.
#:
#: The first two are off the primary assembly — one on a ``_fix`` patch scaffold,
#: one on an ``_alt`` haplotype — and must be skipped by both readers.
#:
#: The third is a 51 bp locus where the two files disagree about the gene's name:
#: the GTF calls it by its Ensembl clone name ``RP11-148K1.15`` and the summary by
#: its NCBI identifier ``LOC128092247``. 47 genes in v1.5 are like this, all of
#: them loci HGNC has not given an approved symbol, and its ``HGNC_ID`` cell is
#: empty. Picking either file's spelling would make the other unfindable, so both
#: resolve. Its 51 bases also make it the smallest interval the index holds.
EXTRA_GENE_IDS: tuple[str, ...] = (
    "ENSG00000293532.1",  # MUC2      chr11_KQ759759v2_fix / NW_015148966.2
    "ENSG00000277656.3",  # GSTT1     chr22_KI270879v1_alt / NT_187633.1
    "ENSG00000288608.1",  # chr7:151,061,928-151,061,978
)

GTF_NAME = "MANE.GRCh38.v1.5.slice.ensembl_genomic.gtf.gz"
SUMMARY_NAME = "MANE.GRCh38.v1.5.slice.summary.txt.gz"

_GENE_ID = re.compile(r'gene_id "([^"]+)"')


def _selected_gene_ids(gtf_path: Path) -> set[str]:
    """Gene identifiers whose ``gene`` row overlaps a region, plus the extras."""
    selected = set(EXTRA_GENE_IDS)
    with gzip.open(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            columns = line.split("\t")
            if len(columns) != 9 or columns[2] != "gene":
                continue
            start, end = int(columns[3]), int(columns[4])
            for contig, region_start, region_end in REGIONS:
                if columns[0] == contig and start <= region_end and region_start <= end:
                    match = _GENE_ID.search(columns[8])
                    if match is not None:
                        selected.add(match.group(1))
                    break
    return selected


def _write_gtf(gtf_path: Path, out_path: Path, gene_ids: set[str]) -> int:
    """Copy every line belonging to a selected gene, in the source's own order."""
    kept: list[str] = []
    with gzip.open(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            match = _GENE_ID.search(line)
            if match is not None and match.group(1) in gene_ids:
                kept.append(line)
    with gzip.open(out_path, "wt", encoding="utf-8") as out:
        out.writelines(kept)
    return len(kept)


def _write_summary(summary_path: Path, out_path: Path, gene_ids: set[str]) -> int:
    """Copy the header plus every summary row for a selected gene."""
    kept: list[str] = []
    with gzip.open(summary_path, "rt", encoding="utf-8") as handle:
        kept.append(next(handle))
        for line in handle:
            fields = line.split("\t")
            if len(fields) > 1 and fields[1] in gene_ids:
                kept.append(line)
    with gzip.open(out_path, "wt", encoding="utf-8") as out:
        out.writelines(kept)
    return len(kept) - 1


def main() -> None:
    parser = argparse.ArgumentParser(description="Cut the committed MANE test slice.")
    parser.add_argument("--gtf", required=True, type=Path, help="MANE.*.ensembl_genomic.gtf.gz")
    parser.add_argument("--summary", required=True, type=Path, help="MANE.*.summary.txt.gz")
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()

    gtf_path: Path = args.gtf
    summary_path: Path = args.summary
    out_dir: Path = args.out_dir

    gene_ids = _selected_gene_ids(gtf_path)
    gtf_out = out_dir / GTF_NAME
    summary_out = out_dir / SUMMARY_NAME

    gtf_lines = _write_gtf(gtf_path, gtf_out, gene_ids)
    summary_rows = _write_summary(summary_path, summary_out, gene_ids)

    print(f"selected {len(gene_ids)} genes")
    print(f"{gtf_out}: {gtf_lines} lines, {gtf_out.stat().st_size} bytes")
    print(f"{summary_out}: {summary_rows} rows, {summary_out.stat().st_size} bytes")


if __name__ == "__main__":
    main()
