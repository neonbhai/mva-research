"""Regenerate ``gnomad.exomes.v4.1.sites.chr21.slice.vcf.bgz`` from the real release.

The gnomAD v4.1 exomes chr21 sites VCF is 2,250,319,809 bytes and is NOT
committed. This script cuts a fixed, documented set of tiny windows out of it so
the adapter's tests run against genuine gnomAD records — real ``AC``/``AN``/``AF``
numbers for all nine genetic ancestry groups, real FILTER combinations, real
indel representations — instead of records someone hand-wrote to match the parser.

Reproduce (the input is the file the acquisition step downloads):

    uv run python tests/fixtures/gnomad/make_fixture.py \
        --source ~/Contri/bio-hackathon/mva-resources/gnomad/v4.1_exomes/\
gnomad.exomes.v4.1.sites.chr21.vcf.bgz

Equivalent with htslib on the PATH (this repo has no bgzip/tabix binary, which
is why the script uses pysam instead):

    { tabix -H "$SRC"; tabix "$SRC" chr21:5031903-5031916 ... ; } > slice.vcf
    bgzip -f slice.vcf && tabix -p vcf slice.vcf.gz

Why windows of a few bases each rather than one contiguous region: a gnomAD v4.1
exomes record carries 413 INFO fields and averages ~8 KB, so even a 6 kb region
is 21 MB. Fourteen windows totalling 154 records compress to ~150 KB, which is
committable, while still giving every shape the adapter has to get right.

The windows, and the record each was chosen for (all coordinates are public
gnomAD reference data, not patient data):

* ``5031903-5031916`` — ``chr21:5031905 C>A`` and ``C>T``: **AC=0 over
  AN=403,848**. A genuine observation of zero carriers in 400,000 chromosomes,
  which is the strongest rarity evidence gnomAD can give, and which must never
  be confused with "gnomAD has no record here" (GP-14). Both carry
  ``FILTER=AC0;AS_VQSR``, which is not incidental: a full scan of chr21 found
  **no** ``PASS`` record with ``AC=0``, because zero post-QC alleles is exactly
  what the ``AC0`` filter marks. An adapter that dropped filtered records would
  therefore discard every one of these, converting the strongest rarity signal
  in the dataset into silence. Two ALTs at one POS on separate lines also make
  "position match is not allele match" testable.
* ``5031933-5031937``, ``5031982-5031993`` — ``5031984 CA>C`` (deletion,
  AN=850,514, PASS) and, at one position, ``5031991 G>GA`` *and* ``5031991
  GA>G``: an insertion and a deletion anchored at the same base, with different
  frequencies. Nothing tests an indel join key harder.
* ``5033362-5033368``, ``5033393-5033397``, ``5033452-5033456`` — indels in
  homopolymer context (``5033364 A>AT``, ``5033365 C>CT``, ``5033395 GT>G``,
  ``5033454 CG>C``). These are where a left-alignment disagreement between the
  proband VCF and gnomAD shows up as a silent join failure.
* ``5035215-5035226`` — ``5035217 T>G``: ``AN=0``. gnomAD emits the record but
  no ``AF`` key at all, because no genotype survived QC. 33,596 of chr21's
  2,188,842 records are this shape. It is absence of information, not an
  observation of zero, and the adapter omits it.
* ``5035655-5035661`` — ``5035658 C>T``: a **common** variant, AF=0.336 over
  AN=759,570, PASS. The end of the frequency range opposite ``5031905``.
* ``5036075-5036081`` — ``5036078 A>C``: global AF=0.0037 but **afr AF=0.170**
  (AN=15,066) against **asj AF=0.0** (AN=4,528). The ADR 0010 case in real data:
  a variant that reads as rare globally and common in one ancestry group. If the
  per-population numbers are missing or misparsed, the population maximum
  silently collapses to the global figure and this variant looks rare.
* ``6086418-6086424`` — ``6086421 G>T``: ``FILTER=AC0;InbreedingCoeff``. Also
  ``6086419 GCGGAG>G``, a five-base deletion.
* ``9027126-9027130``, ``9027273-9027277`` — ``9027128 C>T``:
  ``FILTER=InbreedingCoeff`` on a record with AF=0.512 over AN=1,128,800, and
  ``9027275 G>A``: ``FILTER=AS_VQSR;InbreedingCoeff``. Between these windows and
  the ones above, all four of the release's FILTER IDs appear, alone and in
  every combination chr21 contains.
* ``9810548-9810552``, ``9810833-9810837`` — ``9810550 G>A`` (sas 0.216 vs eas
  0.00017) and ``9810835 GAC>G``, a **common deletion** with divergent
  per-population frequencies: the indel and the ancestry traps at once.

What the slice deliberately CANNOT cover: multi-ALT records. A full scan of
chr21 in this release found **0** in 2,188,842 records — gnomAD ships one ALT
per line. The adapter still splits per-ALT because VCF permits it and because
``Number=A`` INFO values would otherwise be attributed to the wrong allele; that
path is exercised by a hand-built VCF inside the test module rather than faked
here.
"""

from __future__ import annotations

import argparse
import gzip
import os
from pathlib import Path

import pysam

from mva.config import find_repo_root
from mva.privacy.audit import pinned_source


def _default_resource_root() -> Path | None:
    """``$MVA_RESOURCES``, expanded, if it is set.

    Where the acquisition tool puts its downloads is that tool's business
    (``tools.acquire.fetch.resolve_resource_root``); this only reads the same
    environment variable rather than growing a second opinion about the default.
    """
    raw = os.environ.get("MVA_RESOURCES")
    return Path(raw).expanduser() if raw else None


#: The manifest resource this fixture is cut from. Must match the ``resource:``
#: recorded for it in ``fixture_provenance.yaml``.
RESOURCE = "gnomad_exomes_chr21"

#: 1-based inclusive windows on ``chr21``, in ascending order. tabix requires a
#: position-sorted input, and emitting the windows in this order keeps it sorted
#: without a re-sort that could reorder records sharing a position.
REGIONS: tuple[tuple[int, int], ...] = (
    (5_031_903, 5_031_916),
    (5_031_933, 5_031_937),
    (5_031_982, 5_031_993),
    (5_033_362, 5_033_368),
    (5_033_393, 5_033_397),
    (5_033_452, 5_033_456),
    (5_035_215, 5_035_226),
    (5_035_655, 5_035_661),
    (5_036_075, 5_036_081),
    (6_086_418, 6_086_424),
    (9_027_126, 9_027_130),
    (9_027_273, 9_027_277),
    (9_810_548, 9_810_552),
    (9_810_833, 9_810_837),
)

CONTIG = "chr21"

#: Keeps ``GnomadSitesFrequencyAdapter``'s release/subset filename cross-check
#: meaningful on the fixture: the name still says which release it was cut from.
FIXTURE_NAME = "gnomad.exomes.v4.1.sites.chr21.slice.vcf.bgz"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        default=None,
        help=(
            "Path to the pinned release. Optional: resolved from "
            "knowledge/manifests/resources.yaml under --resource-root by default. "
            "Either way its sha256 must match the digest the manifest pins."
        ),
    )
    parser.add_argument(
        "--resource-root",
        type=Path,
        default=_default_resource_root(),
        help="External resource root holding the acquired releases. Defaults to $MVA_RESOURCES.",
    )
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()

    source = pinned_source(
        find_repo_root(),
        RESOURCE,
        source=args.source,
        resource_root=args.resource_root,
    )
    out_dir: Path = args.out_dir
    plain = out_dir / FIXTURE_NAME.removesuffix(".bgz")
    compressed = out_dir / FIXTURE_NAME

    tabix = pysam.TabixFile(str(source))
    try:
        lines: list[str] = list(tabix.header)
        seen: set[str] = set()
        for start, end in REGIONS:
            # pysam.fetch takes 0-based half-open coordinates.
            for line in tabix.fetch(CONTIG, start - 1, end):
                # Windows are disjoint, but a deletion anchored before a window
                # can be returned by it and by its neighbour; dedupe so the
                # fixture never carries a record twice.
                if line not in seen:
                    seen.add(line)
                    lines.append(line)
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
