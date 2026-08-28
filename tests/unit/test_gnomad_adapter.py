"""Unit tests for the real gnomAD sites-VCF frequency adapter.

Everything here runs against
``tests/fixtures/gnomad/gnomad.exomes.v4.1.sites.chr21.slice.vcf.bgz`` — a
154-record slice cut out of the genuine gnomAD v4.1 exomes chr21 sites VCF
(2,250,319,809 bytes) by ``tests/fixtures/gnomad/make_fixture.py``, which
documents the exact windows and why each was chosen. The records are real, so the
AC/AN/AF numbers, the FILTER combinations, the omission of ``AF_<grp>`` where a
group has no called alleles, and the indel representations are all the source's
own, not something written to match the parser.

A frequency adapter's interesting failures are not "wrong float returned". They
are:

* a variant gnomAD has never seen coming back at ``allele_frequency=0.0``, which
  manufactures an ultra-rare candidate out of nothing — worst in exactly the
  ancestries reference panels under-sample (GP-14);
* the mirror image: dropping ``AC0``-filtered records, which on this dataset
  deletes *every* genuine "zero carriers in 400,000 chromosomes" observation,
  because no ``PASS`` record on chr21 has ``AC=0``;
* per-population AC/AN going missing or landing on the wrong allele, which makes
  the ADR 0010 population maximum silently collapse to the global figure;
* an indel representation that disagrees with ``mva.ingestion.normalise``, which
  does not raise — it just fails to join, and a variant with no frequency record
  looks exactly like a rare one. This adapter shipped that defect: a private
  ``minimal_representation`` that trimmed but could not left-align, so a
  left-aligned proband indel and the gnomAD record for the same event were
  looked up under two different keys. ADR 0018 puts the rule in
  ``mva.alleles`` once; the tests under "normalisation / join" below are what
  keep it there.

Each of those has a test below. Tests that need the full 2 GB release are
skipped unless it is present; everything load-bearing runs off the committed
slice.
"""

from __future__ import annotations

import ast
import gzip
import os
import re
import time
import traceback
from collections.abc import Callable, Iterator, Mapping, Sequence
from datetime import UTC, datetime
from pathlib import Path
from typing import Protocol, cast

import pysam
import pytest

from mva.alleles import LeftAlignmentStatus, canonicalise_allele
from mva.annotation import (
    SYNTHETIC_STANDIN_LIMITATION,
    AdapterRole,
    AdapterSet,
    FrequencyAdapter,
    LocalConsequenceAdapter,
    annotate_variants,
    is_synthetic,
    load_default_adapters,
)
from mva.annotation.gnomad_sites import (
    ADAPTER_NAME,
    BGZF_EOF,
    GLOBAL_POPULATION,
    PASS_FILTER,
    UNREPRESENTABLE_GNOMAD_FACTS,
    FrequencyLookup,
    GnomadSitesFrequencyAdapter,
    check_source_complete,
    has_bgzf_eof,
    index_path_for,
    merge_query_regions,
)
from mva.clock import FixedClock
from mva.determinism import stable_hash
from mva.errors import AdapterUnavailableError, GenomeBuildMismatchError, NetworkDeniedError
from mva.ingestion.normalise import normalise_variants, trim_and_left_align
from mva.models.base import AssertionTier
from mva.models.genome import GenomeBuild, GenomicCoordinate
from mva.models.variant import (
    FilterStatus,
    Genotype,
    PopulationFrequency,
    VariantRecord,
    Zygosity,
)

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE = (
    REPO_ROOT / "tests" / "fixtures" / "gnomad" / "gnomad.exomes.v4.1.sites.chr21.slice.vcf.bgz"
)
KNOWLEDGE_ROOT = REPO_ROOT / "knowledge"
ADAPTER_SOURCE = REPO_ROOT / "src" / "mva" / "annotation" / "gnomad_sites.py"

RELEASE = "v4.1"
SUBSET = "exomes"

#: The full release the fixture was cut from. Present on a development machine,
#: absent in CI; the tests that need it skip rather than fail.
FULL_RELEASE_DIR = Path(
    os.environ.get(
        "MVA_GNOMAD_SITES_DIR",
        str(REPO_ROOT.parent / "mva-resources" / "gnomad" / "v4.1_exomes"),
    )
)
FULL_CHR21 = FULL_RELEASE_DIR / "gnomad.exomes.v4.1.sites.chr21.vcf.bgz"

requires_full_release = pytest.mark.skipif(
    not FULL_CHR21.is_file(),
    reason=f"full gnomAD chr21 release not present at {FULL_CHR21}",
)

#: The nine genetic ancestry groups gnomAD v4.1 **exomes** reports. Written down
#: so the header-driven discovery has something to be wrong against. Note the
#: absence of ``ami`` (Amish): that group exists only in the genomes callset, and
#: hard-coding a v4 group list copied from the genomes release would emit a
#: population this file has no numbers for.
EXOME_ANCESTRY_GROUPS = ("afr", "amr", "asj", "eas", "fin", "mid", "nfe", "remaining", "sas")

#: Names a re-introduced private canonicalisation rule would plausibly be given.
#: ``minimal_representation`` is the one that actually shipped here; the rest are
#: the spellings the same mistake wears elsewhere. The repo-wide version of this
#: lint, over every adapter, lives in ``tests/unit/test_normalise_representation.py``.
TRIMMING_FUNCTION_NAMES = frozenset(
    {
        "minimal_representation",
        "_minimal_representation",
        "parsimony_trim",
        "_parsimony_trim",
        "trim_alleles",
        "_trim_alleles",
        "normalise_allele",
        "_normalise_allele",
        "left_align",
        "_left_align",
        "left_shift",
        "_left_shift",
        "canonicalise_allele",
    }
)

# --------------------------------------------------------------------------- anchors
#
# Real records from the slice, each chosen for the trap it exercises. Coordinates
# are public gnomAD reference data, not patient data.

#: ``AC=0`` over ``AN=403,848``, with every one of the nine groups reporting
#: ``AC=0`` over a non-zero AN. A genuine observation of zero carriers in 400,000
#: chromosomes — the strongest rarity evidence gnomAD can give. FILTER is
#: ``AC0;AS_VQSR``, which is not incidental: see ``test_no_pass_record_carries_ac0``.
TRUE_ZERO = "GRCh38:chr21:5031905:C:A"

#: Same POS and REF as ``TRUE_ZERO``, a different ALT that gnomAD also holds.
TRUE_ZERO_SIBLING = "GRCh38:chr21:5031905:C:T"

#: Same POS as ``TRUE_ZERO``, an ALT gnomAD does **not** hold. Position match is
#: not allele match.
UNLISTED_ALT_AT_A_KNOWN_SITE = "GRCh38:chr21:5031905:C:G"

#: Common: AF 0.336 over AN 759,570, PASS. The far end of the range from ``TRUE_ZERO``.
COMMON = "GRCh38:chr21:5035658:C:T"

#: Global AF 0.0037 but ``afr`` AF 0.170 (AN 15,066) against ``asj`` AF 0.0
#: (AN 4,528). The ADR 0010 case in real data.
DIVERGENT = "GRCh38:chr21:5036078:A:C"

#: ``AN=0`` in every cohort: gnomAD emits the record and no ``AF`` key at all.
#: Absence of information, not an observation of zero.
NO_CALLED_ALLELES = "GRCh38:chr21:5035217:T:G"

#: ``AC0;InbreedingCoeff``. Six groups have ``AN=0`` (omitted) and three have
#: ``AN>0`` with ``AC=0`` (kept at 0.0) — both shapes inside one record.
MIXED_ABSENCE_AND_ZERO = "GRCh38:chr21:6086421:G:T"

#: One-base deletion, PASS, AN 850,514.
DELETION = "GRCh38:chr21:5031984:CA:C"

#: An insertion and a deletion anchored at the same base, with frequencies that
#: differ by a factor of 27. Nothing tests an indel join key harder.
INSERTION_AT_5031991 = "GRCh38:chr21:5031991:G:GA"
DELETION_AT_5031991 = "GRCh38:chr21:5031991:GA:G"

#: An insertion of ``C`` immediately before a ``CC`` run (the real bases at
#: 5033364-5033366 are ``A``, ``C``, ``C``), so the same event has three legal
#: right-shifted spellings. gnomAD holds the left-most one, at AF 0.337. This is
#: where a left-alignment disagreement between the proband VCF and gnomAD shows
#: up as a silent join failure.
REPEAT_INSERTION = "GRCh38:chr21:5033364:A:AC"

#: The same single-``C`` insertion written one and two bases to the right. Neither
#: is a key gnomAD holds, so both join only once the adapter can left-align them
#: back onto ``REPEAT_INSERTION`` — which is the whole of ADR 0018 in two strings.
REPEAT_INSERTION_SHIFTED_ONE = "GRCh38:chr21:5033365:C:CC"
REPEAT_INSERTION_SHIFTED_TWO = "GRCh38:chr21:5033366:C:CC"

#: A 19-base deletion with AC=0 over AN=1,195,784 — a long REF span the region
#: query has to reach the end of, and another genuine zero.
LONG_DELETION = "GRCh38:chr21:9027270:TGCTGGCACTGGCTCCACA:T"

#: A common deletion with strongly divergent per-population frequencies
#: (sas 0.171 vs eas 0.0048): the indel trap and the ancestry trap at once.
COMMON_DELETION = "GRCh38:chr21:9810835:GAC:G"

#: ``FILTER=InbreedingCoeff`` on a record at AF 0.512 over AN 1,128,800 — gnomAD
#: distrusts the call and definitely saw it.
INBREEDING_COEFF = "GRCh38:chr21:9027128:C:T"

#: ``FILTER=AS_VQSR;InbreedingCoeff``: two filters, rejoined in file order.
TWO_FILTERS = "GRCh38:chr21:9027275:G:A"

#: Inside a window the slice covers, but no gnomAD record at this allele.
NOT_IN_GNOMAD = "GRCh38:chr21:5031906:A:G"

#: Another miss, at a position the slice does not cover at all.
ANOTHER_ABSENT = "GRCh38:chr21:7000000:G:A"

#: A contig the slice's index does not hold.
ON_AN_UNCOVERED_CONTIG = "GRCh38:chr7:117509123:G:A"


# --------------------------------------------------------------------------- helpers


def mini_header(
    *,
    groups: Sequence[str] = ("afr", "nfe"),
    assembly: str = "gnomAD_GRCh38",
    contigs: Sequence[str] = ("chr21",),
    hail_version: str = "0.2.123-12ebb27db620",
) -> str:
    """A minimal gnomAD-shaped VCF header, for shapes the real release cannot supply."""
    lines = [
        "##fileformat=VCFv4.2",
        f"##hailversion={hail_version}",
        "##vep_version=v105",
        *(f"##contig=<ID={contig},length=46709983,assembly={assembly}>" for contig in contigs),
        '##FILTER=<ID=AC0,Description="Allele count is zero after filtering">',
        '##FILTER=<ID=AS_VQSR,Description="Failed VQSR">',
        '##FILTER=<ID=InbreedingCoeff,Description="Inbreeding coefficient < -0.3">',
        '##FILTER=<ID=PASS,Description="Passed all variant filters">',
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="x">',
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="x">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="x">',
        '##INFO=<ID=nhomalt,Number=A,Type=Integer,Description="x">',
        '##INFO=<ID=AF_grpmax,Number=A,Type=Float,Description="x">',
        '##INFO=<ID=AF_raw,Number=A,Type=Float,Description="x">',
        '##INFO=<ID=AF_XX,Number=A,Type=Float,Description="x">',
        '##INFO=<ID=AF_non_ukb,Number=A,Type=Float,Description="x">',
    ]
    for group in groups:
        lines += [
            f'##INFO=<ID=AC_{group},Number=A,Type=Integer,Description="x">',
            f'##INFO=<ID=AN_{group},Number=1,Type=Integer,Description="x">',
            f'##INFO=<ID=AF_{group},Number=A,Type=Float,Description="x">',
            f'##INFO=<ID=nhomalt_{group},Number=A,Type=Integer,Description="x">',
        ]
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    return "\n".join(lines) + "\n"


def write_sites_vcf(
    directory: Path,
    body: str,
    *,
    header: str | None = None,
    name: str = "gnomad.exomes.v4.1.sites.chr21.mini.vcf",
) -> Path:
    """bgzip + tabix a sites VCF built in-test, returning the compressed path.

    The default filename carries a real release token and the subset, because the
    adapter cross-checks both (the header has no release string to check against).
    """
    plain = directory / name
    plain.write_text((header if header is not None else mini_header()) + body, encoding="utf-8")
    compressed = directory / f"{name}.gz"
    pysam.tabix_compress(str(plain), str(compressed), force=True)
    pysam.tabix_index(str(compressed), preset="vcf", force=True)
    plain.unlink()
    return compressed


def open_adapter(path: Path, **kwargs: object) -> GnomadSitesFrequencyAdapter:
    kwargs.setdefault("release", RELEASE)
    kwargs.setdefault("subset", SUBSET)
    return GnomadSitesFrequencyAdapter(path, **kwargs)  # type: ignore[arg-type]


@pytest.fixture(scope="module")
def adapter() -> Iterator[GnomadSitesFrequencyAdapter]:
    instance = open_adapter(FIXTURE)
    yield instance
    instance.close()


@pytest.fixture(scope="module")
def aligning_adapter() -> Iterator[GnomadSitesFrequencyAdapter]:
    """The same slice, with the reference the adapter needs to left-align.

    Kept as a second fixture rather than replacing the first: the two states are
    both real deployments — a run configured with ``inputs.reference_fasta`` and
    one without — and the difference between them is the thing under test.
    """
    instance = open_adapter(FIXTURE, reference=_SliceReference.from_fixture())
    yield instance
    instance.close()


def fixture_rows() -> list[list[str]]:
    """Every record in the committed slice, as raw column lists."""
    rows: list[list[str]] = []
    with gzip.open(FIXTURE, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                rows.append(line.rstrip("\n").split("\t"))
    return rows


def fixture_header_lines() -> list[str]:
    with gzip.open(FIXTURE, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle if line.startswith("#")]


def info_of(columns: list[str]) -> dict[str, str]:
    """Split a VCF INFO column, by hand, so the adapter's parse has an independent check."""
    fields: dict[str, str] = {}
    for part in columns[7].split(";"):
        key, separator, value = part.partition("=")
        fields[key] = value if separator else ""
    return fields


def row_variant_id(columns: list[str]) -> str:
    return f"GRCh38:{columns[0]}:{columns[1]}:{columns[3]}:{columns[4]}"


def by_population(
    frequencies: Sequence[PopulationFrequency],
) -> dict[str, PopulationFrequency]:
    return {item.population: item for item in frequencies}


def make_record(variant_id: str) -> VariantRecord:
    """A minimal proband-shaped record at a canonical variant ID."""
    build, contig, position, ref, alt = variant_id.split(":")
    return VariantRecord(
        coordinate=GenomicCoordinate(
            build=GenomeBuild(build), contig=contig, position=int(position), ref=ref, alt=alt
        ),
        genotype=Genotype(
            zygosity=Zygosity.HET,
            genotype_string="0/1",
            depth=40,
            ref_reads=20,
            alt_reads=20,
            genotype_quality=99,
        ),
        filter_status=FilterStatus.PASS,
        source_artifact="tests/fixtures/gnomad/gnomad.exomes.v4.1.sites.chr21.slice.vcf.bgz",
    )


def dump(result: Mapping[str, tuple[PopulationFrequency, ...]]) -> str:
    """Canonical hash of a frequency mapping, for byte-identity comparison."""
    return stable_hash(
        {key: [item.model_dump(mode="json") for item in value] for key, value in result.items()}
    )


def real_reader_factory() -> Callable[[Path], _ReaderLike]:
    """The adapter's own reader factory, for spies that wrap it.

    Reached for in exactly one place. The spies below monkeypatch
    ``gnomad_sites._open_reader`` because counting *the adapter's* region queries
    is the only way to assert the batching claim from outside; wrapping pysam
    globally would also catch the fixture helpers.
    """
    import mva.annotation.gnomad_sites as module

    return cast(Callable[[Path], _ReaderLike], module._open_reader)


class _ReaderLike(Protocol):
    """The reader surface the adapter's ``_open_reader`` factory returns.

    Declared here rather than imported so the spies below are fully typed against
    a written-down contract instead of against ``object``.
    """

    @property
    def raw_header(self) -> str: ...

    def __call__(self, region: str) -> Iterator[object]: ...

    def close(self) -> None: ...


class _ConstantReference:
    """A ``ReferenceLookup`` over one hard-coded window, for left-alignment tests."""

    def __init__(self, contig: str, start: int, sequence: str) -> None:
        self._contig = contig
        self._start = start
        self._sequence = sequence

    def fetch(self, contig: str, start: int, end: int) -> str:
        if contig != self._contig:
            raise KeyError(contig)
        low = start - self._start
        high = end - self._start + 1
        if low < 0 or high > len(self._sequence):
            raise IndexError(start)
        return self._sequence[low:high]


class _SliceReference:
    """A ``ReferenceLookup`` whose bases are the slice's own REF columns.

    A real GRCh38 FASTA is 3 GB and absent from CI, and a synthetic one covering
    chr21 out to 9.8 Mb would be a 10 MB test artifact. Every base this adapter's
    left-alignment can need, though, is already in the fixture: a gnomAD record
    states the reference bases under its own REF span, so the 154 committed
    records supply the true GRCh38 base at every position they cover, from the
    release itself rather than from something written to make a test pass.

    Positions outside those spans raise, which is the honest answer and the one
    :func:`mva.alleles.canonicalise_allele` is built for — an unreadable base
    stops the shift rather than inventing one (GP-14). It therefore never claims
    a left-most position it could not prove.
    """

    def __init__(self, contig: str, bases: Mapping[int, str]) -> None:
        self._contig = contig
        self._bases = dict(bases)

    @classmethod
    def from_fixture(cls, contig: str = "chr21") -> _SliceReference:
        bases: dict[int, str] = {}
        for row in fixture_rows():
            start = int(row[1])
            for offset, base in enumerate(row[3]):
                bases[start + offset] = base
        return cls(contig, bases)

    @property
    def known_positions(self) -> int:
        return len(self._bases)

    def fetch(self, contig: str, start: int, end: int) -> str:
        if contig != self._contig:
            raise KeyError(contig)
        try:
            return "".join(self._bases[position] for position in range(start, end + 1))
        except KeyError as exc:
            raise IndexError(start) from exc


# --------------------------------------------------------------------------- premise


def test_fixture_is_a_real_gnomad_release_slice() -> None:
    """Guard the premise of every other test: this is genuine gnomAD v4.1 exomes."""
    header = fixture_header_lines()
    assert "##contig=<ID=chr21,length=46709983,assembly=gnomAD_GRCh38>" in header
    assert "##hailversion=0.2.123-12ebb27db620" in header
    assert any(
        line.startswith('##INFO=<ID=AF_afr,Number=A,Type=Float,Description="Alternate allele ')
        for line in header
    )
    assert len(fixture_rows()) == 154


def test_the_release_header_line_does_not_exist() -> None:
    """The reason ``release`` is a required constructor argument.

    gnomAD's sites VCF names Hail, VEP, dbSNP, GENCODE, CADD, REVEL, SpliceAI and
    a dozen other tools, and never itself. There is no ``##gnomad_version``, no
    ``##source``, no ``##fileDate``. Asserted rather than asserted-in-prose,
    because if a future release *does* add one, this test fails and the adapter
    should start reading it instead of trusting the caller.
    """
    meta = [line for line in fixture_header_lines() if line.startswith("##")]
    keys = {line[2:].split("=", 1)[0] for line in meta}
    assert not keys & {"source", "fileDate", "gnomad_version", "gnomAD_version", "release"}
    release_bearing = [
        line
        for line in meta
        if not line.startswith(("##INFO", "##FILTER", "##contig"))
        and "4.1" in line
        and "version" not in line.lower()
    ]
    assert release_bearing == []


def test_no_pass_record_carries_ac0() -> None:
    """Why dropping filtered records would delete the strongest rarity evidence.

    ``AC=0`` over a large ``AN`` is a real observation of zero carriers. On this
    dataset every such record is FILTERed, because zero post-QC alleles is exactly
    what the ``AC0`` filter marks — so "drop non-PASS" and "lose every genuine
    zero" are the same code change.
    """
    zeroes = [row for row in fixture_rows() if info_of(row).get("AC") == "0"]
    assert zeroes, "fixture drifted: no AC=0 records left"
    assert all(row[6] != "PASS" for row in zeroes)
    assert any(int(info_of(row).get("AN", "0")) > 100_000 for row in zeroes)


def test_the_slice_carries_every_shape_the_adapter_must_handle() -> None:
    """A fixture that quietly stopped covering a case would make its test vacuous."""
    rows = fixture_rows()
    filters = {row[6] for row in rows}
    assert filters >= {"PASS", "AC0", "AS_VQSR", "AC0;AS_VQSR", "InbreedingCoeff"}
    assert sum(1 for row in rows if len(row[3]) != len(row[4])) >= 30, "indels"
    no_af = [row for row in rows if "AF" not in info_of(row)]
    assert len(no_af) == 3, "AN=0 records"
    assert all(info_of(row)["AN"] == "0" for row in no_af)
    assert sum(1 for row in rows if float(info_of(row).get("AF", "0")) > 0.05) >= 1, "common"


def test_the_release_contains_no_multi_allelic_records() -> None:
    """gnomAD ships one ALT per line — 0 multi-ALT in chr21's 2,188,842 records.

    Recorded as a property of the data rather than assumed by the parser: the
    adapter still splits per-ALT (``test_multi_alt_records_are_split_per_allele``),
    because VCF permits several and a ``Number=A`` value attributed to the wrong
    allele would be a wrong frequency rather than a missing one.
    """
    assert all("," not in row[4] for row in fixture_rows())


# --------------------------------------------------------------------------- identity


def test_adapter_identity(adapter: GnomadSitesFrequencyAdapter) -> None:
    assert adapter.name == ADAPTER_NAME
    assert adapter.version == RELEASE
    assert adapter.build is GenomeBuild.GRCH38


def test_source_label_is_built_from_the_header_not_a_display_string(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """``gnomAD`` comes from ``assembly=gnomAD_GRCh38``; ``exomes`` is the declared subset."""
    assert adapter.source_label == "gnomAD_exomes"
    assert adapter.header_facts.dataset == "gnomAD"
    frequency = adapter.frequencies([COMMON])[COMMON][0]
    assert frequency.source == "gnomAD_exomes"
    assert frequency.version == RELEASE


def test_adapter_declares_itself_real(adapter: GnomadSitesFrequencyAdapter) -> None:
    """GP-20: ``is_synthetic`` fails closed, so ``synthetic = False`` is deliberate."""
    assert adapter.synthetic is False
    assert is_synthetic(adapter) is False


def test_satisfies_the_frequency_adapter_protocol(adapter: GnomadSitesFrequencyAdapter) -> None:
    assert isinstance(adapter, FrequencyAdapter)


def test_adapter_set_reports_the_frequency_slot_as_real(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """The descriptor the run manifest and report footers are built from."""
    adapters = AdapterSet(
        consequence=LocalConsequenceAdapter(
            KNOWLEDGE_ROOT / "public" / "consequences.tsv", version="synthetic-v0.0"
        ),
        frequency=adapter,
    )
    frequency = next(d for d in adapters.descriptors() if d.role is AdapterRole.FREQUENCY)
    assert frequency.synthetic is False
    assert frequency.label == f"{ADAPTER_NAME}@{RELEASE}"


def test_populations_are_read_from_the_header_not_hardcoded(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """Ten cohorts: ``global`` plus the nine v4.1 **exome** ancestry groups.

    ``ami`` is deliberately absent — the Amish group exists only in the genomes
    callset, so a hard-coded "gnomAD v4 groups" list copied from the genomes
    release would emit a population this file holds no numbers for. ``grpmax``
    (a derived maximum, which would double-count the group that produced it),
    ``raw`` (pre-QC), the ``XX``/``XY`` sex strata and the ``non_ukb`` subset are
    all excluded too, and all four have ``AF_`` INFO IDs in this header.
    """
    assert adapter.populations == (GLOBAL_POPULATION, *EXOME_ANCESTRY_GROUPS)
    assert "ami" not in adapter.populations
    for excluded in ("grpmax", "raw", "XX", "XY", "non_ukb"):
        assert excluded not in adapter.populations


def test_header_fingerprint_identifies_the_pipeline_run(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """The identity the *data* carries, since the header has no release string."""
    fingerprint = adapter.header_fingerprint
    assert fingerprint.startswith("gnomAD/GRCh38/")
    assert "hailversion=0.2.123-12ebb27db620" in fingerprint
    assert "vep_version=v105" in fingerprint
    assert adapter.header_facts.filter_ids == ("AC0", "AS_VQSR", "InbreedingCoeff", "PASS")


def test_unrepresentable_facts_are_written_down() -> None:
    """An undocumented drop is indistinguishable from a parser that never looked."""
    assert any("AN=0" in item for item in UNREPRESENTABLE_GNOMAD_FACTS)
    assert any("faf95" in item for item in UNREPRESENTABLE_GNOMAD_FACTS)
    assert any("grpmax" in item for item in UNREPRESENTABLE_GNOMAD_FACTS)


# ----------------------------------------------------------- GP-14: absence is absence


def test_a_variant_gnomad_has_never_seen_is_omitted_entirely(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """Not an empty tuple, not AF=0 — no key at all (GP-14).

    This is the single most important property of the adapter. A variant defaulted
    to ``allele_frequency=0.0`` is scored as maximally rare and promoted to the top
    of the candidate list, and it happens most often where the reference panel is
    thinnest — turning under-representation into false ultra-rare candidates.
    """
    result = adapter.frequencies([NOT_IN_GNOMAD, ANOTHER_ABSENT, COMMON])
    assert NOT_IN_GNOMAD not in result
    assert ANOTHER_ABSENT not in result
    assert set(result) == {COMMON}
    assert all(value for value in result.values()), "no key may map to an empty tuple"


def test_a_genuine_zero_is_present_with_allele_frequency_zero(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """The counterpart, and the reason the test above cannot be satisfied by "return nothing".

    ``chr21:5031905 C>A`` is AC=0 over AN=403,848: gnomAD looked at 400,000
    chromosomes and found no carriers. That is the strongest rarity evidence the
    dataset can give and it must arrive, with its allele number attached, so that
    a reader can tell it apart from silence.
    """
    frequencies = adapter.frequencies([TRUE_ZERO])[TRUE_ZERO]
    populations = by_population(frequencies)
    assert populations[GLOBAL_POPULATION].allele_frequency == 0.0
    assert populations[GLOBAL_POPULATION].allele_count == 0
    assert populations[GLOBAL_POPULATION].allele_number == 403_848
    # Every one of the nine groups reports a real cohort here, all with AC=0.
    assert set(populations) == {GLOBAL_POPULATION, *EXOME_ANCESTRY_GROUPS}
    assert all(item.allele_frequency == 0.0 for item in frequencies)
    assert all((item.allele_number or 0) > 0 for item in frequencies)


def test_absence_and_a_genuine_zero_are_distinguishable_in_one_call(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """The whole adapter in one assertion: same call, two different facts."""
    result = adapter.frequencies([NOT_IN_GNOMAD, TRUE_ZERO])
    assert NOT_IN_GNOMAD not in result
    assert result[TRUE_ZERO][0].allele_frequency == 0.0


def test_an_entirely_unknown_batch_returns_an_empty_mapping(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    assert adapter.frequencies([NOT_IN_GNOMAD, ANOTHER_ABSENT]) == {}


def test_no_absent_variant_is_given_a_fabricated_frequency(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """The failure this test exists for would look like a helpful default."""
    for variant_id, frequencies in adapter.frequencies([NOT_IN_GNOMAD, ANOTHER_ABSENT]).items():
        pytest.fail(f"absent variant {variant_id} was given {len(frequencies)} frequency row(s)")


def test_a_group_with_no_called_alleles_is_omitted_while_a_genuine_zero_is_kept(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """The same distinction one level down, in a single real record.

    ``chr21:6086421 G>T`` has ``AN=0`` for six groups — gnomAD omits their
    ``AF_<grp>`` key entirely — and ``AN>0`` with ``AC=0`` for three. Emitting the
    first six at 0.0 would invent six "no carriers in this ancestry" observations
    out of six cohorts that were never sampled here.
    """
    frequencies = adapter.frequencies([MIXED_ABSENCE_AND_ZERO])[MIXED_ABSENCE_AND_ZERO]
    populations = by_population(frequencies)
    assert set(populations) == {GLOBAL_POPULATION, "fin", "nfe", "sas"}
    for missing in ("afr", "amr", "asj", "eas", "mid", "remaining"):
        assert missing not in populations
    assert populations["nfe"].allele_frequency == 0.0
    assert populations["nfe"].allele_number == 22
    assert populations[GLOBAL_POPULATION].allele_number == 26


def test_a_record_with_no_called_alleles_anywhere_is_omitted(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """``AN=0`` in every cohort yields no frequency at all, not a zero.

    33,596 of chr21's 2,188,842 records are this shape. gnomAD attempted the site
    and retained no genotype, so there is no allele frequency to report — and
    ``PopulationFrequency.allele_frequency`` is a required float, so "observed,
    frequency unknown" is not representable. Collapsing it into absence is the
    conservative direction; emitting 0.0 would not be.
    """
    row = next(r for r in fixture_rows() if row_variant_id(r) == NO_CALLED_ALLELES)
    assert info_of(row)["AN"] == "0"
    assert "AF" not in info_of(row), "fixture drifted: this record gained an AF"
    assert adapter.frequencies([NO_CALLED_ALLELES]) == {}


def test_a_query_on_a_contig_with_no_shard_fails_closed(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """ "We could not look" must never be rendered as "we looked and found nothing".

    This is GP-14's core failure mode and it is a live one, not a hypothesis. A
    gnomAD exomes release is ~250 GB of shards that arrive over hours, so a run
    started mid-download has no shard for most chromosomes. If a missing shard
    produced a missing key, a **common** variant on that chromosome would look
    exactly like a variant gnomAD has never seen: it would evade the
    common-frequency down-rank and be promoted as a rare candidate.

    The coverage hole is therefore an error, not an omission. The message names
    the chromosome — that is what the operator has to act on, it is a property of
    the release rather than of the proband, and one of twenty-four chromosome
    names discloses nothing — but never a position.
    """
    assert adapter.contig_map == {"chr21": "chr21"}
    with pytest.raises(AdapterUnavailableError) as excinfo:
        adapter.frequencies([ON_AN_UNCOVERED_CONTIG, COMMON])
    message = str(excinfo.value)
    assert "chr7" in message
    assert "117509123" not in message, "PRIV-09: the position must not be echoed"
    assert "lookup_partial" in message, "the message must name the deliberate way through"


def test_lookup_partial_returns_the_gap_in_a_typed_result(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """The explicit opt-in to an incomplete release.

    The gap is a field of the result, not a property of the adapter, precisely
    because a property is something a caller can forget to read — and forgetting
    it here reintroduces the bug above.
    """
    result = adapter.lookup_partial([ON_AN_UNCOVERED_CONTIG, COMMON, NOT_IN_GNOMAD])
    assert not result.is_complete
    assert result.unqueryable_contigs == ("chr7",)
    assert result.unqueryable_variant_ids == (ON_AN_UNCOVERED_CONTIG,)
    assert set(result.frequencies) == {COMMON}
    # NOT_IN_GNOMAD was looked up and genuinely has no record; the chr7 variant was
    # never looked up at all. Two different facts, and they stay different.
    assert NOT_IN_GNOMAD not in result.unqueryable_variant_ids
    assert "must not be scored as rarity" in result.describe_gap()


def test_lookup_partial_reports_completeness_when_there_is_no_gap(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    result = adapter.lookup_partial([COMMON, NOT_IN_GNOMAD])
    assert result.is_complete
    assert result.unqueryable_contigs == ()
    assert result.unqueryable_variant_ids == ()
    assert result.frequencies == adapter.frequencies([COMMON, NOT_IN_GNOMAD])
    assert "Every queried contig" in result.describe_gap()


def test_position_match_is_not_allele_match(adapter: GnomadSitesFrequencyAdapter) -> None:
    """The join key is all five fields. A neighbouring ALT must not answer for this one."""
    result = adapter.frequencies([TRUE_ZERO, TRUE_ZERO_SIBLING, UNLISTED_ALT_AT_A_KNOWN_SITE])
    assert set(result) == {TRUE_ZERO, TRUE_ZERO_SIBLING}
    assert UNLISTED_ALT_AT_A_KNOWN_SITE not in result


# ------------------------------------------------------------------------ PRIV-05


@pytest.mark.parametrize(
    "url",
    [
        "https://storage.googleapis.com/gcp-public-data--gnomad/x.vcf.bgz",
        "http://example.org/gnomad.exomes.v4.1.sites.chr21.vcf.bgz",
        "s3://gnomad-public/release/4.1/vcf/exomes/x.vcf.bgz",
        "gs://gcp-public-data--gnomad/release/4.1/vcf/exomes/x.vcf.bgz",
        "ftp://ftp.example.org/gnomad.exomes.v4.1.sites.chr21.vcf.bgz",
    ],
)
def test_a_remote_source_is_refused_before_htslib_sees_it(
    url: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """PRIV-05. htslib *can* range-request a remote tabix index; this adapter cannot.

    The guard runs on the raw text and matches a single-slash scheme on purpose:
    ``Path("https://host/x")`` stringifies as ``"https:/host/x"``, so a check
    written against ``"://"`` never fires on a caller who wrapped a URL in a
    ``Path`` — which is exactly what this adapter's signature invites.

    The backends are monkeypatched to explode, so this proves the refusal happens
    *before* any coordinate could reach htslib, not merely that construction fails.
    """
    import mva.annotation.gnomad_sites as module

    def forbidden(path: Path) -> object:
        pytest.fail(f"htslib was handed {path}")

    monkeypatch.setattr(module, "_open_reader", forbidden)
    monkeypatch.setattr(module, "_open_tabix", forbidden)

    with pytest.raises(NetworkDeniedError, match="local files only"):
        open_adapter(Path(url))


def test_the_remote_refusal_does_not_echo_the_path() -> None:
    """PRIV-09: a URL can carry a coordinate in its query string."""
    url = "https://example.org/lookup?chrom=chr21&pos=5031905"
    with pytest.raises(NetworkDeniedError) as excinfo:
        open_adapter(Path(url))
    assert "5031905" not in str(excinfo.value)
    assert "https" in str(excinfo.value)


def test_the_module_imports_no_network_client() -> None:
    """The structural half of PRIV-05, scoped to this file.

    ``tests/unit/test_architecture.py`` enforces this across the whole package;
    duplicating it here means a change to *this* module is caught by *this*
    module's suite rather than by a distant test a reviewer may not run.
    """
    forbidden = {"requests", "httpx", "urllib", "aiohttp", "http", "ftplib", "smtplib", "socket"}
    tree = ast.parse(ADAPTER_SOURCE.read_text(encoding="utf-8"), filename=str(ADAPTER_SOURCE))
    imported: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
            imported.add(node.module.split(".")[0])
    assert not imported & forbidden


# ------------------------------------------------- per-population AC/AN (ADR 0010)


def test_every_population_number_matches_the_vcf_info_exactly(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """The whole slice, field by field, against an independent hand parse.

    Not a spot check: 154 real records, 10 cohorts each, 4 numbers each. A per-allele
    index bug, an ``AN``/``AC`` transposition or a suffix typo would have to be
    consistent across all of them to survive this.
    """
    rows = fixture_rows()
    result = adapter.frequencies([row_variant_id(row) for row in rows])
    checked = 0
    for row in rows:
        info = info_of(row)
        variant_id = row_variant_id(row)
        if "AF" not in info and all(f"AF_{g}" not in info for g in EXOME_ANCESTRY_GROUPS):
            assert variant_id not in result
            continue
        populations = by_population(result[variant_id])
        for label, suffix in [(GLOBAL_POPULATION, "")] + [
            (g, f"_{g}") for g in EXOME_ANCESTRY_GROUPS
        ]:
            raw_af = info.get(f"AF{suffix}")
            raw_an = info.get(f"AN{suffix}")
            if raw_af is None or raw_an == "0":
                assert label not in populations
                continue
            observed = populations[label]
            assert observed.allele_frequency == pytest.approx(float(raw_af), rel=1e-6, abs=1e-12)
            assert observed.allele_number == int(raw_an or 0)
            assert observed.allele_count == int(info[f"AC{suffix}"])
            assert observed.homozygote_count == int(info[f"nhomalt{suffix}"])
            checked += 1
    assert checked > 1_000, f"only {checked} population rows compared"


def test_divergent_population_frequencies_survive(adapter: GnomadSitesFrequencyAdapter) -> None:
    """The ADR 0010 case in real data, and what a missing per-population parse costs.

    ``chr21:5036078 A>C`` reads as rare on the global figure (AF 0.0037) and is
    **common in one ancestry** (afr 0.170 over 15,066 alleles). An adapter that
    emitted only the global row would score this variant as a plausible rare
    candidate; the population maximum exists precisely so it does not.
    """
    populations = by_population(adapter.frequencies([DIVERGENT])[DIVERGENT])
    assert populations[GLOBAL_POPULATION].allele_frequency == pytest.approx(0.00373752, rel=1e-6)
    assert populations["afr"].allele_frequency == pytest.approx(0.169919, rel=1e-6)
    assert populations["afr"].allele_count == 2560
    assert populations["afr"].allele_number == 15_066
    # asj reports a genuine zero over a real cohort, which is not the same as absence.
    assert populations["asj"].allele_frequency == 0.0
    assert populations["asj"].allele_number == 4528
    assert populations["nfe"].allele_frequency == pytest.approx(5.45737e-05, rel=1e-6)


def test_the_allele_number_guard_reads_the_parsed_numbers(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """ADR 0010 end to end: the guard is inert unless per-group AN is right.

    With the configured ``min_allele_number`` of 2000 the ``afr`` row sets the
    maximum, and the variant is correctly seen as common-in-an-ancestry rather
    than rare. Raising the threshold above ``afr``'s cohort moves the maximum,
    which is only observable because the allele numbers are real.
    """
    record = make_record(DIVERGENT).with_annotations(
        population_frequencies=adapter.frequencies([DIVERGENT])[DIVERGENT]
    )
    selection = record.select_max_allele_frequency(min_allele_number=2000)
    assert selection.observed is not None
    assert selection.observed.population == "afr"
    assert selection.observed.allele_frequency == pytest.approx(0.169919, rel=1e-6)
    # Every cohort here clears 2000 except mid (1328) and fin (266).
    assert {item.population for item in selection.excluded} == {"amr", "fin", "mid"}

    strict = record.select_max_allele_frequency(min_allele_number=20_000)
    assert strict.observed is not None
    assert strict.observed.population != "afr"


def test_homozygote_counts_are_carried(adapter: GnomadSitesFrequencyAdapter) -> None:
    """``nhomalt`` is ``Number=A``; a per-allele slip would land it on the wrong allele."""
    populations = by_population(adapter.frequencies([COMMON])[COMMON])
    assert populations[GLOBAL_POPULATION].homozygote_count == 79_929
    assert populations["eas"].homozygote_count == 1003
    assert populations["eas"].allele_frequency == pytest.approx(0.771831, rel=1e-6)


def test_populations_are_emitted_in_a_fixed_order(adapter: GnomadSitesFrequencyAdapter) -> None:
    """GP-30: global first, then ancestry groups A-Z. Never set-iteration order."""
    order = [item.population for item in adapter.frequencies([COMMON])[COMMON]]
    assert order == [GLOBAL_POPULATION, *EXOME_ANCESTRY_GROUPS]


def test_a_missing_allele_number_leaves_a_population_eligible(tmp_path: Path) -> None:
    """ADR 0010: an unreported cohort size is unknown, not small (GP-14).

    Not a shape gnomAD emits — every ``AF_<grp>`` in the release is accompanied by
    an ``AN_<grp>`` — so it is built here rather than faked into the slice.
    """
    path = write_sites_vcf(
        tmp_path,
        "chr21\t100\t.\tA\tG\t.\tPASS\tAC=1;AF=0.01;AF_nfe=0.02;AC_nfe=1\n",
    )
    with open_adapter(path) as instance:
        populations = by_population(
            instance.frequencies(["GRCh38:chr21:100:A:G"])["GRCh38:chr21:100:A:G"]
        )
    assert populations[GLOBAL_POPULATION].allele_number is None
    assert populations["nfe"].allele_number is None
    assert populations["nfe"].allele_frequency == pytest.approx(0.02, rel=1e-6)


# --------------------------------------------------------------------------- FILTER


def test_filter_status_is_carried_verbatim(adapter: GnomadSitesFrequencyAdapter) -> None:
    """gnomAD's own string, including multi-filter combinations, in file order."""
    result = adapter.frequencies([COMMON, TRUE_ZERO, INBREEDING_COEFF, TWO_FILTERS])
    assert result[COMMON][0].filter_status == PASS_FILTER
    assert result[TRUE_ZERO][0].filter_status == "AC0;AS_VQSR"
    assert result[INBREEDING_COEFF][0].filter_status == "InbreedingCoeff"
    assert result[TWO_FILTERS][0].filter_status == "AS_VQSR;InbreedingCoeff"
    # The status is a property of the record, so every population row carries it.
    assert {item.filter_status for item in result[TRUE_ZERO]} == {"AC0;AS_VQSR"}


def test_filtered_records_are_not_dropped(adapter: GnomadSitesFrequencyAdapter) -> None:
    """ "gnomAD distrusted this call" must not become "gnomAD never saw it".

    ``chr21:9027128 C>T`` is ``InbreedingCoeff``-filtered at AF 0.512 over
    AN 1,128,800. Dropping it would make a variant present in half of a million
    chromosomes look novel.
    """
    frequencies = adapter.frequencies([INBREEDING_COEFF])[INBREEDING_COEFF]
    assert frequencies[0].allele_frequency == pytest.approx(0.512108, rel=1e-6)
    assert frequencies[0].allele_number == 1_128_800
    assert frequencies[0].filter_status == "InbreedingCoeff"


def test_every_filter_id_the_release_declares_reaches_a_result(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """All four declared FILTER IDs actually appear in emitted frequencies."""
    rows = fixture_rows()
    result = adapter.frequencies([row_variant_id(row) for row in rows])
    seen: set[str] = set()
    for frequencies in result.values():
        status = frequencies[0].filter_status
        if status:
            seen.update(status.split(";"))
    assert seen == set(adapter.header_facts.filter_ids)


def test_an_unfiltered_record_is_not_reported_as_pass(tmp_path: Path) -> None:
    """``.`` means "no filtering applied"; ``PASS`` means "passed all filters".

    cyvcf2 renders both as ``FILTER = None``, so an adapter reading the scalar
    writes ``"PASS"`` for a record whose file made no quality claim at all —
    inventing a judgement the source never made. Reading ``FILTERS`` keeps them
    apart, and ``filter_status=None`` is the honest rendering of ``.``.
    """
    path = write_sites_vcf(
        tmp_path,
        "chr21\t100\t.\tA\tG\t.\t.\tAC=1;AN=100;AF=0.01\n"
        "chr21\t200\t.\tA\tG\t.\tPASS\tAC=1;AN=100;AF=0.01\n",
    )
    with open_adapter(path) as instance:
        result = instance.frequencies(["GRCh38:chr21:100:A:G", "GRCh38:chr21:200:A:G"])
    assert result["GRCh38:chr21:100:A:G"][0].filter_status is None
    assert result["GRCh38:chr21:200:A:G"][0].filter_status == PASS_FILTER


# ------------------------------------------------------------- normalisation / join


def test_the_adapter_owns_no_canonicalisation_of_its_own() -> None:
    """ADR 0018, asserted against the module rather than trusted (GP-03 does not force a copy).

    This adapter used to define ``minimal_representation``: a private trim that,
    unlike the shared rule, could not left-align. The two agreed on every allele
    they were ever tested with and disagreed on the one class that matters — an
    indel inside a repeat tract — so the join failed silently and the variant was
    scored as novel and ultra-rare. The defect was not the algorithm; it was the
    second copy.

    ``mva.alleles`` sits at layer 1, below both ``ingestion`` and ``annotation``,
    so sharing it costs no GP-03 violation and there is nothing left to justify a
    duplicate.
    """
    import mva.annotation.gnomad_sites as module

    assert not hasattr(module, "minimal_representation")
    assert module.canonicalise_allele is canonicalise_allele

    tree = ast.parse(ADAPTER_SOURCE.read_text(encoding="utf-8"), filename=str(ADAPTER_SOURCE))
    local_functions = {
        node.name
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    assert not local_functions & TRIMMING_FUNCTION_NAMES, (
        "gnomad_sites.py defines its own trimming/canonicalisation function: "
        f"{sorted(local_functions & TRIMMING_FUNCTION_NAMES)}. Call "
        "mva.alleles.canonicalise_allele instead."
    )


def test_the_adapter_and_ingestion_agree_on_every_real_fixture_allele() -> None:
    """The shared rule, exercised over the 154 alleles this adapter will actually see.

    Both sides are asked in both reference states, because the two callers can
    only disagree in the state where one of them can move a coordinate. Comparing
    functions rather than comparing a join result is deliberate: a passing join is
    exactly the evidence that was available while the two implementations
    disagreed.
    """
    reference = _SliceReference.from_fixture()
    with open_adapter(FIXTURE) as plain, open_adapter(FIXTURE, reference=reference) as aligning:
        for lookup, instance in ((None, plain), (reference, aligning)):
            for row in fixture_rows():
                position, ref, alt = int(row[1]), row[3], row[4]
                from_adapter = instance.canonicalise("chr21", position, ref, alt)
                from_ingestion = trim_and_left_align(
                    make_record(f"GRCh38:chr21:{position}:{ref}:{alt}"), lookup
                ).coordinate
                assert (from_adapter.position, from_adapter.ref, from_adapter.alt) == (
                    from_ingestion.position,
                    from_ingestion.ref,
                    from_ingestion.alt,
                )


def test_representation_status_is_typed_and_names_its_own_limitation(
    adapter: GnomadSitesFrequencyAdapter,
    aligning_adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """GP-14: "we could not left-align" is a state a caller receives, not a silent skip.

    Without it, the repeat-tract miss below is indistinguishable from "gnomAD has
    no record", and a report has nothing to put in its footer.
    """
    assert adapter.representation_status is LeftAlignmentStatus.UNAVAILABLE_NO_REFERENCE
    assert aligning_adapter.representation_status is LeftAlignmentStatus.APPLIED
    assert aligning_adapter.representation_limitation is None

    limitation = adapter.representation_limitation
    assert limitation is not None
    assert "not evidence of rarity" in limitation
    # PRIV-09: the sentence reaches a report footer, a log and an agent's context.
    for forbidden in ("chr21", "5033364", "GRCh38", "A:AC", "0/1"):
        assert forbidden not in limitation


def test_real_indels_join_after_the_ingestion_normaliser_has_run(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """The join tested through the pipeline's own normaliser, not by hand.

    Every indel in the slice is turned into a ``VariantRecord``, pushed through
    ``normalise_variants`` (no reference, so trim only — the common configuration),
    and looked up by the ``variant_id`` that comes out. A representation mismatch
    would show up here as a miss, which is the failure mode that looks like
    "novel variant".
    """
    indel_rows = [row for row in fixture_rows() if len(row[3]) != len(row[4])]
    assert len(indel_rows) >= 30
    records = [make_record(row_variant_id(row)) for row in indel_rows]
    normalised = normalise_variants(records)
    ids = [record.coordinate.variant_id for record in normalised.variants]
    result = adapter.frequencies(ids)

    expected = {
        row_variant_id(row)
        for row in indel_rows
        if "AF" in info_of(row) or any(f"AF_{g}" in info_of(row) for g in EXOME_ANCESTRY_GROUPS)
    }
    assert set(result) == expected
    assert len(expected) >= 25


def test_an_insertion_and_a_deletion_at_one_position_do_not_cross_join(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """``chr21:5031991 G>GA`` and ``GA>G`` differ 27-fold. Confusing them is silent."""
    result = adapter.frequencies([INSERTION_AT_5031991, DELETION_AT_5031991])
    insertion = by_population(result[INSERTION_AT_5031991])
    deletion = by_population(result[DELETION_AT_5031991])
    assert insertion[GLOBAL_POPULATION].allele_frequency == pytest.approx(8.70054e-05, rel=1e-6)
    assert deletion[GLOBAL_POPULATION].allele_frequency == pytest.approx(0.00236422, rel=1e-6)
    assert insertion["afr"].allele_count == 24
    assert deletion["afr"].allele_count == 1359


def test_an_untrimmed_representation_of_the_same_indel_still_joins(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """Both sides of the join go through one trim rule, so padding is harmless.

    ``chr21:5031991 GAA>GA`` is the same deletion as ``GA>G`` with an extra shared
    base. A caller that has not run the normaliser still gets the right record.
    """
    padded = "GRCh38:chr21:5031991:GAA:GA"
    result = adapter.frequencies([padded, DELETION_AT_5031991])
    assert result[padded] == result[DELETION_AT_5031991]


def test_a_right_shifted_indel_joins_when_a_reference_is_supplied(
    aligning_adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """**The join this adapter used to lose.** Same event, three spellings, one answer.

    The real bases at chr21:5033364-5033366 are ``A``, ``C``, ``C``, so the
    single-``C`` insertion gnomAD stores as ``5033364 A>AC`` — AF 0.337, one
    chromosome in three — has two other legal spellings. Under the adapter's old
    private ``minimal_representation`` neither of them joined: trimming cannot undo
    a right shift, undoing it needs the reference bases to the *left*, and the
    result was that a variant a third of the population carries came back with no
    frequency data and was scored as novel and ultra-rare.

    Given the same ``ReferenceLookup`` ingestion and the ClinVar adapter already
    take, all three spellings canonicalise to the key gnomAD holds and all three
    get the same frequency. The assertion is on the value, not merely on
    membership: a join that returned *a* record would be just as wrong as no join
    if it returned the neighbouring allele's.
    """
    bases = {
        int(row[1]): row[3][0]
        for row in fixture_rows()
        if int(row[1]) in {5_033_364, 5_033_365, 5_033_366}
    }
    assert bases == {5_033_364: "A", 5_033_365: "C", 5_033_366: "C"}, "the repeat context"

    spellings = [REPEAT_INSERTION, REPEAT_INSERTION_SHIFTED_ONE, REPEAT_INSERTION_SHIFTED_TWO]
    result = aligning_adapter.frequencies(spellings)

    assert set(result) == set(spellings), "every legal spelling of one event must be answered"
    for spelling in spellings:
        assert by_population(result[spelling])[GLOBAL_POPULATION].allele_frequency == pytest.approx(
            0.337086, rel=1e-6
        )
    # Keyed by what the caller passed, not by what canonicalisation turned it into.
    assert result[REPEAT_INSERTION_SHIFTED_TWO] == result[REPEAT_INSERTION]


def test_without_a_reference_the_same_indel_is_an_honest_miss(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """The degraded state, kept explicit rather than quietly fixed away.

    Left-alignment is impossible without the reference genome, and an adapter that
    pretended otherwise would be moving a coordinate on a guess. What must not
    happen is for the resulting miss to look like a clean answer: the variant is
    omitted (never ``allele_frequency=0.0``), and ``representation_status`` says
    the omission may be representational rather than biological.

    This is the half of the old ``test_a_right_shifted_indel_does_not_join`` worth
    keeping. The other half — that a supplied reference could not fix it — was not
    a property of the problem, only of the missing parameter.
    """
    assert adapter.representation_status is LeftAlignmentStatus.UNAVAILABLE_NO_REFERENCE

    # The left-most spelling still joins: trimming is unconditional.
    joined = adapter.frequencies([REPEAT_INSERTION])
    assert joined[REPEAT_INSERTION][0].allele_frequency == pytest.approx(0.337086, rel=1e-6)

    shifted = [REPEAT_INSERTION_SHIFTED_ONE, REPEAT_INSERTION_SHIFTED_TWO]
    assert adapter.frequencies(shifted) == {}, "omitted, not zero-filled (GP-14)"


def test_the_region_window_reaches_a_record_the_release_spells_further_right(
    tmp_path: Path,
) -> None:
    """Fixing the key is half the fix; the fetch window has to cover both spellings.

    gnomAD's own release is left-aligned, so this shape is built rather than found
    — but the failure it guards against is not hypothetical, it is the same
    silent miss one layer down. A query left-aligned to POS 300 and a record
    spelled at POS 302 occupy disjoint spans, so a window of ``[300, 300]`` never
    fetches the record, never canonicalises it, and reports "gnomAD has no record"
    with a join key that was by then perfectly correct.

    The window runs from the left-most position out to
    ``rightmost_equivalent_position``, which is read from the reference rather than
    guessed as a padding constant. Asserted with ``merge_window_bp=0`` so the
    coalescing cannot supply the coverage by accident and hide a broken bound.
    """
    # chr21:296-308 = T G A T A C C C G T T A G, so position 300 is an ``A``
    # followed by a three-base ``C`` run: a single-``C`` insertion after it is
    # legal at 300, 301, 302 and 303 and every spelling is the same event.
    reference = _ConstantReference("chr21", 296, "TGATACCCGTTAG")
    body = "chr21\t302\t.\tC\tCC\t.\tPASS\tAC=10;AN=1000;AF=0.01;AC_afr=10;AN_afr=500;AF_afr=0.02\n"
    path = write_sites_vcf(tmp_path, body)

    left_aligned_query = "GRCh38:chr21:300:A:AC"
    with open_adapter(path, reference=reference, merge_window_bp=0) as instance:
        # The record's own spelling and the left-aligned one are the same event.
        record_side = instance.canonicalise("chr21", 302, "C", "CC")
        assert (record_side.position, record_side.ref, record_side.alt) == (300, "A", "AC")

        result = instance.frequencies([left_aligned_query])
        assert left_aligned_query in result, (
            "the record is spelled two bases right of the query key; a window that "
            "stops at the query position loses it to the fetch, not to the key"
        )
        assert by_population(result[left_aligned_query])["afr"].allele_count == 10

    # Without the reference neither side moves, so the miss is honest and omitted.
    with open_adapter(path, merge_window_bp=0) as plain:
        assert plain.frequencies([left_aligned_query]) == {}


def test_left_alignment_against_a_reference_restores_the_join(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """With a reference, the normaliser produces the key gnomAD actually holds.

    The reference window is the real chr21 sequence taken from the slice's own REF
    columns: 5033362-5033366 are ``A``, ``C``, ``A``, ``C``, ``C``. Given that,
    ``normalise_variants`` rolls both right-shifted spellings back onto
    ``5033364 A>AC`` and the lookup succeeds — which is the whole argument for
    wiring a reference FASTA into ingestion before trusting an indel's rarity.
    """
    reference = _ConstantReference("chr21", 5_033_362, "ACACC")
    records = [
        make_record(REPEAT_INSERTION_SHIFTED_ONE),
        make_record(REPEAT_INSERTION_SHIFTED_TWO),
    ]
    normalised = normalise_variants(records, reference=reference)
    ids = {record.coordinate.variant_id for record in normalised.variants}
    assert ids == {REPEAT_INSERTION}
    assert REPEAT_INSERTION in adapter.frequencies(sorted(ids))


def test_a_reference_only_ever_adds_joins_and_never_moves_one(
    adapter: GnomadSitesFrequencyAdapter,
    aligning_adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """Configuring a reference must be safe, not merely useful.

    Left-alignment now runs on the *record* side too, and a record whose key moved
    when it should not have would silently detach an answer that used to join —
    trading one silent miss for another. Asserted over the whole slice: the
    aligning adapter answers every variant the plain one does, with byte-identical
    frequencies, and the aligning adapter's extra answers are additional spellings
    rather than different numbers.
    """
    ids = [row_variant_id(row) for row in fixture_rows()]
    plain = adapter.frequencies(ids)
    aligned = aligning_adapter.frequencies(ids)
    assert set(plain) <= set(aligned)
    assert dump(plain) == dump({key: aligned[key] for key in plain})


def test_absence_is_still_absence_once_a_reference_is_configured(
    aligning_adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """GP-14 survives the change that made more things join.

    Two ways this could have been blurred, both checked. A widened fetch window
    reads neighbouring records, and a variant gnomAD has never seen must not
    collect one of them; and nothing anywhere may turn a missing record into
    ``allele_frequency=0.0``, which would convert "no evidence" into the strongest
    rarity evidence the dataset can give.
    """
    absent = [NOT_IN_GNOMAD, ANOTHER_ABSENT, UNLISTED_ALT_AT_A_KNOWN_SITE]
    assert aligning_adapter.frequencies(absent) == {}

    mixed = aligning_adapter.frequencies([*absent, TRUE_ZERO])
    assert set(mixed) == {TRUE_ZERO}
    assert by_population(mixed[TRUE_ZERO])[GLOBAL_POPULATION].allele_frequency == 0.0
    assert by_population(mixed[TRUE_ZERO])[GLOBAL_POPULATION].allele_number == 403_848

    # The three ``AN=0`` records stay omitted: observed, but with no frequency to
    # report, which is not the same fact as a genuine zero.
    assert NO_CALLED_ALLELES not in aligning_adapter.frequencies([NO_CALLED_ALLELES])


def test_the_left_aligning_lookup_is_byte_identical_across_runs_and_windows(
    aligning_adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """GP-30, over the path that now reads a reference.

    The merge window is still a query-count knob and not a correctness one, and the
    reference reads it triggers must not make the result depend on how the batch
    happened to be coalesced.
    """
    ids = [REPEAT_INSERTION, REPEAT_INSERTION_SHIFTED_TWO, DELETION, TRUE_ZERO, NOT_IN_GNOMAD]
    reference = _SliceReference.from_fixture()
    digests = {dump(aligning_adapter.frequencies(ids)) for _ in range(2)}
    for window in (0, 1, 250, 10_000):
        with open_adapter(FIXTURE, reference=reference, merge_window_bp=window) as instance:
            digests.add(dump(instance.frequencies(ids)))
    assert len(digests) == 1


def test_a_reference_that_raises_cannot_leak_a_coordinate_or_break_the_lookup(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """PRIV-09 and GP-14 over the surface left-alignment newly added.

    A ``ReferenceLookup`` is caller-supplied and may raise anything — pysam puts
    the region it failed on into its own messages, exactly as cyvcf2 does. The
    shift must therefore treat an unreadable base as absence of information and
    stop where it is, never propagate the exception, and never let its text out.
    The un-shifted representation is still a valid one; a crash, or a coordinate in
    a traceback, is not.
    """

    class _Leaky:
        def fetch(self, contig: str, start: int, end: int) -> str:
            raise RuntimeError(f"{contig}:{start}-{end}")

    with open_adapter(FIXTURE, reference=_Leaky()) as instance:
        try:
            result = instance.frequencies([REPEAT_INSERTION, REPEAT_INSERTION_SHIFTED_ONE])
        except Exception as exc:  # pragma: no cover - the assertion is the message
            rendered = "".join(traceback.format_exception(exc))
            pytest.fail(f"an unreadable reference must not raise:\n{rendered}")

    # Degraded exactly as far as the reference failed: the left-most spelling still
    # joins on the unconditional trim, the shifted one is omitted rather than
    # zero-filled, and neither answer was invented.
    assert set(result) == {REPEAT_INSERTION}
    assert by_population(result[REPEAT_INSERTION])[GLOBAL_POPULATION].allele_frequency == (
        pytest.approx(0.337086, rel=1e-6)
    )


def test_several_indels_at_one_position_stay_apart(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """chr21:5033364 carries nine records; their frequencies span five orders of magnitude.

    ``A>AC`` is 0.337 and ``A>ACC`` is 2.4e-05. Any confusion between two alleles
    anchored at one base is therefore not a rounding error, it is the difference
    between "one in three people" and "one in forty thousand".
    """
    ids = [
        "GRCh38:chr21:5033364:A:AC",
        "GRCh38:chr21:5033364:A:ACC",
        "GRCh38:chr21:5033364:A:AT",
        "GRCh38:chr21:5033364:AC:A",
        "GRCh38:chr21:5033364:A:G",
    ]
    result = adapter.frequencies(ids)
    assert set(result) == set(ids)
    observed = {key: value[0].allele_frequency for key, value in result.items()}
    assert observed["GRCh38:chr21:5033364:A:AC"] == pytest.approx(0.337086, rel=1e-6)
    assert observed["GRCh38:chr21:5033364:A:ACC"] == pytest.approx(2.43633e-05, rel=1e-6)
    assert observed["GRCh38:chr21:5033364:A:AT"] == pytest.approx(7.21842e-06, rel=1e-6)
    assert observed["GRCh38:chr21:5033364:AC:A"] == pytest.approx(2.25594e-05, rel=1e-6)
    assert observed["GRCh38:chr21:5033364:A:G"] == 0.0


def test_two_spellings_of_one_event_both_get_an_answer(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """Distinct caller keys that trim to the same record must not evict each other.

    ``chr21:5031991 GAA>GA`` and ``GA>G`` are one deletion written two ways. A
    key-to-single-ID match table lets the second silently overwrite the first, and
    whichever caller lost the race is told gnomAD has no record — the "rare
    variant" failure again, this time triggered by nothing more than batch order.
    """
    padded = "GRCh38:chr21:5031991:GAA:GA"
    forward = adapter.frequencies([padded, DELETION_AT_5031991])
    backward = adapter.frequencies([DELETION_AT_5031991, padded])
    assert set(forward) == set(backward) == {padded, DELETION_AT_5031991}
    assert forward[padded] == forward[DELETION_AT_5031991]


def test_a_long_deletion_joins(adapter: GnomadSitesFrequencyAdapter) -> None:
    """A 19-base REF span must still be reached by the coalesced region query.

    Also a second genuine zero: AC=0 over AN=1,195,784, FILTER ``AC0;AS_VQSR``.
    """
    populations = by_population(adapter.frequencies([LONG_DELETION])[LONG_DELETION])
    assert populations[GLOBAL_POPULATION].allele_frequency == 0.0
    assert populations[GLOBAL_POPULATION].allele_number == 1_195_784
    assert populations["nfe"].allele_number == 942_034
    assert populations[GLOBAL_POPULATION].filter_status == "AC0;AS_VQSR"


def test_a_common_deletion_keeps_its_per_population_split(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """Indel and ancestry traps at once: ``sas`` 0.171 against ``eas`` 0.0048."""
    populations = by_population(adapter.frequencies([COMMON_DELETION])[COMMON_DELETION])
    assert populations["sas"].allele_frequency == pytest.approx(0.170910, rel=1e-6)
    assert populations["eas"].allele_frequency == pytest.approx(0.00480769, rel=1e-6)
    assert populations["sas"].allele_number == 41_086


# --------------------------------------------------------------------- contig naming


def test_contig_map_is_resolved_from_the_tabix_index(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """Asserted directly, not inferred from a lookup that happened to succeed.

    gnomAD v4 uses UCSC-style contigs, so the mapping is the identity here — but
    it is *resolved* rather than assumed, because getting it wrong fails by
    finding nothing, which reads exactly like "gnomAD has no record" for an entire
    chromosome.
    """
    assert adapter.contig_map == {"chr21": "chr21"}
    assert adapter.available_contigs == ("chr21",)


def test_the_header_declares_every_contig_but_the_index_holds_one() -> None:
    """Why the contig list comes from the index and not from the header.

    Each gnomAD shard's header declares all 24 contigs. An adapter that trusted it
    would have every shard claim the whole genome, and whichever shard was opened
    first would answer for chromosomes it does not contain — silently, with
    nothing.
    """
    header = fixture_header_lines()
    declared = [line for line in header if line.startswith("##contig=")]
    assert len(declared) == 24
    tabix = pysam.TabixFile(str(FIXTURE))
    try:
        assert list(tabix.contigs) == ["chr21"]
    finally:
        tabix.close()


def test_a_bare_ensembl_contig_in_a_variant_id_still_joins(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """``21`` is normalised to ``chr21`` before the query, and the caller's key is returned."""
    bare = "GRCh38:21:5035658:C:T"
    result = adapter.frequencies([bare])
    assert set(result) == {bare}
    assert result[bare] == adapter.frequencies([COMMON])[COMMON]


def test_available_contigs_are_in_karyotype_order(tmp_path: Path) -> None:
    """GP-30: chr2 before chr10, never lexicographic and never dict order."""
    header = mini_header(contigs=("chr2", "chr10", "chrX"))
    body = (
        "".join(
            f"{contig}\t100\t.\tA\tG\t.\tPASS\tAC=1;AN=100;AF=0.01\n"
            for contig in ("chr2", "chr10")
        )
        + "chrX\t100\t.\tA\tG\t.\tPASS\tAC=1;AN=100;AF=0.01\n"
    )
    path = write_sites_vcf(tmp_path, body, header=header)
    with open_adapter(path) as instance:
        assert instance.available_contigs == ("chr2", "chr10", "chrX")


def test_cross_build_variant_ids_are_refused_not_silently_missed(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """A GRCh37 key against a GRCh38 release must raise, not return "no record"."""
    with pytest.raises(GenomeBuildMismatchError, match="GRCh37"):
        adapter.frequencies(["GRCh37:chr21:5035658:C:T"])


def test_malformed_variant_ids_are_refused_without_echoing_them(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    with pytest.raises(ValueError, match="Malformed variant ID") as excinfo:
        adapter.frequencies(["chr21:5035658:C:T"])
    assert "5035658" not in str(excinfo.value)

    with pytest.raises(ValueError, match="non-integer position") as excinfo:
        adapter.frequencies(["GRCh38:chr21:not-a-position:C:T"])
    assert "not-a-position" not in str(excinfo.value)


def test_a_non_canonical_contig_is_refused(adapter: GnomadSitesFrequencyAdapter) -> None:
    """Alt and decoy contigs are out of scope, and saying so beats returning nothing."""
    with pytest.raises(ValueError, match="Non-canonical contig"):
        adapter.frequencies(["GRCh38:chr21_KI270741v1_random:100:A:G"])


# ------------------------------------------------------------------------ determinism


def test_repeat_calls_are_byte_identical(adapter: GnomadSitesFrequencyAdapter) -> None:
    """GP-30. Same inputs, same bytes — including across a freshly opened adapter."""
    ids = [COMMON, TRUE_ZERO, DIVERGENT, NOT_IN_GNOMAD, DELETION]
    first = adapter.frequencies(ids)
    second = adapter.frequencies(ids)
    with open_adapter(FIXTURE) as other:
        third = other.frequencies(ids)
    assert dump(first) == dump(second) == dump(third)


def test_input_order_does_not_change_per_variant_results(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    ids = [COMMON, TRUE_ZERO, DIVERGENT, DELETION, INBREEDING_COEFF]
    assert adapter.frequencies(ids) == adapter.frequencies(list(reversed(ids)))


def test_result_keys_follow_caller_order(adapter: GnomadSitesFrequencyAdapter) -> None:
    """Deterministic iteration order, taken from the caller rather than from a set."""
    ids = [DIVERGENT, NOT_IN_GNOMAD, COMMON, ANOTHER_ABSENT, TRUE_ZERO]
    assert list(adapter.frequencies(ids)) == [DIVERGENT, COMMON, TRUE_ZERO]


def test_duplicate_ids_are_collapsed(adapter: GnomadSitesFrequencyAdapter) -> None:
    assert adapter.frequencies([COMMON]) == adapter.frequencies([COMMON] * 4)


def test_the_merge_window_never_changes_the_answer() -> None:
    """The coalescing window is a query-count knob, not a correctness one.

    Every fetched record is still matched on the exact ``(pos, ref, alt)`` key, so
    a wider window can only cause records to be read and discarded. Proved by
    hashing the whole-slice result at four window sizes.
    """
    ids = [row_variant_id(row) for row in fixture_rows()]
    digests: set[str] = set()
    for window in (0, 1, 250, 10_000):
        with open_adapter(FIXTURE, merge_window_bp=window) as instance:
            digests.add(dump(instance.frequencies(ids)))
    assert len(digests) == 1


def test_merge_query_regions_coalesces_only_within_the_window() -> None:
    assert merge_query_regions([], 100) == ()
    assert merge_query_regions([(5, 5), (5, 5)], 100) == ((5, 5),)
    assert merge_query_regions([(300, 300), (100, 100), (150, 150)], 100) == (
        (100, 150),
        (300, 300),
    )
    assert merge_query_regions([(300, 300), (100, 100), (150, 150)], 200) == ((100, 300),)
    assert merge_query_regions([(100, 100), (150, 150), (400, 400)], 0) == (
        (100, 100),
        (150, 150),
        (400, 400),
    )
    # A long REF extends the span, so the region query still reaches its far end.
    assert merge_query_regions([(100, 500), (520, 520)], 25) == ((100, 520),)


# --------------------------------------------------------------- release identity


def test_a_release_the_filename_contradicts_is_refused() -> None:
    """The header has no release string, so the filename is the only cross-check."""
    with pytest.raises(AdapterUnavailableError, match="whole version token"):
        open_adapter(FIXTURE, release="v2.1.1")


def test_a_release_prefix_is_not_accepted_as_the_release() -> None:
    """``"v4" in "gnomad.exomes.v4.1..."`` is true, and would mislabel every citation.

    Matched as a whole ``v``-dotted token instead, so neither a prefix nor an
    extension of the real release passes.
    """
    with pytest.raises(AdapterUnavailableError, match=r"\['v4.1'\]"):
        open_adapter(FIXTURE, release="v4")
    with pytest.raises(AdapterUnavailableError, match="whole version token"):
        open_adapter(FIXTURE, release="v4.10")


def test_a_wrong_subset_is_refused() -> None:
    """Exome and genome callsets have different cohort sizes; ADR 0010 guards on those."""
    with pytest.raises(AdapterUnavailableError, match="declared subset"):
        open_adapter(FIXTURE, subset="genomes")


def test_an_empty_release_or_subset_is_refused() -> None:
    with pytest.raises(AdapterUnavailableError, match="non-empty"):
        open_adapter(FIXTURE, release="   ")


def test_shards_from_different_pipeline_runs_are_refused(tmp_path: Path) -> None:
    """The only signal that a directory mixes releases, since the header has no version."""
    body = "chr21\t100\t.\tA\tG\t.\tPASS\tAC=1;AN=100;AF=0.01\n"
    write_sites_vcf(tmp_path, body, name="gnomad.exomes.v4.1.sites.chr21.mini.vcf")
    write_sites_vcf(
        tmp_path,
        "chr22\t100\t.\tA\tG\t.\tPASS\tAC=1;AN=100;AF=0.01\n",
        header=mini_header(contigs=("chr22",), hail_version="0.2.99-deadbeef"),
        name="gnomad.exomes.v4.1.sites.chr22.mini.vcf",
    )
    with pytest.raises(AdapterUnavailableError, match="different pipeline run"):
        open_adapter(tmp_path)


def test_two_shards_claiming_one_contig_are_refused(tmp_path: Path) -> None:
    """Order-dependent lookup is not resolved silently (GP-30)."""
    body = "chr21\t100\t.\tA\tG\t.\tPASS\tAC=1;AN=100;AF=0.01\n"
    write_sites_vcf(tmp_path, body, name="gnomad.exomes.v4.1.sites.chr21.a.vcf")
    write_sites_vcf(tmp_path, body, name="gnomad.exomes.v4.1.sites.chr21.b.vcf")
    with pytest.raises(AdapterUnavailableError, match="indexed by both"):
        open_adapter(tmp_path)


def test_a_header_without_a_build_is_refused(tmp_path: Path) -> None:
    """GP-11: guessing GRCh38 would mis-locate every variant by megabases."""
    path = write_sites_vcf(
        tmp_path,
        "chr21\t100\t.\tA\tG\t.\tPASS\tAC=1;AN=100;AF=0.01\n",
        header=mini_header(assembly="gnomAD_hg00"),
    )
    with pytest.raises(AdapterUnavailableError, match="unrecognised genome build"):
        open_adapter(path)


def test_a_source_with_no_per_population_fields_is_refused(tmp_path: Path) -> None:
    """A global-only source makes the ADR 0010 guard silently inert."""
    path = write_sites_vcf(
        tmp_path,
        "chr21\t100\t.\tA\tG\t.\tPASS\tAC=1;AN=100;AF=0.01\n",
        header=mini_header(groups=()),
    )
    with pytest.raises(AdapterUnavailableError, match="AF_<group>"):
        open_adapter(path)


def test_a_grch37_release_joins_grch37_ids(tmp_path: Path) -> None:
    """The build in the join key comes from the release header, never a hardcoded GRCh38.

    A hardcoded build makes the query key and the record key disagree, so every
    lookup returns nothing and the run reports "gnomAD has no record" for the
    whole genome while raising no error at all.
    """
    path = write_sites_vcf(
        tmp_path,
        "chr21\t100\t.\tA\tG\t.\tPASS\tAC=1;AN=100;AF=0.01\n",
        header=mini_header(assembly="gnomAD_GRCh37"),
    )
    with open_adapter(path) as instance:
        assert instance.build is GenomeBuild.GRCH37
        assert set(instance.frequencies(["GRCh37:chr21:100:A:G"])) == {"GRCh37:chr21:100:A:G"}
        with pytest.raises(GenomeBuildMismatchError):
            instance.frequencies(["GRCh38:chr21:100:A:G"])


# ------------------------------------------------------------------ truncation safety


def test_bgzf_eof_detects_a_truncated_stream(tmp_path: Path) -> None:
    """The check that matters: a half-downloaded shard decompresses cleanly.

    A truncated ``.vcf.bgz`` reads perfectly up to its cut point, so every variant
    past it reports as absent — which downstream reads as "novel, therefore
    interesting". The empty-BGZF end-of-file block is what proves the stream is
    whole.
    """
    assert has_bgzf_eof(FIXTURE)
    truncated = tmp_path / "gnomad.exomes.v4.1.sites.chr21.part.vcf.bgz"
    truncated.write_bytes(FIXTURE.read_bytes()[: -len(BGZF_EOF)])
    assert not has_bgzf_eof(truncated)
    assert not has_bgzf_eof(tmp_path / "absent.vcf.bgz")
    tiny = tmp_path / "tiny.vcf.bgz"
    tiny.write_bytes(b"\x1f\x8b")
    assert not has_bgzf_eof(tiny)


def test_check_source_complete_names_each_failure(tmp_path: Path) -> None:
    """Structural facts and file names only — never record content (PRIV-09)."""
    missing = check_source_complete(tmp_path / "absent.vcf.bgz")
    assert not missing.is_complete
    assert missing.reasons == ("file does not exist",)

    orphan = tmp_path / "gnomad.exomes.v4.1.sites.chr21.orphan.vcf.bgz"
    orphan.write_bytes(FIXTURE.read_bytes())
    status = check_source_complete(orphan)
    assert not status.is_complete
    assert status.has_bgzf_eof
    assert status.reasons == ("no .tbi/.csi index beside it",)
    assert status.size_bytes == FIXTURE.stat().st_size
    assert index_path_for(orphan) is None
    assert index_path_for(FIXTURE) == FIXTURE.with_name(FIXTURE.name + ".tbi")

    complete = check_source_complete(FIXTURE)
    assert complete.is_complete
    assert complete.reasons == ()
    assert complete.size_stable is None, "not probed is unknown, not stable"


def test_the_stability_probe_rejects_a_file_that_is_still_growing(tmp_path: Path) -> None:
    """Catches the pathological writer that has flushed a valid-looking tail.

    The sleeper is injected so the probe costs no wall time and reads no clock
    (GP-30); it is also what appends the bytes, standing in for the download.
    """
    growing = tmp_path / "gnomad.exomes.v4.1.sites.chr21.growing.vcf.bgz"
    growing.write_bytes(FIXTURE.read_bytes())
    calls: list[float] = []

    def grow(seconds: float) -> None:
        calls.append(seconds)
        with growing.open("ab") as handle:
            handle.write(b"more")

    status = check_source_complete(growing, stability_probe_seconds=0.5, sleeper=grow)
    assert calls == [0.5]
    assert status.size_stable is False
    assert not status.is_complete
    assert "still being written" in status.reasons[-1]

    stable = check_source_complete(FIXTURE, stability_probe_seconds=0.5, sleeper=lambda _s: None)
    assert stable.size_stable is True
    assert stable.is_complete


def test_the_probe_is_off_by_default_so_construction_never_sleeps() -> None:
    def explode(_seconds: float) -> None:
        pytest.fail("the stability probe must be opt-in")

    with open_adapter(FIXTURE, sleeper=explode) as instance:
        assert instance.available_contigs == ("chr21",)


def test_a_directory_of_incomplete_shards_refuses_to_construct(tmp_path: Path) -> None:
    """Refusing beats annotating against a dataset with a hole in it."""
    truncated = tmp_path / "gnomad.exomes.v4.1.sites.chr21.vcf.bgz"
    truncated.write_bytes(FIXTURE.read_bytes()[: -len(BGZF_EOF)])
    with pytest.raises(AdapterUnavailableError, match="partially written BGZF stream"):
        open_adapter(tmp_path)


def test_incomplete_shards_are_named_and_never_read(tmp_path: Path) -> None:
    """A shard still arriving is excluded *and* reported, not silently skipped."""
    good = tmp_path / FIXTURE.name
    good.write_bytes(FIXTURE.read_bytes())
    index = tmp_path / (FIXTURE.name + ".tbi")
    index.write_bytes((FIXTURE.parent / (FIXTURE.name + ".tbi")).read_bytes())
    arriving = tmp_path / "gnomad.exomes.v4.1.sites.chr22.vcf.bgz"
    arriving.write_bytes(FIXTURE.read_bytes()[:10_000])

    with open_adapter(tmp_path) as instance:
        assert instance.available_contigs == ("chr21",)
        assert len(instance.incomplete_sources) == 1
        reason = instance.incomplete_sources[0]
        assert reason.startswith("gnomad.exomes.v4.1.sites.chr22.vcf.bgz:")
        assert "no .tbi/.csi index beside it" in reason
        assert "BGZF end-of-file marker missing" in reason
        # And querying it is refused, not answered with silence.
        with pytest.raises(AdapterUnavailableError, match=r"contig\(s\) \['chr22'\]"):
            instance.frequencies(["GRCh38:chr22:100:A:G"])


def test_require_contigs_turns_a_still_arriving_shard_into_a_hard_failure(
    tmp_path: Path,
) -> None:
    """ "No frequency data for all of chr22" must be distinguishable from "chr22 is rare".

    Without this the run would silently score every chr22 candidate as having no
    frequency evidence, which is the configured mid-range score rather than a
    missing-resource error.
    """
    good = tmp_path / FIXTURE.name
    good.write_bytes(FIXTURE.read_bytes())
    (tmp_path / (FIXTURE.name + ".tbi")).write_bytes(
        (FIXTURE.parent / (FIXTURE.name + ".tbi")).read_bytes()
    )
    arriving = tmp_path / "gnomad.exomes.v4.1.sites.chr22.vcf.bgz"
    arriving.write_bytes(FIXTURE.read_bytes()[:10_000])

    with pytest.raises(AdapterUnavailableError, match=r"required contig\(s\) \['chr22'\]"):
        open_adapter(tmp_path, require_contigs=["chr22"])
    # chr21 is present, so requiring only it succeeds.
    with open_adapter(tmp_path, require_contigs=["chr21"]) as instance:
        assert instance.available_contigs == ("chr21",)


def test_query_on_incomplete_shard_fails_closed(tmp_path: Path) -> None:
    """One complete shard, one truncated: querying the truncated contig must raise.

    The regression this exists for is the one that matters most in this module.
    Construction succeeds (the complete shard is usable), so nothing warns; the
    lookup then returns a mapping with the chr22 variant simply missing, which is
    byte-for-byte what "gnomAD has never seen this variant" looks like. A common
    chr22 allele would sail past the common-frequency down-rank and be ranked as
    an ultra-rare candidate — on the strength of a download that had not finished.
    """
    good = tmp_path / FIXTURE.name
    good.write_bytes(FIXTURE.read_bytes())
    (tmp_path / (FIXTURE.name + ".tbi")).write_bytes(
        (FIXTURE.parent / (FIXTURE.name + ".tbi")).read_bytes()
    )
    truncated = tmp_path / "gnomad.exomes.v4.1.sites.chr22.vcf.bgz"
    truncated.write_bytes(FIXTURE.read_bytes()[: -len(BGZF_EOF)])

    with open_adapter(tmp_path) as instance:
        assert instance.available_contigs == ("chr21",), "construction still succeeds"
        with pytest.raises(AdapterUnavailableError) as excinfo:
            instance.frequencies(["GRCh38:chr22:42126611:C:T"])
        message = str(excinfo.value)
        assert "chr22" in message
        assert "42126611" not in message, "PRIV-09: the position must not be echoed"
        assert "gnomad.exomes.v4.1.sites.chr22.vcf.bgz" in message, "name the missing shard"

        # The deliberate way through, with the gap in the result rather than hidden.
        partial = instance.lookup_partial(["GRCh38:chr22:42126611:C:T", "GRCh38:chr21:5035658:C:T"])
        assert isinstance(partial, FrequencyLookup)
        assert partial.unqueryable_contigs == ("chr22",)
        assert set(partial.frequencies) == {COMMON}
        assert any("chr22" in item for item in partial.incomplete_sources)


def test_a_backend_failure_never_echoes_the_queried_region(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """PRIV-09 on the error path. cyvcf2 puts the region into what it raises.

    Verified against the real backend behaviour: an unwrapped failure produced a
    bare ``RuntimeError: chr21:5031905-5031905``, i.e. a proband coordinate in the
    terminal, the log, the crash report and an agent's context. Both the region
    call and the iteration are guarded — htslib defers work to the first
    ``next()``, so guarding only the construction would move the leak rather than
    remove it — and the replacement is raised with the context suppressed, because
    chaining would print the original message again.
    """
    import mva.annotation.gnomad_sites as module

    real_open = real_reader_factory()

    class _Exploding:
        def __init__(self, inner: _ReaderLike, *, on_iteration: bool) -> None:
            self._inner = inner
            self._on_iteration = on_iteration

        @property
        def raw_header(self) -> str:
            return self._inner.raw_header

        def __call__(self, region: str) -> Iterator[object]:
            if not self._on_iteration:
                raise RuntimeError(region)

            def _iterate() -> Iterator[object]:
                raise OSError(f"htslib failed reading {region}")
                yield  # pragma: no cover - unreachable, makes this a generator

            return _iterate()

        def close(self) -> None:
            self._inner.close()

    for on_iteration in (False, True):

        def exploding_open(path: Path, *, flag: bool = on_iteration) -> _ReaderLike:
            return _Exploding(real_open(path), on_iteration=flag)

        monkeypatch.setattr(module, "_open_reader", exploding_open)
        with open_adapter(FIXTURE) as instance, pytest.raises(AdapterUnavailableError) as excinfo:
            instance.frequencies([TRUE_ZERO])

        rendered = "".join(
            traceback.format_exception(
                type(excinfo.value), excinfo.value, excinfo.value.__traceback__
            )
        )
        assert "5031905" not in rendered, "the position reached the traceback"
        assert "chr21:" not in rendered, "the region reached the traceback"
        assert "htslib failed reading" not in rendered, "the backend's own text was chained"
        assert "AdapterUnavailableError" in rendered
        # The backend's exception *class* is kept — it is a diagnostic and carries
        # no patient data — while its message and frame are not.
        assert ("RuntimeError" in rendered) is not on_iteration
        assert ("OSError" in rendered) is on_iteration
        assert re.search(r"<region:[0-9a-f]{8}>", rendered), "a correlation handle is still given"


def test_a_shard_with_an_unreadable_index_is_excluded_not_fatal(tmp_path: Path) -> None:
    """The ``.tbi`` may itself be mid-download; that is a missing shard, not a crash."""
    good = tmp_path / FIXTURE.name
    good.write_bytes(FIXTURE.read_bytes())
    (tmp_path / (FIXTURE.name + ".tbi")).write_bytes(
        (FIXTURE.parent / (FIXTURE.name + ".tbi")).read_bytes()
    )
    broken = tmp_path / "gnomad.exomes.v4.1.sites.chr22.vcf.bgz"
    broken.write_bytes(FIXTURE.read_bytes())
    (tmp_path / "gnomad.exomes.v4.1.sites.chr22.vcf.bgz.tbi").write_bytes(b"not-an-index" * 4)

    with open_adapter(tmp_path) as instance:
        assert instance.available_contigs == ("chr21",)
        assert any("index unreadable" in item for item in instance.incomplete_sources)


def test_a_missing_path_is_a_clear_failure(tmp_path: Path) -> None:
    with pytest.raises(AdapterUnavailableError, match="neither a file nor a directory"):
        open_adapter(tmp_path / "nowhere")
    with pytest.raises(AdapterUnavailableError, match=r"no \*\.vcf\.bgz files found"):
        open_adapter(tmp_path)


# ------------------------------------------------------------------- multi-allelic


def test_multi_alt_records_are_split_per_allele(tmp_path: Path) -> None:
    """VCF permits several ALTs per line; gnomAD's v4.1 exomes release uses none.

    A full scan of chr21 found 0 in 2,188,842 records, so this shape cannot be cut
    from the real file — it is built here instead. A reader that assumed one ALT
    per line would attribute allele 1's ``Number=A`` frequency to allele 2, which
    is a *wrong* number rather than a missing one.
    """
    path = write_sites_vcf(
        tmp_path,
        "chr21\t100\t.\tA\tG,T\t.\tPASS\t"
        "AC=1,2;AN=100;AF=0.01,0.02;nhomalt=0,1;"
        "AC_nfe=1,0;AN_nfe=50;AF_nfe=0.02,0.00;nhomalt_nfe=0,0\n",
    )
    with open_adapter(path) as instance:
        result = instance.frequencies(["GRCh38:chr21:100:A:G", "GRCh38:chr21:100:A:T"])
    first = by_population(result["GRCh38:chr21:100:A:G"])
    second = by_population(result["GRCh38:chr21:100:A:T"])
    assert first[GLOBAL_POPULATION].allele_frequency == pytest.approx(0.01, rel=1e-6)
    assert second[GLOBAL_POPULATION].allele_frequency == pytest.approx(0.02, rel=1e-6)
    assert first[GLOBAL_POPULATION].allele_count == 1
    assert second[GLOBAL_POPULATION].allele_count == 2
    assert second[GLOBAL_POPULATION].homozygote_count == 1
    # AN is Number=1: the cohort size is shared by both alleles, and is not indexed.
    assert first[GLOBAL_POPULATION].allele_number == second[GLOBAL_POPULATION].allele_number == 100
    assert first["nfe"].allele_frequency == pytest.approx(0.02, rel=1e-6)
    assert second["nfe"].allele_frequency == 0.0


def test_an_unlisted_alt_at_a_multi_allelic_site_is_still_absent(tmp_path: Path) -> None:
    path = write_sites_vcf(
        tmp_path,
        "chr21\t100\t.\tA\tG,T\t.\tPASS\tAC=1,2;AN=100;AF=0.01,0.02;AF_nfe=0.02,0.0;AN_nfe=50\n",
    )
    with open_adapter(path) as instance:
        assert instance.frequencies(["GRCh38:chr21:100:A:C"]) == {}


# -------------------------------------------------------------------- performance


def test_a_batch_is_one_region_query_per_coalesced_span(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A 5,000-variant batch must not be 5,000 index seeks.

    Measured by counting region queries rather than by timing, which would be
    flaky. The un-coalesced implementation issued one query per distinct position;
    on the real 2.2 GB chr21 shard that is ~56 s for 5,000 clustered variants,
    because each seek decompresses a fresh 64 KB BGZF block for an 8 KB record.
    """
    import mva.annotation.gnomad_sites as module

    real_open = real_reader_factory()
    queries: list[str] = []

    class _Counting:
        def __init__(self, inner: _ReaderLike) -> None:
            self._inner = inner

        @property
        def raw_header(self) -> str:
            return self._inner.raw_header

        def __call__(self, region: str) -> Iterator[object]:
            queries.append(region)
            return self._inner(region)

        def close(self) -> None:
            self._inner.close()

    def counting_open(path: Path) -> _ReaderLike:
        return _Counting(real_open(path))

    monkeypatch.setattr(module, "_open_reader", counting_open)

    ids = [row_variant_id(row) for row in fixture_rows()]
    assert len(ids) == 154

    with open_adapter(FIXTURE, merge_window_bp=0) as instance:
        uncoalesced = instance.frequencies(ids)
    assert len(queries) == 55, "one query per distinct REF span, the old behaviour"

    queries.clear()
    with open_adapter(FIXTURE) as instance:
        result = instance.frequencies(ids)

    # The slice's 154 records sit in 6 clusters spread over 4.8 Mb, so a 1 kb
    # merge window collapses 55 seeks into 6 sequential reads.
    assert len(queries) == 6
    assert all(region.startswith("chr21:") for region in queries)
    assert result == uncoalesced
    # The three ``AN=0`` records in the slice have no allele frequency to report,
    # so they are omitted rather than zero-filled (GP-14).
    assert len(result) == 151


def test_one_reader_is_opened_per_shard(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Handles are opened once at construction and reused, not per lookup."""
    import mva.annotation.gnomad_sites as module

    real_open = real_reader_factory()
    opened: list[Path] = []

    def counting_open(path: Path) -> _ReaderLike:
        opened.append(path)
        return real_open(path)

    monkeypatch.setattr(module, "_open_reader", counting_open)
    with open_adapter(FIXTURE) as instance:
        for _ in range(5):
            instance.frequencies([COMMON, TRUE_ZERO, DIVERGENT])
    assert len(opened) == 1


@requires_full_release
def test_the_full_release_directory_excludes_shards_still_downloading() -> None:
    """The reason :func:`check_source_complete` exists, against the live directory.

    This adapter is pointed at a directory whose 250 GB of shards arrive over
    hours. Every shard without a BGZF end-of-file marker must be named and
    excluded rather than read, because a truncated shard answers "absent" for
    every variant past its cut point.
    """
    with GnomadSitesFrequencyAdapter(
        FULL_RELEASE_DIR, release=RELEASE, subset=SUBSET, require_contigs=["chr21"]
    ) as instance:
        assert "chr21" in instance.available_contigs
        for reason in instance.incomplete_sources:
            assert reason.split(":")[0].endswith((".vcf.bgz", ".vcf.gz"))
            assert "chr21" not in reason.split(":")[0]
        assert instance.header_fingerprint.startswith("gnomAD/GRCh38/")


@requires_full_release
def test_five_thousand_variants_against_the_real_shard() -> None:
    """The performance claim, on the real 2,250,319,809-byte file.

    Asserted as a generous ceiling rather than a benchmark: the point is that the
    batch is a handful of coalesced region queries, not 5,000 seeks. Measured at
    ~0.5 s here against ~56 s before the query coalescing was added.
    """
    tabix = pysam.TabixFile(str(FULL_CHR21))
    try:
        ids: list[str] = []
        for line in tabix.fetch("chr21", 25_880_000, 26_400_000):
            columns = line.split("\t", 6)
            if "," in columns[4]:
                continue
            ids.append(f"GRCh38:chr21:{columns[1]}:{columns[3]}:{columns[4]}")
            if len(ids) >= 5000:
                break
    finally:
        tabix.close()
    assert len(ids) == 5000

    with GnomadSitesFrequencyAdapter(FULL_CHR21, release=RELEASE, subset=SUBSET) as instance:
        started = time.perf_counter()
        result = instance.frequencies(ids)
        elapsed = time.perf_counter() - started
    assert len(result) > 4500
    assert elapsed < 20.0, f"5,000-variant lookup took {elapsed:.1f}s"


# -------------------------------------------------------------- service integration


def test_the_annotation_service_files_this_adapter_s_output_correctly(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """End to end through ``annotate_variants``, with nothing in the service changed.

    Three things are locked at once, all properties of what this adapter emits:

    * a variant gnomAD has no record for gets the "frequency data unavailable"
      evidence item, not a frequency item at 0.0 (GP-14);
    * a genuine ``AC=0`` gets a real frequency item whose numeric value is 0.0 —
      the same number, a completely different claim;
    * ``EvidenceItem`` refuses a DATABASE_ASSERTION whose citation has no version,
      so this proves the release string survives the whole way.
    """
    base = load_default_adapters(KNOWLEDGE_ROOT, KNOWLEDGE_ROOT / "manifests" / "knowledge.yaml")
    # ``clinical`` is supplied because ``annotate_variants`` currently zips two
    # descriptors against three indexes under ``strict=True`` when the clinical
    # slot is unbound; that is a bug in ``annotation.service``, not in this
    # adapter, and it is out of scope for this module to change.
    adapters = AdapterSet(consequence=base.consequence, frequency=adapter, clinical=base.clinical)
    records = [make_record(TRUE_ZERO), make_record(DIVERGENT), make_record(NOT_IN_GNOMAD)]
    result = annotate_variants(
        records, adapters=adapters, clock=FixedClock(datetime(2026, 1, 1, tzinfo=UTC))
    )

    annotated = {record.coordinate.variant_id: record for record in result.variants}
    assert annotated[TRUE_ZERO].has_frequency_data
    assert annotated[DIVERGENT].has_frequency_data
    assert not annotated[NOT_IN_GNOMAD].has_frequency_data

    frequency_items = [
        item
        for item in result.evidence
        if item.tool == ADAPTER_NAME and item.tier is AssertionTier.DATABASE_ASSERTION
    ]
    zero_items = [item for item in frequency_items if item.subject_id == TRUE_ZERO]
    assert zero_items
    assert all(item.numeric_value == 0.0 for item in zero_items)
    citation = zero_items[0].citation
    assert citation is not None
    assert citation.source == "gnomAD_exomes"
    assert citation.version == RELEASE
    assert not any(item.subject_id == NOT_IN_GNOMAD for item in frequency_items)


def test_a_real_adapter_carries_no_synthetic_disclosure(
    adapter: GnomadSitesFrequencyAdapter,
) -> None:
    """GP-20 cuts both ways: a real source must not be labelled a mock either."""
    base = load_default_adapters(KNOWLEDGE_ROOT, KNOWLEDGE_ROOT / "manifests" / "knowledge.yaml")
    adapters = AdapterSet(consequence=base.consequence, frequency=adapter, clinical=base.clinical)
    result = annotate_variants(
        [make_record(COMMON)],
        adapters=adapters,
        clock=FixedClock(datetime(2026, 1, 1, tzinfo=UTC)),
    )
    items = [item for item in result.evidence if item.tool == ADAPTER_NAME]
    assert items
    for item in items:
        assert SYNTHETIC_STANDIN_LIMITATION not in item.limitations
        assert item.limitations, "a real source still has to state its limitations (GP-17)"


# ------------------------------------------------------------------------ lifecycle


def test_close_is_idempotent_and_use_after_close_raises() -> None:
    """A lookup against released handles would return nothing — which reads as absence."""
    instance = open_adapter(FIXTURE)
    assert instance.frequencies([COMMON])
    instance.close()
    instance.close()
    with pytest.raises(AdapterUnavailableError, match="has been closed"):
        instance.frequencies([COMMON])


def test_the_context_manager_closes_on_the_way_out() -> None:
    with open_adapter(FIXTURE) as instance:
        assert instance.frequencies([COMMON])
    with pytest.raises(AdapterUnavailableError, match="has been closed"):
        instance.frequencies([COMMON])


def test_a_failed_construction_leaves_no_open_handles(tmp_path: Path) -> None:
    """A shard opened before the failing one must not leak its htslib handle."""
    good = tmp_path / FIXTURE.name
    good.write_bytes(FIXTURE.read_bytes())
    (tmp_path / (FIXTURE.name + ".tbi")).write_bytes(
        (FIXTURE.parent / (FIXTURE.name + ".tbi")).read_bytes()
    )
    write_sites_vcf(
        tmp_path,
        "chr22\t100\t.\tA\tG\t.\tPASS\tAC=1;AN=100;AF=0.01\n",
        header=mini_header(contigs=("chr22",), hail_version="0.2.99-deadbeef"),
        name="gnomad.exomes.v4.1.sites.chr22.mini.vcf",
    )
    closed: list[str] = []
    import mva.annotation.gnomad_sites as module

    real_open = real_reader_factory()

    class _Tracking:
        def __init__(self, inner: _ReaderLike, name: str) -> None:
            self._inner = inner
            self._name = name

        @property
        def raw_header(self) -> str:
            return self._inner.raw_header

        def __call__(self, region: str) -> Iterator[object]:
            return self._inner(region)

        def close(self) -> None:
            closed.append(self._name)
            self._inner.close()

    def tracking_open(path: Path) -> _ReaderLike:
        return _Tracking(real_open(path), path.name)

    with pytest.MonkeyPatch.context() as patcher:
        patcher.setattr(module, "_open_reader", tracking_open)
        with pytest.raises(AdapterUnavailableError, match="different pipeline run"):
            open_adapter(tmp_path)
    assert FIXTURE.name in closed
