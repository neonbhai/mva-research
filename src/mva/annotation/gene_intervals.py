"""Gene assignment by interval join against a real, local MANE GRCh38 release.

Without this, the pipeline scores zero. ``prioritization.pairing.generate_pairs``
is gene-scoped: it groups variants by ``VariantRecord.gene_symbols``, which is a
derived property reading ``csq.gene_symbol`` off ``VariantRecord.consequences``.
With no adapter populating ``consequences``, every record's ``gene_symbols`` is
``()``, every gene group is empty, and ``generate_pairs`` returns ``()`` for a
whole-genome VCF. A compound heterozygote is by definition two variants *in one
gene*, so a Track 1 answer cannot be formed at all until something assigns genes.

This module does that, and nothing more than that. It answers one question —
*which annotated genes does this coordinate fall inside?* — from the gene rows of
the MANE GTF, joined to the MANE summary for the stable identifiers. It is a
locus lookup, not a variant-effect predictor: it does not know whether a variant
is missense, nonsense or splice-disrupting, and it never guesses. The adapter is
named ``mane-interval-join`` precisely so that no report footer can mistake it
for VEP, SnpEff or Nirvana.

Why local files, permanently (PRIV-05)
--------------------------------------

A proband coordinate is not metadata. Sending ``chr15:40200239`` to Ensembl REST
to ask which gene it is in discloses patient genetic data to a third party,
irreversibly. ``mva.annotation`` may not import a network client, and this module
reads two files a separate, public-only acquisition step already downloaded and
hashed. A pinned local release also gives byte-identical repeat runs (GP-30),
which a live API cannot.

The five traps this module exists to not fall into
--------------------------------------------------

1. **A variant can be in more than one gene, and that must survive.** MANE v1.5
   holds 4,470 overlapping gene pairs on the primary contigs. CEP57 — an MVA
   gene — overlaps FAM76B on one side and MTMR2 on the other, so a variant at
   chr11:95,790,000 is genuinely in two genes on opposite strands. Every
   overlapping gene is returned, in a documented total order. Collapsing to one
   is the same family of data-loss bug as collapsing to the canonical transcript,
   which ``knowledge/adapters/README.md`` forbids outright.

2. **Absence is representable (GP-14).** A variant in an intergenic region has no
   gene. It is *omitted* from the returned mapping — never present with an empty
   tuple, and never given a nearest-gene assignment by default. Nearest-gene is
   available (:class:`GeneRelation.NEAREST`) but is off unless a caller asks for
   it by distance, and every such assignment is labelled as an inference with the
   gap in bases attached, so it can never be read as an overlap.

3. **Contig naming.** The GTF says ``chr15``; the summary says ``NC_000015.10``.
   Neither spelling contains the other, and comparing them raw fails by finding
   nothing — indistinguishable from "this variant is intergenic". The mapping is
   resolved once, explicitly, and exposed as :data:`REFSEQ_ACCESSION_TO_CONTIG`
   so a test can assert it rather than infer it from a miss. Every gene present
   in both files is cross-checked, and a disagreement raises at construction.

4. **Coordinate conventions.** GTF ``start``/``end`` are 1-based and inclusive at
   both ends. VCF ``POS`` is 1-based. This module therefore does **no half-open
   conversion anywhere**: intervals stay 1-based-inclusive from parse to answer,
   which removes the entire off-by-one class rather than managing it. A variant's
   span is ``[POS, POS + len(REF) - 1]``, so a deletion beginning one base before
   a gene still overlaps it.

5. **Symbol drift.** The gene clinicians still call *CENPJ* appears in MANE v1.5
   only as **CPAP** (ENSG00000151849, HGNC:17272). ``CASC5`` is likewise only
   ``KNL1``. Every :class:`GeneInterval` therefore carries the Ensembl gene ID,
   the HGNC ID and the NCBI GeneID alongside the symbol, so a downstream join
   never depends on the symbol string; and :meth:`ManeGeneIndex.resolve_panel`
   *reports* the panel genes it could not locate instead of dropping them.

What this module cannot say
---------------------------

``ConsequenceAnnotation.impact`` is a required, non-optional field, and this
adapter computes no impact. Every value of the enum would be a fabrication with
measurable downstream consequences — see :data:`IMPACT_NOT_ASSESSED_REMEDIATION`
and :class:`ManeGeneAdapter`, which refuses to fabricate one and instead fails
closed with the exact model change needed. There is deliberately no override: a
setting that supplied a stand-in severity would be the same fabrication behind a
switch, and it would enter ranking as evidence while the adapter still declared
``synthetic = False``.
"""

from __future__ import annotations

import gzip
import struct
from bisect import bisect_right
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from re import compile as _re_compile
from types import MappingProxyType
from typing import Final, cast

from pydantic import ValidationError

from mva.annotation.base import ConsequenceAdapter, is_synthetic
from mva.determinism import hash_file
from mva.errors import AdapterUnavailableError, GenomeBuildMismatchError
from mva.models.genome import GenomeBuild, contig_sort_key, normalise_contig
from mva.models.variant import ConsequenceAnnotation, ImpactSeverity

# --------------------------------------------------------------------------- identity

#: Adapter identity, stamped on every EvidenceItem this adapter justifies. It
#: names the *operation* — an interval join against MANE — and deliberately does
#: not resemble the name of a variant-effect predictor, because a reader of a
#: report footer must not mistake a locus lookup for a consequence prediction.
MANE_ADAPTER_NAME: Final = "mane-interval-join"

#: The two MANE_status values a summary row can carry, verbatim.
MANE_SELECT_STATUS: Final = "MANE Select"
MANE_PLUS_CLINICAL_STATUS: Final = "MANE Plus Clinical"

#: SO:0001564. The most specific thing a pure interval join is entitled to say:
#: the variant lies within the span of this gene. It asserts *location* and makes
#: no claim about the molecular effect, which is why it is not one of the terms
#: ``prioritization.scoring.LOF_TERMS`` recognises.
GENE_LOCUS_TERM: Final = "gene_variant"

#: SO:0001628. Used only for a nearest-gene inference, where the variant is by
#: construction *outside* every gene. Pairing this term with a gene symbol is the
#: label that keeps a proximity guess from reading as an overlap.
INTERGENIC_TERM: Final = "intergenic_variant"

#: Facts this module reads out of MANE and then cannot carry through
#: :class:`~mva.models.variant.ConsequenceAnnotation`, recorded here rather than
#: dropped in silence — an undocumented drop is indistinguishable from a parser
#: that never looked.
#:
#: * ``HGNC_ID`` / ``NCBI_GeneID`` — the stable cross-database identifiers.
#:   ``ConsequenceAnnotation`` has one identifier slot (``gene_id``) and the
#:   Ensembl gene ID takes it, because that is the key the GTF and the summary
#:   are joined on. Both other IDs remain reachable on :class:`GeneInterval` via
#:   :meth:`ManeGeneIndex.genes_for_symbol` and :meth:`ManeGeneIndex.assign`.
#: * ``strand`` — MANE states it; the model has no field for it.
#: * ``distance_bp`` — the gap of a nearest-gene inference. Carried on
#:   :class:`GeneAssignment`; in the annotation only the ``intergenic_variant``
#:   term survives to say the assignment is not an overlap.
#: * ``is_canonical`` — MANE does not assert canonicality, so it is left at the
#:   model's ``False`` default rather than being inferred from MANE Select. The
#:   fact MANE *does* assert travels in ``is_mane_select``.
UNREPRESENTABLE_MANE_FIELDS: Final[tuple[str, ...]] = (
    "HGNC_ID",
    "NCBI_GeneID",
    "strand",
    "distance_bp",
    "is_canonical",
)

# --------------------------------------------------------------------------- contigs


def _build_refseq_contig_map() -> Mapping[str, str]:
    """RefSeq assembly accession (version stripped) -> UCSC contig.

    GRCh38's primary assembly units are numbered ``NC_000001`` through
    ``NC_000024``, with 23 and 24 being X and Y; the mitochondrion is the
    separate ``NC_012920``. The accession *version* suffix is deliberately not
    part of the key: it changes when a chromosome is patched, and pinning the
    map to ``NC_000015.10`` would silently stop resolving chr15 on the release
    after the one that bumps it — failing by finding no genes, which reads
    exactly like "this variant is intergenic".

    Anything outside this set is an alt haplotype, a fix patch or an unplaced
    scaffold. MANE v1.5 places 64 genes on such sequences; they are out of scope
    for this pipeline (``normalise_contig`` refuses them too) and are counted
    rather than silently ignored.
    """
    mapping: dict[str, str] = {"NC_012920": "chrM"}
    for number in range(1, 25):
        name = {23: "chrX", 24: "chrY"}.get(number, f"chr{number}")
        mapping[f"NC_{number:06d}"] = name
    return MappingProxyType(mapping)


#: The join key that makes the GTF and the summary describe the same genome.
#: Exposed so a test can assert the mapping directly; a test that only checked
#: lookups would pass against a broken map for every gene the fixture omits.
REFSEQ_ACCESSION_TO_CONTIG: Final[Mapping[str, str]] = _build_refseq_contig_map()


def contig_for_refseq_accession(accession: str) -> str | None:
    """UCSC contig for a RefSeq accession, or ``None`` if it is not a chromosome."""
    return REFSEQ_ACCESSION_TO_CONTIG.get(accession.strip().split(".", 1)[0])


def _primary_contig(raw: str) -> str | None:
    """UCSC contig for a GTF ``seqname``, or ``None`` for alt/fix/unplaced."""
    try:
        return normalise_contig(raw)
    except ValueError:
        return None


# --------------------------------------------------------------------------- models


class GeneRelation(StrEnum):
    """How a gene came to be assigned to a coordinate.

    The two members are not interchangeable and must never be merged in output.
    ``OVERLAP`` is an observation about the MANE gene model; ``NEAREST`` is an
    inference this module made because nothing overlapped, and is only ever
    produced when a caller explicitly asks for it by distance.
    """

    OVERLAP = "overlap"
    NEAREST = "nearest"


@dataclass(frozen=True, slots=True)
class ManeTranscript:
    """One MANE transcript of a gene, as the summary states it.

    A gene can have two: the ``MANE Select`` and, for 74 genes in v1.5, an
    additional ``MANE Plus Clinical``. Both are kept. Reducing a gene to its
    Select transcript is the collapse ``knowledge/adapters/README.md`` rule 3
    forbids, and it is a real loss here: Plus Clinical exists precisely because
    the Select transcript does not carry some clinically reported variants.
    """

    transcript_id: str
    """Versioned Ensembl transcript, e.g. ``ENST00000287598.11``."""

    refseq_nuc: str
    """The matched RefSeq transcript, e.g. ``NM_001211.6``."""

    mane_status: str
    """``MANE Select`` or ``MANE Plus Clinical``, verbatim from the summary."""

    @property
    def is_select(self) -> bool:
        return self.mane_status == MANE_SELECT_STATUS

    def sort_key(self) -> tuple[int, str]:
        """Select first, then Plus Clinical, then by transcript ID."""
        return (0 if self.is_select else 1, self.transcript_id)


@dataclass(frozen=True, slots=True)
class GeneInterval:
    """One MANE gene: its span on the primary assembly, plus its stable IDs.

    ``start`` and ``end`` are 1-based and **inclusive at both ends**, exactly as
    the GTF states them. No half-open form of this interval exists anywhere in
    this module.

    The span is the GTF *gene* row, not the MANE Select transcript span from the
    summary. They differ: BUB1B's gene row is 40,160,984-40,221,137 while its
    Select transcript is 40,161,069-40,221,123. Using the transcript span would
    drop variants in the 85 bases of leading and 14 bases of trailing gene
    sequence — a small window, but a window at exactly the UTR boundary.
    """

    contig: str
    """UCSC-style, e.g. ``chr15``. Always a primary contig."""

    start: int
    end: int
    strand: str
    gene_id: str
    """Versioned Ensembl gene ID, e.g. ``ENSG00000156970.15``. Unique in MANE,
    which is what makes it the tie-breaker of the output's total order."""

    gene_symbol: str
    """The GTF's ``gene_name``, i.e. the Ensembl symbol."""

    summary_symbol: str
    """The summary's ``symbol``, i.e. the NCBI one. Usually identical to
    :attr:`gene_symbol`, but 47 genes in MANE v1.5 disagree — an Ensembl clone
    name against an NCBI ``LOC`` identifier (``RP4-534N18.5`` / ``LOC128031832``),
    for loci HGNC has not given an approved symbol. Both are kept and both resolve
    in :meth:`ManeGeneIndex.genes_for_symbol`, because picking one would make the
    other silently unfindable, and neither file is wrong."""

    gene_biotype: str
    hgnc_id: str | None
    """e.g. ``HGNC:1149``. ``None`` only if the summary omits it."""

    ncbi_gene_id: str | None
    """e.g. ``GeneID:701``."""

    transcripts: tuple[ManeTranscript, ...]
    """Every MANE transcript of this gene, Select first. Never empty."""

    @property
    def symbols(self) -> tuple[str, ...]:
        """Every symbol this release states for the gene, Ensembl first, deduplicated."""
        if self.summary_symbol == self.gene_symbol:
            return (self.gene_symbol,)
        return (self.gene_symbol, self.summary_symbol)

    @property
    def length(self) -> int:
        """Bases spanned, inclusive of both endpoints."""
        return self.end - self.start + 1

    @property
    def mane_select(self) -> ManeTranscript:
        """The MANE Select transcript; guaranteed present and first."""
        return self.transcripts[0]

    def contains(self, position: int) -> bool:
        """Whether a 1-based position lies within the gene, endpoints included."""
        return self.start <= position <= self.end

    def overlaps(self, span_start: int, span_end: int) -> bool:
        """Whether a 1-based inclusive span shares at least one base with the gene."""
        return self.start <= span_end and span_start <= self.end

    def gap_to(self, span_start: int, span_end: int) -> int:
        """Bases separating this gene from a 1-based inclusive span; 0 if they touch."""
        if self.overlaps(span_start, span_end):
            return 0
        return self.start - span_end if self.start > span_end else span_start - self.end

    def sort_key(self) -> tuple[int, int, int, str]:
        """Karyotype contig, then start, then end, then the unique gene ID."""
        return (contig_sort_key(self.contig), self.start, self.end, self.gene_id)


@dataclass(frozen=True, slots=True)
class GeneAssignment:
    """One gene assigned to one coordinate, and the basis on which it was assigned.

    ``relation`` is the whole point of this type. A caller that reads only
    ``gene`` and ignores ``relation`` will treat a proximity guess as a fact, so
    the two are carried together and the guess is never produced unless asked
    for.
    """

    gene: GeneInterval
    relation: GeneRelation
    distance_bp: int
    """Gap in bases. Always 0 for :attr:`GeneRelation.OVERLAP`."""

    @property
    def is_overlap(self) -> bool:
        return self.relation is GeneRelation.OVERLAP

    def sort_key(self) -> tuple[int, int, int, int, str]:
        """The documented total order of this module's output.

        ``(distance, contig, start, end, gene_id)``. Overlaps all have distance
        0, so within one variant's overlaps this reduces to genomic position;
        across a nearest-gene result the closest gene sorts first. ``gene_id`` is
        unique in a MANE release, so the order is total and no tie is ever broken
        by set or dict iteration, by insertion order, or by ``hash()`` (GP-30).
        """
        contig, start, end, gene_id = self.gene.sort_key()
        return (self.distance_bp, contig, start, end, gene_id)


@dataclass(frozen=True, slots=True)
class PanelResolution:
    """What a gene panel resolved to, including what it did not.

    ``missing`` exists so that a panel gene absent from MANE is a statement the
    run can make rather than a row that quietly vanished. ``CENPJ`` is the live
    example: it is a valid, widely used symbol that MANE v1.5 does not contain,
    because HGNC renamed the gene to ``CPAP``.
    """

    found: Mapping[str, tuple[GeneInterval, ...]]
    missing: tuple[str, ...]

    @property
    def is_complete(self) -> bool:
        return not self.missing

    def describe_missing(self) -> str:
        """One clause naming every unlocated symbol, for a warning or a report."""
        if not self.missing:
            return ""
        return (
            f"{len(self.missing)} panel symbol(s) absent from this MANE release: "
            f"{', '.join(self.missing)}. A symbol MANE does not carry is usually "
            "an HGNC rename (CENPJ -> CPAP, CASC5 -> KNL1), not a missing gene; "
            "resolve it by Ensembl or HGNC ID rather than by symbol."
        )


# --------------------------------------------------------------------------- parsing

#: The summary columns this module needs. Checked against the real header rather
#: than assumed by position: MANE has added columns between releases, and reading
#: ``chr_start`` out of whatever happens to be column 12 would mis-locate genes
#: rather than fail.
_REQUIRED_SUMMARY_COLUMNS: Final[tuple[str, ...]] = (
    "NCBI_GeneID",
    "Ensembl_Gene",
    "HGNC_ID",
    "symbol",
    "RefSeq_nuc",
    "Ensembl_nuc",
    "MANE_status",
    "GRCh38_chr",
    "chr_start",
    "chr_end",
    "chr_strand",
)

#: ``MANE.GRCh38.v1.5.<anything>``. A MANE distribution carries its release
#: nowhere else — see :func:`_release_version`.
_RELEASE_RE = _re_compile(r"^MANE\.(?P<build>GRCh\d+)\.v(?P<release>\d+(?:\.\d+)*)\.")

#: Placeholder MANE writes for an absent value.
_MISSING_VALUE: Final = "-"

_GTF_COLUMNS: Final = 9
_GTF_FEATURE_GENE: Final = "gene"


def _gtf_attribute(attributes: str, key: str) -> str | None:
    """Read one quoted GTF attribute without a regex.

    The gene rows are 19,363 of the GTF's 486,796 lines and the attribute column
    is scanned three times per row, so this stays a plain string scan. The
    boundary check keeps ``gene_id`` from matching inside a hypothetical
    ``havana_gene_id``, which a bare ``find`` would do.
    """
    needle = f'{key} "'
    offset = 0
    while True:
        start = attributes.find(needle, offset)
        if start < 0:
            return None
        if start == 0 or attributes[start - 1] in {" ", ";"}:
            value_start = start + len(needle)
            end = attributes.find('"', value_start)
            return None if end < 0 else attributes[value_start:end]
        offset = start + 1


@dataclass(frozen=True, slots=True)
class _GtfGene:
    """A GTF ``gene`` row, before the summary supplies its identifiers."""

    contig: str
    start: int
    end: int
    strand: str
    gene_id: str
    gene_symbol: str
    gene_biotype: str


def _parse_gtf_genes(path: Path) -> tuple[dict[str, _GtfGene], int]:
    """Read every ``gene`` row on a primary contig. Returns (genes, skipped).

    Only ``feature == "gene"`` rows are read. The transcript, exon, CDS and UTR
    rows describe structure this module makes no claim about, and the MANE
    transcript identifiers come from the summary, which states them directly.
    """
    genes: dict[str, _GtfGene] = {}
    skipped = 0
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != _GTF_COLUMNS:
                msg = (
                    f"{path.name!r} line {line_number} has {len(columns)} tab-separated "
                    f"columns, expected {_GTF_COLUMNS}. This is not a GTF; refusing to "
                    "assign genes from it. The line is not echoed (PRIV-09)."
                )
                raise AdapterUnavailableError(msg)
            if columns[2] != _GTF_FEATURE_GENE:
                continue
            contig = _primary_contig(columns[0])
            if contig is None:
                skipped += 1
                continue
            attributes = columns[8]
            gene_id = _gtf_attribute(attributes, "gene_id")
            gene_symbol = _gtf_attribute(attributes, "gene_name")
            gene_biotype = _gtf_attribute(attributes, "gene_type")
            if gene_id is None or gene_symbol is None:
                msg = (
                    f"{path.name!r} line {line_number} is a gene row without a gene_id or "
                    "gene_name attribute. Refusing to build an index in which some genes "
                    "have no identity."
                )
                raise AdapterUnavailableError(msg)
            if gene_id in genes:
                msg = (
                    f"{path.name!r} declares gene {gene_id} on more than one primary "
                    "contig. One gene cannot be in two places; refusing to index a "
                    "contradictory gene model."
                )
                raise AdapterUnavailableError(msg)
            genes[gene_id] = _GtfGene(
                contig=contig,
                start=int(columns[3]),
                end=int(columns[4]),
                strand=columns[6],
                gene_id=gene_id,
                gene_symbol=gene_symbol,
                gene_biotype=gene_biotype or "unknown",
            )
    if not genes:
        msg = (
            f"{path.name!r} yielded no gene rows on any primary contig. An empty gene "
            "model would assign no variant to any gene, which is indistinguishable from "
            "a genome of entirely intergenic variants; refusing to run that way."
        )
        raise AdapterUnavailableError(msg)
    return genes, skipped


@dataclass(frozen=True, slots=True)
class _SummaryGene:
    """The per-gene facts the summary contributes, after its rows are grouped."""

    contig: str | None
    """UCSC contig derived from the RefSeq accession, or ``None`` for a scaffold."""

    hgnc_id: str | None
    ncbi_gene_id: str | None
    symbol: str
    transcripts: tuple[ManeTranscript, ...]


def _parse_summary(path: Path) -> tuple[dict[str, _SummaryGene], int]:
    """Group the summary's transcript rows by Ensembl gene. Returns (genes, skipped).

    A gene has one row per MANE transcript, so the 74 genes with a Plus Clinical
    transcript appear twice. Both rows are kept as transcripts of one gene; a
    parser that keyed a dict on the gene would silently keep whichever came last.
    """
    rows: dict[str, list[Mapping[str, str]]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header_line = handle.readline()
        if not header_line:
            msg = f"{path.name!r} is empty; it carries no MANE summary header."
            raise AdapterUnavailableError(msg)
        header = header_line.lstrip("#").rstrip("\n").split("\t")
        absent = [name for name in _REQUIRED_SUMMARY_COLUMNS if name not in header]
        if absent:
            msg = (
                f"{path.name!r} is missing MANE summary column(s) {', '.join(absent)}. "
                "Columns are read by name because MANE has added columns between "
                "releases; reading them by position would mis-locate genes silently."
            )
            raise AdapterUnavailableError(msg)
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(header):
                msg = (
                    f"{path.name!r} has a row with {len(fields)} fields where the header "
                    f"declares {len(header)}. The row is not echoed (PRIV-09)."
                )
                raise AdapterUnavailableError(msg)
            row = dict(zip(header, fields, strict=True))
            rows.setdefault(row["Ensembl_Gene"], []).append(row)

    genes: dict[str, _SummaryGene] = {}
    skipped = 0
    for gene_id, gene_rows in rows.items():
        first = gene_rows[0]
        contig = contig_for_refseq_accession(first["GRCh38_chr"])
        if contig is None:
            skipped += 1
        transcripts = tuple(
            sorted(
                (
                    ManeTranscript(
                        transcript_id=row["Ensembl_nuc"],
                        refseq_nuc=row["RefSeq_nuc"],
                        mane_status=row["MANE_status"],
                    )
                    for row in gene_rows
                ),
                key=lambda transcript: transcript.sort_key(),
            )
        )
        genes[gene_id] = _SummaryGene(
            contig=contig,
            hgnc_id=_optional(first["HGNC_ID"]),
            ncbi_gene_id=_optional(first["NCBI_GeneID"]),
            symbol=first["symbol"],
            transcripts=transcripts,
        )
    return genes, skipped


def _optional(value: str) -> str | None:
    """Empty cells and MANE's ``-`` placeholder become ``None``, never ``""``."""
    stripped = value.strip()
    return None if stripped in {"", _MISSING_VALUE} else stripped


def _join(
    gtf_genes: Mapping[str, _GtfGene],
    summary_genes: Mapping[str, _SummaryGene],
    *,
    gtf_name: str,
    summary_name: str,
) -> tuple[GeneInterval, ...]:
    """Join the two files on the versioned Ensembl gene ID, checking they agree.

    The contig cross-check is the point of this function. It is the one place
    where the chr-prefixed and RefSeq-accession spellings of the same genome meet,
    and a mismatch there means every variant on that chromosome would be assigned
    to the wrong gene or to none. It is refused loudly rather than resolved by
    preferring one file: a run cannot be trusted if the two halves of its gene
    model disagree about where a gene is. Across MANE v1.5 all 19,299
    primary-contig genes agree, so a mismatch means the two files are from
    different releases or one of them is not what it claims to be.
    """
    intervals: list[GeneInterval] = []
    unmatched: list[str] = []
    contig_conflicts: list[str] = []
    for gene_id, gtf_gene in gtf_genes.items():
        summary_gene = summary_genes.get(gene_id)
        if summary_gene is None:
            unmatched.append(gene_id)
            continue
        if summary_gene.contig != gtf_gene.contig:
            contig_conflicts.append(gene_id)
            continue
        intervals.append(
            GeneInterval(
                contig=gtf_gene.contig,
                start=gtf_gene.start,
                end=gtf_gene.end,
                strand=gtf_gene.strand,
                gene_id=gene_id,
                gene_symbol=gtf_gene.gene_symbol,
                summary_symbol=summary_gene.symbol,
                gene_biotype=gtf_gene.gene_biotype,
                hgnc_id=summary_gene.hgnc_id,
                ncbi_gene_id=summary_gene.ncbi_gene_id,
                transcripts=summary_gene.transcripts,
            )
        )
    if contig_conflicts:
        msg = (
            f"{len(contig_conflicts)} gene(s) sit on different contigs in {gtf_name!r} "
            f"and {summary_name!r}; the first is {contig_conflicts[0]}. The GTF is "
            "chr-prefixed and the summary is RefSeq-accession, so a disagreement here "
            "means the accession map is wrong or the two files are different releases. "
            "Refusing to assign genes from a contradictory gene model."
        )
        raise AdapterUnavailableError(msg)
    if unmatched:
        msg = (
            f"{len(unmatched)} gene(s) in {gtf_name!r} have no row in {summary_name!r}; "
            f"the first is {unmatched[0]}. The two files must be the same MANE release: "
            "the summary supplies the HGNC and NCBI identifiers and the MANE transcript "
            "IDs, without which a gene has no stable identity to join on downstream."
        )
        raise AdapterUnavailableError(msg)
    return tuple(sorted(intervals, key=lambda interval: interval.sort_key()))


# --------------------------------------------------------------------------- integrity


def _verify_integrity(path: Path, *, expected_sha256: str | None, role: str) -> None:
    """Hash the bytes before a single row is read as data.

    Mirrors ``clinvar_vcf._verify_integrity`` and ``local_tables.verify_manifest``:
    hash first, open second. Fails closed when unpinned — an unpinned reference
    release quietly becomes a different file between runs, and every gene
    assignment derived from it would then be unreproducible (GP-30).
    """
    if expected_sha256 is None:
        msg = (
            f"Refusing to open the MANE {role} {path.name!r} without an integrity pin. "
            "Pass the expected sha256 (recorded in the resource manifest). Gene "
            "assignment decides which variants can ever be paired, so an unpinned "
            "gene model makes every candidate pair unreproducible."
        )
        raise AdapterUnavailableError(msg)
    actual = hash_file(path)
    if actual != expected_sha256.strip().lower():
        msg = (
            f"MANE {role} {path.name!r} failed its sha256 integrity check: expected "
            f"{expected_sha256.strip().lower()}, found {actual}. Refusing to assign "
            "genes from an unverified gene model; re-run the acquisition step and "
            "review the manifest diff. The file's contents are not echoed (PRIV-09)."
        )
        raise AdapterUnavailableError(msg)


def gzip_stored_filename(path: Path) -> str | None:
    """The original filename gzip recorded *inside* the compressed member, if any.

    RFC 1952 stores it under the ``FNAME`` flag. It matters here because it is the
    only release identifier that lives within the bytes the sha256 pin covers:
    the on-disk name can be changed by anyone, this cannot be changed without
    breaking the pin. MANE's summary carries it; the GTF, as distributed, does
    not, so it is a corroborating check and never the sole source.
    """
    with path.open("rb") as handle:
        head = handle.read(10)
        if len(head) < 10 or head[0] != 0x1F or head[1] != 0x8B:
            return None
        flags = head[3]
        if not flags & 0x08:
            return None
        if flags & 0x04:
            raw_extra_length = handle.read(2)
            if len(raw_extra_length) < 2:
                return None
            (extra_length,) = struct.unpack("<H", raw_extra_length)
            handle.read(extra_length)
        name = bytearray()
        while (byte := handle.read(1)) not in {b"", b"\x00"}:
            name.extend(byte)
        return name.decode("latin-1") or None


def _release_token(name: str) -> tuple[str, str] | None:
    """``(build, release)`` from a MANE distribution filename, or ``None``."""
    match = _RELEASE_RE.match(name)
    if match is None:
        return None
    return match.group("build"), match.group("release")


def _release_version(gtf_path: Path, summary_path: Path, *, build: GenomeBuild) -> str:
    """The release identity, read from the files rather than invented.

    A MANE distribution states its version in exactly one place: the release
    prefix of its filenames, ``MANE.GRCh38.v1.5.``. The GTF carries no comment or
    header line at all (v1.5 has zero lines beginning ``#`` in 486,796 lines) and
    the summary's single header line names columns, not a version. So the prefix
    is the artifact's own self-description, and it is corroborated three ways
    rather than trusted once:

    * both files must carry the same release prefix — a v1.3 summary paired with
      a v1.5 GTF is caught here rather than surfacing later as missing genes;
    * where gzip stored the original filename inside the compressed bytes, that
      name must agree too, which puts one copy of the release string under the
      sha256 pin;
    * the build the prefix declares must be the build the caller is joining
      against, because a GRCh37 coordinate against a GRCh38 gene model finds the
      wrong gene or no gene, silently (GP-11).
    """
    tokens: dict[str, tuple[str, str]] = {}
    for label, path in (("GTF", gtf_path), ("summary", summary_path)):
        token = _release_token(path.name)
        if token is None:
            msg = (
                f"MANE {label} {path.name!r} does not begin with a MANE release prefix "
                "('MANE.<build>.v<release>.'). That prefix is the only place a MANE "
                "distribution states its own version — the GTF has no header line and "
                "the summary's header names columns — so a renamed file has no "
                "identifiable release and every claim from it would be uncitable."
            )
            raise AdapterUnavailableError(msg)
        tokens[label] = token
        stored = gzip_stored_filename(path)
        if stored is not None:
            stored_token = _release_token(stored)
            if stored_token is not None and stored_token != token:
                msg = (
                    f"MANE {label} {path.name!r} was renamed: gzip records the release "
                    f"'{stored_token[0]}.v{stored_token[1]}' inside the compressed bytes "
                    f"but the filename claims '{token[0]}.v{token[1]}'. Refusing to stamp "
                    "reports with a version the artifact contradicts."
                )
                raise AdapterUnavailableError(msg)

    if tokens["GTF"] != tokens["summary"]:
        gtf_build, gtf_release = tokens["GTF"]
        summary_build, summary_release = tokens["summary"]
        msg = (
            f"MANE release mismatch: the GTF is {gtf_build} v{gtf_release} and the "
            f"summary is {summary_build} v{summary_release}. The two files are joined "
            "on versioned Ensembl gene IDs, which change between releases, so mixing "
            "them drops genes rather than failing."
        )
        raise AdapterUnavailableError(msg)

    declared_build, release = tokens["GTF"]
    try:
        parsed_build = GenomeBuild.parse(declared_build)
    except ValueError as exc:
        msg = (
            f"MANE release declares an unrecognised assembly {declared_build!r}. A "
            "coordinate without an establishable build is invalid (GP-11)."
        )
        raise AdapterUnavailableError(msg) from exc
    if parsed_build is not build:
        msg = (
            f"MANE release is {parsed_build.value} but this index was asked to join "
            f"{build.value} coordinates. Cross-build gene assignment mis-locates every "
            "variant by megabases while looking entirely successful: the GRCh37 "
            "coordinate of TRIP13 is around chr5:600,000, which in GRCh38 is intergenic "
            "and 293 kb from the gene (GP-11)."
        )
        raise AdapterUnavailableError(msg)
    return f"{parsed_build.value}-MANE-v{release}"


# --------------------------------------------------------------------------- index


@dataclass(frozen=True, slots=True)
class _ContigIndex:
    """One contig's genes as parallel arrays, sorted by start.

    ``max_ends[i]`` is the running maximum of ``ends[0..i]``. It is what makes a
    backward scan terminable: once it drops below the query's start, no earlier
    gene can reach the query either, so the scan stops instead of walking the
    contig. Without it a nested gene — ANKRD63 sits entirely inside PLCB2 — would
    force either a full scan or a wrong early exit.
    """

    starts: tuple[int, ...]
    ends: tuple[int, ...]
    max_ends: tuple[int, ...]
    genes: tuple[GeneInterval, ...]


def _build_contig_index(genes: Sequence[GeneInterval]) -> _ContigIndex:
    ordered = sorted(genes, key=lambda gene: (gene.start, gene.end, gene.gene_id))
    starts = tuple(gene.start for gene in ordered)
    ends = tuple(gene.end for gene in ordered)
    running: list[int] = []
    highest = 0
    for end in ends:
        highest = max(highest, end)
        running.append(highest)
    return _ContigIndex(starts=starts, ends=ends, max_ends=tuple(running), genes=tuple(ordered))


class ManeGeneIndex:
    """A per-contig interval index over the gene rows of a MANE GRCh38 release.

    Built once, on this explicitly constructed object. There is deliberately no
    module-level cache and no ``lru_cache`` over a function taking paths: a
    process-global gene model is state one stage can mutate on behalf of another,
    it makes construction order observable in results, and it defeats the layering
    the adapter boundary exists to enforce. Two indexes over two releases can
    coexist in one process, which is what makes a release comparison testable.

    Lookup is a binary search plus a short backward scan, not a linear pass: the
    intended input is a whole-genome VCF of 4-5 million records, where a scan of
    19,299 genes per variant would be roughly 10^11 comparisons.
    """

    def __init__(
        self,
        gtf_path: Path,
        summary_path: Path,
        *,
        expected_gtf_sha256: str | None = None,
        expected_summary_sha256: str | None = None,
        build: GenomeBuild = GenomeBuild.GRCH38,
    ) -> None:
        """Load, verify and index a MANE release.

        Args:
            gtf_path: ``MANE.<build>.v<release>.ensembl_genomic.gtf.gz``. Only its
                ``gene`` feature rows are read.
            summary_path: ``MANE.<build>.v<release>.summary.txt.gz``. Supplies the
                HGNC and NCBI identifiers and the MANE transcript IDs, and is the
                RefSeq-accession cross-check on every gene's contig.
            expected_gtf_sha256: Integrity pin from the resource manifest.
                Required: construction fails closed without it.
            expected_summary_sha256: Integrity pin for the summary. Also required —
                the summary is read as data too, and an unverified identifier table
                would put unverified gene IDs into every downstream join.
            build: The assembly the caller's variant IDs are in. The release's own
                filename prefix must declare the same build (GP-11).

        Raises:
            AdapterUnavailableError: a file is missing, unpinned, fails its hash,
                is not the shape MANE distributes, declares a different release or
                build from its partner, or describes a gene model whose two halves
                disagree about where a gene is.
        """
        for label, path in (("GTF", gtf_path), ("summary", summary_path)):
            if not path.is_file():
                msg = (
                    f"MANE {label} {path.as_posix()!r} not found. This adapter reads a "
                    "pre-downloaded, hash-pinned local release; it never fetches "
                    "anything (PRIV-05)."
                )
                raise AdapterUnavailableError(msg)
        _verify_integrity(gtf_path, expected_sha256=expected_gtf_sha256, role="GTF")
        _verify_integrity(summary_path, expected_sha256=expected_summary_sha256, role="summary")

        self._gtf_path = gtf_path
        self._summary_path = summary_path
        self._build = build
        self._version = _release_version(gtf_path, summary_path, build=build)

        gtf_genes, gtf_skipped = _parse_gtf_genes(gtf_path)
        summary_genes, summary_skipped = _parse_summary(summary_path)
        self._genes = _join(
            gtf_genes,
            summary_genes,
            gtf_name=gtf_path.name,
            summary_name=summary_path.name,
        )
        self._skipped_gtf_genes = gtf_skipped
        self._skipped_summary_genes = summary_skipped

        by_contig: dict[str, list[GeneInterval]] = {}
        by_symbol: dict[str, list[GeneInterval]] = {}
        by_gene_id: dict[str, GeneInterval] = {}
        for gene in self._genes:
            by_contig.setdefault(gene.contig, []).append(gene)
            for symbol in gene.symbols:
                by_symbol.setdefault(symbol, []).append(gene)
            by_gene_id[gene.gene_id] = gene
        self._index: Mapping[str, _ContigIndex] = MappingProxyType(
            {contig: _build_contig_index(genes) for contig, genes in by_contig.items()}
        )
        self._by_symbol: Mapping[str, tuple[GeneInterval, ...]] = MappingProxyType(
            {symbol: tuple(genes) for symbol, genes in by_symbol.items()}
        )
        self._by_gene_id: Mapping[str, GeneInterval] = MappingProxyType(by_gene_id)

    # ------------------------------------------------------------------ identity

    @property
    def version(self) -> str:
        """The release, read out of the distribution: e.g. ``GRCh38-MANE-v1.5``."""
        return self._version

    @property
    def build(self) -> GenomeBuild:
        return self._build

    @property
    def gtf_path(self) -> Path:
        return self._gtf_path

    @property
    def summary_path(self) -> Path:
        return self._summary_path

    @property
    def genes(self) -> tuple[GeneInterval, ...]:
        """Every indexed gene, in the module's documented total order."""
        return self._genes

    @property
    def gene_count(self) -> int:
        return len(self._genes)

    @property
    def contigs(self) -> tuple[str, ...]:
        """Indexed contigs, in karyotype order."""
        return tuple(sorted(self._index, key=contig_sort_key))

    @property
    def skipped_gene_counts(self) -> Mapping[str, int]:
        """Genes each file placed off the primary assembly, counted not hidden.

        MANE v1.5 puts 64 genes on alt haplotypes, fix patches and unplaced
        scaffolds — ``MUC2`` on ``chr11_KQ759759v2_fix``, ``GSTT1`` on
        ``chr22_KI270879v1_alt``. This pipeline does not reason about those
        sequences, and both files must drop the same ones; the counts are exposed
        so a test can assert they agree instead of inferring it from a miss.
        """
        return MappingProxyType(
            {"gtf": self._skipped_gtf_genes, "summary": self._skipped_summary_genes}
        )

    # ------------------------------------------------------------------ lookup

    def genes_at(self, contig: str, start: int, end: int | None = None) -> tuple[GeneInterval, ...]:
        """Every gene overlapping a 1-based inclusive span, in total order.

        ``end`` defaults to ``start``, i.e. a single base. Both endpoints are
        inclusive, matching the GTF and VCF conventions exactly: ``genes_at(c, g.start)``
        and ``genes_at(c, g.end)`` both return the gene, and one base either side
        does not.

        Returns an empty tuple for an intergenic span and for a contig the release
        holds no genes on. The caller distinguishes "no gene here" from "we did not
        look" through :meth:`assign_variants`, which omits the variant entirely.
        """
        span_end = start if end is None else end
        if span_end < start:
            msg = f"span end {span_end} precedes span start {start}"
            raise ValueError(msg)
        contig_index = self._index.get(normalise_contig(contig))
        if contig_index is None:
            return ()
        return tuple(sorted(self._stab(contig_index, start, span_end), key=_gene_order))

    def assign(
        self,
        contig: str,
        start: int,
        end: int | None = None,
        *,
        nearest_within_bp: int | None = None,
    ) -> tuple[GeneAssignment, ...]:
        """Genes for a span, each labelled with how it was assigned.

        Overlapping genes are returned as :attr:`GeneRelation.OVERLAP` with
        ``distance_bp == 0``.

        ``nearest_within_bp`` is **off by default** and changes nothing when a gene
        overlaps. Only when nothing overlaps does it look outward, and then only as
        far as it is told; every gene it finds is returned as
        :attr:`GeneRelation.NEAREST` with the true gap in bases. That labelling is
        not decoration. A GRCh37 TRIP13 coordinate around chr5:600,000 is 12 kb
        from CEP72 in GRCh38, so a nearest-gene lookup there confidently returns
        the wrong gene; only the ``NEAREST`` label and the distance let a reader
        see that the assignment is a guess about a variant that is in no gene at
        all. All genes tied at the minimum distance are returned, never one.
        """
        if nearest_within_bp is not None and nearest_within_bp < 0:
            msg = "nearest_within_bp must not be negative."
            raise ValueError(msg)
        span_end = start if end is None else end
        overlaps = self.genes_at(contig, start, span_end)
        if overlaps:
            return tuple(
                GeneAssignment(gene=gene, relation=GeneRelation.OVERLAP, distance_bp=0)
                for gene in overlaps
            )
        if nearest_within_bp is None:
            return ()
        contig_index = self._index.get(normalise_contig(contig))
        if contig_index is None:
            return ()
        window = self._stab(
            contig_index, max(1, start - nearest_within_bp), span_end + nearest_within_bp
        )
        candidates = [
            GeneAssignment(
                gene=gene,
                relation=GeneRelation.NEAREST,
                distance_bp=gene.gap_to(start, span_end),
            )
            for gene in window
        ]
        candidates = [
            candidate for candidate in candidates if candidate.distance_bp <= nearest_within_bp
        ]
        if not candidates:
            return ()
        closest = min(candidate.distance_bp for candidate in candidates)
        return tuple(
            sorted(
                (candidate for candidate in candidates if candidate.distance_bp == closest),
                key=lambda candidate: candidate.sort_key(),
            )
        )

    def assign_variants(
        self, variant_ids: Sequence[str], *, nearest_within_bp: int | None = None
    ) -> Mapping[str, tuple[GeneAssignment, ...]]:
        """Assign genes to canonical variant IDs, omitting the ones with no gene.

        Keys appear **only** for variants that landed in (or, if asked, near) at
        least one gene. A variant in an intergenic region is absent from the
        mapping rather than present with an empty tuple, so "this variant is in no
        gene" and "this variant was never looked up" stay distinguishable at the
        type level (GP-14).

        A variant's span is its REF span, ``[POS, POS + len(REF) - 1]``, so a
        deletion beginning just outside a gene and reaching into it is assigned to
        that gene.

        Raises:
            GenomeBuildMismatchError: a variant ID names a different build from the
                release. Refused rather than skipped: a silent miss would read as
                "every variant in this batch is intergenic".
            ValueError: an ID is not the canonical
                ``{build}:{contig}:{pos}:{ref}:{alt}`` form, or names a contig this
                pipeline does not reason about.
        """
        assigned: dict[str, tuple[GeneAssignment, ...]] = {}
        for variant_id in _unique(variant_ids):
            contig, position, ref = self._parse_variant_id(variant_id)
            genes = self.assign(
                contig,
                position,
                position + len(ref) - 1,
                nearest_within_bp=nearest_within_bp,
            )
            if genes:
                assigned[variant_id] = genes
        return assigned

    # ------------------------------------------------------------------ by identity

    def genes_for_symbol(self, symbol: str) -> tuple[GeneInterval, ...]:
        """Genes carrying this exact MANE symbol. Empty tuple if MANE has none.

        Both symbols the release states for a gene resolve here — the GTF's
        Ensembl ``gene_name`` and the summary's NCBI ``symbol`` — because 47 genes
        in v1.5 spell themselves differently in the two files and neither spelling
        is wrong. This is not an alias table: every string accepted was read out
        of the pinned release for that exact gene.

        Exact match only. There is no alias table and no case-insensitive fallback:
        guessing that ``CENPJ`` means ``CPAP`` is a curation claim this module has
        no source for, and a wrong guess assigns variants to the wrong gene. Use
        :meth:`gene_for_id` with an Ensembl or HGNC ID when the symbol may have
        drifted.
        """
        return self._by_symbol.get(symbol, ())

    def gene_for_id(self, gene_id: str) -> GeneInterval | None:
        """Look up by versioned Ensembl gene ID, e.g. ``ENSG00000156970.15``."""
        return self._by_gene_id.get(gene_id)

    def resolve_panel(self, symbols: Sequence[str]) -> PanelResolution:
        """Resolve a gene panel, reporting what could not be located.

        The missing list is the product here. A panel symbol MANE does not carry
        is almost always an HGNC rename, and silently dropping it removes the gene
        from the hypothesis space without anything in the run saying so.
        """
        found: dict[str, tuple[GeneInterval, ...]] = {}
        missing: list[str] = []
        for symbol in _unique(symbols):
            genes = self.genes_for_symbol(symbol)
            if genes:
                found[symbol] = genes
            else:
                missing.append(symbol)
        return PanelResolution(found=MappingProxyType(found), missing=tuple(missing))

    # ------------------------------------------------------------------ internals

    @staticmethod
    def _stab(index: _ContigIndex, start: int, end: int) -> list[GeneInterval]:
        """Genes overlapping ``[start, end]``, 1-based inclusive.

        ``bisect_right(starts, end)`` bounds the candidates to genes that begin at
        or before the span's last base. Walking backwards from there, a gene is a
        hit when its own end reaches the span's first base, and the walk stops as
        soon as ``max_ends`` says no earlier gene can.
        """
        hits: list[GeneInterval] = []
        position = bisect_right(index.starts, end) - 1
        max_ends = index.max_ends
        ends = index.ends
        while position >= 0:
            if max_ends[position] < start:
                break
            if ends[position] >= start:
                hits.append(index.genes[position])
            position -= 1
        return hits

    def _parse_variant_id(self, variant_id: str) -> tuple[str, int, str]:
        """Decompose a canonical variant ID into ``(contig, position, ref)``.

        Error text carries field counts and build names only. The coordinate is
        patient data and is never echoed (PRIV-09).
        """
        parts = variant_id.split(":")
        if len(parts) != 5:
            msg = (
                f"Variant ID has {len(parts)} colon-separated fields, expected 5 "
                "({build}:{contig}:{pos}:{ref}:{alt}). The value is not echoed (PRIV-09)."
            )
            raise ValueError(msg)
        build_token, contig_token, position_token, ref, _alt = parts
        try:
            build = GenomeBuild.parse(build_token)
        except ValueError as exc:
            msg = (
                "Variant ID does not begin with a recognised genome build. A coordinate "
                "without a build is invalid (GP-11); the value is not echoed (PRIV-09)."
            )
            raise ValueError(msg) from exc
        if build is not self._build:
            msg = (
                f"Variant ID is {build.value} but this MANE index is "
                f"{self._build.value}. Cross-build gene assignment is refused rather "
                "than silently missed: the same locus differs by megabases between "
                "assemblies, so the join would return the wrong gene or none at all "
                "while looking successful. The coordinate is not echoed (PRIV-09)."
            )
            raise GenomeBuildMismatchError(msg)
        if not position_token.isdigit():
            msg = "Variant ID position is not a positive integer (PRIV-09: not echoed)."
            raise ValueError(msg)
        if not ref:
            msg = "Variant ID carries an empty REF allele (PRIV-09: not echoed)."
            raise ValueError(msg)
        return normalise_contig(contig_token), int(position_token), ref


def _gene_order(gene: GeneInterval) -> tuple[int, int, int, str]:
    return gene.sort_key()


def _unique(values: Iterable[str]) -> tuple[str, ...]:
    """Deduplicate while preserving first-seen order; never a set (GP-30)."""
    seen: dict[str, None] = {}
    for value in values:
        seen.setdefault(value, None)
    return tuple(seen)


# --------------------------------------------------------------------------- adapter

#: Why this adapter cannot produce a ``ConsequenceAnnotation`` today, and the
#: change that would let it. Kept as a constant so the failure message, the module
#: docstring and the test that asserts the refusal all state one thing.
#:
#: ``ConsequenceAnnotation.impact`` is a required ``ImpactSeverity``. An interval
#: join computes no impact, and every member of that enum is a fabrication with a
#: measurable, checkable downstream effect:
#:
#: * ``MODIFIER`` scores 0.05 in ``prioritization.scoring._IMPACT_BASE`` and is a
#:   member of ``prioritization.filters.BENIGN_IMPACTS``, so it would attach a
#:   ``benign_consequence`` flag to a variant nobody assessed — absence of
#:   information rendered as negative information, which is exactly GP-14.
#: * ``LOW`` is also in ``BENIGN_IMPACTS`` and does the same thing at 0.15.
#: * ``MODERATE`` (0.50) and ``HIGH`` (0.90) invent severity, and ``HIGH``
#:   additionally makes every lone heterozygote a single-variant candidate in
#:   ``pairing._wants_single_candidate``.
#:
#: The pipeline already handles the honest answer correctly and needs no new
#: logic for it: ``VariantRecord.worst_impact_for_gene`` returning ``None`` makes
#: ``scoring._variant_consequence`` score ``_IMPACT_UNANNOTATED`` with the
#: rationale "no consequence annotation for this gene (unknown, not benign)", and
#: ``pairing._IMPACT_ORDER_UNKNOWN`` already sorts an unannotated variant with
#: ``LOW`` on the stated grounds that absence is not evidence of benignity.
IMPACT_NOT_ASSESSED_REMEDIATION: Final = (
    "This adapter assigns genes by interval join and computes no molecular "
    "consequence, so it has no impact severity to report — and ConsequenceAnnotation "
    "requires one. Fabricating MODIFIER or LOW would attach a 'benign_consequence' "
    "flag to a variant nobody assessed (GP-14); MODERATE or HIGH would invent "
    "severity. Remediation, in three files: (1) in src/mva/models/variant.py make the "
    "field optional -- `impact: ImpactSeverity | None = None` -- and skip None in "
    "`worst_impact_for_gene` so it keeps returning None for an unassessed gene; "
    "(2) in src/mva/annotation/service.py::_consequence_evidence handle "
    "`annotation.impact is None` (NEUTRAL direction, INSUFFICIENT strength, "
    "'not assessed' in the claim); (3) in src/mva/evidence/store.py the "
    "consequence row must write None rather than `csq.impact.value`. This adapter "
    "offers no way to override the refusal: a configuration knob supplying a "
    "stand-in severity is the same fabrication with a switch in front of it, and it "
    "would reach ranking as evidence while the adapter still declared itself real."
)


def model_accepts_unassessed_impact() -> bool:
    """Whether ``ConsequenceAnnotation`` can represent 'impact not assessed'.

    Public so a composition root can ask before binding: while this is False,
    :class:`ManeGeneAdapter` refuses to construct and the consequence slot needs a
    real effect predictor or nothing.

    Probed by construction rather than by reading pydantic's field internals, so
    that the moment the model is widened this adapter starts emitting the honest
    annotation with no change here. Called once per adapter, not per variant.
    """
    try:
        ConsequenceAnnotation.model_validate(
            {
                "gene_symbol": "PROBE",
                "gene_id": None,
                "transcript_id": "PROBE",
                "consequence_terms": [GENE_LOCUS_TERM],
                "impact": None,
            }
        )
    except ValidationError:
        return False
    return True


class ManeGeneAdapter:
    """A ``ConsequenceAdapter`` that reports gene locus only, from MANE intervals.

    This fills the consequence slot because that is the only place in the adapter
    boundary whose return type carries a gene symbol, and ``gene_symbol`` on a
    ``ConsequenceAnnotation`` is the single field that reaches
    ``VariantRecord.gene_symbols`` and therefore ``generate_pairs``. A separate
    ``GeneLocusAdapter`` Protocol would be the cleaner long-term shape — gene
    assignment is a locational fact, not a prediction — but it would return values
    nothing currently reads, so pairing would stay starved. See the report
    accompanying this module.

    It is a *partial* consequence annotation and says so: the gene, its stable
    identifiers and its MANE transcripts are real facts from a hash-pinned
    release; the molecular effect is simply not assessed. It composes with a real
    effect predictor rather than competing with it — see
    :class:`GeneBackfillConsequenceAdapter`.
    """

    def __init__(self, index: ManeGeneIndex, *, nearest_within_bp: int | None = None) -> None:
        """Bind an index to the consequence slot, or refuse to.

        Args:
            index: The loaded, verified MANE gene model.
            nearest_within_bp: Off by default. When set, a variant in no gene is
                annotated against the nearest gene(s) within this many bases, with
                the ``intergenic_variant`` term marking the assignment as an
                inference rather than an overlap.

        Raises:
            AdapterUnavailableError: ``ConsequenceAnnotation`` cannot represent an
                unassessed impact severity. There is deliberately no way to talk
                this adapter past that: the message carries the exact remediation,
                and a configuration knob that supplied a stand-in severity would be
                a fabricated prediction with a switch in front of it. See
                :data:`IMPACT_NOT_ASSESSED_REMEDIATION`.
        """
        if not model_accepts_unassessed_impact():
            raise AdapterUnavailableError(IMPACT_NOT_ASSESSED_REMEDIATION)
        self._index = index
        self._nearest_within_bp = nearest_within_bp
        self._version = index.version

    # ------------------------------------------------------------------ identity

    @property
    def name(self) -> str:
        return MANE_ADAPTER_NAME

    @property
    def version(self) -> str:
        """The MANE release the assignments came from, e.g. ``GRCh38-MANE-v1.5``."""
        return self._version

    @property
    def synthetic(self) -> bool:
        """False, deliberately (GP-20).

        ``is_synthetic`` fails closed, so every adapter is a mock until it says
        otherwise. This one says otherwise: the gene spans, symbols, Ensembl and
        HGNC identifiers and MANE transcripts all come from an actual MANE release
        whose bytes were checked against a pinned sha256 before the first row was
        read. It is not a claim that the *consequence* is real — no consequence is
        computed, which is the point of :data:`IMPACT_NOT_ASSESSED_REMEDIATION`.
        """
        return False

    @property
    def index(self) -> ManeGeneIndex:
        return self._index

    # ------------------------------------------------------------------ lookup

    def annotate(
        self, variant_ids: Sequence[str]
    ) -> Mapping[str, tuple[ConsequenceAnnotation, ...]]:
        """One annotation per (overlapping gene, MANE transcript), in a total order.

        Both MANE transcripts of a gene are emitted where the release holds two, so
        a Plus Clinical transcript is never collapsed into its Select sibling.

        Keys appear only for variants assigned at least one gene: a variant in an
        intergenic region is missing from the mapping, not present and empty (GP-14).
        """
        annotations: dict[str, tuple[ConsequenceAnnotation, ...]] = {}
        assigned = self._index.assign_variants(
            variant_ids, nearest_within_bp=self._nearest_within_bp
        )
        for variant_id, assignments in assigned.items():
            annotations[variant_id] = tuple(
                self._annotation(assignment, transcript)
                for assignment in assignments
                for transcript in assignment.gene.transcripts
            )
        return annotations

    def _annotation(
        self, assignment: GeneAssignment, transcript: ManeTranscript
    ) -> ConsequenceAnnotation:
        gene = assignment.gene
        return ConsequenceAnnotation(
            gene_symbol=gene.gene_symbol,
            gene_id=gene.gene_id,
            transcript_id=transcript.transcript_id,
            transcript_biotype=gene.gene_biotype,
            # MANE asserts a matched transcript pair, not canonicality, so this is
            # left at the model default rather than inferred from MANE Select.
            is_canonical=False,
            is_mane_select=transcript.is_select,
            consequence_terms=((GENE_LOCUS_TERM,) if assignment.is_overlap else (INTERGENIC_TERM,)),
            # NOT ASSESSED. An interval join computes no molecular consequence, and
            # every member of ImpactSeverity would be a prediction this adapter did
            # not make. Construction already refused unless the model can hold it.
            impact=cast(ImpactSeverity, None),
            source_tool=MANE_ADAPTER_NAME,
            source_tool_version=self._version,
        )


class GeneBackfillConsequenceAdapter:
    """Compose a real effect predictor with this gene-locus adapter.

    ``annotation.service`` binds exactly one consequence adapter and replaces
    ``VariantRecord.consequences`` wholesale from it, so a run cannot have both a
    variant-effect predictor and a gene-locus fallback without composing them into
    one object. This is that composition, and it is a backfill rather than a merge:
    the primary adapter's answer is used verbatim wherever it has one, and this
    adapter contributes only for variants the primary omitted entirely. Merging
    per-gene would create annotations that look like one tool's opinion but are two.

    ``synthetic`` is the OR of both members, not the primary's: a real predictor
    backed by a synthetic fallback is still partly synthetic, and ``is_synthetic``
    fails closed by design (GP-20).
    """

    def __init__(self, primary: ConsequenceAdapter, fallback: ConsequenceAdapter) -> None:
        self._primary = primary
        self._fallback = fallback

    @property
    def name(self) -> str:
        return f"{self._primary.name}+{self._fallback.name}"

    @property
    def version(self) -> str:
        return f"{self._primary.version}+{self._fallback.version}"

    @property
    def synthetic(self) -> bool:
        return is_synthetic(self._primary) or is_synthetic(self._fallback)

    @property
    def primary(self) -> ConsequenceAdapter:
        return self._primary

    @property
    def fallback(self) -> ConsequenceAdapter:
        return self._fallback

    def annotate(
        self, variant_ids: Sequence[str]
    ) -> Mapping[str, tuple[ConsequenceAnnotation, ...]]:
        """Primary answers, backfilled with gene-locus answers where absent."""
        primary = self._primary.annotate(variant_ids)
        uncovered = tuple(
            variant_id for variant_id in _unique(variant_ids) if variant_id not in primary
        )
        if not uncovered:
            return dict(primary)
        merged: dict[str, tuple[ConsequenceAnnotation, ...]] = dict(primary)
        merged.update(self._fallback.annotate(uncovered))
        return merged
