"""Genome build and coordinate primitives.

Core invariant of this project: **a genomic coordinate without a genome build is
invalid**. GRCh37 and GRCh38 positions differ by megabases for the same locus;
silently mixing them produces confidently wrong answers. Every coordinate-bearing
model therefore carries its build, and cross-build comparison raises rather than
guesses.
"""

from __future__ import annotations

import re
from enum import StrEnum
from typing import Self

from pydantic import Field, field_validator, model_validator

from mva.models.base import FrozenModel


class GenomeBuild(StrEnum):
    """Supported reference assemblies.

    The MVA Hackathon 2026 challenge mandates GRCh38 (verified from the Space's
    ``tabs/submit_track1.py``). GRCh37 is modelled so that mismatches can be
    *detected and rejected*, not so that they can be silently accepted.
    """

    GRCH38 = "GRCh38"
    GRCH37 = "GRCh37"

    @property
    def aliases(self) -> frozenset[str]:
        """Accepted spellings seen in real VCF headers."""
        if self is GenomeBuild.GRCH38:
            return frozenset({"grch38", "hg38", "GRCh38", "hg38.p13", "b38", "38"})
        return frozenset({"grch37", "hg19", "GRCh37", "b37", "37"})

    @classmethod
    def parse(cls, raw: str) -> GenomeBuild:
        """Resolve a header string to a build, or raise.

        Deliberately strict: an unrecognised assembly string is an error, because
        the alternative is assuming GRCh38 and mis-locating every variant.
        """
        token = raw.strip().lower()
        for build in cls:
            if token in {alias.lower() for alias in build.aliases}:
                return build
        msg = f"Unrecognised genome build {raw!r}; expected one of {[b.value for b in cls]}"
        raise ValueError(msg)


class ContigStyle(StrEnum):
    """Chromosome naming convention.

    This is not cosmetic. The challenge scorer compares chromosome strings raw
    (``chrom.strip()`` only, no prefix normalisation), so emitting ``1`` where the
    answer key holds ``chr1`` scores zero while looking correct to a human reader.
    """

    UCSC = "ucsc"  # chr1, chr2, chrX, chrM
    ENSEMBL = "ensembl"  # 1, 2, X, MT


_CONTIG_RE = re.compile(r"^(?:chr)?([0-9]{1,2}|[XYxy]|M|MT|mt)$")

#: Canonical nuclear + mitochondrial contigs, in karyotype order.
CANONICAL_CONTIGS: tuple[str, ...] = tuple(
    [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
)
_CONTIG_ORDER: dict[str, int] = {name: i for i, name in enumerate(CANONICAL_CONTIGS)}


def normalise_contig(raw: str, style: ContigStyle = ContigStyle.UCSC) -> str:
    """Convert a contig name to the requested style.

    Raises for non-canonical contigs (alts, decoys, unplaced scaffolds) rather than
    passing them through: a candidate variant on ``chr1_KI270766v1_alt`` is not
    something this pipeline is entitled to reason about.
    """
    match = _CONTIG_RE.match(raw.strip())
    if match is None:
        msg = f"Non-canonical contig {raw!r} (alt/decoy/unplaced contigs are not supported)"
        raise ValueError(msg)
    core = match.group(1).upper()
    if core in {"M", "MT"}:
        core = "M" if style is ContigStyle.UCSC else "MT"
    return f"chr{core}" if style is ContigStyle.UCSC else core


def contig_sort_key(contig: str) -> int:
    """Karyotype ordering (chr1..chr22, X, Y, M) for deterministic output."""
    return _CONTIG_ORDER.get(normalise_contig(contig), len(CANONICAL_CONTIGS))


_ALLELE_RE = re.compile(r"^[ACGTN]+$")


class GenomicCoordinate(FrozenModel):
    """A build-anchored, 1-based, VCF-convention position with alleles."""

    build: GenomeBuild
    contig: str = Field(description="UCSC-style contig, e.g. 'chr15'")
    position: int = Field(gt=0, description="1-based VCF POS")
    ref: str = Field(min_length=1)
    alt: str = Field(min_length=1)

    @field_validator("contig")
    @classmethod
    def _canonical_contig(cls, value: str) -> str:
        return normalise_contig(value, ContigStyle.UCSC)

    @field_validator("ref", "alt")
    @classmethod
    def _valid_allele(cls, value: str) -> str:
        upper = value.strip().upper()
        if upper in {"*", "."}:
            # Spanning-deletion and missing alleles are legal VCF but are not
            # independently rankable candidates; ingestion filters them earlier.
            return upper
        if not _ALLELE_RE.match(upper):
            msg = f"Invalid allele {value!r}: expected IUPAC ACGTN, '*' or '.'"
            raise ValueError(msg)
        return upper

    @model_validator(mode="after")
    def _ref_alt_distinct(self) -> Self:
        if self.ref == self.alt:
            msg = f"REF and ALT are identical ({self.ref!r}) at {self.contig}:{self.position}"
            raise ValueError(msg)
        return self

    @property
    def variant_id(self) -> str:
        """Canonical, build-qualified identifier used as the join key everywhere."""
        return f"{self.build.value}:{self.contig}:{self.position}:{self.ref}:{self.alt}"

    @property
    def end(self) -> int:
        """Inclusive end of the REF span."""
        return self.position + len(self.ref) - 1

    @property
    def is_snv(self) -> bool:
        return len(self.ref) == 1 and len(self.alt) == 1 and self.alt not in {"*", "."}

    @property
    def is_indel(self) -> bool:
        return len(self.ref) != len(self.alt) and self.alt not in {"*", "."}

    def assert_same_build(self, other: GenomicCoordinate) -> None:
        """Guard every cross-variant comparison. Never coerce; always raise."""
        if self.build is not other.build:
            msg = (
                f"Genome build mismatch: {self.variant_id} ({self.build.value}) vs "
                f"{other.variant_id} ({other.build.value}). Cross-build comparison is "
                "refused; lift-over must be an explicit, provenance-tracked stage."
            )
            raise ValueError(msg)

    def sort_key(self) -> tuple[int, int, str, str]:
        return (contig_sort_key(self.contig), self.position, self.ref, self.alt)

    def __str__(self) -> str:
        return self.variant_id
