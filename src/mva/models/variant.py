"""Variant-level domain models."""

from __future__ import annotations

from enum import StrEnum
from typing import Self

from pydantic import Field, computed_field, model_validator

from mva.models.base import FrozenModel
from mva.models.genome import GenomeBuild, GenomicCoordinate


class Zygosity(StrEnum):
    """Called genotype state.

    ``UNKNOWN`` is a first-class value, not an error. A no-call is information
    ("we did not observe this") and must never be silently coerced to hom-ref,
    which would fabricate a negative finding.
    """

    HOM_REF = "hom_ref"
    HET = "het"
    HOM_ALT = "hom_alt"
    HEMIZYGOUS = "hemizygous"
    UNKNOWN = "unknown"


class FilterStatus(StrEnum):
    """Caller FILTER interpretation."""

    PASS = "pass"  # noqa: S105 - VCF FILTER value, not a credential
    FILTERED = "filtered"
    MISSING = "missing"


class ImpactSeverity(StrEnum):
    """Ensembl/VEP-style ordinal impact classes."""

    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    MODIFIER = "modifier"


class Genotype(FrozenModel):
    """The observed call plus the analytical evidence supporting it."""

    zygosity: Zygosity
    genotype_string: str = Field(description="Raw VCF GT, e.g. '0/1', '1|0', './.'")
    phased: bool = Field(default=False, description="True only if the caller emitted '|'.")
    phase_set: int | None = Field(default=None, description="VCF PS tag, when present.")
    depth: int | None = Field(default=None, ge=0, description="DP at this site.")
    ref_reads: int | None = Field(default=None, ge=0)
    alt_reads: int | None = Field(default=None, ge=0)
    genotype_quality: int | None = Field(default=None, ge=0, description="GQ.")

    @computed_field  # type: ignore[prop-decorator]
    @property
    def allele_balance(self) -> float | None:
        """ALT fraction. ``None`` when unknowable — never defaulted to 0.5.

        For a true germline het this sits near 0.5. Marked deviation suggests a
        sequencing artifact, allelic dropout, or — critically for this project —
        **somatic mosaicism**, which in Mosaic Variegated Aneuploidy is exactly the
        kind of signal that must not be flattened away.
        """
        if self.ref_reads is None or self.alt_reads is None:
            return None
        total = self.ref_reads + self.alt_reads
        if total == 0:
            return None
        return self.alt_reads / total

    @property
    def is_heterozygous(self) -> bool:
        return self.zygosity is Zygosity.HET

    @property
    def carries_alt(self) -> bool:
        return self.zygosity in {Zygosity.HET, Zygosity.HOM_ALT, Zygosity.HEMIZYGOUS}


class PopulationFrequency(FrozenModel):
    """An allele frequency observation.

    Every field here is mandatory-by-design except the counts: a bare float with no
    source, no version and no population is unusable. "AF = 0.001" is common in a
    Finnish cohort and vanishingly rare globally.
    """

    source: str = Field(description="e.g. 'gnomAD_genomes'")
    version: str = Field(description="e.g. 'v4.1.0'")
    population: str = Field(description="e.g. 'global', 'nfe', 'afr'")
    allele_frequency: float = Field(ge=0.0, le=1.0)
    allele_count: int | None = Field(default=None, ge=0)
    allele_number: int | None = Field(default=None, ge=0)
    homozygote_count: int | None = Field(default=None, ge=0)
    filter_status: str | None = Field(default=None, description="gnomAD site FILTER.")

    @property
    def provenance_key(self) -> str:
        return f"{self.source}/{self.version}/{self.population}"


class ConsequenceAnnotation(FrozenModel):
    """A predicted molecular consequence on one specific transcript.

    Deliberately transcript-scoped and stored as a *list* on the variant. Collapsing
    to the canonical transcript alone is a known way to lose clinically relevant
    effects (a variant can be benign on MANE-Select and splice-disrupting on the
    tissue-relevant isoform).
    """

    gene_symbol: str
    gene_id: str | None = None
    transcript_id: str
    transcript_biotype: str = "protein_coding"
    is_canonical: bool = False
    is_mane_select: bool = False
    consequence_terms: tuple[str, ...] = Field(
        description="Sequence Ontology terms, most severe first."
    )
    impact: ImpactSeverity
    hgvs_c: str | None = None
    hgvs_p: str | None = None
    exon: str | None = None
    intron: str | None = None
    protein_position: int | None = Field(default=None, ge=1)
    amino_acids: str | None = None
    splice_ai_delta_max: float | None = Field(default=None, ge=0.0, le=1.0)
    pathogenicity_scores: dict[str, float] = Field(
        default_factory=dict, description="e.g. {'CADD_phred': 28.1, 'REVEL': 0.81}"
    )
    source_tool: str = Field(default="unknown")
    source_tool_version: str = Field(default="unknown")

    @property
    def most_severe_term(self) -> str:
        return self.consequence_terms[0] if self.consequence_terms else "unknown"


class ClinicalAssertion(FrozenModel):
    """A curated clinical-significance assertion (ClinVar-style)."""

    source: str = "ClinVar"
    version: str
    accession: str | None = None
    significance: str = Field(description="e.g. 'Pathogenic', 'VUS', 'Benign'")
    review_status: str | None = None
    star_rating: int | None = Field(default=None, ge=0, le=4)
    conditions: tuple[str, ...] = ()


class VariantRecord(FrozenModel):
    """The canonical per-variant record flowing through the pipeline.

    Immutable and additive: annotation stages return a *new* record via
    ``with_annotations``, so the ingested record and the annotated record are
    separately hashable artifacts with distinct provenance.
    """

    coordinate: GenomicCoordinate
    genotype: Genotype
    filter_status: FilterStatus
    raw_filters: tuple[str, ...] = Field(default=(), description="Verbatim VCF FILTER field.")
    quality: float | None = Field(default=None, description="VCF QUAL.")

    consequences: tuple[ConsequenceAnnotation, ...] = ()
    population_frequencies: tuple[PopulationFrequency, ...] = ()
    clinical_assertions: tuple[ClinicalAssertion, ...] = ()

    qc_flags: tuple[str, ...] = Field(
        default=(),
        description="Non-destructive quality markers. Flagging down-ranks; it does not delete.",
    )
    source_artifact: str = Field(description="Reference to the artifact this came from.")
    source_line_index: int | None = Field(
        default=None, ge=0, description="Ordinal within the source, for traceability."
    )
    normalisation_ops: tuple[str, ...] = Field(
        default=(), description="e.g. ('split_multiallelic', 'left_align', 'trim')"
    )

    # ---------------------------------------------------------------- accessors

    @property
    def variant_id(self) -> str:
        return self.coordinate.variant_id

    @property
    def build(self) -> GenomeBuild:
        return self.coordinate.build

    @property
    def gene_symbols(self) -> tuple[str, ...]:
        """All genes this variant is annotated against, deduplicated, ordered."""
        seen: dict[str, None] = {}
        for csq in self.consequences:
            seen.setdefault(csq.gene_symbol, None)
        return tuple(seen)

    def consequences_for_gene(self, gene: str) -> tuple[ConsequenceAnnotation, ...]:
        return tuple(c for c in self.consequences if c.gene_symbol == gene)

    def worst_impact_for_gene(self, gene: str) -> ImpactSeverity | None:
        """Most severe predicted impact across *all* transcripts of a gene.

        Uses the maximum rather than the canonical transcript on purpose; see
        docs/scientific-assumptions.md (ASSUMPTION-TRANSCRIPT-01).
        """
        order = {
            ImpactSeverity.HIGH: 0,
            ImpactSeverity.MODERATE: 1,
            ImpactSeverity.LOW: 2,
            ImpactSeverity.MODIFIER: 3,
        }
        impacts = [c.impact for c in self.consequences_for_gene(gene)]
        return min(impacts, key=lambda i: order[i]) if impacts else None

    def max_allele_frequency(self) -> PopulationFrequency | None:
        """Highest observed AF across all recorded populations.

        Maximum (not global) is the conservative choice for rarity: a variant common
        in any single ancestry is not a plausible ultra-rare cause, and using the
        global AF alone systematically under-estimates frequency for variants
        enriched in populations under-represented in reference cohorts.
        """
        if not self.population_frequencies:
            return None
        return max(self.population_frequencies, key=lambda p: p.allele_frequency)

    @property
    def has_frequency_data(self) -> bool:
        """Absence of frequency data is NOT evidence of rarity. Callers must check."""
        return bool(self.population_frequencies)

    # ---------------------------------------------------------------- builders

    def with_annotations(
        self,
        *,
        consequences: tuple[ConsequenceAnnotation, ...] | None = None,
        population_frequencies: tuple[PopulationFrequency, ...] | None = None,
        clinical_assertions: tuple[ClinicalAssertion, ...] | None = None,
    ) -> VariantRecord:
        return self.model_copy(
            update={
                "consequences": consequences if consequences is not None else self.consequences,
                "population_frequencies": (
                    population_frequencies
                    if population_frequencies is not None
                    else self.population_frequencies
                ),
                "clinical_assertions": (
                    clinical_assertions
                    if clinical_assertions is not None
                    else self.clinical_assertions
                ),
            }
        )

    def with_qc_flags(self, *flags: str) -> VariantRecord:
        """Add QC flags. Additive and order-stable; flags are never removed."""
        merged = list(self.qc_flags)
        for flag in flags:
            if flag not in merged:
                merged.append(flag)
        return self.model_copy(update={"qc_flags": tuple(merged)})

    @model_validator(mode="after")
    def _phase_consistency(self) -> Self:
        if self.genotype.phased and "|" not in self.genotype.genotype_string:
            msg = (
                f"{self.variant_id}: genotype marked phased but GT "
                f"{self.genotype.genotype_string!r} contains no '|' separator."
            )
            raise ValueError(msg)
        return self

    def sort_key(self) -> tuple[int, int, str, str]:
        return self.coordinate.sort_key()
