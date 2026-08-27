"""Variant-level domain models."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Final, Self

from pydantic import Field, computed_field, model_validator

from mva.models.base import FrozenModel
from mva.models.genome import GenomeBuild, GenomicCoordinate

# ---------------------------------------------------------------------------
# Shared vocabulary
#
# These names are the contract between `mva.ingestion`, which emits them, and
# `mva.prioritization`, which reacts to them. They live here because GP-01/GP-03
# forbid either package importing the other, and the alternative — re-typing the
# string literals in both places — is how `INGESTION_QUALITY_FLAGS` came to name
# eight flags nothing emitted while missing three that ingestion really raises.
# One definition, imported by both sides, is the only version a lint can check.
# ---------------------------------------------------------------------------

#: Normalisation operations recorded in ``VariantRecord.normalisation_ops``.
OP_SPLIT_MULTIALLELIC: Final = "split_multiallelic"
OP_TRIM: Final = "trim"
OP_LEFT_ALIGN: Final = "left_align"

FLAG_REF_ALLELE_MISMATCH: Final = "ref_allele_mismatch"
FLAG_FILTERED_BY_CALLER: Final = "filtered_by_caller"
FLAG_NO_CALLER_FILTER: Final = "no_caller_filter"
FLAG_LOW_DEPTH: Final = "low_depth"
FLAG_LOW_GQ: Final = "low_gq"
FLAG_POSSIBLE_MOSAIC: Final = "possible_mosaic"
FLAG_LOW_ALLELE_BALANCE: Final = "low_allele_balance"
FLAG_HIGH_ALLELE_BALANCE: Final = "high_allele_balance"
FLAG_NO_QUALITY_METRICS: Final = "no_quality_metrics"

#: Every QC flag the ingestion stage is entitled to attach to a record. A
#: downstream stage that reacts to ingestion findings by name checks itself
#: against this set, so a flag added on one side cannot go unrecognised on the
#: other.
INGESTION_QC_FLAGS: Final[frozenset[str]] = frozenset(
    {
        FLAG_REF_ALLELE_MISMATCH,
        FLAG_FILTERED_BY_CALLER,
        FLAG_NO_CALLER_FILTER,
        FLAG_LOW_DEPTH,
        FLAG_LOW_GQ,
        FLAG_POSSIBLE_MOSAIC,
        FLAG_LOW_ALLELE_BALANCE,
        FLAG_HIGH_ALLELE_BALANCE,
        FLAG_NO_QUALITY_METRICS,
    }
)

#: The subset of :data:`INGESTION_QC_FLAGS` that argues the genotype call itself
#: may not be real. Membership means "do not trust this call", which is a
#: stronger claim than "something was noted here".
#:
#: Three members of :data:`INGESTION_QC_FLAGS` are deliberately absent:
#:
#: * ``possible_mosaic`` — a skewed het may be the finding itself in a mosaic
#:   aneuploidy disorder (ASSUMPTION-MOSAIC-01). Calling it low quality buries
#:   exactly the signal this pipeline exists to surface.
#: * ``no_quality_metrics`` — no DP, GQ or AD was reported. That is absence of
#:   information, not evidence of poor quality (GP-14).
#: * ``no_caller_filter`` — the caller recorded no FILTER opinion at all. It
#:   never filtered anything, so it never objected to anything.
UNTRUSTED_CALL_FLAGS: Final[frozenset[str]] = frozenset(
    {
        FLAG_REF_ALLELE_MISMATCH,
        FLAG_FILTERED_BY_CALLER,
        FLAG_LOW_DEPTH,
        FLAG_LOW_GQ,
        FLAG_LOW_ALLELE_BALANCE,
        FLAG_HIGH_ALLELE_BALANCE,
    }
)


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
    alt_allele_index: int | None = Field(
        default=None,
        ge=1,
        description=(
            "1-based index of THIS record's ALT within the source site's ALT list. "
            "Set when the record was decomposed from a multiallelic site, where the "
            "raw GT is retained verbatim and so no longer says which of its allele "
            "numbers belongs to this record. Without it a phased '1|2' reads as "
            "'both haplotypes carry an alternate allele' and resolvable phase is "
            "thrown away. ``None`` means the index is unrecorded, and every consumer "
            "must fall back to its pre-existing any-alternate-allele behaviour."
        ),
    )

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


@dataclass(frozen=True, slots=True)
class FrequencySelection:
    """Which population frequency was used for rarity, and which were set aside.

    ``excluded`` exists so that skipping an under-powered cohort is a statement a
    report can make, rather than a silent omission. A reader who disagrees with
    the threshold can see exactly what it cost them.
    """

    observed: PopulationFrequency | None
    excluded: tuple[PopulationFrequency, ...]
    min_allele_number: int

    @property
    def has_exclusions(self) -> bool:
        return bool(self.excluded)

    def describe_exclusions(self) -> str:
        """One clause naming every skipped cohort, for a rationale sentence."""
        if not self.excluded:
            return ""
        rendered = "; ".join(
            f"{frequency.provenance_key} at AF {frequency.allele_frequency:.3g} "
            f"(AC {frequency.allele_count}/AN {frequency.allele_number})"
            for frequency in self.excluded
        )
        return (
            f"excluded from the maximum for reporting fewer than {self.min_allele_number} "
            f"alleles: {rendered}"
        )


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

    @property
    def allele_fraction(self) -> float | None:
        """The ALT read fraction that the heterozygous band may actually be applied to.

        For a biallelic call this is exactly ``Genotype.allele_balance``.

        For a record decomposed from a multiallelic site it is not, and the
        difference decides whether a textbook compound heterozygote survives. After
        splitting ``A>G,GT`` with ``GT=1/2`` and ``AD=2,21,21``, allele 1 holds
        ``ref_reads=2`` and ``alt_reads=21`` — but those two reads are the *site's*
        REF depth and exclude the 21 reads carrying the other ALT. ``alt/(ref+alt)``
        is then 0.91 and looks homozygous, so the call gets flagged as an artifact
        of its own VCF formatting. Where the site depth is known, the site allele
        fraction ``alt/DP`` (0.48 here) is the right quantity; where it is not,
        ``None`` is returned and no caller may raise an allele-balance finding,
        because a confident wrong flag is worse than no flag (GP-14).

        Lives on the record rather than on ``Genotype`` because deciding which
        formula applies needs ``normalisation_ops``, which the genotype cannot see.
        """
        genotype = self.genotype
        balance = genotype.allele_balance
        if OP_SPLIT_MULTIALLELIC not in self.normalisation_ops:
            return balance
        depth, alt_reads, ref_reads = genotype.depth, genotype.alt_reads, genotype.ref_reads
        if depth is None or alt_reads is None:
            return None
        if depth <= 0 or depth < (ref_reads or 0) + alt_reads:
            return balance
        return alt_reads / depth

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

    def select_max_allele_frequency(self, *, min_allele_number: int = 0) -> FrequencySelection:
        """Pick the frequency row rarity should be judged on, and say what was skipped.

        Maximum (not global) AF is the conservative choice: a variant common in any
        single ancestry is not a plausible ultra-rare cause, and the global figure
        systematically under-estimates frequency for alleles enriched in cohorts
        that reference panels under-sample (ASSUMPTION-FREQUENCY-02).

        Taken raw, though, "maximum" hands the decision to whichever population was
        sampled least. One allele seen once in a 40-chromosome subpopulation is an
        AF of 0.025 — above any plausible recessive cut-point — while the same site
        sits at 8e-6 across 125,000 global chromosomes. The maximum then reports a
        genuine founder allele as common on the strength of a single observation.
        gnomAD's grpmax and the ACMG BA1/BS1 rules avoid this with a filtering AF
        (the lower bound of the 95% CI); requiring a minimum ``allele_number``
        before a population may set the maximum is the same guard in its simplest
        form, and it uses a field this model already stores.

        Populations below ``min_allele_number`` are RECORDED in
        :attr:`FrequencySelection.excluded`, never silently dropped: the caller can
        state which cohort was set aside and why. A population with no ``allele_number``
        at all stays eligible — an unreported cohort size is unknown, not small, and
        excluding it would discard the only record a source supplied (GP-14).

        When every population falls below the threshold, no adequately powered
        observation exists and ``observed`` is ``None``. Downstream that reads as
        absence of frequency data, which is scored mid-range rather than as rarity
        — the honest answer, because a frequency measured on 40 chromosomes
        establishes neither commonness nor rarity.

        See docs/decisions/0010-filtering-allele-frequency.md.
        """
        eligible: list[PopulationFrequency] = []
        excluded: list[PopulationFrequency] = []
        for frequency in self.population_frequencies:
            powered = (
                frequency.allele_number is None or frequency.allele_number >= min_allele_number
            )
            (eligible if powered else excluded).append(frequency)
        observed = (
            max(eligible, key=lambda p: (p.allele_frequency, p.provenance_key))
            if eligible
            else None
        )
        return FrequencySelection(
            observed=observed,
            excluded=tuple(sorted(excluded, key=lambda p: (-p.allele_frequency, p.provenance_key))),
            min_allele_number=min_allele_number,
        )

    def max_allele_frequency(self, *, min_allele_number: int = 0) -> PopulationFrequency | None:
        """Highest AF across populations whose cohort is large enough to be believed.

        Thin accessor over :meth:`select_max_allele_frequency`; use that when the
        excluded populations need to be reported. The default of ``0`` applies no
        guard, so a caller with no configured threshold keeps the historical
        behaviour rather than silently acquiring a new one.
        """
        return self.select_max_allele_frequency(min_allele_number=min_allele_number).observed

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
