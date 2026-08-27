"""Gene-level pair construction and phase inference (GP-15).

Two heterozygous variants in one gene are a compound heterozygote **only if
they sit on opposite haplotypes**. Proband-only short-read data usually cannot
tell, and the honest answer — ``PhaseStatus.UNKNOWN`` — is preserved here and
carried all the way into the report. Nothing in this module ever upgrades an
unknown phase to trans; the only upgrades are downgrades to *cis*, which is
near-disqualifying and therefore the one direction where being wrong is safe.

Scope note: candidates are gene-scoped, so a variant with no gene annotation
forms no candidate. That is a limitation of a gene-based hypothesis space, not a
filter — such records remain in the flagged variant set the caller holds.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from itertools import combinations

from mva.models.pair import InheritanceModel, PhaseEvidence, PhaseStatus, make_pair_id
from mva.models.variant import ImpactSeverity, VariantRecord, Zygosity

#: Approximate span a short-read fragment can bridge. Beyond this, read-backed
#: phasing is physically impossible, which is *why* phase is unknown — it says
#: nothing whatsoever about whether the variants are in trans.
READ_BACKED_PHASING_SPAN_BP = 500

PHASE_METHOD_PHASE_SET = "phase_set"
PHASE_METHOD_NONE = "none"

FLAG_PHASE_UNKNOWN = "phase_unknown"
FLAG_PHASE_CIS = "phase_in_cis"
FLAG_PHASE_TRANS = "phase_in_trans"
FLAG_SINGLE_VARIANT = "single_variant_hypothesis"
FLAG_MIXED_ZYGOSITY = "mixed_zygosity_pair"
FLAG_CAP_TRUNCATED = "gene_pair_cap_truncated"

#: Member-variant soft flags promoted onto the candidate so a reviewer sees the
#: caveat on the ranked row without opening the variant records.
PROMOTED_VARIANT_FLAGS: tuple[str, ...] = (
    "common_variant",
    "low_frequency_variant",
    "no_frequency_data",
    "benign_consequence",
    "low_quality_call",
    "possible_mosaic",
    "homozygous_call",
)

#: Sorts before any real coordinate, so single-variant candidates order ahead of
#: two-variant candidates that share the same first variant.
_NO_SECOND_VARIANT: tuple[int, int, str, str] = (-1, -1, "", "")


@dataclass(frozen=True)
class PairCandidate:
    """A gene-scoped hypothesis: one or two variants proposed as jointly causal."""

    pair_id: str
    gene_symbol: str
    variant_a: VariantRecord
    variant_b: VariantRecord | None
    inheritance_model: InheritanceModel
    phase: PhaseEvidence
    flags: tuple[str, ...]

    @property
    def variants(self) -> tuple[VariantRecord, ...]:
        return (self.variant_a,) if self.variant_b is None else (self.variant_a, self.variant_b)

    @property
    def variant_ids(self) -> tuple[str, ...]:
        return tuple(variant.variant_id for variant in self.variants)

    @property
    def is_pair(self) -> bool:
        return self.variant_b is not None

    def has_variant_flag(self, flag: str) -> bool:
        return any(flag in variant.qc_flags for variant in self.variants)

    def sort_key(self) -> tuple[str, tuple[int, int, str, str], tuple[int, int, str, str], str]:
        """Total order: gene, first coordinate, second coordinate, pair id."""
        second = _NO_SECOND_VARIANT if self.variant_b is None else self.variant_b.sort_key()
        return (self.gene_symbol, self.variant_a.sort_key(), second, self.pair_id)


def _alt_haplotype_indices(genotype_string: str) -> frozenset[int] | None:
    """Haplotype slots carrying an alternate allele in a phased diploid GT.

    ``'1|0'`` -> ``{0}``, ``'0|1'`` -> ``{1}``, ``'1|1'`` -> ``{0, 1}``.
    Returns ``None`` when the genotype is unphased, not diploid, or contains a
    missing allele — all cases where nothing may be concluded.
    """
    if "|" not in genotype_string:
        return None
    alleles = genotype_string.split("|")
    if len(alleles) != 2:
        return None
    indices: set[int] = set()
    for slot, allele in enumerate(alleles):
        if allele in {".", ""}:
            return None
        if allele != "0":
            indices.add(slot)
    return frozenset(indices)


def _distance_bp(a: VariantRecord, b: VariantRecord) -> int | None:
    if a.coordinate.contig != b.coordinate.contig:
        return None
    return abs(a.coordinate.position - b.coordinate.position)


def _unknown_phase_note(a: VariantRecord, b: VariantRecord, distance: int | None) -> str:
    """State *why* phase could not be determined. Never why it might be trans."""
    if distance is None:
        return (
            "Variants lie on different contigs, so no read or phase-set evidence can "
            "relate their haplotypes. Phase is unknown, which is not evidence of trans."
        )
    same_set = (
        a.genotype.phase_set is not None
        and b.genotype.phase_set is not None
        and a.genotype.phase_set == b.genotype.phase_set
    )
    if not (a.genotype.phased and b.genotype.phased):
        reason = "at least one genotype was emitted unphased by the caller"
    elif not same_set:
        reason = "the two phased genotypes belong to different phase sets"
    else:
        reason = "the phased genotypes do not resolve to distinct single haplotypes"
    span = (
        f"the sites are {distance} bp apart, beyond the ~{READ_BACKED_PHASING_SPAN_BP} bp "
        "a short-read fragment can bridge, so read-backed phasing cannot reach them"
        if distance > READ_BACKED_PHASING_SPAN_BP
        else (
            f"the sites are {distance} bp apart, within read-backed range, so targeted "
            "re-analysis of the alignments may yet resolve them"
        )
    )
    return (
        f"Phase undetermined because {reason}; {span}. UNKNOWN is preserved rather than "
        "assumed to be trans (GP-15)."
    )


def infer_phase(a: VariantRecord, b: VariantRecord) -> PhaseEvidence:
    """Determine haplotype relationship from caller phase sets only.

    Resolves to ``CIS_CONFIRMED`` when both genotypes are phased into the same
    phase set on the same haplotype slot, and ``TRANS_CONFIRMED`` when they are
    phased into the same set on opposite slots. Every other situation is
    ``UNKNOWN`` with method ``none``, a recorded inter-site distance and a note
    saying what was missing. There is no path in this function that upgrades an
    absent observation into trans.
    """
    a.coordinate.assert_same_build(b.coordinate)
    distance = _distance_bp(a, b)

    same_phase_set = (
        a.genotype.phased
        and b.genotype.phased
        and a.genotype.phase_set is not None
        and a.genotype.phase_set == b.genotype.phase_set
    )
    if same_phase_set:
        slots_a = _alt_haplotype_indices(a.genotype.genotype_string)
        slots_b = _alt_haplotype_indices(b.genotype.genotype_string)
        if slots_a is not None and slots_b is not None and len(slots_a) == 1 == len(slots_b):
            if slots_a == slots_b:
                return PhaseEvidence(
                    status=PhaseStatus.CIS_CONFIRMED,
                    method=PHASE_METHOD_PHASE_SET,
                    distance_bp=distance,
                    notes=(
                        f"Both alternate alleles sit on haplotype {next(iter(slots_a))} of "
                        f"phase set {a.genotype.phase_set}. One gene copy is therefore "
                        "intact, which a compound heterozygote cannot be."
                    ),
                )
            return PhaseEvidence(
                status=PhaseStatus.TRANS_CONFIRMED,
                method=PHASE_METHOD_PHASE_SET,
                distance_bp=distance,
                notes=(
                    f"Alternate alleles occupy opposite haplotypes of phase set "
                    f"{a.genotype.phase_set}; both gene copies are affected."
                ),
            )

    return PhaseEvidence(
        status=PhaseStatus.UNKNOWN,
        method=PHASE_METHOD_NONE,
        distance_bp=distance,
        notes=_unknown_phase_note(a, b, distance),
    )


def _single_variant_phase(variant: VariantRecord) -> PhaseEvidence:
    return PhaseEvidence(
        status=PhaseStatus.UNKNOWN,
        method=PHASE_METHOD_NONE,
        distance_bp=None,
        notes=(
            f"Single-variant hypothesis for {variant.variant_id}: there is no second "
            "allele to phase against. Phase remains UNKNOWN rather than being recorded "
            "as resolved."
        ),
    )


def _phase_flags(phase: PhaseEvidence) -> tuple[str, ...]:
    if phase.status is PhaseStatus.UNKNOWN:
        return (FLAG_PHASE_UNKNOWN,)
    if phase.status in {PhaseStatus.CIS_CONFIRMED, PhaseStatus.CIS_LIKELY}:
        return (FLAG_PHASE_CIS,)
    return (FLAG_PHASE_TRANS,)


def _promoted_flags(variants: tuple[VariantRecord, ...]) -> tuple[str, ...]:
    present = {flag for variant in variants for flag in variant.qc_flags}
    return tuple(flag for flag in PROMOTED_VARIANT_FLAGS if flag in present)


def _pair_inheritance_model(a: VariantRecord, b: VariantRecord) -> tuple[InheritanceModel, bool]:
    """Model for a two-variant candidate, plus whether the zygosities are mixed."""
    if a.genotype.zygosity is Zygosity.HET and b.genotype.zygosity is Zygosity.HET:
        return InheritanceModel.COMPOUND_HETEROZYGOUS, False
    # A homozygous or hemizygous member already accounts for both gene copies on
    # its own, so the two calls are not a compound heterozygote. The combined
    # hypothesis stays representable, but the model is honestly UNKNOWN.
    return InheritanceModel.UNKNOWN, True


def _single_inheritance_model(variant: VariantRecord) -> InheritanceModel:
    if variant.genotype.zygosity is Zygosity.HOM_ALT:
        return InheritanceModel.HOMOZYGOUS_RECESSIVE
    if variant.genotype.zygosity is Zygosity.HEMIZYGOUS:
        return (
            InheritanceModel.X_LINKED_RECESSIVE
            if variant.coordinate.contig == "chrX"
            else InheritanceModel.HOMOZYGOUS_RECESSIVE
        )
    return InheritanceModel.UNKNOWN


def _build_pair_candidate(gene: str, a: VariantRecord, b: VariantRecord) -> PairCandidate:
    model, mixed = _pair_inheritance_model(a, b)
    phase = infer_phase(a, b)
    flags = list(_phase_flags(phase))
    if mixed:
        flags.append(FLAG_MIXED_ZYGOSITY)
    flags.extend(_promoted_flags((a, b)))
    return PairCandidate(
        pair_id=make_pair_id(gene, (a.variant_id, b.variant_id)),
        gene_symbol=gene,
        variant_a=a,
        variant_b=b,
        inheritance_model=model,
        phase=phase,
        flags=tuple(flags),
    )


def _build_single_candidate(gene: str, variant: VariantRecord) -> PairCandidate:
    phase = _single_variant_phase(variant)
    flags = [FLAG_SINGLE_VARIANT, *_phase_flags(phase), *_promoted_flags((variant,))]
    return PairCandidate(
        pair_id=make_pair_id(gene, (variant.variant_id,)),
        gene_symbol=gene,
        variant_a=variant,
        variant_b=None,
        inheritance_model=_single_inheritance_model(variant),
        phase=phase,
        flags=tuple(flags),
    )


def _wants_single_candidate(gene: str, variant: VariantRecord) -> bool:
    """Single-variant hypotheses worth carrying alongside the pairs.

    A homozygous (or hemizygous) call explains both gene copies by itself. A
    lone heterozygote is kept only when its predicted impact is HIGH, so that a
    dominant or de-novo model stays representable instead of being lost to a
    recessive-shaped pipeline.
    """
    zygosity = variant.genotype.zygosity
    if zygosity in {Zygosity.HOM_ALT, Zygosity.HEMIZYGOUS}:
        return True
    return zygosity is Zygosity.HET and variant.worst_impact_for_gene(gene) is ImpactSeverity.HIGH


def _variants_by_gene(
    variants: Sequence[VariantRecord],
) -> dict[str, list[VariantRecord]]:
    """Group alt-carrying variants by gene; a variant may appear under several."""
    grouped: dict[str, list[VariantRecord]] = {}
    for variant in variants:
        if not variant.genotype.carries_alt:
            continue
        for gene in variant.gene_symbols:
            grouped.setdefault(gene, []).append(variant)
    for members in grouped.values():
        members.sort(key=lambda v: (v.sort_key(), v.variant_id))
    return grouped


def generate_pairs(
    variants: Sequence[VariantRecord], *, max_pairs_per_gene: int = 20
) -> tuple[PairCandidate, ...]:
    """Enumerate every gene-scoped candidate hypothesis, deterministically.

    Within each gene: all unordered pairs of alt-carrying variants, plus a
    single-variant candidate for each homozygous/hemizygous call and each
    HIGH-impact heterozygote. Candidates are ordered by genomic position and
    truncated at ``max_pairs_per_gene``; when truncation happens every surviving
    candidate for that gene carries ``gene_pair_cap_truncated``, because a
    silently shortened hypothesis list is indistinguishable from a complete one.
    """
    grouped = _variants_by_gene(variants)
    candidates: list[PairCandidate] = []

    for gene in sorted(grouped):
        members = grouped[gene]
        gene_candidates: list[PairCandidate] = [
            _build_single_candidate(gene, variant)
            for variant in members
            if _wants_single_candidate(gene, variant)
        ]
        gene_candidates.extend(
            _build_pair_candidate(gene, a, b) for a, b in combinations(members, 2)
        )
        gene_candidates.sort(key=lambda candidate: candidate.sort_key())

        if len(gene_candidates) > max_pairs_per_gene:
            gene_candidates = [
                PairCandidate(
                    pair_id=candidate.pair_id,
                    gene_symbol=candidate.gene_symbol,
                    variant_a=candidate.variant_a,
                    variant_b=candidate.variant_b,
                    inheritance_model=candidate.inheritance_model,
                    phase=candidate.phase,
                    flags=(*candidate.flags, FLAG_CAP_TRUNCATED),
                )
                for candidate in gene_candidates[:max_pairs_per_gene]
            ]
        candidates.extend(gene_candidates)

    return tuple(candidates)
