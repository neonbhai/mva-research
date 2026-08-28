"""Candidate selection: the named, counted gate in front of pairing (ADR 0019).

This is the one place in the pipeline entitled to remove a variant that is
*valid* — and it exists because the alternative was worse, not because GP-13 was
relaxed.

The problem it solves
---------------------
``apply_hard_filters`` removes only invalid records, correctly (GP-13, ADR 0005).
So ``generate_pair_candidates`` receives every alt-carrying annotated variant
with a gene symbol: 2,020,500 of them for a whole-genome callset, mean 103.6 per
gene (``docs/scale-report.md`` §5). ``max_pairs_per_gene`` then truncates
**every one of the 19,500 genes**, discards 94.7% of hypotheses, and stamps
``gene_pair_cap_truncated`` on 100% of the output. A flag raised on everything
distinguishes nothing, and the ranking behind it is an artifact of a cap whose
ordering was never meant to be the primary filter.

The cap was doing a selection stage's job, silently, on two million variants.
This module does that job explicitly instead: by a stated rule, with every
dropped variant counted by reason, and with those counts reaching the run
artifacts and the evidence ledger. The cap then fires where it was designed to —
on the handful of hyper-variable genes (TTN, MUC16, OBSCN, NEB, RYR1, HLA) that
genuinely carry more rare coding variants than a hypothesis list can hold.

The three ways this stage could lose the answer, and what stops each
--------------------------------------------------------------------
**Absence of frequency data read as AF = 0 (GP-14).** A variant gnomAD has never
seen is *unknown*, not rare and not common. It is RETAINED by default
(:attr:`SelectionThresholds.retain_unknown_frequency`), because a genuinely novel
pathogenic variant is exactly what a rare-disease search is for, and because
absence correlates with poor reference coverage and under-represented ancestry —
so the "drop the unknown" policy loses the answer preferentially for the patients
least well served already. The same applies when every reporting population is
under-powered (``min_allele_number``): a frequency measured on 40 chromosomes
establishes neither commonness nor rarity (ADR 0010).

**``impact is None`` read as MODIFIER (ADR 0016).** ``None`` means NOT ASSESSED —
the source located the gene and computed no molecular consequence, which is what
a MANE interval join produces for every variant it places. MODIFIER is a positive
prediction of negligible effect. Collapsing the first into the second would drop
every gene-assignment-only variant as predicted-benign. A variant whose impact is
nowhere assessed is RETAINED (:attr:`SelectionThresholds.retain_unassessed_impact`).

**The ranking cut-point reused as a deletion cut-point.** ``FrequencyThresholds``
holds *ranking* cut-points: above ``max_plausible_recessive`` the rarity component
scores ~0 but the candidate survives, flagged. This stage's cut-point is separate
and deliberately looser, so that nothing this stage deletes was still rankable
above zero. A configuration that inverts that relationship is reported as a
warning rather than silently honoured.

Everything else about the stage is ordinary: it is deterministic, it holds one
record at a time, and it never mutates a record — a selected variant is the same
object that arrived, so no artifact downstream changes shape because selection
ran.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from datetime import datetime
from typing import Final

from mva.clock import Clock
from mva.config import FrequencyThresholds, SelectionThresholds
from mva.models.base import AssertionTier
from mva.models.evidence import (
    EvidenceCategory,
    EvidenceDirection,
    EvidenceItem,
    EvidenceStrength,
    EvidenceType,
    make_evidence_id,
)
from mva.models.variant import ConsequenceAnnotation, ImpactSeverity, VariantRecord

# ---------------------------------------------------------------------------
# Consequence vocabulary
# ---------------------------------------------------------------------------

#: Sequence Ontology terms that alter the coding sequence.
#:
#: ``synonymous_variant`` is deliberately ABSENT. It is coding but does not alter
#: the protein, and retaining it would roughly double the surviving set for no
#: gain. A synonymous call that also carries a splice term is still retained,
#: because matching is over a variant's whole term list across every transcript —
#: which is the case ``prioritization.filters`` warns about when it declines to
#: use predicted consequence at all.
CODING_CONSEQUENCE_TERMS: Final[frozenset[str]] = frozenset(
    {
        "coding_sequence_variant",
        "coding_transcript_variant",
        "conservative_inframe_deletion",
        "conservative_inframe_insertion",
        "disruptive_inframe_deletion",
        "disruptive_inframe_insertion",
        "exon_loss_variant",
        "feature_elongation",
        "feature_truncation",
        "frameshift_variant",
        "incomplete_terminal_codon_variant",
        "inframe_deletion",
        "inframe_insertion",
        "initiator_codon_variant",
        "missense_variant",
        "protein_altering_variant",
        "protein_protein_contact",
        "rare_amino_acid_variant",
        "start_lost",
        "start_retained_variant",
        "stop_gained",
        "stop_lost",
        "stop_retained_variant",
        "structural_interaction_variant",
        "transcript_ablation",
        "transcript_amplification",
    }
)

#: Terms that implicate splicing, including the weak ones.
#:
#: ``splice_region_variant`` is a deliberately generous inclusion: it fires on
#: synonymous and intronic calls near an exon boundary, which is precisely the
#: class a naive coding-only rule loses.
SPLICE_CONSEQUENCE_TERMS: Final[frozenset[str]] = frozenset(
    {
        "exonic_splice_region_variant",
        "splice_acceptor_variant",
        "splice_branch_variant",
        "splice_donor_5th_base_variant",
        "splice_donor_region_variant",
        "splice_donor_variant",
        "splice_polypyrimidine_tract_variant",
        "splice_region_variant",
        "splice_site_variant",
    }
)

#: Impacts that retain a variant whatever its terms say.
#:
#: Fail-open on purpose: an unrecognised SO term with a HIGH impact must not be
#: dropped because this module's vocabulary is out of date. A new annotation tool
#: adds terms far more often than it adds impact classes.
SELECTING_IMPACTS: Final[frozenset[ImpactSeverity]] = frozenset(
    {ImpactSeverity.HIGH, ImpactSeverity.MODERATE}
)

# ---------------------------------------------------------------------------
# Reason codes. Both sets partition their side of the decision exactly once, so
# the counts sum to the input and a reader can check the arithmetic.
# ---------------------------------------------------------------------------

KEEP_SELECTION_DISABLED = "retained_selection_disabled"
"""The stage was configured off. Everything through, nothing counted against it."""

KEEP_CLINICAL_ASSERTION = "retained_clinical_assertion"
"""A curated pathogenic assertion. Overrides frequency and consequence alike: a
known pathogenic allele that is common in some cohort is still the answer."""

KEEP_CODING_OR_SPLICE = "retained_coding_or_splice"
"""At least one transcript predicts a coding or splice-relevant consequence."""

KEEP_IMPACT_NOT_ASSESSED = "retained_impact_not_assessed"
"""A gene was located but no molecular consequence was computed on any transcript
(ADR 0016). NOT ASSESSED is not MODIFIER, and this variant has not been judged."""

KEEP_NO_GENE_ASSIGNMENT = "retained_no_gene_assignment"
"""No consequence annotation at all, retained because the policy says so. Pairing
is gene-scoped and will still not see it; the count says how many that is."""

DROP_COMMON_IN_POPULATION = "dropped_common_in_population"
"""An adequately powered population reports an allele frequency above the
selection cut-point. The only reason here that rests on a measured number."""

DROP_FREQUENCY_UNKNOWN = "dropped_frequency_unknown"
"""No usable frequency observation, and the policy was configured to drop those.
Off by default; turning it on is the GP-14 failure this stage is written around."""

DROP_NOT_CODING_OR_SPLICE = "dropped_not_coding_or_splice"
"""Every transcript was assessed and none is coding-altering or splice-relevant."""

DROP_IMPACT_NOT_ASSESSED = "dropped_impact_not_assessed"
"""Impact was never assessed and the policy was configured to drop those. Off by
default; ADR 0016 exists because this is how a not-assessed variant gets read as
predicted-harmless."""

DROP_NO_GENE_ASSIGNMENT = "dropped_no_gene_assignment"
"""No consequence annotation, so no gene, so no gene-scoped hypothesis. Counted
rather than silently ignored: pairing already discards these, invisibly."""

#: Fixed order, so report keys and their JSON serialisation are stable (GP-30).
KEEP_REASONS: Final[tuple[str, ...]] = (
    KEEP_SELECTION_DISABLED,
    KEEP_CLINICAL_ASSERTION,
    KEEP_CODING_OR_SPLICE,
    KEEP_IMPACT_NOT_ASSESSED,
    KEEP_NO_GENE_ASSIGNMENT,
)

DROP_REASONS: Final[tuple[str, ...]] = (
    DROP_COMMON_IN_POPULATION,
    DROP_FREQUENCY_UNKNOWN,
    DROP_NOT_CODING_OR_SPLICE,
    DROP_IMPACT_NOT_ASSESSED,
    DROP_NO_GENE_ASSIGNMENT,
)

NOTE_FREQUENCY_UNKNOWN = "frequency_unknown"
"""Retained with no usable frequency observation at all (GP-14)."""

NOTE_FREQUENCY_UNDERPOWERED = "frequency_underpowered"
"""Every reporting population fell below ``min_allele_number`` (ADR 0010)."""

NOTE_CLINICAL_PATHOGENIC = "clinical_pathogenic_assertion"
"""Carried a curated pathogenic assertion, whether or not that decided the case."""

NOTE_SPLICE_AI = "splice_ai_above_threshold"
"""Retained (at least in part) on a SpliceAI delta rather than on a term."""

#: Fixed order, as above. These annotate a decision; they do not partition it.
NOTES: Final[tuple[str, ...]] = (
    NOTE_FREQUENCY_UNKNOWN,
    NOTE_FREQUENCY_UNDERPOWERED,
    NOTE_CLINICAL_PATHOGENIC,
    NOTE_SPLICE_AI,
)

#: Subject of the aggregate evidence this stage emits. Constant rather than
#: run-derived so the evidence IDs are stable across runs of the same case (GP-30).
SELECTION_SUBJECT_ID: Final[str] = "candidate_selection"

_SELECTION_LIMITATION: Final[str] = (
    "A count of what this pipeline set aside, not a biological finding. Selection "
    "deletes valid candidates in order to make the per-gene pairing cap meaningful "
    "(ADR 0019), and a deleted candidate is invisible and unfalsifiable. It rests on "
    "predicted consequence — the weakest link in the annotation chain — and on "
    "reference-cohort frequencies that are population- and release-specific. Variants "
    "with no frequency data and variants with no assessed impact are retained rather "
    "than judged, so this stage's counts understate how much is genuinely unknown."
)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
#
# :class:`~mva.config.SelectionThresholds` is imported above rather than declared
# here, and is re-exported from this module's ``__all__`` so a reader who arrives
# at the stage finds the knob named where the stage is. It LIVES in config:
# thresholds belong in configuration, never as constants inside a filter (GP-32),
# and this is the one stage entitled to delete a valid record, so its cut-points
# have to be reviewable in the case file rather than in a dataclass default. ADR
# 0019 records why each default is what it is; ``config/default.yaml`` carries the
# values a run actually uses. Nothing in this module reads a cut-point that is not
# on that object or on the ``FrequencyThresholds`` passed alongside it.


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class SelectionDecision:
    """One verdict, with the reason and any qualifying notes."""

    variant: VariantRecord
    retained: bool
    reason: str
    notes: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class SelectionReport:
    """What the stage did, in counts. This is the provenance record.

    Counts, thresholds and warnings only — never a variant ID, so the artifact
    this becomes is classified DERIVED_SAFE alongside the QC report (GP-41).
    """

    input_count: int
    retained_count: int
    retained_by_reason: dict[str, int]
    dropped_by_reason: dict[str, int]
    notes: dict[str, int]
    thresholds: dict[str, bool | float]
    warnings: tuple[str, ...] = ()

    @property
    def dropped_count(self) -> int:
        return self.input_count - self.retained_count

    @property
    def retained_fraction(self) -> float:
        return (self.retained_count / self.input_count) if self.input_count else 0.0

    def as_payload(self) -> dict[str, object]:
        """The report artifact's body. Key order is fixed; values are counts."""
        return {
            "input_count": self.input_count,
            "retained_count": self.retained_count,
            "dropped_count": self.dropped_count,
            "retained_fraction": self.retained_fraction,
            "retained_by_reason": dict(self.retained_by_reason),
            "dropped_by_reason": dict(self.dropped_by_reason),
            "notes": dict(self.notes),
            "thresholds": dict(self.thresholds),
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class SelectionResult:
    """The materialised form: the selected variants, the report and its evidence."""

    variants: tuple[VariantRecord, ...]
    report: SelectionReport
    evidence: tuple[EvidenceItem, ...] = ()
    decisions: tuple[SelectionDecision, ...] = field(default=())
    """Every verdict, retained and dropped alike. Populated only by
    :func:`select_variants`, which is the fixture-scale API; the streaming API
    cannot hold 4.5 M of them and does not try."""


# ---------------------------------------------------------------------------
# The stage
# ---------------------------------------------------------------------------


def select_variants(
    variants: Sequence[VariantRecord],
    *,
    frequency: FrequencyThresholds,
    thresholds: SelectionThresholds | None = None,
    clock: Clock,
) -> SelectionResult:
    """Select plausible candidates, materialising every verdict.

    The fixture-scale API: it holds the input, the selected subset and one
    :class:`SelectionDecision` per input record. :func:`iter_selected` is the same
    computation, one record at a time, and is what a whole-genome caller uses.
    """
    stream = iter_selected(variants, frequency=frequency, thresholds=thresholds, clock=clock)
    decisions = tuple(stream.decisions())
    return SelectionResult(
        variants=tuple(decision.variant for decision in decisions if decision.retained),
        report=stream.report(),
        evidence=stream.evidence(),
        decisions=decisions,
    )


def iter_selected(
    variants: Iterable[VariantRecord],
    *,
    frequency: FrequencyThresholds,
    thresholds: SelectionThresholds | None = None,
    clock: Clock,
) -> SelectionStream:
    """Select plausible candidates from a stream, holding one record at a time.

    Iterating the stream yields the RETAINED records, in input order. Everything
    else is counted; :meth:`SelectionStream.report` and
    :meth:`SelectionStream.evidence` are meaningful once the stream is drained.
    A caller that wants every verdict rather than the survivors drives
    :meth:`SelectionStream.decisions` instead — but only one of the two, because
    the stream is single-pass.
    """
    return SelectionStream(
        variants,
        frequency=frequency,
        thresholds=thresholds or SelectionThresholds(),
        clock=clock,
    )


class SelectionStream:
    """Retained records one at a time; counts, warnings and evidence at the end.

    Construct via :func:`iter_selected`.
    """

    __slots__ = (
        "_dropped",
        "_exhausted",
        "_frequency",
        "_input",
        "_notes",
        "_retained",
        "_source",
        "_started",
        "_thresholds",
        "_timestamp",
    )

    def __init__(
        self,
        variants: Iterable[VariantRecord],
        *,
        frequency: FrequencyThresholds,
        thresholds: SelectionThresholds,
        clock: Clock,
    ) -> None:
        self._source = variants
        self._frequency = frequency
        self._thresholds = thresholds
        # Sampled once, before any record is seen (GP-30).
        self._timestamp = clock.now()
        self._input = 0
        self._retained: dict[str, int] = dict.fromkeys(KEEP_REASONS, 0)
        self._dropped: dict[str, int] = dict.fromkeys(DROP_REASONS, 0)
        self._notes: dict[str, int] = dict.fromkeys(NOTES, 0)
        self._started = False
        self._exhausted = False

    # ------------------------------------------------------------------ driving

    def __iter__(self) -> Iterator[VariantRecord]:
        return (decision.variant for decision in self.decisions() if decision.retained)

    def decisions(self) -> Iterator[SelectionDecision]:
        """Every verdict, retained and dropped alike, in input order."""
        if self._started:
            msg = (
                "This selection stream has already been iterated. Re-iterating would "
                "double every count. Call iter_selected() again against a re-read "
                "source."
            )
            raise ValueError(msg)
        self._started = True
        return self._drive()

    def _drive(self) -> Iterator[SelectionDecision]:
        for variant in self._source:
            self._input += 1
            decision = self._decide(variant)
            if decision.retained:
                self._retained[decision.reason] += 1
            else:
                self._dropped[decision.reason] += 1
            for note in decision.notes:
                self._notes[note] += 1
            yield decision
        self._exhausted = True

    # ------------------------------------------------------------------ the rule

    def _decide(self, variant: VariantRecord) -> SelectionDecision:
        thresholds = self._thresholds
        if not thresholds.enabled:
            return SelectionDecision(variant, retained=True, reason=KEEP_SELECTION_DISABLED)

        notes: list[str] = []

        pathogenic = thresholds.retain_pathogenic_clinical and _has_pathogenic_assertion(variant)
        if pathogenic:
            notes.append(NOTE_CLINICAL_PATHOGENIC)

        frequency_reason = self._frequency_verdict(variant, notes)
        # The clinical rescue runs after the frequency verdict rather than instead
        # of it, so a pathogenic variant that is also common still records
        # `frequency_unknown` / `frequency_underpowered` truthfully in the notes.
        if frequency_reason is not None and not pathogenic:
            return SelectionDecision(
                variant, retained=False, reason=frequency_reason, notes=tuple(notes)
            )

        consequence_reason, consequence_kept, consequence_notes = self._consequence_verdict(variant)
        notes.extend(consequence_notes)
        if pathogenic:
            return SelectionDecision(
                variant, retained=True, reason=KEEP_CLINICAL_ASSERTION, notes=tuple(notes)
            )
        return SelectionDecision(
            variant, retained=consequence_kept, reason=consequence_reason, notes=tuple(notes)
        )

    def _frequency_verdict(self, variant: VariantRecord, notes: list[str]) -> str | None:
        """``None`` when the variant passes the frequency gate, else a drop reason."""
        thresholds = self._thresholds
        selection = variant.select_max_allele_frequency(
            min_allele_number=self._frequency.min_allele_number
        )
        if selection.observed is None:
            # Two different gaps, both of which mean "unknown" and neither of which
            # means "rare" (GP-14). They are counted apart because they have
            # different fixes: one needs a better reference, one a bigger cohort.
            notes.append(
                NOTE_FREQUENCY_UNDERPOWERED if selection.has_exclusions else NOTE_FREQUENCY_UNKNOWN
            )
            return None if thresholds.retain_unknown_frequency else DROP_FREQUENCY_UNKNOWN
        if selection.observed.allele_frequency > thresholds.max_population_frequency:
            return DROP_COMMON_IN_POPULATION
        return None

    def _consequence_verdict(self, variant: VariantRecord) -> tuple[str, bool, tuple[str, ...]]:
        """``(reason, retained, notes)`` from the predicted consequences alone."""
        thresholds = self._thresholds
        consequences = variant.consequences
        if not consequences:
            if thresholds.drop_without_gene_assignment:
                return DROP_NO_GENE_ASSIGNMENT, False, ()
            return KEEP_NO_GENE_ASSIGNMENT, True, ()

        splice_ai = False
        for annotation in consequences:
            if _is_selecting(annotation):
                return KEEP_CODING_OR_SPLICE, True, ()
            delta = annotation.splice_ai_delta_max
            if delta is not None and delta >= thresholds.min_splice_ai_delta:
                splice_ai = True
        if splice_ai:
            return KEEP_CODING_OR_SPLICE, True, (NOTE_SPLICE_AI,)

        # Reached only when nothing argued for the variant. Whether that is a
        # judgement or an absence of one is decided by ADR 0016: an impact of None
        # on every transcript means nobody assessed it.
        if all(annotation.impact is None for annotation in consequences):
            if thresholds.retain_unassessed_impact:
                return KEEP_IMPACT_NOT_ASSESSED, True, ()
            return DROP_IMPACT_NOT_ASSESSED, False, ()
        return DROP_NOT_CODING_OR_SPLICE, False, ()

    # ------------------------------------------------------------- reporting

    @property
    def drained(self) -> bool:
        """Whether the stream ran to completion. Counts are partial until it did."""
        return self._exhausted

    def report(self) -> SelectionReport:
        """The counts, thresholds and warnings this run produced."""
        retained_total = sum(self._retained.values())
        return SelectionReport(
            input_count=self._input,
            retained_count=retained_total,
            retained_by_reason=dict(self._retained),
            dropped_by_reason=dict(self._dropped),
            notes=dict(self._notes),
            thresholds=self._thresholds.as_payload(),
            warnings=self._warnings(retained_total),
        )

    def _warnings(self, retained_total: int) -> tuple[str, ...]:
        thresholds = self._thresholds
        warnings: list[str] = []

        if self._started and not self._exhausted:
            warnings.append(
                f"Selection stream was not drained: {self._input} record(s) were read "
                "before the consumer stopped. Every count below describes that prefix."
            )

        if not thresholds.enabled:
            warnings.append(
                "Candidate selection is DISABLED: every variant was passed to pairing. "
                "The per-gene pairing cap then performs the selection instead, by an "
                "ordering that was never intended as the primary filter, and flags every "
                "candidate it emits (ADR 0019, docs/scale-report.md §5)."
            )
            return tuple(warnings)

        dropped = self._input - retained_total
        if dropped:
            breakdown = ", ".join(
                f"{reason}={self._dropped[reason]}"
                for reason in DROP_REASONS
                if self._dropped[reason]
            )
            warnings.append(
                f"Candidate selection removed {dropped}/{self._input} valid variants before "
                f"pairing ({breakdown}). This is the one stage entitled to delete a valid "
                "record, and a deleted candidate is invisible and unfalsifiable (GP-13, "
                "ADR 0019). Re-run with selection disabled to see what it cost."
            )

        unknown = self._notes[NOTE_FREQUENCY_UNKNOWN] + self._notes[NOTE_FREQUENCY_UNDERPOWERED]
        if unknown:
            verb = "retained" if thresholds.retain_unknown_frequency else "DROPPED"
            warnings.append(
                f"{unknown}/{self._input} variants reaching selection had no usable "
                f"population-frequency observation and were {verb}. Absence of frequency "
                "data is not evidence of rarity and must never be scored as allele "
                "frequency 0 (GP-14)."
            )
        if not thresholds.retain_unknown_frequency:
            warnings.append(
                "retain_unknown_frequency is FALSE: variants with no frequency data were "
                "deleted as if their absence were evidence about them. That is the GP-14 "
                "failure this stage was written to avoid, and it discards novel variants "
                "preferentially for ancestries the reference cohorts under-sample."
            )
        if not thresholds.retain_unassessed_impact:
            warnings.append(
                "retain_unassessed_impact is FALSE: variants whose molecular impact was "
                "never assessed were deleted as if they had been predicted harmless. "
                "`impact is None` means NOT ASSESSED, not MODIFIER (ADR 0016)."
            )

        ranking_cut = self._frequency.max_plausible_recessive
        if thresholds.max_population_frequency <= ranking_cut:
            warnings.append(
                f"Selection cut-point ({thresholds.max_population_frequency}) is not looser "
                f"than the ranking cut-point max_plausible_recessive ({ranking_cut}). "
                "Selection is therefore deleting variants the ranker could still have "
                "scored above zero, which collapses the filter/rank separation this "
                "pipeline is built on (GP-13, ADR 0005)."
            )
        return tuple(warnings)

    def evidence(self) -> tuple[EvidenceItem, ...]:
        """One aggregate item per non-zero reason, plus one for the stage itself.

        Aggregate rather than per-variant, and that is a real concession: GP-19
        wants failed candidates persisted with their reasons, and four million
        individual rejection items is not a thing this pipeline can hold. The
        counts are the persisted form, and ADR 0019 says so in as many words.
        """
        report = self.report()
        items = [
            _evidence_item(
                claim=(
                    f"Candidate selection retained {report.retained_count} of "
                    f"{report.input_count} annotated variants "
                    f"({report.retained_fraction:.6f}) for pairing."
                ),
                numeric_value=float(report.retained_count),
                payload={
                    "input_count": report.input_count,
                    "retained_count": report.retained_count,
                    "dropped_count": report.dropped_count,
                    **report.thresholds,
                },
                timestamp=self._timestamp,
            )
        ]
        items.extend(
            _evidence_item(
                claim=(f"Candidate selection recorded {count} variant(s) under {reason!r}."),
                numeric_value=float(count),
                payload={"reason": reason, "count": count, "input_count": report.input_count},
                timestamp=self._timestamp,
            )
            for reason, count in (
                *((r, report.retained_by_reason[r]) for r in KEEP_REASONS),
                *((r, report.dropped_by_reason[r]) for r in DROP_REASONS),
                *((n, report.notes[n]) for n in NOTES),
            )
            if count
        )
        return tuple(items)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _terms(annotation: ConsequenceAnnotation) -> Iterator[str]:
    """Individual SO terms, splitting the joined forms some tools emit.

    SnpEff writes ``splice_region_variant&synonymous_variant`` in a single ANN
    field, and a VEP-style adapter may pass ``+``-joined terms through. Matching
    the joined string against a set of atomic terms silently fails to recognise
    either, which for this stage means dropping a splice-region call as
    non-coding.
    """
    for term in annotation.consequence_terms:
        for token in term.replace("+", "&").split("&"):
            cleaned = token.strip().lower()
            if cleaned:
                yield cleaned


def _is_selecting(annotation: ConsequenceAnnotation) -> bool:
    """Does this one transcript argue for keeping the variant?"""
    if annotation.impact in SELECTING_IMPACTS:
        return True
    return any(
        term in CODING_CONSEQUENCE_TERMS or term in SPLICE_CONSEQUENCE_TERMS
        for term in _terms(annotation)
    )


def _has_pathogenic_assertion(variant: VariantRecord) -> bool:
    """A curated assertion of pathogenicity, read conservatively.

    ``Conflicting_classifications_of_pathogenicity`` contains the substring but is
    an explicit statement that curators disagree, so it is excluded — and it does
    not need to rescue anything, because a conflicting variant that is also rare
    and coding is retained on its own merits anyway.
    """
    for assertion in variant.clinical_assertions:
        token = assertion.significance.strip().lower()
        if "pathogenic" in token and "conflict" not in token and "benign" not in token:
            return True
    return False


def _evidence_item(
    *,
    claim: str,
    numeric_value: float,
    payload: Mapping[str, str | int | float | bool | None],
    timestamp: datetime,
) -> EvidenceItem:
    return EvidenceItem(
        evidence_id=make_evidence_id(
            subject_id=SELECTION_SUBJECT_ID,
            category=EvidenceCategory.PROVENANCE,
            claim=claim,
            tool="mva-selection",
        ),
        subject_id=SELECTION_SUBJECT_ID,
        subject_kind="run",
        claim=claim,
        category=EvidenceCategory.PROVENANCE,
        # Bookkeeping about what the pipeline did, not an argument about a
        # hypothesis. NEUTRAL is the only honest direction for it.
        direction=EvidenceDirection.NEUTRAL,
        strength=EvidenceStrength.MODERATE,
        evidence_type=EvidenceType.PIPELINE_INFERENCE,
        tier=AssertionTier.INFERENCE,
        citation=None,
        method=(
            "Counted by mva.prioritization.selection over the annotated variant stream, "
            "against the thresholds recorded in the payload."
        ),
        tool="mva-selection",
        tool_version="mva-selection/0.1.0",
        limitations=_SELECTION_LIMITATION,
        timestamp=timestamp,
        numeric_value=numeric_value,
        payload=dict(payload),
    )


__all__ = [
    "CODING_CONSEQUENCE_TERMS",
    "DROP_COMMON_IN_POPULATION",
    "DROP_FREQUENCY_UNKNOWN",
    "DROP_IMPACT_NOT_ASSESSED",
    "DROP_NOT_CODING_OR_SPLICE",
    "DROP_NO_GENE_ASSIGNMENT",
    "DROP_REASONS",
    "KEEP_CLINICAL_ASSERTION",
    "KEEP_CODING_OR_SPLICE",
    "KEEP_IMPACT_NOT_ASSESSED",
    "KEEP_NO_GENE_ASSIGNMENT",
    "KEEP_REASONS",
    "KEEP_SELECTION_DISABLED",
    "NOTES",
    "NOTE_CLINICAL_PATHOGENIC",
    "NOTE_FREQUENCY_UNDERPOWERED",
    "NOTE_FREQUENCY_UNKNOWN",
    "NOTE_SPLICE_AI",
    "SELECTING_IMPACTS",
    "SELECTION_SUBJECT_ID",
    "SPLICE_CONSEQUENCE_TERMS",
    "SelectionDecision",
    "SelectionReport",
    "SelectionResult",
    "SelectionStream",
    "SelectionThresholds",
    "iter_selected",
    "select_variants",
]
