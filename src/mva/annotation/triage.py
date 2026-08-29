"""The relevance bar: which annotated variants are worth the expensive half.

ADR 0028. Consequence annotation runs over the whole callset — it is one JVM
launch and a streaming pass, and it is cheap. What is not cheap is what the
pipeline used to do next for every one of those records: a gnomAD region query, a
ClinVar region query, three ``EvidenceItem``s and ~6 KB of JSON. A measured
whole-genome run was killed after 7 h 30 m having annotated 38.5% of 4,962,060
variants, projecting ~19.5 h and ~95 GB of intermediates onto a 60 GB volume.

Of that callset, 14,228 variants — 0.284% — clear a HIGH/MODERATE impact bar.
The pipeline was computing the selectivity and then throwing it away.

This module is the point where it stops. It reads a variant's consequence
annotations and nothing else, and returns a verdict. That constraint is what makes
the partition possible at all: the verdict has to be reachable **before** the
frequency and clinical adapters are called, because avoiding those calls is the
entire saving. It follows that a rule needing an allele frequency or a ClinVar
assertion cannot live here — see :ref:`the clinical rescue <clinical-rescue>`
below, which is this design's named residual risk.

The bar
-------
A variant is kept when any transcript annotation says any of:

* the impact is at least ``min_impact`` (default LOW — so HIGH, MODERATE and LOW
  all survive, and only a MODIFIER-only variant is at risk);
* any Sequence Ontology term is coding-altering or splice-relevant, whatever the
  impact says — the shared
  :func:`~mva.models.variant.consequence_is_selecting` predicate, which is the
  same one the selection stage applies later;
* SpliceAI reports a delta at or above ``min_splice_ai_delta``;
* the impact was never assessed on any transcript (ADR 0016) — NOT ASSESSED is
  not a prediction of harmlessness;
* the record already carries a clinical assertion from its input.

A variant with no consequence annotation at all is dropped by default, and that
one is provable rather than heuristic: pairing is gene-scoped
(``_variants_by_gene`` iterates ``variant.gene_symbols``), so a variant with no
gene symbol forms no candidate pair and cannot reach the output under any
configuration.

Why the bar is wider than the numbers require
---------------------------------------------
At 14,228 survivors on HIGH/MODERATE there is room to be wrong. The default keeps
the whole LOW impact band on top of that, because impact class is an ordinal
bucket produced by a predictor, and the class this pipeline can least afford to
lose — a synonymous variant that disrupts splicing — is precisely the one a
tight impact bar loses. ``prioritization.filters`` declines to use predicted
consequence at all for the same reason, in as many words.

.. _clinical-rescue:

The residual risk, stated
-------------------------
``SelectionThresholds.retain_pathogenic_clinical`` lets a curated pathogenic
assertion rescue a variant that both the frequency and the consequence gate would
otherwise drop. That rescue needs a ClinVar lookup, and a ClinVar lookup is one of
the two things this bar exists to avoid. So a variant that is MODIFIER on every
transcript, carries no coding or splice term, has no SpliceAI signal, and holds a
pathogenic ClinVar assertion is dropped here where selection would have kept it.

Deep-intronic pathogenic alleles are real and SnpEff calls them
``intron_variant``/MODIFIER, so this class is not empty. Three things bound it,
and none of them makes it disappear:

1. Every such variant is in the discard ledger by id, with ``triage`` as the
   stage and its impact as the value. A targeted ClinVar lookup over the dropped
   ids is a scan of one file, not a re-run.
2. ``triage.enabled: false`` restores the old behaviour exactly, and the run
   warning says so when the risk is live.
3. ``triage.min_impact: modifier`` keeps every annotated record while still
   dropping the ~55% with no gene assignment, which is the halfway position for a
   run that wants the ClinVar rescue and most of the saving.

Nothing here is a scoring decision. The verdict decides how much *work* a variant
gets, never how it ranks, and the composition root proves the point by asserting
that the submission is byte-identical with the bar on and off.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Final

from mva.config import SelectionThresholds, TriageThresholds
from mva.models.discard import DiscardRecord, DiscardStage
from mva.models.variant import (
    ClinicalAssertion,
    ConsequenceAnnotation,
    ImpactSeverity,
    consequence_is_selecting,
    impact_at_least,
    impact_rank,
)

# ---------------------------------------------------------------------------
# Reason codes. Each set partitions its side of the decision exactly once, so the
# counts sum to the input and a reader can check the arithmetic — the same
# property ADR 0019 gave the selection report.
# ---------------------------------------------------------------------------

KEEP_TRIAGE_DISABLED: Final = "kept_triage_disabled"
"""The bar was configured off, or narrowed to nothing. Everything through."""

KEEP_RELEVANT_CONSEQUENCE: Final = "kept_relevant_consequence"
"""A coding-altering or splice-relevant term, or HIGH/MODERATE impact. The same
predicate the selection stage applies, asked earlier."""

KEEP_IMPACT_AT_FLOOR: Final = "kept_impact_at_or_above_floor"
"""Impact clears ``min_impact`` without any term having argued for it — in
practice a LOW-impact call under the default bar."""

KEEP_SPLICE_AI: Final = "kept_splice_ai"
"""A SpliceAI delta at or above the threshold, whatever the terms say."""

KEEP_IMPACT_NOT_ASSESSED: Final = "kept_impact_not_assessed"
"""No transcript had an impact computed. NOT ASSESSED is not MODIFIER (ADR 0016)."""

KEEP_CLINICAL_ASSERTION: Final = "kept_clinical_assertion"
"""The record arrived already carrying a curated assertion. Rare — it means the
input was pre-annotated — but free to honour, and the one case where the clinical
rescue is reachable without a lookup."""

KEEP_NO_GENE_ASSIGNMENT: Final = "kept_no_gene_assignment"
"""No consequence annotation at all, retained because the policy says so."""

DROP_IMPACT_BELOW_FLOOR: Final = "dropped_impact_below_floor"
"""Every transcript was assessed, none is coding or splice-relevant, none carries
a SpliceAI signal, and the most severe impact among them is below the floor."""

DROP_NO_GENE_ASSIGNMENT: Final = "dropped_no_gene_assignment"
"""No consequence annotation, so no gene, so no gene-scoped hypothesis. The only
drop reason here that is provable rather than heuristic."""

#: Fixed order, so report keys and their JSON serialisation are stable (GP-30).
KEEP_REASONS: Final[tuple[str, ...]] = (
    KEEP_TRIAGE_DISABLED,
    KEEP_RELEVANT_CONSEQUENCE,
    KEEP_IMPACT_AT_FLOOR,
    KEEP_SPLICE_AI,
    KEEP_IMPACT_NOT_ASSESSED,
    KEEP_CLINICAL_ASSERTION,
    KEEP_NO_GENE_ASSIGNMENT,
)

DROP_REASONS: Final[tuple[str, ...]] = (
    DROP_IMPACT_BELOW_FLOOR,
    DROP_NO_GENE_ASSIGNMENT,
)

#: Rendered into the ``value`` column of a discard row when nothing was assessed.
#: Never ``"modifier"`` and never ``"0"``: the row must not say the annotator
#: made a prediction it did not make (GP-14).
NOT_ASSESSED_VALUE: Final = "not_assessed"


@dataclass(frozen=True, slots=True)
class TriageVerdict:
    """One decision: keep or drop, the rule, and the value the rule fired on."""

    retained: bool
    reason: str
    value: str = ""

    def as_discard(self, variant_id: str) -> DiscardRecord:
        """The ledger row for a dropped variant. Only valid when not retained."""
        return DiscardRecord(
            variant_id=variant_id,
            stage=DiscardStage.TRIAGE,
            rule=self.reason,
            value=self.value,
        )


#: The verdict every variant gets when the bar is off. A module-level constant
#: rather than a fresh object per record: at 5 M records that is 5 M allocations
#: bought for nothing, in the configuration that is meant to cost nothing.
PASS_THROUGH: Final = TriageVerdict(retained=True, reason=KEEP_TRIAGE_DISABLED)


class VariantTriage:
    """The bar, bound to one configuration. Pure, stateless, and cheap per call.

    Holds no counters: an instance is shared across batches and the counting is
    the stream's job. :meth:`verdict` reads a variant's consequence annotations
    and its clinical assertions and touches nothing else, which is what makes it
    answerable before the frequency and clinical adapters have been called.
    """

    __slots__ = (
        "_clinical_before_bar",
        "_enabled",
        "_floor",
        "_retain_no_gene",
        "_retain_unassessed",
        "_splice_ai",
    )

    def __init__(self, thresholds: TriageThresholds) -> None:
        self._enabled = thresholds.enabled
        self._clinical_before_bar = thresholds.clinical_before_bar
        self._floor = thresholds.min_impact
        self._retain_unassessed = thresholds.retain_unassessed_impact
        self._retain_no_gene = thresholds.retain_without_gene_assignment
        self._splice_ai = thresholds.min_splice_ai_delta

    @property
    def enabled(self) -> bool:
        return self._enabled

    def _unannotated_verdict(
        self, clinical_assertions: Sequence[ClinicalAssertion]
    ) -> TriageVerdict:
        """The verdict for a record the consequence adapter said nothing about.

        Separate from :meth:`verdict` because it is the one branch that is provable
        rather than heuristic, and it deserves to be readable on its own: pairing is
        gene-scoped, so a record with no gene symbol forms no candidate pair under
        any configuration. A curated assertion still keeps it — the assertion is a
        fact about the variant, not about the annotator's coverage.
        """
        if clinical_assertions:
            return TriageVerdict(retained=True, reason=KEEP_CLINICAL_ASSERTION)
        if self._retain_no_gene:
            return TriageVerdict(retained=True, reason=KEEP_NO_GENE_ASSIGNMENT)
        return TriageVerdict(retained=False, reason=DROP_NO_GENE_ASSIGNMENT)

    @property
    def clinical_before_bar(self) -> bool:
        """Whether the clinical adapter is asked about every record, before the bar.

        Read by :class:`~mva.annotation.service.AnnotationStream`, which owns the
        adapters. When True the verdict below receives the adapter's answer and the
        clinical rescue is reachable; when False it receives only whatever the input
        record already carried, and the rescue is not.
        """
        return self._clinical_before_bar

    def verdict(
        self,
        consequences: Sequence[ConsequenceAnnotation],
        *,
        clinical_assertions: Sequence[ClinicalAssertion] = (),
    ) -> TriageVerdict:
        """Decide, from the consequence annotations alone (plus any input assertions).

        Checks run in a fixed order, so a variant that clears the bar on two
        grounds records the same reason on every run (GP-30). The order is by
        strength of the argument, not by cost: a reader of the triage report
        should see ``kept_relevant_consequence`` where a coding term is the reason
        and ``kept_impact_at_or_above_floor`` only where nothing better applied.
        """
        if not self._enabled:
            return PASS_THROUGH
        if not consequences:
            return self._unannotated_verdict(clinical_assertions)

        worst: ImpactSeverity | None = None
        splice_ai = False
        for annotation in consequences:
            if consequence_is_selecting(annotation):
                return TriageVerdict(retained=True, reason=KEEP_RELEVANT_CONSEQUENCE)
            impact = annotation.impact
            if impact is not None and (worst is None or impact_rank(impact) < impact_rank(worst)):
                worst = impact
            delta = annotation.splice_ai_delta_max
            if delta is not None and delta >= self._splice_ai:
                splice_ai = True

        if impact_at_least(worst, self._floor):
            # `worst` is not None here: impact_at_least is False for None.
            return TriageVerdict(retained=True, reason=KEEP_IMPACT_AT_FLOOR, value=str(worst))
        if splice_ai:
            return TriageVerdict(retained=True, reason=KEEP_SPLICE_AI)
        if clinical_assertions:
            return TriageVerdict(retained=True, reason=KEEP_CLINICAL_ASSERTION)
        if worst is None:
            # Every transcript reported `impact is None`. Nobody judged this
            # variant, and a bar is not entitled to read silence as a verdict.
            if self._retain_unassessed:
                return TriageVerdict(
                    retained=True, reason=KEEP_IMPACT_NOT_ASSESSED, value=NOT_ASSESSED_VALUE
                )
            return TriageVerdict(
                retained=False, reason=DROP_IMPACT_BELOW_FLOOR, value=NOT_ASSESSED_VALUE
            )
        return TriageVerdict(retained=False, reason=DROP_IMPACT_BELOW_FLOOR, value=str(worst))


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class TriageReport:
    """What the bar did, in counts. Counts and thresholds only — never an id.

    Classified DERIVED_SAFE alongside the QC and selection reports for exactly
    that reason (GP-41). The identities live in the discard ledger, which is
    SENSITIVE.
    """

    input_count: int
    retained_count: int
    retained_by_reason: dict[str, int] = field(default_factory=dict)
    dropped_by_reason: dict[str, int] = field(default_factory=dict)
    thresholds: dict[str, bool | float | str] = field(default_factory=dict)
    narrowings: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()

    @property
    def dropped_count(self) -> int:
        return self.input_count - self.retained_count

    @property
    def retained_fraction(self) -> float:
        return (self.retained_count / self.input_count) if self.input_count else 0.0

    @property
    def speedup_factor(self) -> float:
        """How many times fewer records reach the expensive half.

        A ratio of record counts, not of measured seconds. It is the honest thing
        this report can state: the per-record cost downstream is a tabix seek and
        a JSON row, and both scale with the count — but the wall-clock ratio also
        depends on how the region merge coalesces, which this stage cannot see.
        """
        return (self.input_count / self.retained_count) if self.retained_count else 0.0

    def as_payload(self) -> dict[str, object]:
        """The report artifact's body. Fixed key order; counts only."""
        return {
            "input_count": self.input_count,
            "retained_count": self.retained_count,
            "dropped_count": self.dropped_count,
            "retained_fraction": self.retained_fraction,
            "record_reduction_factor": self.speedup_factor,
            "retained_by_reason": dict(self.retained_by_reason),
            "dropped_by_reason": dict(self.dropped_by_reason),
            "thresholds": dict(self.thresholds),
            "narrowings": list(self.narrowings),
            "warnings": list(self.warnings),
        }


# ---------------------------------------------------------------------------
# The safety coupling
# ---------------------------------------------------------------------------


def narrow_to_selection(
    triage: TriageThresholds, selection: SelectionThresholds
) -> tuple[TriageThresholds, tuple[str, ...]]:
    """Widen the bar toward selection on six of its seven axes.

    **It does not close the clinical axis, and it does not make the bar a superset
    of a disabled selection stage.** Both exceptions are stated below and both are
    warned about at runtime; neither is a defect. They are in this first sentence
    because a reader who trusts the summary and skips the rest gets a guarantee
    this function does not give — and a wrong model of what selection keeps is the
    error that already cost sixteen ClinVar-rescued variants once.

    The honest invariant is not "the bar is a superset of selection". It is **"the
    bar is a superset of (selection ∩ pairable)"**, and the licence for the
    difference comes from pairing being gene-scoped (`pairing._variants_by_gene`),
    not from selection at all.

    This function is the determinism argument in code. The partition is
    answer-preserving if and only if every variant that could reach the output
    still reaches the frequency and clinical adapters — and the two ways it could
    fail are both configuration, not logic:

    * **Selection is switched off.** Then selection retains every valid record and
      the bar has no downstream gate behind it to be wider than. The only drops
      that stay provable are the ones that rest on pairing rather than on
      selection: pairing is gene-scoped, so a variant with no gene symbol forms no
      candidate pair whatever selection does. So the bar collapses to exactly
      that, by forcing ``min_impact`` to MODIFIER and unassessed-impact retention
      on. The synthetic demo case runs in this state deliberately (ADR 0027 §7),
      which is why the demo still exercises the partition rather than skipping it.

    * **Selection is configured to keep records with no gene assignment.** Then
      the two stages would be reporting about different populations, and the bar
      would be answering a question selection had declined to. Retention is forced
      on to match.

    The SpliceAI threshold is clamped for the same reason: a bar stricter than the
    stage behind it drops variants that stage retains on the identical signal.

    Returns the effective thresholds and a list of the narrowings applied, in a
    fixed order. The caller must surface the list — a bar that silently widened is
    a bar whose configured value is a lie, and the run report has to say which one
    was in force.

    What this function does NOT protect is the clinical rescue. That one needs a
    ClinVar lookup, which is the thing being avoided, and no rearrangement of
    configuration reaches it. It is named in the module docstring, in the run
    warnings and in ADR 0028, and it is the design's known cost.
    """
    if not triage.enabled:
        return triage, ()

    narrowings: list[str] = []
    updates: dict[str, object] = {}

    if not selection.enabled:
        if triage.min_impact is not ImpactSeverity.MODIFIER:
            updates["min_impact"] = ImpactSeverity.MODIFIER
            narrowings.append(
                f"min_impact widened from {triage.min_impact} to modifier: candidate "
                "selection is disabled, so every valid record reaches pairing and the "
                "only variants the bar may drop are the ones pairing cannot use at all "
                "(no gene assignment)."
            )
        if not triage.retain_unassessed_impact:
            updates["retain_unassessed_impact"] = True
            narrowings.append(
                "retain_unassessed_impact forced true: candidate selection is disabled, "
                "so nothing downstream would have dropped a NOT ASSESSED impact "
                "(ADR 0016)."
            )
    elif selection.retain_unassessed_impact and not triage.retain_unassessed_impact:
        updates["retain_unassessed_impact"] = True
        narrowings.append(
            "retain_unassessed_impact forced true: selection retains variants whose "
            "impact was never assessed, so the bar may not drop them first (ADR 0016)."
        )

    if not selection.drop_without_gene_assignment and not triage.retain_without_gene_assignment:
        updates["retain_without_gene_assignment"] = True
        narrowings.append(
            "retain_without_gene_assignment forced true: selection is configured to "
            "keep variants with no gene assignment, so the bar may not drop them first."
        )

    if triage.min_splice_ai_delta > selection.min_splice_ai_delta:
        updates["min_splice_ai_delta"] = selection.min_splice_ai_delta
        narrowings.append(
            f"min_splice_ai_delta lowered from {triage.min_splice_ai_delta} to "
            f"{selection.min_splice_ai_delta} to match selection: a stricter bar would "
            "drop variants selection retains on the identical SpliceAI signal."
        )

    if not updates:
        return triage, ()
    return triage.model_copy(update=updates), tuple(narrowings)


def triage_warnings(
    thresholds: TriageThresholds, selection: SelectionThresholds, narrowings: Sequence[str]
) -> tuple[str, ...]:
    """Run-level caveats about the bar, in a fixed order. These belong in the report."""
    warnings: list[str] = []
    if not thresholds.enabled:
        warnings.append(
            "Annotation triage is DISABLED: every annotated variant received a frequency "
            "lookup, a clinical lookup, a full evidence set and a row in the variants "
            "artifact. This is the pre-ADR-0028 pass and the state a reviewer re-runs in "
            "to see what the bar cost; it is not the state a whole-genome case should run "
            "in, because it is the one that projected ~19.5 h and ~95 GB."
        )
        return tuple(warnings)

    for narrowing in narrowings:
        warnings.append(f"Annotation triage narrowed against candidate selection: {narrowing}")

    if selection.retain_pathogenic_clinical and not thresholds.clinical_before_bar:
        warnings.append(
            "triage.clinical_before_bar is FALSE while selection.retain_pathogenic_clinical "
            "is TRUE. A curated pathogenic assertion therefore cannot rescue a variant the "
            "bar drops, because the rescue needs a lookup the bar skipped. The exposed class "
            "is a variant that is MODIFIER on every transcript, carries no coding or splice "
            "term and no SpliceAI signal, and holds a pathogenic ClinVar record — deep "
            "intronic pathogenic alleles have exactly that shape, and on the scale phantom "
            "this setting changed the selected set from 3,616 to 3,600. Every such variant "
            "is in the discard ledger by id, so a targeted ClinVar scan over the dropped ids "
            "closes it without a re-run; setting clinical_before_bar back to true closes it "
            "outright (ADR 0028)."
        )
    return tuple(warnings)
