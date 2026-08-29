"""Variants the relevance bar must NOT drop, and one it must (ADR 0028).

`test_staged_annotation_determinism.py` proves that switching the bar on does not
reorder or corrupt what survives. It cannot prove the bar is *wide enough*, and
the reason is in its fixtures: the shipped synthetic case is twelve records with
one MODIFIER-only variant, no deep-intronic allele, and
:class:`~mva.annotation.local_tables.NullClinicalAdapter` in the clinical slot —
an adapter that returns nothing for every id.

That gap is not hypothetical. The first implementation of ADR 0028 put the
ClinVar lookup *behind* the bar and lost sixteen variants on the scale phantom —
exactly the ones ``SelectionThresholds.retain_pathogenic_clinical`` exists to
rescue. Because the synthetic case's clinical adapter cannot produce an
assertion, **every byte-identity test stayed green for the whole time that bug
was live.** A partition test whose clinical adapter has nothing to say cannot
observe the clinical rescue at all.

This module supplies the missing half: a small set of synthetic variants, each
constructed so that exactly one of the bar's keep rules stands between it and the
discard ledger, driven through the real annotation stream and the real selection
stage with a clinical adapter that actually answers.

Five of them are the shapes named in the review of this branch:

1. MODIFIER on every transcript, ClinVar-pathogenic — the bug that occurred.
2. Deep intronic, MODIFIER, SpliceAI above threshold — the splice-disrupting
   variant no consequence term flags.
3. MODIFIER on every transcript in a gene on the disease panel — which the bar
   does *not* rescue, and neither does selection. See
   :func:`test_a_modifier_only_variant_in_a_panel_gene_is_dropped_by_both_stages`;
   that parity is the finding, and it is asserted rather than papered over.
4. Impact never assessed on any transcript (ADR 0016, GP-14).
5. MODIFIER on every transcript with nothing arguing for it — the positive
   control. A suite in which everything survives proves nothing about the bar.

Beyond the fixtures, three checks answer the question ``narrow_to_selection``
raises about itself — that it "proves the bar is wider than selection *as I
modelled selection*", by a model that was wrong once:

* a **configuration grid** (7 bar configurations x 9 selection configurations)
  that runs the real selection stage on both sides of the partition and compares
  the selected sets, rather than reading the coupling's model of it;
* a **shape grid** (625 consequence shapes, with and without a curated
  assertion) that does the same over the data rather than over the knobs — which
  matters because only ``consequence_is_selecting`` is actually shared between
  the two stages; "never assessed" and "SpliceAI above threshold" are each
  written twice, once per stage, in different shapes;
* an explicit test that the coupling does **not** close the clinical axis, so the
  one gap it leaves is pinned rather than described.

Nothing here reads a repository knowledge table or a case volume. Every gene is
``SYNTH``-prefixed, every coordinate is fabricated, and the adapters are declared
synthetic (GP-20, GP-40).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import pytest

from mva.annotation.base import AdapterSet
from mva.annotation.service import iter_annotated
from mva.annotation.triage import (
    DROP_IMPACT_BELOW_FLOOR,
    DROP_NO_GENE_ASSIGNMENT,
    KEEP_CLINICAL_ASSERTION,
    KEEP_IMPACT_AT_FLOOR,
    KEEP_IMPACT_NOT_ASSESSED,
    KEEP_NO_GENE_ASSIGNMENT,
    KEEP_RELEVANT_CONSEQUENCE,
    KEEP_SPLICE_AI,
    VariantTriage,
    narrow_to_selection,
)
from mva.clock import Clock, FixedClock
from mva.config import FrequencyThresholds, SelectionThresholds, TriageThresholds
from mva.evidence.discards import DiscardLedger, explain_variant
from mva.models.discard import DiscardStage
from mva.models.genome import GenomeBuild, GenomicCoordinate
from mva.models.variant import (
    ClinicalAssertion,
    ConsequenceAnnotation,
    FilterStatus,
    Genotype,
    ImpactSeverity,
    PopulationFrequency,
    VariantRecord,
    Zygosity,
)
from mva.prioritization.selection import iter_selected

pytestmark = [pytest.mark.integration]

#: Fabricated contig and position block. Canonical contig so the record is
#: well-formed, positions chosen not to collide with the shipped synthetic case.
_CONTIG = "chr21"
_FIRST_POSITION = 10_000_100

#: The demo gene panel's first row (`knowledge/public/gene_panel.tsv`). Used by
#: fixture 3 so that "the gene is on the disease panel" is a fact about this
#: repository's own panel table rather than an assumption in a test.
PANEL_GENE = "SYNTHKIN1"


# ---------------------------------------------------------------------------
# Fixture adapters
#
# The clinical one is the point of this module. `NullClinicalAdapter` returns {}
# for every id, so with it bound the bar's clinical rule and selection's
# `retain_pathogenic_clinical` are both unreachable and both untested.
# ---------------------------------------------------------------------------


class FixtureConsequenceAdapter:
    """Serves a fabricated per-variant consequence table. Declares itself synthetic."""

    def __init__(self, table: Mapping[str, tuple[ConsequenceAnnotation, ...]]) -> None:
        self._table = dict(table)

    @property
    def name(self) -> str:
        return "synthetic-consequence-fixture"

    @property
    def version(self) -> str:
        return "v0-negative-fixtures"

    @property
    def synthetic(self) -> bool:
        return True

    def annotate(
        self, variant_ids: Sequence[str]
    ) -> Mapping[str, tuple[ConsequenceAnnotation, ...]]:
        """Keys only for ids the table knows. A missing key is NOT 'no consequence'."""
        return {vid: self._table[vid] for vid in variant_ids if vid in self._table}


class FixtureFrequencyAdapter:
    """Knows nothing about anything, and records how often it was asked.

    Two jobs. It leaves every fixture variant with *no* frequency data, so
    selection retains them all under ``retain_unknown_frequency`` and the only
    gate that can decide their fate is the consequence rule this module is about
    (GP-14: unknown is not rare and not common). And ``asked`` is the measurement
    ADR 0028 is built on — the count of records that reached the expensive half.
    """

    def __init__(self) -> None:
        self.asked: list[str] = []

    @property
    def name(self) -> str:
        return "synthetic-frequency-fixture"

    @property
    def version(self) -> str:
        return "v0-negative-fixtures"

    @property
    def synthetic(self) -> bool:
        return True

    def frequencies(
        self, variant_ids: Sequence[str]
    ) -> Mapping[str, tuple[PopulationFrequency, ...]]:
        self.asked.extend(variant_ids)
        return {}


class FixtureClinicalAdapter:
    """A clinical adapter that actually returns an assertion.

    The whole reason this file exists. It also records the ids it was asked
    about, which is how the tests below tell "the lookup ran in front of the bar"
    from "the lookup ran behind it".
    """

    def __init__(self, table: Mapping[str, tuple[ClinicalAssertion, ...]]) -> None:
        self._table = dict(table)
        self.asked: list[str] = []

    @property
    def name(self) -> str:
        return "synthetic-clinical-fixture"

    @property
    def version(self) -> str:
        return "v0-negative-fixtures"

    @property
    def synthetic(self) -> bool:
        return True

    def assertions(self, variant_ids: Sequence[str]) -> Mapping[str, tuple[ClinicalAssertion, ...]]:
        self.asked.extend(variant_ids)
        return {vid: self._table[vid] for vid in variant_ids if vid in self._table}


# ---------------------------------------------------------------------------
# The fixtures
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BarFixture:
    """One synthetic variant and the single bar rule that decides it.

    ``expected_triage_reason`` is asserted, not just ``must_survive``: a variant
    that survives for a *different* reason than the one it was built to exercise
    means the fixture stopped testing what it claims to test.
    """

    key: str
    what_it_targets: str
    gene: str | None
    consequences: tuple[ConsequenceAnnotation, ...]
    assertions: tuple[ClinicalAssertion, ...] = ()
    must_survive: bool = True
    #: Whether the variant must also clear the SELECTION stage. Distinct from
    #: ``must_survive`` on purpose: the bar is deliberately WIDER than selection,
    #: so a variant can clear the bar and still be dropped later on its merits.
    #: That is the permitted direction, and separating the two fields is what
    #: stops a selection drop from being read as a partition bug.
    reaches_selection: bool = True
    expected_triage_reason: str = ""
    expected_drop_value: str = ""
    position: int = field(default=0)

    @property
    def variant_id(self) -> str:
        return f"{GenomeBuild.GRCH38}:{_CONTIG}:{self.position}:A:G"


def _csq(
    gene: str,
    terms: tuple[str, ...],
    impact: ImpactSeverity | None,
    *,
    transcript: str = "SYNTHT0001.1",
    splice_ai: float | None = None,
) -> ConsequenceAnnotation:
    return ConsequenceAnnotation(
        gene_symbol=gene,
        transcript_id=transcript,
        consequence_terms=terms,
        impact=impact,
        splice_ai_delta_max=splice_ai,
        source_tool="synthetic-consequence-fixture",
        source_tool_version="v0-negative-fixtures",
    )


def _pathogenic() -> ClinicalAssertion:
    return ClinicalAssertion(
        source="SYNTH_curation",
        version="v0-negative-fixtures",
        accession="SYNTHVAR000001",
        significance="Pathogenic",
        review_status="criteria_provided,_single_submitter",
        star_rating=1,
        conditions=("Synthetic variegated aneuploidy (demo)",),
    )


#: Ordered so that ``position`` is assigned once, below, and never by hand.
_FIXTURE_SPECS: tuple[BarFixture, ...] = (
    BarFixture(
        key="modifier_only_clinvar_pathogenic",
        what_it_targets=(
            "MODIFIER on every transcript, no coding or splice term, no SpliceAI — "
            "rescued only by a curated pathogenic assertion. The exact shape of the "
            "3,616 -> 3,600 bug (ADR 0028 decision 6)."
        ),
        gene="SYNTHBAR1",
        consequences=(
            _csq("SYNTHBAR1", ("intron_variant",), ImpactSeverity.MODIFIER),
            _csq(
                "SYNTHBAR1",
                ("intron_variant",),
                ImpactSeverity.MODIFIER,
                transcript="SYNTHT0002.1",
            ),
        ),
        assertions=(_pathogenic(),),
        expected_triage_reason=KEEP_CLINICAL_ASSERTION,
    ),
    BarFixture(
        key="deep_intronic_high_spliceai",
        what_it_targets=(
            "Deep intronic, MODIFIER on every transcript, no splice SO term — the "
            "splice-disrupting allele that only a SpliceAI delta argues for."
        ),
        gene="SYNTHBAR2",
        consequences=(
            _csq(
                "SYNTHBAR2",
                ("intron_variant",),
                ImpactSeverity.MODIFIER,
                splice_ai=0.62,
            ),
        ),
        expected_triage_reason=KEEP_SPLICE_AI,
    ),
    BarFixture(
        key="modifier_only_panel_gene",
        what_it_targets=(
            "MODIFIER on every transcript, in a gene on knowledge/public/gene_panel.tsv. "
            "Panel membership is not a keep rule in EITHER stage, so this is dropped — "
            "the finding, asserted rather than forced green."
        ),
        gene=PANEL_GENE,
        consequences=(_csq(PANEL_GENE, ("intron_variant",), ImpactSeverity.MODIFIER),),
        must_survive=False,
        expected_triage_reason=DROP_IMPACT_BELOW_FLOOR,
        expected_drop_value="modifier",
    ),
    BarFixture(
        key="impact_never_assessed",
        what_it_targets=(
            "Impact is None on every transcript: the annotator located the gene and "
            "computed nothing. NOT ASSESSED is not MODIFIER (ADR 0016, GP-14)."
        ),
        gene="SYNTHBAR4",
        consequences=(
            _csq("SYNTHBAR4", ("intron_variant",), None),
            _csq("SYNTHBAR4", ("upstream_gene_variant",), None, transcript="SYNTHT0002.1"),
        ),
        expected_triage_reason=KEEP_IMPACT_NOT_ASSESSED,
    ),
    BarFixture(
        key="modifier_only_no_rescue",
        what_it_targets=(
            "The positive control for the bar itself: MODIFIER everywhere, no coding "
            "or splice term, no SpliceAI signal, no assertion. Nothing argues for it "
            "and it must land in the ledger."
        ),
        gene="SYNTHBAR5",
        consequences=(_csq("SYNTHBAR5", ("intron_variant",), ImpactSeverity.MODIFIER),),
        must_survive=False,
        expected_triage_reason=DROP_IMPACT_BELOW_FLOOR,
        expected_drop_value="modifier",
    ),
    BarFixture(
        key="no_consequence_at_all",
        what_it_targets=(
            "The consequence adapter returned no key at all. The one drop that is "
            "provable rather than heuristic — pairing is gene-scoped, so this record "
            "forms no candidate pair under any configuration."
        ),
        gene=None,
        consequences=(),
        must_survive=False,
        expected_triage_reason=DROP_NO_GENE_ASSIGNMENT,
        expected_drop_value="",
    ),
    BarFixture(
        key="splice_region_synonymous",
        what_it_targets=(
            "The class ADR 0028 says the pipeline can least afford to lose: a "
            "synonymous call that also carries a splice term, at LOW impact."
        ),
        gene="SYNTHBAR7",
        consequences=(
            _csq(
                "SYNTHBAR7",
                ("splice_region_variant&synonymous_variant",),
                ImpactSeverity.LOW,
            ),
        ),
        expected_triage_reason=KEEP_RELEVANT_CONSEQUENCE,
    ),
    BarFixture(
        key="low_impact_synonymous",
        what_it_targets=(
            "LOW impact with no coding or splice term. The bar keeps it and selection "
            "drops it, which is the direction the coupling permits: wider is safe."
        ),
        gene="SYNTHBAR8",
        consequences=(_csq("SYNTHBAR8", ("synonymous_variant",), ImpactSeverity.LOW),),
        reaches_selection=False,
        expected_triage_reason=KEEP_IMPACT_AT_FLOOR,
    ),
    BarFixture(
        key="moderate_missense",
        what_it_targets=(
            "An ordinary MODERATE missense. Present so that raising triage.min_impact "
            "to HIGH is a test of something: this variant must survive a floor it does "
            "not clear, on its SO term alone."
        ),
        gene="SYNTHBAR9",
        consequences=(_csq("SYNTHBAR9", ("missense_variant",), ImpactSeverity.MODERATE),),
        expected_triage_reason=KEEP_RELEVANT_CONSEQUENCE,
    ),
    BarFixture(
        key="partially_assessed",
        what_it_targets=(
            "One transcript never assessed, one predicted MODIFIER. Both stages read "
            "'unassessed' as ALL transcripts unassessed, so both drop this. The parity "
            "is what is under test; the reading itself is a shared decision (ADR 0016)."
        ),
        gene="SYNTHBAR10",
        consequences=(
            _csq("SYNTHBAR10", ("intron_variant",), None),
            _csq(
                "SYNTHBAR10",
                ("intron_variant",),
                ImpactSeverity.MODIFIER,
                transcript="SYNTHT0002.1",
            ),
        ),
        must_survive=False,
        expected_triage_reason=DROP_IMPACT_BELOW_FLOOR,
        expected_drop_value="modifier",
    ),
)

FIXTURES: tuple[BarFixture, ...] = tuple(
    BarFixture(
        key=spec.key,
        what_it_targets=spec.what_it_targets,
        gene=spec.gene,
        consequences=spec.consequences,
        assertions=spec.assertions,
        must_survive=spec.must_survive,
        reaches_selection=spec.reaches_selection,
        expected_triage_reason=spec.expected_triage_reason,
        expected_drop_value=spec.expected_drop_value,
        position=_FIRST_POSITION + 100 * index,
    )
    for index, spec in enumerate(_FIXTURE_SPECS)
)

BY_KEY: Mapping[str, BarFixture] = {fixture.key: fixture for fixture in FIXTURES}


def _record(fixture: BarFixture, line_index: int) -> VariantRecord:
    return VariantRecord(
        coordinate=GenomicCoordinate(
            build=GenomeBuild.GRCH38,
            contig=_CONTIG,
            position=fixture.position,
            ref="A",
            alt="G",
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
        source_artifact="tests/integration/test_relevance_bar_negative_fixtures.py",
        source_line_index=line_index,
    )


def _adapters() -> tuple[AdapterSet, FixtureFrequencyAdapter, FixtureClinicalAdapter]:
    frequency = FixtureFrequencyAdapter()
    clinical = FixtureClinicalAdapter(
        {f.variant_id: f.assertions for f in FIXTURES if f.assertions}
    )
    consequence = FixtureConsequenceAdapter(
        {f.variant_id: f.consequences for f in FIXTURES if f.consequences}
    )
    return (
        AdapterSet(consequence=consequence, frequency=frequency, clinical=clinical),
        frequency,
        clinical,
    )


# ---------------------------------------------------------------------------
# The harness: the real annotation stream, then the real selection stage
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Outcome:
    """What one configuration did, end to end."""

    selected: tuple[str, ...]
    annotated: tuple[str, ...]
    ledger: Path
    #: ``keep reason -> count`` from the bar's own report, so a test can assert
    #: WHICH rule kept a variant and not merely that something did.
    triage_keeps: Mapping[str, int]
    frequency_asked: tuple[str, ...]
    clinical_asked: tuple[str, ...]
    narrowings: tuple[str, ...]


def _clock() -> Clock:
    return FixedClock(datetime(2026, 1, 1, tzinfo=UTC))


def _run(
    workdir: Path,
    *,
    bar: bool,
    triage_thresholds: TriageThresholds | None = None,
    selection_thresholds: SelectionThresholds | None = None,
) -> Outcome:
    """Annotate, then select, writing every discard to one ledger.

    Deliberately assembled from the same two public entry points the composition
    root uses — :func:`iter_annotated` and :func:`iter_selected` — rather than
    from the bar's ``verdict`` method. A unit test of ``verdict`` cannot see the
    thing that broke: whether the clinical adapter's answer reaches it.
    """
    selection = selection_thresholds or SelectionThresholds()
    configured = triage_thresholds or TriageThresholds()
    effective, narrowings = narrow_to_selection(configured, selection)

    adapters, frequency, clinical = _adapters()
    records = [_record(fixture, index) for index, fixture in enumerate(FIXTURES)]

    ledger_path = workdir / "discards" / "discard_ledger.tsv"
    ledger = DiscardLedger(ledger_path)
    try:
        stream = iter_annotated(
            records,
            adapters=adapters,
            clock=_clock(),
            batch_size=None,
            triage=VariantTriage(effective) if bar else None,
            discards=ledger,
        )
        annotated = [annotation.variant for annotation in stream]
        selection_stream = iter_selected(
            annotated,
            frequency=FrequencyThresholds(),
            thresholds=selection,
            clock=_clock(),
            discards=ledger,
        )
        selected = tuple(variant.variant_id for variant in selection_stream)
        report = stream.triage_report()
    finally:
        ledger.close()

    return Outcome(
        selected=selected,
        annotated=tuple(variant.variant_id for variant in annotated),
        ledger=ledger_path,
        triage_keeps=dict(report.retained_by_reason),
        frequency_asked=tuple(frequency.asked),
        clinical_asked=tuple(clinical.asked),
        narrowings=narrowings,
    )


def _triage_reason(fixture: BarFixture) -> str:
    """The bar's own verdict for one fixture, under the default thresholds."""
    verdict = VariantTriage(TriageThresholds()).verdict(
        fixture.consequences, clinical_assertions=fixture.assertions
    )
    return verdict.reason


# ---------------------------------------------------------------------------
# 1-4: variants that must survive the bar
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "key",
    [f.key for f in FIXTURES if f.must_survive],
)
def test_a_variant_that_must_survive_the_bar_reaches_the_selected_set(
    key: str, tmp_path: Path
) -> None:
    """Each keep rule, exercised by a variant that depends on it alone.

    Two claims, kept apart. Clearing the bar means reaching the frequency lookup
    and the annotated artifact — that is the bar's whole promise, and it is
    asserted for every fixture here. Reaching the *selected* set is a second,
    stronger claim, and it is asserted only where selection would keep the variant
    on its own merits: the bar is deliberately wider than selection, so a LOW-impact
    synonymous call clearing the bar and then being dropped by
    ``dropped_not_coding_or_splice`` is the coupling working, not failing.

    What is asserted in every case is that no TRIAGE row exists for the variant.
    A loss attributed to the wrong stage is the failure mode a ledger exists to
    prevent.
    """
    fixture = BY_KEY[key]
    outcome = _run(tmp_path, bar=True)

    assert fixture.variant_id in outcome.annotated, (
        f"[{key}] {fixture.what_it_targets}\n"
        "The relevance bar dropped this variant before the frequency lookup, so it "
        "has no evidence, no row in variants/annotated.json and no chance of being "
        "ranked. The bar is too narrow (ADR 0028)."
    )
    reason = _triage_reason(fixture)
    assert reason == fixture.expected_triage_reason, (
        f"[{key}] survived the bar for reason {reason!r}, not the "
        f"{fixture.expected_triage_reason!r} it was built to exercise. The fixture "
        "has stopped testing the rule it claims to test."
    )
    assert outcome.triage_keeps.get(reason, 0) >= 1, (
        f"[{key}] the bar kept it under {reason!r} but the triage report counts none. "
        "The report is what the run publishes; a counter that disagrees with the rule "
        "is a report that describes a run that did not happen."
    )
    rows = explain_variant(outcome.ledger, fixture.variant_id)
    assert not [row for row in rows if row.stage is DiscardStage.TRIAGE], (
        f"[{key}] cleared the bar and still has a triage row in the discard ledger"
    )
    if fixture.reaches_selection:
        assert fixture.variant_id in outcome.selected, (
            f"[{key}] cleared the bar and then vanished at selection, which this "
            "fixture was not built to exercise"
        )
        assert not rows, f"[{key}] is in the discard ledger AND in the selected set"
    else:
        assert [row.stage for row in rows] == [DiscardStage.SELECTION], (
            f"[{key}] is meant to be dropped by SELECTION on its own merits, and the "
            f"ledger instead says {[str(row.stage) for row in rows]}"
        )


def test_a_modifier_only_clinvar_pathogenic_variant_survives_the_bar(tmp_path: Path) -> None:
    """Fixture 1. The bug that occurred, in the configuration that hid it.

    MODIFIER on both transcripts, no coding or splice term, no SpliceAI signal.
    The only thing between it and the ledger is a ClinVar assertion, and the only
    way the bar can see one is if the clinical lookup ran in front of it. With
    `NullClinicalAdapter` bound — as in the shipped synthetic case — this
    assertion cannot fail, because there is no assertion to lose.
    """
    fixture = BY_KEY["modifier_only_clinvar_pathogenic"]
    outcome = _run(tmp_path, bar=True)

    assert fixture.variant_id in outcome.selected, (
        "a MODIFIER-only variant with a curated PATHOGENIC assertion did not reach "
        "the selected set. This is the exact class selection.retain_pathogenic_clinical "
        "exists to rescue, and the 3,616 -> 3,600 regression ADR 0028 decision 6 "
        "reverted."
    )
    assert fixture.variant_id in outcome.clinical_asked, (
        "the clinical adapter was never asked about this variant, so the rescue "
        "cannot have fired on evidence — check triage.clinical_before_bar"
    )


def test_a_deep_intronic_variant_survives_on_its_spliceai_delta_alone(tmp_path: Path) -> None:
    """Fixture 2. No SO term flags it; only the delta does.

    This is the variant class the module docstring of `triage` calls the one the
    pipeline can least afford to lose, and it is absent from the twelve-record
    fixture entirely.
    """
    fixture = BY_KEY["deep_intronic_high_spliceai"]
    outcome = _run(tmp_path, bar=True)
    assert fixture.variant_id in outcome.selected
    assert all(
        annotation.impact is ImpactSeverity.MODIFIER for annotation in fixture.consequences
    ), "the fixture stopped being MODIFIER-only, so it no longer isolates the SpliceAI rule"
    terms = [term for annotation in fixture.consequences for term in annotation.consequence_terms]
    assert not any("splice" in term for term in terms), (
        "the fixture gained a splice SO term, so it now survives on the term rather "
        "than on the delta"
    )


def test_a_modifier_only_variant_in_a_panel_gene_is_dropped_by_both_stages(
    tmp_path: Path,
) -> None:
    """Fixture 3. Panel membership is NOT a keep rule, and this records that.

    `knowledge/public/gene_panel.tsv` lists SYNTHKIN1, and nothing in `src/mva/`
    reads that table: `config.py` names the path and no stage opens it, which
    `knowledge/disease/mva_panel_sources.md` §7 states in as many words. The gene
    table the pipeline DOES consume, `gene_phenotype.tsv`, is loaded by the
    composition root *after* pairing and feeds `score_all_genes` — it is a ranking
    input, strictly downstream of both gates, and a variant dropped at either one
    never reaches it.

    So a MODIFIER-only variant in a panel gene is dropped by the bar, and would
    have been dropped by `selection` under `dropped_not_coding_or_splice` anyway.
    That parity is what makes it *not* a partition bug: the bar is not narrower
    than the stage behind it. It is a gap in both — inherited from ADR 0019, which
    ADR 0028's own Consequences section says in as many words — and it is asserted
    here so that the day someone teaches selection to rescue panel genes, this
    test fails and points at the bar that also needs teaching.
    """
    fixture = BY_KEY["modifier_only_panel_gene"]
    with_bar = _run(tmp_path / "on", bar=True)
    without_bar = _run(tmp_path / "off", bar=False)

    assert fixture.variant_id not in with_bar.selected
    assert fixture.variant_id not in without_bar.selected, (
        "selection now retains a MODIFIER-only variant in a panel gene while the "
        "relevance bar still drops it. The bar has become narrower than the stage "
        "behind it and narrow_to_selection does not cover this axis (ADR 0028)."
    )
    (row,) = explain_variant(with_bar.ledger, fixture.variant_id)
    assert row.stage is DiscardStage.TRIAGE
    assert row.rule == DROP_IMPACT_BELOW_FLOOR
    assert row.value == "modifier"


def test_a_variant_whose_impact_was_never_assessed_survives_the_bar(tmp_path: Path) -> None:
    """Fixture 4. ADR 0016 / GP-14: nobody judged it, so the bar may not judge it.

    `impact is None` on every transcript is the shape a MANE interval join
    produces for every variant it places. Reading that silence as MODIFIER is the
    fabricated-negative failure GP-14 exists to prevent.
    """
    fixture = BY_KEY["impact_never_assessed"]
    outcome = _run(tmp_path, bar=True)
    assert fixture.variant_id in outcome.selected
    assert all(annotation.impact is None for annotation in fixture.consequences)


# ---------------------------------------------------------------------------
# 5: the positive control, and what the ledger says about it
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("key", [f.key for f in FIXTURES if not f.must_survive])
def test_a_variant_that_should_be_dropped_is_named_in_the_ledger_with_its_rule(
    key: str, tmp_path: Path
) -> None:
    """The bar's positive control. A suite where everything survives proves nothing.

    Every assertion here is about the *ledger*, not about the counter: a reviewer
    asking "why isn't variant X in the output?" has only this file, so the rule
    and the value it fired on have to be reconstructable from its bytes.
    """
    fixture = BY_KEY[key]
    outcome = _run(tmp_path, bar=True)

    assert fixture.variant_id not in outcome.annotated, (
        f"[{key}] {fixture.what_it_targets}\n"
        "This variant was supposed to exercise the drop path and did not, so every "
        "other assertion in this module is about a bar that never fires."
    )
    assert fixture.variant_id not in outcome.frequency_asked, (
        f"[{key}] was dropped by the bar and still reached the frequency adapter, "
        "which is the one call ADR 0028 exists to avoid"
    )

    rows = explain_variant(outcome.ledger, fixture.variant_id)
    assert rows, (
        f"[{key}] was dropped and is NOT in the discard ledger. A dropped candidate "
        "that is not accounted for is invisible and unfalsifiable (GP-13 as amended, "
        "ADR 0028)."
    )
    (row,) = rows
    assert row.stage is DiscardStage.TRIAGE
    assert row.rule == fixture.expected_triage_reason, (
        f"[{key}] the ledger names rule {row.rule!r}; the bar decided on "
        f"{fixture.expected_triage_reason!r}. A ledger that names the wrong rule is "
        "worse than no ledger."
    )
    assert row.value == fixture.expected_drop_value, (
        f"[{key}] the ledger's value column says {row.value!r}, expected "
        f"{fixture.expected_drop_value!r}. An absence must render as an absence, "
        "never as a zero or as a prediction nobody made (GP-14)."
    )


def test_the_ledger_never_renders_an_absence_as_a_prediction(tmp_path: Path) -> None:
    """GP-14 in the one column where it is cheap to get wrong.

    A variant with no consequence annotation was never assessed by anything. Its
    row must carry an EMPTY value — not ``modifier``, not ``0``, not ``none``.
    """
    fixture = BY_KEY["no_consequence_at_all"]
    outcome = _run(tmp_path, bar=True)
    (row,) = explain_variant(outcome.ledger, fixture.variant_id)
    assert row.rule == DROP_NO_GENE_ASSIGNMENT
    assert row.value == ""
    assert row.value.lower() not in {"modifier", "0", "none", "null"}


def test_an_absent_ledger_row_is_ambiguous_rather_than_good_news(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """GP-14 on the read side: `explain` returning nothing is not "correctly kept".

    The ledger records removals, not a census. A variant absent from it either
    survived every gate or was never in the input, and this file cannot tell them
    apart — so the CLI must say both, in words, rather than print a reassuring
    silence.
    """
    from mva.cli import discards_explain

    outcome = _run(tmp_path, bar=True)
    survivor = BY_KEY["deep_intronic_high_spliceai"]
    never_seen = f"{GenomeBuild.GRCH38}:{_CONTIG}:99999999:A:G"

    for variant_id in (survivor.variant_id, never_seen):
        assert explain_variant(outcome.ledger, variant_id) == ()
        discards_explain(outcome.ledger, variant_id)
        rendered = " ".join(capsys.readouterr().out.split())
        assert "no discard record" in rendered
        assert "survived every gate" in rendered
        assert "never in the input" in rendered, (
            "an empty explain result must be reported as ambiguous. A variant that "
            "survived and a variant that was never in the callset produce the same "
            "empty result, and the ledger cannot distinguish them (GP-14)."
        )


# ---------------------------------------------------------------------------
# The bar is not vacuous, and it does what ADR 0028 claims it does
# ---------------------------------------------------------------------------


def test_the_bar_keeps_the_frequency_lookup_away_from_what_it_dropped(tmp_path: Path) -> None:
    """The saving, measured on this fixture set rather than asserted.

    ADR 0028's whole claim is that fewer records reach the gnomAD lookup. If the
    bar drops a record but the adapter is still asked about it, the partition is
    cosmetic.
    """
    with_bar = _run(tmp_path / "on", bar=True)
    without_bar = _run(tmp_path / "off", bar=False)

    dropped = set(without_bar.annotated) - set(with_bar.annotated)
    assert dropped, "the bar dropped nothing, so this fixture set proves nothing about it"
    assert not dropped & set(with_bar.frequency_asked)
    assert len(with_bar.frequency_asked) < len(without_bar.frequency_asked)


def test_the_clinical_lookup_runs_in_front_of_the_bar_by_default(tmp_path: Path) -> None:
    """Decision 6 of ADR 0028, asserted on the adapter's own call log.

    Every id — including the ones the bar is about to drop — reaches the clinical
    adapter, because the verdict needs its answer. The frequency adapter sees only
    the survivors. That asymmetry IS the decision.
    """
    outcome = _run(tmp_path, bar=True)
    for fixture in FIXTURES:
        assert fixture.variant_id in outcome.clinical_asked, (
            f"[{fixture.key}] the clinical adapter was not asked about this variant "
            "before the bar decided, so a curated assertion could not have rescued it"
        )
    assert len(outcome.frequency_asked) < len(outcome.clinical_asked)


def test_deferring_the_clinical_lookup_loses_exactly_the_rescued_variant(
    tmp_path: Path,
) -> None:
    """`clinical_before_bar: false` reproduces the 3,616 -> 3,600 regression here.

    This is the assertion the shipped synthetic case cannot make. It is not a
    complaint about the setting — ADR 0028 offers it deliberately, and
    `triage_warnings` says what it costs — but the cost has to be demonstrable,
    and with `NullClinicalAdapter` bound it is invisible.

    Note what is NOT asserted: that `narrow_to_selection` prevents this. It does
    not, and says so. This test is what stands in for that missing coupling.
    """
    deferred = TriageThresholds(clinical_before_bar=False)
    outcome = _run(tmp_path, bar=True, triage_thresholds=deferred)
    rescued = BY_KEY["modifier_only_clinvar_pathogenic"]

    assert rescued.variant_id not in outcome.selected, (
        "triage.clinical_before_bar=false no longer loses the clinically rescued "
        "variant. Either the hole was closed — in which case delete this test and "
        "the warning in triage_warnings — or this fixture stopped depending on the "
        "clinical rule."
    )
    (row,) = explain_variant(outcome.ledger, rescued.variant_id)
    assert row.rule == DROP_IMPACT_BELOW_FLOOR
    assert row.value == "modifier", (
        "the ledger row for a clinically-rescuable variant lost to clinical_before_bar "
        "must still name the value the bar read, so a targeted ClinVar scan over the "
        "dropped ids can close the hole without a re-run (ADR 0028)"
    )


# ---------------------------------------------------------------------------
# Drift: the bar's definition of "relevant" against selection's
# ---------------------------------------------------------------------------

#: Selection configurations to check the coupling against. Each one moves a knob
#: `narrow_to_selection` claims to defend, plus the two that only widen selection
#: (frequency), which the bar cannot see at all and must therefore not depend on.
_SELECTION_GRID: tuple[tuple[str, SelectionThresholds], ...] = (
    ("default", SelectionThresholds()),
    ("selection_disabled", SelectionThresholds(enabled=False)),
    ("keeps_unannotated", SelectionThresholds(drop_without_gene_assignment=False)),
    ("drops_unassessed_impact", SelectionThresholds(retain_unassessed_impact=False)),
    ("no_clinical_rescue", SelectionThresholds(retain_pathogenic_clinical=False)),
    ("looser_spliceai", SelectionThresholds(min_splice_ai_delta=0.05)),
    ("stricter_spliceai", SelectionThresholds(min_splice_ai_delta=0.8)),
    ("drops_unknown_frequency", SelectionThresholds(retain_unknown_frequency=False)),
    (
        "disabled_and_keeps_unannotated",
        SelectionThresholds(enabled=False, drop_without_gene_assignment=False),
    ),
)

#: Bar configurations. `clinical_before_bar` is deliberately absent: it is the ONE
#: axis on which the bar can be narrower than selection, `narrow_to_selection`
#: does not defend it, and it has its own test above
#: (:func:`test_deferring_the_clinical_lookup_loses_exactly_the_rescued_variant`)
#: which asserts the loss rather than hiding it. Every other axis belongs here,
#: including the three the composition root never sets, because a knob nobody
#: proved safe is a knob somebody will turn.
_TRIAGE_GRID: tuple[tuple[str, TriageThresholds], ...] = (
    ("default", TriageThresholds()),
    ("floor_moderate", TriageThresholds(min_impact=ImpactSeverity.MODERATE)),
    ("floor_high", TriageThresholds(min_impact=ImpactSeverity.HIGH)),
    ("floor_modifier", TriageThresholds(min_impact=ImpactSeverity.MODIFIER)),
    ("drops_unassessed", TriageThresholds(retain_unassessed_impact=False)),
    ("keeps_unannotated", TriageThresholds(retain_without_gene_assignment=True)),
    ("stricter_spliceai", TriageThresholds(min_splice_ai_delta=0.9)),
)


@pytest.mark.parametrize(
    ("selection_label", "selection"), _SELECTION_GRID, ids=[g[0] for g in _SELECTION_GRID]
)
@pytest.mark.parametrize(("triage_label", "triage"), _TRIAGE_GRID, ids=[g[0] for g in _TRIAGE_GRID])
def test_the_bar_drops_nothing_pairable_that_selection_would_have_kept(
    selection_label: str,
    selection: SelectionThresholds,
    triage_label: str,
    triage: TriageThresholds,
    tmp_path: Path,
) -> None:
    """The coupling `narrow_to_selection` exists to guarantee, checked end to end.

    `narrow_to_selection` is the branch author's own model of the selection stage,
    and their model was wrong once already: they read `retain_pathogenic_clinical`
    as reachable behind the bar and lost sixteen variants. This test does not read
    the model. It runs the *real* selection stage on both sides of the partition,
    over every bar configuration crossed with every selection configuration, and
    compares the sets.

    The invariant is stated precisely, and the precision is the point. It is NOT
    "the bar is a superset of selection". When `selection.enabled` is false —
    which is the shipped demo's configuration, ADR 0027 §7 — selection retains
    literally everything, records with no gene symbol included, and the bar still
    drops those. The bar's licence to drop them does not come from selection at
    all; it comes from pairing being gene-scoped, so such a record forms no
    candidate pair under any configuration. The honest invariant is therefore
    about the PAIRABLE subset, and the exception is enumerated here rather than
    gestured at in a docstring.
    """
    label = f"{triage_label}/{selection_label}"
    with_bar = _run(
        tmp_path / "on", bar=True, triage_thresholds=triage, selection_thresholds=selection
    )
    without_bar = _run(tmp_path / "off", bar=False, selection_thresholds=selection)

    lost = set(without_bar.selected) - set(with_bar.selected)
    unpairable = {f.variant_id for f in FIXTURES if f.gene is None}
    unexplained = lost - unpairable
    assert not unexplained, (
        f"[{label}] the relevance bar removed {sorted(unexplained)} from the selected "
        "set, and selection would have kept them. narrow_to_selection is supposed to "
        "make that impossible (ADR 0028 decision 5)."
    )
    if lost:
        assert not selection.enabled, (
            f"[{label}] the bar dropped an unpairable record that selection retained "
            "while selection was ENABLED. With selection on, drop_without_gene_assignment "
            "should have made the two stages agree."
        )
    for variant_id in lost:
        assert explain_variant(with_bar.ledger, variant_id), (
            f"[{label}] {variant_id} was removed by the bar and is not in the ledger"
        )


def test_the_bars_residual_exposure_is_exactly_the_unpairable_class(tmp_path: Path) -> None:
    """Names the one gap the coupling leaves, so it cannot widen unnoticed.

    With `selection.enabled: false` — the shipped demo's configuration, ADR 0027
    §7 — selection retains every record and the bar still drops the ones with no
    gene assignment. This test pins that difference to exactly one variant, so
    that if the bar ever starts dropping something else in this configuration it
    fails here rather than in a byte-identity test on a different fixture.
    """
    disabled = SelectionThresholds(enabled=False)
    with_bar = _run(tmp_path / "on", bar=True, selection_thresholds=disabled)
    without_bar = _run(tmp_path / "off", bar=False, selection_thresholds=disabled)

    lost = set(without_bar.selected) - set(with_bar.selected)
    assert lost == {BY_KEY["no_consequence_at_all"].variant_id}, (
        "with selection disabled the bar may drop records with no gene assignment "
        "and nothing else. Pairing is gene-scoped, which is what makes that one "
        "provably answer-preserving; every other drop in this configuration would "
        "not be."
    )
    for variant_id in lost:
        assert explain_variant(with_bar.ledger, variant_id), (
            "the one class the bar is allowed to drop without selection behind it "
            "must still be accounted for by id (GP-13, ADR 0028)"
        )


def test_narrow_to_selection_does_not_close_the_clinical_axis() -> None:
    """The one thing the coupling does NOT do, pinned in a test rather than prose.

    `narrow_to_selection` opens with "Widen the bar until it cannot drop anything
    the selection stage would keep", and on every axis but one it does exactly
    that. On the clinical axis it does not, and it cannot: the rescue needs a
    ClinVar lookup and the lookup is the thing being skipped. The function's own
    closing paragraph says so.

    The gap is real and is left open deliberately (ADR 0028 decision 6 keeps the
    lookup in front of the bar by default and lets an operator move it back). What
    must not happen is that it becomes invisible — so this test asserts three
    things at once: the function leaves the setting alone, it reports no narrowing,
    and the warning that is the entire mitigation is present. If someone later
    teaches the function to force `clinical_before_bar` on, this test fails and
    tells them to delete the warning too.
    """
    from mva.annotation.triage import triage_warnings

    configured = TriageThresholds(clinical_before_bar=False)
    selection = SelectionThresholds(retain_pathogenic_clinical=True)
    effective, narrowings = narrow_to_selection(configured, selection)

    assert effective.clinical_before_bar is False, (
        "narrow_to_selection now forces the clinical lookup in front of the bar. "
        "That closes the gap this branch documents as its named residual risk — "
        "update ADR 0028 decision 6 and triage_warnings in the same change."
    )
    assert narrowings == (), (
        "a narrowing was reported for an axis the function does not actually narrow"
    )
    warnings = triage_warnings(effective, selection, narrowings)
    assert any("clinical_before_bar" in warning for warning in warnings), (
        "the clinical axis is the one axis the coupling leaves open, so the run "
        "warning is the whole mitigation and it has to fire"
    )


def test_narrowing_is_reported_whenever_the_bar_is_widened(tmp_path: Path) -> None:
    """A bar that silently widened is a bar whose configured value is a lie.

    Asserted on the configurations that provoke each narrowing, so the run report
    cannot quietly stop carrying one of them.
    """
    provoking = (
        SelectionThresholds(enabled=False),
        SelectionThresholds(drop_without_gene_assignment=False),
        SelectionThresholds(min_splice_ai_delta=0.05),
    )
    for selection in provoking:
        outcome = _run(tmp_path / str(id(selection)), bar=True, selection_thresholds=selection)
        assert outcome.narrowings, (
            f"narrow_to_selection widened the bar against {selection!r} without "
            "returning a narrowing for the run report to surface"
        )
    unprovoked = _run(tmp_path / "default", bar=True)
    assert unprovoked.narrowings == ()


# ---------------------------------------------------------------------------
# Fixture hygiene
# ---------------------------------------------------------------------------


def test_every_fixture_is_obviously_synthetic() -> None:
    """GP-20 / GP-40. No real gene symbol and no coordinate from any case volume."""
    for fixture in FIXTURES:
        for annotation in fixture.consequences:
            assert annotation.gene_symbol.startswith("SYNTH"), (
                f"[{fixture.key}] gene {annotation.gene_symbol!r} is not SYNTH-prefixed"
            )
            assert annotation.transcript_id.startswith("SYNTHT")
        for assertion in fixture.assertions:
            assert assertion.source.startswith("SYNTH")
            assert (assertion.accession or "").startswith("SYNTHVAR")


def test_every_rule_the_bar_can_fire_is_exercised_by_a_fixture(tmp_path: Path) -> None:
    """Coverage of the bar's own vocabulary, not of this file's line count.

    Every reason code in ``KEEP_REASONS`` and ``DROP_REASONS`` is a branch that
    can decide whether a variant is looked up. A reason no fixture reaches is a
    branch nothing tests, and the two that the DEFAULT bar cannot reach are named
    here with the configuration that reaches them, so "untested" and "unreachable
    by default" stay different statements.
    """
    from_defaults = {_triage_reason(f) for f in FIXTURES}
    assert from_defaults == {
        KEEP_RELEVANT_CONSEQUENCE,
        KEEP_IMPACT_AT_FLOOR,
        KEEP_SPLICE_AI,
        KEEP_IMPACT_NOT_ASSESSED,
        KEEP_CLINICAL_ASSERTION,
        DROP_IMPACT_BELOW_FLOOR,
        DROP_NO_GENE_ASSIGNMENT,
    }, f"the default bar reaches {sorted(from_defaults)} on this fixture set"

    # The two remaining reasons are unreachable under the default thresholds, so
    # each gets the configuration that reaches it rather than being left uncovered.
    unannotated = BY_KEY["no_consequence_at_all"]
    kept_no_gene = VariantTriage(TriageThresholds(retain_without_gene_assignment=True)).verdict(
        unannotated.consequences
    )
    assert kept_no_gene.reason == KEEP_NO_GENE_ASSIGNMENT

    off = _run(tmp_path, bar=False)
    assert set(off.annotated) == {f.variant_id for f in FIXTURES}, (
        "with the bar off every record must reach the annotated set, which is what "
        "makes it the control the byte-identity comparison is taken against"
    )
    assert len({f.position for f in FIXTURES}) == len(FIXTURES)


def test_a_never_assessed_drop_is_distinguishable_in_the_ledger(tmp_path: Path) -> None:
    """The one place the ledger's rule name is weaker than its value column.

    With ``retain_unassessed_impact`` false on BOTH stages, a variant nobody
    assessed is dropped by the bar under ``dropped_impact_below_floor`` — a rule
    name that asserts a comparison against the floor which never happened, because
    there was no impact to compare. Selection names the same event
    ``dropped_impact_not_assessed``, so the two stages spell one fact two ways.

    What rescues it is the value column: it reads ``not_assessed``, never
    ``modifier`` and never ``0`` (GP-14). This test pins that, because the value
    column is the only thing standing between the ledger and a row that says the
    annotator made a prediction it did not make.
    """
    triage = TriageThresholds(retain_unassessed_impact=False)
    selection = SelectionThresholds(retain_unassessed_impact=False)
    outcome = _run(tmp_path, bar=True, triage_thresholds=triage, selection_thresholds=selection)
    fixture = BY_KEY["impact_never_assessed"]

    (row,) = explain_variant(outcome.ledger, fixture.variant_id)
    assert row.stage is DiscardStage.TRIAGE
    assert row.value == "not_assessed", (
        "a variant nobody assessed was recorded with the impact value "
        f"{row.value!r}. NOT ASSESSED is not MODIFIER and the ledger may not say it "
        "is (ADR 0016, GP-14)."
    )
    assert row.rule == DROP_IMPACT_BELOW_FLOOR


# ---------------------------------------------------------------------------
# Drift, part two: the two stages' predicates over the whole consequence shape
# space, rather than over the ten shapes this file happens to name.
# ---------------------------------------------------------------------------

#: Impact values a transcript annotation can carry, NOT ASSESSED included.
_IMPACTS: tuple[ImpactSeverity | None, ...] = (
    None,
    ImpactSeverity.HIGH,
    ImpactSeverity.MODERATE,
    ImpactSeverity.LOW,
    ImpactSeverity.MODIFIER,
)

#: One term from each class the shared predicate distinguishes, plus a joined
#: form (SnpEff emits `a&b`) because the splitter is part of the predicate.
_TERMS: tuple[tuple[str, ...], ...] = (
    ("intron_variant",),
    ("synonymous_variant",),
    ("missense_variant",),
    ("splice_region_variant&synonymous_variant",),
    ("upstream_gene_variant",),
)

_DELTAS: tuple[float | None, ...] = (None, 0.0, 0.1, 0.2, 0.9)


def _shapes() -> list[tuple[ConsequenceAnnotation, ...]]:
    """Every two-transcript combination of impact x terms x SpliceAI delta.

    Two transcripts rather than one on purpose: every rule in both stages
    aggregates ACROSS transcripts, and a one-transcript grid cannot tell an
    aggregation bug from a per-annotation one. The mixed cases — assessed on one
    transcript, never assessed on the other — are the ones where the two stages
    could plausibly have read ADR 0016 differently.
    """
    singles = [
        _csq("SYNTHGRID1", terms, impact, splice_ai=delta)
        for impact in _IMPACTS
        for terms in _TERMS
        for delta in _DELTAS
    ]
    shapes: list[tuple[ConsequenceAnnotation, ...]] = [(one,) for one in singles]
    # A full cross-product of 125 x 125 would be 15,625 pairs for no extra
    # information: pair each annotation with a small fixed set of partners that
    # covers the aggregation boundaries.
    partners = (
        _csq("SYNTHGRID1", ("intron_variant",), None, transcript="SYNTHT0002.1"),
        _csq("SYNTHGRID1", ("intron_variant",), ImpactSeverity.MODIFIER, transcript="SYNTHT0002.1"),
        _csq(
            "SYNTHGRID1", ("missense_variant",), ImpactSeverity.MODERATE, transcript="SYNTHT0002.1"
        ),
        _csq(
            "SYNTHGRID1",
            ("intron_variant",),
            ImpactSeverity.MODIFIER,
            transcript="SYNTHT0002.1",
            splice_ai=0.9,
        ),
    )
    shapes.extend((one, partner) for one in singles for partner in partners)
    return shapes


@pytest.mark.parametrize("assertions", [(), (_pathogenic(),)], ids=["no_clinvar", "pathogenic"])
def test_no_consequence_shape_is_kept_by_selection_and_dropped_by_the_bar(
    assertions: tuple[ClinicalAssertion, ...],
) -> None:
    """The coupling over the input space, not over this file's ten fixtures.

    The configuration grid above varies the *knobs*; this varies the *data*. Both
    are needed, because the two stages do not share their implementations of
    "never assessed" or "SpliceAI above threshold" — only
    ``consequence_is_selecting`` is extracted. Each of the other two is written
    twice, once per stage, in slightly different shapes: the bar tracks the worst
    impact seen and asks ``worst is None``, selection asks
    ``all(annotation.impact is None ...)``. They agree today. Nothing but this
    test makes them keep agreeing.

    125 single-transcript shapes and 500 two-transcript ones, under the default
    thresholds, run once with a curated pathogenic assertion in play and once
    without — 1,250 comparisons of one stage's answer against the other's.
    """
    from mva.prioritization.selection import select_variants

    triage = VariantTriage(TriageThresholds())
    base = _record(FIXTURES[0], 0)
    shapes = _shapes()
    assert len(shapes) == 625, f"the shape grid is {len(shapes)}, not the 625 it describes"

    for consequences in shapes:
        verdict = triage.verdict(consequences, clinical_assertions=assertions)
        if verdict.retained:
            continue
        record = base.with_annotations(consequences=consequences, clinical_assertions=assertions)
        result = select_variants(
            [record],
            frequency=FrequencyThresholds(),
            thresholds=SelectionThresholds(),
            clock=_clock(),
        )
        (decision,) = result.decisions
        shape = [(a.consequence_terms, a.impact, a.splice_ai_delta_max) for a in consequences]
        assert not decision.retained, (
            "the relevance bar dropped a variant the selection stage retains.\n"
            f"  consequences: {shape}\n"
            f"  clinical:     {[a.significance for a in assertions]}\n"
            f"  bar said:     {verdict.reason} (value={verdict.value!r})\n"
            f"  selection:    {decision.reason}\n"
            "Two stages, two implementations of the same rule, and they have drifted "
            "(ADR 0028 decision 5)."
        )
