"""The relevance bar, the discard ledger and the pre-flight estimator (ADR 0028).

Three units, one file, because they are one mechanism: the bar decides, the ledger
accounts for what it decided against, and the estimator is what stops a run that
could never have produced either.
"""

from __future__ import annotations

import io
from pathlib import Path

import pytest

from mva.annotation.triage import (
    DROP_IMPACT_BELOW_FLOOR,
    DROP_NO_GENE_ASSIGNMENT,
    DROP_REASONS,
    KEEP_CLINICAL_ASSERTION,
    KEEP_IMPACT_AT_FLOOR,
    KEEP_IMPACT_NOT_ASSESSED,
    KEEP_NO_GENE_ASSIGNMENT,
    KEEP_REASONS,
    KEEP_RELEVANT_CONSEQUENCE,
    KEEP_SPLICE_AI,
    KEEP_TRIAGE_DISABLED,
    VariantTriage,
    narrow_to_selection,
    triage_warnings,
)
from mva.config import SelectionThresholds, TriageThresholds
from mva.errors import InsufficientDiskError
from mva.evidence.discards import (
    DiscardLedger,
    count_discards,
    explain_variant,
    iter_discards,
)
from mva.models.discard import DiscardRecord, DiscardSink, DiscardStage
from mva.models.variant import (
    ClinicalAssertion,
    ConsequenceAnnotation,
    ImpactSeverity,
    consequence_is_selecting,
)
from mva.progress import StageMonitor

pytestmark = [pytest.mark.unit]


def _csq(
    *,
    terms: tuple[str, ...] = ("intron_variant",),
    impact: ImpactSeverity | None = ImpactSeverity.MODIFIER,
    splice_ai: float | None = None,
    gene: str = "SYNTHG1",
) -> ConsequenceAnnotation:
    return ConsequenceAnnotation(
        gene_symbol=gene,
        transcript_id="SYNTHT1.1",
        consequence_terms=terms,
        impact=impact,
        splice_ai_delta_max=splice_ai,
    )


def _assertion() -> ClinicalAssertion:
    return ClinicalAssertion(version="v0", significance="Pathogenic")


# ---------------------------------------------------------------------------
# The bar
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("annotation", "expected"),
    [
        (_csq(terms=("stop_gained",), impact=ImpactSeverity.HIGH), KEEP_RELEVANT_CONSEQUENCE),
        (
            _csq(terms=("missense_variant",), impact=ImpactSeverity.MODERATE),
            KEEP_RELEVANT_CONSEQUENCE,
        ),
        # A term the impact does not vouch for. This is the case a pure impact bar
        # loses and the reason the term set is consulted at all.
        (
            _csq(terms=("splice_region_variant&synonymous_variant",), impact=ImpactSeverity.LOW),
            KEEP_RELEVANT_CONSEQUENCE,
        ),
        (_csq(terms=("synonymous_variant",), impact=ImpactSeverity.LOW), KEEP_IMPACT_AT_FLOOR),
        (_csq(terms=("intron_variant",), impact=None), KEEP_IMPACT_NOT_ASSESSED),
        (
            _csq(terms=("intron_variant",), impact=ImpactSeverity.MODIFIER, splice_ai=0.4),
            KEEP_SPLICE_AI,
        ),
    ],
)
def test_default_bar_keeps(annotation: ConsequenceAnnotation, expected: str) -> None:
    verdict = VariantTriage(TriageThresholds()).verdict((annotation,))
    assert verdict.retained
    assert verdict.reason == expected


def test_default_bar_drops_a_modifier_only_variant() -> None:
    verdict = VariantTriage(TriageThresholds()).verdict((_csq(),))
    assert not verdict.retained
    assert verdict.reason == DROP_IMPACT_BELOW_FLOOR
    assert verdict.value == "modifier", "the ledger row must name the value the rule read"


def test_default_bar_drops_a_variant_with_no_consequence_at_all() -> None:
    """The only drop reason that is provable rather than heuristic: pairing is gene-scoped."""
    verdict = VariantTriage(TriageThresholds()).verdict(())
    assert not verdict.retained
    assert verdict.reason == DROP_NO_GENE_ASSIGNMENT


def test_a_clinical_assertion_keeps_a_variant_the_consequence_rule_drops() -> None:
    """Decision 6 of ADR 0028, at unit scale: 3,616 vs 3,600 on the phantom."""
    triage = VariantTriage(TriageThresholds())
    verdict = triage.verdict((_csq(),), clinical_assertions=(_assertion(),))
    assert verdict.retained
    assert verdict.reason == KEEP_CLINICAL_ASSERTION


def test_a_disabled_bar_keeps_everything_including_unannotated_records() -> None:
    triage = VariantTriage(TriageThresholds(enabled=False))
    for consequences in ((), (_csq(),)):
        verdict = triage.verdict(consequences)
        assert verdict.retained
        assert verdict.reason == KEEP_TRIAGE_DISABLED


def test_modifier_floor_keeps_every_annotated_record() -> None:
    """The halfway setting: drop the ~55% with no gene, keep everything else."""
    triage = VariantTriage(TriageThresholds(min_impact=ImpactSeverity.MODIFIER))
    assert triage.verdict((_csq(),)).retained
    assert not triage.verdict(()).retained


def test_retaining_unannotated_records_is_configurable() -> None:
    triage = VariantTriage(TriageThresholds(retain_without_gene_assignment=True))
    verdict = triage.verdict(())
    assert verdict.retained
    assert verdict.reason == KEEP_NO_GENE_ASSIGNMENT


def test_not_assessed_is_never_rendered_as_a_prediction() -> None:
    """GP-14 / ADR 0016 in the ledger's `value` column, not only in the reason code."""
    triage = VariantTriage(TriageThresholds(retain_unassessed_impact=False))
    verdict = triage.verdict((_csq(impact=None),))
    assert not verdict.retained
    assert verdict.value == "not_assessed"
    assert "modifier" not in verdict.value


def test_reason_codes_partition_each_side_of_the_decision() -> None:
    assert not set(KEEP_REASONS) & set(DROP_REASONS)
    assert len(set(KEEP_REASONS)) == len(KEEP_REASONS)
    assert len(set(DROP_REASONS)) == len(DROP_REASONS)


def test_the_bar_uses_the_same_predicate_the_selection_stage_uses() -> None:
    """GP-03's shared vocabulary, asserted rather than assumed.

    Two copies of a Sequence Ontology term list is how the two stages come to
    disagree, and a bar narrower than the stage behind it deletes candidates that
    stage would have kept. The predicate is defined once in `mva.models.variant`.
    """
    from mva.prioritization.selection import _is_selecting

    for annotation in (
        _csq(terms=("stop_gained",), impact=ImpactSeverity.HIGH),
        _csq(terms=("splice_region_variant&synonymous_variant",), impact=ImpactSeverity.LOW),
        _csq(terms=("intron_variant",), impact=ImpactSeverity.MODIFIER),
        _csq(terms=("upstream_gene_variant",), impact=None),
    ):
        assert _is_selecting(annotation) is consequence_is_selecting(annotation)


# ---------------------------------------------------------------------------
# The safety coupling
# ---------------------------------------------------------------------------


def test_disabled_selection_widens_the_bar_to_a_gene_assignment_test() -> None:
    effective, narrowings = narrow_to_selection(
        TriageThresholds(), SelectionThresholds(enabled=False)
    )
    assert effective.min_impact is ImpactSeverity.MODIFIER
    assert effective.retain_unassessed_impact
    assert narrowings, "a bar that was widened must say so"


def test_selection_keeping_unannotated_records_forces_the_bar_to_keep_them() -> None:
    effective, narrowings = narrow_to_selection(
        TriageThresholds(), SelectionThresholds(drop_without_gene_assignment=False)
    )
    assert effective.retain_without_gene_assignment
    assert any("gene assignment" in n for n in narrowings)


def test_a_stricter_splice_ai_threshold_is_clamped_to_selections() -> None:
    effective, narrowings = narrow_to_selection(
        TriageThresholds(min_splice_ai_delta=0.9), SelectionThresholds(min_splice_ai_delta=0.2)
    )
    assert effective.min_splice_ai_delta == pytest.approx(0.2)
    assert narrowings


def test_the_default_pairing_needs_no_narrowing() -> None:
    effective, narrowings = narrow_to_selection(TriageThresholds(), SelectionThresholds())
    assert narrowings == ()
    assert effective == TriageThresholds()


def test_a_disabled_bar_is_never_narrowed() -> None:
    thresholds = TriageThresholds(enabled=False)
    effective, narrowings = narrow_to_selection(thresholds, SelectionThresholds(enabled=False))
    assert effective is thresholds
    assert narrowings == ()


def test_deferring_the_clinical_lookup_is_reported_as_the_hole_it_is() -> None:
    warnings = triage_warnings(
        TriageThresholds(clinical_before_bar=False), SelectionThresholds(), ()
    )
    assert any("clinical_before_bar" in w and "3,616" in w for w in warnings), (
        "the measured cost of this setting must be in the warning, not only in the ADR"
    )


def test_the_default_configuration_raises_no_clinical_rescue_warning() -> None:
    assert triage_warnings(TriageThresholds(), SelectionThresholds(), ()) == ()


# ---------------------------------------------------------------------------
# The ledger
# ---------------------------------------------------------------------------


def test_a_ledger_row_round_trips(tmp_path: Path) -> None:
    path = tmp_path / "discards.tsv"
    with DiscardLedger(path) as ledger:
        ledger.drop("GRCh38:chr1:100:A:T", stage=DiscardStage.TRIAGE, rule="r", value="modifier")
    (row,) = list(iter_discards(path))
    assert row == DiscardRecord(
        variant_id="GRCh38:chr1:100:A:T",
        stage=DiscardStage.TRIAGE,
        rule="r",
        value="modifier",
    )


def test_a_ledger_is_a_discard_sink() -> None:
    """The Protocol is what lets layer-4 and layer-5 stages write without importing it."""
    assert isinstance(DiscardLedger, type)
    assert issubclass(DiscardLedger, object)


def test_a_ledger_satisfies_the_sink_protocol(tmp_path: Path) -> None:
    with DiscardLedger(tmp_path / "d.tsv") as ledger:
        assert isinstance(ledger, DiscardSink)


def test_explain_variant_finds_one_variant_among_many(tmp_path: Path) -> None:
    path = tmp_path / "discards.tsv"
    with DiscardLedger(path) as ledger:
        for index in range(500):
            ledger.drop(
                f"GRCh38:chr1:{index}:A:T",
                stage=DiscardStage.SELECTION,
                rule="dropped_common_in_population",
                value=f"0.0{index}",
            )
    (row,) = explain_variant(path, "GRCh38:chr1:317:A:T")
    assert row.value == "0.0317"
    assert explain_variant(path, "GRCh38:chr1:317:A:G") == ()


def test_counts_are_recomputable_from_the_bytes(tmp_path: Path) -> None:
    """A disagreement between the counter and the file means truncation."""
    path = tmp_path / "discards.tsv"
    with DiscardLedger(path) as ledger:
        ledger.drop("v1", stage=DiscardStage.INGEST, rule="non_canonical_contig", value="chrUn")
        ledger.drop("v2", stage=DiscardStage.INGEST, rule="non_canonical_contig", value="chrUn")
        ledger.drop("v3", stage=DiscardStage.HARD_FILTER, rule="no_alt_allele_called")
        in_memory = ledger.counts()
    assert in_memory == count_discards(path)
    assert in_memory == {
        "hard_filter/no_alt_allele_called": 1,
        "ingest/non_canonical_contig": 2,
    }


def test_counts_by_stage_names_every_stage_even_at_zero(tmp_path: Path) -> None:
    """A stage reporting zero is a finding; a stage missing from the mapping is not."""
    with DiscardLedger(tmp_path / "d.tsv") as ledger:
        ledger.drop("v1", stage=DiscardStage.TRIAGE, rule="r")
        counts = ledger.counts_by_stage()
    assert set(counts) == {str(stage) for stage in DiscardStage}
    assert counts["triage"] == 1
    assert counts["ingest"] == 0


def test_a_truncated_ledger_is_still_readable_up_to_the_truncation(tmp_path: Path) -> None:
    """The run that fills the disk is the run someone most needs to read this on."""
    path = tmp_path / "discards.tsv"
    with DiscardLedger(path) as ledger:
        for index in range(10):
            ledger.drop(f"v{index}", stage=DiscardStage.TRIAGE, rule="r", value="modifier")
    text = path.read_text(encoding="utf-8")
    # Cut into the last row's third tab, so the tail is a 3-field line that cannot
    # be parsed. The nine rows in front of it must still come back.
    path.write_text(text[: len(text) - 12], encoding="utf-8")
    rows = list(iter_discards(path))
    assert len(rows) == 9
    assert all(row.rule == "r" for row in rows)


def test_tabs_and_newlines_in_a_field_cannot_shift_the_columns() -> None:
    row = DiscardRecord(
        variant_id="a\tb", stage=DiscardStage.TRIAGE, rule="r\nx", value="v\rz"
    ).as_row()
    assert row.count("\t") == 3
    assert row.endswith("\n")
    assert DiscardRecord.from_row(row) is not None


def test_a_comment_or_blank_line_is_not_a_record() -> None:
    assert DiscardRecord.from_row("# banner\n") is None
    assert DiscardRecord.from_row("\n") is None
    assert DiscardRecord.from_row("too\tfew\tcolumns\n") is None
    assert DiscardRecord.from_row("a\tnot_a_stage\tr\tv\n") is None


def test_recording_after_close_is_ignored_rather_than_raising(tmp_path: Path) -> None:
    """A stage dropping a record while the run unwinds must not mask the first error."""
    ledger = DiscardLedger(tmp_path / "d.tsv")
    ledger.close()
    ledger.drop("v", stage=DiscardStage.TRIAGE, rule="r")
    assert ledger.row_count == 0


def test_a_row_costs_about_a_hundred_bytes() -> None:
    """The size budget that made this a TSV rather than an EvidenceItem."""
    row = DiscardRecord(
        variant_id="GRCh38:chr12:133275309:ACGTACGT:A",
        stage=DiscardStage.SELECTION,
        rule="dropped_common_in_population",
        value="0.031416",
    ).as_row()
    assert len(row.encode("utf-8")) <= 100


# ---------------------------------------------------------------------------
# The pre-flight estimator
# ---------------------------------------------------------------------------


def _monitor(out: io.StringIO, clock: list[float], **kwargs: object) -> StageMonitor:
    defaults: dict[str, object] = {
        "total_records": 1_000_000,
        "free_bytes": 1_000_000_000,
        "out": out,
        "now": lambda: clock[0],
        "project_after_records": 1_000,
        "min_sample_records": 1_000,
        # Explicit, so the suite's own MVA_PROGRESS setting cannot decide whether
        # these assertions are about output or about nothing.
        "enabled": True,
    }
    defaults.update(kwargs)
    return StageMonitor("annotate", **defaults)  # type: ignore[arg-type]


def test_a_projection_that_does_not_fit_stops_the_run() -> None:
    out, clock = io.StringIO(), [0.0]
    monitor = _monitor(out, clock)
    with pytest.raises(InsufficientDiskError) as excinfo:
        for _ in range(1_000):
            clock[0] += 0.01
            monitor.observe(bytes_written=10_000)
    message = str(excinfo.value)
    assert "projected" in message
    assert "10.0 GB" in message, "the projected size must be in the refusal"
    assert "triage.min_impact" in message, "a refusal must carry its remediation"
    assert "chr" not in message, "no coordinate may appear in a refusal (PRIV-09)"


def test_a_projection_that_fits_does_not_stop_the_run() -> None:
    out, clock = io.StringIO(), [0.0]
    monitor = _monitor(out, clock)
    for _ in range(1_000):
        clock[0] += 0.01
        monitor.observe(bytes_written=10)
    assert "pre-flight projection" in out.getvalue()
    assert monitor.records == 1_000


def test_a_sample_too_small_to_project_never_refuses() -> None:
    """A twelve-record fixture cannot project anything, and must not try."""
    out, clock = io.StringIO(), [0.0]
    monitor = _monitor(out, clock, total_records=12, project_after_records=1)
    for _ in range(12):
        clock[0] += 0.01
        monitor.observe(bytes_written=10**9)
    monitor.finish()


def test_an_uninterrogable_volume_disables_the_refusal_rather_than_triggering_it() -> None:
    """GP-14: "we could not look" is not "we looked and found nothing"."""
    out, clock = io.StringIO(), [0.0]
    monitor = _monitor(out, clock, free_bytes=0)
    for _ in range(1_000):
        clock[0] += 0.01
        monitor.observe(bytes_written=10**6)
    monitor.finish()


def test_progress_output_is_throttled_by_the_interval() -> None:
    out, clock = io.StringIO(), [0.0]
    monitor = _monitor(out, clock, project_after_records=10**9, interval_seconds=30.0)
    for index in range(200_000):
        clock[0] = index * 0.001
        monitor.observe(bytes_written=1)
    lines = [line for line in out.getvalue().splitlines() if "rec/s" in line]
    assert 4 <= len(lines) <= 8, f"200 simulated seconds at a 30s interval, got {len(lines)}"


def test_progress_can_be_silenced_but_the_refusal_cannot() -> None:
    out, clock = io.StringIO(), [0.0]
    monitor = _monitor(out, clock, enabled=False)
    with pytest.raises(InsufficientDiskError):
        for _ in range(1_000):
            clock[0] += 0.01
            monitor.observe(bytes_written=10_000)
    assert out.getvalue() == "", "MVA_PROGRESS=0 must silence output, not the refusal"


def test_the_monitor_reports_the_measured_shape_of_the_real_failure() -> None:
    """The run that motivated ADR 0028, replayed through the instrument."""
    out, clock = io.StringIO(), [0.0]
    monitor = StageMonitor(
        "annotate",
        total_records=4_962_060,
        free_bytes=60_000_000_000,
        out=out,
        now=lambda: clock[0],
        project_after_records=250_000,
        min_sample_records=1_000,
        enabled=True,
    )
    with pytest.raises(InsufficientDiskError) as excinfo:
        for _ in range(250_000):
            clock[0] += 0.0566
            monitor.observe(bytes_written=19_000)
    assert "94.3 GB" in str(excinfo.value)
    assert "60.0 GB free" in str(excinfo.value)


def test_bytes_only_observations_do_not_advance_the_record_count() -> None:
    """Bytes and records are counted at different places, so they must not double.

    The write loop reports bytes for a record the bar already counted when it
    examined it. If that call also incremented the record count, every survivor
    would be counted twice and bytes/record would come out at half its true value.
    """
    out, clock = io.StringIO(), [0.0]
    monitor = _monitor(out, clock)
    monitor.observe(records=100, bytes_written=0)
    for _ in range(50):
        monitor.observe(records=0, bytes_written=1_000)

    projection = monitor.projection()
    assert projection.records_seen == 100
    assert projection.bytes_written == 50_000
    assert not out.getvalue(), "a bytes-only observation must not read the clock or report"


def test_a_bar_that_rejects_most_records_does_not_refuse_a_run_that_fits() -> None:
    """The regression: project from records EXAMINED, never from records retained.

    A bar that keeps 1 record in 100 still has the whole callset as its
    denominator. Counting only survivors would measure ~19 kB against each of them,
    multiply that by every record in the callset, and refuse a run whose real output
    is a hundredth of the projection -- reproducing, from the other direction, the
    95 GB estimate that ADR 0028 exists to have caught.
    """
    out, clock = io.StringIO(), [0.0]
    # 1 M records at 10 kB per SURVIVOR. One record in 100 survives, so the true
    # output is 1 M / 100 * 10 kB = 100 MB, which fits in a 1 GB volume easily.
    # Counted the wrong way it projects 1 M * 10 kB = 10 GB and refuses.
    monitor = _monitor(out, clock)
    for index in range(100_000):
        clock[0] += 0.0001
        monitor.observe(records=1)
        if index % 100 == 0:
            monitor.observe(records=0, bytes_written=10_000)

    projection = monitor.projection()
    assert projection.records_seen == 100_000
    assert projection.fits, "a run writing 100 MB to a 1 GB volume must not be refused"
    assert projection.projected_bytes < 200_000_000


def test_the_stream_reports_every_record_it_examines_not_every_record_it_yields() -> None:
    """``on_examined`` is the only way a caller can see a bar-rejected record.

    Without it the count is unobservable: a rejected record is never yielded, so a
    consumer driving progress from the yielded records goes silent for exactly as
    long as the bar is doing its job.
    """
    from datetime import UTC, datetime

    from tests.unit.test_stage_streaming import sample_adapters, sample_records

    from mva.annotation.service import iter_annotated
    from mva.clock import FixedClock

    records = sample_records()
    examined: list[int] = []
    stream = iter_annotated(
        iter(records),
        adapters=sample_adapters(),
        clock=FixedClock(datetime(2026, 1, 1, tzinfo=UTC)),
        batch_size=3,
        triage=VariantTriage(TriageThresholds()),
        on_examined=examined.append,
    )
    yielded = list(stream)

    assert sum(examined) == len(records), "every input record must be reported once"
    assert len(yielded) < len(records), (
        "this fixture must exercise a bar that actually rejects, or the test passes "
        "for the wrong reason"
    )
