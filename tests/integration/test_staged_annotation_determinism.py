"""The staged annotation path does not change the answer (ADR 0028, GP-30).

This is the test the whole partition rests on. Everything else about ADR 0028 is
a throughput argument, and a throughput argument that moves a rank is worthless.

The claim under test is narrow and total: **turning the relevance bar on changes
which variants are annotated in full, and changes nothing that reaches the
output.** Not "the same gene ranks first"; not "the scores are close". The
submission CSV and the ranked-pairs artifact must be byte-identical, because
byte-identity is the only version of this claim that cannot be argued with.

Two configurations are compared, and the second one is the one that matters:

1. The shipped synthetic case, which runs with ``selection.enabled: false``
   (ADR 0027 §7). The bar narrows itself to a gene-assignment test here — see
   :func:`mva.annotation.triage.narrow_to_selection` — so this asserts the
   narrowing works, not that the bar is selective.
2. The same case with selection ENABLED, which is what a real case runs. Here the
   bar is at its configured width and does drop records, so this is the
   configuration in which a partition bug would show up.

A third test asserts the bar is not vacuous: if it dropped nothing in
configuration 2, the first two tests would pass for the wrong reason.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from mva.config import CaseConfig, Workspace
from mva.evidence.discards import count_discards, explain_variant, iter_discards
from mva.models.discard import DiscardStage
from mva.orchestrator import execute_pipeline

pytestmark = [pytest.mark.integration]

#: Artifacts that must not move by one byte when the bar is switched on. The
#: submission is the deliverable; the ranked pairs are what it is composed from,
#: and including it turns "the top ten agree" into "the whole ranking agrees".
INVARIANT_ARTIFACTS: tuple[str, ...] = (
    "submission/track1_submission.csv",
    "candidates/ranked_pairs.json",
)


def _workspace(root: Path, fixtures: Path) -> Workspace:
    inputs = root / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    shutil.copy(fixtures / "synthetic_case.vcf", inputs / "synthetic_case.vcf")
    shutil.copy(fixtures / "synthetic_phenotype.tsv", inputs / "synthetic_phenotype.tsv")
    return Workspace(root=root.resolve(), repo_root=fixtures.parents[2])


@pytest.fixture
def fixtures(repo_root: Path) -> Path:
    return repo_root / "tests" / "fixtures" / "synthetic"


def _run(config: CaseConfig, root: Path, fixtures: Path) -> tuple[Path, dict[str, bytes]]:
    """Execute the pipeline into a fresh workspace and return the invariant bytes."""
    workspace = _workspace(root, fixtures)
    result = execute_pipeline(config, workspace)
    run_dir = workspace.runs_dir / result.run_id
    return run_dir, {
        relative: (run_dir / relative).read_bytes() for relative in INVARIANT_ARTIFACTS
    }


def _with(config: CaseConfig, **triage: object) -> CaseConfig:
    return config.model_copy(update={"triage": config.triage.model_copy(update=triage)})


def _selection_on(config: CaseConfig) -> CaseConfig:
    enabled = config.selection.model_copy(update={"enabled": True})
    return config.model_copy(update={"selection": enabled})


# ---------------------------------------------------------------------------
# The proof
# ---------------------------------------------------------------------------


def test_bar_on_and_off_produce_identical_output_on_the_shipped_case(
    synthetic_config: CaseConfig, tmp_path: Path, fixtures: Path
) -> None:
    """The demo case: selection is off, so the bar narrows to a gene-assignment test."""
    _, without = _run(_with(synthetic_config, enabled=False), tmp_path / "off", fixtures)
    _, with_bar = _run(_with(synthetic_config, enabled=True), tmp_path / "on", fixtures)
    for relative in INVARIANT_ARTIFACTS:
        assert with_bar[relative] == without[relative], (
            f"{relative} differs with the annotation relevance bar on. The partition "
            "changed the answer, which is the one thing ADR 0028 may not do (GP-30)."
        )


def test_bar_on_and_off_produce_identical_output_with_selection_enabled(
    synthetic_config: CaseConfig, tmp_path: Path, fixtures: Path
) -> None:
    """The real-case configuration, where the bar runs at its configured width."""
    config = _selection_on(synthetic_config)
    _, without = _run(_with(config, enabled=False), tmp_path / "off", fixtures)
    _, with_bar = _run(_with(config, enabled=True), tmp_path / "on", fixtures)
    for relative in INVARIANT_ARTIFACTS:
        assert with_bar[relative] == without[relative], (
            f"{relative} differs with the annotation relevance bar on, under "
            "selection.enabled=true. A variant the bar dropped was one selection "
            "would have kept — see triage.narrow_to_selection (ADR 0028, GP-30)."
        )


def test_the_bar_is_not_vacuous_with_selection_enabled(
    synthetic_config: CaseConfig, tmp_path: Path, fixtures: Path
) -> None:
    """A bar that drops nothing would make the two tests above pass for free.

    The synthetic fixture holds one MODIFIER-only record — an intronic call in
    SYNTHSOL5 — precisely so that this assertion has something to find. If the
    fixture changes and this fails, the fix is a new MODIFIER-only record, not a
    weaker assertion.
    """
    config = _selection_on(synthetic_config)
    run_dir, _ = _run(_with(config, enabled=True), tmp_path / "on", fixtures)
    counts = count_discards(run_dir / "discards" / "discard_ledger.tsv")
    prefix = f"{DiscardStage.TRIAGE}/"
    triage_drops = {key: n for key, n in counts.items() if key.startswith(prefix)}
    assert triage_drops, (
        "the relevance bar dropped nothing on the synthetic case, so the "
        "byte-identity tests above prove nothing about the partition"
    )


# ---------------------------------------------------------------------------
# The ledger's own contract
# ---------------------------------------------------------------------------


def test_every_dropped_variant_is_named_in_the_ledger(
    synthetic_config: CaseConfig, tmp_path: Path, fixtures: Path
) -> None:
    """Nothing dropped becomes invisible: the amended GP-13 in one assertion.

    Cross-checked against the *other* run rather than against the pipeline's own
    counters, which is the point — a variant present in the annotated artifact
    when the bar is off and absent when it is on must appear in the ledger by id.
    A counter agreeing with itself would prove nothing.
    """
    import json

    config = _selection_on(synthetic_config)
    off_dir, _ = _run(_with(config, enabled=False), tmp_path / "off", fixtures)
    on_dir, _ = _run(_with(config, enabled=True), tmp_path / "on", fixtures)

    def annotated_ids(run_dir: Path) -> set[str]:
        rows = json.loads((run_dir / "variants" / "annotated.json").read_text(encoding="utf-8"))
        return {
            f"{r['coordinate']['build']}:{r['coordinate']['contig']}:"
            f"{r['coordinate']['position']}:{r['coordinate']['ref']}:{r['coordinate']['alt']}"
            for r in rows
        }

    lost = annotated_ids(off_dir) - annotated_ids(on_dir)
    assert lost, "the bar removed nothing from the annotated artifact"
    ledger = on_dir / "discards" / "discard_ledger.tsv"
    for variant_id in sorted(lost):
        rows = explain_variant(ledger, variant_id)
        assert rows, (
            f"{variant_id} vanished from variants/annotated.json when the relevance "
            "bar was switched on and is NOT in the discard ledger. A dropped candidate "
            "that is not accounted for is invisible and unfalsifiable (GP-13, ADR 0028)."
        )
        assert all(row.rule for row in rows), "a ledger row carries no rule"


def test_ledger_rows_are_compact_enough_to_hold_a_whole_genome(
    synthetic_config: CaseConfig, tmp_path: Path, fixtures: Path
) -> None:
    """~100 bytes a row, so 5 M rows cost ~500 MB rather than ~65 GB.

    The size budget is the reason this is a TSV and not an ``EvidenceItem``, and a
    budget nothing measures is a wish. Asserted on the mean over real rows rather
    than on a constructed one, so a field added later that happens to be long
    fails here rather than on the run that needed the space.
    """
    config = _selection_on(synthetic_config)
    run_dir, _ = _run(_with(config, enabled=True), tmp_path / "on", fixtures)
    ledger = run_dir / "discards" / "discard_ledger.tsv"
    rows = list(iter_discards(ledger))
    assert rows, "no discard rows to measure"
    mean = sum(len(row.as_row().encode("utf-8")) for row in rows) / len(rows)
    assert mean <= 128, (
        f"discard rows average {mean:.0f} bytes. At 5 M rows that is "
        f"{mean * 5_000_000 / 1e9:.1f} GB, and the whole point of not using "
        "EvidenceItems was that 5 M of those cost ~65 GB (ADR 0028)."
    )


def test_ledger_holds_no_row_for_a_variant_that_survived(
    synthetic_config: CaseConfig, tmp_path: Path, fixtures: Path
) -> None:
    """The answer's own variants are never in the ledger, under any configuration."""
    config = _selection_on(synthetic_config)
    workspace = _workspace(tmp_path / "on", fixtures)
    result = execute_pipeline(_with(config, enabled=True), workspace)
    assert result.ranked_pairs
    ledger = workspace.runs_dir / result.run_id / "discards" / "discard_ledger.tsv"
    for variant_id in result.ranked_pairs[0].variant_ids:
        assert not explain_variant(ledger, variant_id), (
            f"{variant_id} is in the top-ranked pair AND in the discard ledger"
        )
