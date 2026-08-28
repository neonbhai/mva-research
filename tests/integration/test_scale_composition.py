"""The whole-genome composition, proved equivalent to the one it replaces.

``docs/handoff-scale.md`` asks the composition root to swap a five-line batch
sequence for a single streamed pass. That is only a safe trade if the pass emits
the *same artifact* and reaches pairing with the *same records*. This test runs
both compositions over the synthetic fixture with the real local adapters and
compares them byte for byte, so the handoff is a checked claim rather than a
described one.

It also pins what selection does to the demo case, because that is the number the
composition root will be surprised by: 18 candidates become 9, and the golden
rank-1 pair survives untouched.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import pytest

from mva.annotation.local_tables import load_default_adapters
from mva.annotation.service import annotate_variants, iter_annotated
from mva.clock import demo_clock
from mva.config import CaseConfig, FrequencyThresholds, QualityThresholds, Workspace
from mva.evidence.ledger import EvidenceLedger
from mva.ingestion.normalise import normalise_variants
from mva.ingestion.qc import assess_quality
from mva.ingestion.reader import read_vcf
from mva.models.genome import GenomeBuild
from mva.models.provenance import ArtifactKind
from mva.models.variant import VariantRecord
from mva.pipeline import validate_case
from mva.prioritization.filters import (
    apply_hard_filters,
    apply_soft_flags,
    iter_hard_filtered,
    iter_soft_flagged,
)
from mva.prioritization.pairing import generate_pair_candidates
from mva.prioritization.selection import SelectionThresholds, iter_selected

REPO_ROOT = Path(__file__).resolve().parents[2]
KNOWLEDGE_ROOT = REPO_ROOT / "knowledge"
MANIFEST_PATH = KNOWLEDGE_ROOT / "manifests" / "knowledge.yaml"

pytestmark = pytest.mark.integration


def _qc_variants(workspace: Workspace, config: CaseConfig) -> tuple[VariantRecord, ...]:
    ingestion = read_vcf(
        workspace.path(config.inputs.vcf),
        expected_build=config.genome_build,
        source_artifact="input_vcf",
    )
    normalised = normalise_variants(ingestion.variants)
    return assess_quality(
        normalised.variants, thresholds=config.quality, clock=demo_clock()
    ).variants


def test_the_streamed_composition_writes_the_same_annotated_artifact(
    synthetic_config: CaseConfig, synthetic_workspace: Workspace
) -> None:
    """Same bytes, same evidence, same records into pairing — one pass instead of five."""
    context = validate_case(synthetic_config, synthetic_workspace, clock=demo_clock())
    variants = _qc_variants(synthetic_workspace, synthetic_config)
    adapters = load_default_adapters(KNOWLEDGE_ROOT, MANIFEST_PATH)

    # ---- the composition as orchestrator.py has it today
    batch = annotate_variants(variants, adapters=adapters, clock=context.clock)
    before_ledger = EvidenceLedger(run_id=context.run_id)
    before_ledger.extend(batch.evidence)
    before_art = context.write_json_artifact(
        "before/annotated.json",
        [v.model_dump(mode="json") for v in batch.variants],
        kind=ArtifactKind.ANNOTATED_VARIANTS,
        stage="annotate",
        row_count=len(batch.variants),
    )
    before_filtered = apply_hard_filters(
        batch.variants, expected_build=synthetic_config.genome_build
    )
    before_flagged = apply_soft_flags(
        before_filtered.retained,
        frequency=synthetic_config.frequency,
        quality=synthetic_config.quality,
    )

    # ---- the composition docs/handoff-scale.md asks for
    after_ledger = EvidenceLedger(run_id=context.run_id)
    annotated = iter_annotated(iter(variants), adapters=adapters, clock=context.clock)
    with context.open_json_rows_artifact(
        "after/annotated.json",
        kind=ArtifactKind.ANNOTATED_VARIANTS,
        stage="annotate",
    ) as sink:

        def _recorded() -> Iterator[VariantRecord]:
            for item in annotated:
                after_ledger.extend(item.evidence)
                sink.write(item.variant.model_dump(mode="json"))
                yield item.variant

        after_filtered = iter_hard_filtered(
            _recorded(), expected_build=synthetic_config.genome_build
        )
        after_flagged = iter_soft_flagged(
            after_filtered,
            frequency=synthetic_config.frequency,
            quality=synthetic_config.quality,
        )
        selection = iter_selected(
            after_flagged,
            frequency=synthetic_config.frequency,
            thresholds=SelectionThresholds(enabled=False),
            clock=context.clock,
        )
        selected = list(selection)

    after_art = sink.provenance
    assert after_art is not None

    root = synthetic_workspace.root
    assert (root / after_art.relative_path).read_bytes() == (
        root / before_art.relative_path
    ).read_bytes()
    assert after_art.content_hash == before_art.content_hash
    assert after_art.row_count == before_art.row_count

    assert after_ledger.items() == before_ledger.items()
    assert annotated.coverage() == batch.coverage
    assert annotated.warnings() == batch.warnings
    assert after_filtered.counts() == before_filtered.counts
    # Selection disabled, so the streamed path reaches pairing with exactly the
    # records the batch path does.
    assert tuple(selected) == before_flagged


def test_selection_shrinks_the_demo_candidate_set_without_moving_rank_one(
    synthetic_config: CaseConfig, synthetic_workspace: Workspace
) -> None:
    """The before/after GP-32 asks for, on the case the golden files lock.

    Selection removes four of twelve variants and halves the candidate list. The
    golden rank-1 pair is built from two variants that both survive, so
    `tests/golden/expected_ranking.tsv` needs no re-baselining — which is the
    thing to check before wiring this in, not after.
    """
    variants = _qc_variants(synthetic_workspace, synthetic_config)
    adapters = load_default_adapters(KNOWLEDGE_ROOT, MANIFEST_PATH)
    annotated = annotate_variants(variants, adapters=adapters, clock=demo_clock()).variants
    flagged = apply_soft_flags(
        apply_hard_filters(annotated, expected_build=GenomeBuild.GRCH38).retained,
        frequency=FrequencyThresholds(),
        quality=QualityThresholds(),
    )

    before = generate_pair_candidates(
        flagged,
        max_pairs_per_gene=synthetic_config.max_pairs_per_gene,
        frequency=synthetic_config.frequency,
    )
    stream = iter_selected(iter(flagged), frequency=FrequencyThresholds(), clock=demo_clock())
    selected = list(stream)
    report = stream.report()
    after = generate_pair_candidates(
        selected,
        max_pairs_per_gene=synthetic_config.max_pairs_per_gene,
        frequency=synthetic_config.frequency,
    )

    assert report.input_count == 12
    assert report.retained_count == 8
    assert report.dropped_by_reason["dropped_common_in_population"] == 3
    assert report.dropped_by_reason["dropped_not_coding_or_splice"] == 1
    # One variant has no frequency data at all and is RETAINED (GP-14).
    assert report.notes["frequency_unknown"] == 1

    assert len(before.candidates) == 18
    assert len(after.candidates) == 9

    golden = ("GRCh38:chr15:40200000:C:T", "GRCh38:chr15:40210500:G:A")
    survivors = {v.variant_id for v in selected}
    assert set(golden) <= survivors, "the golden rank-1 pair must survive selection"
    assert any(tuple(candidate.variant_ids) == golden for candidate in after.candidates), (
        "the golden rank-1 candidate must still be enumerated"
    )
