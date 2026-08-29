"""The clinical rescue, exercised through the real pipeline (ADR 0028 decision 6).

`test_staged_annotation_determinism.py` runs the shipped synthetic case end to
end and compares the submission byte for byte. It cannot see the clinical rescue,
because the case binds :class:`~mva.annotation.local_tables.NullClinicalAdapter`
— an adapter that returns ``{}`` for every id. Every branch that depends on a
curated assertion is dead code for the whole of that test:
``AnnotationStream._triage_batch`` never receives one, ``VariantTriage``'s
clinical rule never fires, and ``selection.retain_pathogenic_clinical`` never
rescues anything.

That is not a hypothetical blind spot. ADR 0028's own decision 6 records that the
first implementation deferred the ClinVar lookup behind the bar and lost sixteen
variants on the scale phantom — and the byte-identity tests stayed green for the
entire time that bug was live, because with no assertions in play there was
nothing to lose.

This module binds a clinical adapter that answers, on the case's one MODIFIER-only
variant, and asks the three questions the shipped fixture cannot:

1. does the partition still preserve the submission byte for byte when a clinical
   adapter is live;
2. does the assertion actually rescue the variant the bar would otherwise drop;
3. does turning the lookup back off lose it again, visibly and by id.

The assertion is fabricated for a fictional gene in a case whose own VCF header
says ``mva_synthetic=true``. It is a stand-in for ClinVar, is labelled as one, and
never leaves this file (GP-20, GP-40).
"""

from __future__ import annotations

import json
import shutil
from collections.abc import Mapping, Sequence
from pathlib import Path

import pytest

from mva.annotation.base import AdapterSet
from mva.annotation.triage import DROP_IMPACT_BELOW_FLOOR
from mva.config import CaseConfig, Workspace
from mva.evidence.discards import explain_variant
from mva.models.discard import DiscardStage
from mva.models.variant import ClinicalAssertion
from mva.orchestrator import PipelineResult, execute_pipeline

pytestmark = [pytest.mark.integration]

#: The synthetic case's only MODIFIER-on-every-transcript record: an intronic call
#: in the fictional SYNTHSOL5, SpliceAI 0.01, no coding or splice term. It is the
#: single variant in this fixture that the relevance bar drops, which makes it the
#: only one on which the clinical rescue is observable at all.
MODIFIER_ONLY_VARIANT = "GRCh38:chr3:10000000:G:A"

#: Artifacts that must not move by one byte when the bar is switched on — the same
#: two `test_staged_annotation_determinism.py` locks, so the two claims are
#: directly comparable and only the adapter set differs.
INVARIANT_ARTIFACTS: tuple[str, ...] = (
    "submission/track1_submission.csv",
    "candidates/ranked_pairs.json",
)


def _assertion() -> ClinicalAssertion:
    """A fabricated PATHOGENIC record, labelled as fabricated in every field."""
    return ClinicalAssertion(
        source="SYNTH_curation",
        version="v0.0-synthetic",
        accession="SYNTHVAR000042",
        significance="Pathogenic",
        review_status="criteria_provided,_single_submitter",
        star_rating=1,
        conditions=("Synthetic transport disorder (demo)",),
    )


class SpeakingClinicalAdapter:
    """A clinical adapter with exactly one thing to say, and a record of being asked.

    ``asked`` is the measurement that separates "the lookup ran in front of the
    bar" from "the lookup ran behind it", which is the whole of ADR 0028 decision
    6 and is invisible from the outputs alone.
    """

    def __init__(self, table: Mapping[str, tuple[ClinicalAssertion, ...]]) -> None:
        self._table = dict(table)
        self.asked: list[str] = []

    @property
    def name(self) -> str:
        return "synthetic-clinical-standin"

    @property
    def version(self) -> str:
        return "v0.0-synthetic"

    @property
    def synthetic(self) -> bool:
        """True. GP-20: a mock that does not declare itself is the failure."""
        return True

    def assertions(self, variant_ids: Sequence[str]) -> Mapping[str, tuple[ClinicalAssertion, ...]]:
        self.asked.extend(variant_ids)
        return {vid: self._table[vid] for vid in variant_ids if vid in self._table}


@pytest.fixture
def fixtures(repo_root: Path) -> Path:
    return repo_root / "tests" / "fixtures" / "synthetic"


def _workspace(root: Path, fixtures: Path) -> Workspace:
    inputs = root / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    shutil.copy(fixtures / "synthetic_case.vcf", inputs / "synthetic_case.vcf")
    shutil.copy(fixtures / "synthetic_phenotype.tsv", inputs / "synthetic_phenotype.tsv")
    return Workspace(root=root.resolve(), repo_root=fixtures.parents[2])


@pytest.fixture
def speaking_clinical_adapter(monkeypatch: pytest.MonkeyPatch) -> SpeakingClinicalAdapter:
    """Swap the run's clinical slot for one that answers, leaving the rest alone.

    Patched at the composition root rather than by editing a knowledge table: the
    tables are hash-pinned in `knowledge/manifests/knowledge.yaml` and locked by
    `tests/golden/test_locked_files.py`, and inventing clinical significance for
    fictional variants inside the shipped knowledge base is exactly what
    `NullClinicalAdapter`'s docstring refuses to do. It is legitimate in one test's
    process; it would not be legitimate on disk.
    """
    from mva import orchestrator

    adapter = SpeakingClinicalAdapter({MODIFIER_ONLY_VARIANT: (_assertion(),)})
    original = orchestrator.load_default_adapters

    def _load_with_clinical(knowledge_root: Path, manifest_path: Path) -> AdapterSet:
        shipped = original(knowledge_root, manifest_path)
        return AdapterSet(
            consequence=shipped.consequence,
            frequency=shipped.frequency,
            clinical=adapter,
        )

    monkeypatch.setattr(orchestrator, "load_default_adapters", _load_with_clinical)
    return adapter


def _with(config: CaseConfig, **triage: object) -> CaseConfig:
    return config.model_copy(update={"triage": config.triage.model_copy(update=triage)})


def _selection_on(config: CaseConfig) -> CaseConfig:
    enabled = config.selection.model_copy(update={"enabled": True})
    return config.model_copy(update={"selection": enabled})


def _run(config: CaseConfig, root: Path, fixtures: Path) -> tuple[Path, PipelineResult]:
    workspace = _workspace(root, fixtures)
    result = execute_pipeline(config, workspace)
    return workspace.runs_dir / result.run_id, result


def _annotated_ids(run_dir: Path) -> set[str]:
    rows = json.loads((run_dir / "variants" / "annotated.json").read_text(encoding="utf-8"))
    return {
        f"{r['coordinate']['build']}:{r['coordinate']['contig']}:"
        f"{r['coordinate']['position']}:{r['coordinate']['ref']}:{r['coordinate']['alt']}"
        for r in rows
    }


# ---------------------------------------------------------------------------


def test_the_partition_is_byte_identical_with_a_live_clinical_adapter(
    synthetic_config: CaseConfig,
    tmp_path: Path,
    fixtures: Path,
    speaking_clinical_adapter: SpeakingClinicalAdapter,
) -> None:
    """The determinism claim, re-made in the configuration that can break it.

    The existing byte-identity tests cannot exercise the bar's clinical rule or
    selection's clinical rescue, because their clinical adapter returns nothing.
    This runs the same comparison with one live assertion in play, under
    ``selection.enabled: true`` where the bar is at its configured width.
    """
    config = _selection_on(synthetic_config)
    off_dir, _ = _run(_with(config, enabled=False), tmp_path / "off", fixtures)
    on_dir, _ = _run(_with(config, enabled=True), tmp_path / "on", fixtures)

    for relative in INVARIANT_ARTIFACTS:
        assert (on_dir / relative).read_bytes() == (off_dir / relative).read_bytes(), (
            f"{relative} differs with the relevance bar on, while a clinical adapter "
            "is live. The partition changed the answer (ADR 0028, GP-30)."
        )


def test_a_curated_assertion_rescues_the_one_variant_the_bar_would_drop(
    synthetic_config: CaseConfig,
    tmp_path: Path,
    fixtures: Path,
    speaking_clinical_adapter: SpeakingClinicalAdapter,
) -> None:
    """Decision 6, on the shipped case: the rescue fires, and it fires on evidence.

    Without the assertion this variant is the case's one triage casualty — the
    fixture that makes `test_the_bar_is_not_vacuous_with_selection_enabled` pass.
    With it, it must survive, and the adapter must have been asked about it BEFORE
    the verdict, or the rescue was luck rather than evidence.
    """
    config = _selection_on(synthetic_config)
    run_dir, _ = _run(_with(config, enabled=True), tmp_path / "on", fixtures)

    assert MODIFIER_ONLY_VARIANT in speaking_clinical_adapter.asked, (
        "the clinical adapter was never asked about the variant it holds an assertion "
        "for, so triage.clinical_before_bar is not doing what ADR 0028 decision 6 says"
    )
    assert MODIFIER_ONLY_VARIANT in _annotated_ids(run_dir), (
        "a MODIFIER-only variant with a curated PATHOGENIC assertion was dropped by "
        "the relevance bar. This is the 3,616 -> 3,600 regression, on the shipped case."
    )
    ledger = run_dir / "discards" / "discard_ledger.tsv"
    assert not explain_variant(ledger, MODIFIER_ONLY_VARIANT), (
        "the rescued variant is in the discard ledger and in the annotated artifact"
    )


def test_deferring_the_clinical_lookup_loses_it_again_visibly(
    synthetic_config: CaseConfig,
    tmp_path: Path,
    fixtures: Path,
    speaking_clinical_adapter: SpeakingClinicalAdapter,
) -> None:
    """`clinical_before_bar: false` reproduces the bug on the shipped case.

    The setting exists deliberately and `triage_warnings` states its cost, but a
    cost nothing demonstrates is a claim. This is the demonstration, and it also
    checks the mitigation ADR 0028 offers in its place: the lost variant is in the
    ledger by id, with the value the rule read, so a targeted ClinVar scan over the
    dropped ids closes the hole without a re-run.
    """
    config = _selection_on(synthetic_config)
    run_dir, _ = _run(
        _with(config, enabled=True, clinical_before_bar=False), tmp_path / "deferred", fixtures
    )

    assert MODIFIER_ONLY_VARIANT not in _annotated_ids(run_dir), (
        "triage.clinical_before_bar=false no longer loses the clinically rescuable "
        "variant. Either the hole was closed — in which case the warning in "
        "triage_warnings is now false — or this fixture stopped depending on the rule."
    )
    (row,) = explain_variant(run_dir / "discards" / "discard_ledger.tsv", MODIFIER_ONLY_VARIANT)
    assert row.stage is DiscardStage.TRIAGE
    assert row.rule == DROP_IMPACT_BELOW_FLOOR
    assert row.value == "modifier", (
        "the mitigation for this hole is a targeted scan over the dropped ids, which "
        "needs the row to name the value the bar read (ADR 0028)"
    )


def test_the_run_warns_when_the_rescue_is_unreachable(
    synthetic_config: CaseConfig,
    tmp_path: Path,
    fixtures: Path,
    speaking_clinical_adapter: SpeakingClinicalAdapter,
) -> None:
    """A hole the run does not mention is a hole nobody will look for.

    `narrow_to_selection` silently widens the bar on every other axis where it
    could be narrower than selection. This is the one axis it does NOT close, so
    the warning is the entire mitigation and it has to be present in the artifact
    a reviewer reads.
    """
    config = _selection_on(synthetic_config)
    _, result = _run(
        _with(config, enabled=True, clinical_before_bar=False), tmp_path / "deferred", fixtures
    )
    warnings = " ".join(result.warnings)
    assert "clinical_before_bar" in warnings, (
        "the run does not warn that the clinical rescue is unreachable. The setting is "
        "allowed; leaving it unannounced is not (ADR 0028)."
    )
