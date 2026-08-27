"""End-to-end stage wiring.

The stage graph lives here and nowhere else (GP-03). Stages are pure functions
over typed models; this module sequences them, writes their artifacts with
provenance, and accumulates the evidence ledger.

Ordering:

    validate -> ingest -> annotate -> prioritise -> mechanism -> drugs -> report
             -> evidence persistence -> provenance -> privacy audit
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from mva.annotation.local_tables import load_default_adapters
from mva.annotation.service import annotate_variants
from mva.config import CaseConfig, Workspace
from mva.determinism import canonical_json
from mva.errors import ConfigError
from mva.evidence.ledger import AssertionResolver, EvidenceLedger
from mva.evidence.store import EvidenceStore, GraphEdge
from mva.ingestion.normalise import normalise_variants
from mva.ingestion.qc import assess_quality
from mva.ingestion.reader import read_vcf
from mva.interventions.catalog import DrugCatalog
from mva.interventions.generate import generate_drug_hypotheses
from mva.mechanisms.builder import build_mechanism, mechanism_relevance_score
from mva.mechanisms.library import MechanismLibrary
from mva.models.evidence import EvidenceItem
from mva.models.mechanism import MechanismHypothesis
from mva.models.pair import CandidatePair
from mva.models.provenance import ArtifactKind, RunManifest
from mva.phenotype.hpo import GenePhenotypeIndex
from mva.phenotype.loader import load_phenotype_profile
from mva.phenotype.scoring import score_all_genes
from mva.pipeline import (
    RunContext,
    artifact_digest,
    build_run_manifest,
    validate_case,
    write_provenance_manifest,
)
from mva.prioritization.filters import apply_hard_filters, apply_soft_flags
from mva.prioritization.pairing import generate_pairs
from mva.prioritization.ranking import assign_discriminating_experiments, rank_pairs
from mva.prioritization.scoring import score_pair
from mva.reporting.dossier import build_candidate_dossier
from mva.reporting.track1 import build_submission_rows, render_submission_csv
from mva.reporting.track2 import (
    build_drug_report,
    build_mechanism_report,
    build_rejection_record,
    build_track2_report,
)

#: Stage names in execution order, for `--stop-after`.
STAGES: tuple[str, ...] = (
    "validate",
    "ingest",
    "annotate",
    "prioritise",
    "mechanism",
    "drugs",
    "report",
)

#: The knowledge root is repo-relative and PUBLIC. It never contains patient data,
#: which is why it may live inside the repository while the workspace may not.
KNOWLEDGE_ROOT_NAME = "knowledge"


@dataclass
class PipelineResult:
    """Summary of a completed run. Counts and identifiers only — never records."""

    run_id: str
    manifest: RunManifest
    digest: dict[str, str]
    pair_count: int
    evidence_count: int
    drugs_accepted: int
    drugs_rejected: int
    top_gene: str | None
    warnings: tuple[str, ...] = ()
    ranked_pairs: tuple[CandidatePair, ...] = field(default_factory=tuple)
    mechanism: MechanismHypothesis | None = None


def _should_run(stage: str, stop_after: str | None) -> bool:
    if stop_after is None:
        return True
    if stop_after not in STAGES:
        msg = f"Unknown stage {stop_after!r}. Valid stages: {', '.join(STAGES)}"
        raise ConfigError(msg)
    return STAGES.index(stage) <= STAGES.index(stop_after)


def execute_pipeline(  # noqa: PLR0915 - the composition root is legitimately long
    config: CaseConfig,
    workspace: Workspace,
    *,
    allow_workspace_in_repo: bool = False,
    stop_after: str | None = None,
) -> PipelineResult:
    """Run the pipeline and write every artifact."""
    repo_root = workspace.repo_root
    knowledge_root = repo_root / KNOWLEDGE_ROOT_NAME

    # ---------------------------------------------------------------- validate
    context = validate_case(config, workspace, allow_workspace_in_repo=allow_workspace_in_repo)
    ledger = EvidenceLedger(run_id=context.run_id)

    # ---------------------------------------------------------------- ingest
    ingestion = read_vcf(
        workspace.path(config.inputs.vcf),
        expected_build=config.genome_build,
        source_artifact="input_vcf",
    )
    normalised = normalise_variants(ingestion.variants)
    qc = assess_quality(normalised.variants, thresholds=config.quality, clock=context.clock)
    ledger.extend(qc.evidence)
    context.warnings.extend(ingestion.warnings)
    context.warnings.extend(normalised.warnings)

    normalised_art = context.write_json_artifact(
        "variants/normalised.json",
        [v.model_dump(mode="json") for v in qc.variants],
        kind=ArtifactKind.NORMALISED_VARIANTS,
        stage="ingest",
        row_count=len(qc.variants),
    )
    context.write_json_artifact(
        "qc/qc_report.json",
        {
            "metrics": qc.metrics,
            "normalisation_operations": normalised.operations_applied,
            "skipped_count": ingestion.skipped_count,
            "skipped_reasons": list(ingestion.skipped_reasons),
        },
        kind=ArtifactKind.QC_REPORT,
        stage="ingest",
        upstream=[normalised_art.artifact_id],
    )
    if not _should_run("annotate", stop_after):
        return _finish(context, repo_root=repo_root, ledger=ledger, ranked=[])

    # ---------------------------------------------------------------- annotate
    adapters = load_default_adapters(
        knowledge_root, knowledge_root / "manifests" / "knowledge.yaml"
    )
    annotation = annotate_variants(qc.variants, adapters=adapters, clock=context.clock)
    ledger.extend(annotation.evidence)
    context.warnings.extend(annotation.warnings)

    annotated_art = context.write_json_artifact(
        "variants/annotated.json",
        [v.model_dump(mode="json") for v in annotation.variants],
        kind=ArtifactKind.ANNOTATED_VARIANTS,
        stage="annotate",
        upstream=[normalised_art.artifact_id],
        row_count=len(annotation.variants),
    )
    if not _should_run("prioritise", stop_after):
        return _finish(context, repo_root=repo_root, ledger=ledger, ranked=[])

    # ---------------------------------------------------------------- phenotype
    profile = load_phenotype_profile(
        workspace.path(config.inputs.phenotype),
        subject_id=config.proband_id,
        hpo_version="synthetic-v0.0",
        source_artifact="input_phenotype",
    )
    gene_index = GenePhenotypeIndex.from_tsv(
        knowledge_root / "public" / "gene_phenotype.tsv", version="synthetic-v0.0"
    )

    # ---------------------------------------------------------------- prioritise
    filtered = apply_hard_filters(annotation.variants, expected_build=config.genome_build)
    flagged = apply_soft_flags(
        filtered.retained, frequency=config.frequency, quality=config.quality
    )
    candidates = generate_pairs(flagged, max_pairs_per_gene=config.max_pairs_per_gene)

    gene_symbols = sorted({c.gene_symbol for c in candidates})
    phenotype_matches = score_all_genes(
        gene_symbols, profile=profile, index=gene_index, clock=context.clock
    )
    for match in phenotype_matches.values():
        ledger.extend(match.evidence)

    library = MechanismLibrary.from_tsv(
        knowledge_root / "public" / "mechanisms.tsv",
        knowledge_root / "public" / "mechanism_meta.tsv",
        version="synthetic-v0.0",
    )

    scored = []
    for candidate in candidates:
        phenotype_score = (
            phenotype_matches[candidate.gene_symbol].score
            if candidate.gene_symbol in phenotype_matches
            else 0.5
        )
        gene_mechanism = library.for_gene(candidate.gene_symbol)
        scored_pair = score_pair(
            candidate,
            phenotype_score=phenotype_score,
            mechanism_score=mechanism_relevance_score(gene_mechanism),
            weights=config.weights,
            phase_weights=config.phase_weights,
            frequency=config.frequency,
            quality=config.quality,
            clock=context.clock,
        )
        scored.append(scored_pair)
        ledger.extend(scored_pair.supporting_evidence)
        ledger.extend(scored_pair.contradicting_evidence)

    ranked = assign_discriminating_experiments(rank_pairs(scored, clock=context.clock))
    pairs_art = context.write_json_artifact(
        "candidates/ranked_pairs.json",
        [p.model_dump(mode="json") for p in ranked],
        kind=ArtifactKind.CANDIDATE_PAIRS,
        stage="prioritise",
        upstream=[annotated_art.artifact_id],
        row_count=len(ranked),
    )
    top_gene = ranked[0].gene_symbol if ranked else None

    # ---------------------------------------------------------------- track 1
    rows = build_submission_rows(
        ranked, proband_id=config.proband_id, max_rows=config.max_submission_rows
    )
    context.write_text_artifact(
        "submission/track1_submission.csv",
        render_submission_csv(rows),
        kind=ArtifactKind.SUBMISSION,
        stage="report",
        upstream=[pairs_art.artifact_id],
        row_count=len(rows),
    )
    resolver = AssertionResolver(ledger)
    context.write_text_artifact(
        "reports/candidate_dossier.md",
        build_candidate_dossier(ranked, resolver=resolver, clock=context.clock),
        kind=ArtifactKind.DOSSIER,
        stage="report",
        upstream=[pairs_art.artifact_id],
    )

    if not _should_run("mechanism", stop_after):
        return _finish(context, repo_root=repo_root, ledger=ledger, ranked=ranked)

    # ---------------------------------------------------------------- mechanism
    mechanism_result = build_mechanism(
        top_gene or "",
        pair_id=ranked[0].pair_id if ranked else None,
        library=library,
        clock=context.clock,
    )
    ledger.extend(mechanism_result.evidence)
    context.warnings.extend(mechanism_result.warnings)
    mechanism = mechanism_result.hypothesis

    if mechanism is not None:
        context.write_text_artifact(
            "reports/mechanism_report.md",
            build_mechanism_report(mechanism, resolver=resolver, clock=context.clock),
            kind=ArtifactKind.MECHANISM_REPORT,
            stage="mechanism",
            upstream=[pairs_art.artifact_id],
        )

    if mechanism is None or not _should_run("drugs", stop_after):
        return _finish(
            context, repo_root=repo_root, ledger=ledger, ranked=ranked, mechanism=mechanism
        )

    # ---------------------------------------------------------------- drugs
    catalog = DrugCatalog.from_tsv(
        knowledge_root / "public" / "drug_catalog.tsv", version="synthetic-v0.0"
    )
    triage = generate_drug_hypotheses(mechanism=mechanism, catalog=catalog, clock=context.clock)
    ledger.extend(triage.evidence)

    context.write_text_artifact(
        "reports/drug_hypotheses.md",
        build_drug_report(
            triage.accepted,
            triage.rejected,
            mechanism=mechanism,
            resolver=resolver,
            clock=context.clock,
        ),
        kind=ArtifactKind.DRUG_HYPOTHESES,
        stage="drugs",
        row_count=len(triage.accepted),
    )
    context.write_text_artifact(
        "reports/rejection_record.md",
        build_rejection_record(triage.rejected, clock=context.clock),
        kind=ArtifactKind.REJECTION_RECORD,
        stage="drugs",
        row_count=len(triage.rejected),
    )
    context.write_text_artifact(
        "reports/track2_report.md",
        build_track2_report(
            mechanism,
            triage.accepted,
            triage.rejected,
            pair=ranked[0] if ranked else None,
            resolver=resolver,
            clock=context.clock,
        ),
        kind=ArtifactKind.TRACK2_REPORT,
        stage="report",
    )

    return _finish(
        context,
        repo_root=repo_root,
        ledger=ledger,
        ranked=ranked,
        mechanism=mechanism,
        accepted=list(triage.accepted),
        rejected=list(triage.rejected),
        profile=profile,
    )


def _finish(
    context: RunContext,
    *,
    repo_root: Path,
    ledger: EvidenceLedger,
    ranked: list[CandidatePair] | tuple[CandidatePair, ...],
    mechanism: MechanismHypothesis | None = None,
    accepted: list[object] | None = None,
    rejected: list[object] | None = None,
    profile: object | None = None,
) -> PipelineResult:
    """Persist the evidence store, write provenance, and summarise."""
    accepted = accepted or []
    rejected = rejected or []
    _persist_evidence(
        context,
        ledger=ledger,
        ranked=ranked,
        mechanism=mechanism,
        accepted=accepted,
        rejected=rejected,
        profile=profile,
    )

    manifest = build_run_manifest(context, repo_root=repo_root)
    write_provenance_manifest(context, manifest)
    # Rebuild after the provenance artifact registers itself, so the manifest we
    # return matches what is on disk.
    manifest = build_run_manifest(context, repo_root=repo_root)

    return PipelineResult(
        run_id=context.run_id,
        manifest=manifest,
        digest=artifact_digest(context),
        pair_count=len(ranked),
        evidence_count=len(ledger),
        drugs_accepted=len(accepted),
        drugs_rejected=len(rejected),
        top_gene=ranked[0].gene_symbol if ranked else None,
        warnings=tuple(context.warnings),
        ranked_pairs=tuple(ranked),
        mechanism=mechanism,
    )


def _persist_evidence(
    context: RunContext,
    *,
    ledger: EvidenceLedger,
    ranked: list[CandidatePair] | tuple[CandidatePair, ...],
    mechanism: MechanismHypothesis | None,
    accepted: list[object],
    rejected: list[object],
    profile: object | None,
) -> None:
    """Write the evidence database and its Parquet export."""
    db_path = context.artifact_path("evidence/evidence.duckdb")
    with EvidenceStore(db_path) as store:
        store.initialise()
        store.write_evidence(ledger.items())
        if ranked:
            store.write_pairs(list(ranked))
            store.write_variants([v for pair in ranked for v in pair.variants])
            store.write_consequences([v for pair in ranked for v in pair.variants])
            store.write_frequencies([v for pair in ranked for v in pair.variants])
        if profile is not None:
            store.write_phenotypes(profile)  # type: ignore[arg-type]
        if mechanism is not None:
            store.write_mechanism(mechanism)
        drugs = [*accepted, *rejected]
        if drugs:
            store.write_drugs(drugs)  # type: ignore[arg-type]
        store.write_edges(
            _graph_edges(ranked=ranked, mechanism=mechanism, accepted=accepted, rejected=rejected)
        )
        counts = store.counts()

    context.register_artifact(
        kind=ArtifactKind.EVIDENCE_DB,
        path=db_path,
        stage="evidence",
        row_count=sum(counts.values()),
        notes=canonical_json(counts),
    )


def _graph_edges(
    *,
    ranked: list[CandidatePair] | tuple[CandidatePair, ...],
    mechanism: MechanismHypothesis | None,
    accepted: list[object],
    rejected: list[object],
) -> list[GraphEdge]:
    """Relational subject-predicate-object edges (ADR 0002).

    The domain is a graph; the storage is not. Keeping edges relational avoids a
    second system to secure and delete, and leaves the export path open.
    """
    edges: list[GraphEdge] = []
    for pair in ranked:
        for variant in pair.variants:
            edges.append(
                GraphEdge(
                    subject_id=variant.variant_id,
                    subject_kind="variant",
                    predicate="participates_in",
                    object_id=pair.pair_id,
                    object_kind="pair",
                )
            )
        edges.append(
            GraphEdge(
                subject_id=pair.pair_id,
                subject_kind="pair",
                predicate="implicates_gene",
                object_id=pair.gene_symbol,
                object_kind="gene",
                confidence=pair.composite_score,
            )
        )
    if mechanism is not None:
        edges.append(
            GraphEdge(
                subject_id=mechanism.gene_symbol,
                subject_kind="gene",
                predicate="has_mechanism",
                object_id=mechanism.mechanism_id,
                object_kind="mechanism",
            )
        )
        for drug in [*accepted, *rejected]:
            drug_id = getattr(drug, "drug_id", None)
            if not isinstance(drug_id, str):
                continue
            predicate = "rejected_for_mechanism" if drug in rejected else "targets_mechanism"
            edges.append(
                GraphEdge(
                    subject_id=drug_id,
                    subject_kind="drug",
                    predicate=predicate,
                    object_id=mechanism.mechanism_id,
                    object_kind="mechanism",
                )
            )
    return edges


def collect_evidence(items: list[EvidenceItem]) -> tuple[EvidenceItem, ...]:
    return tuple(items)
