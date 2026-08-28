"""End-to-end stage wiring.

The stage graph lives here and nowhere else (GP-03). Stages are pure functions
over typed models; this module sequences them, writes their artifacts with
provenance, and accumulates the evidence ledger.

Ordering:

    validate -> ingest -> annotate -> select -> prioritise -> mechanism -> drugs
             -> report -> evidence persistence -> provenance -> privacy audit
"""

from __future__ import annotations

from collections.abc import Generator, Iterator
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass, field
from pathlib import Path

from mva.alleles import ReferenceLookup
from mva.annotation.binding import (
    BoundAdapters,
    ResolvedResources,
    build_real_adapter_set,
    resolve_real_resources,
)
from mva.annotation.local_tables import load_default_adapters
from mva.annotation.service import iter_annotated
from mva.config import CaseConfig, NetworkProfile, Workspace
from mva.determinism import canonical_json
from mva.errors import ConfigError, ExportBlockedError
from mva.evidence.ledger import AssertionResolver, EvidenceLedger
from mva.evidence.store import EvidenceStore, GraphEdge
from mva.ingestion.normalise import normalise_variants, open_reference_fasta
from mva.ingestion.qc import assess_quality
from mva.ingestion.reader import read_vcf
from mva.interventions.catalog import DrugCatalog
from mva.interventions.generate import generate_drug_hypotheses
from mva.mechanisms.builder import build_mechanism, mechanism_relevance_score
from mva.mechanisms.library import MechanismLibrary
from mva.models.base import Sensitivity
from mva.models.drug import DrugHypothesis
from mva.models.evidence import EvidenceItem
from mva.models.mechanism import MechanismHypothesis
from mva.models.pair import CandidatePair
from mva.models.provenance import ArtifactKind, ArtifactProvenance, RunManifest
from mva.models.variant import VariantRecord
from mva.phenotype.hpo import GenePhenotypeIndex
from mva.phenotype.loader import load_phenotype_profile
from mva.phenotype.scoring import score_all_genes
from mva.pipeline import (
    RunContext,
    artifact_digest,
    build_run_manifest,
    reference_versions_from_manifest,
    validate_case,
    write_provenance_manifest,
)
from mva.prioritization.filters import iter_hard_filtered, iter_soft_flagged
from mva.prioritization.pairing import generate_pair_candidates
from mva.prioritization.ranking import assign_discriminating_experiments, rank_pairs
from mva.prioritization.scoring import score_pair
from mva.prioritization.selection import iter_selected
from mva.privacy.audit import run_audit
from mva.privacy.export import PUBLIC_EXPORT_ALLOWLIST as _PUBLIC_EXPORT_ALLOWLIST
from mva.privacy.export import gate_public_export
from mva.privacy.netguard import OfflineProfile
from mva.reporting.dossier import build_candidate_dossier
from mva.reporting.track1 import (
    build_submission_rows,
    render_submission_csv,
    validate_submission,
)
from mva.reporting.track2 import (
    build_drug_report,
    build_mechanism_report,
    build_rejection_record,
    build_track2_report,
)
from mva.resources import (
    ResourceRoot,
    ResourceRootError,
    load_resource_manifest,
    reference_fasta_path,
    resolve_resource_root,
)

#: Stage names in execution order, for `--stop-after`.
STAGES: tuple[str, ...] = (
    "validate",
    "ingest",
    "annotate",
    "select",
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
    #: Every drug hypothesis considered, accepted and rejected alike. Exposed so
    #: acceptance tests can assert against what the pipeline DID, rather than
    #: against the expectation file they were given (GP-19).
    accepted_drugs: tuple[DrugHypothesis, ...] = field(default_factory=tuple)
    rejected_drugs: tuple[DrugHypothesis, ...] = field(default_factory=tuple)

    def drug(self, drug_id: str) -> DrugHypothesis | None:
        """Look up a drug hypothesis by ID across both accepted and rejected."""
        for hypothesis in (*self.accepted_drugs, *self.rejected_drugs):
            if hypothesis.drug_id == drug_id:
                return hypothesis
        return None


#: Re-exported, not redeclared. This module previously carried its own byte-equal
#: copy of the list in `mva.privacy.export`, which is the arrangement that lets the
#: two drift silently — and the policy is supposed to have exactly one home (see
#: that module's docstring). The name stays importable from here because the
#: composition root is where a reader looks for what this run may publish.
PUBLIC_EXPORT_ALLOWLIST: tuple[str, ...] = _PUBLIC_EXPORT_ALLOWLIST


def _gate_public_artifact(context: RunContext, artifact: ArtifactProvenance) -> None:
    """Run the public-export gate over an artifact classified PUBLIC.

    Wired to `RunContext.on_register`, so it runs on EVERY artifact the run
    produces, at the moment it is produced. It previously ran on the submission
    and on nothing else, while `MECHANISM_REPORT`, `DRUG_HYPOTHESES`,
    `REJECTION_RECORD` and `TRACK2_REPORT` were all classified PUBLIC and all on
    the export allowlist — so four of the five artifacts a reader would actually
    publish were gated zero times, and `docs/architecture.md`'s "gated twice"
    described the submission alone.

    Non-PUBLIC artifacts no-op here, which is what makes gating everything cheap
    enough to do unconditionally: the scan cost is paid only by the files the
    classification already claims are publishable.

    Classification is a claim; the allowlist match and the content re-scan are the
    verification (GP-43). A refusal raises: an invalid public artifact sitting on
    disk is one someone will find and upload.
    """
    if artifact.sensitivity is not Sensitivity.PUBLIC:
        return
    path = context.workspace.root / artifact.relative_path
    decision = gate_public_export(
        path, declared=artifact.sensitivity, allowlist=PUBLIC_EXPORT_ALLOWLIST
    )
    if not decision.allowed:
        msg = (
            f"Artifact {artifact.relative_path!r} is classified PUBLIC but failed the "
            f"export gate: {'; '.join(decision.reasons)}"
        )
        raise ExportBlockedError(msg)


def _resolve_case_resources(config: CaseConfig, workspace: Workspace) -> ResolvedResources | None:
    """Which reference releases this case runs against, or ``None`` for synthetic.

    **This is the policy the whole real-data path turns on, so it is stated here,
    at the composition root, rather than buried in an adapter factory (GP-03).**

    * ``synthetic: true`` -> ``None``. The case runs on the hash-pinned demo tables
      under ``knowledge/public/``. That is what `just demo`, `just
      demo-determinism` and every test do, and it must stay independent of whether
      the machine happens to have 200 GB of reference data: a suite whose result
      depends on ``$MVA_RESOURCES`` is a suite that passes in one place and fails
      in another. It also keeps the fictional ``SYNTH*`` genes away from a real
      gene model that has never heard of them.
    * ``synthetic: false`` -> the real releases, **required**. No resource root, an
      unregistered file, an unfetched entry or a failed integrity pin all raise.

    There is deliberately no third branch, and in particular no silent fallback for
    a real case with missing resources. The synthetic tables hold fictional genes
    and invented allele frequencies; a real proband ranked against them would
    produce a submission, a dossier and a provenance manifest that all looked
    entirely healthy and were entirely fabricated (GP-20). ADR 0027 records the
    trade — a loud failure that stops a run is recoverable in minutes; a quiet one
    that ships is not recoverable at all.
    """
    if config.synthetic:
        return None

    root = _require_resource_root(config, workspace)
    manifest_path = workspace.repo_root / config.resources.manifest
    manifest = load_resource_manifest(manifest_path)
    return resolve_real_resources(config, resource_root=root, manifest=manifest)


def _require_resource_root(config: CaseConfig, workspace: Workspace) -> ResourceRoot:
    """The external reference-data root, or a refusal that says why it matters."""
    try:
        return resolve_resource_root(repo_root=workspace.repo_root)
    except ResourceRootError as exc:
        msg = (
            f"Case {config.case_id!r} declares synthetic=false, so it must be annotated "
            f"against the real public reference releases — but: {exc}\n"
            "Refusing to fall back to the synthetic tables under knowledge/public/. They "
            "hold fictional genes and invented allele frequencies, so this run would rank a "
            "real proband against fabricated evidence and every artifact it produced would "
            "look healthy (GP-20, ADR 0027). Set MVA_RESOURCES, or set synthetic=true if "
            "this case really is a demo."
        )
        raise ConfigError(msg) from exc


def _open_reference(
    config: CaseConfig,
    workspace: Workspace,
    *,
    resolved: ResolvedResources | None,
    stack: ExitStack,
) -> ReferenceLookup | None:
    """Open the GRCh38 FASTA that left-alignment and both joins need.

    For a real case this is the verified release :func:`_resolve_case_resources`
    already required, so it is never ``None``. For a synthetic case it is the
    per-case ``inputs.reference_fasta`` override if one is set, and otherwise
    ``None`` — which is the state every fixture and the demo run in, and which
    ``normalise_variants`` reports as
    :attr:`~mva.alleles.LeftAlignmentStatus.UNAVAILABLE_NO_REFERENCE` rather than
    passing over in silence.
    """
    path = (
        resolved.reference_fasta
        if resolved is not None
        else reference_fasta_path(config, workspace=workspace, resource_root=None)
    )
    if path is None:
        return None
    return stack.enter_context(_closing_reference(path))


@contextmanager
def _closing_reference(path: Path) -> Generator[ReferenceLookup]:
    """Open an indexed FASTA and guarantee the htslib handle is released."""
    handle = open_reference_fasta(path)
    try:
        yield handle
    finally:
        handle.close()


def _bind_adapters(
    *,
    knowledge_root: Path,
    manifest_path: Path,
    resolved: ResolvedResources | None,
    reference: ReferenceLookup | None,
    stack: ExitStack,
) -> BoundAdapters:
    """Bind the annotation adapter set this case is entitled to.

    The synthetic branch is byte-for-byte what the pipeline did before the real
    adapters existed. The real branch passes ``reference=`` to **both** joining
    adapters, which is the whole point: constructed without it they still work,
    still pass their tests, and silently lose every repeat-tract indel join
    (ADR 0018). ``build_real_adapter_set`` takes the keyword without a default so
    it cannot be forgotten here, and reports the adapters' own
    ``representation_limitation`` when it is ``None``.
    """
    bound = (
        BoundAdapters(adapters=load_default_adapters(knowledge_root, manifest_path))
        if resolved is None
        else build_real_adapter_set(resolved, reference=reference)
    )
    stack.callback(bound.close)
    return bound


def _should_run(stage: str, stop_after: str | None) -> bool:
    if stop_after is None:
        return True
    if stop_after not in STAGES:
        msg = f"Unknown stage {stop_after!r}. Valid stages: {', '.join(STAGES)}"
        raise ConfigError(msg)
    return STAGES.index(stage) <= STAGES.index(stop_after)


def execute_pipeline(
    config: CaseConfig,
    workspace: Workspace,
    *,
    allow_workspace_in_repo: bool = False,
    stop_after: str | None = None,
) -> PipelineResult:
    """Run the pipeline and write every artifact, with the offline guard ARMED.

    The arming lives here and not only in `mva.cli` because `execute_pipeline` is
    what actually reads the VCF. Every other entry point — Snakemake, a test, a
    notebook, a future service — reaches patient data through this function, and a
    guard that only the CLI installs is a guard the next entry point does not have.
    `mva.cli._install_privacy_guards` still arms for the whole process; the two
    nest, which `OfflineProfile` is explicitly built for.

    `strict=False` is deliberate and not a weakening: the run shells out to `git`
    for the provenance manifest and the privacy audit, and `strict` blocks
    `subprocess.Popen`. Strict mode would make every run fail at the manifest step,
    which is how a control gets switched off. What is denied here is Python-level
    outbound network — the realistic accident of a stage calling `requests.get`.
    The honest limits (C extensions, spawned subprocesses, `ctypes`) are in
    `mva.privacy.netguard`, and closing them needs an OS control (TD-06).
    """
    if config.network_profile is NetworkProfile.ONLINE:
        return _execute_stages(
            config,
            workspace,
            allow_workspace_in_repo=allow_workspace_in_repo,
            stop_after=stop_after,
        )
    with OfflineProfile(workspace_root=workspace.root):
        return _execute_stages(
            config,
            workspace,
            allow_workspace_in_repo=allow_workspace_in_repo,
            stop_after=stop_after,
        )


def _execute_stages(  # noqa: PLR0915 - the composition root is legitimately long
    config: CaseConfig,
    workspace: Workspace,
    *,
    allow_workspace_in_repo: bool = False,
    stop_after: str | None = None,
) -> PipelineResult:
    """The stage graph itself. Call `execute_pipeline`, which arms the guard first."""
    repo_root = workspace.repo_root
    knowledge_root = repo_root / KNOWLEDGE_ROOT_NAME

    # ---------------------------------------------------------------- validate
    # Reference hashes are resolved BEFORE the run id is derived: two runs over
    # different knowledge tables must not collide in the same run directory
    # (they previously did, silently overwriting each other's submission).
    manifest_path = knowledge_root / "manifests" / "knowledge.yaml"
    reference_versions = reference_versions_from_manifest(manifest_path)

    context = validate_case(
        config,
        workspace,
        allow_workspace_in_repo=allow_workspace_in_repo,
        reference_versions=reference_versions,
    )
    # Wired BEFORE the first artifact is written, so there is no window in which an
    # artifact is registered ungated. See `_gate_public_artifact`.
    context.on_register = _gate_public_artifact
    # The spill directory is inside the workspace: evidence subjects are variant
    # IDs, so the file carries proband coordinates and is SENSITIVE (GP-40).
    # `close()` in `_finish` unlinks it. Without a spill_dir the ledger keeps its
    # pre-existing all-in-memory behaviour, which is what every test gets.
    ledger = EvidenceLedger(run_id=context.run_id, spill_dir=workspace.tmp_dir)

    # Resolved BEFORE the VCF is opened. A real case whose reference releases are
    # missing must fail before it has read one patient record, not after (ADR 0027).
    resolved = _resolve_case_resources(config, workspace)

    # `ExitStack` spans ingest and annotate because the SAME reference FASTA handle
    # serves both: normalisation left-aligns against it, and the ClinVar and gnomAD
    # adapters reconcile a shifted indel against it. Two handles would be two
    # caches over one immutable file. Everything is released as the block exits,
    # including on an early `return _finish(...)` inside it.
    with ExitStack() as stack:
        reference = _open_reference(config, workspace, resolved=resolved, stack=stack)

        # ---------------------------------------------------------------- ingest
        ingestion = read_vcf(
            workspace.path(config.inputs.vcf),
            expected_build=config.genome_build,
            source_artifact="input_vcf",
        )
        normalised = normalise_variants(ingestion.variants, reference=reference)
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
                # The typed degraded state, not only the prose in `warnings`. A
                # reader deciding whether an absent ClinVar assertion means
                # anything has to be able to branch on this (GP-14, ADR 0018).
                "left_alignment": normalised.left_alignment.as_dict(),
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
        bound = _bind_adapters(
            knowledge_root=knowledge_root,
            manifest_path=manifest_path,
            resolved=resolved,
            reference=reference,
            stack=stack,
        )
        annotated = iter_annotated(qc.variants, adapters=bound.adapters, clock=context.clock)

        # ONE pass drives everything: annotation feeds the ledger and the artifact,
        # then the hard filter, the soft flags and selection. Nothing but the
        # selected variants accumulates, and that is a few hundred records rather
        # than 4.5 M (docs/scale-report.md, docs/handoff-scale.md §1.3).
        #
        # The artifact is written from the point BEFORE hard filtering, so
        # `variants/annotated.json` still holds every annotated record. That is why
        # this uses the PUSH sink: a puller cannot both write every record and hand
        # a filtered subset onward in the same pass.
        with context.open_json_rows_artifact(
            "variants/annotated.json",
            kind=ArtifactKind.ANNOTATED_VARIANTS,
            stage="annotate",
            upstream=[normalised_art.artifact_id],
        ) as sink:

            def _recorded() -> Iterator[VariantRecord]:
                for item in annotated:
                    ledger.extend(item.evidence)
                    sink.write(item.variant.model_dump(mode="json"))
                    yield item.variant

            filtered = iter_hard_filtered(_recorded(), expected_build=config.genome_build)
            flagged = iter_soft_flagged(
                filtered, frequency=config.frequency, quality=config.quality
            )
            selection = iter_selected(
                flagged,
                frequency=config.frequency,
                thresholds=config.selection,
                clock=context.clock,
            )
            selected = list(selection)

        annotated_art = sink.provenance
        if annotated_art is None:  # pragma: no cover - set by the context manager
            msg = "The annotated-variants artifact was not registered on exit."
            raise ConfigError(msg)
        context.warnings.extend(annotated.warnings())
        # Coverage holes are known only once the adapters have been asked, so this
        # is read after the pass, not at bind time.
        context.warnings.extend(bound.run_warnings())

        # ---------------------------------------------------------------- select
        selection_report = selection.report()
        ledger.extend(selection.evidence())
        context.warnings.extend(selection_report.warnings)
        context.write_json_artifact(
            "selection/selection_report.json",
            {
                **selection_report.as_payload(),
                "hard_filter": filtered.counts(),
            },
            kind=ArtifactKind.SELECTION_REPORT,
            stage="select",
            upstream=[annotated_art.artifact_id],
            row_count=selection_report.input_count,
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
    # Hard filtering and soft flagging happened inside the annotate pass above; the
    # records in `selected` have been through both and carry identical flags.
    #
    # ADR 0013: bound the hypothesis space by plausibility, never by coordinate, and
    # surface it when a cap fires. A cap that silently deletes the correct pair is the
    # worst failure this pipeline can have, so the warning is part of the contract.
    pairing = generate_pair_candidates(
        selected,
        max_pairs_per_gene=config.max_pairs_per_gene,
        frequency=config.frequency,
    )
    candidates = pairing.candidates
    context.warnings.extend(pairing.warnings)

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
    submission_text = render_submission_csv(rows)

    # Self-check BEFORE the bytes reach disk: a malformed submission that exists is
    # worse than one that does not, because it will be uploaded.
    valid, submission_errors = validate_submission(submission_text)
    if not valid:
        msg = "Rendered Track 1 submission failed its own contract check: " + "; ".join(
            submission_errors
        )
        raise ExportBlockedError(msg)

    context.write_text_artifact(
        "submission/track1_submission.csv",
        submission_text,
        kind=ArtifactKind.SUBMISSION,
        stage="report",
        upstream=[pairs_art.artifact_id],
        row_count=len(rows),
    )
    # No explicit gate call here any more: `context.on_register` ran it as
    # `submission_art` was registered, along with every other artifact.
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

    _write_privacy_audit(context, repo_root=repo_root)
    # The spill file carries proband coordinates and the run is over. `close()` is a
    # no-op on a ledger that never spilled, and `len(ledger)` keeps working after it;
    # reading the ledger's CONTENTS afterwards raises rather than returning an empty
    # ledger, which would read as a run that produced no evidence.
    ledger.close()

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
        accepted_drugs=tuple(d for d in accepted if isinstance(d, DrugHypothesis)),
        rejected_drugs=tuple(d for d in rejected if isinstance(d, DrugHypothesis)),
    )


def _write_privacy_audit(context: RunContext, *, repo_root: Path) -> None:
    """Run the privacy audit over the repo and the workspace, and record it.

    Written as part of every run rather than as a separate manual step: an audit
    you have to remember to run is one you will forget on the run that mattered.
    The report contains paths, line numbers and counts only — never matched
    content (GP-41).
    """
    report = run_audit(repo_root, workspace=context.workspace.root)
    context.write_text_artifact(
        "privacy/privacy_audit.md",
        report.to_markdown(),
        kind=ArtifactKind.PRIVACY_AUDIT,
        stage="privacy",
        row_count=len(report.results),
    )
    if not report.passed:
        context.warnings.append(
            "Privacy audit FAILED for this run: "
            f"{', '.join(report.failed_checks)}. See privacy/privacy_audit.md. "
            "Do not export or share any artifact from this run until it passes."
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
        # `items()` materialises; a spilled ledger holds more items than the process
        # has memory for, which is why it spilled.
        store.write_evidence_stream(ledger.iter_items())
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
        # The Parquet export is the artifact determinism is asserted on; the
        # DuckDB file is a storage-engine container. See verify_determinism.
        exported = store.export_parquet(db_path.parent / "parquet")

    context.register_artifact(
        kind=ArtifactKind.EVIDENCE_DB,
        path=db_path,
        stage="evidence",
        row_count=sum(counts.values()),
        notes=canonical_json(counts),
    )
    for table in sorted(exported):
        context.register_artifact(
            kind=ArtifactKind.EVIDENCE_DB,
            path=exported[table],
            stage="evidence",
            row_count=counts.get(table),
            notes=f"Parquet export of table {table!r}.",
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
