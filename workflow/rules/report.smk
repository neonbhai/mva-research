# ---------------------------------------------------------------------------
# Reporting, provenance and the privacy gate
# ---------------------------------------------------------------------------
# The two report rules are the DAG's two exits (docs/architecture.md). Both
# invoke the SAME `mva run report` command — the CLI surface has one report
# subcommand, which renders the full report set for the run — and they differ
# only in which artifacts each one declares, and therefore in where each sits in
# the DAG. That is safe precisely because repeat runs are byte-identical (GP-30):
# the second invocation rewrites the Track 1 files with the same bytes.
# `report_track2` takes the Track 1 artifacts as an input so the two are ordered
# rather than racing.
#
# Rendering is gated twice (GP-43): an allowlist decides what may be emitted, and
# the rendered bytes are re-scanned afterwards. Classification is a claim; the
# scan is the verification.


rule report_track1:
    """Ranked pairs -> the Track 1 submission CSV and its candidate dossier.

    The submission format is a verified contract, not a guess: exact column
    order, at most 10 rows, one row per compound-het pair rather than two, and
    chromosome names compared raw — so `chr15` and `15` are different variants to
    the scorer. See docs/references/track1-submission-contract.md.
    """
    input:
        ranked=RANKED_PAIRS,
    output:
        submission=TRACK1_SUBMISSION,
        dossier=TRACK1_DOSSIER,
    log:
        f"{LOG_DIR}/report-track1.log",
    shell:
        mva_run("report")


rule report_track2:
    """Drug hypotheses -> the Track 2 mechanism-and-intervention report."""
    input:
        accepted=DRUG_HYPOTHESES,
        rejected=REJECTION_RECORD,
        track1=TRACK1_ARTIFACTS,
    output:
        report=TRACK2_REPORT,
    log:
        f"{LOG_DIR}/report-track2.log",
    shell:
        mva_run("report")


rule provenance:
    """The join point: the run manifest exists and the artifact set is complete.

    There is deliberately no `mva run provenance` subcommand to call here. Each
    stage registers its own ArtifactProvenance as it writes (GP-31), and the
    manifest is rewritten by whichever stage finishes the run; a manifest
    assembled after the fact by a separate pass could only describe what it
    found, not witness what happened. So this rule asserts rather than produces:
    it is the DAG node that says the artifact set is complete, so the privacy
    audit runs against a finished run instead of a partial one.
    """
    input:
        VALIDATED_OK,
        NORMALISED_VARIANTS,
        QC_REPORT,
        ANNOTATED_VARIANTS,
        RANKED_PAIRS,
        TRACK1_ARTIFACTS,
        MECHANISM_REPORT,
        DRUG_HYPOTHESES,
        REJECTION_RECORD,
        TRACK2_ARTIFACTS,
    output:
        touch(PROVENANCE_OK),
    params:
        manifest=PROVENANCE_MANIFEST,
        evidence_db=EVIDENCE_DB,
    shell:
        "test -s {params.manifest} && test -s {params.evidence_db}"


rule privacy_audit:
    """Blocking gate: the repository holds no patient data (PRIV-01, ADR 0009).

    Last in the DAG on purpose. It scans the REPOSITORY, not the workspace, and
    reports paths, line numbers, span lengths and rule IDs — never the matched
    bytes, because a scanner that echoes what it finds is itself the leak
    (GP-41). Nothing downstream consumes its output; it exists to fail the run.
    """
    input:
        PROVENANCE_OK,
    output:
        touch(PRIVACY_AUDIT_OK),
    log:
        f"{LOG_DIR}/privacy-audit.log",
    params:
        repo=str(REPO_ROOT),
        workspace=WORKSPACE,
    shell:
        "{MVA} privacy audit --repo {params.repo:q} --workspace {params.workspace:q}"
        " > {log} 2>&1"
