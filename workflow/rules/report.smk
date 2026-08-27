# ---------------------------------------------------------------------------
# Reporting, provenance and the privacy gate
# ---------------------------------------------------------------------------
# The two report rules are the DAG's two exits. Both invoke the SAME
# `mva run report` command, which renders the full report set for the run; they
# differ only in which artifacts each declares and therefore in where they sit in
# the DAG. That is safe precisely because repeat runs are byte-identical
# (GP-30): the second invocation rewrites the Track 1 files with the same bytes.
# `report_track2` takes the Track 1 artifacts as input so the two are ordered
# rather than racing.
#
# Rendering is gated twice (GP-43): an allowlist decides what may be emitted, and
# the rendered bytes are re-scanned afterwards. Classification is a claim; the
# scan is the verification.


rule report_track1:
    """Ranked pairs -> the Track 1 submission CSV and its dossier.

    The submission format is a verified contract, not a guess: exact column
    order, at most 10 rows, one row per compound-het pair, and chromosome names
    compared raw so `chr15` and `15` are different variants to the scorer. See
    docs/references/track1-submission-contract.md.
    """
    input:
        pairs=CANDIDATE_PAIRS,
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
        rejected=REJECTED_DRUGS,
        track1=TRACK1_ARTIFACTS,
    output:
        report=TRACK2_REPORT,
    log:
        f"{LOG_DIR}/report-track2.log",
    shell:
        mva_run("report")


rule provenance:
    """The join point: every declared artifact of this run exists.

    There is deliberately no `mva run provenance` subcommand. Each stage
    registers its own ArtifactProvenance as it writes (GP-31); a manifest
    assembled after the fact could only describe what it found, not witness what
    happened. This rule is the DAG node that says the artifact set is complete,
    so the privacy audit runs against a finished run rather than a partial one.
    """
    input:
        VALIDATED_OK,
        NORMALISED_VARIANTS,
        QC_REPORT,
        ANNOTATED_VARIANTS,
        PHENOTYPE_PROFILE,
        CANDIDATE_PAIRS,
        RANKED_PAIRS,
        TRACK1_ARTIFACTS,
        MECHANISM_REPORT,
        DRUG_HYPOTHESES,
        REJECTED_DRUGS,
        TRACK2_ARTIFACTS,
    output:
        touch(PROVENANCE_OK),
    shell:
        "true"


rule privacy_audit:
    """Blocking gate: the repository holds no patient data (PRIV-01, ADR 0009).

    Last in the DAG on purpose. It scans the REPOSITORY, not the workspace, and
    reports paths, line numbers, span lengths and rule IDs — never the matched
    bytes, because a scanner that echoes what it finds is itself the leak
    (GP-41).
    """
    input:
        PROVENANCE_OK,
    output:
        touch(PRIVACY_AUDIT_OK),
    log:
        f"{LOG_DIR}/privacy-audit.log",
    shell:
        f"{MVA} privacy audit --repo {shlex.quote(str(REPO_ROOT))}"
        f" --workspace {shlex.quote(WORKSPACE)} > {{log}} 2>&1"
