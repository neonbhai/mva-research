# ---------------------------------------------------------------------------
# Stage 1 — validate and ingest
# ---------------------------------------------------------------------------
# Thin rules only (ADR 0001): each one is a single `mva` subcommand. All output
# paths are rooted at WORKSPACE (ADR 0006); see the header of the Snakefile.


rule validate:
    """Fail before the first patient byte is read.

    Config validation, the workspace boundary check and the input inventory all
    happen here, so a misconfigured run stops while it has still written nothing
    that matters. Modelled as a touched marker because a check produces a
    verdict, not an artifact — the run manifest it opens is asserted by the
    `provenance` rule at the other end of the DAG.
    """
    input:
        case_config=CASE_CONFIG,
        defaults=DEFAULTS_CONFIG,
        vcf=VCF_INPUT,
        phenotype=PHENOTYPE_INPUT,
    output:
        touch(VALIDATED_OK),
    log:
        f"{LOG_DIR}/validate.log",
    shell:
        mva_run("validate")


rule ingest:
    """VCF -> normalised, QC-annotated variant records.

    Reads, normalises (multiallelic split, trim, left-align) and applies
    analytical-validity checks. Nothing is deleted here except records that are
    invalid or impossible; weak calls are flagged and survive into ranking
    (GP-13, ADR 0005), which is why the QC report is an artifact in its own right
    rather than a log line.
    """
    input:
        validated=VALIDATED_OK,
        vcf=VCF_INPUT,
    output:
        variants=NORMALISED_VARIANTS,
        qc=QC_REPORT,
    log:
        f"{LOG_DIR}/ingest.log",
    conda:
        "../envs/mva.yaml"
    shell:
        mva_run("ingest")
