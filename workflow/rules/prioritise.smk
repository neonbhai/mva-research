# ---------------------------------------------------------------------------
# Stage 3 — prioritise
# ---------------------------------------------------------------------------
# Filtering and ranking are separate (GP-13, ADR 0005). Hard filters remove only
# invalid or impossible records; everything else is flagged and down-ranked. The
# score is a vector, not a scalar: the component scores stay visible into the
# report so "why is this ranked first?" is answerable without re-running.


rule prioritise:
    """Annotated variants -> candidate pairs, scored and ranked.

    Phase is never assumed: two heterozygous variants in one gene are a compound
    heterozygote only if in trans, and PhaseStatus.UNKNOWN survives end to end
    rather than being silently upgraded (GP-15, ASSUMPTION-PHASE-01).
    """
    input:
        variants=ANNOTATED_VARIANTS,
    output:
        ranked=RANKED_PAIRS,
    log:
        f"{LOG_DIR}/prioritise.log",
    conda:
        "../envs/mva.yaml"
    shell:
        mva_run("prioritise")
