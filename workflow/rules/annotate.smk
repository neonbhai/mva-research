# ---------------------------------------------------------------------------
# Stage 2 — annotate
# ---------------------------------------------------------------------------
# Annotation is LOCAL by construction (PRIV-05). Consequence, frequency and
# gene-phenotype tables are pre-downloaded, hash-pinned files under knowledge/;
# no proband coordinate ever leaves the machine. A structural test forbids
# importing a network client anywhere under src/mva/annotation.


rule annotate:
    """Normalised variants + phenotype -> annotated variant records.

    Consequences and frequencies are LISTS on the record, not a single chosen
    transcript (ASSUMPTION-TRANSCRIPT-01), and a variant with no frequency data
    is recorded as having none rather than as rare (GP-14). The phenotype profile
    is built in the same stage and travels in the evidence store.
    """
    input:
        variants=NORMALISED_VARIANTS,
        phenotype=INPUT_PHENOTYPE,
    output:
        variants=ANNOTATED_VARIANTS,
    log:
        f"{LOG_DIR}/annotate.log",
    conda:
        "../envs/mva.yaml"
    shell:
        mva_run("annotate")
