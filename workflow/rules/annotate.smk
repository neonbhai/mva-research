# ---------------------------------------------------------------------------
# Stage 2 — annotate
# ---------------------------------------------------------------------------
# Annotation is LOCAL by construction (PRIV-05). Consequence, frequency and
# gene-phenotype tables are pre-downloaded, hash-pinned files under knowledge/;
# no proband coordinate ever leaves the machine. A structural test forbids
# importing a network client anywhere under src/mva/annotation.


rule annotate:
    """Normalised variants + phenotype -> annotated variants and a HPO profile.

    Two artifacts because they are separately hashable and separately classified:
    the annotated variants are SENSITIVE, the phenotype profile is the term set
    the ranking stage consumes.
    """
    input:
        variants=NORMALISED_VARIANTS,
        phenotype=INPUT_PHENOTYPE,
    output:
        variants=ANNOTATED_VARIANTS,
        profile=PHENOTYPE_PROFILE,
    log:
        f"{LOG_DIR}/annotate.log",
    conda:
        "../envs/mva.yaml"
    shell:
        mva_run("annotate")
