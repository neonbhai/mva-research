# ---------------------------------------------------------------------------
# Track 2 — mechanism and intervention
# ---------------------------------------------------------------------------
# The second exit from the one pipeline (docs/architecture.md). Track 2 starts
# from the ranked pair rather than from a fresh analysis, because a drug
# hypothesis that is not anchored to the ranked variant pair is not grounded in
# anything.


rule mechanism:
    """Ranked pair -> a graded mechanistic chain.

    Every link in the chain is graded individually (ASSUMPTION-MECHANISM-01); the
    chain is never asserted at the strength of its strongest link, and an
    inferred link is marked as inferred. Direction of effect is mandatory and
    signed, and "unknown" is a third state that never counts as agreement
    (GP-16).
    """
    input:
        ranked=RANKED_PAIRS,
    output:
        mechanism=MECHANISM_REPORT,
    log:
        f"{LOG_DIR}/mechanism.log",
    conda:
        "../envs/mva.yaml"
    shell:
        mva_run("mechanism")


rule drugs:
    """Mechanism -> drug hypotheses, accepted and rejected.

    The rejection record is an output, not a side effect: contradicting evidence
    and failed candidates are persisted with their reasons (GP-19). A drug whose
    observed direction disagrees with the required correction cannot be
    constructed as an accepted hypothesis at all (ASSUMPTION-DRUG-01), so the
    rejections are where most of the reasoning is visible.
    """
    input:
        mechanism=MECHANISM_REPORT,
    output:
        accepted=DRUG_HYPOTHESES,
        rejected=REJECTION_RECORD,
    log:
        f"{LOG_DIR}/drugs.log",
    conda:
        "../envs/mva.yaml"
    shell:
        mva_run("drugs")
