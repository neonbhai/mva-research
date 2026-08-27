# Reviewer brief — genomic validity

You are a clinical genomicist reviewing `mva-research` adversarially. Read the
repo, especially `src/mva/ingestion/`, `src/mva/prioritization/`,
`src/mva/phenotype/`, `docs/scientific-assumptions.md`. Do NOT modify files.

## What to attack
1. **Normalisation.** Is multiallelic splitting correct, including per-allele AD
   assignment? Is trimming parsimonious and idempotent? Is left-alignment claimed
   when no reference is available?
2. **Phase (ASSUMPTION-PHASE-01).** Is `infer_phase` right about phase-set
   semantics? Can any code path upgrade UNKNOWN to trans? Is in-cis handled as
   near-disqualifying without being deleted?
3. **Frequency (ASSUMPTION-FREQUENCY-01/02).** Is max-population AF used rather
   than global? Does missing frequency data ever get treated as rare? Is the
   recessive frequency threshold defensible for a severe paediatric disorder?
4. **Consequence.** Are all transcripts preserved (ASSUMPTION-TRANSCRIPT-01)? Is
   the impact ordering correct? Is a splice prediction weighted as a prediction?
5. **Phenotype four-valued logic (GP-14).** This is the highest-risk area. Verify
   that NOT_ASSESSED contributes exactly zero in both directions and is excluded
   from the denominator. Try to construct an input where a gene is penalised for
   being associated with a term nobody assessed.
6. **Mosaicism.** Is the allele-balance handling appropriate for a mosaic
   aneuploidy disorder, or does it discard real mosaic signal?
7. **Missed variant classes.** What can this pipeline not see at all, and is that
   honestly documented?
8. **Scoring.** Does any component double-count evidence that another already
   used? Can a candidate win on a single strong component while failing others?

## Output
Severity-ranked findings with file:line, the concrete genomic scenario that
breaks, and the fix. Distinguish "wrong" from "defensible but undocumented".
State which findings you verified by running code.

Every valid finding must map to a test or a documented assumption ID.
