# Reviewer brief — Track 2 pharmacology

You are a pharmacologist reviewing the drug-repurposing reasoning in
`mva-research` adversarially. Read `src/mva/mechanisms/`,
`src/mva/interventions/`, `knowledge/public/drug_catalog.tsv`,
`knowledge/public/mechanism_meta.tsv`, `docs/scientific-assumptions.md`
(ASSUMPTION-DRUG-01..07). Do NOT modify files.

## The central claim to attack
The pipeline claims it cannot recommend a compound that acts on the right target
in the wrong direction. **Try to break that.** Construct a catalogue entry or a
mechanism that gets a contraindicated agent through. Specifically:
- Can an unsigned direction (`UNKNOWN` / `CONTEXT_DEPENDENT`) be laundered into
  agreement anywhere?
- Can a scoring weight change resurrect a wrong-direction drug? (The claim is
  that a model validator prevents this — verify it.)
- Is `required_correction` derived correctly from the mechanism's
  `disease_direction`, or can a sign error invert the whole check?
- Can a drug acting on a node that is downstream of the real defect appear to
  correct it?

## Also examine
1. **Evidence tiering.** Is binding-vs-cells-vs-animals-vs-humans ranked honestly?
   Is indirect evidence penalised relative to direct?
2. **Concentration.** Is `concentration_achievable` used, and is `None` handled as
   a gap rather than a pass (ASSUMPTION-DRUG-05)?
3. **Paediatric evidence.** Is the caveat that oncology-population tolerability
   does not transfer to a germline-instability population actually applied, or
   just stated?
4. **Symptomatic vs disease-modifying (ASSUMPTION-DRUG-04).** Can a symptomatic
   agent be presented as mechanism correction anywhere in the rendered output?
5. **Approval status (ASSUMPTION-DRUG-03).** Is a tool compound ever described as
   a repurposing candidate?
6. **Oncogenic risk.** Is `worsens_chromosomal_instability = None` treated as a
   blocking gap in this disease context, or silently as fine?
7. **Developmental window (ASSUMPTION-MECHANISM-02).** Does the report surface
   that post-natal correction cannot reverse established structural findings, or
   is it buried?
8. **The rejection record.** Are rejections preserved with reasons (GP-19), and
   are the reasons correctly distinguished (wrong-direction vs not-approved vs
   target-not-in-mechanism)?

## Output
Severity-ranked findings, file:line, the concrete scenario, the fix. Be
especially harsh about anything that could read as clinical recommendation.
Every valid finding must map to a test.
