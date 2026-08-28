"""Regenerate the BUB1B drug catalogue TSV with guaranteed column alignment."""

from pathlib import Path

COLUMNS = (
    "drug_id",
    "name",
    "approved_name",
    "approval_status",
    "intervention_class",
    "target",
    "target_node_id",
    "mechanism_of_action",
    "observed_direction",
    "is_direct_evidence",
    "strongest_evidence_type",
    "has_pediatric_exposure",
    "youngest_age_studied",
    "pediatric_indication",
    "route",
    "cns_penetrant",
    "achievable_plasma_um",
    "required_effective_um",
    "half_life_hours",
    "worsens_cin",
    "validation_experiment",
    "notes",
)

HEADER = """\
# LITERATURE-DERIVED drug catalogue for the BUB1B / BubR1 mechanism. NOT synthetic; see SOURCES.md.
# Every agent is real and every direction cell is sourced to a named publication. NOTHING HERE IS
# MEDICAL ADVICE: each row is a research hypothesis requiring pre-clinical validation, and several
# rows exist precisely to be REJECTED, because the rejections are the most useful output.
#
# Two deliberate uniformities. (1) Every concentration cell is blank, because no publication
# reports an effective concentration at the BubR1-abundance node for any of these agents;
# ASSUMPTION-DRUG-05 therefore fires on every candidate, which is the honest state of the field.
# (2) `worsens_cin` is blank wherever the question has not actually been asked of the agent.
# Blank is unassessed, and in a chromosomal-instability disorder unassessed BLOCKS (ADR 0011,
# ASSUMPTION-DRUG-07). It is not a clean bill of health and must never be filled in optimistically.
#
# GENERATED from tools/build_knowledge/bub1b_catalog.py. Edit that file, not this one.
"""

ROWS: list[dict[str, str]] = [
    {
        "drug_id": "DRUG-NMN",
        "name": "Nicotinamide mononucleotide",
        "approval_status": "investigational_phase_2",
        "intervention_class": "disease_modifying",
        "target": "NAD+ pool feeding SIRT2-dependent deacetylation of BubR1 lysine-668",
        "target_node_id": "N3_bubr1",
        "mechanism_of_action": (
            "Raises the NAD+ pool. SIRT2 uses NAD+ to hold BubR1 lysine-668 deacetylated against "
            "CBP-dependent acetylation, which slows BubR1 turnover. NMN administration raised "
            "BubR1 abundance in mice in vivo."
        ),
        "observed_direction": "increase",
        "is_direct_evidence": "1",
        "strongest_evidence_type": "animal_model",
        "has_pediatric_exposure": "0",
        "route": "oral",
        "validation_experiment": (
            "Quantify total BubR1 by immunoblot and targeted mass spectrometry in BUB1B-MVA "
            "patient-derived fibroblasts across a clinically plausible exposure range, with "
            "intracellular NAD+ and BubR1 K668 acetylation as intermediate readouts and blinded "
            "scoring of premature chromatid separation and micronucleus frequency as the "
            "phenotype. Include a SIRT2-knockdown arm for mechanism specificity and an "
            "ectopic-BubR1 arm as the ceiling control."
        ),
        "notes": (
            "Direction AGREES at the therapeutic target on the strongest direct evidence in this "
            "catalogue. Not repurposable: there is no approved medicine. Note the seam in the "
            "evidence: NMN was shown to raise BubR1 in mice, while the lifespan rescue in the "
            "BubR1-hypomorphic genotype used SIRT2 overexpression, so the agent and the disease "
            "genotype have never been combined in one experiment. The NAD+ decline that motivates "
            "the axis is an ageing phenomenon, and whether a young child has any headroom on it "
            "is unknown."
        ),
    },
    {
        "drug_id": "DRUG-NR",
        "name": "Nicotinamide riboside",
        "approval_status": "investigational_phase_2",
        "intervention_class": "disease_modifying",
        "target": "NAD+ pool feeding SIRT2-dependent deacetylation of BubR1 lysine-668",
        "target_node_id": "N3_bubr1",
        "mechanism_of_action": (
            "Orally bioavailable NAD+ precursor. Raises blood NAD+ in humans. Its effect on BubR1 "
            "abundance has never been measured."
        ),
        "observed_direction": "unknown",
        "is_direct_evidence": "0",
        "strongest_evidence_type": "human_trial",
        "has_pediatric_exposure": "0",
        "route": "oral",
        "validation_experiment": (
            "Before any efficacy question, measure whether a clinically achievable exposure raises "
            "BubR1 protein at all in BUB1B-MVA patient cells. The human data establish that this "
            "agent raises NAD+; nothing establishes that raising NAD+ in a child raises BubR1."
        ),
        "notes": (
            "Sold as a food ingredient, not authorised as a medicine, so it is not a repurposing "
            "candidate however attractive the mechanism. Direction on the target node is UNKNOWN, "
            "which is neither agreement nor disagreement."
        ),
    },
    {
        "drug_id": "DRUG-NAM",
        "name": "Nicotinamide",
        "approved_name": "nicotinamide",
        "approval_status": "approved_other_indication",
        "intervention_class": "disease_modifying",
        "target": "NAD+ salvage pathway; also a product inhibitor of sirtuins including SIRT2",
        "target_node_id": "N3_bubr1",
        "mechanism_of_action": (
            "Feeds the NAD+ salvage pathway, which should support SIRT2-dependent deacetylation of "
            "BubR1 K668 and slow BubR1 turnover. It is simultaneously a non-competitive product "
            "inhibitor of sirtuins, which would act in the opposite direction on the same axis."
        ),
        "observed_direction": "context_dependent",
        "is_direct_evidence": "0",
        "strongest_evidence_type": "biochemical_binding",
        "has_pediatric_exposure": "1",
        "youngest_age_studied": "children (five-year randomised exposure at 1.2 g/m2 daily)",
        "pediatric_indication": "Type 1 diabetes prevention trial; pellagra",
        "route": "oral",
        "worsens_cin": "0",
        "validation_experiment": (
            "Resolve the SIGN before anything else. Dose-response measurement of intracellular "
            "NAD+, SIRT2 activity, BubR1 K668 acetylation and total BubR1 in BUB1B-MVA patient "
            "fibroblasts, across the plasma range achieved by the doses already given to children, "
            "with blinded premature-chromatid-separation scoring as the phenotypic readout."
        ),
        "notes": (
            "Direction CANNOT BE DETERMINED and must not be rounded to agreement. Nicotinamide is "
            "both an NAD+ precursor and an inhibitor of the very sirtuin the mechanism depends on, "
            "so the net sign at the BubR1 node is genuinely bidirectional. worsens_cin is recorded "
            "as assessed-and-negative on the strength of a phase 3 randomised trial with a "
            "cancer-incidence endpoint, which found fewer cancers; that trial was in adults with "
            "UV-driven skin cancer and does not transfer to a germline chromosomal-instability "
            "population."
        ),
    },
    {
        "drug_id": "DRUG-NA",
        "name": "Nicotinic acid",
        "approved_name": "nicotinic acid",
        "approval_status": "approved_other_indication",
        "intervention_class": "disease_modifying",
        "target": "NAD+ pool via the Preiss-Handler pathway",
        "target_node_id": "N3_bubr1",
        "mechanism_of_action": (
            "Raises the NAD+ pool through the Preiss-Handler pathway. In adults with mitochondrial "
            "myopathy it raised blood NAD+ several-fold and improved muscle performance. It has no "
            "sirtuin-inhibiting liability. Its effect on BubR1 abundance has never been measured."
        ),
        "observed_direction": "unknown",
        "is_direct_evidence": "0",
        "strongest_evidence_type": "human_trial",
        "has_pediatric_exposure": "1",
        "youngest_age_studied": "children",
        "pediatric_indication": "Familial hypercholesterolaemia",
        "route": "oral",
        "validation_experiment": (
            "Measure BubR1 abundance and missegregation in patient cells at the plasma NAD+ "
            "elevations already documented in humans; separately, interrogate the existing large "
            "niacin cardiovascular outcome datasets for cancer incidence, which would close the "
            "blocking oncogenic-risk gap without a new study."
        ),
        "notes": (
            "The only approved agent here with a clean mechanistic story and no known "
            "counter-direction. It is nonetheless rejected, because whether it changes aneuploidy "
            "or cancer susceptibility has never been asked. That gap is blocking in this disease "
            "context by design, and it is unusually cheap to close."
        ),
    },
    {
        "drug_id": "DRUG-ATA",
        "name": "Ataluren",
        "approval_status": "withdrawn",
        "intervention_class": "disease_modifying",
        "target": "Ribosomal readthrough of premature termination codons",
        "target_node_id": "N3_bubr1",
        "mechanism_of_action": (
            "Promotes ribosomal readthrough of premature termination codons, which would restore "
            "full-length protein from a nonsense allele. Never tested on BUB1B."
        ),
        "observed_direction": "unknown",
        "is_direct_evidence": "0",
        "strongest_evidence_type": "cell_line",
        "has_pediatric_exposure": "1",
        "youngest_age_studied": "2 years",
        "pediatric_indication": "Nonsense-mutation Duchenne muscular dystrophy",
        "route": "oral",
        "validation_experiment": (
            "Only meaningful after two prior determinations: that the patient carries a nonsense "
            "allele at all, and that the allele produces a transcript for a ribosome to read "
            "through. Assay order is transcript quantification first, readthrough second."
        ),
        "notes": (
            "Genotype-conditional and mechanistically undercut here. Readthrough needs a "
            "transcript, and truncating BUB1B alleles in MVA are reported to yield no detectable "
            "transcript, so nonsense-mediated decay removes the substrate before the ribosome "
            "reaches the stop codon. Regulatory status is also against it: the EU conditional "
            "marketing authorisation was not renewed and the European Commission decision took "
            "effect in 2025. The approval_status enum has no per-jurisdiction granularity, so "
            "`withdrawn` records the EU position and should not be read as a global statement; "
            "national authorisations elsewhere have differed."
        ),
    },
    {
        "drug_id": "DRUG-GEN",
        "name": "Gentamicin",
        "approved_name": "gentamicin",
        "approval_status": "approved",
        "intervention_class": "disease_modifying",
        "target": "Ribosomal readthrough of premature termination codons",
        "target_node_id": "N3_bubr1",
        "mechanism_of_action": (
            "Aminoglycoside that promotes readthrough of premature termination codons. "
            "Demonstrated to restore a full-length protein in a different nonsense-mutation "
            "disorder. Never tested on BUB1B."
        ),
        "observed_direction": "unknown",
        "is_direct_evidence": "0",
        "strongest_evidence_type": "human_trial",
        "has_pediatric_exposure": "1",
        "youngest_age_studied": "neonates",
        "pediatric_indication": "Serious gram-negative infection",
        "route": "intravenous",
        "cns_penetrant": "0",
        "validation_experiment": (
            "The same two determinations as ataluren, plus a chronic-exposure toxicity plan that "
            "does not currently exist for this indication."
        ),
        "notes": (
            "The same nonsense-mediated-decay objection applies. Separately, cumulative "
            "ototoxicity and nephrotoxicity make chronic paediatric dosing implausible for a "
            "lifelong condition. This catalogue schema has no field for an agent-specific "
            "toxicity, so that hazard is recorded here and in the report rather than scored, "
            "which is a limitation of the model and not an absence of risk."
        ),
    },
    {
        "drug_id": "DRUG-MPS1",
        "name": "Reversine",
        "approval_status": "tool_compound",
        "intervention_class": "disease_modifying",
        "target": "TTK/MPS1 kinase",
        "target_node_id": "N5_sac",
        "mechanism_of_action": (
            "Inhibits MPS1, the kinase that builds spindle assembly checkpoint signalling at "
            "kinetochores. Abrogates the checkpoint and causes chromosome missegregation. This is "
            "the intended effect in oncology: push chromosomally unstable tumour cells past their "
            "aneuploidy-tolerance ceiling."
        ),
        "observed_direction": "decrease",
        "is_direct_evidence": "1",
        "strongest_evidence_type": "cell_line",
        "has_pediatric_exposure": "0",
        "route": "none",
        "cns_penetrant": "0",
        "worsens_cin": "1",
        "validation_experiment": (
            "None. This agent must not be advanced for the soma. If an MPS1 inhibitor is ever "
            "considered it is as a tumour-directed agent under an entirely separate risk "
            "assessment, and the experiment is a tumour-versus-normal-tissue therapeutic-index "
            "study, not a rescue study."
        ),
        "notes": (
            "THE NAIVE TOP HIT, and the reason this pipeline exists. Maximal target proximity: it "
            "binds the kinase that builds the exact checkpoint named in the mechanism. The sign is "
            "inverted. This child is not a tumour, and the non-tumour cells already sit at the "
            "aneuploidy-tolerance ceiling with no reserve. Any ranking by pathway membership or "
            "target proximity that cannot see signs will put this compound first."
        ),
    },
    {
        "drug_id": "DRUG-HDAC",
        "name": "Vorinostat",
        "approved_name": "vorinostat",
        "approval_status": "approved_adult_only",
        "intervention_class": "disease_modifying",
        "target": "Class I and II histone deacetylases",
        "target_node_id": "N5_sac",
        "mechanism_of_action": (
            "Inhibits class I and II histone deacetylases. Measured at the mitotic readout, the "
            "class induces premature sister-chromatid separation and overrides the spindle "
            "assembly checkpoint."
        ),
        "observed_direction": "decrease",
        "is_direct_evidence": "0",
        "strongest_evidence_type": "cell_line",
        "has_pediatric_exposure": "1",
        "youngest_age_studied": "children (phase I in paediatric oncology)",
        "pediatric_indication": "Relapsed or refractory paediatric solid and CNS tumours",
        "route": "oral",
        "worsens_cin": "1",
        "validation_experiment": (
            "None as mechanism correction. If the acetylation hypothesis is to be tested at all, "
            "the experiment is a targeted one on BubR1 K250 acetylation and abundance, not a "
            "class-wide HDAC inhibitor, and it has to be read out against the checkpoint rather "
            "than against the acetylation mark."
        ),
        "notes": (
            "Here because a plausible reader would propose it: BubR1 is stabilised by acetylation "
            "at K250, and HDAC inhibitors raise acetylation, which looks like the right direction "
            "at the abundance node. Measured at the checkpoint, the class does the opposite. A "
            "worked example of a mechanism-shaped argument that reverses when someone measures the "
            "actual readout. The same class objection applies by INFERENCE to valproate, an "
            "anticonvulsant a child with this syndrome could plausibly be prescribed; that "
            "inference is untested, is not clinical advice, and does not weigh the seizure "
            "indication."
        ),
    },
    {
        "drug_id": "DRUG-CQ",
        "name": "Hydroxychloroquine",
        "approved_name": "hydroxychloroquine",
        "approval_status": "approved",
        "intervention_class": "disease_modifying",
        "target": "Lysosomal acidification and autophagic flux",
        "target_node_id": "N9_stress",
        "mechanism_of_action": (
            "Blocks lysosomal acidification and autophagic flux. Autophagy inhibition is "
            "selectively lethal to aneuploid cells, which already depend on autophagic degradation "
            "to manage their protein burden. In the soma this raises, rather than lowers, the "
            "proteotoxic stress node."
        ),
        "observed_direction": "increase",
        "is_direct_evidence": "0",
        "strongest_evidence_type": "cell_line",
        "has_pediatric_exposure": "1",
        "youngest_age_studied": "infants",
        "pediatric_indication": "Juvenile idiopathic arthritis; malaria",
        "route": "oral",
        "validation_experiment": (
            "None for the soma. As a tumour-directed agent it would need a therapeutic-index study "
            "against the patient's own non-tumour cells, because those cells share the aneuploidy "
            "that makes the tumour sensitive."
        ),
        "notes": (
            "THE THERAPEUTIC-WINDOW INVERSION, made mechanical. Autophagy inhibition is "
            "selectively lethal to aneuploid cells, which is exactly why it is attractive against "
            "the rhabdomyosarcoma and exactly why it is dangerous everywhere else: the same "
            "property that kills the tumour depletes a child whose every tissue is mosaically "
            "aneuploid and already growth-restricted. Rejected for the soma on the same evidence "
            "that recommends it for the tumour. Chloroquine, not hydroxychloroquine, was the "
            "screened compound, so the direction is carried across a congener."
        ),
    },
    {
        "drug_id": "DRUG-HSP90",
        "name": "Tanespimycin",
        "approval_status": "investigational_phase_2",
        "intervention_class": "disease_modifying",
        "target": "HSP90 chaperone",
        "target_node_id": "N9_stress",
        "mechanism_of_action": (
            "Inhibits HSP90-dependent protein folding. Aneuploid cells already have impaired "
            "HSP90-dependent folding and are selectively killed by this class, which means the "
            "agent raises the proteotoxic stress node rather than lowering it."
        ),
        "observed_direction": "increase",
        "is_direct_evidence": "1",
        "strongest_evidence_type": "cell_line",
        "has_pediatric_exposure": "0",
        "route": "intravenous",
        "validation_experiment": (
            "None for the soma. As a tumour-directed agent it would need a therapeutic-index study "
            "against the patient's own non-tumour cells, because those cells share the aneuploidy "
            "that makes the tumour sensitive."
        ),
        "notes": (
            "Second aneuploidy-selective antiproliferative from the same screen as "
            "hydroxychloroquine, and the same tumour-versus-soma inversion. There is a "
            "BUB1B-specific reason to be more concerned than for a generic cytotoxic: MVA missense "
            "BubR1 is degraded faster than wild type, so impairing chaperone-dependent folding "
            "could plausibly lower the very quantity this mechanism needs raised. That is an "
            "INFERENCE and not a measurement; BubR1 is not an established HSP90 client."
        ),
    },
    {
        "drug_id": "DRUG-SEN",
        "name": "Dasatinib",
        "approved_name": "dasatinib",
        "approval_status": "approved",
        "intervention_class": "disease_modifying",
        "target": "SRC and ABL kinases; senescent-cell anti-apoptotic pathways",
        "target_node_id": "N10_clearance",
        "mechanism_of_action": (
            "Used intermittently with quercetin as a senolytic, clearing senescent cells. Genetic "
            "clearance of p16-positive senescent cells delayed age-related deterioration in "
            "BubR1-hypomorphic mice, the mouse model of this gene."
        ),
        "observed_direction": "decrease",
        "is_direct_evidence": "0",
        "strongest_evidence_type": "animal_model",
        "has_pediatric_exposure": "1",
        "youngest_age_studied": "under 18 years (phase II, n=113)",
        "pediatric_indication": "Philadelphia-chromosome-positive chronic myeloid leukaemia",
        "route": "oral",
        "cns_penetrant": "0",
        "validation_experiment": (
            "Before any efficacy question, determine whether senescence in this disorder is on "
            "balance a burden or a barrier: measure tumour incidence, not only healthspan, in "
            "BubR1-hypomorphic mice under senolytic dosing, powered for the cancer endpoint. A "
            "healthspan-only readout cannot answer the question that matters in a "
            "cancer-predisposition syndrome."
        ),
        "notes": (
            "The best-motivated consequence-directed hypothesis here, and the clearest "
            "illustration of why a deviation from wild type is not the same as a pathology. The "
            "p16 and p53 response that drives the progeroid phenotype in BubR1-hypomorphic mice is "
            "also what suppresses aneuploidy-driven tumorigenesis. Removing it in a "
            "cancer-predisposed child may remove a barrier rather than a burden. The node is "
            "therefore flagged compensatory, no corrective direction is derived from its state, "
            "and the direction check returns cannot-determine rather than guessing. One further "
            "signal the schema has no field for: the paediatric CML trial reported bone growth "
            "and development events in 4 percent of patients, which is a specific concern in a "
            "child whose defining feature is growth restriction."
        ),
    },
    {
        "drug_id": "DRUG-METF",
        "name": "Metformin",
        "approved_name": "metformin",
        "approval_status": "approved",
        "intervention_class": "disease_modifying",
        "target": "AMP-activated protein kinase; mitochondrial complex I",
        "target_node_id": "N9_stress",
        "mechanism_of_action": (
            "Activates AMPK and imposes an energy-stress state. Energy stress induced by AICAR is "
            "selectively lethal to aneuploid cells, so the direction at the stress node is upward."
        ),
        "observed_direction": "increase",
        "is_direct_evidence": "0",
        "strongest_evidence_type": "cell_line",
        "has_pediatric_exposure": "1",
        "youngest_age_studied": "10 years",
        "pediatric_indication": "Type 2 diabetes",
        "route": "oral",
        "validation_experiment": (
            "None as mechanism correction. If the aneuploidy-selective-lethality logic is to be "
            "pursued at all it belongs in the tumour-directed arm, with a therapeutic-index "
            "readout against the patient's own non-tumour cells."
        ),
        "notes": (
            "An approved, cheap, extremely well-tolerated paediatric drug that a repurposing "
            "search would surface as obviously safe, and would be wrong to. AICAR, not metformin, "
            "was the screened compound, so the direction is carried across the AMPK-activator "
            "class rather than measured for this agent."
        ),
    },
    {
        "drug_id": "DRUG-RAPA",
        "name": "Sirolimus",
        "approved_name": "sirolimus",
        "approval_status": "approved",
        "intervention_class": "disease_modifying",
        "target": "mTORC1",
        "target_node_id": "N9_stress",
        "mechanism_of_action": (
            "Inhibits mTORC1, lowering translational load and raising autophagic flux, which is "
            "the obvious proteostasis move against aneuploidy-associated stress. Its net direction "
            "on that stress node in an aneuploid background has not been measured."
        ),
        "observed_direction": "unknown",
        "is_direct_evidence": "0",
        "strongest_evidence_type": "cell_line",
        "has_pediatric_exposure": "1",
        "youngest_age_studied": "infants",
        "pediatric_indication": "Complicated lymphatic anomalies",
        "route": "oral",
        "validation_experiment": (
            "Measure proteotoxic-stress markers and viability in BUB1B-MVA patient cells under "
            "mTOR inhibition before anything else, because the sign is genuinely unknown: reducing "
            "translational load could relieve the burden, or removing the mTOR-dependent survival "
            "programme could kill the cells that are carrying it."
        ),
        "notes": (
            "Two objections this catalogue cannot express as fields, and which the report "
            "therefore states in prose. It is anti-growth in a child whose defining feature is "
            "growth restriction, and it is immunosuppressive in a child whose aneuploid cells are "
            "partly held in check by immune clearance. Both act at nodes other than the one this "
            "agent targets, which is a general limitation of scoring a drug against a single node."
        ),
    },
]


def render() -> str:
    lines = [HEADER.rstrip("\n"), "\t".join(COLUMNS)]
    for row in ROWS:
        cells = []
        for column in COLUMNS:
            value = row.get(column, "")
            if "\t" in value or "\n" in value:
                msg = f"{row['drug_id']}: cell {column!r} contains a tab or newline"
                raise ValueError(msg)
            cells.append(value)
        lines.append("\t".join(cells))
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    import sys

    Path(sys.argv[1]).write_text(render(), encoding="utf-8")
    print(  # noqa: T201 -- this is the generator's CLI report
        f"wrote {len(ROWS)} rows x {len(COLUMNS)} columns to {sys.argv[1]}"
    )
