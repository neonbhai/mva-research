# ADR 0025 — The Track 2 target node is BubR1 abundance, not the checkpoint

**Status:** accepted · **Date:** 2026-08-28

## Context

`MechanismHypothesis.therapeutic_target_node_id` names one node on the chain. Every
drug's direction is checked against that node's required correction, so the choice
decides which agents can be accepted at all. For BUB1B-MVA there were two candidates.

**The checkpoint node** (`N5_sac`, spindle assembly checkpoint signalling) is the
obvious one. It is what the disease is named after, it is where the pathology is most
vividly described, and it is the node a reader's eye goes to.

**The protein-abundance node** (`N3_bubr1`, BubR1 protein abundance) sits two rungs
upstream and is where the patient's deviation is a *quantity*.

## Decision

The therapeutic target is `N3_bubr1`. Three reasons, in order of weight.

**1. The checkpoint node has no direction-correct pharmacological handle, and the
compounds that bind it all push the wrong way.** Everything developed against
MPS1/TTK, Aurora B or the checkpoint kinases is an *inhibitor*, because the oncology
programme wants chromosomally unstable tumour cells pushed past their
aneuploidy-tolerance ceiling. BubR1 itself is a **pseudokinase** (Suijkerbuijk et al.,
*Dev Cell* 2012, PMID 22698286), so there is no catalytic activity to agonise even in
principle. Designating the checkpoint as the target creates a target node whose entire
accessible pharmacology is contraindicated, and then measures every candidate against
it. That is not a target; it is a trap.

**2. BubR1 reaches missegregation by two independent routes, and only the protein is
upstream of both.** The checkpoint arm restrains APC/C-Cdc20. The second arm recruits
PP2A-B56 to kinetochores through the KARD domain to stabilise microtubule attachments
(Suijkerbuijk et al., *Dev Cell* 2012, PMID 23079597; Kruse et al., *J Cell Sci* 2013,
PMID 23345399; Lampson & Kapoor, *Nat Cell Biol* 2005, PMID 15592459). MVA patient
cells show *both* an impaired checkpoint and chromosome alignment defects
(Suijkerbuijk et al., *Cancer Res* 2010, PMID 20516114). An intervention at the
checkpoint node addresses one arm and leaves the other untouched. The protein node is
the last point at which one quantity governs both.

**3. The deviation at the protein node is quantitative, graded, and demonstrated to
be reversible in the right direction.** Residual BubR1 in MVA is reduced, not absent.
Graded knockdown places the emergence of aneuploidy near 50% residual protein (Bohers
et al., *Hum Genet* 2008, PMID 18932004), so there is a threshold to climb back over.
Ectopic BubR1 expression restored checkpoint activity **in MVA patient cells**
(Suijkerbuijk 2010), and sustained elevated BubR1 in mice reduced aneuploidy, reduced
tumorigenesis and extended healthy lifespan (Baker et al., *Nat Cell Biol* 2013,
PMID 23242215). MVA can even be caused by a purely regulatory allele that only lowers
expression (Ochiai et al., *PNAS* 2014, PMID 24344301). This is a dosage disease, and
dosage is the thing a drug could move.

The generalised rule is written up as ASSUMPTION-MECHANISM-06: **choose the target at
the node where the deviation is a movable quantity, not at the node that best names
the disease.**

## Consequences

- `disease_direction = decrease`, `required_correction = increase`, and the target
  node's `state_in_patient = decrease`. The `MechanismHypothesis` validator checks all
  three agree at construction (ASSUMPTION-MECHANISM-05).
- Agents acting at the checkpoint node are still evaluated — they are not filtered out
  (GP-13) — and are scored against that node's own corrective direction (`restore`).
  Reversine and vorinostat are consequently rejected with `WRONG_DIRECTION` rather
  than being silently absent.
- Agents acting on downstream consequence nodes carry `MECHANISM_MISMATCH` as a
  non-fatal concern, printed beside them (ASSUMPTION-DRUG-04).
- The choice is not free of cost. Raising BubR1 has an upper bound nobody has
  characterised: Mad2 overexpression is tumorigenic in mice (Sotillo et al., *Cancer
  Cell* 2007;11:9-23, PMID 17189715), so "more checkpoint protein" is not
  monotonically safe,
  even though BubR1 overexpression specifically was protective in Baker 2013. Any
  dose-finding work has to look for a ceiling, not only a floor.

## Alternative rejected

**Target the transcript node (`N2_transcript`) instead.** Also quantitative, and
arguably closer to the lesion. Rejected because two of the three published mechanisms
of BubR1 loss in MVA do not act on transcript level at all — accelerated turnover of
the missense protein is post-translational (Suijkerbuijk 2010) — so a transcript-level
target would be blind to the allele class the pharmacology is most likely to reach.
