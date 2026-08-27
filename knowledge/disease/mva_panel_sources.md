# MVA panel — sources, method, and every place the curated record disagrees

Companion to `knowledge/disease/mva_panel.tsv`. Panel version `2026-08-28-a`, **116 genes**.

This file exists so a reviewer can audit every claim in the panel without trusting
the author. **Nothing in the TSV came from model recall.** Every gene–disease link,
every OMIM/MONDO/HGNC/Ensembl identifier and every PMID was read out of a curated
file on this machine by a script, and the script refuses to emit an identifier it
cannot find in one of those files. Where recall and the curated record disagreed,
the curated record won and the disagreement is written down below rather than
silently resolved.

## 1. Sources of record

| Source | File | Snapshot | Used for |
|---|---|---|---|
| EBI Gene2Phenotype (DDG2P) | `mva-resources/genepanels/DDG2P.csv` | downloaded 2026-08-27, 2 869 records | gene–disease link, disease name, disease OMIM + MONDO, allelic requirement, confidence, molecular mechanism, **publication PMIDs**, gene MIM, previous symbols |
| ClinGen Gene–Disease Validity | `mva-resources/genepanels/clingen_gene_validity.csv` | file created 2026-08-27, 3 668 curations | independent classification, MOI, MONDO, expert panel, classification date |
| HPO annotations | `mva-resources/hpo/phenotype.hpoa` | v2026-06-23 (8 574 OMIM, 4 337 ORPHANET) | disease names, per-disease **reference PMIDs**, curated inheritance terms, phenotype profile |
| HPO gene↔phenotype | `mva-resources/hpo/phenotype_to_genes.txt`, `genes_to_phenotype.txt` | v2026-06-23 | gene→disease links absent from DDG2P/ClinGen (this is how MVA4, MVA7 and the Orphanet MVA set were found) |
| HPO ontology | `mva-resources/hpo/hp.obo` | v2026-06-23 | every `HP:` id in §4 checked for existence and non-obsolescence |
| MANE | `mva-resources/mane/MANE.GRCh38.v1.5.summary.txt.gz` | v1.5, 19 363 MANE Select genes | **HGNC id, Ensembl gene id, NCBI GeneID, GRCh38 contig, HGNC-approved current symbol** |
| ClinVar | `mva-resources/clinvar/clinvar.vcf.gz` | `fileDate=2026-08-22`, GRCh38, 4 467 990 records | `clinvar_plp_n` = records with `CLNSIG` in {Pathogenic, Likely_pathogenic, Pathogenic/Likely_pathogenic, ±low_penetrance} whose `GENEINFO` names the gene |

`DDG2P.csv.gz` in the same directory is **not** a gzip file — it is an HTML error page
from a failed download (1 120 bytes). The uncompressed `DDG2P.csv` beside it is the
real one and is what was used. Anyone re-running this should ignore the `.gz`.

## 2. Method, and the checks that make the citations trustworthy

1. A seed list names each gene, its panel tier, its mechanism class, and *optionally*
   the OMIM disease id the curator believes is the relevant one.
2. The builder joins that seed against all six sources above.
3. **Every seeded OMIM id is verified twice before it is allowed into the file:**
   it must (a) exist in DDG2P's disease-MIM column or in `phenotype.hpoa`, and
   (b) be linked to *that gene* by DDG2P or by `phenotype_to_genes.txt`. An id that
   fails either check is blanked and reported.
4. Citations are only ever copied — from DDG2P's `publications` column, from
   `phenotype.hpoa`'s `reference` column, or from an ORPHA/OMIM disease id. They are
   never composed.
5. Symbols, HGNC ids, Ensembl gene ids and contigs come from MANE Select only. A gene
   absent from MANE Select would have been reported and dropped; none were.

**The verification caught a real error.** The seed proposed `OMIM:614928` for
*NSMCE2* from recall. The check found that this id exists in HPOA but names
*Ectodermal dysplasia 6, hair/nail type* and is not linked to *NSMCE2*. The correct
id, read from `phenotype_to_genes.txt`, is `OMIM:617253` *Seckel syndrome 10*
(PMID:25105364). This is exactly the failure mode a fabricated citation produces, and
it is why the check exists rather than a proofread.

Reproduce: the builder and seed used are in the session scratchpad; the joins are
five `csv`/`awk` reads and are re-derivable from §1 in a few minutes.

## 3. Panel composition

| | count |
|---|---|
| Total genes | 116 |
| `panel_tier` 1 (MVA phenotype series + Orphanet MVA gene set) | 7 |
| `panel_tier` 2 (CIN / repair / centrosome / microcephaly / growth differential) | 109 |
| `evidence_tier` definitive | 52 |
| `evidence_tier` strong | 22 |
| `evidence_tier` moderate | 11 |
| `evidence_tier` limited | 20 |
| `evidence_tier` **candidate** (no curated gene–disease validity anywhere) | 11 |
| Recessive-compatible (`autosomal_recessive` or `x_linked`) | 92 |
| Absent from DDG2P | 29 |
| Absent from ClinGen | 52 |

### Tier 1 — the genes this case is actually about

| Gene | Disease | OMIM | MONDO | Inheritance | Evidence tier | Citation | ClinVar P/LP |
|---|---|---|---|---|---|---|---|
| **BUB1B** | BUB1B-related mosaic variegated aneuploidy syndrome (MVA1) | OMIM:257300 | MONDO:0009759 | AR, biallelic | **definitive** (DDG2P definitive + ClinGen Definitive 2019-11-22) | PMID:9916837, 11169558, 15475955, 16411201, 21190457 | 81 |
| **CEP57** | CEP57-related mosaic variegated aneuploidy syndrome (MVA2) | OMIM:614114 | MONDO:0013582 | AR, biallelic | **definitive** (DDG2P definitive + ClinGen Definitive 2019-11-22) | PMID:12116237, 21552266, 24259107 | 19 |
| **TRIP13** | TRIP13-related mosaic variegated aneuploidy and Wilms tumour (MVA3) | OMIM:617598 | MONDO:0054736 | AR, biallelic | **strong** (DDG2P strong; **no ClinGen curation exists**) | PMID:28553959 | 9 |
| **CENATAC** | Mosaic variegated aneuploidy syndrome 4 | OMIM:620153 | — | AR (HPOA HP:0000007) | limited | PMID:34009673 | 2 |
| **MAD1L1** | Mosaic variegated aneuploidy syndrome 7 with inflammation and tumour predisposition | OMIM:620189 | — | AR (HPOA HP:0000007) | limited | PMID:36322655 | 8 |
| **BUB1** | BUB1-related microcephaly and developmental disorder (MCPH30); also in the Orphanet MVA gene set | OMIM:620183 | MONDO:0859342 | AR, biallelic | moderate (DDG2P) | PMID:35044816 | 6 |
| **BUB3** | Mosaic variegated aneuploidy syndrome (Orphanet gene set only) | — | — | not curated | limited | ORPHA:1052 | **0** |

Tier 1 is 5 genes with an OMIM MVA phenotype-series entry plus the 2 extra genes
(*BUB1*, *BUB3*) that Orphanet places under `ORPHA:1052 Mosaic variegated aneuploidy
syndrome`. The full Orphanet MVA gene set is exactly {BUB1, BUB1B, BUB3, CEP57,
TRIP13}, read from `phenotype_to_genes.txt`.

`BUB1B` additionally carries `OMIM:176430 Premature chromatid separation trait`,
which is annotated in HPOA as autosomal **dominant** with the single phenotype
HP:0200024 — a heterozygous-carrier trait, not the syndrome. Do not confuse the two
when a single *BUB1B* het turns up.

## 4. HPO terms for this phenotype — all verified against `hp.obo` v2026-06-23

Every id below exists in the downloaded ontology and none is obsolete. The right-hand
column is which MVA phenotype-series entries HPOA actually annotates the term to
(`-` = the term is real and relevant but is not currently annotated to any MVA entry).

### The cytogenetic hallmark

| HP id | Label | Annotated to |
|---|---|---|
| `HP:0200024` | Premature chromatid separation | MVA1, MVA3 |
| `HP:0003220` | Abnormality of chromosome stability | MVA4 |
| `HP:0040012` | Chromosome breakage | — |
| `HP:0003221` | Chromosomal breakage induced by crosslinking agents | — |
| `HP:0010997` | Chromosomal breakage induced by ionizing radiation | — |
| `HP:0010998` | Increased susceptibility to spontaneous sister chromatid exchange | — |

> **The HPO has no term for aneuploidy, mosaic aneuploidy, or mosaic variegated
> aneuploidy.** Searching `hp.obo` for `aneuploid` returns nothing (only unrelated
> `aneurysm` string matches). The disease's defining cytogenetic finding is
> represented as `HP:0200024` plus the breakage terms above. Any phenotype-matching
> code that expects to find an "aneuploidy" HPO term for this proband will silently
> match nothing. `HP:0200024` is the term to key on, and note that HPOA does **not**
> annotate it to MVA2 (*CEP57*) or MVA4 (*CENATAC*) — absence there is an annotation
> gap, not evidence that those patients lack PCS.

### Growth and neurodevelopment

| HP id | Label | Annotated to |
|---|---|---|
| `HP:0000252` | Microcephaly | MVA1, MVA2, MVA3, MVA4, MVA7 |
| `HP:0011451` | Primary microcephaly | — |
| `HP:0001511` | Intrauterine growth retardation | MVA1, MVA2 |
| `HP:0001518` | Small for gestational age | MVA1, MVA2 |
| `HP:0008846` | Severe intrauterine growth retardation | MVA2 |
| `HP:0001510` | Growth delay | MVA2, MVA3 |
| `HP:0008897` | Postnatal growth retardation | MVA1 |
| `HP:0004322` | Short stature | MVA1, MVA2, MVA3 |
| `HP:0030674` | Antenatal onset | MVA1 |
| `HP:0001263` | Global developmental delay | MVA3, MVA7 |
| `HP:0011344` | Severe global developmental delay | MVA1 |
| `HP:0002187` | Profound intellectual disability | MVA1 |
| `HP:0001250` | Seizure | MVA1, MVA3 |
| `HP:0001290` | Generalized hypotonia | MVA1 |

### Brain structure

| HP id | Label | Annotated to |
|---|---|---|
| `HP:0001305` | Dandy-Walker malformation | MVA1 |
| `HP:0001321` | Cerebellar hypoplasia | MVA1 |
| `HP:0001274` | Agenesis of corpus callosum | MVA1 |
| `HP:0002119` | Ventriculomegaly | MVA1 |
| `HP:0000238` | Hydrocephalus | MVA1 |

### Tumour predisposition — the discriminating feature

| HP id | Label | Annotated to |
|---|---|---|
| `HP:0002667` | Nephroblastoma (Wilms tumour) | MVA1, MVA3 |
| `HP:0006743` | Embryonal rhabdomyosarcoma | MVA1, MVA7 |
| `HP:0001909` | Leukemia | MVA1 |
| `HP:0002664` | Neoplasm | — |

`HP:0002667` is annotated to MVA3 at frequency **6/6** — every reported *TRIP13*
patient in PMID:28553959 had a Wilms tumour. If the proband has a nephroblastoma,
*TRIP13* and *BUB1B* move ahead of *CEP57*, which has no tumour annotation at all.

### Other MVA features worth carrying

| HP id | Label | Annotated to |
|---|---|---|
| `HP:0005387` | Combined immunodeficiency | MVA1 |
| `HP:0000518` | Cataract | MVA1 |
| `HP:0000639` | Nystagmus | MVA1, MVA3, MVA7 |
| `HP:0000365` | Hearing impairment | MVA2 |
| `HP:0000821` | Hypothyroidism | MVA2 |

### Inheritance terms used to fill `inheritance` where DDG2P/ClinGen were silent

`HP:0000007` autosomal recessive · `HP:0000006` autosomal dominant ·
`HP:0001417`/`HP:0001419`/`HP:0001423` X-linked · `HP:0001427` mitochondrial.

## 5. Where a curated source contradicts common recall — read this before ranking

Each of these is a place where the panel deliberately does **not** say what a
textbook or a model would say. `evidence_tier` in the TSV follows the curated
source; the `notes` column repeats the short version.

### 5.1 Two curated sources contradict each other

| Gene | DDG2P | ClinGen | What the panel does |
|---|---|---|---|
| **FANCM** | `FANCM-related Fanconi anemia`, **strong** — with *no PMIDs and no OMIM id in the record* | **Refuted** for Fanconi anaemia (MONDO:0019391, 2024-12-10) **and** Definitive for a *separate* entity, `FANCM Fanconi-like genomic instability disorder` (MONDO:0100578, 2024-11-12) | Kept, `evidence_tier=strong` from DDG2P, with the conflict in `notes`. **A biallelic *FANCM* hit must not be reported as Fanconi anaemia.** ClinGen is the more recent and better-evidenced call here. |
| **CDC6** | Meier-Gorlin syndrome 5, **definitive** | **Limited** (2025-01-17, Syndromic Disorders GCEP) | Kept at DDG2P's definitive per the panel's "strongest classification" rule, conflict flagged. Treat as *limited* when ranking. |
| **ORC4** | Meier-Gorlin syndrome 2, **definitive** | **Moderate** (2025-06-30) | Same treatment. |
| **RAD51C** | FA complementation group O, **strong** | **Limited** for FA-O (2023-12-20) — while Definitive for the *dominant* cancer-predisposition phenotype | Same treatment. The recessive FA-O claim is the weak one. |
| **RECQL4** | disease name recorded as **Baller-Gerold syndrome**, OMIM:218600 | **Rothmund-Thomson syndrome** (MONDO:0010002), Definitive | The TSV carries DDG2P's Baller-Gerold label because that is the row that was joined. *RECQL4* causes Rothmund-Thomson, RAPADILINO **and** Baller-Gerold; do not read the single label as the full allelic series. |
| **DNA2** | microcephalic primordial dwarfism ±poikiloderma/cataracts, **limited** | **Moderate**, but for *mitochondrial disease*, **AD** (2026-02-09) | Panel keeps the recessive dwarfism link at `limited`. The two sources are describing different diseases. |
| **NSMCE3** | distinct DNA breakage syndrome, **limited** | **Limited** (2026-01-15) | Consistent; both weak. |

### 5.2 A curated source contradicts what a clinician or model would assume

| Gene | Common assumption | What the curated record says |
|---|---|---|
| **CENPE** | An established autosomal-recessive primary microcephaly gene (MCPH13, OMIM:616051) | ClinGen classifies it **Limited** for AR primary microcephaly (2023-12-19, Brain Malformations GCEP), and it is **absent from DDG2P entirely**. Only 5 ClinVar P/LP records. Treat a *CENPE* biallelic hit as a hypothesis needing functional work, not an answer. |
| **NBN** | A hereditary breast-cancer gene | ClinGen **Refuted** the AD breast-carcinoma link (2023-03-14) while keeping AR Nijmegen breakage syndrome **Definitive**. Only the recessive claim is supportable. |
| **BRIP1**, **SLX4**, **XRCC2**, **RAD50** | Breast-cancer predisposition genes | All four are **Refuted** for hereditary breast carcinoma by ClinGen. Their real paediatric relevance is recessive (FA-J, FA-P, FA-U, NBS-like). |
| **RAD51** | Fanconi anaemia complementation group R | HPO links *RAD51* to `OMIM:617244 FA-R` — but the mechanism is **monoallelic dominant-negative**, not biallelic. DDG2P curates only a *limited* monoallelic mirror-movements phenotype. **A compound-het model does not apply to this gene.** |
| **GMNN** | Recessive Meier-Gorlin | DDG2P: `monoallelic_autosomal`, **gain of function**. Not a compound-het candidate. |
| **POLE** | Dominant polymerase-proofreading cancer syndrome | DDG2P curates a *separate* **biallelic** IMAGE-I syndrome with immunodeficiency (OMIM:618336, strong). Both are real; only the biallelic one fits this case. |
| **NIPBL**, **SMC3**, **RAD21**, **WT1**, **KIF11**, **TUBB** | Panel genes for a recessive workup | All curated **monoallelic/dominant** and usually *de novo*. They are on the panel as differentials, and a single het in them is a legitimate finding — but it is not a compound-het answer. `inheritance` in the TSV is the field to filter on. |
| **CENPJ** | The MCPH6/Seckel-4 gene symbol | HGNC-approved current symbol is **CPAP** (HGNC:17272, ENSG00000151849). MANE and DDG2P both use `CPAP`; ClinGen and most annotation tools still emit `CENPJ`. ClinVar has **zero** records under `CENPJ`. The alias is in `gene_symbol_aliases` — join on HGNC or Ensembl id, not on symbol. |
| **KNL1** | `CASC5` | Renamed. Alias carried. |
| **ATRIP** | Seckel syndrome gene | Orphanet lists it under `ORPHA:808`, but the **only ClinGen curation is a Limited AD breast-carcinoma link** and it is absent from DDG2P. Its ClinVar P/LP count (81) is inflated: *ATRIP* shares its chr3 locus with *TREX1*, so `GENEINFO` names both on the same records. |

### 5.3 Genes on the panel with no curated gene–disease validity at all

Eleven genes carry `evidence_tier=candidate` and an **empty citation field**, which
is deliberate: no PMID, OMIM id or ORPHA id for them exists in any source in §1.

`TTK` · `AURKB` · `MAD2L1` · `ESPL1` · `ZW10` · `ZWILCH` · `KNTC1` · `NDC80` ·
`NUF2` · `SPC24` · `SPC25`

They are here on mechanism alone — they are the core mitotic-checkpoint kinase, the
chromosomal passenger complex, the mitotic checkpoint complex, separase, the RZZ
complex and the NDC80 kinetochore–microtubule attachment complex, i.e. the machinery
whose failure produces exactly the missegregation phenotype this proband has. Their
purpose is to stop a novel biallelic hit in the pathway being invisible.
`clinvar_plp_n` is 0 or 1 for every one of them.

**None of these may be reported as a causal gene–disease finding.** A biallelic hit
here is a research observation requiring functional validation, and the panel says so
in the `citation_source` column, verbatim:
`none: absent from DDG2P 2026-08-27, ClinGen 2026-08-27, HPO 2026-06-23`.

Also note *CDC20*: it has an OMIM link (`OMIM:620276`), but to
*oocyte/zygote/embryo maturation arrest 14* — a female-infertility phenotype, not a
paediatric somatic one. It is kept for mechanism; its curated disease does not fit
this proband.

## 6. Known gaps and things this panel deliberately does not claim

- **MVA5 and MVA6 are not in the panel.** The OMIM phenotype-series numbering implies
  they exist, but `phenotype.hpoa` v2026-06-23 contains only MVA1, MVA2, MVA3, MVA4
  and MVA7. No gene is asserted for numbers that could not be read from a file.
- **`clinvar_plp_n` is a record count, not a gene-level evidence score.** It counts
  submitted records, so a gene with many submissions for one recurrent allele scores
  high, and `GENEINFO` names every overlapping gene (see *ATRIP*). Use it to sanity-
  check that a gene has *any* known pathogenic variation, not to rank genes.
- **`evidence_tier` takes the strongest available classification.** Where DDG2P and
  ClinGen disagree (§5.1) that means the panel is optimistic by one tier. The
  disagreement is in `ddg2p_confidence` and `clingen_classification`, both of which
  are carried per row so a ranker can implement a stricter rule (e.g. take the
  *weakest*) without re-curating.
- **One row per gene.** Genes with several curated diseases (e.g. *SMC1A*: CdLS2 and
  developmental epileptic encephalopathy; *POLE*; *BRCA2*) carry the row most relevant
  to this phenotype. The others are still in DDG2P/ClinGen and are not lost, just not
  in this file.
- **No copy-number, repeat-expansion or mosaicism content.** The panel is a gene list
  for SNV/indel reasoning. Repo tech debt TD-13 (repeat expansions invisible) and
  TD-03 apply unchanged.
- **`gene_symbol_aliases` comes from DDG2P's `previous gene symbols` column**, so it
  is empty for the 29 genes absent from DDG2P — including *CENATAC*, *MAD1L1* and
  *BUB3*. Join on `hgnc_id` or `ensembl_gene_id` in preference to the symbol.

## 7. Relationship to `knowledge/public/gene_panel.tsv`

`knowledge/public/gene_panel.tsv` is a **synthetic, fictional** five-gene table
(`SYNTHKIN1`…`SYNTHSOL5`). It is referenced by `src/mva/config.py:299` as
`paths.gene_panel` but — as of this snapshot — **no adapter or pipeline stage reads
it**; grep finds no consumer for `mechanism_class`, `full_name` or
`disease_association` anywhere in `src/mva/`. Replacing it is therefore not enough on
its own; a consumer has to exist.

This file is a *different artifact* and does not overwrite it. Projection onto that
schema, if a loader is written, is column-for-column:

| `gene_panel.tsv` column | from `mva_panel.tsv` |
|---|---|
| `gene_symbol` | `gene_symbol` |
| `gene_id` | `ensembl_gene_id` |
| `full_name` | `disease_name` (no gene full-name column is carried; MANE column 5 has one if needed) |
| `inheritance` | `inheritance` |
| `mechanism_class` | `mechanism_class` |
| `disease_association` | `disease_name` + `disease_omim` |
| `source` | `citation_source` |
| `version` | `2026-08-28-a` |

The extra columns that have no home in that schema — `evidence_tier`, `panel_tier`,
`ddg2p_confidence`, `clingen_classification`, `citation`, `clinvar_plp_n`,
`gene_symbol_aliases` — are the ones ranking actually needs, so the better move is to
teach the loader this schema rather than to down-project.
