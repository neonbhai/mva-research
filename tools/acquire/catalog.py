"""The resource registry: a declarative table of public reference resources.

This module names *what* this tool knows how to fetch and *where it comes from* --
nothing here reads the filesystem or the network. ``tools.acquire.inspect`` is the
module that cross-references this catalog against what is actually on disk to
produce the record that gets written to ``knowledge/manifests/resources.yaml``.

Every entry declared here is deliberately hand-written, not glob-discovered the way
``mva.annotation.local_tables.compute_manifest`` discovers ``knowledge/public/*.tsv``:
a license string, a real release version and an upstream URL are facts about the
world that cannot be inferred from a filename, and inventing them would be worse
than an admittedly manual table.
"""

from __future__ import annotations

from typing import Final

from tools.acquire.models import ResourceEntry

# --------------------------------------------------------------------------- licenses
#
# Kept as named constants because getting these wrong in the committed manifest is
# the whole failure mode this file exists to avoid.

_CLINVAR_LICENSE: Final = (
    "Public domain (NCBI / U.S. Government work); no usage restriction. Citation "
    "requested: Landrum MJ et al., 'ClinVar: improving access to variant "
    "interpretations and supporting evidence', Nucleic Acids Res. 2018."
)

_GNOMAD_LICENSE: Final = (
    "CC0 1.0 Universal (public domain dedication). No usage restriction, but the "
    "gnomAD consortium requests citation of the relevant flagship paper (Chen S et "
    "al., 'A genomic mutational constraint map using variation in 76,156 human "
    "genomes', Nature 2024, for v4; Karczewski KJ et al., Nature 2020, for the "
    "original gnomAD paper) -- see https://gnomad.broadinstitute.org/about."
)

_HPO_LICENSE: Final = (
    "CC BY 4.0 (Human Phenotype Ontology Consortium). Attribution required; see "
    "https://hpo.jax.org/app/license."
)

_G2P_LICENSE: Final = (
    "EMBL-EBI terms of use (https://www.ebi.ac.uk/about/terms-of-use): EMBL-EBI "
    "imposes no additional restriction beyond the data owner's own terms, and "
    "expects attribution in line with good scientific practice. Cite: Thormann A "
    "et al., 'Flexible and scalable diagnostic filtering of genomic variants using "
    "G2P with Ensembl VEP', Nat Commun 2019, and the DD (Developmental Disorders) "
    "panel curators at https://www.ebi.ac.uk/gene2phenotype."
)

_CLINGEN_LICENSE: Final = (
    "CC0 1.0 Universal (public domain dedication) -- ClinGen Gene-Disease Validity "
    "curations; see https://search.clinicalgenome.org/kb/gene-validity."
)

# --------------------------------------------------------------------------- fixed entries

_CLINVAR: Final[tuple[ResourceEntry, ...]] = (
    ResourceEntry(
        name="clinvar_vcf",
        url="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz",
        version="2026-08-22",
        path="clinvar/clinvar.vcf.gz",
        license=_CLINVAR_LICENSE,
        description="ClinVar clinical variant assertions, GRCh38 weekly VCF release.",
    ),
    ResourceEntry(
        name="clinvar_vcf_tbi",
        url="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi",
        version="2026-08-22",
        path="clinvar/clinvar.vcf.gz.tbi",
        license=_CLINVAR_LICENSE,
        description="Tabix index for the ClinVar GRCh38 weekly VCF.",
    ),
)

_HPO_RELEASE_VERSION: Final = "2026-06-23"
_HPO_BASE_URL: Final = (
    "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download"
)

_HPO: Final[tuple[ResourceEntry, ...]] = (
    ResourceEntry(
        name="hpo_ontology",
        url=f"{_HPO_BASE_URL}/hp.obo",
        version=_HPO_RELEASE_VERSION,
        path="hpo/hp.obo",
        license=_HPO_LICENSE,
        description="Human Phenotype Ontology, OBO format.",
    ),
    ResourceEntry(
        name="hpo_disease_annotations",
        url=f"{_HPO_BASE_URL}/phenotype.hpoa",
        version=_HPO_RELEASE_VERSION,
        path="hpo/phenotype.hpoa",
        license=_HPO_LICENSE,
        description="HPO disease-to-phenotype annotation table (phenotype.hpoa).",
    ),
    ResourceEntry(
        name="hpo_genes_to_phenotype",
        url=f"{_HPO_BASE_URL}/genes_to_phenotype.txt",
        version=_HPO_RELEASE_VERSION,
        path="hpo/genes_to_phenotype.txt",
        license=_HPO_LICENSE,
        description="HPO gene-to-phenotype association table.",
    ),
    ResourceEntry(
        name="hpo_phenotype_to_genes",
        url=f"{_HPO_BASE_URL}/phenotype_to_genes.txt",
        version=_HPO_RELEASE_VERSION,
        path="hpo/phenotype_to_genes.txt",
        license=_HPO_LICENSE,
        description="HPO phenotype-to-gene association table.",
    ),
)

_GNOMAD_CONSTRAINT: Final[tuple[ResourceEntry, ...]] = (
    ResourceEntry(
        name="gnomad_constraint_metrics",
        url=(
            "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/"
            "constraint/gnomad.v4.1.constraint_metrics.tsv"
        ),
        version="4.1",
        path="gnomad/gnomad.v4.1.constraint_metrics.tsv",
        license=_GNOMAD_LICENSE,
        description="gnomAD v4.1 per-transcript constraint metrics (pLI, LOEUF, missense o/e).",
    ),
)

_GENEPANELS: Final[tuple[ResourceEntry, ...]] = (
    ResourceEntry(
        name="ddg2p_dd_panel",
        # The historical bulk-download URL (.../downloads/DDG2P.csv.gz) now serves an
        # HTML app shell instead of data -- Gene2Phenotype moved bulk export behind a
        # per-panel API. This is the corrected endpoint, confirmed to return
        # `content-type: text/csv` for the DD (Developmental Disorders) panel, which is
        # the successor to the DDG2P dataset this tool was asked to register.
        url="https://www.ebi.ac.uk/gene2phenotype/api/panel/DD/download/",
        version="DD-panel-2026-08-27",
        path="genepanels/DDG2P.csv",
        license=_G2P_LICENSE,
        description=(
            "Gene2Phenotype Developmental Disorders (DD) gene panel, CSV export -- the "
            "successor to the DDG2P dataset (gene-disease associations with allelic "
            "requirement, mutation consequence and confidence)."
        ),
    ),
    ResourceEntry(
        name="clingen_gene_validity",
        url="https://search.clinicalgenome.org/kb/gene-validity/download",
        version="2026-08-27",
        path="genepanels/clingen_gene_validity.csv",
        license=_CLINGEN_LICENSE,
        description="ClinGen Gene-Disease Validity curations, all genes and classifications.",
    ),
)

#: Chromosomes gnomAD v4.1 exomes ships sites-only VCFs for. Declared as a fixed plan
#: because that is what was actually asked for (see gnomad.sh's CHROMS list); which of
#: these are FETCHED vs NOT_FETCHED is a question for tools.acquire.inspect, not this file.
GNOMAD_EXOME_CHROMOSOMES: Final[tuple[str, ...]] = (*(str(n) for n in range(1, 23)), "X", "Y")

_GNOMAD_EXOME_BASE_URL: Final = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes"
)


def _gnomad_exome_entries(chromosomes: tuple[str, ...]) -> tuple[ResourceEntry, ...]:
    """One VCF + one tabix-index entry per chromosome, built from a single URL/path pattern.

    This is the "register them by pattern" the acquisition tool was asked for: one
    template, expanded over the declared chromosome list, rather than 48 hand-typed
    near-duplicate entries.
    """
    entries: list[ResourceEntry] = []
    for chrom in chromosomes:
        filename = f"gnomad.exomes.v4.1.sites.chr{chrom}.vcf.bgz"
        entries.append(
            ResourceEntry(
                name=f"gnomad_exomes_chr{chrom}",
                url=f"{_GNOMAD_EXOME_BASE_URL}/{filename}",
                version="4.1",
                path=f"gnomad/v4.1_exomes/{filename}",
                license=_GNOMAD_LICENSE,
                description=f"gnomAD v4.1 exome sites-only VCF (no genotypes), chromosome {chrom}.",
            )
        )
        entries.append(
            ResourceEntry(
                name=f"gnomad_exomes_chr{chrom}_tbi",
                url=f"{_GNOMAD_EXOME_BASE_URL}/{filename}.tbi",
                version="4.1",
                path=f"gnomad/v4.1_exomes/{filename}.tbi",
                license=_GNOMAD_LICENSE,
                description=(
                    f"Tabix index for the gnomAD v4.1 exome sites-only VCF, chromosome {chrom}."
                ),
            )
        )
    return tuple(entries)


#: The complete, declarative resource registry. Nullable-hash `ResourceEntry` objects
#: in their as-declared state (status=NOT_FETCHED, sha256=None) -- see
#: `tools.acquire.inspect.inspect_all` for how these get filled in from real disk state.
KNOWN_RESOURCES: Final[tuple[ResourceEntry, ...]] = (
    _CLINVAR
    + _HPO
    + _GNOMAD_CONSTRAINT
    + _GENEPANELS
    + _gnomad_exome_entries(GNOMAD_EXOME_CHROMOSOMES)
)

_names = [entry.name for entry in KNOWN_RESOURCES]
assert len(_names) == len(set(_names)), (
    "duplicate resource name in KNOWN_RESOURCES; every manifest entry must be unique"
)
