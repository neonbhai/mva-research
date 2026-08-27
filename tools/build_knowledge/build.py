"""Build ``knowledge/real/*.tsv`` from downloaded, real gene-disease resources.

Deterministic and offline: reads local files only, under ``--resources-root`` (default
``$MVA_RESOURCES``), never touches the network, and writes no timestamp except the
fixed :data:`tools.build_knowledge.tsv_io.RETRIEVED_DATE` recorded in every table's
header comment. Re-running against the same input bytes reproduces byte-identical
output (GP-30); ``tests/unit/test_knowledge_tables.py`` proves this on a small fixture.

Usage::

    uv run python -m tools.build_knowledge.build \\
        [--resources-root PATH] [--output-dir PATH]

Expected layout under the resources root (matching what was actually downloaded)::

    genepanels/clingen_gene_validity.csv
    genepanels/DDG2P.csv              (optional, gzip or plain -- see gene_disease.try_parse_ddg2p)
    gnomad/gnomad.v4.1.constraint_metrics.tsv
    hpo/genes_to_phenotype.txt
    hpo/phenotype.hpoa
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from tools.build_knowledge.gene_disease import (
    CLINGEN_MOI_GLOSSARY,
    GENE_DISEASE_COLUMNS,
    AssociationRow,
    ClinGenResult,
    DDG2PResult,
    parse_clingen,
    try_parse_ddg2p,
)
from tools.build_knowledge.gene_phenotype import (
    GENE_PHENOTYPE_COLUMNS,
    HpoGenePhenotypeResult,
    parse_gene_to_phenotype,
    read_hpo_release,
)
from tools.build_knowledge.gnomad_constraint import (
    GENE_PANEL_COLUMNS,
    GnomadConstraintResult,
    parse_gnomad_constraint,
)
from tools.build_knowledge.tsv_io import RETRIEVED_DATE, format_cell, format_float, write_tsv

#: Real definitions, quoted directly from the downloaded hp.obo (HPO's frequency
#: subontology, term IDs HP:0040280-HP:0040284). Documentation only: association_strength
#: always carries the source's own token verbatim, never this expansion.
_HPO_FREQUENCY_GLOSSARY: Final[tuple[tuple[str, str, str], ...]] = (
    ("HP:0040280", "Obligate", "Always present, i.e. in 100% of the cases."),
    ("HP:0040281", "Very frequent", "Present in 80% to 99% of the cases."),
    ("HP:0040282", "Frequent", "Present in 30% to 79% of the cases."),
    ("HP:0040283", "Occasional", "Present in 5% to 29% of the cases."),
    ("HP:0040284", "Very rare", "Present in 1% to 4% of the cases."),
)

_SAMPLE_SIZE: Final = 20


def _fmt_sample(values: tuple[str, ...], *, limit: int = _SAMPLE_SIZE) -> str:
    if not values:
        return "(none)"
    shown = ", ".join(values[:limit])
    if len(values) > limit:
        return f"{shown}, ... ({len(values) - limit} more)"
    return shown


@dataclass(frozen=True, slots=True)
class ResourcePaths:
    clingen: Path
    ddg2p: Path
    gnomad_constraint: Path
    hpo_genes_to_phenotype: Path
    hpo_phenotype_hpoa: Path

    @staticmethod
    def under(resources_root: Path) -> ResourcePaths:
        return ResourcePaths(
            clingen=resources_root / "genepanels" / "clingen_gene_validity.csv",
            ddg2p=resources_root / "genepanels" / "DDG2P.csv",
            gnomad_constraint=resources_root / "gnomad" / "gnomad.v4.1.constraint_metrics.tsv",
            hpo_genes_to_phenotype=resources_root / "hpo" / "genes_to_phenotype.txt",
            hpo_phenotype_hpoa=resources_root / "hpo" / "phenotype.hpoa",
        )


@dataclass(frozen=True, slots=True)
class BuildReport:
    gnomad: GnomadConstraintResult
    clingen: ClinGenResult
    ddg2p: DDG2PResult
    hpo: HpoGenePhenotypeResult
    gene_panel_path: Path
    gene_disease_path: Path
    gene_phenotype_path: Path


def _build_gene_panel(result: GnomadConstraintResult, *, output_path: Path) -> None:
    unmatched = result.genes_without_ensembl_row
    comments = [
        "Real gene-level constraint metrics from gnomAD.",
        "",
        "Source: gnomAD gene constraint metrics "
        "(gnomad.v4.1.constraint_metrics.tsv), downloaded from",
        "  https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/"
        "gnomad.v4.1.constraint_metrics.tsv",
        f"Version: {result.version} (parsed from the source filename; the TSV body "
        "carries no version column)",
        f"Retrieved: {RETRIEVED_DATE}",
        "Licence: public domain (CC0). gnomAD requests attribution/citation to the "
        "flagship publication wherever the data are used; no other restriction applies.",
        "  See https://gnomad.broadinstitute.org/policies",
        "Generator: tools/build_knowledge/build.py (uv run python -m tools.build_knowledge.build)",
        "",
        "One row per gene: gnomAD reports constraint per TRANSCRIPT, so for each gene "
        "the MANE Select transcript is preferred, then the canonical transcript, then "
        "the lexicographically smallest transcript ID (deterministic tie-break; mirrors "
        "mva.annotation.local_tables._consequence_sort_key). Only transcripts carrying "
        "an Ensembl gene ID (gene_id starting 'ENSG') are used, so ensembl_gene_id is "
        "always populated for an included gene; the gnomAD 'gene' value literally equal "
        "to 'NA' (its own missing-gene sentinel) is excluded entirely, not treated as a "
        "gene symbol.",
        "",
        "pli = lof.pLI, loeuf = lof.oe_ci.upper, mis_z = mis.z_score -- gnomAD's own "
        "high-confidence-pLoF-only ('lof.', not the merged 'lof_hc_lc.' set) and "
        "sample-size-calibrated ('z_score', not 'z_raw') columns: the same metrics "
        "gnomAD's browser shows on a gene's page. A gene with an insufficient expected-"
        "variant count for a metric ('NA' in the source) has an EMPTY cell here, NEVER "
        '0 -- an empty LOEUF means unmeasured, not "extremely constrained" (GP-14).',
        "constraint_flags carries gnomAD's own QC flags for the selected transcript "
        "(e.g. no_exp_lof, outlier_mis), comma-separated; empty means no flag was set.",
        "",
        "KNOWN COVERAGE GAP in this download (found by inspecting the file, not "
        "assumed): it contains autosomal genes (chr1-chr22) ONLY. No chrX or chrY gene "
        "appears anywhere in the source (verified directly: e.g. DMD, AR, ATRX, ATP7A, "
        "BTK are entirely absent), and the mitochondrial genome has no nuclear-LOFTEE-"
        "style constraint model at all. A gene absent from this table may be real and "
        'clinically important -- absence here means "not covered by this gnomAD '
        'release/download," never "unconstrained."',
        "",
        f"Accounting: {result.total_transcript_rows} transcript rows read, "
        f"{result.na_gene_rows} with gnomAD's own 'NA' gene sentinel (excluded), "
        f"{result.genes_seen} distinct real gene symbols seen, "
        f"{len(unmatched)} of those have no Ensembl-tagged transcript row at all and "
        f"are excluded from this table: {_fmt_sample(unmatched)}.",
        f"{result.genes_included} genes included below.",
    ]

    rows: list[dict[str, str]] = []
    for gene_row in result.rows:
        rows.append(
            {
                "gene_symbol": gene_row.gene_symbol,
                "ensembl_gene_id": gene_row.ensembl_gene_id,
                "transcript_id": gene_row.transcript_id,
                "canonical": "true" if gene_row.canonical else "false",
                "mane_select": "true" if gene_row.mane_select else "false",
                "pli": format_float(gene_row.pli),
                "loeuf": format_float(gene_row.loeuf),
                "mis_z": format_float(gene_row.mis_z),
                "constraint_flags": ",".join(gene_row.constraint_flags),
                "source": "gnomAD",
                "version": result.version,
            }
        )

    write_tsv(output_path, comment_lines=comments, columns=GENE_PANEL_COLUMNS, rows=rows)


def _build_gene_disease(
    clingen: ClinGenResult,
    ddg2p: DDG2PResult,
    gnomad_genes: frozenset[str],
    *,
    output_path: Path,
) -> None:
    unmatched_in_gnomad = tuple(sorted(g for g in clingen.unique_genes if g not in gnomad_genes))
    unmatched_xl = tuple(g for g in unmatched_in_gnomad if "XL" in _moi_codes_for_gene(clingen, g))
    moi_glossary = "; ".join(f"{code}={name}" for code, name in CLINGEN_MOI_GLOSSARY.items())
    classification_summary = "; ".join(
        f"{name}={count}" for name, count in sorted(clingen.classification_counts.items())
    )

    ddg2p_block: list[str]
    if ddg2p.available:
        ddg2p_missing_disease_id = sum(1 for row in ddg2p.rows if row.disease_id is None)
        ddg2p_unmatched = tuple(sorted(g for g in ddg2p.unique_genes if g not in gnomad_genes))
        ddg2p_unmatched_xl = tuple(
            g
            for g in ddg2p_unmatched
            if any(
                code.startswith(("monoallelic_X", "monoallelic_Y"))
                for code in _inheritance_codes_for_gene(ddg2p, g)
            )
        )
        ddg2p_confidence_summary = "; ".join(
            f"{name}={count}" for name, count in sorted((ddg2p.confidence_counts or {}).items())
        )
        ddg2p_block = [
            "Source 2: EBI Gene2Phenotype (G2P) developmental-disorder panel export "
            "(DDG2P.csv), downloaded from",
            "  https://www.ebi.ac.uk/gene2phenotype/api/panel/DD/download/",
            f"Version: {RETRIEVED_DATE} (accessed-on date; this G2P export does not embed "
            "a dataset release/version string of its own)",
            "Licence: EBI Gene2Phenotype terms -- freely available for use, with "
            "attribution requested. See https://www.ebi.ac.uk/gene2phenotype/about/project",
            f"{ddg2p.total_rows} associations across {len(ddg2p.unique_genes)} genes. "
            f"confidence values seen: {ddg2p_confidence_summary}.",
            "hgnc_id is DDG2P's own 'hgnc id' column, prefixed 'HGNC:'; "
            "mutation_consequence is DDG2P's 'variant consequence'; inheritance is "
            "DDG2P's own 'allelic requirement' vocabulary (biallelic_autosomal, "
            "monoallelic_autosomal, monoallelic_X_hemizygous, monoallelic_X_heterozygous, "
            "monoallelic_X, monoallelic_Y_hemizygous, mitochondrial) -- a DIFFERENT code "
            "space from ClinGen's MOI column above; the two are never merged into one "
            "vocabulary, only stacked as separate rows distinguished by `source`. "
            "classification_date is DDG2P's 'date of last review'; gcep is always empty "
            "for DDG2P rows (ClinGen-specific concept); panel is DDG2P's own "
            "semicolon-joined sub-panel membership (e.g. 'DD; Eye'), always empty for "
            "ClinGen rows.",
            f"disease_id prefers DDG2P's 'disease MONDO' column, falls back to "
            f"'disease mim' (prefixed 'OMIM:') when MONDO is blank, and is left EMPTY "
            f"when both are blank ({ddg2p_missing_disease_id} of {ddg2p.total_rows} rows) "
            "-- never a fabricated identifier.",
            f"Gene-symbol cross-check against knowledge/real/gene_panel.tsv (gnomAD): "
            f"{len(ddg2p_unmatched)} of {len(ddg2p.unique_genes)} DDG2P gene symbols "
            f"were not found there ({len(ddg2p_unmatched_xl)} carry an X-linked/Y-linked "
            "allelic requirement and are explained by gene_panel.tsv's chrX/Y coverage "
            f"gap); remaining sample: "
            f"{_fmt_sample(tuple(g for g in ddg2p_unmatched if g not in ddg2p_unmatched_xl))}.",
        ]
    else:
        ddg2p_block = [
            "Source 2: EBI Gene2Phenotype DDG2P developmental-disorder panel -- UNAVAILABLE.",
            f"  Reason: {ddg2p.reason}",
            "  This table therefore carries ClinGen associations only. "
            "tools/build_knowledge/gene_disease.py:try_parse_ddg2p will parse a "
            "corrected download automatically on the next run -- no code change needed.",
        ]

    comments = [
        "Real gene-disease clinical validity associations.",
        "",
        "Source 1: ClinGen Gene-Disease Validity curations "
        "(clingen_gene_validity.csv), downloaded from",
        "  https://search.clinicalgenome.org/kb/gene-validity/download",
        f"Version: {clingen.version} (ClinGen's own 'FILE CREATED:' date; ClinGen is a "
        "living curation, not a tagged release)",
        "Licence: CC0 1.0 Universal (public domain); ClinGen requests attribution and "
        "the date accessed.",
        "  See https://clinicalgenome.org/docs/terms-of-use/",
        *ddg2p_block,
        f"Retrieved: {RETRIEVED_DATE}",
        "Generator: tools/build_knowledge/build.py (uv run python -m tools.build_knowledge.build)",
        "",
        "Every classification from BOTH sources is kept VERBATIM in confidence, "
        "including ClinGen's Disputed, Refuted and 'No Known Disease Relationship' and "
        "DDG2P's (lower-case) disputed/refuted: curated-and-refuted is real evidence, "
        "and dropping it would make absence indistinguishable from refutation (GP-14, "
        f"generalised). ClinGen values seen in this download: {classification_summary}.",
        f"inheritance carries ClinGen's own MOI code unmodified for ClinGen rows "
        f"({moi_glossary}); mutation_consequence is always empty for ClinGen rows -- "
        "ClinGen's curation model does not capture this concept.",
        "",
        f"Gene-symbol cross-check against knowledge/real/gene_panel.tsv (gnomAD): "
        f"{len(unmatched_in_gnomad)} of {len(clingen.unique_genes)} ClinGen gene symbols "
        "were not found there.",
        f"  {len(unmatched_xl)} of those carry an X-linked (XL) MOI record and are "
        "explained by gene_panel.tsv's documented chrX/Y coverage gap, not a symbol "
        "mismatch.",
        f"  The remaining {len(unmatched_in_gnomad) - len(unmatched_xl)} are mostly "
        "mitochondrial genes (gnomAD's nuclear constraint model does not cover mtDNA) "
        "and a handful of recent HGNC renames / non-protein-coding loci: "
        f"{_fmt_sample(tuple(g for g in unmatched_in_gnomad if g not in unmatched_xl))}.",
    ]

    all_rows: tuple[AssociationRow, ...] = clingen.rows + ddg2p.rows
    formatted: list[dict[str, str]] = [
        {
            "gene_symbol": row.gene_symbol,
            "hgnc_id": format_cell(row.hgnc_id),
            "disease_id": format_cell(row.disease_id),
            "disease_name": row.disease_name,
            "inheritance": format_cell(row.inheritance),
            "mutation_consequence": format_cell(row.mutation_consequence),
            "confidence": row.confidence,
            "classification_date": format_cell(row.classification_date),
            "gcep": format_cell(row.gcep),
            "panel": format_cell(row.panel),
            "source": row.source,
            "version": row.version,
        }
        for row in sorted(all_rows, key=_association_sort_key)
    ]
    write_tsv(output_path, comment_lines=comments, columns=GENE_DISEASE_COLUMNS, rows=formatted)


def _association_sort_key(row: AssociationRow) -> tuple[str, str, str, str, str]:
    return (
        row.gene_symbol,
        row.disease_id or "",
        row.source,
        row.classification_date or "",
        row.gcep or "",
    )


def _moi_codes_for_gene(clingen: ClinGenResult, gene: str) -> frozenset[str]:
    return frozenset(row.inheritance or "" for row in clingen.rows if row.gene_symbol == gene)


def _inheritance_codes_for_gene(ddg2p: DDG2PResult, gene: str) -> frozenset[str]:
    return frozenset(row.inheritance or "" for row in ddg2p.rows if row.gene_symbol == gene)


def _build_gene_phenotype(result: HpoGenePhenotypeResult, *, output_path: Path) -> None:
    glossary_lines = [
        f'  {term_id} {name:<14s} "{definition}"'
        for term_id, name, definition in _HPO_FREQUENCY_GLOSSARY
    ]
    comments = [
        "Real gene-to-phenotype (HPO term) annotations, restricted to genes curated in",
        "knowledge/real/gene_disease.tsv.",
        "",
        "Source: HPO gene-to-phenotype annotations (genes_to_phenotype.txt), downloaded",
        "  from https://github.com/obophenotype/human-phenotype-ontology/releases/"
        "latest/download/genes_to_phenotype.txt",
        f"Version: {result.version} (HPO's own release date, read from "
        "phenotype.hpoa's '#version:' line)",
        "Licence: CC BY 4.0 (Human Phenotype Ontology Consortium). The HPO content and "
        "its logical relationships may not be altered; cite the HPO project.",
        f"Retrieved: {RETRIEVED_DATE}",
        "Generator: tools/build_knowledge/build.py (uv run python -m tools.build_knowledge.build)",
        "",
        f"Restriction and size: the full source file has {result.total_rows_read} rows "
        f"across many more genes than are curated anywhere in this pipeline. "
        f"Restricted here to the {result.restrict_gene_count} gene symbols in "
        f"knowledge/real/gene_disease.tsv: {result.genes_matched} of them have HPO "
        f"annotations ({result.rows_kept} matching source rows, collapsed to "
        f"{len(result.rows)} unique (gene, HPO term) pairs below). "
        f"{len(result.genes_unmatched)} genes from gene_disease.tsv have ZERO rows "
        "here -- this may be a real absence of curated phenotype annotation, or "
        "gene-symbol drift between sources; this generator does not attempt to tell the "
        f"two apart. Unmatched sample: {_fmt_sample(result.genes_unmatched)}.",
        "",
        "One row per (gene_symbol, hpo_id): the source has one row per (gene, HPO term, "
        "DISEASE), since the same phenotype can be annotated via several diseases with "
        "different evidence. Collapsed here by preferring a row that carries real "
        "frequency data over one recorded as '-', tie-broken deterministically on "
        "disease_id -- never averaged or otherwise invented.",
        "",
        "association_strength is HPO's OWN phenotype-FREQUENCY vocabulary (how often "
        "the phenotype occurs among cases of the linked disease), carried through "
        "VERBATIM -- NOT the definitive/strong/moderate/supporting curation-confidence "
        "scale the SYNTHETIC demo table (knowledge/public/gene_phenotype.tsv) invented "
        "for its fictional genes. These are different scientific concepts and must not "
        "be conflated. Real HPO frequency terms (quoted from hp.obo):",
        *glossary_lines,
        "A row may instead carry a raw fraction ('n/m' cases) or a percentage recorded "
        "directly in the source. An EMPTY cell means the source recorded '-' (frequency "
        'not stated) -- unmeasured, not "never occurs" (GP-14). HP:0040285 (Excluded, '
        "0% -- i.e. NOT part of the disease) does not appear in this file: HPO's own "
        "genes_to_phenotype.txt build already excludes negated annotations, so no "
        "exclusion filtering was implemented here.",
        "",
        "Downstream code that expects the closed definitive/strong/moderate/supporting "
        "vocabulary (src/mva/phenotype/hpo.py's STRENGTH_WEIGHTS) will need adapting "
        "before it can read this file directly; that adaptation is out of scope for this "
        "generator (src/mva/phenotype/ is owned elsewhere).",
    ]

    rows = [
        {
            "gene_symbol": row.gene_symbol,
            "hpo_id": row.hpo_id,
            "label": row.label,
            "association_strength": format_cell(row.association_strength),
            "source": "HPO",
            "version": result.version,
        }
        for row in result.rows
    ]
    write_tsv(output_path, comment_lines=comments, columns=GENE_PHENOTYPE_COLUMNS, rows=rows)


def build_all(resources_root: Path, output_dir: Path) -> BuildReport:
    """Parse every source once and write the three real knowledge tables.

    Pure aside from the three file writes at the end: parsing happens first and in
    full, so a malformed source fails before anything is written (no half-written
    table on error).
    """
    paths = ResourcePaths.under(resources_root)

    gnomad = parse_gnomad_constraint(paths.gnomad_constraint)
    clingen = parse_clingen(paths.clingen)
    ddg2p = try_parse_ddg2p(paths.ddg2p, version=RETRIEVED_DATE)
    hpo_version = read_hpo_release(paths.hpo_phenotype_hpoa)
    # Restricted to genes curated by EITHER clinical-validity source -- not an
    # arbitrary truncation, and not ClinGen alone: when DDG2P is available it widens
    # coverage exactly as this generator's brief asks for ("genes present in DDG2P or
    # ClinGen"). When DDG2P is unavailable, ddg2p.unique_genes is simply empty and
    # this reduces to ClinGen alone.
    disease_genes = frozenset(clingen.unique_genes) | frozenset(ddg2p.unique_genes)
    hpo = parse_gene_to_phenotype(
        paths.hpo_genes_to_phenotype,
        restrict_to_genes=disease_genes,
        version=hpo_version,
    )

    gene_panel_path = output_dir / "gene_panel.tsv"
    gene_disease_path = output_dir / "gene_disease.tsv"
    gene_phenotype_path = output_dir / "gene_phenotype.tsv"

    gnomad_genes = frozenset(row.gene_symbol for row in gnomad.rows)
    _build_gene_panel(gnomad, output_path=gene_panel_path)
    _build_gene_disease(clingen, ddg2p, gnomad_genes, output_path=gene_disease_path)
    _build_gene_phenotype(hpo, output_path=gene_phenotype_path)

    return BuildReport(
        gnomad=gnomad,
        clingen=clingen,
        ddg2p=ddg2p,
        hpo=hpo,
        gene_panel_path=gene_panel_path,
        gene_disease_path=gene_disease_path,
        gene_phenotype_path=gene_phenotype_path,
    )


def _print_report(report: BuildReport) -> None:
    def _size(path: Path) -> str:
        return f"{path.stat().st_size:,} bytes"

    print("== knowledge/real build report ==")  # noqa: T201 -- this is the CLI's report
    print(  # noqa: T201
        f"gene_panel.tsv:     {report.gnomad.genes_included:>6} rows, "
        f"{_size(report.gene_panel_path)}  ({report.gene_panel_path})"
    )
    print(  # noqa: T201
        f"gene_disease.tsv:   {len(report.clingen.rows) + len(report.ddg2p.rows):>6} rows, "
        f"{_size(report.gene_disease_path)}  ({report.gene_disease_path})"
    )
    print(  # noqa: T201
        f"gene_phenotype.tsv: {len(report.hpo.rows):>6} rows, "
        f"{_size(report.gene_phenotype_path)}  ({report.gene_phenotype_path})"
    )
    print()  # noqa: T201
    print(  # noqa: T201
        f"gnomAD: {report.gnomad.total_transcript_rows} transcript rows -> "
        f"{report.gnomad.genes_included} genes; "
        f"{len(report.gnomad.genes_without_ensembl_row)} genes excluded (no Ensembl-"
        "tagged transcript)."
    )
    print(  # noqa: T201
        f"ClinGen: {report.clingen.total_rows} associations across "
        f"{len(report.clingen.unique_genes)} genes."
    )
    if report.ddg2p.available:
        print(  # noqa: T201
            f"DDG2P: {report.ddg2p.total_rows} associations across "
            f"{len(report.ddg2p.unique_genes)} genes."
        )
    else:
        print(f"DDG2P: UNAVAILABLE -- {report.ddg2p.reason}")  # noqa: T201
    print(  # noqa: T201
        f"HPO gene_phenotype: restricted to {report.hpo.restrict_gene_count} genes "
        "curated by ClinGen and/or DDG2P; "
        f"{report.hpo.genes_matched} matched, "
        f"{len(report.hpo.genes_unmatched)} unmatched."
    )


def _resolve_resources_root(cli_value: Path | None) -> Path:
    if cli_value is not None:
        return cli_value
    env_value = os.environ.get("MVA_RESOURCES")
    if env_value:
        return Path(env_value)
    msg = (
        "No resources root given: pass --resources-root or set the MVA_RESOURCES "
        "environment variable to the directory containing genepanels/, gnomad/ and "
        "hpo/ subdirectories."
    )
    raise SystemExit(msg)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--resources-root",
        type=Path,
        default=None,
        help="Directory containing genepanels/, gnomad/ and hpo/ (default: $MVA_RESOURCES).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory to write gene_panel.tsv / gene_disease.tsv / gene_phenotype.tsv into "
        "(default: <repo root>/knowledge/real).",
    )
    args = parser.parse_args(argv)

    resources_root = _resolve_resources_root(args.resources_root)
    if args.output_dir is not None:
        output_dir = args.output_dir
    else:
        repo_root = Path(__file__).resolve().parents[2]
        output_dir = repo_root / "knowledge" / "real"

    report = build_all(resources_root, output_dir)
    _print_report(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
