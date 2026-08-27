"""Parse HPO's own flat gene-to-phenotype table into one row per (gene, HPO term).

``genes_to_phenotype.txt`` is already the flat gene -> HPO-term view the HPO project
builds from ``phenotype.hpoa``; this module does no ontology-graph work of its own (no
term propagation, no parent/child reasoning, no disease-annotation qualifier logic) --
it only reads what HPO already flattened, per this task's explicit boundary with
``src/mva/phenotype/`` (owned elsewhere).

The source has one row per (gene, HPO term, **disease**): the same phenotype can be
annotated on a gene through several diseases, each with its own frequency evidence.
Collapsing to one row per (gene, HPO term) prefers whichever disease-row carries real
frequency data over one recorded as ``"-"`` (not stated), tie-broken deterministically
on ``disease_id`` -- never averaged or otherwise invented.

``association_strength`` carries HPO's OWN phenotype-**frequency** vocabulary verbatim
(one of five real HPO frequency-subontology term IDs, a fraction of annotated cases, or
a percentage), NOT the ``definitive``/``strong``/``moderate``/``supporting``
curation-confidence scale the SYNTHETIC demo table
(``knowledge/public/gene_phenotype.tsv``) invented for the hackathon's fictional genes.
Those are different scientific concepts: population frequency of a symptom within a
disease is not the same claim as curated confidence that a gene causes a phenotype at
all, and collapsing one onto the other would fabricate precision this project's rules
forbid (GP-14). See ``build.py`` for the real term definitions, quoted from ``hp.obo``,
recorded in the rendered table's header comment.
"""

from __future__ import annotations

import re
from collections.abc import Collection
from dataclasses import dataclass
from pathlib import Path
from typing import Final

GENE_PHENOTYPE_COLUMNS: Final[tuple[str, ...]] = (
    "gene_symbol",
    "hpo_id",
    "label",
    "association_strength",
    "source",
    "version",
)

#: HPO's own "not stated" sentinel for the frequency column in this file.
_NOT_STATED: Final = "-"

_REQUIRED_SOURCE_COLUMNS: Final[tuple[str, ...]] = (
    "gene_symbol",
    "hpo_id",
    "hpo_name",
    "frequency",
    "disease_id",
)

_HPOA_VERSION_RE: Final[re.Pattern[str]] = re.compile(r"^#version:\s*(?P<version>\S+)")


@dataclass(frozen=True, slots=True)
class GenePhenotypeRow:
    gene_symbol: str
    hpo_id: str
    label: str
    association_strength: str | None


@dataclass(frozen=True, slots=True)
class HpoGenePhenotypeResult:
    rows: tuple[GenePhenotypeRow, ...]
    version: str
    total_rows_read: int
    rows_kept: int
    restrict_gene_count: int
    genes_matched: int
    genes_unmatched: tuple[str, ...]


def read_hpo_release(phenotype_hpoa_path: Path) -> str:
    """The real HPO release string, read from ``phenotype.hpoa``'s own ``#version:`` line."""
    if not phenotype_hpoa_path.is_file():
        msg = f"HPO annotation file not found: {phenotype_hpoa_path}"
        raise FileNotFoundError(msg)
    with phenotype_hpoa_path.open(encoding="utf-8") as handle:
        for _ in range(10):
            line = handle.readline()
            if not line:
                break
            match = _HPOA_VERSION_RE.match(line.strip())
            if match:
                return match.group("version")
    msg = (
        f"{phenotype_hpoa_path.name}: no '#version:' header line found in the first 10 "
        "lines. Refusing to invent an HPO release string (GP-18)."
    )
    raise ValueError(msg)


def parse_gene_to_phenotype(
    path: Path, *, restrict_to_genes: Collection[str], version: str
) -> HpoGenePhenotypeResult:
    """Parse ``genes_to_phenotype.txt``, kept to genes in ``restrict_to_genes``.

    The restriction exists for size, not correctness: the full file annotates
    thousands of genes with no curated clinical-validity record anywhere in this
    pipeline (see ``build.py`` for the measured before/after row counts). A gene in
    ``restrict_to_genes`` with zero matching rows is reported in
    :attr:`HpoGenePhenotypeResult.genes_unmatched`, never silently dropped.
    """
    if not path.is_file():
        msg = f"HPO gene-to-phenotype file not found: {path}"
        raise FileNotFoundError(msg)
    restrict_set = frozenset(restrict_to_genes)

    best_by_pair: dict[tuple[str, str], tuple[tuple[int, str], GenePhenotypeRow]] = {}
    total_rows = 0
    rows_kept = 0
    genes_matched: set[str] = set()

    with path.open(encoding="utf-8") as handle:
        header_line = handle.readline()
        header = header_line.rstrip("\n").split("\t")
        idx = {name: position for position, name in enumerate(header)}
        missing_columns = [c for c in _REQUIRED_SOURCE_COLUMNS if c not in idx]
        if missing_columns:
            msg = f"{path.name} is missing column(s) {missing_columns}."
            raise ValueError(msg)

        for line in handle:
            total_rows += 1
            fields = line.rstrip("\n").split("\t")
            gene = fields[idx["gene_symbol"]]
            if gene not in restrict_set:
                continue
            rows_kept += 1
            genes_matched.add(gene)

            hpo_id = fields[idx["hpo_id"]]
            label = fields[idx["hpo_name"]]
            frequency = fields[idx["frequency"]]
            disease_id = fields[idx["disease_id"]]
            strength = None if frequency == _NOT_STATED else frequency

            candidate = GenePhenotypeRow(
                gene_symbol=gene, hpo_id=hpo_id, label=label, association_strength=strength
            )
            # Prefer a row with real frequency data over one without; break remaining
            # ties on disease_id for a total, deterministic order (GP-30).
            sort_key = (0 if strength is not None else 1, disease_id)
            pair = (gene, hpo_id)
            existing = best_by_pair.get(pair)
            if existing is None or sort_key < existing[0]:
                best_by_pair[pair] = (sort_key, candidate)

    rows = tuple(
        sorted(
            (entry[1] for entry in best_by_pair.values()),
            key=lambda r: (r.gene_symbol, r.hpo_id),
        )
    )
    genes_unmatched = tuple(sorted(restrict_set - genes_matched))
    return HpoGenePhenotypeResult(
        rows=rows,
        version=version,
        total_rows_read=total_rows,
        rows_kept=rows_kept,
        restrict_gene_count=len(restrict_set),
        genes_matched=len(genes_matched),
        genes_unmatched=genes_unmatched,
    )
