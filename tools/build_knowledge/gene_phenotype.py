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

Two qualifiers, two columns (ADR 0021)
--------------------------------------

This generator used to write HPO's frequency vocabulary into a column named
``association_strength``, which ``src/mva/phenotype/hpo.py`` validates as curated
gene-disease clinical validity. Both modules were internally consistent and the
table was unloadable. The two quantities now have their own columns, their own
vocabularies and their own parsers:

``hpo_frequency``
    HPO's OWN phenotype-frequency vocabulary, carried through VERBATIM: one of the
    five HP:00402xx frequency-subontology terms, an ``n/m`` fraction of annotated
    cases, or a percentage. Empty where the source recorded ``-``. This says how
    often a feature occurs *among cases of the linked disease*.

``association_strength``
    Curated gene-disease **clinical validity** -- how confident an expert panel is
    that variation in the gene causes disease at all -- taken from ClinGen
    Gene-Disease Validity and EBI Gene2Phenotype (DDG2P) via
    :func:`curated_strength_by_gene`. HPO's ``genes_to_phenotype.txt`` contains no
    such judgement, so it cannot be the source of this column. Empty where no
    source classifies the gene: absent, never a default (GP-14).

Population frequency of a symptom within a disease and curated confidence that a
gene causes disease at all are different scientific claims, and collapsing one onto
the other fabricates precision this project's rules forbid. See ``build.py`` for the
real term definitions, quoted from ``hp.obo``, recorded in the rendered table's
header comment.
"""

from __future__ import annotations

import re
from collections.abc import Collection, Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Final

GENE_PHENOTYPE_COLUMNS: Final[tuple[str, ...]] = (
    "gene_symbol",
    "hpo_id",
    "label",
    "association_strength",
    "association_strength_source",
    "hpo_frequency",
    "source",
    "version",
)

#: HPO's own "not stated" sentinel for the frequency column in this file.
_NOT_STATED: Final = "-"

#: HPO's "Excluded" (0% of cases) frequency term: a NEGATED annotation saying the
#: feature is *not* part of the disease. A gene-phenotype association table asserts
#: the opposite claim, so such a row is excluded here and counted in the rendered
#: header rather than emitted with an empty or a zero frequency. HPO's own
#: genes_to_phenotype.txt build already drops negated annotations, so the count is
#: expected to be 0 -- it is measured, not assumed.
EXCLUDED_FREQUENCY_TERM: Final = "HP:0040285"

#: ClinGen Gene-Disease Validity and DDG2P confidence classifications, strongest
#: first. Case-folded: ClinGen writes ``Definitive`` and DDG2P ``definitive`` for one
#: concept, and folding case is not a semantic remap. The order is the SOURCES' own
#: published ladder, not a scale invented here; :func:`curated_strength_by_gene`
#: reads it as a total order and nothing else.
#:
#: An unlisted classification raises rather than sorting to the bottom: a new ClinGen
#: tier is a curation-model change that must be looked at, not silently ranked.
CURATED_VALIDITY_RANK: Final[tuple[str, ...]] = (
    "definitive",
    "strong",
    "moderate",
    "limited",
    "disputed",
    "refuted",
    "no known disease relationship",
)

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
    """Curated gene-disease clinical validity; ``None`` when no source classifies
    this gene. Never derived from the HPO file -- see the module docstring."""
    association_strength_source: str | None
    """Which curation panel supplied ``association_strength``. ``None`` exactly when
    the strength is absent. Separate from ``source`` (which names where the
    gene->term ANNOTATION came from, always HPO) because the two really are
    different provenances, and one ``source`` column covering both would recreate
    the conflation ADR 0021 removed."""
    hpo_frequency: str | None
    """HPO's own frequency token, verbatim; ``None`` where the source said ``-``."""


@dataclass(frozen=True, slots=True)
class HpoGenePhenotypeResult:
    rows: tuple[GenePhenotypeRow, ...]
    version: str
    total_rows_read: int
    rows_kept: int
    restrict_gene_count: int
    genes_matched: int
    genes_unmatched: tuple[str, ...]
    genes_without_curated_strength: tuple[str, ...]
    """Genes emitted below with an EMPTY association_strength, because no curation
    source classifies them. Reported so the gap is visible and countable, never
    filled with a default."""
    rows_without_curated_strength: int
    excluded_frequency_rows: int
    """Source rows dropped for carrying HP:0040285 (Excluded, 0% of cases)."""


def curated_strength_by_gene(
    records: Iterable[tuple[str, str, str]],
) -> dict[str, tuple[str, str]]:
    """Strongest curated clinical-validity classification per gene, case-folded.

    ``records`` is ``(gene_symbol, classification, source)`` triples from any curation
    source; ``build.py`` passes ClinGen's and DDG2P's rows together, so the result is
    the strongest statement made about the gene by *either* panel, mapped to
    ``(classification, source)`` so the borrowed value keeps its provenance (GP-31).
    When both panels give the strongest tier, the source name that sorts first wins --
    a total, deterministic tie-break (GP-30), not a judgement about which panel is
    better.

    **Granularity, stated plainly.** ClinGen and DDG2P curate the (gene, DISEASE)
    pair; this reduces to the gene. A gene curated ``definitive`` for one disease and
    ``limited`` for another is recorded here as ``definitive``, which over-states the
    validity of the second disease's phenotype annotations. The alternative -- joining
    HPO's annotation to the curation for its *own* disease -- was measured and reaches
    only 75,379 of 275,046 source rows (27%): HPO annotates against OMIM and ORPHA,
    ClinGen curates against MONDO, and no crosswalk between them exists on disk. A
    column that is right for a quarter of rows and empty for the rest is not more
    honest than a gene-level one whose granularity is documented; see ADR 0021.

    Raises ``ValueError`` on a classification outside :data:`CURATED_VALIDITY_RANK`.
    """
    rank = {name: position for position, name in enumerate(CURATED_VALIDITY_RANK)}
    best: dict[str, tuple[str, str]] = {}
    for gene_symbol, classification, source in records:
        gene = gene_symbol.strip()
        token = " ".join(classification.strip().lower().split())
        panel = source.strip()
        if not gene or not token:
            continue
        if token not in rank:
            allowed = ", ".join(CURATED_VALIDITY_RANK)
            msg = (
                f"Unknown gene-disease validity classification {classification.strip()!r} "
                f"for {gene}. Known classifications are: {allowed}. Refusing to rank an "
                "unrecognised curation tier: a new tier is a change in the source's "
                "curation model and must be placed in CURATED_VALIDITY_RANK deliberately, "
                "not sorted silently to the bottom."
            )
            raise ValueError(msg)
        current = best.get(gene)
        if current is None or (rank[token], panel) < (rank[current[0]], current[1]):
            best[gene] = (token, panel)
    return best


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
    path: Path,
    *,
    restrict_to_genes: Collection[str],
    version: str,
    curated_strengths: Mapping[str, tuple[str, str]],
) -> HpoGenePhenotypeResult:
    """Parse ``genes_to_phenotype.txt``, kept to genes in ``restrict_to_genes``.

    The restriction exists for size, not correctness: the full file annotates
    thousands of genes with no curated clinical-validity record anywhere in this
    pipeline (see ``build.py`` for the measured before/after row counts). A gene in
    ``restrict_to_genes`` with zero matching rows is reported in
    :attr:`HpoGenePhenotypeResult.genes_unmatched`, never silently dropped.

    ``curated_strengths`` supplies the ``association_strength`` column from
    :func:`curated_strength_by_gene`. It is a **separate** argument from
    ``restrict_to_genes`` on purpose: today the two happen to hold the same genes, but
    conflating them would make "which genes do we keep" and "which genes are curated"
    one fact, and a widened restriction would then silently invent curation for the
    genes it added. A gene absent from this mapping is emitted with an EMPTY strength
    and counted in :attr:`HpoGenePhenotypeResult.genes_without_curated_strength`.
    """
    if not path.is_file():
        msg = f"HPO gene-to-phenotype file not found: {path}"
        raise FileNotFoundError(msg)
    restrict_set = frozenset(restrict_to_genes)

    best_by_pair: dict[tuple[str, str], tuple[tuple[int, str], GenePhenotypeRow]] = {}
    total_rows = 0
    rows_kept = 0
    excluded_frequency_rows = 0
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

            frequency_cell = fields[idx["frequency"]]
            if frequency_cell.strip().upper() == EXCLUDED_FREQUENCY_TERM:
                # A negated annotation ("0% of cases") is not an association. Dropped
                # and counted, never rendered as an empty frequency, which would be
                # indistinguishable from "not stated" (GP-14).
                excluded_frequency_rows += 1
                continue

            rows_kept += 1
            genes_matched.add(gene)

            hpo_id = fields[idx["hpo_id"]]
            label = fields[idx["hpo_name"]]
            disease_id = fields[idx["disease_id"]]
            frequency = None if frequency_cell == _NOT_STATED else frequency_cell

            curated = curated_strengths.get(gene)
            candidate = GenePhenotypeRow(
                gene_symbol=gene,
                hpo_id=hpo_id,
                label=label,
                association_strength=None if curated is None else curated[0],
                association_strength_source=None if curated is None else curated[1],
                hpo_frequency=frequency,
            )
            # Prefer a row with real frequency data over one without; break remaining
            # ties on disease_id for a total, deterministic order (GP-30).
            sort_key = (0 if frequency is not None else 1, disease_id)
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
    genes_without_curated_strength = tuple(
        sorted({row.gene_symbol for row in rows if row.association_strength is None})
    )
    rows_without_curated_strength = sum(1 for row in rows if row.association_strength is None)
    return HpoGenePhenotypeResult(
        rows=rows,
        version=version,
        total_rows_read=total_rows,
        rows_kept=rows_kept,
        restrict_gene_count=len(restrict_set),
        genes_matched=len(genes_matched),
        genes_unmatched=genes_unmatched,
        genes_without_curated_strength=genes_without_curated_strength,
        rows_without_curated_strength=rows_without_curated_strength,
        excluded_frequency_rows=excluded_frequency_rows,
    )
