"""Gene-to-phenotype knowledge: association records and the lookup index.

This module owns three things:

* :func:`normalise_hpo_id` — the single place in this package that decides whether
  a string is a usable HPO identifier. It mirrors the validator on
  :class:`mva.models.phenotype.PhenotypeObservation` so that a malformed term is
  rejected at the *file* boundary (GP-02) with a message that names the file and
  line, rather than deep inside Pydantic with no positional context.
* :func:`read_tsv_rows` — the shared reader for the comment-prefixed TSV files this
  package consumes. It lives here rather than in ``loader`` so that the dependency
  between the two modules points one way only (``loader`` -> ``hpo``).
* :class:`GenePhenotypeIndex` — an immutable, deterministically ordered index over
  :class:`GeneAssociation` records.

**Maturity (GP-20):** the shipped ``knowledge/public/gene_phenotype.tsv`` is a
*synthetic substitute*. The gene symbols are fictional and the association
strengths were written for the demo case. The parsing and indexing logic here is
real; the biology it is pointed at is not.

That claim is now narrower than it used to be, and the narrowing matters. The
*ontology* side of this package — :mod:`mva.phenotype.ontology`,
:mod:`mva.phenotype.corpus`, :mod:`mva.phenotype.similarity` — reads the real
published HPO release and its real annotation corpus. What remains synthetic is
this file's gene→term table. No output of the scorer may be described as
biologically valid until a real gene-phenotype knowledge base replaces it, but
the reason is the gene symbols, not the ontology.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from mva.errors import IngestionError

#: Mirrors the pattern enforced by ``PhenotypeObservation.hpo_id``.
HPO_ID_PATTERN: Final[re.Pattern[str]] = re.compile(r"^HP:\d{7}$")

#: Numeric weight per curated association strength.
#:
#: **These are heuristics, not calibrated parameters (GP-32).** They encode one
#: claim only: that a definitively curated gene-phenotype link should count for
#: more than a single supporting case report. The ratios were chosen to be
#: readable (1.0 / 0.8 / 0.6 / 0.4), not fitted to any labelled dataset. Changing
#: them requires a decision record and a before/after comparison.
#:
#: Note the floor is 0.4, not 0.0: a "supporting" association is still a real
#: curated statement. Zero would silently erase it.
STRENGTH_WEIGHTS: dict[str, float] = {
    "definitive": 1.0,
    "strong": 0.8,
    "moderate": 0.6,
    "supporting": 0.4,
}

#: Columns every gene-phenotype knowledge file must provide.
GENE_PHENOTYPE_COLUMNS: Final[tuple[str, ...]] = (
    "gene_symbol",
    "hpo_id",
    "label",
    "association_strength",
    "source",
    "version",
)


def normalise_hpo_id(raw: str, *, context: str) -> str:
    """Return a canonical ``HP:0000000`` identifier, or raise.

    Accepts the two spellings that appear in real exports (``HP:0000252`` and the
    OBO-style ``HP_0000252``) and nothing else. In particular a bare integer or a
    truncated identifier such as ``HP:123`` is an error: silently zero-padding it
    would fabricate a term that the clinician never recorded.

    ``context`` is a caller-supplied locator (file and line, or gene symbol) used
    only to make the message actionable. It must contain identifiers and counts
    only, never record text (PRIV-09).
    """
    token = raw.strip().upper().replace("HP_", "HP:")
    if not HPO_ID_PATTERN.match(token):
        msg = (
            f"Invalid HPO identifier {raw.strip()!r} in {context}; expected the form "
            "'HP:0001250' — the prefix 'HP:' (or 'HP_') followed by exactly seven digits."
        )
        raise IngestionError(msg)
    return token


def read_tsv_rows(
    path: Path,
    *,
    required_columns: Sequence[str],
) -> list[tuple[int, dict[str, str]]]:
    """Read a ``#``-commented TSV into ``(line_number, row)`` pairs.

    Deliberately strict, because the failure mode of a lenient tabular reader is a
    column silently shifting by one and every downstream status being wrong:

    * the first non-comment, non-blank line is the header;
    * every required column must be present, or the file is rejected by name;
    * a row whose field count differs from the header is rejected with its line
      number rather than padded or truncated.

    Extra columns are tolerated and preserved — a knowledge file gaining an
    annotation column should not break a consumer that does not read it.
    """
    if not path.is_file():
        msg = f"Phenotype input not found: {path.as_posix()}"
        raise IngestionError(msg)

    header: list[str] | None = None
    rows: list[tuple[int, dict[str, str]]] = []

    text = path.read_text(encoding="utf-8-sig")
    for lineno, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if header is None:
            header = [name.strip() for name in fields]
            _assert_columns_present(path, header, required_columns)
            continue
        if len(fields) != len(header):
            msg = (
                f"{path.name} line {lineno}: expected {len(header)} tab-separated fields "
                f"(from the header), found {len(fields)}. Refusing to guess which column "
                "is missing; a shifted column silently corrupts every status in the file."
            )
            raise IngestionError(msg)
        rows.append(
            (lineno, {name: value.strip() for name, value in zip(header, fields, strict=True)})
        )

    if header is None:
        msg = f"{path.name} contains no header row (only comments or blank lines)."
        raise IngestionError(msg)
    return rows


def _assert_columns_present(
    path: Path, header: Sequence[str], required_columns: Sequence[str]
) -> None:
    missing = [name for name in required_columns if name not in header]
    if missing:
        msg = (
            f"{path.name} is missing required column(s): {', '.join(sorted(missing))}. "
            f"Required columns are: {', '.join(required_columns)}."
        )
        raise IngestionError(msg)


@dataclass(frozen=True, slots=True)
class GeneAssociation:
    """One curated gene -> HPO term association.

    Self-validating and self-normalising, so an association can never enter the
    index with an unrecognised strength (which would have no weight) or a
    malformed HPO identifier (which would never match a patient observation and
    would therefore quietly deflate the gene's score).
    """

    gene_symbol: str
    hpo_id: str
    label: str
    association_strength: str
    """One of ``definitive``, ``strong``, ``moderate``, ``supporting``."""
    source: str
    version: str

    def __post_init__(self) -> None:
        symbol = self.gene_symbol.strip()
        if not symbol:
            msg = "Gene-phenotype association has an empty gene_symbol."
            raise IngestionError(msg)

        strength = self.association_strength.strip().lower()
        if strength not in STRENGTH_WEIGHTS:
            allowed = ", ".join(sorted(STRENGTH_WEIGHTS))
            msg = (
                f"Unknown association_strength {self.association_strength.strip()!r} for "
                f"{symbol}; allowed values are: {allowed}. An unrecognised strength has no "
                "weight, so accepting it would silently drop the association."
            )
            raise IngestionError(msg)

        object.__setattr__(self, "gene_symbol", symbol)
        object.__setattr__(self, "association_strength", strength)
        object.__setattr__(self, "hpo_id", normalise_hpo_id(self.hpo_id, context=f"gene {symbol}"))
        object.__setattr__(self, "label", self.label.strip())
        object.__setattr__(self, "source", self.source.strip())
        object.__setattr__(self, "version", self.version.strip())

    @property
    def weight(self) -> float:
        """Numeric weight for this association's curated strength."""
        return STRENGTH_WEIGHTS[self.association_strength]

    @property
    def sort_key(self) -> tuple[str, str, str]:
        """Total order used everywhere this record is emitted (GP-30)."""
        return (self.gene_symbol.upper(), self.hpo_id, self.gene_symbol)


class GenePhenotypeIndex:
    """Immutable lookup over gene-phenotype associations.

    Every accessor returns a ``tuple`` in a total, documented order. Nothing here
    iterates a ``set`` or a ``dict`` whose insertion order depends on the input
    file's order, because artifacts derived from this index must be byte-identical
    across runs (GP-30).
    """

    __slots__ = ("_associations", "_by_gene", "_by_term", "_version")

    def __init__(self, associations: Sequence[GeneAssociation], *, version: str) -> None:
        ordered = tuple(sorted(associations, key=lambda assoc: assoc.sort_key))

        by_gene: dict[str, list[GeneAssociation]] = {}
        by_term: dict[str, list[str]] = {}
        seen: set[tuple[str, str]] = set()

        for assoc in ordered:
            key = (assoc.gene_symbol.upper(), assoc.hpo_id)
            if key in seen:
                msg = (
                    f"Duplicate gene-phenotype association {assoc.gene_symbol}/{assoc.hpo_id}. "
                    "Duplicates would be counted twice in the score's weighted denominator; "
                    "de-duplicate the knowledge file or merge the sources explicitly."
                )
                raise IngestionError(msg)
            seen.add(key)
            by_gene.setdefault(key[0], []).append(assoc)
            by_term.setdefault(assoc.hpo_id, []).append(assoc.gene_symbol)

        self._associations: tuple[GeneAssociation, ...] = ordered
        self._by_gene: dict[str, tuple[GeneAssociation, ...]] = {
            gene: tuple(items) for gene, items in by_gene.items()
        }
        self._by_term: dict[str, tuple[str, ...]] = {
            term: tuple(sorted(set(genes), key=lambda name: (name.upper(), name)))
            for term, genes in by_term.items()
        }
        self._version: str = version

    @classmethod
    def from_tsv(cls, path: Path, *, version: str) -> GenePhenotypeIndex:
        """Load associations from a ``#``-commented TSV (GP-02: typed, not dicts)."""
        associations = [
            GeneAssociation(
                gene_symbol=row["gene_symbol"],
                hpo_id=normalise_hpo_id(row["hpo_id"], context=f"{path.name} line {lineno}"),
                label=row["label"],
                association_strength=row["association_strength"],
                source=row["source"],
                version=row["version"],
            )
            for lineno, row in read_tsv_rows(path, required_columns=GENE_PHENOTYPE_COLUMNS)
        ]
        return cls(associations, version=version)

    def terms_for_gene(self, gene_symbol: str) -> tuple[GeneAssociation, ...]:
        """Associations for one gene, ordered by HPO id. Empty tuple if unknown.

        An unknown gene yields ``()`` rather than raising: "this pipeline holds no
        phenotype knowledge about the gene" is a normal state, and the scorer turns
        it into a neutral score, never into a penalty (GP-14).
        """
        return self._by_gene.get(gene_symbol.strip().upper(), ())

    def genes_for_term(self, hpo_id: str) -> tuple[str, ...]:
        """Gene symbols associated with one HPO term, ordered. Empty if unknown."""
        return self._by_term.get(hpo_id.strip().upper().replace("HP_", "HP:"), ())

    @property
    def version(self) -> str:
        """Release identifier for this knowledge snapshot (cited in evidence)."""
        return self._version

    @property
    def associations(self) -> tuple[GeneAssociation, ...]:
        """Every association, in total order."""
        return self._associations

    @property
    def gene_symbols(self) -> tuple[str, ...]:
        """Distinct gene symbols in the index, in total order."""
        return tuple(
            sorted(
                {assoc.gene_symbol for assoc in self._associations},
                key=lambda name: (name.upper(), name),
            )
        )

    def __len__(self) -> int:
        return len(self._associations)

    def __repr__(self) -> str:
        return (
            f"GenePhenotypeIndex(version={self._version!r}, "
            f"genes={len(self._by_gene)}, associations={len(self._associations)})"
        )
