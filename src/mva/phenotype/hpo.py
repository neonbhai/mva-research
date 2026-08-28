"""Gene-to-phenotype knowledge: association records and the lookup index.

This module owns four things:

* :func:`normalise_hpo_id` — the single place in this package that decides whether
  a string is a usable HPO identifier. It mirrors the validator on
  :class:`mva.models.phenotype.PhenotypeObservation` so that a malformed term is
  rejected at the *file* boundary (GP-02) with a message that names the file and
  line, rather than deep inside Pydantic with no positional context.
* :func:`read_tsv_rows` — the shared reader for the comment-prefixed TSV files this
  package consumes. It lives here rather than in ``loader`` so that the dependency
  between the two modules points one way only (``loader`` -> ``hpo``).
* :func:`parse_hpo_frequency` — the parser for HPO's phenotype-**frequency**
  vocabulary, which is a different quantity from association strength and now has
  its own column, its own type and its own tests (ADR 0021).
* :class:`GenePhenotypeIndex` — an immutable, deterministically ordered index over
  :class:`GeneAssociation` records.

**Two quantities, two columns (ADR 0021).** ``knowledge/real/gene_phenotype.tsv``
used to write HPO's frequency vocabulary (``HP:0040283``, ``12/45``, ``50%``) into
a column called ``association_strength``, which this reader validates as curated
gene-disease clinical validity. Two modules agreed on a name and disagreed about
the meaning, and the reader failed closed on the whole file — correctly. The
columns are now separate and mean what they say:

``association_strength``
    Curated gene-disease **clinical validity**: how confident an expert panel is
    that variation in this gene causes disease *at all*. Sourced from ClinGen
    Gene-Disease Validity and EBI Gene2Phenotype (DDG2P), whose vocabularies
    :data:`STRENGTH_WEIGHTS` carries verbatim (case-folded). May be **absent**,
    which means no source classifies the gene — never "weak" (GP-14). It is a
    GENE-level judgement and therefore cannot discriminate between one gene's
    terms; see ASSUMPTION-PHENOTYPE-07 for what that costs and what carries the
    within-gene signal instead.

``hpo_frequency``
    HPO's phenotype frequency: how often the feature occurs *among cases of the
    linked disease*. One of five HP:00402xx terms, an ``n/m`` fraction of
    annotated cases, or a percentage. Optional column; absent means HPO recorded
    ``-`` (not stated), which is unmeasured, not "never occurs".

A frequency token supplied as an ``association_strength`` is refused with a
message that names the confusion — that refusal is the regression guard for the
defect this contract replaced.

**Maturity (GP-20):** the shipped ``knowledge/public/gene_phenotype.tsv`` is a
*synthetic substitute*. The gene symbols are fictional and the association
strengths were written for the demo case. The parsing and indexing logic here is
real; the biology it is pointed at is not.

That claim is now narrower than it used to be, and the narrowing matters. The
*ontology* side of this package — :mod:`mva.phenotype.ontology`,
:mod:`mva.phenotype.corpus`, :mod:`mva.phenotype.similarity` — reads the real
published HPO release and its real annotation corpus, and
``knowledge/real/gene_phenotype.tsv`` is a real HPO table joined to real ClinGen
and DDG2P curations. What remains synthetic is the *demo* gene→term table. No
output of the scorer may be described as biologically valid while it is pointed
at the synthetic file, but the reason is the gene symbols, not the ontology.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Final

from mva.errors import IngestionError

#: Mirrors the pattern enforced by ``PhenotypeObservation.hpo_id``.
HPO_ID_PATTERN: Final[re.Pattern[str]] = re.compile(r"^HP:\d{7}$")

#: Numeric weight per curated gene-disease clinical-validity classification.
#:
#: **These are heuristics, not calibrated parameters (GP-32).** They encode one
#: claim only: that a definitively curated gene-disease link should count for more
#: than a disputed or refuted one. The ratios were chosen to be readable, not
#: fitted to any labelled dataset. Changing them requires a decision record and a
#: before/after comparison; the current values and their before/after are recorded
#: in ADR 0021.
#:
#: The vocabulary is the union of the two real curation sources and the synthetic
#: demo table's own scale, case-folded so that ClinGen's ``Definitive`` and DDG2P's
#: ``definitive`` are one token (a case fold is not a semantic remap):
#:
#: * ClinGen Gene-Disease Validity: definitive, strong, moderate, limited,
#:   disputed, refuted, no known disease relationship.
#: * DDG2P confidence: definitive, strong, moderate, limited, disputed, refuted.
#: * ``supporting`` exists only in ``knowledge/public/gene_phenotype.tsv``, the
#:   synthetic demo table; no real source emits it. It keeps its original 0.4 so
#:   the golden expectation is untouched by ADR 0021.
#:
#: Note the floor is 0.1, not 0.0: a refuted gene-disease claim is still a real
#: curated statement and the HPO annotation hanging off it is still real. Zero
#: would silently erase both, and GP-13/GP-19 say contradicted candidates are
#: down-ranked and kept, never deleted.
STRENGTH_WEIGHTS: dict[str, float] = {
    "definitive": 1.0,
    "strong": 0.8,
    "moderate": 0.6,
    "limited": 0.4,
    "supporting": 0.4,
    "disputed": 0.2,
    "refuted": 0.1,
    "no known disease relationship": 0.1,
}

#: Weight for an association whose gene carries NO curated validity classification.
#:
#: Absence of curation is not a weak curation (GP-14), so this value is deliberately
#: **not on the curated ladder above**: reading 0.5 in an evidence item's
#: ``numeric_value`` is itself the signal that nobody has classified this gene.
#: ``test_uncurated_weight_is_not_mistakable_for_a_curated_one`` enforces that.
#:
#: The two tempting defaults are both wrong and both silent. Zero erases a real HPO
#: annotation (and would drop the gene out of the score's denominator entirely);
#: 1.0 invents a definitive expert classification that nobody made. Every call site
#: names this constant explicitly via :meth:`GeneAssociation.weight_or` — there is
#: no default argument, so the choice cannot be made by omission. See ADR 0021.
UNCURATED_ASSOCIATION_WEIGHT: Final[float] = 0.5

#: Columns every gene-phenotype knowledge file must provide.
GENE_PHENOTYPE_COLUMNS: Final[tuple[str, ...]] = (
    "gene_symbol",
    "hpo_id",
    "label",
    "association_strength",
    "source",
    "version",
)

#: Columns a gene-phenotype knowledge file MAY provide, and what absence means.
#:
#: ``hpo_frequency`` is optional rather than required so that a table carrying only
#: curated strength (the synthetic demo table) and a table carrying both (the real
#: HPO-derived table) are both readable by one reader. A file without the column is
#: not a file whose phenotypes never occur; it is a file that does not state how
#: often they do.
OPTIONAL_GENE_PHENOTYPE_COLUMNS: Final[tuple[str, ...]] = (
    "association_strength_source",
    "hpo_frequency",
)


class HpoFrequencyKind(StrEnum):
    """Which of HPO's three frequency spellings a value used.

    Kept distinct rather than normalised away to a single number: ``HP:0040282``
    ("Frequent", 30-79% of cases) and ``12/45`` (a counted 26.7% in one annotated
    series) are different epistemic objects, and collapsing them would hide that
    one is a curator's band and the other is a tally.
    """

    TERM = "term"
    """One of the five HP:00402xx frequency-subontology terms."""

    FRACTION = "fraction"
    """``n/m``: n of m annotated cases had the feature."""

    PERCENTAGE = "percentage"
    """A percentage recorded directly in the source."""


#: HPO's frequency subontology, quoted from ``hp.obo``, with the proportion bounds
#: those definitions state. Bounds are the SOURCE's own, never interpolated.
HPO_FREQUENCY_TERMS: Final[Mapping[str, tuple[str, float, float]]] = {
    "HP:0040280": ("Obligate", 1.0, 1.0),
    "HP:0040281": ("Very frequent", 0.80, 0.99),
    "HP:0040282": ("Frequent", 0.30, 0.79),
    "HP:0040283": ("Occasional", 0.05, 0.29),
    "HP:0040284": ("Very rare", 0.01, 0.04),
}

#: HPO's "Excluded" (0% of cases) frequency term.
#:
#: Deliberately NOT in :data:`HPO_FREQUENCY_TERMS`. It is a *negated* annotation —
#: "this feature is not part of this disease" — and a gene-phenotype association
#: table states the opposite claim. Carrying it here would turn an exclusion into
#: an association weighted like any other, so it is refused by name with a message
#: saying where such a row belongs instead.
EXCLUDED_FREQUENCY_TERM: Final = "HP:0040285"

_FRACTION_PATTERN: Final[re.Pattern[str]] = re.compile(r"^(?P<numerator>\d+)/(?P<denominator>\d+)$")
_PERCENTAGE_PATTERN: Final[re.Pattern[str]] = re.compile(r"^(?P<value>\d+(?:\.\d+)?)%$")


@dataclass(frozen=True, slots=True)
class HpoFrequency:
    """How often a phenotype occurs among cases of the linked disease.

    ``lower_bound``/``upper_bound`` are proportions in ``[0, 1]``. For a term they
    are the band its ``hp.obo`` definition states; for a fraction or a percentage
    they are equal, because a counted value is a point, not a range. ``raw`` is
    kept so the source's own spelling survives into provenance unaltered.
    """

    raw: str
    kind: HpoFrequencyKind
    label: str
    lower_bound: float
    upper_bound: float


def looks_like_hpo_frequency(token: str) -> bool:
    """True if ``token`` is HPO frequency vocabulary rather than a validity class.

    Exists so that the one confusion this contract was built to prevent produces a
    *diagnostic* refusal instead of a generic "unknown value" — the reader may be
    an agent whose only view of the rule is the error string.
    """
    candidate = token.strip().upper().replace("HP_", "HP:")
    if candidate in HPO_FREQUENCY_TERMS or candidate == EXCLUDED_FREQUENCY_TERM:
        return True
    stripped = token.strip()
    return bool(_FRACTION_PATTERN.match(stripped) or _PERCENTAGE_PATTERN.match(stripped))


def parse_hpo_frequency(raw: str, *, context: str) -> HpoFrequency | None:
    """Parse one HPO frequency cell, or ``None`` when the source did not state one.

    Accepts exactly what HPO emits and nothing else:

    * one of the five terms in :data:`HPO_FREQUENCY_TERMS`;
    * ``n/m`` — n of m annotated cases, with ``m > 0`` and ``n <= m``;
    * ``p%`` — a percentage in ``[0, 100]``.

    An empty cell (HPO's ``-``, rendered empty by the generator) is ``None``: not
    stated is unmeasured, not zero (GP-14). Anything else raises rather than
    degrading to ``None``, because a silently-dropped frequency is indistinguishable
    from one the source never recorded, and the two mean different things.
    """
    token = raw.strip()
    if not token:
        return None

    identifier = token.upper().replace("HP_", "HP:")
    if identifier == EXCLUDED_FREQUENCY_TERM:
        msg = (
            f"HPO frequency {token!r} in {context} is HP:0040285 (Excluded, 0% of cases): "
            "a NEGATED annotation stating the feature is not part of the disease. A "
            "gene-phenotype association table asserts the opposite, so this row must be "
            "dropped by the generator (and counted in the table's header), not carried "
            "here where it would be weighted as an association (GP-14, GP-19)."
        )
        raise IngestionError(msg)

    term = HPO_FREQUENCY_TERMS.get(identifier)
    if term is not None:
        label, lower, upper = term
        return HpoFrequency(
            raw=token,
            kind=HpoFrequencyKind.TERM,
            label=label,
            lower_bound=lower,
            upper_bound=upper,
        )

    fraction = _FRACTION_PATTERN.match(token)
    if fraction is not None:
        numerator = int(fraction.group("numerator"))
        denominator = int(fraction.group("denominator"))
        if denominator == 0:
            msg = (
                f"HPO frequency {token!r} in {context} has a zero denominator; "
                "'n of 0 cases' is not a proportion."
            )
            raise IngestionError(msg)
        if numerator > denominator:
            msg = (
                f"HPO frequency {token!r} in {context} counts more affected cases than "
                "annotated cases; refusing to clamp a value the source got wrong."
            )
            raise IngestionError(msg)
        value = numerator / denominator
        return HpoFrequency(
            raw=token,
            kind=HpoFrequencyKind.FRACTION,
            label=f"{numerator} of {denominator} annotated cases",
            lower_bound=value,
            upper_bound=value,
        )

    percentage = _PERCENTAGE_PATTERN.match(token)
    if percentage is not None:
        percent = float(percentage.group("value"))
        if percent > 100.0:
            msg = f"HPO frequency {token!r} in {context} exceeds 100%."
            raise IngestionError(msg)
        value = percent / 100.0
        return HpoFrequency(
            raw=token,
            kind=HpoFrequencyKind.PERCENTAGE,
            label=f"{percentage.group('value')}% of cases",
            lower_bound=value,
            upper_bound=value,
        )

    allowed = ", ".join(sorted(HPO_FREQUENCY_TERMS))
    msg = (
        f"Unrecognised hpo_frequency {token!r} in {context}. Allowed: one of the HPO "
        f"frequency terms ({allowed}), an 'n/m' fraction of annotated cases, a "
        "percentage such as '50%', or an EMPTY cell meaning the source recorded '-' "
        "(not stated). Note this column carries FREQUENCY, not curated validity — a "
        f"value such as 'definitive' belongs in association_strength."
    )
    raise IngestionError(msg)


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
    """One gene -> HPO term association, with its two independent qualifiers.

    Self-validating and self-normalising, so an association can never enter the
    index with an unrecognised strength (which would have no weight) or a
    malformed HPO identifier (which would never match a patient observation and
    would therefore quietly deflate the gene's score).
    """

    gene_symbol: str
    hpo_id: str
    label: str
    association_strength: str | None
    """Curated gene-disease clinical validity — a key of :data:`STRENGTH_WEIGHTS`.

    ``None`` means no curation source classifies this gene. That is a real and
    common state, and it is not the same claim as ``limited`` (GP-14): see
    :data:`UNCURATED_ASSOCIATION_WEIGHT` for what it is worth and why the caller
    has to say so.
    """
    source: str
    version: str
    hpo_frequency: HpoFrequency | None = None
    """How often the feature occurs among cases of the linked disease, when the
    source stated it. A *different quantity* from ``association_strength``; ADR
    0021 exists because the two once shared a column."""
    association_strength_source: str | None = None
    """Which curation panel supplied ``association_strength`` (ClinGen, DDG2P, ...).

    Distinct from ``source``, which names where the gene->term ANNOTATION came from.
    In the real table those are genuinely different bodies -- HPO annotates the term,
    ClinGen and DDG2P classify the gene -- and one provenance field covering both
    would misattribute the validity call to HPO (GP-31)."""

    def __post_init__(self) -> None:
        symbol = self.gene_symbol.strip()
        if not symbol:
            msg = "Gene-phenotype association has an empty gene_symbol."
            raise IngestionError(msg)

        object.__setattr__(self, "gene_symbol", symbol)
        object.__setattr__(self, "association_strength", _normalise_strength(self, symbol=symbol))
        object.__setattr__(self, "hpo_id", normalise_hpo_id(self.hpo_id, context=f"gene {symbol}"))
        object.__setattr__(self, "label", self.label.strip())
        object.__setattr__(self, "source", self.source.strip())
        object.__setattr__(self, "version", self.version.strip())
        strength_source = (self.association_strength_source or "").strip()
        if strength_source and self.association_strength is None:
            msg = (
                f"{symbol}/{self.hpo_id} names a curation panel "
                f"({strength_source!r}) but carries no association_strength. A "
                "provenance for a value that does not exist is not a partial record, "
                "it is a contradiction: either the classification was dropped in "
                "transit, or the panel name was. Supply both, or neither."
            )
            raise IngestionError(msg)
        object.__setattr__(self, "association_strength_source", strength_source or None)

    @property
    def weight(self) -> float | None:
        """Numeric weight for this association's curated strength, or ``None``.

        Optional on purpose. ``None`` is not a number a caller can accidentally add
        to a sum — pyright refuses it — so an uncurated association cannot slip into
        a score with an implied value. Use :meth:`weight_or` and name the value.
        """
        if self.association_strength is None:
            return None
        return STRENGTH_WEIGHTS[self.association_strength]

    def weight_or(self, uncurated: float) -> float:
        """Weight, with the caller stating what an *uncurated* association is worth.

        There is no default argument, and that is the point (ADR 0021): the value
        an unknown curation status contributes is a scientific choice, so it is
        made visibly at the call site — ``assoc.weight_or(UNCURATED_ASSOCIATION_WEIGHT)``
        — and is greppable, rather than being decided by omission inside a property.
        """
        weight = self.weight
        if weight is None:
            return uncurated
        return weight

    @property
    def sort_key(self) -> tuple[str, str, str]:
        """Total order used everywhere this record is emitted (GP-30)."""
        return (self.gene_symbol.upper(), self.hpo_id, self.gene_symbol)


def _normalise_strength(assoc: GeneAssociation, *, symbol: str) -> str | None:
    """Case-fold and validate one clinical-validity token, or accept its absence.

    Empty is accepted and means *absent*; anything non-empty must be a curated
    classification. The frequency-vocabulary branch is the regression guard for the
    defect ADR 0021 replaced: it names the conflation instead of reporting a
    generic unknown value, because that message is the only view of this rule a
    reader debugging a rejected table will get.
    """
    raw = (assoc.association_strength or "").strip()
    if not raw:
        return None

    strength = " ".join(raw.lower().split())
    if strength in STRENGTH_WEIGHTS:
        return strength

    allowed = ", ".join(sorted(STRENGTH_WEIGHTS))
    if looks_like_hpo_frequency(raw):
        msg = (
            f"association_strength {raw!r} for {symbol} is HPO phenotype-FREQUENCY "
            "vocabulary, not gene-disease clinical validity. Frequency ('how often does "
            "this feature occur in cases of this disease') and validity ('how confident "
            "are we that this gene causes disease at all') are different quantities and "
            "must not share a column (ADR 0021). Put it in the 'hpo_frequency' column; "
            f"association_strength takes one of: {allowed}, or an EMPTY cell when no "
            "curation source classifies the gene."
        )
        raise IngestionError(msg)

    msg = (
        f"Unknown association_strength {raw!r} for {symbol}; allowed values are: "
        f"{allowed}, or an EMPTY cell meaning no curation source classifies this gene. "
        "An unrecognised strength has no weight, so accepting it would silently drop "
        "the association."
    )
    raise IngestionError(msg)


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
        """Load associations from a ``#``-commented TSV (GP-02: typed, not dicts).

        ``hpo_frequency`` is read when the file provides it and left ``None`` when it
        does not (:data:`OPTIONAL_GENE_PHENOTYPE_COLUMNS`), so one reader serves both
        the synthetic demo table and the real HPO-derived one.
        """
        associations = [
            GeneAssociation(
                gene_symbol=row["gene_symbol"],
                hpo_id=normalise_hpo_id(row["hpo_id"], context=f"{path.name} line {lineno}"),
                label=row["label"],
                association_strength=row["association_strength"],
                source=row["source"],
                version=row["version"],
                hpo_frequency=parse_hpo_frequency(
                    row.get("hpo_frequency", ""), context=f"{path.name} line {lineno}"
                ),
                association_strength_source=row.get("association_strength_source", ""),
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
