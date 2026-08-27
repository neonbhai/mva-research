"""The Human Phenotype Ontology as a directed acyclic graph.

This module exists because HPO is **a DAG, not a tree**, and almost every cheap
implementation gets that wrong. ``HP:0000252`` (Microcephaly) reaches the root
through several distinct paths of different lengths; ``HP:0001305``
(Dandy-Walker malformation) has more than one ``is_a`` parent. Code that reads a
single parent, or that assumes one path implies one depth, silently loses most of
the ontology's structure and then reports the loss as "no match".

Three consequences are designed for here:

* **Multiple parents.** :meth:`HpoOntology.parents` returns every ``is_a`` target,
  and the closure walks all of them.
* **Multiple paths / re-convergence.** The closure is a visited-set traversal, so
  a term reachable by five paths is visited once and the cost stays linear in the
  reachable subgraph rather than exponential in the path count.
* **Cycles.** The release is expected to be acyclic, but a malformed or
  hand-edited file is not, and an ontology loader that recurses off a cycle takes
  the whole pipeline with it. The traversal is iterative with an explicit visited
  set, so a cycle terminates instead of exhausting the stack.

**Direction matters and is not symmetric.** :meth:`ancestor_closure` walks *up*
(toward the root) and :meth:`descendant_closure` walks *down*. The phenotype
scorer uses them for opposite purposes — see :mod:`mva.phenotype.propagation` —
and swapping them inverts the clinical meaning of every negative finding.

**No module-level state.** Every cache lives on an explicitly constructed
:class:`HpoOntology` instance. A module-level dict or an ``lru_cache`` on a
parsing function would make the ontology a hidden global, break the layering rule
this repository tests for, and leak one case's data into the next.

**No network (PRIV-05).** This module parses a local file that the composition
root supplies by path. It never fetches anything.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from mva.determinism import hash_file
from mva.errors import IngestionError
from mva.phenotype.hpo import HPO_ID_PATTERN

#: OBO header key carrying the ontology release. The value is recorded verbatim in
#: provenance: "the HPO release we scored against" is not a claim anyone should
#: have to reconstruct from a file modification time.
DATA_VERSION_KEY: Final[str] = "data-version"

#: Stanza header that opens a term. Every other stanza — ``[Typedef]`` in current
#: releases — is skipped: they define relation types, not phenotype terms, and
#: admitting them would put non-terms in the graph.
_TERM_STANZA: Final[str] = "[Term]"


@dataclass(frozen=True, slots=True)
class HpoTerm:
    """One HPO term and its outgoing ``is_a`` edges.

    ``parents`` is a tuple, not a single value, because a term routinely has
    several. It is stored sorted so that anything derived from it is byte-stable
    across runs (GP-30).
    """

    hpo_id: str
    name: str
    parents: tuple[str, ...]
    alt_ids: tuple[str, ...]
    is_obsolete: bool
    replaced_by: str | None


class HpoOntology:
    """Immutable, deterministically ordered view of one HPO release.

    Construction parses the OBO file once; the ancestor closure is memoised on the
    instance as it is requested. The descendant closure is memoised too but is
    computed only on demand, because the reflexive descendant set of a term near
    the root is most of the ontology and materialising all of them would cost far
    more memory than the handful of excluded terms a real profile contains.

    Obsolete terms are **not** part of the graph. They are kept in a side table so
    that :meth:`resolve` can forward an obsoleted identifier to its ``replaced_by``
    successor where the release states one, and return ``None`` where it does not.
    Returning ``None`` — rather than inventing a placeholder node — is what lets
    the caller report "this term is not in this release" instead of scoring it as
    a term that matches nothing.
    """

    __slots__ = (
        "_alt_to_primary",
        "_ancestor_cache",
        "_children",
        "_data_version",
        "_descendant_cache",
        "_obsolete",
        "_source_path",
        "_source_sha256",
        "_terms",
    )

    def __init__(
        self,
        terms: Sequence[HpoTerm],
        *,
        data_version: str,
        source_path: Path | None = None,
        source_sha256: str | None = None,
    ) -> None:
        primary: dict[str, HpoTerm] = {}
        obsolete: dict[str, HpoTerm] = {}
        for term in terms:
            target = obsolete if term.is_obsolete else primary
            if term.hpo_id in target:
                msg = (
                    f"Duplicate HPO term {term.hpo_id} in the ontology source. A term "
                    "defined twice has two parent sets, and no reader can pick between "
                    "them; fix the source file rather than merging them silently."
                )
                raise IngestionError(msg)
            target[term.hpo_id] = term

        # Alternate identifiers are how the release records historical merges. A
        # profile written against an older release will contain them, and dropping
        # them would silently delete real observations.
        alt_to_primary: dict[str, str] = {}
        for term in primary.values():
            for alt in term.alt_ids:
                alt_to_primary.setdefault(alt, term.hpo_id)

        children: dict[str, list[str]] = {}
        for term in primary.values():
            for parent in term.parents:
                if parent in primary:
                    children.setdefault(parent, []).append(term.hpo_id)

        self._terms: dict[str, HpoTerm] = primary
        self._obsolete: dict[str, HpoTerm] = obsolete
        self._alt_to_primary: dict[str, str] = alt_to_primary
        self._children: dict[str, tuple[str, ...]] = {
            parent: tuple(sorted(kids)) for parent, kids in children.items()
        }
        self._data_version: str = data_version
        self._source_path: Path | None = source_path
        self._source_sha256: str | None = source_sha256
        self._ancestor_cache: dict[str, tuple[str, ...]] = {}
        self._descendant_cache: dict[str, tuple[str, ...]] = {}

    # -- construction -------------------------------------------------------

    @classmethod
    def from_obo(
        cls,
        path: Path,
        *,
        expected_sha256: str | None = None,
    ) -> HpoOntology:
        """Parse an ``hp.obo`` release.

        ``expected_sha256`` is a caller-supplied pin. It is a constructor argument
        rather than a constant in this module because the composition root owns
        resource provenance; library code that hardcodes an absolute path or a
        digest cannot be pointed at a different release without an edit.

        The digest is always computed and recorded even when nothing is pinned, so
        provenance can state exactly which bytes were scored against.
        """
        if not path.is_file():
            msg = (
                f"HPO ontology file not found: {path.as_posix()}. The phenotype stage "
                "reads a local release only (PRIV-05); download it in a separate "
                "offline acquisition step and pass the path in from the composition root."
            )
            raise IngestionError(msg)

        digest = hash_file(path)
        if expected_sha256 is not None and digest != expected_sha256.strip().lower():
            msg = (
                f"HPO ontology {path.name} has sha256 {digest}, but {expected_sha256.strip()} "
                "was pinned. Scoring against an unexpected ontology release silently "
                "changes every phenotype score; re-pin the digest deliberately or restore "
                "the pinned release."
            )
            raise IngestionError(msg)

        terms, data_version = _parse_obo(path)
        if not terms:
            msg = f"{path.name} contains no [Term] stanzas; it is not an OBO ontology."
            raise IngestionError(msg)
        return cls(
            terms,
            data_version=data_version,
            source_path=path,
            source_sha256=digest,
        )

    # -- identity and provenance -------------------------------------------

    @property
    def data_version(self) -> str:
        """The release string from the OBO header, verbatim (e.g. ``hp/releases/2026-06-23``)."""
        return self._data_version

    @property
    def release(self) -> str:
        """Trailing component of :attr:`data_version` (e.g. ``2026-06-23``).

        Convenience only. :attr:`data_version` is the value that goes into
        provenance, because it is the value the release actually declares.
        """
        return self._data_version.rsplit("/", 1)[-1]

    @property
    def source_sha256(self) -> str | None:
        """Digest of the parsed file, for provenance. ``None`` for in-memory graphs."""
        return self._source_sha256

    @property
    def source_path(self) -> Path | None:
        return self._source_path

    # -- lookups ------------------------------------------------------------

    def resolve(self, hpo_id: str) -> str | None:
        """Canonical primary identifier for ``hpo_id``, or ``None`` if this release has none.

        Handles the three cases a real profile produces:

        * a current primary identifier — returned unchanged;
        * an ``alt_id`` from a historical merge — forwarded to its primary term;
        * an obsoleted term — forwarded to ``replaced_by`` when the release states
          one, otherwise ``None``.

        ``None`` means "this release does not contain this term". It deliberately
        does not mean "no match": the caller must report an unresolvable term as an
        information gap, never as a failed comparison (GP-14).
        """
        token = hpo_id.strip().upper().replace("HP_", "HP:")
        if token in self._terms:
            return token
        primary = self._alt_to_primary.get(token)
        if primary is not None:
            return primary
        dead = self._obsolete.get(token)
        if dead is not None and dead.replaced_by is not None:
            return self.resolve(dead.replaced_by)
        return None

    def term(self, hpo_id: str) -> HpoTerm | None:
        """The term record for ``hpo_id`` after alias resolution, or ``None``."""
        resolved = self.resolve(hpo_id)
        return self._terms.get(resolved) if resolved is not None else None

    def label(self, hpo_id: str) -> str | None:
        term = self.term(hpo_id)
        return term.name if term is not None else None

    def contains(self, hpo_id: str) -> bool:
        return self.resolve(hpo_id) is not None

    def parents(self, hpo_id: str) -> tuple[str, ...]:
        """Direct ``is_a`` parents, sorted. Empty for a root or an unknown term.

        A tuple because HPO terms have several parents; a scalar accessor here is
        the single commonest way ontology code loses most of the graph.
        """
        term = self.term(hpo_id)
        if term is None:
            return ()
        return tuple(parent for parent in term.parents if parent in self._terms)

    def children(self, hpo_id: str) -> tuple[str, ...]:
        """Direct ``is_a`` children, sorted. Empty for a leaf or an unknown term."""
        resolved = self.resolve(hpo_id)
        if resolved is None:
            return ()
        return self._children.get(resolved, ())

    # -- closures -----------------------------------------------------------

    def ancestor_closure(self, hpo_id: str) -> tuple[str, ...]:
        """Every term reachable by following ``is_a`` **upward**, including ``hpo_id`` itself.

        Reflexive on purpose: the two consumers — "does this gene's term subsume an
        observed feature?" and "what is the most informative common ancestor?" —
        both need the term itself to be a candidate, and a non-reflexive closure
        makes ``sim(t, t)`` a special case that is easy to forget.

        Returns ``()`` for a term this release does not contain, which is
        distinguishable from ``(hpo_id,)`` for a known isolated root.

        Iterative with a visited set: multi-parent DAGs re-converge constantly, and
        a naive recursive walk both re-explores exponentially and never terminates
        on a cyclic file.
        """
        resolved = self.resolve(hpo_id)
        if resolved is None:
            return ()
        cached = self._ancestor_cache.get(resolved)
        if cached is not None:
            return cached
        closure = self._walk(resolved, self._parents_of)
        self._ancestor_cache[resolved] = closure
        return closure

    def descendant_closure(self, hpo_id: str) -> tuple[str, ...]:
        """Every term reachable by following ``is_a`` **downward**, including ``hpo_id``.

        Computed on demand rather than eagerly for every term: the reflexive
        descendant set of a term near the root is most of the ontology, so
        materialising all 20k of them would cost orders of magnitude more memory
        than the few excluded terms a clinical profile actually contains.
        """
        resolved = self.resolve(hpo_id)
        if resolved is None:
            return ()
        cached = self._descendant_cache.get(resolved)
        if cached is not None:
            return cached
        closure = self._walk(resolved, self._children_of)
        self._descendant_cache[resolved] = closure
        return closure

    def _walk(self, start: str, step: Callable[[str], tuple[str, ...]]) -> tuple[str, ...]:
        """Iterative reflexive closure over one edge direction.

        ``step`` is passed as a bound method rather than a boolean flag, so the two
        directions cannot be confused by a stray argument at the call site.
        """
        seen: set[str] = {start}
        stack: list[str] = [start]
        while stack:
            current = stack.pop()
            for neighbour in step(current):
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
        return tuple(sorted(seen))

    def _parents_of(self, hpo_id: str) -> tuple[str, ...]:
        term = self._terms.get(hpo_id)
        if term is None:
            return ()
        return tuple(parent for parent in term.parents if parent in self._terms)

    def _children_of(self, hpo_id: str) -> tuple[str, ...]:
        return self._children.get(hpo_id, ())

    def is_ancestor_of(self, candidate: str, term: str) -> bool:
        """True when ``candidate`` subsumes ``term`` (reflexively).

        Named from the ontology's point of view, not the caller's, because
        ``is_ancestor_of(a, b)`` and ``is_descendant_of(a, b)`` are exactly the
        confusion that inverts a negative finding.
        """
        resolved = self.resolve(candidate)
        return resolved is not None and resolved in self.ancestor_closure(term)

    # -- bulk views ---------------------------------------------------------

    @property
    def term_ids(self) -> tuple[str, ...]:
        """Every non-obsolete term identifier, sorted (GP-30)."""
        return tuple(sorted(self._terms))

    @property
    def roots(self) -> tuple[str, ...]:
        """Terms with no ``is_a`` parent inside this release, sorted."""
        return tuple(
            sorted(
                hpo_id for hpo_id, term in self._terms.items() if not self._parents_of(term.hpo_id)
            )
        )

    @property
    def obsolete_ids(self) -> tuple[str, ...]:
        return tuple(sorted(self._obsolete))

    def resolve_all(self, hpo_ids: Iterable[str]) -> tuple[tuple[str, ...], tuple[str, ...]]:
        """Split identifiers into ``(resolved, unresolvable)``, both sorted and deduplicated.

        Two return values rather than one filtered list, because the terms this
        release cannot resolve are exactly the ones a caller must report as an
        information gap rather than drop.
        """
        resolved: set[str] = set()
        missing: set[str] = set()
        for raw in hpo_ids:
            primary = self.resolve(raw)
            if primary is None:
                missing.add(raw.strip().upper().replace("HP_", "HP:"))
            else:
                resolved.add(primary)
        return tuple(sorted(resolved)), tuple(sorted(missing))

    def __len__(self) -> int:
        return len(self._terms)

    def __repr__(self) -> str:
        return (
            f"HpoOntology(data_version={self._data_version!r}, terms={len(self._terms)}, "
            f"obsolete={len(self._obsolete)})"
        )


# ---------------------------------------------------------------------------
# OBO parsing
# ---------------------------------------------------------------------------


class _StanzaAccumulator:
    """Mutable scratch state for one ``[Term]`` stanza.

    A class rather than a pile of ``nonlocal`` variables so that the tag dispatch
    lives next to the fields it writes, and so the parse loop stays short enough to
    read in one screen. Instances are per-stanza and never escape :func:`_parse_obo`.
    """

    __slots__ = ("alt_ids", "hpo_id", "name", "obsolete", "parents", "replaced_by")

    def __init__(self) -> None:
        self.hpo_id: str = ""
        self.name: str = ""
        self.parents: list[str] = []
        self.alt_ids: list[str] = []
        self.obsolete: bool = False
        self.replaced_by: str | None = None

    def consume(self, line: str, *, path: Path) -> None:
        """Apply one ``tag: value`` line. Unknown tags are ignored."""
        tag, _, value = line.partition(":")
        if tag == "id":
            self.hpo_id = _term_id(_before_comment(value), path=path)
        elif tag == "name":
            self.name = value.strip()
        elif tag == "is_a":
            self.parents.append(_term_id(_before_comment(value), path=path))
        elif tag == "alt_id":
            self.alt_ids.append(_term_id(_before_comment(value), path=path))
        elif tag == "is_obsolete":
            self.obsolete = value.strip().lower() == "true"
        elif tag == "replaced_by":
            self.replaced_by = _term_id(_before_comment(value), path=path)

    def build(self) -> HpoTerm | None:
        """The accumulated term, or ``None`` for an empty or non-term stanza."""
        if not self.hpo_id:
            return None
        return HpoTerm(
            hpo_id=self.hpo_id,
            name=self.name,
            parents=tuple(sorted(set(self.parents))),
            alt_ids=tuple(sorted(set(self.alt_ids))),
            is_obsolete=self.obsolete,
            replaced_by=self.replaced_by,
        )


def _parse_obo(path: Path) -> tuple[list[HpoTerm], str]:
    """Read the subset of OBO 1.2 this package needs.

    Deliberately not a general OBO parser. It reads ``id``, ``name``, ``is_a``,
    ``alt_id``, ``is_obsolete`` and ``replaced_by`` and ignores every other tag,
    which keeps a 11 MB file to a single linear pass with no intermediate objects
    per tag. ``[Typedef]`` stanzas are skipped entirely — they describe relation
    types, and admitting them would insert non-phenotype nodes into the graph.

    Only ``is_a`` is treated as a subsumption edge. ``relationship:`` lines (none
    are present in current releases, but the tag is legal) are deliberately not
    followed: ``part_of`` and ``has_part`` are not ``is_a``, and treating them as
    subsumption would make anatomically-related terms score as semantically
    equivalent.
    """
    terms: list[HpoTerm] = []
    data_version = ""
    stanza = _StanzaAccumulator()
    in_term = False
    in_header = True

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("["):
                built = stanza.build()
                if built is not None:
                    terms.append(built)
                stanza = _StanzaAccumulator()
                in_term = line.strip() == _TERM_STANZA
                in_header = False
                continue
            if in_header:
                if line.startswith(DATA_VERSION_KEY + ":"):
                    data_version = line.split(":", 1)[1].strip()
                continue
            if in_term:
                stanza.consume(line, path=path)
    built = stanza.build()
    if built is not None:
        terms.append(built)

    if not data_version:
        msg = (
            f"{path.name} has no '{DATA_VERSION_KEY}:' header. The ontology release is "
            "mandatory provenance: a phenotype score is only interpretable against a "
            "named release, and reconstructing it from a file timestamp is a guess."
        )
        raise IngestionError(msg)
    return terms, data_version


def _before_comment(value: str) -> str:
    """Strip an OBO tag value down to its identifier.

    Two suffixes are legal on an ``is_a`` line and both appear in real releases::

        is_a: HP:0001507 ! Growth abnormality
        is_a: HP:0030991 {xref="PMID:3972225"} ! Abnormal cardiac ventricle morphology

    The trailing-modifier block ``{...}`` precedes the ``!`` comment in the OBO
    grammar, so cutting at whichever delimiter appears first is correct for both
    and for the two combined. Cutting only at ``!`` was a real bug: it left the
    modifier attached and rejected the whole release.
    """
    for delimiter in ("{", "!"):
        value = value.split(delimiter, 1)[0]
    return value


def _term_id(value: str, *, path: Path) -> str:
    """Validate one identifier token, naming the file when it is malformed."""
    token = value.strip().upper().replace("HP_", "HP:")
    if not HPO_ID_PATTERN.match(token):
        msg = (
            f"{path.name}: malformed HPO identifier {value.strip()!r}. Expected the form "
            "'HP:0001250'. Accepting it would put a node in the graph that no observation "
            "can ever match."
        )
        raise IngestionError(msg)
    return token


def ontology_provenance(ontology: HpoOntology) -> Mapping[str, str]:
    """Provenance fields for the ontology, as a plain sorted mapping.

    Returned as data rather than written anywhere: the phenotype stage does not
    own artifact writing, and the composition root decides where these land.
    """
    fields: dict[str, str] = {
        "hpo_data_version": ontology.data_version,
        "hpo_release": ontology.release,
        "hpo_term_count": str(len(ontology)),
    }
    if ontology.source_sha256 is not None:
        fields["hpo_sha256"] = ontology.source_sha256
    if ontology.source_path is not None:
        fields["hpo_source_file"] = ontology.source_path.name
    return dict(sorted(fields.items()))
