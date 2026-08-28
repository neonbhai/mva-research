"""Regenerate the committed HPO test subgraph from a real HPO release.

Run it with::

    MVA_RESOURCES=~/Contri/bio-hackathon/mva-resources \
        uv run python tests/fixtures/hpo/build_fixture.py

Each of ``hp.obo``, ``phenotype.hpoa`` and ``genes_to_phenotype.txt`` (from
https://hpo.jax.org/data/annotations) is located through
``knowledge/manifests/resources.yaml`` and its sha256 is checked against the
digest that file pins, so this script cannot build a "real HPO slice" out of
whatever happens to be in a directory someone typed. That is the generator half
of ADR 0012 condition 2; ``fixture_provenance.yaml`` in this directory is the
half the privacy audit reads.

The release is 67 MB and is **not** committed; this script cuts a ~600-term slice
of the real graph, with real annotation rows, that is small enough to live in the
repository and still exercises everything the scorer depends on.

Why a real slice rather than a hand-written toy ontology
--------------------------------------------------------
A synthetic graph tests the code against the author's mental model of HPO, which
is exactly where the bugs are. The three structural facts that break naive
implementations are all preserved here, and each is asserted by
``tests/unit/test_hpo_ontology.py``:

* **Multi-parent terms.** ``HP:0001305`` (Dandy-Walker malformation) has four real
  ``is_a`` parents, and roughly 30 other terms in the slice have more than one.
* **Unequal path lengths to the root.** ``HP:0007370`` (Aplasia/Hypoplasia of the
  corpus callosum) reaches ``HP:0000001`` in as few as 7 and as many as 12 steps.
  Any measure that uses depth as a proxy for specificity gets a different answer
  depending on which path it walked; this is why information content comes from
  the annotation corpus instead.
* **Re-convergence.** The slice is closed upward, so the same ancestor is reached
  by many distinct paths and a closure that does not memoise blows up.

Selection rule (deterministic, and the whole rule)
--------------------------------------------------
1. Start from :data:`SEED_TERMS` — the terms the synthetic case uses, plus two
   broad terms so downward exclusion propagation has something to propagate over.
2. Close **upward**: add every ancestor, so every path to the root is intact.
3. Add one **downward** layer: every direct child of every term so far, so
   ``descendant_closure`` has real structure to walk.
4. Close upward again, so no ``is_a`` edge points outside the slice.

Annotation rows are the real rows from the release, filtered to terms in the
slice and capped at :data:`MAX_DISEASES` / :data:`MAX_GENES` entities. Entities
annotating a seed term are taken first, then the rest, both in sorted identifier
order — see :func:`_slice_rows`. Nothing is invented, and nothing is sampled
randomly.
"""

from __future__ import annotations

import argparse
import os
from collections.abc import Iterable, Sequence
from pathlib import Path

from mva.config import find_repo_root
from mva.phenotype.ontology import HpoOntology, HpoTerm
from mva.privacy.audit import pinned_source


def _default_resource_root() -> Path | None:
    """``$MVA_RESOURCES``, expanded, if it is set.

    Where the acquisition tool puts its downloads is that tool's business
    (``tools.acquire.fetch.resolve_resource_root``); this only reads the same
    environment variable rather than growing a second opinion about the default.
    """
    raw = os.environ.get("MVA_RESOURCES")
    return Path(raw).expanduser() if raw else None


#: The manifest resource each input comes from. Every name here must match the
#: ``resource:`` recorded for the corresponding fixture in
#: ``fixture_provenance.yaml``, which is what the privacy audit reads.
ONTOLOGY_RESOURCE = "hpo_ontology"
HPOA_RESOURCE = "hpo_disease_annotations"
G2P_RESOURCE = "hpo_genes_to_phenotype"

#: Terms the slice is grown from.
#:
#: The first eight are the synthetic case's own terms, so the fixture can score the
#: same profile the rest of the suite uses. ``HP:0000707`` and ``HP:0000478`` are
#: deliberately broad: excluding a whole organ system is the case that exercises
#: downward propagation, and without a broad term in the slice there is nothing to
#: exclude.
SEED_TERMS: tuple[str, ...] = (
    "HP:0000252",  # Microcephaly - two is_a parents
    "HP:0000478",  # Abnormality of the eye - broad, for downward exclusion
    "HP:0000518",  # Cataract
    "HP:0000707",  # Abnormality of the nervous system - broad
    "HP:0000822",  # Hypertension
    "HP:0001250",  # Seizure
    "HP:0001305",  # Dandy-Walker malformation - four is_a parents
    "HP:0001511",  # Intrauterine growth retardation
    "HP:0002667",  # Nephroblastoma
    "HP:0200024",  # Premature chromatid separation
)

#: Entity caps. Large enough that information content is non-degenerate, small
#: enough that the committed files stay well under a megabyte.
MAX_DISEASES: int = 300
MAX_GENES: int = 200

_OBO_NAME = "hp_subset.obo"
_HPOA_NAME = "phenotype_subset.hpoa"
_G2P_NAME = "genes_to_phenotype_subset.txt"

_REGENERATION_NOTE = (
    "Generated by tests/fixtures/hpo/build_fixture.py from a real HPO release. "
    "Regenerate with: MVA_RESOURCES=<resource-root> uv run python "
    "tests/fixtures/hpo/build_fixture.py"
)


def select_terms(ontology: HpoOntology, seeds: Sequence[str]) -> tuple[str, ...]:
    """Apply the four-step selection rule. Sorted, so the output is byte-stable."""
    upward: set[str] = set()
    for seed in sorted(seeds):
        upward.update(ontology.ancestor_closure(seed))

    with_children: set[str] = set(upward)
    for term in sorted(upward):
        with_children.update(ontology.children(term))

    closed: set[str] = set(with_children)
    for term in sorted(with_children):
        closed.update(ontology.ancestor_closure(term))
    return tuple(sorted(closed))


def render_obo(ontology: HpoOntology, terms: Sequence[str]) -> str:
    """Emit the slice as OBO, keeping only ``is_a`` edges that stay inside it."""
    included = frozenset(terms)
    lines: list[str] = [
        f"! {_REGENERATION_NOTE}",
        f"! {len(terms)} terms, sliced from the release named in data-version below.",
        "! Seed terms: " + ", ".join(sorted(SEED_TERMS)),
        "format-version: 1.2",
        f"data-version: {ontology.data_version}",
        "default-namespace: human_phenotype",
        "ontology: hp.obo",
        "",
    ]
    for term_id in terms:
        term: HpoTerm | None = ontology.term(term_id)
        if term is None:  # pragma: no cover - selection only yields known terms
            continue
        lines.append("[Term]")
        lines.append(f"id: {term.hpo_id}")
        lines.append(f"name: {term.name}")
        for alt in term.alt_ids:
            lines.append(f"alt_id: {alt}")
        for parent in term.parents:
            if parent in included:
                lines.append(f"is_a: {parent} ! {ontology.label(parent) or ''}".rstrip())
        lines.append("")
    return "\n".join(lines) + "\n"


def _read_rows(path: Path) -> tuple[list[str], list[str], list[list[str]]]:
    """Split a ``#``-commented TSV into (comments, header, data rows)."""
    comments: list[str] = []
    header: list[str] = []
    rows: list[list[str]] = []
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        if line.startswith("#"):
            comments.append(line)
            continue
        if not line.strip():
            continue
        fields = line.split("\t")
        if not header:
            header = fields
            continue
        rows.append(fields)
    return comments, header, rows


def _slice_rows(
    rows: Iterable[Sequence[str]],
    *,
    entity_column: int,
    term_column: int,
    included: frozenset[str],
    max_entities: int,
) -> list[Sequence[str]]:
    """Keep rows whose term is in the slice, for at most ``max_entities`` entities.

    Entities that annotate a :data:`SEED_TERMS` member are taken first, then the
    rest; both groups in sorted identifier order. The seed-first rule is what keeps
    the rare seed terms — the ones the scoring tests actually compare — inside the
    information-content corpus. A plain alphabetical cut dropped every disease
    annotated with ``HP:0200024``, leaving that term with no IC and the fixture
    unable to exercise the measure it exists to test.

    Ordering by identifier rather than by file position means the fixture does not
    change when the release reorders its rows.
    """
    by_entity: dict[str, list[Sequence[str]]] = {}
    seeded: set[str] = set()
    seeds = frozenset(SEED_TERMS)
    for row in rows:
        if len(row) <= max(entity_column, term_column):
            continue
        term = row[term_column].strip()
        if term not in included:
            continue
        entity = row[entity_column].strip()
        by_entity.setdefault(entity, []).append(row)
        if term in seeds:
            seeded.add(entity)
    ordered = sorted(seeded) + sorted(set(by_entity) - seeded)
    kept: list[Sequence[str]] = []
    for entity in ordered[:max_entities]:
        kept.extend(by_entity[entity])
    return kept


def _render_tsv(
    comments: Sequence[str], header: Sequence[str], rows: Sequence[Sequence[str]]
) -> str:
    lines = [f"# {_REGENERATION_NOTE}", *comments, "\t".join(header)]
    lines.extend("\t".join(row) for row in rows)
    return "\n".join(lines) + "\n"


def build(out: Path, *, resource_root: Path | None = None) -> None:
    """Write all three fixture files from the manifest-pinned releases.

    Each input is hash-verified against ``knowledge/manifests/resources.yaml``
    before a byte of it is read, so a directory of look-alike files cannot be
    passed off as an HPO release (ADR 0012 condition 2).
    """
    repo_root = find_repo_root()
    obo_source = pinned_source(repo_root, ONTOLOGY_RESOURCE, resource_root=resource_root)
    hpoa_source = pinned_source(repo_root, HPOA_RESOURCE, resource_root=resource_root)
    g2p_source = pinned_source(repo_root, G2P_RESOURCE, resource_root=resource_root)

    ontology = HpoOntology.from_obo(obo_source)
    terms = select_terms(ontology, SEED_TERMS)
    included = frozenset(terms)
    out.mkdir(parents=True, exist_ok=True)

    (out / _OBO_NAME).write_text(render_obo(ontology, terms), encoding="utf-8")

    comments, header, rows = _read_rows(hpoa_source)
    hpoa_rows = _slice_rows(
        rows,
        entity_column=header.index("database_id"),
        term_column=header.index("hpo_id"),
        included=included,
        max_entities=MAX_DISEASES,
    )
    (out / _HPOA_NAME).write_text(_render_tsv(comments, header, hpoa_rows), encoding="utf-8")

    comments, header, rows = _read_rows(g2p_source)
    gene_rows = _slice_rows(
        rows,
        entity_column=header.index("gene_symbol"),
        term_column=header.index("hpo_id"),
        included=included,
        max_entities=MAX_GENES,
    )
    (out / _G2P_NAME).write_text(_render_tsv(comments, header, gene_rows), encoding="utf-8")

    print(f"{len(terms)} terms -> {out / _OBO_NAME}")
    print(f"{len(hpoa_rows)} disease rows -> {out / _HPOA_NAME}")
    print(f"{len(gene_rows)} gene rows -> {out / _G2P_NAME}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--resource-root",
        type=Path,
        default=_default_resource_root(),
        help="External resource root holding the acquired releases. Defaults to $MVA_RESOURCES.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory to write the fixture into (default: this directory).",
    )
    args = parser.parse_args()
    build(args.out, resource_root=args.resource_root)


if __name__ == "__main__":
    main()
