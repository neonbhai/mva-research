"""The evidence store: a DuckDB database plus a deterministic Parquet export.

Why DuckDB and not a service
----------------------------
The whole pipeline runs on one machine, offline, over patient data that may never
leave the workspace. A single embedded file gives transactional writes, real SQL
over list and JSON columns, and a zero-dependency reader — no daemon to secure, no
port to accidentally expose, no second copy of the genotypes.

Why relational and not a graph database
---------------------------------------
``graph_edges`` stores subject-predicate-object triples in an ordinary table. The
graph is small, traversal is one or two hops from a known subject, and the edges
must share a transaction, a determinism guarantee and a Parquet export with the
tabular evidence they cite. A second engine would buy traversals this workload
never performs at the price of a second integrity boundary. See schema.sql and the
corresponding decision record.

Determinism (GP-30)
-------------------
Every writer upserts on a primary key, so a re-run converges instead of appending.
Every export sorts by that primary key, pins the Parquet codec, level, version and
row-group size, and embeds no timestamp of its own. Two exports of the same data
are byte-identical, which ``tests/unit/test_evidence_store.py`` asserts with
sha256.

Privacy
-------
This database holds genotypes. Exceptions raised from here name identifiers and
counts only — never record content — because a traceback travels to every log the
process can reach.
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Final

import duckdb
import pyarrow.parquet as pq

from mva.determinism import canonical_json, hash_file, short_hash
from mva.errors import EvidenceError
from mva.models import (
    ArtifactKind,
    ArtifactProvenance,
    CandidatePair,
    Citation,
    DrugHypothesis,
    EvidenceItem,
    MechanismHypothesis,
    PhenotypeProfile,
    RunManifest,
    Sensitivity,
    VariantRecord,
    contig_sort_key,
)

_SCHEMA_PATH: Final[Path] = Path(__file__).with_name("schema.sql")

#: Every table in the store, in the order the schema declares them. This tuple is
#: the export manifest and the key set of :meth:`EvidenceStore.counts`.
TABLES: Final[tuple[str, ...]] = (
    "variants",
    "consequences",
    "frequencies",
    "clinical_assertions",
    "phenotype_observations",
    "genes",
    "candidate_pairs",
    "evidence_items",
    "mechanism_nodes",
    "mechanism_links",
    "mechanisms",
    "drugs",
    "drug_rejections",
    "citations",
    "pipeline_runs",
    "artifact_provenance",
    "graph_edges",
)

#: Primary key of each table. Doubles as the export sort key: sorting by the PK is
#: what makes a Parquet file a function of its content rather than of insertion
#: order (GP-30).
PRIMARY_KEYS: Final[Mapping[str, tuple[str, ...]]] = {
    "variants": ("variant_id",),
    "consequences": ("variant_id", "gene_symbol", "transcript_id"),
    "frequencies": ("variant_id", "source", "version", "population"),
    "clinical_assertions": ("assertion_key",),
    "phenotype_observations": ("subject_id", "hpo_id"),
    "genes": ("gene_symbol",),
    "candidate_pairs": ("pair_id",),
    "evidence_items": ("evidence_id",),
    "mechanism_nodes": ("mechanism_id", "node_id"),
    "mechanism_links": ("mechanism_id", "link_id"),
    "mechanisms": ("mechanism_id",),
    "drugs": ("drug_id",),
    "drug_rejections": ("drug_id", "reason"),
    "citations": ("citation_key",),
    "pipeline_runs": ("run_id",),
    "artifact_provenance": ("run_id", "artifact_id"),
    "graph_edges": ("edge_id",),
}

# --------------------------------------------------------------------- parquet
# Pinned so the bytes are a pure function of the data. Changing any of these
# constants changes every exported file's hash and therefore every golden
# expectation, so they are named rather than inlined.

#: zstd at a fixed level. The default level is a library choice that may move
#: between releases; pinning it keeps old exports reproducible.
PARQUET_COMPRESSION: Final[str] = "zstd"
PARQUET_COMPRESSION_LEVEL: Final[int] = 3
#: Fixed row-group size. Left to the writer, this varies with input size and
#: silently changes the file layout for identical data.
PARQUET_ROW_GROUP_SIZE: Final[int] = 122_880
PARQUET_DATA_PAGE_SIZE: Final[int] = 1 << 20
#: Format version pinned: 2.6 fixes the physical encoding of the logical types
#: used here (timestamps as int64 micros, lists as 3-level LIST).
PARQUET_VERSION: Final[str] = "2.6"


@dataclass(frozen=True)
class GraphEdge:
    """One subject-predicate-object assertion.

    An edge is a claim like any other, so it carries ``evidence_ids`` (GP-10).
    ``confidence`` is optional and is deliberately not defaulted to 1.0: an edge
    with no stated confidence means "unscored", not "certain".
    """

    subject_id: str
    subject_kind: str
    predicate: str
    object_id: str
    object_kind: str
    evidence_ids: tuple[str, ...] = ()
    confidence: float | None = None

    @property
    def edge_id(self) -> str:
        """Content-derived identity: the same triple asserted twice is one edge."""
        return "EDG-" + short_hash([self.subject_id, self.predicate, self.object_id], 16)


# ------------------------------------------------------------------- utilities


def _evidence_id_list(values: Iterable[str]) -> list[str]:
    """Normalise an evidence-ID list column: sorted and deduplicated.

    Citation lists are sets in meaning, so they are stored as sets in bytes. This
    removes the last place where caller iteration order could leak into an
    exported file.
    """
    return sorted(set(values))


def _text_list(values: Iterable[str]) -> list[str]:
    """A list column whose order is meaningful (e.g. consequence terms, ordered
    most-severe first) and is therefore preserved exactly as given."""
    return list(values)


def _naive_utc(value: datetime) -> datetime:
    """The instant as a naive UTC datetime, for the queryable TIMESTAMP column."""
    if value.tzinfo is None:
        return value
    return value.astimezone(UTC).replace(tzinfo=None)


def _iso(value: datetime | None) -> str | None:
    return None if value is None else value.isoformat()


class EvidenceStore:
    """Transactional, idempotent persistence for every domain model.

    Usage::

        with EvidenceStore(workspace / "evidence.duckdb") as store:
            store.initialise()
            store.write_evidence(ledger.items())

    Every ``write_*`` method is an upsert on the table's primary key and returns
    the number of primary records persisted. Calling one twice with the same input
    leaves the database bit-for-bit identical.
    """

    def __init__(self, db_path: Path, *, read_only: bool = False) -> None:
        self._db_path = Path(db_path)
        self._read_only = read_only
        if not read_only:
            self._db_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            self._connection: duckdb.DuckDBPyConnection | None = duckdb.connect(
                str(self._db_path), read_only=read_only
            )
        except duckdb.Error as exc:  # pragma: no cover - environment dependent
            msg = f"Could not open evidence store at {self._db_path.name!r}: {type(exc).__name__}"
            raise EvidenceError(msg) from exc
        # Single-threaded execution. At this data scale the cost is negligible and
        # it removes an entire class of ordering nondeterminism from scans that
        # feed the Parquet export.
        self._connection.execute("SET threads TO 1")

    # ------------------------------------------------------------- lifecycle

    @property
    def db_path(self) -> Path:
        return self._db_path

    @property
    def read_only(self) -> bool:
        return self._read_only

    def __enter__(self) -> EvidenceStore:
        return self

    def __exit__(
        self,
        *exc: object,
    ) -> None:
        del exc
        self.close()

    def close(self) -> None:
        """Close the connection. Idempotent."""
        if self._connection is not None:
            self._connection.close()
            self._connection = None

    def initialise(self) -> None:
        """Apply ``schema.sql``.

        Idempotent: every statement in the schema is ``CREATE ... IF NOT EXISTS``,
        so opening an existing workspace re-validates rather than recreating.
        """
        connection = self._require_connection()
        try:
            connection.execute(_SCHEMA_PATH.read_text(encoding="utf-8"))
        except duckdb.Error as exc:
            msg = f"Evidence schema could not be applied: {type(exc).__name__}"
            raise EvidenceError(msg) from exc

    def _require_connection(self) -> duckdb.DuckDBPyConnection:
        if self._connection is None:
            msg = f"Evidence store {self._db_path.name!r} is closed."
            raise EvidenceError(msg)
        return self._connection

    # ----------------------------------------------------------- low-level IO

    def _upsert(self, table: str, rows: Sequence[Mapping[str, object]]) -> int:
        """``INSERT ... ON CONFLICT DO UPDATE`` for one table.

        Rows are deduplicated on the primary key before execution (last writer
        wins) and then ordered by it. DuckDB refuses to update the same row twice
        in one command, and an undeduplicated batch is also a determinism hazard.
        """
        if not rows:
            return 0
        key_columns = PRIMARY_KEYS[table]
        columns = list(rows[0].keys())

        unique: dict[tuple[object, ...], Mapping[str, object]] = {}
        for row in rows:
            unique[tuple(row[key] for key in key_columns)] = row
        ordered = [unique[key] for key in sorted(unique, key=lambda k: tuple(str(v) for v in k))]

        placeholders = ", ".join("?" for _ in columns)
        updates = ", ".join(f"{c} = excluded.{c}" for c in columns if c not in key_columns)
        action = f"DO UPDATE SET {updates}" if updates else "DO NOTHING"
        # Table and column names come from PRIMARY_KEYS and the schema, never from
        # caller input; every value travels as a bound parameter.
        column_sql = ", ".join(columns)
        conflict_sql = ", ".join(key_columns)
        statement = (
            f"INSERT INTO {table} ({column_sql}) VALUES ({placeholders}) "  # noqa: S608
            f"ON CONFLICT ({conflict_sql}) {action}"
        )
        connection = self._require_connection()
        try:
            connection.executemany(statement, [[row[c] for c in columns] for row in ordered])
        except duckdb.Error as exc:
            msg = f"Write to {table!r} failed for {len(ordered)} row(s): {type(exc).__name__}"
            raise EvidenceError(msg) from exc
        return len(ordered)

    # ---------------------------------------------------------------- writers

    def write_variants(
        self, variants: Sequence[VariantRecord], *, run_id: str | None = None
    ) -> int:
        """Persist variant rows and their curated clinical assertions.

        Clinical assertions travel with the variant because they are part of the
        same ingested record; there is no separate stage that produces them alone.
        Returns the number of variant rows.
        """
        rows: list[Mapping[str, object]] = []
        for variant in variants:
            coordinate = variant.coordinate
            genotype = variant.genotype
            rows.append(
                {
                    "variant_id": variant.variant_id,
                    "build": coordinate.build.value,
                    "contig": coordinate.contig,
                    "contig_order": contig_sort_key(coordinate.contig),
                    "position": coordinate.position,
                    "ref": coordinate.ref,
                    "alt": coordinate.alt,
                    "zygosity": genotype.zygosity.value,
                    "genotype_string": genotype.genotype_string,
                    "phased": genotype.phased,
                    "phase_set": genotype.phase_set,
                    "depth": genotype.depth,
                    "ref_reads": genotype.ref_reads,
                    "alt_reads": genotype.alt_reads,
                    "genotype_quality": genotype.genotype_quality,
                    "allele_balance": genotype.allele_balance,
                    "filter_status": variant.filter_status.value,
                    "raw_filters": _text_list(variant.raw_filters),
                    "quality": variant.quality,
                    "qc_flags": _text_list(variant.qc_flags),
                    "normalisation_ops": _text_list(variant.normalisation_ops),
                    "source_artifact": variant.source_artifact,
                    "source_line_index": variant.source_line_index,
                    "evidence_ids": [],
                    "run_id": run_id,
                }
            )
        written = self._upsert("variants", rows)
        self.write_clinical_assertions(variants, run_id=run_id)
        return written

    def write_clinical_assertions(
        self, variants: Sequence[VariantRecord], *, run_id: str | None = None
    ) -> int:
        """Persist ClinVar-style assertions attached to the given variants."""
        rows: list[Mapping[str, object]] = []
        for variant in variants:
            for assertion in variant.clinical_assertions:
                accession = assertion.accession or ""
                key = f"{variant.variant_id}|{assertion.source}|{assertion.version}|{accession}"
                rows.append(
                    {
                        "assertion_key": key,
                        "variant_id": variant.variant_id,
                        "source": assertion.source,
                        "version": assertion.version,
                        "accession": assertion.accession,
                        "significance": assertion.significance,
                        "review_status": assertion.review_status,
                        "star_rating": assertion.star_rating,
                        "conditions": _text_list(assertion.conditions),
                        "evidence_ids": [],
                        "run_id": run_id,
                    }
                )
        return self._upsert("clinical_assertions", rows)

    def write_consequences(
        self, variants: Sequence[VariantRecord], *, run_id: str | None = None
    ) -> int:
        """Persist every transcript-scoped consequence, then rebuild ``genes``."""
        rows: list[Mapping[str, object]] = []
        for variant in variants:
            for csq in variant.consequences:
                rows.append(
                    {
                        "variant_id": variant.variant_id,
                        "gene_symbol": csq.gene_symbol,
                        "transcript_id": csq.transcript_id,
                        "gene_id": csq.gene_id,
                        "transcript_biotype": csq.transcript_biotype,
                        "is_canonical": csq.is_canonical,
                        "is_mane_select": csq.is_mane_select,
                        "consequence_terms": _text_list(csq.consequence_terms),
                        "most_severe_term": csq.most_severe_term,
                        "impact": None if csq.impact is None else csq.impact.value,
                        "hgvs_c": csq.hgvs_c,
                        "hgvs_p": csq.hgvs_p,
                        "exon": csq.exon,
                        "intron": csq.intron,
                        "protein_position": csq.protein_position,
                        "amino_acids": csq.amino_acids,
                        "splice_ai_delta_max": csq.splice_ai_delta_max,
                        "pathogenicity_scores": canonical_json(csq.pathogenicity_scores),
                        "source_tool": csq.source_tool,
                        "source_tool_version": csq.source_tool_version,
                        "evidence_ids": [],
                        "run_id": run_id,
                    }
                )
        written = self._upsert("consequences", rows)
        self._rebuild_genes()
        return written

    def _rebuild_genes(self) -> None:
        """Refresh the ``genes`` projection wholesale from ``consequences``.

        Rebuilt rather than incrementally merged: merging list columns across
        successive writes is exactly how an "idempotent" writer stops being one.
        """
        connection = self._require_connection()
        connection.execute("DELETE FROM genes")
        connection.execute(
            """
            INSERT INTO genes (
                gene_symbol, gene_id, transcript_ids, variant_ids, worst_impact, evidence_ids
            )
            SELECT
                gene_symbol,
                min(gene_id),
                list_sort(list_distinct(list(transcript_id))),
                list_sort(list_distinct(list(variant_id))),
                CASE min(
                    CASE impact
                        WHEN 'high' THEN 0
                        WHEN 'moderate' THEN 1
                        WHEN 'low' THEN 2
                        ELSE 3
                    END
                )
                    WHEN 0 THEN 'high'
                    WHEN 1 THEN 'moderate'
                    WHEN 2 THEN 'low'
                    ELSE 'modifier'
                END,
                list_sort(list_distinct(flatten(list(evidence_ids))))
            FROM consequences
            GROUP BY gene_symbol
            ORDER BY gene_symbol
            """
        )

    def write_frequencies(
        self, variants: Sequence[VariantRecord], *, run_id: str | None = None
    ) -> int:
        """Persist population allele frequencies (source + version + population)."""
        rows: list[Mapping[str, object]] = []
        for variant in variants:
            for freq in variant.population_frequencies:
                rows.append(
                    {
                        "variant_id": variant.variant_id,
                        "source": freq.source,
                        "version": freq.version,
                        "population": freq.population,
                        "allele_frequency": freq.allele_frequency,
                        "allele_count": freq.allele_count,
                        "allele_number": freq.allele_number,
                        "homozygote_count": freq.homozygote_count,
                        "filter_status": freq.filter_status,
                        "evidence_ids": [],
                        "run_id": run_id,
                    }
                )
        return self._upsert("frequencies", rows)

    def write_phenotypes(self, profile: PhenotypeProfile, *, run_id: str | None = None) -> int:
        """Persist a phenotype profile as one row per HPO observation.

        ``status`` is stored verbatim: NOT_ASSESSED and EXCLUDED are different
        facts and are never collapsed to a boolean (GP-14).
        """
        rows: list[Mapping[str, object]] = [
            {
                "subject_id": profile.subject_id,
                "hpo_id": obs.hpo_id,
                "label": obs.label,
                "status": obs.status.value,
                "onset": obs.onset.value,
                "provenance": obs.provenance,
                "extraction_confidence": obs.extraction_confidence,
                "source_excerpt_hash": obs.source_excerpt_hash,
                "notes": obs.notes,
                "source_artifact": profile.source_artifact,
                "hpo_version": profile.hpo_version,
                "evidence_ids": [],
                "run_id": run_id,
            }
            for obs in profile.observations
        ]
        return self._upsert("phenotype_observations", rows)

    def write_pairs(self, pairs: Sequence[CandidatePair], *, run_id: str | None = None) -> int:
        """Persist candidate pairs with their full score vector, not just the
        composite, so a ranking can be explained from the database alone."""
        rows: list[Mapping[str, object]] = []
        for pair in pairs:
            scores = pair.scores
            coordinate = pair.variant_a.coordinate
            rows.append(
                {
                    "pair_id": pair.pair_id,
                    "gene_symbol": pair.gene_symbol,
                    "variant_a_id": pair.variant_a.variant_id,
                    "variant_b_id": (None if pair.variant_b is None else pair.variant_b.variant_id),
                    "is_pair": pair.is_pair,
                    "inheritance_model": pair.inheritance_model.value,
                    "phase_status": pair.phase.status.value,
                    "phase_method": pair.phase.method,
                    "phase_supporting_reads": pair.phase.supporting_reads,
                    "phase_distance_bp": pair.phase.distance_bp,
                    "phase_notes": pair.phase.notes,
                    "score_analytical_validity": scores.analytical_validity,
                    "score_rarity": scores.rarity,
                    "score_molecular_consequence": scores.molecular_consequence,
                    "score_inheritance_consistency": scores.inheritance_consistency,
                    "score_phenotype_similarity": scores.phenotype_similarity,
                    "score_mechanistic_relevance": scores.mechanistic_relevance,
                    "score_evidence_quality": scores.evidence_quality,
                    "score_contradiction_penalty": scores.contradiction_penalty,
                    "composite_score": pair.composite_score,
                    "rank": pair.rank,
                    "contig_order": contig_sort_key(coordinate.contig),
                    "position": coordinate.position,
                    "supporting_evidence_ids": _evidence_id_list(pair.supporting_evidence_ids),
                    "contradicting_evidence_ids": _evidence_id_list(
                        pair.contradicting_evidence_ids
                    ),
                    "missing_evidence": canonical_json(pair.missing_evidence),
                    "blocking_question_count": len(pair.blocking_questions),
                    "recommended_next_test": pair.recommended_next_test,
                    "discriminating_experiment": pair.discriminating_experiment,
                    "rank_rationale": pair.rank_rationale,
                    "flags": _text_list(pair.flags),
                    "run_id": run_id,
                }
            )
        return self._upsert("candidate_pairs", rows)

    def write_evidence(self, items: Sequence[EvidenceItem]) -> int:
        """Persist evidence items and their citations.

        Contradictions are written exactly like support (GP-19): there is no
        filter on ``direction`` anywhere on this path.
        """
        rows: list[Mapping[str, object]] = []
        citations: list[Mapping[str, object]] = []
        for item in items:
            citation = item.citation
            rows.append(
                {
                    "evidence_id": item.evidence_id,
                    "subject_id": item.subject_id,
                    "subject_kind": item.subject_kind,
                    "claim": item.claim,
                    "category": item.category.value,
                    "direction": item.direction.value,
                    "strength": item.strength.value,
                    "evidence_type": item.evidence_type.value,
                    "tier": item.tier.value,
                    "method": item.method,
                    "tool": item.tool,
                    "tool_version": item.tool_version,
                    "limitations": item.limitations,
                    "citation_source": None if citation is None else citation.source,
                    "citation_identifier": None if citation is None else citation.identifier,
                    "citation_version": None if citation is None else citation.version,
                    "citation_url": None if citation is None else citation.url,
                    "citation_title": None if citation is None else citation.title,
                    "citation_key": None if citation is None else citation.key,
                    "timestamp": _naive_utc(item.timestamp),
                    "timestamp_iso": item.timestamp.isoformat(),
                    "run_id": item.run_id,
                    "numeric_value": item.numeric_value,
                    "payload": canonical_json(item.payload),
                }
            )
            if citation is not None:
                citations.append(
                    {
                        "citation_key": citation.key,
                        "source": citation.source,
                        "identifier": citation.identifier,
                        "version": citation.version,
                        "url": citation.url,
                        "title": citation.title,
                    }
                )
        written = self._upsert("evidence_items", rows)
        self._upsert("citations", citations)
        return written

    def write_mechanism(self, mechanism: MechanismHypothesis, *, run_id: str | None = None) -> int:
        """Persist a mechanism hypothesis with its nodes and links.

        Returns 1: the nodes and links are parts of the same record, written in
        the same call so a mechanism can never be half-persisted.
        """
        self._upsert(
            "mechanisms",
            [
                {
                    "mechanism_id": mechanism.mechanism_id,
                    "gene_symbol": mechanism.gene_symbol,
                    "pair_id": mechanism.pair_id,
                    "summary": mechanism.summary,
                    "disease_direction": mechanism.disease_direction.value,
                    "therapeutic_target_node_id": mechanism.therapeutic_target_node_id,
                    "required_correction": mechanism.required_correction.value,
                    "node_count": len(mechanism.nodes),
                    "link_count": len(mechanism.links),
                    "inferred_link_count": len(mechanism.inferred_links),
                    "is_fully_demonstrated": mechanism.is_fully_demonstrated,
                    "supporting_evidence_ids": _evidence_id_list(mechanism.supporting_evidence_ids),
                    "contradicting_evidence_ids": _evidence_id_list(
                        mechanism.contradicting_evidence_ids
                    ),
                    "uncertainties": _text_list(mechanism.uncertainties),
                    "discriminating_experiments": canonical_json(
                        mechanism.discriminating_experiments
                    ),
                    "developmental_window_caveat": mechanism.developmental_window_caveat,
                    "run_id": run_id,
                }
            ],
        )
        self._upsert(
            "mechanism_nodes",
            [
                {
                    "mechanism_id": mechanism.mechanism_id,
                    "node_id": node.node_id,
                    "kind": node.kind.value,
                    "label": node.label,
                    "identifier": node.identifier,
                    "state_in_patient": node.state_in_patient.value,
                    "description": node.description,
                }
                for node in mechanism.nodes
            ],
        )
        self._upsert(
            "mechanism_links",
            [
                {
                    "mechanism_id": mechanism.mechanism_id,
                    "link_id": link.link_id,
                    "source_node_id": link.source_node_id,
                    "target_node_id": link.target_node_id,
                    "relation": link.relation,
                    "direction": link.direction.value,
                    "tier": link.tier.value,
                    "strength": link.strength.value,
                    "is_directly_demonstrated": link.is_directly_demonstrated,
                    "uncertainty": link.uncertainty,
                    "evidence_ids": _evidence_id_list(link.evidence_ids),
                    "contradicting_evidence_ids": _evidence_id_list(
                        link.contradicting_evidence_ids
                    ),
                }
                for link in mechanism.links
            ],
        )
        return 1

    def write_drugs(
        self,
        drugs: Sequence[DrugHypothesis],
        *,
        mechanism_id: str | None = None,
        run_id: str | None = None,
    ) -> int:
        """Persist drug hypotheses, accepted and rejected alike.

        Rejections are additionally exploded into ``drug_rejections``, one row per
        reason, so "what did we consider and why did we drop it?" is answerable
        without parsing prose (GP-19).
        """
        rows: list[Mapping[str, object]] = []
        rejections: list[Mapping[str, object]] = []
        for drug in drugs:
            pediatric = drug.pediatric_evidence
            pk = drug.pharmacokinetics
            rows.append(
                {
                    "drug_id": drug.drug_id,
                    "name": drug.name,
                    "approved_name": drug.approved_name,
                    "approval_status": drug.approval_status.value,
                    "intervention_class": drug.intervention_class.value,
                    "target": drug.target,
                    "target_node_id": drug.target_node_id,
                    "mechanism_id": mechanism_id,
                    "mechanism_of_action": drug.mechanism_of_action,
                    "required_direction": drug.required_direction.value,
                    "observed_direction": drug.observed_direction.value,
                    "directions_agree": drug.directions_agree,
                    "is_direct_evidence": drug.is_direct_evidence,
                    "strongest_evidence_type": drug.strongest_evidence_type.value,
                    "has_in_vivo_evidence": drug.has_in_vivo_evidence,
                    "is_repurposable": drug.is_repurposable,
                    "pediatric_has_exposure": pediatric.has_pediatric_exposure,
                    "pediatric_youngest_age_studied": pediatric.youngest_age_studied,
                    "pediatric_indication": pediatric.indication,
                    "pediatric_tolerability_summary": pediatric.tolerability_summary,
                    "pediatric_caveat": pediatric.caveat,
                    "pediatric_evidence_ids": _evidence_id_list(pediatric.evidence_ids),
                    "pk_route": pk.route,
                    "pk_cns_penetrant": pk.cns_penetrant,
                    "pk_achievable_plasma_um": pk.achievable_plasma_concentration_um,
                    "pk_required_effective_um": pk.required_effective_concentration_um,
                    "pk_concentration_achievable": pk.concentration_achievable,
                    "pk_half_life_hours": pk.half_life_hours,
                    "pk_notes": pk.notes,
                    "safety_concerns": canonical_json(drug.safety_concerns),
                    "disqualifying_safety_count": len(drug.disqualifying_safety),
                    "worsens_chromosomal_instability": drug.worsens_chromosomal_instability,
                    "proposed_validation_experiment": drug.proposed_validation_experiment,
                    "evidence_ids": _evidence_id_list(drug.evidence_ids),
                    "contradicting_evidence_ids": _evidence_id_list(
                        drug.contradicting_evidence_ids
                    ),
                    "score": drug.score,
                    "rank": drug.rank,
                    "rejected": drug.rejected,
                    "rejection_rationale": drug.rejection_rationale,
                    "run_id": run_id,
                }
            )
            rejections.extend(
                {
                    "drug_id": drug.drug_id,
                    "reason": reason.value,
                    "drug_name": drug.name,
                    "rationale": drug.rejection_rationale,
                    "required_direction": drug.required_direction.value,
                    "observed_direction": drug.observed_direction.value,
                    "run_id": run_id,
                }
                for reason in drug.rejection_reasons
            )
        written = self._upsert("drugs", rows)
        self._upsert("drug_rejections", rejections)
        return written

    def write_run_manifest(self, manifest: RunManifest) -> int:
        """Persist a run manifest and every artifact it registered (GP-31)."""
        self._upsert(
            "pipeline_runs",
            [
                {
                    "run_id": manifest.run_id,
                    "case_id": manifest.case_id,
                    "genome_build": manifest.genome_build,
                    "started_at": _naive_utc(manifest.started_at),
                    "started_at_iso": manifest.started_at.isoformat(),
                    "completed_at": (
                        None if manifest.completed_at is None else _naive_utc(manifest.completed_at)
                    ),
                    "completed_at_iso": _iso(manifest.completed_at),
                    "config_hash": manifest.config_hash,
                    "config_snapshot": canonical_json(manifest.config_snapshot),
                    "git_commit": manifest.git_commit,
                    "git_dirty": manifest.git_dirty,
                    "is_reproducible": manifest.is_reproducible,
                    "inputs": canonical_json(manifest.inputs),
                    "commands": canonical_json(manifest.commands),
                    "tool_versions": canonical_json(manifest.tool_versions),
                    "reference_versions": canonical_json(manifest.reference_versions),
                    "python_version": manifest.python_version,
                    "platform": manifest.platform,
                    "network_profile": manifest.network_profile,
                    "synthetic": manifest.synthetic,
                    "warnings": _text_list(manifest.warnings),
                }
            ],
        )
        self._upsert(
            "artifact_provenance",
            [
                {
                    "run_id": manifest.run_id,
                    "artifact_id": artifact.artifact_id,
                    "kind": artifact.kind.value,
                    "relative_path": artifact.relative_path,
                    "sensitivity": artifact.sensitivity.value,
                    "is_exportable": artifact.is_exportable,
                    "content_hash": artifact.content_hash,
                    "size_bytes": artifact.size_bytes,
                    "produced_by_stage": artifact.produced_by_stage,
                    "upstream_artifact_ids": _text_list(artifact.upstream_artifact_ids),
                    "tool_versions": canonical_json(artifact.tool_versions),
                    "created_at": _naive_utc(artifact.created_at),
                    "created_at_iso": artifact.created_at.isoformat(),
                    "row_count": artifact.row_count,
                    "notes": artifact.notes,
                }
                for artifact in manifest.artifacts
            ],
        )
        return 1

    def write_edges(self, edges: Sequence[GraphEdge], *, run_id: str | None = None) -> int:
        """Persist subject-predicate-object edges into the relational graph table."""
        rows: list[Mapping[str, object]] = [
            {
                "edge_id": edge.edge_id,
                "subject_id": edge.subject_id,
                "subject_kind": edge.subject_kind,
                "predicate": edge.predicate,
                "object_id": edge.object_id,
                "object_kind": edge.object_kind,
                "evidence_ids": _evidence_id_list(edge.evidence_ids),
                "confidence": edge.confidence,
                "run_id": run_id,
            }
            for edge in edges
        ]
        return self._upsert("graph_edges", rows)

    # ---------------------------------------------------------------- readers

    def evidence_for(self, subject_id: str) -> tuple[EvidenceItem, ...]:
        """Every evidence item about ``subject_id``, in evidence-ID order.

        Support and contradiction are returned together, on purpose: a caller that
        has to opt in to seeing the counter-evidence will forget to.
        """
        rows = self.query(
            "SELECT * FROM evidence_items WHERE subject_id = ? ORDER BY evidence_id",
            (subject_id,),
        )
        return tuple(_row_to_evidence(row) for row in rows)

    def contradictions_for(self, subject_id: str) -> tuple[EvidenceItem, ...]:
        """Only the items that argue against ``subject_id`` (GP-19)."""
        rows = self.query(
            "SELECT * FROM evidence_items WHERE subject_id = ? AND direction = 'contradicts' "
            "ORDER BY evidence_id",
            (subject_id,),
        )
        return tuple(_row_to_evidence(row) for row in rows)

    def ranked_pairs(self, limit: int | None = None) -> list[dict[str, object]]:
        """Candidate pairs ordered by composite score, highest first.

        The tiebreak is genomic position then pair ID, giving a total order: equal
        scores must not be resolved by storage order, or the ranking stops being
        reproducible (GP-30).
        """
        if limit is not None and limit < 0:
            msg = f"ranked_pairs limit must be non-negative, got {limit}."
            raise EvidenceError(msg)
        statement = (
            "SELECT pair_id, gene_symbol, variant_a_id, variant_b_id, inheritance_model, "
            "phase_status, composite_score, rank, blocking_question_count, "
            "supporting_evidence_ids, contradicting_evidence_ids, flags "
            "FROM candidate_pairs "
            "ORDER BY composite_score DESC, contig_order ASC, position ASC, pair_id ASC"
        )
        if limit is not None:
            return self.query(f"{statement} LIMIT {int(limit)}")
        return self.query(statement)

    def query(self, sql: str, params: Sequence[object] = ()) -> list[dict[str, object]]:
        """Run an arbitrary read query and return rows as dicts."""
        connection = self._require_connection()
        try:
            cursor = connection.execute(sql, list(params)) if params else connection.execute(sql)
            description = cursor.description
            if not description:
                return []
            columns = [str(column[0]) for column in description]
            return [dict(zip(columns, row, strict=True)) for row in cursor.fetchall()]
        except duckdb.Error as exc:
            msg = f"Evidence query failed ({len(params)} parameter(s)): {type(exc).__name__}"
            raise EvidenceError(msg) from exc

    def counts(self) -> dict[str, int]:
        """Row count of every table, in schema order."""
        connection = self._require_connection()
        result: dict[str, int] = {}
        for table in TABLES:
            row = connection.execute(f"SELECT count(*) FROM {table}").fetchone()  # noqa: S608 - name from TABLES
            result[table] = 0 if row is None else int(row[0])
        return result

    # ---------------------------------------------------------------- parquet

    def export_parquet(self, out_dir: Path, tables: Sequence[str] | None = None) -> dict[str, Path]:
        """Export tables to Parquet, byte-identically on repeat runs.

        The guarantee rests on four things and nothing else: rows are sorted by
        primary key, the codec and its level are pinned, the row-group and page
        sizes are fixed, and the writer embeds no clock reading. Empty tables are
        exported too, so the file set is a function of the schema rather than of
        which stages happened to run.
        """
        selected = TABLES if tables is None else tuple(tables)
        unknown = [name for name in selected if name not in PRIMARY_KEYS]
        if unknown:
            msg = f"Unknown table(s) requested for export: {sorted(unknown)}"
            raise EvidenceError(msg)

        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        connection = self._require_connection()
        written: dict[str, Path] = {}
        for table in selected:
            order_by = ", ".join(PRIMARY_KEYS[table])
            result = connection.execute(
                f"SELECT * FROM {table} ORDER BY {order_by}"  # noqa: S608 - names from PRIMARY_KEYS
            )
            # DuckDB 1.5 renamed ``fetch_arrow_table`` to ``to_arrow_table`` and
            # deprecated the old spelling. This project turns DeprecationWarnings
            # from mva.* into errors, so prefer the new name and fall back only on
            # the older DuckDB releases the dependency range still allows.
            to_arrow_table = getattr(result, "to_arrow_table", None)
            arrow_table = (
                to_arrow_table() if to_arrow_table is not None else result.fetch_arrow_table()
            )
            path = out_dir / f"{table}.parquet"
            pq.write_table(
                arrow_table,
                path,
                compression=PARQUET_COMPRESSION,
                compression_level=PARQUET_COMPRESSION_LEVEL,
                row_group_size=PARQUET_ROW_GROUP_SIZE,
                data_page_size=PARQUET_DATA_PAGE_SIZE,
                version=PARQUET_VERSION,
                use_dictionary=True,
                write_statistics=True,
                store_schema=True,
                coerce_timestamps="us",
                allow_truncated_timestamps=False,
            )
            written[table] = path
        return written

    def export_provenance(
        self,
        exported: Mapping[str, Path],
        *,
        created_at: datetime,
        run_id: str,
        workspace_root: Path | None = None,
        produced_by_stage: str = "evidence.export_parquet",
    ) -> tuple[ArtifactProvenance, ...]:
        """Describe an export as :class:`ArtifactProvenance` records (GP-31).

        Returned rather than written: inserting provenance for an export into the
        same database would change the next export's bytes and break GP-30. The
        composition root decides where these records land.

        ``created_at`` comes from the injected clock, never from a wall-clock read.
        Sensitivity is SENSITIVE unconditionally — this database holds genotypes,
        and promotion to PUBLIC is an explicit gate, never a default.
        """
        counts = self.counts()
        records: list[ArtifactProvenance] = []
        for table in sorted(exported):
            path = exported[table]
            relative = (
                path.name if workspace_root is None else str(path.relative_to(workspace_root))
            )
            records.append(
                ArtifactProvenance(
                    artifact_id=f"{run_id}:evidence:{table}",
                    kind=ArtifactKind.EVIDENCE_DB,
                    relative_path=str(relative),
                    sensitivity=Sensitivity.SENSITIVE,
                    content_hash=hash_file(path),
                    size_bytes=path.stat().st_size,
                    produced_by_stage=produced_by_stage,
                    created_at=created_at,
                    row_count=counts.get(table, 0),
                    notes=f"Parquet export of table {table!r}.",
                )
            )
        return tuple(records)


def _row_to_evidence(row: Mapping[str, Any]) -> EvidenceItem:
    """Rebuild an :class:`EvidenceItem` from a database row.

    The timestamp is reconstructed from the stored ISO text rather than from the
    TIMESTAMP column, so an aware datetime survives the round trip exactly.
    Pydantic re-validates on construction, which means a row that violates a model
    invariant (a database assertion with no versioned citation, say) fails loudly
    on read instead of quietly re-entering the pipeline.
    """
    citation: Citation | None = None
    if row["citation_source"] is not None:
        citation = Citation(
            source=str(row["citation_source"]),
            identifier=str(row["citation_identifier"]),
            version=row["citation_version"],
            url=row["citation_url"],
            title=row["citation_title"],
        )
    payload_text = row["payload"]
    payload: dict[str, str | int | float | bool | None] = (
        json.loads(payload_text) if payload_text else {}
    )
    return EvidenceItem(
        evidence_id=str(row["evidence_id"]),
        subject_id=str(row["subject_id"]),
        subject_kind=str(row["subject_kind"]),
        claim=str(row["claim"]),
        category=row["category"],
        direction=row["direction"],
        strength=row["strength"],
        evidence_type=row["evidence_type"],
        tier=row["tier"],
        citation=citation,
        method=str(row["method"]),
        tool=str(row["tool"]),
        tool_version=str(row["tool_version"]),
        limitations=str(row["limitations"]),
        timestamp=datetime.fromisoformat(str(row["timestamp_iso"])),
        run_id=row["run_id"],
        numeric_value=row["numeric_value"],
        payload=payload,
    )


__all__ = [
    "PARQUET_COMPRESSION",
    "PARQUET_COMPRESSION_LEVEL",
    "PARQUET_ROW_GROUP_SIZE",
    "PARQUET_VERSION",
    "PRIMARY_KEYS",
    "TABLES",
    "EvidenceStore",
    "GraphEdge",
]
