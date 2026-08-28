-- ---------------------------------------------------------------------------
-- mva evidence store — DuckDB DDL
-- ---------------------------------------------------------------------------
-- This file is the single source of truth for the persistence layer. It is
-- applied verbatim by EvidenceStore.initialise() and is written to be *idempotent*
-- (every statement is CREATE ... IF NOT EXISTS), so re-opening a workspace is
-- always safe.
--
-- Cross-cutting design choices, stated once here rather than repeated per table:
--
--  * EXPLICIT TYPES EVERYWHERE. No inferred columns. A silently widened type is
--    how a position becomes a float and a coordinate stops joining.
--
--  * EVERY TABLE HAS A PRIMARY KEY, and every writer upserts on it. Re-running a
--    stage must converge, not accumulate: a pipeline that appends duplicates on
--    retry cannot satisfy GP-30 (byte-identical repeat runs).
--
--  * `evidence_ids VARCHAR[]` ON EVERY CLAIM-BEARING TABLE. DuckDB's native LIST
--    type is used in preference to a join table because these lists are short,
--    always read whole, and never joined against on their own. GP-10 ("no claim
--    without evidence") is only enforceable if the citation travels *with* the
--    claim rather than in a table someone can forget to query. Lists are stored
--    sorted and deduplicated by the writer so the bytes are stable (GP-30).
--
--  * TEXT-ENCODED ENUMS. Every domain enum is a Python StrEnum; storing the
--    string keeps the database readable by a human with a SQL prompt, keeps
--    Parquet exports self-describing, and avoids a schema migration every time an
--    enum gains a member. Validation lives in the Pydantic models, which are the
--    only writers.
--
--  * TIMESTAMPS ARE STORED TWICE, deliberately: `*_at` as DuckDB TIMESTAMP (UTC,
--    naive) so SQL predicates and range scans work, and `*_at_iso` as the exact
--    ISO-8601 text of the original aware datetime so a read reconstructs the
--    model field byte-for-byte. DuckDB's TIMESTAMPTZ fetch path requires `pytz`,
--    which is not a dependency of this project, so TIMESTAMPTZ is avoided.
--
--  * SORT KEYS ARE MATERIALISED (e.g. `contig_order`). Karyotype order is not
--    lexicographic; storing the ordinal means every reader gets the same total
--    order without re-deriving it, which is what makes Parquet exports stable.
--
--  * JSON COLUMNS hold structured detail that is read as a whole and never
--    filtered on (payloads, score snapshots, config). They are always written as
--    canonical JSON (sorted keys, no insignificant whitespace) — see
--    mva.determinism.canonical_json — so identical content yields identical bytes.
-- ---------------------------------------------------------------------------


-- ------------------------------------------------------------------ variants

-- One row per called variant. `variant_id` is the build-qualified canonical ID
-- (GenomicCoordinate.variant_id) and is the join key for the whole database.
CREATE TABLE IF NOT EXISTS variants (
    variant_id          VARCHAR PRIMARY KEY,
    build               VARCHAR   NOT NULL,
    contig              VARCHAR   NOT NULL,
    -- Karyotype ordinal (chr1..chr22, X, Y, M). Materialised so that ORDER BY is
    -- a total order over integers rather than a string sort that puts chr10 first.
    contig_order        INTEGER   NOT NULL,
    position            BIGINT    NOT NULL,
    ref                 VARCHAR   NOT NULL,
    alt                 VARCHAR   NOT NULL,
    zygosity            VARCHAR   NOT NULL,
    genotype_string     VARCHAR   NOT NULL,
    phased              BOOLEAN   NOT NULL,
    phase_set           BIGINT,
    depth               INTEGER,
    ref_reads           INTEGER,
    alt_reads           INTEGER,
    genotype_quality    INTEGER,
    -- Computed on the model, persisted here so mosaicism screens are a SQL query.
    -- NULL means "unknowable", never 0.5 (see Genotype.allele_balance).
    allele_balance      DOUBLE,
    filter_status       VARCHAR   NOT NULL,
    raw_filters         VARCHAR[] NOT NULL,
    quality             DOUBLE,
    -- QC flags down-rank; they never delete (GP-13). Kept as a list so the
    -- flagged record survives with its reasons attached.
    qc_flags            VARCHAR[] NOT NULL,
    normalisation_ops   VARCHAR[] NOT NULL,
    source_artifact     VARCHAR   NOT NULL,
    source_line_index   BIGINT,
    evidence_ids        VARCHAR[] NOT NULL,
    run_id              VARCHAR
);

CREATE INDEX IF NOT EXISTS idx_variants_locus ON variants (contig_order, position);


-- Transcript-scoped consequence predictions. Deliberately NOT collapsed to the
-- canonical transcript: a variant can be benign on MANE-Select and splice
-- disrupting on the tissue-relevant isoform, and the canonical-only shortcut is a
-- known way to lose the real finding.
CREATE TABLE IF NOT EXISTS consequences (
    variant_id            VARCHAR   NOT NULL,
    gene_symbol           VARCHAR   NOT NULL,
    transcript_id         VARCHAR   NOT NULL,
    gene_id               VARCHAR,
    transcript_biotype    VARCHAR   NOT NULL,
    is_canonical          BOOLEAN   NOT NULL,
    is_mane_select        BOOLEAN   NOT NULL,
    consequence_terms     VARCHAR[] NOT NULL,
    most_severe_term      VARCHAR   NOT NULL,
    -- NULLABLE ON PURPOSE (ADR 0016). NULL means NOT ASSESSED: the source located
    -- the gene -- a MANE interval join answers "which gene is this variant in?" --
    -- but computed no molecular consequence. NULL is emphatically not 'modifier':
    -- MODIFIER is a positive prediction of negligible effect and
    -- prioritization.filters.BENIGN_IMPACTS treats it as one, so writing 'modifier'
    -- here would file a variant nobody assessed as predicted-benign (GP-14).
    impact                VARCHAR,
    hgvs_c                VARCHAR,
    hgvs_p                VARCHAR,
    exon                  VARCHAR,
    intron                VARCHAR,
    protein_position      INTEGER,
    amino_acids           VARCHAR,
    splice_ai_delta_max   DOUBLE,
    -- {'CADD_phred': 28.1, 'REVEL': 0.81}. JSON rather than columns because the
    -- set of scorers is configuration, not schema.
    pathogenicity_scores  JSON      NOT NULL,
    source_tool           VARCHAR   NOT NULL,
    source_tool_version   VARCHAR   NOT NULL,
    evidence_ids          VARCHAR[] NOT NULL,
    run_id                VARCHAR,
    PRIMARY KEY (variant_id, gene_symbol, transcript_id)
);

CREATE INDEX IF NOT EXISTS idx_consequences_gene ON consequences (gene_symbol);


-- Population allele frequencies. The primary key includes source, version AND
-- population because "AF = 0.001" without all three is not a fact (GP-18).
CREATE TABLE IF NOT EXISTS frequencies (
    variant_id        VARCHAR   NOT NULL,
    source            VARCHAR   NOT NULL,
    version           VARCHAR   NOT NULL,
    population        VARCHAR   NOT NULL,
    allele_frequency  DOUBLE    NOT NULL,
    allele_count      BIGINT,
    allele_number     BIGINT,
    homozygote_count  BIGINT,
    filter_status     VARCHAR,
    evidence_ids      VARCHAR[] NOT NULL,
    run_id            VARCHAR,
    PRIMARY KEY (variant_id, source, version, population)
);


-- Curated clinical-significance assertions (ClinVar-style).
-- `assertion_key` is a synthetic, deterministic primary key: `accession` is
-- legitimately NULL for some submitters, and DuckDB primary keys are NOT NULL.
CREATE TABLE IF NOT EXISTS clinical_assertions (
    assertion_key   VARCHAR PRIMARY KEY,
    variant_id      VARCHAR   NOT NULL,
    source          VARCHAR   NOT NULL,
    version         VARCHAR   NOT NULL,
    accession       VARCHAR,
    significance    VARCHAR   NOT NULL,
    review_status   VARCHAR,
    star_rating     INTEGER,
    conditions      VARCHAR[] NOT NULL,
    evidence_ids    VARCHAR[] NOT NULL,
    run_id          VARCHAR
);

CREATE INDEX IF NOT EXISTS idx_clinical_assertions_variant
    ON clinical_assertions (variant_id);


-- ---------------------------------------------------------------- phenotype

-- Four-valued phenotype logic (GP-14). `status` is never collapsed to a boolean:
-- NOT_ASSESSED and EXCLUDED are stored distinctly because only the latter is
-- usable negative evidence.
--
-- Note the deliberate absence of any free-text clinical note column. Only
-- `source_excerpt_hash` is persisted; the excerpt itself is directly identifying
-- and must not enter a derived artifact (PRIV).
CREATE TABLE IF NOT EXISTS phenotype_observations (
    subject_id             VARCHAR   NOT NULL,
    hpo_id                 VARCHAR   NOT NULL,
    label                  VARCHAR   NOT NULL,
    status                 VARCHAR   NOT NULL,
    onset                  VARCHAR   NOT NULL,
    provenance             VARCHAR   NOT NULL,
    extraction_confidence  DOUBLE    NOT NULL,
    source_excerpt_hash    VARCHAR,
    notes                  VARCHAR,
    source_artifact        VARCHAR   NOT NULL,
    hpo_version            VARCHAR   NOT NULL,
    evidence_ids           VARCHAR[] NOT NULL,
    run_id                 VARCHAR,
    PRIMARY KEY (subject_id, hpo_id)
);


-- --------------------------------------------------------------------- genes

-- A materialised projection of `consequences`, refreshed wholesale by
-- EvidenceStore.write_consequences. It is a convenience index for gene-level
-- questions, never an independent source of truth — which is why it is rebuilt
-- rather than incrementally merged (incremental merge of list columns is the
-- classic way an "idempotent" writer stops being idempotent).
CREATE TABLE IF NOT EXISTS genes (
    gene_symbol     VARCHAR PRIMARY KEY,
    gene_id         VARCHAR,
    transcript_ids  VARCHAR[] NOT NULL,
    variant_ids     VARCHAR[] NOT NULL,
    -- Most severe *assessed* impact across the gene's transcripts, or NULL when no
    -- transcript carries an assessed impact at all. Unassessed rows are excluded
    -- from the minimum rather than bucketed as least-severe, so a
    -- gene-assignment-only annotation can neither dilute a real prediction from
    -- another adapter nor manufacture a 'modifier' verdict of its own (ADR 0016).
    worst_impact    VARCHAR,
    evidence_ids    VARCHAR[] NOT NULL
);


-- ----------------------------------------------------------- candidate pairs

-- The Track 1 unit of prediction. The score is stored as a VECTOR, not just the
-- composite: "why is this ranked first?" must be answerable from the database
-- alone, without re-running the pipeline. `variant_b_id` is NULL for
-- single-variant (dominant / de-novo / homozygous) hypotheses, which share the
-- ranked list with compound-het pairs.
CREATE TABLE IF NOT EXISTS candidate_pairs (
    pair_id                     VARCHAR PRIMARY KEY,
    gene_symbol                 VARCHAR   NOT NULL,
    variant_a_id                VARCHAR   NOT NULL,
    variant_b_id                VARCHAR,
    is_pair                     BOOLEAN   NOT NULL,
    inheritance_model           VARCHAR   NOT NULL,
    -- Phase is never assumed (GP-15). UNKNOWN is stored and preserved, and the
    -- method/distance columns record *why* it could not be determined.
    phase_status                VARCHAR   NOT NULL,
    phase_method                VARCHAR   NOT NULL,
    phase_supporting_reads      INTEGER,
    phase_distance_bp           BIGINT,
    phase_notes                 VARCHAR,
    score_analytical_validity   DOUBLE    NOT NULL,
    score_rarity                DOUBLE    NOT NULL,
    score_molecular_consequence DOUBLE    NOT NULL,
    score_inheritance_consistency DOUBLE  NOT NULL,
    score_phenotype_similarity  DOUBLE    NOT NULL,
    score_mechanistic_relevance DOUBLE    NOT NULL,
    score_evidence_quality      DOUBLE    NOT NULL,
    score_contradiction_penalty DOUBLE    NOT NULL,
    composite_score             DOUBLE    NOT NULL,
    rank                        INTEGER,
    -- Positional tiebreak for a total ordering; see CandidatePair.sort_key.
    contig_order                INTEGER   NOT NULL,
    position                    BIGINT    NOT NULL,
    supporting_evidence_ids     VARCHAR[] NOT NULL,
    -- GP-19: what argues AGAINST this candidate is a first-class column, not a
    -- footnote, and is never pruned when a pair is de-prioritised.
    contradicting_evidence_ids  VARCHAR[] NOT NULL,
    -- Enumerated gaps (OpenQuestion). A candidate ranked first with four open
    -- questions is a different object from one with none, and the report says so.
    missing_evidence            JSON      NOT NULL,
    blocking_question_count     INTEGER   NOT NULL,
    recommended_next_test       VARCHAR   NOT NULL,
    discriminating_experiment   VARCHAR,
    rank_rationale              VARCHAR   NOT NULL,
    flags                       VARCHAR[] NOT NULL,
    run_id                      VARCHAR
);

CREATE INDEX IF NOT EXISTS idx_candidate_pairs_gene ON candidate_pairs (gene_symbol);
CREATE INDEX IF NOT EXISTS idx_candidate_pairs_score ON candidate_pairs (composite_score);


-- ------------------------------------------------------------ evidence items

-- The append-only spine of the whole system. Everything else cites into it.
--
-- `limitations` is NOT NULL by design (GP-17): an evidence item that cannot say
-- what it fails to establish is how a prediction gets mistaken for proof.
-- `direction` is indexed because "show me everything that contradicts this" is a
-- first-class query, not an exception path (GP-19).
CREATE TABLE IF NOT EXISTS evidence_items (
    evidence_id          VARCHAR PRIMARY KEY,
    subject_id           VARCHAR NOT NULL,
    subject_kind         VARCHAR NOT NULL,
    claim                VARCHAR NOT NULL,
    category             VARCHAR NOT NULL,
    direction            VARCHAR NOT NULL,
    strength             VARCHAR NOT NULL,
    evidence_type        VARCHAR NOT NULL,
    tier                 VARCHAR NOT NULL,
    method               VARCHAR NOT NULL,
    tool                 VARCHAR NOT NULL,
    tool_version         VARCHAR NOT NULL,
    limitations          VARCHAR NOT NULL,
    -- Citation flattened into columns (rather than a foreign key) so that a single
    -- SELECT answers "what is this claim and who says so". The normalised copy
    -- lives in `citations` for bibliography rendering.
    citation_source      VARCHAR,
    citation_identifier  VARCHAR,
    citation_version     VARCHAR,
    citation_url         VARCHAR,
    citation_title       VARCHAR,
    citation_key         VARCHAR,
    -- UTC, naive: for SQL predicates. See the header note on TIMESTAMPTZ.
    timestamp            TIMESTAMP NOT NULL,
    -- Exact ISO-8601 of the original aware datetime: the read path reconstructs
    -- EvidenceItem.timestamp from this, so a round trip is lossless.
    timestamp_iso        VARCHAR NOT NULL,
    run_id               VARCHAR,
    numeric_value        DOUBLE,
    -- Structured machine-readable detail. Canonical JSON; never patient identifiers.
    payload              JSON    NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_evidence_subject ON evidence_items (subject_id);
CREATE INDEX IF NOT EXISTS idx_evidence_direction ON evidence_items (direction);
CREATE INDEX IF NOT EXISTS idx_evidence_category ON evidence_items (category);


-- Normalised bibliography. One row per distinct source pointer, so a report can
-- render references without walking every evidence row.
CREATE TABLE IF NOT EXISTS citations (
    citation_key  VARCHAR PRIMARY KEY,
    source        VARCHAR NOT NULL,
    identifier    VARCHAR NOT NULL,
    version       VARCHAR,
    url           VARCHAR,
    title         VARCHAR
);


-- ---------------------------------------------------------------- mechanisms

-- A mechanism is stored as an explicit chain of nodes and links rather than a
-- prose paragraph. That is what makes Track 2 auditable: the database can be
-- asked "which link is only inferred?" and "in which direction must the target
-- node be pushed?", questions a paragraph cannot answer.
CREATE TABLE IF NOT EXISTS mechanisms (
    mechanism_id                 VARCHAR PRIMARY KEY,
    gene_symbol                  VARCHAR   NOT NULL,
    pair_id                      VARCHAR,
    summary                      VARCHAR   NOT NULL,
    -- Signed (GP-16). The therapeutic requirement is the inverse of this, and the
    -- drug direction check compares against `required_correction`.
    disease_direction            VARCHAR   NOT NULL,
    therapeutic_target_node_id   VARCHAR   NOT NULL,
    required_correction          VARCHAR   NOT NULL,
    node_count                   INTEGER   NOT NULL,
    link_count                   INTEGER   NOT NULL,
    inferred_link_count          INTEGER   NOT NULL,
    is_fully_demonstrated        BOOLEAN   NOT NULL,
    supporting_evidence_ids      VARCHAR[] NOT NULL,
    contradicting_evidence_ids   VARCHAR[] NOT NULL,
    uncertainties                VARCHAR[] NOT NULL,
    discriminating_experiments   JSON      NOT NULL,
    developmental_window_caveat  VARCHAR   NOT NULL,
    run_id                       VARCHAR
);


CREATE TABLE IF NOT EXISTS mechanism_nodes (
    mechanism_id      VARCHAR NOT NULL,
    node_id           VARCHAR NOT NULL,
    kind              VARCHAR NOT NULL,
    label             VARCHAR NOT NULL,
    identifier        VARCHAR,
    -- How this node deviates from wild type in the affected individual. Signed.
    state_in_patient  VARCHAR NOT NULL,
    description       VARCHAR NOT NULL,
    PRIMARY KEY (mechanism_id, node_id)
);


CREATE TABLE IF NOT EXISTS mechanism_links (
    mechanism_id                VARCHAR   NOT NULL,
    link_id                     VARCHAR   NOT NULL,
    source_node_id              VARCHAR   NOT NULL,
    target_node_id              VARCHAR   NOT NULL,
    relation                    VARCHAR   NOT NULL,
    direction                   VARCHAR   NOT NULL,
    tier                        VARCHAR   NOT NULL,
    strength                    VARCHAR   NOT NULL,
    -- FALSE means "inferred by analogy or pathway membership". Persisted so the
    -- weak points of a chain can be listed mechanically rather than remembered.
    is_directly_demonstrated    BOOLEAN   NOT NULL,
    uncertainty                 VARCHAR   NOT NULL,
    evidence_ids                VARCHAR[] NOT NULL,
    contradicting_evidence_ids  VARCHAR[] NOT NULL,
    PRIMARY KEY (mechanism_id, link_id)
);

CREATE INDEX IF NOT EXISTS idx_mechanism_links_source ON mechanism_links (source_node_id);
CREATE INDEX IF NOT EXISTS idx_mechanism_links_target ON mechanism_links (target_node_id);


-- --------------------------------------------------------------------- drugs

-- Drug-repurposing hypotheses, accepted AND rejected. Both directions are stored
-- in the same table so that a query for "what did we consider?" cannot
-- accidentally return only the survivors.
CREATE TABLE IF NOT EXISTS drugs (
    drug_id                          VARCHAR PRIMARY KEY,
    name                             VARCHAR   NOT NULL,
    approved_name                    VARCHAR,
    approval_status                  VARCHAR   NOT NULL,
    intervention_class               VARCHAR   NOT NULL,
    target                           VARCHAR   NOT NULL,
    target_node_id                   VARCHAR   NOT NULL,
    mechanism_id                     VARCHAR,
    mechanism_of_action              VARCHAR   NOT NULL,
    -- The three columns the whole module exists to protect. `directions_agree` is
    -- TRI-STATE: NULL means "undeterminable", which is NOT agreement (GP-16).
    required_direction               VARCHAR   NOT NULL,
    observed_direction               VARCHAR   NOT NULL,
    directions_agree                 BOOLEAN,
    is_direct_evidence               BOOLEAN   NOT NULL,
    strongest_evidence_type          VARCHAR   NOT NULL,
    has_in_vivo_evidence             BOOLEAN   NOT NULL,
    is_repurposable                  BOOLEAN   NOT NULL,
    pediatric_has_exposure           BOOLEAN   NOT NULL,
    pediatric_youngest_age_studied   VARCHAR,
    pediatric_indication             VARCHAR,
    pediatric_tolerability_summary   VARCHAR,
    pediatric_caveat                 VARCHAR   NOT NULL,
    pediatric_evidence_ids           VARCHAR[] NOT NULL,
    pk_route                         VARCHAR,
    -- NULL means unknown, which is a blocking gap for a neurological phenotype —
    -- not a permissive default.
    pk_cns_penetrant                 BOOLEAN,
    pk_achievable_plasma_um          DOUBLE,
    pk_required_effective_um         DOUBLE,
    pk_concentration_achievable      BOOLEAN,
    pk_half_life_hours               DOUBLE,
    pk_notes                         VARCHAR   NOT NULL,
    safety_concerns                  JSON      NOT NULL,
    disqualifying_safety_count       INTEGER   NOT NULL,
    -- NULL = unassessed. For a chromosomal-instability disorder that is itself a
    -- blocking gap, so it is stored tri-state rather than defaulted to FALSE.
    worsens_chromosomal_instability  BOOLEAN,
    proposed_validation_experiment   VARCHAR   NOT NULL,
    evidence_ids                     VARCHAR[] NOT NULL,
    contradicting_evidence_ids       VARCHAR[] NOT NULL,
    score                            DOUBLE    NOT NULL,
    rank                             INTEGER,
    rejected                         BOOLEAN   NOT NULL,
    rejection_rationale              VARCHAR   NOT NULL,
    run_id                           VARCHAR
);

CREATE INDEX IF NOT EXISTS idx_drugs_rejected ON drugs (rejected);
CREATE INDEX IF NOT EXISTS idx_drugs_target_node ON drugs (target_node_id);


-- GP-19 made physical. One row per (drug, reason): a rejected candidate keeps its
-- full reason set forever, so a later reviewer can ask "what did the pipeline
-- throw away, and was it right to?". Nothing in this codebase ever DELETEs from
-- here; the store has no API that can.
CREATE TABLE IF NOT EXISTS drug_rejections (
    drug_id     VARCHAR NOT NULL,
    reason      VARCHAR NOT NULL,
    drug_name   VARCHAR NOT NULL,
    rationale   VARCHAR NOT NULL,
    -- Denormalised copies of the direction pair: the single most important thing
    -- to be able to read off a rejection row without a join.
    required_direction  VARCHAR NOT NULL,
    observed_direction  VARCHAR NOT NULL,
    run_id      VARCHAR,
    PRIMARY KEY (drug_id, reason)
);


-- ---------------------------------------------------------------- provenance

-- One row per pipeline execution (GP-31). `config_snapshot` is stored whole so a
-- run can be reproduced from the database alone.
CREATE TABLE IF NOT EXISTS pipeline_runs (
    run_id            VARCHAR PRIMARY KEY,
    case_id           VARCHAR   NOT NULL,
    genome_build      VARCHAR   NOT NULL,
    started_at        TIMESTAMP NOT NULL,
    started_at_iso    VARCHAR   NOT NULL,
    completed_at      TIMESTAMP,
    completed_at_iso  VARCHAR,
    config_hash       VARCHAR   NOT NULL,
    config_snapshot   JSON      NOT NULL,
    git_commit        VARCHAR,
    -- A dirty tree is not reproducible from the commit alone, and is marked as
    -- such rather than quietly recorded as if it were.
    git_dirty         BOOLEAN   NOT NULL,
    is_reproducible   BOOLEAN   NOT NULL,
    inputs            JSON      NOT NULL,
    commands          JSON      NOT NULL,
    tool_versions     JSON      NOT NULL,
    reference_versions JSON     NOT NULL,
    python_version    VARCHAR   NOT NULL,
    platform          VARCHAR   NOT NULL,
    network_profile   VARCHAR   NOT NULL,
    synthetic         BOOLEAN   NOT NULL,
    warnings          VARCHAR[] NOT NULL
);


-- One row per artifact a run produced. `sensitivity` is the machine-readable
-- basis for the public-export gate: classification is a claim, and the export
-- scanner is the verification.
CREATE TABLE IF NOT EXISTS artifact_provenance (
    run_id                VARCHAR   NOT NULL,
    artifact_id           VARCHAR   NOT NULL,
    kind                  VARCHAR   NOT NULL,
    -- Relative to the workspace root, never absolute: an absolute path leaks the
    -- location and naming of patient files into a publishable manifest.
    relative_path         VARCHAR   NOT NULL,
    sensitivity           VARCHAR   NOT NULL,
    is_exportable         BOOLEAN   NOT NULL,
    content_hash          VARCHAR   NOT NULL,
    size_bytes            BIGINT    NOT NULL,
    produced_by_stage     VARCHAR   NOT NULL,
    upstream_artifact_ids VARCHAR[] NOT NULL,
    tool_versions         JSON      NOT NULL,
    created_at            TIMESTAMP NOT NULL,
    created_at_iso        VARCHAR   NOT NULL,
    row_count             BIGINT,
    notes                 VARCHAR   NOT NULL,
    PRIMARY KEY (run_id, artifact_id)
);


-- --------------------------------------------------------------- graph edges

-- The subject–predicate–object table.
--
-- This project deliberately does NOT adopt a graph database. The graph here is
-- small (thousands of edges), is always traversed one or two hops from a known
-- subject, and — decisively — must live in the same transactional store and the
-- same Parquet export as the tabular evidence it cites. A second engine would
-- mean a second determinism story, a second provenance story and a second
-- deployment dependency, in exchange for traversals this workload never performs.
-- Recording the graph relationally keeps one storage engine, one export format
-- and one integrity boundary. (An ADR records this decision in full.)
--
-- `edge_id` is a deterministic hash of (subject_id, predicate, object_id): the
-- same relation asserted twice is one edge, not two, so writes converge.
CREATE TABLE IF NOT EXISTS graph_edges (
    edge_id       VARCHAR PRIMARY KEY,
    subject_id    VARCHAR   NOT NULL,
    subject_kind  VARCHAR   NOT NULL,
    predicate     VARCHAR   NOT NULL,
    object_id     VARCHAR   NOT NULL,
    object_kind   VARCHAR   NOT NULL,
    -- GP-10: an edge is a claim, so it carries its citations like any other.
    evidence_ids  VARCHAR[] NOT NULL,
    confidence    DOUBLE,
    run_id        VARCHAR
);

-- Both directions are indexed: traversal from a subject ("what does this variant
-- affect?") and from an object ("what implicates this phenotype?") are equally
-- common, and an un-indexed reverse lookup is a full scan.
CREATE INDEX IF NOT EXISTS idx_graph_edges_subject ON graph_edges (subject_id);
CREATE INDEX IF NOT EXISTS idx_graph_edges_object ON graph_edges (object_id);
CREATE INDEX IF NOT EXISTS idx_graph_edges_predicate ON graph_edges (predicate);
