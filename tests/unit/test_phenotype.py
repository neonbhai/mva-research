"""Unit tests for phenotype parsing and gene-phenotype scoring.

The centre of gravity here is GP-14. Most of these tests are not checking that a
number is right; they are checking that a *distinction survives*: that
``NOT_ASSESSED`` never behaves like ``EXCLUDED``, that ``UNCERTAIN`` never behaves
like either, and that neither can move a score. The equality assertion in
``test_not_assessed_association_scores_identically_to_no_association`` is the
regression guard the module exists for — if it ever fails, the scorer has started
penalising genes for questions the clinic did not ask.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pytest

from mva.clock import Clock, demo_clock
from mva.determinism import stable_hash
from mva.errors import IngestionError
from mva.models.base import AssertionTier
from mva.models.evidence import EvidenceCategory, EvidenceDirection, EvidenceType
from mva.models.phenotype import (
    NEGATIVE_EVIDENCE_STATUSES,
    UNINFORMATIVE_STATUSES,
    ObservationStatus,
    Onset,
    PhenotypeProfile,
)
from mva.phenotype import (
    NEUTRAL_SCORE,
    STATUS_ALIASES,
    STRENGTH_WEIGHTS,
    GeneAssociation,
    GenePhenotypeIndex,
    PhenotypeMatch,
    load_phenotype_profile,
    normalise_hpo_id,
    parse_onset,
    parse_status,
    score_all_genes,
    score_gene_phenotype,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
PHENOTYPE_FIXTURE = REPO_ROOT / "tests" / "fixtures" / "synthetic" / "synthetic_phenotype.tsv"
GENE_PHENOTYPE_KNOWLEDGE = REPO_ROOT / "knowledge" / "public" / "gene_phenotype.tsv"

HPO_VERSION = "v0.0-synthetic-hpo"
KNOWLEDGE_VERSION = "v0.0-synthetic"
SOURCE_ARTIFACT = "synthetic_phenotype.tsv"

MICROCEPHALY = "HP:0000252"
IUGR = "HP:0001511"
PREMATURE_CHROMATID_SEPARATION = "HP:0200024"
DANDY_WALKER = "HP:0001305"
CATARACT = "HP:0000518"
NEPHROBLASTOMA = "HP:0002667"
SEIZURE = "HP:0001250"

TSV_HEADER = ("hpo_id", "label", "status", "onset", "provenance", "extraction_confidence", "notes")


# ---------------------------------------------------------------------------
# Helpers and fixtures
# ---------------------------------------------------------------------------


def _assoc(gene: str, hpo_id: str, label: str, strength: str) -> GeneAssociation:
    return GeneAssociation(
        gene_symbol=gene,
        hpo_id=hpo_id,
        label=label,
        association_strength=strength,
        source="unit_test_curation",
        version="v0",
    )


def _index(associations: Sequence[GeneAssociation]) -> GenePhenotypeIndex:
    return GenePhenotypeIndex(associations, version=KNOWLEDGE_VERSION)


def _write_phenotype_tsv(path: Path, rows: Sequence[Sequence[str]]) -> Path:
    lines = ["# mva_synthetic=true  fabricated unit-test profile", "\t".join(TSV_HEADER)]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _row(hpo_id: str, label: str, status: str) -> tuple[str, ...]:
    return (hpo_id, label, status, "unknown", "unit_test", "1.0", "")


def _fingerprint(matches: dict[str, PhenotypeMatch]) -> str:
    """Total, order-independent content hash of a scoring result."""
    return stable_hash(
        {
            gene: {
                "score": match.score,
                "matched": list(match.matched_terms),
                "contradicted": list(match.contradicted_terms),
                "unassessed": list(match.unassessed_terms),
                "unexplained": list(match.unexplained_terms),
                "rationale": match.rationale,
                "evidence": [item.model_dump(mode="json") for item in match.evidence],
            }
            for gene, match in matches.items()
        }
    )


@pytest.fixture
def clock() -> Clock:
    return demo_clock()


@pytest.fixture
def profile() -> PhenotypeProfile:
    return load_phenotype_profile(
        PHENOTYPE_FIXTURE,
        subject_id="PROBAND01",
        hpo_version=HPO_VERSION,
        source_artifact=SOURCE_ARTIFACT,
    )


@pytest.fixture
def knowledge_index() -> GenePhenotypeIndex:
    return GenePhenotypeIndex.from_tsv(GENE_PHENOTYPE_KNOWLEDGE, version=KNOWLEDGE_VERSION)


# ---------------------------------------------------------------------------
# 1. Negation parsing
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize(
    "spelling",
    [
        "excluded",
        "EXCLUDED",
        "  Excluded  ",
        "absent",
        "Absent",
        "negated",
        "negative",
        "no",
        "not_present",
        "not present",
        "NOT-PRESENT",
        "not_observed",
        "ruled out",
        "ruled_out",
        "denied",
    ],
)
def test_negation_spellings_parse_as_excluded(spelling: str) -> None:
    """Real extracts spell negation many ways; all of them mean EXCLUDED."""
    status = parse_status(spelling)
    assert status is ObservationStatus.EXCLUDED
    assert status in NEGATIVE_EVIDENCE_STATUSES
    assert status not in UNINFORMATIVE_STATUSES


@pytest.mark.unit
@pytest.mark.parametrize(
    "spelling",
    ["not_assessed", "Not Assessed", "NOT-ASSESSED", "not evaluated", "unknown", "n/a", "na", "nd"],
)
def test_no_information_spellings_never_become_excluded(spelling: str) -> None:
    """The dangerous confusion, asserted directly: no information is not negation."""
    status = parse_status(spelling)
    assert status is ObservationStatus.NOT_ASSESSED
    assert status not in NEGATIVE_EVIDENCE_STATUSES
    assert status in UNINFORMATIVE_STATUSES


@pytest.mark.unit
def test_negation_spelling_survives_the_loader(tmp_path: Path) -> None:
    """The mapping is not just in `parse_status`; it reaches the profile."""
    path = _write_phenotype_tsv(
        tmp_path / "negation.tsv",
        [
            _row(MICROCEPHALY, "Microcephaly", "present"),
            _row(CATARACT, "Cataract", "NOT-PRESENT"),
            _row(SEIZURE, "Seizure", "equivocal"),
            _row(NEPHROBLASTOMA, "Nephroblastoma", "n/a"),
        ],
    )
    loaded = load_phenotype_profile(
        path, subject_id="S1", hpo_version=HPO_VERSION, source_artifact="negation.tsv"
    )
    assert loaded.observed_terms == (MICROCEPHALY,)
    assert loaded.excluded_terms == (CATARACT,)
    assert loaded.not_assessed_terms == (NEPHROBLASTOMA,)
    assert loaded.status_of(SEIZURE) is ObservationStatus.UNCERTAIN


@pytest.mark.unit
def test_every_alias_maps_to_a_real_status() -> None:
    """The alias table cannot drift away from the four-valued enum."""
    assert set(STATUS_ALIASES.values()) == set(ObservationStatus)


# ---------------------------------------------------------------------------
# 2. THE regression guard: missing vs explicitly absent
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_not_assessed_association_scores_identically_to_no_association(
    profile: PhenotypeProfile, clock: Clock
) -> None:
    """A gene must not be penalised for a term nobody assessed.

    This is the core guarantee of the module. If the denominator were the gene's
    total associated weight, adding a NOT_ASSESSED association would drag the
    score down; here it must change nothing at all.
    """
    base = [
        _assoc("SYNTHKIN1", MICROCEPHALY, "Microcephaly", "strong"),
        _assoc("SYNTHKIN1", IUGR, "Intrauterine growth retardation", "strong"),
        _assoc(
            "SYNTHKIN1", PREMATURE_CHROMATID_SEPARATION, "Premature chromatid sep.", "definitive"
        ),
        _assoc("SYNTHKIN1", DANDY_WALKER, "Dandy-Walker malformation", "moderate"),
    ]
    with_gap = [*base, _assoc("SYNTHKIN1", NEPHROBLASTOMA, "Nephroblastoma", "definitive")]

    without = score_gene_phenotype("SYNTHKIN1", profile=profile, index=_index(base), clock=clock)
    with_unassessed = score_gene_phenotype(
        "SYNTHKIN1", profile=profile, index=_index(with_gap), clock=clock
    )

    assert with_unassessed.score == without.score
    assert with_unassessed.matched_terms == without.matched_terms
    assert with_unassessed.contradicted_terms == without.contradicted_terms
    assert with_unassessed.unexplained_terms == without.unexplained_terms

    # The gap is not silently dropped: it is recorded, just not scored.
    assert without.unassessed_terms == ()
    assert with_unassessed.unassessed_terms == (NEPHROBLASTOMA,)


@pytest.mark.unit
def test_many_unassessed_associations_still_do_not_move_the_score(
    profile: PhenotypeProfile, clock: Clock
) -> None:
    """The invariant must not decay as the number of unassessed terms grows."""
    base = [_assoc("GENEX", MICROCEPHALY, "Microcephaly", "strong")]
    padded = [
        *base,
        *(
            _assoc("GENEX", f"HP:900000{n}", f"Unassessed feature {n}", "definitive")
            for n in range(1, 8)
        ),
    ]
    lean = score_gene_phenotype("GENEX", profile=profile, index=_index(base), clock=clock)
    padded_match = score_gene_phenotype("GENEX", profile=profile, index=_index(padded), clock=clock)

    assert padded_match.score == lean.score
    assert len(padded_match.unassessed_terms) == 7


@pytest.mark.unit
def test_unassessed_term_emits_a_neutral_information_gap_item(
    profile: PhenotypeProfile, clock: Clock
) -> None:
    """A gap is reported as an actionable NEUTRAL item, never as absence."""
    match = score_gene_phenotype(
        "SYNTHKIN1",
        profile=profile,
        index=_index(
            [
                _assoc("SYNTHKIN1", MICROCEPHALY, "Microcephaly", "strong"),
                _assoc("SYNTHKIN1", NEPHROBLASTOMA, "Nephroblastoma", "moderate"),
            ]
        ),
        clock=clock,
    )
    gaps = [item for item in match.evidence if NEPHROBLASTOMA in item.claim]
    assert len(gaps) == 1
    gap = gaps[0]
    assert gap.direction is EvidenceDirection.NEUTRAL
    assert gap.numeric_value == 0.0
    assert gap.category is EvidenceCategory.PHENOTYPE
    assert "not absence of the feature" in gap.limitations


# ---------------------------------------------------------------------------
# 3. Excluded terms are real negative evidence
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_excluded_term_lowers_score_and_records_a_contradiction(
    profile: PhenotypeProfile, clock: Clock
) -> None:
    base = [_assoc("GENEY", MICROCEPHALY, "Microcephaly", "strong")]
    with_excluded = [*base, _assoc("GENEY", CATARACT, "Cataract", "strong")]

    without = score_gene_phenotype("GENEY", profile=profile, index=_index(base), clock=clock)
    with_contra = score_gene_phenotype(
        "GENEY", profile=profile, index=_index(with_excluded), clock=clock
    )

    assert with_contra.score < without.score
    assert with_contra.contradicted_terms == (CATARACT,)
    assert with_contra.has_contradiction

    contradictions = [
        item
        for item in with_contra.evidence
        if item.direction is EvidenceDirection.CONTRADICTS and CATARACT in item.claim
    ]
    assert len(contradictions) == 1
    item = contradictions[0]
    assert item.category is EvidenceCategory.PHENOTYPE
    assert item.tier is AssertionTier.OBSERVED_DATA
    assert item.evidence_type is EvidenceType.DIRECT_MEASUREMENT
    assert item.is_contradiction
    assert item.limitations.strip()


@pytest.mark.unit
def test_excluded_is_the_only_status_that_can_push_below_neutral(
    profile: PhenotypeProfile, clock: Clock
) -> None:
    """Only an assessed absence may drag a gene below the no-information value."""
    contradicted = score_gene_phenotype(
        "GENEZ",
        profile=profile,
        index=_index([_assoc("GENEZ", CATARACT, "Cataract", "definitive")]),
        clock=clock,
    )
    unassessed = score_gene_phenotype(
        "GENEZ",
        profile=profile,
        index=_index([_assoc("GENEZ", NEPHROBLASTOMA, "Nephroblastoma", "definitive")]),
        clock=clock,
    )
    uncertain = score_gene_phenotype(
        "GENEZ",
        profile=profile,
        index=_index([_assoc("GENEZ", SEIZURE, "Seizure", "definitive")]),
        clock=clock,
    )

    assert contradicted.score < NEUTRAL_SCORE
    assert unassessed.score == NEUTRAL_SCORE
    assert uncertain.score == NEUTRAL_SCORE


# ---------------------------------------------------------------------------
# 4. Uncertain contributes zero
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_uncertain_term_contributes_zero(profile: PhenotypeProfile, clock: Clock) -> None:
    """An equivocal finding is neither presence nor absence."""
    base = [_assoc("SYNTHMET2", MICROCEPHALY, "Microcephaly", "supporting")]
    with_uncertain = [*base, _assoc("SYNTHMET2", SEIZURE, "Seizure", "definitive")]

    without = score_gene_phenotype("SYNTHMET2", profile=profile, index=_index(base), clock=clock)
    with_equivocal = score_gene_phenotype(
        "SYNTHMET2", profile=profile, index=_index(with_uncertain), clock=clock
    )

    assert with_equivocal.score == without.score
    assert SEIZURE not in with_equivocal.matched_terms
    assert SEIZURE not in with_equivocal.contradicted_terms
    assert SEIZURE not in with_equivocal.unassessed_terms
    assert SEIZURE not in with_equivocal.unexplained_terms
    assert "uncertain" in with_equivocal.rationale.lower()
    assert SEIZURE in with_equivocal.rationale

    neutral_items = [
        item
        for item in with_equivocal.evidence
        if SEIZURE in item.claim and item.direction is EvidenceDirection.NEUTRAL
    ]
    assert len(neutral_items) == 1
    assert neutral_items[0].numeric_value == 0.0


# ---------------------------------------------------------------------------
# 5. Unknown terms default to NOT_ASSESSED
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_status_of_unknown_term_is_not_assessed(profile: PhenotypeProfile) -> None:
    """A term nobody recorded is a term nobody assessed."""
    assert profile.status_of("HP:9999999") is ObservationStatus.NOT_ASSESSED
    assert profile.status_of("HP:9999999") not in NEGATIVE_EVIDENCE_STATUSES
    assert profile.status_of(NEPHROBLASTOMA) is ObservationStatus.NOT_ASSESSED
    assert profile.status_of(CATARACT) is ObservationStatus.EXCLUDED
    assert profile.status_of(MICROCEPHALY) is ObservationStatus.OBSERVED
    assert profile.status_of(SEIZURE) is ObservationStatus.UNCERTAIN


# ---------------------------------------------------------------------------
# 6. Invalid HPO identifiers
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize(
    "bad", ["HP:123", "HP:00002520", "0000252", "HPO:0000252", "", "HP:abcdefg"]
)
def test_invalid_hpo_id_raises(bad: str) -> None:
    with pytest.raises(IngestionError) as excinfo:
        normalise_hpo_id(bad, context="unit test")
    assert "HP:0001250" in str(excinfo.value)


@pytest.mark.unit
def test_invalid_hpo_id_in_profile_file_raises(tmp_path: Path) -> None:
    path = _write_phenotype_tsv(
        tmp_path / "bad_hpo.tsv", [_row("HP:123", "Microcephaly", "observed")]
    )
    with pytest.raises(IngestionError) as excinfo:
        load_phenotype_profile(
            path, subject_id="S1", hpo_version=HPO_VERSION, source_artifact="bad_hpo.tsv"
        )
    assert "bad_hpo.tsv" in str(excinfo.value)


@pytest.mark.unit
def test_invalid_hpo_id_in_association_raises() -> None:
    with pytest.raises(IngestionError):
        _assoc("GENEA", "HP:123", "Microcephaly", "strong")


@pytest.mark.unit
def test_obo_style_hpo_id_is_normalised() -> None:
    assert normalise_hpo_id("hp_0000252", context="unit test") == MICROCEPHALY


# ---------------------------------------------------------------------------
# 7. Unrecognised status is a loud error
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("bad", ["mostly there", "", "   ", "definitely-maybe", "1"])
def test_unrecognised_status_raises_with_helpful_message(bad: str) -> None:
    with pytest.raises(IngestionError) as excinfo:
        parse_status(bad)
    message = str(excinfo.value)
    for canonical in ("observed", "excluded", "uncertain", "not_assessed"):
        assert canonical in message
    assert "GP-14" in message
    assert "absent" in message  # the accepted-spellings listing is included


@pytest.mark.unit
def test_unrecognised_status_in_file_names_the_file_and_survives_to_the_caller(
    tmp_path: Path,
) -> None:
    path = _write_phenotype_tsv(
        tmp_path / "bad_status.tsv", [_row(MICROCEPHALY, "Microcephaly", "probably-ish")]
    )
    with pytest.raises(IngestionError) as excinfo:
        load_phenotype_profile(
            path, subject_id="S1", hpo_version=HPO_VERSION, source_artifact="bad_status.tsv"
        )
    assert "probably-ish" in str(excinfo.value)


@pytest.mark.unit
def test_unrecognised_onset_raises() -> None:
    with pytest.raises(IngestionError):
        parse_onset("perinatal-ish")
    assert parse_onset("") is Onset.UNKNOWN
    assert parse_onset("Congenital") is Onset.CONGENITAL


@pytest.mark.unit
def test_unknown_association_strength_raises() -> None:
    with pytest.raises(IngestionError) as excinfo:
        _assoc("GENEB", MICROCEPHALY, "Microcephaly", "quite good")
    message = str(excinfo.value)
    for allowed in sorted(STRENGTH_WEIGHTS):
        assert allowed in message


@pytest.mark.unit
def test_duplicate_hpo_term_in_profile_raises(tmp_path: Path) -> None:
    path = _write_phenotype_tsv(
        tmp_path / "dupe.tsv",
        [
            _row(MICROCEPHALY, "Microcephaly", "observed"),
            _row(MICROCEPHALY, "Microcephaly", "excluded"),
        ],
    )
    with pytest.raises(IngestionError) as excinfo:
        load_phenotype_profile(
            path, subject_id="S1", hpo_version=HPO_VERSION, source_artifact="dupe.tsv"
        )
    assert "duplicate" in str(excinfo.value).lower()


@pytest.mark.unit
def test_missing_required_column_raises(tmp_path: Path) -> None:
    path = tmp_path / "no_status.tsv"
    path.write_text("hpo_id\tlabel\nHP:0000252\tMicrocephaly\n", encoding="utf-8")
    with pytest.raises(IngestionError) as excinfo:
        load_phenotype_profile(
            path, subject_id="S1", hpo_version=HPO_VERSION, source_artifact="no_status.tsv"
        )
    assert "status" in str(excinfo.value)


# ---------------------------------------------------------------------------
# 8. Determinism (GP-30)
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_scoring_is_deterministic_across_runs(
    profile: PhenotypeProfile, knowledge_index: GenePhenotypeIndex
) -> None:
    genes = ["SYNTHOTH3", "SYNTHKIN1", "SYNTHMET2", "SYNTHMUL4", "SYNTHSOL5", "NOTINDEXED"]
    first = score_all_genes(genes, profile=profile, index=knowledge_index, clock=demo_clock())
    second = score_all_genes(
        list(reversed(genes)), profile=profile, index=knowledge_index, clock=demo_clock()
    )

    assert list(first) == sorted(genes)
    assert list(first) == list(second)
    assert _fingerprint(first) == _fingerprint(second)


@pytest.mark.unit
def test_index_accessors_are_ordered(knowledge_index: GenePhenotypeIndex) -> None:
    terms = knowledge_index.terms_for_gene("SYNTHKIN1")
    assert [assoc.hpo_id for assoc in terms] == sorted(assoc.hpo_id for assoc in terms)
    assert knowledge_index.terms_for_gene("synthkin1") == terms  # case-insensitive lookup
    assert knowledge_index.terms_for_gene("NOT_IN_INDEX") == ()

    genes = knowledge_index.genes_for_term(SEIZURE)
    assert genes == ("SYNTHMET2", "SYNTHMUL4")
    assert knowledge_index.genes_for_term("HP:9999999") == ()
    assert knowledge_index.version == KNOWLEDGE_VERSION


@pytest.mark.unit
def test_duplicate_association_is_rejected() -> None:
    duplicated = [
        _assoc("GENEC", MICROCEPHALY, "Microcephaly", "strong"),
        _assoc("GENEC", MICROCEPHALY, "Microcephaly", "moderate"),
    ]
    with pytest.raises(IngestionError):
        _index(duplicated)


# ---------------------------------------------------------------------------
# 9. Fixture end-to-end
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_fixture_profile_parses_with_all_four_statuses(profile: PhenotypeProfile) -> None:
    assert profile.subject_id == "PROBAND01"
    assert len(profile.observations) == 7
    assert profile.observed_terms == (
        MICROCEPHALY,
        DANDY_WALKER,
        IUGR,
        PREMATURE_CHROMATID_SEPARATION,
    )
    assert profile.excluded_terms == (CATARACT,)
    assert profile.not_assessed_terms == (NEPHROBLASTOMA,)
    assert profile.status_of(SEIZURE) is ObservationStatus.UNCERTAIN
    assert [obs.hpo_id for obs in profile.observations] == sorted(
        obs.hpo_id for obs in profile.observations
    )


@pytest.mark.unit
def test_fixture_ranks_kin1_above_met2_above_oth3(
    profile: PhenotypeProfile, knowledge_index: GenePhenotypeIndex, clock: Clock
) -> None:
    """The headline expectation for the synthetic case."""
    matches = score_all_genes(
        ["SYNTHKIN1", "SYNTHMET2", "SYNTHOTH3"],
        profile=profile,
        index=knowledge_index,
        clock=clock,
    )
    kin1, met2, oth3 = matches["SYNTHKIN1"], matches["SYNTHMET2"], matches["SYNTHOTH3"]

    assert kin1.score > met2.score > oth3.score
    assert oth3.score >= 0.0
    assert kin1.score <= 1.0

    # SYNTHKIN1 explains the whole observed picture and is not dragged down by the
    # nephroblastoma nobody assessed.
    assert kin1.matched_terms == (MICROCEPHALY, DANDY_WALKER, IUGR, PREMATURE_CHROMATID_SEPARATION)
    assert kin1.unassessed_terms == (NEPHROBLASTOMA,)
    assert kin1.unexplained_terms == ()
    assert kin1.score > 0.9

    # SYNTHMET2 matches microcephaly only; its seizure association is UNCERTAIN.
    assert met2.matched_terms == (MICROCEPHALY,)
    assert met2.contradicted_terms == ()
    assert NEUTRAL_SCORE < met2.score < kin1.score

    # SYNTHOTH3's only association was explicitly excluded after assessment.
    assert oth3.contradicted_terms == (CATARACT,)
    assert oth3.score < NEUTRAL_SCORE


@pytest.mark.unit
def test_gene_with_no_phenotype_knowledge_is_neutral_not_zero(
    profile: PhenotypeProfile, knowledge_index: GenePhenotypeIndex, clock: Clock
) -> None:
    """Ignorance about a gene must not rank below evidence against a gene."""
    unknown = score_gene_phenotype(
        "GENE_NOT_IN_INDEX", profile=profile, index=knowledge_index, clock=clock
    )
    contradicted = score_gene_phenotype(
        "SYNTHOTH3", profile=profile, index=knowledge_index, clock=clock
    )
    assert unknown.score == NEUTRAL_SCORE
    assert unknown.score > contradicted.score
    assert "absence of information is not evidence" in unknown.rationale.lower()


# ---------------------------------------------------------------------------
# Evidence hygiene (GP-10 / GP-17)
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_every_evidence_item_is_phenotype_categorised_and_states_limitations(
    profile: PhenotypeProfile, knowledge_index: GenePhenotypeIndex, clock: Clock
) -> None:
    matches = score_all_genes(
        list(knowledge_index.gene_symbols), profile=profile, index=knowledge_index, clock=clock
    )
    seen_ids: set[str] = set()
    for match in matches.values():
        assert match.evidence, f"{match.gene_symbol} produced a score with no evidence (GP-10)"
        for item in match.evidence:
            assert item.category is EvidenceCategory.PHENOTYPE
            assert item.subject_id == match.gene_symbol
            assert item.subject_kind == "gene"
            assert len(item.limitations.strip()) > 3
            assert item.tool_version
            assert item.citation is not None
            assert item.evidence_id not in seen_ids
            seen_ids.add(item.evidence_id)
            if item.tier is AssertionTier.OBSERVED_DATA:
                assert item.evidence_type is EvidenceType.DIRECT_MEASUREMENT


@pytest.mark.unit
def test_computed_match_is_an_inference_tier_item(
    profile: PhenotypeProfile, knowledge_index: GenePhenotypeIndex, clock: Clock
) -> None:
    match = score_gene_phenotype("SYNTHKIN1", profile=profile, index=knowledge_index, clock=clock)
    inferences = [
        item
        for item in match.evidence
        if item.tier is AssertionTier.INFERENCE
        and item.evidence_type is EvidenceType.PIPELINE_INFERENCE
        and item.numeric_value == match.score
    ]
    assert len(inferences) == 1
    assert inferences[0].direction is EvidenceDirection.SUPPORTS
    assert "curation judgement" in inferences[0].limitations
