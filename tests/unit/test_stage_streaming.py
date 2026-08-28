"""Tests for the three changes that let this pipeline survive a WGS callset.

Streaming annotation, a spilling evidence ledger and a named selection stage are
all *equivalence* claims: the new path must produce what the old path produced,
byte for byte, or the change is a silent behaviour edit wearing a performance
costume. Almost every test here is therefore a comparison rather than an
assertion about a literal, and the ones that are not are the GP-14 / ADR 0016
doors that must stay shut.

See ``docs/decisions/0019-selection-before-pairing.md`` and
``docs/scale-report.md`` §9.
"""

from __future__ import annotations

import json
import sqlite3
import subprocess
import sys
from collections.abc import Iterator, Mapping, Sequence
from datetime import UTC, datetime
from pathlib import Path

import pytest

from mva.annotation.base import AdapterSet
from mva.annotation.service import annotate_variants, iter_annotated
from mva.clock import Clock, FixedClock
from mva.config import FrequencyThresholds, QualityThresholds, Workspace
from mva.determinism import canonical_json, stable_hash, write_canonical_json_rows
from mva.errors import AnnotationError, EvidenceError
from mva.evidence.ledger import EvidenceLedger
from mva.evidence.spill import (
    DICTIONARY_BYTES,
    SqliteEvidenceSpill,
    build_dictionary,
    compress,
    decode,
    decompress,
    encode,
)
from mva.evidence.store import EvidenceStore
from mva.models.base import AssertionTier
from mva.models.evidence import (
    Citation,
    EvidenceCategory,
    EvidenceDirection,
    EvidenceItem,
    EvidenceStrength,
    EvidenceType,
    make_evidence_id,
)
from mva.models.genome import GenomeBuild, GenomicCoordinate
from mva.models.provenance import ArtifactKind
from mva.models.variant import (
    ClinicalAssertion,
    ConsequenceAnnotation,
    FilterStatus,
    Genotype,
    ImpactSeverity,
    PopulationFrequency,
    VariantRecord,
    Zygosity,
)
from mva.pipeline import JsonRowsSink, RunContext
from mva.prioritization.filters import (
    apply_hard_filters,
    apply_soft_flags,
    iter_hard_filtered,
    iter_soft_flagged,
)
from mva.prioritization.selection import (
    DROP_COMMON_IN_POPULATION,
    DROP_NO_GENE_ASSIGNMENT,
    DROP_NOT_CODING_OR_SPLICE,
    DROP_REASONS,
    KEEP_CLINICAL_ASSERTION,
    KEEP_CODING_OR_SPLICE,
    KEEP_IMPACT_NOT_ASSESSED,
    KEEP_REASONS,
    NOTE_FREQUENCY_UNDERPOWERED,
    NOTE_FREQUENCY_UNKNOWN,
    SelectionThresholds,
    iter_selected,
    select_variants,
)

REPO_ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture
def clock() -> Clock:
    return FixedClock(datetime(2026, 1, 1, tzinfo=UTC))


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------


def record(
    contig: str = "chr1",
    position: int = 1_000,
    *,
    ref: str = "A",
    alt: str = "G",
    consequences: Sequence[ConsequenceAnnotation] = (),
    frequencies: Sequence[PopulationFrequency] = (),
    assertions: Sequence[ClinicalAssertion] = (),
    zygosity: Zygosity = Zygosity.HET,
) -> VariantRecord:
    return VariantRecord(
        coordinate=GenomicCoordinate(
            build=GenomeBuild.GRCH38, contig=contig, position=position, ref=ref, alt=alt
        ),
        genotype=Genotype(
            zygosity=zygosity,
            genotype_string="1/1" if zygosity is Zygosity.HOM_ALT else "0/1",
            depth=40,
            ref_reads=20,
            alt_reads=20,
            genotype_quality=99,
        ),
        filter_status=FilterStatus.PASS,
        consequences=tuple(consequences),
        population_frequencies=tuple(frequencies),
        clinical_assertions=tuple(assertions),
        source_artifact="tests/unit/test_stage_streaming.py",
    )


def consequence(
    *,
    gene: str = "SYNTHG1",
    terms: Sequence[str] = ("missense_variant",),
    impact: ImpactSeverity | None = ImpactSeverity.MODERATE,
    transcript: str = "ENST00000000001.1",
    splice_ai: float | None = None,
) -> ConsequenceAnnotation:
    return ConsequenceAnnotation(
        gene_symbol=gene,
        transcript_id=transcript,
        consequence_terms=tuple(terms),
        impact=impact,
        splice_ai_delta_max=splice_ai,
        source_tool="test",
        source_tool_version="0",
    )


def frequency(value: float, *, allele_number: int | None = 152_312) -> PopulationFrequency:
    return PopulationFrequency(
        source="test_gnomad",
        version="v0",
        population="global",
        allele_frequency=value,
        allele_number=allele_number,
    )


def evidence_item(
    subject_id: str,
    *,
    claim: str = "A claim about a subject.",
    direction: EvidenceDirection = EvidenceDirection.NEUTRAL,
    category: EvidenceCategory = EvidenceCategory.ANALYTICAL,
    subject_kind: str = "variant",
    run_id: str | None = None,
) -> EvidenceItem:
    return EvidenceItem(
        evidence_id=make_evidence_id(
            subject_id=subject_id, category=category, claim=claim, tool="test"
        ),
        subject_id=subject_id,
        subject_kind=subject_kind,
        claim=claim,
        category=category,
        direction=direction,
        strength=EvidenceStrength.MODERATE,
        evidence_type=EvidenceType.PIPELINE_INFERENCE,
        tier=AssertionTier.INFERENCE,
        citation=Citation(source="test", identifier=subject_id, version="v0"),
        method="Constructed by a test.",
        tool="test",
        tool_version="0",
        limitations="Fabricated for a test; establishes nothing.",
        timestamp=datetime(2026, 1, 1, tzinfo=UTC),
        run_id=run_id,
        numeric_value=1.5,
        payload={"subject": subject_id, "flag": True, "empty": None},
    )


# ---------------------------------------------------------------------------
# Fake adapters, exercising every branch the service has
# ---------------------------------------------------------------------------


class FakeConsequenceAdapter:
    name = "fake-consequence"
    version = "0.1"
    synthetic = True

    def __init__(self, index: Mapping[str, tuple[ConsequenceAnnotation, ...]]) -> None:
        self._index = index
        self.calls: list[tuple[str, ...]] = []

    def annotate(
        self, variant_ids: Sequence[str]
    ) -> Mapping[str, tuple[ConsequenceAnnotation, ...]]:
        self.calls.append(tuple(variant_ids))
        return {k: v for k, v in self._index.items() if k in set(variant_ids)}


class FakeFrequencyAdapter:
    name = "fake-frequency"
    version = "0.1"
    synthetic = False

    def __init__(self, index: Mapping[str, tuple[PopulationFrequency, ...]]) -> None:
        self._index = index
        self.calls: list[tuple[str, ...]] = []

    def frequencies(
        self, variant_ids: Sequence[str]
    ) -> Mapping[str, tuple[PopulationFrequency, ...]]:
        self.calls.append(tuple(variant_ids))
        return {k: v for k, v in self._index.items() if k in set(variant_ids)}


class FakeClinicalAdapter:
    name = "fake-clinical"
    version = "0.1"
    synthetic = True

    def __init__(self, index: Mapping[str, tuple[ClinicalAssertion, ...]]) -> None:
        self._index = index

    def assertions(self, variant_ids: Sequence[str]) -> Mapping[str, tuple[ClinicalAssertion, ...]]:
        return {k: v for k, v in self._index.items() if k in set(variant_ids)}


def sample_records() -> list[VariantRecord]:
    """A callset covering every branch: multi-transcript, no-frequency, pre-annotated."""
    return [
        record("chr1", 100),
        record("chr1", 200),
        record("chr2", 300),
        record("chr2", 400),
        record("chr3", 500, consequences=[consequence(gene="PREANNOTATED")]),
        record("chr4", 600),
        record("chr5", 700),
    ]


def sample_adapters(*, clinical: bool = True) -> AdapterSet:
    records = sample_records()
    ids = [r.variant_id for r in records]
    return AdapterSet(
        consequence=FakeConsequenceAdapter(
            {
                ids[0]: (
                    consequence(transcript="ENST1"),
                    consequence(transcript="ENST2", impact=ImpactSeverity.HIGH),
                ),
                ids[2]: (consequence(gene="SYNTHG2", impact=None),),
                ids[4]: (consequence(gene="SYNTHG3"),),
            }
        ),
        frequency=FakeFrequencyAdapter(
            {
                ids[0]: (frequency(0.0),),
                ids[1]: (frequency(0.3), frequency(0.02)),
                ids[3]: (frequency(1e-6),),
            }
        ),
        clinical=(
            FakeClinicalAdapter(
                {
                    ids[2]: (
                        ClinicalAssertion(
                            source="test", version="v0", significance="Pathogenic", star_rating=2
                        ),
                    )
                }
            )
            if clinical
            else None
        ),
    )


# ---------------------------------------------------------------------------
# Blocker 1: streaming annotation
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize("batch_size", [1, 2, 3, 7, 100])
@pytest.mark.parametrize("clinical", [True, False])
def test_streaming_annotation_equals_the_batch_path(
    batch_size: int, clinical: bool, clock: Clock
) -> None:
    """The equivalence the whole change rests on: same records, same evidence, same report.

    Batching is a memory decision. If it were also a semantics decision — different
    coverage denominators, a different evidence order, a warning that appears only
    at one batch size — the streaming path would be a different stage wearing the
    same name, and no measurement of it would mean anything.
    """
    records = sample_records()
    expected = annotate_variants(records, adapters=sample_adapters(clinical=clinical), clock=clock)

    stream = iter_annotated(
        iter(records),
        adapters=sample_adapters(clinical=clinical),
        clock=clock,
        batch_size=batch_size,
    )
    variants: list[VariantRecord] = []
    evidence: list[EvidenceItem] = []
    for annotated in stream:
        variants.append(annotated.variant)
        evidence.extend(annotated.evidence)

    assert tuple(variants) == expected.variants
    assert tuple(evidence) == expected.evidence
    assert stream.coverage() == expected.coverage
    assert stream.warnings() == expected.warnings
    assert stable_hash([v.model_dump(mode="json") for v in variants]) == stable_hash(
        [v.model_dump(mode="json") for v in expected.variants]
    )


@pytest.mark.unit
def test_streaming_annotation_calls_the_adapters_in_batches_not_per_variant(
    clock: Clock,
) -> None:
    """Adapters are bulk lookups. 4.5 M round trips is not a smaller version of 900."""
    records = sample_records()
    adapters = sample_adapters()
    consequence_adapter = adapters.consequence
    assert isinstance(consequence_adapter, FakeConsequenceAdapter)

    list(iter_annotated(iter(records), adapters=adapters, clock=clock, batch_size=3))

    assert [len(call) for call in consequence_adapter.calls] == [3, 3, 1]


@pytest.mark.unit
def test_streaming_annotation_never_holds_the_callset(clock: Clock) -> None:
    """The source is consumed lazily: nothing pulls more than one batch ahead."""
    pulled = 0

    def counted() -> Iterator[VariantRecord]:
        nonlocal pulled
        for item in sample_records():
            pulled += 1
            yield item

    stream = iter_annotated(counted(), adapters=sample_adapters(), clock=clock, batch_size=2)
    iterator = iter(stream)
    next(iterator)
    assert pulled == 2, "one batch was read, not the whole callset"


@pytest.mark.unit
def test_re_iterating_an_annotation_stream_raises(clock: Clock) -> None:
    """A second pass would re-call every adapter and double every counter."""
    stream = iter_annotated(iter(sample_records()), adapters=sample_adapters(), clock=clock)
    list(stream)
    with pytest.raises(AnnotationError, match="already been iterated"):
        list(stream)


@pytest.mark.unit
def test_a_zero_batch_size_is_refused(clock: Clock) -> None:
    with pytest.raises(AnnotationError, match="at least 1"):
        iter_annotated(iter(()), adapters=sample_adapters(), clock=clock, batch_size=0)


@pytest.mark.unit
def test_a_partially_consumed_stream_says_so_in_its_warnings(clock: Clock) -> None:
    """Counts for a prefix must not be readable as counts for the run."""
    stream = iter_annotated(iter(sample_records()), adapters=sample_adapters(), clock=clock)
    iterator = iter(stream)
    next(iterator)

    assert stream.drained is False
    assert any("not drained" in warning for warning in stream.warnings())


@pytest.mark.unit
def test_a_duplicate_variant_id_across_a_batch_boundary_is_counted_once(clock: Clock) -> None:
    """Coverage denominators are DISTINCT ids, batched or not.

    A coordinate-sorted stream can only repeat an id adjacently, so the carry of
    the previous batch's last id makes per-batch deduplication exactly the global
    deduplication the whole-callset path performs.
    """
    duplicated = record("chr1", 100)
    records = [duplicated, duplicated, duplicated, record("chr2", 200)]
    adapters = AdapterSet(
        consequence=FakeConsequenceAdapter({duplicated.variant_id: (consequence(),)}),
        frequency=FakeFrequencyAdapter({}),
        clinical=None,
    )
    batched = iter_annotated(iter(records), adapters=adapters, clock=clock, batch_size=2)
    list(batched)

    whole = annotate_variants(
        records,
        adapters=AdapterSet(
            consequence=FakeConsequenceAdapter({duplicated.variant_id: (consequence(),)}),
            frequency=FakeFrequencyAdapter({}),
            clinical=None,
        ),
        clock=clock,
    )
    assert batched.coverage() == whole.coverage
    assert batched.coverage()["fake-consequence"] == pytest.approx(0.5)


@pytest.mark.unit
def test_streaming_annotation_is_identical_across_repeat_runs(clock: Clock) -> None:
    """GP-30 at the stage level."""
    digests: list[str] = []
    for _ in range(2):
        stream = iter_annotated(
            iter(sample_records()), adapters=sample_adapters(), clock=clock, batch_size=3
        )
        payload = [
            [a.variant.model_dump(mode="json"), [e.model_dump(mode="json") for e in a.evidence]]
            for a in stream
        ]
        digests.append(stable_hash([payload, stream.coverage(), list(stream.warnings())]))
    assert digests[0] == digests[1]


# ---------------------------------------------------------------------------
# Blocker 2: the spilling ledger
# ---------------------------------------------------------------------------


def in_memory_items(
    corpus: Sequence[EvidenceItem], *, run_id: str = "RUN-1"
) -> tuple[EvidenceItem, ...]:
    """What an ordinary in-memory ledger holds for this corpus, in its order."""
    ledger = EvidenceLedger(run_id=run_id)
    ledger.extend(corpus)
    return ledger.items()


def in_memory_order(corpus: Sequence[EvidenceItem], *, run_id: str = "RUN-1") -> list[str]:
    """The order an ordinary in-memory ledger puts these items in.

    Used as the expectation rather than re-implementing ``_sort_key`` here: the
    claim under test is that the spill agrees with the ledger, and a test that
    re-derived the key could agree with itself while both drifted from the code.
    """
    ledger = EvidenceLedger(run_id=run_id)
    ledger.extend(corpus)
    return [item.evidence_id for item in ledger.items()]


def ledger_corpus(count: int = 40) -> list[EvidenceItem]:
    """Items whose sort keys interleave, so ordering is actually exercised."""
    directions = (
        EvidenceDirection.SUPPORTS,
        EvidenceDirection.CONTRADICTS,
        EvidenceDirection.NEUTRAL,
    )
    categories = (
        EvidenceCategory.ANALYTICAL,
        EvidenceCategory.POPULATION,
        EvidenceCategory.CONSEQUENCE,
    )
    kinds = ("variant", "pair", "gene")
    return [
        evidence_item(
            f"SUBJ-{(i * 7) % count:03d}",
            claim=f"Claim number {i} about subject {(i * 7) % count}.",
            direction=directions[i % 3],
            category=categories[i % 3],
            subject_kind=kinds[i % 3],
        )
        for i in range(count)
    ]


@pytest.mark.unit
def test_a_spilled_ledger_is_observationally_identical_to_an_in_memory_one(
    tmp_path: Path,
) -> None:
    """The whole claim of the spill: nothing about the ledger's behaviour changes."""
    corpus = ledger_corpus()

    memory = EvidenceLedger(run_id="RUN-1")
    memory.extend(corpus)

    spilled = EvidenceLedger(run_id="RUN-1", spill_dir=tmp_path, spill_threshold=5)
    spilled.extend(corpus)

    assert spilled.spilled is True
    assert memory.spilled is False
    assert len(spilled) == len(memory)
    assert spilled.items() == memory.items()
    assert list(spilled.iter_items()) == list(memory.iter_items())
    assert spilled.contradictions() == memory.contradictions()
    for subject in {item.subject_id for item in corpus}:
        assert spilled.for_subject(subject) == memory.for_subject(subject)
    for item in corpus:
        assert spilled.get(item.evidence_id) == memory.get(item.evidence_id)
        assert item.evidence_id in spilled
    assert "EV-NOPE-0000000000000000" not in spilled
    assert spilled.get("EV-NOPE-0000000000000000") is None

    spilled.close()


@pytest.mark.unit
def test_the_spill_orders_exactly_as_python_sorts_including_non_ascii(tmp_path: Path) -> None:
    """SQLite's BINARY collation vs Python's code-point order, asserted not assumed.

    UTF-8 is order-preserving with respect to code points, so a byte-wise SQL sort
    and a Python string sort agree — but the ledger's total order (GP-30) rests on
    that being true, and 'well-known' is not the same as 'checked'.
    """
    subjects = ["a", "A", "z", "Z", "é", "e", "ß", "ss", "中", "🧬", "a\x7f", "a "]
    corpus = [evidence_item(subject, claim=f"Claim about {subject}.") for subject in subjects]

    spilled = EvidenceLedger(run_id="RUN-1", spill_dir=tmp_path, spill_threshold=1)
    spilled.extend(corpus)
    from_disk = [item.subject_id for item in spilled.iter_items()]
    spilled.close()

    memory = EvidenceLedger(run_id="RUN-1")
    memory.extend(corpus)
    assert from_disk == [item.subject_id for item in memory.items()]


@pytest.mark.unit
def test_the_spill_round_trips_an_item_byte_for_byte(tmp_path: Path) -> None:
    """Encode/decode fidelity, because a lossy round trip breaks GP-30 silently."""
    original = evidence_item("SUBJ-1", claim="A claim with a float and a null in its payload.")
    assert decode(encode(original)) == original
    # The stored form is Pydantic's serialiser, not canonical JSON — the spill is
    # scratch and is never hashed. What must survive is the canonical form of the
    # ROUND-TRIPPED item, because that is what reaches the evidence store.
    assert canonical_json(decode(encode(original)).model_dump(mode="json")) == canonical_json(
        original.model_dump(mode="json")
    )

    spilled = EvidenceLedger(run_id="RUN-1", spill_dir=tmp_path, spill_threshold=1)
    stored = spilled.add(original)
    (recovered,) = spilled.items()
    assert recovered == stored.model_copy(update={"run_id": "RUN-1"})
    spilled.close()


@pytest.mark.unit
def test_a_spilled_ledger_still_refuses_a_reused_evidence_id(tmp_path: Path) -> None:
    """Integrity survives the move to disk; only the moment of detection moves."""
    original = evidence_item("SUBJ-1")
    forged = original.model_copy(update={"claim": "A different claim under the same ID."})

    spilled = EvidenceLedger(run_id="RUN-1", spill_dir=tmp_path, spill_threshold=1)
    spilled.add(original)
    with pytest.raises(EvidenceError, match="reused for a different claim"):
        spilled.add(forged)
        spilled.items()  # forces the flush that detects a cross-batch collision
    spilled.close()


@pytest.mark.unit
def test_re_adding_an_identical_item_to_a_spilled_ledger_is_a_no_op(tmp_path: Path) -> None:
    item = evidence_item("SUBJ-1")
    spilled = EvidenceLedger(run_id="RUN-1", spill_dir=tmp_path, spill_threshold=1)
    spilled.add(item)
    spilled.add(item)
    spilled.add(item)
    assert len(spilled) == 1
    spilled.close()


@pytest.mark.unit
def test_insertion_order_does_not_change_a_spilled_ledger(tmp_path: Path) -> None:
    """GP-30: two runs that produce the same evidence in a different order agree."""
    corpus = ledger_corpus(20)

    forward = EvidenceLedger(run_id="RUN-1", spill_dir=tmp_path / "a", spill_threshold=3)
    forward.extend(corpus)
    forward_digest = stable_hash([item.model_dump(mode="json") for item in forward.iter_items()])
    forward.close()

    backward = EvidenceLedger(run_id="RUN-1", spill_dir=tmp_path / "b", spill_threshold=3)
    backward.extend(reversed(corpus))
    backward_digest = stable_hash([item.model_dump(mode="json") for item in backward.iter_items()])
    backward.close()

    assert forward_digest == backward_digest


@pytest.mark.unit
def test_closing_a_spilled_ledger_removes_the_file(tmp_path: Path) -> None:
    """The spill carries proband coordinates. It does not outlive the run (GP-40)."""
    spilled = EvidenceLedger(run_id="RUN-1", spill_dir=tmp_path, spill_threshold=1)
    spilled.extend(ledger_corpus(5))
    path = tmp_path / "evidence-ledger-RUN-1.sqlite"
    assert path.is_file()

    spilled.close()
    assert not path.exists()
    spilled.close()  # idempotent


@pytest.mark.unit
def test_a_ledger_without_a_spill_dir_never_touches_the_filesystem(tmp_path: Path) -> None:
    """The default is unchanged behaviour, and that has to be structural."""
    ledger = EvidenceLedger(run_id="RUN-1", spill_threshold=1)
    ledger.extend(ledger_corpus(50))
    assert ledger.spilled is False
    assert list(tmp_path.iterdir()) == []


@pytest.mark.unit
def test_a_stale_spill_file_is_replaced_not_appended_to(tmp_path: Path) -> None:
    """A re-run of the same case derives the same run id, so the same filename."""
    first = EvidenceLedger(run_id="RUN-1", spill_dir=tmp_path, spill_threshold=1)
    first.extend(ledger_corpus(6))
    path = tmp_path / "evidence-ledger-RUN-1.sqlite"
    # Simulate a crashed run: the file is left behind rather than closed.
    first._spill._connection.close()  # type: ignore[union-attr]
    assert path.is_file()

    second = EvidenceLedger(run_id="RUN-1", spill_dir=tmp_path, spill_threshold=1)
    second.extend(ledger_corpus(3))
    assert len(second) == 3
    second.close()


@pytest.mark.unit
def test_the_store_writes_an_evidence_stream_exactly_as_it_writes_a_sequence(
    tmp_path: Path,
) -> None:
    """Batching the write must not change a byte of the Parquet export (GP-30)."""
    corpus = [item.model_copy(update={"run_id": "RUN-1"}) for item in ledger_corpus(30)]

    with EvidenceStore(tmp_path / "batch.duckdb") as store:
        store.initialise()
        assert store.write_evidence(corpus) == len(corpus)
        batch_files = store.export_parquet(tmp_path / "batch-parquet")
        batch_counts = store.counts()

    with EvidenceStore(tmp_path / "stream.duckdb") as store:
        store.initialise()
        assert store.write_evidence_stream(iter(corpus), batch_size=4) == len(corpus)
        stream_files = store.export_parquet(tmp_path / "stream-parquet")
        stream_counts = store.counts()

    assert batch_counts == stream_counts
    for table, path in batch_files.items():
        assert path.read_bytes() == stream_files[table].read_bytes(), table


@pytest.mark.unit
def test_a_spilled_ledger_persists_through_the_store_unchanged(tmp_path: Path) -> None:
    """The composition the orchestrator will use: spill in, stream out, same rows."""
    corpus = ledger_corpus(25)

    memory = EvidenceLedger(run_id="RUN-1")
    memory.extend(corpus)
    with EvidenceStore(tmp_path / "a.duckdb") as store:
        store.initialise()
        store.write_evidence(memory.items())
        expected = store.export_parquet(tmp_path / "a-parquet")["evidence_items"].read_bytes()

    spilled = EvidenceLedger(run_id="RUN-1", spill_dir=tmp_path, spill_threshold=4)
    spilled.extend(corpus)
    with EvidenceStore(tmp_path / "b.duckdb") as store:
        store.initialise()
        store.write_evidence_stream(spilled.iter_items())
        observed = store.export_parquet(tmp_path / "b-parquet")["evidence_items"].read_bytes()
    spilled.close()

    assert observed == expected


@pytest.mark.unit
def test_the_spill_can_be_driven_directly(tmp_path: Path) -> None:
    """The backend is usable on its own, so it can be tested without a ledger."""
    spill = SqliteEvidenceSpill(tmp_path / "direct.sqlite", flush_batch=2)
    corpus = ledger_corpus(9)
    spill.extend(corpus)
    assert len(spill) == 9
    assert [item.evidence_id for item in spill.iter_items()] == in_memory_order(corpus)
    spill.close()


# ---------------------------------------------------------------------------
# Blocker 3: selection
# ---------------------------------------------------------------------------

THRESHOLDS = SelectionThresholds()
FREQUENCY = FrequencyThresholds()


def decide(
    variant: VariantRecord, clock: Clock, **overrides: object
) -> tuple[bool, str, tuple[str, ...]]:
    thresholds = SelectionThresholds(**overrides)  # type: ignore[arg-type]
    stream = iter_selected([variant], frequency=FREQUENCY, thresholds=thresholds, clock=clock)
    (decision,) = list(stream.decisions())
    return decision.retained, decision.reason, decision.notes


@pytest.mark.unit
def test_a_variant_with_no_frequency_data_is_retained(clock: Clock) -> None:
    """GP-14, and the failure mode that loses the whole challenge.

    A variant no reference cohort has recorded is UNKNOWN. Treating it as AF = 0
    would be the strongest possible rarity signal; treating it as common enough to
    delete is the same error with the opposite sign. It is retained, and the gap is
    recorded as a note so the report can state it.
    """
    novel = record(consequences=[consequence()])
    assert novel.has_frequency_data is False

    retained, reason, notes = decide(novel, clock)
    assert retained is True
    assert reason == KEEP_CODING_OR_SPLICE
    assert NOTE_FREQUENCY_UNKNOWN in notes


@pytest.mark.unit
def test_an_underpowered_frequency_is_unknown_not_rare_and_not_common(clock: Clock) -> None:
    """AC=1 in AN=40 is an AF of 0.025 and is evidence of nothing (ADR 0010)."""
    tiny_cohort = record(
        consequences=[consequence()], frequencies=[frequency(0.025, allele_number=40)]
    )
    retained, reason, notes = decide(tiny_cohort, clock)
    assert retained is True
    assert reason == KEEP_CODING_OR_SPLICE
    assert NOTE_FREQUENCY_UNDERPOWERED in notes


@pytest.mark.unit
def test_a_not_assessed_impact_is_never_dropped_as_benign(clock: Clock) -> None:
    """ADR 0016. `impact is None` means NOT ASSESSED, emphatically not MODIFIER.

    A MANE interval join produces exactly this shape for every variant it places:
    a gene, a transcript, and no computed consequence. Reading it as a prediction
    of negligible effect would delete every gene-assignment-only candidate.
    """
    unassessed = record(consequences=[consequence(terms=("transcript_variant",), impact=None)])
    retained, reason, _ = decide(unassessed, clock)
    assert retained is True
    assert reason == KEEP_IMPACT_NOT_ASSESSED


@pytest.mark.unit
def test_an_explicit_modifier_prediction_is_dropped(clock: Clock) -> None:
    """The other half of ADR 0016: MODIFIER *is* a positive prediction, and counts."""
    modifier = record(
        consequences=[consequence(terms=("intron_variant",), impact=ImpactSeverity.MODIFIER)],
        frequencies=[frequency(1e-6)],
    )
    retained, reason, _ = decide(modifier, clock)
    assert retained is False
    assert reason == DROP_NOT_CODING_OR_SPLICE


@pytest.mark.unit
def test_a_common_variant_is_dropped_on_a_measured_number(clock: Clock) -> None:
    common = record(consequences=[consequence()], frequencies=[frequency(0.3)])
    retained, reason, _ = decide(common, clock)
    assert retained is False
    assert reason == DROP_COMMON_IN_POPULATION


@pytest.mark.unit
def test_a_variant_between_the_selection_and_ranking_cut_points_survives(clock: Clock) -> None:
    """Nothing this stage deletes was still rankable above zero on rarity."""
    borderline = record(
        consequences=[consequence()],
        frequencies=[frequency((FREQUENCY.max_plausible_recessive + 0.02) / 2)],
    )
    retained, _, _ = decide(borderline, clock)
    assert retained is True


@pytest.mark.unit
def test_a_splice_region_synonymous_variant_is_retained(clock: Clock) -> None:
    """The exact case `prioritization.filters` warns a consequence filter would lose.

    Both the joined SnpEff spelling and separate terms have to work: matching a
    combined `a&b` string against a set of atomic terms silently fails, and for
    this stage 'silently fails' means deleting a splice candidate.
    """
    joined = record(
        consequences=[
            consequence(
                terms=("splice_region_variant&synonymous_variant",), impact=ImpactSeverity.LOW
            )
        ],
        frequencies=[frequency(1e-6)],
    )
    separate = record(
        consequences=[
            consequence(
                terms=("synonymous_variant", "splice_region_variant"), impact=ImpactSeverity.LOW
            )
        ],
        frequencies=[frequency(1e-6)],
    )
    for variant in (joined, separate):
        retained, reason, _ = decide(variant, clock)
        assert retained is True
        assert reason == KEEP_CODING_OR_SPLICE


@pytest.mark.unit
def test_a_plain_synonymous_variant_is_dropped(clock: Clock) -> None:
    """Deliberate: including it roughly doubles the surviving set for no gain."""
    synonymous = record(
        consequences=[consequence(terms=("synonymous_variant",), impact=ImpactSeverity.LOW)],
        frequencies=[frequency(1e-6)],
    )
    retained, reason, _ = decide(synonymous, clock)
    assert retained is False
    assert reason == DROP_NOT_CODING_OR_SPLICE


@pytest.mark.unit
def test_any_transcript_can_rescue_a_variant(clock: Clock) -> None:
    """GP-20 / ASSUMPTION-TRANSCRIPT-01: never collapse to the canonical transcript."""
    mixed = record(
        consequences=[
            consequence(terms=("intron_variant",), impact=ImpactSeverity.MODIFIER, transcript="T1"),
            consequence(terms=("stop_gained",), impact=ImpactSeverity.HIGH, transcript="T2"),
        ],
        frequencies=[frequency(1e-6)],
    )
    retained, reason, _ = decide(mixed, clock)
    assert retained is True
    assert reason == KEEP_CODING_OR_SPLICE


@pytest.mark.unit
def test_a_high_impact_with_an_unrecognised_term_is_retained(clock: Clock) -> None:
    """Fail open on vocabulary: a new tool adds terms faster than impact classes."""
    exotic = record(
        consequences=[
            consequence(terms=("some_future_so_term",), impact=ImpactSeverity.HIGH),
        ],
        frequencies=[frequency(1e-6)],
    )
    retained, reason, _ = decide(exotic, clock)
    assert retained is True
    assert reason == KEEP_CODING_OR_SPLICE


@pytest.mark.unit
def test_a_splice_ai_delta_retains_a_variant_with_no_splice_term(clock: Clock) -> None:
    predicted = record(
        consequences=[
            consequence(terms=("intron_variant",), impact=ImpactSeverity.MODIFIER, splice_ai=0.7)
        ],
        frequencies=[frequency(1e-6)],
    )
    retained, reason, _ = decide(predicted, clock)
    assert retained is True
    assert reason == KEEP_CODING_OR_SPLICE


@pytest.mark.unit
def test_a_curated_pathogenic_assertion_overrides_both_gates(clock: Clock) -> None:
    """A known pathogenic allele that is common in some cohort is still the answer."""
    known = record(
        consequences=[consequence(terms=("intron_variant",), impact=ImpactSeverity.MODIFIER)],
        frequencies=[frequency(0.4)],
        assertions=[
            ClinicalAssertion(source="ClinVar", version="v0", significance="Likely_pathogenic")
        ],
    )
    retained, reason, _ = decide(known, clock)
    assert retained is True
    assert reason == KEEP_CLINICAL_ASSERTION


@pytest.mark.unit
def test_a_conflicting_clinvar_classification_does_not_rescue(clock: Clock) -> None:
    """`Conflicting_classifications_of_pathogenicity` contains the substring and is
    an explicit statement that curators disagree."""
    conflicting = record(
        consequences=[consequence(terms=("intron_variant",), impact=ImpactSeverity.MODIFIER)],
        frequencies=[frequency(0.4)],
        assertions=[
            ClinicalAssertion(
                source="ClinVar",
                version="v0",
                significance="Conflicting_classifications_of_pathogenicity",
            )
        ],
    )
    retained, reason, _ = decide(conflicting, clock)
    assert retained is False
    assert reason == DROP_COMMON_IN_POPULATION


@pytest.mark.unit
def test_a_variant_with_no_gene_assignment_is_dropped_and_counted(clock: Clock) -> None:
    """Pairing already ignores these. Counting them makes the loss visible."""
    intergenic = record(frequencies=[frequency(1e-6)])
    retained, reason, _ = decide(intergenic, clock)
    assert retained is False
    assert reason == DROP_NO_GENE_ASSIGNMENT

    retained, _, _ = decide(intergenic, clock, drop_without_gene_assignment=False)
    assert retained is True


@pytest.mark.unit
def test_selection_never_mutates_a_record(clock: Clock) -> None:
    """A selected variant is the same object, so no downstream artifact changes."""
    original = record(consequences=[consequence()], frequencies=[frequency(1e-6)])
    result = select_variants([original], frequency=FREQUENCY, clock=clock)
    (selected,) = result.variants
    assert selected is original
    assert selected.qc_flags == ()


@pytest.mark.unit
def test_the_counts_partition_the_input(clock: Clock) -> None:
    """Arithmetic a reader can check: keeps plus drops equal what went in."""
    variants = [
        record("chr1", 1, consequences=[consequence()]),
        record("chr1", 2, consequences=[consequence()], frequencies=[frequency(0.9)]),
        record("chr1", 3, frequencies=[frequency(1e-6)]),
        record("chr1", 4, consequences=[consequence(impact=None)]),
        record(
            "chr1",
            5,
            consequences=[consequence(terms=("intron_variant",), impact=ImpactSeverity.MODIFIER)],
            frequencies=[frequency(1e-6)],
        ),
    ]
    result = select_variants(variants, frequency=FREQUENCY, clock=clock)
    report = result.report

    assert report.input_count == len(variants)
    assert sum(report.retained_by_reason.values()) == report.retained_count
    assert sum(report.dropped_by_reason.values()) == report.dropped_count
    assert report.retained_count + report.dropped_count == report.input_count
    assert list(report.retained_by_reason) == list(KEEP_REASONS)
    assert list(report.dropped_by_reason) == list(DROP_REASONS)
    assert len(result.decisions) == len(variants)


@pytest.mark.unit
def test_the_streaming_and_materialising_selection_paths_agree(clock: Clock) -> None:
    variants = [
        record(
            "chr1",
            i,
            consequences=[consequence()] if i % 2 else [],
            frequencies=[frequency(i / 100)],
        )
        for i in range(1, 30)
    ]
    materialised = select_variants(variants, frequency=FREQUENCY, clock=clock)

    stream = iter_selected(iter(variants), frequency=FREQUENCY, clock=clock)
    streamed = list(stream)

    assert streamed == list(materialised.variants)
    assert stream.report().as_payload() == materialised.report.as_payload()
    assert stream.evidence() == materialised.evidence


@pytest.mark.unit
def test_disabling_selection_retains_everything_and_says_so(clock: Clock) -> None:
    variants = [record("chr1", i, frequencies=[frequency(0.9)]) for i in range(1, 6)]
    result = select_variants(
        variants, frequency=FREQUENCY, thresholds=SelectionThresholds(enabled=False), clock=clock
    )
    assert len(result.variants) == len(variants)
    assert any("DISABLED" in warning for warning in result.report.warnings)


@pytest.mark.unit
def test_dropping_unknown_frequencies_is_possible_and_loudly_warned(clock: Clock) -> None:
    """The switch exists so the policy is arguable, not so it is used (ADR 0019)."""
    novel = record(consequences=[consequence()])
    result = select_variants(
        [novel],
        frequency=FREQUENCY,
        thresholds=SelectionThresholds(retain_unknown_frequency=False),
        clock=clock,
    )
    assert result.variants == ()
    assert any("GP-14" in warning for warning in result.report.warnings)
    assert any("retain_unknown_frequency is FALSE" in w for w in result.report.warnings)


@pytest.mark.unit
def test_a_selection_cut_point_tighter_than_the_ranking_cut_point_is_reported(
    clock: Clock,
) -> None:
    """Collapsing the two would delete candidates the ranker could still score."""
    result = select_variants(
        [record(consequences=[consequence()], frequencies=[frequency(1e-6)])],
        frequency=FREQUENCY,
        thresholds=SelectionThresholds(max_population_frequency=FREQUENCY.max_plausible_recessive),
        clock=clock,
    )
    assert any("not looser than the ranking cut-point" in w for w in result.report.warnings)


@pytest.mark.unit
def test_selection_evidence_is_aggregate_deterministic_and_limited(clock: Clock) -> None:
    """GP-10 and GP-17: the counts are citable, and the item states what it is not."""
    variants = [record("chr1", i, consequences=[consequence()]) for i in range(1, 6)]
    first = select_variants(variants, frequency=FREQUENCY, clock=clock)
    second = select_variants(variants, frequency=FREQUENCY, clock=clock)

    assert first.evidence == second.evidence
    assert first.evidence, "the stage must be able to justify what it removed"
    assert len(first.evidence) < len(variants), "aggregate, not one item per variant"
    for item in first.evidence:
        assert item.limitations
        assert item.direction is EvidenceDirection.NEUTRAL
        assert item.category is EvidenceCategory.PROVENANCE
        assert item.subject_kind == "run"


@pytest.mark.unit
def test_a_selection_report_carries_no_variant_identifier(clock: Clock) -> None:
    """It is classified DERIVED_SAFE, so it must be counts only (GP-41)."""
    variants = [record("chr7", 123_456, consequences=[consequence()]) for _ in range(1)]
    result = select_variants(variants, frequency=FREQUENCY, clock=clock)
    rendered = canonical_json(result.report.as_payload())
    assert variants[0].variant_id not in rendered
    assert "123456" not in rendered
    assert "chr7" not in rendered


@pytest.mark.unit
def test_an_impossible_selection_threshold_is_refused() -> None:
    """The bounds are the config model's now, so the message is Pydantic's.

    `SelectionThresholds` moved from a dataclass in this module to
    `mva.config.SelectionThresholds` when the field landed on `CaseConfig`
    (docs/handoff-scale.md §2, ADR 0019). `ValidationError` is a `ValueError`, so
    what changes is only the wording: the bound is now declared on the field and
    is therefore visible in the config schema rather than only in a `__post_init__`.
    Zero remains refused for the same reason — it would delete every variant any
    cohort has ever seen once.
    """
    with pytest.raises(ValueError, match="greater than 0"):
        SelectionThresholds(max_population_frequency=0.0)
    with pytest.raises(ValueError, match="less than or equal to 1"):
        SelectionThresholds(min_splice_ai_delta=1.5)


@pytest.mark.unit
def test_re_iterating_a_selection_stream_raises(clock: Clock) -> None:
    stream = iter_selected([record()], frequency=FREQUENCY, clock=clock)
    list(stream)
    with pytest.raises(ValueError, match="already been iterated"):
        list(stream)


# ---------------------------------------------------------------------------
# GP-30 across the whole chain, in separate processes
# ---------------------------------------------------------------------------


@pytest.mark.integration
@pytest.mark.slow
def test_the_whole_streamed_chain_is_byte_identical_under_two_hash_seeds(
    tmp_path: Path,
) -> None:
    """Annotate -> ledger -> select, digested end to end, under two PYTHONHASHSEEDs.

    In-process determinism tests cannot see hash-seed dependence, because the seed
    is fixed for the life of the interpreter. Set iteration order, dict ordering
    from a `frozenset` and any accidental reliance on `hash()` only diverge across
    processes (GP-30).
    """
    digests: list[dict[str, object]] = []
    for seed in ("0", "987654321"):
        proc = subprocess.run(  # noqa: S603 - fixed argv, no shell
            [
                sys.executable,
                "-m",
                "tools.scale.stage_harness",
                "digest",
                "--records",
                "5000",
                "--tmp",
                str(tmp_path / seed),
                "--rss-cap",
                "0",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
            env={
                "PYTHONHASHSEED": seed,
                "PATH": "/usr/bin:/bin",
                "HOME": str(tmp_path),
                "MVA_SCALE_TMP": str(tmp_path / seed),
            },
            timeout=600,
        )
        # The harness prints one INDENTED object, so parse all of stdout: its last
        # line is "}".
        digests.append(json.loads(proc.stdout.strip())["notes"])

    first, second = digests
    assert first["hash_seed"] != second["hash_seed"]
    for key in (
        "variants_sha256",
        "verdicts_sha256",
        "ledger_sha256",
        "report_sha256",
        "annotation_warnings_sha256",
    ):
        assert first[key] == second[key], key


# ---------------------------------------------------------------------------
# The two stream shims the composition needs: hard filter, soft flags, artifact
# ---------------------------------------------------------------------------


def filterable_records() -> list[VariantRecord]:
    """A set covering every hard-filter reason plus records that survive."""
    return [
        record("chr1", 100, consequences=[consequence()], frequencies=[frequency(1e-6)]),
        record("chr1", 200, zygosity=Zygosity.HOM_REF),
        record("chr1", 300, alt="*"),
        record("chr1", 400, zygosity=Zygosity.UNKNOWN),
        record("chr2", 500, consequences=[consequence()], frequencies=[frequency(0.4)]),
        record("chr2", 600, zygosity=Zygosity.HOM_ALT, frequencies=[frequency(1e-6)]),
    ]


@pytest.mark.unit
def test_the_streaming_hard_filter_matches_the_batch_one() -> None:
    """Same records, same counts — only the memory profile differs.

    ``retained``, not ``flagged``: the composition root consumes ``.retained`` and
    soft-flags it, and the two differ in flag ORDER even where they agree on the
    flag set (see ``HardFilterStream``). Comparing against ``.flagged`` here would
    lock in a silent change to every downstream artifact's bytes.
    """
    records = filterable_records()
    batch = apply_hard_filters(records, expected_build=GenomeBuild.GRCH38)

    stream = iter_hard_filtered(iter(records), expected_build=GenomeBuild.GRCH38)
    streamed = list(stream)

    assert tuple(streamed) == batch.retained
    assert stream.counts() == batch.counts


@pytest.mark.unit
def test_soft_flagging_retained_and_flagged_records_differs_in_order() -> None:
    """The discrepancy the streaming shim had to be careful about, pinned.

    ``apply_hard_filters`` attaches configuration-free flags to ``.flagged``;
    ``with_qc_flags`` appends rather than re-orders; so soft-flagging an
    already-flagged record yields the same flag SET in a different ORDER. This is
    not a bug being fixed here — it is the reason ``iter_hard_filtered`` yields
    ``.retained`` — but an undocumented difference like this is how a refactor
    changes results, so it gets a test rather than a comment.
    """
    # A homozygous call that is also common: `homozygous_call` comes from the
    # configuration-free pass and is LAST in SOFT_FLAGS, while `common_variant`
    # comes from the threshold pass and is FIRST. Flag them in either order and
    # the set matches; the tuple does not.
    records = [
        record(
            "chr2",
            600,
            zygosity=Zygosity.HOM_ALT,
            consequences=[consequence()],
            frequencies=[frequency(0.4)],
        )
    ]
    batch = apply_hard_filters(records, expected_build=GenomeBuild.GRCH38)

    from_retained = apply_soft_flags(
        batch.retained, frequency=FREQUENCY, quality=QualityThresholds()
    )
    from_flagged = apply_soft_flags(batch.flagged, frequency=FREQUENCY, quality=QualityThresholds())

    assert [set(v.qc_flags) for v in from_retained] == [set(v.qc_flags) for v in from_flagged]
    assert from_retained[0].qc_flags != from_flagged[0].qc_flags, (
        "if these ever become equal the ordering caveat in HardFilterStream is stale"
    )


@pytest.mark.unit
def test_the_streaming_hard_filter_keeps_no_per_variant_removal_list() -> None:
    """Counts, not coordinates (GP-41), and not millions of tuples either."""
    stream = iter_hard_filtered(iter(filterable_records()), expected_build=GenomeBuild.GRCH38)
    list(stream)
    rendered = canonical_json(stream.counts())
    assert "chr1" not in rendered
    assert stream.counts()["removed"] == 3


@pytest.mark.unit
def test_the_streaming_soft_flags_match_the_batch_ones() -> None:
    records = [r for r in filterable_records() if r.genotype.carries_alt]
    batch = apply_soft_flags(records, frequency=FREQUENCY, quality=QualityThresholds())
    streamed = tuple(
        iter_soft_flagged(iter(records), frequency=FREQUENCY, quality=QualityThresholds())
    )
    assert streamed == batch


@pytest.mark.unit
def test_re_iterating_a_hard_filter_stream_raises() -> None:
    stream = iter_hard_filtered(iter(filterable_records()), expected_build=GenomeBuild.GRCH38)
    list(stream)
    with pytest.raises(ValueError, match="already been iterated"):
        list(stream)


@pytest.mark.unit
def test_the_pushed_json_rows_artifact_is_byte_identical_to_the_pulled_one(
    tmp_path: Path,
) -> None:
    """A mid-chain artifact must not be a different artifact (GP-30).

    ``write_json_rows_artifact`` pulls, ``open_json_rows_artifact`` pushes. They
    have to produce the same bytes or the choice between them becomes a choice
    about output.
    """
    rows: list[dict[str, object]] = [
        {"b": 2, "a": 1, "nested": {"z": None, "y": [1, 2]}},
        {"b": 3, "a": 2, "nested": {"z": "x", "y": []}},
    ]
    pulled = tmp_path / "pulled.json"
    write_canonical_json_rows(pulled, rows)

    pushed = tmp_path / "pushed.json"
    with pushed.open("w", encoding="utf-8") as handle:
        sink = JsonRowsSink(handle)
        for row in rows:
            sink.write(row)
        sink.close()

    assert pushed.read_bytes() == pulled.read_bytes()
    assert pulled.read_text(encoding="utf-8") == canonical_json(rows) + "\n"


@pytest.mark.unit
def test_an_empty_pushed_artifact_is_still_a_valid_empty_array(tmp_path: Path) -> None:
    path = tmp_path / "empty.json"
    with path.open("w", encoding="utf-8") as handle:
        sink = JsonRowsSink(handle)
        sink.close()
    assert path.read_text(encoding="utf-8") == "[]\n"
    assert json.loads(path.read_text(encoding="utf-8")) == []


@pytest.mark.unit
def test_the_spill_compresses_against_a_dictionary_it_stores(tmp_path: Path) -> None:
    """The size fix, locked as behaviour rather than left as a comment.

    ``limitations`` is 38% of an evidence document and is byte-identical across
    millions of items. Stored verbatim that measured 4,933 bytes per item on disk;
    deflated against 4 KiB of its own boilerplate it measures ~236. The dictionary
    lives in the file, so a row is always readable by whatever wrote it.
    """
    # Stamped up front: SqliteEvidenceSpill stores what it is given, while
    # EvidenceLedger stamps run_id on the way in. Comparing the two needs both
    # sides to have been stamped the same way.
    corpus = [item.model_copy(update={"run_id": "RUN-1"}) for item in ledger_corpus(200)]
    path = tmp_path / "compressed.sqlite"
    spill = SqliteEvidenceSpill(path, flush_batch=100)
    spill.extend(corpus)
    stored = list(spill.iter_items())

    connection = sqlite3.connect(path)
    (dictionary,) = connection.execute("SELECT value FROM meta WHERE key='dictionary'").fetchone()
    blobs = [bytes(row[0]) for row in connection.execute("SELECT document FROM ledger")]
    connection.close()
    spill.close()

    assert 0 < len(dictionary) <= DICTIONARY_BYTES
    assert tuple(stored) == in_memory_items(corpus)

    plain = sum(len(encode(item).encode()) for item in corpus)
    packed = sum(len(blob) for blob in blobs)
    assert packed * 2 < plain, f"compression achieved only {plain / packed:.2f}x"
    # And the raw bytes really are compressed, not merely re-encoded text.
    assert decompress(blobs[0], dictionary).startswith("{")


@pytest.mark.unit
def test_compression_round_trips_every_document_exactly() -> None:
    """A lossy round trip here would corrupt evidence, silently, at scale only."""
    corpus = ledger_corpus(40)
    documents = [encode(item) for item in corpus]
    dictionary = build_dictionary(documents)
    for document in documents:
        assert decompress(compress(document, dictionary), dictionary) == document


@pytest.mark.unit
def test_a_closed_ledger_still_reports_its_count_but_refuses_its_contents(
    tmp_path: Path,
) -> None:
    """The composition root reports `evidence_count` after persisting and closing.

    Making that number depend on whether someone had already closed the ledger
    would be a trap. Reading the *contents* after close must raise instead, because
    an empty ledger would read as a run that produced no evidence.
    """
    corpus = ledger_corpus(12)
    ledger = EvidenceLedger(run_id="RUN-1", spill_dir=tmp_path, spill_threshold=3)
    ledger.extend(corpus)
    count = len(ledger)
    ledger.close()

    assert len(ledger) == count == len(corpus)
    with pytest.raises(EvidenceError, match="was closed"):
        ledger.items()
    with pytest.raises(EvidenceError, match="was closed"):
        ledger.get(corpus[0].evidence_id)


@pytest.mark.unit
def test_closing_an_unspilled_ledger_changes_nothing(tmp_path: Path) -> None:
    ledger = EvidenceLedger(run_id="RUN-1")
    ledger.extend(ledger_corpus(4))
    ledger.close()
    assert len(ledger) == 4
    assert len(ledger.items()) == 4


@pytest.mark.unit
def test_the_ledger_works_as_a_context_manager(tmp_path: Path) -> None:
    with EvidenceLedger(run_id="RUN-1", spill_dir=tmp_path, spill_threshold=2) as ledger:
        ledger.extend(ledger_corpus(9))
        assert ledger.spilled is True
    assert not list(tmp_path.glob("*.sqlite"))


@pytest.mark.unit
def test_a_failed_pushed_artifact_registers_no_provenance(tmp_path: Path) -> None:
    """A half-written artifact with provenance claiming it is complete is worse
    than one that is obviously unfinished."""
    context = RunContext(
        config=None,  # type: ignore[arg-type] - unused by the paths under test
        workspace=Workspace(root=tmp_path, repo_root=tmp_path),
        clock=FixedClock(datetime(2026, 1, 1, tzinfo=UTC)),
        run_id="RUN-1",
        run_dir=tmp_path / "runs" / "RUN-1",
        started_at=datetime(2026, 1, 1, tzinfo=UTC),
    )
    context.run_dir.mkdir(parents=True, exist_ok=True)

    with (
        pytest.raises(RuntimeError, match="stage blew up"),
        context.open_json_rows_artifact(
            "variants/annotated.json",
            kind=ArtifactKind.ANNOTATED_VARIANTS,
            stage="annotate",
        ) as sink,
    ):
        sink.write({"a": 1})
        msg = "the stage blew up"
        raise RuntimeError(msg)

    assert context.artifacts == []
    assert sink.provenance is None


@pytest.mark.unit
def test_a_pushed_artifact_hands_back_its_provenance(tmp_path: Path) -> None:
    context = RunContext(
        config=None,  # type: ignore[arg-type] - unused by the paths under test
        workspace=Workspace(root=tmp_path, repo_root=tmp_path),
        clock=FixedClock(datetime(2026, 1, 1, tzinfo=UTC)),
        run_id="RUN-1",
        run_dir=tmp_path / "runs" / "RUN-1",
        started_at=datetime(2026, 1, 1, tzinfo=UTC),
    )
    context.run_dir.mkdir(parents=True, exist_ok=True)

    with context.open_json_rows_artifact(
        "variants/annotated.json",
        kind=ArtifactKind.ANNOTATED_VARIANTS,
        stage="annotate",
    ) as sink:
        for index in range(3):
            sink.write({"index": index})

    assert sink.provenance is not None
    assert sink.provenance.row_count == 3
    assert context.artifacts == [sink.provenance]


@pytest.mark.unit
def test_reading_an_empty_spill_before_writing_does_not_disable_compression(
    tmp_path: Path,
) -> None:
    """A read before the first write must not pin an empty dictionary."""
    spill = SqliteEvidenceSpill(tmp_path / "empty-first.sqlite", flush_batch=2)
    assert len(spill) == 0
    assert list(spill.iter_items()) == []
    assert spill.get("EV-NOPE-0000000000000000") is None

    corpus = [item.model_copy(update={"run_id": "RUN-1"}) for item in ledger_corpus(80)]
    spill.extend(corpus)
    assert tuple(spill.iter_items()) == in_memory_items(corpus)

    connection = sqlite3.connect(tmp_path / "empty-first.sqlite")
    (dictionary,) = connection.execute("SELECT value FROM meta WHERE key='dictionary'").fetchone()
    connection.close()
    spill.close()
    assert len(dictionary) > 0, "the dictionary was never built"
