"""Annotation orchestration: adapters in, annotated records + evidence out.

This module is the only place that knows *what annotation means* for this
pipeline. The adapters know how to look things up; the service decides what a
lookup result licenses you to say. Three of those decisions carry most of the
weight:

**Absence is recorded, never imputed (GP-14).** A variant the frequency adapter
has never heard of keeps an empty ``population_frequencies`` tuple — so
``VariantRecord.has_frequency_data`` is False — and gains a NEUTRAL evidence item
that says, in words, that frequency data is unavailable and that this is not
evidence of rarity. The tempting alternative, defaulting to AF = 0, would turn a
gap in reference coverage into the strongest possible rarity signal, promoting
exactly those variants the reference cohorts cover worst. That is a fabricated
negative, and it is how under-represented ancestries get spurious "ultra-rare"
candidates.

**Prediction is filed as prediction (GP-12).** Consequence annotations are
``COMPUTATIONAL_PREDICTION`` / ``IN_SILICO_PREDICTION``. Frequencies are
``DATABASE_ASSERTION`` / ``CURATED_DATABASE`` and therefore carry a versioned
``Citation``, which ``EvidenceItem`` enforces rather than trusts.

**All transcripts survive (GP-20 / ASSUMPTION-TRANSCRIPT-01).** Whatever the
adapter returns is attached in full. Nothing here collapses a variant to its
canonical transcript.

Everything is deterministic: no wall clock (the ``Clock`` is injected), no set or
dict iteration in output ordering, and evidence emitted in input-record order with
a fixed per-record sequence.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from datetime import datetime
from itertools import islice

from mva.annotation.base import (
    SYNTHETIC_STANDIN_LIMITATION,
    AdapterDescriptor,
    AdapterSet,
)
from mva.annotation.triage import (
    DROP_REASONS as TRIAGE_DROP_REASONS,
)
from mva.annotation.triage import (
    KEEP_REASONS as TRIAGE_KEEP_REASONS,
)
from mva.annotation.triage import (
    KEEP_TRIAGE_DISABLED,
    TriageReport,
    TriageVerdict,
    VariantTriage,
)
from mva.clock import Clock
from mva.errors import AnnotationError
from mva.models.base import AssertionTier
from mva.models.discard import DiscardSink
from mva.models.evidence import (
    Citation,
    EvidenceCategory,
    EvidenceDirection,
    EvidenceItem,
    EvidenceStrength,
    EvidenceType,
    make_evidence_id,
)
from mva.models.variant import (
    ClinicalAssertion,
    ConsequenceAnnotation,
    ImpactSeverity,
    PopulationFrequency,
    VariantRecord,
)

#: Ordinal impact -> evidence strength. A prediction never rises above MODERATE:
#: "predicted stop-gained" is a strong hint, not a demonstrated loss of function.
_IMPACT_STRENGTH: Mapping[ImpactSeverity, EvidenceStrength] = {
    ImpactSeverity.HIGH: EvidenceStrength.MODERATE,
    ImpactSeverity.MODERATE: EvidenceStrength.SUPPORTING,
    ImpactSeverity.LOW: EvidenceStrength.WEAK,
    ImpactSeverity.MODIFIER: EvidenceStrength.INSUFFICIENT,
}

#: Impacts whose prediction argues *for* a variant mattering. LOW/MODIFIER are
#: neutral rather than contradicting: a synonymous call is uninformative, not
#: exculpatory (it can still disrupt splicing or codon usage).
_SUPPORTING_IMPACTS: frozenset[ImpactSeverity] = frozenset(
    {ImpactSeverity.HIGH, ImpactSeverity.MODERATE}
)

_CONSEQUENCE_LIMITATION = (
    "A predicted molecular consequence, not an observed one: no functional assay, no "
    "RNA evidence and no segregation data support it. Transcript-level predictions "
    "disagree between tools and between transcripts of the same gene, and impact class "
    "is an ordinal bucket, not a measurement of severity."
)

_FREQUENCY_LIMITATION = (
    "An observation from a reference cohort, not an interpretation. Whether this "
    "frequency is compatible with the disease model is decided by the prioritisation "
    "stage against configured thresholds, never here. Frequencies are population- and "
    "release-specific, and cohort composition biases them for under-represented "
    "ancestries."
)

_NO_FREQUENCY_LIMITATION = (
    "This item records a GAP, not a finding. Absence from a reference cohort commonly "
    "reflects low coverage at the site, an ancestry not represented in the panel, an "
    "indel-representation mismatch or a failed liftover — not true rarity. It must "
    "never be scored as allele frequency 0, and it is not evidence that the variant is "
    "rare (GP-14)."
)

_CLINICAL_LIMITATION = (
    "A third-party curated assertion, inheriting that submitter's evidence, review "
    "depth and date. Star ratings compress conflicting submissions, and re-curation "
    "changes significance over time; the assertion is only as current as its release."
)


@dataclass(frozen=True)
class AnnotationResult:
    """Everything one annotation pass produced.

    ``coverage`` and ``warnings`` are first-class results rather than log output:
    "the frequency adapter knew nothing about 30% of these variants" is a fact the
    report must be able to state, and a log line cannot be cited.
    """

    variants: tuple[VariantRecord, ...]
    evidence: tuple[EvidenceItem, ...]
    coverage: dict[str, float] = field(default_factory=dict)
    """Fraction of distinct input variants for which each adapter returned data,
    keyed by adapter name. 0.0 means "asked, told nothing" — not "not asked"."""

    warnings: tuple[str, ...] = ()


#: Records held in memory at once by :func:`iter_annotated`, and therefore the
#: number of IDs one adapter call receives.
#:
#: Adapters are batch lookups by contract (``Sequence[str]`` in,
#: ``Mapping[str, ...]`` out), so the choice is a trade between round trips and
#: resident memory, not between batching and not batching.
#:
#: **Why 250,000 and not 5,000.** The dominant adapter here is SnpEff, and its
#: cost is ~100% fixed start-up: measured on this machine, one variant takes
#: 35.5 s and 5,000 take 31.6 s, because both are a GRCh38.115 database load and
#: the marginal per-variant cost is indistinguishable from zero. At 5,000 the
#: 4,962,060-variant callset pays 993 JVM launches — about 9.8 h of nothing but
#: reloading the same database, inside a run that projected to ~30 h. At 250,000
#: it pays 20, and SnpEff's share drops from ~16.5 h to ~25-30 min.
#:
#: The bound is memory. 250,000 records is ~918 MB of ``VariantRecord`` (measured
#: at 3,672 bytes each, ``docs/scale-report.md`` §2) plus the adapter results for
#: them, against 24 GB of RAM and a 6 GB JVM heap. It is deliberately not 500,000:
#: the scale sweep puts batch-1M at 7.79 GiB against a 12 GiB watchdog, and there
#: is no throughput left to buy — at 20 launches the fixed cost is already amortised.
#:
#: ``snpeff_local.DEFAULT_BATCH_SIZE`` is kept equal to this. When the SnpEff
#: adapter's own chunk size is the smaller of the two it re-chunks what the service
#: already batched, and the launches come back: 5,000 against 25,000 meant the
#: adapter's chunking never triggered at all, and 250,000 against 25,000 would mean
#: 10 launches per service batch, i.e. 199 of the 993 reloads surviving the fix.
#: Batch size cannot change the answer either way (GP-30,
#: ``test_batch_size_does_not_change_the_answer``); it only decides how many times
#: the database is read off disk.
DEFAULT_ANNOTATION_BATCH_SIZE: int = 250_000


@dataclass(frozen=True, slots=True)
class AnnotatedVariant:
    """One annotated record, with exactly the evidence its annotation licensed.

    The pairing of record and evidence is what lets the caller stream: a consumer
    can write the record to an artifact and push the evidence into a ledger in the
    same step, and then drop both. In the batch API the same values arrive as two
    parallel tuples on :class:`AnnotationResult`, in the same order.
    """

    variant: VariantRecord
    evidence: tuple[EvidenceItem, ...]


def annotate_variants(
    variants: Sequence[VariantRecord],
    *,
    adapters: AdapterSet,
    clock: Clock,
) -> AnnotationResult:
    """Annotate records against the bound adapter set, materialising the result.

    Every input record appears in the output, in input order, annotated or not:
    filtering is a separate concern from annotation (GP-13). Adapters are called
    once with the full deduplicated ID list rather than per variant.

    This stage is the *producer* of ``consequences``, ``population_frequencies``
    and ``clinical_assertions``; it replaces rather than merges those fields, so a
    result always reflects exactly one adapter set. Re-annotating already-annotated
    records is reported in ``warnings``.

    **This is the fixture-scale API.** It holds the input sequence, the annotated
    copy of it and every evidence item simultaneously — three full copies of a
    callset, ~16 GB of ``VariantRecord`` alone at whole-genome scale
    (``docs/scale-report.md`` §2). :func:`iter_annotated` is the same computation
    with the same output, one record at a time; a caller whose record count scales
    with the callset must use it. This function is now a thin drain of that stream
    with ``batch_size=None`` (one batch, one adapter call, exactly as before), so
    the two paths cannot drift.
    """
    stream = iter_annotated(variants, adapters=adapters, clock=clock, batch_size=None)
    annotated: list[VariantRecord] = []
    evidence: list[EvidenceItem] = []
    for annotation in stream:
        annotated.append(annotation.variant)
        evidence.extend(annotation.evidence)
    return AnnotationResult(
        variants=tuple(annotated),
        evidence=tuple(evidence),
        coverage=stream.coverage(),
        warnings=stream.warnings(),
    )


def iter_annotated(
    variants: Iterable[VariantRecord],
    *,
    adapters: AdapterSet,
    clock: Clock,
    batch_size: int | None = DEFAULT_ANNOTATION_BATCH_SIZE,
    triage: VariantTriage | None = None,
    discards: DiscardSink | None = None,
) -> AnnotationStream:
    """Annotate a stream of records without holding the callset.

    Emits one :class:`AnnotatedVariant` per input record, in input order, with the
    same annotations and the same evidence in the same order as
    :func:`annotate_variants` produces for the same input. What changes is only
    what is alive at once: one batch of records and the adapter results for that
    batch, rather than the whole callset three times over.

    Adapters are still called in bulk — ``batch_size`` IDs at a time, deduplicated
    within the batch exactly as the whole-callset call deduplicates globally. A
    per-variant call would be far worse than the memory it saves: a real adapter is
    a tabix seek, an index probe or a subprocess, and 4.5 M of those is not a
    smaller version of 911 of them.

    ``batch_size=None`` means "one batch": every record is read into memory and one
    adapter call is made, which is precisely the pre-streaming behaviour and is what
    :func:`annotate_variants` uses.

    Coverage denominators are *distinct* variant IDs, as before. Duplicates are
    collapsed within a batch, and a duplicate spanning a batch boundary is collapsed
    too (the last ID of the previous batch is carried). For a coordinate-sorted
    stream — which is what ingestion emits, and the only order in which two records
    can share an ID — that is exactly the global deduplication the batch path does.
    An unsorted stream with far-apart duplicates would count a repeat twice, which
    can only *understate* coverage; it can never invent it.

    **``triage`` partitions the pass (ADR 0028).** With a bar bound, the consequence
    adapter is still asked about every record, and the frequency and clinical
    adapters are asked only about the records that clear it. Records that do not
    are yielded to ``discards`` instead of being yielded to the caller: no evidence,
    no artifact row, no lookup. Without a bar this is exactly the pre-ADR-0028 pass,
    which is what every existing caller gets.
    """
    return AnnotationStream(
        variants,
        adapters=adapters,
        clock=clock,
        batch_size=batch_size,
        triage=triage,
        discards=discards,
    )


class AnnotationStream:
    """Annotated records one at a time; coverage and warnings at the end.

    Construct via :func:`iter_annotated`. Single-pass: re-iterating would call every
    adapter a second time and double every counter, so it raises instead.
    """

    __slots__ = (
        "_adapters",
        "_asked",
        "_batch_size",
        "_carry",
        "_clinical_before_bar",
        "_clinical_desc",
        "_consequence_desc",
        "_descriptors",
        "_discards",
        "_distinct",
        "_exhausted",
        "_frequency_desc",
        "_hits",
        "_pre_annotated",
        "_records",
        "_source",
        "_started",
        "_timestamp",
        "_triage",
        "_triage_dropped",
        "_triage_kept",
    )

    def __init__(
        self,
        variants: Iterable[VariantRecord],
        *,
        adapters: AdapterSet,
        clock: Clock,
        batch_size: int | None,
        triage: VariantTriage | None = None,
        discards: DiscardSink | None = None,
    ) -> None:
        if batch_size is not None and batch_size < 1:
            msg = (
                f"batch_size={batch_size} must be at least 1, or None for a single "
                "whole-input batch. A batch size of zero would call the adapters "
                "forever without ever emitting a record."
            )
            raise AnnotationError(msg)
        self._source = variants
        self._adapters = adapters
        self._batch_size = batch_size
        descriptors = adapters.descriptors()
        self._descriptors = descriptors
        self._consequence_desc = descriptors[0]
        self._frequency_desc = descriptors[1]
        self._clinical_desc = descriptors[2] if len(descriptors) > 2 else None
        # Sampled once, before any record is seen, so every evidence item in a run
        # carries the same timestamp however long the stream takes (GP-30).
        self._timestamp = clock.now()
        self._records = 0
        self._pre_annotated = 0
        self._distinct = 0
        #: Distinct IDs the consequence / frequency / clinical adapter answered for.
        self._hits = [0, 0, 0]
        self._carry: str | None = None
        self._started = False
        self._exhausted = False
        # A bar that is present but disabled is treated as absent, so the two
        # spellings of "off" cannot diverge and the hot path has one branch.
        self._triage = triage if (triage is not None and triage.enabled) else None
        self._clinical_before_bar = triage.clinical_before_bar if triage is not None else True
        self._discards = discards
        self._triage_kept: dict[str, int] = dict.fromkeys(TRIAGE_KEEP_REASONS, 0)
        self._triage_dropped: dict[str, int] = dict.fromkeys(TRIAGE_DROP_REASONS, 0)
        #: Distinct IDs each adapter was ASKED about. Equal to `_distinct` for the
        #: consequence adapter always, and for the other two only when no bar is
        #: bound — which is why coverage divides by this and not by `_distinct`.
        #: Dividing a triaged run's frequency hits by the whole callset would
        #: report a coverage collapse that never happened.
        self._asked = [0, 0, 0]

    # ------------------------------------------------------------------ driving

    def __iter__(self) -> Iterator[AnnotatedVariant]:
        if self._started:
            msg = (
                "This annotation stream has already been iterated. Re-iterating would "
                "call every adapter a second time and double every coverage counter. "
                "Call iter_annotated() again against a re-read source."
            )
            raise AnnotationError(msg)
        self._started = True
        return self._drive()

    def _drive(self) -> Iterator[AnnotatedVariant]:
        adapters = self._adapters
        clinical_adapter = adapters.clinical
        triage = self._triage
        for batch in self._batches():
            variant_ids = _unique_ids(record.variant_id for record in batch)
            # The consequence adapter always sees the whole batch. It is one JVM
            # launch per batch and a streaming pass, it is what makes the bar
            # answerable, and it is not the cost this partition exists to avoid.
            consequence_index = adapters.consequence.annotate(variant_ids)

            # The clinical adapter runs BEFORE the bar by default, over the whole
            # batch, and its answer is an input to the verdict. Without that, a
            # curated pathogenic assertion cannot rescue a variant the consequence
            # rule drops, because the rescue needs the lookup the bar skipped —
            # measured on the scale phantom as 16 of 3,616 selected variants lost.
            # ClinVar is one tabix file of ~100 MB against gnomAD's 22 multi-gigabyte
            # shards, so this is the cheap half of what was deferred and the
            # frequency lookup, which is the expensive half, stays behind the bar.
            clinical_upfront = (
                clinical_adapter is not None and triage is not None and self._clinical_before_bar
            )
            early_clinical: Mapping[str, tuple[ClinicalAssertion, ...]] = (
                clinical_adapter.assertions(variant_ids)
                if clinical_upfront and clinical_adapter is not None
                else {}
            )

            verdicts = self._triage_batch(batch, variant_ids, consequence_index, early_clinical)
            # The ONLY ids that reach the frequency and clinical adapters. This is
            # the line the whole of ADR 0028 is about: a gnomAD lookup is a region
            # query against a multi-gigabyte shard, and asking it about 14 thousand
            # coordinates instead of 4.96 million collapses the merged span count
            # by three orders of magnitude.
            lookup_ids = (
                variant_ids
                if triage is None
                else tuple(vid for vid in variant_ids if verdicts[vid].retained)
            )

            frequency_index = adapters.frequency.frequencies(lookup_ids)
            if clinical_upfront:
                clinical_index: Mapping[str, tuple[ClinicalAssertion, ...]] = early_clinical
            else:
                clinical_index = (
                    clinical_adapter.assertions(lookup_ids) if clinical_adapter is not None else {}
                )
            self._observe_batch(
                variant_ids,
                lookup_ids,
                clinical_ids=variant_ids if clinical_upfront else lookup_ids,
                consequence_index=consequence_index,
                frequency_index=frequency_index,
                clinical_index=clinical_index,
            )
            for record in batch:
                self._records += 1
                if (
                    record.consequences
                    or record.population_frequencies
                    or record.clinical_assertions
                ):
                    self._pre_annotated += 1
                if triage is not None:
                    verdict = verdicts[record.variant_id]
                    if not verdict.retained:
                        # Accounted for, not retained (GP-13 as amended, ADR 0028).
                        # One ~75-byte row instead of ~6 KB of JSON and three
                        # EvidenceItems, and the variant stays answerable by id.
                        if self._discards is not None:
                            self._discards.record(verdict.as_discard(record.variant_id))
                        continue
                yield self._annotate(
                    record,
                    consequence_index=consequence_index,
                    frequency_index=frequency_index,
                    clinical_index=clinical_index,
                )
        self._exhausted = True

    def _triage_batch(
        self,
        batch: Sequence[VariantRecord],
        variant_ids: Sequence[str],
        consequence_index: Mapping[str, tuple[ConsequenceAnnotation, ...]],
        early_clinical: Mapping[str, tuple[ClinicalAssertion, ...]],
    ) -> dict[str, TriageVerdict]:
        """One verdict per DISTINCT id in the batch, and the running counts.

        Keyed by id rather than computed per record so that two records sharing an
        id — the only way a coordinate-sorted stream repeats — get one verdict and
        are counted once, exactly as coverage deduplicates them. Returns an empty
        mapping when no bar is bound, which is the allocation the disabled path
        must not pay.
        """
        triage = self._triage
        if triage is None:
            return {}
        # Assertions the ADAPTER just returned, plus any the input record already
        # carried. Both are curated-assertion evidence and either is enough to keep
        # a variant, which makes the bar's clinical test at least as wide as
        # selection's `retain_pathogenic_clinical` — that one fires on PATHOGENIC
        # assertions only, so keeping on any assertion is a strict superset.
        assertions_by_id: dict[str, tuple[ClinicalAssertion, ...]] = dict(early_clinical)
        for record in batch:
            if record.clinical_assertions:
                assertions_by_id.setdefault(record.variant_id, record.clinical_assertions)
        verdicts: dict[str, TriageVerdict] = {}
        for variant_id in variant_ids:
            verdict = triage.verdict(
                consequence_index.get(variant_id, ()),
                clinical_assertions=assertions_by_id.get(variant_id, ()),
            )
            verdicts[variant_id] = verdict
            if verdict.retained:
                self._triage_kept[verdict.reason] += 1
            else:
                self._triage_dropped[verdict.reason] += 1
        return verdicts

    def _batches(self) -> Iterator[list[VariantRecord]]:
        """Successive record batches. ``batch_size=None`` yields exactly one.

        The single-batch case yields even when the input is empty, because the
        pre-streaming implementation called every adapter with an empty ID tuple and
        an adapter is entitled to notice that it was asked. The chunked case skips
        the empty tail instead, so a drained stream makes no pointless final call.
        """
        if self._batch_size is None:
            yield list(self._source)
            return
        iterator = iter(self._source)
        while batch := list(islice(iterator, self._batch_size)):
            yield batch

    def _observe_batch(
        self,
        variant_ids: Sequence[str],
        lookup_ids: Sequence[str],
        *,
        clinical_ids: Sequence[str],
        consequence_index: Mapping[str, Sequence[object]],
        frequency_index: Mapping[str, Sequence[object]],
        clinical_index: Mapping[str, Sequence[object]],
    ) -> None:
        """Fold this batch into the running coverage counters.

        Two denominators, not one. ``_distinct`` counts every id the batch held and
        is the denominator for the consequence adapter, which saw all of them.
        ``_asked`` counts the ids each adapter was actually given, which for the
        frequency and clinical adapters under a triage bar is the surviving subset.

        Dividing a triaged run's frequency hits by the whole callset would report
        that gnomAD coverage collapsed by 350x on the day the bar was switched on.
        It did not: the same fraction of the variants that were *asked about* came
        back. A coverage figure whose denominator is a population the adapter was
        never shown is not a measurement of the adapter (GP-14).
        """
        carry = self._carry
        for variant_id in variant_ids:
            if variant_id == carry:
                # Counted as the tail of the previous batch; counting it again would
                # inflate the denominator and understate coverage.
                continue
            self._distinct += 1
            self._asked[0] += 1
            if consequence_index.get(variant_id):
                self._hits[0] += 1
        for variant_id in lookup_ids:
            if variant_id == carry:
                continue
            self._asked[1] += 1
            if frequency_index.get(variant_id):
                self._hits[1] += 1
        for variant_id in clinical_ids:
            if variant_id == carry:
                continue
            self._asked[2] += 1
            if clinical_index.get(variant_id):
                self._hits[2] += 1
        if variant_ids:
            # The carry is the last id of the batch, whether or not it was looked
            # up: it is the tail-deduplication handle for the NEXT batch, and a
            # triaged-out id can still be the one that repeats across the boundary.
            self._carry = variant_ids[-1]

    def _annotate(
        self,
        record: VariantRecord,
        *,
        consequence_index: Mapping[str, tuple[ConsequenceAnnotation, ...]],
        frequency_index: Mapping[str, tuple[PopulationFrequency, ...]],
        clinical_index: Mapping[str, tuple[ClinicalAssertion, ...]],
    ) -> AnnotatedVariant:
        variant_id = record.variant_id
        consequences = tuple(consequence_index.get(variant_id, ()))
        frequencies = tuple(frequency_index.get(variant_id, ()))
        assertions = tuple(clinical_index.get(variant_id, ()))
        timestamp = self._timestamp

        annotated = record.with_annotations(
            consequences=consequences,
            # Explicitly empty when unknown. `with_annotations` treats None as
            # "leave alone", and leaving stale frequencies alone would be worse
            # than either alternative.
            population_frequencies=frequencies,
            clinical_assertions=assertions,
        )

        evidence: list[EvidenceItem] = [
            _consequence_evidence(
                variant_id, annotation, adapter=self._consequence_desc, timestamp=timestamp
            )
            for annotation in consequences
        ]
        if frequencies:
            evidence.extend(
                _frequency_evidence(
                    variant_id, frequency, adapter=self._frequency_desc, timestamp=timestamp
                )
                for frequency in frequencies
            )
        else:
            evidence.append(
                _missing_frequency_evidence(
                    variant_id, adapter=self._frequency_desc, timestamp=timestamp
                )
            )
        if self._clinical_desc is not None:
            evidence.extend(
                _clinical_evidence(
                    variant_id, assertion, adapter=self._clinical_desc, timestamp=timestamp
                )
                for assertion in assertions
            )
        return AnnotatedVariant(variant=annotated, evidence=tuple(evidence))

    # ------------------------------------------------------------- reporting

    @property
    def drained(self) -> bool:
        """Whether the stream ran to completion. Counters are partial until it did."""
        return self._exhausted

    def coverage(self) -> dict[str, float]:
        """Fraction of the variants each adapter was ASKED about that it answered for.

        Reported even when it is 0.0, because "we asked and got nothing" is a
        finding. Meaningful only once the stream has been consumed; a caller that
        stops early gets the coverage of the prefix it read, and
        :meth:`warnings` says so in words.

        The denominator is per adapter (see :meth:`_observe_batch`). Without a
        triage bar all three denominators are the distinct input count, which is
        the pre-ADR-0028 definition exactly. With one, the frequency and clinical
        denominators are the surviving subset — because a coverage number computed
        over variants an adapter was never shown says nothing about the adapter.
        :meth:`triage_report` carries how many that was.
        """
        return {
            descriptor.name: (self._hits[index] / self._asked[index]) if self._asked[index] else 0.0
            for index, descriptor in enumerate(self._descriptors)
        }

    def triage_report(self) -> TriageReport:
        """What the relevance bar did, in counts. Counts only — never an id (GP-41).

        Meaningful once the stream is drained, like every other counter here. With
        no bar bound the report says so: input equals retained, and the single
        keep reason is ``kept_triage_disabled``.
        """
        kept = sum(self._triage_kept.values())
        dropped = sum(self._triage_dropped.values())
        if self._triage is None:
            return TriageReport(
                input_count=self._records,
                retained_count=self._records,
                retained_by_reason={KEEP_TRIAGE_DISABLED: self._records},
                dropped_by_reason=dict.fromkeys(TRIAGE_DROP_REASONS, 0),
                thresholds={"enabled": False},
            )
        return TriageReport(
            input_count=kept + dropped,
            retained_count=kept,
            retained_by_reason=dict(self._triage_kept),
            dropped_by_reason=dict(self._triage_dropped),
        )

    def warnings(self) -> tuple[str, ...]:
        """Run-level caveats, in a fixed order. These belong in the report, not a log."""
        total = self._distinct
        warnings: list[str] = []

        if self._started and not self._exhausted:
            warnings.append(
                f"Annotation stream was not drained: {self._records} record(s) were "
                "read before the consumer stopped. Every count below describes that "
                "prefix only, and coverage computed from it is not the run's coverage."
            )

        mocked = [descriptor.label for descriptor in self._descriptors if descriptor.synthetic]
        if mocked:
            warnings.append(
                "GP-20: annotation ran against SYNTHETIC stand-in adapter(s) "
                f"{', '.join(mocked)}. The resulting consequences and frequencies are "
                "fabricated demo data and are NOT biologically valid."
            )

        missing_consequences = total - self._hits[0]
        if missing_consequences:
            warnings.append(
                f"{missing_consequences}/{total} variants received no consequence annotation "
                f"from '{self._consequence_desc.name}'. They are retained unannotated rather "
                "than dropped (GP-13); an unannotated variant is not a benign variant."
            )

        asked_frequency = self._asked[1]
        missing_frequencies = asked_frequency - self._hits[1]
        if missing_frequencies:
            warnings.append(
                f"{missing_frequencies}/{asked_frequency} variants have no population-frequency "
                f"data from '{self._frequency_desc.name}'. Their frequencies are recorded as "
                "empty with a neutral no-data evidence item; they must not be scored as allele "
                "frequency 0 (GP-14)."
            )

        # Emitted only when the bar actually removed something. A bar that ran and
        # dropped nothing is a non-event, and `annotation/triage_report.json` says
        # so unconditionally for anyone who wants to check. The clinical-rescue
        # caveat is NOT repeated here: it depends on `selection.retain_pathogenic_
        # clinical`, which this stage cannot see, so `triage.triage_warnings` owns
        # it at the composition root and states it once.
        dropped = sum(self._triage_dropped.values())
        if self._triage is not None and dropped:
            kept = sum(self._triage_kept.values())
            breakdown = ", ".join(
                f"{reason}={self._triage_dropped[reason]}"
                for reason in TRIAGE_DROP_REASONS
                if self._triage_dropped[reason]
            )
            warnings.append(
                f"Annotation triage sent {kept}/{kept + dropped} variants to the frequency and "
                f"clinical adapters and set {dropped} aside on consequence alone"
                + (f" ({breakdown})" if breakdown else "")
                + ". Every one of those is in the discard ledger by id with the rule and the "
                "value that dropped it — accounted for, not retained (GP-13 as amended, "
                "ADR 0028). A variant this bar drops has no frequency lookup, no clinical "
                "lookup and no row in the variants artifact, so absence of any of those for it "
                "is a statement about this pipeline and not about the variant (GP-14)."
            )

        clinical = self._adapters.clinical
        if clinical is None:
            warnings.append(
                "No clinical-assertion adapter was configured, so no curated significance was "
                "sought. Absence of a clinical assertion is not evidence of benignity."
            )
        elif self._hits[2] == 0:
            warnings.append(
                f"Clinical adapter '{clinical.name}' returned no assertions for any of "
                f"the {self._asked[2]} variants it was asked about. Absence of a curated "
                "assertion is not evidence of benignity (GP-14)."
            )

        if self._pre_annotated:
            warnings.append(
                f"{self._pre_annotated}/{self._records} input records already carried "
                "annotations. This stage replaces rather than merges them, so the result "
                "reflects exactly one adapter set."
            )

        return tuple(warnings)


# --------------------------------------------------------------------------- evidence


def _limitations(specific: str, *, adapter: AdapterDescriptor) -> str:
    """Per-claim limitations, plus the mock disclosure when the adapter is a mock."""
    if adapter.synthetic:
        return f"{specific} {SYNTHETIC_STANDIN_LIMITATION}"
    return specific


def _consequence_evidence(
    variant_id: str,
    annotation: ConsequenceAnnotation,
    *,
    adapter: AdapterDescriptor,
    timestamp: datetime,
) -> EvidenceItem:
    """One item per transcript. Two transcripts means two claims, never one merged."""
    claim = (
        f"{variant_id} is predicted to cause {annotation.most_severe_term} in "
        f"{annotation.gene_symbol} on transcript {annotation.transcript_id} "
        f"(predicted impact: {_impact_text(annotation.impact)})."
    )
    payload: dict[str, str | int | float | bool | None] = {
        "gene_symbol": annotation.gene_symbol,
        "transcript_id": annotation.transcript_id,
        "transcript_biotype": annotation.transcript_biotype,
        "is_canonical": annotation.is_canonical,
        "is_mane_select": annotation.is_mane_select,
        "consequence_terms": ",".join(annotation.consequence_terms),
        "impact": _impact_text(annotation.impact),
        "hgvs_c": annotation.hgvs_c,
        "hgvs_p": annotation.hgvs_p,
        "splice_ai_delta_max": annotation.splice_ai_delta_max,
    }
    for score_name in sorted(annotation.pathogenicity_scores):
        payload[f"score_{score_name}"] = annotation.pathogenicity_scores[score_name]

    return _evidence_item(
        subject_id=variant_id,
        category=EvidenceCategory.CONSEQUENCE,
        claim=claim,
        # GP-12: a tool's opinion about a transcript is a prediction. Pairing this
        # with OBSERVED_DATA is rejected by EvidenceItem itself.
        direction=(
            EvidenceDirection.SUPPORTS
            if annotation.impact is not None and annotation.impact in _SUPPORTING_IMPACTS
            else EvidenceDirection.NEUTRAL
        ),
        strength=(
            EvidenceStrength.INSUFFICIENT
            if annotation.impact is None
            else _IMPACT_STRENGTH[annotation.impact]
        ),
        evidence_type=EvidenceType.IN_SILICO_PREDICTION,
        tier=AssertionTier.COMPUTATIONAL_PREDICTION,
        citation=Citation(
            source=annotation.source_tool,
            identifier=f"{variant_id}/{annotation.transcript_id}",
            version=annotation.source_tool_version,
        ),
        method=(
            "Transcript-level consequence lookup via the bound consequence adapter; all "
            "returned transcripts retained, none collapsed to canonical."
        ),
        adapter=adapter,
        specific_limitation=_CONSEQUENCE_LIMITATION,
        timestamp=timestamp,
        # No single scalar magnitude: CADD/REVEL/SpliceAI are in the payload under
        # their own names rather than flattened into one ambiguous number.
        numeric_value=None,
        payload=payload,
    )


def _frequency_evidence(
    variant_id: str,
    frequency: PopulationFrequency,
    *,
    adapter: AdapterDescriptor,
    timestamp: datetime,
) -> EvidenceItem:
    """One item per (source, version, population). Never merged into a single AF."""
    claim = (
        f"{variant_id} has allele frequency {frequency.allele_frequency:.7f} in the "
        f"{frequency.population} population of {frequency.source} {frequency.version}."
    )
    return _evidence_item(
        subject_id=variant_id,
        category=EvidenceCategory.POPULATION,
        claim=claim,
        # Neutral by construction: this stage reports what the cohort holds. Turning a
        # frequency into "too common to be causal" is a thresholded judgement and
        # belongs to prioritisation, where the thresholds are configuration (GP-32).
        direction=EvidenceDirection.NEUTRAL,
        strength=EvidenceStrength.MODERATE,
        evidence_type=EvidenceType.CURATED_DATABASE,
        tier=AssertionTier.DATABASE_ASSERTION,
        # GP-18/GP-12: source and version come from the frequency record itself, and
        # EvidenceItem refuses a DATABASE_ASSERTION whose citation has no version.
        citation=Citation(
            source=frequency.source,
            identifier=f"{variant_id}/{frequency.population}",
            version=frequency.version,
        ),
        method="Keyed lookup in the bound population-frequency adapter.",
        adapter=adapter,
        specific_limitation=_FREQUENCY_LIMITATION,
        timestamp=timestamp,
        numeric_value=frequency.allele_frequency,
        payload={
            "source": frequency.source,
            "version": frequency.version,
            "population": frequency.population,
            "allele_frequency": frequency.allele_frequency,
            "allele_count": frequency.allele_count,
            "allele_number": frequency.allele_number,
            "homozygote_count": frequency.homozygote_count,
            "filter_status": frequency.filter_status,
        },
    )


def _missing_frequency_evidence(
    variant_id: str, *, adapter: AdapterDescriptor, timestamp: datetime
) -> EvidenceItem:
    """The GP-14 item: an explicit, citable record that there is no data.

    Emitting nothing would be indistinguishable from "we never looked", and would
    leave downstream scoring free to read the empty tuple as AF = 0. This item makes
    the gap addressable: it has an ID, it is NEUTRAL, and its claim says in plain
    words that absence is not rarity.
    """
    claim = (
        f"No population-frequency record exists for {variant_id} in {adapter.name} "
        f"{adapter.version}: frequency data is UNAVAILABLE for this variant. Absence of "
        "frequency data is not evidence of rarity and must not be scored as allele "
        "frequency 0 (GP-14)."
    )
    return _evidence_item(
        subject_id=variant_id,
        category=EvidenceCategory.POPULATION,
        claim=claim,
        direction=EvidenceDirection.NEUTRAL,
        strength=EvidenceStrength.INSUFFICIENT,
        # A statement about this pipeline's own lookup, not a database assertion: the
        # cohort asserted nothing, so no assertion may be attributed to it.
        evidence_type=EvidenceType.PIPELINE_INFERENCE,
        tier=AssertionTier.INFERENCE,
        citation=Citation(
            source=adapter.name,
            identifier=f"{variant_id}/no-record",
            version=adapter.version,
        ),
        method=(
            "Keyed lookup in the bound population-frequency adapter returned no rows for "
            "this variant ID."
        ),
        adapter=adapter,
        specific_limitation=_NO_FREQUENCY_LIMITATION,
        timestamp=timestamp,
        # Deliberately None, not 0.0. The whole point of this item is that there is no
        # number to report.
        numeric_value=None,
        payload={
            "frequency_data_available": False,
            "queried_source": adapter.name,
            "queried_version": adapter.version,
        },
    )


def _clinical_evidence(
    variant_id: str,
    assertion: ClinicalAssertion,
    *,
    adapter: AdapterDescriptor,
    timestamp: datetime,
) -> EvidenceItem:
    """A curated significance assertion, filed with its review depth and release."""
    claim = (
        f"{variant_id} is asserted {assertion.significance} by {assertion.source} "
        f"{assertion.version}"
        + (f" (review status: {assertion.review_status})." if assertion.review_status else ".")
    )
    return _evidence_item(
        subject_id=variant_id,
        category=EvidenceCategory.CONSEQUENCE,
        claim=claim,
        direction=_significance_direction(assertion.significance),
        strength=_significance_strength(assertion.star_rating),
        evidence_type=EvidenceType.CURATED_DATABASE,
        tier=AssertionTier.DATABASE_ASSERTION,
        citation=Citation(
            source=assertion.source,
            identifier=assertion.accession or variant_id,
            version=assertion.version,
        ),
        method="Keyed lookup in the bound clinical-assertion adapter.",
        adapter=adapter,
        specific_limitation=_CLINICAL_LIMITATION,
        timestamp=timestamp,
        numeric_value=float(assertion.star_rating) if assertion.star_rating is not None else None,
        payload={
            "source": assertion.source,
            "version": assertion.version,
            "accession": assertion.accession,
            "significance": assertion.significance,
            "review_status": assertion.review_status,
            "star_rating": assertion.star_rating,
            "conditions": ",".join(assertion.conditions),
        },
    )


def _significance_direction(significance: str) -> EvidenceDirection:
    """Map a significance string to a direction, conservatively.

    Anything conflicting, uncertain or unrecognised is NEUTRAL. A VUS is not weak
    support; it is an explicit statement that the curator could not decide.
    """
    token = significance.strip().lower()
    if "conflict" in token:
        return EvidenceDirection.NEUTRAL
    if "pathogenic" in token:
        return EvidenceDirection.SUPPORTS
    if "benign" in token:
        return EvidenceDirection.CONTRADICTS
    return EvidenceDirection.NEUTRAL


def _significance_strength(star_rating: int | None) -> EvidenceStrength:
    """Review depth, not submission count, drives strength."""
    if star_rating is None:
        return EvidenceStrength.WEAK
    if star_rating >= 3:
        return EvidenceStrength.STRONG
    if star_rating == 2:
        return EvidenceStrength.MODERATE
    if star_rating == 1:
        return EvidenceStrength.SUPPORTING
    return EvidenceStrength.WEAK


def _evidence_item(
    *,
    subject_id: str,
    category: EvidenceCategory,
    claim: str,
    direction: EvidenceDirection,
    strength: EvidenceStrength,
    evidence_type: EvidenceType,
    tier: AssertionTier,
    citation: Citation,
    method: str,
    adapter: AdapterDescriptor,
    specific_limitation: str,
    timestamp: datetime,
    numeric_value: float | None,
    payload: dict[str, str | int | float | bool | None],
) -> EvidenceItem:
    """Assemble an EvidenceItem with a content-derived, reproducible ID."""
    return EvidenceItem(
        evidence_id=make_evidence_id(
            subject_id=subject_id, category=category, claim=claim, tool=adapter.name
        ),
        subject_id=subject_id,
        subject_kind="variant",
        claim=claim,
        category=category,
        direction=direction,
        strength=strength,
        evidence_type=evidence_type,
        tier=tier,
        citation=citation,
        method=method,
        tool=adapter.name,
        tool_version=adapter.version,
        limitations=_limitations(specific_limitation, adapter=adapter),
        timestamp=timestamp,
        numeric_value=numeric_value,
        payload=payload,
    )


# --------------------------------------------------------------------------- reporting


def _unique_ids(variant_ids: Iterable[str]) -> tuple[str, ...]:
    """Deduplicate variant IDs while preserving first-seen order.

    Order-preserving rather than sorted: adapters receive IDs in the order the
    caller supplied them, which keeps a debug trace aligned with the input file.
    """
    seen: dict[str, None] = {}
    for variant_id in variant_ids:
        seen.setdefault(variant_id, None)
    return tuple(seen)


def _impact_text(impact: ImpactSeverity | None) -> str:
    """Render an impact for human-facing evidence text.

    ``None`` means the source located the gene but predicted no consequence. It
    renders as "not assessed" rather than as a severity, so a reader can never
    mistake a gene-assignment-only annotation for a prediction of harmlessness
    (GP-14).
    """
    return "not assessed" if impact is None else impact.value
