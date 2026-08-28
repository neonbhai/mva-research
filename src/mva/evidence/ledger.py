"""The in-memory evidence ledger and the GP-10 assertion gate.

Two objects live here, and the split is deliberate.

:class:`EvidenceLedger` is what pipeline stages write into: an append-only,
deduplicating accumulator that never loses an item and never reorders one run
relative to the next. It exists so that a stage can record *why* it did something
at the moment it does it, rather than reconstructing a justification afterwards.

:class:`AssertionResolver` is what the reporting layer reads through. It is the
mechanical form of GP-10 — "no claim without evidence". A renderer calls
:meth:`AssertionResolver.require` before emitting a sentence; if the claim cites
nothing, or cites an ID the ledger has never seen, the sentence does not get
written. The gate is a hard failure rather than a warning because a warning in a
clinical-adjacent report is a warning nobody reads.

Neither class touches the evidence *store*. The ledger is the write buffer, the
store is the durable copy, and keeping them apart means a stage can be
unit-tested without a filesystem.

**Scale.** The buffer is a dict, which is right for a case and wrong for a
genome: at 4.5 M records, ingestion and annotation together emit on the order of
10-20 M items, tens of gigabytes of hydrated Pydantic models. A ledger given a
``spill_dir`` therefore migrates to a SQLite file once it passes
:data:`DEFAULT_SPILL_THRESHOLD` items and keeps writing there
(:mod:`mva.evidence.spill`). Without a ``spill_dir`` it stays in memory forever,
which is exactly the pre-existing behaviour and is what every test and every
fixture-scale caller gets by default. Spilling changes nothing observable: the
same items, the same total order, the same bytes on a repeat run (GP-30). What it
does change is that :meth:`EvidenceLedger.items` — which materialises a tuple —
stops being a safe call, so :meth:`EvidenceLedger.iter_items` exists beside it and
is what a whole-genome caller must use.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path
from types import TracebackType
from typing import Final

from mva.errors import EvidenceError, UnsourcedAssertionError
from mva.evidence.spill import DEFAULT_FLUSH_BATCH, SqliteEvidenceSpill
from mva.models import EvidenceDirection, EvidenceItem

#: Items held in memory before a ledger with a ``spill_dir`` moves to disk.
#:
#: Chosen so that every fixture-scale run, the synthetic demo and the whole test
#: suite stay entirely in memory — the synthetic case produces low hundreds of
#: items — while a whole-genome run crosses it within the first percent of
#: ingestion. A hydrated ``EvidenceItem`` measures ~3.5 KB resident, so this
#: buffer is the ledger's memory ceiling before the spill takes over: ~175 MB.
#: (It was 200,000, which measured 1.08 GB peak and is more headroom than any
#: caller needs.)
DEFAULT_SPILL_THRESHOLD: Final[int] = 50_000

#: How much of a claim is echoed in an error message. Claims are scientific
#: sentences, not record content, but they are still truncated: an exception
#: message is a log line, and log lines should be scannable.
_CLAIM_EXCERPT = 120


def _sort_key(item: EvidenceItem) -> tuple[str, str, str, str, str]:
    """A total order over evidence.

    Sorting by content rather than by insertion means two runs that produce the
    same evidence in a different order still serialise identically (GP-30). The
    evidence ID is last and is unique, so the order is total and no tie is ever
    resolved by chance.
    """
    return (
        item.subject_kind,
        item.subject_id,
        item.category.value,
        item.direction.value,
        item.evidence_id,
    )


class EvidenceLedger:
    """Append-only, deduplicating accumulator used by pipeline stages.

    Deduplication is by ``evidence_id``, which is content-derived
    (:func:`mva.models.make_evidence_id`): the same conclusion drawn twice is one
    piece of evidence, not two, and must not be able to inflate a score by being
    recorded from two call sites.

    The ledger deliberately has no ``remove`` and no ``clear``. Contradictions and
    failed hypotheses stay (GP-19); a stage that wants to discount evidence
    down-weights it in scoring, where the decision is visible.
    """

    def __init__(
        self,
        *,
        run_id: str,
        spill_dir: Path | None = None,
        spill_threshold: int = DEFAULT_SPILL_THRESHOLD,
        flush_batch: int = DEFAULT_FLUSH_BATCH,
    ) -> None:
        """Open a ledger for one run.

        ``spill_dir`` is the only new decision. ``None`` — the default — keeps
        every item in memory for the life of the ledger, which is what a
        fixture-scale run wants and what every existing caller already gets. A
        path turns on overflow: once ``spill_threshold`` items have accumulated
        the ledger migrates to a SQLite file in that directory and keeps writing
        there, bounding resident memory at roughly the threshold.

        The directory MUST be inside the external workspace. Evidence subjects are
        variant IDs, so the spill file carries proband coordinates and is SENSITIVE
        (GP-40); ``Workspace.tmp_dir`` is the intended home, and
        :meth:`close` unlinks the file.
        """
        if spill_threshold < 1:
            msg = f"spill_threshold={spill_threshold} must be at least 1."
            raise EvidenceError(msg)
        self._run_id = run_id
        self._items: dict[str, EvidenceItem] = {}
        self._spill_dir = spill_dir
        self._spill_threshold = spill_threshold
        self._flush_batch = flush_batch
        self._spill: SqliteEvidenceSpill | None = None
        self._closed_count: int | None = None

    @property
    def run_id(self) -> str:
        return self._run_id

    @property
    def spilled(self) -> bool:
        """Whether this ledger has moved to disk. Reported, never inferred."""
        return self._spill is not None

    def add(self, item: EvidenceItem) -> EvidenceItem:
        """Record one item and return the ledger's copy of it.

        An item with no ``run_id`` is stamped with the ledger's. An item carrying a
        *different* ``run_id`` is rejected: mixing runs in one ledger would make
        the provenance of every downstream score ambiguous.

        Re-adding an identical item is a no-op and returns the stored copy.
        Re-using an ID for different content is an integrity violation and raises,
        because content-derived IDs mean a collision is a bug, not a coincidence.
        Once spilled, that collision is detected at the next flush rather than on
        the call itself — see :meth:`mva.evidence.spill.SqliteEvidenceSpill.add`.
        """
        if item.run_id is None:
            item = item.model_copy(update={"run_id": self._run_id})
        elif item.run_id != self._run_id:
            msg = (
                f"Evidence {item.evidence_id!r} belongs to run {item.run_id!r} but was "
                f"added to the ledger for run {self._run_id!r}."
            )
            raise EvidenceError(msg)

        spill = self._spill
        if spill is not None:
            return spill.add(item)

        existing = self._items.get(item.evidence_id)
        if existing is None:
            self._items[item.evidence_id] = item
            spill_dir = self._spill_dir
            if spill_dir is not None and len(self._items) >= self._spill_threshold:
                self._begin_spill(spill_dir)
            return item
        if existing != item:
            msg = (
                f"Evidence ID {item.evidence_id!r} was reused for a different claim. "
                "Evidence IDs are content-derived; a collision means the ID was "
                "constructed from the wrong inputs."
            )
            raise EvidenceError(msg)
        return existing

    def extend(self, items: Iterable[EvidenceItem]) -> None:
        """Add many items. Fails on the first integrity violation."""
        for item in items:
            self.add(item)

    def _begin_spill(self, spill_dir: Path) -> None:
        """Move the in-memory buffer to disk and keep writing there."""
        spill = SqliteEvidenceSpill(
            spill_dir / f"evidence-ledger-{self._run_id}.sqlite",
            flush_batch=self._flush_batch,
        )
        for item in self._items.values():
            spill.add(item)
        self._items = {}
        self._spill = spill

    def items(self) -> tuple[EvidenceItem, ...]:
        """Every item, in a deterministic content-derived order.

        Materialises. Safe at case scale and nowhere else: a spilled ledger holds
        more items than the process has memory for, which is why it spilled. Use
        :meth:`iter_items` wherever the count scales with the callset.
        """
        return tuple(self.iter_items())

    def iter_items(self) -> Iterator[EvidenceItem]:
        """Every item, in the same order as :meth:`items`, one at a time.

        The spilled path reads through a composite index whose columns are exactly
        the sort key, so this is an ordered scan rather than a sort of everything.
        """
        self._require_open()
        spill = self._spill
        if spill is not None:
            yield from spill.iter_items()
        else:
            yield from sorted(self._items.values(), key=_sort_key)

    def for_subject(self, subject_id: str) -> tuple[EvidenceItem, ...]:
        """Everything known about one subject, supporting and contradicting alike."""
        self._require_open()
        spill = self._spill
        if spill is not None:
            return tuple(spill.iter_for_subject(subject_id))
        # Sort the matching subset, not the whole ledger: the answer is identical
        # because _sort_key is total, and the cost stops scaling with the callset.
        return tuple(
            sorted(
                (item for item in self._items.values() if item.subject_id == subject_id),
                key=_sort_key,
            )
        )

    def contradictions(self) -> tuple[EvidenceItem, ...]:
        """Every item whose direction is ``contradicts`` (GP-19)."""
        self._require_open()
        spill = self._spill
        if spill is not None:
            return tuple(spill.iter_contradictions())
        return tuple(
            sorted(
                (
                    item
                    for item in self._items.values()
                    if item.direction is EvidenceDirection.CONTRADICTS
                ),
                key=_sort_key,
            )
        )

    def get(self, evidence_id: str) -> EvidenceItem | None:
        """Look up one item, or ``None``. Raises only if the ledger was closed."""
        self._require_open()
        spill = self._spill
        if spill is not None:
            return spill.get(evidence_id)
        return self._items.get(evidence_id)

    def close(self) -> None:
        """Release the spill file, if one was opened. Idempotent.

        A ledger that never spilled has nothing to release, so this is a no-op for
        every case-scale caller.

        For one that did, the file is unlinked: it held proband-derived evidence
        and the run is over. ``len()`` keeps working afterwards, from a count taken
        on the way out — the composition root reports ``evidence_count`` in its
        result, and making that depend on whether someone had already closed the
        ledger would be a trap rather than a contract. Reading the *contents* after
        close raises, because returning an empty ledger would look like a run that
        produced no evidence.
        """
        spill = self._spill
        if spill is not None:
            self._closed_count = len(spill)
            self._spill = None
            spill.close()

    def _require_open(self) -> None:
        if self._closed_count is not None:
            msg = (
                f"Evidence ledger for run {self._run_id!r} was closed; its spill file "
                "has been removed. Read the evidence before closing, or persist it to "
                "the evidence store first. (len() still reports the final count.)"
            )
            raise EvidenceError(msg)

    def __enter__(self) -> EvidenceLedger:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self.close()

    def __contains__(self, evidence_id: object) -> bool:
        if not isinstance(evidence_id, str):
            return False
        self._require_open()
        spill = self._spill
        if spill is not None:
            return evidence_id in spill
        return evidence_id in self._items

    def __len__(self) -> int:
        if self._closed_count is not None:
            return self._closed_count
        spill = self._spill
        return len(spill) if spill is not None else len(self._items)

    def __iter__(self) -> Iterator[EvidenceItem]:
        return self.iter_items()


class AssertionResolver:
    """GP-10 enforcement: nothing is stated in a report unless it resolves here.

    The resolver reads a ledger; it never writes to one. That asymmetry is the
    point — a renderer that could add evidence could manufacture its own
    justification.
    """

    def __init__(self, ledger: EvidenceLedger) -> None:
        self._ledger = ledger

    @property
    def ledger(self) -> EvidenceLedger:
        return self._ledger

    def resolve(self, evidence_ids: Sequence[str]) -> tuple[EvidenceItem, ...]:
        """Best-effort lookup: the items that exist, in the order requested.

        Unknown IDs are silently skipped. Use this only where a partial answer is
        meaningful (rendering the citations a claim *does* have). Anything that
        gates output must use :meth:`require`.
        """
        seen: set[str] = set()
        resolved: list[EvidenceItem] = []
        for evidence_id in evidence_ids:
            if evidence_id in seen:
                continue
            seen.add(evidence_id)
            item = self._ledger.get(evidence_id)
            if item is not None:
                resolved.append(item)
        return tuple(resolved)

    def require(self, claim: str, evidence_ids: Sequence[str]) -> tuple[EvidenceItem, ...]:
        """Resolve every ID or refuse the claim.

        Raises :class:`~mva.errors.UnsourcedAssertionError` when the citation list
        is empty or names an ID the ledger does not hold. Both failures are the
        same failure in practice: a sentence about a patient with nothing behind
        it.

        Error messages carry the claim, the missing IDs and counts — never record
        content, because a traceback goes everywhere the process's logs go.
        """
        if not evidence_ids:
            msg = (
                f"Unsourced assertion (GP-10): {_excerpt(claim)!r} cites no evidence. "
                "Every rendered statement must cite at least one EvidenceItem ID; "
                "if the claim is genuinely unsupported, record it as an OpenQuestion "
                "instead of stating it."
            )
            raise UnsourcedAssertionError(msg)

        missing = [eid for eid in evidence_ids if self._ledger.get(eid) is None]
        if missing:
            msg = (
                f"Unsourced assertion (GP-10): {_excerpt(claim)!r} cites "
                f"{len(missing)} of {len(evidence_ids)} evidence ID(s) that do not "
                f"resolve in the ledger for run {self._ledger.run_id!r}: "
                f"{sorted(missing)}."
            )
            raise UnsourcedAssertionError(msg)

        return self.resolve(evidence_ids)

    def contradictions_for(self, evidence_ids: Sequence[str]) -> tuple[EvidenceItem, ...]:
        """The subset of the cited evidence that argues against the claim.

        A renderer calls this so that a hedge can be attached automatically: a
        statement cited by evidence that partly contradicts it must not read as
        settled (GP-19).
        """
        return tuple(
            item
            for item in self.resolve(evidence_ids)
            if item.direction is EvidenceDirection.CONTRADICTS
        )


def _excerpt(claim: str) -> str:
    stripped = claim.strip()
    if len(stripped) <= _CLAIM_EXCERPT:
        return stripped
    return stripped[: _CLAIM_EXCERPT - 1] + "…"


__all__ = ["DEFAULT_SPILL_THRESHOLD", "AssertionResolver", "EvidenceLedger"]
