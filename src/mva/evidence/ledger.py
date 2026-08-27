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

Neither class touches the database. The ledger is the write buffer, the store is
the durable copy, and keeping them apart means a stage can be unit-tested without
a filesystem.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence

from mva.errors import EvidenceError, UnsourcedAssertionError
from mva.models import EvidenceDirection, EvidenceItem

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

    def __init__(self, *, run_id: str) -> None:
        self._run_id = run_id
        self._items: dict[str, EvidenceItem] = {}

    @property
    def run_id(self) -> str:
        return self._run_id

    def add(self, item: EvidenceItem) -> EvidenceItem:
        """Record one item and return the ledger's copy of it.

        An item with no ``run_id`` is stamped with the ledger's. An item carrying a
        *different* ``run_id`` is rejected: mixing runs in one ledger would make
        the provenance of every downstream score ambiguous.

        Re-adding an identical item is a no-op and returns the stored copy.
        Re-using an ID for different content is an integrity violation and raises,
        because content-derived IDs mean a collision is a bug, not a coincidence.
        """
        if item.run_id is None:
            item = item.model_copy(update={"run_id": self._run_id})
        elif item.run_id != self._run_id:
            msg = (
                f"Evidence {item.evidence_id!r} belongs to run {item.run_id!r} but was "
                f"added to the ledger for run {self._run_id!r}."
            )
            raise EvidenceError(msg)

        existing = self._items.get(item.evidence_id)
        if existing is None:
            self._items[item.evidence_id] = item
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

    def items(self) -> tuple[EvidenceItem, ...]:
        """Every item, in a deterministic content-derived order."""
        return tuple(sorted(self._items.values(), key=_sort_key))

    def for_subject(self, subject_id: str) -> tuple[EvidenceItem, ...]:
        """Everything known about one subject, supporting and contradicting alike."""
        return tuple(item for item in self.items() if item.subject_id == subject_id)

    def contradictions(self) -> tuple[EvidenceItem, ...]:
        """Every item whose direction is ``contradicts`` (GP-19)."""
        return tuple(
            item for item in self.items() if item.direction is EvidenceDirection.CONTRADICTS
        )

    def get(self, evidence_id: str) -> EvidenceItem | None:
        """Look up one item, or ``None``. Never raises; see AssertionResolver."""
        return self._items.get(evidence_id)

    def __contains__(self, evidence_id: object) -> bool:
        return isinstance(evidence_id, str) and evidence_id in self._items

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self) -> Iterator[EvidenceItem]:
        return iter(self.items())


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


__all__ = ["AssertionResolver", "EvidenceLedger"]
