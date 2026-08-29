"""The discard ledger: every dropped variant, by id, with the rule that dropped it.

ADR 0028. This is the mechanism that discharges the amended GP-13 — the
requirement is that a dropped candidate be *accounted for*, not that it be
retained — and it is the reviewer-facing artifact `docs/repo-followup-backlog.md`
§1.3 asked for.

Why a TSV and not the evidence store
------------------------------------
The evidence store is DuckDB and Parquet behind an ``EvidenceLedger`` that spills
to SQLite past a threshold (ADR 0002, ADR 0004). It is the right home for a claim
that a report will cite. It is the wrong home for five million rows nothing will
cite, and ADR 0019 already declined to put them there: "four million items of the
exact kind the ledger spill exists to survive, describing candidates nobody will
read".

A discard row is 4 fields and about 75 bytes. Five million of them is ~375 MB,
written straight through with no index, no schema and no transaction. The same
information as ``EvidenceItem``s would be ~65 GB — on a run whose whole problem
was that it projected 95 GB of intermediates onto a 60 GB volume.

The format is the query interface. ``grep`` answers "why isn't variant X in the
output?" against a 375 MB file in a second, and :func:`explain_variant` is the
same scan with the reasons spelled out. Neither loads the file.

What this is not
----------------
It is **not** a ranking of near-misses. ADR 0019 was right that the honest way to
see what a stage cost you is to re-run with the stage off; this file tells you
*which* variants and *on what value*, which is enough to decide whether re-running
is worth it and enough to check a specific variant by name. It carries no score,
no evidence and no claim.

Privacy
-------
Every row holds a build-qualified proband coordinate, so the file is SENSITIVE and
lives in the workspace beside the evidence spill (GP-40, ADR 0006). It is
registered as ``ArtifactKind.DISCARD_LEDGER`` and is not on the export allowlist.
:meth:`DiscardLedger.counts` returns counts only and is what a report may state
(GP-41).
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from pathlib import Path
from types import TracebackType
from typing import IO, Final

from mva.models.discard import (
    DISCARD_LEDGER_COLUMNS,
    FIELD_SEPARATOR,
    DiscardRecord,
    DiscardStage,
)

#: Written above the header. Says what the file is, that it is not evidence, and
#: that reading it means reading proband coordinates — a reader who opens it by
#: accident should learn all three from the first screen.
BANNER: Final[tuple[str, ...]] = (
    "# mva discard ledger — one row per variant this run dropped.",
    "# GP-13 (amended, ADR 0028): a dropped candidate is ACCOUNTED FOR, not retained.",
    "# Columns: variant_id, stage, rule, value. `value` is the observation the rule",
    "# fired on; empty means the rule fired on an ABSENCE and no value exists (GP-14).",
    "# This is NOT evidence and carries no claim. It holds proband coordinates and is",
    "# classified SENSITIVE: it never leaves the workspace (GP-40, GP-41).",
)


class DiscardLedger:
    """An open discard ledger. One row per dropped variant, written straight through.

    Streaming and append-only by construction: nothing is held, nothing is sorted,
    nothing is rewritten. Rows land in the order the pipeline dropped them, which
    for a coordinate-sorted callset is coordinate order within each stage — and
    which is deterministic for a fixed input either way (GP-30).

    Used as a context manager, or opened and closed explicitly by the composition
    root, which must close it before the artifact is hashed.
    """

    __slots__ = ("_closed", "_counts", "_handle", "_path", "_rows")

    def __init__(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        self._path = path
        self._handle: IO[str] = path.open("w", encoding="utf-8")
        self._counts: dict[str, int] = {}
        self._rows = 0
        self._closed = False
        for line in BANNER:
            self._handle.write(line + "\n")
        self._handle.write("#" + FIELD_SEPARATOR.join(DISCARD_LEDGER_COLUMNS) + "\n")

    # ------------------------------------------------------------------ writing

    def record(self, discard: DiscardRecord) -> None:
        """Account for one dropped variant.

        Silently ignored after :meth:`close`. A stage that drops a record while
        the run is unwinding should not raise a second exception on top of the
        first one, and an un-ledgered discard during teardown is a smaller
        problem than a lost traceback.
        """
        if self._closed:
            return
        self._handle.write(discard.as_row())
        self._rows += 1
        key = discard.counter_key
        self._counts[key] = self._counts.get(key, 0) + 1

    def drop(self, variant_id: str, *, stage: DiscardStage, rule: str, value: str = "") -> None:
        """:meth:`record` without constructing the dataclass at the call site."""
        self.record(DiscardRecord(variant_id=variant_id, stage=stage, rule=rule, value=value))

    # ------------------------------------------------------------------ reading

    @property
    def path(self) -> Path:
        return self._path

    @property
    def row_count(self) -> int:
        return self._rows

    def counts(self) -> dict[str, int]:
        """``stage/rule -> count``, in a fixed order. Counts only — no ids (GP-41).

        Sorted by key rather than by insertion, so the mapping this becomes in an
        artifact is byte-stable between runs whatever order the drops arrived in.
        """
        return {key: self._counts[key] for key in sorted(self._counts)}

    def counts_by_stage(self) -> dict[str, int]:
        """Total discards per stage, in :class:`DiscardStage` declaration order."""
        totals: dict[str, int] = {str(stage): 0 for stage in DiscardStage}
        for key, count in self._counts.items():
            stage = key.split("/", 1)[0]
            if stage in totals:
                totals[stage] += count
        return totals

    # ------------------------------------------------------------------ lifecycle

    def close(self) -> None:
        """Flush and close. Idempotent."""
        if self._closed:
            return
        self._closed = True
        self._handle.close()

    @property
    def closed(self) -> bool:
        return self._closed

    def __enter__(self) -> DiscardLedger:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self.close()


# ---------------------------------------------------------------------------
# The read side
# ---------------------------------------------------------------------------


def iter_discards(path: Path) -> Iterator[DiscardRecord]:
    """Every parseable row, in file order. One row resident at a time.

    Unparseable rows are skipped rather than raised on — see
    :meth:`DiscardRecord.from_row`. A ledger truncated by the disk exhaustion the
    pre-flight estimator exists to prevent is still readable up to the truncation,
    and that is the run where someone most needs to read it.
    """
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            record = DiscardRecord.from_row(line)
            if record is not None:
                yield record


def explain_variant(path: Path, variant_id: str) -> tuple[DiscardRecord, ...]:
    """Every discard row for one variant. Empty means it was never dropped.

    "Empty" is genuinely ambiguous and the caller must say so rather than
    resolving it: a variant absent from this file either survived to the output or
    was never in the input at all. The ledger is a record of removals, not a
    census (GP-14). :func:`mva.cli` renders that distinction in words.

    Matching is exact on the canonical id. A partial coordinate finds nothing,
    which is the right answer — a pipeline that guessed which of two variants at
    one position you meant would be worse than one that said no.
    """
    return tuple(record for record in iter_discards(path) if record.variant_id == variant_id)


def count_discards(path: Path) -> Mapping[str, int]:
    """``stage/rule -> count`` read back from a finished ledger.

    The same mapping :meth:`DiscardLedger.counts` returns, recomputed from the
    bytes. Two of them disagreeing means the file was truncated, which is worth
    knowing and is why this exists separately from the in-memory counter.
    """
    counts: dict[str, int] = {}
    for record in iter_discards(path):
        key = record.counter_key
        counts[key] = counts.get(key, 0) + 1
    return {key: counts[key] for key in sorted(counts)}
