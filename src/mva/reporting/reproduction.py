"""Does a produced ranking reproduce the submitted one, without anyone reading either?

`docs/track1-method.md` §4 states the gap plainly: the submitted Track 1 ranking
came out of a scripted run, `iter_selected` still applies ADR 0019's rules rather
than the rule that run applied, and nothing had been executed end to end by the
library against the real callset. Once it has been, one question is left — does
the produced CSV propose the same candidates in the same order? — and this module
answers exactly that.

**Both inputs are patient data.** A Track 1 submission is a list of one child's
variant coordinates. The submitted file is not in this repository and never will
be (GP-40), and an agent runs this comparison, so its output lands in that agent's
context (PRIV-03, ADR 0008). The module therefore holds one invariant, which is
stronger than an instruction to be careful:

    **No byte read from either file reaches a return value, a log record or an
    exception message.**

Everything emitted is a count, a boolean, a correlation coefficient, or a verdict
from a fixed vocabulary. Files are named by ROLE (``submitted`` / ``produced``)
and by basename; the directory is withheld because it is the workspace location,
which :mod:`mva.cli` already refuses to print. Malformed input is reported as a
line NUMBER, never a line (GP-41, PRIV-09). Every exception raised while parsing
is re-raised ``from None``, because a chained :class:`UnicodeDecodeError` carries
the offending bytes in its own message and a chained ``csv.Error`` can carry a
field — the two ways this invariant would otherwise be lost to a traceback.

The comparison itself is the scorer's. Rows are keyed with
:func:`~mva.reporting.track1.submission_row_key` and ranked the way
``evaluation.py`` ranks them, ``sorted(enumerate(rows), key=(-epcr, file_index))``
— so "same order" means the order the challenge would score, not the order the
bytes happen to sit in.
"""

from __future__ import annotations

import csv
import io
import math
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

from mva.errors import SubmissionFormatError
from mva.reporting.track1 import MAX_SUBMISSION_ROWS, TRACK1_COLUMNS, submission_row_key

#: One variant as the scorer keys it: ``(chrom, pos, ref, alt)``.
VariantKey = tuple[str, int, str, str]

#: One submission row. A frozenset because the scorer matches rows by frozenset
#: equality, so a compound-het pair is order-free within its row.
RowKey = frozenset[VariantKey]

#: The mandatory first variant. Read from this constant and never from the parsed
#: record's own keys: iterating the file's header would put file bytes into a
#: message, which is the one thing this module may not do.
_FIRST_VARIANT_COLUMNS = ("chrom_1", "pos_1", "ref_1", "alt_1")

#: The k values top-k agreement is reported at.
#:
#: Not a round-numbers choice. They are the challenge scorer's own rank-point tier
#: boundaries, ``[(1,100),(3,50),(5,25),(10,10)]``: rank 1 scores 100, ranks 2-3
#: score 50, 4-5 score 25, 6-10 score 10, and everything below scores nothing. So
#: agreement at each k is agreement about a scoring tier, and a disagreement that
#: shows up first at k=5 is worth strictly less than one that shows up at k=1.
TOP_K_TIERS: tuple[int, ...] = (1, 3, 5, 10)


class ReproductionVerdict(StrEnum):
    """The five ways two rankings can relate. Closed vocabulary, by design."""

    IDENTICAL = "identical"
    SAME_SET_DIFFERENT_ORDER = "same-set-different-order"
    SUPERSET = "superset"
    SUBSET = "subset"
    DIVERGENT = "divergent"


#: What each verdict means for the question "did the library reproduce the run?".
#:
#: Held here rather than in the CLI because the verdict word alone is thin, and a
#: reader who gets `superset` without being told which direction it points has
#: been given a token rather than an answer. A test asserts the map is total.
VERDICT_MEANINGS: dict[ReproductionVerdict, str] = {
    ReproductionVerdict.IDENTICAL: (
        "same rows, same scored order — the produced ranking reproduces the submission"
    ),
    ReproductionVerdict.SAME_SET_DIFFERENT_ORDER: (
        "same rows, different scored order — the candidates agree, the ranking does not"
    ),
    ReproductionVerdict.SUPERSET: ("produced proposes every submitted row and at least one more"),
    ReproductionVerdict.SUBSET: ("produced proposes only a strict subset of the submitted rows"),
    ReproductionVerdict.DIVERGENT: ("each file proposes at least one row the other does not"),
}


@dataclass(frozen=True, slots=True)
class TopKAgreement:
    """How many of the produced top-k rows are in the submitted top-k.

    ``produced_window`` and ``submitted_window`` are ``min(k, rows)`` for each
    file, and are carried rather than recomputed because they are the difference
    between "0 of 3 agree" and "0 of 0 rows exist to agree" (GP-14). A ratio alone
    cannot tell those apart, and the second is not a disagreement.
    """

    k: int
    agreeing: int
    produced_window: int
    submitted_window: int


@dataclass(frozen=True, slots=True)
class ReproductionReport:
    """Agreement statistics for two submissions. Counts, flags and a verdict only.

    Nothing on this dataclass can hold a coordinate: every field is an ``int``, a
    ``bool``, a bounded ``float`` or a member of a closed enum. That is deliberate
    — a report object that *could* carry a row is one refactor away from printing
    one, and the caller here is an agent (ADR 0008).
    """

    submitted_rows: int
    produced_rows: int
    intersection: int
    submitted_only: int
    produced_only: int
    top1_identical: bool
    top_k: tuple[TopKAgreement, ...]
    spearman: float | None
    verdict: ReproductionVerdict

    @property
    def reproduced(self) -> bool:
        """True only for :attr:`ReproductionVerdict.IDENTICAL`.

        Deliberately strict, because the question this module was written to
        answer is "does it reproduce the submitted ranking?" and every other
        verdict is a no with a reason attached. The CLI exits non-zero on this.
        """
        return self.verdict is ReproductionVerdict.IDENTICAL


@dataclass(frozen=True, slots=True)
class _RankedRow:
    """One parsed row, with everything needed to rank it and nothing else."""

    line: int
    order: int
    epcr: float
    key: RowKey


def _where(role: str, path: Path) -> str:
    """How a file is named in a message: its role and its basename.

    Not its directory. The submitted CSV lives on the case volume, and
    :mod:`mva.cli` states as a load-bearing property that it prints no workspace
    absolute path (GP-44). The role is what a reader actually needs anyway — with
    two files both called ``track1_submission.csv``, the parent directory
    disambiguates nothing that ``submitted`` / ``produced`` has not already said.
    """
    return f"{role} file {path.name!r}"


def _malformed(role: str, path: Path, detail: str) -> SubmissionFormatError:
    """Build the refusal. **Raise it outside the `except` block that produced it.**

    ``raise ... from None`` is not enough here, and the difference is a disclosure.
    It sets ``__suppress_context__``, which stops the standard traceback formatter
    printing the original — but ``__context__`` still *holds* it, and the original
    is ``ValueError: could not convert string to float: '<the row's epcr>'``.
    Anything that reads the attribute rather than asking the formatter — pytest's
    failure repr does, and so does any log handler that walks the chain — prints
    the field verbatim (PRIV-09). This was not hypothetical: the leak test in
    ``tests/unit/test_track1_reproduction.py`` found it, with the fixture token in
    the pytest output.

    The interpreter sets ``__context__`` at raise time, so no factory can clear it.
    The only reliable fix is to leave the handler first: an exception raised after
    the ``try`` statement finishes has no exception being handled to chain to. Each
    caller below therefore records what went wrong and raises on the way out.
    """
    return SubmissionFormatError(f"{_where(role, path)}: {detail}")


def _read_text(path: Path, role: str) -> str:
    """Read the file as UTF-8, converting a decode failure into a clean error.

    ``UnicodeDecodeError`` puts the offending bytes in ``str(exc)``, so it may
    neither propagate nor be chained. The byte OFFSET is reportable and the byte
    is not, which is what the message says.
    """
    detail = ""
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        detail = (
            f"is not valid UTF-8 (first bad byte at offset {exc.start}). A Track 1 "
            "submission is a UTF-8 CSV; the byte itself is withheld under GP-41."
        )
    except OSError as exc:
        detail = f"could not be read (errno {exc.errno})."
    raise _malformed(role, path, detail)


def _header_of(reader: csv.DictReader[str], role: str, path: Path) -> tuple[str, ...]:
    """The header row, validated against the contract without echoing it."""
    unparseable = False
    try:
        fieldnames = tuple(reader.fieldnames or ())
    except csv.Error:
        fieldnames = ()
        unparseable = True
    if unparseable:
        raise _malformed(role, path, "line 1: header row is not parseable as CSV.")

    if not fieldnames:
        raise _malformed(
            role,
            path,
            "is empty: no header row. The Track 1 contract requires the 12-column "
            "header, so a file without one is not a submission.",
        )
    if fieldnames != TRACK1_COLUMNS:
        raise _malformed(
            role,
            path,
            f"line 1: header has {len(fieldnames)} column(s) and does not match the "
            f"Track 1 contract's {len(TRACK1_COLUMNS)} columns in order (see "
            "mva.reporting.track1.TRACK1_COLUMNS). The columns found are withheld: "
            "a header can be edited to hold data, so it is treated as content.",
        )
    return fieldnames


def _records(
    reader: csv.DictReader[str], role: str, path: Path
) -> list[tuple[int, dict[str, str | None]]]:
    """Every data row paired with its 1-based line number in the file.

    ``reader.line_num`` rather than an enumerate counter: ``notes`` is a quoted
    field and may legally hold a newline, in which case the two disagree — and a
    diagnostic that names the wrong line costs a reader more than it saves.
    """
    records: list[tuple[int, dict[str, str | None]]] = []
    unparseable = False
    try:
        for record in reader:
            records.append((reader.line_num, record))
    except csv.Error:
        unparseable = True
    if unparseable:
        raise _malformed(role, path, f"line {reader.line_num}: row is not parseable as CSV.")

    if len(records) > MAX_SUBMISSION_ROWS:
        raise _malformed(
            role,
            path,
            f"has {len(records)} data rows; the scorer raises ValueError above "
            f"{MAX_SUBMISSION_ROWS}. This file would not have been accepted, so "
            "comparing against it would compare against something unsubmittable.",
        )
    return records


def _epcr_of(record: dict[str, str | None], role: str, path: Path, line: int) -> float:
    """The row's EPCR, which is what ranks it. Never reported, only ordered by.

    An unparseable value and a NaN collapse into one message so the refusal can be
    raised outside the handler — see :func:`_malformed`. They also deserve the same
    answer: neither can be ordered, and the scorer rejects both.
    """
    try:
        epcr = float((record.get("epcr") or "").strip())
    except ValueError:
        epcr = math.nan
    if not math.isfinite(epcr):
        raise _malformed(
            role,
            path,
            f"line {line}: epcr is not a finite number, so the row cannot be ranked. "
            "The value is withheld (GP-41).",
        )
    return epcr


def _key_of(record: dict[str, str | None], role: str, path: Path, line: int) -> RowKey:
    """The scorer's row key, or a line-numbered refusal."""
    # Iterate the CONTRACT's column names, never the record's own keys: the
    # record's keys come from the file's header, and naming them would put file
    # bytes into a message.
    blank = [column for column in _FIRST_VARIANT_COLUMNS if not (record.get(column) or "").strip()]
    if blank:
        raise _malformed(
            role,
            path,
            f"line {line}: {len(blank)} of the {len(_FIRST_VARIANT_COLUMNS)} mandatory "
            f"first-variant columns are blank ({sorted(blank)}). A row with no first "
            "variant has no comparable key.",
        )
    key = submission_row_key(record)
    if key is None:
        raise _malformed(
            role,
            path,
            f"line {line}: the variant columns do not parse into a scorer key "
            "(pos must be a plain integer). The values are withheld (GP-41).",
        )
    return key


def read_ranked_keys(path: Path, *, role: str) -> tuple[RowKey, ...]:
    """Row keys for one submission, in the order the challenge scorer would rank them.

    ``role`` is ``"submitted"`` or ``"produced"`` and appears in every message
    raised about this file; it is how a reader tells two files apart without being
    told where either one lives.

    Ordering is the scorer's, ``sorted(enumerate(rows), key=(-epcr, file_index))``,
    which is total and so leaves no set- or dict-order dependence anywhere in this
    module (GP-30). Pre-sorting is not required of a submission, and a file that
    is not pre-sorted still has a well-defined ranking; comparing file order
    directly would report a difference the scorer cannot see.
    """
    reader = csv.DictReader(io.StringIO(_read_text(path, role), newline=""))
    _header_of(reader, role, path)

    rows: list[_RankedRow] = []
    first_seen: dict[RowKey, int] = {}
    for order, (line, record) in enumerate(_records(reader, role, path)):
        key = _key_of(record, role, path, line)
        seen_at = first_seen.setdefault(key, line)
        if seen_at != line:
            raise _malformed(
                role,
                path,
                f"lines {seen_at} and {line} propose the same variant set. A repeated "
                "row has no single rank, so rank agreement against it is undefined.",
            )
        rows.append(
            _RankedRow(line=line, order=order, epcr=_epcr_of(record, role, path, line), key=key)
        )

    rows.sort(key=lambda row: (-row.epcr, row.order))
    return tuple(row.key for row in rows)


def spearman_over_intersection(
    submitted: Sequence[RowKey], produced: Sequence[RowKey]
) -> float | None:
    """Spearman's rho between the two rankings, over the rows both propose.

    Implemented here rather than pulled in: scipy is not a dependency of this
    project and adding one to compute a closed-form statistic over at most ten
    ranks would be a poor trade (see ``pyproject.toml``).

    Rows are **re-ranked within the intersection**. Using each row's rank in its
    own whole file instead would let a row the two files disagree about shift
    every rank beneath it, so the coefficient would measure set disagreement — a
    thing already reported exactly, as three counts — rather than order
    disagreement, which is what nothing else here measures.

    Ranks are distinct integers by construction (the scorer's sort is total), so
    there are no ties and the closed form ``1 - 6·Σd² / (n(n²-1))`` is exact
    rather than an approximation of the Pearson-on-ranks definition.

    Returns ``None`` when fewer than two rows are shared. Rank correlation is
    undefined there — with one shared row there is nothing to be in the wrong
    order — and ``None`` says so. Returning 0.0 would state "no correlation",
    which is a measurement, and absence of a measurement is not one (GP-14).
    """
    shared = frozenset(submitted) & frozenset(produced)
    if len(shared) < 2:
        return None

    # Both comprehensions walk an ORDERED sequence and use the set for membership
    # only, so nothing here depends on set iteration order (GP-30).
    produced_rank = {key: rank for rank, key in enumerate(k for k in produced if k in shared)}
    n = len(shared)
    squared = sum(
        (rank - produced_rank[key]) ** 2
        for rank, key in enumerate(k for k in submitted if k in shared)
    )
    return 1.0 - (6.0 * squared) / (n * (n * n - 1))


def _top_k_agreement(
    submitted: Sequence[RowKey], produced: Sequence[RowKey]
) -> tuple[TopKAgreement, ...]:
    """Overlap between the two top-k windows, at each scoring tier boundary."""
    agreements: list[TopKAgreement] = []
    for k in TOP_K_TIERS:
        window = frozenset(submitted[:k])
        agreements.append(
            TopKAgreement(
                k=k,
                agreeing=sum(1 for key in produced[:k] if key in window),
                produced_window=min(k, len(produced)),
                submitted_window=min(k, len(submitted)),
            )
        )
    return tuple(agreements)


def _verdict(submitted: Sequence[RowKey], produced: Sequence[RowKey]) -> ReproductionVerdict:
    """Classify the relationship. Order is only consulted once the sets agree."""
    submitted_set = frozenset(submitted)
    produced_set = frozenset(produced)
    if submitted_set == produced_set:
        if tuple(submitted) == tuple(produced):
            return ReproductionVerdict.IDENTICAL
        return ReproductionVerdict.SAME_SET_DIFFERENT_ORDER
    if submitted_set < produced_set:
        return ReproductionVerdict.SUPERSET
    if produced_set < submitted_set:
        return ReproductionVerdict.SUBSET
    return ReproductionVerdict.DIVERGENT


def compare_rankings(submitted: Sequence[RowKey], produced: Sequence[RowKey]) -> ReproductionReport:
    """Agreement statistics for two rank-ordered key sequences.

    Separate from :func:`compare_submission_files` so the statistics can be
    exercised without a file on disk — which is what lets the test suite cover
    this at all, given that no real submission may enter the repository (GP-40).
    """
    submitted_set = frozenset(submitted)
    produced_set = frozenset(produced)
    return ReproductionReport(
        submitted_rows=len(submitted),
        produced_rows=len(produced),
        intersection=len(submitted_set & produced_set),
        submitted_only=len(submitted_set - produced_set),
        produced_only=len(produced_set - submitted_set),
        # False when either file has no rank-1 row at all. The row counts, reported
        # alongside, are what distinguish "the top rows differ" from "there is no
        # top row to compare" (GP-14).
        top1_identical=bool(submitted) and bool(produced) and submitted[0] == produced[0],
        top_k=_top_k_agreement(submitted, produced),
        spearman=spearman_over_intersection(submitted, produced),
        verdict=_verdict(submitted, produced),
    )


def compare_submission_files(submitted: Path, produced: Path) -> ReproductionReport:
    """Compare two Track 1 submission CSVs and return counts only.

    ``submitted`` is the file that was uploaded; ``produced`` is what a run
    generated. The order matters to the verdict and to nothing else:
    :attr:`ReproductionVerdict.SUPERSET` means the *produced* file proposes
    everything the submitted one did and more.
    """
    return compare_rankings(
        read_ranked_keys(submitted, role="submitted"),
        read_ranked_keys(produced, role="produced"),
    )


__all__ = [
    "TOP_K_TIERS",
    "VERDICT_MEANINGS",
    "ReproductionReport",
    "ReproductionVerdict",
    "RowKey",
    "TopKAgreement",
    "VariantKey",
    "compare_rankings",
    "compare_submission_files",
    "read_ranked_keys",
    "spearman_over_intersection",
]
