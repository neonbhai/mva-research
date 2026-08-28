"""The committed record of what PubMed returned, and the comparison against it.

Offline by construction: this module opens files, not sockets.

The ledger is what makes the citation check runnable in CI on the sensitive path.
Without it the only way to know a PMID is real is to ask NCBI, and PRIV-05 forbids
asking from anywhere the pipeline runs. With it, the always-on test compares our
prose against a checked-in snapshot of the authoritative record, and the network
appears exactly once, in an opt-in refresh.

The trade is stated rather than hidden: a ledger is a *cache*, so an always-on pass
means "our document agrees with PubMed as of `resolved_utc`", not "as of now". The
online test closes that gap when it is run. A record that is never refreshed is a
record that ages, and `resolved_utc` is in the file so a reviewer can see how much.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any, Final

from tools.citations.extract import (
    CitationForm,
    WrittenCitation,
    normalise_journal,
    normalise_surname,
    normalise_title,
)

__all__ = [
    "LEDGER_PATH",
    "LEDGER_SCHEMA",
    "Ledger",
    "LedgerRecord",
    "Mismatch",
    "MismatchKind",
    "compare",
    "format_mismatches",
    "load_ledger",
    "write_ledger",
]

LEDGER_SCHEMA: Final[str] = "mva.citations.ledger/1"

#: Relative to the repository root.
LEDGER_PATH: Final[str] = "knowledge/literature/citations.verified.json"


@dataclass(frozen=True, slots=True)
class LedgerRecord:
    """One PubMed esummary record, reduced to the fields a citation asserts."""

    pmid: str
    first_author: str
    first_author_surname: str
    journal: str
    pubdate: str
    epubdate: str
    year: int | None
    epub_year: int | None
    title: str


@dataclass(frozen=True, slots=True)
class Ledger:
    """A dated snapshot of PubMed for a set of PMIDs."""

    schema: str
    resolved_utc: str
    source: str
    records: dict[str, LedgerRecord]


class MismatchKind(StrEnum):
    """What disagreed. Each is a different bug with a different fix."""

    #: A PMID we cite is not in the ledger at all. Refresh the ledger.
    NOT_IN_LEDGER = "not_in_ledger"
    #: The ledger carries a PMID nothing cites any more. Refresh the ledger.
    ORPHANED_IN_LEDGER = "orphaned_in_ledger"
    #: THE ONE THIS EXISTS FOR: the PMID resolves to a paper by someone else.
    FIRST_AUTHOR = "first_author"
    #: The year we wrote is neither the print year nor the electronic year.
    YEAR = "year"
    #: The journal we wrote is not the journal of record.
    JOURNAL = "journal"
    #: The title we wrote is not the title of record. The strongest signal here:
    #: two papers can share an author, a year and a journal; none share a title.
    TITLE = "title"


@dataclass(frozen=True, slots=True)
class Mismatch:
    """One disagreement between our prose and the record, in both directions."""

    kind: MismatchKind
    pmid: str
    source: str
    line: int
    wrote: str
    pubmed: str
    context: str


def load_ledger(path: Path) -> Ledger:
    """Read the ledger, refusing an unknown schema rather than guessing at it."""
    payload: Any = json.loads(path.read_text(encoding="utf-8"))
    schema = payload.get("schema")
    if schema != LEDGER_SCHEMA:
        msg = (
            f"{path} declares schema {schema!r}, expected {LEDGER_SCHEMA!r}. "
            f"Refusing to read it: a ledger read under the wrong schema would "
            f"silently verify the wrong fields."
        )
        raise ValueError(msg)
    return Ledger(
        schema=schema,
        resolved_utc=str(payload["resolved_utc"]),
        source=str(payload["source"]),
        records={
            pmid: LedgerRecord(**record) for pmid, record in sorted(payload["records"].items())
        },
    )


def write_ledger(path: Path, ledger: Ledger) -> None:
    """Write the ledger deterministically: sorted keys, one trailing newline.

    GP-30. Two refreshes that resolve the same records must produce the same bytes
    apart from ``resolved_utc``, or the diff is unreviewable.
    """
    payload = {
        "schema": ledger.schema,
        "resolved_utc": ledger.resolved_utc,
        "source": ledger.source,
        "records": {pmid: asdict(record) for pmid, record in sorted(ledger.records.items())},
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _year_matches(written: int, record: LedgerRecord) -> bool:
    """Accept either the print year or the electronic year.

    Not laxity. A 2013 epub with a 2014 print issue is cited both ways in the real
    literature and by different reference managers, and calling one of them a
    citation error would train a reader to ignore this check. Both years come from
    the same record, so a paper published in neither year still fails.
    """
    return written in {record.year, record.epub_year}


def compare(citations: list[WrittenCitation], ledger: Ledger) -> list[Mismatch]:
    """Every disagreement between what we wrote and what PubMed returned.

    GP-14 is load-bearing here. A citation written without a journal (a ``Baker
    2004`` back-reference) is not compared on its journal, because we made no
    journal claim there -- absence of a written journal is not a wrong journal.
    Only fields actually written are checked, and the form of every citation is
    recorded so that "we checked nothing" cannot look like "everything passed".
    """
    mismatches: list[Mismatch] = []
    cited: set[str] = set()

    for citation in citations:
        cited.add(citation.pmid)
        record = ledger.records.get(citation.pmid)
        if record is None:
            mismatches.append(
                Mismatch(
                    kind=MismatchKind.NOT_IN_LEDGER,
                    pmid=citation.pmid,
                    source=citation.source,
                    line=citation.line,
                    wrote=citation.context,
                    pubmed="(absent from the ledger)",
                    context=citation.context,
                )
            )
            continue

        if citation.first_author_surname is not None and normalise_surname(
            citation.first_author_surname
        ) != normalise_surname(record.first_author_surname):
            mismatches.append(
                Mismatch(
                    kind=MismatchKind.FIRST_AUTHOR,
                    pmid=citation.pmid,
                    source=citation.source,
                    line=citation.line,
                    wrote=citation.first_author_surname,
                    pubmed=f"{record.first_author} — {record.title}",
                    context=citation.context,
                )
            )

        if citation.year is not None and not _year_matches(citation.year, record):
            mismatches.append(
                Mismatch(
                    kind=MismatchKind.YEAR,
                    pmid=citation.pmid,
                    source=citation.source,
                    line=citation.line,
                    wrote=str(citation.year),
                    pubmed=f"pubdate {record.pubdate!r}, epubdate {record.epubdate!r}",
                    context=citation.context,
                )
            )

        if citation.title is not None and normalise_title(citation.title) != normalise_title(
            record.title
        ):
            mismatches.append(
                Mismatch(
                    kind=MismatchKind.TITLE,
                    pmid=citation.pmid,
                    source=citation.source,
                    line=citation.line,
                    wrote=citation.title,
                    pubmed=record.title,
                    context=citation.context,
                )
            )

        if citation.journal is not None and normalise_journal(
            citation.journal
        ) != normalise_journal(record.journal):
            mismatches.append(
                Mismatch(
                    kind=MismatchKind.JOURNAL,
                    pmid=citation.pmid,
                    source=citation.source,
                    line=citation.line,
                    wrote=citation.journal,
                    pubmed=record.journal,
                    context=citation.context,
                )
            )

    for pmid in sorted(set(ledger.records) - cited, key=int):
        record = ledger.records[pmid]
        mismatches.append(
            Mismatch(
                kind=MismatchKind.ORPHANED_IN_LEDGER,
                pmid=pmid,
                source=LEDGER_PATH,
                line=0,
                wrote="(nothing cites this PMID any more)",
                pubmed=f"{record.first_author} — {record.title}",
                context="",
            )
        )
    return mismatches


def format_mismatches(mismatches: list[Mismatch]) -> str:
    """A failure message that names what we wrote and what PubMed returned.

    Written for the reader of a red CI log, who has no context and needs the whole
    finding in the assertion text.
    """
    lines: list[str] = []
    for mismatch in mismatches:
        lines.append(
            f"  [{mismatch.kind.value}] PMID {mismatch.pmid} "
            f"({mismatch.source}:{mismatch.line})\n"
            f"      we wrote : {mismatch.wrote}\n"
            f"      PubMed   : {mismatch.pubmed}"
        )
        if mismatch.context:
            lines.append(f"      context  : ...{mismatch.context[-140:]}")
    return "\n".join(lines)


def form_counts(citations: list[WrittenCitation]) -> dict[str, int]:
    """How many citations were checked at each level of written detail."""
    counts = dict.fromkeys((form.value for form in CitationForm), 0)
    for citation in citations:
        counts[citation.form.value] += 1
    return counts
