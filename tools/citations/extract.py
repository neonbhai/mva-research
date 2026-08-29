"""Pull every PMID -- and the author, year and journal written beside it -- out of our own prose.

Offline by construction. This module imports nothing that can reach a network and
is the half of the citation check that runs on every ``pytest`` invocation.

Why parse prose at all, rather than keeping a hand-maintained citation database?
Because the failure this exists to catch is a *disagreement between the two*: a
PMID that resolves perfectly well but resolves to a different paper than the
sentence around it describes. Checking a database against PubMed proves the
database is internally consistent and proves nothing about the document a reviewer
actually reads. So the document is the input.

Three written forms appear in this repository, and each PMID is classified into
exactly one of them:

``REFERENCE_ROW``
    ``| R01 | Hanks S, et al. (2004) Title. *Nat Genet* 36(11):1159-61. PMID 15475955 |``
    -- the ``knowledge/literature/bub1b/SOURCES.md`` table. Year in parentheses
    *before* the journal, and a **title**, which is the single most discriminating
    field there is: two papers can share an author, a year and a journal, and no
    two share a title.

``PROSE``
    ``Suijkerbuijk et al., *Cancer Res* 2010;70:4891-900 (PMID 20516114)``
    -- the inline form used throughout ``docs/track2-hypothesis.md``. Year *after*
    the journal.

``BARE``
    ``Baker 2004 (PMID 15208629)`` -- a back-reference to a citation given in full
    elsewhere. Surname and year are checkable; no journal was written.

``UNANNOTATED``
    A PMID with no parseable surname beside it. Its existence is still checked; its
    metadata is not, because there is no written claim to contradict. Recording it
    as a *fourth category* rather than silently skipping it is GP-14: absence of a
    written author is not evidence that the author matches.

The category counts are asserted by the test suite. A regex regression that quietly
demoted every citation to ``UNANNOTATED`` would otherwise turn this file into a
checker that checks nothing and still passes -- the same failure mode
:func:`mva.privacy.audit.check_detector_liveness` exists to refuse.
"""

from __future__ import annotations

import re
import unicodedata
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Final

__all__ = [
    "DOCUMENT_SOURCES",
    "JOURNAL_ABBREVIATION_ALIASES",
    "CitationForm",
    "WrittenCitation",
    "collect_written_citations",
    "extract_written_citations",
    "normalise_journal",
    "normalise_surname",
    "normalise_title",
    "unique_pmids",
]

#: Files whose PMIDs are verified. Relative to the repository root.
#:
#: ``knowledge/literature/`` is swept whole (markdown *and* TSV) rather than named
#: file by file, so a new literature table is covered the day it lands instead of
#: the day someone remembers to add it here.
#:
#: NOT in scope, deliberately: ``knowledge/disease/mva_panel.tsv`` carries ~350
#: gene-disease PMIDs with no written author/year/journal beside them, so only the
#: existence half of this check could ever apply to it, at ~10x the request budget.
#: That is a real uncovered surface and it is named here rather than left implicit.
#: Every file whose written citations are checked against PubMed.
#:
#: A document left off this list is not "unchecked but fine" — it is unchecked and
#: silent about it, which is the exact shape of the failure this whole tool exists
#: to catch. The claims ledger and the quantified-uncertainty note were both added
#: here after review found they carried PMIDs no gate had ever looked at.
#:
#: Adding a file is cheap. Leaving one off is invisible. Prefer over-inclusion.
DOCUMENT_SOURCES: Final[tuple[str, ...]] = (
    "docs/track2-hypothesis.md",
    "docs/track2-claims-ledger.md",
    "docs/track2-quantified-uncertainty.md",
    "docs/track2-drugs-to-avoid.md",
    "knowledge/literature/",
)

_SOURCE_SUFFIXES: Final[frozenset[str]] = frozenset({".md", ".tsv"})


class CitationForm(StrEnum):
    """How much of a citation was written next to a PMID."""

    REFERENCE_ROW = "reference_row"
    PROSE = "prose"
    BARE = "bare"
    UNANNOTATED = "unannotated"


@dataclass(frozen=True, slots=True)
class WrittenCitation:
    """One PMID as it appears in one of our files, with whatever we claimed about it."""

    pmid: str
    form: CitationForm
    first_author_surname: str | None
    year: int | None
    journal: str | None
    title: str | None
    source: str
    line: int
    context: str

    @property
    def has_metadata(self) -> bool:
        """True when there is a written claim to contradict."""
        return self.first_author_surname is not None


# --------------------------------------------------------------------------- regexes

_PMID = re.compile(r"PMID[:\s]\s*(?P<pmid>\d{6,9})")

#: A markdown italic span. Journal names are the only italics that ever sit
#: immediately before a citation's volume/page string.
#:
#: The ``(?<!\*)``/``(?!\*)`` guards make a lone asterisk the delimiter, so that
#: ``**bold**`` is skipped rather than consumed. Without them the asterisk parity of
#: a sentence like ``**context_dependent** ... *and* a ... (Bitterman et al., *J Biol
#: Chem* 2002`` shifts by one and the journal span is silently mispaired -- which is
#: how seven correctly written citations first came back as ``UNANNOTATED``.
_ITALIC = re.compile(r"(?<!\*)\*(?!\*)(?P<journal>[^*\n]{2,200})(?<!\*)\*(?!\*)")

#: Longest string accepted *as a journal name*. The span regex above scans much
#: wider than this on purpose: a cap applied during scanning does not reject a long
#: italic, it mispairs the asterisks after it and swallows the real journal. So the
#: scan is generous and the acceptance is strict.
_MAX_JOURNAL_CHARS: Final[int] = 60

#: An alphabetic run long enough to be a word. Used to reject an italic span that
#: is prose rather than a journal: between a journal name and its PMID there is
#: only a volume/page/year string, whose only letters are fragments -- ``e5``,
#: ``ra15``, ``Pt``, the ``E`` of ``E4134-41``.
_WORDISH = re.compile(r"[A-Za-z]{3,}")

#: An author list ending immediately before ``, *Journal*``. Anchored at the end so
#: that leftmost-match cannot start inside the preceding sentence.
#: ``\u2019`` rather than a literal typographic apostrophe: the ``re`` engine
#: understands the escape, and the source file stays free of characters that look
#: like ASCII and are not.
_AUTHORS_TAIL = re.compile(
    r"""
    (?P<first>[A-ZÀ-Þ][A-Za-zÀ-ɏ'\u2019-]{1,30})
    (?:\s+[A-Z]{1,3}\b)?
    (?:
        \s*(?:,|&|and|et\s+al\.?)\s*
        (?:[A-ZÀ-Þ][A-Za-zÀ-ɏ'\u2019-]{1,30}(?:\s+[A-Z]{1,3}\b)?)?
    )*
    \s*,?\s*$
    """,
    re.VERBOSE,
)

#: ``(2004)`` -- the reference-table year, which precedes the journal.
_PARENTHESISED_YEAR = re.compile(r"\((?P<year>(?:19|20)\d{2})\)")

#: ``2010;70:4891-900`` -- the prose year, which follows the journal.
_TRAILING_YEAR = re.compile(r"(?<!\d)(?P<year>(?:19|20)\d{2})(?!\d)")

#: The opening of a reference-table cell: ``Hanks S, et al. (2004)``.
_ROW_FIRST_AUTHOR = re.compile(
    r"^[\s|]*(?:R\d+\s*\|\s*)?(?P<first>[A-ZÀ-Þ][A-Za-zÀ-ɏ'\u2019-]{1,30})"
)

#: ``Baker 2004 (`` / ``Baker 2011,`` -- a back-reference with no journal.
_BARE = re.compile(
    r"(?P<first>[A-ZÀ-Þ][A-Za-zÀ-ɏ'\u2019-]{1,30})"
    r"\s+(?P<year>(?:19|20)\d{2})\s*[,;(]?\s*$"
)

#: Journal names we abbreviate differently from PubMed's ``source`` field.
#:
#: Kept explicit and tiny on purpose. Every entry is a place where a human chose a
#: shorter name than the record carries; a fuzzy matcher here would also swallow a
#: genuinely wrong journal, which is the thing being looked for.
JOURNAL_ABBREVIATION_ALIASES: Final[dict[str, str]] = {
    "pnas": "Proc Natl Acad Sci U S A",
}

_MAX_WINDOW = 320


def normalise_surname(surname: str) -> str:
    """Casefold, strip diacritics, keep letters and internal hyphens.

    ``García-Castillo`` and ``Garcia-Castillo`` are the same author written by
    two people; ``Baker`` and ``Bakker`` are not, and nothing here conflates them.
    """
    decomposed = unicodedata.normalize("NFKD", surname)
    stripped = "".join(ch for ch in decomposed if not unicodedata.combining(ch))
    kept = [ch for ch in stripped.casefold() if ch.isalpha() or ch == "-"]
    return "".join(kept).strip("-")


def normalise_title(title: str) -> str:
    """Casefold, strip diacritics, drop everything non-alphanumeric.

    Markdown emphasis, a trailing full stop, ``NAD+`` versus ``NAD(+)`` and British
    versus American spelling of a Greek-derived word are all differences in how a
    title is *typeset*. A different paper is a different sequence of words.
    """
    decomposed = unicodedata.normalize("NFKD", title)
    stripped = "".join(ch for ch in decomposed if not unicodedata.combining(ch))
    return "".join(ch for ch in stripped.casefold() if ch.isalnum())


def normalise_journal(journal: str) -> str:
    """Casefold and drop everything that is not alphanumeric.

    ``Proc Natl Acad Sci USA`` and ``Proc Natl Acad Sci U S A`` differ only in
    spacing, and that difference is not a citation error.
    """
    alias = JOURNAL_ABBREVIATION_ALIASES.get(
        "".join(ch for ch in journal.casefold() if ch.isalnum())
    )
    resolved = alias if alias is not None else journal
    return "".join(ch for ch in resolved.casefold() if ch.isalnum())


def _window(flat: str, start: int, previous_end: int) -> tuple[str, int, bool]:
    """The text a citation could plausibly occupy, where it begins, and how it began.

    Bounded three ways, tightest wins: the end of the previous PMID (two citations
    in one sentence must not borrow each other's authors), a fixed character budget,
    and the last markdown table pipe (a table cell is a hard citation boundary).

    The third return value says whether the window starts at a *pipe* rather than at
    the character budget. It matters: the reference-row parse anchors the first
    author to the start of the window with ``^``, and that anchor is only meaningful
    if the start is a real boundary. A window truncated mid-sentence by the budget
    would make ``^`` point at an arbitrary word, so the row form is refused there.
    """
    begin = max(previous_end, start - _MAX_WINDOW, 0)
    text = flat[begin:start]
    pipe = text.rfind("|")
    if pipe != -1:
        begin += pipe + 1
        return flat[begin:start], begin, True
    return text, begin, False


def _journal_in(window: str) -> re.Match[str] | None:
    """The rightmost italic span that is followed only by a volume/page string."""
    for match in reversed(list(_ITALIC.finditer(window))):
        if len(match.group("journal")) > _MAX_JOURNAL_CHARS:
            continue
        if not _WORDISH.search(window[match.end() :]):
            return match
    return None


@dataclass(frozen=True, slots=True)
class _Parsed:
    form: CitationForm
    surname: str | None
    year: int | None
    journal: str | None
    title: str | None


_UNPARSED = _Parsed(CitationForm.UNANNOTATED, None, None, None, None)


def _parse_window(window: str, *, at_cell_boundary: bool) -> _Parsed:
    journal_match = _journal_in(window)
    if journal_match is not None:
        head = window[: journal_match.start()]
        tail = window[journal_match.end() :]
        journal = " ".join(journal_match.group("journal").split())
        trailing = _TRAILING_YEAR.search(tail)
        if trailing is not None:
            authors = _AUTHORS_TAIL.search(head.rstrip().rstrip(","))
            if authors is not None:
                return _Parsed(
                    CitationForm.PROSE,
                    authors.group("first"),
                    int(trailing.group("year")),
                    journal,
                    None,
                )
        parenthesised = list(_PARENTHESISED_YEAR.finditer(head))
        if at_cell_boundary and parenthesised:
            row_author = _ROW_FIRST_AUTHOR.match(head)
            if row_author is not None:
                title = head[parenthesised[-1].end() :].strip().strip("*").strip()
                return _Parsed(
                    CitationForm.REFERENCE_ROW,
                    row_author.group("first"),
                    int(parenthesised[-1].group("year")),
                    journal,
                    " ".join(title.rstrip(".").split()) or None,
                )
    bare = _BARE.search(window)
    if bare is not None:
        return _Parsed(CitationForm.BARE, bare.group("first"), int(bare.group("year")), None, None)
    return _UNPARSED


def extract_written_citations(text: str, *, source: str) -> list[WrittenCitation]:
    """Every PMID in one file, with whatever was written beside it.

    Newlines are replaced by spaces rather than collapsed, so offsets into the
    flattened text still index the original and line numbers stay exact. A citation
    wrapped across two lines is then just a citation.
    """
    flat = text.replace("\n", " ")
    citations: list[WrittenCitation] = []
    previous_end = 0
    for match in _PMID.finditer(flat):
        window, begin, at_cell_boundary = _window(flat, match.start(), previous_end)
        previous_end = match.end()
        parsed = _parse_window(window, at_cell_boundary=at_cell_boundary)
        citations.append(
            WrittenCitation(
                pmid=match.group("pmid"),
                form=parsed.form,
                first_author_surname=parsed.surname,
                year=parsed.year,
                journal=parsed.journal,
                title=parsed.title,
                source=source,
                line=text.count("\n", 0, match.start()) + 1,
                context=" ".join(flat[begin : match.end()].split())[-200:],
            )
        )
    return citations


def _source_files(repo_root: Path) -> list[Path]:
    files: list[Path] = []
    for entry in DOCUMENT_SOURCES:
        target = repo_root / entry
        if target.is_dir():
            files.extend(
                path
                for path in sorted(target.rglob("*"))
                if path.is_file() and path.suffix in _SOURCE_SUFFIXES
            )
        elif target.is_file():
            files.append(target)
    return files


def collect_written_citations(repo_root: Path) -> list[WrittenCitation]:
    """Every citation in :data:`DOCUMENT_SOURCES`, sorted for a stable diff."""
    citations: list[WrittenCitation] = []
    for path in _source_files(repo_root):
        citations.extend(
            extract_written_citations(
                path.read_text(encoding="utf-8"),
                source=path.relative_to(repo_root).as_posix(),
            )
        )
    return sorted(citations, key=lambda c: (c.pmid, c.source, c.line))


def unique_pmids(citations: list[WrittenCitation]) -> list[str]:
    """Deduplicated PMIDs in numeric order -- the batch list for a resolver."""
    return sorted({c.pmid for c in citations}, key=int)
