"""The opt-in, ONLINE half of the citation check.

Skipped by default. Runs only under ``pytest -m network`` or with
``MVA_VERIFY_CITATIONS=1`` set (see ``tests/conftest.py``), and never on the
patient-data path.

What it adds over the always-on offline test: the offline test proves our prose
agrees with a *snapshot*, and a snapshot can be stale, or wrong in the same way the
document is wrong if someone hand-edited both. This one asks NCBI. It is the thing
that would actually catch a hallucinated PMID the first time it is run after the
citation is written.

Two assertions, in this order, because they fail differently:

1. every cited PMID resolves at all -- an unresolvable identifier raises out of
   :func:`tools.citations.resolve.resolve_pmids` rather than being skipped;
2. the live records equal the committed ledger, field for field. A difference here
   is not necessarily an error -- PubMed does correct records -- but it means the
   committed file no longer says what NCBI says, and nobody should find that out by
   accident.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from mva.config import find_repo_root

REPO_ROOT = find_repo_root(Path(__file__))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools.citations.extract import collect_written_citations, unique_pmids  # noqa: E402
from tools.citations.ledger import (  # noqa: E402
    LEDGER_PATH,
    Ledger,
    LedgerRecord,
    compare,
    format_mismatches,
    load_ledger,
)

pytestmark = [pytest.mark.network, pytest.mark.integration, pytest.mark.slow]


@pytest.fixture(scope="module")
def live_records() -> dict[str, LedgerRecord]:
    from tools.citations.resolve import resolve_pmids

    citations = collect_written_citations(REPO_ROOT)
    return resolve_pmids(unique_pmids(citations))


def test_every_cited_pmid_resolves_against_pubmed(live_records: dict[str, LedgerRecord]) -> None:
    """Existence, live. ``resolve_pmids`` raises on a PMID NCBI does not know."""
    citations = collect_written_citations(REPO_ROOT)
    expected = set(unique_pmids(citations))
    assert set(live_records) == expected
    assert all(record.first_author for record in live_records.values()), (
        "A resolved record with no author is a record we cannot check a citation "
        "against, and must not be treated as verified."
    )


def test_written_citations_match_the_live_records(live_records: dict[str, LedgerRecord]) -> None:
    """The real check: author, year, journal and title, against PubMed itself."""
    citations = collect_written_citations(REPO_ROOT)
    ledger = load_ledger(REPO_ROOT / LEDGER_PATH)
    live = Ledger(
        schema=ledger.schema,
        resolved_utc="(live)",
        source=ledger.source,
        records=live_records,
    )
    mismatches = compare(citations, live)
    assert not mismatches, (
        f"{len(mismatches)} citation(s) disagree with PubMed RIGHT NOW.\n"
        + format_mismatches(mismatches)
        + "\n\nThis is the finding this test exists for. Fix the document."
    )


def test_the_committed_ledger_still_matches_pubmed(live_records: dict[str, LedgerRecord]) -> None:
    """Staleness detector. Refresh and commit if this fails for a benign reason."""
    ledger = load_ledger(REPO_ROOT / LEDGER_PATH)
    drifted = [
        f"  PMID {pmid}\n    ledger: {ledger.records.get(pmid)}\n    pubmed: {record}"
        for pmid, record in sorted(live_records.items())
        if ledger.records.get(pmid) != record
    ]
    assert not drifted, (
        f"{LEDGER_PATH} no longer matches PubMed:\n"
        + "\n".join(drifted)
        + "\n\nRemediation: `just verify-citations` rewrites the ledger; review the "
        "diff before committing it. PubMed does correct records, so a difference is "
        "not automatically an error -- but it is always something a human should see."
    )
