"""The ONLINE half: resolve PMIDs against NCBI E-utilities.

Nothing on the always-on test path imports this module, and nothing on the
patient-data path may import it at all. It exists so that the offline half has a
committed record to check against; it is not itself part of any pipeline run.

PRIV-05. ``urllib`` is imported at module scope, not inside a function, precisely so
that ``tests/unit/test_architecture.py::test_no_network_clients_in_sensitive_stages``
can see it with an AST lint. Hiding a network client inside a function body to slip
past that lint would be the exact failure the lint is for. This module lives under
``tools/`` -- outside ``src/mva`` entirely -- for the same reason
``tools.acquire.fetch`` does: public reference data is fetched by a separate,
explicitly online step that never sees a proband coordinate.

What is sent: a comma-separated list of integers. Nothing else. The PMIDs of
published papers are not patient data, but the request is still worth stating in
full, because "we only send X" is checkable and "we are careful" is not.
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from collections.abc import Iterable, Sequence
from typing import Any, Final

from tools.citations.ledger import LedgerRecord

__all__ = [
    "BATCH_SIZE",
    "ESUMMARY_URL",
    "MIN_REQUEST_INTERVAL_S",
    "CitationResolutionError",
    "parse_esummary",
    "resolve_pmids",
    "surname_of",
]

ESUMMARY_URL: Final[str] = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

#: NCBI documents 3 requests/second without an API key. 0.34 s between requests
#: keeps us under it with one whole batch of headroom, and the whole job is a
#: handful of requests, so there is nothing to gain by going faster.
MIN_REQUEST_INTERVAL_S: Final[float] = 0.34

#: E-utilities accepts several hundred UIDs per GET. 100 keeps the URL short enough
#: to stay well inside any proxy's line limit and makes a partial failure cheap to
#: attribute.
BATCH_SIZE: Final[int] = 100

_REQUEST_TIMEOUT_S: Final[float] = 30.0

#: Sent so NCBI can identify the traffic, as their usage policy asks.
_USER_AGENT: Final[str] = "mva-research-citation-check/0.1 (+https://github.com/mva-research)"


class CitationResolutionError(RuntimeError):
    """A PMID batch could not be resolved. Never swallowed: a silent skip here
    would let an unverifiable citation pass as a verified one."""


def surname_of(author: str) -> str:
    """``"Lara-Gonzalez P"`` -> ``"Lara-Gonzalez"``.

    PubMed writes ``Surname Initials``. The initials are a trailing run of one to
    three capitals; everything before them is the surname, particles included, so
    ``van der Meer J`` keeps all three words.
    """
    parts = author.split()
    if len(parts) > 1 and 1 <= len(parts[-1]) <= 3 and parts[-1].isupper():
        parts = parts[:-1]
    return " ".join(parts)


def _year_of(pubdate: str) -> int | None:
    token = pubdate.strip()[:4]
    return int(token) if token.isdigit() else None


def parse_esummary(payload: dict[str, Any], requested: Sequence[str]) -> dict[str, LedgerRecord]:
    """Turn one esummary response into records, failing loudly on anything odd.

    A PMID that does not exist comes back as a summary carrying an ``error`` key
    rather than as an HTTP failure, so the absent case has to be detected here or
    it is not detected at all.
    """
    result = payload.get("result")
    if not isinstance(result, dict):
        msg = f"esummary response has no 'result' object for {list(requested)!r}."
        raise CitationResolutionError(msg)

    records: dict[str, LedgerRecord] = {}
    for pmid in requested:
        summary = result.get(pmid)
        if not isinstance(summary, dict):
            msg = (
                f"PMID {pmid} is absent from the esummary response. It was requested "
                f"and NCBI returned nothing for it: the identifier does not resolve."
            )
            raise CitationResolutionError(msg)
        if summary.get("error"):
            msg = f"PMID {pmid} does not resolve: NCBI returned error {summary['error']!r}."
            raise CitationResolutionError(msg)

        authors = summary.get("authors") or []
        first = ""
        for entry in authors:
            if isinstance(entry, dict) and entry.get("authtype") == "Author":
                first = str(entry.get("name", ""))
                break
        if not first and authors and isinstance(authors[0], dict):
            first = str(authors[0].get("name", ""))

        records[pmid] = LedgerRecord(
            pmid=pmid,
            first_author=first,
            first_author_surname=surname_of(first),
            journal=str(summary.get("source", "")),
            pubdate=str(summary.get("pubdate", "")),
            epubdate=str(summary.get("epubdate", "")),
            year=_year_of(str(summary.get("pubdate", ""))),
            epub_year=_year_of(str(summary.get("epubdate", ""))),
            title=str(summary.get("title", "")).strip(),
        )
    return records


def _fetch(url: str) -> dict[str, Any]:
    request = urllib.request.Request(  # noqa: S310 - fixed https host, see ESUMMARY_URL
        url, method="GET", headers={"User-Agent": _USER_AGENT}
    )
    with urllib.request.urlopen(request, timeout=_REQUEST_TIMEOUT_S) as response:  # noqa: S310
        body = response.read()
    parsed = json.loads(body.decode("utf-8"))
    if not isinstance(parsed, dict):
        msg = "esummary returned a non-object JSON body."
        raise CitationResolutionError(msg)
    return parsed


def resolve_pmids(
    pmids: Iterable[str],
    *,
    batch_size: int = BATCH_SIZE,
    sleep: float = MIN_REQUEST_INTERVAL_S,
) -> dict[str, LedgerRecord]:
    """Resolve every PMID, batched and rate-limited. Raises on the first failure."""
    ordered = sorted({str(p) for p in pmids}, key=int)
    records: dict[str, LedgerRecord] = {}
    for index in range(0, len(ordered), batch_size):
        batch = ordered[index : index + batch_size]
        query = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(batch), "retmode": "json"})
        if index:
            time.sleep(sleep)
        try:
            payload = _fetch(f"{ESUMMARY_URL}?{query}")
        except CitationResolutionError:
            raise
        except OSError as exc:
            msg = f"Could not reach NCBI E-utilities for PMIDs {batch[0]}..{batch[-1]}: {exc!r}"
            raise CitationResolutionError(msg) from exc
        records.update(parse_esummary(payload, batch))
    return records
