"""The always-on, zero-network half of the citation check.

``docs/track2-hypothesis.md`` §0 says "every citation was resolved against PubMed
before being written down." This file is what turns that sentence from a claim into
a property. It never opens a socket: it compares our prose against
``knowledge/literature/citations.verified.json``, a committed snapshot of what NCBI
returned, so it runs under the same PRIV-05 offline rules as the patient-data path.

The staleness of that snapshot is the honest cost, and it is bounded by the opt-in
online test in ``tests/integration/test_citations_online.py``.

Two of these tests are positive controls rather than assertions about the document,
and they are the reason this file can be trusted. A citation checker that has
silently stopped parsing emits exactly the report a clean document does -- the same
indistinguishability :func:`mva.privacy.audit.check_detector_liveness` exists to
refuse. So one test corrupts a ledger record and requires the comparison to catch
it, and another pins the number of citations parsed at each level of written detail.
"""

from __future__ import annotations

import dataclasses
import json
import sys
from pathlib import Path

import pytest

from mva.config import find_repo_root

# `tools/` is not part of the installed `mva` distribution, so it is importable
# only with the repo root on `sys.path`. Same bootstrap, same reason, as
# `tests/unit/test_knowledge_tables.py`.
REPO_ROOT = find_repo_root(Path(__file__))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tools.citations.extract import (  # noqa: E402
    CitationForm,
    WrittenCitation,
    collect_written_citations,
    extract_written_citations,
    normalise_journal,
    normalise_surname,
    normalise_title,
    unique_pmids,
)
from tools.citations.ledger import (  # noqa: E402
    LEDGER_PATH,
    LEDGER_SCHEMA,
    Ledger,
    LedgerRecord,
    MismatchKind,
    compare,
    form_counts,
    format_mismatches,
    load_ledger,
    write_ledger,
)

#: The floor for each written form, checked by
#: :func:`test_citation_parse_coverage_has_not_regressed`. Raise these when the
#: document grows; never lower one to make a red test green, because lowering one
#: is exactly how this check would be switched off without anybody deciding to.
MINIMUM_FORM_COUNTS = {
    CitationForm.REFERENCE_ROW.value: 54,
    CitationForm.PROSE.value: 47,
    CitationForm.BARE.value: 5,
}


@pytest.fixture(scope="module")
def citations() -> list[WrittenCitation]:
    return collect_written_citations(REPO_ROOT)


@pytest.fixture(scope="module")
def ledger() -> Ledger:
    return load_ledger(REPO_ROOT / LEDGER_PATH)


@pytest.mark.unit
def test_every_cited_pmid_matches_the_verified_record(
    citations: list[WrittenCitation], ledger: Ledger
) -> None:
    """THE GATE. Every PMID resolves, and to the paper the sentence around it names.

    A PMID that resolves to a *different* paper is the failure this exists for: it
    passes every existence check, reads correctly to a human skimming the reference
    list, and is wrong. So author, year, journal and -- where we wrote one -- title
    are compared, not merely the identifier.
    """
    mismatches = compare(citations, ledger)
    assert not mismatches, (
        f"{len(mismatches)} citation(s) disagree with the verified PubMed record.\n"
        + format_mismatches(mismatches)
        + f"\n\nRemediation: if the document is wrong, fix the document. If the ledger "
        f"is stale, refresh it with `just verify-citations` (ONLINE) and commit "
        f"{LEDGER_PATH}. Do not edit the ledger by hand: it is a record of what NCBI "
        f"returned, and a hand-edited record verifies nothing."
    )


@pytest.mark.unit
def test_citation_parse_coverage_has_not_regressed(citations: list[WrittenCitation]) -> None:
    """A positive control for the parser, not for the document.

    If a regex change quietly demoted every citation to ``unannotated``, the gate
    above would still be green -- there would be no written claim left to
    contradict. That is the shape of a check that has stopped checking, and this
    test is what makes it visible.
    """
    counts = form_counts(citations)
    assert counts[CitationForm.UNANNOTATED.value] == 0, (
        "Some PMIDs have no parseable author beside them, so only their existence "
        "is verified:\n"
        + "\n".join(
            f"  {c.pmid} at {c.source}:{c.line} — ...{c.context[-120:]}"
            for c in citations
            if c.form is CitationForm.UNANNOTATED
        )
    )
    for form, floor in MINIMUM_FORM_COUNTS.items():
        assert counts[form] >= floor, (
            f"Citation parsing regressed: {counts[form]} citation(s) in {form!r} form, "
            f"expected at least {floor}. The document did not shrink; the parser did."
        )


@pytest.mark.unit
def test_comparison_catches_a_citation_pointing_at_the_wrong_paper(
    citations: list[WrittenCitation], ledger: Ledger
) -> None:
    """The canary. Repoint one ledger record at a different paper; demand a failure.

    Modelled on :func:`mva.privacy.patterns.detector_canaries`: every rule carries a
    specimen it MUST match, pushed through the function the real check calls rather
    than through a reimplementation of it. ``J Irreproducible Results`` and
    ``Nobody Q`` are fabricated and are not a journal or a person.
    """
    pmid = "15475955"
    assert pmid in ledger.records, "fixture PMID is no longer cited; pick another"
    corrupted = dict(ledger.records)
    corrupted[pmid] = dataclasses.replace(
        corrupted[pmid],
        first_author="Nobody Q",
        first_author_surname="Nobody",
        journal="J Irreproducible Results",
        pubdate="1999 Jan",
        epubdate="",
        year=1999,
        epub_year=1999,
        title="A fabricated title that is not the paper we cite",
    )
    mismatches = compare(citations, dataclasses.replace(ledger, records=corrupted))
    kinds = {m.kind for m in mismatches}
    assert MismatchKind.FIRST_AUTHOR in kinds
    assert MismatchKind.YEAR in kinds
    assert MismatchKind.JOURNAL in kinds
    assert MismatchKind.TITLE in kinds
    assert all(m.pmid == pmid for m in mismatches)
    rendered = format_mismatches(mismatches)
    assert "Nobody Q" in rendered and "Hanks" in rendered, (
        "The failure message must name BOTH what we wrote and what PubMed returned; "
        "a message that names only one leaves the reader unable to tell which is wrong."
    )


@pytest.mark.unit
def test_a_pmid_missing_from_the_ledger_is_a_failure_not_a_pass(
    citations: list[WrittenCitation], ledger: Ledger
) -> None:
    """GP-14, applied to the ledger. An unresolved PMID is not a verified one."""
    dropped = {k: v for k, v in ledger.records.items() if k != "15475955"}
    mismatches = compare(citations, dataclasses.replace(ledger, records=dropped))
    assert {m.kind for m in mismatches} == {MismatchKind.NOT_IN_LEDGER}


@pytest.mark.unit
def test_a_ledger_record_nothing_cites_is_reported(
    citations: list[WrittenCitation], ledger: Ledger
) -> None:
    """Orphans are surfaced, not deleted (GP-19). They usually mean a stale refresh."""
    extra = dict(ledger.records)
    extra["99999999"] = LedgerRecord(
        pmid="99999999",
        first_author="Canaryson A",
        first_author_surname="Canaryson",
        journal="J Irreproducible Results",
        pubdate="1999 Jan",
        epubdate="",
        year=1999,
        epub_year=1999,
        title="A fabricated record that nothing in this repository cites",
    )
    mismatches = compare(citations, dataclasses.replace(ledger, records=extra))
    assert [m.kind for m in mismatches] == [MismatchKind.ORPHANED_IN_LEDGER]


@pytest.mark.unit
def test_ledger_file_is_byte_identical_to_a_rewrite(ledger: Ledger, tmp_path: Path) -> None:
    """GP-30. A ledger that reformats on every refresh has an unreviewable diff."""
    rewritten = tmp_path / "citations.verified.json"
    write_ledger(rewritten, ledger)
    assert rewritten.read_bytes() == (REPO_ROOT / LEDGER_PATH).read_bytes()


@pytest.mark.unit
def test_ledger_refuses_an_unknown_schema(tmp_path: Path) -> None:
    path = tmp_path / "ledger.json"
    path.write_text(
        json.dumps({"schema": "something/else", "resolved_utc": "", "source": "", "records": {}}),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Refusing to read it"):
        load_ledger(path)


@pytest.mark.unit
def test_ledger_declares_when_it_was_resolved(ledger: Ledger) -> None:
    """A cache with no date is a cache nobody can judge the age of."""
    assert ledger.schema == LEDGER_SCHEMA
    assert ledger.resolved_utc.endswith("Z")
    assert "eutils.ncbi.nlm.nih.gov" in ledger.source


# --------------------------------------------------------------- parser unit tests
#
# Fabricated documents. `Canaryson`, `Q Fabricated` and `J Irreproducible Results`
# are not people or journals, and the PMIDs below are out of the assigned range.


@pytest.mark.unit
def test_parses_the_inline_prose_form() -> None:
    text = (
        "Checkpoint failure (Suijkerbuijk et al., *Cancer Res* 2010;70:4891-900 (PMID 20516114))."
    )
    (citation,) = extract_written_citations(text, source="x.md")
    assert citation.form is CitationForm.PROSE
    assert (citation.first_author_surname, citation.year, citation.journal) == (
        "Suijkerbuijk",
        2010,
        "Cancer Res",
    )


@pytest.mark.unit
def test_parses_a_multi_author_prose_form_and_takes_the_first_author() -> None:
    text = "Kops, Weaver & Cleveland, *Nat Rev Cancer* 2005;5:773-85 (PMID 16195750)."
    (citation,) = extract_written_citations(text, source="x.md")
    assert citation.first_author_surname == "Kops"


@pytest.mark.unit
def test_bold_markers_do_not_capture_the_journal() -> None:
    """The regression that made seven real citations look unparseable.

    ``**bold**`` shifts asterisk parity, so a naive italic regex pairs the wrong
    asterisks and swallows the journal name into a prose span.
    """
    text = (
        "**context_dependent** - a *precursor* and "
        "(Bitterman et al., *J Biol Chem* 2002;277:45099-107, PMID 12297502)"
    )
    (citation,) = extract_written_citations(text, source="x.md")
    assert (citation.first_author_surname, citation.journal) == ("Bitterman", "J Biol Chem")


@pytest.mark.unit
def test_a_citation_wrapped_across_lines_is_still_one_citation() -> None:
    text = "Kalish et al., *Clin Cancer Res*\n  2024;30:5260-5269, PMID 39320341"
    (citation,) = extract_written_citations(text, source="x.md")
    assert (citation.first_author_surname, citation.year, citation.journal) == (
        "Kalish",
        2024,
        "Clin Cancer Res",
    )
    # The line number points at the PMID, which is where a reader goes to fix it.
    assert citation.line == 2


@pytest.mark.unit
def test_two_citations_in_one_sentence_do_not_borrow_each_others_authors() -> None:
    text = (
        "Aa et al., *J Alpha* 2001;1:1-2 (PMID 10000001) and then "
        "Bb et al., *J Beta* 2002;2:3-4 (PMID 10000002)."
    )
    first, second = extract_written_citations(text, source="x.md")
    assert (first.first_author_surname, first.journal) == ("Aa", "J Alpha")
    assert (second.first_author_surname, second.journal) == ("Bb", "J Beta")


@pytest.mark.unit
def test_parses_the_reference_table_row_form_including_the_title() -> None:
    text = (
        "| R01 | Canaryson A, et al. (2004) A fabricated study of nothing at all. "
        "*J Irreproducible Results* 36(11):1159-61. PMID 10000003 |"
    )
    (citation,) = extract_written_citations(text, source="x.md")
    assert citation.form is CitationForm.REFERENCE_ROW
    assert citation.first_author_surname == "Canaryson"
    assert citation.year == 2004
    assert citation.journal == "J Irreproducible Results"
    assert citation.title == "A fabricated study of nothing at all"


@pytest.mark.unit
def test_a_back_reference_claims_no_journal_and_is_not_checked_for_one() -> None:
    """GP-14. ``Baker 2004`` asserts an author and a year, and nothing else."""
    (citation,) = extract_written_citations(
        "...deleting p16 attenuated it (Baker 2008, PMID 18516091)", source="x.md"
    )
    assert citation.form is CitationForm.BARE
    assert (citation.first_author_surname, citation.year, citation.journal) == (
        "Baker",
        2008,
        None,
    )


@pytest.mark.unit
def test_a_pmid_with_nothing_written_beside_it_is_recorded_as_unannotated() -> None:
    """Not skipped. A fourth category, so "checked nothing" cannot read as "passed"."""
    (citation,) = extract_written_citations("see also PMID 10000004.", source="x.md")
    assert citation.form is CitationForm.UNANNOTATED
    assert not citation.has_metadata


@pytest.mark.unit
@pytest.mark.parametrize(
    ("written", "record"),
    [
        ("García-Castillo", "Garcia-Castillo"),
        ("Lara-Gonzalez", "lara-gonzalez"),
        ("Méhes", "Mehes"),
    ],
)
def test_surname_normalisation_ignores_diacritics_and_case(written: str, record: str) -> None:
    assert normalise_surname(written) == normalise_surname(record)


@pytest.mark.unit
def test_surname_normalisation_does_not_conflate_different_people() -> None:
    assert normalise_surname("Baker") != normalise_surname("Bakker")


@pytest.mark.unit
def test_journal_normalisation_resolves_the_declared_abbreviations() -> None:
    assert normalise_journal("PNAS") == normalise_journal("Proc Natl Acad Sci U S A")
    assert normalise_journal("Proc Natl Acad Sci USA") == normalise_journal(
        "Proc Natl Acad Sci U S A"
    )
    assert normalise_journal("Nat Genet") != normalise_journal("Nat Cell Biol")


@pytest.mark.unit
def test_title_normalisation_ignores_typesetting_only() -> None:
    assert normalise_title("A study of NAD+ decline.") == normalise_title("A study of NAD+ decline")
    assert normalise_title("Alpha beta") != normalise_title("Alpha gamma")


@pytest.mark.unit
def test_every_source_file_contributes_or_is_deliberately_empty(
    citations: list[WrittenCitation],
) -> None:
    """Guards the sweep itself: a glob that silently matched nothing would pass."""
    assert len(unique_pmids(citations)) >= 54
    sources = {c.source for c in citations}
    assert "docs/track2-hypothesis.md" in sources
    assert "knowledge/literature/bub1b/SOURCES.md" in sources


@pytest.mark.unit
@pytest.mark.privacy
def test_the_offline_half_imports_no_network_client() -> None:
    """PRIV-05, structurally. The modules CI runs must not be able to reach out.

    An AST lint over imports, reusing the forbidden set that guards the
    patient-data path. `tools.citations.resolve` is deliberately excluded: it is
    the online half, it imports `urllib` at module scope so that a lint can see it,
    and nothing on this path imports it.
    """
    from tests.unit.test_architecture import forbidden_imports

    for module in ("extract.py", "ledger.py"):
        path = REPO_ROOT / "tools" / "citations" / module
        found = forbidden_imports(path.read_text(encoding="utf-8"))
        assert not found, f"tools/citations/{module} imports a network route: {found!r}"

    resolver = (REPO_ROOT / "tools" / "citations" / "resolve.py").read_text(encoding="utf-8")
    assert any(name == "urllib" for name, _ in forbidden_imports(resolver)), (
        "tools/citations/resolve.py must import urllib at module scope. Moving the "
        "import into a function body would hide the only network client in this "
        "package from every AST lint that looks for one."
    )
