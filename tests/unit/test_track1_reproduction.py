"""Comparing two Track 1 submissions without either of them becoming readable.

The comparison itself is arithmetic over row keys and is easy to get right. The
property that is *not* easy — and that the second half of this module exists to
prove — is that nothing the comparator emits is derived from the bytes it read.
Both inputs are lists of one child's variant coordinates (GP-40), the caller is
an agent, and its output goes straight into that agent's context (PRIV-03,
ADR 0008). So a leak here is not a cosmetic defect; it is the disclosure the
whole privacy model is built to prevent.

The leak tests therefore work the way the PRIV-05 import lint does: rather than
asserting that a tree which already passes still passes, they feed the tool input
constructed to be recognisable — every one of the twelve columns filled with a
token no legitimate output could contain — and then assert those tokens appear
nowhere in stdout, stderr, the exception, or the report object. Each case also
runs a positive control: the tokens must be present in the fixture, or the test
would pass for the uninteresting reason that there was nothing to leak.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path

import pytest
from typer.testing import CliRunner

from mva.cli import app
from mva.errors import SubmissionFormatError
from mva.reporting.reproduction import (
    TOP_K_TIERS,
    VERDICT_MEANINGS,
    ReproductionVerdict,
    compare_rankings,
    compare_submission_files,
    read_ranked_keys,
    spearman_over_intersection,
)
from mva.reporting.track1 import MAX_SUBMISSION_ROWS, TRACK1_COLUMNS

pytestmark = [pytest.mark.unit, pytest.mark.privacy]

Variant = tuple[str, int, str, str]

#: Synthetic coordinates. `SYNTH`-style gene names are not available in a CSV
#: that has no gene column, so distinctness is carried by position instead; the
#: contig is a real one because the scorer compares contig strings raw and a
#: fixture using a fictional contig would not exercise the same code path.
_CONTIG = "chr15"


def _variant(ordinal: int) -> Variant:
    return (_CONTIG, 40_200_000 + ordinal, "C", "T")


def _row(variants: Sequence[Variant], epcr: float) -> list[str]:
    """One contract row. `epcr` is what ranks it, so callers set it explicitly."""
    first = variants[0]
    second = variants[1] if len(variants) > 1 else ("", "", "", "")
    return [
        "PROBAND01",
        *[str(cell) for cell in first],
        *[str(cell) for cell in second],
        f"{epcr:.4f}",
        "primary",
        "",
    ]


def _csv(rows: Sequence[Sequence[str]], *, header: Sequence[str] = TRACK1_COLUMNS) -> str:
    return "\n".join([",".join(header), *[",".join(row) for row in rows]]) + "\n"


def _ranking(*ordinals: int, descending: bool = True) -> str:
    """A submission proposing one single-variant row per ordinal, in that order.

    EPCRs are laid out so the file's own order IS the scored order, which is what
    a rendered submission looks like; `descending=False` inverts them to exercise
    the fact that rank comes from EPCR rather than from file position.
    """
    step = -0.01 if descending else 0.01
    start = 0.9 if descending else 0.9 - step * (len(ordinals) - 1)
    return _csv([_row([_variant(o)], start + step * i) for i, o in enumerate(ordinals)])


def _write(directory: Path, name: str, text: str) -> Path:
    path = directory / name
    path.write_text(text, encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# The statistics
# ---------------------------------------------------------------------------


def test_identical_files_reproduce(tmp_path: Path) -> None:
    text = _ranking(1, 2, 3)
    report = compare_submission_files(
        _write(tmp_path, "submitted.csv", text), _write(tmp_path, "produced.csv", text)
    )

    assert report.verdict is ReproductionVerdict.IDENTICAL
    assert report.reproduced
    assert (report.submitted_rows, report.produced_rows) == (3, 3)
    assert (report.intersection, report.submitted_only, report.produced_only) == (3, 0, 0)
    assert report.top1_identical
    assert report.spearman == 1.0
    assert [a.agreeing for a in report.top_k] == [1, 3, 3, 3]
    # The windows shrink with the file, so "3 of 3" at k=5 is not a claim that
    # five rows agreed (GP-14).
    assert [a.produced_window for a in report.top_k] == [1, 3, 3, 3]


def test_reordering_the_same_rows_is_not_a_reproduction(tmp_path: Path) -> None:
    """Same candidates, different ranking. The set counts cannot see this."""
    report = compare_submission_files(
        _write(tmp_path, "submitted.csv", _ranking(1, 2, 3)),
        _write(tmp_path, "produced.csv", _ranking(3, 2, 1)),
    )

    assert report.verdict is ReproductionVerdict.SAME_SET_DIFFERENT_ORDER
    assert not report.reproduced
    assert (report.intersection, report.submitted_only, report.produced_only) == (3, 0, 0)
    assert not report.top1_identical
    assert report.spearman == -1.0


def test_a_produced_superset_is_reported_as_a_superset(tmp_path: Path) -> None:
    report = compare_submission_files(
        _write(tmp_path, "submitted.csv", _ranking(1, 2)),
        _write(tmp_path, "produced.csv", _ranking(1, 2, 3, 4)),
    )

    assert report.verdict is ReproductionVerdict.SUPERSET
    assert (report.intersection, report.submitted_only, report.produced_only) == (2, 0, 2)
    assert report.top1_identical


def test_a_produced_subset_is_reported_as_a_subset(tmp_path: Path) -> None:
    report = compare_submission_files(
        _write(tmp_path, "submitted.csv", _ranking(1, 2, 3, 4)),
        _write(tmp_path, "produced.csv", _ranking(1, 2)),
    )

    assert report.verdict is ReproductionVerdict.SUBSET
    assert (report.intersection, report.submitted_only, report.produced_only) == (2, 2, 0)


def test_disjoint_files_are_divergent(tmp_path: Path) -> None:
    report = compare_submission_files(
        _write(tmp_path, "submitted.csv", _ranking(1, 2, 3)),
        _write(tmp_path, "produced.csv", _ranking(7, 8, 9)),
    )

    assert report.verdict is ReproductionVerdict.DIVERGENT
    assert (report.intersection, report.submitted_only, report.produced_only) == (0, 3, 3)
    assert not report.top1_identical
    assert all(agreement.agreeing == 0 for agreement in report.top_k)
    assert report.spearman is None


def test_partial_overlap_in_both_directions_is_divergent(tmp_path: Path) -> None:
    """Neither containment holds, so neither `superset` nor `subset` is honest."""
    report = compare_submission_files(
        _write(tmp_path, "submitted.csv", _ranking(1, 2, 3)),
        _write(tmp_path, "produced.csv", _ranking(2, 3, 4)),
    )

    assert report.verdict is ReproductionVerdict.DIVERGENT
    assert (report.intersection, report.submitted_only, report.produced_only) == (2, 1, 1)


def test_two_empty_rankings_are_identical() -> None:
    """Vacuously true, and reported as such rather than as a failure to compare."""
    report = compare_rankings((), ())

    assert report.verdict is ReproductionVerdict.IDENTICAL
    assert not report.top1_identical  # there is no rank-1 row to be identical to
    assert report.spearman is None


def test_spearman_is_none_below_two_shared_rows() -> None:
    """GP-14: an undefined statistic is reported as undefined, not as zero.

    Returning 0.0 would state "the two rankings are uncorrelated", which is a
    measurement. One shared row admits no measurement at all — there is nothing
    for it to be in the wrong order relative to.
    """
    shared_one = compare_rankings(
        (frozenset({_variant(1)}), frozenset({_variant(2)})),
        (frozenset({_variant(1)}), frozenset({_variant(9)})),
    )

    assert shared_one.intersection == 1
    assert shared_one.spearman is None


@pytest.mark.parametrize(
    ("submitted", "produced", "expected"),
    [
        # Perfect agreement, then perfect reversal: the two ends of the range.
        ((1, 2, 3, 4), (1, 2, 3, 4), 1.0),
        ((1, 2, 3, 4), (4, 3, 2, 1), -1.0),
        # One adjacent swap at the bottom: d^2 = 1 + 1 over n = 4, so
        # rho = 1 - 6*2 / (4 * 15) = 0.8.
        ((1, 2, 3, 4), (1, 2, 4, 3), 0.8),
        # Rows the two files do not share are excluded and the rest re-ranked, so
        # this is the swap above with two strangers added on either side.
        ((1, 2, 3, 4, 8), (7, 1, 2, 4, 3), 0.8),
    ],
)
def test_spearman_matches_the_closed_form(
    submitted: tuple[int, ...], produced: tuple[int, ...], expected: float
) -> None:
    rho = spearman_over_intersection(
        tuple(frozenset({_variant(o)}) for o in submitted),
        tuple(frozenset({_variant(o)}) for o in produced),
    )

    assert rho is not None
    assert rho == pytest.approx(expected)


def test_rank_order_comes_from_epcr_not_file_order(tmp_path: Path) -> None:
    """The scorer sorts by EPCR, so a file that is not pre-sorted still ranks.

    Comparing file order directly would report a difference the challenge cannot
    see, and would call a correct reproduction divergent.
    """
    descending = _write(tmp_path, "submitted.csv", _ranking(1, 2, 3))
    ascending = _write(tmp_path, "produced.csv", _ranking(3, 2, 1, descending=False))

    assert read_ranked_keys(descending, role="submitted") == read_ranked_keys(
        ascending, role="produced"
    )
    assert compare_submission_files(descending, ascending).reproduced


def test_a_pair_row_is_one_key_whichever_way_round_it_is_written(tmp_path: Path) -> None:
    """Rows are matched by frozenset equality, exactly as the scorer matches them."""
    forward = _csv([_row([_variant(1), _variant(2)], 0.9)])
    reversed_pair = _csv([_row([_variant(2), _variant(1)], 0.9)])

    report = compare_submission_files(
        _write(tmp_path, "submitted.csv", forward),
        _write(tmp_path, "produced.csv", reversed_pair),
    )

    assert report.reproduced


def test_a_pair_and_its_two_halves_are_different_rows(tmp_path: Path) -> None:
    """A split pair is a real divergence, not a formatting difference (ADR 0023)."""
    report = compare_submission_files(
        _write(tmp_path, "submitted.csv", _csv([_row([_variant(1), _variant(2)], 0.9)])),
        _write(tmp_path, "produced.csv", _ranking(1, 2)),
    )

    assert report.verdict is ReproductionVerdict.DIVERGENT
    assert report.intersection == 0


def test_comparison_is_deterministic(tmp_path: Path) -> None:
    """GP-30: no set-ordering or hash-ordering reaches the result."""
    submitted = _write(tmp_path, "submitted.csv", _ranking(1, 2, 3, 4))
    produced = _write(tmp_path, "produced.csv", _ranking(2, 1, 5, 3))

    first = compare_submission_files(submitted, produced)
    second = compare_submission_files(submitted, produced)

    assert first == second


def test_every_verdict_has_a_meaning() -> None:
    """A verdict word with no gloss is a token, not an answer."""
    assert set(VERDICT_MEANINGS) == set(ReproductionVerdict)
    assert all(len(meaning) > 20 for meaning in VERDICT_MEANINGS.values())


def test_top_k_tiers_are_the_scorers_rank_point_boundaries() -> None:
    """If the contract's tiers move, these k values are wrong and must move too."""
    assert TOP_K_TIERS == (1, 3, 5, 10)
    assert TOP_K_TIERS[-1] == MAX_SUBMISSION_ROWS


# ---------------------------------------------------------------------------
# Malformed input
# ---------------------------------------------------------------------------


def test_an_empty_file_is_refused(tmp_path: Path) -> None:
    path = _write(tmp_path, "produced.csv", "")

    with pytest.raises(SubmissionFormatError) as excinfo:
        read_ranked_keys(path, role="produced")

    assert "no header row" in str(excinfo.value)


def test_a_header_only_file_is_read_as_zero_rows(tmp_path: Path) -> None:
    """Zero rows is legal — the scorer sets no minimum — and is not a reproduction."""
    report = compare_submission_files(
        _write(tmp_path, "submitted.csv", _ranking(1, 2, 3)),
        _write(tmp_path, "produced.csv", _csv([])),
    )

    assert report.produced_rows == 0
    assert report.verdict is ReproductionVerdict.SUBSET
    assert not report.top1_identical
    assert [agreement.produced_window for agreement in report.top_k] == [0, 0, 0, 0]


def test_a_mismatched_header_is_refused(tmp_path: Path) -> None:
    path = _write(tmp_path, "produced.csv", _csv([], header=TRACK1_COLUMNS[:-1]))

    with pytest.raises(SubmissionFormatError) as excinfo:
        read_ranked_keys(path, role="produced")

    assert "line 1" in str(excinfo.value)


def test_more_rows_than_the_scorer_accepts_are_refused(tmp_path: Path) -> None:
    """A file the challenge would reject is not a thing to be measured against."""
    path = _write(tmp_path, "produced.csv", _ranking(*range(MAX_SUBMISSION_ROWS + 1)))

    with pytest.raises(SubmissionFormatError) as excinfo:
        read_ranked_keys(path, role="produced")

    assert str(MAX_SUBMISSION_ROWS) in str(excinfo.value)


def test_a_non_integer_position_is_refused_by_line_number(tmp_path: Path) -> None:
    rows = [_row([_variant(1)], 0.9), _row([_variant(2)], 0.8)]
    rows[1][2] = "40,200,002"
    path = _write(tmp_path, "produced.csv", _csv(rows))

    with pytest.raises(SubmissionFormatError) as excinfo:
        read_ranked_keys(path, role="produced")

    # Line 3, because the header is line 1. The number is reportable; the line
    # is not (GP-41).
    assert "line 3" in str(excinfo.value)
    assert "40,200,002" not in str(excinfo.value)


def test_a_blank_first_variant_is_refused(tmp_path: Path) -> None:
    rows = [_row([_variant(1)], 0.9)]
    rows[0][1] = ""
    path = _write(tmp_path, "produced.csv", _csv(rows))

    with pytest.raises(SubmissionFormatError) as excinfo:
        read_ranked_keys(path, role="produced")

    assert "line 2" in str(excinfo.value)
    assert "chrom_1" in str(excinfo.value)  # a contract column name, not file content


def test_a_non_numeric_epcr_is_refused(tmp_path: Path) -> None:
    rows = [_row([_variant(1)], 0.9)]
    rows[0][9] = "high"
    path = _write(tmp_path, "produced.csv", _csv(rows))

    with pytest.raises(SubmissionFormatError) as excinfo:
        read_ranked_keys(path, role="produced")

    assert "epcr" in str(excinfo.value)
    assert "high" not in str(excinfo.value)
    assert excinfo.value.__context__ is None


def test_a_repeated_row_is_refused_naming_both_line_numbers(tmp_path: Path) -> None:
    """A duplicated row has no single rank, so rank agreement against it is undefined."""
    row = _row([_variant(1)], 0.9)
    path = _write(tmp_path, "produced.csv", _csv([row, _row([_variant(2)], 0.8), row]))

    with pytest.raises(SubmissionFormatError) as excinfo:
        read_ranked_keys(path, role="produced")

    assert "lines 2 and 4" in str(excinfo.value)


def test_undecodable_bytes_are_refused_without_the_bytes(tmp_path: Path) -> None:
    """`UnicodeDecodeError` carries the offending bytes in its own message.

    Letting it propagate, or chaining it with `raise ... from exc`, would put
    patient bytes into a traceback that reaches the terminal and the agent's
    context (PRIV-09). The offset is reportable; the byte is not.
    """
    path = tmp_path / "produced.csv"
    path.write_bytes(b",".join(c.encode() for c in TRACK1_COLUMNS) + b"\nPROBAND01,\xff\xfe\n")

    with pytest.raises(SubmissionFormatError) as excinfo:
        read_ranked_keys(path, role="produced")

    assert "not valid UTF-8" in str(excinfo.value)
    assert excinfo.value.__cause__ is None
    assert excinfo.value.__context__ is None


def test_the_error_names_the_role_and_the_basename_but_not_the_directory(
    tmp_path: Path,
) -> None:
    """GP-44: the directory is the workspace location, which the CLI never prints."""
    case_volume = tmp_path / "encrypted-case-volume"
    case_volume.mkdir()
    path = _write(case_volume, "track1_submission.csv", "")

    with pytest.raises(SubmissionFormatError) as excinfo:
        read_ranked_keys(path, role="submitted")

    message = str(excinfo.value)
    assert "submitted" in message
    assert "track1_submission.csv" in message
    assert "encrypted-case-volume" not in message


# ---------------------------------------------------------------------------
# The leak tests. This is the point of the module.
# ---------------------------------------------------------------------------

#: Anything matching this in the tool's output came out of a fixture file, and a
#: fixture file stands in for one child's variant coordinates. Two shapes: an
#: alphabetic token (`ZQX...`, which no label, verdict or remediation contains)
#: and a numeric one (`91234` followed by at least three digits, which no count,
#: line number or correlation this tool prints can reach).
_SENTINEL = re.compile(r"ZQX[A-Z0-9]+|91234\d{3,}")

#: A directory named to be recognisable, standing in for the case volume.
_SENTINEL_DIR = "ZQXCASEVOLUME"


def _leaked(*texts: str) -> list[str]:
    found: set[str] = set()
    for text in texts:
        found.update(_SENTINEL.findall(text))
    return sorted(found)


def _sentinel_row(tag: str, ordinal: int) -> list[str]:
    """One row with a distinct, recognisable token in every one of the 12 columns.

    Values still parse — the position columns are integers and the EPCR is a
    float in (0, 1] — because a fixture the comparator rejects at the first field
    would never exercise the success path, which is where a leak is most likely.
    """
    return [
        f"ZQXPROBAND{tag}",
        f"chrZQXONE{tag}{ordinal}",
        f"91234{ordinal:04d}1",
        f"ZQXREFONE{tag}{ordinal}",
        f"ZQXALTONE{tag}{ordinal}",
        f"chrZQXTWO{tag}{ordinal}",
        f"91234{ordinal:04d}2",
        f"ZQXREFTWO{tag}{ordinal}",
        f"ZQXALTTWO{tag}{ordinal}",
        f"0.91234{ordinal:04d}",
        f"ZQXFINDING{tag}{ordinal}",
        f"ZQXNOTE{tag}{ordinal}",
    ]


def _sentinel_csv(tag: str, count: int = 3) -> str:
    return _csv([_sentinel_row(tag, ordinal) for ordinal in range(1, count + 1)])


def _mutated(tag: str, column: int, value: str) -> str:
    rows = [_sentinel_row(tag, ordinal) for ordinal in range(1, 4)]
    rows[1][column] = value
    return _csv(rows)


#: `(name, submitted text, produced text, expected exit code)`. Every parsing
#: refusal is here as well as the success path: an error message is the likeliest
#: place for a record to escape, because it is written while something is wrong.
_LEAK_SCENARIOS: tuple[tuple[str, str, str, int], ...] = (
    ("identical", _sentinel_csv("A"), _sentinel_csv("A"), 0),
    ("divergent", _sentinel_csv("A"), _sentinel_csv("B"), 1),
    ("subset", _sentinel_csv("A"), _sentinel_csv("A", count=1), 1),
    ("superset", _sentinel_csv("A", count=1), _sentinel_csv("A"), 1),
    ("empty-produced", _sentinel_csv("A"), "", 1),
    ("header-only-produced", _sentinel_csv("A"), _sentinel_csv("A", count=0), 1),
    ("bad-position", _sentinel_csv("A"), _mutated("A", 2, "ZQXNOTANUMBER"), 1),
    ("bad-epcr", _sentinel_csv("A"), _mutated("A", 9, "ZQXNOTAFLOAT"), 1),
    ("blank-first-variant", _sentinel_csv("A"), _mutated("A", 1, ""), 1),
    (
        "duplicate-row",
        _sentinel_csv("A"),
        _csv([_sentinel_row("A", 1), _sentinel_row("A", 2), _sentinel_row("A", 1)]),
        1,
    ),
    (
        "too-many-rows",
        _sentinel_csv("A"),
        _sentinel_csv("A", count=MAX_SUBMISSION_ROWS + 1),
        1,
    ),
)


def _header_mismatch_fixture() -> str:
    """A header whose column NAMES are sentinels, to prove the header is content too.

    A header can be edited to hold data. Echoing "expected these columns, found
    those" would be a natural, helpful diagnostic and would print whatever the
    file's first line happens to say.
    """
    return _csv([], header=[f"ZQXCOLUMN{index}" for index in range(len(TRACK1_COLUMNS))])


@pytest.mark.parametrize(
    ("name", "submitted_text", "produced_text", "exit_code"),
    [
        pytest.param(*scenario, id=scenario[0])
        for scenario in (
            *_LEAK_SCENARIOS,
            ("header-mismatch", _sentinel_csv("A"), _header_mismatch_fixture(), 1),
        )
    ],
)
def test_the_cli_never_prints_a_byte_it_read(
    tmp_path: Path, name: str, submitted_text: str, produced_text: str, exit_code: int
) -> None:
    """Run the real command and assert no fixture token reaches stdout or stderr.

    `COLUMNS` is set wide on purpose: rich truncates a table cell that does not
    fit, and a truncated leak is a leak this assertion cannot see. A narrow
    console would make the test pass for the wrong reason.
    """
    case_volume = tmp_path / _SENTINEL_DIR
    case_volume.mkdir()
    submitted = _write(case_volume, "submitted.csv", submitted_text)
    produced = _write(case_volume, "produced.csv", produced_text)

    # Positive control: the detector is live and the fixtures really are loaded
    # with tokens. Without this the assertion below passes on empty input.
    assert _leaked(submitted_text, produced_text), f"{name}: fixture carries no sentinel"

    result = CliRunner().invoke(
        app,
        ["verify", "reproduces", str(submitted), str(produced)],
        env={"COLUMNS": "300", "TERM": "dumb"},
    )

    assert result.exit_code == exit_code, f"{name}: unexpected exit code"
    leaked = _leaked(result.stdout, result.stderr, repr(result.exception))
    assert not leaked, f"{name}: the CLI printed content read from a submission: {leaked}"
    # The directory is withheld as well: it is the case-volume path (GP-44).
    assert _SENTINEL_DIR not in result.stdout + result.stderr


@pytest.mark.parametrize(
    ("name", "submitted_text", "produced_text", "exit_code"),
    [pytest.param(*scenario, id=scenario[0]) for scenario in _LEAK_SCENARIOS],
)
def test_the_api_never_returns_or_raises_a_byte_it_read(
    tmp_path: Path, name: str, submitted_text: str, produced_text: str, exit_code: int
) -> None:
    """The same property one level down, where a traceback would carry it.

    The CLI catches `MvaError` and prints the message; an exception that escapes
    a library caller instead carries its message, its cause and its context into
    a traceback, into the log, and into a pytest failure diff (PRIV-09).
    """
    submitted = _write(tmp_path, "submitted.csv", submitted_text)
    produced = _write(tmp_path, "produced.csv", produced_text)
    assert _leaked(submitted_text, produced_text), f"{name}: fixture carries no sentinel"

    try:
        report = compare_submission_files(submitted, produced)
    except SubmissionFormatError as exc:
        assert exc.__cause__ is None and exc.__context__ is None, (
            f"{name}: a chained exception carries the original message with it"
        )
        assert not _leaked(str(exc), repr(exc)), f"{name}: the exception message leaks"
        return

    assert exit_code == 0 or not report.reproduced
    assert not _leaked(repr(report), str(report)), f"{name}: the report object leaks"


def test_the_report_carries_no_field_that_could_hold_a_row() -> None:
    """Structural, because a message can be reviewed and a data model cannot.

    A `ReproductionReport` that merely *happens* not to contain a coordinate today
    is one refactor from printing one. Every field is an int, a bool, a float, or
    a member of a closed enum, and the only nested type is `TopKAgreement`, whose
    fields are all ints.
    """
    from dataclasses import fields

    from mva.reporting.reproduction import ReproductionReport, TopKAgreement

    allowed = {"int", "bool", "float | None", "ReproductionVerdict"}
    offending = [
        f"{field.name}: {field.type}"
        for field in fields(ReproductionReport)
        if str(field.type) not in allowed | {"tuple[TopKAgreement, ...]"}
    ]
    offending += [
        f"TopKAgreement.{field.name}: {field.type}"
        for field in fields(TopKAgreement)
        if str(field.type) != "int"
    ]

    assert not offending, (
        "ReproductionReport gained a field that can hold submission content:\n  "
        + "\n  ".join(offending)
        + "\n\nRemediation: the report is read by an agent (ADR 0008). Report a "
        "count, a flag or a verdict. If a per-row handle is genuinely needed, it "
        "must be a truncated digest and never a coordinate (GP-41)."
    )
