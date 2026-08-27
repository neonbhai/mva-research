"""Exception messages must not carry patient data (PRIV-09).

An adversarial review found variant coordinates, alleles, genotype strings and
HPO terms interpolated into `ValueError` messages across the model and ingestion
layers — and pydantic additionally appending `input_value=<the whole record>`.
These are ordinary exceptions: nothing catches them, so the traceback reaches the
terminal, the log file, any crash report, and an AI agent's context.

The rule is not "be careful". It is: a message may name the field and the
constraint, and may carry a `mva.models.base.error_token` handle for
correlation, but it may never echo the value.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from mva.config import find_repo_root
from mva.models.base import error_token
from mva.models.genome import GenomeBuild, GenomicCoordinate
from mva.models.phenotype import PhenotypeObservation

pytestmark = [pytest.mark.unit, pytest.mark.privacy]

REPO = find_repo_root(Path(__file__))
SRC = REPO / "src" / "mva"

REMEDIATION = (
    "\n\nPRIV-09 remediation: an exception message may state the field name and "
    "the constraint that failed. It may carry a correlation handle from "
    "`mva.models.base.error_token(value)`. It may NOT interpolate the value. "
    "Debuggability is preserved by the token: two messages about the same record "
    "share it within a run, and nothing outside the run can reverse it."
)


def test_pydantic_does_not_echo_the_input_record() -> None:
    """`hide_input_in_errors` must stay on: the default echoes the whole record."""
    with pytest.raises(ValueError) as excinfo:
        GenomicCoordinate(
            build=GenomeBuild.GRCH38, contig="chr12", position=9_999_999, ref="C", alt="C"
        )
    message = str(excinfo.value)
    assert "input_value" not in message, (
        "pydantic is echoing the input record into the error message" + REMEDIATION
    )


@pytest.mark.parametrize(
    ("factory", "secret"),
    [
        (
            lambda: GenomicCoordinate(
                build=GenomeBuild.GRCH38,
                contig="chr12",
                position=9_999_999,
                ref="C",
                alt="C",
            ),
            "9999999",
        ),
        (
            lambda: GenomicCoordinate(
                build=GenomeBuild.GRCH38,
                contig="chr9_KI270717v1_random",
                position=1,
                ref="A",
                alt="T",
            ),
            "KI270717",
        ),
        (
            lambda: PhenotypeObservation(
                hpo_id="HP:123",
                label="x",
                status="observed",  # type: ignore[arg-type]
                provenance="test",
                extraction_confidence=1.0,
            ),
            "HP:123",
        ),
    ],
)
def test_error_messages_do_not_echo_the_offending_value(factory: object, secret: str) -> None:
    with pytest.raises(ValueError) as excinfo:
        factory()  # type: ignore[operator]
    message = str(excinfo.value)
    assert secret not in message, (
        f"the value {secret!r} appears in the exception message" + REMEDIATION
    )


def test_cross_build_comparison_does_not_echo_coordinates() -> None:
    a = GenomicCoordinate(
        build=GenomeBuild.GRCH38, contig="chr12", position=9_999_999, ref="C", alt="T"
    )
    b = GenomicCoordinate(
        build=GenomeBuild.GRCH37, contig="chr12", position=8_888_888, ref="A", alt="G"
    )
    with pytest.raises(ValueError) as excinfo:
        a.assert_same_build(b)
    message = str(excinfo.value)
    assert "9999999" not in message and "8888888" not in message, (
        "a build-mismatch error echoed both coordinates" + REMEDIATION
    )
    # The builds themselves are safe and necessary to state.
    assert "GRCh38" in message and "GRCh37" in message


def test_error_token_is_stable_within_a_run_and_not_the_value() -> None:
    """The handle must correlate without disclosing."""
    value = "GRCh38:chr15:40200000:C:T"
    token = error_token(value)
    assert token == error_token(value), "token is not stable within the run"
    assert value not in token
    assert "40200000" not in token
    assert re.fullmatch(r"[0-9a-f]{8}", token), "token should be a short hex handle"
    assert error_token(value) != error_token(value + "A"), "token does not discriminate"


def test_no_module_interpolates_a_genotype_string_into_a_message() -> None:
    """Static guard: `genotype_string` must not reach an exception message."""
    offenders: list[str] = []
    pattern = re.compile(r"msg\s*=.*genotype_string", re.DOTALL)
    for path in sorted(SRC.rglob("*.py")):
        if "__pycache__" in path.parts:
            continue
        text = path.read_text(encoding="utf-8")
        for lineno, line in enumerate(text.splitlines(), start=1):
            if pattern.search(line):
                offenders.append(f"  {path.relative_to(REPO)}:{lineno}")
    assert not offenders, (
        "a genotype string is interpolated into an error message:\n"
        + "\n".join(offenders)
        + REMEDIATION
    )
