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
    """Static guard: `genotype_string` must not reach an exception message.

    **Line-scoped, and therefore partial.** The pattern is applied per line, so a
    `msg = (` spanning several lines never matches — which is exactly how
    `VariantRecord._phase_consistency` interpolated a genotype underneath this
    test for as long as it did. Kept because it sweeps all of `src/mva` and costs
    nothing; the guard that actually catches the shape is
    `test_no_model_validator_interpolates_a_record_field_into_a_message` below,
    which parses instead of grepping but is scoped to `mva.models`. Neither is
    complete on its own, and saying so is cheaper than a reader assuming either is.
    """
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


def test_gnomad_query_failure_redacts_region() -> None:
    """A gnomAD backend failure must not put the queried coordinate in the traceback.

    Found by adversarial review and reproduced before it was fixed: cyvcf2 puts
    the region string into whatever it raises, and the adapter let that propagate
    unwrapped. The result was a bare ``RuntimeError: chr21:5031905-5031905`` — a
    proband coordinate reaching the terminal, the log file, a crash report and an
    agent's context through the *error* path, which no amount of care on the
    success path prevents.

    The whole traceback is rendered, not just ``str(exc)``: chaining with
    ``raise ... from exc`` would keep the original message and frame visible even
    though the new message is clean, so the fix has to suppress the context and
    this test has to look at what a reader would actually see.
    """
    import traceback
    from collections.abc import Iterator

    import mva.annotation.gnomad_sites as gnomad
    from mva.errors import AdapterUnavailableError

    fixture = (
        REPO / "tests" / "fixtures" / "gnomad" / "gnomad.exomes.v4.1.sites.chr21.slice.vcf.bgz"
    )
    if not fixture.is_file():  # pragma: no cover - fixture is committed
        pytest.skip(f"gnomAD fixture not present at {fixture}")

    real_open = gnomad._open_reader

    class _Exploding:
        """A reader that raises the region it was handed, as cyvcf2 does."""

        def __init__(self, inner: object) -> None:
            self._inner = inner

        @property
        def raw_header(self) -> str:
            return str(getattr(self._inner, "raw_header", ""))

        def __call__(self, region: str) -> Iterator[object]:
            raise RuntimeError(region)

        def close(self) -> None:
            close = getattr(self._inner, "close", None)
            if callable(close):
                close()

    # Built rather than written as a literal. A traceback prints each frame's
    # SOURCE LINE, so a hard-coded coordinate in this test's own call would show up
    # in the rendering and fail the assertion for a reason that has nothing to do
    # with the adapter. Real callers pass a variable, which is what this mimics.
    queried = ":".join(["GRCh38", "chr" + "21", str(5_031_905), "C", "A"])

    monkeypatch = pytest.MonkeyPatch()
    try:

        def exploding_open(path: Path) -> _Exploding:
            return _Exploding(real_open(path))

        monkeypatch.setattr(gnomad, "_open_reader", exploding_open)
        adapter = gnomad.GnomadSitesFrequencyAdapter(fixture, release="v4.1", subset="exomes")
        try:
            with pytest.raises(AdapterUnavailableError) as excinfo:
                adapter.frequencies([queried])
        finally:
            adapter.close()
    finally:
        monkeypatch.undo()

    rendered = "".join(
        traceback.format_exception(type(excinfo.value), excinfo.value, excinfo.value.__traceback__)
    )
    assert "5031905" not in rendered, "the queried position reached the traceback" + REMEDIATION
    assert "chr21" not in rendered, "the queried contig reached the traceback" + REMEDIATION
    assert "AdapterUnavailableError" in rendered
    assert re.search(r"<region:[0-9a-f]{8}>", rendered), (
        "the sanitized message must still carry a correlation handle" + REMEDIATION
    )


# ---------------------------------------------------------------------------
# Model validators (PRIV-09)
#
# `hide_input_in_errors=True` stops pydantic appending `input_value=<the whole
# record>`. It does NOT touch a string a validator built by hand, and a
# `model_validator(mode="after")` runs with the fully-constructed record in scope,
# so it is exactly where an f-string reaches for the value. Two of them did.
# ---------------------------------------------------------------------------


def _phased_without_pipe() -> object:
    """A record marked phased whose GT has no '|'. Real: a caller that sets PS
    without emitting a phased separator, which is what makes this reachable."""
    from mva.models.variant import FilterStatus, Genotype, VariantRecord, Zygosity

    return VariantRecord(
        coordinate=GenomicCoordinate(
            build=GenomeBuild.GRCH38,
            contig="chr" + "15",
            position=40_200_000,
            ref="C",
            alt="T",
        ),
        genotype=Genotype(
            zygosity=Zygosity.HET,
            genotype_string="0" + "/" + "1",
            phased=True,
        ),
        filter_status=FilterStatus.PASS,
        source_artifact="tests/unit/test_error_message_privacy.py",
    )


def test_phase_inconsistency_message_carries_neither_coordinate_nor_genotype() -> None:
    """GP-15's own validator was disclosing the genotype it refused to trust.

    The trigger is a routine callset defect, so this is not an exotic path: the
    message reached the terminal, the log file and an agent's context carrying a
    proband coordinate *and* the call at it. Under this project's threat model a
    handful of rare coordinates identifies an individual and their parents, and a
    genotype identifies their parents' too.
    """
    with pytest.raises(ValueError) as excinfo:
        _phased_without_pipe()
    message = str(excinfo.value)

    assert "40200000" not in message, "the position reached the message" + REMEDIATION
    assert "chr15" not in message, "the contig reached the message" + REMEDIATION
    assert "0/1" not in message, "the genotype reached the message" + REMEDIATION

    # Still debuggable without the data: the rule, the two fields that disagree,
    # and a correlation handle.
    assert "phased" in message and "genotype_string" in message, (
        "the message must name the fields that disagree, or it just moves the cost"
    )
    assert re.search(r"<variant:[0-9a-f]{8}>", message), (
        "the message must carry a correlation handle" + REMEDIATION
    )


def _self_paired() -> object:
    """A pairing defect: one variant offered as both halves of a compound het."""
    from mva.models.pair import (
        CandidatePair,
        ComponentScores,
        InheritanceModel,
        PhaseEvidence,
        PhaseStatus,
        make_pair_id,
    )
    from mva.models.variant import FilterStatus, Genotype, VariantRecord, Zygosity

    variant = VariantRecord(
        coordinate=GenomicCoordinate(
            build=GenomeBuild.GRCH38,
            contig="chr" + "15",
            position=40_200_000,
            ref="C",
            alt="T",
        ),
        genotype=Genotype(zygosity=Zygosity.HET, genotype_string="0" + "/" + "1"),
        filter_status=FilterStatus.PASS,
        source_artifact="tests/unit/test_error_message_privacy.py",
    )
    return CandidatePair(
        pair_id=make_pair_id("SYNTH1", (variant.variant_id, variant.variant_id)),
        gene_symbol="SYNTH1",
        variant_a=variant,
        variant_b=variant,
        inheritance_model=InheritanceModel.COMPOUND_HETEROZYGOUS,
        phase=PhaseEvidence(status=PhaseStatus.UNKNOWN, method="none"),
        scores=ComponentScores(
            analytical_validity=0.9,
            rarity=0.9,
            molecular_consequence=0.9,
            inheritance_consistency=0.5,
            phenotype_similarity=0.8,
            mechanistic_relevance=0.7,
            evidence_quality=0.6,
            contradiction_penalty=0.0,
        ),
        composite_score=0.8,
        recommended_next_test="Parental segregation testing.",
    )


def test_self_pairing_message_does_not_echo_the_variant_id() -> None:
    """A pairing bug must not pay for itself with a disclosure.

    `pair_id` is safe and stays: it is a blake2b digest of the gene symbol and the
    member IDs, it is already written into every artifact, and without it the
    message cannot be tied to the candidate it is about. The full variant ID is a
    coordinate and an allele pair, and it is what had to go.
    """
    with pytest.raises(ValueError) as excinfo:
        _self_paired()
    message = str(excinfo.value)

    assert "40200000" not in message, "the position reached the message" + REMEDIATION
    assert "GRCh38:chr15" not in message, "the coordinate reached the message" + REMEDIATION
    assert "PAIR-SYNTH1-" in message, "the pair identifier is the safe handle and must stay"
    assert "HOMOZYGOUS_RECESSIVE" in message, "the message must still say what to do instead"
    assert re.search(r"<variant:[0-9a-f]{8}>", message), (
        "the message must carry a correlation handle" + REMEDIATION
    )


def test_no_model_validator_interpolates_a_record_field_into_a_message() -> None:
    """The static half, over `mva.models`, and AST-based rather than line-based.

    The previous guard was a per-line regex for `msg = ... genotype_string`. Both
    defects it was meant to catch were written as multi-line parenthesised
    assignments, so the pattern never matched a single line and the guard passed
    over the exact code it was written for. Matching the *shape* of the expression
    instead of its formatting is the difference between a lint and a decoration.

    Values wrapped in `error_token(...)` are pruned: that is the sanctioned way to
    put a record in a message, and flagging it would push authors back to raw
    interpolation.
    """
    import ast

    forbidden = {"genotype_string", "variant_id", "position", "ref", "alt"}

    def unprotected(node: ast.AST) -> list[ast.AST]:
        found: list[ast.AST] = []
        stack: list[ast.AST] = [node]
        while stack:
            current = stack.pop()
            if isinstance(current, ast.Call):
                func = current.func
                name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", None)
                if name == "error_token":
                    continue
            found.append(current)
            stack.extend(ast.iter_child_nodes(current))
        return found

    offenders: list[str] = []
    for path in sorted((SRC / "models").rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign) and any(
                isinstance(target, ast.Name) and target.id in {"msg", "message"}
                for target in node.targets
            ):
                value: ast.AST | None = node.value
            elif isinstance(node, ast.Raise) and node.exc is not None:
                value = node.exc
            else:
                continue
            names = {
                sub.attr
                for sub in unprotected(value)
                if isinstance(sub, ast.Attribute) and sub.attr in forbidden
            }
            if names:
                offenders.append(
                    f"  {path.relative_to(REPO)}:{node.lineno} — interpolates "
                    f"{', '.join(sorted(names))}"
                )
    assert not offenders, (
        "a record field is interpolated into an exception message:\n"
        + "\n".join(offenders)
        + REMEDIATION
    )
