"""Does the composition root actually hand the reference genome to the adapters?

This file exists because of one specific, measured, invisible failure. Both
`ClinvarVcfAdapter` and `GnomadSitesFrequencyAdapter` take `reference=` as an
*optional* keyword. Constructed without it they still work, still pass every one
of their own tests, and silently run in a trim-only mode that cannot reconcile a
right-shifted indel — and that mode is the default.

A test that merely asserted "the binding constructed two adapters" would pass in
exactly the state this whole exercise exists to prevent. So the tests below assert
the reference **reached** them, two ways:

* through the adapters' own typed self-report (`representation_status`, which is
  derived from whether a reference was supplied), and
* behaviourally, with a recording reference that logs every `fetch`: the join
  path is driven with a real indel from the fixture release and the reference has
  to have been consulted.

The counterfactual is asserted too, because a degraded run must be *loud*: with
`reference=None` both adapters report `UNAVAILABLE_NO_REFERENCE` and
`representation_warnings` states the cost in the run manifest.
"""

from __future__ import annotations

import os
from collections.abc import Iterator
from pathlib import Path

import pytest

from mva.alleles import LeftAlignmentStatus
from mva.annotation.binding import (
    AdapterBindingError,
    BoundAdapters,
    ResolvedResources,
    build_clinical_adapter,
    build_frequency_adapter,
    representation_warnings,
)
from mva.annotation.clinvar_vcf import ClinvarVcfAdapter
from mva.annotation.gnomad_sites import GnomadSitesFrequencyAdapter
from mva.config import CaseConfig, Workspace, find_repo_root
from mva.determinism import hash_file
from mva.errors import ConfigError
from mva.models.genome import GenomeBuild
from mva.orchestrator import _bind_adapters, _resolve_case_resources

pytestmark = [pytest.mark.unit]

REPO_ROOT = find_repo_root(Path(__file__))
FIXTURES = REPO_ROOT / "tests" / "fixtures"

CLINVAR_SLICE = FIXTURES / "clinvar" / "clinvar_slice.vcf.gz"
GNOMAD_SLICE_DIR = FIXTURES / "gnomad"

#: A real deletion from the ClinVar slice (CFTR region). Chosen because it is an
#: INDEL: the reference is consulted only for alleles that can be shifted, so a
#: SNV would prove nothing about whether it was supplied.
CLINVAR_INDEL_ID = "GRCh38:chr7:117509031:CA:C"

#: A real deletion from the gnomAD chr21 slice, for the same reason.
GNOMAD_INDEL_ID = "GRCh38:chr21:5031914:GT:G"

#: gnomAD v4.1 exomes ships no chrM shard *at all* — not a failed download, a
#: property of the release. Any mitochondrial variant therefore reaches the
#: coverage-hole path even against a complete install.
MITOCHONDRIAL_ID = "GRCh38:chrM:8993:T:G"


class RecordingReference:
    """A `ReferenceLookup` that records every read and answers with valid bases.

    Returning a constant nucleotide is safe and deliberate: `canonicalise_allele`
    stops shifting the moment the base it reads does not match, so this reference
    is consulted, bounded, and never able to move a coordinate — which keeps the
    test about *whether the adapter asked* rather than about alignment arithmetic,
    which `tests/unit/test_normalise_representation.py` already covers.
    """

    def __init__(self) -> None:
        self.calls: list[tuple[str, int, int]] = []

    def fetch(self, contig: str, start: int, end: int) -> str:
        self.calls.append((contig, start, end))
        return "A" * (end - start + 1)


def _resolved(reference_fasta: Path | None = None) -> ResolvedResources:
    """A resolution over the committed public-reference fixtures (ADR 0012).

    `ResolvedResources` is constructed directly rather than through
    `resolve_real_resources` because the two halves are deliberately separate
    functions: resolution answers "are these releases present, registered and
    intact?" and construction opens handles. Building the value here exercises the
    second half against fixtures without a `verify=False` switch existing anywhere
    — there is no way to talk the production path out of verifying, because
    verification is a different function that it always calls.
    """
    unused = FIXTURES / "mane"
    return ResolvedResources(
        reference_fasta=reference_fasta or (unused / "not-used-by-these-two-slots.fa"),
        clinvar_vcf=CLINVAR_SLICE,
        clinvar_sha256=hash_file(CLINVAR_SLICE),
        gnomad_dir=GNOMAD_SLICE_DIR,
        gnomad_release="v4.1",
        gnomad_contigs=("chr21",),
        mane_gtf=unused / "MANE.GRCh38.v1.5.slice.ensembl_genomic.gtf.gz",
        mane_gtf_sha256="",
        mane_summary=unused / "MANE.GRCh38.v1.5.slice.summary.txt.gz",
        mane_summary_sha256="",
        snpeff_jar=unused / "snpEff.jar",
        snpeff_config=unused / "snpEff.config",
        snpeff_data_dir=unused,
        snpeff_genome="GRCh38.115",
        snpeff_java_binary=unused / "java",
        snpeff_pins=unused / "snpeff_pins.json",
        build=GenomeBuild.GRCH38,
    )


@pytest.fixture
def clinvar_with_reference() -> Iterator[tuple[ClinvarVcfAdapter, RecordingReference]]:
    reference = RecordingReference()
    adapter = build_clinical_adapter(_resolved(), reference=reference)
    yield adapter, reference
    adapter.close()


@pytest.fixture
def gnomad_with_reference() -> Iterator[tuple[GnomadSitesFrequencyAdapter, RecordingReference]]:
    reference = RecordingReference()
    wrapper = build_frequency_adapter(_resolved(), reference=reference)
    yield wrapper.inner, reference
    wrapper.close()


# ---------------------------------------------------------------------------
# The reference reaches BOTH joining adapters
# ---------------------------------------------------------------------------


def test_the_binding_gives_clinvar_the_reference_and_clinvar_uses_it(
    clinvar_with_reference: tuple[ClinvarVcfAdapter, RecordingReference],
) -> None:
    """ClinVar's half of the trap, asserted on behaviour and not only on shape.

    Measured cost of getting this wrong, on the real 2026-08-22 release over
    chr17:43,000,000-43,520,000: of 2,215 indel ALT alleles in repeat tracts, 0
    right-shifted spellings join without a reference against 2,211 with one — and
    1,761 of those are assertions ClinVar calls Pathogenic or Likely pathogenic.
    """
    adapter, reference = clinvar_with_reference

    assert adapter.representation_status is LeftAlignmentStatus.APPLIED
    assert adapter.representation_limitation is None

    adapter.assertions([CLINVAR_INDEL_ID])
    assert reference.calls, (
        "ClinVar was constructed but never consulted the reference on an indel lookup: "
        "the keyword was accepted and dropped, which is the failure this test exists for"
    )


def test_the_binding_gives_gnomad_the_reference_and_gnomad_uses_it(
    gnomad_with_reference: tuple[GnomadSitesFrequencyAdapter, RecordingReference],
) -> None:
    """gnomAD's half of the same trap.

    Measured on the real chr21 shard, a 520 kb exonic window: of 1,029 indel
    records in repeat tracts, 0 right-shifted spellings join without a reference
    against 989 with one. Thirty of the recovered records are variants gnomAD
    itself calls common, 22 of them above 5% AF — scored as novel and ultra-rare
    without it, which is the strongest promoting signal this pipeline has.
    """
    adapter, reference = gnomad_with_reference

    assert adapter.representation_status is LeftAlignmentStatus.APPLIED
    assert adapter.representation_limitation is None

    adapter.frequencies([GNOMAD_INDEL_ID])
    assert reference.calls, (
        "gnomAD was constructed but never consulted the reference on an indel lookup"
    )


def test_a_reference_less_binding_says_so_out_loud_for_both_slots() -> None:
    """The counterfactual. A degraded run must not be able to look healthy (GP-14).

    `reference=None` remains *representable* — refusing to model a degraded state
    only pushes the degradation somewhere less visible — but it is reported by
    both adapters and surfaced into the run warnings, naming what it costs.
    """
    clinical = build_clinical_adapter(_resolved(), reference=None)
    wrapper = build_frequency_adapter(_resolved(), reference=None)
    try:
        assert clinical.representation_status is LeftAlignmentStatus.UNAVAILABLE_NO_REFERENCE
        assert wrapper.inner.representation_status is LeftAlignmentStatus.UNAVAILABLE_NO_REFERENCE
        assert clinical.representation_limitation is not None
        assert wrapper.inner.representation_limitation is not None

        warnings = representation_warnings(clinical=clinical, frequency=wrapper.inner)
        assert len(warnings) == 2, "a degraded run reported fewer than both degraded slots"
        assert any("ClinVar" in warning for warning in warnings)
        assert any("gnomAD" in warning for warning in warnings)
        assert all("ADR 0018" in warning for warning in warnings)
    finally:
        clinical.close()
        wrapper.close()


def test_a_bound_real_set_with_a_reference_raises_no_representation_warning() -> None:
    """The healthy state is silent, so the degraded one is legible when it appears."""
    reference = RecordingReference()
    clinical = build_clinical_adapter(_resolved(), reference=reference)
    wrapper = build_frequency_adapter(_resolved(), reference=reference)
    try:
        assert representation_warnings(clinical=clinical, frequency=wrapper.inner) == ()
    finally:
        clinical.close()
        wrapper.close()


# ---------------------------------------------------------------------------
# The coverage hole that is a property of the release, not a bug
# ---------------------------------------------------------------------------


def test_a_contig_the_release_does_not_cover_is_reported_not_raised(
    gnomad_with_reference: tuple[GnomadSitesFrequencyAdapter, RecordingReference],
) -> None:
    """gnomAD v4.1 exomes has no chrM shard. One chrM call must not abort a run.

    `frequencies()` on the raw adapter fails closed on a contig it cannot query,
    which is right for that API and wrong for a whole-callset run: a proband VCF
    routinely holds a few dozen mitochondrial calls. The wrapper answers what is
    answerable and NAMES the hole, so the variant is recorded as having no
    frequency data — absence of a resource, never evidence of rarity (GP-14).
    """
    adapter, reference = gnomad_with_reference
    wrapper = build_frequency_adapter(_resolved(), reference=reference)
    try:
        with pytest.raises(Exception, match="no complete shard"):
            adapter.frequencies([MITOCHONDRIAL_ID])

        result = wrapper.frequencies([MITOCHONDRIAL_ID])
        assert MITOCHONDRIAL_ID not in result
        assert wrapper.unqueryable_contigs == ("chrM",)
        assert wrapper.unqueryable_count == 1

        (warning,) = wrapper.warnings()
        assert "chrM=1" in warning
        assert "GP-14" in warning
        assert "8993" not in warning, "a coordinate was echoed into a warning (PRIV-09)"
    finally:
        wrapper.close()


def test_the_coverage_warning_is_empty_until_a_hole_is_actually_hit(
    gnomad_with_reference: tuple[GnomadSitesFrequencyAdapter, RecordingReference],
) -> None:
    """A gap is reported when it happens, not declared in advance from the shard list."""
    _, reference = gnomad_with_reference
    wrapper = build_frequency_adapter(_resolved(), reference=reference)
    try:
        wrapper.frequencies([GNOMAD_INDEL_ID])
        assert wrapper.warnings() == ()
    finally:
        wrapper.close()


# ---------------------------------------------------------------------------
# The composition root's policy: a real case never falls back to fabricated data
# ---------------------------------------------------------------------------


def test_a_synthetic_case_binds_the_synthetic_tables_even_with_resources_present(
    synthetic_config: CaseConfig,
    synthetic_workspace: Workspace,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The demo must not change behaviour on a machine that happens to hold 200 GB.

    A suite whose result depends on `$MVA_RESOURCES` passes in one place and fails
    in another, and the synthetic case's `SYNTH*` genes do not exist in any real
    gene model anyway.
    """
    monkeypatch.setenv("MVA_RESOURCES", str(REPO_ROOT.parent))
    assert synthetic_config.synthetic is True
    assert _resolve_case_resources(synthetic_config, synthetic_workspace) is None


def test_a_real_case_without_a_resource_root_refuses_rather_than_faking_it(
    synthetic_config: CaseConfig,
    synthetic_workspace: Workspace,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The single most consequential branch in the composition root (ADR 0027).

    A real proband ranked against `knowledge/public/` — fictional genes, invented
    allele frequencies — would produce a submission, a dossier and a provenance
    manifest that all looked entirely healthy and were entirely fabricated. So the
    absence of the real releases is a refusal, and the message says why rather
    than only what.
    """
    monkeypatch.delenv("MVA_RESOURCES", raising=False)
    real_case = synthetic_config.model_copy(update={"synthetic": False})

    with pytest.raises(ConfigError) as excinfo:
        _resolve_case_resources(real_case, synthetic_workspace)

    message = str(excinfo.value)
    assert "MVA_RESOURCES" in message
    assert "knowledge/public/" in message, "the refusal does not name what it refused to use"
    assert "ADR 0027" in message


def test_a_real_case_with_an_empty_resource_root_refuses_at_resolution(
    synthetic_config: CaseConfig,
    synthetic_workspace: Workspace,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A root that exists but holds nothing is the same refusal, one step later.

    This is the shape of a half-finished acquisition: the directory is there, so
    root resolution succeeds, and every release inside it is absent. Failing here
    rather than at the first lookup means the run has read no patient record yet.
    """
    empty = tmp_path / "resources"
    empty.mkdir()
    monkeypatch.setenv("MVA_RESOURCES", str(empty))
    real_case = synthetic_config.model_copy(update={"synthetic": False})

    with pytest.raises(AdapterBindingError) as excinfo:
        _resolve_case_resources(real_case, synthetic_workspace)
    assert "knowledge/public/" in str(excinfo.value)


def test_the_synthetic_branch_still_binds_exactly_what_it_bound_before(
    synthetic_workspace: Workspace,
) -> None:
    """Part 3: the synthetic path is unchanged, adapters and all."""
    from contextlib import ExitStack

    knowledge_root = REPO_ROOT / "knowledge"
    with ExitStack() as stack:
        bound = _bind_adapters(
            knowledge_root=knowledge_root,
            manifest_path=knowledge_root / "manifests" / "knowledge.yaml",
            resolved=None,
            reference=None,
            stack=stack,
        )
        assert isinstance(bound, BoundAdapters)
        assert bound.adapters.consequence.name == "local-tsv-consequence"
        assert bound.adapters.frequency.name == "local-tsv-frequency"
        assert bound.coverage is None
        assert bound.run_warnings() == ()


# ---------------------------------------------------------------------------
# Against the real releases, when they are installed
# ---------------------------------------------------------------------------

_RESOURCES = Path(os.environ.get("MVA_RESOURCES") or REPO_ROOT.parent / "mva-resources")


@pytest.mark.slow
@pytest.mark.integration
def test_resolution_against_the_real_releases_binds_the_real_reference(
    synthetic_config: CaseConfig, synthetic_workspace: Workspace, monkeypatch: pytest.MonkeyPatch
) -> None:
    """End to end on the actual installed releases, or skipped with the reason.

    A representational failure that has never been checked against real data must
    not be able to hide behind a green suite that quietly checked nothing, so this
    skips loudly rather than passing vacuously.
    """
    if not (_RESOURCES / "clinvar" / "clinvar.vcf.gz").is_file():
        pytest.skip(f"no installed reference releases under {_RESOURCES}")

    monkeypatch.setenv("MVA_RESOURCES", str(_RESOURCES))
    real_case = synthetic_config.model_copy(update={"synthetic": False})
    resolved = _resolve_case_resources(real_case, synthetic_workspace)

    assert resolved is not None
    assert resolved.reference_fasta.is_file(), "the GRCh38 FASTA is required, not optional"
    assert resolved.gnomad_release.startswith("v")
    assert "chr21" in resolved.gnomad_contigs
    assert "chrM" not in resolved.gnomad_contigs, (
        "gnomAD v4.1 exomes has no chrM shard; requiring one would fail every complete install"
    )
    assert resolved.verified_count > 0
