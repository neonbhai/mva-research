"""Unit tests for the real, SnpEff-backed consequence adapter.

Most of these run against a **stub SnpEff**
(``tests/fixtures/synthetic/consequence/stub_snpeff.sh``)
rather than the real ~1 GB installation, and that is deliberate rather than a
compromise. The adapter under test is entirely real; what is faked is the
counterparty, so that the properties which actually go wrong in a subprocess
integration — an offline flag quietly dropped, a chromosome name silently
unmatched, a transcript list quietly collapsed, an input VCF echoed into a
traceback — are locked by tests that run in a second on any machine, with no JRE.

The real tool is exercised by :func:`test_real_snpeff_annotates_a_known_variant`,
which skips when the installation produced by ``tools/setup/install_snpeff.sh``
is not present.
"""

from __future__ import annotations

import json
import math
import shutil
import stat
import subprocess
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

import pytest

from mva.annotation import snpeff_local
from mva.annotation.base import ConsequenceAdapter, is_synthetic
from mva.annotation.service import DEFAULT_ANNOTATION_BATCH_SIZE
from mva.annotation.snpeff_local import (
    DEFAULT_BATCH_SIZE,
    JAVA_HOME_ENV,
    SNPEFF_ADAPTER_NAME,
    SNPEFF_OFFLINE_FLAGS,
    SnpEffArtifactPins,
    SnpEffConsequenceAdapter,
    SnpEffRunReport,
    consequence_sort_key,
    genome_is_declared,
    load_mane_select_ids,
    parse_ann_entries,
    plan_batch,
    render_input_vcf,
    resolve_java_binary,
)
from mva.errors import AdapterUnavailableError, GenomeBuildMismatchError
from mva.models.genome import ContigStyle, GenomeBuild
from mva.models.variant import ConsequenceAnnotation, ImpactSeverity

REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURES = REPO_ROOT / "tests" / "fixtures" / "synthetic" / "consequence"

#: Where tools/setup/install_snpeff.sh puts the real installation. Kept out of the
#: repository: the database is ~1 GB.
REAL_INSTALL_ROOT = REPO_ROOT.parent / "mva-resources" / "snpeff"
REAL_DATABASE = "GRCh38.115"

#: Three transcripts in the canned table, one of them a retained-intron isoform of
#: a different gene. The whole point of `annotate` is that all three survive.
MULTI_TRANSCRIPT = "GRCh38:chr15:40200000:C:T"
SINGLE_TRANSCRIPT = "GRCh38:chr7:117559590:A:G"
INTRONIC = "GRCh38:chr11:5227002:T:A"
INTERGENIC = "GRCh38:chr19:1000000:C:G"
FRAMESHIFT = "GRCh38:chr2:47403000:G:GA"
#: Present in no canned row: the stub returns the record with no ANN at all.
UNKNOWN_TO_TOOL = "GRCh38:chr12:999999:A:C"
#: Both canned as ERROR_CHROMOSOME_NOT_FOUND.
UNKNOWN_CONTIG_A = "GRCh38:chr1:100:A:T"
UNKNOWN_CONTIG_B = "GRCh38:chr3:200:A:T"

AdapterFactory = Callable[..., SnpEffConsequenceAdapter]

#: Database names the stub SnpEff treats as behaviour switches rather than genomes.
#: The last four exit ZERO while returning output the adapter must refuse.
STUB_DATABASES = (
    "SLEEPDB",
    "FAILDB",
    "TRUNCDB",
    "UNKNOWNIDDB",
    "DUPIDDB",
    "MISSINGDB",
)


# --------------------------------------------------------------------------- rig


@pytest.fixture
def adapter_factory(tmp_path: Path) -> AdapterFactory:
    """A SnpEff 'installation' whose java is the stub script.

    Everything the adapter validates at construction is present, so the test
    exercises the real constructor path rather than bypassing it.
    """
    stub = tmp_path / "bin" / "stub_snpeff.sh"
    stub.parent.mkdir(parents=True)
    shutil.copy(FIXTURES / "stub_snpeff.sh", stub)
    shutil.copy(FIXTURES / "canned_snpeff_ann.tsv", stub.parent / "canned_snpeff_ann.tsv")
    stub.chmod(stub.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    install = tmp_path / "snpEff"
    install.mkdir()
    (install / "snpEff.jar").write_text("not a real jar", encoding="utf-8")
    # Shaped like the real 5.4c config: SnpEff resolves a genome name through this
    # file, so a database whose files exist but whose name is undeclared is a
    # construction error. The stub databases are declared for the same reason.
    (install / "snpEff.config").write_text(
        "data.dir = ./data/\n"
        + f"{REAL_DATABASE}.genome: Homo_sapiens\n"
        + "".join(f"{name}.genome: stub\n" for name in STUB_DATABASES),
        encoding="utf-8",
    )

    data_dir = tmp_path / "data"
    for database in (REAL_DATABASE, *STUB_DATABASES):
        predictor = data_dir / database / "snpEffectPredictor.bin"
        predictor.parent.mkdir(parents=True)
        predictor.write_bytes(b"stub predictor")

    def build(**overrides: Any) -> SnpEffConsequenceAdapter:
        kwargs: dict[str, Any] = {
            "java_binary": stub,
            "jar_path": install / "snpEff.jar",
            "data_dir": data_dir,
            "genome_database": REAL_DATABASE,
            "contig_style": ContigStyle.ENSEMBL,
        }
        kwargs.update(overrides)
        # Pins are mandatory, so the rig measures the stub install unless a test is
        # deliberately exercising a mismatch. Measuring here is legitimate precisely
        # because these bytes were written by this fixture moments ago.
        if "pins" not in kwargs:
            database = str(kwargs["genome_database"])
            mane = kwargs.get("mane_summary")
            predictor = Path(str(kwargs["data_dir"])) / database / "snpEffectPredictor.bin"
            if not predictor.is_file():
                # The test is exercising a missing database; _require_installed runs
                # before pin verification and is what should report it.
                kwargs["pins"] = SnpEffArtifactPins(
                    jar="0" * 64, config="0" * 64, predictor="0" * 64
                )
                return SnpEffConsequenceAdapter(**kwargs)
            kwargs["pins"] = SnpEffArtifactPins.measure(
                jar_path=Path(str(kwargs["jar_path"])),
                config_path=Path(str(kwargs.get("config_path", install / "snpEff.config"))),
                predictor_path=Path(str(kwargs["data_dir"])) / database / "snpEffectPredictor.bin",
                mane_summary=Path(str(mane)) if mane is not None else None,
            )
        return SnpEffConsequenceAdapter(**kwargs)

    return build


def _snapshot(result: Mapping[str, tuple[ConsequenceAnnotation, ...]]) -> str:
    """Byte-exact rendering of a result, including mapping and tuple order."""
    return json.dumps(
        [
            [variant_id, [annotation.model_dump(mode="json") for annotation in annotations]]
            for variant_id, annotations in result.items()
        ],
        sort_keys=False,
    )


# --------------------------------------------------------------- offline guarantees


@pytest.mark.unit
def test_argv_pins_snpeff_offline(adapter_factory: AdapterFactory) -> None:
    """PRIV-05: the offline flags are in the argv, not in a comment.

    ``-noLog`` in particular is load-bearing: SnpEff posts a usage record to its
    own server on every run unless it is passed, which would be an outbound
    connection from the annotation stage.
    """
    argv = adapter_factory().build_argv()
    for flag in SNPEFF_OFFLINE_FLAGS:
        assert flag in argv, f"{flag} missing; SnpEff may reach the network or emit a timestamp"
    assert "-noLog" in SNPEFF_OFFLINE_FLAGS
    assert "-nodownload" in SNPEFF_OFFLINE_FLAGS
    assert "-download" not in argv
    assert all(isinstance(item, str) for item in argv)
    # No shell: nothing in argv may be a shell fragment, and the database name is
    # the final positional rather than something spliced into a command string.
    assert not any(any(ch in item for ch in ";|&$`") for item in argv)
    assert argv[-1] == REAL_DATABASE


@pytest.mark.unit
def test_argv_never_asks_snpeff_to_collapse_transcripts(
    adapter_factory: AdapterFactory,
) -> None:
    """SnpEff has a `-canon` switch. Passing it would be the data-loss bug, in argv."""
    argv = adapter_factory().build_argv()
    assert "-canon" not in argv
    assert "-canonList" not in argv


@pytest.mark.unit
def test_annotation_module_imports_no_network_client() -> None:
    """The structural rule, asserted locally so this module owns its own guarantee."""
    import ast

    source = (REPO_ROOT / "src" / "mva" / "annotation" / "snpeff_local.py").read_text(
        encoding="utf-8"
    )
    forbidden = {"requests", "httpx", "urllib", "aiohttp", "http", "ftplib", "smtplib", "socket"}
    imported: set[str] = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            imported.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
            imported.add(node.module.split(".")[0])
    assert not (imported & forbidden)


# --------------------------------------------------------------- adapter identity


@pytest.mark.unit
def test_satisfies_the_consequence_adapter_protocol(adapter_factory: AdapterFactory) -> None:
    assert isinstance(adapter_factory(), ConsequenceAdapter)


@pytest.mark.unit
def test_declares_itself_non_synthetic(adapter_factory: AdapterFactory) -> None:
    """GP-20: `is_synthetic` fails closed, so this must be declared, not defaulted."""
    adapter = adapter_factory()
    assert adapter.synthetic is False
    assert is_synthetic(adapter) is False
    assert adapter.name == SNPEFF_ADAPTER_NAME


@pytest.mark.unit
def test_version_names_the_database_that_produced_the_answer(
    adapter_factory: AdapterFactory,
) -> None:
    """The version is read off the tool and the database, never invented."""
    adapter = adapter_factory()
    expected = f"5.2/{REAL_DATABASE}+{adapter.pins.composite_digest}"
    assert adapter.version == expected
    assert REAL_DATABASE in adapter.version
    # The digest is what makes this a version rather than a label: the release and
    # the database name together do NOT identify the bytes that produced an answer.
    assert len(adapter.pins.composite_digest) == 12
    # And it is stamped onto every annotation the adapter emits.
    annotations = adapter.annotate([SINGLE_TRANSCRIPT])[SINGLE_TRANSCRIPT]
    assert {a.source_tool_version for a in annotations} == {expected}
    assert {a.source_tool for a in annotations} == {SNPEFF_ADAPTER_NAME}


@pytest.mark.unit
def test_missing_database_fails_at_construction(adapter_factory: AdapterFactory) -> None:
    """An absent database is a wiring error, not a run-time surprise."""
    with pytest.raises(AdapterUnavailableError, match=r"GRCh38\.999"):
        adapter_factory(genome_database="GRCh38.999")


# ------------------------------------------------------------------- annotate()


@pytest.mark.unit
def test_every_transcript_survives(adapter_factory: AdapterFactory) -> None:
    """Rule 3 of the adapter contract: collapsing to MANE-Select is a data-loss bug."""
    annotations = adapter_factory().annotate([MULTI_TRANSCRIPT])[MULTI_TRANSCRIPT]
    assert len(annotations) == 3
    assert [a.transcript_id for a in annotations] == [
        "ENST00000287598",
        "ENST00000412359",
        "ENST00000672080",
    ]
    # Including the isoform that disagrees with the others about both gene and
    # severity -- exactly the row a "just take the canonical one" shortcut loses.
    assert {a.gene_symbol for a in annotations} == {"BUB1B", "BUB1B-PAK6"}
    assert {a.impact for a in annotations} == {ImpactSeverity.HIGH, ImpactSeverity.MODIFIER}
    assert {a.transcript_biotype for a in annotations} == {"protein_coding", "retained_intron"}


@pytest.mark.unit
def test_unknown_variant_is_omitted_not_empty(adapter_factory: AdapterFactory) -> None:
    """GP-14: absence is a missing key. An empty tuple would assert 'no consequence'.

    Read through `annotate_with_report`, which is the explicit opt-in to a batch
    the tool did not fully answer: `annotate` refuses a record with no ANN field
    rather than returning it as an absent key. See
    `test_a_record_with_no_ann_field_fails_closed_through_the_protocol`.
    """
    result, _ = adapter_factory().annotate_with_report([SINGLE_TRANSCRIPT, UNKNOWN_TO_TOOL])
    assert SINGLE_TRANSCRIPT in result
    assert UNKNOWN_TO_TOOL not in result
    assert result.get(UNKNOWN_TO_TOOL) is None


@pytest.mark.unit
def test_typed_fields_are_parsed_at_the_boundary(adapter_factory: AdapterFactory) -> None:
    """GP-02: models out, not dicts; and every optional field is genuinely optional."""
    adapter = adapter_factory()
    (missense,) = adapter.annotate([SINGLE_TRANSCRIPT])[SINGLE_TRANSCRIPT]
    assert missense.gene_symbol == "CFTR"
    assert missense.gene_id == "ENSG00000001626"
    assert missense.consequence_terms == ("missense_variant", "splice_region_variant")
    assert missense.most_severe_term == "missense_variant"
    assert missense.impact is ImpactSeverity.MODERATE
    assert missense.hgvs_c == "c.1521A>G"
    assert missense.hgvs_p == "p.Ile507Met"
    assert missense.protein_position == 507
    assert missense.exon == "11/27"
    assert missense.intron is None
    # Base SnpEff computes neither, so neither is invented (GP-14).
    assert missense.splice_ai_delta_max is None
    assert missense.pathogenicity_scores == {}
    assert missense.amino_acids is None

    (intronic,) = adapter.annotate([INTRONIC])[INTRONIC]
    assert intronic.exon is None
    assert intronic.intron == "1/2"
    assert intronic.protein_position is None
    assert intronic.hgvs_p is None

    (frameshift,) = adapter.annotate([FRAMESHIFT])[FRAMESHIFT]
    assert frameshift.impact is ImpactSeverity.HIGH
    assert frameshift.consequence_terms == ("frameshift_variant",)


@pytest.mark.unit
def test_non_transcript_features_are_kept_rather_than_dropped(
    adapter_factory: AdapterFactory,
) -> None:
    """An intergenic call is an answer. Dropping it would fabricate an absence."""
    (intergenic,) = adapter_factory().annotate([INTERGENIC])[INTERGENIC]
    assert intergenic.transcript_biotype == "intergenic_region"
    assert intergenic.impact is ImpactSeverity.MODIFIER


@pytest.mark.unit
def test_result_order_follows_the_caller(adapter_factory: AdapterFactory) -> None:
    """The subprocess input is coordinate-sorted; the result is not reordered by it."""
    requested = [INTRONIC, MULTI_TRANSCRIPT, SINGLE_TRANSCRIPT]
    assert list(adapter_factory().annotate(requested)) == requested


@pytest.mark.unit
def test_duplicate_ids_are_requested_once(adapter_factory: AdapterFactory) -> None:
    result = adapter_factory().annotate([SINGLE_TRANSCRIPT, SINGLE_TRANSCRIPT])
    assert list(result) == [SINGLE_TRANSCRIPT]


# ------------------------------------------------------------------ determinism


@pytest.mark.unit
def test_repeat_runs_are_byte_identical(adapter_factory: AdapterFactory) -> None:
    """GP-30. A subprocess is a determinism hazard; this is the acceptance criterion."""
    variants = [MULTI_TRANSCRIPT, SINGLE_TRANSCRIPT, INTRONIC, UNKNOWN_TO_TOOL, FRAMESHIFT]
    adapter = adapter_factory()
    # `annotate_with_report` because the batch deliberately contains a record the
    # stub returns with no ANN field, which `annotate` refuses; determinism has to
    # hold on the partial answer too, since that is what a caller who opted in gets.
    first = _snapshot(adapter.annotate_with_report(variants)[0])
    second = _snapshot(adapter.annotate_with_report(variants)[0])
    assert first == second
    # A freshly constructed adapter, i.e. a fresh JVM and a fresh scratch
    # directory, must agree too: no temp path or start time may leak into output.
    assert _snapshot(adapter_factory().annotate_with_report(variants)[0]) == first
    assert "/tmp" not in first  # noqa: S108 - asserting a temp path is ABSENT
    assert "mva-snpeff-" not in first
    assert "2024-09-24" not in first


@pytest.mark.unit
def test_batch_size_does_not_change_the_answer(adapter_factory: AdapterFactory) -> None:
    """Chunking is a memory bound, not a semantic one (GP-30)."""
    variants = [MULTI_TRANSCRIPT, SINGLE_TRANSCRIPT, INTRONIC, FRAMESHIFT, INTERGENIC]
    whole = _snapshot(adapter_factory().annotate(variants))
    assert _snapshot(adapter_factory(batch_size=2).annotate(variants)) == whole
    assert _snapshot(adapter_factory(batch_size=1).annotate(variants)) == whole


@pytest.mark.unit
def test_batch_size_must_be_positive(adapter_factory: AdapterFactory) -> None:
    with pytest.raises(AdapterUnavailableError, match="batch_size"):
        adapter_factory(batch_size=0)


@pytest.mark.unit
def test_transcript_order_is_a_total_documented_key(adapter_factory: AdapterFactory) -> None:
    """Ordering is presentation. It must still be total, or repeat runs diverge."""
    annotations = adapter_factory(mane_summary=FIXTURES / "mane_summary_excerpt.tsv").annotate(
        [MULTI_TRANSCRIPT]
    )[MULTI_TRANSCRIPT]
    keys = [consequence_sort_key(a) for a in annotations]
    assert keys == sorted(keys)
    assert len(set(keys)) == len(keys), "sort key is not total; ties are nondeterministic"
    # MANE-Select sorts first -- and is still only first, never alone.
    assert annotations[0].is_mane_select is True
    assert len(annotations) == 3


# ------------------------------------------------------------------ contig names


@pytest.mark.unit
def test_contig_style_is_mapped_explicitly() -> None:
    """The repo key is `chr`-prefixed; SnpEff's Ensembl databases are not."""
    (ensembl_site,) = plan_batch(
        [MULTI_TRANSCRIPT], build=GenomeBuild.GRCH38, contig_style=ContigStyle.ENSEMBL
    )
    assert ensembl_site.contig == "15"
    assert ensembl_site.variant_id == MULTI_TRANSCRIPT
    (ucsc_site,) = plan_batch(
        [MULTI_TRANSCRIPT], build=GenomeBuild.GRCH38, contig_style=ContigStyle.UCSC
    )
    assert ucsc_site.contig == "chr15"
    assert "\n15\t40200000\t" in render_input_vcf([ensembl_site])
    assert "\nchr15\t40200000\t" in render_input_vcf([ucsc_site])


@pytest.mark.unit
def test_mitochondrial_contig_maps_to_the_ensembl_spelling() -> None:
    (site,) = plan_batch(
        ["GRCh38:chrM:8993:T:G"], build=GenomeBuild.GRCH38, contig_style=ContigStyle.ENSEMBL
    )
    assert site.contig == "MT"


@pytest.mark.unit
def test_a_whole_run_of_unknown_contigs_is_raised_not_returned_empty(
    adapter_factory: AdapterFactory,
) -> None:
    """Every variant unplaceable is a configuration failure, and must not look clean.

    A run where the database recognises no contig returns zero annotations for
    everything, which is indistinguishable from "this exome has nothing to
    report". Measured against SnpEff 5.4c, the tool strips a leading `chr` itself,
    so this fires on a database that genuinely lacks the contigs (a RefSeq
    accession naming scheme, a missing MT) rather than on `chr15` alone.
    """
    adapter = adapter_factory()
    with pytest.raises(AdapterUnavailableError) as excinfo:
        adapter.annotate([UNKNOWN_CONTIG_A, UNKNOWN_CONTIG_B])
    message = str(excinfo.value)
    assert "contig_style" in message
    assert "ensembl" in message
    # PRIV-09: the diagnosis names the setting, never a coordinate.
    assert "100" not in message.replace("PRIV-09", "")
    assert "chr1" not in message


@pytest.mark.unit
def test_sites_the_tool_cannot_place_do_not_poison_the_whole_batch(
    adapter_factory: AdapterFactory,
) -> None:
    """One unknown contig among real hits is dropped, not escalated."""
    result = adapter_factory().annotate([UNKNOWN_CONTIG_A, SINGLE_TRANSCRIPT])
    assert list(result) == [SINGLE_TRANSCRIPT]


@pytest.mark.unit
def test_cross_build_annotation_is_refused() -> None:
    """GRCh37 coordinates against a GRCh38 database return confident nonsense."""
    with pytest.raises(GenomeBuildMismatchError):
        plan_batch(
            ["GRCh37:chr15:40200000:C:T"],
            build=GenomeBuild.GRCH38,
            contig_style=ContigStyle.ENSEMBL,
        )


@pytest.mark.unit
def test_unannotatable_alleles_are_skipped_not_sent() -> None:
    """Spanning deletions and missing alleles are legal VCF but not annotatable."""
    assert (
        plan_batch(
            ["GRCh38:chr15:40200000:C:*"],
            build=GenomeBuild.GRCH38,
            contig_style=ContigStyle.ENSEMBL,
        )
        == ()
    )


@pytest.mark.unit
def test_malformed_variant_id_is_tokenised_not_echoed() -> None:
    with pytest.raises(AdapterUnavailableError) as excinfo:
        plan_batch(
            ["chr15-40200000-C-T"], build=GenomeBuild.GRCH38, contig_style=ContigStyle.ENSEMBL
        )
    assert "40200000" not in str(excinfo.value)


@pytest.mark.unit
def test_input_vcf_carries_no_timestamp_and_no_genotype() -> None:
    sites = plan_batch(
        [MULTI_TRANSCRIPT, SINGLE_TRANSCRIPT],
        build=GenomeBuild.GRCH38,
        contig_style=ContigStyle.ENSEMBL,
    )
    vcf = render_input_vcf(sites)
    assert "##fileDate" not in vcf
    assert "FORMAT" not in vcf
    assert vcf.startswith("##fileformat=VCFv4.2\n")
    # Coordinate-sorted, so the bytes depend on the set of variants, not the order
    # the caller iterated in.
    assert (
        render_input_vcf(
            plan_batch(
                [SINGLE_TRANSCRIPT, MULTI_TRANSCRIPT],
                build=GenomeBuild.GRCH38,
                contig_style=ContigStyle.ENSEMBL,
            )
        )
        == vcf
    )


# ------------------------------------------------------------- subprocess hygiene


@pytest.mark.unit
def test_failure_names_the_tool_and_leaks_nothing(adapter_factory: AdapterFactory) -> None:
    """PRIV-09. SnpEff echoes offending VCF records on failure; we must not relay them."""
    adapter = adapter_factory(genome_database="FAILDB")
    with pytest.raises(AdapterUnavailableError) as excinfo:
        adapter.annotate([MULTI_TRANSCRIPT, SINGLE_TRANSCRIPT])
    message = str(excinfo.value)
    assert "SnpEff" in message
    assert "exit code 3" in message
    for secret in ("40200000", "117559590", "mva0", "\tC\tT", "chr15", "ENST"):
        assert secret not in message, f"{secret!r} leaked into the exception message"


@pytest.mark.unit
def test_timeout_is_bounded_and_raises(adapter_factory: AdapterFactory) -> None:
    adapter = adapter_factory(genome_database="SLEEPDB", timeout_seconds=1.0)
    with pytest.raises(AdapterUnavailableError, match="timeout"):
        adapter.annotate([SINGLE_TRANSCRIPT])


@pytest.mark.unit
def test_missing_java_runtime_is_a_named_failure(
    adapter_factory: AdapterFactory, tmp_path: Path
) -> None:
    with pytest.raises(AdapterUnavailableError, match="Java runtime"):
        adapter_factory(java_binary=tmp_path / "no" / "such" / "java")


# --------------------------------------------------------------- locating the JVM
#
# The machine this was written on has a real JDK *and* macOS's /usr/bin/java stub,
# which is exactly the configuration that makes a PATH lookup dangerous: the stub
# is found first, passes every existence check, and fails only when executed.


@pytest.mark.unit
def test_java_home_is_the_only_fallback_and_path_is_never_consulted(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """PATH is not a JVM resolution mechanism here, even when it holds a working java.

    A `java` that works perfectly is placed on PATH and JAVA_HOME is cleared. The
    resolver must still refuse: on macOS the name `java` resolves to a stub, and a
    resolver that trusted PATH would silently pick whichever came first.
    """
    on_path = tmp_path / "pathbin"
    on_path.mkdir()
    shutil.copy(FIXTURES / "stub_snpeff.sh", on_path / "java")
    (on_path / "java").chmod(0o755)
    monkeypatch.setenv("PATH", str(on_path))
    monkeypatch.delenv(JAVA_HOME_ENV, raising=False)

    with pytest.raises(AdapterUnavailableError) as excinfo:
        resolve_java_binary()
    message = str(excinfo.value)
    assert JAVA_HOME_ENV in message
    assert "PATH" in message
    assert resolve_java_binary(on_path / "java") == (on_path / "java").resolve()


@pytest.mark.unit
def test_java_home_supplies_the_runtime_when_none_is_passed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The documented fallback: $JAVA_HOME/bin/java, by absolute path."""
    home = tmp_path / "jdk"
    (home / "bin").mkdir(parents=True)
    shutil.copy(FIXTURES / "stub_snpeff.sh", home / "bin" / "java")
    (home / "bin" / "java").chmod(0o755)
    monkeypatch.setenv(JAVA_HOME_ENV, str(home))
    assert resolve_java_binary() == (home / "bin" / "java").resolve()

    monkeypatch.setenv(JAVA_HOME_ENV, str(tmp_path / "not-a-jdk"))
    with pytest.raises(AdapterUnavailableError, match=JAVA_HOME_ENV):
        resolve_java_binary()


@pytest.mark.unit
def test_adapter_uses_java_home_without_an_explicit_binary(
    adapter_factory: AdapterFactory, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Constructing with no java_binary at all must still work off JAVA_HOME."""
    home = tmp_path / "jdk-home"
    (home / "bin").mkdir(parents=True)
    shutil.copy(FIXTURES / "stub_snpeff.sh", home / "bin" / "java")
    shutil.copy(FIXTURES / "canned_snpeff_ann.tsv", home / "bin" / "canned_snpeff_ann.tsv")
    (home / "bin" / "java").chmod(0o755)
    monkeypatch.setenv(JAVA_HOME_ENV, str(home))

    adapter = adapter_factory(java_binary=None)
    assert adapter.java_binary == (home / "bin" / "java").resolve()
    assert adapter.annotate([SINGLE_TRANSCRIPT])[SINGLE_TRANSCRIPT]


@pytest.mark.unit
def test_a_java_that_is_not_a_runtime_is_diagnosed_as_such(
    adapter_factory: AdapterFactory, tmp_path: Path
) -> None:
    """The macOS stub case: present, executable, not a JRE.

    This is the failure the previous shape of this adapter turned into "SnpEff
    version probe failed with exit code 1", which sends the operator to look at
    SnpEff. The message must name the JVM, and it must say how to fix it.
    """
    fake = tmp_path / "stubjava" / "java"
    fake.parent.mkdir()
    shutil.copy(FIXTURES / "stub_java_not_a_runtime.sh", fake)
    fake.chmod(0o755)

    with pytest.raises(AdapterUnavailableError) as excinfo:
        adapter_factory(java_binary=fake)
    message = str(excinfo.value)
    assert "not a working Java runtime" in message
    assert "Unable to locate a Java Runtime" in message, "the actual diagnosis is not relayed"
    assert JAVA_HOME_ENV in message
    # It must NOT be misattributed to SnpEff, which is what makes it a bug report
    # against the wrong tool.
    assert "version probe failed" not in message


# ------------------------------------------------------- the genome must be declared


@pytest.mark.unit
def test_genome_absent_from_the_config_is_refused_at_construction(
    adapter_factory: AdapterFactory, tmp_path: Path
) -> None:
    """Files present + name undeclared is a naming bug SnpEff reports as a download bug.

    With `-nodownload` set, SnpEff answers an unknown genome name by refusing to
    fetch it, which reads as a network failure in an adapter whose entire premise
    is that it never uses the network.
    """
    config = tmp_path / "bare.config"
    config.write_text("data.dir = ./data/\n", encoding="utf-8")
    with pytest.raises(AdapterUnavailableError) as excinfo:
        adapter_factory(config_path=config)
    message = str(excinfo.value)
    assert f"{REAL_DATABASE}.genome" in message
    assert "-nodownload" in message


@pytest.mark.unit
def test_genome_is_declared_reads_the_config_verbatim(tmp_path: Path) -> None:
    config = tmp_path / "snpEff.config"
    config.write_text(
        "# comment\nGRCh38.115.genome: Homo_sapiens\nGRCh37.75.genome : Human\n",
        encoding="utf-8",
    )
    assert genome_is_declared(config, "GRCh38.115") is True
    assert genome_is_declared(config, "GRCh37.75") is True
    assert genome_is_declared(config, "GRCh38.99") is False
    assert genome_is_declared(tmp_path / "absent.config", "GRCh38.115") is False


# ------------------------------------------------------- failing closed on output
#
# The dangerous case is not a crash. It is SnpEff exiting ZERO having written an
# answer that is missing records, because a missing annotation is read downstream
# as "this variant has no consequence" -- and since VariantRecord.gene_symbols is
# derived from consequences alone, that removes the variant from gene grouping and
# therefore from pairing. A truncated run would delete candidate pairs silently.


@pytest.mark.unit
@pytest.mark.parametrize(
    ("database", "expected"),
    [
        ("TRUNCDB", "truncated mid-record"),
        ("UNKNOWNIDDB", "never sent"),
        ("DUPIDDB", "came back twice"),
        ("MISSINGDB", "no record at all"),
    ],
)
def test_missing_or_malformed_output_record_fails_closed(
    adapter_factory: AdapterFactory, database: str, expected: str
) -> None:
    """Every way SnpEff can succeed while returning an unusable answer must raise."""
    adapter = adapter_factory(genome_database=database)
    with pytest.raises(AdapterUnavailableError) as excinfo:
        adapter.annotate([MULTI_TRANSCRIPT, SINGLE_TRANSCRIPT])
    message = str(excinfo.value)
    assert expected in message
    # It must explain the consequence, because "returned fewer rows" does not
    # obviously mean "candidate pairs were deleted".
    assert "pairing" in message
    # PRIV-09: the diagnosis describes the shape of the failure, never the record.
    for secret in ("40200000", "117559590", "chr15", "\tC\tT"):
        assert secret not in message, f"{secret!r} leaked into the exception message"


@pytest.mark.unit
def test_a_silently_partial_run_is_never_returned_as_an_empty_result(
    adapter_factory: AdapterFactory,
) -> None:
    """The regression this fails-closed behaviour exists for, stated directly.

    Before this, a SnpEff run that returned nothing produced `{}` -- byte-identical
    to a clean run over variants with no consequence. The two must not be the same
    value.
    """
    adapter = adapter_factory(genome_database="MISSINGDB")
    with pytest.raises(AdapterUnavailableError):
        adapter.annotate([SINGLE_TRANSCRIPT])


@pytest.mark.unit
def test_every_variant_lands_in_exactly_one_outcome_bucket(
    adapter_factory: AdapterFactory,
) -> None:
    """Omission is classified, never inferred from an absent key (GP-14)."""
    requested = [
        SINGLE_TRANSCRIPT,  # annotated
        UNKNOWN_TO_TOOL,  # a record with no ANN at all
        UNKNOWN_CONTIG_A,  # SnpEff cannot place it
        "GRCh38:chr15:40200000:C:*",  # never sent: spanning deletion
    ]
    result, report = adapter_factory().annotate_with_report(requested)

    assert report.requested == tuple(requested)
    assert report.annotated == (SINGLE_TRANSCRIPT,)
    assert report.without_ann == (UNKNOWN_TO_TOOL,)
    assert report.unplaceable == (UNKNOWN_CONTIG_A,)
    assert report.skipped_unannotatable == ("GRCh38:chr15:40200000:C:*",)
    assert set(report.unannotated) == set(requested) - {SINGLE_TRANSCRIPT}
    # The partition is checked by the adapter itself, not only by this test.
    report.assert_partitions()
    # And the mapping still omits rather than empties (GP-14).
    assert set(result) == {SINGLE_TRANSCRIPT}
    for omitted in report.unannotated:
        assert omitted not in result


@pytest.mark.unit
def test_run_report_is_exposed_after_a_plain_annotate_call(
    adapter_factory: AdapterFactory,
) -> None:
    adapter = adapter_factory()
    assert adapter.last_run_report is None
    adapter.annotate([SINGLE_TRANSCRIPT, MULTI_TRANSCRIPT])
    report = adapter.last_run_report
    assert isinstance(report, SnpEffRunReport)
    assert report.records_returned == 2
    assert report.without_ann == ()

    # And it survives the fail-closed path, so the refusal is diagnosable rather
    # than leaving the caller with an exception and no accounting.
    with pytest.raises(AdapterUnavailableError):
        adapter.annotate([SINGLE_TRANSCRIPT, UNKNOWN_TO_TOOL])
    after = adapter.last_run_report
    assert after is not None
    assert after.without_ann == (UNKNOWN_TO_TOOL,)


@pytest.mark.unit
def test_partition_check_catches_an_inconsistent_report() -> None:
    """The accounting guard is itself tested, or it is just a comment."""
    broken = SnpEffRunReport(
        requested=("a", "b"),
        skipped_unannotatable=(),
        annotated=("a",),
        without_ann=(),
        unplaceable=(),
        incomplete=(),
        records_returned=1,
    )
    with pytest.raises(AdapterUnavailableError, match="accounting is inconsistent"):
        broken.assert_partitions()


# --------------------------------------------------------------- artifact pinning


@pytest.mark.unit
def test_unpinned_installation_is_refused(adapter_factory: AdapterFactory) -> None:
    """`snpEff_latest_core.zip` is a moving target by name; an unpinned run is not
    reproducible and its provenance string is actively misleading."""
    with pytest.raises(AdapterUnavailableError) as excinfo:
        adapter_factory(pins=None)
    message = str(excinfo.value)
    assert "SnpEffArtifactPins" in message
    assert "snpEff_latest_core.zip" in message


@pytest.mark.unit
@pytest.mark.parametrize("artifact", ["snpEff.jar", "snpEff.config", "snpEffectPredictor.bin"])
def test_snpeff_rejects_unpinned_or_mismatched_artifacts(
    adapter_factory: AdapterFactory, tmp_path: Path, artifact: str
) -> None:
    """Mutate each pinned artifact independently; construction must fail each time."""
    good = adapter_factory()
    pins = good.pins

    target = {
        "snpEff.jar": good.jar_path,
        "snpEff.config": good.config_path,
        "snpEffectPredictor.bin": good.predictor_path,
    }[artifact]
    target.write_bytes(target.read_bytes() + b"\n# tampered")

    with pytest.raises(AdapterUnavailableError) as excinfo:
        adapter_factory(pins=pins)
    message = str(excinfo.value)
    assert artifact in message
    assert "failed its sha256 pin" in message
    # The remediation has to say why a changed gene model matters, not just that a
    # hash differs.
    assert "different transcripts" in message


@pytest.mark.unit
def test_mane_summary_must_be_pinned_when_it_is_used(
    adapter_factory: AdapterFactory,
) -> None:
    """A file that decides which transcripts are MANE Select is an input to the answer."""
    unpinned = SnpEffArtifactPins(
        jar=adapter_factory().pins.jar,
        config=adapter_factory().pins.config,
        predictor=adapter_factory().pins.predictor,
        mane_summary=None,
    )
    with pytest.raises(AdapterUnavailableError, match="pinned or not used"):
        adapter_factory(mane_summary=FIXTURES / "mane_summary_excerpt.tsv", pins=unpinned)


@pytest.mark.unit
def test_the_composite_digest_changes_when_any_artifact_changes() -> None:
    """Two different installations must not be able to cite the same provenance."""
    base = SnpEffArtifactPins(jar="a" * 64, config="b" * 64, predictor="c" * 64)
    assert (
        base.composite_digest
        == SnpEffArtifactPins(jar="a" * 64, config="b" * 64, predictor="c" * 64).composite_digest
    )
    for changed in (
        SnpEffArtifactPins(jar="z" * 64, config="b" * 64, predictor="c" * 64),
        SnpEffArtifactPins(jar="a" * 64, config="z" * 64, predictor="c" * 64),
        SnpEffArtifactPins(jar="a" * 64, config="b" * 64, predictor="z" * 64),
        SnpEffArtifactPins(
            jar="a" * 64, config="b" * 64, predictor="c" * 64, mane_summary="d" * 64
        ),
    ):
        assert changed.composite_digest != base.composite_digest


# ------------------------------------------------------------------------- MANE


@pytest.mark.unit
def test_mane_select_flag_comes_from_the_mane_release(adapter_factory: AdapterFactory) -> None:
    """SnpEff's ANN carries no MANE flag; it is joined from NCBI's summary or left unset."""
    without = adapter_factory()
    assert without.mane_select_count == 0
    assert all(not a.is_mane_select for a in without.annotate([MULTI_TRANSCRIPT])[MULTI_TRANSCRIPT])

    with_mane = adapter_factory(mane_summary=FIXTURES / "mane_summary_excerpt.tsv")
    assert with_mane.mane_select_count > 0
    annotations = with_mane.annotate([MULTI_TRANSCRIPT])[MULTI_TRANSCRIPT]
    flagged = {a.transcript_id for a in annotations if a.is_mane_select}
    assert flagged == {"ENST00000287598"}


@pytest.mark.unit
def test_mane_plus_clinical_is_not_mane_select() -> None:
    """The two statuses are different claims and the join must not conflate them."""
    selected = load_mane_select_ids(FIXTURES / "mane_summary_excerpt.tsv")
    assert "ENST00000287598" in selected
    assert "NM_001211" in selected
    assert "ENST00000233146" not in selected  # MSH2 row is MANE Plus Clinical


@pytest.mark.unit
def test_mane_summary_with_unexpected_columns_is_refused(tmp_path: Path) -> None:
    path = tmp_path / "bad.tsv"
    path.write_text("#gene\ttranscript\nA\tB\n", encoding="utf-8")
    with pytest.raises(AdapterUnavailableError, match="expected columns"):
        load_mane_select_ids(path)


# ------------------------------------------------------------------ ANN parsing


@pytest.mark.unit
def test_ann_percent_escapes_are_decoded() -> None:
    ann = (
        "T|missense_variant|MODERATE|GENE%3DX|ENSG1|transcript|ENST1|protein_coding|"
        "2/5|c.10A>T|p.Lys4Ter|100/900|10/800|4/266||"
    )
    parsed = parse_ann_entries(
        f"DP=10;ANN={ann};AF=0.5", tool="snpeff", tool_version="5.2/GRCh38.110"
    )
    assert parsed.contig_unknown is False
    assert parsed.present is True
    assert parsed.annotations[0].gene_symbol == "GENE=X"


@pytest.mark.unit
def test_record_without_ann_yields_nothing() -> None:
    """No ANN field at all is reported as `present=False`, not as an empty answer.

    The distinction is the whole of GP-14 at this level: an empty tuple with
    `present=True` would mean "SnpEff looked and found nothing", which is a claim
    this adapter is not entitled to make on the strength of a missing field.
    """
    parsed = parse_ann_entries(".", tool="snpeff", tool_version="v")
    assert parsed.annotations == ()
    assert parsed.present is False
    assert parsed.contig_unknown is False


@pytest.mark.unit
def test_truncated_ann_entry_is_refused_not_guessed() -> None:
    with pytest.raises(AdapterUnavailableError, match="sub-fields"):
        parse_ann_entries("ANN=T|missense_variant|MODERATE", tool="snpeff", tool_version="v")


@pytest.mark.unit
def test_reference_mismatch_warning_does_not_discard_the_annotation() -> None:
    """A REF-mismatch entry still describes a real transcript. Dropping it is data loss."""
    ann = (
        "T|missense_variant|MODERATE|GENE|ENSG1|transcript|ENST1|protein_coding|"
        "2/5|c.10A>T|p.Lys4Ter|100/900|10/800|4/266||WARNING_REF_DOES_NOT_MATCH_GENOME"
    )
    parsed = parse_ann_entries(f"ANN={ann}", tool="snpeff", tool_version="v")
    assert len(parsed.annotations) == 1
    assert parsed.contig_unknown is False


@pytest.mark.unit
def test_unplaceable_entries_are_not_returned_as_annotations() -> None:
    """An out-of-range placeholder has no transcript behind it; it is a failure notice."""
    ann = "T||MODIFIER|||||||||||||ERROR_OUT_OF_CHROMOSOME_RANGE"
    parsed = parse_ann_entries(f"ANN={ann}", tool="snpeff", tool_version="v")
    assert parsed.annotations == ()
    assert parsed.unplaceable == 1
    # Out-of-range is not a contig-naming problem, so it must not trip that alarm.
    assert parsed.contig_unknown is False


@pytest.mark.unit
def test_impact_is_always_asserted_never_left_not_assessed(
    adapter_factory: AdapterFactory,
) -> None:
    """`ConsequenceAnnotation.impact` is nullable; this adapter never uses that.

    `None` means NOT ASSESSED -- the source found the gene but computed no
    consequence -- and it is emphatically not MODIFIER. SnpEff always states an
    impact class, and an entry whose class this adapter cannot recognise is
    refused rather than degraded to None, so every annotation returned here
    carries a positive severity claim.
    """
    adapter = adapter_factory()
    for variant in (MULTI_TRANSCRIPT, SINGLE_TRANSCRIPT, INTRONIC, INTERGENIC, FRAMESHIFT):
        for annotation in adapter.annotate([variant])[variant]:
            assert annotation.impact is not None, variant


@pytest.mark.unit
def test_sort_key_ranks_not_assessed_apart_from_modifier() -> None:
    """Ordering must not collapse 'nobody predicted' into 'predicted harmless'."""

    def make(impact: ImpactSeverity | None) -> ConsequenceAnnotation:
        return ConsequenceAnnotation(
            gene_symbol="G",
            transcript_id="ENST1",
            consequence_terms=("x",),
            impact=impact,
        )

    not_assessed = consequence_sort_key(make(None))
    modifier = consequence_sort_key(make(ImpactSeverity.MODIFIER))
    high = consequence_sort_key(make(ImpactSeverity.HIGH))
    assert not_assessed != modifier, "NOT ASSESSED was conflated with MODIFIER"
    assert high < modifier < not_assessed


@pytest.mark.unit
def test_unrecognised_impact_class_is_refused() -> None:
    ann = (
        "T|missense_variant|CATASTROPHIC|GENE|ENSG1|transcript|ENST1|protein_coding|"
        "2/5|c.10A>T|p.Lys4Ter|100/900|10/800|4/266||"
    )
    with pytest.raises(AdapterUnavailableError, match="impact class"):
        parse_ann_entries(f"ANN={ann}", tool="snpeff", tool_version="v")


# -------------------------------------------------------------- the real tool
#
# Everything below runs the genuine SnpEff 5.4c against the genuine GRCh38.115
# gene model, and skips when that ~790 MB installation is absent.
#
# The two variants are REAL, PUBLIC, PATHOGENIC ClinVar records -- not patient
# data and not invented coordinates. Each was verified independently in the NCBI
# ClinVar GRCh38 release before being written down here, and each is asserted
# against ClinVar's own molecular-consequence term, so the test fails if SnpEff
# and ClinVar disagree rather than merely if SnpEff changes its output format:
#
#   chr15:40165186 C>T  BUB1B  Pathogenic  MC=SO:0001587|nonsense
#                       ALLELEID 1361261, Mosaic variegated aneuploidy syndrome 1
#   chr11:95812970 C>T  CEP57  Pathogenic  MC=SO:0001587|nonsense
#                       ALLELEID 39649, OMIM 607951.0003, MVA syndrome 2
#
# A previous revision of this file asserted against chr15:40224704 C>T, which is
# not a coding change at all -- it produced no HGVS.p and the test failed the
# moment a real database was present to run it. Coordinates here are checked
# against a source, never recalled.

#: Digests of the installation these expectations were measured against.
#: Hardcoded deliberately: if the database is rebuilt or re-downloaded, the pin
#: fails loudly rather than the assertions below drifting silently to whatever
#: the new gene model says.
REAL_PINS = SnpEffArtifactPins(
    jar="5e8f75cbf908a33c6fb2e65c81e66fe31236cb21bb0195541c4703d8202c22b3",
    config="777768cee885c91c7396ca716c03abb254ef4cde129de2a10513e2844895df5a",
    predictor="f318b79b26ced7e8a44ac18e749ade83bef642e74387163d71885127d85b357a",
    mane_summary="d10ace2720681a3b2e0eefd9da4f551274a6b4141ac9bfd6a2565dfb6e9ad55c",
)

REAL_DATA_DIR = REAL_INSTALL_ROOT / "data"
REAL_JAR = REAL_INSTALL_ROOT / "snpEff" / "snpEff.jar"
REAL_MANE = REAL_INSTALL_ROOT / "mane" / "MANE.GRCh38.v1.5.summary.txt.gz"

BUB1B_NONSENSE = "GRCh38:chr15:40165186:C:T"
CEP57_NONSENSE = "GRCh38:chr11:95812970:C:T"
#: ClinVar Pathogenic 2 bp deletion, MC=SO:0001589|frameshift_variant. Included
#: because a real short-read callset is indel-heavy — a large minority of records
#: carry an indel — which the SNV-only fixtures would not have exercised.
BUB1B_FRAMESHIFT = "GRCh38:chr15:40170578:AGG:A"

REAL_VARIANTS = [BUB1B_NONSENSE, CEP57_NONSENSE, BUB1B_FRAMESHIFT]

requires_real_snpeff = pytest.mark.skipif(
    not (REAL_DATA_DIR / REAL_DATABASE / "snpEffectPredictor.bin").is_file(),
    reason="real SnpEff installation absent; run tools/setup/install_snpeff.sh",
)


def _build_real_adapter(**overrides: Any) -> SnpEffConsequenceAdapter:
    kwargs: dict[str, Any] = {
        "jar_path": REAL_JAR,
        "data_dir": REAL_DATA_DIR,
        "genome_database": REAL_DATABASE,
        "pins": REAL_PINS,
        "mane_summary": REAL_MANE,
    }
    kwargs.update(overrides)
    return SnpEffConsequenceAdapter(**kwargs)


@pytest.fixture(scope="module")
def real_adapter() -> SnpEffConsequenceAdapter:
    """One real adapter for the module: a JVM start plus a 228 MB gene model load.

    No java_binary is passed, so this also proves the $JAVA_HOME resolution path
    works against a genuine JDK -- the machine's `java` on PATH is the macOS stub
    and would fail.
    """
    return _build_real_adapter()


@pytest.fixture(scope="module")
def real_result(
    real_adapter: SnpEffConsequenceAdapter,
) -> tuple[Mapping[str, tuple[ConsequenceAnnotation, ...]], SnpEffRunReport]:
    """One SnpEff invocation shared by every assertion that does not need its own."""
    return real_adapter.annotate_with_report(REAL_VARIANTS)


def _for(annotations: tuple[ConsequenceAnnotation, ...], transcript: str) -> ConsequenceAnnotation:
    """The annotation on one specific transcript, which must be present exactly once."""
    matches = [a for a in annotations if a.transcript_id.split(".", 1)[0] == transcript]
    assert len(matches) == 1, f"expected exactly one {transcript} annotation, got {len(matches)}"
    return matches[0]


@pytest.mark.slow
@requires_real_snpeff
def test_real_snpeff_reports_the_installed_release_and_database(
    real_adapter: SnpEffConsequenceAdapter,
) -> None:
    """The version is read off the artifacts, and it is not just a label."""
    assert real_adapter.version == f"5.4c/{REAL_DATABASE}+{REAL_PINS.composite_digest}"
    assert real_adapter.synthetic is False
    assert is_synthetic(real_adapter) is False
    assert real_adapter.pins == REAL_PINS
    # Resolved from $JAVA_HOME, never from PATH.
    assert real_adapter.java_binary.is_file()
    assert real_adapter.mane_select_count > 15_000, "the real MANE release was not loaded"


@pytest.mark.slow
@requires_real_snpeff
def test_real_snpeff_bub1b_nonsense_hgvs_protein_position_and_exon(
    real_result: tuple[Mapping[str, tuple[ConsequenceAnnotation, ...]], SnpEffRunReport],
) -> None:
    """The fields this adapter exists to deliver, on a real coding variant.

    ClinVar independently classifies this record `MC=SO:0001587|nonsense`, so
    `stop_gained` here is agreement between two sources, not a snapshot of one.
    """
    annotations = real_result[0][BUB1B_NONSENSE]
    mane = _for(annotations, "ENST00000287598")

    assert mane.gene_symbol == "BUB1B"
    assert mane.gene_id == "ENSG00000156970"
    assert mane.consequence_terms == ("stop_gained",)
    assert mane.impact is ImpactSeverity.HIGH
    assert mane.hgvs_c == "c.169C>T"
    assert mane.hgvs_p == "p.Gln57*"
    assert mane.protein_position == 57
    assert mane.exon == "2/23"
    assert mane.intron is None
    assert mane.transcript_biotype == "protein_coding"
    # MANE Select comes from the NCBI release, joined on the unversioned accession:
    # SnpEff emits ENST00000287598.11 and the MANE summary must still match it.
    assert mane.is_mane_select is True
    assert mane.transcript_id.startswith("ENST00000287598.")
    # Base SnpEff computes neither; neither is invented (GP-14).
    assert mane.splice_ai_delta_max is None
    assert mane.pathogenicity_scores == {}


@pytest.mark.slow
@requires_real_snpeff
def test_real_snpeff_cep57_nonsense_hgvs_protein_position_and_exon(
    real_result: tuple[Mapping[str, tuple[ConsequenceAnnotation, ...]], SnpEffRunReport],
) -> None:
    """The second real coding variant, on a different chromosome and gene."""
    annotations = real_result[0][CEP57_NONSENSE]
    mane = _for(annotations, "ENST00000325542")

    assert mane.gene_symbol == "CEP57"
    assert mane.gene_id == "ENSG00000166037"
    assert mane.consequence_terms == ("stop_gained",)
    assert mane.impact is ImpactSeverity.HIGH
    assert mane.hgvs_c == "c.241C>T"
    assert mane.hgvs_p == "p.Arg81*"
    assert mane.protein_position == 81
    assert mane.exon == "3/11"
    assert mane.intron is None
    assert mane.is_mane_select is True


@pytest.mark.slow
@requires_real_snpeff
def test_real_snpeff_annotates_a_real_frameshift_indel(
    real_result: tuple[Mapping[str, tuple[ConsequenceAnnotation, ...]], SnpEffRunReport],
) -> None:
    """Indels are ~19% of the real proband VCF, so they are not an edge case here."""
    annotations = real_result[0][BUB1B_FRAMESHIFT]
    mane = _for(annotations, "ENST00000287598")
    assert mane.gene_symbol == "BUB1B"
    # ClinVar: MC=SO:0001589|frameshift_variant.
    assert "frameshift_variant" in mane.consequence_terms
    assert mane.impact is ImpactSeverity.HIGH
    assert mane.hgvs_c is not None and mane.hgvs_c.startswith("c.")
    assert mane.hgvs_p is not None and mane.hgvs_p.startswith("p.")
    assert mane.protein_position is not None and mane.protein_position >= 1
    assert mane.exon is not None


@pytest.mark.slow
@requires_real_snpeff
def test_real_snpeff_gene_symbol_is_set_on_every_annotation(
    real_result: tuple[Mapping[str, tuple[ConsequenceAnnotation, ...]], SnpEffRunReport],
) -> None:
    """The field the whole pipeline is starved without.

    `VariantRecord.gene_symbols` is derived from these, and `prioritization.pairing`
    groups candidate pairs by them: an annotation with a blank gene symbol removes
    a variant from pairing as effectively as no annotation at all.
    """
    result, report = real_result
    assert set(result) == set(REAL_VARIANTS), "every real variant must be annotated"
    assert report.unannotated == ()
    for variant_id, per_variant in result.items():
        assert per_variant, variant_id
        for annotation in per_variant:
            assert annotation.gene_symbol
            assert annotation.gene_symbol.strip() == annotation.gene_symbol
            assert annotation.transcript_id
            assert annotation.source_tool_version.startswith("5.4c/")
    assert {a.gene_symbol for a in result[BUB1B_NONSENSE]} == {"BUB1B"}
    assert {a.gene_symbol for a in result[CEP57_NONSENSE]} == {"CEP57"}


@pytest.mark.slow
@requires_real_snpeff
def test_real_snpeff_preserves_every_transcript(
    real_result: tuple[Mapping[str, tuple[ConsequenceAnnotation, ...]], SnpEffRunReport],
) -> None:
    """No collapsing to MANE-Select. A real gene has many isoforms and all survive."""
    result = real_result[0]
    bub1b = result[BUB1B_NONSENSE]
    cep57 = result[CEP57_NONSENSE]
    assert len(bub1b) > 5, "transcripts were collapsed"
    assert len(cep57) > 5, "transcripts were collapsed"
    # Exactly one MANE Select each, and the rest are kept alongside it rather than
    # discarded in its favour.
    assert sum(a.is_mane_select for a in bub1b) == 1
    assert sum(a.is_mane_select for a in cep57) == 1
    assert bub1b[0].is_mane_select is True, "MANE Select sorts first"
    # Non-MANE isoforms disagree with MANE about the effect, which is the entire
    # reason collapsing is a data-loss bug.
    assert len({a.most_severe_term for a in cep57}) > 1

    # SnpEff emits more entries than there are transcripts, and that is correct
    # rather than a duplicate: this real variant lands in three NMD transcripts
    # that are reported TWICE each, once as `3_prime_UTR_variant` and once as
    # `non_coding_transcript_exon_variant`, agreeing on transcript, gene, impact,
    # HGVS.c and exon and differing only in the effect term. Both are kept -- the
    # tool said both -- which is why the ordering key cannot stop at the
    # transcript ID.
    assert len({a.transcript_id for a in cep57}) < len(cep57), (
        "expected at least one transcript reported under two effect terms"
    )
    assert len({(a.transcript_id, a.consequence_terms) for a in cep57}) == len(cep57)

    keys = [consequence_sort_key(a) for a in cep57]
    assert keys == sorted(keys)
    assert len(set(keys)) == len(keys), "sort key is not total on real output"


@pytest.mark.slow
@requires_real_snpeff
def test_real_snpeff_repeat_runs_are_byte_identical(
    real_adapter: SnpEffConsequenceAdapter,
    real_result: tuple[Mapping[str, tuple[ConsequenceAnnotation, ...]], SnpEffRunReport],
) -> None:
    """GP-30 against the real subprocess, which is where determinism actually risks.

    A second run through a freshly constructed adapter: a new JVM, a new scratch
    directory, a new temp path, a different wall-clock time. None of that may reach
    the result.
    """
    first = _snapshot(real_result[0])
    assert _snapshot(real_adapter.annotate(REAL_VARIANTS)) == first
    assert _snapshot(_build_real_adapter().annotate(REAL_VARIANTS)) == first
    assert "/tmp" not in first  # noqa: S108 - asserting a temp path is ABSENT
    assert "mva-snpeff-" not in first
    assert "2026-02-23" not in first, "the SnpEff build date leaked into the result"


@pytest.mark.slow
@requires_real_snpeff
def test_real_database_contig_naming_is_verified_not_assumed() -> None:
    """What GRCh38.115 actually calls its chromosomes, checked rather than assumed.

    Measured against this database, not inferred:

    * ``15`` annotates. The database is Ensembl-named.
    * ``chr15`` ALSO annotates -- SnpEff strips a leading ``chr`` before lookup, so
      the autosomes forgive a mismatch.
    * ``chrM`` does **not**. Stripping ``chr`` yields ``M``; the database calls it
      ``MT``. This is the case that makes the explicit mapping load-bearing rather
      than defence in depth, and it fails by silently returning no annotation.

    So `ContigStyle.ENSEMBL` is the correct default, and the reason is the
    mitochondrion, not the autosomes.
    """
    adapter = _build_real_adapter(contig_style=ContigStyle.ENSEMBL)
    mito = "GRCh38:chrM:8993:T:G"
    result, report = adapter.annotate_with_report([BUB1B_NONSENSE, mito])

    # The Ensembl mapping sends MT, and MT-ATP6 comes back. Sent as chrM it would
    # have been ERROR_CHROMOSOME_NOT_FOUND and reported as 'no consequence'.
    assert mito in result, "chrM was not mapped to the database's MT spelling"
    assert any(a.gene_symbol == "MT-ATP6" for a in result[mito])
    assert report.unplaceable == ()

    # And the mapping is asserted directly, not only through a lookup.
    (site,) = plan_batch([mito], build=GenomeBuild.GRCH38, contig_style=ContigStyle.ENSEMBL)
    assert site.contig == "MT"
    (autosome,) = plan_batch(
        [BUB1B_NONSENSE], build=GenomeBuild.GRCH38, contig_style=ContigStyle.ENSEMBL
    )
    assert autosome.contig == "15"


@pytest.mark.slow
@requires_real_snpeff
def test_real_snpeff_argv_is_offline_and_lossless(
    real_adapter: SnpEffConsequenceAdapter,
) -> None:
    """The offline pins, asserted on the argv that really ran (PRIV-05)."""
    argv = real_adapter.build_argv()
    assert "-nodownload" in argv, "SnpEff may fetch a database over the network"
    assert "-noLog" in argv, "SnpEff posts a usage record to its own server"
    assert "-noStats" in argv, "the HTML summary embeds a wall-clock timestamp"
    assert "-canon" not in argv
    assert argv[-1] == REAL_DATABASE
    assert str(real_adapter.java_binary) == argv[0]


# ------------------------------------------------- adversarial review, finding 5
#
# Two halves, both reproduced by execution before they were fixed.
#
# (a) SnpEff can exit ZERO having written a record with no ``ANN`` field at all.
#     ``parse_ann_entries`` already refuses an ANN whose *shape* it cannot read,
#     but a record carrying none was classified into ``SnpEffRunReport.without_ann``
#     and then omitted from the mapping. ``annotate`` is the method the
#     ``ConsequenceAdapter`` Protocol exposes and the only one
#     ``annotation.service`` calls, so that classification was invisible: the
#     caller saw a key missing, which it reads as "this variant has no
#     consequence" -- and since ``VariantRecord.gene_symbols`` is derived from
#     ``consequences`` alone, the variant leaves gene grouping and pairing
#     entirely. Same shape as the gnomAD defect ADR 0018 closed: a report property
#     is something a caller can forget to read.
#
# (b) The installation was not content-pinned in any way a run could actually use.
#     ``SnpEffArtifactPins`` existed and was mandatory, but the only way to obtain
#     one was ``measure()``, which hashes whatever happens to be installed -- so
#     the pin verified the install against itself. ``tools/setup/install_snpeff.sh``
#     writes ``snpeff_pins.json`` *after* checking every digest against the
#     reviewed EXPECT_* constants in the script, and nothing in ``src`` could read
#     it.


@pytest.mark.unit
def test_a_record_with_no_ann_field_fails_closed_through_the_protocol(
    adapter_factory: AdapterFactory,
) -> None:
    """`annotate` refuses output it could not interpret; `annotate_with_report` opts in.

    The stub returns `UNKNOWN_TO_TOOL` as a well-formed VCF row whose INFO carries
    no `ANN=` at all. SnpEff annotates even an intergenic site, so a record with no
    ANN is output this adapter cannot read -- not a finding that the variant has no
    consequence.
    """
    adapter = adapter_factory()
    with pytest.raises(AdapterUnavailableError) as excinfo:
        adapter.annotate([SINGLE_TRANSCRIPT, UNKNOWN_TO_TOOL])
    message = str(excinfo.value)
    assert "no ANN" in message
    assert "annotate_with_report" in message, "the deliberate way through must be named"
    assert "pairing" in message, "the consequence of the silence must be stated"
    # PRIV-09: the shape of the failure, never the record.
    for secret in ("999999", "chr12", "\tA\tC"):
        assert secret not in message, f"{secret!r} leaked into the exception message"

    # The report is still exposed, so the failure is diagnosable.
    report = adapter.last_run_report
    assert report is not None
    assert report.without_ann == (UNKNOWN_TO_TOOL,)

    # And the explicit opt-in still answers, with the gap named in the result.
    result, report = adapter_factory().annotate_with_report([SINGLE_TRANSCRIPT, UNKNOWN_TO_TOOL])
    assert set(result) == {SINGLE_TRANSCRIPT}
    assert report.without_ann == (UNKNOWN_TO_TOOL,)


@pytest.mark.unit
def test_a_clean_batch_is_not_escalated(adapter_factory: AdapterFactory) -> None:
    """The fail-closed check fires on unreadable output, not on a variant SnpEff placed."""
    result = adapter_factory().annotate([SINGLE_TRANSCRIPT, MULTI_TRANSCRIPT])
    assert set(result) == {SINGLE_TRANSCRIPT, MULTI_TRANSCRIPT}


@pytest.mark.unit
def test_pins_can_be_loaded_from_the_installer_manifest(tmp_path: Path) -> None:
    """The reviewed manifest is loadable, so a run can pin to it rather than to itself."""
    manifest = tmp_path / "snpeff_pins.json"
    manifest.write_text(
        json.dumps(
            {
                "snpeff_release": "5.4c",
                "genome_database": "GRCh38.115",
                "jar": "a" * 64,
                "config": "b" * 64,
                "predictor": "c" * 64,
                "mane_summary": "d" * 64,
            }
        ),
        encoding="utf-8",
    )
    pins = SnpEffArtifactPins.from_manifest(manifest)
    assert pins.jar == "a" * 64
    assert pins.config == "b" * 64
    assert pins.predictor == "c" * 64
    assert pins.mane_summary == "d" * 64
    assert pins.snpeff_release == "5.4c"
    assert pins.genome_database == "GRCh38.115"
    # The declared identity is not part of the artifact digest: that is over bytes.
    assert pins.as_rows() == (
        ("jar", "a" * 64),
        ("config", "b" * 64),
        ("predictor", "c" * 64),
        ("mane_summary", "d" * 64),
    )


def _drop_jar(payload: dict[str, Any]) -> None:
    del payload["jar"]


def _drop_database(payload: dict[str, Any]) -> None:
    del payload["genome_database"]


def _corrupt_predictor(payload: dict[str, Any]) -> None:
    payload["predictor"] = "not-a-digest"


def _empty_config(payload: dict[str, Any]) -> None:
    payload["config"] = ""


@pytest.mark.unit
@pytest.mark.parametrize(
    ("mutate", "expected"),
    [
        (_drop_jar, "jar"),
        (_corrupt_predictor, "64 lowercase hex"),
        (_drop_database, "genome_database"),
        (_empty_config, "64 lowercase hex"),
    ],
)
def test_a_manifest_that_pins_less_than_it_claims_is_refused(
    tmp_path: Path, mutate: Callable[[dict[str, Any]], None], expected: str
) -> None:
    """A manifest missing a role would verify three artifacts and report success."""
    payload: dict[str, Any] = {
        "snpeff_release": "5.4c",
        "genome_database": "GRCh38.115",
        "jar": "a" * 64,
        "config": "b" * 64,
        "predictor": "c" * 64,
    }
    mutate(payload)
    manifest = tmp_path / "snpeff_pins.json"
    manifest.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(AdapterUnavailableError) as excinfo:
        SnpEffArtifactPins.from_manifest(manifest)
    assert expected in str(excinfo.value)


@pytest.mark.unit
def test_a_manifest_for_another_database_is_refused_at_construction(
    adapter_factory: AdapterFactory, tmp_path: Path
) -> None:
    """Pins that describe a different genome verify bytes that answer a different question.

    The digests would still match if the operator happened to install both, and
    the run would then stamp one database's provenance onto the other's answers.
    """
    good = adapter_factory().pins
    mismatched = SnpEffArtifactPins(
        jar=good.jar,
        config=good.config,
        predictor=good.predictor,
        snpeff_release=good.snpeff_release,
        genome_database="GRCh37.75",
    )
    with pytest.raises(AdapterUnavailableError) as excinfo:
        adapter_factory(pins=mismatched)
    message = str(excinfo.value)
    assert "GRCh37.75" in message
    assert REAL_DATABASE in message


@pytest.mark.unit
def test_a_manifest_naming_another_snpeff_release_is_refused(
    adapter_factory: AdapterFactory,
) -> None:
    """The JAR that actually runs must be the release the manifest was reviewed for."""
    good = adapter_factory().pins
    mismatched = SnpEffArtifactPins(
        jar=good.jar,
        config=good.config,
        predictor=good.predictor,
        snpeff_release="4.3t",
        genome_database=REAL_DATABASE,
    )
    with pytest.raises(AdapterUnavailableError) as excinfo:
        adapter_factory(pins=mismatched)
    message = str(excinfo.value)
    assert "4.3t" in message


@pytest.mark.unit
def test_the_shipped_pins_manifest_matches_the_expectations_in_this_file() -> None:
    """The manifest the installer wrote and the digests this file asserts on are one thing.

    Skips where the ~1 GB installation is absent. Where it is present, this is the
    check that stops the two drifting: `REAL_PINS` is what the real-SnpEff
    assertions below were measured against, and `snpeff_pins.json` is what a
    composition root will actually load.
    """
    manifest = REAL_INSTALL_ROOT / "snpeff_pins.json"
    if not manifest.is_file():
        pytest.skip(f"no installer manifest at {manifest}")
    pins = SnpEffArtifactPins.from_manifest(manifest, genome_database=REAL_DATABASE)
    assert pins.jar == REAL_PINS.jar
    assert pins.config == REAL_PINS.config
    assert pins.predictor == REAL_PINS.predictor
    assert pins.mane_summary == REAL_PINS.mane_summary
    assert pins.genome_database == REAL_DATABASE


# ------------------------------------------------------ finding 6: build confusion


@pytest.mark.unit
def test_cross_build_annotation_is_refused_through_the_protocol(
    adapter_factory: AdapterFactory,
) -> None:
    """The guard is on the method `annotation.service` calls, not only on `plan_batch`.

    A GRCh37 batch that came back as `{}` would read as "none of these variants has
    a consequence", which removes every one of them from `gene_symbols`, from gene
    grouping and from pairing. One stray build among valid ids must poison the
    batch rather than be dropped from it: a plausible partial answer is the outcome
    a caller cannot detect.
    """
    adapter = adapter_factory()
    with pytest.raises(GenomeBuildMismatchError) as excinfo:
        adapter.annotate(["GRCh37:chr15:40200000:C:T"])
    assert "GRCh37" in str(excinfo.value)

    with pytest.raises(GenomeBuildMismatchError):
        adapter.annotate([SINGLE_TRANSCRIPT, "GRCh37:chr15:40200000:C:T"])

    # PRIV-09: build names, never the coordinate.
    assert "40200000" not in str(excinfo.value)


@pytest.mark.slow
@requires_real_snpeff
def test_the_real_installation_is_pinned_from_the_manifest_it_shipped_with() -> None:
    """End to end: the adapter verifies the bytes it is actually about to run.

    Not `measure()`, which hashes whatever is installed and therefore pins the
    installation to itself. `tools/setup/install_snpeff.sh` checks every artifact
    against the reviewed EXPECT_* digests in the script and only then writes
    `snpeff_pins.json`, so loading that file is a verification rather than a
    recording. The release string in the manifest is checked against what the JAR
    says about itself, which is the half a digest cannot cover.
    """
    manifest = REAL_INSTALL_ROOT / "snpeff_pins.json"
    if not manifest.is_file():
        pytest.skip(f"no installer manifest at {manifest}")
    pins = SnpEffArtifactPins.from_manifest(manifest, genome_database=REAL_DATABASE)
    adapter = _build_real_adapter(pins=pins)
    assert adapter.pins.genome_database == REAL_DATABASE
    assert adapter.pins.snpeff_release == "5.4c"
    assert adapter.version == f"5.4c/{REAL_DATABASE}+{pins.composite_digest}"
    # The digests are over bytes, so a manifest-loaded pin and a hardcoded one
    # produce the identical provenance string.
    assert pins.composite_digest == REAL_PINS.composite_digest
    result = adapter.annotate([BUB1B_NONSENSE])
    assert result[BUB1B_NONSENSE]


# --------------------------------------------------------- one launch per chunk
#
# SnpEff's cost is ~100% fixed start-up: measured on this machine, 35.5 s to
# annotate ONE variant and 31.6 s to annotate 5,000, because both are a
# GRCh38.115 database load. The marginal per-variant cost is indistinguishable
# from zero, so the number of JVM launches IS the runtime, and the batch size is
# the only lever on it.
#
# The 4,962,060-variant callset projected to ~30 h. `service.DEFAULT_ANNOTATION_
# BATCH_SIZE` was 5,000 while `snpeff_local.DEFAULT_BATCH_SIZE` was 25,000, so the
# adapter's own chunking never triggered and the run paid one launch per 5,000
# variants: 993 launches, ~9.8 h of nothing but reloading the same database.
#
# These tests guard the two halves of the fix that are silent when broken -- a
# mismatched pair of constants, and a return to buffering the output in memory.


def test_the_two_batch_sizes_are_equal() -> None:
    """A mismatch silently keeps most of the reloads, and nothing else notices.

    Whichever constant is smaller decides the launch count. 5,000 against 25,000
    meant the adapter's chunking never fired at all; 250,000 against a stale
    25,000 would mean ten launches per service batch and 199 of the original 993
    surviving a change that looks, from the outside, like it worked.
    """
    assert DEFAULT_ANNOTATION_BATCH_SIZE == DEFAULT_BATCH_SIZE, (
        "the service batch and the SnpEff chunk disagree; the smaller one wins and "
        "the JVM reloads come back"
    )


def test_the_batch_size_is_the_one_the_machine_can_afford() -> None:
    """Both bounds, asserted, because both are real and they point opposite ways.

    Too small and the run pays for a database reload it does not need; too large
    and it trips the memory watchdog. 500,000 is deliberately not taken: the scale
    sweep puts batch-1M at 7.79 GiB against a 12 GiB limit, and at 20 launches the
    fixed cost is already amortised, so there is no throughput left to buy.
    """
    callset = 4_962_060
    launches = math.ceil(callset / DEFAULT_ANNOTATION_BATCH_SIZE)
    assert launches <= 25, (
        f"{launches} JVM launches over the real callset; each is a full GRCh38.115 "
        "database load and the run has no time for them"
    )
    assert DEFAULT_ANNOTATION_BATCH_SIZE <= 250_000, (
        "24 GB of RAM and a 6 GB JVM heap; the scale sweep's batch-1M measured "
        "7.79 GiB against a 12 GiB watchdog"
    )


def test_each_chunk_costs_exactly_one_subprocess(
    adapter_factory: AdapterFactory, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Chunk count and launch count are the same number, which is the whole premise."""
    adapter = adapter_factory(batch_size=2)
    variants = [MULTI_TRANSCRIPT, SINGLE_TRANSCRIPT, INTRONIC, FRAMESHIFT, INTERGENIC]

    launches: list[tuple[object, ...]] = []
    real_run = subprocess.run

    def _counting(*args: Any, **kwargs: Any) -> Any:
        launches.append(args)
        return real_run(*args, **kwargs)

    # Patched after construction so the version probe's own launch is not counted.
    monkeypatch.setattr(snpeff_local.subprocess, "run", _counting)
    adapter.annotate(variants)

    assert len(launches) == math.ceil(len(variants) / 2) == 3


def test_annotation_output_is_spooled_to_a_file_not_held_in_memory(
    adapter_factory: AdapterFactory, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The reason the 25,000 ceiling existed, removed rather than raised.

    `_invoke` used `capture_output=True`, which accumulates the whole annotated VCF
    as one `bytes` object and then again as the `str` that `.splitlines()` was
    called on. At ~700 bytes per annotated record that is ~350 MB per 250,000
    variants, before the list of lines -- and that memory, not any property of
    SnpEff, is what the old batch ceiling encoded. Raising the batch size without
    removing it would trade JVM reloads for an OOM.

    Asserted structurally: a test that only measured wall time would pass against a
    reintroduced buffer on a five-variant fixture.
    """
    adapter = adapter_factory()
    captured: list[dict[str, Any]] = []
    real_run = subprocess.run

    def _recording(*args: Any, **kwargs: Any) -> Any:
        captured.append(kwargs)
        return real_run(*args, **kwargs)

    monkeypatch.setattr(snpeff_local.subprocess, "run", _recording)
    adapter.annotate([MULTI_TRANSCRIPT, SINGLE_TRANSCRIPT])

    (kwargs,) = captured
    assert not kwargs.get("capture_output"), (
        "stdout is being buffered in memory again; the batch ceiling this removed "
        "would have to come back with it"
    )
    stdout = kwargs["stdout"]
    assert hasattr(stdout, "fileno"), "stdout is not a real file the child writes into"
    # Inside the run's own scratch directory: the spool holds annotated proband
    # coordinates and must not outlive the call or land somewhere shared (GP-40).
    assert "mva-snpeff-" in stdout.name
    assert kwargs["stderr"] is not stdout, "stderr merged into the annotated output stream"


def test_the_scratch_spool_does_not_outlive_the_call(
    adapter_factory: AdapterFactory, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The annotated VCF carries proband coordinates; it is unlinked with the scratch."""
    adapter = adapter_factory()
    paths: list[Path] = []
    real_run = subprocess.run

    def _recording(*args: Any, **kwargs: Any) -> Any:
        paths.append(Path(kwargs["stdout"].name))
        return real_run(*args, **kwargs)

    monkeypatch.setattr(snpeff_local.subprocess, "run", _recording)
    adapter.annotate([MULTI_TRANSCRIPT])

    (spool,) = paths
    assert not spool.exists(), "the annotated output was left on disk after the call"
    assert not spool.parent.exists(), "the scratch directory outlived the call"
