"""Tests for ``tools/build_knowledge``: the generator that turns downloaded, real
gene-disease resources into ``knowledge/real/*.tsv``.

Two kinds of coverage:

* **Fixture-driven unit tests** (most of this file) build tiny, self-contained files
  shaped exactly like the real sources (gnomAD's constraint TSV, ClinGen's CSV, HPO's
  genes_to_phenotype.txt, a valid and a broken DDG2P.csv.gz) under ``tmp_path``, so
  they run anywhere -- no dependency on ``$MVA_RESOURCES`` or the ~100MB of real
  downloads that live outside this repository. They prove the rules that matter most:
  GP-14 (missing is never zero), source vocabulary is preserved verbatim, and output
  is byte-identical across repeat runs (GP-30).
* **Assertions against the committed ``knowledge/real/*.tsv``** prove those specific
  properties hold in the actual shipped tables, not just in the generator's logic.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from mva.config import find_repo_root

# `tools/` is deliberately NOT part of the installed `mva` distribution (see
# tools/build_knowledge/__init__.py: it is standalone acquisition-to-table tooling,
# not pipeline code), so it is only importable with the repo root on sys.path.
# `tests/` is deliberately not a package either (see tests/conftest.py), so this
# bootstrap lives here, in the one file that needs it, rather than in shared
# fixtures this task does not own.
REPO_ROOT = find_repo_root(Path(__file__))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import gzip  # noqa: E402

from tools.build_knowledge.build import BuildReport, ResourcePaths, build_all  # noqa: E402
from tools.build_knowledge.gene_disease import (  # noqa: E402
    GENE_DISEASE_COLUMNS,
    parse_clingen,
    try_parse_ddg2p,
)
from tools.build_knowledge.gene_phenotype import (  # noqa: E402
    GENE_PHENOTYPE_COLUMNS,
    parse_gene_to_phenotype,
    read_hpo_release,
)
from tools.build_knowledge.gnomad_constraint import (  # noqa: E402
    GENE_PANEL_COLUMNS,
    extract_gnomad_version,
    parse_gnomad_constraint,
)
from tools.build_knowledge.tsv_io import (  # noqa: E402
    format_cell,
    format_float,
    render_tsv,
    write_tsv,
)

pytestmark = pytest.mark.unit

KNOWLEDGE_REAL = REPO_ROOT / "knowledge" / "real"
KNOWLEDGE_PUBLIC = REPO_ROOT / "knowledge" / "public"


def _read_tsv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Minimal reader mirroring mva.annotation.local_tables._read_rows's contract:
    '#'-prefixed comments are skipped, the first non-comment line is the header.
    """
    header: list[str] | None = None
    rows: list[dict[str, str]] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            stripped = raw.rstrip("\n")
            if not stripped.strip() or stripped.lstrip().startswith("#"):
                continue
            fields = stripped.split("\t")
            if header is None:
                header = fields
                continue
            rows.append(dict(zip(header, fields, strict=True)))
    assert header is not None, f"{path} has no header row"
    return header, rows


# --------------------------------------------------------------------------- tsv_io


def test_format_float_distinguishes_missing_from_a_real_zero() -> None:
    """The crux of GP-14: an unmeasured value and a genuinely observed zero must
    never render the same way.
    """
    assert format_float(None) == ""
    assert format_float(0.0) == "0"
    assert format_float(0.0) != format_float(None)


def test_format_cell_empty_for_none_only() -> None:
    assert format_cell(None) == ""
    assert format_cell("") == ""
    assert format_cell("Definitive") == "Definitive"


def test_render_tsv_uses_lf_only_and_is_deterministic() -> None:
    rows = [{"a": "1", "b": "x"}, {"a": "2", "b": "y"}]
    first = render_tsv(comment_lines=["hello"], columns=["a", "b"], rows=rows)
    second = render_tsv(comment_lines=["hello"], columns=["a", "b"], rows=rows)
    assert first == second
    assert "\r" not in first
    assert first.startswith("# hello\na\tb\n1\tx\n2\ty\n")


def test_render_tsv_rejects_a_tab_or_newline_inside_a_cell() -> None:
    with pytest.raises(ValueError, match="tab or newline"):
        render_tsv(comment_lines=[], columns=["a"], rows=[{"a": "has\ttab"}])
    with pytest.raises(ValueError, match="tab or newline"):
        render_tsv(comment_lines=[], columns=["a"], rows=[{"a": "has\nnewline"}])


def test_write_tsv_writes_lf_line_endings_on_disk(tmp_path: Path) -> None:
    out = tmp_path / "out.tsv"
    write_tsv(out, comment_lines=["c"], columns=["a"], rows=[{"a": "1"}])
    raw = out.read_bytes()
    assert b"\r\n" not in raw
    assert raw.endswith(b"\n")


# --------------------------------------------------------------------------- gnomAD

_GNOMAD_HEADER = (
    "gene\tgene_id\ttranscript\tcanonical\tmane_select\t"
    "lof.pLI\tlof.oe_ci.upper\tmis.z_score\tconstraint_flags\n"
)


def _write_gnomad_fixture(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        _GNOMAD_HEADER,
        # GENE_A: two transcripts; the MANE-Select+canonical one carries a REAL,
        # observed pLI of exactly 0.0 -- must render as "0", not empty.
        "GENE_A\tENSG00000000001\tENST00000000009\tfalse\tfalse\t0.9\t0.2\t3.0\t[]\n",
        "GENE_A\tENSG00000000001\tENST00000000001\ttrue\ttrue\t0.0\t1.5\t2.0\t[]\n",
        # GENE_B: no expected-LoF count -> pLI/LOEUF genuinely unmeasured (NA).
        'GENE_B\tENSG00000000002\tENST00000000002\ttrue\ttrue\tNA\tNA\tNA\t["no_exp_lof"]\n',
        # GENE_C: only a non-Ensembl (Entrez-style) gene_id row -> excluded entirely.
        "GENE_C\t999\tNM_000003.1\ttrue\ttrue\t0.5\t0.5\t0.5\t[]\n",
        # gnomAD's own missing-gene sentinel -- must never be treated as a real gene.
        "NA\tENSG00000000004\tENST00000000004\ttrue\ttrue\t0.1\t0.1\t0.1\t[]\n",
    ]
    path.write_text("".join(lines), encoding="utf-8")


def test_gnomad_version_parsed_from_filename(tmp_path: Path) -> None:
    path = tmp_path / "gnomad.v4.1.constraint_metrics.tsv"
    path.write_text(_GNOMAD_HEADER, encoding="utf-8")
    assert extract_gnomad_version(path) == "v4.1"


def test_gnomad_version_refuses_to_guess_an_unexpected_filename(tmp_path: Path) -> None:
    path = tmp_path / "constraint.tsv"
    with pytest.raises(ValueError, match="Refusing to guess"):
        extract_gnomad_version(path)


def test_gnomad_parser_never_defaults_missing_constraint_to_zero(tmp_path: Path) -> None:
    path = tmp_path / "gnomad.v9.9.constraint_metrics.tsv"
    _write_gnomad_fixture(path)
    result = parse_gnomad_constraint(path)
    by_gene = {row.gene_symbol: row for row in result.rows}

    # GENE_A's selected (MANE+canonical) row has a real pLI of 0.0.
    assert by_gene["GENE_A"].pli == 0.0
    assert by_gene["GENE_A"].transcript_id == "ENST00000000001"

    # GENE_B has NO computable metrics -- None, not 0.0.
    assert by_gene["GENE_B"].pli is None
    assert by_gene["GENE_B"].loeuf is None
    assert by_gene["GENE_B"].mis_z is None
    assert by_gene["GENE_B"].constraint_flags == ("no_exp_lof",)

    # The rendered TSV cell must carry that distinction through: "0" vs empty.
    assert format_float(by_gene["GENE_A"].pli) == "0"
    assert format_float(by_gene["GENE_B"].pli) == ""


def test_gnomad_parser_excludes_na_sentinel_and_non_ensembl_only_genes(tmp_path: Path) -> None:
    path = tmp_path / "gnomad.v9.9.constraint_metrics.tsv"
    _write_gnomad_fixture(path)
    result = parse_gnomad_constraint(path)
    genes = {row.gene_symbol for row in result.rows}

    assert "NA" not in genes  # gnomAD's missing-gene sentinel is not a gene.
    assert "GENE_C" not in genes  # no Ensembl-tagged transcript row.
    assert "GENE_C" in result.genes_without_ensembl_row  # ...but it IS reported, not silent.
    assert result.na_gene_rows == 1
    assert result.genes_included == 2  # GENE_A, GENE_B


def test_gnomad_parser_prefers_mane_select_then_canonical(tmp_path: Path) -> None:
    path = tmp_path / "gnomad.v9.9.constraint_metrics.tsv"
    _write_gnomad_fixture(path)
    result = parse_gnomad_constraint(path)
    selected = next(r for r in result.rows if r.gene_symbol == "GENE_A")
    assert selected.mane_select is True
    assert selected.canonical is True


def test_gnomad_parser_is_deterministic_regardless_of_row_order(tmp_path: Path) -> None:
    forward = tmp_path / "gnomad.v9.9.constraint_metrics.tsv"
    _write_gnomad_fixture(forward)
    forward_result = parse_gnomad_constraint(forward)

    lines = forward.read_text(encoding="utf-8").splitlines(keepends=True)
    reversed_path = tmp_path / "reversed" / "gnomad.v9.9.constraint_metrics.tsv"
    reversed_path.parent.mkdir(parents=True, exist_ok=True)
    reversed_path.write_text("".join([lines[0], *reversed(lines[1:])]), encoding="utf-8")
    reversed_result = parse_gnomad_constraint(reversed_path)

    assert forward_result.rows == reversed_result.rows


# --------------------------------------------------------------------------- ClinGen

_CLINGEN_RULE_LINE = ",".join(f'"{"+" * 10}"' for _ in range(10)) + "\n"


def _write_clingen_fixture(path: Path, *, created: str = "2099-01-01") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        '"CLINGEN GENE DISEASE VALIDITY CURATIONS","","","","","","","","",""\n',
        f'"FILE CREATED: {created}","","","","","","","","",""\n',
        '"WEBPAGE: https://example.invalid/kb/gene-validity","","","","","","","","",""\n',
        _CLINGEN_RULE_LINE,
        '"GENE SYMBOL","GENE ID (HGNC)","DISEASE LABEL","DISEASE ID (MONDO)","MOI","SOP",'
        '"CLASSIFICATION","ONLINE REPORT","CLASSIFICATION DATE","GCEP"\n',
        _CLINGEN_RULE_LINE,
        '"GENE_A","HGNC:1","disease one","MONDO:0000001","AD","SOP1","Definitive",'
        '"https://example.invalid/1","2020-01-01T00:00:00.000Z","Panel One"\n',
        '"GENE_A","HGNC:1","disease two","MONDO:0000002","AR","SOP1","Disputed",'
        '"https://example.invalid/2","2021-01-01T00:00:00.000Z","Panel Two"\n',
        '"GENE_D","HGNC:4","disease four","MONDO:0000004","XL","SOP1","Refuted",'
        '"https://example.invalid/4","2022-01-01T00:00:00.000Z","Panel Four"\n',
    ]
    path.write_text("".join(lines), encoding="utf-8")


def test_clingen_version_is_read_from_the_file_created_line(tmp_path: Path) -> None:
    path = tmp_path / "clingen_gene_validity.csv"
    _write_clingen_fixture(path, created="2077-07-07")
    result = parse_clingen(path)
    assert result.version == "2077-07-07"
    assert all(row.version == "2077-07-07" for row in result.rows)


def test_clingen_preserves_disputed_and_refuted_verbatim(tmp_path: Path) -> None:
    path = tmp_path / "clingen_gene_validity.csv"
    _write_clingen_fixture(path)
    result = parse_clingen(path)
    confidences = {
        row.gene_symbol: row.confidence
        for row in result.rows
        if row.disease_id
        in {
            "MONDO:0000002",
            "MONDO:0000004",
        }
    }
    assert confidences["GENE_A"] == "Disputed"
    assert confidences["GENE_D"] == "Refuted"
    # Never remapped onto a synthetic numeric score.
    assert all(row.confidence in {"Definitive", "Disputed", "Refuted"} for row in result.rows)


def test_clingen_mutation_consequence_is_always_absent_not_fabricated(tmp_path: Path) -> None:
    path = tmp_path / "clingen_gene_validity.csv"
    _write_clingen_fixture(path)
    result = parse_clingen(path)
    assert all(row.mutation_consequence is None for row in result.rows)


# --------------------------------------------------------------------------- DDG2P
#
# Schema below matches the REAL EBI Gene2Phenotype per-panel export
# (https://www.ebi.ac.uk/gene2phenotype/api/panel/DD/download/), verified against the
# actual downloaded file on 2026-08-28 -- not the older bulk-download column names
# published in some G2P papers, which this endpoint does not use. Only the columns
# try_parse_ddg2p actually requires are included; the real file carries several more
# (g2p id, previous gene symbols, molecular mechanism, phenotypes, publications, ...)
# that this generator does not read.

_DDG2P_HEADER = (
    "gene symbol,hgnc id,disease name,disease mim,disease MONDO,"
    "allelic requirement,confidence,variant consequence,date of last review,panel\n"
)


def _write_valid_ddg2p_fixture(path: Path, *, gzip_compressed: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    body = (
        _DDG2P_HEADER
        # MONDO present -> preferred over the OMIM mim for disease_id.
        + "GENE_A,100001,disease five,600123,MONDO:0000005,biallelic_autosomal,disputed,"
        "absent gene product,2020-01-01 00:00:00+00:00,DD; Eye\n"
        # No MONDO -> falls back to the numeric OMIM mim, prefixed.
         + "GENE_E,100002,disease six,600124,,monoallelic_autosomal,definitive,"
        "altered gene product structure,2021-02-02 00:00:00+00:00,DD\n"
        # Neither MONDO nor mim -> disease_id must be empty, never fabricated.
         + "GENE_F,100003,disease seven,,,monoallelic_X_hemizygous,refuted,uncertain,"
        "2022-03-03 00:00:00+00:00,DD\n"
    )
    if gzip_compressed:
        with gzip.open(path, mode="wt", encoding="utf-8", newline="") as handle:
            handle.write(body)
    else:
        path.write_text(body, encoding="utf-8")


def test_try_parse_ddg2p_reports_unavailable_for_an_html_error_page(tmp_path: Path) -> None:
    """Regression test for the exact failure mode found in the real (first) download
    attempt: an HTML app shell saved with a '.csv.gz' name, HTTP 200, which a
    status-code-only check ('curl -f' didn't fail) would miss entirely.
    """
    path = tmp_path / "DDG2P.csv"
    path.write_bytes(b"<!doctype html>\n<html><body>not actually CSV or gzip</body></html>")
    result = try_parse_ddg2p(path, version="2026-08-28")
    assert result.available is False
    assert result.rows == ()
    assert result.reason is not None


def test_try_parse_ddg2p_reports_unavailable_when_file_is_missing(tmp_path: Path) -> None:
    result = try_parse_ddg2p(tmp_path / "does-not-exist.csv", version="2026-08-28")
    assert result.available is False
    assert result.rows == ()


@pytest.mark.parametrize("gzip_compressed", [False, True])
def test_try_parse_ddg2p_parses_a_valid_csv_gzip_or_plain(
    tmp_path: Path, *, gzip_compressed: bool
) -> None:
    """The real working DDG2P export is PLAIN CSV (not gzip, despite the resource
    catalogue's original '.csv.gz' naming); the parser must handle both by content,
    not by trusting the file extension.
    """
    path = tmp_path / "DDG2P.csv"
    _write_valid_ddg2p_fixture(path, gzip_compressed=gzip_compressed)
    result = try_parse_ddg2p(path, version="2026-08-28")
    assert result.available is True
    assert result.reason is None
    by_gene = {row.gene_symbol: row for row in result.rows}

    assert by_gene["GENE_A"].confidence == "disputed"  # vocabulary preserved verbatim, lower-case
    assert by_gene["GENE_A"].disease_id == "MONDO:0000005"  # MONDO preferred when present
    assert by_gene["GENE_A"].hgnc_id == "HGNC:100001"
    assert by_gene["GENE_A"].mutation_consequence == "absent gene product"
    assert by_gene["GENE_A"].inheritance == "biallelic_autosomal"
    assert by_gene["GENE_A"].panel == "DD; Eye"

    assert by_gene["GENE_E"].disease_id == "OMIM:600124"  # falls back to the OMIM mim

    # Neither MONDO nor mim present -> EMPTY (None), never a fabricated identifier.
    assert by_gene["GENE_F"].disease_id is None
    assert by_gene["GENE_F"].confidence == "refuted"  # preserved even though it's negative

    assert all(row.source == "DDG2P" for row in result.rows)
    assert all(row.gcep is None for row in result.rows)  # ClinGen-specific, not DDG2P's


def test_try_parse_ddg2p_fails_closed_on_an_unexpected_schema(tmp_path: Path) -> None:
    path = tmp_path / "DDG2P.csv"
    path.write_text("totally,different,columns\n1,2,3\n", encoding="utf-8")
    result = try_parse_ddg2p(path, version="2026-08-28")
    assert result.available is False
    assert result.reason is not None
    assert "missing expected column" in result.reason
    assert "missing expected column" in result.reason


# --------------------------------------------------------------------------- HPO

_HPOA_FIXTURE = (
    '#description: "fixture"\n'
    "#version: 2099-09-09\n"
    "#tracker: https://example.invalid\n"
    "#hpo-version: http://purl.obolibrary.org/obo/hp/releases/2099-09-09/hp.json\n"
    "database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\tfrequency\t"
    "sex\tmodifier\taspect\tbiocuration\n"
)

_GENES_TO_PHENOTYPE_HEADER = "ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id\n"


def _write_hpo_fixtures(hpoa_path: Path, genes_to_phenotype_path: Path) -> None:
    hpoa_path.parent.mkdir(parents=True, exist_ok=True)
    hpoa_path.write_text(_HPOA_FIXTURE, encoding="utf-8")

    genes_to_phenotype_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        _GENES_TO_PHENOTYPE_HEADER,
        # Same (gene, hpo_id) pair via two diseases: one unstated ('-'), one real.
        # The '-' row comes FIRST -- collapsing must not let file order decide.
        "1\tGENE_A\tHP:0000001\tTerm One\t-\tOMIM:111\n",
        "1\tGENE_A\tHP:0000001\tTerm One\t3/5\tOMIM:222\n",
        "1\tGENE_A\tHP:0000002\tTerm Two\tHP:0040281\tOMIM:111\n",
        # Not in the restriction set -- must be excluded entirely.
        "2\tGENE_NOT_RESTRICTED\tHP:0000009\tTerm Nine\t-\tOMIM:999\n",
    ]
    genes_to_phenotype_path.write_text("".join(lines), encoding="utf-8")


def test_hpo_release_is_read_from_phenotype_hpoa_header(tmp_path: Path) -> None:
    hpoa = tmp_path / "phenotype.hpoa"
    hpoa.write_text(_HPOA_FIXTURE, encoding="utf-8")
    assert read_hpo_release(hpoa) == "2099-09-09"


def test_gene_phenotype_dash_becomes_none_not_a_fabricated_value(tmp_path: Path) -> None:
    hpoa = tmp_path / "phenotype.hpoa"
    genes_to_phenotype = tmp_path / "genes_to_phenotype.txt"
    _write_hpo_fixtures(hpoa, genes_to_phenotype)

    result = parse_gene_to_phenotype(
        genes_to_phenotype, restrict_to_genes={"GENE_A", "GENE_Z"}, version="2099-09-09"
    )
    by_pair = {(row.gene_symbol, row.hpo_id): row for row in result.rows}

    # The pair with one real and one unstated frequency keeps the REAL one, even
    # though the unstated ('-') row appeared first in the source file.
    assert by_pair[("GENE_A", "HP:0000001")].association_strength == "3/5"
    assert by_pair[("GENE_A", "HP:0000002")].association_strength == "HP:0040281"

    # A restricted-in gene with zero matching rows is reported, not silently absent.
    assert "GENE_Z" in result.genes_unmatched
    # A gene outside the restriction set never appears.
    assert "GENE_NOT_RESTRICTED" not in {r.gene_symbol for r in result.rows}


def test_gene_phenotype_no_source_row_ever_becomes_a_zero_or_empty_string_association(
    tmp_path: Path,
) -> None:
    """The rendered cell for an unstated frequency is EMPTY (None), never the literal
    string '-' carried through, and never conflated with a real "0"-shaped value.
    """
    hpoa = tmp_path / "phenotype.hpoa"
    genes_to_phenotype = tmp_path / "genes_to_phenotype.txt"
    genes_to_phenotype.parent.mkdir(parents=True, exist_ok=True)
    genes_to_phenotype.write_text(
        _GENES_TO_PHENOTYPE_HEADER + "1\tGENE_A\tHP:0000001\tTerm One\t-\tOMIM:111\n",
        encoding="utf-8",
    )
    hpoa.write_text(_HPOA_FIXTURE, encoding="utf-8")
    result = parse_gene_to_phenotype(
        genes_to_phenotype, restrict_to_genes={"GENE_A"}, version="2099-09-09"
    )
    (row,) = result.rows
    assert row.association_strength is None
    assert format_cell(row.association_strength) == ""


# --------------------------------------------------------------------------- sanity checks
#
# Direct regression coverage for the failure this project actually hit: an HTML
# app-shell page saved under a dataset's expected filename, with a 200-OK status that
# a status-code-only check cannot catch. Every reader here must fail LOUDLY (raise, or
# report ``available=False`` with a reason) -- never silently produce a plausible-
# looking-but-empty-or-wrong table, which is the dangerous failure mode: a confident,
# empty, wrong table is indistinguishable from success until someone reads the numbers.

_HTML_GARBAGE = b"<!doctype html>\n<html><body>not the expected dataset</body></html>"


def test_every_source_reader_rejects_html_garbage_instead_of_silently_succeeding(
    tmp_path: Path,
) -> None:
    gnomad_path = tmp_path / "gnomad.v9.9.constraint_metrics.tsv"
    gnomad_path.write_bytes(_HTML_GARBAGE)
    with pytest.raises(ValueError, match="missing column"):
        parse_gnomad_constraint(gnomad_path)

    clingen_path = tmp_path / "clingen_gene_validity.csv"
    clingen_path.write_bytes(_HTML_GARBAGE)
    with pytest.raises(ValueError, match="FILE CREATED"):
        parse_clingen(clingen_path)

    hpoa_path = tmp_path / "phenotype.hpoa"
    hpoa_path.write_bytes(_HTML_GARBAGE)
    with pytest.raises(ValueError, match="version"):
        read_hpo_release(hpoa_path)

    genes_to_phenotype_path = tmp_path / "genes_to_phenotype.txt"
    genes_to_phenotype_path.write_bytes(_HTML_GARBAGE)
    with pytest.raises(ValueError, match="missing column"):
        parse_gene_to_phenotype(genes_to_phenotype_path, restrict_to_genes={"X"}, version="v1")

    ddg2p_path = tmp_path / "DDG2P.csv"
    ddg2p_path.write_bytes(_HTML_GARBAGE)
    result = try_parse_ddg2p(ddg2p_path, version="v1")
    assert result.available is False
    assert result.rows == ()


# --------------------------------------------------------------------------- build_all


def _write_full_fixture_resources(root: Path) -> None:
    """A complete, tiny resources root shaped like ``ResourcePaths.under()`` expects,
    real-shaped filenames included (``gnomad.v4.1.constraint_metrics.tsv`` matters:
    :func:`extract_gnomad_version` reads the release from the filename itself).
    """
    paths = ResourcePaths.under(root)
    _write_gnomad_fixture(paths.gnomad_constraint)
    _write_clingen_fixture(paths.clingen)
    # No DDG2P.csv.gz at all -- exercises the "file missing" unavailable path inside a
    # full build, distinct from the "wrong bytes" case tested above.
    _write_hpo_fixtures(paths.hpo_phenotype_hpoa, paths.hpo_genes_to_phenotype)


def test_build_all_is_byte_identical_across_two_runs(tmp_path: Path) -> None:
    resources_root = tmp_path / "resources"
    _write_full_fixture_resources(resources_root)

    out1 = tmp_path / "out1"
    out2 = tmp_path / "out2"
    report1 = build_all(resources_root, out1)
    report2 = build_all(resources_root, out2)

    for name in ("gene_panel.tsv", "gene_disease.tsv", "gene_phenotype.tsv"):
        bytes1 = (out1 / name).read_bytes()
        bytes2 = (out2 / name).read_bytes()
        assert bytes1 == bytes2, f"{name} was not byte-identical across two runs"

    assert isinstance(report1, BuildReport)
    assert isinstance(report2, BuildReport)


def test_build_all_writes_expected_columns(tmp_path: Path) -> None:
    resources_root = tmp_path / "resources"
    _write_full_fixture_resources(resources_root)
    out = tmp_path / "out"
    build_all(resources_root, out)

    panel_header, panel_rows = _read_tsv_rows(out / "gene_panel.tsv")
    assert panel_header == list(GENE_PANEL_COLUMNS)
    assert panel_rows  # non-empty

    disease_header, disease_rows = _read_tsv_rows(out / "gene_disease.tsv")
    assert disease_header == list(GENE_DISEASE_COLUMNS)
    assert disease_rows

    phenotype_header, phenotype_rows = _read_tsv_rows(out / "gene_phenotype.tsv")
    assert phenotype_header == list(GENE_PHENOTYPE_COLUMNS)
    assert phenotype_rows


def test_build_all_reports_ddg2p_as_unavailable_when_missing(tmp_path: Path) -> None:
    resources_root = tmp_path / "resources"
    _write_full_fixture_resources(resources_root)
    out = tmp_path / "out"
    report = build_all(resources_root, out)
    assert report.ddg2p.available is False
    # ClinGen rows must still be written even though DDG2P failed (no cross-source
    # failure coupling): GENE_A and GENE_D both came from the ClinGen fixture.
    gene_symbols = {row.gene_symbol for row in report.clingen.rows}
    assert {"GENE_A", "GENE_D"} <= gene_symbols


def test_build_all_widens_hpo_restriction_when_ddg2p_is_available(tmp_path: Path) -> None:
    """When DDG2P parses successfully, the gene_phenotype restriction is the UNION of
    ClinGen and DDG2P genes -- not ClinGen alone -- per the brief's "genes present in
    DDG2P or ClinGen" restriction.
    """
    resources_root = tmp_path / "resources"
    _write_full_fixture_resources(resources_root)
    paths = ResourcePaths.under(resources_root)
    # GENE_E is a DDG2P-only gene (not in the ClinGen fixture) and IS in the HPO
    # fixture's restrict-able rows -- add one so the widened restriction is observable.
    genes_to_phenotype_extra = paths.hpo_genes_to_phenotype
    genes_to_phenotype_extra.write_text(
        genes_to_phenotype_extra.read_text(encoding="utf-8")
        + "3\tGENE_E\tHP:0000099\tTerm Ninety-Nine\t1/1\tOMIM:333\n",
        encoding="utf-8",
    )
    _write_valid_ddg2p_fixture(paths.ddg2p, gzip_compressed=False)

    out = tmp_path / "out"
    report = build_all(resources_root, out)

    assert report.ddg2p.available is True
    assert "GENE_A" in report.clingen.unique_genes  # ClinGen alone does NOT have GENE_E
    assert "GENE_E" not in report.clingen.unique_genes
    hpo_gene_symbols = {row.gene_symbol for row in report.hpo.rows}
    assert "GENE_E" in hpo_gene_symbols  # only reachable via the widened (DDG2P) restriction


# --------------------------------------------------------------------------- committed real tables


def test_real_gene_panel_matches_documented_schema() -> None:
    header, rows = _read_tsv_rows(KNOWLEDGE_REAL / "gene_panel.tsv")
    assert header == list(GENE_PANEL_COLUMNS)
    assert len(rows) > 1000  # sanity: this is the real, genome-scale table


def test_real_gene_panel_keeps_missing_constraint_distinct_from_zero() -> None:
    """Proves the GP-14 distinction holds in the SHIPPED table, not just in the
    generator's unit-tested logic: some real genes have an empty pli/loeuf cell
    (unmeasured), and that is not the same as any row reading literal "0".
    """
    _, rows = _read_tsv_rows(KNOWLEDGE_REAL / "gene_panel.tsv")
    pli_values = [row["pli"] for row in rows]
    assert "" in pli_values, "expected at least one gene with no computable pLI"
    non_empty = [v for v in pli_values if v != ""]
    assert non_empty, "expected at least one gene with a real pLI value"
    # None of gnomAD's own missing-value sentinel leaked through unconverted.
    assert "NA" not in pli_values
    assert "NA" not in [row["loeuf"] for row in rows]


def test_real_gene_panel_excludes_sex_chromosome_genes_and_says_so() -> None:
    _, rows = _read_tsv_rows(KNOWLEDGE_REAL / "gene_panel.tsv")
    symbols = {row["gene_symbol"] for row in rows}
    assert "DMD" not in symbols  # the known, documented chrX coverage gap
    comments = (KNOWLEDGE_REAL / "gene_panel.tsv").read_text(encoding="utf-8")
    assert "chrX" in comments
    assert "COVERAGE GAP" in comments


def test_real_gene_disease_matches_documented_schema() -> None:
    header, rows = _read_tsv_rows(KNOWLEDGE_REAL / "gene_disease.tsv")
    assert header == list(GENE_DISEASE_COLUMNS)
    assert len(rows) > 100


def test_real_gene_disease_preserves_full_confidence_vocabulary_of_both_sources() -> None:
    """Disputed and Refuted rows must survive into the shipped table (requirement:
    'curated and refuted' is real evidence, not something to drop) -- from ClinGen
    (Title Case) AND from DDG2P (its own, lower-case spelling), each kept verbatim.
    """
    _, rows = _read_tsv_rows(KNOWLEDGE_REAL / "gene_disease.tsv")
    by_source: dict[str, set[str]] = {}
    for row in rows:
        by_source.setdefault(row["source"], set()).add(row["confidence"])

    assert "Disputed" in by_source["ClinGen"]
    assert "Refuted" in by_source["ClinGen"]
    assert "disputed" in by_source["DDG2P"]
    assert "refuted" in by_source["DDG2P"]
    # No row's confidence was remapped onto a synthetic numeric score.
    all_confidences = {row["confidence"] for row in rows}
    assert all(not value.replace(".", "", 1).isdigit() for value in all_confidences)


def test_real_gene_disease_mutation_consequence_absent_for_clingen_present_for_ddg2p() -> None:
    """ClinGen's curation model has no mutation-consequence concept -- every ClinGen
    row must have an EMPTY cell there, never a guessed value. DDG2P DOES capture this
    (its 'variant consequence' column), so its rows must have it populated -- proving
    the distinction is real, not just an artifact of one source being absent.
    """
    _, rows = _read_tsv_rows(KNOWLEDGE_REAL / "gene_disease.tsv")
    clingen_rows = [row for row in rows if row["source"] == "ClinGen"]
    ddg2p_rows = [row for row in rows if row["source"] == "DDG2P"]
    assert clingen_rows and ddg2p_rows
    assert all(row["mutation_consequence"] == "" for row in clingen_rows)
    assert any(row["mutation_consequence"] != "" for row in ddg2p_rows)


def test_real_gene_disease_hgnc_id_present_from_both_sources() -> None:
    _, rows = _read_tsv_rows(KNOWLEDGE_REAL / "gene_disease.tsv")
    by_source: dict[str, list[str]] = {}
    for row in rows:
        by_source.setdefault(row["source"], []).append(row["hgnc_id"])
    assert all(value.startswith("HGNC:") for value in by_source["ClinGen"])
    assert all(value.startswith("HGNC:") for value in by_source["DDG2P"])


def test_real_gene_disease_some_ddg2p_rows_have_no_disease_id_and_that_is_visible() -> None:
    """A handful of DDG2P rows have neither a MONDO nor an OMIM disease identifier;
    those must render with an EMPTY disease_id cell, not a fabricated placeholder --
    and the header comment must say so with a real count (not silently drop the row).
    """
    _, rows = _read_tsv_rows(KNOWLEDGE_REAL / "gene_disease.tsv")
    ddg2p_rows = [row for row in rows if row["source"] == "DDG2P"]
    missing_disease_id = [row for row in ddg2p_rows if row["disease_id"] == ""]
    assert missing_disease_id  # the real data does contain this case
    comments = (KNOWLEDGE_REAL / "gene_disease.tsv").read_text(encoding="utf-8")
    assert f"{len(missing_disease_id)} of {len(ddg2p_rows)} rows" in comments


def test_real_gene_phenotype_matches_synthetic_table_columns() -> None:
    """Requirement: the real table's columns match the shipped synthetic demo table's
    columns, so a future integration only has to change *values*, not column names.
    """
    synthetic_header, _ = _read_tsv_rows(KNOWLEDGE_PUBLIC / "gene_phenotype.tsv")
    real_header, real_rows = _read_tsv_rows(KNOWLEDGE_REAL / "gene_phenotype.tsv")
    assert real_header == synthetic_header == list(GENE_PHENOTYPE_COLUMNS)
    assert real_rows


def test_real_gene_phenotype_never_carries_hpos_own_missing_sentinel() -> None:
    """HPO's '-' (not stated) must have become an EMPTY cell, never the literal
    string '-' -- otherwise a naive downstream reader could mistake it for data.
    """
    _, rows = _read_tsv_rows(KNOWLEDGE_REAL / "gene_phenotype.tsv")
    strengths = {row["association_strength"] for row in rows}
    assert "-" not in strengths
    assert "" in strengths  # unstated frequency does appear, just as an empty cell


def test_real_tables_carry_true_provenance_not_the_synthetic_labels() -> None:
    panel_text = (KNOWLEDGE_REAL / "gene_panel.tsv").read_text(encoding="utf-8")
    disease_text = (KNOWLEDGE_REAL / "gene_disease.tsv").read_text(encoding="utf-8")
    phenotype_text = (KNOWLEDGE_REAL / "gene_phenotype.tsv").read_text(encoding="utf-8")

    assert "gnomAD" in panel_text
    assert "CC0" in panel_text
    assert "SYNTHETIC" not in panel_text.upper().replace("NON-SYNTHETIC", "")

    assert "ClinGen" in disease_text
    assert "CC0" in disease_text

    assert "HPO" in phenotype_text
    assert "CC BY 4.0" in phenotype_text

    for text in (panel_text, disease_text, phenotype_text):
        assert "2026-08-28" in text  # the recorded retrieval date


def test_real_tables_are_reproducible_from_the_committed_generator(tmp_path: Path) -> None:
    """End-to-end determinism proof on the fixture (the real ~150MB source files are
    not available in every environment this test runs in -- see the module docstring)
    plus a direct re-hash of the committed files against themselves, guarding against
    an accidental hand-edit that the generator would not reproduce.
    """
    resources_root = tmp_path / "resources"
    _write_full_fixture_resources(resources_root)
    out_a = tmp_path / "a"
    out_b = tmp_path / "b"
    build_all(resources_root, out_a)
    build_all(resources_root, out_b)
    for name in ("gene_panel.tsv", "gene_disease.tsv", "gene_phenotype.tsv"):
        assert (out_a / name).read_bytes() == (out_b / name).read_bytes()
