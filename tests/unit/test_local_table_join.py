"""The default annotation adapters must join on the ONE representation rule (ADR 0018).

`mva.alleles.canonicalise_allele` exists because `chr1:100 AT>AG` and `chr1:101 T>G`
are the same substitution written two ways, and every join in this pipeline is a
string comparison on `GenomicCoordinate.variant_id`. When the two sides disagree the
join does not raise. It returns nothing — and "no frequency record" is scored by this
pipeline as *novel and ultra-rare*, the strongest promoting signal the ranker has,
while "no consequence record" can get the allele deleted outright by selection.

The gnomAD and ClinVar adapters were fixed under ADR 0018 and are covered by
`tests/unit/test_gnomad_adapter.py` and `tests/unit/test_clinvar_adapter.py`. The
local TSV adapters were not: they stored the table's `variant_id` column verbatim
and looked it up by direct string membership. That matters more than "they are
synthetic" suggests, because they are the **default executable adapter path** — the
one an organizer reruns, and the one `just demo` uses.

These tests assert the agreement *directly against `canonicalise_allele`*, not by
observing that a join happened to succeed. Agreement inferred from a passing join is
exactly the evidence that was available while the two implementations disagreed.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mva.alleles import LeftAlignmentStatus, canonicalise_allele
from mva.annotation.local_tables import (
    LocalConsequenceAdapter,
    LocalFrequencyAdapter,
    _key_is_indel,
)
from mva.errors import AdapterUnavailableError

pytestmark = pytest.mark.unit

_CONSEQUENCE_HEADER = "\t".join(
    (
        "variant_id",
        "gene_symbol",
        "gene_id",
        "transcript_id",
        "transcript_biotype",
        "is_canonical",
        "is_mane_select",
        "consequence_terms",
        "impact",
        "hgvs_c",
        "hgvs_p",
        "exon",
        "protein_position",
        "amino_acids",
        "splice_ai_delta_max",
        "cadd_phred",
        "revel",
    )
)

_FREQUENCY_HEADER = "\t".join(
    (
        "variant_id",
        "source",
        "version",
        "population",
        "allele_frequency",
        "allele_count",
        "allele_number",
        "homozygote_count",
        "filter_status",
    )
)


def _consequence_row(variant_id: str, *, gene: str = "SYNTHKIN1") -> str:
    return "\t".join(
        (
            variant_id,
            gene,
            "SYNTHG0001",
            f"SYNTHT-{variant_id.replace(':', '-')}.1",
            "protein_coding",
            "1",
            "1",
            "missense_variant",
            "moderate",
            "c.1A>G",
            "p.Met1Val",
            "1/2",
            "1",
            "M/V",
            "0.01",
            "22.0",
            "0.40",
        )
    )


def _frequency_row(variant_id: str, *, allele_frequency: str = "0.25") -> str:
    return "\t".join(
        (
            variant_id,
            "SYNTHETIC-gnomAD-substitute",
            "synthetic-v0.0",
            "global",
            allele_frequency,
            "5",
            "20",
            "0",
            "PASS",
        )
    )


def _consequence_table(tmp_path: Path, *variant_ids: str) -> LocalConsequenceAdapter:
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "consequences.tsv"
    body = "\n".join(_consequence_row(vid) for vid in variant_ids)
    path.write_text(f"# SYNTHETIC fixture\n{_CONSEQUENCE_HEADER}\n{body}\n", encoding="utf-8")
    return LocalConsequenceAdapter(path, version="synthetic-v0.0")


def _frequency_table(tmp_path: Path, *variant_ids: str) -> LocalFrequencyAdapter:
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "frequencies.tsv"
    body = "\n".join(_frequency_row(vid) for vid in variant_ids)
    path.write_text(f"# SYNTHETIC fixture\n{_FREQUENCY_HEADER}\n{body}\n", encoding="utf-8")
    return LocalFrequencyAdapter(path, version="synthetic-v0.0")


#: The ADR 0018 worked example, verbatim. Same substitution, two legal spellings.
NON_MINIMAL = "GRCh38:chr1:100:AT:AG"
MINIMAL = "GRCh38:chr1:101:T:G"


def test_the_two_spellings_really_are_one_variant() -> None:
    """Pin the premise before testing the join, so a failure below is unambiguous."""
    canonical = canonicalise_allele(contig="chr1", position=100, ref="AT", alt="AG")
    assert (canonical.position, canonical.ref, canonical.alt) == (101, "T", "G")
    assert canonical.trimmed
    assert not canonical.left_aligned, "no reference was supplied; nothing may claim otherwise"


# ---------------------------------------------------------------------------
# The join, both directions
# ---------------------------------------------------------------------------


def test_a_non_minimal_table_key_is_found_by_a_minimal_query(tmp_path: Path) -> None:
    """The reproduction. The table spells it `100 AT>AG`; ingestion produced `101 T>G`.

    Before the fix this returned ``{}`` — indistinguishable from "this adapter has
    never heard of that variant", which is what a genuinely novel causal allele
    also looks like.
    """
    adapter = _consequence_table(tmp_path, NON_MINIMAL)
    found = adapter.annotate([MINIMAL])

    assert MINIMAL in found, "the minimal query missed the non-minimal table row"
    assert found[MINIMAL][0].gene_symbol == "SYNTHKIN1"


def test_a_minimal_table_key_is_found_by_a_non_minimal_query(tmp_path: Path) -> None:
    """Both sides go through the rule, so the failure is symmetric and so is the fix."""
    adapter = _consequence_table(tmp_path, MINIMAL)
    found = adapter.annotate([NON_MINIMAL])

    assert NON_MINIMAL in found
    assert found[NON_MINIMAL][0].gene_symbol == "SYNTHKIN1"


def test_the_frequency_adapter_joins_the_same_way(tmp_path: Path) -> None:
    """A missed frequency is the expensive half: absence reads as ultra-rare (GP-14)."""
    adapter = _frequency_table(tmp_path, NON_MINIMAL)
    found = adapter.frequencies([MINIMAL])

    assert MINIMAL in found, (
        "a common variant spelled non-minimally in the table returned no frequency at "
        "all, which this pipeline scores as novel and ultra-rare"
    )
    assert found[MINIMAL][0].allele_frequency == pytest.approx(0.25)


def test_a_bare_contig_in_the_table_still_joins(tmp_path: Path) -> None:
    """`15` and `chr15` are one contig; only the challenge SCORER compares them raw."""
    adapter = _consequence_table(tmp_path, "GRCh38:15:40200000:C:T")
    found = adapter.annotate(["GRCh38:chr15:40200000:C:T"])
    assert "GRCh38:chr15:40200000:C:T" in found


def test_a_different_build_never_joins(tmp_path: Path) -> None:
    """The same locus differs by megabases between assemblies (GP-11)."""
    adapter = _consequence_table(tmp_path, "GRCh37:chr1:101:T:G")
    assert adapter.annotate([MINIMAL]) == {}


# ---------------------------------------------------------------------------
# The contract the callers depend on is unchanged
# ---------------------------------------------------------------------------


def test_results_are_keyed_by_the_caller_s_own_id_not_the_canonical_one(
    tmp_path: Path,
) -> None:
    """`mva.annotation.service` looks the result up by ``record.variant_id``.

    Re-keying the mapping to the canonical form would turn a fixed join into a
    silent miss one layer up, which is the same defect wearing a different hat.
    """
    adapter = _consequence_table(tmp_path, MINIMAL)
    found = adapter.annotate([NON_MINIMAL])

    assert list(found) == [NON_MINIMAL]
    assert MINIMAL not in found


def test_an_unknown_variant_is_omitted_and_not_defaulted(tmp_path: Path) -> None:
    """GP-14 survives the change: absent stays absent, never an empty tuple."""
    adapter = _frequency_table(tmp_path, MINIMAL)
    assert adapter.frequencies(["GRCh38:chr9:1234:A:T"]) == {}


def test_an_unparseable_key_is_matched_verbatim_rather_than_guessed_at(
    tmp_path: Path,
) -> None:
    """An ID that is not ``build:contig:pos:ref:alt`` is not silently reinterpreted.

    Refusing to guess is the point: a coordinate moved on a guess is exactly what
    `canonicalise_allele` documents itself as refusing to do for symbolic alleles.
    """
    adapter = _consequence_table(tmp_path, "not-a-variant-id")
    assert "not-a-variant-id" in adapter.annotate(["not-a-variant-id"])


def test_two_spellings_of_one_variant_in_the_table_are_merged_not_shadowed(
    tmp_path: Path,
) -> None:
    """Both rows describe the same allele, so both annotations belong to it.

    Keying on the raw string kept them as two unrelated entries and a query found
    at most one of them, depending on how the caller happened to spell it.
    """
    adapter = _consequence_table(tmp_path, NON_MINIMAL, MINIMAL)
    found = adapter.annotate([MINIMAL])
    assert len(found[MINIMAL]) == 2


def test_the_index_is_deterministic_under_table_row_order(tmp_path: Path) -> None:  # GP-30
    forward = _consequence_table(tmp_path / "a", NON_MINIMAL, MINIMAL).annotate([MINIMAL])
    backward = _consequence_table(tmp_path / "b", MINIMAL, NON_MINIMAL).annotate([MINIMAL])
    assert forward == backward


# ---------------------------------------------------------------------------
# GP-14 / ADR 0018: what these adapters may NOT claim
# ---------------------------------------------------------------------------


def test_an_indel_query_declares_that_it_was_never_left_aligned(
    tmp_path: Path,
) -> None:
    """Trimming needs no reference; left-alignment does, and a TSV has none.

    The honest output is not "left-alignment skipped" in a log. It is a typed
    statement, carried on the adapter, that every indel join in this run may be
    missing for representational reasons — the same `LeftAlignmentReport` the
    ingestion stage and the real adapters use.

    The report is about the *run*, so it is empty until the run asks something.
    """
    adapter = _frequency_table(tmp_path, "GRCh38:chr1:100:A:AT")
    assert adapter.left_alignment.status is LeftAlignmentStatus.NOT_REQUIRED

    adapter.frequencies(["GRCh38:chr1:100:A:AT"])
    report = adapter.left_alignment

    assert report.status is LeftAlignmentStatus.UNAVAILABLE_NO_REFERENCE
    assert report.reference_available is False
    assert report.indel_count == 1
    assert report.shifted_count == 0
    assert report.is_degraded
    assert "left-align" in report.describe().lower()


def test_a_run_of_snvs_only_says_left_alignment_was_not_required(
    tmp_path: Path,
) -> None:
    """ "Could not" and "did not need to" are opposite claims about how much to trust
    the rarity of every indel in the run, so they must not share a value."""
    adapter = _consequence_table(tmp_path, MINIMAL, "GRCh38:chr2:500:G:C")
    adapter.annotate([MINIMAL, "GRCh38:chr2:500:G:C"])
    report = adapter.left_alignment

    assert report.status is LeftAlignmentStatus.NOT_REQUIRED
    assert not report.is_degraded


# --------------------------------------------------------------------------- the
# count is on the query side of the join, not the table side


def test_an_snv_only_table_does_not_deny_the_runs_indels(tmp_path: Path) -> None:
    """The defect, in the direction that hides the damage.

    `_left_alignment_report` counted `sum(1 for key in index if _key_is_indel(key))`
    — indels in the LOOKUP TABLE. An SNV-only table therefore reported
    `NOT_REQUIRED`, whose `describe()` states "this run contains no indel records",
    over a run whose indels are exactly the ones silently returning absence: a
    trim-only key cannot match a left-aligned source, the miss is read as "no
    frequency data", and that is scored novel and ultra-rare (GP-14).

    The table's composition is not a fact about the run. The query set is.
    """
    adapter = _frequency_table(tmp_path, MINIMAL, "GRCh38:chr2:500:G:C")
    assert not any(_key_is_indel(key) for key in adapter._index), "premise: no indels in the table"

    # A run made of indels, none of which the table holds.
    assert adapter.frequencies(["GRCh38:chr1:100:A:AT", "GRCh38:chr3:700:GT:G"]) == {}

    report = adapter.left_alignment
    assert report.status is LeftAlignmentStatus.UNAVAILABLE_NO_REFERENCE, (
        "an SNV-only table denied that the run contained any indels, over the very "
        "indels that had just come back as absent"
    )
    assert report.indel_count == 2
    assert report.is_degraded
    assert "no indel records" not in report.describe()


def test_a_table_full_of_indels_does_not_degrade_an_snv_only_run(tmp_path: Path) -> None:
    """The mirror. A table's indels say nothing about a run that asked about none.

    Counting the table also over-reported: a run of pure SNVs against an
    indel-heavy table was told its indel joins might be wrong when it had made no
    indel joins at all. A warning a reader cannot act on crowds out the ones they can.
    """
    adapter = _consequence_table(tmp_path, "GRCh38:chr1:100:A:AT", "GRCh38:chr3:700:GT:G")
    adapter.annotate([MINIMAL, "GRCh38:chr2:500:G:C"])

    report = adapter.left_alignment
    assert report.status is LeftAlignmentStatus.NOT_REQUIRED
    assert report.indel_count == 0
    assert not report.is_degraded


def test_a_missed_indel_query_is_counted_exactly_like_a_hit(tmp_path: Path) -> None:
    """A miss is the shape a representation mismatch takes, so it must be counted.

    Counting hits only would reproduce the defect one layer down: the indels that
    fail to join would be invisible to the report that exists to qualify them.
    """
    adapter = _frequency_table(tmp_path, "GRCh38:chr1:100:A:AT")
    adapter.frequencies(["GRCh38:chr1:100:A:AT", "GRCh38:chr9:900:C:CT"])
    assert adapter.left_alignment.indel_count == 2


def test_an_indel_written_at_another_position_in_a_repeat_tract_still_misses(
    tmp_path: Path,
) -> None:
    """The residual limitation, asserted rather than hoped about.

    `chr1:100 A>AT` and `chr1:104 T>TT` are the same insertion into a homopolymer,
    and only a reference FASTA can prove it. Trim-only canonicalisation cannot, so
    the join fails — and the adapter's own report is what says so. Locking this in
    stops a future reader concluding from the tests above that the local adapters
    are representation-complete.
    """
    adapter = _frequency_table(tmp_path, "GRCh38:chr1:100:A:AT")
    assert adapter.frequencies(["GRCh38:chr1:104:T:TT"]) == {}
    assert adapter.left_alignment.is_degraded, (
        "the miss above is invisible unless the adapter declares the degradation that causes it"
    )


# --------------------------------------------------------------------------- an
# empty table is refused, not hash-verified into uselessness


def test_a_header_only_consequence_table_is_refused(tmp_path: Path) -> None:
    """`_read_rows` refused a missing header and accepted zero data rows.

    A header-only `consequences.tsv` therefore built a hash-verified `AdapterSet`
    over an empty index. Every lookup missed, and a miss in these adapters is the
    *correct* answer for a variant the table does not hold -- so nothing anywhere
    could tell the two apart. Every variant lost its `gene_symbols`, and
    compound-heterozygous pairing groups by gene, so the run produced zero
    candidate pairs and reported success.

    Both sibling modules already guarded this: `gene_intervals._read_gtf` refuses a
    gene model with no rows, and `GnomadSitesFrequencyAdapter` refuses a release
    with no complete shard, both for the same reason.
    """
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "consequences.tsv"
    path.write_text(f"# SYNTHETIC fixture\n{_CONSEQUENCE_HEADER}\n", encoding="utf-8")

    with pytest.raises(AdapterUnavailableError) as excinfo:
        LocalConsequenceAdapter(path, version="synthetic-v0.0")
    message = str(excinfo.value)
    assert "no data rows" in message
    assert "consequences.tsv" in message
    # The remediation is in the message: a reader may be an agent whose only view
    # of the rule is this string.
    assert "Re-run the knowledge-table build" in message


def test_a_header_only_frequency_table_is_refused(tmp_path: Path) -> None:
    """The same guard on the frequency slot, which feeds the rarity signal."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "frequencies.tsv"
    path.write_text(f"# SYNTHETIC fixture\n{_FREQUENCY_HEADER}\n", encoding="utf-8")

    with pytest.raises(AdapterUnavailableError, match="no data rows"):
        LocalFrequencyAdapter(path, version="synthetic-v0.0")


def test_a_comments_only_table_is_refused_too(tmp_path: Path) -> None:
    """Comment lines are skipped before the header, so they must not count as rows."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "consequences.tsv"
    path.write_text(
        f"# SYNTHETIC fixture\n{_CONSEQUENCE_HEADER}\n# a trailing note\n\n",
        encoding="utf-8",
    )
    with pytest.raises(AdapterUnavailableError, match="no data rows"):
        LocalConsequenceAdapter(path, version="synthetic-v0.0")


def test_a_single_row_table_still_builds(tmp_path: Path) -> None:
    """The guard is 'no rows', not 'few rows'. One real row is a valid table."""
    adapter = _consequence_table(tmp_path, MINIMAL)
    assert adapter.annotate([MINIMAL])
