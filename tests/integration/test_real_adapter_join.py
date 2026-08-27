"""Do the real files actually join? Measured against the real files, not fixtures.

Every representation test elsewhere in the suite runs on constructed records, so
each one proves the *code* is self-consistent. None of them proves the thing that
matters: that the proband's spelling of a variant and the spelling in the NCBI
ClinVar release resolve to the same string. That question can only be answered by
the actual release and the actual reference genome, so it is answered here.

Two sources are needed and neither ships with the repository:

* the ClinVar GRCh38 release (~185 MB, hash-pinned at open time), and
* the GRCh38 no-alt analysis-set FASTA (~3.1 GB decompressed, ``chr``-prefixed).

They are located under a resources root outside the repo (public reference data,
never patient data — GP-40 is about the workspace, which nothing here touches).
When one is missing the tests **skip with the reason**, because a
representational failure that has never been checked against real data must not
be able to hide behind a green suite that quietly checked nothing.
"""

from __future__ import annotations

import os
from collections.abc import Iterator
from pathlib import Path

import pytest

from mva.alleles import LeftAlignmentStatus, canonicalise_allele, rightmost_equivalent_position
from mva.annotation.clinvar_vcf import ClinvarVcfAdapter
from mva.config import find_repo_root
from mva.determinism import hash_file, stable_hash
from mva.errors import AdapterUnavailableError
from mva.ingestion.normalise import FastaReference, normalise_variants, open_reference_fasta
from mva.models.genome import GenomeBuild, GenomicCoordinate
from mva.models.variant import OP_LEFT_ALIGN, FilterStatus, Genotype, VariantRecord, Zygosity

pytestmark = [pytest.mark.integration, pytest.mark.slow]

REPO_ROOT = find_repo_root(Path(__file__))
RESOURCES = Path(os.environ.get("MVA_RESOURCES") or REPO_ROOT.parent / "mva-resources")

#: GRCh38 chr1 length. A FASTA still being written has a short or absent chr1, so
#: this is the completeness check: an index built over a partial download would
#: answer every fetch with silence or the wrong base, and left-alignment would then
#: be confidently wrong rather than merely absent.
CHR1_LENGTH_GRCH38 = 248_956_422

#: BUB1B. The MVA gene, so the region this project actually reasons about is the
#: region the join is measured on.
BUB1B_REGION = ("chr15", 40_160_000, 40_230_000)


def _clinvar_release() -> Path | None:
    path = RESOURCES / "clinvar" / "clinvar.vcf.gz"
    return path if path.is_file() and path.with_name(path.name + ".tbi").is_file() else None


def _reference_fasta() -> Path | None:
    """The first candidate FASTA that exists and is demonstrably complete."""
    override = os.environ.get("MVA_REFERENCE_FASTA")
    candidates = [Path(override)] if override else []
    candidates += [
        RESOURCES / "reference" / "GRCh38_no_alt.fa",
        RESOURCES / "reference" / "GRCh38_full_analysis_set_plus_decoy_hla.fa",
    ]
    for path in candidates:
        if not path.is_file():
            continue
        try:
            reference = open_reference_fasta(path)
        except AdapterUnavailableError:
            continue  # still downloading, or not indexable yet
        try:
            if len(reference.fetch("chr1", CHR1_LENGTH_GRCH38, CHR1_LENGTH_GRCH38)) == 1:
                return path
        except (KeyError, ValueError, IndexError, OSError):
            continue
        finally:
            reference.close()
    return None


@pytest.fixture(scope="module")
def clinvar() -> Iterator[ClinvarVcfAdapter]:
    path = _clinvar_release()
    if path is None:
        pytest.skip(f"ClinVar release not present under {RESOURCES / 'clinvar'}")
    adapter = ClinvarVcfAdapter(path, expected_sha256=hash_file(path))
    yield adapter
    adapter.close()


@pytest.fixture(scope="module")
def reference() -> Iterator[FastaReference]:
    path = _reference_fasta()
    if path is None:
        pytest.skip(
            "No complete GRCh38 reference FASTA found. Left-alignment cannot be "
            f"measured against real data. Looked under {RESOURCES / 'reference'}; "
            "set MVA_REFERENCE_FASTA to override."
        )
    handle = open_reference_fasta(path)
    yield handle
    handle.close()


def make_record(contig: str, position: int, ref: str, alt: str) -> VariantRecord:
    return VariantRecord(
        coordinate=GenomicCoordinate(
            build=GenomeBuild.GRCH38, contig=contig, position=position, ref=ref, alt=alt
        ),
        genotype=Genotype(
            zygosity=Zygosity.HET,
            genotype_string="0/1",
            depth=40,
            ref_reads=20,
            alt_reads=20,
            genotype_quality=99,
        ),
        filter_status=FilterStatus.PASS,
        source_artifact="tests/integration/test_real_adapter_join.py",
    )


def _release_records(adapter: ClinvarVcfAdapter) -> list[list[str]]:
    """Raw columns for every ClinVar record over BUB1B, read independently.

    Opened separately rather than through the adapter's own handle: this is the
    control side of the comparison, and it must not share the code path it exists
    to check.
    """
    import pysam

    contig, start, end = BUB1B_REGION
    index_contig = adapter.contig_map[contig]
    handle = pysam.TabixFile(str(adapter.vcf_path))
    try:
        return [line.split("\t", 8) for line in handle.fetch(index_contig, start - 1, end)]
    finally:
        handle.close()


# ---------------------------------------------------------------------------
# What the real release actually looks like
# ---------------------------------------------------------------------------


def test_the_real_release_is_already_minimally_represented(clinvar: ClinvarVcfAdapter) -> None:
    """Measured, not assumed — and it decides how much each half of the fix buys.

    A full pass over the pinned release found 0 non-minimal ALT entries in
    4,467,990 records, and this re-checks the MVA region of whatever release is
    installed. So the trimming half of the fix currently buys robustness rather
    than recovered assertions: nothing in *this* release needs it. The
    left-alignment half is the opposite — the release is left-aligned and a
    right-shifted proband indel misses every time — which is why the reference
    FASTA is not optional for a real run.
    """
    non_minimal = [
        (row[1], row[3], alt)
        for row in _release_records(clinvar)
        for alt in row[4].split(",")
        if alt not in {".", "*", ""}
        and canonicalise_allele(contig="chr15", position=int(row[1]), ref=row[3], alt=alt).changed
    ]
    assert non_minimal == [], "the release is no longer minimal; the trimming path now matters"


def test_a_real_pathogenic_record_joins_through_the_canonical_key(
    clinvar: ClinvarVcfAdapter,
) -> None:
    """The whole point, on real records: keys built by the shared rule find them.

    Every germline-classified record in the region is looked up under the key the
    shared canonicalisation produces from its own columns. All of them must come
    back. A single miss is a variant that would be reported as having no ClinVar
    record — which is not evidence of benignity, and reads downstream as novelty.
    """
    rows = _release_records(clinvar)
    germline = [row for row in rows if "CLNSIG=" in row[7]]
    assert germline, "the installed release has no germline classifications over BUB1B"

    # ALT "." means the record states no alternate allele, so it yields no
    # assertion at all. Excluded here because it is the source saying nothing, not
    # a join that failed — the distinction GP-14 exists to keep.
    keys = [
        f"GRCh38:chr15:{canonical.position}:{canonical.ref}:{canonical.alt}"
        for row in germline
        for alt in row[4].split(",")
        if alt not in {".", "*", ""}
        for canonical in [
            canonicalise_allele(contig="chr15", position=int(row[1]), ref=row[3], alt=alt)
        ]
    ]
    found = clinvar.assertions(keys)
    assert set(found) == set(keys), "a real ClinVar record did not answer to its own key"


def test_real_lookups_are_byte_identical_across_runs(clinvar: ClinvarVcfAdapter) -> None:
    """GP-30 on the real release, over the new canonicalising lookup path."""
    rows = _release_records(clinvar)[:200]
    keys = [f"GRCh38:chr15:{row[1]}:{row[3]}:{row[4].split(',')[0]}" for row in rows]

    def digest() -> str:
        result = clinvar.assertions(keys)
        return stable_hash(
            {key: [item.model_dump(mode="json") for item in value] for key, value in result.items()}
        )

    assert digest() == digest()


# ---------------------------------------------------------------------------
# The join that needs the reference genome
# ---------------------------------------------------------------------------


def _homopolymer_insertion(
    clinvar: ClinvarVcfAdapter, reference: FastaReference
) -> tuple[list[str], str] | None:
    """A real ClinVar single-base insertion sitting in a real repeat tract.

    Returns every legal spelling of it, plus the canonical key, or ``None`` when
    the region holds no such record. Spellings are derived from the reference
    rather than invented: for a one-base insertion of ``b`` anchored at ``p``, every
    position from ``p`` to the end of the run of ``b`` that follows is an equally
    valid VCF representation of the identical event.
    """
    for row in _release_records(clinvar):
        position, ref, alts = int(row[1]), row[3], row[4]
        for alt in alts.split(","):
            if len(ref) != 1 or len(alt) != 2 or alt[0] != ref:
                continue
            canonical = canonicalise_allele(
                contig="chr15", position=position, ref=ref, alt=alt, reference=reference
            )
            rightmost = rightmost_equivalent_position(
                contig="chr15",
                position=canonical.position,
                ref=canonical.ref,
                alt=canonical.alt,
                reference=reference,
            )
            if rightmost <= canonical.position + 1:
                continue  # not in a repeat: nothing to shift, nothing to prove
            inserted = alt[1]
            spellings = [f"{canonical.position}:{canonical.ref}:{canonical.alt}"]
            spellings.extend(
                f"{shifted}:{reference.fetch('chr15', shifted, shifted)}:"
                f"{reference.fetch('chr15', shifted, shifted)}{inserted}"
                for shifted in range(canonical.position + 1, rightmost)
                if reference.fetch("chr15", shifted, shifted) == inserted
            )
            if len(spellings) < 2:
                continue
            key = f"GRCh38:chr15:{canonical.position}:{canonical.ref}:{canonical.alt}"
            return [f"GRCh38:chr15:{item}" for item in spellings], key
    return None


def test_a_right_shifted_proband_indel_joins_the_real_release(
    clinvar: ClinvarVcfAdapter, reference: FastaReference
) -> None:
    """The finding, measured end to end on real data.

    A proband VCF that spells an insertion at the right-hand end of a repeat has a
    join key ClinVar does not hold. Normalising against the real GRCh38 FASTA moves
    it to the left-most position — the one ClinVar and gnomAD store — and the
    assertion appears. Without the FASTA it does not, and the variant is scored as
    novel: this test asserts both halves, because the second is the measurement of
    what the missing reference costs.
    """
    found = _homopolymer_insertion(clinvar, reference)
    if found is None:
        pytest.skip("no ClinVar single-base insertion in a repeat tract over BUB1B")
    spellings, canonical_key = found

    for spelling in spellings:
        _, contig, position, ref, alt = spelling.split(":")
        record = make_record(contig, int(position), ref, alt)

        aligned = normalise_variants([record], reference=reference)
        assert aligned.left_alignment.status is LeftAlignmentStatus.APPLIED
        assert aligned.variants[0].variant_id == canonical_key, (
            "a legal spelling of a real ClinVar insertion did not reach its key"
        )
        assert clinvar.assertions([aligned.variants[0].variant_id]), (
            "the left-aligned key found no ClinVar record"
        )

        if spelling == canonical_key:
            continue
        # The counterfactual: the same record with no reference available.
        degraded = normalise_variants([record])
        assert degraded.left_alignment.status is LeftAlignmentStatus.UNAVAILABLE_NO_REFERENCE
        assert degraded.variants[0].variant_id != canonical_key
        assert clinvar.assertions([degraded.variants[0].variant_id]) == {}, (
            "the un-aligned key must miss — that miss is what the reference buys back"
        )


def test_the_reference_confirms_the_ref_alleles_of_the_real_release(
    clinvar: ClinvarVcfAdapter, reference: FastaReference
) -> None:
    """Guards the premise of every left-alignment: the FASTA is the right assembly.

    If the FASTA and the release disagree about REF, the alignment is being
    computed against a different genome and every shifted coordinate is wrong while
    looking entirely successful. ClinVar's REF alleles are the cheapest independent
    check of that available.
    """
    rows = _release_records(clinvar)[:500]
    assert rows, "the installed release has no records over BUB1B"
    mismatches = [
        row[1]
        for row in rows
        if reference.fetch("chr15", int(row[1]), int(row[1]) + len(row[3]) - 1) != row[3].upper()
    ]
    assert mismatches == [], "reference FASTA and ClinVar release disagree on REF"


def test_left_alignment_moves_real_coordinates_and_records_that_it_did(
    clinvar: ClinvarVcfAdapter, reference: FastaReference
) -> None:
    """GP-31 on real coordinates: a POS that moved is a POS with a provenance entry.

    The Track 1 submission is scored on exact coordinates. A silently moved
    coordinate is both unauditable and, if the move was wrong, unrecoverable.
    """
    found = _homopolymer_insertion(clinvar, reference)
    if found is None:
        pytest.skip("no ClinVar single-base insertion in a repeat tract over BUB1B")
    spellings, canonical_key = found
    shifted = [spelling for spelling in spellings if spelling != canonical_key]

    for spelling in shifted:
        _, contig, position, ref, alt = spelling.split(":")
        result = normalise_variants(
            [make_record(contig, int(position), ref, alt)], reference=reference
        )
        assert OP_LEFT_ALIGN in result.variants[0].normalisation_ops
        assert result.operations_applied[OP_LEFT_ALIGN] == 1
