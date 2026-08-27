"""Clinical assertions from the real NCBI ClinVar GRCh38 release VCF, read locally.

This is the first *real* annotation source in the repository. It replaces
``NullClinicalAdapter`` for the clinical slot (TD-01) and therefore declares
``synthetic = False`` — a deliberate opt-out of the GP-20 mock disclosure, made
only because this genuinely parses a hash-pinned NCBI release rather than a
fabricated table.

Why a local file and tabix, permanently (PRIV-05):

* A proband coordinate is not metadata. Handing ``chr15:40200239 A>G`` to a
  remote variant endpoint discloses patient genetic data to a third party,
  irreversibly and without consent. ``mva.annotation`` is structurally forbidden
  from importing a network client, and this module reaches the ClinVar release
  the same way ``local_tables`` reaches its TSVs: through a file that a separate,
  public-only acquisition step already downloaded and hashed.
* A pinned local file gives byte-identical repeat runs (GP-30). A live API
  cannot.

Four ClinVar-specific traps this module exists to not fall into:

1. **Absence is not benignity (GP-14).** A variant with no ClinVar record is
   *missing* from the returned mapping — never present with an empty tuple, and
   never given a fabricated "not pathogenic" call. The same rule applies one
   level down: a record that carries only an *oncogenicity* (``ONC``) or somatic
   (``SCI``) classification and no germline ``CLNSIG`` has no germline
   classification on record, so it yields no assertion at all. 308 of the 1425
   records in the committed test fixture are exactly that shape.
2. **Review status is not significance.** ``CLNSIG=Pathogenic`` with
   ``CLNREVSTAT=no_assertion_criteria_provided`` (0 stars, one unreviewed
   submitter) and ``CLNSIG=Pathogenic`` with
   ``CLNREVSTAT=reviewed_by_expert_panel`` (3 stars) are the same string and
   wildly different evidence. The star rating is carried through to
   ``ClinicalAssertion.star_rating``, which is what
   ``annotation.service._significance_strength`` grades the evidence on. Both
   cases are in the fixture, at ``chr15:40200239`` and ``chr7:117509033``.
3. **Conflict is a finding, not a tie to break.**
   ``Conflicting_classifications_of_pathogenicity`` is passed through verbatim as
   the significance, so ``_significance_direction`` reads it as NEUTRAL. This
   adapter never resolves a conflict to the most severe or the most common member
   call. (What it cannot yet carry is the ``CLNSIGCONF`` breakdown — see
   :data:`UNREPRESENTABLE_CLINVAR_FIELDS`.)
4. **Contig naming.** ClinVar's VCF uses bare Ensembl-style contigs (``15``,
   ``X``, ``MT``); this pipeline's join key is UCSC-style (``chr15``, ``chrX``,
   ``chrM``). Comparing those raw is the single most dangerous detail called out
   in ``docs/references/track1-submission-contract.md``: it fails by silently
   finding nothing, which is indistinguishable from "ClinVar has no record".
   The mapping is therefore resolved once, explicitly, against the contig names
   the tabix index actually holds, and is exposed as :attr:`contig_map` so a test
   can assert it rather than infer it from a miss.
5. **Allele representation.** ClinVar writes ``100 AT>AG`` where this pipeline's
   normalised proband record holds the minimal ``101 T>G``. Same base change, same
   locus, different string — so the pathogenic assertion vanishes and the variant
   is scored as having no ClinVar record, which is one of the strongest promoting
   signals the ranker has. This adapter used to build its key from the raw VCF
   columns and so could not join those two. It now canonicalises **both** the
   caller's query and every record it reads through
   :func:`mva.alleles.canonicalise_allele` — the same function
   :mod:`mva.ingestion.normalise` uses, not a second copy of the rule, because two
   implementations of one representation rule is precisely what caused the miss.
"""

from __future__ import annotations

import hashlib
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Final, Protocol, cast

from mva.alleles import (
    CanonicalAllele,
    LeftAlignmentStatus,
    ReferenceLookup,
    canonicalise_allele,
    rightmost_equivalent_position,
)
from mva.determinism import hash_file
from mva.errors import AdapterUnavailableError, GenomeBuildMismatchError
from mva.models.genome import CANONICAL_CONTIGS, GenomeBuild, normalise_contig
from mva.models.variant import ClinicalAssertion

# --------------------------------------------------------------------------- identity

#: Adapter identity, stamped on every EvidenceItem this adapter justifies. Names
#: the *format* rather than the database, because the version string already
#: names the release and the two together must be unambiguous in a report footer.
CLINVAR_ADAPTER_NAME: Final = "clinvar-vcf"

#: ``ClinicalAssertion.source``. ClinVar's own name for itself, taken from the
#: ``##source`` header line rather than hardcoded as a display string.
CLINVAR_SOURCE: Final = "ClinVar"

#: Header lines this adapter refuses to run without. ``##fileDate`` is the weekly
#: release date and is the only version ClinVar's VCF carries; ``##reference`` is
#: what makes the coordinates joinable (GP-11).
_HEADER_SOURCE: Final = "##source="
_HEADER_FILE_DATE: Final = "##fileDate="
_HEADER_REFERENCE: Final = "##reference="

# --------------------------------------------------------------------------- semantics

#: ClinVar's review status -> star rating, keyed on the RAW ``CLNREVSTAT`` token
#: exactly as it appears in the VCF (underscores and all), so that the lookup
#: cannot be broken by a change to how the text is prettified for display.
#:
#: This table is the whole of trap 2 above. It is public so a reviewer can check
#: it against ClinVar's published definition without reading the parser, and so a
#: test can assert every value the current release actually contains.
#:
#: The pre-2024 spellings are kept because a run may legitimately be pinned to an
#: older release; ClinVar renamed "interpretation" to "classification" in its
#: 2024 schema change and both forms are in the wild.
CLINVAR_STAR_RATINGS: Final[Mapping[str, int]] = MappingProxyType(
    {
        "practice_guideline": 4,
        "reviewed_by_expert_panel": 3,
        "criteria_provided,_multiple_submitters,_no_conflicts": 2,
        "criteria_provided,_single_submitter": 1,
        "criteria_provided,_conflicting_classifications": 1,
        "criteria_provided,_conflicting_interpretations": 1,  # pre-2024 spelling
        "no_assertion_criteria_provided": 0,
        "no_classification_provided": 0,
        "no_assertion_provided": 0,  # pre-2024 spelling
        "no_classification_for_the_single_variant": 0,
        "no_interpretation_for_the_single_variant": 0,  # pre-2024 spelling
        "no_classifications_from_unflagged_records": 0,
        "no_assertion_for_the_individual_variant": 0,  # pre-2024 spelling
    }
)

#: INFO fields this adapter parses and then has to drop, because
#: ``mva.models.variant.ClinicalAssertion`` has no field that can hold them
#: without lying about what they are. Recorded here rather than left implicit:
#: an undocumented drop is indistinguishable from a parser that never looked.
#:
#: The consequential one is ``CLNSIGCONF``. The *fact* of conflict does survive,
#: in ``significance``; the breakdown ("Likely_pathogenic(1)|Uncertain_significance(3)")
#: does not, and a 1-vs-3 split reads very differently from a 3-vs-1 split. It is
#: deliberately NOT smuggled into ``conditions`` (those are diseases) or emitted
#: as separate assertions (downstream would then see a standalone "Likely
#: pathogenic" call and a standalone "Uncertain significance" call at one site and
#: file SUPPORTS and NEUTRAL evidence for each, which is the conflict being
#: silently resolved by a different route).
UNREPRESENTABLE_CLINVAR_FIELDS: Final[tuple[str, ...]] = (
    "CLNSIGCONF",
    "ALLELEID",
    "CLNHGVS",
    "CLNVI",
    "MC",
    "GENEINFO",
    "ONC",
    "SCI",
)

#: Percent escapes are legal in VCF INFO values and ClinVar uses them (``%3D``
#: and ``%3B`` both appear in the 2026-08-22 release). Decoded with a hand-rolled
#: scanner because ``urllib.parse.unquote`` is unavailable here by design: the
#: whole ``urllib`` package is a forbidden import on the patient-data path
#: (PRIV-05), and importing it for a string helper would trip the architecture
#: test for no benefit.
_HEX_DIGITS: Final[frozenset[str]] = frozenset("0123456789abcdefABCDEF")

#: ClinVar's placeholder for "no value here" in a ``|``-delimited INFO list.
_MISSING_VALUE: Final = "."

#: Positions within this many bases are fetched in one tabix query rather than
#: one each. Exome candidate sets cluster, so this roughly halves the wall time
#: for a realistic batch; a small window keeps the number of ClinVar records read
#: but discarded bounded even in a mutation hotspot.
DEFAULT_MERGE_WINDOW_BP: Final = 1000


class _TabixHandle(Protocol):
    """Exactly the pysam surface this module uses.

    pysam ships no ``py.typed`` marker, so everything it returns is Unknown to
    pyright. Narrowing it to a Protocol at the single point of construction means
    the rest of the module is fully typed against a contract that is written down
    here, instead of against whatever the C extension happens to return.
    """

    @property
    def contigs(self) -> list[str]: ...

    @property
    def header(self) -> Iterator[str]: ...

    def fetch(self, reference: str, start: int, end: int) -> Iterator[str]: ...

    def close(self) -> None: ...


@dataclass(frozen=True, slots=True)
class _Query:
    """One caller-supplied variant ID, decomposed into a joinable coordinate.

    ``position``/``ref``/``alt`` are the **canonical** form, not the raw text of the
    ID: they have been through :func:`mva.alleles.canonicalise_allele`, so they are
    the same representation the records read out of the release are reduced to.
    Both sides of the join therefore agree by construction rather than by luck.
    """

    variant_id: str
    """The caller's original string, which is what the result is keyed by."""

    build: GenomeBuild
    """Carried rather than assumed. A hardcoded ``GRCh38`` here would make the
    query key and the record key disagree the moment the adapter is pointed at a
    GRCh37 release — and disagree silently, returning nothing for every variant."""

    contig: str
    """UCSC-style, e.g. ``chr15``."""

    position: int
    ref: str
    alt: str

    search_end: int
    """Right-most POS at which an equivalent spelling of this variant could sit.

    Equal to ``position`` unless a reference is configured and the variant sits in
    a repeat. An index query finds records by the span they occupy, so a source
    record that spells this same insertion further right occupies a disjoint span
    and is never even fetched — a miss indistinguishable from "ClinVar has no
    record". See :func:`mva.alleles.rightmost_equivalent_position`."""

    @property
    def canonical_key(self) -> str:
        """The build-qualified join key, rebuilt from the normalised parts."""
        return f"{self.build.value}:{self.contig}:{self.position}:{self.ref}:{self.alt}"


# --------------------------------------------------------------------------- integrity


def read_shipped_md5(md5_path: Path) -> str:
    """Parse NCBI's ``clinvar.vcf.gz.md5`` sidecar into a bare hex digest.

    The sidecar's second field is the *remote* build path on NCBI's filesystem,
    which is why only the first whitespace-delimited token is taken.

    A digest that travelled alongside the file it describes proves the download
    was not truncated or corrupted. It does not prove the file is the one anyone
    reviewed — an attacker or a mistake that replaces the VCF can replace the
    sidecar too. Prefer a sha256 recorded in the repository's own resource
    manifest; this exists so a run can be pinned before that manifest exists.
    """
    if not md5_path.is_file():
        msg = (
            f"ClinVar md5 sidecar {md5_path.as_posix()!r} not found. Pass "
            "expected_sha256 instead, or re-run the acquisition step."
        )
        raise AdapterUnavailableError(msg)
    first_line = md5_path.read_text(encoding="utf-8").strip().split("\n", 1)[0]
    digest = first_line.split()[0] if first_line.split() else ""
    if len(digest) != 32 or any(char not in _HEX_DIGITS for char in digest):
        msg = (
            f"ClinVar md5 sidecar {md5_path.name} does not begin with a 32-character hex "
            "digest. The file is not echoed (PRIV-09)."
        )
        raise AdapterUnavailableError(msg)
    return digest.lower()


def _md5_file(path: Path) -> str:
    """Streamed md5, to check a shipped sidecar. Never used as a security control."""
    digest = hashlib.new("md5", usedforsecurity=False)
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def _verify_integrity(path: Path, *, expected_sha256: str | None, expected_md5: str | None) -> None:
    """Check the release's bytes before a single record is read as data.

    Mirrors ``local_tables.verify_manifest``: hash first, open second. An adapter
    built over unverified bytes is annotating against an unknown file, and the
    resulting report would carry a release string it cannot substantiate.

    Fails closed when neither hash is supplied — an unpinned 185 MB download is
    exactly the artifact that quietly becomes a different file between runs.
    Messages name the file and the digests; they never echo its contents
    (PRIV-09), and a hash of public reference data is not patient-derived.
    """
    if expected_sha256 is None and expected_md5 is None:
        msg = (
            f"Refusing to open ClinVar release {path.name!r} without an integrity pin. "
            "Pass expected_sha256 (preferred; recorded in the resource manifest) or "
            "expected_md5 (read from NCBI's shipped .md5 sidecar via read_shipped_md5). "
            "An unpinned resource makes every claim derived from it unreproducible."
        )
        raise AdapterUnavailableError(msg)
    if expected_sha256 is not None:
        actual = hash_file(path)
        if actual != expected_sha256.strip().lower():
            msg = (
                f"ClinVar release {path.name!r} failed its sha256 integrity check: expected "
                f"{expected_sha256.strip().lower()}, found {actual}. Refusing to annotate "
                "against an unpinned resource; re-run the acquisition step and review the "
                "manifest diff."
            )
            raise AdapterUnavailableError(msg)
    if expected_md5 is not None:
        actual_md5 = _md5_file(path)
        if actual_md5 != expected_md5.strip().lower():
            msg = (
                f"ClinVar release {path.name!r} failed the md5 check shipped with it: "
                f"expected {expected_md5.strip().lower()}, found {actual_md5}. The download "
                "is truncated or corrupted; refusing to annotate against it."
            )
            raise AdapterUnavailableError(msg)


# --------------------------------------------------------------------------- parsing


def _percent_decode(value: str) -> str:
    """Decode VCF ``%XX`` escapes in a single left-to-right pass.

    One pass, not a sequence of ``str.replace`` calls: replacing ``%25`` at any
    point other than last would let ``%2525`` decode twice and turn an escaped
    literal into a delimiter.
    """
    if "%" not in value:
        return value
    out: list[str] = []
    index = 0
    length = len(value)
    while index < length:
        char = value[index]
        if char == "%" and index + 2 < length:
            pair = value[index + 1 : index + 3]
            if pair[0] in _HEX_DIGITS and pair[1] in _HEX_DIGITS:
                out.append(chr(int(pair, 16)))
                index += 3
                continue
        out.append(char)
        index += 1
    return "".join(out)


def _decode_text(value: str) -> str:
    """Restore ClinVar's display text from its VCF-safe encoding.

    ClinVar writes spaces as ``_`` because a VCF INFO value may not contain one,
    so ``Uncertain_significance`` *is* "Uncertain significance" and
    ``criteria_provided,_single_submitter`` *is* "criteria provided, single
    submitter". Decoding is restoring the source's own string, not reformatting
    it. Percent escapes are undone first so that a ``%5F`` cannot become a space.
    """
    return _percent_decode(value).replace("_", " ").strip()


def _parse_info(raw: str) -> dict[str, str]:
    """Split a VCF INFO column into keys and raw (still-encoded) values.

    Flag keys (no ``=``) map to the empty string. Values are left encoded so that
    each consumer decodes at the granularity it needs — splitting on ``|`` must
    happen before ``%7C`` is turned back into a literal bar.
    """
    fields: dict[str, str] = {}
    for part in raw.split(";"):
        if not part:
            continue
        key, separator, value = part.partition("=")
        fields[key] = value if separator else ""
    return fields


def _split_encoded_list(value: str) -> tuple[str, ...]:
    """Decode a ``|``-delimited ClinVar INFO list, dropping ``.`` placeholders."""
    items: list[str] = []
    for token in value.split("|"):
        if not token or token == _MISSING_VALUE:
            continue
        decoded = _decode_text(token)
        if decoded and decoded != _MISSING_VALUE:
            items.append(decoded)
    return tuple(items)


def _star_rating(review_status_raw: str | None) -> int | None:
    """Stars for a raw ``CLNREVSTAT`` token, or ``None`` if the token is unknown.

    ``None`` rather than ``0``. ClinVar adds review-status strings over time, and
    an unrecognised one means "we could not grade this review depth" — which is
    absence of information, not a zero-star review (GP-14). ``None`` grades the
    downstream evidence WEAK, which is the conservative reading; a fabricated 0
    would be indistinguishable from a genuinely unreviewed submission.
    """
    if review_status_raw is None:
        return None
    return CLINVAR_STAR_RATINGS.get(review_status_raw)


def _assertion_sort_key(assertion: ClinicalAssertion) -> tuple[str, str, str]:
    """Total order over the assertions attached to one variant (GP-30).

    ``(accession, significance, review_status)``. The accession — ClinVar's
    Variation ID, prefixed — is unique per record and so decides the order on its
    own in every real case; the other two are tie-breakers that keep the order
    total even for a release that somehow repeats one. Lexicographic on the
    accession string, which is stable but not numeric: ``VariationID:100`` sorts
    before ``VariationID:99``. Stability is the requirement; readability of the
    ordering is not.
    """
    return (assertion.accession or "", assertion.significance, assertion.review_status or "")


def _unique_ids(variant_ids: Sequence[str]) -> tuple[str, ...]:
    """Deduplicate while preserving caller order, so the result order is stable."""
    seen: dict[str, None] = {}
    for variant_id in variant_ids:
        seen.setdefault(variant_id, None)
    return tuple(seen)


def merge_query_regions(positions: Sequence[int], window: int) -> tuple[tuple[int, int], ...]:
    """Collapse sorted positions into ``(start, end)`` 1-based inclusive spans.

    Purely a query-count optimisation: every fetched record is still matched on
    the exact canonical key, so a wider window can only cause records to be read
    and discarded, never to be joined to the wrong variant.
    """
    return merge_query_spans([(position, position) for position in positions], window)


def merge_query_spans(spans: Sequence[tuple[int, int]], window: int) -> tuple[tuple[int, int], ...]:
    """Collapse 1-based inclusive spans, merging any pair closer than ``window``.

    The span form exists because a query is no longer a point. With a reference
    configured, a variant in a repeat must be searched for from its left-most to
    its right-most legal position, since a source record anywhere in that tract
    canonicalises to the same key.
    """
    if not spans:
        return ()
    ordered = sorted((min(a, b), max(a, b)) for a, b in spans)
    merged: list[tuple[int, int]] = []
    low, high = ordered[0]
    for start, end in ordered[1:]:
        if start - high <= window:
            high = max(high, end)
        else:
            merged.append((low, high))
            low, high = start, end
    merged.append((low, high))
    return tuple(merged)


def _resolve_contig_map(index_contigs: Sequence[str]) -> dict[str, str]:
    """Map each canonical UCSC contig to the name the tabix index actually uses.

    Resolved from the index rather than assumed, and resolved once, because the
    failure mode of getting it wrong is silence: ``fetch("chr15", ...)`` against
    an Ensembl-named index raises, and a defensively-swallowed exception would
    report "ClinVar has nothing on chr15" for an entire chromosome.

    Mitochondria carry three spellings in the wild (``MT``, ``M``, ``chrM``) and
    ClinVar uses ``MT``; the candidate order below is preference order, not a
    guess. Contigs absent from the index are absent from the map, and variants on
    them are simply not looked up — ClinVar genuinely holds nothing there.
    """
    available = set(index_contigs)
    resolved: dict[str, str] = {}
    for ucsc in CANONICAL_CONTIGS:
        bare = ucsc.removeprefix("chr")
        candidates = ("MT", "M", "chrM", "chrMT") if ucsc == "chrM" else (bare, ucsc)
        for candidate in candidates:
            if candidate in available:
                resolved[ucsc] = candidate
                break
    return resolved


# --------------------------------------------------------------------------- adapter


class ClinvarVcfAdapter:
    """Curated clinical significance from a local, hash-pinned ClinVar release VCF.

    Constructed against a file the acquisition step downloaded; never fetches
    anything itself. The tabix index is opened once and reused for the life of
    the adapter, and lookups are region queries, so a 5,000-variant batch costs
    thousands of index seeks rather than one 185 MB scan per variant.
    """

    def __init__(
        self,
        vcf_path: Path,
        *,
        expected_sha256: str | None = None,
        expected_md5: str | None = None,
        build: GenomeBuild = GenomeBuild.GRCH38,
        merge_window_bp: int = DEFAULT_MERGE_WINDOW_BP,
        reference: ReferenceLookup | None = None,
    ) -> None:
        """Open a verified ClinVar release for region queries.

        Args:
            vcf_path: The bgzipped ClinVar VCF. Its ``.tbi`` must sit beside it.
            expected_sha256: Pin from the repository's resource manifest. Preferred.
            expected_md5: Pin from NCBI's shipped ``.md5`` sidecar, via
                :func:`read_shipped_md5`. Detects a corrupt download, not a
                substituted one. At least one of the two pins is required.
            build: The assembly the caller's variant IDs are in. The release's
                ``##reference`` header must agree, or construction fails: joining
                GRCh37 coordinates to a GRCh38 release mis-locates every variant
                by megabases while looking entirely successful (GP-11).
            merge_window_bp: Query-coalescing window; see :func:`merge_query_regions`.
            reference: optional 1-based inclusive reference accessor. Trimming to a
                minimal representation never needs one and is always applied, so the
                common non-minimal mismatch joins with or without it. Left-alignment
                does need one: without it, an indel this release spells at a
                different position inside the same repeat cannot be reconciled, and
                :attr:`representation_status` says so rather than leaving the miss to
                read as "ClinVar has no record" (GP-14).
        """
        if merge_window_bp < 0:
            msg = "merge_window_bp must not be negative."
            raise AdapterUnavailableError(msg)
        if not vcf_path.is_file():
            msg = (
                f"ClinVar release {vcf_path.as_posix()!r} not found. This adapter reads a "
                "pre-downloaded, hash-pinned local file; it never fetches anything (PRIV-05)."
            )
            raise AdapterUnavailableError(msg)
        index_path = vcf_path.with_name(vcf_path.name + ".tbi")
        if not index_path.is_file():
            msg = (
                f"ClinVar tabix index {index_path.as_posix()!r} not found. Without it every "
                "lookup would become a full scan of the release; refusing to run that way."
            )
            raise AdapterUnavailableError(msg)

        _verify_integrity(vcf_path, expected_sha256=expected_sha256, expected_md5=expected_md5)

        self._vcf_path = vcf_path
        self._index_path = index_path
        self._build = build
        self._merge_window_bp = merge_window_bp
        self._reference = reference
        self._tabix = _open_tabix(vcf_path)
        self._version = _release_version(self._tabix, path=vcf_path, build=build)
        self._contig_map = _resolve_contig_map(self._tabix.contigs)

    # ------------------------------------------------------------------ identity

    @property
    def name(self) -> str:
        return CLINVAR_ADAPTER_NAME

    @property
    def version(self) -> str:
        """The release, read out of the file: ``<##reference>-<##fileDate>``.

        e.g. ``GRCh38-2026-08-22``. Both halves are needed and neither is
        invented: the date identifies which weekly ClinVar release produced the
        classification (they change), and the assembly identifies which of the
        two coordinate systems it is stated in. This string is stamped onto every
        Citation, and ``EvidenceItem`` refuses a DATABASE_ASSERTION without it.
        """
        return self._version

    @property
    def synthetic(self) -> bool:
        """False, deliberately (GP-20).

        ``is_synthetic`` fails closed, so every adapter is treated as a mock until
        it says otherwise. This one says otherwise because it parses an actual
        NCBI ClinVar release whose bytes were verified against a pinned hash
        before the first record was read. It is *not* a claim that the underlying
        assertions are correct — ClinVar submissions conflict, get reclassified
        and carry review depths from "one unreviewed submitter" to "practice
        guideline", which is exactly why the star rating is preserved.
        """
        return False

    @property
    def vcf_path(self) -> Path:
        return self._vcf_path

    @property
    def index_path(self) -> Path:
        return self._index_path

    @property
    def build(self) -> GenomeBuild:
        return self._build

    @property
    def representation_status(self) -> LeftAlignmentStatus:
        """Whether this adapter can reconcile a *shifted* indel spelling, typed.

        ``APPLIED`` when a reference was supplied, ``UNAVAILABLE_NO_REFERENCE``
        otherwise. Trimming to the minimal representation is unconditional either
        way, so the ``100 AT>AG`` versus ``101 T>G`` class of mismatch joins in both
        states; what the degraded state costs is the repeat-tract class, where the
        release and the proband place the same insertion at different positions.
        Typed rather than logged so a report can state the limitation instead of a
        reader having to infer it from an absent assertion.
        """
        if self._reference is None:
            return LeftAlignmentStatus.UNAVAILABLE_NO_REFERENCE
        return LeftAlignmentStatus.APPLIED

    @property
    def representation_limitation(self) -> str | None:
        """One sentence for a report footer, or ``None`` when nothing is degraded."""
        if self._reference is not None:
            return None
        return (
            "ClinVar lookups were canonicalised by trimming only: no reference FASTA "
            "was supplied to the adapter, so an indel that ClinVar places at a "
            "different position within the same repeat tract could not be matched. "
            "Such a variant is reported as having no ClinVar record, which is absence "
            "of information and not evidence of benignity."
        )

    @property
    def contig_map(self) -> Mapping[str, str]:
        """UCSC contig -> the name this release's index uses, e.g. ``chr15 -> 15``.

        Exposed so the prefix mapping is assertable directly. A test that only
        checks lookups would pass against a broken map for any variant ClinVar
        happens not to hold.
        """
        return MappingProxyType(self._contig_map)

    # ------------------------------------------------------------------ lookup

    def assertions(self, variant_ids: Sequence[str]) -> Mapping[str, tuple[ClinicalAssertion, ...]]:
        """Return ClinVar's germline classifications for the variants it holds.

        Keys appear **only** for variants with at least one germline
        classification on record. A missing key means "ClinVar has nothing to say
        here", which is not evidence of benignity and must not be rendered as any
        (GP-14). A variant whose only ClinVar record carries an oncogenicity or
        somatic classification is missing for the same reason.

        Raises:
            GenomeBuildMismatchError: a variant ID is in a different build from
                the release. Refused rather than skipped — a silent miss would
                read as "ClinVar has no record" for every variant in the batch.
            ValueError: a variant ID is not the canonical
                ``{build}:{contig}:{pos}:{ref}:{alt}`` form, or names a contig
                this pipeline does not reason about.
        """
        ordered = _unique_ids(variant_ids)
        queries = tuple(self._parse_query(variant_id) for variant_id in ordered)

        # Canonical key -> every caller spelling that reduces to it, so the result is
        # keyed by exactly what was passed in even where canonicalisation changed the
        # text. A list, not a single ID: two spellings of one variant in the same
        # batch now share a key, and mapping the key to only the last of them would
        # answer one caller and silently drop the other.
        by_key: dict[str, list[str]] = {}
        for query in queries:
            by_key.setdefault(query.canonical_key, []).append(query.variant_id)

        by_contig: dict[str, list[tuple[int, int]]] = {}
        for query in queries:
            by_contig.setdefault(query.contig, []).append((query.position, query.search_end))

        found: dict[str, list[ClinicalAssertion]] = {}
        # tabix returns every record that *overlaps* a region, and the regions no
        # longer partition the queried positions, so one release record can be
        # handed back by two fetches. Identity is (join key, ClinVar record), which
        # keeps two genuinely distinct submissions at one allele and collapses the
        # same submission seen twice. The previous guard — drop any record whose raw
        # POS falls outside the span — cannot be used any more: a record at POS 100
        # spelled ``AT>AG`` canonicalises into the span queried for POS 101, and
        # discarding it on its raw position is the very miss being fixed.
        seen: set[tuple[str, tuple[str, str, str]]] = set()
        for ucsc_contig in sorted(by_contig):
            index_contig = self._contig_map.get(ucsc_contig)
            if index_contig is None:
                # This release's index holds nothing for that chromosome at all.
                continue
            for span in merge_query_spans(by_contig[ucsc_contig], self._merge_window_bp):
                for line in self._tabix.fetch(index_contig, span[0] - 1, span[1]):
                    for key, assertion in self._parse_record(line, ucsc_contig):
                        for variant_id in by_key.get(key, ()):
                            identity = (variant_id, _assertion_sort_key(assertion))
                            if identity in seen:
                                continue
                            seen.add(identity)
                            found.setdefault(variant_id, []).append(assertion)

        return {
            variant_id: tuple(sorted(found[variant_id], key=_assertion_sort_key))
            for variant_id in ordered
            if variant_id in found
        }

    def close(self) -> None:
        """Release the htslib file handle. Idempotent."""
        self._tabix.close()

    # ------------------------------------------------------------------ internals

    def canonicalise(self, contig: str, position: int, ref: str, alt: str) -> CanonicalAllele:
        """The single entry point both sides of the join go through.

        Delegates to :func:`mva.alleles.canonicalise_allele` — the same function
        :mod:`mva.ingestion.normalise` calls. There is deliberately no trimming or
        shifting logic in this module: a second implementation of the rule is what
        made ``100 AT>AG`` and ``101 T>G`` fail to join in the first place.

        Public so that a test can compare this adapter's representation against the
        ingestion stage's *directly*, rather than inferring that they agree from a
        join that happened to succeed. Agreement inferred from a passing join is
        exactly the evidence that was available while the two disagreed.
        """
        return canonicalise_allele(
            contig=contig,
            position=position,
            ref=ref,
            alt=alt,
            reference=self._reference,
        )

    def _parse_query(self, variant_id: str) -> _Query:
        """Decompose a canonical variant ID, refusing anything ambiguous.

        The caller's ID is canonicalised here rather than trusted. Ingestion already
        normalises proband records, so in the pipeline this is a no-op — but the
        adapter is also called with hand-written and third-party IDs, and an
        adapter that only joins correctly when its caller happened to normalise
        first is an adapter whose correctness lives somewhere else.

        Error text carries the field count and the build names only. The
        coordinate itself is patient data and is never echoed (PRIV-09).
        """
        parts = variant_id.split(":")
        if len(parts) != 5:
            msg = (
                f"Variant ID has {len(parts)} colon-separated fields, expected 5 "
                "({build}:{contig}:{pos}:{ref}:{alt}). The value is not echoed (PRIV-09)."
            )
            raise ValueError(msg)
        build_token, contig_token, position_token, ref, alt = parts
        build = GenomeBuild.parse(build_token)
        if build is not self._build:
            msg = (
                f"A {build.value} variant ID was passed to a {self._build.value} ClinVar "
                "adapter. Cross-build lookup is refused: the same locus differs by megabases "
                "between assemblies, so the join would silently return nothing or, worse, the "
                "wrong record. Lift-over must be an explicit, provenance-tracked stage. The "
                "coordinate is not echoed (PRIV-09)."
            )
            raise GenomeBuildMismatchError(msg)
        try:
            position = int(position_token)
        except ValueError as exc:
            msg = "Variant ID position is not an integer. The value is not echoed (PRIV-09)."
            raise ValueError(msg) from exc
        contig = normalise_contig(contig_token)
        canonical = self.canonicalise(contig, position, ref.strip().upper(), alt.strip().upper())
        search_end = canonical.position
        if self._reference is not None:
            search_end = max(
                search_end,
                rightmost_equivalent_position(
                    contig=contig,
                    position=canonical.position,
                    ref=canonical.ref,
                    alt=canonical.alt,
                    reference=self._reference,
                ),
            )
        return _Query(
            variant_id=variant_id,
            build=build,
            contig=contig,
            position=canonical.position,
            ref=canonical.ref,
            alt=canonical.alt,
            search_end=search_end,
        )

    def _parse_record(self, line: str, ucsc_contig: str) -> Iterator[tuple[str, ClinicalAssertion]]:
        """Turn one ClinVar VCF line into ``(canonical key, assertion)`` pairs.

        One pair per ALT allele. ClinVar writes one ALT per line in practice — a
        full scan of the 2026-08-22 release found 0 multi-ALT records in
        4,467,990 — but VCF permits several and a reader that assumed otherwise
        would attach a pathogenic call to the first allele and drop the rest.

        On a multi-ALT record every split allele receives the *same* assertion,
        because ClinVar's INFO fields are declared ``Number=1`` or ``Number=.``,
        never ``Number=A``: the file simply does not say which ALT a
        classification belongs to. That is a limitation of the source recorded
        faithfully, not a judgement this adapter is making.

        The key is built from the **canonicalised** allele, never from the raw VCF
        columns. ClinVar's release contains genuinely non-minimal spellings — a
        pathogenic ``100 AT>AG`` is the same substitution as the proband's
        ``101 T>G`` — and keying on the raw columns loses the assertion entirely
        while looking exactly like "ClinVar has no record here" (GP-14).

        De-duplication of a record returned by two overlapping fetches is the
        caller's job, in :meth:`assertions`, on the join key rather than on the raw
        position: after canonicalisation a record's POS is no longer the thing that
        decides which query it answers.
        """
        columns = line.split("\t", 7)
        if len(columns) < 8:
            return
        _, position_text, record_id, ref, alts, _qual, _filter, info_text = columns
        position = int(position_text)
        info = _parse_info(info_text)

        significance_raw = info.get("CLNSIG")
        if not significance_raw:
            # Germline classification absent. The record may still carry ONC
            # (oncogenicity) or SCI (somatic clinical impact), which are separate
            # classification axes that ClinicalAssertion cannot represent and
            # which are emphatically not germline evidence. Nothing on record.
            return

        review_status_raw = info.get("CLNREVSTAT") or None
        assertion = ClinicalAssertion(
            source=CLINVAR_SOURCE,
            version=self._version,
            accession=_accession(record_id, info.get("ALLELEID")),
            # Verbatim structure: ClinVar's own ``|`` between multiple aggregate
            # classifications is preserved, so "Pathogenic|drug_response" stays
            # visibly two calls rather than collapsing into one.
            significance="|".join(_decode_text(token) for token in significance_raw.split("|")),
            review_status=_decode_text(review_status_raw) if review_status_raw else None,
            star_rating=_star_rating(review_status_raw),
            conditions=_split_encoded_list(info.get("CLNDN", "")),
        )

        ref_allele = ref.strip().upper()
        for alt in alts.split(","):
            allele = alt.strip().upper()
            if not allele or allele == _MISSING_VALUE:
                continue
            canonical = self.canonicalise(ucsc_contig, position, ref_allele, allele)
            key = (
                f"{self._build.value}:{ucsc_contig}:"
                f"{canonical.position}:{canonical.ref}:{canonical.alt}"
            )
            yield key, assertion


def _accession(record_id: str, allele_id: str | None) -> str | None:
    """ClinVar's identifier for the record, prefixed so it is self-describing.

    The VCF ID column is the ClinVar Variation ID (``##ID=<Description="ClinVar
    Variation ID">``), which is what resolves at
    ``ncbi.nlm.nih.gov/clinvar/variation/<id>``. A bare integer in a Citation
    would be unresolvable by a reader, so the kind is named. The VCV accession is
    deliberately not synthesised from the Variation ID: it is derivable, but
    printing an accession the file never contained is inventing provenance.
    """
    if record_id and record_id != _MISSING_VALUE:
        return f"VariationID:{record_id.strip()}"
    if allele_id:
        return f"AlleleID:{allele_id.strip()}"
    return None


def _open_tabix(path: Path) -> _TabixHandle:
    """Open the release for region queries, or explain why we cannot.

    pysam is imported here rather than at module scope so that importing this
    module — which the composition root and the architecture tests both do —
    does not require the optional ``genomics`` extra to be installed.
    """
    try:
        import pysam  # noqa: PLC0415 - native backend, imported on demand
    except ImportError as exc:  # pragma: no cover - guarded by the genomics extra
        msg = (
            "The pysam backend is not installed; install the 'genomics' extra. Tabix region "
            "queries are how this adapter avoids scanning a 185 MB release once per variant."
        )
        raise AdapterUnavailableError(msg) from exc
    try:
        return cast(_TabixHandle, pysam.TabixFile(str(path)))
    except (OSError, ValueError) as exc:
        msg = (
            f"htslib could not open ClinVar release {path.name!r} with its tabix index; the "
            "file may be truncated, not bgzipped, or indexed against different bytes."
        )
        raise AdapterUnavailableError(msg) from exc


def _release_version(tabix: _TabixHandle, *, path: Path, build: GenomeBuild) -> str:
    """Read the release identity out of the VCF header, and check it is joinable.

    Three header lines matter and all three are required. ``##source`` proves the
    file is ClinVar and not some other VCF that happens to be bgzipped;
    ``##reference`` proves the coordinates are in the build the caller's variant
    IDs use; ``##fileDate`` is the only release version ClinVar's VCF carries,
    and without it every assertion would be uncitable.
    """
    source: str | None = None
    file_date: str | None = None
    reference: str | None = None
    for line in tabix.header:
        if line.startswith(_HEADER_SOURCE):
            source = line[len(_HEADER_SOURCE) :].strip()
        elif line.startswith(_HEADER_FILE_DATE):
            file_date = line[len(_HEADER_FILE_DATE) :].strip()
        elif line.startswith(_HEADER_REFERENCE):
            reference = line[len(_HEADER_REFERENCE) :].strip()
        elif not line.startswith("##"):
            break

    if source is None or source.lower() != CLINVAR_SOURCE.lower():
        msg = (
            f"{path.name!r} does not declare '##source=ClinVar'. This adapter parses "
            "ClinVar's INFO vocabulary (CLNSIG/CLNREVSTAT/CLNDN); pointing it at another "
            "VCF would produce assertions attributed to ClinVar that ClinVar never made."
        )
        raise AdapterUnavailableError(msg)
    if file_date is None:
        msg = (
            f"{path.name!r} carries no '##fileDate' header, so the release cannot be "
            "identified. Every ClinVar claim is only as current as its release; an "
            "unversioned assertion is unreproducible and EvidenceItem refuses one."
        )
        raise AdapterUnavailableError(msg)
    if reference is None:
        msg = (
            f"{path.name!r} carries no '##reference' header, so the assembly its "
            "coordinates are stated in is unknown. A coordinate without a build is "
            "invalid (GP-11)."
        )
        raise AdapterUnavailableError(msg)
    try:
        declared = GenomeBuild.parse(reference)
    except ValueError as exc:
        msg = (
            f"{path.name!r} declares an unrecognised '##reference={reference}'. Refusing to "
            "join coordinates whose assembly cannot be established (GP-11)."
        )
        raise AdapterUnavailableError(msg) from exc
    if declared is not build:
        msg = (
            f"{path.name!r} is a {declared.value} release but this adapter was constructed "
            f"for {build.value}. The same locus differs by megabases between assemblies, so "
            "the join would return wrong records rather than no records (GP-11)."
        )
        raise AdapterUnavailableError(msg)

    return f"{declared.value}-{file_date}"
