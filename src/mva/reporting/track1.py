"""The Track 1 submission CSV — the artifact the challenge actually scores.

The contract is vendored at ``docs/references/track1-submission-contract.md`` and
was derived by reading the challenge Space's source. Everything here implements
that file; where the two disagree, the contract wins and this module is the bug.

Three details are worth more than the rest of this module combined:

1. **Chromosome strings are compared raw.** The scorer builds its variant key as
   ``(chrom.strip(), int(pos), ref.strip().upper(), alt.strip().upper())`` — there
   is no chr-prefix normalisation on either side. Emitting ``15`` where the answer
   key holds ``chr15`` scores exactly zero while looking perfectly correct to a
   human reader. Ensembl-style GRCh38 VCFs use bare contig names, so this is a
   live hazard. Every row is forced through
   :func:`~mva.models.genome.normalise_contig` and
   :func:`validate_submission` re-checks the prefix on the rendered bytes.
2. **A compound-heterozygous pair is ONE row**, using the ``_2`` columns. Two rows
   would be read as two independent single-variant proposals and would burn two of
   the ten slots.
3. **``epcr`` must satisfy ``0 < epcr <= 1``.** Zero is rejected outright, so the
   floor of :func:`composite_to_epcr` is positive by construction: a candidate we
   scored 0.0 is a weak prediction, not an absent one, and it still deserves its
   rank-ordered slot.

**Privacy.** This file is a PUBLIC artifact that carries genomic coordinates by
necessity — that is the challenge's own format, and the coordinates are the
prediction. Everything *else* is discretionary, so ``notes`` is constructed from
structured, gene-level and mechanism-class values only. Free text from the
pipeline (``rank_rationale``, clinical narrative, phenotype terms) never reaches
it; see :func:`_safe_note`.
"""

from __future__ import annotations

import csv
import importlib
import io
import logging
import math
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from mva.models import ContigStyle, Sensitivity, normalise_contig

if TYPE_CHECKING:  # pragma: no cover - typing only
    from mva.models import CandidatePair

_LOG = logging.getLogger(__name__)

#: The exact 12 columns, in the exact order, from the verified contract.
TRACK1_COLUMNS: tuple[str, ...] = (
    "proband_id",
    "chrom_1",
    "pos_1",
    "ref_1",
    "alt_1",
    "chrom_2",
    "pos_2",
    "ref_2",
    "alt_2",
    "epcr",
    "finding_type",
    "notes",
)

#: Hard limit. ``if len(rows) > 10: raise ValueError`` in the scorer.
MAX_SUBMISSION_ROWS = 10

#: The only proband ID the scorer accepts; anything else hard-fails the whole
#: submission with "Unknown proband_id".
ACCEPTED_PROBAND_ID = "PROBAND01"

#: EPCR bounds. The floor is strictly positive because the scorer validates
#: ``0 < epcr <= 1`` and rejects a zero outright.
EPCR_FLOOR = 0.01
EPCR_CEILING = 1.0

#: Decimals in the rendered EPCR. Four is enough to keep every distinct composite
#: score distinguishable at this scale while staying readable; ties that survive
#: rounding fall back to file order, which is why rows are pre-sorted.
EPCR_DECIMALS = 4

#: Accepted values for ``finding_type``. Blank is also legal (the scorer reads the
#: column with ``.get()``), but we always write one for legibility.
FINDING_TYPES: frozenset[str] = frozenset({"primary", "secondary"})

#: Flag on a :class:`~mva.models.CandidatePair` that marks it as an incidental
#: finding unrelated to the primary phenotype.
INCIDENTAL_FLAG = "incidental_finding"

#: ``notes`` is never read by the scorer, so its only job is to be safe and
#: legible. Anything longer is a paragraph, and paragraphs are where clinical
#: detail leaks in.
NOTE_MAX_CHARS = 120

#: Characters permitted in ``notes``. A whitelist, not a blacklist: the field is
#: assembled from structured values, so anything outside this set means something
#: unexpected reached it.
_NOTE_ALLOWED = re.compile(r"[^A-Za-z0-9 ;=,._+/-]")
_DIGITS_ONLY = re.compile(r"^[0-9]+$")
_ALLELE_RE = re.compile(r"^[ACGTN]+$")
_SECOND_VARIANT_COLUMNS = ("chrom_2", "pos_2", "ref_2", "alt_2")


def composite_to_epcr(composite: float) -> float:
    """Map a composite score on ``[0, 1]`` to an EPCR on ``[0.01, 1.0]``.

    ``epcr = 0.01 + 0.99 * clamp(composite, 0, 1)``, rounded to four decimals.

    Why affine and why a floor:

    * The scorer validates ``0 < epcr <= 1`` and **rejects zero**. A composite of
      0.0 is a weak candidate, not a withdrawn one, so it must still map to a
      positive value or the row is discarded and its rank slot wasted.
    * The map is strictly increasing, so it preserves the ranking exactly. Rank
      points are derived from EPCR order alone, and F-max is swept over EPCR
      thresholds; any monotone map gives identical rank points, and an affine one
      keeps the threshold sweep evenly spread rather than bunching candidates.

    **This is a rank-ordering confidence, not a calibrated probability.** Nothing
    in this pipeline is fitted to labelled outcomes (ASSUMPTION-SCORING-01), so
    ``epcr = 0.60`` does not mean "60% chance this pair is causal". It means "we
    rank this above everything with a lower number". The reports say so in words.
    """
    if math.isnan(composite):
        msg = "composite score is NaN; a candidate without a score cannot be ranked."
        raise ValueError(msg)
    clamped = min(max(composite, 0.0), 1.0)
    epcr = EPCR_FLOOR + clamped * (EPCR_CEILING - EPCR_FLOOR)
    return min(max(round(epcr, EPCR_DECIMALS), EPCR_FLOOR), EPCR_CEILING)


@dataclass(frozen=True)
class SubmissionRow:
    """One CSV row: either a variant pair (``_2`` populated) or a single variant.

    ``pos_2`` is a string rather than an ``int | None`` because the column is
    genuinely empty for single-variant proposals and the CSV is the source of
    truth for what "empty" looks like. ``pos_1`` is an ``int`` because it is never
    optional.
    """

    proband_id: str
    chrom_1: str
    pos_1: int
    ref_1: str
    alt_1: str
    chrom_2: str
    pos_2: str
    ref_2: str
    alt_2: str
    epcr: float
    finding_type: str = "primary"
    notes: str = ""

    @property
    def is_pair(self) -> bool:
        return bool(self.chrom_2)

    def as_mapping(self) -> dict[str, str]:
        """The row as CSV-ready strings, keyed by column name."""
        return {
            "proband_id": self.proband_id,
            "chrom_1": self.chrom_1,
            "pos_1": str(self.pos_1),
            "ref_1": self.ref_1,
            "alt_1": self.alt_1,
            "chrom_2": self.chrom_2,
            "pos_2": self.pos_2,
            "ref_2": self.ref_2,
            "alt_2": self.alt_2,
            "epcr": f"{self.epcr:.{EPCR_DECIMALS}f}",
            "finding_type": self.finding_type,
            "notes": self.notes,
        }

    def variant_keys(self) -> frozenset[tuple[str, int, str, str]]:
        """The scorer's own key form, used here to detect duplicate rows."""
        keys = {(self.chrom_1.strip(), self.pos_1, self.ref_1.upper(), self.alt_1.upper())}
        if self.is_pair:
            keys.add(
                (
                    self.chrom_2.strip(),
                    int(self.pos_2),
                    self.ref_2.upper(),
                    self.alt_2.upper(),
                )
            )
        return frozenset(keys)


def _safe_note(pair: CandidatePair) -> str:
    """Build ``notes`` from structured, gene-level values only.

    Deliberately assembled rather than copied: gene symbol, inheritance model and
    phase status are enum values and identifiers, none of which is clinical detail
    about the child. ``rank_rationale`` and anything phenotype-derived is excluded
    by construction, not by filtering — the free-text fields are simply never read
    here. The result is then whitelisted and truncated as a second line of defence
    against an unexpected value arriving in a gene symbol.
    """
    parts = [
        pair.gene_symbol,
        pair.inheritance_model.value,
        f"phase={pair.phase.status.value}",
    ]
    note = _NOTE_ALLOWED.sub("", "; ".join(parts)).strip()
    return note[:NOTE_MAX_CHARS].rstrip("; ")


def _row_for_pair(pair: CandidatePair, proband_id: str) -> SubmissionRow:
    first = pair.variant_a.coordinate
    second = pair.variant_b.coordinate if pair.variant_b is not None else None
    return SubmissionRow(
        proband_id=proband_id,
        # Forced, not trusted: see the module docstring. A bare '15' here is a
        # silent zero score.
        chrom_1=normalise_contig(first.contig, ContigStyle.UCSC),
        pos_1=int(first.position),
        ref_1=first.ref.strip().upper(),
        alt_1=first.alt.strip().upper(),
        chrom_2=normalise_contig(second.contig, ContigStyle.UCSC) if second else "",
        pos_2=str(int(second.position)) if second else "",
        ref_2=second.ref.strip().upper() if second else "",
        alt_2=second.alt.strip().upper() if second else "",
        epcr=composite_to_epcr(pair.composite_score),
        finding_type="secondary" if INCIDENTAL_FLAG in pair.flags else "primary",
        notes=_safe_note(pair),
    )


def truncation_notice(total_candidates: int, *, max_rows: int = MAX_SUBMISSION_ROWS) -> str | None:
    """A one-line record of truncation, or ``None`` if nothing was dropped.

    Returned rather than hidden in a log so the dossier can state it: a reader
    comparing the dossier's ranked list against the submission must be able to see
    why the list is shorter, instead of inferring that the missing candidates were
    rejected on their merits.
    """
    if total_candidates <= max_rows:
        return None
    return (
        f"Submission truncated: {total_candidates} ranked candidates, "
        f"{max_rows} submitted, {total_candidates - max_rows} omitted solely because "
        f"the challenge format accepts at most {max_rows} rows. Omission here is a "
        "format limit, not a scientific rejection; the full ranking is in the dossier."
    )


def build_submission_rows(
    pairs: Sequence[CandidatePair],
    *,
    proband_id: str,
    max_rows: int = MAX_SUBMISSION_ROWS,
) -> tuple[SubmissionRow, ...]:
    """Rank, convert and truncate candidate pairs into submission rows.

    Ordering is by EPCR descending. The scorer derives rank exactly this way
    (``sorted(enumerate(rows), key=lambda x: (-x[1][1], x[0]))``), so ties break by
    file order — which means our file order is a scientific decision, not a
    formatting one. Ties therefore fall back to the candidate sort key (composite
    score, then genomic position), giving a total, reproducible order (GP-30).
    """
    if not 1 <= max_rows <= MAX_SUBMISSION_ROWS:
        msg = (
            f"max_rows={max_rows} is outside 1..{MAX_SUBMISSION_ROWS}. The challenge "
            f"scorer raises ValueError above {MAX_SUBMISSION_ROWS} rows."
        )
        raise ValueError(msg)

    ordered = sorted(pairs, key=lambda pair: pair.sort_key())
    rows = [_row_for_pair(pair, proband_id) for pair in ordered]
    # Stable sort: candidates whose EPCR collides after rounding keep the
    # score-then-position order established above.
    rows.sort(key=lambda row: -row.epcr)

    if len(rows) > max_rows:
        # Counts only, never coordinates: an exception or log line travels
        # everywhere the process's logs travel (GP-41).
        _LOG.warning(
            "Track 1 submission truncated to the format limit: %d candidate rows, "
            "%d submitted, %d omitted.",
            len(rows),
            max_rows,
            len(rows) - max_rows,
        )
    return tuple(rows[:max_rows])


def render_submission_csv(rows: Sequence[SubmissionRow]) -> str:
    """Render rows to the exact CSV the challenge expects, header included."""
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=list(TRACK1_COLUMNS),
        lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    for row in rows:
        writer.writerow(row.as_mapping())
    return buffer.getvalue()


def write_submission(rows: Sequence[SubmissionRow], path: Path) -> Path:
    """Render, self-validate, write, and gate the submission as a public export.

    The self-check runs *before* the bytes hit disk: a malformed submission that
    exists is worse than one that does not, because it will be uploaded. The
    export gate then runs on the written file (GP-43) — classification is a claim,
    the re-scan is the verification — and a blocked file is removed rather than
    left for someone to find and submit.
    """
    text = render_submission_csv(rows)
    ok, errors = validate_submission(text)
    if not ok:
        msg = (
            "Refusing to write a submission that fails its own contract check "
            f"({len(errors)} problem(s)): " + "; ".join(errors)
        )
        raise ValueError(msg)

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")
    _gate_public_export(path)
    return path


def _gate_public_export(path: Path) -> None:
    """Run the public-export gate over the written file, if it is available.

    Bound dynamically because ``mva.privacy.export`` is a peer deliverable that may
    not be present in every checkout; its absence degrades to "ungated", which is
    recorded in the log rather than silently assumed to be fine. When the gate is
    present and blocks, the file is deleted and the error propagates: fail closed.
    """
    try:
        module = importlib.import_module("mva.privacy.export")
    except ImportError:
        _LOG.info(
            "mva.privacy.export is unavailable; submission %s was written without the "
            "public-export gate (GP-43). The composition root must gate it before upload.",
            path.name,
        )
        return
    gate = getattr(module, "gate_public_export", None)
    if gate is None:  # pragma: no cover - defensive against interface drift
        _LOG.info("mva.privacy.export exposes no gate_public_export; submission left ungated.")
        return
    try:
        gate(path, declared=Sensitivity.PUBLIC, allowlist=(path.name,))
    except Exception:
        path.unlink(missing_ok=True)
        raise


def _validate_row(index: int, row: dict[str, str | None], errors: list[str]) -> None:
    """Apply every per-row rule from the contract. ``index`` is 1-based."""
    where = f"row {index}"
    missing = [column for column in TRACK1_COLUMNS if row.get(column) is None]
    if missing:
        errors.append(f"{where}: missing value(s) for {missing}; every column must be present.")
        return
    values = {column: (row[column] or "").strip() for column in TRACK1_COLUMNS}

    if values["proband_id"] != ACCEPTED_PROBAND_ID:
        errors.append(
            f"{where}: proband_id {values['proband_id']!r} is not {ACCEPTED_PROBAND_ID!r}; "
            "the scorer hard-fails the whole submission with 'Unknown proband_id'."
        )

    _validate_variant(where, "1", values, errors, required=True)
    populated = [column for column in _SECOND_VARIANT_COLUMNS if values[column]]
    if populated and len(populated) != len(_SECOND_VARIANT_COLUMNS):
        errors.append(
            f"{where}: second variant is partially populated ({sorted(populated)}). It must "
            "be entirely blank (single-variant proposal) or entirely present (one pair row)."
        )
    elif populated:
        _validate_variant(where, "2", values, errors, required=False)
        if (values["chrom_1"], values["pos_1"], values["ref_1"], values["alt_1"]) == (
            values["chrom_2"],
            values["pos_2"],
            values["ref_2"],
            values["alt_2"],
        ):
            errors.append(
                f"{where}: both variants are identical; a variant cannot pair with itself."
            )

    _validate_epcr(where, values["epcr"], errors)

    if values["finding_type"] and values["finding_type"] not in FINDING_TYPES:
        errors.append(
            f"{where}: finding_type {values['finding_type']!r} must be one of "
            f"{sorted(FINDING_TYPES)} (or blank)."
        )

    note = row["notes"] or ""
    if len(note) > NOTE_MAX_CHARS:
        errors.append(
            f"{where}: notes is {len(note)} characters (limit {NOTE_MAX_CHARS}). The "
            "submission is a public artifact; notes carries mechanism-class and "
            "gene-level statements only, never clinical narrative."
        )
    if _NOTE_ALLOWED.search(note):
        errors.append(
            f"{where}: notes contains characters outside the permitted set. It is "
            "assembled from structured values; free text must not reach it."
        )


def _validate_contig(where: str, suffix: str, chrom: str, errors: list[str]) -> None:
    """The one check that silently costs every point if it is wrong."""
    if not chrom.startswith("chr"):
        errors.append(
            f"{where}: chrom_{suffix}={chrom!r} lacks the UCSC 'chr' prefix. The scorer "
            "compares chromosome strings raw, so this row would silently score zero."
        )
        return
    try:
        canonical = normalise_contig(chrom, ContigStyle.UCSC)
    except ValueError:
        errors.append(f"{where}: chrom_{suffix}={chrom!r} is not a canonical contig.")
        return
    if canonical != chrom:
        errors.append(
            f"{where}: chrom_{suffix}={chrom!r} is not canonical; expected {canonical!r}."
        )


def _validate_alleles(where: str, suffix: str, ref: str, alt: str, errors: list[str]) -> None:
    for label, allele in ((f"ref_{suffix}", ref), (f"alt_{suffix}", alt)):
        if not allele:
            errors.append(f"{where}: {label} is empty.")
        elif allele != allele.upper():
            errors.append(f"{where}: {label}={allele!r} must be uppercase.")
        elif not _ALLELE_RE.match(allele):
            errors.append(f"{where}: {label}={allele!r} is not an IUPAC ACGTN allele string.")
    if ref and ref == alt:
        errors.append(f"{where}: ref_{suffix} and alt_{suffix} are identical ({ref!r}).")


def _validate_variant(
    where: str,
    suffix: str,
    values: dict[str, str],
    errors: list[str],
    *,
    required: bool,
) -> None:
    chrom = values[f"chrom_{suffix}"]
    pos = values[f"pos_{suffix}"]

    if not chrom:
        if required:
            errors.append(f"{where}: chrom_{suffix} is empty; the first variant is mandatory.")
        return

    _validate_contig(where, suffix, chrom, errors)

    if not _DIGITS_ONLY.match(pos):
        errors.append(
            f"{where}: pos_{suffix}={pos!r} must be a plain integer — no commas, no "
            "decimal point, no scientific notation. The scorer calls int() on it."
        )
    elif int(pos) <= 0:
        errors.append(f"{where}: pos_{suffix}={pos!r} must be a 1-based position greater than 0.")

    _validate_alleles(where, suffix, values[f"ref_{suffix}"], values[f"alt_{suffix}"], errors)


def _validate_epcr(where: str, raw: str, errors: list[str]) -> None:
    try:
        epcr = float(raw)
    except ValueError:
        errors.append(f"{where}: epcr={raw!r} is not a float; the scorer raises ValueError.")
        return
    if math.isnan(epcr) or not 0.0 < epcr <= 1.0:
        errors.append(
            f"{where}: epcr={raw!r} violates 0 < epcr <= 1. Zero is rejected outright, "
            "which is why composite_to_epcr has a positive floor."
        )


def validate_submission(csv_text: str) -> tuple[bool, tuple[str, ...]]:
    """Re-parse rendered CSV and re-apply every rule in the contract.

    This is a self-check against *our own* reading of the challenge scorer, not a
    substitute for it. It exists because the failure mode it guards against is
    silent: a submission with a bare contig or a zero EPCR uploads cleanly, scores
    nothing, and burns one of six attempts.
    """
    errors: list[str] = []
    reader = csv.DictReader(io.StringIO(csv_text, newline=""))
    fieldnames = tuple(reader.fieldnames or ())
    if fieldnames != TRACK1_COLUMNS:
        errors.append(
            f"header mismatch: got {list(fieldnames)}, expected {list(TRACK1_COLUMNS)}. The "
            "header row is required and the column order is fixed by the contract."
        )
        return (False, tuple(errors))

    rows = list(reader)
    if len(rows) > MAX_SUBMISSION_ROWS:
        errors.append(
            f"{len(rows)} data rows exceed the {MAX_SUBMISSION_ROWS}-row limit; the scorer "
            "raises ValueError rather than truncating."
        )

    probands: list[str] = []
    epcrs: list[float] = []
    seen_keys: dict[frozenset[tuple[str, int, str, str]], int] = {}
    for index, row in enumerate(rows, start=1):
        _validate_row(index, row, errors)
        probands.append((row.get("proband_id") or "").strip())
        try:
            epcrs.append(float((row.get("epcr") or "").strip()))
        except ValueError:
            epcrs.append(math.nan)
        key = _row_key(row)
        if key is not None:
            first_seen = seen_keys.setdefault(key, index)
            if first_seen != index:
                errors.append(
                    f"row {index}: duplicates the variant set already proposed in row "
                    f"{first_seen}; a repeated row wastes one of the ten slots."
                )

    distinct = sorted(set(probands))
    if len(distinct) > 1:
        errors.append(
            f"submission mixes {len(distinct)} proband IDs ({distinct}); the scorer "
            "evaluates only the first proband encountered and hard-fails on a stray ID."
        )

    ordered = [value for value in epcrs if not math.isnan(value)]
    if ordered != sorted(ordered, reverse=True):
        errors.append(
            "rows are not sorted by epcr descending. Rank is derived from epcr order with "
            "ties broken by file order, so ordering is part of the prediction."
        )

    return (not errors, tuple(errors))


def _row_key(row: dict[str, str | None]) -> frozenset[tuple[str, int, str, str]] | None:
    """The scorer's variant-key set for a parsed row, or ``None`` if unparseable."""
    try:
        keys = {
            (
                (row["chrom_1"] or "").strip(),
                int((row["pos_1"] or "").strip()),
                (row["ref_1"] or "").strip().upper(),
                (row["alt_1"] or "").strip().upper(),
            )
        }
        if (row.get("chrom_2") or "").strip():
            keys.add(
                (
                    (row["chrom_2"] or "").strip(),
                    int((row["pos_2"] or "").strip()),
                    (row["ref_2"] or "").strip().upper(),
                    (row["alt_2"] or "").strip().upper(),
                )
            )
    except (KeyError, TypeError, ValueError):
        return None
    return frozenset(keys)


__all__ = [
    "ACCEPTED_PROBAND_ID",
    "EPCR_CEILING",
    "EPCR_FLOOR",
    "MAX_SUBMISSION_ROWS",
    "TRACK1_COLUMNS",
    "SubmissionRow",
    "build_submission_rows",
    "composite_to_epcr",
    "render_submission_csv",
    "truncation_notice",
    "validate_submission",
    "write_submission",
]
