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
from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING, cast

from mva.errors import ExportBlockedError
from mva.models import ContigStyle, Sensitivity, normalise_contig

if TYPE_CHECKING:  # pragma: no cover - typing only
    from mva.models import CandidatePair
    from mva.privacy.export import ExportDecision

_LOG = logging.getLogger(__name__)

#: Used only when ``mva.privacy.export`` is present but predates
#: ``PUBLIC_EXPORT_ALLOWLIST``. It is a fixed name, not ``path.name``: an
#: allowlist derived from the file under inspection is not an allowlist.
TRACK1_EXPORT_ALLOWLIST_FALLBACK: tuple[str, ...] = ("track1_submission.csv",)

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

#: Minimum gap between the EPCRs of two adjacent submitted rows.
#:
#: The scorer can see exactly two properties of the EPCR vector: the ORDER of the
#: rows and which rows are TIED. Rank is ``sorted(enumerate(rows), key=(-epcr,
#: file_index))``, and F-max sweeps thresholds over the emitted values, unioning
#: ``predicted_variants |= row.variants`` down to each threshold. Two rows sharing
#: an EPCR are therefore inseparable at every threshold: no sweep can ever predict
#: one without the other. If the true row is one of a tied pair, the partner's
#: variants are false positives at the true row's own threshold — worth up to
#: -0.333 F-max — and if file order puts the wrong row first, the true row is also
#: demoted one rank, which across a tier boundary is a further -50 rank points.
#:
#: One rounding unit (1e-4) would be enough to break a tie: the contract's
#: re-verification against the executed scorer confirms the sweep runs over the
#: unique values WE emit, compared with ``>=``, and not over a fixed grid — so
#: magnitude carries no information and only injectivity matters. 0.01 is used
#: anyway because it costs nothing to buy the margin: ten rows need 0.09 of the
#: 0.99 available range, and a gap a human can see in the CSV is a gap no
#: re-rounding, spreadsheet round-trip or future grid-based sweep can close.
MIN_EPCR_SEPARATION = 0.01

#: EPCR arithmetic is done in integer units of ``10 ** -EPCR_DECIMALS`` so that
#: the separation pass is exact and byte-reproducible (GP-30). Repeated float
#: subtraction of 0.01 is neither.
_EPCR_UNIT_SCALE = 10**EPCR_DECIMALS

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
    """Rank, convert, truncate and compose candidate pairs into submission rows.

    Ordering is by EPCR descending. The scorer derives rank exactly this way
    (``sorted(enumerate(rows), key=lambda x: (-x[1][1], x[0]))``), so ties break by
    file order — which means our file order is a scientific decision, not a
    formatting one. Ties therefore fall back to
    :meth:`~mva.models.pair.CandidatePair.sort_key`, which is total: composite
    score, then both variants' coordinates, then the content-derived ``pair_id``.
    No row's position depends on the order the caller happened to supply (GP-30).

    Three composition passes then run over that order, in this sequence:

    1. :func:`_drop_subsumed` removes a candidate already fully covered by a
       *higher*-ranked one — it can contribute no variant the F-max union does not
       already hold, and rank points go to the first full match.
    2. :func:`_promote_pairs_above_subsets` handles the mirror case, where the
       subset outranks its own superset. Inside the submitted window it reorders;
       a superset ranked *below* the window is exchanged with its own subset. It
       runs over the WHOLE ranked list, **before** truncation (ADR 0023).
    3. :func:`_enforce_epcr_separation` renders that settled order as strictly
       decreasing EPCRs. It runs LAST on purpose: it is order-preserving by
       construction, so it must be handed the order we mean to submit.

    **Truncation happens between (2) and (3), and the order of (2) and truncation
    is worth 50 rank points.** Slicing to ``max_rows`` first deletes any pair that
    happens to rank below the cut — including the pair whose own subset is sitting
    at rank 1, which is precisely the pair promotion exists to rescue. The scorer
    then sees a partial match where a full match was available (100 -> 50), and if
    the pair carried the only copy of the second causal allele, F-max is capped
    below 1.0 as well. See ADR 0023.
    """
    if not 1 <= max_rows <= MAX_SUBMISSION_ROWS:
        msg = (
            f"max_rows={max_rows} is outside 1..{MAX_SUBMISSION_ROWS}. The challenge "
            f"scorer raises ValueError above {MAX_SUBMISSION_ROWS} rows."
        )
        raise ValueError(msg)

    ordered = sorted(pairs, key=lambda pair: pair.sort_key())
    ordered, subsumed = _drop_subsumed(ordered)
    if subsumed:
        # Counts only, never coordinates (GP-41).
        _LOG.info(
            "Dropped %d candidate row(s) whose variants are already covered by a "
            "higher-ranked candidate.",
            subsumed,
        )
    # Each row is carried with the candidate key that produced it so the EPCR sort
    # below is a TOTAL order in its own right. Sorting on ``-epcr`` alone would be
    # correct only by stability — it would inherit the preceding order for rows
    # whose EPCRs collide after rounding, which quietly makes this pass depend on
    # the one above it. Composite score still leads (EPCR is a strictly increasing
    # function of it), so this appends tiebreakers and re-ranks nothing.
    keyed = [(_row_for_pair(pair, proband_id), pair.sort_key()) for pair in ordered]
    keyed.sort(key=lambda item: (-item[0].epcr, item[1]))
    rows = [row for row, _ in keyed]

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
    # Compose, THEN truncate, THEN render (ADR 0023). Promotion is what decides
    # which rows belong in the window, so it cannot be handed a window that has
    # already been cut. Separation is still applied after truncation, because only
    # the submitted rows are scored and reserving range for rows that will be
    # dropped would push the kept ones needlessly close to the floor.
    composed, promoted, exchanged = _promote_pairs_above_subsets(rows, window=max_rows)
    if promoted:
        # Counts only, never coordinates (GP-41).
        _LOG.info(
            "Promoted %d pair row(s) above a single-variant row whose variant they already carry.",
            promoted,
        )
    if exchanged:
        _LOG.info(
            "Exchanged %d single-variant row(s) for a candidate pair that carries the same "
            "variant and had ranked below the %d-row limit. A one-variant row cannot equal a "
            "two-variant answer as a frozenset, so the pair is the only one of the two that "
            "can score a full match (ADR 0023).",
            exchanged,
            max_rows,
        )
    return _enforce_epcr_separation(composed[:max_rows])


def _promote_pairs_above_subsets(
    rows: Sequence[SubmissionRow], *, window: int = MAX_SUBMISSION_ROWS
) -> tuple[tuple[SubmissionRow, ...], int, int]:
    """Move a pair row above any earlier row whose variants it already carries.

    :func:`_drop_subsumed` removes a single-variant candidate that falls *below*
    the pair it was carved out of. The mirror case survives it: when the single
    outscores its own parent pair, the single is kept first and the pair lands
    beneath it. Nothing is wrong with either row — the defect is the order.

    **The cost is one-sided, so this is a straight bet.** Rank points come only
    from the best *full* match, and a match is frozenset equality on the row's
    variants. If the answer is the pair, pair-above-single scores 100 and
    single-above-pair scores 50; no arrangement puts both at rank 1. F-max is
    identical either way, because the pair re-emits the single's variant, so the
    predicted-variant union at the top threshold is the same set.

    We bet on the pair, because the challenge says so: the Track 1 answer is
    stated to be a clinically validated compound heterozygote, the scorer's own
    ``is_compound_het = len(true_variants) == 2`` gates the only partial-credit
    branch on it, and the repository's local fallback ground truth holds two
    variants.

    **The single is kept, immediately below — while keeping it is free.** A row
    ranked below the answer costs nothing (verified by executing the scorer), so
    dropping it would forfeit a full match in the unlikely world where the truth is
    a single variant and gain exactly nothing in the likely one. That argument
    holds only inside the ten-row window; see the exchange rule below for what
    happens when keeping it would cost another row its slot.

    This is composition, not scoring. Nothing here touches a composite score: the
    ranking is a scientific judgement and this is a statement about how to lay
    that judgement out in a twelve-column CSV. Implementing it as a score
    adjustment would disguise a scoring change as a rendering fix.

    **``window`` is the row limit the caller will truncate to, and this pass runs
    BEFORE that truncation (ADR 0023).** It used to be handed ``rows[:max_rows]``,
    so a pair ranked below the cut had already been deleted by the time promotion
    looked for it — and the pair most likely to sit below the cut is exactly the
    one this pass exists to rescue, because the compound-het hypothesis carries a
    phase penalty its own halves do not (GP-15). A single at rank 1 with its parent
    pair at rank 11 scored a partial match: 50 rank points instead of 100, plus an
    F-max capped below 1.0 whenever the pair held the only copy of the second
    causal allele.

    Reaching past the cut changes what promotion *costs*, so the pass is in two
    parts and they are justified differently:

    1. **Reordering inside the window is free.** No row enters or leaves, so every
       superset already submitted is lifted above every subset it covers — the
       original ADR 0015 rule, unchanged, including the case of several pairs
       carrying the same variant.
    2. **Crossing the cut costs a row, so it is rationed.** A subset row with no
       superset inside the window at all is **exchanged** with the highest-ranked
       superset outside it. The pair takes the subset's slot and the subset takes
       the pair's; exactly one row enters the submission and exactly one leaves,
       and the row that leaves is the subset itself. No unrelated row moves by a
       single position.

    The exchange is the part that needed a decision, and it inverts one sentence of
    ADR 0015. That ADR keeps the outranked single *because keeping it is free* — a
    row below the answer cannot lower either metric. Here it is not free: the only
    way to keep it is to evict something else. And a single-variant row is worth
    nothing in the world we are betting on, since ``row.variants == true_variants``
    cannot hold for a one-variant row against a two-variant answer, while the row
    it would displace is another pair that can. Its variant is not lost either —
    the promoted pair re-emits it, so the F-max union at that threshold is
    unchanged. Keep the single when it is free; drop it when it costs a slot.

    Sliding the pair up and letting the last row fall off the end was the other
    candidate rule and is worse: it spends a submitted pair to keep a single that
    can never full-match. Promoting *every* out-of-window superset is worse still —
    it demotes higher-scoring unrelated hypotheses on the strength of a
    subset relation, which is a re-ranking, and re-ranking is a scientific
    judgement that does not belong in a rendering pass.

    Terminates and is deterministic: pass 1 strictly decreases the number of
    subset-above-superset inversions inside a fixed window, and pass 2 is a single
    forward scan in which each step swaps one in-window index with one out-of-window
    index and never revisits either.
    """
    result = list(rows)
    limit = min(window, len(result))

    promoted = 0
    index = 0
    while index < limit:
        variants = result[index].variant_keys()
        superset = next(
            (later for later in range(index + 1, limit) if variants < result[later].variant_keys()),
            None,
        )
        if superset is None:
            index += 1
            continue
        result.insert(index, result.pop(superset))
        promoted += 1

    exchanged = 0
    for index in range(limit):
        variants = result[index].variant_keys()
        if any(variants < result[other].variant_keys() for other in range(limit) if other != index):
            # Already covered by a pair inside the window; pass 1 put it above.
            continue
        outside = next(
            (
                later
                for later in range(limit, len(result))
                if variants < result[later].variant_keys()
            ),
            None,
        )
        if outside is None:
            continue
        result[index], result[outside] = result[outside], result[index]
        exchanged += 1

    return tuple(result), promoted, exchanged


def _epcr_units(value: float) -> int:
    """EPCR as an exact integer count of ``10 ** -EPCR_DECIMALS``."""
    return round(value * _EPCR_UNIT_SCALE)


def _enforce_epcr_separation(
    rows: Sequence[SubmissionRow], *, separation: float = MIN_EPCR_SEPARATION
) -> tuple[SubmissionRow, ...]:
    """Force strictly-decreasing, ``separation``-apart EPCRs down the file order.

    **Order-preserving by construction, not by care.** The pass walks the rows in
    the order it is given — which is already the submitted order — and assigns
    each row a value strictly below its predecessor's. The output is therefore
    strictly decreasing *in file order* whatever the inputs were, so for any two
    rows ``i < j`` the emitted EPCRs satisfy ``epcr[i] > epcr[j]``. No row can
    overtake another; the only thing this transformation can do is turn a
    collision into a separation. That is what makes it a formatting fix rather
    than a scientific one: it cannot change which hypothesis we rank above which,
    and by the invariance of both scorer metrics under any strictly increasing
    reparameterisation of EPCR, it cannot change the score except by removing
    ties.

    Each row is assigned ``max(min(own value, predecessor - step), reserved)``
    where ``reserved = EPCR_FLOOR + step * (rows below it)``. The ``min`` enforces
    the separation; the ``reserved`` term keeps enough range for every row still
    to come, so the pass can never run out of room and re-tie two rows at the
    floor. ``reserved`` decreases by exactly ``step`` per row, so it can never
    break the separation it is protecting.

    Arithmetic is integer, in units of ``10 ** -EPCR_DECIMALS``: the rendered CSV
    is byte-identical across runs and platforms (GP-30), which repeated float
    subtraction of 0.01 does not guarantee.
    """
    if len(rows) < 2:
        return tuple(rows)

    step = _epcr_units(separation)
    if step < 1:
        msg = (
            f"separation={separation} is smaller than one rendered EPCR unit "
            f"(1e-{EPCR_DECIMALS}); it cannot separate two rounded values."
        )
        raise ValueError(msg)

    floor = _epcr_units(EPCR_FLOOR)
    ceiling = _epcr_units(EPCR_CEILING)
    if floor + step * (len(rows) - 1) > ceiling:
        msg = (
            f"{len(rows)} rows separated by {separation} do not fit in "
            f"({EPCR_FLOOR}, {EPCR_CEILING}]. Reduce the separation or the row count."
        )
        raise ValueError(msg)

    units = [_epcr_units(row.epcr) for row in rows]
    for index in range(len(units)):
        reserved = floor + step * (len(units) - 1 - index)
        upper = ceiling if index == 0 else units[index - 1] - step
        units[index] = max(min(units[index], upper), reserved)

    return tuple(
        row if unit == _epcr_units(row.epcr) else replace(row, epcr=unit / _EPCR_UNIT_SCALE)
        for row, unit in zip(rows, units, strict=True)
    )


def _drop_subsumed(
    ordered: Sequence[CandidatePair],
) -> tuple[list[CandidatePair], int]:
    """Remove candidates already fully covered by a higher-ranked candidate.

    The ranked list legitimately contains single-variant hypotheses carved out of
    a pair — a lone high-impact heterozygote is a real alternative under a
    dominant or de-novo model, so prioritisation is right to emit it. But once the
    *pair* has been submitted, the single-variant row adds nothing the scorer can
    use: F-max unions variants across rows (``predicted_variants |= row.variants``),
    so a subset row contributes no new variant, and rank points go to the first
    full match. It only consumes one of the ten available rows, which would
    otherwise carry a different gene's hypothesis.

    So: drop a candidate whose variant set is a subset of an already-kept
    candidate's. Strictly a submission-shaping step — nothing is removed from the
    ranked list, the dossier or the evidence store (GP-13, GP-19).
    """
    kept: list[CandidatePair] = []
    kept_sets: list[frozenset[str]] = []
    dropped = 0
    for pair in ordered:
        variants = frozenset(pair.variant_ids)
        if any(variants <= seen for seen in kept_sets):
            dropped += 1
            continue
        kept.append(pair)
        kept_sets.append(variants)
    return kept, dropped


def _write_rows(rows: Sequence[SubmissionRow]) -> str:
    """Serialise rows to the exact 12 columns. No checking of any kind."""
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


def render_submission_csv(rows: Sequence[SubmissionRow]) -> str:
    """Render rows to the exact CSV the challenge expects, header included.

    **This function validates. That is a deliberate public-API decision, taken
    because the alternative was an exported validation bypass.**

    ``SubmissionRow`` is a plain dataclass with no invariants of its own: it will
    hold a bare contig, a proband ID the scorer hard-fails on, an EPCR outside
    ``(0, 1]``, or a variant repeated across two rows, and a caller can hand
    eleven of them to a serialiser. Every one of those states is silent at upload
    time and fatal at scoring time. While the raw serialiser was the module's
    public renderer, `validate_submission` was reachable only by callers who
    remembered to call it, which is the same as saying the contract was optional.

    Three shapes were considered:

    * *Stop exporting the raw renderer.* Rejected: `validate_submission` has to be
      testable against deliberately-invalid bytes, and a test that reaches in for a
      private name is a bypass with extra steps.
    * *Rename it to something frightening and leave it public.* Rejected on its
      own, because the frightening name is still the shortest path from
      ``SubmissionRow`` to a file.
    * *Make the safe spelling the obvious one.* Adopted. The plain name validates
      and raises; :func:`render_submission_csv_unvalidated` keeps the raw path for
      the tests that need it, is named for what it does not do, and is excluded
      from ``mva.reporting``'s package-level exports so it cannot be reached by
      ``from mva.reporting import ...``.

    Raises:
        ValueError: listing every failing check. Errors are collected rather than
            raised one at a time so a caller sees the whole problem at once, which
            is also why :func:`validate_submission` returns a tuple rather than
            raising.
    """
    text = _write_rows(rows)
    ok, errors = validate_submission(text)
    if not ok:
        msg = (
            "Refusing to render a Track 1 submission that fails its own contract check "
            f"({len(errors)} problem(s)): " + "; ".join(errors)
        )
        raise ValueError(msg)
    return text


def render_submission_csv_unvalidated(rows: Sequence[SubmissionRow]) -> str:
    """Serialise rows with **no contract check whatsoever**. Not a submission path.

    The name is the warning, and it is the whole point of this function existing
    separately: anything produced here may be rejected outright by the scorer, or
    accepted and scored zero, and nothing in this module will have objected.

    It exists for exactly two callers:

    * tests that must feed :func:`validate_submission` bytes it is supposed to
      reject — a validator can only be shown to work on input that fails it;
    * tests that price a *hypothetical* submission shape against the scorer replica
      (a split pair, a tie, a bare contig) to record what that shape costs.

    Anything that writes a file uses :func:`write_submission`. Anything that
    produces text for a caller uses :func:`render_submission_csv`. This function is
    absent from ``mva.reporting``'s package exports on purpose.
    """
    return _write_rows(rows)


def write_submission(rows: Sequence[SubmissionRow], path: Path) -> Path:
    """Render, self-validate, write, and gate the submission as a public export.

    The self-check runs *before* the bytes hit disk: a malformed submission that
    exists is worse than one that does not, because it will be uploaded. The
    export gate then runs on the written file (GP-43) — classification is a claim,
    the re-scan is the verification — and a blocked file is removed rather than
    left for someone to find and submit.
    """
    text = _write_rows(rows)
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
    """Run the public-export gate over the written file, and act on its verdict.

    Bound dynamically because ``mva.privacy.export`` is a peer deliverable that may
    not be present in every checkout; its absence degrades to "ungated", which is
    recorded in the log rather than silently assumed to be fine.

    Two things used to make this function decorative, and both were exercised in
    review:

    * ``gate_public_export`` RETURNS an :class:`~mva.privacy.export.ExportDecision`
      and never raises for a refusal — that is its documented contract, so the
      caller can see every failing check at once. Wrapping it in ``try/except``
      and ignoring the return value meant a verdict of ``allowed=False`` deleted
      nothing, raised nothing and left the submission on disk, while the docstring
      claimed the opposite.
    * ``allowlist=(path.name,)`` builds the allowlist out of the file being
      checked, so the allowlist check could not fail. The deny-by-default list
      is :data:`mva.privacy.export.PUBLIC_EXPORT_ALLOWLIST`, which is what makes
      it a decision taken in advance rather than a tautology.

    The ``try/except`` is kept, narrowed to the call itself, because a gate that
    raises for an unexpected reason must still delete the file before the error
    propagates. Fail closed means the artifact does not survive the refusal.
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
    declared_allowlist = getattr(
        module, "PUBLIC_EXPORT_ALLOWLIST", TRACK1_EXPORT_ALLOWLIST_FALLBACK
    )
    allowlist = cast("Sequence[str]", declared_allowlist)
    try:
        decision = cast(
            "ExportDecision", gate(path, declared=Sensitivity.PUBLIC, allowlist=allowlist)
        )
    except Exception:
        path.unlink(missing_ok=True)
        raise
    if not decision.allowed:
        path.unlink(missing_ok=True)
        msg = (
            f"Submission {path.name!r} was refused by the public-export gate (GP-43) and "
            f"has been deleted. Failing checks: {'; '.join(decision.reasons)}. "
            "Content is withheld under GP-41."
        )
        raise ExportBlockedError(msg)


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
    singles: list[_SingleRow] = []
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
            if len(key) == 1:
                gene = _note_gene(row.get("notes") or "")
                finding = (row.get("finding_type") or "").strip()
                if gene and finding != "secondary":
                    singles.append(_SingleRow(index=index, variant=next(iter(key)), gene=gene))
    _validate_split_pairs(singles, frozenset(seen_keys), errors)

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
    _validate_epcr_ties(ordered, errors)

    return (not errors, tuple(errors))


def _validate_epcr_ties(epcrs: Sequence[float], errors: list[str]) -> None:
    """Reject a repeated EPCR value.

    A tie is the one property of the EPCR vector, other than its order, that the
    scorer can see, and it is pure loss. Thresholds in the F-max sweep are the
    unique emitted values, so two rows sharing one can never be separated: at
    every threshold the sweep either predicts both or neither, and the tied
    partner's variants are false positives at the true row's own threshold. If
    file order also puts the wrong row first, the true row is demoted a rank.
    Nothing is bought in exchange, which is why this is an error and not a note.
    """
    seen: dict[str, int] = {}
    for index, value in enumerate(epcrs, start=1):
        rendered = f"{value:.{EPCR_DECIMALS}f}"
        first = seen.setdefault(rendered, index)
        if first != index:
            errors.append(
                f"row {index}: epcr={rendered} repeats the value already used by row "
                f"{first}. Tied rows are inseparable at every F-max threshold and the "
                "later row is demoted by file order; emit strictly decreasing values "
                f"at least {MIN_EPCR_SEPARATION} apart."
            )


@dataclass(frozen=True)
class _SingleRow:
    """A parsed single-variant row, with the gene the pipeline attributed it to."""

    index: int
    variant: tuple[str, int, str, str]
    gene: str


def _note_gene(note: str) -> str:
    """The gene symbol :func:`_safe_note` writes first, or ``""`` if there is none.

    ``notes`` is the only place a rendered row carries gene identity: the twelve
    columns the challenge specifies have no gene field, so a validator working on
    CSV text has nothing else to key on. :func:`_safe_note` always writes
    ``"<gene>; <inheritance>; phase=<status>"``, so the first ``;``-separated token
    is the symbol.

    Deliberately conservative. A token containing ``=`` is one of the later fields
    arriving first, and an empty ``notes`` (a hand-built row) yields ``""``. Both
    mean "no attribution", and a row with no attribution is skipped by the
    split-pair check rather than guessed at — the check is best-effort on an
    optional field, and a false positive there would refuse a legitimate
    submission.
    """
    token = note.split(";", 1)[0].strip()
    return "" if "=" in token else token


def _validate_split_pairs(
    singles: Sequence[_SingleRow],
    proposed: frozenset[frozenset[tuple[str, int, str, str]]],
    errors: list[str],
) -> None:
    """Reject one compound-het hypothesis taken apart into two single-variant rows.

    **The gap this closes.** The other cross-row checks look for rows that are too
    similar: the same variant set twice, the same EPCR twice. This one looks for
    the opposite failure — two rows that are each individually perfect and are
    jointly a downgrade. Right proband, ``chr``-prefixed contigs, in-range and
    untied EPCRs, distinct variant sets: nothing per-row is wrong, and the
    submission still cannot score what it should.

    **What it costs, verified.** The scorer matches with
    ``row.variants == true_variants`` — frozenset equality. A one-variant row can
    never equal a two-variant answer, so once the pair is split there is no full
    match to be had anywhere in the file and the best remaining outcome is the
    partial-credit branch: **100 -> 50 rank points**. F-max is untouched (it unions
    variants across rows, so both alleles are still predicted), which is precisely
    what makes the defect invisible without this check — the submission looks right
    on one metric while silently halving the other. Contract, "Never split a
    candidate pair across two rows".

    **Why the gene and not the coordinates.** Two variants are a compound
    heterozygote when they sit in the *same gene*; the check therefore fires on two
    single-variant rows the pipeline attributed to one gene and no row proposing
    them together. Requiring every unordered pair of single-variant rows to be
    joined would demand 45 pair rows for ten singles, which is neither possible nor
    a claim we would want to make. Gene identity is read from ``notes`` (see
    :func:`_note_gene`), the only column that carries it; rows with no attribution
    are skipped rather than guessed at.

    **Secondary findings are exempt.** ``finding_type=secondary`` states that the
    row is an incidental finding *unrelated to the primary phenotype*. Joining it
    to a primary hypothesis in the same gene would assert a compound heterozygote
    the pipeline explicitly did not claim.
    """
    for position, left in enumerate(singles):
        for right in singles[position + 1 :]:
            if left.gene != right.gene:
                continue
            if frozenset({left.variant, right.variant}) in proposed:
                continue
            errors.append(
                f"rows {left.index} and {right.index}: two single-variant rows in "
                f"{left.gene} are one compound-heterozygous hypothesis split across two "
                "rows, and no row proposes the two variants together. The scorer matches "
                "by frozenset equality, so neither row can equal a two-variant answer and "
                "no full match exists anywhere in the file: the best available outcome "
                "drops to partial credit, 100 -> 50 rank points, for a prediction that is "
                "otherwise exactly right (F-max is unaffected, which is why this is "
                "silent). A compound-het pair is ONE row, using the _2 columns."
            )


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
    "MIN_EPCR_SEPARATION",
    "TRACK1_COLUMNS",
    "SubmissionRow",
    "build_submission_rows",
    "composite_to_epcr",
    "render_submission_csv",
    # Deliberately module-level only: `mva.reporting` does not re-export it, so the
    # unvalidated path cannot be reached by `from mva.reporting import ...`.
    "render_submission_csv_unvalidated",
    "truncation_notice",
    "validate_submission",
    "write_submission",
]
