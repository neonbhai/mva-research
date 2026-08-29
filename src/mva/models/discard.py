"""The record of a variant that did not make it, and the contract for writing one.

GP-13's promise is that a candidate never becomes invisible. Until ADR 0028 the
pipeline kept that promise by *retaining* records — everything that was not
invalid travelled all the way to pairing, carrying flags. That is affordable at
fixture scale and not at whole-genome scale: 4,962,060 variants each cost a
gnomAD lookup, a ClinVar lookup, ~6 KB of JSON and roughly three ``EvidenceItem``
rows, and a measured run projected ~19.5 h and ~95 GB of intermediates on a 60 GB
volume for a callset in which 14,228 variants — 0.284% — could possibly matter.

The amended GP-13 asks for something weaker and sufficient: every dropped record
is **accounted for**, by id, with the rule that dropped it and the value that
triggered it. That is what a :class:`DiscardRecord` is, and it costs about 75
bytes rather than about 13 KB.

Three properties are deliberate:

**It is not an ``EvidenceItem``.** An evidence item carries a claim, a tier, a
strength, a direction, limitations prose and a citation, and it exists to be
cited by a report (GP-10). Nothing cites a discard. ADR 0019 already refused to
emit four million rejection items for exactly this reason and settled for
per-reason counts; this is the middle term it could not afford — per-variant, but
not per-variant *prose*.

**It is a slots dataclass, not a Pydantic model.** ``VariantRecord`` costs 3,672
bytes resident (``docs/scale-report.md`` §2). At 5 M discards the model tax alone
would be the thing this record exists to avoid. Nothing validates a discard row
because nothing constructs one from untrusted input: the writer is the pipeline.

**The sink is a Protocol declared here, in layer 0.** ``mva.annotation`` (layer 4)
and ``mva.prioritization`` (layer 5) both drop records; the ledger that persists
them is ``mva.evidence`` (layer 6). Neither stage may import it (GP-01), and
neither should know it exists. They take a :class:`DiscardSink` and the
composition root passes the real one.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Final, Protocol, runtime_checkable

#: Column separator and row terminator. A discard row is a TSV line because the
#: format has to survive `grep`: the question this ledger answers — "why isn't
#: variant X in the output?" — is one a reviewer should be able to ask of a 500 MB
#: file without loading it, and `mva discards --variant` is a convenience over
#: exactly the same scan.
FIELD_SEPARATOR: Final = "\t"

#: Written once at the top of a ledger file. Comment-prefixed so that a reader
#: splitting on tabs can skip it on the same rule that skips the provenance
#: banner above it.
DISCARD_LEDGER_COLUMNS: Final[tuple[str, ...]] = ("variant_id", "stage", "rule", "value")


class DiscardStage(StrEnum):
    """Which stage dropped the record.

    Ordered by where they sit in the pass, so a ledger read top to bottom reads
    in pipeline order for any one variant — though a variant appears at most
    once, because the first stage to drop it is the last one to see it.
    """

    INGEST = "ingest"
    """The reader refused the record before it became a ``VariantRecord`` at all:
    a non-canonical contig, a symbolic or structural allele, a spanning deletion,
    a missing ALT, a malformed line, an uninterpretable coordinate.

    This stage exists because a discard audit that begins downstream of the first
    gate is complete about everything except the largest single loss. On this
    callset that gate removed 116,485 decoy and unplaced records — more than eight
    times the number of variants that survive to the end — and until now they were
    a count in ``qc_report.json`` with no identities behind it. "Never extracted"
    and "extracted, then dropped" are different facts about a variant and stay
    separately answerable here.

    A record refused at this gate may have no canonical id, because forming one is
    the thing that failed. The ``variant_id`` field then carries the best handle
    that exists — a build-qualified coordinate where the coordinate parsed, and
    ``line:<index>`` where it did not — which is enough to find the line in the
    source VCF and is never presented as a canonical id."""

    TRIAGE = "triage"
    """The annotation stage's relevance bar (ADR 0028). The record was annotated
    for consequence and judged not worth a frequency lookup, a clinical lookup,
    an evidence row or a line in the variants artifact."""

    HARD_FILTER = "hard_filter"
    """Invalid or impossible: wrong build, non-canonical contig, placeholder
    allele, hom-ref, no-call (ADR 0005). These were always deleted; they were
    previously counted and not named."""

    SELECTION = "selection"
    """The candidate-selection stage in front of pairing (ADR 0019). The one
    stage entitled to delete a *valid* record, and until now the one whose
    casualties were counted in aggregate only."""


@dataclass(frozen=True, slots=True)
class DiscardRecord:
    """One dropped variant: what it was, what dropped it, and on what value.

    ``value`` is the observation that triggered the rule, rendered as a short
    string — an impact class, an allele frequency, a contig name. It is the
    difference between a ledger that says "dropped by selection" and one that
    says "dropped because an adequately powered population reports AF 0.031", and
    the second is the only one a reviewer can argue with.

    ``value`` may be empty. Some rules fire on an absence (no consequence
    annotation at all), and an empty string is the honest rendering of that —
    never ``"0"``, never ``"none"`` (GP-14).
    """

    variant_id: str
    stage: DiscardStage
    rule: str
    value: str = ""

    def as_row(self) -> str:
        """One TSV line, terminated. Fields are sanitised, never quoted.

        A tab or a newline inside a field would silently shift every later column
        and make the file unparseable at exactly the row a reviewer is looking
        for. None of the four fields can legitimately contain either — a variant
        id is ``build:contig:pos:ref:alt``, a rule is an identifier, a value is a
        number or a short token — so the whitespace is collapsed rather than
        escaped, which keeps the reader a `str.split` for a 5 M-row file.
        """
        fields = (self.variant_id, str(self.stage), self.rule, self.value)
        return FIELD_SEPARATOR.join(_sanitise(field) for field in fields) + "\n"

    @classmethod
    def from_row(cls, line: str) -> DiscardRecord | None:
        """Parse one line, or ``None`` for a comment, a blank or a malformed row.

        Lenient on purpose. This is the read side of a file that may have been
        truncated by the very disk exhaustion the stage was trying to avoid, and
        a reader that raises on the last half-written line answers nothing about
        the 4,999,999 rows in front of it.
        """
        stripped = line.rstrip("\n")
        if not stripped or stripped.startswith("#"):
            return None
        parts = stripped.split(FIELD_SEPARATOR)
        if len(parts) != len(DISCARD_LEDGER_COLUMNS):
            return None
        variant_id, stage, rule, value = parts
        try:
            parsed_stage = DiscardStage(stage)
        except ValueError:
            return None
        return cls(variant_id=variant_id, stage=parsed_stage, rule=rule, value=value)

    @property
    def counter_key(self) -> str:
        """``stage/rule``, the key a per-reason count is kept under."""
        return f"{self.stage}/{self.rule}"


def _sanitise(field: str) -> str:
    """Collapse anything that would break the row shape into a single space."""
    if "\t" not in field and "\n" not in field and "\r" not in field:
        return field
    return field.replace("\t", " ").replace("\r", " ").replace("\n", " ")


@runtime_checkable
class DiscardSink(Protocol):
    """Anywhere a dropped record can be accounted for.

    Deliberately one method and no lifecycle. A stage that drops a record should
    not be able to open, flush, close or query the ledger; it should be able to
    say what it dropped and why, and nothing else. ``mva.evidence.discards``
    supplies the implementation, and ``None`` is always a legal value for a caller
    that has no ledger — every stage's pre-ADR-0028 behaviour.
    """

    def record(self, discard: DiscardRecord) -> None:
        """Account for one dropped variant. Must not raise for a well-formed record."""
        ...
