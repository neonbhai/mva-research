"""GP-10 in executable form: no claim reaches a report without evidence behind it.

An :class:`Assertion` is a *scientific statement about the patient or the
biology* — the kind of sentence a clinician would underline. Every one carries
the tier it was believed at and the evidence IDs it rests on, and every one is
run past an :class:`AssertionChecker` before it can be rendered.

What is deliberately **not** an Assertion: structural readouts of the pipeline's
own state ("the composite score is 0.71", "phase is UNKNOWN", "three open
questions are blocking"). Those are facts about this program, not claims about
biology, and the report labels them as such rather than laundering them through
the evidence gate. The distinction is the whole point — if everything is an
assertion, the marker stops meaning anything.

Two visible outcomes:

* An assertion citing nothing, or citing an ID the ledger does not hold, is
  refused with :class:`~mva.errors.UnsourcedAssertionError`. It is not downgraded,
  hedged or emitted with a warning; a warning inside a clinical-adjacent report
  is a warning nobody reads.
* An assertion at an unproven tier (``INFERENCE``, ``SPECULATION``) is rendered
  with its marker **leading** the sentence, so a reader skimming the left margin
  cannot miss it. Proven tiers carry a trailing marker, which labels without
  shouting.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING

from mva.errors import UnsourcedAssertionError
from mva.models import UNPROVEN_TIERS, AssertionTier

if TYPE_CHECKING:  # pragma: no cover - typing only, no runtime dependency
    from mva.evidence.ledger import AssertionResolver
    from mva.models import EvidenceItem

#: Visible label attached to every rendered claim. Tiers in
#: :data:`~mva.models.UNPROVEN_TIERS` get a marker that survives skim-reading;
#: ``SPECULATION`` is shouted because it is the one tier that must never be
#: mistaken for a finding.
TIER_MARKERS: dict[AssertionTier, str] = {
    AssertionTier.OBSERVED_DATA: "[observed]",
    AssertionTier.DATABASE_ASSERTION: "[database]",
    AssertionTier.COMPUTATIONAL_PREDICTION: "[predicted]",
    AssertionTier.LITERATURE_MECHANISM: "[literature]",
    AssertionTier.INFERENCE: "[inferred]",
    AssertionTier.SPECULATION: "[SPECULATIVE]",
}

#: Appended when the cited evidence itself contains a contradiction (GP-19). A
#: claim supported by evidence that partly argues against it must not read as
#: settled.
CONTRADICTED_MARKER = "[contradicted]"

#: The epistemic ladder, in declaration order: observed -> speculation.
_TIER_LADDER: tuple[AssertionTier, ...] = tuple(AssertionTier)


def weakest_tier(items: Sequence[EvidenceItem]) -> AssertionTier:
    """The lowest rung any cited evidence stands on.

    A claim is only as strong as its weakest support, so a sentence resting on one
    curated database record and one pipeline inference is labelled ``[inferred]``,
    not ``[database]``. Taking the maximum here instead would let a single strong
    citation launder a chain of guesses.
    """
    if not items:
        msg = "weakest_tier requires at least one evidence item."
        raise ValueError(msg)
    return max((item.tier for item in items), key=_TIER_LADDER.index)


@dataclass(frozen=True)
class Assertion:
    """One evidence-backed sentence destined for a report."""

    text: str
    tier: AssertionTier
    evidence_ids: tuple[str, ...] = ()
    contradicted: bool = False
    """Set by :meth:`AssertionChecker.check` when the cited evidence includes an
    item whose direction is ``CONTRADICTS``. Never set by hand."""

    def rendered(self) -> str:
        """The sentence with its tier marker, positioned by how much it needs one."""
        marker = TIER_MARKERS[self.tier]
        body = f"{marker} {self.text}" if self.tier in UNPROVEN_TIERS else f"{self.text} {marker}"
        return f"{body} {CONTRADICTED_MARKER}" if self.contradicted else body

    def citation(self) -> str:
        """Evidence IDs in a stable, sorted order for rendering next to the claim."""
        return ", ".join(sorted(self.evidence_ids)) if self.evidence_ids else "none"

    @property
    def is_unproven(self) -> bool:
        return self.tier in UNPROVEN_TIERS


class AssertionChecker:
    """The gate every rendered claim passes through.

    ``strict=True`` (the default, and the only mode used for artifacts that leave
    the machine) refuses uncited claims and claims citing IDs the ledger has never
    seen. ``strict=False`` exists for drafting and for diagnostics: it resolves
    what it can, marks contradictions, and never raises. A non-strict render is
    not a publishable artifact and the report templates say so.
    """

    def __init__(self, resolver: AssertionResolver, *, strict: bool = True) -> None:
        self._resolver = resolver
        self._strict = strict

    @property
    def strict(self) -> bool:
        return self._strict

    def check(self, assertion: Assertion) -> Assertion:
        """Verify one assertion and return it, annotated.

        Raises :class:`~mva.errors.UnsourcedAssertionError` in strict mode when the
        claim cites nothing or cites an unresolvable ID.
        """
        if not assertion.evidence_ids:
            if not self._strict:
                return assertion
            msg = (
                f"Unsourced assertion (GP-10): {_excerpt(assertion.text)!r} cites no "
                f"evidence and cannot be rendered. Either attach the EvidenceItem ID "
                f"that supports it, or record it as an OpenQuestion / uncertainty so "
                f"the report states the gap instead of stating the claim."
            )
            raise UnsourcedAssertionError(msg)

        items = (
            self._resolver.require(assertion.text, assertion.evidence_ids)
            if self._strict
            else self._resolver.resolve(assertion.evidence_ids)
        )
        contradicted = any(item.is_contradiction for item in items)
        if contradicted == assertion.contradicted:
            return assertion
        return replace(assertion, contradicted=contradicted)

    def check_all(self, assertions: Sequence[Assertion]) -> tuple[Assertion, ...]:
        """Check every assertion, in order. Fails on the first violation.

        Fail-fast rather than collect-and-report: a report with one unsourced
        sentence is not 90% publishable, it is unpublishable, and the first failure
        names the sentence that has to change.
        """
        return tuple(self.check(assertion) for assertion in assertions)

    def rendered_all(self, assertions: Sequence[Assertion]) -> tuple[str, ...]:
        """Check, then render. The only path templates use to obtain claim text."""
        return tuple(checked.rendered() for checked in self.check_all(assertions))


def _excerpt(claim: str, limit: int = 120) -> str:
    stripped = claim.strip()
    return stripped if len(stripped) <= limit else stripped[: limit - 1] + "…"


__all__ = [
    "CONTRADICTED_MARKER",
    "TIER_MARKERS",
    "Assertion",
    "AssertionChecker",
    "weakest_tier",
]
