"""Four-valued observation logic propagated over the HPO DAG.

This is the module where the two easiest ways to break a phenotype pipeline live
next to each other, so both are spelled out.

**The directions are opposite, and that is not a symmetry to tidy away.**

``OBSERVED`` propagates **upward**.
    Recording "Microcephaly" entails "Abnormality of skull size" and every other
    ``is_a`` ancestor: if the specific thing is present, the general thing is
    present. It entails nothing about siblings and nothing about children —
    microcephaly does not imply *Severe microcephaly*, and it says nothing about
    *Macrocephaly*. Propagating downward would credit a gene for a specific
    feature that was never seen, which inflates every score in the run.

``EXCLUDED`` propagates **downward**.
    Recording "no abnormality of the eye" entails "no cataract" and every other
    descendant: if the general thing is absent, every specific form of it is
    absent. It entails nothing about ancestors — "no cataract" plainly does not
    mean "no abnormality of the eye". Propagating a negative upward would
    manufacture a broad negative finding out of one narrow one and would then
    penalise every gene under that branch.

``UNCERTAIN`` and ``NOT_ASSESSED`` propagate in **neither** direction.
    An equivocal finding does not entail its ancestors, and "nobody looked" cannot
    entail anything at all. They stay exactly where they were recorded (GP-14).

**Contradictory records are reported, not resolved.** There are two shapes and
both land in :attr:`PropagatedProfile.conflicts`:

*Across two terms* — a profile states "Cataract: observed" and "Abnormality of the
eye: excluded". Those cannot both be true. The explicit clinical assertion is kept
for the term it was made about, and the disagreement is reported.

*Within one term, under two spellings* — a profile states ``HP:0000002``
(current id) EXCLUDED and ``HP:0009999`` (its retired ``alt_id``) OBSERVED. Both
name the same term. This one is not resolved at all: :func:`_reconcile` demotes
the term to ``UNCERTAIN`` so it contributes nothing in either direction and
propagates nowhere. Ranking the two statuses and keeping the winner — the
behaviour this module used to have — converted a genuine contradiction into
positive evidence *and* made it vanish from the conflict report, which is the
precise inversion the four-valued logic exists to prevent.

Silently picking a winner would hide a data-entry error that a human needs to see.

**Unresolvable terms are an information gap.** A term the loaded release does not
contain — obsoleted without a successor, or newer than the ontology — is recorded
in :attr:`PropagatedProfile.unresolved_terms` and contributes nothing. It is not
scored as a failed match.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum

from mva.models.phenotype import ObservationStatus, PhenotypeProfile
from mva.phenotype.ontology import HpoOntology


class InferenceBasis(StrEnum):
    """Why a term carries the status it does. Always reported alongside the status.

    A match inferred through the graph and a match recorded by a clinician are not
    the same quality of evidence, and a report that cannot tell them apart is
    overstating its case.
    """

    RECORDED = "recorded"
    """The profile states this exact term's status."""

    ANCESTOR_OF_OBSERVED = "ancestor_of_observed"
    """Entailed upward: a descendant of this term was recorded as OBSERVED."""

    DESCENDANT_OF_EXCLUDED = "descendant_of_excluded"
    """Entailed downward: an ancestor of this term was recorded as EXCLUDED."""

    UNRECORDED = "unrecorded"
    """Neither recorded nor entailed. Always NOT_ASSESSED (GP-14)."""

    UNRESOLVED = "unresolved"
    """Not present in the loaded ontology release; no inference is possible."""


@dataclass(frozen=True, slots=True)
class InferredObservation:
    """The status of one term for one subject, with its justification."""

    hpo_id: str
    status: ObservationStatus
    basis: InferenceBasis
    source_term: str | None
    """The recorded term the inference came from. ``None`` when ``basis`` is
    ``RECORDED``, ``UNRECORDED`` or ``UNRESOLVED``."""

    @property
    def is_inferred(self) -> bool:
        return self.basis in (
            InferenceBasis.ANCESTOR_OF_OBSERVED,
            InferenceBasis.DESCENDANT_OF_EXCLUDED,
        )


@dataclass(frozen=True, slots=True)
class PropagationSummary:
    """Counts describing what propagation added, for rationale and evidence."""

    recorded_observed: int
    recorded_excluded: int
    recorded_uncertain: int
    recorded_not_assessed: int
    implied_observed: int
    """Ancestors gained by upward closure, excluding the recorded terms themselves."""
    implied_excluded: int
    """Descendants gained by downward closure, excluding the recorded terms themselves."""
    unresolved: int
    conflicts: int
    """Distinct terms in either conflict category below."""
    direct_conflicts: int
    """Terms the record asserts BOTH present and absent under two spellings of the
    same identifier. Never resolved by precedence; see :func:`_reconcile`."""
    entailed_conflicts: int
    """Terms the two closures disagree about: an observed finding sitting under an
    ancestor that was assessed and excluded."""


class PropagatedProfile:
    """A :class:`PhenotypeProfile` with its entailments materialised over the DAG.

    Built once per subject per run and passed to the scorer. All accessors return
    sorted tuples, and nothing iterates a set in output position (GP-30).
    """

    __slots__ = (
        "_conflicts",
        "_direct_conflicts",
        "_entailed_conflicts",
        "_excluded_closure",
        "_observed_closure",
        "_ontology",
        "_recorded",
        "_summary",
        "_unresolved",
    )

    def __init__(self, profile: PhenotypeProfile, *, ontology: HpoOntology) -> None:
        # Group by RESOLVED primary identifier before reconciling anything. The
        # loader already refuses two rows carrying the same literal HPO id, so the
        # only way one term reaches here twice is under two spellings — a primary
        # id and one of its alt_ids. Reconciling those by "most informative status
        # wins" silently turned "EXCLUDED under the retired id, OBSERVED under the
        # current one" into plain OBSERVED, which is the four-valued logic
        # inverting itself on exactly the input alt_id exists to describe.
        by_primary: dict[str, list[ObservationStatus]] = {}
        unresolved: list[str] = []
        for observation in profile.observations:
            resolved = ontology.resolve(observation.hpo_id)
            if resolved is None:
                unresolved.append(observation.hpo_id)
                continue
            by_primary.setdefault(resolved, []).append(observation.status)

        recorded: dict[str, ObservationStatus] = {}
        direct_conflicts: list[str] = []
        for primary, statuses in by_primary.items():
            status, conflicted = _reconcile(statuses)
            recorded[primary] = status
            if conflicted:
                direct_conflicts.append(primary)

        observed_seeds = sorted(
            term for term, status in recorded.items() if status is ObservationStatus.OBSERVED
        )
        excluded_seeds = sorted(
            term for term, status in recorded.items() if status is ObservationStatus.EXCLUDED
        )

        # UP for observed, DOWN for excluded. Swapping these two calls is the single
        # highest-impact bug this module can contain, and it is tested in both
        # directions in tests/unit/test_phenotype_propagation.py.
        observed_closure: dict[str, str] = {}
        for seed in observed_seeds:
            for ancestor in ontology.ancestor_closure(seed):
                observed_closure.setdefault(ancestor, seed)

        excluded_closure: dict[str, str] = {}
        for seed in excluded_seeds:
            for descendant in ontology.descendant_closure(seed):
                excluded_closure.setdefault(descendant, seed)

        entailed_conflicts = sorted(set(observed_closure) & set(excluded_closure))

        self._ontology: HpoOntology = ontology
        self._recorded: dict[str, ObservationStatus] = dict(sorted(recorded.items()))
        self._observed_closure: dict[str, str] = dict(sorted(observed_closure.items()))
        self._excluded_closure: dict[str, str] = dict(sorted(excluded_closure.items()))
        self._unresolved: tuple[str, ...] = tuple(sorted(set(unresolved)))
        self._direct_conflicts: tuple[str, ...] = tuple(sorted(direct_conflicts))
        self._entailed_conflicts: tuple[str, ...] = tuple(entailed_conflicts)
        self._conflicts: tuple[str, ...] = tuple(
            sorted(set(direct_conflicts) | set(entailed_conflicts))
        )
        self._summary: PropagationSummary = PropagationSummary(
            recorded_observed=len(observed_seeds),
            recorded_excluded=len(excluded_seeds),
            recorded_uncertain=sum(
                1 for status in recorded.values() if status is ObservationStatus.UNCERTAIN
            ),
            recorded_not_assessed=sum(
                1 for status in recorded.values() if status is ObservationStatus.NOT_ASSESSED
            ),
            implied_observed=len(observed_closure) - len(observed_seeds),
            implied_excluded=len(excluded_closure) - len(excluded_seeds),
            unresolved=len(self._unresolved),
            conflicts=len(self._conflicts),
            direct_conflicts=len(self._direct_conflicts),
            entailed_conflicts=len(self._entailed_conflicts),
        )

    # -- status -------------------------------------------------------------

    def status_of(self, hpo_id: str) -> InferredObservation:
        """The subject's status for ``hpo_id``, recorded or entailed.

        Precedence, in order, and each step is a claim about evidence quality:

        0. **Reconciliation first.** Rows are grouped by resolved primary id and
           collapsed by :func:`_reconcile` before any of the steps below run, so
           "recorded" below always means one reconciled status per term. A term
           the record asserts both present and absent is ``UNCERTAIN`` by the time
           it reaches step 1, and therefore falls through to step 4.
        1. **A recorded OBSERVED or EXCLUDED** — a clinician's direct assertion
           about this exact term outranks anything derived from another term.
        2. **Entailed OBSERVED** — a descendant was seen, so this term is present by
           the meaning of ``is_a``. This outranks a recorded ``UNCERTAIN`` or
           ``NOT_ASSESSED``: if the specific finding was made, the general one is
           not "unassessed" whatever the form said.
        3. **Entailed EXCLUDED** — an ancestor was assessed and found absent.
           Ranked below entailed OBSERVED, so a directly seen finding is never
           overturned by a broader negative.
        4. **A recorded UNCERTAIN**, then **a recorded NOT_ASSESSED**.
        5. **NOT_ASSESSED**, unrecorded. The safe default: a term nobody wrote down
           is a term nobody assessed.

        A term outside the ontology release returns ``NOT_ASSESSED`` with basis
        ``UNRESOLVED`` — no information, in either direction.
        """
        resolved = self._ontology.resolve(hpo_id)
        if resolved is None:
            return InferredObservation(
                hpo_id=hpo_id.strip().upper().replace("HP_", "HP:"),
                status=ObservationStatus.NOT_ASSESSED,
                basis=InferenceBasis.UNRESOLVED,
                source_term=None,
            )

        recorded = self._recorded.get(resolved)
        if recorded in (ObservationStatus.OBSERVED, ObservationStatus.EXCLUDED):
            assert recorded is not None
            return InferredObservation(
                hpo_id=resolved,
                status=recorded,
                basis=InferenceBasis.RECORDED,
                source_term=None,
            )

        seed = self._observed_closure.get(resolved)
        if seed is not None:
            return InferredObservation(
                hpo_id=resolved,
                status=ObservationStatus.OBSERVED,
                basis=InferenceBasis.ANCESTOR_OF_OBSERVED,
                source_term=seed,
            )

        seed = self._excluded_closure.get(resolved)
        if seed is not None:
            return InferredObservation(
                hpo_id=resolved,
                status=ObservationStatus.EXCLUDED,
                basis=InferenceBasis.DESCENDANT_OF_EXCLUDED,
                source_term=seed,
            )

        if recorded is not None:
            return InferredObservation(
                hpo_id=resolved,
                status=recorded,
                basis=InferenceBasis.RECORDED,
                source_term=None,
            )
        return InferredObservation(
            hpo_id=resolved,
            status=ObservationStatus.NOT_ASSESSED,
            basis=InferenceBasis.UNRECORDED,
            source_term=None,
        )

    # -- views --------------------------------------------------------------

    @property
    def recorded_observed(self) -> tuple[str, ...]:
        """Terms the clinician recorded as present, sorted. No entailments."""
        return tuple(
            term for term, status in self._recorded.items() if status is ObservationStatus.OBSERVED
        )

    @property
    def recorded_excluded(self) -> tuple[str, ...]:
        """Terms the clinician recorded as absent, sorted. No entailments."""
        return tuple(
            term for term, status in self._recorded.items() if status is ObservationStatus.EXCLUDED
        )

    @property
    def observed_closure(self) -> tuple[str, ...]:
        """Recorded OBSERVED terms **plus all their ancestors**, sorted."""
        return tuple(self._observed_closure)

    @property
    def excluded_closure(self) -> tuple[str, ...]:
        """Recorded EXCLUDED terms **plus all their descendants**, sorted."""
        return tuple(self._excluded_closure)

    @property
    def unresolved_terms(self) -> tuple[str, ...]:
        """Recorded terms this ontology release does not contain, sorted."""
        return self._unresolved

    @property
    def conflicts(self) -> tuple[str, ...]:
        """Every term the record contradicts itself about, sorted.

        The union of :attr:`direct_conflicts` and :attr:`entailed_conflicts`.
        Non-empty means the source record is internally inconsistent. Reported, not
        silently repaired.
        """
        return self._conflicts

    @property
    def direct_conflicts(self) -> tuple[str, ...]:
        """Terms recorded BOTH observed and excluded, under two spellings of one id.

        These never reach :attr:`observed_closure` or :attr:`excluded_closure`:
        :func:`_reconcile` demotes them to ``UNCERTAIN``, so they contribute
        nothing in either direction and propagate nowhere until a human resolves
        the record.
        """
        return self._direct_conflicts

    @property
    def entailed_conflicts(self) -> tuple[str, ...]:
        """Terms the upward and downward closures disagree about, sorted."""
        return self._entailed_conflicts

    @property
    def summary(self) -> PropagationSummary:
        return self._summary

    @property
    def ontology(self) -> HpoOntology:
        return self._ontology

    def __repr__(self) -> str:
        return (
            f"PropagatedProfile(observed={self._summary.recorded_observed}"
            f"+{self._summary.implied_observed} implied, "
            f"excluded={self._summary.recorded_excluded}"
            f"+{self._summary.implied_excluded} implied, "
            f"unresolved={self._summary.unresolved}, conflicts={self._summary.conflicts})"
        )


def _reconcile(statuses: Sequence[ObservationStatus]) -> tuple[ObservationStatus, bool]:
    """Collapse the statuses recorded for **one** primary term.

    Returns ``(status, is_conflicted)``.

    A term reaches here more than once only when the profile spelled it two ways —
    a current identifier and one of its ``alt_id`` aliases — because
    :func:`mva.phenotype.loader.load_phenotype_profile` already refuses two rows
    with the same literal id. ``alt_id`` exists precisely because real annotations
    accumulate under retired identifiers, so "EXCLUDED under the old spelling,
    OBSERVED under the new one" is a realistic export, not a curiosity.

    The rule, and why each branch is what it is:

    * **Unanimous** — collapse. Nothing is being decided.
    * **More than one distinct informative status** (``OBSERVED`` *and*
      ``EXCLUDED``) — an irreconcilable clinical contradiction about a single
      term. The result is ``UNCERTAIN`` and ``is_conflicted`` is ``True``.
      ``UNCERTAIN`` is the honest bucket: the feature *was* assessed, twice, and
      the assessments disagree. It is also the only outcome that contributes
      nothing in either direction and propagates nowhere, so a contradictory
      record can neither support nor argue against a gene until a human resolves
      it. An earlier version ranked ``OBSERVED`` above ``EXCLUDED`` and kept it;
      that turned a contradiction into positive evidence and made the conflict
      disappear from the report entirely.
    * **Exactly one informative status, plus uninformative duplicates** — the
      informative one wins and this is *not* a conflict. "Assessed and present"
      together with "nobody assessed it" is a redundant row, not a disagreement:
      nothing contradicts the finding, so discarding it would lose real evidence
      to a clerical artifact.
    * **All uninformative** — ``UNCERTAIN`` beats ``NOT_ASSESSED``. Both carry
      zero weight in both directions, so no evidence can be created or destroyed
      here; the more specific description of the gap is kept.

    Only identical or uninformative duplicates are therefore ever collapsed
    silently. Anything else is either reported as a conflict or is a strict
    superset of the information in the rows it replaces.
    """
    distinct = set(statuses)
    if len(distinct) == 1:
        return statuses[0], False

    informative = distinct & _INFORMATIVE_STATUSES
    if len(informative) > 1:
        return ObservationStatus.UNCERTAIN, True
    if len(informative) == 1:
        return next(iter(informative)), False
    if ObservationStatus.UNCERTAIN in distinct:
        return ObservationStatus.UNCERTAIN, False
    return ObservationStatus.NOT_ASSESSED, False


#: The two statuses that can move a score. A disagreement *between* these two for
#: one term is the only duplicate this module refuses to resolve.
_INFORMATIVE_STATUSES: frozenset[ObservationStatus] = frozenset(
    {ObservationStatus.OBSERVED, ObservationStatus.EXCLUDED}
)
