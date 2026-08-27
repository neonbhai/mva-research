"""Adapter contracts for the annotation stage (GP-02, GP-20).

The annotation stage never talks to an annotation *tool*; it talks to these
Protocols. That indirection is what lets the repository ship a deterministic,
fully local, obviously-synthetic substitute today and swap in a real VEP /
gnomAD / ClinVar adapter later without touching a line of orchestration code.
See `knowledge/adapters/README.md` for the slot-in procedure.

Three rules are baked into the signatures rather than left to convention:

* **Typed models only (GP-02).** Every adapter returns `mva.models.variant`
  records. Parsing happens once, at the boundary, inside the adapter.
* **Absence is representable (GP-14).** The return type is a `Mapping` keyed by
  variant ID, and an adapter omits variants it knows nothing about. It never
  invents an empty-but-present answer, and it certainly never returns AF = 0 for
  a variant it has never seen. "Not in the table" and "in the table with value
  zero" are different facts and must stay different in the type.
* **Every adapter is signed (GP-20).** `name` and `version` are mandatory and
  flow into every `EvidenceItem` and into the run manifest, so a reader can
  always tell which tool at which version produced a claim — and, for the
  adapters shipped here, that the tool was a synthetic stand-in.

Batch (`Sequence[str]` in, `Mapping[str, ...]` out) rather than per-variant
calls: a real adapter is a bulk lookup or a subprocess, and a per-variant API
would push callers into N round-trips.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from typing import Protocol, runtime_checkable

from mva.models.variant import ClinicalAssertion, ConsequenceAnnotation, PopulationFrequency

#: Appended to the ``limitations`` of every EvidenceItem produced from a
#: synthetic-substitute adapter. GP-20: a mock is labelled as a mock everywhere it
#: is visible, not only in the maturity ledger.
SYNTHETIC_STANDIN_LIMITATION: str = (
    "Produced by a SYNTHETIC STAND-IN adapter reading a fabricated local table, not by "
    "a real annotation tool (VEP/SnpEff/gnomAD/ClinVar/SpliceAI). The values are invented "
    "for a fictional demo case, are NOT biologically valid, and must never be used for "
    "clinical interpretation or presented as real annotation."
)


class AdapterRole(StrEnum):
    """Which slot in the annotation stage an adapter fills."""

    CONSEQUENCE = "consequence"
    FREQUENCY = "frequency"
    CLINICAL = "clinical"


@runtime_checkable
class ConsequenceAdapter(Protocol):
    """Predicts molecular consequences, per transcript.

    Implementations must return **all** transcript annotations they hold for a
    variant, in a deterministic order. Collapsing to the canonical transcript is a
    data-loss bug, not an optimisation: a variant can be benign on MANE-Select and
    splice-disrupting on the tissue-relevant isoform.
    """

    @property
    def name(self) -> str:
        """Stable adapter identity; recorded on every EvidenceItem it justifies."""
        ...

    @property
    def version(self) -> str:
        """Version of the tool or table behind this adapter."""
        ...

    def annotate(
        self, variant_ids: Sequence[str]
    ) -> Mapping[str, tuple[ConsequenceAnnotation, ...]]:
        """Annotate canonical variant IDs.

        Keys appear only for variants the adapter actually knows about; a missing
        key means "not annotated", never "no consequence".
        """
        ...


@runtime_checkable
class FrequencyAdapter(Protocol):
    """Looks up population allele frequencies.

    A variant missing from the result has **no frequency data**. That is not the
    same as, and may never be rendered as, an allele frequency of zero (GP-14).
    """

    @property
    def name(self) -> str:
        """Stable adapter identity; recorded on every EvidenceItem it justifies."""
        ...

    @property
    def version(self) -> str:
        """Version of the dataset release behind this adapter."""
        ...

    def frequencies(
        self, variant_ids: Sequence[str]
    ) -> Mapping[str, tuple[PopulationFrequency, ...]]:
        """Return per-variant frequency observations.

        Each observation carries its own source, version and population (GP-18);
        those come from the underlying dataset and are never supplied by the caller.
        """
        ...


@runtime_checkable
class ClinicalAdapter(Protocol):
    """Looks up curated clinical-significance assertions (ClinVar-style)."""

    @property
    def name(self) -> str:
        """Stable adapter identity; recorded on every EvidenceItem it justifies."""
        ...

    @property
    def version(self) -> str:
        """Version of the curated release behind this adapter."""
        ...

    def assertions(self, variant_ids: Sequence[str]) -> Mapping[str, tuple[ClinicalAssertion, ...]]:
        """Return per-variant clinical assertions.

        An empty result means "this source holds nothing on record", which is not
        evidence of benignity.
        """
        ...


@runtime_checkable
class MaturityAware(Protocol):
    """Optional mixin protocol: an adapter that declares its own maturity grade."""

    @property
    def synthetic(self) -> bool:
        """True when this adapter fabricates data rather than reading a real source."""
        ...


def is_synthetic(adapter: object) -> bool:
    """Whether an adapter must be labelled as a mock (GP-20).

    **Fails closed.** An adapter that does not declare a ``synthetic`` property is
    treated as synthetic, so the disclosure appears by default and has to be
    switched off deliberately by a real adapter declaring ``synthetic = False``.
    The reverse default would let an undeclared mock quietly present itself as a
    real annotation source, which is the exact failure GP-20 exists to prevent.
    """
    if isinstance(adapter, MaturityAware):
        return adapter.synthetic
    return True


@dataclass(frozen=True)
class AdapterDescriptor:
    """Identity of one bound adapter: coverage key, warning subject, provenance row."""

    role: AdapterRole
    name: str
    version: str
    synthetic: bool

    @property
    def label(self) -> str:
        """``name@version``, the form used in warnings and report footers."""
        return f"{self.name}@{self.version}"


@dataclass(frozen=True)
class AdapterSet:
    """The adapters bound for one annotation run.

    ``clinical`` is optional and defaults to *absent* rather than to a silently
    empty stub, so that "no clinical source was configured" stays distinguishable
    from "the configured clinical source returned nothing".
    """

    consequence: ConsequenceAdapter
    frequency: FrequencyAdapter
    clinical: ClinicalAdapter | None = None

    def descriptors(self) -> tuple[AdapterDescriptor, ...]:
        """Bound adapters in a fixed role order, for deterministic reporting."""
        bound = [
            AdapterDescriptor(
                role=AdapterRole.CONSEQUENCE,
                name=self.consequence.name,
                version=self.consequence.version,
                synthetic=is_synthetic(self.consequence),
            ),
            AdapterDescriptor(
                role=AdapterRole.FREQUENCY,
                name=self.frequency.name,
                version=self.frequency.version,
                synthetic=is_synthetic(self.frequency),
            ),
        ]
        if self.clinical is not None:
            bound.append(
                AdapterDescriptor(
                    role=AdapterRole.CLINICAL,
                    name=self.clinical.name,
                    version=self.clinical.version,
                    synthetic=is_synthetic(self.clinical),
                )
            )
        return tuple(bound)
