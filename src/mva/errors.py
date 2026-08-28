"""Typed exceptions.

Exception messages are constructed from identifiers and counts only, never from
record text (PRIV-09). A traceback that escapes to a terminal enters the model's
context and any log aggregator; a genotype embedded in an exception message would
travel everywhere the traceback travels.
"""

from __future__ import annotations


class MvaError(Exception):
    """Base for every error this package raises."""


class ConfigError(MvaError):
    """Invalid or unsafe configuration."""


class WorkspaceError(ConfigError):
    """The configured workspace violates a privacy boundary."""


class GenomeBuildMismatchError(MvaError):
    """Two records from different assemblies were compared."""


class IngestionError(MvaError):
    """A source artifact could not be read or is malformed."""


class ReferenceMismatchError(IngestionError):
    """A record's REF allele disagrees with the reference genome."""


class ReferenceUnusableError(MvaError):
    """The reference genome could not supply a base the allele rule needed.

    Distinct from :class:`ReferenceMismatchError`, which is a statement *about the
    record* made by a reference that worked. This one says the reference itself is
    broken — it raised, it is missing the contig, or it returned something that is
    not a nucleotide — so nothing is known about the record at all.

    Also distinct from "there is genuinely no base at this position", which is the
    reference being complete and correct at the edge of a contig. Collapsing the
    two is what let a failed FASTA read degrade a run to trim-only join keys while
    the adapters went on reporting that left-alignment had been applied (ADR 0018).

    Raised only by the reference-consuming primitives in :mod:`mva.alleles` whose
    return type has nowhere to carry the degraded state.
    """


class AdapterUnavailableError(MvaError):
    """A required external tool or dataset is not present."""


class AnnotationError(MvaError):
    """The annotation stage was driven in a way it cannot honour.

    Raised for misuse of the stage rather than for adapter failure: re-iterating a
    single-pass stream, or asking for a batch size that cannot bound anything. An
    adapter that is missing raises :class:`AdapterUnavailableError` instead.
    """


class EvidenceError(MvaError):
    """Evidence-store integrity violation."""


class UnsourcedAssertionError(EvidenceError):
    """A report tried to state a claim with no resolvable evidence (GP-10)."""


class ReportCompletenessError(MvaError):
    """A renderer would have emitted fewer hypotheses than it was given.

    Report sections are chosen by predicate. A predicate set that does not
    partition the input does not shorten the report — it *erases* the rows it
    misses, silently, and the row most likely to be missed is a contraindicated
    compound in a state no section expects.
    """


class PrivacyViolationError(MvaError):
    """A privacy boundary would be crossed. Always fail closed."""


class NetworkDeniedError(PrivacyViolationError):
    """Outbound network was attempted while the offline profile was armed."""


class ExportBlockedError(PrivacyViolationError):
    """An artifact failed the public-export gate."""


class DeterminismError(MvaError):
    """A repeat run produced different bytes (GP-30)."""
