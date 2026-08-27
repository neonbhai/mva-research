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


class AdapterUnavailableError(MvaError):
    """A required external tool or dataset is not present."""


class EvidenceError(MvaError):
    """Evidence-store integrity violation."""


class UnsourcedAssertionError(EvidenceError):
    """A report tried to state a claim with no resolvable evidence (GP-10)."""


class PrivacyViolationError(MvaError):
    """A privacy boundary would be crossed. Always fail closed."""


class NetworkDeniedError(PrivacyViolationError):
    """Outbound network was attempted while the offline profile was armed."""


class ExportBlockedError(PrivacyViolationError):
    """An artifact failed the public-export gate."""


class DeterminismError(MvaError):
    """A repeat run produced different bytes (GP-30)."""
