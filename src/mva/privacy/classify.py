"""Path-based sensitivity classification.

Classification here is a *cheap prior*, not a proof. It answers "how should this
path be treated before anyone opens it", which is the only question you can ask
about a 60 GiB CRAM without paying to read it. The verification is always the
content scan (GP-43); see :mod:`mva.privacy.export`, where both must agree.

The default is :attr:`~mva.models.base.Sensitivity.SENSITIVE`. An unrecognised
path in a pipeline that handles patient genomes is far more likely to be an
unanticipated output than a safe one, and the cost asymmetry is total: a
misclassified-public file is one export away from being irreversible, while a
misclassified-sensitive file costs someone an explicit allowlist entry.
"""

from __future__ import annotations

from pathlib import Path
from typing import Final

from mva.models.base import Sensitivity

#: Extensions that are patient-derived by construction, or that routinely carry a
#: genotype matrix. Kept in sync with the corresponding sections of ``.gitignore``
#: — the audit's ``gitignore_effectiveness`` check is what proves they agree.
SENSITIVE_EXTENSIONS: Final[frozenset[str]] = frozenset(
    {
        # Variant calls
        ".vcf",
        ".gvcf",
        ".bcf",
        ".tbi",
        ".csi",
        ".idx",
        # Alignments
        ".bam",
        ".bai",
        ".cram",
        ".crai",
        ".sam",
        ".sai",
        # Raw reads
        ".fastq",
        ".fq",
        ".ubam",
        ".sra",
        ".bcl",
        # Reference / interval / array / pedigree
        ".fa",
        ".fasta",
        ".fai",
        ".2bit",
        ".dict",
        ".ped",
        ".fam",
        ".bim",
        ".bed",
        ".pgen",
        ".psam",
        ".pvar",
        ".cnv",
        ".seg",
        ".idat",
        ".gtc",
        # Clinical
        ".phenopacket",
        ".mrn",
        ".phi",
        # Derived analytical stores that can hold genotypes
        ".duckdb",
        ".db",
        ".sqlite",
        ".sqlite3",
        ".parquet",
        ".feather",
        ".arrow",
        ".h5",
        ".hdf5",
        ".npz",
        ".pkl",
        ".pickle",
    }
)

#: Compression suffixes that are stripped before extension matching, so that
#: ``proband.vcf.gz`` and ``sample_R1.fastq.gz`` classify like their plain forms.
_COMPRESSION_SUFFIXES: Final[frozenset[str]] = frozenset({".gz", ".bgz", ".bz2", ".zst", ".xz"})

#: Directory names that mean "patient workspace" wherever they appear in a path.
_SENSITIVE_DIR_NAMES: Final[frozenset[str]] = frozenset(
    {
        "workspace",
        "work",
        ".work",
        "private",
        "sensitive",
        "patient",
        "case_data",
        "real-case",
        "runs",
    }
)

#: Directory sequences whose contents are public by construction and are audited
#: as such (``synthetic_fixtures_marked`` for the fixtures, curation review for
#: the knowledge tables).
_PUBLIC_PREFIXES: Final[tuple[tuple[str, ...], ...]] = (
    ("knowledge", "public"),
    ("knowledge", "manifests"),
    ("tests", "fixtures", "synthetic"),
    ("tests", "golden"),
    ("docs",),
    ("templates",),
    ("config",),
    ("prompts",),
)


def _stripped_suffix(path: Path) -> str:
    """The meaningful extension, with any compression suffix removed."""
    suffixes = [s.lower() for s in path.suffixes]
    while suffixes and suffixes[-1] in _COMPRESSION_SUFFIXES:
        suffixes.pop()
    return suffixes[-1] if suffixes else ""


def is_sensitive_extension(path: Path) -> bool:
    """Whether the path's extension marks it as patient-derived.

    Handles multi-part suffixes: ``.vcf.gz``, ``.fastq.gz``, ``.g.vcf.gz`` and
    ``.phenopacket.json`` all resolve to their genomic component. A bare ``.gz``
    with no inner extension is *not* treated as sensitive by extension alone — the
    content scan's magic-byte sniffer is the control for that case.
    """
    suffix = _stripped_suffix(path)
    if suffix in SENSITIVE_EXTENSIONS:
        return True
    # `.phenopacket.json` and friends: the payload extension is the penultimate one.
    lowered = [s.lower() for s in path.suffixes]
    return any(s in SENSITIVE_EXTENSIONS for s in lowered[:-1]) if len(lowered) > 1 else False


def _has_prefix(parts: tuple[str, ...], prefix: tuple[str, ...]) -> bool:
    span = len(prefix)
    return any(parts[i : i + span] == prefix for i in range(len(parts) - span + 1))


def classify_path(path: Path) -> Sensitivity:
    """Classify a path without reading it.

    Order matters and is deliberate:

    1. A sensitive extension wins outright, wherever the file sits. A ``.vcf``
       under ``docs/`` is a mistake, not a document.
    2. A sensitive directory component wins next, so every untyped artifact in a
       run directory inherits the run's sensitivity.
    3. Only then can an audited public prefix grant ``PUBLIC``.
    4. Anything else is ``SENSITIVE``. There is no "probably fine" tier.

    ``DERIVED_SAFE`` is never returned here: it is a claim about *content*
    (aggregate counts, no genotypes) that no path inspection can support. It is
    assigned by the stage that produces the artifact and re-checked at export.
    """
    if is_sensitive_extension(path):
        return Sensitivity.SENSITIVE

    parts = tuple(path.parts)
    if any(part in _SENSITIVE_DIR_NAMES for part in parts):
        return Sensitivity.SENSITIVE

    if any(_has_prefix(parts, prefix) for prefix in _PUBLIC_PREFIXES):
        return Sensitivity.PUBLIC

    return Sensitivity.SENSITIVE
