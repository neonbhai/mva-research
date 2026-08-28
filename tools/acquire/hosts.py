"""The public-reference host allowlist (PRIV-05 boundary, acquisition side).

``src/mva/annotation`` (and the rest of the patient-data path) is structurally
forbidden from importing a network client at all -- see
``tests/unit/test_architecture.py::test_no_network_clients_in_sensitive_stages``
and the "Why no adapter here may touch the network" section of
``knowledge/adapters/README.md``. This tool is the other half of that design: it
is the *only* place in the project allowed to open a socket, and in exchange it
carries a structural constraint of its own -- it may only ever talk to a small,
named set of public reference hosts.

This is not a courtesy check. It is the difference between "this tool can only
ever download ClinVar/gnomAD/HPO/etc." and "this tool can be pointed at
literally any URL", and the latter is one typo away from becoming a general
network client with patient-adjacent code nearby. A closed allowlist makes that
mistake a loud, immediate exception instead of a silent request.
"""

from __future__ import annotations

from typing import Final
from urllib.parse import urlsplit

from tools.acquire.errors import DisallowedHostError

#: The only hosts this tool will ever open a connection to. Every host here
#: serves a named, versioned, public reference dataset used by this project
#: (ClinVar, gnomAD, HPO, Gene2Phenotype/DDG2P, ClinGen gene validity, Ensembl,
#: MANE, OBO Foundry). Adding a host here is a deliberate, reviewable decision --
#: never widen this set to make a one-off fetch succeed.
ALLOWED_HOSTS: Final[frozenset[str]] = frozenset(
    {
        "ftp.ncbi.nlm.nih.gov",
        "storage.googleapis.com",
        "github.com",
        "objects.githubusercontent.com",
        "ftp.ensembl.org",
        "purl.obolibrary.org",
        "www.ebi.ac.uk",
        "ftp.ebi.ac.uk",
        "search.clinicalgenome.org",
        # Added deliberately, 2026-08-28, for the SnpEff consequence adapter. SnpEff
        # publishes its core jar and its per-genome databases only from its own
        # bucket -- there is no mirror on any host already listed here -- and the
        # GRCh38.115 database is what gene assignment runs against, which is the
        # difference between the pipeline emitting candidate pairs and emitting
        # none. It serves fixed, named archives over https and is content-pinned in
        # knowledge/manifests/resources.yaml, so a silent republish under an
        # unchanged name (which SnpEff does do: the core archive is literally called
        # `snpEff_latest_core.zip`) fails the manifest check rather than changing
        # transcripts underneath an existing report.
        "snpeff-public.s3.amazonaws.com",
    }
)


def assert_allowed_host(url: str) -> None:
    """Raise :class:`DisallowedHostError` unless ``url`` is https and its host is allowed.

    Called before any connection is opened -- both when a :class:`ResourceEntry`
    is constructed (so a bad URL can never even enter the registry) and again
    immediately before a fetch (so nothing can reach this point by constructing
    the object some other way). Plain ``http`` is refused outright: every host on
    this list serves https, and accepting cleartext for "just this one resource"
    is exactly the kind of quiet downgrade this allowlist exists to prevent.
    """
    parsed = urlsplit(url)
    if parsed.scheme != "https":
        msg = (
            f"Refusing to fetch {url!r}: only https URLs are permitted for public "
            f"reference acquisition (got scheme {parsed.scheme!r})."
        )
        raise DisallowedHostError(msg)

    host = (parsed.hostname or "").lower()
    if host not in ALLOWED_HOSTS:
        msg = (
            f"Refusing to fetch from host {host!r}: not on the public-reference "
            f"allowlist ({', '.join(sorted(ALLOWED_HOSTS))}). This tool fetches named, "
            "versioned public reference datasets only. If this is a genuine new public "
            "source, add it to ALLOWED_HOSTS deliberately, with a reason -- never widen "
            "this set to make one fetch succeed, and never point this tool at anything "
            "that could echo back patient data."
        )
        raise DisallowedHostError(msg)
