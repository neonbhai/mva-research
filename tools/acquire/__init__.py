"""The public-reference acquisition tool.

``src/mva/annotation`` (and every other module on the patient-data path) is
structurally forbidden from importing a network client -- see
``tests/unit/test_architecture.py::test_no_network_clients_in_sensitive_stages``.
This package is the other half of that design: the one place in the project
allowed to fetch reference data over the network, kept deliberately outside
``src/mva`` so it can never be imported by (or accidentally entangled with) code
that also sees patient data.

See ``knowledge/adapters/README.md``, "Why no adapter here may touch the network
(PRIV-05)", for the full contract this package implements:

* it runs before, and separately from, any patient data being loaded;
* it sends only public reference identifiers -- a dataset name, a release tag, a
  URL -- to a small, closed allowlist of public reference hosts
  (:mod:`tools.acquire.hosts`), and refuses (loudly) to fetch anything else;
* it writes hash-pinned artifacts to an external resource root that is asserted to
  be outside this repository (:mod:`tools.acquire.fetch`), and records what it
  fetched -- and what it did not -- in a committed, deterministic manifest
  (:mod:`tools.acquire.manifest`) that ``src/mva/annotation`` adapters can later
  read and verify offline.

Run it via ``uv run python -m tools.acquire <command>`` (see :mod:`tools.acquire.cli`).
"""

from __future__ import annotations

from tools.acquire.catalog import KNOWN_RESOURCES
from tools.acquire.errors import (
    AcquisitionError,
    DisallowedHostError,
    ResourceFetchError,
    ResourceRootError,
    ResourceVerificationError,
)
from tools.acquire.fetch import (
    DEFAULT_RESOURCE_ROOT,
    fetch_resource,
    is_download_stable,
    resolve_resource_root,
    sniff_content_mismatch,
)
from tools.acquire.hosts import ALLOWED_HOSTS, assert_allowed_host
from tools.acquire.manifest import (
    ResourceManifest,
    load_resources_manifest,
    render_resources_yaml,
    write_resources_manifest,
)
from tools.acquire.models import ResourceEntry, ResourceStatus
from tools.acquire.survey import survey_all, survey_resource
from tools.acquire.verify import (
    VerificationResult,
    VerificationStatus,
    assert_verified,
    verify_all,
    verify_resource,
)

__all__ = [
    "ALLOWED_HOSTS",
    "DEFAULT_RESOURCE_ROOT",
    "KNOWN_RESOURCES",
    "AcquisitionError",
    "DisallowedHostError",
    "ResourceEntry",
    "ResourceFetchError",
    "ResourceManifest",
    "ResourceRootError",
    "ResourceStatus",
    "ResourceVerificationError",
    "VerificationResult",
    "VerificationStatus",
    "assert_allowed_host",
    "assert_verified",
    "fetch_resource",
    "is_download_stable",
    "load_resources_manifest",
    "render_resources_yaml",
    "resolve_resource_root",
    "sniff_content_mismatch",
    "survey_all",
    "survey_resource",
    "verify_all",
    "verify_resource",
    "write_resources_manifest",
]
