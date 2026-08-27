"""Offline enforcement for the patient-data path — a tripwire, not a boundary.

What this actually is
---------------------
A ``sys.addaudithook`` that raises :class:`~mva.errors.NetworkDeniedError` when
CPython emits one of a small set of network audit events while the profile is
armed. It catches the realistic failure mode: someone adds ``requests.get(...)``
to an annotation step, or a library silently phones home for a reference, and the
run stops loudly instead of shipping a proband's coordinates to a third party.

**Honest limits. Read these before relying on it.**

* **C extensions bypass it entirely.** ``pysam``/``htslib`` and ``cyvcf2`` call
  ``connect(2)``, ``getaddrinfo(3)`` and libcurl from C. No Python audit event is
  emitted, so nothing here fires. This is not a corner case — it is the primary
  genomics I/O path in this project.
* **Subprocesses are unpoliced.** Audit hooks are per-interpreter. Once
  ``bcftools``/``samtools``/``git`` is spawned it has an unrestricted network.
  ``strict=True`` blocks the *spawn*, which is a different and much blunter
  control.
* **``ctypes`` bypasses everything.** ``ctypes.CDLL("libc").connect(...)`` reaches
  the kernel without touching the socket module.
* **Audit hooks cannot be removed.** ``sys.addaudithook`` is append-only by design
  (a removable hook would be a trivial bypass). Hence the module-level ``_armed``
  flag: we gate inside a permanently installed hook rather than pretending to
  uninstall one.
* **In-process bypass is trivial for hostile code.** Any code that can set
  ``mva.privacy.netguard._armed = False`` is past the guard. This defends against
  mistakes, not adversaries.
* **Only the events in :data:`BLOCKED_EVENTS` are seen.** CPython emits audit
  events for a specific, finite list of socket operations, and anything not on
  that list is invisible here regardless of what it does. Raw ``socket.send`` on
  an already-connected socket, ``sendmsg``, ``os.write`` to a socket file
  descriptor and an ``AF_UNIX`` hop to a local relay all reach the network
  without emitting anything this hook can block. ``socket.connect``,
  ``socket.sendto`` and ``socket.bind`` cover how a datagram or a stream actually
  gets started from Python, which is where a mistake shows up; they are not a
  proof that nothing left.

The real boundary is at the OS: ``sandbox-exec``/Seatbelt with a no-network
profile, a ``pf`` deny rule for the run user, a network namespace on Linux, or —
the control that never has a bug — the Wi-Fi switch. ``NetworkProfile.OFFLINE_ENFORCED``
in :mod:`mva.config` means "this hook is armed **and** an OS control is asserted";
``OFFLINE_BEST_EFFORT`` means only this hook, and the run manifest records which.

The CRAM reference trap
-----------------------
htslib, when it opens a CRAM whose reference sequences are not local, fetches them
from ``https://www.ebi.ac.uk/ena/cram/md5/<md5>`` — by default, silently, with no
Python-visible network call. The MD5s it sends are of reference contigs rather
than patient sequence, but the request itself discloses that a specific assembly
is being processed from this IP at this time, and the code path is one htslib
version away from carrying more. Entering the profile therefore also sets
``REF_PATH=/dev/null`` (never look anything up remotely) and points ``REF_CACHE``
at the workspace. This is an environment control, not a hook, and it is the only
thing in this module that works against C code.
"""

from __future__ import annotations

import contextlib
import os
import sys
import tempfile
from pathlib import Path
from types import TracebackType
from typing import Any, Final, Self

from mva.config import resolve_workspace
from mva.errors import ConfigError, NetworkDeniedError

#: Audit events that constitute outbound network from Python.
#:
#: ``socket.sendto`` is the one that matters most and was missing. A connectionless
#: UDP socket never emits ``socket.connect``: ``sendto()`` puts a datagram on the
#: wire in a single call, so a few lines of pure Python — no library, no C
#: extension — could carry a phenotype profile straight out of an armed profile
#: with nothing firing. ``socket.bind`` covers the listening half (a reverse
#: channel is still a channel) and ``socket.gethostbyaddr`` the reverse-DNS lookup
#: that ``getaddrinfo`` blocking left open.
BLOCKED_EVENTS: Final[frozenset[str]] = frozenset(
    {
        "socket.connect",
        "socket.sendto",
        "socket.bind",
        "socket.getaddrinfo",
        "socket.gethostbyname",
        "socket.gethostbyaddr",
        "urllib.Request",
    }
)

#: Additional events blocked under ``strict=True``. These are the documented
#: bypass routes, not network operations in themselves — blocking them is coarse
#: and will break any stage that shells out (the privacy audit itself shells out
#: to ``git``, so it must never run inside a strict profile).
STRICT_EXTRA_EVENTS: Final[frozenset[str]] = frozenset(
    {
        "subprocess.Popen",
        "os.system",
        "os.exec",
        "os.posix_spawn",
        "ctypes.dlopen",
    }
)

_armed: bool = False
_strict_mode: bool = False
_hook_installed_flag: bool = False


def _audit_hook(event: str, args: tuple[Any, ...]) -> None:
    """The installed hook. Inert until a profile arms it.

    ``args`` is accepted because the audit-hook protocol requires it and is then
    immediately discarded: for ``socket.connect`` it is the destination address,
    for ``urllib.Request`` the full URL. Including any of it in the raised error
    would push a hostname (and, for a leaking call, potentially a query string
    built from patient data) into the traceback, the log and the model context —
    exactly the GP-41 failure this package exists to prevent.
    """
    del args
    if not _armed:
        return
    if event in BLOCKED_EVENTS or (_strict_mode and event in STRICT_EXTRA_EVENTS):
        msg = (
            f"{event}: blocked by the armed offline profile. Event arguments are "
            "withheld under GP-41. The patient-data path must not make outbound "
            "connections; pre-download public reference data in a separate, "
            "explicitly online acquisition step."
        )
        raise NetworkDeniedError(msg)


def arm_audit_hook() -> None:
    """Install the audit hook exactly once for this interpreter.

    Named for its purpose rather than its effect: it makes arming *possible*. The
    hook stays inert until an :class:`OfflineProfile` sets ``_armed``. Installation
    is unconditional and permanent because ``sys.addaudithook`` has no removal
    counterpart — so the flag, not the hook, is what a context manager toggles.
    """
    global _hook_installed_flag  # noqa: PLW0603 - audit hooks are permanent; see docstring
    if _hook_installed_flag:
        return
    sys.addaudithook(_audit_hook)
    _hook_installed_flag = True


def is_armed() -> bool:
    """Whether outbound network is currently being denied at the Python level."""
    return _armed


def hook_installed() -> bool:
    """Whether the audit hook has been installed in this interpreter."""
    return _hook_installed_flag


def _ref_cache_root() -> Path:
    """Where htslib may cache reference sequences: inside the workspace if we have one."""
    try:
        return resolve_workspace().root
    except (ConfigError, OSError):
        # No workspace configured yet (WorkspaceError subclasses ConfigError). A
        # process-local temp cache still keeps htslib away from the EBI endpoint.
        return Path(tempfile.gettempdir()) / "mva-ref-cache"


def reference_cache_env(workspace_root: Path) -> dict[str, str]:
    """The htslib environment that keeps CRAM reference resolution local.

    ``REF_PATH=/dev/null`` is the important half. htslib's *default* REF_PATH ends
    in the EBI URL template, so leaving it unset means "look locally, then ask the
    internet"; setting it to a path that resolves nothing means "look locally, then
    give up", which is the behaviour we want on the patient-data path.

    ``REF_CACHE`` uses htslib's ``%2s/%2s/%s`` expansion over the MD5 hex digest, so
    any sequence that *is* resolved locally is cached inside the workspace rather
    than in ``~/.cache``, where it would sit outside the privacy boundary.
    """
    cache = workspace_root / "ref_cache"
    # A read-only cache root is a degradation, not a privacy failure: REF_PATH alone
    # already suppresses the remote lookup.
    with contextlib.suppress(OSError):
        cache.mkdir(parents=True, exist_ok=True)
    return {"REF_PATH": "/dev/null", "REF_CACHE": f"{cache.as_posix()}/%2s/%2s/%s"}


def configure_reference_cache(workspace_root: Path) -> None:
    """Apply :func:`reference_cache_env` to this process, permanently.

    For the composition root, which sets it once at startup and never unwinds it.
    :class:`OfflineProfile` applies the same variables but restores the previous
    values on exit, because a context manager that silently mutated the process
    environment for good would be a trap.

    This is the only control in this module that constrains C code, and it is
    therefore the only one that constrains ``pysam``/``htslib`` — the audit hook
    does not see anything libhts does.
    """
    os.environ.update(reference_cache_env(workspace_root))


class OfflineProfile:
    """Context manager that arms the hook and neutralises htslib's remote fetch.

    Re-entrant and nesting-safe: the previous armed/strict state is restored on
    exit, so an inner profile cannot silently disarm an outer one.

    ``strict=True`` additionally blocks process spawning and ``ctypes.dlopen``.
    Use it only around pure-Python analysis; it will break any stage that invokes
    an external binary.
    """

    def __init__(self, *, strict: bool = False) -> None:
        self._strict = strict
        self._previous_armed = False
        self._previous_strict = False
        self._previous_env: dict[str, str | None] = {}

    def _set_env(self, key: str, value: str) -> None:
        self._previous_env[key] = os.environ.get(key)
        os.environ[key] = value

    def __enter__(self) -> Self:
        global _armed, _strict_mode  # noqa: PLW0603 - the module flag IS the arming mechanism
        arm_audit_hook()
        self._previous_armed = _armed
        self._previous_strict = _strict_mode

        for key, value in reference_cache_env(_ref_cache_root()).items():
            self._set_env(key, value)

        _armed = True
        _strict_mode = self._strict
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        global _armed, _strict_mode  # noqa: PLW0603 - the module flag IS the arming mechanism
        _armed = self._previous_armed
        _strict_mode = self._previous_strict
        for key, value in self._previous_env.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
        self._previous_env.clear()


# Installed at import: the hook must exist before any stage runs, and it costs
# nothing while disarmed.
arm_audit_hook()
