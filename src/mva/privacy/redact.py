"""Redaction of text and of the logging pipeline (GP-42).

Why a Filter and not a Formatter
--------------------------------
A ``logging.Formatter`` is the wrong instrument for this job, in four independent
ways:

1. **It runs too late and too narrowly.** A Formatter is consulted by exactly one
   handler, at emit time. Any handler configured without our formatter — a
   third-party library's ``StreamHandler``, a ``QueueHandler``, an APM/Sentry
   integration, ``logging.lastResort`` — emits the raw record. A Filter mutates
   the ``LogRecord`` *object*, so every consumer downstream of that point sees the
   scrubbed version, whether or not it uses our formatter.
2. **Structured handlers never call it.** JSON and OTLP handlers serialise
   ``record.msg`` and ``record.args`` directly. There is no format string to
   intercept.
3. **``Formatter.format`` caches.** It writes ``record.exc_text`` back onto the
   record; a second handler then reuses that cached, unscrubbed traceback.
4. **Tracebacks are not part of the message.** ``exc_info`` carries frame objects
   whose locals routinely include the parsed VCF row that caused the exception.
   Only record-level surgery removes them.

Why the filter goes on every HANDLER, not on the logger
-------------------------------------------------------
This is the trap that makes most "we redact our logs" claims false. ``Logger.filter``
is consulted only for records logged *directly to that logger*. A record created
by ``logging.getLogger("mva.ingestion.vcf")`` and propagated up to the root
logger's handlers is filtered by the child logger's filters and then handed
straight to each ancestor's *handlers* — ancestor **logger** filters are never
consulted (see ``Logger.callHandlers`` in CPython: it iterates ``c.handlers`` and
calls ``hdlr.handle(record)``, which applies handler filters only). So a filter on
the root logger protects nothing that any submodule logs.

Handler filters plus a ``LogRecordFactory`` give two independent layers: the
factory scrubs at construction (covering handlers we never see), and the handler
filters scrub again at emit (covering records built by other means, e.g.
``logging.makeLogRecord`` in a ``QueueListener``).
"""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sized
from typing import Any, Final, cast

from mva.privacy.patterns import REDACTION_RULES, placeholder

#: Stand-in for a suppressed traceback. The exception TYPE is retained because a
#: type name is a shape, not a value, and without it a redacted log is undebuggable.
_TRACEBACK_PLACEHOLDER: Final = "<REDACTED:traceback>"
_STACK_PLACEHOLDER: Final = "<REDACTED:stack_info>"

#: Scalars are passed through unchanged so that ``%d``/``%f`` format specifiers in
#: existing log statements keep working. Their textual form is still redacted when
#: it reaches ``record.msg``, and a bare integer carries no rule-detectable payload.
_PASSTHROUGH_TYPES: Final = (bool, int, float, type(None))

#: Attributes the stdlib itself puts on every ``LogRecord``. Everything on the
#: record that is NOT in this set arrived through ``extra=``, which is exactly the
#: channel structured/JSON/OTLP handlers serialise — and exactly the channel that
#: used to survive untouched, because ``Logger.makeRecord`` writes ``extra`` into
#: ``record.__dict__`` *after* the record factory has run.
_RESERVED_RECORD_ATTRS: Final[frozenset[str]] = frozenset(
    {
        "args",
        "asctime",
        "created",
        "exc_info",
        "exc_text",
        "filename",
        "funcName",
        "levelname",
        "levelno",
        "lineno",
        "message",
        "module",
        "msecs",
        "msg",
        "name",
        "pathname",
        "process",
        "processName",
        "relativeCreated",
        "stack_info",
        "taskName",
        "thread",
        "threadName",
    }
)


# ---------------------------------------------------------------------------
# Text redaction
# ---------------------------------------------------------------------------


def redact_bytes(data: bytes) -> bytes:
    """Replace every rule match with a fixed placeholder.

    Overlapping matches from different rules are merged into a single placeholder
    labelled with the first (lowest-offset) rule, because emitting two nested
    placeholders would leak the *structure* of the overlap — and, worse, could
    leave an unredacted sliver between two spans that were computed independently.
    """
    spans: list[tuple[int, int, str]] = []
    for rule in REDACTION_RULES:
        for match in rule.pattern.finditer(data):
            start, end = match.span()
            if end > start:
                spans.append((start, end, rule.rule_id))
    if not spans:
        return data

    spans.sort()
    merged: list[tuple[int, int, str]] = []
    cur_start, cur_end, cur_id = spans[0]
    for start, end, rule_id in spans[1:]:
        if start <= cur_end:
            cur_end = max(cur_end, end)
        else:
            merged.append((cur_start, cur_end, cur_id))
            cur_start, cur_end, cur_id = start, end, rule_id
    merged.append((cur_start, cur_end, cur_id))

    parts: list[bytes] = []
    previous = 0
    for start, end, rule_id in merged:
        parts.append(data[previous:start])
        parts.append(placeholder(rule_id, end - start).encode("ascii"))
        previous = end
    parts.append(data[previous:])
    return b"".join(parts)


def redact_text(text: str) -> str:
    """Apply the whole rule battery to a string.

    Detection happens in bytes (the rules are byte patterns), so the string is
    encoded with ``surrogateescape`` to survive lone surrogates that arrive from
    ``os`` APIs, and decoded back with ``errors="replace"`` because a merged span
    can legitimately cut a multi-byte character in half. Lossy round-tripping is
    acceptable here and lossless round-tripping is not: this function exists to
    destroy information.
    """
    return redact_bytes(text.encode("utf-8", errors="surrogateescape")).decode(
        "utf-8", errors="replace"
    )


def safe_repr(obj: object) -> str:
    """Describe an object by type and shape, never by contents.

    Intended for exception messages and debug logs on the patient-data path, where
    ``repr(variant_row)`` would print the genotype. What survives is what a
    developer actually needs at 2am — "it was a list of 4310 things, not None".
    """
    if obj is None:
        return "None"
    type_name = type(obj).__name__
    shape = getattr(obj, "shape", None)
    if isinstance(shape, tuple):
        return f"<{type_name} shape={shape!r}>"
    if isinstance(obj, Mapping):
        return f"<{type_name} n_keys={len(obj)}>"
    if isinstance(obj, str | bytes | bytearray):
        return f"<{type_name} len={len(obj)}>"
    if isinstance(obj, Sized):
        return f"<{type_name} len={len(obj)}>"
    return f"<{type_name}>"


def _scrub_arg(value: object) -> object:
    if isinstance(value, _PASSTHROUGH_TYPES):
        return value
    if isinstance(value, str):
        return redact_text(value)
    return safe_repr(value)


def scrub_record(record: logging.LogRecord) -> None:
    """Rewrite a record in place so no rule-detectable content remains.

    Idempotent: re-running it over an already-scrubbed record is a no-op, which
    matters because the record factory and the handler filters both call it.
    """
    if isinstance(record.msg, str):
        record.msg = redact_text(record.msg)
    elif isinstance(record.msg, _PASSTHROUGH_TYPES):
        record.msg = str(record.msg)
    else:
        record.msg = safe_repr(record.msg)

    args = record.args
    if isinstance(args, Mapping):
        record.args = {key: _scrub_arg(value) for key, value in args.items()}
    elif isinstance(args, tuple):
        record.args = tuple(_scrub_arg(value) for value in args)

    if record.exc_info is not None:
        exc_type = record.exc_info[0]
        name = exc_type.__name__ if exc_type is not None else "Exception"
        # exc_info must stay a valid tuple-or-None for the stdlib, so the placeholder
        # goes into exc_text: Formatter.formatException is skipped whenever exc_text
        # is already set, which is precisely the interception point we want.
        record.exc_info = None
        record.exc_text = f"{_TRACEBACK_PLACEHOLDER}:{name}"
    elif record.exc_text:
        record.exc_text = _TRACEBACK_PLACEHOLDER

    if record.stack_info:
        record.stack_info = _STACK_PLACEHOLDER

    _scrub_extras(record)


def _scrub_extras(record: logging.LogRecord) -> None:
    """Scrub every ``extra=`` key the caller attached to the record.

    ``Logger.makeRecord`` copies ``extra`` into ``record.__dict__`` after the
    record factory returns, so the factory can never see it. A JSON, OTLP or
    logfmt handler serialises precisely those keys. Skipping them meant the
    module docstring's claim about structured handlers was false.

    Keys are left alone and values are scrubbed with the same ``_scrub_arg`` used
    for ``%``-args, so an ``extra={"n_variants": 4310}`` counter still formats as
    an integer while an ``extra={"row": <VcfRecord>}`` becomes a shape.
    """
    for key, value in list(record.__dict__.items()):
        if key in _RESERVED_RECORD_ATTRS or key.startswith("_"):
            continue
        record.__dict__[key] = _scrub_arg(value)


class GenomicRedactionFilter(logging.Filter):
    """A ``logging.Filter`` that never filters — it launders.

    ``filter()`` always returns ``True``; the value is the in-place mutation of the
    record. Dropping records instead would hide operational failures, and a
    privacy control that makes the system less observable gets removed.
    """

    def filter(self, record: logging.LogRecord) -> bool:
        scrub_record(record)
        return True


# ---------------------------------------------------------------------------
# Installation
# ---------------------------------------------------------------------------

_factory_installed = False


def _install_record_factory() -> None:
    global _factory_installed  # noqa: PLW0603 - the record factory is process-global
    if _factory_installed:
        return
    previous = logging.getLogRecordFactory()

    def _redacting_factory(*args: Any, **kwargs: Any) -> logging.LogRecord:
        record = previous(*args, **kwargs)
        scrub_record(record)
        return record

    logging.setLogRecordFactory(_redacting_factory)
    _factory_installed = True


def _handlers_for(logger: logging.Logger | None) -> list[logging.Handler]:
    """Every handler a record could reach.

    With no logger given this walks the entire manager registry, because a library
    that configured its own handler on its own logger is exactly the leak we are
    trying to close. With a logger given it walks that logger and its ancestors,
    which is the set ``Logger.callHandlers`` will visit.
    """
    handlers: list[logging.Handler] = []
    seen: set[int] = set()

    def _collect(candidate: logging.Logger) -> None:
        for handler in candidate.handlers:
            if id(handler) not in seen:
                seen.add(id(handler))
                handlers.append(handler)

    if logger is None:
        _collect(logging.getLogger())
        for existing in list(logging.getLogger().manager.loggerDict.values()):
            if isinstance(existing, logging.Logger):
                _collect(existing)
    else:
        current: logging.Logger | None = logger
        while current is not None:
            _collect(current)
            if not current.propagate:
                break
            current = current.parent
    return handlers


_handler_hook_installed = False


def _install_handler_hook() -> None:
    """Cover handlers that do not exist yet.

    Two independent holes are closed here, and both were live:

    1. ``install_redaction()`` armed the handlers present *at that moment*. A
       library that adds its own handler afterwards — and libraries do, on
       import, on first use, on ``basicConfig`` — got an unfiltered path
       straight out of the process. Wrapping ``Handler.__init__`` means a
       handler is born filtered.
    2. A handler constructed *before* the hook (or one whose ``__init__`` never
       chains to ``logging.Handler.__init__``) still has no filter. Wrapping
       ``Handler.handle`` scrubs at the last point every record must pass
       through, whatever filters the handler ended up with. ``handle`` is also
       the entry point used by records that never met the record factory —
       ``logging.makeLogRecord`` in a ``QueueListener``, or a record replayed by
       a test harness.

    Both patches are idempotent and process-global, matching the audit hook and
    the record factory: there is no correct way to make a redaction control
    removable, because removability is itself the bypass.
    """
    global _handler_hook_installed  # noqa: PLW0603 - process-global, like the record factory
    if _handler_hook_installed:
        return

    original_init = logging.Handler.__init__
    original_handle = logging.Handler.handle

    def _redacting_init(self: logging.Handler, *args: Any, **kwargs: Any) -> None:
        original_init(self, *args, **kwargs)
        if not any(isinstance(existing, GenomicRedactionFilter) for existing in self.filters):
            self.addFilter(GenomicRedactionFilter())

    def _redacting_handle(self: logging.Handler, record: logging.LogRecord) -> bool:
        scrub_record(record)
        return original_handle(self, record)

    cast(Any, logging.Handler).__init__ = _redacting_init
    cast(Any, logging.Handler).handle = _redacting_handle
    _handler_hook_installed = True


def redaction_installed() -> bool:
    """Whether GP-42 is armed in this interpreter.

    Exists so the audit's live probe can assert the state it found *before*
    arming anything. A probe that installs the control it is measuring can only
    ever report success.
    """
    return _factory_installed and _handler_hook_installed


def unfiltered_handlers() -> list[str]:
    """Type names of handlers, anywhere in the registry, carrying no redaction filter.

    Walks ``Logger.manager.loggerDict`` rather than the root logger's handlers:
    the handler a library attached to its own logger is never consulted by the
    root, so looking only at the root reports "clean" for exactly the
    configuration this control exists to catch.
    """
    return [
        type(handler).__name__
        for handler in _handlers_for(None)
        if not any(isinstance(f, GenomicRedactionFilter) for f in handler.filters)
    ]


def install_redaction(logger: logging.Logger | None = None) -> None:
    """Arm GP-42: scrub at record construction, at every handler, and at emit.

    Safe to call repeatedly (startup, then again after a library reconfigures
    logging); duplicate filters, factories and hooks are all suppressed.
    """
    _install_record_factory()
    _install_handler_hook()
    for handler in _handlers_for(logger):
        if not any(isinstance(existing, GenomicRedactionFilter) for existing in handler.filters):
            handler.addFilter(GenomicRedactionFilter())
