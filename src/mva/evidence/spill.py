"""Disk-backed overflow for :class:`~mva.evidence.ledger.EvidenceLedger`.

The ledger is an in-memory dict, and that is the right shape for a case: a few
thousand claims, all of them cheap to hold. It stops being the right shape at
whole-genome scale. Ingestion emits one analytical `EvidenceItem` per record and
annotation emits at least one more (the GP-14 "no frequency data" item) plus one
per transcript and per population, so a 4.5 M-record callset produces on the
order of 10-20 M items. At the ~1 KB a hydrated `EvidenceItem` costs that is tens
of gigabytes, and the process dies somewhere in the annotate stage.

This module is the overflow. When the ledger crosses its spill threshold it
migrates into a SQLite file and keeps writing there. **Nothing observable
changes**: the same items come back, in the same total order, and a repeat run
produces the same bytes (GP-30).

Why SQLite rather than the DuckDB evidence store
------------------------------------------------
The store is the *durable, queryable* copy of a finished run, with a schema, a
Parquet export and a byte-identity guarantee. This is a scratch buffer for a run
in progress: it is created, written once, read once and deleted. Making the store
serve both roles would mean writing 20 M rows into the artifact that
`verify_determinism` hashes, purely to have somewhere to put them. SQLite is in
the standard library, needs no dependency, and gives exactly the two operations
required — an ordered scan and a primary-key lookup.

Size, and why the documents are compressed
------------------------------------------
Measured on 762,705 annotation evidence items (``tools/scale/stage_harness.py
ledger``): the canonical-JSON form of one item is **1,640 bytes**, of which
``limitations`` alone is **624 bytes — 38%** — and ``method`` a further 93. Those
two fields are byte-identical across millions of items, because they are module
constants: the same four sentences about what a frequency lookup does not
establish, written out once per variant.

Stored verbatim in SQLite that came to **4,933 bytes per item on disk** (3.0x the
payload: rows that large overflow the page in a ``WITHOUT ROWID`` table, and each
secondary index re-carries the text primary key) and wrote at 4,726 items/s. At
whole-genome scale that is 34 GB and 24 minutes for a scratch file.

So documents are deflated against a **shared dictionary** built from the first
flush's own documents and stored in the file's ``meta`` table. Compressing each
row independently against 4 KiB of the boilerplate takes it to **236 bytes, 6.9x
smaller**, at 32,400 items/s — about what encoding the JSON costs anyway. The
dictionary is derived from the data rather than hard-coded so this module does
not have to know what any other stage's limitation text says, and it is stored in
the file so a row can always be read back by whatever wrote it.

Ordering
--------
The scan order is ``ORDER BY subject_kind, subject_id, category, direction,
evidence_id``, which is the column-by-column form of ``ledger._sort_key``. SQLite
compares ``TEXT`` with the ``BINARY`` collation, which is a byte-wise comparison
of the UTF-8 encoding; UTF-8 is order-preserving with respect to code points, so
that is the same order Python's ``sorted`` produces. A test asserts it on
non-ASCII subjects rather than leaving the reader to trust the argument.

Privacy
-------
Evidence subjects are variant IDs, so the spill file carries build-qualified
proband coordinates and is SENSITIVE. It must be created inside the external
workspace (``Workspace.tmp_dir``), never in the repository, and
:meth:`SqliteEvidenceSpill.close` unlinks it. A run that dies before close leaves
it behind, which is the same exposure as the run's own artifacts and is covered
by the same workspace boundary (GP-40, ADR 0006).
"""

from __future__ import annotations

import hashlib
import sqlite3
import zlib
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Final

from mva.errors import EvidenceError
from mva.models import EvidenceDirection, EvidenceItem

#: Items buffered in memory before a round trip to SQLite. Sized so the integrity
#: check (one ``SELECT`` per flush, not one per item) stays a single query while
#: the buffer itself stays small: 10,000 hydrated items is roughly 35 MB.
DEFAULT_FLUSH_BATCH: Final[int] = 10_000

#: Bytes of shared deflate dictionary. Measured sensitivity on real evidence
#: documents, ratio / throughput: 1 KiB 3.46x @ 32,400/s; 2 KiB 4.91x @ 33,900/s;
#: **4 KiB 6.94x @ 32,400/s**; 8 KiB 7.00x @ 26,100/s; 32 KiB 13.95x @ 13,400/s.
#: 4 KiB is the knee -- 8 KiB buys 1% more compression for 20% less throughput,
#: and 32 KiB doubles the ratio at less than half the speed. Priming the deflate
#: state is what costs, so a bigger dictionary is not free.
DICTIONARY_BYTES: Final[int] = 4096

#: Documents sampled to build that dictionary. zlib matches against the TAIL of
#: the dictionary, so the sample is truncated from the right.
_DICTIONARY_SAMPLE: Final[int] = 64

#: Raw deflate, no zlib header -- the header carries a dictionary checksum this
#: file does not need, since the dictionary travels in the same file.
_WBITS: Final[int] = -15

#: Deflate level. 1 rather than 6: with a shared dictionary the ratio difference
#: is small and the throughput difference is not.
_LEVEL: Final[int] = 1

#: Column list, in the order ``ledger._sort_key`` compares them. The composite
#: index below is declared on exactly this prefix, so the ordered scan reads the
#: index rather than sorting a temporary table.
_ORDER_COLUMNS: Final[str] = "subject_kind, subject_id, category, direction, evidence_id"

_SCHEMA: Final[str] = """
CREATE TABLE IF NOT EXISTS meta (
    key   TEXT PRIMARY KEY,
    value BLOB NOT NULL
);
CREATE TABLE IF NOT EXISTS ledger (
    evidence_id  TEXT PRIMARY KEY,
    subject_kind TEXT NOT NULL,
    subject_id   TEXT NOT NULL,
    category     TEXT NOT NULL,
    direction    TEXT NOT NULL,
    content_sha  TEXT NOT NULL,
    document     BLOB NOT NULL
);
CREATE INDEX IF NOT EXISTS ledger_order ON ledger (subject_kind, subject_id, category, direction, evidence_id);
CREATE INDEX IF NOT EXISTS ledger_subject ON ledger (subject_id);
"""  # noqa: E501 - one DDL statement per line reads better than a wrapped index definition


def encode(item: EvidenceItem) -> str:
    """The stored form of an item.

    Pydantic's own serialiser, **not** :func:`~mva.determinism.canonical_json`,
    and that is deliberate. Canonical JSON exists so that artifact bytes are
    reproducible (GP-30); this file is not an artifact — it is scratch, it is
    never hashed, and it is unlinked when the run ends. What it needs is an exact
    round trip, which a test asserts. Measured on real evidence items,
    ``model_dump_json`` encodes at 69,500/s against ``canonical_json``'s 37,400
    and decodes at 35,700 against 20,300, for the same 1,638 bytes -- which is
    ~90 seconds off a whole-genome ledger write, for no change to anything
    anyone can observe.

    Determinism is unaffected because nothing downstream sees this encoding:
    :meth:`~mva.evidence.ledger.EvidenceLedger.iter_items` hands back hydrated
    models, and :meth:`~mva.evidence.store.EvidenceStore.write_evidence_stream`
    re-serialises the payload canonically on its way into DuckDB.
    """
    return item.model_dump_json()


def decode(document: str) -> EvidenceItem:
    """Rehydrate a stored item. Round-trips exactly; a test holds it to that."""
    return EvidenceItem.model_validate_json(document)


def build_dictionary(documents: Sequence[str]) -> bytes:
    """A shared deflate dictionary derived from real documents.

    Truncated from the RIGHT because zlib matches against the tail of the
    dictionary, so the material nearest the end is the material that pays.
    """
    blob = "".join(documents[:_DICTIONARY_SAMPLE]).encode("utf-8")
    return blob[-DICTIONARY_BYTES:]


def compress(document: str, dictionary: bytes) -> bytes:
    compressor = zlib.compressobj(_LEVEL, zlib.DEFLATED, _WBITS, zdict=dictionary)
    return compressor.compress(document.encode("utf-8")) + compressor.flush()


def decompress(blob: bytes, dictionary: bytes) -> str:
    decompressor = zlib.decompressobj(_WBITS, zdict=dictionary)
    return (decompressor.decompress(blob) + decompressor.flush()).decode("utf-8")


def _fingerprint(document: str) -> str:
    """A short content fingerprint, for the ID-reuse check.

    64 bits, not 256. The check is "did the same content-derived ID arrive with
    different content", and over the ~10^7 items a whole-genome run produces the
    chance of a 64-bit collision hiding a real one is about 10^-6 -- against a
    stored cost of 16 bytes per row instead of 64, on rows that are now ~236
    bytes. blake2b at this size is also faster than sha256.
    """
    return hashlib.blake2b(
        document.encode("utf-8"), digest_size=8, usedforsecurity=False
    ).hexdigest()


class SqliteEvidenceSpill:
    """A single-run, write-then-read evidence buffer on disk.

    Writes are buffered and flushed in batches. A read of any kind flushes first,
    so a caller never sees a partial view; in a pipeline run there is exactly one
    such flush, because every write happens before any read.
    """

    __slots__ = ("_connection", "_count", "_dictionary", "_flush_batch", "_path", "_pending")

    def __init__(self, path: Path, *, flush_batch: int = DEFAULT_FLUSH_BATCH) -> None:
        if flush_batch < 1:
            msg = f"flush_batch={flush_batch} must be at least 1."
            raise EvidenceError(msg)
        self._path = path
        self._flush_batch = flush_batch
        self._pending: dict[str, EvidenceItem] = {}
        self._count = 0
        self._dictionary: bytes | None = None

        path.parent.mkdir(parents=True, exist_ok=True)
        # A previous run of the same case derives the same run id and therefore the
        # same filename. Starting from whatever it left behind would silently merge
        # two runs' evidence, so the file is replaced, never appended to.
        path.unlink(missing_ok=True)
        self._connection = sqlite3.connect(path, isolation_level=None)
        # page_size must be set BEFORE the first table exists or it is ignored.
        # 8 KiB rather than the 4 KiB default: it halves the B-tree depth for the
        # same data and leaves room for a compressed document to sit inline.
        self._connection.execute("PRAGMA page_size = 8192")
        # Durability is worthless here: the file is deleted at the end of the run,
        # and a crash mid-run invalidates it whatever the journal says. Turning both
        # off is what makes a multi-million-row write tolerable.
        self._connection.execute("PRAGMA journal_mode = OFF")
        self._connection.execute("PRAGMA synchronous = OFF")
        self._connection.execute("PRAGMA temp_store = MEMORY")
        self._connection.execute("PRAGMA cache_size = -65536")  # 64 MiB page cache
        self._connection.executescript(_SCHEMA)

    # ------------------------------------------------------------------ writing

    @property
    def path(self) -> Path:
        return self._path

    def add(self, item: EvidenceItem) -> EvidenceItem:
        """Buffer one item, returning the ledger's copy of it.

        A collision with something already buffered is detected here and raises
        immediately. A collision with something already flushed is detected at the
        next flush — at most ``flush_batch`` additions later, always within the same
        run, and with the same message. Detecting it immediately would need every
        ID ever written to stay resident, which is the memory this class exists to
        avoid.
        """
        existing = self._pending.get(item.evidence_id)
        if existing is not None:
            if existing != item:
                raise _collision(item.evidence_id)
            return existing
        self._pending[item.evidence_id] = item
        if len(self._pending) >= self._flush_batch:
            self.flush()
        return item

    def extend(self, items: Sequence[EvidenceItem]) -> None:
        for item in items:
            self.add(item)

    def _ensure_dictionary(self, documents: Sequence[str]) -> bytes:
        """The shared deflate dictionary, built once and stored in the file.

        Derived from the first flush's own documents, so this module never has to
        know what any other stage's boilerplate says, and persisted so that a row
        can always be read back by the file that holds it.
        """
        if self._dictionary is not None:
            return self._dictionary
        row = self._connection.execute("SELECT value FROM meta WHERE key = 'dictionary'").fetchone()
        dictionary = bytes(row[0]) if row is not None else build_dictionary(documents)
        if row is None:
            self._connection.execute(
                "INSERT INTO meta (key, value) VALUES ('dictionary', ?)", (dictionary,)
            )
        self._dictionary = dictionary
        return dictionary

    def flush(self) -> None:
        """Write the buffer, refusing any ID reused for different content."""
        if not self._pending:
            return
        pending = self._pending
        self._pending = {}
        connection = self._connection

        documents = [encode(item) for item in pending.values()]
        dictionary = self._ensure_dictionary(documents)

        rows: list[tuple[str, str, str, str, str, str, bytes]] = [
            (
                item.evidence_id,
                item.subject_kind,
                item.subject_id,
                item.category.value,
                item.direction.value,
                _fingerprint(document),
                compress(document, dictionary),
            )
            for item, document in zip(pending.values(), documents, strict=True)
        ]
        # Order the batch by primary key before writing. SQLite does not care, but a
        # deterministic write order keeps two runs' page layouts comparable when
        # someone inevitably diffs the file while debugging.
        rows.sort(key=lambda row: row[0])

        known = {row[0]: row[1] for row in self._existing_hashes([row[0] for row in rows])}
        for row in rows:
            stored = known.get(row[0])
            if stored is not None and stored != row[5]:
                raise _collision(row[0])

        connection.execute("BEGIN")
        connection.executemany(
            "INSERT OR IGNORE INTO ledger "
            "(evidence_id, subject_kind, subject_id, category, direction, content_sha, document) "
            "VALUES (?, ?, ?, ?, ?, ?, ?)",
            rows,
        )
        connection.execute("COMMIT")
        self._count += sum(1 for row in rows if row[0] not in known)

    def _dictionary_for_read(self) -> bytes:
        """The stored dictionary, or empty when nothing has been written yet.

        The empty case is deliberately NOT cached: a read before the first write
        would otherwise pin ``b""`` as this spill's dictionary and every later row
        would be compressed against nothing. Correct, but 7x larger, and silently
        so.
        """
        if self._dictionary is not None:
            return self._dictionary
        row = self._connection.execute("SELECT value FROM meta WHERE key = 'dictionary'").fetchone()
        if row is None:
            return b""
        self._dictionary = bytes(row[0])
        return self._dictionary

    def _existing_hashes(self, ids: Sequence[str]) -> list[tuple[str, str]]:
        """``(evidence_id, content_sha)`` for the subset of ``ids`` already stored."""
        found: list[tuple[str, str]] = []
        # SQLITE_MAX_VARIABLE_NUMBER is 999 on older builds; chunk rather than
        # depend on the local library's limit.
        for start in range(0, len(ids), 500):
            chunk = ids[start : start + 500]
            placeholders = ",".join("?" for _ in chunk)
            cursor = self._connection.execute(
                f"SELECT evidence_id, content_sha FROM ledger WHERE evidence_id IN ({placeholders})",  # noqa: E501,S608 - placeholders only; every value is bound
                tuple(chunk),
            )
            found.extend((str(row[0]), str(row[1])) for row in cursor.fetchall())
        return found

    # ------------------------------------------------------------------ reading

    def __len__(self) -> int:
        self.flush()
        return self._count

    def get(self, evidence_id: str) -> EvidenceItem | None:
        pending = self._pending.get(evidence_id)
        if pending is not None:
            return pending
        row = self._connection.execute(
            "SELECT document FROM ledger WHERE evidence_id = ?", (evidence_id,)
        ).fetchone()
        if row is None:
            return None
        return decode(decompress(bytes(row[0]), self._dictionary_for_read()))

    def __contains__(self, evidence_id: str) -> bool:
        if evidence_id in self._pending:
            return True
        row = self._connection.execute(
            "SELECT 1 FROM ledger WHERE evidence_id = ?", (evidence_id,)
        ).fetchone()
        return row is not None

    def iter_items(self) -> Iterator[EvidenceItem]:
        """Every item, in ``ledger._sort_key`` order, one at a time."""
        self.flush()
        dictionary = self._dictionary_for_read()
        cursor = self._connection.execute(f"SELECT document FROM ledger ORDER BY {_ORDER_COLUMNS}")  # noqa: S608 - constant column list
        for row in cursor:
            yield decode(decompress(bytes(row[0]), dictionary))

    def iter_for_subject(self, subject_id: str) -> Iterator[EvidenceItem]:
        self.flush()
        dictionary = self._dictionary_for_read()
        cursor = self._connection.execute(
            f"SELECT document FROM ledger WHERE subject_id = ? ORDER BY {_ORDER_COLUMNS}",  # noqa: S608 - constant column list
            (subject_id,),
        )
        for row in cursor:
            yield decode(decompress(bytes(row[0]), dictionary))

    def iter_contradictions(self) -> Iterator[EvidenceItem]:
        self.flush()
        dictionary = self._dictionary_for_read()
        cursor = self._connection.execute(
            f"SELECT document FROM ledger WHERE direction = ? ORDER BY {_ORDER_COLUMNS}",  # noqa: S608 - constant column list
            (EvidenceDirection.CONTRADICTS.value,),
        )
        for row in cursor:
            yield decode(decompress(bytes(row[0]), dictionary))

    # ------------------------------------------------------------------ closing

    def close(self) -> None:
        """Close the connection and remove the file. Idempotent."""
        try:
            self._connection.close()
        finally:
            self._pending = {}
            self._path.unlink(missing_ok=True)


def _collision(evidence_id: str) -> EvidenceError:
    return EvidenceError(
        f"Evidence ID {evidence_id!r} was reused for a different claim. "
        "Evidence IDs are content-derived; a collision means the ID was "
        "constructed from the wrong inputs."
    )


__all__ = [
    "DEFAULT_FLUSH_BATCH",
    "DICTIONARY_BYTES",
    "SqliteEvidenceSpill",
    "build_dictionary",
    "compress",
    "decode",
    "decompress",
    "encode",
]
