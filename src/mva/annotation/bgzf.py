"""Is this bgzipped file, and the tabix index beside it, safe to read as data?

Both local-VCF adapters in this package reach their release through htslib, and
both can be handed a file pair that htslib will read *successfully* and
*incompletely*. That combination is the dangerous one here: a region query against
a short file or a short index does not raise, it returns fewer records, and a
variant with no record is scored by this pipeline as novel and ultra-rare —
the strongest promoting signal the ranker has (GP-14).

The rule therefore lives here once and both adapters call it, for the same reason
:mod:`mva.alleles` exists: the previous shape of this repository had one adapter
checking a property and the other not, which is indistinguishable from neither
checking it until someone reads both files.

Two independent things can be short, and each hides the other:

* **The data.** A partially downloaded ``.vcf.bgz`` decompresses perfectly up to
  its truncation point. :func:`has_bgzf_eof` looks for the empty end-of-file block
  every complete bgzip stream ends with.
* **The index.** An index built from a shorter file names only the blocks it saw.
  htslib will query a *complete* data file through it without complaint and answer
  "no record" for everything past its reach. :func:`index_covers_data` measures
  that reach out of the index itself.

Nothing here uses mtime. htslib warns when an index is older than its data, and on
the real gnomAD v4.1 exomes release it warns on all 24 shards — the ``.tbi`` files
finished downloading before the multi-gigabyte shards did. Every one of those
indexes is complete. An mtime gate would reject the entire real release while
still missing an index rebuilt from a truncated file, so the check is on content.
"""

from __future__ import annotations

import gzip
import struct
import zlib
from pathlib import Path
from typing import Final

__all__ = [
    "BGZF_EOF",
    "has_bgzf_eof",
    "index_covers_data",
    "index_path_for",
    "tabix_index_reach",
]

#: The 28-byte empty-BGZF block every complete bgzip stream ends with. Its absence
#: means the file is truncated — mid-download, interrupted, or corrupt.
BGZF_EOF: Final[bytes] = (
    b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00"
    b"\x00\x00\x00\x00\x00\x00\x00\x00"
)

#: Index suffixes tabix/htslib will accept for a bgzipped VCF.
_INDEX_SUFFIXES: Final[tuple[str, ...]] = (".tbi", ".csi")

#: Leading bytes of a tabix index. A CSI index carries ``CSI\1`` and a different
#: layout; :func:`tabix_index_reach` measures only the ``.tbi`` form and says so.
_TBI_MAGIC: Final[bytes] = b"TBI\x01"


def index_path_for(path: Path) -> Path | None:
    """The tabix/CSI index beside ``path``, or ``None`` when neither exists."""
    for suffix in _INDEX_SUFFIXES:
        candidate = path.with_name(path.name + suffix)
        if candidate.is_file():
            return candidate
    return None


def has_bgzf_eof(path: Path) -> bool:
    """Whether the file ends with the empty-BGZF block that terminates every bgzip stream.

    A partially downloaded ``.vcf.bgz`` decompresses perfectly up to its truncation
    point, so a reader that does not look for the EOF marker gets a *silently short*
    dataset: every variant past the truncation point reports as absent, which
    downstream reads as "novel, therefore interesting".
    """
    try:
        size = path.stat().st_size
    except OSError:
        return False
    if size < len(BGZF_EOF):
        return False
    with path.open("rb") as handle:
        handle.seek(size - len(BGZF_EOF))
        return handle.read(len(BGZF_EOF)) == BGZF_EOF


def tabix_index_reach(index_path: Path) -> int | None:
    """The furthest byte offset in the data file that this ``.tbi`` index can reach.

    Parsed straight out of the index rather than inferred, because the question it
    answers has no other source. Every chunk and linear-index entry in a tabix index
    is a *virtual offset*: the high 48 bits are the start of a BGZF block in the
    data file, the low 16 bits an offset inside it (htslib's tabix specification).
    The maximum block offset across all of them is the last place in the data file
    the index knows about.

    Returns ``None`` — unknown, not zero — when the index is a CSI (a different
    layout, needed only for contigs over 512 Mb and shipped by neither gnomAD nor
    ClinVar) or cannot be parsed. ``None`` is deliberately not an accusation.
    """
    try:
        raw = gzip.decompress(index_path.read_bytes())
    except (OSError, EOFError, gzip.BadGzipFile, zlib.error):
        return None
    if not raw.startswith(_TBI_MAGIC):
        return None
    try:
        n_ref, _fmt, _sq, _bg, _en, _meta, _skip, name_bytes = struct.unpack_from("<8i", raw, 4)
        offset = 4 + 32 + name_bytes
        reach = 0
        for _ in range(n_ref):
            (bin_count,) = struct.unpack_from("<i", raw, offset)
            offset += 4
            for _ in range(bin_count):
                _bin, chunk_count = struct.unpack_from("<Ii", raw, offset)
                offset += 8
                for _ in range(chunk_count):
                    _begin, end = struct.unpack_from("<QQ", raw, offset)
                    offset += 16
                    reach = max(reach, end >> 16)
            (interval_count,) = struct.unpack_from("<i", raw, offset)
            offset += 4
            for _ in range(interval_count):
                (virtual,) = struct.unpack_from("<Q", raw, offset)
                offset += 8
                reach = max(reach, virtual >> 16)
    except (struct.error, ValueError):
        return None
    return reach


def index_covers_data(path: Path, index_path: Path | None) -> bool | None:
    """Whether the index beside ``path`` reaches the end of ``path``'s data.

    This is the second door into the failure :func:`has_bgzf_eof` closes. The EOF
    marker proves the *data* stream is whole and says nothing about the *index*, and
    htslib will region-query a complete data file through an index built from a
    shorter one perfectly happily — returning nothing for every record past the
    point the index reaches.

    The bound is exact rather than a tolerance. A complete index reaches either the
    file size (gnomAD's own indexes, whose last chunk ends past the final data
    block) or ``size - len(BGZF_EOF)`` (an index built locally by ``tabix``, whose
    last chunk ends at the start of the empty end-of-file block). Anything short of
    that leaves data blocks the index cannot name. Measured against all 24 real
    v4.1 exome shards and the 193 MB ClinVar release: every one reaches exactly one
    of those two values.

    Returns ``None`` when :func:`tabix_index_reach` could not measure the index.
    """
    if index_path is None:
        return None
    reach = tabix_index_reach(index_path)
    if reach is None:
        return None
    try:
        size = path.stat().st_size
    except OSError:
        return None
    return reach >= size - len(BGZF_EOF)
