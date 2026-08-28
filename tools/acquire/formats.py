"""Proving a downloaded file is the format its name claims.

A sha256 pins whatever arrived. It does not notice that what arrived is an HTML
error page — which happened here for real: a stale Gene2Phenotype bulk-download
URL returned a JavaScript app shell with HTTP 200, and ``file -b`` said "HTML
document text". Hashing that page produces a perfectly valid hash of the wrong
thing, and every downstream check then passes forever.

``tools.acquire.fetch.sniff_content_mismatch`` is the cheap tripwire for that: 64
bytes, magic-number only, run on every fetch. This module is the deep version,
run once at registration, that actually **opens the file as the format it claims
to be** — a BGZF VCF whose tabix index loads and whose first record parses, a
FASTA whose ``.fai`` really indexes it, a GTF that decompresses into nine
tab-separated columns. The result is recorded in the manifest
(:class:`mva.resources.IntegrityRecord`) so a reader can see what was proven
rather than inferring it from the presence of a hash.

Every detail string here describes **public reference data only**. No patient
coordinate has a code path into this module: its only input is a path under the
resource root.

``pysam`` is an optional extra (``uv sync --extra genomics``). Without it the
structural checks that need htslib degrade to magic-number checks and say so in
the recorded detail, rather than silently reporting a weaker check as a stronger
one.
"""

from __future__ import annotations

import csv
import gzip
import io
import struct
import zipfile
from collections.abc import Callable
from pathlib import Path
from typing import Any, Final

from mva.resources import FormatCheck, IndexCheck, ResourceKind

#: Read from the head of a text/gzip member when checking column shape.
_PROBE_LINES: Final = 40

#: How many contigs of a FASTA to actually seek into when cross-checking its .fai.
#: Three rather than all 195: the check is looking for an index built against a
#: different file, which is a whole-file property, not a per-contig one.
_FAIDX_PROBE_CONTIGS: Final = 3

_TBI_MAGIC: Final = b"TBI\x01"
_CSI_MAGIC: Final = b"CSI\x01"
_BGZF_EOF: Final = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")

_HTML_SNIFFS: Final[tuple[bytes, ...]] = (b"<!doctype html", b"<html", b"<!DOCTYPE HTML")


class FormatProbeResult:
    """What a probe concluded. Deliberately not frozen-model machinery: this type
    never leaves the tool, and the manifest stores the strings, not the object."""

    __slots__ = ("check", "detail", "index_check", "index_detail", "problem")

    def __init__(
        self,
        check: FormatCheck,
        detail: str,
        problem: str | None = None,
        *,
        index_check: IndexCheck = IndexCheck.NOT_APPLICABLE,
        index_detail: str = "",
    ) -> None:
        self.check = check
        self.detail = detail
        self.problem = problem
        self.index_check = index_check
        self.index_detail = index_detail

    @property
    def ok(self) -> bool:
        return self.problem is None

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"FormatProbeResult({self.check.value!r}, {self.detail!r}, {self.problem!r})"


def _pysam() -> Any | None:
    try:
        import pysam  # noqa: PLC0415 - optional extra, imported only when a probe needs it
    except ImportError:
        return None
    return pysam


def _bgzf_eof_present(path: Path) -> bool:
    """Does the file end with BGZF's mandatory 28-byte end-of-file block?

    This is the single most useful cheap check on a multi-gigabyte bgzipped file:
    the block is written last, so its presence means the writer finished. A
    stalled 12 GB download has a perfect head, a plausible size and no EOF block.
    """
    size = path.stat().st_size
    if size < len(_BGZF_EOF):
        return False
    with path.open("rb") as handle:
        handle.seek(size - len(_BGZF_EOF))
        return handle.read(len(_BGZF_EOF)) == _BGZF_EOF


def _looks_like_html(path: Path) -> bool:
    with path.open("rb") as handle:
        head = handle.read(64).lower()
    return any(head.startswith(sniff.lower()) for sniff in _HTML_SNIFFS)


def _gzip_head_lines(path: Path, count: int = _PROBE_LINES) -> list[str]:
    lines: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            lines.append(line.rstrip("\n"))
            if len(lines) >= count:
                break
    return lines


def _text_head_lines(path: Path, count: int = _PROBE_LINES) -> list[str]:
    lines: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            lines.append(line.rstrip("\n"))
            if len(lines) >= count:
                break
    return lines


# --------------------------------------------------------------------------- probes


class _TabixIndex:
    """The parts of a ``.tbi`` that say which bytes of which file it describes.

    Parsed here rather than through pysam because pysam exposes the index as an
    *accessor* -- contig names and a fetch API -- and deliberately hides the raw
    virtual offsets. Those offsets are the evidence: a BGZF virtual offset is
    ``(block_offset << 16) | offset_within_block``, so the largest one an index
    references is a hard statement about how big the file it was built from was.
    An index whose deepest reference lies beyond the end of the data file cannot
    possibly be an index OF that file, and no amount of successful header-reading
    will reveal that.

    Format: SAMtools tabix specification, section 'The Tabix index format'.
    """

    __slots__ = ("contigs", "linear_intervals", "max_block_offset")

    def __init__(
        self, contigs: list[str], max_block_offset: int, linear_intervals: dict[str, int]
    ) -> None:
        self.contigs = contigs
        self.max_block_offset = max_block_offset
        self.linear_intervals = linear_intervals


#: The tabix linear index buckets records into fixed 16 kb windows.
_TBI_LINEAR_SHIFT: Final = 14


def _parse_tabix_index(index_path: Path) -> _TabixIndex | None:
    """Read a ``.tbi``. Returns ``None`` for anything that is not one."""
    try:
        with gzip.open(index_path, "rb") as handle:
            blob = handle.read()
    except (OSError, EOFError):
        return None
    if len(blob) < 36 or blob[:4] != _TBI_MAGIC:
        return None

    cursor = 4
    (n_ref,) = struct.unpack_from("<i", blob, cursor)
    cursor += 4
    cursor += 24  # format, col_seq, col_beg, col_end, meta, skip
    (l_nm,) = struct.unpack_from("<i", blob, cursor)
    cursor += 4
    names = blob[cursor : cursor + l_nm].split(b"\x00")
    contigs = [name.decode("ascii", errors="replace") for name in names if name]
    cursor += l_nm

    max_block_offset = 0
    linear_intervals: dict[str, int] = {}
    try:
        for ref in range(n_ref):
            (n_bin,) = struct.unpack_from("<i", blob, cursor)
            cursor += 4
            for _ in range(n_bin):
                cursor += 4  # bin id
                (n_chunk,) = struct.unpack_from("<i", blob, cursor)
                cursor += 4
                for _ in range(n_chunk):
                    _, chunk_end = struct.unpack_from("<QQ", blob, cursor)
                    cursor += 16
                    max_block_offset = max(max_block_offset, chunk_end >> 16)
            (n_intv,) = struct.unpack_from("<i", blob, cursor)
            cursor += 4
            cursor += 8 * n_intv
            if ref < len(contigs):
                linear_intervals[contigs[ref]] = n_intv
    except struct.error:
        return None
    return _TabixIndex(contigs, max_block_offset, linear_intervals)


def _index_describes_this_file(
    data: Path, index: Path, contigs: list[str], variants: object
) -> tuple[IndexCheck, str, str | None]:
    """Prove the index was built from ``data``, or say why that could not be shown.

    Returns ``(verdict, detail, problem)``. ``problem`` is non-None only for a
    definitive STALE finding, which fails the resource closed rather than pinning
    an index that does not describe its data.

    Two independent lines of evidence, because they fail differently:

    1. **Structural.** The deepest BGZF block offset the index references must lie
       inside the data file. An index for a larger or different file blows past the
       end, and this is arithmetic on ~100 kB of index -- no seeking, no htslib.
    2. **Behavioural.** A targeted fetch, aimed at the deepest region the index's
       own linear index claims to cover, must return records inside the window that
       was asked for. This is what catches an index of the *same size* built from
       different data.

    mtime is deliberately absent from both. htslib warns when an index is older
    than its data, which is true for every gnomAD shard here purely because each
    small ``.tbi`` finished downloading before its multi-gigabyte ``.bgz`` did.
    Timestamps are not evidence about content (ADR 0020).
    """
    data_size = data.stat().st_size
    skew = index.stat().st_mtime < data.stat().st_mtime
    skew_note = (
        "; index mtime precedes data mtime (htslib warns about this) -- ordinary "
        "download ordering, and not evidence either way, so it was checked by access"
        if skew
        else ""
    )

    parsed = _parse_tabix_index(index)
    if parsed is None:
        return (
            IndexCheck.NOT_CHECKED,
            f"index {index.name} could not be parsed as a tabix index{skew_note}",
            None,
        )

    # STRICTLY greater. A chunk end whose block offset equals the file size is the
    # canonical "the last record runs to the end" terminator, and gnomAD's shards
    # all use it -- an earlier `>=` here flagged every one of the 25 as stale,
    # which would have failed 184.8 GB of good data closed.
    if parsed.max_block_offset > data_size:
        return (
            IndexCheck.STALE,
            f"index references block offset {parsed.max_block_offset} past EOF ({data_size})",
            (
                f"the index {index.name!r} references BGZF block offset "
                f"{parsed.max_block_offset}, which lies beyond the end of the {data_size}-byte "
                "data file. It cannot be an index of this file -- it was built from a larger "
                "or different one. Re-fetch the pair; do not pin them together."
            ),
        )

    reach = parsed.max_block_offset / data_size if data_size else 0.0
    detail = (
        f"index {index.name} references up to block offset {parsed.max_block_offset} "
        f"({reach:.1%} into the {data_size}-byte file)"
    )

    probed = _deep_fetch_probe(parsed, contigs, variants)
    if probed is None:
        return (IndexCheck.NOT_CHECKED, f"{detail}; no deep region to probe{skew_note}", None)
    ok, probe_detail = probed
    if not ok:
        return (
            IndexCheck.STALE,
            f"{detail}; {probe_detail}",
            (
                f"seeking through {index.name!r} returned records outside the window that was "
                f"requested ({probe_detail}). The index does not describe this data file."
            ),
        )
    return (IndexCheck.CONSISTENT, f"{detail}; {probe_detail}{skew_note}", None)


def _deep_fetch_probe(
    parsed: _TabixIndex, contigs: list[str], variants: object
) -> tuple[bool, str] | None:
    """Fetch at the deepest coordinate the linear index claims, and check the answer.

    Aimed deep on purpose. A fetch at the start of a contig is satisfied by the
    first BGZF block, which is where a wrong index is most likely to be
    accidentally right; the far end of the linear index is the offset a mismatched
    index gets wrong.
    """
    for contig in contigs:
        intervals = parsed.linear_intervals.get(contig, 0)
        if intervals < 2:
            continue
        start = (intervals - 1) << _TBI_LINEAR_SHIFT
        end = start + (4 << _TBI_LINEAR_SHIFT)
        try:
            records = list(variants.fetch(contig, start, end))  # type: ignore[attr-defined]
        except (OSError, ValueError) as exc:
            return (False, f"deep fetch on {contig}:{start}-{end} raised {exc.__class__.__name__}")
        if not records:
            # An empty window is not a contradiction: the linear index sizes a
            # contig by its last record's bucket, and the final buckets of a
            # sites-only VCF are legitimately sparse. Try the next contig.
            continue
        outside = [r for r in records if not start < r.pos <= end + 1]
        if outside:
            return (False, f"deep fetch on {contig}:{start}-{end} returned out-of-window records")
        return (
            True,
            f"deep fetch at {contig}:{start}-{end} returned {len(records)} in-window record(s)",
        )
    return None


def _probe_bgzf_vcf(path: Path) -> FormatProbeResult:
    """A BGZF VCF: complete stream, loadable tabix index, parseable first record."""
    if not _bgzf_eof_present(path):
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "BGZF end-of-file block absent",
            "the BGZF EOF block is missing, so this stream was never finished — an "
            "interrupted download, not a complete file",
        )
    index = path.with_suffix(path.suffix + ".tbi")
    if not index.is_file():
        index = path.with_suffix(path.suffix + ".csi")
    if not index.is_file():
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "BGZF EOF present; no .tbi/.csi index alongside",
            "no tabix index found next to the VCF; random access is impossible without it",
        )

    pysam = _pysam()
    if pysam is None:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            f"BGZF EOF present, index {index.name} present; pysam unavailable so the "
            "index was NOT opened",
            None,
            index_check=IndexCheck.NOT_CHECKED,
            index_detail="pysam unavailable; install the 'genomics' extra to cross-check",
        )

    try:
        with pysam.TabixFile(str(path)) as tabix:
            contigs = list(tabix.contigs)
        with pysam.VariantFile(str(path)) as variants:
            record = next(iter(variants.fetch()), None)
            samples = len(variants.header.samples)
            index_verdict, index_detail, index_problem = _index_describes_this_file(
                path, index, contigs, variants
            )
    except (OSError, ValueError) as exc:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "index open failed",
            f"htslib could not open the file and its index ({exc.__class__.__name__})",
        )

    if not contigs:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "tabix index opened but lists no contigs",
            "the tabix index is empty; it indexes no reference sequence at all",
        )
    if record is None:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            f"tabix index opened; {len(contigs)} contig(s); no data records",
            "the VCF carries a header but not a single data record",
        )
    if index_problem is not None:
        # Fails CLOSED. A stale index is not a warning to record next to a hash --
        # it silently returns the wrong records, or none, for exactly the rare
        # coordinates this pipeline cares about.
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            index_detail,
            index_problem,
            index_check=index_verdict,
            index_detail=index_detail,
        )

    first = f"{record.chrom}:{record.pos}"
    return FormatProbeResult(
        FormatCheck.BGZF_TABIX_INDEXED,
        (
            f"BGZF EOF present; tabix index {index.name} opened; {len(contigs)} contig(s) "
            f"({contigs[0]}...); sites-only={samples == 0}; first record {first}"
        ),
        index_check=index_verdict,
        index_detail=index_detail,
    )


def _probe_tabix_index(path: Path) -> FormatProbeResult:
    try:
        with gzip.open(path, "rb") as handle:
            magic = handle.read(4)
    except OSError as exc:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "not a readable gzip member",
            f"could not decompress the index ({exc.__class__.__name__})",
        )
    if magic == _TBI_MAGIC:
        return FormatProbeResult(FormatCheck.TABIX_INDEX, "gzip member carrying TBI\\1 magic")
    if magic == _CSI_MAGIC:
        return FormatProbeResult(FormatCheck.TABIX_INDEX, "gzip member carrying CSI\\1 magic")
    return FormatProbeResult(
        FormatCheck.NOT_CHECKED,
        "gzip member with unexpected magic",
        "decompresses, but does not begin with the TBI\\1 or CSI\\1 magic an index must",
    )


def _probe_fasta(path: Path) -> FormatProbeResult:
    """A FASTA cross-checked against its own ``.fai``.

    The interesting failure is an index built against a *different* FASTA — the
    two files are downloaded and generated separately, and a mismatched pair
    yields sequence that is silently offset. Reading a slice through the index and
    checking it against the same slice read directly is what catches it.
    """
    index = Path(f"{path}.fai")
    if not index.is_file():
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "no .fai alongside",
            "the FASTA has no faidx index; random reference access is impossible",
        )
    with path.open("rb") as handle:
        if handle.read(1) != b">":
            return FormatProbeResult(
                FormatCheck.NOT_CHECKED,
                "first byte is not '>'",
                "does not begin with a FASTA header line",
            )

    records = _parse_fai(index)
    if not records:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED, ".fai parsed to no records", "the faidx index is empty"
        )

    pysam = _pysam()
    if pysam is None:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            f"FASTA header present; .fai lists {len(records)} contig(s); pysam unavailable "
            "so offsets were NOT cross-checked",
            None,
        )

    try:
        with pysam.FastaFile(str(path)) as fasta:
            indexed = list(fasta.references)
            probes: list[str] = []
            for name, length, offset, linebases in records[:_FAIDX_PROBE_CONTIGS]:
                width = min(60, length)
                through_index = fasta.fetch(name, 0, width).upper()
                direct = _read_raw_bases(path, offset, width, linebases).upper()
                if through_index != direct:
                    return FormatProbeResult(
                        FormatCheck.NOT_CHECKED,
                        f"faidx offset for {name} does not land on its sequence",
                        f"the .fai is not an index of this FASTA: contig {name}'s declared "
                        "offset reads different bases than the index returns. Re-run "
                        "`samtools faidx` on this exact file.",
                        index_check=IndexCheck.STALE,
                        index_detail=f"contig {name}: indexed read disagrees with a raw read",
                    )
                probes.append(f"{name}:{length}")
    except (OSError, ValueError) as exc:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "faidx open failed",
            f"htslib could not open the FASTA with its index ({exc.__class__.__name__})",
        )

    total = sum(length for _, length, _, _ in records)
    cross_read = (
        f"faidx offsets cross-checked against raw reads for {', '.join(probes)}: each "
        "declared offset lands on that contig's own sequence"
    )
    return FormatProbeResult(
        FormatCheck.FASTA_FAIDX_CONSISTENT,
        f"{len(indexed)} contig(s), {total} bases indexed; {cross_read}",
        index_check=IndexCheck.CONSISTENT,
        index_detail=cross_read,
    )


def _parse_fai(index: Path) -> list[tuple[str, int, int, int]]:
    """``(name, length, offset, linebases)`` per faidx record; empty if malformed."""
    records: list[tuple[str, int, int, int]] = []
    for line in _text_head_lines(index, count=1_000_000):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            return []
        try:
            records.append((fields[0], int(fields[1]), int(fields[2]), int(fields[3])))
        except ValueError:
            return []
    return records


def _read_raw_bases(path: Path, offset: int, count: int, linebases: int) -> str:
    """``count`` bases starting at faidx ``offset``, skipping the line breaks."""
    # Over-read to cover the newlines interleaved in the requested span.
    span = count + (count // max(linebases, 1)) + 2
    with path.open("rb") as handle:
        handle.seek(offset)
        blob = handle.read(span).decode("ascii", errors="replace")
    return "".join(blob.split())[:count]


def _probe_fasta_index(path: Path) -> FormatProbeResult:
    records = _parse_fai(path)
    if not records:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "not a parseable faidx index",
            "does not parse as five tab-separated faidx columns",
        )
    return FormatProbeResult(
        FormatCheck.FASTA_INDEX,
        f"{len(records)} faidx record(s); first {records[0][0]} length {records[0][1]}",
    )


def _probe_gzip_gtf(path: Path) -> FormatProbeResult:
    try:
        lines = _gzip_head_lines(path)
    except (OSError, EOFError) as exc:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "gzip member unreadable",
            f"could not decompress ({exc.__class__.__name__})",
        )
    data = [line for line in lines if line and not line.startswith("#")]
    if not data:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "decompressed to comments only",
            f"no data line in the first {_PROBE_LINES} decompressed lines",
        )
    fields = data[0].split("\t")
    if len(fields) != 9:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            f"first data line has {len(fields)} columns",
            "a GTF data line must have exactly 9 tab-separated columns",
        )
    try:
        start, end = int(fields[3]), int(fields[4])
    except ValueError:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "GTF start/end are not integers",
            "columns 4 and 5 of a GTF line must be integer coordinates",
        )
    header = len(lines) - len(data)
    return FormatProbeResult(
        FormatCheck.GZIP_TEXT,
        (
            f"gzip GTF; {header} comment line(s); first feature {fields[0]}:{start}-{end} "
            f"type={fields[2]} source={fields[1]}"
        ),
    )


def _probe_gzip_table(path: Path) -> FormatProbeResult:
    try:
        lines = _gzip_head_lines(path)
    except (OSError, EOFError) as exc:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "gzip member unreadable",
            f"could not decompress ({exc.__class__.__name__})",
        )
    if not lines:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED, "decompressed to nothing", "the gzip member is empty"
        )
    columns = len(lines[0].split("\t"))
    if columns < 2:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            f"header row has {columns} column(s)",
            "the first decompressed line does not split into a tab-separated header",
        )
    return FormatProbeResult(
        FormatCheck.GZIP_TEXT,
        f"gzip TSV; {columns}-column header ({lines[0].split(chr(9))[0]}...); "
        f"{len(lines)} line(s) probed",
    )


def _probe_gzip_fasta(path: Path) -> FormatProbeResult:
    """A gzipped FASTA. Distinct from :func:`_probe_gzip_table` because a FASTA is
    not a table and checking it as one reports a perfectly good file as broken --
    which is exactly what this probe did to the GRCh38 archive on its first run."""
    try:
        lines = _gzip_head_lines(path, count=8)
    except (OSError, EOFError) as exc:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "gzip member unreadable",
            f"could not decompress ({exc.__class__.__name__})",
        )
    if not lines or not lines[0].startswith(">"):
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "first decompressed line is not a FASTA header",
            "decompresses, but does not begin with a '>' header line",
        )
    sequence = "".join(lines[1:5]).upper()
    if sequence and not set(sequence) <= set("ACGTNRYKMSWBDHV"):
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "non-IUPAC characters after the header",
            "the lines after the FASTA header are not nucleotide sequence",
        )
    return FormatProbeResult(
        FormatCheck.GZIP_TEXT,
        f"gzip FASTA; first record {lines[0][:60]}; {len(sequence)} sequence bases probed",
    )


def _probe_properties(path: Path) -> FormatProbeResult:
    """A ``key : value`` / ``key=value`` configuration file.

    For ``snpEff.config`` this is load-bearing rather than cosmetic: SnpEff
    resolves a genome name THROUGH this file, so a config that parses but is
    missing the genome stanza makes an installed database report as absent. The
    probe therefore reports which ``*.genome`` entries it found.
    """
    # Whole file, not a head: snpEff.config is 19 MB and the genome stanza that
    # makes it load-bearing can be anywhere in it. A 19 MB text read is ~50 ms,
    # and a probe that reports "no genome declared" because it stopped reading
    # early would be worse than no probe.
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        return FormatProbeResult(FormatCheck.NOT_CHECKED, "empty file", "the file has no content")
    settings = [
        line
        for line in lines
        if line.strip() and not line.lstrip().startswith("#") and (":" in line or "=" in line)
    ]
    if not settings:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "no key/value lines",
            "no line parses as 'key : value' or 'key=value'",
        )
    genomes = sorted(
        line.split(".genome")[0].strip() for line in settings if ".genome" in line.split(":")[0]
    )
    detail = f"{len(settings)} setting(s) over {len(lines)} line(s)"
    if genomes:
        detail += f"; {len(genomes)} genome stanza(s) declared"
    return FormatProbeResult(FormatCheck.PLAIN_TEXT, detail)


def _probe_obo(path: Path) -> FormatProbeResult:
    lines = _text_head_lines(path)
    if not lines or not lines[0].startswith("format-version:"):
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "no OBO format-version header",
            "an .obo file must begin with a 'format-version:' line",
        )
    version = next((line for line in lines if line.startswith("data-version:")), "")
    return FormatProbeResult(
        FormatCheck.PLAIN_TEXT, f"OBO {lines[0]}" + (f"; {version}" if version else "")
    )


def _probe_text_table(path: Path) -> FormatProbeResult:
    lines = _text_head_lines(path)
    if not lines:
        return FormatProbeResult(FormatCheck.NOT_CHECKED, "empty file", "the file has no content")
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    header = next((line for line in lines if line.strip() and not line.startswith("#")), lines[0])
    columns = len(next(csv.reader(io.StringIO(header), delimiter=delimiter), []))
    # A single-column file is REPORTED, not refused. What this probe can honestly
    # prove is that the bytes decode as text and carry content; "how many columns"
    # is a shape a reviewer reads off the manifest, not a rule the probe should
    # enforce -- a plain one-column list is a legitimate resource, and rejecting it
    # would mark a perfectly good file as unfetchable.
    shape = f"{columns}-column {delimiter!r}-delimited table" if columns > 1 else "single-column"
    return FormatProbeResult(FormatCheck.PLAIN_TEXT, f"{shape}; {len(lines)} line(s) probed")


def _probe_zip(path: Path) -> FormatProbeResult:
    try:
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
            broken = archive.testzip()
    except (OSError, zipfile.BadZipFile) as exc:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "not a readable zip",
            f"the archive did not open ({exc.__class__.__name__})",
        )
    if broken is not None:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "zip central directory lists a corrupt member",
            f"member {broken!r} fails its CRC",
        )
    return FormatProbeResult(FormatCheck.JAVA_ARCHIVE, f"zip/jar with {len(names)} entries")


def _probe_snpeff_database(path: Path) -> FormatProbeResult:
    predictor = path / "snpEffectPredictor.bin"
    if not predictor.is_file() or predictor.stat().st_size == 0:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "snpEffectPredictor.bin missing or empty",
            "a SnpEff genome database without a non-empty snpEffectPredictor.bin is not a "
            "database; SnpEff would report the genome as not installed",
        )
    sequences = sorted(path.glob("sequence*.bin"))
    empty = [p.name for p in path.rglob("*") if p.is_file() and p.stat().st_size == 0]
    detail = (
        f"snpEffectPredictor.bin {predictor.stat().st_size} bytes; "
        f"{len(sequences)} sequence*.bin file(s)"
    )
    if empty:
        detail += f"; {len(empty)} zero-length member(s): {', '.join(sorted(empty))}"
    if not sequences:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            detail,
            "no sequence*.bin files: the database cannot supply codons, so HGVS.c/HGVS.p "
            "would be absent from every annotation",
        )
    return FormatProbeResult(FormatCheck.SNPEFF_DATABASE, detail)


def _probe_opaque(path: Path) -> FormatProbeResult:
    size = path.stat().st_size
    if size == 0:
        return FormatProbeResult(FormatCheck.NOT_CHECKED, "zero bytes", "the file is empty")
    return FormatProbeResult(
        FormatCheck.OPAQUE_BINARY,
        f"{size} bytes of binary with no independently checkable structure",
    )


def _probe_json(path: Path) -> FormatProbeResult:
    import json  # noqa: PLC0415 - only this probe needs it

    try:
        payload: Any = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED, "not valid JSON", f"did not parse ({exc.__class__.__name__})"
        )
    keys = sorted(payload) if isinstance(payload, dict) else []
    return FormatProbeResult(
        FormatCheck.PLAIN_TEXT, f"JSON object with keys: {', '.join(keys)}" if keys else "JSON"
    )


# --------------------------------------------------------------------------- dispatch


#: Suffix -> probe, longest suffix first. A table rather than a chain of ``if``s
#: because the ordering IS the logic here: ``.gtf.gz`` and ``.fna.gz`` must both be
#: consulted before the generic ``.gz``, and checking a gzipped FASTA as a TSV
#: reports a perfectly good file as broken (it did, on this table's first run).
_Probe = Callable[[Path], FormatProbeResult]

_PROBES_BY_SUFFIX: Final[tuple[tuple[tuple[str, ...], _Probe], ...]] = (
    ((".vcf.gz", ".vcf.bgz"), _probe_bgzf_vcf),
    ((".tbi", ".csi"), _probe_tabix_index),
    ((".fa", ".fasta", ".fna"), _probe_fasta),
    ((".fai",), _probe_fasta_index),
    ((".gtf.gz",), _probe_gzip_gtf),
    ((".fa.gz", ".fna.gz", ".fasta.gz"), _probe_gzip_fasta),
    ((".gz",), _probe_gzip_table),
    ((".obo",), _probe_obo),
    ((".jar", ".zip"), _probe_zip),
    ((".json",), _probe_json),
    ((".config", ".properties"), _probe_properties),
    ((".tsv", ".csv", ".txt", ".hpoa"), _probe_text_table),
)


def _dispatch(path: Path, kind: ResourceKind) -> _Probe:
    if kind is ResourceKind.DIRECTORY:
        return _probe_snpeff_database
    name = path.name.lower()
    for suffixes, probe in _PROBES_BY_SUFFIX:
        if name.endswith(suffixes):
            return probe
    return _probe_opaque


def probe_format(path: Path, kind: ResourceKind = ResourceKind.FILE) -> FormatProbeResult:
    """Open ``path`` as the format its name claims, and report what was proven.

    Never raises for a bad file: a malformed resource is a *result* to record, not
    an exception to propagate, because the registration pass must be able to
    report every problem in one run rather than stopping at the first.
    """
    exists = path.is_dir() if kind is ResourceKind.DIRECTORY else path.is_file()
    if not exists:
        return FormatProbeResult(FormatCheck.NOT_CHECKED, "absent", "nothing at this path")

    if kind is ResourceKind.FILE and _looks_like_html(path):
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "starts with an HTML document",
            "the file begins with an HTML document rather than the promised data — almost "
            "certainly a server error or landing page saved under the right filename. This "
            "is the failure that produced a valid sha256 of the wrong thing once already; "
            "re-fetch, do not re-pin.",
        )

    try:
        return _dispatch(path, kind)(path)
    except (OSError, ValueError, UnicodeError) as exc:
        return FormatProbeResult(
            FormatCheck.NOT_CHECKED,
            "probe raised",
            f"the format probe failed with {exc.__class__.__name__}",
        )
