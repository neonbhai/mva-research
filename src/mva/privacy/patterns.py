"""The shared detection battery: regexes over BYTES, plus magic-byte sniffers.

Why bytes and not text
----------------------
Genomic files are frequently mis-encoded, truncated, or binary with textual
islands (a BAM header inside BGZF). Decoding first would either raise (losing the
detection) or silently mangle the very bytes we are trying to recognise. Every
rule therefore matches ``bytes``; decoding happens only at the reporting edge and
only for content we never emit.

GP-41 — the cardinal rule
-------------------------
This module is a *detector*, not a *reporter*. Nothing here returns, logs or
stores ``match.group()``. The only permitted derivations of a match are its span,
its length, the rule ID, the line number and the file path. The scanner's output
is itself a leak vector: an agent runs the audit and the result lands in a model
context window, a CI log and a terminal scrollback.

The one concession is :func:`correlation_id`, which lets a caller count *distinct*
matches (e.g. "three different HPO terms") without ever holding onto them. It is
keyed by a per-process random salt that is never persisted, because a plain
truncated hash of a low-entropy value — a seven-digit HPO term, a short MRN — is
brute-forceable in milliseconds. That salt is also why the value must stay in
memory: rendering it into a file makes the file a random function of a secret, and
every artifact here is inside the GP-30 byte-identity claim.

False positives are a design input, not an afterthought
-------------------------------------------------------
A privacy scanner that cries wolf gets disabled, which is strictly worse than one
that is narrower. Each rule therefore carries an explicit
``false_positive_risk`` note, and the two rules with genuinely high collision
rates against legitimate public data — ``vcf_data_line`` (public coordinate
tables) and ``hpo_term`` (the public ontology, and code constants) — are ``warn``
or threshold-gated rather than unconditional failures. Promotion logic lives in
:mod:`mva.privacy.audit`, which has the whole-file context needed to apply it.
"""

from __future__ import annotations

import hmac
import re
import secrets
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Final, Literal

from mva.errors import PrivacyViolationError

#: Rule severities. `fail` blocks; `warn` is recorded and surfaced but does not
#: block unless the audit is run in strict mode.
type Severity = Literal["fail", "warn"]

# ---------------------------------------------------------------------------
# Correlation identifiers (GP-41)
# ---------------------------------------------------------------------------

#: Random, per-process, NEVER written to disk and never included in any artifact.
#: Regenerated on every interpreter start, so a correlation ID is meaningful only
#: within one run and cannot be used to link findings across runs or machines.
_RUN_SALT: Final[bytes] = secrets.token_bytes(16)


def correlation_id(matched: bytes) -> str:
    """A within-PROCESS-stable, non-invertible tag for a matched byte string.

    **In-memory only. Never render this into an artifact.** The key is
    ``secrets.token_bytes(16)`` regenerated at import, so the tag is a random
    function of a secret from one process to the next. Anything that reaches a
    file inherits that randomness, and a file that changes between two identical
    runs is a GP-30 violation.

    That is not hypothetical: ``mva.privacy.audit.redact_path`` used this to label
    redacted path components, ``privacy/privacy_audit.md`` is a registered artifact
    inside the byte-identity claim and is not in ``verify_determinism``'s skip set,
    and the two disagreed across processes. The report now labels by first
    appearance (``mva.privacy.audit.PathRedactor``), which is deterministic AND
    less invertible than this.

    Use this — and only this — when you need to know whether two matches were the
    *same* value (counting distinct HPO terms, deduplicating findings) without
    retaining the value, and the answer lives only in a local variable.

    A bare ``sha256(value)[:8]`` would be useless here: HPO terms, MRNs and
    genomic positions all live in tiny keyspaces, so an unsalted digest is a
    lookup table away from plaintext. The HMAC key is per-process and unpersisted,
    which makes the tag meaningless outside the process that produced it.
    """
    return hmac.new(_RUN_SALT, matched, "sha256").hexdigest()[:8]


def placeholder(rule_id: str, length: int) -> str:
    """The single approved rendering of a match: rule ID and length, nothing else."""
    return f"<REDACTED:{rule_id}:len={length}>"


# ---------------------------------------------------------------------------
# Rule type
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Rule:
    """One detection rule.

    ``severity`` is the rule's severity *in isolation*. Some rules are refined by
    whole-file context in :mod:`mva.privacy.audit` (see ``vcf_data_line`` and
    ``hpo_term``); the value here is the floor, never a promotion.
    """

    rule_id: str
    pattern: re.Pattern[bytes]
    severity: Severity
    description: str
    false_positive_risk: str


# ---------------------------------------------------------------------------
# Patterns
#
# NOTE ON SELF-MATCHING: this file necessarily contains the *source text* of every
# pattern. It does not match itself, and that is deliberate rather than lucky:
# the line-structured rules are anchored with ``^`` under re.MULTILINE (a pattern
# written inside a quoted Python expression is never at the start of a line), and
# the field-structured rules require literal TAB and NEWLINE bytes, which a Python
# source file writes as the two-character escapes ``\t`` and ``\n``.
# ---------------------------------------------------------------------------

_VCF_HEADER: Final = re.compile(rb"^##fileformat=VCFv4\.\d", re.MULTILINE)

_VCF_CHROM_LINE: Final = re.compile(
    rb"^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    re.MULTILINE,
)

# Five tab-separated VCF columns followed by a sixth field boundary. The trailing
# \t is load-bearing: it forces at least CHROM/POS/ID/REF/ALT/QUAL, which a
# two-column coordinate list cannot satisfy.
#
# CHROM is a free-form contig NAME, not a human chromosome number. Restricting it
# to `(?:chr)?(1..22|X|Y|MT)` meant a VCF called against RefSeq accessions
# (`NC_000012.12`), against a patched assembly (`chr12_KI270834v1_alt`) or against
# any non-human reference was invisible to this rule. The discriminating power of
# the rule was never in CHROM: it is in REF/ALT being nucleotide alphabets in
# fixed tab-delimited positions.
_VCF_DATA_LINE: Final = re.compile(
    rb"^[A-Za-z0-9_.\-]{1,32}\t[0-9]{1,12}\t"
    rb"(?:\.|rs[0-9]+|[A-Za-z0-9_:.\-]{1,64})\t"
    rb"[ACGTNacgtn]{1,200}\t"
    rb"(?:[ACGTNacgtn,.*]{1,200}|<[A-Za-z0-9:_]{1,20}>)\t",
    re.MULTILINE,
)

# FORMAT keys are alphanumeric, not alphabetic: Mutect2 emits `GT:AD:AF:F1R2:F2R1`
# and `[A-Z]{1,4}` rejected `F1R2`, so a somatic call set matched nothing. The
# allele pair is also optional-repeating rather than mandatory-paired, because a
# haploid call (chrY, chrM, a hemizygous male X) is written as a single allele
# with no `/` or `|` at all and was therefore invisible.
_GENOTYPE_FIELD: Final = re.compile(rb"\tGT(?::[A-Z0-9]{1,4})*\t[0-9.]+(?:[/|][0-9.]+)*")

# A FASTQ record is recognised structurally, not by extension: header line,
# >=25nt of a homogeneous nucleotide alphabet, separator, >=25 quality chars.
# The homogeneity requirement is what stops ordinary prose from matching.
_FASTQ_RECORD: Final = re.compile(
    rb"^@[!-~][ -~]{0,300}\r?\n[ACGTNacgtn]{25,}\r?\n\+[ -~]{0,300}\r?\n[!-~]{25,}",
    re.MULTILINE,
)

# @RG SM: carries the sample name a sequencing centre actually used, which in
# practice is a hospital accession or the proband's initials plus a date.
# The field separator is `[ \t]` rather than `\t`: the SAM spec is tab-delimited,
# but header lines survive round-trips through editors, `tr`, log excerpts and
# `samtools view -H | column` as space-delimited text, and the sample name is
# just as disclosive there.
_SAM_RG_SAMPLE: Final = re.compile(rb"^@RG[ \t](?:[!-~]+[ \t])*SM:([!-~]+)", re.MULTILINE)

# Both renderings of an HPO identifier: the canonical `HP:0000001` and the OBO/OWL
# underscore form `HP_0000001` used in ontology dumps, RDF exports and file names.
_HPO_TERM: Final = re.compile(rb"\bHP[:_]\d{7}\b")

# A FASTA record: a `>` description line followed by a long homogeneous nucleotide
# run. Patient consensus sequences, extracted read sets and reference slices all
# arrive in this shape, and no rule recognised it at all before.
_FASTA_RECORD: Final = re.compile(
    rb"^>[ -~]{0,300}\r?\n(?:[ACGTNUacgtnu]{40,}\r?\n?)",
    re.MULTILINE,
)

# A PLINK .ped body line: FID IID PAT MAT SEX PHENO then allele pairs. Whitespace
# delimited (PLINK accepts spaces or tabs), so the VCF rules never saw it, and the
# genotype matrix it carries is exactly as disclosive as a VCF.
_PLINK_PED_LINE: Final = re.compile(
    rb"^[!-~]{1,64}[ \t]+[!-~]{1,64}[ \t]+[!-~]{1,64}[ \t]+[!-~]{1,64}[ \t]+"
    rb"[012][ \t]+(?:-9|[012])(?:[ \t]+[ACGTNacgtn0]){6,}",
    re.MULTILINE,
)

# Keyword-anchored ONLY. A bare \d{7,10} run is deliberately NOT matched: it
# collides with every genomic POS in every VCF, every gnomAD allele count and
# every byte offset, which would make the rule useless and then ignored.
# The negative lookahead rejects a purely alphabetic follower, so prose such as
# "the patient id field" does not match, while the same keyword followed by a
# separator and an alphanumeric token does. (No worked example is spelled out
# here: a live one would make this file trip its own scanner, which is both
# correct behaviour and a permanent false alarm.)
#: The separator between a PHI keyword and its value. `,` and TAB are in the class
#: because the commonest way a phenotype/identifier table reaches this repository
#: is as a CSV or TSV, where the delimiter IS the separator: a two-column
#: keyword/value row was invisible while only `[:=#]` was accepted, which made the
#: keyed rules blind to the single most likely carrier format. (As elsewhere in
#: this file, no worked example is written out — a live one would make the module
#: trip its own scanner, which is correct behaviour and a permanent false alarm.)
_KEYED_SEP: Final[bytes] = rb"[ \t]{0,4}[:=#,\t][ \t]{0,4}"

_MRN: Final = re.compile(
    rb"(?i:\bMRN\b"
    rb"|\bmedical[ _-]record[ _-](?:number|no\.?|id)\b"
    rb"|\bNHS[ _-]number\b"
    rb"|\bpatient[ _-]?id\b"
    rb"|\bhospital[ _-]?(?:number|no\.?|id)\b)"
    + _KEYED_SEP
    + rb"(?![A-Za-z]+\b)[A-Za-z0-9][A-Za-z0-9\-]{3,31}\b"
)

_DOB: Final = re.compile(
    rb"(?i:\bDOB\b|\bdate[ _-]of[ _-]birth\b|\bbirth[ _-]?date\b|\bdate[ _-]born\b)"
    + _KEYED_SEP
    + rb"[0-9]{1,4}[/\-.][0-9]{1,2}[/\-.][0-9]{1,4}"
)

# A bare ISO date is WARN and nothing more: every provenance manifest, changelog
# and lockfile in the repository is full of them. Escalating this to `fail` would
# make the audit unrunnable within a day.
# A trailing \b would miss the commonest form of all, `2026-01-01T00:00:00Z`,
# because `1` and `T` are both word characters. The boundaries are digit-aware
# instead, so a date embedded in a timestamp still matches.
_ISO_DATE_BARE: Final = re.compile(
    rb"(?<![\d-])(?:19|20)\d{2}-(?:0[1-9]|1[0-2])-(?:0[1-9]|[12]\d|3[01])(?!\d)"
)

_PERSON_NAME_KEYED: Final = re.compile(
    rb"(?i:(?:patient|proband|child|subject|mother|father)"
    rb"[ _-]?(?:surname|first[ _-]?name|last[ _-]?name|name))"
    + _KEYED_SEP
    + rb"[A-Z][A-Za-z'\-]{1,30}\b"
)

# ---------------------------------------------------------------------------
# Redaction-only patterns
#
# These are NOT in RULES and therefore never run in the audit's content scan.
# They exist because :func:`mva.privacy.redact.redact_text` was blind to this
# project's OWN canonical renderings — the exact strings other modules
# interpolate into log lines and exception messages — while being unusable as
# audit rules: a bare `0/1` token occurs in every ratio, fraction and code
# comment in the tree, so scanning for it would bury the report.
# ---------------------------------------------------------------------------

#: `GRCh38:chr12:9999999:C:T` — the canonical variant_id this project builds and
#: then puts in messages ("no evidence for {variant_id}"). Position plus alleles
#: is an identifying genotype observation.
_VARIANT_ID: Final = re.compile(
    rb"GRCh3[78]:chr[^:\s]{1,32}:\d{1,12}:[ACGTN*.]{1,64}:[ACGTN*.]{1,64}"
)

#: A bare called genotype token (`0/1`, `1|0`, `./.`). The lookarounds keep it to
#: a single allele on each side, so `50/50` and `1/2/2020` do not match.
_GENOTYPE_TOKEN: Final = re.compile(rb"(?<![\d/|.])[0-9.][/|][0-9.](?![\d/|])")


RULES: Final[tuple[Rule, ...]] = (
    Rule(
        rule_id="vcf_header",
        pattern=_VCF_HEADER,
        severity="fail",
        description="VCF format declaration at the start of a line.",
        false_positive_risk=(
            "Near zero. Only a real VCF (or a document quoting one at column 0) "
            "carries this token line-anchored."
        ),
    ),
    Rule(
        rule_id="vcf_chrom_line",
        pattern=_VCF_CHROM_LINE,
        severity="fail",
        description="VCF column header row; everything after it is genotype data.",
        false_positive_risk="Near zero. Requires the exact eight-column tab-separated header.",
    ),
    Rule(
        rule_id="vcf_data_line",
        pattern=_VCF_DATA_LINE,
        severity="warn",
        description="Tab-structured CHROM/POS/ID/REF/ALT record line.",
        false_positive_risk=(
            "HIGH in isolation, and higher since CHROM was widened to any contig "
            "name: public coordinate tables (ClinVar exports, gnomAD extracts, our "
            "own knowledge/public/*.tsv) match legitimately, and so does any TSV "
            "whose 4th and 5th columns are nucleotide alphabets. This is why the "
            "rule is `warn`; mva.privacy.audit promotes it to `fail` only when the "
            "same file also matches vcf_header or genotype_field."
        ),
    ),
    Rule(
        rule_id="genotype_field",
        pattern=_GENOTYPE_FIELD,
        severity="fail",
        description="A FORMAT/GT column paired with an actual called genotype.",
        false_positive_risk=(
            "Low. Requires literal tabs around a GT-led FORMAT string followed by "
            "a called allele. Accepting alphanumeric FORMAT keys (F1R2, F2R1) and "
            "a single haploid allele widens it slightly: a TSV column literally "
            "named GT followed by a numeric column now matches. That is a rare "
            "shape, and the missed somatic and haploid call sets were not."
        ),
    ),
    Rule(
        rule_id="fastq_record",
        pattern=_FASTQ_RECORD,
        severity="fail",
        description="Four-line FASTQ record with >=25nt of homogeneous nucleotide alphabet.",
        false_positive_risk=(
            "Low. A 25nt run drawn only from ACGTN, sandwiched between an '@' "
            "header and a '+' separator, does not occur in prose or code."
        ),
    ),
    Rule(
        rule_id="sam_rg_sample",
        pattern=_SAM_RG_SAMPLE,
        severity="fail",
        description="SAM/BAM @RG read-group line carrying an SM: sample identifier.",
        false_positive_risk=(
            "Very low, and the payload is severe: SM: values are real hospital "
            "sample or accession IDs written by the sequencing centre. Accepting a "
            "space as the field separator adds only lines that begin `@RG ` at "
            "column 0 and carry an SM: token."
        ),
    ),
    Rule(
        rule_id="hpo_term",
        pattern=_HPO_TERM,
        severity="fail",
        description="Human Phenotype Ontology term identifier.",
        false_positive_risk=(
            "HIGH in isolation, for both the `HP:` and the OBO `HP_` rendering. "
            "The public ontology tables under knowledge/public/, "
            "test fixtures and code constants all contain HPO IDs legitimately. "
            "mva.privacy.audit only fails a NON-allowlisted file carrying >=3 "
            "DISTINCT terms, on the reasoning that three co-occurring phenotypes "
            "start to be a clinical profile rather than a reference."
        ),
    ),
    Rule(
        rule_id="mrn",
        pattern=_MRN,
        severity="fail",
        description="Keyword-anchored medical record / patient / NHS identifier.",
        false_positive_risk=(
            "Low, and deliberately so. The rule is anchored on the keyword and "
            "MUST NOT match a bare \\d{7,10} run, which would collide with every "
            "genomic POS and make the scanner unusable. `,` and TAB are accepted "
            "as separators so a two-column CSV/TSV row matches; the cost is that "
            "the same keyword followed by a comma and a value in prose matches "
            "too, which is the correct trade."
        ),
    ),
    Rule(
        rule_id="dob",
        pattern=_DOB,
        severity="fail",
        description="Keyword-anchored date of birth.",
        false_positive_risk=(
            "Low; requires an explicit DOB-style keyword, a separator (now "
            "including `,` and TAB, so a CSV/TSV column pair matches) and a "
            "day/month/year triple."
        ),
    ),
    Rule(
        rule_id="iso_date_bare",
        pattern=_ISO_DATE_BARE,
        severity="warn",
        description="An unkeyed ISO-8601 calendar date.",
        false_positive_risk=(
            "Very high by construction: provenance timestamps, lockfiles and "
            "changelogs are made of these. Kept as `warn` for visibility only, and "
            "exempted from log redaction (see REDACTION_EXEMPT_RULES)."
        ),
    ),
    Rule(
        rule_id="person_name_keyed",
        pattern=_PERSON_NAME_KEYED,
        severity="fail",
        description="A relationship/role keyword bound to a capitalised personal name.",
        false_positive_risk=(
            "Low. Requires keyword + explicit separator (including `,` and TAB, so "
            "a CSV/TSV column pair matches) + a capitalised token, so 'patient "
            "name field' in prose does not match."
        ),
    ),
    Rule(
        rule_id="fasta_record",
        pattern=_FASTA_RECORD,
        severity="fail",
        description="FASTA description line followed by >=40nt of nucleotide alphabet.",
        false_positive_risk=(
            "Low. A `>`-led line at column 0 followed by a 40-character run drawn "
            "only from ACGTNU does not occur in prose, code or markdown. A PUBLIC "
            "reference FASTA matches too, and that is intended: a reference slice "
            "belongs in the workspace or a cache, never committed to this repo."
        ),
    ),
    Rule(
        rule_id="plink_ped_line",
        pattern=_PLINK_PED_LINE,
        severity="fail",
        description="PLINK .ped body line: six pedigree fields followed by allele pairs.",
        false_positive_risk=(
            "Low. Requires six whitespace-separated fields with a PLINK sex code "
            "and phenotype code in positions 5 and 6, then at least six "
            "single-character allele tokens. The payload is a family structure "
            "plus a genotype matrix, which is as disclosive as a VCF."
        ),
    ),
)

_RULES_BY_ID: Final[dict[str, Rule]] = {rule.rule_id: rule for rule in RULES}


# ---------------------------------------------------------------------------
# Positive controls
#
# A detector that has silently stopped detecting reports the same "clean" as a
# repository that is genuinely clean, and the two are indistinguishable from the
# outside. That is not hypothetical here: this project's audit passed over an IRB
# protocol number, a specimen accession, a sequencing flowcell ID and several
# derived callset statistics, because no rule was looking for them. It reported
# success and was believed.
#
# So every rule carries a specimen it MUST match. The battery is re-proved against
# these on every audit, and a rule that fails to fire on its own canary aborts the
# run instead of contributing a passing check.
#
# The specimens are ASSEMBLED rather than written as literals wherever a literal
# would make this file match its own battery. ``hpo_term`` and ``iso_date_bare``
# are unanchored and would otherwise turn the positive controls into findings — the
# same self-matching hazard documented above the patterns, arriving from the other
# direction.
# ---------------------------------------------------------------------------


def detector_canaries() -> dict[str, bytes]:
    """One specimen per rule, each of which that rule must match.

    Fabricated throughout. ``CANARY`` and ``Canaryson`` are not a real sample name
    or a real person, the coordinates are not a real locus, and the read is four
    bases repeated -- a specimen has to be recognisable to the rule, not plausible
    as data.
    """
    return {
        "vcf_header": b"##fileformat=VCFv4.2\n",
        "vcf_chrom_line": b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        "vcf_data_line": b"chrCANARY\t1\t.\tA\tG\t.\tPASS\t.\n",
        "genotype_field": b"canary\tGT:AD:DP\t0/1\n",
        "fastq_record": b"@CANARY\n" + b"ACGT" * 10 + b"\n+\n" + b"I" * 40 + b"\n",
        "sam_rg_sample": b"@RG\tID:1\tSM:CANARY\n",
        # assembled: a literal would make this file match hpo_term
        "hpo_term": b" HP" + b":" + b"0000118 ",
        "mrn": b"MRN: 123456789\n",
        "dob": b"DOB: 01/02/1990\n",
        # assembled: a literal would make this file match iso_date_bare
        "iso_date_bare": b" 2026" + b"-01-02 ",
        "person_name_keyed": b"patient name: Canaryson\n",
        "fasta_record": b">canary\n" + b"ACGT" * 15 + b"\n",
        "plink_ped_line": b"F1\tI1\tP1\tM1\t1\t-9\tA\tG\tA\tG\tA\tG\n",
    }


def rule_by_id(rule_id: str) -> Rule:
    """Look up a rule, failing loudly on an unknown ID."""
    try:
        return _RULES_BY_ID[rule_id]
    except KeyError:
        known = ", ".join(sorted(_RULES_BY_ID))
        msg = f"Unknown privacy rule {rule_id!r}. Known rules: {known}."
        raise KeyError(msg) from None


#: Rules whose matches are NOT stripped from log records. Redacting every ISO date
#: would erase the provenance timestamps that make a log readable, for a rule that
#: is warn-only precisely because it is dominated by false positives.
REDACTION_EXEMPT_RULES: Final[frozenset[str]] = frozenset({"iso_date_bare"})

#: Rules that redact but never appear in an audit finding.
#:
#: The split exists because the two jobs have opposite cost functions. An audit
#: rule that fires on a legitimate file costs a red build and eventually gets the
#: audit switched off, so audit rules must be specific. A redaction rule that
#: fires spuriously costs one unreadable token in a log line, so redaction rules
#: should be greedy. These two are unusable as audit rules — every fraction and
#: ratio in the tree is a `0/1` — and mandatory as redaction rules, because they
#: are the canonical forms THIS project builds and then interpolates into log
#: messages and exception text.
REDACTION_ONLY_RULES: Final[tuple[Rule, ...]] = (
    Rule(
        rule_id="variant_id",
        pattern=_VARIANT_ID,
        severity="fail",
        description="This project's canonical build:chrom:pos:ref:alt variant identifier.",
        false_positive_risk=(
            "Not applicable: redaction-only. A build-qualified coordinate with "
            "both alleles is an identifying observation wherever it appears."
        ),
    ),
    Rule(
        rule_id="genotype_token",
        pattern=_GENOTYPE_TOKEN,
        severity="fail",
        description="A bare called genotype token such as 0/1, 1|0 or ./.",
        false_positive_risk=(
            "HIGH, which is why it is redaction-only and never an audit rule: a "
            "prose fraction like 1/2 matches. In a log line that is an acceptable "
            "loss; in the audit report it would be unusable noise."
        ),
    ),
)

#: Rules applied by :func:`mva.privacy.redact.redact_text`.
REDACTION_RULES: Final[tuple[Rule, ...]] = (
    tuple(rule for rule in RULES if rule.rule_id not in REDACTION_EXEMPT_RULES)
    + REDACTION_ONLY_RULES
)

#: Rules that are threshold- or context-gated by the audit rather than absolute.
CONTEXT_GATED_RULES: Final[frozenset[str]] = frozenset({"vcf_data_line", "hpo_term"})

#: Distinct HPO terms in one non-allowlisted file before the finding becomes `fail`.
HPO_DISTINCT_FAIL_THRESHOLD: Final[int] = 3


# ---------------------------------------------------------------------------
# Magic-byte sniffing (NOT regex)
# ---------------------------------------------------------------------------

#: BAM/CRAM/BGZF are containers. A regex over compressed bytes finds nothing, and
#: an extension check is trivially defeated by a rename, so the container is
#: identified from its header bytes.
_CRAM_MAGIC: Final = b"CRAM"
_BGZF_MAGIC: Final = b"\x1f\x8b\x08\x04"
_BGZF_EXTRA: Final = b"BC\x02\x00"
_BAM_MAGIC: Final = b"BAM\x01"

#: Severity for each sniffed container type. BGZF alone is `warn` because every
#: bgzipped PUBLIC resource (reference FASTA, GTF, dbSNP) has exactly this header.
BINARY_SEVERITY: Final[dict[str, Severity]] = {
    "cram": "fail",
    "bam": "fail",
    "bgzf": "warn",
}


def _first_bgzf_block_is_bam(head: bytes) -> bool:
    """Inflate the first BGZF block and look for the BAM magic.

    Uses the BSIZE field from the BGZF extra subfield to bound the block, then
    inflates with a gzip-aware window. Any zlib failure is treated as "not BAM":
    a truncated head is a normal condition when only the first few KiB were read.
    """
    if len(head) < 18:
        return False
    bsize = int.from_bytes(head[16:18], "little") + 1
    block = head[:bsize] if 0 < bsize <= len(head) else head
    try:
        decompressed = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(block, 64)
    except zlib.error:
        return False
    return decompressed[:4] == _BAM_MAGIC


def sniff_binary(head: bytes) -> str | None:
    """Identify a genomic container from its leading bytes.

    Returns ``"cram"``, ``"bam"``, ``"bgzf"`` or ``None``. Never returns any of
    the bytes it inspected.
    """
    if head[:4] == _CRAM_MAGIC:
        return "cram"
    if head[:4] == _BGZF_MAGIC and head[12:16] == _BGZF_EXTRA:
        return "bam" if _first_bgzf_block_is_bam(head) else "bgzf"
    return None


# ---------------------------------------------------------------------------
# Safe file reading
# ---------------------------------------------------------------------------

#: Files larger than this are not content-scanned. A 60 GiB CRAM is identified by
#: its magic bytes and its extension; reading it would be pointless and slow.
MAX_SCAN_BYTES: Final[int] = 8 * 1024 * 1024

#: Enough for every container header we sniff, and for a BGZF first block.
HEAD_BYTES: Final[int] = 65536


def read_capped(path: Path, limit: int = MAX_SCAN_BYTES) -> bytes:
    """Read at most ``limit`` bytes. Always binary; never decodes."""
    with path.open("rb") as handle:
        return handle.read(limit)


#: Bound on how many gzip members :func:`gunzip_capped` will walk. A BGZF file is
#: a chain of ~64 KiB members; the header and the first records — which is all the
#: rules need — are in the first few.
_MAX_GZIP_MEMBERS: Final[int] = 64


def gunzip_capped(data: bytes, limit: int = MAX_SCAN_BYTES) -> bytes | None:
    """Inflate a gzip/BGZF byte string, capped, or return ``None``.

    Every rule in this module matches plaintext, so a gzipped VCF — the form
    variant callers actually emit — matched nothing at all: not the header, not
    the ``#CHROM`` line, not a single genotype. Magic-byte sniffing caught BGZF
    and rated it ``warn`` (every bgzipped public resource has that header), and
    plain gzip was not caught at all.

    Returns ``None`` for non-gzip input and for a stream that fails to inflate; a
    truncated or corrupt member is a normal condition when only a capped prefix
    was read, and is not something to raise about. Output is bounded by ``limit``
    and by :data:`_MAX_GZIP_MEMBERS`, so a decompression bomb cannot exhaust
    memory here.
    """
    if data[:2] != b"\x1f\x8b":
        return None
    out = bytearray()
    remaining = data
    for _ in range(_MAX_GZIP_MEMBERS):
        obj = zlib.decompressobj(16 + zlib.MAX_WBITS)
        try:
            out += obj.decompress(remaining, max(0, limit - len(out)))
        except zlib.error:
            break
        if len(out) >= limit or not obj.eof or not obj.unused_data:
            break
        remaining = obj.unused_data
    return bytes(out) if out else None


def decode_scrubbed(data: bytes, *, path: Path) -> str:
    """Strict-decode bytes, converting a decode failure into a content-free error.

    ``UnicodeDecodeError`` embeds the offending bytes in ``str(exc)``. Letting one
    escape would defeat every other control in this package, because a traceback
    travels to the terminal, the CI log and the model context. So the original is
    raised ``from None`` — chaining it with ``from exc`` would re-attach the very
    message we are suppressing to the ``__cause__`` of the new one.
    """
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError as exc:
        msg = (
            f"Undecodable bytes in {path.as_posix()} at offset {exc.start} "
            f"(run of {exc.end - exc.start} bytes). Content withheld under GP-41."
        )
        raise PrivacyViolationError(msg) from None


def decode_lossy(data: bytes) -> str:
    """Decode for scanning purposes only. Never raises; never round-trips."""
    return data.decode("utf-8", errors="replace")
