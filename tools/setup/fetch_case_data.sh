#!/bin/bash
# Fetch the GATED PROBAND DATA into the encrypted case volume, and nowhere else.
#
# This is real patient data under a gated IRB protocol. No resharing.
# Deletion within 30 days of challenge close (2026-10-24) from ALL environments
# including derived datasets, confirmed by email. See docs/privacy-model.md.
#
# The case file names are themselves patient identifiers -- they carry the
# specimen accession and the sequencer flowcell ID -- so they are NOT written
# down in this public repository. They are DISCOVERED from the gated dataset
# listing at run time, using the same token the download already needs. Anyone
# with dataset access can therefore still run this in one command:
#
#   bash tools/setup/fetch_case_data.sh          # VCF + index + phenotype
#   MVA_FETCH_FASTQ=1 bash tools/setup/fetch_case_data.sh   # also the raw FASTQs (tens of GB)
#
# Discovery can be overridden without ever naming a file in-tree:
#
#   MVA_CASE_MANIFEST=/path/to/manifest   one file name per line, `#` comments
#                                         allowed. Keep it outside the repo --
#                                         on the encrypted volume is best.
#   MVA_CASE_FILES="name1 name2 ..."      explicit names, whitespace-separated.
set -euo pipefail

MOUNT="${MVA_MOUNT:-/Volumes/MVACASE}"
DEST="$MOUNT/raw"
REPO="SageBio/mva-hackathon-2026-data"
BASE="https://huggingface.co/datasets/$REPO/resolve/main"
API="https://huggingface.co/api/datasets/$REPO/tree/main?recursive=1"

# --- refuse to write patient data anywhere but the encrypted volume ---
if ! mount | grep -q " $MOUNT "; then
  echo "FATAL: $MOUNT is not mounted. Run tools/setup/mount_case.sh first." >&2
  echo "Patient data must never land on the plain disk or inside the repo." >&2
  exit 1
fi
case "$(cd "$DEST" && pwd -P)" in
  "$MOUNT"/*) ;;
  *) echo "FATAL: destination $DEST resolves outside $MOUNT. Refusing." >&2; exit 1 ;;
esac

TOK="${HF_TOKEN:-$(cat "$HOME/.cache/huggingface/token" 2>/dev/null || true)}"
[ -n "$TOK" ] || { echo "FATAL: no HF token. Set HF_TOKEN or run 'hf auth login'." >&2; exit 1; }

# --- where the file names come from: the gated listing, never this repo ---
list_dataset() {
  curl -fsSL -H "Authorization: Bearer $TOK" "$API" | python3 -c '
import json, sys
for entry in json.load(sys.stdin):
    if entry.get("type") == "file":
        print(entry["path"])
'
}

case_files() {
  if [ -n "${MVA_CASE_FILES:-}" ]; then
    # shellcheck disable=SC2086  # deliberate word splitting: one name per field
    printf '%s\n' ${MVA_CASE_FILES}
  elif [ -n "${MVA_CASE_MANIFEST:-}" ]; then
    [ -f "$MVA_CASE_MANIFEST" ] || {
      echo "FATAL: MVA_CASE_MANIFEST=$MVA_CASE_MANIFEST does not exist." >&2
      exit 1
    }
    sed -e 's/#.*//' -e 's/[[:space:]]*$//' "$MVA_CASE_MANIFEST" | grep -v '^$'
  else
    echo ">> Listing the gated dataset to discover the case file names." >&2
    list_dataset
  fi
}

FILES=$(case_files)
[ -n "$FILES" ] || {
  echo "FATAL: no case files found. Check dataset access, or set MVA_CASE_MANIFEST." >&2
  exit 1
}

# Select by FILE TYPE, which is not identifying, rather than by name.
matching() {
  printf '%s\n' "$FILES" | grep -E "$1" || true
}

get() {
  local f="$1"
  local enc
  enc=$(python3 -c 'import urllib.parse,sys; print(urllib.parse.quote(sys.argv[1]))' "$f")
  echo ">> $f"
  mkdir -p "$(dirname "$DEST/$f")"
  curl -fL --progress-bar -C - -H "Authorization: Bearer $TOK" -o "$DEST/$f" "$BASE/$enc"
}

CORE=$(matching '\.vcf\.gz$|\.vcf\.gz\.tbi$|\.vcf\.gz\.csi$|\.docx$|(^|/)README\.md$')
[ -n "$CORE" ] || {
  echo "FATAL: the listing contains no VCF/phenotype files. Refusing to guess." >&2
  exit 1
}
while IFS= read -r f; do get "$f"; done <<<"$CORE"

if [ "${MVA_FETCH_FASTQ:-0}" = "1" ]; then
  FASTQ=$(matching '\.fastq\.gz$|\.fq\.gz$')
  if [ -n "$FASTQ" ]; then
    echo ">> Fetching the raw FASTQs, tens of GB (only needed for realignment / SV calling)."
    while IFS= read -r f; do get "$f"; done <<<"$FASTQ"
  else
    echo ">> MVA_FETCH_FASTQ=1 but the listing has no FASTQs. Skipping." >&2
  fi
fi

chmod -R 700 "$DEST"
echo
echo ">> Done. Checksums:"
shasum -a 256 "$DEST"/*.vcf.gz "$DEST"/*.tbi 2>/dev/null || true
echo
echo ">> These files are patient data. They live only inside $MOUNT."
echo ">> Their NAMES are patient data too -- do not paste them into the repo, an"
echo ">> issue, a commit message or a log. .gitignore denies genomic formats by default."
