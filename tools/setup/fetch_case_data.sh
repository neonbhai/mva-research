#!/bin/bash
# Fetch the GATED PROBAND DATA into the encrypted case volume, and nowhere else.
#
# This is real patient data under WCG IRB protocol #REDACTED-PROTOCOL. No resharing.
# Deletion within 30 days of challenge close (2026-10-24) from ALL environments
# including derived datasets, confirmed by email. See docs/privacy-model.md.
#
#   bash tools/setup/fetch_case_data.sh          # VCF + index + phenotype (303 MB)
#   MVA_FETCH_FASTQ=1 bash tools/setup/fetch_case_data.sh   # also 79 GB of FASTQs
set -euo pipefail

MOUNT="${MVA_MOUNT:-/Volumes/MVACASE}"
DEST="$MOUNT/raw"
REPO="SageBio/mva-hackathon-2026-data"
BASE="https://huggingface.co/datasets/$REPO/resolve/main"

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

get() {
  local f="$1"
  local enc; enc=$(python3 -c 'import urllib.parse,sys; print(urllib.parse.quote(sys.argv[1]))' "$f")
  echo ">> $f"
  curl -fL --progress-bar -C - -H "Authorization: Bearer $TOK" -o "$DEST/$f" "$BASE/$enc"
}

get "WGS_SPECIMEN_FLOWCELL.vcf.gz"
get "WGS_SPECIMEN_FLOWCELL.vcf.gz.tbi"
get "Challenge_Clinical_Phenotype_1.docx"
get "README.md"

if [ "${MVA_FETCH_FASTQ:-0}" = "1" ]; then
  echo ">> Fetching 79 GB of FASTQs (only needed for realignment / SV calling)."
  for L in 1 2 3 4; do for R in 1 2; do
    get "WGS_SPECIMEN_FLOWCELL_S16_L00${L}_R${R}_001.fastq.gz"
  done; done
fi

chmod -R 700 "$DEST"
echo
echo ">> Done. Checksums:"
shasum -a 256 "$DEST"/*.vcf.gz "$DEST"/*.tbi 2>/dev/null || true
echo
echo ">> These files are patient data. They live only inside $MOUNT."
echo ">> Do NOT copy them into the repo. .gitignore denies genomic formats by default."
