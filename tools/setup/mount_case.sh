#!/bin/bash
# Create + mount the encrypted case volume for REAL PATIENT DATA.
# Implements docs/privacy-model.md. Run this yourself: it needs a passphrase and sudo.
#
#   bash tools/setup/mount_case.sh
#
# Teardown (cryptographic erasure) is in docs/privacy-model.md. Do not skip it:
# the challenge requires deletion within 30 days of close (2026-10-24).
set -euo pipefail

BUNDLE="${MVA_BUNDLE:-$HOME/private/mvacase.sparsebundle}"
MOUNT="${MVA_MOUNT:-/Volumes/MVACASE}"
SIZE="${MVA_SIZE:-60g}"   # sparse: consumes only what is written. VCF is 300MB.

if [ ! -d "$BUNDLE" ]; then
  mkdir -p "$(dirname "$BUNDLE")"
  echo ">> Creating encrypted sparse bundle at $BUNDLE ($SIZE max, AES-256)."
  echo ">> Choose a passphrase and record it OUT OF BAND -- not the login keychain,"
  echo ">> which is itself backed up and snapshotted."
  hdiutil create -size "$SIZE" -type SPARSEBUNDLE -fs APFS \
    -encryption AES-256 -stdinpass -volname MVACASE "$BUNDLE"
  echo ">> Excluding from Time Machine (needs sudo)."
  sudo tmutil addexclusion -p "$BUNDLE"
  tmutil isexcluded "$BUNDLE"          # expect [Excluded]
else
  echo ">> Reusing existing bundle at $BUNDLE"
fi

if ! mount | grep -q " $MOUNT "; then
  echo ">> Attaching (passphrase again):"
  hdiutil attach -stdinpass -nobrowse -noindex -mountpoint "$MOUNT" "$BUNDLE"
fi

sudo mdutil -i off -d "$MOUNT" || true       # Spotlight off
touch "$MOUNT/.metadata_never_index"
mkdir -p "$MOUNT/case01" "$MOUNT/tmp" "$MOUNT/raw"
chmod 700 "$MOUNT/case01" "$MOUNT/tmp" "$MOUNT/raw"

cat <<EOF

>> Mounted at $MOUNT
>> Now export these in the shell that will run the pipeline:

     export MVA_WORKSPACE=$MOUNT/case01
     export TMPDIR=$MOUNT/tmp

>> Then fetch the proband data:  bash tools/setup/fetch_case_data.sh
EOF
