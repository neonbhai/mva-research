#!/usr/bin/env bash
# Install SnpEff + a GRCh38 gene model for the real consequence adapter.
#
#   tools/setup/install_snpeff.sh [install_root]
#
# Default install_root: ../mva-resources/snpeff (deliberately OUTSIDE the repo --
# the database is ~790 MB compressed, ~1.6 GB unpacked, and must never be
# committed).
#
# WHAT THIS IS, AND WHY IT IS NOT IN src/mva/annotation/
# ------------------------------------------------------
# This is the "separate, public-only acquisition step" the adapter README requires
# (PRIV-05). It downloads PUBLIC reference data only -- a SnpEff release, an
# Ensembl-derived gene model and the NCBI MANE summary -- and it runs before, and
# separately from, any patient data being loaded. It sends dataset names and URLs
# and nothing else. No proband coordinate, genotype, sample ID or pedigree detail
# is transmitted, by this script or by anything downstream of it: annotation
# itself runs offline against what this script leaves on disk.
#
# Everything is idempotent: re-running skips artifacts already present and
# re-verifies the ones that are.
set -euo pipefail

INSTALL_ROOT="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/../mva-resources/snpeff}"
mkdir -p "$INSTALL_ROOT"
INSTALL_ROOT="$(cd "$INSTALL_ROOT" && pwd)"

# Pinned so that a re-run a month from now reproduces today's answer rather than
# silently upgrading the gene model underneath an existing report (GP-30).
SNPEFF_CORE_URL="https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip"
SNPEFF_DB_VERSION="v5_4"
SNPEFF_GENOME="GRCh38.115"
SNPEFF_DB_URL="https://snpeff-public.s3.amazonaws.com/databases/${SNPEFF_DB_VERSION}/snpEff_${SNPEFF_DB_VERSION}_${SNPEFF_GENOME}.zip"
MANE_RELEASE="1.5"
MANE_URL="https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_${MANE_RELEASE}/MANE.GRCh38.v${MANE_RELEASE}.summary.txt.gz"
JDK_FORMULA="openjdk@21"
DOWNLOAD_PARTS=8

# ---------------------------------------------------------------- EXPECTED PINS
#
# sha256 of the INSTALLED (unpacked) artifacts, not of the archives: these are the
# exact bytes SnpEffConsequenceAdapter verifies before it annotates anything, and
# the adapter refuses to run without them.
#
# Why this matters more here than for a normal dependency: the core archive is
# published as `snpEff_latest_core.zip`. The name is a promise that the bytes
# behind that URL WILL change. SnpEff also rebuilds a genome under an unchanged
# name when the upstream Ensembl annotation is corrected. Without a pin, two runs
# months apart can produce different transcripts and different HGVS while
# stamping the identical `5.4c/GRCh38.115` provenance onto every EvidenceItem --
# which is exactly the failure that recording a version is supposed to prevent.
#
# Measured 2026-08-28 against SnpEff 5.4c (build 2026-02-23) and the Ensembl-115
# database (retrieval_date 2025-12-01 per snpEff.config).
#
# Override any of these with an environment variable to install a different
# pinned release; set to the empty string to record-without-enforcing on a first
# install of a version nobody has pinned yet.
EXPECT_JAR="${EXPECT_JAR-5e8f75cbf908a33c6fb2e65c81e66fe31236cb21bb0195541c4703d8202c22b3}"
EXPECT_CONFIG="${EXPECT_CONFIG-777768cee885c91c7396ca716c03abb254ef4cde129de2a10513e2844895df5a}"
EXPECT_PREDICTOR="${EXPECT_PREDICTOR-f318b79b26ced7e8a44ac18e749ade83bef642e74387163d71885127d85b357a}"
EXPECT_MANE="${EXPECT_MANE-d10ace2720681a3b2e0eefd9da4f551274a6b4141ac9bfd6a2565dfb6e9ad55c}"

say() { printf '\n==> %s\n' "$*"; }
die() { printf '\nFAILED: %s\n' "$*" >&2; exit 1; }

sha256_of() { shasum -a 256 "$1" | awk '{print $1}'; }

# check_pin <label> <path> <expected>
# Empty expected records instead of enforcing, and says so loudly.
check_pin() {
  local label="$1" path="$2" expected="$3" actual
  [ -f "$path" ] || die "$label missing at $path"
  actual="$(sha256_of "$path")"
  if [ -z "$expected" ]; then
    printf '  %-22s %s  (UNPINNED - recorded, not enforced)\n' "$label" "$actual"
    return 0
  fi
  if [ "$actual" != "$expected" ]; then
    die "$label sha256 mismatch.
    expected $expected
    found    $actual
  The bytes behind these URLs change. Review the difference, then update the
  EXPECT_* pins in this script AND the pins passed to SnpEffConsequenceAdapter.
  Refusing to leave an installation whose provenance string would be misleading."
  fi
  printf '  %-22s %s  OK\n' "$label" "$actual"
}

# fetch_parallel <url> <output> <parts>
# Resumable, ranged, parallel download. Ranged parts are used because these are
# ~800 MB single-connection transfers that a laptop's link will otherwise stall
# on; -C - makes an interrupted run resume rather than restart.
fetch_parallel() {
  local url="$1" out="$2" parts="$3" total part_size i start end pids=()
  [ -f "$out" ] && return 0
  total="$(curl -fsSLI "$url" | awk 'tolower($1) ~ /^content-length:/ {print $2}' | tr -d '\r' | tail -1)"
  if [ -z "$total" ] || [ "$total" -le 0 ] 2>/dev/null; then
    say "content-length unavailable; falling back to a single stream"
    curl -fL --retry 8 --retry-all-errors -C - -o "$out.tmp" "$url"
    mv "$out.tmp" "$out"
    return 0
  fi
  part_size=$(( (total + parts - 1) / parts ))
  say "downloading $((total / 1024 / 1024)) MB in $parts parallel parts"
  for ((i = 0; i < parts; i++)); do
    start=$((i * part_size))
    end=$((start + part_size - 1))
    [ "$end" -ge "$total" ] && end=$((total - 1))
    curl -sS -L --retry 8 --retry-all-errors -C - -r "${start}-${end}" -o "$out.part$i" "$url" &
    pids+=($!)
  done
  for pid in "${pids[@]}"; do wait "$pid" || die "a download part failed"; done

  # Every part must be exactly its requested length before anything is
  # concatenated. Joining a part that is still being written produces a zip that
  # unpacks *partially* -- which is far worse than one that fails outright.
  for ((i = 0; i < parts; i++)); do
    start=$((i * part_size))
    end=$((start + part_size - 1))
    [ "$end" -ge "$total" ] && end=$((total - 1))
    local want=$((end - start + 1))
    local got
    got="$(wc -c < "$out.part$i" | tr -d ' ')"
    [ "$got" = "$want" ] || die "part $i is $got bytes, expected $want; re-run to resume"
  done

  cat "$out".part* > "$out.tmp"
  local assembled
  assembled="$(wc -c < "$out.tmp" | tr -d ' ')"
  [ "$assembled" = "$total" ] || die "assembled $assembled bytes, expected $total"
  mv "$out.tmp" "$out"
  rm -f "$out".part*
}

# --------------------------------------------------------------------- 1. a JRE
# macOS ships a /usr/bin/java *stub* that reports "Unable to locate a Java
# Runtime". SnpEff needs a real one; Homebrew's openjdk is keg-only, so it is
# symlinked into the install root rather than onto PATH (non-destructive: nothing
# outside INSTALL_ROOT and Homebrew's own Cellar is touched).
#
# The adapter never resolves `java` through PATH for exactly this reason -- see
# resolve_java_binary() -- so the path below is the one to hand it.
say "Java runtime"
if [ ! -x "$INSTALL_ROOT/jdk/bin/java" ]; then
  if [ -n "${JAVA_HOME:-}" ] && [ -x "$JAVA_HOME/bin/java" ]; then
    ln -sfn "$JAVA_HOME" "$INSTALL_ROOT/jdk"
  else
    if ! brew --prefix "$JDK_FORMULA" >/dev/null 2>&1; then
      brew install "$JDK_FORMULA"
    fi
    ln -sfn "$(brew --prefix "$JDK_FORMULA")/libexec/openjdk.jdk/Contents/Home" "$INSTALL_ROOT/jdk"
  fi
fi
JAVA="$INSTALL_ROOT/jdk/bin/java"
"$JAVA" -version 2>&1 | head -1

# ------------------------------------------------------------------ 2. the tool
say "SnpEff core"
if [ ! -f "$INSTALL_ROOT/snpEff/snpEff.jar" ]; then
  fetch_parallel "$SNPEFF_CORE_URL" "$INSTALL_ROOT/snpEff_core.zip" 4
  unzip -q -o "$INSTALL_ROOT/snpEff_core.zip" -d "$INSTALL_ROOT"
fi
"$JAVA" -jar "$INSTALL_ROOT/snpEff/snpEff.jar" -version

# -------------------------------------------------------------- 3. the gene model
# ~790 MB compressed. The database VERSION must match the jar's: SnpEff refuses a
# predictor built by an incompatible release, and the genome name must appear in
# snpEff.config, which is why the two are pinned together above.
say "SnpEff database ${SNPEFF_GENOME} (~790 MB download)"
if [ ! -f "$INSTALL_ROOT/data/${SNPEFF_GENOME}/snpEffectPredictor.bin" ]; then
  fetch_parallel "$SNPEFF_DB_URL" "$INSTALL_ROOT/snpEff_db.zip" "$DOWNLOAD_PARTS"
  unzip -q -o "$INSTALL_ROOT/snpEff_db.zip" -d "$INSTALL_ROOT"
fi
grep -q "^${SNPEFF_GENOME}.genome" "$INSTALL_ROOT/snpEff/snpEff.config" \
  || die "${SNPEFF_GENOME} is not declared in snpEff.config; SnpEff resolves a genome
  through the config, so it would report the database as missing and -- with
  -nodownload set -- refuse to fetch it."

# ------------------------------------------------------------- 4. the MANE flags
# SnpEff's ANN field carries no MANE-Select flag. Without this file the adapter
# leaves is_mane_select unset for every transcript rather than guessing.
say "MANE summary v${MANE_RELEASE}"
mkdir -p "$INSTALL_ROOT/mane"
MANE_FILE="$INSTALL_ROOT/mane/MANE.GRCh38.v${MANE_RELEASE}.summary.txt.gz"
[ -f "$MANE_FILE" ] || curl -fL --retry 8 --retry-all-errors -C - -o "$MANE_FILE" "$MANE_URL"

# --------------------------------------------------------------- 5. verify pins
say "Artifact pins (sha256 of the bytes the adapter will annotate against)"
PREDICTOR="$INSTALL_ROOT/data/${SNPEFF_GENOME}/snpEffectPredictor.bin"
CONFIG="$INSTALL_ROOT/snpEff/snpEff.config"
JAR="$INSTALL_ROOT/snpEff/snpEff.jar"
check_pin "snpEff.jar"          "$JAR"       "$EXPECT_JAR"
check_pin "snpEff.config"       "$CONFIG"    "$EXPECT_CONFIG"
check_pin "snpEffectPredictor"  "$PREDICTOR" "$EXPECT_PREDICTOR"
check_pin "MANE summary"        "$MANE_FILE" "$EXPECT_MANE"

PINS_FILE="$INSTALL_ROOT/snpeff_pins.json"
cat > "$PINS_FILE" <<PINS
{
  "snpeff_release": "$("$JAVA" -jar "$JAR" -version 2>&1 | head -1 | awk '{print $2}')",
  "genome_database": "$SNPEFF_GENOME",
  "jar": "$(sha256_of "$JAR")",
  "config": "$(sha256_of "$CONFIG")",
  "predictor": "$(sha256_of "$PREDICTOR")",
  "mane_summary": "$(sha256_of "$MANE_FILE")"
}
PINS

# ---------------------------------------------------------------- 6. prove it works
# A real, PUBLIC, pathogenic ClinVar record -- NOT patient data -- annotated with
# the same offline flags the adapter uses, so an install that would fail at run
# time fails here instead.
#
#   chr15:40165186 C>T  BUB1B  ClinVar Pathogenic, ALLELEID 1361261,
#                       MC=SO:0001587|nonsense, MVA syndrome 1.
#
# Asserted against ClinVar's own molecular consequence: a smoke test that only
# checked "some ANN came back" would pass against a database annotating the wrong
# assembly. Bare `15` is used because this database is Ensembl-named.
say "Offline smoke test (public ClinVar variant, not patient data)"
SMOKE_OUT="$(printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n15\t40165186\tsmoke\tC\tT\t.\t.\t.\n' \
  | "$JAVA" -Xmx6g -jar "$JAR" ann \
      -noLog -nodownload -noStats \
      -config "$CONFIG" -dataDir "$INSTALL_ROOT/data" \
      -i vcf -o vcf "$SNPEFF_GENOME" 2>/dev/null | grep -v '^#')"
printf '%s\n' "$SMOKE_OUT" | cut -c1-200
case "$SMOKE_OUT" in
  *ERROR_CHROMOSOME_NOT_FOUND*) die "contig naming mismatch: the database did not recognise '15'" ;;
esac
case "$SMOKE_OUT" in
  *"stop_gained|HIGH|BUB1B"*) ;;
  *) die "expected stop_gained/HIGH/BUB1B (ClinVar calls this variant a nonsense change)
  and did not get it. The database may be built against the wrong assembly." ;;
esac
case "$SMOKE_OUT" in
  *"p.Gln57*"*) ;;
  *) die "expected HGVS.p p.Gln57* and did not get it" ;;
esac

# The mitochondrion is the case where contig mapping is load-bearing: SnpEff
# strips a leading 'chr' itself, so chr15 forgives a mismatch, but 'chrM' becomes
# 'M' and this database calls it 'MT'. Verified rather than assumed.
MT_OUT="$(printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nMT\t8993\tmt\tT\tG\t.\t.\t.\n' \
  | "$JAVA" -Xmx6g -jar "$JAR" ann -noLog -nodownload -noStats \
      -config "$CONFIG" -dataDir "$INSTALL_ROOT/data" -i vcf -o vcf "$SNPEFF_GENOME" 2>/dev/null | grep -v '^#')"
case "$MT_OUT" in
  *MT-ATP6*) ;;
  *) die "the database did not recognise 'MT'; ContigStyle.ENSEMBL would silently
  return no annotation for every mitochondrial variant." ;;
esac

# ------------------------------------------------------------------- 7. provenance
say "Installed"
cat <<SUMMARY
  install root : $INSTALL_ROOT
  java         : $("$JAVA" -version 2>&1 | head -1)
  snpEff       : $("$JAVA" -jar "$JAR" -version 2>&1 | head -1 | tr '\t' ' ')
  database     : $SNPEFF_GENOME ($(du -sh "$INSTALL_ROOT/data/$SNPEFF_GENOME" | cut -f1))
  mane         : $(basename "$MANE_FILE")
  pins         : $PINS_FILE

Wire it up (composition root):

    from mva.annotation.snpeff_local import SnpEffArtifactPins, SnpEffConsequenceAdapter

    SnpEffConsequenceAdapter(
        jar_path=Path("$JAR"),
        data_dir=Path("$INSTALL_ROOT/data"),
        genome_database="$SNPEFF_GENOME",
        mane_summary=Path("$MANE_FILE"),
        pins=SnpEffArtifactPins(
            jar="$(sha256_of "$JAR")",
            config="$(sha256_of "$CONFIG")",
            predictor="$(sha256_of "$PREDICTOR")",
            mane_summary="$(sha256_of "$MANE_FILE")",
        ),
        # java_binary defaults to \$JAVA_HOME/bin/java; pass it explicitly to use
        # this install's JDK regardless of the environment:
        java_binary=Path("$JAVA"),
    )
SUMMARY
