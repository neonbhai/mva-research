#!/bin/sh
# mva_synthetic=true
#
# A stand-in for `java -jar snpEff.jar` used ONLY by tests/unit/test_consequence_adapter.py.
#
# It exists so the adapter's real behaviour -- argv construction, stdin plumbing,
# ANN parsing, transcript preservation, ordering, omission, timeout handling and
# the no-leak rule on failure -- is testable in under a second, with no JRE and no
# 600 MB database on the machine. It is NOT an annotator: it replays canned ANN
# strings from canned_snpeff_ann.tsv next to it. The adapter under test is real;
# only this counterparty is faked, and the real tool is exercised separately by
# the integration test that skips when SnpEff is not installed.
#
# Behaviour is selected by the LAST argument (the genome database name), matching
# how SnpEff itself is invoked:
#   -version anywhere  -> print a SnpEff-shaped version line, exit 0
#   ...    SLEEPDB     -> sleep past any sane timeout, exit 0
#   ...    FAILDB      -> echo stdin (i.e. the input VCF) to stderr, exit 3
#   ...    <anything>  -> annotate stdin from the canned table
#
# Four modes exit ZERO while returning output the adapter must refuse. They exist
# because that is the dangerous case: a non-zero exit is already handled, whereas
# a successful-looking run that returns a partial answer is read downstream as
# "these variants have no consequence" and deletes them from gene grouping.
#   ...    TRUNCDB     -> emit a data row cut short mid-record, exit 0
#   ...    UNKNOWNIDDB -> emit a record under an ID that was never sent, exit 0
#   ...    DUPIDDB     -> emit the first record twice, exit 0
#   ...    MISSINGDB   -> drop the last record, exit 0
set -e
here=$(dirname "$0")

for arg in "$@"; do
  if [ "$arg" = "-version" ]; then
    # Real SnpEff prints its version on stderr, not stdout. Faithfully so here,
    # because the adapter's probe has to cope with exactly that.
    printf 'SnpEff\t5.2\t2024-09-24\n' >&2
    exit 0
  fi
done

db=""
for arg in "$@"; do db="$arg"; done

case "$db" in
  SLEEPDB)
    sleep 30
    exit 0
    ;;
  FAILDB)
    # Deliberately the worst case: SnpEff really does echo offending VCF records
    # on failure, so the adapter must not put this stream into its exception.
    cat >&2
    exit 3
    ;;
  TRUNCDB)
    # A record cut short mid-write, which is what a killed pipe or a full disk
    # leaves behind. SnpEff still exits 0.
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    printf '15\t40200000\tmva0\tC\n'
    exit 0
    ;;
  UNKNOWNIDDB)
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    printf '15\t40200000\tnot-an-id-we-sent\tC\tT\t.\t.\t.\n'
    exit 0
    ;;
  DUPIDDB)
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    printf '15\t40200000\tmva0\tC\tT\t.\t.\t.\n'
    printf '15\t40200000\tmva0\tC\tT\t.\t.\t.\n'
    exit 0
    ;;
  MISSINGDB)
    # The run stops after the header: every planned site is unaccounted for.
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    exit 0
    ;;
esac

awk -F '\t' -v canned="$here/canned_snpeff_ann.tsv" '
BEGIN {
  OFS = "\t"
  while ((getline line < canned) > 0) {
    if (line ~ /^#/) continue
    n = split(line, f, "\t")
    if (n < 5 || f[1] == "contig") continue
    ann[f[1] ":" f[2] ":" f[3] ":" f[4]] = f[5]
  }
}
/^##/ { print; next }
/^#CHROM/ {
  print "##SnpEffVersion=\"5.2 (build 2024-09-24 10:11), by Pablo Cingolani\""
  print "##SnpEffCmd=\"SnpEff  -i vcf -o vcf GRCh38.110 /tmp/mva-scratch-9f2a/input.vcf \""
  print
  next
}
{
  key = $1 ":" $2 ":" $4 ":" $5
  if (key in ann) { $8 = "ANN=" ann[key] }
  print
}
'
