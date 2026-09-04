#!/usr/bin/env bash
# USAGE: tsd_annotate_vcf.sh <vcf_in> <TSD_summary.txt> <vcf_out>
#
# Writes INFO/TSD onto every record whose TSD search passed. The value is the
# 5' copy and the 3' copy of the duplication, comma-separated, as
# TSD_Match_v2.sh reports them.
#
# Records are matched to summary rows by ID in awk. The join this replaces
# needed both streams sorted the same way, which depends on sort and join
# agreeing on a locale, and it matched any row containing PASS.
set -euo pipefail

VCF_IN=$1
SUMMARY=$2
VCF_OUT=$3

awk -F'\t' -v OFS='\t' '
  NR == FNR { if ($NF == "PASS") tsd[$1] = toupper($(NF-2)) "," toupper($(NF-1)); next }
  /^#/      { next }
  ($3 in tsd) { print $1, $2, $3, $4, $5, tsd[$3] }
' "${SUMMARY}" "${VCF_IN}" | LC_ALL=C sort -k1,1 -k2,2n > TSD_annotation

n_pass=$(awk -F'\t' '$NF == "PASS"' "${SUMMARY}" | wc -l)
n_annot=$(wc -l < TSD_annotation)
echo "TSD: ${n_pass} PASS in ${SUMMARY}, ${n_annot} matched to records of ${VCF_IN}"

echo '##INFO=<ID=TSD,Number=1,Type=String,Description="Target site duplication sequence passing filters">' > tsd_header
bgzip -f TSD_annotation
tabix -f -s1 -b2 -e2 TSD_annotation.gz
bcftools annotate -a TSD_annotation.gz -h tsd_header \
    -c CHROM,POS,~ID,REF,ALT,INFO/TSD -Ov -o "${VCF_OUT}" "${VCF_IN}"
