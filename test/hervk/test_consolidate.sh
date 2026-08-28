#!/usr/bin/env bash
# Regression test for `hervk_reconcile.py consolidate`.
#
# Real CaG data: the ten member records of the five flagged loci, from the
# unfiltered graph-genotyped VCF, with the assembly callset alongside. These
# five cover every shape the consolidator has to handle -- two loci where both
# members describe the same allele, one genuinely triallelic locus, one where a
# member is a tandem duplication whose genotypes must be withheld, and one where
# a member has no usable allele state.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
BIN=../../bin
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

python3 $BIN/hervk_reconcile.py consolidate \
    --genotyped-vcf hervk_graph_fixture.vcf \
    --loci          hervk_loci_fixture.tsv \
    --calls         hervk_calls_fixture.tsv \
    --discovery-vcf hervk_discovery_fixture.vcf \
    --out-vcf       "$TMP/out.vcf" \
    --out-archive   "$TMP/archive.vcf" \
    --report        "$TMP/report.md" 2>"$TMP/err"
cat "$TMP/err"

fail=0
chk(){ if [[ "$2" == "$3" ]]; then echo "  [ ok ] $1"; else echo "  [FAIL] $1: got '$2', want '$3'"; fail=1; fi; }

# accounting: 10 members in, 6 archived, 3 consolidated records emitted
chk "input records"      "$(grep -vc '^#' hervk_graph_fixture.vcf)" "10"
chk "output records"     "$(grep -vc '^#' "$TMP/out.vcf")"          "7"
chk "archived records"   "$(grep -vc '^#' "$TMP/archive.vcf")"      "6"
# 3 merged into new multi-allelic records, 2 left in place but given their
# locus identity (one member masked, one member unresolved)
chk "merged records"     "$(awk -F'\t' '!/^#/ && $3 ~ /^HERVK_/' "$TMP/out.vcf" | wc -l | tr -d ' ')" "3"
chk "loci annotated"     "$(grep -c 'HERVK_LOCUS=HERVK_' "$TMP/out.vcf" | tr -d ' ')" "5"

get(){ grep -m1 "$1" "$TMP/out.vcf" | tr '\t' '\n' | sed -n 8p | tr ';' '\n' | grep -m1 "^$2=" | cut -d= -f2; }

# chr11: both members are the same allele; graph and assemblies agree exactly
chk "chr11 AC"           "$(get HERVK_chr11_101704640 HERVK_AC)"      "23"
chk "chr11 AN"           "$(get HERVK_chr11_101704640 HERVK_AN)"      "40"
chk "chr11 discovery AC" "$(get HERVK_chr11_101704640 HERVK_AC_DISC)" "23"
grep -q 'HERVK_LOCUS=HERVK_chr11_101704640.*HERVK_DISC_CONCORDANT' "$TMP/out.vcf" \
  && echo "  [ ok ] chr11 flagged concordant with the assemblies" \
  || { echo "  [FAIL] chr11 should be flagged HERVK_DISC_CONCORDANT"; fail=1; }

# chr12: genuinely triallelic, and the two ploidy violations are the third
# allele flattened by bcftools norm -m-
chk "chr12 alleles"      "$(get HERVK_chr12_55299985 HERVK_ALLELE)"   "null,provirus"
chk "chr12 AC"           "$(get HERVK_chr12_55299985 HERVK_AC)"       "6,21"
chk "chr12 AN"           "$(get HERVK_chr12_55299985 HERVK_AN)"       "33"
chk "chr12 ploidy exceeded" "$(get HERVK_chr12_55299985 HERVK_N_PLOIDY_EXCEEDED)" "2"

# chr6: the tandem member is withheld, its partner keeps the locus identity
# match on the ID column: the partner record names the tandem member in its
# INFO, so a whole-line grep would pick that up too
tandem(){ awk -F'\t' '$3=="chr6-78894876-INS-8465_106221"' "$TMP/out.vcf"; }
tandem | grep -q 'HERVK_GT_MASKED' \
  && echo "  [ ok ] chr6 tandem record genotypes withheld" \
  || { echo "  [FAIL] chr6 tandem record should carry HERVK_GT_MASKED"; fail=1; }
chk "chr6 tandem has no called GT" \
    "$(tandem | cut -f10- | tr '\t' '\n' | cut -d: -f1 | grep -cv '^\.[/|]\?\.\?$' | tr -d ' ')" "0"
chk "chr6 partner keeps its genotypes" \
    "$(awk -F'\t' '$3=="chr6-78894317-DEL-8465_106220"' "$TMP/out.vcf" | cut -f10- | tr '\t' '\n' | cut -d: -f1 | grep -c '1' | tr -d ' ')" "8"

[[ $fail -eq 0 ]] && echo "PASS" || { echo "FAIL"; exit 1; }
