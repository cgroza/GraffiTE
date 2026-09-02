#!/usr/bin/env bash
# Regression test for `hervk_reconcile.py consolidate`.
#
# Real CaG data: the ten member records of the five flagged loci, from the
# unfiltered graph-genotyped VCF, with the assembly callset alongside. These
# five cover every shape the consolidator has to handle -- two loci where both
# members describe the same allele, one genuinely triallelic locus, one where a
# member carries a copy-number allele the graph cannot count, and one where a
# member has no usable allele state.
#
# The copy-number case is the one that changed shape: masking used to drop the
# member from the merge, leaving one usable member and no consolidation, so the
# locus stayed as two separate records. It consolidates now, and the masking is
# per allele -- the clean allele keeps its graph genotypes and the copy-number
# one is named in HERVK_ALLELE_NOGT with its count in HERVK_AC_DISC.
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
chk "output records"     "$(grep -vc '^#' "$TMP/out.vcf")"          "6"
chk "archived records"   "$(grep -vc '^#' "$TMP/archive.vcf")"      "8"
# 4 merged into new multi-allelic records, 1 left in place and given its locus
# identity (its partner has no usable allele state)
chk "merged records"     "$(awk -F'\t' '!/^#/ && $3 ~ /^HERVK_/' "$TMP/out.vcf" | wc -l | tr -d ' ')" "4"
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

# chr6: both alleles are on one record now. The solo allele is genotyped from
# the graph; the copy-number allele is not, and says so rather than reporting a
# 0 that would read as absent.
chk "chr6 consolidates"  "$(awk -F'\t' '!/^#/ && $3=="HERVK_chr6_78894316"' "$TMP/out.vcf" | wc -l | tr -d ' ')" "1"
chk "chr6 carries both alleles" "$(get HERVK_chr6_78894316 HERVK_ALLELE)" "solo,prov_x2"
chk "chr6 graph counts the solo allele only" "$(get HERVK_chr6_78894316 HERVK_AC)" "8,0"
chk "chr6 names the allele it cannot genotype" \
    "$(get HERVK_chr6_78894316 HERVK_ALLELE_NOGT)" "prov_x2"
chk "chr6 discovery counts both alleles" \
    "$(get HERVK_chr6_78894316 HERVK_AC_DISC)" "8,1"
chk "chr6 keeps the clean allele's carriers" \
    "$(awk -F'\t' '!/^#/ && $3=="HERVK_chr6_78894316"' "$TMP/out.vcf" | cut -f10- | tr '\t' '\n' | cut -d: -f1 | grep -c '1' | tr -d ' ')" "8"
chk "chr6 is not flagged concordant" \
    "$(awk -F'\t' '!/^#/ && $3=="HERVK_chr6_78894316"' "$TMP/out.vcf" | grep -c 'HERVK_DISC_CONCORDANT' | tr -d ' ')" "0"

# Header validity. The input fixture defines SVTYPE (Number=1) and SVLEN
# (Number=.), and the consolidated records carry one value per ALT, so our
# Number=A definitions have to REPLACE those rather than sit beside them --
# two ##INFO lines with one ID is invalid VCF and readers disagree on the winner.
chk "no INFO id defined twice" \
    "$(grep '^##INFO=<ID=' "$TMP/out.vcf" | sed 's/##INFO=<ID=\([^,>]*\).*/\1/' | sort | uniq -d | wc -l | tr -d ' ')" "0"
chk "SVTYPE defined once"  "$(grep -c '##INFO=<ID=SVTYPE,' "$TMP/out.vcf")" "1"
chk "SVLEN defined once"   "$(grep -c '##INFO=<ID=SVLEN,'  "$TMP/out.vcf")" "1"
chk "SVTYPE is per-ALT"    "$(grep -c '##INFO=<ID=SVTYPE,Number=A,' "$TMP/out.vcf")" "1"
chk "SVLEN is per-ALT"     "$(grep -c '##INFO=<ID=SVLEN,Number=A,'  "$TMP/out.vcf")" "1"
# the multiallelic locus is what makes Number=A load-bearing
chk "multiallelic SVLEN kept per-ALT" \
    "$(awk -F'\t' '$3=="HERVK_chr12_55299985"' "$TMP/out.vcf" | grep -o 'SVLEN=[^;]*' | head -1)" "SVLEN=-974,8212"

[[ $fail -eq 0 ]] && echo "PASS" || { echo "FAIL"; exit 1; }
