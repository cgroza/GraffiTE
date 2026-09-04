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

has(){ awk -F'\t' -v id="$2" '!/^#/ && $3==id' "$TMP/out.vcf" \
       | grep -c "$1" | tr -d ' '; }

# The locus flags travel from hervk_loci.tsv onto the record, so the two say
# the same thing and neither recomputes it. chr6 is the case the single
# locus_type label cannot express: solo, provirus and a two-unit allele.
chk "chr6 is labelled copy_number" "$(get HERVK_chr6_78894316 HERVK_LOCUS_TYPE)" "copy_number"
chk "chr6 is flagged solo/provirus too" "$(has 'HERVK_SOLO_PROV' HERVK_chr6_78894316)" "1"
chk "chr6 is flagged copy-number"       "$(has 'HERVK_CNV' HERVK_chr6_78894316)"       "1"
chk "chr6 is not an insertion polymorphism" "$(has 'HERVK_MEI' HERVK_chr6_78894316)"   "0"
chk "chr12 is an insertion polymorphism"    "$(has 'HERVK_MEI' HERVK_chr12_55299985)"  "1"
chk "chr12 is flagged solo/provirus too"    "$(has 'HERVK_SOLO_PROV' HERVK_chr12_55299985)" "1"
chk "chr11 is not an insertion polymorphism" "$(has 'HERVK_MEI' HERVK_chr11_101704640)" "0"

# A locus the --human filter cut in half. chr8 keeps one of two members, so its
# counts cannot cover the locus and the record has to say so rather than read
# as a complete allele frequency.
chk "chr8 says a member is missing" \
    "$(grep -c 'HERVK_LOCUS=HERVK_chr8_7552031.*HERVK_LOCUS_INCOMPLETE' "$TMP/out.vcf")" "1"
chk "chr8 names the missing member" \
    "$(grep -m1 'HERVK_LOCUS=HERVK_chr8_7552031' "$TMP/out.vcf" | grep -o 'HERVK_MEMBERS_ABSENT=[^;]*' | cut -d= -f2)" \
    "chr8-7552072-INS-218_114261"
chk "chr8 still reports the whole allele set" \
    "$(grep -m1 'HERVK_LOCUS=HERVK_chr8_7552031' "$TMP/out.vcf" | grep -o 'HERVK_ALLELE_SET=[^;]*' | cut -d= -f2)" \
    "null,solo"
chk "a complete locus is not flagged incomplete" \
    "$(has 'HERVK_LOCUS_INCOMPLETE' HERVK_chr6_78894316)" "0"

# --- the same machinery over the discovery callset --------------------------
# Assembly alignments resolve a tandem array directly, so nothing is masked
# there: chr6's two-unit allele is counted rather than withheld. The genotyper
# guard does not apply, and neither does the discovery fallback, whose source
# would be this file.
python3 $BIN/hervk_reconcile.py consolidate --source discovery \
    --vcf-in hervk_discovery_fixture.vcf \
    --loci   hervk_loci_fixture.tsv \
    --calls  hervk_calls_fixture.tsv \
    --out-vcf "$TMP/disc.vcf" 2>"$TMP/derr"

d(){ grep -m1 "$1" "$TMP/disc.vcf" | tr '\t' '\n' | sed -n 8p | tr ';' '\n' \
     | grep -m1 "^$2=" | cut -d= -f2; }

chk "discovery consolidates chr6" \
    "$(awk -F'\t' '!/^#/ && $3=="HERVK_chr6_78894316"' "$TMP/disc.vcf" | wc -l | tr -d ' ')" "1"
chk "discovery counts the two-unit allele" "$(d HERVK_chr6_78894316 HERVK_AC)" "8,1"
chk "discovery masks nothing at chr6" \
    "$(grep -m1 'HERVK_LOCUS=HERVK_chr6_78894316' "$TMP/disc.vcf" | grep -c 'HERVK_ALLELE_NOGT')" "0"
chk "discovery keeps both flags at chr6" \
    "$(grep -m1 'HERVK_LOCUS=HERVK_chr6_78894316' "$TMP/disc.vcf" | grep -c 'HERVK_SOLO_PROV.*HERVK_CNV')" "1"
chk "discovery records its source" \
    "$(grep -c '^##hervk_consolidation=source:discovery' "$TMP/disc.vcf")" "1"
# No graph to compare against, so no ploidy claim either way. The flag used to
# fire on every record of a run without --discovery-vcf, reading "nothing to
# check" as "the two disagree".
chk "discovery claims no ploidy mismatch" \
    "$(awk -F'\t' '!/^#/' "$TMP/disc.vcf" | grep -c 'HERVK_DISC_PLOIDY_MISMATCH')" "0"

# --- a haploid discovery caller ---------------------------------------------
# SVIM-asm run per haplotype emits GT 1 or 0, not 0/1. Falling back to those
# counts on a diploid denominator would halve every frequency, so the fallback
# is refused and says why. AN comes from the GT, never from an assumed 2.
sed 's/\([0-9]\)|\([0-9]\)/\1/g; s/\([0-9]\)\/\([0-9]\)/\1/g' \
    hervk_discovery_fixture.vcf > "$TMP/haploid.vcf"
python3 $BIN/hervk_reconcile.py consolidate \
    --genotyped-vcf hervk_graph_fixture.vcf \
    --loci   hervk_loci_fixture.tsv \
    --calls  hervk_calls_fixture.tsv \
    --discovery-vcf "$TMP/haploid.vcf" \
    --out-vcf "$TMP/hap.vcf" 2>"$TMP/haperr"

chk "haploid discovery is refused as a fallback" \
    "$(grep -c 'HERVK_LOCUS=HERVK_chr11_101704640.*HERVK_DISC_PLOIDY_MISMATCH' "$TMP/hap.vcf")" "1"
chk "no discovery counts are reported from it" \
    "$(awk -F'\t' '!/^#/' "$TMP/hap.vcf" | grep -c 'HERVK_AC_DISC')" "0"
chk "graph counts are untouched by it" \
    "$(grep -m1 'HERVK_LOCUS=HERVK_chr11_101704640' "$TMP/hap.vcf" | grep -o 'HERVK_AN=[0-9]*' | cut -d= -f2)" "40"
# Consolidating that haploid callset on its own is fine: AN follows its ploidy.
python3 $BIN/hervk_reconcile.py consolidate --source discovery \
    --vcf-in "$TMP/haploid.vcf" --loci hervk_loci_fixture.tsv \
    --calls hervk_calls_fixture.tsv --out-vcf "$TMP/haploid_cons.vcf" 2>/dev/null
chk "a haploid callset consolidates at AN=20" \
    "$(grep -m1 'HERVK_LOCUS=HERVK_chr11_101704640' "$TMP/haploid_cons.vcf" | grep -o 'HERVK_AN=[0-9]*' | cut -d= -f2)" "20"
# --discovery-vcf would name this same file, which cross-checks nothing.
if python3 $BIN/hervk_reconcile.py consolidate --source discovery \
     --vcf-in hervk_discovery_fixture.vcf --loci hervk_loci_fixture.tsv \
     --calls hervk_calls_fixture.tsv --discovery-vcf hervk_discovery_fixture.vcf \
     --out-vcf "$TMP/x.vcf" 2>/dev/null; then
  echo "  [FAIL] --source discovery should refuse --discovery-vcf"; fail=1
else
  echo "  [ ok ] --source discovery refuses --discovery-vcf"
fi

[[ $fail -eq 0 ]] && echo "PASS" || { echo "FAIL"; exit 1; }
