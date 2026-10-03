#!/usr/bin/env bash
# Regression test for add_polyA.py when the reported TSD is part of the tail.
#
# add_polyA.py trims the TSD copy from the end it scans before looking for the
# tail, which exposes a tail that sits behind the copy. But the TSD search
# sometimes reports a run of the tail base as the duplication, and trimming it
# removed the tail. The two real Alu insertions below, both from chm13v2.0, have
# a tail at the expected end and came back polyA=FALSE: chr2:212,323,901 ends in
# 9 A and its TSD is AAAAA, which leaves 4; chr22:25,461,790 starts with 20 T
# and its TSD is 19 T, which leaves 1. The tail now counts with or without the
# trim. The two synthetic cases check that a tail found only once its TSD copy
# is trimmed is still found, and that a sequence without a tail still is not.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export PATH="$(cd ../../bin && pwd):$PATH"

fail=0
chk(){ if [[ "$2" == "$3" ]]; then echo "  [ ok ] $1"; else
       echo "  [FAIL] $1: got '$2', want '$3'"; fail=1; fi; }

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

ALU_PLUS=AAAAATGTTTTTCTAGACATCCTTGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCCGGCTAAAACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGTAGTGGCGGGCGCCTGTAGTCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCCCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAA
ALU_MINUS=TTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGCTGGAGTGCAGTGGCGGGATCTCGGCTCACTGCAAGCTCCGCCTCCCGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAGGCGCGCGCCACTACGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTTTTAGCCGGGATGGTCTCGATCTCCTGACCTCGTGATCCGCCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGCC
T19=TTTTTTTTTTTTTTTTTTT
BEHIND=GTCAGCGTCAGTGCATGCAGTCAAAAAAAAAAGCTCGTCCGC  # 10 A, then the 3' TSD copy
NOTAIL=GTCAGCGTCAGTGCATGCAGTCGCTGCAGGTCGATC

{
  printf '##fileformat=VCFv4.2\n'
  printf '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type">\n'
  printf '##INFO=<ID=n_hits,Number=1,Type=Integer,Description="Number of repeats">\n'
  printf '##INFO=<ID=RM_hit_strands,Number=.,Type=String,Description="RepeatMasker hit strands">\n'
  printf '##INFO=<ID=TSD,Number=2,Type=String,Description="Target site duplication">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  printf 'chr2\t212323900\talu_plus\tT\tT%s\t.\tPASS\tSVTYPE=INS;n_hits=1;RM_hit_strands=+;TSD=AAAAA,AAAAA\n' "$ALU_PLUS"
  printf 'chr22\t25461789\talu_minus\tT\tT%s\t.\tPASS\tSVTYPE=INS;n_hits=1;RM_hit_strands=C;TSD=%s,%s\n' "$ALU_MINUS" "$T19" "$T19"
  printf 'chrS\t100\tbehind\tC\tC%s\t.\tPASS\tSVTYPE=INS;n_hits=1;RM_hit_strands=+;TSD=GCTCGTCCGC,GCTCGTCCGC\n' "$BEHIND"
  printf 'chrS\t200\tnotail\tC\tC%s\t.\tPASS\tSVTYPE=INS;n_hits=1;RM_hit_strands=+;TSD=GATC,GATC\n' "$NOTAIL"
} > "$tmp/in.vcf"

add_polyA.py "$tmp/in.vcf" -o "$tmp/out.vcf"
call(){ awk -F'\t' -v id="$1" '$3 == id { n = split($8, kv, ";"); for (i = 1; i <= n; i++) if (kv[i] ~ /^polyA=/) print substr(kv[i], 7) }' "$tmp/out.vcf"; }

echo "a TSD that is a run of the tail base"
chk "chr2:212,323,901, + strand, TSD AAAAA"  "$(call alu_plus)"  "TRUE"
chk "chr22:25,461,790, C strand, TSD 19 T"   "$(call alu_minus)" "TRUE"
echo "controls"
chk "a tail behind its TSD copy is found"    "$(call behind)"    "TRUE"
chk "no tail, no call"                       "$(call notail)"    "FALSE"

exit ${fail}
