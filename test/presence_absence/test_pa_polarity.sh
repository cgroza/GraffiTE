#!/usr/bin/env bash
# vcf_to_pa_tsv.py has to read insertion polarity without INFO/SVTYPE.
#
# truvari_merge strips every upstream INFO field before the collapse and puts
# only SVLEN back, so on any run with two or more caller VCFs pangenome.vcf
# carries no SVTYPE. gt_presence() keyed on SVTYPE alone and returned NA for
# every sample, which made the presence-absence TSVs carry no information at
# all on the pipeline's main path. Measured before the fix on a 29-record run:
# 116 of 116 sample values NA.
#
# The first block feeds records that DO carry SVTYPE, so a fix that ignored it
# and guessed from lengths would still have to agree with the field.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
BIN="$(cd ../../bin && pwd)"

fail=0
chk(){ if [[ "$2" == "$3" ]]; then echo "  [ ok ] $1"; else
       echo "  [FAIL] $1: got '$2', want '$3'"; fail=1; fi; }

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

hdr() {
  printf '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n'
  printf '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="t">\n'
  printf '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="l">\n'
  printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tpresent\tabsent\n'
}

# With SVTYPE: an insertion and a deletion, each carried by one sample.
{ hdr
  printf 'chr1\t100\twith_ins\tA\tACCCCCCCCC\t.\tPASS\tSVTYPE=INS;SVLEN=9\tGT\t1/1\t0/0\n'
  printf 'chr1\t200\twith_del\tACCCCCCCCC\tA\t.\tPASS\tSVTYPE=DEL;SVLEN=-9\tGT\t1/1\t0/0\n'
} > "$tmp/with.vcf"

# The same records with SVTYPE removed, which is what the merge leaves.
{ hdr
  printf 'chr1\t100\tno_ins\tA\tACCCCCCCCC\t.\tPASS\tSVLEN=9\tGT\t1/1\t0/0\n'
  printf 'chr1\t200\tno_del\tACCCCCCCCC\tA\t.\tPASS\tSVLEN=-9\tGT\t1/1\t0/0\n'
  printf 'chr1\t300\tsymbolic\tA\t<INV>\t.\tPASS\tSVLEN=100\tGT\t1/1\t0/0\n'
  printf 'chr1\t400\tsame_len\tACGT\tACGA\t.\tPASS\t.\tGT\t1/1\t0/0\n'
} > "$tmp/without.vcf"

pa(){ python3 "$BIN/vcf_to_pa_tsv.py" "$1" -o "$2" >/dev/null 2>&1; }
col(){ awk -F'\t' -v id="$2" 'NR>1 && $4==id {print $(NF-1)"/"$NF}' "$1"; }

pa "$tmp/with.vcf" "$tmp/with.tsv"
chk "SVTYPE=INS: carrier present, other absent"  "$(col "$tmp/with.tsv" with_ins)" "1/0"
chk "SVTYPE=DEL: carrier absent, other present"  "$(col "$tmp/with.tsv" with_del)" "0/1"

pa "$tmp/without.vcf" "$tmp/without.tsv"
chk "no SVTYPE, longer ALT: reads as an insertion" "$(col "$tmp/without.tsv" no_ins)" "1/0"
chk "no SVTYPE, longer REF: reads as a deletion"   "$(col "$tmp/without.tsv" no_del)" "0/1"
chk "no SVTYPE, symbolic ALT: still NA"            "$(col "$tmp/without.tsv" symbolic)" "NA/NA"
chk "no SVTYPE, equal lengths: still NA"           "$(col "$tmp/without.tsv" same_len)" "NA/NA"

chk "no row is NA just because SVTYPE is gone" \
    "$(awk -F'\t' 'NR>1 && $4 ~ /^no_/ {print $(NF-1)$NF}' "$tmp/without.tsv" | grep -c NA)" "0"

exit $fail
