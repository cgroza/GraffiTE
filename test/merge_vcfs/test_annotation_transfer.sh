#!/usr/bin/env bash
# Regression test for the INFO transfer in merge_VCFs.
#
# merge_VCFs copies the annotation off pangenome.vcf onto the genotyped records
# with `bcftools annotate -a pangenome.vcf -c CHROM,POS,ID,REF,ALT,INFO`. With a
# VCF as -a that matches on CHROM, POS, REF and a shared ALT; ID is transferred,
# not matched.
#
# The pangenie path ran `bcftools norm -f <ref> -m-` to split the multi-allelic
# records PanGenie writes, and -f turns on left-alignment as well. pangenome.vcf
# is left-aligned only when two or more caller VCFs went through the truvari
# merge, so on a --vcf, single-caller or --graffite_vcf run an insertion inside a
# homopolymer came back at a different POS and matched nothing: no repeat
# annotation, no TSD, no polyA, and an empty ID. -N keeps the position.
#
# Also pinned here: the match is case-insensitive over alleles. bin/fix_vcf.py
# leaves REF soft-masked lowercase and bin/merge_vcfs.py upper-cases it for the
# graph, so every TE record crosses that boundary.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

for t in bcftools bgzip tabix samtools; do
  command -v $t >/dev/null || { echo "  [skip] $t not on PATH"; exit 0; }
done

fail=0
chk(){ if [[ "$2" == "$3" ]]; then echo "  [ ok ] $1"; else
       echo "  [FAIL] $1: got '$2', want '$3'"; fail=1; fi; }

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
cd "$tmp"

# 100 bp of C, a 10 bp A homopolymer at 101-110, then G. An inserted A anywhere
# in the run is the same variant; only POS 100 is the left-aligned form.
python3 - <<'PY'
s = "C" * 100 + "A" * 10 + "G" * 90
open("ref.fa", "w").write(">ctg1\n" + "\n".join(s[i:i+60] for i in range(0, len(s), 60)) + "\n")
PY
samtools faidx ref.fa

# pangenome.vcf as GraffiTE writes it: not left-aligned, REF soft-masked.
{ echo '##fileformat=VCFv4.2'
  echo '##contig=<ID=ctg1,length=200>'
  echo '##INFO=<ID=n_hits,Number=1,Type=Integer,Description="Number of repeats found in insertion">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  printf 'ctg1\t110\tctg1-110-INS-1\ta\taa\t.\tPASS\tn_hits=1\n'
} > pangenome.vcf
bgzip -c pangenome.vcf > pangenome.vcf.gz && tabix -p vcf pangenome.vcf.gz

# what PanGenie hands back: multi-allelic, alleles upper-cased by the graph.
{ echo '##fileformat=VCFv4.2'
  echo '##contig=<ID=ctg1,length=200>'
  echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n'
  printf 'ctg1\t110\t.\tA\tAA,AAA\t.\t.\t.\tGT\t1/2\n'
} > called.vcf

transfer() {   # $1 = extra norm flags; echoes the INFO of the matched record
  bcftools norm -f ref.fa $1 -m- -Oz -o n.vcf.gz called.vcf 2>/dev/null
  tabix -f -p vcf n.vcf.gz
  bcftools annotate -a pangenome.vcf.gz -c CHROM,POS,ID,REF,ALT,INFO n.vcf.gz 2>/dev/null \
    | grep -v '^#' | head -1
}

old=$(transfer "")
chk "without -N the record is realigned off its position" "$(echo "$old" | cut -f2)" "100"
chk "without -N it loses its ID"                          "$(echo "$old" | cut -f3)" "."
chk "without -N it loses the annotation"                  "$(echo "$old" | cut -f8)" "."

new=$(transfer "-N")
chk "with -N the position is kept"       "$(echo "$new" | cut -f2)" "110"
chk "with -N the ID is transferred"      "$(echo "$new" | cut -f3)" "ctg1-110-INS-1"
chk "with -N the annotation arrives"     "$(echo "$new" | cut -f8)" "n_hits=1"
chk "with -N the record is still split"  \
    "$(bcftools view -H n.vcf.gz | wc -l | tr -d ' ')" "2"
chk "lowercase REF matches the upper-cased graph allele" \
    "$(echo "$new" | cut -f4)" "a"

if [[ $fail -eq 0 ]]; then echo PASS; else echo FAIL; exit 1; fi
