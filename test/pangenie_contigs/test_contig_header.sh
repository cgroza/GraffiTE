#!/usr/bin/env bash
# Regression test for the contig headers on PanGenie's output.
#
# PanGenie writes <sample>_genotyping.vcf with no ##contig lines. The pangenie
# process then ran `bcftools norm -Oz` straight off it, and bcftools needs a
# header contig to BCF-encode a record, so it stopped on the first one with
# "Invalid BCF, CONTIG id=0 not present in the header" and left a zero-byte
# .vcf.gz. Every --graph_method pangenie run died there.
#
# The fixture is a PanGenie-shaped VCF over two contigs. The first check runs
# the old chain and requires it to fail: without that, a fixture that had
# contig lines in it would let the second check pass while testing nothing.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

command -v bcftools >/dev/null || { echo "  [skip] bcftools not on PATH"; exit 0; }
command -v samtools >/dev/null || { echo "  [skip] samtools not on PATH"; exit 0; }
command -v tabix    >/dev/null || { echo "  [skip] tabix not on PATH";    exit 0; }

fail=0
chk(){ if [[ "$2" == "$3" ]]; then echo "  [ ok ] $1"; else
       echo "  [FAIL] $1: got '$2', want '$3'"; fail=1; fi; }

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

python3 - "$tmp" <<'EOF'
import random, sys
tmp = sys.argv[1]
random.seed(93)
contigs = {'chr1': 4000, 'X': 3000}
seqs = {c: ''.join(random.choice('ACGT') for _ in range(n)) for c, n in contigs.items()}
with open(f'{tmp}/ref.fa', 'w') as fh:
    for c, s in seqs.items():
        fh.write(f'>{c}\n' + '\n'.join(s[i:i+60] for i in range(0, len(s), 60)) + '\n')

# PanGenie's header, as it writes it: no ##contig anywhere.
hdr = ['##fileformat=VCFv4.2',
       '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
       '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
       '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
       '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2']
def ins(c, pos, n, gt):
    r = seqs[c][pos-1]
    alt = r + ''.join(random.choice('ACGT') for _ in range(n))
    return f'{c}\t{pos}\t.\t{r}\t{alt}\t.\tPASS\tAF=0.5\tGT:GQ\t{gt}:60'
recs = [ins('chr1', 500, 300, '1/1'), ins('chr1', 1500, 150, '0/1'),
        ins('X', 800, 250, '1/1'),    ins('X', 2000, 90, '0/0')]
open(f'{tmp}/S2_genotyping.vcf', 'w').write('\n'.join(hdr + recs) + '\n')
EOF

chk "fixture carries no contig lines" \
    "$(grep -c '^##contig' "$tmp/S2_genotyping.vcf")" "0"

# 1. the chain as it was: must fail, or this fixture proves nothing
bcftools norm -f "$tmp/ref.fa" -N -m- -Oz \
         -o "$tmp/old.vcf.gz" "$tmp/S2_genotyping.vcf" >/dev/null 2>&1
chk "norm alone still fails on a contig-less header" "$?" "255"

# 2. the chain as it stands in module/main.nf
samtools faidx "$tmp/ref.fa" >/dev/null 2>&1
bcftools reheader -f "$tmp/ref.fa.fai" "$tmp/S2_genotyping.vcf" \
         > "$tmp/S2_genotyping.contigs.vcf" 2>/dev/null
chk "reheader carries the reference contigs over" \
    "$(grep -c '^##contig' "$tmp/S2_genotyping.contigs.vcf")" "2"

bcftools norm -f "$tmp/ref.fa" -N -m- -Oz \
         -o "$tmp/new.vcf.gz" "$tmp/S2_genotyping.contigs.vcf" >/dev/null 2>&1
chk "norm writes the genotypes out" "$?" "0"
tabix -p vcf "$tmp/new.vcf.gz" >/dev/null 2>&1
chk "tabix indexes the result" "$?" "0"

chk "every record survives" "$(bcftools view -H "$tmp/new.vcf.gz" | wc -l)" "4"
# LC_ALL=C or the answer depends on the collation the host happens to set:
# a UTF-8 locale sorts this chr1,X and the C locale sorts it X,chr1.
chk "both contigs still carry records" \
    "$(bcftools view -H "$tmp/new.vcf.gz" | cut -f1 | LC_ALL=C sort -u | paste -sd, -)" "X,chr1"
chk "the genotypes are the ones PanGenie called" \
    "$(bcftools query -f '[%GT]\n' "$tmp/new.vcf.gz" | paste -sd, -)" "1/1,0/1,1/1,0/0"

exit $fail
