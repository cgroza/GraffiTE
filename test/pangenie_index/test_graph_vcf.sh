#!/usr/bin/env bash
# Regression test for the pangenie_index graph input (#93).
#
# merge_vcfs.py raised "An ID needs to be provided for each individual ID"
# when `bcftools norm -m+` had joined records at one position into a record
# with a different number of IDs and ALT alleles. That happens with two
# records of the same CHROM/POS/REF/ALT and different IDs, two records sharing
# an ID, and a record with ID "." beside another at its position. The test
# runs the process script, minus PanGenie-index, on those records and checks
# the tracking table.
#
# The fixture also has a record with two ALT alleles, an ID containing ":",
# three records at one position and a deletion that overlaps an insertion.
# merge_vcfs.py drops alleles from the last two, and the table must mark
# them in_graph=no.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export PATH="$(cd ../../bin && pwd):$PATH"

command -v bcftools >/dev/null || { echo "  [skip] bcftools not on PATH"; exit 0; }
python3 -c 'import pyfaidx' 2>/dev/null || { echo "  [skip] pyfaidx not importable"; exit 0; }

fail=0
chk(){ if [[ "$2" == "$3" ]]; then echo "  [ ok ] $1"; else
       echo "  [FAIL] $1: got '$2', want '$3'"; fail=1; fi; }

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

python3 - "$tmp" <<'EOF'
import random, sys
tmp = sys.argv[1]
random.seed(93)
seq = ''.join(random.choice('ACGT') for _ in range(3000))
open(f'{tmp}/ref.fa', 'w').write('>chr1\n' + '\n'.join(seq[i:i+60] for i in range(0, len(seq), 60)) + '\n')
def s(n): return ''.join(random.choice('ACGT') for _ in range(n))
I1, I2, I3 = s(60), s(70), s(80)
def ins(pos, vid, *alts):
    r = seq[pos-1]
    return f"chr1\t{pos}\t{vid}\t{r}\t{','.join(r + a for a in alts)}\t.\tPASS\tSVTYPE=INS\tGT\t1/1"
def dl(pos, vid, n):
    return f"chr1\t{pos}\t{vid}\t{seq[pos-1:pos+n]}\t{seq[pos-1]}\t.\tPASS\tSVTYPE=DEL\tGT\t0/1"
recs = [
    ins(100, 'dupA', I1), ins(100, 'dupB', I1),         # 1-2 same allele, two IDs
    ins(400, 'same', I1), ins(400, 'same', I2),         # 3-4 one ID, two alleles
    ins(700, '.', I1), ins(700, 'dotmate', I2),         # 5-6 missing ID
    ins(1000, 'multi', I1, I2),                         # 7   two ALT alleles
    ins(1300, 'chr1:1300:INS', I1),                     # 8   ID with ':'
    ins(1600, 't1', I1), ins(1600, 't2', I2), ins(1600, 't3', I3),  # 9-11
    dl(1900, 'ovDel', 200), ins(2000, 'ovIns', I1),     # 12-13 overlap
    ins(2500, 'plain', I1),                             # 14
]
hdr = ['##fileformat=VCFv4.2', '##contig=<ID=chr1,length=3000>',
       '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">',
       '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
       '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1']
open(f'{tmp}/pangenome.vcf', 'w').write('\n'.join(hdr + recs) + '\n')
EOF

cd "$tmp"
# A copy of the pangenie_index script in module/main.nf, minus PanGenie-index.
# Keep the two in step.
bcftools view -G -Ov -o sites.vcf pangenome.vcf
pangenie_graph_vcf.py prepare sites.vcf graph_input.vcf pangenie_graph_variants.tsv 2> prepare.log
bcftools sort graph_input.vcf 2>/dev/null | bcftools norm -m+ -Ov -o graph.vcf 2>/dev/null
# merge_vcfs.py names /usr/bin/python3, which may not be the python3 that has pyfaidx.
rc=0; python3 "$(command -v merge_vcfs.py)" merge -r ref.fa -v graph.vcf -ploidy 2 > graph_merged.vcf 2> merge.log || rc=$?
chk "merge_vcfs.py exits 0" "$rc" "0"
pangenie_graph_vcf.py report pangenie_graph_variants.tsv graph_merged.vcf 2> report.log

row(){ awk -F'\t' -v r="$1" -v a="${2:-1}" '$1==r && $5==a {print $6, $7, $8}' pangenie_graph_variants.tsv; }
chk "one row per ALT allele"          "$(tail -n +2 pangenie_graph_variants.tsv | wc -l | tr -d ' ')" "15"
chk "first of a duplicate pair"       "$(row 1)" "dupA yes ."
chk "second of a duplicate pair"      "$(row 2)" "dupA yes duplicate_of_record_1"
chk "shared ID, first"                "$(row 3)" "same yes ."
chk "shared ID, second renamed"       "$(row 4)" "graffite_rec4 yes id_replaced"
chk "missing ID renamed"              "$(row 5)" "graffite_rec5 yes id_replaced"
chk "multi-ALT, allele 2"             "$(row 7 2)" "graffite_rec7_2 yes id_replaced"
chk "ID with ':' renamed"             "$(row 8)" "graffite_rec8 yes id_replaced"
# Which of the three is left out depends on the ALT order bcftools norm writes.
chk "three alleles at one position"   "$(for r in 9 10 11; do row $r; done | awk '$2=="no"' | wc -l | tr -d ' ')" "1"
chk "overlapping deletion"            "$(row 12)" "ovDel no ."
chk "overlapping insertion"           "$(row 13)" "ovIns no ."
chk "plain record keeps its ID"       "$(row 14)" "plain yes ."
chk "report counts"                   "$(cat report.log)" \
    "pangenie_graph_vcf.py: 3 of 15 pangenome.vcf alleles are not in the PanGenie graph and get no genotype"

# The graph IDs marked in_graph=yes are the same set as the IDs in INFO/ID of
# the graph, which PanGenie copies to the genotyped VCFs.
ids=$(grep -v '^#' graph_merged.vcf | cut -f8 | sed 's/.*ID=//' | tr ',:' '\n\n' | sort -u)
want=$(awk -F'\t' 'NR>1 && $7=="yes" {print $6}' pangenie_graph_variants.tsv | sort -u)
chk "table matches graph INFO/ID"     "$ids" "$want"

exit $fail
