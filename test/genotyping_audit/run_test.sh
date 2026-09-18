#!/usr/bin/env bash
# Regression test for bin/genotyping_audit.py.
#
# Four ways a pangenome.vcf allele fails to reach the merged genotypes, one
# fixture record each, plus one that arrives intact. The PanGenie fixture
# carries INFO/ID and a graph table so the audit can separate "never entered the
# graph" from "genotyped but unmatched by merge_VCFs"; the vg fixture has
# neither, and the audit has to say so rather than guess.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="$(cd "${HERE}/../../bin" && pwd):$PATH"

fail=0
chk(){ if [[ "$2" == "$3" ]]; then echo "  [ ok ] $1"; else
       echo "  [FAIL] $1: got '$2', want '$3'"; fail=1; fi; }

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

{ echo '##fileformat=VCFv4.2'
  echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">'
  echo '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="len">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  printf 'chr1\t100\tkept\tA\tAGGG\t.\tPASS\tSVTYPE=INS;SVLEN=3\n'
  printf 'chr1\t200\tdupA\tA\tAGGG\t.\tPASS\tSVTYPE=INS;SVLEN=3\n'
  printf 'chr1\t200\tdupB\tA\tAGGG\t.\tPASS\tSVTYPE=INS;SVLEN=3\n'
  printf 'chr1\t300\thasN\tA\tANNN\t.\tPASS\tSVTYPE=INS;SVLEN=3\n'
  printf 'chr1\t400\tshifted\tA\tAGGG\t.\tPASS\tSVTYPE=INS;SVLEN=3\n'
} > "$tmp/pangenome.vcf"

# pangenie_graph_vcf.py output: hasN never reached the graph, dupB shares dupA's
# allele, the rest each got their own graph variant.
{ printf 'record\tCHROM\tPOS\tpangenome_ID\tallele\tgraph_ID\tin_graph\tnote\n'
  printf '1\tchr1\t100\tkept\t1\tkept\tyes\t.\n'
  printf '2\tchr1\t200\tdupA\t1\tdupA\tyes\t.\n'
  printf '3\tchr1\t200\tdupB\t1\tdupA\tyes\tduplicate_of_record_2\n'
  printf '4\tchr1\t300\thasN\t1\thasN\tno\t.\n'
  printf '5\tchr1\t400\tshifted\t1\tshifted\tyes\t.\n'
} > "$tmp/graph_table.tsv"

# PanGenie's merged output. "shifted" was genotyped -- INFO/ID names it -- but
# left-alignment moved it, so merge_VCFs' annotate found no match and it has no
# ID and no annotation.
{ echo '##fileformat=VCFv4.2'
  echo '##source=PanGenie'
  echo '##INFO=<ID=ID,Number=A,Type=String,Description="graph variant">'
  echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">'
  echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n'
  printf 'chr1\t100\tkept\tA\tAGGG\t.\t.\tID=kept;SVTYPE=INS\tGT\t0/1\t1/1\n'
  printf 'chr1\t200\tdupA\tA\tAGGG\t.\t.\tID=dupA;SVTYPE=INS\tGT\t0/1\t0/0\n'
  printf 'chr1\t397\t.\tC\tCGGG\t.\t.\tID=shifted\tGT\t0/1\t./.\n'
} > "$tmp/merged_pangenie.vcf"

printf 'kept\ndupA\n' > "$tmp/trusted.ids"

genotyping_audit.py --pangenome "$tmp/pangenome.vcf" --merged "$tmp/merged_pangenie.vcf" \
  --graph-table "$tmp/graph_table.tsv" --trusted-ids "$tmp/trusted.ids" \
  -o "$tmp/audit.tsv" 2>/dev/null

at(){ awk -F'\t' -v id="$1" -v col="$2" '
  /^#/ {next} NR_h==0 && /^pangenome_ID/ {for(i=1;i<=NF;i++) h[$i]=i; NR_h=1; next}
  $1==id {print $(h[col])}' "$tmp/audit.tsv"; }

chk "backend named in the summary" \
    "$(grep -c '^# genotyper: pangenie (join key: INFO/ID)' "$tmp/audit.tsv")" "1"
chk "one row per ALT allele" \
    "$(grep -vc '^#' "$tmp/audit.tsv")" "6"   # 5 records + the column header

chk "a matched record is genotyped"        "$(at kept lost_at)"    "genotyped"
chk "its samples are counted"              "$(at kept n_samples_genotyped)" "2"
chk "the trusted flag is carried"          "$(at kept trusted)"    "yes"
chk "an untrusted record says so"          "$(at shifted trusted)" "no"

chk "the N-containing allele never entered the graph" "$(at hasN lost_at)" "not_in_graph"
chk "and is not reported as annotated"                "$(at hasN annotated)" "no"

chk "the first of a duplicate pair is genotyped" "$(at dupA lost_at)" "genotyped"
chk "the second names the record it merged into" \
    "$(at dupB lost_at)" "duplicate_of_record_2"

chk "a genotyped record the merge could not match is flagged" \
    "$(at shifted lost_at)" "no_match_in_merge"
chk "its calls are still counted through INFO/ID" \
    "$(at shifted n_samples_genotyped)" "1"
chk "it is recorded as genotyped"       "$(at shifted genotyped)" "yes"
chk "but not as annotated"              "$(at shifted annotated)" "no"
chk "the summary counts it"             \
    "$(grep -c '^# genotyped but unmatched by merge_VCFs: 1' "$tmp/audit.tsv")" "1"

# vg: no INFO/ID, no graph table. The audit falls back on the ID column and says
# it cannot separate the two loss modes.
{ echo '##fileformat=VCFv4.2'
  echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  echo '##FORMAT=<ID=MAD,Number=1,Type=Integer,Description="Minimum site allele depth">'
  echo '##FORMAT=<ID=XD,Number=1,Type=Float,Description="eXpected Depth, background coverage as used for the Poisson model">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n'
  printf 'chr1\t100\tkept\tA\tAGGG\t.\t.\t.\tGT:MAD:XD\t0/1:9:30\n'
  printf 'chr1\t900\t>12345>678\tA\tAT\t.\t.\t.\tGT:MAD:XD\t0/1:4:30\n'
} > "$tmp/merged_vg.vcf"

genotyping_audit.py --pangenome "$tmp/pangenome.vcf" --merged "$tmp/merged_vg.vcf" \
  -o "$tmp/audit_vg.tsv" 2>/dev/null
chk "vg is detected and the ID column is used" \
    "$(grep -c '^# genotyper: giraffe (join key: ID)' "$tmp/audit_vg.tsv")" "1"
chk "without a graph table the report says what it cannot tell apart" \
    "$(grep -c 'not separable' "$tmp/audit_vg.tsv")" "1"
chk "a snarl record in the merged VCF adds no row" \
    "$(grep -vc '^#' "$tmp/audit_vg.tsv")" "6"
chk "the ID-column caveat is stated" \
    "$(grep -c 'only the PanGenie path carries INFO/ID' "$tmp/audit_vg.tsv")" "1"

if [[ $fail -eq 0 ]]; then echo PASS; else echo FAIL; exit 1; fi
