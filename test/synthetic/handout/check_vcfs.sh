#!/usr/bin/env bash
# The vcfs cell hands the spine's own per-sample vg call VCFs back through
# --vcfs. With no alignment and no vg_call, the run has nothing to recompute, so
# its merged genotypes must match the spine's record for record.
#
#   check_vcfs.sh <spine run dir> <vcfs run dir>
set -uo pipefail
S=${1:?spine run dir}; V=${2:?vcfs run dir}
fail=0
ok()  { echo "  [ ok ] $1"; }
bad() { echo "  [FAIL] $1"; fail=1; }

M=4_Genotyping/GraffiTE.merged.genotypes.vcf.gz
for d in "$S" "$V"; do
  [[ -s "$d/$M" ]] || { bad "no $M in $d"; echo "vcfs fail=$fail"; exit 1; }
done

# what the run executed, from its trace: the skipped stages must be absent
ran=$(awk -F'\t' 'NR>1 {sub(/ \(.*/, "", $4); print $4}' "$V/nextflow_trace.txt" | sort -u | tr '\n' ' ')
if grep -qwE "vg_call|graph_align_reads|make_graph" <<<"$ran"; then
  bad "the run aligned or called anyway: $ran"
else
  ok "no make_graph, graph_align_reads or vg_call in the trace"
fi

q() { bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$1/$M"; }
samples_s=$(bcftools query -l "$S/$M" | tr '\n' ' '); samples_v=$(bcftools query -l "$V/$M" | tr '\n' ' ')
[[ "$samples_s" == "$samples_v" ]] && ok "same samples: $samples_v" \
                                   || bad "samples differ: spine $samples_s, vcfs $samples_v"
n_s=$(q "$S" | wc -l); n_v=$(q "$V" | wc -l)
called=$(q "$V" | awk -F'\t' '{for (i = 5; i <= NF; i++) if ($i != "0/0" && $i != "./." && $i != ".") {n++; break}} END {print n + 0}')
[[ $n_s -gt 0 && $called -gt 0 ]] && ok "$n_v records, $called with a non-reference genotype" \
                                  || bad "nothing to compare: $n_s records in the spine, $called called here"
d=$(diff <(q "$S") <(q "$V") | grep -c '^[<>]')
[[ $d -eq 0 && $n_s -eq $n_v ]] && ok "every record and genotype matches the spine's" \
                                || bad "$d lines differ ($n_s records in the spine, $n_v here)"

echo "vcfs fail=$fail"; exit $fail
