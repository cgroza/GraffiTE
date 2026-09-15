#!/usr/bin/env bash
# Pack the #93 PanGenie test outputs into one archive to bring back.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
# shellcheck disable=SC1091
[[ -f INPUTS.env ]] && source ./INPUTS.env
OUTDIR="${1:-${OUTDIR:-pangenie_93_run}}"
BUNDLE="pangenie_93_results_$(date +%Y%m%d)"
rm -rf "$BUNDLE" && mkdir -p "$BUNDLE"

for f in "$OUTDIR"/4_Genotyping/pangenie_graph_variants.tsv \
         "$OUTDIR"/4_Genotyping/GraffiTE.merged.genotypes.vcf.gz \
         "$OUTDIR"/4_Genotyping/*_genotyping.vcf.gz \
         "$OUTDIR"/pangenie_93_assertions.log \
         "$OUTDIR"/nextflow_report.html \
         "$OUTDIR"/nextflow_trace.txt \
         "$OUTDIR"/reads.csv; do
  [[ -f "$f" ]] && cp "$f" "$BUNDLE/" || echo "  (missing: $f)"
done

# The pangenie_index task logs hold the pangenie_graph_vcf.py counts (records,
# duplicates merged, IDs replaced, alleles left out of the graph) and the
# merge_vcfs.py messages. Find the task directory through the trace.
WORK="${NXF_WORK:-work}"
if [[ -f "$OUTDIR/nextflow_trace.txt" ]]; then
  awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) col[$i]=i; next}
              $col["name"] ~ /^pangenie_index/ || $col["name"] ~ /^pangenie / {print $col["name"] "\t" $col["hash"] "\t" $col["status"] "\t" $col["exit"]}' \
      "$OUTDIR/nextflow_trace.txt" > "$BUNDLE/pangenie_tasks.tsv" || true
  while IFS=$'\t' read -r name hash status rc; do
    d=$(ls -d "$WORK/$hash"* 2>/dev/null | head -1) || true
    [[ -n "$d" ]] || continue
    tag=$(echo "$name" | tr -c 'A-Za-z0-9_\n' '_')
    for x in .command.sh .command.err .command.out .exitcode; do
      [[ -f "$d/$x" ]] && cp "$d/$x" "$BUNDLE/${tag}${x}"
    done
  done < "$BUNDLE/pangenie_tasks.tsv"
fi

GT_DIR="${NXF_ASSETS:-$HOME/.nextflow/assets}/${PROJECT:-cgroza/GraffiTE}"
{
  echo "date        : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host        : $(hostname)"
  echo "nextflow    : $(nextflow -v 2>/dev/null | head -1)"
  echo "project     : ${PROJECT:-cgroza/GraffiTE} -r ${REVISION:-?}"
  echo "commit      : $(cd "$GT_DIR" 2>/dev/null && git rev-parse HEAD 2>/dev/null || echo '?')"
  echo "graffite_vcf: ${GRAFFITE_VCF:-<unset>}"
  echo "reference   : ${REFERENCE:-<unset>}"
  echo "reads_1     : ${READS_1:-<unset>}"
  echo "reads_2     : ${READS_2:-<same as reads_1>}"
  echo "sif         : ${GRAFFITE_SIF:-<pulled by nextflow>}"
} > "$BUNDLE/RUN_INFO.txt"

tar czf "${BUNDLE}.tar.gz" "$BUNDLE" && rm -rf "$BUNDLE"
echo "wrote ${BUNDLE}.tar.gz ($(du -h "${BUNDLE}.tar.gz" | cut -f1))"
