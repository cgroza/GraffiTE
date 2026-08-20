#!/usr/bin/env bash
# Pack the HERV-K outputs into one archive to bring back for analysis.
set -euo pipefail
OUTDIR="${1:-${OUTDIR:-hervk_v2_run}}"
STAMP="$(date +%Y%m%d)"
BUNDLE="hervk_v2_results_${STAMP}"
rm -rf "$BUNDLE" && mkdir -p "$BUNDLE"

for f in hervk_calls.tsv hervk_loci.tsv hervk_arch.tsv hervk_refstate.tsv \
         hervk_polymorphism_summary.md pangenome.human.vcf \
         pangenome.presence-absence_human.tsv human_filter_summary.txt; do
  [[ -f "$OUTDIR/3_TSD_search/$f" ]] && cp "$OUTDIR/3_TSD_search/$f" "$BUNDLE/" \
    || echo "  (missing: $f)"
done
for f in hervk_assertions.log nextflow_report.html nextflow_trace.txt; do
  [[ -f "$OUTDIR/$f" ]] && cp "$OUTDIR/$f" "$BUNDLE/"
done

# Provenance: without this the numbers cannot be tied to a commit.
{
  echo "date        : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host        : $(hostname)"
  GT_DIR="${GT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
  echo "graffite    : $GT_DIR"
  echo "branch      : $(cd "$GT_DIR" && git rev-parse --abbrev-ref HEAD 2>/dev/null || echo '?')"
  echo "commit      : $(cd "$GT_DIR" && git rev-parse HEAD 2>/dev/null || echo '?')"
  echo "dirty       : $(cd "$GT_DIR" && git status --porcelain 2>/dev/null | wc -l | tr -d ' ') modified files"
  echo "pav_vcf     : ${PAV_VCF:-<unset>}"
  echo "reference   : ${REFERENCE:-<unset>}"
  echo "te_library  : ${TE_LIBRARY:-<unset>}"
} > "$BUNDLE/RUN_INFO.txt"

# gzip the VCF only if it is big enough to matter
[[ -f "$BUNDLE/pangenome.human.vcf" ]] && gzip -f "$BUNDLE/pangenome.human.vcf"

tar czf "${BUNDLE}.tar.gz" "$BUNDLE" && rm -rf "$BUNDLE"
echo "wrote ${BUNDLE}.tar.gz ($(du -h "${BUNDLE}.tar.gz" | cut -f1))"
