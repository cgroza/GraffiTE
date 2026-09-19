#!/usr/bin/env bash
# Pack what came back into one archive.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
# shellcheck disable=SC1091
[[ -f INPUTS.env ]] && source ./INPUTS.env
WORKDIR="${WORKDIR:-$PWD}"
BUNDLE="synthetic_results_$(date +%Y%m%d_%H%M)"
rm -rf "$BUNDLE" && mkdir -p "$BUNDLE"

# Per-run: the published outputs that are small, plus the logs.
for d in "$WORKDIR"/runs/*/; do
  [[ -d "$d" ]] || continue
  r="$(basename "$d")"
  mkdir -p "$BUNDLE/$r"
  for f in 3_TSD_search/pangenome.vcf 3_TSD_search/pangenome.trusted.vcf \
           3_TSD_search/pangenome.human.vcf 3_TSD_search/human_filter_summary.txt \
           3_TSD_search/TSD_summary.txt \
           4_Genotyping/genotyping_record_audit.tsv \
           4_Genotyping/pangenie_graph_variants.tsv \
           4_Genotyping/GraffiTE.merged.genotypes.presence-absence_trusted.tsv \
           hervk_loci.tsv hervk_reconciliation_report.md \
           nextflow_trace.txt nextflow_report.html run.log published.txt; do
    [[ -f "$d/$f" ]] && cp "$d/$f" "$BUNDLE/$r/$(basename "$f")"
  done
  # Headers only for the big VCFs; the records are not what a reviewer reads.
  for v in "$d"/4_Genotyping/*.vcf.gz; do
    [[ -f "$v" ]] || continue
    { echo "== $(basename "$v")  records=$(bcftools view -H "$v" 2>/dev/null | wc -l)"
      bcftools view -h "$v" 2>/dev/null | grep -E '^##(INFO|FORMAT|source|GraffiTE)' || true
    } >> "$BUNDLE/$r/vcf_headers.txt"
  done
  ( cd "$d" && find . -type f | sort ) > "$BUNDLE/$r/published_files.txt" 2>/dev/null || true
done

cp OBSERVED.tsv CALIBRATION.log ASSERTIONS.log MATRIX.log "$BUNDLE/" 2>/dev/null || true
[[ -f "$WORKDIR/build/truth.tsv" ]] && cp "$WORKDIR/build/truth.tsv" "$BUNDLE/"
[[ -f "$WORKDIR/build/MANIFEST.sha256" ]] && cp "$WORKDIR/build/MANIFEST.sha256" "$BUNDLE/"

# What the tools actually were. container.md currently guesses at versions.
if [[ -n "${GRAFFITE_SIF:-}" && -f "${GRAFFITE_SIF}" ]]; then
  R=$(command -v apptainer || command -v singularity || true)
  if [[ -n "$R" ]]; then
    { for t in RepeatMasker ultra minimap2 samtools bcftools vg sniffles svim-asm truvari GraphAligner PanGenie winnowmap; do
        printf '%-14s ' "$t"
        "$R" exec "$GRAFFITE_SIF" bash -lc "command -v $t >/dev/null && ($t --version 2>&1 || $t -v 2>&1 || true) | head -1" 2>/dev/null || echo "(absent)"
      done; } > "$BUNDLE/VERSIONS.txt"
  fi
fi

{
  echo "date     : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host     : $(hostname)"
  echo "nextflow : $(nextflow -v 2>/dev/null | head -1)"
  echo "project  : ${PROJECT:-?} -r ${REVISION:-?}"
  echo "commit   : $(cd "${NXF_ASSETS:-$HOME/.nextflow/assets}/${PROJECT:-cgroza/GraffiTE}" 2>/dev/null && git rev-parse HEAD 2>/dev/null || echo '?')"
  echo "profile  : ${PROFILE:-?}  cpus=${CPUS:-?}  mem=${MEM_GB:-?}G"
  echo "sif      : ${GRAFFITE_SIF:-<pulled by nextflow>}"
  echo "seed     : ${SEED:-?}"
} > "$BUNDLE/RUN_INFO.txt"

tar czf "${BUNDLE}.tar.gz" "$BUNDLE" && rm -rf "$BUNDLE"
echo "wrote ${BUNDLE}.tar.gz ($(du -h "${BUNDLE}.tar.gz" | cut -f1))"
