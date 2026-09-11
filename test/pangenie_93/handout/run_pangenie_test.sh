#!/usr/bin/env bash
# #93 PanGenie test: genotype the reporter's pangenome.vcf with PanGenie,
# starting from --graffite_vcf (discovery is skipped), then run the assertions.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
[[ -f INPUTS.env ]] || { echo "INPUTS.env not found -- run ./bootstrap.sh first" >&2; exit 1; }
# shellcheck disable=SC1091
source ./INPUTS.env

for v in GRAFFITE_VCF REFERENCE READS_1; do
  [[ -n "${!v:-}" ]] || { echo "$v is empty in INPUTS.env" >&2; exit 1; }
done
OUTDIR="${OUTDIR:-pangenie_93_run}"
PROFILE="${PROFILE:-cluster}"
CPUS="${CPUS:-16}"
PANGENIE_MEMORY="${PANGENIE_MEMORY:-120G}"
PANGENIE_TIME="${PANGENIE_TIME:-24h}"
SAMPLE_1="${SAMPLE_1:-S1}"
SAMPLE_2="${SAMPLE_2:-S2}"
REVISION="${REVISION:-fix/pangenie-dup-ids-93}"
PROJECT="${PROJECT:-cgroza/GraffiTE}"
mkdir -p "$OUTDIR"
OUTDIR_ABS="$(cd "$OUTDIR" && pwd)"

# Two samples so that the run shows pangenie starting once per sample. With no
# READS_2 the same reads go in twice; GraffiTE sees two samples either way.
READS_CSV="$OUTDIR_ABS/reads.csv"
{
  echo "path,sample,type"
  echo "${READS_1},${SAMPLE_1},short"
  echo "${READS_2:-$READS_1},${SAMPLE_2},short"
} > "$READS_CSV"
echo "reads    : $READS_CSV"
sed 's/^/           /' "$READS_CSV"

# Resume safety. Nextflow's task hash does not cover bin/, so after a new
# commit -resume could reuse a pangenie_index built by the old
# pangenie_graph_vcf.py. Drop -resume whenever the revision moves.
ASSET="${NXF_ASSETS:-$HOME/.nextflow/assets}/$PROJECT"
COMMIT="$(cd "$ASSET" 2>/dev/null && git rev-parse HEAD 2>/dev/null || echo unknown)"
STAMP="$OUTDIR/.last_commit"
RESUME_ARG=(-resume)
if [[ -f "$STAMP" && "$(cat "$STAMP")" != "$COMMIT" ]]; then
  RESUME_ARG=()
  echo "!! pipeline moved $(cut -c1-8 "$STAMP") -> ${COMMIT:0:8}: running WITHOUT -resume"
elif [[ ! -f "$STAMP" && -d "$OUTDIR/4_Genotyping" ]]; then
  RESUME_ARG=()
  echo "!! existing output with no revision stamp: running without -resume"
fi

echo "project  : $PROJECT -r $REVISION (${COMMIT:0:8})"
echo "nextflow : $(nextflow -v 2>/dev/null | head -1)"
echo "out      : $OUTDIR"

# -latest re-pulls the branch tip so a fix pushed upstream is picked up.
nextflow run "$PROJECT" -r "$REVISION" -latest \
    --graffite_vcf    "$GRAFFITE_VCF" \
    --reference       "$REFERENCE" \
    --graph_method    pangenie \
    --genotype_with   "$READS_CSV" \
    --cores           "$CPUS" \
    --pangenie_memory "$PANGENIE_MEMORY" \
    --pangenie_time   "$PANGENIE_TIME" \
    --out             "$OUTDIR" \
    -profile          "$PROFILE" \
    ${GRAFFITE_SIF:+-with-singularity "$GRAFFITE_SIF"} \
    -with-report "$OUTDIR/nextflow_report.html" \
    -with-trace  "$OUTDIR/nextflow_trace.txt" \
    "${RESUME_ARG[@]}"

echo "$COMMIT" > "$STAMP"

echo
echo "== assertions =="
python3 ./assert_pangenie_test.py \
    --outdir "$OUTDIR" \
    --graffite-vcf "$GRAFFITE_VCF" \
    --samples "$SAMPLE_1,$SAMPLE_2" \
    ${READS_2:+--distinct-reads} \
  | tee "$OUTDIR/pangenie_93_assertions.log"
exit "${PIPESTATUS[0]}"
