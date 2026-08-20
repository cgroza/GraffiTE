#!/usr/bin/env bash
# HERV-K v2 discovery test: annotate an existing PAV call set, no genotyping.
#
# Graph genotyping is NOT run and CANNOT be run for this cohort -- there are no
# raw reads. --genotype false is deliberate; do not remove it.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
[[ -f INPUTS.env ]] || { echo "INPUTS.env not found — run ./bootstrap.sh first" >&2; exit 1; }
# shellcheck disable=SC1091
source ./INPUTS.env

for v in REFERENCE TE_LIBRARY; do
  [[ -n "${!v:-}" ]] || { echo "$v is empty in INPUTS.env" >&2; exit 1; }
done
if [[ -n "${RM_DIR:-}" ]]; then
  ENTRY=(--RM_dir "$RM_DIR")
  echo "entry    : --RM_dir $RM_DIR  (RepeatMasker skipped)"
elif [[ -n "${PAV_VCF:-}" ]]; then
  ENTRY=(--vcf "$PAV_VCF")
  echo "entry    : --vcf $PAV_VCF  (full re-mask, slow path)"
else
  echo "set RM_DIR or PAV_VCF in INPUTS.env" >&2; exit 1
fi
OUTDIR="${OUTDIR:-hervk_v2_run}"
PROFILE="${PROFILE:-cluster}"
CPUS="${CPUS:-8}"
REVISION="${REVISION:-v1.1dev-hervk-v2}"
PROJECT="${PROJECT:-cgroza/GraffiTE}"

mkdir -p "$OUTDIR"
echo "project  : $PROJECT -r $REVISION"
echo "out      : $OUTDIR"
echo "profile  : $PROFILE"

# -latest re-pulls the branch tip, so a fix pushed upstream is picked up
# instead of silently running a stale cached copy.
nextflow run "$PROJECT" -r "$REVISION" -latest \
    "${ENTRY[@]}" \
    --reference  "$REFERENCE" \
    --TE_library "$TE_LIBRARY" \
    --human      true \
    --genotype   false \
    --out        "$OUTDIR" \
    --cores      "$CPUS" \
    -profile     "$PROFILE" \
    ${GRAFFITE_SIF:+-with-singularity "$GRAFFITE_SIF"} \
    -with-report "$OUTDIR/nextflow_report.html" \
    -with-trace  "$OUTDIR/nextflow_trace.txt" \
    -resume

echo
echo "== assertions =="
python3 ./assert_hervk_test.py --outdir "$OUTDIR" | tee "$OUTDIR/hervk_assertions.log"
