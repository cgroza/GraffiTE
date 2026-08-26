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

# Stage E runs against an existing genotyped VCF rather than re-genotyping.
STAGE_E=()
if [[ -n "${GENOTYPED_VCF:-}" ]]; then
  [[ -f "$GENOTYPED_VCF" ]] || { echo "GENOTYPED_VCF=$GENOTYPED_VCF not found" >&2; exit 1; }
  STAGE_E=(--hervk_reconcile_vcf "$GENOTYPED_VCF")
  echo "stage E   : consolidating against $GENOTYPED_VCF"
else
  echo "stage E   : skipped (GENOTYPED_VCF not set)"
fi

mkdir -p "$OUTDIR"

# Resume safety.
#
# Nextflow's task hash covers the process script text and the input files. It
# does NOT cover the contents of bin/, which is staged onto PATH. So when the
# HERV-K scripts change but the process block does not, -resume happily reuses
# the previous task output and the run reports stale results while appearing to
# succeed. Stamp the revision and drop -resume whenever it moves.
ASSET="${NXF_ASSETS:-$HOME/.nextflow/assets}/$PROJECT"
COMMIT="$(cd "$ASSET" 2>/dev/null && git rev-parse HEAD 2>/dev/null || echo unknown)"
STAMP="$OUTDIR/.last_commit"
RESUME_ARG=(-resume)
if [[ -f "$STAMP" ]]; then
  PREV="$(cat "$STAMP")"
  if [[ "$PREV" != "$COMMIT" ]]; then
    RESUME_ARG=()
    echo "!! pipeline moved ${PREV:0:8} -> ${COMMIT:0:8}"
    echo "!! running WITHOUT -resume: Nextflow does not hash bin/, so resuming"
    echo "!! here would reuse the previous HERV-K output and hide the change."
  fi
elif [[ -d "$OUTDIR/3_TSD_search" ]]; then
  RESUME_ARG=()
  echo "!! existing output with no revision stamp — running without -resume"
fi

echo "project  : $PROJECT -r $REVISION (${COMMIT:0:8})"
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
    "${STAGE_E[@]}" \
    --out        "$OUTDIR" \
    --cores      "$CPUS" \
    -profile     "$PROFILE" \
    ${GRAFFITE_SIF:+-with-singularity "$GRAFFITE_SIF"} \
    -with-report "$OUTDIR/nextflow_report.html" \
    -with-trace  "$OUTDIR/nextflow_trace.txt" \
    "${RESUME_ARG[@]}"

echo "$COMMIT" > "$STAMP"

echo
echo "== assertions =="
python3 ./assert_hervk_test.py --outdir "$OUTDIR" | tee "$OUTDIR/hervk_assertions.log"
