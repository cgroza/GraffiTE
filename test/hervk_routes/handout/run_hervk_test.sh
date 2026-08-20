#!/usr/bin/env bash
# HERV-K v2 discovery test: annotate an existing PAV call set, no genotyping.
#
# Graph genotyping is NOT run and CANNOT be run for this cohort -- there are no
# raw reads. --genotype false is deliberate; do not remove it.
set -euo pipefail

: "${PAV_VCF:?set PAV_VCF to the merged PAV call set}"
: "${REFERENCE:?set REFERENCE to the CHM13v2 FASTA}"
: "${TE_LIBRARY:?set TE_LIBRARY to the RepeatMasker library FASTA}"
OUTDIR="${OUTDIR:-hervk_v2_run}"
GT_DIR="${GT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
PROFILE="${PROFILE:-standard}"
CPUS="${CPUS:-8}"

echo "GraffiTE : $GT_DIR  ($(cd "$GT_DIR" && git rev-parse --abbrev-ref HEAD 2>/dev/null || echo '?'))"
echo "out      : $OUTDIR"

nextflow run "$GT_DIR/main.nf" \
    --vcf        "$PAV_VCF" \
    --reference  "$REFERENCE" \
    --TE_library "$TE_LIBRARY" \
    --human      true \
    --genotype   false \
    --out        "$OUTDIR" \
    --cores      "$CPUS" \
    -profile     "$PROFILE" \
    -with-report "$OUTDIR/nextflow_report.html" \
    -with-trace  "$OUTDIR/nextflow_trace.txt" \
    -resume

echo
echo "== assertions =="
python3 "$(dirname "${BASH_SOURCE[0]}")/assert_hervk_test.py" --outdir "$OUTDIR" \
    | tee "$OUTDIR/hervk_assertions.log"
