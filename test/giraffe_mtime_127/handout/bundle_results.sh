#!/usr/bin/env bash
# Pack the #127 test results into one archive to bring back.
set -uo pipefail
# shellcheck disable=SC1091
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"
BUNDLE="$HANDOUT/issue127_results_$(date +%Y%m%d)"
rm -rf "$BUNDLE"; mkdir -p "$BUNDLE"

[[ -f "$WORKDIR/VERDICT.txt" ]] && cp "$WORKDIR/VERDICT.txt" "$BUNDLE/"
if [[ -d "$WORKDIR/partA" ]]; then
  mkdir -p "$BUNDLE/partA"
  cp "$WORKDIR"/partA/{results.tsv,inodes.tsv,meta.txt} "$BUNDLE/partA/" 2>/dev/null
  cp -r "$WORKDIR/partA/logs" "$BUNDLE/partA/" 2>/dev/null
fi
if [[ -d "$WORKDIR/partB" ]]; then
  mkdir -p "$BUNDLE/partB"
  cp "$WORKDIR"/partB/{results.tsv,index_mtimes.tsv} "$BUNDLE/partB/" 2>/dev/null
  for d in "$WORKDIR"/partB/*/; do
    c=$(basename "$d"); [[ "$c" == work_* ]] && continue
    mkdir -p "$BUNDLE/partB/$c"
    find "$d" -maxdepth 1 -type f \( -name '*.txt' -o -name '*.log' -o -name '*.err' -o -name '*.sh' -o -name '*.exitcode' \) \
      -exec cp {} "$BUNDLE/partB/$c/" \;
  done
fi
{
  echo "date     : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host     : $(hostname)"
  echo "nextflow : $(nextflow -v 2>/dev/null | head -1)"
  echo "project  : $PROJECT -r $REVISION"
  echo "commit   : $(cd "${NXF_ASSETS:-$HOME/.nextflow/assets}/$PROJECT" 2>/dev/null && git rev-parse HEAD 2>/dev/null || echo '?')"
  echo "image    : $IMG"
  echo "inputs   : $REFERENCE | $GRAFFITE_VCF | $READS_CSV"
  echo "scratch  : A=$SCRATCH_DIRS  B=$NF_SCRATCH  SLURM_TMPDIR=${SLURM_TMPDIR:-<unset>}"
} > "$BUNDLE/RUN_INFO.txt"
tar czf "$BUNDLE.tar.gz" -C "$HANDOUT" "$(basename "$BUNDLE")" && echo "bundle: $BUNDLE.tar.gz"
