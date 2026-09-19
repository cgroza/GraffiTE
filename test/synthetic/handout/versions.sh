#!/usr/bin/env bash
# Which versions of the tools the image carried.
#
# bundle_results.sh writes VERSIONS.txt only when GRAFFITE_SIF names a file, and
# INPUTS.env leaves GRAFFITE_SIF empty on purpose, so a correctly configured run
# produces no VERSIONS.txt at all. The image Nextflow pulled sits in the
# singularity cache. Probe that, and leave INPUTS.env alone.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
# shellcheck disable=SC1091
source ./INPUTS.env
H="$PWD"
OUT="${1:-$H/VERSIONS.txt}"

RUNNER="$(command -v apptainer || command -v singularity || true)"
[[ -n "$RUNNER" ]] || { echo "no apptainer/singularity on PATH -- run this inside an allocation" >&2; exit 1; }

SIF="${GRAFFITE_SIF:-}"
if [[ -z "$SIF" || ! -f "$SIF" ]]; then
  for d in "${NXF_SINGULARITY_CACHEDIR:-$H/singularity_cache}" "$H/work/singularity" \
           "$H/singularity_cache"; do
    [[ -d "$d" ]] || continue
    cand=$(find "$d" -maxdepth 1 -name '*graffite*' \( -name '*.img' -o -name '*.sif' \) 2>/dev/null | head -1)
    [[ -n "$cand" ]] && { SIF="$cand"; break; }
  done
fi
[[ -n "$SIF" && -f "$SIF" ]] || { echo "no graffite image found in the singularity cache" >&2; exit 1; }

{
  echo "image     : $SIF"
  echo "size      : $(du -h "$SIF" | cut -f1)"
  echo "mtime     : $(date -u -r "$SIF" +%Y-%m-%dT%H:%M:%SZ)"
  echo "sha256    : $(sha256sum "$SIF" | cut -c1-64)"
  echo "runner    : $("$RUNNER" --version 2>&1 | head -1)"
  echo "probed    : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo
  for t in RepeatMasker ProcessRepeats ultra minimap2 samtools bcftools bedtools vg \
           sniffles svim-asm truvari GraphAligner PanGenie PanGenie-index winnowmap \
           tabix bgzip python3 Rscript blastn makeblastdb; do
    printf '%-16s ' "$t"
    "$RUNNER" exec "$SIF" bash -lc "command -v $t >/dev/null 2>&1 || { echo '(absent)'; exit 0; }
      ( $t --version 2>&1 || $t -v 2>&1 || $t -h 2>&1 || true ) | grep -m1 -iE '[0-9]+\.[0-9]' || echo '(no version string)'" 2>/dev/null \
      || echo "(probe failed)"
  done
  echo
  echo "-- R packages --"
  "$RUNNER" exec "$SIF" Rscript -e 'for (p in c("optparse","dplyr","stringr","tidyr","readr","vcfR")) cat(sprintf("%-10s %s\n", p, tryCatch(as.character(packageVersion(p)), error=function(e) "(absent)")))' 2>/dev/null
  echo
  echo "-- python packages --"
  "$RUNNER" exec "$SIF" python3 -c 'import importlib.metadata as m
for p in ("pysam","vcfpy","truvari","pandas","numpy"):
    try: print("%-10s %s" % (p, m.version(p)))
    except Exception: print("%-10s (absent)" % p)' 2>/dev/null
} > "$OUT"
echo "wrote $OUT"
