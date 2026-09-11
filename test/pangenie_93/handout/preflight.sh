#!/usr/bin/env bash
# Pre-flight for the #93 PanGenie test. Checks inputs, tools and the cached
# revision before committing a cluster allocation.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
if [[ -f INPUTS.env ]]; then
  # shellcheck disable=SC1091
  source ./INPUTS.env
else
  echo "INPUTS.env not found -- run ./bootstrap.sh first" >&2; exit 1
fi
fail=0
ok(){   printf '  [ ok ] %s\n' "$1"; }
bad(){  printf '  [FAIL] %s\n' "$1"; fail=1; }
warn(){ printf '  [warn] %s\n' "$1"; }
cat_any(){ case "$1" in *.gz) gzip -cd "$1" ;; *) cat "$1" ;; esac; }

echo "== inputs =="
for var in GRAFFITE_VCF REFERENCE READS_1; do
  val="${!var:-}"
  if   [[ -z "$val"   ]]; then bad "$var is not set"
  elif [[ ! -f "$val" ]]; then bad "$var=$val does not exist"
  elif [[ "$val" != /* ]]; then bad "$var=$val is relative -- use an absolute path"
  else ok "$var=$val ($(du -h "$val" | cut -f1))"; fi
done
if [[ -n "${READS_2:-}" ]]; then
  [[ -f "$READS_2" && "$READS_2" == /* ]] && ok "READS_2=$READS_2" || bad "READS_2=$READS_2 missing or relative"
else
  warn "READS_2 empty -- READS_1 is genotyped twice, as ${SAMPLE_1:-S1} and ${SAMPLE_2:-S2}"
fi
[[ -n "${REFERENCE:-}" && -f "${REFERENCE}.fai" ]] \
  && ok "reference index ${REFERENCE}.fai" \
  || warn "no ${REFERENCE:-<REFERENCE>}.fai -- contig check skipped; the pipeline builds one"

echo "== pangenome VCF =="
if [[ -f "${GRAFFITE_VCF:-}" ]]; then
  # Expected for the #93 file: 51910 records, 31 duplicate CHROM/POS/REF/ALT
  # keys, 0 missing IDs. A deduplicated copy (0 duplicates) makes the test
  # vacuous for the fix, so it is flagged.
  read -r n_rec n_dup n_dot n_multi < <(cat_any "$GRAFFITE_VCF" | awk -F'\t' '
    /^#/ {next}
    {n++; k=$1"\t"$2"\t"$4"\t"$5; c[k]++; if(c[k]==2) d++; if($3==".") dot++; if($5 ~ /,/) m++}
    END {print n+0, d+0, dot+0, m+0}')
  ok "$n_rec records, $n_dup duplicate CHROM/POS/REF/ALT keys, $n_dot missing IDs, $n_multi multi-ALT records"
  [[ "$n_dup" -gt 0 ]] || warn "no duplicate records -- is this the deduplicated copy? The fix is only exercised by the original."
  [[ "$n_rec" -eq 51910 && "$n_dup" -eq 31 ]] || warn "counts differ from the #93 file (51910 records, 31 duplicates) -- assertions derived from it will say so"
  if [[ -f "${REFERENCE:-}.fai" ]]; then
    missing=$(cat_any "$GRAFFITE_VCF" | awk -F'\t' '!/^#/ {print $1}' | sort -u \
              | comm -23 - <(cut -f1 "${REFERENCE}.fai" | sort -u) | wc -l | tr -d ' ')
    [[ "$missing" -eq 0 ]] && ok "every VCF contig is in the reference" \
      || bad "$missing VCF contigs are not in ${REFERENCE}.fai -- wrong reference?"
  fi
fi

echo "== reads =="
for var in READS_1 READS_2; do
  val="${!var:-}"; [[ -z "$val" || ! -f "$val" ]] && continue
  first=$(cat_any "$val" 2>/dev/null | head -c1)
  [[ "$first" == "@" ]] && ok "$var looks like FASTQ" || bad "$var does not start with '@' -- not a FASTQ?"
done

echo "== container =="
PROFILE="${PROFILE:-cluster}"
RUNNER=""
command -v apptainer   >/dev/null && RUNNER=apptainer
[[ -z "$RUNNER" ]] && command -v singularity >/dev/null && RUNNER=singularity
if [[ -n "${GRAFFITE_SIF:-}" ]]; then
  [[ -f "$GRAFFITE_SIF" ]] && ok "container $GRAFFITE_SIF" || bad "GRAFFITE_SIF=$GRAFFITE_SIF missing"
  if [[ -n "$RUNNER" && -f "$GRAFFITE_SIF" ]]; then
    for t in PanGenie PanGenie-index bcftools python3; do
      "$RUNNER" exec "$GRAFFITE_SIF" which "$t" >/dev/null 2>&1 \
        && ok "container has $t" || bad "container is missing $t"
    done
    "$RUNNER" exec "$GRAFFITE_SIF" python3 -c 'import pyfaidx' >/dev/null 2>&1 \
      && ok "container python3 has pyfaidx (merge_vcfs.py)" || bad "container python3 lacks pyfaidx"
  fi
else
  warn "GRAFFITE_SIF not set -- Nextflow will pull the container for -profile $PROFILE"
  [[ -n "$RUNNER" ]] && ok "$RUNNER present" \
    || warn "no apptainer/singularity on this host -- fine if the compute nodes have it"
fi

echo "== pipeline revision =="
PROJECT="${PROJECT:-cgroza/GraffiTE}"
REVISION="${REVISION:-fix/pangenie-dup-ids-93}"
GT_DIR="${NXF_ASSETS:-$HOME/.nextflow/assets}/$PROJECT"
if [[ ! -d "$GT_DIR" ]]; then
  bad "$GT_DIR not cached -- run ./bootstrap.sh"
else
  BR="$(cd "$GT_DIR" && git rev-parse --abbrev-ref HEAD 2>/dev/null || echo '?')"
  [[ "$BR" == "$REVISION" ]] && ok "cached revision $BR ($(cd "$GT_DIR" && git rev-parse --short HEAD))" \
    || bad "cached revision is $BR, expected $REVISION -- rerun ./bootstrap.sh"
  for f in main.nf module/main.nf bin/pangenie_graph_vcf.py bin/merge_vcfs.py test/pangenie_index/test_graph_vcf.sh; do
    [[ -f "$GT_DIR/$f" ]] && ok "$f" || bad "$GT_DIR/$f missing"
  done
  [[ -f "$GT_DIR/panmethyl/module/main.nf" ]] && ok "panmethyl submodule present" \
    || bad "panmethyl/module/main.nf missing -- the include in main.nf will fail (nextflow pull should fetch submodules)"
  # The graph-input unit test needs bcftools and pyfaidx. Run it in the
  # container when one is given, otherwise on the host (it skips itself when
  # a tool is missing).
  if [[ -n "$RUNNER" && -f "${GRAFFITE_SIF:-}" ]]; then
    out=$("$RUNNER" exec "$GRAFFITE_SIF" bash "$GT_DIR/test/pangenie_index/test_graph_vcf.sh" 2>&1); rc=$?
  else
    out=$(bash "$GT_DIR/test/pangenie_index/test_graph_vcf.sh" 2>&1); rc=$?
  fi
  if [[ "$out" == *"[skip]"* ]]; then warn "graph-input unit test skipped: $(echo "$out" | grep skip | head -1)"
  elif [[ $rc -eq 0 ]]; then ok "graph-input unit test passes ($(echo "$out" | grep -c '\[ ok \]') checks)"
  else bad "graph-input unit test FAILED -- do not run the pipeline, report this:"; echo "$out" | grep FAIL; fi
fi

echo "== nextflow =="
if command -v nextflow >/dev/null; then
  NXF_V=$(nextflow -v 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
  ok "nextflow $NXF_V"
  # Nextflow 26.04 parses with the strict syntax by default; GraffiTE did not
  # compile under it before this branch. nextflow lint exists from 25.04.
  if [[ -d "$GT_DIR" ]] && nextflow lint -h >/dev/null 2>&1; then
    lint=$(cd "$GT_DIR" && nextflow lint -o concise main.nf module/main.nf nextflow.config 2>&1)
    nerr=$(echo "$lint" | grep -c '^Error')
    [[ "$nerr" -eq 0 ]] && ok "nextflow lint: 0 errors in main.nf, module/main.nf, nextflow.config" \
      || { bad "nextflow lint: $nerr errors"; echo "$lint" | grep '^Error' | head -10; }
  else
    warn "nextflow lint not available in $NXF_V -- strict-syntax check skipped"
  fi
else
  bad "nextflow not on PATH (module load nextflow?)"
fi

echo
[[ $fail -eq 0 ]] && echo "pre-flight PASSED" || echo "pre-flight FAILED -- fix the above before running"
exit $fail
