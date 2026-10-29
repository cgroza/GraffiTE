#!/usr/bin/env bash
# Part B: the pipeline itself, --graph_method giraffe from --graffite_vcf with
# long reads, as in issue #127. Three kinds of cell, each in its own work
# directory and without -resume:
#
#   noscratch        no scratch (what the synthetic suite has always run)
#   scratch_<n>      process.scratch as in -profile cluster, REPS_B times
#   scratch_move     the same, with stageOutMode 'move' on make_graph
#
# Writes $WORKDIR/partB/{results.tsv,index_mtimes.tsv} and one directory per
# cell with its run.log, nextflow.log and trace.
set -uo pipefail
export LC_ALL=C
# shellcheck disable=SC1091
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

for f in "$REFERENCE" "$GRAFFITE_VCF" "$READS_CSV"; do
  [[ -f "$f" ]] || { echo "missing input: $f" >&2; exit 1; }
done
command -v nextflow >/dev/null || { echo "nextflow not on PATH (module load nextflow)" >&2; exit 1; }
ensure_image || exit 1

OUT="$WORKDIR/partB"
mkdir -p "$OUT"
RESULTS="$OUT/results.tsv"
MTIMES="$OUT/index_mtimes.tsv"
printf 'cell\toverlay\tnf_rc\ttask\tstatus\texit\tpredicted\tobserved\tmessage\n' > "$RESULTS"
printf 'cell\tfile\tmtime\tsize\n' > "$MTIMES"

export GT_CPUS="$CPUS" GT_MEM_GB="$MEM_GB" MTIME_NF_SCRATCH="$NF_SCRATCH"
SIF_ARG=(); [[ -n "${GRAFFITE_SIF:-}" ]] && SIF_ARG=(-with-singularity "$GRAFFITE_SIF")
TMP_ARG=(); [[ -n "${CONTAINER_TMP:-}" ]] && TMP_ARG=(--container_tmp "$CONTAINER_TMP")

workdir_of() {  # workdir_of <work> <trace hash, e.g. ab/123456>
  ls -d "$1/$2"* 2>/dev/null | head -1
}

run_cell() {  # run_cell <cell> <overlay config or "none">
  local cell="$1" overlay="$2"
  local out="$OUT/$cell" work="$OUT/work_$cell"
  rm -rf "$out" "$work"; mkdir -p "$out"
  local cfg=(-c "$HANDOUT/local.config")
  [[ "$overlay" != none ]] && cfg+=(-c "$HANDOUT/$overlay")
  echo "== $cell ($overlay) =="
  ( set -x
    nextflow -log "$out/nextflow.log" run "$PROJECT" -r "$REVISION" \
      -profile standard "${cfg[@]}" -w "$work" \
      --graffite_vcf  "$GRAFFITE_VCF" \
      --reference     "$REFERENCE" \
      --genotype_with "$READS_CSV" \
      --graph_method  giraffe \
      --cores         "$CPUS" \
      --out           "$out" \
      "${SIF_ARG[@]}" "${TMP_ARG[@]}" \
      -with-trace "$out/nextflow_trace.txt"
  ) > "$out/run.log" 2>&1
  local nf_rc=$?

  local trace="$out/nextflow_trace.txt"
  if [[ ! -f "$trace" ]]; then
    printf '%s\t%s\t%s\t-\tno_trace\t-\t-\tnot_run\t%s\n' "$cell" "$overlay" "$nf_rc" \
      "$(grep -m1 -iE 'error' "$out/run.log" | tr '\t' ' ')" >> "$RESULTS"
    echo "  no trace; see $out/run.log"; return
  fi

  # make_graph: the index mtimes as graph_align_reads will see them, and
  # whether vg's rule predicts a failure from them.
  local pred=unknown mg
  mg=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} $c["name"]~/^make_graph/{print $c["hash"]; exit}' "$trace")
  local mgdir; mgdir=$(workdir_of "$work" "$mg")
  if [[ -n "$mgdir" && -d "$mgdir/index" ]]; then
    find "$mgdir/index" -maxdepth 1 -type f -printf '%f\t%T@\t%s\n' | sort \
      | awk -v c="$cell" 'BEGIN{OFS="\t"}{print c, $0}' >> "$MTIMES"
    pred=$(awk -F'\t' -v c="$cell" '$1==c {t[$2]=$3}
      END {d=t["index.dist"]; m=t["index.longread.withzip.min"]; z=t["index.longread.zipcodes"]
           if (d=="" || m=="" || z=="") print "unknown"
           else print ((d+0 <= m+0 && d+0 <= z+0) ? "pass" : "fail_mtime")}' "$MTIMES")
    grep -m1 'NXF_SCRATCH=' "$mgdir/.command.run" > "$out/make_graph_scratch_line.txt" 2>/dev/null || true
    grep -A8 '^nxf_unstage_outputs()' "$mgdir/.command.run" > "$out/make_graph_unstage.txt" 2>/dev/null || true
  fi

  # graph_align_reads, one row per sample
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next}
              $c["name"]~/^(make_graph|graph_align_reads)/{print $c["name"] "\t" $c["hash"] "\t" $c["status"] "\t" $c["exit"]}' "$trace" \
  | while IFS=$'\t' read -r name hash status rc; do
      local d obs msg="" p="-"
      d=$(workdir_of "$work" "$hash")
      if [[ "$status" == COMPLETED ]]; then obs=pass
      elif [[ -n "$d" ]] && grep -q 'is newer than' "$d/.command.err" 2>/dev/null; then obs=fail_mtime
      else obs=fail_other; fi
      [[ -n "$d" ]] && msg=$(grep -m1 -E 'is newer than|rror' "$d/.command.err" 2>/dev/null | tr '\t' ' ')
      [[ "$name" == graph_align_reads* ]] && p="$pred"
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$cell" "$overlay" "$nf_rc" "$name" "$status" "$rc" "$p" "$obs" "$msg" >> "$RESULTS"
      printf '  %-28s %-10s predicted=%-11s observed=%s\n' "$name" "$status" "$p" "$obs"
      if [[ -n "$d" ]]; then
        local tag; tag=$(echo "$name" | tr -c 'A-Za-z0-9_\n' '_')
        for x in .command.sh .command.err .exitcode; do
          [[ -f "$d/$x" ]] && cp "$d/$x" "$out/${tag}${x}"
        done
      fi
    done
  echo "  nextflow rc=$nf_rc"
}

SEL=("$@")
if [[ ${#SEL[@]} -eq 0 ]]; then
  SEL=(noscratch)
  for r in $(seq 1 "$REPS_B"); do SEL+=("scratch_$r"); done
  SEL+=(scratch_move)
fi
for cell in "${SEL[@]}"; do
  case "$cell" in
    noscratch)    run_cell "$cell" none ;;
    scratch_*)    [[ "$cell" == scratch_move ]] && run_cell "$cell" scratch_move.config || run_cell "$cell" scratch.config ;;
    *) echo "unknown cell $cell" >&2 ;;
  esac
done
echo
echo "Part B done: $RESULTS"
