#!/usr/bin/env bash
# Tally Parts A and B and apply the outcome labels fixed in README_HPC_CLAUDE.md
# before the run. Writes $WORKDIR/VERDICT.txt. Reports; does not interpret.
set -uo pipefail
export LC_ALL=C
# shellcheck disable=SC1091
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

A="$WORKDIR/partA/results.tsv"
B="$WORKDIR/partB/results.tsv"
V="$WORKDIR/VERDICT.txt"

{
echo "# issue #127 mtime test -- $(date -u +%Y-%m-%dT%H:%M:%SZ) on $(hostname)"
echo

if [[ -f "$A" ]]; then
  echo "## Part A: observed per source and variant"
  printf "%-24s %-12s %5s %5s %10s %10s %10s\n" "source (fs)" variant n pass fail_mtime fail_other =predicted
  awk -F'\t' 'NR>1 {k=$1" ("$3")\t"$6; n[k]++; o[k,$8]++; if($7==$8) agree[k]++; keys[k]=1}
    END {for (k in keys) {split(k, a, "\t")
           printf "%-24s %-12s %5d %5d %10d %10d %10d\n", a[1], a[2], n[k], o[k,"pass"], o[k,"fail_mtime"], o[k,"fail_other"], agree[k]}}' "$A" | sort
  echo
else
  echo "## Part A: no results ($A)"; echo
fi

if [[ -f "$B" ]]; then
  echo "## Part B: per task"
  awk -F'\t' 'NR>1 {printf "%-14s %-28s %-10s predicted=%-11s observed=%s\n", $1, $4, $5, $7, $8}' "$B"
  echo
else
  echo "## Part B: no results ($B)"; echo
fi

echo "## Labels"
# Validity: the controls have to behave, or nothing else in the run counts.
inv=""
if [[ -f "$A" ]]; then
  bad=$(awk -F'\t' 'NR>1 && $6=="forced_bad" && $8!="fail_mtime"' "$A" | wc -l | tr -d " ")
  [[ $bad -gt 0 ]] && inv+=" forced_bad_did_not_fail($bad)"
  bad=$(awk -F'\t' 'NR>1 && ($6=="forced_good" || $6=="control") && $8!="pass"' "$A" | wc -l | tr -d " ")
  [[ $bad -gt 0 ]] && inv+=" control_or_forced_good_failed($bad)"
  [[ $(awk -F'\t' 'NR>1 && $6=="forced_bad"' "$A" | wc -l | tr -d " ") -eq 0 ]] && inv+=" no_forced_bad_rows"
fi
if [[ -f "$B" ]]; then
  bad=$(awk -F'\t' 'NR>1 && $1=="noscratch" && $8!="pass"' "$B" | wc -l | tr -d " ")
  [[ $bad -gt 0 ]] && inv+=" noscratch_failed($bad)"
fi
echo "VALIDITY: $([[ -z "$inv" ]] && echo OK || echo "INVALID:$inv")"

# Reproduction: the #127 error from a plain copy (A cp) or from the pipeline
# under the cluster scratch setting (B scratch_<n>).
a_nf=0; a_extra=0; b_fail=0
[[ -f "$A" ]] && a_nf=$(awk -F'\t' 'NR>1 && $1=="nf_scratch" && $6=="cp" && $8=="fail_mtime"' "$A" | wc -l | tr -d " ")
[[ -f "$A" ]] && a_extra=$(awk -F'\t' 'NR>1 && $1!="nf_scratch" && $6=="cp" && $8=="fail_mtime"' "$A" | wc -l | tr -d " ")
[[ -f "$B" ]] && b_fail=$(awk -F'\t' 'NR>1 && $1 ~ /^scratch_[0-9]+$/ && $4 ~ /^graph_align_reads/ && $8=="fail_mtime"' "$B" | wc -l | tr -d " ")
if   [[ $a_nf -gt 0 || $b_fail -gt 0 ]]; then repro=REPRODUCED
elif [[ $a_extra -gt 0 ]];               then repro=REPRODUCED_OTHER_FS_ONLY
else                                          repro=NOT_REPRODUCED; fi
echo "REPRODUCTION: $repro (A cp nf_scratch fail_mtime=$a_nf, A cp other fs fail_mtime=$a_extra, B scratch fail_mtime=$b_fail)"

# Mechanism: did the inode-order prediction match the cp outcome?
if [[ -f "$A" ]]; then
  read -r n_cp n_agree < <(awk -F'\t' 'NR>1 && $6=="cp" {n++; if($7==$8) a++} END{print n+0, a+0}' "$A")
  echo "INODE_ORDER_PREDICTS_CP: $n_agree of $n_cp"
fi

# Candidate fix: mv in A, stageOutMode 'move' in B.
mv_bad=0; mv_n=0
[[ -f "$A" ]] && mv_bad=$(awk -F'\t' 'NR>1 && $6=="mv" && $8!="pass"' "$A" | wc -l | tr -d " ")
[[ -f "$A" ]] && mv_n=$(awk -F'\t' 'NR>1 && $6=="mv"' "$A" | wc -l | tr -d " ")
bm_bad=0; bm_n=0
[[ -f "$B" ]] && bm_bad=$(awk -F'\t' 'NR>1 && $1=="scratch_move" && $4 ~ /^graph_align_reads/ && $8!="pass"' "$B" | wc -l | tr -d " ")
[[ -f "$B" ]] && bm_n=$(awk -F'\t' 'NR>1 && $1=="scratch_move" && $4 ~ /^graph_align_reads/' "$B" | wc -l | tr -d " ")
if [[ $((mv_n + bm_n)) -eq 0 ]]; then mv_label=NOT_RUN
elif [[ $((mv_bad + bm_bad)) -eq 0 ]]; then mv_label=MOVE_OK
else mv_label=MOVE_FAILED; fi
echo "CANDIDATE_FIX: $mv_label (A mv: $((mv_n - mv_bad))/$mv_n pass, B scratch_move: $((bm_n - bm_bad))/$bm_n pass)"
} | tee "$V"
