#!/usr/bin/env bash
# What each run of the ison cell has to leave behind. Every run passed one switch
# as a string from a params file, the only way a string reaches isOn(): the
# command line turns "false" into a Boolean before the pipeline sees it.
#
#   check_ison.sh <runs dir>
#
# Off means the subset and the stages the switch gates are absent; on means they
# are there. "True" reads as on to isOn(), so it has to behave like "true".
set -uo pipefail
R=${1:?runs dir}
fail=0

ok()   { echo "  [ ok ] $1"; }
bad()  { echo "  [FAIL] $1"; fail=1; }
rc()   { cat "$R/ison_$1.rc" 2>/dev/null || echo missing; }
have() { [[ -s "$R/ison_$1/$2" || -d "$R/ison_$1/$2" ]]; }

human_off() {
  local v=$1
  [[ $(rc "$v") == 0 ]] || { bad "$v: run exited $(rc "$v")"; return; }
  if have "$v" 3_TSD_search/pangenome.trusted.vcf && ! have "$v" 3_TSD_search/pangenome.human.vcf \
     && ! have "$v" 3_TSD_search/hervk_loci.tsv; then
    ok "$v: trusted subset, no human subset, no HERV-K tables"
  else
    bad "$v: expected the trusted subset only; 3_TSD_search has: $(ls "$R/ison_$v/3_TSD_search" 2>/dev/null | tr '\n' ' ')"
  fi
}

human_on() {
  local v=$1
  [[ $(rc "$v") == 0 ]] || { bad "$v: run exited $(rc "$v")"; return; }
  if have "$v" 3_TSD_search/pangenome.human.vcf && have "$v" 3_TSD_search/human_filter_summary.txt \
     && have "$v" 3_TSD_search/hervk_loci.tsv && ! have "$v" 3_TSD_search/pangenome.trusted.vcf; then
    ok "$v: human subset and HERV-K tables, no trusted subset"
  else
    bad "$v: expected the human subset; 3_TSD_search has: $(ls "$R/ison_$v/3_TSD_search" 2>/dev/null | tr '\n' ' ')"
  fi
}

for v in human_str_false human_str_FALSE human_str_empty; do human_off "$v"; done
for v in human_str_true human_str_True; do human_on "$v"; done

v=genotype_str_false
if [[ $(rc $v) != 0 ]]; then
  bad "$v: run exited $(rc $v)"
elif ! have $v 4_Genotyping && ! have $v GraffiTE_graph && have $v 3_TSD_search/pangenome.vcf; then
  ok "$v: discovery ran, no graph and no genotyping"
else
  bad "$v: expected no 4_Genotyping and no GraffiTE_graph; the run has: $(ls "$R/ison_$v" | tr '\n' ' ')"
fi

echo "ison fail=$fail"; exit $fail
