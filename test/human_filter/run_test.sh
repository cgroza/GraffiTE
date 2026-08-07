#!/usr/bin/env bash
# Regression test for the --human pME filter (module/main.nf, process
# concat_repeatmask). Builds the same bcftools expression from the defaults in
# nextflow.config and checks which fixture records survive.
#
# The expression assembly below mirrors the Groovy in module/main.nf and must be
# kept in sync with it.
#
# Usage: bash test/human_filter/run_test.sh
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VCF="${HERE}/human_filter_fixture.vcf"
CONFIG="${HERE}/../../nextflow.config"

command -v bcftools >/dev/null || { echo "bcftools not found"; exit 1; }

# --- defaults from nextflow.config -------------------------------------------
cfg() { sed -n "s/^[[:space:]]*$1[[:space:]]*=[[:space:]]*\"\{0,1\}\([^\"\/]*\)\"\{0,1\}.*/\1/p" "$CONFIG" | head -1 | sed 's/[[:space:]]*$//'; }
ALU=$(cfg human_alu_ids);        L1=$(cfg human_l1_ids)
SVA=$(cfg human_sva_ids);        HERVK=$(cfg human_hervk_ids)
MINSV=$(cfg human_min_svlen);    MAXTR=$(cfg human_max_ultra_span)
IGNORE=$(cfg human_ignore_filter); PAIR=$(cfg hervk_sva_pair)
PAIRMAX=$(cfg hervk_pair_max_svlen)

# --- expression assembly (mirrors module/main.nf) ----------------------------
or_ids() { # comma-separated regex list -> OR'd repeat_ids clauses
  local out="" re
  IFS=',' read -ra RE <<< "$1"
  for re in "${RE[@]}"; do
    re="$(echo "$re" | tr -d '[:space:]')"
    [ -n "$out" ] && out+=" | "
    out+="repeat_ids~\"${re}\""
  done
  echo "(${out})"
}
grp() { # class, regex list
  if [ -n "$2" ]; then echo "(matching_classes=\"$1\" & $(or_ids "$2"))"
  else echo "matching_classes=\"$1\""; fi
}

build_filter() {
  local ids size single pair hits base
  ids="$(grp 'SINE/Alu' "$ALU") | $(grp 'LINE/L1' "$L1") | $(grp 'Retroposon/SVA' "$SVA") | $(grp 'Simple_repeat' "$SVA") | $(grp 'LTR/ERVK' "$HERVK")"
  size="abs(SVLEN)>=${MINSV} & (ULTRA_TR_span<${MAXTR} | matching_classes=\"Simple_repeat\")"
  single="n_hits==1 & (matching_classes=\"LTR/ERVK\" | matching_classes=\"Simple_repeat\" | polyA=\"TRUE\")"
  pair="n_hits==2 & matching_classes=\"LTR/ERVK\" & matching_classes=\"Retroposon/SVA\" & repeat_ids~\"^HERVK-int\" & abs(SVLEN)<=${PAIRMAX}"
  if [ "$PAIR" = "true" ]; then hits="((${single}) | (${pair}))"; else hits="(${single})"; fi
  base="(${ids}) & ${size} & ${hits}"
  if [ "$IGNORE" = "true" ]; then echo "$base"; else echo "(${base}) & FILTER=\"PASS\""; fi
}

kept() { bcftools query -i "$(build_filter)" -f '%ID\n' "$VCF" 2>/dev/null | grep -v '^pass=' | sort | tr '\n' ' ' | sed 's/ $//'; }

fail=0
check() { # name, expected
  local got; got="$(kept)"
  if [ "$got" = "$2" ]; then
    printf 'ok   %s\n' "$1"
  else
    printf 'FAIL %s\n       expected: %s\n       got     : %s\n' "$1" "$2" "$got"; fail=1
  fi
}

# --- cases -------------------------------------------------------------------
check "defaults" \
  "AluY HERVKint L1HS L1HSdel LTR5Hs PAIRfwd PAIRrev SVAD SVAF"

L1="^L1HS,^L1PA2"
check "human_l1_ids adds L1PA2" \
  "AluY HERVKint L1HS L1HSdel L1PA2 LTR5Hs PAIRfwd PAIRrev SVAD SVAF"
L1="$(cfg human_l1_ids)"

MINSV=100
check "human_min_svlen=100 admits the 140 bp AluY" \
  "AluY AluYshort HERVKint L1HS L1HSdel LTR5Hs PAIRfwd PAIRrev SVAD SVAF"
MINSV=$(cfg human_min_svlen)

PAIR=false
check "hervk_sva_pair=false drops the HERVK+SVA pairs" \
  "AluY HERVKint L1HS L1HSdel LTR5Hs SVAD SVAF"
PAIR=$(cfg hervk_sva_pair)

IGNORE=true
check "human_ignore_filter=true admits the non-PASS record" \
  "AluY AluYnopass HERVKint L1HS L1HSdel LTR5Hs PAIRfwd PAIRrev SVAD SVAF"
IGNORE=$(cfg human_ignore_filter)

HERVK="^HERVK-int,^LTR5"
check "human_hervk_ids=^LTR5 admits the ancestral LTR5" \
  "AluY HERVKint L1HS L1HSdel LTR5 LTR5Hs PAIRfwd PAIRrev SVAD SVAF"
HERVK=$(cfg human_hervk_ids)

SVA=""
check "empty human_sva_ids keeps the whole SVA/Simple_repeat classes" \
  "AluY HERVKint L1HS L1HSdel LTR5Hs PAIRfwd PAIRrev SVAA SVAD SVAF"
SVA=$(cfg human_sva_ids)

echo
if [ "$fail" -eq 0 ]; then echo "all human filter tests passed"; else echo "human filter tests FAILED"; fi
exit "$fail"
