#!/usr/bin/env bash
# Regression test for hervk_classify.py's candidate gate.
#
# is_candidate() used to admit an LTR/ERVK record only at n_hits == 1, or
# n_hits == 2 alongside SVA. That is a proxy for "this SV is HML-2 and nothing
# else", and a poor one: whether RepeatMasker calls a split terminal LTR
# LTR5_Hs or SVA_A changes n_hits without changing a base. Two of the three
# records at chr7:4.70 Mb carry n_hits = 3 for that reason alone and were
# dropped, including the one describing the common allele at the locus.
#
# The architecture table is asked first now. The hit-count rule has to stay as
# the fallback for records with no architecture row, and neither route may
# admit an SV that is mostly not HML-2.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

fail=0
chk(){ if [[ "$2" == "$3" ]]; then echo "  [ ok ] $1"; else
       echo "  [FAIL] $1: got '$2', want '$3'"; fail=1; fi; }

out=$(python3 - <<'EOF'
import importlib.util
spec = importlib.util.spec_from_file_location(
    'hervk_classify', '../../bin/hervk_classify.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
cfg = m.DEFAULTS

def info(svlen, n_hits, classes='LTR/ERVK,Retroposon/SVA'):
    return {'matching_classes': classes, 'n_hits': str(n_hits),
            '_svlen': str(svlen)}

def arch(ltr, integ):
    return {'ltr_bp': str(ltr), 'int_bp': str(integ)}

def y(b):
    return 'yes' if b else 'no'

# chr7-4706809-DEL-8503: n_hits 3, and 8503 of its 8503 bp are HML-2.
print('n3_with_arch', y(m.is_candidate(info(-8503, 3), cfg, arch(968, 7535))))
# chr7-4699715-INS-8504: the other n_hits 3 record at the same locus.
print('n3_ins_with_arch', y(m.is_candidate(info(8504, 3), cfg, arch(968, 7536))))
# Same record with no architecture row: the hit-count fallback still refuses.
print('n3_no_arch', y(m.is_candidate(info(-8503, 3), cfg, None)))
# The two shapes the old rule allowed must keep working without an arch row.
print('n1_no_arch', y(m.is_candidate(info(8504, 1, 'LTR/ERVK'), cfg, None)))
print('n2_sva_no_arch', y(m.is_candidate(info(8504, 2), cfg, None)))
# An SV that is only 6% HML-2 is not a candidate however few hits it has.
print('mostly_not_hml2', y(m.is_candidate(info(8000, 3), cfg, arch(200, 300))))
print('mostly_not_hml2_n1',
      y(m.is_candidate(info(8000, 1, 'LTR/ERVK'), cfg, arch(200, 300))))
# The size cap is checked before anything else.
print('over_cap', y(m.is_candidate(info(4000000, 1), cfg, arch(2000000, 2000000))))
# A record with no LTR/ERVK at all never reaches the HML-2 tests.
print('not_ervk', y(m.is_candidate(info(8504, 1, 'LINE/L1'), cfg, arch(968, 7536))))
EOF
)
g(){ echo "$out" | awk -v k="$1" '$1==k{print $2}'; }

chk "n_hits=3 admitted when the architecture is HML-2"  "$(g n3_with_arch)"      "yes"
chk "the other n_hits=3 chr7 record is admitted too"    "$(g n3_ins_with_arch)"  "yes"
chk "n_hits=3 with no architecture row still refused"   "$(g n3_no_arch)"        "no"
chk "n_hits=1 keeps working with no architecture row"   "$(g n1_no_arch)"        "yes"
chk "n_hits=2 with SVA keeps working"                   "$(g n2_sva_no_arch)"    "yes"
chk "mostly non-HML-2 refused at n_hits=3"              "$(g mostly_not_hml2)"   "no"
chk "mostly non-HML-2 refused at n_hits=1 as before"    "$(g mostly_not_hml2_n1)" "yes"
chk "over the SVLEN cap is refused"                     "$(g over_cap)"          "no"
chk "a non-ERVK record is never a candidate"            "$(g not_ervk)"          "no"

[[ $fail -eq 0 ]] && echo "PASS" || { echo "FAIL"; exit 1; }
