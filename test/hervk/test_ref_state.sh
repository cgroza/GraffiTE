#!/usr/bin/env bash
# Regression test for hervk_ref_state.py's array handling.
#
# Two things are pinned down here.
#
# unit_structure() counts proviral units in a clustered reference element and
# measures their period. Every geometry below is real, taken from chm13v2.0 and
# checked by aligning the chr6 provirus back onto each window, so the expected
# periods are measurements rather than round numbers: chr6 8465, chr7 8504,
# chr12 4935. The chr12 case matters because RepeatMasker splits its degraded
# internal region into four fragments, and that is still one unit.
#
# truncated_by_window() used to exempt state == 'provirus' from the rescue
# pass. That hid the case the rescue exists for: a tandem array reads
# 'provirus' off its first unit and runs into the window edge, so the rest of
# the array is never masked. chr7 did that at flank=12000, reporting an element
# that ended at POS+12000 to the base.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

fail=0
chk(){ if [[ "$2" == "$3" ]]; then echo "  [ ok ] $1"; else
       echo "  [FAIL] $1: got '$2', want '$3'"; fail=1; fi; }

out=$(python3 - <<'EOF'
import importlib.util
spec = importlib.util.spec_from_file_location(
    'hervk_ref_state', '../../bin/hervk_ref_state.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

def L(qs, qe):
    return {'name': 'LTR5_Hs', 'qstart': qs, 'qend': qe, 'bp': qe - qs + 1}

def I(qs, qe):
    return {'name': 'HERVK-int', 'qstart': qs, 'qend': qe, 'bp': qe - qs + 1}

def show(tag, el):
    n, p = m.unit_structure(el)
    print(tag, 'None' if n is None else n, 'None' if p is None else p)

# chr6:78,894,317-78,903,741 -- one complete provirus, period 8465.
show('chr6', [L(78894317, 78895276), I(78895277, 78902781), L(78902782, 78903741)])

# chr7:4,699,540-4,717,514 (7p22.1a) -- two proviruses sharing the middle LTR,
# period 8504.
show('chr7', [L(4699540, 4700507), I(4700508, 4708043), L(4708044, 4709011),
              I(4709012, 4716546), L(4716547, 4717514)])

# chr12:133,147,852-133,153,798 -- one unit whose internal region survives as
# four fragments. Grouping has to collapse them, or this reads as four units.
show('chr12', [L(133147852, 133148864),
               I(133148865, 133149571), I(133149572, 133150255),
               I(133150256, 133150531), I(133150532, 133152481),
               L(133152787, 133153798)])

# A solo LTR is not an array.
show('solo', [L(1000, 1968)])

# chr6-78894317-DEL-8465 reads INT/LTR: a truncated view, not a unit.
show('truncated', [I(1000, 8536), L(8537, 9504)])

# --- truncated_by_window: the edge decides, not the state ---
cfg = dict(m.DEFAULTS)
fp = ('sv', 'chr7', 4700333, 4700333, 8504)          # insertion: a point
# element runs to exactly POS+12000, i.e. hard against the window edge
res = {'state': 'provirus', 'elem_start': 4699540, 'elem_end': 4700333 + 12000}
print('edge_provirus', m.truncated_by_window(res, fp, cfg, 12000))
# same element, comfortably inside a wider window
print('inside_provirus', m.truncated_by_window(res, fp, cfg, 24000))
EOF
)
g(){ echo "$out" | awk -v k="$1" '$1==k{print $2, $3}'; }
g1(){ echo "$out" | awk -v k="$1" '$1==k{print $2}'; }

chk "chr6 provirus is 1 unit, period 8465"      "$(g chr6)"       "1 8465"
chk "chr7 array is 2 units, period 8504"        "$(g chr7)"       "2 8504"
chk "chr12 fragmented INT is still 1 unit"      "$(g chr12)"      "1 4935"
chk "solo LTR is not an array"                  "$(g solo)"       "None None"
chk "truncated INT/LTR read is not an array"    "$(g truncated)"  "None None"
chk "provirus at the window edge is rescued"    "$(g1 edge_provirus)"   "True"
chk "provirus inside its window is left alone"  "$(g1 inside_provirus)" "False"

[[ $fail -eq 0 ]] && echo "PASS" || { echo "FAIL"; exit 1; }
