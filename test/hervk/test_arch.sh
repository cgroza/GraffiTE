#!/usr/bin/env bash
# Regression test for hervk_arch.py's SINE-R reassignment.
#
# reassign_sine_r() decides whether an SVA hit beside an HML-2 element is really
# the LTR-derived SINE-R domain, and it decides that from the SVA consensus
# START coordinate. Hits that come from a BED4 annotation -- the
# --hervk_ref_annotation path in hervk_ref_state.py -- have no consensus
# coordinates at all, so that field is None and the comparison used to raise
# TypeError. The fix must skip those hits WITHOUT disabling the reassignment for
# real RepeatMasker hits, which is what the second case here pins down.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

fail=0
chk(){ if [[ "$2" == "$3" ]]; then echo "  [ ok ] $1"; else
       echo "  [FAIL] $1: got '$2', want '$3'"; fail=1; fi; }

out=$(python3 - <<'EOF'
import importlib.util, sys
spec = importlib.util.spec_from_file_location('hervk_arch', '../../bin/hervk_arch.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

def frag(name, klass, cs, ce, qs, qe, sw):
    return {'name': name, 'klass': klass, 'cons_start': cs, 'cons_end': ce,
            'cons_len': (ce or 0), 'qstart': qs, 'qend': qe, 'bp': qe - qs + 1,
            'sw': sw, 'strand': '+', 'link': name}

INT = lambda: frag('HERVK-int', 'LTR/ERVK', 1, 7536, 310, 7845, 5000)

# All three cases share one geometry -- the SVA hit abuts the internal region
# with a 9 bp gap, well inside sine_r_max_gap -- so consensus START is the only
# variable between them.
# 1. BED4-derived: no consensus coordinates anywhere.
try:
    r = m.reassign_sine_r([frag('SVA_A', 'Retroposon/SVA', None, None, 1, 300, 1000),
                           INT()], m.DEFAULTS)
    print('bed_ok', 'yes')
    print('bed_reassigned', 'yes' if not any(f['name'].startswith('SVA') for f in r) else 'no')
except TypeError:
    print('bed_ok', 'no'); print('bed_reassigned', 'crash')

# 2. A real RepeatMasker hit in the SINE-R range still gets reassigned --
#    the None guard must not silently switch the whole mechanism off.
r = m.reassign_sine_r([frag('SVA_A', 'Retroposon/SVA', 951, 1113, 1, 300, 1000),
                       INT()], m.DEFAULTS)
print('rm_reassigned', 'yes' if not any(f['name'].startswith('SVA') for f in r) else 'no')

# 3. A real hit BELOW the SINE-R threshold is left alone.
r = m.reassign_sine_r([frag('SVA_A', 'Retroposon/SVA', 100, 400, 1, 300, 1000),
                       INT()], m.DEFAULTS)
print('rm_low_reassigned', 'yes' if not any(f['name'].startswith('SVA') for f in r) else 'no')
EOF
)
g(){ echo "$out" | awk -v k="$1" '$1==k{print $2}'; }

chk "BED4 hit without consensus coords does not crash" "$(g bed_ok)"            "yes"
chk "BED4 hit is left as SVA (cannot identify SINE-R)" "$(g bed_reassigned)"    "no"
chk "real SINE-R hit is still reassigned to LTR"       "$(g rm_reassigned)"     "yes"
chk "hit below the SINE-R threshold is left alone"     "$(g rm_low_reassigned)" "no"

[[ $fail -eq 0 ]] && echo "PASS" || { echo "FAIL"; exit 1; }
