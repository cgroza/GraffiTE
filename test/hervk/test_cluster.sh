#!/usr/bin/env bash
# Regression test for hervk_reconcile.py's locus clustering.
#
# cluster() joined two records when their reference-element tuples matched to
# the base, or when their footprints were within --window. Both tests fail at a
# copy-number array. hervk_ref_state masks a window per record, so the element
# span it reports depends on that window: the two chr7 insertions came back
# 4699540-4712333 and 4699540-4711714 for one array and never compared equal.
# And the deletion at the same locus sits 6475 bp from the nearer insertion,
# well outside the 1200 bp default window but well inside the 17,975 bp array.
#
# Overlap decides now, and a record inside a known element joins it.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

fail=0
chk(){ if [[ "$2" == "$3" ]]; then echo "  [ ok ] $1"; else
       echo "  [FAIL] $1: got '$2', want '$3'"; fail=1; fi; }

out=$(python3 - <<'EOF'
import importlib.util
spec = importlib.util.spec_from_file_location(
    'hervk_reconcile', '../../bin/hervk_reconcile.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

def rec(vid, chrom, start, end):
    return {'id': vid, 'chrom': chrom, 'start': start, 'end': end,
            'ref_allele': 'prov_x2', 'alt_allele': 'prov_x3',
            'cls': 'copy_number', 'evidence': 'CNV_PERIOD', 'k': '',
            'arch': ''}

def ref(chrom, start, end):
    return {'ref_elem_chrom': chrom, 'ref_elem_start': str(start),
            'ref_elem_end': str(end)}

# The real chr7 locus: two insertions 619 bp apart and a deletion 6475 bp on,
# with the element span reported differently for each record because each was
# masked in its own window.
recs = [rec('ins_a', 'chr7', 4699714, 4699714),
        rec('ins_b', 'chr7', 4700333, 4700333),
        rec('del_c', 'chr7', 4706808, 4715311)]
ref_tbl = {'ins_a': ref('chr7', 4699540, 4711714),
           'ins_b': ref('chr7', 4699540, 4712333),
           'del_c': ref('chr7', 4699540, 4717514)}
_, table = m.cluster(recs, ref_tbl, 1200)
print('chr7_loci', len(table))
print('chr7_members', table[0]['n_records'] if table else 0)

# Control: two records far apart with no element information must stay apart.
far = [rec('x', 'chr9', 1000, 1000), rec('y', 'chr9', 90000, 90000)]
_, t2 = m.cluster(far, {}, 1200)
print('far_loci', len(t2))

# Control: a record beyond the array does not get pulled in.
recs3 = recs + [rec('outside', 'chr7', 4760000, 4760000)]
_, t3 = m.cluster(recs3, ref_tbl, 1200)
print('outside_loci', len(t3))

# Records at one element whose spans differ still merge, which is the case
# exact-tuple matching missed.
two = [rec('p', 'chr6', 78894316, 78902781), rec('q', 'chr6', 78894875, 78894875)]
ref2 = {'p': ref('chr6', 78893849, 78903320), 'q': ref('chr6', 78894317, 78903741)}
_, t4 = m.cluster(two, ref2, 1200)
print('chr6_loci', len(t4))
EOF
)
g(){ echo "$out" | awk -v k="$1" '$1==k{print $2}'; }

chk "the three chr7 records form one locus"        "$(g chr7_loci)"    "1"
chk "that locus holds all three records"           "$(g chr7_members)" "3"
chk "distant records with no element stay apart"   "$(g far_loci)"     "2"
chk "a record beyond the array is not pulled in"   "$(g outside_loci)" "2"
chk "same element, different spans, still merges"  "$(g chr6_loci)"    "1"

[[ $fail -eq 0 ]] && echo "PASS" || { echo "FAIL"; exit 1; }
