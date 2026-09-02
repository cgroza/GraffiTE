#!/usr/bin/env bash
# Regression test for hervk_reconcile.py's build_locus_record().
#
# Three shapes reach it. (a) and (b) worked already; (c) returned an error and
# the locus was skipped, so a multi-allelic locus whose members are insertions
# with no deletion spanning them could never consolidate. chr7:4,699,714 is that
# shape -- two insertions at 4,699,714 and 4,700,333 and a deletion starting at
# 4,706,808, downstream of both -- which is why it was the one CaG locus that
# never merged.
#
# The geometry here is the chr7 one in miniature, against a synthetic contig so
# the test needs no reference genome. Case (b) is re-checked alongside, because
# a member's own REF field must keep winning over anything fetched.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

command -v samtools >/dev/null || { echo "  [skip] samtools not on PATH"; exit 0; }

fail=0
chk(){ if [[ "$2" == "$3" ]]; then echo "  [ ok ] $1"; else
       echo "  [FAIL] $1: got '$2', want '$3'"; fail=1; fi; }

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

python3 - "$tmp" <<'EOF'
import sys, os, random
tmp = sys.argv[1]
random.seed(7)
seq = ''.join(random.choice('ACGT') for _ in range(200))
with open(os.path.join(tmp, 'ref.fa'), 'w') as fh:
    fh.write('>t1\n')
    for i in range(0, len(seq), 60):
        fh.write(seq[i:i+60] + '\n')
open(os.path.join(tmp, 'seq.txt'), 'w').write(seq)
EOF
samtools faidx "$tmp/ref.fa"

out=$(python3 - "$tmp" <<'EOF'
import importlib.util, sys, os
tmp = sys.argv[1]
spec = importlib.util.spec_from_file_location(
    'hervk_reconcile', '../../bin/hervk_reconcile.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
seq = open(os.path.join(tmp, 'seq.txt')).read().strip()
ref = os.path.join(tmp, 'ref.fa')

def rec(pos, r, a):                       # a VCF record as the code sees it
    return ['t1', str(pos), f'v{pos}', r, a, '.', 'PASS', '.', 'GT']

def at(pos):                              # 1-based base at pos
    return seq[pos - 1]

INS_A = 'A' * 20
INS_B = 'C' * 20

# (c) two insertions and a deletion starting downstream of both: the chr7 shape
a = rec(10, at(10), at(10) + INS_A)
b = rec(30, at(30), at(30) + INS_B)
d = rec(100, seq[99:149], at(100))        # REF 50 bp, ALT 1 bp -> removes 49
mem = [a, b, d]
keep = [0, 1, 2]
idx  = {0: 1, 1: 2, 2: 3}
states = ['prov_x3', 'prov_x2', 'provirus']

chrom, pos, r, alt, err = m.build_locus_record(mem, keep, states, idx, ref)
print('c_err', err or 'none')
if not err:
    alts = alt.split(',')
    span = seq[9:149]                     # 10..149, 140 bp
    print('c_pos', pos)
    print('c_reflen', len(r))
    print('c_ref_exact', 'yes' if r == span else 'no')
    print('c_altlens', ','.join(str(len(x)) for x in alts))
    print('c_altA_exact', 'yes' if alts[0] == span[:1] + INS_A + span[1:] else 'no')
    print('c_altB_exact', 'yes' if alts[1] == span[:21] + INS_B + span[21:] else 'no')
    print('c_altD_exact', 'yes' if alts[2] == span[:91] else 'no')

# (c) without a reference: must refuse rather than invent a span
_, _, _, _, err2 = m.build_locus_record(mem, keep, states, idx, None)
print('c_norefs', 'refused' if err2 else 'built')

# (b) a deletion that DOES span the others: its own REF must be used, and the
# result must not depend on the FASTA being available
big = rec(5, seq[4:160], at(5))           # spans 5..160, covers both insertions
mem2 = [a, b, big]
chrom2, pos2, r2, alt2, err3 = m.build_locus_record(
    mem2, [0, 1, 2], states, idx, ref)
print('b_err', err3 or 'none')
print('b_pos', pos2)
print('b_ref_is_member', 'yes' if r2 == seq[4:160] else 'no')
_, pos2b, r2b, alt2b, _ = m.build_locus_record(mem2, [0, 1, 2], states, idx, None)
print('b_noref_same', 'yes' if (pos2b, r2b, alt2b) == (pos2, r2, alt2) else 'no')

# (a) one allele state: the leftmost member is emitted unchanged
_, posa, ra, alta, erra = m.build_locus_record(
    [a, b], [0, 1], ['prov_x2'], {0: 1, 1: 1}, ref)
print('a_pos', posa); print('a_ref', ra); print('a_alt_is_A', 'yes' if alta == a[4] else 'no')
EOF
)
g(){ echo "$out" | awk -v k="$1" '$1==k{$1=""; sub(/^ /,""); print}'; }

chk "case c builds without a spanning deletion" "$(g c_err)"        "none"
chk "case c anchors at the leftmost member"     "$(g c_pos)"        "10"
chk "case c REF spans every member"             "$(g c_reflen)"     "140"
chk "case c REF is the reference span verbatim" "$(g c_ref_exact)"  "yes"
chk "case c ALT lengths are +20, +20, -49"      "$(g c_altlens)"    "160,160,91"
chk "case c splices insertion A exactly"        "$(g c_altA_exact)" "yes"
chk "case c splices insertion B exactly"        "$(g c_altB_exact)" "yes"
chk "case c splices the deletion exactly"       "$(g c_altD_exact)" "yes"
chk "case c refuses when no reference is given" "$(g c_norefs)"     "refused"
chk "case b still builds"                       "$(g b_err)"        "none"
chk "case b anchors on the spanning deletion"   "$(g b_pos)"        "5"
chk "case b uses the member's own REF"          "$(g b_ref_is_member)" "yes"
chk "case b ignores the reference either way"   "$(g b_noref_same)" "yes"
chk "case a emits the leftmost member"          "$(g a_pos)"        "10"
chk "case a keeps that member's ALT"            "$(g a_alt_is_A)"   "yes"

[[ $fail -eq 0 ]] && echo "PASS" || { echo "FAIL"; exit 1; }
