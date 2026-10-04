#!/usr/bin/env python3
"""What the hervkref cell has to show.

Each run gave hervk_annotate a reference annotation through
--hervk_ref_annotation instead of letting it mask windows of the reference:
a RepeatMasker .out of the whole reference, the same hits as a BED, and a BED
with no hit at all. The first two hold what the in-pipeline masking finds, so
their reference states must match the spine's. The empty one holds nothing, so
every locus the spine found in the reference must come back null; that is what
shows the annotation was read and the masking skipped.

usage: check_hervkref.py <runs dir> <spine run dir> <work dir>
"""
import csv
import glob
import os
import re
import sys


def states(run):
    path = os.path.join(run, '3_TSD_search', 'hervk_refstate.tsv')
    with open(path) as fh:
        return {r['id']: r for r in csv.DictReader(fh, delimiter='\t')}


def annotation_seen(run, work):
    """The value hervk_annotate's script tested to choose its branch. The script
    holds both branches, so their text proves nothing; the rendered test does:
    `if [[ -n "<annotation>" ]]` is empty when the reference was masked."""
    trace = os.path.join(run, 'nextflow_trace.txt')
    for r in csv.DictReader(open(trace), delimiter='\t'):
        if r['name'].startswith('hervk_annotate'):
            for cmd in glob.glob(os.path.join(work, r['hash'] + '*', '.command.sh')):
                m = re.search(r'if \[\[ -n "([^"]*)" \]\]; then\n\s*hervk_ref_state\.py', open(cmd).read())
                return m.group(1) if m else None
    return None


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    runs, spine, work = sys.argv[1:]
    fail = 0

    def check(label, good, detail=''):
        nonlocal fail
        print(f"  [{' ok ' if good else 'FAIL'}] {label}" + (f": {detail}" if detail and not good else ''))
        fail += not good

    s = states(spine)
    found = sorted(k for k, r in s.items() if r['ref_state'] not in ('null', 'unknown'))
    check("the spine finds HERV-K sequence in the reference", bool(found), "every candidate is null")

    for mode in ('out', 'bed', 'empty'):
        run = os.path.join(runs, f'hervkref_{mode}')
        seen = annotation_seen(run, work)
        want = os.path.join(runs, 'hervkref_inputs', {'out': 'ref.out', 'bed': 'ref.bed', 'empty': 'ref.empty.bed'}[mode])
        check(f"{mode}: hervk_annotate took the annotation branch, not the masking one",
              seen == want, f"its branch test read {seen!r}")
        r = states(run)
        if mode in ('out', 'bed'):
            # A BED holds no consensus coordinates, so hervk_arch.py writes the
            # architecture as a length (LTR:~968bp, not LTR:1-968); every other
            # column has to match.
            skip = {'ref_arch'} if mode == 'bed' else set()
            cut = lambda d: {c: v for c, v in d.items() if c not in skip} if d else d
            diff = sorted(k for k in s.keys() | r.keys() if cut(s.get(k)) != cut(r.get(k)))
            what = "reference state" + (", ref_arch aside," if skip else "")
            check(f"{mode}: the same {what} as the spine on all {len(s)} candidates", not diff,
                  f"{len(diff)} differ, first {diff[0]}: {s.get(diff[0])} vs {r.get(diff[0])}" if diff else '')
            for k in sorted(found):
                if skip and s[k]['ref_arch'] != r.get(k, {}).get('ref_arch'):
                    print(f"  [info] {k}: ref_arch {s[k]['ref_arch']} from the .out, {r[k]['ref_arch']} from the BED")
        else:
            still = [k for k in found if r.get(k, {}).get('ref_state') != 'null']
            check(f"empty: the {len(found)} reference loci come back null", not still,
                  f"{still} still found")
            others = sorted(k for k in s if k not in found and s[k] != r.get(k))
            check(f"empty: the other {len(s) - len(found)} candidates unchanged", not others,
                  f"{len(others)} differ, first {others[0]}" if others else '')

    print(f"hervkref fail={fail}")
    sys.exit(1 if fail else 0)


if __name__ == '__main__':
    main()
