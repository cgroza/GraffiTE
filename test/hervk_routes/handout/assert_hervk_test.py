#!/usr/bin/env python3
"""
Check a GraffiTE --human discovery run against the HERV-K v2 predictions.

Every expectation here was derived from the raw RepeatMasker tables of the
existing CaG run, so this is a genuine test: the architecture was measured
first and the classifier is being asked to reproduce it. A failure means the
new code disagrees with the RepeatMasker evidence, not that a number moved.

Exits non-zero on any failure.
"""
import argparse
import os
import sys


def load_tsv(path, key):
    rows = {}
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        ki = header.index(key)
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            f += [''] * (len(header) - len(f))
            rows[f[ki]] = dict(zip(header, f))
    return rows, header


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--outdir', required=True, help='GraffiTE --out directory')
    ap.add_argument('--expected', default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'EXPECTED.tsv'))
    args = ap.parse_args()

    d = os.path.join(args.outdir, '3_TSD_search')
    calls_p = os.path.join(d, 'hervk_calls.tsv')
    loci_p = os.path.join(d, 'hervk_loci.tsv')
    for p in (calls_p, loci_p):
        if not os.path.exists(p):
            sys.exit(f'FATAL: missing {p}\n'
                     '  The hervk_annotate process did not run or did not publish.\n'
                     '  Check that --human was set and that the run reached 3_TSD_search.')

    calls, _ = load_tsv(calls_p, 'id')
    loci, _ = load_tsv(loci_p, 'locus_id')
    expected, _ = load_tsv(args.expected, 'id')

    failures, warnings, checked = [], [], 0

    # ---- 1. per-record expectations -------------------------------------
    for vid, want in expected.items():
        got = calls.get(vid)
        if got is None:
            failures.append(f'{vid}: absent from hervk_calls.tsv')
            continue
        checked += 1
        if want['expect_evidence'] and got['evidence'] != want['expect_evidence']:
            failures.append(
                f'{vid}: evidence {got["evidence"]}, expected '
                f'{want["expect_evidence"]} ({want["note"]})')
        wc = want['expect_class']
        if wc and wc != '*_prov' and got['class'] != wc:
            failures.append(f'{vid}: class {got["class"]}, expected {wc}')
        if wc == '*_prov' and not got['class'].endswith('_prov'):
            failures.append(f'{vid}: class {got["class"]}, expected a proviral class')
        if want['expect_k'] and got['k'] != want['expect_k']:
            failures.append(f'{vid}: k={got["k"]}, expected {want["expect_k"]}')

    # ---- 2. the headline results ----------------------------------------
    null_prov = sorted(v for v, r in calls.items() if r['class'] == 'null_prov')
    if not null_prov:
        failures.append(
            'null_prov count is 0. This is the whole point of the rewrite: a '
            'complete provirus into an empty site was unreachable in v1 '
            'because it needed ~1936 bp of LTR that annotate_vcf.R collapses '
            'away. Zero here means the architecture layer is not reaching the '
            'classifier -- check hervk_arch.tsv is populated.')

    merge_loci = {r['locus_id']: r for r in loci.values()
                  if 'MERGE_CANDIDATE' in r['flags']}
    for want_locus, why in (
            ('HERVK_chr11_101704640', 'the pair snarl-based grouping misses '
                                      '(574 bp apart, k=574)'),
            ('HERVK_chr12_55299985', 'the three-allele locus'),
            ('HERVK_chr6_78894316', 'the polarity contradiction')):
        if want_locus not in merge_loci:
            failures.append(f'{want_locus} not flagged MERGE_CANDIDATE — {why}')

    # ---- 3. things worth seeing but not failing on -----------------------
    unknown_ref = [v for v, r in calls.items()
                   if r['evidence'] == 'REF_ANNOT' and r['ref_state'] == 'unknown']
    if unknown_ref:
        warnings.append(
            f'{len(unknown_ref)} records reached REF_ANNOT with '
            'ref_state=unknown — the reference masking produced nothing for '
            'them. Check hervk_refstate.tsv and that the TE library contains '
            'LTR5_Hs / HERVK-int.')
    unresolved = [v for v, r in calls.items() if r['evidence'] == 'UNRESOLVED']
    if unresolved:
        warnings.append(f'{len(unresolved)} records UNRESOLVED (kept and flagged, '
                        'not dropped): ' + ', '.join(sorted(unresolved)[:5]))
    split = [r['locus_id'] for r in loci.values()
             if 'LOCUS_SPLIT_BY_HUMAN_FILTER' in r['flags']]
    if split:
        warnings.append(f'{len(split)} loci split by the --human FILTER="PASS" '
                        'requirement: ' + ', '.join(split[:5]))

    # ---- report ----------------------------------------------------------
    print(f'HERV-K v2 assertions — {checked} records checked against '
          f'{len(expected)} expectations')
    print(f'  candidates classified : {len(calls)}')
    print(f'  loci                  : {len(loci)}  '
          f'({len(merge_loci)} flagged for merge)')
    print(f'  null_prov             : {len(null_prov)}  {", ".join(null_prov)}')
    by_ev = {}
    for r in calls.values():
        by_ev[r['evidence']] = by_ev.get(r['evidence'], 0) + 1
    print('  evidence              : ' +
          ', '.join(f'{k}={v}' for k, v in sorted(by_ev.items())))

    for w in warnings:
        print(f'\n  NOTE: {w}')
    if failures:
        print(f'\n  {len(failures)} FAILURE(S):')
        for f in failures:
            print(f'    - {f}')
        sys.exit(1)
    print('\n  PASS — all expectations met.')


if __name__ == '__main__':
    main()
