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
import gzip
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


def read_consolidated(path):
    """Locus records and per-locus INFO from a consolidated VCF.

    Returns (by_locus, n_records, header_ids). by_locus is keyed on
    HERVK_LOCUS, so it holds merged records and unmerged ones that carry the
    locus annotation alike, which is what lets the same checks run over the
    discovery and the graph file.
    """
    opener = gzip.open if path.endswith('.gz') else open
    by_locus, n, header_ids = {}, 0, set()
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('##INFO=<ID='):
                header_ids.add(line.split('ID=', 1)[1].split(',', 1)[0])
                continue
            if line.startswith('#'):
                continue
            n += 1
            f = line.rstrip('\n').split('\t')
            info = dict(kv.split('=', 1) if '=' in kv else (kv, '')
                        for kv in f[7].split(';'))
            if 'HERVK_LOCUS' in info:
                by_locus[info['HERVK_LOCUS']] = info
    return by_locus, n, header_ids


def check_locus_layer(path, label, failures, warnings):
    """Assertions that must hold in every consolidated VCF, whichever stage.

    A user filtering on HERVK_MEI has to get the same set of loci from either
    file. The counts here are the CaG cohort: 31 HERV-K loci, 21 of them
    insertion polymorphisms.
    """
    by_locus, n_records, header_ids = read_consolidated(path)

    # Every field the records use has to be declared. An undeclared flag makes
    # bcftools guess Type=String, and `-i 'INFO/HERVK_MEI=1'` then matches
    # nothing without erroring.
    for tag in ('HERVK_LOCUS', 'HERVK_LOCUS_N', 'HERVK_MEI', 'HERVK_SOLO_PROV',
                'HERVK_CNV', 'HERVK_LOCUS_TYPE'):
        if tag not in header_ids:
            failures.append(f'{label}: {tag} is not declared in the header')

    mei = {l for l, i in by_locus.items() if 'HERVK_MEI' in i}
    sp = {l for l, i in by_locus.items() if 'HERVK_SOLO_PROV' in i}
    cnv = {l for l, i in by_locus.items() if 'HERVK_CNV' in i}
    counts = {'loci': len(by_locus), 'mei': len(mei),
              'solo_prov': len(sp), 'cnv': len(cnv),
              'incomplete': sum(1 for i in by_locus.values()
                                if 'HERVK_LOCUS_INCOMPLETE' in i)}
    for key, want in (('loci', 31), ('mei', 21), ('solo_prov', 9), ('cnv', 3),
                      ('incomplete', 3)):
        if counts[key] != want:
            warnings.append(f'{label}: {key}={counts[key]}, expected {want} '
                            'on the CaG cohort')

    # chr6 is the locus a single category cannot describe.
    m = by_locus.get('HERVK_chr6_78894316')
    if m is None:
        failures.append(f'{label}: HERVK_chr6_78894316 carries no locus id')
    else:
        for tag in ('HERVK_SOLO_PROV', 'HERVK_CNV'):
            if tag not in m:
                failures.append(f'{label}: chr6 should carry {tag}; it '
                                'segregates a solo LTR, a provirus and a '
                                'two-unit allele')
        if 'HERVK_MEI' in m:
            failures.append(f'{label}: chr6 has no null allele and must not '
                            'carry HERVK_MEI')

    # The two loci the --human filter cuts in half have to say so, or their
    # allele frequencies read as covering the locus.
    for lid in ('HERVK_chr7_4699714', 'HERVK_chr8_7552031'):
        i = by_locus.get(lid)
        if i is None:
            failures.append(f'{label}: {lid} carries no locus id')
        elif 'HERVK_LOCUS_INCOMPLETE' not in i:
            failures.append(f'{label}: {lid} lost members to the --human '
                            'filter and must carry HERVK_LOCUS_INCOMPLETE')
        elif not i.get('HERVK_ALLELE_SET'):
            failures.append(f'{label}: {lid} is incomplete but names no '
                            'HERVK_ALLELE_SET, so a reader cannot see what '
                            'else segregates there')
    return counts, n_records


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
    disc_counts = graph_counts = None

    # hervk_loci.tsv gained the locus properties, and both consolidated VCFs
    # read them from here rather than recomputing, so a missing column means
    # every downstream flag is missing too.
    for col in ('locus_type', 'mei', 'solo_prov', 'cnv'):
        if loci and col not in next(iter(loci.values())):
            failures.append(f'hervk_loci.tsv has no {col} column')

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
        # j, and the allele states either side. The states are the point of the
        # copy-number work: an array locus has to say how many units each
        # allele carries, not just that something proviral is segregating.
        if want.get('expect_j') and got.get('j', '') != want['expect_j']:
            failures.append(f'{vid}: j={got.get("j")}, expected {want["expect_j"]}')
        for col, field in (('expect_allele_ref', 'allele_ref'),
                           ('expect_allele', 'allele')):
            if want.get(col) and got.get(field, '') != want[col]:
                failures.append(f'{vid}: {field}={got.get(field)}, '
                                f'expected {want[col]}')

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

    # ---- 3. copy-number loci --------------------------------------------
    # Discovery genotypes stay on these records: they come from
    # haplotype-resolved alignments and carry the allele frequencies. It is the
    # GRAPH genotypes that get withheld, in hervk_reconcile.py consolidate,
    # because the ALT path repeats sequence the reference already carries.
    tandem = sorted(v for v, r in calls.items() if r['class'] == 'copy_number')
    for vid in tandem:
        ref_st = calls[vid].get('allele_ref', '')
        alt_st = calls[vid].get('allele', '')
        if not (ref_st.startswith('prov_x') or alt_st.startswith('prov_x')):
            failures.append(f'{vid}: classed copy_number but neither allele is '
                            f'a prov_xN state ({ref_st} -> {alt_st})')

    # ---- 3b. the consolidated discovery VCF ------------------------------
    disc_cons = os.path.join(d, 'pangenome.human.consolidated.vcf')
    if os.path.exists(disc_cons):
        disc_counts, disc_n = check_locus_layer(
            disc_cons, 'pangenome.human.consolidated.vcf', failures, warnings)
        # The unmerged file has to keep every record: it induces the graph.
        raw = os.path.join(d, 'pangenome.human.vcf')
        if os.path.exists(raw):
            with open(raw) as fh:
                raw_n = sum(1 for l in fh if not l.startswith('#'))
            if disc_n >= raw_n:
                failures.append(
                    f'the consolidated discovery VCF has {disc_n} records and '
                    f'the unmerged one {raw_n}; merging should reduce the count')
    else:
        failures.append(f'missing {disc_cons} -- hervk_annotate did not write '
                        'the consolidated discovery VCF')

    # ---- 4. stage E, when it ran ----------------------------------------
    gt_dir = os.path.join(args.outdir, '4_Genotyping')
    cons = os.path.join(gt_dir, 'GraffiTE.merged.genotypes.human.vcf.gz')
    cons_plain = cons[:-3]
    cons_path = cons if os.path.exists(cons) else (
        cons_plain if os.path.exists(cons_plain) else None)
    stage_e = {}
    if cons_path:
        import gzip
        opener = gzip.open if cons_path.endswith('.gz') else open
        merged, masked_ok, gt_masked = {}, True, []
        with opener(cons_path, 'rt') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                f = line.rstrip('\n').split('\t')
                info = dict(kv.split('=', 1) if '=' in kv else (kv, '')
                            for kv in f[7].split(';'))
                if f[2].startswith('HERVK_'):
                    merged[f[2]] = info
                if 'HERVK_GT_MASKED' in info:
                    gt_masked.append(f[2])
                    called = [g.split(':')[0] for g in f[9:]]
                    if any(c.replace('|', '/').replace('.', '').replace('/', '')
                           for c in called):
                        masked_ok = False
        stage_e = {'merged': merged, 'gt_masked': gt_masked}

        graph_counts, _ = check_locus_layer(
            cons_path, 'GraffiTE.merged.genotypes.human.vcf.gz',
            failures, warnings)
        if disc_counts and graph_counts and disc_counts != graph_counts:
            failures.append(
                'the two delivered VCFs disagree about the locus layer: '
                f'discovery {disc_counts}, graph {graph_counts}. A user '
                'filtering on HERVK_MEI must get the same loci from either.')

        for lid in ('HERVK_chr11_101704640', 'HERVK_chr12_55299985',
                    'HERVK_chr6_78894316'):
            if lid not in merged:
                failures.append(f'{lid} was not consolidated in {cons_path}')

        # chr6 is where masking per allele instead of per locus shows up: the
        # solo allele keeps its graph genotypes and only the two-unit allele
        # is withheld.
        if 'HERVK_chr6_78894316' in merged:
            m6 = merged['HERVK_chr6_78894316']
            if m6.get('HERVK_ALLELE_NOGT') != 'prov_x2':
                failures.append(
                    f"chr6 HERVK_ALLELE_NOGT = {m6.get('HERVK_ALLELE_NOGT')}, "
                    'expected prov_x2 -- the two-unit allele is the one the '
                    'graph cannot count')
            if not (m6.get('HERVK_AC', '').startswith('8,')):
                failures.append(
                    f"chr6 HERVK_AC = {m6.get('HERVK_AC')}, expected the solo "
                    'allele at 8; the graph finds every solo carrier')
            if m6.get('HERVK_AC_DISC') != '8,1':
                warnings.append(
                    f"chr6 HERVK_AC_DISC = {m6.get('HERVK_AC_DISC')}, expected "
                    '8,1 from the assemblies')
        if 'HERVK_chr11_101704640' in merged:
            m = merged['HERVK_chr11_101704640']
            if m.get('HERVK_AC') != '23' or m.get('HERVK_AN') != '40':
                failures.append(
                    f"chr11 AC/AN = {m.get('HERVK_AC')}/{m.get('HERVK_AN')}, "
                    'expected 23/40 -- graph and assemblies agree exactly here, '
                    'so a mismatch means the dosage resolution changed')
            if 'HERVK_DISC_CONCORDANT' not in m:
                warnings.append('chr11 not flagged HERVK_DISC_CONCORDANT; the '
                                'graph and the assemblies disagree there now')
        if 'HERVK_chr12_55299985' in merged:
            m = merged['HERVK_chr12_55299985']
            if m.get('HERVK_N_PLOIDY_EXCEEDED') != '2':
                warnings.append(
                    f"chr12 ploidy violations = {m.get('HERVK_N_PLOIDY_EXCEEDED')}, "
                    'expected 2 (the third allele flattened by bcftools norm -m-)')
        if gt_masked and not masked_ok:
            failures.append('a HERVK_GT_MASKED record still carries called '
                            'genotypes')

    # ---- 5. things worth seeing but not failing on -----------------------
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
    print(f'HERV-K v3 assertions — {checked} records checked against '
          f'{len(expected)} expectations')
    print(f'  copy_number           : {len(tandem)}')
    for v in tandem:
        r = calls[v]
        print(f'      {v}  {r.get("allele_ref")} -> {r.get("allele")}  '
              f'(ref_units={r.get("n_units_ref")}, period={r.get("unit_bp")}, '
              f'k={r.get("k")}, j={r.get("j")})')
    if disc_counts:
        print('  stage 3: consolidated : ' +
              '  '.join(f'{k}={v}' for k, v in disc_counts.items()))
    if stage_e:
        print(f'  stage E: consolidated : {len(stage_e["merged"])}  '
              f'({", ".join(sorted(stage_e["merged"]))})')
        print(f'  stage E: GT withheld  : {len(stage_e["gt_masked"])}')
    if graph_counts:
        print('  stage E: locus layer  : ' +
              '  '.join(f'{k}={v}' for k, v in graph_counts.items()))
    else:
        print('  stage E               : not run (set GENOTYPED_VCF to test it)')
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
