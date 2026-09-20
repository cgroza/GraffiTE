#!/usr/bin/env python3
"""Measure one cell's output against truth.tsv. Writes OBSERVED.tsv and CALIBRATION.log.

The measurement that has to happen before any expectation is frozen. Nothing
here asserts on pipeline behaviour. A planted site that never reached
pangenome.vcf is written down as such, and a designed field that did not fire is
written down as None.

  observe.py --workdir <WORKDIR> --cell spine
"""

import argparse
import glob
import gzip
import os
import re
import sys
from collections import defaultdict

POS_TOL = 200        # bcftools norm -f left-aligns the truvari output, so a
                     # record can sit left of where the site was planted. The
                     # width is a guess; widen it if sites go unmatched.
SVLEN_REL_TOL = 0.30

SVA_E_VNTR = (428, 864)   # annotate_vcf.R:133-134; admissible band is 429..863
SPAN_CUTOFF = 0.80        # nextflow.config: repeat_span_cutoff, applied with '>'


def op(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def read_vcf(path):
    """[(chrom, pos, id, ref, alt, {info}, line)]"""
    recs = []
    if not os.path.exists(path):
        return recs
    with op(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 8:
                continue
            info = {}
            for kv in f[7].split(';'):
                if '=' in kv:
                    k, v = kv.split('=', 1)
                    info[k] = v
                elif kv:
                    info[kv] = 'TRUE'
            recs.append((f[0], int(f[1]), f[2], f[3], f[4], info, f))
    return recs


def svlen_of(ref, alt, info):
    if 'SVLEN' in info:
        try:
            return abs(int(info['SVLEN']))
        except ValueError:
            pass
    if alt.startswith('<'):
        return 0
    return abs(len(alt) - len(ref))


def read_rm_out(path):
    """RepeatMasker .out -> list of dicts. Columns per annotate_vcf.R:10-14."""
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i < 3:
                continue
            p = line.split()
            if len(p) < 15:
                continue
            def num(x):
                return int(re.sub(r'[()]', '', x))
            rows.append(dict(
                sw=int(p[0]), div=float(p[1]), qry_id=p[4],
                qry_start=int(p[5]), qry_end=int(p[6]), strand=p[8],
                repeat_id=p[9], matching_class=p[10],
                in_repeat_start=num(p[11]), in_repeat_end=num(p[12]),
                in_repeat_left=num(p[13]), link_id=p[14],
                star=(len(p) > 15 and p[15] == '*'),
            ))
            # annotate_vcf.R:50-51
            rows[-1]['target_start'] = (rows[-1]['in_repeat_start'] if p[8] == '+'
                                        else rows[-1]['in_repeat_left'])
            rows[-1]['target_end'] = rows[-1]['in_repeat_end']
    return rows


def read_truth(path):
    rows = []
    with open(path) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            v = line.rstrip('\n').split('\t')
            r = dict(zip(hdr, v))
            r['pos'] = int(r['pos'])
            r['svlen'] = int(r['svlen'])
            rows.append(r)
    return rows


def match(truth_row, recs):
    """Best VCF record for a planted site: same contig, near position, near length."""
    best, best_key = None, None
    for r in recs:
        chrom, pos, vid, ref, alt, info, _ = r
        if chrom != truth_row['contig']:
            continue
        d = abs(pos - truth_row['pos'])
        if d > POS_TOL:
            continue
        L = svlen_of(ref, alt, info)
        if truth_row['svlen'] and L:
            if abs(L - truth_row['svlen']) / truth_row['svlen'] > SVLEN_REL_TOL:
                continue
        key = (d, abs(L - truth_row['svlen']))
        if best_key is None or key < best_key:
            best, best_key = r, key
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workdir', required=True)
    ap.add_argument('--cell', default='spine')
    ap.add_argument('--out', default=None)
    ap.add_argument('--log', default=None)
    a = ap.parse_args()

    W = a.workdir
    RUN = os.path.join(W, 'runs', a.cell)
    here = os.path.dirname(os.path.abspath(__file__))
    out_tsv = a.out or os.path.join(here, 'OBSERVED.tsv')
    out_log = a.log or os.path.join(here, 'CALIBRATION.log')

    truth = read_truth(os.path.join(W, 'build', 'truth.tsv'))

    # ---- sources -------------------------------------------------------
    pangenome = read_vcf(os.path.join(RUN, '3_TSD_search', 'pangenome.vcf'))
    trusted   = read_vcf(os.path.join(RUN, '3_TSD_search', 'pangenome.trusted.vcf'))
    human     = read_vcf(os.path.join(RUN, '3_TSD_search', 'pangenome.human.vcf'))
    svs       = read_vcf(os.path.join(RUN, '1_SV_search', 'SVs.vcf'))

    # every chunk's pre-cutoff annotated VCF, so a record the span filter dropped
    # is still visible
    prefilter = []
    for p in sorted(glob.glob(os.path.join(RUN, '2_Repeat_Filtering', '*',
                                           'genotypes_repmasked.vcf.gz'))):
        prefilter += read_vcf(p)
    postfilter = []
    for p in sorted(glob.glob(os.path.join(RUN, '2_Repeat_Filtering', '*',
                                           'genotypes_repmasked_filtered.vcf'))):
        postfilter += read_vcf(p)

    rm_rows = []
    for p in sorted(glob.glob(os.path.join(RUN, '2_Repeat_Filtering', '*',
                                           'repeatmasker_dir', '*.out'))):
        rm_rows += read_rm_out(p)
    rm_by_qry = defaultdict(list)
    for r in rm_rows:
        rm_by_qry[r['qry_id']].append(r)

    all_sw = [r['sw'] for r in rm_rows]
    min_sw = min(all_sw) if all_sw else None

    cols = ['name', 'kind', 'contig', 'pos', 'truth_svlen', 'truth_tsd', 'carriers',
            'in_svs', 'in_prefilter', 'in_postfilter', 'in_pangenome',
            'in_trusted', 'in_human',
            'vcf_id', 'vcf_pos', 'vcf_svlen',
            'n_hits', 'repeat_ids', 'matching_classes', 'fragmts', 'strands',
            'RM_hit_IDs', 'L1_5PINV', 'total_repeat_span', 'ULTRA_TR_span',
            'TSD', 'polyA',
            'rm_n_rows', 'rm_link_ids', 'rm_best_sw', 'rm_min_sw', 'rm_repeat_ids',
            'rm_target_start', 'rm_target_end', 'rm_strands']

    rows = []
    for t in truth:
        rec_pan = match(t, pangenome)
        rec_pre = match(t, prefilter)
        rec_post = match(t, postfilter)
        rec_svs = match(t, svs)
        rec_tru = match(t, trusted)
        rec_hum = match(t, human)
        src = rec_pan or rec_post or rec_pre
        info = src[5] if src else {}

        # RepeatMasker rows are keyed by the VCF ID the chunk carried
        qid = None
        for cand in (rec_pre, rec_post, rec_pan):
            if cand and cand[2] in rm_by_qry:
                qid = cand[2]
                break
        hits = rm_by_qry.get(qid, []) if qid else []

        def g(k, d='None'):
            return info.get(k, d)

        rows.append({
            'name': t['name'], 'kind': t['kind'], 'contig': t['contig'],
            'pos': t['pos'], 'truth_svlen': t['svlen'], 'truth_tsd': t['tsd'],
            'carriers': t['carriers'],
            'in_svs': 'yes' if rec_svs else 'no',
            'in_prefilter': 'yes' if rec_pre else 'no',
            'in_postfilter': 'yes' if rec_post else 'no',
            'in_pangenome': 'yes' if rec_pan else 'no',
            'in_trusted': 'yes' if rec_tru else 'no',
            'in_human': 'yes' if rec_hum else 'no',
            'vcf_id': src[2] if src else 'None',
            'vcf_pos': src[1] if src else 'None',
            'vcf_svlen': svlen_of(src[3], src[4], info) if src else 'None',
            'n_hits': g('n_hits'), 'repeat_ids': g('repeat_ids'),
            'matching_classes': g('matching_classes'), 'fragmts': g('fragmts'),
            'strands': g('RM_hit_strands', g('strands')),
            'RM_hit_IDs': g('RM_hit_IDs'), 'L1_5PINV': g('L1_5PINV'),
            'total_repeat_span': g('total_repeat_span'),
            'ULTRA_TR_span': g('ULTRA_TR_span'),
            'TSD': g('TSD'), 'polyA': g('polyA'),
            'rm_n_rows': len(hits),
            'rm_link_ids': ','.join(sorted({h['link_id'] for h in hits})) or 'None',
            'rm_best_sw': max((h['sw'] for h in hits), default='None'),
            'rm_min_sw': min((h['sw'] for h in hits), default='None'),
            'rm_repeat_ids': ','.join(h['repeat_id'] for h in hits) or 'None',
            'rm_target_start': ','.join(str(h['target_start']) for h in hits) or 'None',
            'rm_target_end': ','.join(str(h['target_end']) for h in hits) or 'None',
            'rm_strands': ''.join(h['strand'] for h in hits) or 'None',
        })

    with open(out_tsv, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r[c]) for c in cols) + '\n')

    by_name = {r['name']: r for r in rows}

    # ---- calibration ---------------------------------------------------
    L = []
    def say(s=''):
        L.append(s)

    say(f'cell            : {a.cell}')
    say(f'run dir         : {RUN}')
    say(f'planted sites   : {len(truth)}')
    say(f'SVs.vcf         : {len(svs)} records')
    say(f'pre-cutoff      : {len(prefilter)} records over '
        f'{len(glob.glob(os.path.join(RUN, "2_Repeat_Filtering", "*")))} chunks')
    say(f'post-cutoff     : {len(postfilter)} records')
    say(f'pangenome.vcf   : {len(pangenome)} records')
    say(f'trusted.vcf     : {len(trusted)} records')
    say(f'human.vcf       : {len(human)} records')
    say(f'RepeatMasker    : {len(rm_rows)} hit rows, min SW {min_sw}')
    say()

    verdict = []

    # 1. every designed element produced at least one RepeatMasker hit
    no_hit = [r['name'] for r in rows if r['rm_n_rows'] == 0]
    if no_hit:
        verdict.append(('FIXTURE', 'every designed element has >=1 RepeatMasker hit',
                        f'{len(no_hit)} with none: ' + ','.join(no_hit)))
    else:
        verdict.append(('ok', 'every designed element has >=1 RepeatMasker hit', ''))

    # 2. no planted copy within 15 % of the reporting floor.
    #
    # Read literally against "the lowest observed SW score" the check cannot
    # fail usefully. Every row of indels.fa is a planted copy, so the weakest
    # planted copy is the minimum by construction and sits within 15 % of
    # itself. The floor that decides anything is the score RepeatMasker will
    # not report below. bin/repmask_vcf.sh invokes RepeatMasker with no
    # -cutoff, which I checked; that the resulting default is 225 comes from
    # RepeatMasker's documentation and is not verified against this container.
    # Confirm it against VERSIONS.txt before freezing. Both readings are
    # printed.
    RM_DEFAULT_CUTOFF = 225
    best = [(r['name'], r['rm_best_sw']) for r in rows
            if isinstance(r['rm_best_sw'], int)]
    if best:
        weakest = min(best, key=lambda x: x[1])
        say(f'lowest SW over all RepeatMasker hit rows : {min_sw}')
        say(f'weakest planted copy (best hit)          : {weakest[0]} = {weakest[1]}')
        say(f'RepeatMasker reporting cutoff (default)  : {RM_DEFAULT_CUTOFF}')
        say()
        close = [(n, v) for n, v in best if v < RM_DEFAULT_CUTOFF * 1.15]
        if close:
            verdict.append(('FIXTURE',
                            'no planted copy within 15 % of the reporting cutoff (225)',
                            'within 15 %: ' + ','.join(f'{n}={v}' for n, v in close)))
        else:
            verdict.append(('ok',
                            'no planted copy within 15 % of the reporting cutoff (225)',
                            f'weakest is {weakest[0]}={weakest[1]}'))
    else:
        verdict.append(('None',
                        'no planted copy within 15 % of the reporting cutoff (225)',
                        'no RepeatMasker hits at all'))

    # 3. the span ladder has exactly one cut point, bracketing 0.80
    ladder = [r for r in rows if r['kind'] == 'span_ladder']
    ladder.sort(key=lambda r: -int(re.search(r'_span_(\d+)', r['name']).group(1)))
    say('span ladder (designed fraction -> observed total_repeat_span, survival):')
    for r in ladder:
        say(f"  {r['name']:<18} span={r['total_repeat_span']:<8} "
            f"ultra={r['ULTRA_TR_span']:<8} postfilter={r['in_postfilter']:<4} "
            f"pangenome={r['in_pangenome']}")
    say()
    surv = [(r['name'], r['in_pangenome'] == 'yes') for r in ladder]
    cuts = sum(1 for i in range(len(surv) - 1) if surv[i][1] != surv[i + 1][1])
    if cuts == 1:
        i = next(i for i in range(len(surv) - 1) if surv[i][1] != surv[i + 1][1])
        lo, hi = ladder[i], ladder[i + 1]
        try:
            a_s, b_s = float(lo['total_repeat_span']), float(hi['total_repeat_span'])
            brackets = a_s > SPAN_CUTOFF >= b_s
        except ValueError:
            a_s = b_s = None
            brackets = False
        verdict.append(('ok' if brackets else 'FIXTURE',
                        'span ladder has exactly one cut point bracketing 0.80',
                        f'cut between {lo["name"]}({a_s}) and {hi["name"]}({b_s})'))
    else:
        verdict.append(('FIXTURE', 'span ladder has exactly one cut point bracketing 0.80',
                        f'{cuts} cut points'))

    # 4. hervk_loci.tsv has a row with n_records > 1
    hl = os.path.join(RUN, '3_TSD_search', 'hervk_loci.tsv')
    if not os.path.exists(hl):
        hl_alt = glob.glob(os.path.join(RUN, '**', 'hervk_loci.tsv'), recursive=True)
        hl = hl_alt[0] if hl_alt else hl
    if os.path.exists(hl):
        with open(hl) as fh:
            hdr = fh.readline().rstrip('\n').split('\t')
            hrows = [dict(zip(hdr, l.rstrip('\n').split('\t'))) for l in fh if l.strip()]
        say(f'hervk_loci.tsv  : {len(hrows)} rows, columns {hdr}')
        key = 'n_records' if 'n_records' in hdr else None
        multi = [r for r in hrows if key and r.get(key, '0').isdigit()
                 and int(r[key]) > 1]
        if key is None:
            verdict.append(('FIXTURE', 'hervk_loci.tsv has a row with n_records>1',
                            f'no n_records column; columns are {hdr}'))
        elif multi:
            verdict.append(('ok', 'hervk_loci.tsv has a row with n_records>1',
                            f'{len(multi)} such rows'))
        else:
            verdict.append(('FIXTURE', 'hervk_loci.tsv has a row with n_records>1',
                            'truvari collapse merged the designed multi-record locus'))
    else:
        verdict.append(('None', 'hervk_loci.tsv has a row with n_records>1',
                        'hervk_loci.tsv not produced'))
    say()

    # ---- the four items on the README's calibration list ----------------
    say('-- L1_5PINV locus (A09_l1_twin_primed_Cplus) --')
    r = by_name.get('A09_l1_twin_primed_Cplus')
    if r:
        say(f"  reached pangenome : {r['in_pangenome']}")
        say(f"  RM link IDs       : {r['rm_link_ids']}  (rows={r['rm_n_rows']})")
        say(f"  RM strands        : {r['rm_strands']}")
        say(f"  RM repeat ids     : {r['rm_repeat_ids']}")
        say(f"  INFO/RM_hit_IDs   : {r['RM_hit_IDs']}")
        say(f"  INFO/L1_5PINV     : {r['L1_5PINV']}")
        say(f"  VERDICT           : " +
            ('set' if r['L1_5PINV'] not in ('None', '.', '') else 'None'))
    else:
        say('  None -- site absent from truth.tsv')
    r2 = by_name.get('A10_l1_plusC_negative')
    if r2:
        say(f"  negative control A10: strands={r2['rm_strands']} "
            f"L1_5PINV={r2['L1_5PINV']}")
    say()

    say('-- SVA VNTR-only locus (A12_sva_vntr_only) --')
    r = by_name.get('A12_sva_vntr_only')
    if r:
        say(f"  reached pangenome : {r['in_pangenome']}")
        say(f"  target_start      : {r['rm_target_start']}")
        say(f"  target_end        : {r['rm_target_end']}")
        say(f"  repeat_ids        : {r['repeat_ids']}")
        say(f"  matching_classes  : {r['matching_classes']}")
        ok = None
        try:
            ts = [int(x) for x in str(r['rm_target_start']).split(',') if x.strip('-').isdigit()]
            te = [int(x) for x in str(r['rm_target_end']).split(',') if x.strip('-').isdigit()]
            if ts and te:
                ok = (min(ts) > SVA_E_VNTR[0] and max(te) < SVA_E_VNTR[1])
        except Exception:
            ok = None
        say(f"  strictly inside SVA_E 429..863 : "
            f"{'yes' if ok else ('no' if ok is False else 'None')}")
        say(f"  (VNTR_only) in repeat_ids      : "
            f"{'yes' if 'VNTR_only' in str(r['repeat_ids']) else 'no'}")
    else:
        say('  None -- site absent from truth.tsv')
    say()

    say('-- polyA and TSD probes --')
    for n in ['A01_alu_tsd_polyA', 'A02_alu_no_polyA_no_tsd', 'A03_alu_polyA_too_short',
              'A04_alu_polyA_beyond_slack', 'A05_alu_revcomp', 'A07_alu_full_length_clean',
              'A08_alu_no_tsd']:
        r = by_name.get(n)
        if r:
            say(f"  {n:<28} TSD={r['TSD']:<12} polyA={r['polyA']:<6} "
                f"n_hits={r['n_hits']:<3} strands={r['strands']}")
    say()

    say('-- subsets --')
    for n, r in by_name.items():
        if r['in_pangenome'] == 'no':
            continue
    miss = [r['name'] for r in rows if r['in_pangenome'] == 'no']
    say(f"  planted sites that never reached pangenome.vcf ({len(miss)}): "
        + (','.join(miss) if miss else 'none'))
    say()

    say('== verdict ==')
    for status, what, detail in verdict:
        say(f'  [{status:^7}] {what}' + (f' -- {detail}' if detail else ''))

    text = '\n'.join(L) + '\n'
    with open(out_log, 'w') as fh:
        fh.write(text)
    sys.stdout.write(text)
    print(f'\nwrote {out_tsv} and {out_log}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
