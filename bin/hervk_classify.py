#!/usr/bin/env python3
"""
HERV-K (HML-2) SV polymorphism classifier for GraffiTE --human runs.

Classifies each TE-annotated SV into one of five hypotheses based on its
size and HERV-K family content:

    H_C : null     <-> solo-LTR        (canonical |SVLEN| ~  968 bp)
    H_T : truncated proviral           (1500 <= |SVLEN| <= 8000)
    H_B : solo-LTR <-> proviral        (canonical |SVLEN| ~ 8504 bp)
    H_A : null     <-> proviral        (canonical |SVLEN| ~ 9472 bp)
    H_X : non-transposition / other    (flat background)

Inputs/outputs are GraffiTE-style VCF and presence-absence TSV. The tool
adds INFO/HERVK_CLASS, INFO/HERVK_PMAP, INFO/HERVK_LAMBDA, INFO/HERVK_NU
to the VCF, plus an optional FORMAT/HERVK_AS field giving the per-sample
per-haplotype allelic state (e.g. "solo|null"). When --strict is set,
candidate rows whose MAP class is "other" or whose MAP posterior is below
the threshold are dropped.

Usage examples:

    # Annotate the main VCF + TSV; emit summary
    hervk_classify.py \\
        --vcf-in pangenome.vcf --vcf-out pangenome.hervk.vcf \\
        --tsv-in pangenome.presence-absence.tsv \\
        --tsv-out pangenome.presence-absence.hervk.tsv \\
        --summary hervk_polymorphism_summary.md

    # Strict-filter the trusted/human VCF + TSV
    hervk_classify.py --strict \\
        --vcf-in pangenome.trusted.human.vcf \\
        --vcf-out pangenome.trusted.human.hervk.vcf \\
        --tsv-in pangenome.presence-absence_human.tsv \\
        --tsv-out pangenome.presence-absence_human.hervk.tsv
"""

import argparse
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict


# -------- Reference architecture --------
LTR_LEN, INT_LEN = 968, 7536
SOLO_PROV = LTR_LEN + INT_LEN          # 8504
NULL_PROV = 2 * LTR_LEN + INT_LEN      # 9472

EXPECTED = {
    'C': {'s': LTR_LEN,    'lam': LTR_LEN,   'nu': 0       },
    'B': {'s': SOLO_PROV,  'lam': LTR_LEN,   'nu': INT_LEN },
    'A': {'s': NULL_PROV,  'lam': 2*LTR_LEN, 'nu': INT_LEN },
}

LTR_FAMILY = {'LTR5_Hs', 'LTR5A', 'LTR5B'}
INT_FAMILY = {'HERVK-int'}
SVA_FAMILIES = {'SVA_A', 'SVA_B', 'SVA_C', 'SVA_D', 'SVA_E', 'SVA_F'}

# Defaults — see HERVK.config.json for runtime overrides.
DEFAULTS = {
    "sigmas": {
        "s_C": 30.0,    # solo-LTRs are very length-uniform
        "s_B": 800.0,   # allow ~10% INT truncation
        "s_A": 800.0,
        "lam": 300.0,   # absorbs (x)-merged LTR/INT calls
        "nu":  1000.0,
        "t":   200.0,   # SV must be mostly HERV-K
    },
    "priors": {"C": 0.55, "T": 0.08, "B": 0.05, "A": 0.02, "X": 0.30},
    "t_min": 1500,
    "t_max": 8000,
    # Background flat density support for H_X.
    "s_range": 30000.0,
    "lam_range": 5000.0,
    "nu_range": 30000.0,
    # SVA-as-LTR-mimic window for n_hits==2 HERVK+SVA cases.
    # Empirical observation: SVA "steals" ~328 bp of LTR5_Hs annotation due
    # to LTR5_Hs/SVA homology; treat SVA bp as LTR-equivalent in this window.
    "sva_mimic_min": 250,
    "sva_mimic_max": 400,
    # Strict-mode threshold (only applied when --strict).
    "pmap_min": 0.90,
    # Polyallelic flag window (bp).
    "polyallelic_window": 100,
    # FORMAT/HERVK_AS: emit "?|?" when MAP posterior below this.
    "as_min_posterior": 0.90,
}

ALLELE_INTERPRETATION = {
    'C': {'ref': 'null',  'alt': 'solo'           },
    'T': {'ref': 'null',  'alt': 'truncated_prov' },
    'B': {'ref': 'solo',  'alt': 'prov'           },
    'A': {'ref': 'null',  'alt': 'prov'           },
    'X': {'ref': '?',     'alt': '?'              },
}

CLASS_LABEL = {
    'C': 'null_solo',
    'T': 'truncated_prov',
    'B': 'solo_prov',
    'A': 'null_prov',
    'X': 'other',
}


# -------- Config loading --------
def load_config(path):
    cfg = {k: (dict(v) if isinstance(v, dict) else v)
           for k, v in DEFAULTS.items()}
    if path:
        with open(path) as fh:
            user = json.load(fh)
        for k, v in user.items():
            if k in cfg and isinstance(cfg[k], dict) and isinstance(v, dict):
                cfg[k].update(v)
            else:
                cfg[k] = v
    return cfg


# -------- Annotation parsing --------
def parse_hits(match_lengths, repeat_ids, matching_classes, n_hits, cfg):
    """Return (lambda, nu, sva_neutral_bp).

    lambda  = bp matching LTR family (+ SVA bp when SVA mimics LTR5_Hs)
    nu      = bp matching INT family
    sva_neutral_bp = SVA bp that is *not* counted as LTR-mimic; subtracted
                     from s for the coverage term so a true HERV-K event
                     accompanied by a separate large SVA insertion is not
                     incorrectly pushed into H_X.
    """
    if not match_lengths or not repeat_ids:
        return 0.0, 0.0, 0.0
    try:
        lens = [float(x) for x in str(match_lengths).split(',')]
    except ValueError:
        return 0.0, 0.0, 0.0
    ids = [r.strip().replace('(x)', '') for r in str(repeat_ids).split(',')]
    if len(lens) != len(ids):
        return 0.0, 0.0, 0.0

    lam = sum(L for L, r in zip(lens, ids) if r in LTR_FAMILY)
    nu  = sum(L for L, r in zip(lens, ids) if r in INT_FAMILY)
    sva_neutral = 0.0

    # Special case: n_hits == 2 with HERVK + SVA. Add SVA bp to lambda when
    # in the empirical mimic window; otherwise treat as neutral background
    # (subtracted from s for coverage term only).
    try:
        n_hits_int = int(float(n_hits))
    except (TypeError, ValueError):
        n_hits_int = -1
    classes = set()
    if matching_classes:
        classes = {c.strip() for c in str(matching_classes).split(',')}
    is_hervk_sva = (n_hits_int == 2
                    and 'LTR/ERVK' in classes
                    and 'Retroposon/SVA' in classes)

    if is_hervk_sva:
        mimic_lo = cfg['sva_mimic_min']
        mimic_hi = cfg['sva_mimic_max']
        for L, r in zip(lens, ids):
            if r in SVA_FAMILIES:
                if mimic_lo <= L <= mimic_hi:
                    lam += L
                else:
                    sva_neutral += L

    return lam, nu, sva_neutral


# -------- Likelihoods --------
def log_likelihood_gaussian(s, lam, nu, hyp, cfg):
    sigmas = cfg['sigmas']
    e = EXPECTED[hyp]
    sig_s = sigmas[f's_{hyp}']
    t = lam + nu
    return (
        -0.5 * ((s   - e['s'])   / sig_s        ) ** 2
        -0.5 * ((lam - e['lam']) / sigmas['lam']) ** 2
        -0.5 * ((nu  - e['nu'])  / sigmas['nu'] ) ** 2
        -0.5 * ((s   - t)        / sigmas['t']  ) ** 2
    )


def log_likelihood_truncated(s, lam, nu, cfg):
    if s < cfg['t_min'] or s > cfg['t_max']:
        return -1e6
    sigmas = cfg['sigmas']
    t = lam + nu
    return (
        -math.log(cfg['t_max'] - cfg['t_min'])
        - 0.5 * ((lam - LTR_LEN) / sigmas['lam']) ** 2
        - 0.5 * ((s   - t)       / sigmas['t']  ) ** 2
    )


def log_background(cfg):
    return -math.log(cfg['s_range'] * cfg['lam_range'] * cfg['nu_range'])


def classify(s, lam, nu, cfg):
    """Return posterior dict over {'C','T','B','A','X'}."""
    priors = cfg['priors']
    lp = {k: log_likelihood_gaussian(s, lam, nu, k, cfg) + math.log(priors[k])
          for k in 'ABC'}
    lp['T'] = log_likelihood_truncated(s, lam, nu, cfg) + math.log(priors['T'])
    lp['X'] = log_background(cfg) + math.log(priors['X'])
    m = max(lp.values())
    e = {k: math.exp(v - m) for k, v in lp.items()}
    Z = sum(e.values())
    return {k: v / Z for k, v in e.items()}


def map_class(post):
    k = max(post, key=post.get)
    return k, post[k]


# -------- VCF I/O --------
INFO_HEADERS = [
    '##INFO=<ID=HERVK_CLASS,Number=1,Type=String,'
    'Description="HERV-K (HML-2) MAP class: null_solo|truncated_prov|'
    'solo_prov|null_prov|other|NA. Computed only for SVs in the HERV-K '
    'candidate set (LTR/ERVK with HML-2 family bp, or n_hits==2 with '
    'LTR/ERVK+Retroposon/SVA).">',
    '##INFO=<ID=HERVK_PMAP,Number=1,Type=Float,'
    'Description="Posterior probability of the HERV-K MAP class.">',
    '##INFO=<ID=HERVK_LAMBDA,Number=1,Type=Float,'
    'Description="bp matching HML-2 LTR family (LTR5_Hs/LTR5A/LTR5B); '
    'includes SVA bp when SVA mimics LTR5_Hs in n_hits==2 HERVK+SVA case.">',
    '##INFO=<ID=HERVK_NU,Number=1,Type=Float,'
    'Description="bp matching HML-2 internal family (HERVK-int).">',
]
FORMAT_HEADER = (
    '##FORMAT=<ID=HERVK_AS,Number=1,Type=String,'
    'Description="Per-haplotype HERV-K allelic state derived from GT and '
    'INFO/HERVK_CLASS (e.g. solo|null, prov|solo). Set to ?|? when MAP '
    'posterior is below threshold or class is other.">'
)


def parse_info(info):
    d = {}
    if not info or info == '.':
        return d
    for kv in info.split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k] = v
        else:
            d[kv] = ''
    return d


def info_to_str(d):
    parts = []
    for k, v in d.items():
        parts.append(k if v == '' else f'{k}={v}')
    return ';'.join(parts) if parts else '.'


def is_candidate(info_d, cfg):
    """Return True if SV qualifies for HERV-K classification."""
    matching_classes = info_d.get('matching_classes', '')
    if not matching_classes or matching_classes == 'NA':
        return False
    classes = {c.strip() for c in matching_classes.split(',')}
    try:
        n_hits = int(float(info_d.get('n_hits', '0')))
    except ValueError:
        return False
    if 'LTR/ERVK' not in classes:
        return False
    if n_hits == 1:
        return True
    if n_hits == 2 and 'Retroposon/SVA' in classes:
        return True
    return False


def interpret_haplotype_allele(allele, klass, svtype):
    """Map a single allele char ('0','1','.') to an allelic state."""
    if allele in ('.', ''):
        return '?'
    interp = ALLELE_INTERPRETATION[klass]
    if svtype == 'DEL':
        interp = {'ref': interp['alt'], 'alt': interp['ref']}
    return interp['alt'] if allele not in ('0',) else interp['ref']


def build_hervk_as(gt, klass, posterior, svtype, cfg):
    if klass == 'X' or posterior < cfg['as_min_posterior']:
        # Preserve haplotype shape (slash vs pipe, ploidy).
        sep_match = re.search(r'[/|]', gt or '')
        if sep_match:
            sep = sep_match.group(0)
            n = len(re.split(r'[/|]', gt))
            return sep.join(['?'] * n)
        return '?'
    alleles = re.split(r'([/|])', gt)  # keep separators
    out = []
    for tok in alleles:
        if tok in ('/', '|'):
            out.append(tok)
        else:
            out.append(interpret_haplotype_allele(tok, klass, svtype))
    return ''.join(out)


def process_vcf(vcf_in, vcf_out, cfg, strict, classifications):
    """Stream-rewrite the VCF, adding HERVK INFO and FORMAT fields.

    classifications: dict to be populated with vid -> dict for downstream
    TSV annotation and summary.
    """
    in_close = vcf_in != '-'
    out_close = vcf_out != '-'
    fin = open(vcf_in) if in_close else sys.stdin
    fout = open(vcf_out, 'w') if out_close else sys.stdout

    samples = []
    header_done = False

    for line in fin:
        if line.startswith('##'):
            fout.write(line)
            continue

        if line.startswith('#CHROM'):
            for h in INFO_HEADERS:
                fout.write(h + '\n')
            fout.write(FORMAT_HEADER + '\n')
            cols = line.rstrip('\n').split('\t')
            samples = cols[9:] if len(cols) > 9 else []
            fout.write(line)
            header_done = True
            continue

        if not header_done or not line.strip():
            fout.write(line)
            continue

        fields = line.rstrip('\n').split('\t')
        if len(fields) < 8:
            fout.write(line)
            continue

        chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]
        info_d = parse_info(info)

        if not is_candidate(info_d, cfg):
            fout.write(line)
            continue

        # Compute lambda/nu/coverage-adjusted s.
        lam, nu, sva_neutral = parse_hits(
            info_d.get('match_lengths', ''),
            info_d.get('repeat_ids', ''),
            info_d.get('matching_classes', ''),
            info_d.get('n_hits', '0'),
            cfg,
        )
        if lam + nu <= 0:
            fout.write(line)
            continue

        try:
            svlen = abs(float(info_d.get('SVLEN', '0')))
        except ValueError:
            fout.write(line)
            continue

        s_for_coverage = max(0.0, svlen - sva_neutral)
        post = classify(s_for_coverage, lam, nu, cfg)
        klass, pmap = map_class(post)

        # Stash for TSV/summary.
        classifications[vid] = {
            'chrom': chrom, 'pos': int(pos),
            'svtype': info_d.get('SVTYPE', ''),
            'svlen': svlen,
            'class': klass, 'pmap': pmap,
            'lambda': lam, 'nu': nu,
            'has_x_merge': '(x)' in str(info_d.get('repeat_ids', '')),
        }

        # Strict-mode exclusion.
        if strict and (klass == 'X' or pmap < cfg['pmap_min']):
            continue

        # Update INFO.
        info_d['HERVK_CLASS'] = CLASS_LABEL[klass]
        info_d['HERVK_PMAP']  = f'{pmap:.4f}'
        info_d['HERVK_LAMBDA'] = f'{lam:.0f}'
        info_d['HERVK_NU']     = f'{nu:.0f}'
        fields[7] = info_to_str(info_d)

        # Update FORMAT/per-sample with HERVK_AS.
        if samples and len(fields) > 9:
            svtype = info_d.get('SVTYPE', '')
            fmt_keys = fields[8].split(':')
            if 'HERVK_AS' not in fmt_keys:
                fmt_keys.append('HERVK_AS')
                fields[8] = ':'.join(fmt_keys)
                try:
                    gt_idx = fmt_keys.index('GT')
                except ValueError:
                    gt_idx = 0
                for i in range(9, len(fields)):
                    parts = fields[i].split(':') if fields[i] else ['.']
                    gt = parts[gt_idx] if gt_idx < len(parts) else '.'
                    parts.append(build_hervk_as(gt, klass, pmap, svtype, cfg))
                    fields[i] = ':'.join(parts)

        fout.write('\t'.join(fields) + '\n')

    if in_close:
        fin.close()
    if out_close:
        fout.close()


# -------- TSV annotation --------
def annotate_tsv(tsv_in, tsv_out, classifications, strict, cfg):
    """Append HERVK columns to a presence-absence TSV (matched by ID).

    When strict, drop rows whose HERVK class is 'other' or whose pmap is
    below cfg['pmap_min']. Non-candidate rows always pass through.
    """
    new_cols = ['HERVK_class', 'HERVK_pmap', 'HERVK_lambda', 'HERVK_nu']
    with open(tsv_in) as fin, open(tsv_out, 'w') as fout:
        header = fin.readline().rstrip('\n').split('\t')
        try:
            id_idx = header.index('ID')
        except ValueError:
            raise SystemExit(f'TSV {tsv_in} has no ID column')
        fout.write('\t'.join(header + new_cols) + '\n')
        for line in fin:
            row = line.rstrip('\n').split('\t')
            vid = row[id_idx] if id_idx < len(row) else ''
            c = classifications.get(vid)
            if c:
                if strict and (c['class'] == 'X'
                               or c['pmap'] < cfg['pmap_min']):
                    continue
                row += [CLASS_LABEL[c['class']],
                        f"{c['pmap']:.4f}",
                        f"{c['lambda']:.0f}",
                        f"{c['nu']:.0f}"]
            else:
                row += ['NA', 'NA', 'NA', 'NA']
            fout.write('\t'.join(row) + '\n')


# -------- Polyallelic flagging --------
def find_polyallelic(classifications, cfg):
    """Pairs of nearby SVs where one is H_C and the other is H_T or H_B."""
    by_chrom = defaultdict(list)
    for vid, c in classifications.items():
        by_chrom[c['chrom']].append((c['pos'], vid, c))
    flagged = []
    window = cfg['polyallelic_window']
    for chrom, items in by_chrom.items():
        items.sort()
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                if items[j][0] - items[i][0] > window:
                    break
                ka, kb = items[i][2]['class'], items[j][2]['class']
                pair = {ka, kb}
                if 'C' in pair and ('B' in pair or 'T' in pair):
                    flagged.append((chrom, items[i][1], items[j][1],
                                    items[i][2]['class'], items[j][2]['class'],
                                    items[i][0], items[j][0]))
    return flagged


# -------- Summary --------
def write_summary(path, classifications, cfg, polyallelic):
    if not classifications:
        with open(path, 'w') as fh:
            fh.write('# HERV-K polymorphism summary\n\nNo HERV-K candidate '
                     'SVs found.\n')
        return

    classes = [c['class'] for c in classifications.values()]
    counts = Counter(classes)
    total = len(classifications)
    pmaps = sorted(c['pmap'] for c in classifications.values())
    median_pmap = pmaps[len(pmaps) // 2]
    confident = sum(1 for p in pmaps if p >= 0.90)
    ambiguous = total - confident

    def median(xs):
        xs = sorted(xs)
        return xs[len(xs) // 2] if xs else float('nan')

    lines = []
    lines.append('# HERV-K polymorphism summary\n')
    lines.append(f'**Total candidate SVs:** {total}\n')
    lines.append('## Per-class counts\n')
    lines.append('| Class | Label | Count | Median |SVLEN| |')
    lines.append('|---|---|---|---|')
    for k in ('C', 'T', 'B', 'A', 'X'):
        n = counts.get(k, 0)
        med = median([c['svlen'] for c in classifications.values()
                      if c['class'] == k]) if n else float('nan')
        lines.append(f'| H_{k} | {CLASS_LABEL[k]} | {n} | '
                     f'{med:.0f} |' if n else
                     f'| H_{k} | {CLASS_LABEL[k]} | 0 | — |')

    lines.append('\n## Posterior confidence\n')
    lines.append(f'- Confident (pmap >= 0.90): {confident}')
    lines.append(f'- Ambiguous  (pmap <  0.90): {ambiguous}')
    lines.append(f'- Median pmap: {median_pmap:.3f}\n')

    n_xmerge = sum(1 for c in classifications.values() if c['has_x_merge'])
    lines.append('## (x)-merged annotations\n')
    lines.append(f'- Total: {n_xmerge}')
    if n_xmerge:
        xm_by_class = Counter(c['class'] for c in classifications.values()
                              if c['has_x_merge'])
        for k in ('C', 'T', 'B', 'A', 'X'):
            if xm_by_class.get(k):
                lines.append(f'  - H_{k} ({CLASS_LABEL[k]}): '
                             f'{xm_by_class[k]}')
    lines.append('')

    lines.append('## Polyallelic candidates\n')
    if polyallelic:
        lines.append(f'{len(polyallelic)} site(s) flagged within '
                     f'{cfg["polyallelic_window"]} bp:\n')
        lines.append('| chrom | id_a | class_a | pos_a | id_b | class_b | pos_b |')
        lines.append('|---|---|---|---|---|---|---|')
        for chrom, va, vb, ka, kb, pa, pb in polyallelic:
            lines.append(f'| {chrom} | {va} | H_{ka} | {pa} | {vb} | '
                         f'H_{kb} | {pb} |')
    else:
        lines.append('None flagged.')
    lines.append('')

    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


# -------- Self-test on the regression fixture --------
def selftest(tsv_path):
    """Replicates §7 of the plan against the toy fixture (TSV-only)."""
    cfg = load_config(None)
    counts = Counter()
    confident = 0
    with open(tsv_path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        idx = {k: header.index(k) for k in
               ('SVLEN', 'match_lengths', 'repeat_ids',
                'matching_classes', 'n_hits')}
        for line in fh:
            row = line.rstrip('\n').split('\t')
            try:
                svlen = abs(float(row[idx['SVLEN']]))
            except ValueError:
                continue
            lam, nu, sva = parse_hits(
                row[idx['match_lengths']], row[idx['repeat_ids']],
                row[idx['matching_classes']], row[idx['n_hits']], cfg)
            post = classify(max(0.0, svlen - sva), lam, nu, cfg)
            k, p = map_class(post)
            counts[k] += 1
            if p >= 0.90:
                confident += 1
    print('Self-test on', tsv_path)
    for k in 'CTBAX':
        print(f'  H_{k} ({CLASS_LABEL[k]:<14}): {counts.get(k,0)}')
    print(f'  confident (>=0.90): {confident}')


# -------- Main --------
def main():
    ap = argparse.ArgumentParser(
        description='HERV-K (HML-2) SV polymorphism classifier.')
    ap.add_argument('--vcf-in', required=False)
    ap.add_argument('--vcf-out', required=False)
    ap.add_argument('--tsv-in')
    ap.add_argument('--tsv-out')
    ap.add_argument('--summary')
    ap.add_argument('--config')
    ap.add_argument('--strict', action='store_true',
                    help='Drop HERV-K candidate rows whose MAP class is '
                         '"other" or whose pmap < pmap_min.')
    ap.add_argument('--selftest',
                    help='Run regression on the toy fixture TSV and exit.')
    args = ap.parse_args()

    if args.selftest:
        selftest(args.selftest)
        return

    if not args.vcf_in or not args.vcf_out:
        ap.error('--vcf-in and --vcf-out required (unless --selftest)')

    cfg = load_config(args.config)
    classifications = {}
    process_vcf(args.vcf_in, args.vcf_out, cfg, args.strict, classifications)

    if args.tsv_in and args.tsv_out:
        annotate_tsv(args.tsv_in, args.tsv_out, classifications,
                     args.strict, cfg)

    if args.summary:
        polyallelic = find_polyallelic(classifications, cfg)
        write_summary(args.summary, classifications, cfg, polyallelic)

    sys.stderr.write(
        f'[hervk_classify] candidates classified: {len(classifications)}\n')


if __name__ == '__main__':
    main()
