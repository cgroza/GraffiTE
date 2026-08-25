#!/usr/bin/env python3
"""
HERV-K (HML-2) allele-state classifier for GraffiTE --human runs.

Version 2 -- evidence first. The previous version decided REF and ALT states
from expected allele-size arithmetic alone (a Gaussian MAP over |SVLEN|, LTR bp
and internal bp against 968 / 8504 / 9472). That is not sufficient: an 8.5 kb
insertion is equally consistent with LTR+INT entering a solo LTR and with a
complete provirus entering an empty site whose internal region carries a
deletion. Size cannot separate them; the reference and the element's own
architecture can.

Resolution order -- the first rule that fires wins, and every call records the
evidence that produced it:

    ARCH_2LTR    two full-length terminal LTRs on the SV allele
                 -> nothing was consumed by the alignment, so REF = null
    ARCH_PERM    one LTR split across the termini, consensus intervals
                 complementary -> the SV sits inside a solo LTR, REF = solo
    ARCH_SOLO    a lone LTR with no internal region
    REF_ANNOT    architecture is degenerate (k = 0, or no terminal LTR);
                 the masked reference window decides
    DENOVO_LTR   reference unavailable/ambiguous; terminal direct repeat scan
    UNRESOLVED   nothing resolved it -- kept and flagged, never dropped

HERVK_PMAP is now a *confidence* derived from the size model, reported and used
by --strict. It never decides the class.

Usage:
    hervk_classify.py \\
        --vcf-in pangenome.human.vcf --vcf-out pangenome.human.hervk.vcf \\
        --arch hervk_arch.tsv --ref-state hervk_refstate.tsv \\
        --tsv-in pangenome.presence-absence_human.tsv \\
        --tsv-out pangenome.presence-absence_human.hervk.tsv \\
        --summary hervk_polymorphism_summary.md
"""

import argparse
import json
import math
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hervk_arch import INT_CONSENSUS_LEN, ltr_len

# -------- Reference architecture --------
LTR_LEN, INT_LEN = 968, INT_CONSENSUS_LEN
SOLO_PROV = LTR_LEN + INT_LEN          # 8504
NULL_PROV = 2 * LTR_LEN + INT_LEN      # 9472

LTR_FAMILY = {'LTR5_Hs', 'LTR5A', 'LTR5B', 'LTR5'}
INT_FAMILY = {'HERVK-int'}
SVA_FAMILIES = {'SVA_A', 'SVA_B', 'SVA_C', 'SVA_D', 'SVA_E', 'SVA_F'}

DEFAULTS = {
    "sigmas": {"solo": 60.0, "prov": 800.0},
    "priors": {"null_solo": 0.55, "solo_prov": 0.20, "null_prov": 0.08,
               "truncated_prov": 0.07, "other": 0.10},
    "t_min": 1500,
    "t_max": 9000,
    "s_range": 30000.0,
    # Below this fraction of the full internal consensus, a proviral allele is
    # reported as truncated_prov rather than solo_prov / null_prov.
    "int_full_frac": 0.80,
    # Strict-mode threshold (only applied when --strict).
    "pmap_min": 0.90,
    # Minimum HML-2 bp for a candidate to be classified at all.
    "min_hml2_bp": 50,
    # Upper |SVLEN| bound for candidacy. The gate is matching_classes=LTR/ERVK
    # with no size limit, which let a 25.3 Mb deletion into the HERV-K set --
    # it carried 38430 RepeatMasker hits and dominated the cost of the whole
    # stage. A whole provirus is 9472 bp; nothing plausible needs 25 kb.
    "max_svlen": 25000,
    # Minimum LTR bp before the SV allele counts as carrying a whole solo LTR
    # (half a consensus). Below this it is an LTR fragment, not an allele.
    "min_solo_bp": 484,
    # Minimum internal bp before the SV allele counts as proviral.
    "min_prov_int_bp": 500,
}

CLASSES = ('null_solo', 'solo_prov', 'truncated_prov', 'null_prov', 'other')


# -------- Config --------
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


# -------- Side-table loading --------
def load_table(path, key='id'):
    if not path:
        return {}
    rows = {}
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        ki = header.index(key)
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < len(header):
                f += [''] * (len(header) - len(f))
            rows[f[ki]] = dict(zip(header, f))
    return rows


# -------- INFO helpers --------
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
    return ';'.join(k if v == '' else f'{k}={v}' for k, v in d.items()) or '.'


def is_candidate(info_d, cfg=None):
    """HERV-K candidate gate: an LTR/ERVK SV, alone or paired with SVA.

    Kept identical to the bcftools --human carve-out in module/main.nf so the
    two gates cannot disagree.
    """
    matching_classes = info_d.get('matching_classes', '')
    if not matching_classes or matching_classes == 'NA':
        return False
    classes = {c.strip() for c in matching_classes.split(',')}
    if 'LTR/ERVK' not in classes:
        return False
    try:
        n_hits = int(float(info_d.get('n_hits', '0')))
    except ValueError:
        return False
    if not (n_hits == 1 or (n_hits == 2 and 'Retroposon/SVA' in classes)):
        return False
    cap = (cfg or DEFAULTS).get('max_svlen')
    try:
        if cap and abs(int(float(info_d.get('_svlen', 0)))) > cap:
            return False
    except (TypeError, ValueError):
        pass
    return True


def hml2_bp_from_info(info_d):
    """Fallback lambda/nu when no architecture table row exists.

    Only used when the raw RepeatMasker table is unavailable for a record;
    it inherits the (x)-collapse problem and so cannot resolve architecture.
    """
    lens, ids = info_d.get('match_lengths', ''), info_d.get('repeat_ids', '')
    if not lens or not ids:
        return 0.0, 0.0
    try:
        lv = [float(x) for x in lens.split(',')]
    except ValueError:
        return 0.0, 0.0
    iv = [r.strip().replace('(x)', '').replace('(VNTR_only)', '')
          for r in ids.split(',')]
    if len(lv) != len(iv):
        return 0.0, 0.0
    lam = sum(v for v, r in zip(lv, iv) if r in LTR_FAMILY)
    nu = sum(v for v, r in zip(lv, iv) if r in INT_FAMILY)
    return lam, nu


# -------- Allele-state resolution --------
def alt_state_from_content(lam, nu, cfg):
    """What kind of HML-2 allele does the SV sequence itself represent?

    Fragments below half an LTR are not an allele of anything -- a 141 bp
    LTR5B piece is a fragment, and calling it a solo LTR would invent a
    polymorphism that is not there.
    """
    if nu >= cfg['min_prov_int_bp']:
        return 'provirus'
    if lam >= cfg['min_solo_bp']:
        return 'solo'
    return None


def classify_pair(ref_state, alt_state, int_bp, cfg):
    """Map a (REF, ALT) allele-state pair onto a HERVK_CLASS label."""
    pair = {ref_state, alt_state}
    if pair == {'null', 'solo'}:
        return 'null_solo'
    if pair == {'null', 'provirus'}:
        return ('truncated_prov'
                if int_bp < cfg['int_full_frac'] * INT_LEN else 'null_prov')
    if pair == {'solo', 'provirus'}:
        return ('truncated_prov'
                if int_bp < cfg['int_full_frac'] * INT_LEN else 'solo_prov')
    return 'other'


def resolve(arch, ref, svlen, lam, nu, cfg):
    """Return (ref_state, alt_state, evidence, notes).

    `arch` is a row from hervk_arch.py; `ref` a row from hervk_ref_state.py.
    Either may be missing.
    """
    notes = []
    sig = (arch or {}).get('signature', '')
    is_ins = svlen >= 0
    observed_ref = (ref or {}).get('ref_state', '')

    def check(arch_ref, alt, code):
        """Architecture resolved it -- but say so if the reference disagrees.

        The two are independent measurements of the same thing, so a
        disagreement is information, not noise. chr6-78894876 reads ARCH_PERM
        (REF = solo) at a locus where the masked reference holds a whole
        provirus; its DEL partner reads the reference correctly. One of the two
        records is a mis-polarised representation of the other, and the
        reference is what settles which.
        """
        if observed_ref in ('null', 'solo', 'provirus') and observed_ref != arch_ref:
            notes.append(f'REF_ARCH_CONFLICT:arch={arch_ref},ref={observed_ref}')
        return arch_ref, alt, code, notes

    # 1-2. Architecture settles it outright -- but which side of the pair the
    # SV *carries* depends on its polarity. For an insertion the architecture
    # describes the allele being added; for a deletion it describes reference
    # sequence being removed, so the same signature means the opposite thing.
    #   ARCH_2LTR  INS: null -> provirus     DEL: provirus -> null
    #   ARCH_PERM  INS: solo -> provirus     DEL: provirus -> solo
    # (The class label is the same either way, since it names the pair, but
    # HERVK_ALLELE_REF / HERVK_ALLELE are what consolidation and the figures
    # read, and those were inverted for deletions.)
    if sig == 'ARCH_2LTR':
        return check('null', 'provirus', 'ARCH_2LTR') if is_ins \
            else check('provirus', 'null', 'ARCH_2LTR')
    if sig == 'ARCH_PERM':
        return check('solo', 'provirus', 'ARCH_PERM') if is_ins \
            else check('provirus', 'solo', 'ARCH_PERM')
    if sig == 'ARCH_SOLO':
        return check('null', 'solo', 'ARCH_SOLO') if is_ins \
            else check('solo', 'null', 'ARCH_SOLO')

    # 4. Degenerate architecture: the reference decides.
    ref_state = observed_ref
    sv_allele = alt_state_from_content(lam, nu, cfg)

    if ref_state in ('null', 'solo', 'provirus', 'partial'):
        if is_ins:
            if ref_state == 'null':
                return 'null', sv_allele, 'REF_ANNOT', notes
            if ref_state == 'solo':
                # A solo LTR plus an inserted LTR+INT block is a provirus.
                return 'solo', 'provirus', 'REF_ANNOT', notes
            notes.append('INS_INTO_NONEMPTY_REF')
            return ref_state, sv_allele, 'REF_ANNOT', notes
        # Deletion: the SV sequence is reference sequence being removed.
        if ref_state == 'provirus':
            # Removing LTR+INT out of LTR-INT-LTR leaves one LTR behind.
            remaining = 'solo' if sv_allele == 'provirus' else 'null'
            return 'provirus', remaining, 'REF_ANNOT', notes
        if ref_state == 'solo':
            return 'solo', 'null', 'REF_ANNOT', notes
        if ref_state == 'partial':
            return 'partial', 'null', 'REF_ANNOT', notes
        notes.append('DEL_FROM_EMPTY_REF')
        return ref_state, None, 'REF_ANNOT', notes

    # 5-6. Nothing resolved it. Keep the record, say so plainly.
    return None, sv_allele, 'UNRESOLVED', notes


# -------- Confidence (never decides the class) --------
def _gauss(x, mu, sigma):
    return math.exp(-0.5 * ((x - mu) / sigma) ** 2) / (sigma * math.sqrt(2 * math.pi))


def size_confidence(svlen, resolved_class, cfg):
    """Posterior of the resolved class under a properly normalised size model.

    Every component is a real probability density here -- the v1 model mixed
    unnormalised Gaussians of different sigma with normalised flat densities,
    which silently favoured the narrow solo-LTR hypothesis by ~26x.
    """
    s = abs(svlen)
    sig, pri = cfg['sigmas'], cfg['priors']
    t_lo, t_hi = cfg['t_min'], cfg['t_max']
    lik = {
        'null_solo': _gauss(s, LTR_LEN, sig['solo']),
        'solo_prov': _gauss(s, SOLO_PROV, sig['prov']),
        'null_prov': _gauss(s, NULL_PROV, sig['prov']),
        'truncated_prov': (1.0 / (t_hi - t_lo)) if t_lo <= s <= t_hi else 0.0,
        'other': 1.0 / cfg['s_range'],
    }
    post = {k: lik[k] * pri.get(k, 0.0) for k in lik}
    z = sum(post.values())
    if z <= 0:
        return 0.0
    return post.get(resolved_class, post['other']) / z


# -------- VCF I/O --------
INFO_HEADERS = [
    '##INFO=<ID=HERVK_CLASS,Number=1,Type=String,Description="HERV-K (HML-2) '
    'polymorphism class: null_solo|solo_prov|truncated_prov|null_prov|other.">',
    '##INFO=<ID=HERVK_ALLELE_REF,Number=1,Type=String,Description="HERV-K state '
    'of the REF allele: null|solo|provirus|partial|. (unresolved).">',
    '##INFO=<ID=HERVK_ALLELE,Number=.,Type=String,Description="HERV-K state of '
    'each ALT allele, in ALT order: null|solo|provirus|partial|.">',
    '##INFO=<ID=HERVK_EVIDENCE,Number=1,Type=String,Description="Evidence that '
    'resolved the allele states: ARCH_2LTR|ARCH_PERM|ARCH_SOLO|REF_ANNOT|'
    'DENOVO_LTR|UNRESOLVED|NON_HML2.">',
    '##INFO=<ID=HERVK_ARCH,Number=1,Type=String,Description="Element 5-prime to '
    '3-prime architecture of the SV allele with consensus intervals, e.g. '
    'LTR:575-968/INT:1-7536/LTR:1-574.">',
    '##INFO=<ID=HERVK_K,Number=1,Type=Integer,Description="LTR permutation '
    'point: alignment breakpoint inside the reference solo LTR. An alignment '
    'property, not a biological one; do not key on its value.">',
    '##INFO=<ID=HERVK_REF_STATE,Number=1,Type=String,Description="HML-2 state '
    'of the masked reference window: null|solo|provirus|partial|unknown.">',
    '##INFO=<ID=HERVK_LAMBDA,Number=1,Type=Float,Description="bp of HML-2 LTR '
    'sequence on the SV allele, from the tiled RepeatMasker fragments.">',
    '##INFO=<ID=HERVK_NU,Number=1,Type=Float,Description="bp of HML-2 internal '
    '(HERVK-int) sequence on the SV allele.">',
    '##INFO=<ID=HERVK_COV,Number=1,Type=Float,Description="Fraction of the SV '
    'allele that is HML-2 sequence.">',
    '##INFO=<ID=HERVK_PMAP,Number=1,Type=Float,Description="Confidence in the '
    'resolved class under the size model. Reporting only -- it does not '
    'determine the class.">',
    '##INFO=<ID=HERVK_NOTE,Number=.,Type=String,Description="Diagnostics for '
    'this call, e.g. REF_ARCH_CONFLICT when the architecture and the masked '
    'reference imply different REF allele states.">',
]

TSV_COLUMNS = ['HERVK_class', 'HERVK_allele_ref', 'HERVK_allele',
               'HERVK_evidence', 'HERVK_k', 'HERVK_ref_state',
               'HERVK_lambda', 'HERVK_nu', 'HERVK_pmap', 'HERVK_arch']

# Per-candidate call table. Written for *every* HERV-K candidate in the input,
# including records the --human filter will later drop (a non-PASS FILTER from
# the SV caller is enough to remove one: chr15-2092086-DEL-8221 carries TRIM).
# The locus layer groups from this table rather than from the human VCF, so a
# merge partner cannot go missing just because it failed an unrelated filter.
CALLS_COLUMNS = ['id', 'chrom', 'pos', 'svlen', 'class', 'allele_ref',
                 'allele', 'evidence', 'k', 'ref_state', 'lambda', 'nu',
                 'cov', 'pmap', 'arch', 'notes']


def write_calls(results, path):
    with open(path, 'w') as fh:
        fh.write('\t'.join(CALLS_COLUMNS) + '\n')
        for vid in sorted(results):
            r = results[vid]
            fh.write('\t'.join([
                vid, r['chrom'], str(r['pos']), str(r['svlen']), r['cls'],
                r['ref_allele'] or '.', r['alt_allele'] or '.', r['evidence'],
                str(r['k']) if r['k'] not in ('', None) else '.',
                r['ref_state'] or 'unknown',
                f"{r['lambda']:.0f}", f"{r['nu']:.0f}",
                f"{r['cov']:.4f}", f"{r['pmap']:.4f}", r['arch'] or '.',
                ','.join(r.get('notes') or []) or '.',
            ]) + '\n')


def classify_record(vid, info_d, svlen, arch_tbl, ref_tbl, cfg):
    """Resolve one candidate SV. Returns a result dict."""
    arch = arch_tbl.get(vid)
    ref = ref_tbl.get(vid)

    if arch:
        lam = float(arch.get('ltr_bp') or 0)
        nu = float(arch.get('int_bp') or 0)
    else:
        lam, nu = hml2_bp_from_info(info_d)

    result = {
        'chrom': info_d.get('_chrom', ''), 'pos': info_d.get('_pos', 0),
        'svlen': svlen, 'lambda': lam, 'nu': nu,
        'arch': (arch or {}).get('arch', ''),
        'k': (arch or {}).get('k', ''),
        'ref_state': (ref or {}).get('ref_state', 'unknown'),
        'notes': [],
    }

    # Non-HML-2 LTR/ERVK (HERVK9-int, MER11A, LTR13 ...): say so explicitly.
    # v1 skipped these silently, which also let them slip through --strict.
    if lam + nu < cfg['min_hml2_bp']:
        result.update(cls='other', ref_allele=None, alt_allele=None,
                      evidence='NON_HML2', pmap=0.0, cov=0.0)
        return result

    ref_state, alt_state, evidence, notes = resolve(arch, ref, svlen, lam, nu, cfg)
    cls = classify_pair(ref_state, alt_state, nu, cfg) \
        if (ref_state and alt_state) else 'other'

    result.update(cls=cls, ref_allele=ref_state, alt_allele=alt_state,
                  evidence=evidence, notes=notes,
                  pmap=size_confidence(svlen, cls, cfg),
                  cov=((lam + nu) / abs(svlen)) if svlen else 0.0)
    return result


def process_vcf(vcf_in, vcf_out, cfg, strict, arch_tbl, ref_tbl, results):
    fin = open(vcf_in) if vcf_in != '-' else sys.stdin
    if vcf_out is None:
        fout = open(os.devnull, 'w')
    else:
        fout = open(vcf_out, 'w') if vcf_out != '-' else sys.stdout
    header_done = False
    try:
        for line in fin:
            if line.startswith('##'):
                fout.write(line)
                continue
            if line.startswith('#CHROM'):
                for h in INFO_HEADERS:
                    fout.write(h + '\n')
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

            info_d = parse_info(fields[7])
            svlen = len(fields[4].split(',')[0]) - len(fields[3])
            info_d['_svlen'] = str(svlen)
            if not is_candidate(info_d, cfg):
                info_d.pop('_svlen', None)
                fout.write(line)
                continue
            info_d.pop('_svlen', None)

            vid = fields[2]
            info_d['_chrom'], info_d['_pos'] = fields[0], int(fields[1])
            r = classify_record(vid, info_d, svlen, arch_tbl, ref_tbl, cfg)
            r['chrom'], r['pos'] = fields[0], int(fields[1])
            results[vid] = r

            if strict and (r['cls'] == 'other' or r['pmap'] < cfg['pmap_min']):
                continue

            for key in ('_chrom', '_pos'):
                info_d.pop(key, None)
            info_d['HERVK_CLASS'] = r['cls']
            info_d['HERVK_ALLELE_REF'] = r['ref_allele'] or '.'
            info_d['HERVK_ALLELE'] = r['alt_allele'] or '.'
            info_d['HERVK_EVIDENCE'] = r['evidence']
            if r['arch']:
                info_d['HERVK_ARCH'] = r['arch']
            if r['k'] not in ('', None):
                info_d['HERVK_K'] = str(r['k'])
            info_d['HERVK_REF_STATE'] = r['ref_state'] or 'unknown'
            info_d['HERVK_LAMBDA'] = f"{r['lambda']:.0f}"
            info_d['HERVK_NU'] = f"{r['nu']:.0f}"
            info_d['HERVK_COV'] = f"{r['cov']:.4f}"
            info_d['HERVK_PMAP'] = f"{r['pmap']:.4f}"
            if r.get('notes'):
                info_d['HERVK_NOTE'] = ','.join(r['notes'])
            fields[7] = info_to_str(info_d)
            fout.write('\t'.join(fields) + '\n')
    finally:
        if vcf_in != '-':
            fin.close()
        if vcf_out != '-':
            fout.close()



def annotate_tsv(tsv_in, tsv_out, results, strict, cfg):
    with open(tsv_in) as fin, open(tsv_out, 'w') as fout:
        header = fin.readline().rstrip('\n').split('\t')
        try:
            id_idx = header.index('ID')
        except ValueError:
            raise SystemExit(f'TSV {tsv_in} has no ID column')
        fout.write('\t'.join(header + TSV_COLUMNS) + '\n')
        for line in fin:
            row = line.rstrip('\n').split('\t')
            vid = row[id_idx] if id_idx < len(row) else ''
            r = results.get(vid)
            if r:
                if strict and (r['cls'] == 'other' or r['pmap'] < cfg['pmap_min']):
                    continue
                row += [r['cls'], r['ref_allele'] or 'NA', r['alt_allele'] or 'NA',
                        r['evidence'], str(r['k']) if r['k'] not in ('', None) else 'NA',
                        r['ref_state'] or 'unknown',
                        f"{r['lambda']:.0f}", f"{r['nu']:.0f}",
                        f"{r['pmap']:.4f}", r['arch'] or 'NA']
            else:
                row += ['NA'] * len(TSV_COLUMNS)
            fout.write('\t'.join(row) + '\n')


# -------- Summary --------
def write_summary(path, results, cfg):
    with open(path, 'w') as fh:
        fh.write('# HERV-K (HML-2) polymorphism summary\n\n')
        if not results:
            fh.write('No HERV-K candidate SVs found.\n')
            return
        fh.write(f'Candidates classified: **{len(results)}** '
                 '(human subset of the pangenome VCF).\n\n')

        fh.write('## Classes\n\n| Class | Count | Median abs(SVLEN) |\n')
        fh.write('|---|---|---|\n')
        by_cls = defaultdict(list)
        for r in results.values():
            by_cls[r['cls']].append(abs(r['svlen']))
        for cls in CLASSES:
            sizes = sorted(by_cls.get(cls, []))
            med = sizes[len(sizes) // 2] if sizes else '--'
            fh.write(f'| {cls} | {len(sizes)} | {med} |\n')

        fh.write('\n## Evidence\n\n| Evidence | Count |\n|---|---|\n')
        ev = defaultdict(int)
        for r in results.values():
            ev[r['evidence']] += 1
        for k in sorted(ev, key=lambda x: -ev[x]):
            fh.write(f'| {k} | {ev[k]} |\n')

        perm = sorted((int(r['k']), vid) for vid, r in results.items()
                      if r['evidence'] == 'ARCH_PERM' and r['k'] not in ('', None))
        fh.write(f'\n## LTR permutation points ({len(perm)} records)\n\n')
        if perm:
            fh.write('`k` is where the aligner broke the reference solo LTR. It is an '
                     'alignment property, varies between haplotypes and callers, and '
                     'nothing downstream may key on it.\n\n| k | Record |\n|---|---|\n')
            for k, vid in perm:
                fh.write(f'| {k} | {vid} |\n')
        else:
            fh.write('None.\n')

        unres = [vid for vid, r in results.items()
                 if r['evidence'] in ('UNRESOLVED',) or r['cls'] == 'other']
        fh.write(f'\n## Unresolved / other ({len(unres)})\n\n')
        fh.write('These are retained, never dropped: a locus that cannot be resolved '
                 'is a locus to look at, not one to discard.\n\n')
        for vid in sorted(unres):
            r = results[vid]
            fh.write(f'- `{vid}` — {r["cls"]}, {r["evidence"]}, '
                     f'ref_state={r["ref_state"]}\n')


# -------- Self-test --------
def selftest(arch_path, ref_path, expect_path):
    """Assert the architecture + reference tables reproduce known calls."""
    cfg = load_config(None)
    arch_tbl, ref_tbl = load_table(arch_path), load_table(ref_path)
    failures, checked = [], 0
    with open(expect_path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            row = dict(zip(header, line.rstrip('\n').split('\t')))
            vid = row['id']
            info_d = {'matching_classes': row.get('matching_classes', 'LTR/ERVK'),
                      'n_hits': row.get('n_hits', '1')}
            r = classify_record(vid, info_d, int(row['svlen']),
                                arch_tbl, ref_tbl, cfg)
            checked += 1
            for field, got in (('class', r['cls']), ('evidence', r['evidence']),
                               ('allele_ref', r['ref_allele'] or ''),
                               ('allele', r['alt_allele'] or '')):
                want = row.get(field, '')
                if want and want != got:
                    failures.append(f'{vid}: {field} expected {want}, got {got}')
            if row.get('k'):
                got_k = str(r['k'])
                if got_k != row['k']:
                    failures.append(f"{vid}: k expected {row['k']}, got {got_k}")

    print(f'hervk_classify selftest: {checked} records checked')
    for f in failures:
        print('  FAIL', f)
    if failures:
        sys.exit(1)
    print('  all expectations met')


def main():
    ap = argparse.ArgumentParser(
        description='HERV-K (HML-2) allele-state classifier (evidence-first).')
    ap.add_argument('--vcf-in')
    ap.add_argument('--vcf-out',
                    help='annotated VCF; omit to classify without writing one')
    ap.add_argument('--calls-out',
                    help='per-candidate call table for the locus layer')
    ap.add_argument('--arch', help='architecture TSV from hervk_arch.py')
    ap.add_argument('--ref-state', help='reference-state TSV from hervk_ref_state.py')
    ap.add_argument('--tsv-in')
    ap.add_argument('--tsv-out')
    ap.add_argument('--summary')
    ap.add_argument('--config')
    ap.add_argument('--max-svlen', type=int,
                    help='|SVLEN| cap for candidacy; must match the cap used to '
                         'build the candidate list the reference masking ran on')
    ap.add_argument('--strict', action='store_true',
                    help='drop candidates classed "other" or below pmap_min '
                         '(off by default: dropping records is what hid the '
                         'null_prov failure in v1)')
    ap.add_argument('--selftest', nargs=3,
                    metavar=('ARCH_TSV', 'REFSTATE_TSV', 'EXPECT_TSV'))
    args = ap.parse_args()

    if args.selftest:
        selftest(*args.selftest)
        return

    if not args.vcf_in:
        ap.error('--vcf-in is required (unless --selftest)')
    if not args.vcf_out and not args.calls_out:
        ap.error('at least one of --vcf-out or --calls-out is required')

    cfg = load_config(args.config)
    if args.max_svlen:
        cfg['max_svlen'] = args.max_svlen
    arch_tbl = load_table(args.arch)
    ref_tbl = load_table(args.ref_state)
    results = {}
    process_vcf(args.vcf_in, args.vcf_out, cfg, args.strict,
                arch_tbl, ref_tbl, results)

    if args.calls_out:
        write_calls(results, args.calls_out)
    if args.tsv_in and args.tsv_out:
        annotate_tsv(args.tsv_in, args.tsv_out, results, args.strict, cfg)
    if args.summary:
        write_summary(args.summary, results, cfg)

    counts = defaultdict(int)
    for r in results.values():
        counts[r['cls']] += 1
    detail = ', '.join(f'{k}={counts[k]}' for k in CLASSES if counts[k])
    sys.stderr.write(f'[hervk_classify] {len(results)} candidates ({detail})\n')


if __name__ == '__main__':
    main()
