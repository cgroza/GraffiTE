#!/usr/bin/env python3
"""
HERV-K (HML-2) locus layer.

`flag` (stage 3, this phase)
    Group HERV-K candidate records into loci and mark them in place. Nothing is
    merged or dropped here: the discovery VCF is what induces the graph, and
    the human subset is what the analysis reads, so both must keep their record
    structure exactly. The locus grouping is written out as INFO tags plus a
    dedicated table for downstream use.

    Grouping is anchored on the reference HML-2 element when one is known, and
    otherwise on proximity within one LTR length. That radius is derivable
    rather than arbitrary: two records describing the same insertion can be
    placed anywhere inside the shared reference LTR, so they land up to one LTR
    apart. In the CaG cohort the observed offsets are 574 bp (chr11), 559 bp
    (chr6), 107 bp (chr12) and 871 bp (chr1) -- all within one LTR, and the
    first three sit outside truvari's default 500 bp refdist, which is why
    truvari did not collapse them.

`consolidate` (stage 4)
    Not implemented in this phase -- see the phase plan. Genotyping cannot be
    re-run for this cohort, so consolidation is built and validated against a
    supplied giraffe VCF in phase 2.

Usage:
    hervk_reconcile.py flag --vcf-in in.vcf --vcf-out out.vcf \
        --loci-out hervk_loci.tsv [--ref-state refstate.tsv] [--window 1200]
"""

import argparse
import sys
from collections import defaultdict

AUTOSOME_PREFIXES = tuple(f'chr{i}' for i in range(1, 23))

INFO_HEADERS = [
    '##INFO=<ID=HERVK_LOCUS,Number=1,Type=String,Description="HERV-K locus '
    'identifier grouping records that describe the same element.">',
    '##INFO=<ID=HERVK_LOCUS_N,Number=1,Type=Integer,Description="Number of '
    'HERV-K records assigned to this locus.">',
    '##INFO=<ID=HERVK_MERGE_FLAG,Number=0,Type=Flag,Description="This locus '
    'holds more than one record describing the same element. Flagged only -- '
    'records are never merged at this stage, because this VCF must keep its '
    'structure for graph induction and downstream ID matching.">',
    '##INFO=<ID=HERVK_POLARITY_CONFLICT,Number=0,Type=Flag,Description="Records '
    'at this locus imply different REF allele states. The masked reference is '
    'authoritative; consolidation resolves this.">',
]

LOCI_COLUMNS = ['locus_id', 'chrom', 'start', 'end', 'n_records', 'record_ids',
                'n_in_human', 'records_not_in_human', 'ref_state', 'allele_set',
                'per_record_class', 'per_record_evidence', 'per_record_k',
                'arch', 'flags']


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


def read_calls(path):
    """Load every HERV-K candidate from the classifier's call table.

    Clustering must see all candidates, not only the ones in the human subset.
    The --human filter requires FILTER="PASS", and PAV emits TRIM/COMPOUND on
    perfectly real HML-2 records (chr15-2092086-DEL-8221 is one), so grouping
    from the human VCF alone could split a locus whose partner was removed for
    an unrelated reason.
    """
    recs = []
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            if not line.strip():
                continue
            r = dict(zip(header, line.rstrip('\n').split('\t')))
            # Non-HML-2 LTR/ERVK (HERVK9-int, MER11A, LTR13 ...) are reported
            # as `other` by the classifier but must not form loci: they are a
            # different lineage and grouping them would invent merge flags.
            if r['evidence'] == 'NON_HML2':
                continue
            pos, svlen = int(r['pos']), int(r['svlen'])
            recs.append({
                'id': r['id'], 'chrom': r['chrom'], 'pos': pos, 'svlen': svlen,
                'start': pos, 'end': pos if svlen >= 0 else pos + abs(svlen),
                'cls': r['class'], 'evidence': r['evidence'],
                'k': '' if r['k'] == '.' else r['k'],
                'ref_allele': r['allele_ref'], 'alt_allele': r['allele'],
                'ref_state': r['ref_state'], 'arch': r['arch'],
            })
    return recs


def vcf_ids(vcf_path):
    ids = set()
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t', 3)
            if len(f) > 2:
                ids.add(f[2])
    return ids


def cluster(recs, ref_tbl, window, human_ids=None):
    """Assign each record to a locus.

    Two records share a locus when they touch the same reference HML-2 element,
    or when their footprints lie within `window` bp of each other.
    """
    by_chrom = defaultdict(list)
    for r in recs:
        by_chrom[r['chrom']].append(r)

    loci = []
    for chrom in sorted(by_chrom):
        items = sorted(by_chrom[chrom], key=lambda r: (r['start'], r['end']))
        current = []
        cur_end = None
        cur_elems = set()
        for r in items:
            ref = ref_tbl.get(r['id'], {})
            elem = None
            if ref.get('ref_elem_start') and ref.get('ref_elem_end'):
                elem = (ref['ref_elem_chrom'], ref['ref_elem_start'],
                        ref['ref_elem_end'])
            joins = bool(current) and (
                (elem is not None and elem in cur_elems)
                or r['start'] - cur_end <= window)
            if not joins and current:
                loci.append(current)
                current, cur_elems, cur_end = [], set(), None
            current.append(r)
            if elem is not None:
                cur_elems.add(elem)
            cur_end = max(cur_end or r['end'], r['end'])
        if current:
            loci.append(current)

    assignment = {}
    table = []
    for members in loci:
        chrom = members[0]['chrom']
        start = min(m['start'] for m in members)
        end = max(m['end'] for m in members)
        locus_id = f'HERVK_{chrom}_{start}'

        ref_states = {m['ref_allele'] for m in members if m['ref_allele'] not in ('.', '')}
        alleles = set()
        for m in members:
            for a in (m['ref_allele'], m['alt_allele']):
                if a not in ('.', '', 'None'):
                    alleles.add(a)

        in_human = [m['id'] for m in members
                    if human_ids is None or m['id'] in human_ids]
        absent = [m['id'] for m in members if m['id'] not in in_human]

        flags = []
        if len(members) > 1:
            flags.append('MERGE_CANDIDATE')
        # Only interesting when the locus is *split* by the human filter --
        # some members annotated, some not. A locus wholly outside the human
        # subset is simply out of scope, not a problem.
        if absent and in_human:
            flags.append('LOCUS_SPLIT_BY_HUMAN_FILTER')
        if len(ref_states) > 1:
            flags.append('POLARITY_CONFLICT')
        if len(alleles) > 2:
            flags.append('MULTIALLELIC')
        if not chrom.startswith(AUTOSOME_PREFIXES):
            flags.append('PLOIDY_UNVERIFIED')

        for m in members:
            assignment[m['id']] = (locus_id, len(members), flags)

        table.append({
            'locus_id': locus_id, 'chrom': chrom,
            'start': str(start), 'end': str(end),
            'n_records': str(len(members)),
            'record_ids': ','.join(m['id'] for m in members),
            'n_in_human': str(len(in_human)),
            'records_not_in_human': ','.join(absent) or '.',
            'ref_state': ','.join(sorted(ref_states)) or '.',
            'allele_set': ','.join(sorted(alleles)) or '.',
            'per_record_class': ','.join(m['cls'] for m in members),
            'per_record_evidence': ','.join(m['evidence'] for m in members),
            'per_record_k': ','.join(m['k'] or '.' for m in members),
            'arch': ';'.join(m['arch'] or '.' for m in members),
            'flags': ','.join(flags) or '.',
        })
    return assignment, table


def write_vcf(vcf_in, vcf_out, assignment):
    with open(vcf_in) as fin, open(vcf_out, 'w') as fout:
        header_done = False
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
            f = line.rstrip('\n').split('\t')
            if len(f) < 8 or f[2] not in assignment:
                fout.write(line)
                continue
            locus_id, n, flags = assignment[f[2]]
            info = parse_info(f[7])
            info['HERVK_LOCUS'] = locus_id
            info['HERVK_LOCUS_N'] = str(n)
            if 'MERGE_CANDIDATE' in flags:
                info['HERVK_MERGE_FLAG'] = ''
            if 'POLARITY_CONFLICT' in flags:
                info['HERVK_POLARITY_CONFLICT'] = ''
            f[7] = info_to_str(info)
            fout.write('\t'.join(f) + '\n')


def write_loci(table, path):
    with open(path, 'w') as fh:
        fh.write('\t'.join(LOCI_COLUMNS) + '\n')
        for row in table:
            fh.write('\t'.join(row[c] for c in LOCI_COLUMNS) + '\n')


def cmd_flag(args):
    recs = read_calls(args.calls)
    ref_tbl = load_table(args.ref_state)
    human = vcf_ids(args.vcf_in)
    assignment, table = cluster(recs, ref_tbl, args.window, human)
    write_vcf(args.vcf_in, args.vcf_out, assignment)
    write_loci(table, args.loci_out)

    multi = [r for r in table if int(r['n_records']) > 1]
    conflict = [r for r in table if 'POLARITY_CONFLICT' in r['flags']]
    partial = [r for r in table if 'LOCUS_SPLIT_BY_HUMAN_FILTER' in r['flags']]
    sys.stderr.write(
        f'[hervk_reconcile] {len(recs)} candidates -> {len(table)} loci '
        f'({len(human & {r["id"] for r in recs})} annotated in the human VCF); '
        f'{len(multi)} flagged for merge, {len(conflict)} polarity conflicts, '
        f'{len(partial)} loci split by the human filter\n')
    for r in multi:
        sys.stderr.write(f'    {r["locus_id"]}: {r["record_ids"]} [{r["flags"]}]\n')


def cmd_consolidate(args):
    sys.exit(
        'hervk_reconcile consolidate: not implemented in this phase.\n'
        'Stage-4 consolidation (literal multi-allelic records in the human '
        'merged genotypes VCF) is phase 2, and is built against a supplied '
        'giraffe VCF because graph genotyping cannot be re-run for this cohort.')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    f = sub.add_parser('flag', help='group candidates into loci and flag them')
    f.add_argument('--calls', required=True,
                   help='call table from hervk_classify.py --calls-out, '
                        'covering every candidate in the discovery VCF')
    f.add_argument('--vcf-in', required=True,
                   help='human-subset VCF to annotate (the discovery VCF is '
                        'never written to)')
    f.add_argument('--vcf-out', required=True)
    f.add_argument('--loci-out', required=True)
    f.add_argument('--ref-state')
    f.add_argument('--window', type=int, default=1200,
                   help='max footprint gap within one locus (default: one LTR '
                        'plus tolerance)')
    f.set_defaults(func=cmd_flag)

    c = sub.add_parser('consolidate', help='(phase 2) collapse flagged loci')
    c.add_argument('--genotyped-vcf')
    c.add_argument('--loci')
    c.add_argument('--human-vcf')
    c.add_argument('--reference')
    c.add_argument('--genotyper', default='giraffe')
    c.add_argument('--out-vcf')
    c.add_argument('--out-archive')
    c.add_argument('--report')
    c.set_defaults(func=cmd_consolidate)

    args = ap.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
