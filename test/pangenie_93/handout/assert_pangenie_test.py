#!/usr/bin/env python3
"""Assertions for the #93 PanGenie test. Standard library only.

FAIL lines decide the result (exit 1 on any). NOTE lines are measurements to
report back verbatim; they are not pass/fail.
"""
import argparse
import collections
import gzip
import os
import sys

# Measured on the #93 pangenome.vcf before this test was written, with
# merge_vcfs.py run against a stand-in reference (every base A): 1,550 alleles
# with N in ALT and 19 in overlapping clusters. The real reference is expected
# to give the same number, but that is a prediction, so a mismatch is a NOTE.
ISSUE93 = {'records': 51910, 'duplicates': 31, 'not_in_graph': 1569}
RESERVED = set(';,:=| \t')

fails = 0


def ok(msg):
    print('  [ ok ] ' + msg)


def fail(msg):
    global fails
    fails += 1
    print('  [FAIL] ' + msg)


def note(msg):
    print('  NOTE: ' + msg)


def open_any(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def read_vcf(path):
    """Return (sample names, list of (chrom, pos, id, ref, alt, info, gts))."""
    samples, recs = [], []
    with open_any(path) as f:
        for line in f:
            if line.startswith('##'):
                continue
            fields = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                samples = fields[9:]
                continue
            gts = []
            if len(fields) > 9:
                fmt = fields[8].split(':')
                gi = fmt.index('GT') if 'GT' in fmt else None
                gts = [s.split(':')[gi] if gi is not None else '.' for s in fields[9:]]
            recs.append((fields[0], fields[1], fields[2], fields[3], fields[4], fields[7], gts))
    return samples, recs


def info_ids(info):
    for field in info.split(';'):
        if field.startswith('ID='):
            out = []
            for allele in field[3:].split(','):
                out.extend(allele.split(':'))
            return out
    return []


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--outdir', required=True)
    p.add_argument('--graffite-vcf', required=True)
    p.add_argument('--samples', required=True, help='comma-separated')
    p.add_argument('--distinct-reads', action='store_true',
                   help='the samples have different reads (skips the same-reads concordance note)')
    a = p.parse_args()
    samples = a.samples.split(',')
    gdir = os.path.join(a.outdir, '4_Genotyping')

    print('== outputs ==')
    table_path = os.path.join(gdir, 'pangenie_graph_variants.tsv')
    merged_path = os.path.join(gdir, 'GraffiTE.merged.genotypes.vcf.gz')
    per_sample = {s: os.path.join(gdir, s + '_genotyping.vcf.gz') for s in samples}
    missing = [x for x in [table_path, merged_path] + list(per_sample.values()) if not os.path.exists(x)]
    for x in [table_path, merged_path] + list(per_sample.values()):
        (ok if os.path.exists(x) else fail)(('found ' if os.path.exists(x) else 'missing ') + x)
    if missing:
        print('\nstopping: outputs missing (did the run finish? see nextflow_trace.txt)')
        return 1

    print('== tracking table vs pangenome.vcf ==')
    _, pan = read_vcf(a.graffite_vcf)
    alleles, first_of, dup_pairs = [], {}, []
    for n, (c, pos, vid, ref, alts, _, _) in enumerate(pan, start=1):
        for i, alt in enumerate(alts.split(','), start=1):
            key = (c, pos, ref, alt)
            if key in first_of:
                dup_pairs.append(((str(n), str(i)), first_of[key]))
            else:
                first_of[key] = (str(n), str(i))
            alleles.append((n, vid, key))
    dup_expect = len(dup_pairs)
    with open(table_path) as t:
        header = t.readline().rstrip('\n').split('\t')
        rows = [dict(zip(header, l.rstrip('\n').split('\t'))) for l in t]
    want_cols = ['record', 'CHROM', 'POS', 'pangenome_ID', 'allele', 'graph_ID', 'in_graph', 'note']
    (ok if header == want_cols else fail)('table columns: ' + ','.join(header))
    (ok if len(rows) == len(alleles) else fail)(
        '%d table rows for %d ALT alleles in pangenome.vcf' % (len(rows), len(alleles)))

    dups = [r for r in rows if r['note'].startswith('duplicate_of_record_')]
    (ok if len(dups) == dup_expect else fail)(
        '%d duplicate rows, %d duplicate alleles counted in pangenome.vcf' % (len(dups), dup_expect))
    by_rec = {(r['record'], r['allele']): r for r in rows}
    bad_dup = [d for d, first in dup_pairs
               if by_rec.get(d, {}).get('graph_ID') is None
               or by_rec[d]['graph_ID'] != by_rec.get(first, {}).get('graph_ID')]
    (ok if not bad_dup else fail)(
        'every duplicate allele shares the graph_ID of its first occurrence (%d do not)' % len(bad_dup))

    replaced = [r for r in rows if r['note'] == 'id_replaced']
    note('%d IDs replaced (0 expected for the #93 file, which has no missing or repeated IDs)' % len(replaced))
    bad_ids = [r['graph_ID'] for r in rows if set(r['graph_ID']) & RESERVED or r['graph_ID'] in ('', '.')]
    (ok if not bad_ids else fail)('no graph_ID contains a reserved character or is empty (%d do)' % len(bad_ids))

    yes = {r['graph_ID'] for r in rows if r['in_graph'] == 'yes'}
    no_rows = [r for r in rows if r['in_graph'] == 'no']
    note('%d of %d alleles in_graph=no; %d distinct graph IDs in the graph'
         % (len(no_rows), len(rows), len(yes)))
    if len(pan) == ISSUE93['records'] and dup_expect == ISSUE93['duplicates']:
        msg = 'in_graph=no is %d; predicted %d from the stand-in reference run' % (len(no_rows), ISSUE93['not_in_graph'])
        note(msg + (' (matches)' if len(no_rows) == ISSUE93['not_in_graph'] else ' (DIFFERS -- report this)'))

    print('== every sample genotyped (the "1 of 6 samples" half of #93) ==')
    msamples, merged = read_vcf(merged_path)
    for s in samples:
        (ok if s in msamples else fail)('%s is a sample column of the merged VCF' % s)
    (ok if len(msamples) == len(samples) else fail)(
        'merged VCF has %d sample columns, %d samples in reads.csv' % (len(msamples), len(samples)))

    print('== graph IDs reach the genotyped VCFs ==')
    for s, path in per_sample.items():
        _, recs = read_vcf(path)
        got = set()
        for r in recs:
            got.update(info_ids(r[5]))
        (ok if got == yes else fail)(
            '%s: INFO/ID holds %d graph IDs, table has %d in_graph=yes (%d only in VCF, %d only in table)'
            % (s, len(got), len(yes), len(got - yes), len(yes - got)))
    got = set()
    for r in merged:
        got.update(info_ids(r[5]))
    (ok if got == yes else fail)(
        'merged VCF: INFO/ID holds %d graph IDs, table has %d in_graph=yes (%d only in VCF, %d only in table)'
        % (len(got), len(yes), len(got - yes), len(yes - got)))
    dup_gids = {r['graph_ID'] for r in dups if r['in_graph'] == 'yes'}
    (ok if dup_gids <= got else fail)(
        '%d of %d duplicate groups in the graph have a genotyped record' % (len(dup_gids & got), len(dup_gids)))

    print('== measurements ==')
    pan_ids_of = collections.defaultdict(set)
    for r in rows:
        pan_ids_of[r['graph_ID']].add(r['pangenome_ID'])
    match = sum(1 for r in merged if any(r[2] in pan_ids_of.get(g, ()) for g in info_ids(r[5])))
    note('merged VCF: %d records; %d have an ID column that is one of the pangenome IDs of their graph_ID'
         % (len(merged), match))
    for i, s in enumerate(msamples):
        c = collections.Counter(r[6][i].replace('|', '/') for r in merged)
        note('%s genotypes: %s' % (s, ', '.join('%s=%d' % kv for kv in c.most_common())))
    if len(msamples) == 2 and not a.distinct_reads:
        same = sum(1 for r in merged if r[6][0].replace('|', '/') == r[6][1].replace('|', '/'))
        note('same reads under two names: %d of %d genotypes identical' % (same, len(merged)))

    print()
    print('RESULT: %s (%d failure%s)' % ('PASS' if fails == 0 else 'FAIL', fails, '' if fails == 1 else 's'))
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
