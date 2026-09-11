#!/usr/bin/env python3
"""Build the PanGenie graph input from pangenome.vcf and track each variant.

prepare: write the one-sample VCF that merge_vcfs.py turns into the graph, and
         a table with one row per ALT allele of pangenome.vcf.
report:  fill the table's in_graph column from merge_vcfs.py output.

merge_vcfs.py raises an error unless each record has one ID per ALT allele
after `bcftools norm -m+` joins the records at a position (#93). norm breaks
that count for two records with the same CHROM/POS/REF/ALT and different IDs
(it keeps one ALT and both IDs), for two records with one ID and different
ALTs (it keeps the ID once), and for a record with ID "." beside another at
its position (it drops the "."). prepare makes one graph variant of records
with the same CHROM/POS/REF/ALT, and replaces an ID that is missing, already
used, or contains a character merge_vcfs.py or PanGenie splits on.

PanGenie writes the graph ID of each allele to INFO/ID of the genotyped VCFs.
To find the genotype of a pangenome.vcf record, look up its graph_ID in the
table and match it against INFO/ID.
"""

import argparse
import sys

RESERVED = set(';,:=| \t')
COLUMNS = ['record', 'CHROM', 'POS', 'pangenome_ID', 'allele', 'graph_ID', 'in_graph', 'note']


def safe_id(vid):
    return vid != '.' and vid != '' and not (set(vid) & RESERVED)


def prepare(vcf_in, vcf_out, table_out):
    used = set()
    graph_id_of = {}    # (CHROM, POS, REF, ALT) -> (graph ID, record number)
    rows = []
    n_dup = n_renamed = 0
    record = 0

    with open(vcf_in) as fin, open(vcf_out, 'w') as fout:
        has_gt = False
        for line in fin:
            if line.startswith('##'):
                if line.startswith('##FORMAT=<ID=GT,'):
                    has_gt = True
                fout.write(line)
                continue
            if line.startswith('#'):
                if not has_gt:
                    fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                fout.write('\t'.join(line.rstrip('\n').split('\t')[:8] + ['FORMAT', 'ref']) + '\n')
                continue

            record += 1
            f = line.rstrip('\n').split('\t')
            chrom, pos, vid, ref, alts = f[0], f[1], f[2], f[3], f[4].split(',')
            # There is no reliable way to pair several IDs (a;b) with ALT
            # alleles, so safe_id() rejects them and the alleles get new IDs.
            for i, alt in enumerate(alts, start=1):
                key = (chrom, pos, ref, alt)
                row = {'record': record, 'CHROM': chrom, 'POS': pos, 'pangenome_ID': vid,
                       'allele': i, 'in_graph': 'no', 'note': '.'}
                if key in graph_id_of:
                    gid, first = graph_id_of[key]
                    row['graph_ID'] = gid
                    row['note'] = 'duplicate_of_record_%d' % first
                    n_dup += 1
                    rows.append(row)
                    continue

                if len(alts) == 1 and safe_id(vid) and vid not in used:
                    gid = vid
                else:
                    gid = 'graffite_rec%d' % record + ('_%d' % i if len(alts) > 1 else '')
                    while gid in used:
                        gid += '_x'
                    row['note'] = 'id_replaced'
                    n_renamed += 1
                used.add(gid)
                graph_id_of[key] = (gid, record)
                row['graph_ID'] = gid
                rows.append(row)
                fout.write('\t'.join([chrom, pos, gid, ref, alt, f[5], f[6], '.', 'GT', '1|0']) + '\n')

    with open(table_out, 'w') as t:
        t.write('\t'.join(COLUMNS) + '\n')
        for r in rows:
            t.write('\t'.join(str(r[c]) for c in COLUMNS) + '\n')

    sys.stderr.write('pangenie_graph_vcf.py: %d records, %d ALT alleles, %d merged into an '
                     'identical earlier allele, %d IDs replaced\n'
                     % (record, len(rows), n_dup, n_renamed))


def report(table, graph_vcf):
    in_graph = set()
    with open(graph_vcf) as g:
        for line in g:
            if line.startswith('#'):
                continue
            info = line.split('\t')[7]
            for field in info.split(';'):
                if field.startswith('ID='):
                    for allele_ids in field[3:].split(','):
                        in_graph.update(allele_ids.split(':'))

    with open(table) as t:
        header = t.readline().rstrip('\n').split('\t')
        rows = [dict(zip(header, l.rstrip('\n').split('\t'))) for l in t]
    missing = 0
    for r in rows:
        r['in_graph'] = 'yes' if r['graph_ID'] in in_graph else 'no'
        missing += r['in_graph'] == 'no'
    with open(table, 'w') as t:
        t.write('\t'.join(header) + '\n')
        for r in rows:
            t.write('\t'.join(r[c] for c in header) + '\n')

    # merge_vcfs.py leaves out overlapping alleles (every ALT sits on the same
    # haplotype of the graph sample, so an overlap drops both), all but two
    # ALT alleles at one position, and alleles with a base other than A, C, G
    # or T.
    sys.stderr.write('pangenie_graph_vcf.py: %d of %d pangenome.vcf alleles are not in the '
                     'PanGenie graph and get no genotype\n' % (missing, len(rows)))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    a = sub.add_parser('prepare')
    a.add_argument('vcf_in', help='sites-only pangenome VCF')
    a.add_argument('vcf_out', help='graph input VCF, one sample with GT 1|0')
    a.add_argument('table', help='tracking table to write')
    b = sub.add_parser('report')
    b.add_argument('table', help='tracking table written by prepare, updated in place')
    b.add_argument('graph_vcf', help='merge_vcfs.py output')
    args = p.parse_args()
    if args.cmd == 'prepare':
        prepare(args.vcf_in, args.vcf_out, args.table)
    else:
        report(args.table, args.graph_vcf)
