#!/usr/bin/env python3
"""Say what became of every pangenome.vcf allele once the graph had genotyped it.

A record can leave the callset at four points between pangenome.vcf and
GraffiTE.merged.genotypes.vcf.gz, and until now only the first two were written
down anywhere, for one back end:

  not_in_graph          the allele never entered the graph. On the PanGenie path
                        merge_vcfs.py leaves out alleles that overlap another at
                        the same site, all but two ALTs at one position, and
                        anything holding a base other than A, C, G or T
                        (see issue #102). pangenie_graph_variants.tsv marks
                        these in_graph=no.
  duplicate_of_record_N two pangenome.vcf records with the same CHROM, POS, REF
                        and ALT became one graph allele, so one genotyped record
                        stands for both and only one of the two IDs survives.
  no_match_in_merge     merge_VCFs transfers the annotation with
                        `bcftools annotate -a pangenome.vcf`, which matches on
                        CHROM, POS, REF and a shared ALT. A record the genotyper
                        moved or rewrote matches nothing and arrives with no ID
                        and no annotation.
  genotyped             the record is in the merged VCF under its own ID.

Nothing equivalent exists for the vg back ends, and vg construct names its alt
paths from a hash of the coordinates and alleles rather than from the VCF ID, so
position is the only link back. The two back ends therefore need different join
keys, and the report names which one it used:

  pangenie   INFO/ID of the genotyped records, via the graph_ID column of
             pangenie_graph_variants.tsv. Number=A, so it survives
             `bcftools norm -m-` onto each split record, and it does not depend
             on the position holding still.
  giraffe    the ID column after merge_VCFs' annotate step, which is the
             pangenome.vcf ID wherever the match succeeded.

Counting records in the two files does not answer this. The merged VCF is a
superset of pangenome.vcf on the vg path, because `vg call -a -A` calls every
snarl including ones that were never GraffiTE variants, and a subset on the
PanGenie path.

Usage:
  genotyping_audit.py --pangenome pangenome.vcf \
                      --merged GraffiTE.merged.genotypes.vcf.gz \
                      [--graph-table pangenie_graph_variants.tsv] \
                      [--trusted-ids trusted.ids] \
                      -o genotyping_record_audit.tsv
"""

import argparse
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hervk_reconcile import detect_genotyper  # noqa: E402

COLUMNS = ['pangenome_ID', 'CHROM', 'POS', 'SVTYPE', 'SVLEN', 'allele', 'trusted',
           'in_graph', 'graph_note', 'graph_ID', 'genotyped', 'annotated',
           'n_samples_genotyped', 'lost_at']

MISSING_GT = {'.', './.', '.|.', '', None}


def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def info_get(info, key):
    for field in info.split(';'):
        if field.startswith(key + '='):
            return field[len(key) + 1:]
    return 'NA'


def read_pangenome(path):
    """One row per ALT allele, in file order."""
    rows = []
    with opener(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            for i, _alt in enumerate(f[4].split(','), start=1):
                rows.append({'pangenome_ID': f[2], 'CHROM': f[0], 'POS': f[1],
                             'SVTYPE': info_get(f[7], 'SVTYPE'),
                             'SVLEN': info_get(f[7], 'SVLEN'),
                             'allele': i})
    return rows


def read_graph_table(path):
    """(CHROM, POS, pangenome_ID, allele) -> the prepare/report row."""
    table = {}
    with opener(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            r = dict(zip(header, line.rstrip('\n').split('\t')))
            table[(r['CHROM'], r['POS'], r['pangenome_ID'], r['allele'])] = r
    return table


def read_merged(path):
    """IDs, INFO/ID values, and how many samples carry a genotype for each."""
    ids, graph_ids, called = set(), set(), {}
    with opener(path) as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            vid = f[2]
            n = 0
            if len(f) > 9:
                for sample in f[9:]:
                    gt = sample.split(':')[0]
                    if gt not in MISSING_GT:
                        n += 1
            if vid != '.':
                ids.add(vid)
                called[vid] = max(called.get(vid, 0), n)
            gid = info_get(f[7], 'ID')
            if gid != 'NA':
                for one in gid.replace(':', ',').split(','):
                    if one:
                        graph_ids.add(one)
                        called[one] = max(called.get(one, 0), n)
    return ids, graph_ids, called


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--pangenome', required=True)
    ap.add_argument('--merged', required=True)
    ap.add_argument('--graph-table', help='pangenie_graph_variants.tsv; pangenie runs only')
    ap.add_argument('--trusted-ids', help='one ID per line, as trusted_genotypes writes')
    ap.add_argument('-o', '--output', default='-')
    args = ap.parse_args()

    backend = detect_genotyper(args.merged) or 'unknown'
    # The graph table is what makes the INFO/ID join possible; without it even a
    # PanGenie file has to fall back on the ID column.
    key = 'INFO/ID' if (backend == 'pangenie' and args.graph_table) else 'ID'

    rows = read_pangenome(args.pangenome)
    graph = read_graph_table(args.graph_table) if args.graph_table else {}
    merged_ids, merged_graph_ids, called = read_merged(args.merged)
    trusted = set()
    if args.trusted_ids:
        with open(args.trusted_ids) as fh:
            trusted = {l.strip() for l in fh if l.strip()}

    counts = {'genotyped': 0, 'not_in_graph': 0, 'duplicate': 0,
              'no_match_in_merge': 0, 'not_genotyped': 0}
    out = []
    for r in rows:
        g = graph.get((r['CHROM'], r['POS'], r['pangenome_ID'], str(r['allele'])), {})
        in_graph = g.get('in_graph', 'NA')
        note = g.get('note', 'NA')
        graph_id = g.get('graph_ID', 'NA')

        # Two different questions. INFO/ID says whether the graph produced a
        # call; the ID column says whether merge_VCFs' annotate found the record
        # and transferred the annotation onto it. A record can be genotyped and
        # still arrive anonymous, which is the case worth surfacing.
        annotated = r['pangenome_ID'] in merged_ids
        if key == 'INFO/ID':
            genotyped = graph_id in merged_graph_ids
            n_called = called.get(graph_id, 0)
        else:
            genotyped = annotated   # no INFO/ID carrier: the two cannot be told apart
            n_called = called.get(r['pangenome_ID'], 0)

        if in_graph == 'no':
            lost = 'not_in_graph'
            counts['not_in_graph'] += 1
        elif note.startswith('duplicate_of_record_'):
            # The allele is in the graph and genotyped; the record that shares it
            # carries the call, and this ID does not appear on its own.
            lost = note
            counts['duplicate'] += 1
        elif genotyped and annotated:
            lost = 'genotyped'
            counts['genotyped'] += 1
        elif genotyped:
            lost = 'no_match_in_merge'
            counts['no_match_in_merge'] += 1
        else:
            lost = 'not_genotyped'
            counts['not_genotyped'] += 1

        r.update({'trusted': 'yes' if r['pangenome_ID'] in trusted else
                             ('no' if trusted else 'NA'),
                  'in_graph': in_graph, 'graph_note': note, 'graph_ID': graph_id,
                  'genotyped': ('NA' if key == 'ID' else ('yes' if genotyped else 'no')),
                  'annotated': 'yes' if annotated else 'no',
                  'n_samples_genotyped': n_called, 'lost_at': lost})
        out.append(r)

    summary = [
        '# genotyper: %s (join key: %s)' % (backend, key),
        '# pangenome.vcf alleles: %d' % len(out),
        '# genotyped: %d' % counts['genotyped'],
        '# never entered the graph: %d' % counts['not_in_graph'],
        '# merged into a duplicate record: %d' % counts['duplicate'],
        '# genotyped but unmatched by merge_VCFs: %d' % counts['no_match_in_merge'],
        '# absent from the merged VCF: %d' % counts['not_genotyped'],
    ]
    if key == 'ID':
        summary.append('# joined on the ID column, so a record the graph '
                       'genotyped and the merge then failed to match counts as '
                       'absent; only the PanGenie path carries INFO/ID')
    if not args.graph_table:
        summary.append('# no graph table: "never entered the graph" is not '
                       'separable from the rest on this back end')

    fh = sys.stdout if args.output == '-' else open(args.output, 'w')
    for line in summary:
        fh.write(line + '\n')
    fh.write('\t'.join(COLUMNS) + '\n')
    for r in out:
        fh.write('\t'.join(str(r[c]) for c in COLUMNS) + '\n')
    if fh is not sys.stdout:
        fh.close()
    sys.stderr.write('\n'.join(s[2:] for s in summary) + '\n')


if __name__ == '__main__':
    main()
