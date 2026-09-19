#!/usr/bin/env python3
"""Append the records the genotyping audit's lost_at categories need.

run_matrix.sh's `pangenie` cell feeds the spine's pangenome.vcf in unaltered, so
every record is well formed and genotyping_record_audit.tsv comes back with zero
not_in_graph and zero duplicates. A working audit and a no-op audit look the
same in that file. docs/design-notes/synthetic-test-set.md specifies a doctored
VCF for this cell and the handout does not build one, so this does.

The records go on the end rather than in sorted position, which leaves the file
unsorted. That is why it goes only to --graph_method pangenie; the design note
says vg autoindex refuses an unsorted VCF, which I have not tested.

Three additions, one per category in bin/pangenie_graph_vcf.py:
  - an exact CHROM/POS/REF/ALT copy of an existing record, new ID  -> duplicate_of_record_N
  - a copy with N inside the ALT                                   -> never enters the graph
  - a copy with ID "."                                             -> id_replaced

  doctor_pangenome.py pangenome.vcf doctored.vcf
"""

import sys


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, dst = sys.argv[1], sys.argv[2]

    header, records = [], []
    with open(src) as fh:
        for line in fh:
            (header if line.startswith('#') else records).append(line.rstrip('\n'))
    if not records:
        sys.exit('doctor_pangenome.py: no records in %s' % src)

    def fields(i):
        return records[i].split('\t')

    # pick donors with a long enough ALT that an N is not the whole allele
    order = sorted(range(len(records)), key=lambda i: -len(fields(i)[4]))
    donors = order[:3] if len(order) >= 3 else order * 3

    added = []

    f = fields(donors[0])[:]          # exact duplicate, new ID
    f[2] = f[2] + '_DUPE'
    added.append(('duplicate', '\t'.join(f)))

    f = fields(donors[1])[:]          # N inside the ALT
    alt = f[4]
    mid = len(alt) // 2
    f[4] = alt[:mid] + 'N' * 10 + alt[mid + 10:]
    f[2] = f[2] + '_NALT'
    added.append(('N_in_ALT', '\t'.join(f)))

    f = fields(donors[2])[:]          # missing ID
    f[1] = str(int(f[1]) + 1)
    f[2] = '.'
    added.append(('no_ID', '\t'.join(f)))

    with open(dst, 'w') as fh:
        for h in header:
            fh.write(h + '\n')
        for r in records:
            fh.write(r + '\n')
        for _, r in added:
            fh.write(r + '\n')

    for what, r in added:
        g = r.split('\t')
        print('added %-10s %s:%s ID=%s len(ALT)=%d' % (what, g[0], g[1], g[2], len(g[4])))
    print('wrote %s: %d original + %d doctored' % (dst, len(records), len(added)))


if __name__ == '__main__':
    main()
