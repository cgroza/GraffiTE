#!/usr/bin/env python3
"""What the breakscaf cell has to show.

The haplotype assemblies carry the reference's 120 bp run of N on chr1, and
single N bases inside the planted insertions that copy them from their Dfam
consensus. --break_scaffolds has to cut at the run and leave the insertions
whole, so discovery with it must find the same records as discovery without.

usage: check_breakscaf.py <runs dir> <assemblies.csv> <work dir> <min gap>
"""
import csv
import glob
import gzip
import os
import re
import sys


def fasta(path):
    seqs, name = {}, None
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: ''.join(v) for k, v in seqs.items()}


def records(vcf):
    """(CHROM, POS, length change) -> True when the longer allele carries an N."""
    out = {}
    with open(vcf) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t')
            longer = f[4] if len(f[4]) > len(f[3]) else f[3]
            out[(f[0], int(f[1]), len(f[4]) - len(f[3]))] = 'N' in longer[1:].upper()
    return out


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    runs, asm_csv, work, min_gap = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    gap = re.compile(r'N{%d,}' % min_gap)
    fail = 0

    def check(label, good, detail=''):
        nonlocal fail
        print(f"  [{' ok ' if good else 'FAIL'}] {label}" + (f": {detail}" if detail and not good else ''))
        fail += not good

    # what each assembly should become, from the input itself
    expected = {}
    for row in csv.DictReader(open(asm_csv)):
        seqs = fasta(row['path'])
        pieces = sum(1 for s in seqs.values() for p in gap.split(s) if p)
        short_n = sum(len(m.group()) for s in seqs.values() for m in re.finditer(r'N+', s)
                      if len(m.group()) < min_gap)
        expected[os.path.basename(row['path'])] = (len(seqs), pieces, short_n)

    trace = os.path.join(runs, 'breakscaf_break', 'nextflow_trace.txt')
    hashes = [r['hash'] for r in csv.DictReader(open(trace), delimiter='\t')
              if r['name'].startswith('break_scaffold')]
    check(f"break_scaffold ran once per assembly ({len(expected)})", len(hashes) == len(expected),
          f"{len(hashes)} tasks")
    for h in hashes:
        for out in glob.glob(os.path.join(work, h + '*', 'broken', '*.fa.gz')):
            asm = os.path.basename(out)[:-len('.fa.gz')]
            n_in, want, short_n = expected.get(asm, (0, -1, -1))
            seqs = fasta(out)
            got_short = sum(s.count('N') for s in seqs.values())
            left = sum(1 for s in seqs.values() if gap.search(s))
            check(f"{asm}: {n_in} sequences cut into {want}, {short_n} N in shorter runs kept",
                  len(seqs) == want and got_short == short_n and left == 0,
                  f"{len(seqs)} pieces, {got_short} N kept, {left} pieces still holding a gap")

    a = records(os.path.join(runs, 'breakscaf_nobreak', '3_TSD_search', 'pangenome.vcf'))
    b = records(os.path.join(runs, 'breakscaf_break', '3_TSD_search', 'pangenome.vcf'))
    with_n = sum(a.values())
    check("discovery without the break found records, some carrying an N", len(a) > 0 and with_n > 0,
          f"{len(a)} records, {with_n} with an N")
    check(f"the same {len(a)} records with --break_scaffolds ({with_n} with an N)", a.keys() == b.keys(),
          f"{len(a.keys() - b.keys())} lost, {len(b.keys() - a.keys())} new")

    print(f"breakscaf fail={fail}")
    sys.exit(1 if fail else 0)


if __name__ == '__main__':
    main()
