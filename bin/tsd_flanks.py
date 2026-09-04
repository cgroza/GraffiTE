#!/usr/bin/env python3
"""Extract the reference flanks the TSD search compares against.

For every insertion or deletion in the VCF, write two records to a FASTA:

  >ID__L   the WINDOW bases ending at POS, the anchor base included
  >ID__R   the WINDOW bases after the variant: from POS+1 for an insertion,
           from POS+len(REF)+1 for a deletion

These are the windows prepTSD.sh has always used. TSD_Match_v2.sh scores
offsets against them, so they are not to be moved. Sequences are written on one
line because the matcher reads each with `grep -A 1 | tail -n 1`.

Reads plain or BGZF FASTA through pysam. A window running off the start of a
contig is clamped; one running off the end is truncated by htslib.
"""

import argparse
import sys

import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--reference', required=True, help='plain or BGZF FASTA')
    ap.add_argument('--window', type=int, required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    fa = pysam.FastaFile(args.reference)
    contigs = set(fa.references)
    win = args.window

    n_records = 0
    n_written = 0
    with open(args.vcf) as vcf, open(args.out, 'w') as out:
        for line in vcf:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 5:
                continue
            n_records += 1
            chrom, pos, vid, ref, alt = f[0], int(f[1]), f[2], f[3], f[4]
            if chrom not in contigs:
                sys.exit(f'tsd_flanks.py: contig {chrom} of record {vid} '
                         f'is not in {args.reference}')
            if len(ref) < len(alt):
                r_start = pos
            elif len(ref) > len(alt):
                r_start = pos + len(ref)
            else:
                continue
            # pysam.fetch takes 0-based half-open coordinates.
            l_seq = fa.fetch(chrom, max(0, pos - win), pos)
            r_seq = fa.fetch(chrom, r_start, r_start + win)
            out.write(f'>{vid}__L\n{l_seq}\n>{vid}__R\n{r_seq}\n')
            n_written += 1

    print(f'tsd_flanks.py: {n_records} records, flanks written for {n_written}',
          file=sys.stderr)
    if n_records and not n_written:
        sys.exit('tsd_flanks.py: no flank could be extracted; the reference '
                 'and the VCF do not describe the same sequence')


if __name__ == '__main__':
    main()
