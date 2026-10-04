#!/usr/bin/python3

# USAGE: breakgaps.py <assembly.fa> [min_gap]
#
# Splits each scaffold into contigs at runs of at least min_gap N. A shorter
# run is an unknown base inside a contig, not a gap between two, and splitting
# there cuts any insertion that carries it. Without min_gap every run splits.

import pysam
import sys
import re

fasta = pysam.FastaFile(sys.argv[1])
min_gap = int(sys.argv[2]) if len(sys.argv) > 2 else 1
gap = re.compile(r'N{%d,}' % min_gap)

for scaffold in fasta.references:
    scaffold_seq = fasta[scaffold]
    contig_count = 0
    contigs = gap.split(scaffold_seq)

    for contig in contigs:
        if len(contig) == 0:
            continue
        contig_name = ">" + scaffold + "_" + str(contig_count)
        print(contig_name)
        print(contig)
        contig_count = contig_count + 1
