import pysam
import sys
import re



fasta = pysam.FastaFile(sys.argv[1])

for scaffold in fasta.references:
    scaffold_seq = fasta[scaffold]
    contig_count = 0
    contigs = re.split(r'N+', scaffold_seq)

    for contig in contigs:
        contig_name = ">" + scaffold + "_" + str(contig_count)
        print(contig_name)
        print(contig)
        contig_count = contig_count + 1
