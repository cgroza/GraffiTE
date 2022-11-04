#!/usr/bin/python3

import vcfpy
import pysam
import argparse
import sys

parser = argparse.ArgumentParser("Replace deletion REF sequences in VCF with interval from FASTA.")
parser.add_argument("--ref", metavar="ref", type = str, nargs = 1,
                    help = "Reference genome FASTA.")

parser.add_argument("--vcf_in", metavar="vcf_in", type = str, nargs = 1,
                    help = "Input VCF to correct REF field.")

parser.add_argument("--vcf_out", metavar="vcf_out", type = str, nargs = 1,
                    help = "Output VCF file.")

args = parser.parse_args()

fasta = pysam.FastaFile(args.ref[0])
reader = vcfpy.Reader.from_path(args.vcf_in[0])
writer = vcfpy.Writer.from_path(args.vcf_out[0], reader.header)

i = 1
for record in reader:
    print("Record", i, file = sys.stderr)
    ref = record.REF
    alt = record.ALT[0].value
    chrom = record.CHROM
    pos = record.POS
    svlen = len(ref) - len(alt)

    # insertions
    if(svlen < 0):
        new_ref = fasta.fetch(chrom, pos - 1, pos)
        record.REF = new_ref
        new_alt = new_ref + alt[1:]
        record.ALT[0].value = new_alt
    # deletions, SNVs
    else:
        new_ref = fasta.fetch(chrom, pos - 1, pos + svlen)
        new_alt = new_ref[0]
        record.REF = new_ref
        record.ALT[0].value = new_alt

    writer.write_record(record)
    i = i + 1
