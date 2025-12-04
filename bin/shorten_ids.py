#!/usr/bin/python3

import vcfpy
import argparse
import sys

parser = argparse.ArgumentParser("Replace long truvari variant IDs with shorter ones.")
parser.add_argument("--vcf_in", metavar="vcf_in", type = str, nargs = 1,
                    help = "Input VCF to shorten ID field.")

parser.add_argument("--vcf_out", metavar="vcf_out", type = str, nargs = 1,
                    help = "Output VCF file.")

args = parser.parse_args()

reader = vcfpy.Reader.from_path(args.vcf_in[0])
writer = vcfpy.Writer.from_path(args.vcf_out[0], reader.header)

i = 0
for record in reader:
    record.ID = [record.ID[0] + '_' + str(i)]
    writer.write_record(record)
    i += 1
