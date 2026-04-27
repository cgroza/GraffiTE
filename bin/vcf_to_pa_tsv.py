#!/usr/bin/env python3
"""Convert a GraffiTE VCF to a flat presence/absence TSV.

Fixed output schema (21 columns + one column per sample):
  CHROM POS END ID SVTYPE SVLEN n_hits match_lengths repeat_ids matching_classes
  fragmts RM_hit_strands RM_hit_IDs total_match_length total_match_span L1_5PINV
  ULTRA_TR ULTRA_TR_span total_repeat_span TSD polyA <sample1> <sample2> ...

- Missing INFO fields are written as NA.
- FORMAT column is dropped.
- Per-sample values: TE presence based on SVTYPE:
    INS: any ALT allele in GT -> 1 ; all ref -> 0
    DEL: any ALT allele in GT -> 0 ; all ref -> 1
  Missing genotypes (./., .) -> NA.
"""

import argparse
import re
import sys

INFO_COLS = [
    'END', 'SVTYPE', 'SVLEN', 'n_hits', 'match_lengths', 'repeat_ids',
    'matching_classes', 'fragmts', 'RM_hit_strands', 'RM_hit_IDs',
    'total_match_length', 'total_match_span', 'L1_5PINV', 'ULTRA_TR',
    'ULTRA_TR_span', 'total_repeat_span', 'TSD', 'polyA',
]


def parse_info(info):
    d = {}
    if info == '.' or not info:
        return d
    for kv in info.split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k] = v
        else:
            d[kv] = ''
    return d


def gt_presence(gt, svtype):
    if not gt or gt in ('.', './.', '.|.'):
        return 'NA'
    alleles = re.split(r'[/|]', gt)
    if all(a == '.' for a in alleles):
        return 'NA'
    has_alt = any(a not in ('0', '.', '') for a in alleles)
    if svtype == 'INS':
        return '1' if has_alt else '0'
    if svtype == 'DEL':
        return '0' if has_alt else '1'
    return 'NA'


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('vcf', nargs='?', default='-')
    ap.add_argument('-o', '--output', default='-')
    args = ap.parse_args()

    fin = sys.stdin if args.vcf == '-' else open(args.vcf)
    fout = sys.stdout if args.output == '-' else open(args.output, 'w')

    samples = []
    header_written = False

    for line in fin:
        line = line.rstrip('\n')
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            cols = line.split('\t')
            samples = cols[9:] if len(cols) > 9 else []
            header = ['CHROM', 'POS'] + INFO_COLS[:1] + ['ID'] + INFO_COLS[1:] + samples
            # Re-order: CHROM POS END ID SVTYPE ... (END before ID, per spec)
            header = ['CHROM', 'POS', 'END', 'ID'] + INFO_COLS[1:] + samples
            fout.write('\t'.join(header) + '\n')
            header_written = True
            continue
        if not line or not header_written:
            continue
        fields = line.split('\t')
        if len(fields) < 8:
            continue
        chrom, pos, vid, _ref, _alt, _qual, _filt, info = fields[:8]
        info_d = parse_info(info)
        svtype = info_d.get('SVTYPE', '')
        end_val = info_d.get('END', 'NA')
        other_vals = []
        for k in INFO_COLS[1:]:  # skip END (already placed)
            if k in info_d:
                v = info_d[k]
                other_vals.append('1' if v == '' else v)
            else:
                other_vals.append('NA')
        sample_vals = []
        if samples and len(fields) > 9:
            fmt = fields[8].split(':')
            try:
                gt_idx = fmt.index('GT')
            except ValueError:
                gt_idx = 0
            for s_data in fields[9:]:
                gt = s_data.split(':')[gt_idx] if s_data else '.'
                sample_vals.append(gt_presence(gt, svtype))
        row = [chrom, pos, end_val, vid] + other_vals + sample_vals
        fout.write('\t'.join(row) + '\n')

    if fin is not sys.stdin:
        fin.close()
    if fout is not sys.stdout:
        fout.close()


if __name__ == '__main__':
    main()
