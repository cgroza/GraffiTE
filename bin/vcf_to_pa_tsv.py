#!/usr/bin/env python3
"""Convert a VCF to a flat presence/absence TSV.

- INFO fields are flattened to columns (ordered by ##INFO header lines).
  Missing values are written as NA. Flag-type (no '=') INFO entries become '1'.
- FORMAT column is dropped.
- Per-sample columns report TE presence based on SVTYPE:
    INS: any ALT allele in GT -> 1 ; all ref -> 0
    DEL: any ALT allele in GT -> 0 ; all ref -> 1
  Missing genotypes (./., .) -> NA.
"""

import argparse
import re
import sys

INFO_HEADER_RE = re.compile(r'^##INFO=<ID=([^,]+),')


def parse_info(info):
    d = {}
    if info == '.' or not info:
        return d
    for kv in info.split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k] = v
        else:
            d[kv] = ''  # flag
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

    info_ids = []
    seen_info = set()
    samples = []
    header_written = False

    for line in fin:
        line = line.rstrip('\n')
        if line.startswith('##'):
            m = INFO_HEADER_RE.match(line)
            if m and m.group(1) not in seen_info:
                info_ids.append(m.group(1))
                seen_info.add(m.group(1))
            continue
        if line.startswith('#CHROM'):
            cols = line.split('\t')
            samples = cols[9:] if len(cols) > 9 else []
            header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER'] + info_ids + samples
            fout.write('\t'.join(header) + '\n')
            header_written = True
            continue
        if not line or not header_written:
            continue
        fields = line.split('\t')
        if len(fields) < 8:
            continue
        chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]
        info_d = parse_info(info)
        info_vals = []
        for k in info_ids:
            if k in info_d:
                v = info_d[k]
                info_vals.append('1' if v == '' else v)
            else:
                info_vals.append('NA')
        svtype = info_d.get('SVTYPE', '')
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
        fout.write('\t'.join([chrom, pos, vid, ref, alt, qual, filt] + info_vals + sample_vals) + '\n')

    if fin is not sys.stdin:
        fin.close()
    if fout is not sys.stdout:
        fout.close()


if __name__ == '__main__':
    main()
