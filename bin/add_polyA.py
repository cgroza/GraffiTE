#!/usr/bin/env python3
"""
Annotate a GraffiTE VCF with polyA=TRUE/FALSE in the INFO column.

Rules:
  - n_hits > 1           -> polyA=NA (skip search)
  - RM_hit_strands='+'   -> look for polyA run near the 3' end of the variant sequence
  - RM_hit_strands='C'   -> look for polyT run near the 5' end of the variant sequence
  - INS: variant sequence = ALT[1:]; DEL: variant sequence = REF[1:]
  - TSD may be absent or truncated at the breakpoint, so we scan a window
    of (len(TSD) + FLANK) bp at the relevant end rather than hard-trimming.

A tail is called if a sliding window of length >= MIN_LEN within the scan
region has A (or T) fraction >= MIN_PURITY.
"""

import argparse
import re
import sys

MIN_LEN = 8          # minimum tail length
MIN_PURITY = 0.8     # minimum A/T fraction within the tail window
FLANK = 50           # bp to scan past the reported TSD length

INFO_HEADER = ('##INFO=<ID=polyA,Number=1,Type=String,'
               'Description="TRUE if an imperfect polyA (or polyT on minus '
               'strand) tail is detected near the breakpoint of a single-hit '
               'TE insertion/deletion; FALSE if no tail found; NA if n_hits>1. '
               f'Params: min_len={MIN_LEN}, '
               f'min_purity={MIN_PURITY}, flank={FLANK}.">')


def parse_info(info):
    d = {}
    for kv in info.split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k] = v
        else:
            d[kv] = True
    return d


def has_tail(seq, base):
    """Return True if seq contains a window of length >= MIN_LEN with
    >= MIN_PURITY fraction of `base`."""
    if len(seq) < MIN_LEN:
        return False
    seq = seq.upper()
    base = base.upper()
    # Try every window size from MIN_LEN up to len(seq); accept if any window matches.
    # Efficient check: for each window of length MIN_LEN, compute purity.
    # Additionally try to extend: scan for maximal runs of base allowing mismatches.
    n = len(seq)
    # Prefix sum of base occurrences
    ps = [0] * (n + 1)
    for i, c in enumerate(seq):
        ps[i + 1] = ps[i] + (1 if c == base else 0)
    # Check any window of length >= MIN_LEN with purity >= MIN_PURITY.
    for L in range(MIN_LEN, n + 1):
        need = L * MIN_PURITY
        for i in range(0, n - L + 1):
            if ps[i + L] - ps[i] >= need:
                return True
        # early exit: if max possible count in any window of length L is below need,
        # longer windows will also fail at that start, but could still succeed
        # elsewhere. Simpler to just keep looping; sequences here are short.
    return False


def detect_polyA(variant_seq, strand, tsd_len):
    """variant_seq: inserted/deleted sequence (without the anchor base).
       strand: '+' or 'C'.
       tsd_len: reported TSD length.
       Returns True if a polyA (or polyT) tail is called."""
    window = tsd_len + FLANK
    if strand == '+':
        region = variant_seq[-window:] if window < len(variant_seq) else variant_seq
        return has_tail(region, 'A')
    elif strand == 'C':
        region = variant_seq[:window] if window < len(variant_seq) else variant_seq
        return has_tail(region, 'T')
    return False


def annotate_record(fields):
    info_str = fields[7]
    info = parse_info(info_str)

    svtype = info.get('SVTYPE', '')
    n_hits = info.get('n_hits')
    strands = info.get('RM_hit_strands', '')

    try:
        n_hits_i = int(n_hits) if n_hits is not None else 0
    except ValueError:
        n_hits_i = 0

    if n_hits_i > 1:
        call = 'NA'
    else:
        call = 'FALSE'

    if n_hits_i == 1 and strands in ('+', 'C'):
        ref = fields[3]
        alt = fields[4]
        if svtype == 'INS':
            variant_seq = alt[1:]
        elif svtype == 'DEL':
            variant_seq = ref[1:]
        else:
            variant_seq = ''

        tsd = info.get('TSD', '')
        tsd_len = len(tsd) if isinstance(tsd, str) and tsd not in ('.', '') else 0

        if variant_seq and detect_polyA(variant_seq, strands, tsd_len):
            call = 'TRUE'

    # Append polyA to INFO (replace if somehow already present)
    new_info = re.sub(r'(^|;)polyA=[^;]*', '', info_str).strip(';')
    new_info = f'{new_info};polyA={call}' if new_info else f'polyA={call}'
    fields[7] = new_info
    return fields


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('vcf', nargs='?', default='-',
                    help='Input VCF (default: stdin)')
    ap.add_argument('-o', '--output', default='-',
                    help='Output VCF (default: stdout)')
    args = ap.parse_args()

    fin = sys.stdin if args.vcf == '-' else open(args.vcf)
    fout = sys.stdout if args.output == '-' else open(args.output, 'w')

    header_info_injected = False
    saw_header = False

    for line in fin:
        if line.startswith('##'):
            saw_header = True
            # Inject our INFO line before the #CHROM line; track last INFO seen.
            fout.write(line)
            continue
        if line.startswith('#CHROM'):
            saw_header = True
            if not header_info_injected:
                fout.write(INFO_HEADER + '\n')
                header_info_injected = True
            fout.write(line)
            continue

        if not line.strip():
            fout.write(line)
            continue

        # Data line. If no header was present (headerless test VCF), emit nothing extra.
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 8:
            fout.write(line)
            continue
        fields = annotate_record(fields)
        fout.write('\t'.join(fields) + '\n')

    if fin is not sys.stdin:
        fin.close()
    if fout is not sys.stdout:
        fout.close()


if __name__ == '__main__':
    main()
