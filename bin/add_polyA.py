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
MAX_SLACK = 5        # max bp between tail end and the (TSD-trimmed) terminus

INFO_HEADER = ('##INFO=<ID=polyA,Number=1,Type=String,'
               'Description="TRUE if an imperfect polyA (or polyT on minus '
               'strand) tail is detected anchored to the 3\' (or 5\') end of a '
               'single-hit TE insertion/deletion after trimming any exact TSD '
               'suffix/prefix; FALSE if no tail found; NA if n_hits>1. '
               f'Params: min_len={MIN_LEN}, '
               f'min_purity={MIN_PURITY}, max_slack={MAX_SLACK}.">')


def parse_info(info):
    d = {}
    for kv in info.split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k] = v
        else:
            d[kv] = True
    return d


def has_anchored_tail(seq, base):
    """Return True iff seq contains a window of length >= MIN_LEN with
    >= MIN_PURITY fraction of `base` whose END lies within MAX_SLACK bp
    of the 3' terminus of seq."""
    n = len(seq)
    if n < MIN_LEN:
        return False
    seq = seq.upper()
    base = base.upper()
    ps = [0] * (n + 1)
    for i, c in enumerate(seq):
        ps[i + 1] = ps[i] + (1 if c == base else 0)
    # window [i, i+L) must have j = i+L in [n-MAX_SLACK, n]
    for L in range(MIN_LEN, n + 1):
        need = L * MIN_PURITY
        i_lo = max(0, n - MAX_SLACK - L)
        i_hi = n - L  # inclusive
        for i in range(i_lo, i_hi + 1):
            if ps[i + L] - ps[i] >= need:
                return True
    return False


def detect_polyA(variant_seq, strand, tsd):
    """variant_seq: inserted/deleted sequence (without the anchor base).
       strand: '+' or 'C'.
       tsd: reported TSD sequence (may be empty/'.')."""
    seq = variant_seq.upper()
    tsd_up = (tsd or '').upper()
    if strand == '+':
        # polyA at 3' end: trim an exact TSD suffix if present, then scan
        # anchored to the new 3' terminus.
        if tsd_up and seq.endswith(tsd_up):
            seq = seq[:-len(tsd_up)]
        return has_anchored_tail(seq, 'A')
    elif strand == 'C':
        # polyT at 5' end: trim an exact TSD prefix, reverse, then scan
        # anchored to the (reversed) 3' terminus for T.
        if tsd_up and seq.startswith(tsd_up):
            seq = seq[len(tsd_up):]
        return has_anchored_tail(seq[::-1], 'T')
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
        if not isinstance(tsd, str) or tsd in ('.', ''):
            tsd = ''

        if variant_seq and detect_polyA(variant_seq, strands, tsd):
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
