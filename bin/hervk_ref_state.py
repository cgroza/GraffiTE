#!/usr/bin/env python3
"""
Determine the HERV-K (HML-2) state of the *reference* allele at each candidate
locus.

The architecture of the SV sequence alone cannot always say what the reference
carries. Two full terminal LTRs prove the reference was empty; one LTR split
across the termini proves it held a solo LTR. But at the degenerate breakpoint
k = 0 the SV shows one whole LTR and nothing at the other end, which looks
identical whether the reference held a solo LTR (insertion right at the LTR
boundary) or nothing at all (a one-LTR-truncated provirus into an empty site).
About half the CaG candidates land in that degenerate form, so this is the
workhorse of the classifier and not a rare fallback.

The check is direct: cut a window out of the reference around the SV footprint,
mask it, and read off what HML-2 sequence is actually there.

    null      no HML-2 sequence at the locus
    solo      a single LTR, no internal region
    provirus  internal region flanked by ~2 LTRs worth of LTR sequence
    partial   internal region with only ~1 LTR, or a fragmentary LTR
    unknown   window could not be evaluated

Usage:
    hervk_ref_state.py --vcf in.vcf --reference ref.fa --te-library lib.fa \
        --out refstate.tsv [--ids ids.txt] [--flank 1500] [--threads 4]
    hervk_ref_state.py --vcf in.vcf --reference ref.fa \
        --rm-annotation rmsk.bed --out refstate.tsv
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hervk_arch import (DEFAULTS as ARCH_DEFAULTS, INT_FAMILIES, arch_string,
                        is_ltr_family, ltr_len, parse_rm_out, reassign_sine_r,
                        tile_hits)

DEFAULTS = {
    "flank": 1500,
    # Max gap (bp) between HML-2 fragments still counted as one element.
    "element_gap": 1000,
    # Min HML-2 bp in the window before the locus is anything but `null`.
    "min_hml2_bp": 100,
    # Min internal-region bp before the element counts as carrying an INT.
    "min_int_bp": 200,
    # Tolerance (bp) when matching total LTR bp against 1x or 2x consensus.
    "ltr_tol": 250,
}


# -------- Candidate footprints from the VCF --------
def read_footprints(vcf_path, wanted=None):
    """Return [(id, chrom, fp_start, fp_end, svlen)] with 1-based inclusive
    footprints: a point at POS for insertions, POS..POS+|SVLEN| for deletions.

    Polarity comes from len(ALT) - len(REF), never from INFO/SVTYPE: the
    multi-sample merge strips INFO (module/main.nf `bcftools annotate -x INFO`)
    and only SVLEN is added back, so SVTYPE is not reliably present.
    """
    out = []
    opener = open
    if vcf_path.endswith('.gz'):
        import gzip
        opener = gzip.open
    with opener(vcf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 8:
                continue
            sv_id = f[2]
            if wanted is not None and sv_id not in wanted:
                continue
            pos = int(f[1])
            svlen = len(f[4].split(',')[0]) - len(f[3])
            if svlen >= 0:
                out.append((sv_id, f[0], pos, pos, svlen))
            else:
                out.append((sv_id, f[0], pos, pos + abs(svlen), svlen))
    return out


# -------- Reference window extraction + masking --------
def extract_windows(footprints, reference, flank, workdir):
    """Cut one FASTA record per candidate; return {id: (chrom, win_start)}."""
    regions, offsets = [], {}
    for sv_id, chrom, fp_start, fp_end, _ in footprints:
        start = max(1, fp_start - flank)
        end = fp_end + flank
        regions.append(f'{chrom}:{start}-{end}')
        offsets[sv_id] = (chrom, start)

    region_file = os.path.join(workdir, 'regions.txt')
    with open(region_file, 'w') as fh:
        fh.write('\n'.join(regions) + '\n')

    raw = subprocess.run(['samtools', 'faidx', reference, '-r', region_file],
                         capture_output=True, text=True, check=True).stdout

    # Rename each record from chrom:start-end to the SV id so the RepeatMasker
    # query ids come back keyed the way everything else in GraffiTE is keyed.
    fasta = os.path.join(workdir, 'ref_windows.fa')
    ids = [fp[0] for fp in footprints]
    with open(fasta, 'w') as fh:
        i = -1
        for line in raw.splitlines():
            if line.startswith('>'):
                i += 1
                fh.write(f'>{ids[i]}\n')
            else:
                fh.write(line + '\n')
    return fasta, offsets


def run_repeatmasker(fasta, te_library, threads, workdir):
    rm_dir = os.path.join(workdir, 'rm')
    os.makedirs(rm_dir, exist_ok=True)
    subprocess.run(['RepeatMasker', '-lib', te_library, '-s',
                    '-dir', rm_dir, '-pa', str(max(1, threads)), fasta],
                   check=True, capture_output=True, text=True)
    out = os.path.join(rm_dir, os.path.basename(fasta) + '.out')
    if not os.path.exists(out):
        # RepeatMasker omits the .out when nothing was masked anywhere.
        return None
    return out


def hits_from_annotation(path, footprints, flank):
    """Load reference HML-2 hits from a precomputed track instead of masking.

    Accepts a RepeatMasker .out over the reference, or a BED4
    (chrom, start, end, repeat_name). Coordinates are converted to
    window-relative so the tiling code sees the same shape either way.
    """
    windows = {}
    for sv_id, chrom, fp_start, fp_end, _ in footprints:
        windows.setdefault(chrom, []).append(
            (sv_id, max(1, fp_start - flank), fp_end + flank))

    hits = {}

    def add(chrom, start, end, name, klass, sw, cons_start, cons_end, cons_left, strand):
        for sv_id, wstart, wend in windows.get(chrom, []):
            if end < wstart or start > wend:
                continue
            hits.setdefault(sv_id, []).append({
                'sw': sw, 'qstart': max(start, wstart) - wstart + 1,
                'qend': min(end, wend) - wstart + 1, 'strand': strand,
                'name': name, 'klass': klass,
                'cons_start': cons_start, 'cons_end': cons_end,
                'cons_len': (cons_end or 0) + (cons_left or 0), 'link': name,
            })

    with open(path) as fh:
        if path.endswith('.out'):
            for line in fh:
                f = line.split()
                if len(f) < 15 or not f[0].isdigit():
                    continue
                nums = [int(x.strip('()')) for x in (f[11], f[12], f[13])]
                if f[8] == 'C':
                    cons_left, cons_end, cons_start = nums
                else:
                    cons_start, cons_end, cons_left = nums
                add(f[4], int(f[5]), int(f[6]), f[9], f[10], int(f[0]),
                    cons_start, cons_end, cons_left, f[8])
        else:
            for line in fh:
                if line.startswith(('#', 'track', 'browser')):
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) < 4:
                    continue
                name = f[3]
                klass = f[4] if len(f) > 4 and '/' in f[4] else 'LTR/ERVK'
                # BED is half-open 0-based; RM-style coords are 1-based inclusive.
                add(f[0], int(f[1]) + 1, int(f[2]), name, klass, 1000,
                    None, None, None, f[5] if len(f) > 5 else '+')
    return hits


# -------- State calling --------
def cluster_elements(frags, cfg):
    """Merge HML-2 fragments separated by <= element_gap into single elements,
    so a neighbouring unrelated HML-2 copy in the window is not pooled in."""
    hml2 = sorted((f for f in frags
                   if is_ltr_family(f['name']) or f['name'] in INT_FAMILIES),
                  key=lambda f: f['qstart'])
    elements, current = [], []
    for f in hml2:
        if current and f['qstart'] - current[-1]['qend'] > cfg['element_gap']:
            elements.append(current)
            current = []
        current.append(f)
    if current:
        elements.append(current)
    return elements


def call_state(element, cfg):
    """Reference state for one clustered HML-2 element."""
    if not element:
        return 'null', 0, 0, ''

    ltr_bp = sum(f['bp'] for f in element if is_ltr_family(f['name']))
    int_bp = sum(f['bp'] for f in element if f['name'] in INT_FAMILIES)
    arch = arch_string(element)
    if ltr_bp + int_bp < cfg['min_hml2_bp']:
        return 'null', ltr_bp, int_bp, arch

    ltrs = [f for f in element if is_ltr_family(f['name'])]
    fam = max(ltrs, key=lambda f: f['bp'])['name'] if ltrs else 'LTR5_Hs'
    L = ltr_len(fam)

    if int_bp < cfg['min_int_bp']:
        state = 'solo' if abs(ltr_bp - L) <= cfg['ltr_tol'] else 'partial'
    else:
        state = 'provirus' if ltr_bp >= 1.5 * L else 'partial'
    return state, ltr_bp, int_bp, arch


def evaluate(footprints, rm_hits, offsets, cfg):
    """Pick the HML-2 element at (or nearest) each footprint and call its state."""
    arch_cfg = dict(ARCH_DEFAULTS)
    results = {}
    for sv_id, chrom, fp_start, fp_end, svlen in footprints:
        hits = rm_hits.get(sv_id, [])
        if not hits:
            results[sv_id] = {'state': 'null', 'ltr_bp': 0, 'int_bp': 0,
                              'dist': '', 'arch': 'NONE', 'chrom': chrom,
                              'elem_start': '', 'elem_end': ''}
            continue
        frags = reassign_sine_r(tile_hits(hits, arch_cfg), arch_cfg)
        elements = cluster_elements(frags, cfg)
        if not elements:
            results[sv_id] = {'state': 'null', 'ltr_bp': 0, 'int_bp': 0,
                              'dist': '', 'arch': 'NONE', 'chrom': chrom,
                              'elem_start': '', 'elem_end': ''}
            continue

        win_start = offsets.get(sv_id, (chrom, max(1, fp_start - cfg['flank'])))[1]
        fp_lo = fp_start - win_start + 1
        fp_hi = fp_end - win_start + 1

        def distance(el):
            lo = min(f['qstart'] for f in el)
            hi = max(f['qend'] for f in el)
            if hi < fp_lo:
                return fp_lo - hi
            if lo > fp_hi:
                return lo - fp_hi
            return 0

        best = min(elements, key=distance)
        state, ltr_bp, int_bp, arch = call_state(best, cfg)
        results[sv_id] = {'state': state, 'ltr_bp': ltr_bp, 'int_bp': int_bp,
                          'dist': distance(best), 'arch': arch,
                          'chrom': chrom,
                          'elem_start': min(f['qstart'] for f in best) + win_start - 1,
                          'elem_end': max(f['qend'] for f in best) + win_start - 1}
    return results


COLUMNS = ['id', 'ref_state', 'ref_ltr_bp', 'ref_int_bp', 'ref_dist',
           'ref_elem_chrom', 'ref_elem_start', 'ref_elem_end', 'ref_arch']


def write_tsv(results, path):
    with open(path, 'w') as fh:
        fh.write('\t'.join(COLUMNS) + '\n')
        for sv_id in sorted(results):
            r = results[sv_id]
            fh.write('\t'.join([sv_id, r['state'], str(int(r['ltr_bp'])),
                                str(int(r['int_bp'])), str(r['dist']),
                                r.get('chrom', ''), str(r.get('elem_start', '')),
                                str(r.get('elem_end', '')), r['arch']]) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--reference', required=True)
    ap.add_argument('--te-library')
    ap.add_argument('--rm-annotation',
                    help='precomputed reference annotation (.out or BED) instead of masking')
    ap.add_argument('--out', required=True)
    ap.add_argument('--ids', help='optional file of SV IDs to restrict to')
    ap.add_argument('--flank', type=int, default=DEFAULTS['flank'])
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--keep-temp', action='store_true')
    args = ap.parse_args()

    if not args.rm_annotation and not args.te_library:
        sys.exit('ERROR: one of --te-library (to mask) or --rm-annotation is required')

    cfg = dict(DEFAULTS)
    cfg['flank'] = args.flank

    wanted = None
    if args.ids:
        with open(args.ids) as fh:
            wanted = {ln.strip() for ln in fh if ln.strip()}

    footprints = read_footprints(args.vcf, wanted)
    if not footprints:
        write_tsv({}, args.out)
        sys.stderr.write('hervk_ref_state: no candidates; empty table written\n')
        return

    if args.rm_annotation:
        rm_hits = hits_from_annotation(args.rm_annotation, footprints, args.flank)
        offsets = {fp[0]: (fp[1], max(1, fp[2] - args.flank)) for fp in footprints}
    else:
        workdir = tempfile.mkdtemp(prefix='hervk_ref_')
        try:
            fasta, offsets = extract_windows(footprints, args.reference,
                                             args.flank, workdir)
            rm_out = run_repeatmasker(fasta, args.te_library, args.threads, workdir)
            rm_hits = parse_rm_out([rm_out]) if rm_out else {}
        finally:
            if not args.keep_temp:
                shutil.rmtree(workdir, ignore_errors=True)

    results = evaluate(footprints, rm_hits, offsets, cfg)
    write_tsv(results, args.out)
    counts = {}
    for r in results.values():
        counts[r['state']] = counts.get(r['state'], 0) + 1
    summary = ', '.join(f'{k}={v}' for k, v in sorted(counts.items()))
    sys.stderr.write(f'hervk_ref_state: {len(results)} loci ({summary})\n')


if __name__ == '__main__':
    main()
