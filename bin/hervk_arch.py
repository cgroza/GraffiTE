#!/usr/bin/env python3
"""
HERV-K (HML-2) architecture extraction from raw RepeatMasker output.

GraffiTE's `annotate_vcf.R` groups RepeatMasker fragments by link ID and
collapses each group to a single name (top SW score) plus an "(x)" marker when
the fragments disagree. That is the right summary for most TE classes, but for
HML-2 it destroys exactly the information the classifier needs: a full provirus
collapses to `HERVK-int(x)` spanning the whole SV, with zero LTR bp.

This module reads the raw `indels.fa.out` instead and recovers, per SV:

  * a winner-take-all tiling of hits on the query axis (so overlapping
    LTR5_Hs / SVA_A calls are counted once rather than summed twice),
  * SVA SINE-R reassignment -- an SVA hit at SVA-consensus >= ~900 abutting an
    HML-2 hit is LTR5-derived sequence, not a separate element,
  * fragments ordered 5'->3' along the *element* (not the query), and
  * the two diagnostic signatures:

      ARCH_2LTR   two full-length terminal LTRs (both ~1..L)
                  -> the SV carries a complete provirus
      ARCH_PERM   one LTR split across the termini, consensus intervals
                  complementary: 5' fragment [k+1..L], 3' fragment [1..k]
                  -> the SV was placed inside a pre-existing solo LTR

`k` is the alignment breakpoint inside the reference LTR. It is a property of
the alignment, not of the biology: the inserted length is the same for every k,
so the gap penalty cancels and the placement is decided by the handful of
substitutions separating the two LTR copies. k therefore varies between
haplotypes and between callers, and nothing downstream may key on its value.
At k = 0 the permutation is degenerate -- one full LTR and nothing at the other
terminus -- which is indistinguishable from a one-LTR-truncated provirus, and
must fall through to the reference check.

Usage:
    hervk_arch.py --rm-out DIR_OR_FILE [...] --out arch.tsv [--ids ids.txt]
"""

import argparse
import glob
import os
import re
import sys
from collections import Counter

# -------- HML-2 reference architecture --------
LTR_CONSENSUS_LEN = {'LTR5_Hs': 968, 'LTR5A': 1033, 'LTR5B': 968, 'LTR5': 968}
INT_CONSENSUS_LEN = 7536

# The HML-2 internal region is named differently by different libraries:
# Dfam calls it HERVK, RepBase-derived sets HERVK-int, and HERVK_int / HERVKint
# also occur. Match all of them, but only them -- HERVK9-int, HERVK11-int and
# HERVK14-int are separate ERV lineages and a loose prefix would sweep them in.
# Requiring the name to end right after the optional "int" is what keeps the
# digit-suffixed families out.
INT_RE = re.compile(r'^HERVK[-_]?(int(ernal)?)?$', re.IGNORECASE)
INT_FAMILIES = {'HERVK-int', 'HERVK'}   # kept for callers that want a literal set
SKIP_CLASSES = {'Simple_repeat', 'Low_complexity'}

DEFAULTS = {
    # Consensus start at/after which an SVA hit is SINE-R (HERV-K LTR derived).
    "sine_r_min": 900,
    # Max query gap (bp) for an SVA SINE-R hit to count as abutting HML-2.
    "sine_r_max_gap": 50,
    # Tolerance (bp) on consensus-boundary and complementarity comparisons.
    "perm_tol": 50,
    # An LTR terminus must start within this many bp of the SV end.
    "terminus_tol": 60,
    # Minimum INT bp for the SV to be considered to carry an internal region.
    "min_int_bp": 200,
    # Minimum bp for a tiled fragment to be reported at all.
    "min_frag_bp": 20,
}


def is_int_family(name):
    """HML-2 internal region, whatever the library calls it."""
    return bool(INT_RE.match(name))


def is_ltr_family(name):
    """HML-2 LTRs only.

    This was a startswith('LTR5') prefix test, which is wrong: a Dfam human
    library carries ~20 families on that prefix and only these four are HML-2.
    LTR57-int and LTR53-int are not even LTRs -- they are the internal regions
    of other ERV lineages -- and were being counted as LTR bp.
    """
    return name in LTR_CONSENSUS_LEN


def is_sva_family(name):
    return name.startswith('SVA_')


def ltr_len(name):
    return LTR_CONSENSUS_LEN.get(name, 968)


# -------- RepeatMasker .out parsing --------
def parse_rm_out(paths):
    """Return {qry_id: [hit, ...]} from one or more RepeatMasker .out files.

    Fields follow the cross_match layout that `annotate_vcf.R:read_rm_custom`
    also parses. Consensus coordinates are orientation-corrected here: for a
    'C' hit RepeatMasker prints "(left) end start", so the consensus interval
    is (field 14, field 13); for '+' it is (field 12, field 13).
    """
    hits = {}
    for path in paths:
        with open(path) as fh:
            for line in fh:
                f = line.split()
                if len(f) < 15 or not f[0].isdigit():
                    continue  # header, blank, or malformed
                try:
                    sw = int(f[0])
                    qstart, qend = int(f[5]), int(f[6])
                except ValueError:
                    continue
                strand = f[8]
                name, klass = f[9], f[10]
                if klass in SKIP_CLASSES:
                    continue
                nums = [int(x.strip('()')) for x in (f[11], f[12], f[13])]
                if strand == 'C':
                    cons_left, cons_end, cons_start = nums
                else:
                    cons_start, cons_end, cons_left = nums
                hits.setdefault(f[4], []).append({
                    'sw': sw, 'qstart': qstart, 'qend': qend, 'strand': strand,
                    'name': name, 'klass': klass,
                    'cons_start': cons_start, 'cons_end': cons_end,
                    'cons_len': cons_end + cons_left,
                    'link': f[14],
                })
    return hits


def expand_rm_paths(inputs):
    """Accept .out files, or directories/globs containing indels.fa.out."""
    paths = []
    for item in inputs:
        for match in sorted(glob.glob(item)) or [item]:
            if os.path.isdir(match):
                paths.extend(sorted(glob.glob(
                    os.path.join(match, '**', 'indels.fa.out'), recursive=True)))
            elif os.path.exists(match):
                paths.append(match)
    return paths


# -------- Query-axis tiling --------
def tile_hits(hits, cfg):
    """Winner-take-all on the query axis, highest SW score first.

    Summing `match_lengths` double-counts wherever RepeatMasker reports an
    SVA_A hit on top of an LTR5_Hs hit (chr6-78894317 does exactly that), so
    every bp is awarded to one hit only before any bp is totted up.
    """
    if not hits:
        return []
    lo = min(h['qstart'] for h in hits)
    hi = max(h['qend'] for h in hits)
    owner = [None] * (hi - lo + 1)
    for idx, h in enumerate(sorted(range(len(hits)),
                                   key=lambda i: -hits[i]['sw'])):
        for p in range(hits[h]['qstart'] - lo, hits[h]['qend'] - lo + 1):
            if owner[p] is None:
                owner[p] = h

    # One pass over `owner` instead of one per hit. The previous form,
    # `sum(1 for o in owner if o == i)` inside this loop, is O(n_hits * span):
    # chr1-120594342-DEL-25264467 is a 25.3 Mb DEL carrying 38430 hits, i.e.
    # ~9.7e11 comparisons, which ran 45 min at ~4% before being killed.
    # `owner[p]` already holds the winning index for every base, so counting it
    # once is the same arithmetic.
    claims = Counter(o for o in owner if o is not None)

    tiled = []
    for i, h in enumerate(hits):
        claimed = claims.get(i, 0)
        if claimed < cfg['min_frag_bp']:
            continue
        frag = dict(h)
        frag['bp'] = claimed
        tiled.append(frag)
    tiled.sort(key=lambda h: h['qstart'])
    return tiled


def reassign_sine_r(frags, cfg):
    """Relabel SVA SINE-R fragments abutting HML-2 sequence as LTR.

    SVA's SINE-R domain is HERV-K LTR derived, so RepeatMasker will sometimes
    win the terminal LTR of a provirus for SVA. Every SVA hit seen beside an
    HML-2 element in the CaG set sits at SVA consensus >= ~900, and its length
    completes the LTR exactly (chr12-58305931: 303 bp SVA + LTR5_Hs 304-968).
    """
    hml2 = [f for f in frags
            if is_ltr_family(f['name']) or is_int_family(f['name'])]
    if not hml2:
        return frags
    for f in frags:
        if not is_sva_family(f['name']):
            continue
        # A hit from a BED annotation has no consensus coordinates, so there
        # is no way to tell the SINE-R domain from the rest of SVA. Leave it
        # alone rather than guess -- and never compare None to an int.
        if f['cons_start'] is None or f['cons_start'] < cfg['sine_r_min']:
            continue
        gap = min(max(f['qstart'] - o['qend'], o['qstart'] - f['qend'], 0)
                  for o in hml2)
        if gap <= cfg['sine_r_max_gap']:
            f['reassigned_from'] = f['name']
            f['name'] = 'LTR5_Hs'
            f['klass'] = 'LTR/ERVK'
            # SVA consensus coordinates do not translate to LTR consensus
            # coordinates; mark them unknown and let the length carry it.
            f['cons_start'] = None
            f['cons_end'] = None
    return frags


def element_order(frags):
    """Order fragments 5'->3' along the element rather than along the query.

    A 'C'-strand element runs opposite to the query, so its terminal fragments
    are swapped relative to query order.
    """
    hml2 = [f for f in frags
            if is_ltr_family(f['name']) or is_int_family(f['name'])]
    pool = hml2 or frags
    rev = sum(f['sw'] for f in pool if f['strand'] == 'C') > \
          sum(f['sw'] for f in pool if f['strand'] != 'C')
    ordered = sorted(frags, key=lambda f: f['qstart'], reverse=rev)
    return ordered, ('C' if rev else '+')


# -------- Architecture --------
def _close(a, b, tol):
    return a is not None and b is not None and abs(a - b) <= tol


def architecture(frags, svlen, cfg):
    """Classify the internal architecture of one SV allele.

    Returns a dict with the tiled bp totals, the ordered architecture string,
    the number of full/partial terminal LTRs, the permutation point k when the
    ARCH_PERM signature holds, and a signature in
    {ARCH_2LTR, ARCH_PERM, ARCH_SOLO, ARCH_NONE}.
    """
    ordered, strand = element_order(frags)
    ltr_bp = sum(f['bp'] for f in ordered if is_ltr_family(f['name']))
    int_bp = sum(f['bp'] for f in ordered if is_int_family(f['name']))
    other_bp = sum(f['bp'] for f in ordered
                   if not is_ltr_family(f['name'])
                   and f['name'] not in INT_FAMILIES)

    out = {
        'strand': strand, 'n_frag': len(ordered),
        'ltr_bp': ltr_bp, 'int_bp': int_bp, 'other_bp': other_bp,
        'k': None, 'signature': 'ARCH_NONE', 'ltr_family': None,
        'arch': arch_string(ordered),
        'n_ltr_termini': 0,
        'int_gaps': int_consensus_gaps(ordered),
    }

    ltrs = [f for f in ordered if is_ltr_family(f['name'])]
    if not ltrs:
        return out
    fam = max(ltrs, key=lambda f: f['bp'])['name']
    L = ltr_len(fam)
    out['ltr_family'] = fam

    has_int = int_bp >= cfg['min_int_bp']
    first, last = ordered[0], ordered[-1]
    span_lo = min(f['qstart'] for f in ordered)
    span_hi = max(f['qend'] for f in ordered)
    first_terminal = is_ltr_family(first['name']) and (
        (first['qstart'] - span_lo if strand == '+' else span_hi - first['qend'])
        <= cfg['terminus_tol'])
    last_terminal = is_ltr_family(last['name']) and (
        (span_hi - last['qend'] if strand == '+' else last['qstart'] - span_lo)
        <= cfg['terminus_tol'])
    out['n_ltr_termini'] = int(first_terminal) + int(last_terminal)

    # Lone LTR, no internal region: a solo LTR allele.
    if not has_int and len(ltrs) == 1 and _close(ltr_bp, L, 4 * cfg['perm_tol']):
        out['signature'] = 'ARCH_SOLO'
        return out

    if not (has_int and first_terminal and last_terminal):
        return out

    tol = cfg['perm_tol']
    first_full = _close(first['cons_start'], 1, tol) and _close(first['cons_end'], L, tol)
    last_full = _close(last['cons_start'], 1, tol) and _close(last['cons_end'], L, tol)

    # Two whole LTRs: nothing was consumed by the alignment, so the SV carries
    # a complete provirus and the reference must have been empty.
    if first_full and last_full:
        out['signature'] = 'ARCH_2LTR'
        return out

    # One LTR split across the two termini. The 5' piece must run to the end of
    # the consensus; the 3' piece supplies the rest. Where the 3' piece was won
    # by SVA its consensus coordinates are unknown, so its length carries the
    # test instead -- that is the chr7 and chr12 case.
    if _close(first['cons_end'], L, tol) and first['cons_start'] and first['cons_start'] > 1 + tol:
        k = first['cons_start'] - 1
        tail_ok = _close(last['bp'], k, tol)
        if last['cons_start'] is not None:
            tail_ok = tail_ok or (_close(last['cons_start'], 1, tol)
                                  and _close(last['cons_end'], k, tol))
        if tail_ok and _close(ltr_bp, L, 2 * tol):
            out['signature'] = 'ARCH_PERM'
            out['k'] = k
    return out


def int_consensus_gaps(ordered):
    """Internal-region consensus gaps, e.g. the recurrent 291 bp HERVK-int
    5535-5825 deletion seen at chr3/chr5/chr12/chr19."""
    ints = [f for f in ordered
            if is_int_family(f['name']) and f['cons_start'] is not None]
    if len(ints) < 2:
        return ''
    ints.sort(key=lambda f: f['cons_start'])
    gaps = []
    for a, b in zip(ints, ints[1:]):
        if b['cons_start'] - a['cons_end'] > 1:
            gaps.append(f"{a['cons_end'] + 1}-{b['cons_start'] - 1}")
    return ','.join(gaps)


def arch_string(ordered):
    """Compact 5'->3' architecture, e.g. LTR:575-968/INT:1-7536/LTR:1-574."""
    parts = []
    for f in ordered:
        if is_ltr_family(f['name']):
            tag = 'LTR'
        elif is_int_family(f['name']):
            tag = 'INT'
        else:
            tag = f['name']
        if f['cons_start'] is None:
            parts.append(f"{tag}:~{f['bp']}bp")
        else:
            parts.append(f"{tag}:{f['cons_start']}-{f['cons_end']}")
    return '/'.join(parts) if parts else 'NONE'


def analyse(rm_hits, cfg, wanted=None):
    """Return {sv_id: architecture dict} for every SV in the RM tables."""
    out = {}
    for sv_id, hits in rm_hits.items():
        if wanted is not None and sv_id not in wanted:
            continue
        frags = reassign_sine_r(tile_hits(hits, cfg), cfg)
        out[sv_id] = architecture(frags, None, cfg)
    return out


COLUMNS = ['id', 'signature', 'k', 'ltr_family', 'n_ltr_termini', 'n_frag',
           'ltr_bp', 'int_bp', 'other_bp', 'strand', 'int_gaps', 'arch']


def write_tsv(results, path):
    with open(path, 'w') as fh:
        fh.write('\t'.join(COLUMNS) + '\n')
        for sv_id in sorted(results):
            r = results[sv_id]
            row = [sv_id, r['signature'],
                   '' if r['k'] is None else str(r['k']),
                   r['ltr_family'] or '', str(r['n_ltr_termini']),
                   str(r['n_frag']), str(int(r['ltr_bp'])),
                   str(int(r['int_bp'])), str(int(r['other_bp'])),
                   r['strand'], r['int_gaps'], r['arch']]
            fh.write('\t'.join(row) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--rm-out', nargs='+', required=True,
                    help='RepeatMasker .out files, or repeatmasker_dir paths/globs')
    ap.add_argument('--out', required=True, help='output architecture TSV')
    ap.add_argument('--ids', help='optional file of SV IDs to restrict to')
    args = ap.parse_args()

    paths = expand_rm_paths(args.rm_out)
    if not paths:
        sys.exit(f'ERROR: no RepeatMasker .out files found in {args.rm_out}')

    wanted = None
    if args.ids:
        with open(args.ids) as fh:
            wanted = {ln.strip() for ln in fh if ln.strip()}

    results = analyse(parse_rm_out(paths), dict(DEFAULTS), wanted)
    write_tsv(results, args.out)
    sys.stderr.write(f'hervk_arch: {len(results)} SVs from {len(paths)} RM table(s)\n')


if __name__ == '__main__':
    main()
