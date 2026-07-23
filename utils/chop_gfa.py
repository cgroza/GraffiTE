#!/usr/bin/env python3
"""
chop_gfa.py  –  Chop every GFA node to unit length (1 base per node).

Produces:
  1. A new GFA where every segment has sequence length 1.
  2. A TSV translation table: (original_id, offset) -> new_id.

Supports GFA 1 (S / L / P lines).  Segments with unknown sequence ('*')
are kept as single placeholder nodes.  All other record types are passed
through verbatim.

Usage
-----
    python chop_gfa.py  input.gfa  output.gfa  [-t translation.tsv]
"""

import re
import sys
import argparse
from collections import OrderedDict


# ── helpers ──────────────────────────────────────────────────────────────────

def natural_key(s: str):
    """Numeric-aware sort key ("2" sorts before "10")."""
    return [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', str(s))]


def sub_id(orig: str, offset: int) -> str:
    """Canonical name for the chopped sub-node at *offset* inside *orig*."""
    return f"{orig}_{offset}"


def seg_len(seq: str) -> int:
    """Effective length; unknown ('*') or empty sequences count as 1."""
    return len(seq) if seq not in ('', '*') else 1


# ── parsing ───────────────────────────────────────────────────────────────────

def parse_gfa(path: str):
    """
    Parse a GFA 1 file.

    Returns
    -------
    headers  : list[str]    – raw H lines
    segments : OrderedDict  – id -> {'seq': str, 'tags': list[str]}
    links    : list[dict]   – L records
    paths    : list[dict]   – P records
    other    : list[str]    – every other record type (passed through)
    """
    headers, other = [], []
    segments = OrderedDict()
    links, paths = [], []

    with open(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip('\n')
            if not line:
                continue
            cols = line.split('\t')
            rt = cols[0]

            if rt == 'H':
                headers.append(line)

            elif rt == 'S':
                if len(cols) < 3:
                    sys.exit(f"ERROR line {lineno}: malformed S record: {line!r}")
                segments[cols[1]] = {'seq': cols[2], 'tags': cols[3:]}

            elif rt == 'L':
                if len(cols) < 6:
                    sys.exit(f"ERROR line {lineno}: malformed L record: {line!r}")
                links.append({
                    'from_id':     cols[1],
                    'from_orient': cols[2],
                    'to_id':       cols[3],
                    'to_orient':   cols[4],
                    'overlap':     cols[5],
                    'tags':        cols[6:],
                })

            elif rt == 'P':
                if len(cols) < 3:
                    sys.exit(f"ERROR line {lineno}: malformed P record: {line!r}")
                paths.append({
                    'name':     cols[1],
                    'segments': cols[2],
                    'overlaps': cols[3] if len(cols) > 3 else '*',
                    'tags':     cols[4:],
                })

            else:
                other.append(line)

    return headers, segments, links, paths, other


# ── chopping ──────────────────────────────────────────────────────────────────

def chop_segments(segments: OrderedDict):
    """
    Split every segment into one node per base.

    Returns
    -------
    new_segs    : OrderedDict  new_id -> {'seq': base, 'tags': []}
    int_links   : list[dict]   internal forward-chain links
    translation : dict         (original_id, offset) -> new_id
    """
    new_segs = OrderedDict()
    int_links = []
    translation = {}

    for sid, data in segments.items():
        seq = data['seq']

        if seq in ('', '*'):
            # Unknown sequence: single placeholder node
            nid = sub_id(sid, 0)
            new_segs[nid] = {'seq': '*', 'tags': []}
            translation[(sid, 0)] = nid
            continue

        n = len(seq)
        for i, base in enumerate(seq):
            nid = sub_id(sid, i)
            new_segs[nid] = {'seq': base, 'tags': []}
            translation[(sid, i)] = nid

        # internal chain: sid_0 + -> sid_1 + -> ... -> sid_{n-1} +
        for i in range(n - 1):
            int_links.append({
                'from_id':     sub_id(sid, i),
                'from_orient': '+',
                'to_id':       sub_id(sid, i + 1),
                'to_orient':   '+',
                'overlap':     '0M',
                'tags':        [],
            })

    return new_segs, int_links, translation


# ── link update ───────────────────────────────────────────────────────────────

def update_links(links: list, segments: dict) -> list:
    """
    Re-attach every external link to the correct terminal sub-node.

    For a link  A <oA>  B <oB>  (with lengths n_A, n_B):
      oA '+' -> last  sub-node of A  (A_{n_A-1}, orient '+')
      oA '-' -> first sub-node of A  (A_0,       orient '-')
      oB '+' -> first sub-node of B  (B_0,       orient '+')
      oB '-' -> last  sub-node of B  (B_{n_B-1}, orient '-')
    """
    result = []
    for lk in links:
        fid, fo = lk['from_id'], lk['from_orient']
        tid, to = lk['to_id'],   lk['to_orient']

        fs = segments.get(fid)
        ts = segments.get(tid)
        if fs is None or ts is None:
            sys.stderr.write(
                f"WARNING: unknown node in link {fid}{fo}->{tid}{to}; skipped\n")
            continue

        n_f = seg_len(fs['seq'])
        n_t = seg_len(ts['seq'])

        new_from = sub_id(fid, n_f - 1) if fo == '+' else sub_id(fid, 0)
        new_to   = sub_id(tid, 0)       if to == '+' else sub_id(tid, n_t - 1)

        result.append({
            'from_id':     new_from,
            'from_orient': fo,
            'to_id':       new_to,
            'to_orient':   to,
            'overlap':     '0M',
            'tags':        lk['tags'],
        })

    return result


# ── path update ───────────────────────────────────────────────────────────────

def update_paths(paths: list, segments: dict) -> list:
    """
    Expand each path segment step into unit-length sub-node steps.

    X+  expands to  X_0+, X_1+, ..., X_{n-1}+
    X-  expands to  X_{n-1}-, X_{n-2}-, ..., X_0-
    """
    result = []
    for path in paths:
        expanded = []
        for item in path['segments'].split(','):
            item = item.strip()
            if item and item[-1] in '+-':
                sid, orient = item[:-1], item[-1]
            else:
                sid, orient = item, '+'

            seg = segments.get(sid)
            if seg is None:
                sys.stderr.write(
                    f"WARNING: unknown segment {sid!r} in path "
                    f"{path['name']!r}\n")
                continue

            n = seg_len(seg['seq'])
            rng = range(n) if orient == '+' else range(n - 1, -1, -1)
            for i in rng:
                expanded.append(f"{sub_id(sid, i)}{orient}")

        n_steps = len(expanded)
        overlaps = ','.join(['0M'] * (n_steps - 1)) if n_steps > 1 else '*'
        result.append({
            'name':     path['name'],
            'segments': ','.join(expanded),
            'overlaps': overlaps,
            'tags':     path['tags'],
        })

    return result


# ── output ────────────────────────────────────────────────────────────────────

def write_gfa(out_path, headers, new_segs, all_links, new_paths, other):
    with open(out_path, 'w') as fh:
        for line in headers:
            fh.write(line + '\n')
        for sid, d in new_segs.items():
            fh.write('\t'.join(['S', sid, d['seq']] + d['tags']) + '\n')
        for lk in all_links:
            fh.write('\t'.join([
                'L',
                lk['from_id'], lk['from_orient'],
                lk['to_id'],   lk['to_orient'],
                lk['overlap'],
            ] + lk['tags']) + '\n')
        for p in new_paths:
            fh.write('\t'.join(
                ['P', p['name'], p['segments'], p['overlaps']] + p['tags']
            ) + '\n')


def write_translation(out_path, translation: dict, segments: dict):
    """Write the (original_id, offset) -> new_id mapping as a TSV."""
    rows = sorted(translation.items(),
                  key=lambda kv: (natural_key(kv[0][0]), kv[0][1]))
    with open(out_path, 'w') as fh:
        fh.write('original_id\toriginal_length\toffset\tnew_id\n')
        for (oid, off), nid in rows:
            olen = seg_len(segments[oid]['seq'])
            fh.write(f"{oid}\t{olen}\t{off}\t{nid}\n")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description='Chop every GFA segment to unit length (1 base per node).')
    ap.add_argument('input',  help='Input GFA file')
    ap.add_argument('output', help='Output (chopped) GFA file')
    ap.add_argument('-t', '--translation', default='translation.tsv',
                    help='Translation table (default: translation.tsv)')
    args = ap.parse_args()

    print(f"[1/5] Parsing {args.input} ...")
    headers, segments, links, paths, other = parse_gfa(args.input)
    total_bp = sum(seg_len(d['seq']) for d in segments.values())
    print(f"      {len(segments):,} segments · {total_bp:,} bp · "
          f"{len(links):,} links · {len(paths):,} paths")

    print("[2/5] Chopping segments ...")
    new_segs, int_links, translation = chop_segments(segments)

    print("[3/5] Re-wiring external links ...")
    ext_links = update_links(links, segments)

    print("[4/5] Expanding paths ...")
    new_paths = update_paths(paths, segments)

    all_links = int_links + ext_links

    print("[5/5] Writing output ...")
    write_gfa(args.output, headers, new_segs, all_links, new_paths, other)
    write_translation(args.translation, translation, segments)

    print(f"\n  New segments    : {len(new_segs):>10,}")
    print(f"  Internal links  : {len(int_links):>10,}")
    print(f"  External links  : {len(ext_links):>10,}")
    print(f"  GFA  -> {args.output}")
    print(f"  Table-> {args.translation}")


if __name__ == '__main__':
    main()
