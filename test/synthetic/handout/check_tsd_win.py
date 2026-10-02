#!/usr/bin/env python3
"""Compare a --tsd_win run's TSD search against the spine's, which used 30.

The window sets the flank width prepTSD.sh cuts and the junction column
TSD_Match_v2.sh scores against, and the two have to agree. Flanks of 40 scored
against a junction at 30 made the matcher pick other hits: all 29 SVs of the
test set changed, and PASS went from 27 to 16. When both move, the offsets and
scores are measured from the junction, so each SV should come back with the
same verdict as in the spine, and each PASS with the same TSD, offsets and score.

A FAIL can report a different best hit. Since 170c788 the matcher ranks a hit
under 6 bp that is off the junction below any longer one, so a wider window can
bring in a longer, farther hit and report it instead. Both fail, and neither
reaches INFO/TSD; the script lists them without counting them as a failure.

usage: check_tsd_win.py <spine run dir> <run dir> <window> [spine window]
"""
import os
import sys


def fragment_widths(path):
    """Lengths of the [flank][SV] fragments TSD_Match_v2.sh logs, one per end."""
    widths = set()
    with open(path) as fh:
        lines = iter(fh)
        for line in lines:
            if line.startswith('>L|5P_end') or line.startswith('>R|3P_end'):
                widths.add(len(next(lines, '').rstrip('\n')))
    return widths


def summary(path):
    """SV id -> (verdict, L_TSD, R_TSD, score, four junction offsets)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 16:
                continue
            if f[-2] == 'no_hit':
                out[f[0]] = (f[-1], 'no_hit', 'no_hit', 'NA', ())
            else:
                # f[13:17] are win-R_start, win-L_start, win-R_end, win-L_end
                out[f[0]] = (f[-1], f[-3], f[-2], f[-4], tuple(f[13:17]))
    return out


def vcf_tsd(path):
    """Record id -> INFO/TSD, or None when the record has none."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t', 8)
            tsd = None
            for kv in f[7].split(';'):
                if kv.startswith('TSD='):
                    tsd = kv[4:]
            out[f[2]] = tsd
    return out


def main():
    if len(sys.argv) not in (4, 5):
        sys.exit(__doc__)
    spine, run, win = sys.argv[1], sys.argv[2], int(sys.argv[3])
    spine_win = int(sys.argv[4]) if len(sys.argv) == 5 else 30
    fail = 0

    def check(label, good, detail=''):
        nonlocal fail
        print(f"  [{' ok ' if good else 'FAIL'}] {label}" + (f": {detail}" if detail and not good else ''))
        fail += not good

    s_dir, r_dir = (os.path.join(d, '3_TSD_search') for d in (spine, run))
    need = [os.path.join(d, f) for d in (s_dir, r_dir)
            for f in ('TSD_full_log.txt', 'TSD_summary.txt', 'pangenome.vcf')]
    missing = [p for p in need if not os.path.isfile(p)]
    check("TSD outputs present", not missing, ', '.join(missing))
    if missing:
        print(f"tsd_win fail={fail}")
        sys.exit(1)

    # Without these two, a run that ignored --tsd_win would match the spine
    # record for record and pass.
    w = fragment_widths(os.path.join(r_dir, 'TSD_full_log.txt'))
    check(f"fragments are 2 x {win} bp", w == {2 * win}, f"widths seen {sorted(w)}")
    w = fragment_widths(os.path.join(s_dir, 'TSD_full_log.txt'))
    check(f"spine fragments are 2 x {spine_win} bp", w == {2 * spine_win}, f"widths seen {sorted(w)}")

    a = summary(os.path.join(s_dir, 'TSD_summary.txt'))
    b = summary(os.path.join(r_dir, 'TSD_summary.txt'))
    check("same SVs searched", a.keys() == b.keys(),
          f"{len(a.keys() - b.keys())} only in spine, {len(b.keys() - a.keys())} only here")
    n_pass = sum(v[0] == 'PASS' for v in a.values())
    check("spine has TSDs to compare", n_pass > 0, "no PASS in the spine")
    both = a.keys() & b.keys()
    diff = sorted(k for k in both if a[k][0] != b[k][0])
    check(f"same verdict on all {len(both)} SVs", not diff,
          f"{len(diff)} differ, first {diff[0]}: {a[diff[0]]} vs {b[diff[0]]}" if diff else '')
    passed = [k for k in both if a[k][0] == 'PASS' or b[k][0] == 'PASS']
    diff = sorted(k for k in passed if a[k] != b[k])
    check(f"same TSD, offsets and score on all {len(passed)} PASS SVs", not diff,
          f"{len(diff)} differ, first {diff[0]}: {a[diff[0]]} vs {b[diff[0]]}" if diff else '')
    moved = sorted(k for k in both if k not in passed and a[k] != b[k])
    for k in moved:
        print(f"  [info] FAIL in both, best hit moved: {k}: {a[k][1]} -> {b[k][1]}")

    a = vcf_tsd(os.path.join(s_dir, 'pangenome.vcf'))
    b = vcf_tsd(os.path.join(r_dir, 'pangenome.vcf'))
    n_tsd = sum(v is not None for v in a.values())
    check("pangenome.vcf has the same records", a.keys() == b.keys(),
          f"{len(a)} in spine, {len(b)} here")
    check("spine pangenome.vcf carries TSDs", n_tsd > 0, "no record has INFO/TSD")
    diff = sorted(k for k in a.keys() & b.keys() if a[k] != b[k])
    check(f"same INFO/TSD on all {len(a.keys() & b.keys())} records ({n_tsd} with a TSD)",
          not diff, f"{len(diff)} differ, first {diff[0]}: {a[diff[0]]} vs {b[diff[0]]}" if diff else '')

    print(f"tsd_win fail={fail}")
    sys.exit(1 if fail else 0)


if __name__ == '__main__':
    main()
