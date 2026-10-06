#!/usr/bin/env python3
"""What the human_graphaligner and human_precomputed cells have to show.

Both repeat the spine's discovery with --human and change only the genotyping
back end. hervk_reconcile used to accept giraffe alone, so every other back end
stopped after genotyping (issue #100); it now accepts the three that genotype
through vg call. Each cell has to show that the reconciler ran on its back end
and consolidated what it consolidates on the spine:

  precomputed   vg call reads the spine's own graph and packs, so the
                consolidated VCF must hold the spine's records, INFO and
                genotypes included, and the report must be the spine's line for
                line. The records may come in another order: they follow the
                ##contig lines of the genotyped VCF, and this cell's header has
                listed the contigs in a different order from the spine's.
  graphaligner  the long reads come off h2 and h4, which carry the same sites
                as h1 and h3, the haplotypes behind the spine's S1 and S2. The
                records and the HERV-K fields discovery sets must be the
                spine's, and at every HERV-K locus the genotypes and the dosage
                counts must match the spine's, L1 against S1 and L2 against S2.
                Other genotypes are listed where they differ but do not fail
                the cell: a difference there says something about
                GraphAligner, not about the consolidation.

usage: check_human_backends.py <spine run> <cell run> <graphaligner|precomputed> <work dir>
"""
import csv
import glob
import gzip
import os
import re
import sys

OUT = '4_Genotyping'
VCF = 'GraffiTE.merged.genotypes.human.vcf.gz'
FILES = (VCF, 'hervk_unconsolidated_records.vcf', 'hervk_reconciliation_report.md')
# Set from the genotypes, so they can differ between back ends that call
# different genotypes; every other HERVK_* field comes from discovery.
FROM_GENOTYPES = {'HERVK_AC', 'HERVK_AN', 'HERVK_N_RESOLVED', 'HERVK_N_PARTIAL',
                  'HERVK_N_PLOIDY_EXCEEDED', 'HERVK_DISC_CONCORDANT', 'HERVK_GT_MASKED',
                  'HERVK_MEMBERS_MASKED', 'HERVK_DISC_PLOIDY_MISMATCH'}


def read_vcf(path):
    """(sample names, [record fields]) below the ## header."""
    samples, records = [], []
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            f = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                samples = f[9:]
            else:
                records.append(f)
    return samples, records


def info(field):
    return dict(kv.split('=', 1) if '=' in kv else (kv, True) for kv in field.split(';') if kv)


def reconcile_task(run, work):
    """The hervk_reconcile row of the run's trace and the --genotyper its script passed."""
    for r in csv.DictReader(open(os.path.join(run, 'nextflow_trace.txt')), delimiter='\t'):
        if r['name'].startswith('hervk_reconcile'):
            seen = None
            for cmd in glob.glob(os.path.join(work, r['hash'] + '*', '.command.sh')):
                m = re.search(r'--genotyper\s+(\S+)', open(cmd).read())
                seen = m.group(1) if m else None
            return r, seen
    return None, None


def main():
    if len(sys.argv) != 5 or sys.argv[3] not in ('graphaligner', 'precomputed'):
        sys.exit(__doc__)
    spine, run, mode, work = sys.argv[1:]
    fail = 0

    def check(label, good, detail=''):
        nonlocal fail
        print(f"  [{' ok ' if good else 'FAIL'}] {label}" + (f": {detail}" if detail and not good else ''))
        fail += not good

    row, seen = reconcile_task(run, work)
    check("hervk_reconcile ran and exited 0", row is not None and row['status'] == 'COMPLETED',
          f"trace says {row['status']} (exit {row['exit']})" if row else "not in the trace")
    check(f"its script passed --genotyper {mode}", seen == mode, f"it passed {seen!r}")
    missing = [f for f in FILES if not os.path.isfile(os.path.join(run, OUT, f))]
    check("it published the consolidated VCF, the archive and the report", not missing,
          f"missing {', '.join(missing)}")
    if row is None or missing:
        print(f"human_{mode} fail={fail}")
        sys.exit(1)

    s_samples, s_recs = read_vcf(os.path.join(spine, OUT, VCF))
    c_samples, c_recs = read_vcf(os.path.join(run, OUT, VCF))
    loci = [r for r in s_recs if 'HERVK_LOCUS' in info(r[7])]
    check("the spine has HERV-K locus records to compare against", bool(loci), "none in its VCF")

    if mode == 'precomputed':
        check(f"the same samples as the spine ({', '.join(s_samples)})", c_samples == s_samples,
              f"{', '.join(c_samples)}")
        a_set, b_set = set(map(tuple, s_recs)), set(map(tuple, c_recs))
        only_s = [r[2] for r in s_recs if tuple(r) not in b_set]
        only_c = [r[2] for r in c_recs if tuple(r) not in a_set]
        same = sorted(map(tuple, s_recs)) == sorted(map(tuple, c_recs))
        check(f"the spine's {len(s_recs)} records, INFO and genotypes included, in any order", same,
              f"{len(c_recs)} records here; differing records: {only_s[:3]} in the spine, {only_c[:3]} here")
        if same and s_recs != c_recs:
            order = lambda recs: ' '.join(dict.fromkeys(r[0] for r in recs))
            print(f"  [info] the same records in another contig order: spine {order(s_recs)}, here {order(c_recs)}")
        a = open(os.path.join(spine, OUT, FILES[2])).read().splitlines()
        b = open(os.path.join(run, OUT, FILES[2])).read().splitlines()
        first = next((i for i, (x, y) in enumerate(zip(a, b)) if x != y), None)
        check("the spine's reconciliation report", a == b,
              f"line {first + 1}: {b[first]!r}" if first is not None else f"{len(b)} lines against {len(a)}")
    else:
        key = lambda r: tuple(r[:5])
        check(f"the spine's {len(s_recs)} records (CHROM, POS, ID, REF, ALT)",
              sorted(map(key, s_recs)) == sorted(map(key, c_recs)),
              f"{len(c_recs)} records here")
        pairs = list(zip(c_samples, s_samples))
        print(f"  [info] samples compared: {', '.join(f'{c} against {s}' for c, s in pairs)}")
        theirs = {key(r): r for r in c_recs}
        struct, locus_gt, other_gt = [], [], []
        for r in s_recs:
            c = theirs.get(key(r))
            if c is None:
                continue
            si, ci = info(r[7]), info(c[7])
            hk = lambda d: {k: v for k, v in d.items() if k.startswith('HERVK_') and k not in FROM_GENOTYPES}
            if hk(si) != hk(ci):
                struct.append(r[2])
            gts = [(x.split(':')[0], y.split(':')[0]) for x, y in zip(r[9:], c[9:])]
            if 'HERVK_LOCUS' in si:
                counts = [(si.get(k), ci.get(k)) for k in ('HERVK_AC', 'HERVK_AN')]
                if any(x != y for x, y in gts) or any(x != y for x, y in counts):
                    locus_gt.append(f"{r[2]} GT {[x for x, _ in gts]} vs {[y for _, y in gts]}, "
                                    f"AC/AN {[x for x, _ in counts]} vs {[y for _, y in counts]}")
            elif any(x != y for x, y in gts):
                other_gt.append(f"{r[2]} {[x for x, _ in gts]} vs {[y for _, y in gts]}")
        check("the HERV-K fields discovery sets, on every record", not struct,
              f"{len(struct)} differ, first {struct[0]}" if struct else '')
        check(f"genotypes and HERVK_AC/HERVK_AN at the spine's {len(loci)} HERV-K locus records",
              not locus_gt, '; '.join(locus_gt))
        for d in other_gt:
            print(f"  [info] genotype differs outside the HERV-K loci: {d}")

    print(f"human_{mode} fail={fail}")
    sys.exit(1 if fail else 0)


if __name__ == '__main__':
    main()
