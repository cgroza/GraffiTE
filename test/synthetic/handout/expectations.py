#!/usr/bin/env python3
"""Freeze an OBSERVED.tsv, and check a later one against it.

Only run the freeze once calibration is green. A blessed wrong expectation reads
exactly like a right one.

  expectations.py freeze  OBSERVED.tsv EXPECTED.tsv
  expectations.py check   OBSERVED.tsv EXPECTED.tsv [--log ASSERTIONS.log]

Exit codes follow docs/design-notes/synthetic-test-set.md:
  0  everything matches
  1  a structural invariant broke: a site stopped reaching a file it reached
     before, a subset lost a member, the span ladder changed shape
  2  a tool heuristic moved: a score, a reported edge, a link ID, an annotation
     string. Review and re-freeze as a committed diff.

The split is by column, not by size of change. Which file a record reaches is
the pipeline's decision and ours to assert. What RepeatMasker scored it is
RepeatMasker's, and a container rebuild is allowed to move it, but not silently.
"""

import argparse
import hashlib
import os
import re
import sys
import time

STRUCTURAL = [
    'in_svs', 'in_prefilter', 'in_postfilter', 'in_pangenome',
    'in_trusted', 'in_human', 'n_hits',
]
HEURISTIC = [
    'repeat_ids', 'matching_classes', 'fragmts', 'strands', 'L1_5PINV',
    'total_repeat_span', 'ULTRA_TR_span', 'TSD', 'polyA',
    'rm_n_rows', 'rm_link_ids', 'rm_best_sw', 'rm_repeat_ids',
    'rm_target_start', 'rm_target_end', 'rm_strands', 'vcf_svlen',
]
KEY = 'name'


def read(path):
    meta = {}
    with open(path) as fh:
        line = fh.readline()
        while line.startswith('#'):
            k, _, v = line[1:].strip().partition(':')
            meta[k.strip()] = v.strip()
            line = fh.readline()
        hdr = line.rstrip('\n').split('\t')
        rows = [dict(zip(hdr, l.rstrip('\n').split('\t'))) for l in fh if l.strip()]
    return hdr, rows, meta


def provenance(observed):
    """What the measurement was of. Without this a freeze and its own source are
    indistinguishable from a freeze and a fresh run, and `check` reports a match
    either way. That is how a matrix that failed every cell at launch, leaving
    the previous run's output in place, produced a green freeze."""
    here = os.path.dirname(os.path.abspath(observed))
    out = {'observed': os.path.abspath(observed),
           'observed_sha256': sha(observed),
           'observed_mtime': stamp(observed)}
    env = os.path.join(here, 'INPUTS.env')
    workdir = image = None
    if os.path.exists(env):
        for line in open(env):
            m = re.match(r'\s*(WORKDIR|GRAFFITE_SIF)="([^"]*)"', line)
            if m and m.group(1) == 'WORKDIR':
                workdir = m.group(2)
            elif m:
                image = m.group(2)
    if workdir:
        mani = os.path.join(workdir, 'build', 'MANIFEST.sha256')
        out['build_manifest_sha256'] = sha(mani)
        run = os.path.join(workdir, 'runs', 'spine', '3_TSD_search', 'pangenome.vcf')
        out['spine_pangenome_vcf'] = run if os.path.exists(run) else '(absent)'
        out['spine_mtime'] = stamp(run)
    if image:
        out['image'] = image
        out['image_sha256'] = sha(image)
    return out


def sha(path):
    if not path or not os.path.exists(path):
        return '(absent)'
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()[:24]


def stamp(path):
    if not path or not os.path.exists(path):
        return '(absent)'
    return time.strftime('%Y-%m-%dT%H:%M:%S', time.localtime(os.path.getmtime(path)))


def freeze(src, dst):
    hdr, rows, _ = read(src)
    cols = [KEY, 'kind', 'contig', 'pos', 'truth_svlen'] + \
           [c for c in STRUCTURAL + HEURISTIC if c in hdr]
    prov = provenance(src)
    with open(dst, 'w') as fh:
        for k, v in prov.items():
            fh.write('# %s: %s\n' % (k, v))
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(r.get(c, '') for c in cols) + '\n')
    print('froze %d sites, %d columns into %s' % (len(rows), len(cols), dst))
    for k, v in prov.items():
        print('  %-22s %s' % (k, v))
    print('commit this file: a re-freeze has to be reviewable in git log.')
    return 0


def check(obs_path, exp_path, log_path):
    ohdr, orows, _ = read(obs_path)
    ehdr, erows, emeta = read(exp_path)
    now = provenance(obs_path)
    same = (emeta.get('observed_sha256') not in (None, '', '(absent)')
            and emeta.get('observed_sha256') == now.get('observed_sha256'))
    obs = {r[KEY]: r for r in orows}
    out, structural, heuristic = [], [], []

    missing = [r[KEY] for r in erows if r[KEY] not in obs]
    extra = [k for k in obs if k not in {r[KEY] for r in erows}]
    for m in missing:
        structural.append('%s: site is in EXPECTED.tsv and not in OBSERVED.tsv' % m)
    for m in extra:
        structural.append('%s: site is in OBSERVED.tsv and not in EXPECTED.tsv' % m)

    for e in erows:
        o = obs.get(e[KEY])
        if not o:
            continue
        for c in ehdr:
            if c in (KEY, 'kind', 'contig', 'pos', 'truth_svlen'):
                continue
            want, got = e.get(c, ''), o.get(c, '')
            if want == got:
                continue
            line = '%-28s %-18s expected %-22s got %s' % (e[KEY], c, want, got)
            (structural if c in STRUCTURAL else heuristic).append(line)

    if same:
        out.append('!! this OBSERVED.tsv is byte-identical to the one EXPECTED.tsv')
        out.append('!! was frozen from, so the comparison below cannot fail. Re-run')
        out.append('!! the pipeline and measure again before trusting a match.')
        out.append('')
    for k in ('build_manifest_sha256', 'image_sha256', 'spine_mtime'):
        if emeta.get(k) and emeta[k] != now.get(k):
            out.append('note: %s moved since the freeze: %s -> %s'
                       % (k, emeta[k], now.get(k)))
    out.append('observed : %s' % obs_path)
    out.append('expected : %s' % exp_path)
    out.append('sites    : %d expected, %d observed' % (len(erows), len(orows)))
    out.append('')
    out.append('== structural (exit 1) == %d' % len(structural))
    out += ['  ' + s for s in structural] or ['  none']
    out.append('')
    out.append('== tool heuristics moved (exit 2) == %d' % len(heuristic))
    out += ['  ' + s for s in heuristic] or ['  none']
    out.append('')
    rc = 1 if structural else (2 if heuristic else 0)
    out.append('verdict  : %s' % {0: 'match', 1: 'STRUCTURAL BREAK',
                                  2: 'heuristics moved, review and re-freeze'}[rc])
    text = '\n'.join(out) + '\n'
    if log_path:
        open(log_path, 'w').write(text)
    sys.stdout.write(text)
    return rc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['freeze', 'check'])
    ap.add_argument('observed')
    ap.add_argument('expected')
    ap.add_argument('--log', default=None)
    a = ap.parse_args()
    if a.mode == 'freeze':
        return freeze(a.observed, a.expected)
    return check(a.observed, a.expected, a.log)


if __name__ == '__main__':
    sys.exit(main())
