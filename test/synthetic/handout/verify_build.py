#!/usr/bin/env python3
"""Check the generated test set against itself before any of it reaches GraffiTE.

The design note asks for the fixture to be asserted before the pipeline is, so
that a fixture that drifted does not read as a classifier regression. Everything
here is arithmetic on the files the generator wrote; nothing runs the pipeline.

  verify_build.py <WORKDIR>/build

Checks:
  1. every haplotype carrying a site holds the planted sequence at the reference
     position truth.tsv records, and every non-carrier does not
  2. the anchor base of each --svs record is the reference base at its POS
  3. each --svs ALT is the anchor base plus the planted insert
  4. sample names are distinct across --assemblies and --svs
  5. every variant ID is at most 50 characters (bin/repmask_vcf.sh:43 exits 1 above that)
  6. contigs are in ASCII-lexicographic order, which vg autoindex requires
"""

import gzip
import pathlib
import sys


def read_fasta(path):
    seqs, name, buf = {}, None, []
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(buf)
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = ''.join(buf)
    return seqs


def main():
    B = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else 'build')
    bad = []
    note = []

    rows = [l.rstrip('\n').split('\t') for l in open(B / 'truth.tsv')]
    hdr, sites = rows[0], [dict(zip(rows[0], r)) for r in rows[1:]]

    ref = read_fasta(B / 'ref' / 'synth.fa')
    haps = {h: read_fasta(B / 'hap' / f'{h}.fa') for h in ['h1', 'h2', 'h3', 'h4']}

    # the inserted sequence is not in truth.tsv, so recover it from the --svs fixture
    inserts = {}
    for hap in ['h1', 'h2', 'h3', 'h4']:
        p = B / 'vcf' / f'svs_{hap}.vcf.gz'
        if not p.exists():
            continue
        with gzip.open(p, 'rt') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                f = line.rstrip('\n').split('\t')
                inserts[f[2]] = (f[0], int(f[1]), f[3], f[4])

    # 1 + 2 + 3
    #
    # `pos` in truth.tsv is a REFERENCE coordinate. A haplotype carrying several
    # sites on one contig holds the nth one at pos plus the total length of the
    # inserts placed before it, which is the offset checked here. The reference
    # position is what svim-asm reports and what the --svs fixture must agree with.
    by_contig = {}
    for s_ in sites:
        by_contig.setdefault(s_['contig'], []).append(s_)

    for s_ in sites:
        name, contig, pos = s_['name'], s_['contig'], int(s_['pos'])
        carriers = s_['carriers'].split('|')
        if name not in inserts:
            bad.append(f'{name}: no record in any svs_h*.vcf.gz')
            continue
        c, p, refb, alt = inserts[name]
        if (c, p) != (contig, pos):
            bad.append(f'{name}: svs record at {c}:{p}, truth.tsv says {contig}:{pos}')
        if refb != ref[contig][pos - 1]:
            bad.append(f'{name}: anchor base {refb} but reference has '
                       f'{ref[contig][pos - 1]} at {contig}:{pos}')
        ins = alt[1:]
        if len(ins) != int(s_['svlen']):
            bad.append(f'{name}: ALT carries {len(ins)} bp, truth.tsv says {s_["svlen"]}')
        for h in ['h1', 'h2', 'h3', 'h4']:
            shift = sum(len(inserts[o['name']][3]) - 1
                        for o in by_contig[contig]
                        if int(o['pos']) < pos and h in o['carriers'].split('|')
                        and o['name'] in inserts)
            off = pos + shift
            got = haps[h][contig][off:off + len(ins)]
            if h in carriers and got != ins:
                bad.append(f'{name}: {h} carries it but {contig} offset {off} '
                           f'(reference {pos} + {shift}) does not hold the insert')
            if h not in carriers and got == ins:
                bad.append(f'{name}: {h} is not a carrier yet holds the insert')

    # every haplotype must be the reference with exactly its carrier inserts
    # spliced in at their reference positions, and nothing else
    for h in ['h1', 'h2', 'h3', 'h4']:
        for contig in ref:
            want = ref[contig]
            for s_ in sorted(by_contig.get(contig, []), key=lambda r: -int(r['pos'])):
                if h in s_['carriers'].split('|') and s_['name'] in inserts:
                    q = int(s_['pos'])
                    want = want[:q] + inserts[s_['name']][3][1:] + want[q:]
            if haps[h][contig] != want:
                first = next((i for i, (a, b) in enumerate(zip(haps[h][contig], want))
                              if a != b), min(len(haps[h][contig]), len(want)))
                bad.append(f'{h} {contig}: {len(haps[h][contig])} bp, expected '
                           f'{len(want)} bp, first difference at {first}')

    # 4
    asm = [l.split(',')[0] for l in (B / 'assemblies.csv').read_text().splitlines()[1:] if l]
    svs_samples = []
    for hap in ['h1', 'h2', 'h3', 'h4']:
        p = B / 'vcf' / f'svs_{hap}.vcf.gz'
        if not p.exists():
            continue
        with gzip.open(p, 'rt') as fh:
            for line in fh:
                if line.startswith('#CHROM'):
                    svs_samples += line.rstrip('\n').split('\t')[9:]
                    break
    clash = sorted(set(asm) & set(svs_samples))
    if clash:
        bad.append(f'sample names shared by --assemblies and --svs: {clash} '
                   f'(module/main.nf:223 runs bcftools merge with no --force-samples)')

    # 5
    long_ids = [i for i in inserts if len(i) > 50]
    if long_ids:
        bad.append(f'variant IDs over 50 characters: {long_ids}')

    # 6
    order = list(ref)
    if order != sorted(order):
        bad.append(f'contigs are not in ASCII-lexicographic order: {order}')

    note.append(f'sites            : {len(sites)}')
    note.append(f'contigs          : {len(ref)}, {sum(len(v) for v in ref.values())} bp reference')
    note.append(f'--assemblies     : {asm}')
    note.append(f'--svs samples    : {svs_samples}')
    for h in haps:
        note.append(f'{h} length        : {sum(len(v) for v in haps[h].values())} bp')
    print('\n'.join(note))
    print()
    if bad:
        print(f'verify_build: {len(bad)} problem(s)')
        for b in bad:
            print('  [FAIL]', b)
        return 1
    print('verify_build: the test set is internally consistent')
    return 0


if __name__ == '__main__':
    sys.exit(main())
