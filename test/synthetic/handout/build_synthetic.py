#!/usr/bin/env python3
"""Build the synthetic GraffiTE test set: reference, library, haplotypes, inputs.

Deterministic. Same --seed and same Dfam source give the same bytes, and
MANIFEST.sha256 is the guard. Uses only the standard library; minimap2 and
samtools are needed for the BAM inputs and are skipped with a warning when
absent.

Consensus sequences come from the human_DFAM3.6.fasta inside the repository's
test/human_test_set.tar.gz, except L1HS: that library ships L1HS_5end and
L1HS_3end as separate entries and no full-length copy, and two library entries
would almost certainly get separate RepeatMasker link IDs, which is what
INFO/L1_5PINV keys on. The full-length L1HS here is built by joining the two
halves, so one library entry covers the whole element.

The HERV-K lengths are read out of bin/hervk_arch.py rather than restated: the
classifier compares observed hit lengths against them, so a library that
disagrees produces allele states that look like classifier bugs.

  build_synthetic.py --dfam human_DFAM3.6.fasta --out build/
"""

import argparse
import gzip
import hashlib
import os
import pathlib
import random
import re
import shutil
import subprocess
import sys

# ---------------------------------------------------------------- consensus

def read_fasta(path):
    seqs, name, buf = {}, None, []
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = ''.join(buf)
    return seqs


def pick(seqs, want):
    """Exact Dfam entry whose name before '#' is `want`."""
    for k in sorted(seqs):
        if k.split('#')[0] == want:
            return seqs[k].upper()
    sys.exit(f'build_synthetic.py: {want} not found in the Dfam library')


def hervk_constants(repo):
    """INT_CONSENSUS_LEN and LTR_CONSENSUS_LEN['LTR5_Hs'], from the classifier."""
    src = (repo / 'bin' / 'hervk_arch.py').read_text()
    m_int = re.search(r'^INT_CONSENSUS_LEN\s*=\s*(\d+)', src, re.M)
    m_ltr = re.search(r"^LTR_CONSENSUS_LEN\s*=\s*\{[^}]*'LTR5_Hs'\s*:\s*(\d+)", src, re.M)
    if not (m_int and m_ltr):
        sys.exit('build_synthetic.py: could not read the HERV-K constants from bin/hervk_arch.py')
    return int(m_int.group(1)), int(m_ltr.group(1))


# ---------------------------------------------------------------- sequence

COMP = str.maketrans('ACGTN', 'TGCAN')


def rc(s):
    return s.translate(COMP)[::-1]


def background(rng, n):
    """Order-1 Markov background. Uniform random sequence is unrealistically
    even and gives RepeatMasker an easier discrimination problem than a real
    genome does."""
    # Mildly AT-rich with a CpG deficit, which is what makes a background that
    # RepeatMasker has to work against rather than through.
    table = {
        'A': 'AAAAAAACCCGGGTTTTTTT', 'C': 'AAAAAACCCCGTTTTTTTTT',
        'G': 'AAAAAACCCCCGGGTTTTTT', 'T': 'AAAAAACCCGGGGTTTTTTT',
    }
    out = [rng.choice('ACGT')]
    for _ in range(n - 1):
        out.append(rng.choice(table[out[-1]]))
    return ''.join(out)


def mutate(rng, s, pct):
    """Substitute pct% of bases. Divergence is what decides whether
    RepeatMasker reports a copy at all, so it is a per-site knob."""
    if pct <= 0 or not s:
        return s
    out = list(s)
    for i in rng.sample(range(len(out)), max(1, int(len(out) * pct / 100.0))):
        out[i] = rng.choice([b for b in 'ACGT' if b != out[i]])
    return ''.join(out)


def tsd(rng, n=8):
    return ''.join(rng.choice('ACGT') for _ in range(n))


# ---------------------------------------------------------------- the design
#
# One row per planted variant. `kind` is what the site is for; the assertions
# on the cluster are written against these names, and truth.tsv carries them.
#
#   contig, name, kind, builder(rng, cons, const) -> inserted sequence
#
# Everything here is a NON-REFERENCE insertion unless kind starts with 'del_',
# in which case the element is in the reference and absent from the haplotype.

def plan_sites(cons, const):
    INT_LEN, LTR_LEN = const
    # L1HS here is L1HS_5end + L1HS_3end joined; the slices below must stay
    # inside it, and a silent empty slice would plant nothing.
    if len(cons['L1HS']) < 3038:
        sys.exit(f"build_synthetic.py: joined L1HS is {len(cons['L1HS'])} bp, "
                 "too short for the twin-priming slices")
    aluY, aluSx, sva_e, l1hs = cons['AluY'], cons['AluSx'], cons['SVA_E'], cons['L1HS']
    ltr5, hervk = cons['LTR5_Hs'], cons['HERVK']

    def ins(seq):
        return lambda r: seq(r) if callable(seq) else seq

    S = []
    # -- chr1: the annotation zone -------------------------------------------
    S += [
        ('chr1', 'A01_alu_tsd_polyA', 'canonical',
         lambda r: mutate(r, aluY, 5) + 'A' * 40, 8),
        ('chr1', 'A02_alu_no_polyA_no_tsd', 'trusted_noop_probe',
         lambda r: mutate(r, aluY, 5), 0),
        ('chr1', 'A03_alu_polyA_too_short', 'polyA_below_min',
         lambda r: mutate(r, aluY, 5) + 'A' * 6, 8),
        ('chr1', 'A04_alu_polyA_beyond_slack', 'polyA_beyond_slack',
         lambda r: mutate(r, aluY, 5) + 'A' * 15 + background(r, 10), 8),
        ('chr1', 'A05_alu_revcomp', 'strand_C_polyT',
         lambda r: rc(mutate(r, aluY, 5) + 'A' * 40), 8),
        ('chr1', 'A06_alusx_trusted_not_human', 'subfamily_gate',
         lambda r: mutate(r, aluSx, 8) + 'A' * 40, 8),
        ('chr1', 'A07_alu_full_length_clean', 'canonical_low_div',
         lambda r: mutate(r, aluY, 1) + 'A' * 30, 8),
        ('chr1', 'A08_alu_no_tsd', 'tsd_fail',
         lambda r: mutate(r, aluY, 5) + 'A' * 40, 0),
        # L1 twin priming: an inverted 5' piece then a forward 3' piece, both
        # from one consensus so ProcessRepeats can link them into one C+ group.
        ('chr1', 'A09_l1_twin_primed_Cplus', 'L1_5PINV_positive',
         lambda r: rc(mutate(r, l1hs[200:1400], 3)) + mutate(r, l1hs[2200:3038], 3) + 'A' * 30, 8),
        ('chr1', 'A10_l1_plusC_negative', 'L1_5PINV_negative',
         lambda r: mutate(r, l1hs[2200:3038], 3) + rc(mutate(r, l1hs[200:1400], 3)) + 'A' * 30, 8),
        ('chr1', 'A11_l1_full', 'L1_full',
         lambda r: mutate(r, l1hs[:3000], 4) + 'A' * 40, 8),
        # SVA_E VNTR-only. The window is 429..863 (435 bp wide) and the test in
        # annotate_vcf.R is strict on both sides; SVA_D's 433..688 band clears
        # the 250 bp floor by six bases, which is not a margin when
        # RepeatMasker decides the reported edges.
        ('chr1', 'A12_sva_vntr_only', 'SVA_VNTR_only',
         lambda r: mutate(r, sva_e[440:850], 4), 8),
        ('chr1', 'A13_sva_full', 'SVA_full',
         lambda r: mutate(r, sva_e, 4) + 'A' * 40, 8),
    ]
    # Eight-rung total_repeat_span ladder straddling --repeat_span_cutoff 0.80.
    # Each rung is a fixed TE core padded with unique sequence to hit a target
    # repeat fraction. Exactly one cut point is the assertion.
    for i, frac in enumerate([0.95, 0.90, 0.85, 0.82, 0.78, 0.70, 0.60, 0.45], start=14):
        core = 300
        total = int(round(core / frac))
        S.append(('chr1', f'A{i}_span_{int(frac*100)}', 'span_ladder',
                  (lambda c, t: (lambda r: mutate(r, aluY[:c], 6) + background(r, t - c)))(core, total), 8))

    # -- chr2: the HML-2 zone ------------------------------------------------
    S += [
        ('chr2', 'H01_solo_ltr', 'hervk_solo',
         lambda r: mutate(r, ltr5, 2), 8),
        ('chr2', 'H02_provirus', 'hervk_provirus',
         lambda r: mutate(r, ltr5, 2) + mutate(r, hervk, 3) + mutate(r, ltr5, 2), 8),
        ('chr2', 'H03_provirus_two_unit', 'hervk_cnv',
         lambda r: (mutate(r, ltr5, 2) + mutate(r, hervk, 3)) * 2 + mutate(r, ltr5, 2), 8),
        ('chr2', 'H04_truncated_int', 'hervk_partial',
         lambda r: mutate(r, ltr5, 2) + mutate(r, hervk[:3000], 4), 8),
        ('chr2', 'H05_hervk9_decoy', 'non_hml2_decoy',
         lambda r: mutate(r, cons['HERVK9'], 6), 8),
        ('chr2', 'H06_ltr5_sva_pair', 'hervk_sva_pair',
         lambda r: mutate(r, ltr5, 2) + mutate(r, hervk, 3) + mutate(r, sva_e[:400], 6) + mutate(r, ltr5, 2), 8),
    ]
    # -- the ploidy probes ---------------------------------------------------
    # Same element, same offset, different divergence so the alignments stay
    # separable against min_support = '2,4'.
    for contig, div in [('X', 4), ('chrX', 7), ('chrX_alt', 9), ('chrY', 4)]:
        S.append((contig, f'P_{contig}_aluY', 'ploidy_probe',
                  (lambda d: (lambda r: mutate(r, aluY, d) + 'A' * 30))(div), 8))
    return S


# ---------------------------------------------------------------- assembly

CONTIGS = [
    # name, background length. ASCII-lexicographic order is what vg autoindex
    # wants, and it is also the order bcftools sorts into.
    ('X', 4000), ('chr1', 60000), ('chr2', 60000), ('chr3_quiet', 20000),
    ('chr4_narrow', 20000), ('chrX', 4000), ('chrX_alt', 4000), ('chrY', 4000),
]
PITCH = 2500      # spacing between planted sites on chr1/chr2
MARGIN = 900      # keep every site clear of a contig end; tsd_win is 30 but
                  # the TSD flank extraction clamps near the start


def build(args):
    out = pathlib.Path(args.out)
    for d in ['ref', 'lib', 'hap', 'reads', 'vcf', 'bam']:
        (out / d).mkdir(parents=True, exist_ok=True)
    repo = pathlib.Path(args.repo).resolve()
    const = hervk_constants(repo)
    dfam = read_fasta(args.dfam)

    cons = {k: pick(dfam, k) for k in
            ['AluY', 'AluSx', 'SVA_E', 'LTR5_Hs', 'HERVK', 'HERVK9']}
    # One full-length L1HS out of the two halves the library ships.
    cons['L1HS'] = pick(dfam, 'L1HS_5end') + pick(dfam, 'L1HS_3end')

    if len(cons['HERVK']) != const[0]:
        sys.exit(f"build_synthetic.py: HERVK consensus is {len(cons['HERVK'])} bp, "
                 f"bin/hervk_arch.py expects {const[0]}")
    if len(cons['LTR5_Hs']) != const[1]:
        sys.exit(f"build_synthetic.py: LTR5_Hs is {len(cons['LTR5_Hs'])} bp, "
                 f"bin/hervk_arch.py expects {const[1]}")

    rng = random.Random(args.seed)
    scale = max(1, args.scale)
    # Sites keep their pitch; scaling adds unique background around them. That
    # is the knob for k-mer uniqueness, which is what decides whether PanGenie
    # can genotype at all.
    ref = {name: background(rng, n * scale) for name, n in CONTIGS}

    # Reference-resident copies. These are what hervk_ref_state.py masks and
    # what makes a deletion polymorphism possible at all.
    ref_copies = []
    for i, (elem, div) in enumerate([('AluY', 6), ('AluSx', 10), ('L1HS', 8), ('SVA_E', 6)]):
        pos = (3000 + i * 4000) * scale
        seq = mutate(rng, cons[elem][:600], div)
        ref['chr3_quiet'] = ref['chr3_quiet'][:pos] + seq + ref['chr3_quiet'][pos:]
        ref_copies.append(('chr3_quiet', pos, elem))
    # A reference solo LTR on chr2, so at least one HERV-K locus has a
    # non-empty reference state.
    solo_pos = 50000 * scale
    ref['chr2'] = (ref['chr2'][:solo_pos] + mutate(rng, cons['LTR5_Hs'], 3)
                   + ref['chr2'][solo_pos:])
    ref_copies.append(('chr2', solo_pos, 'LTR5_Hs'))

    # Plant a run of N in chr1 for --break_scaffolds.
    npos = 58000 * scale
    ref['chr1'] = ref['chr1'][:npos] + 'N' * 120 + ref['chr1'][npos + 120:]

    sites = plan_sites(cons, const)
    # Assign a position per site, per contig, on a fixed pitch.
    used = {}
    placed = []
    for contig, name, kind, builder, tsd_len in sites:
        idx = used.get(contig, 0)
        used[contig] = idx + 1
        pos = MARGIN + idx * PITCH * scale
        if pos + PITCH > len(ref[contig]):
            sys.exit(f'build_synthetic.py: {contig} too short for site {name}')
        placed.append((contig, pos, name, kind, builder, tsd_len))

    # Haplotypes. h1 and h2 carry the odd-indexed sites, h3 and h4 the even
    # ones, so every site is polymorphic across the four and none is fixed.
    haps = {h: {c: s for c, s in ref.items()} for h in ['h1', 'h2', 'h3', 'h4']}
    truth = []
    for i, (contig, pos, name, kind, builder, tsd_len) in enumerate(placed):
        carriers = ['h1', 'h2'] if i % 2 == 0 else ['h3', 'h4']
        if kind == 'ploidy_probe':
            carriers = ['h1', 'h2', 'h3', 'h4']
        seq = builder(rng)
        dup = tsd(rng, tsd_len) if tsd_len else ''
        insert = dup + seq + dup if dup else seq
        for h in carriers:
            s = haps[h][contig]
            haps[h][contig] = s[:pos] + insert + s[pos:]
        truth.append({'contig': contig, 'pos': pos, 'name': name, 'kind': kind,
                      'svlen': len(insert), 'tsd': dup or 'NONE',
                      'carriers': '|'.join(carriers), 'seq': insert})

    # ---- write ----------------------------------------------------------
    def write_fa(path, d):
        with open(path, 'w') as fh:
            for k, _ in CONTIGS:
                fh.write(f'>{k}\n')
                s = d[k]
                for j in range(0, len(s), 60):
                    fh.write(s[j:j + 60] + '\n')

    write_fa(out / 'ref' / 'synth.fa', ref)
    for h, d in sorted(haps.items()):
        write_fa(out / 'hap' / f'{h}.fa', d)

    # The library. Names carry the RepeatMasker "#Class/Family" convention,
    # which is what annotate_vcf.R reads matching_classes out of.
    lib_entries = [
        ('AluY', 'SINE/Alu', cons['AluY']),
        ('AluSx', 'SINE/Alu', cons['AluSx']),
        ('L1HS', 'LINE/L1', cons['L1HS']),
        ('SVA_E', 'Retroposon/SVA', cons['SVA_E']),
        ('LTR5_Hs', 'LTR/ERVK', cons['LTR5_Hs']),
        ('HERVK9-int', 'LTR/ERVK', cons['HERVK9']),
    ]
    for tag, intname in [('synth_TE', 'HERVK-int'), ('synth_TE_bare_HERVK', 'HERVK')]:
        with open(out / 'lib' / f'{tag}.fasta', 'w') as fh:
            for n, c, s in lib_entries + [(intname, 'LTR/ERVK', cons['HERVK'])]:
                fh.write(f'>{n}#{c}\n')
                for j in range(0, len(s), 60):
                    fh.write(s[j:j + 60] + '\n')

    with open(out / 'truth.tsv', 'w') as fh:
        cols = ['contig', 'pos', 'name', 'kind', 'svlen', 'tsd', 'carriers']  # 'seq' stays out: it is the insert, not a fact to assert
        fh.write('\t'.join(cols) + '\n')
        for r in truth:
            fh.write('\t'.join(str(r[c]) for c in cols) + '\n')

    # Samplesheets. Paths are absolute: nextflow.config runs --contain.
    ab = out.resolve()
    (out / 'assemblies.csv').write_text(
        'sample,path\n' + ''.join(f'{h},{ab}/hap/{h}.fa\n' for h in sorted(haps)))

    write_reads(rng, ref, haps, out, args)
    write_vcf_fixtures(out, ref, placed, truth)
    make_bams(out, args)
    manifest(out)
    print(f'built {len(placed)} sites across {len(CONTIGS)} contigs in {out}')
    return 0


# ---------------------------------------------------------------- reads

def write_reads(rng, ref, haps, out, args):
    """Short and long reads straight off the haplotypes. No simulator is in the
    container, and error-free reads are the right starting point: a genotyping
    failure should be the graph's, not the read set's. Errors are a knob to
    turn once the matrix is green."""
    def fq(path, recs):
        # mtime=0 and no embedded filename: gzip's default header carries both
        # and would make two identical runs differ.
        with open(path, 'wb') as raw:
            with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as gz:
                for name, seq in recs:
                    gz.write(f'@{name}\n{seq}\n+\n{"I" * len(seq)}\n'.encode())

    def sample_reads(d, rlen, depth, tag, paired):
        recs = []
        for c, _ in CONTIGS:
            s = d[c]
            n = max(1, int(len(s) * depth / rlen))
            for i in range(n):
                p = rng.randrange(0, max(1, len(s) - rlen))
                r = s[p:p + rlen]
                if len(r) < rlen:
                    continue
                if paired:
                    recs.append((f'{tag}_{c}_{i}/1', r))
                    m = max(0, p + 300 - rlen)
                    recs.append((f'{tag}_{c}_{i}/2', rc(d[c][m:m + rlen])))
                else:
                    recs.append((f'{tag}_{c}_{i}', r if i % 2 else rc(r)))
        rng.shuffle(recs)
        return recs

    short_depth = args.short_depth
    long_depth = args.long_depth
    # Two short-read samples for pangenie/giraffe, two long-read for
    # graphaligner and the sniffles entry points.
    fq(out / 'reads' / 'S1.short.fq.gz', sample_reads(haps['h1'], 150, short_depth, 'S1', True))
    fq(out / 'reads' / 'S2.short.fq.gz', sample_reads(haps['h3'], 150, short_depth, 'S2', True))
    fq(out / 'reads' / 'L1.long.fq.gz', sample_reads(haps['h2'], 5000, long_depth, 'L1', False))
    fq(out / 'reads' / 'L2.long.fq.gz', sample_reads(haps['h4'], 5000, long_depth, 'L2', False))

    ab = out.resolve()
    (out / 'reads.csv').write_text(
        'path,sample,type\n'
        f'{ab}/reads/S1.short.fq.gz,S1,short\n'
        f'{ab}/reads/S2.short.fq.gz,S2,short\n')
    (out / 'reads_long.csv').write_text(
        'path,sample,type\n'
        f'{ab}/reads/L1.long.fq.gz,L1,hifi\n'
        f'{ab}/reads/L2.long.fq.gz,L2,ont\n')
    (out / 'longreads.csv').write_text(
        'sample,path,type\n'
        f'L1,{ab}/reads/L1.long.fq.gz,hifi\n'
        f'L2,{ab}/reads/L2.long.fq.gz,ont\n')


def make_bams(out, args):
    """--bams needs coordinate-sorted BAMs. sniffles takes the sample name from
    @RG SM when it is there, so set it: the sample column of the output depends
    on it."""
    if not (shutil.which('minimap2') and shutil.which('samtools')):
        print('  [skip] minimap2/samtools not on PATH: no BAMs, --bams entry unavailable')
        return
    ab = out.resolve()
    rows = []
    for sample, fq in [('B1', 'L1.long.fq.gz'), ('B2', 'L2.long.fq.gz')]:
        bam = out / 'bam' / f'{sample}.bam'
        sam = subprocess.run(
            ['minimap2', '-ax', 'map-hifi', '-t', '2',
             '-R', f'@RG\\tID:{sample}\\tSM:{sample}',
             str(out / 'ref' / 'synth.fa'), str(out / 'reads' / fq)],
            capture_output=True, text=True)
        if sam.returncode != 0:
            print('  [warn] minimap2 failed:', sam.stderr.strip().splitlines()[-1:])
            return
        p = subprocess.run(['samtools', 'sort', '-o', str(bam)],
                           input=sam.stdout, text=True, capture_output=True)
        if p.returncode != 0:
            print('  [warn] samtools sort failed:', p.stderr.strip()[-200:])
            return
        subprocess.run(['samtools', 'index', str(bam)], check=True)
        rows.append(f'{sample},{ab}/bam/{sample}.bam\n')
    (out / 'bams.csv').write_text('sample,path\n' + ''.join(rows))


# ---------------------------------------------------------------- vcf inputs

VCF_HDR = """##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Variant length">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##FILTER=<ID=PASS,Description="passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
"""


def write_vcf_fixtures(out, ref, placed, truth):
    """Two hand-written entry points.

    svs_*.vcf.gz feed --svs, which applies no filter and no process, so these
    reach the truvari collapse verbatim. merged.vcf.gz feeds --vcf, the branch
    that keeps INFO and does NOT recompute SVLEN, so its SVLEN must be right in
    the file.

    chr4_narrow carries one symbolic <INV> and no indel. bcftools index -s
    lists only contigs holding records, so chr4_narrow DOES become a chunk,
    `bcftools view --types indels` returns nothing for it, and the no-indel
    exit in bin/repmask_vcf.sh fires. That is the only end-to-end route to
    that branch. Every INFO tag it uses is declared above; bcftools exits 255
    on an undeclared one.
    """
    contigs = ''.join(f'##contig=<ID={c},length={len(ref[c])}>\n' for c, _ in CONTIGS)

    # --svs: one VCF per haplotype, bgzipped, since the multi-VCF path tabixes
    # every staged input and then globs *.vcf.gz.
    for hap in ['h1', 'h2', 'h3', 'h4']:
        lines = [VCF_HDR + contigs + f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{hap}']
        for (contig, pos, name, kind, _b, _t), tr in zip(placed, truth):
            if hap not in tr['carriers'].split('|'):
                continue
            refb = ref[contig][pos - 1]
            alt = refb + tr['seq']
            lines.append(f'{contig}\t{pos}\t{name}\t{refb}\t{alt}\t.\tPASS\t'
                         f'SVTYPE=INS;SVLEN={len(tr["seq"])}\tGT\t1')
        p = out / 'vcf' / f'svs_{hap}.vcf'
        p.write_text('\n'.join(lines) + '\n')
        subprocess.run(['bgzip', '-f', str(p)], check=False)
        subprocess.run(['tabix', '-f', '-p', 'vcf', str(p) + '.gz'], check=False)

    # --vcf: one merged VCF, plus the chr4_narrow <INV> that has no indel.
    lines = [VCF_HDR + contigs + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tm1']
    for (contig, pos, name, kind, _b, _t), tr in zip(placed, truth):
        refb = ref[contig][pos - 1]
        alt = refb + tr['seq']
        lines.append(f'{contig}\t{pos}\t{name}\t{refb}\t{alt}\t.\tPASS\t'
                     f'SVTYPE=INS;SVLEN={len(tr["seq"])}\tGT\t1')
    lines.append('chr4_narrow\t5000\tchr4_inv_only\tA\t<INV>\t.\tPASS\t'
                 'SVTYPE=INV;END=5400;SVLEN=400\tGT\t1')
    p = out / 'vcf' / 'merged.vcf'
    p.write_text('\n'.join(lines) + '\n')
    subprocess.run(['bgzip', '-f', str(p)], check=False)
    subprocess.run(['tabix', '-f', '-p', 'vcf', str(p) + '.gz'], check=False)

    ab = out.resolve()
    (out / 'svs.csv').write_text(
        'sample,path\n' + ''.join(f'{h},{ab}/vcf/svs_{h}.vcf.gz\n' for h in ['h1', 'h2', 'h3', 'h4']))
    # One-VCF sheet: truvari_merge's num_files==1 branch with from_vcf=false,
    # which is the single-input path that DOES tabix its input. The --vcf entry
    # takes the other single-input branch, so both need covering.
    (out / 'svs_one.csv').write_text(f'sample,path\nh1,{ab}/vcf/svs_h1.vcf.gz\n')

    # --epigenomes with --lifted: the CSV of already-lifted methylation calls.
    # --lifted short-circuits bamtags_to_BED and lift_epigenome, so this needs
    # no MM/ML BAM and no methylation simulator. Columns follow the graph node
    # schema annotate_vcf.py reads; two rows per sample is enough to show the
    # fields reach the merged VCF.
    (out / 'lifted.csv').write_text(
        'sample,path\n' + ''.join(f'{s},{ab}/meth/{s}.methylation.csv\n' for s in ['S1', 'S2']))
    (out / 'epigenomes.csv').write_text(
        'sample,path\n' + ''.join(f'{s},{ab}/bam/B1.bam\n' for s in ['S1', 'S2']))
    meth = out / 'meth'
    meth.mkdir(exist_ok=True)
    for smp in ['S1', 'S2']:
        (meth / f'{smp}.methylation.csv').write_text(
            'node,levels,coverage\n1,0.8,10\n2,0.2,12\n')


# ---------------------------------------------------------------- manifest

# The samplesheets and the BAMs embed the absolute output path, so they differ
# between two builds in different directories even when the data is identical.
# The manifest covers the data; PATH_DEPENDENT is listed separately so that a
# mismatch there is read as "different --out", not as a broken generator.
PATH_DEPENDENT = ('.csv', '.bam', '.bai')


def manifest(out):
    out = pathlib.Path(out)
    lines, skipped = [], []
    for p in sorted(out.rglob('*')):
        if not p.is_file() or p.name == 'MANIFEST.sha256':
            continue
        rel = p.relative_to(out)
        if p.suffix in PATH_DEPENDENT:
            skipped.append(str(rel))
            continue
        lines.append(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {rel}')
    body = '\n'.join(lines) + '\n'
    body += ''.join(f'# path-dependent, not hashed: {s}\n' for s in skipped)
    (out / 'MANIFEST.sha256').write_text(body)
    print(f'  manifest: {len(lines)} hashed, {len(skipped)} path-dependent')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dfam', required=True, help='human_DFAM3.6.fasta from test/human_test_set.tar.gz')
    ap.add_argument('--repo', default='.', help='GraffiTE checkout, for bin/hervk_arch.py')
    ap.add_argument('--out', default='build')
    ap.add_argument('--seed', type=int, default=20260918)
    ap.add_argument('--short-depth', type=int, default=30)
    ap.add_argument('--long-depth', type=int, default=15)
    ap.add_argument('--scale', type=int, default=1,
                    help='multiply every contig length. Raise it when a k-mer '
                         'or alignment-based stage has too little unique '
                         'sequence to work with; the planted sites do not move.')
    sys.exit(build(ap.parse_args()))


if __name__ == '__main__':
    main()
