#!/usr/bin/env python3
"""
HERV-K (HML-2) locus layer.

`flag` (stage 3, this phase)
    Group HERV-K candidate records into loci and mark them in place. Nothing is
    merged or dropped here: the discovery VCF is what induces the graph, and
    the human subset is what the analysis reads, so both must keep their record
    structure exactly. The locus grouping is written out as INFO tags plus a
    dedicated table for downstream use.

    Grouping is anchored on the reference HML-2 element when one is known, and
    otherwise on proximity within one LTR length. That radius is derivable
    rather than arbitrary: two records describing the same insertion can be
    placed anywhere inside the shared reference LTR, so they land up to one LTR
    apart. In the CaG cohort the observed offsets are 574 bp (chr11), 559 bp
    (chr6), 107 bp (chr12) and 871 bp (chr1) -- all within one LTR, and the
    first three sit outside truvari's default 500 bp refdist, which is why
    truvari did not collapse them.

`consolidate` (stage 4)
    Collapse each flagged locus in the graph-genotyped VCF into one
    multi-allelic record.

    Genotypes are resolved by *allele dosage* across the member records, not
    by taking each record's call at face value. A locus is one place with one
    allele set; the members are partial views of it, and `vg call` frequently
    leaves one member uncalled because the reads took another member's path.
    Measured on the CaG set: that structural missingness (no DP at all, not
    low coverage) affects 9 samples at chr6 and 11 at chr12.

    Summing dosages fixes most of it. Where the called members already account
    for every haplotype, an uncalled member is pinned to zero -- arithmetic,
    not inference. At chr12 that resolves 17 of 20 samples where treating any
    missing member as fatal resolves 7.

    It also gives a free consistency check: dosages cannot exceed ploidy.
    Three CaG samples violate it (chr12 x2, chr8 x1), and at chr12 those two
    samples are exactly the gap between the graph and the assemblies. They are
    the third allele being flattened by `bcftools norm -m-`.

    No genotype is ever invented. Where dosage leaves haplotypes unaccounted
    the record is missing there, and AC and AN are reported separately so the
    shortfall is visible: at chr6 the graph finds every solo carrier (AC=8,
    matching the assemblies exactly) and only fails to confirm non-carriers
    (AN=22 of 40). Reporting AF alone would hide that; reporting AC/AN does
    not.

Usage:
    hervk_reconcile.py flag --vcf-in in.vcf --vcf-out out.vcf \
        --loci-out hervk_loci.tsv [--ref-state refstate.tsv] [--window 1200]
"""

import argparse
import gzip
import re
import sys
from collections import defaultdict

AUTOSOME_PREFIXES = tuple(f'chr{i}' for i in range(1, 23))

INFO_HEADERS = [
    '##INFO=<ID=HERVK_LOCUS,Number=1,Type=String,Description="HERV-K locus '
    'identifier grouping records that describe the same element.">',
    '##INFO=<ID=HERVK_LOCUS_N,Number=1,Type=Integer,Description="Number of '
    'HERV-K records assigned to this locus.">',
    '##INFO=<ID=HERVK_MERGE_FLAG,Number=0,Type=Flag,Description="This locus '
    'holds more than one record describing the same element. Flagged only -- '
    'records are never merged at this stage, because this VCF must keep its '
    'structure for graph induction and downstream ID matching.">',
    '##INFO=<ID=HERVK_POLARITY_CONFLICT,Number=0,Type=Flag,Description="Records '
    'at this locus imply different REF allele states. The masked reference is '
    'authoritative; consolidation resolves this.">',
]

LOCI_COLUMNS = ['locus_id', 'chrom', 'start', 'end', 'n_records', 'record_ids',
                'n_in_human', 'records_not_in_human', 'ref_state', 'allele_set',
                'per_record_class', 'per_record_evidence', 'per_record_k',
                'arch', 'flags']


def parse_info(info):
    d = {}
    if not info or info == '.':
        return d
    for kv in info.split(';'):
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k] = v
        else:
            d[kv] = ''
    return d


def info_to_str(d):
    return ';'.join(k if v == '' else f'{k}={v}' for k, v in d.items()) or '.'


def load_table(path, key='id'):
    if not path:
        return {}
    rows = {}
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        ki = header.index(key)
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < len(header):
                f += [''] * (len(header) - len(f))
            rows[f[ki]] = dict(zip(header, f))
    return rows


def read_calls(path):
    """Load every HERV-K candidate from the classifier's call table.

    Clustering must see all candidates, not only the ones in the human subset.
    The --human filter requires FILTER="PASS", and PAV emits TRIM/COMPOUND on
    perfectly real HML-2 records (chr15-2092086-DEL-8221 is one), so grouping
    from the human VCF alone could split a locus whose partner was removed for
    an unrelated reason.
    """
    recs = []
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            if not line.strip():
                continue
            r = dict(zip(header, line.rstrip('\n').split('\t')))
            # Non-HML-2 LTR/ERVK (HERVK9-int, MER11A, LTR13 ...) are reported
            # as `other` by the classifier but must not form loci: they are a
            # different lineage and grouping them would invent merge flags.
            if r['evidence'] == 'NON_HML2':
                continue
            pos, svlen = int(r['pos']), int(r['svlen'])
            recs.append({
                'id': r['id'], 'chrom': r['chrom'], 'pos': pos, 'svlen': svlen,
                'start': pos, 'end': pos if svlen >= 0 else pos + abs(svlen),
                'cls': r['class'], 'evidence': r['evidence'],
                'k': '' if r['k'] == '.' else r['k'],
                'ref_allele': r['allele_ref'], 'alt_allele': r['allele'],
                'ref_state': r['ref_state'], 'arch': r['arch'],
            })
    return recs


def vcf_ids(vcf_path):
    ids = set()
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t', 3)
            if len(f) > 2:
                ids.add(f[2])
    return ids


def cluster(recs, ref_tbl, window, human_ids=None):
    """Assign each record to a locus.

    Two records share a locus when they touch the same reference HML-2 element,
    or when their footprints lie within `window` bp of each other.
    """
    by_chrom = defaultdict(list)
    for r in recs:
        by_chrom[r['chrom']].append(r)

    loci = []
    for chrom in sorted(by_chrom):
        items = sorted(by_chrom[chrom], key=lambda r: (r['start'], r['end']))
        current = []
        cur_end = None
        cur_elems = set()
        for r in items:
            ref = ref_tbl.get(r['id'], {})
            elem = None
            if ref.get('ref_elem_start') and ref.get('ref_elem_end'):
                elem = (ref['ref_elem_chrom'], ref['ref_elem_start'],
                        ref['ref_elem_end'])
            joins = bool(current) and (
                (elem is not None and elem in cur_elems)
                or r['start'] - cur_end <= window)
            if not joins and current:
                loci.append(current)
                current, cur_elems, cur_end = [], set(), None
            current.append(r)
            if elem is not None:
                cur_elems.add(elem)
            cur_end = max(cur_end or r['end'], r['end'])
        if current:
            loci.append(current)

    assignment = {}
    table = []
    for members in loci:
        chrom = members[0]['chrom']
        start = min(m['start'] for m in members)
        end = max(m['end'] for m in members)
        locus_id = f'HERVK_{chrom}_{start}'

        ref_states = {m['ref_allele'] for m in members if m['ref_allele'] not in ('.', '')}
        alleles = set()
        for m in members:
            for a in (m['ref_allele'], m['alt_allele']):
                if a not in ('.', '', 'None'):
                    alleles.add(a)

        in_human = [m['id'] for m in members
                    if human_ids is None or m['id'] in human_ids]
        absent = [m['id'] for m in members if m['id'] not in in_human]

        flags = []
        if len(members) > 1:
            flags.append('MERGE_CANDIDATE')
        # Only interesting when the locus is *split* by the human filter --
        # some members annotated, some not. A locus wholly outside the human
        # subset is simply out of scope, not a problem.
        if absent and in_human:
            flags.append('LOCUS_SPLIT_BY_HUMAN_FILTER')
        if len(ref_states) > 1:
            flags.append('POLARITY_CONFLICT')
        if len(alleles) > 2:
            flags.append('MULTIALLELIC')
        if not chrom.startswith(AUTOSOME_PREFIXES):
            flags.append('PLOIDY_UNVERIFIED')

        for m in members:
            assignment[m['id']] = (locus_id, len(members), flags)

        table.append({
            'locus_id': locus_id, 'chrom': chrom,
            'start': str(start), 'end': str(end),
            'n_records': str(len(members)),
            'record_ids': ','.join(m['id'] for m in members),
            'n_in_human': str(len(in_human)),
            'records_not_in_human': ','.join(absent) or '.',
            'ref_state': ','.join(sorted(ref_states)) or '.',
            'allele_set': ','.join(sorted(alleles)) or '.',
            'per_record_class': ','.join(m['cls'] for m in members),
            'per_record_evidence': ','.join(m['evidence'] for m in members),
            'per_record_k': ','.join(m['k'] or '.' for m in members),
            'arch': ';'.join(m['arch'] or '.' for m in members),
            'flags': ','.join(flags) or '.',
        })
    return assignment, table


def write_vcf(vcf_in, vcf_out, assignment):
    with open(vcf_in) as fin, open(vcf_out, 'w') as fout:
        header_done = False
        for line in fin:
            if line.startswith('##'):
                fout.write(line)
                continue
            if line.startswith('#CHROM'):
                for h in INFO_HEADERS:
                    fout.write(h + '\n')
                fout.write(line)
                header_done = True
                continue
            if not header_done or not line.strip():
                fout.write(line)
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 8 or f[2] not in assignment:
                fout.write(line)
                continue
            locus_id, n, flags = assignment[f[2]]
            info = parse_info(f[7])
            info['HERVK_LOCUS'] = locus_id
            info['HERVK_LOCUS_N'] = str(n)
            if 'MERGE_CANDIDATE' in flags:
                info['HERVK_MERGE_FLAG'] = ''
            if 'POLARITY_CONFLICT' in flags:
                info['HERVK_POLARITY_CONFLICT'] = ''
            f[7] = info_to_str(info)
            fout.write('\t'.join(f) + '\n')


def write_loci(table, path):
    with open(path, 'w') as fh:
        fh.write('\t'.join(LOCI_COLUMNS) + '\n')
        for row in table:
            fh.write('\t'.join(row[c] for c in LOCI_COLUMNS) + '\n')


def cmd_flag(args):
    recs = read_calls(args.calls)
    ref_tbl = load_table(args.ref_state)
    human = vcf_ids(args.vcf_in)
    assignment, table = cluster(recs, ref_tbl, args.window, human)
    write_vcf(args.vcf_in, args.vcf_out, assignment)
    write_loci(table, args.loci_out)

    multi = [r for r in table if int(r['n_records']) > 1]
    conflict = [r for r in table if 'POLARITY_CONFLICT' in r['flags']]
    partial = [r for r in table if 'LOCUS_SPLIT_BY_HUMAN_FILTER' in r['flags']]
    sys.stderr.write(
        f'[hervk_reconcile] {len(recs)} candidates -> {len(table)} loci '
        f'({len(human & {r["id"] for r in recs})} annotated in the human VCF); '
        f'{len(multi)} flagged for merge, {len(conflict)} polarity conflicts, '
        f'{len(partial)} loci split by the human filter\n')
    for r in multi:
        sys.stderr.write(f'    {r["locus_id"]}: {r["record_ids"]} [{r["flags"]}]\n')


SUPPORTED_GENOTYPERS = ('giraffe',)


def detect_genotyper(vcf_path):
    """Identify which back end wrote a genotyped VCF, from its header alone.

    Needed because --hervk_reconcile_vcf points at a VCF from an *earlier* run,
    whose back end is a property of that file and not of this run's
    params.graph_method (which is moot anyway -- that path pairs with
    --genotype false, so nothing is being genotyped here).

    vg call declares FORMAT/MAD ("Minimum site allele depth") and FORMAT/XD
    ("eXpected Depth ... Poisson model"); PanGenie declares neither and names
    itself in ##source or ##commandline.

    giraffe, graphaligner and precomputed all genotype *through* vg call, so a
    vg-shaped file is reported as 'giraffe'. That names the validated VCF
    shape, not a claim about which aligner produced the GAM -- and since
    --genotyper only selects the guard, the distinction has no effect on the
    consolidation itself.

    Returns a genotyper name, or None when the header settles nothing.
    """
    fmt_ids, other = set(), []
    opener = gzip.open if vcf_path.endswith('.gz') else open
    with opener(vcf_path, 'rt') as fh:
        for line in fh:
            if not line.startswith('##'):
                break  # #CHROM or the first record: header is done
            if line.startswith('##FORMAT=<ID='):
                fmt_ids.add(line.split('##FORMAT=<ID=', 1)[1].split(',', 1)[0])
            else:
                other.append(line.lower())

    if 'pangenie' in ''.join(other):
        return 'pangenie'
    if {'MAD', 'XD'} <= fmt_ids:
        return 'giraffe'
    if 'KC' in fmt_ids:          # PanGenie's kmer count, if it named nothing
        return 'pangenie'
    return None


ID_RE = re.compile(r'##(?:INFO|FORMAT)=<ID=([^,>]+)')


def merge_headers(head, new_lines):
    """Drop any INFO/FORMAT definition from `head` that `new_lines` redefines.

    Two ##INFO lines with the same ID is invalid VCF and readers disagree about
    which one wins, so ours must replace rather than append. Number=A is the
    correct shape for both kinds of record here: a consolidated locus carries
    one SVTYPE/SVLEN per ALT, and a record carried over unchanged is biallelic,
    so one ALT means one value. Filtering by ID also makes a second
    consolidation pass over our own output idempotent.
    """
    ids = set()
    for h in new_lines:
        m = ID_RE.match(h)
        if m:
            ids.add(m.group(1))
    keep = []
    for h in head:
        m = ID_RE.match(h)
        if m and m.group(1) in ids:
            continue
        keep.append(h)
    return keep


CONSOLIDATED_HEADERS = [
    '##INFO=<ID=SVTYPE,Number=A,Type=String,Description="Variant type per ALT.">',
    '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Variant length per ALT.">',
    '##INFO=<ID=HERVK_LOCUS,Number=1,Type=String,Description="HERV-K locus id.">',
    '##INFO=<ID=HERVK_ALLELE_REF,Number=1,Type=String,Description="HERV-K state '
    'of the REF allele.">',
    '##INFO=<ID=HERVK_ALLELE,Number=.,Type=String,Description="HERV-K state of '
    'each ALT allele, in ALT order.">',
    '##INFO=<ID=HERVK_MEMBERS,Number=.,Type=String,Description="Record IDs '
    'consolidated into this locus.">',
    '##INFO=<ID=HERVK_AC,Number=.,Type=Integer,Description="Allele count per '
    'ALT, after dosage resolution across the member records.">',
    '##INFO=<ID=HERVK_AN,Number=1,Type=Integer,Description="Total alleles '
    'called at this locus. Compare against 2N: a shortfall means members were '
    'structurally uncalled, not that the alleles are absent.">',
    '##INFO=<ID=HERVK_N_RESOLVED,Number=1,Type=Integer,Description="Samples '
    'fully resolved by dosage.">',
    '##INFO=<ID=HERVK_N_PARTIAL,Number=1,Type=Integer,Description="Samples with '
    'some haplotypes unaccounted; those haplotypes are reported missing.">',
    '##INFO=<ID=HERVK_N_PLOIDY_EXCEEDED,Number=1,Type=Integer,Description='
    '"Samples whose summed ALT dosage exceeds ploidy -- a third allele lost to '
    'bcftools norm -m-. Their genotypes are set missing.">',
    '##INFO=<ID=HERVK_AC_DISC,Number=.,Type=Integer,Description="Allele count '
    'from the discovery (assembly) callset, per ALT.">',
    '##INFO=<ID=HERVK_AN_DISC,Number=1,Type=Integer,Description="Total alleles '
    'in the discovery callset.">',
    '##INFO=<ID=HERVK_DISC_CONCORDANT,Number=0,Type=Flag,Description="Graph and '
    'discovery agree on every ALT count at this locus.">',
    '##INFO=<ID=HERVK_GT_MASKED,Number=0,Type=Flag,Description="Genotypes '
    'withheld: a tandem duplication, which is neither transposition nor '
    'intra-element recombination and so does not belong in allele frequencies. '
    'The call itself is kept in the HERV-K tables.">',
    '##INFO=<ID=HERVK_MEMBERS_MASKED,Number=.,Type=String,Description="Locus '
    'members whose genotypes were withheld.">',
    '##INFO=<ID=HERVK_MEMBERS_UNRESOLVED,Number=.,Type=String,Description="Locus '
    'members with no usable allele state.">',
    '##INFO=<ID=HERVK_POLARITY_FLIPPED,Number=.,Type=String,Description="Members '
    'whose own polarity disagreed with the locus REF state; re-expressed '
    'against it.">',
]


def parse_gt(gt):
    """Return (alleles, ploidy). alleles is None when the call is missing."""
    field = gt.split(':', 1)[0]
    a = field.replace('|', '/').split('/')
    ploidy = len(a)
    if any(x == '.' for x in a):
        return None, ploidy
    return [int(x) for x in a], ploidy


def structurally_missing(sample_field, fmt):
    """Missing with no depth reported at all -- the record was never evaluated,
    as opposed to evaluated and found ambiguous. Every missing member call at
    chr6 and chr12 is of this kind."""
    d = dict(zip(fmt, sample_field.split(':')))
    return d.get('DP', '.') in ('.', '', None)


def resolve_sample(member_gts, fmt_list, ploidy_hint):
    """Resolve one sample at one locus by ALT dosage across member records.

    Returns (dosages, status) where dosages maps member index -> ALT dosage,
    and status is one of resolved | partial | ploidy_exceeded.
    """
    dos, ploidy, missing = {}, ploidy_hint, []
    for i, (gt, fmt) in enumerate(zip(member_gts, fmt_list)):
        alleles, p = parse_gt(gt)
        ploidy = max(ploidy, p)
        if alleles is None:
            missing.append(i)
        else:
            dos[i] = sum(1 for x in alleles if x > 0)

    known = sum(dos.values())
    if known > ploidy:
        return dos, 'ploidy_exceeded'
    if not missing:
        return dos, 'resolved'
    if known >= ploidy:
        # Every haplotype is already spoken for, so the uncalled members carry
        # nothing. This is arithmetic, not an assumption.
        for i in missing:
            dos[i] = 0
        return dos, 'resolved'
    return dos, 'partial'


def build_gt(dos, allele_index, ploidy, status):
    """Compose the consolidated GT from per-member ALT dosages."""
    if status in ('ploidy_exceeded',):
        return '/'.join(['.'] * ploidy)
    call = []
    for i, d in sorted(dos.items()):
        call.extend([str(allele_index[i])] * d)
    if status == 'partial':
        # Report what is known and leave the rest missing rather than filling
        # it with reference.
        call.extend(['.'] * (ploidy - len(call)))
        return '/'.join(call[:ploidy])
    call.extend(['0'] * (ploidy - len(call)))
    return '/'.join(sorted(call[:ploidy]))


def read_vcf_records(path, wanted):
    """Return (header_lines, samples, {id: fields}) for the wanted IDs."""
    opener = gzip.open if path.endswith('.gz') else open
    head, samples, recs = [], [], {}
    with opener(path, 'rt') as fh:
        for line in fh:
            if line.startswith('##'):
                head.append(line.rstrip('\n')); continue
            if line.startswith('#CHROM'):
                samples = line.rstrip('\n').split('\t')[9:]; continue
            f = line.rstrip('\n').split('\t')
            if f[2] in wanted:
                recs[f[2]] = f
    return head, samples, recs


def discovery_counts(path, members, samples):
    """Per-member ALT counts from the assembly-based discovery callset.

    An independent measurement of the same haplotypes. It is what caught the
    flattening at chr12, so it is reported beside the graph counts rather than
    used to correct them.
    """
    if not path:
        return None, None
    _, dsamples, recs = read_vcf_records(path, set(members))
    if not recs:
        return None, None
    idx = {s: i for i, s in enumerate(dsamples)}
    counts, an = {}, 0
    for m in members:
        counts[m] = 0
    for s in samples:
        if s not in idx:
            continue
        j = idx[s]
        seen = False
        for m in members:
            if m not in recs:
                continue
            alleles, _ = parse_gt(recs[m][9 + j])
            if alleles is None:
                continue
            seen = True
            counts[m] += sum(1 for x in alleles if x > 0)
        if seen:
            an += 2
    return counts, an


def cmd_consolidate(args):
    if args.genotyper == 'auto':
        detected = detect_genotyper(args.genotyped_vcf)
        if detected is None:
            sys.exit('hervk_reconcile consolidate: --genotyper auto could not '
                     f'identify the back end that wrote {args.genotyped_vcf} '
                     '(no vg call FORMAT/MAD+XD, no PanGenie marker). Pass '
                     '--genotyper explicitly.')
        sys.stderr.write(f'[hervk_reconcile] --genotyper auto: header says '
                         f'{detected}\n')
        args.genotyper = detected

    if args.genotyper not in SUPPORTED_GENOTYPERS:
        sys.exit(f'hervk_reconcile consolidate: --genotyper {args.genotyper} is '
                 f'not supported yet (only {", ".join(SUPPORTED_GENOTYPERS)}). '
                 'The internals are back-end agnostic; the guard is here so an '
                 'unvalidated back end fails loudly instead of quietly.')

    loci = [r for r in load_table(args.loci, 'locus_id').values()]
    flagged = [r for r in loci if int(r['n_records']) > 1]
    if not flagged:
        sys.stderr.write('[hervk_reconcile] no multi-record loci to consolidate\n')

    members_all = set()
    for r in flagged:
        members_all.update(r['record_ids'].split(','))

    calls = load_table(args.calls, 'id') if args.calls else {}
    # tandem records may sit outside any flagged locus, so they have to be in
    # the read set too
    members_all |= {vid for vid, c in calls.items()
                    if c.get('class') == 'tandem_prov'}
    head, samples, recs = read_vcf_records(args.genotyped_vcf, members_all)

    report, consolidated, archived, dropped = [], {}, [], set()
    annotate_only = {}

    # Every tandem-duplication record gets its genotypes withheld, not just the
    # ones that happen to sit in a multi-record locus. Two of the three in CaG
    # (chr7-4700334, chr12-133148145) are single-record loci and so are never
    # reached by the loop below -- they were being masked in the discovery VCF
    # and left callable in the genotyped one, which is exactly the
    # inconsistency the masking exists to prevent.
    mask_gt = {vid for vid, c in calls.items()
               if c.get('class') == 'tandem_prov'}
    for locus in flagged:
        mem = [m for m in locus['record_ids'].split(',') if m in recs]
        if len(mem) < 2:
            report.append((locus['locus_id'], 'skipped',
                           f'{len(mem)} of {len(locus["record_ids"].split(","))} '
                           'members present in the genotyped VCF'))
            continue

        # Locus REF state is authoritative and comes from the masked reference.
        ref_state = (locus['ref_state'].split(',')[0]
                     if locus['ref_state'] not in ('.', '') else '.')
        alt_states, allele_index, flipped, masked, unresolved = [], {}, [], [], []
        for i, m in enumerate(mem):
            c = calls.get(m, {})
            st = c.get('allele', '.')
            if st in ('.', '', 'NA') or c.get('class') in ('other', ''):
                # No allele state means no allele. Counting "." as a distinct
                # ALT would turn a biallelic locus into a spurious multiallelic
                # one and then demand a spanning deletion to build it.
                unresolved.append(m)
                continue
            if c.get('class') == 'tandem_prov':
                # Withheld here as well as in the discovery VCF: a tandem
                # duplication is not an ERV life-cycle event and must not enter
                # allele frequencies from either callset.
                masked.append(m)
            if c.get('allele_ref', '.') not in ('.', '') and \
               ref_state not in ('.', '') and c['allele_ref'] != ref_state:
                flipped.append(m)
            if st in alt_states:
                allele_index[i] = alt_states.index(st) + 1
            else:
                alt_states.append(st)
                allele_index[i] = len(alt_states)

        keep = [i for i, m in enumerate(mem)
                if m not in masked and m not in unresolved]
        if len(keep) < 2:
            # Nothing to merge, but the locus is still real and the surviving
            # member should carry its identity rather than look like a lone
            # unrelated record.
            why = ('one usable member (others masked or unresolved)'
                   if keep else 'no member carries a usable allele state')
            if keep:
                annotate_only[mem[keep[0]]] = (locus, masked, unresolved)
            report.append((locus['locus_id'], 'annotated' if keep else 'skipped', why))
            continue

        fmt_list = [recs[mem[i]][8].split(':') for i in keep]
        n_res = n_part = n_viol = 0
        gts, struct_miss = [], 0
        ac = {i: 0 for i in keep}
        an = 0
        for si in range(len(samples)):
            gt_fields = [recs[mem[i]][9 + si] for i in keep]
            for gf, fm in zip(gt_fields, fmt_list):
                a, _ = parse_gt(gf)
                if a is None and structurally_missing(gf, fm):
                    struct_miss += 1
            _, ploidy = parse_gt(gt_fields[0])
            dos, status = resolve_sample(gt_fields, fmt_list, ploidy)
            local = {keep[k]: v for k, v in
                     zip(range(len(keep)), [dos.get(k, 0) for k in range(len(keep))])}
            if status == 'resolved':
                n_res += 1
                for i, d in local.items():
                    ac[i] += d
                an += ploidy
            elif status == 'partial':
                n_part += 1
                for i, d in local.items():
                    ac[i] += d
                an += sum(local.values())
            else:
                n_viol += 1
            gts.append(build_gt({k: dos.get(k, 0) for k in range(len(keep))},
                                {k: allele_index[keep[k]] for k in range(len(keep))},
                                ploidy, status))

        alt_ac = {}
        for i in keep:
            alt_ac.setdefault(allele_index[i], 0)
            alt_ac[allele_index[i]] += ac[i]

        dcounts, dan = discovery_counts(args.discovery_vcf,
                                        [mem[i] for i in keep], samples)
        disc_ac = {}
        if dcounts:
            for i in keep:
                disc_ac.setdefault(allele_index[i], 0)
                disc_ac[allele_index[i]] += dcounts.get(mem[i], 0)

        consolidated[locus['locus_id']] = {
            'locus': locus, 'mem': mem, 'keep': keep, 'gts': gts,
            'alt_states': alt_states, 'allele_index': allele_index,
            'ref_state': ref_state, 'ac': alt_ac, 'an': an,
            'n_res': n_res, 'n_part': n_part, 'n_viol': n_viol,
            'struct_miss': struct_miss, 'flipped': flipped, 'masked': masked,
            'disc_ac': disc_ac, 'disc_an': dan,
        }
        dropped.update(mem)
        archived.extend(recs[m] for m in mem)
        report.append((locus['locus_id'], 'consolidated', ''))

    write_consolidated(args, head, samples, consolidated, dropped,
                       archived, report, recs, mask_gt, annotate_only)


def build_locus_record(mem_recs, keep, alt_states, allele_index):
    """Build one multi-allelic record with literal REF/ALT sequence.

    Three shapes occur, and only the first two are reachable without the
    reference FASTA:

      a) every member describes the same allele (chr1, chr8, chr11) -- emit the
         leftmost member's own REF/ALT unchanged; it already is that allele.
      b) members describe different alleles and one is a deletion spanning the
         locus (chr12) -- its REF field *is* the reference sequence over the
         span, so the other alleles can be spliced into it exactly.
      c) different alleles with no spanning deletion -- needs the FASTA.

    Nothing is reconstructed from coordinates alone; every base emitted here
    comes from a REF or ALT field that the caller already wrote.
    """
    kept = [mem_recs[i] for i in keep]
    anchor = min(kept, key=lambda r: int(r[1]))

    if len(alt_states) == 1:
        return anchor[0], anchor[1], anchor[3], anchor[4], None

    spanning = None
    for r in kept:
        span_end = int(r[1]) + len(r[3]) - 1
        if len(r[3]) > len(r[4]) and all(int(o[1]) >= int(r[1]) and
                                         int(o[1]) <= span_end for o in kept):
            spanning = r
            break
    if spanning is None:
        return None, None, None, None, (
            'members describe different alleles and no deletion spans the '
            'locus, so the reference sequence over the span is not available '
            'from the records; pass --reference to splice it from the FASTA')

    ref_seq = spanning[3]
    base = int(spanning[1])
    alts = [None] * len(alt_states)
    for i, r in zip(keep, kept):
        ai = allele_index[i] - 1
        if r is spanning:
            alts[ai] = r[4]
        else:
            off = int(r[1]) - base
            ins = r[4][1:] if len(r[4]) > len(r[3]) else ''
            skip = len(r[3]) - 1 if len(r[3]) > len(r[4]) else 0
            alts[ai] = ref_seq[:off + 1] + ins + ref_seq[off + 1 + skip:]
    if any(a is None for a in alts):
        return None, None, None, None, 'could not build every ALT allele'
    return spanning[0], spanning[1], ref_seq, ','.join(alts), None


def mask_genotypes(fields):
    """Blank every genotype, preserving ploidy."""
    for i in range(9, len(fields)):
        parts = fields[i].split(':')
        n = len(parts[0].replace('|', '/').split('/'))
        parts[0] = '/'.join(['.'] * n)
        fields[i] = ':'.join(parts)
    return fields


def write_consolidated(args, head, samples, consolidated, dropped, archived,
                       report, recs, mask_gt=frozenset(), annotate_only=None):
    out_by_id, skipped = {}, []
    for lid, c in consolidated.items():
        mem_recs = [recs[m] for m in c['mem']]
        chrom, pos, ref, alt, err = build_locus_record(
            mem_recs, c['keep'], c['alt_states'], c['allele_index'])
        if err:
            skipped.append((lid, err))
            report.append((lid, 'skipped', err))
            for m in c['mem']:
                dropped.discard(m)
            continue

        # SVTYPE/SVLEN per ALT. A consolidated record without them is badly
        # formed for anything downstream that keys off variant type, and the
        # values are not recoverable from the members once the alleles have
        # been re-expressed against a common REF.
        alt_list = alt.split(',')
        info = {
            'SVTYPE': ','.join('DEL' if len(a) < len(ref) else 'INS'
                               for a in alt_list),
            'SVLEN': ','.join(str(len(a) - len(ref)) for a in alt_list),
            'HERVK_LOCUS': lid,
            'HERVK_ALLELE_REF': c['ref_state'] or '.',
            'HERVK_ALLELE': ','.join(c['alt_states']),
            'HERVK_MEMBERS': ','.join(c['mem']),
            'HERVK_AC': ','.join(str(c['ac'].get(i + 1, 0))
                                 for i in range(len(c['alt_states']))),
            'HERVK_AN': str(c['an']),
            'HERVK_N_RESOLVED': str(c['n_res']),
            'HERVK_N_PARTIAL': str(c['n_part']),
            'HERVK_N_PLOIDY_EXCEEDED': str(c['n_viol']),
        }
        if c['flipped']:
            info['HERVK_POLARITY_FLIPPED'] = ','.join(c['flipped'])
        if c['masked']:
            info['HERVK_MEMBERS_MASKED'] = ','.join(c['masked'])
        if c['disc_an']:
            info['HERVK_AC_DISC'] = ','.join(str(c['disc_ac'].get(i + 1, 0))
                                             for i in range(len(c['alt_states'])))
            info['HERVK_AN_DISC'] = str(c['disc_an'])
            if all(c['ac'].get(k, 0) == c['disc_ac'].get(k, 0)
                   for k in set(c['ac']) | set(c['disc_ac'])):
                info['HERVK_DISC_CONCORDANT'] = ''
        fields = [chrom, pos, lid, ref, alt, '.', 'PASS', info_to_str(info), 'GT']
        fields.extend(c['gts'])
        out_by_id[c['mem'][0]] = fields

    with open(args.out_vcf, 'w') as fh:
        for h in merge_headers(head, CONSOLIDATED_HEADERS):
            fh.write(h + '\n')
        for h in CONSOLIDATED_HEADERS:
            fh.write(h + '\n')
        fh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
                 + '\t'.join(samples) + '\n')
        annotate_only = annotate_only or {}
        opener = gzip.open if args.genotyped_vcf.endswith('.gz') else open
        with opener(args.genotyped_vcf, 'rt') as gin:
            for line in gin:
                if line.startswith('#'):
                    continue
                vid = line.split('\t', 3)[2]
                if vid in out_by_id:
                    fh.write('\t'.join(out_by_id[vid]) + '\n')
                    continue
                if vid in dropped:
                    continue
                f = line.rstrip('\n').split('\t')
                if vid in mask_gt:
                    f = mask_genotypes(f)
                    info = parse_info(f[7])
                    info['HERVK_GT_MASKED'] = ''
                    f[7] = info_to_str(info)
                if vid in annotate_only:
                    locus, masked, unresolved = annotate_only[vid]
                    info = parse_info(f[7])
                    info['HERVK_LOCUS'] = locus['locus_id']
                    info['HERVK_LOCUS_N'] = locus['n_records']
                    if masked:
                        info['HERVK_MEMBERS_MASKED'] = ','.join(masked)
                    if unresolved:
                        info['HERVK_MEMBERS_UNRESOLVED'] = ','.join(unresolved)
                    f[7] = info_to_str(info)
                fh.write('\t'.join(f) + '\n')

    if args.out_archive:
        with open(args.out_archive, 'w') as fh:
            for h in head:
                fh.write(h + '\n')
            fh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
                     + '\t'.join(samples) + '\n')
            for r in archived:
                fh.write('\t'.join(r) + '\n')

    if args.report:
        write_report(args.report, consolidated, report, samples)

    late = {l for l, _ in skipped}
    n_ok = sum(1 for lid, st, _ in report if st == 'consolidated' and lid not in late)
    n_ann = sum(1 for _, st, _ in report if st == 'annotated')
    all_skips = [(l, w) for l, st, w in report if st == 'skipped'] + skipped
    sys.stderr.write(f'[hervk_reconcile] consolidated {n_ok} loci; '
                     f'{n_ann} annotated in place; {len(all_skips)} skipped\n')
    for lid, why in all_skips:
        sys.stderr.write(f'    SKIP {lid}: {why}\n')
    if mask_gt:
        sys.stderr.write(f'    genotypes withheld on {len(mask_gt)} tandem '
                         f'record(s): {", ".join(sorted(mask_gt))}\n')


def write_report(path, consolidated, report, samples):
    n = len(samples)
    with open(path, 'w') as fh:
        fh.write('# HERV-K locus consolidation\n\n')
        fh.write('Genotypes are resolved by ALT dosage across the member records. '
                 'Where the called members already account for every haplotype, an '
                 'uncalled member is pinned to zero -- arithmetic, not inference. '
                 'Nothing is filled in beyond that: unaccounted haplotypes stay '
                 'missing, which is why AC and AN are reported separately.\n\n')
        fh.write('| locus | alleles | AC | AN | 2N | resolved | partial | '
                 'ploidy exceeded | discovery AC/AN |\n')
        fh.write('|---|---|---|---|---|---|---|---|---|\n')
        for lid, c in consolidated.items():
            ac = ','.join(str(c['ac'].get(i + 1, 0))
                          for i in range(len(c['alt_states'])))
            disc = (','.join(str(c['disc_ac'].get(i + 1, 0))
                             for i in range(len(c['alt_states'])))
                    + f"/{c['disc_an']}") if c['disc_an'] else 'n/a'
            fh.write(f"| {lid} | {c['ref_state']} -> {','.join(c['alt_states'])} "
                     f"| {ac} | {c['an']} | {2*n} | {c['n_res']} | {c['n_part']} "
                     f"| {c['n_viol']} | {disc} |\n")
        viol = [l for l, c in consolidated.items() if c['n_viol']]
        if viol:
            fh.write('\n## Ploidy exceeded\n\nSummed ALT dosage above ploidy '
                     'means a third allele was flattened by `bcftools norm -m-`. '
                     'Those samples are set missing rather than guessed at.\n\n')
            for l in viol:
                fh.write(f"- `{l}`: {consolidated[l]['n_viol']} sample(s)\n")
        short = [l for l, c in consolidated.items() if c['an'] < 2 * n]
        if short:
            fh.write('\n## AN below 2N\n\nMembers were structurally uncalled '
                     '(no depth reported at all) for some samples, so those '
                     'haplotypes are missing. This depresses AN without '
                     'affecting AC: at chr6 the graph recovers every carrier and '
                     'only fails to confirm non-carriers.\n\n')
            for l in short:
                c = consolidated[l]
                fh.write(f"- `{l}`: AN={c['an']} of {2*n}, "
                         f"{c['struct_miss']} structurally-missing member calls\n")
        skipped = [(l, w) for l, st, w in report if st == 'skipped']
        if skipped:
            fh.write('\n## Skipped\n\n')
            for l, w in skipped:
                fh.write(f'- `{l}`: {w}\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    f = sub.add_parser('flag', help='group candidates into loci and flag them')
    f.add_argument('--calls', required=True,
                   help='call table from hervk_classify.py --calls-out, '
                        'covering every candidate in the discovery VCF')
    f.add_argument('--vcf-in', required=True,
                   help='human-subset VCF to annotate (the discovery VCF is '
                        'never written to)')
    f.add_argument('--vcf-out', required=True)
    f.add_argument('--loci-out', required=True)
    f.add_argument('--ref-state')
    f.add_argument('--window', type=int, default=1200,
                   help='max footprint gap within one locus (default: one LTR '
                        'plus tolerance)')
    f.set_defaults(func=cmd_flag)

    c = sub.add_parser('consolidate',
                       help='collapse flagged loci in the graph-genotyped VCF')
    c.add_argument('--genotyped-vcf', required=True)
    c.add_argument('--loci', required=True,
                   help='hervk_loci.tsv from the flag step')
    c.add_argument('--calls',
                   help='hervk_calls.tsv, for per-member allele states')
    c.add_argument('--discovery-vcf',
                   help='pangenome.human.vcf -- assembly genotypes, reported '
                        'beside the graph counts as an independent check')
    c.add_argument('--reference',
                   help='FASTA, needed only for loci with several alleles and '
                        'no deletion spanning the locus')
    c.add_argument('--genotyper', default='giraffe',
                   help="back end that wrote --genotyped-vcf; 'auto' "
                        'reads it from the VCF header, which is what a '
                        'run consolidating against an existing VCF wants')
    c.add_argument('--out-vcf', required=True)
    c.add_argument('--out-archive',
                   help='the member records removed, kept verbatim')
    c.add_argument('--report')
    c.set_defaults(func=cmd_consolidate)

    args = ap.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
