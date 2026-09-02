# HERV-K classifier v3: copy-number alleles

Plan, 2026-09-02. Successor to `HERVK_CLASSIFIER_V2.md`. Nothing here is
implemented yet.

## Context

The v2 classifier resolves three allele states at an HML-2 locus: `null`,
`solo` and `provirus`. We audited the CaG results for the paper and found that
state space too small. One locus in the call set is a copy-number array the
classifier cannot see.

chr7:4,699,540–4,717,514 in CHM13v2 is **7p22.1a / HERV-K108**, and the
reference carries two complete proviruses in tandem sharing a central LTR:
`LTR-INT-LTR-INT-LTR`, units 99.778% identical, period 8,504 bp. Three alleles
segregate in the 20 CaG samples:

| record | state | haplotypes | what v2 does |
|---|---|---|---|
| `chr7-4706809-DEL-8503_108500` | 1 unit | 26/40, in 16 samples | never classified |
| (reference) | 2 units | 12/40 | |
| `chr7-4700334-INS-8504_108499` | 3 units | 1/40 | `tandem_prov`, genotypes masked |
| `chr7-4699715-INS-8504_108498` | 3 units | 1/40 | never classified |

So the locus contributes nothing: its common allele never reaches
`hervk_calls.tsv`, and the record that does carries the rarest allele and is
then masked. Pasternack, Paulsen & Nath 2025 (*Front Genet* 16:1498978) report the same
loci at the same frequencies: an 8,501 bp deletion in ~80% of 222 individuals
at 7p22.1a, against 16 of 20 here, and a 3-LTR tandem at 6q14.1 in ~5%, against
1 of 20 at our chr6 locus. Their 6q14.1 element is 9,422 bp against our 9,425,
so it is almost certainly the same locus. We did not lift over to confirm.

The same architecture appears at chr6:78,894,316 (1 unit → 2) and
chr12:133,148,144 (1 unit → 2). Full evidence in
`paper/HERVK_FINAL/NOTE_tandem_alleles.md`.

**Goal:** make copy number a first-class allele state, so the classifier
describes an HML-2 array by its unit count instead of forcing it into
`provirus` plus a `tandem_prov` escape hatch.

## Constraint: the re-run has no graph genotyping

Re-runs use the existing handout path (`test/hervk_routes/handout/`), which
already supports this. `--RM_dir` supplies a published `2_Repeat_Filtering`,
`--genotype false` skips graph construction and `vg call` entirely, and
`--hervk_reconcile_vcf` points stage E at a genotyped VCF from an earlier run
(`main.nf:268-279`, outside the `if(params.genotype)` block). Inputs stay the
same: `RM_DIR`, `REFERENCE`, `TE_LIBRARY`, `GENOTYPED_VCF`.

Two consequences for the design:

1. **Discovery genotypes are the deliverable.** They come from haplotype-resolved
   assembly alignments and are sound at these records. The graph genotypes are
   not: at a copy-number locus, reads from the pre-existing reference copy
   traverse the duplicated path, so non-carriers acquire ALT support. At chr6
   the ALT fraction tracks provirus dosage rather than carriage: provirus
   homozygotes 0.13 to 0.42, solo-LTR carriers 0.00. Move masking from the
   classifier to the reconciler, and justify it on mapping ambiguity.
2. **`hervk_ref_state.py` re-masks the reference on every run.** Widening its
   windows is affordable, so the fixes in change 2 cost minutes, not hours.

`bin/` is not part of Nextflow's task hash, so `-resume` silently replays old
HERV-K output when only a script changes. `run_hervk_test.sh:41-64` already
guards this with a `.last_commit` stamp. Keep using it.

## What we found, and where

Everything below was verified against the CaG run and `chm13v2.0.fa`.

### 1. The candidacy gate drops records on a RepeatMasker naming quirk

`bin/hervk_classify.py:161` requires `n_hits == 1`, or `n_hits == 2` with
`Retroposon/SVA`. Both missing chr7 records have `n_hits = 3`. The third hit is
an artefact of how RepeatMasker named a split LTR fragment: 792 bp came back
`LTR5_Hs` where the surviving record's 172 bp fragment came back `SVA_A`, which
`hervk_arch.py` reassigns as SINE-R, after the gate has already run.

This gate is separate from the `--human` filter. The candidate ID list at
`module/main.nf:298` is far looser (`matching_classes="LTR/ERVK" &
abs(SVLEN)<=25000`) and admits all three chr7 records, which is why they appear
in `hervk_arch.tsv` and `hervk_refstate.tsv` but not in `hervk_calls.tsv`.
**The tables already hold the evidence; only `is_candidate` refuses it.**

How much the candidate set grows, we cannot say in advance. The local discovery
VCF holds 208 records carrying `LTR/ERVK`, against 131 candidates in v2's
assertion log, and we could not reconcile those two counts against each other.

### 2. The reference window truncates the array, and a guard suppresses the fix

`truncated_by_window` (`bin/hervk_ref_state.py:257-274`) ends with

```python
return at_edge and result['state'] != 'provirus'
```

An element that reads `provirus` is never re-cut, even when it runs into the
window edge. chr7 reads `provirus` at `flank=12000`, touches the edge, and is
not rescued. The reported element ends at exactly `POS + 12000` for both
records, 4,712,333 and 4,711,714, so the third LTR and the second unit were
never masked.

The rescue is also single-pass (`:398-414`), and 12 kb cannot hold a 17,975 bp
array anchored 794 bp inside it.

### 3. Locus clustering compares reference elements by exact coordinates

`cluster` (`bin/hervk_reconcile.py:158-189`) joins two records when their
`(ref_elem_chrom, ref_elem_start, ref_elem_end)` tuples are **equal**, or when
their footprints lie within `--window` (default 1200,
`nextflow.config:99-102`). Because the element span is window-dependent, the
two chr7 insertions get different tuples and join only by the 619 bp window
test. The deletion sits 6,475 bp away and would form its own locus even if it
were classified.

### 4. The deletion carries its own diagnostic signature, currently unused

`hervk_arch.py` already tiled the chr7 deletion correctly:

```
chr7-4706809-DEL-8503_108500   ARCH_NONE   ltr_bp 968   int_bp 7535
arch: INT:1237-7536/LTR:1-968/INT:1-1236
```

The internal region is split across the termini with complementary consensus
intervals (`1237-7536` and `1-1236`, meeting at j = 1236) and one complete LTR
sits between them. That is a whole proviral unit in circular permutation, with
the permutation point in the **internal region** rather than the LTR.

The genome confirms it. The deletion starts 6,301 bp into the first internal
region, and the element is on the minus strand, so the consensus position is
7536 − 6301 + 1 = **1236**. The same arithmetic that validates `k` for
`ARCH_PERM` validates `j` here.

This signature is worth more than the record it came from. `ARCH_PERM` and this
one are the two crossover positions of the same reaction:

| signature | split across termini | whole in the middle | crossover in |
|---|---|---|---|
| `ARCH_PERM` | one LTR | internal region | the LTR |
| `ARCH_INT_PERM` (new) | internal region | one LTR | the internal region |

An `ARCH_INT_PERM` signature is only possible when the reference already holds
two or more units, because it needs internal-region homology. It therefore
identifies a multi-unit array **from the SV sequence alone**, with no reference
masking.

## The changes

In dependency order. Base the work on **`v1.1dev`**: PR #94 merged
`v1.1dev-hervk-v2` into it on 2026-08-28 (merge commit `a8d7cc2`), so all four
HERV-K scripts and the handout are on the mainline. Line numbers below are
against `origin/v1.1dev`.

**Two fixes never made it into that merge.** Both were written in response to
the Copilot review on #94 and committed to `v1.1dev-hervk-v2` about six hours
after the merge landed. No follow-up PR exists:

| commit | fix | consequence on `v1.1dev` today |
|---|---|---|
| `1415680` | guard SINE-R reassignment against missing consensus coordinates | `hervk_ref_state.py --rm-annotation` (the BED4 path) raises `TypeError` on a `None < int` comparison |
| `a487582` | replace `SVTYPE`/`SVLEN` header definitions instead of appending them | every consolidated VCF carries two conflicting definitions of each, which is invalid VCF and makes vcfR warn on every read |

They also carry `test/hervk/test_arch.sh` and six assertions in
`test_consolidate.sh`. Land them first, as a small PR of their own, before any
v3 work. Nothing in this plan depends on them, but the v3 branch should not be
cut from a base that is missing them.

### C1. Let the architecture decide candidacy (`bin/hervk_classify.py`)

Replace the `n_hits` rule in `is_candidate` (`:145-170`) with a test on the
architecture table, which `process_vcf` already has in scope as `arch_tbl`.

Admit an SV when `LTR/ERVK` is among its classes, `|SVLEN|` is within
`max_svlen`, and its `hervk_arch.tsv` row shows `ltr_bp + int_bp >=
min_hml2_bp` and covers at least `hml2_frac_min` (propose 0.80) of `|SVLEN|`.
Keep the current `n_hits` rule as the fallback for records with no arch row.

Pass `arch_tbl` into `is_candidate`; update the docstring, which currently
claims the gate is identical to the `--human` carve-out.

**Leave `module/main.nf:429-447` alone.** The `--human` pME filter defines the
paper's TE set and its counts, and widening it would move numbers that have
nothing to do with HERV-K. The calls and loci tables are already built from the
full discovery VCF, so they become complete on their own. `n_in_human` and
`LOCUS_SPLIT_BY_HUMAN_FILTER` (`hervk_reconcile.py:217-220`) already record the
difference.

Add `--vcf-out hervk_candidates.vcf` to the full-candidate classify call
(`module/main.nf:315-318`) so the annotated HERV-K records and their discovery
genotypes land in one file, instead of only reaching the human subset.

### C2. Recover the whole array (`bin/hervk_ref_state.py`)

- Drop the `and result['state'] != 'provirus'` clause from
  `truncated_by_window` (`:274`). A provirus at a window edge may be a
  truncated array, which is exactly the case that matters.
- Make the rescue iterative, or size the flank from the element found in the
  first pass, such as `max(rescue_flank, 2 * elem_span + flank)`, and re-check
  the edge. Cap the number of rounds.
- `max_element_span` (`:66`, 10500) currently flags a merged array as
  `OVERSIZE_ELEMENT`. Keep the flag for genuinely suspect merges, but do not
  raise it when the span is a clean multiple of a unit period (C3).
  `cluster_elements`' `element_gap = 1000` is correct here and should stay: it
  is what merges array units into one element.

### C3. Count units (`bin/hervk_ref_state.py`)

Parse the ordered fragment list already in `ref_arch` and emit four columns:

| column | meaning |
|---|---|
| `ref_n_units` | complete `LTR-INT` units, counting a shared LTR once |
| `ref_unit_bp` | LTR-to-LTR period, internal region plus one LTR |
| `ref_array_start` / `ref_array_end` | full array extent |

A clean array has `n_LTR = n_INT + 1`; the period is the element span minus one
LTR consensus length, divided by `n_INT`. Set `ref_state = provirus` with
`ref_n_units = N`, keeping `null`, `solo` and `partial` unchanged. Expected for
chr7 once C2 lands: `ref_n_units = 2`, `ref_unit_bp = 8504`, array
4,699,540–4,717,514.

### C4. Copy-number classification (`hervk_arch.py`, `hervk_classify.py`)

**New signature in `hervk_arch.py`.** `ARCH_INT_PERM`: internal-region
fragments at both termini with complementary consensus intervals
(`[j+1..7536]` and `[1..j]`, within `perm_tol`) and one complete LTR between
them. Emit `j` alongside `k`.

**New evidence code in `hervk_classify.py`.** `CNV_PERIOD`: `|SVLEN|` is within
tolerance of an integer multiple of `ref_unit_bp` (propose ±50 bp or ±0.5%,
whichever is larger). This is the general test, and it covers deletions, which
`ARCH_PERM` alone cannot. Verified against all four records: chr6 8,465 against
a period of 8,465; chr7 8,504 and 8,503 against 8,504; chr12 4,933 against
4,935.

**Allele states.** Extend the state space from `provirus` to `prov_xN`.
Insertions raise the unit count, deletions lower it:

| record | REF | ALT |
|---|---|---|
| `chr7-4706809-DEL-8503` | `prov_x2` | `prov_x1` |
| `chr7-4700334-INS-8504` | `prov_x2` | `prov_x3` |
| `chr6-78894876-INS-8465` | `prov_x1` | `prov_x2` |
| `chr12-133148145-INS-4933` | `prov_x1` | `prov_x2` |

`prov_x1` is today's `provirus`; keep that spelling so existing consumers and
fixtures do not break, and introduce `prov_xN` only for N ≥ 2. Retire
`tandem_prov` as a class, replacing it with `prov_x2`. Add
`HERVK_LOCUS_TYPE=copy_number` for loci where unit count varies.

The evidence ladder gains one rung, above the reference check because it does
not depend on masking:

1. `ARCH_2LTR`, `ARCH_PERM`, `ARCH_INT_PERM`: architecture
2. `CNV_PERIOD`: period arithmetic against `ref_unit_bp`
3. `REF_ANNOT`: the masked reference
4. size confidence, which never decides

Record `HERVK_K` or `HERVK_J` and the predicted junction offset in the
reference unit, so the genome cross-check that validated all four records
becomes an automatic assertion rather than a manual one.

### C5. Cluster by overlap (`bin/hervk_reconcile.py`)

In `cluster` (`:158-189`), replace tuple equality on the reference element with
**interval overlap**, and take the join window from the array span
(`ref_array_end - ref_array_start`) when one is known, falling back to
`--window` otherwise. This is what puts all three chr7 records in one locus.
Raise `hervk_locus_window` only if C2 leaves a case that overlap does not
catch.

### C6. Move masking from the classifier to the reconciler

`hervk_classify.py:493-503` blanks genotypes for `tandem_prov`. That is the
wrong place: it blanks the discovery genotypes too, and those are the only ones
a no-genotyping re-run produces.

- Report discovery AC/AN for every copy-number allele.
- Mask graph genotypes at copy-number loci inside `hervk_reconcile.py
  consolidate`, and emit the diagnostic that justifies it: the spread of ALT
  fraction across samples the discovery call says are non-carriers.
- Rename `hervk_mask_tandem` (`nextflow.config:92-97`) to
  `hervk_mask_graph_gt_at_cnv`, keeping the old name as a deprecated alias.

### C7. Wording

The `mask_tandem` comment (`hervk_classify.py:70-76`) and
`paper/HERVK_FINAL/HERVK_RESULTS_CaG.md` both say a tandem unit is "neither
transposition nor intra-element recombination", and both call the tandem
records singletons at all three loci. Unequal LTR–LTR exchange **is**
intra-element recombination, and chr7 is 2/40 at a locus where a third allele
sits at 26/40. Rewrite against `NOTE_tandem_alleles.md` §5 and §7, and cite
Hughes & Coffin 2004 (*PNAS* 101:1668-1672) rather than presenting the
mechanism as new. Their wording is still unverified; check it before quoting.

### C8. Fixtures

`test/hervk/` and `test/hervk_routes/handout/assert_hervk_test.py` already hold
the harness. Add:

- an `ARCH_INT_PERM` fixture built from the chr7 deletion's real tiling
- a three-state chr7 locus fixture, asserting one locus and three records
- a `ref_n_units = 2` refstate fixture
- a regression asserting `truncated_by_window` fires on an edge-touching
  `provirus`

## What the re-run should produce

Falsifiable expectations, measured from the current data.

| check | expected |
|---|---|
| chr7 locus | 1 locus, 3 records, states `prov_x1`/`prov_x2`/`prov_x3` |
| chr7 discovery AC/AN | 1 unit 26/40, 2 units 12/40, 3 units 2/40 |
| chr7 `ref_n_units` | 2, `ref_unit_bp` 8504, array 4,699,540–4,717,514 |
| chr6 locus | provirus 31/40, solo 8/40, `prov_x2` 1/40, carrier sets disjoint |
| chr12 locus | `prov_x1` 39/40, `prov_x2` 1/40 |
| `ARCH_INT_PERM` records | at least 1, `j = 1236` at chr7 |
| `CNV_PERIOD` agreement | `|SVLEN|` within 50 bp of the period at all four records |
| junction cross-check | predicted offset within 20 bp of observed at every `ARCH_PERM`/`ARCH_INT_PERM` record |
| `ARCH_PERM` k values | 102, 174, 303, 410, 574, 707, 735, 868 unchanged |
| `null_prov` | still 3 |
| `partial` count | no higher than re-run 2 |
| loci total | 55 or more; report any that merged or split against re-run 2 |

Two things to watch and report rather than tune:

- **Loci split by the `--human` filter.** C1 makes the calls table complete
  while the human VCF stays as it was, so `LOCUS_SPLIT_BY_HUMAN_FILTER` should
  appear at chr7 and possibly elsewhere. That is the correct outcome, not a
  fault.
- **New `prov_xN` calls at loci we have not looked at.** 208 records carry `LTR/ERVK` in the local discovery VCF, against 131 classified in v2. Any new
  copy-number locus needs the same genome cross-check before we believe it.

## Verification

Sequence-level checks, independent of the pipeline:

```bash
# duplicated unit is a copy of the reference starting at the breakpoint
blastn -query <alt_allele.fa> -db <ref_window> \
       -outfmt "6 qstart qend sstart send pident length mismatch gapopen"

# array architecture: LTR and internal blocks over the locus
blastn -query <chr6_LTR.fa|chr6_INT.fa> -db <ref_window> -outfmt "6 sstart send pident length"
```

Pipeline run, unchanged from the handout except for the revision:

```bash
cd /xdisk/cgoubert/cgoubert/GraffiTE1.1/CaG
./bootstrap.sh                 # pull the v3 revision; INPUTS.env is preserved
unset HERVK_REF_RM_OUT         # windows change in C2; a stale .out cannot describe them
./preflight.sh
./run_hervk_test.sh            # confirm it prints the "pipeline moved" line
./bundle_results.sh
```

Report the assertion log, every `HERVK_NOTE` that is not `.`, the `ref_state`
and `ref_n_units` distributions, how many windows the rescue re-cut, and the
`hervk_annotate` wall time from `nextflow_trace.txt`.

## Out of scope

- Re-genotyping in the graph. A copy-number-aware genotyper for these loci is a
  separate problem, and the paper does not need one.
- The `--human` pME filter. C1 deliberately routes around it.
- Proviral arrays with no null, solo or copy-number variation. Clément's
  standing call: of no interest without a segregating state.
- Deciding assembly artefact against real allele. That needs per-base read
  coverage over the duplicated unit and junction-spanning reads in the carrier
  assemblies, neither of which this pipeline produces.
