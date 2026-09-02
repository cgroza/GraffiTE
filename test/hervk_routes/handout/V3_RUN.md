# v3 run: copy-number alleles

Read `README_HPC_CLAUDE.md` for the general setup and `E2E_TEST.md` for how
stage E works. This file covers only what is different.

Pipeline is at **`feat/hervk-copy-number`** (PR #96), branched off `v1.1dev`
after #95 landed. Five commits.

## What changed

v2 could describe null, solo and provirus. It could not describe an HML-2 locus
that carries more than one proviral unit, and chr7:4,699,540-4,717,514
(7p22.1a, HERV-K108) is one: CHM13 holds two proviruses in tandem sharing a
central LTR, and across the 20 CaG samples one, two and three units all
segregate. v2 classified one of the three records at that locus, the one
carrying the rarest allele, and then withheld its genotypes.

| # | change | changes output? |
|---|---|---|
| 1 | `is_candidate()` asks the architecture table before the `n_hits` rule | **yes**, two chr7 records enter |
| 2 | `truncated_by_window()` no longer exempts `provirus`; the rescue loops, doubling the flank | **yes**, chr7 reads as a 2-unit array |
| 3 | `unit_structure()` counts units and measures the period | **yes**, new columns |
| 4 | `ARCH_INT_PERM` and `CNV_PERIOD`; states gain `prov_xN` | **yes**, `tandem_prov` is gone |
| 5 | `cluster()` joins by element overlap, not exact coordinates | **yes**, chr7 becomes one locus |
| 6 | genotype masking moves to `hervk_reconcile.py consolidate` | **yes**, discovery GTs survive |
| 7 | `hervk_candidates.vcf`, the full candidate set with genotypes | new file |

## Leave `HERVK_REF_RM_OUT` unset

Same instruction as re-run 2, different reason. Change 2 makes the rescue
iterative, so windows are cut at 12000, then 24000, then 48000 for anything
still touching an edge. A `.out` computed over the old windows cannot describe
the new ones, and reusing it would quietly reproduce the truncated chr7 array
this run exists to fix.

```bash
unset HERVK_REF_RM_OUT
```

## Cost

More reference masking than re-run 2, still far less than run 1.

Run 1 sent 26.4 Mb to RepeatMasker and died on the wall clock. Re-run 2 sent
about 1.1 Mb plus roughly 30 rescue windows. This run adds the `provirus`
records that were exempt before, and gives anything still at an edge a second
and third round.

Rough arithmetic, not a measurement: 1.1 Mb first pass, plus perhaps 1 to 2 Mb
of rescue windows, assuming most records settle after one or two rounds. We do
not know how many reach round three. If the masking phase runs past an hour,
stop and say so.

Report the number of windows re-cut per round. `hervk_ref_state` prints one
line per round to stderr, so there should be up to three:

```
hervk_ref_state: re-cutting N window(s) at flank=12000 (element reached the window edge)
hervk_ref_state: re-cutting N window(s) at flank=24000 (element reached the window edge)
```

## Run it

```bash
cd /xdisk/cgoubert/cgoubert/GraffiTE1.1/CaG
REVISION=feat/hervk-copy-number ./bootstrap.sh
unset HERVK_REF_RM_OUT
./preflight.sh
./run_hervk_test.sh
./bundle_results.sh
```

**Edit two lines in `INPUTS.env` before `preflight.sh`.** `bootstrap.sh` never
overwrites a populated `INPUTS.env`, by design, so the copy on the cluster keeps
re-run 2's values and drops the new template beside it as `INPUTS.env.new`:

```bash
REVISION="feat/hervk-copy-number"
OUTDIR="hervk_v3_run"
```

`RM_DIR`, `REFERENCE`, `TE_LIBRARY` and `GENOTYPED_VCF` stay as they are.
The separate `OUTDIR` keeps this run from colliding with re-run 2.

The `REVISION=` in front of `bootstrap.sh` is needed on the first run only. The
copy of `bootstrap.sh` already on the cluster still defaults to
`v1.1dev-hervk-v2` and does not read `INPUTS.env`; the version it pulls fixes
both.

The `-resume` guard still applies: you should see the "pipeline moved" line. If
you do not, and `$OUTDIR` already holds output, stop and say so.

## Nextflow was never executed against these changes

There is no cluster on the machine the code was written on, so we reviewed
`module/main.nf` and `nextflow.config` by eye and nothing else. A first-launch
failure most likely sits in one of:

- `params.hervk_mask_graph_gt_at_cnv`, new, replacing `hervk_mask_tandem`
  (kept as a `null` alias and read in the `hervk_reconcile` process)
- the new `path("hervk_candidates.vcf")` output on `hervk_annotate`
- the `--vcf-out-candidates-only` flag on the first `hervk_classify.py` call

If it dies inside a minute with a Groovy or parameter error, that is where to
look. Report the stack trace verbatim rather than patching around it.

## What must be true

**`chr7-4706809-DEL-8503_108500` must appear in `hervk_calls.tsv`.** It is the
common allele at 7p22.1a, 26 of 40 haplotypes in 16 of 20 samples, and no
previous run has classified it.

| check | expected |
|---|---|
| `chr7-4706809-DEL-8503_108500` | present, `copy_number`, `prov_x2` -> `provirus`, `j=1236` |
| `chr7-4699715-INS-8504_108498` | present, `copy_number`, `prov_x2` -> `prov_x3`, `k=794` |
| `chr7-4700334-INS-8504_108499` | `copy_number`, `prov_x2` -> `prov_x3`, `k=174` |
| chr7 locus | all three records in one locus, flagged `MERGE_CANDIDATE` |
| chr7 `ref_elem_end` | 4717514, **not** POS+12000 (4712333 or 4711714) |
| chr7 `ref_n_units` / `ref_unit_bp` | 2 / 8504 |
| chr6:78,894,316 `ref_n_units` / `ref_unit_bp` | 1 / 8465 |
| chr12:133,148,144 `ref_n_units` / `ref_unit_bp` | 1 / ~4935 |
| chr6-78894876, chr12-133148145 | `copy_number`, `provirus` -> `prov_x2` |
| chr6-78894317-DEL-8465 | unchanged: `REF_ANNOT`, `provirus` -> `solo` |
| `copy_number` records | at least 5 |
| `null_prov` | still 3 |
| `ARCH_PERM` k values | 102, 174, 303, 410, 574, 707, 735, 868 unchanged, with 794 joining |
| `OVERSIZE_ELEMENT` on the chr7 array | absent |
| discovery genotypes on copy-number records | **present, not `./.`** |
| `3_TSD_search/hervk_candidates.vcf` | exists, holds the chr7 records with genotypes |
| chr7 locus flags | `LOCUS_SPLIT_BY_HUMAN_FILTER` |

That last one is correct behaviour, not a fault. The `--human` pME filter is
deliberately narrower than the HERV-K candidate list and was left alone, so two
of the three chr7 records stay out of `pangenome.human.vcf`.
`hervk_candidates.vcf` holds their annotation and genotypes, and the flag
records that the locus lost members.

The discovery-genotype row is the one worth checking by hand. In v2 those
records were `./.` in `pangenome.human.vcf`; masking now happens only in stage
E, on the graph calls. If they still come back `./.`, change 6 did not take.

## What will move, and should

Do not treat these as failures.

- `tandem_prov` is gone from every class count. `copy_number` replaces it.
- `candidates classified` rises above re-run 2's 131. The `n_hits` gate no
  longer drops records, and how many that is worth is not known in advance.
  Report the new number.
- Two new evidence codes appear, `CNV_PERIOD` and `ARCH_INT_PERM`.
- Locus and record counts shift at chr7: two more records, still one locus.
- `partial` may drop again, since change 2 rescues windows re-run 2 left
  alone. It may also not move at all; re-run 2 already cleared most of them.

## Report back

The bundle, plus:

- the assertion log
- the per-round rescue lines from `hervk_ref_state` (stderr, in the Nextflow log)
- `ref_state` and `ref_n_units` distributions
- the full `copy_number` block the assertion script prints, with unit counts,
  periods, k and j for each record
- `candidates classified`, `loci`, and the evidence distribution
- any `HERVK_NOTE` that is not `.`, in full
- wall time for `hervk_annotate` from `nextflow_trace.txt`

Do not tune parameters to make assertions pass. A disagreement is the result,
and the three that matter most (`j=1236`, `ref_unit_bp=8504`,
`ref_elem_end=4717514`) were each measured off chm13v2.0 by hand before the
code was written.
