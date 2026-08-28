# HERV-K v2 test — run state as of 2026-08-28 (E2E, stage E)

Status: **COMPLETE, PASS.** Pipeline commit `22c7395`, driver job 23713404,
exit 0, 33m44s. Bundle `hervk_v2_results_20260827.tar.gz` (7.7 MB, not in git)
is the 2026-08-27 run; its numbers are identical to this one.

This supersedes the re-run 2 state. **The three `REF_ARCH_CONFLICT` that the
previous version of this file flagged as unresolved are resolved** — see below.

## Result

```
HERV-K v2 assertions — 41 records checked against 41 expectations
  tandem_prov           : 3  chr12-133148145-INS-4933_23163,
                             chr6-78894876-INS-8465_106221,
                             chr7-4700334-INS-8504_108499
  stage E: consolidated : 3  (HERVK_chr11_101704640, HERVK_chr12_55299985,
                              HERVK_chr1_75219429)
  stage E: GT withheld  : 1
  candidates classified : 131
  loci                  : 55  (5 flagged for merge)
  null_prov             : 3  chr19-21797327-INS-9478_60578,
                             chr19-22370220-INS-8229_60593,
                             chr8-7226885-INS-9468_114234
  evidence              : ARCH_2LTR=2, ARCH_PERM=8, ARCH_SOLO=30,
                          NON_HML2=71, REF_ANNOT=20
  NOTE: 1 loci split by the --human FILTER="PASS" requirement: HERVK_chr8_7552031
  PASS — all expectations met.
```

## Resolved: the three REF_ARCH_CONFLICT were tandem duplications

Re-run 2 left three records whose architecture said REF = solo while the masked
reference read `provirus`. RERUN_2.md predicted zero such conflicts and carved
out only chr6; the other two had no DEL partner to corroborate them, so the
output could not settle whether they were the same phenomenon or a genuine
disagreement.

`cacf5f0` (call tandem duplications and withhold their genotypes) answers it.
All three are tandem duplications, and all three now carry
`TANDEM_DUP,GT_MASKED`:

```
chr12-133148145-INS-4933_23163   TANDEM_DUP,GT_MASKED   (k=735, LTR5A)
chr6-78894876-INS-8465_106221    TANDEM_DUP,GT_MASKED   (k=410)
chr7-4700334-INS-8504_108499     TANDEM_DUP,GT_MASKED   (k=174)
```

**`REF_ARCH_CONFLICT` count is now 0**, which is what RERUN_2.md expected. The
remaining notes are three `INS_INTO_NONEMPTY_REF` (`chr21-13434-INS-318_75608`,
`chr5-26290893-INS-60_99699`, `chr8-7510479-INS-69_114254`).

Class counts moved accordingly — the 3 tandem records came out of `solo_prov`
(2) and `truncated_prov` (1); the total is still 131:

| Class | re-run 2 | now |
|---|---|---|
| null_solo | 31 | 31 |
| solo_prov | 11 | 9 |
| truncated_prov | 8 | 7 |
| null_prov | 3 | 3 |
| tandem_prov | — | 3 |
| other | 78 | 78 |

`ref_state` unchanged: null 123, solo 43, partial 5, provirus 15.

## Stage E — all seven E2E_TEST.md checks match

| check | expected | observed |
|---|---|---|
| consolidated | 3 (chr1, chr11, chr12) | 3 |
| annotated in place | 1 (chr6) | 1 |
| skipped | 1 (chr8) | 1, `1 of 2 members present` |
| chr11 | `AC=23 AN=40`, `HERVK_DISC_CONCORDANT` | both |
| chr12 | `AC=6,21 AN=33`, ploidy exceeded 2 | plus `N_RESOLVED=15`, `N_PARTIAL=3`, discovery `8,21/40` |
| chr6 tandem | `HERVK_GT_MASKED`, all `./.` | 20 genotypes, 0 non-missing |
| accounting | in − 6 + 3 = out | 5808 − 6 + 3 = 5805 |

Consolidation table:

```
HERVK_chr1_75219429   solo -> provirus        AC=18   AN=40  2N=40  discovery 20/40
HERVK_chr11_101704640 solo -> provirus        AC=23   AN=40  2N=40  discovery 23/40
HERVK_chr12_55299985  solo -> null,provirus   AC=6,21 AN=33  2N=40  discovery 8,21/40
```

Outputs in `hervk_v2_run/4_Genotyping/`:
`GraffiTE.merged.genotypes.human.vcf.gz` (+ `.tbi`),
`hervk_unconsolidated_records.vcf`, `hervk_reconciliation_report.md`.

## Three fixes this round, all pushed

| commit | what it unblocked |
|---|---|
| `04d8970` | `human_hervk_ids = "^HERVK$..."` was double-quoted, so Groovy read `$` as an interpolation and the whole pipeline config failed to parse. **Every** run on the branch died at launch in ~6s. Single-quoted. |
| `a535dd4` | `bundle_results.sh` predated stage E and collected only `3_TSD_search/` and `$OUTDIR/`, so the bundle silently omitted all three stage-E outputs — the entire point of the test. |
| `22c7395` | `main.nf` passed `params.graph_method` to `hervk_reconcile` on the external-VCF path too. That param describes *this* run's genotyping, which `--genotype false` disabled, so its `pangenie` default rejected a valid `vg call` VCF. The path now passes `auto` and `detect_genotyper()` reads the back end from the header. |

`22c7395` was validated by reverting the local `graph_method="giraffe"` workaround
back to the stock `pangenie` and rerunning, so job 23713404 reproduced the exact
configuration that failed. It logged
`[hervk_reconcile] --genotyper auto: header says giraffe` and passed, with every
stage-E number identical to the run before it. `--genotyper` feeds only the
guard (it appears nowhere else in `hervk_reconcile.py`), so it cannot change
consolidation output.

## Local environment — not in git, differs from a fresh checkout

- `INPUTS.env` carries `GENOTYPED_VCF=.../june_out/4_Genotyping/GraffiTE.merged.genotypes.vcf.gz`.
  `bootstrap.sh` preserves a populated `INPUTS.env` and writes new fields only to
  `INPUTS.env.new`, so **anyone upgrading from an earlier run has no
  `GENOTYPED_VCF` at all** — and empty means "skip stage E entirely", i.e. the
  run passes while testing nothing. Worth a preflight warning.
  `june_out` is deliberate: same run as `RM_DIR`, so discovery and genotyped
  candidate sets match. E2E_TEST.md's expectations only hold for a matched pair.
- `nextflow.config`: `executor.$local.cpus/memory` pinned to 16 / 80 GB
  (Nextflow otherwise sizes the local pool from the node's 94 cores / 470 G, not
  the job cgroup). `env.HERVK_REF_RM_OUT` commented out. `graph_method` is back
  at the stock `pangenie` and no longer matters after `22c7395`.
- `PROFILE="standard"` (local executor), therefore **run via
  `sbatch submit_hervk_driver.sh`, not `./run_hervk_test.sh`**. The handout's
  foreground form assumes the cluster profile; with the local executor it would
  run every task on the login node, which has no `singularity`. The cluster
  profile is not an option here: Puma charges ~85 min of queue latency *per
  task*, and there are 1337 `tsd_search` tasks.

## History

| job | outcome |
|---|---|
| 23613713 | TIMEOUT 12 h — cluster profile, ~85 min queue wait per task |
| 23667884 | killed — `tile_hits` O(n_hits*span) blowup |
| 23668829 | failed 3 s — Nextflow refuses `-latest` on a dirty asset repo |
| 23668844 | TIMEOUT 4 h — inside RepeatMasker ProcessRepeats |
| 23670310 | ProcessRepeats rerun standalone, 2h14m, COMPLETED |
| 23671057 | run 1 PASS (`28109e7`) |
| 23684967 | cancelled 2h02m — `max_svlen` capped at classify, not candidacy |
| 23687175 | run 2 PASS (`a507cc8`), 34m13s |
| 23703904 | FAILED 6 s — unparseable pipeline config |
| 23703909 | FAILED 32m — `--genotyper pangenie` rejected |
| 23704225 | E2E PASS (`04d8970` + local giraffe override), 5m41s |
| 23713404 | **E2E PASS (`22c7395`), 33m44s — auto-detection, no override** |

## Leftovers, safe to delete

`hervk_resume/` (64.5 MB, run 1's `ref_windows.fa.cat.gz` and derived `.out`)
and `work/hervk_resume/`. Only needed if re-enabling `HERVK_REF_RM_OUT`.
