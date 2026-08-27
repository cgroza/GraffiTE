# HERV-K v2 test — run state as of 2026-08-25 (re-run 2)

> Historical record of re-run 2, kept for its measurements and job history.
> Its open question — the three `REF_ARCH_CONFLICT` records — has been
> answered: they are tandem duplications, now classified `tandem_prov` with
> genotypes withheld. For what to run and expect now, see `E2E_TEST.md`.

Status: **COMPLETE, PASS.** Pipeline commit `a507cc8`.
Bundle: `hervk_v2_results_20260825.tar.gz` (2.7 MB, not in git).
Driver job 23687175, exit 0, 34m13s total; `hervk_annotate` 1m45s / 264 MB.
Results under `hervk_v2_run/3_TSD_search/`.

Run 1 (2026-08-24, commit `28109e7`) also passed; its bundle is
`hervk_v2_results_20260824.tar.gz` and its assertion log is kept at
`hervk_v2_run/hervk_assertions.log.run1`. Run 1's numbers are superseded —
three of the five fixes in `0681c33` change output.

## Result

```
HERV-K v2 assertions — 41 records checked against 41 expectations
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

`ref_state` distribution, rescue pass re-cut **13** windows at `flank=12000`:

| ref_state | run 1 | run 2 |
|---|---|---|
| null | 134 | 123 |
| solo | 45 | 43 |
| partial | 22 | **5** |
| provirus | 7 | **15** |

Classes: `null_solo` 31, `solo_prov` 11, `truncated_prov` 8 (was 7),
`null_prov` 3, `other` 78 (was 79).

## THE THING TO ANALYSE: three REF_ARCH_CONFLICT, not zero

`RERUN_2.md` predicts **zero** `REF_ARCH_CONFLICT` and carves out only
`chr6-78894876` as expected. There are **three**. Its own instruction: "Any
conflict remaining is a genuine disagreement and is the most interesting thing
in the output — report it in full."

```
chr12-133148145-INS-4933_23163  k=735  REF_ARCH_CONFLICT:arch=solo,ref=provirus
  arch=LTR:736-1033/INT:1-707/INT:3007-3690/INT:5187-5462/INT:5587-7536/LTR:1-735
chr7-4700334-INS-8504_108499    k=174  REF_ARCH_CONFLICT:arch=solo,ref=provirus
  arch=LTR:175-968/INT:1-7536/LTR:~174bp
chr6-78894876-INS-8465_106221   k=410  REF_ARCH_CONFLICT:arch=solo,ref=provirus   <- expected
  arch=LTR:411-968/INT:1-7536/LTR:1-410
```

All three share one signature: `ARCH_PERM` architecture (implying REF = solo)
against a masked reference reading `provirus`. What separates them is
corroboration:

- `chr6-78894876` sits in a **2-record locus** with its DEL partner
  (`HERVK_chr6_78894316`, flags `MERGE_CANDIDATE,POLARITY_CONFLICT`), so
  `chr6-78894317-DEL-8465` independently confirms the reference.
- `chr12-133148145` and `chr7-4700334` are **single-record loci, no flags** —
  no partner to cross-check. Both appear in `EXPECTED.tsv`.

Unresolved: whether these two are the same mis-polarised-INS phenomenon as
chr6, or a real disagreement. Without a DEL partner the two cannot be
distinguished from the output alone. Note `chr12-133148145` is on **LTR5A**
(1033 bp consensus, not 968) and is one of the records behind `truncated_prov`
going 7 -> 8.

Other notes (3, all `INS_INTO_NONEMPTY_REF`): `chr21-13434-INS-318_75608`,
`chr5-26290893-INS-60_99699`, `chr8-7510479-INS-69_114254`.

Also weaker than predicted: `other` fell only 79 -> 78, not to zero.

## Why a507cc8 exists — fix 4 shipped in the wrong stage

`0681c33` added `max_svlen = 25000` to `bin/hervk_classify.py` DEFAULTS
(line 75, used at 153). `RERUN_2.md` describes fix 4 as a cap "on candidacy".
It was not: `module/main.nf` built `hervk_candidate.ids` with an uncapped
`bcftools view -H -i 'matching_classes="LTR/ERVK"'`, and **that file is what
`hervk_ref_state.py` masks against**. So the cap dropped the record from the
calls only after its reference window had already been masked.

Consequence: re-run 2's first attempt (job 23684967) masked the same 26.4 Mb as
run 1 and was killed at 2h02m on the same trajectory that timed out job
23668844. `RERUN_2.md`'s "under 2 Mb / minutes" depends entirely on the half
that did not ship.

`a507cc8` moves the cap into the candidacy query. Measured effect:

| | before | after |
|---|---|---|
| candidates | 208 | 186 |
| masking input | 26,714,672 B | 865,486 B |
| largest window | 25,267,468 bp | 26,391 bp |
| `hervk_annotate` | timed out (>2 h) | 1m45s |

Two caveats, both in the commit message: it drops **22** candidates, not just
the one artifact (everything over 25 kb — a complete HML-2 provirus is 9,472 bp,
so this is defensible but larger than "drops junk" implies); and 25000 is now
duplicated in `module/main.nf` and `hervk_classify.py`, which is exactly how it
drifted. Upstream should promote it to one shared param.

## The record that caused all of it

`chr1-120594342-DEL-25264467` — a PAV alignment artifact, not a variant:
`FILTER=COMPOUND`, `CALL_SOURCE=ALNTRUNC`, **zero-length** `QRY_REGION`,
REF->ALT 25,264,468 bp -> 1 bp, window spanning the chr1 centromere, 84.3%
satellite. It qualified because `matching_classes` is `Number=.` and bcftools
`=` matches if ANY element does — 12 of its 18,667 fragments are LTR/ERVK
(0.064%). Now excluded at candidacy.

## Local environment — differs from the handout on purpose

These are in the working directory, not in git:

- `INPUTS.env`: `PROFILE="standard"` (local executor). The cluster profile costs
  ~85 min of queue latency **per task** on Puma; 1337 `tsd_search` tasks made a
  12 h driver die at 57% (job 23613713). Local executor in one 16-core
  allocation does the same stage in ~30 min.
- Therefore **run via `sbatch submit_hervk_driver.sh`, not `./run_hervk_test.sh`**.
  The handout's foreground invocation assumes the cluster profile; with the
  local executor it would run every task on the login node, which has no
  `singularity`.
- `nextflow.config`: `executor.$local.cpus/memory` pinned to 16 / 80 GB
  (Nextflow otherwise sizes the local pool from the node's 94 cores / 470 G, not
  the cgroup). `env.HERVK_REF_RM_OUT` is **commented out** — RERUN_2.md requires
  clean masking, and it was an `env.*` entry, so the handout's
  `unset HERVK_REF_RM_OUT` would not have reached it.
- `preflight.sh`: patched so `-profile standard` is treated as containerised
  (`singularity.enabled` is set outside the profiles block, and the standard
  profile does set `process.container`), and the login-node singularity check is
  a warning rather than a failure. **`bootstrap.sh` reverts this** — it copies
  the handout over local files, protecting only `INPUTS.env`. Re-apply after any
  bootstrap or preflight will fail on `RepeatMasker not on PATH`.

Backups: `*.bak.cluster`, `main.nf.bak.nocap`, `hervk_arch.py.bak.quadratic`.

## History

| job | outcome |
|---|---|
| 23613713 | TIMEOUT 12 h — cluster profile, ~85 min queue wait per task |
| 23667884 | killed — `tile_hits` O(n_hits*span) blowup (now upstream in `0681c33`) |
| 23668829 | failed 3 s — Nextflow refuses `-latest` on a dirty asset repo |
| 23668844 | TIMEOUT 4 h — inside RepeatMasker ProcessRepeats |
| 23670310 | ProcessRepeats rerun standalone, 2h14m, COMPLETED |
| 23671057 | **run 1 PASS** (`28109e7`) |
| 23684967 | cancelled 2h02m — fix 4 gap, masking 26.4 Mb again |
| 23687175 | **run 2 PASS** (`a507cc8`), 34m13s |

Run 1's two local patches (`2c6d36d`, `28109e7`) are now redundant — both are in
upstream `0681c33`. They survive as tag `hervk-local-patches-20260824` and
`hervk_resume/graffite-local-patches.bundle`.

## Leftovers, safe to delete

`hervk_resume/` (64.5 MB — run 1's `ref_windows.fa.cat.gz` and the `.out`
derived from it) and `work/hervk_resume/`. Only needed if re-enabling
`HERVK_REF_RM_OUT`, which RERUN_2.md says not to do.
