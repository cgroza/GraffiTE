# Brief for the SSH agent: issue #127 mtime test on Puma

You are on the University of Arizona HPC (Puma). Your job is to run one test,
make sure it ran properly, and bring back the results. You do not decide what
the results mean: the person who sent you will. Report facts, quote files, and
separate what you observed from anything you infer.

Full reference: `README_HPC_CLAUDE.md` in this directory. This brief is the
order of work and the things that are likely to go wrong on a first run.

## Context in five lines

- Issue: https://github.com/cgroza/GraffiTE/issues/127. `graph_align_reads`
  fails with `index/index.dist is newer than index/index.longread.withzip.min`
  under `-profile cluster` with `--graph_method giraffe`.
- vg giraffe checks index file modification times. Under the cluster profile,
  `make_graph` runs in scratch and Nextflow copies `index/` back with
  `cp -fRL`, which resets them. GNU `cp` copies in inode order.
- The question is whether this reproduces here, or depends on the reporter's
  filesystem. It decides whether GraffiTE gets patched or the reporter gets a
  workaround.
- Branch `test/giraffe-mtime-127` = `main` (v1.1.0.0) + this handout. No
  pipeline code differs from the release.
- Nothing in this handout has run yet. Its tally logic was dry-run on fake
  data only. Expect first-run bugs in the scripts.

## Order of work

### 1. Set up (login node)

```bash
mkdir -p /xdisk/cgoubert/cgoubert/GraffiTE1.1/issue127 && cd $_
module load nextflow
nextflow pull cgroza/GraffiTE -r test/giraffe-mtime-127
cp ~/.nextflow/assets/cgroza/GraffiTE/test/giraffe_mtime_127/handout/* .
./bootstrap.sh
```

Check `bootstrap.sh` printed a revision on `test/giraffe-mtime-127`.

### 2. Find the inputs

This test reuses the synthetic test set. Its `WORKDIR` was most likely
`/xdisk/cgoubert/cgoubert/GraffiTE1.1/synthetic`; confirm by looking for its
`INPUTS.env` and `build/`. You need, under it:

| file | used by |
|---|---|
| `build/ref/synth.fa` | both parts |
| `build/reads_long.csv` and the `build/reads/L1.long.fq.gz` it names | Part B, Part A |
| `runs/spine/3_TSD_search/pangenome.vcf` | both parts |

Check that the paths inside `reads_long.csv` exist: the CSV holds absolute
paths written at build time, so a moved or rebuilt set breaks it silently.

If `build/` is missing, run that set's `./preflight.sh` (it builds inputs and
runs nothing else). If the spine's `pangenome.vcf` is missing, look for any
other GraffiTE `pangenome.vcf` called against `synth.fa` (other cells under
`runs/`), and set `GRAFFITE_VCF` to it. If there is none, stop and report: do
not run the spine just for this.

Fill `INPUTS.env`: `WORKDIR` (this directory is fine), `SYNTH_DIR`, and
`CONTAINER_TMP` (a directory on `/xdisk`, as the synthetic set uses). Leave
the rest at their defaults.

### 3. Find the sbatch account and partition

`submit.sh` carries no `--account` or `--partition`. Use what the synthetic
driver jobs used: look at earlier `synth_driver_*.out` files, or
`sacct -u $USER -S 2026-09-01 --format=JobID,JobName,Account,Partition | grep graffite`.
Pass them on the command line (`sbatch --account=... --partition=... submit.sh`);
do not edit them into `submit.sh`.

### 4. Smoke test in a short interactive allocation

The login node has no `apptainer`. Before the full job, take a short
interactive allocation (8 cores, 32 GB, 1 h) and run one Part A build:

```bash
REPS_A=1 ./part_a.sh
```

`REPS_A` on the command line overrides `INPUTS.env`, so nothing needs
editing. `submit.sh` later runs with the `INPUTS.env` value (10).

What should come out: the image pulls once into `$WORKDIR/singularity`
(several GB, minutes), the index builds, and five lines per source directory:
`control`, `cp`, `forced_bad`, `forced_good`, `mv`. `forced_bad` must say
`observed=fail_mtime`. `control` and `forced_good` must say `observed=pass`.
If they do, run `./summarize.sh` and check `VALIDITY: OK`.

Also note from `partA/meta.txt` whether `SLURM_TMPDIR` is set on Puma. If it
is unset, the Nextflow scratch falls back to `/tmp`, which is what the test
then measures; record it either way.

### 5. Full run

```bash
sbatch --account=... --partition=... submit.sh
```

Parts A and B, the tally, the bundle. The synthetic genome is ~200 kb, so
expect well under the 6 h wall; if it gets close, report the stage it was in.

### 6. Check Part B actually used scratch

For one `scratch_<n>` cell and for `scratch_move`, open
`partB/<cell>/make_graph_unstage.txt` and `make_graph_scratch_line.txt`.
The scratch cells must show `nxf_fs_copy` / `cp -fRL` and an `NXF_SCRATCH=`
line; `scratch_move` must show `nxf_fs_move` / `mv -f`. If a scratch cell has
no `NXF_SCRATCH=` line, that cell tested nothing; say so in the report.

## What you may and may not change

You **may** fix bugs in the handout scripts (`part_a.sh`, `part_b.sh`,
`summarize.sh`, `bundle_results.sh`, `common.sh`, `submit.sh`) so that they do
what `README_HPC_CLAUDE.md` says they do: a wrong path, a bind mount, a
quoting error, a GNU-vs-other tool flag. Keep every change as a diff
(`diff -u original fixed`) and put the diffs in the report. Do not push them.

You **may not**:

- change the pipeline, its parameters, or the GraffiTE revision
- change the variants in Part A, the cells in Part B, or the label rules in
  `summarize.sh`
- rerun a cell to turn a failure into a pass. A `fail_mtime` in a `cp` row or
  a `scratch_<n>` cell is a result
- reuse work directories or add `-resume`

If the VALIDITY line says `INVALID`, stop after the bundle and report: the
rest of the run does not count.

## What to bring back

1. The bundle path: `issue127_results_<date>.tar.gz` in this directory.
2. `VERDICT.txt`, verbatim.
3. From `partA/meta.txt`: vg version, `cp` and `mv` versions, kernel,
   filesystem type of the work directory, `SLURM_TMPDIR`.
4. For each Part A source directory: its filesystem type, and in how many of
   the builds `index.dist` had a lower inode number than both
   `index.longread.withzip.min` and `index.longread.zipcodes` (from
   `partA/inodes.tsv`, stage `src`).
5. For Part B: the unstage block of one scratch cell and of `scratch_move`,
   and the per-cell index mtimes from `partB/index_mtimes.tsv`.
6. Any failure outside `make_graph` and `graph_align_reads`, quoted from its
   `.command.err`.
7. Every change you made to the scripts, as diffs, with one line each on why.
8. Wall time of the job (`sacct -j <id> --format=Elapsed,MaxRSS`).

Mark each statement as observed (with the file it came from) or inferred. Do
not explain why the result came out as it did, and do not recommend a fix.
