# Issue #127 giraffe index mtime test — handout

You are running one test on a cluster and reporting what came back. Everything
needed is in this directory. Read this file top to bottom before starting.

## What is being tested

Issue #127 (https://github.com/cgroza/GraffiTE/issues/127): with
`--graph_method giraffe` and `-profile cluster`, `graph_align_reads` stops at
once with

    error[vg giraffe] index/index.dist is newer than index/index.longread.withzip.min which depends on it

The hypothesis, with what is verified and what is not:

- **Verified (source).** vg giraffe compares file modification times:
  `IndexRegistry::predates()` takes `std::filesystem::last_write_time` of each
  index and requires the distance index to be no newer than the minimizer and
  zipcode indexes (`src/index_registry.cpp`, `src/subcommand/giraffe_main.cpp`
  in vg v1.77.0, the version the container installs).
- **Verified (source).** The cluster profile sets
  `process.scratch = '$SLURM_TMPDIR'`. `make_graph` then builds `index/` in
  scratch and Nextflow copies it back with `cp -fRL` (`nxf_fs_copy` in
  `.command.run`). Plain `cp` does not keep modification times: every file gets
  the time it was written in the work directory.
- **Verified (source).** GNU `cp -R` copies a directory's entries in inode
  number order (gnulib `savedir`, `SAVEDIR_SORT_FASTREAD`). GNU `mv` keeps
  modification times, also across filesystems (`mv.c`,
  `preserve_timestamps = true`).
- **Not verified.** That this is what happened to the reporter. If the scratch
  filesystem gives `index.dist` a higher inode number than the long-read
  minimizer or zipcode file, `cp` writes it later, vg sees it as newer, and the
  run fails. Whether that happens depends on how the filesystem hands out inode
  numbers, so it may happen on one cluster and never on another.

The test answers two questions: does the failure reproduce here, under
GraffiTE's own cluster scratch setting; and does `stageOutMode 'move'` on
`make_graph` avoid it. The pipeline code on this branch is v1.1.0.0
(`main`) plus this handout; nothing in `main.nf`, `module/` or `bin/` is
changed.

## The one thing not to do

**Do not change the pipeline or its parameters to get a pass.** A failure in
the `cp` or `scratch_<n>` rows is a result, not a problem with the run. The only
pass/fail that says the test itself is broken is the VALIDITY line (below).

## Inputs

This test reuses the synthetic test set (`test/synthetic/handout`): its
reference, its long reads, and the spine's `pangenome.vcf`. Point `SYNTH_DIR`
at that set's `WORKDIR`. If the spine was never run there, set `GRAFFITE_VCF`
to any GraffiTE `pangenome.vcf` called against the same reference; if there is
no `build/` either, run the synthetic handout's `./preflight.sh` once, which
builds the inputs and nothing else.

| variable | what |
|---|---|
| `WORKDIR` | absolute, on `/xdisk`. Parts A and B write here. |
| `SYNTH_DIR` | the synthetic set's `WORKDIR` |
| `REFERENCE`, `GRAFFITE_VCF`, `READS_CSV`, `READS_FASTQ` | default to files under `SYNTH_DIR`; override only if they are elsewhere |
| `GRAFFITE_SIF` | leave empty: Part A pulls `docker://cgroza/graffite:1.1.0.0` into `$WORKDIR/singularity` and Part B reuses it |
| `CONTAINER_TMP` | as in the synthetic handout |
| `SCRATCH_DIRS` | Part A source directories. Empty means `${SLURM_TMPDIR:-/tmp}` (what Nextflow uses) plus `$WORKDIR/partA_src` (a second filesystem) |
| `REPS_A`, `REPS_B` | index builds per source in A (10), scratch repeats in B (3) |
| `NF_SCRATCH` | Part B's `process.scratch`; default is the cluster profile's literal `$SLURM_TMPDIR` |
| `CPUS`, `MEM_GB` | must match `--cpus-per-task` and `--mem` in `submit.sh` |

## Run it

```bash
mkdir -p /xdisk/cgoubert/cgoubert/GraffiTE1.1/issue127 && cd $_
module load nextflow

nextflow pull cgroza/GraffiTE -r test/giraffe-mtime-127
cp ~/.nextflow/assets/cgroza/GraffiTE/test/giraffe_mtime_127/handout/* .

./bootstrap.sh                  # pull + refresh this handout
$EDITOR INPUTS.env              # WORKDIR, SYNTH_DIR
sbatch submit.sh                # Part A, Part B, tally, bundle
```

`sbatch submit.sh A` or `B` runs one part. Inside an allocation you can call
`./part_a.sh`, `./part_b.sh [cell ...]`, `./summarize.sh`, `./bundle_results.sh`
directly. The login node has no `apptainer`; run everything inside the
allocation.

Part B uses `-profile standard` (local executor) in one allocation, as the
synthetic suite does on Puma, plus a config overlay that sets
`process.scratch` exactly as the cluster profile does. The scratch staging
code in `.command.run` is written by the same wrapper for every executor, so
this exercises the copy step `-profile cluster` runs without paying Puma's
per-task queue latency. `make_graph_unstage.txt` in each Part B cell shows the
`nxf_unstage_outputs` block that ran; check it says `cp -fRL` for the scratch
cells and `mv -f` for `scratch_move`.

## What each part does

**Part A** (`part_a.sh`, vg only, no Nextflow). For each source directory and
each of `REPS_A` builds: create a task-style scratch directory with
`mktemp -d -t nxf.XXXXXXXXXX`, build the index there with `make_graph`'s recipe
copied verbatim, record every index file's inode number, mtime and readdir
position, then run `vg giraffe --parameter-preset hifi` (the reporter's preset)
against five variants:

| variant | what | expected |
|---|---|---|
| `control` | in place, no copy | pass |
| `cp` | `cp -fRL index <dest>/`, as `nxf_fs_copy` does | predicted from inode order |
| `forced_bad` | every file, then `index.dist` last after `sleep 1` | `fail_mtime` (positive control) |
| `forced_good` | `index.dist` first, `sleep 1`, then the rest | pass |
| `mv` | `mv -f index <dest>/` | pass |

**Part B** (`part_b.sh`, the pipeline). `--graffite_vcf` + `--graph_method
giraffe` + long reads, each cell in its own work directory without `-resume`:

| cell | overlay | expected |
|---|---|---|
| `noscratch` | none | pass |
| `scratch_1..REPS_B` | `scratch.config` | unknown: this is the test |
| `scratch_move` | `scratch_move.config` | pass |

Only `make_graph` and `graph_align_reads` are tallied. If a later process
(`vg_call`, `merge_VCFs`, ...) fails, record it in the summary, but it does not
enter the labels.

## Labels

`summarize.sh` writes `$WORKDIR/VERDICT.txt` with these lines. They were fixed
before the run; do not redefine them afterwards.

- `VALIDITY`: `OK`, or `INVALID` with the reason. Invalid if any `forced_bad`
  does not fail with `is newer than`, any `control`/`forced_good` fails, or
  `noscratch` fails. If invalid, stop and report: the rest means nothing.
- `REPRODUCTION`:
  - `REPRODUCED`: at least one `fail_mtime` in Part A `cp` from the Nextflow
    scratch source, or in a Part B `scratch_<n>` `graph_align_reads`
  - `REPRODUCED_OTHER_FS_ONLY`: `fail_mtime` only in Part A `cp` from the extra
    source directory
  - `NOT_REPRODUCED`: none
- `INODE_ORDER_PREDICTS_CP`: in how many Part A `cp` rows the inode-order
  prediction matched the outcome
- `CANDIDATE_FIX`: `MOVE_OK` if every Part A `mv` and every Part B
  `scratch_move` `graph_align_reads` passed, else `MOVE_FAILED`

## Reading a failure

| symptom | what it means |
|---|---|
| `index build failed` in Part A | read `partA/logs/<src>_r<n>_build.log`. Most likely the VCF and reference disagree on contigs, or the container could not see a path (all paths absolute, all under bound directories) |
| every giraffe run `fail_other` | read one `partA/logs/*.err`; usually the reads file or a bind path |
| Part B `no_trace` | Nextflow did not start; `partB/<cell>/run.log` |
| `scratch_move` unstage block says `cp` | the overlay did not apply; report `nextflow -v` and the block |
| `make_graph_scratch_line.txt` missing for a scratch cell | scratch was not used; the cell does not test anything, report it |

## What to report back

`bundle_results.sh` packs `issue127_results_<date>.tar.gz`. Send it with a short
summary:

- `VERDICT.txt` verbatim
- `partA/meta.txt` (vg, cp, mv versions, kernel, filesystem types,
  `SLURM_TMPDIR`)
- for the Part A source directories: the filesystem type, and whether
  `index.dist` came before or after the two long-read files in inode order in
  each build (from `partA/inodes.tsv`, stage `src`)
- for Part B: the `make_graph_unstage.txt` of one scratch cell and of
  `scratch_move`, and the index mtimes per cell from `partB/index_mtimes.tsv`
- anything that failed outside `make_graph` / `graph_align_reads`, quoted

Report what the files say. Do not explain why a result came out as it did
beyond what the files show.
