# Issue #127: giraffe index mtimes after a scratch copy

`vg giraffe` exits with an error when `index.dist` is newer than the minimizer
or zipcode index built from it. It checks this by comparing file modification
times. Under `-profile cluster`, `make_graph` builds the index in
`$SLURM_TMPDIR`, and Nextflow copies it back with `cp -fRL`, which does not
keep modification times. GNU `cp` copies the files of a directory in inode
order, so `index.dist` ends up newer whenever the scratch filesystem gave it a
higher inode number than the long-read index files. We think this is what
happened in #127, but have not reproduced it.

The scripts in `handout/` test this on a cluster in two parts:

- Part A runs `make_graph`'s recipe and `vg giraffe` without Nextflow. It
  copies the index with `cp` and with `mv`, adds a positive control that
  copies `index.dist` last, and records inode numbers to check that inode
  order predicts the `cp` outcome.
- Part B runs the pipeline on the synthetic test set three ways: without
  scratch, with the cluster profile's scratch setting, and with that setting
  plus `stageOutMode 'move'` on `make_graph`.

The pipeline code is unchanged. `handout/README_HPC_CLAUDE.md` explains how to
run the test and how results are labelled.
