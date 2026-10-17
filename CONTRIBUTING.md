# Contributing to GraffiTE

This page is for people who change the code. The user documentation is at
[cgroza.github.io/GraffiTE](https://cgroza.github.io/GraffiTE/).

## Branches

- `main` holds releases, and each release is a tag on it.
- `dev` collects the work for the next release. Open pull requests against `dev`.
- A fix to a released version can branch from `main`, be tagged there, and be merged back into
  `dev`.

## Version numbers

A GraffiTE version has four numbers, W.X.Y.Z. The first release under this scheme is 1.1.0.0.

| Part | Goes up when | Examples |
|---|---|---|
| W | GraffiTE is redesigned: how it runs or what it produces changes wholesale, usually with a new paper. No compatibility promise holds across W versions. | GraffiTE 2 |
| X | By design, a call an earlier version made changes on the same inputs, or a mainstay behaves differently: the trusted or `--human` subsets, a default, a parameter renamed or removed. Adding something new is not an X. | a new rule in the trusted subset; a removed parameter |
| Y | Anything else that changes the pipeline: a bug fix, or an addition that leaves earlier calls as they were (new calls, a new parameter, a new output). When a fix changes results, the changelog says which records and conditions. | BC200 calls added to `--human`, planned for 1.1.2.0; `--human` given as a string (#121) |
| Z | The pipeline logic does not change: a rebuilt container image, a pinned tool or submodule, the test suite, CI. | a rebuilt image with a newer bcftools; the panmethyl pin (#123) |

A change to the documentation alone does not change the version. The site redeploys from `main`.

Every release has:

- a git tag, `vW.X.Y.Z`, on `main`;
- the same tag on the container image, `cgroza/graffite:W.X.Y.Z`, named in `nextflow.config`, so
  that `nextflow run cgroza/GraffiTE -r vW.X.Y.Z` runs the image built for it, while `latest`
  follows the newest release;
- the version in `version.txt` and `manifest.version`; `version.txt` is what the VCF headers carry
  as `##GraffiTE_version`;
- a GitHub release whose notes come from the changelog;
- documentation stamped with the release.

## Synthetic test suite

The scripts in `test/synthetic/handout/` build a small synthetic genome with planted transposable
element insertions, run GraffiTE on it in 20 configurations (the matrix), and record what each run
produced in `OBSERVED.tsv`. `expectations.py check OBSERVED.tsv EXPECTED.tsv` compares a run
with the expectations frozen at the last release:

- exit 0: everything matches;
- exit 1: a structural result moved, such as a planted site that stopped reaching a file it
  reached before, or a subset that lost a member;
- exit 2: a tool heuristic moved: a score, a coordinate, a TSD sequence, an annotation string.

`README_HPC_CLAUDE.md` in that directory says how to build the test set and run the matrix on a
SLURM cluster, and `EXPECTED.tsv` beside it holds the frozen expectations.

## What a release must pass

- Z: exit 0, or exit 2 with the tool versions that moved named in the release notes. Never exit 1.
- Y: structural results move only for the records the change targets, and the notes list them.
- X: a full matrix run, `EXPECTED.tsv` re-frozen and committed with the release, and a run on real
  data.
- W: all of the above, and a "vX vs vW" page in the documentation.

## Release checklist

1. On the release candidate: the full matrix green, and for X and W `EXPECTED.tsv` re-frozen and
   committed.
2. The image built from `Dockerfile` and pushed with the version tag; the tool versions read back
   into `docs/reference/container.md`; the matrix run again on that image.
3. `version.txt` and `manifest.version` set, and `nextflow.config` naming the image tag.
4. Docs: stamps, `-r` commands and blob links pointing at the tag, and the changelog entry dated.
5. Merge into `main`, tag, move `latest`, publish the GitHub release, and check the live site.
6. From a clean account, `nextflow run cgroza/GraffiTE -r vW.X.Y.Z` on the quick install test.
7. Close the issues the release fixes, with a note naming the version.
