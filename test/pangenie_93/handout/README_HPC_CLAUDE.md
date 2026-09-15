# Issue #93 PanGenie test — handout

You are running one test on a cluster and reporting what came back. Everything
needed is in this directory. Read this file top to bottom before starting.

## What is being tested

Issue #93 (https://github.com/cgroza/GraffiTE/issues/93) reported two failures
on the v1.1dev PanGenie path. Branch `fix/pangenie-dup-ids-93` fixes both:

1. **`pangenie_index` crashed** in `merge_vcfs.py` with
   `An ID needs to be provided for each individual ID`. The reporter's
   `pangenome.vcf` has 31 pairs of records with the same CHROM/POS/REF/ALT and
   different IDs; `bcftools norm -m+` joined each pair into one ALT with two
   IDs. The new `bin/pangenie_graph_vcf.py` merges such records into one graph
   variant before `norm`, and writes `4_Genotyping/pangenie_graph_variants.tsv`:
   one row per ALT allele of `pangenome.vcf`, with its graph ID (which PanGenie
   copies to INFO/ID of the genotyped VCFs) and whether it made it into the
   graph.
2. **Only 1 of 6 samples was genotyped.** `pangenie` received the reference as
   a one-item queue channel, so it ran once. It now gets a value channel.

Also on the branch, not exercised by this run's inputs except where noted:

- `truvari_merge` with two or more discovery VCFs read a file that no longer
  existed (not exercised: this run starts from `--graffite_vcf`).
- The workflow did not compile under Nextflow 26's strict syntax.
  **Exercised** if the site's Nextflow is 26.x (strict syntax is its default),
  and `preflight.sh` runs `nextflow lint` when the site's Nextflow has it.

Until now, `PanGenie-index` and `PanGenie` have not run on this branch: the
fixes were tested locally without the container. This run is the first time
the whole PanGenie path executes on the fixed code.

## The one thing not to do

**Do not deduplicate or otherwise edit `pangenome.vcf`.** The reporter worked
around the crash by deduplicating it; the point here is that the pipeline
handles the original file. `preflight.sh` warns if the file has no duplicates.

## Inputs

Everything is driven by **`INPUTS.env`**. Populate the paths at the top:

| variable | what |
|---|---|
| `GRAFFITE_VCF` | the reporter's `pangenome.vcf` (original, 51,910 records) |
| `REFERENCE` | the reference it was called against (EptFus1.0, `mEF_genomic.fa` in the issue log) |
| `READS_1` | one short-read FASTQ (interleaved is fine) |
| `READS_2` | optional second read set; empty reuses `READS_1` as a second sample |
| `GRAFFITE_SIF` | optional local `.sif`; empty lets Nextflow pull the image |
| `OUTDIR`, `PROFILE`, `CPUS` | defaults `pangenie_93_run`, `cluster`, `16` |
| `PANGENIE_MEMORY`, `PANGENIE_TIME` | defaults `120G`, `24h` (the config leaves memory unset) |
| `REVISION` | `fix/pangenie-dup-ids-93` |

All paths must be **absolute** (`--contain` in `nextflow.config`).

## Run it

```bash
mkdir -p /xdisk/cgoubert/cgoubert/GraffiTE1.1/issue93 && cd $_
module load nextflow            # whatever the site provides

# first time only: fetch the handout out of the Nextflow asset cache
nextflow pull cgroza/GraffiTE -r fix/pangenie-dup-ids-93
cp ~/.nextflow/assets/cgroza/GraffiTE/test/pangenie_93/handout/* .

./bootstrap.sh                  # nextflow pull + refresh this handout
$EDITOR INPUTS.env              # fill in the paths
./preflight.sh                  # stop here if it fails
./run_pangenie_test.sh          # pipeline, then assertions
./bundle_results.sh             # -> pangenie_93_results_<date>.tar.gz
```

Run under `tmux`/`screen` or as a batch job: with `-profile cluster` Nextflow
submits its own slurm jobs and the driver must stay alive.

Expected work: `pangenie_index` once (graph of ~50,000 variants on a ~2 Gb
genome), then `pangenie` once **per sample** — two tasks. If the trace shows
one `pangenie` task, the second fix did not take; report it.

`run_pangenie_test.sh` passes `-latest` and `-resume`, and drops `-resume` when
the pipeline commit moves (Nextflow does not hash `bin/`, so a resumed
`pangenie_index` could come from the old `pangenie_graph_vcf.py`).

## What must pass

`assert_pangenie_test.py` exits non-zero on any `[FAIL]`:

- the table, the merged VCF and both `<sample>_genotyping.vcf.gz` exist
- one table row per ALT allele of `pangenome.vcf`; 31 duplicate rows for the
  #93 file, each sharing the graph ID of its first occurrence
- no graph ID empty or holding `; , : = |`
- both samples are columns of `GraffiTE.merged.genotypes.vcf.gz`
- the graph IDs in INFO/ID of each per-sample VCF and of the merged VCF are
  exactly the table's `in_graph=yes` graph IDs

## Reading a failure

| symptom | what it means |
|---|---|
| `pangenie_index` fails with `An ID needs to be provided` | the old script ran: check `REVISION`, the cached commit, and that `-resume` was dropped |
| `pangenie_graph_vcf.py: command not found` or permission denied | `bin/pangenie_graph_vcf.py` lost its executable bit; report the `ls -l` of it in the asset dir |
| only one `<sample>_genotyping.vcf.gz` / one `pangenie` task | the value-channel fix did not take |
| PanGenie-index killed / exit 137 / 140 | memory or time; raise `PANGENIE_MEMORY` / `PANGENIE_TIME` and rerun (resume is fine) |
| INFO/ID set differs from the table | report the counts in that line and the first few IDs on each side; do not edit anything |
| compile error mentioning syntax | report `nextflow -v` and the error verbatim; this is the Nextflow 26 fix |

## What to report back

Send the bundle plus a short summary:

- the assertion log verbatim, including every `NOTE:` line
- `nextflow -v`, and whether it ran with the strict parser (26.x default) or not
- the `pangenie_graph_vcf.py:` lines from the `pangenie_index` `.command.err`
  in the bundle (records, duplicates merged, IDs replaced, alleles not in graph)
- the `in_graph=no` count against the predicted 1,569 (`NOTE` line). This is a
  prediction from a stand-in reference, not a test; report it either way
- wall time and peak memory of `pangenie_index` and each `pangenie` task
  (from `nextflow_trace.txt` / the report)

Do not tune parameters to make assertions pass. If something disagrees, that is
the result.
