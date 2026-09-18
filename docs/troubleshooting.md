---
title: Troubleshooting
description: >-
  Messages GraffiTE prints on purpose, failures users have reported, and what to
  do about each.
---

# Troubleshooting

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `1a050f9`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](getting-started/v1.0-vs-v1.1.md).

---

## Nextflow 26 rejects main.nf with "Unexpected input"

```
Error main.nf:158:27: Unexpected input: 'splitText'
```

You are running a checkout older than commit `4605630` (2026-09-14) on Nextflow 26, whose
default parser is the strict one. Update to a newer commit, or select the legacy parser for that
checkout:

```bash
export NXF_SYNTAX_PARSER=v1
```

One more thing to know about Nextflow 26: its parser passes a command-line `--flag false` as the
string `"false"`, which Groovy treats as true. Before commit `cfaff1e` that made `--genotype false`
run genotyping anyway on Nextflow 26. Since that commit every boolean parameter is read through
`isOn()` <span class="src">`module/main.nf:1-9`</span>, and both parsers agree.

---

## Messages the workflow prints on purpose

Each of these stops the run before any process starts. The fix is in the message.

| Message | Cause | What to do |
|---|---|---|
| `No input given. Pass one of --longreads, --bams, --assemblies, --pav, --svs (discovery), --vcf (a merged SV VCF), --RM_dir (RepeatMasker output of an earlier run) or --graffite_vcf (a pangenome.vcf from an earlier run).` <span class="src">`main.nf:143`</span> | No entry flag on the command line. | Pick one on [Choosing your inputs](getting-started/choosing-your-inputs.md). |
| `--vcf cannot be combined with --assemblies. Pass --vcf alone, or drop it and use --svs to add your own per-sample VCFs to the discovery merge.` <span class="src">`main.nf:46-49`</span> | `--vcf` replaces the discovery merge, so a discovery flag beside it has nothing to feed. | Use `--svs` if you want your VCF merged with caller output. |
| `--graph_method precomputed builds nothing itself: it needs --graph (...) and either --vcfs (...) or --graph_alignments (...).` <span class="src">`main.nf:53-56`</span> | `precomputed` reuses a graph and alignments from an earlier run and was given neither. | Pass `--graph` plus `--vcfs` or `--graph_alignments`, or choose another method. |
| `Unsupported --graph_method. --graph_method must be pangenie, giraffe, graphaligner or precomputed.` <span class="src">`main.nf:255`</span> | A misspelt method. | One of the four. |
| `--graffite_vcf skips discovery, and HERV-K reconciliation needs the hervk_annotate outputs discovery produces. Pass --hervk_reconcile false, or start from --RM_dir instead of --graffite_vcf.` <span class="src">`main.nf:61-63`</span> | `--human` consolidation reads files that only the annotation stage writes. | Either option in the message. |
| `Error main.nf:12:1: Invalid include source: '<dir>/panmethyl/module.nf'` <span class="src">`main.nf:12`</span> | A clone made without `--recurse-submodules`; the `panmethyl` submodule directory is empty. This one is Nextflow's own message. | `git submodule update --init` in the clone. `nextflow pull` clones do not hit this. |

---

## A process cannot see an input file

Every container runs with `--contain --bind $(pwd):/tmp` <span class="src">`nextflow.config:5`</span>.
Files outside the launch directory are only visible if given by absolute path, and `/tmp`
inside the container is the launch directory. Symptoms are `No such file` inside a process for
a file that exists on the host, or a full filesystem when the launch directory is small. Use
absolute paths, and bind a larger scratch directory to `/tmp` if needed; see
[Installation](getting-started/installation.md).

---

## `TSD_summary.txt` is empty and no record has `INFO/TSD`

On `v1.1dev` before commit `884afa8` the flank extraction failed without an error on a
gzip-compressed reference, and a missing `exact_match.py` read as no hit on every variant. Both
now stop the run instead. If your copy is older than that commit, update it
([Installation](getting-started/installation.md#updating)); `test/tsd/test_tsd_chain.sh`
checks the chain.

---

## `variant IDs must be no greater than 50 characters`

The annotation script stops on any variant ID longer than 50 characters
<span class="src">`bin/repmask_vcf.sh:41-45`</span>, the length RepeatMasker keeps intact in
its output. IDs written by the pipeline's own callers are short. The message appears with
`--vcf`, whose IDs pass through unchanged <span class="src">`module/main.nf:196-198`</span>,
or with `--svs`: rename them before the run.

---

## RepeatMasker runs out of memory or time

Stage B is split by contig <span class="src">`module/main.nf:472-485`</span> and each piece runs
with `--repeatmasker_memory` (default `10G`) and `--repeatmasker_time` (default `12h`)
<span class="src">`nextflow.config:146-148`</span>. Raise them, and raise
`--repeatmasker_threads`, for large contigs or large libraries. Nextflow's message when a job is
killed by the scheduler can be misleading; check the `.command.log` in the failing work
directory. For maize-sized, repeat-rich genomes users have needed up to 120 h and 400 GB per
process, with long-read alignment the most demanding step. See
[Resources and scaling](guides/resources.md).

---

## The run finishes but there is no `3_TSD_search/pangenome.vcf`

Check whether `tsd_search` and `tsd_report` ran. When the container cannot write to `/tmp` the
pipeline skips them and still reports success. `/tmp` inside the container is the launch
directory by default; point `--container_tmp` at a writable directory instead, see
[Paths and `/tmp`](getting-started/installation.md#paths-and-tmp).

The other way to get here is a run where no chunk survived annotation. Stage B is split one piece
per contig <span class="src">`module/main.nf:472-485`</span>, and a chunk whose records all fail
`--repeat_span_cutoff` drops out by design. If every chunk drops out, `concat_repeatmask` never
runs and nothing says so. `2_Repeat_Filtering/*/genotypes_repmasked_filtered.vcf` holds the
per-chunk record counts.

---

## `cgroza/GraffiTE contains uncommitted changes -- cannot pull from repository`

The pipeline cached under `~/.nextflow/assets/cgroza/GraffiTE/` was edited in place. Delete it
and pull again; see [Installation](getting-started/installation.md#updating).

Usually the edit was to `nextflow.config`, to bind a scratch directory to `/tmp`. That is what
`--container_tmp` is for; use it and leave the cached copy alone.

---

## A full-length LTR element shows `n_hits` greater than 1

RepeatMasker annotates the LTR and the internal region of an LTR retrotransposon as separate
hits when they are separate entries in the library, and the pipeline reports what RepeatMasker
finds. Such an element has `n_hits=2` or `3` and is left out of the subsets that require a
single hit (`pangenome.trusted.vcf`, and `pangenome.human.vcf` apart from the HERV-K pair rule
<span class="src">`module/main.nf:492,514`</span>). It stays in `pangenome.vcf` with all its
hits listed. A library entry that holds the full element as one consensus gives one hit.

---

## The run is slow

Two steps dominate: `samtools sort` after minimap2, whose parallelism is set by `--stSort_t`
and `--stSort_m`, and RepeatMasker with a large library and many SVs. `-with-report report.html`
shows where the time went and which requests to change.

---

## Getting help

Open an issue at <https://github.com/cgroza/GraffiTE/issues>. It is the only channel the
authors watch, and someone may have hit the same thing. Include the commit you run
(`nextflow info cgroza/GraffiTE`), the command line, and the `.command.log` of the failing task.
