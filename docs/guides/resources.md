---
title: Resources and scaling
description: >-
  Execution profiles, the --cores shortcut, per-process CPU, memory and time
  allocation, and advice for large repetitive genomes.
---

# Resources and scaling

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `1a050f9`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](../getting-started/v1.0-vs-v1.1.md).

## Execution profiles

Pass one with `-profile` (single dash). All three run the same container.
<span class="src">`nextflow.config:7-24`</span>

| Profile | Executor | Notes |
|---|---|---|
| `standard` | `local` | The default. Every task runs on the machine you launched from. |
| `cluster` | `slurm` | One SLURM job per task. Tasks run in `$SLURM_TMPDIR` and copy results back. |
| `cloud` | `aws` | Does not launch: Nextflow 26.04.6 stops with `Unknown executor name: aws`. The AWS Batch executor is called `awsbatch`, and needs your own AWS settings in a config passed with `-c`. |

```bash
nextflow run cgroza/GraffiTE -r v1.1dev -profile cluster ...
```

Site-specific settings (a partition, an account, a container cache) go in a config file of your
own, passed with `-c site.config`. The `cluster` profile sets the executor, the container and
`process.scratch`, and nothing else. It sets no queue and no retry policy, and most SLURM sites
need at least a queue, plus an account where the site bills one.

```groovy title="site.config"
process {
    queue          = 'nocona'                  // your partition
    clusterOptions = '--nodes=1 --ntasks=1'    // add --account=... if your site bills one
    errorStrategy  = { task.exitStatus in [137, 140, 143] ? 'retry' : 'terminate' }
    maxRetries     = 3
}
```

```bash
nextflow run cgroza/GraffiTE -r v1.1dev -profile cluster -c site.config ...
```

137 and 143 are SIGKILL and SIGTERM, which is how SLURM ends a job that ran past its time or
memory limit, and some sites report the time limit as 140. Retrying anything else spends an
allocation on a failure that will repeat.

Two defaults to check against your site:

- `process.scratch = '$SLURM_TMPDIR'` <span class="src">`nextflow.config:16`</span> assumes the
  scheduler defines that variable, and many sites do not. Where it is unset, set `process.scratch`
  to a path of your own or to `false`.
- Several `*_memory` parameters default to `null`, so the pipeline requests no `--mem`. On a
  scheduler that needs one, set them; see [Per-process allocation](#per-process-allocation).

If the container cannot write to `/tmp` on your nodes, point `--container_tmp` at a directory it
can; see [Paths and `/tmp`](../getting-started/installation.md#paths-and-tmp).

See also the
[Nextflow executor documentation](https://www.nextflow.io/docs/latest/executor.html).

## The `--cores` shortcut

`--cores N` sets `cpus` to `N` for every process whose CPU count is configurable, in place of all
the `*_threads` parameters, and replaces the `32` that `pav_asm` asks for. Set it for a quick run
on one machine; set the individual parameters for a cluster, where a lone `bcftools` step does
not need the 40 CPUs an alignment does. <span class="src">`nextflow.config:50,176-326`</span>

Processes pinned to one CPU ignore it: `break_scaffold`, `tsd_prep`, `tsd_search`, `tsd_report`,
`merge_VCFs`, `hervk_reconcile`, and the eight methylation processes.
<span class="src">`nextflow.config:179-181,227-241,273-277,293-337`</span>

## Per-process allocation

Each `withName` block in `nextflow.config` reads a `*_threads`, `*_memory` and `*_time`
parameter. Memory and time are Nextflow strings such as `"40G"` and `"12h"`. A memory default of
`null` means the process declares no memory requirement at all; on a scheduler that needs one, set
it. <span class="src">`nextflow.config:123-161,171-327`</span>

| Process | CPUs | Memory | Time | Source |
|---|---|---|---|---|
| `break_scaffold` | `1` | none | none | <span class="src">`nextflow.config:179-181`</span> |
| `map_asm` | `--map_asm_threads` `1` | `--map_asm_memory` `null` | `--map_asm_time` `"3h"` | <span class="src">`nextflow.config:135-137,182-186`</span> |
| `map_longreads` | `--map_longreads_threads` `1` | `--map_longreads_memory` `null` | `--map_longreads_time` `"12h"` | <span class="src">`nextflow.config:138-140,187-191`</span> |
| `sniffles_sample_call`, `sniffles_population_call` | `--sniffles_threads` `1` | `--sniffles_memory` `null` | `--sniffles_time` `"12h"` | <span class="src">`nextflow.config:149-151,192-201`</span> |
| `svim_asm`, `truvari_merge` | `--svim_asm_threads` `1` | `--svim_asm_memory` `null` | `--svim_asm_time` `"12h"` | <span class="src">`nextflow.config:154-156,202-211`</span> |
| `pav_asm` | `32` | `--pav_memory` `"120G"` | `--pav_time` `"12h"` | <span class="src">`nextflow.config:152-153,338-343`</span> |
| `split_repeatmask`, `repeatmask_VCF`, `concat_repeatmask` | `--repeatmasker_threads` `1` | `--repeatmasker_memory` `"10G"` | `--repeatmasker_time` `"12h"` | <span class="src">`nextflow.config:146-148,212-226`</span> |
| `tsd_prep`, `tsd_search`, `tsd_report` | `1` | `--tsd_memory` `"10G"` | `--tsd_time` `"1h"` | <span class="src">`nextflow.config:157-158,227-241`</span> |
| `hervk_annotate` | `--hervk_annotate_threads` `1` | `--hervk_annotate_memory` `"10G"` | `--hervk_annotate_time` `"12h"` | <span class="src">`nextflow.config:127-129,288-292`</span> |
| `hervk_reconcile` | `1` | `--hervk_reconcile_memory` `"10G"` | `--hervk_reconcile_time` `"1h"` | <span class="src">`nextflow.config:130-131,293-297`</span> |
| `pangenie_index`, `pangenie` | `--pangenie_threads` `1` | `--pangenie_memory` `null` | `--pangenie_time` `"12h"` | <span class="src">`nextflow.config:143-145,242-251`</span> |
| `make_graph` | `--make_graph_threads` `1` | `--make_graph_memory` `"40G"` | `--make_graph_time` `"6h"` | <span class="src">`nextflow.config:132-134,252-256`</span> |
| `bam_to_fastq`, `graph_align_reads` | `--graph_align_threads` `1` | `--graph_align_memory` `null` | `--graph_align_time` `"12h"` | <span class="src">`nextflow.config:124-126,257-267`</span> |
| `vg_call` | `--vg_call_threads` `1` | `--vg_call_memory` `null` | `--vg_call_time` `"2h"` | <span class="src">`nextflow.config:159-161,268-272`</span> |
| `merge_VCFs` | `1` | `--merge_vcf_memory` `"10G"` | `--merge_vcf_time` `"1h"` | <span class="src">`nextflow.config:141-142,273-277`</span> |
| `bamtags_to_BED` | `2` | `50 GB` | `6 h` | <span class="src">`nextflow.config:298-302`</span> |
| `lift_epigenome`, `merge_CSV` | `1` | `60 GB` | `6 h` | <span class="src">`nextflow.config:303-312`</span> |
| `index_graph`, `annotate_VCF`, `annotate_BED`, `BED_to_graph`, `merge_BED` | `1` | `40 GB` | `6 h` | <span class="src">`nextflow.config:313-337`</span> |

Three things the table does not show:

- `graph_align_reads` has `errorStrategy = 'finish'`: when one sample's alignment fails, the
  samples already running finish and the run then stops, instead of being killed at once.
  <span class="src">`nextflow.config:266`</span>
- The default thread count is `1` everywhere. A run with defaults aligns each assembly on one CPU.
  Raise `--map_asm_threads`, `--map_longreads_threads`, `--pangenie_threads` and
  `--graph_align_threads` first; those are the steps that scale with CPUs.
- `repmask_vcf.sh` sizes RepeatMasker and ULTRA from `nproc`, not from the `cpus` you allocate.
  On a shared node that reports every core, request a whole node for `repeatmask_VCF` or accept
  the oversubscription. <span class="src">`bin/repmask_vcf.sh:67,92`</span>

What drives each step, as the v1.0 README put it: discovery scales with genome size, the merge
with the number of assemblies, and genotyping with genome size and the size of the read sets. We
have no measurements for RepeatMasker or the graph steps.

`-with-report` writes an HTML report with the peak CPU and memory of every task. Run once with
generous values, read the report, then trim.

## Large, complex and highly repetitive genomes

The defaults, in particular for memory and wall time, can be too small for large, repeat-rich
genomes such as maize. Nextflow's error message for a job that ran out of memory or time is often
about something else, so a run that fails with a confusing message on a large genome is worth
retrying with more of both before looking for another cause. Maize runs have been reported to
complete with requests of the order of 120 h and 400 GB per process, the long-read alignment being
the most demanding step; those are requested values, not measured peaks.

## Measured resource usage

The measured table from the paper (human, *Cannabis sativa* and *Zea mays*, on v1.0) lives in the
[README on the `main` branch](https://github.com/cgroza/GraffiTE/blob/main/README.md#resource-usage-examples).
It was not re-run for v1.1.
