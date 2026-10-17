---
title: Resources and scaling
description: >-
  Execution profiles, the --cores shortcut, per-process CPU, memory and time
  allocation, and advice for large repetitive genomes.
---

# Resources and scaling

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1.0.0` at commit `bf51777`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](../getting-started/v1.0-vs-v1.1.md).

## Execution profiles

Pass one with `-profile` (single dash). All three run the same image, `cgroza/graffite:1.1.0.0`:
`standard` and `cluster` with Apptainer, `cloud` with Docker.
<span class="src">`nextflow.config:7-29`</span>

| Profile | Executor | Notes |
|---|---|---|
| `standard` | `local` | The default. Every task runs on the machine you launched from. |
| `cluster` | `slurm` | One SLURM job per task. Tasks run in `$SLURM_TMPDIR` and copy results back. |
| `cloud` | `awsbatch` | One AWS Batch job per task, with Docker. Needs your AWS settings in a config passed with `-c`: a Batch queue, a region and an S3 work directory. `pav_asm` cannot run. Not tested. |

```bash
nextflow run cgroza/GraffiTE -r v1.1.0.0 -profile cluster ...
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
nextflow run cgroza/GraffiTE -r v1.1.0.0 -profile cluster -c site.config ...
```

137 and 143 are SIGKILL and SIGTERM, which is how SLURM ends a job that ran past its time or
memory limit, and some sites report the time limit as 140. Retrying anything else spends an
allocation on a failure that will repeat.

Two defaults to check against your site:

- `process.scratch = '$SLURM_TMPDIR'` <span class="src">`nextflow.config:16`</span> names a
  variable that many sites do not define. Where it is unset, set `process.scratch` to a path of
  your own or to `false`.
- Several `*_memory` parameters default to `null`, so the pipeline requests no `--mem`. On a
  scheduler that needs one, set them; see [Per-process allocation](#per-process-allocation).

If the container cannot write to `/tmp` on your nodes, point `--container_tmp` at a directory it
can; see [Paths and `/tmp`](../getting-started/installation.md#paths-and-tmp).

See also the
[Nextflow executor documentation](https://www.nextflow.io/docs/latest/executor.html).

## `--cores` shortcut

`--cores N` sets `cpus` to `N` for every process whose CPU count is configurable, in place of all
the `*_threads` parameters, and replaces `pav_asm`'s default of `32`. Set it for a quick run
on one machine; set the individual parameters for a cluster, where a lone `bcftools` step does
not need the 40 CPUs an alignment does. <span class="src">`nextflow.config:55,184-350`</span>

Processes with a hardcoded `cpus` ignore it: `break_scaffold`, `tsd_prep`, `tsd_search`,
`tsd_report`, `merge_VCFs`, `trusted_genotypes`, `genotyping_audit`, `hervk_reconcile`, and the
eight methylation processes. Each gets one CPU except `bamtags_to_BED`, which gets two.
<span class="src">`nextflow.config:185-187,233-247,279-293,299-343`</span>

## Per-process allocation

The `withName` blocks in `nextflow.config` read the `*_threads`, `*_memory` and `*_time`
parameters. Memory and time are Nextflow strings such as `"40G"` and `"12h"`. A memory default of
`null` means the process declares no memory requirement at all; on a scheduler that needs one, set
it. <span class="src">`nextflow.config:128-166,184-350`</span>

| Process | CPUs | Memory | Time | Source |
|---|---|---|---|---|
| `break_scaffold` | `1` | none | none | <span class="src">`nextflow.config:185-187`</span> |
| `map_asm` | `--map_asm_threads` `1` | `--map_asm_memory` `null` | `--map_asm_time` `"3h"` | <span class="src">`nextflow.config:140-142,188-192`</span> |
| `map_longreads` | `--map_longreads_threads` `1` | `--map_longreads_memory` `null` | `--map_longreads_time` `"12h"` | <span class="src">`nextflow.config:143-145,193-197`</span> |
| `sniffles_sample_call`, `sniffles_population_call` | `--sniffles_threads` `1` | `--sniffles_memory` `null` | `--sniffles_time` `"12h"` | <span class="src">`nextflow.config:154-156,198-207`</span> |
| `svim_asm`, `truvari_merge` | `--svim_asm_threads` `1` | `--svim_asm_memory` `null` | `--svim_asm_time` `"12h"` | <span class="src">`nextflow.config:159-161,208-217`</span> |
| `pav_asm` | `32` | `--pav_memory` `"120G"` | `--pav_time` `"12h"` | <span class="src">`nextflow.config:157-158,344-349`</span> |
| `split_repeatmask`, `repeatmask_VCF`, `concat_repeatmask` | `--repeatmasker_threads` `1` | `--repeatmasker_memory` `"10G"` | `--repeatmasker_time` `"12h"` | <span class="src">`nextflow.config:151-153,218-232`</span> |
| `tsd_prep`, `tsd_search`, `tsd_report` | `1` | `--tsd_memory` `"10G"` | `--tsd_time` `"1h"` | <span class="src">`nextflow.config:162-163,233-247`</span> |
| `hervk_annotate` | `--hervk_annotate_threads` `1` | `--hervk_annotate_memory` `"10G"` | `--hervk_annotate_time` `"12h"` | <span class="src">`nextflow.config:132-134,294-298`</span> |
| `hervk_reconcile` | `1` | `--hervk_reconcile_memory` `"10G"` | `--hervk_reconcile_time` `"1h"` | <span class="src">`nextflow.config:135-136,299-303`</span> |
| `pangenie_index`, `pangenie` | `--pangenie_threads` `1` | `--pangenie_memory` `"32G"` | `--pangenie_time` `"12h"` | <span class="src">`nextflow.config:148-150,248-257`</span> |
| `make_graph` | `--make_graph_threads` `1` | `--make_graph_memory` `"40G"` | `--make_graph_time` `"6h"` | <span class="src">`nextflow.config:137-139,258-262`</span> |
| `bam_to_fastq`, `graph_align_reads` | `--graph_align_threads` `1` | `--graph_align_memory` `null` | `--graph_align_time` `"12h"` | <span class="src">`nextflow.config:129-131,263-273`</span> |
| `vg_call` | `--vg_call_threads` `1` | `--vg_call_memory` `null` | `--vg_call_time` `"2h"` | <span class="src">`nextflow.config:164-166,274-278`</span> |
| `merge_VCFs` | `1` | `--merge_vcf_memory` `"10G"` | `--merge_vcf_time` `"1h"` | <span class="src">`nextflow.config:146-147,279-283`</span> |
| `bamtags_to_BED` | `2` | `50 GB` | `6 h` | <span class="src">`nextflow.config:304-308`</span> |
| `lift_epigenome`, `merge_CSV` | `1` | `60 GB` | `6 h` | <span class="src">`nextflow.config:309-318`</span> |
| `index_graph`, `annotate_VCF`, `annotate_BED`, `BED_to_graph`, `merge_BED` | `1` | `40 GB` | `6 h` | <span class="src">`nextflow.config:319-343`</span> |

Three things the table does not show:

- `graph_align_reads` has `errorStrategy = 'finish'`: when one sample's alignment fails, the
  samples already running finish and the run then stops, instead of being killed at once.
  <span class="src">`nextflow.config:272`</span>
- Every `*_threads` parameter defaults to `1`. A run with defaults aligns each assembly on one CPU.
  Raise `--map_asm_threads`, `--map_longreads_threads`, `--pangenie_threads` and
  `--graph_align_threads` first; those are the steps that scale with CPUs.
- `repmask_vcf.sh` sizes RepeatMasker and ULTRA from `nproc`, not from the `cpus` you allocate.
  On a shared node that reports every core, request a whole node for `repeatmask_VCF` or accept
  the oversubscription. <span class="src">`bin/repmask_vcf.sh:67,92`</span>

As the v1.0 README put it, discovery scales with genome size, the merge with the number of
assemblies, and genotyping with genome size and the size of the read sets. We
have no measurements for RepeatMasker or the graph steps.

`-with-report` writes an HTML report with the peak CPU and memory of every task. Run once with
generous values, read the report, then trim.

## Large, complex and highly repetitive genomes

The defaults, in particular for memory and wall time, can be too small for large, repeat-rich
genomes such as maize. Nextflow's error message for a job that ran out of memory or time is often
about something else, so if a run on a large genome fails with a confusing message, retry it with
more of both before looking for another cause. The v1.0 README reports maize runs that completed
with requests of the order of 120 h and 400 GB per process, the long-read alignment being the most
demanding step; those are requested values, not measured peaks.

## Measured resource usage

The measured table from the paper (human, *Cannabis sativa* and *Zea mays*, on v1.0) lives in the
[README on the `main` branch](https://github.com/cgroza/GraffiTE/blob/main/README.md#resource-usage-examples).
It was not re-run for v1.1.
