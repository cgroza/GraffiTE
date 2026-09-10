---
title: Parameters
description: >-
  Every GraffiTE parameter, its exact default, what it controls, and where it is
  declared in the source.
---

# Parameters

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `18a76d9`. Every default in this page was read from
    `nextflow.config` or from the code that consumes it; the **Source** column gives the line so
    you can check.

## How to pass parameters

Nextflow distinguishes two kinds of option by the number of leading dashes, and the distinction is
load-bearing:

| Form | Belongs to | Example |
|---|---|---|
| `--name value` | **GraffiTE** — anything on this page | `--assemblies assemblies.csv` |
| `-name value` | **Nextflow itself** | `-resume`, `-profile cluster`, `-with-report` |

```bash
nextflow run cgroza/GraffiTE -r v1.1dev -latest \
  -profile cluster \
  --reference hs37d5.fa \
  --assemblies assemblies.csv \
  --TE_library human_DFAM3.6.fasta \
  --genotype_with reads.csv
```

Defaults live in [`nextflow.config`](https://github.com/cgroza/GraffiTE/blob/v1.1dev/nextflow.config).
You can override them on the command line as above, or edit a local copy of the config, or supply a
`-params-file`.

!!! warning "`-params-file` will not see every parameter"
    Six parameters are read by the workflow but never declared in `nextflow.config` — see
    [Undeclared parameters](#undeclared-parameters). They work on the command line but are invisible
    to schema-based tooling.

!!! note "Underscores and hyphens"
    This page uses the names exactly as declared in `nextflow.config`, which is the form guaranteed
    to work. The current README writes `--genotype-with` (hyphen) for what is declared as
    `genotype_with` (underscore). Prefer the underscore form.

---

## Required inputs

These three are required for essentially every run, and all three have defaults that are *filenames
rather than sensible values* — if you do not set them, GraffiTE looks for `reference.fa` and
`TE_library.fa` in the launch directory and fails if they are absent.

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--reference` | `"reference.fa"` | Reference genome FASTA. Everything is called relative to this. Existence is checked at launch. May be plain or bgzip-compressed; a non-BGZF gzip is transparently re-compressed. | <span class="src">`nextflow.config:41`</span> |
| `--TE_library` | `"TE_library.fa"` | FASTA of repeat consensus sequences passed to RepeatMasker as `-lib`. Required unless you skip Stage B with `--RM_dir` or `--graffite_vcf`. | <span class="src">`nextflow.config:42`</span> |
| `--genotype_with` | `"reads.csv"` | Samplesheet of read sets to genotype. Required when `--genotype true` (the default). See [Samplesheets](samplesheets.md). | <span class="src">`nextflow.config:36`</span> |

---

## Stage A — discovery inputs

At least one of these must be supplied, unless you enter the pipeline further downstream with
`--vcf`, `--RM_dir` or `--graffite_vcf`. They are **additive**: supply several and all of their
variant calls are merged into one non-redundant set.

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--assemblies` | `false` | Samplesheet of genome assemblies. Each is aligned with minimap2 (or winnowmap) and called with `svim-asm haploid`, keeping only `INS`/`DEL` ≥ 100 bp. | <span class="src">`nextflow.config:38`</span> |
| `--longreads` | `false` | Samplesheet of raw long reads. Aligned, then called per-sample and jointly with Sniffles2 at `--minsvlen 100`. | <span class="src">`nextflow.config:37`</span> |
| `--bams` | `false` | Samplesheet of **already aligned** long-read BAMs. Skips the alignment step and goes straight to Sniffles2. Can be combined with `--longreads`. | <span class="src">`nextflow.config:29`</span> |
| `--pav` | `false` | Samplesheet of phased assemblies to call with [PAV](https://github.com/EichlerLab/pav). Runs in its own container. Keeps variants with \|SVLEN\| > 50. | <span class="src">`nextflow.config:39`</span> |
| `--svs` | *(undeclared)* | Samplesheet of per-sample SV VCFs you called yourself. No discovery process runs; the VCFs go straight into the merge. | <span class="src">`main.nf:90`</span> |

See [Stage A — discovery](../guides/discovery.md) for what each backend does and how to choose.

### Discovery tuning

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--aligner` | `"minimap2"` | Aligner for assemblies and long reads. The only other accepted value is `"winnowmap"`. **Any other value produces a process with no script body and the run fails.** | <span class="src">`nextflow.config:70`</span> |
| `--asm_divergence` | `"asm5"` | minimap2 `-x` preset for assembly alignment. Use `asm10` or `asm20` for more divergent assemblies. Applies to assemblies only, not long reads. | <span class="src">`nextflow.config:69`</span> |
| `--break_scaffolds` | `false` | Split input assemblies into contigs at runs of `N` before aligning. Use when your input is scaffolded. | <span class="src">`nextflow.config:40`</span> |
| `--mini_K` | `"500M"` | minimap2/winnowmap `-K` — the minibatch size. Raise it to improve throughput at the cost of memory. | <span class="src">`nextflow.config:47`</span> |
| `--stSort_m` | `"4G"` | `samtools sort -m`, memory **per thread**. Total sort memory is roughly `stSort_m × stSort_t`. | <span class="src">`nextflow.config:48`</span> |
| `--stSort_t` | `4` | `samtools sort -@`, sort threads. Note this is independent of the process `cpus`. | <span class="src">`nextflow.config:49`</span> |

---

## Stage B — repeat annotation

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--repeat_span_cutoff` | `0.80` | Minimum `total_repeat_span` — the fraction of the variant's sequence covered by the non-redundant union of RepeatMasker hits and ULTRA tandem repeats — for a variant to be kept. The single most consequential filter in the pipeline. Applied **twice**: once per contig chunk, once after concatenation. | <span class="src">`nextflow.config:51`</span> |
| `--tsd_win` | `30` | Width in bp of the flanking window searched for target site duplications. | <span class="src">`nextflow.config:44`</span> |
| `--tsd_batch_size` | `100` | Number of variants per TSD-search task. Lower it to increase parallelism, raise it to reduce scheduler overhead. | <span class="src">`nextflow.config:50`</span> |

!!! danger "`--tsd_win` is not safely tunable"
    `bin/TSD_Match_v2.sh` hardcodes the value `30` as the reference point in its TSD scoring
    arithmetic, independently of what `--tsd_win` is set to. Setting `--tsd_win` to anything other
    than `30` changes the window that is searched but **not** the scoring, silently corrupting the
    TSD score and therefore the `PASS`/`FAIL` decision. Leave it at the default until this is fixed.

### The trusted subset

With default settings GraffiTE writes `pangenome.trusted.vcf`, a conservative subset of
`pangenome.vcf`. These parameters define it. They have **no effect when `--human` is set**, because
`--human` replaces the trusted subset rather than refining it.

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--trusted_min_svlen` | `250` | Minimum \|SVLEN\| in bp. | <span class="src">`nextflow.config:52`</span> |
| `--trusted_max_ultra_span` | `0.6` | Maximum `ULTRA_TR_span` — reject variants that are mostly tandem repeat. Bypassed for records whose class is `Simple_repeat`. | <span class="src">`nextflow.config:53`</span> |
| `--trusted_ignore_filter` | `false` | When `true`, drop the requirement that the record already carries `FILTER=PASS` from the upstream caller. | <span class="src">`nextflow.config:54`</span> |

The full expression also requires `n_hits==1`, and requires `polyA="TRUE"` for LINE, SINE and
Retroposon classes. See [Stage B — annotation](../guides/annotation.md#the-trusted-subset).

---

## The `--human` pME subset

Setting `--human` swaps the trusted subset for a subset restricted to **recent, still-active human
mobile element subfamilies**. The output is `pangenome.human.vcf`; `pangenome.trusted.vcf` is not
written.

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--human` | `false` | Emit `pangenome.human.vcf` instead of `pangenome.trusted.vcf`, and run the HERV-K classifier. | <span class="src">`nextflow.config:55`</span> |
| `--human_alu_ids` | `"^AluY"` | Alu subfamilies to admit. The default keeps all `AluY*` (AluYa5, AluYb8, …) and excludes the older `AluS*` and `AluJ*`. | <span class="src">`nextflow.config:60`</span> |
| `--human_l1_ids` | `"^L1HS"` | L1 subfamilies. Add `^L1PA2` to relax by one subfamily. | <span class="src">`nextflow.config:61`</span> |
| `--human_sva_ids` | `"^SVA_[DEF]"` | SVA subfamilies. Also gates the `Simple_repeat` records that arise from VNTR-only SVA variants. | <span class="src">`nextflow.config:62`</span> |
| `--human_hervk_ids` | `"^HERVK-int,^LTR5_Hs,^LTR5A,^LTR5B"` | HML-2 lineage families. | <span class="src">`nextflow.config:63`</span> |
| `--human_min_svlen` | `250` | Minimum \|SVLEN\| in bp. | <span class="src">`nextflow.config:64`</span> |
| `--human_max_ultra_span` | `0.6` | Maximum `ULTRA_TR_span`. | <span class="src">`nextflow.config:65`</span> |
| `--human_ignore_filter` | `false` | Drop the `FILTER=PASS` requirement. | <span class="src">`nextflow.config:66`</span> |
| `--hervk_sva_pair` | `true` | Admit two-hit records that are a HERVK-int + SVA pair, which arise from LTR5_Hs/SVA sequence homology rather than from a genuine composite element. | <span class="src">`nextflow.config:67`</span> |
| `--hervk_pair_max_svlen` | `10500` | \|SVLEN\| ceiling for the above carve-out. Sized for the 9472 bp proviral element plus tolerance. | <span class="src">`nextflow.config:68`</span> |
| `--hervk_pair_max_hits` | `3` | Hit ceiling for the same carve-out. RepeatMasker splits the internal region of a degraded provirus as well as the LTR, giving three hits rather than two. Set to `2` for pre-v1.1 behaviour. | <span class="src">`nextflow.config:69`</span> |
| `--hervk_config` | `null` | Path to a JSON file overriding the HERV-K classifier's priors, sigmas and thresholds. Template at [`utils/HERVK.config.json`](https://github.com/cgroza/GraffiTE/blob/v1.1dev/utils/HERVK.config.json). | <span class="src">`nextflow.config:56`</span> |

!!! note "Whitelist syntax"
    The `*_ids` parameters are **comma-separated lists of bcftools regular expressions** matched
    against the `repeat_ids` INFO field. bcftools regexes support `^`, `$`, `.`, `*` and `[...]`
    but **not alternation** — which is exactly why these are lists rather than single patterns.
    An empty string keeps the whole class unfiltered.

See [Human mobile element insertions](../guides/human-mei.md) for the assembled filter expression
and the reasoning behind each gate.

---

## Stage C — genotyping

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--genotype` | `true` | Run Stage C at all. Set `false` to stop after annotation. | <span class="src">`nextflow.config:34`</span> |
| `--graph_method` | `"pangenie"` | One of `pangenie`, `giraffe`, `graphaligner`, `precomputed`. Determines how the graph is built and how reads are mapped onto it. | <span class="src">`nextflow.config:35`</span> |
| `--min_mapq` | `0` | `vg pack -Q` — minimum mapping quality for a read to contribute coverage. Ignored by the `pangenie` method. | <span class="src">`nextflow.config:71`</span> |
| `--min_support` | `"2,4"` | `vg call -m` — minimum support to call an allele, as `min_support_for_ref,min_support_for_alt`. Ignored by the `pangenie` method. | <span class="src">`nextflow.config:72`</span> |

!!! warning "`precomputed` is under-documented in the code"
    `--graph_method precomputed` is accepted by the branch test in `main.nf`, but it is absent from
    the error message listing valid methods and has no case in the `make_graph` or
    `graph_align_reads` switches. It is only usable in combination with `--graph` and either
    `--vcfs` or `--graph_alignments`. The comment in `nextflow.config:35` still lists only three
    methods.

### Reusing existing intermediates

All four are undeclared in `nextflow.config` but read by the workflow.

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--graffite_vcf` | `false` | Skip Stages A **and** B entirely and genotype this already-annotated GraffiTE VCF. | <span class="src">`nextflow.config:31`</span> |
| `--vcf` | `false` | Skip Stage A. Annotate this single external SV VCF. | <span class="src">`nextflow.config:32`</span> |
| `--RM_dir` | `false` | Skip RepeatMasker. Reuse an existing `2_Repeat_Filtering/` directory. Each subdirectory must contain `genotypes_repmasked_filtered.vcf` and `repeatmasker_dir`. Labelled "mainly for debug" in the config. | <span class="src">`nextflow.config:33`</span> |
| `--graph` | *(undeclared)* | Path to a pre-built graph index directory. Skips `make_graph`. | <span class="src">`main.nf:168`</span> |
| `--graph_alignments` | *(undeclared)* | Samplesheet of precomputed graph alignments (`sample,gaf,pack`). Skips `graph_align_reads`. | <span class="src">`main.nf:181`</span> |
| `--vcfs` | *(undeclared)* | Samplesheet of precomputed per-sample `vg call` VCFs. Skips both alignment and calling. | <span class="src">`main.nf:176`</span> |

!!! bug "`--vcf` cannot be combined with a discovery flag"
    When `--vcf` is set, Stage A is skipped — but the code then unconditionally reads the channel
    that Stage A would have created if any discovery flag is also present, and the run dies with
    `No such variable: sv_variants_ch`. Use `--vcf` alone, or use `--svs` instead if you want your
    VCFs merged with other discovery output. <span class="src">`main.nf:109`</span>

---

## Methylation (`--epigenomes`)

Requires the `panmethyl` submodule to be initialised.

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--epigenomes` | `false` | Project per-read methylation calls onto the graph during genotyping. Only compatible with the `giraffe` / `graphaligner` / `precomputed` methods. | <span class="src">`nextflow.config:30`</span> |
| `--code` | `"C+m"` | SAM `MM`/`ML` tag code identifying the modification to extract. | <span class="src">`nextflow.config:110`</span> |
| `--motif` | `"CG"` | Sequence motif indexed on the graph. | <span class="src">`nextflow.config:111`</span> |
| `--lifted` | *(undeclared)* | Samplesheet of already-lifted methylation CSVs, skipping extraction and lifting. | <span class="src">`main.nf:197`</span> |
| `--bed` | *(undeclared)* | BED file to project onto the graph and annotate with methylation. | <span class="src">`main.nf:211`</span> |

---

## Resources

### The global override

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--cores` | `false` | When set to an integer, overrides the `cpus` of **every** process that has a configurable thread count — ignoring all the individual `*_threads` parameters below. The simplest way to scale a run. | <span class="src">`nextflow.config:45`</span> |
| `--out` | `"out"` | Root directory for all published output. | <span class="src">`nextflow.config:43`</span> |

### Per-process allocation

Each process reads a `*_threads`, `*_memory` and `*_time` triple. `--cores`, if set, wins over
every `*_threads` value.

!!! warning "A `null` memory default means *no memory directive at all*"
    Seven of these default to `null`, which is not "unlimited" but "the process declares no memory
    requirement". On a scheduler that requires one, or for a memory-hungry step like `pangenie`,
    you must set it explicitly.

| Process | Threads | Memory | Time | Source |
|---|---|---|---|---|
| `map_asm` | `--map_asm_threads` `1` | `--map_asm_memory` **`null`** | `--map_asm_time` `"3h"` | <span class="src">`nextflow.config:81-83`</span> |
| `map_longreads` | `--map_longreads_threads` `1` | `--map_longreads_memory` **`null`** | `--map_longreads_time` `"12h"` | <span class="src">`nextflow.config:84-86`</span> |
| `sniffles_sample_call`, `sniffles_population_call` | `--sniffles_threads` `1` | `--sniffles_memory` **`null`** | `--sniffles_time` `"12h"` | <span class="src">`nextflow.config:95-97`</span> |
| `svim_asm`, `truvari_merge` | `--svim_asm_threads` `1` | `--svim_asm_memory` **`null`** | `--svim_asm_time` `"12h"` | <span class="src">`nextflow.config:100-102`</span> |
| `pav_asm` | `32` (literal; only `--cores` overrides) | `--pav_memory` `"120G"` | `--pav_time` `"12h"` | <span class="src">`nextflow.config:98-99, 257`</span> |
| `split_repeatmask`, `repeatmask_VCF`, `concat_repeatmask` | `--repeatmasker_threads` `1` | `--repeatmasker_memory` `"10G"` | `--repeatmasker_time` `"12h"` | <span class="src">`nextflow.config:92-94`</span> |
| `tsd_prep`, `tsd_search`, `tsd_report` | `1` (fixed) | `--tsd_memory` `"10G"` | `--tsd_time` `"1h"` | <span class="src">`nextflow.config:103-104`</span> |
| `pangenie_index`, `pangenie` | `--pangenie_threads` `1` | `--pangenie_memory` **`null`** | `--pangenie_time` `"12h"` | <span class="src">`nextflow.config:89-91`</span> |
| `make_graph` | `--make_graph_threads` `1` | `--make_graph_memory` `"40G"` | `--make_graph_time` `"6h"` | <span class="src">`nextflow.config:78-80`</span> |
| `bam_to_fastq`, `graph_align_reads` | `--graph_align_threads` `1` | `--graph_align_memory` **`null`** | `--graph_align_time` `"12h"` | <span class="src">`nextflow.config:75-77`</span> |
| `vg_call` | `--vg_call_threads` `1` | `--vg_call_memory` **`null`** | `--vg_call_time` `"2h"` | <span class="src">`nextflow.config:105-107`</span> |
| `merge_VCFs` | `1` (fixed) | `--merge_vcf_memory` `"10G"` | `--merge_vcf_time` `"1h"` | <span class="src">`nextflow.config:87-88`</span> |

The eight methylation processes carry fixed allocations (40–60 GB, 6 h) that are not parameterised.
<span class="src">`nextflow.config:215-253`</span>

!!! note "`truvari_merge` borrows the svim-asm knobs"
    `truvari_merge` has no parameters of its own — it reuses `--svim_asm_threads`, `--svim_asm_memory`
    and `--svim_asm_time`. This matters because the merge is now internally parallel: it shards with
    `truvari divide` and runs `truvari collapse` across `task.cpus` workers, so raising the svim-asm
    thread count speeds up merging too. <span class="src">`nextflow.config:144-148`</span>

### RepeatMasker threading

`bin/repmask_vcf.sh` does **not** use `task.cpus`. It computes its own thread count as
`nproc / 4`, floored at 1, and passes that to RepeatMasker's `-pa`. On a node where `nproc`
reports all host cores rather than the cores allocated to your job, this can oversubscribe badly.

---

## Undeclared parameters

Six parameters are read by `main.nf` but never declared in the `params { }` block. They work —
Groovy treats an unset property as null, which is falsy, so the `if (params.x)` guards behave
correctly — but they do not appear in `nextflow.config`, are invisible to `-params-file` schema
tooling, and are easy to miss when reading the config.

| Parameter | Read at | Purpose |
|---|---|---|
| `--svs` | `main.nf:90` | Samplesheet of external per-sample SV VCFs |
| `--graph` | `main.nf:168` | Pre-built graph index directory |
| `--vcfs` | `main.nf:176` | Precomputed per-sample `vg call` VCFs |
| `--graph_alignments` | `main.nf:181` | Precomputed graph alignments (`sample,gaf,pack`) |
| `--lifted` | `main.nf:197` | Pre-lifted methylation CSVs |
| `--bed` | `main.nf:211` | BED to project onto the graph |

---

## Deprecated and inert

| Parameter | Default | Status | Source |
|---|---|---|---|
| `--mammal` | `false` | **Inert since v1.1.** Still plumbed through to `bin/repmask_vcf.sh`, but that script now only prints a deprecation notice. L1 5′ inversions and SVA VNTR-only polymorphisms are reported unconditionally instead, via the `L1_5PINV` INFO field and the `(VNTR_only)` suffix on `repeat_ids`. The `mam_filter_1` and `mam_filter_2` fields no longer exist. | <span class="src">`nextflow.config:46`</span> |

---

## Nextflow's own options

Not GraffiTE parameters, but the ones you will use constantly. Single dash.

| Option | Effect |
|---|---|
| `-profile` | `standard` (local), `cluster` (SLURM), or `cloud` (AWS). See [Resources](../guides/resources.md#execution-profiles). |
| `-resume` | Reuse cached results from the previous run of the same pipeline in the same directory. |
| `-r` | Git revision — branch, tag or commit — when running straight from GitHub. |
| `-latest` | Pull the newest commit for `-r` before running. |
| `-with-report` | Write an HTML execution report, including per-process CPU and memory high-water marks. Useful for tuning the table above. |
| `-with-trace` | Write a tab-delimited trace of every task. |

Full list in the [Nextflow CLI documentation](https://www.nextflow.io/docs/latest/cli.html).
