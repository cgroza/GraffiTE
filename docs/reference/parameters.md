---
title: Parameters
description: >-
  Every GraffiTE parameter, its exact default, what it controls, and the line in
  the source where it is declared or read.
---

# Parameters

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `4c8e385`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](../getting-started/v1.0-vs-v1.1.md).

Every default on this page was read from `nextflow.config` or from the code that consumes it. The
**Source** column gives the file and line so you can check.

## How to pass parameters

Nextflow distinguishes two kinds of option by the number of leading dashes:

| Form | Belongs to | Example |
|---|---|---|
| `--name value` | GraffiTE, anything on this page | `--assemblies assemblies.csv` |
| `-name value` | Nextflow itself | `-resume`, `-profile cluster`, `-with-report` |

```bash
nextflow run cgroza/GraffiTE -r v1.1dev -latest \
  -profile cluster \
  --reference hs37d5.fa \
  --assemblies assemblies.csv \
  --TE_library human_DFAM3.6.fasta \
  --genotype_with reads.csv
```

Defaults live in [`nextflow.config`](https://github.com/cgroza/GraffiTE/blob/v1.1dev/nextflow.config).
Override them on the command line, in a local copy of the config passed with `-c`, or in a
`-params-file`.

!!! warning "Use underscores, not hyphens"
    `--genotype-with x.csv` sets a parameter named `genotypeWith` and leaves `genotype_with` at
    its default, `reads.csv`; `--genotype_with x.csv` sets `genotype_with`. The pipeline ignores
    the hyphenated form without a message and reads the default samplesheet instead, or stops
    because `reads.csv` is missing from the launch directory. Older versions of the README wrote
    the hyphenated form. Verified with Nextflow 26.04.6.

---

## Required inputs

Almost every run needs three parameters, and all three default to a filename rather than a
value. If you do not set them, GraffiTE looks for `reference.fa` and `TE_library.fa` in the
launch directory and stops if they are absent.

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--reference` | `"reference.fa"` | Reference genome FASTA. Everything is called relative to it. Its existence is checked at launch. Plain or BGZF-compressed; a plain gzip is re-compressed with bgzip where an index is needed. | <span class="src">`nextflow.config:46`, `module/main.nf:534-541`</span> |
| `--TE_library` | `"TE_library.fa"` | FASTA of repeat consensus sequences, passed to RepeatMasker as `-lib`. Needed for Stage B, and again by the HERV-K step under `--human`, even with `--RM_dir`. | <span class="src">`nextflow.config:47`, `main.nf:144,175`</span> |
| `--genotype_with` | `"reads.csv"` | Samplesheet of the read sets to genotype. Read when `--genotype` is true (the default). See [Samplesheets](samplesheets.md). | <span class="src">`nextflow.config:40`, `main.nf:189`</span> |

---

## Stage A: discovery inputs

Supply at least one of these unless you enter further downstream with `--vcf`, `--RM_dir` or
`--graffite_vcf`. They add up: pass several and all of their calls go into one truvari merge.

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--assemblies` | `false` | Samplesheet of genome assemblies. Each is aligned with minimap2 (or winnowmap) and called with `svim-asm haploid`, keeping `INS` and `DEL` of 100 bp or more. | <span class="src">`nextflow.config:42`, `module/main.nf:153`</span> |
| `--longreads` | `false` | Samplesheet of unaligned long reads. Aligned, then called per sample and jointly with Sniffles2 at `--minsvlen 100`. | <span class="src">`nextflow.config:41`, `module/main.nf:81,97`</span> |
| `--bams` | `false` | Samplesheet of long-read BAMs that are already aligned to `--reference`. Skips alignment and goes straight to Sniffles2. Combines with `--longreads`. | <span class="src">`nextflow.config:29`, `main.nf:92-101`</span> |
| `--pav` | `false` | Samplesheet of phased assemblies to call with [PAV](https://github.com/EichlerLab/pav), which runs in its own container. Keeps variants with \|SVLEN\| above 50 bp. | <span class="src">`nextflow.config:43`, `module/main.nf:136`</span> |
| `--svs` | `false` | Samplesheet of per-sample SV VCFs you called yourself. No caller runs; the files go straight into the merge. | <span class="src">`nextflow.config:44`, `main.nf:120-123`</span> |

`--vcf` is the one input that does not combine. Passing it beside any of the five above stops the
run at launch with a message naming the flags. <span class="src">`main.nf:54-57`</span>

See [Stage A: discovery](../guides/discovery.md) for what each backend does.

### Discovery tuning

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--aligner` | `"minimap2"` | Aligner for assemblies and long reads. The only other accepted value is `"winnowmap"`. Any other value leaves `map_asm` and `map_longreads` with no script and the run fails. | <span class="src">`nextflow.config:119`, `module/main.nf:24-37`</span> |
| `--asm_divergence` | `"asm5"` | minimap2 `-x` preset for assembly alignment. Use `asm10` or `asm20` for assemblies further from the reference. Assemblies only. | <span class="src">`nextflow.config:118`, `module/main.nf:26`</span> |
| `--break_scaffolds` | `false` | Split each assembly into contigs at runs of `N` before aligning. For scaffolded input. | <span class="src">`nextflow.config:45`, `main.nf:107-109`</span> |
| `--mini_K` | `"500M"` | minimap2 and winnowmap `-K`, the number of bases loaded per batch. Larger is faster and uses more memory. | <span class="src">`nextflow.config:52`, `module/main.nf:26`</span> |
| `--stSort_m` | `"4G"` | `samtools sort -m`, memory per sort thread. Total sort memory is about `stSort_m` times `stSort_t`, on top of the aligner. | <span class="src">`nextflow.config:53`, `module/main.nf:27`</span> |
| `--stSort_t` | `4` | `samtools sort -@`, sort threads. Independent of the process `cpus`. | <span class="src">`nextflow.config:54`, `module/main.nf:27`</span> |

---

## Stage B: repeat annotation

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--repeat_span_cutoff` | `0.80` (fraction of variant length) | The filter keeps a variant when `total_repeat_span`, the fraction of its sequence covered by the union of RepeatMasker hits and ULTRA tandem repeats, is above this value. Applied twice: per contig chunk, and again after concatenation. | <span class="src">`nextflow.config:56`, `module/main.nf:533,615`</span> |
| `--tsd_win` | `30` (bp) | Width of the flank on each side of the variant, and of the variant end trimmed for the search, when looking for target site duplications. Sizes the sequences and the scoring alike. | <span class="src">`nextflow.config:49`, `module/main.nf:630,646`</span> |
| `--tsd_batch_size` | `100` (variants) | Variants per TSD-search task. Lower for more parallel tasks, higher for fewer. | <span class="src">`nextflow.config:55`, `main.nf:158`</span> |

### The trusted subset

Without `--human`, GraffiTE writes `pangenome.trusted.vcf`, a conservative subset of
`pangenome.vcf`. These parameters define it. They do nothing when `--human` is set, because
`--human` replaces the trusted subset instead of narrowing it.

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--trusted_min_svlen` | `250` (bp) | Minimum \|SVLEN\|. | <span class="src">`nextflow.config:57`, `module/main.nf:482`</span> |
| `--trusted_max_ultra_span` | `0.6` (fraction of variant length) | Maximum `ULTRA_TR_span`; rejects variants that are mostly tandem repeat. Records whose class is `Simple_repeat` bypass it. | <span class="src">`nextflow.config:58`, `module/main.nf:482`</span> |
| `--trusted_ignore_filter` | `false` | When true, the record no longer needs `FILTER=PASS` from the upstream caller. | <span class="src">`nextflow.config:59`, `module/main.nf:483`</span> |

The full expression also requires `n_hits==1`, and requires `polyA="TRUE"` for the LINE, SINE and
Retroposon classes. See [Stage B: annotation](../guides/annotation.md).

---

## The `--human` pME subset

`--human` swaps the trusted subset for one restricted to recent human mobile element subfamilies.
The output is `pangenome.human.vcf`; `pangenome.trusted.vcf` is not written. It also turns on the
HERV-K steps below.

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--human` | `false` | Write `pangenome.human.vcf` instead of `pangenome.trusted.vcf`, run `hervk_annotate`, and run `hervk_reconcile` after genotyping. | <span class="src">`nextflow.config:60`, `main.nf:171,282`</span> |
| `--human_alu_ids` | `"^AluY"` | Alu subfamilies to keep. The default keeps every `AluY*` and drops `AluS*` and `AluJ*`. | <span class="src">`nextflow.config:65`</span> |
| `--human_l1_ids` | `"^L1HS"` | L1 subfamilies. Add `^L1PA2` to widen by one subfamily. | <span class="src">`nextflow.config:66`</span> |
| `--human_sva_ids` | `"^SVA_[DEF]"` | SVA subfamilies. The same list gates the `Simple_repeat` records that come from VNTR-only SVA variants. | <span class="src">`nextflow.config:67`, `module/main.nf:495-496`</span> |
| `--human_hervk_ids` | `'^HERVK-int,^HERVK$,^LTR5_Hs,^LTR5A,^LTR5B'` | HML-2 lineage names. Both `HERVK-int` and a bare `HERVK` are listed because libraries differ; `^HERVK$` is anchored at both ends so that HERVK9, HERVK11 and HERVK14 stay out. Single-quoted in the config because `$` inside double quotes is a Groovy interpolation. | <span class="src">`nextflow.config:68-76`</span> |
| `--human_min_svlen` | `250` (bp) | Minimum \|SVLEN\|. | <span class="src">`nextflow.config:77`, `module/main.nf:498`</span> |
| `--human_max_ultra_span` | `0.6` (fraction of variant length) | Maximum `ULTRA_TR_span`; `Simple_repeat` records bypass it. | <span class="src">`nextflow.config:78`, `module/main.nf:498`</span> |
| `--human_ignore_filter` | `false` | When true, drop the `FILTER=PASS` requirement. | <span class="src">`nextflow.config:79`, `module/main.nf:527`</span> |
| `--hervk_sva_pair` | `true` | Also admit multi-hit records that pair `HERVK-int` with an SVA hit. RepeatMasker assigns part of the LTR5_Hs sequence to SVA, so a provirus arrives as two or three hits. | <span class="src">`nextflow.config:80`, `module/main.nf:524-525`</span> |
| `--hervk_pair_max_svlen` | `10500` (bp) | \|SVLEN\| ceiling for that carve-out: a 9472 bp provirus plus tolerance. | <span class="src">`nextflow.config:81`, `module/main.nf:524`</span> |
| `--hervk_pair_max_hits` | `3` (hits) | `n_hits` ceiling for the same carve-out. Three admits a provirus whose internal region RepeatMasker split in two; `2` restores the earlier behaviour. | <span class="src">`nextflow.config:82`, `module/main.nf:509-524`</span> |

!!! note "Whitelist syntax"
    The `*_ids` parameters are comma-separated lists of bcftools regular expressions matched against
    `repeat_ids`. bcftools regexes support `^`, `$`, `.`, `*` and `[...]` and have no alternation,
    which is why they are lists. `~` is applied to each element of the field, so `^` anchors per
    element. An empty string keeps the whole class. <span class="src">`nextflow.config:62-64`, `module/main.nf:487-491`</span>

See [Human mobile element insertions](../guides/human-mei.md) for the assembled expression.

### HERV-K classifier and locus reconciliation

All `--human` only. `hervk_annotate` runs after `concat_repeatmask` and reads the raw RepeatMasker
tables; `hervk_reconcile` runs after `merge_VCFs`.

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--hervk_config` | `null` | JSON file overriding the classifier's thresholds. Template at [`utils/HERVK.config.json`](https://github.com/cgroza/GraffiTE/blob/v1.1dev/utils/HERVK.config.json). | <span class="src">`nextflow.config:61`, `module/main.nf:276`</span> |
| `--hervk_max_svlen` | `25000` (bp) | \|SVLEN\| ceiling for HERV-K candidacy. Applied when the candidate list is built, before the reference windows are masked, so an oversized record costs no masking. | <span class="src">`nextflow.config:84`, `module/main.nf:300`</span> |
| `--hervk_ref_flank` | `1500` (bp) | Reference sequence masked on each side of a candidate footprint to establish the REF state. | <span class="src">`nextflow.config:108`, `module/main.nf:308,312`</span> |
| `--hervk_ref_annotation` | `null` | A precomputed RepeatMasker `.out` or BED for the reference. When set, the in-pipeline masking is skipped. | <span class="src">`nextflow.config:113`, `module/main.nf:305-308`</span> |
| `--hervk_locus_window` | `1200` (bp) | Largest gap between two record footprints that still counts as one locus: one LTR plus tolerance. | <span class="src">`nextflow.config:109`, `module/main.nf:340`</span> |
| `--hervk_strict` | `false` | Drop candidates classed `other` or below the confidence floor from `pangenome.human.vcf`. Off by default because dropping records is what hid a classifier failure before. | <span class="src">`nextflow.config:115`, `module/main.nf:277,330`</span> |
| `--hervk_reconcile` | `true` | Consolidate flagged HERV-K loci in the human subset of the genotyped calls. Only the giraffe back end is validated; the reconciler refuses others. | <span class="src">`nextflow.config:94`, `main.nf:282,298`</span> |
| `--hervk_reconcile_vcf` | `null` | Consolidate against this genotyped VCF from an earlier run instead of one produced now. Pair it with `--genotype false`. | <span class="src">`nextflow.config:90`, `main.nf:298-310`</span> |
| `--hervk_mask_graph_gt_at_cnv` | `true` | Withhold the graph genotypes at copy-number loci, where reads from the pre-existing reference copy give non-carriers ALT support. The calls are kept either way. | <span class="src">`nextflow.config:98`, `module/main.nf:394-396`</span> |

`--graffite_vcf` with `--human` and the default `--hervk_reconcile true` stops at launch, because
the reconciler needs the outputs of `hervk_annotate` and that process only runs during discovery.
Pass `--hervk_reconcile false`, or enter from `--RM_dir`. <span class="src">`main.nf:69-71`</span>

---

## Stage C: genotyping

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--genotype` | `true` | Run Stage C. Set `false` to stop after annotation. | <span class="src">`nextflow.config:34`, `main.nf:188`</span> |
| `--graph_method` | `"pangenie"` | One of `pangenie`, `giraffe`, `graphaligner`, `precomputed`. See below. | <span class="src">`nextflow.config:35`, `main.nf:214-217,274`</span> |
| `--min_mapq` | `0` | `vg pack -Q`, the lowest mapping quality a read needs to contribute coverage. Ignored by `pangenie`. | <span class="src">`nextflow.config:120`, `module/main.nf:782,790`</span> |
| `--min_support` | `"2,4"` | `vg call -m`, minimum support to call an allele, as `ref,alt`. Ignored by `pangenie`. | <span class="src">`nextflow.config:121`, `module/main.nf:811`</span> |

| `--graph_method` | Graph | Read mapping | Notes |
|---|---|---|---|
| `pangenie` | `PanGenie-index` on the annotated VCF | k-mer counting, no alignment | Written for short reads. The read preset, `--min_mapq` and `--min_support` do nothing. |
| `giraffe` | `vg autoindex` | `vg giraffe` | Short reads interleaved (`-i`) unless the samplesheet `type` selects a long-read preset. |
| `graphaligner` | `vg construct` | `GraphAligner` | Long reads. |
| `precomputed` | supplied with `--graph` | supplied with `--graph_alignments`, or skipped with `--vcfs` | Builds nothing. Without `--graph` and one of the two, the run stops at launch. |

<span class="src">`main.nf:61-64,214-243`, `module/main.nf:723-741,778-797`</span>

### Reusing existing intermediates

Each of these replaces the process that would produce it.

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--graffite_vcf` | `false` | Skip Stages A and B and genotype this `pangenome.vcf` from an earlier run. | <span class="src">`nextflow.config:31`, `main.nf:183-186`</span> |
| `--vcf` | `false` | Skip Stage A. Annotate this one merged SV VCF. Cannot be combined with a discovery flag. | <span class="src">`nextflow.config:32`, `main.nf:54-57,148-149`</span> |
| `--RM_dir` | `false` | Skip RepeatMasker. Reuse a `2_Repeat_Filtering/` directory whose subdirectories each hold `genotypes_repmasked_filtered.vcf` and `repeatmasker_dir/`. | <span class="src">`nextflow.config:33`, `main.nf:133-142`</span> |
| `--graph` | `false` | A graph index directory as `make_graph` writes it (`index.gfa`, `index.pb`, and `index.giraffe.gbz` for giraffe). Skips `make_graph`. | <span class="src">`nextflow.config:37`, `main.nf:221-225`</span> |
| `--graph_alignments` | `false` | Samplesheet of graph alignments (`sample,gaf,pack`). Skips `graph_align_reads`. | <span class="src">`nextflow.config:38`, `main.nf:234-236`</span> |
| `--vcfs` | `false` | Samplesheet of per-sample `vg call` VCFs. Skips alignment and `vg_call`. | <span class="src">`nextflow.config:39`, `main.nf:229-231`</span> |

See [Resuming and skipping work](../guides/skipping-work.md).

---

## Methylation (`--epigenomes`)

Needs the `panmethyl` submodule and one of the `giraffe`, `graphaligner` or `precomputed` methods;
the branch sits inside that block. <span class="src">`main.nf:38-42,245`</span>

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--epigenomes` | `false` | Lift per-read modification calls from the genotyping BAMs onto the graph and annotate the genotyped VCFs with them. | <span class="src">`nextflow.config:30`, `main.nf:245-269`</span> |
| `--code` | `"C+m"` | The SAM `MM`/`ML` modification code to extract. | <span class="src">`nextflow.config:164`, `main.nf:256`</span> |
| `--motif` | `"CG"` | Motif indexed on the graph. | <span class="src">`nextflow.config:165`, `main.nf:247`</span> |
| `--lifted` | `false` | Samplesheet of modification tables already lifted onto the graph; skips extraction and lifting. | <span class="src">`nextflow.config:166`, `main.nf:250-252`</span> |
| `--bed` | `false` | A BED file to project onto the graph and annotate with methylation levels. | <span class="src">`nextflow.config:167`, `main.nf:264-267`</span> |

See [Methylation](../guides/methylation.md).

---

## Resources

### Global switches

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--cores` | `false` | An integer here overrides the `cpus` of every process that reads a `*_threads` parameter, and the `32` of `pav_asm`. The processes fixed at one CPU are unaffected. | <span class="src">`nextflow.config:50,176-326`</span> |
| `--out` | `"out"` | Root of the published output tree. | <span class="src">`nextflow.config:48`</span> |

### Per-process allocation

Memory and time values are Nextflow strings (`"10G"`, `"12h"`). Where one parameter serves several processes, they are listed once.

!!! warning "A `null` memory default means no memory directive"
    Seven of these default to `null`, which means the process declares no memory requirement. A
    scheduler that needs one, or a memory-hungry step such as `pangenie`, needs the value set.

| Parameter | Default | Processes | Source |
|---|---|---|---|
| `--map_asm_threads` | `1` | `map_asm` | <span class="src">`nextflow.config:136,176`</span> |
| `--map_asm_memory` | `null` | `map_asm` | <span class="src">`nextflow.config:135,177`</span> |
| `--map_asm_time` | `"3h"` | `map_asm` | <span class="src">`nextflow.config:137,178`</span> |
| `--map_longreads_threads` | `1` | `map_longreads` | <span class="src">`nextflow.config:139,181`</span> |
| `--map_longreads_memory` | `null` | `map_longreads` | <span class="src">`nextflow.config:138,182`</span> |
| `--map_longreads_time` | `"12h"` | `map_longreads` | <span class="src">`nextflow.config:140,183`</span> |
| `--sniffles_threads` | `1` | `sniffles_sample_call`, `sniffles_population_call` | <span class="src">`nextflow.config:150,186,191`</span> |
| `--sniffles_memory` | `null` | same | <span class="src">`nextflow.config:149,187,192`</span> |
| `--sniffles_time` | `"12h"` | same | <span class="src">`nextflow.config:151,188,193`</span> |
| `--svim_asm_threads` | `1` | `svim_asm`, `truvari_merge` | <span class="src">`nextflow.config:155,196,201`</span> |
| `--svim_asm_memory` | `null` | same | <span class="src">`nextflow.config:154,197,202`</span> |
| `--svim_asm_time` | `"12h"` | same | <span class="src">`nextflow.config:156,198,203`</span> |
| `--pav_memory` | `"120G"` | `pav_asm` (its `cpus` is the literal `32` unless `--cores` is set) | <span class="src">`nextflow.config:152,323-324`</span> |
| `--pav_time` | `"12h"` | `pav_asm` | <span class="src">`nextflow.config:153,325`</span> |
| `--repeatmasker_threads` | `1` | `split_repeatmask`, `repeatmask_VCF`, `concat_repeatmask` | <span class="src">`nextflow.config:147,206,211,216`</span> |
| `--repeatmasker_memory` | `"10G"` | same | <span class="src">`nextflow.config:146,207,212,217`</span> |
| `--repeatmasker_time` | `"12h"` | same | <span class="src">`nextflow.config:148,208,213,218`</span> |
| `--tsd_memory` | `"10G"` | `tsd_prep`, `tsd_search`, `tsd_report` (one CPU each) | <span class="src">`nextflow.config:157,222,227,232`</span> |
| `--tsd_time` | `"1h"` | same | <span class="src">`nextflow.config:158,223,228,233`</span> |
| `--hervk_annotate_threads` | `1` | `hervk_annotate` | <span class="src">`nextflow.config:128,272`</span> |
| `--hervk_annotate_memory` | `"10G"` | `hervk_annotate` | <span class="src">`nextflow.config:127,273`</span> |
| `--hervk_annotate_time` | `"12h"` | `hervk_annotate` | <span class="src">`nextflow.config:129,274`</span> |
| `--hervk_reconcile_memory` | `"10G"` | `hervk_reconcile` (one CPU) | <span class="src">`nextflow.config:130,278`</span> |
| `--hervk_reconcile_time` | `"1h"` | `hervk_reconcile` | <span class="src">`nextflow.config:131,279`</span> |
| `--pangenie_threads` | `1` | `pangenie_index`, `pangenie` | <span class="src">`nextflow.config:144,236,241`</span> |
| `--pangenie_memory` | `null` | same | <span class="src">`nextflow.config:143,237,242`</span> |
| `--pangenie_time` | `"12h"` | same | <span class="src">`nextflow.config:145,238,243`</span> |
| `--make_graph_threads` | `1` | `make_graph` | <span class="src">`nextflow.config:133,246`</span> |
| `--make_graph_memory` | `"40G"` | `make_graph` | <span class="src">`nextflow.config:132,247`</span> |
| `--make_graph_time` | `"6h"` | `make_graph` | <span class="src">`nextflow.config:134,248`</span> |
| `--graph_align_threads` | `1` | `bam_to_fastq`, `graph_align_reads` | <span class="src">`nextflow.config:125,251,256`</span> |
| `--graph_align_memory` | `null` | same | <span class="src">`nextflow.config:124,252,257`</span> |
| `--graph_align_time` | `"12h"` | same | <span class="src">`nextflow.config:126,253,258`</span> |
| `--vg_call_threads` | `1` | `vg_call` | <span class="src">`nextflow.config:160,262`</span> |
| `--vg_call_memory` | `null` | `vg_call` | <span class="src">`nextflow.config:159,263`</span> |
| `--vg_call_time` | `"2h"` | `vg_call` | <span class="src">`nextflow.config:161,264`</span> |
| `--merge_vcf_memory` | `"10G"` | `merge_VCFs` (one CPU) | <span class="src">`nextflow.config:141,268`</span> |
| `--merge_vcf_time` | `"1h"` | `merge_VCFs` | <span class="src">`nextflow.config:142,269`</span> |

`break_scaffold` runs on one CPU with no memory or time directive. The eight methylation processes
have fixed allocations of 40 to 60 GB and 6 h that no parameter changes. See
[Resources and scaling](../guides/resources.md) for the full table.
<span class="src">`nextflow.config:172-174,281-320`</span>

!!! note "`truvari_merge` shares the svim-asm parameters"
    The merge has no parameters of its own. It shards the merged VCF with `truvari divide` and runs
    `truvari collapse` on `task.cpus` shards at once, so `--svim_asm_threads` also speeds up
    merging. <span class="src">`module/main.nf:202-209`</span>

### RepeatMasker threading

`bin/repmask_vcf.sh` ignores `task.cpus`. It passes `nproc / 4` (at least 1) to RepeatMasker's `-pa`
and `nproc` to ULTRA's `-t`. On a node where `nproc` reports the host's cores rather than the cores
allocated to the job, both run wider than the allocation. <span class="src">`bin/repmask_vcf.sh:22,30,39`</span>

---

## Deprecated and inert

| Parameter | Default | Status | Source |
|---|---|---|---|
| `--mammal` | `false` | Inert since v1.1. Still passed to `bin/repmask_vcf.sh`, which prints a notice and does nothing else. L1 5′ inversions and SVA VNTR-only polymorphisms are always reported, through the `L1_5PINV` INFO field and the `(VNTR_only)` suffix on `repeat_ids`. The `mam_filter_1` and `mam_filter_2` fields no longer exist. | <span class="src">`nextflow.config:51`, `bin/repmask_vcf.sh:112-116`</span> |
| `--hervk_mask_tandem` | `null` | Old name of `--hervk_mask_graph_gt_at_cnv`. When set, its value wins over the new name. | <span class="src">`nextflow.config:107`, `module/main.nf:394-395`</span> |

---

## Nextflow's own options

Single dash. Not GraffiTE parameters, but the ones you will use most.

| Option | Effect |
|---|---|
| `-profile` | `standard` (local), `cluster` (SLURM) or `cloud` (AWS). See [Resources and scaling](../guides/resources.md). |
| `-resume` | Reuse cached results from the previous run of the same pipeline in the same directory. |
| `-r` | Git revision (branch, tag or commit) when running straight from GitHub. |
| `-latest` | Pull the newest commit of `-r` before running. |
| `-with-report` | Write an HTML execution report with per-process CPU and memory peaks. Useful for tuning the table above. |
| `-with-trace` | Write a tab-delimited trace of every task. |
| `-c` | Add a config file, for site-specific executor settings. |

Full list in the [Nextflow CLI documentation](https://www.nextflow.io/docs/latest/cli.html).
