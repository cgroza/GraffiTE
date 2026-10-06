---
title: Parameters
description: >-
  Every GraffiTE parameter, its exact default, what it controls, and the line in
  the source where it is declared or read.
---

# Parameters

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `25e417a`. The
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
| `--reference` | `"reference.fa"` | Reference genome FASTA. Everything is called relative to it. Its existence is checked at launch. Plain or BGZF-compressed; a plain gzip is re-compressed with bgzip where an index is needed. | <span class="src">`nextflow.config:46`, `module/main.nf:578-585`</span> |
| `--TE_library` | `"TE_library.fa"` | FASTA of repeat consensus sequences, passed to RepeatMasker as `-lib`. Needed for Stage B, and again by the HERV-K step under `--human`, even with `--RM_dir`. | <span class="src">`nextflow.config:48`, `main.nf:144,175`</span> |
| `--genotype_with` | `"reads.csv"` | Samplesheet of the read sets to genotype. Read when `--genotype` is true (the default). See [Samplesheets](samplesheets.md). | <span class="src">`nextflow.config:39`, `main.nf:189`</span> |

---

## Stage A: discovery inputs

Supply at least one of these unless you enter further downstream with `--vcf`, `--RM_dir` or
`--graffite_vcf`. They add up: pass several and all of their calls go into one truvari merge.

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--assemblies` | `false` | Samplesheet of genome assemblies. Each is aligned with minimap2 (or winnowmap) and called with `svim-asm haploid`, keeping `INS` and `DEL` of 100 bp or more. | <span class="src">`nextflow.config:41`, `module/main.nf:184`</span> |
| `--longreads` | `false` | Samplesheet of unaligned long reads. Aligned, then called per sample and jointly with Sniffles2 at `--minsvlen 100`. | <span class="src">`nextflow.config:40`, `module/main.nf:110,128`</span> |
| `--bams` | `false` | Samplesheet of long-read BAMs that are already aligned to `--reference`. Skips alignment and goes straight to Sniffles2. Combines with `--longreads`. | <span class="src">`nextflow.config:28`, `main.nf:92-101`</span> |
| `--pav` | `false` | Samplesheet of phased assemblies to call with [PAV](https://github.com/EichlerLab/pav), which runs in its own container. Keeps variants with \|SVLEN\| above 50 bp. | <span class="src">`nextflow.config:42`, `module/main.nf:167`</span> |
| `--svs` | `false` | Samplesheet of per-sample SV VCFs you called yourself. No caller runs; the files go straight into the merge. | <span class="src">`nextflow.config:43`, `main.nf:120-123`</span> |

`--vcf` is the one input that does not combine. Passing it beside any of the five above stops the
run at launch with a message naming the flags. <span class="src">`main.nf:46-49`</span>

See [Stage A: discovery](../guides/discovery.md) for what each backend does.

### Discovery tuning

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--aligner` | `"minimap2"` | Aligner for assemblies and long reads. The only other accepted value is `"winnowmap"`. Any other value leaves `map_asm` and `map_longreads` with no script and the run fails. | <span class="src">`nextflow.config:120`, `module/main.nf:48-61`</span> |
| `--asm_divergence` | `"asm5"` | minimap2 `-x` preset for assembly alignment. Use `asm10` or `asm20` for assemblies further from the reference. Assemblies only. | <span class="src">`nextflow.config:119`, `module/main.nf:50`</span> |
| `--break_scaffolds` | `false` | Split each assembly into contigs at runs of at least `--break_scaffolds_min_gap` `N` before aligning. For scaffolded input. | <span class="src">`nextflow.config:44`, `main.nf:107-109`</span> |
| `--break_scaffolds_min_gap` | `10` | Shortest run of `N` that `--break_scaffolds` treats as a gap. A shorter run is an unknown base and stays inside its contig, so an insertion that carries one is not cut. `1` splits at every `N`. | <span class="src">`nextflow.config:45`, `module/main.nf:36`, `bin/breakgaps.py:14-15`</span> |
| `--mini_K` | `"500M"` | minimap2 and winnowmap `-K`, the number of bases loaded per batch. Larger is faster and uses more memory. | <span class="src">`nextflow.config:53`, `module/main.nf:50`</span> |
| `--stSort_m` | `"4G"` | `samtools sort -m`, memory per sort thread. Total sort memory is about `stSort_m` times `stSort_t`, on top of the aligner. | <span class="src">`nextflow.config:54`, `module/main.nf:51`</span> |
| `--stSort_t` | `4` | `samtools sort -@`, sort threads. Independent of the process `cpus`. | <span class="src">`nextflow.config:55`, `module/main.nf:51`</span> |

---

## Stage B: repeat annotation

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--repeat_span_cutoff` | `0.80` (fraction of variant length) | The filter keeps a variant when `total_repeat_span`, the fraction of its sequence covered by the union of RepeatMasker hits and ULTRA tandem repeats, is above this value. Applied twice: per contig chunk, and again after concatenation. | <span class="src">`nextflow.config:57`, `module/main.nf:577,659`</span> |
| `--tsd_win` | `30` (bp) | Width of the flank on each side of the variant, and of the variant end trimmed for the search, when looking for target site duplications. Sizes the sequences and the scoring alike. | <span class="src">`nextflow.config:50`, `module/main.nf:673,688`</span> |
| `--tsd_batch_size` | `100` (variants) | Variants per TSD-search task. Lower for more parallel tasks, higher for fewer. | <span class="src">`nextflow.config:56`, `main.nf:158`</span> |

### The trusted subset

Without `--human`, GraffiTE writes `pangenome.trusted.vcf`, a conservative subset of
`pangenome.vcf`. These parameters define it. With `--human`, GraffiTE writes neither
`pangenome.trusted.vcf` nor `GraffiTE.merged.genotypes.trusted.vcf.gz`, because `--human` replaces
the trusted subset instead of narrowing it. `genotyping_audit` still reads all three on every
genotyping run to fill the `trusted` column of `genotyping_record_audit.tsv`.

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--trusted_min_svlen` | `250` (bp) | Minimum \|SVLEN\|. Sized for insertions; it truncates the `SVA_*(VNTR_only)` records, whose unit is about 49 bp. See [SVA VNTR polymorphisms](../background/sva-vntr.md#the-subsets-truncate-this-set). | <span class="src">`nextflow.config:58`, `module/main.nf:14`</span> |
| `--trusted_max_ultra_span` | `0.6` (fraction of variant length) | Maximum `ULTRA_TR_span`; rejects variants that are mostly tandem repeat. Records whose class is `Simple_repeat` bypass it. | <span class="src">`nextflow.config:59`, `module/main.nf:14`</span> |
| `--trusted_ignore_filter` | `false` | When true, the record no longer needs `FILTER=PASS` from the upstream caller. | <span class="src">`nextflow.config:60`, `module/main.nf:22`</span> |

The full expression also requires `n_hits==1`, and requires `polyA="TRUE"` for the LINE, SINE and
Retroposon classes. See [Stage B: annotation](../guides/annotation.md).

---

## The `--human` pME subset

`--human` swaps the trusted subset for one restricted to recent human mobile element subfamilies.
The output is `pangenome.human.vcf`; `pangenome.trusted.vcf` is not written. It also turns on the
HERV-K steps below.

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--human` | `false` | Write `pangenome.human.vcf` instead of `pangenome.trusted.vcf`, run `hervk_annotate`, and run `hervk_reconcile` after genotyping. | <span class="src">`nextflow.config:61`, `main.nf:171,286`</span> |
| `--human_alu_ids` | `"^AluY"` | Alu subfamilies to keep. The default keeps every `AluY*` and drops `AluS*` and `AluJ*`. | <span class="src">`nextflow.config:66`</span> |
| `--human_l1_ids` | `"^L1HS"` | L1 subfamilies. Add `^L1PA2` to widen by one subfamily. | <span class="src">`nextflow.config:67`</span> |
| `--human_sva_ids` | `"^SVA_[DEF]"` | SVA subfamilies. The same list gates the `Simple_repeat` records that come from VNTR-only SVA variants, where it matches the consensus the VNTR sequence scored against rather than the host element’s subfamily. See [SVA VNTR polymorphisms](../background/sva-vntr.md#the-subsets-truncate-this-set). | <span class="src">`nextflow.config:68`, `module/main.nf:539-540`</span> |
| `--human_hervk_ids` | `'^HERVK-int,^HERVK$,^LTR5_Hs,^LTR5A,^LTR5B'` | HML-2 lineage names. Both `HERVK-int` and a bare `HERVK` are listed because libraries differ; `^HERVK$` is anchored at both ends so that HERVK9, HERVK11 and HERVK14 stay out. Single-quoted in the config because `$` inside double quotes is a Groovy interpolation. | <span class="src">`nextflow.config:69-77`</span> |
| `--human_min_svlen` | `250` (bp) | Minimum \|SVLEN\|. Sized for insertions; it truncates the `SVA_*(VNTR_only)` records, whose unit is about 49 bp. See [SVA VNTR polymorphisms](../background/sva-vntr.md#the-subsets-truncate-this-set). | <span class="src">`nextflow.config:78`, `module/main.nf:542`</span> |
| `--human_max_ultra_span` | `0.6` (fraction of variant length) | Maximum `ULTRA_TR_span`; `Simple_repeat` records bypass it. | <span class="src">`nextflow.config:79`, `module/main.nf:542`</span> |
| `--human_ignore_filter` | `false` | When true, drop the `FILTER=PASS` requirement. | <span class="src">`nextflow.config:80`, `module/main.nf:571`</span> |
| `--hervk_sva_pair` | `true` | Also admit multi-hit records that pair `HERVK-int` with an SVA hit. RepeatMasker assigns part of the LTR5_Hs sequence to SVA, so a provirus arrives as two or three hits. | <span class="src">`nextflow.config:81`, `module/main.nf:568-569`</span> |
| `--hervk_pair_max_svlen` | `10500` (bp) | \|SVLEN\| ceiling for that carve-out: a 9472 bp provirus plus tolerance. | <span class="src">`nextflow.config:82`, `module/main.nf:568`</span> |
| `--hervk_pair_max_hits` | `3` (hits) | `n_hits` ceiling for the same carve-out. Three admits a provirus whose internal region RepeatMasker split in two; `2` restores the earlier behaviour. | <span class="src">`nextflow.config:83`, `module/main.nf:553-568`</span> |

!!! note "Whitelist syntax"
    The `*_ids` parameters are comma-separated lists of bcftools regular expressions matched against
    `repeat_ids`. bcftools regexes support `^`, `$`, `.`, `*` and `[...]` and have no alternation,
    which is why they are lists. `~` is applied to each element of the field, so `^` anchors per
    element. An empty string keeps the whole class. <span class="src">`nextflow.config:63-65`, `module/main.nf:531-535`</span>

See [Human mobile element insertions](../guides/human-mei.md) for the assembled expression.

### HERV-K classifier and locus reconciliation

All `--human` only. `hervk_annotate` runs after `concat_repeatmask` and reads the raw RepeatMasker
tables; `hervk_reconcile` runs after `merge_VCFs`.

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--hervk_config` | `null` | JSON file overriding the classifier's thresholds. Template at [`utils/HERVK.config.json`](https://github.com/cgroza/GraffiTE/blob/v1.1dev/utils/HERVK.config.json). | <span class="src">`nextflow.config:62`, `module/main.nf:321`</span> |
| `--hervk_max_svlen` | `25000` (bp) | \|SVLEN\| ceiling for HERV-K candidacy. Applied when the candidate list is built, before the reference windows are masked, so an oversized record costs no masking. | <span class="src">`nextflow.config:85`, `module/main.nf:345`</span> |
| `--hervk_ref_flank` | `1500` (bp) | Reference sequence masked on each side of a candidate footprint to establish the REF state. | <span class="src">`nextflow.config:109`, `module/main.nf:353,357`</span> |
| `--hervk_ref_annotation` | `null` | A precomputed RepeatMasker `.out` or BED for the reference. When set, the in-pipeline masking is skipped. | <span class="src">`nextflow.config:114`, `module/main.nf:350-353`</span> |
| `--hervk_locus_window` | `1200` (bp) | Largest gap between two record footprints that still counts as one locus: one LTR plus tolerance. | <span class="src">`nextflow.config:110`, `module/main.nf:385`</span> |
| `--hervk_strict` | `false` | Drop candidates classed `other` or below the confidence floor from `pangenome.human.vcf`. Off by default because dropping records is what hid a classifier failure before. | <span class="src">`nextflow.config:116`, `module/main.nf:322,375`</span> |
| `--hervk_reconcile` | `true` | Consolidate flagged HERV-K loci in the human subset of the genotyped calls. Reads `vg call` genotypes: giraffe, graphaligner or precomputed. | <span class="src">`nextflow.config:95`, `main.nf:286,302`, `bin/hervk_reconcile.py:434`</span> |
| `--hervk_reconcile_vcf` | `null` | Consolidate against this genotyped VCF from an earlier run instead of one produced now. Pair it with `--genotype false`. | <span class="src">`nextflow.config:91`, `main.nf:302-314`</span> |
| `--hervk_mask_graph_gt_at_cnv` | `true` | Withhold the graph genotypes at copy-number loci, where reads from the pre-existing reference copy give non-carriers ALT support. The calls are kept either way. | <span class="src">`nextflow.config:99`, `module/main.nf:439-441`</span> |

`--graffite_vcf` with `--human` and the default `--hervk_reconcile true` stops at launch, because
the reconciler needs the outputs of `hervk_annotate` and that process only runs during discovery.
Pass `--hervk_reconcile false`, or enter from `--RM_dir`. <span class="src">`main.nf:61-63`</span>

A `--human` run that genotypes with `--graph_method pangenie` also stops at launch while
`--hervk_reconcile` is `true`, because the reconciler reads `vg call` genotypes and PanGenie does
not produce them. Pass `--hervk_reconcile false`, or genotype with giraffe, graphaligner or
precomputed. <span class="src">`main.nf:65-71`</span>

---

## Stage C: genotyping

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--genotype` | `true` | Run Stage C. Set `false` to stop after annotation. | <span class="src">`nextflow.config:33`, `main.nf:188`</span> |
| `--graph_method` | `"pangenie"` | One of `pangenie`, `giraffe`, `graphaligner`, `precomputed`. See below. | <span class="src">`nextflow.config:34`, `main.nf:200-206,263`</span> |
| `--min_mapq` | `0` | `vg pack -Q`, the lowest mapping quality a read needs to contribute coverage. Ignored by `pangenie`. | <span class="src">`nextflow.config:121`, `module/main.nf:843,851`</span> |
| `--min_support` | `"2,4"` | `vg call -m`, minimum support to call an allele, as `ref,alt`. Ignored by `pangenie`. | <span class="src">`nextflow.config:122`, `module/main.nf:872`</span> |

| `--graph_method` | Graph | Read mapping | Notes |
|---|---|---|---|
| `pangenie` | `PanGenie-index` on the annotated VCF | k-mer counting, no alignment | Written for short reads. The read preset, `--min_mapq` and `--min_support` do nothing. |
| `giraffe` | `vg autoindex` | `vg giraffe` | Short reads interleaved (`-i`) unless the samplesheet `type` selects a long-read preset. |
| `graphaligner` | `vg construct` | `GraphAligner` | Long reads. |
| `precomputed` | supplied with `--graph` | supplied with `--graph_alignments`, or skipped with `--vcfs` | Builds nothing. Without `--graph` and one of the two, the run stops at launch. |

<span class="src">`main.nf:53-56,200-232`, `module/main.nf:786-799,840-858`</span>

### Reusing existing intermediates

Each of these replaces the process that would produce it.

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--graffite_vcf` | `false` | Skip Stages A and B and genotype this `pangenome.vcf` from an earlier run. | <span class="src">`nextflow.config:30`, `main.nf:183-186`</span> |
| `--vcf` | `false` | Skip Stage A. Annotate this one merged SV VCF. Cannot be combined with a discovery flag. | <span class="src">`nextflow.config:31`, `main.nf:46-49,148-149`</span> |
| `--RM_dir` | `false` | Skip RepeatMasker. Reuse a `2_Repeat_Filtering/` directory whose subdirectories each hold `genotypes_repmasked_filtered.vcf` and `repeatmasker_dir/`. | <span class="src">`nextflow.config:32`, `main.nf:133-142`</span> |
| `--graph` | `false` | A graph index directory as `make_graph` writes it (`index.gfa`, `index.pb`, and `index.giraffe.gbz` for giraffe). Skips `make_graph`. | <span class="src">`nextflow.config:36`, `main.nf:210-214`</span> |
| `--graph_alignments` | `false` | Samplesheet of graph alignments (`sample,gaf,pack`). Skips `graph_align_reads`. | <span class="src">`nextflow.config:37`, `main.nf:223-225`</span> |
| `--vcfs` | `false` | Samplesheet of per-sample `vg call` VCFs. Skips alignment and `vg_call`. | <span class="src">`nextflow.config:38`, `main.nf:218-220`</span> |

See [Resuming and skipping work](../guides/skipping-work.md).

---

## Methylation (`--epigenomes`)

Needs the `panmethyl` submodule and one of the `giraffe`, `graphaligner` or `precomputed` methods;
the branch sits inside that block. <span class="src">`main.nf:12,234`</span>

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--epigenomes` | `false` | Lift per-read modification calls from the genotyping BAMs onto the graph and annotate the genotyped VCFs with them. | <span class="src">`nextflow.config:29`, `main.nf:234-258`</span> |
| `--code` | `"C+m"` | The SAM `MM`/`ML` modification code to extract. | <span class="src">`nextflow.config:165`, `main.nf:245`</span> |
| `--motif` | `"CG"` | Motif indexed on the graph. | <span class="src">`nextflow.config:166`, `main.nf:236`</span> |
| `--lifted` | `false` | Samplesheet of modification tables already lifted onto the graph; skips extraction and lifting. | <span class="src">`nextflow.config:167`, `main.nf:239-241`</span> |
| `--bed` | `false` | A BED file to project onto the graph and annotate with methylation levels. | <span class="src">`nextflow.config:168`, `main.nf:253-256`</span> |

See [Methylation](../guides/methylation.md).

---

## Resources

### Global switches

| Parameter | Default | Effect | Source |
|---|---|---|---|
| `--cores` | `false` | An integer here overrides the `cpus` of every process that reads a `*_threads` parameter, and the `32` of `pav_asm`. The processes fixed at one CPU are unaffected. | <span class="src">`nextflow.config:51,179-345`</span> |
| `--out` | `"out"` | Root of the published output tree. | <span class="src">`nextflow.config:49`</span> |
| `--container_tmp` | `false` | Directory bound to `/tmp` inside the container. The launch directory when unset. Set it when the launch filesystem is small, slow, or unwritable from the container. | <span class="src">`nextflow.config:47,176`</span> |

### Per-process allocation

Memory and time values are Nextflow strings (`"10G"`, `"12h"`). Where one parameter serves several processes, they are listed once.

!!! warning "A `null` memory default means no memory directive"
    Seven of these default to `null`, which means the process declares no memory requirement. A
    scheduler that needs one, or a memory-hungry step such as `pangenie`, needs the value set.

| Parameter | Default | Processes | Source |
|---|---|---|---|
| `--map_asm_threads` | `1` | `map_asm` | <span class="src">`nextflow.config:137,184`</span> |
| `--map_asm_memory` | `null` | `map_asm` | <span class="src">`nextflow.config:136,185`</span> |
| `--map_asm_time` | `"3h"` | `map_asm` | <span class="src">`nextflow.config:138,186`</span> |
| `--map_longreads_threads` | `1` | `map_longreads` | <span class="src">`nextflow.config:140,189`</span> |
| `--map_longreads_memory` | `null` | `map_longreads` | <span class="src">`nextflow.config:139,190`</span> |
| `--map_longreads_time` | `"12h"` | `map_longreads` | <span class="src">`nextflow.config:141,191`</span> |
| `--sniffles_threads` | `1` | `sniffles_sample_call`, `sniffles_population_call` | <span class="src">`nextflow.config:151,194,199`</span> |
| `--sniffles_memory` | `null` | same | <span class="src">`nextflow.config:150,195,200`</span> |
| `--sniffles_time` | `"12h"` | same | <span class="src">`nextflow.config:152,196,201`</span> |
| `--svim_asm_threads` | `1` | `svim_asm`, `truvari_merge` | <span class="src">`nextflow.config:156,204,209`</span> |
| `--svim_asm_memory` | `null` | same | <span class="src">`nextflow.config:155,205,210`</span> |
| `--svim_asm_time` | `"12h"` | same | <span class="src">`nextflow.config:157,206,211`</span> |
| `--pav_memory` | `"120G"` | `pav_asm` (its `cpus` is the literal `32` unless `--cores` is set) | <span class="src">`nextflow.config:153,341-342`</span> |
| `--pav_time` | `"12h"` | `pav_asm` | <span class="src">`nextflow.config:154,343`</span> |
| `--repeatmasker_threads` | `1` | `split_repeatmask`, `repeatmask_VCF`, `concat_repeatmask` | <span class="src">`nextflow.config:148,214,219,224`</span> |
| `--repeatmasker_memory` | `"10G"` | same | <span class="src">`nextflow.config:147,215,220,225`</span> |
| `--repeatmasker_time` | `"12h"` | same | <span class="src">`nextflow.config:149,216,221,226`</span> |
| `--tsd_memory` | `"10G"` | `tsd_prep`, `tsd_search`, `tsd_report` (one CPU each) | <span class="src">`nextflow.config:158,230,235,240`</span> |
| `--tsd_time` | `"1h"` | same | <span class="src">`nextflow.config:159,231,236,241`</span> |
| `--hervk_annotate_threads` | `1` | `hervk_annotate` | <span class="src">`nextflow.config:129,290`</span> |
| `--hervk_annotate_memory` | `"10G"` | `hervk_annotate` | <span class="src">`nextflow.config:128,291`</span> |
| `--hervk_annotate_time` | `"12h"` | `hervk_annotate` | <span class="src">`nextflow.config:130,292`</span> |
| `--hervk_reconcile_memory` | `"10G"` | `hervk_reconcile` (one CPU) | <span class="src">`nextflow.config:131,296`</span> |
| `--hervk_reconcile_time` | `"1h"` | `hervk_reconcile` | <span class="src">`nextflow.config:132,297`</span> |
| `--pangenie_threads` | `1` | `pangenie_index`, `pangenie` | <span class="src">`nextflow.config:145,244,249`</span> |
| `--pangenie_memory` | `null` | same | <span class="src">`nextflow.config:144,245,250`</span> |
| `--pangenie_time` | `"12h"` | same | <span class="src">`nextflow.config:146,246,251`</span> |
| `--make_graph_threads` | `1` | `make_graph` | <span class="src">`nextflow.config:134,254`</span> |
| `--make_graph_memory` | `"40G"` | `make_graph` | <span class="src">`nextflow.config:133,255`</span> |
| `--make_graph_time` | `"6h"` | `make_graph` | <span class="src">`nextflow.config:135,256`</span> |
| `--graph_align_threads` | `1` | `bam_to_fastq`, `graph_align_reads` | <span class="src">`nextflow.config:126,259,264`</span> |
| `--graph_align_memory` | `null` | same | <span class="src">`nextflow.config:125,260,265`</span> |
| `--graph_align_time` | `"12h"` | same | <span class="src">`nextflow.config:127,261,266`</span> |
| `--vg_call_threads` | `1` | `vg_call` | <span class="src">`nextflow.config:161,270`</span> |
| `--vg_call_memory` | `null` | `vg_call` | <span class="src">`nextflow.config:160,271`</span> |
| `--vg_call_time` | `"2h"` | `vg_call` | <span class="src">`nextflow.config:162,272`</span> |
| `--merge_vcf_memory` | `"10G"` | `merge_VCFs`, `trusted_genotypes`, `genotyping_audit` (one CPU each) | <span class="src">`nextflow.config:142,276,281,286`</span> |
| `--merge_vcf_time` | `"1h"` | same | <span class="src">`nextflow.config:143,277,282,287`</span> |

`break_scaffold` runs on one CPU with no memory or time directive. The eight methylation processes
have fixed allocations of 40 to 60 GB and 6 h that no parameter changes. See
[Resources and scaling](../guides/resources.md) for the full table.
<span class="src">`nextflow.config:180-182,299-338`</span>

!!! note "`truvari_merge` shares the svim-asm parameters"
    The merge has no parameters of its own. It shards the merged VCF with `truvari divide` and runs
    `truvari collapse` on `task.cpus` shards at once, so `--svim_asm_threads` also speeds up
    merging. <span class="src">`module/main.nf:247-254`</span>

### RepeatMasker threading

`bin/repmask_vcf.sh` ignores `task.cpus`. It passes `nproc / 4` (at least 1) to RepeatMasker's `-pa`
and `nproc` to ULTRA's `-t`. On a node where `nproc` reports the host's cores rather than the cores
allocated to the job, both run wider than the allocation. <span class="src">`bin/repmask_vcf.sh:67,75,92`</span>

---

## Deprecated and inert

| Parameter | Default | Status | Source |
|---|---|---|---|
| `--mammal` | `false` | Inert since v1.1. Still passed to `bin/repmask_vcf.sh`, which prints a notice and does nothing else. L1 5′ inversions and SVA VNTR-only polymorphisms are always reported, through the `L1_5PINV` INFO field and the `(VNTR_only)` suffix on `repeat_ids`. The `mam_filter_1` and `mam_filter_2` fields no longer exist. | <span class="src">`nextflow.config:52`, `bin/repmask_vcf.sh:165-169`</span> |
| `--hervk_mask_tandem` | `null` | Old name of `--hervk_mask_graph_gt_at_cnv`. When set, its value wins over the new name. | <span class="src">`nextflow.config:108`, `module/main.nf:439-440`</span> |

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
