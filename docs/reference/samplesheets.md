---
title: Samplesheet formats
description: The exact CSV columns each GraffiTE samplesheet parameter reads, with a copyable example for each.
---

# Samplesheet formats

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `4c8e385`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](../getting-started/v1.0-vs-v1.1.md).

Every samplesheet is a CSV read by `main.nf`. With one exception (`--pav`) the file has a header row
and columns are picked by name, so their order does not matter and extra columns are ignored. The
column names below are the ones the code reads; with a misspelt header the value is empty and the
run stops when it tries to open the file.

!!! warning "Use absolute paths"
    Nextflow resolves a relative path against the directory you launched from, whatever directory
    the CSV sits in. Absolute paths remove the question.

---

## `--assemblies`

Columns `sample`, `path`. One row per haploid assembly; a diploid sample gets two rows with two
names. <span class="src">`main.nf:105-106`</span>

```csv title="assemblies.csv"
sample,path
HG002_mat,/data/HG002.mat.fa.gz
HG002_pat,/data/HG002.pat.fa.gz
```

The `sample` name becomes the sample column of the svim-asm VCF and, after the merge, the
haplotype name in `pangenome.vcf`. <span class="src">`module/main.nf:153-154`</span>

## `--longreads`

Columns `sample`, `path`, `type`. Unaligned reads in FASTQ or FASTA, plain or gzipped. `type` picks
the minimap2 (or winnowmap) preset: it becomes `map-<type>`, except `lr:hq` which is passed as is.
So `ont`, `hifi` and `pb` are the usual values. <span class="src">`main.nf:87-88`, `module/main.nf:48-51`</span>

```csv title="longreads.csv"
sample,path,type
HG002,/data/HG002.hifi.fastq.gz,hifi
HG005,/data/HG005.ont.fastq.gz,ont
```

## `--bams`

Columns `sample`, `path`. Long reads already aligned to `--reference`, coordinate-sorted; the
process indexes the BAM before calling. <span class="src">`main.nf:93-94`, `module/main.nf:80`</span>

```csv title="bams.csv"
sample,path
HG002,/data/HG002.hifi.sorted.bam
```

## `--pav`

No header names are read. The first line is skipped as a header, the first column is the sample
name, and every further column is one haplotype FASTA of that sample. All haplotypes of a sample
share one row. <span class="src">`main.nf:115-116`, `module/main.nf:120-132`</span>

```csv title="pav.csv"
sample,hap1,hap2
HG002,/data/HG002.mat.fa.gz,/data/HG002.pat.fa.gz
```

## `--svs`

Columns `sample`, `path`. Per-sample SV VCFs from your own caller, bgzip-compressed and named
`*.vcf.gz`: the merge indexes each file with `tabix` and then loops over `*.vcf.gz`. `sample` is
read but not used. <span class="src">`main.nf:120-123`, `module/main.nf:181-194`</span>

```csv title="svs.csv"
sample,path
HG002,/data/HG002.sv.vcf.gz
```

Nothing on this path filters by `SVTYPE`, so anything other than insertions and deletions
reaches RepeatMasker. The other discovery backends keep `INS` and `DEL` only.
<span class="src">`module/main.nf:98,153`</span>

## `--genotype_with`

Columns `sample`, `path`, `type`. One read set per row, FASTQ (plain or gzipped) or BAM. `type` maps
to a read preset and anything else, including an empty value, means short reads:

| `type` | Preset | Effect |
|---|---|---|
| `pb`, `hifi` | `hifi` | `vg giraffe --parameter-preset hifi` |
| `ont` | `r10` | `vg giraffe --parameter-preset r10` |
| anything else | `default` | `vg giraffe` with `-i`, interleaved paired-end short reads |

<span class="src">`main.nf:189-205`, `module/main.nf:773-781`</span>

`pangenie` and `graphaligner` ignore the preset. A `.bam` path is converted to FASTQ first by
`bam_to_fastq`; a BAM is also what `--epigenomes` reads its modification tags from.
<span class="src">`main.nf:206-211,255`</span>

```csv title="reads.csv"
sample,path,type
HG002,/data/HG002.illumina.interleaved.fastq.gz,
HG005,/data/HG005.hifi.bam,hifi
```

## `--vcfs`

Columns `sample`, `path`. Per-sample `vg call` VCFs from an earlier run, bgzipped with the
`.vcf.gz` suffix, because `merge_VCFs` collects every `*vcf.gz` it is given. `path` has to be a
glob such as `HG002.vcf.gz*`: the code sorts the files it matches, and a plain path is taken apart
into its directory components instead (checked with Nextflow 26.04.6). The glob also brings the
`.tbi` index along. <span class="src">`main.nf:229-231`, `module/main.nf:830`</span>

```csv title="vcfs.csv"
sample,path
HG002,/data/earlier_run/HG002.vcf.gz*
```

## `--graph_alignments`

Columns `sample`, `gaf`, `pack`. The two files `graph_align_reads` publishes to
`GraffiTE_alignments/` for each sample. Needs `--graph` pointing at the index they were made
against. <span class="src">`main.nf:234-236`, `module/main.nf:763-769`</span>

```csv title="graph_alignments.csv"
sample,gaf,pack
HG002,/data/earlier_run/GraffiTE_alignments/HG002.gaf.gz,/data/earlier_run/GraffiTE_alignments/HG002.pack
```

## `--lifted`

Columns `sample`, `path`. Modification tables already lifted onto the graph, as `lift_epigenome`
writes them. Skips `bamtags_to_BED` and `lift_epigenome`. <span class="src">`main.nf:250-252`</span>

```csv title="lifted.csv"
sample,path
HG002,/data/earlier_run/lifted/HG002.csv.gz
```
