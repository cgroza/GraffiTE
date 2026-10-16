---
title: Quickstart
description: >-
  Run GraffiTE end to end on the bundled human chromosome 22 test set, and know
  which files to look for when it finishes.
---

# Quickstart

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `e687712`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](v1.0-vs-v1.1.md).

The repository ships a small human dataset that exercises every stage: one haplotype assembly of
chromosome 22, the matching reference, a Dfam library, and a short-read set to genotype.

---

## Test dataset

```bash
tar xzf test/human_test_set.tar.gz
cd GraffiTE_testset
```

| File | Role |
|---|---|
| `HG002.mat.cur.20211005_chr22.fasta.gz` | The maternal HG002 assembly, chromosome 22 only. Stage 1 input. |
| `hs37d5.chr22.fa` | Reference, chromosome 22 of hs37d5. |
| `human_DFAM3.6.fasta` | TE consensus library for RepeatMasker. |
| `HG002.100k.set1.Illumina.fastq.gz` | 100 k Illumina reads from HG002. Stage 3 input. |
| `assemblies.csv` | `path,sample` samplesheet naming the assembly `HG002.mat`. |
| `reads.csv` | `path,sample` samplesheet naming the read set `short_test1`. |
| `out/` | The output of the command below, made with GraffiTE 1.1.0.0, to compare with yours. |

`reads.csv` has no `type` column. The workflow reads `row.type`, finds nothing, and falls to
the `default` preset <span class="src">`main.nf:190`</span>; PanGenie, the default
genotyper, does not use the preset.

!!! note "Compare your run with the bundled `out/`"
    `out/` holds the output of the command below, made with GraffiTE 1.1.0.0 on 2026-10-16.
    Your `out_v1.1/` should hold the same records; the `##fileDate` header lines differ.

---

## Run it

```bash
nextflow run cgroza/GraffiTE -r v1.1dev \
  --assemblies assemblies.csv \
  --reference hs37d5.chr22.fa \
  --TE_library human_DFAM3.6.fasta \
  --genotype_with reads.csv \
  --out out_v1.1
```

Pass the four data arguments as shown; the defaults for `--reference` and `--TE_library` are
placeholder filenames that do not exist <span class="src">`nextflow.config:46-48`</span>.
`--genotype_with` is the parameter name as declared <span class="src">`nextflow.config:39`</span>.
The hyphenated form does nothing: Nextflow turns
`--genotype-with` into a parameter named `genotypeWith`, and the pipeline reads the default
`reads.csv` instead (checked with Nextflow 26.04.6). Add `-profile cluster` on SLURM, and
`-with-singularity /abs/path/graffite_1.1.0.0.sif` if you pulled the image by hand
([Installation](installation.md)).

`--out out_v1.1` keeps the results apart from the bundled `out/`.

---

## What you should see

When the run completes, `out_v1.1/` holds the files below. With GraffiTE 1.1.0.0,
`1_SV_search/SVs.vcf` has 263 records, `pangenome.vcf` 221, `pangenome.trusted.vcf` 25 and
`GraffiTE.merged.genotypes.vcf.gz` 218. The run took 152 seconds with `--cores 8`.

```
out_v1.1/
├── 1_SV_search/
│   ├── svim-asm_individual_VCFs/HG002.mat.vcf.gz
│   └── SVs.vcf
├── 2_Repeat_Filtering/
│   └── 1/                       one directory per contig, here only chr22
│       ├── genotypes_repmasked_filtered.vcf
│       ├── repeatmasker_dir/
│       └── ultra_out.bed ...
├── 3_TSD_search/
│   ├── pangenome.vcf
│   ├── pangenome.trusted.vcf
│   ├── pangenome.presence-absence.tsv
│   ├── pangenome.presence-absence_trusted.tsv
│   ├── TSD_summary.txt
│   └── TSD_full_log.txt
└── 4_Genotyping/
    ├── pangenie_graph_variants.tsv
    ├── short_test1_genotyping.vcf.gz
    ├── short_test1_genotyping.vcf.gz.tbi
    ├── genotyping_record_audit.tsv
    ├── GraffiTE.merged.genotypes.vcf.gz
    ├── GraffiTE.merged.genotypes.trusted.vcf.gz
    ├── GraffiTE.merged.genotypes.trusted.vcf.gz.tbi
    └── GraffiTE.merged.genotypes.presence-absence_trusted.tsv
```

`pangenome.vcf` is the annotated callset,
`pangenome.trusted.vcf` its conservative subset, and `GraffiTE.merged.genotypes.vcf.gz` the
genotypes of `short_test1`. That last file does not carry every `pangenome.vcf` allele: on the
default `--graph_method pangenie`, `merge_vcfs.py` leaves out alleles that overlap another at the
same site, all but two ALT alleles at one position, and alleles holding a base other than A, C, G
or T. Those alleles get no genotype <span class="src">`bin/pangenie_graph_vcf.py:119-124`</span>.
`pangenie_graph_variants.tsv` says which alleles reached the graph
<span class="src">`module/main.nf:713`</span>, and `genotyping_record_audit.tsv` lists what became
of each one <span class="src">`module/main.nf:907-927`</span>.
`GraffiTE.merged.genotypes.trusted.vcf.gz` is that conservative subset carried through to the
genotypes, written whenever `--human` is off <span class="src">`module/main.nf:929-959`</span>.
The presence-absence TSVs say, per sample, whether the TE is present, whichever way the VCF
record points. Every file is described on [Output files](../reference/outputs.md), and every
INFO field on [VCF fields](../reference/vcf-fields.md).

Nextflow also leaves a `work/` directory in the launch directory with every intermediate
file. It is safe to delete once you have what you need, and it is what `-resume` reads.

---

## Next steps

- [Choosing your inputs](choosing-your-inputs.md) maps your own data to the entry flags.
- [Stage 1: SV search](../guides/discovery.md) explains what happened to the assembly.
- [Resources and scaling](../guides/resources.md) covers whole genomes and clusters.
- [Human MEIs](../guides/human-mei.md) covers the `--human` filter for human data.
