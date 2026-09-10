---
title: Quickstart
description: >-
  Run GraffiTE end to end on the bundled human chromosome 22 test set, and know
  which files to look for when it finishes.
---

# Quickstart

!!! info "Applies to GraffiTE v1.1"
    Verified against `v1.1dev` at commit `4c8e385`. The
    [2024 paper](https://www.nature.com/articles/s41467-024-53294-2) describes v1.0, which
    differs in places; see [v1.0 vs v1.1](v1.0-vs-v1.1.md).

The repository ships a small human dataset that exercises every stage: one haplotype assembly of
chromosome 22, the matching reference, a Dfam library, and a short-read set to genotype.

---

## The test dataset

```bash
tar xzf test/human_test_set.tar.gz
cd GraffiTE_testset
```

| File | Role |
|---|---|
| `HG002.mat.cur.20211005_chr22.fasta.gz` | The maternal HG002 assembly, chromosome 22 only. Stage A input. |
| `hs37d5.chr22.fa` | Reference, chromosome 22 of hs37d5. |
| `human_DFAM3.6.fasta` | TE consensus library for RepeatMasker. |
| `HG002.100k.set1.Illumina.fastq.gz` | 100 k Illumina reads from HG002. Stage C input. |
| `assemblies.csv` | `path,sample` samplesheet naming the assembly `HG002.mat`. |
| `reads.csv` | `path,sample` samplesheet naming the read set `short_test1`. |
| `out/` | Output of a run made before v1.1. See the warning below. |

`reads.csv` has no `type` column. The workflow reads `row.type`, finds nothing, and falls to
the `default` preset <span class="src">`main.nf:190-204`</span>; PanGenie, the default
genotyper, does not use the preset.

!!! warning "The bundled `out/` is not a comparison target"
    That directory was produced by a pre-1.1 release. It contains `svim-asm_variants.vcf`,
    `vcfs.txt` and the OneCode files under `repeatmasker_dir/`, none of which v1.1 writes, and
    it has no `pangenome.trusted.vcf`, no presence-absence TSVs and no ULTRA or polyA fields.
    Use it to see the shape of the older output, not to diff against your run.

---

## Run it

```bash
export NXF_SYNTAX_PARSER=v1
nextflow run cgroza/GraffiTE -r v1.1dev \
  --assemblies assemblies.csv \
  --reference hs37d5.chr22.fa \
  --TE_library human_DFAM3.6.fasta \
  --genotype_with reads.csv \
  --out out_v1.1
```

Pass the four data arguments as shown; the defaults for `--reference` and `--TE_library` are
placeholder filenames that do not exist <span class="src">`nextflow.config:46-47`</span>.
`--genotype_with` is the parameter name as declared <span class="src">`nextflow.config:40`</span>.
The README on `main` writes it with a hyphen, and that form does nothing: Nextflow turns
`--genotype-with` into a parameter named `genotypeWith`, and the pipeline reads the default
`reads.csv` instead (checked with Nextflow 26.04.6). Add `-profile cluster` on SLURM, and
`-with-singularity /abs/path/graffite_latest.sif` if you pulled the image by hand
([Installation](installation.md)).

`--out out_v1.1` keeps the results apart from the bundled `out/`.

---

## What you should see

This page was written from the code, not from a run, so it lists files and not record counts.
When the run completes, `out_v1.1/` holds:

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
    ├── short_test1_genotyping.vcf.gz
    ├── short_test1_genotyping.vcf.gz.tbi
    └── GraffiTE.merged.genotypes.vcf.gz
```

Each directory is one stage. `pangenome.vcf` is the annotated callset,
`pangenome.trusted.vcf` its conservative subset, and `GraffiTE.merged.genotypes.vcf.gz` the
genotypes of `short_test1` at every polymorphism. The presence-absence TSVs say, per sample,
whether the TE is present, whichever way the VCF record points. Every file is described on
[Output files](../reference/outputs.md), and every INFO field on [VCF fields](../reference/vcf-fields.md).

Nextflow also leaves a `work/` directory in the launch directory with every intermediate
file. It is safe to delete once you have what you need, and it is what `-resume` reads.

---

## Next steps

- [Choosing your inputs](choosing-your-inputs.md) maps your own data to the entry flags.
- [Stage A: discovery](../guides/discovery.md) explains what happened to the assembly.
- [Resources and scaling](../guides/resources.md) covers whole genomes and clusters.
- [Human MEIs](../guides/human-mei.md) is the `--human` filter you will want on human data.
