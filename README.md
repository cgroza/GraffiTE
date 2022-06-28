# GraffiTE

GraffiTE is a pipeline that finds polymorphic transposable elements in genome assemblies and genotypes the discovered polymorphisms in read sets using a pangenomic approach.

To download and run `GraffiTE`, install `nextflow` and use the following command line:

```
nextflow run cgroza/GraffiTE --assemblies assemblies.csv --reads reads.csv --TE_library library.fa --reference reference.fa
```

This will download and cache the `GraffiTE` pipeline and Singularity image for local use. Later iterations will run faster.

## Parameters
- `--assemblies`: a CSV file that lists the genome assemblies and sample names from which polymorphisms are to be discovered.

Example `assemblies.csv`:
```
path,sample
/path/to/assembly/sampleA.fa,sampleA_name
/path/to/assembly/sampleB.fa,sampleB_name
/path/to/assembly/sampleZ.fa,sampleZ_name

```

- `--reads`: a CSV file that lists the read sets (FASTQs) and sample names from which polymorphisms are to be genotyped. These samples may be different than the genome assemblies.

Example `reads.csv`:
```
path,sample
/path/to/reads/sample1.fastq,sample1_name
/path/to/reads/sample2.fastq,sample2_name
/path/to/reads/sampleN.fastq,sampleN_name

```

- `--TE_library`: a FASTA file that lists the consensus sequences of the transposable elements to be discovered. Must be compatible with `RepeatMasker`.

-- `--reference`: a reference genome of the species being studied. All assemblies are compared to this reference genome.

## Additional parameters

- `--vcf`: a *fully phased* VCF file. Use this if you already have a *phased* VCF file that was produced by GraffiTE, or from a difference source and would like to use the graph genotyping step.
- `--out`: if you would like to change the default output directory (`out/`).
- `--genotype`: true or false. Use this if you would like to discover polymorphisms in assemblies but you would like to skip genotyping polymorphisms from reads.
