# `GraffiTE`


## Description

`GraffiTE` is a pipeline that finds polymorphic transposable elements in genome assemblies and genotypes the discovered polymorphisms in read sets using a pangenomic approach.

```mermaid
graph TD;
    A-->B;
    A-->C;
    B-->D;
    C-->D;
```

## Installation


### Prerequisites

`GraffiTE` is a `Nextflow` pipeline, with all the dependencies wrapped in a `Singularity` image. It is thus compatible with any Linux system including HPCs.

- install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- install [Singularity](https://apptainer.org/admin-docs/master/installation.html)

### GraffiTE install

To download and run `GraffiTE`, install `nextflow` and use the following command line:

```
nextflow run cgroza/GraffiTE \
   --assemblies assemblies.csv \
   --reads reads.csv \
   --TE_library library.fa \
   --reference reference.fa
```

This will download and cache the `GraffiTE` pipeline and Singularity image for local use. Later runs will skip the slow download step.

## Running GraffiTE

### Parameters

- `--assemblies`: a CSV file that lists the genome assemblies and sample names from which polymorphisms are to be discovered. One assembly per sample and sample names must be unique.

Example `assemblies.csv`:
```
path,sample
/path/to/assembly/sampleA.fa,sampleA_name
/path/to/assembly/sampleB.fa,sampleB_name
/path/to/assembly/sampleZ.fa,sampleZ_name

```

- `--reads`: a CSV file that lists the read sets (FASTQs) and sample names from which polymorphisms are to be genotyped. These samples may be different than the genome assemblies. Only one FASTQ per sample, and sample names must be unique. Paired-end reads must be concatenated in the same file.
> Note that the current genotyper, `PanGenie` is optimized for short-reads. Long-read support will be available soon!

Example `reads.csv`:
```
path,sample
/path/to/reads/sample1.fastq,sample1_name
/path/to/reads/sample2.fastq,sample2_name
/path/to/reads/sampleN.fastq,sampleN_name

```

- `--TE_library`: a FASTA file that lists the consensus sequences of the transposable elements to be discovered. Must be compatible with `RepeatMasker`, i.e. with header in the format: `>TEname#type/subtype` for example `AluY#SINE/Alu`
   - From [DFAM](https://dfam.org/releases/current/families/) (open access): download the latest DFAM release (`Dfam.h5` or `Dfam_curatedonly.h5` files) and use the tool [FamDB](https://github.com/Dfam-consortium/FamDB) to extract the consensus for your model: `famdb.py -i <Dfam.h5> families -f fasta_name -a <taxa> --include-class-in-name > TE_library.fasta`
   - From [Repbase](https://www.girinst.org/server/RepBase/index.php) (paid subscription): use the "RepeatMasker Edition" libraries

- `--reference`: a reference genome of the species being studied. All assemblies are compared to this reference genome.

### Additional parameters

- `--vcf`: a *fully phased* VCF file. Use this if you already have a *phased* VCF file that was produced by GraffiTE, or from a difference source and would like to use the graph genotyping step.
- `--out`: if you would like to change the default output directory (`out/`).
- `--genotype`: true or false. Use this if you would like to discover polymorphisms in assemblies but you would like to skip genotyping polymorphisms from reads.


### Execution profiles
By default, the pipeline will inherit the your `nextflow` configuration and run accordingly.
To execute locally, on SLURM, or AWS, pass one of the `-profile` provided with the `GraffiTE`:
- `standard`
- `cluster`
- `cloud`

For example,
```
nextflow run cgroza/GraffiTE -profile cluster ...
```
will run on SLURM.

## Changing the number of CPUs and memory required by each step
You may alter the following parameters on the command line or in your own `nextflow` configuration file to change how many CPUs and how much memory will be required by each step.

- Step 1, polymorphisms discovery. The memory requirement depends on the genome size of the species. More cores is faster.
```
params.svim_asm_memory
params.svim_asm_threads
```

- Step 2, merging polymorphisms. The requirements depends on the number of assemblies.
```
params.make_vcf_memory
params.make_vcf_threads
```
- Step 3, genotyping polymorphisms from reads. The memory requirements depend on the genome size and size of the read sets. More cores is faster.
```
params.pangenie_memory
params.pangenie_threads
```

The requirements are numbers or strings accepted by `nextflow`. For example, 40 for number of CPUs and '100G' for memory.
