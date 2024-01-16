![](https://i.imgur.com/jvprOAS.png)

[![status](https://img.shields.io/badge/status:-v0.2.5_beta-orange)]() [![status: support](https://img.shields.io/badge/support:-yes-green)]() [![status: preprint](https://img.shields.io/badge/preprint:-BioRxiv:doi.org/10.1101/2023.09.11.557209-red)](https://doi.org/10.1101/2023.09.11.557209)

## Description

`GraffiTE` is a pipeline that finds polymorphic transposable elements in genome assemblies or long read datasets and genotypes the discovered polymorphisms in read sets using a pangenomic approach. `GraffiTE` is developed by **Cristian Groza** and **Clément Goubert** in [Guillaume Bourque's group](https://computationalgenomics.ca/BourqueLab/) at the [Genome Centre of McGill University](https://www.mcgillgenomecentre.ca/) (Montréal, Canada). `GraffiTE` is based on the concept developped in [Groza et al., 2022](https://link.springer.com/protocol/10.1007/978-1-0716-2883-6_5).

1. First, each genome assembly or long read dataset is aligned to the reference genome with [`minimap2`](https://github.com/lh3/minimap2). For each sample considered, structural variants (SVs) are called with [`svim-asm`](https://github.com/eldariont/svim-asm) if using assemblies or [`sniffles2`](https://github.com/fritzsedlazeck/Sniffles) if using long reads and only insertions and deletions relative to the reference genome are kept.
![](https://i.imgur.com/V5NHK3G.png)
2. Candidate SVs (INS and DEL) are scanned with [`RepeatMasker`](https://www.repeatmasker.org/), using a user-provided library of repeats of interest (.fasta). SVs covered ≥80% by repeats are kept. At this step, target site duplications (TSDs) are searched for SVs representing a single TE family.
![](https://i.imgur.com/2qRpojE.png)

3. Each candidate repeat polymorphism is induced in a graph-genome where TEs and repeats are represented as bubbles, allowing reads to be mapped on either presence of absence alleles with [`Pangenie`](https://github.com/eblerjana/pangenie), [`Giraffe`](https://www.science.org/doi/10.1126/science.abg8871) or  [`GraphAligner`](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02157-2).
![](https://i.imgur.com/UyT62yp.png)

🗞️ GraffiTE preprint now on [BioRxiv](https://doi.org/10.1101/2023.09.11.557209)!

----

⚠️ **This is a beta version, with no guarantees! Bug/issues as well as comments and suggestions are welcomed in the [Issue](https://github.com/cgroza/GraffiTE/issues) section of this Github.**

----

## Changelog

**beta 0.2.5 (09-11-23):**
- :beetle: bug fix: fix a VCF annotation issue that was happening when two distinct variants shared the same VCF POS field. Annotations are now distinct depending on the variant sequence.
- cleanup GraphAligner VCF outputs for clarity.


<details><summary>beta 0.2.4 (06-27-23)**:</summary>
<p>

- Refactored `GraffiTE` to use the DSL2 Nextflow syntax.

</p>
</details>

<details><summary>beta 0.2.3 (02-21-22):</summary>
<p>

- :new: feature: You can now perform the initial SV search from both assemblies and long-read together. The variants discovered with each method will be merged together for the filtering and genotyping.
- :new: parameters with defaults added to control time, cpu and memory for each process. This is useful to manage cluster requests when `-profile cluster` is used.
- :beetle: bug fix: merging of variant now only occurs for the same SVTYPE flag (INS or DEL).

</p>
</details>

<details><summary>beta 0.2.2 (02-01-22):</summary>
<p>

- :new: feature: adds `sniffles2` as an alternative to `svim-asm` in order to start SV search from long reads (instead of a genomic assembly).
   - Using the parameter `--longreads` instead of `--assembly` (see inputs) will prompt `GraffiTE` to use `sniffles2`
   - For now, `svim-asm` and `sniffles2` pipeline are separated (either `--longreads` or `--assembly`. We will soon allow to merge the findings of both callers before filtering for repeats.
- :new: feature: adds a divergence preset option to `minimap2` ahead of `svim-asm`. Use the flag `--asm_divergence <asm5/asm10/asm20>`. Defaults is `asm5` (< 5% expected divergence between assembly and reference genome). [See minimap2 documentation](https://lh3.github.io/minimap2/minimap2.html).
- :new: `time`, `cpu` and `memory` directives options added to control the resources needed for each `GraffiTE` process. Useful to optimize scheduler requests while using the `cluster` profile of `GraffiTE`. See details here.

</p>
</details>

<details><summary>beta 0.2.1 (11-30-22 - click to drop-down details):</summary>
<p>

- :new: feature: adds `--RM_vcf` and `--RM_dir` input options. Allows to start a run directly at the TSD search step by providing the VCF and `repeatmasker_dir` produced by the processes `repeatmasker` or `repeatmasker_fromVCF` (found in the output folder `2_Repeat_Filtering`). This is useful if a run crashed during any of the TSD search processes and the job is not recoverable by Nextflow. Providing `--RM_vcf` and `--RM_dir` will bypass SV calling with `minimap2/svim_asm` (`svim_asm` process) and `repeatmasker/repeatmasker_fromVCF` processes.
- :beetle: bug fix: TSD search is now performed by batches of 100 variants, which will reduce by a factor 100 the number of temporary working directories (which can cause storage to run over inodes' quota). If more than 100 variants are present, TSDs will be searched in parallel batches (up to the number of available CPUs).

</p>
</details>

<details><summary>beta 0.2 (11-11-22 - click to drop-down details):</summary>
<p>

- :new: feature:  adds two new read aligners: [`giraffe`](https://github.com/vgteam/vg#mapping) (short read optimized, works also with long-reads) and [`graphAligner`](https://github.com/maickrau/GraphAligner) (long-read, error-prone compliant). 
   - usage: `--graph_method [pangenie/giraffe/graphaligner]` default: `pangenie` (short accurate reads)
- :new: feature:  adds `--vcf` input option: requires a sequence resolved (REF and ALT allele sequences in VCF). Will bypass genome alignments and proceed with repeat annotations, TSD search, and reads mapping (optional).
- :new: feature:  adds `--graffite_vcf` input option: requires a VCF created by `GraffiTE` (in the outputs `3_TSD_search/pangemome.vcf`). Will skip all steps but read mapping.
- :beetle: bug fix: remove the dependency to `biomartr`

</p>
</details>



<details><summary>beta 0.1 (11-02-22 - click to drop-down details):</summary>
<p>

- first release

</p>
</details>

*It is required to update both the repository (`git pull`) and image to see changes*

----

## Workflow

![](https://i.imgur.com/X0jOkVn.png)

## Installation

### Prerequisites

`GraffiTE` is a `Nextflow` pipeline, with all the dependencies wrapped in an `Apptainer` image. It is thus compatible with any Linux system including HPCs.

- install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- install [Apptainer](https://apptainer.org/admin-docs/master/installation.html)

### GraffiTE install

- If an internet connection is accessible from the compute nodes, the general command shown in the next section will download and cache the `GraffiTE` pipeline and Apptainer image for local use. Later runs will skip the slow download step. It is however required to add the repository to the apptainer list, by typing:

```
apptainer remote add --no-login SylabsCloud cloud.sycloud.io
apptainer remote use SylabsCloud
```

- Alternatively, this repository can be cloned and the apptainer image downloaded at a specific location:
   - 1. Clone the Github repository
   ```
   git clone https://github.com/cgroza/GraffiTE.git
   ```
   - 2. Pull the apptainer image (this is long but only required once)
   ```
   apptainer remote add --no-login SylabsCloud cloud.sycloud.io
   apptainer remote use SylabsCloud
   apptainer pull --arch amd64 graffite_latest.sif library://cgroza/collection/graffite:latest
   ```
   - 3. Override the default image path in the file `nextflow.config` from `library://cgroza/collection/graffite:latest` to `<your-path>/graffite_latest.sif`. Alternatively, the `Nextflow` command `-with-singularity <your-path>/graffite_latest.sif` can be used when running `GraffiTE` (it will override the presets in `nextflow.config`).

## Important note

We are aware of a common issue araising when the pipeline call a temporary directory (/tmp). The most common symptom is that though the program may complete without error, it skips over "tsd_search" and "tsd_report". The program wont produce a vcf file (`3_TSD_search/pangenome.vcf`) and the vcf in `2_Repeat_Filter` has no variants. While we will try to fix this in a next update, an easy fix is to ammend the `nextflow.config` file as follow. 

1. Locate the file:
   - Either in `~/.nextflow/assets/cgroza/GraffiTE/nextflow.config`
   - or in the cloned GitHub repository.

2. Ammend the file:

replace:
```
singularity.runOptions = '--contain'
```

with

```
singularity.runOptions = '--contain -B <path-to-writable-dir>/:/tmp'
```
> replace `<path-to-writable-dir>` with any writable path on your host machine

## Running GraffiTE

- The general command to run `GraffiTE` is as follow:

```
nextflow run cgroza/GraffiTE \
   --assemblies assemblies.csv \
   --TE_library library.fa \
   --reference reference.fa \
   --graph_method pangenie \
   --reads reads.csv
```

- If using from a local apptainer image and with a clone of the Github repository:

```
nextflow run <path-to-install>/GraffiTE/main.nf \
   --assemblies assemblies.csv \
   --TE_library library.fa \
   --reference reference.fa \
   --reads reads.csv [-with-singularity <your-path>/graffite_latest.sif]
```
> As a `Nextflow` pipeline, commad line arguments for `GraffiTE` can be distinguished between pipeline-related commands, prefixed with `--` such as `--reference` and `Nextflow`-specific commands, prefixed with `-` such as `-resume` (see [`Nextflow` documentation](https://www.nextflow.io/docs/latest/index.html)).

A small test set is included in the `test/human_test_set.tar.gz` file.
Download and decompress the file and run:
```
nextflow run https://github.com/cgroza/GraffiTE --reference hs37d5.chr22.fa --assemblies assemblies.csv --reads reads.csv --TE_library human_DFAM3.6.fasta
```
This will show a complete run of the GraffiTE pipeline, with the output stored in `out`.

## Parameters

#### Input files

- `--assemblies`: a CSV file that lists the genome assemblies and sample names from which polymorphisms are to be discovered. One assembly per sample and sample names must be unique. **The header is required**. 

   Example `assemblies.csv`:
   ```
   path,sample
   /path/to/assembly/sampleA.fa,sampleA_name
   /path/to/assembly/sampleB.fa,sampleB_name
   /path/to/assembly/sampleZ.fa,sampleZ_name
   ```

AND/OR

- `--longreads`: a CSV file that lists the longreads FASTQ, sample names, and type of longreads (hifi/pb/ont) from which polymorphisms are to be discovered. One FASTQ per sample and sample names must be unique. **The header is required**.

   Example `longreads.csv`:
   ```
   path,sample,type
   /path/to/reads/sampleA.fq.gz,sampleA_name,pb
   /path/to/reads/sampleB.fq.gz,sampleB_name,hifi
   /path/to/reads/sampleZ.fq.gz,sampleZ_name,ont
   ```

AND (always required)

- `--TE_library`: a FASTA file that lists the consensus sequences (models) of the transposable elements to be discovered. Must be compatible with `RepeatMasker`, i.e. with header in the format: `>TEname#type/subtype` for example `AluY#SINE/Alu`. The library can include a single repeat model or all the known repeat models of your species of interest.
   - From [DFAM](https://dfam.org/releases/current/families/) (open access): download the latest DFAM release (`Dfam.h5` or `Dfam_curatedonly.h5` files) and use the tool [FamDB](https://github.com/Dfam-consortium/FamDB) to extract the consensus for your model: `famdb.py -i <Dfam.h5> families -f fasta_name -a <taxa> --include-class-in-name > TE_library.fasta`
   - From [Repbase](https://www.girinst.org/server/RepBase/index.php) (paid subscription): use the "RepeatMasker Edition" libraries

- `--reference`: a reference genome of the species being studied. All assemblies or long-reads  in input are compared to this reference genome.

- `--graph_method`: can be `pangenie`, `giraffe` or `graphaligner`, select which graph method will be used to genotyped TEs. Default is `pangenie` and it is optimized for short-reads. `giraffe` can handle both short and long reads, and `graphaligner` is optimized for long reads. 
>Note that both `giraffe` and `graphaligner` will spawn a process called `graphAlignReads`, while `pangenie` will spawn a process called `pangenie`.

- `--reads`: a CSV file that lists the read sets (FASTQ/FASTQ.GZ) and sample names from which polymorphisms are to be genotyped. These samples may be different than the genome assemblies. **The header is required**. Only one FASTQ/FASTQ.GZ per sample, and sample names must be unique. **Paired-end reads must be concatenated into a single file (`Pangenie`)**. In case `--longreads` is used as input, the same table can be used for `--longreads` and `--reads` (but not the opposite: `type` column is needed in `--longreads`, optional for `--reads`).

   Example `reads.csv`:
   ```
   path,sample
   /path/to/reads/sample1.fastq,sample1_name
   /path/to/reads/sample2.fastq,sample2_name
   /path/to/reads/sampleN.fastq,sampleN_name
   ```
   or
   ```
   path,sample,type
   /path/to/reads/sampleA.fq.gz,sampleA_name,pb
   /path/to/reads/sampleB.fq.gz,sampleB_name,hifi
   /path/to/reads/sampleZ.fq.gz,sampleZ_name,ont
   ```

#### Additional parameters

- `--out`: if you would like to change the default output directory (`out/`).
- `--genotype`: true or false. Use this if you would like to discover polymorphisms in assemblies but you would like to skip genotyping polymorphisms from reads.
- `--tsd_win`: the length (in bp) of flanking region (5' and 3' ends) for Target Site Duplication (TSD) search. Default 30bp. By default, 30bp upstream and downstream each variant will be added to search for TSD. (see also [TSD section](#tsd-module))
- `--cores`: global CPU parameter. Will apply the chosen integer to all multi-threaded processes. See [here](#changing-the-number-of-cpus-and-memory-required-by-each-step) for more customization.
- `--mammal`: Apply mammal-specific annotation filters (see [Mammal filter section](#mammalian-filters---mammal) for more details). 
   - (i) will search for LINE1 5' inversion (due to Twin Priming or similar mechanisms). Will call 5' inversion if (and only if) the variant has two RepeatMasker hits on the same L1 model (for example L1HS, L1HS) with the same hit ID, and a `C,+` strand pattern. 
   - (ii) will search for VNTR polymorphism between orthologous SVA elements.

#### Pipeline Shortcuts

These parameters can be used to bypass different steps of the pipeline.

- `--vcf`: a *sequence resolved* VCF containing both REF and ALT variants sequences. This option will bypass the SV discovery and will proceed to annotate and filter the input VCF for repeats and TSD, as well as genoyping (unless `--genotype false` is set)
- `--RM_vcf`+`--RM_dir`: bypasses SV discovery and filtering (RepeatMasker) and starts at the TSD search process. `--RM_vcf` can be found in the outputs: `2_Repeat_Filtering/genotypes_repmasked_filtered.vcf` and `--RM_dir` in `2_Repeat_Filtering/repeatmasker_dir`
- `--graffite_vcf`: Use this if you already have a VCF file that was produced by GraffiTE (see output: `3_TSD_Search/pangenome.vcf`), or from a difference source and would like to use the graph genotyping step. The file must be a [fully-phased](https://github.com/eblerjana/pangenie#input-variants) VCF. Note that TE annotation won't be performed on this file (see `--vcf` instead), and only genotyping will be performed.

#### Process-specific parameters

##### SV detection with `svim-asm` (from assemblies)

- `--map_asm_threads`: number of cores allocated to `minimap2` or `winnowmap` in the `map_asm` process. Overrides `--cores`.
- `--map_asm_memory`: RAM limit for the `minimap2` or `winnowmap` in the `map_asm` process. Default is unset.
- `--map_asm_time`: for `cluster` profile, max time for the scheduler for the `map_asm` process. Default is 1h.
- `--svim_asm_threads`: number of cores allocated to `svim_asm` process. Overrides `--cores`.
- `--svim_asm_memory`: RAM limit for the SV search in `svim_asm` process. Default is unset.
- `--svim_asm_time`: for `cluster` profile, max time for the scheduler for this process. Default is 1h.
- `--asm_divergence`: divergence preset option for `minimap2` ahead of `svim-asm`. Use the flag . `asm5`/`asm10`/`asm20`  Defaults is `asm5` (< 5% expected divergence between assembly and reference genome). [See minimap2 documentation](https://lh3.github.io/minimap2/minimap2.html).


- `--mini_K`: `minimap2` parameter `-K`. *Number of bases loaded into memory to process in a mini-batch. Similar to option -I, K/M/G/k/m/g suffix is accepted. A large NUM helps load balancing in the multi-threading mode, at the cost of increased memory.* Default 500M.
- `--stSort_m`: `samtools sort` parameter `-m` (for each alternative assembly, post-`minimap2`): *Approximately the maximum required memory per thread, specified either in bytes or with a K, M, or G suffix.* Default in `GraffiTE` is 4G.
- `--stSort_t`: `samtools sort` parameter `@` (for each alternative assembly, post-`minimap2`): *Set number of sorting and compression threads.* Default in `GraffiTE` is 4 threads. 

##### SV detection with `sniffles2` (from long reads)

- `--map_longreads_threads`: number of cores allocated to `minimap2` or `winnowmap` in the `map_longreads` process. Overrides `--cores`.
- `--map_longreads_memory`: RAM limit for the `minimap2` or `winnowmap` in the `map_longreads` process. Default is unset.
- `--map_longreads_time`: for `cluster` profile, max time for the scheduler for the `map_longreads` process. Default is 1h.
-`--sniffles_threads`:  number of `minimap2` threads (parameter `-t` in `minimap2`). Overrides `--cores`.
-`--sniffles_memory`: RAM limit for the SV search (`minimap2`+`sniffles2`) process. Default is unset.
-`--sniffles_time`: for `cluster` profile, max time for the scheduler for this process. Default is 2h.
- `--stSort_m`: `samtools sort` parameter `-m` (for each long-read alignment, post-`minimap2`): *Approximately the maximum required memory per thread, specified either in bytes or with a K, M, or G suffix.* Default in `GraffiTE` is 4G.
- `--stSort_t`: `samtools sort` parameter `@` (for each long-read alignment, post-`minimap2`): *Set number of sorting and compression threads.* Default in `GraffiTE` is 4 threads. 

##### SV Annotation (RepeatMasker)
- `--repeatmasker_threads`: number of RepeatMasker threads. Overrides `--cores`.
- `--repeatmasker_memory`: RAM limit for the RepeatMasker (annotation) process. Default is unset.
- `--repeatmasker_time`: for `cluster` profile, max time for the scheduler for this process. Default is 2h.

##### Genotyping with Pangenie
- `--pangenie_threads`: number of `Pangenie` threads. Overrides `--cores`.
- `--pangenie_memory`: RAM limit for the Pangenie (genotyping) process. Default is unset.
- `--pangenie_time`: for `cluster` profile, max time for the scheduler for this process. Default is 2h.

##### Genotyping with Giraffe, GraphAligner and `vg call`
- `--make_graph_threads`: threads for creating the graph with `vg autoindex` (Giraffe) or `vg construct` (GraphAligner). Default is 1.
- `--make_graph_memory`: RAM limit for creating the graph with `vg autoindex` (Giraffe) or `vg construct` (GraphAligner). Default is unset.

- `--graph_align_theads`: threads for aligning reads to the graph with `vg giraffe` or `GraphAligner`. Default is 1.
- `--graph_align_memory`: RAM limit for aligning reads to the graph with `vg giraffe` or `GraphAligner`. Default is unset.
- `--graph_align_time`: for `cluster` profile, max time for the scheduler for this process. Default is 12h.

- `--vg_call_threads`: threads for calling genotypes with `vg call` on graph alignments. Default is 1.
- `--vg_call_memory`: RAM limit for calling genotypes with `vg call` on graph alignments. Default is unset.

- `--min_mapq`: Minimum mapping quality to consider when counting read depth on nodes. Default is 0.
- `--min_support`: Minimum required read depth on `allele,bubble` to consider for genotyping. The first number is the minimum read depth on allele, and the second is the minimum depth on the entire bubble/locus. Default is `2,4`.

#### `Nextflow` parameters

`Nextflow`-specific parameters can be passed in addition to those presented above. These parameters can be distinguished by the use of a single `-`, such as `-resume`. See `Nextflow` documentation for more details.

- `-resume`: if nothing is changed in the command line and the `/work` folder created by `Nextflow`, the pipeline will resume after the last chached process.
- `-with-singularity`: if a local apptainer image is used, this parameter will override the default image path given in `nextflow.config`.
- `-with-report report.html`: for a Nextflow report on resource usage to help tune the CPU and memory parameters for your genome/species.

## Outputs

The results of `GraffiTE` will be produced in a designated folder with the option `--out`. The output folder contains up to 4 sub-folders (3 if `--genotype false` is set). Below is an example of the output folder using two alternative assemblies of the human chromosome 1 (maternal and paternal haplotypes of HG002) and two read-sets from HG002 for genotyping.

```
OUTPUT_FOLDER/
├── 1_SV_search
│   ├── HG002_mat.vcf
│   └── HG002_pat.vcf
├── 2_Repeat_Filtering
│   ├── genotypes_repmasked_filtered.vcf
│   └── repeatmasker_dir
│       ├── ALL.onecode.elem_sorted.bak
│       ├── indels.fa.cat.gz
│       ├── indels.fa.masked
│       ├── indels.fa.onecode.out
│       ├── indels.fa.out
│       ├── indels.fa.out.length
│       ├── indels.fa.out.log.txt
│       ├── indels.fa.tbl
│       ├── onecode.log
│       └── OneCode_LTR.dic
├── 3_TSD_search
│   ├── pangenome.vcf
│   ├── TSD_full_log.txt
│   └── TSD_summary.txt
└── 4_Genotyping
    ├── GraffiTE.merged.genotypes.vcf
    ├── HG002_s1_10X_genotyping.vcf.gz
    ├── HG002_s1_10X_genotyping.vcf.gz.tbi
    ├── HG002_s2_10X_genotyping.vcf.gz
    └── HG002_s2_10X_genotyping.vcf.gz.tbi
```

- `1_SV_search/`
   - This folder will contain 1 VCF file per alternative assembly. The format is `[assembly_name].vcf` with `[assembly_name]` as set in the file `assemblies.csv`
- `2_Repeat_Filtering/`
   - `genotypes_repmasked_filtered.vcf` a vcf file with the merged variants detected in each alternative assembly. The merge is made with [`SURVIVOR`](https://github.com/fritzsedlazeck/SURVIVOR) with the parameters `SURVIVOR merge vcfs.txt 0.1 0 0 0 0 100`. Details about the vcf annotation can be found in the [VCF section](#output-vcfs) of the manual. This VCF contains only variants for witch repeats in the `--TE_library` file span more than 80% of the sequence (from 1 or more repeat models).
   - `repeatmasker_dir/`:
      - `indels.fa.*`: `RepeatMasker` output files. `indels.fa` represents all SV sequences queried to `RepeatMasker`. See the [RepeatMasker documentation](https://www.repeatmasker.org/webrepeatmaskerhelp.html) for more information. 
      - `ALL.onecode.elem_sorted.bak`: original `OneCodeToFindThemAll` outputs. see [here](https://mobilednajournal.biomedcentral.com/articles/10.1186/1759-8753-5-13) fore more details.
      - `OneCode_LTR.dic`: `OneCodeToFindThemAll` LTR dictionary automatically produced from `--TE_library` see [here](https://mobilednajournal.biomedcentral.com/articles/10.1186/1759-8753-5-13) fore more details.
      - `onecode.log`: log file for `OneCodeToFindThemAll` process.
- `3_TSD_Search/` (see [TSD section](#tsd-module))
   - `pangenome.vcf` final VCF containing all retained repeat variants and annotations (with TSD if passing the TSD filters). This file is used later by `Pangenie`,`Giraffe` or `graphAligner` to create the genome-graph onto which reads are mapped for genotyping. (example [here](#output-vcfs)). Can be re-used for genotyping only with `--graffite_vcf pangenome.vcf`
   - `TSD_summary.txt`: tab delimited output of the TSD search module. 1 line per variant. See [TSD section](#tsd-module) for more information. "PASS" entries are reported in the `pangenie.vcf` and final (with genotypes) VCF.
   - `TSD_full_log.txt:`detailed (verbose rich) report of TSD search for each SV (see [TSD section](#tsd-module)).
- `4_Genotyping/`
   - `GraffiTE.merged.genotypes.vcf`: final mutli-sample VCF with the genotypes for each sample present in the `--reads` file. See [VCF section](#output-vcfs) for more details.
   - `*.vcf.gz` individual genotypes (do not contain TE annotation)
   - `*.vcf.gz.tbi` index for individual VCFs.

> Note that intermediate files will be written in the `./work` folder created by `Nextflow`. Each `Nextflow` process is run in a separate working directory. If an error occurs, `Nextflow` will points to the specific working directory. Moreover, it is possible to resume interrupted jobs if the `./work` folder is intact and you use the same command, plus the `-resume` (1 single `-`) tag after your command. It is recommended to delete the `./work` folder regularly to avoid storage issues (more than space, it can aggregate a LOT of files through time). More info about `Nextflow` usage can be found [here](https://www.nextflow.io/docs/latest/index.html).

#### Output VCFs

`GraffiTE` outputs variants in the [VCF 4.2 format](https://samtools.github.io/hts-specs/VCFv4.2.pdf). Additional fields are added in the INFO column of the VCF to annotate SVs containing TEs and other repeats (`3_TSD_Search/pangenie.vcf` [do not contain individual genotypes, only the list of variants] and `4_Genotyping/GraffiTE.merged.genotypes.vcf` which contains a genotype column for each reads-set).

- `3_TSD_Search/pangenie.vcf`
```
1       8501990 HG002_mat.svim_asm.INS.94       T       TCAATACACACACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGCCGGACTGCGGACTGCAGTGGCGCAATCTCGGCTCACTGCAAGCTCCGCTTCCCGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCCAGTAGCTGGGACTACAGGCGCCCGCCACCGCGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTCGATCTCCTGACCTCATGATCCACCCGCCTCGGCCTCCCAAAGTGCTGGGACTACAGGCGTGAGCCACCGCGCCCGGC        .       PASS    SUPP=1;SUPP_VEC=10;SVLEN=345;SVTYPE=INS;SVMETHOD=SURVIVOR1.0.7;CHR2=1;END=8501990;CIPOS=0,0;CIEND=0,0;STRANDS=+-;n_hits=1;fragmts=1;match_lengths=316;repeat_ids=AluYb9;matching_classes=SINE/Alu;RM_hit_strands=C;RM_hit_IDs=15016;total_match_length=316;total_match_span=0.913295;mam_filter_1=None;mam_filter_2=None;TSD=AATACACACACTTTTT,AATACACACACTTTTT    GT  1|0
```
> An example of AluYb9 insertion relative to the reference genome (hg19 was used for this example). The genotype is always heterozygous in order to create both allele in the graph used for genotyping

- `4_Genotyping/GraffiTE.merged.genotypes.vcf`
```
1       8501990 HG002_mat.svim_asm.INS.94       T       TCAATACACACACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGCCGGACTGCGGACTGCAGTGGCGCAATCTCGGCTCACTGCAAGCTCCGCTTCCCGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCCAGTAGCTGGGACTACAGGCGCCCGCCACCGCGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCGTGTTAGCCAGGATGGTCTCGATCTCCTGACCTCATGATCCACCCGCCTCGGCCTCCCAAAGTGCTGGGACTACAGGCGTGAGCCACCGCGCCCGGC        .       PASS    UK=51;MA=0;AF=0.5;AK=13,38;CIEND=0,0;CIPOS=0,0;CHR2=1;END=8501990;SVLEN=345;SVMETHOD=SURVIVOR1.0.7;SVTYPE=INS;SUPP_VEC=10;SUPP=1;STRANDS=+-;n_hits=1;match_lengths=316;repeat_ids=AluYb9;matching_classes=SINE/Alu;fragmts=1;RM_hit_strands=C;RM_hit_IDs=15016;total_match_length=316;total_match_span=0.913295;mam_filter_1=None;mam_filter_2=None;TSD=AATACACACACTTTTT,AATACACACACTTTTT      GT:GQ:GL:KC     0/1:10000:-81.8909,0,-64.99:7   0/1:10000:-81.8909,0,-64.99:7
```
> An example of AluYb9 insertion relative to the reference genome (hg19 was used for this example). Genotypes are based on read mapping for each individual.

```
1  33108378 HG002_pat.svim_asm.INS.206 T  TTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCACCAGACTGGAGTACAATGGCACAATCTCGGCTTACTGCAACTTCCGCCTCCTGGGTTCAAGCAATTCCCCTGCCTCAGCCTCCTGAGTAGCTGGGATTACAGACGTGTGCCACCATGCCTGGCTAATTTTTTGTATTTTA
GCAGAGACGGAGTTTCACCATGTTGGCCAGGATGCTCTCAATCTCCTTACCTCATGATCCGCCAGCCTCGGCCTCCCAAAGTGCTGGGATTATTACAGGCATGAGCCACAGTCCCAGGTCTTTAGACAAACTCAACCCATTATCAATCAAAAAATGTTTAAATTCACTTATAGCATGGAAGCTACCCCACCCCTCCCCCCTCCCCCCTCCCGCCCCCCCCAGCTTTGAGTTGTCCCACCTTTCTGGACCAAAGCA ATGTATTTCTTAAACTTAATTGATTAATGTCTCATGCCTCTCTGAAATGTATAAAACCAAACTGTGCCCTGACCACCTTGGGCACACTGAGCACATGTTCTCAGGATCTCCAGAGGGCTGTGTCAGGGGCCATGGTCACATTTGGCTCAGAATACATCTCTTCAAATATTTTATAGAGTTCGACTATTTTGTCAACAATTAAAAAGGCACCTATTCAGAAT
ATTAAAAGTTAAGATTTAATAACATCAACAGTTCTTACTGATTCATCAAATATTTTTTTTTTTGAGACCGAGTCTCGCTCTATCGCCCAGGCTGGAGGGCAGTGGCACAATCTCTGTTCACTGCAACCTCCGCCTCCCGGGTTCAAGCGATTCTCCTGCCTCAGCCTCCCGAATAGCTGGGACTACATGCGCGTGCCACCACGCCTGGCTAATTTTTGTATTTTTAGTAGAGACGGAGTTTCACAACGTTGGCCAGGATGGTCTCGATCCCTTGACCTCATGATCCGCCTGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGTGTGAGCCACCGGCGCCTGGCCAAAACAAAA  .PASS K=301;MA=0;AF=0.5;AK=2,299;CIEND=0,1;CIPOS=0,0;CHR2=1;END=33108378;SVLEN=1002;SVMETHOD=SURVIVOR1.0.7;SVTYPE=INS;SUPP_VEC=11;SUPP=2;STRANDS=+-;n_hits=4;match_lengths=293,331,80,291;repeat_ids=AluSc8,MER4E1,Charlie1a,AluSc;
matching_classes=SINE/Alu,LTR/ERV1,DNA/hAT-Charlie,SINE/Alu;fragmts=1,1,1,1;RM_hit_strands=C,+,C,C;RM_hit_IDs=28269,28270,28271,28272;total_match_length=991;total_match_span=0.988036;mam_filter_1=None;mam_filter_2=None   GT:GQ:GL:KC 1/1:10000:-450.343,-147.4,0:4 1/1:10000:-450.343,-147.4,0:4
```
> A more complex example with `n_hit=4`

VCF column:
- `(1) CHROM`: chromosome/scaffold/contig 
- `(2) POS`: position (in bp) of the SV start, relative to the reference genome
- `(3) ID`: variant name
- `(4) REF`: reference allele
- `(5) ALT`: alternative allele
- `(6) QUAL`: not used
- `(7) FILTER`: currently not used. "PASS" is used by default but does not inform about variant quality (for now!)
- `(8) INFO`:
   - `UK` (`4_Genotyping/GraffiTE.merged.genotypes.vcf` only): [`Pangenie`] Total number of unique kmers
   - `MA` (`4_Genotyping/GraffiTE.merged.genotypes.vcf` only): [`Pangenie`] Number of alleles missing in panel haplotypes
   - `AF` (`4_Genotyping/GraffiTE.merged.genotypes.vcf` only): [`Pangenie`] Allele Frequency
   - `AK` (`4_Genotyping/GraffiTE.merged.genotypes.vcf` only): [`Pangenie`] Number of unique kmers per allele. Will be -1 for alleles not covered by any input haplotype path
   - `CIEND` (ignore)
   - `CIPOS` (ignore)
   - `CHR2` (ignore)
   - `END`: End position of the SV on the reference genome
   - `SVLEN`: Length of the SV (bp), can be negative
   - `SVMETHOD=SURVIVOR1.0.7;` (ignore)
   - `SVTYPE`: Type of SV (can be INS or DEL)
   - `SUPP_VEC`: Support Vector from SURVIVOR (merge of individual loci). SUPP_VEC=01 means two alternative assemblies were used, the SV is absent from the first one and present in the second one.
   - `SUPP`: Number of assemblies with the variant
   - `STRANDS=+-;` (ignore)
   - `n_hits`: number of distinct RepeatMasker hits on the SV
   - `match_lengths`: length of each RepeatMasker hit. If `n_hits` > 1, lengths of each hit are comma separated
   - `repeat_ids`: target name of each RepeatMasker hit. If `n_hits` > 1, names for each hit are comma separated
   - `matching_classes`: classification of each RepeatMasker hit. If `n_hits` > 1, classification for each hit are comma separated
   - `fragmts`: number of fragments stitched together for each RepeatMasker hit. If `n_hits` > 1, the number of stitched fragments for each hit are comma separated
   - `RM_hit_strands`: strands for each RepeatMasker hit. If `n_hits` > 1, the strands of each hit are comma separated. Can be `+` or `C` (complement)
   - `RM_hit_IDs`: unique RepeatMasker hit ID (last column of the `.out` file of repeatmasker). If `n_hits` > 1, hit IDs are comma separated. Fragments stitched with `OneCodeToFindThemAll` are shown separated with `/`.
   - `total_match_length`: total number of bp covered by repeats in the SV
   - `total_match_span`: proportion of the SV covered by repeats (minimum is 0.8)
   - `mam_filter_1`: `5P_INV` will be shown if the SV is a LINE1 with a 5' inversion; Null otherwise; (only present if `--mammal` is set)
   - `mam_filter_2`: `SVA_VNTR` if the SV is a length polymorphism of the VNTR region of an SVA element; Null otherwise; (only present if `--mammal` is set)
   - `TSD`: Target Site Duplication (left_TSD,right_TSD); only present if TSD passes filters (see TSD section)
- `(9) FORMAT` and `(10) GENOTYPE`
   - `GT`: Genotype (0=reference allele, 1=alternative allele, .=missing)
   - `GQ`: (`4_Genotyping/GraffiTE.merged.genotypes.vcf` only): [`Pangenie`] Genotype quality: phred scaled probability that the genotype is wrong.
   - `GL`: (`4_Genotyping/GraffiTE.merged.genotypes.vcf` only): [`Pangenie`] Comma-separated log10-scaled genotype likelihoods for absent, heterozygous, homozygous.
   - `KC`: (`4_Genotyping/GraffiTE.merged.genotypes.vcf` only): [`Pangenie`] Local kmer coverage.

When using `Giraffe` and `GraphAligner` with `vg call`, the following fields are also present:
- `AT`: Allele traversal as path in graph
- `DP`: Total Depth
- `AD`: Allelic depths for the ref and alt alleles in the order listed">
- `MAD`: Minimum site allele depth
- `GL`: Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy
- `GQ`: Genotype Quality, the Phred-scaled probability estimate of the called genotype
- `GP`: Genotype Probability, the log-scaled posterior probability of the called genotype
- `XD`: eXpected Depth, background coverage as used for the Poisson model

## TSD module

For SVs with a single TE insertion detected (`n_hits=1`, and LINE1s with the flag `mam_filter_1=5P_INV`) target site duplication are searched by comparing the flanking regions following this workflow:

- 1 extract the flanking sequences of each filtered SV: 
   - 1.1 extract the bases not identified as repeat by RepeatMasker in the 5' and 3' end of the SV (these regions will often include one TSD, or a partial sequence of the TSD)
   - 1.2 extract an additional (by default 30) bp on each side of the SV from the reference genome.
- 2 perform the TSD search:
   - Combine the extracted flanking and create the L (5') and R (3') fragments for each SV. 
   - If present, trim 5' poly-A or 3' poly-T (leaves only 3 As or Ts) before alignments but keep track of the poly-A/T length.
   - Call `blastn` to align with a seed of 4 bp 
   - Applies PASS filters and return summary files. PASS is currently given if:
      - L and R flanks match within +/- 5 bp of the TE ends (as defined by `RepeatMasker`, "Ns" nucleotides)
      - tolerate (TE hit divergence to consensus x alignment length) mismatches+gaps or 1 mismatch+gap if (TE hit divergence to consensus x alignment length) < 1
      - tolerate offset of +/- poly-A/T length

![](https://i.imgur.com/ZzO1ZcQ.png)

The script also account for the presence of poly-A/T

![](https://i.imgur.com/ejDKo5x.png)

- `TSD_summary.txt` output file (The header is not present in the real file).
   ```
   SV_name                          RM_family_name    RM_hit_strand  RM_hit_divergence TSD_length  Mismatches  Gaps    5P_TSD_end   5P_offset      3P_TSD_start    3P_offset     5P_TSD            3P_TSD            FILTER
   HG002_mat.svim_asm.DEL.1014      AluY              C              2.2               10          0           0       -1           0              1               0             ATTATTATTA        ATTATTATTA        PASS
   HG002_mat.svim_asm.DEL.1013      L1HS              C              1.3               16          0           0       -15          3              1               0             AGTATTCTGGATTTTT  AGTATTCTGGATTTTT  FAIL
   G002_mat.svim_asm.DEL.1015       L1HS              +              1.0738            4           0           0       -9           0              1               0             AAAG              AAAG              FAIL
   HG002_mat.svim_asm.DEL.102       AluYa5            C              0.3               11          0           0       -1           0              1               0             CTGCATACTTT       CTGCATACTTT       PASS
   HG002_mat.svim_asm.DEL.1011      L1P2              C              6.9               4           0           0       -21          0              1               0             CATC              CATC              FAIL
   HG002_mat.svim_asm.DEL.1005      AluY              C              1.0               12          0           0       -1           0              1               0             CCAGAAGTCTTT      CCAGAAGTCTTT      PASS
   HG002_mat.svim_asm.DEL.1010      AluYh3            +              2.4               12          0           0       -1           0              1               0             AATTTCTATCTC      AATTTCTATCTC      PASS
   ```
 - `TSD_full_log.txt:`detailed (verbose rich) report of TSD search for each SV.
   ```
      --- TSD search for HG002_mat.svim_asm.DEL.1014 ---

   >L|5P_end
   ACAGGCGTGAGCCTCCACGCCTGGCCTAGATATTATTATTATTATTATTA
   ||||||||||||||||||||||||||||||||||||||||||||||||||
   1   5    10   15   20   25   30   35   40   45   50
   >R|3P_end
   ATTATTATTAACCTATTTTACAGATGAGGG
   ||||||||||||||||||||||||||||||||||||||||||||||||||
   1   5    10   15   20   25   30   35   40   45   50

   3' poly_A: element is in C orientation, will not search for poly_A
   5' poly_T: 0 bp, will not remove anything for alignment


   Building a new DB, current time: 11/02/2022 22:27:12
   New DB name:   /scratch/cgoubert/GraffiTE/work/d1/3d8805a29e13fad52ed5aa1e7a9e76/L.short.fasta
   New DB title:  L.short.fasta
   Sequence type: Nucleotide
   Keep MBits: T
   Maximum file size: 1000000000B
   Adding sequences from FASTA; added 1 sequences in 0.000507116 seconds.

   candidate hits from blastn:
   R|3P_end        L|5P_end        100.000 10      0       0       1       10      41      50      0.001   19.6
   R|3P_end        L|5P_end        100.000 4       0       0       1       4       47      50      3.1      8.5
   R|3P_end        L|5P_end        100.000 10      0       0       1       10      38      47      0.001   19.6
   R|3P_end        L|5P_end        100.000 10      0       0       1       10      35      44      0.001   19.6
   R|3P_end        L|5P_end        100.000 10      0       0       1       10      32      41      0.001   19.6
   R|3P_end        L|5P_end        100.000 8       0       0       3       10      31      38      0.018   15.9
   R|3P_end        L|5P_end        100.000 4       0       0       12      15      25      28      3.1      8.5
   R|3P_end        L|5P_end        87.500  8       0       1       14      20      37      44      3.1      8.5
   R|3P_end        L|5P_end        87.500  8       0       1       14      20      31      38      3.1      8.5
   R|3P_end        L|5P_end        100.000 4       0       0       20      23      1       4       3.1      8.5
   R|3P_end        L|5P_end        100.000 4       0       0       22      25      28      31      3.1      8.5
   R|3P_end        L|5P_end        100.000 4       0       0       25      28      8       11      3.1      8.5

   candidate TSDs:
   ACAGGCGTGAGCCTCCACGCCTGGCCTAGATATTATTATTATTATTATTA[ <<< AluY C <<< ]ATTATTATTAACCTATTTTACAGATGAGGG
                                           ‾‾‾‾‾‾‾‾‾‾                  ‾‾‾‾‾‾‾‾‾‾

   PASS

   3' end: nothing to extend
   5' end: nothing to extend
   SVname  TEname  Strand  Div     AlnLen  MM      Gaps    5P_TSD_end      5P_offset       3P_TSD_start    3P_offset       5P_TSD  3P_TSD
   HG002_mat.svim_asm.DEL.1014     AluY    C       2.2     10      0       0       -1      0       1       0       ATTATTATTA      ATTATTATTA      PASS
   ```

## Mammalian filters `--mammal`

In order to account for the particularities of several TE families, we have introduced a `--mammal` flag that will search for specific features associated with mammalian TEs. So far we are accounting for two particular cases: 5' Inversion of L1 elements and VNTR polymorphism between orthologous SVA insertions. We will try to add more of these filters, for example to detect solo vs full-length LTR polymorphisms. If you would like to see more of these filters, please share your suggestions on the [Issue](https://github.com/cgroza/GraffiTE/issues) page!

## L1 5' inversion

SV detected by GraffiTE and corresponding to non-canonical TPRT (Twin Priming Reverse Transcription), such as Twin Priming (see [here](https://genome.cshlp.org/content/11/12/2059.long) and [here](https://mobilednajournal.biomedcentral.com/articles/10.1186/1759-8753-1-7)) may be skipped by the TSD script because it artificially creates 2 hits instead of one for a single TE insert. 

![](https://i.imgur.com/YfukCpL.png)

Whether or not the L1 is inserted on the + or - strand, at Twin-Primed L1 will have the same pattern with RepeatMasker:
- hit 1 = C
- hit2 = + 

![](https://i.imgur.com/NfyCXZd.png)
> This is because an inversion on the - strand feature will look like + on the consensus (`(-)*(-) = (+)` or a "reverted reverse")

However, we can differentiate the two based on the coordinates of the hit on the TE consensus (cartoon not to scale to compare two L1 insertions with the same consensus):

![](https://i.imgur.com/XtS5FGQ.png)

For each pair (C,+) of hits, we look at the target hit coordinates:
- if hit 1 ( C ) coordinates are < hit 2 (+), the TE inserted on the + strand (top, blue example)
- if hit 1 ( C ) coordinates are > hit 2 (+), the TE inserted on the - strand (bottom, orange example)

L1 inversions will be reported with the flag `mam_filter_1=5P_INV` in the INFO field of the VCFs.

## VNTR polymorphisms in SVA elements

<img src="https://i.imgur.com/l0DPyRL.png" alt="drawing" width="500"/>

If `GraffiTE` detects:
- SV annotated as SVA **and**,
- RepeatMasker hit corresponding only to the VNTR region of these elements **and**,
- If the flanking is an SVA in the same orientation

The variant will be flagged with `mam_filter_2=VNTR_ONLY:SVA_F:544:855` with `SVA_F:544:855` varying according to the element family and VNTR region:

| SVA model   | VNTR period size | Repeat # | start | end |
| ----------- | ----------- | ----------- |----------- |----------- |
| SVA_A       | 37          | 10.5        | 436 |855 |
| SVA_B       | 37          | 10.8        | 431| 867 |
| SVA_C       | 37          | 10.5        | 432| 851|
| SVA_D       | 37          | 6.4         | 432| 689|
| SVA_E       | 37          | 10.8        | 428| 864|
| SVA_F       | 37          | 10.5        | 435| 857|


## `GraffiTE` execution profiles
By default, the pipeline will inherit the `nextflow` configuration and run accordingly.
To execute locally, on SLURM, or AWS, pass one of the `-profile` provided with the `GraffiTE`:
- `standard`
- `cluster`
- `cloud`

For example,
```
nextflow run cgroza/GraffiTE -profile cluster ...
```
will run on SLURM.

## Specifying memory and CPU allocation at each step
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

## Resource usage example:

(this section will be updated based on our ongoing tests)

- Human chromosome 1: 10 cpu, 80Gb RAM. SV discovery \~30mn to 1h per genome, but can be improved by fine tuning the process-specific parameters.
- *Drosophila melanogaster* full genomes: 4 cpu, 40Gb RAM. SV discovery \~15mn per genome.

## Known Issues / Notes / FAQ

- The "stitching" method to identify unique TE insertion from fragmented hits has some degree of limitation. This can be flagrant for full-length LTR insertion, which can show `n_hits` > 1, and thus wont be recognized as a "single" element insertion, nor run through the TSD module. For now, names between LTR and I(nternal) sequences much match in the header name (e.g. TIRANT_LTR and TIRANT_I) to be automatically recognized as a single hit. We will make use of the RepeatMasker hit ID in order to improve this stitching procedure. In the meantime, we recommend to check/rename your LTR of interest in the `--TE_library` file. 

- As mentioned above, in order to improve runtime, the TSD module is only run for SVs with a single TE hit. We will improve this feature in order to be able to run the module on all SVs.

- The TSD module will currently spawn one process per TSD, which can create a lot of folders and files. Make sure to delete the `work/` folder regularly to stay below quotas!

- There are currently several bottlenecks in the pipeline: `samtools sort` can be tricky to parallelize properly (piped from `minimap2` alignments, which are often fast) and the performance will depends on the genomes size, complexity and the parameter used. `RepeatMasker` can be slow with a large number of SVs and a large library, hang-on! If you find satisfactory combinations of parameters for your model, please share them in the issues section! Thanks!
