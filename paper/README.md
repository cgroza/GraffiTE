# GraffiTE paper // supporting data and code

This directory contains analysis scripts associated with the GraffiTE paper. 
Associated datasets can be downloaded from [Zenondo]() (the archive uses the same directory topology).

>note (09-11-24): I am currently awaiting the BioRxiv DOI to push the data to Zenondo, for specifics, please [email me](mailto:goubert.clement@gmail.com) directly -- CG.

## Simulations

- simulator: https://github.com/clemgoub/pMEsim
- R analyses: `GraffiTE_simulations.R`
- data: `sim_output/` directory (Zenondo archive): VCF from the simulations (see paper's Methods).

## GIAB/HG002 Benchmark

- R analyses: `GraffiTE.GIAB.sveval.R`
- data: `sveval_VCFS` directory (Zenondo archive): the directory contains the filtered (see paper's method) VCF files for each method and different coverage levels. See also R code file.

## Application Examples

### Human (HPRC)

- R analyses: `GraffiTE_HPRC_pangenome.R`
- data: `HPRC_pangenome.vcf` (Zenondo archive): Raw GraffiTE VCF after step 1 and 2 (no graph genotyping)

### Drosophila melanogaster (Rech. 2022)

- R analyses: `GraffiTE_DM30.pME_Figures`
- data:
	- `DM30_sv-sn_graphaligner_1hits_08212023_withCounts_withHeader.tsv` (Zenondo archive): annotated main output in TSV format
	- `DM30_sv-sn_graphaligner_1hits_08212023_GA_SUPP.recode.vcf` (Zenondo archive): Filtered GraffiTE VCF (see paper's Methods)

### Cannabis sativa

- R analyses: `Cs_GraffiTE.R`
- data:
	- `N3_200_40000bp_clusters_pangenome.vcf.recode.vcf` (Zenondo archive): Filtered GraffiTE VCF, retaining pMEs from clusters with N = 3 pMEs, and variants length:  200 <= L <= 40000 (bp) -- see also paper's Methods
	- `Cs_GraffiTE_280423_pangenome.variants_clustered_mode0_with_ANNOTATION` (Zenondo archive): Clustering of pME locus sequences detected by GraffiTE, tsv format.
