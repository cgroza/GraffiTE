# GraffiTE paper // supporting data and code

This directory contains analysis scripts associated with the GraffiTE paper. 
Associated datasets can be downloaded from [zenodo](https://zenodo.org/record/8400868) (the archive uses the same directory topology).

## Simulations

- simulator: https://github.com/clemgoub/pMEsim
- R analyses: `GraffiTE_simulations.R`
- data: `sim_output/` directory (zenodo archive): VCF from the simulations (see paper's Methods).

## GIAB/HG002 Benchmark

- R analyses: `GraffiTE.GIAB.sveval.R`
- data: `sveval_VCFS` directory ([zenodo archive](https://zenodo.org/record/8400868)): the directory contains the filtered (see paper's method) VCF files for each method and different coverage levels. See also R code file.
	- `GIAB.Tier1.Alu.L1.SVA.250bp.vcf` ([zenodo archive](https://zenodo.org/record/8400868)): "Truth" VCF conaining high-quality calls for Alu, L1 and SVA elements in the individual HG002. Derived from [Genome In A Bottle](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/).
 
## Application Examples

### Human (HPRC)

- R analyses: `GraffiTE_HPRC_pangenome.R`
- data: `HPRC_pangenome.vcf` ([zenodo archive](https://zenodo.org/record/8400868)): Raw GraffiTE VCF after step 1 and 2 (no graph genotyping)

### Drosophila melanogaster (Rech. 2022)

- R analyses: `GraffiTE_DM30.pME_Figures`
- data:
	- `DM30_sv-sn_graphaligner_1hits_08212023_withCounts_withHeader.tsv` ([zenodo archive](https://zenodo.org/record/8400868)): annotated main output in TSV format
	- `DM30_sv-sn_graphaligner_1hits_08212023_GA_SUPP.recode.vcf` ([zenodo archive](https://zenodo.org/record/8400868)): Filtered GraffiTE VCF (see paper's Methods)

### Cannabis sativa

- R analyses: `Cs_GraffiTE.R`
- data:
	- `N3_200_40000bp_clusters_pangenome.vcf.recode.vcf` ([zenodo archive](https://zenodo.org/record/8400868)): Filtered GraffiTE VCF, retaining pMEs from clusters with N = 3 pMEs, and variants length:  200 <= L <= 40000 (bp) -- see also paper's Methods
	- `Cs_GraffiTE_280423_pangenome.variants_clustered_mode0_with_ANNOTATION` ([zenodo archive](https://zenodo.org/record/8400868)): Clustering of pME locus sequences detected by GraffiTE, tsv format.
