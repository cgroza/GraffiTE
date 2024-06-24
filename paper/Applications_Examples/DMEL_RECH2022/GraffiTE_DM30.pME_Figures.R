### Analysis of Rech's (DM30) pangenome with GraffiTE
## We will use here the full GT-sv-sn-GA and filter for 1 hits 

################################################################################
# Data prep (bash commands used to prepare the data for R)
# The starting file is the output if Graffite mode GT-svsn-GA DM30_sv-sn_graphaligner.vcf
# 1. list variants with a single TE hit using the pre-genotyping VCF:
# grep -v '#' DM30_sv-sn_pangenome.vcf | grep 'n_hits=1;' | cut -f 3 > DM30_sv-sn_pangenome_1hits_IDs
# 2. filter for this variants and exclude fixed loci (same genotype as reference genome for all samples)
# vcftools --vcf  DM30_sv-sn_graphaligner.vcf --snps DM30_sv-sn_pangenome_1hits_IDs --non-ref-ac-any 1 --recode --recode-INFO-all --out DM30_sv-sn_GA_1hits_03142024
# 3. convert filtered VCF into tsv for R analyses
# cat GA_tsv_head <(paste <(grep -v '#' DM30_sv-sn_GA_1hits_03142024.recode.vcf | cut -f 1-7) <(grep -v '#' DM30_sv-sn_GA_1hits_03142024.recode.vcf | cut -f 8 | sed 's/;/\t/g;s/[A-Za-z0-9_-]*=//g;s/\tlowdepth//g' | awk '$2~/>/ {if (NF == 29) {print $3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29"\tNA\t"$1"\t"$2"\t"} else {print $3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30"\t"$1"\t"$2}; next}  {if (NF == 29) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\tNA\t"$28"\t"$29} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30}}') <(grep -v '#' DM30_sv-sn_GA_1hits_03142024.recode.vcf | cut -f 10- | sed 's/\.\/\.:/N\/A:/g' | awk '{ for (i=1; i<=NF; i++) printf "%s ", substr($i, 1, 3); print ""; }') <(grep -v '#' DM30_sv-sn_GA_1hits_03142024.recode.vcf | cut -f 10- | sed 's/\.\/\.:/N\/A:/g' | awk '{ for (i=1; i<=NF; i++) printf "%s ", substr($i, 1, 3); print ""; }' | awk '{ count = 0; for (i=1; i<=NF; i++) { if ($i == "0/0" || $i == "N/A") count++}; {print 30-count}}') | awk -v OFS="\t" '$1=$1') | awk '/POS/ {print $0"\thasMissing"; next} /N\/A/ {print $0"\tTRUE"; next} {print $0"\tFALSE"}' > DM30_sv-sn_GA_1hits_03142024.recode.NAomitCounts.tsv
# 4. header for this file is:
# CHROM	POS	ID	REF	ALT	QUAL	FILTER	CIEND	CIPOS	CHR2	END	SVLEN	SVMETHOD	SVTYPE	SUPP_VEC	SUPP	STRANDS	sniffles2_SUPP	sniffles2_SVLEN	sniffles2_SVTYPE	sniffles2_ID	svim-asm_SUPP	svim-asm_SVLEN	svim-asm_SVTYPE	svim-asm_ID	n_hits	match_lengths	TEname	matching_classes	fragmts	RM_hit_strands	RM_hit_IDs	total_match_length	total_match_span	TSD	DP	AT	AKA-017	AKA-018	COR-014	COR-018	COR-023	JUT-008	KIE-094	RAL-059	COR-025	GIM-012	GIM-024	JUT-011	LUN-004	LUN-007	MUN-008	MUN-009	MUN-013	MUN-015	MUN-020	RAL-091	RAL-176	RAL-177	RAL-375	RAL-426	RAL-737	RAL-855	SLA-001	STO-022	TEN-015	TOM-008 GraphAligner_genome_counts
################################################################################

### R analyses start here:

# Load R libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(grDevices)

# load main table converted from the GraffiTE VCF
#GTDM<-read.table("DM30_sv-sn_graphaligner_1hits_08212023_withCounts_withHeader.tsv", 
#                 h = T)
options(scipen=999) # prevent scientific notation of genomic position and other large integers
GTDM<-read.table("DM30_sv-sn_GA_1hits_03142024.recode.NAomitCounts.tsv", 
                 h = T)
# check that the input file is already filtered for 1 hits only
# length(GTDM[,1])
# length(subset(GTDM, (GTDM$n_hits == 1))[,1])
# length(subset(GTDM, (GTDM$n_hits == 1 & GTDM$GraphAligner_genome_counts != 0))[,1])
# for code legacy and to avoid bug, rename input with .f
GTDM.f<-GTDM
# load and join the TE classification
MCTE<-read.table("Rech_TE_list.txt", h = T)
GTDM.f<-merge(GTDM.f, MCTE, by = "TEname")
# list the chr we will use
DMchrs<-c("NC_004354.4", "NT_033779.5","NT_033778.4", "NT_037436.4", "NT_033777.3")
GTDM.f<-GTDM.f[GTDM.f$CHROM %in% DMchrs,]
# remove loci with missing genotypes
GTDM.f<-GTDM.f[GTDM.f$hasMissing == FALSE,]

# annotate euchromatic and heterochromatic pME
Euch.bed<-read.table("Euchromatic_regions_RECH2022.bed")
GTDM.f.bed<-data.frame(cbind(GTDM.f$CHROM, GTDM.f$POS, GTDM.f$POS+1,GTDM.f$ID))
Euch.GTDM.f<-bedtoolsr::bt.intersect(GTDM.f.bed, Euch.bed)[,4]
GTDM.f$region<-ifelse(GTDM.f$ID %in% Euch.GTDM.f, "euch", "heter")
DMcsize<-data.frame(DMchrs, chr=c("X", "2L", "2R", "3L", "3R"), start=c(0,0,0,0,0), end=c(23542271, 23513712, 25286936, 28110227, 32079331))
names(DMcsize)<-c("CHROM", "chr", "start", "end")
#GTDM.f$CHROM<-factor(GTDM.f$CHROM, levels = c("NC_004354.4", "NT_033779.5","NT_033778.4", "NT_037436.4", "NT_033777.3", "NC_004353.4"),
#                     labels = c("X", "2L", "2R", "3L", "3R", "4"))
GTDM.f<-merge(GTDM.f, DMcsize, by = "CHROM")
# Heterochromatic regions defined in Rech 2022
Het<-data.frame(chr = c("2L", "2R", "3L", "3R", "X"),
                X1_s = rep(0,5),
                X1_e = c(530000, 5982495, 750000, 6754278, 1325967),
                X2_s = c(18870000, 24972477, 19026900, 31614278, 21338973),
                X2_e = c(23513712, 25286936, 28110227, 32079331, 23542271))
# get dm6 annotation for background densities
rpmsk<-read.table("dm6_rpmsk_nochr.bed")
names(rpmsk)<-c("chr", "start", "end", "strand", "name", "Order", "SuperFam")
rpmsk.TEs<-rpmsk[!rpmsk$Order %in% c("Simple_repeat", "Low_complexity", "Satellite"),]
rpmsk.reps<-rpmsk[rpmsk$Order %in% c("Simple_repeat", "Low_complexity", "Satellite"),]
dm6.500kbW<- read.table("dm6.500kbW.bed")
cov_reps<-bedtoolsr::bt.coverage(a = dm6.500kbW, b = rpmsk.reps)
cov_reps$mean_pos<-(cov_reps$V3+cov_reps$V2)/2
names(cov_reps)<-c("chr", "w_start", "w_end", "V4", "V5", "w_len", "density", "mean_pos")
cov_reps$type<-rep("simple_repeats", length(cov_reps$chr))
cov_tes<-bedtoolsr::bt.coverage(a = dm6.500kbW, b = rpmsk.TEs)
cov_tes$mean_pos<-(cov_tes$V3+cov_tes$V2)/2
names(cov_tes)<-c("chr", "w_start", "w_end", "V4", "V5", "w_len", "density", "mean_pos")
cov_tes$type<-rep("TEs", length(cov_tes$chr))
cov_all<-rbind(cov_reps, cov_tes)
cov_all$panel<-rep("dm6", length(cov_all$chr))

# counts
length(GTDM.f$ID)
length(GTDM.f[GTDM.f$region == "euch",]$ID)
length(GTDM.f[GTDM.f$region == "heter",]$ID)

### Barplot comparing E/H pME content
# reorder factors
GTDM.f$panel<-rep("pME", length(GTDM.f$CHROM))
GTDM.f$Order<-factor(GTDM.f$Order, levels = c("LTR", "LINE", "TIR", "Helitron", "LARD", "TRIM"))
# normalized counts in loci #, separated by INS/DEL and HET/EUC
DMpME<-ggplot(GTDM.f[GTDM.f$SuperFamily != "NewFam",])+
  geom_rect(dat = Het, aes(xmin = X1_s/1000000, xmax = X1_e/1000000, ymin = 0, ymax = 150), fill = "grey90")+
  geom_rect(dat = Het, aes(xmin = X2_s/1000000, xmax = X2_e/1000000, ymin = 0, ymax = 150), fill = "grey90")+
  geom_histogram(aes(x = POS/1000000, fill = Order), binwidth = 0.5)+
  #geom_line(data = cov_all, aes(x = mean_pos/1000000, y = density, col = type), show.legend = T)+
  scale_fill_manual(values = c("green3","royalblue", "salmon", "gold", "purple"))+
  xlab("chromosome position (Mbp)")+
  ylab("pME count")+
  facet_grid(~chr, scales = "free_x", space = "free_x")+
  #theme_minimal()+
  theme(#strip.background =element_rect(fill="grey10"),
        #strip.text = element_text(colour = 'white'),
    panel.border = element_blank(),
    strip.text = element_blank(),    
    panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        legend.justification = "left",
        )

cov_all$type<-factor(cov_all$type, levels = c("TEs", "simple_repeats"), labels = c("TEs", "other repeats"))
#label(cov_all$type)
DM_reps<-ggplot(GTDM.f[GTDM.f$SuperFamily != "NewFam",])+
  geom_rect(dat = Het, aes(xmin = X1_s/1000000, xmax = X1_e/1000000, ymin = 0, ymax = 1), fill = "grey90")+
  geom_rect(dat = Het, aes(xmin = X2_s/1000000, xmax = X2_e/1000000, ymin = 0, ymax = 1), fill = "grey90")+
  geom_line(data = cov_all, aes(x = mean_pos/1000000, y = density, col = type, linetype = type))+
  xlab("")+
  #theme_minimal()+
  scale_color_manual( values = c("black", "grey50"))+
  facet_grid(~chr, scales = "free_x", space = "free_x")+
  theme(strip.background = element_rect(fill="white"),
        #strip.text = element_text(colour = 'white'),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        legend.justification = "left",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  guides(color=guide_legend(title="dm6 annotations"))

cowplot::plot_grid(DM_reps, DMpME, ncol = 1, align = "hv", rel_heights = c(0.5,1))

levels(GTDM.f$SuperFamily)[levels(GTDM.f$SuperFamily)=="Gypsy"] <- "Ty3"
GTDM.f$SuperFamily<-factor(plyr::revalue(GTDM.f$SuperFamily,c("Gypsy"="Ty3")))
GTDM.f$SuperFamily == "Ty3"

## Count of loci per Subfamily and frequency of loci in pop
# order pop frequency
GTDM.f$genome_count<-factor(GTDM.f$GraphAligner_genome_counts, levels = c(30:1))
# create ramp palette
colfuncHIGH<-colorRampPalette(c("#9d0000", "#e94a29"))
colfuncMID<-colorRampPalette(c("#ffad00","#fed6ac"))
colfuncLOW<-colorRampPalette(c("#9fa2ff", "#6a6fff"))
# group and count per TE order
df2<-GTDM.f[GTDM.f$SuperFamily != "NewFam",] %>%
  group_by(SuperFamily) %>%
  summarise(Count = length(ID))
df2$SuperFamily<-factor(df2$SuperFamily, 
                        levels = c("Copia", "Gypsy", "Pao",
                                   "CR1", "I", "Jockey", "R1",
                                   "CMC-Transib", "hAT-hobo", "MULE-NOF", "P", "PiggyBac", "TcMar-Pogo", "TcMar-Tc1",
                                   "INE-1",
                                   "LARD"))
# plot
bpDM<-ggplot(GTDM.f[GTDM.f$SuperFamily != "NewFam",], aes(x = factor(SuperFamily), fill = factor(genome_count)))+
  geom_bar(stat = "count")+
  scale_fill_manual(values = c(colfuncHIGH(2), colfuncMID(26), colfuncLOW(2)))+
  #geom_text(data = df2, aes(x = SuperFamily, y = Count))+
  facet_grid(~Order, scales = "free", space = "free")+
  #theme_minimal()+
  ylab("pME count")+
  xlab("pME Superfamily")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "right",
        strip.background =element_rect(fill="white"),
        strip.text = element_text(colour = 'black'),
        #legend.text = element_text(size=.5),
        legend.key.size = unit(0.5, "cm"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  labs(fill = "strains \n")+
  guides(fill = guide_legend(ncol = 2, label.position = "right",
                             reverse = T)
         )

# export list of analysed variants
write.table(as.data.frame(GTDM.f$ID), "analyzed_variants", quote = F, row.names = F, col.names = F)

### Saturation curves:

################################################################################
# for this analysis, we are using a tool we created early on that require the
# variants to be encoded in the SUPP_VEC of the VCF. Since we are working on the
# genotyped variants (not the detected variants), we needed to recode the SUPP_VEC
# accordingly
#
# # first filter VCF to match the filtered list of variant
# vcftools --vcf DM30_sv-sn_GA_1hits_03142024.recode.vcf --snps analyzed_variants --recode --recode-INFO-all --out DM30_sv-sn_GA_1hits_03142024__analyzed
# # second replace the SUPP_VEC using GA genotypes
# cat <(grep '#' DM30_sv-sn_GA_1hits_03142024__analyzed.recode.vcf) <(paste <(grep -v '#' DM30_sv-sn_GA_1hits_03142024__analyzed.recode.vcf | sed 's/SUPP_VEC=/SUPP_VEC=\t/g' | cut -f 1-8) <(grep -v '#' DM30_sv-sn_GA_1hits_03142024__analyzed.recode.vcf | cut -f 10- | sed 's/\.\/\.:/N\/A:/g' | awk '{ for (i=1; i<=NF; i++) printf "%s ", substr($i, 1, 3); print ""; }' | sed -e 's/0\/0/0/g' -e 's/N\/A/0/g' -e 's/.\/./1/g' -e 's/ //g') <(grep -v '#' DM30_sv-sn_GA_1hits_03142024__analyzed.recode.vcf | sed 's/[0-1]*;SUPP=/\t;SUPP=/g' | cut -f 9-) | sed 's/SUPP_VEC=\t/SUPP_VEC=/g;s/\t;SUPP/;SUPP/g') > DM30_sv-sn_GA_1hits_03142024__analyzed_GA_SUPP.recode.vcf
################################################################################

library(vcfR)
library(ggrepel)
library(tidyverse)
library(optparse)
library(ggplot2)
library(reshape2)
library(ade4)
library(stringr)

theme_set(theme_bw(base_size = 14))

# functions to read VCF's support vector and plot the saturation curves
discovery_curve <- function(ref_vcf, type = c("polymorphism", "TE"), nperm = 100) {
  vcf_matrix <- lapply(str_split(extract.info(ref_vcf, "SUPP_VEC"), ""), as.numeric)
  na_entries <- is.na(vcf_matrix)
  vcf_matrix <- vcf_matrix[!na_entries] %>%
    simplify2array() %>%
    t()
  
  pca_matrix <- t(vcf_matrix)
  
  if (type == "TE") {
    indels <- str_length(getALT(ref_vcf)) - str_length(getREF(ref_vcf))
    indels <- indels[!na_entries]
    flip_ref_TE <- indels < 0
    vcf_matrix[flip_ref_TE, ] <- as.numeric(!vcf_matrix[flip_ref_TE, ])
  }
  
  permute_curve <- function(i) {
    vcf_m <- vcf_matrix[, sample(ncol(vcf_matrix))]
    
    vcf_vec_found <- vector(length = nrow(vcf_m))
    vcf_vec_count <- vector(mode = "numeric", length = ncol(vcf_m))
    
    for (j in 1:ncol(vcf_m)) {
      vcf_vec_found[as.logical(vcf_m[, j])] <- T
      vcf_vec_count[j] <- sum(vcf_vec_found)
    }
    
    tibble(
      TEs = vcf_vec_count,
      genome = 1:ncol(vcf_m),
      permutation = i
    )
  }
  return(list(
    curve = bind_rows(lapply(1:nperm, permute_curve)) %>%
      group_by(genome) %>%
      summarise(
        MeanTEs = mean(TEs),
        SdTEs = sd(TEs),
        nperm = n(),
        lowCI = MeanTEs - qnorm(0.975) * (SdTEs / sqrt(nperm)),
        highCI = MeanTEs + qnorm(0.975) * (SdTEs / sqrt(nperm)),
      ),
    pca = pca_matrix
  ))
}
filter_vcf_freq <- function(ref_vcf, min_freq = 0, max_freq = Inf, nhits = c("any", "single")) {
  vcf_freq <- lapply(str_split(extract.info(ref_vcf, "SUPP_VEC"), ""), as.numeric) %>%
    simplify2array() %>%
    t() %>%
    rowSums()
  
  if(nhits == "single") {
    vcf_hits <- sapply(extract.info(ref_vcf, "n_hits"), as.numeric) == 1
    vcf_filter <- str_detect(extract.info(ref_vcf, "mam_filter_1"), "5P_INV")
    
    ref_vcf[(vcf_freq >= min_freq & vcf_freq <= max_freq) & (vcf_hits | vcf_filter)]
    
  } else {
    ref_vcf[vcf_freq >= min_freq & vcf_freq <= max_freq]
  }
}
plot_curves <- function(vcf, rare_freq, fixed_freq, nperm, type = c("polymorphism", "TE"), mode = c("cristian", "clem")) {
  
  n_genomes <- str_length(extract.info(vcf, "SUPP_VEC")[1])
  rare_freq_break <- ceiling(n_genomes * rare_freq)
  fixed_freq_break <- ceiling(n_genomes * fixed_freq)
  
  if (mode == "cristian") {
    
    df <- bind_rows(
      discovery_curve(filter_vcf_freq(vcf, 1, rare_freq_break, nhits = "any"), type, nperm)$curve %>%
        mutate(Frequency = "Rare", Variants = "All"),
      discovery_curve(filter_vcf_freq(vcf, rare_freq_break, fixed_freq_break, nhits = "any"), type, nperm)$curve %>%
        mutate(Frequency = "Common", Variants = "All"),
      discovery_curve(filter_vcf_freq(vcf, fixed_freq_break, n_genomes, nhits = "any"), type, nperm)$curve %>%
        mutate(Frequency = "Fixed", Variants = "All"),
      discovery_curve(filter_vcf_freq(vcf, 1, rare_freq_break, nhits = "single"), type, nperm)$curve %>%
        mutate(Frequency = "Rare", Variants = "nhits=1"),
      discovery_curve(filter_vcf_freq(vcf, rare_freq_break, fixed_freq_break, nhits = "single"), type, nperm)$curve %>%
        mutate(Frequency = "Common", Variants = "nhits=1"),
      discovery_curve(filter_vcf_freq(vcf, fixed_freq_break, n_genomes, nhits = "single"), type, nperm)$curve %>%
        mutate(Frequency = "Fixed", Variants = "nhits=1")
    ) %>%
      mutate(
        Frequency = factor(Frequency, levels = c("Rare", "Common", "Fixed")),
        Variants = factor(Variants, levels = c("All", "nhits=1")),
        group = paste(Frequency, Variants)
      )
    
    ggplot(df) +
      geom_point(aes(x = genome, y = MeanTEs, color = Frequency, shape = Variants), alpha = 0.7) +
      geom_line(aes(x = genome, y = MeanTEs, color = Frequency, linetype = Variants, group = group), alpha = 0.7) +
      geom_ribbon(aes(x = genome, ymin = lowCI, ymax = highCI, fill = Frequency, group = group), alpha = 0.1) +
      geom_text_repel(data = df %>% filter(genome %in% c(1, ceiling(n_genomes/2), n_genomes)),
                      mapping = aes(x = genome, y = MeanTEs, color = Frequency, label = ceiling(MeanTEs))) +
      labs(title = paste0("Number of ", type, "s (permutation)"),
           subtitle = "Stratified in mutually exclusive frequency bins",
           x = "Nth additional genome", y = paste0("Cumulative ", type, "s")) +
      scale_y_log10()
  } 
  else if (mode == "clem") {
    df <- bind_rows(
      discovery_curve(filter_vcf_freq(vcf, 1, rare_freq_break, nhits = "any"), type, nperm)$curve %>%
        mutate(Frequency = "Rare", Variants = "All"),
      discovery_curve(filter_vcf_freq(vcf, rare_freq_break+1, fixed_freq_break-1, nhits = "any"), type, nperm)$curve %>%
        mutate(Frequency = "Common", Variants = "All"),
      discovery_curve(filter_vcf_freq(vcf, fixed_freq_break, n_genomes, nhits = "any"), type, nperm)$curve %>%
        mutate(Frequency = "Fixed", Variants = "All"),
    ) %>%
      mutate(
        Frequency = factor(Frequency, levels = c("Rare", "Common", "Fixed")),
        # Variants = factor(Variants, levels = c("All", "nhits=1")),
        group = paste(Frequency, Variants)
      )
    
    ggplot(df) +
      geom_point(aes(x = genome, y = MeanTEs, color = Frequency), alpha = 0.7) +
      geom_line(aes(x = genome, y = MeanTEs, color = Frequency, group = group), alpha = 0.7) +
      geom_ribbon(aes(x = genome, ymin = lowCI, ymax = highCI, fill = Frequency, group = group), alpha = 0.1) +
      geom_text_repel(data = df %>% filter(genome %in% c(1, ceiling(n_genomes/2), n_genomes)),
                      mapping = aes(x = genome, y = MeanTEs, color = Frequency, label = ceiling(MeanTEs))) +
      labs(title = paste0("Number of pME loci (permutation)"),
           subtitle = "Stratified in mutually exclusive frequency bins",
           x = "Nth additional genome", y = paste0("Cumulative pME")) +
      scale_y_log10(breaks = c(300, 1000, 3000, 8000))+
      scale_x_continuous(breaks = seq(0,95, by = 10))+
      scale_color_manual(values = c("#6a6fff", "#ffad00", "#9d0000"))+
      theme_classic()+
      theme(legend.position = c(60, 1000))
  }
}
# load the VCF
vcfDM<-read.vcfR("DM30_sv-sn_GA_1hits_03142024__analyzed_GA_SUPP.recode.vcf")
t_2_29<-plot_curves(vcfDM, rare_freq = 2/30, fixed_freq = 29/30, nperm = 100, type = "polymorphism", mode = "clem")


CP1<-cowplot::plot_grid(DM_reps, DMpME, ncol = 1, align = "hv", rel_heights = c(0.5,1), axis = "r")
CP2<-cowplot::plot_grid(t_2_29, bpDM,
                   rel_widths = c(0.3, 0.7),
                   ncol = 2
                   )

# final plot
cowplot::plot_grid(CP1, CP2, ncol = 1)

### Comparison with short reads:

################################################################################
# data prep to match long-reads analyses
#
# 1. First, I need to filter the pangenie VCF in the same fashion: 1 hits, no fixed
# ```sh
# vcftools --vcf  DM30_sv-sn_pangenie.vcf --snps DM30_sv-sn_pangenome_1hits_IDs --non-ref-ac-any 1 --recode --recode-INFO-all --out DM30_sv-sn_PA_1hits_03192024
# ```
# After filtering, kept 10950 out of a possible 32909 Site
# 
# 2. Identify missing data -- we don't need to parse the info field this time (we just need to positions and the genotypes)
# ```sh
# paste <(grep -v '#' DM30_sv-sn_PA_1hits_03192024.recode.vcf | cut -f 1-7)  <(grep -v '#' DM30_sv-sn_PA_1hits_03192024.recode.vcf | cut -f 10- | sed 's/\.\/\.:/N\/A:/g' | awk '{ for (i=1; i<=NF; i++) printf "%s ", substr($i, 1, 3); print ""; }') <(grep -v '#' DM30_sv-sn_PA_1hits_03192024.recode.vcf | cut -f 10- | sed 's/\.\/\.:/N\/A:/g' | awk '{ for (i=1; i<=NF; i++) printf "%s ", substr($i, 1, 3); print ""; }' | awk '{ count = 0; for (i=1; i<=NF; i++) { if ($i == "0/0" || $i == "N/A") count++}; {print 30-count}}') | awk -v OFS="\t" '$1=$1' | awk '/POS/ {print $0"\thasMissing"; next} /N\/A/ {print $0"\tTRUE"; next} {print $0"\tFALSE"}' | awk '{print $1"\t"$2"\t"$3"\t"$(NF-1)"\t"$NF}' > DM30_sv-sn_PG_1hits_03192024.simple.NAomitCounts.tsv
# ```
################################################################################

GTDMshort<-read.table("DM30_sv-sn_PG_1hits_03192024.simple.NAomitCounts.tsv")
names(GTDMshort)<-c("CHROM", "POS", "ID", "PGcount", "hasMissing")
# filter for main chromosomes
GTDMshort<-GTDMshort[GTDMshort$CHROM %in% DMchrs,]
# remove loci with missing genotypes
GTDMshort<-GTDMshort[GTDMshort$hasMissing == FALSE,]

library(ggVennDiagram)
vlist<-list(
  GTlong=GTDM.f$ID,
  GTshort=GTDMshort$ID)
ggVennDiagram(vlist) + scale_fill_gradient(low="grey99",high = "red")

length(GTDMshort$ID)
