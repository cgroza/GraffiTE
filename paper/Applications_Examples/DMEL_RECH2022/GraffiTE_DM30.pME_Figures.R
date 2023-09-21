### Analysis of Rech's (DM30) pangenome with GraffiTE
## We will use here the full GT-sv-sn-GA and filter for 1 hits 
library(ggplot2)
library(ggrepel)
library(dplyr)
library(grDevices)

# load main table converted from the GraffiTE VCF
GTDM<-read.table("DM30_sv-sn_graphaligner_1hits_08212023_withCounts_withHeader.tsv", 
                 h = T)
length(subset(GTDM, (GTDM$n_hits == 1))[,1])
# # subset to trusted pME (n_hits = 1 and not fixed relative to dm6 after GraphAligner)
GTDM.f<-subset(GTDM, (GTDM$n_hits == 1 & GTDM$GraphAligner_genome_counts != 0))

# join the TE classification
MCTE<-read.table("Rech_TE_list.txt", h = T)
colnames(GTDM.f)[31]<-"TEname"
GTDM.f<-merge(GTDM.f, MCTE, by = "TEname")
# list the chr we will use
DMchrs<-c("NC_004354.4", "NT_033779.5","NT_033778.4", "NT_037436.4", "NT_033777.3")
GTDM.f<-GTDM.f[GTDM.f$CHROM %in% DMchrs,]
# annotate euchromatic and heterochromatic pME
Euch.bed<-read.table("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/Other\ computers/Zephyrantes/GraffiTE/Comparison_GraffiTE_Rech_DM30_plus_ISO-1/GraffiTE_VCFS/Euchromatic_regions_RECH2022.bed")
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
rpmsk<-read.table("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/Other\ computers/Zephyrantes/GraffiTE/Comparison_GraffiTE_Rech_DM30_plus_ISO-1/dm6_rpmsk_nochr.bed")
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
label(cov_all$type)
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

######### Saturation curves ##################

library(vcfR)
library(ggrepel)
library(tidyverse)
library(optparse)
library(ggplot2)
library(reshape2)
library(ade4)
library(stringr)

# we need to load these functions first
theme_set(theme_bw(base_size = 14))

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


vcfDM<-read.vcfR("DM30_sv-sn_graphaligner_1hits_08212023_GA_SUPP.recode.vcf")
t_2_29<-plot_curves(vcfDM, rare_freq = 2/30, fixed_freq = 29/30, nperm = 100, type = "polymorphism", mode = "clem")


CP1<-cowplot::plot_grid(DM_reps, DMpME, ncol = 1, align = "hv", rel_heights = c(0.5,1), axis = "r")
CP2<-cowplot::plot_grid(t_2_29, bpDM,
                   rel_widths = c(0.3, 0.7),
                   ncol = 2
                   )
cowplot::plot_grid(CP1, CP2, ncol = 1)
