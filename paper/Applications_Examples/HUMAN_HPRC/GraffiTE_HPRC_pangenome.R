# Human HPRC / GraffiTE analyses

library(vcfR)
library(ggrepel)
library(tidyverse)
library(optparse)
library(ggplot2)
library(reshape2)
library(ade4)
library(ggnewscale)
library(ggbeeswarm)
library(stringr)

# Reusing Cristian's code to make PCA and sat curves // this needs to be loaded first
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
      labs(title = paste0("Number of ", type, "s (permutation)"),
           subtitle = "Stratified in mutually exclusive frequency bins",
           x = "Nth additional genome", y = paste0("Cumulative ", type, "s")) +
      scale_y_log10(breaks = c(300, 1000, 3000, 8000))+
      scale_x_continuous(breaks = seq(0,95, by = 10))+
      theme_classic()+
      theme(legend.position = c(60, 1000))
  }
}

#############################
### Analyses
#############################

HPRCvcf<-read.vcfR("HPRC_pangenome.vcf")
Lfilter<-abs(as.numeric(extract.info(HPRCvcf, "SVLEN"))) >= 250
HitsFilter<-as.numeric(extract.info(HPRCvcf, "n_hits")) == 1 | extract.info(HPRCvcf, "mam_filter_1") != "None"
SVAVNTF_filter<-extract.info(HPRCvcf, "mam_filter_2") == "None"
TEclassFilter<-grepl(paste(c("Alu", "L1", "SVA"), collapse = "|"), extract.info(HPRCvcf, "repeat_ids"))
HPRCvcf_F<-HPRCvcf[Lfilter & HitsFilter & TEclassFilter & SVAVNTF_filter]

## Saturation
Sat<-plot_curves(HPRCvcf_F, rare_freq = 0.05, fixed_freq = 0.95, nperm = 100, type = "polymorphism", mode = "clem")
## PCoA
HPRCmeta<-read.table("~/Documents/GraffiTE/HPRC_metadata.txt", h = T, sep = "\t")
HPRCmeta2<-HPRCmeta[rep(seq_len(nrow(HPRCmeta)), each = 2), ]
HPRCmeta2$ID<-paste(HPRCmeta2$Sample, rep(c("mat", "pat"), 47), sep = "_")
HPRCmeta2<-HPRCmeta2[order(HPRCmeta2$Sample),]
HPRCmeta2$Subpopulation<-factor(HPRCmeta2$Subpopulation,
                                levels = c("ACB", "ESN", "GWD", "MSL", "KEN", "YRI", "ASW",
                                           "CLM", "PEL", "PUR",
                                           "ASJ",
                                           "CHI", "CHS", "KHV",
                                           "PJL"),
                                labels = c("African Caribbean, Barbados (ACB)", "Esan, Nigeria (ESN)", "Gambian, Western Division (GWD)", "Mende, Sierra Leone (MSL)", "Maasai, Kenya (KEN)", "Yoruba, Nigeria (YRI)", "African Ancestry, South West US (ASW)",
                                           "Colombians, Medelin (CLM)", "Peruvian, Lima (PEL)", "Puerto Rican, Puerto Rico (PUR)", 
                                           "Ashkenazim Jewish, Europe, GIAB (ASJ)", 
                                           "Chinese, GIAB (CHI))", "Southern Han Chinese (CHS)", "Kinh Vietnamese (KHV)", 
                                           "Punjabi (PJL)"))
HPRCcolors<-c("#F5D470","darkorange", "#EAAA51", "#DC5936", "#E08244", "gold", "firebrick",
              "#ad22f2", "#8959B9", "#BBA9C4",
              "#E69FC6",
              "#3971C6", "royalblue", "#829CB3",
              "#39bf44")
pca_matrix_H <- discovery_curve(HPRCvcf_F, "polymorphism", 1)$pca
HPinds<-HPRCmeta2$ID

# for counting purposing only, count per individual instead that per haplotype
HPdiplo<-str_split(HPinds, "_", n = 2,simplify = T)[,1]
IndVarCount<-aggregate(pca_matrix_H, data.frame(HPdiplo), sum)
row.names(IndVarCount)<-IndVarCount$HPdiplo
IndVarCount<-IndVarCount[,-1]
IndVarCount[IndVarCount == 2]<-1
## Total
mean(rowSums(IndVarCount))
sd(rowSums(IndVarCount))
## Per TE class
ALUfilter<-grepl(paste(c("Alu"), collapse = "|"), extract.info(HPRCvcf_F, "repeat_ids"))
L1filter<-grepl(paste(c("L1"), collapse = "|"), extract.info(HPRCvcf_F, "repeat_ids"))
L1REVfilter<-grepl(paste(c("5P_INV"), collapse = "|"), extract.info(HPRCvcf_F, "mam_filter_1"))
SVAfilter<-grepl(paste(c("SVA"), collapse = "|"), extract.info(HPRCvcf_F, "repeat_ids"))
VNTRfilter<-grepl(paste(c("VNTR"), collapse = "|"), extract.info(HPRCvcf, "mam_filter_2"))
table(ALUfilter)
table(L1filter)
table(L1REVfilter)
table(SVAfilter)
table(VNTRfilter)
## For each individual (two haplotypes)
mean(rowSums(IndVarCount[,ALUfilter]))
sd(rowSums(IndVarCount[,ALUfilter]))

mean(rowSums(IndVarCount[,L1filter]))
sd(rowSums(IndVarCount[,L1filter]))

mean(rowSums(IndVarCount[,L1REVfilter]))
sd(rowSums(IndVarCount[,L1REVfilter]))
mean(rowSums(IndVarCount[,L1REVfilter]))/mean(rowSums(IndVarCount[,L1filter]))

mean(rowSums(IndVarCount[,SVAfilter]))
sd(rowSums(IndVarCount[,SVAfilter]))

VNTR_Vcf<-HPRCvcf[VNTRfilter]
VNTRmat<-discovery_curve(VNTR_Vcf, "polymorphism", 1)$pca
VNTRmatCount<-aggregate(VNTRmat, data.frame(HPdiplo), sum)
row.names(VNTRmatCount)<-VNTRmatCount$HPdiplo
VNTRmatCount<-VNTRmatCount[,-1]
VNTRmatCount[VNTRmatCount == 2]<-1
mean(rowSums(VNTRmatCount))
sd(rowSums(VNTRmatCount))

# back to the binary matrix for PCoA using 1/0 calls in each haplotype
Jmat_H<-dist.binary(pca_matrix_H[,colSums(pca_matrix_H) > 0.05*94 & colSums(pca_matrix_H) < 0.95*94], method = 5)
pco_H<-dudi.pco(na.omit(Jmat_H), nf = 3, scannf = F)
insH<-inertia.dudi(pco_H)
pco_H_table<-cbind(pco_H$li, HPinds, HPRCmeta2$Subpopulation)
# PCoA
HPRCmeta2$SubCode<-factor(HPRCmeta2$Subpopulation, levels = c("African Caribbean, Barbados (ACB)", "Esan, Nigeria (ESN)", "Gambian, Western Division (GWD)", "Mende, Sierra Leone (MSL)", "Maasai, Kenya (KEN)", "Yoruba, Nigeria (YRI)", "African Ancestry, South West US (ASW)",
                                                              "Colombians, Medelin (CLM)", "Peruvian, Lima (PEL)", "Puerto Rican, Puerto Rico (PUR)", 
                                                              "Ashkenazim Jewish, Europe, GIAB (ASJ)", 
                                                              "Chinese, GIAB (CHI))", "Southern Han Chinese (CHS)", "Kinh Vietnamese (KHV)", 
                                                              "Punjabi (PJL)"),
                          labels = c("ACB", "ESN", "GWD", "KEN", "MSL", "YRI", "ASW",
                                              "CLM", "PEL", "PUR",
                                              "ASJ",
                                              "CHI", "CHS", "KHV",
                                              "PJL")
)
PCoA<-ggplot(pco_H_table, aes(x = A1, y = A2, label = HPinds, col = HPRCmeta2$Subpopulation))+
  geom_point()+
  #geom_label_repel()+
  xlab(paste("Axis 1 (", round(insH$tot.inertia$`cum(%)`[1], 1), "%)"))+
  ylab(paste("Axis 2 (", round(insH$tot.inertia$`cum(%)`[2]-insH$tot.inertia$`cum(%)`[1], 1), "%)"))+
  scale_color_manual(values = HPRCcolors)+
  theme_classic()

## Count per TE and Pop

TEmat<-data.frame(cbind(t(pca_matrix_H), unlist((lapply(str_split(extract.info(HPRCvcf_F, "matching_classes"), ","), `[[`, 1))), getID(HPRCvcf_F)))
names(TEmat)<-c(HPinds, c("pME_Class", "pME_ID"))
TElong<-melt(TEmat, id.vars = c("pME_Class", "pME_ID")) %>%
  group_by(variable, pME_Class) %>%
    summarize(
      Count = sum(as.numeric(value))
      )
names(TElong)<-c("ID", "pME_Class", "Count")
TElonger<-merge(TElong, HPRCmeta2, by = "ID")
# ggplot(TElonger, aes(x = Region, y = Count, fill = pME_Class))+
#   geom_boxplot(outlier.colour = "white")+
#   theme_classic()
TElonger$Region<-factor(TElonger$Region, levels = c("AFR", "AMR", "ASJ", "EAS", "SAS"), labels = c("Africa", "Americas", "Europe", "Eastern Asia", "Southern Asia"))
TElonger$pME_Class<-factor(TElonger$pME_Class, levels = c("SINE/Alu", "LINE/L1", "Retroposon/SVA"))
RegionColor<-c("darkorange", "#ad22f2", "#E69FC6", "royalblue", "#39bf44")
Counts<-ggplot(TElonger, aes(x = SubCode, y = Count))+
  geom_point(aes(colour = Subpopulation), position = position_jitterdodge(jitter.width = 1), size = .6, show.legend = F)+
  #geom_beeswarm(aes(colour = Subpopulation), dodge.width = 1, cex = .5, size = .5)+#, position = position_jitterdodge(jitter.width = .15), size = .5)+
  #geom_boxplot(aes(colour = Subpopulation), outlier.colour = "white", lwd = 0.2, fatten = 3, fill = "white", alpha = .5)+
  scale_color_manual(values = HPRCcolors)+
  #scale_fill_manual(values = HPRCcolors)+
  #new_scale_color()+
  geom_boxplot(aes(colour = Subpopulation), fill="white", outlier.alpha = 0, outlier.size = .5, lwd = 0.3, fatten = 3, alpha = .8, show.legend = F)+
  #scale_color_manual(values = RegionColor)+
  scale_fill_manual(values = HPRCcolors)+
  facet_wrap(~pME_Class, scales = "free")+
  #scale_color_manual(values = RegionColor)+
  #scale_y_log10()+
  theme_classic()+
  xlab("Subpopulation")+
  ylab("pME count")+
  theme(axis.text.x = element_text(angle = 90))

Aplot<-cowplot::plot_grid(Sat, PCoA, rel_widths = c(0.35,0.65))
Bplot<-cowplot::plot_grid( NULL, Counts, rel_widths = c(0,1))
cowplot::plot_grid(Aplot, Bplot, nrow = 2)