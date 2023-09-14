# GraffiTE C sativa analyses

library(vcfR)
library(ggrepel)
library(tidyverse)
library(optparse)
library(ggplot2)

# Reusing Cristian's code to make PCA and sat curves // need to be loaded before anything
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
      scale_y_log10()+
      scale_x_continuous(breaks = c(1:9))+
      theme_classic()
  }
}

#########################################################
# Analyses
#########################################################

#### CLUSTERS ####

ClusAll<-read.table("~/Documents/Csativa/GraffiTE/Cs_GraffiTE_280423_pangenome.variants_clustered_mode0_with_ANNOTATION", h = F)
#head(ClusAll)
names(ClusAll)<-c("cluster", "member", "SUPP", "SUPP_VEC", "SVLEN", "n_hits", "n_frags", "hit_length", "hit_names", "hit_class", "hit_strands", "frags", "cum_hits_len", "SV_TE_cov", "TSD")
Ncount<-tapply(ClusAll$member, ClusAll$cluster, length)
NcumL<-tapply(ClusAll$SVLEN, ClusAll$cluster, function (x) sum(abs(x)))

# summarize

ClusN3<-ClusAll %>%
  group_by(cluster) %>%
  summarise(
    Count = length(member),
    mSVLEN = mean(abs(SVLEN)),
    sdSVLEN = sd(abs(SVLEN)),
    cumSVLEN = sum(abs(SVLEN))
    ) %>%
  filter(Count > 2) %>%
  arrange(desc(Count), desc(cumSVLEN))
write.table(ClusN3, file = "~/Documents/Csativa/GraffiTE/ClusN3.tsv", quote = F, row.names = F)

# All N >= 3
barplot(height = ClusN3$Count, width = ClusN3$cumSVLEN)
# M >= 100
barplot(height = ClusN3[ClusN3$Count > 99,]$Count, width = ClusN3[ClusN3$Count > 99,]$cumSVLEN)
# top 100
barplot(height = ClusN3[1:100,]$Count, width = ClusN3[1:100,]$cumSVLEN/1e6, space = 0, xlab = "cummulative pME sequences length (Mbp)", ylab = "copies",
        ylim = c(0,800))
axis(side = 1, at = seq(0,250, by = 1))


#### VCF ####
vcf_N3<-read.vcfR("/Users/cgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/Other\ computers/Zephyrantes/N3_200_40000bp_clusters_pangenome.vcf.recode.vcf")

# Saturation
plot_curves(vcf_N3, 1/9, 9/9, 100, "polymorphism", "clem")