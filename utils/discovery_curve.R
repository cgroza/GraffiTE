#!/usr/bin/env Rscript

library(vcfR)
library(ggrepel)
library(tidyverse)
library(optparse)

discovery_curve <- function(ref_vcf, type = c("polymorphism", "TE"), nperm = 100) {
  vcf_matrix <- lapply(str_split(extract.info(ref_vcf, "SUPP_VEC"), ""), as.numeric)
  na_entries <- is.na(vcf_matrix)
  vcf_matrix <- vcf_matrix[!na_entries] %>%
    simplify2array() %>%
    t()

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
  bind_rows(lapply(1:1000, permute_curve)) %>%
    group_by(genome) %>%
    summarise(
      MeanTEs = mean(TEs),
      SdTEs = sd(TEs),
      nperm = n(),
      lowCI = MeanTEs - qnorm(0.975) * (SdTEs / sqrt(nperm)),
      highCI = MeanTEs + qnorm(0.975) * (SdTEs / sqrt(nperm)),
      )
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

plot_curves <- function(vcf, rare_freq, fixed_freq, nperm, type = c("polymorphism", "TE")) {
  theme_set(theme_bw(base_size = 14))

  n_genomes <- str_length(extract.info(vcf, "SUPP_VEC")[1])
  rare_freq_break <- ceiling(n_genomes * rare_freq)
  fixed_freq_break <- ceiling(n_genomes * fixed_freq)

  df <- bind_rows(
    discovery_curve(filter_vcf_freq(vcf, 1, rare_freq_break, nhits = "any"), type, nperm) %>%
    mutate(Frequency = "Rare", Variants = "All"),
    discovery_curve(filter_vcf_freq(vcf, rare_freq_break, fixed_freq_break, nhits = "any"), type, nperm) %>%
    mutate(Frequency = "Common", Variants = "All"),
    discovery_curve(filter_vcf_freq(vcf, fixed_freq_break, n_genomes, nhits = "any"), type, nperm) %>%
    mutate(Frequency = "Fixed", Variants = "All"),
    discovery_curve(filter_vcf_freq(vcf, 1, rare_freq_break, nhits = "single"), type, nperm) %>%
    mutate(Frequency = "Rare", Variants = "nhits=1"),
    discovery_curve(filter_vcf_freq(vcf, rare_freq_break, fixed_freq_break, nhits = "single"), type, nperm) %>%
    mutate(Frequency = "Common", Variants = "nhits=1"),
    discovery_curve(filter_vcf_freq(vcf, fixed_freq_break, n_genomes, nhits = "single"), type, nperm) %>%
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

option_list <- list(
  make_option(c("-p", "--nperm"),
              type = "numeric", help = "Number of genome order permutations", metavar = "number",
              default = 100
  ),
  make_option(c("-v", "--vcf"),
    type = "character", help = "GraffiTE VCF file with TE genotypes.", metavar = "file"
  ),
  make_option(c("-o", "--out"),
    type = "character", help = "Prefix to write the output plots (PDF).", metavar = "file"
  ),
  make_option(c("-r", "--rare"),
              type = "numeric", help = "Frequency (0-1) breakpoint for rare alleles.", metavar = "value",
              default = 0.10
  ),
  make_option(c("-f", "--fixed"),
              type = "numeric", help = "Frequency (0-1) breakpoint for fixed alleles.", metavar = "value",
              default = 0.95
  )

)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (any(is.null(opt$vcf), is.null(opt$out))) {
  print_help(opt_parser)
  stop()
}

vcf <- read.vcfR(opt$vcf)
plot_curves(vcf, opt$rare, opt$fixed, opt$nperm, "TE")
ggsave(paste0(opt$out, "_TEs.pdf"))
plot_curves(vcf, opt$rare, opt$fixed, opt$nperm, "polymorphism")
ggsave(paste0(opt$out, "_polymorphisms.pdf"))
