#!/usr/bin/env Rscript

library(biomartr)
library(optparse)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(vcfR)

option_list <- list(
  make_option(c("-r", "--dotout"),
    type = "character", help = "RepeatMasker annotation .out file.", metavar = "file"
  ),
  make_option(c("-v", "--vcf"),
    type = "character", help = "VCF file to be annotated.", metavar = "file"
    ),
  make_option(c("-a", "--annotation"),
    type = "character", help = "File to write the final annotation file for bcftools.", metavar = "file"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


if (any(is.null(opt$dotout), is.null(opt$vcf), is.null(opt$annotation))) {
  print_help(opt_parser)
  stop()
}

rep_mask <- read_rm(opt$dotout) %>%
  mutate(match_len = qry_end - qry_start) %>%
  group_by(qry_id) %>%
  filter(match_len == max(match_len)) %>%
  filter(sw_score == max(sw_score)) %>%
  select(matching_class, repeat_id, match_len)

vcf <- read.vcfR(opt$vcf)
vcf_df <- tibble(
  CHROM = getCHROM(vcf),
  POS = getPOS(vcf),
  qry_length = abs(str_length(getALT(vcf)) - str_length(getREF(vcf))),
  qry_id = getID(vcf)
)

annot <- left_join(vcf_df, rep_mask, by = "qry_id") %>%
  replace_na(list(matching_class = "None", repeat_id = "None", match_len = 0)) %>%
  mutate(match_span = match_len / qry_length) %>%
  select(-c(qry_id, qry_length, match_len)) %>%
  arrange(CHROM, POS)

write_tsv(annot, file = opt$annotation, col_names = F)
print(colnames(annot))
