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
  # merge overlapping TE annotations in the same query
  mutate(match_len = qry_end - qry_start) %>%
  arrange(qry_start) %>%
  group_by(qry_id, g = cumsum(cummax(lag(qry_end, default = first(qry_end))) < qry_start)) %>%
  summarise(
    start = first(qry_start), stop = max(qry_end),
    match_lengths = paste0(match_len, collapse = ","),
    repeat_ids = paste0(repeat_id, collapse = ","),
    matching_classes = paste0(matching_class, collapse = ","),
    n_hits = n()
  ) %>%
  # merge non-overlapping TE annotation in the same query
  mutate(total_match_length = stop - start) %>%
  group_by(qry_id) %>%
    summarise(
      n_hits = sum(n_hits),
      total_match_length = sum(total_match_length),
      match_lengths = paste0(match_lengths, collapse = ","),
      repeat_ids = paste0(repeat_ids, collapse = ","),
      matching_classes = paste0(matching_classes, collapse = ",")
    )

vcf <- read.vcfR(opt$vcf)
vcf_df <- tibble(
  CHROM = getCHROM(vcf),
  POS = getPOS(vcf),
  REF = getREF(vcf),
  ALT = getALT(vcf),
  qry_length = abs(str_length(getALT(vcf)) - str_length(getREF(vcf))),
  qry_id = getID(vcf)
)

annot <- left_join(vcf_df, rep_mask, by = "qry_id") %>%
  replace_na(list(matching_classes = "None",
                  repeat_ids = "None",
                  n_hits = 0,
                  total_match_length = 0,
                  match_lengths = "0")) %>%
  mutate(total_match_span = total_match_length / qry_length) %>%
    select(-c(qry_id, qry_length)) %>%
    arrange(CHROM, POS) %>%
    select(CHROM, POS, n_hits, match_lengths, repeat_ids, matching_classes, total_match_length, total_match_span)

write_tsv(annot, file = opt$annotation, col_names = F)
print(colnames(annot))
