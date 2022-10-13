#!/usr/bin/env Rscript

library(biomartr)
library(optparse)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(vcfR)

read_rm_custom <- function(file) {
    rm_file <- readr::read_lines(file = file, skip = 3)
    rm_file <- lapply(rm_file, function(x) {
        str.res <- unlist(stringr::str_split(x, "\\s+"))#[-1]
        str.res <- str.res[1:16]
        return(str.res)
    })
    rm_file <- tibble::as_tibble(do.call(rbind, rm_file))
    colnames(rm_file) <- c("sw_score", "perc_div", "perc_del",
        "perc_insert", "qry_id", "qry_start", "qry_end", "qry_left", "strand",
        "repeat_id", "matching_class", "in_repeat_start", "in_repeat_end", "in_repeat_left", "ID", "fragmts")
    qry_end <- qry_start <- NULL
    nrow_before_filtering <- nrow(rm_file)
    suppressWarnings(rm_file <- dplyr::mutate(rm_file,
      qry_start = as.integer(qry_start),
      qry_end = as.integer(qry_end), fragmts = as.integer(fragmts)
    ))
    rm_file <- dplyr::filter(rm_file, !is.na(qry_start), !is.na(qry_end))
    rm_file <- dplyr::mutate(rm_file, qry_width = as.integer(qry_end -
        qry_start + 1L))
    nrow_after_filtering <- nrow(rm_file)
    if ((nrow_before_filtering - nrow_after_filtering) > 0)
        message((nrow_before_filtering - nrow_after_filtering) +
            1, " out of ", nrow_before_filtering, " rows ~ ",
            round(((nrow_before_filtering - nrow_after_filtering) +
                1) / nrow_before_filtering, 3), "% were removed from the imported RepeatMasker file, ",
            "because they contained 'NA' values in either 'qry_start' or 'qry_end'.")
    return(rm_file)
}


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

rep_mask <- read_rm_custom(opt$dotout) %>%
  # compile RM hit stats for each SV
  mutate(match_len = qry_end - qry_start) %>%
  arrange(qry_start) %>%
  group_by(qry_id) %>%
  summarise(
    start = first(qry_start), stop = max(qry_end),
    match_lengths = paste0(match_len, collapse = ","),
    fragmts = paste0(fragmts, collapse = ","),
    repeat_ids = paste0(repeat_id, collapse = ","),
    matching_classes = paste0(matching_class, collapse = ","),
    strands = paste0(strand, collapse = ","),
    RM_id = paste0(ID, collapse = ","),
    n_hits = n()
  )
  # merge non-overlapping TE annotation in the same query
  # mutate(total_match_length = stop - start) %>%
  # group_by(qry_id) %>%
  #   summarise(
  #     n_hits = sum(n_hits),
  #     total_match_length = sum(total_match_length),
  #     match_lengths = paste0(match_lengths, collapse = ","),
  #     fragmts = paste0(fragmts, collapse = ","),
  #     repeat_ids = paste0(repeat_ids, collapse = ","),
  #     matching_classes = paste0(matching_classes, collapse = ",")
  #   )

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
                  #total_match_length = 0,
                  fragmts = "0",
                  match_lengths = "0",
                  strands = "None",
                  RM_id = "None")) %>%
  #mutate(total_match_span = total_match_length / qry_length) %>%
    select(-c(qry_id, qry_length)) %>%
    arrange(CHROM, POS) %>%
    select(CHROM, POS, n_hits, fragmts, match_lengths, repeat_ids, matching_classes, strands, RM_id) #, total_match_length, total_match_span)

write_tsv(annot, file = opt$annotation, col_names = F)
print(colnames(annot))
