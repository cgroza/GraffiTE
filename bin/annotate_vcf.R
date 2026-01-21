#!/usr/bin/env Rscript

library(optparse)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(vcfR)

read_rm_custom <- function(file) {
  rm_file <- readr::read_lines(file = file, skip = 3)
  rm_file <- lapply(rm_file, function(x) {
    str.res <- unlist(stringr::str_split(stringr::str_trim(x), "\\s+"))#[-1]
    str.res <- str.res[1:15]
    return(str.res)
  })
  rm_file <- tibble::as_tibble(do.call(rbind, rm_file))
  colnames(rm_file) <- c("sw_score", "perc_div", "perc_del",
                         "perc_insert", "qry_id", "qry_start", "qry_end", "qry_left", "strand",
                         "repeat_id", "matching_class", "in_repeat_start", "in_repeat_end", "in_repeat_left", "ID")
  # fix the cross_match confusing columns for the target
  rm_file$in_repeat_start <- as.integer(stringr::str_remove_all(rm_file$in_repeat_start, "\\(|\\)"))
  rm_file$in_repeat_end <- as.integer(stringr::str_remove_all(rm_file$in_repeat_end, "\\(|\\)"))
  rm_file$in_repeat_left <- as.integer(stringr::str_remove_all(rm_file$in_repeat_left, "\\(|\\)"))
  rm_file$target_start <- ifelse(rm_file$strand == "+", rm_file$in_repeat_start, rm_file$in_repeat_left)
  rm_file$target_end <- rm_file$in_repeat_end

  # commented section below since V1.1
  # qry_end <- qry_start <- NULL
  # nrow_before_filtering <- nrow(rm_file)
  # # suppressWarnings(rm_file <- dplyr::mutate(rm_file,
  # #                                           qry_start = as.integer(qry_start)-1, # add -1 to be 0-based and calculate length
  # #                                           qry_end = as.integer(qry_end), fragmts = as.integer(fragmts)
  # # ))
  # rm_file <- dplyr::filter(rm_file, !is.na(qry_start), !is.na(qry_end))
  # rm_file <- dplyr::mutate(rm_file, qry_width = as.integer(qry_end -
  #                                                            qry_start + 1L))
  # nrow_after_filtering <- nrow(rm_file)
  # if ((nrow_before_filtering - nrow_after_filtering) > 0)
  #   message((nrow_before_filtering - nrow_after_filtering) +
  #             1, " out of ", nrow_before_filtering, " rows ~ ",
  #           round(((nrow_before_filtering - nrow_after_filtering) +
  #                    1) / nrow_before_filtering, 3), "% were removed from the imported RepeatMasker file, ",
  #           "because they contained 'NA' values in either 'qry_start' or 'qry_end'.")
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

### V1.1 new code below 

# first parse per hit ID (replace OneCodeToFindThemAll)
read_rm_custom(opt$dotout) %>% group_by(ID) %>%
  summarize(
    qry_id = unique(qry_id),
    hit_qry_start = min(as.integer(qry_start)),
    hit_qry_end = min(as.integer(qry_end)),
    hit_strand = paste0(unique(strand), collapse =  ""),
    frg_targ_start = paste0(unique(target_start), collapse =  ","),
    frg_targ_end = paste0(unique(target_end), collapse =  ","),
    hit_name = repeat_id[which.max(as.integer(sw_score))],
    fragments_names = paste0(repeat_id, collapse =  ","),
    hit_matching_class = matching_class[which.max(as.integer(sw_score))],
    is_confused = ifelse(length(unique(repeat_id)) > 1, "yes", "no"),
    final_hit_name = ifelse(is_confused == "no", hit_name, paste(hit_name, "(x)", sep = "")),
    hit_perc_div = sum(as.numeric(perc_div) * (as.numeric(qry_end) - as.numeric(qry_start)))/sum(as.numeric(qry_end) - as.numeric(qry_start)),
    fragments = n()
  ) -> dotout_parsed_per_ID

# L1 5' INVERSIONS (Twin Priming) see: https://genome.cshlp.org/content/11/12/2059.long & https://link.springer.com/article/10.1186/1759-8753-1-7
# assess L1 inversions: must be a single L1 hit, 2 fragments, with C+ strand pattern.
# positive strand L1 with 5P inversions have the frag_targ_start of the first fragments (C) smaller than the second fragment (+).
dotout_parsed_per_ID$is_L15PINV <- ifelse(dotout_parsed_per_ID$hit_matching_class == "LINE/L1" & dotout_parsed_per_ID$hit_strand == "C+", "yes", "no")
dotout_parsed_per_ID$final_hit_strand <- ifelse(
  dotout_parsed_per_ID$hit_matching_class == "LINE/L1" & dotout_parsed_per_ID$hit_strand == "C+",
  sapply(dotout_parsed_per_ID$frg_targ_start, function(x) {
    parts <- as.numeric(unlist(str_split(x, ",")))
    ifelse(parts[1] < parts[2], "+", "-") # assess orientation of the L1 5' inversion   
  }),
  dotout_parsed_per_ID$hit_strand         # otherwise report hit orientation
)                                                                 

# SVA VNTR polymorphism
# Assess whether the variant annotated as SVA is in fact only the VNTR region on an SVA
# create a table with the VNTR coordinate of each SVA family
SVA<-data.frame(family = c("SVA_A", "SVA_B", "SVA_C", "SVA_D", "SVA_E", "SVA_F"),
           VNTR_start = c(436, 431, 432, 432, 428, 435),
           VNTR_end = c(855, 867, 851, 689, 864, 857))
# check the hits
dotout_parsed_per_ID <- dotout_parsed_per_ID %>%
  left_join(SVA, by = c("hit_name" = "family")) %>%
  mutate(
    in_vntr = case_when(
      hit_matching_class == "Retroposon/SVA" & 
        fragments == 1 & 
        !is.na(VNTR_start) & 
        !is.na(VNTR_end) &
        suppressWarnings(as.numeric(frg_targ_start) - 1) > VNTR_start - 1 & 
        suppressWarnings(as.numeric(frg_targ_end)) < VNTR_end ~ "yes",
      TRUE ~ "no"
    )
  ) %>%
  select(-VNTR_start, -VNTR_end)  # remove the joined columns
# update the finale repeat name
dotout_parsed_per_ID$final_hit_name <- ifelse(dotout_parsed_per_ID$hit_matching_class == "Retroposon/SVA" & dotout_parsed_per_ID$in_vntr == "yes", 
                                              paste(dotout_parsed_per_ID$final_hit_name, "(VNTR_only)", sep = ""),
                                              dotout_parsed_per_ID$final_hit_name)
dotout_parsed_per_ID$hit_matching_class <- ifelse(dotout_parsed_per_ID$hit_matching_class == "Retroposon/SVA" & dotout_parsed_per_ID$in_vntr == "yes", 
                                              "Simple_repeat", # Switch from Retroposon/SVA to Simple_repeat to avoid accidental count of these as SVA MEI
                                              dotout_parsed_per_ID$hit_matching_class)

# compile RM hit stats for each SV
dotout_parsed_per_ID %>%
  mutate(match_len = hit_qry_end - hit_qry_start - 1) %>%
  arrange(hit_qry_start) %>%
  group_by(qry_id) %>%
  summarise(
    start = first(hit_qry_start), stop = max(hit_qry_end),
    match_lengths = paste0(match_len, collapse = ","),
    fragmts = paste0(fragments, collapse = ","),
    repeat_ids = paste0(final_hit_name, collapse = ","),
    matching_classes = paste0(hit_matching_class, collapse = ","),
    strands = paste0(hit_strand, collapse = ","),
    RM_id = paste0(ID, collapse = ","),
    L1_5PINV = paste0(ID[is_L15PINV == "yes"], collapse = ","),
    n_hits = n()
  ) -> rep_mask 

### legacy code below 
# we pull the VCF to make an annotation table to use later with bcftools.
vcf <- read.vcfR(opt$vcf)
vcf_df <- tibble(
  CHROM = getCHROM(vcf),
  POS = getPOS(vcf),
  REF = getREF(vcf),
  ALT = getALT(vcf),
  qry_length = abs(str_length(getALT(vcf)) - str_length(getREF(vcf))),
  qry_id = getID(vcf)
)
# create the annotation table for each SV with default
annot <- left_join(vcf_df, rep_mask, by = "qry_id") %>%
  replace_na(list(matching_classes = "None",
                  repeat_ids = "None",
                  n_hits = 0,
                  fragmts = "0",
                  match_lengths = "0",
                  strands = "None",
                  RM_id = "None",
                  L1_5PINV = "None")) %>%
    select(-c(qry_length)) %>%
    arrange(CHROM, POS, qry_id) %>%
    select(CHROM, POS, qry_id, REF, ALT, n_hits, fragmts, match_lengths, repeat_ids, matching_classes, strands, RM_id, L1_5PINV)

write_tsv(annot, file = opt$annotation, col_names = F)
print(colnames(annot))
