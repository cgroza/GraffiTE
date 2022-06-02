#!/bin/bash

VCF=$1
OUT_VCF=$2
FASTA_LIB=$3

FASTA_FILE=$(mktemp)
cat ${VCF} | grep -v "#" | grep "svim_asm.INS" | awk '{print(sprintf(">%s\n%s", $3, $5))}' > ${FASTA_FILE}
cat ${VCF} | grep -v "#" | grep "svim_asm.DEL" | awk '{print(sprintf(">%s\n%s", $3, $4))}' >> ${FASTA_FILE}

REPMASK_DIR=$(mktemp -d)

RepeatMasker -lib ${FASTA_LIB} -s -dir ${REPMASK_DIR} -pa 16 ${FASTA_FILE}

REPMASK_OUT=${REPMASK_DIR}/$(basename ${FASTA_FILE}).out

ANNOT_FILE=$(mktemp)

Rscript - <<SCRIPT
library(biomartr)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(vcfR)

repmask_out <- "$REPMASK_OUT"
annot_vcf <- "$VCF"
annot_file <- "$ANNOT_FILE"

rep_mask <- read_rm(repmask_out) %>%
  mutate(match_len = qry_end - qry_start) %>%
  group_by(qry_id) %>%
  filter(match_len == max(match_len)) %>%
  filter(sw_score == max(sw_score)) %>%
  select(matching_class, repeat_id, match_len)

vcf <- read.vcfR(annot_vcf)
vcf_df <- tibble(CHROM = getCHROM(vcf),
                 POS = getPOS(vcf),
                 qry_length = abs(str_length(getALT(vcf)) - str_length(getREF(vcf))),
                 qry_id = getID(vcf))

annot <- left_join(vcf_df, rep_mask, by = "qry_id") %>%
  replace_na(list(matching_class = "None", repeat_id = "None", match_len = 0)) %>%
  mutate(match_span = match_len / qry_length) %>%
  select(-c(qry_id, qry_length, match_len)) %>%
  arrange(CHROM, POS)

write_tsv(annot, file = annot_file, col_names = F)
print(colnames(annot))
SCRIPT

bgzip ${ANNOT_FILE}
tabix -s1 -b2 -e2 ${ANNOT_FILE}.gz

HDR_FILE=$(mktemp)

echo -e '##INFO=<ID=repFamily,Number=1,Type=String,Description="Repeat family">' > ${HDR_FILE}
echo -e '##INFO=<ID=repClass,Number=1,Type=String,Description="Repeat name">' >> ${HDR_FILE}
echo -e '##INFO=<ID=match_span,Number=1,Type=Float,Description="Insertion length spaned by repeat">' >> ${HDR_FILE}

bcftools annotate -a ${ANNOT_FILE}.gz -h ${HDR_FILE} \
         -c CHROM,POS,INFO/repClass,INFO/repFamily,INFO/match_span $VCF | \
    bcftools view -Oz -o ${OUT_VCF}

