#!/bin/bash

VCF=$1
OUT_VCF=$2
FASTA_LIB=$3

FASTA_FILE=indels.fa
cat ${VCF} | grep -v "#" | grep "SVTYPE=INS" | awk '{print(sprintf(">%s\n%s", $3, $5))}' > ${FASTA_FILE}
cat ${VCF} | grep -v "#" | grep "SVTYPE=DEL" | awk '{print(sprintf(">%s\n%s", $3, $4))}' >> ${FASTA_FILE}

mkdir repeatmasker_dir
REPMASK_DIR=repeatmasker_dir

RepeatMasker -nolow -lib ${FASTA_LIB} -s -dir ${REPMASK_DIR} -pa 16 ${FASTA_FILE}

REPMASK_OUT=${REPMASK_DIR}/$(basename ${FASTA_FILE}).out

ANNOT_FILE=vcf_annotation

annotate_vcf.R --dotout ${REPMASK_OUT} --vcf ${VCF} --annotation ${ANNOT_FILE}

bgzip ${ANNOT_FILE}
tabix -s1 -b2 -e2 ${ANNOT_FILE}.gz

HDR_FILE=$(mktemp)

echo -e '##INFO=<ID=n_hits,Number=1,Type=Float,Description="Number of repeats found in insertion">' >> ${HDR_FILE}
echo -e '##INFO=<ID=match_lengths,Number=1,Type=String,Description="Insertion lengths spanned by each repeat">' >> ${HDR_FILE}
echo -e '##INFO=<ID=repeat_ids,Number=1,Type=String,Description="Repeat family IDs">' >> ${HDR_FILE}
echo -e '##INFO=<ID=matching_classes,Number=1,Type=String,Description="Repeat class names">' >> ${HDR_FILE}
echo -e '##INFO=<ID=total_match_length,Number=1,Type=Float,Description="Insertion length spanned by repeats">' >> ${HDR_FILE}
echo -e '##INFO=<ID=total_match_span,Number=1,Type=Float,Description="Insertion span spanned by repeats">' >> ${HDR_FILE}
echo -e '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> ${HDR_FILE}

bcftools annotate -a ${ANNOT_FILE}.gz -h ${HDR_FILE} \
         -c CHROM,POS,INFO/n_hits,INFO/match_lengths,INFO/repeat_ids,INFO/matching_classes,INFO/total_match_length,INFO/total_match_span $VCF | \
    bcftools view -Oz -o ${OUT_VCF}

