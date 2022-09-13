#! /bin/bash

# USAGE: ./findTSD.sh <REF_GENOME>

# variable list
VCF="genotypes_repmasked_filtered.vcf" # filtered vcf with repeatmasker
REF=$1 # ref genome
WIN=$2 # windows size in flanking to search TSD
MSK="indels.fa.masked" # masked SVs from repeatmasker
OUT_VCF="pangenie.vcf"

#######################
# Step 3: search TSDs #
#######################
echo "searching TSDs..."
TSD_Match_v2.sh SV_sequences_L_R_trimmed_WIN.fa flanking_sequences.fasta

# ########################
# # Step 4: annotate VCF #
# ########################

# # create TSD_annotation from TSD_summary.txt
# join -13 -21 <(grep -v "#" ${VCF} | cut -f 1-3 | sort -k3,3) <(grep 'PASS' TSD_summary.txt | awk '{print $1"\t"$(NF-2)","$(NF-1)}' | sort -k1,1) | awk '{print $2"\t"$3"\t"$4}' | sort -k1,1 -k2,2n > TSD_annotation
# # add new header line to describe it
# HDR_FILE=$(mktemp)
# echo -e '##INFO=<ID=TSD,Number=1,Type=String,Description="Target site duplication sequence passing filters">' >> ${HDR_FILE}
# # bgzip it and annotate
# TSD_FILE=TSD_annotation 
# bgzip ${TSD_FILE}
# tabix -s1 -b2 -e2 ${TSD_FILE}.gz 
# bcftools annotate -a ${TSD_FILE}.gz -h ${HDR_FILE} -c CHROM,POS,INFO/TSD ${VCF} | bcftools view > ${OUT_VCF}


# # clean 
# rm oneHit_SV_coordinates.bed
# rm oneHit_SV_coordinates_win.bed