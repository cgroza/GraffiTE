#! /bin/bash

# USAGE: ./prepTSD.sh <REF_GENOME> <WINDOW_SIZE> [THREADS]
#
# Prepares the two FASTA files TSD_Match_v2.sh compares: the reference flanks
# of every indel and the two ends of every indel sequence.

set -euo pipefail

VCF="genotypes_repmasked_filtered.vcf" # filtered vcf with repeatmasker
REF=$1 # ref genome
WIN=$2 # windows size in flanking to search TSD
THREADS=${3:-1}
FASTA_FILE=indels.fa

n_records=$(bcftools view -H ${VCF} | wc -l)
if [[ ${n_records} -eq 0 ]]; then
    echo "prepTSD.sh: ${VCF} has no records, nothing to prepare"
    : > flanking_sequences.fasta
    : > ${FASTA_FILE}
    : > SV_sequences_L_R_trimmed_WIN.fa
    : > indels.txt
    exit 0
fi

###################################################
# Step 1: extract flanking of each retained TE SV #
###################################################
echo "extracting flanking..."

# htslib reads plain and BGZF FASTA and refuses gzip. A gzip reference used to
# reach bedtools getfasta here, which failed the same way, and the failure was
# silent: an empty flank file, and a search of each SV's two ends against each
# other. Same re-compression as concat_repeatmask.
if [[ "${REF}" == *.gz ]] && ! (file -L "${REF}" | grep -q "BGZF"); then
    BGZF="bgzf_$(basename "${REF}")"
    echo "re-compressing ${REF} to ${BGZF}"
    gzip -dc "${REF}" | bgzip -@ "${THREADS}" -c > "${BGZF}"
    REF=${BGZF}
fi
samtools faidx "${REF}"

tsd_flanks.py --vcf ${VCF} --reference "${REF}" --window ${WIN} \
    --out flanking_sequences.fasta

##################################################
# Step 2: extract 5' and 3' of each masked TE SV #
##################################################
echo "extracting SVs' 5' and 3' ends..."

bcftools view -H --types indels --include 'ILEN>0' ${VCF} | awk '{print(sprintf(">%s\n%s", $3, $5))}' > ${FASTA_FILE}
bcftools view -H --types indels --include 'ILEN<0' ${VCF} | awk '{print(sprintf(">%s\n%s", $3, $4))}' >> ${FASTA_FILE}

# linearize fasta, then trim and split in two seq (L and R)
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' indels.fa | \
awk -v len=${WIN} -F '\t' '{x=len;L=length($2);printf("%s\n%s\n%s\n%s\n",$1"__L",(L<=x?$2:substr($2,2,x)),$1"__R",(L<=x?$2:substr($2,1+L-x,x)));}' > SV_sequences_L_R_trimmed_WIN.fa
# export the list of SV to search TSD for next process parallelization
grep '>' SV_sequences_L_R_trimmed_WIN.fa | sed 's/>//g;s/__/\t/g' | cut -f 1 | sort | uniq > indels.txt

n_sv=$(wc -l < indels.txt)
n_flank=$(grep -c '^>' flanking_sequences.fasta)
echo "${n_sv} indels, ${n_flank} flank sequences"
if [[ $((n_sv * 2)) -ne ${n_flank} ]]; then
    echo "prepTSD.sh: expected two flanks per indel" >&2
    exit 1
fi
