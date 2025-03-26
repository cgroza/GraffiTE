#! /bin/bash

# USAGE: ./prepTSD.sh <REF_GENOME> <WINDOW_SIZE>

# variable list
VCF="genotypes_repmasked_filtered.vcf" # filtered vcf with repeatmasker
REF=$1 # ref genome
WIN=$2 # windows size in flanking to search TSD
MSK="indels.fa.masked" # masked SVs from repeatmasker
OUT_VCF="pangenie.vcf"

###################################################
# Step 1: extract flanking of each retained TE SV #
###################################################
echo "extracting flanking..."

# get contig length for bedtools
grep "contig=" ${VCF} | sed 's/\#\#contig=<ID=//g;s/,length=/\t/g;s/>//g' > gLength.txt

# create a bed with vcf entries
grep -v '#' ${VCF} | \
 grep 'n_hits=1;\|n_hits=2;' | \
 grep -v 'mam_filter_2=VNTR_ONLY' | \
 awk '/n_hits=1/ && length($4) < length($5) {print $1"\t"$2"\t"($2)+1"\t"$3; next} /n_hits=1/ && length($4) > length($5) {print $1"\t"$2"\t"($2+length($4))"\t"$3; next} /n_hits=2/ && length($4) < length($5) && /5P_INV/ {print $1"\t"$2"\t"($2)+1"\t"$3; next} /n_hits=2/ && length($4) > length($5) && /5P_INV/ {print $1"\t"$2"\t"($2+length($4))"\t"$3}' > oneHit_SV_coordinates.bed

# extend +/- ${WIN} bp in two entries per SV
cat <(bedtools slop -i oneHit_SV_coordinates.bed -g gLength.txt -l 30 -r 0 | awk '{print $0"__L"}') \
<(bedtools slop -i oneHit_SV_coordinates.bed -g gLength.txt -l 0 -r 30 | awk '{print $0"__R"}') | \
sort -k1,1 -k2,2n -k3,3n | awk '/__L/ {print $1":"$2"-"($2+30); next} /__R/ {print $1":"($3-30)"-"$3}' > oneHit_SV_coordinates_win.regions
# extract fasta from flanking
samtools faidx -r oneHit_SV_coordinates_win.regions -o flanking_sequences.fasta ${REF}

##################################################
# Step 2: extract 5' and 3' of each masked TE SV #
##################################################
echo "extracting SVs' 5' and 3' ends..."

# filter the indels.fa.masked to keep only single RM hits (1 TE per SV)
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(cut -f 4 oneHit_SV_coordinates.bed) ${MSK} > oneHit_indels.fa.masked
# linearize fasta, then trim and split in two seq (L and R)
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' oneHit_indels.fa.masked | \
awk -v len=${WIN} -F '\t' '{x=len;L=length($2);printf("%s\n%s\n%s\n%s\n",$1"__L",(L<=x?$2:substr($2,2,x+1)),$1"__R",(L<=x?$2:substr($2,1+L-x,x)));}' > SV_sequences_L_R_trimmed_WIN.fa
# export the list of SV to search TSD for next process parallelization
grep '>' SV_sequences_L_R_trimmed_WIN.fa | sed 's/>//g;s/__/\t/g' | cut -f 1 | sort | uniq > indels.txt
