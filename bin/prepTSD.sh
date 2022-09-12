#! /bin/bash

# USAGE: ./findTSD.sh <REF_GENOME>

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
## USE THIS FOR NOW TO FILTER REAL 1 hits

awk 'NR > 3' indels.fa.onecode.out | cut -f 5 | sort | uniq -d > remove

grep -v '#' ${VCF} | \
 grep -vwf remove | \
 awk '/INS/ {print $1"\t"$2"\t"($2)+1"\t"$3; next} /DEL/ {print $1"\t"$2"\t"($2+length($4))"\t"$3}' > oneHit_SV_coordinates.bed

## REVERSE TO THIS WHEN CRISTIAN HAS FIXED THE RM FILTER SCRIPT

# grep -v '#' ${VCF} | \
#  grep 'n_hits=1;' | \
#  awk '/INS/ {print $1"\t"$2"\t"($2)+1"\t"$3; next} /DEL/ {print $1"\t"$2"\t"($2+length($4))"\t"$3}' > oneHit_SV_coordinates.bed

# extend +/- ${WIN} bp in two entries per SV
cat <(bedtools slop -i oneHit_SV_coordinates.bed -g gLength.txt -l 30 -r 0 | awk '{print $0"__L"}') \
<(bedtools slop -i oneHit_SV_coordinates.bed -g gLength.txt -l 0 -r 30 | awk '{print $0"__R"}') | \
sort -k1,1 -k2,2n -k3,3n | awk '/__L/ {print $1"\t"$2"\t"($2+30)"\t"$4; next} /__R/ {print $1"\t"($3-30)"\t"$3"\t"$4}' > oneHit_SV_coordinates_win.bed
# extract fasta from flanking
bedtools getfasta -fi ${REF} -bed oneHit_SV_coordinates_win.bed -name > flanking_sequences.fasta

##################################################
# Step 2: extract 5' and 3' of each masked TE SV #
##################################################
echo "extracting SVs' 5' and 3' ends..."

# filter the indels.fa.masked to keep only single RM hits (1 TE per SV)
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' <(cut -f 4 oneHit_SV_coordinates.bed) ${MSK} > oneHit_indels.fa.masked
# linearize fasta, then trim and split in two seq (L and R)
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' oneHit_indels.fa.masked | \
awk -v len=${WIN} -F '\t' '{x=len;L=length($2);printf("%s\n%s\n%s\n%s\n",$1"__L",(L<=x?$2:substr($2,2,x+1)),$1"__R",(L<=x?$2:substr($2,1+L-x,x)));}' > SV_sequences_L_R_trimmed_WIN.fa