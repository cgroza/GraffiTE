#!/bin/bash

VCF=$1
OUT_VCF=$2
FASTA_LIB=$3
MAM=$4

FASTA_FILE=indels.fa
bcftools view -H --types indels --include 'ILEN>0' ${VCF} | awk '{print(sprintf(">%s\n%s", $3, $5))}' >> ${FASTA_FILE}
bcftools view -H --types indels --include 'ILEN<0' ${VCF} | awk '{print(sprintf(">%s\n%s", $3, $4))}' >> ${FASTA_FILE}

# verify that indels.fa ids are not longer than 50 characters
if grep ">" indels.fa | awk 'length > 50' | grep -q .; then
    echo "variant IDs must be no greater than 50 characters"
    exit 1
fi

mkdir repeatmasker_dir
REPMASK_DIR=repeatmasker_dir

# detect number of cores allocated to this nextflow process
repmask_cores="$(($(nproc)/4))"

# guard against number of cores smaller than 1
if [ "${repmask_cores}" -lt "1" ]; then
    repmask_cores="1"
fi

# run RepeatMasker
RepeatMasker -lib ${FASTA_LIB} -s -dir ${REPMASK_DIR} -pa ${repmask_cores} ${FASTA_FILE}

REPMASK_OUT=${REPMASK_DIR}/$(basename ${FASTA_FILE}).out

# run ULTRA to detect tandem repeat the RM may have annotated as TE

mkdir ultra_temp
ULTRA_DIR=ultra_temp
ULTRA_OUT=ultra_out.bed
ultra --bed -t $(nproc) -o ${ULTRA_DIR}/ultra_out ${FASTA_FILE}
# ULTRA's -o naming varies: single-format mode may write the literal filename;
# multi-format mode treats -o as a prefix and appends .bed/.tsv.
if [ ! -f ${ULTRA_DIR}/ultra_out.bed ] && [ -f ${ULTRA_DIR}/ultra_out ]; then
    mv ${ULTRA_DIR}/ultra_out ${ULTRA_DIR}/ultra_out.bed
fi
cp ${ULTRA_DIR}/ultra_out.bed ./ultra_out.bed
# Non-redundant ULTRA-annotated bases per SV (each SV is a "chr" in the bed).
# bedtools merge collapses overlapping intervals within each SV; we then sum
# the merged interval widths to get one row per SV: <SV_id>\t<non_redundant_bp>
sort -k1,1 -k2,2n ultra_out.bed | bedtools merge -i - | \
  awk 'BEGIN{OFS="\t"} {sum[$1]+=$3-$2} END{for (i in sum) print i, sum[i]}' | \
  sort -k1,1 > ultra_out.span

ANNOT_FILE=vcf_annotation

annotate_vcf.R --dotout ${REPMASK_OUT} --vcf ${VCF} --annotation ${ANNOT_FILE}_1 # v1.1 changed to pure RM .out (remove OneCode)

# calculate the total span of TEs on the SV without overlap
echo "compute repeat proportion for each SVs..."
samtools faidx indels.fa
awk '{print $1"\t"$2}' indels.fa.fai > indels.length
grep -v 'Simple_repeat\|Low_complexity' ${REPMASK_OUT} | awk 'NR > 3 {print $5"\t"$6-1"\t"$7"\t"$10}' | sort -k1,1 -k2,2n -k3,3n | bedtools merge > merge.bed # add -1 to start to meet .bed format
rm -rf span &> /dev/null # clean in case there is a "span" file already
rm -rf ${ANNOT_FILE}.gz &> /dev/null # clean in case there was a ${ANNOT_FILE}.gz file already
# join allele length with TE span
# input 1: sum of TE length group by SV ID
# input 2: SV length
# output:  id, sum of TE length, SV length, TE span; sorted by SV ID
join -11 -21 \
<(awk '{sum[$1] += $3-$2} END {for (i in sum) print i"\t"sum[i]}' merge.bed | sort -k1,1) \
<(sort -k1,1 indels.length) | \
awk '{print $1"\t"$2"\t"$3"\t"($2/$3)}' > span

# merge with ${ANNOT_FILE}_1
join -13 -21 -a1 <(sort -k3,3 ${ANNOT_FILE}_1)  <(sort -k1,1 span) | sed 's/ /\t/g' # this is to see how it looks
join -13 -21 -a1 <(sort -k3,3 ${ANNOT_FILE}_1)  <(sort -k1,1 span) | sed 's/ /\t/g' | \
 awk '{print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$16}' | \
  awk '{if (NF == 13) {print $0"\t0\t0"} else {print $0}}' > vcf_annotation.tmp
# left-join ULTRA non-redundant span onto the annotation (key = SV ID, col 3)
# unmatched SVs (no tandem repeat found by ULTRA) get ULTRA_TR=0
# combine ULTRA non-redundant bp with variant sequence length to get span ratio
# (capped at 1 in case ULTRA annotates 1bp more than the ALT/REF sequence length)
join -11 -21 <(sort -k1,1 ultra_out.span) <(sort -k1,1 indels.length) | \
 awk 'BEGIN{OFS="\t"} {r=$2/$3; if(r>1)r=1; print $1, $2, r}' | \
 sort -k1,1 > ultra_out.stats
# left-join ULTRA stats onto the annotation (key = SV ID, col 3)
# unmatched SVs (no tandem repeat found by ULTRA) get ULTRA_TR=0, ULTRA_TR_span=0
join -13 -21 -a1 <(sort -k3,3 vcf_annotation.tmp) ultra_out.stats | sed 's/ /\t/g' | \
 awk 'BEGIN{OFS="\t"} {u_bp=(NF>=17)?$16:0; u_sp=(NF>=17)?$17:0; print $2,$3,$1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,u_bp,u_sp}' | \
  sort -k1,1 -k2,2n > vcf_annotation #${ANNOT_FILE}
# copy for dev
cp vcf_annotation vcf_annotation.bak.txt

# if --mammal if set, search for L1 5' inversion (Twin Priming and similar) and if SVA hits are within VNTR only (non retrotransposition polymorphism)
if [[ ${MAM} == "MAM" ]]
then
    echo "[INFO] Mammalian filters (--mammal) option has been discontinued since version 1.1"
    echo "       L1 5' inversions and SVA VNTR only polymorphisms are now systematically reported"
fi

echo "writing vcf..."
#head ${ANNOT_FILE}
bgzip vcf_annotation #${ANNOT_FILE}
tabix -s1 -b2 -e2 vcf_annotation.gz #${ANNOT_FILE}.gz

HDR_FILE=hdr_file

echo -e '##INFO=<ID=n_hits,Number=1,Type=Integer,Description="Number of repeats found in insertion">' >> ${HDR_FILE}
echo -e '##INFO=<ID=match_lengths,Number=.,Type=Integer,Description="Insertion lengths spanned by each repeat">' >> ${HDR_FILE}
echo -e '##INFO=<ID=repeat_ids,Number=.,Type=String,Description="Repeat family IDs">' >> ${HDR_FILE}
echo -e '##INFO=<ID=matching_classes,Number=.,Type=String,Description="Repeat class names">' >> ${HDR_FILE}
echo -e '##INFO=<ID=fragmts,Number=.,Type=Integer,Description="Number of fragments merged into one by one code">' >> ${HDR_FILE}
echo -e '##INFO=<ID=RM_hit_strands,Number=.,Type=String,Description="RepeatMasker hit strands">' >> ${HDR_FILE}
echo -e '##INFO=<ID=RM_hit_IDs,Number=.,Type=String,Description="RepeatMasker hit IDs">' >> ${HDR_FILE}
echo -e '##INFO=<ID=total_match_length,Number=1,Type=Integer,Description="Insertion length spanned by repeats">' >> ${HDR_FILE}
echo -e '##INFO=<ID=total_match_span,Number=1,Type=Float,Description="Insertion span spanned by repeats">' >> ${HDR_FILE}
echo -e '##INFO=<ID=L1_5PINV,Number=.,Type=String,Description="RM hit ID in this SV flagged as LINE1 with 5-prime inversion">' >> ${HDR_FILE}
echo -e '##INFO=<ID=ULTRA_TR,Number=1,Type=Integer,Description="Non-redundant bases of tandem repeats annotated by ULTRA within the insertion (bedtools-merged)">' >> ${HDR_FILE}
echo -e '##INFO=<ID=ULTRA_TR_span,Number=1,Type=Float,Description="Fraction of the variant sequence spanned by ULTRA tandem repeats (ULTRA_TR / variant length, capped at 1)">' >> ${HDR_FILE}
echo -e '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> ${HDR_FILE}

cat <(bcftools view -h ${VCF}) <(bcftools view -H ${VCF} | sort -k1,1 -k2,2n) > genotypes.sorted.vcf
bcftools annotate -a ${ANNOT_FILE}.gz -h ${HDR_FILE} \
-c CHROM,POS,~ID,REF,ALT,INFO/n_hits,INFO/fragmts,INFO/match_lengths,INFO/repeat_ids,INFO/matching_classes,INFO/RM_hit_strands,INFO/RM_hit_IDs,INFO/L1_5PINV,INFO/total_match_length,INFO/total_match_span,INFO/ULTRA_TR,INFO/ULTRA_TR_span genotypes.sorted.vcf | \
bcftools view -Oz -o ${OUT_VCF}
