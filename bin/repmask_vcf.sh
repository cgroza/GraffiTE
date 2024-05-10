#!/bin/bash

VCF=$1
OUT_VCF=$2
FASTA_LIB=$3
MAM=$4

FASTA_FILE=indels.fa
bcftools view --types indels --include 'ILEN>0' ${VCF} | grep -v "#" | awk '{print(sprintf(">%s\n%s", $3, $5))}' >> ${FASTA_FILE}
bcftools view --types indels --include 'ILEN<0' ${VCF} | grep -v "#" | awk '{print(sprintf(">%s\n%s", $3, $4))}' >> ${FASTA_FILE}

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
RepeatMasker -nolow -lib ${FASTA_LIB} -s -dir ${REPMASK_DIR} -pa $(nproc) ${FASTA_FILE}

REPMASK_OUT=${REPMASK_DIR}/$(basename ${FASTA_FILE}).out
REPMASK_ONECODE_OUT=${REPMASK_DIR}/$(basename ${FASTA_FILE}).onecode.out

### Step 1 make the dictionary file if LTRs are split between _LTR and _I
echo "Building onecode LTR dictionary..."
build_dictionary.pl --rm ${REPMASK_OUT} --unknown > ${REPMASK_DIR}/OneCode_LTR.dic

### Step 2 run OneCodeToFindThemAll
echo "Running onecode..."
one_code_to_find_them_all.pl --rm ${REPMASK_OUT} --unknown --ltr ${REPMASK_DIR}/OneCode_LTR.dic > ${REPMASK_DIR}/onecode.log 2>&1

### Step 3 parse the output
# remove output if previous run
echo "Concatenate outputs..."
rm ${REPMASK_DIR}/ALL.elem_sorted.csv 2>/dev/null
# loop over each individual output and concatenate in a single file (can't cat, often there are too many files)
for i in ${REPMASK_DIR}/*.elem_sorted.csv
do
	cat $i >> ${REPMASK_DIR}/ALL.elem_sorted.csv
done

# parse the concatenated outputs and save the raw to recover stitched copies
echo "Parse outputs..."
cat <(head -n 3 ${REPMASK_OUT} | sed 's/ID/ID\tfragmts/g') <(grep '###' ${REPMASK_DIR}/ALL.elem_sorted.csv | cut -f 1-16 | sed 's/###//g') > ${REPMASK_ONECODE_OUT}
mv ${REPMASK_DIR}/ALL.elem_sorted.csv ${REPMASK_DIR}/ALL.onecode.elem_sorted.bak

# clean-up the mess (BE SURE YOU HAVE RENAMED THE OUTPUT)
echo "Cleanup..."
find ${REPMASK_DIR} -name "*.copynumber.csv" -delete
find ${REPMASK_DIR} -name "*.elem_sorted.csv" -delete
find ${REPMASK_DIR} -name "*.ltr.csv" -delete
find ${REPMASK_DIR} -name "*.transposons.csv" -delete

ANNOT_FILE=vcf_annotation

annotate_vcf.R --dotout ${REPMASK_ONECODE_OUT} --vcf ${VCF} --annotation ${ANNOT_FILE}_1

# calculate the total span of TEs on the SV without overlap
echo "compute repeat proportion for each SVs..."
samtools faidx indels.fa
awk '{print $1"\t"$2}' indels.fa.fai > indels.length
awk 'NR > 3 {print $5"\t"$6"\t"$7"\t"$10}' ${REPMASK_ONECODE_OUT} | bedtools merge > merge.bed
RMQUERIES=$(awk 'NR > 3 {print $5}' ${REPMASK_ONECODE_OUT} | sort | uniq)
rm -rf span &> /dev/null # clean in case there is a "span" file already
rm -rf ${ANNOT_FILE}.gz &> /dev/null # clean in case there was a ${ANNOT_FILE}.gz file already
for i in ${RMQUERIES}
do
paste -d "\t" <(echo -e "${i}") <(grep -w "${i}" merge.bed | awk '{print ($3-$2)}' | paste -sd+ | bc) <(grep -w "${i}" indels.length) | awk '{print $1"\t"$2"\t"$4"\t"($2/$4)}' >> span
done
# merge with ${ANNOT_FILE}_1
join -13 -21 -a1 <(sort -k3,3 ${ANNOT_FILE}_1)  <(sort -k1,1 span) | sed 's/ /\t/g' | \
 awk '{print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$15}' | \
 awk '{if (NF == 10) {print $0"\t0\t0"} else {print $0}}' | \
 sort -k1,1 -k2,2n > ${ANNOT_FILE}

# if --mammal if set, search for L1 5' inversion (Twin Priming and similar) and if SVA hits are within VNTR only (non retrotransposition polymorphism)
if [[ ${MAM} == "MAM" ]]
then
    echo "Mammalian filters ON. Filtering..."
    # FILTER 1: L1 Twin Priming and similar
    TwP=$(awk '{if ($6 == 2 && $11 == "C,+" && $10~/LINE/) {names = split($10,n,",",seps); ids = split($12, i,",",seps); if (names n[1] == names n[2] && ids i[1] == ids i[2]) {print $0"\tTwP"}}}' vcf_annotation | awk '{print $3}')
    rm TwP.txt &> /dev/null
    for i in $TwP 
    do
    grep -w $(grep -w "${i}" <(grep -v "#" genotypes.vcf | awk '{print $1"_"$2"\t"$3}') | cut -f 2) repeatmasker_dir/indels.fa.onecode.out | \
    awk 'getline second {line=split(second, a, "\t", sep); print $12"\t"$13"\t"$14"\t"a[12]"\t"a[13]"\t"a[14]}' | \
    awk -v coord=${i} '{if ($2<$5) {print coord"\t5P_INV:plus"} else {print coord"\t5P_INV:minus"}}' >> TwP.txt
    done

    # FILTER 2: SVA VNTR
    # get coordinates of each putative TE in TE (single TE hits)
    awk '$6 == 1 && $14 > 0.9 {print $1"_"$2"\t"$1"\t"$2"\t"($2+1)"\t"$8"\t"$10}' vcf_annotation | grep 'SVA' | sort -k1,1 > SVA_candidates
    join -17 -25 <(join -11 -21 SVA_candidates <(awk 'NR > 3 && !/^#/ {print $1"_"$2"\t"$3}' genotypes.vcf | \
     sort -k1,1) | sort -k7,7) <(sort -k5,5 repeatmasker_dir/indels.fa.onecode.out) | \
     sed 's/ /\t/g' | \
     awk '{if ($15 == "+") {print $16"\t"$18"\t"$19"\t"$1} else {print $16"\t"$20"\t"$19"\t"$1}}' > SVA_candidates.bed
    # create bed with SVAs' VNTR position according to TRF on DFAM3.6 models (2022)
    rm SVA_VNTR.bed &> /dev/null
    echo -e "SVA_A\t436\t855\t+" >> SVA_VNTR.bed
    echo -e "SVA_B\t431\t867\t+" >> SVA_VNTR.bed
    echo -e "SVA_C\t432\t851\t+" >> SVA_VNTR.bed
    echo -e "SVA_D\t432\t689\t+" >> SVA_VNTR.bed
    echo -e "SVA_E\t428\t864\t+" >> SVA_VNTR.bed
    echo -e "SVA_F\t435\t857\t+" >> SVA_VNTR.bed
    bedtools intersect -wao -a <(sort -k1,1 -k2,2n SVA_candidates.bed) -b SVA_VNTR.bed | awk '{if ($9/($3-$2) > 0.9) {print $4"\tVNTR_ONLY:"$5":"$2":"$3}}' | sort -k1,1 > SVA_VNTR.txt

    join -a1 -13 -21 <(sort -k3,3 vcf_annotation) <(sort -k1,1 TwP.txt) | sed 's/ /\t/g' | awk '{if (NF == 14) {print $0"\tNone"} else {print $0}}' > vcf_annotation.temp
    join -a1 -11 -21 <(sort -k1,1 vcf_annotation.temp) SVA_VNTR.txt | sed 's/ /\t/g' | awk '{if (NF == 15) {print $0"\tNone"} else {print $0}}' | \
     awk '{print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16}' | sort -k1,1 -k2,2n > ${ANNOT_FILE}

    echo "writing vcf..."
    bgzip ${ANNOT_FILE}
    tabix -s1 -b2 -e2 ${ANNOT_FILE}.gz

    HDR_FILE=$(mktemp)

    echo -e '##INFO=<ID=n_hits,Number=1,Type=Integer,Description="Number of repeats found in insertion">' >> ${HDR_FILE}
    echo -e '##INFO=<ID=match_lengths,Number=.,Type=Integer,Description="Insertion lengths spanned by each repeat">' >> ${HDR_FILE}
    echo -e '##INFO=<ID=repeat_ids,Number=.,Type=String,Description="Repeat family IDs">' >> ${HDR_FILE}
    echo -e '##INFO=<ID=matching_classes,Number=.,Type=String,Description="Repeat class names">' >> ${HDR_FILE}
    echo -e '##INFO=<ID=fragmts,Number=.,Type=Integer,Description="Number of fragments merged into one by one code">' >> ${HDR_FILE}
    echo -e '##INFO=<ID=RM_hit_strands,Number=.,Type=String,Description="RepeatMasker hit strands">' >> ${HDR_FILE}
    echo -e '##INFO=<ID=RM_hit_IDs,Number=.,Type=String,Description="RepeatMasker hit IDs">' >> ${HDR_FILE}
    echo -e '##INFO=<ID=total_match_length,Number=1,Type=Integer,Description="Insertion length spanned by repeats">' >> ${HDR_FILE}
    echo -e '##INFO=<ID=total_match_span,Number=1,Type=Float,Description="Insertion span spanned by repeats">' >> ${HDR_FILE}
    echo -e '##INFO=<ID=mam_filter_1,Number=1,Type=String,Description="L1 with 5P inversion; strand">' >> ${HDR_FILE}
    echo -e '##INFO=<ID=mam_filter_2,Number=1,Type=String,Description=">90% of the SV is a SVA VNTR; SVA_family:consensus_hit_coordinates">' >> ${HDR_FILE}
    echo -e '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> ${HDR_FILE}

    cat <(grep '#' ${VCF}) <(grep -v '#' ${VCF} | sort -k1,1 -k2,2n) > genotypes.sorted.vcf
    bcftools annotate -a ${ANNOT_FILE}.gz -h ${HDR_FILE} \
    -c CHROM,POS,~ID,REF,ALT,INFO/n_hits,INFO/fragmts,INFO/match_lengths,INFO/repeat_ids,INFO/matching_classes,INFO/RM_hit_strands,INFO/RM_hit_IDs,INFO/total_match_length,INFO/total_match_span,INFO/mam_filter_1,INFO/mam_filter_2 genotypes.sorted.vcf | \
    bcftools view -Oz -o ${OUT_VCF}

else
    echo "Mammalian filters OFF, writing vcf..."
    bgzip ${ANNOT_FILE}
    tabix -s1 -b2 -e2 ${ANNOT_FILE}.gz

    HDR_FILE=$(mktemp)

    echo -e '##INFO=<ID=n_hits,Number=1,Type=Integer,Description="Number of repeats found in insertion">' >> ${HDR_FILE}
    echo -e '##INFO=<ID=match_lengths,Number=.,Type=Integer,Description="Insertion lengths spanned by each repeat">' >> ${HDR_FILE}
    echo -e '##INFO=<ID=repeat_ids,Number=.,Type=String,Description="Repeat family IDs">' >> ${HDR_FILE}
    echo -e '##INFO=<ID=matching_classes,Number=.,Type=String,Description="Repeat class names">' >> ${HDR_FILE}
    echo -e '##INFO=<ID=fragmts,Number=.,Type=Integer,Description="Number of fragments merged into one by one code">' >> ${HDR_FILE}
    echo -e '##INFO=<ID=RM_hit_strands,Number=.,Type=String,Description="RepeatMasker hit strands">' >> ${HDR_FILE}
    echo -e '##INFO=<ID=RM_hit_IDs,Number=.,Type=String,Description="RepeatMasker hit IDs">' >> ${HDR_FILE}
    echo -e '##INFO=<ID=total_match_length,Number=1,Type=Integer,Description="Insertion length spanned by repeats">' >> ${HDR_FILE}
    echo -e '##INFO=<ID=total_match_span,Number=1,Type=Float,Description="Insertion span spanned by repeats">' >> ${HDR_FILE}
    echo -e '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> ${HDR_FILE}

    cat <(grep '#' ${VCF}) <(grep -v '#' ${VCF} | sort -k1,1 -k2,2n) > genotypes.sorted.vcf
    bcftools annotate -a ${ANNOT_FILE}.gz -h ${HDR_FILE} \
    -c CHROM,POS,~ID,REF,ALT,INFO/n_hits,INFO/fragmts,INFO/match_lengths,INFO/repeat_ids,INFO/matching_classes,INFO/RM_hit_strands,INFO/RM_hit_IDs,INFO/total_match_length,INFO/total_match_span genotypes.sorted.vcf | \
    bcftools view -Oz -o ${OUT_VCF}
fi
