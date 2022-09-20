#!/bin/bash

VCF=$1
OUT_VCF=$2
FASTA_LIB=$3

FASTA_FILE=indels.fa
cat ${VCF} | grep -v "#" | grep "SVTYPE=INS" | awk '{print(sprintf(">%s\n%s", $3, $5))}' > ${FASTA_FILE}
cat ${VCF} | grep -v "#" | grep "SVTYPE=DEL" | awk '{print(sprintf(">%s\n%s", $3, $4))}' >> ${FASTA_FILE}

mkdir repeatmasker_dir
REPMASK_DIR=repeatmasker_dir

# detect number of cores allocated to this nextflow process
repmask_cores="$(($(nproc)/4))"

# guard against number of cores smaller than 1
if [ "${repmask_cores}" -lt "1" ]; then
    repmask_cores="1"
fi


RepeatMasker -nolow -lib ${FASTA_LIB} -s -dir ${REPMASK_DIR} -pa ${repmask_cores} ${FASTA_FILE}

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

annotate_vcf.R --dotout ${REPMASK_ONECODE_OUT} --vcf ${VCF} --annotation ${ANNOT_FILE}

bgzip ${ANNOT_FILE}
tabix -s1 -b2 -e2 ${ANNOT_FILE}.gz

HDR_FILE=$(mktemp)

echo -e '##INFO=<ID=n_hits,Number=1,Type=Integer,Description="Number of repeats found in insertion">' >> ${HDR_FILE}
echo -e '##INFO=<ID=match_lengths,Number=.,Type=Integer,Description="Insertion lengths spanned by each repeat">' >> ${HDR_FILE}
echo -e '##INFO=<ID=repeat_ids,Number=.,Type=String,Description="Repeat family IDs">' >> ${HDR_FILE}
echo -e '##INFO=<ID=matching_classes,Number=.,Type=String,Description="Repeat class names">' >> ${HDR_FILE}
echo -e '##INFO=<ID=fragmts,Number=.,Type=Integer,Description="Number of fragments merged into one by one code">' >> ${HDR_FILE}
echo -e '##INFO=<ID=total_match_length,Number=1,Type=Integer,Description="Insertion length spanned by repeats">' >> ${HDR_FILE}
echo -e '##INFO=<ID=total_match_span,Number=1,Type=Float,Description="Insertion span spanned by repeats">' >> ${HDR_FILE}
echo -e '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> ${HDR_FILE}

bcftools annotate -a ${ANNOT_FILE}.gz -h ${HDR_FILE} \
         -c CHROM,POS,INFO/n_hits,INFO/fragmts,INFO/match_lengths,INFO/repeat_ids,INFO/matching_classes,INFO/total_match_length,INFO/total_match_span $VCF | \
    bcftools view -Oz -o ${OUT_VCF}

