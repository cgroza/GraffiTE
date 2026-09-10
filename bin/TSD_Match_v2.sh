#! /bin/bash
# This is V2 of this script, with main diff being BLASTN to replace EMBOSS:WATER
# Author: Clément Goubert, 2022
# email: goubert.clement@gmail.com

SVSEQ=$1 # SV_sequences_L_R_trimmed_WIN.fa
FLANK=$2 # flanking_sequences.fasta
indels=$(cat $3)
WIN=${4:-30} # flank width prepTSD.sh used (--tsd_win); the junction is column WIN of each fragment
VERBOSE=$5   # for debug only
name=$(cat $3 | head -n 1)
[[ "${WIN}" =~ ^[1-9][0-9]*$ ]] || { echo "TSD_Match_v2.sh: window must be a positive integer, got '${WIN}'" >&2; exit 1; }
# rulers for the log, one column per base of a [WIN][WIN] fragment
RULER_BAR=$(printf '%*s' $((2 * WIN)) '' | tr ' ' '|')
RULER_NUM=$(awk -v n=$((2 * WIN)) 'BEGIN { s = "1"; for (i = 5; i <= n; i += 5) { pad = i - length(s) - 1; s = s sprintf("%*s%d", pad, "", i) } print s }')

# A missing matcher would otherwise read as "no_hit" on every variant.
command -v exact_match.py > /dev/null || { echo "TSD_Match_v2.sh: exact_match.py is not on PATH" >&2; exit 1; }
for f in "${SVSEQ}" "${FLANK}"; do
	[[ -s "${f}" ]] || { echo "TSD_Match_v2.sh: ${f} is empty" >&2; exit 1; }
done

# clean the summary file if exists
rm $name.TSD_summary.txt 2> /dev/null

# def main function
function tsdfind {
# loop over each TE/SV
while IFS= read -r i
do

echo ""
echo ""
echo ""
echo "--- TSD search for ${i} ---"
echo ""
# create 5' and 3' fragments: L = [WIN bp 5' flank][WIN bp 5' SV] R = [WIN bp 3' SV][WIN bp 3' flank] 
cat <(echo ">L|5P_end") <(paste -d '\0' <(grep -A 1 "${i}__L" ${FLANK} | tail -n 1) <(grep -A 1 "${i}__L" ${SVSEQ} | tail -n 1) <(echo " ") <(grep -A 1 "${i}__R" ${SVSEQ} | tail -n 1) <(grep -A 1 "${i}__R" ${FLANK} | tail -n 1) | awk '{print $1}') > L.fasta
cat <(echo ">R|3P_end") <(paste -d '\0' <(grep -A 1 "${i}__L" ${FLANK} | tail -n 1) <(grep -A 1 "${i}__L" ${SVSEQ} | tail -n 1) <(echo " ") <(grep -A 1 "${i}__R" ${SVSEQ} | tail -n 1) <(grep -A 1 "${i}__R" ${FLANK} | tail -n 1) | awk '{print $2}') > R.fasta

# print 5' and 3' fragments to compare with ruler
cat <(awk -v bar="${RULER_BAR}" -v num="${RULER_NUM}" 'getline seq {print $0"\n"seq"\n"bar"\n"num}' L.fasta) <(awk -v bar="${RULER_BAR}" -v num="${RULER_NUM}" 'getline seq {print $0"\n"seq"\n"bar"\n"num}' R.fasta)
echo ""

# now we "blast" no matter what and I will just save the table for now and apply TSD length filters to avoid longer matches in repetitive regions
# TSD_MIN: default = 4 ; min = 4 max = 30
# TSD_MAX: default = 20; min = 4 max = 30
TSD_MIN=4
TSD_MAX=20
exact_match.py -word_size 4 -query R.fasta -db L.fasta -outfmt 6 -strand plus | awk -v tsdmin=${TSD_MIN} -v tsdmax=${TSD_MAX} -v win=${WIN} 'function abs(x) { return x < 0 ? -x : x } function min(x,y) { return x < y ? x : y } { a = (abs(win-$7) + abs(win-$9)) / 2; b = (abs(win-$8) + abs(win-$10)) / 2; score = min(a, b);if($4 >= tsdmin && $4 <= tsdmax){print $0"\t"(win-$7)"\t"(win-$9)"\t"(win-$8)"\t"(win-$10)"\t"score} }' > blastout 2>&1


if ! [[ -s blastout ]]

then # no match detected

	output=$(echo -e "${i}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tno_hit\tno_hit\tFAIL" | sed 's/\n//g;s/ /\t/g')
	# test="FAIL"
	# # debug
	# echo ""
	# echo "$test"
	# echo ""

else # match found

	# get best hit using lowest (best) TSD score (= how close to the edges of the SV/flank the TSD are -- it should be ideally: [TSD----]TSD or TSD[-----TSD])
	# example of good scoring: TSD are snug with the breakpoints, one on the genome, one on the SV:
	# best hit: R|3P_end	L|5P_end	100.000	23	0	0	31	53	31	53	1.00e-20	46.		-1	-1	-23	-23	1 <---- score is 1 (should be 0, need to fix 31 not 30 as reference point)
	# 
	# >L|5P_end                     ***********************       
	# acataaaatatcaaagtacccaaactatacATTATATACTGTACATAAAATATAAAATTA
	# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
	# 1   5    10   15   20   25   30   35   40   45   50   55   60
	#
	# >R|3P_end                     ***********************
	# TGTACATAAAATAAAGTACACAAACTATAAattatatactgtacataaaatatgaaatta
	# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
	# 1   5    10   15   20   25   30   35   40   45   50   55   60
	#
	#
	# candidate TSDs:
	# acataaaatatcaaagtacccaaactatac [ATTATATACTGTACATAAAATATAAAATTA---SV/TE(s)---TGTACATAAAATAAAGTACACAAACTATAA] attatatactgtacataaaatatgaaatta
	#                                 ^^^^^^^^^^^^^^^^^^^^^^^                                                     ^^^^^^^^^^^^^^^^^^^^^^^
	
	# get best hit
	sort -k17,17n -k4,4nr blastout | head -n 1 > best_hit
	# print the candidate hits and best hit
	echo -e "R_query\tL_target\tidty\tmatch_len\tMM\tgaps\tR_start\tR_end\tL_start\tL_end\te-value\tR_start_offset\tL_start_offset\tR_end_offset\tL_end_offset\tTSD_score" > header
	echo ""
	echo "candidate hits:"
	cat header blastout
	echo ""
	echo "best hit:"
	cat header best_hit
	echo ""

	# display TSD model

	Lstart=$(awk '{print $9}' best_hit)
	Rstart=$(awk '{print $7}' best_hit)
	Lend=$(awk '{print $10}' best_hit)
	Rend=$(awk '{print $8}' best_hit)

	echo "candidate TSDs:"
	paste -d '\0' <(echo -e $(awk -v Lstart=${Lstart} -v Lend=${Lend} 'getline seq {printf substr(seq, 1,((Lstart-1))) "\\e[4m"substr(seq, ((Lstart)),((Lend-Lstart+1)))"\\e[0m" substr(seq, ((Lend+1)))}' L.fasta)) <(echo -e "---SV/TE(s)---") <(echo -e $(awk -v Rstart=$((${Rstart})) -v Rend=$((${Rend})) 'getline seq {printf substr(seq, 1,((Rstart-1))) "\\e[4m"substr(seq, ((Rstart)),((Rend-Rstart+1)))"\\e[0m" substr(seq, ((Rend+1)))}' R.fasta))

	# grab the sequences
	L_TSD=$(awk -v Lstart=${Lstart} -v Lend=${Lend} 'getline seq {printf substr(seq, Lstart, Lend-Lstart+1)}' L.fasta)
	R_TSD=$(awk -v Rstart=${Rstart} -v Rend=${Rend} 'getline seq {printf substr(seq, Rstart, Rend-Rstart+1)}' R.fasta)

	# add sequences to summary report and assign PASS/FAIL - We have an opportunity here to offer user parameters
	output=$(awk -v sv=${i} -v ltsd=${L_TSD} -v rtsd=${R_TSD} '{ if($NF <= 5) { print sv"\t"$0"\t"ltsd"\t"rtsd"\tPASS" } else { print sv"\t"$0"\t"ltsd"\t"rtsd"\tFAIL" } }' best_hit)
fi

echo ""
echo "$output" | tee -a $name.TSD_summary.txt

done <<< "$indels"
}

# exec and print output according to verbose option
if [[ ${VERBOSE} == "V" ]]
then
	tsdfind | tee $name.TSD_full_log.txt
else
	tsdfind > $name.TSD_full_log.txt
fi