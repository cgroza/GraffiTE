#! /bin/bash

SVSEQ=$1 # SV_sequences_L_R_trimmed_WIN.fa
FLANK=$2 # flanking_sequences.fasta
VERBOSE=$3 # for debug only

# get each TE/SV in the variable indels
indels=$(grep '>' ${SVSEQ} | sed 's/>//g;s/__/\t/g' | cut -f 1 | sort | uniq)
#indels=""
# clean the summary file if exists
rm TSD_summary.txt 2> /dev/null

# def main function
function tsdfind {
# loop over each TE/SV
while IFS= read -r i
do
# get the strand from the RepeatMasker file
strand=$(grep -w "${i}" indels.fa.onecode.out | cut -f 9 | sort | uniq)
TE=$(grep -w "${i}" indels.fa.onecode.out | cut -f 10 | sort | uniq)
DIV=$(grep -w "${i}" indels.fa.onecode.out | cut -f 2 | sort | uniq)
echo ""
echo ""
echo ""
echo "--- TSD search for ${i} ---"
echo ""
# create 5' and 3' fragments: L = [WIN bp 5' flank][WIN bp 5' SV] R = [WIN bp 3' SV][WIN bp 3' flank]
cat <(echo ">L (5')") <(paste -d "" <(grep -A 1 "${i}__L" ${FLANK} | tail -n 1) <(grep -A 1 "${i}__L" ${SVSEQ} | tail -n 1) <(grep -A 1 "${i}__R" ${SVSEQ} | tail -n 1) <(grep -A 1 "${i}__R" ${FLANK} | tail -n 1) | sed 's/N./ /g' | awk '{print $1}') > L.fasta
cat <(echo ">R (3')") <(paste -d "" <(grep -A 1 "${i}__L" ${FLANK} | tail -n 1) <(grep -A 1 "${i}__L" ${SVSEQ} | tail -n 1) <(grep -A 1 "${i}__R" ${SVSEQ} | tail -n 1) <(grep -A 1 "${i}__R" ${FLANK} | tail -n 1) | sed 's/N./ /g' | awk '{print $2}') > R.fasta
# print 5' and 3' fragments to compare with ruler
cat <(awk 'getline seq {print $0"\n"seq"\n||||||||||||||||||||||||||||||||||||||||||||||||||\n1   5    10   15   20   25   30   35   40   45   50"}' L.fasta) <(awk 'getline seq {print $0"\n"seq"\n||||||||||||||||||||||||||||||||||||||||||||||||||\n1   5    10   15   20   25   30   35   40   45   50"}' R.fasta)
echo ""
# get 5' (L) SV sequence length
Llen=$(awk 'BEGIN {OFS = "\n"}; /^>/ {print(substr(sequence_id, 2)" "sequence_length); sequence_length = 0; sequence_id = $0}; /^[^>]/ {sequence_length += length($0)}; END {print(sequence_length)}' L.fasta | tail -n 1)
# get buffer size for prefix or trailing A or T
# 5' (L) --> [XXXXXXXXXXTTTTTTTTTT] gather trailing T at the end of the left (5') fragment (reverse TE insertion)
poly_L=$(grep "[T]*TTT\b" <(tail -n 1 L.fasta) -o)
# 3' (R) --> [AAAAAAAAAAXXXXXXXXXX] gather prefixes A at the beginning of the right (3') fragment (forward TE insertion) 
poly_R=$(grep "^AAA[A]*" <(tail -n 1 R.fasta) -o)
echo ""
if (( ${#poly_L}-3 < 1 )); then echo "5' poly_T: ${#poly_L} bp, will not remove anything for alignment"; offsetL=0; else echo "5' poly_T: ${#poly_L} bp, will remove $((${#poly_L}-3)) trailing T for alignment"; offsetL=$((${#poly_L}-3));fi
if (( ${#poly_R}-3 < 1 )); then echo "3' poly_A: ${#poly_R} bp, will not remove anything for alignment"; offsetR=0; else echo "3' poly_A: ${#poly_R} bp, will remove $((${#poly_R}-3)) starting A for alignment"; offsetR=$((${#poly_R}-3));fi
echo ""
# trim poly-A or -T for alignments
awk '/>/{getline seq; sub(/T+TT$/,"TTT", seq); print $0"\n"seq}' L.fasta > L.short.fasta
awk '/>/{getline seq; sub(/^AA+A/,"AAA", seq); print $0"\n"seq}' R.fasta > R.short.fasta
# get short sequences length
Lshort=$(awk 'BEGIN {OFS = "\n"}; /^>/ {print(substr(sequence_id, 2)" "sequence_length); sequence_length = 0; sequence_id = $0}; /^[^>]/ {sequence_length += length($0)}; END {print(sequence_length)}' L.short.fasta | tail -n 1)
Rshort=$(awk 'BEGIN {OFS = "\n"}; /^>/ {print(substr(sequence_id, 2)" "sequence_length); sequence_length = 0; sequence_id = $0}; /^[^>]/ {sequence_length += length($0)}; END {print(sequence_length)}' R.short.fasta | tail -n 1)

# check if no short sequence is empty, in which case call bypass alignment and return no TSD
if (( $Lshort == 0 ))
then
	echo -e "$i" "$TE" "$strand" "$DIV" "NA\tNA\tNA\tNA\tNA\tNA\tNA\tpoly-T\tNA\tFAIL" | sed 's/\n//g;s/ /\t/g' | tee -a TSD_summary.txt
elif (( $Rshort == 0 ))
	then
		echo -e "$i" "$TE" "$strand" "$DIV" "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tpoly-A\tFAIL" | sed 's/\n//g;s/ /\t/g' | tee -a TSD_summary.txt
else
	# water to generate match
	# print the match
	water -auto -stdout -gapextend 10 -gapopen 50 -asequence L.short.fasta -bsequence R.short.fasta | grep -v "#\|^$" | awk -v offsetR=${offsetR} '/^R/ {$2=(($2+offsetR)); if ($2 < 10) {print $1"                 "$2"  "$3"     "$4} else {print $1"                 "$2" "$3"     "$4};next} !/^R/ {print $0}'
	echo ""
	# store match info to variables
	eval $(water -auto -stdout -gapextend 10 -gapopen 50 -asequence L.short.fasta -bsequence R.short.fasta | awk '/Length/ {print "length="$3; next} /Identity/ {gsub("/.*","", $3); print "matches="$3; next} /Gaps/ {gsub("/.*","", $3); print "gaps="$3; next} /^L/ {print "Lstart="$2"\nLend="$4; next} /^R/ {print "Rstart="$2"\nRend="$4; next}')
	# print 5' and 3' with underlined hits
	echo "candidate TSDs:"
	paste <(echo -e $(awk -v Lstart=${Lstart} -v Lend=${Lend} 'getline seq {printf substr(seq, 1,((Lstart-1))) "\\e[4m"substr(seq, ((Lstart)),((Lend-Lstart+1)))"\\e[0m" substr(seq, ((Lend+1)))}' L.fasta)) <(echo -e $(awk -v Rstart=$((${Rstart}+${offsetR})) -v Rend=$((${Rend}+${offsetR})) 'getline seq {printf substr(seq, 1,((Rstart-1))) "\\e[4m"substr(seq, ((Rstart)),((Rend-Rstart+1)))"\\e[0m" substr(seq, ((Rend+1)))}' R.fasta))
	echo ""
	# print match info, add >>> <<< if selected
	#echo "Alignment=$length," "Mismatches=$(($length-$matches))," "Gaps=$gaps" | awk -v LL="$Llen" -v MM=$(($length-$matches)) -v GP="$gaps" -v RS="$Rstart" -v LE="$Lend" '{if (RS <= 5 && LE > LL-5 && MM+GP < 2) {print $0"\tFILTER PASS!!!"} else {print $0"\tFILTER FAIL"}}' 
	#echo "R" "L" "$length" "$(($length-$matches))" "$gaps" "$Rstart" "$Rend" "$Lstart" "$Lend" | awk -v size="$Llen" -v pL="${#poly_L}" -v pR="${#poly_R}" '{if ($6 <= 5+pR && $9 > size-pL && ($4+$5) < 2 ) {print ">>>"$0"<<<\tTSD FOUND!!!"} else {print $0"\tNO TSD"}}' 
	# [[ADD FILE!!!]] write output to simple file
	TSDs=$(awk -v Lstart=${Lstart} -v Lend=${Lend} 'getline seq {print substr(seq,Lstart,Lend-Lstart+1)}' L.fasta)" "$(awk -v Rstart=$((${Rstart}+${offsetR})) -v Rend=$((${Rend}+${offsetR})) 'getline seq {print substr(seq,Rstart,Rend-Rstart+1)}' R.fasta)
	# split TSDs in L and R TSDs:
	LTSD=$(echo "$TSDs" | awk '{print $1}')
	RTSD=$(echo "$TSDs" | awk '{print $2}')

	# create test variable to see if the candidate TSDs are PASS: must be in the +/- 5bp of the flanking (+/- offset for poly-A/T) + having either <= 1 GAP+MM or (GAP+MM)/L <= DIV of the TE
	test=$(echo -e "$i" "$TE" "$strand" "$DIV" "$length" "$(($length-$matches))" "$gaps" "$(($Lend-$Llen-1))" "$offsetL" "$((Rstart+$offsetR))" "$offsetR" "$TSDs" | sed 's/\n//g;s/ /\t/g' | awk '{if ($8+$9 <= 0 && $8+$9 > -6 && $10-$11 <= 5 && (($6+$7)/$5 <= $4/100 || ($6+$7 <= 1))) {print $0"\tPASS"} else {print $0"\tFAIL"}}' | awk '{print $NF}')
	output=$(echo -e "$i" "$TE" "$strand" "$DIV" "$length" "$(($length-$matches))" "$gaps" "$(($Lend-$Llen-1))" "$offsetL" "$((Rstart+$offsetR))" "$offsetR" "$TSDs" | sed 's/\n//g;s/ /\t/g' | awk '{if ($8+$9 <= 0 && $8+$9 > -6 && $10-$11 <= 5 && (($6+$7)/$5 <= $4/100 || ($6+$7 <= 1))) {print $0"\tPASS"} else {print $0"\tFAIL"}}')
	# debug
	echo ""
	echo "$test"
	echo ""
	# poly-A/T elongation
	if [[ "$test" == "PASS" ]]
		then
			#echo -e "SVname\tTEname\tStrand\tDiv\tAlnLen\tMM\tGaps\t5P_TSD_end\t5P_offset\t3P_TSD_start\t3P_offset\t5P_TSD\t3P_TSD"
			#echo "$output"
			if (( $offsetL > 0 && $Lend == $Lshort ))
				then
				echo "3' end: can extend poly T"
				# count how many T are left in the L (5') sequence (in the poly-T)
				# take the L short sequence
				LSseq=$(tail -n 1 L.short.fasta)
				# grep it out of the long sequence and count the T
				Lleft=$(sed "s/$LSseq//g" <(tail -n 1 L.fasta))
				LLcount=$(echo ${#Lleft})
				# count how many T are left in the R (3') sequence on the 3' of the last matching residue
				Rleft=$(sed "s/$RTSD/\t/g" <(tail -n 1 R.fasta) | cut -f 2 | grep "^[T]*T" -o)
				RLcount=$(echo ${#Rleft})
				# now find the smallest and at those Ts to the TSDs
				TtoAdd=$(echo $((LLcount<RLcount ? LLcount : RLcount)))
				echo "we could add $TtoAdd T to the 3' of the TSDs"
				# adds the T to the TSDs
				LTSDext=$(awk -v Ttime=${TtoAdd} -v LTSD=${LTSD} 'BEGIN{for(c=0;c<Ttime;c++) printf LTSD; printf "T"; printf "\n"}')
				RTSDext=$(awk -v Ttime=${TtoAdd} -v RTSD=${RTSD} 'BEGIN{for(c=0;c<Ttime;c++) printf RTSD; printf "T"; printf "\n"}')
			else
				echo "3' end: nothing to extend"
				LTSDext=${LTSD}
				RTSDext=${RTSD}
			fi
			if (( $offsetR > 0 && $Rstart == 1 ))
				then
				echo "5' end: can extend poly A"
				# count how many A are left in the R (3') sequence (in the poly-A)
				# take the R short sequence
				RSseq=$(tail -n 1 R.short.fasta)
				# grep it out of the long sequence and count the A
				Rleft2=$(sed "s/$RSseq//g" <(tail -n 1 R.fasta))
				LRcount=$(echo ${#Rleft2})
				# count how many A are left in the L (5') sequence on the 5' of the last matching residue
				Lleft2=$(sed "s/$LTSD/\t/g" <(tail -n 1 L.fasta) | cut -f 1 | grep "[A]*A\b" -o)
				RLcount2=$(echo ${#Lleft2})
				AtoAdd=$(echo $((LRcount<RLcount2 ? LRcount : RLcount2)))
				echo "we could add $AtoAdd A to the 5' of the TSDs"
				# adds the A to the TSDs
				LTSDextF=$(awk -v Ttime=${AtoAdd} -v LTSD=${LTSDext} 'BEGIN{for(c=0;c<Ttime;c++) printf "A"; printf LTSD; printf "\n"}')
				RTSDextF=$(awk -v Ttime=${AtoAdd} -v RTSD=${RTSDext} 'BEGIN{for(c=0;c<Ttime;c++) printf "A"; printf RTSD; printf "\n"}')
			else
				echo "5' end: nothing to extend"
				LTSDextF=${LTSDext}
				RTSDextF=${RTSDext}
			fi
			echo -e "SVname\tTEname\tStrand\tDiv\tAlnLen\tMM\tGaps\t5P_TSD_end\t5P_offset\t3P_TSD_start\t3P_offset\t5P_TSD\t3P_TSD"
			echo -e "$i" "$TE" "$strand" "$DIV" "$length" "$(($length-$matches))" "$gaps" "$(($Lend-$Llen-1))" "$offsetL" "$((Rstart+$offsetR))" "$offsetR" "$LTSDextF" "$RTSDextF" "PASS" | sed 's/\n//g;s/ /\t/g' | tee -a TSD_summary.txt
		else
			echo -e "SVname\tTEname\tStrand\tDiv\tAlnLen\tMM\tGaps\t5P_TSD_end\t5P_offset\t3P_TSD_start\t3P_offset\t5P_TSD\t3P_TSD"
			echo "$output" | tee -a TSD_summary.txt
	fi # close loop for elongation
fi # close loop that check if one end is poly-A or poly-T and skip
# close loop and feed it with each selected (1 TE hit) SV name
done <<< "$indels" 
}

# exec and print output according to verbose option
if [[ ${VERBOSE} == "V" ]]
then
	tsdfind | tee TSD_full_log.txt
else
	tsdfind > TSD_full_log.txt
fi