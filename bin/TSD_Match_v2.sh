#! /bin/bash
# This is V2 of this script, with main diff being BLASTN to replace EMBOSS:WATER
# Author: Clément Goubert, 2022
# email: goubert.clement@gmail.com

SVSEQ=$1 # SV_sequences_L_R_trimmed_WIN.fa
FLANK=$2 # flanking_sequences.fasta
indels=$(cat $3)  
VERBOSE=$4 # for debug only
name=$(cat $3 | head -n 1)

# clean the summary file if exists
rm $name.TSD_summary.txt 2> /dev/null

# def main function
function tsdfind {
# loop over each TE/SV
while IFS= read -r i
do
# check if the variant is a L1 with 5P inversion (OBSOLETE)
# L15P=$(grep -w ${i} genotypes_repmasked_filtered.vcf | sed 's/mam_filter_1=/\t/g;s/;mam_filter_2=/\t/g;s/5P_INV:plus/+/g;s/5P_INV:minus/C/g' | awk '{print $9}')
# if [[ $L15P == "None" || $L15P == "GT" ]] # GT appears if params.mammal is false in the Nextflow pipeline.
# then
# 	# get all info from the RepeatMasker file
# 	strand=$(grep -w "${i}" indels.fa.onecode.out | cut -f 9 | sort | uniq)
# 	TE=$(grep -w "${i}" indels.fa.onecode.out | cut -f 10 | sort | uniq)
# 	DIV=$(grep -w "${i}" indels.fa.onecode.out | cut -f 2 | sort | uniq)
# else
# 	# get the strand from the vcf
# 	strand=$L15P
# 	# get the TE name from the RepeatMasker file
# 	TE=$(grep -w "${i}" indels.fa.onecode.out | cut -f 10 | sort | uniq)
# 	# get the divergence by averaging the two hits in the RepeatMasker file (weighted average by hit length)
# 	DIV=$(grep -w "${i}" indels.fa.onecode.out | awk '{print $2"\t"($7-$6)}' | awk 'getline second {print $0"\t"second}' | awk '{print ($1*$2+$3*$4)/($2+$4)}')
# fi

echo ""
echo ""
echo ""
echo "--- TSD search for ${i} ---"
echo ""
# create 5' and 3' fragments: L = [WIN bp 5' flank][WIN bp 5' SV] R = [WIN bp 3' SV][WIN bp 3' flank] 
# cat <(echo ">L|5P_end") <(paste -d "" <(grep -A 1 "${i}__L" ${FLANK} | tail -n 1) <(grep -A 1 "${i}__L" ${SVSEQ} | tail -n 1 | sed 's/N*$/ /g') <(grep -A 1 "${i}__R" ${SVSEQ} | tail -n 1 | sed 's/^N*N/ /g') <(grep -A 1 "${i}__R" ${FLANK} | tail -n 1) | awk '{print $1}') > L.fasta
# cat <(echo ">R|3P_end") <(paste -d "" <(grep -A 1 "${i}__L" ${FLANK} | tail -n 1) <(grep -A 1 "${i}__L" ${SVSEQ} | tail -n 1 | sed 's/N*$/ /g') <(grep -A 1 "${i}__R" ${SVSEQ} | tail -n 1 | sed 's/^N*N/ /g') <(grep -A 1 "${i}__R" ${FLANK} | tail -n 1) | awk '{print $2}') > R.fasta
cat <(echo ">L|5P_end") <(paste -d "" <(grep -A 1 "${i}__L" ${FLANK} | tail -n 1) <(grep -A 1 "${i}__L" ${SVSEQ} | tail -n 1) <(echo " ") <(grep -A 1 "${i}__R" ${SVSEQ} | tail -n 1) <(grep -A 1 "${i}__R" ${FLANK} | tail -n 1) | awk '{print $1}') > L.fasta
cat <(echo ">R|3P_end") <(paste -d "" <(grep -A 1 "${i}__L" ${FLANK} | tail -n 1) <(grep -A 1 "${i}__L" ${SVSEQ} | tail -n 1) <(echo " ") <(grep -A 1 "${i}__R" ${SVSEQ} | tail -n 1) <(grep -A 1 "${i}__R" ${FLANK} | tail -n 1) | awk '{print $2}') > R.fasta

# print 5' and 3' fragments to compare with ruler
cat <(awk 'getline seq {print $0"\n"seq"\n||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n1   5    10   15   20   25   30   35   40   45   50   55   60"}' L.fasta) <(awk 'getline seq {print $0"\n"seq"\n||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n1   5    10   15   20   25   30   35   40   45   50   55   60"}' R.fasta)
echo ""

# get 5' (L) SV sequence length (old routine)
#Llen=$(awk 'BEGIN {OFS = "\n"}; /^>/ {print(substr(sequence_id, 2)" "sequence_length); sequence_length = 0; sequence_id = $0}; /^[^>]/ {sequence_length += length($0)}; END {print(sequence_length)}' L.fasta | tail -n 1)

# FOR NOW, WE REMOVE THIS AS WE IMPLEMENT THE NAIVE APPROCH FOR ALL SV, THUS WE DON'T HAVE A PRIOR OF TE ORIENTATION
# get offset size for prefix or trailing A or T
# if [[ $strand == "+" ]]
# then
# 	# do nothing on 5'
# 	echo "5' poly_T: element is in ${strand} orientation, will not search for poly_T"
#  	offsetL=0
#  	cp L.fasta L.short.fasta
#  	#3' (R) --> [AAAAAAAAAAXXXXXXXXXX] gather prefixes A at the beginning of the right (3') fragment (forward TE insertion) 
#  	poly_R=$(grep "^AAA[A]*" <(tail -n 1 R.fasta) -o)
# 	if (( ${#poly_R}-3 < 1 ))
# 	then
# 		echo "3' poly_A: ${#poly_R} bp, will not remove anything for alignment"
# 		offsetR=0
# 		cp R.fasta R.short.fasta
# 	else
# 		echo "3' poly_A: ${#poly_R} bp, will remove $((${#poly_R}-3)) starting A for alignment"
# 		offsetR=$((${#poly_R}-3))
# 		awk '/>/{getline seq; sub(/^AA+A/,"AAA", seq); print $0"\n"seq}' R.fasta > R.short.fasta
# 	fi
# else
# 	# do nothing on 3'
# 	echo "3' poly_A: element is in ${strand} orientation, will not search for poly_A"
# 	offsetR=0
# 	cp R.fasta R.short.fasta
# 	# 5' (L) --> [XXXXXXXXXXTTTTTTTTTT] gather trailing T at the end of the left (5') fragment (reverse TE insertion)
# 	poly_L=$(grep "[T]*TTT$" <(tail -n 1 L.fasta) -o)
# 	if (( ${#poly_L}-3 < 1 ))
# 	then
# 		echo "5' poly_T: ${#poly_L} bp, will not remove anything for alignment"
# 		offsetL=0
# 		cp L.fasta L.short.fasta
# 	else
# 		echo "5' poly_T: ${#poly_L} bp, will remove $((${#poly_L}-3)) starting A for alignment"
# 		offsetL=$((${#poly_L}-3))
# 		awk '/>/{getline seq; sub(/T+TT$/,"TTT", seq); print $0"\n"seq}' L.fasta > L.short.fasta
# 	fi
# fi

# # get short sequences length (old routine)
# Lshort=$(awk 'BEGIN {OFS = "\n"}; /^>/ {print(substr(sequence_id, 2)" "sequence_length); sequence_length = 0; sequence_id = $0}; /^[^>]/ {sequence_length += length($0)}; END {print(sequence_length)}' L.short.fasta | tail -n 1)
# Rshort=$(awk 'BEGIN {OFS = "\n"}; /^>/ {print(substr(sequence_id, 2)" "sequence_length); sequence_length = 0; sequence_id = $0}; /^[^>]/ {sequence_length += length($0)}; END {print(sequence_length)}' R.short.fasta | tail -n 1)

# # check if no short sequence is empty, in which case call bypass alignment and return no TSD
# if (( $Lshort == 0 ))
# then
# 	echo -e "$i" "$TE" "$strand" "$DIV" "NA\tNA\tNA\tNA\tNA\tNA\tNA\tpoly-T\tNA\tFAIL" | sed 's/\n//g;s/ /\t/g' | tee -a $name.TSD_summary.txt
# elif (( $Rshort == 0 ))
# 	then
# 		echo -e "$i" "$TE" "$strand" "$DIV" "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tpoly-A\tFAIL" | sed 's/\n//g;s/ /\t/g' | tee -a $name.TSD_summary.txt
# else
# 	#make blast db for short L sequence
# 	#makeblastdb -in L.short.fasta -out L.short.fasta -dbtype="nucl" &> /dev/null

# now we directly make a db out of the left fragment no matter what
makeblastdb -in L.fasta -out L.fasta -dbtype="nucl"

# 	makeblastdb -in L.short.fasta -out L.short.fasta -dbtype="nucl"
# 	#blast, keep outputs with less than 1 or DIV MM+GAP, sort the edge-most hit on top
# 	blastn -word_size 4 -query R.short.fasta -db L.short.fasta -outfmt 6 -strand plus | awk -v div={$DIV} '$5+$6 <= 1 || ($5+$6)/$4 <= div/100' | sort -k7,7n -k10,10nr -k4,4nr > blastout 2>&1

# now we blast no matter what and I will just save the table for now
exact_match.py -word_size 4 -query R.fasta -db L.fasta -outfmt 6 -strand plus | awk 'function abs(x) { return x < 0 ? -x : x } function min(x,y) { return x < y ? x : y } { a = (abs(30-$7) + abs(30-$9)) / 2; b = (abs(30-$8) + abs(30-$10)) / 2; score = min(a, b); print $0"\t"(30-$7)"\t"(30-$9)"\t"(30-$8)"\t"(30-$10)"\t"score }' > blastout 2>&1
sort -k17,17n -k4,4nr blastout | head -n 1 > best_hit

# 	if ! [[ -s blastout ]]
# 	then
# 		output=$(echo -e "$i" "$TE" "$strand" "$DIV" "NA\tNA\tNA\tNA\tNA\tNA\tNA\tno_hit\tno_hit\tFAIL" | sed 's/\n//g;s/ /\t/g')
# 		test="FAIL"
# 		# debug
# 		echo ""
# 		echo "$test"
# 		echo ""
# 	else

#print the candidate hits
header=$(echo -e "R_query\tL_target\tidty\tmatch_len\tMM\tgaps\tR_start\tR_end\tL_start\tL_end\tbitScore\tR_start_offset\tL_start_offset\tR_end_offset\tL_end_offset\tTSD_score")
echo ""
echo "candidate hits:"
cat <(echo ${header}) blastout
echo ""
echo "best hit:"
cat <(${header}) best_hit
echo ""

Lstart=$(awk '{print $9}' best_hit)
Rstart=$(awk '{print $7}' best_hit)
Lend=$(awk '{print $10}' best_hit)
Rend=$(awk '{print $8}' best_hit)


echo "candidate TSDs:"
paste -d "" <(echo -e $(awk -v Lstart=${Lstart} -v Lend=${Lend} 'getline seq {printf substr(seq, 1,((Lstart-1))) "\\e[4m"substr(seq, ((Lstart)),((Lend-Lstart+1)))"\\e[0m" substr(seq, ((Lend+1)))}' L.fasta)) <(echo -e "---SV/TE(s)---") <(echo -e $(awk -v Rstart=$((${Rstart})) -v Rend=$((${Rend})) 'getline seq {printf substr(seq, 1,((Rstart-1))) "\\e[4m"substr(seq, ((Rstart)),((Rend-Rstart+1)))"\\e[0m" substr(seq, ((Rend+1)))}' R.fasta))
#paste -d "" <(echo -e $(awk -v Lstart=${Lstart} -v Lend=${Lend} 'getline seq {printf substr(seq, 1,((Lstart-1))) "\\e[4m"substr(seq, ((Lstart)),((Lend-Lstart+1)))"\\e[0m" substr(seq, ((Lend+1)))}' L.fasta)) <(if [[ ${strand} == "C" ]]; then echo -e "[ <<< ${TE} ${strand} <<< ]"; else echo "[ >>> ${TE} ${strand} >>> ]";fi) <(echo -e $(awk -v Rstart=$((${Rstart}+${offsetR})) -v Rend=$((${Rend}+${offsetR})) 'getline seq {printf substr(seq, 1,((Rstart-1))) "\\e[4m"substr(seq, ((Rstart)),((Rend-Rstart+1)))"\\e[0m" substr(seq, ((Rend+1)))}' R.fasta))

# 		#take best hit and export variables to match previous water format
# 		eval $(head -n 1 blastout | awk '{print "length="$4; print "MM="$5; print "gaps="$6; print "Lstart="$9; print "Lend="$10; print "Rstart="$7; print "Rend="$8 }')
# 		# print 5' and 3' with underlined hits
# 		echo "candidate TSDs:"
# 		paste -d "" <(echo -e $(awk -v Lstart=${Lstart} -v Lend=${Lend} 'getline seq {printf substr(seq, 1,((Lstart-1))) "\\e[4m"substr(seq, ((Lstart)),((Lend-Lstart+1)))"\\e[0m" substr(seq, ((Lend+1)))}' L.fasta)) <(if [[ ${strand} == "C" ]]; then echo -e "[ <<< ${TE} ${strand} <<< ]"; else echo "[ >>> ${TE} ${strand} >>> ]";fi) <(echo -e $(awk -v Rstart=$((${Rstart}+${offsetR})) -v Rend=$((${Rend}+${offsetR})) 'getline seq {printf substr(seq, 1,((Rstart-1))) "\\e[4m"substr(seq, ((Rstart)),((Rend-Rstart+1)))"\\e[0m" substr(seq, ((Rend+1)))}' R.fasta))
# 		echo ""
# 		# print match info, add >>> <<< if selected
# 		#echo "Alignment=$length," "Mismatches=$(($length-$matches))," "Gaps=$gaps" | awk -v LL="$Llen" -v MM=$(($length-$matches)) -v GP="$gaps" -v RS="$Rstart" -v LE="$Lend" '{if (RS <= 5 && LE > LL-5 && MM+GP < 2) {print $0"\tFILTER PASS!!!"} else {print $0"\tFILTER FAIL"}}' 
# 		#echo "R" "L" "$length" "$(($length-$matches))" "$gaps" "$Rstart" "$Rend" "$Lstart" "$Lend" | awk -v size="$Llen" -v pL="${#poly_L}" -v pR="${#poly_R}" '{if ($6 <= 5+pR && $9 > size-pL && ($4+$5) < 2 ) {print ">>>"$0"<<<\tTSD FOUND!!!"} else {print $0"\tNO TSD"}}' 
# 		# [[ADD FILE!!!]] write output to simple file
# 		TSDs=$(awk -v Lstart=${Lstart} -v Lend=${Lend} 'getline seq {print substr(seq,Lstart,Lend-Lstart+1)}' L.fasta)" "$(awk -v Rstart=$((${Rstart}+${offsetR})) -v Rend=$((${Rend}+${offsetR})) 'getline seq {print substr(seq,Rstart,Rend-Rstart+1)}' R.fasta)
# 		# split TSDs in L and R TSDs:
# 		LTSD=$(echo "$TSDs" | awk '{print $1}')
# 		RTSD=$(echo "$TSDs" | awk '{print $2}')

# 		# create test variable to see if the candidate TSDs are PASS: must be in the +/- 5bp of the flanking (+/- offset for poly-A/T) + having either <= 1 GAP+MM or (GAP+MM)/L <= DIV of the TE
# 		test=$(echo -e "$i" "$TE" "$strand" "$DIV" "$length" "$MM" "$gaps" "$(($Lend-$Llen-1))" "$offsetL" "$((Rstart+$offsetR))" "$offsetR" "$TSDs" | sed 's/\n//g;s/ /\t/g' | awk '{if ($8+$9 <= 0 && $8+$9 > -6 && $10-$11 <= 5 && (($6+$7)/$5 <= $4/100 || ($6+$7 <= 1))) {print $0"\tPASS"} else {print $0"\tFAIL"}}' | awk '{print $NF}')
# 		output=$(echo -e "$i" "$TE" "$strand" "$DIV" "$length" "$MM" "$gaps" "$(($Lend-$Llen-1))" "$offsetL" "$((Rstart+$offsetR))" "$offsetR" "$TSDs" | sed 's/\n//g;s/ /\t/g' | awk '{if ($8+$9 <= 0 && $8+$9 > -6 && $10-$11 <= 5 && (($6+$7)/$5 <= $4/100 || ($6+$7 <= 1))) {print $0"\tPASS"} else {print $0"\tFAIL"}}')
# 		# debug
# 		echo ""
# 		echo "$test"
# 		echo ""
# 	fi
# 	# poly-A/T elongation
# 	if [[ "$test" == "PASS" ]]
# 		then
# 			#echo -e "SVname\tTEname\tStrand\tDiv\tAlnLen\tMM\tGaps\t5P_TSD_end\t5P_offset\t3P_TSD_start\t3P_offset\t5P_TSD\t3P_TSD"
# 			#echo "$output"
# 			if (( $offsetL > 0 && $Lend == $Lshort ))
# 				then
# 				echo "3' end: can extend poly T"
# 				# count how many T are left in the L (5') sequence (in the poly-T)
# 				# take the L short sequence
# 				LSseq=$(tail -n 1 L.short.fasta)
# 				# grep it out of the long sequence and count the T
# 				Lleft=$(sed "s/$LSseq//g" <(tail -n 1 L.fasta))
# 				LLcount=$(echo ${#Lleft})
# 				# count how many T are left in the R (3') sequence on the 3' of the last matching residue
# 				Rleft=$(sed "s/$RTSD/\t/g" <(tail -n 1 R.fasta) | cut -f 2 | grep "^[T]*T" -o)
# 				RLcount=$(echo ${#Rleft})
# 				# now find the smallest and at those Ts to the TSDs
# 				TtoAdd=$(echo $((LLcount<RLcount ? LLcount : RLcount)))
# 				echo "we could add $TtoAdd T to the 3' of the TSDs"
# 				# adds the T to the TSDs
# 				LTSDext=$(paste -d "" <(echo $LTSD) <(awk -v Ttime=${TtoAdd} -v LTSD=${LTSD} 'BEGIN {for(c=0;c<Ttime;c++) printf "T"; printf "\n"}'))
# 				RTSDext=$(paste -d "" <(echo $RTSD) <(awk -v Ttime=${TtoAdd} -v LTSD=${LTSD} 'BEGIN {for(c=0;c<Ttime;c++) printf "T"; printf "\n"}'))
# 			else
# 				echo "3' end: nothing to extend"
# 				LTSDext=${LTSD}
# 				RTSDext=${RTSD}
# 			fi
# 			if (( $offsetR > 0 && $Rstart == 1 ))
# 				then
# 				echo "5' end: can extend poly A"
# 				# count how many A are left in the R (3') sequence (in the poly-A)
# 				# take the R short sequence
# 				RSseq=$(tail -n 1 R.short.fasta)
# 				# grep it out of the long sequence and count the A
# 				Rleft2=$(sed "s/$RSseq//g" <(tail -n 1 R.fasta))
# 				LRcount=$(echo ${#Rleft2})
# 				# count how many A are left in the L (5') sequence on the 5' of the last matching residue
# 				Lleft2=$(sed "s/$LTSD/\t/g" <(tail -n 1 L.fasta) | cut -f 1 | grep "[A]*A\b" -o)
# 				RLcount2=$(echo ${#Lleft2})
# 				AtoAdd=$(echo $((LRcount<RLcount2 ? LRcount : RLcount2)))
# 				echo "we could add $AtoAdd A to the 5' of the TSDs"
# 				# adds the A to the TSDs
# 				#LTSDextF=$(awk -v Ttime=${AtoAdd} -v LTSD=${LTSDext} 'BEGIN{for(c=0;c<Ttime;c++) printf "A"; printf LTSD; printf "\n"}')
# 				LTSDextF=$(paste -d "" <(awk -v Ttime=${AtoAdd} -v LTSD=${LTSDext} 'BEGIN {for(c=0;c<Ttime;c++) printf "A"}') <(echo $LTSDext)) 
# 				#RTSDextF=$(awk -v Ttime=${AtoAdd} -v RTSD=${RTSDext} 'BEGIN{for(c=0;c<Ttime;c++) printf "A"; printf RTSD; printf "\n"}')
# 				RTSDextF=$(paste -d "" <(awk -v Ttime=${AtoAdd} -v LTSD=${RTSDext} 'BEGIN {for(c=0;c<Ttime;c++) printf "A"}') <(echo $RTSDext)) 
# 			else
# 				echo "5' end: nothing to extend"
# 				LTSDextF=${LTSDext}
# 				RTSDextF=${RTSDext}
# 			fi
# 			echo -e "SVname\tTEname\tStrand\tDiv\tAlnLen\tMM\tGaps\t5P_TSD_end\t5P_offset\t3P_TSD_start\t3P_offset\t5P_TSD\t3P_TSD"
# 			echo -e "$i" "$TE" "$strand" "$DIV" "$length" "$MM" "$gaps" "$(($Lend-$Llen-1))" "$offsetL" "$((Rstart+$offsetR))" "$offsetR" "$LTSDextF" "$RTSDextF" "PASS" | sed 's/\n//g;s/ /\t/g' | tee -a $name.TSD_summary.txt
# 		else
# 			echo -e "SVname\tTEname\tStrand\tDiv\tAlnLen\tMM\tGaps\t5P_TSD_end\t5P_offset\t3P_TSD_start\t3P_offset\t5P_TSD\t3P_TSD"
# 			echo "$output" | tee -a $name.TSD_summary.txt
# 	fi # close loop for elongation
# fi # close loop that check if one end is poly-A or poly-T and skip
# close loop and feed it with each selected (1 TE hit) SV name
done <<< "$indels" 
}

# exec and print output according to verbose option
if [[ ${VERBOSE} == "V" ]]
then
	tsdfind | tee $name.TSD_full_log.txt
else
	tsdfind > $name.TSD_full_log.txt
fi