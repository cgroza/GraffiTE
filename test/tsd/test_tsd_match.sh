#!/usr/bin/env bash
# Regression test for how TSD_Match_v2.sh picks the best TSD candidate.
#
# The search keeps every exact match of 4 to 20 bp between the two junction
# fragments, L = [5' flank][SV 5' end] and R = [SV 3' end][3' flank], scores
# each by its distance from the junction, and keeps the best. A TSD can sit in
# either of two layouts:
#
#   [TSD---SV---]TSD   copies start the SV and start the 3' flank: both start
#                      at column WIN+1 of their fragment
#   TSD[---SV---TSD]   copies end the 5' flank and end the SV: both end at
#                      column WIN
#
# Start offsets used to be taken from WIN rather than WIN+1, so a real TSD in
# the first layout scored 1 instead of 0. It then lost to any shorter match
# scoring 0.5, and TPRT insertions often have one: a TSD opening with AAAA (the
# L1 endonuclease cuts at TTTT/AA) sits against the element's poly(A) tail,
# and the tail's last A plus the TSD's first three As make a one-base-shifted
# AAAA. In a sample of 150 insertions the search had called at 4 bp, 83% had a
# longer candidate that passed and lost, and in 68% the winner was this
# shifted AAAA.
#
# With only the offset fixed, a real TSD one base off the junction loses to a
# short homopolymer sitting exactly on it. So candidates scoring 1.5 or less
# are ties and the longest wins. At 1.5 the two copies can sit up to three
# bases off the junction between them.
#
# Case 1 is the insertion above. Case 2 is the end-anchored layout, which scored
# 0 before and must still. Case 3 is the reverse fault, taken from chm13v2.0: the
# Alu at chr1:116,545,968, whose 17 bp TSD opens with As that run into the
# poly(A) tail, so its junction is ambiguous by a base. Case 4 is the Alu at
# chr16:11,028,501, whose two TSD copies sit two bases and one base off the
# junction, a score of 1.5.
#
# Candidates shorter than 6 bp must also score 0.5 or less (TSD_SHORT and
# TSD_SHORT_SCORE in TSD_Match_v2.sh). Case 5 is the Alu at chr12:66,597,281:
# its two TSD copies differ at column 32, so the exact match is 18 bp starting
# two bases in on both sides (score 2), and it lost to a 4 bp match scoring
# 1.5. Case 6 is synthetic: L holds only C and G, R only A and T, except for
# one planted 4 bp word 2.5 bases off the junction. It passed before and must
# now fail.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export PATH="$(cd ../../bin && pwd):${PATH}"

fail=0
chk(){ if [[ "$2" == "$3" ]]; then echo "  [ ok ] $1"; else
       echo "  [FAIL] $1: got '$2', want '$3'"; fail=1; fi; }

work=$(mktemp -d); trap 'rm -rf "${work}"' EXIT
cd "${work}"

TSD=AAAAGTTATGTCT     # 13 bp, opens with AAAA
# case 1, [TSD---SV---]TSD, the SV ending in a poly(A) tail
L1_FLANK=GCTAGCTAGGCTCCAGTGCATGCCTGACGT          # 30 bp, ends in T
L1_SV=${TSD}GGCCGGGCGCGGTGGCTCA                   # 30 bp: TSD + element start
R1_SV=CTGTAATCCCAGCACTTTGGAAAAAAAAAA               # 30 bp: element end + poly(A)
R1_FLANK=${TSD}CAGGTCTGCATCACC                     # 30 bp: TSD + flank
# case 2, TSD[---SV---TSD]
TSD2=GATTCAGTTC
L2_FLANK=CCTGATGCATGCAGCATCAC${TSD2}             # 30 bp ending in the TSD
L2_SV=GGCCGGGCGCGGTGGCTCACGCCTGTAATC
R2_SV=CCAGCACTTTGGGAGGCCGA${TSD2}                # 30 bp ending in the TSD
R2_FLANK=CGTCAAGCTTGCATGCAGGTCTGCCTTAAC
# case 3, chr1-116545968-DEL-321_2222 as the search sees it
TSD3=AAAAATTGTAACTTGTT
L3_FLANK=tagctTTGTATCTATAATGGGGTTTCTGTT
L3_SV=AAAAATTGTAacttgttggccgggcgcggt
R3_SV=ctgtctcaaaaaaaaaaaaaaaaaaaaaaa
R3_FLANK=aaaattgtaacttgttttacatttcaaaat
# case 4, chr16-11028501-DEL-315_39892
TSD4=TAAAAAATAAAAAAGAA
L4_FLANK=ataaccatgtacaattataatgcatccatt
L4_SV=aaaaaataaaaaagaaggccgggcgcggtg
R4_SV=agactccgtctcaaaaaaaaaaaaaaaata
R4_FLANK=aaaaataaaaaagaaaaagatagcatTAAA
# case 5, chr12-66597281-INS-337_21819
TSD5=AAGAAATGCATATTAAAG
L5_FLANK=tcatgaaagaaaagatgttcaacttcactc
L5_SV=AAAAGAAATGCATATTAAAGGCCGGGCGCG
R5_SV=AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
R5_FLANK=ataagaaatgcatattaaagccacaggaaa
# case 6, synthetic: the only shared word is GATC, at L 33-36 and R 34-37
L6_FLANK=CGGCGCCGCGGCGCCGGCGCGGCCGCGCGG
L6_SV=CCGATCGGCCGCGGCGCCGGCGCGCCGGCG
R6_SV=TTATAATTTAATATTATAATTTATATTAAT
R6_FLANK=ATTGATCATATTTAATATTTATATTTAATA

cat > flanking_sequences.fasta <<EOF
>case5__L
${L5_FLANK}
>case5__R
${R5_FLANK}
>case6__L
${L6_FLANK}
>case6__R
${R6_FLANK}
>case1__L
${L1_FLANK}
>case1__R
${R1_FLANK}
>case2__L
${L2_FLANK}
>case2__R
${R2_FLANK}
>case3__L
${L3_FLANK}
>case3__R
${R3_FLANK}
>case4__L
${L4_FLANK}
>case4__R
${R4_FLANK}
EOF
cat > SV_sequences_L_R_trimmed_WIN.fa <<EOF
>case5__L
${L5_SV}
>case5__R
${R5_SV}
>case6__L
${L6_SV}
>case6__R
${R6_SV}
>case1__L
${L1_SV}
>case1__R
${R1_SV}
>case2__L
${L2_SV}
>case2__R
${R2_SV}
>case3__L
${L3_SV}
>case3__R
${R3_SV}
>case4__L
${L4_SV}
>case4__R
${R4_SV}
EOF
printf 'case1\ncase2\ncase3\ncase4\ncase5\ncase6\n' > indels.txt

TSD_Match_v2.sh SV_sequences_L_R_trimmed_WIN.fa flanking_sequences.fasta indels.txt 30 > /dev/null 2>&1
S=case1.TSD_summary.txt
# columns counted from the end, as tsd_annotate_vcf.sh reads them:
# status, 3' copy, 5' copy, score
# uppercased, as tsd_annotate_vcf.sh writes them: the matcher keeps the
# reference's soft-masking
field(){ awk -F'\t' -v id="$1" -v k="$2" '$1 == id { print toupper($(NF - k)) }' "${S}"; }

echo "case 1: start-anchored TSD against a poly(A) tail"
chk "reported TSD is the 13 bp duplication" "$(field case1 2)" "${TSD}"
chk "it scores 0"                            "$(field case1 3)" "0"
chk "it passes"                              "$(field case1 0)" "PASS"
echo "case 2: end-anchored TSD"
chk "reported TSD is the 10 bp duplication"  "$(field case2 2)" "${TSD2}"
chk "it scores 0"                            "$(field case2 3)" "0"
chk "it passes"                              "$(field case2 0)" "PASS"
echo "case 3: TSD one base off the junction, against a poly(A) tail"
chk "reported TSD is the 17 bp duplication"  "$(field case3 2)" "${TSD3}"
chk "it passes"                              "$(field case3 0)" "PASS"
echo "case 4: TSD copies two bases and one base off the junction"
chk "reported TSD is the 17 bp duplication"  "$(field case4 2)" "${TSD4}"
chk "it passes"                              "$(field case4 0)" "PASS"
echo "case 5: 18 bp TSD two bases off, against a 4 bp match scoring 1.5"
chk "reported TSD is the 18 bp duplication"  "$(field case5 2)" "${TSD5}"
chk "it passes"                              "$(field case5 0)" "PASS"
echo "case 6: a lone 4 bp match off the junction"
chk "it fails"                               "$(field case6 0)" "FAIL"

exit ${fail}
