#!/usr/bin/env bash
# Regression test for the TSD chain: prepTSD.sh, TSD_Match_v2.sh,
# tsd_annotate_vcf.sh, fix_vcf.py and add_polyA.py, from a VCF and a reference
# to INFO/TSD and INFO/polyA in the VCF, in the order concat_repeatmask runs
# them.
#
# The chain could fail without an error. bedtools getfasta cannot read a gzip
# reference, and prepTSD.sh went on with an empty flank file, so the matcher
# compared each SV's two ends against each other. A missing exact_match.py
# read as no hit on every variant. On macOS the chain failed on a plain
# reference too, because `paste -d ""` is GNU-only. And add_polyA.py read the
# two-copy TSD value as one string, so it never trimmed the TSD before scanning
# for a tail.
#
# 884afa8 fixed that read, but in a pipeline run add_polyA.py still trimmed
# nothing, and this test passed because it skipped fix_vcf.py. With INFO/TSD
# declared Number=1, fix_vcf.py's vcfpy writer turned the comma into %2C, and
# add_polyA.py's split on ',' returned one piece again. The polyA call on insA
# is TRUE with or without the trim, since its 15 bp tail passes with the TSD
# left on. insC's 10 bp tail passes only once its TSD is trimmed.
#
# Three insertions and a deletion with a planted 8 bp duplication on a
# synthetic contig, run against the reference as plain FASTA, gzip and BGZF.
# One insertion sits closer to the contig start than the window is wide.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
export PATH="$(cd ../../bin && pwd):$PATH"

for t in samtools bcftools bgzip tabix; do
  command -v $t >/dev/null || { echo "  [skip] $t not on PATH"; exit 0; }
done
python3 -c 'import pysam' 2>/dev/null || { echo "  [skip] pysam not importable"; exit 0; }

fail=0
chk(){ if [[ "$2" == "$3" ]]; then echo "  [ ok ] $1"; else
       echo "  [FAIL] $1: got '$2', want '$3'"; fail=1; fi; }

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

# Without vcfpy the chain runs without fix_vcf.py, as it did when the %2C bug
# got past this test, and the test prints a [skip] line saying so.
# fix_vcf.py names /usr/bin/python3, which may not be the python3 that has
# vcfpy; a shim first on PATH redirects it, as in
# test/truvari_merge/test_multi_input.sh.
if python3 -c 'import vcfpy' 2>/dev/null; then
  have_vcfpy=1
  mkdir "$tmp/shim"
  printf '#!/usr/bin/env bash\nexec python3 %s "$@"\n' "$(command -v fix_vcf.py)" > "$tmp/shim/fix_vcf.py"
  chmod +x "$tmp/shim/fix_vcf.py"
  export PATH="$tmp/shim:$PATH"
else
  have_vcfpy=0
  echo "  [skip] vcfpy not importable: fix_vcf.py is left out and its checks do not run"
fi

python3 - "$tmp" <<'EOF'
import random, sys, os
tmp = sys.argv[1]
random.seed(11)
def rnd(n): return ''.join(random.choice('ACGT') for _ in range(n))
TSD = 'GATTACAG'
TSD_C = 'CGTCTGAC'                     # no A in its first three bases, see insC

seq = list(rnd(400))
def plant(pos, tsd=TSD):               # the duplication occupies pos-7..pos, 1-based
    seq[pos-8:pos] = list(tsd)
    # The matcher extends an exact match as far as it goes, so the bases on
    # either side of the planted copies must differ between the two copies:
    # C before and G after in the reference, G before and C after in the SV.
    seq[pos-9] = 'C'
    seq[pos] = 'G'
    seq[pos+1] = 'G'
plant(150)                             # insertion A
plant(12)                              # insertion B, inside the first 30 bp
plant(250)                             # deletion
plant(350, TSD_C)                      # insertion C, at 418 once the deletion is in
te = 'C' + rnd(58) + 'G'
del_seq = te + TSD                     # reference reads TSD [te TSD] after this
seq = ''.join(seq[:250]) + del_seq + ''.join(seq[250:])

with open(os.path.join(tmp, 'ref.fa'), 'w') as fh:
    fh.write('>t1\n')
    for i in range(0, len(seq), 60):
        fh.write(seq[i:i+60] + '\n')

ins_a = 'C' + rnd(59) + 'A' * 15 + TSD # polyA tail, then the 3' copy of the TSD
ins_b = 'C' + rnd(48) + 'G' + TSD
# A 10 bp tail after CG, then the 3' copy. add_polyA.py lets a tail end up to
# 5 bp short of the terminus, so with the TSD left on, every window it scans
# includes the CGT that opens the TSD. The best of them, the tail plus CGT,
# holds 10 A in 13 bp, under the 80% the call needs.
ins_c = 'C' + rnd(59) + 'A' * 10 + TSD_C
recs = [
  ('t1', 12,  'insB', seq[11],           seq[11] + ins_b),
  ('t1', 150, 'insA', seq[149],          seq[149] + ins_a),
  ('t1', 250, 'del1', seq[249] + del_seq, seq[249]),
  ('t1', 418, 'insC', seq[417],          seq[417] + ins_c),
]
with open(os.path.join(tmp, 'genotypes_repmasked_filtered.vcf'), 'w') as fh:
    fh.write('##fileformat=VCFv4.2\n')
    fh.write(f'##contig=<ID=t1,length={len(seq)}>\n')
    fh.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="x">\n')
    fh.write('##INFO=<ID=n_hits,Number=1,Type=Integer,Description="x">\n')
    fh.write('##INFO=<ID=RM_hit_strands,Number=.,Type=String,Description="x">\n')
    fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="x">\n')
    fh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n')
    for c, p, i, r, a in recs:
        t = 'INS' if len(a) > len(r) else 'DEL'
        fh.write(f'{c}\t{p}\t{i}\t{r}\t{a}\t.\tPASS\t'
                 f'SVTYPE={t};n_hits=1;RM_hit_strands=+\tGT\t0/1\n')
EOF

# insB's row is forced to FAIL before annotation. Its flank is 12 bp, and a
# chance 4-mer can still score inside the threshold there, so the search
# outcome is not deterministic; what the test checks on insB is the clamp and
# that a FAIL row leaves the record without a TSD.
run_chain(){                           # $1 workdir, $2 reference basename, $3 window (default 30)
  local win=${3:-30}
  ( cd "$1" \
    && prepTSD.sh "$2" "$win" 1 > prep.log 2>&1 \
    && TSD_Match_v2.sh SV_sequences_L_R_trimmed_WIN.fa flanking_sequences.fasta indels.txt "$win" > /dev/null 2>&1 \
    && cat ./*.TSD_summary.txt > TSD_summary.raw.txt \
    && awk -F'\t' -v OFS='\t' '$1 == "insB" { $NF = "FAIL" } 1' TSD_summary.raw.txt > TSD_summary.txt \
    && tsd_annotate_vcf.sh genotypes_repmasked_filtered.vcf TSD_summary.txt pangenome.vcf > annot.log 2>&1 \
    && fix_step "$2" \
    && add_polyA.py pangenome_nopa.vcf -o pangenome.polyA.vcf )
}

# concat_repeatmask's fix_vcf.py step (module/main.nf:578-586), run in the
# workdir. Like concat_repeatmask, it re-compresses a gzip reference to BGZF
# first so that pysam can open it.
fix_step(){                            # $1 reference basename
  local ref=$1
  if (( ! have_vcfpy )); then cp pangenome.vcf pangenome_nopa.vcf; return; fi
  if [[ $ref == *.gz ]] && ! (file -L "$ref" | grep -q BGZF); then
    gzip -dc "$ref" | bgzip -c > fix_ref.fa.gz || return 1
    ref=fix_ref.fa.gz
  fi
  fix_vcf.py --ref "$ref" --vcf_in pangenome.vcf --vcf_out pangenome_nopa.vcf 2> fix.log
}

for enc in plain gzip bgzf; do
  w="$tmp/$enc"; mkdir -p "$w"; cp "$tmp/genotypes_repmasked_filtered.vcf" "$w/"
  case $enc in
    plain) cp "$tmp/ref.fa" "$w/ref.fa";            ref=ref.fa ;;
    gzip)  gzip  -c "$tmp/ref.fa" > "$w/ref.fa.gz"; ref=ref.fa.gz ;;
    bgzf)  bgzip -c "$tmp/ref.fa" > "$w/ref.fa.gz"; ref=ref.fa.gz ;;
  esac
  rc=0; run_chain "$w" "$ref" || rc=$?
  chk "$enc: chain exits 0" "$rc" "0"
  [[ $rc -eq 0 ]] || { cat "$w/prep.log" "$w/annot.log" "$w/fix.log" 2>/dev/null; continue; }

  chk "$enc: two flanks per indel" "$(grep -c '^>' "$w/flanking_sequences.fasta")" "8"
  chk "$enc: flank clamped at the contig start" \
      "$(grep -A1 '^>insB__L$' "$w/flanking_sequences.fasta" | tail -1 | tr -d '\n' | wc -c | tr -d ' ')" "12"
  chk "$enc: one summary row per indel" "$(wc -l < "$w/TSD_summary.txt" | tr -d ' ')" "4"
  chk "$enc: insA passes with the planted TSD" \
      "$(awk -F'\t' '$1=="insA"{print $(NF-2)","$(NF-1)","$NF}' "$w/TSD_summary.txt")" "GATTACAG,GATTACAG,PASS"
  chk "$enc: del1 passes with the planted TSD" \
      "$(awk -F'\t' '$1=="del1"{print $(NF-2)","$(NF-1)","$NF}' "$w/TSD_summary.txt")" "GATTACAG,GATTACAG,PASS"
  chk "$enc: insC passes with the planted TSD" \
      "$(awk -F'\t' '$1=="insC"{print $(NF-2)","$(NF-1)","$NF}' "$w/TSD_summary.txt")" "CGTCTGAC,CGTCTGAC,PASS"
  chk "$enc: TSD is declared in the header, with two values" \
      "$(grep -c '^##INFO=<ID=TSD,Number=2,' "$w/pangenome.vcf")" "1"
  chk "$enc: insA carries INFO/TSD" \
      "$(bcftools query -i 'ID="insA"' -f '%INFO/TSD\n' "$w/pangenome.vcf")" "GATTACAG,GATTACAG"
  chk "$enc: del1 carries INFO/TSD" \
      "$(bcftools query -i 'ID="del1"' -f '%INFO/TSD\n' "$w/pangenome.vcf")" "GATTACAG,GATTACAG"
  chk "$enc: a FAIL row leaves TSD unset" \
      "$(bcftools query -i 'ID="insB"' -f '%INFO/TSD\n' "$w/pangenome.vcf")" "."
  if (( have_vcfpy )); then
    chk "$enc: fix_vcf.py writes the TSD comma as a comma" \
        "$(awk -F'\t' '$3=="insA"{print $8}' "$w/pangenome_nopa.vcf" | tr ';' '\n' | grep '^TSD=')" "TSD=GATTACAG,GATTACAG"
  fi
  chk "$enc: polyA is found on insA" \
      "$(bcftools query -i 'ID="insA"' -f '%INFO/polyA\n' "$w/pangenome.polyA.vcf")" "TRUE"
  chk "$enc: polyA is found on insC once its TSD is trimmed" \
      "$(bcftools query -i 'ID="insC"' -f '%INFO/polyA\n' "$w/pangenome.polyA.vcf")" "TRUE"
done

# Two checks on insC, on the plain run. With its TSD removed from INFO,
# add_polyA.py does not call the tail, so the TRUE above comes from the trim.
# And add_polyA.py gives the same call on a pangenome.vcf from before TSD was
# Number=2, where the comma is %2C.
w="$tmp/plain"
awk -F'\t' -v OFS='\t' '$3 == "insC" { sub(/;?TSD=[^;]*/, "", $8) } 1' \
    "$w/pangenome_nopa.vcf" > "$w/no_tsd.vcf"
add_polyA.py "$w/no_tsd.vcf" -o "$w/no_tsd.polyA.vcf"
chk "insC without its TSD: the tail alone is not called" \
    "$(bcftools query -i 'ID="insC"' -f '%INFO/polyA\n' "$w/no_tsd.polyA.vcf")" "FALSE"
sed -e 's/^##INFO=<ID=TSD,Number=2,/##INFO=<ID=TSD,Number=1,/' \
    -e '/^#/!s/\(TSD=[ACGT]*\),/\1%2C/' "$w/pangenome_nopa.vcf" > "$w/encoded.vcf"
chk "encoded fixture: three TSDs carry %2C" "$(grep -c 'TSD=[ACGT]*%2C' "$w/encoded.vcf")" "3"
add_polyA.py "$w/encoded.vcf" -o "$w/encoded.polyA.vcf"
chk "a %2C TSD from an earlier run is still trimmed" \
    "$(bcftools query -i 'ID="insC"' -f '%INFO/polyA\n' "$w/encoded.polyA.vcf")" "TRUE"

# A window other than 30. The matcher scored hit offsets against a literal 30,
# so at 40 a snug TSD scored 3 instead of 0, and a wider window pushes it past
# the PASS threshold.
w="$tmp/win40"; mkdir -p "$w"; cp "$tmp/genotypes_repmasked_filtered.vcf" "$tmp/ref.fa" "$w/"
rc=0; run_chain "$w" ref.fa 40 || rc=$?
chk "win40: chain exits 0" "$rc" "0"
chk "win40: flank is 40 bp" \
    "$(grep -A1 '^>insA__L$' "$w/flanking_sequences.fasta" | tail -1 | tr -d '\n' | wc -c | tr -d ' ')" "40"
chk "win40: insA passes with the planted TSD" \
    "$(awk -F'\t' '$1=="insA"{print $(NF-2)","$(NF-1)","$NF}' "$w/TSD_summary.txt")" "GATTACAG,GATTACAG,PASS"
chk "win40: del1 passes with the planted TSD" \
    "$(awk -F'\t' '$1=="del1"{print $(NF-2)","$(NF-1)","$NF}' "$w/TSD_summary.txt")" "GATTACAG,GATTACAG,PASS"
chk "win40: insA scores 0 against the junction" \
    "$(awk -F'\t' '$1=="insA"{print $(NF-3)}' "$w/TSD_summary.txt")" "0"

# Wrong reference: the contig is not there. This has to stop the run rather
# than write an empty flank file.
w="$tmp/wrong"; mkdir -p "$w"; cp "$tmp/genotypes_repmasked_filtered.vcf" "$w/"
sed 's/^>t1$/>t2/' "$tmp/ref.fa" > "$w/ref.fa"
rc=0; ( cd "$w" && prepTSD.sh ref.fa 30 1 > prep.log 2>&1 ) || rc=$?
chk "wrong reference: prepTSD.sh exits non-zero" "$([[ $rc -ne 0 ]] && echo yes || echo no)" "yes"
chk "wrong reference: the message names the contig" "$(grep -c 'contig t1' "$w/prep.log")" "1"

# A chromosome with no records after filtering is not an error.
w="$tmp/empty"; mkdir -p "$w"; cp "$tmp/ref.fa" "$w/"
grep '^#' "$tmp/genotypes_repmasked_filtered.vcf" > "$w/genotypes_repmasked_filtered.vcf"
rc=0; ( cd "$w" && prepTSD.sh ref.fa 30 1 > prep.log 2>&1 ) || rc=$?
chk "empty VCF: prepTSD.sh exits 0" "$rc" "0"
chk "empty VCF: no indel listed" "$(wc -c < "$w/indels.txt" | tr -d ' ')" "0"

# No INFO/SVTYPE, which is what the multi-VCF path produces. truvari_merge
# strips every upstream INFO field before the collapse (module/main.nf:216) and
# puts only SVLEN back, so a run with two or more caller VCFs reaches add_polyA.py
# with no SVTYPE at all. Keyed on SVTYPE alone it scanned an empty string and
# answered FALSE for every record, and the polyA="TRUE" clause of the --human
# filter then dropped every Alu, L1 and SVA. insA is the same record as above and
# still carries a 15 bp A tail behind its TSD.
w="$tmp/nosvtype"; mkdir -p "$w"; cp "$tmp/ref.fa" "$w/"
grep -v '^##INFO=<ID=SVTYPE' "$tmp/genotypes_repmasked_filtered.vcf" \
  | awk -F'\t' -v OFS='\t' '/^#/ {print; next} {sub(/SVTYPE=[^;]*;/, "", $8); print}' \
  > "$w/genotypes_repmasked_filtered.vcf"
chk "no-SVTYPE fixture really has none" \
    "$(grep -c 'SVTYPE' "$w/genotypes_repmasked_filtered.vcf")" "0"
rc=0; run_chain "$w" ref.fa || rc=$?
chk "no SVTYPE: chain exits 0" "$rc" "0"
chk "no SVTYPE: polyA still found on insA" \
    "$(bcftools query -i 'ID="insA"' -f '%INFO/polyA\n' "$w/pangenome.polyA.vcf")" "TRUE"
chk "no SVTYPE: insB has no tail and stays FALSE" \
    "$(bcftools query -i 'ID="insB"' -f '%INFO/polyA\n' "$w/pangenome.polyA.vcf")" "FALSE"
chk "no SVTYPE: the deletion is still read from REF" \
    "$(bcftools query -i 'ID="del1"' -f '%INFO/polyA\n' "$w/pangenome.polyA.vcf")" "FALSE"

if [[ $fail -eq 0 ]]; then echo PASS; else echo FAIL; exit 1; fi
