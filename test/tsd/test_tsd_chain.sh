#!/usr/bin/env bash
# Regression test for the TSD chain: prepTSD.sh, TSD_Match_v2.sh,
# tsd_annotate_vcf.sh and add_polyA.py, from a VCF and a reference to INFO/TSD
# and INFO/polyA in the VCF.
#
# The chain could fail without an error. bedtools getfasta cannot read a gzip
# reference, and prepTSD.sh went on with an empty flank file, so the matcher
# compared each SV's two ends against each other. A missing exact_match.py
# read as no hit on every variant. On macOS the chain failed on a plain
# reference too, because `paste -d ""` is GNU-only. And add_polyA.py read the
# two-copy TSD value as one string, so it never trimmed the TSD before scanning
# for a tail.
#
# Two insertions and a deletion with a planted 8 bp duplication on a synthetic
# contig, run against the reference as plain FASTA, gzip and BGZF. One
# insertion sits closer to the contig start than the window is wide.
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

python3 - "$tmp" <<'EOF'
import random, sys, os
tmp = sys.argv[1]
random.seed(11)
def rnd(n): return ''.join(random.choice('ACGT') for _ in range(n))
TSD = 'GATTACAG'

seq = list(rnd(400))
def plant(pos):                        # the duplication occupies pos-7..pos, 1-based
    seq[pos-8:pos] = list(TSD)
    # The matcher extends an exact match as far as it goes, so the bases on
    # either side of the planted copies must differ between the two copies:
    # C before and G after in the reference, G before and C after in the SV.
    seq[pos-9] = 'C'
    seq[pos] = 'G'
    seq[pos+1] = 'G'
plant(150)                             # insertion A
plant(12)                              # insertion B, inside the first 30 bp
plant(250)                             # deletion
te = 'C' + rnd(58) + 'G'
del_seq = te + TSD                     # reference reads TSD [te TSD] after this
seq = ''.join(seq[:250]) + del_seq + ''.join(seq[250:])

with open(os.path.join(tmp, 'ref.fa'), 'w') as fh:
    fh.write('>t1\n')
    for i in range(0, len(seq), 60):
        fh.write(seq[i:i+60] + '\n')

ins_a = 'C' + rnd(59) + 'A' * 15 + TSD # polyA tail, then the 3' copy of the TSD
ins_b = 'C' + rnd(48) + 'G' + TSD
recs = [
  ('t1', 12,  'insB', seq[11],           seq[11] + ins_b),
  ('t1', 150, 'insA', seq[149],          seq[149] + ins_a),
  ('t1', 250, 'del1', seq[249] + del_seq, seq[249]),
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
    && add_polyA.py pangenome.vcf -o pangenome.polyA.vcf )
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
  [[ $rc -eq 0 ]] || { cat "$w/prep.log" "$w/annot.log" 2>/dev/null; continue; }

  chk "$enc: two flanks per indel" "$(grep -c '^>' "$w/flanking_sequences.fasta")" "6"
  chk "$enc: flank clamped at the contig start" \
      "$(grep -A1 '^>insB__L$' "$w/flanking_sequences.fasta" | tail -1 | tr -d '\n' | wc -c | tr -d ' ')" "12"
  chk "$enc: one summary row per indel" "$(wc -l < "$w/TSD_summary.txt" | tr -d ' ')" "3"
  chk "$enc: insA passes with the planted TSD" \
      "$(awk -F'\t' '$1=="insA"{print $(NF-2)","$(NF-1)","$NF}' "$w/TSD_summary.txt")" "GATTACAG,GATTACAG,PASS"
  chk "$enc: del1 passes with the planted TSD" \
      "$(awk -F'\t' '$1=="del1"{print $(NF-2)","$(NF-1)","$NF}' "$w/TSD_summary.txt")" "GATTACAG,GATTACAG,PASS"
  chk "$enc: TSD is declared in the header" "$(grep -c '^##INFO=<ID=TSD,' "$w/pangenome.vcf")" "1"
  chk "$enc: insA carries INFO/TSD" \
      "$(bcftools query -i 'ID="insA"' -f '%INFO/TSD\n' "$w/pangenome.vcf")" "GATTACAG,GATTACAG"
  chk "$enc: del1 carries INFO/TSD" \
      "$(bcftools query -i 'ID="del1"' -f '%INFO/TSD\n' "$w/pangenome.vcf")" "GATTACAG,GATTACAG"
  chk "$enc: a FAIL row leaves TSD unset" \
      "$(bcftools query -i 'ID="insB"' -f '%INFO/TSD\n' "$w/pangenome.vcf")" "."
  chk "$enc: polyA is found once the TSD is trimmed" \
      "$(bcftools query -i 'ID="insA"' -f '%INFO/polyA\n' "$w/pangenome.polyA.vcf")" "TRUE"
done

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

if [[ $fail -eq 0 ]]; then echo PASS; else echo FAIL; exit 1; fi
