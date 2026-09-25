#!/usr/bin/env bash
# Regression test for a RepeatMasker chunk with nothing to annotate (#93).
#
# split_repeatmask cuts the VCF one piece per contig, so a contig whose records
# are all non-indel gave RepeatMasker an empty FASTA. RepeatMasker failed,
# repmask_vcf.sh had no `set -e` and carried on, repeatmasker_dir was left
# empty, and `cp repeatmasker_dir/repeatmasker_dir/*` killed tsd_prep two
# processes later. Turning on `set -e` then turned a second silent failure into
# a loud one: annotate_vcf.R aborted on a .out holding only its three header
# lines, which is what RepeatMasker writes when it finds nothing.
#
# Checked here: the no-indel early exit produces every output repeatmask_VCF
# declares, its VCF still parses under the total_repeat_span filter that runs
# next, and read_rm_custom returns an empty table rather than aborting.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
BIN="$(cd ../../bin && pwd)"
export PATH="${BIN}:$PATH"

command -v bcftools >/dev/null || { echo "  [skip] bcftools not on PATH"; exit 0; }

fail=0
chk(){ if [[ "$2" == "$3" ]]; then echo "  [ ok ] $1"; else
       echo "  [FAIL] $1: got '$2', want '$3'"; fail=1; fi; }

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

cat > "$tmp/no_indels.vcf" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=ctg1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1
ctg1	100	v1	A	T	.	PASS	.	GT	1
ctg1	200	v2	C	G	.	PASS	.	GT	0
EOF

w="$tmp/chunk"; mkdir -p "$w"; cp "$tmp/no_indels.vcf" "$w/"
rc=0
( cd "$w" && repmask_vcf.sh no_indels.vcf genotypes_repmasked.vcf.gz lib.fa > run.log 2>&1 ) || rc=$?
chk "no indel: repmask_vcf.sh exits 0" "$rc" "0"
chk "no indel: it says why" \
    "$(grep -c 'nothing to annotate' "$w/run.log" || true)" "1"
chk "no indel: RepeatMasker was not run" \
    "$(grep -c 'RepeatMasker' "$w/run.log" || true)" "0"

# repeatmask_VCF declares these, and Nextflow fails the task if one is absent.
missing=""
for f in genotypes_repmasked.vcf.gz ultra_out.bed ultra_out.span ultra_out.stats \
         union.bp total_repeat_span.tsv combined.stats vcf_annotation.bak.txt; do
  [[ -e "$w/$f" ]] || missing="${missing} $f"
done
chk "no indel: every declared output exists" "${missing}" ""
chk "no indel: repeatmasker_dir exists" "$([[ -d "$w/repeatmasker_dir" ]] && echo yes)" "yes"

# The line repeatmask_VCF runs straight after the script. bcftools refuses an
# expression naming a tag the header does not declare, so this is the check
# that the early exit still wrote the INFO definitions.
rc=0
( cd "$w" && bcftools view -Ov -o filtered.vcf \
    -i 'INFO/total_repeat_span > 0.80' genotypes_repmasked.vcf.gz ) 2>/dev/null || rc=$?
chk "no indel: the span filter parses" "$rc" "0"
chk "no indel: no record survives it" \
    "$(bcftools view -H "$w/filtered.vcf" | wc -l | tr -d ' ')" "0"
chk "no indel: the two input records are kept upstream of it" \
    "$(bcftools view -H "$w/genotypes_repmasked.vcf.gz" | wc -l | tr -d ' ')" "2"

# A .out with only the three lines read_rm_custom skips: what RepeatMasker
# leaves when it finds nothing, and what the missing-.out guard synthesises.
if command -v Rscript >/dev/null && \
   Rscript -e 'quit(status = !all(c("dplyr","stringr","tibble","readr") %in% rownames(installed.packages())))' 2>/dev/null; then
  printf '   SW   perc perc perc  query\nscore   div. del. ins.  sequence\n\n' > "$tmp/empty.out"
  n=$(Rscript -e '
    suppressMessages({library(dplyr);library(stringr);library(tibble);library(readr)})
    src <- readLines(file.path(commandArgs(TRUE)[1], "annotate_vcf.R"))
    eval(parse(text = paste(src[1:(grep("^option_list", src)[1] - 1)], collapse = "\n")))
    cat(nrow(read_rm_custom(commandArgs(TRUE)[2])))
  ' "$BIN" "$tmp/empty.out" 2>/dev/null)
  chk "header-only .out: read_rm_custom returns no row" "$n" "0"

  # read_rm_custom returning no row was as far as this test used to go, and the
  # crash is further down: summarise() types L1_5PINV from ifelse(), which gives
  # logical(0) when there are no groups, and replace_na() cannot then put "None"
  # into it. Run the whole script, which is where a user meets it.
  if Rscript -e 'quit(status = !all(c("optparse","vcfR","tidyr") %in% rownames(installed.packages())))' 2>/dev/null; then
    cat > "$tmp/one.vcf" <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	ins1	A	ACGTACGTACGT	.	PASS	.
VCF
    ( cd "$tmp" && Rscript "$BIN/annotate_vcf.R" --dotout empty.out --vcf one.vcf \
        --annotation empty_annot.tsv ) >/dev/null 2>&1
    chk "header-only .out: annotate_vcf.R exits 0" "$?" "0"
    chk "header-only .out: it still writes one row" \
        "$(wc -l < "$tmp/empty_annot.tsv" 2>/dev/null | tr -d ' ')" "1"
    chk "header-only .out: the row reads n_hits 0 and None" \
        "$(cut -f6,9,10,13 "$tmp/empty_annot.tsv" 2>/dev/null)" \
        "$(printf '0\tNone\tNone\tNone')"

    # and the same script on a .out that does have a hit, so the fix cannot pass
    # by turning every row into None.
    printf '  300   5.0  0.0  0.0  ins1              2      11    (0) +  AluY           SINE/Alu        1    10   (1)      1\n' \
      >> "$tmp/empty.out"
    ( cd "$tmp" && Rscript "$BIN/annotate_vcf.R" --dotout empty.out --vcf one.vcf \
        --annotation hit_annot.tsv ) >/dev/null 2>&1
    chk "one hit: the row reads n_hits 1 and the family" \
        "$(cut -f6,9,10 "$tmp/hit_annot.tsv" 2>/dev/null)" \
        "$(printf '1\tAluY\tSINE/Alu')"

    # One variant carrying an inverted L1 hit beside an ordinary hit. summarise()
    # used to see a length-2 condition here and write the variant out twice, once
    # with the hit IDs and once with "None".
    printf '   SW   perc perc perc  query\nscore   div. del. ins.  sequence\n\n' > "$tmp/mixed.out"
    printf '  500   2.0  0.0  0.0  ins1    1   10   (0) C  L1HS   LINE/L1   (0)  100   50      1\n' >> "$tmp/mixed.out"
    printf '  500   2.0  0.0  0.0  ins1   11   20   (0) +  L1HS   LINE/L1    200  300  (0)      1\n' >> "$tmp/mixed.out"
    printf '  300   5.0  0.0  0.0  ins1   30   40   (0) +  AluY   SINE/Alu     1   10  (1)      2\n' >> "$tmp/mixed.out"
    ( cd "$tmp" && Rscript "$BIN/annotate_vcf.R" --dotout mixed.out --vcf one.vcf \
        --annotation mixed_annot.tsv ) >/dev/null 2>&1
    chk "inverted L1 beside another hit: one row, not two" \
        "$(wc -l < "$tmp/mixed_annot.tsv" 2>/dev/null | tr -d ' ')" "1"
    chk "inverted L1 beside another hit: both hits and the link ID kept" \
        "$(cut -f6,9,13 "$tmp/mixed_annot.tsv" 2>/dev/null)" \
        "$(printf '2\tL1HS,AluY\t1')"
  else
    echo "  [skip] optparse, vcfR or tidyr not available; full-script check skipped"
  fi
else
  echo "  [skip] Rscript or its libraries not available"
fi

if [[ $fail -eq 0 ]]; then echo PASS; else echo FAIL; exit 1; fi
