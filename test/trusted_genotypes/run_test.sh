#!/usr/bin/env bash
# Regression test for the trusted subset of the genotyped calls
# (module/main.nf, process trusted_genotypes).
#
# The trusted expression reads FILTER, and merge_VCFs transfers INFO but not
# FILTER (-c CHROM,POS,ID,REF,ALT,INFO). The FILTER column of
# GraffiTE.merged.genotypes.vcf.gz is therefore the genotyper's, so running the
# expression against that file answers a different question. trusted_genotypes
# evaluates it on pangenome.vcf and subsets the genotyped calls by the IDs it
# returns; this test pins that difference.
#
# The expression assembly below mirrors trustedFilterFull() in module/main.nf
# and must be kept in sync with it.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${HERE}/../../nextflow.config"
export PATH="$(cd "${HERE}/../../bin" && pwd):$PATH"

for t in bcftools bgzip tabix; do
  command -v $t >/dev/null || { echo "  [skip] $t not on PATH"; exit 0; }
done

fail=0
chk(){ if [[ "$2" == "$3" ]]; then echo "  [ ok ] $1"; else
       echo "  [FAIL] $1: got '$2', want '$3'"; fail=1; fi; }

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

# --- defaults from nextflow.config -------------------------------------------
cfg() { sed -n "s/^[[:space:]]*$1[[:space:]]*=[[:space:]]*[\"']\{0,1\}\([^\"'\/]*\)[\"']\{0,1\}.*/\1/p" "$CONFIG" | head -1 | sed 's/[[:space:]]*$//'; }
MINSV=$(cfg trusted_min_svlen)
MAXTR=$(cfg trusted_max_ultra_span)
chk "min SVLEN read from the config" "$MINSV" "250"
chk "max ULTRA span read from the config" "$MAXTR" "0.6"

# --- expression assembly (mirrors trustedFilter/trustedFilterFull) -----------
TRUSTED="n_hits==1 & abs(SVLEN)>=${MINSV} & (ULTRA_TR_span<${MAXTR} | matching_classes=\"Simple_repeat\") & ((matching_classes!~\"LINE\" & matching_classes!~\"SINE\" & matching_classes!~\"Retroposon\") | polyA=\"TRUE\")"
TRUSTED_FULL="(${TRUSTED}) & FILTER=\"PASS\""

# --- fixture -----------------------------------------------------------------
hdr() {
  echo '##fileformat=VCFv4.2'
  echo '##contig=<ID=chr1,length=100000>'
  echo '##FILTER=<ID=PASS,Description="passed">'
  echo '##FILTER=<ID=LowQual,Description="caller flagged">'
  echo '##FILTER=<ID=lowad,Description="genotyper flagged">'
  echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">'
  echo '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="length">'
  echo '##INFO=<ID=n_hits,Number=1,Type=Integer,Description="hits">'
  echo '##INFO=<ID=matching_classes,Number=.,Type=String,Description="classes">'
  echo '##INFO=<ID=ULTRA_TR_span,Number=1,Type=Float,Description="tandem span">'
  echo '##INFO=<ID=polyA,Number=1,Type=String,Description="polyA">'
  echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n'
}
# id pos filter info                                                      -> trusted?
rows='t1 1000 PASS    SVTYPE=INS;SVLEN=500;n_hits=1;matching_classes=LINE/L1;ULTRA_TR_span=0.10;polyA=TRUE  yes
t2 2000 PASS    SVTYPE=INS;SVLEN=300;n_hits=1;matching_classes=DNA/hAT;ULTRA_TR_span=0.20;polyA=FALSE       yes
u1 3000 PASS    SVTYPE=INS;SVLEN=900;n_hits=2;matching_classes=LINE/L1,LTR/ERVK;ULTRA_TR_span=0.10;polyA=TRUE no
u2 4000 PASS    SVTYPE=INS;SVLEN=100;n_hits=1;matching_classes=LINE/L1;ULTRA_TR_span=0.10;polyA=TRUE        no
u3 5000 PASS    SVTYPE=INS;SVLEN=400;n_hits=1;matching_classes=LINE/L1;ULTRA_TR_span=0.90;polyA=TRUE        no
f1 6000 LowQual SVTYPE=INS;SVLEN=600;n_hits=1;matching_classes=LINE/L1;ULTRA_TR_span=0.10;polyA=TRUE        no'

{ hdr
  while read -r id pos filt info _; do
    printf 'chr1\t%s\t%s\tA\tAGGGG\t.\t%s\t%s\tGT\t1\t0\n' "$pos" "$id" "$filt" "$info"
  done <<< "$rows"
} > "$tmp/pangenome.vcf"

# The merged genotypes: same records, same INFO, but FILTER is the genotyper's
# and the calls are diploid. This is what merge_VCFs produces.
{ hdr
  while read -r id pos _ info _; do
    printf 'chr1\t%s\t%s\tA\tAGGGG\t.\tlowad\t%s\tGT\t0/1\t0/0\n' "$pos" "$id" "$info"
  done <<< "$rows"
} > "$tmp/merged.vcf"
bgzip -c "$tmp/merged.vcf" > "$tmp/merged.vcf.gz" && tabix -p vcf "$tmp/merged.vcf.gz"

# --- what trusted_genotypes does ---------------------------------------------
# bcftools prints a "pass=1 [...]" debug line on stderr for each !~ clause on a
# Number=. field; it is noise, and the expression is evaluated on stdout.
bcftools view -H -i "$TRUSTED_FULL" "$tmp/pangenome.vcf" 2>/dev/null | cut -f3 | sort -u > "$tmp/trusted.ids"
chk "the trusted ids come off pangenome.vcf" \
    "$(tr '\n' ',' < "$tmp/trusted.ids")" "t1,t2,"

bcftools view -i 'ID=@'"$tmp/trusted.ids" -Oz -o "$tmp/trusted.vcf.gz" "$tmp/merged.vcf.gz"
tabix -p vcf "$tmp/trusted.vcf.gz"
chk "the genotyped subset holds those records" \
    "$(bcftools query -f '%ID,' "$tmp/trusted.vcf.gz" | tr -d '\n')" "t1,t2,"
chk "every record in it has one hit" \
    "$(bcftools query -f '%INFO/n_hits\n' "$tmp/trusted.vcf.gz" | sort -u | tr '\n' ',')" "1,"
chk "the sample columns survive" \
    "$(bcftools query -l "$tmp/trusted.vcf.gz" | tr '\n' ',')" "s1,s2,"
chk "the genotypes survive" \
    "$(bcftools query -f '[%GT ]\n' "$tmp/trusted.vcf.gz" | head -1 | tr -d ' ')" "0/10/0"

# --- why it is not run against the merged VCF --------------------------------
direct=$(bcftools view -H -i "$TRUSTED_FULL" "$tmp/merged.vcf.gz" 2>/dev/null | wc -l | tr -d ' ')
chk "running the expression on the merged VCF drops everything" "$direct" "0"

# --- --trusted_ignore_filter ------------------------------------------------
bcftools view -H -i "$TRUSTED" "$tmp/pangenome.vcf" 2>/dev/null | cut -f3 | sort -u > "$tmp/ignore.ids"
chk "ignoring FILTER admits the caller-flagged record" \
    "$(tr '\n' ',' < "$tmp/ignore.ids")" "f1,t1,t2,"

# --- the presence-absence table ----------------------------------------------
if python3 -c 'import sys' 2>/dev/null; then
  bcftools view "$tmp/trusted.vcf.gz" | vcf_to_pa_tsv.py -o "$tmp/pa.tsv"
  chk "the TSV has one row per record plus a header" \
      "$(wc -l < "$tmp/pa.tsv" | tr -d ' ')" "3"
  chk "the TSV carries both samples" \
      "$(head -1 "$tmp/pa.tsv" | tr '\t' '\n' | tail -2 | tr '\n' ',')" "s1,s2,"
  chk "an insertion called 0/1 reads as present" \
      "$(awk -F'\t' 'NR==2{print $NF"|"$(NF-1)}' "$tmp/pa.tsv")" "0|1"
fi

if [[ $fail -eq 0 ]]; then echo PASS; else echo FAIL; exit 1; fi
