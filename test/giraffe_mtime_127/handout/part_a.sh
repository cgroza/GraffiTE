#!/usr/bin/env bash
# Part A: the copy step on its own, no Nextflow.
#
# Builds the giraffe index with make_graph's exact recipe in a scratch
# directory, then hands it to vg giraffe five ways:
#
#   control      in place, no copy
#   cp           cp -fRL index <dest>/      what Nextflow's nxf_fs_copy runs
#   forced_bad   every file, then index.dist last (positive control: must fail)
#   forced_good  index.dist first, then the rest
#   mv           mv -f index <dest>/        stageOutMode 'move'
#
# GNU cp -R copies a directory's entries in inode order (gnulib savedir,
# SAVEDIR_SORT_FASTREAD), so for the cp variant the outcome is predicted from
# the inode numbers the scratch filesystem gave the files. Both the prediction
# and the observation are recorded.
#
# Writes $WORKDIR/partA/{results.tsv,inodes.tsv,meta.txt,logs/}.
set -uo pipefail
export LC_ALL=C   # sort and join below must agree on collation
# shellcheck disable=SC1091
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

for f in "$REFERENCE" "$GRAFFITE_VCF" "$READS_FASTQ"; do
  [[ -f "$f" ]] || { echo "missing input: $f" >&2; exit 1; }
done
ensure_image || exit 1

OUT="$WORKDIR/partA"
rm -rf "$OUT"; mkdir -p "$OUT/logs" "$OUT/dest"
RESULTS="$OUT/results.tsv"
INODES="$OUT/inodes.tsv"
printf 'src_label\tsrc_dir\tsrc_fs\tdest_fs\trep\tvariant\tpredicted\tobserved\trc\tmessage\n' > "$RESULTS"
printf 'src_label\trep\tstage\tfile\tinode\tmtime\tsize\treaddir_pos\n' > "$INODES"

# A few hundred reads are plenty: the age check runs when the indexes load,
# before the first read is mapped.
READS="$OUT/reads_subset.fq.gz"
gzip -cdf "$READS_FASTQ" | head -n 800 | gzip > "$READS"

{
  echo "date      : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host      : $(hostname)"
  echo "kernel    : $(uname -r)"
  echo "image     : $IMG"
  echo "vg        : $(ctr "$OUT" 'vg version' 2>/dev/null | head -1)"
  echo "host cp   : $(cp --version | head -1)"
  echo "host mv   : $(mv --version | head -1)"
  echo "workdir fs: $(fs_type "$WORKDIR")"
  echo "SLURM_TMPDIR: ${SLURM_TMPDIR:-<unset>}"
} > "$OUT/meta.txt"
cat "$OUT/meta.txt"

# record <src_label> <rep> <stage> <dir>: one row per index file, with its
# position in raw readdir order.
record() {
  local label="$1" rep="$2" stage="$3" dir="$4"
  ( cd "$dir/index" && ls -1 -U ) | awk '{print NR "\t" $0}' > "$OUT/.pos"
  find "$dir/index" -maxdepth 1 -type f -printf '%f\t%i\t%T@\t%s\n' | sort \
    | join -t $'\t' -1 1 -2 2 - <(sort -t $'\t' -k2,2 "$OUT/.pos") \
    | awk -v l="$label" -v r="$rep" -v s="$stage" 'BEGIN{OFS="\t"} {print l, r, s, $1, $2, $3, $4, $5}' \
    >> "$INODES"
}

# giraffe_in <dir> <logtag>: sets RC, OBS and MSG.
giraffe_in() {
  local dir="$1" log="$OUT/logs/$2.err"
  ctr "$dir" "vg giraffe --parameter-preset hifi -o gam -t $CPUS --index-basename index/index -f $READS > /dev/null" 2> "$log"
  RC=$?
  MSG="$(grep -m1 -E 'is newer than|error' "$log" | tr '\t' ' ')"
  if [[ $RC -eq 0 ]]; then OBS=pass
  elif grep -q 'is newer than' "$log"; then OBS=fail_mtime
  else OBS=fail_other; fi
}

row() {  # row <label> <src> <dest_fs> <rep> <variant> <predicted>
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$(fs_type "$2")" "$3" "$4" "$5" "$6" "$OBS" "$RC" "$MSG" >> "$RESULTS"
  printf '  %-10s rep %-3s %-12s predicted=%-11s observed=%s\n' "$1" "$4" "$5" "$6" "$OBS"
}

read -r -a SRC_BASES <<< "$SCRATCH_DIRS"
i=0
for base in "${SRC_BASES[@]}"; do
  i=$((i + 1))
  label=$([[ $i -eq 1 ]] && echo nf_scratch || echo "extra$((i - 1))")
  mkdir -p "$base" || { echo "cannot create $base" >&2; continue; }
  echo "== $label: $base ($(fs_type "$base")) =="
  for rep in $(seq 1 "$REPS_A"); do
    # nxf_mktemp, as .command.run creates the task's scratch directory
    src="$(TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX)"
    ln -s "$REFERENCE" "$src/$(basename "$REFERENCE")"
    ln -s "$GRAFFITE_VCF" "$src/$(basename "$GRAFFITE_VCF")"
    # make_graph, giraffe branch, verbatim from module/main.nf
    ctr "$src" "set -e
      mkdir index
      bcftools +setGT $(basename "$GRAFFITE_VCF") -- -t a -n u > unphased.vcf
      vg autoindex --tmp-dir \$PWD  -p index/index -w sr-giraffe -w lr-giraffe -v unphased.vcf -r $(basename "$REFERENCE")
      vg convert --vg-algorithm -f index/index.giraffe.gbz > index/index.gfa
      vg snarls index/index.giraffe.gbz > index/index.pb" \
      > "$OUT/logs/${label}_r${rep}_build.log" 2>&1 \
      || { echo "  index build failed, see logs/${label}_r${rep}_build.log" >&2; rm -rf "$src"; continue; }
    record "$label" "$rep" src "$src"

    # The cp prediction: index.dist must come before both long-read files in
    # inode order, or cp writes it later and vg finds it newer.
    pred=$(awk -F'\t' -v l="$label" -v r="$rep" '
      $1==l && $2==r && $3=="src" {ino[$4]=$5}
      END {d=ino["index.dist"]; m=ino["index.longread.withzip.min"]; z=ino["index.longread.zipcodes"]
           if (d=="" || m=="" || z=="") print "unknown"
           else print ((d+0 < m+0 && d+0 < z+0) ? "pass" : "fail_mtime")}' "$INODES")

    giraffe_in "$src" "${label}_r${rep}_control"
    row "$label" "$src" "$(fs_type "$src")" "$rep" control pass

    d="$OUT/dest/${label}_r${rep}"
    mkdir -p "$d/cp" "$d/forced_bad/index" "$d/forced_good/index" "$d/mv"

    ( cd "$src" && cp -fRL index "$d/cp/" )
    record "$label" "$rep" cp "$d/cp"
    giraffe_in "$d/cp" "${label}_r${rep}_cp"
    row "$label" "$src" "$(fs_type "$d")" "$rep" cp "$pred"

    for f in "$src"/index/*; do [[ "$(basename "$f")" == index.dist ]] || cp -fL "$f" "$d/forced_bad/index/"; done
    sleep 1; cp -fL "$src/index/index.dist" "$d/forced_bad/index/"
    giraffe_in "$d/forced_bad" "${label}_r${rep}_forced_bad"
    row "$label" "$src" "$(fs_type "$d")" "$rep" forced_bad fail_mtime

    cp -fL "$src/index/index.dist" "$d/forced_good/index/"; sleep 1
    for f in "$src"/index/*; do [[ "$(basename "$f")" == index.dist ]] || cp -fL "$f" "$d/forced_good/index/"; done
    giraffe_in "$d/forced_good" "${label}_r${rep}_forced_good"
    row "$label" "$src" "$(fs_type "$d")" "$rep" forced_good pass

    ( cd "$src" && mv -f index "$d/mv/" )
    record "$label" "$rep" mv "$d/mv"
    giraffe_in "$d/mv" "${label}_r${rep}_mv"
    row "$label" "$src" "$(fs_type "$d")" "$rep" mv pass

    rm -rf "$src" "$d"/*/index
  done
done
rm -f "$OUT/.pos"
echo
echo "Part A done: $RESULTS"
