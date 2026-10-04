#!/usr/bin/env bash
# Run one or more cells of the synthetic matrix. Each cell is independent
# except that everything after "spine" reuses the spine's outputs, so run
# spine first and let it finish.
#
#   ./run_matrix.sh -l            list the cells
#   ./run_matrix.sh spine         run one
#   ./run_matrix.sh               run $RUNS from INPUTS.env
#
# A failing cell does not stop the others; the per-cell status lands in
# MATRIX.log and the exit code is non-zero if any cell failed.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
[[ -f INPUTS.env ]] || { echo "INPUTS.env not found -- run ./bootstrap.sh first" >&2; exit 1; }
# shellcheck disable=SC1091
source ./INPUTS.env

PROFILE="${PROFILE:-standard}"
CPUS="${CPUS:-16}"
PROJECT="${PROJECT:-cgroza/GraffiTE}"
REVISION="${REVISION:-v1.1dev}"
HANDOUT="$PWD"

CELLS=(spine pangenie graphaligner precomputed longreads bams vcf svs duallib guards epi epi_bam winnowmap tsd_win40 ison vcfs breakscaf hervkref)
describe() { case "$1" in
  spine)        echo "--assemblies x4 + --svs + --human + giraffe + genotyping. Pays discovery, RepeatMasker, TSD and the HERV-K stack once; publishes the graph, the alignments and the RM_dir for everything after it.";;
  pangenie)     echo "--graffite_vcf on the spine's pangenome.vcf, pangenie. The -N left-alignment guard and the allele-drop categories in the audit.";;
  graphaligner) echo "--graffite_vcf, graphaligner, long reads.";;
  precomputed)  echo "--graph and --graph_alignments from the spine. Builds nothing.";;
  longreads)    echo "--longreads, two presets (hifi, ont), sniffles per sample then the joint call.";;
  bams)         echo "--bams. No type column; sample names come from @RG SM.";;
  vcf)          echo "--vcf, the from_vcf=true branch: no tabix, INFO kept, SVLEN not recomputed. Carries the chr4_narrow <INV> that has no indel.";;
  svs)          echo "--svs alone, one VCF: the single-input branch that DOES tabix its input.";;
  duallib)      echo "The spine's inputs against the library that spells the internal region HERVK rather than HERVK-int. The pair-rule records should vanish with no error.";;
  guards)       echo "Launch-time guards only. No container, seconds.";;
  epi)          echo "--epigenomes via a hand-written --lifted CSV. Tier 2.";;
  epi_bam)      echo "--epigenomes with no --lifted: the bamtags_to_BED and lift_epigenome branch, off MM/ML-tagged BAMs.";;
  winnowmap)    echo "--aligner winnowmap over the spine's assemblies. Tier 2.";;
  tsd_win40)    echo "The spine's discovery at --tsd_win 40, no genotyping. Every verdict, and every PASS's TSD, must match the spine's at 30.";;
  ison)         echo "--human and --genotype passed as strings from a params file (\"false\", \"FALSE\", \"\", \"true\", \"True\"). Each has to switch its stages off or on as isOn() reads it.";;
  vcfs)         echo "--graph_method precomputed with --vcfs: the spine's own vg call VCFs handed back. No alignment, no vg_call, and the spine's merged genotypes.";;
  breakscaf)    echo "--break_scaffolds on the four haplotypes: cut at the 120 bp N run, keep the insertions that carry single N bases, and find what discovery without it finds.";;
  hervkref)     echo "--hervk_ref_annotation as a RepeatMasker .out, a BED and an empty BED: no masking in hervk_annotate, the spine's reference states from the first two, none from the third.";;
esac; }

if [[ "${1:-}" == "-l" ]]; then
  for c in "${CELLS[@]}"; do printf '  %-13s %s\n' "$c" "$(describe "$c")"; done; exit 0
fi

WORKDIR="${WORKDIR:?set WORKDIR in INPUTS.env}"
B="$WORKDIR/build"
RUNS_DIR="$WORKDIR/runs"

SEL=("$@"); [[ ${#SEL[@]} -eq 0 ]] && read -r -a SEL <<< "${RUNS:-spine}"

[[ -d "$B" ]] || { echo "no build at $B -- run ./preflight.sh, which builds it" >&2; exit 1; }
mkdir -p "$RUNS_DIR"

export GT_CPUS="$CPUS" GT_MEM_GB="${MEM_GB:-80}"
SIF_ARG=(); [[ -n "${GRAFFITE_SIF:-}" ]] && SIF_ARG=(-with-singularity "$GRAFFITE_SIF")
TMP_ARG=(); [[ -n "${CONTAINER_TMP:-}" ]] && TMP_ARG=(--container_tmp "$CONTAINER_TMP")

nf() {  # nf <cell> <extra args...>
  local cell="$1"; shift
  local out="$RUNS_DIR/$cell"
  # Clear the cell's own directory first. Nextflow refuses to start when
  # -with-trace names a file that exists, so re-running a cell that failed
  # dies on "Trace file already exists" before it reaches the pipeline, and
  # the log then describes the previous run rather than this one.
  rm -rf "$out"
  mkdir -p "$out"
  echo "=== $cell ===" | tee -a "$HANDOUT/MATRIX.log"
  ( set -x
    nextflow -log "$out/nextflow.log" run "$PROJECT" -r "$REVISION" \
      -profile "$PROFILE" -c "$HANDOUT/local.config" \
      --reference "$B/ref/synth.fa" \
      --TE_library "$B/lib/synth_TE.fasta" \
      --cores "$CPUS" \
      --out "$out" \
      "${SIF_ARG[@]}" "${TMP_ARG[@]}" \
      -with-report "$out/nextflow_report.html" \
      -with-trace  "$out/nextflow_trace.txt" \
      "$@"
  ) >"$out/run.log" 2>&1
  local rc=$?
  ( cd "$out" && find . -type f | sort ) > "$out/published.txt" 2>/dev/null
  echo "$cell rc=$rc" | tee -a "$HANDOUT/MATRIX.log"
  return $rc
}

SPINE="$RUNS_DIR/spine"
fail=0
for cell in "${SEL[@]}"; do
  case "$cell" in
  spine)
    nf spine --assemblies "$B/assemblies.csv" --svs "$B/svs.csv" \
       --human --genotype_with "$B/reads.csv" --graph_method giraffe \
       --repeatmasker_memory "${REPEATMASKER_MEMORY:-16G}" || fail=1 ;;
  pangenie)
    nf pangenie --graffite_vcf "$SPINE/3_TSD_search/pangenome.vcf" \
       --genotype_with "$B/reads.csv" --graph_method pangenie \
       --pangenie_memory "${PANGENIE_MEMORY:-40G}" --pangenie_time "${PANGENIE_TIME:-4h}" || fail=1 ;;
  graphaligner)
    nf graphaligner --graffite_vcf "$SPINE/3_TSD_search/pangenome.vcf" \
       --genotype_with "$B/reads_long.csv" --graph_method graphaligner || fail=1 ;;
  precomputed)
    # main.nf:54 requires --graph AND one of --vcfs / --graph_alignments. The cell
    # as shipped passed only --graph, so it stopped at the launch guard every time
    # (measured: rc=1 in 6 s). run_matrix.sh -l already describes this cell as
    # "--graph and --graph_alignments from the spine"; the CSV is built here from
    # what the spine published, since the generator does not write one.
    GA="$RUNS_DIR/precomputed_alignments.csv"
    { echo "sample,gaf,pack"
      for g in "$SPINE"/GraffiTE_alignments/*.gaf.gz; do
        [ -e "$g" ] || continue
        s=$(basename "$g" .gaf.gz)
        [ -f "$SPINE/GraffiTE_alignments/$s.pack" ] && \
          echo "$s,$g,$SPINE/GraffiTE_alignments/$s.pack"
      done
    } > "$GA"
    nf precomputed --graffite_vcf "$SPINE/3_TSD_search/pangenome.vcf" \
       --graph_method precomputed --graph "$SPINE/GraffiTE_graph/index" \
       --graph_alignments "$GA" \
       --genotype_with "$B/reads.csv" || fail=1 ;;
  longreads)
    nf longreads --longreads "$B/longreads.csv" --genotype false || fail=1 ;;
  bams)
    nf bams --bams "$B/bams.csv" --genotype false || fail=1 ;;
  vcf)
    nf vcf --vcf "$B/vcf/merged.vcf.gz" --genotype false || fail=1 ;;
  svs)
    nf svs --svs "$B/svs_one.csv" --genotype false || fail=1 ;;
  duallib)
    # Same inputs, the other spelling of the internal region. This arm calls
    # nextflow itself rather than nf(), so it clears its directory and keeps
    # its exit code on its own. It used to do neither: a rerun stopped at
    # launch on "Trace file already exists", and the rc=$? printed after
    # "|| fail=1" said 0 either way.
    local_out="$RUNS_DIR/duallib"; rm -rf "$local_out"; mkdir -p "$local_out"
    echo "=== duallib ===" | tee -a "$HANDOUT/MATRIX.log"
    ( set -x
      nextflow -log "$local_out/nextflow.log" run "$PROJECT" -r "$REVISION" \
        -profile "$PROFILE" -c "$HANDOUT/local.config" \
        --assemblies "$B/assemblies.csv" --human --genotype false \
        --reference "$B/ref/synth.fa" --TE_library "$B/lib/synth_TE_bare_HERVK.fasta" \
        --cores "$CPUS" --out "$local_out" \
        "${SIF_ARG[@]}" "${TMP_ARG[@]}" -with-trace "$local_out/nextflow_trace.txt"
    ) >"$local_out/run.log" 2>&1
    rc=$?; [[ $rc -eq 0 ]] || fail=1
    echo "duallib rc=$rc" | tee -a "$HANDOUT/MATRIX.log" ;;
  guards)
    bash "$HANDOUT/check_guards.sh" 2>&1 | tee -a "$HANDOUT/MATRIX.log" || fail=1 ;;
  epi)
    nf epi --graffite_vcf "$SPINE/3_TSD_search/pangenome.vcf" --graph_method giraffe \
       --genotype_with "$B/reads.csv" --epigenomes "$B/epigenomes.csv" \
       --lifted "$B/lifted.csv" || fail=1 ;;
  epi_bam)
    # No --lifted, so main.nf:235-239 runs bamtags_to_BED and lift_epigenome
    # instead of reading a CSV. --genotype_with must name BAMs: only a .bam row
    # reaches reads_input_ch.bam (main.nf:185), and the generator writes those
    # with MM/ML tags for tagtobed to read.
    nf epi_bam --graffite_vcf "$SPINE/3_TSD_search/pangenome.vcf" --graph_method giraffe \
       --genotype_with "$B/reads_bam.csv" --epigenomes true || fail=1 ;;
  winnowmap)
    nf winnowmap --assemblies "$B/assemblies.csv" --aligner winnowmap --genotype false || fail=1 ;;
  tsd_win40)
    # The window reaches prepTSD.sh and TSD_Match_v2.sh (module/main.nf:673,688).
    # The same inputs as the spine at a different window should give the same
    # TSD on every SV; check_tsd_win.py compares the two and confirms that the
    # fragments are 80 bp here.
    nf tsd_win40 --assemblies "$B/assemblies.csv" --svs "$B/svs.csv" --human --genotype false \
       --tsd_win 40 --repeatmasker_memory "${REPEATMASKER_MEMORY:-16G}" || fail=1
    python3 "$HANDOUT/check_tsd_win.py" "$SPINE" "$RUNS_DIR/tsd_win40" 40 2>&1 \
      | tee -a "$HANDOUT/MATRIX.log" || fail=1 ;;
  ison)
    # A string reaches isOn() only from a params or config file; the command line
    # turns "false" into a Boolean first. Each run passes one switch that way, on
    # the short --svs input, with -resume so the runs share the masking.
    for v in human_str_false:human:false human_str_FALSE:human:FALSE human_str_empty:human: \
             human_str_true:human:true human_str_True:human:True genotype_str_false:genotype:false; do
      name=${v%%:*}; rest=${v#*:}; key=${rest%%:*}; val=${rest#*:}
      pf="$RUNS_DIR/ison_$name.yaml"; printf '%s: "%s"\n' "$key" "$val" > "$pf"
      args=(--svs "$B/svs_one.csv")
      # genotyping would run if the string were read as on: give it reads
      if [[ $key == genotype ]]; then args+=(--genotype_with "$B/reads.csv" --graph_method giraffe); else args+=(--genotype false); fi
      nf "ison_$name" "${args[@]}" -params-file "$pf" -resume
      echo $? > "$RUNS_DIR/ison_$name.rc"
    done
    bash "$HANDOUT/check_ison.sh" "$RUNS_DIR" 2>&1 | tee -a "$HANDOUT/MATRIX.log" || fail=1 ;;
  vcfs)
    # The spine does not publish its per-sample vg call VCFs, so copy them out of
    # the work directories its trace names. path has to be a glob, as
    # docs/reference/samplesheets.md says, so the .tbi comes along.
    VI="$RUNS_DIR/vcfs_inputs"; mkdir -p "$VI"; VC="$RUNS_DIR/vcfs.csv"
    { echo "sample,path"
      for h in $(awk -F'\t' 'NR>1 && $4 ~ /^vg_call/ {print $2}' "$SPINE/nextflow_trace.txt"); do
        for v in "$HANDOUT"/work/${h}*/*.vcf.gz; do
          [ -f "$v.tbi" ] || continue
          s=$(basename "$v" .vcf.gz); cp "$v" "$v.tbi" "$VI/"
          echo "$s,$VI/$s.vcf.gz*"
        done
      done
    } > "$VC"
    nf vcfs --graffite_vcf "$SPINE/3_TSD_search/pangenome.vcf" \
       --graph_method precomputed --graph "$SPINE/GraffiTE_graph/index" \
       --vcfs "$VC" --genotype_with "$B/reads.csv" || fail=1
    bash "$HANDOUT/check_vcfs.sh" "$SPINE" "$RUNS_DIR/vcfs" 2>&1 | tee -a "$HANDOUT/MATRIX.log" || fail=1 ;;
  breakscaf)
    # The haplotypes carry the reference's 120 bp N run on chr1, and single N
    # bases inside planted insertions. Discovery with --break_scaffolds has to cut
    # at the run only, so it must find what discovery without it finds.
    for mode in nobreak break; do
      args=(--assemblies "$B/assemblies.csv" --genotype false)
      [[ $mode == break ]] && args+=(--break_scaffolds)
      nf "breakscaf_$mode" "${args[@]}" -resume || fail=1
    done
    python3 "$HANDOUT/check_breakscaf.py" "$RUNS_DIR" "$B/assemblies.csv" "$HANDOUT/work" \
      "${BREAK_MIN_GAP:-10}" 2>&1 | tee -a "$HANDOUT/MATRIX.log" || fail=1 ;;
  hervkref)
    # --hervk_ref_annotation gives hervk_annotate a repeat track for the reference
    # in place of its own masking. The test set has none, so mask the whole
    # reference in the pipeline's image with the options hervk_ref_state.py uses,
    # and pass the result as a .out, as a BED, and as a BED with no hit.
    HI="$RUNS_DIR/hervkref_inputs"; mkdir -p "$HI/rm"
    IMG="${GRAFFITE_SIF:-${NXF_SINGULARITY_CACHEDIR:-}/cgroza-graffite-latest.img}"
    "$(command -v apptainer || command -v singularity)" exec -B "$WORKDIR" "$IMG" \
      RepeatMasker -lib "$B/lib/synth_TE.fasta" -s -dir "$HI/rm" -pa "$CPUS" "$B/ref/synth.fa" > "$HI/rm.log" 2>&1
    cp "$HI/rm/synth.fa.out" "$HI/ref.out"
    awk -v OFS='\t' 'NR > 3 {print $5, $6 - 1, $7, $10, $11, ($9 == "C" ? "-" : "+")}' "$HI/ref.out" > "$HI/ref.bed"
    : > "$HI/ref.empty.bed"
    for mode in out bed empty; do
      case $mode in out) f=$HI/ref.out ;; bed) f=$HI/ref.bed ;; empty) f=$HI/ref.empty.bed ;; esac
      nf "hervkref_$mode" --assemblies "$B/assemblies.csv" --svs "$B/svs.csv" --human --genotype false \
         --repeatmasker_memory "${REPEATMASKER_MEMORY:-16G}" --hervk_ref_annotation "$f" -resume || fail=1
    done
    python3 "$HANDOUT/check_hervkref.py" "$RUNS_DIR" "$SPINE" "$HANDOUT/work" 2>&1 | tee -a "$HANDOUT/MATRIX.log" || fail=1 ;;
  *) echo "unknown cell: $cell (see ./run_matrix.sh -l)" >&2; fail=1 ;;
  esac
done
echo "matrix done, fail=$fail"
exit $fail
