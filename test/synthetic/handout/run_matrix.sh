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
REVISION="${REVISION:-test/synthetic-end-to-end}"
HANDOUT="$PWD"

CELLS=(spine pangenie graphaligner precomputed longreads bams vcf svs duallib guards epi winnowmap)
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
  winnowmap)    echo "--aligner winnowmap over the spine's assemblies. Tier 2.";;
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
    nf precomputed --graffite_vcf "$SPINE/3_TSD_search/pangenome.vcf" \
       --graph_method precomputed --graph "$SPINE/GraffiTE_graph/index" \
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
    # Same inputs, the other spelling of the internal region.
    local_out="$RUNS_DIR/duallib"; mkdir -p "$local_out"
    echo "=== duallib ===" | tee -a "$HANDOUT/MATRIX.log"
    ( set -x
      nextflow -log "$local_out/nextflow.log" run "$PROJECT" -r "$REVISION" \
        -profile "$PROFILE" -c "$HANDOUT/local.config" \
        --assemblies "$B/assemblies.csv" --human --genotype false \
        --reference "$B/ref/synth.fa" --TE_library "$B/lib/synth_TE_bare_HERVK.fasta" \
        --cores "$CPUS" --out "$local_out" \
        "${SIF_ARG[@]}" "${TMP_ARG[@]}" -with-trace "$local_out/nextflow_trace.txt"
    ) >"$local_out/run.log" 2>&1 || fail=1
    echo "duallib rc=$?" | tee -a "$HANDOUT/MATRIX.log" ;;
  guards)
    bash "$HANDOUT/check_guards.sh" 2>&1 | tee -a "$HANDOUT/MATRIX.log" || fail=1 ;;
  epi)
    nf epi --graffite_vcf "$SPINE/3_TSD_search/pangenome.vcf" --graph_method giraffe \
       --genotype_with "$B/reads.csv" --epigenomes "$B/epigenomes.csv" \
       --lifted "$B/lifted.csv" || fail=1 ;;
  winnowmap)
    nf winnowmap --assemblies "$B/assemblies.csv" --aligner winnowmap --genotype false || fail=1 ;;
  *) echo "unknown cell: $cell (see ./run_matrix.sh -l)" >&2; fail=1 ;;
  esac
done
echo "matrix done, fail=$fail"
exit $fail
