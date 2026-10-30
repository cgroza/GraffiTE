# Sourced by the other scripts: INPUTS.env, its defaults, and the container.
# shellcheck shell=bash
HANDOUT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ -f "$HANDOUT/INPUTS.env" ]] || { echo "INPUTS.env not found -- run ./bootstrap.sh first" >&2; exit 1; }
# REPS_A and REPS_B set on the command line win over INPUTS.env, so a smoke
# test (REPS_A=1 ./part_a.sh) needs no edit.
_cli_reps_a="${REPS_A:-}"; _cli_reps_b="${REPS_B:-}"
# shellcheck disable=SC1091
source "$HANDOUT/INPUTS.env"
[[ -n "$_cli_reps_a" ]] && REPS_A="$_cli_reps_a"
[[ -n "$_cli_reps_b" ]] && REPS_B="$_cli_reps_b"

: "${WORKDIR:?WORKDIR is empty in INPUTS.env}"
: "${SYNTH_DIR:?SYNTH_DIR is empty in INPUTS.env}"
REFERENCE="${REFERENCE:-$SYNTH_DIR/build/ref/synth.fa}"
GRAFFITE_VCF="${GRAFFITE_VCF:-$SYNTH_DIR/runs/spine/3_TSD_search/pangenome.vcf}"
READS_CSV="${READS_CSV:-$SYNTH_DIR/build/reads_long.csv}"
READS_FASTQ="${READS_FASTQ:-$SYNTH_DIR/build/reads/L1.long.fq.gz}"
SCRATCH_DIRS="${SCRATCH_DIRS:-${SLURM_TMPDIR:-/tmp} $WORKDIR/partA_src}"
REPS_A="${REPS_A:-10}"
REPS_B="${REPS_B:-3}"
NF_SCRATCH="${NF_SCRATCH:-\$SLURM_TMPDIR}"
CPUS="${CPUS:-8}"
MEM_GB="${MEM_GB:-32}"
REVISION="${REVISION:-test/giraffe-mtime-127}"
PROJECT="${PROJECT:-cgroza/GraffiTE}"
IMAGE_TAG="1.1.0.0"

RUNNER=""
command -v apptainer   >/dev/null && RUNNER=apptainer
[[ -z "$RUNNER" ]] && command -v singularity >/dev/null && RUNNER=singularity

# The image both parts use. Nextflow names a pulled docker://cgroza/graffite:X
# image cgroza-graffite-X.img in its cache, so pulling it there once serves
# Part A directly and Part B through NXF_SINGULARITY_CACHEDIR.
export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-$WORKDIR/singularity}"
IMG="${GRAFFITE_SIF:-$NXF_SINGULARITY_CACHEDIR/cgroza-graffite-$IMAGE_TAG.img}"

ensure_image() {
  [[ -n "$RUNNER" ]] || { echo "no apptainer/singularity on this host -- run inside the allocation" >&2; return 1; }
  [[ -f "$IMG" ]] && return 0
  mkdir -p "$(dirname "$IMG")"
  echo "pulling docker://cgroza/graffite:$IMAGE_TAG -> $IMG"
  "$RUNNER" pull "$IMG" "docker://cgroza/graffite:$IMAGE_TAG"
}

# ctr <workdir> <command>: run a command in the container, the way nextflow.config
# does (--contain, /tmp bound to the task directory), with the input and output
# trees bound at the same paths.
ctr() {
  local wd="$1"; shift
  local binds=(--bind "$wd:/tmp" --bind "$wd" --bind "$WORKDIR")
  local d
  for d in "$(dirname "$REFERENCE")" "$(dirname "$GRAFFITE_VCF")" "$(dirname "$READS_FASTQ")"; do
    binds+=(--bind "$d")
  done
  "$RUNNER" exec --contain "${binds[@]}" --pwd "$wd" "$IMG" bash -c "$*"
}

fs_type() { stat -f -c %T "$1" 2>/dev/null || echo '?'; }
