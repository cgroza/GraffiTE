#!/usr/bin/env bash
# Check the environment, then build the synthetic inputs. Stop here if it fails.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
[[ -f INPUTS.env ]] || { echo "INPUTS.env not found -- run ./bootstrap.sh first" >&2; exit 1; }
# shellcheck disable=SC1091
source ./INPUTS.env

nbad=0
ok()   { echo "  [ ok ] $*"; }
warn() { echo "  [warn] $*"; }
bad()  { echo "  [FAIL] $*"; nbad=$((nbad+1)); }

echo "== paths =="
[[ -n "${WORKDIR:-}" ]] || bad "WORKDIR is empty in INPUTS.env"
if [[ -n "${WORKDIR:-}" ]]; then
  [[ "$WORKDIR" = /* ]] || bad "WORKDIR must be absolute (--contain)"
  mkdir -p "$WORKDIR" 2>/dev/null && ok "WORKDIR $WORKDIR" || bad "cannot create $WORKDIR"
  case "$WORKDIR" in
    /home/*|"$HOME"/*) warn "WORKDIR is under \$HOME. work/ will outgrow a home quota; /xdisk is the place for it.";;
  esac
  avail=$(df -Pk "$WORKDIR" 2>/dev/null | awk 'NR==2{print int($4/1048576)}')
  [[ -n "$avail" && "$avail" -ge 20 ]] && ok "${avail} GB free on WORKDIR" \
    || warn "only ${avail:-?} GB free on WORKDIR; the spine's work/ wants room"
fi
if [[ -n "${CONTAINER_TMP:-}" ]]; then
  mkdir -p "$CONTAINER_TMP" 2>/dev/null && ok "CONTAINER_TMP $CONTAINER_TMP" \
    || bad "cannot create CONTAINER_TMP=$CONTAINER_TMP"
else
  warn "CONTAINER_TMP empty: /tmp inside the container is the launch directory. If it fills, the run SKIPS tsd_search and tsd_report and still reports success."
fi

echo "== tools on the driver =="
command -v nextflow >/dev/null && ok "nextflow $(nextflow -v 2>/dev/null | head -1)" \
  || bad "nextflow not on PATH (module load it)"
command -v python3 >/dev/null && ok "python3 $(python3 -V 2>&1)" || bad "python3 not on PATH"
RUNNER="$(command -v apptainer || command -v singularity || true)"
[[ -n "$RUNNER" ]] && ok "container runtime $RUNNER" \
  || bad "no apptainer/singularity on PATH. The login node usually has none -- run this inside the allocation."

echo "== pipeline =="
PROJECT="${PROJECT:-cgroza/GraffiTE}"
GT_DIR="${GRAFFITE_REPO:-${NXF_ASSETS:-$HOME/.nextflow/assets}/$PROJECT}"
if [[ -d "$GT_DIR" ]]; then
  ok "checkout $GT_DIR ($(cd "$GT_DIR" && git rev-parse --short HEAD 2>/dev/null || echo '?'))"
  for f in main.nf module/main.nf bin/hervk_arch.py bin/genotyping_audit.py test/human_test_set.tar.gz; do
    [[ -e "$GT_DIR/$f" ]] || bad "missing $f in the checkout"
  done
  [[ -f "$GT_DIR/panmethyl/module/main.nf" ]] && ok "panmethyl submodule present" \
    || bad "panmethyl/module/main.nf missing -- the include in main.nf will fail"
  # -latest refuses to run against a dirty asset repo, and the failure is a
  # three-second exit that reads like a config error.
  if [[ -n "$(cd "$GT_DIR" && git status --porcelain 2>/dev/null)" ]]; then
    warn "the cached checkout has local modifications; nextflow will refuse -latest"
  fi
else
  bad "no checkout at $GT_DIR -- run ./bootstrap.sh"
fi

echo "== profile =="
case "${PROFILE:-}" in
  standard) ok "PROFILE=standard (local executor), the right choice inside an sbatch allocation" ;;
  cluster)  warn "PROFILE=cluster submits one slurm job per task. On Puma that has cost ~85 min of queue latency PER TASK and timed out a 12 h run. Use standard inside sbatch." ;;
  *)        bad "PROFILE='${PROFILE:-}' is not standard or cluster" ;;
esac

echo "== build the inputs =="
DFAM_TAR="$GT_DIR/test/human_test_set.tar.gz"
if [[ -f "$DFAM_TAR" && -n "${WORKDIR:-}" ]]; then
  DF="$WORKDIR/dfam/human_DFAM3.6.fasta"
  if [[ ! -f "$DF" ]]; then
    mkdir -p "$WORKDIR/dfam"
    tar xzf "$DFAM_TAR" -C "$WORKDIR/dfam" --strip-components=1 GraffiTE_testset/human_DFAM3.6.fasta \
      && ok "extracted the Dfam library" || bad "could not extract human_DFAM3.6.fasta"
  else
    ok "Dfam library already at $DF"
  fi
  if [[ -f "$DF" ]]; then
    python3 ./build_synthetic.py --dfam "$DF" --repo "$GT_DIR" --out "$WORKDIR/build" \
      --seed "${SEED:-20260918}" --short-depth "${SHORT_DEPTH:-30}" --long-depth "${LONG_DEPTH:-15}" \
      && ok "inputs built in $WORKDIR/build" || bad "build_synthetic.py failed"
  fi
else
  bad "no $DFAM_TAR"
fi

echo
if [[ $nbad -gt 0 ]]; then
  echo "preflight: $nbad problem(s). Fix them before running the matrix."; exit 1
fi
echo "preflight OK. Next: sbatch submit_driver.sh   (or ./run_matrix.sh spine if you are already inside an allocation)"
