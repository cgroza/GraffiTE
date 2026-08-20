#!/usr/bin/env bash
# Pre-flight for the HERV-K v2 discovery test. Checks inputs and tools before
# committing a cluster allocation to a run that would fail an hour in.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
if [[ -f INPUTS.env ]]; then
  # shellcheck disable=SC1091
  source ./INPUTS.env
else
  echo "INPUTS.env not found — run ./bootstrap.sh first" >&2; exit 1
fi
fail=0
ok(){   printf '  [ ok ] %s\n' "$1"; }
bad(){  printf '  [FAIL] %s\n' "$1"; fail=1; }
warn(){ printf '  [warn] %s\n' "$1"; }

echo "== inputs =="
for var in PAV_VCF REFERENCE TE_LIBRARY; do
  val="${!var:-}"
  if   [[ -z "$val"    ]]; then bad "$var is not set"
  elif [[ ! -f "$val"  ]]; then bad "$var=$val does not exist"
  else ok "$var=$val"; fi
done
[[ -n "${REFERENCE:-}" && -f "${REFERENCE:-}.fai" ]] \
  && ok "reference index ${REFERENCE}.fai" \
  || warn "no ${REFERENCE:-<REFERENCE>}.fai — the run will build one (needs a writable directory)"

if [[ -n "${TE_LIBRARY:-}" && -f "${TE_LIBRARY:-}" ]]; then
  for fam in LTR5_Hs HERVK-int; do
    grep -q ">$fam" "$TE_LIBRARY" \
      && ok "TE library contains $fam" \
      || bad "TE library has no >$fam — reference masking cannot call HML-2 states"
  done
fi

echo "== container =="
PROFILE="${PROFILE:-cluster}"
if [[ -n "${GRAFFITE_SIF:-}" ]]; then
  [[ -f "$GRAFFITE_SIF" ]] && ok "container $GRAFFITE_SIF" || bad "GRAFFITE_SIF=$GRAFFITE_SIF missing"
  RUNNER=""
  command -v apptainer   >/dev/null && RUNNER=apptainer
  [[ -z "$RUNNER" ]] && command -v singularity >/dev/null && RUNNER=singularity
  if [[ -n "$RUNNER" ]]; then
    ok "$RUNNER available"
    for t in RepeatMasker samtools bcftools python3; do
      "$RUNNER" exec "$GRAFFITE_SIF" which "$t" >/dev/null 2>&1 \
        && ok "container has $t" || bad "container is missing $t"
    done
  else
    bad "neither apptainer nor singularity on PATH"
  fi
elif [[ "$PROFILE" == "cluster" || "$PROFILE" == "aws" ]]; then
  # Tools live in library://cgroza/collection/graffite:latest, which Nextflow
  # pulls on first use. Nothing to check on the host.
  warn "GRAFFITE_SIF not set — Nextflow will pull the container for -profile $PROFILE."
  warn "Set GRAFFITE_SIF to a local .sif to check its contents here instead."
  command -v apptainer >/dev/null || command -v singularity >/dev/null \
    && ok "container runtime present" || bad "no apptainer/singularity on PATH"
else
  warn "-profile $PROFILE with no container — tools must be on PATH"
  for t in RepeatMasker samtools bcftools python3; do
    command -v "$t" >/dev/null && ok "$t" || bad "$t not on PATH"
  done
fi

echo "== data visibility =="
# nextflow.config sets singularity.runOptions = "--contain --bind $(pwd):/tmp",
# so inputs outside the launch directory need autoMounts to bind them. It is
# on by default, but a path on a filesystem the node cannot see fails late.
for var in PAV_VCF REFERENCE TE_LIBRARY; do
  val="${!var:-}"
  [[ -z "$val" ]] && continue
  case "$val" in
    /*) ok "$var is an absolute path" ;;
    *)  warn "$var=$val is relative — use an absolute path so the container resolves it" ;;
  esac
done

echo "== pipeline revision =="
PROJECT="${PROJECT:-cgroza/GraffiTE}"
REVISION="${REVISION:-v1.1dev-hervk-v2}"
GT_DIR="${NXF_ASSETS:-$HOME/.nextflow/assets}/$PROJECT"
if [[ ! -d "$GT_DIR" ]]; then
  bad "$GT_DIR not cached — run ./bootstrap.sh (or: nextflow pull $PROJECT -r $REVISION)"
else
  BR="$(cd "$GT_DIR" && git rev-parse --abbrev-ref HEAD 2>/dev/null || echo '?')"
  [[ "$BR" == "$REVISION" ]] && ok "cached revision $BR ($(cd "$GT_DIR" && git rev-parse --short HEAD))" \
    || bad "cached revision is $BR, expected $REVISION — rerun ./bootstrap.sh"
  for f in main.nf bin/hervk_arch.py bin/hervk_ref_state.py bin/hervk_classify.py bin/hervk_reconcile.py; do
    [[ -f "$GT_DIR/$f" ]] && ok "$f" || bad "$GT_DIR/$f missing"
  done
  if [[ -f "$GT_DIR/bin/hervk_arch.py" ]]; then
    ( cd "$GT_DIR" && python3 bin/hervk_classify.py --selftest \
        test/hervk/hervk_arch_fixture.tsv \
        test/hervk/hervk_refstate_fixture.tsv \
        test/hervk/hervk_expect.tsv >/dev/null 2>&1 ) \
      && ok "classifier selftest passes" \
      || bad "classifier selftest FAILED — do not run the pipeline, report this"
  fi
fi

echo "== nextflow =="
command -v nextflow >/dev/null && ok "nextflow $(nextflow -v 2>/dev/null)" \
  || bad "nextflow not on PATH (module load nextflow?)"

echo
[[ $fail -eq 0 ]] && echo "pre-flight PASSED" || echo "pre-flight FAILED — fix the above before running"
exit $fail
