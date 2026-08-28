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

echo "== entry point =="
if [[ -n "${RM_DIR:-}" ]]; then
  if [[ ! -d "$RM_DIR" ]]; then
    bad "RM_DIR=$RM_DIR is not a directory"
  else
    ok "RM_DIR=$RM_DIR (RepeatMasker will be skipped)"
    n_sh=$(find "$RM_DIR" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l | tr -d ' ')
    n_vcf=$(find "$RM_DIR" -mindepth 2 -maxdepth 2 -name genotypes_repmasked_filtered.vcf 2>/dev/null | wc -l | tr -d ' ')
    n_rm=$(find "$RM_DIR" -mindepth 2 -maxdepth 2 -type d -name repeatmasker_dir 2>/dev/null | wc -l | tr -d ' ')
    n_out=$(find "$RM_DIR" -mindepth 3 -maxdepth 3 -name indels.fa.out 2>/dev/null | wc -l | tr -d ' ')
    ok "$n_sh shard directories"
    [[ "$n_vcf" -eq "$n_sh" && "$n_sh" -gt 0 ]] \
      && ok "$n_vcf genotypes_repmasked_filtered.vcf" \
      || bad "$n_vcf of $n_sh shards have genotypes_repmasked_filtered.vcf — is this a published 2_Repeat_Filtering output rather than a work/ directory?"
    [[ "$n_rm" -eq "$n_sh" ]] && ok "$n_rm repeatmasker_dir" \
      || bad "$n_rm of $n_sh shards have repeatmasker_dir"
    [[ "$n_out" -eq "$n_sh" ]] && ok "$n_out indels.fa.out (the HERV-K architecture source)" \
      || bad "$n_out of $n_sh shards have repeatmasker_dir/indels.fa.out — without these the classifier cannot read LTR architecture"
    hk=$(grep -lE "[[:space:]](HERVK[-_]?(int(ernal)?)?|LTR5_Hs|LTR5A|LTR5B)[[:space:]]" "$RM_DIR"/*/repeatmasker_dir/indels.fa.out 2>/dev/null | wc -l | tr -d ' ')
    [[ "$hk" -gt 0 ]] && ok "HML-2 hits present in $hk shard(s)" \
      || bad "no HML-2 hits (HERVK*/LTR5*) anywhere in the RepeatMasker output"
  fi
elif [[ -n "${PAV_VCF:-}" ]]; then
  [[ -f "$PAV_VCF" ]] && ok "PAV_VCF=$PAV_VCF (full re-mask, slow path)" \
    || bad "PAV_VCF=$PAV_VCF does not exist"
else
  bad "neither RM_DIR nor PAV_VCF is set in INPUTS.env"
fi

echo "== stage E =="
if [[ -n "${GENOTYPED_VCF:-}" ]]; then
  if [[ ! -f "$GENOTYPED_VCF" ]]; then
    bad "GENOTYPED_VCF=$GENOTYPED_VCF does not exist"
  else
    ok "GENOTYPED_VCF=$GENOTYPED_VCF"
    n=$(bcftools query -l "$GENOTYPED_VCF" 2>/dev/null | wc -l | tr -d ' ')
    [[ "${n:-0}" -gt 0 ]] && ok "$n samples in the genotyped VCF" \
      || warn "could not read samples from it (bcftools missing on this host?)"
  fi
else
  warn "GENOTYPED_VCF not set — stage E will be skipped, only discovery runs"
fi

echo "== inputs =="
for var in REFERENCE TE_LIBRARY; do
  val="${!var:-}"
  if   [[ -z "$val"    ]]; then bad "$var is not set"
  elif [[ ! -f "$val"  ]]; then bad "$var=$val does not exist"
  else ok "$var=$val"; fi
done
[[ -n "${REFERENCE:-}" && -f "${REFERENCE:-}.fai" ]] \
  && ok "reference index ${REFERENCE}.fai" \
  || warn "no ${REFERENCE:-<REFERENCE>}.fai — the run will build one (needs a writable directory)"

if [[ -n "${TE_LIBRARY:-}" && -f "${TE_LIBRARY:-}" ]]; then
  # The HML-2 internal region is named HERVK by Dfam and HERVK-int by other
  # sets; both are accepted, so no library needs renaming. HERVK9/11/14-int are
  # different lineages and do not count.
  int_name=$(grep -oE "^>(HERVK[-_]?(int(ernal)?)?)\b" "$TE_LIBRARY" 2>/dev/null \
             | sed 's/^>//' | sort -u | head -1)
  if [[ -n "$int_name" ]]; then
    ok "TE library internal region: $int_name"
  else
    bad "TE library has no HML-2 internal region (looked for HERVK, HERVK-int, HERVK_int)"
  fi
  ltr_found=$(grep -cE "^>(LTR5_Hs|LTR5A|LTR5B|LTR5)\b" "$TE_LIBRARY" 2>/dev/null || echo 0)
  [[ "${ltr_found:-0}" -gt 0 ]] \
    && ok "TE library HML-2 LTRs: $ltr_found (LTR5*)" \
    || bad "TE library has no LTR5_Hs/LTR5A/LTR5B — reference masking cannot call HML-2 states"
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
elif [[ "$PROFILE" == "cluster" || "$PROFILE" == "aws" || "$PROFILE" == "standard" ]]; then
  # Tools live in library://cgroza/collection/graffite:latest, which Nextflow
  # pulls on first use. Nothing to check on the host.
  # Every profile in nextflow.config sets process.container, and
  # singularity.enabled is set outside the profiles block -- so `standard`
  # means "local executor", not "no container". The tools are in the image,
  # never on the host PATH.
  warn "GRAFFITE_SIF not set — Nextflow will pull the container for -profile $PROFILE."
  warn "Set GRAFFITE_SIF to a local .sif to check its contents here instead."
  if command -v apptainer >/dev/null || command -v singularity >/dev/null; then
    ok "container runtime present"
  else
    # A login node often has no runtime while the compute nodes do.
    warn "no apptainer/singularity on this host — fine if you submit to nodes that have it"
  fi
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
for var in RM_DIR PAV_VCF REFERENCE TE_LIBRARY GENOTYPED_VCF; do
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
