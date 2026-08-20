#!/usr/bin/env bash
# Pre-flight for the HERV-K v2 discovery test. Checks inputs and tools before
# committing a cluster allocation to a run that would fail an hour in.
set -uo pipefail
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

echo "== tools =="
for t in nextflow; do
  command -v "$t" >/dev/null && ok "$t $(nextflow -v 2>/dev/null)" || bad "$t not on PATH"
done
if [[ -n "${GRAFFITE_SIF:-}" ]]; then
  [[ -f "$GRAFFITE_SIF" ]] && ok "container $GRAFFITE_SIF" || bad "GRAFFITE_SIF=$GRAFFITE_SIF missing"
  for t in RepeatMasker samtools bcftools python3; do
    singularity exec "$GRAFFITE_SIF" which "$t" >/dev/null 2>&1 \
      && ok "container has $t" || bad "container is missing $t"
  done
else
  warn "GRAFFITE_SIF not set — assuming tools are on PATH (-profile local)"
  for t in RepeatMasker samtools bcftools python3; do
    command -v "$t" >/dev/null && ok "$t" || bad "$t not on PATH"
  done
fi

echo "== pipeline =="
GT_DIR="${GT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
for f in main.nf bin/hervk_arch.py bin/hervk_ref_state.py bin/hervk_classify.py bin/hervk_reconcile.py; do
  [[ -f "$GT_DIR/$f" ]] && ok "$f" || bad "$GT_DIR/$f missing — wrong branch? expected v1.1dev-hervk-v2"
done
if [[ -f "$GT_DIR/bin/hervk_arch.py" ]]; then
  ( cd "$GT_DIR" && python3 bin/hervk_classify.py --selftest \
      test/hervk/hervk_arch_fixture.tsv \
      test/hervk/hervk_refstate_fixture.tsv \
      test/hervk/hervk_expect.tsv >/dev/null 2>&1 ) \
    && ok "classifier selftest passes" \
    || bad "classifier selftest FAILED — do not run the pipeline, report this"
fi

echo
[[ $fail -eq 0 ]] && echo "pre-flight PASSED" || echo "pre-flight FAILED — fix the above before running"
exit $fail
