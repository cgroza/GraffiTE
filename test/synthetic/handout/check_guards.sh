#!/usr/bin/env bash
# The launch-time guards. No container, no data: each should stop at launch
# with its own message, and none of them is tested by anything else.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
# shellcheck disable=SC1091
source ./INPUTS.env
PROJECT="${PROJECT:-cgroza/GraffiTE}"; REVISION="${REVISION:-test/synthetic-end-to-end}"
B="${WORKDIR:?}/build"
fail=0
expect() {  # expect <label> <substring> <args...>
  local label="$1" want="$2"; shift 2
  local out; out=$(nextflow run "$PROJECT" -r "$REVISION" -preview "$@" 2>&1)
  if grep -qF -- "$want" <<<"$out"; then echo "  [ ok ] $label"
  else echo "  [FAIL] $label: did not say '$want'"; echo "$out" | tail -3 | sed 's/^/         /'; fail=1; fi
}
expect "no input at all" "No input given" --reference "$B/ref/synth.fa"
expect "--graffite_vcf with --human and reconciliation" "HERV-K reconciliation needs" \
  --graffite_vcf "$B/vcf/merged.vcf.gz" --reference "$B/ref/synth.fa" --human
expect "unsupported --graph_method" "graph_method must be" \
  --graffite_vcf "$B/vcf/merged.vcf.gz" --reference "$B/ref/synth.fa" --graph_method nope
echo "guards fail=$fail"; exit $fail
