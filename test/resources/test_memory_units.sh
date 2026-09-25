#!/usr/bin/env bash
# Every *_memory and *_time param must carry a unit, or be null.
#
# Nextflow reads a bare number in a memory directive as BYTES. `--graph_align_memory 400`
# asks for 400 B, not 400 GB, and nothing warns: the requirement is trivially
# satisfiable, so the scheduler starts every task at once and the host runs out of
# memory. That is what issue #100 looks like from the outside, where the reporter
# raised the number twice and the run died the same way each time.
#
# The pipeline cannot stop a user writing `--graph_align_memory 400` on the command
# line. It can stop the defaults drifting into the same shape, which is what this
# checks. If a default ever loses its unit, every run inherits 400-bytes behaviour
# with no flag passed at all.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
CFG=nextflow.config
[[ -f "$CFG" ]] || { echo "  [FAIL] $CFG not found"; exit 1; }

fail=0
chk(){ if [[ "$2" == "$3" ]]; then echo "  [ ok ] $1"; else
       echo "  [FAIL] $1: got '$2', want '$3'"; fail=1; fi; }

# Defaults: null, or a quoted string ending in a memory/time unit.
bad_mem=""
while IFS= read -r line; do
  name=$(sed -E 's/^[[:space:]]*([A-Za-z_][A-Za-z0-9_]*)[[:space:]]*=.*/\1/' <<<"$line")
  val=$(sed -E 's/^[^=]*=[[:space:]]*//; s#[[:space:]]*//.*$##; s/[[:space:]]*$//' <<<"$line")
  [[ "$val" == "null" ]] && continue
  # a unit-carrying string: '40G', "10G", '2.5 GB'
  if ! grep -qiE "^['\"][0-9]+(\.[0-9]+)?[[:space:]]*(B|KB|MB|GB|TB|K|M|G|T)['\"]$" <<<"$val"; then
    bad_mem+="${name}=${val} "
  fi
done < <(grep -E "^[[:space:]]*[A-Za-z_][A-Za-z0-9_]*_memory[[:space:]]*=" "$CFG")
chk "every *_memory default is null or carries a unit" "${bad_mem:-none}" "none"

bad_time=""
while IFS= read -r line; do
  name=$(sed -E 's/^[[:space:]]*([A-Za-z_][A-Za-z0-9_]*)[[:space:]]*=.*/\1/' <<<"$line")
  val=$(sed -E 's/^[^=]*=[[:space:]]*//; s#[[:space:]]*//.*$##; s/[[:space:]]*$//' <<<"$line")
  [[ "$val" == "null" ]] && continue
  if ! grep -qiE "^['\"][0-9]+(\.[0-9]+)?[[:space:]]*(s|m|h|d|sec|min|hour|day)s?['\"]$" <<<"$val"; then
    bad_time+="${name}=${val} "
  fi
done < <(grep -E "^[[:space:]]*[A-Za-z_][A-Za-z0-9_]*_time[[:space:]]*=" "$CFG")
chk "every *_time default is null or carries a unit" "${bad_time:-none}" "none"

# Each params.*_memory referenced in a withName block must be declared in params.
missing=""
for p in $(grep -oE "params\.[A-Za-z_][A-Za-z0-9_]*_memory" "$CFG" | sort -u | sed 's/params\.//'); do
  grep -qE "^[[:space:]]*${p}[[:space:]]*=" "$CFG" || missing+="$p "
done
chk "every params.*_memory used by a process is declared" "${missing:-none}" "none"

# Nextflow's own reading of a bare number, when we can run it. A control run
# goes first: if Nextflow here will not even take `--m 1.GB`, this harness
# cannot test anything and says so rather than passing quietly. Only once the
# control works does the bare number become an assertion.
if command -v nextflow >/dev/null 2>&1; then
  tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
  cat > "$tmp/m.nf" <<'EOF'
params.m = null
process show { memory = params.m
  "echo \"MEM=${task.memory}\"" }
workflow { show() }
EOF
  runmem(){ rm -rf "$tmp/work" "$tmp/.nextflow"*;
            ( cd "$tmp" && NXF_ANSI_LOG=false timeout 240 nextflow run m.nf --m "$1" >/dev/null 2>&1 );
            find "$tmp/work" -name .command.out -exec cat {} \; 2>/dev/null | grep -o "MEM=.*" | head -1; }
  control=$(runmem 1.GB)
  if [[ "$control" == "MEM=1 GB" ]]; then
    chk "a bare number is read as bytes, not GB" "$(runmem 400)" "MEM=400 B"
  else
    echo "  [SKIP] this Nextflow did not honour --m 1.GB (got '${control:-nothing}'),"
    echo "         so the bare-number check cannot run here. Static checks above still apply."
  fi
else
  echo "  [SKIP] nextflow not on PATH; the bare-number check cannot run here"
fi

exit $fail
