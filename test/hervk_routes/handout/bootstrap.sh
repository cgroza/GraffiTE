#!/usr/bin/env bash
# Fetch the GraffiTE revision and drop this handout into the current directory.
# Run this first, from the directory you want to work in.
set -euo pipefail
# Read INPUTS.env when it is already here, so the revision does not have to be
# passed twice. An explicit REVISION= in the environment still wins, which is
# what a first bootstrap onto a new branch needs.
_env_rev=""; _env_proj=""
if [[ -f ./INPUTS.env ]]; then
  # shellcheck disable=SC1091
  _env_rev="$(source ./INPUTS.env >/dev/null 2>&1; echo "${REVISION:-}")"
  _env_proj="$(source ./INPUTS.env >/dev/null 2>&1; echo "${PROJECT:-}")"
fi
PROJECT="${PROJECT:-${_env_proj:-cgroza/GraffiTE}}"
REVISION="${REVISION:-${_env_rev:-fix/hervk-pair-rule-n-hits}}"

command -v nextflow >/dev/null || {
  echo "nextflow not on PATH — module load it first (e.g. 'module load nextflow')" >&2
  exit 1; }

echo "pulling $PROJECT -r $REVISION ..."
nextflow pull "$PROJECT" -r "$REVISION"

ASSET="${NXF_ASSETS:-$HOME/.nextflow/assets}/$PROJECT"
SRC="$ASSET/test/hervk_routes/handout"
[[ -d "$SRC" ]] || { echo "handout not found at $SRC" >&2; exit 1; }

# Never clobber a populated INPUTS.env on re-run.
for f in "$SRC"/*; do
  b="$(basename "$f")"
  if [[ "$b" == "INPUTS.env" && -f ./INPUTS.env ]]; then
    echo "  keeping your existing INPUTS.env (new template at INPUTS.env.new)"
    cp "$f" ./INPUTS.env.new
  elif [[ -f "./$b" ]] && ! cmp -s "$f" "./$b"; then
    # Do not silently discard a local edit. Run 1 patched preflight.sh for the
    # local-executor setup and bootstrap reverted it, so the next preflight
    # failed on "RepeatMasker not on PATH" for no visible reason.
    cp "./$b" "./$b.local.bak"
    cp "$f" ./
    echo "  $b differed — yours saved as $b.local.bak"
  else
    cp "$f" ./
  fi
done
chmod +x ./*.sh ./*.py 2>/dev/null || true

echo
echo "handout ready in $PWD"
echo "revision: $(cd "$ASSET" && git rev-parse --short HEAD 2>/dev/null || echo '?') on $REVISION"
echo
echo "next: edit INPUTS.env, then ./preflight.sh"
