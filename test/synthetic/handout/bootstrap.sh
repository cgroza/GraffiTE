#!/usr/bin/env bash
# Fetch the GraffiTE revision and drop this handout into the current directory.
# Run this first, from the directory you want to work in.
set -euo pipefail
_env_rev=""; _env_proj=""
if [[ -f ./INPUTS.env ]]; then
  # shellcheck disable=SC1091
  _env_rev="$(source ./INPUTS.env >/dev/null 2>&1; echo "${REVISION:-}")"
  _env_proj="$(source ./INPUTS.env >/dev/null 2>&1; echo "${PROJECT:-}")"
fi
PROJECT="${PROJECT:-${_env_proj:-cgroza/GraffiTE}}"
REVISION="${REVISION:-${_env_rev:-test/synthetic-end-to-end}}"

command -v nextflow >/dev/null || {
  echo "nextflow not on PATH -- module load it first" >&2; exit 1; }

echo "pulling $PROJECT -r $REVISION ..."
nextflow pull "$PROJECT" -r "$REVISION"

ASSET="${NXF_ASSETS:-$HOME/.nextflow/assets}/$PROJECT"
SRC="$ASSET/test/synthetic/handout"
[[ -d "$SRC" ]] || { echo "handout not found at $SRC" >&2; exit 1; }

for f in "$SRC"/*; do
  b="$(basename "$f")"
  if [[ "$b" == "INPUTS.env" && -f ./INPUTS.env ]]; then
    echo "  keeping your existing INPUTS.env (new template at INPUTS.env.new)"
    cp "$f" ./INPUTS.env.new
  elif [[ -f "./$b" ]] && ! cmp -s "$f" "./$b"; then
    cp "./$b" "./$b.local.bak"
    cp "$f" ./
    echo "  $b differed -- yours saved as $b.local.bak"
  else
    cp "$f" ./
  fi
done
chmod +x ./*.sh ./*.py 2>/dev/null || true

# The Dfam library the generator reads its consensus sequences from lives in a
# tarball in the repo, not in the handout: 36 MB is too much to copy per run.
echo "  dfam source: $ASSET/test/human_test_set.tar.gz"

echo
echo "handout ready in $PWD"
echo "revision: $(cd "$ASSET" && git rev-parse --short HEAD 2>/dev/null || echo '?') on $REVISION"
echo
echo "next: edit INPUTS.env, then ./preflight.sh"
