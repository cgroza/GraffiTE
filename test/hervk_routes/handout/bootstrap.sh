#!/usr/bin/env bash
# Fetch the GraffiTE revision and drop this handout into the current directory.
# Run this first, from the directory you want to work in.
set -euo pipefail
PROJECT="${PROJECT:-cgroza/GraffiTE}"
REVISION="${REVISION:-v1.1dev-hervk-v2}"

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
