#!/usr/bin/env bash
#SBATCH --job-name=graffite-127
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=issue127_%j.out
#
# Parts A and B in one allocation, then the tally and the bundle. Keep
# --cpus-per-task and --mem in step with CPUS and MEM_GB in INPUTS.env.
#
#   sbatch submit.sh          # both parts
#   sbatch submit.sh A        # Part A only
#   sbatch submit.sh B        # Part B only
set -uo pipefail
cd "$SLURM_SUBMIT_DIR" 2>/dev/null || cd "$(dirname "${BASH_SOURCE[0]}")"
module load nextflow 2>/dev/null || true
what="${1:-AB}"
[[ "$what" == *A* ]] && ./part_a.sh
[[ "$what" == *B* ]] && ./part_b.sh
./summarize.sh
./bundle_results.sh || true
