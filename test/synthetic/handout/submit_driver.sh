#!/usr/bin/env bash
#SBATCH --job-name=graffite-synth
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=8:00:00
#SBATCH --output=synth_driver_%j.out
#
# The Nextflow driver runs INSIDE one allocation with the LOCAL executor.
#
# Not -profile cluster: on Puma that submits one slurm job per task and charges
# roughly 85 minutes of queue latency each. A 12 h job died that way (23613713)
# having finished a fraction of the work. One allocation, local executor, and
# local.config sized to the allocation is what has actually worked here.
#
# Adjust --cpus-per-task/--mem/--time together with CPUS/MEM_GB in INPUTS.env;
# local.config reads the latter and they must agree or Nextflow oversubscribes
# the allocation until the scheduler kills it.
#
#   sbatch submit_driver.sh            # the cells named in RUNS
#   sbatch submit_driver.sh spine      # one cell
set -uo pipefail
cd "$SLURM_SUBMIT_DIR" 2>/dev/null || cd "$(dirname "${BASH_SOURCE[0]}")"
module load nextflow 2>/dev/null || true
./run_matrix.sh "$@"
rc=$?
./bundle_results.sh || true
exit $rc
