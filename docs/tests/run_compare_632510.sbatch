#!/bin/bash
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --job-name=cusp632510
#SBATCH --output=/scratch/users/nberrios/3dhbi/fix_632510/compare_%j.log

# Compare the recorded pre-fix 632510 against the post-fix rerun.
# Submit with --dependency=afterok:<rerun jobid>.

set -euo pipefail

ml purge
ml load math devel
ml load py-numpy/1.24.2_py39
ml load viz py-matplotlib/3.7.1_py39

python3 /home/users/nberrios/3dhbi/hbi_git/docs/compare_632510_cusp.py

# durable record, so the result is on disk whether or not anyone is watching
printf '| %s | cusp comparison | see compare log | - | `%s` |\n' \
  "$SLURM_JOB_ID" "/scratch/users/nberrios/3dhbi/fix_632510/compare_${SLURM_JOB_ID}.log" \
  >> /home/users/nberrios/3dhbi/JOB_STATUS.md

echo "COMPARISON DONE"
