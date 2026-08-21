#!/bin/bash
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --job-name=triage
#SBATCH --output=/scratch/users/nberrios/3dhbi/fix_632510/triage_%j.log

# Which recorded runs did the limitsigma ratchet actually corrupt?
# Part B streams every monitor*.dat under /scratch/users/nberrios/3dhbi/output
# (~333 run dirs, several MB each), hence the 2 h ask.

set -euo pipefail

ml purge
ml load math devel
ml load py-numpy/1.24.2_py39

python3 /home/users/nberrios/3dhbi/hbi_git/docs/triage_ratchet_impact.py

printf '| %s | ratchet triage across recorded runs | - | - | `%s` |\n' \
  "$SLURM_JOB_ID" "/scratch/users/nberrios/3dhbi/fix_632510/triage_${SLURM_JOB_ID}.log" \
  >> /home/users/nberrios/3dhbi/JOB_STATUS.md

echo "TRIAGE DONE"
