#!/bin/bash
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=48G
#SBATCH --job-name=prepostfigs
#SBATCH --output=/scratch/users/nberrios/3dhbi/fix_632510/prepost_figs_%j.log

# Pre-fix vs post-fix overlay figures for 632510.
# Fields are memmapped, but figures 5 and 7 stream every frame of sigma and slip
# for both runs (~16 GB of reads), hence the 2 h / 48 G ask rather than the 30 min
# the cusp comparison needed.

set -euo pipefail

ml purge
ml load math devel
ml load py-numpy/1.24.2_py39
ml load viz py-matplotlib/3.7.1_py39

python3 /home/users/nberrios/3dhbi/hbi_git/docs/compare_prepost_632510.py

printf '| %s | pre/post overlay figures | docs/figs/prepost_632510/ | - | `%s` |\n' \
  "$SLURM_JOB_ID" "/scratch/users/nberrios/3dhbi/fix_632510/prepost_figs_${SLURM_JOB_ID}.log" \
  >> /home/users/nberrios/3dhbi/JOB_STATUS.md

echo "PREPOST FIGURES DONE"
