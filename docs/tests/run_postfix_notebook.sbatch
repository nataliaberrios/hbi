#!/bin/bash
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --job-name=postfixfigs
#SBATCH --output=/scratch/users/nberrios/3dhbi/fix_632510/postfix_notebook_%j.log

# The run_batch_figs.sbatch figure set, but against the POST-FIX rerun of 632510.
#
# The stock driver (batch_run_notebook.py) hardcodes
# DATA_DIR = /scratch/users/nberrios/3dhbi/output/{JOBID}. The post-fix rerun keeps
# filenumber 632510 while living elsewhere, so the stock driver would re-read the
# pre-fix data and clobber batch_figures/632510/. batch_run_postfix.py takes the
# data dir and an output label explicitly instead.
#
# 96 G rather than the usual 64 G: the notebook load_simulation()s six ~4 GB fields
# and several cells then take full copies (pf_grid, slip_grid_cm, slip_anim).

set -euo pipefail

cd /home/users/nberrios/3dhbi/hbi_analysis/notebooks

/home/users/nberrios/miniconda3/bin/python batch_run_postfix.py \
    632510 \
    /scratch/users/nberrios/3dhbi/fix_632510/output \
    632510_postfix

printf '| %s | post-fix notebook figures | notebooks/batch_figures/632510_postfix/ | - | `%s` |\n' \
  "$SLURM_JOB_ID" "/scratch/users/nberrios/3dhbi/fix_632510/postfix_notebook_${SLURM_JOB_ID}.log" \
  >> /home/users/nberrios/3dhbi/JOB_STATUS.md

echo "POSTFIX NOTEBOOK FIGURES DONE"
