#!/bin/bash
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --job-name=frontdp
#SBATCH --output=/scratch/users/nberrios/3dhbi/fix_632510/frontdp_%j.log

# Contour the 632507 pressure field at each model's failure margin, to test
# whether the slip-front gap vs Taiyi is a prestress difference rather than an
# elasticity-solver difference.  Memmaps pf632507.dat frame by frame, so 16 G is
# plenty for the 728 MB field.

set -euo pipefail

/home/users/nberrios/miniconda3/bin/python \
    /home/users/nberrios/3dhbi/hbi_git/docs/front_vs_dp_contour.py

echo "FRONT-VS-DP SBATCH DONE"
