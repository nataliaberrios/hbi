#!/bin/bash
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --job-name=analyse_ab
#SBATCH --output=/scratch/users/nberrios/3dhbi/ab_test/analyse_%j.log

# Run the T1-T4 checks on the output of run_ab_test.sbatch.

set -euo pipefail

ml purge
ml load math devel
ml load py-numpy/1.24.2_py39
ml load viz py-matplotlib/3.7.1_py39

python3 /home/users/nberrios/3dhbi/hbi_git/docs/analyse_ab_test.py
echo "ANALYSIS DONE"
