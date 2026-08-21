#!/bin/bash
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --job-name=hbibuild
#SBATCH --output=/scratch/users/nberrios/3dhbi/build_%j.log

# Rebuild ./lhbiem in the repo from the current working tree.
#
# The binary sitting in hbi_git predates the limitsigma fix (built 2026-06-25,
# main_LH.f90 last touched 2026-07-31), so running it silently reproduces the
# bug. This job replaces it and asserts the fix is actually in the built source.
#
# Use `make clean`, NOT `make rmod`, on an already-built tree. rmod is
# `rm -f m_*.o *.mod`: it removes every .mod but only the m_*.o objects, so
# okada.o and TDstressFS.o survive, make considers them current, and okada.mod /
# dtriangular.mod are never regenerated -- m_HACApK_calc_entry_ij.f90 then dies on
# "Error in opening the compiled module file [DTRIANGULAR]". rmod is only safe on
# a fresh checkout where there are no stale objects to begin with.

set -euo pipefail

ml purge
ml load ifort
ml load impi

REPO=/home/users/nberrios/3dhbi/hbi_git
cd "$REPO"

if ! grep -q "sigmat(i)=sigmat(i)+(sigma(i)-sigmac(i))" main_LH.f90; then
  echo "FIX NOT PRESENT in main_LH.f90 -- refusing to build"; exit 1
fi

make clean
make 2>&1 | tail -20
test -x ./lhbiem || { echo "BUILD FAILED"; exit 1; }

echo "--- built ---"
ls -l --time-style=+'%Y-%m-%d %H:%M' ./lhbiem

printf '| %s | hbibuild | COMPLETED | - | `%s` — rebuilt `hbi_git/lhbiem` from the working tree |\n' \
  "$SLURM_JOB_ID" "/scratch/users/nberrios/3dhbi/build_${SLURM_JOB_ID}.log" \
  >> /home/users/nberrios/3dhbi/JOB_STATUS.md

echo "BUILD DONE"
