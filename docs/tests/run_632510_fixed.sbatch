#!/bin/bash
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --constraint=CPU_GEN:RME
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --job-name=fix632510
#SBATCH --output=/scratch/users/nberrios/3dhbi/fix_632510/run_%j.log

# Rerun job 632510's own deck with the post-fix binary.
#
# The A/B test (run_ab_test.sbatch) proves the normal-stress ratchet is gone,
# but its deck is 2 days at 151x151 and never grows the injector slip cusp, so
# it cannot say whether the ratchet is what produced the cusp in 632510. This
# job is the end-to-end test of that: same deck, same inputs, only the binary
# differs, so cumulative slip through the injector can be compared directly
# against the recorded pre-fix run already on disk at
#   /scratch/users/nberrios/3dhbi/output/632510/
#
# The original took 15h08m on these same resources (job 36097851, 1 MPI rank x
# 16 OpenMP threads), hence the 24h request.
#
# The deck is copied byte-for-byte, filenumber 632510 included, so the ONLY
# difference from the recorded run is the executable. Output lands in this job's
# own directory, so the recorded run is not touched.

set -euo pipefail

ml purge
ml load ifort
ml load impi

REPO=/home/users/nberrios/3dhbi/hbi_git
ORIG=/scratch/users/nberrios/3dhbi/runs/36097851     # the original 632510 run
RUN=/scratch/users/nberrios/3dhbi/fix_632510
FIX_COMMIT=${FIX_COMMIT:-94da536}                    # "Stop the limitsigma clamp..."

mkdir -p "$RUN"
SRC="$RUN/src"
rm -rf "$SRC"; mkdir -p "$SRC"

echo "=============== building the fixed binary from $FIX_COMMIT ==============="
( cd "$REPO" && git archive "$FIX_COMMIT" ) | tar -x -C "$SRC"
# the Sherlock/Intel compiler selection is an uncommitted working-tree edit
cp "$REPO/Makefile" "$SRC/Makefile"
cd "$SRC"
make rmod >/dev/null 2>&1 || true
make 2>&1 | tail -5
test -x ./lhbiem || { echo "BUILD FAILED"; exit 1; }

# confirm the built source really carries the fix
if grep -q "sigmat(i)=sigmat(i)+(sigma(i)-sigmac(i))" main_LH.f90; then
  echo "fix present in main_LH.f90"
else
  echo "FIX NOT PRESENT in the built source -- wrong commit?"; exit 1
fi

echo "=============== staging 632510's deck ==============="
cd "$RUN"
mkdir -p output
cp "$SRC/lhbiem" .
for f in res632510.in perm_632510_ds5.txt june_clean.txt; do
  cp "$ORIG/$f" .
done
# prove the deck is unmodified relative to the recorded run
diff "$ORIG/res632510.in" ./res632510.in && echo "deck identical to the recorded run"

echo "=============== running ==============="
export OMP_NUM_THREADS=16
time ./lhbiem res632510.in

echo "--- output ---"
ls -la output/
echo "ALL DONE"
