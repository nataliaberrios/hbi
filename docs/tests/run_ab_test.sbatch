#!/bin/bash
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --constraint=CPU_GEN:RME
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --job-name=hbi_ab
#SBATCH --output=/scratch/users/nberrios/3dhbi/ab_test/ab_%j.log

# A/B test for the limitsigma normal-stress ratchet.
#
# Builds two binaries from the same tree -- one with the fix, one with the
# pre-fix main_LH.f90 restored -- runs the identical input deck through both,
# and leaves the output in
#   $AB/buggy/output/   and   $AB/fixed/output/
# for docs/analyse_ab_test.py to compare.

set -euo pipefail

ml purge
ml load ifort
ml load impi

REPO=/home/users/nberrios/3dhbi/hbi_git
AB=/scratch/users/nberrios/3dhbi/ab_test
FIX_COMMIT=${FIX_COMMIT:-HEAD}          # commit holding the fix
BASE_COMMIT=${BASE_COMMIT:-HEAD~1}      # commit holding the pre-fix main_LH.f90

mkdir -p "$AB"

for VARIANT in buggy fixed; do
  SRC="$AB/$VARIANT"
  echo "=============== building $VARIANT in $SRC ==============="
  rm -rf "$SRC"
  mkdir -p "$SRC"

  # same tree for both, then restore the pre-fix main_LH.f90 for the buggy one
  git -C "$REPO" archive "$FIX_COMMIT" | tar -x -C "$SRC"
  if [ "$VARIANT" = "buggy" ]; then
    git -C "$REPO" show "$BASE_COMMIT:main_LH.f90" > "$SRC/main_LH.f90"
  fi
  # the Sherlock/Intel compiler selection is an uncommitted working-tree edit,
  # so take the Makefile from the live tree rather than from the commit
  cp "$REPO/Makefile" "$SRC/Makefile"

  cd "$SRC"
  make rmod >/dev/null 2>&1 || true
  make 2>&1 | tail -5
  test -x ./lhbiem || { echo "BUILD FAILED for $VARIANT"; exit 1; }

  mkdir -p output
  cp "$REPO/docs/tests/ab_test.in" .
  cp "$REPO/docs/tests/inj_ab_test.txt" .
done

# confirm the two binaries really differ only in the intended way
echo "=============== source diff (buggy -> fixed) ==============="
diff "$AB/buggy/main_LH.f90" "$AB/fixed/main_LH.f90" || true

for VARIANT in buggy fixed; do
  echo "=============== running $VARIANT ==============="
  cd "$AB/$VARIANT"
  time ./lhbiem ab_test.in
  echo "--- output written ---"
  ls -la output/ | head
done

echo "ALL DONE"
