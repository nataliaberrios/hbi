#!/bin/bash
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --constraint=CPU_GEN:RME
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --job-name=hbi_ab_ctl
#SBATCH --output=/scratch/users/nberrios/3dhbi/ab_test/ctl_%j.log

# Null test: run ab_control.in (sigmainit 500 MPa, so the clamp never fires)
# through both binaries built by run_ab_test.sbatch. With no clamp corrections
# to accumulate, the two must agree bit-for-bit.
#
# Run run_ab_test.sbatch first -- this reuses the binaries it built.

set -euo pipefail

ml purge
ml load ifort
ml load impi

REPO=/home/users/nberrios/3dhbi/hbi_git
AB=/scratch/users/nberrios/3dhbi/ab_test

for VARIANT in buggy fixed; do
  test -x "$AB/$VARIANT/lhbiem" || { echo "missing $AB/$VARIANT/lhbiem -- run run_ab_test.sbatch first"; exit 1; }
  DIR="$AB/ctl_$VARIANT"
  rm -rf "$DIR"; mkdir -p "$DIR/output"
  cp "$AB/$VARIANT/lhbiem" "$DIR/"
  cp "$REPO/docs/tests/ab_control.in" "$DIR/"
  cp "$REPO/docs/tests/inj_ab_test.txt" "$DIR/"
  echo "=============== running control, $VARIANT ==============="
  cd "$DIR"
  time ./lhbiem ab_control.in
done

echo "=============== byte-for-byte comparison ==============="
STATUS=0
for f in slip sigma pf tau vel kp; do
  A="$AB/ctl_buggy/output/${f}990001.dat"
  B="$AB/ctl_fixed/output/${f}990001.dat"
  if [ -f "$A" ] && [ -f "$B" ]; then
    if cmp -s "$A" "$B"; then
      echo "  $f: IDENTICAL  ($(stat -c%s "$A") bytes)"
    else
      echo "  $f: DIFFER"; STATUS=1
    fi
  fi
done

if [ $STATUS -eq 0 ]; then
  echo "CONTROL PASS -- the fix is a no-op when the clamp never fires"
else
  echo "CONTROL FAIL -- the fix changed a run it should not have"
fi
exit $STATUS
