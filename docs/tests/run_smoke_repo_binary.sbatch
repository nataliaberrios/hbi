#!/bin/bash
#SBATCH --partition=serc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --constraint=CPU_GEN:RME
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --job-name=smokerepo
#SBATCH --output=/scratch/users/nberrios/3dhbi/smoke_repo_%j.log

# Regression test on THE BINARY THE GRID-SEARCH HARNESS ACTUALLY RUNS.
#
# Verifying "the source has the fix and the binary is newer" is not the same as
# verifying the binary behaves. This runs hbi_git/lhbiem -- the exact file
# EXEC_PATH points at in march26_submit_hbi_git.sh -- on the A/B control deck.
#
# Sizing, from the recorded A/B pair (ab_test/{buggy,fixed}/output/monitor990000.dat):
#   step 1114  the clamp starts binding (minnorm reaches minsig = 1)
#   step 8309  the buggy binary's maxnorm first exceeds sigmainit (30.032)
#   step 14078 the fixed binary is still at exactly 30.000
# So nstep 10000 clears the detection threshold with margin.
#
# PASS  = maxnorm stays at 30.000 for all steps AND the clamp actually bound
# FAIL  = maxnorm exceeds sigmainit  -> the harness is running pre-fix code
# INCONCLUSIVE = clamp never bound   -> deck did not exercise the bug

set -euo pipefail
ml purge; ml load ifort; ml load impi

REPO=/home/users/nberrios/3dhbi/hbi_git
RUN=/scratch/users/nberrios/3dhbi/smoke_repo
rm -rf "$RUN"; mkdir -p "$RUN/output"; cd "$RUN"

cp "$REPO/lhbiem" ./lhbiem_under_test
cp "$REPO/res990009_smoke.in" .
cp "$REPO/docs/tests/inj_ab_test.txt" .

echo "binary under test: $REPO/lhbiem"
ls -l --time-style=+'%Y-%m-%d %H:%M' "$REPO/lhbiem"
sha256sum ./lhbiem_under_test | cut -c1-12

export OMP_NUM_THREADS=16
time ./lhbiem_under_test smoke_repo_binary.in | tail -5

M=output/monitor990009.dat
echo "=============================== VERDICT ==============================="
awk 'NR>1{
  if($7<=1.0001) bound=1
  if($6>30.001 && !bad){bad=1; badstep=$1; badval=$6}
  last=$1; lm=$6; ln=$7
}
END{
  printf "steps run: %d   final maxnorm: %.4f   final minnorm: %.4f\n", last, lm, ln
  printf "clamp bound at some point: %s\n", (bound?"YES":"NO")
  if(!bound) print "INCONCLUSIVE -- the clamp never bound, deck did not exercise the bug"
  else if(bad) printf "FAIL -- maxnorm reached %.4f at step %d; this binary RATCHETS\n", badval, badstep
  else print "PASS -- maxnorm never left sigmainit; the harness binary carries the fix"
}' $M
echo "======================================================================="

printf '| %s | smokerepo | COMPLETED | - | `%s` — regression test of `hbi_git/lhbiem` (the harness binary) |\n' \
  "$SLURM_JOB_ID" "/scratch/users/nberrios/3dhbi/smoke_repo_${SLURM_JOB_ID}.log" \
  >> /home/users/nberrios/3dhbi/JOB_STATUS.md
