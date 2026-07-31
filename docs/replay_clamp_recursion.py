#!/usr/bin/env python3
"""Test 1 -- identification by replay.

Takes the pore-pressure history p(t) actually recorded by a finished HBI run
and pushes it through two implementations of the operator-split normal-stress
update:

  OLD (main_LH.f90 before the fix)
      Sigma_n   = sigma_{n-1} + p_{n-1}          <- rebuilt from a clamped sigma
      sigma_n   = clamp(Sigma_n - p_n)

  NEW (after the fix)
      Sigma_n   = Sigma_{n-1}                    <- carried as a state
      sigma_n   = clamp(Sigma_n - p_n)

Neither replay reads the recorded sigma. If the OLD recursion reproduces the
recorded sigma trace and the NEW one does not, the mechanism is identified: the
only thing that produced the recorded stress history is the clamp feeding back
into the stored total normal stress.

Both replays assume sigmaconst (problem 3dp), i.e. no elastic normal-stress
change, so the correct answer is exactly max(sigmainit - p, minsig).

Usage:  python replay_clamp_recursion.py [jobid]
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

JOBID = int(sys.argv[1]) if len(sys.argv) > 1 else 632510
BASE = Path(f"/scratch/users/nberrios/3dhbi/output/{JOBID}")
OUT = Path(__file__).resolve().parent / "figs"
OUT.mkdir(parents=True, exist_ok=True)

IMAX = JMAX = 601
SIGMAINIT, MINSIG, MAXSIG = 30.0, 1.0, 300.0
IWELL = JWELL = 300

NCELL = IMAX * JMAX
time_days = np.loadtxt(BASE / f"time{JOBID}.dat")[:, 1] / 86400.0
NT = len(time_days)
row = np.arange(IWELL * JMAX, IWELL * JMAX + JMAX)


def read_row(name):
    m = np.memmap(BASE / f"{name}{JOBID}.dat", np.float64, "r", shape=(NT, NCELL))
    return np.array(m[:, row])


pf = read_row("pf")
sigma_rec = read_row("sigma")          # ground truth: what the run actually stored
clamp = lambda s: np.clip(s, MINSIG, MAXSIG)


def replay(p, mode):
    """p: (nt, ncell). Returns sigma history under the given update rule."""
    nt = p.shape[0]
    sig = np.empty_like(p)
    sig[0] = clamp(SIGMAINIT - p[0])
    Sigma = sig[0] + p[0]              # both rules initialise sigmat the same way
    for n in range(1, nt):
        if mode == "old":
            Sigma = sig[n - 1] + p[n - 1]   # l.1124: rebuilt from a clamped sigma
        # mode == "new": Sigma is a state, and with sigmaconst it never moves
        sig[n] = clamp(Sigma - p[n])
    return sig


sig_old = replay(pf, "old")
sig_new = replay(pf, "new")
sig_ideal = np.maximum(SIGMAINIT - pf, MINSIG)

c = JWELL
print(f"job {JOBID}: replaying {NT} frames at {pf.shape[1]} cells\n")
print("agreement with the RECORDED sigma (which neither replay was shown):")
for name, s in (("OLD rule", sig_old), ("NEW rule", sig_new)):
    d = np.abs(s - sigma_rec)
    print(f"  {name}:  max |diff| = {d.max():8.3f} MPa   "
          f"mean = {d.mean():7.4f}   at injector = {np.abs(s[:,c]-sigma_rec[:,c]).max():8.3f}")

print("\nagreement with the CORRECT answer max(sigmainit - p, minsig):")
for name, s in (("recorded ", sigma_rec), ("OLD rule ", sig_old), ("NEW rule ", sig_new)):
    d = np.abs(s - sig_ideal)
    print(f"  {name}:  max |diff| = {d.max():8.3f} MPa   mean = {d.mean():7.4f}")

print("\nNEW rule invariant  sigma + p == sigmainit  wherever unclamped:")
free = sig_new > MINSIG + 1e-12
print(f"  max |sigma + p - sigmainit| = {np.abs(sig_new[free] + pf[free] - SIGMAINIT).max():.3e} MPa")

# ------------------------------------------------------------------- figure
fig, ax = plt.subplots(2, 1, figsize=(7.5, 6.0), sharex=True, constrained_layout=True)

ax[0].plot(time_days, sigma_rec[:, c], lw=2.6, color="0.75", label="recorded by the run")
ax[0].plot(time_days, sig_old[:, c], lw=1.0, ls="-", color="k", label="OLD rule, replayed")
ax[0].plot(time_days, sig_new[:, c], lw=1.4, ls="--", color="C0", label="NEW rule, replayed")
ax[0].plot(time_days, sig_ideal[:, c], lw=1.0, ls=":", color="C3", label=r"$\max(\sigma_0-p,\sigma_{\min})$")
ax[0].set_ylabel(r"$\bar\sigma$ at injector (MPa)")
ax[0].legend(fontsize=8, loc="upper left")
ax[0].set_title(f"job {JOBID}: replaying the recorded p(t) through both update rules")

ax[1].semilogy(time_days, np.abs(sig_old[:, c] - sigma_rec[:, c]) + 1e-16, lw=1.2,
               color="k", label="|OLD replay - recorded|")
ax[1].semilogy(time_days, np.abs(sig_new[:, c] - sigma_rec[:, c]) + 1e-16, lw=1.2,
               color="C0", label="|NEW replay - recorded|")
ax[1].set_ylabel("discrepancy (MPa)")
ax[1].set_xlabel("time (days)")
ax[1].legend(fontsize=8, loc="center right")

stem = OUT / f"clamp_replay_{JOBID}"
for ext in ("pdf", "png"):
    fig.savefig(f"{stem}.{ext}", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"\nwrote {stem}.pdf/.png")
