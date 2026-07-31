#!/usr/bin/env python3
"""Reproduce the limitsigma normal-stress ratchet from HBI output.

Reads the binary field output of a `3dp` run with `pressurediffusion T` and
checks the stored effective normal stress against the value the effective
stress principle requires, sigma = max(sigmainit - pf, minsig).

Usage:  python verify_clamp_ratchet.py [jobid]

Writes docs/figs/clamp_ratchet_*.{pdf,png} next to this script.
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

# from res632510.in
IMAX = JMAX = 601
DS_KM = 0.005
SIGMAINIT = 30.0
MINSIG = 1.0
RIGID_GPA = 24.0
IWELL = JWELL = 300          # 0-based; injection file gives Fortran (301,301)

NCELL = IMAX * JMAX
time_days = np.loadtxt(BASE / f"time{JOBID}.dat")[:, 1] / 86400.0
NT = len(time_days)

# one grid row through the injector: flat index = i*jmax + j
row = np.arange(IWELL * JMAX, IWELL * JMAX + JMAX)


def read_row(name):
    m = np.memmap(BASE / f"{name}{JOBID}.dat", np.float64, "r", shape=(NT, NCELL))
    return np.array(m[:, row])


def read_frame(name, k):
    m = np.memmap(BASE / f"{name}{JOBID}.dat", np.float64, "r", shape=(NT, NCELL))
    return np.array(m[k]).reshape(IMAX, JMAX)


pf, sigma, slip, vel, tau = (read_row(n) for n in ("pf", "sigma", "slip", "vel", "tau"))

# what the effective stress principle requires, given the floor
sigma_true = np.maximum(SIGMAINIT - pf, MINSIG)
err = sigma - sigma_true

# the total normal stress the code is carrying in `sigmat`
sigmat_final = sigma[-1] + pf[-1]
pf_running_max = pf.max(axis=0)

dist_m = (np.arange(JMAX) - JWELL) * DS_KM * 1e3
c, ref = JWELL, JWELL + 17          # injector and the r = 85 m peak

# ---------------------------------------------------------------- reporting
print(f"job {JOBID}: {NT} frames, 0-{time_days[-1]:.2f} d, grid {IMAX}x{JMAX}")
print("\nratchet identity  sigmat_final == minsig + max_t pf   (MPa)")
print(f"{'j':>5}{'r (m)':>8}{'max pf':>9}{'sigmat':>9}{'minsig+maxpf':>14}")
for j in [300, 305, 310, 317, 325, 340, 360, 400, 450]:
    print(f"{j:5d}{dist_m[j]:8.0f}{pf_running_max[j]:9.2f}"
          f"{sigmat_final[j]:9.2f}{MINSIG + pf_running_max[j]:14.2f}")

clamped = np.where(pf_running_max > SIGMAINIT - MINSIG)[0]
resid = np.abs(sigmat_final[clamped] - (MINSIG + pf_running_max[clamped])).max()
print(f"\nmax |sigmat - (minsig + max pf)| over clamped cells: {resid:.2e} MPa")
print(f"sigmat at a never-clamped cell (j=450):             {sigmat_final[450]:.4f} MPa")
print(f"clamped region: |r| <= {dist_m[clamped].max():.0f} m")

# slip-deficit accounting: injector vs the r = 85 m cell
deficit = (slip[:, ref] - slip[:, c]) * 100.0            # cm
d_def = np.diff(deficit)
bad = err[:-1, c] > 0.5
print(f"\nslip deficit gained while sigma is corrupted : {d_def[bad].sum():+.3f} cm")
print(f"slip deficit gained while sigma is correct   : {d_def[~bad].sum():+.3f} cm")
print(f"net deficit (= depth of the cusp)            : {deficit[-1]:+.3f} cm")

# elastic back-stress implied by that deficit, dtau ~ mu * delta / a
a_km = 0.1
print(f"\nimplied back-stress mu*delta/a for delta={deficit[-1]/100:.4f} m, a={a_km} km:"
      f" {RIGID_GPA * (deficit[-1] / 100) / a_km:.2f} MPa")

# ------------------------------------------------------------------ fig. 1
fig, ax = plt.subplots(3, 1, figsize=(7.0, 8.0), sharex=True, constrained_layout=True)

s_final = read_frame("slip", -1)[IWELL] * 100.0
ax[0].plot(dist_m, s_final, lw=1.8, color="k")
ax[0].axvline(0, color="0.6", lw=0.8, ls=":")
ax[0].set_ylabel("cumulative slip (cm)")
ax[0].set_title(f"job {JOBID}: profile through the injector, final state")
ax[0].set_ylim(bottom=0)

ax[1].plot(dist_m, sigmat_final, lw=1.8, color="k", label=r"$\Sigma$ carried by the code")
ax[1].plot(dist_m, MINSIG + pf_running_max, lw=1.0, ls="--", color="C3",
           label=r"$\sigma_{\min}+\max_t p$")
ax[1].axhline(SIGMAINIT, color="C0", lw=1.0, ls=":", label=r"$\sigma_{\rm init}$ (correct)")
ax[1].set_ylabel("total normal stress (MPa)")
ax[1].legend(fontsize=8)

ax[2].plot(dist_m, err.max(axis=0), lw=1.8, color="C3")
ax[2].axhline(0, color="0.6", lw=0.8)
ax[2].set_ylabel(r"max$_t\;(\sigma_{\rm code}-\sigma_{\rm true})$ (MPa)")
ax[2].set_xlabel("distance from injector along the profile (m)")
ax[2].set_xlim(-600, 600)

f1 = OUT / f"clamp_ratchet_profile_{JOBID}"
for ext in ("pdf", "png"):
    fig.savefig(f"{f1}.{ext}", dpi=200, bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------------ fig. 2
fig, ax = plt.subplots(3, 1, figsize=(7.0, 7.5), sharex=True, constrained_layout=True)

ax[0].plot(time_days, pf[:, c], lw=1.2, color="C0", label=r"$p$ at injector")
ax[0].plot(time_days, sigma[:, c], lw=1.4, color="k", label=r"$\sigma$ stored by the code")
ax[0].plot(time_days, sigma_true[:, c], lw=1.2, ls="--", color="C3",
           label=r"$\max(\sigma_{\rm init}-p,\ \sigma_{\min})$")
ax[0].set_ylabel("stress (MPa)")
ax[0].legend(fontsize=8, loc="upper left")
ax[0].set_title(f"job {JOBID}: injector cell ({IWELL+1},{JWELL+1}) in Fortran indexing")

ax[1].semilogy(time_days, np.maximum(vel[:, c], 1e-30), lw=1.0, color="k", label="injector")
ax[1].semilogy(time_days, np.maximum(vel[:, ref], 1e-30), lw=1.0, color="C1",
               label=f"r = {dist_m[ref]:.0f} m")
ax[1].set_ylabel("slip rate (m/s)")
ax[1].set_ylim(1e-28, 1e-4)
ax[1].legend(fontsize=8, loc="lower left")

ax[2].plot(time_days, deficit, lw=1.4, color="k")
ax[2].axhline(0, color="0.6", lw=0.8)
ax[2].fill_between(time_days[:-1], deficit.min(), deficit.max(), where=bad,
                   color="C3", alpha=0.15, lw=0, label=r"$\sigma$ corrupted at injector")
ax[2].set_ylabel("slip deficit,\ninjector vs r=85 m (cm)")
ax[2].set_xlabel("time (days)")
ax[2].legend(fontsize=8, loc="upper left")

f2 = OUT / f"clamp_ratchet_history_{JOBID}"
for ext in ("pdf", "png"):
    fig.savefig(f"{f2}.{ext}", dpi=200, bbox_inches="tight")
plt.close(fig)

print(f"\nwrote {f1}.pdf/.png\nwrote {f2}.pdf/.png")
