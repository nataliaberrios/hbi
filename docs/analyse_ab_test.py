#!/usr/bin/env python3
"""Test 2 -- end-to-end A/B comparison of the pre-fix and post-fix binaries.

Reads the output of docs/tests/run_ab_test.sbatch, which builds two binaries
differing only in the limitsigma normal-stress update and runs them on the
identical input deck (docs/tests/ab_test.in).

Checks, for problem 3dp where sigmaconst is .true. and the total normal stress
is therefore forbidden to change:

  T1  invariant   sigma == max(sigmainit - pf, minsig)   at every cell and step
  T2  total normal stress sigma + pf == sigmainit wherever the clamp is inactive
  T3  the two binaries agree bit-for-bit until the clamp first fires
  T4  cumulative slip is maximised at the injector, not on a ring around it

Usage:  python analyse_ab_test.py [ab_dir]
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

AB = Path(sys.argv[1] if len(sys.argv) > 1 else "/scratch/users/nberrios/3dhbi/ab_test")
OUT = Path(__file__).resolve().parent / "figs"
OUT.mkdir(parents=True, exist_ok=True)

JOBID = 990000
IMAX = JMAX = 151
DS_KM = 0.005
SIGMAINIT, MINSIG = 30.0, 1.0
IWELL = JWELL = 75            # 0-based; deck puts the well at Fortran (76,76)
NCELL = IMAX * JMAX
WELL = IWELL * JMAX + JWELL


def load(variant):
    d = AB / variant / "output"
    t = np.loadtxt(d / f"time{JOBID}.dat")
    nt = len(t)
    out = {"t_days": t[:, 1] / 86400.0, "nt": nt}
    for name in ("pf", "sigma", "slip", "vel", "tau"):
        a = np.fromfile(d / f"{name}{JOBID}.dat", np.float64)
        n = a.size // NCELL
        out[name] = a[: n * NCELL].reshape(n, NCELL)
        out["nt"] = min(out["nt"], n)
    n = out["nt"]
    for name in ("pf", "sigma", "slip", "vel", "tau"):
        out[name] = out[name][:n]
    out["t_days"] = out["t_days"][:n]
    return out


runs = {v: load(v) for v in ("buggy", "fixed")}
fail = []
print(f"grid {IMAX}x{JMAX}, sigmainit={SIGMAINIT} MPa, minsig={MINSIG} MPa")
for v, r in runs.items():
    print(f"  {v:5s}: {r['nt']} frames, 0-{r['t_days'][-1]:.4f} d, "
          f"max pf = {r['pf'].max():.2f} MPa, max slip = {r['slip'][-1].max()*100:.4f} cm")

n = min(r["nt"] for r in runs.values())
print(f"\ncomparing over the {n} frames both runs share")

# ---------------------------------------------------------------- T1 and T2
print("\nT1  sigma == max(sigmainit - pf, minsig)")
for v, r in runs.items():
    ideal = np.maximum(SIGMAINIT - r["pf"][:n], MINSIG)
    e = np.abs(r["sigma"][:n] - ideal).max()
    ok = e < 1e-9
    print(f"    {v:5s}: max |error| = {e:10.4f} MPa   {'PASS' if ok else 'FAIL'}")
    if v == "fixed" and not ok:
        fail.append("T1")

print("\nT2  sigma + pf == sigmainit wherever the clamp is inactive")
for v, r in runs.items():
    free = r["sigma"][:n] > MINSIG + 1e-12
    e = np.abs(r["sigma"][:n][free] + r["pf"][:n][free] - SIGMAINIT).max()
    ok = e < 1e-9
    print(f"    {v:5s}: max |sigma + pf - sigmainit| = {e:10.4f} MPa   "
          f"{'PASS' if ok else 'FAIL'}")
    if v == "fixed" and not ok:
        fail.append("T2")

# ------------------------------------------------------------------- T3
print("\nT3  the binaries agree until the clamp first fires")
clamped = runs["buggy"]["sigma"][:n] <= MINSIG + 1e-12
if clamped.any():
    first = int(np.argmax(clamped.any(axis=1)))
    print(f"    clamp first fires at frame {first} (t = {runs['buggy']['t_days'][first]:.5f} d)")
    agree = True
    for name in ("pf", "sigma", "slip", "tau"):
        d = np.abs(runs["buggy"][name][:first] - runs["fixed"][name][:first]).max() if first else 0.0
        agree &= d == 0.0
        print(f"      {name:6s} max |buggy - fixed| before first clamp = {d:.3e}")
    print(f"    {'PASS -- identical before the clamp fires' if agree else 'FAIL -- differ beforehand'}")
    if not agree:
        fail.append("T3")
else:
    print("    clamp never fired; this deck does not exercise the bug")

# ------------------------------------------------------------------- T4
print("\nT4  cumulative slip peaks at the injector")
for v, r in runs.items():
    s = r["slip"][n - 1].reshape(IMAX, JMAX)
    mi, mj = np.unravel_index(np.argmax(s), s.shape)
    off = np.hypot(mi - IWELL, mj - JWELL) * DS_KM * 1e3
    # within one cell of the injector: a 1-cell offset is discretisation, not a
    # ring.  The pre-fix cusp put the peak 17 cells (85 m) out, so this still
    # separates the two cases.
    ok = off <= DS_KM * 1e3 * 1.001
    print(f"    {v:5s}: max slip at ({mi},{mj}), {off:6.1f} m from the injector; "
          f"injector {s[IWELL,JWELL]*100:.4f} cm vs peak {s[mi,mj]*100:.4f} cm   "
          f"{'PASS' if ok else 'FAIL'}")
    if v == "fixed" and not ok:
        fail.append("T4")

print("\n" + ("ALL CHECKS PASS for the fixed binary"
              if not fail else f"FAILURES in the fixed binary: {sorted(set(fail))}"))

# ---------------------------------------------------------------- figures
dist = (np.arange(JMAX) - JWELL) * DS_KM * 1e3
fig, ax = plt.subplots(1, 3, figsize=(14, 3.8), constrained_layout=True)
col = {"buggy": "C3", "fixed": "C0"}

for v, r in runs.items():
    s = r["slip"][n - 1].reshape(IMAX, JMAX)[IWELL] * 100
    ax[0].plot(dist, s, lw=1.8, color=col[v], label=v)
ax[0].axvline(0, color="0.6", lw=0.8, ls=":")
ax[0].set(xlabel="distance from injector (m)", ylabel="cumulative slip (cm)",
          title="final slip through the injector")
ax[0].set_xlim(-250, 250)
ax[0].legend(fontsize=9)

for v, r in runs.items():
    ax[1].plot(r["t_days"][:n], r["sigma"][:n, WELL], lw=1.6, color=col[v], label=v)
ideal = np.maximum(SIGMAINIT - runs["fixed"]["pf"][:n, WELL], MINSIG)
ax[1].plot(runs["fixed"]["t_days"][:n], ideal, lw=1.0, ls="--", color="k",
           label=r"$\max(\sigma_0-p,\sigma_{\min})$")
ax[1].set(xlabel="time (days)", ylabel=r"$\bar\sigma$ at injector (MPa)",
          title="effective normal stress")
ax[1].legend(fontsize=9)

for v, r in runs.items():
    ax[2].plot(r["t_days"][:n], r["sigma"][:n, WELL] + r["pf"][:n, WELL], lw=1.6,
               color=col[v], label=v)
ax[2].axhline(SIGMAINIT, color="k", lw=1.0, ls="--", label=r"$\sigma_{\rm init}$")
ax[2].set(xlabel="time (days)", ylabel=r"$\Sigma=\bar\sigma+p$ at injector (MPa)",
          title="total normal stress (must not move)")
ax[2].legend(fontsize=9)

stem = OUT / "ab_test_comparison"
for ext in ("pdf", "png"):
    fig.savefig(f"{stem}.{ext}", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"\nwrote {stem}.pdf/.png")

sys.exit(1 if fail else 0)
