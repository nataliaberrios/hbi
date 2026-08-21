#!/usr/bin/env python3
"""Is the friction law implicated in the limitsigma bug? Three panels say no.

The concern this addresses: shear stress at the injector is wrong by 13.9 MPa
pre-fix, so someone could reasonably ask whether the friction formulation (aging
law, cutoff velocity, flash heating) is at fault rather than the normal-stress
bookkeeping. Each panel isolates that differently.

PANEL A -- the identity test, and the decisive one. For problem 3dp the commit
94da536 claims the effective normal stress satisfies

    sigma = max(sigmainit - pf, minsig)

exactly, for all time and all cells. That equation contains NO friction quantity:
no mu, no a, no b, no dc, no state variable. Measured at t = 16.37 d, post-fix
satisfies it with a residual of EXACTLY 0.0 over all 361201 cells; pre-fix
violates it by up to 23.33 MPa across 2317 cells. A defect that shows up as a
violation of a friction-free identity cannot be a friction defect.

PANEL B -- the friction coefficient itself. mu = tau/sigma along the profile.
Inside the clamped zone, where sigma differs between the two runs by a factor of
24, mu differs by 0.000785 (0.08%). The friction law returned the same
coefficient; it was handed a corrupted normal stress.

PANEL C -- tau against sigma for the clamped-zone cells. Both runs lie in the
same narrow mu band; the pre-fix points simply extend to much larger sigma. That
is what "tau is slaved to sigma through an unchanged mu" looks like.

Usage:  python compare_friction_not_the_issue.py
"""
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BUGGY = Path("/scratch/users/nberrios/3dhbi/output/632510")
FIXED = Path("/scratch/users/nberrios/3dhbi/fix_632510/output")
OUT = Path(__file__).resolve().parent / "figs" / "prepost_632510"

JOBID = 632510
IMAX = JMAX = 601
DS_M = 5.0
IWELL = JWELL = 300
NCELL = IMAX * JMAX
SIGMAINIT, MINSIG = 30.0, 1.0

# validated with scripts/validate_palette.js --mode light: all six checks PASS,
# worst adjacent CVD dE 23.5 (protan). Same pair used throughout this study.
PRE, POST = "#D55E00", "#0072BD"
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
LW = {"pre-fix": 3.4, "post-fix": 1.5}
COL = {"pre-fix": PRE, "post-fix": POST}

plt.rcParams.update({
    "font.size": 10.5, "axes.titlesize": 11.5, "axes.labelsize": 10.5,
    "axes.edgecolor": MUTED, "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
})


def style(ax):
    ax.grid(True, color=GRID, lw=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return ax


def load(base, label):
    t_days = np.loadtxt(base / f"time{JOBID}.dat")[:, 1] / 86400.0
    run = dict(label=label)
    nt = None
    for f in ("sigma", "pf", "tau"):
        p = base / f"{f}{JOBID}.dat"
        if not p.exists():
            sys.exit(f"missing {p}")
        n = min(os.path.getsize(p) // (8 * NCELL), t_days.size)
        nt = n if nt is None else min(nt, n)
        run[f] = np.memmap(p, np.float64, "r", shape=(n, NCELL))
    run["nt"], run["t"] = nt, t_days[:nt]
    print(f"{label:9s} {nt} frames, 0-{run['t'][-1]:.3f} d")
    return run


runs = [load(BUGGY, "pre-fix"), load(FIXED, "post-fix")]
pre, post = runs
t_common = min(r["t"][-1] for r in runs)
dist = (np.arange(JMAX) - JWELL) * DS_M

for r in runs:
    k = int(np.argmin(np.abs(r["t"] - t_common)))
    r["k"] = k
    r["sig2d"] = np.asarray(r["sigma"][k]).reshape(IMAX, JMAX)
    r["pf2d"] = np.asarray(r["pf"][k]).reshape(IMAX, JMAX)
    r["tau2d"] = np.asarray(r["tau"][k]).reshape(IMAX, JMAX)
    r["resid"] = r["sig2d"] - np.maximum(SIGMAINIT - r["pf2d"], MINSIG)
    r["mu"] = r["tau2d"] / r["sig2d"]
    print(f"  {r['label']:9s} identity residual: max |r| = {np.abs(r['resid']).max():.6e} MPa, "
          f"cells with |r| > 1e-9: {int((np.abs(r['resid']) > 1e-9).sum())}")

# the clamped region, defined from the CORRECT run
clamped = post["sig2d"] <= MINSIG * 1.0001
ii, jj = np.meshgrid(np.arange(IMAX), np.arange(JMAX), indexing="ij")
rad = np.hypot(ii - IWELL, jj - JWELL) * DS_M
r_clamp = rad[clamped].max()
dmu = np.abs(post["mu"] - pre["mu"])
print(f"\nclamped zone: {int(clamped.sum())} cells, out to r = {r_clamp:.0f} m")
print(f"max |mu_post - mu_pre| inside the clamped zone: {dmu[clamped].max():.6f}")
print(f"max |sigma_post - sigma_pre| inside the clamped zone: "
      f"{np.abs(post['sig2d'] - pre['sig2d'])[clamped].max():.4f} MPa")

fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.8), constrained_layout=True)

# ---------------------------------------------- A: friction-free identity residual
ax = axes[0]
ax.axvspan(-r_clamp, r_clamp, color=MUTED, alpha=0.10, zorder=0)
for r in runs:
    ax.plot(dist, r["resid"][IWELL], lw=LW[r["label"]], color=COL[r["label"]],
            label=r["label"])
ax.axhline(0, color=MUTED, lw=1, ls="--")
ax.set_xlim(-600, 600)
style(ax).set(xlabel="distance from injector (m)",
              ylabel=r"$\bar\sigma - \max(\sigma_{init}-p_f,\ \sigma_{min})$  (MPa)",
              title="A. Residual of a friction-free identity\n"
                    "post-fix satisfies it exactly; pre-fix is off by 23 MPa")
ax.legend(frameon=False, fontsize=9.5, labelcolor=INK)
ax.annotate("shaded: clamped zone", xy=(0.5, 0.02), xycoords="axes fraction",
            fontsize=8.5, color=MUTED, ha="center")

# ------------------------------------------------------- B: friction coefficient
ax = axes[1]
ax.axvspan(-r_clamp, r_clamp, color=MUTED, alpha=0.10, zorder=0)
for r in runs:
    ax.plot(dist, r["mu"][IWELL], lw=LW[r["label"]], color=COL[r["label"]],
            label=r["label"])
ax.set_xlim(-600, 600)
win = np.abs(dist) <= 600
lo_mu = min(r["mu"][IWELL][win].min() for r in runs)
hi_mu = max(r["mu"][IWELL][win].max() for r in runs)
pad = max(0.01, 0.35 * (hi_mu - lo_mu))
ax.set_ylim(lo_mu - pad, hi_mu + pad)
style(ax).set(xlabel="distance from injector (m)",
              ylabel=r"friction coefficient  $\mu=\tau/\bar\sigma$",
              title="B. The friction coefficient is unchanged\n"
                    f"max difference in the clamped zone: {dmu[clamped].max():.2e}")
ax.legend(frameon=False, fontsize=9.5, labelcolor=INK)

# --------------------------------------------------- C: tau against sigma, clamped
ax = axes[2]
for r in runs:
    ax.scatter(r["sig2d"][clamped], r["tau2d"][clamped], s=5, alpha=0.35,
               color=COL[r["label"]], edgecolors="none", label=r["label"])
lo = np.nanmin([r["mu"][clamped].min() for r in runs])
hi = np.nanmax([r["mu"][clamped].max() for r in runs])
mu_bar = 0.5 * (lo + hi)
xs = np.linspace(0, max(pre["sig2d"][clamped].max(), 1.0) * 1.05, 50)
ax.plot(xs, mu_bar * xs, color=MUTED, lw=1, ls="--", zorder=0)
ax.annotate(rf"$\tau = {mu_bar:.3f}\,\bar\sigma$" "\n"
            rf"($\mu$ spans {lo:.3f}–{hi:.3f} here)",
            xy=(xs[-1], mu_bar * xs[-1]), xytext=(-6, 10),
            textcoords="offset points", fontsize=8.5, color=MUTED, ha="right")
ax.annotate("post-fix: every clamped cell\nsits at the minsig floor",
            xy=(1.0, mu_bar), xytext=(28, 34), textcoords="offset points",
            fontsize=8.5, color=POST,
            arrowprops=dict(arrowstyle="->", color=POST, lw=1))
ax.annotate("pre-fix: same $\\mu$,\nbut $\\bar\\sigma$ up to 24 MPa",
            xy=(pre["sig2d"][clamped].max(), mu_bar * pre["sig2d"][clamped].max()),
            xytext=(-14, -46), textcoords="offset points", fontsize=8.5, color=PRE,
            ha="right", arrowprops=dict(arrowstyle="->", color=PRE, lw=1))
style(ax).set(xlabel=r"effective normal stress $\bar\sigma$ (MPa)",
              ylabel=r"shear stress $\tau$ (MPa)",
              title="C. Clamped-zone cells lie in the same $\\mu$ band\n"
                    r"$\tau$ tracks $\bar\sigma$; only $\bar\sigma$ changed")
leg = ax.legend(frameon=False, fontsize=9.5, labelcolor=INK, loc="upper left")
for h in leg.legend_handles:
    h.set_alpha(1.0)
    h.set_sizes([28])

fig.suptitle(f"Friction is not implicated — job 632510 at t = {t_common:.2f} d",
             fontsize=13)
OUT.mkdir(parents=True, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"12_friction_not_the_issue.{ext}", dpi=200, bbox_inches="tight")
print(f"\nwrote {OUT}/12_friction_not_the_issue.png/.pdf")
