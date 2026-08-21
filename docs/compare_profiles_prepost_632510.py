#!/usr/bin/env python3
"""Profiles through the injector, pre-fix vs post-fix: sigma_eff, pf and tau.

The companion to 06_slip_profiles_vs_time (cumulative slip). Same construction --
a cut through the injector row at matched times, both binaries overlaid -- but for
the three stress/pressure fields:

    effective normal stress   sigma632510.dat
    pore fluid pressure       pf632510.dat
    shear stress              tau632510.dat

  pre-fix  (recorded)  /scratch/users/nberrios/3dhbi/output/632510/      1402 frames, 0-30.66 d
  post-fix (rerun)     /scratch/users/nberrios/3dhbi/fix_632510/output/  1378 frames, 0-16.38 d

Frames are selected by TIME in each run separately, never by shared index -- the
two runs write frames on different step sequences, so frame k is not the same
instant in both (pre-fix frame 1378 is at ~30 d, post-fix at 16.4 d).

Drawn pre-fix wide / post-fix narrow, so where the runs agree the blue core sits
inside an orange rim rather than one curve hiding the other.

Usage:  python compare_profiles_prepost_632510.py
"""
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

SIGMAINIT = 30.0
MINSIG = 1.0

# validated with scripts/validate_palette.js --mode light: all six checks PASS
PRE, POST = "#D55E00", "#0072BD"
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
LW = {"pre-fix": 3.4, "post-fix": 1.5}
COL = {"pre-fix": PRE, "post-fix": POST}

FIELDS = [
    dict(key="sigma", label=r"effective normal stress $\bar\sigma$ (MPa)",
         short="sigma_eff", refs=[(SIGMAINIT, "sigmainit"), (MINSIG, "minsig")]),
    dict(key="pf", label=r"pore fluid pressure $p_f$ (MPa)",
         short="pf", refs=[(SIGMAINIT, "sigmainit")]),
    dict(key="tau", label=r"shear stress $\tau$ (MPa)",
         short="tau", refs=[]),
]

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


def save(fig, stem):
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {stem}.png/.pdf", flush=True)


def load(base, label):
    tfile = base / f"time{JOBID}.dat"
    if not tfile.exists():
        sys.exit(f"missing {tfile}")
    t_days = np.loadtxt(tfile)[:, 1] / 86400.0
    run = dict(label=label)
    nt = None
    for f in ("sigma", "pf", "tau"):
        p = base / f"{f}{JOBID}.dat"
        if not p.exists():
            sys.exit(f"missing {p}")
        n = min(np.memmap(p, np.float64, "r").size // NCELL, t_days.size)
        nt = n if nt is None else min(nt, n)
        run[f] = np.memmap(p, np.float64, "r", shape=(n, NCELL))
    run["nt"] = nt
    run["t"] = t_days[:nt]
    print(f"{label:9s} {nt} frames, 0-{run['t'][-1]:.3f} d")
    return run


def profile(run, field, t_target):
    """Cut through the injector row at the frame nearest t_target (by time)."""
    k = int(np.argmin(np.abs(run["t"] - t_target)))
    row = np.asarray(run[field][k]).reshape(IMAX, JMAX)[IWELL]
    return row, run["t"][k]


runs = [load(BUGGY, "pre-fix"), load(FIXED, "post-fix")]
pre, post = runs
t_common = min(r["t"][-1] for r in runs)
dist = (np.arange(JMAX) - JWELL) * DS_M
print(f"\ncommon window 0-{t_common:.3f} d\n")

TIMES = [t for t in (1.0, 4.0, 8.0, 12.0, t_common) if t <= t_common + 1e-9]

# ------------------------------------------------- figure 1: fields x times grid
fig, axes = plt.subplots(len(FIELDS), len(TIMES),
                         figsize=(3.15 * len(TIMES), 3.0 * len(FIELDS)),
                         sharex=True, constrained_layout=True)
for ri, fld in enumerate(FIELDS):
    for ci, tt in enumerate(TIMES):
        ax = axes[ri, ci]
        for r in runs:
            prof, tact = profile(r, fld["key"], tt)
            ax.plot(dist, prof, lw=LW[r["label"]], color=COL[r["label"]],
                    label=r["label"])
        for yval, name in fld["refs"]:
            ax.axhline(yval, color=MUTED, lw=0.9, ls="--")
        ax.axvline(0, color=MUTED, lw=0.8, ls=":")
        ax.set_xlim(-900, 900)
        style(ax)
        if ri == 0:
            ax.set_title(f"t = {tt:.2f} d")
        if ri == len(FIELDS) - 1:
            ax.set_xlabel("distance from injector (m)")
        if ci == 0:
            ax.set_ylabel(fld["label"])
            ax.legend(frameon=False, fontsize=8.5, labelcolor=INK)
fig.suptitle("Profiles through the injector, pre-fix vs post-fix (job 632510)",
             fontsize=13)
save(fig, "10_profiles_sigma_pf_tau")

# ------------------------------------- figure 2: at t_common, full span and zoom
fig, axes = plt.subplots(len(FIELDS), 2, figsize=(11.5, 3.1 * len(FIELDS)),
                         constrained_layout=True)
sel = np.abs(dist) <= 300
print(f"at t = {t_common:.2f} d, profile through the injector:")
for ri, fld in enumerate(FIELDS):
    for r in runs:
        prof, _ = profile(r, fld["key"], t_common)
        r["_p"] = prof
        axes[ri, 0].plot(dist, prof, lw=LW[r["label"]], color=COL[r["label"]],
                         label=r["label"])
        axes[ri, 1].plot(dist[sel], prof[sel], lw=LW[r["label"]],
                         color=COL[r["label"]], marker="o", ms=2.6,
                         label=r["label"])
    d = post["_p"] - pre["_p"]
    denom = np.maximum(np.abs(pre["_p"]), 1e-12)
    print(f"  {fld['short']:9s} injector pre {pre['_p'][JWELL]:9.4f}  "
          f"post {post['_p'][JWELL]:9.4f}   max |post-pre| over the profile "
          f"{np.abs(d).max():8.4f} MPa ({100*np.abs(d/denom).max():6.2f}%)")
    for c in (0, 1):
        for yval, name in fld["refs"]:
            axes[ri, c].axhline(yval, color=MUTED, lw=0.9, ls="--")
            if c == 0:
                axes[ri, c].annotate(name, xy=(0.985, yval),
                                     xycoords=("axes fraction", "data"),
                                     xytext=(0, 4), textcoords="offset points",
                                     fontsize=8.5, color=MUTED, ha="right")
        axes[ri, c].axvline(0, color=MUTED, lw=0.8, ls=":")
        style(axes[ri, c])
        axes[ri, c].set_xlabel("distance from injector (m)")
    axes[ri, 0].set_ylabel(fld["label"])
    axes[ri, 0].legend(frameon=False, fontsize=9, labelcolor=INK)
    axes[ri, 0].set_title("full profile" if ri == 0 else "")
    axes[ri, 1].set_title("zoom on the injector" if ri == 0 else "")
fig.suptitle(f"Profiles through the injector at t = {t_common:.2f} d, "
             "pre-fix vs post-fix (job 632510)", fontsize=13)
save(fig, "11_profiles_at_common_time")

print(f"\nfigures -> {OUT}")
