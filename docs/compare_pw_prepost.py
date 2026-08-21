#!/usr/bin/env python3
"""Wellbore pressure (pw) time series, pre-fix vs post-fix, job 632510.

pw is written by output_monitor (main_LH.f90:1337) as  step, time_s, pw_MPa  --
one row per monitor step, so ~40k rows rather than the ~1400 field frames. It is
the wellbore pressure from the well model (rw, skin, Sw_fwid), which sits above
the injector cell's pf.

Why this needs care: the two runs end at different times (pre-fix 30.66 d,
post-fix 16.37 d), so the last rows are NOT comparable -- 10.08 MPa vs 163.05
MPa is almost entirely the 14 d of post-injection decay the post-fix run never
reached. Everything here is interpolated onto a common time grid. Comparing row
index to row index on this pair produces nonsense.

Usage:  python compare_pw_prepost.py
"""
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BUGGY = Path("/scratch/users/nberrios/3dhbi/output/632510")
FIXED = Path("/scratch/users/nberrios/3dhbi/fix_632510/output")
OUT = Path(__file__).resolve().parent / "figs" / "prepost_632510"

JOBID = 632510
INJ_END_D = 17.15                # last node with a nonzero injection rate

# validated with scripts/validate_palette.js --mode light: all six checks PASS
PRE, POST = "#D55E00", "#0072BD"
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
# pre-fix drawn wide, post-fix narrow on top: where only the blue core shows
# inside an orange rim, the two runs agree
LW = {"pre-fix": 3.4, "post-fix": 1.5}
COL = {"pre-fix": PRE, "post-fix": POST}

plt.rcParams.update({
    "font.size": 11, "axes.titlesize": 12.5, "axes.labelsize": 11,
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
    a = np.loadtxt(base / f"pw{JOBID}.dat")
    t = a[:, 1] / 86400.0
    pw = a[:, 2]
    print(f"{label:9s} {len(t):6d} rows, 0-{t[-1]:.3f} d, "
          f"pw range [{pw.min():.3f}, {pw.max():.3f}] MPa, final {pw[-1]:.3f}")
    return dict(label=label, t=t, pw=pw)


runs = [load(BUGGY, "pre-fix"), load(FIXED, "post-fix")]
pre, post = runs

t_common = min(r["t"][-1] for r in runs)
grid = np.linspace(0.0, t_common, 4000)
pw_pre = np.interp(grid, pre["t"], pre["pw"])
pw_post = np.interp(grid, post["t"], post["pw"])
d = pw_post - pw_pre

print(f"\ncommon window 0-{t_common:.3f} d, {grid.size} interpolated points")
sel = np.abs(pw_pre) > 0.01 * np.abs(pw_pre).max()
print(f"max |post - pre|      {np.abs(d).max():.4f} MPa")
print(f"max relative diff     {100*np.abs(d[sel]/pw_pre[sel]).max():.3f}%")
print(f"RMS difference        {np.sqrt(np.mean(d**2)):.4f} MPa")
print(f"peak pw, pre  {pw_pre.max():.3f} MPa   post {pw_post.max():.3f} MPa")
print("\nfor reference, the naive row-index comparison would report "
      f"{abs(post['pw'][-1] - pre['pw'][-1]):.1f} MPa -- all of it the "
      "different end times")

# ------------------------------------------------------------------- figure
fig, axes = plt.subplots(3, 1, figsize=(9, 10.5), constrained_layout=True)

# full span, on each run's own clock
for r in runs:
    axes[0].plot(r["t"], r["pw"], lw=LW[r["label"]], color=COL[r["label"]],
                 label=r["label"])
axes[0].axvline(INJ_END_D, color=MUTED, lw=1, ls=":")
axes[0].annotate("injection rate → 0", xy=(INJ_END_D, 0.97),
                 xycoords=("data", "axes fraction"), xytext=(-4, 0),
                 textcoords="offset points", fontsize=9, color=MUTED,
                 rotation=90, ha="right", va="top")
style(axes[0]).set(xlabel="time (days)", ylabel="wellbore pressure (MPa)",
                   title="pw over each run's full span — the pre-fix tail is post-injection decay")
axes[0].legend(frameon=False, fontsize=9, labelcolor=INK)

# first two days, where the injection transients live
m_pre = pre["t"] <= 2.0
m_post = post["t"] <= 2.0
axes[1].plot(pre["t"][m_pre], pre["pw"][m_pre], lw=LW["pre-fix"], color=PRE, label="pre-fix")
axes[1].plot(post["t"][m_post], post["pw"][m_post], lw=LW["post-fix"], color=POST, label="post-fix")
style(axes[1]).set(xlabel="time (days)", ylabel="wellbore pressure (MPa)",
                   title="First 2 days — the startup transient")
axes[1].legend(frameon=False, fontsize=9, labelcolor=INK)

# difference on the common grid
axes[2].axhline(0, color=MUTED, lw=1, ls="--")
axes[2].plot(grid, d, lw=1.6, color=INK)
style(axes[2]).set(xlabel="time (days)", ylabel="post − pre (MPa)",
                   title=f"Difference on a common time grid — max {np.abs(d).max():.3f} MPa "
                         f"({100*np.abs(d[sel]/pw_pre[sel]).max():.2f}% peak relative)")

OUT.mkdir(parents=True, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"09_pw_prepost.{ext}", dpi=200, bbox_inches="tight")
print(f"\nwrote {OUT}/09_pw_prepost.png/.pdf")
