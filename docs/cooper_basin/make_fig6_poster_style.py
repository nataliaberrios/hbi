#!/usr/bin/env python3
"""Fig 6 in the poster's own layout and styling, with updated runs.

The poster does not stack three panels. It has two images:

  image16  a 1x2 side-by-side, observed on the left and the simulation on the
           right, figsize (16, 3) at dpi 300, label fontsize 13, title 14,
           legend 9, both panels forced to a shared y limit. Those numbers come
           straight from cell 107 of
           notebooks/cooper_basin_plots-27_abs_pressure.ipynb, which is the cell
           that produced it.
  image18  a single wide panel for the second simulation, larger fonts.

This reproduces both with the new runs: 632894 (Job 807's physics on current
code) in the pair, and 632895 (632510 with the limitsigma cusp gone, now
reaching a true 17 d) as the wide panel.

Kept from the original: viridis over np.linspace(0, 1, 8), the 3-17 day set in
2-day steps, xlim +/-1.5 km, cm on the y axis, the shared y limit across the
pair.

Changed, deliberately, and both are flagged in the printed summary:
  * the simulated profiles are along STRIKE, frame[:, c]. The original used
    frame[c, :], which is along DIP -- j is the innermost loop in
    coordinate3ddip (main_LH.f90:1577) and controls y,z. The observed panel is a
    real geographic along-strike profile, so the original compared two different
    directions. --axis dip reproduces the poster exactly.
  * no time extrapolation. The original interpolated with
    fill_value="extrapolate", so a run ending at 16.37 d still yielded a curve
    labelled "17 days". Both runs now reach 17.000 d, asserted here.

Usage:  python make_fig6_poster_style.py [--axis strike|dip]
"""
import argparse
import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OBS = Path("/home/users/nberrios/3dhbi/hbi/slip_profiles_strike.txt")
SCRATCH = Path("/scratch/users/nberrios/3dhbi/output")
OUT = Path("/home/users/nberrios/3dhbi/hbi_analysis/figures/fig6_remake")
PAIR_SIM, WIDE_SIM = 632894, 632895
IM = JM = 601
DS_M = 5.0
TIMES_D = [3, 5, 7, 9, 11, 13, 15, 17]


def sim_profiles(job, axis):
    base = SCRATCH / str(job)
    p = base / f"slip{job}.dat"
    t = np.atleast_2d(np.loadtxt(base / f"time{job}.dat"))[:, 1] / 86400.0
    nt = min(os.path.getsize(p) // (8 * IM * JM), len(t))
    arr = np.memmap(p, np.float64, "r", shape=(nt, IM * JM))
    c = (IM - 1) // 2
    x = (np.arange(IM) - c) * DS_M / 1000.0
    cols = []
    for td in TIMES_D:
        k = int(np.argmin(np.abs(t[:nt] - td)))
        assert abs(t[k] - td) < 0.05, (
            f"run {job}: nearest frame to {td} d is {t[k]:.3f} d; the original "
            f"figure would have extrapolated here")
        f = np.asarray(arr[k]).reshape(IM, JM)
        cols.append((f[:, c] if axis == "strike" else f[c, :]) * 100.0)
    return x, np.column_stack(cols)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--axis", choices=("strike", "dip"), default="strike")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    xlab = f"Distance along-{a.axis} (km)"

    o = np.loadtxt(OBS)
    x_obs, obs = o[:, 0], o[:, 1:]
    x_p, pair = sim_profiles(PAIR_SIM, a.axis)
    x_w, wide = sim_profiles(WIDE_SIM, a.axis)
    colors = plt.cm.viridis(np.linspace(0, 1, len(TIMES_D)))

    # ---- the side-by-side pair, image16's geometry
    fig, (ax_obs, ax_sim) = plt.subplots(1, 2, figsize=(16, 3), dpi=300)
    for i, td in enumerate(TIMES_D):
        ax_obs.plot(x_obs, obs[:, i], lw=2, color=colors[i], label=f"{td} days")
        ax_sim.plot(x_p, pair[:, i], "-", lw=2, color=colors[i], alpha=0.8,
                    label=f"{td} days")
    ax_obs.set_xlabel(xlab, fontsize=13)
    ax_obs.set_ylabel("Cumulative slip (cm)", fontsize=13)
    ax_obs.set_xlim([-1.5, 1.5])
    ax_obs.legend(fontsize=9, loc="best")
    ax_obs.set_title("Observed, 20 m west of well", fontsize=14)
    ax_sim.set_xlabel(xlab, fontsize=13)
    ax_sim.set_ylabel("Cumulative slip (cm)", fontsize=13)
    ax_sim.set_xlim([-1.5, 1.5])
    ax_sim.legend(fontsize=9, loc="upper right")
    ax_sim.set_title(f"Simulation (Job {PAIR_SIM})", fontsize=14)
    ymax = max(ax_obs.get_ylim()[1], ax_sim.get_ylim()[1])
    ax_obs.set_ylim([0, ymax])
    ax_sim.set_ylim([0, ymax])
    plt.tight_layout()
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"fig6_pair_{a.axis}.{e}", bbox_inches="tight")
    plt.close(fig)

    # ---- the wide single panel, image18's geometry
    fig, ax = plt.subplots(figsize=(10.5, 3.9), dpi=300)
    for i, td in enumerate(TIMES_D):
        ax.plot(x_w, wide[:, i], lw=2.2, color=colors[i], label=f"{td} days")
    ax.set_xlabel(xlab, fontsize=17)
    ax.set_ylabel("Cumulative slip (cm)", fontsize=17)
    ax.set_xlim([-1.5, 1.5])
    ax.set_ylim(bottom=0)
    ax.tick_params(labelsize=14)
    ax.legend(fontsize=13, loc="upper right")
    plt.tight_layout()
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"fig6_wide_{a.axis}.{e}", bbox_inches="tight")
    plt.close(fig)

    print(f"axis: along-{a.axis}")
    print(f"wrote {OUT}/fig6_pair_{a.axis}.png   (16x3, matches image16)")
    print(f"wrote {OUT}/fig6_wide_{a.axis}.png   (matches image18)")
    print(f"\nshared y limit on the pair: 0 to {ymax:.2f} cm")
    c = len(x_w) // 2
    print(f"\n{'t (d)':>6s} {'observed':>9s} {f'sim {PAIR_SIM}':>12s} "
          f"{f'sim {WIDE_SIM}':>12s} {'cusp dip %':>11s}")
    for i, td in enumerate(TIMES_D):
        print(f"{td:>6d} {obs[:, i].max():>9.3f} {pair[:, i].max():>12.3f} "
              f"{wide[:, i].max():>12.3f} "
              f"{100*(wide[:,i].max()-wide[c,i])/wide[:,i].max():>11.2f}")


if __name__ == "__main__":
    main()
