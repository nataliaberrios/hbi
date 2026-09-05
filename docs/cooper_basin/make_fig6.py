#!/usr/bin/env python3
"""Fig 6 remade: observed, Job 807's physics on current code, and 632510 to 17 d.

Three stacked panels, as the poster caption describes:

  top     observed cumulative slip 20 m west of the well, from the seismicity
          catalog under a circular crack model. Read from slip_profiles_strike.txt,
          the same precomputed file the original figure used, so this panel is
          unchanged.
  middle  632894 -- Job 807's physics on the current code. res807.in could not be
          rerun: shear_mod is no longer a recognised key (so 807 ran at
          rigid = 32.04 rather than 24), Sw_fwid predates the wellbore model and
          had defaulted to 1.0, and oct_clean_300.txt is in the older injection
          format the current reader would misparse.
  bottom  632895 -- 632510 with the limitsigma cusp gone. The poster's version is
          pre-fix, where the ratchet dimples the injector by up to 7.8%.

TWO DEPARTURES FROM THE ORIGINAL, both deliberate.

1. AXIS. The simulated panels here are along STRIKE, frame[:, c]. The original
   plotted frame[c, :], which is along DIP -- in coordinate3ddip
   (main_LH.f90:1577) j is innermost and controls y,z, so reshape(IM,JM)[a,b] has
   a = STRIKE and b = DIP. Verified against the poster: image18 matches the dip
   profile (20.487 cm, 770 m at 17 d) rather than strike (20.426 cm, 840 m).
   Since the observed panel is a real geographic along-strike profile, comparing
   it against a simulated dip profile is inconsistent, so the simulated panels
   are switched to strike. It costs little: peaks agree to three decimals
   between the two axes and only the extent moves, 840 vs 775 m for 632895.

2. NO TIME EXTRAPOLATION. The original interpolated with
   fill_value="extrapolate", so a run ending at 16.37 d still produced a curve
   labelled "17 days". Both runs here reach 17.000 d exactly, and this script
   asserts that rather than extrapolating.

Usage:  python make_fig6.py [--axis strike|dip]
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
MID, BOT = 632894, 632895
IM = JM = 601
DS_M = 5.0
TIMES_D = [3, 5, 7, 9, 11, 13, 15, 17]
XLIM = 1.5


def sim_profiles(job, axis):
    """(x_km, slip_cm[len(x), len(TIMES_D)]) at the frames nearest TIMES_D."""
    base = SCRATCH / str(job)
    p = base / f"slip{job}.dat"
    t = np.atleast_2d(np.loadtxt(base / f"time{job}.dat"))[:, 1] / 86400.0
    nt = min(os.path.getsize(p) // (8 * IM * JM), len(t))
    arr = np.memmap(p, np.float64, "r", shape=(nt, IM * JM))
    c = (IM - 1) // 2
    x = (np.arange(IM) - c) * DS_M / 1000.0
    cols, actual = [], []
    for td in TIMES_D:
        k = int(np.argmin(np.abs(t[:nt] - td)))
        # Assert rather than extrapolate. The original figure used
        # fill_value="extrapolate", which drew a "17 days" curve from a run that
        # stopped at 16.37 d.
        assert abs(t[k] - td) < 0.05, (
            f"run {job} has no frame within 0.05 d of {td} d "
            f"(nearest {t[k]:.3f} d, run ends {t[nt-1]:.3f} d)")
        frame = np.asarray(arr[k]).reshape(IM, JM)
        cols.append((frame[:, c] if axis == "strike" else frame[c, :]) * 100.0)
        actual.append(t[k])
    return x, np.column_stack(cols), np.array(actual)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--axis", choices=("strike", "dip"), default="strike")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    o = np.loadtxt(OBS)
    x_obs, obs = o[:, 0], o[:, 1:]
    x_mid, mid, t_mid = sim_profiles(MID, a.axis)
    x_bot, bot, t_bot = sim_profiles(BOT, a.axis)
    print(f"axis: along-{a.axis}")
    print(f"  frames used, middle: {np.round(t_mid,3)}")
    print(f"  frames used, bottom: {np.round(t_bot,3)}")

    cmap = plt.get_cmap("viridis")
    cols = [cmap(v) for v in np.linspace(0.0, 1.0, len(TIMES_D))]
    panels = [("Observed, 20 m west of well", x_obs, obs),
              (f"Simulation, Job 807 physics on current code (run {MID})",
               x_mid, mid),
              (f"Simulation, cusp fixed (run {BOT})", x_bot, bot)]

    fig, axes = plt.subplots(3, 1, figsize=(9.0, 10.2), dpi=200,
                             constrained_layout=True)
    for ax, (title, x, Y) in zip(axes, panels):
        for i, td in enumerate(TIMES_D):
            ax.plot(x, Y[:, i], lw=1.9, color=cols[i], label=f"{td} days")
        ax.set(ylabel="Cumulative slip (cm)", xlim=(-XLIM, XLIM), title=title)
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=8, ncol=2, loc="upper right", frameon=True)
    # top two share a y-limit, as in the original; the bottom is 2.5x larger and
    # would flatten them
    ymax = max(axes[0].get_ylim()[1], axes[1].get_ylim()[1])
    axes[0].set_ylim(0, ymax)
    axes[1].set_ylim(0, ymax)
    axes[2].set_xlabel(f"Distance along-{a.axis} (km)")
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"fig6_remake_{a.axis}.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/fig6_remake_{a.axis}.png")

    c = len(x_mid) // 2
    print(f"\n{'t (d)':>6s} {'observed':>9s} {'632894':>9s} {'632895':>9s} "
          f"{'894 cusp%':>10s} {'895 cusp%':>10s}")
    for i, td in enumerate(TIMES_D):
        print(f"{td:>6d} {obs[:, i].max():>9.3f} {mid[:, i].max():>9.3f} "
              f"{bot[:, i].max():>9.3f} "
              f"{100*(mid[:,i].max()-mid[c,i])/mid[:,i].max():>10.2f} "
              f"{100*(bot[:,i].max()-bot[c,i])/bot[:,i].max():>10.2f}")


if __name__ == "__main__":
    main()
