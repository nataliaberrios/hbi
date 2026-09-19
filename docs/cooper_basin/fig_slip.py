#!/usr/bin/env python3
"""PAPER FIGURE: cumulative slip, in the fig6_remake_strike_notitle format.

Stacked full-width panels, one per dataset, each showing the whole time family
in viridis with a two-column legend, boxed axes, both axes labelled on every
panel, no titles. That is make_fig6.py's layout and constants -- FS_LABEL 13,
FS_LEGEND 9, lw 1.9, xlim +-1.5 km, shared y-limit -- reused here rather than
inventing a second slip-profile style.

    top     observed, from the seismicity catalogue
    bottom  the reference simulation

THE OBSERVED PANEL IS AN INCREMENT SINCE t0, WHICH make_fig6.py's IS NOT, and
that is the one deliberate departure. main_LH.f90:923 sets slip = 0, so a
cycle-2 run starts with no slip history while the real fault had already
accumulated slip by data-day 4.300. Plotting the absolute observed profile
against a model that starts from zero compares two different quantities and
makes the model look about twice as bad as it is. The increment is formed
PER-X -- the observed profile is interpolated in time at t0 at every position
and subtracted -- because the pre-t0 slip is not spatially uniform and removing
a single scalar would distort the shape. Pass --absolute for the untouched
profiles.

THE TIMES ARE DATA-DAYS IN BOTH PANELS and the colours line up, so one colour
is one instant in the field. The simulation is sampled at td - t0_of(deck).
Data-day 3 is dropped because it precedes t0.

BOTH PANELS SHARE A Y-LIMIT. Autoscaling each would put the two families at the
same apparent height and hide the amplitude shortfall, which is the entire
content of the comparison.

WHAT IT SHOWS. The model is too LOW, not too NARROW: its profiles reach roughly
the right distance at about half the amplitude. The observed lobe is visibly
asymmetric, peaking near x = -96 m, while the model is symmetric about the
injector -- which is why a scalar peak ratio is slippery here and runs from
0.35x to 0.49x depending on whether peaks or the injector are compared.

THE CATALOGUE SLIP IS AN INVERSION, NOT A MEASUREMENT: cumulative slip from the
seismicity under a circular crack model, scaling as sum(M_0)/(mu pi r^2). It is
sensitive to the assumed rupture radius and is a LOWER bound on total slip,
since the catalogue holds only what radiated. The ragged curves are that
inversion, not gauge noise.

Usage:
  python fig_slip.py                     # 633001 against the catalogue
  python fig_slip.py --runs 633001 633110
  python fig_slip.py --absolute
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = Path("/home/users/nberrios/3dhbi/hbi_analysis")
sys.path.insert(0, str(H / "notebooks"))
sys.path.insert(0, str(H))
from sim_curves import load_slip, _deck, _ff
from injection_t0 import t0_of

FIG = H / "figures" / "cycle2"
OBS = Path("/home/users/nberrios/3dhbi/hbi/slip_profiles_strike.txt")
OT = np.array([3, 5, 7, 9, 11, 13, 15, 17], float)   # data-days in the file
XLIM = 1.5                                           # km, as make_fig6.py
FS_LABEL, FS_LEGEND = 13, 9


def observed():
    o = np.loadtxt(OBS)
    return o[:, 0], o[:, 1:]        # x in km, profiles in cm


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", type=int, default=[633001])
    ap.add_argument("--absolute", action="store_true",
                    help="observed absolute, not the increment since t0")
    ap.add_argument("--out", default="fig_slip")
    a = ap.parse_args(argv)

    x_obs, prof = observed()
    t0 = t0_of(_deck(a.runs[0]))
    times = [td for td in OT if td > t0 + 0.2]
    dropped = [td for td in OT if td not in times]
    if dropped:
        print(f"data-days before t0 = {t0:g}, dropped: "
              f"{[f'{d:.0f}' for d in dropped]}")

    at0 = np.array([np.interp(t0, OT, prof[i, :]) for i in range(len(x_obs))])
    Yobs = np.column_stack([
        np.array([np.interp(td, OT, prof[i, :]) for i in range(len(x_obs))])
        - (0.0 if a.absolute else at0) for td in times])
    if not a.absolute:
        Yobs = np.maximum(Yobs, 0.0)

    panels = [(x_obs, Yobs)]
    for n in a.runs:
        tn = t0_of(_deck(n))
        stack, xs = [], None
        for td in times:
            r, sl, ta = load_slip(n, td - tn, how="strike")
            rk = np.asarray(r) / 1000.0                  # metres -> km
            xs = np.concatenate([-rk[::-1], rk[1:]])     # mirror the model
            stack.append(np.concatenate([sl[::-1], sl[1:]]) * 100.0)
        panels.append((xs, np.column_stack(stack)))

    cmap = plt.get_cmap("viridis")
    cols = [cmap(v) for v in np.linspace(0.0, 1.0, len(times))]
    nrow = len(panels)
    fig, axes = plt.subplots(nrow, 1, figsize=(9.0, 3.8 * nrow), dpi=300,
                             constrained_layout=True)
    axes = np.atleast_1d(axes)
    fig.get_layout_engine().set(h_pad=0.18, hspace=0.10)

    for ax, (x, Y) in zip(axes, panels):
        for i, td in enumerate(times):
            ax.plot(x, Y[:, i], lw=1.9, color=cols[i], label=f"{td:.0f} days")
        ax.set_xlabel("Distance along-strike (km)", fontsize=FS_LABEL)
        ax.set_ylabel("Cumulative slip (cm)", fontsize=FS_LABEL)
        ax.set_xlim(-XLIM, XLIM)
        ax.tick_params(labelsize=FS_LABEL - 2)
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=FS_LEGEND, ncol=2, loc="upper right", frameon=True)

    ymax = max(ax.get_ylim()[1] for ax in axes)
    for ax in axes:
        ax.set_ylim(0, ymax)

    for e in ("png", "pdf"):
        fig.savefig(FIG / f"{a.out}.{e}")
    plt.close(fig)

    kind = "ABSOLUTE" if a.absolute else f"increment since data-day {t0:g}"
    print(f"observed panel: {kind}")
    print(f"{'data-day':>9} {'obs peak':>9} " +
          " ".join(f"{n:>9}" for n in a.runs) + f" {'ratio':>7}")
    for i, td in enumerate(times):
        row = f"{td:>9.0f} {Yobs[:, i].max():>9.3f}"
        for j in range(len(a.runs)):
            row += f" {panels[j + 1][1][:, i].max():>9.3f}"
        row += f" {panels[1][1][:, i].max() / Yobs[:, i].max():>7.2f}"
        print(row)
    print(f"wrote {FIG}/{a.out}.png")


if __name__ == "__main__":
    main()
