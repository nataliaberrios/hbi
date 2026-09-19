#!/usr/bin/env python3
"""PAPER FIGURE: cumulative slip, modelled against the catalogue.

(a) slip against distance along strike, at three times
(b) slip at the injector against time

This is the figure section 8 rests on, and the set did not have it: Figs. 4 and
5 carry slip only as a ratio, and the one model-vs-observed slip comparison
that existed for the reference run (obsoverlay_633001) is a diagnostic drawn in
the per-run house style.

THE OBSERVED CURVE IS AN INCREMENT SINCE t0. main_LH.f90:923 sets slip = 0
unconditionally, so the model starts with no slip history while the real fault
had already accumulated slip by data-day 4.300. The increment is formed per-x
-- the observed profile is interpolated in time at t0 at EVERY position and
subtracted -- because the pre-t0 slip is not spatially uniform and subtracting
a single scalar would distort the shape.

AND THE SCALAR RATIO DEPENDS ON A CHOICE WORTH 40%. The observed profile peaks
at x = -96 m, not at the well, while the model peaks at x = 0 by symmetry, so
"peak" and "at the injector" are not the same comparison. At data-day 13:

    peak of the per-x difference     3.061 cm   ->  0.35x
    difference of the peaks          2.598 cm   ->  0.42x
    at the injector, x = 0           2.214 cm   ->  0.49x

The middle one is what score_cycle2.py uses and it is the least defensible of
the three: it subtracts the peak at t0 from the peak at t when those sit at
different positions, which is not an operation on a field. Panel (b) uses the
last, both curves at x = 0, and panel (a) shows the profiles so the offset peak
is visible rather than buried in a scalar. Comparing against the ABSOLUTE
observed peak, as the per-run overlay does, gives 0.20x.

BOTH SIDES ARE PLOTTED AS THEY ARE. The observed lobe is genuinely asymmetric
-- it runs about -400 to +80 m, so its front is on the negative side -- while
the model is symmetric about the injector. Folding either one to compare
half-widths would hide that, so the model is mirrored and the observation is
left alone, and the asymmetry is visible.

WHAT THE FIGURE SHOWS. The model reaches 0.49x of the observed increment at
the injector at 8.7 d, and the shortfall is not uniform in shape: panel (a) shows the modelled
profile is too low everywhere rather than too narrow, and panel (b) shows why
-- the model slips in two bursts, at disc breakout near 3.5 d and at the day-9
rate step-up, and is flat between 4 and 9 d while the observation climbs
steadily. A model that under-predicts slip by 2x to 2.9x, depending on the
definition, while matching the seismicity front and the wellhead pressure is
the open question of section 8,
and this is where a reader should be able to see it.

THE CATALOGUE SLIP IS AN INVERSION, NOT A MEASUREMENT -- cumulative slip 20 m
west of the well, from the seismicity catalogue under a circular crack model
(see make_fig6.py). It scales as sum(M_0)/(mu pi r^2), so it is sensitive to
the assumed rupture radius, and it is a lower bound on TOTAL slip because the
catalogue contains only what radiated. Both belong in the caption: the quantity
being under-predicted is itself uncertain, and is the lone dissenter among the
three observables.

Usage:  python fig_slip.py [--run 633001]
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
OT = np.array([3, 5, 7, 9, 11, 13, 15, 17], float)   # data-days
SHOW = (9.0, 13.0, 17.0)                             # data-days to draw in (a)
T_EVAL = 8.7                                         # sim days, the scoring time
INK = "#141413"

plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 9, "legend.fontsize": 7.8,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "axes.edgecolor": INK, "axes.linewidth": 0.7,
    "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": INK, "ytick.color": INK,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 3.0, "ytick.major.size": 3.0,
    "axes.grid": False, "legend.frameon": False,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def panel(ax, letter):
    ax.text(-0.17, 1.02, f"({letter})", transform=ax.transAxes,
            fontsize=10, fontweight="bold", va="bottom", ha="left")


def observed():
    """(x_m, profiles_cm[nx, nt]) of absolute cumulative slip."""
    o = np.loadtxt(OBS)
    return o[:, 0] * 1000.0, o[:, 1:]


def obs_increment(xo, prof, td, t0):
    """Observed slip at data-day td MINUS its value at t0, per position.

    Per-x rather than a single scalar: the slip already on the fault at t0 is
    not spatially uniform, so subtracting its peak would distort the shape.
    """
    at_t = np.array([np.interp(td, OT, prof[i, :]) for i in range(len(xo))])
    at_0 = np.array([np.interp(t0, OT, prof[i, :]) for i in range(len(xo))])
    return np.maximum(at_t - at_0, 0.0)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", type=int, default=633001)
    ap.add_argument("--out", default="fig_slip")
    a = ap.parse_args(argv)

    d = _deck(a.run)
    t0 = t0_of(d)
    t_end = _ff(d["tmax"]) * 365.0
    xo, prof = observed()

    fig, (a0, a1) = plt.subplots(1, 2, figsize=(7.48, 3.25), dpi=400,
                                 constrained_layout=True)
    cols = plt.cm.viridis(np.linspace(0.08, 0.78, len(SHOW)))
    tlab = []

    print(f"run {a.run}, t0 = data-day {t0:g}, tmax {t_end:.2f} d")
    print(f"{'data-day':>9} {'sim d':>7} {'obs peak':>9} {'model':>8} {'ratio':>7}")
    for c, td in zip(cols, SHOW):
        ts = td - t0
        if ts <= 0 or ts > t_end + 0.3:
            print(f"{td:>9.1f}  out of range"); continue
        oi = obs_increment(xo, prof, td, t0)
        a0.plot(xo, oi, ls=(0, (3, 2)), lw=1.6, color=c)
        r, sl, ta = load_slip(a.run, ts, how="strike")
        rm = np.asarray(r)
        xs = np.concatenate([-rm[::-1], rm[1:]])          # mirror the model
        ss = np.concatenate([sl[::-1], sl[1:]]) * 100.0
        a0.plot(xs, ss, "-", lw=1.7, color=c)
        tlab.append((c, f"data-day {td:.0f}  (sim {ts:.1f} d)"))
        print(f"{td:>9.1f} {ta:>7.2f} {oi.max():>9.3f} {ss.max():>8.3f} "
              f"{ss.max()/oi.max():>7.2f}")

    # ONE legend. Two overlapping legends collided in the upper-left, where
    # the data also lives.
    hs = [plt.Line2D([], [], color=c, lw=1.7) for c, _ in tlab] + \
         [plt.Line2D([], [], color="none"),
          plt.Line2D([], [], color=INK, lw=1.7, ls="-"),
          plt.Line2D([], [], color=INK, lw=1.6, ls=(0, (3, 2)))]
    ls_ = [t for _, t in tlab] + ["", "model", "catalogue increment"]
    a0.legend(hs, ls_, loc="upper right", handlelength=2.0,
              labelspacing=0.28, borderaxespad=0.2)
    a0.set(xlabel="distance along strike from the injector (m)",
           ylabel="cumulative slip (cm)", xlim=(-1000, 1000))

    # --- (b) at the injector, against time
    ts_, ss_ = [], []
    for t_ in np.linspace(0.1, t_end, 60):
        try:
            _, sl, ta = load_slip(a.run, t_, how="strike")
        except Exception:
            continue
        if abs(ta - t_) < 0.3:
            ts_.append(ta); ss_.append(sl[0] * 100.0)
    ts_, ss_ = np.array(ts_), np.array(ss_)
    # AT THE INJECTOR FOR BOTH CURVES, x = 0. The observed profile peaks at
    # x = -96 m, not at the well, so a "peak" comparison is not like-for-like
    # while the model's peak is at x = 0 by symmetry. Three definitions of the
    # observed increment are defensible and they span 40%:
    #   peak of the per-x difference        3.061 cm   -> 0.35x
    #   difference of the peaks (scorer)    2.598 cm   -> 0.42x
    #   at the injector                     2.214 cm   -> 0.49x
    # The middle one subtracts peaks measured at different locations, which is
    # not a field operation; this panel uses the last, and panel (a) shows the
    # profiles so the offset peak is visible rather than hidden in a scalar.
    tdg = np.linspace(t0, OT.max(), 200)
    at_well = np.array([np.interp(0.0, xo, prof[:, j])
                        for j in range(prof.shape[1])])      # obs at x = 0
    so = np.interp(tdg, OT, at_well) - float(np.interp(t0, OT, at_well))
    a1.plot(tdg - t0, so, ls=(0, (3, 2)), lw=1.8, color=INK,
            label="catalogue increment")
    a1.plot(ts_, ss_, "-", lw=1.8, color="#D55E00", label=f"{a.run}")
    rat = float(np.interp(T_EVAL, ts_, ss_)) / float(np.interp(T_EVAL, tdg - t0, so))
    a1.set(xlabel=f"days since injection resumed (data-day {t0:g})",
           ylabel="cumulative slip at the injector (cm)",
           xlim=(0, max(t_end, tdg.max() - t0)), ylim=(0, None))
    a1.legend(loc="upper left", handlelength=2.0, labelspacing=0.3)

    for x, l in ((a0, "a"), (a1, "b")):
        panel(x, l)
    for e in ("png", "pdf"):
        fig.savefig(FIG / f"{a.out}.{e}")
    plt.close(fig)
    print(f"\nat sim {T_EVAL} d the model is {rat:.2f}x the catalogue increment")
    print(f"wrote {FIG}/{a.out}.png")


if __name__ == "__main__":
    main()
