#!/usr/bin/env python3
"""Presentation cut of the r-t front figure. Same fits, no words on it.

This is classic_front_fit.py's figure redrawn for a slide. Two differences,
both requested:

  1. NO LEGEND AND NO IN-PANEL TEXT. Not a legend, not curve labels, not the
     event count, not the shut-in markers' names -- every word is added in
     PowerPoint, where it can be placed and sized against the rest of the
     slide. Only axis and tick labels remain, and the graphical marks that
     cannot be redrawn in a slide editor over the right x values: the dashed
     verticals at injection resumption (4.300 d) and shut-in (17.455 d), and
     the shaded first shut-in, 1.582-3.558 d.

     THE SHADED BAND ENDS AT 3.558 d, NOT 4.300. The rate is zero only to
     3.558 d; from there to 4.298 d a 1.6-2.5 L/s trickle runs, and 4.300 d is
     where it steps to ~8 L/s and climbs. Earlier versions shaded through to
     4.300 and so called 0.74 d of low-rate injection a shut-in. The interval
     now comes from front_backfront.zero_rate_interval(), read off the rate
     file. 4.300 d remains the cycle-2 origin and its dashed marker stays --
     sustained injection does resume there -- it is just not where the shut-in
     ended.

     No grid either, and the dashed verticals are drawn in a grey that no data
     series uses. They were the same INK as the Q line, which made a time
     marker look like part of the rate history.

     What this costs: the figure no longer states what its own curves are, so
     everything below has to be said out loud. main() prints all of it on every
     run, which is the intended source for the slide's text.

  2. LOWER PANEL REDRAWN AFTER WANG & DUNHAM figure 2a: injection rate as a
     LINE on the left axis and the wellhead pressure CHANGE on a twin right
     axis, rather than rate alone as a filled area. The axis labels and ticks
     are coloured to match their curves -- with no legend and no annotation,
     colour is the only thing left to identify them by.

WHAT IS DRAWN, since nothing on the figure says so. Grey dots: all 20 734
located events of the November 2012 stimulation, distance from the injector
against time. Red: the triggering front r = sqrt(4 pi D t), t measured from
injection resumption, D fitted so 95% of events fall below. Blue: Parotidis's
back front, t from first injection and t_s the injection duration, D fitted so
5% of post-shut-in events fall below.

THE FITS ARE NOT REDONE HERE. D_trig and D_back come from calling
classic_front_fit's own fit_quantile/trig_front/back_front, so this figure and
the notebook's cannot drift apart; if the fit changes there it changes here.

dp CONVENTION. dp is measured from the first sample of the wellhead record,
34.412 MPa, which is Wang & Dunham's convention (see
front_backfront.pressure_history). Two consequences leave the axis, and with
the annotations gone NEITHER IS MARKED:

  - during the first shut-in the gauge bleeds to atmospheric, so dp falls to
    -35.07 MPa and leaves the bottom of the axis. The shaded band is the only
    remaining indication. The axis is clipped at -2 as in their figure 2a
    rather than rescaled to hold the gauge emptying, which is not the
    reservoir.
  - a burst at 14.2 d reaches 87.76 MPa absolute, dp = +53.34, which is 92% of
    the vertical stress and almost certainly water hammer. 1036 samples clear
    +25 over 78 minutes, only 2 above +50; away from 14.20-14.38 d nothing in
    the record exceeds +20, which is what sets the axis at +22. The green trace
    simply runs off the top there.

Usage:  python slide_front_fit.py
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import importlib.util as iu


def _load(name):
    s = iu.spec_from_file_location(name, str(Path(__file__).with_name(name + ".py")))
    m = iu.module_from_spec(s)
    s.loader.exec_module(m)
    return m


fb = _load("front_backfront")
cf = _load("classic_front_fit")

T_ON, T_SHUT, TS = cf.T_ON, cf.T_SHUT, cf.TS
T_RESUME = cf.T_RESUME
T_SHUT1, T_TRICKLE = fb.zero_rate_interval()     # 1.582 -> 3.558 d
OUT = Path("/home/users/nberrios/3dhbi/hbi_analysis/figures/postshutin")
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
RED, BLUE, GRN = "#a8071a", "#1d4ed8", "#009E73"
EVERY = 30                                       # pressure decimation, ~35 s
# The time markers get a colour NO DATA SERIES USES. They were INK, which is
# also the Q line, so a dashed vertical read as part of the rate history.
MARK = "#8E44AD"

plt.rcParams.update({"font.size": 12.5, "axes.titlesize": 14,
                     "axes.labelsize": 13.5, "axes.edgecolor": MUTED,
                     "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED,
                     # tick NUMBERS black on every axis; xtick.color/ytick.color
                     # above now govern only the tick marks themselves
                     "xtick.labelcolor": INK, "ytick.labelcolor": INK,
                     "xtick.labelsize": 12, "ytick.labelsize": 12})


def main(argv=None):
    """argv accepted and ignored, matching the other two scripts here."""
    t_abs, r, n_raw = fb.catalogue()
    ti, q = fb.rate_history()
    tp, p_abs, dp = fb.pressure_history(every=EVERY)
    tp_full, p_full, dp_full = fb.pressure_history()

    tt = t_abs - T_RESUME
    kt = tt > 0
    D_trig = cf.fit_quantile(tt[kt], r[kt], cf.trig_front, 0.95)
    tb = t_abs - T_ON
    post = tb > TS
    D_back = cf.fit_quantile(tb[post], r[post], cf.back_front, 0.05)

    i_pk = int(np.argmax(dp_full))
    print(f"{len(t_abs)} events of {n_raw}; D_trig = {D_trig:.4f}, "
          f"D_back = {D_back:.4f} m^2/s  (fits from classic_front_fit)")
    print(f"dp measured from p[0] = {p_full[0]:.3f} MPa (Wang & Dunham's "
          f"convention); dp spans {dp_full.min():+.2f} to {dp_full.max():+.2f}")
    print(f"  peak dp {dp_full.max():+.2f} MPa (p = {p_full[i_pk]:.2f}) at "
          f"{tp_full[i_pk]:.2f} d -- {(dp_full > 25).sum()} samples above +25, "
          f"off the plotted axis")
    print(f"  minimum dp {dp_full.min():+.2f} MPa, gauge bled off during the "
          f"first shut-in ({T_SHUT1:.3f}-{T_TRICKLE:.3f} d)")
    print(f"  pressure decimated {EVERY}x for drawing "
          f"({len(dp_full)} -> {len(dp)} samples); extrema above are undecimated")

    fig, (ax, axq) = plt.subplots(
        2, 1, figsize=(12.0, 8.6), dpi=200, sharex=True,
        gridspec_kw=dict(height_ratios=[2.9, 1.45], hspace=0.08))

    # ------------------------------------------------------------ upper: r-t
    ax.scatter(t_abs, r, s=4.0, alpha=0.22, color=MUTED, lw=0)
    g = np.linspace(T_RESUME + 1e-3, t_abs.max(), 800)
    rg = cf.trig_front(g - T_RESUME, D_trig)
    ax.plot(g, rg, "-", lw=3.2, color=RED)
    gb = np.linspace(T_SHUT + 1e-3, t_abs.max(), 800)
    rb = cf.back_front(gb - T_ON, D_back)
    ax.plot(gb, rb, "-", lw=3.2, color=BLUE)

    # NO IN-PANEL TEXT AT ALL -- every word goes on in PowerPoint. The
    # graphical marks stay, because they are not words and cannot be redrawn in
    # a slide editor over the right x values: the two dashed verticals at
    # injection resumption and shut-in, and the shaded first shut-in.
    for x in (T_RESUME, T_SHUT):
        ax.axvline(x, color=MARK, lw=1.7, ls="--")
    ax.axvspan(T_SHUT1, T_TRICKLE, color=GRID, alpha=0.65, zorder=0)
    ax.set(ylabel="Distance from injection point (m)",
           xlim=(0, t_abs.max()), ylim=(0, 1750))

    # ------------------------------------------- lower: Q and dp, Taiyi 2a
    axq.plot(ti, q, "-", lw=2.0, color=INK, zorder=4)
    axq.set(xlabel="Days since 2012-11-13", xlim=(0, t_abs.max()),
            ylim=(0, 70))
    axq.set_ylabel("Q (L/s)", color=INK)
    axq.tick_params(axis="y", color=MUTED, labelcolor=INK)
    axq.set_axisbelow(True)

    axp = axq.twinx()
    axp.plot(tp, dp, "-", lw=1.5, color=GRN, alpha=0.9, zorder=3)
    # Spelled out, not "$\Delta p$ at wellhead" -- on a slide the axis has to
    # name the quantity without the audience decoding a symbol first.
    # Two lines: spelled out in full it is longer than the 2.9 in panel is
    # tall, so on one line it overran into the panel above.
    axp.set_ylabel("Wellhead pressure\nchange (MPa)", color=GRN,
                   linespacing=1.5)
    axp.tick_params(axis="y", color=GRN, labelcolor=INK)
    axp.set_ylim(-2, 22)
    axp.spines["right"].set_color(GRN)
    axp.spines["top"].set_visible(False)
    # Both off-axis features are now UNLABELLED, so they need saying out loud
    # when the slide goes up. The gauge bleeds to atmospheric during the shaded
    # interval, dp reaching -35.07 MPa, and the burst at 14.2 d leaves the top:
    # 1036 samples clear +25 over 78 minutes, the tallest at 87.8 MPa absolute,
    # only 2 above +50. Away from 14.20-14.38 d nothing exceeds +20, which is
    # what sets the axis at +22. main() prints all of it every run.

    for a in (axq, axp):
        a.axvspan(T_SHUT1, T_TRICKLE, color=GRID, alpha=0.65, zorder=0)
        for x in (T_RESUME, T_SHUT):
            a.axvline(x, color=MARK, lw=1.7, ls="--", zorder=1)

    OUT.mkdir(parents=True, exist_ok=True)
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"slide_front_fit.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {OUT}/slide_front_fit.png and .pdf")
    return D_trig, D_back


if __name__ == "__main__":
    main()
