#!/usr/bin/env python3
"""Presentation cut of the r-t front figure. Same fits, no words on it.

This is classic_front_fit.py's figure redrawn for a slide. Two differences,
both requested:

  0. TWO CUTS ARE WRITTEN. slide_front_fit.* is the plain figure;
     slide_front_fit_M3.* adds magenta dots on the four Mw >= 3 events, same
     marker as the rest of the catalogue and 2.5x the diameter. Both come from
     one pass over the data, so they cannot diverge.

  1. NO LEGEND AND NO IN-PANEL TEXT. Not a legend, not curve labels, not the
     event count, not the shut-in markers' names -- every word is added in
     PowerPoint, where it can be placed and sized against the rest of the
     slide. Only axis and tick labels remain, plus the one graphical mark that
     cannot be redrawn in a slide editor over the right x values.

  1b. ONE MARK, ONE MEANING. A purple dashed vertical wherever injection
     stops, and nothing else -- no end-of-shut-in, no resumption marker at
     4.300 d, no seismicity onset, no shaded band, no grid. Earlier versions
     mixed all of these into the same dashed-line style, so the reader had to
     be told case by case what each line meant.

     front_backfront.shutin_starts() reads them off the rate file: 1.582 and
     17.151 d. It thresholds on duration, because the record has 27 zero-rate
     runs and 24 of them are under half an hour -- dropouts, not decisions. The
     cut at 0.1 d falls in a wide gap; the longest run it excludes is 0.035 d.

     THE FINAL SHUT-IN IS AT 17.151 d, NOT T_SHUT = 17.455. Injection stops at
     17.095, restarts for 31 minutes, and stops for good at 17.151; 17.455 is
     the rate file's LAST SAMPLE. So the blue curve, which starts at T_SHUT
     because that is the t_s the fit uses, begins 0.30 d to the right of the
     purple line it belongs to. That gap is real and visible, not a drawing
     error, and re-fitting on t_s = 16.650 d rather than 16.954 d would close
     it.

     The purple is #8E44AD. The marks were INK, the same colour as the Q line,
     which made a time marker look like part of the rate history.

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
    -35.07 MPa and leaves the bottom of the axis. With the band gone there is
    now NO indication of this at all; the green trace simply drops out between
    1.582 and 3.558 d. The axis is clipped at -2 as in their figure 2a rather
    than rescaled to hold the gauge emptying, which is not the reservoir.
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
SHUTINS = fb.shutin_starts()                     # [1.582, 17.151] d
OUT = Path("/home/users/nberrios/3dhbi/hbi_analysis/figures/postshutin")
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
RED, BLUE, GRN = "#a8071a", "#1d4ed8", "#009E73"
EVERY = 30                                       # pressure decimation, ~35 s
# The time markers get a colour NO DATA SERIES USES. They were INK, which is
# also the Q line, so a dashed vertical read as part of the rate history.
MARK = "#8E44AD"
MAG = "#e6007e"                                  # the Mw >= MBIG dots
MBIG = 3.0

# "Days since injection began" would be WRONG on the current axis: injection
# starts at T_ON = 0.501 d, so t = 0 here is half a day BEFORE it. Flipping
# this to True subtracts T_ON from every plotted time and relabels, which is
# self-consistent but renumbers the figure against every other one in the
# project -- the shut-ins become 1.081 and 16.650 d, resumption 3.799 d.
X_FROM_INJECTION = False

plt.rcParams.update({"font.size": 12.5, "axes.titlesize": 14,
                     "axes.labelsize": 13.5, "axes.edgecolor": MUTED,
                     "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED,
                     # tick NUMBERS black on every axis; xtick.color/ytick.color
                     # above now govern only the tick marks themselves
                     "xtick.labelcolor": INK, "ytick.labelcolor": INK,
                     "xtick.labelsize": 12, "ytick.labelsize": 12})


def main(argv=None):
    """Writes BOTH cuts: slide_front_fit.* and slide_front_fit_M3.*.

    argv accepted and ignored, matching the other two scripts here. The starred
    version is a variant to choose between, not a replacement, so both are
    built from one pass over the data.
    """
    t_abs, r, ml, mw = fb.catalogue_mag()
    _, _, n_raw = fb.catalogue()
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
          f"first shut-in")
    print(f"  shut-ins marked at " +
          ", ".join(f"{x:.3f}" for x in SHUTINS) + " d; the blue curve starts "
          f"at T_SHUT = {T_SHUT} d, the rate file's last sample, "
          f"{T_SHUT - SHUTINS[-1]:.3f} d later")
    print(f"  pressure decimated {EVERY}x for drawing "
          f"({len(dp_full)} -> {len(dp)} samples); extrema above are undecimated")
    hi = mw >= MBIG
    print(f"\nM3 version: {int(hi.sum())} events with Mw >= {MBIG} "
          f"(ML >= {MBIG} would give {int((ml >= MBIG).sum())}, "
          f"{int(((mw >= MBIG) & (ml >= MBIG)).sum())} in both)")
    idx = np.where(hi)[0]
    for j in idx[np.argsort(-mw[idx])]:
        print(f"    Mw {mw[j]:.2f}  ML {ml[j]:.2f}  t = {t_abs[j]:7.3f} d  "
              f"r = {r[j]:6.0f} m")

    for big in (False, True):
        _draw(big, t_abs, r, mw, ti, q, tp, dp, D_trig, D_back)
    return D_trig, D_back


def _draw(big, t_abs, r, mw, ti, q, tp, dp, D_trig, D_back):
    """One figure. `big` adds the Mw >= MBIG overlay and changes the stem."""
    off = T_ON if X_FROM_INJECTION else 0.0
    xlab = ("Days since injection began" if X_FROM_INJECTION
            else "Days since 2012-11-13")
    fig, (ax, axq) = plt.subplots(
        2, 1, figsize=(12.0, 8.6), dpi=200, sharex=True,
        gridspec_kw=dict(height_ratios=[2.9, 1.45], hspace=0.08))

    # ------------------------------------------------------------ upper: r-t
    ax.scatter(t_abs - off, r, s=4.0, alpha=0.22, color=MUTED, lw=0)
    g = np.linspace(T_RESUME + 1e-3, t_abs.max(), 800)
    rg = cf.trig_front(g - T_RESUME, D_trig)
    ax.plot(g - off, rg, "-", lw=3.2, color=RED)
    gb = np.linspace(T_SHUT + 1e-3, t_abs.max(), 800)
    rb = cf.back_front(gb - T_ON, D_back)
    ax.plot(gb - off, rb, "-", lw=3.2, color=BLUE)

    # ONE MARK, ONE MEANING: a purple dashed vertical wherever injection stops.
    # Nothing else is marked -- not the end of a shut-in, not the resumption at
    # 4.300 d, not the onset of seismicity -- and the shaded band is gone with
    # them, since its right edge marked the end of the first shut-in.
    for x in SHUTINS:
        ax.axvline(x - off, color=MARK, lw=1.7, ls="--")

    # The one optional overlay: same marker as every other event, magenta and
    # a little bigger. s = 4 -> 26 is 2.5x the DIAMETER, which is as small as
    # it can be and still be found on a projected slide. Mw, not ML -- see
    # fb.catalogue_mag(); at the top of the catalogue the two scales pick
    # different events.
    if big:
        k = mw >= MBIG
        ax.scatter(t_abs[k] - off, r[k], s=26, c=MAG, lw=0, zorder=6)

    ax.set(ylabel="Distance from injection point (m)",
           xlim=(0, t_abs.max() - off), ylim=(0, 1750))

    # ------------------------------------------- lower: Q and dp, Taiyi 2a
    axq.plot(ti - off, q, "-", lw=2.0, color=INK, zorder=4)
    axq.set(xlabel=xlab, xlim=(0, t_abs.max() - off), ylim=(0, 70))
    axq.set_ylabel("Q (L/s)", color=INK)
    axq.tick_params(axis="y", color=MUTED, labelcolor=INK)
    axq.set_axisbelow(True)

    axp = axq.twinx()
    axp.plot(tp - off, dp, "-", lw=1.5, color=GRN, alpha=0.9, zorder=3)
    # Spelled out, not "$\Delta p$ at wellhead" -- on a slide the axis has to
    # name the quantity without the audience decoding a symbol first.
    # Two lines: spelled out in full it is longer than the 2.9 in panel is
    # tall, so on one line it overran into the panel above.
    axp.set_ylabel("Wellhead pressure\nchange (MPa)", color=GRN,
                   linespacing=1.5)
    # Right-hand numbers GREEN, matching their own label and curve. The left
    # axes stay black. This is the one place colour earns its keep: two
    # quantities share the panel, so the numbers have to say which scale they
    # belong to.
    axp.tick_params(axis="y", color=GRN, labelcolor=GRN)
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
        for x in SHUTINS:
            a.axvline(x - off, color=MARK, lw=1.7, ls="--", zorder=1)

    stem = "slide_front_fit_M3" if big else "slide_front_fit"
    OUT.mkdir(parents=True, exist_ok=True)
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"{stem}.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/{stem}.png and .pdf")


if __name__ == "__main__":
    main()
