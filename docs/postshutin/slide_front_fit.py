#!/usr/bin/env python3
"""Presentation cut of the r-t front figure. Same fits, no legend.

This is classic_front_fit.py's figure redrawn for a slide. Two differences,
both requested:

  1. NO LEGEND. On a slide the speaker says what the curves are, and the box
     sat over the empty upper-left corner that the eye wants for the front's
     early rise. The two curves are instead labelled INLINE at their right-hand
     ends, where the axes are empty, each in its own colour.

  2. LOWER PANEL REDRAWN AFTER WANG & DUNHAM figure 2a: injection rate as a
     LINE on the left axis and the wellhead pressure CHANGE on a twin right
     axis, rather than rate alone as a filled area. Axis labels are coloured to
     match their curves, which is what carries the identification with no
     legend present.

THE FITS ARE NOT REDONE HERE. D_trig and D_back come from calling
classic_front_fit's own fit_quantile/trig_front/back_front, so this figure and
the notebook's cannot drift apart; if the fit changes there it changes here.

dp CONVENTION. dp is measured from the first sample of the wellhead record,
34.412 MPa, which is Wang & Dunham's convention (see
front_backfront.pressure_history). Consequences visible on the panel:

  - during the first shut-in the gauge bleeds to atmospheric, so dp falls to
    -35 MPa and leaves the bottom of the axis. The interval is shaded, and the
    axis is clipped at -2 as in their figure 2a rather than rescaled to hold an
    excursion that is the gauge emptying, not the reservoir.
  - a two-sample transient at 14.23 d reaches 87.76 MPa absolute, dp = +53.3,
    which is 92% of the vertical stress and almost certainly water hammer. It
    is off the top of the axis and annotated as such, not deleted.

Usage:  python slide_front_fit.py
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

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
T_SHUT1 = 1.582                                  # first shut-in begins
OUT = Path("/home/users/nberrios/3dhbi/hbi_analysis/figures/postshutin")
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
RED, BLUE, GRN = "#a8071a", "#1d4ed8", "#009E73"
EVERY = 30                                       # pressure decimation, ~35 s
HALO = [pe.withStroke(linewidth=3.4, foreground="white")]

plt.rcParams.update({"font.size": 12.5, "axes.titlesize": 14,
                     "axes.labelsize": 13.5, "axes.edgecolor": MUTED,
                     "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED,
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
          f"first shut-in ({T_SHUT1}-{T_RESUME} d)")
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

    # Inline curve labels replace the legend, each in the empty region its own
    # curve bounds: the triggering front's is up and to its left, where no
    # event can be because none has arrived yet, and the back front's is inside
    # the expanding quiet wedge, which is the feature it exists to mark. The
    # formulae are deliberately not on the figure -- they are in the caption and
    # in what the speaker says, and on a projector they only cost legibility.
    ax.text(5.55, 760, "triggering front\n" f"$D$ = {D_trig:.3f} m$^2$/s",
            ha="left", va="bottom", color=RED, fontsize=13.5,
            linespacing=1.45, path_effects=HALO, zorder=7)
    ax.text(18.15, 75, "back front\n" f"$D$ = {D_back:.3f} m$^2$/s",
            ha="left", va="bottom", color=BLUE, fontsize=13.5,
            linespacing=1.45, path_effects=HALO, zorder=7)
    ax.text(0.35, 1690, f"{len(t_abs):,} located events".replace(",", " "),
            ha="left", va="top", color=MUTED, fontsize=12,
            path_effects=HALO, zorder=7)

    for x, lab in ((T_RESUME, "injection resumes"), (T_SHUT, "shut-in, $t_s$")):
        ax.axvline(x, color=INK, lw=1.7, ls="--")
        ax.text(x + 0.17, 900, lab, rotation=90, ha="left", va="top",
                fontsize=11.5, color=INK, path_effects=HALO, zorder=6)
    ax.axvspan(T_SHUT1, T_RESUME, color=GRID, alpha=0.65, zorder=0)
    ax.text(2.94, 900, "1st shut-in", rotation=90, ha="center", va="top",
            fontsize=11.5, color=MUTED, path_effects=HALO, zorder=6)
    ax.set(ylabel="Distance from injection point (m)",
           xlim=(0, t_abs.max()), ylim=(0, 1750))
    ax.grid(alpha=0.28, color=GRID)

    # ------------------------------------------- lower: Q and dp, Taiyi 2a
    axq.plot(ti, q, "-", lw=2.0, color=INK, zorder=4)
    axq.set(xlabel="Days since 2012-11-13", xlim=(0, t_abs.max()),
            ylim=(0, 70))
    axq.set_ylabel("Q (L/s)", color=INK)
    axq.tick_params(axis="y", colors=INK)
    axq.grid(alpha=0.28, color=GRID)
    axq.set_axisbelow(True)

    axp = axq.twinx()
    axp.plot(tp, dp, "-", lw=1.5, color=GRN, alpha=0.9, zorder=3)
    axp.set_ylabel(r"$\Delta p$ at wellhead (MPa)", color=GRN)
    axp.tick_params(axis="y", colors=GRN)
    axp.set_ylim(-2, 22)
    axp.spines["right"].set_color(GRN)
    axp.spines["top"].set_visible(False)
    # The excursion is a 78-minute BURST of spikes, not one sample: 1036
    # samples clear +25 between 14.205 and 14.260 d, the tallest reaching
    # 87.8 MPa absolute. Away from 14.20-14.38 d nothing exceeds +20, which is
    # what sets the axis at +22.
    axp.text(14.23 + 0.25, 20.6,
             f"spike burst to {p_full[i_pk]:.1f} MPa\n"
             r"($\Delta p$ = " f"{dp_full.max():+.0f}), off scale",
             ha="left", va="top", color=GRN, fontsize=10.5,
             linespacing=1.35, path_effects=HALO, zorder=7)
    axp.annotate("", xy=(14.23, 22), xytext=(14.23, 18.9),
                 arrowprops=dict(arrowstyle="-|>", color=GRN, lw=1.5))
    axp.text(2.94, 1.2, "gauge bled off", rotation=90, ha="center", va="bottom",
             fontsize=10.5, color=MUTED, path_effects=HALO, zorder=7)

    for a in (axq, axp):
        a.axvspan(T_SHUT1, T_RESUME, color=GRID, alpha=0.65, zorder=0)
        for x in (T_RESUME, T_SHUT):
            a.axvline(x, color=INK, lw=1.7, ls="--", zorder=1)

    OUT.mkdir(parents=True, exist_ok=True)
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"slide_front_fit.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {OUT}/slide_front_fit.png and .pdf")
    return D_trig, D_back


if __name__ == "__main__":
    main()
