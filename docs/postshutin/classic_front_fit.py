#!/usr/bin/env python3
"""The classic r-t plot: sqrt(t) triggering front and back front, fitted.

The standard figure in the induced-seismicity literature (Shapiro et al. 1997,
2002; Parotidis et al. 2004): every event as a dot in distance-vs-time, with two
smooth analytical curves fitted as ENVELOPES of the cloud.

  TRIGGERING FRONT       r(t) = sqrt(4 pi D t)
    Shapiro's result: the leading edge of seismicity triggered by pore-pressure
    diffusion from a point source. Fitted as an UPPER envelope.

  BACK FRONT             r(t) = sqrt( 4 D t (t - t_s)/t_s * ln(t/(t - t_s)) )
    Parotidis's result, derived for 2D radial (Theis/E1) diffusion, which is the
    right geometry for a fault plane -- a 3D point source gives a 6D prefactor
    instead of 4D. After shut-in at t_s the near-well pressure falls while the
    perturbation keeps spreading, so seismicity ceases near the well first and
    the quiet region expands. Fitted as a LOWER envelope of the post-shut-in
    events.

D IS FITTED BY QUANTILE, NOT BY EYE. Each curve belongs to a one-parameter
family in D, so the fit is: choose D such that a stated fraction of events falls
on the correct side. That is reproducible and reportable, unlike drawing the
envelope by hand. Both fractions are printed.

TIME ORIGIN. Both formulae assume injection starting at t = 0 and running
continuously. The record does not: injection begins at 0.501 d, stops at
1.582 d, and only resumes continuously at 4.300 d. Anchoring the triggering
front at 4.300 d rather than 0.501 d therefore fits far better, measured against
the per-day observed maximum:

    origin                       D          rms dev   worst under
    first injection  0.501 d     0.0803     207 m       -404 m
    injection resumes 4.300 d    0.1018     147 m       -122 m

so 4.300 d is used, and the x axis stays on the full record so cycle 1 and the
first shut-in remain visible. Both fits are printed.

Usage:  python classic_front_fit.py
"""
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import importlib.util as iu
_s = iu.spec_from_file_location(
    "fb", str(Path(__file__).with_name("front_backfront.py")))
fb = iu.module_from_spec(_s); _s.loader.exec_module(fb)

T_ON, T_SHUT = fb.T_ON, fb.T_SHUT          # 0.501 d, 17.455 d after 2012-11-13
TS = T_SHUT - T_ON                          # injection duration, days
OUT = Path("/home/users/nberrios/3dhbi/hbi_analysis/figures/postshutin")
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
RED, BLUE = "#a8071a", "#1d4ed8"
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12,
                     "axes.labelsize": 11.5, "axes.edgecolor": MUTED,
                     "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED,
                     "legend.fontsize": 10})


T_RESUME = 4.300                            # injection resumes, cycle-2 zero


def trig_front(t, D):
    """Shapiro: r = sqrt(4 pi D t). t in days SINCE ITS OWN ORIGIN."""
    return np.sqrt(4 * np.pi * D * np.maximum(t, 0.0) * 86400.0)


def back_front(t, D, ts=TS):
    """Parotidis, 2D radial. NaN before shut-in."""
    t = np.atleast_1d(np.asarray(t, float))
    out = np.full_like(t, np.nan)
    k = t > ts
    tt = t[k]
    out[k] = np.sqrt(4 * D * (tt * (tt - ts) / ts)
                     * np.log(tt / (tt - ts)) * 86400.0)
    return out


def fit_quantile(t, r, curve, frac, lo=1e-4, hi=1e3, **kw):
    """D such that `frac` of the events lie BELOW curve(t, D).

    frac = 0.95 gives an upper envelope containing 95% of events; frac = 0.05
    gives a lower envelope with 5% below it. Monotone in D, so bisection.
    """
    def g(D):
        c = curve(t, D, **kw)
        m = np.isfinite(c)
        return (r[m] < c[m]).mean() - frac
    return brentq(g, lo, hi, xtol=1e-8, rtol=1e-10)


def main():
    t_abs, r, _ = fb.catalogue()
    ti, q = fb.rate_history()
    ETA, PHI, BETA = 0.89e-3, 0.01, 2.25e-8

    # triggering front, anchored at injection RESUMPTION (see the docstring)
    tt = t_abs - T_RESUME
    kt = tt > 0
    D_trig = fit_quantile(tt[kt], r[kt], trig_front, 0.95)
    # the conventional alternative, for the record
    ta = t_abs - T_ON
    ka = ta > 0
    D_alt = fit_quantile(ta[ka], r[ka], trig_front, 0.95)

    # back front, anchored at shut-in
    tb_all = t_abs - T_ON
    post = tb_all > TS
    D_back = fit_quantile(tb_all[post], r[post], back_front, 0.05)

    print(f"{len(t_abs)} events; {kt.sum()} after injection resumes "
          f"({T_RESUME} d), {post.sum()} after shut-in ({T_SHUT} d)\n")
    print(f"TRIGGERING FRONT   r = sqrt(4 pi D t),  t from {T_RESUME} d")
    print(f"  D = {D_trig:.4f} m^2/s  (95% of events below)"
          f"  ->  k = {D_trig*ETA*PHI*BETA:.2e} m^2")
    print(f"  [alternative, t from {T_ON} d: D = {D_alt:.4f}, a worse envelope]")
    print(f"BACK FRONT         r = sqrt(4 D t(t-t_s)/t_s ln(t/(t-t_s)))")
    print(f"  D = {D_back:.4f} m^2/s  (5% of post-shut-in events below)"
          f"  ->  k = {D_back*ETA*PHI*BETA:.2e} m^2")
    print(f"\nratio D_back/D_trig = {D_back/D_trig:.1f}x")
    print(f"model bounds: kpmin -> {1e-15/(ETA*PHI*BETA):.4f}, "
          f"kpmax -> {2.5e-13/(ETA*PHI*BETA):.3f} m^2/s")

    fig, (ax, axq) = plt.subplots(
        2, 1, figsize=(11.5, 7.6), dpi=200, sharex=True,
        gridspec_kw=dict(height_ratios=[3.4, 1.0], hspace=0.06))

    ax.scatter(t_abs, r, s=3.5, alpha=0.22, color=MUTED, lw=0,
               label=f"{len(t_abs)} events")
    g = np.linspace(T_RESUME + 1e-3, t_abs.max(), 800)
    ax.plot(g, trig_front(g - T_RESUME, D_trig), "-", lw=2.8, color=RED,
            label=r"triggering front  $r=\sqrt{4\pi D t}$,  "
                  f"$D$ = {D_trig:.3f} m$^2$/s")
    gb = np.linspace(T_SHUT + 1e-3, t_abs.max(), 800)
    ax.plot(gb, back_front(gb - T_ON, D_back), "-", lw=2.8, color=BLUE,
            label=r"back front  $r=\sqrt{4Dt\frac{(t-t_s)}{t_s}"
                  r"\ln\frac{t}{t-t_s}}$,  " f"$D$ = {D_back:.3f} m$^2$/s")
    for x, lab, yl in ((T_RESUME, "injection resumes  ($t=0$ for $r=\\sqrt{4\\pi Dt}$)",
                        640),
                       (T_SHUT, "shut-in  ($t_s$)", 1700)):
        ax.axvline(x, color=INK, lw=1.6, ls="--")
        ax.text(x - 0.18, yl, lab, rotation=90, ha="right", va="top",
                fontsize=9.5, color=INK)
    ax.axvspan(1.582, T_RESUME, color=GRID, alpha=0.65, zorder=0)
    ax.text(2.94, 60, "1st shut-in", rotation=90, ha="center", va="bottom",
            fontsize=9.5, color=MUTED)
    ax.set(ylabel="Distance from injection point (m)",
           xlim=(0, t_abs.max()), ylim=(0, 1750))
    ax.set_title("Cooper Basin Habanero 4, November 2012 — "
                 "triggering front and back front")
    ax.legend(loc="upper left", framealpha=0.92)
    ax.grid(alpha=0.28, color=GRID)

    axq.fill_between(ti, 0, q, color=MUTED, alpha=0.45, lw=0)
    for x in (T_RESUME, T_SHUT):
        axq.axvline(x, color=INK, lw=1.6, ls="--")
    axq.axvspan(1.582, T_RESUME, color=GRID, alpha=0.65, zorder=0)
    axq.set(xlabel="Days since 2012-11-13", ylabel="q (L/s)",
            xlim=(0, t_abs.max()), ylim=(0, 70))
    axq.grid(alpha=0.28, color=GRID)

    OUT.mkdir(parents=True, exist_ok=True)
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"classic_front_fit.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {OUT}/classic_front_fit.png")
    return D_trig, D_back


if __name__ == "__main__":
    main()
