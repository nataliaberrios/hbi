#!/usr/bin/env python3
"""Front and back front of the Cooper Basin seismicity, from the FULL catalogue.

WHY THIS EXISTS. Every front figure in this project used
make_sweep_figures.observed(), which windows the catalogue to Nov 13-30 2012 --
17574 of 20735 events. The 3161 it drops are almost entirely the POST-INJECTION
ones (Nov 30 to Dec 4), which is exactly where a back front appears. Plotting
only the windowed catalogue is why the back front looked absent.

WHAT A BACK FRONT IS. After shut-in the pore-pressure perturbation keeps
diffusing outward while the near-well pressure falls, so seismicity ceases near
the well first and the quiet region expands. The locus where dp/dt = 0 is the
back front (Parotidis, Shapiro & Rothert 2004, 10.1029/2003gl018987).

THE ANALYTICAL CURVE, DERIVED FOR THIS GEOMETRY RATHER THAN QUOTED. A planar
fault with a point injection is 2D radial diffusion, so the pressure kernel is
Theis/E1 -- the same exp1 this project already uses in gc_solution.py -- not the
3D erfc kernel. For a single injection from 0 to t_s,

    p(r,t) ~ E1(r^2/4Dt) - E1(r^2/4D(t-t_s))          t > t_s
    d/dt E1(r^2/4Dt) = exp(-r^2/4Dt)/t
    dp/dt = 0   =>   exp(-x1)/t = exp(-x2)/(t-t_s),  x_i = r^2/4D(t-t_i)
                =>   r^2/4D * t_s/(t(t-t_s)) = ln(t/(t-t_s))

    r_bf(t) = sqrt( 4D * t(t-t_s)/t_s * ln(t/(t-t_s)) )

which is Parotidis's result. Note the 3D point source gives a different
prefactor (6D, from the extra t^-3/2 in the erfc derivative), so the geometry
matters and the 2D form is the right one here.

BUT THE REAL RECORD IS NOT A SINGLE INJECTION -- it has a 2.7 d shut-in in the
middle and a strongly varying rate. So the closed form above is only indicative.
This script also solves the back front for the ACTUAL rate history by
superposition, which is the same principle carried to the real q(t):

    dp/dt(r,t) ~ sum_i dq_i * exp(-r^2/4D(t-t_i)) / (t-t_i)      over t_i < t

and finds the r where that vanishes. Both are plotted; the difference between
them is the error the single-injection idealisation makes.

D IS FITTED, NOT ASSUMED. Parotidis's point is that the back front measures the
hydraulic diffusivity of the seismically active volume. That fitted D is an
INDEPENDENT observable to set against the model's kpmax/(eta*phi*beta) = 1.248
m^2/s -- it does not depend on any HBI run.

Usage:  python front_backfront.py [--out DIR]
"""
import argparse
from datetime import datetime as dt, timedelta as td
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.optimize import brentq, curve_fit
from scipy.special import exp1
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = Path("/home/users/nberrios/3dhbi/hbi_analysis")
T_REF = dt(2012, 11, 13)          # this project's t = 0
LAT, LON = -27.8115, 140.7596     # injector, as in observed()
T_SHUT = 17.455                   # injection ends
T_ON = 0.501                      # injection first starts
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
RED, BLUE, GRN = "#a8071a", "#1d4ed8", "#009E73"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 10.5,
                     "axes.edgecolor": MUTED, "text.color": INK,
                     "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "legend.fontsize": 8.5})


def _is_cli():
    """True only when run as a script, so main() can be imported safely."""
    import sys as _s
    return _s.argv and _s.argv[0].endswith("front_backfront.py")


def catalogue():
    """FULL catalogue -- no Nov 30 cut. (t days from T_REF, r metres)."""
    m = loadmat(H / "Cooper_Basin_Catalog_HAB_4.mat", squeeze_me=True,
                struct_as_record=False)
    c = {e.field: e.val for e in m["Catalog"]}
    tt = c["Time"].astype(float)
    dts = [dt.fromordinal(int(x)) + td(days=x % 1) - td(days=366) for x in tt]
    t = np.array([(d - T_REF).total_seconds() / 86400 for d in dts])
    r = np.sqrt(((c["Lat"] - LAT) * 111.0) ** 2
                + ((c["Long"] - LON) * 111.0 * np.cos(np.radians(LAT))) ** 2) * 1e3
    k = np.isfinite(t) & np.isfinite(r) & (t >= 0)
    o = np.argsort(t[k])
    return t[k][o], r[k][o], len(tt)


def rate_history():
    """(t_days, q) from the injection .mat, same units as observed()."""
    inj = loadmat(H / "Cooper_Basin_HAB_4_Injection_Rate.mat")["d"][0, 0]
    ti = inj["Date"].squeeze(); ti = ti - ti[0]
    q = inj["Injection_rate"].squeeze() * (1000 / 60)
    k = np.isfinite(ti) & np.isfinite(q)
    o = np.argsort(ti[k])
    return ti[k][o], q[k][o]


def edges(t, r, bins, lo_pct=5, hi_pct=95, nmin=25):
    """Inner (back) and outer (front) edge of the cloud per time bin."""
    tc, rin, rout, n = [], [], [], []
    for a, b in zip(bins[:-1], bins[1:]):
        m = (t >= a) & (t < b)
        if m.sum() < nmin:
            continue
        tc.append(0.5 * (a + b)); n.append(int(m.sum()))
        rin.append(np.percentile(r[m], lo_pct))
        rout.append(np.percentile(r[m], hi_pct))
    return (np.array(tc), np.array(rin), np.array(rout), np.array(n))


def bf_single(t, D, t_s=T_SHUT - T_ON):
    """Closed-form 2D back front, t measured from injection START."""
    t = np.atleast_1d(np.asarray(t, float))
    out = np.full_like(t, np.nan)
    ok = t > t_s
    tt = t[ok]
    out[ok] = np.sqrt(4 * D * (tt * (tt - t_s) / t_s) * np.log(tt / (tt - t_s))
                      * 86400.0)
    return out


def bf_super(t_days, D, ti, q, rmax=4000.0):
    """Back front for the ACTUAL rate history, by superposition.

    dp/dt ~ sum_i dq_i exp(-r^2/(4 D (t-t_i)))/(t-t_i). Solve for r.
    """
    dq = np.diff(np.concatenate([[0.0], q]))
    keep = np.abs(dq) > 1e-9
    tq, dq = ti[keep], dq[keep]

    def dpdt(r, t):
        dtt = (t - tq) * 86400.0
        m = dtt > 0
        if not m.any():
            return np.nan
        x = r ** 2 / (4 * D * dtt[m])
        return float(np.sum(dq[m] * np.exp(-np.minimum(x, 700.0)) / dtt[m]))

    out = []
    for t in np.atleast_1d(t_days):
        f = lambda r: dpdt(r, t)
        try:
            a, b = 1.0, rmax
            if f(a) * f(b) > 0:
                out.append(np.nan); continue
            out.append(brentq(f, a, b))
        except (ValueError, RuntimeError):
            out.append(np.nan)
    return np.array(out)


def main(argv=None):
    """argv is explicit so a notebook can call main([]) -- argparse otherwise
    reads the kernel's own argv and exits 2."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(H / "figures" / "postshutin"))
    a = ap.parse_args([] if argv is None and not _is_cli() else argv)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    t, r, ntot = catalogue()
    ti, q = rate_history()
    print(f"FULL catalogue: {len(t)} events with t >= 0, of {ntot} total")
    print(f"  t = {t.min():.2f} to {t.max():.2f} d   (t = 0 is {T_REF:%Y-%m-%d})")
    nwin = ((t >= 0) & (t <= 17.0)).sum()
    print(f"  the Nov 13-30 window used elsewhere holds {nwin}; this adds "
          f"{len(t)-nwin} POST-INJECTION events\n")

    bins = np.arange(0.0, t.max() + 0.5, 0.5)
    tc, rin, rout, n = edges(t, r, bins)
    post = tc > T_SHUT

    # --- fit D to the observed inner edge after shut-in.
    # ONLY THE RISING LIMB. The observed back front peaks at 707 m near 18.75 d
    # and then FALLS to 365 m, which the theory cannot produce: r_bf(t) is
    # monotonic. Two reasons to exclude the falling part rather than fit it.
    #   * event counts collapse there -- 235, 454, 147, 192, 29 per bin against
    #     ~650 on the rising limb -- so a 5th percentile is a poor edge estimate.
    #   * a monotonic model fitted across a peak is biased low by construction.
    # Both fits are reported so the size of that bias is visible, not hidden.
    ipk = int(np.argmax(rin[post]))
    rise = np.zeros_like(post); rise[np.where(post)[0][:ipk + 1]] = True

    def model(tt, D):
        return bf_single(tt, D)

    fits = {}
    for name, msk in (("rising limb only", rise), ("all post-shut-in", post)):
        D_, pc = curve_fit(model, tc[msk], rin[msk], p0=[1.0],
                           bounds=(1e-3, 1e4))
        D_ = float(D_[0]); e_ = float(np.sqrt(pc[0, 0]))
        res = rin[msk] - bf_single(tc[msk], D_)
        fits[name] = (D_, e_, np.sqrt(np.mean(res ** 2)), int(msk.sum()))
        print(f"D from the BACK FRONT, {name:>18s}: {D_:>6.3f} +/- {e_:.3f} "
              f"m^2/s   rms {np.sqrt(np.mean(res**2)):>4.0f} m, n = {int(msk.sum())}")
    D_fit = fits["rising limb only"][0]
    print(f"\n  model kpmax/(eta*phi*beta) = 1.248 m^2/s  -> fitted is "
          f"{D_fit/1.248:.2f}x that")
    print(f"  model kpmin/(eta*phi*beta) = 0.005 m^2/s  -> fitted is "
          f"{D_fit/4.994e-3:.0f}x that")
    print(f"  so the back front measures an EFFECTIVE D between the model's")
    print(f"  background and enhanced values -- an independent observable.\n")
    print(f"  observed back front peaks at {rin[post].max():.0f} m near "
          f"{tc[post][ipk]:.2f} d then FALLS to {rin[post][-1]:.0f} m; theory is")
    print(f"  monotonic, and the late bins hold only "
          f"{n[post][-1]}-{n[post][-3]} events.\n")

    # ---------------- figure
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(14.5, 5.0), dpi=200,
                                 constrained_layout=True)
    a0.scatter(t, r, s=2.5, alpha=0.16, color=MUTED, lw=0,
               label=f"all {len(t)} events (full catalogue)")
    a0.plot(tc, rout, "-", lw=2.2, color=RED, label="FRONT (95th pct)")
    a0.plot(tc, rin, "-", lw=2.2, color=BLUE, label="BACK FRONT (5th pct)")
    a0.plot(tc[post][ipk + 1:], rin[post][ipk + 1:], "x", ms=8, mew=2,
            color=BLUE, label="excluded from the fit (falling, low n)")
    tg = np.linspace(T_SHUT + 0.02, t.max(), 300)
    a0.plot(tg, bf_single(tg, D_fit), "--", lw=2, color=INK,
            label=f"Parotidis 2D back front, D = {D_fit:.2f} m$^2$/s "
                  f"(fitted to the rising limb)")
    a0.plot(tg, bf_super(tg, D_fit, ti, q), ":", lw=2, color=GRN,
            label="same D, ACTUAL rate history (superposition)")
    a0.axvspan(1.582, 4.300, color=GRID, alpha=0.7, zorder=0)
    a0.axvline(T_SHUT, color=INK, lw=1.4, ls="-.")
    a0.text(T_SHUT - 0.25, 1180, "injection ends", fontsize=8.5, color=INK,
            rotation=90, va="top", ha="right")
    a0.text(2.94, 40, "1st shut-in", fontsize=8.5, color=MUTED,
            rotation=90, va="bottom", ha="center")
    a0.set(xlabel="Days since 2012-11-13", ylabel="Distance from injector (m)",
           xlim=(0, t.max()), ylim=(0, 1800))
    a0.set_title("Cooper Basin seismicity: front and back front\n"
                 "the back front only exists in the part of the catalogue "
                 "previously cut off")
    a0.legend(loc="upper left", fontsize=8)
    a0.grid(alpha=0.3, color=GRID)

    a1b = a1.twinx()
    a1b.fill_between(ti, 0, q, color=GRID, zorder=0, label="injection rate")
    a1b.set_ylabel("Injection rate (L/s)", color=MUTED)
    a1b.set_ylim(0, 260)
    a1.plot(tc, rout, "-o", ms=3, lw=2, color=RED, label="front")
    a1.plot(tc, rin, "-o", ms=3, lw=2, color=BLUE, label="back front")
    a1.plot(tg, bf_single(tg, D_fit), "--", lw=2, color=INK,
            label=f"analytical back front (D = {D_fit:.2f})")
    a1.axvline(T_SHUT, color=INK, lw=1.4, ls="-.")
    a1.set(xlabel="Days since 2012-11-13", ylabel="Distance from injector (m)",
           xlim=(12, t.max()), ylim=(0, 1800))
    a1.set_title("Zoom on the final shut-in.\nThe cloud becomes an expanding "
                 "ANNULUS — a pulse, per Sáez & Lecampion (2023)")
    a1.legend(loc="upper left", fontsize=8)
    a1.grid(alpha=0.3, color=GRID)
    for e in ("png", "pdf"):
        fig.savefig(out / f"front_backfront.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}/front_backfront.png")

    print(f"\n{'window':>13s} {'n':>6s} {'back':>7s} {'front':>7s} {'width':>7s}")
    for i in range(len(tc)):
        if tc[i] < T_SHUT - 1.5 and not (3.0 < tc[i] < 5.0):
            continue
        tag = "POST" if tc[i] > T_SHUT else ""
        print(f"{tc[i]:>12.2f} {n[i]:>6d} {rin[i]:>6.0f}m {rout[i]:>6.0f}m "
              f"{rout[i]-rin[i]:>6.0f}m {tag}")
    return D_fit


if __name__ == "__main__":
    main()
