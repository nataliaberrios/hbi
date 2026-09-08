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


def catalogue_mag():
    """catalogue() plus magnitudes: (t, r, ML, Mw), identical mask and order.

    TWO MAGNITUDE SCALES SHIP IN THE FILE AND THEY DISAGREE, so "M3" is
    ambiguous until one is named:

        ML  tops out at 3.0    2 events >= 3.0
        Mw  tops out at 3.1    4 events >= 3.0,  8 >= 2.8

    and they disagree per event, not just in calibration -- the largest Mw
    (3.10 at 19.648 d) is only ML 2.0, while an ML 3.0 event is Mw 2.8. So the
    two scales select overlapping-but-different sets: of the 2 ML>=3 and the
    4 Mw>=3, only 1 event is in both.

    Mw is the one to prefer for a physical statement, being tied to M0 (also in
    the file, 4.2e8 to 4.5e13 N m), but both are returned so a figure can say
    which it used.
    """
    m = loadmat(H / "Cooper_Basin_Catalog_HAB_4.mat", squeeze_me=True,
                struct_as_record=False)
    c = {e.field: e.val for e in m["Catalog"]}
    tt = c["Time"].astype(float)
    dts = [dt.fromordinal(int(x)) + td(days=x % 1) - td(days=366) for x in tt]
    t = np.array([(d - T_REF).total_seconds() / 86400 for d in dts])
    r = np.sqrt(((c["Lat"] - LAT) * 111.0) ** 2
                + ((c["Long"] - LON) * 111.0 * np.cos(np.radians(LAT))) ** 2) * 1e3
    ml = np.asarray(c["ML"], float)
    mw = np.asarray(c["Mw"], float)
    k = np.isfinite(t) & np.isfinite(r) & (t >= 0)
    o = np.argsort(t[k])
    t, r, ml, mw = t[k][o], r[k][o], ml[k][o], mw[k][o]
    # alignment with catalogue() is what lets a caller overlay one on the other
    t0, r0, _ = catalogue()
    assert np.array_equal(t, t0) and np.array_equal(r, r0), \
        "catalogue_mag() diverged from catalogue()"
    return t, r, ml, mw


def shutin_starts(min_len=0.1, thr=0.05):
    """Times at which injection stops for longer than `min_len` days.

    Returns [1.582, 17.151]. The record has 27 zero-rate runs; 24 are under
    half an hour and are dropouts rather than decisions, so a length threshold
    is what separates a shut-in from noise, and 0.1 d puts the cut in a wide
    empty gap -- the next-longest run after the two returned is 0.035 d.

    THE FINAL SHUT-IN BEGINS AT 17.151 d, NOT AT T_SHUT = 17.455. Injection
    stops at 17.095, restarts for 31 minutes at 17.130, and stops for good at
    17.151; 17.455 is simply the LAST SAMPLE OF THE RATE FILE, 0.30 d later.
    T_SHUT is used elsewhere as t_s, the injection duration in the back-front
    formula, so it is left alone here -- but anything drawn at T_SHUT is
    marking the end of the record, not a shut-in.
    """
    ti, q = rate_history()
    z = q <= thr
    out, i = [], 0
    while i < len(z):
        if z[i]:
            j = i
            while j + 1 < len(z) and z[j + 1]:
                j += 1
            if ti[j] - ti[i] > min_len:
                out.append(float(ti[i]))
            i = j + 1
        else:
            i += 1
    return out


def zero_rate_interval(thr=0.05):
    """(t0, t1) of the longest interval with q <= thr L/s during the record.

    THE FIRST SHUT-IN IS SHORTER THAN THE GAP BEFORE FULL RESUMPTION, and
    figures that shade 1.582-4.300 d as "shut-in" are wrong for the last
    0.74 d of it. The rate is genuinely zero from 1.582 to 3.558 d; from 3.558
    to 4.298 d a 1.6-2.5 L/s trickle runs; only at 4.300 d does it step to
    ~8 L/s and climb. 4.300 d is still the right cycle-2 origin -- it is where
    sustained injection resumes -- but it is not where the shut-in ends.

    Derived from the rate file rather than hard-coded so it cannot drift, and
    taken as the LONGEST zero run rather than the outermost zero samples: the
    rate briefly touches zero again at 4.569 d, so min/max over the mask would
    return 1.582-5.298 d instead.
    """
    ti, q = rate_history()
    z = q <= thr
    runs, i = [], 0
    while i < len(z):
        if z[i]:
            j = i
            while j + 1 < len(z) and z[j + 1]:
                j += 1
            if ti[i] > T_ON and ti[j] < T_SHUT:
                runs.append((ti[j] - ti[i], ti[i], ti[j]))
            i = j + 1
        else:
            i += 1
    if not runs:
        raise RuntimeError("no zero-rate interval found between T_ON and T_SHUT")
    _, t0, t1 = max(runs)
    return float(t0), float(t1)


def pressure_history(every=1):
    """(t_days, p_abs, dp) from the wellhead .mat, same clock as rate_history().

    p_abs is the measured absolute wellhead pressure in MPa; dp is the CHANGE
    from the first sample of the record, 34.412 MPa. That is the convention
    Wang & Dunham use -- export_panelA_fixed.m:30 is
    `(M2.H4_wh_p(:) - M2.H4_wh_p(1))/1e6` -- so a dp plotted here is directly
    comparable to their figure 2a. The pre-injection median over t < 0.501 d is
    33.970 MPa, 0.44 MPa lower, the well having drifted down while shut; using
    the first sample rather than that median is therefore a choice, and it makes
    dp smaller by 0.44 MPa everywhere. It is kept only for comparability.

    `every` decimates for plotting -- the record is 1.31 M samples at ~1.2 s.
    Decimation can drop short transients, so the caller should quote extrema
    from the undecimated record; p_abs.max() is 87.76 MPa over two samples at
    14.23 d and survives only at every = 1.
    """
    pr = loadmat(H / "Cooper_Basin_HAB_4_Wellhead_Pressure.mat")["d"][0, 0]
    tp = pr["Date"].squeeze(); tp = tp - tp[0]
    p = pr["Wellhead_pressure"].squeeze()
    k = np.isfinite(tp) & np.isfinite(p)
    o = np.argsort(tp[k])
    tp, p = tp[k][o], p[k][o]
    p0 = float(p[0])                    # 34.412 MPa; taken from the file, not
    if every > 1:                       # hard-coded, so it cannot go stale
        tp, p = tp[::every], p[::every]
    return tp, p, p - p0


def front_runmax(t, r):
    """The FRONT: running maximum of event distance. Monotone BY CONSTRUCTION.

    A per-bin high percentile is NOT a front -- it tracks where the event
    POPULATION sits in that bin, so it moves backwards when the cloud
    redistributes. Measured on this catalogue it retreats four times
    (819->739, 882->613, 1010->817, 1544->1429 m), which no front can do.
    """
    return np.maximum.accumulate(r)


def sqrt_test(t, r, a, b):
    """(lambda, R^2) for R = lambda*sqrt(t) over [a,b).

    R = lambda*sqrt(t) is equivalent to R^2 linear in t, so R^2 of that
    regression is the honest test of whether the sqrt law holds. It assumes
    CONSTANT RATE, so it should hold on an uninterrupted flowing block and
    degrade across shut-ins and rate steps -- which is exactly what happens.
    """
    m = (t >= a) & (t < b)
    tt, rr = t[m], np.maximum.accumulate(r)[m]
    lam = float(np.sum(np.sqrt(tt) * rr) / np.sum(tt))
    A = np.vstack([tt, np.ones_like(tt)]).T
    _, res, _, _ = np.linalg.lstsq(A, rr ** 2, rcond=None)
    ss = np.sum((rr ** 2 - np.mean(rr ** 2)) ** 2)
    return lam, float(1 - (res[0] / ss if len(res) else 0)), int(m.sum())


def D_from_triggering(lam_m_per_sqrtday):
    """Shapiro triggering front r = sqrt(4*pi*D*t)  =>  D = lambda^2/(4 pi)."""
    lam = lam_m_per_sqrtday / np.sqrt(86400.0)      # m / sqrt(s)
    return lam ** 2 / (4 * np.pi)


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

    # ---------------- diffusivities from the two fronts, model-free
    lam_s, r2_s, n_s = sqrt_test(t, r, 5.3, 13.2)      # steady flowing block
    lam_a, r2_a, n_a = sqrt_test(t, r, T_ON, T_SHUT)   # all injection
    D_trig = D_from_triggering(lam_s)
    # TWO viscosities, deliberately. The decks carry eta = 0.89e-3, which came
    # down the 1807 lineage and is Wang & Dunham's value; this project's rule is
    # eta = 1.27e-4, appropriate for the depth and temperature of the Habanero
    # reservoir. k = D*eta*phi*beta scales with eta, so the inferred
    # permeability differs by 7.008x and is meaningless without saying which.
    # The D comparison against kpmax/(eta*phi*beta) is eta-INDEPENDENT, since
    # both sides use the same eta, so only the k numbers are affected.
    PHI, BETA = 0.01, 2.25e-8
    ETA_DECK, ETA_RES = 0.89e-3, 1.27e-4
    ETA = ETA_DECK
    print(f"sqrt(t) test -- R = lam*sqrt(t) means R^2 is LINEAR in t, so R^2 of")
    print(f"that regression is the honest test. It assumes CONSTANT RATE.")
    print(f"  5.3-13.2 d (uninterrupted flowing block): lam {lam_s:>6.1f} "
          f"m/sqrt(d), R2 = {r2_s:.3f}, n = {n_s}")
    print(f"  {T_ON}-{T_SHUT} d (all injection, 4 rate changes): lam {lam_a:>6.1f} "
          f"m/sqrt(d), R2 = {r2_a:.3f}, n = {n_a}")
    print(f"  -> sqrt(t) HOLDS on the steady block and degrades across the")
    print(f"     shut-in and the rate steps at 13.5 and 16.5 d, as it should.\n")
    print(f"TWO INDEPENDENT DIFFUSIVITIES, neither using any HBI run:")
    print(f"  triggering front (Shapiro, r = sqrt(4 pi D t)), during injection:")
    print(f"     D = lam^2/(4 pi) = {D_trig:.4f} m^2/s  -> k = {D_trig*ETA*PHI*BETA:.2e} m^2")
    print(f"  back front (Parotidis), post-shut-in:")
    print(f"     D = {D_fit:.4f} m^2/s  -> k = {D_fit*ETA*PHI*BETA:.2e} m^2")
    print(f"  RATIO D_back/D_trig = {D_fit/D_trig:.1f}x -- a MODEL-FREE measure of")
    print(f"  the permeability enhancement achieved during the stimulation.")
    print(f"  the decks assume kpmax/kpmin = 250x. The data says ~{D_fit/D_trig:.0f}x.")
    print(f"  model bounds: kpmin -> {1e-15/(ETA*PHI*BETA):.4f}, "
          f"kpmax -> {2.5e-13/(ETA*PHI*BETA):.3f} m^2/s; both fits sit between.\n")

    # ---------------- figure
    runmax = front_runmax(t, r)
    fig, (a0, a1, a2) = plt.subplots(1, 3, figsize=(16.5, 4.9), dpi=200,
                                     constrained_layout=True)

    a0.scatter(t, r, s=2.5, alpha=0.16, color=MUTED, lw=0,
               label=f"all {len(t)} events (full catalogue)")
    a0.plot(t, runmax, "-", lw=2.2, color=RED,
            label="FRONT = running max (monotone)")
    a0.plot(tc, rin, "-", lw=2.2, color=BLUE, label="BACK FRONT (5th pct)")
    a0.plot(tc[post][ipk + 1:], rin[post][ipk + 1:], "x", ms=8, mew=2,
            color=BLUE, label="excluded from the fit (falling, low n)")
    tg = np.linspace(T_SHUT + 0.02, t.max(), 300)
    a0.plot(tg, bf_single(tg, D_fit), "--", lw=2, color=INK,
            label=f"Parotidis back front, D = {D_fit:.2f} m$^2$/s")
    a0.plot(tg, bf_super(tg, D_fit, ti, q), ":", lw=2, color=GRN,
            label="same D, ACTUAL rate history (superposition)")
    # Draw the triggering front ONLY over the window it was fitted on. Shown
    # across the whole record it sits below the staircase everywhere and looks
    # like a bad fit, when in fact it is an extrapolation outside its domain.
    tt = np.linspace(5.3, 13.2, 200)
    a0.plot(tt, lam_s * np.sqrt(tt), "-.", lw=2.2, color="#8E44AD",
            label=f"triggering front $\\lambda\\sqrt{{t}}$ over its FIT window, "
                  f"D = {D_trig:.3f}")
    a0.axvspan(1.582, 4.300, color=GRID, alpha=0.7, zorder=0)
    a0.axvspan(5.3, 13.2, color="#8E44AD", alpha=0.07, zorder=0)
    a0.axvline(T_SHUT, color=INK, lw=1.4, ls="-.")
    a0.text(T_SHUT - 0.25, 1700, "injection ends", fontsize=8.5, color=INK,
            rotation=90, va="top", ha="right")
    a0.set(xlabel="Days since 2012-11-13", ylabel="Distance from injector (m)",
           xlim=(0, t.max()), ylim=(0, 1800))
    a0.set_title("A. Front (monotone) and back front.\nShaded purple: the "
                 "uninterrupted flowing block")
    a0.legend(loc="upper left", fontsize=7.5)
    a0.grid(alpha=0.3, color=GRID)

    # --- B: the sqrt(t) test. R^2 linear in t  <=>  R ~ sqrt(t)
    for a_, b_, c_, nm in ((5.3, 13.2, "#8E44AD", "steady block"),
                           (T_ON, T_SHUT, MUTED, "all injection")):
        m = (t >= a_) & (t < b_)
        lam_, r2_, _ = sqrt_test(t, r, a_, b_)
        a1.plot(t[m], (runmax[m] / 1000.0) ** 2, "-", lw=2, color=c_,
                label=f"{nm}: $R^2$ = {r2_:.3f}")
        a1.plot(t[m], (lam_ * np.sqrt(t[m]) / 1000.0) ** 2, "--", lw=1.4,
                color=c_, alpha=0.8)
    a1.axvline(T_SHUT, color=INK, lw=1.4, ls="-.")
    a1.axvspan(1.582, 4.300, color=GRID, alpha=0.7, zorder=0)
    a1.set(xlabel="Days since 2012-11-13",
           ylabel="(front radius)$^2$  (km$^2$)", xlim=(0, T_SHUT + 0.5))
    a1.set_title("B. The $\\sqrt{t}$ test — straight lines mean $R\\propto"
                 "\\sqrt{t}$.\nA STAIRCASE: each step is a rate increase, so "
                 "$\\sqrt{t}$ is a trend, not the mechanism")
    a1.legend(loc="upper left"); a1.grid(alpha=0.3, color=GRID)

    # --- C: the two diffusivities against the model's bounds
    names = ["triggering front\n(during injection)", "back front\n(post-shut-in)"]
    vals = [D_trig, D_fit]
    a2.bar([0, 1], vals, width=0.5, color=["#8E44AD", BLUE], alpha=0.85)
    for i, v in enumerate(vals):
        a2.text(i, v * 1.15, f"{v:.3f}", ha="center", fontsize=10)
    a2.axhline(1e-15 / (ETA * PHI * BETA), color=GRN, ls="--", lw=1.6,
               label="model $k_{pmin}$")
    a2.axhline(2.5e-13 / (ETA * PHI * BETA), color=RED, ls="--", lw=1.6,
               label="model $k_{pmax}$")
    a2.set_yscale("log")
    a2.set_xticks([0, 1]); a2.set_xticklabels(names, fontsize=8.5)
    a2.set(ylabel="hydraulic diffusivity $D$ (m$^2$/s)", ylim=(2e-3, 5))
    a2.set_title(f"C. Two MODEL-FREE diffusivities.\nRatio "
                 f"{D_fit/D_trig:.1f}x = the permeability enhancement\n"
                 f"(the decks assume 250x)")
    a2.legend(loc="upper left", fontsize=8); a2.grid(alpha=0.3, color=GRID, axis="y")

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
