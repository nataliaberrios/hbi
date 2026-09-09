#!/usr/bin/env python3
"""One figure per cycle-2 arm: front-vs-time, front-vs-volume, seismicity, dp.

WHAT THIS IS FOR. Across all four arms, the front and the wellhead pressure
respond to tau_0 with OPPOSITE SIGN: raising tau_0 grows the front and shrinks
dp. Since the front is already at or above the observed value where dp is worst,
there is no direction to move that improves both -- every sweep slides along a
trade-off instead of closing on a fit. That claim is easier to see than to read,
which is what these figures are: within one arm the only thing that changes is
tau_0, so each panel shows one parameter's whole effect.

FOUR PANELS PER ARM

  (a) front against time, and (b) front against cumulative injected volume.
      BOTH USE THE PROJECT'S OWN DEFINITIONS, taken from
      hbi_analysis/notebooks/cooper_basin_plots-30.cleaned.ipynb rather than
      reinvented here:

        observed front   cell 49/51 -- events sorted by time, binned 100 at a
                         time, keeping those between the 90th and 95th
                         percentile of each bin's distance, with the initial
                         cloud radius (median of the first ten) subtracted.
        model front      cell 50 -- one point per cell, at the first time that
                         cell's slip exceeds the threshold, plotted against
                         |x|. This is what run_data already returns.
        fit              cell 55 -- lambda = sum(sqrt(x) R) / sum(x), least
                         squares THROUGH THE ORIGIN, which is right because the
                         cloud radius has already been subtracted.

      Both fronts are SCATTERS of points, so both are drawn as points with one
      lambda*sqrt line through them. The comparison is the ratio of lambdas,
      printed and in the legend.

      An earlier version of this figure used a running maximum and a
      free-offset fit instead. That changed the quantity every number was
      measured against, and the free offset re-solved a problem the origin
      subtraction had already solved.

  (c) the seismicity itself -- every event after t0 as a dot, with the same
      model fronts over it. This is the panel that says whether a front is
      "too big": a model front well above the cloud is claiming slip where
      nothing was recorded.

  (d) wellhead PRESSURE CHANGE. dp, not absolute pressure, and this is not a
      presentational choice -- the model's absolute wellhead carries a +-2 MPa
      uncertainty from the assumed density of the 4077 m water column alone
      (33.8 MPa at 1000 kg/m3, 31.8 at 1050, 35.8 at 950), against a mean
      observed dp of 7.9 MPa. Absolute agreement inside 2 MPa is unfalsifiable.
      Model dp is ppw - P_STATIC; observed dp is pm - P_REF. See
      score_cycle2.py's docstring for the three-baseline trap.

READ (a) AND (d) TOGETHER, one arm at a time. The curves fan the same way in
(a) and the opposite way in (d), and no single tau_0 sits in the middle of both.

Usage:
  python compare_arms.py              # all four arms
  python compare_arms.py --arms 2 4
"""
import argparse
import importlib.util as iu
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = "/home/users/nberrios/3dhbi/hbi_analysis"
_s = iu.spec_from_file_location("sf", H + "/make_sweep_figures.py")
sf = iu.module_from_spec(_s); _s.loader.exec_module(sf)
sys.path.insert(0, H + "/notebooks")
from sim_curves import _deck, _ff
_f = iu.spec_from_file_location(
    "fb", "/home/users/nberrios/3dhbi/hbi_git/docs/postshutin/front_backfront.py")
fb = iu.module_from_spec(_f); _f.loader.exec_module(fb)

T0 = 4.300
TMAX = 13.2
P_STATIC = sf.P0 - sf.RHO * sf.G * sf.HW / 1e6      # 33.805 MPa
P_REF = 34.412
OUT = Path(H) / "figures" / "cycle2"
INK, MUTED, GRID, OBSC = "#1a1a19", "#6b6b66", "#d8d8d4", "#a8071a"

ARMS = {
    1: (range(632950, 632955), r"arm 1:  $\eta$ = 0.89e-3 Pa s,  "
                               r"$\beta$ = 2.25e-8 Pa$^{-1}$"),
    2: (range(632960, 632965), r"arm 2:  $\eta$ = 1.27e-4 Pa s,  "
                               r"$\beta$ = 2.25e-8 Pa$^{-1}$"),
    3: (range(632970, 632975), r"arm 3:  $\eta$ = 1.27e-4 Pa s,  "
                               r"$\beta$ = 1.58e-7 Pa$^{-1}$"),
    4: (range(632980, 632985), r"arm 4:  $\eta$ = 1.27e-4 Pa s,  "
                               r"$\beta$ = 7.72e-9 Pa$^{-1}$"),
}

plt.rcParams.update({"font.size": 10.5, "axes.titlesize": 11.5,
                     "axes.labelsize": 11, "axes.edgecolor": MUTED,
                     "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED,
                     "xtick.labelcolor": INK, "ytick.labelcolor": INK,
                     "legend.fontsize": 8.5})


def seismicity_front(t, d, bin_size=100, lower=90, upper=95):
    """THE PROJECT'S OWN FRONT DEFINITION, not a new one.

    Verbatim from calculate_seismicity_front_percentiles in
    hbi_analysis/notebooks/cooper_basin_plots-30.cleaned.ipynb, cell 49, called
    at cell 51 with bin_size=100, lower=90, upper=95: sort events by time, take
    them in bins of 100, and keep those whose distance lies between the 90th
    and 95th percentile OF THAT BIN. The kept (t, d) pairs ARE the front. It is
    a scatter, not a curve.

    An earlier version of this figure replaced it with a running maximum on the
    grounds that a percentile front can retreat. That was solving a problem the
    definition does not have -- the front is a cloud of points to be fitted, so
    a non-monotone sequence of points is not a defect -- and it silently
    changed the quantity every number was measured against.
    """
    o = np.argsort(t)
    t, d = np.asarray(t)[o], np.asarray(d)[o]
    ft, fd = [], []
    for i in range(0, len(t), bin_size):
        bt, bd = t[i:i + bin_size], d[i:i + bin_size]
        if not len(bt):
            continue
        lo, hi = np.percentile(bd, lower), np.percentile(bd, upper)
        m = (bd >= lo) & (bd <= hi)
        ft.extend(bt[m]); fd.extend(bd[m])
    return np.array(ft), np.array(fd)


def lam(x, R):
    """lambda in R = lambda*sqrt(x), least squares THROUGH THE ORIGIN.

    fit_sqrt_front from the same notebook, cell 55:
    lambda = sum(sqrt(x)*R) / sum(x).

    Through the origin is correct here and I had misread why. The observed
    front distances have origin_distance_km subtracted first -- the median of
    the first ten event distances -- so the initial cloud radius is already
    removed and R(0) = 0 is the right constraint, not an approximation. Fitting
    a free offset instead, as an earlier version did, re-solves a problem the
    subtraction has already solved and makes lambda incomparable with every
    lambda in the notebooks.
    """
    x, R = np.asarray(x, float), np.asarray(R, float)
    k = np.isfinite(x) & np.isfinite(R) & (x > 0)
    if k.sum() < 3:
        return np.nan
    b = np.sqrt(x[k])
    return float(np.sum(b * R[k]) / np.sum(b ** 2))


def volume_since_t0(obs):
    """(t_sim, cumulative ML) from the injection record, integrated from t0.

    q is L/s and dt is in days, so the trapezoidal integral is LITRES; 1 ML is
    1e6 L. Getting this wrong by 1e6 was a real bug in an earlier figure.
    """
    ti, q = obs["ti"], obs["q"]
    k = ti >= T0
    v = np.concatenate([[0.0], np.cumsum(np.diff(ti[k]) * 86400.0 * 0.5
                                        * (q[k][1:] + q[k][:-1]))]) / 1e6
    return ti[k] - T0, v


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", nargs="+", type=int, default=[1, 2, 3, 4])
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)

    obs = sf.observed()
    tcat, rcat, _ = fb.catalogue()
    kc = (tcat - T0) > 0
    tvol, vol = volume_since_t0(obs)
    Vobs = np.interp(tcat[kc] - T0, tvol, vol)

    for arm in a.arms:
        runs, title = ARMS[arm]
        runs = list(runs)
        cols = plt.cm.viridis(np.linspace(0.05, 0.85, len(runs)))

        fig, ax = plt.subplots(2, 2, figsize=(15.0, 9.6), dpi=200,
                               constrained_layout=True)
        (art, arv), (asz, adp) = ax

        # ---------------------------------------------------------- observed
        # The front is the project's own percentile definition, on the shifted
        # clock, with the initial cloud radius subtracted exactly as cell 51
        # does: origin = median of the first ten event distances. Plotted as
        # the SCATTER it is, with one lambda*sqrt(t) line through it.
        te, re_ = tcat[kc] - T0, rcat[kc]
        org = float(np.median(re_[np.argsort(te)][:10]))
        ft, fd = seismicity_front(te, re_ - org)
        Lo = lam(ft, fd)
        tf = np.linspace(0.0, TMAX, 300)
        art.scatter(ft, fd + org, s=13, color=OBSC, alpha=0.75, lw=0,
                    label=f"observed seismicity front ({len(ft)} points)")
        art.plot(tf, Lo * np.sqrt(tf) + org, "-", lw=2.4, color=OBSC,
                 label=f"   $\\lambda$ = {Lo:.1f} m/$\\sqrt{{d}}$")
        Vf = np.interp(ft, tvol, vol)
        LoV = lam(Vf, fd)
        vf = np.linspace(0.0, vol.max(), 300)
        arv.scatter(Vf, fd + org, s=13, color=OBSC, alpha=0.75, lw=0,
                    label="observed seismicity front")
        arv.plot(vf, LoV * np.sqrt(vf) + org, "-", lw=2.4, color=OBSC,
                 label=f"   $\\lambda_V$ = {LoV:.1f} m/$\\sqrt{{ML}}$")
        print(f"  observed front: {len(ft)} points, org {org:.0f} m, "
              f"lambda {Lo:.2f} m/sqrt(d), lambda_V {LoV:.2f} m/sqrt(ML)")
        asz.scatter(te, re_, s=3.0, alpha=0.18, color=MUTED, lw=0,
                    label=f"{int(kc.sum())} events after $t_0$")
        adp.plot(obs["tp"] - T0, obs["pm"] - P_REF, lw=1.2, color=MUTED,
                 alpha=0.85, label=r"measured $\Delta p$")

        print(f"\n{title}")
        print(f"  {'run':>7} {'tau_0':>6} {'front@8.7d':>11} {'dp bias':>9} "
              f"{'dp %':>7}")
        for c, n in zip(cols, runs):
            dk = sf.deck(n)
            d = sf.run_data(n, dk)
            if d is None:
                print(f"  {n}  no output")
                continue
            tau = _ff(dk["muinit"]) * _ff(dk["sigmainit"])
            lab = r"$\tau_0$ = " f"{tau:.2f} MPa"
            Rm = d["R"] * 1000.0                 # run_data returns KILOMETRES
            fr = dpb = dpp = np.nan
            if len(d["T"]) > 2:
                # SORT BY T -- run_data returns (T, R) indexed by cell, with
                # R = |x[i]|, so the raw arrays trace out and back. Ordering
                # also matters for np.interp; see score_cycle2.front_at.
                _o = np.argsort(d["T"])
                Tm, Rm = np.asarray(d["T"])[_o], np.asarray(d["R"])[_o] * 1000.0
                Vs = np.interp(Tm, tvol, vol)
                # The model slip front is the same kind of object as the
                # observed one -- one (first-exceedance time, |x|) point per
                # cell, calculate_slip_front in the notebook's cell 50 -- so it
                # is plotted the same way, as points with one lambda*sqrt line.
                Lr, LrV = lam(Tm, Rm), lam(Vs, Rm)
                art.scatter(Tm, Rm, s=9, color=c, alpha=0.7, lw=0,
                            label=lab + f"   |   $\\lambda$ = {Lr:.1f}"
                                  f"  ({Lr/Lo:.2f}$\\times$ observed)")
                art.plot(tf, Lr * np.sqrt(tf), "-", lw=1.4, color=c, alpha=0.9)
                arv.scatter(Vs, Rm, s=9, color=c, alpha=0.7, lw=0,
                            label=lab + f"   |   $\\lambda_V$ = {LrV:.1f}"
                                  f"  ({LrV/LoV:.2f}$\\times$)")
                arv.plot(vf, LrV * np.sqrt(vf), "-", lw=1.4, color=c, alpha=0.9)
                asz.plot(Tm, Rm, lw=1.9, color=c, label=lab)
                print(f"    {n}  lambda {Lr:7.2f} ({Lr/Lo:.2f}x obs)   "
                      f"lambda_V {LrV:7.2f} ({LrV/LoV:.2f}x)")
            else:
                # arm 3's lowest tau_0 never slips, so it has no front at all;
                # say so on the figure rather than leaving a silent gap
                art.plot([], [], lw=1.9, color=c, label=lab + "  (no slip)")
            if d.get("tpw") is not None:
                adp.plot(d["tpw"], d["ppw"] - P_STATIC, lw=1.9, color=c,
                         label=lab)
                hi = min(d["tpw"][-1], (obs["tp"] - T0).max(), 8.7)
                if hi > 0.2:
                    gr = np.linspace(0.05, hi, 2000)
                    ps = np.interp(gr, d["tpw"], d["ppw"]) - P_STATIC
                    ob = np.interp(gr, obs["tp"] - T0, obs["pm"]) - P_REF
                    qg = np.interp(gr, obs["ti"] - T0, obs["q"])
                    fl = (qg > 0.25 * np.nanmax(obs["q"])) & (ob > 5.0 - P_REF)
                    dpb = float(np.mean(ps[fl] - ob[fl]))
                    dpp = 100.0 * dpb / float(np.mean(ob[fl]))
            print(f"  {n:>7} {tau:6.2f} {fr:11.2f} {dpb:+9.2f} {dpp:+7.1f}")

        art.set(xlabel="Days since injection resumed (data-day 4.300)",
                ylabel="Front radius (m)", xlim=(0, TMAX), ylim=(0, 3000))
        art.set_title("(a)  Front radius vs time")
        art.legend(loc="upper left", framealpha=0.93)

        arv.set(xlabel="Cumulative injected volume since $t_0$ (ML)",
                ylabel="Front radius (m)", ylim=(0, 3000))
        arv.set_title("(b)  Front radius vs injected volume")
        arv.legend(loc="upper left", framealpha=0.93)

        asz.set(xlabel="Days since injection resumed",
                ylabel="Distance from injection point (m)",
                xlim=(0, TMAX), ylim=(0, 3000))
        asz.set_title("(c)  Seismicity, with the same model fronts")
        asz.legend(loc="upper left", framealpha=0.93)

        adp.set(xlabel="Days since injection resumed",
                ylabel="Wellhead pressure change (MPa)",
                xlim=(0, TMAX), ylim=(-2, 30))
        adp.set_title("(d)  Wellhead pressure change")
        adp.legend(loc="upper left", framealpha=0.93)

        fig.suptitle(title, fontsize=13)
        for e in ("png", "pdf"):
            fig.savefig(OUT / f"arm{arm}_compare.{e}", bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {OUT}/arm{arm}_compare.png")


if __name__ == "__main__":
    main()
