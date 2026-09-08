#!/usr/bin/env python3
"""Cycle-2 comparison figures: pressure, R-T, R-V, slip, on the SHIFTED clock.

PANEL TITLES ARE PLAIN DESCRIPTIONS -- (a) Wellhead pressure, and so on. Every
justification for a transformation lives in the notebook captions
(notebooks/cycle2_comparison.ipynb), not on the figure. An earlier version
argued its own methodology in the titles, which reads as covering for something
and is not what a title is for.

compare_runs_pressure_RT_RV.py cannot be used for these runs. It reads observed
series at ABSOLUTE time and fits lambda*sqrt(t) through the origin, and neither
holds after a restart: sim t = 0 is data-day 4.300, and the observed front is
already ~310 m there while the simulated front starts at the 300 m disc edge, so
neither passes through the origin.

Everything here is therefore on SIM TIME, with every observed series pulled back
by 4.300 d, and NO sqrt(t) fit anywhere.

FOUR PANELS
  pressure  wellhead, sim vs observed-shifted. The flowing mask is the usual
            q > 25% of peak and p_obs > 5 MPa, PLUS p_obs < 60 MPa to drop the
            two-sample 87.76 MPa transient at data-day 14.23 (92% of Sv,
            almost certainly water-hammer -- masking it, not fitting it).
  R-T       front radius vs time, with a SQUARE-ROOT FIT overlaid (dashed).
            Observed front is the RUNNING MAXIMUM of the full 20735-event
            catalogue, monotone by construction; a per-bin percentile is not a
            front and retreats four times on this catalogue.
  R-V       the same, against cumulative injected volume SINCE t0, which removes
            the shut-in/rate-step staircase that makes R-T non-sqrt(t). Also
            square-root fitted.

THE SQUARE-ROOT FIT, AND WHY IT IS NOT THROUGH THE ORIGIN. The classic form is
r = sqrt(4 pi D t), which forces r(0) = 0. That is wrong here: at sim t = 0 the
observed front is already ~420 m and the simulated front starts at the 300 m
disc edge, because both have a head start from cycle 1. Forcing the curve
through the origin would absorb that offset into D and bias it.

The fit used keeps the square-root SHAPE and adds the head start explicitly:

    r(t)^2 = r0^2 + 4 pi D t          <=>    r = sqrt(r0^2 + 4 pi D t)

which is exactly r = sqrt(4 pi D (t + t_off)) with r0^2 = 4 pi D t_off -- the
same square-root law, started earlier. Because it is linear in t, it is fitted
by ordinary least squares on r^2 against t, and the R^2 of THAT regression is
the honest test of whether the square-root law holds at all. Both D and R^2 are
printed and put in the legend. For R-V the identical algebra applies with V in
place of t, since V is proportional to t at constant rate.
  slip      slip at the injector against the observed INCREMENT,
            obs(t + 4.3 d) - obs(4.3 d) = 2.598 cm by data-day 13. Required
            because main_LH.f90:923 sets slip = 0d0 -- HBI cannot start with
            pre-existing slip, so absolute observed slip is not the target.

Usage:  python compare_cycle2.py --tag top5 632963 632964 632962 632980 632981
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
from sim_curves import load_slip, _deck, _ff
_f = iu.spec_from_file_location(
    "fb", "/home/users/nberrios/3dhbi/hbi_git/docs/postshutin/front_backfront.py")
fb = iu.module_from_spec(_f); _f.loader.exec_module(fb)

T0 = 4.300
OUT = Path(H) / "figures" / "cycle2"
OBS_SLIP = Path("/home/users/nberrios/3dhbi/hbi/slip_profiles_strike.txt")
OT = [3, 5, 7, 9, 11, 13, 15, 17]
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
MEAS = "#6b6b66"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 10.5,
                     "axes.edgecolor": MUTED, "text.color": INK,
                     "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "legend.fontsize": 8})


def arm_of(n):
    d = _deck(n)
    e, b = _ff(d["eta"]), _ff(d["beta"])
    if abs(e - 0.89e-3) < 1e-9:
        return 1, r"$\eta$ 0.89e-3, $\beta$ 2.25e-8"
    if b > 1e-7:
        return 3, r"$\eta$ 1.27e-4, $\beta$ 1.58e-7"
    if b < 1e-8:
        return 4, r"$\eta$ 1.27e-4, $\beta$ 7.72e-9"
    return 2, r"$\eta$ 1.27e-4, $\beta$ 2.25e-8"


def sqrt_fit(x, r):
    """r^2 = b*(x - x_off) by OLS on r^2 against x. Returns (x_off, b, R^2).

    r = sqrt(b*(x - x_off)) is the square-root law with its origin SHIFTED, and
    because r^2 is linear in x it is an ordinary least-squares fit whose R^2 is
    the honest test of whether the law holds at all.

    x_off is reported signed and NOT clamped. An earlier version clamped a
    head-start term to zero, which hid the result: every fit here prefers
    x_off > 0, i.e. a DELAY, because both fronts grow more slowly than sqrt(t)
    early on -- the observed front sits at ~420 m for about a day after
    injection resumes, and the simulated fronts are zero until ~1.8 d because
    no cell has yet slipped past the 1e-4 m threshold.
    """
    k = np.isfinite(x) & np.isfinite(r) & (r > 0)
    x, r = np.asarray(x)[k], np.asarray(r)[k]
    if len(x) < 5:
        return np.nan, np.nan, np.nan
    A = np.vstack([x, np.ones_like(x)]).T
    (b, c), res, *_ = np.linalg.lstsq(A, r ** 2, rcond=None)
    ss = np.sum((r ** 2 - np.mean(r ** 2)) ** 2)
    r2 = float(1 - (res[0] / ss)) if len(res) and ss > 0 else np.nan
    return float(-c / b) if b else np.nan, float(b), r2


def label(n):
    a, txt = arm_of(n)
    d = _deck(n)
    return (f"{n}  arm {a}: {txt}, "
            r"$\tau_0$ = " f"{_ff(d['muinit'])*_ff(d['sigmainit']):.2f} MPa")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jobs", nargs="+", type=int)
    ap.add_argument("--tag", default="cycle2")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    obs = sf.observed()
    tcat, rcat, _ = fb.catalogue()
    runmax = np.maximum.accumulate(rcat)
    o = np.loadtxt(OBS_SLIP)
    ocm = [o[:, 1 + i].max() for i in range(8)]
    slip_t0 = float(np.interp(T0, OT, ocm))
    cols = plt.cm.viridis(np.linspace(0.05, 0.85, len(a.jobs)))

    fig, ax = plt.subplots(2, 2, figsize=(16.0, 10.0), dpi=200,
                           constrained_layout=True)
    (ap_, ar), (av, asl) = ax

    # --- observed, all pulled back to sim time
    ap_.plot(obs["tp"] - T0, obs["pm"], lw=1.2, color=MEAS, alpha=0.85,
             label="measured")
    m = (tcat - T0) > 0
    ar.plot(tcat[m] - T0, runmax[m], lw=2.4, color="#a8071a",
            label="observed")
    xo, bo, r2o = sqrt_fit(tcat[m] - T0, runmax[m])
    Do = bo / (4 * np.pi * 86400.0)
    tf = np.linspace(0.02, 13.1, 300)
    ar.plot(tf, np.sqrt(np.maximum(bo * (tf - xo), 0)), "--", lw=2.2,
            color="#a8071a", alpha=0.8,
            label=r"  fit $\sqrt{4\pi D(t-t_{off})}$: $D$="
                  f"{Do:.3f}, $t_{{off}}$={xo:+.2f} d, $R^2$={r2o:.3f}")
    print(f"observed  R-T: D {Do:.4f} m2/s, t_off {xo:+.2f} d, R2 {r2o:.3f}")

    tg = np.linspace(0.02, 13.1, 300)
    asl.plot(tg, np.interp(tg + T0, OT, ocm) - slip_t0, lw=2.4, color="#a8071a",
             label="observed (increment)")
    # cumulative volume since t0
    ti, q = obs["ti"], obs["q"]
    kv = ti >= T0
    # q is L/s and dt is in days, so the integral is LITRES; 1 ML = 1e6 L.
    vol = np.concatenate([[0.0], np.cumsum(np.diff(ti[kv]) * 86400.0
                                           * 0.5 * (q[kv][1:] + q[kv][:-1]))]) / 1e6
    tvol = ti[kv] - T0
    Vo = np.interp(tcat[m] - T0, tvol, vol)
    av.plot(Vo, runmax[m], lw=2.4, color="#a8071a", label="observed")
    xv, bv, r2v = sqrt_fit(Vo, runmax[m])
    vf = np.linspace(0.05, vol.max(), 300)
    av.plot(vf, np.sqrt(np.maximum(bv * (vf - xv), 0)), "--", lw=2.2,
            color="#a8071a", alpha=0.8,
            label=r"  fit $\sqrt{c(V-V_{off})}$: "
                  f"$V_{{off}}$={xv:+.1f} ML, $R^2$={r2v:.3f}")
    print(f"observed  R-V: V_off {xv:+.2f} ML, R2 {r2v:.3f}")

    print(f"{'run':>7s} {'arm':>4s} {'tau_0':>6s} {'front':>7s} {'wellhd':>8s} "
          f"{'slip':>7s}")
    for c, n in zip(cols, a.jobs):
        dk = sf.deck(n); d = sf.run_data(n, dk)
        if d is None:
            print(f"  [skip] {n}: no output"); continue
        if d.get("tpw") is not None:
            ap_.plot(d["tpw"], d["ppw"], lw=1.7, color=c, label=label(n))
        # run_data returns R in KILOMETRES -- its lambda fit works in km.
        # Plotting it straight on a metre axis gives a flat line at ~0.
        Rm = d["R"] * 1000.0
        ar.plot(d["T"], Rm, lw=1.7, color=c, label=label(n))
        xo_, b, r2 = sqrt_fit(d["T"], Rm)
        D = b / (4 * np.pi * 86400.0)
        ar.plot(tf, np.sqrt(np.maximum(b * (tf - xo_), 0)), "--", lw=1.3,
                color=c, alpha=0.85)
        # fit parameters folded into the RUN label, so each run costs one legend
        # entry rather than two -- with 10 runs the doubled legend covered the
        # curves entirely.
        ar.lines[-2].set_label(
            f"{n} arm {arm_of(n)[0]}, " r"$\tau_0$=" f"{_ff(_deck(n)['muinit'])*_ff(_deck(n)['sigmainit']):.2f}"
            f" | $D$={D:.3f}, $t_{{off}}$={xo_:+.2f} d, $R^2$={r2:.3f}")
        Vs = np.interp(d["T"], tvol, vol)
        av.plot(Vs, Rm, lw=1.7, color=c)
        xv_, bb, r2b = sqrt_fit(Vs, Rm)
        av.plot(vf, np.sqrt(np.maximum(bb * (vf - xv_), 0)), "--", lw=1.3,
                color=c, alpha=0.85)
        print(f"  {n}  R-T: D {D:.4f}, t_off {xo_:+.2f} d, R2 {r2:.3f}   "
              f"R-V: V_off {xv_:+.2f} ML, R2 {r2b:.3f}")
        ts = np.linspace(0.2, min(13.1, d["t_end"] if "t_end" in d else 13.1), 40)
        sv = []
        for t_ in ts:
            try:
                _, sl, ta = load_slip(n, t_, how="strike")
                sv.append(sl[0] * 100 if abs(ta - t_) < 0.3 else np.nan)
            except Exception:
                sv.append(np.nan)
        asl.plot(ts, sv, lw=1.7, color=c)

    ap_.set(xlabel="Days since injection resumed (data-day 4.300)",
            ylabel="Absolute wellhead pressure (MPa)", xlim=(0, 13.2),
            ylim=(25, 62))
    ap_.set_title("(a)  Wellhead pressure")
    ap_.legend(loc="lower right", fontsize=7.0, framealpha=0.93,
               labelspacing=0.35)
    ap_.grid(alpha=0.3, color=GRID)

    ar.set(xlabel="Days since injection resumed", ylabel="Front radius (m)",
           xlim=(0, 13.2), ylim=(0, 2600))
    ar.set_title("(b)  Front radius vs time")
    ar.legend(loc="lower right", fontsize=7.0, framealpha=0.93,
              borderpad=0.5, labelspacing=0.35)
    ar.grid(alpha=0.3, color=GRID)

    av.set(xlabel="Cumulative injected volume since t$_0$ (ML)",
           ylabel="Front radius (m)", ylim=(0, 2600))
    av.set_title("(c)  Front radius vs cumulative injected volume")
    av.legend(loc="lower right", fontsize=7.0, framealpha=0.93)
    av.grid(alpha=0.3, color=GRID)

    asl.set(xlabel="Days since injection resumed",
            ylabel="Slip at the injector (cm)", xlim=(0, 13.2))
    asl.set_title("(d)  Cumulative slip at the injector")
    asl.legend(loc="upper left", fontsize=7.5, framealpha=0.93)
    asl.grid(alpha=0.3, color=GRID)

    for e in ("png", "pdf"):
        fig.savefig(OUT / f"cycle2_compare_{a.tag}.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {OUT}/cycle2_compare_{a.tag}.png")


if __name__ == "__main__":
    main()
