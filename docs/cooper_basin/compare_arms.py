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

  (a) front radius vs time. Model fronts are the ds-resolved slip contour at
      FRONT_THR = 1e-4 m; the observed front is the running maximum of event
      distance over the full catalogue, monotone by construction.

  (b) front radius vs cumulative injected volume since t0. Volume rather than
      time removes the rate steps and the shut-ins, which is what makes the
      raw r-t curve depart from sqrt(t) even when the physics is diffusive.

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
    runmax = np.maximum.accumulate(rcat)
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
        art.plot(tcat[kc] - T0, runmax[kc], lw=2.6, color=OBSC,
                 label="observed front (running max)")
        arv.plot(Vobs, runmax[kc], lw=2.6, color=OBSC, label="observed front")
        asz.scatter(tcat[kc] - T0, rcat[kc], s=3.0, alpha=0.18, color=MUTED,
                    lw=0, label=f"{int(kc.sum())} events after $t_0$")
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
                art.plot(d["T"], Rm, lw=1.9, color=c, label=lab)
                asz.plot(d["T"], Rm, lw=1.9, color=c, label=lab)
                arv.plot(np.interp(d["T"], tvol, vol), Rm, lw=1.9, color=c,
                         label=lab)
                if d["t_end"] >= 8.7:
                    fr = (float(np.interp(8.7, d["T"], d["R"])) * 1000.0
                          / float(np.interp(8.7 + T0, tcat, runmax)))
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
