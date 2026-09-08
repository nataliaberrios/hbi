#!/usr/bin/env python3
"""Cycle-2 comparison figures: pressure, R-T, R-V, slip, on the SHIFTED clock.

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
  R-T       front radius vs time. Observed front is the RUNNING MAXIMUM of the
            full 20735-event catalogue, which is monotone by construction; a
            per-bin percentile is not a front and retreats four times on this
            catalogue.
  R-V       the same front against cumulative injected volume SINCE t0, which
            removes the shut-in/rate-step staircase that makes R-T non-sqrt(t).
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

    fig, ax = plt.subplots(2, 2, figsize=(15.0, 9.4), dpi=200,
                           constrained_layout=True)
    (ap_, ar), (av, asl) = ax

    # --- observed, all pulled back to sim time
    ap_.plot(obs["tp"] - T0, obs["pm"], lw=1.2, color=MEAS, alpha=0.85,
             label="measured wellhead (shifted)")
    m = (tcat - T0) > 0
    ar.plot(tcat[m] - T0, runmax[m], lw=2.4, color="#a8071a",
            label="observed front (running max)")
    tg = np.linspace(0.02, 13.1, 300)
    asl.plot(tg, np.interp(tg + T0, OT, ocm) - slip_t0, lw=2.4, color="#a8071a",
             label="observed slip INCREMENT")
    # cumulative volume since t0
    ti, q = obs["ti"], obs["q"]
    kv = ti >= T0
    # q is L/s and dt is in days, so the integral is LITRES; 1 ML = 1e6 L.
    vol = np.concatenate([[0.0], np.cumsum(np.diff(ti[kv]) * 86400.0
                                           * 0.5 * (q[kv][1:] + q[kv][:-1]))]) / 1e6
    tvol = ti[kv] - T0
    av.plot(np.interp(tcat[m] - T0, tvol, vol), runmax[m], lw=2.4,
            color="#a8071a", label="observed front")

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
        ar.plot(d["T"], d["R"] * 1000.0, lw=1.7, color=c, label=label(n))
        av.plot(np.interp(d["T"], tvol, vol), d["R"] * 1000.0, lw=1.7, color=c)
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
    ap_.set_title("Wellhead pressure. Observed pulled back 4.300 d;\n"
                  "the 87.8 MPa water-hammer transient is off-scale by design")
    ap_.legend(loc="lower right", fontsize=7.5); ap_.grid(alpha=0.3, color=GRID)

    ar.set(xlabel="Days since injection resumed", ylabel="Front radius (m)",
           xlim=(0, 13.2), ylim=(0, 2600))
    ar.set_title("R–T. Observed is a RUNNING MAXIMUM (monotone);\n"
                 "no $\\sqrt{t}$ fit — neither front starts at the origin")
    ar.legend(loc="upper left", fontsize=7.5); ar.grid(alpha=0.3, color=GRID)

    av.set(xlabel="Cumulative injected volume since t$_0$ (ML)",
           ylabel="Front radius (m)", ylim=(0, 2600))
    av.set_title("R–V. Volume removes the rate-step staircase\n"
                 "that makes R–T non-$\\sqrt{t}$")
    av.legend(loc="upper left", fontsize=7.5)
    av.grid(alpha=0.3, color=GRID)

    asl.set(xlabel="Days since injection resumed",
            ylabel="Slip at the injector (cm)", xlim=(0, 13.2))
    asl.set_title("Slip vs the observed INCREMENT since t$_0$.\n"
                  "HBI zeroes slip at restart, so absolute slip is not the target")
    asl.legend(loc="upper left", fontsize=7.5); asl.grid(alpha=0.3, color=GRID)

    for e in ("png", "pdf"):
        fig.savefig(OUT / f"cycle2_compare_{a.tag}.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {OUT}/cycle2_compare_{a.tag}.png")


if __name__ == "__main__":
    main()
