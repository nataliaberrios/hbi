#!/usr/bin/env python3
"""The front constrains a STRENGTH MARGIN, not a stress state.

The argument this figure makes, in three panels:

  A  wellhead pressure. The measured record fixes the pressure amplitude
     A = Q0*eta/(4*pi*k*w) = 0.572 MPa. Both candidate models share it, so it is
     not what distinguishes them.
  B  front radius against time, from Wang & Dunham's OWN closed-form criterion
     (source_code/cmp_seis_extent.m): Dp(r,t) = A*E1(phi*eta*beta*r^2/(4*k*t)),
     front where Dp = Dtauc = f0*sigmabar_0 - tau_0. Three curves: their
     calibrated tau_0 = 15.0 at lab friction (matches), the resolved
     tau_0 = 10.26 at lab friction (7 m at 5 d, sixty times too small), and the
     resolved stress at f0 = 0.433 (matches). The data cannot tell the first
     from the third.
  C  the Dtauc = 1.86 MPa locus in (tau_0, f0). Both models sit on it. Shaded:
     the published friction range for chlorite-bearing granitic gouge, and the
     resolved-stress band. The figure's point is that the observation pins the
     COMBINATION, and the stress measurement survives if f0 is near 0.43.

WHY THIS RATHER THAN "TAIYI DOES NOT MATCH THE SLIP". He does match the front --
196/340/438 m against an observed 187/323/417, about 5%. A figure claiming
otherwise collapses the moment anyone opens cmp_seis_extent.m. What is true is
that his framework predicts no slip AMPLITUDE on the main fault at all, the
seismicity coming from separate off-fault spring-sliders, so there is no
cumulative-slip curve to compare against the observed 2.81 -> 9.17 cm. That is
"does not predict", not "fails to match", and it is a footnote here rather than
the headline.

Dp falls off logarithmically in r, so the front radius depends EXPONENTIALLY on
Dtauc -- a 0.2 MPa change in tau_0 moves the 5 d front by 20%. That sharpness is
why the locus is narrow and why it is worth plotting as a locus.

Usage:  python make_strength_margin_figure.py
"""
import importlib.util as iu
import sys
from pathlib import Path

import numpy as np
from scipy.special import exp1
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = "/home/users/nberrios/3dhbi/hbi_analysis"
_s = iu.spec_from_file_location("sf", H + "/make_sweep_figures.py")
sf = iu.module_from_spec(_s); _s.loader.exec_module(sf)

OUT = Path(H) / "figures" / "strength_margin"
# Wang & Dunham Table 1
Q0, ETA, K, W, PHI, BETA = 20e-3, 8.9e-4, 4e-13, 6.0, 0.01, 1e-8
SIG, DZ, STD = 28.0e6, 0.5, 2.0
A = Q0 * ETA / (4 * np.pi * K * W) * np.exp(-DZ ** 2 / (2 * STD ** 2))
TAU_RES, F0_ALT, F0_LAB, TAU_WD = 10.26, 0.433, 0.60, 15.00
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 10.5,
                     "axes.edgecolor": MUTED, "text.color": INK,
                     "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "legend.fontsize": 8.5})


def front(t_s, dtauc_pa):
    f = lambda r: A * exp1(PHI * ETA * BETA * r ** 2 / (4 * K * t_s)) - dtauc_pa
    try:
        return brentq(f, 1e-4, 8000.0)
    except ValueError:
        return np.nan


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    obs = sf.observed()
    fig, (a0, a1, a2) = plt.subplots(1, 3, figsize=(15.5, 4.5), dpi=200,
                                     constrained_layout=True)

    # --- A: the wellhead pins the amplitude
    a0.plot(obs["tp"], obs["pm"], lw=1.3, color=MUTED,
            label="measured wellhead (Habanero 4)")
    dk = sf.deck(632881); d = sf.run_data(632881, dk)
    if d is not None and d.get("tpw") is not None:
        a0.plot(d["tpw"], d["ppw"], lw=2, color="#1d4ed8",
                label="Wang & Dunham parameters in HBI (632881), +8.2%")
    a0.set(xlabel="Days since injection began",
           ylabel="Absolute wellhead pressure (MPa)", xlim=(0, 5))
    a0.set_title("A. The wellhead fixes the pressure amplitude\n"
                 f"$A = Q_0\\eta/4\\pi k w$ = {A/1e6:.3f} MPa — shared by both models")
    a0.legend(loc="lower right"); a0.grid(alpha=0.3, color=GRID)

    # --- B: the front, from their own criterion
    t = np.linspace(0.15, 6.0, 260)
    cases = [(TAU_WD, F0_LAB, "#1d4ed8",
              f"Wang & Dunham: $\\tau_0$={TAU_WD:.2f}, $f_0$={F0_LAB:.2f}"),
             (TAU_RES, F0_LAB, "#a8071a",
              f"resolved stress at lab friction: $\\tau_0$={TAU_RES:.2f}, "
              f"$f_0$={F0_LAB:.2f}"),
             (TAU_RES, F0_ALT, "#009E73",
              f"resolved stress at $f_0$={F0_ALT:.3f}: $\\tau_0$={TAU_RES:.2f}")]
    for tau0, f0, c, lb in cases:
        dt = (f0 * SIG / 1e6 - tau0) * 1e6
        r = np.array([front(x * 86400, dt) for x in t])
        a1.plot(t, r, lw=2.2, color=c,
                label=lb + f"   ($\\Delta\\tau_c$={dt/1e6:.2f} MPa)")
    # The front distances carry a separate ORIGIN offset (obs["org"], 90 m --
    # the injector's own offset in the catalogue frame). Omitting it plots the
    # observed front 90 m too shallow: 323 m instead of 413 m at 5 d, which
    # made the curves look 40% high when they are within 6%.
    mo = obs["ft"] <= 6
    a1.scatter(obs["ft"][mo], (obs["fd"][mo] + obs["org"]) * 1000.0,
               s=13, color=MUTED, zorder=1, label="observed seismicity front")
    a1.set(xlabel="Days since injection began",
           ylabel="Front radius (m)", xlim=(0, 6), ylim=(0, 700))
    a1.set_title("B. Their own front criterion.\nThe data cannot distinguish "
                 "blue from green.")
    a1.legend(loc="upper left"); a1.grid(alpha=0.3, color=GRID)

    # --- C: the locus
    dt_req = (F0_ALT * SIG / 1e6 - TAU_RES)
    f0g = np.linspace(0.35, 0.70, 300)
    a2.plot(f0g * SIG / 1e6 - dt_req, f0g, lw=2.4, color=INK,
            label=f"$\\Delta\\tau_c$ = {dt_req:.2f} MPa (what the front requires)")
    a2.axhspan(0.37, 0.50, color="#009E73", alpha=0.16,
               label="chlorite-bearing granitic gouge, $f_0$ 0.37–0.50")
    a2.axhspan(0.60, 0.71, color="#1d4ed8", alpha=0.13,
               label="unaltered granite / feldspar, $f_0$ 0.60–0.71")
    a2.axvspan(10.26, 10.36, color="#a8071a", alpha=0.20,
               label="resolved $\\tau_0$ from the stress measurements")
    for tau0, f0, c, nm in ((TAU_WD, F0_LAB, "#1d4ed8", "Wang & Dunham"),
                            (TAU_RES, F0_ALT, "#009E73", "this model")):
        a2.scatter([tau0], [f0], s=150, color=c, edgecolors=INK, zorder=5)
        a2.annotate(nm, (tau0, f0), fontsize=9, color=c, fontweight="bold",
                    textcoords="offset points", xytext=(-10, 12), ha="center")
    a2.set(xlabel=r"$\tau_0$ (MPa)", ylabel="$f_0$", xlim=(9, 17),
           ylim=(0.35, 0.72))
    a2.set_title("C. The front pins the COMBINATION\n"
                 "$\\Delta\\tau_c = f_0\\bar\\sigma_0 - \\tau_0$, not $\\tau_0$")
    a2.legend(loc="lower right", fontsize=7.8); a2.grid(alpha=0.3, color=GRID)

    for e in ("png", "pdf"):
        fig.savefig(OUT / f"strength_margin.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/strength_margin.png")
    print(f"\nA = {A/1e6:.3f} MPa;  required Dtauc = {dt_req:.2f} MPa")
    print(f"{'route':>30s} {'tau_0':>6s} {'f0':>6s} {'dtauc':>7s} "
          f"{'1d':>6s} {'3d':>6s} {'5d':>6s}")
    print(f"{'OBSERVED':>30s} {'':>6s} {'':>6s} {'':>7s} "
          f"{187:>5d}m {323:>5d}m {417:>5d}m")
    for tau0, f0, c, lb in cases:
        dt = (f0 * SIG / 1e6 - tau0) * 1e6
        print(f"{lb.split(':')[0][:30]:>30s} {tau0:>6.2f} {f0:>6.3f} "
              f"{dt/1e6:>6.2f}M "
              + " ".join(f"{front(x*86400, dt):>5.0f}m" for x in (1, 3, 5)))


if __name__ == "__main__":
    main()
