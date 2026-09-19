#!/usr/bin/env python3
"""PAPER FIGURE: what the analytical crack says before any simulation is run.

(a) The stress-injection parameter T against the crack-growth factor
    lambda = a(t)/sqrt(4 alpha t), for the uniform-permeability solution
    (Saez et al. 2022) and the permeability-enhancement solution
    (Dunham 2024), with this project's own quadrature over the same relation
    as a check on the implementation.

(b) Both observables as model/observed against a UNIFORM fault permeability,
    evaluated at t = 8.7 d with the record's mean flowing rate.

WHAT (a) SHOWS. The two solutions are indistinguishable at small T -- a
critically stressed fault propagates a crack whether or not permeability
evolves, so enhancement is irrelevant there. They separate as T grows: at
T = 2 the enhanced crack grows 1.7x faster, at T = 5 it grows 5.2x faster, and
above T ~ 8 the uniform solution HAS NO SOLUTION AT ALL while the enhanced one
still returns a finite lambda. Enhancement matters precisely and only in the
understressed regime this paper is about.

WHAT (b) SHOWS, AND THE SUBTLETY THAT MAKES IT WORK. The overpressure scale
deltaP = q eta / (4 pi k) does not depend on which crack solution is used, so
THERE IS ONE PRESSURE CURVE AND TWO FRONT CURVES. The pressure pins k tightly:
it passes through the observed value near k = 5e-15 and is 4.7x too high at
1e-15 and 20x too low at 2.5e-13. The front is where the two models disagree,
and at the k the pressure selects the enhanced solution gives 0.93x against the
uniform solution's 0.80x.

WHAT IS *NOT* CLAIMED. An earlier version of this figure asserted that the
uniform solution cannot propagate a crack at Cooper Basin at all -- T = 27 to
257, no solution. That was wrong: it came from evaluating the
UNIFORM-permeability solution at kpmax, which is incoherent. Evaluated at the
permeability it actually assumes, the uniform solution is perfectly healthy
(lambda = 1.87 at k = 1e-15). The real result is the trade-off in (b), not an
impossibility, and it is weaker but true.

CAVEATS THAT BELONG IN THE CAPTION
  - the solution assumes a CONSTANT rate; the record is variable, so the mean
    flowing rate is used and T would span roughly 0.08 to 0.8 over the history
  - the solution is for a circular crack at nu = 0; HBI runs nu = 0.25
    (m_const.f90:4), which makes its crack elliptical
  - the analytical front peaks at 0.82x (uniform) and 1.06x (enhanced) rather
    than passing cleanly through 1, so this figure motivates the two-zone
    structure, it does not calibrate it

The curve labels come from the column names in Tlambda3D.csv and the legend of
2_revised_analytical_sol.ipynb (T3_SL = Saez, T3 = Dunham), not from reading
the two papers.

Usage:  python fig_analytical.py
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import exp1
from scipy.integrate import quad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = Path("/home/users/nberrios/3dhbi/hbi_analysis")
FIG = H / "figures" / "cycle2"
CSV = Path("/home/users/nberrios/3dhbi/hbi/Tlambda3D.csv")
IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")

ETA, PHI, BETA = 1.27e-4, 0.01, 2.25e-8
F0, SIG0, TAU0 = 0.60, 27.99e6, 11.53e6      # res633001.in
T_EVAL_D = 8.7
R_OBS, DP_OBS = 763.0, 11.7                  # m, MPa

INK = "#141413"
C_UNIF, C_ENH, C_PRESS = "#8E44AD", "#1d4ed8", "#0e7a63"

plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 9, "legend.fontsize": 7.8,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "axes.edgecolor": INK, "axes.linewidth": 0.7,
    "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": INK, "ytick.color": INK,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 3.0, "ytick.major.size": 3.0,
    "axes.grid": False, "legend.frameon": False,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def panel(ax, letter):
    ax.text(-0.19, 1.03, f"({letter})", transform=ax.transAxes,
            fontsize=10, fontweight="bold", va="bottom", ha="left")


def T_quad(lam, eps=0.0):
    """T(lambda) by direct quadrature -- the implementation check for (a)."""
    def g(xi):
        return exp1(lam**2 * xi**2) - exp1(lam**2) + lam**-2*np.exp(-lam**2) - eps
    return quad(lambda xi: g(xi) * xi / np.sqrt(1 - xi**2), 0, 1)[0]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="fig_analytical")
    a = ap.parse_args(argv)

    df = pd.read_csv(CSV).sort_values("lambda")
    lam_q = np.logspace(-2, 1.0, 60)
    T_q = np.array([T_quad(l) for l in lam_q])

    def lam_at(T, col):
        d = df[["lambda", col]].dropna().sort_values(col)
        if not (d[col].min() <= T <= d[col].max()):
            return np.nan
        return float(np.interp(T, d[col], d["lambda"]))

    q = np.array(open(IN / "june_clean_from_d4300.txt").read()
                 .splitlines()[4].split(), float)
    qbar = float(np.mean(q[q > 0]))
    strength = F0 * SIG0 - TAU0
    t = T_EVAL_D * 86400.0

    ks = np.logspace(np.log10(4e-16), np.log10(4e-13), 160)
    alpha = ks / (ETA * PHI * BETA)
    dP = qbar * ETA / (4 * np.pi * ks)
    Ts = strength / (F0 * dP)
    fr = {c: np.array([lam_at(T, c) for T in Ts]) * np.sqrt(4 * alpha * t) / R_OBS
          for c in ("T3_SL", "T3")}
    pr = dP / 1e6 / DP_OBS

    fig, (aa, ab) = plt.subplots(1, 2, figsize=(7.48, 3.25), dpi=400,
                                 constrained_layout=True)

    aa.loglog(df["lambda"], df["T3_SL"], "-", lw=1.8, color=C_UNIF,
              label="uniform permeability (Sáez et al. 2022)")
    aa.loglog(df["lambda"], df["T3"], "-", lw=1.8, color=C_ENH,
              label="permeability enhancement (Dunham 2024)")
    aa.loglog(lam_q, T_q, ls=(0, (2, 2)), lw=1.2, color=INK,
              label="this implementation")
    aa.set(xlabel=r"crack-growth factor  $\lambda = a(t)/\sqrt{4\alpha t}$",
           ylabel="stress-injection parameter, $T$",
           xlim=(1e-2, 10), ylim=(1e-1, 1e2))
    aa.legend(loc="lower left", handlelength=1.9, labelspacing=0.35)

    ab.axhline(1.0, color=INK, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ab.loglog(ks, fr["T3_SL"], "-", lw=1.8, color=C_UNIF, zorder=3,
              label="front, uniform permeability")
    ab.loglog(ks, fr["T3"], "-", lw=1.8, color=C_ENH, zorder=3,
              label="front, with enhancement")
    ab.loglog(ks, pr, "-", lw=1.8, color=C_PRESS, zorder=3,
              label="wellhead pressure (both)")
    ab.set(xlabel=r"uniform fault permeability, $k$ (m$^2$)",
           ylabel="model / observed", xlim=(ks[0], ks[-1]), ylim=(0.02, 20))
    ab.legend(loc="lower left", handlelength=1.9, labelspacing=0.35)

    for x, l in ((aa, "a"), (ab, "b")):
        panel(x, l)
    for e in ("png", "pdf"):
        fig.savefig(FIG / f"{a.out}.{e}")
    plt.close(fig)

    # --- the numbers the caption needs
    med = np.median(np.abs(np.interp(df["lambda"], lam_q, T_q) - df["T3"])
                    / df["T3"])
    print(f"implementation check vs Dunham 2024: median |rel diff| {med*100:.2f}%")
    kp1 = float(np.interp(0.0, np.log(pr[::-1]), np.log(ks[::-1])))
    print(f"pressure matches at k = {np.exp(kp1):.2e} m^2")
    for c, nm in (("T3_SL", "uniform"), ("T3", "enhanced")):
        v = fr[c][np.isfinite(fr[c])]
        kk = ks[np.isfinite(fr[c])]
        print(f"  front, {nm:8s}: peaks {v.max():.2f}x at k = "
              f"{kk[v.argmax()]:.1e}; at the pressure-matching k it is "
              f"{np.interp(np.exp(kp1), ks, fr[c]):.2f}x"
              + (f"; no crack above k = {kk.max():.1e}" if kk.max() < ks[-1]
                 else ""))
    print(f"wrote {FIG}/{a.out}.png")


if __name__ == "__main__":
    main()
