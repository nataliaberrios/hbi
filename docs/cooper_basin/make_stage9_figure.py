#!/usr/bin/env python3
"""Stage 9: what kpmax actually controls, and why it cannot deliver both targets.

Three panels:

  left    pressure profiles at 5 d for the four kpmax values. Shows the
          structure directly -- a plateau with a logarithmic peak on top -- and
          that BOTH fall as kpmax rises. The peak falls as 1/kpmax, which is the
          prediction that motivated the sweep; the plateau falling with it is
          what was not predicted and is what binds.
  middle  the two decompositions against kpmax on log axes, with the 1/kpmax
          reference line through the peak, plus dp_crit as a horizontal line.
          Where the plateau crosses dp_crit is where slip switches off.
  right   peak slip against plateau - dp_crit, with the observed 2.81 cm at 5 d
          marked. This is the relation Stage 10 exploits: to move right along
          this axis at fixed kpmax, lower dp_crit by raising tau_0.

Usage:  python make_stage9_figure.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/home/users/nberrios/3dhbi/hbi_analysis/notebooks")
from sim_curves import load_dp, load_slip, _deck, _ff

OUT = Path("/home/users/nberrios/3dhbi/hbi_analysis/figures/stage9")
SCORES = Path("/home/users/nberrios/3dhbi/hbi_git/docs/figs/cooper_grid/grid_scores.json")
RUNS = [632875, 632896, 632897, 632898]
OBS_SLIP_CM = 2.81      # observed cumulative slip at the injector, 5 d
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 10.5,
                     "axes.edgecolor": MUTED, "text.color": INK,
                     "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "legend.fontsize": 8.5})


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sc = {x["n"]: x for x in json.load(open(SCORES))}
    rows = []
    for n in RUNS:
        d = _deck(n)
        sig, f0, mu0 = _ff(d["sigmainit"]), _ff(d["f0"]), _ff(d["muinit"])
        dpc = sig - mu0 * sig / f0
        r, dp, ta = load_dp(n, 5.0, how="strike")
        rs, sl, _ = load_slip(n, 5.0, how="strike")
        b = np.where(sl < 1e-4)[0]
        R = rs[b[0]] if len(b) else rs[-1]
        m = (r >= 50) & (r <= max(0.7 * R, 60))
        rows.append(dict(n=n, kpmax=_ff(d["kpmax"]), r=r, dp=dp, rs=rs, sl=sl,
                         R=R, plat=np.median(dp[m]) / 1e6, wl=dp[0] / 1e6,
                         dpc=dpc, pk=sl[0] * 100.0,
                         wh=sc[n]["p_pct"], fr=sc[n]["lam_ratio"]))
    dpc = rows[0]["dpc"]
    cols = plt.cm.viridis(np.linspace(0.05, 0.85, len(rows)))

    fig, (a0, a1, a2) = plt.subplots(1, 3, figsize=(15.0, 4.3), dpi=200,
                                     constrained_layout=True)

    for w, c in zip(rows, cols):
        a0.plot(w["r"], w["dp"] / 1e6, lw=1.9, color=c,
                label=f"kpmax {w['kpmax']:.1e}")
    a0.axhline(dpc, color=INK, ls="--", lw=1.2)
    a0.text(1050, dpc + 0.25, f"$\\Delta p_{{crit}}$ = {dpc:.2f} MPa",
            fontsize=8.5, color=INK)
    a0.set(xlabel="Distance along strike from injector (m)",
           ylabel="Overpressure $\\Delta p$ (MPa)", xlim=(0, 1400))
    a0.set_title("A plateau with a peak on top — and\nBOTH fall as kpmax rises")
    a0.legend(loc="upper right")
    a0.grid(alpha=0.3, color=GRID)

    k = np.array([w["kpmax"] for w in rows])
    plat = np.array([w["plat"] for w in rows])
    peak = np.array([w["wl"] - w["plat"] for w in rows])
    a1.loglog(k, peak, "o-", lw=1.9, ms=7, color="#a8071a",
              label="peak above plateau")
    a1.loglog(k, peak[0] * k[0] / k, ":", lw=1.4, color="#a8071a",
              label="$1/k_{pmax}$ reference")
    a1.loglog(k, plat, "s-", lw=1.9, ms=7, color="#1d4ed8", label="plateau")
    a1.axhline(dpc, color=INK, ls="--", lw=1.2, label="$\\Delta p_{crit}$")
    a1.set(xlabel="kpmax (m$^2$)", ylabel="MPa")
    a1.set_title("The peak obeys $1/k_{pmax}$ as predicted.\n"
                 "The plateau falls too — that is the trap.")
    a1.legend(loc="lower left")
    a1.grid(alpha=0.3, which="both", color=GRID)

    exc = plat - dpc
    pk = np.array([w["pk"] for w in rows])
    a2.plot(exc, pk, "o-", lw=1.9, ms=8, color="#1a1a19", zorder=3)
    for w, e, p in zip(rows, exc, pk):
        a2.annotate(f"{w['n']}\n{w['wh']:+.1f}%", (e, p),
                    textcoords="offset points", xytext=(8, -4), fontsize=8,
                    color=MUTED)
    a2.axhline(OBS_SLIP_CM, color="#a8071a", ls="--", lw=1.4)
    a2.text(-2.0, OBS_SLIP_CM + 0.2, f"observed {OBS_SLIP_CM:.2f} cm at 5 d",
            fontsize=8.5, color="#a8071a")
    a2.axvline(0, color=MUTED, lw=1.0, ls=":")
    a2.set(xlabel="plateau $-\\ \\Delta p_{crit}$ (MPa)",
           ylabel="Cumulative slip at injector, 5 d (cm)")
    a2.set_title("Slip tracks the excess over $\\Delta p_{crit}$.\n"
                 "Stage 10 moves right by lowering $\\Delta p_{crit}$.")
    a2.grid(alpha=0.3, color=GRID)

    for e in ("png", "pdf"):
        fig.savefig(OUT / f"stage9_mechanism.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/stage9_mechanism.png")
    print(f"\n{'run':>7s} {'kpmax':>9s} {'plateau':>8s} {'peak':>7s} "
          f"{'dp(rw)':>7s} {'excess':>7s} {'R':>6s} {'slip':>8s} "
          f"{'front':>6s} {'wellhd':>8s}")
    for w, e in zip(rows, exc):
        print(f"{w['n']:>7d} {w['kpmax']:>9.1e} {w['plat']:>8.2f} "
              f"{w['wl']-w['plat']:>7.2f} {w['wl']:>7.2f} {e:>+7.2f} "
              f"{w['R']:>5.0f}m {w['pk']:>7.3f}cm {w['fr']:>6.2f} "
              f"{w['wh']:>+7.1f}%")


if __name__ == "__main__":
    main()
