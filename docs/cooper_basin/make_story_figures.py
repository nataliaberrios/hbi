#!/usr/bin/env python3
"""The two figures the advisor story needs that did not exist yet.

FIGURE A, story_search_summary.png -- what was searched and what came of it.
Every other figure in this project shows a handful of runs. This one shows all
of them, because the argument is that the search was EXHAUSTIVE, and that claim
needs the whole population visible.

  A1  what was swept: one row per parameter, its range, and how many runs.
  A2  the outcome: front vs wellhead for every scored run, coloured by slip
      ratio, with the +/-15% bands drawn. The target corner is empty.
  A3  THE INVARIANT: slip ratio against wellhead. Every run that reproduces the
      observed slip sits at +70 to +80% wellhead, across different tau_0,
      different maps, and enhancement on and off. This is the single most
      robust result in the project.

FIGURE B, story_activation_pressure.png -- the punchline, from the field data.

  B1  dp_crit: the MEASURED activation pressure (Habanero FDP sec 4 / Holl 2015
      sec 8.4 -- 2.6 MPa for H04 local Oct 2012, ~0.4 MPa for the extended
      Nov 2012 stimulation this project models) against what each model needs.
      The model requires 4-27x the overpressure at which the fault actually
      failed.
  B2  f0 implied by that measurement at the MEASURED tau_0 and sigmabar_0:
      0.375-0.408. f0 = 0.60 is off the scale.
  B3  the stress resolution, so the argument cannot be dismissed as ignoring the
      measurements: Holl & Barton's ratios with Sv ~ 95.3 MPa put tau_0 = 10.36
      and sigmabar_0 = 27.99 squarely inside the admissible band for an
      18-20 deg thrust -- and put Wang & Dunham's 15.0 inside it too.

All scores here are on the 0-5 d window so all 82 runs are comparable; the
figures that use the honest 0-18 d window are the separate Stage 11 ones.

Usage:  python make_story_figures.py
"""
import importlib.util as iu
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = "/home/users/nberrios/3dhbi/hbi_analysis"
_s = iu.spec_from_file_location("sf", H + "/make_sweep_figures.py")
sf = iu.module_from_spec(_s); _s.loader.exec_module(sf)
sys.path.insert(0, H + "/notebooks")
from sim_curves import load_slip, _deck, _ff

OUT = Path(H) / "figures" / "story"
SCORES = Path("/home/users/nberrios/3dhbi/hbi_git/docs/figs/cooper_grid/grid_scores.json")
OBS5 = 2.806
SIG, TAU0 = 27.99, 10.36
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
RED, BLUE, GRN = "#a8071a", "#1d4ed8", "#009E73"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 10.5,
                     "axes.edgecolor": MUTED, "text.color": INK,
                     "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "legend.fontsize": 8.5})

SWEPT = [
    ("initial shear stress $\\tau_0$", "11.0 - 15.0 MPa", 36),
    ("permeability ceiling $k_{pmax}$", "5e-14 - 2.5e-10 (4 decades)", 12),
    ("permeability structure", "uniform / 2-zone, 50x & 250x", 30),
    ("enhancement on/off", "permev T / F", 82),
    ("fluid viscosity $\\eta$ + storage $\\beta$", "8.9e-4 & 1.27e-4, D-matched", 20),
    ("porosity grading $\\phi(r)$", "0.005-0.020, 3 gradings", 6),
    ("friction $a$, $b$", "a 0.015-0.022, b 0.005-0.012", 4),
    ("effective normal stress $\\bar\\sigma_0$", "27.99 & 30.0 MPa", 82),
]


def collect():
    sc = json.load(open(SCORES))
    rows = []
    for r in sc:
        if r.get("lam_ratio") is None or r.get("p_pct") is None:
            continue
        try:
            _, sl, ta = load_slip(r["n"], 5.0, how="strike")
            if abs(ta - 5.0) > 0.3:
                continue
        except Exception:
            continue
        rows.append(dict(n=r["n"], f=r["lam_ratio"], p=r["p_pct"],
                         s=sl[0] * 100 / OBS5, tau0=r["tau0"],
                         pev=str(r.get("permev", "?")).upper().startswith("T")))
    return rows


def fig_a(rows):
    fig = plt.figure(figsize=(15.5, 4.6), dpi=200)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 1], wspace=0.28)
    a1, a2, a3 = (fig.add_subplot(gs[i]) for i in range(3))

    a1.axis("off")
    a1.set_title(f"A. What was searched — {len(rows)} runs with all three\n"
                 "quantities measurable", loc="left")
    y = 0.93
    for name, rng, n in SWEPT:
        a1.text(0.0, y, name, fontsize=9, color=INK, transform=a1.transAxes)
        a1.text(0.60, y, rng, fontsize=8.2, color=MUTED, transform=a1.transAxes)
        y -= 0.115
    a1.text(0.0, y - 0.02, "never varied:  $d_c$, $v_{init}$, rigidity, "
            "$S_{w}$, skin, dilatancy,\n                        viscous flow, "
            "anisotropy", fontsize=8.2, color=RED, transform=a1.transAxes)

    sc = a2.scatter([r["f"] for r in rows], [r["p"] for r in rows],
                    c=[min(r["s"], 2.2) for r in rows], s=42, cmap="viridis",
                    edgecolors=INK, linewidths=0.4, zorder=3)
    a2.axhspan(-15, 15, color=GRID, alpha=0.6, zorder=0)
    a2.axvspan(0.85, 1.15, color=GRID, alpha=0.6, zorder=0)
    a2.set(xlabel="slip front  $\\lambda/\\lambda_{obs}$",
           ylabel="wellhead bias (%)", xlim=(0, 1.3), ylim=(-30, 330))
    a2.set_title("B. Front vs wellhead, every run\n(colour = slip / observed)")
    a2.plot([1], [0], marker="*", ms=22, color=RED, zorder=5)
    a2.annotate("the target", (1, 0), fontsize=9, color=RED, fontweight="bold",
                textcoords="offset points", xytext=(-16, 16), ha="center")
    fig.colorbar(sc, ax=a2, label="slip / observed (capped 2.2)")
    a2.grid(alpha=0.3, color=GRID)

    for r in rows:
        a3.scatter(r["s"], r["p"], s=42, zorder=3,
                   color=(BLUE if r["pev"] else GRN), edgecolors=INK,
                   linewidths=0.4)
    hit = [r for r in rows if abs(r["s"] - 1) <= 0.15]
    a3.axvspan(0.85, 1.15, color=RED, alpha=0.10, zorder=0)
    a3.axhspan(-15, 15, color=GRID, alpha=0.6, zorder=0)
    if hit:
        lo, hi = min(r["p"] for r in hit), max(r["p"] for r in hit)
        # State the enhancement status from the data rather than asserting it:
        # all six slip-matching runs turn out to have permev T, which is a
        # STRONGER statement than "on and off" and was mis-annotated at first.
        non = sum(1 for r in hit if not r["pev"])
        enh = ("all with enhancement ON" if non == 0
               else f"{len(hit)-non} with enhancement on, {non} off")
        a3.annotate(f"all {len(hit)} runs that match the slip sit at\n"
                    f"+{lo:.0f} to +{hi:.0f}% wellhead\n"
                    f"($\\tau_0$ {min(r['tau0'] for r in hit):.1f}–"
                    f"{max(r['tau0'] for r in hit):.1f} MPa, uniform and\n"
                    f"2-zone maps, {enh})",
                    (1.0, hi), fontsize=8.5, color=RED, ha="center",
                    textcoords="offset points", xytext=(28, 96),
                    arrowprops=dict(arrowstyle="->", color=RED, lw=1.2))
    a3.scatter([], [], color=BLUE, label="enhancement on", edgecolors=INK)
    a3.scatter([], [], color=GRN, label="enhancement off", edgecolors=INK)
    a3.set(xlabel="slip / observed at 5 d", ylabel="wellhead bias (%)",
           xlim=(0, 2.3), ylim=(-30, 330))
    a3.set_title("C. THE INVARIANT — matching the slip\ncosts +70 to +80% "
                 "on the wellhead")
    a3.legend(loc="center right"); a3.grid(alpha=0.3, color=GRID)
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"story_search_summary.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/story_search_summary.png   ({len(rows)} runs, "
          f"{len(hit)} matching slip)")
    return hit


def fig_b():
    fig, (b1, b2, b3) = plt.subplots(1, 3, figsize=(15.5, 4.4), dpi=200,
                                     constrained_layout=True)
    names = ["measured\n(H04, 2012)", "Wang & Dunham\n$\\tau_0$ 15.0, $f_0$ 0.60",
             "this model\n$\\tau_0$ 10.36, $f_0$ 0.60", "632913's actual\nplateau"]
    vals = [None, 3.0, 10.72, 11.9]
    cols = [RED, BLUE, "#8E44AD", MUTED]
    b1.bar([0], [2.6 - 0.4], bottom=[0.4], width=0.55, color=RED, alpha=0.35,
           edgecolor=RED, linewidth=1.6)
    b1.text(0, 3.4, "0.4 – 2.6 MPa\nMEASURED", ha="center", fontsize=8.5,
            color=RED, fontweight="bold")
    for i, (v, c) in enumerate(zip(vals[1:], cols[1:]), start=1):
        b1.bar([i], [v], width=0.55, color=c, alpha=0.85)
        b1.text(i, v + 0.4, f"{v:.1f}", ha="center", fontsize=9, color=c)
    b1.axhspan(0.4, 2.6, color=RED, alpha=0.10, zorder=0)
    b1.set_xticks(range(4)); b1.set_xticklabels(names, fontsize=8)
    b1.set(ylabel="$\\Delta p_{crit}$ — overpressure needed to slip (MPa)",
           ylim=(0, 14))
    b1.set_title("A. The fault slipped at ~1 MPa.\nThe model needs ~11.")
    b1.grid(alpha=0.3, color=GRID, axis="y")

    dp = np.linspace(0.05, 12, 400)
    b2.plot(dp, TAU0 / (SIG - dp), lw=2.4, color=INK)
    b2.axvspan(0.4, 2.6, color=RED, alpha=0.18,
               label="measured $\\Delta p_{crit}$, 0.4–2.6 MPa")
    for d, c in ((0.4, RED), (2.6, RED)):
        b2.scatter([d], [TAU0 / (SIG - d)], s=90, color=c, zorder=5,
                   edgecolors=INK)
    b2.axhline(0.60, color=BLUE, ls="--", lw=1.4,
               label="$f_0$ = 0.60 (lab granite, every deck)")
    b2.annotate(f"$f_0$ = {TAU0/(SIG-0.4):.3f} – {TAU0/(SIG-2.6):.3f}",
                (1.5, TAU0 / (SIG - 1.5)), fontsize=10, color=RED,
                fontweight="bold", textcoords="offset points", xytext=(36, 10))
    b2.set(xlabel="$\\Delta p_{crit}$ (MPa)", ylabel="required $f_0$",
           xlim=(0, 12), ylim=(0.3, 0.72))
    b2.set_title("B. $f_0$ implied by the measurement,\nat the MEASURED "
                 "$\\tau_0$ and $\\bar\\sigma_0$")
    b2.legend(loc="upper left"); b2.grid(alpha=0.3, color=GRID)

    PP, SV = 72.7, 95.3
    dips = np.linspace(5, 30, 200)
    for R, c, ls in ((1.35, GRN, "-"), (1.45, GRN, "--")):
        s1, s3 = R * SV, SV
        d = np.radians(dips)
        tau = 0.5 * (s1 - s3) * np.sin(2 * d)
        b3.plot(dips, tau, lw=2, color=c, ls=ls,
                label=f"$\\tau$ on the fault, $S_{{Hmax}}/S_v$ = {R}")
    b3.axhspan(TAU0 - 0.05, TAU0 + 0.05, color="#8E44AD", alpha=0.9)
    b3.annotate(f"this model, $\\tau_0$ = {TAU0}", (6, TAU0), fontsize=9,
                color="#8E44AD", fontweight="bold",
                textcoords="offset points", xytext=(0, 6))
    b3.axhline(15.0, color=BLUE, lw=2.2)
    b3.annotate("Wang & Dunham, $\\tau_0$ = 15.0", (6, 15.0), fontsize=9,
                color=BLUE, fontweight="bold", textcoords="offset points",
                xytext=(0, 6))
    b3.set(xlabel="fault dip (deg)", ylabel="$\\tau$ resolved on the fault (MPa)",
           xlim=(5, 30), ylim=(0, 22))
    b3.set_title("C. Both stress states are admissible.\n"
                 "Holl & Barton cannot distinguish them.")
    b3.legend(loc="upper left"); b3.grid(alpha=0.3, color=GRID)
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"story_activation_pressure.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/story_activation_pressure.png")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    hit = fig_a(collect())
    fig_b()
    print("\nruns matching the observed slip to +/-15%:")
    for r in sorted(hit, key=lambda x: x["p"]):
        print(f"  {r['n']}  tau0 {r['tau0']:>5.2f}  front {r['f']:.2f}  "
              f"wellhead {r['p']:+.1f}%  slip {r['s']:.2f}x  "
              f"permev {'T' if r['pev'] else 'F'}")
