#!/usr/bin/env python3
"""Stage 11's two findings, one figure each. No log axes on slip.

FIGURE 1, stage11_front_window.png -- the 5 d front match does not survive.
Every "matches the front" claim in this project was made on a 0-5 d window.
Left: front radius against sqrt(t) with the observed points, out to 18 d, so
the divergence is visible rather than inferred. Right: lambda_sim/lambda_obs
against the fitting window, with the observed front REFIT on each window so the
comparison is like-for-like. The tuned runs cross 1.0 near 5 d and keep going.

FIGURE 2, stage11_ab.png -- a-b is not the controlling parameter.
632915 and 632917 share a-b = 0.010 with different `a` and differ by 2.3x in
slip. Left: slip at 5 d against a-b, showing two runs at the same a-b with
different answers. Right: the same points on an (a, b) grid, which is where the
pattern actually lives -- amplitude tracks `a`, and b's effect changes sign
with it.

Usage:  python make_stage11_figures.py
"""
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

OUT = Path(H) / "figures" / "stage11"
RUNS = [(632913, "#a8071a"), (632916, "#1d4ed8"), (632917, "#009E73"),
        (632915, "#8E44AD")]
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 10.5,
                     "axes.edgecolor": MUTED, "text.color": INK,
                     "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "legend.fontsize": 8.5})


def lab(n):
    d = _deck(n)
    return f"{n}: a {_ff(d['a']):.3f}, b {_ff(d['b']):.3f}"


def fig_front(obs):
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(11.5, 4.3), dpi=200,
                                 constrained_layout=True)
    a0.scatter(np.sqrt(obs["ft"]), obs["fd"], s=11, color=MUTED,
               label="observed seismicity front", zorder=1)
    wins = np.arange(2.0, 18.5, 0.5)
    for n, c in RUNS:
        dk = sf.deck(n); d = sf.run_data(n, dk)
        if d is None:
            continue
        m = d["T"] <= 18.0
        a0.plot(np.sqrt(d["T"][m]), d["R"][m], lw=2, color=c, label=lab(n))
        rs = []
        for w in wins:
            mm = d["T"] <= w
            mo = obs["ft"] <= w
            if mm.sum() < 3 or mo.sum() < 3:
                rs.append(np.nan); continue
            rs.append(sf.fit(d["T"][mm], d["R"][mm])
                      / sf.fit(obs["ft"][mo], obs["fd"][mo]))
        a1.plot(wins, rs, lw=2, color=c, label=lab(n))
    a0.set(xlabel=r"$\sqrt{t}$  ($\sqrt{\mathrm{days}}$)",
           ylabel="Front radius (km)")
    a0.set_title("Front growth to 18 d — a straight line would be "
                 "$R\\propto\\sqrt{t}$;\nneither is, both step with the "
                 "injection cycles")
    a1.axhline(1.0, color=INK, ls="--", lw=1.2)
    a1.axhspan(0.85, 1.15, color=GRID, alpha=0.55, zorder=0)
    a1.axvline(5.0, color=MUTED, ls=":", lw=1.2)
    a1.text(5.3, 0.50, "the 0–5 d window\nevery earlier claim used",
            fontsize=8, color=MUTED)
    a1.set(xlabel="Fitting window (days, from t = 0)",
           ylabel=r"$\lambda_{sim}/\lambda_{obs}$  (both refit per window)",
           ylim=(0.4, 1.9))
    a1.set_title("The match is a 5 d crossing, not a match")
    a1.text(11.5, 0.50, "shaded: $\\pm$15% band", fontsize=8, color=MUTED)
    for ax in (a0, a1):
        ax.grid(alpha=0.3, color=GRID); ax.legend(loc="upper left")
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"stage11_front_window.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/stage11_front_window.png")


def fig_ab():
    pts = []
    for n in (632901, 632915, 632916, 632917):
        d = _deck(n)
        s = load_slip(n, 5.0, how="strike")[1][0] * 100
        pts.append((n, _ff(d["a"]), _ff(d["b"]), s))
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(11.5, 4.3), dpi=200,
                                 constrained_layout=True)
    for n, a, b, s in pts:
        a0.scatter(a - b, s, s=90, color=INK, zorder=3)
        dx = (9, 6) if n != 632901 else (-14, 12)
        a0.annotate(f"{n}\na {a:.3f}, b {b:.3f}", (a - b, s), fontsize=8,
                    color=MUTED, textcoords="offset points", xytext=dx)
    same = [(p[1] - p[2], p[3]) for p in pts if abs(p[1] - p[2] - 0.010) < 1e-9]
    if len(same) == 2:
        a0.plot([same[0][0]] * 2, [same[0][1], same[1][1]], color="#a8071a",
                lw=2.5, zorder=2)
        a0.annotate("same $a-b$,\n2.3x apart", (same[0][0], np.mean(
            [same[0][1], same[1][1]])), fontsize=9, color="#a8071a",
            ha="right", textcoords="offset points", xytext=(-14, 14))
    a0.axhline(2.806, color="#a8071a", ls="--", lw=1.4)
    a0.text(0.0035, 2.9, "observed 2.81 cm at 5 d", fontsize=8.5,
            color="#a8071a")
    a0.set(xlabel="$a-b$", ylabel="Cumulative slip at injector, 5 d (cm)",
           xlim=(0.0, 0.020), ylim=(0, 3.2))
    a0.set_title("$a-b$ does NOT set the amplitude")
    sc = a1.scatter([p[1] for p in pts], [p[2] for p in pts],
                    c=[p[3] for p in pts], s=430, cmap="viridis",
                    edgecolors=INK, linewidths=0.8, zorder=3)
    for n, a, b, s in pts:
        a1.annotate(f"{s:.2f}", (a, b), fontsize=8.5, ha="center",
                    va="center", color="w", fontweight="bold", zorder=4)
        a1.annotate(str(n), (a, b), fontsize=7.5, ha="center",
                    textcoords="offset points", xytext=(0, 17), color=MUTED)
    a1.set(xlabel="$a$", ylabel="$b$", xlim=(0.012, 0.025),
           ylim=(0.002, 0.015))
    a1.set_title("Amplitude tracks $a$; $b$'s effect flips sign with it\n"
                 "(labels are slip at 5 d in cm)")
    fig.colorbar(sc, ax=a1, label="slip at 5 d (cm)")
    for ax in (a0, a1):
        ax.grid(alpha=0.3, color=GRID)
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"stage11_ab.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/stage11_ab.png")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    fig_front(sf.observed())
    fig_ab()
