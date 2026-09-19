#!/usr/bin/env python3
"""PAPER FIGURE: the joint constraint is a crossing of two families, not one
lucky run.

The complete (tau_0 x disc radius) grid -- 5 stress states by 7 initial
enhanced-zone radii, all 35 cells run with everything else held fixed. Each
observable is shown as model/observed over that plane, so a perfect match is
the contour at 1.

    (a)  seismicity front       lambda/lambda_obs
    (b)  wellhead pressure      dp/dp_obs
    (c)  total slip             delta/delta_obs
    (d)  the three unity contours on one axis

PANEL (d) IS THE FIGURE. The lambda = 1 and dp = 1 contours are not parallel:
the front contour runs nearly horizontally (the front is set by tau_0 and
barely notices the disc) while the pressure contour runs diagonally (pressure
responds to both). Two non-parallel curves cross at a point, and that point is
the joint constraint. It is a crossing of two one-parameter families, which is
a far stronger statement than "this cell fits" -- a single best-fitting run
could be coincidence; an intersection cannot.

AND IT SHOWS THE OPEN QUESTION GEOMETRICALLY. The slip = 1 contour does not
pass through that intersection. Whether it is absent from the plotted range
entirely or merely offset, the reader sees immediately that two observables
agree about the fault and the third does not -- which is exactly the state of
section 8 and is more honest than a sentence in the discussion.

WHY THE DISC MATTERS LESS AS STRESS RISES, visible in (b): at tau_0 10.36 the
disc spans ~29 points of dp; at tau_0 15.00 only ~3. More slip means more
permeability enhancement, which swamps the initial disc, so the disc is a lever
only at low prestress. That is why the pressure contour is diagonal rather than
vertical.

RATIOS ARE PLOTTED AS log2, so 2x and 0.5x are equally far from the unity
contour and a diverging colormap centred on zero is honest about direction.
Colorbar ticks are relabelled back to ratios.

DEFINITIONS ARE IMPORTED, NOT REIMPLEMENTED -- seismicity_front and lam from
compare_arms, the datum from fault_pressure, grid membership from
compare_disc.discover(). This figure cannot disagree with Fig. 4 or the scorer.

Usage:  python fig_grid.py
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
_a = iu.spec_from_file_location(
    "ca", "/home/users/nberrios/3dhbi/hbi_git/docs/cooper_basin/compare_arms.py")
ca = iu.module_from_spec(_a); _a.loader.exec_module(ca)
_f = iu.spec_from_file_location(
    "fp", "/home/users/nberrios/3dhbi/hbi_git/docs/cooper_basin/fault_pressure.py")
fp = iu.module_from_spec(_f); _f.loader.exec_module(fp)
_d = iu.spec_from_file_location(
    "cd", "/home/users/nberrios/3dhbi/hbi_git/docs/cooper_basin/compare_disc.py")
cd = iu.module_from_spec(_d); _d.loader.exec_module(cd)
sys.path.insert(0, H + "/notebooks")
from sim_curves import load_slip

FIG = Path(H) / "figures" / "cycle2"
OBS_SLIP = Path("/home/users/nberrios/3dhbi/hbi/slip_profiles_strike.txt")
OT = [3, 5, 7, 9, 11, 13, 15, 17]
T0, T_EVAL = 4.300, 8.7
INK = "#141413"
C_FRONT, C_PRESS, C_SLIP = "#1d4ed8", "#0e7a63", "#D55E00"

plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 9, "legend.fontsize": 8,
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
    ax.text(-0.17, 1.02, f"({letter})", transform=ax.transAxes,
            fontsize=10, fontweight="bold", va="bottom", ha="left")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="fig_grid")
    a = ap.parse_args(argv)

    grid = cd.discover()
    taus = sorted(grid)
    discs = sorted({d for v in grid.values() for d, _ in v})

    obs = sf.observed()
    tcat, rcat, _ = ca.fb.catalogue()
    k = (tcat - T0) > 0
    te, re_ = tcat[k] - T0, rcat[k]
    org = float(np.median(re_[np.argsort(te)][:10]))
    ft, fd = ca.seismicity_front(te, re_ - org)
    LOBS = ca.lam(ft, fd)
    to, dpo = fp.dp_observed(obs)
    o = np.loadtxt(OBS_SLIP)
    peak = [o[:, 1 + i].max() for i in range(8)]
    s0 = float(np.interp(T0, OT, peak))
    S_OBS = float(np.interp(T_EVAL + T0, OT, peak)) - s0
    P_STATIC = sf.P0 - sf.RHO * sf.G * sf.HW / 1e6

    shape = (len(taus), len(discs))
    L = np.full(shape, np.nan)
    P = np.full(shape, np.nan)
    S = np.full(shape, np.nan)

    print(f"observed: lambda {LOBS:.1f}, slip {S_OBS:.3f} cm, "
          f"datum {fp.datum():.3f} MPa")
    print(f"{'tau_0':>6} {'disc':>5} {'run':>7} {'lam':>6} {'dp':>6} {'slip':>6}")
    for i, tau in enumerate(taus):
        cells = dict(grid[tau])
        for j, disc in enumerate(discs):
            n = cells.get(disc)
            if n is None:
                continue
            d = sf.run_data(n, sf.deck(n))
            if d is None or len(d["T"]) < 3 or d.get("tpw") is None:
                continue
            oo = np.argsort(d["T"])
            L[i, j] = ca.lam(np.asarray(d["T"])[oo],
                             np.asarray(d["R"])[oo] * 1000.0) / LOBS

            dpm = d["ppw"] - P_STATIC
            gr = np.linspace(0.05, min(d["tpw"][-1], to.max(), T_EVAL), 2000)
            ps = np.interp(gr, d["tpw"], dpm)
            ob = np.interp(gr, to, dpo)
            qg = np.interp(gr, obs["ti"] - T0, obs["q"])
            fl = (qg > 0.25 * np.nanmax(obs["q"])) \
                & (np.interp(gr, obs["tp"] - T0, obs["pm"]) > 5.0)
            P[i, j] = float(np.mean(ps[fl])) / float(np.mean(ob[fl]))

            try:                       # one sample, at the scoring time only
                _, sl, ta = load_slip(n, T_EVAL, how="strike")
                if abs(ta - T_EVAL) < 0.4:
                    S[i, j] = sl[0] * 100 / S_OBS
            except Exception:
                pass
            print(f"{tau:>6.2f} {disc:>5} {n:>7} {L[i,j]:>6.2f} "
                  f"{P[i,j]:>6.2f} {S[i,j]:>6.2f}")

    X, Y = np.meshgrid(np.array(discs, float), np.array(taus, float))
    fig, ax = plt.subplots(2, 2, figsize=(7.48, 5.9), dpi=400,
                           constrained_layout=True)
    (af, apz), (asl, aj) = ax

    # ONE SHARED COLOUR SCALE ACROSS (a)-(c). Letting each panel autoscale made
    # the same colour mean 2.0x in one panel and 1.3x in the next, which is
    # actively misleading when the three are read side by side.
    V = float(np.nanmax([np.nanmax(np.abs(np.log2(Z))) for Z in (L, P, S)]))
    for axx, Z in ((af, L), (apz, P), (asl, S)):
        pc = axx.pcolormesh(X, Y, np.log2(Z), cmap="RdBu_r", shading="nearest",
                            vmin=-V, vmax=V)
        if np.nanmin(Z) < 1 < np.nanmax(Z):
            axx.contour(X, Y, np.log2(Z), levels=[0.0], colors=[INK],
                        linewidths=1.8)
        axx.scatter(X, Y, s=5, color=INK, alpha=0.35, lw=0, zorder=4)
    cb = fig.colorbar(pc, ax=[af, apz, asl], location="right", shrink=0.55,
                      aspect=26, pad=0.015)
    # TICK AT ROUND RATIOS, NOT ROUND log2 VALUES. Ticking at +-0.5, +-1, +-1.5
    # in log2 and relabelling gave "2.82843" and "0.707107" on the bar.
    rt = [r for r in (0.33, 0.5, 1, 2, 3) if abs(np.log2(r)) <= V]
    cb.set_ticks([np.log2(r) for r in rt])
    cb.set_ticklabels([f"{r:g}" for r in rt])
    cb.ax.tick_params(length=2.5, width=0.6)
    cb.outline.set_linewidth(0.6)
    cb.set_label("model / observed", fontsize=8.5)

    # --- panel (d): the unity contours together
    hs, ls = [], []
    for Z, col, lab in ((L, C_FRONT, "seismicity front"),
                        (P, C_PRESS, "wellhead pressure"),
                        (S, C_SLIP, "total slip")):
        if np.nanmin(Z) < 1 < np.nanmax(Z):
            aj.contour(X, Y, np.log2(Z), levels=[0.0], colors=[col],
                       linewidths=2.0)
            hs.append(plt.Line2D([], [], color=col, lw=2.0)); ls.append(lab)
        else:
            lo, hi = np.nanmin(Z), np.nanmax(Z)
            hs.append(plt.Line2D([], [], color=col, lw=2.0, ls=":"))
            ls.append(f"{lab} ({lo:.2f}–{hi:.2f}$\\times$, no match)")
    aj.scatter(X, Y, s=6, color=INK, alpha=0.35, lw=0)
    aj.legend(hs, ls, loc="center", handlelength=1.8, labelspacing=0.35)
    # the band around tau_0 13-14 is the only empty region: the slip
    # contour sits at 14.9 and the other two cross near 11.5

    for axx, l, t in ((af, "a", "Seismicity front"),
                      (apz, "b", "Wellhead pressure"),
                      (asl, "c", "Total slip"),
                      (aj, "d", "Unity contours")):
        axx.set(xlabel="Initial enhanced-zone radius (m)",
                ylabel=r"$\tau_0$ (MPa)",
                xlim=(min(discs) - 25, max(discs) + 25),
                ylim=(min(taus) - 0.3, max(taus) + 0.3))
        panel(axx, l)

    for e in ("png", "pdf"):
        fig.savefig(FIG / f"{a.out}.{e}")
    plt.close(fig)
    for nm, Z in (("front", L), ("pressure", P), ("slip", S)):
        print(f"{nm:>9}: {np.nanmin(Z):.2f} to {np.nanmax(Z):.2f}"
              + ("" if np.nanmin(Z) < 1 < np.nanmax(Z)
                 else "   <-- never reaches 1 on this grid"))
    print(f"wrote {FIG}/{a.out}.png")


if __name__ == "__main__":
    main()
