#!/usr/bin/env python3
"""Slip distributions: the field the slip front is derived from, plotted directly.

The R-T and R-V figures reduce the whole slip field to one number per time (where
slip first crosses dc). That hides shape: two runs can share a front position while
one has a broad low-amplitude bowl and the other a narrow deep patch. These figures
show the field itself.

Per run, into figures/<job>/:
    slipprof_<job>.png   slip along dip and along strike, at a series of times
    slipmap_<job>.png    2D cumulative slip at the final time, with the dc contour
                         (= the front the R-T figures track) drawn on it

Per sweep, into figures/sweeps/:
    slip_<sweepkey>.png  every run's dip and strike profile at one shared time,
                         so the sweep's effect on slip SHAPE is visible

Conventions follow notebooks/slip_profiles_dip_max.txt: distance along dip in km,
slip in m, profiles taken through the injector row/column.

The dc contour matters: the front is defined as slip > dc, so a run whose slip
never reaches dc anywhere has no front at all, and that shows up here as a map with
no contour rather than as a missing point in an R-T plot.

Usage:  python make_slip_figures.py 632721 632750 ...   # per-run
        python make_slip_figures.py --sweeps            # sweep overlays
        python make_slip_figures.py --all               # both, all sweep members
"""
import argparse
import glob
import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = Path("/home/users/nberrios/3dhbi/hbi_analysis")
FIG = H / "figures"
IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")

# validated with scripts/validate_palette.js --mode light: all six checks PASS
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
DCC = "#a8071a"

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "axes.edgecolor": MUTED, "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "legend.fontsize": 8.5,
})

# candidate snapshot times in days; filtered to those inside each run
TIMES = [0.5, 1, 2, 4, 6, 8, 12, 17, 24, 30]


def style(ax):
    ax.grid(True, color=GRID, lw=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return ax


def ffloat(x):
    try:
        return float(str(x).replace("d", "e").replace("D", "e"))
    except (TypeError, ValueError):
        return None


def deck(n):
    d = {}
    for ln in open(IN / f"res{n}.in"):
        w = ln.split()
        if len(w) >= 2 and not ln.startswith("!"):
            d[w[0]] = w[1].strip('"')
    return d


def load(n, dk):
    g = (glob.glob(f"/scratch/users/nberrios/3dhbi/runs/*/output/slip{n}.dat")
         + glob.glob(f"/scratch/users/nberrios/3dhbi/output/{n}/slip{n}.dat")
         + glob.glob(f"/oak/stanford/groups/edunham/nberrios/3doutput/{n}/slip{n}.dat"))
    if not g:
        return None
    p = g[0]
    IM, JM = int(dk["imax"]), int(dk["jmax"])
    ds, NC = ffloat(dk["ds"]), IM * JM
    t = np.atleast_2d(np.loadtxt(p.replace("slip", "time")))[:, 1] / 86400.0
    nt = min(os.path.getsize(p) // (8 * NC), len(t))
    if nt < 3:
        return None
    return dict(mm=np.memmap(p, np.float64, "r", shape=(nt, NC)), t=t[:nt],
                IM=IM, JM=JM, ds=ds, dc=ffloat(dk["dc"]))


def frame(d, k):
    return np.asarray(d["mm"][k]).reshape(d["IM"], d["JM"])


def at_time(d, day):
    """Nearest recorded frame to a wall-clock day, with its actual time."""
    k = int(np.argmin(np.abs(d["t"] - day)))
    return k, d["t"][k]


def axes_km(d):
    x = (np.arange(d["IM"]) - d["IM"] // 2) * d["ds"]      # dip, km
    y = (np.arange(d["JM"]) - d["JM"] // 2) * d["ds"]      # strike, km
    return x, y


def label(n, dk):
    pb = (ffloat(dk.get("phi")) or 0) * (ffloat(dk.get("beta")) or 0)
    bits = [f"$\\mu_0$={ffloat(dk['muinit']):.2f}",
            f"$\\bar\\sigma_0$={ffloat(dk['sigmainit']):.2f} MPa",
            f"$\\phi\\beta$={pb:.1e}", f"dc={ffloat(dk['dc']):.2e}"]
    if dk.get("permev", "F").upper().startswith("T"):
        bits.append(f"permev T, kpmax={dk.get('kpmax','unset')}, kL={dk.get('kL','?')}")
    else:
        bits.append("permev F")
    return f"{n}:  " + ",  ".join(bits)


def per_run(n):
    dk = deck(n)
    d = load(n, dk)
    if d is None:
        print(f"  {n}: no usable slip output")
        return
    folder = FIG / str(n)
    folder.mkdir(parents=True, exist_ok=True)
    x, y = axes_km(d)
    times = [td for td in TIMES if td <= d["t"][-1] + 1e-9][-6:]
    if not times:
        times = [d["t"][-1]]
    cmap = plt.get_cmap("viridis")
    cols = [cmap(v) for v in np.linspace(0.06, 0.88, len(times))]

    # ---- profiles along dip and strike
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), constrained_layout=True)
    for ax, which in ((axes[0], "dip"), (axes[1], "strike")):
        for td, c in zip(times, cols):
            k, ta = at_time(d, td)
            f = frame(d, k)
            prof = f[:, d["JM"] // 2] if which == "dip" else f[d["IM"] // 2, :]
            ax.plot(x if which == "dip" else y, prof, lw=1.8, color=c,
                    label=f"{ta:.2f} d")
        ax.axhline(d["dc"], color=DCC, lw=1.3, ls="--")
        ax.annotate(f"dc = {d['dc']:.1e} m  (front threshold)",
                    xy=(0.99, d["dc"]), xycoords=("axes fraction", "data"),
                    xytext=(0, 4), textcoords="offset points", ha="right",
                    fontsize=8.5, color=DCC)
        ax.set_yscale("log")
        style(ax).set(xlabel=f"distance along {which} from injector (km)",
                      ylabel="cumulative slip (m)",
                      title=f"Slip profile along {which} (through the injector)")
        ax.legend(frameon=False, labelcolor=INK, ncol=2)
    fig.suptitle(label(n, dk) + "\nlog scale, so the dc crossing that defines the "
                 "front is readable", fontsize=11)
    for e in ("png", "pdf"):
        fig.savefig(folder / f"slipprof_{n}.{e}", dpi=180, bbox_inches="tight")
    plt.close(fig)

    # ---- 2D map at the final frame
    k = len(d["t"]) - 1
    f = frame(d, k)
    fig, ax = plt.subplots(figsize=(7.4, 6.1), constrained_layout=True)
    mx = float(f.max())
    if mx <= 0:
        ax.text(0.5, 0.5, "no slip anywhere", transform=ax.transAxes,
                ha="center", fontsize=14, color=DCC)
        ax.set(xticks=[], yticks=[])
    else:
        im = ax.pcolormesh(y, x, f, shading="auto", cmap="magma",
                           vmin=0, vmax=mx)
        cb = fig.colorbar(im, ax=ax, pad=0.02)
        cb.set_label("cumulative slip (m)")
        cb.outline.set_edgecolor(MUTED)
        if mx > d["dc"]:
            ax.contour(y, x, f, levels=[d["dc"]], colors=[DCC], linewidths=1.8)
        ax.set_aspect("equal")
    ax.set(xlabel="along strike (km)", ylabel="along dip (km)")
    ax.set_title(f"{label(n, dk)}\ncumulative slip at {d['t'][k]:.2f} d; "
                 f"red contour = dc = {d['dc']:.1e} m (the front)", fontsize=9.5)
    for e in ("png", "pdf"):
        fig.savefig(folder / f"slipmap_{n}.{e}", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  {n}: wrote slipprof_{n}, slipmap_{n}  (max slip {mx:.3e} m "
          f"at {d['t'][k]:.2f} d)")


def sweep_overlay(sw):
    """Every run in a sweep, profiles at one shared time."""
    runs = []
    for n in sw["runs"]:
        if not (IN / f"res{n}.in").exists():
            continue
        dk = deck(n)
        d = load(n, dk)
        if d is not None:
            runs.append((n, dk, d))
    if len(runs) < 2:
        print(f"  {sw['key']}: fewer than 2 runs, skipping")
        return
    tcut = min(d["t"][-1] for _, _, d in runs)
    cmap = plt.get_cmap("viridis")
    cols = [cmap(v) for v in np.linspace(0.06, 0.86, len(runs))]
    if sw.get("pairs"):
        cmap = plt.get_cmap("tab10")
        cols = [cmap((i // 2) % 10) for i in range(len(runs))]

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.4), constrained_layout=True)
    for ax, which in ((axes[0], "dip"), (axes[1], "strike")):
        for (n, dk, d), c in zip(runs, cols):
            x, y = axes_km(d)
            k, ta = at_time(d, tcut)
            f = frame(d, k)
            prof = f[:, d["JM"] // 2] if which == "dip" else f[d["IM"] // 2, :]
            ls = "--" if (sw.get("pairs")
                          and dk.get("permev", "F").upper().startswith("T")) else "-"
            ax.plot(x if which == "dip" else y, prof, lw=1.9, ls=ls, color=c,
                    label=f"{n}  {sw['param'](dk)}")
            if which == "dip":
                ax.axhline(d["dc"], color=DCC, lw=1.0, ls=":", zorder=0)
        ax.axhline(runs[0][2]["dc"], color=DCC, lw=1.3, ls="--")
        ax.annotate("dc (front threshold)", xy=(0.99, runs[0][2]["dc"]),
                    xycoords=("axes fraction", "data"), xytext=(0, 4),
                    textcoords="offset points", ha="right", fontsize=8.5, color=DCC)
        ax.set_yscale("log")
        style(ax).set(xlabel=f"distance along {which} from injector (km)",
                      ylabel="cumulative slip (m)",
                      title=f"Along {which}")
        ax.legend(frameon=False, labelcolor=INK)
    fig.suptitle(f"Slip distribution — {sw['title']}\n{sw['fixed']}\n"
                 f"all profiles at the shared time {tcut:.2f} d", fontsize=11.5)
    out = FIG / "sweeps"
    out.mkdir(parents=True, exist_ok=True)
    for e in ("png", "pdf"):
        fig.savefig(out / f"slip_{sw['key']}.{e}", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  {sw['key']}: wrote sweeps/slip_{sw['key']}.png  (at {tcut:.2f} d)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jobs", nargs="*", type=int)
    ap.add_argument("--sweeps", action="store_true")
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()

    import importlib.util as iu
    spec = iu.spec_from_file_location("sf", str(H / "make_sweep_figures.py"))
    sf = iu.module_from_spec(spec)
    spec.loader.exec_module(sf)

    jobs = list(a.jobs)
    if a.all:
        jobs = sorted({n for sw in sf.SWEEPS for n in sw["runs"]})
    for n in jobs:
        per_run(n)
    if a.sweeps or a.all:
        for sw in sf.SWEEPS:
            sweep_overlay(sw)


if __name__ == "__main__":
    main()
