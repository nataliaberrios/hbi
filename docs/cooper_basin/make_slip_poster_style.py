#!/usr/bin/env python3
"""Cumulative slip in the poster's visual language: observed | simulated, linear.

The comparison figure this replaces used a log y axis and a normalised-shape
panel. Those answer "where does the front cross a threshold" and "is it the same
shape", but they are not how this figure is read for a talk, and the log axis in
particular makes a factor-of-six amplitude gap look like a small offset. LINEAR
ONLY here -- no log panel.

Conventions taken from the poster (cell 107 of
notebooks/cooper_basin_plots-27_abs_pressure.ipynb, the cell that produced
image16), so these drop into the same slide:

  * observed on the LEFT, simulation on the RIGHT, figsize (16, 3) at dpi 300
  * viridis over np.linspace(0, 1, n_times), one colour per time
  * xlim +/-1.5 km, cumulative slip in CM
  * label fontsize 13, title 14, legend 9
  * BOTH panels forced to a shared y limit, so the amplitude comparison is
    honest -- panels on independent scales hide exactly the discrepancy this
    figure exists to show

Simulated profiles are along STRIKE and mirrored about the injector (the model
is symmetric by construction); the observed profile is genuinely two-sided and
asymmetric, spanning -1.00 to +1.39 km with its peak at -0.136 km. Observed slip
is ALREADY IN CENTIMETRES in slip_profiles_strike.txt; simulated slip is in
metres and gets the factor of 100.

TIME COVERAGE. The observed file has 3,5,...,17 d. A 5 d run therefore
contributes only the 3 and 5 d curves, against the poster's eight, and the
observed panel is trimmed to the same times so the two panels are comparable.
Dropped times are printed and named in the title rather than left implied.
(In the data the 3 d and 5 d columns are byte-identical, so those two observed
curves coincide -- that is the file, not a plotting error.)

Usage:  python make_slip_poster_style.py 632901 [632896 ...]
        python make_slip_poster_style.py --grid --tag story5 632901 632896 ...
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/home/users/nberrios/3dhbi/hbi_analysis/notebooks")
from sim_curves import load_slip, _deck, _ff, _path

OBS = Path("/home/users/nberrios/3dhbi/hbi/slip_profiles_strike.txt")
OBS_TIMES = [3, 5, 7, 9, 11, 13, 15, 17]
XLIM_KM = 1.5
OUT = Path("/home/users/nberrios/3dhbi/hbi_analysis/figures/slip_poster_style")
FS_LABEL, FS_TITLE, FS_LEGEND = 13, 14, 9


def sim_prof(job, td):
    """(x_km, slip_cm) along strike, mirrored, or None if td not reached."""
    r, sl, ta = load_slip(job, td, how="strike")
    d = _deck(job)
    tol = max(0.5 * _ff(d.get("dtout", "0.0002")) * 365.0, 1e-6)
    if abs(ta - td) > max(tol, 0.05):
        return None
    sl_cm = sl * 100.0
    x = np.concatenate([-r[::-1], r[1:]]) / 1000.0
    return x, np.concatenate([sl_cm[::-1], sl_cm[1:]])


def times_for(job):
    d = _deck(job)
    p = _path(job, "slip")
    t = np.atleast_2d(np.loadtxt(p.replace("slip", "time", 1)))[:, 1] / 86400.0
    NC = int(d["imax"]) * int(d["jmax"])
    nt = min(os.path.getsize(p) // (8 * NC), len(t))
    tend = float(t[nt - 1])
    tol = max(0.5 * _ff(d.get("dtout", "0.0002")) * 365.0, 1e-6)
    keep = [td for td in OBS_TIMES if td <= tend + tol]
    return keep, [td for td in OBS_TIMES if td not in keep], tend


def label(job):
    d = _deck(job)
    return (f"Job {job}: $\\tau_0$ = "
            f"{_ff(d['muinit'])*_ff(d['sigmainit']):.2f} MPa, "
            f"kpmax {d.get('kpmax', 'fixed perm')}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jobs", nargs="+", type=int)
    ap.add_argument("--grid", action="store_true",
                    help="one observed panel plus one panel per run")
    ap.add_argument("--tag", default="grid")
    ap.add_argument("--no-titles", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    o = np.loadtxt(OBS)
    x_obs, obs_cm = o[:, 0], o[:, 1:]        # already CENTIMETRES

    def draw_obs(ax, times, colors):
        for c, td in zip(colors, times):
            ax.plot(x_obs, obs_cm[:, OBS_TIMES.index(td)], lw=2, color=c,
                    label=f"{td} days")

    def dress(ax, ymax, title):
        ax.set_xlabel("Distance along-strike (km)", fontsize=FS_LABEL)
        ax.set_ylabel("Cumulative Slip (cm)", fontsize=FS_LABEL)
        ax.set_xlim(-XLIM_KM, XLIM_KM)
        ax.set_ylim(0, ymax)
        ax.tick_params(labelsize=FS_LABEL - 2)
        ax.legend(fontsize=FS_LEGEND, loc="upper right")
        if not a.no_titles:
            ax.set_title(title, fontsize=FS_TITLE)

    if a.grid:
        # shared times = those EVERY run reaches, so all panels are comparable
        per = {j: times_for(j) for j in a.jobs}
        times = [td for td in OBS_TIMES if all(td in per[j][0] for j in a.jobs)]
        if not times:
            sys.exit("no observed time is reached by every run")
        colors = plt.cm.viridis(np.linspace(0, 1, len(times)))
        # Wrap into rows of at most 3. A single row of 6 squashes each panel
        # to ~1/6 of the width and the profiles become unreadable slivers --
        # the poster's pair is (16, 3) for TWO panels, i.e. 8 in wide each,
        # so keep roughly that per panel.
        n = len(a.jobs) + 1
        ncol = min(3, n)
        nrow = -(-n // ncol)
        fig, axg = plt.subplots(nrow, ncol, figsize=(5.3 * ncol, 3.4 * nrow),
                                dpi=300, squeeze=False)
        axes = axg.ravel()
        for ax in axes[n:]:
            ax.axis("off")
        draw_obs(axes[0], times, colors)
        prof = {}
        for ax, j in zip(axes[1:], a.jobs):
            for c, td in zip(colors, times):
                pr = sim_prof(j, td)
                if pr is None:
                    continue
                prof[(j, td)] = pr
                ax.plot(pr[0], pr[1], lw=2, color=c, label=f"{td} days")
        ymax = max([obs_cm[:, [OBS_TIMES.index(t) for t in times]].max()]
                   + [v[1].max() for v in prof.values()]) * 1.08
        dress(axes[0], ymax, "Observed, 20 m west of well")
        for ax, j in zip(axes[1:], a.jobs):
            dress(ax, ymax, label(j))
        plt.tight_layout()
        for e in ("png", "pdf"):
            fig.savefig(OUT / f"slip_poster_{a.tag}"
                        f"{'_notitle' if a.no_titles else ''}.{e}",
                        bbox_inches="tight")
        plt.close(fig)
        print(f"wrote slip_poster_{a.tag}.png   times {times}, shared y 0-"
              f"{ymax:.2f} cm")
        drop = sorted({td for j in a.jobs for td in per[j][1]})
        if drop:
            print(f"  observed times not reached by every run, DROPPED: {drop}")
        print(f"\n  {'run':>8s} " + " ".join(f"{t:>7d}d" for t in times))
        print(f"  {'observed':>8s} " + " ".join(
            f"{obs_cm[:, OBS_TIMES.index(t)].max():>8.3f}" for t in times))
        for j in a.jobs:
            print(f"  {j:>8d} " + " ".join(
                f"{prof[(j,t)][1].max():>8.3f}" if (j, t) in prof else
                f"{'--':>8s}" for t in times))
        return

    for j in a.jobs:
        times, drop, tend = times_for(j)
        if not times:
            print(f"  {j}: ends at {tend:.2f} d, before {OBS_TIMES[0]} d")
            continue
        colors = plt.cm.viridis(np.linspace(0, 1, len(times)))
        fig, (ax_o, ax_s) = plt.subplots(1, 2, figsize=(16, 3), dpi=300)
        draw_obs(ax_o, times, colors)
        pk = []
        for c, td in zip(colors, times):
            pr = sim_prof(j, td)
            if pr is None:
                continue
            ax_s.plot(pr[0], pr[1], lw=2, color=c, label=f"{td} days")
            pk.append(pr[1].max())
        ymax = max([obs_cm[:, [OBS_TIMES.index(t) for t in times]].max()]
                   + pk) * 1.08
        dress(ax_o, ymax, "Observed, 20 m west of well")
        dress(ax_s, ymax, label(j))
        plt.tight_layout()
        for e in ("png", "pdf"):
            fig.savefig(OUT / f"slip_poster_{j}"
                        f"{'_notitle' if a.no_titles else ''}.{e}",
                        bbox_inches="tight")
        plt.close(fig)
        print(f"  {j}: wrote slip_poster_{j}.png   times {times}, "
              f"shared y 0-{ymax:.2f} cm, run ends {tend:.2f} d")
        if drop:
            print(f"      observed times NOT reached, dropped: {drop}")


if __name__ == "__main__":
    main()
