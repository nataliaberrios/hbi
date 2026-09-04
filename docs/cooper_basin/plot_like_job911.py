#!/usr/bin/env python3
"""Cumulative slip in the presentation of the Job 911 figure, for direct overlay.

The reference figure (figures/taiyi_validation/, "constant rate injection
simulation (Job 911)") plots cumulative slip in CM against distance along-strike
over a FIXED +/-1.5 km, at 3,5,...,17 days, with an analytical solution overlaid.
Reproducing those choices exactly is the difference between a comparison and two
pictures that happen to be of the same quantity:

  - cm, not m
  - fixed +/-1.5 km. This matters more than it sounds. The generic slip figure
    auto-zooms to where the slip is, so the same profile renders as a sharp cusp
    on a +/-1.5 km axis and as a rounded dome on a +/-0.6 km one. Nearly all of
    the apparent SHAPE disagreement between the two figures was that zoom, not
    the physics: measured, the profile's r(50% of peak)/r(1% of peak) is 0.471,
    where a triangle gives 0.50.
  - 3 to 17 days in 2-day steps, eight curves, not the generic figure's six
    logarithmically-spaced times.

Usage:  python plot_like_job911.py 632888 [632889 ...]
"""
import argparse
import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import importlib.util as iu

H = Path("/home/users/nberrios/3dhbi/hbi_analysis")
OUT = H / "figures" / "taiyi_validation"
SCRATCH = "/scratch/users/nberrios/3dhbi/output"
# The harness rsyncs output to SCRATCH only when a run FINISHES (with
# --remove-source-files), so a run in flight has its .dat files under
# runs/<jobid>/output/ instead. Look there too, so a partial figure can be made
# while the run is still going -- and label it as partial rather than pretending
# the last curve is the final state.
RUNS = "/scratch/users/nberrios/3dhbi/runs"

TIMES_D = [3, 5, 7, 9, 11, 13, 15, 17]      # as in the reference figure
XLIM_KM = 1.5

_s = iu.spec_from_file_location("sf", str(H / "make_sweep_figures.py"))
sf = iu.module_from_spec(_s)
_s.loader.exec_module(sf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jobs", nargs="+", type=int)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    for n in a.jobs:
        dk = sf.deck(n)
        IM, JM = int(dk["imax"]), int(dk["jmax"])
        ds_m = sf.ffloat(dk["ds"]) * 1000.0
        NC = IM * JM
        import glob as _g
        cands = ([f"{SCRATCH}/{n}/slip{n}.dat"]
                 + sorted(_g.glob(f"{RUNS}/*/output/slip{n}.dat")))
        p = next((q for q in cands if os.path.exists(q) and os.path.getsize(q)), None)
        if p is None:
            print(f"  {n}: no slip output yet")
            continue
        in_flight = "/runs/" in p
        t = np.atleast_2d(np.loadtxt(p.replace("slip", "time")))[:, 1] / 86400.0
        nt = min(os.path.getsize(p) // (8 * NC), len(t))
        arr = np.memmap(p, np.float64, "r", shape=(nt, NC))
        c = (IM - 1) // 2
        # Full along-strike line through the injector, both sides, in km.
        xs = (np.arange(IM) - c) * ds_m / 1000.0

        # Accept a reference time if some snapshot lands within half an output
        # interval of it. A bare `td <= t[-1] + 1e-9` silently DROPPED the 17 d
        # curve here: the run ends at 16.999999100694446 d, which is 9e-7 d
        # (0.08 s) short of 17, six orders outside a 1e-9 tolerance. The figure
        # then looked complete while missing its most important curve.
        dtout_d = sf.ffloat(dk.get("dtout", "0.0002")) * 365.0
        tol = max(0.5 * dtout_d, 1e-6)
        avail = [td for td in TIMES_D
                 if np.min(np.abs(t[:nt] - td)) <= tol]
        if not avail:
            print(f"  {n}: reached only {t[nt-1]:.2f} d, before the first "
                  f"reference time {TIMES_D[0]} d")
            continue
        dropped = [td for td in TIMES_D if td not in avail]
        if dropped:
            print(f"  {n}: reference times not reached: {dropped} "
                  f"(run ends at {t[nt-1]:.4f} d)")

        cmap = plt.get_cmap("viridis")
        cols = [cmap(v) for v in np.linspace(0.0, 0.95, len(avail))]
        fig, ax = plt.subplots(figsize=(11.0, 3.6), constrained_layout=True)
        for td, col in zip(avail, cols):
            k = int(np.argmin(abs(t[:nt] - td)))
            prof = np.asarray(arr[k]).reshape(IM, JM)[c, :]      # along strike
            ax.plot(xs, prof * 100.0, lw=1.9, color=col, label=f"{td} days")
        ax.set(xlabel="Distance along-strike (km)",
               ylabel="Cumulative Slip (cm)", xlim=(-XLIM_KM, XLIM_KM))
        ax.set_ylim(bottom=0)
        pev = dk.get("permev", "F")
        kx, km = sf.ffloat(dk.get("kpmax")), sf.ffloat(dk.get("kpmin"))
        rng = f", enhancement range {kx/km:.0f}x" if pev.upper() == "T" else ""
        partial = "  [PARTIAL -- run still in progress]" if in_flight else ""
        ax.set_title(f"constant rate injection simulation (Job {n}){partial}\n"
                     f"kp {sf.ffloat(dk['kp']):.0e}, ds {ds_m:.0f} m, "
                     f"permev {pev}{rng}", fontsize=11)
        ax.legend(frameon=True, fontsize=9, loc="upper right")
        for e in ("png", "pdf"):
            fig.savefig(OUT / f"job911style_{n}.{e}", dpi=180, bbox_inches="tight")
        plt.close(fig)

        # Numbers for the comparison, at the reference figure's last time.
        k = int(np.argmin(abs(t[:nt] - avail[-1])))
        prof = np.asarray(arr[k]).reshape(IM, JM)[c, c:]
        rr = np.arange(len(prof)) * ds_m
        pk = prof[0]

        def rad(frac):
            b = np.where(prof < frac * pk)[0]
            return rr[b[0]] if len(b) else rr[-1]

        print(f"  {n}: wrote job911style_{n}.png   "
              f"t={t[k]:.2f} d  peak {pk*100:.2f} cm  "
              f"r(1% of peak) {rad(0.01):.0f} m  r50/r1 {rad(0.5)/max(rad(0.01),1):.3f}")


if __name__ == "__main__":
    main()
