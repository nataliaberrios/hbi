#!/usr/bin/env python3
"""Simulated slip against the OBSERVED slip profiles, for any run.

The per-run suite (make_run_figures + make_slip_figures) shows a run against
itself: pressure history, front R(t), R(V), slip profile, slip map. None of it
puts the run next to the data. This does, and it is the figure that shows the
current state of the match honestly, because the three things being matched
separate cleanly into three panels:

  A  slip vs distance, linear. The AMPLITUDE comparison. 632901 sits 5.7x below
     the observed curves here, which is the open problem.
  B  the same on a log y axis with the 1e-4 m front threshold drawn. The FRONT
     comparison. The curves cross the threshold at nearly the same radius even
     though they are a factor of 6 apart in height -- that is the cliff at the
     disc edge, and it is why the front matches while the amplitude does not.
  C  each profile divided by its own peak. The SHAPE comparison, which asks
     whether the simulated crack is the same object scaled down or a different
     one. Answer for 632901: same shape to about 10% over the inner half.

Observed data is slip_profiles_strike.txt -- cumulative slip 20 m west of the
well, inferred from the seismicity catalog under a circular crack model, at
3,5,...,17 d.

TWO PROPERTIES OF THAT FILE THAT BIT AN EARLIER VERSION OF THIS SCRIPT.
  * Its slip columns are ALREADY IN CENTIMETRES, not metres. Multiplying by 100
    reported the 5 d peak as 51 cm instead of 2.81. Simulated slip IS in metres
    and still needs the factor -- the two sources are in different units.
  * The observed profile is two-sided and ASYMMETRIC, spanning -1.00 to
    +1.39 km and peaking at x = -0.136 km, not at the injector. Folding to
    x >= 0 therefore discards the peak. Both sides are plotted, as Fig 6 does,
    and the simulated profile (symmetric by construction) is mirrored.
  * The 3 d and 5 d columns are byte-identical in the file, so those two
    observed curves coincide. That is the data, not a plotting error. Only the times the run actually reaches are drawn; a run that
stops at 5 d gets 3 and 5 d and the rest are reported as dropped rather than
extrapolated.

DIRECTION. Simulated profiles are along STRIKE, frame[:, c]. In coordinate3ddip
(main_LH.f90:1577) j is innermost and controls y,z, so reshape(IM,JM)[a,b] has
a = strike and b = dip. The observed file is a real geographic along-strike
profile, so strike is the consistent choice. There is a separate
slip_profiles_dip.txt on a DIFFERENT time set (3,6,...,21 d); it is not
interchangeable, and pairing simulated dip against observed strike was a bug
fixed earlier in this project.

COMPARE MODE (--compare --tag NAME --at DAYS) puts several runs on shared axes
against the observed profile at ONE time, which is the slip-distribution
counterpart of pressure_compare_*/RT_compare_*/RV_compare_*. Per-run mode varies
time for one run; compare mode varies run at one time.

Usage:  python make_obs_overlay.py 632901 [632900 ...]
        python make_obs_overlay.py --compare --tag story5 --at 5 632901 632896 ...
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
FRONT_THR = 1e-4          # the same fixed threshold score_grid.py uses
OUTROOT = Path("/home/users/nberrios/3dhbi/hbi_analysis/figures")
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 10.5,
                     "axes.edgecolor": MUTED, "text.color": INK,
                     "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "legend.fontsize": 8.5})


def front(r, s, thr=FRONT_THR):
    b = np.where(s < thr)[0]
    return r[b[0]] if len(b) else r[-1]


def obs_at(x_km, obs_cm, td):
    """(x_m, slip_cm, front_m) for the observed profile at time td."""
    i = OBS_TIMES.index(td)
    xo_m, so = x_km * 1000.0, obs_cm[:, i]
    pos = xo_m >= 0
    Ro = max(front(xo_m[pos], so[pos] / 100.0),
             front(-xo_m[~pos][::-1], so[~pos][::-1] / 100.0))
    return xo_m, so, Ro


def compare(jobs, tag, td, x_km, obs_cm):
    """Several runs against the observed profile at ONE time, shared axes."""
    xo_m, so_cm, Ro = obs_at(x_km, obs_cm, td)
    cols = plt.cm.viridis(np.linspace(0.05, 0.85, len(jobs)))
    fig, (a0, a1, a2) = plt.subplots(1, 3, figsize=(15.0, 4.3), dpi=200,
                                     constrained_layout=True)
    for ax in (a0, a1):
        ax.plot(xo_m, so_cm, lw=2.6, color="#a8071a", zorder=5,
                label=f"OBSERVED {td} d")
    a2.plot(xo_m, so_cm / so_cm.max(), lw=2.6, color="#a8071a", zorder=5,
            label=f"OBSERVED {td} d")
    rows = []
    for c, job in zip(cols, jobs):
        d = _deck(job)
        try:
            r, sl, ta = load_slip(job, td, how="strike")
        except FileNotFoundError:
            print(f"  [skip] {job}: no slip output"); continue
        if abs(ta - td) > 0.05:
            print(f"  [skip] {job}: nearest frame to {td} d is {ta:.2f} d")
            continue
        sl_cm = sl * 100.0
        xs = np.concatenate([-r[::-1], r[1:]])
        ss = np.concatenate([sl_cm[::-1], sl_cm[1:]])
        lab = (f"{job}: tau_0 {_ff(d['muinit'])*_ff(d['sigmainit']):.2f} MPa, "
               f"kpmax {d.get('kpmax','-')}")
        for ax in (a0, a1):
            ax.plot(xs, ss, lw=1.9, color=c, label=lab)
        a2.plot(xs, ss / max(ss.max(), 1e-30), lw=1.9, color=c, label=lab)
        rows.append((job, sl_cm[0], front(r, sl)))
    a0.set(xlabel="Distance along strike from injector (m)",
           ylabel="Cumulative slip (cm)", xlim=(-700, 700))
    a0.set_ylim(bottom=0)
    a0.set_title(f"A. Amplitude at {td} d")
    a1.set(xlabel="Distance along strike from injector (m)",
           ylabel="Cumulative slip (cm)", xlim=(-700, 700), yscale="log",
           ylim=(1e-4, 20))
    a1.axhline(FRONT_THR * 100, color=INK, ls=":", lw=1.5)
    a1.text(-660, FRONT_THR * 100 * 1.35,
            f"front threshold {FRONT_THR:.0e} m", fontsize=8, color=INK)
    a1.set_title("B. Front — where each profile crosses the threshold")
    a2.set(xlabel="Distance along strike from injector (m)",
           ylabel="Slip / peak slip", xlim=(-700, 700), ylim=(0, 1.05))
    a2.set_title("C. Shape — each profile over its own peak")
    for ax in (a0, a1, a2):
        ax.grid(alpha=0.3, color=GRID)
        ax.legend(loc="upper right", fontsize=7.5)
    fig.suptitle(f"Slip distribution against the seismicity-derived observation "
                 f"at {td} d.  The observed lobe is ONE-SIDED "
                 f"(-400 to +80 m, peak at -150 m); the simulations are "
                 f"symmetric about the injector.", fontsize=10)
    OUTROOT.mkdir(parents=True, exist_ok=True)
    for e in ("png", "pdf"):
        fig.savefig(OUTROOT / f"cooper_basin_calibration" /
                    f"slip_compare_{tag}.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote slip_compare_{tag}.png/.pdf   (observed peak "
          f"{so_cm.max():.3f} cm, observed R {Ro:.0f} m)")
    print(f"  {'run':>8s} {'peak slip':>10s} {'/obs':>6s} {'R':>6s} {'/obs':>6s}")
    for job, pk, R in rows:
        print(f"  {job:>8d} {pk:>9.3f}cm {pk/so_cm.max():>6.2f} "
              f"{R:>5.0f}m {R/Ro:>6.2f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jobs", nargs="+", type=int)
    ap.add_argument("--compare", action="store_true",
                    help="several runs on shared axes at one time")
    ap.add_argument("--tag", default="compare")
    ap.add_argument("--at", type=int, default=5, choices=OBS_TIMES,
                    help="observed time to compare at (default 5 d)")
    a = ap.parse_args()

    o = np.loadtxt(OBS)
    # column 0 is km; the slip columns are already CENTIMETRES
    x_obs_km, obs_cm = o[:, 0], o[:, 1:]

    if a.compare:
        (OUTROOT / "cooper_basin_calibration").mkdir(parents=True, exist_ok=True)
        compare(a.jobs, a.tag, a.at, x_obs_km, obs_cm)
        return

    for job in a.jobs:
        d = _deck(job)
        p = _path(job, "slip")
        t = np.atleast_2d(np.loadtxt(p.replace("slip", "time", 1)))[:, 1] / 86400.0
        NC = int(d["imax"]) * int(d["jmax"])
        nt = min(os.path.getsize(p) // (8 * NC), len(t))
        tend = float(t[nt - 1])
        dtout_d = _ff(d.get("dtout", "0.0002")) * 365.0
        tol = max(0.5 * dtout_d, 1e-6)
        avail = [i for i, td in enumerate(OBS_TIMES) if td <= tend + tol]
        dropped = [OBS_TIMES[i] for i in range(len(OBS_TIMES)) if i not in avail]
        if not avail:
            print(f"  {job}: ends at {tend:.2f} d, before the first observed "
                  f"time {OBS_TIMES[0]} d -- nothing to compare")
            continue
        times = [OBS_TIMES[i] for i in avail]
        cols = plt.cm.viridis(np.linspace(0.05, 0.8, len(times)))

        fig, (a0, a1, a2) = plt.subplots(1, 3, figsize=(15.0, 4.3), dpi=200,
                                         constrained_layout=True)
        rows = []
        for k, (i, td) in enumerate(zip(avail, times)):
            r, sl, ta = load_slip(job, td, how="strike")
            xo_m, so_cm = x_obs_km * 1000.0, obs_cm[:, i]
            sl_cm = sl * 100.0                      # sim is metres
            # simulated profile mirrored about the injector for a like-for-like
            # two-sided picture; the observed one is genuinely two-sided.
            xs = np.concatenate([-r[::-1], r[1:]])
            ss = np.concatenate([sl_cm[::-1], sl_cm[1:]])
            # front = max extent on EITHER side, in metres
            Rs = front(r, sl)
            pos = xo_m >= 0
            Ro = max(front(xo_m[pos], so_cm[pos] / 100.0),
                     front(-xo_m[~pos][::-1], so_cm[~pos][::-1] / 100.0))
            rows.append((td, ta, sl_cm[0], so_cm.max(), Rs, Ro))
            for ax in (a0, a1):
                ax.plot(xo_m, so_cm, lw=2.2, color=cols[k],
                        label=f"observed {td} d")
                ax.plot(xs, ss, "--", lw=1.9, color=cols[k],
                        label=f"{job} {td} d")
            a2.plot(xo_m, so_cm / max(so_cm.max(), 1e-30), lw=2.2,
                    color=cols[k], label=f"observed {td} d")
            a2.plot(xs, ss / max(ss.max(), 1e-30), "--", lw=1.9,
                    color=cols[k], label=f"{job} {td} d")

        a0.set(xlabel="Distance along strike from injector (m)",
               ylabel="Cumulative slip (cm)", xlim=(-700, 700))
        a0.set_ylim(bottom=0)
        a0.set_title("A. Amplitude — solid observed, dashed simulated")
        a1.set(xlabel="Distance along strike from injector (m)",
               ylabel="Cumulative slip (cm)", xlim=(-700, 700), yscale="log",
               ylim=(1e-4, 20))
        a1.axhline(FRONT_THR * 100, color="#a8071a", ls=":", lw=1.5)
        a1.text(-660, FRONT_THR * 100 * 1.35,
                f"front threshold {FRONT_THR:.0e} m", fontsize=8,
                color="#a8071a")
        a1.set_title("B. Front — the threshold crossings nearly coincide")
        a2.set(xlabel="Distance along strike from injector (m)",
               ylabel="Slip / peak slip", xlim=(-700, 700), ylim=(0, 1.05))
        a2.set_title("C. Shape — each profile over its own peak")
        for ax in (a0, a1, a2):
            ax.grid(alpha=0.3, color=GRID)
            ax.legend(loc="upper right", ncol=1)
        fig.suptitle(f"Run {job}: tau_0 = "
                     f"{_ff(d['muinit'])*_ff(d['sigmainit']):.2f} MPa, "
                     f"kpmax {d.get('kpmax')}, a-b = "
                     f"{_ff(d['a'])-_ff(d['b']):+.3f}   —   observed slip from "
                     f"the seismicity catalog, 20 m west of the well",
                     fontsize=10.5)

        out = OUTROOT / str(job)
        out.mkdir(parents=True, exist_ok=True)
        for e in ("png", "pdf"):
            fig.savefig(out / f"obsoverlay_{job}.{e}", bbox_inches="tight")
        plt.close(fig)
        print(f"  {job}: wrote obsoverlay_{job}  (run ends {tend:.2f} d)")
        if dropped:
            print(f"    observed times not reached, DROPPED: {dropped}")
        print(f"    {'t':>4s} {'sim slip':>9s} {'obs peak':>9s} {'ratio':>7s} "
              f"{'sim R':>7s} {'obs R':>7s} {'R ratio':>8s}")
        for td, ta, ss, so, Rs, Ro in rows:
            print(f"    {td:>3d}d {ss:>8.3f}cm {so:>8.3f}cm "
                  f"{ss/max(so,1e-30):>7.2f} {Rs:>6.0f}m {Ro:>6.0f}m "
                  f"{Rs/max(Ro,1e-30):>8.2f}")


if __name__ == "__main__":
    main()
