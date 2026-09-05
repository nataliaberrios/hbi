#!/usr/bin/env python3
"""Fig 6 bottom panel with the injector cusp removed.

The poster's bottom panel (swarm26_4.pptx, ppt/media/image18.png) is run 632510
on PRE-FIX code, where the limitsigma ratchet produced a dimple at the injector
that deepens with time: 7.58% below the peak at 9 d, 4.94% at 17 d. The post-fix
rerun removes it -- 0.15% at 9 d -- and lives in a separate tree,
/scratch/users/nberrios/3dhbi/fix_632510/output/, which is why it does not turn
up when searching decks in grid_search_inputs.

TWO THINGS ABOUT THE ORIGINAL WORTH KNOWING BEFORE REUSING IT.

1. The axis label is wrong. In coordinate3ddip (main_LH.f90:1577) j is the
   innermost loop and controls y,z (down-dip), so reshape(IM,JM)[a,b] has
   a = i = STRIKE and b = j = DIP. image18 matches frame[c,:], the DIP profile:
   at 17 d that gives peak 20.487 cm and extent 770 m against the figure's
   ~20.5 cm and ~780 m, where the true strike profile gives 20.426 cm and 840 m.
   So the poster plots along-dip under an "along-strike" label, and the Fig 6
   caption repeats it. This script plots the same data -- dip, so the curves are
   comparable to the poster -- and labels it correctly.

   Dip is also the physically meaningful direction here: the fault is reverse, so
   slip is down-dip, along-dip growth is mode II and along-strike is mode III.

2. The post-fix run stops at 16.37 d, not 17. It hit nstep 40000, not tmax. The
   last curve is therefore labelled 16.4 d rather than silently drawn as "17
   days".

Usage:  python fig6_bottom_panel.py            # fixed panel alone
        python fig6_bottom_panel.py --compare  # also a pre/post overlay
"""
import argparse
import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PRE = Path("/scratch/users/nberrios/3dhbi/output/632510")
POST = Path("/scratch/users/nberrios/3dhbi/fix_632510/output")
OUT = Path("/home/users/nberrios/3dhbi/hbi_analysis/figures/fig6_remake")
JOB = 632510
IM = JM = 601
DS_M = 5.0
TIMES_D = [3, 5, 7, 9, 11, 13, 15, 17]      # as in the poster
XLIM_KM = 1.5


def load(base):
    p = base / f"slip{JOB}.dat"
    t = np.atleast_2d(np.loadtxt(base / f"time{JOB}.dat"))[:, 1] / 86400.0
    nt = min(os.path.getsize(p) // (8 * IM * JM), len(t))
    return np.memmap(p, np.float64, "r", shape=(nt, IM * JM)), t, nt


def profile(arr, k):
    """Along DIP through the injector, mirrored to a full diameter.

    frame[c, :] is the dip line -- see the module docstring. Returned in km and
    cm to match the poster's axes.
    """
    c = (IM - 1) // 2
    line = np.asarray(arr[k]).reshape(IM, JM)[c, :]
    x = (np.arange(JM) - c) * DS_M / 1000.0
    return x, line * 100.0


def times_available(t, nt):
    """Reference times a snapshot actually covers, plus the run's true end.

    Never silently promote the last frame to the requested time: the post-fix
    run ends at 16.37 d and drawing that as "17 days" would misstate it.
    """
    tend = t[nt - 1]
    keep = [td for td in TIMES_D if td <= tend + 1e-6]
    dropped = [td for td in TIMES_D if td > tend + 1e-6]
    if dropped:
        keep = keep + [tend]
    return keep, dropped, tend


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--compare", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    A_post, t_post, nt_post = load(POST)
    keep, dropped, tend = times_available(t_post, nt_post)
    print(f"post-fix: {nt_post} frames, 0-{tend:.2f} d")
    if dropped:
        print(f"  requested times not reached: {dropped}  "
              f"-> last curve drawn at the run's actual end, {tend:.2f} d")

    cmap = plt.get_cmap("viridis")
    cols = [cmap(v) for v in np.linspace(0.0, 0.95, len(keep))]

    fig, ax = plt.subplots(figsize=(10.5, 4.2), dpi=200, constrained_layout=True)
    for td, col in zip(keep, cols):
        k = int(np.argmin(np.abs(t_post[:nt_post] - td)))
        x, sl = profile(A_post, k)
        lab = (f"{td:.0f} days" if abs(td - round(td)) < 1e-6
               else f"{td:.1f} days (run end)")
        ax.plot(x, sl, lw=2.0, color=col, label=lab)
    ax.set(xlabel="Distance along-dip (km)", ylabel="Cumulative slip (cm)",
           xlim=(-XLIM_KM, XLIM_KM))
    ax.set_ylim(bottom=0)
    ax.legend(frameon=True, fontsize=10, loc="upper right")
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"fig6_bottom_postfix.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/fig6_bottom_postfix.png")

    # cusp depth, the quantity the fix is about
    A_pre, t_pre, nt_pre = load(PRE)
    print(f"\n{'t (d)':>7s} {'pre peak':>9s} {'pre dip%':>9s} "
          f"{'post peak':>10s} {'post dip%':>10s}")
    for td in keep:
        kp_ = int(np.argmin(np.abs(t_pre[:nt_pre] - td)))
        ko = int(np.argmin(np.abs(t_post[:nt_post] - td)))
        _, sp = profile(A_pre, kp_)
        _, so = profile(A_post, ko)
        c = len(sp) // 2
        print(f"{td:>7.2f} {sp.max():>9.3f} {100*(sp.max()-sp[c])/sp.max():>9.2f} "
              f"{so.max():>10.3f} {100*(so.max()-so[c])/so.max():>10.2f}")

    if a.compare:
        fig, ax = plt.subplots(figsize=(10.5, 4.2), dpi=200,
                               constrained_layout=True)
        for td in (9, 15, keep[-1]):
            kp_ = int(np.argmin(np.abs(t_pre[:nt_pre] - td)))
            ko = int(np.argmin(np.abs(t_post[:nt_post] - td)))
            x, sp = profile(A_pre, kp_)
            _, so = profile(A_post, ko)
            ax.plot(x, sp, lw=2.4, color="#D55E00", alpha=0.85,
                    label="pre-fix (poster)" if td == 9 else None)
            ax.plot(x, so, lw=1.6, color="#0072BD",
                    label="post-fix" if td == 9 else None)
        ax.set(xlabel="Distance along-dip (km)",
               ylabel="Cumulative slip (cm)", xlim=(-0.9, 0.9))
        ax.set_ylim(bottom=0)
        ax.set_title("Injector cusp, pre- vs post-fix, at 9 / 15 / "
                     f"{keep[-1]:.1f} d", fontsize=11)
        ax.legend(frameon=False, fontsize=10)
        for e in ("png", "pdf"):
            fig.savefig(OUT / f"fig6_bottom_prepost.{e}", bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {OUT}/fig6_bottom_prepost.png")


if __name__ == "__main__":
    main()
