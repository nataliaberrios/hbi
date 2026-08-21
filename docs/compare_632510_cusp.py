#!/usr/bin/env python3
"""Test 3 -- does removing the normal-stress ratchet remove the injector slip cusp?

The A/B test settles that the ratchet is gone. It cannot settle what the ratchet
*caused*, because its deck runs 2 days at 151x151 and no cusp forms in either
binary. This script compares job 632510's recorded pre-fix run against a rerun of
the identical deck with the post-fix binary, which is the direct test.

  pre-fix  (recorded)  /scratch/users/nberrios/3dhbi/output/632510/
  post-fix (rerun)     /scratch/users/nberrios/3dhbi/fix_632510/output/

The claim under test, from docs/normal_stress_clamp_ratchet.tex: under the bug the
injector cell is spuriously locked whenever pf recedes from a peak, leaving a local
MINIMUM of cumulative slip at the injection point with a ring of higher slip around
it. If that diagnosis is right, the rerun's slip profile should peak at the
injector instead.

Reported for each run, at a common time:
  - where the cumulative-slip maximum sits relative to the injector
  - the depth of the cusp: (ring peak - injector) / ring peak
  - the profile through the injector, plotted

Usage:  python compare_632510_cusp.py
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BUGGY = Path("/scratch/users/nberrios/3dhbi/output/632510")
FIXED = Path("/scratch/users/nberrios/3dhbi/fix_632510/output")
OUT = Path(__file__).resolve().parent / "figs"

JOBID = 632510
IMAX = JMAX = 601
DS_M = 5.0                      # ds 0.005 km
IWELL = JWELL = 300
NCELL = IMAX * JMAX

PRE, POST = "#D55E00", "#0072BD"     # validated pair, distinct under deuteranopia
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"


def load(base, label):
    tfile = base / f"time{JOBID}.dat"
    if not tfile.exists():
        sys.exit(f"missing {tfile}\nHas the rerun (job 36918035) finished?")
    t_days = np.loadtxt(tfile)[:, 1] / 86400.0
    slip = np.memmap(base / f"slip{JOBID}.dat", np.float64, "r")
    nt = min(slip.size // NCELL, t_days.size)
    slip = np.memmap(base / f"slip{JOBID}.dat", np.float64, "r", shape=(nt, NCELL))
    print(f"{label:8s}: {nt} frames, 0-{t_days[nt-1]:.3f} d")
    return dict(t=t_days[:nt], slip=slip, nt=nt, label=label)


def frame_at(run, t_target):
    """Index of the frame nearest t_target."""
    return int(np.argmin(np.abs(run["t"] - t_target)))


def cusp_metrics(field2d):
    """Where the slip maximum sits, and how deep the dip at the injector is."""
    mi, mj = np.unravel_index(np.argmax(field2d), field2d.shape)
    off_m = np.hypot(mi - IWELL, mj - JWELL) * DS_M
    inj = field2d[IWELL, JWELL]
    peak = field2d[mi, mj]
    depth_pct = 100.0 * (peak - inj) / peak if peak > 0 else 0.0
    return dict(mi=mi, mj=mj, off_m=off_m, inj=inj, peak=peak, depth_pct=depth_pct)


runs = [load(BUGGY, "pre-fix"), load(FIXED, "post-fix")]
t_common = min(r["t"][r["nt"] - 1] for r in runs)
print(f"\ncomparing at t = {t_common:.3f} d (the later time both runs reach)\n")

for r in runs:
    k = frame_at(r, t_common)
    r["frame"] = k
    s2d = np.array(r["slip"][k]).reshape(IMAX, JMAX)
    r["s2d"] = s2d
    m = cusp_metrics(s2d)
    r["m"] = m
    print(f"{r['label']:8s} (frame {k}, t = {r['t'][k]:.3f} d)")
    print(f"    slip max at ({m['mi']},{m['mj']}), {m['off_m']:.1f} m from the injector")
    print(f"    injector {m['inj']*100:.4f} cm   vs   max {m['peak']*100:.4f} cm")
    print(f"    cusp depth = {m['depth_pct']:.2f}% of the maximum\n")

pre, post = runs[0]["m"], runs[1]["m"]
print("-" * 68)
if pre["depth_pct"] > 1.0 and post["depth_pct"] < 0.25 * pre["depth_pct"]:
    print("CUSP REMOVED -- the pre-fix run has a real dip at the injector and the\n"
          "post-fix run does not. The ratchet is what produced it.")
elif pre["depth_pct"] <= 1.0:
    print("NO CUSP IN THE PRE-FIX RUN at this time -- this comparison cannot\n"
          "settle the question; check a later time or the profile figure.")
else:
    print(f"CUSP PERSISTS -- pre-fix dip {pre['depth_pct']:.2f}%, post-fix dip "
          f"{post['depth_pct']:.2f}%.\nRemoving the ratchet did NOT remove the cusp, so "
          "the cusp has another cause.")
print("-" * 68)

# ------------------------------------------------------------------ figure
dist = (np.arange(JMAX) - JWELL) * DS_M
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 4.4), constrained_layout=True)

for r, c in zip(runs, (PRE, POST)):
    ax.plot(dist, r["s2d"][IWELL] * 100, lw=2.0, color=c, label=r["label"])
ax.axvline(0, color=MUTED, lw=0.8, ls=":")
ax.set(xlabel="distance from injector (m)", ylabel="cumulative slip (cm)",
       title=f"Slip through the injector at t = {t_common:.2f} d")
ax.legend(frameon=False, fontsize=9, labelcolor=INK)

# zoom on the injector, where the cusp lives
sel = np.abs(dist) <= 300
for r, c in zip(runs, (PRE, POST)):
    ax2.plot(dist[sel], r["s2d"][IWELL][sel] * 100, lw=2.0, color=c,
             marker="o", ms=3, label=r["label"])
ax2.axvline(0, color=MUTED, lw=0.8, ls=":")
ax2.set(xlabel="distance from injector (m)", ylabel="cumulative slip (cm)",
        title="Zoom: is there a dip at the injection point?")
ax2.legend(frameon=False, fontsize=9, labelcolor=INK)

for a in (ax, ax2):
    a.grid(True, color=GRID, lw=0.6, alpha=0.7)
    a.set_axisbelow(True)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)

OUT.mkdir(parents=True, exist_ok=True)
stem = OUT / "cusp_632510_prefix_vs_postfix"
for ext in ("pdf", "png"):
    fig.savefig(f"{stem}.{ext}", dpi=200, bbox_inches="tight")
print(f"\nwrote {stem}.pdf/.png")
