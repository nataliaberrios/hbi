#!/usr/bin/env python3
"""Pre-fix vs post-fix overlays for job 632510.

compare_632510_cusp.py answers one question (did the injector cusp go away).
This is the broader set: every field worth seeing, with the recorded pre-fix run
and the post-fix rerun on the same axes.

  pre-fix  (recorded)  /scratch/users/nberrios/3dhbi/output/632510/       1402 frames, 0-30.66 d
  post-fix (rerun)     /scratch/users/nberrios/3dhbi/fix_632510/output/   1378 frames, 0-16.38 d

The two runs cover different spans, and the reason is NOT that the fix costs more
per step. While injection is on the two take near-identical timesteps (~60 s at
16 d). Injection stops at 17.15 d; past that the fault quiesces, steps grow to
~1e4 s, and the pre-fix run coasts 13.5 d to tmax on ~2000 steps. The post-fix run
exhausted nstep=40000 at 16.4 d, i.e. ~0.75 d short of that cheap tail. On top of
the physical wind-down the pre-fix run also OVER-locks: the ratchet drives normal
stress to 135 MPa, far above sigmainit=30, which the post-fix run never does.
Several figures exist to show this, so DO NOT "fix" the differing x-ranges.

Fields are memmapped, not loaded -- each of slip/sigma/pf/kp is ~4 GB per run.

Usage:  python compare_prepost_632510.py
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

BUGGY = Path("/scratch/users/nberrios/3dhbi/output/632510")
FIXED = Path("/scratch/users/nberrios/3dhbi/fix_632510/output")
OUT = Path(__file__).resolve().parent / "figs" / "prepost_632510"

JOBID = 632510
IMAX = JMAX = 601
DS_M = 5.0                       # ds 0.005 km
IWELL = JWELL = 300              # injector, Fortran (301,301) -> python (300,300)
NCELL = IMAX * JMAX
WELL = IWELL * JMAX + JWELL

SIGMAINIT = 30.0                 # MPa, from the deck
MINSIG = 1.0                     # MPa, the limitsigma floor
# june_clean.txt has 425 time nodes running to 31.39 d, but the rate column is
# ZERO from node 423 onward -- real injection stops at 17.15 d and the schedule
# is just padded out. Do not read the last time node as the end of injection.
INJ_END_D = 17.15                # last node with a nonzero injection rate

# validated with scripts/validate_palette.js --mode light: all six checks PASS
PRE, POST = "#D55E00", "#0072BD"
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
SEQ = "viridis"                  # magnitude: one hue ramp
DIV = "RdBu_r"                   # polarity: two poles + neutral midpoint

plt.rcParams.update({
    "font.size": 11, "axes.titlesize": 12.5, "axes.labelsize": 11,
    "axes.edgecolor": MUTED, "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
})


def style(ax):
    ax.grid(True, color=GRID, lw=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return ax


def ref_line(ax, y, text, side="left"):
    """Horizontal reference line with a label that does not sit on the data."""
    ax.axhline(y, color=MUTED, lw=1, ls="--")
    x, ha = (0.015, "left") if side == "left" else (0.985, "right")
    ax.annotate(text, xy=(x, y), xycoords=("axes fraction", "data"),
                xytext=(0, 6), textcoords="offset points", fontsize=9, color=MUTED,
                ha=ha, va="bottom",
                bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.75))


def save(fig, stem):
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {stem}.png/.pdf", flush=True)


def load(base, label):
    tfile = base / f"time{JOBID}.dat"
    if not tfile.exists():
        sys.exit(f"missing {tfile}")
    t_days = np.loadtxt(tfile)[:, 1] / 86400.0

    mon = np.loadtxt(base / f"monitor{JOBID}.dat")
    # cols (main_LH.f90:1328): k, x, log10(mvelG), meanslip, meanmu,
    #                          maxnorm, minnorm, errmax, dtdid, nrjct, walltime
    run = dict(label=label, t=t_days,
               mk=mon[:, 0], mt=mon[:, 1] / 86400.0, lvmax=mon[:, 2],
               meanslip=mon[:, 3], meanmu=mon[:, 4],
               maxnorm=mon[:, 5], minnorm=mon[:, 6], dt=mon[:, 8])

    for f in ("slip", "sigma", "pf", "kp"):
        p = base / f"{f}{JOBID}.dat"
        if not p.exists():
            print(f"  [warn] {label}: no {f} file")
            run[f] = None
            continue
        nt = min(np.memmap(p, np.float64, "r").size // NCELL, t_days.size)
        run[f] = np.memmap(p, np.float64, "r", shape=(nt, NCELL))
        run["nt"] = nt
    run["t"] = t_days[:run["nt"]]
    print(f"{label:8s}: {run['nt']} frames, 0-{run['t'][-1]:.3f} d, "
          f"{len(run['mk'])} monitor rows", flush=True)
    return run


def frame_at(run, t):
    return int(np.argmin(np.abs(run["t"] - t)))


def label_end(ax, x, y, text, color):
    ax.annotate(text, xy=(x, y), xytext=(5, 0), textcoords="offset points",
                color=color, fontsize=9, va="center", clip_on=False)


runs = [load(BUGGY, "pre-fix"), load(FIXED, "post-fix")]
pre, post = runs
COL = {"pre-fix": PRE, "post-fix": POST}
# The two runs coincide almost exactly over most of the common window. Drawn at
# equal width the post-fix line simply hides the pre-fix one and it reads as
# missing data, so draw pre-fix wide and post-fix narrow on top of it: where only
# the blue core shows with an orange rim, the runs agree.
LW = {"pre-fix": 3.4, "post-fix": 1.5}
t_common = min(r["t"][-1] for r in runs)
print(f"\ncommon window: 0-{t_common:.3f} d\n", flush=True)


# ----------------------------------------------------------------- 1. the ratchet
fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True, constrained_layout=True)
for r in runs:
    c = COL[r["label"]]
    axes[0].plot(r["mt"], r["maxnorm"], lw=LW[r["label"]], color=c, label=r["label"])
    axes[1].plot(r["mt"], r["minnorm"], lw=LW[r["label"]], color=c, label=r["label"])
ref_line(axes[0], SIGMAINIT, f"sigmainit = {SIGMAINIT:g} MPa", side="right")
axes[0].set(ylabel="max normal stress (MPa)",
            title="Normal-stress drift: post-fix never leaves its initial value")
ref_line(axes[1], MINSIG, f"minsig = {MINSIG:g} MPa", side="right")
axes[1].set(xlabel="time (days)", ylabel="min normal stress (MPa)",
            title="Minimum effective normal stress: pinned on the minsig floor until the pre-fix lockup")
for a in axes:
    style(a).legend(frameon=False, fontsize=9, labelcolor=INK)
save(fig, "01_normal_stress_drift")


# ------------------------------------------------------- 2. slip rate and lockup
fig, ax = plt.subplots(figsize=(9, 4.6), constrained_layout=True)
for r in runs:
    c = COL[r["label"]]
    ax.plot(r["mt"], r["lvmax"], lw=LW[r["label"]], color=c, label=r["label"])
ax.axvline(INJ_END_D, color=MUTED, lw=1, ls=":")
ax.annotate("injection rate -> 0", xy=(INJ_END_D, 0.02), xycoords=("data", "axes fraction"),
            xytext=(-4, 0), textcoords="offset points", fontsize=9, color=MUTED,
            rotation=90, ha="right")
style(ax).set(xlabel="time (days)", ylabel="log$_{10}$ peak slip rate (m/s)",
              title="Peak slip rate: the two runs agree; the tail is post-injection relaxation")
ax.legend(frameon=False, fontsize=9, labelcolor=INK, loc="lower right")
save(fig, "02_peak_slip_rate")


# ---------------------------------------------------------------- 3. the timestep
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
for r in runs:
    c = COL[r["label"]]
    axes[0].semilogy(r["mt"], r["dt"], lw=1.1, color=c, label=r["label"])
    axes[1].plot(r["mt"], r["mk"], lw=LW[r["label"]], color=c, label=r["label"])
ref_line(axes[0], 43200, "dtmax = 43200 s")
style(axes[0]).set(xlabel="time (days)", ylabel="timestep taken (s)",
                   title="Timestep: identical while injecting; pre-fix then coasts post-injection")
style(axes[1]).set(xlabel="time (days)", ylabel="step number",
                   title="Step budget: the pre-fix run's last 13.5 d cost ~2000 steps")
for a in axes:
    a.legend(frameon=False, fontsize=9, labelcolor=INK)
save(fig, "03_timestep_and_step_budget")


# ----------------------------------------------- 4. pressure and stress at injector
fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True, constrained_layout=True)
for r in runs:
    c = COL[r["label"]]
    t = r["t"]
    axes[0].plot(t, np.asarray(r["pf"][:, WELL]), lw=LW[r["label"]], color=c, label=r["label"])
    axes[1].plot(t, np.asarray(r["sigma"][:, WELL]), lw=LW[r["label"]], color=c, label=r["label"])
    if r["kp"] is not None:
        # kp at the injector is pinned at kpmax for the entire run in both cases --
        # the cell slips hard enough that the slip-driven term in the kp ODE
        # (main_LH.f90:2413) saturates immediately. Left to autoscale, matplotlib
        # subtracts an offset and magnifies float noise in the 9th significant
        # digit into a plausible-looking sawtooth, so pin the axis instead.
        kpw = np.asarray(r["kp"][:, WELL]) / 1e-14
        axes[2].plot(t, kpw, lw=LW[r["label"]], color=c, label=r["label"])
        print(f"  {r['label']} kp at injector: [{kpw.min():.6f}, {kpw.max():.6f}] x1e-14")
ref_line(axes[0], SIGMAINIT, f"sigmainit = {SIGMAINIT:g} MPa", side="right")
axes[0].set(ylabel="pore pressure (MPa)",
            title="Pore pressure at the injector -- the overshoot past sigmainit survives the fix")
ref_line(axes[1], MINSIG, f"minsig = {MINSIG:g} MPa", side="right")
axes[1].set(ylabel="normal stress (MPa)",
            title="Effective normal stress at the injector")
axes[2].set(xlabel="time (days)", ylabel="permeability (10$^{-14}$ m$^2$)",
            title="Permeability at the injector: pinned at kpmax throughout, in both runs",
            ylim=(0, 6.5))
axes[2].ticklabel_format(axis="y", useOffset=False, style="plain")
ref_line(axes[2], 5.0, "kpmax = 5e-14 m$^2$", side="right")
for a in axes:
    style(a).legend(frameon=False, fontsize=9, labelcolor=INK)
save(fig, "04_injector_pf_sigma_kp")


# --------------------------------------- 5. how much of the fault sits on the floor
print("scanning sigma fields for cells on the minsig floor ...", flush=True)
fig, ax = plt.subplots(figsize=(9, 4.6), constrained_layout=True)
for r in runs:
    c = COL[r["label"]]
    frac = np.empty(r["nt"])
    ratchet = np.empty(r["nt"])
    for k in range(r["nt"]):
        s = np.asarray(r["sigma"][k])
        frac[k] = np.mean(s <= MINSIG * 1.0001)
        ratchet[k] = np.mean(s > SIGMAINIT * 1.0001)
    r["frac_floor"], r["frac_ratchet"] = frac, ratchet
    ax.plot(r["t"], 100 * frac, lw=2, color=c, label=f"{r['label']} — on minsig floor")
    ax.plot(r["t"], 100 * ratchet, lw=1.6, color=c, ls="--",
            label=f"{r['label']} — drifted above sigmainit")
style(ax).set(xlabel="time (days)", ylabel="fault area (%)",
              title="Fraction of the fault clamped at minsig, and fraction drifted above sigmainit")
ax.legend(frameon=False, fontsize=9, labelcolor=INK)
save(fig, "05_clamp_and_drift_area")


# ------------------------------------------------- 6. slip profiles at matched times
times = [t for t in (1.0, 4.0, 8.0, 12.0, t_common) if t <= t_common]
fig, axes = plt.subplots(1, len(times), figsize=(3.4 * len(times), 3.8),
                         sharey=True, constrained_layout=True)
dist = (np.arange(JMAX) - JWELL) * DS_M
for a, tt in zip(np.atleast_1d(axes), times):
    for r in runs:
        k = frame_at(r, tt)
        prof = np.asarray(r["slip"][k]).reshape(IMAX, JMAX)[IWELL] * 100
        a.plot(dist, prof, lw=LW[r["label"]], color=COL[r["label"]], label=r["label"])
    a.axvline(0, color=MUTED, lw=0.8, ls=":")
    style(a).set(xlabel="distance from injector (m)", title=f"t = {tt:.2f} d")
    a.set_xlim(-900, 900)
np.atleast_1d(axes)[0].set_ylabel("cumulative slip (cm)")
np.atleast_1d(axes)[0].legend(frameon=False, fontsize=9, labelcolor=INK)
fig.suptitle("Slip through the injector -- the pre-fix cusp deepens with time", fontsize=13)
save(fig, "06_slip_profiles_vs_time")


# ------------------------------------------------------------- 7. slip front radius
print("computing slip fronts ...", flush=True)
ii, jj = np.meshgrid(np.arange(IMAX), np.arange(JMAX), indexing="ij")
rad_m = np.hypot(ii - IWELL, jj - JWELL) * DS_M
THRESH = 1e-4                                    # m, = dc
fig, ax = plt.subplots(figsize=(9, 4.6), constrained_layout=True)
for r in runs:
    c = COL[r["label"]]
    R = np.empty(r["nt"])
    for k in range(r["nt"]):
        m = np.asarray(r["slip"][k]).reshape(IMAX, JMAX) > THRESH
        R[k] = rad_m[m].max() if m.any() else 0.0
    ax.plot(r["t"], R, lw=LW[r["label"]], color=c, label=r["label"])
    label_end(ax, r["t"][-1], R[-1], r["label"], c)
style(ax).set(xlabel="time (days)", ylabel="slip-front radius (m)",
              title=f"Slip front (slip > {THRESH:g} m = dc)")
ax.legend(frameon=False, fontsize=9, labelcolor=INK, loc="lower right")
save(fig, "07_slip_front_radius")


# ------------------------------------------------------ 8. field maps at t_common
kpre, kpost = frame_at(pre, t_common), frame_at(post, t_common)
half = 250                                        # cells to either side of the well
sl = slice(IWELL - half, IWELL + half + 1)
ext = [-half * DS_M, half * DS_M, -half * DS_M, half * DS_M]

for fld, unit, scale in (("sigma", "MPa", 1.0), ("slip", "cm", 100.0), ("pf", "MPa", 1.0)):
    a_pre = np.asarray(pre[fld][kpre]).reshape(IMAX, JMAX)[sl, sl] * scale
    a_post = np.asarray(post[fld][kpost]).reshape(IMAX, JMAX)[sl, sl] * scale
    diff = a_post - a_pre
    vmin, vmax = min(a_pre.min(), a_post.min()), max(a_pre.max(), a_post.max())
    lim = np.abs(diff).max() or 1.0

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    for a, arr, ttl in ((axes[0], a_pre, "pre-fix"), (axes[1], a_post, "post-fix")):
        im = a.imshow(arr, origin="lower", extent=ext, cmap=SEQ, vmin=vmin, vmax=vmax)
        a.set(title=f"{ttl}  ({fld}, {unit})", xlabel="m from injector")
    fig.colorbar(im, ax=axes[1], label=f"{fld} ({unit})")
    imd = axes[2].imshow(diff, origin="lower", extent=ext, cmap=DIV,
                         norm=TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=lim))
    axes[2].set(title="post-fix - pre-fix", xlabel="m from injector")
    fig.colorbar(imd, ax=axes[2], label=f"difference ({unit})")
    axes[0].set_ylabel("m from injector")
    for a in axes:
        a.plot(0, 0, marker="+", ms=9, color="w", mew=1.6)
    fig.suptitle(f"{fld} at t = {t_common:.2f} d", fontsize=13)
    save(fig, f"08_{fld}_map")


# ------------------------------------------------------------------- 9. summary
print("\n" + "=" * 70)
print(f"{'':10s} {'end t (d)':>10s} {'steps':>8s} {'max sigma':>10s} {'log vmax':>9s} {'floor %':>8s}")
for r in runs:
    print(f"{r['label']:10s} {r['t'][-1]:10.2f} {int(r['mk'][-1]):8d} "
          f"{r['maxnorm'].max():10.2f} {r['lvmax'][-1]:9.2f} {100*r['frac_floor'][-1]:8.2f}")
print("=" * 70)
print(f"\nfigures -> {OUT}")
