#!/usr/bin/env python3
"""Which existing runs does the limitsigma ratchet actually change?

Two questions, because "rerun everything" is expensive and probably wrong:

PART A -- for the one run where we have both binaries (632510), how much did the
bug move the quantities being calibrated against? Slip front R(t) and injector
pore pressure are the observables; near-injector cumulative slip is the thing the
cusp lives in.

PART B -- across every recorded run, did the ratchet fire at all? The clamp can
only corrupt a run if it pushed total normal stress above its initial value.
monitor col 6 is maxnorm and col 6 at k=0 is sigmainit, so
max(maxnorm)/maxnorm[0] > 1 is a deck-independent detector. Runs that never
exceed 1.0 are untouched by this bug and need no rerun.

Usage:  python triage_ratchet_impact.py
"""
import os
import sys
from pathlib import Path

import numpy as np

BUGGY = Path("/scratch/users/nberrios/3dhbi/output/632510")
FIXED = Path("/scratch/users/nberrios/3dhbi/fix_632510/output")
OUTPUT_ROOT = Path("/scratch/users/nberrios/3dhbi/output")

IMAX = JMAX = 601
NCELL = IMAX * JMAX
IWELL = JWELL = 300
WELL = IWELL * JMAX + JWELL
DS_M = 5.0
THRESH = 1e-4               # m, = dc; same slip-front definition as the figures


def field(base, name, nt=None):
    p = base / f"{name}632510.dat"
    n = os.path.getsize(p) // (8 * NCELL)
    return np.memmap(p, np.float64, "r", shape=(n, NCELL))


print("=" * 74)
print("PART A -- how far did the bug move the calibration targets? (job 632510)")
print("=" * 74)

t_pre = np.loadtxt(BUGGY / "time632510.dat")[:, 1] / 86400.0
t_post = np.loadtxt(FIXED / "time632510.dat")[:, 1] / 86400.0

# Frames are written every `interval` steps OR on a tout crossing, so frame k in
# one run is NOT at the same physical time as frame k in the other -- pre-fix
# frame 1378 sits at ~30 d while post-fix frame 1378 sits at 16.4 d. Comparing by
# index compares different times and manufactures huge spurious differences
# (it reported R(t) off by 141 m when the true answer is far smaller).
# Everything below is interpolated onto a common TIME grid.
t_common = min(t_pre[-1], t_post[-1])
grid = np.linspace(0.0, t_common, 1200)
print(f"common window 0-{t_common:.3f} d, compared on {grid.size} interpolated times")
print(f"  (pre-fix {len(t_pre)} frames to {t_pre[-1]:.2f} d, "
      f"post-fix {len(t_post)} frames to {t_post[-1]:.2f} d)\n")
n_pre, n_post = len(t_pre), len(t_post)

ii, jj = np.meshgrid(np.arange(IMAX), np.arange(JMAX), indexing="ij")
rad = np.hypot(ii - IWELL, jj - JWELL) * DS_M

slip_pre, slip_post = field(BUGGY, "slip"), field(FIXED, "slip")
pf_pre, pf_post = field(BUGGY, "pf"), field(FIXED, "pf")

def front(slip, nt):
    R = np.empty(nt)
    for k in range(nt):
        m = np.asarray(slip[k]).reshape(IMAX, JMAX) > THRESH
        R[k] = rad[m].max() if m.any() else 0.0
    return R


R_pre = np.interp(grid, t_pre[:n_pre], front(slip_pre, n_pre))
R_post = np.interp(grid, t_post[:n_post], front(slip_post, n_post))
pf_w_pre = np.interp(grid, t_pre[:n_pre], np.asarray(pf_pre[:n_pre, WELL]))
pf_w_post = np.interp(grid, t_post[:n_post], np.asarray(pf_post[:n_post, WELL]))
sl_w_pre = np.interp(grid, t_pre[:n_pre], np.asarray(slip_pre[:n_pre, WELL]))
sl_w_post = np.interp(grid, t_post[:n_post], np.asarray(slip_post[:n_post, WELL]))


def report(name, a, b, unit, scale=1.0):
    d = np.abs(b - a) * scale
    denom = np.maximum(np.abs(a) * scale, 1e-12)
    rel = 100.0 * d / denom
    # ignore the first frames where both are ~0 and relative error is meaningless
    sel = (np.abs(a) * scale) > 0.01 * np.abs(a * scale).max()
    print(f"{name:34s} max abs diff {d.max():9.4f} {unit:5s}"
          f"   max rel diff {rel[sel].max() if sel.any() else 0.0:7.3f}%")


report("slip-front radius R(t)", R_pre, R_post, "m")
report("pore pressure at injector", pf_w_pre, pf_w_post, "MPa")
report("cumulative slip at injector", sl_w_pre, sl_w_post, "cm", 100.0)

print("\nslip-front position at selected times (m):")
print(f"  {'t (d)':>8s} {'pre-fix':>10s} {'post-fix':>10s} {'post-pre':>10s}")
for tt in (1.0, 2.0, 4.0, 8.0, 12.0, 16.0, t_common):
    if tt > t_common:
        continue
    a = np.interp(tt, grid, R_pre)
    b = np.interp(tt, grid, R_post)
    print(f"  {tt:8.2f} {a:10.1f} {b:10.1f} {b-a:+10.2f}")

# nearest frame to t_common in each run, again by TIME not index
s_pre = np.asarray(slip_pre[int(np.argmin(np.abs(t_pre[:n_pre] - t_common)))]).reshape(IMAX, JMAX)
s_post = np.asarray(slip_post[int(np.argmin(np.abs(t_post[:n_post] - t_common)))]).reshape(IMAX, JMAX)
print(f"\nat t = {t_common:.2f} d, whole-field cumulative slip:")
print(f"  max |post - pre| = {np.abs(s_post - s_pre).max()*100:.4f} cm "
      f"({100*np.abs(s_post-s_pre).max()/s_pre.max():.2f}% of peak slip)")
ring = rad > 300
print(f"  beyond 300 m of the injector: max |post - pre| = "
      f"{np.abs(s_post - s_pre)[ring].max()*100:.4f} cm")

print("\n" + "=" * 74)
print("PART B -- did the ratchet fire in each recorded run?")
print("=" * 74)
print("ratio = max(maxnorm) / maxnorm[0]; 1.00 means normal stress never rose")
print("above its initial value, i.e. the bug could not have altered that run.\n")

rows = []
for d in sorted(OUTPUT_ROOT.iterdir()):
    if not d.is_dir():
        continue
    mons = list(d.glob("monitor*.dat"))
    if not mons:
        continue
    m = mons[0]
    if os.path.getsize(m) == 0:
        continue
    try:
        # col 6 (1-indexed) = maxnorm; stream it, these files reach ~6 MB
        vals = []
        with open(m) as fh:
            for line in fh:
                p = line.split()
                if len(p) >= 7:
                    vals.append(float(p[5]))
        if not vals:
            continue
        v = np.array(vals)
        rows.append((d.name, v[0], v.max(), v.max() / v[0] if v[0] else np.nan, len(v)))
    except Exception as e:
        print(f"  [skip] {d.name}: {e}")

rows.sort(key=lambda r: -(r[3] if r[3] == r[3] else 0))
print(f"{'job':>10s} {'sigmainit':>10s} {'max sigma':>10s} {'ratio':>8s} {'steps':>8s}")
fired = 0
for name, v0, vmax, ratio, ns in rows:
    flag = ""
    if ratio > 1.001:
        fired += 1
        flag = "  <-- ratchet fired"
    print(f"{name:>10s} {v0:10.3f} {vmax:10.3f} {ratio:8.3f} {ns:8d}{flag}")

print(f"\n{fired} of {len(rows)} runs show the ratchet firing.")
print(f"{len(rows)-fired} runs never exceeded sigmainit and are unaffected by this bug.")
