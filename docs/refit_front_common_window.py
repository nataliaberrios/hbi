#!/usr/bin/env python3
"""Is the R-T fit coefficient change real, or a fitting-window artifact?

The notebook's R-T panel reports, for job 632510:

    pre-fix   downdip R = 0.185 sqrt(t)   perpendicular 0.170
    post-fix  downdip R = 0.173 sqrt(t)   perpendicular 0.159
    observed          R = 0.171 sqrt(t)

Read naively that says the fix moved the downdip fit onto the observed value.
But `fit_sqrt_front` least-squares fits over EVERY available point, so the fit
window is just however far the run got: 30.66 d pre-fix versus 16.37 d post-fix.
A sqrt-t coefficient is not invariant to that. Independent evidence says the
underlying fronts barely moved -- interpolated onto a common time grid the two
runs' front radii differ by under 1 m.

So: refit the PRE-fix run restricted to the post-fix window and compare like
with like. If pre-fix(truncated) ~ post-fix, the apparent improvement is an
artifact and the calibration conclusion is unchanged.

`extract_slip_cross_section` and `calculate_slip_front` are copied verbatim from
cell 76 of cooper_basin_validation_stage_june15.ipynb, and `fit_sqrt_front` from
cell 81, so the numbers are directly comparable to the notebook's.

Usage:  python refit_front_common_window.py
"""
import os
from pathlib import Path

import numpy as np

BUGGY = Path("/scratch/users/nberrios/3dhbi/output/632510")
FIXED = Path("/scratch/users/nberrios/3dhbi/fix_632510/output")

IMAX = JMAX = 601
NCELL = IMAX * JMAX
DS_KM = 0.005
DC = 1e-4


# ---- verbatim from the notebook, cell 76 -------------------------------------
def extract_slip_cross_section(slip, nt, imax, jmax, axis, ds_km):
    if axis == 'downdip':
        cs_index = int(jmax / 2)
        slip_cross_section = np.zeros((imax, nt))
        for t in range(nt):
            grid = slip[:, t].reshape(imax, jmax)
            slip_cross_section[:, t] = grid[:, cs_index]
        x_coords = np.linspace(-ds_km * imax / 2, ds_km * imax / 2, imax)
        n_points = imax
    else:
        cs_index = int(imax / 2)
        slip_cross_section = np.zeros((jmax, nt))
        for t in range(nt):
            grid = slip[:, t].reshape(imax, jmax)
            slip_cross_section[:, t] = grid[cs_index, :]
        x_coords = np.linspace(-ds_km * jmax / 2, ds_km * jmax / 2, jmax)
        n_points = jmax

    return slip_cross_section, x_coords, n_points, cs_index


def calculate_slip_front(slip_cross_section, x_coords, time_years, Dc, n_points):
    cross_spatial = []
    cross_temporal = []

    for i in range(n_points):
        above_dc = slip_cross_section[i, :] > Dc
        if np.any(above_dc):
            first_cross = np.argmax(above_dc)
            cross_spatial.append(abs(x_coords[i]))
            cross_temporal.append(time_years[first_cross])

    cross_spatial = np.array(cross_spatial)
    cross_temporal = np.array(cross_temporal) * 365.0  # days

    return cross_spatial, cross_temporal


# ---- verbatim from the notebook, cell 81 -------------------------------------
def fit_sqrt_front(x, R):
    if len(x) < 2:
        return np.nan, None
    b = np.sqrt(x)
    lam = np.sum(b * R) / np.sum(b**2)
    return lam, lam * b


def load(base, label):
    t_days = np.loadtxt(base / "time632510.dat")[:, 1] / 86400.0
    p = base / "slip632510.dat"
    nt = min(os.path.getsize(p) // (8 * NCELL), len(t_days))
    # the notebook works with slip as (ncell, nt); the file is (nt, ncell)
    slip = np.memmap(p, np.float64, "r", shape=(nt, NCELL)).T
    print(f"{label:18s} {nt} frames, 0-{t_days[nt-1]:.3f} d")
    return slip, nt, t_days[:nt] / 365.0      # time_years, as the notebook uses


def fronts(slip, nt, time_years):
    out = {}
    for axis in ("downdip", "perpendicular"):
        cs, x, n, _ = extract_slip_cross_section(slip, nt, IMAX, JMAX, axis, DS_KM)
        R, T = calculate_slip_front(cs, x, time_years, DC, n)
        out[axis] = (R, T)
    return out


print("=" * 72)
slip_pre, nt_pre, ty_pre = load(BUGGY, "pre-fix")
slip_post, nt_post, ty_post = load(FIXED, "post-fix")
t_cut = ty_post[-1] * 365.0
print(f"\npost-fix window ends at {t_cut:.3f} d -- refitting pre-fix to the same cut\n")

f_pre = fronts(slip_pre, nt_pre, ty_pre)
f_post = fronts(slip_post, nt_post, ty_post)

print(f"{'axis':>15s} {'pre-fix full':>14s} {'pre-fix <=cut':>14s} {'post-fix':>12s}")
print("-" * 72)
for axis in ("downdip", "perpendicular"):
    R, T = f_pre[axis]
    lam_full, _ = fit_sqrt_front(T, R)
    m = T <= t_cut
    lam_cut, _ = fit_sqrt_front(T[m], R[m])
    Rp, Tp = f_post[axis]
    lam_post, _ = fit_sqrt_front(Tp, Rp)
    print(f"{axis:>15s} {lam_full:14.4f} {lam_cut:14.4f} {lam_post:12.4f}")
    print(f"{'':>15s} {'':>14s} {'-> post/pre_cut = ':>14s}{lam_post/lam_cut:6.4f}")

print("\nobserved (from the notebook): 0.171")
print("\nIf pre-fix<=cut is close to post-fix, the apparent improvement from 0.185")
print("to 0.173 is the shorter fitting window, not the limitsigma fix.")
