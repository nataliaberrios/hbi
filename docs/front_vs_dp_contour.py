#!/usr/bin/env python3
"""Is the 632507 slip-front deficit just the failure threshold, not elasticity?

The pressure solvers already agree with Taiyi to ~0.1 MPa RMS, so dp(r,t) is
effectively the same field in both models.  If that is true, the slip front can
only differ through the condition that turns dp into slip.

Both decks set tau0 = 15 MPa.  They differ in effective normal stress:

    Taiyi   sigma_eff = 27.99 MPa  (s_v 100, s_Hmax 160, 10 deg dip, p_pore 73.82)
    632507  sigma_eff = 30.00 MPa  (sigmainit 30.0, muinit 0.5)

so the fault reaches f0 = 0.6 at dp = 2.99 MPa for Taiyi and dp = 5.00 MPa here.
The prediction is that the dp = 5 MPa contour tracks the simulated slip front
(0.079 sqrt(t)) while the dp = 3 MPa contour sits much further out, near the
observed 0.171 sqrt(t) -- which would make the gap a prestress difference and
exonerate the elasticity kernel.

Reads pf632507.dat by memmap, one frame at a time; no rerun needed.
"""
import os
from pathlib import Path

import numpy as np

RUN = Path("/scratch/users/nberrios/3dhbi/output/632507")
JOB = "632507"

IMAX = JMAX = 601
NCELL = IMAX * JMAX
IWELL = JWELL = 300
DS_M = 20.0                     # ds 0.020 km
PA_PER_MPA = 1e6

# dp thresholds to contour, MPa.  2.99 and 5.00 are the two models' failure
# margins; the rest map out the sensitivity.
LEVELS = [1.0, 2.0, 2.99, 4.0, 5.0, 6.0, 8.0]

t = np.loadtxt(RUN / f"time{JOB}.dat")[:, 1] / 86400.0     # days
pf = np.memmap(RUN / f"pf{JOB}.dat", np.float64, "r",
               shape=(os.path.getsize(RUN / f"pf{JOB}.dat") // (8 * NCELL), NCELL))
nt = min(len(t), len(pf))
t = t[:nt]

# radial distance of every cell from the injector, in metres
ii, jj = np.meshgrid(np.arange(IMAX), np.arange(JMAX), indexing="ij")
rad = np.hypot(ii - IWELL, jj - JWELL) * DS_M

print(f"{nt} frames, 0-{t[-1]:.3f} d;  grid {IMAX}x{JMAX} at ds {DS_M:.0f} m")
print("pf is dp directly (deck sets pfinit 0)\n")

# ---------------------------------------------------------------- contours
# R(t) for each level = furthest cell still above that dp.
R = np.zeros((len(LEVELS), nt))
for k in range(nt):
    frame = np.asarray(pf[k]).reshape(IMAX, JMAX) / PA_PER_MPA
    for li, lev in enumerate(LEVELS):
        m = frame > lev
        R[li, k] = rad[m].max() if m.any() else 0.0

print(f"peak dp at injector: {max(np.asarray(pf[k]).reshape(IMAX,JMAX)[IWELL,JWELL] for k in range(nt))/PA_PER_MPA:.3f} MPa\n")

# ------------------------------------------------------------------- fits
# R = C sqrt(t), C in km/sqrt(day), matching the notebook's R-T panel units.
print("  dp level     C in R = C sqrt(t)      R at 15 d")
print("   (MPa)          (km/sqrt(d))            (km)")
print("-" * 52)
for li, lev in enumerate(LEVELS):
    r_km = R[li] / 1000.0
    ok = (t > 0.05) & (r_km > 0)
    if ok.sum() < 5:
        print(f"  {lev:6.2f}        (never reached)")
        continue
    C = np.sum(np.sqrt(t[ok]) * r_km[ok]) / np.sum(t[ok])
    tag = ""
    if abs(lev - 5.00) < 1e-9:
        tag = "  <-- 632507 failure margin"
    if abs(lev - 2.99) < 1e-9:
        tag = "  <-- Taiyi failure margin"
    print(f"  {lev:6.2f}          {C:8.4f}            {C*np.sqrt(15):6.3f}{tag}")

print("\nfor reference: simulated slip front 0.079, observed seismicity 0.171")

np.savez(RUN.parent.parent / "fix_632510" / "front_vs_dp_632507.npz",
         t=t, R=R, levels=np.array(LEVELS))
print("\nsaved front_vs_dp_632507.npz")
print("FRONT-VS-DP DONE")
