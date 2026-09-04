#!/usr/bin/env python3
"""Measure alpha and deltaP from a constant-rate HBI run, for the GC solution.

The Gauss-Chebyshev slip solution takes six physical inputs. Four are read
straight off the deck and already agree with res632888.in exactly:

    f = f0 = 0.6            sigma_0 = sigmainit = 30e6
    mu = rigid = 24e9       tau_0 = muinit*sigmainit = 11.1e6

The other two are properties of the pressure field rather than deck entries, and
have to be measured:

    alpha   hydraulic diffusivity. Exact from the deck: kp/(eta*phi*beta).
            This is only well defined for permev F -- with enhancement, k and
            hence alpha vary in space AND time, and a fixed-alpha point-source
            solution has nothing to match. That is why res632888.in (permev F)
            is the right target and Job 911 (permev T, 250x range) is not.

    deltaP  amplitude of the point-source solution
                dp(r,t) = deltaP * E1(r^2 / (4*alpha*t))
            Analytically deltaP = (Q/w)*eta/(4*pi*k) for a fault zone of width w
            fed at rate Q, with the injection file holding Q/w directly. That
            prediction is printed, but the value reported for use is FITTED to
            the simulated field, so any discretisation or wellbore effect is
            included rather than assumed away.

The fit is a single linear least squares in deltaP at fixed alpha, over an
annulus that excludes both the injector cell (where the Peaceman well model, not
the point-source solution, sets the pressure) and the far field (where the value
is at the solver floor).

Usage:  python fit_pointsource.py 632888
"""
import argparse
import os
from pathlib import Path

import numpy as np
from scipy.special import exp1
import importlib.util as iu

H = Path("/home/users/nberrios/3dhbi/hbi_analysis")
OUT = H / "figures" / "taiyi_validation"
SCRATCH = "/scratch/users/nberrios/3dhbi/output"

_s = iu.spec_from_file_location("sf", str(H / "make_sweep_figures.py"))
sf = iu.module_from_spec(_s)
_s.loader.exec_module(sf)


def radial(arr, IM, JM, c, k):
    """Azimuthally averaged radial profile, which is what a point-source
    solution is a function of. Averaging over angle rather than taking one
    line also gives a direct check on axisymmetry."""
    f = np.asarray(arr[k]).reshape(IM, JM)
    ii, jj = np.mgrid[0:IM, 0:JM]
    rc = np.hypot(ii - c, jj - c)
    nb = int(rc.max()) + 1
    idx = rc.astype(int)
    tot = np.bincount(idx.ravel(), weights=f.ravel(), minlength=nb)
    cnt = np.bincount(idx.ravel(), minlength=nb)
    return tot[:nb] / np.maximum(cnt[:nb], 1), cnt[:nb]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("job", type=int)
    ap.add_argument("--times", type=float, nargs="*",
                    default=[3, 5, 7, 9, 11, 13, 15, 17])
    a = ap.parse_args()
    n = a.job
    dk = sf.deck(n)

    if dk.get("permev", "F").upper() == "T":
        print("WARNING: permev T. alpha is not constant in this run, so a "
              "fixed-alpha point-source solution does not apply. Numbers below "
              "are not usable for that comparison.")

    IM, JM = int(dk["imax"]), int(dk["jmax"])
    ds_m = sf.ffloat(dk["ds"]) * 1000.0
    NC = IM * JM
    kp, eta = sf.ffloat(dk["kp"]), sf.ffloat(dk["eta"])
    phi, beta = sf.ffloat(dk["phi"]), sf.ffloat(dk["beta"])
    alpha = kp / (eta * phi * beta)

    qfile = dk.get("injection_file")
    q_over_w = float(Path(f"/home/groups/edunham/nberrios/3dhbi/examples/"
                          f"grid_search_inputs/{qfile}").read_text()
                     .splitlines()[4].split()[0])
    dP_pred = q_over_w * eta / (4 * np.pi * kp)

    print(f"run {n}   permev {dk.get('permev','F')}   ds {ds_m:.0f} m   "
          f"kp {kp:.2e}")
    print(f"  f       = {sf.ffloat(dk['f0'])}")
    print(f"  sigma_0 = {sf.ffloat(dk['sigmainit'])*1e6:.4g} Pa")
    print(f"  tau_0   = {sf.ffloat(dk['muinit'])*sf.ffloat(dk['sigmainit'])*1e6:.4g} Pa")
    print(f"  mu      = {sf.ffloat(dk['rigid'])*1e9:.4g} Pa")
    print(f"  alpha   = {alpha:.6f} m^2/s          <- kp/(eta*phi*beta), exact")
    print(f"  q/w     = {q_over_w:.4e} m^2/s  from {qfile}")
    print(f"  deltaP predicted = (q/w)*eta/(4*pi*kp) = {dP_pred:.4e} Pa "
          f"= {dP_pred/1e6:.4f} MPa")
    print()

    p = f"{SCRATCH}/{n}/pf{n}.dat"
    t = np.atleast_2d(np.loadtxt(p.replace("pf", "time")))[:, 1] / 86400.0
    nt = min(os.path.getsize(p) // (8 * NC), len(t))
    pf = np.memmap(p, np.float64, "r", shape=(nt, NC))
    c = (IM - 1) // 2

    print(f"{'t (d)':>7s} {'deltaP fit MPa':>15s} {'rel resid':>10s} "
          f"{'fit/pred':>9s} {'n pts':>6s}")
    fits = []
    for td in a.times:
        if td > t[nt - 1] + 1e-9:
            continue
        k = int(np.argmin(abs(t[:nt] - td)))
        prof, cnt = radial(pf, IM, JM, c, k)
        rr = np.arange(len(prof)) * ds_m
        # dp is written in MPa; E1 argument in SI
        x = rr ** 2 / (4 * alpha * t[k] * 86400.0)
        basis = exp1(np.maximum(x, 1e-300))
        # Annulus: skip the injector cell and anything at the solver floor.
        m = (rr > 2 * ds_m) & (prof > 1e-3 * prof.max()) & np.isfinite(basis)
        if m.sum() < 10:
            continue
        dP = float(np.dot(basis[m], prof[m]) / np.dot(basis[m], basis[m]))
        resid = float(np.linalg.norm(prof[m] - dP * basis[m])
                      / np.linalg.norm(prof[m]))
        fits.append(dP)
        print(f"{t[k]:>7.2f} {dP:>15.4f} {resid:>10.4f} "
              f"{dP*1e6/dP_pred:>9.3f} {m.sum():>6d}")

    if fits:
        dP_MPa = float(np.median(fits))
        f0 = sf.ffloat(dk["f0"])
        sig0 = sf.ffloat(dk["sigmainit"]) * 1e6
        tau0 = sf.ffloat(dk["muinit"]) * sf.ffloat(dk["sigmainit"]) * 1e6
        print()
        print("USE THESE IN THE GC SCRIPT")
        print(f"  f       = {f0}")
        print(f"  sigma_0 = {sig0:.6g}")
        print(f"  tau_0   = {tau0:.6g}")
        print(f"  mu      = {sf.ffloat(dk['rigid'])*1e9:.6g}")
        print(f"  alpha   = {alpha:.6f}")
        print(f"  deltaP  = {dP_MPa*1e6:.6g}      # {dP_MPa:.4f} MPa, "
              f"median of the fits above")
        T = (f0 * sig0 - tau0) / (f0 * dP_MPa * 1e6)
        print(f"  -> T_final = (f*sigma_0 - tau_0)/(f*deltaP) = {T:.4f}")
        print(f"     (your script had T_final = 9.583 at deltaP = 1.2 MPa)")


if __name__ == "__main__":
    main()
