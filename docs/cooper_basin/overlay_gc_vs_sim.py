#!/usr/bin/env python3
"""The 17-day curve from job911style_632892.png, overlaid with the GC solution.

Plots exactly the last curve of that figure -- the along-strike line through the
injector at t = 17 d, in cm over +/-1.5 km, which is what the figure draws -- and
the Gauss-Chebyshev point-source solution on the same axes.

The GC solution is reproduced from the notebook cell verbatim: same g(), same
T(lambda) inversion, same n = 200 Gauss-Chebyshev grid, same
dim_factor = R_t*f*deltaP/mu. Parameters as the notebook sets them.

Usage:  python overlay_gc_vs_sim.py [--job 632892] [--days 17] [--n 200]
"""
import argparse
import os
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.special import ellipe as E, ellipk as K, exp1
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("/home/users/nberrios/3dhbi/hbi_analysis/figures/taiyi_validation")
SCRATCH = "/scratch/users/nberrios/3dhbi/output"
DECKS = "/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs"

# Notebook values.
F, SIGMA_0, TAU_0, DELTAP, ALPHA, MU = 0.6, 30e6, 11.1e6, 1.2e6, 1.25, 24e9


def g(xi, lam, eps=0.0):
    return exp1(lam ** 2 * xi ** 2) - exp1(lam ** 2) + lam ** -2 * np.exp(-lam ** 2) - eps


def lambda_from_T_factory():
    lams = np.logspace(-2, 2, 600)
    Ts = np.array([quad(lambda xi: g(xi, l) * xi / np.sqrt(1 - xi ** 2), 0, 1)[0]
                   for l in lams])
    o = np.argsort(Ts)
    inv = interp1d(np.log(Ts[o]), np.log(lams[o]), bounds_error=True)
    return lambda T: float(np.exp(inv(np.log(T))))


def gc_solution(t_sec, n=200):
    """Returns (r_full_m, slip_full_m, lam, R_t). Same algebra as the notebook."""
    a_idx, b_idx = np.arange(1, n), np.arange(1, n + 1)
    rt = np.cos(np.pi * a_idx / n)
    st = np.cos(np.pi * (b_idx - 0.5) / n)
    rbar, sbar = 0.5 * (rt + 1.0), 0.5 * (st + 1.0)
    ros = (rt[:, None] + 1.0) / (st[None, :] + 1.0)
    m = np.minimum(4.0 * ros / (1.0 + ros) ** 2, np.nextafter(1.0, 0.0))
    A = (-1.0 / n) * (K(m) / (st[None, :] + rt[:, None] + 2.0)
                      + E(m) / (st[None, :] - rt[:, None]))
    T_final = (F * SIGMA_0 - TAU_0) / (F * DELTAP)
    lam = lambda_from_T_factory()(T_final)
    gamma_i = g(rbar, lam, 0.0) - T_final
    R_t = lam * np.sqrt(4.0 * ALPHA * t_sec)
    dim_factor = R_t * F * DELTAP / MU
    orow = np.zeros((1, n))
    orow[0, -1] = 1.0 / np.sqrt(1.0 - st[-1] ** 2)
    phi = np.linalg.solve(np.vstack([A, orow]),
                          np.concatenate([gamma_i, [-2.0]]))
    dt_ = -(np.pi / n) * np.concatenate([[0.0], np.cumsum(phi)[:-1]])
    r_phys, slip = sbar * R_t, dim_factor * dt_
    return (np.concatenate([-r_phys, r_phys[::-1]]),
            np.concatenate([slip, slip[::-1]]), lam, R_t, T_final)


def sim_curve(job, days, how="azimuthal"):
    """Radial slip profile. how = "azimuthal" | "strike" | "dip".

    WHICH INDEX IS WHICH. coordinate3ddip (main_LH.f90:1577) fills
        do i=1,imax; do j=1,jmax; k=k+1
          xcol(k)=(i-imax/2-0.5)*ds0        <- x depends on i alone
          ycol(k),zcol(k) from yr(j),zr(j)  <- y,z depend on j alone, down-dip
    with j innermost, so reshape(IM,JM)[a,b] has a = i = STRIKE and
    b = j = DIP. frame[c,:] is therefore along DIP, and frame[:,c] along
    STRIKE. An earlier version of this script plotted frame[c,:] and labelled
    the axis "along-strike", which was backwards, and make_slip_figures.py has
    the same two labels swapped.

    It is not a cosmetic mislabel: the crack is elliptical, 875 m along strike
    against 810 m along dip at 17 d, so the two directions give lambda 0.3229
    and 0.2989. The azimuthal mean, 0.3192, is the right comparison for an
    axisymmetric solution and is the default here.
    """
    d = {}
    for line in open(f"{DECKS}/res{job}.in"):
        if line.startswith("!") or not line.strip():
            continue
        w = line.split()
        if len(w) >= 2:
            d[w[0]] = " ".join(w[1:]).strip('"')
    IM, JM = int(d["imax"]), int(d["jmax"])
    ds_m = float(d["ds"]) * 1000.0
    NC = IM * JM
    p = f"{SCRATCH}/{job}/slip{job}.dat"
    t = np.atleast_2d(np.loadtxt(p.replace("slip", "time")))[:, 1] / 86400.0
    nt = min(os.path.getsize(p) // (8 * NC), len(t))
    arr = np.memmap(p, np.float64, "r", shape=(nt, NC))
    c = (IM - 1) // 2
    k = int(np.argmin(np.abs(t[:nt] - days)))
    frame = np.asarray(arr[k]).reshape(IM, JM)
    if how == "strike":
        half = frame[c:, c]
    elif how == "dip":
        half = frame[c, c:]
    elif how == "azimuthal":
        ii, jj = np.mgrid[0:IM, 0:JM]
        idx = np.hypot(ii - c, jj - c).astype(int)
        nb = idx.max() + 1
        cnt = np.bincount(idx.ravel(), minlength=nb).astype(float)
        half = (np.bincount(idx.ravel(), weights=frame.ravel(), minlength=nb)
                / np.maximum(cnt, 1))
    else:
        raise ValueError(how)
    r_half = np.arange(len(half)) * ds_m
    return (np.concatenate([-r_half[::-1], r_half]),
            np.concatenate([half[::-1], half]), float(t[k]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--job", type=int, default=632892)
    ap.add_argument("--days", type=float, default=17.0)
    ap.add_argument("--n", type=int, default=200)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    r_gc, slip_gc, lam, R_t, T_final = gc_solution(a.days * 86400.0, a.n)
    curves = {h: sim_curve(a.job, a.days, h) for h in ("azimuthal", "strike", "dip")}
    r_sim, slip_sim, t_act = curves["azimuthal"]

    fig, ax = plt.subplots(figsize=(9.0, 5.2), dpi=150, constrained_layout=True)
    # strike and dip bracket the azimuthal mean: the crack is elliptical, so
    # showing only one direction hides a 1.9%/-5.7% spread in lambda.
    ax.fill_between(curves["strike"][0] / 1000.0,
                    np.interp(curves["strike"][0], curves["dip"][0],
                              curves["dip"][1]) * 100.0,
                    curves["strike"][1] * 100.0,
                    color="#1a1a19", alpha=0.15, lw=0,
                    label="HBI, along-strike to along-dip")
    ax.plot(r_sim / 1000.0, slip_sim * 100.0, "-", color="#1a1a19", lw=2.2,
            label=f"HBI {a.job}, azimuthal mean, t = {t_act:.2f} d")
    ax.scatter(r_gc / 1000.0, slip_gc * 100.0, s=13, facecolors="#a8071a",
               edgecolors="none", zorder=3,
               label=f"GC analytical, $\\lambda$ = {lam:.4f}")
    ax.axvline(0, color="b", ls="--", lw=1, label="injection point")
    ax.set(xlabel="Distance from injector (km)", ylabel="Cumulative Slip (cm)",
           xlim=(-1.5, 1.5))
    ax.set_ylim(bottom=0)
    ax.set_title(f"GC analytical solution vs HBI {a.job} at t = {a.days:g} days\n"
                 f"T = {T_final:.3f}, $\\lambda$ = {lam:.4f}, "
                 f"$R_t$ = {R_t:.0f} m,  $\\alpha$ = {ALPHA}, "
                 f"$\\Delta P$ = {DELTAP/1e6:g} MPa", fontsize=11)
    ax.legend(frameon=True, fontsize=9.5)
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"gc_vs_sim_{a.job}_{a.days:g}d.{e}",
                    bbox_inches="tight")
    plt.close(fig)

    pk_s, pk_g = slip_sim.max() * 100, np.nanmax(slip_gc) * 100
    # crack radius: GC's is R_t by construction; the sim's is where the plotted
    # line falls below a threshold, reported at two so it is not hidden
    print(f"wrote {OUT}/gc_vs_sim_{a.job}_{a.days:g}d.png")
    print(f"  t actual        {t_act:.4f} d")
    print(f"  peak slip  sim  {pk_s:.4f} cm   GC {pk_g:.4f} cm   "
          f"sim/GC {pk_s/pk_g:.4f}")
    print(f"  GC   lambda {lam:.4f}   R_t {R_t:.0f} m")
    al = 1.25
    for h in ("strike", "azimuthal", "dip"):
        rr, sl, _ = curves[h]
        half = sl[len(sl) // 2:]
        rh = rr[len(rr) // 2:]
        b = np.where(half < 1e-4)[0]
        R = rh[b[0]] if len(b) else rh[-1]
        lm = R / np.sqrt(4 * al * t_act * 86400)
        print(f"  sim  {h:<10s} R {R:>4.0f} m  lambda {lm:.4f}  "
              f"lam/GC {lm/lam:.4f}")


if __name__ == "__main__":
    main()
