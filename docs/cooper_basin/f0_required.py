#!/usr/bin/env python3
"""What f0 would each run need to put its failure contour at the observed front?

Friction has never been swept: f0 = 0.6 in all 723 decks, (a,b) = (0.015,0.012)
in 721 of them. Every "friction sweep" in this project was actually a muinit
sweep, i.e. a sweep of initial SHEAR STRESS, which is a different object --
muinit is an initial condition applied once (main_LH.f90:861, tau = sigma*muinit)
while f0 is a material parameter passed into deriv() every timestep
(main_LH.f90:1908).

Why f0 is worth screening. The wall is that the wellhead and the front demand
opposite near-well pressures, and every lever tried so far (map contrast,
porosity grading, fluid/storage) moves BOTH targets, so they trade along a line
rather than spanning the plane. f0 is different: it does not appear anywhere in
m_diffusion.f90, so it cannot change the pressure field at all. It only changes
how much pressure is needed to fail. That makes it the one candidate for an
orthogonal lever -- IF the pressure reaches the observed front radius at all.

The screen. Under a Mohr-Coulomb reading, failure at radius r needs

    dp(r) >= dtauc = f0*sigmabar0 - tau0        =>   f0 <= (tau0 + dp(r))/sigmabar0

so f0_req = (tau0 + dp(R_OBS))/sigmabar0 is the LARGEST f0 that still fails out
to the observed front in that run's already-computed pressure field.

IMPORTANT CAVEAT, and the reason this is a screen and not a prediction. HBI's
regularised rate-and-state law has NO failure threshold:

    veltmp = 2*vref*exp(-psi/a)*sinh(tau/sigma/a)

is positive for any tau > 0, and f0 enters through the state variable psi, not as
a yield stress. So f0_req does not predict where HBI will put the front. It is
used here only to CHOOSE which f0 values are worth running; HBI's actual friction
law then decides. Importing Mohr-Coulomb into this model has already produced one
wrong conclusion in this project and must not be repeated.

Hard floor: f0 must exceed muinit. At f0 = muinit the fault is already at failure
with zero overpressure and the initial condition is meaningless.

Defensibility bands for f0, from F0_LITERATURE_REVIEW.md:
    0.60         unaltered granite (Taiyi's assumption; feldspar is 0.71)
    0.55         lightly chloritized
    0.50         ~20-35 wt.% chlorite -- the range MEASURED in the Pohang and
                 Gonghe EGS granite reservoirs (Zhang et al. 2022)
    0.45         strongly chloritized, beyond any measured EGS reservoir
    0.37         pure chlorite gouge -- absolute floor

Usage:  python f0_required.py
"""
import glob
import json
import os
from pathlib import Path

import importlib.util as iu
import numpy as np

H = Path("/home/users/nberrios/3dhbi/hbi_analysis")
OUT = Path("/home/users/nberrios/3dhbi/hbi_git/docs/figs/cooper_grid")

_s = iu.spec_from_file_location("sf", str(H / "make_sweep_figures.py"))
sf = iu.module_from_spec(_s)
_s.loader.exec_module(sf)

R_OBS = 417.0        # observed seismicity front radius at 5 d, m
WINDOW_D = 5.0
PASS_LO, PASS_HI = -15.0, 15.0        # wellhead bias band, percent
F0_DEFENSIBLE = 0.50                  # measured in real EGS granite reservoirs


def pf_path(n):
    for pat in (f"/scratch/users/nberrios/3dhbi/runs/*/output/pf{n}.dat",
                f"/scratch/users/nberrios/3dhbi/output/{n}/pf{n}.dat",
                f"/oak/stanford/groups/edunham/nberrios/3doutput/{n}/pf{n}.dat"):
        g = glob.glob(pat)
        if g:
            return g[0]
    return None


def analyse(n):
    dk = sf.deck(n)
    p = pf_path(n)
    if p is None:
        return None
    IM, JM = int(dk["imax"]), int(dk["jmax"])
    ds_m, NC = sf.ffloat(dk["ds"]) * 1000.0, IM * JM
    sig, mu, f0 = (sf.ffloat(dk["sigmainit"]), sf.ffloat(dk["muinit"]),
                   sf.ffloat(dk["f0"]))
    tau0 = mu * sig
    t = np.atleast_2d(np.loadtxt(p.replace("pf", "time")))[:, 1] / 86400.0
    nt = min(os.path.getsize(p) // (8 * NC), len(t))
    if nt < 5:
        return None
    pf = np.memmap(p, np.float64, "r", shape=(nt, NC))
    c = (IM - 1) // 2
    rr = np.arange(IM - c) * ds_m
    # Take the frame that maximises dp at the observed radius within the window:
    # the front is a running maximum, so the best chance of failing out to R_OBS
    # is whenever dp(R_OBS) peaks, not necessarily at the final time.
    kk = [k for k in range(nt) if t[k] <= WINDOW_D]
    best = None
    for k in kk:
        prof = np.asarray(pf[k]).reshape(IM, JM)[c, c:]
        dpr = float(np.interp(R_OBS, rr, prof))
        if best is None or dpr > best[1]:
            best = (k, dpr, float(prof[0]))
    k, dp_r, dp_w = best
    dtauc = f0 * sig - tau0
    prof = np.asarray(pf[k]).reshape(IM, JM)[c, c:]
    below = np.where(prof < dtauc)[0]

    def reach(level):
        """Radius at which dp falls below `level` -- how far the pressure gets."""
        b = np.where(prof < level)[0]
        return float(rr[b[0]]) if len(b) else float(rr[-1])

    return dict(n=n, sig=sig, mu=mu, tau0=tau0, dtauc=dtauc,
                t_best=float(t[k]), dp_well=dp_w, dp_robs=dp_r,
                f0_req=(tau0 + dp_r) / sig,
                r_fail=float(rr[below[0]]) if len(below) else float(rr[-1]),
                reach_1MPa=reach(1.0), reach_0p1MPa=reach(0.1),
                dp_max_domain=float(prof.max()))


def main():
    scores = {r["n"]: r for r in json.load(open(OUT / "grid_scores.json"))}
    rows = []
    for n in sorted(scores):
        a = analyse(n)
        if a is None:
            continue
        s = scores[n]
        a["p_pct"] = s.get("p_pct")
        a["lam_ratio"] = s.get("lam_ratio")
        rows.append(a)

    print(f"observed front radius R_OBS = {R_OBS:.0f} m;  window 0-{WINDOW_D:.0f} d")
    print("dp columns are at the frame that MAXIMISES dp(417 m) within the window.")
    print("f0_req = (tau0 + dp(417m))/sigmabar0, the largest f0 whose Mohr-Coulomb")
    print("  contour still reaches R_OBS. NOT a prediction -- HBI has no threshold.")
    print("reach = radius where dp drops below 1.0 / 0.1 MPa: profile steepness.\n")
    hdr = (f"{'run':>7s} {'tau0':>6s} {'dtauc':>6s} {'dp(rw)':>7s} {'dp(417)':>8s} "
           f"{'1MPa':>6s} {'0.1MPa':>7s} {'dtauc r':>8s} {'f0_req':>7s} "
           f"{'-mu':>6s} {'wellhd%':>8s} {'lam':>5s}")
    print(hdr)
    print("-" * len(hdr))
    for r in sorted(rows, key=lambda x: -x["dp_robs"]):
        pp = f"{r['p_pct']:+8.1f}" if r["p_pct"] is not None else "      --"
        lr = f"{r['lam_ratio']:5.2f}" if r["lam_ratio"] is not None else "   --"
        print(f"{r['n']:>7d} {r['tau0']:>6.2f} {r['dtauc']:>6.2f} {r['dp_well']:>7.2f} "
              f"{r['dp_robs']:>8.3f} {r['reach_1MPa']:>5.0f}m {r['reach_0p1MPa']:>6.0f}m "
              f"{r['r_fail']:>7.0f}m {r['f0_req']:>7.3f} "
              f"{r['f0_req'] - r['mu']:>+6.3f} {pp} {lr}")

    print(f"\n{len(rows)} runs with a pressure field.\n")

    # The headline. f0 can only matter where there IS overpressure; the margin it
    # would have to buy is dtauc - dp(R_OBS), and the room it has is
    # (f0 - muinit)*sigmabar0. Compare them.
    best = max(rows, key=lambda x: x["dp_robs"])
    print(f"MAX overpressure at 417 m over ALL {len(rows)} runs: "
          f"{best['dp_robs']:.3f} MPa (run {best['n']}, wellhead "
          f"{best['p_pct']:+.1f}%, lambda {best['lam_ratio']:.2f}).")
    print(f"  Its dp is {best['dp_well']:.1f} MPa at the well and "
          f"{best['dp_robs']:.3f} MPa at 417 m -- a factor of "
          f"{best['dp_well']/max(best['dp_robs'],1e-9):.0f} over "
          f"{R_OBS:.0f} m.")
    n_tiny = sum(1 for r in rows if r["dp_robs"] < 0.5)
    print(f"  {n_tiny}/{len(rows)} runs have dp(417 m) < 0.5 MPa.")
    print(f"  f0_req exceeds muinit by at most "
          f"{max(r['f0_req']-r['mu'] for r in rows):+.4f} across the whole grid.")
    print()
    print("READ THIS BEFORE USING f0_req. It collapses onto muinit because")
    print("dp(417 m) ~ 0 everywhere, so the Mohr-Coulomb screen degenerates: it")
    print("says the fault would have to be AT failure with no overpressure. That")
    print("is a statement about the SCREEN, not about f0. HBI's regularised")
    print("rate-and-state law has no threshold, and the runs prove the point --")
    print("632875 fits the front (lambda 1.03, slip front ~431 m at 5 d) while its")
    print("dtauc contour sits at 365 m and dp(417 m) is 0.34 MPa. Its front is")
    print("carried by ELASTIC stress transfer from the slipping patch, not by")
    print("pressure at the front. So f0 has to be swept in HBI, not screened here.")
    json.dump(rows, open(OUT / "f0_required.json", "w"), indent=1)
    print(f"\nwrote {OUT / 'f0_required.json'}")


if __name__ == "__main__":
    main()
