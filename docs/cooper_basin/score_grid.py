#!/usr/bin/env python3
"""Score every grid run on one window and emit a single results table.

One place, one window, one set of conventions, so results arrive as a filled-in
grid rather than as a sequence of ad-hoc readings. Every convention here exists
because getting it wrong has already produced a wrong answer once:

  WINDOW  0-5 d for every run, with the observed front REFIT on that same window.
    lambda depends on the window: the four Stage 2 runs score 0.97-0.99 over
    0-8 d and 0.57-0.73 over 0-5 d, because lambda_obs is 0.1412 over 8 d and
    0.1866 over 5 d. A front ratio without its window is meaningless.

  PRESSURE  mean(sim)/mean(measured) - 1, over FLOWING periods only, plus RMS in
    MPa. Not mean(sim/measured): the measured wellhead falls through ~0 during
    shut-ins, so that ratio divides by near-zero and read 2.17 and 0.32 on two
    samplings of the same run. Flowing = injecting above 25% of peak rate, which
    is only ~14% of a 5 d window -- most of it is the 1.58-3.56 d shut-in, which
    HBI cannot follow at all because it has no wellbore bleed-off.

  PEAK SLIP / dc  reported alongside lambda. The front is the dc contour, so a
    run whose peak slip is only ~1.6x dc has a razor-thin ring and its lambda is
    near-noise. Below 1.0 there is no front at all -- reported as "no slip",
    which is a result, not missing data.

  CODE CHECK  max(maxnorm)/maxnorm[0] from monitor column 6 must be 1. Anything
    above means the run used the pre-fix limitsigma ratchet. bool() matters:
    numpy.False_ is False evaluates False, which previously mislabelled every
    buggy run as merely too short.

  BOUNDS CHECK  kpmax == map max and kpmin == map min. Of 112 archived runs, the
    35 obeying this had 29 reach dc; the 77 that did not had 2. Any run failing
    it is flagged, not silently scored.

Usage:  python score_grid.py 632800 632801 ...
        python score_grid.py --stage 1
        python score_grid.py --all          # every 6328xx deck with output
"""
import argparse
import glob
import importlib.util as iu
import json
import os
from pathlib import Path

import numpy as np

H = Path("/home/users/nberrios/3dhbi/hbi_analysis")
IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
OUT = Path("/home/users/nberrios/3dhbi/hbi_git/docs/figs/cooper_grid")
SCRATCH = "/scratch/users/nberrios/3dhbi/output"
OAK = "/oak/stanford/groups/edunham/nberrios/3doutput"
RUNS = "/scratch/users/nberrios/3dhbi/runs"

WINDOW_D = 5.0
F0 = 0.6
DP_OBS = 10.92          # peak measured downhole overpressure over 0-5 d, MPa

_spec = iu.spec_from_file_location("sf", str(H / "make_sweep_figures.py"))
sf = iu.module_from_spec(_spec)
_spec.loader.exec_module(sf)

STAGES = {1: list(range(632800, 632808)), 2: list(range(632810, 632822))}
EXTRA_STAGE2 = [632522, 632523, 632524, 632525]


def outdir(n):
    for b in (f"{SCRATCH}/{n}", f"{OAK}/{n}"):
        if os.path.exists(f"{b}/time{n}.dat"):
            return b, "final"
    g = glob.glob(f"{RUNS}/*/output/time{n}.dat")
    return (os.path.dirname(g[0]), "in progress") if g else (None, "no output")


def code_is_fixed(n, base):
    p = f"{base}/monitor{n}.dat"
    if not (os.path.exists(p) and os.path.getsize(p)):
        return None
    try:
        a = np.atleast_2d(np.loadtxt(p))
        # bool() is required: numpy.False_ is False evaluates False
        return bool(a[:, 5].max() / a[0, 5] <= 1.0001) if a.shape[1] >= 6 else None
    except Exception:
        return None


def bounds_ok(dk):
    """kpmax == map max and kpmin == map min, or a uniform field matching kp."""
    kx, km = sf.ffloat(dk.get("kpmax")), sf.ffloat(dk.get("kpmin"))
    if dk.get("parameterfromfile", "F").upper() not in ("T", "TRUE", ".TRUE."):
        kp = sf.ffloat(dk.get("kp"))
        if None in (kp, km):
            return None
        return abs(km - kp) / kp < 0.01          # uniform field starts at the floor
    pf = dk.get("parameter_file")
    if not pf or not (IN / pf).exists() or None in (kx, km):
        return None
    k = np.loadtxt(IN / pf, skiprows=1)
    return bool(abs(k.max() - kx) / kx < 0.01 and abs(k.min() - km) / km < 0.01)


def peak_slip(n, base, dk):
    sp = f"{base}/slip{n}.dat"
    if not (os.path.exists(sp) and os.path.getsize(sp)):
        return None
    IM, JM = int(dk["imax"]), int(dk["jmax"])
    NC = IM * JM
    nt = os.path.getsize(sp) // (8 * NC)
    if nt < 1:
        return None
    with open(sp, "rb") as f:               # last frame only: cheap
        f.seek((nt - 1) * NC * 8)
        return float(np.frombuffer(f.read(NC * 8), np.float64).max())


def score(n, obs, lam_obs):
    dk = sf.deck(n)
    base, where = outdir(n)
    r = dict(n=n, where=where,
             sigma=sf.ffloat(dk.get("sigmainit")), mu=sf.ffloat(dk.get("muinit")),
             eta=sf.ffloat(dk.get("eta")), beta=sf.ffloat(dk.get("beta")),
             phi=sf.ffloat(dk.get("phi")), permev=dk.get("permev", "F"),
             kpmax=dk.get("kpmax"), kpmin=dk.get("kpmin"),
             perm=("uniform" if dk.get("parameterfromfile", "F").upper()
                   not in ("T", "TRUE", ".TRUE.") else dk.get("parameter_file")))
    r["tau0"] = (r["mu"] or 0) * (r["sigma"] or 0)
    r["dp_crit"] = r["sigma"] * (1 - r["mu"] / F0)
    r["can_slip"] = r["dp_crit"] <= DP_OBS
    r["bounds_ok"] = bounds_ok(dk)
    if base is None:
        return r
    r["code_fixed"] = code_is_fixed(n, base)
    d = sf.run_data(n, dk)
    if d is None:
        return r
    r["reached"] = float(d["t_end"])
    pk = peak_slip(n, base, dk)
    if pk is not None:
        r["peak_over_dc"] = pk / sf.ffloat(dk["dc"])
    m = d["T"] <= WINDOW_D
    if m.sum() >= 2:
        lam = sf.fit(d["T"][m], d["R"][m])
        if np.isfinite(lam):
            r["lam"] = float(lam)
            r["lam_ratio"] = float(lam / lam_obs)
    if d["tpw"] is not None and dk.get("injection_file") == "june_clean.txt":
        hi = min(d["tpw"][-1], obs["tp"].max(), WINDOW_D)
        if hi > 0.2:
            gr = np.linspace(0.05, hi, 2000)
            ps = np.interp(gr, d["tpw"], d["ppw"])
            ob = np.interp(gr, obs["tp"], obs["pm"])
            qg = np.interp(gr, obs["ti"], obs["q"])
            fl = (qg > 0.25 * np.nanmax(obs["q"])) & (ob > 5.0)
            if fl.sum() >= 30:
                r["p_pct"] = float(100 * (ps[fl].mean() / ob[fl].mean() - 1))
                r["p_rms"] = float(np.sqrt(np.mean((ps[fl] - ob[fl]) ** 2)))
                r["p_sim"] = float(ps[fl].mean())
                r["p_meas"] = float(ob[fl].mean())
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jobs", nargs="*", type=int)
    ap.add_argument("--stage", type=int)
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()

    jobs = list(a.jobs)
    if a.stage:
        jobs = STAGES[a.stage] + (EXTRA_STAGE2 if a.stage == 2 else [])
    if a.all:
        jobs = sorted(set(STAGES[1] + STAGES[2] + EXTRA_STAGE2))
    if not jobs:
        raise SystemExit("give job numbers, --stage N, or --all")

    obs = sf.observed()
    mo = obs["ft"] <= WINDOW_D
    lam_obs = sf.fit(obs["ft"][mo], obs["fd"][mo])
    print(f"window 0-{WINDOW_D:.0f} d   lambda_obs = {lam_obs:.4f}   "
          f"peak measured overpressure {DP_OBS:.2f} MPa\n")

    rows = [score(n, obs, lam_obs) for n in jobs if (IN / f"res{n}.in").exists()]

    print(f"{'run':>7s} {'sigma':>6s} {'mu':>6s} {'tau0':>6s} {'eta':>8s} "
          f"{'beta':>10s} {'permev':>7s} {'perm':>30s} {'reach':>6s} {'lam/obs':>8s} "
          f"{'p %':>8s} {'RMS':>6s} {'pk/dc':>7s} {'code':>6s} {'bnds':>5s}")
    print("-" * 152)
    for r in rows:
        lam = ("-" if "lam_ratio" not in r else
               f"{r['lam_ratio']:.2f}" if r.get("peak_over_dc", 9) >= 1
               else "no slip")
        if "lam_ratio" not in r and r.get("peak_over_dc", 9) < 1:
            lam = "no slip"
        pm = (f"{os.path.basename(str(r['perm']))[:34]}" if r["perm"] != "uniform"
              else "uniform")
        print(f"{r['n']:>7d} {r['sigma']:>6.2f} {r['mu']:>6.3f} {r['tau0']:>6.2f} "
              f"{r['eta']:>8.2e} {r['beta']:>10.3e} "
              f"{r['permev']:>7s} {pm[:30]:>30s} "
              f"{r.get('reached', float('nan')):>6.2f} {lam:>8s} "
              f"{r.get('p_pct', float('nan')):>+8.1f} "
              f"{r.get('p_rms', float('nan')):>6.1f} "
              f"{r.get('peak_over_dc', float('nan')):>7.0f} "
              f"{('fixed' if r.get('code_fixed') else 'BUGGY' if r.get('code_fixed') is False else '?'):>6s} "
              f"{('ok' if r.get('bounds_ok') else 'BAD' if r.get('bounds_ok') is False else '?'):>5s}")

    usable = [r for r in rows
              if r.get("code_fixed") and r.get("bounds_ok")
              and r.get("reached", 0) >= 0.95 * WINDOW_D]
    print(f"\n{len(usable)}/{len(rows)} usable (fixed code, bounds match the map, "
          f"reached >={0.95*WINDOW_D:.2f} d)")
    hits = [r for r in usable if "lam_ratio" in r and "p_pct" in r
            and abs(r["lam_ratio"] - 1) <= 0.15 and abs(r["p_pct"]) <= 15]
    if hits:
        print("\nINSIDE BOTH BANDS (+/-15%):")
        for r in sorted(hits, key=lambda z: abs(z["lam_ratio"] - 1)
                        + abs(z["p_pct"]) / 100):
            print(f"  {r['n']}  tau0 {r['tau0']:.2f} MPa  front {r['lam_ratio']:.2f}x  "
                  f"pressure {r['p_pct']:+.1f}%  permev {r['permev']}")
    else:
        print("\nnothing inside both bands yet")

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "grid_scores.json").write_text(json.dumps(rows, indent=1) + "\n")
    print(f"\nwrote {OUT}/grid_scores.json")


if __name__ == "__main__":
    main()
