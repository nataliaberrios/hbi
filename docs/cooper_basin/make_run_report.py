#!/usr/bin/env python3
"""One folder per simulation: its figures plus a params.txt describing the run.

    figures/<jobid>/params.txt      every deck value, plus derived quantities
    figures/<jobid>/*.png|pdf       figures belonging to that run alone
    figures/comparisons/            multi-run overlays, which belong to no single run

params.txt records more than the deck, because the deck alone does not tell you
what the run actually did:

  - the PERMEABILITY MAP is a separate file, so a deck reading "kp 1e-15" can be
    running a two-zone map with 1.1e-12 near the well. The near/far values and the
    disc radius are read out of the map itself.
  - hydraulic diffusivity D = kp/(eta phi beta), the quantity that actually sets
    how fast pressure spreads, appears nowhere in the deck.
  - the failure threshold: with tau0 = muinit*sigmainit fixed, slip needs
    dp > sigmainit - tau0/f0 of overpressure. That number decides whether the run
    slips at all and is likewise not written down anywhere.
  - tmax in DAYS. HBI's year is 365 days (main_LH.f90:211), not 365.25, so
    converting with 365.25 disagrees with the deck by ~17 min over 17 d.
  - how the run actually ended, from convergence_log.csv. A run that exhausted
    nstep is logged CONVERGED by the older harness, indistinguishable from one
    that reached tmax -- 83 of 343 CONVERGED rows had in fact truncated.

Usage:
    python make_run_report.py 632522 632526 ...
    python make_run_report.py --all          # every res6325xx with output
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path("/home/users/nberrios/3dhbi/hbi_analysis")
FIGROOT = HERE / "figures"
INPUT_DIR = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
SCRATCH = Path("/scratch/users/nberrios/3dhbi/output")
RUNS = Path("/scratch/users/nberrios/3dhbi/runs")
OAK = Path("/oak/stanford/groups/edunham/nberrios/3doutput")
# the 632510 post-fix rerun reuses filenumber 632510 in its own directory
EXTRA = {632510: Path("/scratch/users/nberrios/3dhbi/fix_632510/output")}
CONV = INPUT_DIR / "convergence_log.csv"
CONV2 = INPUT_DIR / "convergence_log_v2.csv"

DAYS_PER_YEAR = 365.0        # HBI's convention, main_LH.f90:211


def ffloat(x):
    try:
        return float(str(x).replace("d", "e").replace("D", "e"))
    except ValueError:
        return None


def read_deck(num):
    p = INPUT_DIR / f"res{num}.in"
    if not p.exists():
        return None, None
    order, d = [], {}
    for line in p.read_text().splitlines():
        if line.startswith("!") or not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 2:
            d[parts[0]] = parts[1].strip('"')
            order.append(parts[0])
    return d, order


def perm_map_summary(deck):
    """The permeability actually in force, which the deck does not state."""
    if deck.get("parameterfromfile", "F").upper() not in ("T", "TRUE", ".TRUE."):
        kp = ffloat(deck.get("kp"))
        return f"uniform {kp:.3e} m^2 (parameterfromfile F)" if kp else "unknown"
    pf = deck.get("parameter_file")
    if not pf:
        return "parameterfromfile T but no parameter_file given"
    p = INPUT_DIR / pf
    if not p.exists():
        return f"{pf} (MISSING)"
    k = np.loadtxt(p, skiprows=1)
    lo, hi = k.min(), k.max()
    # NOT np.isclose: its default atol=1e-8 dwarfs permeabilities of order 1e-15,
    # so every map compared equal and two-zone maps were reported as uniform.
    if lo > 0 and hi / lo < 1.001:
        return f"{pf}: uniform {lo:.3e} m^2"
    n_hi = int((k > (lo + hi) / 2).sum())
    ds_m = ffloat(deck.get("ds", 0)) * 1000.0
    # area of the high-perm patch -> equivalent disc radius
    r_eq = np.sqrt(n_hi * ds_m**2 / np.pi) if ds_m else float("nan")
    return (f"{pf}\n    near-well {hi:.3e} m^2 over {n_hi} cells "
            f"(equivalent disc radius {r_eq:.0f} m)\n"
            f"    far field {lo:.3e} m^2 over {k.size - n_hi} cells\n"
            f"    contrast  {hi/lo:.0f}x")


def conv_rows(num):
    out = []
    for f, label in ((CONV, "convergence_log.csv"), (CONV2, "convergence_log_v2.csv")):
        if not f.exists():
            continue
        head = None
        for ln in f.read_text().splitlines():
            if head is None:
                head = ln.split(","); continue
            c = ln.split(",")
            if len(c) > 3 and c[3] == str(num):
                out.append((label, dict(zip(head, c))))
    return out


def run_state(num):
    """Where the output is and how far it got."""
    for base, where in ((SCRATCH / str(num), "final"),
                        (OAK / str(num), "final (OAK)")):
        t = base / f"time{num}.dat"
        if t.exists():
            a = np.loadtxt(t)
            return where, str(base), len(a), a[-1, 1] / 86400.0
    # still running: output sits in the per-job rundir
    for d in sorted(RUNS.glob("*/output")):
        t = d / f"time{num}.dat"
        if t.exists() and os.path.getsize(t):
            a = np.atleast_2d(np.loadtxt(t))
            return "in progress", str(d), len(a), a[-1, 1] / 86400.0
    return "no output", "-", 0, 0.0


def extra_state(num):
    """A second run sharing this filenumber, e.g. the 632510 post-fix rerun."""
    d = EXTRA.get(num)
    if not d:
        return None
    t = d / f"time{num}.dat"
    if not t.exists():
        return None
    a = np.loadtxt(t)
    return str(d), len(a), a[-1, 1] / 86400.0


def write_params(num):
    deck, order = read_deck(num)
    if deck is None:
        print(f"  {num}: no deck, skipping"); return None
    folder = FIGROOT / str(num)
    folder.mkdir(parents=True, exist_ok=True)

    where, path, nframes, t_end = run_state(num)
    tmax_yr = ffloat(deck.get("tmax"))
    mu, sig, f0 = (ffloat(deck.get(k)) for k in ("muinit", "sigmainit", "f0"))
    phi, beta, eta = (ffloat(deck.get(k)) for k in ("phi", "beta", "eta"))
    imax, jmax = deck.get("imax"), deck.get("jmax")
    ds_km = ffloat(deck.get("ds"))

    L = [f"run {num}", "=" * 60, ""]
    L += ["STATE", f"  output       {where}: {path}",
          f"  frames       {nframes}",
          f"  reached      {t_end:.4f} d"]
    if tmax_yr:
        L += [f"  tmax         {tmax_yr:g} yr = {tmax_yr*DAYS_PER_YEAR:.4f} d "
              f"(HBI year = 365 d)",
              f"  completed    {'YES' if t_end >= tmax_yr*DAYS_PER_YEAR - 1e-3 else 'NO - stopped early'}"]
    for label, row in conv_rows(num):
        L += [f"  {label}: status={row.get('status','?')}"
              + (f", end_days={row['end_days']}" if 'end_days' in row else "")
              + (f", binary={row['exe_sha256']} ({row.get('git_branch','?')})"
                 if 'exe_sha256' in row else "")]
    ex = extra_state(num)
    if ex:
        L += ["", f"  SECOND RUN with this filenumber (post-limitsigma-fix rerun):",
              f"    output     {ex[0]}", f"    frames     {ex[1]}",
              f"    reached    {ex[2]:.4f} d"]
    L += [""]

    L += ["GEOMETRY",
          f"  grid         {imax} x {jmax} at ds {ds_km:g} km "
          f"({ffloat(imax)*ds_km:.3f} x {ffloat(jmax)*ds_km:.3f} km)" if ds_km else "",
          f"  problem      {deck.get('problem','?')}", ""]

    L += ["MECHANICS"]
    if None not in (mu, sig, f0):
        tau0 = mu * sig
        L += [f"  sigmainit    {sig:g} MPa",
              f"  muinit       {mu:g}   -> tau0 = {tau0:.3f} MPa",
              f"  f0           {f0:g}",
              f"  dp to fail   {sig - tau0/f0:.3f} MPa   (= sigmainit - tau0/f0)"]
    for k in ("a", "b", "dc", "limitsigma", "minsig", "velinit", "velmin"):
        if k in deck:
            L += [f"  {k:<12s} {deck[k]}"]
    L += [""]

    L += ["HYDRAULICS",
          f"  permeability {perm_map_summary(deck)}",
          f"  permev       {deck.get('permev','?')}   (does perm evolve with slip)"]
    for k in ("kpmin", "kpmax", "kL", "kT", "phi", "beta", "eta", "pfinit"):
        if k in deck:
            L += [f"  {k:<12s} {deck[k]}"]
    if None not in (phi, beta, eta):
        pfl = INPUT_DIR / deck["parameter_file"] if deck.get("parameter_file") else None
        if pfl and pfl.exists():
            k = np.loadtxt(pfl, skiprows=1)
            for nm, kv in (("near-well", k.max()), ("far field", k.min())):
                L += [f"  D ({nm})   {kv/(eta*phi*beta):.4g} m^2/s"]
        else:
            kp = ffloat(deck.get("kp"))
            if kp:
                L += [f"  D            {kp/(eta*phi*beta):.4g} m^2/s"]
    L += [f"  injection    {deck.get('injection_file','?')}", ""]

    L += ["FULL DECK (verbatim)", "-" * 60]
    L += [f"  {k:<22s} {deck[k]}" for k in order]

    (folder / "params.txt").write_text("\n".join(x for x in L if x is not None) + "\n")
    print(f"  {num}: wrote {folder}/params.txt  ({where}, {t_end:.2f} d)")
    return folder


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jobs", nargs="*", type=int)
    ap.add_argument("--all", action="store_true")
    a = ap.parse_args()
    jobs = a.jobs
    if a.all:
        jobs = sorted({int(p.stem.replace("res", ""))
                       for p in INPUT_DIR.glob("res6325*.in")
                       if p.stem.replace("res", "").isdigit()})
    if not jobs:
        sys.exit("give job numbers or --all")
    for n in jobs:
        write_params(n)


if __name__ == "__main__":
    main()
