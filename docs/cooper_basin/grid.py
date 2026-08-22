#!/usr/bin/env python3
"""The Cooper Basin exploration, stated as the grid search it actually is.

It grew stage by stage, which hid two things: how many cells the axes imply, and
which cells are already covered by the archive. This defines the coordinate
system once, enumerates it, and matches every cell against the existing decks and
output so coverage is a fact rather than a memory.

AXES, and why each one is an axis rather than a fixed value:

  sigmainit   30.0 is what res1807/res1808 use -- a round number. 27.99 is
              derived from the stress data and verified against Taiyi's
              setup_model.m: |snn| - p_pore, snn from rotate_stress(-160e6, 0,
              -100e6, 10 deg) = -101.8092 MPa, minus 73.82 -> 27.9892 MPa.
              It matters because dp_crit = sigmainit*(1 - muinit/f0) scales with
              it, so it decides whether the fault can fail at the measured
              overpressure at all.
  eta         8.9e-4 is Taiyi's, and his own comment says it is water at 25 C --
              wrong for a 4.1 km geothermal well. 1.27e-4 is the correction.
  beta        entangled with eta through D = k/(eta*phi*beta). "D-matched" means
              beta scaled by 8.9e-4/1.27e-4 = 7.007874 so D is unchanged. But eta
              and beta do NOT enter only through the product: the Peaceman well
              index T = 2*pi*k/eta (m_diffusion.f90:669) depends on eta alone, so
              the wellbore coupling gamma is not preserved. Hence beta is its own
              axis, not a slaved quantity.
  kpmax       2.5e-13 (res1807) and 5e-14 (res1808). With kpmin = kp = 1e-15
              these set the enhancement contrast, 250x and 50x.
  permev      enhancement on or off. Off means kpmax/kpmin are inert and only the
              initial field acts -- a real and separate physical model, not a
              degenerate case.
  perm        uniform initial kp, or a two-zone map whose near-well disc IS kpmax
              and background IS kpmin. Never anything else -- of 112 archived
              runs with enhancement and a nonuniform map, the 35 obeying that
              rule had 29 reach dc, the 77 that did not had 2.
  muinit      the strength axis. Floored by the data: slip needs
              dp_crit <= peak measured overpressure, so muinit >= f0*(1 -
              dp_obs/sigmainit).

A full factorial over these is ~800 runs, which is not the plan. The plan takes
SLICES, and this script names which.

Usage:  python grid.py                 # coverage report
        python grid.py --pending       # only cells with no run yet
"""
import argparse
import glob
import os
from pathlib import Path

import numpy as np

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
SCRATCH = "/scratch/users/nberrios/3dhbi/output"
OAK = "/oak/stanford/groups/edunham/nberrios/3doutput"
RUNS = "/scratch/users/nberrios/3dhbi/runs"

ETA_TAIYI, ETA_CORRECT = 8.9e-4, 1.27e-4
RATIO = ETA_TAIYI / ETA_CORRECT

# parent -> (kpmax, beta at Taiyi viscosity)
PARENTS = {1807: ("2.5e-13", 2.25e-8), 1808: ("5e-14", 2.005e-8)}
SIGMAS = ["30.0", "27.99"]
DP_OBS_5D = 10.92          # peak measured downhole overpressure, first 5 d, MPa
WINDOW_D = 5.0             # the comparison window; a run shorter than this is not coverage
F0 = 0.6


def ff(x):
    try:
        return float(str(x).replace("d", "e").replace("D", "e"))
    except (TypeError, ValueError):
        return None


def deck(n):
    d = {}
    p = IN / f"res{n}.in"
    if not p.exists():
        return None
    for line in p.read_text().splitlines():
        if line.startswith("!") or not line.strip():
            continue
        w = line.split()
        if len(w) >= 2:
            d[w[0]] = " ".join(w[1:]).strip('"')
    return d


def perm_of(dk):
    """('uniform', kp) or ('twozone', near, far)."""
    if dk.get("parameterfromfile", "F").upper() not in ("T", "TRUE", ".TRUE."):
        return ("uniform", ff(dk.get("kp")), ff(dk.get("kp")))
    pf = dk.get("parameter_file")
    if not pf or not (IN / pf).exists():
        return ("missing", None, None)
    k = np.loadtxt(IN / pf, skiprows=1)
    lo, hi = float(k.min()), float(k.max())
    return (("uniform", hi, lo) if hi / lo < 1.001 else ("twozone", hi, lo))


def state(n):
    """('done', days) / ('running', days) / ('none', 0)."""
    for base in (f"{SCRATCH}/{n}", f"{OAK}/{n}"):
        t = f"{base}/time{n}.dat"
        if os.path.exists(t) and os.path.getsize(t):
            try:
                return "done", float(np.atleast_2d(np.loadtxt(t))[-1, 1]) / 86400.0
            except Exception:
                return "done", 0.0
    g = glob.glob(f"{RUNS}/*/output/time{n}.dat")
    for t in g:
        if os.path.getsize(t):
            try:
                return "running", float(np.atleast_2d(np.loadtxt(t))[-1, 1]) / 86400.0
            except Exception:
                return "running", 0.0
    return "none", 0.0


def code_is_fixed(n):
    """The limitsigma ratchet detector: max(maxnorm)/maxnorm[0] must be 1."""
    for base in (f"{SCRATCH}/{n}", f"{OAK}/{n}"):
        p = f"{base}/monitor{n}.dat"
        if os.path.exists(p) and os.path.getsize(p):
            try:
                a = np.atleast_2d(np.loadtxt(p))
                if a.shape[1] >= 6:
                    # bool() is essential: a bare numpy comparison yields
                    # numpy.bool_, and `numpy.False_ is False` is False, so an
                    # `fx is False` test silently missed every buggy-code run and
                    # reported it as merely too short.
                    return bool(a[:, 5].max() / a[0, 5] <= 1.0001)
            except Exception:
                return None
    return None


def cell_of(n):
    """Grid coordinates of an existing run, or None if it is off-grid."""
    dk = deck(n)
    if dk is None or dk.get("problem") != "3dp":
        return None
    if dk.get("injection_file") != "june_clean.txt":
        return None                      # not comparable to the June data
    kx, km = dk.get("kpmax"), ff(dk.get("kpmin"))
    if kx is None or km is None or abs(km - 1e-15) > 1e-17:
        return None                      # not on the 1807/1808 kpmin
    parent = next((p for p, (k, _) in PARENTS.items() if k == kx), None)
    if parent is None:
        return None
    ps, near, far = perm_of(dk)
    if ps == "twozone" and not (abs(near - ff(kx)) / ff(kx) < 0.01
                                and abs(far - km) < 1e-17):
        return None                      # bounds do not match the map
    eta, beta = ff(dk.get("eta")), ff(dk.get("beta"))
    pev = "T" if dk.get("permev", "F").upper().startswith("T") else "F"
    return dict(n=n, parent=parent, sigma=dk.get("sigmainit"), eta=eta, beta=beta,
                mu=ff(dk.get("muinit")), perm=ps, permev=pev, dc=ff(dk.get("dc")),
                tmax_d=ff(dk.get("tmax")) * 365.0 if ff(dk.get("tmax")) else None)


def fluid_label(eta, beta, parent):
    b0 = PARENTS[parent][1]
    if abs(eta - ETA_TAIYI) / ETA_TAIYI < 0.02:
        return "A" if abs(beta - b0) / b0 < 0.02 else f"eta_T, beta x{beta/b0:.2f}"
    if abs(eta - ETA_CORRECT) / ETA_CORRECT < 0.02:
        if abs(beta - b0 * RATIO) / (b0 * RATIO) < 0.02:
            return "B (D-matched)"
        return f"eta_c, beta x{beta/b0:.2f}"
    return f"eta {eta:.2e}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pending", action="store_true")
    a = ap.parse_args()

    # every June-driven, kpmin-1e-15, consistent-bounds run in the archive
    found = []
    for p in sorted(IN.glob("res*.in")):
        st = p.stem[3:]
        if st.isdigit():
            c = cell_of(int(st))
            if c:
                found.append(c)

    print(f"{len(found)} archived runs sit on this grid "
          f"(June injection, kpmin 1e-15, bounds matching the map)\n")
    print(f"{'run':>7s} {'parent':>7s} {'sigma':>7s} {'fluid':>16s} {'perm':>9s} "
          f"{'permev':>7s} {'muinit':>7s} {'tau0':>6s} {'dp_crit':>8s} {'can slip?':>10s} "
          f"{'tmax d':>7s} {'state':>8s} {'days':>6s} {'code':>6s}")
    print("-" * 130)
    for c in sorted(found, key=lambda z: (z["parent"], z["sigma"], z["n"])):
        s, days = state(c["n"])
        fx = code_is_fixed(c["n"])
        sig = ff(c["sigma"])
        dpc = sig * (1 - c["mu"] / F0)
        print(f"{c['n']:>7d} {c['parent']:>7d} {c['sigma']:>7s} "
              f"{fluid_label(c['eta'], c['beta'], c['parent']):>16s} {c['perm']:>9s} "
              f"{c['permev']:>7s} {c['mu']:>7.4f} {c['mu']*sig:>6.2f} {dpc:>8.2f} "
              f"{('yes' if dpc <= DP_OBS_5D else 'NO'):>10s} "
              f"{c['tmax_d']:>7.2f} {s:>8s} {days:>6.2f} "
              f"{('fixed' if fx else 'BUGGY' if fx is False else '?'):>6s}")

    # ---- the slices the plan actually takes
    print("\n" + "=" * 130)
    print("PLANNED SLICES  (a full factorial over these axes would be ~800 runs)")
    print("=" * 130)
    def key(c):
        return (c["parent"], c["sigma"],
                fluid_label(c["eta"], c["beta"], c["parent"]),
                c["perm"], c["permev"], round(c["mu"], 4))

    # A cell counts as COVERED only if a run finished on FIXED code. A deck that
    # exists but has not run is "queued", which is not the same thing -- the
    # earlier version of this script conflated them and reported 8/8 for a slice
    # where nothing had actually run.
    have, queued = {}, {}
    for c in found:
        st, days = state(c["n"])
        fx = code_is_fixed(c["n"])
        # Reaching the window matters: res1703 stopped at 0.06 d and res532440 at
        # 0.46 d, so neither covers a 5 d cell even though both "finished".
        if st == "done" and fx and days >= 0.95 * WINDOW_D:
            have.setdefault(key(c), []).append((c["n"], days))
        else:
            why = ("BUGGY code" if fx is False else
                   "no monitor -> code unverified" if fx is None and st == "done"
                   else st if st != "done" else f"only {days:.2f} d")
            queued.setdefault(key(c), []).append((c["n"], why))
    slices = [
        ("S1 baseline, uniform perm, enhancement ON",
         [(p, s, f, "uniform", "T", 0.37) for p in PARENTS for s in SIGMAS
          for f in ("A", "B (D-matched)")]),
        ("S2 two-zone perm, consistent bounds, enhancement ON",
         [(p, s, f, "twozone", "T", 0.37) for p in PARENTS for s in SIGMAS
          for f in ("A", "B (D-matched)")]),
        ("S2b two-zone perm, enhancement OFF (the control for S2)",
         [(p, s, f, "twozone", "F", 0.37) for p in PARENTS for s in SIGMAS
          for f in ("A", "B (D-matched)")]),
    ]
    for name, cells in slices:
        cov = [c for c in cells if c in have]
        que = [c for c in cells if c not in have and c in queued]
        miss = [c for c in cells if c not in have and c not in queued]
        print(f"\n{name}: {len(cells)} cells | {len(cov)} DONE (fixed code, "
              f">={0.95*WINDOW_D:.1f} d), {len(que)} deck exists but not usable yet, "
              f"{len(miss)} no deck")
        for c in cov:
            rr = ", ".join(f"{n} ({d:.2f} d)" for n, d in have[c])
            print(f"   DONE     p{c[0]} sigma {c[1]:>6s} {c[2]:<16s} permev {c[4]}"
                  f"  -> {rr}")
        for c in que:
            rr = ", ".join(f"{n} [{w}]" for n, w in queued[c])
            print(f"   not yet  p{c[0]} sigma {c[1]:>6s} {c[2]:<16s} permev {c[4]}"
                  f"  -> {rr}")
        for c in miss:
            sig = ff(c[1])
            dpc = sig * (1 - c[5] / F0)
            print(f"   NO DECK  p{c[0]} sigma {c[1]:>6s} {c[2]:<16s} permev {c[4]}"
                  f"  (dp_crit {dpc:.2f} -> "
                  f"{'can slip' if dpc <= DP_OBS_5D else 'NO SLIP'})")

    print("\nmuinit axis floor, from the data alone (muinit >= f0*(1 - dp_obs/sigma)):")
    for s in SIGMAS:
        sig = ff(s)
        mm = F0 * (1 - DP_OBS_5D / sig)
        print(f"   sigma {s:>6s} -> muinit >= {mm:.4f}, tau0 >= {mm*sig:.2f} MPa")


if __name__ == "__main__":
    main()
