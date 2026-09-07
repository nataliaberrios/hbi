#!/usr/bin/env python3
"""Stage 16: restart the clock where injection resumes. Day 0 = day 4.300.

THE HYPOTHESIS. 86 runs have failed to reproduce the Cooper Basin slip front,
wellhead pressure and slip magnitude simultaneously. Maybe no single model has
to: cycle 1 changes the fault irreversibly (kT = 1e15, so permeability
enhancement never heals), so cycle 1's end state is a legitimate free initial
condition rather than something a cold-start model must predict. Match the
post-shut-in period on its own terms.

t0 = 4.300 d, WHERE INJECTION RESUMES. The data offers four candidate "end of
the shut-in" times and they are not interchangeable:

    1.582 d   injection stops                              (injection record)
    1.590 d   wellhead gauge bled to ~0                    (pressure record)
    3.555 d   well closed in, gauge recovers to 32.50 MPa
    4.300 d   injection resumes                            <- t0

Seismicity also has a gap from 1.6 to 4.3 d, so nothing observable is lost.
Observed state at t0: front 313 m, cumulative slip 2.81 cm, wellhead 33.27 MPa
against a 34.41 MPa reservoir baseline -- about ZERO overpressure. The pressure
has no memory of cycle 1; only the permeability does. pfinit = 0d0 is therefore
already correct in the parent and is not changed.

ONE PARAMETER PER SIMULATION. An earlier draft of this stage changed three
things at once (injection record, disc radius, f0) and was discarded for exactly
that reason; res632940.in and res632941.in are its dead decks, kept on disk but
referenced by nothing. Here:

  632942  = res632913.in with ONLY injection_file and tmax changed
  632943  = res632942.in with ONLY parameter_file changed (disc 150 -> 300 m)

tmax is a duration, not physics, so changing it alongside the injection record
does not break one-at-a-time. Everything else -- f0 0.60, muinit 0.37,
sigmainit 27.99, a, b, dc, ds, imax, kpmax, kpmin, pfinit -- is inherited.

WHY res632913.in AS PARENT. It is the project's best full-record run: front 1.08
and slip 0.95x on 0-18 d, wellhead +58%. So 632942 answers directly whether the
model that best fits the whole record does better when asked only about the
post-shut-in period, with the restart as the sole cause.

WHY THE DISC RADIUS AS THE FIRST VARIATION. It is the one parameter specific to
this idea: it encodes how much enhancement cycle 1 left behind, and it is the
one thing a cold-start model cannot represent. 300 m is the observed front just
after the shut-in (313 m at 4.3 d). For scale, at 3.6 d res632913 has an
enhanced radius of 245 m and res632924 (f0 0.40) has 570 m, bracketing it.

SCORING IS NOT COMPARABLE TO THE OTHER 86 RUNS. Use score_restart.py, not
score_grid.py:
  * every observed series must be read at t_obs = t_sim + 4.300 d
  * slip must be compared against the observed INCREMENT, obs(t+4.3) - 2.81 cm,
    because main_LH.f90:923 sets slip = 0d0 unconditionally -- HBI cannot start
    with pre-existing slip. Target 6.36 cm by 17 d, i.e. 12.7 d of sim time.
  * lambda from R = lambda*sqrt(t) through the origin is INVALID: the observed
    front is already 313 m at t0 and the simulated front starts at the disc
    edge, so neither passes through the origin, and lambda = sum(sqrt(t)*R)/sum(t)
    weights points by t. Compare R(t) directly at matched times.

Usage:  python build_stage16.py [--write]
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import exp1
from scipy.optimize import brentq

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
PARENT = 632913
QFILE = "june_clean_from_d4300.txt"
T0_D = 4.300
DAYS = 13.1                 # observed record ends 13.155 d after t0
DISC_NEW = 300.0            # m, from the observed 313 m front at t0
# E1 front law anchored to res632913's own measured 890 m at 18 d
A_FRONT, ALPHA = 3.56, 1.248
SAFE = 0.8


def read_deck(p):
    out = []
    for line in Path(p).read_text().splitlines():
        if line.startswith("!") or not line.strip():
            continue
        w = line.split()
        if len(w) >= 2:
            out.append((w[0], " ".join(w[1:])))
    return out


def ff(x):
    return float(str(x).replace("d", "e").replace("D", "e"))


def front(t_d, dtauc):
    f = lambda r: A_FRONT * exp1(r ** 2 / (4 * ALPHA * t_d * 86400)) - dtauc
    try:
        return brentq(f, 1e-3, 3e4)
    except ValueError:
        return float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()

    assert (IN / QFILE).exists(), f"{QFILE} not built"
    pairs = read_deck(IN / f"res{PARENT}.in")
    base = dict(pairs)
    sig = ff(base["sigmainit"]); tau0 = ff(base["muinit"]) * sig
    f0 = ff(base["f0"]); dtauc = f0 * sig - tau0
    ds = ff(base["ds"]) * 1000.0; IMAX = int(base["imax"]); half = IMAX * ds / 2
    kx = ff(base["kpmax"]); km = ff(base["kpmin"])
    D = ff(base["kp"]) / (ff(base["eta"]) * ff(base["phi"]) * ff(base["beta"]))
    Ld = math.sqrt(4 * D * DAYS * 86400)
    src_map = base["parameter_file"].strip('"')
    old_disc = 150.0
    newmap = f"perm_2zone_{IMAX}_ds{int(ds)}_disc{int(DISC_NEW)}_kmax{kx:.1e}.txt"

    Rc = front(DAYS, dtauc)
    usable = SAFE * half - Ld
    print(f"parent res{PARENT}.in: ds {ds:.0f} m, imax {IMAX}, f0 {f0:.2f}, "
          f"tau_0 {tau0:.2f} MPa, pfinit {base.get('pfinit')}")
    print(f"  t0 {T0_D} d, injection {QFILE}, tmax {DAYS} d "
          f"(observed ends {17.455-T0_D:.3f} d after t0)\n")
    print(f"{'run':>7s}  {'changes from':>14s}  {'what changes':<34s} "
          f"{'front est':>10s} {'L/half':>7s}")
    rows = [
        (632942, f"res{PARENT}.in", "injection_file, tmax", Rc),
        (632943, "res632942.in", f"parameter_file (disc {old_disc:.0f}->"
                                 f"{DISC_NEW:.0f} m)", Rc + (DISC_NEW - old_disc)),
    ]
    bad = 0
    for n, frm, what, R in rows:
        ratio = (R + Ld) / half
        ok = ratio < SAFE
        bad += not ok
        print(f"{n:>7d}  {frm:>14s}  {what:<34s} {R:>9.0f}m {ratio:>7.2f}"
              + ("" if ok else "   <-- BREACH"))
    print(f"\n  far-field diffusion L({DAYS} d) = {Ld:.0f} m; usable front "
          f"<= {usable:.0f} m at ds {ds:.0f} m")
    print(f"  front estimate: A*E1(R^2/(4*alpha*t)) = Dtauc {dtauc:.2f} MPa with")
    print(f"  A = {A_FRONT} MPa, alpha = {ALPHA} m^2/s, calibrated to "
          f"res{PARENT}'s measured 890 m at 18 d")
    if bad:
        sys.exit(f"\n{bad} deck(s) would breach L/half = {SAFE}. Nothing written.")

    if not a.write:
        print("\ndry run -- nothing written. re-run with --write")
        return

    if not (IN / newmap).exists():
        c = (IMAX - 1) // 2
        rr = np.hypot(*(np.mgrid[0:IMAX, 0:IMAX] - c)) * ds
        kp = np.where(rr.ravel() <= DISC_NEW + 1e-9, kx, km)
        (IN / newmap).write_text("kp\n" + "\n".join(f"{v:.6e}" for v in kp) + "\n")
        print(f"\nwrote {newmap}  ({int((kp == kx).sum())} disc cells)")

    common = [
        f"! STAGE 16 -- RESTART WHERE INJECTION RESUMES. day 0 = day {T0_D}.",
        "!",
        "! 86 runs have failed to reproduce the slip front, wellhead pressure and",
        "! slip magnitude together. Maybe no single model has to: cycle 1 changes the",
        "! fault irreversibly (kT = 1e15, enhancement never heals), so cycle 1's end",
        "! state is a legitimate free initial condition rather than something a",
        "! cold-start model must predict.",
        "!",
        "! t0 = 4.300 d, where INJECTION RESUMES. The data gives four candidates:",
        "!   1.582 d  injection stops",
        "!   1.590 d  wellhead gauge bled to ~0",
        "!   3.555 d  well closed in, gauge recovers to 32.50 MPa",
        "!   4.300 d  injection resumes                                  <- t0",
        "! Seismicity has a gap from 1.6 to 4.3 d, so nothing observable is lost.",
        "!",
        "! Observed at t0: front 313 m, cumulative slip 2.81 cm, wellhead 33.27 MPa",
        "! against a 34.41 MPa reservoir baseline -- about ZERO overpressure. The",
        "! pressure has no memory of cycle 1; only the permeability does. pfinit = 0d0",
        "! is therefore already correct and is NOT changed.",
        "!",
        f"! INJECTION {QFILE}: june_clean.txt's times shifted by",
        f"! -{T0_D} d, keeping the 284 points at or after it plus an interpolated",
        "! t = 0 point (q = 8.4264e-04). june_clean.txt is UNTOUCHED (md5",
        "! 603285de6a5b1145d61bea05d97ebbf7). Volume per unit width is conserved to",
        "! 0.00000% against the record integrated from t0. The peak rate is LOWER",
        "! than the full record's, 7.99e-3 vs 8.74e-3, because the record's peak",
        "! falls BEFORE the shut-in -- a property of the experiment, not an error.",
        "! Injection is ON at t0, so the run starts on a ramp rather than from rest.",
        "!",
        "! SCORING: use score_restart.py, NOT score_grid.py. Observed series must be",
        f"! read at t_obs = t_sim + {T0_D} d; slip must be compared against the observed",
        "! INCREMENT (obs - 2.81 cm) because main_LH.f90:923 sets slip = 0d0 and HBI",
        "! cannot start with pre-existing slip; and lambda from R = lambda*sqrt(t)",
        "! through the origin is INVALID because neither front starts at the origin.",
        "!",
        "! ONE PARAMETER PER SIMULATION. An earlier draft changed three things at once",
        "! (injection record, disc radius, f0) and was discarded; res632940.in and",
        "! res632941.in are its dead decks, referenced by nothing.",
    ]
    hdrs = {
        632942: common + [
            "!",
            f"! THIS RUN: res{PARENT}.in with ONLY injection_file and tmax changed.",
            f"! res{PARENT}.in is the project's best full-record run -- front 1.08 and",
            "! slip 0.95x on 0-18 d, wellhead +58% -- so this asks directly whether the",
            "! model that best fits the whole record does better when asked only about",
            "! the post-shut-in period, with the restart as the sole cause.",
            "!",
            f"! Its initial disc stays at {old_disc:.0f} m ({src_map}), which is BELOW the",
            "! observed 313 m front at t0. If this run simply trails the observation by",
            "! ~160 m for its whole length, the disc radius is doing all the work and",
            "! res632943.in is the real test.",
            "!",
            f"! Domain: front estimate {Rc:.0f} m at {DAYS} d plus {Ld:.0f} m far-field",
            f"! diffusion against a {half:.0f} m half-domain, L/half {(Rc+Ld)/half:.2f}, "
            f"limit {usable:.0f} m.",
        ],
        632943: common + [
            "!",
            "! THIS RUN: res632942.in with ONLY parameter_file changed, disc",
            f"! {old_disc:.0f} -> {DISC_NEW:.0f} m. That is the one parameter specific to this idea:",
            "! how much permeability enhancement cycle 1 left behind, which a",
            f"! cold-start model cannot represent. {DISC_NEW:.0f} m is the observed front just",
            "! after the shut-in (313 m at 4.3 d). For scale, at 3.6 d res632913 has an",
            "! enhanced radius of 245 m and res632924 (f0 0.40) has 570 m.",
            "!",
            "! This is the original 'nonuniform initial permeability' hypothesis, but",
            "! with a physical origin -- the product of cycle 1 -- rather than an",
            "! arbitrary map.",
            "!",
            f"! Domain: front estimate {Rc+150:.0f} m (cold front plus the extra 150 m of",
            f"! pre-enhanced radius) plus {Ld:.0f} m diffusion, L/half "
            f"{(Rc+150+Ld)/half:.2f}, limit {usable:.0f} m.",
        ],
    }
    ovs = {
        632942: {"injection_file": f'"{QFILE}"',
                 "tmax": f"{DAYS/365.0:.8f}", "nstep": "400000"},
        632943: {"injection_file": f'"{QFILE}"',
                 "tmax": f"{DAYS/365.0:.8f}", "nstep": "400000",
                 "parameter_file": f'"{newmap}"'},
    }
    for n in (632942, 632943):
        lines = [f"filenumber {n}"] + hdrs[n]
        for k, v in pairs:
            if k == "filenumber":
                continue
            lines.append(f"{k} {ovs[n].get(k, v)}")
        (IN / f"res{n}.in").write_text("\n".join(lines) + "\n")

    print("\nverifying -- each deck differs from ITS OWN parent by one thing")
    nbad = 0
    for n, parent_file, allowed in (
            (632942, f"res{PARENT}.in",
             {"filenumber", "injection_file", "tmax", "nstep"}),
            (632943, "res632942.in", {"filenumber", "parameter_file"})):
        p = dict(read_deck(IN / parent_file)); d = dict(read_deck(IN / f"res{n}.in"))
        diff = {k for k in set(d) | set(p) if d.get(k) != p.get(k)}
        mp = np.loadtxt(IN / d["parameter_file"].strip('"'), skiprows=1)
        ck = {
            f"only {sorted(allowed)} vs {parent_file}": not (diff - allowed),
            "shifted injection": d["injection_file"].strip('"') == QFILE,
            "pfinit == 0": ff(d["pfinit"]) == 0.0,
            "tau_0 unchanged": (d["muinit"] == base["muinit"]
                                and d["sigmainit"] == base["sigmainit"]),
            "f0 unchanged": d["f0"] == base["f0"],
            "ds/imax unchanged": (d["ds"] == base["ds"]
                                  and d["imax"] == base["imax"]),
            "kpmax == map max": abs(mp.max() - kx) / kx < 1e-6,
            "kpmin == map min": abs(mp.min() - km) / km < 1e-6,
            "map size == imax^2": mp.size == IMAX ** 2,
            "tmax": abs(ff(d["tmax"]) * 365 - DAYS) < 1e-4,
        }
        ok = all(ck.values()); nbad += not ok
        print(f"  res{n}.in: " + ("OK" if ok else "PROBLEM "
              + str([k for k, v in ck.items() if not v])))
        print(f"      changed vs {parent_file}: {sorted(diff)}")
    print(f"\n  {2-nbad}/2 clean")
    sys.exit(0 if nbad == 0 else 1)


if __name__ == "__main__":
    main()
