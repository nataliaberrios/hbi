#!/usr/bin/env python3
"""Stage 14: f0 = 0.40 -- the value the MEASURED activation pressure requires.

THE MEASUREMENT THIS TARGETS. The Habanero Field Development Plan (sec 4,
Fig 4-3) and Holl (2015, sec 8.4) both report the downhole overpressure at which
seismicity began on the Habanero Fault: dp = 2.6 MPa for the H04 local
stimulation (Oct 2012) and dp ~ 0.4 MPa for the extended stimulation (Nov 2012),
which is the record this project models. That is dp_crit measured directly -- no
stress tensor, no assumed friction, no fault dip.

    measured, H04 2012                     0.4-2.6 MPa
    Wang & Dunham (tau_0 15.0, f0 0.60)          3.0
    res632913     (tau_0 10.36, f0 0.60)        10.72
    res632913's measured plateau                 11.9

The model needs 4-27x the overpressure at which the fault actually failed. That
is its +58% wellhead restated as a measurement, and it explains all 82 runs at
once: anything that reproduced the observed slip had to over-pressurise, because
the threshold was set an order of magnitude too high. Requiring the measured
threshold at this project's tau_0 and sigmabar_0 gives f0 = 0.375 (dp 0.4) to
0.408 (dp 2.6). f0 = 0.40 is the middle.

This is NOT the retracted f0 = 0.433 argument. That was retracted because 0.433
needs ~82 wt.% chlorite against a documented maximum of 35 -- correct on its own
terms, but it answered "what gouge composition gives f0 = 0.433?" The field data
asks "what f0 is consistent with a fault that slipped at 0.4 MPa?", which does
not go through friction databases at all.

The stress state itself is sound and unchanged here. Holl & Barton's ratios
(SHmax/Shmin/Sv 1.35-1.45/1.10-1.25/1.0, Pp 72.7 MPa) with Sv ~ 95.3 MPa -- from
Holl sec 8.4, first seismicity 2.6 MPa above Pp and ~20 MPa below the overburden
-- give tau 9.8-12.6 MPa and sigmabar_n 25.8-26.7 MPa on an 18-20 deg thrust.
tau_0 = 10.36 / sigmabar_0 = 27.99 sits inside that band.

WHY ds CHANGES, WHICH IS THE AWKWARD PART. A weaker fault slips further. Fronts
projected from the E1 law R: A*E1(R^2/(4*alpha*t)) = Dtauc, with
alpha = kpmax/(eta*phi*beta) = 1.248 m^2/s from physics and A = 3.56 MPa
calibrated to res632913's own measured 890 m at 18 d (the same law over-predicts
its 5 d front by 1.23x, so this sizing is conservative):

      f0   Dtauc    R(3d)   R(5d)   R(9d)   R(18d)
    0.60   6.43M     363m    469m    629m     890m
    0.40   0.84M    1115m   1439m   1931m    2731m

Against a 1502 m half-domain at ds 5 m (usable limit 1026 m after allowing 176 m
of far-field diffusion), f0 = 0.40 breaches at 2.5 d. ds 5 m simply cannot test
it. Enlarging imax instead is not available: imax 701 and 901 both SEGFAULT
inside HACApK's HACApK_generate_frame_LH, right after the MPI communicator
split, at --mem=64G and 220G alike and with ulimit -s unlimited -- the H-matrix
is only 17.5 GB, so it is not memory. That is vendored code and is not to be
modified.

     ds    half    usable   disc cells   T vs ds5   safe tmax, f0 0.40
     5m   1502m    1026m         30.0      1.00x          2.5 d
    10m   3005m    2228m         15.0      0.78x         11.8 d
    20m   6010m    4632m          7.5      0.64x         18.0 d

Both ds are run, for the trade they represent: 10 m keeps the 150 m disc at 15
cells and the Peaceman index within 22% of the parent, but stops at 11 d; 20 m
reaches the full 18 d with the disc at 7.5 cells and T 36% off.

EACH ds GETS ITS OWN f0 = 0.60 CONTROL, and the sweep is uninterpretable
without them. ds enters the Peaceman well index
T = 2*pi*k/eta/(log(0.2*ds/rw)+skin) directly, and T sets the wellhead -- one of
the quantities under test. A bare f0 = 0.40 run at a new ds would confound the
friction change with a mesh change. The controls isolate it.

  632924  ds 10 m  f0 0.40  11 d
  632925  ds 10 m  f0 0.60  11 d   <- control
  632926  ds 20 m  f0 0.40  18 d
  632927  ds 20 m  f0 0.60  18 d   <- control

Usage:  python build_stage14.py [--write]
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
DISC_R = 150.0
IMAX = 601
# (new, ds_m, f0, days)
JOBS = [(632924, 10.0, 0.4000, 11.0),
        (632925, 10.0, 0.6000, 11.0),
        (632926, 20.0, 0.4000, 18.0),
        (632927, 20.0, 0.6000, 18.0)]
ALPHA = 1.248           # kpmax/(eta*phi*beta) for this configuration
A_FRONT = 3.56          # MPa, calibrated to 632913's measured 890 m at 18 d


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
        return brentq(f, 1e-3, 5e4)
    except ValueError:
        return float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()

    pairs = read_deck(IN / f"res{PARENT}.in")
    base = dict(pairs)
    sig = ff(base["sigmainit"]); tau0 = ff(base["muinit"]) * sig
    kx = ff(base["kpmax"]); km = ff(base["kpmin"])
    D_far = ff(base["kp"]) / (ff(base["eta"]) * ff(base["phi"]) * ff(base["beta"]))
    rw = ff(base["rw"]); den5 = math.log(0.2 * 5.0 / rw)

    print(f"parent res{PARENT}.in  tau_0 {tau0:.2f} MPa (MEASURED, unchanged), "
          f"sigmabar_0 {sig:.2f}, f0 {ff(base['f0']):.2f}, ds "
          f"{ff(base['ds'])*1000:.0f} m\n")
    print(f"{'new':>7s} {'ds':>5s} {'f0':>7s} {'Dtauc':>7s} {'dp_crit':>8s} "
          f"{'days':>5s} {'R(end)':>8s} {'+diff':>7s} {'half':>7s} {'L/half':>7s} "
          f"{'T/T5':>6s} {'role':>8s}")
    rows, bad = [], 0
    for new, ds, f0, days in JOBS:
        dt = f0 * sig - tau0
        dpc = sig - tau0 / f0
        half = IMAX * ds / 2
        Ld = math.sqrt(4 * D_far * days * 86400)
        R = front(days, dt)
        ratio = (R + Ld) / half
        den = math.log(0.2 * ds / rw)
        ok = ratio < 0.8
        bad += not ok
        rows.append((new, ds, f0, days, dt, dpc, half, R, Ld, ratio))
        print(f"{new:>7d} {ds:>4.0f}m {f0:>7.4f} {dt:>6.2f}M {dpc:>7.2f}M "
              f"{days:>4.0f}d {R:>7.0f}m {Ld:>6.0f}m {half:>6.0f}m {ratio:>7.2f} "
              f"{den5/den:>5.2f}x "
              + f"{'CONTROL' if f0 > 0.5 else 'test':>8s}"
              + ("" if ok else "   <-- BREACH"))
    if bad:
        sys.exit(f"\n{bad} deck(s) would breach L/half = 0.8. Nothing written.")
    print(f"\nall {len(rows)} within L/half < 0.8 "
          f"(front law over-predicts 632913's 5 d front by 1.23x, so conservative)")
    print(f"dp_crit target from the field measurement: 0.4-2.6 MPa")

    if not a.write:
        print("\ndry run -- nothing written. re-run with --write")
        return

    made = []
    for new, ds, f0, days, dt, dpc, half, R, Ld, ratio in rows:
        mapname = f"perm_2zone_{IMAX}_ds{int(ds)}_kmax{kx:.1e}.txt"
        if not (IN / mapname).exists():
            c = (IMAX - 1) // 2
            rr = np.hypot(*(np.mgrid[0:IMAX, 0:IMAX] - c)) * ds
            kp = np.where(rr.ravel() <= DISC_R + 1e-9, kx, km)
            (IN / mapname).write_text("kp\n"
                                      + "\n".join(f"{v:.6e}" for v in kp) + "\n")
            print(f"  wrote {mapname}  ({int((kp==kx).sum())} disc cells)")
        ndisc = int((np.loadtxt(IN / mapname, skiprows=1) == kx).sum())
        role = "CONTROL for the ds change" if f0 > 0.5 else \
               "TEST: f0 from the measured activation pressure"
        hdr = [
            f"! STAGE 14 -- {role}.",
            f"! Built from res{PARENT}.in; f0, ds, the map and tmax change.",
            "!",
            "! THE MEASUREMENT. Habanero FDP sec 4 (Fig 4-3) and Holl (2015) sec 8.4",
            "! both report the overpressure at which seismicity began on the Habanero",
            "! Fault: dp = 2.6 MPa (H04 local, Oct 2012) and dp ~ 0.4 MPa (extended,",
            "! Nov 2012 -- the record modelled here). That is dp_crit, MEASURED.",
            "!",
            "!     measured, H04 2012                    0.4-2.6 MPa",
            "!     Wang & Dunham (tau_0 15.0, f0 0.60)         3.0",
            f"!     res{PARENT} (tau_0 {tau0:.2f}, f0 0.60)           10.72",
            f"!     res{PARENT}'s measured plateau                11.9",
            "!",
            "! The model needs 4-27x the overpressure at which the fault actually",
            "! failed -- its +58% wellhead restated as a measurement. Requiring the",
            "! measured threshold gives f0 = 0.375 (dp 0.4) to 0.408 (dp 2.6).",
            "!",
            f"! f0 {ff(base['f0']):.2f} -> {f0:.4f} at the MEASURED tau_0 = {tau0:.2f} MPa, so",
            f"! dp_crit = sigmabar_0 - tau_0/f0 = {dpc:.2f} MPa and",
            f"! Dtauc = f0*sigmabar_0 - tau_0 = {dt:.2f} MPa.",
            "!",
            "! NOT the retracted f0 = 0.433 argument. That needed ~82 wt.% chlorite",
            "! against a documented max of 35 -- correct, but it asked what gouge",
            "! composition gives 0.433. The field data asks what f0 is consistent",
            "! with a fault that slipped at 0.4 MPa, which needs no friction table.",
            "!",
            "! The stress state is unchanged and is SOUND: Holl & Barton's ratios with",
            "! Sv ~ 95.3 MPa (Holl sec 8.4) give tau 9.8-12.6 and sigmabar_n 25.8-26.7",
            "! on an 18-20 deg thrust, and tau_0 10.36 / sigmabar_0 27.99 sits inside.",
            "!",
            f"! ds {ff(base['ds'])*1000:.0f} -> {ds:.0f} m. NOT a free choice: at ds 5 m the",
            f"! f0 = 0.40 front breaches the 1502 m half-domain at 2.5 d, and enlarging",
            "! imax is unavailable because imax 701 and 901 both segfault inside",
            "! HACApK's HACApK_generate_frame_LH after the MPI split, at 64G and 220G",
            "! alike and with ulimit -s unlimited (the H-matrix is only 17.5 GB, so it",
            "! is not memory). HACApK is vendored and not to be modified.",
            "!",
            f"! Projected front at {days:.0f} d is {R:.0f} m from",
            f"! A*E1(R^2/(4*alpha*t)) = Dtauc with alpha = kpmax/(eta*phi*beta) =",
            f"! {ALPHA:.3f} m^2/s and A = {A_FRONT:.2f} MPa calibrated to res{PARENT}'s own",
            f"! measured 890 m at 18 d. Plus {Ld:.0f} m of far-field diffusion, that is",
            f"! {R+Ld:.0f} m against a {half:.0f} m half-domain, L/half = {ratio:.2f}. The same",
            f"! law over-predicts res{PARENT}'s 5 d front by 1.23x, so this is",
            "! conservative.",
            "!",
            f"! {mapname}: the same {DISC_R:.0f} m disc at kpmax {kx:.1e} over background",
            f"! kpmin {km:.0e}, {ndisc} cells at ds {ds:.0f} m.",
            "!",
            f"! PAIRED with the f0 = {'0.60' if f0 < 0.5 else '0.40'} run at this same ds. ds enters the",
            "! Peaceman index T = 2*pi*k/eta/(log(0.2*ds/rw)+skin) directly, and T sets",
            f"! the wellhead -- here T is {den5/math.log(0.2*ds/rw):.2f}x its ds = 5 m value. Without the",
            "! pair the friction change and the mesh change cannot be separated.",
        ]
        ov = {"f0": f"{f0:.4f}", "ds": f"{ds/1000.0:.6f}",
              "parameter_file": f'"{mapname}"',
              "tmax": f"{days/365.0:.8f}", "nstep": "400000"}
        lines = [f"filenumber {new}"] + hdr
        for k, v in pairs:
            if k == "filenumber":
                continue
            lines.append(f"{k} {ov.get(k, v)}")
        (IN / f"res{new}.in").write_text("\n".join(lines) + "\n")
        made.append((new, ds, f0, days, mapname))

    print("\nverifying")
    nbad = 0
    for new, ds, f0, days, mapname in made:
        d = dict(read_deck(IN / f"res{new}.in"))
        mp = np.loadtxt(IN / mapname, skiprows=1)
        allowed = {"filenumber", "f0", "ds", "parameter_file", "tmax", "nstep"}
        diff = {k for k in set(d) | set(base) if d.get(k) != base.get(k)}
        ck = {
            "only allowed keys": not (diff - allowed),
            "f0": abs(ff(d["f0"]) - f0) < 1e-9,
            "ds": abs(ff(d["ds"]) * 1000 - ds) < 1e-6,
            "tau_0 UNCHANGED": (d["muinit"] == base["muinit"]
                                and d["sigmainit"] == base["sigmainit"]),
            "a,b,dc unchanged": all(d[k] == base[k] for k in ("a", "b", "dc")),
            "imax unchanged": d["imax"] == base["imax"] == str(IMAX),
            "map size": mp.size == IMAX ** 2,
            "kpmax == map max": abs(mp.max() - kx) / kx < 1e-6,
            "kpmin == map min": abs(mp.min() - km) / km < 1e-6,
            "injection unchanged": (d.get("injection_file")
                                    == base.get("injection_file")),
            "tmax": abs(ff(d["tmax"]) * 365 - days) < 1e-4,
        }
        ok = all(ck.values()); nbad += not ok
        print(f"  res{new}.in ds {ds:.0f} m f0 {f0:.2f}: "
              + ("OK" if ok else "PROBLEM "
                 + str([k for k, v in ck.items() if not v])))
    print(f"\n  {len(made)-nbad}/{len(made)} clean")
    sys.exit(0 if nbad == 0 else 1)


if __name__ == "__main__":
    main()
