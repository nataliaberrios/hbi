#!/usr/bin/env python3
"""Stage 13: the f0 sweep on 632913, each run in a domain sized for its front.

WHY f0. Across 82 runs, six reproduce the observed slip to within 3% and all six
sit at +70 to +80% wellhead; nothing matches all three targets at any tolerance
down to +/-50%. The model needs ~1.75x the measured pressure to slip. The only
quantity that sets how much pressure slip requires is the strength margin
Dtauc = f0*sigmabar_0 - tau_0, and tau_0 is measured, so it is f0 or nothing.

Wang & Dunham's own front criterion (source_code/cmp_seis_extent.m) makes this
quantitative and was re-derived from their Table 1 here:

    route                       tau_0     f0   Dtauc     1d     3d     5d
    OBSERVED                                           187m   323m   417m
    Wang & Dunham               15.00  0.600   1.80M   196m   340m   438m
    resolved stress, lab f0     10.26  0.600   6.54M     3m     5m     7m
    resolved stress, f0=0.433   10.26  0.433   1.86M   185m   321m   414m

The front pins Dtauc = 1.86 MPa, not tau_0. Their calibration is one of two
points on that locus; the other keeps the stress measurement and needs f0 near
0.43 -- inside the published range for chlorite-bearing granitic gouge (pure
chlorite 0.37, unaltered granite 0.60-0.71), and the value their own off-fault
spring-sliders run at. They already argue for phyllosilicates in this fault zone
to justify velocity-strengthening.

WHAT THIS TESTS THAT THE FORMULA CANNOT. Their criterion is a Coulomb threshold
on secondary faults. Whether HBI's rate-and-state main fault, with permeability
enhancement and a real injection record, behaves the same way is unchecked --
CONCLUSIONS.md flags it explicitly. Predicted, from Stage 10's finding that the
plateau tracks dp_crit: lowering f0 should drop the pressure while leaving the
slip roughly intact, i.e. 632913's front (1.08) and slip (0.95x) survive while
its +58% wellhead falls.

  632920  f0 0.5500  Dtauc 5.03 MPa  dp_crit 9.15 MPa
  632921  f0 0.5000  Dtauc 3.64 MPa  dp_crit 7.27 MPa
  632922  f0 0.4366  Dtauc 1.86 MPa  dp_crit 4.26 MPa   <- the addendum's value
                                                           for tau_0 = 10.36

DOMAIN, WHICH DRIVES THE ONLY REAL DESIGN CHOICE HERE. A weaker fault slips
further, and at imax 601 (half-domain 1502 m) ALL THREE would breach. Sizing
from the empirical front-vs-Dtauc scaling across runs sharing 632913's map:

     f0      Dtauc   projected front at 18 d   needs half-domain
   0.5500   5.03M            1124 m                 1625 m
   0.5000   3.64M            1532 m                 2135 m
   0.4366   1.86M            2900 m                 3845 m

CAVEAT ON THAT SCALING, stated because it drives the cost: the five runs sharing
632913's map span only two Dtauc values (6.44 and 6.90 MPa) and their fronts
scatter 0.59-0.70 AT FIXED Dtauc -- as much as the difference between the two
groups. The fitted exponent (front ~ Dtauc^-0.95) is therefore barely
constrained, and these projections could be wrong by a factor of two either way.
They are used only to size domains generously, never to predict a result.

ds STAYS 5 m in every run. Enlarging via imax rather than ds keeps the cell size,
the near-well resolution and the Peaceman well index
T = 2*pi*k/eta/(log(0.2*ds/rw)+skin) identical to 632913 -- and T sets the
wellhead, which is the quantity under test. Coarsening ds would change the
answer while appearing to change only the mesh.

632922 IS NOT LAUNCHED BY THIS SCRIPT. Its 3845 m half-domain means imax 1601,
2.56M cells against 632913's 361k -- roughly 7x, so a 10-15x runtime on an 11 h
baseline, past the 2 d partition limit. It is written so the deck exists and can
be sized properly once 632920/632921 measure the real sensitivity, and printed
with --write so the decision is visible rather than buried.

Usage:  python build_stage13.py [--write]
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
PARENT = 632913
DAYS = 18.0
DISC_R = 150.0          # m, from perm_2zone_601_ds5_kmax2.5e-13.txt
# (new, f0, imax, launch now?)
JOBS = [(632920, 0.5500, 701, True),
        (632921, 0.5000, 901, True),
        (632922, 0.4366, 1601, False)]


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()

    pairs = read_deck(IN / f"res{PARENT}.in")
    base = dict(pairs)
    sig, tau0 = ff(base["sigmainit"]), ff(base["muinit"]) * ff(base["sigmainit"])
    ds_m = ff(base["ds"]) * 1000.0
    kx = ff(base["kpmax"]); km = ff(base["kpmin"])
    D = ff(base["kp"]) / (ff(base["eta"]) * ff(base["phi"]) * ff(base["beta"]))
    L = math.sqrt(4 * D * DAYS * 86400)
    # empirical front-vs-Dtauc exponent, from runs sharing this map
    EXP, F_REF, DT_REF = -0.95, 890.0, 0.60 * sig - tau0

    print(f"parent res{PARENT}.in  tau_0 {tau0:.2f} MPa, sigmabar_0 {sig:.2f}, "
          f"f0 {ff(base['f0']):.2f}, ds {ds_m:.0f} m")
    print(f"  its 18 d front is {F_REF:.0f} m at Dtauc {DT_REF:.2f} MPa; "
          f"far-field diffusion L(18 d) = {L:.0f} m\n")
    print(f"{'new':>7s} {'f0':>7s} {'Dtauc':>7s} {'dp_crit':>8s} {'imax':>5s} "
          f"{'half':>6s} {'front?':>7s} {'L/half':>7s} {'cells':>9s} {'run':>5s}")
    rows = []
    for new, f0, imax, go in JOBS:
        dt = f0 * sig - tau0
        dpc = sig - tau0 / f0
        half = imax * ds_m / 2
        fr = F_REF * (dt / DT_REF) ** EXP
        ratio = (fr + L) / half
        rows.append((new, f0, imax, dt, dpc, half, fr, ratio, go))
        print(f"{new:>7d} {f0:>7.4f} {dt:>6.2f}M {dpc:>7.2f}M {imax:>5d} "
              f"{half:>5.0f}m {fr:>6.0f}m {ratio:>7.2f} "
              f"{imax**2:>9,d} {'yes' if go else 'HOLD':>5s}")
        if ratio >= 0.8:
            print(f"          ^ L/half {ratio:.2f} exceeds 0.8 -- domain too small")
    bad = [r for r in rows if r[7] >= 0.8]
    if bad:
        sys.exit(f"\n{len(bad)} deck(s) exceed L/half = 0.8. Not writing.")
    print(f"\nall {len(rows)} projected within L/half < 0.8")

    if not a.write:
        print("\ndry run -- nothing written. re-run with --write")
        return

    made = []
    for new, f0, imax, dt, dpc, half, fr, ratio, go in rows:
        # map for this imax: same 150 m disc at kpmax, background kpmin
        mapname = f"perm_2zone_{imax}_ds{ds_m:.0f}_kmax{kx:.1e}.txt".replace(
            "ds5.", "ds5.")
        mapname = f"perm_2zone_{imax}_ds{int(ds_m)}_kmax{kx:.1e}.txt"
        c = (imax - 1) // 2
        rr = np.hypot(*(np.mgrid[0:imax, 0:imax] - c)) * ds_m
        kp = np.where(rr.ravel() <= DISC_R + 1e-9, kx, km)
        (IN / mapname).write_text("kp\n"
                                  + "\n".join(f"{v:.6e}" for v in kp) + "\n")
        ndisc = int((kp == kx).sum())
        hdr = [
            f"! STAGE 13 -- f0 sweep on res{PARENT}.in. f0 and the DOMAIN change.",
            "!",
            "! Across 82 runs, six reproduce the observed slip to within 3% and all",
            "! six sit at +70 to +80% wellhead; nothing matches all three targets at",
            "! any tolerance down to +/-50%. The model needs ~1.75x the measured",
            "! pressure to slip, and the only quantity setting how much pressure slip",
            "! requires is Dtauc = f0*sigmabar_0 - tau_0. tau_0 is measured.",
            "!",
            "! Wang & Dunham's own front criterion (cmp_seis_extent.m) pins",
            "! Dtauc = 1.86 MPa, not tau_0. Their calibrated tau_0 = 15.0 at",
            "! f0 = 0.60 and a resolved tau_0 = 10.26 at f0 = 0.433 both reproduce",
            "! the observed 187/323/417 m front to ~5%; the data cannot distinguish",
            "! them. f0 = 0.433 is inside the published range for chlorite-bearing",
            "! granitic gouge and is the value their own off-fault spring-sliders",
            "! use, and they already argue for phyllosilicates in this fault zone.",
            "!",
            f"! f0 {ff(base['f0']):.2f} -> {f0:.4f}, so Dtauc {DT_REF:.2f} -> {dt:.2f} MPa",
            f"! and dp_crit = sigmabar_0 - tau_0/f0 goes "
            f"{sig - tau0/ff(base['f0']):.2f} -> {dpc:.2f} MPa.",
            "! tau_0 is UNCHANGED at the measured value.",
            "!",
            "! PREDICTION, from Stage 10's finding that the plateau tracks dp_crit:",
            f"! lowering f0 should drop the pressure while leaving slip roughly",
            f"! intact, so res{PARENT}'s front (1.08 on 0-18 d) and slip (0.95x at",
            "! 17 d) survive while its +58% wellhead falls. If instead the front",
            "! runs away and the slip collapses, pressure and strength are not",
            "! separable and that is the answer.",
            "!",
            "! WHAT THE CLOSED FORM CANNOT TEST: it is a Coulomb threshold on",
            "! secondary faults. Whether HBI's rate-and-state main fault with",
            "! enhancement and a real injection record behaves the same way is",
            "! unchecked -- CONCLUSIONS.md flags it.",
            "!",
            f"! DOMAIN. A weaker fault slips further, and at the parent's imax 601",
            f"! (half 1502 m) this would breach. Sizing from the empirical",
            f"! front-vs-Dtauc scaling over runs sharing this map (front ~",
            f"! Dtauc^{EXP:.2f}) projects a {fr:.0f} m front at {DAYS:.0f} d, so",
            f"! imax {imax} gives a {half:.0f} m half-domain and L/half = {ratio:.2f}.",
            "! That scaling is BARELY CONSTRAINED -- the five runs span only two",
            "! Dtauc values and scatter 0.59-0.70 at fixed Dtauc -- so it is used to",
            "! size generously, never to predict a result.",
            "!",
            f"! ds STAYS {ds_m:.0f} m. Enlarging via imax rather than ds keeps the",
            "! cell size, near-well resolution and Peaceman well index",
            "! T = 2*pi*k/eta/(log(0.2*ds/rw)+skin) identical to the parent -- and T",
            "! sets the wellhead, the quantity under test.",
            "!",
            f"! {mapname}: same {DISC_R:.0f} m disc at kpmax {kx:.1e} over background",
            f"! kpmin {km:.0e}, {ndisc} cells, regenerated for imax {imax}.",
        ]
        if not go:
            hdr += [
                "!",
                "! *** NOT LAUNCHED WITH THE REST OF STAGE 13. ***",
                f"! imax {imax} is {imax**2/601**2:.1f}x the parent's cells, so roughly a",
                "! 10-15x runtime on an 11 h baseline -- past the 2 d partition limit.",
                "! This deck exists so it can be sized properly once 632920/632921",
                "! measure the real front sensitivity. Do not submit as written.",
            ]
        ov = {"f0": f"{f0:.4f}", "imax": str(imax), "jmax": str(imax),
              "parameter_file": f'"{mapname}"',
              "tmax": f"{DAYS/365.0:.8f}", "nstep": "400000"}
        lines = [f"filenumber {new}"] + hdr
        for k, v in pairs:
            if k == "filenumber":
                continue
            lines.append(f"{k} {ov.get(k, v)}")
        (IN / f"res{new}.in").write_text("\n".join(lines) + "\n")
        made.append((new, f0, imax, mapname, go))

    print("\nverifying")
    nbad = 0
    for new, f0, imax, mapname, go in made:
        d = dict(read_deck(IN / f"res{new}.in"))
        mp = np.loadtxt(IN / mapname, skiprows=1)
        allowed = {"filenumber", "f0", "imax", "jmax", "parameter_file",
                   "tmax", "nstep"}
        diff = {k for k in set(d) | set(base) if d.get(k) != base.get(k)}
        ck = {
            "only allowed keys": not (diff - allowed),
            "f0 set": abs(ff(d["f0"]) - f0) < 1e-9,
            "tau_0 unchanged": (d["muinit"] == base["muinit"]
                                and d["sigmainit"] == base["sigmainit"]),
            "a,b,dc unchanged": all(d[k] == base[k] for k in ("a", "b", "dc")),
            "kpmax == map max": abs(mp.max() - kx) / kx < 1e-6,
            "kpmin == map min": abs(mp.min() - km) / km < 1e-6,
            "map size == imax^2": mp.size == imax ** 2,
            "ds unchanged": d["ds"] == base["ds"],
            "injection unchanged": (d.get("injection_file")
                                    == base.get("injection_file")),
        }
        ok = all(ck.values()); nbad += not ok
        print(f"  res{new}.in f0 {f0:.4f} imax {imax} "
              f"({'launch' if go else 'HOLD'}): "
              + ("OK" if ok else "PROBLEM "
                 + str([k for k, v in ck.items() if not v])))
    print(f"\n  {len(made)-nbad}/{len(made)} clean")
    print("\nSUBMIT ONLY: " + ", ".join(f"res{n}.in" for n, _, _, _, g in made if g))
    sys.exit(0 if nbad == 0 else 1)


if __name__ == "__main__":
    main()
