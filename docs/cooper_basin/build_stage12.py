#!/usr/bin/env python3
"""Stage 12: the permev F control on 632901 -- does enhancement do the work?

THE QUESTION. The hypothesis this whole grid was built to test is that matching
both Cooper Basin targets needs permeability ENHANCEMENT plus a NONUNIFORM
initial permeability, on an UNDERSTRESSED fault. Two of those three are now
isolated:

  understressed      confirmed. All 3 runs inside both bands sit at
                     tau_0 = 10.4-11.1 MPa against Wang & Dunham's 15.0, and
                     nothing at 15.0 reaches even 0.35 on the front.
  nonuniform map     confirmed. Every matching run has a two-zone map; the
                     uniform-permeability runs top out at 0.81.
  enhancement        NOT ISOLATED. Correlated only.

The correlation is strong -- all 3 matching runs are permev T, the best permev F
front across 8 runs is 0.64 and that one has a +312% wellhead, and the two
permev F runs with acceptable wellheads (632880/632881, Taiyi's own parameters)
reach only 0.14 and 0.22. There is also a clear mechanism: 632901's map disc is
150 m in radius, while its pressure plateau, slip front and permeability
transition all sit at ~380 m at 5 d, coinciding to within 15-20 m. Nothing in
the initial condition put high permeability at 380 m, so evolution carried it
there.

But none of that is a paired control, and the obvious objection is that the
nonuniform map alone did the work while enhancement came along for the ride.

THE RUN. 632901 with permev F and nothing else touched -- same map, same tau_0,
same kpmax (which becomes inert once k cannot evolve). Predictions, so this is
falsifiable rather than illustrative:

  * if enhancement is doing the work, the high-k zone stays pinned at the map's
    150 m, the sealed-disc plateau shrinks to that radius, and the front should
    collapse from 1.06x toward roughly 150/380 = 0.4x
  * if the map alone was enough, the front stays near 1.0 and the enhancement
    part of the hypothesis is wrong as stated

Either outcome is worth having, and it is the control a reviewer will ask for.

DOMAIN. With permev F the disc cannot grow, so the extent is the map's 150 m
plus diffusion through the kpmin background: D_far = 4.994e-3 m^2/s gives
L(18 d) = 176 m, so 326 m against a 1502 m half-domain, L/half = 0.22. This is
the SAFEST run of the set -- the permev T 18 d runs are 0.48-0.67 precisely
because their enhanced zone grows.

18 d to match Stage 11, so it can go straight into the same figures.

Usage:  python build_stage12.py [--write]
"""
import argparse
import math
import sys
from pathlib import Path

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
NEW, PARENT, DAYS = 632918, 632901, 18.0
NSTEP = "400000"


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
    assert str(base.get("permev", "F")).upper().startswith("T"), \
        f"res{PARENT}.in is not permev T; there is nothing to control for"
    D = ff(base["kp"]) / (ff(base["eta"]) * ff(base["phi"]) * ff(base["beta"]))
    L = math.sqrt(4 * D * DAYS * 86400)
    half = int(base["imax"]) * ff(base["ds"]) * 1000 / 2
    names = (IN / base["parameter_file"].strip('"')).read_text().split("\n", 1)[0].split()
    import numpy as np
    arr = np.loadtxt(IN / base["parameter_file"].strip('"'), skiprows=1)
    kp = arr[:, names.index("kp")]
    im = int(base["imax"])
    c = (im - 1) // 2
    rr = np.hypot(*(np.mgrid[0:im, 0:im] - c)) * ff(base["ds"]) * 1000
    r_disc = rr.ravel()[kp == kp.max()].max()

    print(f"res{NEW}.in  <-  res{PARENT}.in,  permev T -> F only")
    print(f"  map {base['parameter_file']}")
    print(f"  disc {kp.max():.1e} out to {r_disc:.0f} m, background {kp.min():.0e}")
    print(f"  kpmax {base['kpmax']} becomes INERT with permev F")
    print(f"  tau_0 {ff(base['muinit'])*ff(base['sigmainit']):.2f} MPa, unchanged")
    print(f"\n  domain: disc {r_disc:.0f} m + L({DAYS:.0f} d) {L:.0f} m = "
          f"{r_disc+L:.0f} m vs half {half:.0f} m = {(r_disc+L)/half:.2f}")
    assert (r_disc + L) / half < 0.8, "domain too small"
    print(f"  prediction: front collapses from 1.06x toward "
          f"~{r_disc/380:.1f}x if enhancement is what widens the zone")

    if not a.write:
        print("\ndry run -- nothing written. re-run with --write")
        return

    hdr = [
        "! STAGE 12 -- the permev F CONTROL on res632901.in. permev T -> F only.",
        "!",
        "! Two of the hypothesis's three parts are isolated. Understressed:",
        "! confirmed, all 3 runs inside both bands sit at tau_0 = 10.4-11.1 MPa",
        "! against Wang & Dunham's 15.0, and nothing at 15.0 reaches 0.35 on the",
        "! front. Nonuniform map: confirmed, every matching run has a two-zone map",
        "! while the uniform runs top out at 0.81. ENHANCEMENT IS NOT ISOLATED --",
        "! only correlated.",
        "!",
        "! The correlation is strong: all 3 matching runs are permev T, the best",
        "! permev F front across 8 runs is 0.64 and that run has a +312% wellhead,",
        "! and the permev F runs with acceptable wellheads (632880/632881, Taiyi's",
        "! own parameters) reach only 0.14 and 0.22. And there is a mechanism --",
        f"! this map's disc is {r_disc:.0f} m in radius while 632901's pressure",
        "! plateau, slip front and permeability transition all sit at ~380 m at",
        "! 5 d, coinciding to within 15-20 m. Nothing in the initial condition put",
        "! high permeability at 380 m, so evolution carried it there.",
        "!",
        "! But that is not a paired control, and the obvious objection is that the",
        "! nonuniform map alone did the work. This run removes only the evolution.",
        "!",
        "! FALSIFIABLE EITHER WAY:",
        f"!   enhancement matters  -> zone pinned at {r_disc:.0f} m, front collapses",
        f"!                           from 1.06x toward ~{r_disc/380:.1f}x",
        "!   map alone sufficed  -> front stays near 1.0, and the enhancement part",
        "!                          of the hypothesis is wrong as stated",
        "!",
        f"! kpmax {base['kpmax']} is left in place but is INERT with permev F.",
        "!",
        f"! Domain: the disc cannot grow, so extent is {r_disc:.0f} m + "
        f"L({DAYS:.0f} d) = {L:.0f} m",
        f"! = {r_disc+L:.0f} m against a {half:.0f} m half-domain, "
        f"L/half = {(r_disc+L)/half:.2f}. The",
        "! safest run of the set: the permev T 18 d runs are 0.48-0.67 precisely",
        "! because their enhanced zone grows.",
        "!",
        f"! {DAYS:.0f} d to match Stage 11 so it drops into the same figures.",
    ]
    ov = {"permev": "F", "tmax": f"{DAYS/365.0:.8f}", "nstep": NSTEP}
    lines = [f"filenumber {NEW}"] + hdr
    for k, v in pairs:
        if k == "filenumber":
            continue
        lines.append(f"{k} {ov.get(k, v)}")
    (IN / f"res{NEW}.in").write_text("\n".join(lines) + "\n")

    d = dict(read_deck(IN / f"res{NEW}.in"))
    diff = {k for k in set(d) | set(base) if d.get(k) != base.get(k)}
    ck = {
        "only permev/tmax/nstep/filenumber": not (
            diff - {"filenumber", "permev", "tmax", "nstep"}),
        "permev F": d["permev"] == "F",
        "map unchanged": d["parameter_file"] == base["parameter_file"],
        "kpmax unchanged": d["kpmax"] == base["kpmax"],
        "muinit unchanged": d["muinit"] == base["muinit"],
        "a, b unchanged": d["a"] == base["a"] and d["b"] == base["b"],
        "18 d": abs(ff(d["tmax"]) * 365 - DAYS) < 1e-4,
    }
    ok = all(ck.values())
    print(f"\n  res{NEW}.in: " + ("OK" if ok else "PROBLEM "
          + str([k for k, v in ck.items() if not v])))
    print(f"  changed {sorted(diff)}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
