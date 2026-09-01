#!/usr/bin/env python3
"""Stage 4: spatially graded POROSITY, the one untested direction that does not
require overriding a measurement.

Why this and not leakoff or hydraulic fracturing:

  The blocking constraint is the SHAPE of the pressure profile. The injection
  rate is measured (36 ML) and the wellhead is measured (44.73 MPa peak), which
  together pin the transmissivity and hence Dp everywhere. For a UNIFORM medium
  the whole profile is one logarithm, so the ratio Dp(r_w)/Dp(417 m) is fixed by
  geometry: solving for the k and beta that give 11 MPa at the well and 6.54 MPa
  at 417 m returns beta = 6.6e-18, eight orders below the compressibility of
  water. No uniform medium can do it.

  A GRADED medium is not subject to that constraint -- with phi(r) varying, the
  ratio becomes a functional of the grading rather than a fixed logarithm. This
  is therefore a different problem, not another point in the same family.

  Leakoff would steepen the profile (fluid leaves the fault zone, so less Dp at
  large r) -- wrong direction, though worth doing for correctness since Wang &
  Dunham neglect it. Hydraulic fracturing produces the right kind of flattening
  but needs downhole >= 101.8 MPa, i.e. a wellhead of 61.8 against a measured
  44.73 (+38%), so it cannot be invoked without also revising the stress state by
  more than the ~6 MPa that would fix the front on its own.

How it works in HBI: beta is a SCALAR in t_params and cannot vary in space, but
phi CAN -- phiG(:) is a field, it is one of the per-cell quantities the parameter
file supplies (main_LH.f90 case('phi') -> param_diff%phiG), and the solver uses
phiG in both places that matter:

    str   = beta*phiG                     m_diffusion.f90:444   storage
    cdiff = kpG/(eta*beta*phiG)           m_diffusion.f90:445   diffusivity

So LOWERING phi far from the well simultaneously lowers storage and RAISES
diffusivity there, both of which carry pressure further with less drop. That is
the flattening we need. Note they cannot be separated: phi enters both.

Physical motivation: a fault damage zone with porosity concentrated toward the
core is standard fault-zone architecture. phi is held inside 0.005-0.02
throughout, the range already recorded as defensible.

Design: 2 bases x 3 gradings = 6 runs.

  bases     632810  the best PRESSURE match in the grid (wellhead -0.5%, no slip)
            632812  the best FRONT in the grid (255 m, but wellhead 82.5 MPa)
            Both use perm_2zone_601_ds5_kmax2.5e-13.txt, so only phi changes.

  gradings  two-zone at the same 150 m radius as the permeability disc
            G1  phi 0.020 near / 0.005 far   4x, the full physical range
            G2  phi 0.010 near / 0.005 far   2x
            G3  phi 0.005 near / 0.020 far   INVERSE -- a control that should
                                             make the profile steeper and the
                                             front worse. If it does not, the
                                             mechanism is not what we think.

Usage:  python build_stage4.py [--check]
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")

BASES = [632810, 632812]
GRADINGS = [("G1", 0.020, 0.005), ("G2", 0.010, 0.005), ("G3", 0.005, 0.020)]
R_DISC = 150.0          # m, same radius as the permeability disc
START = 632870


def read_deck(p):
    out = []
    for line in Path(p).read_text().splitlines():
        if line.startswith("!") or not line.strip():
            continue
        w = line.split()
        if len(w) >= 2:
            out.append((w[0], " ".join(w[1:])))
    return out


def write_2col(outname, permfile, imax, ds_m, phi_near, phi_far):
    """Two-column parameter file: kp taken from the existing map, phi graded.

    Format the reader expects (main_LH.f90): a header line of ncol names, then
    NCELLg rows of ncol values, i outer / j inner. Column order follows the
    header, so 'kp phi' maps column 1 -> kp, column 2 -> phi.
    """
    kp = np.loadtxt(IN / permfile, skiprows=1)
    assert kp.size == imax * imax, f"{permfile} has {kp.size}, expected {imax**2}"
    c = (imax + 1) / 2.0
    lines = ["kp phi"]
    n_near = 0
    for i in range(1, imax + 1):
        for j in range(1, imax + 1):
            r = ds_m * math.hypot(i - c, j - c)
            ph = phi_near if r <= R_DISC else phi_far
            if r <= R_DISC:
                n_near += 1
            lines.append(f"{kp[(i-1)*imax + (j-1)]:.6e} {ph:.6e}")
    (IN / outname).write_text("\n".join(lines) + "\n")
    return n_near


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    n = START
    made = []
    print(f"{'run':>7s} {'base':>7s} {'grade':>6s} {'phi near':>9s} {'phi far':>8s} "
          f"{'beta':>11s} {'phi*beta near':>14s} {'phi*beta far':>13s} {'D far ratio':>12s}")
    print("-" * 104)
    for base in BASES:
        pairs = read_deck(IN / f"res{base}.in")
        d = dict(pairs)
        imax = int(d["imax"])
        ds_m = float(d["ds"]) * 1000.0
        permfile = d["parameter_file"].strip('"')
        beta = float(d["beta"])
        phi0 = float(d["phi"])
        for gname, ph_near, ph_far in GRADINGS:
            tag = permfile.replace("perm_", "").replace(".txt", "")
            outname = f"permphi_{tag}_{gname}.txt"
            print(f"{n:>7d} {base:>7d} {gname:>6s} {ph_near:>9.3f} {ph_far:>8.3f} "
                  f"{beta:>11.4e} {ph_near*beta:>14.3e} {ph_far*beta:>13.3e} "
                  f"{phi0/ph_far:>11.2f}x")
            if not a.check:
                nn = write_2col(outname, permfile, imax, ds_m, ph_near, ph_far)
                r_eq = math.sqrt(nn * ds_m**2 / math.pi)
                changes = {"filenumber": str(n),
                           "parameter_file": f'"{outname}"',
                           "parameter_file_ncol": "2"}
                hdr = [
                    "! STAGE 4 -- spatially GRADED POROSITY.",
                    f"! Built from res{base}.in; only the parameter file and ncol change.",
                    "!",
                    "! WHY: the blocking constraint is the SHAPE of the pressure profile.",
                    "! The injection rate and the wellhead are both measured, which pins the",
                    "! transmissivity and hence Dp everywhere. For a UNIFORM medium the whole",
                    "! profile is one logarithm, so Dp(r_w)/Dp(417 m) is fixed by geometry --",
                    "! solving for the k and beta giving 11 MPa at the well and 6.54 MPa at",
                    "! 417 m returns beta = 6.6e-18, eight orders below water. A GRADED medium",
                    "! is not subject to that constraint.",
                    "!",
                    "! beta is a scalar in HBI and cannot vary in space, but phi can:",
                    "!   str   = beta*phiG              m_diffusion.f90:444   storage",
                    "!   cdiff = kpG/(eta*beta*phiG)    m_diffusion.f90:445   diffusivity",
                    "! so lowering phi far from the well lowers storage AND raises diffusivity",
                    "! there, both carrying pressure further with less drop. They cannot be",
                    "! separated -- phi enters both.",
                    "!",
                    f"! GRADING {gname}: phi = {ph_near:.3f} within {R_DISC:.0f} m of the well,",
                    f"!   {ph_far:.3f} beyond ({nn} cells near, equivalent radius {r_eq:.0f} m).",
                    "!   Both values are inside the 0.005-0.02 range recorded as defensible."
                    + ("" if gname != "G3" else
                       "\n! G3 is the INVERSE grading, included as a CONTROL: it should steepen"
                       "\n! the profile and make the front worse. If it does not, the mechanism"
                       "\n! is not what we think it is."),
                    "!",
                    f"! phi*beta ranges {ph_far*beta:.3e} (far) to {ph_near*beta:.3e} (near),",
                    f"! against the scalar {phi0*beta:.3e} in the base run.",
                    "!",
                    f"! kp column is copied verbatim from {permfile}, so the permeability",
                    "! field is IDENTICAL to the base run and only porosity differs.",
                ]
                lines = [f"filenumber {n}"] + hdr
                for k, v in pairs:
                    if k == "filenumber":
                        continue
                    lines.append(f"{k} {changes.get(k, v)}")
                    if k == "parameter_file" and "parameter_file_ncol" not in dict(pairs):
                        lines.append("parameter_file_ncol 2")
                (IN / f"res{n}.in").write_text("\n".join(lines) + "\n")
                made.append((n, base, outname))
            n += 1

    if a.check:
        print("\n--check: nothing written")
        return

    print(f"\nwrote {len(made)} decks and {len(made)} two-column parameter files\n")
    print("verifying: kp column identical to the base map, phi as designed,")
    print("and the deck changed only in filenumber / parameter_file / ncol")
    bad = 0
    for num, base, pf in made:
        bd = dict(read_deck(IN / f"res{base}.in"))
        nd = dict(read_deck(IN / f"res{num}.in"))
        diff = {k for k in set(bd) | set(nd) if bd.get(k) != nd.get(k)}
        extra = diff - {"filenumber", "parameter_file", "parameter_file_ncol"}
        v = np.loadtxt(IN / pf, skiprows=1)
        kp_base = np.loadtxt(IN / bd["parameter_file"].strip('"'), skiprows=1)
        kp_ok = np.allclose(v[:, 0], kp_base, rtol=1e-6)
        ncol_ok = nd.get("parameter_file_ncol") == "2"
        ok = (not extra) and kp_ok and ncol_ok
        if not ok:
            bad += 1
        print(f"  res{num}.in  changed {sorted(diff)}  kp identical {kp_ok}  "
              f"ncol=2 {ncol_ok}  phi {v[:,1].min():.3f}-{v[:,1].max():.3f}  "
              f"{'OK' if ok else 'PROBLEM ' + str(sorted(extra))}")
    print(f"\n  {len(made)-bad}/{len(made)} clean")
    sys.exit(0 if bad == 0 else 1)


if __name__ == "__main__":
    main()
