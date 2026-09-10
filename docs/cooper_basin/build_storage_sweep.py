#!/usr/bin/env python3
"""CYCLE 2, Round 4: the STORAGE sweep. Why the first two days build too fast.

THE OBSERVATION THIS ADDRESSES. Measured against the model, at tau_0 = 11.53
(632961, the best joint cell):

    t        q L/s    observed dp    model dp    model/obs
    0.25 d    13.2       0.16          1.51        9.4x
    1.00 d    18.5       1.61          4.07        2.5x
    2.50 d    40.3       9.55         12.28        1.29x
    4-8 d     23-27     10.81         13.4         1.25x

If the model's INJECTIVITY were wrong it would be off by a constant factor at
all times. It is not -- by the plateau it is only 1.25x high, which is the +25%
already known. The entire extra error lives in the first two days, so it is a
TRANSIENT: the model pressurises too fast. Fitted shapes over 0.15-2.5 d:
observed dp ~ t^1.72 (R2 0.986), arm 2 dp ~ t^0.91. The data starts nearly flat
and accelerates; the model climbs immediately.

Too-fast pressurisation with the right steady state means TOO LITTLE STORAGE.
Early on the data absorbs roughly ten times more fluid per MPa than the model.

WHY NOT SWEEP beta, WHICH IS THE OBVIOUS CHOICE. beta does two jobs at once:
storage as phi*beta, and diffusivity as D = kp/(eta phi beta). Raising it slows
the buildup, which is wanted, AND slows the front, which is not. That is exactly
what happened to arm 3: phi*beta x7 brought dp at 1 d to 2.24 against the
observed 1.61 -- the right direction -- while lambda collapsed to 0.08-0.56.
So arm 3 has already established the direction and the cost. Repeating it under
another name would add nothing.

TWO SWEEPS INSTEAD, each ONE KEY off its own parent, and each run at THREE
tau_0 -- 10.36, 11.53 and 12.71 MPa (parents 632960/632961/632962). Three
rather than one because storage and tau_0 both move the front, and only a grid
separates them: if the same storage fixes the early curvature at every tau_0,
the mechanism is storage; if the storage needed varies with tau_0, it is not.
18 runs, 633060-633077.

  A. phi (x3 per tau_0), 0.01 -> 0.02 / 0.05 / 0.10, i.e. phi*beta x2/x5/x10.
     phi enters the governing equations ONLY as the product phi*beta
     (setup_model.m:101 notes this), so for storage it is interchangeable with
     beta -- but it leaves a well-known fluid property alone and varies the one
     that is genuinely uncertain. The porosity of a cycle-1 stimulated fracture
     network is not 1% in any meaningful sense.

     Storage changes the TRANSIENT and not the steady plateau, which is set by
     k. So this should slow the early rise while leaving the 4-8 d value near
     its current 1.25x. It will also slow the front, since D falls as 1/phi;
     lambda goes roughly as 1/sqrt(phi*beta), so x2 costs about 0.71x of
     lambda. THAT IS ACCEPTED HERE ON PURPOSE: this round asks whether storage
     reproduces the early curvature at all. If it does, tau_0 is re-balanced
     afterwards to recover the front. Mechanism first, fit second.

     phi = 0.10 is not a granite matrix porosity and is not offered as one. It
     is an effective storativity for a stimulated fracture volume, and if the
     answer turns out to need it, that is a statement about the model's
     near-well representation rather than about rock.

  B. Sw_fwid (x3 per tau_0), 7.4e-9 -> x10 / x100 / x1000. THE CLEANEST
     EXPERIMENT AVAILABLE HERE, because it cannot touch the front. Sw_fwid is
     wellbore storage over fault width and enters only the well-cell weighting
     gamma = (Sw_fwid/h)/((Sw_fwid/h) + T) at m_diffusion.f90:671. The
     permeability field, the diffusivity and the stress state are untouched.

     Its present value is V_w*beta_w/fwid with V_w = pi R_w^2 H_w = 101 m^3 and
     beta_w = 4.41e-10 -- WATER COMPRESSIBILITY ONLY. But the well had just been
     vented: through the first shut-in the wellhead fell from 39.83 to
     -0.23 MPa in 26 minutes. A partly drained wellbore and near-well fracture
     network refills at almost no pressure rise, which is a free-surface
     capacity rather than a compressibility, and can exceed the compressive
     value by orders of magnitude. Hence x1000 as the upper bracket.

     If the first two days are the WELL refilling rather than the ROCK
     pressurising, this isolates it and costs nothing elsewhere.

A CAVEAT ON THE DATA, so this is not over-fitted. Observed dp of 0.16 MPa at
6 h while taking 13.2 L/s is very low. The gauge had just come off a vent, so
the first ~0.25 d may be recording the wellbore refilling rather than the
reservoir -- which would also be storage, but in the well rather than the rock,
i.e. sweep B rather than sweep A. The two sweeps are distinguishable precisely
because B cannot move the front and A must.

DOMAIN. Both sweeps are safe. Raising phi lowers D and therefore SHRINKS the
front, so each parent's MEASURED 13.1 d front bounds its own block -- 831 m at
tau_0 10.36, 932 at 11.53, 1082 at 12.71. Sw_fwid does not enter the
diffusivity at all, so those runs inherit the parent's front exactly. Worst
cell is L/half 0.49.

Usage:
  python build_storage_sweep.py            # dry run
  python build_storage_sweep.py --write
"""
import argparse
import math
import sys
from pathlib import Path

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
DAYS, SAFE = 13.1, 0.8
# (parent, tau_0, measured 13.1 d front in m, first filenumber for its block).
# Run at THREE tau_0 rather than one, because storage and tau_0 both move the
# front and only a grid separates them: if the early curvature is fixed at
# every tau_0 by the same storage, the mechanism is storage; if the needed
# storage varies with tau_0, it is not.
PARENTS = [(632960, 10.36, 831.0, 633060),
           (632961, 11.53, 932.0, 633066),
           (632962, 12.71, 1082.0, 633072)]
PHI_VALS = [0.02, 0.05, 0.10]        # phi*beta x2 / x5 / x10
SWF_FACS = [10.0, 100.0, 1000.0]


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


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args(argv)

    rows, half = [], None
    print(f"{'run':>7} {'parent':>7} {'tau_0':>6} {'key':>9} {'value':>10} "
          f"{'phi*beta':>10} {'D_disc':>8} {'lam f':>6} {'R bound':>8} "
          f"{'L/half':>7}")
    for PARENT, TAU, R_MAX_PARENT, FIRST in PARENTS:
        pairs = read_deck(IN / f"res{PARENT}.in")
        base = dict(pairs)
        ds = ff(base["ds"]) * 1000.0
        IMAX = int(base["imax"]); half = IMAX * ds / 2
        phi0, beta = ff(base["phi"]), ff(base["beta"])
        eta, kx = ff(base["eta"]), ff(base["kpmax"])
        swf0 = ff(base["Sw_fwid"])
        sig, f0, mu = ff(base["sigmainit"]), ff(base["f0"]), ff(base["muinit"])
        assert abs(f0 - 0.60) < 1e-9, f"res{PARENT}: f0 is {f0}"
        assert abs(mu * sig - TAU) < 0.02, \
            f"res{PARENT}: tau_0 is {mu*sig:.2f}, expected {TAU}"
        n = FIRST
        for phi in PHI_VALS:
            pb = phi * beta
            D = kx / (eta * pb)
            lf = math.sqrt(phi0 / phi)     # lambda ~ sqrt(D)
            R = R_MAX_PARENT * lf
            rows.append((n, PARENT, TAU, pairs, "phi", f"{phi:g}", pb, D, lf,
                         R, phi0, beta, swf0, half, R_MAX_PARENT, mu, sig))
            print(f"{n:>7} {PARENT:>7} {TAU:>6.2f} {'phi':>9} {phi:>10g} "
                  f"{pb:>10.3e} {D:>8.3f} {lf:>6.2f} {R:>7.0f}m "
                  f"{(R+398)/half:>7.2f}")
            n += 1
        for fac in SWF_FACS:
            v = swf0 * fac
            rows.append((n, PARENT, TAU, pairs, "Sw_fwid", f"{v:.2e}",
                         phi0*beta, kx/(eta*phi0*beta), 1.0, R_MAX_PARENT,
                         phi0, beta, swf0, half, R_MAX_PARENT, mu, sig))
            print(f"{n:>7} {PARENT:>7} {TAU:>6.2f} {'Sw_fwid':>9} {v:>10.2e} "
                  f"{phi0*beta:>10.3e} {kx/(eta*phi0*beta):>8.3f} {1.0:>6.2f} "
                  f"{R_MAX_PARENT:>7.0f}m {(R_MAX_PARENT+398)/half:>7.2f}")
            n += 1
    bad = sum(1 for r in rows if (r[9] + 398) / r[13] >= SAFE)
    if bad:
        sys.exit(f"\n{bad} deck(s) breach L/half {SAFE}. Nothing written.")
    print(f"\n  all {len(rows)} within L/half < {SAFE}; the phi runs SHRINK the "
          f"front so the\n  parent's measured 932 m bounds them, and Sw_fwid "
          f"does not enter D at all")

    if not a.write:
        print("\ndry run -- nothing written. re-run with --write")
        return

    made = []
    for (n, PARENT, TAU, pairs, key, val, pb, D, lf, R, phi0, beta, swf0,
         half, R_MAX_PARENT, mu, sig) in rows:
        if key == "phi":
            what = [
                f"! CYCLE 2 ROUND 4, SWEEP A: phi = {val}. ONE KEY off "
                f"res{PARENT}.in.",
                "!",
                "! THE MODEL PRESSURISES TOO FAST IN THE FIRST TWO DAYS. Against",
                f"! res{PARENT}.in (tau_0 {mu*sig:.2f}, the best joint cell):",
                "!     t        q L/s   obs dp   sim dp   ratio",
                "!     0.25 d    13.2     0.16     1.51    9.4x",
                "!     1.00 d    18.5     1.61     4.07    2.5x",
                "!     2.50 d    40.3     9.55    12.28    1.29x",
                "!     4-8 d    23-27    10.81    13.4     1.25x",
                "! A wrong INJECTIVITY would be off by a constant factor. This is",
                "! not: by the plateau it is only 1.25x high. The error is a",
                "! TRANSIENT, and observed dp ~ t^1.72 against the model's t^0.91.",
                "! Too-fast pressurisation at the right steady state means too",
                "! little STORAGE -- early on the data absorbs ~10x more fluid",
                "! per MPa.",
                "!",
                "! phi AND NOT beta. phi enters only as the product phi*beta, so",
                "! for storage the two are interchangeable, but phi is the one",
                "! that is genuinely uncertain -- the porosity of a cycle-1",
                "! stimulated fracture network is not 1% in any useful sense --",
                "! while beta is a known fluid property.",
                f"!   phi {phi0} -> {val}, so phi*beta {phi0*beta:.3e} -> {pb:.3e}",
                "!",
                "! WHAT IS EXPECTED, both ways. Storage changes the TRANSIENT and",
                "! not the steady plateau, which is set by k -- so the early rise",
                "! should slow while the 4-8 d value stays near 1.25x. It also",
                f"! lowers D to {D:.3f} m2/s, so lambda falls by about "
                f"{lf:.2f}x.",
                "! THAT COST IS ACCEPTED DELIBERATELY. This round asks whether",
                "! storage reproduces the early curvature at all; if it does,",
                "! tau_0 is re-balanced afterwards to recover the front.",
                "! Mechanism first, fit second.",
                "!",
                "! ARM 3 ALREADY BOUNDS THIS. It ran phi*beta x7 and brought dp at",
                "! 1 d to 2.24 against the observed 1.61 -- the right direction --",
                "! while lambda collapsed to 0.08-0.56. So the direction and the",
                "! cost are both known; what is not known is the magnitude needed.",
                "!",
                f"! phi = 0.10 is NOT offered as a granite matrix porosity. It is an",
                "! effective storativity for a stimulated fracture volume, and if",
                "! the answer needs it, that is a statement about the model's",
                "! near-well representation rather than about rock.",
                "!",
                f"! DOMAIN: raising phi SHRINKS the front. The parent measured",
                f"! {R_MAX_PARENT:.0f} m at {DAYS} d, so {R:.0f} m bounds this run;",
                f"! plus 398 m far-field against {half:.0f} m half-domain gives",
                f"! L/half {(R+398)/half:.2f}.",
            ]
            ov = {"phi": val}
        else:
            what = [
                f"! CYCLE 2 ROUND 4, SWEEP B: Sw_fwid = {val}. ONE KEY off "
                f"res{PARENT}.in.",
                "!",
                "! Same problem as sweep A addresses -- see res633060.in for the",
                "! measured discrepancy. This is the CLEANEST TEST OF IT, because",
                "! Sw_fwid CANNOT TOUCH THE FRONT. It is wellbore storage over",
                "! fault width and enters only the well-cell weighting",
                "!     gamma = (Sw_fwid/h) / ((Sw_fwid/h) + T)",
                "! at m_diffusion.f90:671. The permeability field, the diffusivity",
                "! and the stress state are all untouched.",
                "!",
                f"! Its present value {swf0:.2e} is V_w*beta_w/fwid with",
                "! V_w = pi R_w^2 H_w = 101 m^3 and beta_w = 4.41e-10, i.e. WATER",
                "! COMPRESSIBILITY ONLY. But the well had just been vented: through",
                "! the first shut-in the wellhead fell from 39.83 to -0.23 MPa in",
                "! 26 minutes. A partly drained wellbore and near-well fracture",
                "! network refills at almost no pressure rise -- a free-surface",
                "! capacity, not a compressibility -- and can exceed the",
                f"! compressive value by orders of magnitude. Hence up to x1000.",
                "!",
                "! IF THE FIRST TWO DAYS ARE THE WELL REFILLING RATHER THAN THE",
                "! ROCK PRESSURISING, this isolates it at no cost elsewhere. The",
                "! two sweeps are distinguishable precisely because B cannot move",
                "! the front and A must.",
                "!",
                "! A CAVEAT ON THE DATA, so this is not over-fitted: observed dp of",
                "! 0.16 MPa at 6 h while taking 13.2 L/s is very low, and the gauge",
                "! had just come off a vent. The first ~0.25 d may be recording the",
                "! wellbore refilling rather than the reservoir.",
                "!",
                f"! DOMAIN: unchanged from the parent, {R_MAX_PARENT:.0f} m at "
                f"{DAYS} d, L/half {(R_MAX_PARENT+398)/half:.2f}.",
            ]
            ov = {"Sw_fwid": val}

        lines = [f"filenumber {n}"] + what
        for k, v in pairs:
            if k == "filenumber":
                continue
            lines.append(f"{k} {ov.get(k, v)}")
        (IN / f"res{n}.in").write_text("\n".join(lines) + "\n")
        made.append((n, PARENT, key, val))

    print("\nverifying each deck against ITS OWN parent")
    nbad = 0
    for n, PARENT, key, val in made:
        ref = dict(read_deck(IN / f"res{PARENT}.in"))
        d = dict(read_deck(IN / f"res{n}.in"))
        diff = sorted(k for k in set(ref) | set(d)
                      if ref.get(k) != d.get(k) and k != "filenumber")
        ok = diff == [key]
        nbad += not ok
        print(f"  {n}  vs res{PARENT}  changed {diff}  "
              + ("OK" if ok else f"EXPECTED ['{key}']"))
    if nbad:
        sys.exit(f"\n{nbad} problem(s). Nothing submitted.")
    print(f"\nall {len(made)} verified. Submit from {IN}:")
    for n, PARENT, key, val in made:
        d = dict(read_deck(IN / f"res{n}.in"))
        print(f"  sbatch march26_submit_hbi_git_scratch.sh -i res{n}.in "
              f"-w {d['injection_file'].strip(chr(34))} "
              f"-p {d['parameter_file'].strip(chr(34))}")


if __name__ == "__main__":
    main()
