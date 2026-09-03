#!/usr/bin/env python3
"""Stage 7: CONSTANT-RATE injection on corrected code, for analytical comparison.

This is the run res911.in was meant to be. res911.in cannot simply be rerun,
because it has three silent defects -- two mis-specified keys and one that means
it never injected at all.

DEFECT 1 -- res911.in never injected. It sets `injectionfromfile F`,
`injection flowrate`, `qinj 4.2e-3`. For `problem 3dp` the pressure solver is
Beuler2d, and its constant-rate branch is COMMENTED OUT and marked
"No longer correct!":

    ! else if(param_diff%injection=='flowrate' .and. time<...) then
      ! b(imax/2,jmax/2)=...+h*param_diff%qinj/str(...)*1e-12/ds0/ds0
    ! end if                                   m_diffusion.f90:677-680

Only Beuler1d/Beuler1db, used for 2D-fault problems, still have it active
(m_diffusion.f90:234 and :347). So a 3dp run with injectionfromfile F gets no
source term whatsoever: pf stays at pfinit and nothing drives slip.

DEFECT 2 -- `shear_mod 24` is not a recognised key. read_inputfile has only
`case('rigid')` (main_LH.f90:2695). Unrecognised keys are silently ignored, so
res911.in ran at the default rigid = 32.04 GPa (m_const.f90:5) rather than the
24 it asked for -- a 33% error in shear modulus, and 24 is Taiyi's value.

DEFECT 3 -- `Sw_fwid` is absent, and it defaults to 1.0 (main_LH.f90:180)
against the 7.4e-9 every grid deck uses. It sets the wellbore coupling
gamma = Sw_fwid/h/(Sw_fwid/h + T) (m_diffusion.f90:671) and divides the injection
term (m_diffusion.f90:675, :763), so eight orders of magnitude there is not a
detail. `rw` is also absent but its default, 0.089, happens to be the value in
use, so that one is harmless.

THE FIX, WITH NO CODE CHANGE. Drive a genuinely constant rate through the
injection-FILE path, which is the tested one and which carries the Peaceman well
model. The reader (input_well, m_diffusion.f90:49) wants

    nwell
    npoint
    times(1:npoint)                     seconds
    iwell jwell                         per well
    qvals(1:npoint)                     per well

and the interpolation is linear between points, holding the last value beyond the
end (m_diffusion.f90:650-660). So two points with EQUAL rates give an exactly
constant rate for all time. june_clean.txt is left untouched; this writes a new
file.

WHY TWO RUNS.

  permev F  Linear diffusion, constant permeability. This is the one that maps
            onto a closed-form solution -- with permev T the diffusivity evolves
            with slip and no analytical form applies. This is the primary run.
  permev T  res911.in's own setting, so the pair isolates what enhancement does
            to the clean constant-rate problem.

A CAVEAT ON kp = 1e-15, which is res911.in's value and is kept here for
faithfulness. It is very tight, and this lineage is on record producing an
injector overpressure of 144 MPa against sigmainit 30 (documented for 632510,
traced to kpmin sitting 20-1100x below the paper's calibration). With
`limitsigma T` and `minsig 1` the effective normal stress will clamp, which is a
nonlinearity a linear-diffusion analytical solution does not contain. If the
comparison needs to stay in the linear regime, rerun at a larger kp -- 4e-13, the
paper's far-field value, is the obvious choice, and needs only the kp line
changed.

Domain is safe either way: at kp 1e-15, D = kp/(eta*phi*beta) = 4.99e-3 m^2/s, so
the 30.66 d diffusion length is 230 m against a 1.5 km half-domain (15%). Under
enhancement the far field stays at kpmin, so it does not move.

Usage:  python build_stage7.py [--rate 5.0e-3] [--write]
"""
import argparse
import sys
from pathlib import Path

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
PARENT = 911
DAYS = 30.66                     # res911.in's tmax 0.084 yr, HBI's year = 365 d
IWELL = JWELL = 301              # exact centre of a 601 grid; res911.in used 600,
                                 # which has no centre cell -- see below


def read_deck(p):
    out = []
    for line in Path(p).read_text().splitlines():
        if line.startswith("!") or not line.strip():
            continue
        w = line.split()
        if len(w) >= 2:
            out.append((w[0], " ".join(w[1:])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rate", type=float, default=5.0e-3,
                    help="constant rate in june_clean.txt's units (per unit fault "
                         "width). Default 5.0e-3 ~ the measured record's mean of "
                         "4.937e-3, so magnitudes stay comparable to the grid.")
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()

    tmax_s = DAYS * 86400.0
    qfile = f"const_rate_{a.rate:.1e}.txt"
    tmax_yr = DAYS / 365.0

    print(f"parent res{PARENT}.in -> res632886 (permev F) and res632887 (permev T)")
    print(f"constant rate {a.rate:.4e} for {DAYS:.2f} d "
          f"(tmax {tmax_yr:.8f} yr), well at ({IWELL},{JWELL})")
    print(f"injection file: {qfile}")
    if not a.write:
        print("\ndry run -- nothing written. re-run with --write")
        return

    # Two points with equal rates -> exactly constant. Never overwrite june_clean.
    assert qfile != "june_clean.txt"
    (IN / qfile).write_text(
        "1\n2\n"
        f"0 {tmax_s:.0f}\n"
        f"{IWELL} {JWELL}\n"
        f"{a.rate:.8f} {a.rate:.8f}\n")
    print(f"wrote {IN / qfile}")

    pairs = read_deck(IN / f"res{PARENT}.in")
    base = dict(pairs)
    assert base["injectionfromfile"].upper() == "F", "parent already uses a file"
    assert base["imax"] == "600", f"parent imax is {base['imax']}, expected 600"

    common = {
        "imax": "601", "jmax": "601",       # 601 gives an exact centre cell
        "tmax": f"{tmax_yr:.8f}",
        "dtout": "0.0002",                  # ~1.75 h, as the grid uses
        "nstep": "400000",                  # headroom for 30.66 d
        "injectionfromfile": "T",
        "injection_file": f'"{qfile}"',
        "rigid": "24",                      # replaces the ignored shear_mod
        "Sw_fwid": "7.4e-9",
        "rw": "8.9e-2",
    }
    drop = {"injection", "qinj", "shear_mod"}   # superseded by the file path

    made = []
    for n, pev in ((632886, "F"), (632887, "T")):
        hdr = [
            "! STAGE 7 -- CONSTANT-RATE injection on corrected code.",
            f"! The run res{PARENT}.in was meant to be. res{PARENT}.in has three",
            "! silent defects and cannot just be rerun:",
            "!",
            "!  1. IT NEVER INJECTED. injectionfromfile F + problem 3dp routes to",
            "!     Beuler2d, whose constant-rate branch is commented out and marked",
            '!     "No longer correct!" (m_diffusion.f90:677-680). Only the 1D',
            "!     solvers still have it (:234, :347). No source term at all.",
            "!  2. shear_mod is not a recognised key -- read_inputfile has only",
            "!     case('rigid') (main_LH.f90:2695) -- so it ran at the default",
            "!     rigid = 32.04 (m_const.f90:5), not the 24 it asked for.",
            "!  3. Sw_fwid was absent and defaults to 1.0 (main_LH.f90:180) against",
            "!     the 7.4e-9 in use, eight orders out. It sets the wellbore",
            "!     coupling gamma (m_diffusion.f90:671) and divides the injection",
            "!     term (:675, :763).",
            "!",
            f"! Constant rate is delivered through the injection FILE path instead,",
            f"! which is tested and carries the Peaceman well model: {qfile} holds",
            f"! two points with equal rates, and the interpolation holds the last",
            f"! value beyond the end (m_diffusion.f90:650-660), so the rate is",
            f"! exactly {a.rate:.4e} throughout. june_clean.txt is untouched.",
            "!",
            f"! imax/jmax 600 -> 601 so there IS an exact centre cell for the well;",
            f"! an even grid has no centre, which matters for comparison against an",
            f"! axisymmetric analytical solution.",
            "!",
        ]
        if pev == "F":
            hdr += [
                "! permev F -- CONSTANT permeability, so this is LINEAR diffusion and",
                "! the run maps onto a closed-form solution. This is the primary run",
                "! for the analytical comparison; with permev T the diffusivity",
                "! evolves with slip and no analytical form applies.",
            ]
        else:
            hdr += [
                f"! permev T -- res{PARENT}.in's own setting, kp = kpmin growing toward",
                "! kpmax. Paired with res632886.in so the effect of enhancement on the",
                "! clean constant-rate problem is isolated. NOT the analytical case.",
            ]
        hdr += [
            "!",
            "! CAVEAT ON kp 1e-15, kept from the parent: it is very tight, and this",
            "! lineage is on record producing 144 MPa of injector overpressure",
            "! against sigmainit 30 (632510). limitsigma T / minsig 1 will then clamp",
            "! the effective normal stress -- a nonlinearity no linear-diffusion",
            "! solution contains. If the comparison must stay linear, raise kp to",
            "! 4e-13 (the paper's far-field value); only that one line need change.",
            "!",
            f"! Domain: D = kp/(eta*phi*beta) = 4.99e-3 m^2/s, so the {DAYS:.2f} d",
            "! diffusion length is 230 m against a 1.5 km half-domain (15%).",
        ]
        lines = [f"filenumber {n}"] + hdr
        seen = set()
        for k, v in pairs:
            if k == "filenumber" or k in drop:
                continue
            if k == "permev":
                lines.append(f"permev {pev}")
                seen.add(k)
                continue
            lines.append(f"{k} {common.get(k, v)}")
            seen.add(k)
        for k, v in common.items():
            if k not in seen:
                lines.append(f"{k} {v}")
        (IN / f"res{n}.in").write_text("\n".join(lines) + "\n")
        made.append((n, pev))

    print("\nverifying")
    bad = 0
    for n, pev in made:
        d = dict(read_deck(IN / f"res{n}.in"))
        checks = {
            "permev": d.get("permev") == pev,
            "injectionfromfile T": d.get("injectionfromfile") == "T",
            "injection_file": d.get("injection_file", "").strip('"') == qfile,
            "no stale qinj/injection/shear_mod": not (drop & set(d)),
            "rigid 24": d.get("rigid") == "24",
            "Sw_fwid": d.get("Sw_fwid") == "7.4e-9",
            "imax 601": d.get("imax") == "601" and d.get("jmax") == "601",
            "physics == parent": all(
                d.get(k) == base.get(k) for k in
                ("a", "b", "dc", "f0", "muinit", "sigmainit", "kp", "kpmax",
                 "kpmin", "kL", "kT", "phi", "beta", "eta", "ds", "problem",
                 "limitsigma", "minsig", "parameterfromfile")),
        }
        ok = all(checks.values())
        bad += not ok
        print(f"  res{n}.in permev {pev}: "
              + ("OK" if ok else "PROBLEM " + str([k for k, v in checks.items()
                                                   if not v])))
    print(f"\n  {len(made)-bad}/{len(made)} clean")
    sys.exit(0 if bad == 0 else 1)


if __name__ == "__main__":
    main()
