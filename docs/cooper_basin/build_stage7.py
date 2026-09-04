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

kp CANNOT STAY AT res911.in's 1e-15 -- established empirically, not argued.
Built verbatim and submitted as 632886/632887, it drove the adaptive timestep down
to 0.01 s and held it there. Reaching tmax would have taken 252 MILLION steps
against an nstep budget of 400000, stopping at 0.05 d of 30.66. Both were
cancelled at step 240. A nearly impermeable medium under constant injection builds
well pressure faster than the integrator can follow. The same tightness is
separately on record producing 144 MPa of injector overpressure against
sigmainit 30 (632510, traced to kpmin sitting 20-1100x below the paper's
calibration), which limitsigma T / minsig 1 would then clamp -- a nonlinearity no
linear-diffusion solution contains.

So kp is raised 100x to 1e-13 and ds from 5 m to 20 m. Both are forced by the
same arithmetic: at kp 1e-13, D = kp/(eta*phi*beta) = 0.499 m^2/s, so the 30.66 d
diffusion length is 2.30 km -- 38% of the 6.01 km half-domain at ds 20 m, but 153%
of the 1.50 km half-domain at ds 5 m. The script ASSERTS L/half < 0.8, so a
combination that would run off the edge fails loudly instead of producing a
boundary-contaminated result.

Usage:  python build_stage7.py [--kp 1e-13] [--ds 0.020] [--days 30.66]
                              [--rate 5.0e-3] [--start 632888] [--write]
"""
import argparse
import sys
from pathlib import Path

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
PARENT = 911
IWELL = JWELL = 301              # exact centre of a 601 grid; res911.in used 600,
                                 # which has no centre cell -- see below

# res911.in's kp = 1e-15 is COMPUTATIONALLY INFEASIBLE here, not merely a
# linearity concern. Submitted as 632886/632887 it drove the adaptive timestep
# down to 0.01 s and held it there: 252 MILLION steps would have been needed to
# reach tmax against an nstep budget of 400000, stopping at 0.05 d of 30.66.
# Both were cancelled at step 240. A nearly impermeable medium under constant
# injection builds well pressure faster than the integrator can follow.
#
# kp is therefore raised 100x to 1e-13, which is still a low permeability and
# keeps 911's full duration inside a safe domain. Grid moves to ds 20 m for the
# same reason: at kp 1e-13, D = kp/(eta*phi*beta) = 0.499 m^2/s, so the 30.66 d
# diffusion length is 2.30 km -- 38% of a 6.01 km half-domain at ds 20 m, but
# 153% of the 1.50 km half-domain at ds 5 m.


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
    ap.add_argument("--kp", type=float, default=1e-13,
                    help="uniform initial permeability, m^2. Default 1e-13; "
                         "res911.in's 1e-15 is infeasible, see the note above.")
    ap.add_argument("--ds", type=float, default=0.020, help="cell size, km")
    ap.add_argument("--days", type=float, default=30.66,
                    help="run duration; res911.in's tmax 0.084 yr is 30.66 d")
    ap.add_argument("--start", type=int, default=632888)
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()

    DAYS = a.days
    tmax_s = DAYS * 86400.0
    qfile = f"const_rate_{a.rate:.1e}.txt"
    tmax_yr = DAYS / 365.0

    print(f"parent res{PARENT}.in -> res{a.start} (permev F) and res{a.start+1} (permev T)")
    print(f"constant rate {a.rate:.4e} for {DAYS:.2f} d "
          f"(tmax {tmax_yr:.8f} yr), well at ({IWELL},{JWELL})")
    print(f"injection file: {qfile}")

    pairs = read_deck(IN / f"res{PARENT}.in")
    base = dict(pairs)
    assert base["injectionfromfile"].upper() == "F", "parent already uses a file"
    assert base["imax"] == "600", f"parent imax is {base['imax']}, expected 600"
    eta, phi, beta = (float(base[k].replace("d", "e")) for k in ("eta", "phi", "beta"))
    D = a.kp / (eta * phi * beta)
    half_km = a.ds * 601 / 2.0
    L_km = (4 * D * tmax_s) ** 0.5 / 1000.0
    print(f"kp {a.kp:.0e} -> D {D:.4f} m^2/s;  L({DAYS:.2f} d) = {L_km:.2f} km "
          f"vs half-domain {half_km:.3f} km  = {L_km/half_km:.2f}")
    assert L_km / half_km < 0.8, (
        f"domain too small: L/half = {L_km/half_km:.2f}. Raise --ds or lower "
        f"--kp/--days.")

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

    common = {
        "imax": "601", "jmax": "601",       # 601 gives an exact centre cell
        "ds": f"{a.ds:g}",
        "kp": f"{a.kp:.1e}",
        "kpmin": f"{a.kp:.1e}",             # uniform field starts AT the floor
        "tmax": f"{tmax_yr:.8f}",
        "dtout": "0.0002",                  # ~1.75 h, as the grid uses
        "nstep": "400000",
        "injectionfromfile": "T",
        "injection_file": f'"{qfile}"',
        "rigid": "24",                      # replaces the ignored shear_mod
        "Sw_fwid": "7.4e-9",
        "rw": "8.9e-2",
    }
    drop = {"injection", "qinj", "shear_mod"}   # superseded by the file path

    made = []
    for n, pev in ((a.start, "F"), (a.start + 1, "T")):
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
                f"! kpmax. Paired with res{a.start}.in so the effect of enhancement on the",
                "! clean constant-rate problem is isolated. NOT the analytical case.",
            ]
        hdr += [
            "!",
            f"! kp RAISED from the parent's 1e-15 to {a.kp:.1e}, and ds from 5 m to",
            f"! {a.ds*1000:.0f} m. This is not a preference -- the parent's value is",
            "! COMPUTATIONALLY INFEASIBLE. Submitted verbatim as 632886/632887 it drove",
            "! the adaptive timestep to 0.01 s and held it there: 252 million steps",
            "! would have been needed to reach tmax against an nstep budget of 400000,",
            "! stopping at 0.05 d of 30.66. Both were cancelled at step 240. A nearly",
            "! impermeable medium under constant injection builds well pressure faster",
            "! than the integrator can follow. The same tightness is on record giving",
            "! 144 MPa of injector overpressure against sigmainit 30 (632510), which",
            "! limitsigma T / minsig 1 would then clamp -- a nonlinearity no",
            "! linear-diffusion solution contains.",
            "!",
            f"! Domain: D = kp/(eta*phi*beta) = {D:.4f} m^2/s, so the {DAYS:.2f} d",
            f"! diffusion length is {L_km:.2f} km against a {half_km:.3f} km",
            f"! half-domain = {100*L_km/half_km:.0f}%. At ds 5 m the half-domain would",
            f"! be 1.50 km and this would be {100*L_km/1.5025:.0f}% -- off the edge.",
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
            # kp, kpmin and ds are DELIBERATELY different from the parent -- the
            # parent's kp = 1e-15 is infeasible and forces ds too. Assert the new
            # values instead of exempting them, and assert they really did move,
            # so a silently-unchanged kp cannot slip through as "clean".
            "kp == kpmin == requested": (
                float(d.get("kp", "nan")) == a.kp
                and float(d.get("kpmin", "nan")) == a.kp),
            "kp actually raised from parent": (
                float(d.get("kp", "nan")) > float(base["kp"].replace("d", "e"))),
            "ds as requested": float(d.get("ds", "nan")) == a.ds,
            "everything else == parent": all(
                d.get(k) == base.get(k) for k in
                ("a", "b", "dc", "f0", "muinit", "sigmainit", "kpmax",
                 "kL", "kT", "phi", "beta", "eta", "problem",
                 "limitsigma", "minsig", "parameterfromfile", "backslip",
                 "vpl", "velinit", "velmin", "pfinit", "pressurediffusion",
                 "bcl", "bcr", "bct", "bcb", "pbcl", "pbcr", "pbct", "pbcb")),
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
