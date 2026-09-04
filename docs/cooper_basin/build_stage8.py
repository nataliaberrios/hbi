#!/usr/bin/env python3
"""Stage 8: reproduce Job 911 by making the WELLBORE compliant, not by editing code.

Job 911's figure (figures/taiyi_validation/) shows, at 17 d, a peak of ~6 cm and
slip out to ~+/-0.9 km with a sharp cusp at the injector. Stage 7's 632888
reproduces the peak (6.00 cm at 17 d) but is 2.6x too narrow (+/-0.34 km) and
rounded at the apex. Two reasons, both introduced by Stage 7 and neither physical:

  * 911 has permev T with kpmax 2.5e-13 over kp = kpmin = 1e-15 -- a 250x
    enhancement range. 632888 is permev F, so it has no enhancement at all, and
    632889 has only 2.5x because kpmin was raised to 1e-13 and kpmax left alone.
    Enhancement propagates a permeable conduit outward, which is the obvious
    candidate for the missing width.
  * ds went 5 m -> 20 m, so the apex is averaged over 4x wider cells and rounds.

Both followed from raising kp 100x, which Stage 7 did because kp = 1e-15 drove the
timestep to 0.01 s. That stiffness is NOT intrinsic to kp = 1e-15 -- 911 ran at
that permeability. It comes from the wellbore model, which 911 predates: at
kp = 1e-15,

    gamma = (Sw_fwid/h)/((Sw_fwid/h) + T) = 0.977        m_diffusion.f90:671

so the formation barely damps the well and the pressure ramps at essentially the
full q/Sw_fwid (m_diffusion.f90:763), 0.66 MPa/s, against 0.27 for the case that
ran fine.

THE FIX IS AN INPUT, NOT A CODE CHANGE. Sw_fwid is well storage per unit fault
width; raising it makes the wellbore more compliant and the ramp is q/Sw_fwid,
so the reduction is linear:

    Sw_fwid 7.4e-9  ->  ramp 0.660 MPa/s      (what failed)
    Sw_fwid 7.4e-7  ->  ramp 0.007 MPa/s      (40x below what worked)
    Sw_fwid 7.4e-5  ->  ramp 0.000 MPa/s

THE RISK, AND WHY THIS IS A BRACKET RATHER THAN ONE RUN. A more compliant well is
a bigger buffer. The formation source term is proportional to gamma*T*pw
(m_diffusion.f90:675) and T is tiny at kp = 1e-15, so if Sw_fwid is too large the
injected fluid fills the wellbore instead of entering the rock: numerically
healthy, but little pressure in the formation and little slip. There is therefore
a window, and its width is not known a priori. Two runs bracket it by 100x. Read
BOTH the timestep and whether slip develops -- a run that is fast and produces
nothing has failed, not succeeded.

Everything else is 911's: kp = kpmin = 1e-15, kpmax 2.5e-13, permev T, ds 5 m,
imax 601 (911 used 600, which has no centre cell). tmax 17 d to match the
reference figure rather than 911's 30.66 d.

Domain: the far field stays at kpmin, so D = 4.99e-3 m^2/s and the 17 d diffusion
length is 171 m against a 1.50 km half-domain.

Usage:  python build_stage8.py [--write]
"""
import argparse
import math
import sys
from pathlib import Path

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
PARENT = 911
DAYS = 17.0
QFILE = "const_rate_5.0e-03_17d.txt"
RATE = 5.0e-3
IWELL = JWELL = 301
SW_VARIANTS = [(632890, "7.4e-7"), (632891, "7.4e-5")]


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
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()

    pairs = read_deck(IN / f"res{PARENT}.in")
    base = dict(pairs)
    eta, phi, beta = (float(base[k].replace("d", "e"))
                      for k in ("eta", "phi", "beta"))
    kp = float(base["kp"].replace("d", "e"))
    kpmax = float(base["kpmax"].replace("d", "e"))
    ds_km = float(base["ds"])
    D = kp / (eta * phi * beta)
    L = math.sqrt(4 * D * DAYS * 86400)
    half = ds_km * 601 / 2 * 1000
    den = math.log(0.2 * ds_km * 1e3 / 0.089)
    T = 2 * math.pi * kp / eta / den

    print(f"911's physics kept: kp = kpmin = {kp:.0e}, kpmax = {kpmax:.1e} "
          f"({kpmax/kp:.0f}x range), permev {base['permev']}, ds {ds_km*1000:.0f} m")
    print(f"D(far) {D:.4f} m^2/s -> L({DAYS:.0f} d) {L:.0f} m vs half-domain "
          f"{half:.0f} m = {L/half:.2f}")
    assert L / half < 0.8, "domain too small"
    print(f"{'run':>8s} {'Sw_fwid':>10s} {'gamma':>8s} {'ramp MPa/s':>11s}")
    for n, sw in SW_VARIANTS:
        s = float(sw)
        g = (s / 60.0) / ((s / 60.0) + T)
        print(f"{n:>8d} {sw:>10s} {g:>8.4f} {g*RATE/s*1e-6:>11.4f}")

    if not a.write:
        print("\ndry run -- nothing written. re-run with --write")
        return

    assert QFILE != "june_clean.txt"
    (IN / QFILE).write_text(
        f"1\n2\n0 {DAYS*86400:.0f}\n{IWELL} {JWELL}\n{RATE:.8f} {RATE:.8f}\n")
    print(f"\nwrote {IN / QFILE}")

    common = {
        "imax": "601", "jmax": "601",
        "tmax": f"{DAYS/365.0:.8f}",
        "dtout": "0.0002",
        "nstep": "400000",
        "injectionfromfile": "T",
        "injection_file": f'"{QFILE}"',
        "rigid": "24",                 # replaces 911's unrecognised shear_mod
        "rw": "8.9e-2",
    }
    drop = {"injection", "qinj", "shear_mod"}

    made = []
    for n, sw in SW_VARIANTS:
        hdr = [
            "! STAGE 8 -- reproduce Job 911 via a COMPLIANT WELLBORE, no code change.",
            f"! 911's physics kept exactly: kp = kpmin = {kp:.0e}, kpmax {kpmax:.1e}",
            f"! ({kpmax/kp:.0f}x enhancement range), permev T, ds {ds_km*1000:.0f} m.",
            "!",
            "! WHY STAGE 7 MISSED. 632888 matches 911's peak (6.00 cm vs ~6 cm at",
            "! 17 d) but is 2.6x too narrow and rounded at the apex, because it is",
            "! permev F (no enhancement at all) at ds 20 m. 632889 has enhancement",
            "! but only 2.5x, kpmin having been raised to 1e-13. Both followed from",
            "! raising kp 100x for numerical reasons.",
            "!",
            "! THAT WAS THE WRONG FIX. kp = 1e-15 is not intrinsically stiff -- 911",
            "! ran at it. The stiffness is the WELLBORE, which 911 predates: at",
            "! kp 1e-15, gamma = (Sw_fwid/h)/((Sw_fwid/h)+T) = 0.977, so the",
            "! formation barely damps the well and pressure ramps at essentially the",
            "! full q/Sw_fwid = 0.66 MPa/s, against 0.27 for the run that completed.",
            "!",
            f"! Sw_fwid raised from 7.4e-9 to {sw}, making the wellbore more",
            "! compliant. The ramp is q/Sw_fwid, so the reduction is linear.",
            "!",
            "! RISK, STATED BEFORE THE FACT: a more compliant well is a bigger",
            "! buffer, and the formation source is proportional to gamma*T*pw with T",
            "! tiny at this kp. Too much storage and the fluid fills the wellbore",
            "! instead of the rock -- numerically healthy, but no pressure in the",
            "! formation and no slip. This deck is one of a 100x bracket",
            f"! (res632890 7.4e-7, res632891 7.4e-5) precisely because the window is",
            "! not known in advance. A fast run that produces no slip has FAILED.",
            "!",
            f"! Domain: far field stays at kpmin, so D = {D:.4f} m^2/s and the",
            f"! {DAYS:.0f} d diffusion length is {L:.0f} m against a {half:.0f} m",
            f"! half-domain ({100*L/half:.0f}%).",
            "!",
            f"! tmax {DAYS:.0f} d to match the reference figure, not 911's 30.66 d.",
            "! imax 600 -> 601 for an exact centre cell.",
        ]
        lines = [f"filenumber {n}"] + hdr
        seen = set()
        for k, v in pairs:
            if k == "filenumber" or k in drop:
                continue
            lines.append(f"{k} {common.get(k, v)}")
            seen.add(k)
        for k, v in list(common.items()) + [("Sw_fwid", sw)]:
            if k not in seen:
                lines.append(f"{k} {v}")
        (IN / f"res{n}.in").write_text("\n".join(lines) + "\n")
        made.append((n, sw))

    print("\nverifying 911's physics is untouched")
    bad = 0
    keep = ("a", "b", "dc", "f0", "muinit", "sigmainit", "kp", "kpmin", "kpmax",
            "kL", "kT", "phi", "beta", "eta", "ds", "problem", "permev",
            "limitsigma", "minsig", "parameterfromfile", "backslip", "vpl",
            "velinit", "velmin", "pfinit", "pressurediffusion",
            "bcl", "bcr", "bct", "bcb", "pbcl", "pbcr", "pbct", "pbcb")
    for n, sw in made:
        d = dict(read_deck(IN / f"res{n}.in"))
        ck = {
            "911 physics identical": all(d.get(k) == base.get(k) for k in keep),
            "permev T": d.get("permev") == "T",
            f"Sw_fwid {sw}": d.get("Sw_fwid") == sw,
            "no stale keys": not (drop & set(d)),
            "rigid 24": d.get("rigid") == "24",
        }
        ok = all(ck.values())
        bad += not ok
        print(f"  res{n}.in Sw {sw}: "
              + ("OK" if ok else "PROBLEM " + str([k for k, v in ck.items() if not v])))
    print(f"\n  {len(made)-bad}/{len(made)} clean")
    sys.exit(0 if bad == 0 else 1)


if __name__ == "__main__":
    main()
