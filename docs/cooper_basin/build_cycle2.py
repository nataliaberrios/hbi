#!/usr/bin/env python3
"""CYCLE 2 as a new project: t = 0 is data-day 4.300. Round 1, the muinit sweep.

86 runs failed to reproduce the Cooper Basin slip front, wellhead pressure and
slip magnitude together over the full 17 d record. Rather than keep perturbing
that, this treats the post-shut-in period as a NEW problem: t = 0 is data-day
4.300 where injection resumes, cycle 1 is DISCARDED rather than inherited, and
the model's own parameters are re-swept against the post-4.3 data.

WHY THE POST-4.3 RECORD IS THE BETTER TARGET. Cycle 1 is all transient -- two
wellhead spikes at 43.15 and 44.73 MPa, a bled-off gauge from 1.59 to 3.56 d and
a seismicity gap. Post-4.3 is a smooth ramp from 32.6 MPa to ~43.7 by 8 d, a
FOUR-DAY QUASI-STEADY PLATEAU at 43.7 MPa (8-12 d), then a step to ~48-52 MPa.
A sustained plateau constrains injectivity directly; two spikes never could.

It is NOT the gentler half -- it sustains higher pressure than cycle 1, plus a
two-sample 87.76 MPa transient at 14.23 d (92% of Sv, almost certainly
water-hammer; to be masked, not fitted).

CUTTING HERE IS PHYSICALLY JUSTIFIED. Through the shut-in (1.58-4.30 d) the
catalogue drops from 615 to 9.6 events/day and the front does NOT advance --
318 m during, against 418 m already reached. Saez & Lecampion (2023,
10.1098/rspa.2022.0810) predict post-injection pulses reaching ~2x the shut-in
size for CRITICALLY STRESSED faults; this record shows arrest instead.

FIXED, ON STATED GROUNDS:
  f0 = 0.60          lab-constrained. Granite gouge is 0.69-0.74 (Zhang et al.
                     2022), so 0.60 already assumes weakening, and
                     CONCLUSIONS.md retracted f0 = 0.433 on exactly these
                     grounds. NOT swept.
  sigmainit = 27.99  because tau_0 = muinit*sigmainit, so varying sigmainit
                     cannot be done one-key-at-a-time.
  a, b, dc           Stage 11 showed `a` (not a-b) drives amplitude, but only
                     ~2x and at the cost of the front.
  ds 10 m, imax 601  set by the domain constraint -- see below.

SWEPT: muinit alone, hence tau_0. TAU_0 IS NOW A RESULT, NOT AN ASSUMPTION. The
old project fixed it at the measured 10.36 MPa; here it is swept exactly as Wang
& Dunham did, and whether the post-4.3 data needs 10.36 or 15.0 is the finding.

A RESULT TO KNOW BEFORE READING THE SWEEP. With f0 and sigmabar_0 both fixed,
dp_crit = sigmabar_0*(1 - muinit/f0), so matching the MEASURED 0.4-2.6 MPa
activation pressure (Habanero FDP sec 4 / Holl 2015 sec 8.4) requires
muinit = 0.544-0.591, i.e. tau_0 = 15.2-16.5 MPa -- AT OR ABOVE Taiyi's 15.0.
Hard ceiling at muinit = 0.60, where tau_0 = f0*sigmabar_0 and the fault fails
at zero overpressure. So with lab friction the activation pressure argues FOR a
critically stressed fault. That is a cross-check, not the target.

DOMAIN IS WHAT SETS ds. A weaker fault runs further, and at ds 5 m this sweep
breaches by muinit 0.42 (L/half 0.83). ds 10 m covers muinit up to ~0.54, which
is exactly the measured-stress-to-Taiyi range. Going past 0.54 needs ds 20 m and
is a Round-2 decision.

Parent is res632925.in (ds 10 m, f0 0.60, disc 150 m) because that mesh is
already characterised: 632925 gives front 0.69 / wellhead +70.5% / slip 0.99x
against res632913's 0.70 / +70.7% / 1.03x at ds 5 m, so the ds change is small.

Usage:  python build_cycle2.py [--write]
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import exp1
from scipy.optimize import brentq

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
PARENT = 632925
QFILE = "june_clean_from_d4300.txt"
T0_D, DAYS, DISC_R = 4.300, 13.1, 300.0
MUINIT = [(632950, 0.3700), (632951, 0.4120), (632952, 0.4540),
          (632953, 0.4950), (632954, 0.5360)]
A_FRONT, ALPHA, SAFE = 3.56, 1.248, 0.8


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


def front(dtauc):
    f = lambda r: A_FRONT * exp1(r ** 2 / (4 * ALPHA * DAYS * 86400)) - dtauc
    try:
        return brentq(f, 1e-3, 5e4)
    except ValueError:
        return float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()
    assert (IN / QFILE).exists(), f"{QFILE} missing"

    pairs = read_deck(IN / f"res{PARENT}.in")
    base = dict(pairs)
    sig, f0 = ff(base["sigmainit"]), ff(base["f0"])
    ds = ff(base["ds"]) * 1000.0
    IMAX = int(base["imax"]); half = IMAX * ds / 2
    kx, km = ff(base["kpmax"]), ff(base["kpmin"])
    D = ff(base["kp"]) / (ff(base["eta"]) * ff(base["phi"]) * ff(base["beta"]))
    Ld = math.sqrt(4 * D * DAYS * 86400)
    mapname = f"perm_2zone_{IMAX}_ds{int(ds)}_disc{int(DISC_R)}_kmax{kx:.1e}.txt"

    assert abs(f0 - 0.60) < 1e-9, f"parent f0 is {f0}, expected 0.60 (FIXED)"
    print(f"parent res{PARENT}.in: ds {ds:.0f} m, imax {IMAX}, f0 {f0:.2f} FIXED, "
          f"sigmainit {sig:.2f} FIXED")
    print(f"  t0 {T0_D} d, injection {QFILE}, tmax {DAYS} d, disc {DISC_R:.0f} m")
    print(f"  ceiling: tau_0 < f0*sigmabar_0 = {f0*sig:.2f} MPa, muinit < {f0:.2f}")
    print(f"  measured activation pressure 0.4-2.6 MPa needs muinit "
          f"{(sig-2.6)*f0/sig:.3f}-{(sig-0.4)*f0/sig:.3f}\n")
    print(f"{'run':>7s} {'muinit':>7s} {'tau_0':>7s} {'dp_crit':>8s} {'Dtauc':>7s} "
          f"{'R est':>7s} {'L/half':>7s}  note")
    rows, bad = [], 0
    for n, mu in MUINIT:
        tau = mu * sig; dpc = sig - tau / f0; dt = f0 * sig - tau
        R = front(dt) + (DISC_R - 150.0)
        ratio = (R + Ld) / half
        ok = ratio < SAFE
        bad += not ok
        note = ("the measured stress" if abs(mu - 0.370) < 1e-9
                else "Wang & Dunham's ratio" if abs(mu - 0.536) < 1e-9 else "")
        rows.append((n, mu, tau, dpc, dt, R, ratio))
        print(f"{n:>7d} {mu:>7.4f} {tau:>6.2f}M {dpc:>7.2f}M {dt:>6.2f}M "
              f"{R:>6.0f}m {ratio:>7.2f}  {note}"
              + ("" if ok else "   <-- BREACH"))
    print(f"\n  far-field L({DAYS} d) = {Ld:.0f} m; usable front <= "
          f"{SAFE*half-Ld:.0f} m at ds {ds:.0f} m, half {half:.0f} m")
    if bad:
        sys.exit(f"\n{bad} deck(s) breach L/half = {SAFE}. Nothing written.")
    print(f"  all {len(rows)} within L/half < {SAFE}")

    if not a.write:
        print("\ndry run -- nothing written. re-run with --write")
        return

    if not (IN / mapname).exists():
        c = (IMAX - 1) // 2
        rr = np.hypot(*(np.mgrid[0:IMAX, 0:IMAX] - c)) * ds
        kp = np.where(rr.ravel() <= DISC_R + 1e-9, kx, km)
        (IN / mapname).write_text("kp\n" + "\n".join(f"{v:.6e}" for v in kp) + "\n")
        print(f"\nwrote {mapname}  ({int((kp == kx).sum())} disc cells)")
    ndisc = int((np.loadtxt(IN / mapname, skiprows=1) == kx).sum())

    made = []
    for n, mu, tau, dpc, dt, R, ratio in rows:
        hdr = [
            f"! CYCLE 2 AS A NEW PROJECT -- t = 0 IS DATA-DAY {T0_D}.",
            f"! Round 1 of the muinit sweep. Built from res{PARENT}.in.",
            "!",
            "! 86 runs failed to reproduce the front, wellhead and slip together over",
            "! the full 17 d record. This treats the post-shut-in period as a NEW",
            "! problem: cycle 1 is DISCARDED, not inherited, and the model's own",
            "! parameters are re-swept against the post-4.3 data.",
            "!",
            "! WHY THIS WINDOW IS THE BETTER TARGET. Cycle 1 is all transient -- two",
            "! wellhead spikes at 43.15 and 44.73 MPa, a bled-off gauge from 1.59 to",
            "! 3.56 d, a seismicity gap. Post-4.3 is a smooth ramp from 32.6 MPa to",
            "! ~43.7 by 8 d, a FOUR-DAY QUASI-STEADY PLATEAU at 43.7 MPa (8-12 d),",
            "! then a step to ~48-52. A plateau constrains injectivity; spikes cannot.",
            "! It is NOT the gentler half: it sustains HIGHER pressure than cycle 1,",
            "! plus a 2-sample 87.76 MPa transient at 14.23 d (92% of Sv, almost",
            "! certainly water-hammer -- mask it, do not fit it).",
            "!",
            "! CUTTING HERE IS PHYSICALLY JUSTIFIED. Through the shut-in the catalogue",
            "! drops from 615 to 9.6 events/day and the front does NOT advance (318 m",
            "! during, vs 418 m already reached). Saez & Lecampion 2023",
            "! (10.1098/rspa.2022.0810) predict post-injection pulses reaching ~2x the",
            "! shut-in size for CRITICALLY STRESSED faults; this record shows arrest.",
            "!",
            "! FIXED, NOT SWEPT:",
            f"!   f0 = {f0:.2f} -- lab-constrained. Granite gouge is 0.69-0.74 (Zhang et al.",
            "!     2022), so 0.60 already assumes weakening; CONCLUSIONS.md retracted",
            "!     f0 = 0.433 on exactly these grounds.",
            f"!   sigmainit = {sig:.2f} -- because tau_0 = muinit*sigmainit, varying it",
            "!     cannot be done one-key-at-a-time.",
            "!   a, b, dc -- Stage 11 showed `a` drives amplitude only ~2x and at the",
            "!     cost of the front.",
            "!",
            f"! SWEPT: muinit alone. THIS RUN muinit = {mu:.4f}, so",
            f"! tau_0 = {tau:.2f} MPa, dp_crit = sigmabar_0(1 - muinit/f0) = {dpc:.2f} MPa,",
            f"! Dtauc = f0*sigmabar_0 - tau_0 = {dt:.2f} MPa.",
            "! TAU_0 IS A RESULT HERE, NOT AN ASSUMPTION -- the old project fixed it at",
            "! the measured 10.36 MPa; this sweeps it as Wang & Dunham did.",
            "!",
            f"! Ceiling: tau_0 < f0*sigmabar_0 = {f0*sig:.2f} MPa, i.e. muinit < {f0:.2f}.",
            f"! The MEASURED 0.4-2.6 MPa activation pressure needs muinit",
            f"! {(sig-2.6)*f0/sig:.3f}-{(sig-0.4)*f0/sig:.3f}, i.e. tau_0 15.2-16.5 MPa -- at or ABOVE",
            "! Taiyi's 15.0. With lab friction the activation pressure argues FOR a",
            "! critically stressed fault. Cross-check, not target.",
            "!",
            f"! INITIAL DISC IS A FREE PARAMETER, here {DISC_R:.0f} m ({ndisc} cells) so the",
            "! simulated front starts near the observed 313 m at t0. Its radius is",
            "! Round 2's variable.",
            "!",
            f"! INJECTION {QFILE}: june_clean.txt shifted by -{T0_D} d,",
            "! 285 points, volume conserved to 0.00000% against the record integrated",
            "! from t0. june_clean.txt UNTOUCHED (md5",
            "! 603285de6a5b1145d61bea05d97ebbf7).",
            "!",
            "! SCORING: score_cycle2.py, NOT score_grid.py. Observed series read at",
            f"! t_obs = t_sim + {T0_D} d; slip compared against the observed INCREMENT",
            "! (obs - 2.81 cm, target 6.36 cm by data-day 17) because main_LH.f90:923",
            "! sets slip = 0d0; front compared as R(t) DIRECTLY with no lambda*sqrt(t)",
            "! fit, since neither front starts at the origin.",
            "!",
            f"! DOMAIN: ds {ds:.0f} m is forced -- at ds 5 m this sweep breaches by",
            f"! muinit 0.42 (L/half 0.83). Front estimate {R:.0f} m at {DAYS} d plus",
            f"! {Ld:.0f} m far-field diffusion against a {half:.0f} m half-domain,",
            f"! L/half = {ratio:.2f}, usable limit {SAFE*half-Ld:.0f} m.",
        ]
        ov = {"muinit": f"{mu:.4f}", "injection_file": f'"{QFILE}"',
              "parameter_file": f'"{mapname}"',
              "tmax": f"{DAYS/365.0:.8f}", "nstep": "400000"}
        lines = [f"filenumber {n}"] + hdr
        for k, v in pairs:
            if k == "filenumber":
                continue
            lines.append(f"{k} {ov.get(k, v)}")
        (IN / f"res{n}.in").write_text("\n".join(lines) + "\n")
        made.append((n, mu))

    print("\nverifying")
    nbad = 0
    ref = dict(read_deck(IN / "res632950.in"))
    for n, mu in made:
        d = dict(read_deck(IN / f"res{n}.in"))
        mp = np.loadtxt(IN / mapname, skiprows=1)
        # vs the FIRST sibling: only filenumber and muinit may differ
        dsib = {k for k in set(d) | set(ref) if d.get(k) != ref.get(k)}
        # vs the parent: the four intended keys plus muinit
        dpar = {k for k in set(d) | set(base) if d.get(k) != base.get(k)}
        ck = {
            "vs sibling: only filenumber/muinit":
                not (dsib - {"filenumber", "muinit"}),
            "vs parent: only the intended keys":
                not (dpar - {"filenumber", "muinit", "injection_file",
                             "parameter_file", "tmax", "nstep"}),
            "muinit": abs(ff(d["muinit"]) - mu) < 1e-9,
            "f0 == 0.60 FIXED": abs(ff(d["f0"]) - 0.60) < 1e-9,
            "sigmainit FIXED": d["sigmainit"] == base["sigmainit"],
            "a,b,dc FIXED": all(d[k] == base[k] for k in ("a", "b", "dc")),
            "tau_0 below ceiling": ff(d["muinit"]) < 0.60,
            "pfinit == 0": ff(d["pfinit"]) == 0.0,
            "shifted injection": d["injection_file"].strip('"') == QFILE,
            "kpmax == map max": abs(mp.max() - kx) / kx < 1e-6,
            "kpmin == map min": abs(mp.min() - km) / km < 1e-6,
            "map size == imax^2": mp.size == IMAX ** 2,
            "ds/imax unchanged": (d["ds"] == base["ds"]
                                  and d["imax"] == base["imax"]),
            "tmax": abs(ff(d["tmax"]) * 365 - DAYS) < 1e-4,
        }
        ok = all(ck.values()); nbad += not ok
        print(f"  res{n}.in muinit {mu:.4f} (tau_0 {mu*sig:>5.2f} MPa): "
              + ("OK" if ok else "PROBLEM "
                 + str([k for k, v in ck.items() if not v])))
    print(f"\n  {len(made)-nbad}/{len(made)} clean")
    print(f"  each sibling vs res632950.in differs only in filenumber + muinit")
    sys.exit(0 if nbad == 0 else 1)


if __name__ == "__main__":
    main()
