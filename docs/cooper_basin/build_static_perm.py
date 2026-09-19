#!/usr/bin/env python3
"""CYCLE 2, Round 7: turn permeability EVOLUTION off and let the background
diffusivity carry the front. One run, off res633001.in.

WHY, AND WHY IT IS NOT A TUNING MOVE

633001 matches the wellhead plateau to +2.4% and then misses everything after
day 9. The miss is not an offset, it is a slope: across the day-9 rate step the
observed dp rises +4.64 MPa and the model +1.68, so the model's response to a
rate increase is 2.75x too small. Sharper still, between 9.6 and 12.5 d the
injection rate FALLS from 44.5 to 37.8 L/s while the observed dp RISES
14.63 -> 16.29 and the model's FALLS 14.80 -> 13.44.

Measured on sustained-flow windows, the well is LINEAR: dp/q is 439 kPa per L/s
at 4-8 d and 431 at 10.2-12.5 d across a 1.4x change in rate. The model is 443
then 355 -- it is the model whose injectivity changes, by 20%, not the well's.

THE CAUSE IS THE MECHANISM THAT MAKES THE FRONT. Inside the disc kp already
starts AT kpmax, so there is no local enhancement left to happen; the model
softens because the enhanced ZONE GROWS OUTWARD as slip spreads. That growing
zone is also what propagates the front. So sweeping kL cannot separate the two
-- it slows the front by exactly the factor it slows the softening. That
objection is why this run changes the structure instead.

A CLOSED-SYSTEM EXPLANATION WAS TESTED AND REJECTED. If a bounded compartment
were filling, dp would go linear in cumulative volume (pseudo-steady state).
Observed d(dp)/dV is 0.015 MPa/ML over 4-8 d, 1.13 over 8-10 d and -0.24 over
10-12.6 d, and a straight line through the late half leaves a 1.49 MPa rms
residual. Not pseudo-steady state.

WHAT THIS RUN DOES INSTEAD. Constant injectivity and a propagating front are
not in conflict -- that is ordinary Theis behaviour. Take the front from a
STATIC background diffusivity rather than from a growing zone, and nothing about
the well changes with time, so the linearity the data show comes out by
construction rather than by fitting.

    front:        the observed percentile front reaches R = 763.1 m at 8.7 d.
                  Shapiro's triggering front r = sqrt(4 pi D t) needs
                  D = R^2/(4 pi t) = 0.0616 m2/s.
    that D needs  kp = D * eta * phi * beta = 0.0616 * 2.8575e-14
                     = 1.76e-15 m^2,
                  against the present background of 1.0e-15 (D 0.0350). The
                  background is 1.8x too low, which is the whole reason the
                  front has to be manufactured by enhancement.
    injectivity:  unchanged -- the disc stays at kpmax 2.5e-13, 142x the new
                  background. That gap is real and is what justifies having two
                  zones at all.

THE DISC IS NOT "ENHANCEMENT REMOVED", IT IS ENHANCEMENT ALREADY DONE. The
simulation starts at data-day 4.300, after cycle 1 injected ~8 ML and slipped
the fault. The near-well permeability was enhanced during cycle 1; by t = 0 that
is finished. A static disc IS that enhanced zone, which is the same logic as the
cycle-2 reframing and makes the disc a physical initial condition rather than a
fitted knob.

THE BACKGROUND MUST GO IN THE MAP, NOT THE kpmin KEY. main_LH.f90:2413 is the
only place kpmin is read, inside the evolution ODE, so with permev F the kpmin
key is INERT. The per-cell map is applied at main_LH.f90:864-886, which is
outside every permev gate (the gates at :749, :802, :961 and :1419 are kp
restart/output I/O only), so the disc survives with evolution off. Hence a new
map file. kpmin and kp are set to the new background anyway so that the deck
cannot disagree with the map about what the background is -- one physical
number, written in three places, asserted equal below.

kL and kT are left at 1e-3 and 1e15. Both are inert with permev F; removing
them would make the diff against res633001.in read as four changes rather than
the two that matter.

WHAT WOULD FALSIFY THIS. HBI puts the front where SLIP crosses a threshold, not
where pressure does, so the sqrt(4 pi D t) inversion above is approximate and
the static front may land short of 763 m. If lambda collapses, permeability
evolution is load-bearing for the front and that is the finding; 633001 remains
the reference either way.

Usage:
  python build_static_perm.py            # dry run, prints the diff and checks
  python build_static_perm.py --write
"""
import argparse
import sys
from pathlib import Path

import numpy as np

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
PARENT = 633001
NEW_DEFAULT = 633120
DAYS, SAFE = 13.1, 0.8
R_OBS, T_OBS = 763.1, 8.7          # percentile front, metres, at sim days


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
    # ROUND 7b. 633120 came out +28.3% on pressure and 0.79x on radius, so the
    # sqrt(4 pi D t) inversion below UNDERESTIMATED the background -- exactly
    # the failure flagged in the header, since HBI puts the front where slip
    # crosses a threshold, not where pressure does. Raising the background
    # lowers dp (~1/k) and lengthens the front (~sqrt k), so ONE knob moves both
    # errors the right way: k x 1.283 predicts dp ~0% and R/R_obs 0.89x.
    ap.add_argument("--background", type=float, default=None,
                    help="override the inverted background, m^2")
    ap.add_argument("--filenumber", type=int, default=None)
    a = ap.parse_args(argv)

    pairs = read_deck(IN / f"res{PARENT}.in")
    base = dict(pairs)
    ds = ff(base["ds"]) * 1000.0
    IMAX = int(base["imax"])
    half = IMAX * ds / 2
    den = ff(base["eta"]) * ff(base["phi"]) * ff(base["beta"])
    kx = ff(base["kpmax"])
    km_old = ff(base["kpmin"])
    disc = 350

    # --- invert the observed front for the background it implies
    D_new = R_OBS ** 2 / (4 * np.pi * T_OBS * 86400.0)
    km_new = float(f"{a.background if a.background else D_new * den:.3g}")
    D_chk = km_new / den

    print(f"parent res{PARENT}.in: ds {ds:.0f} m, imax {IMAX}, "
          f"half-domain {half:.0f} m")
    print(f"  eta*phi*beta = {den:.4e}")
    print(f"  observed front {R_OBS:.1f} m at {T_OBS} d  ->  "
          f"D = R^2/(4 pi t) = {D_new:.4f} m2/s")
    print(f"  background kp: {km_old:.2e} (D {km_old/den:.4f}) -> "
          f"{km_new:.2e} (D {D_chk:.4f}),  {km_new/km_old:.2f}x")
    print(f"  disc unchanged at kpmax {kx:.2e}, contrast "
          f"{kx/km_new:.0f}x (was {kx/km_old:.0f}x)")

    # --- DOMAIN. With evolution off the pressure front is the static Shapiro
    # front through the new background, plus the disc radius it crosses first.
    # This is the CONSERVATIVE bound: the real slip front cannot outrun the
    # pressure that drives it.
    r_press = np.sqrt(4 * np.pi * D_chk * DAYS * 86400.0)
    reach = disc + r_press
    ratio = reach / half
    print(f"\nDOMAIN at tmax {DAYS} d")
    print(f"  static pressure front sqrt(4 pi D t) = {r_press:.0f} m, "
          f"+ disc {disc} m  ->  {reach:.0f} m")
    print(f"  L/half = {ratio:.2f}" + ("  OK" if ratio < SAFE else "  BREACH"))
    if ratio >= SAFE:
        sys.exit(f"L/half {ratio:.2f} >= {SAFE}. Nothing written.")

    mn = (f"perm_2zone_{IMAX}_ds{int(ds)}_disc{disc}_kmax{kx:.2e}"
          f"_kmin{km_new:.2e}.txt")
    NEW = a.filenumber or NEW_DEFAULT
    over = {"filenumber": str(NEW),
            "permev": "F",
            "kpmin": f"{km_new:.2e}",
            "kp": f"{km_new:.2e}",
            "parameter_file": f'"{mn}"'}

    # --- the diff must be exactly the five keys above, and no others
    changed = [(k, base[k], over[k]) for k in over if base.get(k) != over[k]]
    print(f"\nDIFF against res{PARENT}.in -- {len(changed)} keys")
    for k, o, n in changed:
        print(f"  {k:16s} {o:44s} -> {n}")
    assert set(k for k, _, _ in changed) == set(over), \
        "a key in the override dict did not actually change"
    assert len(over) == 5, "expected exactly 5 keys to move"

    # --- the background is written in three places; they must agree
    assert abs(ff(over["kpmin"]) - km_new) < 1e-20
    assert abs(ff(over["kp"]) - km_new) < 1e-20
    assert km_new > km_old, "this run RAISES the background"
    assert kx > 100 * km_new, "the disc must stay well above the background"

    hdr = [
        f"! CYCLE 2 ROUND 7: PERMEABILITY EVOLUTION OFF, background raised to",
        f"! carry the front. ONE STRUCTURAL CHANGE off res{PARENT}.in.",
        "!",
        f"! 633001 matches the wellhead plateau to +2.4% and then misses",
        "! everything after day 9. The miss is a SLOPE, not an offset: across the",
        "! day-9 rate step the observed dp rises +4.64 MPa and the model +1.68.",
        "! Between 9.6 and 12.5 d the rate FALLS 44.5 -> 37.8 L/s while the",
        "! observed dp RISES 14.63 -> 16.29 and the model's FALLS 14.80 -> 13.44.",
        "!",
        "! MEASURED, THE WELL IS LINEAR. dp/q is 439 kPa per L/s at 4-8 d and 431",
        "! at 10.2-12.5 d, across a 1.4x change in rate. The model is 443 then",
        "! 355 -- the model's injectivity improves by 20%, the well's does not.",
        "!",
        "! THE CAUSE IS THE MECHANISM THAT MAKES THE FRONT. Inside the disc kp",
        "! already starts AT kpmax, so no local enhancement remains; the model",
        "! softens because the enhanced ZONE GROWS as slip spreads, and that same",
        "! growth propagates the front. Sweeping kL therefore cannot separate",
        "! them -- it slows the front by the factor it slows the softening.",
        "!",
        "! A CLOSED COMPARTMENT WAS TESTED AND REJECTED: observed d(dp)/dV is",
        "! 0.015 MPa/ML over 4-8 d, 1.13 over 8-10 d, -0.24 over 10-12.6 d, and a",
        "! line through the late half leaves 1.49 MPa rms. Not pseudo-steady.",
        "!",
        "! WHAT THIS RUN DOES. Constant injectivity and a propagating front are",
        "! not in conflict -- that is ordinary Theis behaviour. Take the front",
        "! from a STATIC background instead of a growing zone and nothing about",
        "! the well changes with time, so the observed linearity follows by",
        "! construction rather than by fitting.",
        "!",
        f"!   the observed percentile front is R = {R_OBS:.1f} m at {T_OBS} d;",
        f"!   Shapiro r = sqrt(4 pi D t) needs D = R^2/(4 pi t) = {D_new:.4f} m2/s;",
        f"!   D = kp/(eta phi beta) with eta phi beta = {den:.4e} needs",
        f"!   kp = {km_new:.2e} m^2, against the present {km_old:.1e}",
        f"!   (D {km_old/den:.4f}). The background is {km_new/km_old:.2f}x too low,",
        "!   which is the whole reason the front had to be manufactured.",
        "!",
        f"! THE DISC IS UNCHANGED at kpmax {kx:.2e}, now {kx/km_new:.0f}x the",
        f"! background (was {kx/km_old:.0f}x). The well's injectivity needs that",
        "! gap; it is what justifies two zones at all.",
        "!",
        "! THE DISC IS NOT 'ENHANCEMENT REMOVED', IT IS ENHANCEMENT ALREADY DONE.",
        "! This run starts at data-day 4.300, after cycle 1 injected ~8 ML and",
        "! slipped the fault. The near-well permeability was enhanced THEN. A",
        "! static disc is that enhanced zone, which makes it a physical initial",
        "! condition rather than a fitted knob.",
        "!",
        "! THE BACKGROUND GOES IN THE MAP, NOT THE kpmin KEY. main_LH.f90:2413 is",
        "! the only place kpmin is read, inside the evolution ODE, so with",
        "! permev F the kpmin key is INERT. The per-cell map is applied at",
        "! main_LH.f90:864-886, outside every permev gate (:749, :802, :961 and",
        "! :1419 are kp restart/output I/O only), so the disc survives with",
        "! evolution off. kpmin and kp are set to the new background anyway so the",
        "! deck cannot disagree with the map; the builder asserts all three.",
        "!",
        "! kL 1e-3 and kT 1e15 are left in place. Both are inert with permev F.",
        "!",
        f"! DOMAIN: static front sqrt(4 pi D t) = {r_press:.0f} m at {DAYS} d plus",
        f"! the {disc} m disc gives {reach:.0f} m against a {half:.0f} m",
        f"! half-domain, L/half {ratio:.2f}. The slip front cannot outrun the",
        "! pressure that drives it, so this bound is conservative.",
        "!",
        "! WHAT WOULD FALSIFY IT. HBI puts the front where SLIP crosses a",
        "! threshold, not where pressure does, so the inversion above is",
        "! approximate and the static front may land short. If lambda collapses,",
        "! permeability evolution is load-bearing for the front and that is the",
        f"! result; res{PARENT}.in remains the reference either way.",
        "!",
    ]

    if not a.write:
        print("\ndry run. re-run with --write")
        return

    # --- the map: disc at kpmax, everything else at the new background
    c_ = (IMAX - 1) // 2
    rr = np.hypot(*(np.mgrid[0:IMAX, 0:IMAX] - c_)) * ds
    kp = np.where(rr.ravel() <= disc + 1e-9, kx, km_new)
    if not (IN / mn).exists():
        (IN / mn).write_text("kp\n" + "\n".join(f"{v:.6e}" for v in kp) + "\n")
        print(f"\nwrote {mn}")
    chk = np.loadtxt(IN / mn, skiprows=1)
    assert chk.size == IMAX ** 2, f"map has {chk.size} cells, expected {IMAX**2}"
    assert abs(chk.max() - kx) < 1e-20, "map max != kpmax"
    assert abs(chk.min() - km_new) < 1e-20, "map min != the new background"
    frac = (chk == chk.max()).mean()
    exp = np.pi * disc ** 2 / (IMAX * ds) ** 2
    assert abs(frac - exp) < 0.1 * exp, \
        f"disc covers {frac:.4%} of cells, expected ~{exp:.4%}"
    print(f"  map verified: {chk.size} cells, {frac:.4%} at kpmax "
          f"(pi r^2 predicts {exp:.4%})")

    # FILENUMBER MUST BE THE FIRST LINE, for two independent reasons, and
    # writing hdr first violated both.
    #
    # 1. march26_submit_hbi_git_scratch.sh:85 does
    #    FILE_NUM=$(awk 'NR==1 {print $2}'), so a comment on line 1 made
    #    FILE_NUM="CYCLE" and the post-run rsync globbed *CYCLE*.dat, moving
    #    nothing.
    # 2. WORSE, and the real bug: main_LH.f90:2674 is
    #    read(33,*,iostat=ios) param,pvalue -- LIST-DIRECTED, so it needs TWO
    #    values and CONTINUES INTO THE NEXT RECORD to get them. A bare "!"
    #    header line supplies only one, so the read consumed "filenumber" as
    #    its pvalue and discarded 633120. `number` kept its default 0 and the
    #    run wrote slip0.dat, pf0.dat, ... Every other key parsed fine, since
    #    the next read resumed cleanly, so the run was scientifically valid and
    #    only misnamed -- but silently, which is the dangerous part.
    #
    # A scan of all ~100 res6*.in decks found this deck was the only one
    # affected: the others put filenumber on line 1 and none has a bare "!"
    # immediately preceding a key.
    out = [f"{k} {over.get(k, v)}" for k, v in pairs]
    assert out[0].startswith("filenumber "), \
        f"filenumber must be the first deck line, got {out[0]!r}"
    # Moving filenumber to line 1 protects filenumber and nothing else: hdr
    # still ENDED in a bare "!", which now sits immediately before out[1] and
    # would swallow 'problem' by the same list-directed read. Drop any trailing
    # bare comment so no key is ever the record after a lone "!".
    while hdr and hdr[-1].strip() == "!":
        hdr.pop()
    body = out[:1] + hdr + out[1:]
    for i, (a_, b_) in enumerate(zip(body, body[1:])):
        assert not (a_.strip() == "!" and not b_.startswith("!")), (
            f"line {i+1} is a bare '!' and line {i+2} is the key "
            f"{b_.split()[0]!r}; Fortran's list-directed read will swallow it")
    dp = IN / f"res{NEW}.in"
    assert not dp.exists(), f"{dp} already exists -- pick a new filenumber"
    dp.write_text("\n".join(body) + "\n")
    print(f"wrote {dp}")

    # --- read it back and re-diff, so what is on disk is what was checked
    rb = dict(read_deck(dp))
    d2 = [k for k in set(rb) | set(base) if rb.get(k) != base.get(k)]
    assert sorted(d2) == sorted(over), f"on-disk diff is {sorted(d2)}"
    assert rb["parameter_file"].strip('"') == mn
    assert rb["permev"] == "F"
    print(f"  on-disk diff against res{PARENT}.in is exactly {sorted(d2)}")


if __name__ == "__main__":
    main()
