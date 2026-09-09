#!/usr/bin/env python3
"""CYCLE 2, Round 2: break the front/dp trade-off at the MEASURED stress.

WHAT ROUND 1 ESTABLISHED. Twenty runs, four fluid parameterisations x five
tau_0. Arm 2 (eta 1.27e-4, beta 2.25e-8) is the only live one, and within it:

    tau_0    front@8.7d    dp error
    10.36       0.97         +55%          <- front matches, dp does not
    11.53       1.08         +37%
    12.71       1.26         +18%
    13.86       1.53         -0.2%         <- dp matches, front does not
    15.00       1.99         -18%

The front and the wellhead dp respond to tau_0 with OPPOSITE SIGN, so no tau_0
fits both; the best compromise is ~tau_0 12.1 at front 1.17 / dp +28%. The same
opposition holds in all four arms.

WHY, AND THIS IS WHAT ROUND 2 ACTS ON. Arm 3 gives the mechanism away: its dp
barely moves across the whole tau_0 sweep (-25.8% to -29.9%) because almost
nothing slips, while arms 1, 2 and 4 swing 70-90 points over the same range. So
the coupling runs through permev -- slip raises kp toward kpmax, injectivity
rises, wellhead dp falls, AND the front grows, all through one pathway. Sweeping
anything that changes how much slips therefore trades one target against the
other by construction.

The way out is to raise near-well injectivity WITHOUT requiring slip to do it,
at tau_0 = 10.36 where the front already matches at 0.97. Two levers, tested
separately:

  LEVER 1, skin (632990-632992). Enters ONLY the Peaceman well-cell
    transmissivity, m_diffusion.f90:669:

        T = 2 pi kp / eta / (log(0.2 ds / rw) + skin)

    With ds = 10 m and rw = 0.089 m, log(0.2*10/0.089) = 3.112, so skin shifts
    the well-cell coupling and nothing else -- the far-field diffusivity, the
    permeability map and the stress state are untouched, so the front should
    barely move. That is the prediction being tested.

    Habanero 4 is a hydraulically stimulated injector and NEGATIVE SKIN IS WHAT
    A STIMULATED WELL HAS. skin = 0 was the code default (main_LH.f90:179),
    never a measurement.

    HOW MUCH IS NEEDED, and why it is only an estimate: if the whole dp were
    the Peaceman drop, killing +55% would need (3.112+skin)/3.112 = 1/1.553,
    i.e. skin = -1.11. It is not -- dp also contains the radial drop through
    the 300 m disc and, more importantly, through the kpmin = 1e-15 background
    beyond it, which is 250x less permeable. So -1.11 is an upper bound on the
    effect and the response has to be measured. Hence three values.

    HARD BOUND: skin > -3.112, or the denominator changes sign and T is
    unphysical. -2.0 leaves 1.112, i.e. T x 2.80, and is deliberately the most
    aggressive value tried.

  LEVER 2, the initial high-k disc radius (632993-632994), 300 m -> 450, 600 m.
    Replaces kpmin background with kpmax over a wider annulus. Whether this
    helps dp at all is genuinely open:

      - AGAINST: the radial pressure drop is logarithmic, so the annulus
        300-600 m contributes log(2) = 0.69 out of log(300/0.089) = 8.1 of the
        drop inside the disc -- about 8%.
      - FOR: the drop through the kpmin background beyond the disc is not
        logarithmically small, because kpmin is 250x below kpmax. Pushing that
        boundary outward removes the expensive part.

    UNLIKE skin, THIS LEVER IS NOT CLEAN. Enlarging the high-k region raises
    the effective diffusivity out to that radius, so the front will grow -- the
    front is already 0.97, so growth is the wrong direction. The test is
    quantitative: if dp falls a lot and the front grows a little, the lever is
    still useful; if both move together, it is just another point on the same
    trade-off and the disc is not the answer. Predicted fronts are asserted
    against the domain below either way.

    600 m IS THE LARGEST THE DOMAIN ALLOWS, and this was not a judgement call.
    A 900 m disc was in the first draft and the domain check REFUSED TO WRITE
    IT: on the pessimistic estimate below it reaches L/half = 0.96 against the
    0.8 limit, i.e. into the boundary. The alternatives -- ds 20 m for a 6010 m
    half-domain, or a shorter tmax -- each change a second key and would stop
    this being a one-parameter test, so the disc is capped instead. If the
    lever earns it, 900 m can follow on a wider domain, sized from a MEASURED
    front rather than from this estimate.

EVERYTHING ELSE IS 632960's, unchanged, and each deck differs from it in
exactly ONE parameter -- verified by diff at the end, not by assertion in the
prose. tau_0 stays at the measured 10.36 MPa in all five.

Usage:
  python build_cycle2_round2.py            # dry run: the table and the checks
  python build_cycle2_round2.py --write    # write decks and maps
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
PARENT = 632960
DAYS, SAFE = 13.1, 0.8
SKINS = [(632990, -0.5), (632991, -1.0), (632992, -2.0)]
DISCS = [(632993, 450.0), (632994, 600.0)]
# 632960 measured: front 831 m at 13.1 d with a 300 m disc. The disc runs are
# scaled from that rather than from the E1 estimator, because the estimator was
# calibrated at eta = 0.89e-3 and over-predicts these fronts by ~100x.
R_MEASURED_300 = 831.0


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

    pairs = read_deck(IN / f"res{PARENT}.in")
    base = dict(pairs)
    ds = ff(base["ds"]) * 1000.0
    IMAX = int(base["imax"]); half = IMAX * ds / 2
    kx, km = ff(base["kpmax"]), ff(base["kpmin"])
    rw = ff(base["rw"])
    sig, f0, mu = ff(base["sigmainit"]), ff(base["f0"]), ff(base["muinit"])
    D_far = ff(base["kp"]) / (ff(base["eta"]) * ff(base["phi"]) * ff(base["beta"]))
    Ld = math.sqrt(4 * D_far * DAYS * 86400)
    peac = math.log(0.2 * ds / rw)

    assert abs(f0 - 0.60) < 1e-9, f"parent f0 is {f0}, expected 0.60 (FIXED)"
    assert abs(mu * sig - 10.36) < 0.02, \
        f"parent tau_0 is {mu*sig:.2f}, expected the measured 10.36 MPa"
    assert "skin" not in base, \
        "parent already sets skin -- it would no longer be a one-key change"

    print(f"parent res{PARENT}.in: arm 2, tau_0 = {mu*sig:.2f} MPa (MEASURED), "
          f"ds {ds:.0f} m, imax {IMAX}, disc 300 m")
    print(f"  front measured 831 m at {DAYS} d, dp error +55.3%")
    print(f"  Peaceman log(0.2*ds/rw) = log(0.2*{ds:.0f}/{rw:g}) = {peac:.3f}, "
          f"so skin > {-peac:.3f} is required")
    print(f"  far-field L({DAYS} d) = {Ld:.0f} m, half-domain {half:.0f} m, "
          f"usable front <= {SAFE*half-Ld:.0f} m\n")

    rows, bad = [], 0
    print(f"{'run':>7s} {'lever':>7s} {'value':>7s} {'T factor':>9s} "
          f"{'R pred':>7s} {'L/half':>7s}")
    for n, sk in SKINS:
        fac = peac / (peac + sk)
        # skin does not touch the permeability field or the stress state, so
        # the predicted front is the measured one. If it moves, the premise of
        # the lever is wrong, which is itself the result.
        R = R_MEASURED_300
        ratio = (R + Ld) / half
        bad += ratio >= SAFE
        rows.append((n, "skin", sk, fac, R, ratio))
        print(f"{n:>7d} {'skin':>7s} {sk:>7.2f} {fac:>8.2f}x {R:>6.0f}m "
              f"{ratio:>7.2f}" + ("" if ratio < SAFE else "   <-- BREACH"))
    for n, dr in DISCS:
        # Worst case: the front scales with the disc radius, since inside the
        # disc the diffusivity is uniform at kpmax/(eta phi beta). That is an
        # OVER-estimate -- the front is set by the strength margin, not only by
        # how far the high-k region reaches -- and it is used deliberately so
        # the domain check errs toward refusing to write.
        R = R_MEASURED_300 * dr / 300.0
        ratio = (R + Ld) / half
        bad += ratio >= SAFE
        rows.append((n, "disc", dr, 1.0, R, ratio))
        print(f"{n:>7d} {'disc':>7s} {dr:>6.0f}m {'-':>9s} {R:>6.0f}m "
              f"{ratio:>7.2f}" + ("" if ratio < SAFE else "   <-- BREACH"))

    if bad:
        sys.exit(f"\n{bad} deck(s) breach L/half = {SAFE}. Nothing written.")
    print(f"\n  all {len(rows)} within L/half < {SAFE} on the pessimistic "
          f"front estimate")

    if not a.write:
        print("\ndry run -- nothing written. re-run with --write")
        return

    # ------------------------------------------------------------ perm maps
    c = (IMAX - 1) // 2
    rr = np.hypot(*(np.mgrid[0:IMAX, 0:IMAX] - c)) * ds
    maps = {}
    for _, dr in DISCS:
        mn = f"perm_2zone_{IMAX}_ds{int(ds)}_disc{int(dr)}_kmax{kx:.1e}.txt"
        if not (IN / mn).exists():
            kp = np.where(rr.ravel() <= dr + 1e-9, kx, km)
            (IN / mn).write_text("kp\n" + "\n".join(f"{v:.6e}" for v in kp) + "\n")
            print(f"wrote {mn}  ({int((kp == kx).sum())} disc cells)")
        maps[dr] = mn

    # ---------------------------------------------------------------- decks
    made = []
    for n, lever, val, fac, R, ratio in rows:
        if lever == "skin":
            what = [
                f"! ROUND 2, LEVER 1: skin = {val:.2f}. ONE KEY off res{PARENT}.in.",
                "!",
                "! Round 1 left the front and the wellhead dp irreconcilable: in arm 2",
                "! the front matches at tau_0 10.36 (0.97) where dp is +55%, and dp",
                "! matches at tau_0 13.86 where the front is 1.53. They respond to",
                "! tau_0 with OPPOSITE SIGN, so no tau_0 fits both.",
                "!",
                "! Arm 3 shows why: its dp moves only -25.8 to -29.9% across the whole",
                "! tau_0 sweep because almost nothing slips, while arms 1/2/4 swing",
                "! 70-90 points. The coupling is THROUGH permev -- slip raises kp",
                "! toward kpmax, injectivity rises, dp falls, and the front grows, one",
                "! pathway for both. So anything changing how much slips trades the two",
                "! against each other by construction.",
                "!",
                "! skin does not. m_diffusion.f90:669:",
                "!     T = 2 pi kp / eta / (log(0.2 ds / rw) + skin)",
                f"! It enters ONLY the Peaceman well-cell transmissivity. Here",
                f"! log(0.2*{ds:.0f}/{rw:g}) = {peac:.3f}, so skin = {val:.2f} gives",
                f"! T x {fac:.2f}. The permeability map, the far-field diffusivity and",
                "! the stress state are untouched.",
                "!",
                "! PREDICTION ON RECORD: dp falls, the front stays near 0.97. If the",
                "! front moves appreciably the premise is wrong, and that is the result.",
                "!",
                f"! Killing +55% would need skin = -1.11 IF dp were entirely the",
                "! Peaceman drop. It is not -- dp also contains the radial drop through",
                "! the 300 m disc and through the kpmin = 1e-15 background beyond it,",
                "! 250x less permeable. -1.11 is therefore an upper bound on the",
                f"! effect; three values ({', '.join(f'{s:.1f}' for _, s in SKINS)}) map",
                "! the actual response.",
                "!",
                f"! HARD BOUND skin > {-peac:.3f}, else the denominator changes sign and",
                f"! T is unphysical. This run leaves {peac + val:.3f}.",
                "!",
                "! JUSTIFIED, not fitted: Habanero 4 is a hydraulically stimulated",
                "! injector, and negative skin is what a stimulated well has. skin = 0",
                "! is the code default at main_LH.f90:179, never a measurement.",
                "!",
                f"! tau_0 stays at the MEASURED {mu*sig:.2f} MPa. Domain: front predicted",
                f"! {R:.0f} m (= 632960's measured value, since skin should not move it)",
                f"! plus {Ld:.0f} m far-field against {half:.0f} m half-domain, "
                f"L/half {ratio:.2f}.",
            ]
            ov = {"rw": base["rw"]}          # unchanged; skin is appended below
            extra = [f"skin {val}"]
        else:
            ov = {"parameter_file": f'"{maps[val]}"'}
            extra = []
            what = [
                f"! ROUND 2, LEVER 2: initial high-k disc {val:.0f} m, from 300 m.",
                f"! ONE KEY off res{PARENT}.in (parameter_file).",
                "!",
                "! Same problem as lever 1 addresses -- see res632991.in for the full",
                "! statement. Round 1's arm 2 matches the front at tau_0 10.36 and dp",
                "! at 13.86, and the two cannot be reconciled by tau_0 because both",
                "! respond to it through the same slip -> kp -> injectivity pathway.",
                "!",
                "! This lever replaces kpmin background with kpmax over a wider",
                "! annulus. WHETHER IT HELPS dp IS GENUINELY OPEN:",
                "!   AGAINST: the radial drop is logarithmic, so 300->600 m adds",
                f"!     log(2) = 0.69 of {math.log(300/rw):.1f} inside the disc, about 8%.",
                "!   FOR: the drop through the kpmin background beyond the disc is NOT",
                "!     logarithmically small, kpmin being 250x below kpmax. Moving that",
                "!     boundary outward removes the expensive part.",
                "!",
                "! THIS LEVER IS NOT CLEAN, unlike skin. A wider high-k region raises",
                "! the effective diffusivity out to that radius, so the front will grow,",
                "! and the front is already 0.97 -- growth is the wrong direction. The",
                "! test is quantitative: a large dp drop for a small front increase",
                "! makes it useful; both moving together makes it one more point on the",
                "! same trade-off.",
                "!",
                f"! DOMAIN, on a deliberately pessimistic estimate: the front is taken to",
                f"! scale with the disc radius, {R_MEASURED_300:.0f} x {val:.0f}/300 = "
                f"{R:.0f} m, which over-predicts",
                "! because the front is set by the strength margin and not only by the",
                f"! reach of the high-k region. Plus {Ld:.0f} m far-field against",
                f"! {half:.0f} m half-domain gives L/half {ratio:.2f}. If the measured",
                f"! front exceeds the {SAFE*half-Ld:.0f} m usable limit, late times are",
                "! not quotable.",
                "!",
                f"! tau_0 stays at the MEASURED {mu*sig:.2f} MPa; kpmax and kpmin",
                "! unchanged, so only WHERE kpmax applies differs.",
            ]

        lines = [f"filenumber {n}"] + what
        for k, v in pairs:
            if k == "filenumber":
                continue
            lines.append(f"{k} {ov.get(k, v)}")
        lines += extra
        (IN / f"res{n}.in").write_text("\n".join(lines) + "\n")
        made.append((n, lever, val))

    # -------------------------------------------------------------- verify
    print("\nverifying that each deck differs from the parent in ONE key")
    nbad = 0
    ref = dict(pairs)
    for n, lever, val in made:
        d = dict(read_deck(IN / f"res{n}.in"))
        diff = sorted(set(ref) ^ set(d)) + \
            sorted(k for k in set(ref) & set(d)
                   if ref[k] != d[k] and k != "filenumber")
        want = {"skin"} if lever == "skin" else {"parameter_file"}
        ok = set(diff) == want
        nbad += not ok
        print(f"  {n}  changed: {diff or ['nothing']}  "
              + ("OK" if ok else f"EXPECTED {sorted(want)}  <-- WRONG"))
        if lever == "disc":
            mp = np.loadtxt(IN / d["parameter_file"].strip('"'), skiprows=1)
            got = float(mp.max())
            ncell = int((mp == got).sum())
            rad = math.sqrt(ncell / math.pi) * ds
            print(f"        map max {got:.3e} (kpmax {kx:.3e}), min {mp.min():.3e} "
                  f"(kpmin {km:.3e}), {ncell} cells -> r_eff {rad:.0f} m "
                  f"(target {val:.0f})")
            if abs(got - kx) > 1e-20 or abs(rad - val) > 1.5 * ds:
                nbad += 1
                print("        <-- MAP DOES NOT MATCH THE DECK")
    if nbad:
        sys.exit(f"\n{nbad} problem(s). Fix before submitting.")
    # -w MUST BE THE DECK'S OWN injection_file, read from the deck rather than
    # named here. That flag is the only thing that copies it into the rundir,
    # and HBI's open(77, iostat=) defaults to status='unknown', so a missing
    # injection file is CREATED EMPTY, returns ios = 0, and then dies with
    # "end-of-file during read, unit 77" at m_diffusion.f90:59 -- ten seconds
    # in, with no output and nothing pointing at the real cause. The first
    # submission of these five runs was lost exactly this way, to
    # -w wells_so.dat, which is res3000.in's three-well file.
    print(f"\nall {len(made)} decks verified. Submit with, from {IN}:")
    for n, lever, val in made:
        d = dict(read_deck(IN / f"res{n}.in"))
        inj = d["injection_file"].strip('"')
        perm = d["parameter_file"].strip('"')
        assert (IN / inj).exists() and (IN / inj).stat().st_size > 0, \
            f"{inj} missing or empty in {IN}"
        print(f"  sbatch march26_submit_hbi_git_scratch.sh -i res{n}.in "
              f"-w {inj} -p {perm}")


if __name__ == "__main__":
    main()
