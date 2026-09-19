#!/usr/bin/env python3
"""CYCLE 2, Round 5: sweep kpmax at tau_0 = 10.36, kpmin fixed. Contrast 2.75x
to 400x.

WHY. kpmax = 2.5e-13 has been the value in EVERY run of this project, ~90 of
them, and it is one of the parameters with no justification recorded anywhere
in the repository. CONCLUSIONS.md records that the front is insensitive to it
across a 100x range while the wellhead plateau is not -- so it is constrained
by pressure alone, and the pressure calibration has been corrected twice since
that was written (the Holl datum, and the segmented friction). It therefore
needs re-testing rather than inheriting.

THE RANGE IS SET BY THE COMPARISON WITH WANG & DUNHAM. With kpmin held at
1e-15, the contrast kpmax/kpmin spans:

    2.75x   Taiyi's own contrast -- his fault zone is k = 1.1e-12 within 150 m
            of the well and 4.0e-13 beyond (Table 1), i.e. nearly uniform
    400x    the ratio of his far-field 4.0e-13 to our kpmin

which brackets the present 250x. Note where this puts the two models relative
to each other: OUR kpmax OF 2.5e-13 IS BELOW HIS FAR-FIELD 4.0e-13. His
reservoir is a nearly uniform, far more permeable fault zone; ours is
low-permeability rock with a small sharp conductive patch. The sweep asks
whether the data can tell them apart.

SEVEN VALUES, logarithmically spaced, because the range is a factor of 145 and
linear spacing would put five of seven points in the top decade.

kpmin AND kp STAY AT 1e-15 by instruction. So this varies the contrast and the
near-well diffusivity together, and leaves the far field alone -- D_far stays
0.0350 m^2/s in every run while D_disc spans 0.096 to 14.0 m^2/s.

TWO KEYS CHANGE PER DECK, kpmax and parameter_file, because the permeability
MAP encodes kpmax inside the disc: a deck whose kpmax disagrees with its map's
maximum is a silent inconsistency, and presubmit.py checks for exactly that. So
the pair is one physical change expressed in two places, and the verification
below requires both to move together and to agree.

STRESS STATES. Run first at tau_0 = 10.36 MPa, the MEASURED shear stress, then
extended to the next two up, 11.53 and 12.71, on the reasoning that kpmax
behaves OPPOSITELY to the disc radius as stress rises. The disc sets WHERE kp
starts high; kpmax sets the CEILING that permev grows toward, and kp relaxes
toward it at a rate proportional to slip velocity (main_LH.f90:2426) with no
healing at kT = 1e15. So at higher tau_0, where more of the fault slips, more
of it reaches kpmax -- kpmax matters MORE exactly where the disc matters less
(29 points of dp range at tau_0 10.36, only 3 at 15.00).

The two extra states are also where the answer would be most useful: 11.53 is
the current best joint cell and the sweep tests how ROBUST that match is to a
parameter with no recorded justification; 12.71 is where dp crosses zero on the
disc sweep, and if lambda is kpmax-sensitive then lowering kpmax there should
move lambda toward 1 AND dp toward 0 together.

A CAVEAT ON THE PREMISE. CONCLUSIONS.md records that the front holds at
0.98-1.03x across a 100x change in kpmax, i.e. that lambda is INSENSITIVE to
it. If that survives the corrected datum and the lambda definition, then kpmax
is a clean dp-only lever and the reasoning above about lambda is wrong. The
tau_0 10.36 sweep measures lambda(kpmax) directly and settles it. The
contrast values do not depend on the answer, which is why all three states are
built without waiting.

Parents are res632960/1/2, disc 300 m, arm 2.

WHAT TO WATCH FOR. At low contrast the disc is barely more permeable than the
background, so injectivity collapses and dp should rise steeply; that may drive
MORE slip, not less, since dp_crit is unchanged. The 2.75x and 6x runs are the
ones that might fail to reach tmax on nstep rather than on time.

Usage:
  python build_kmax_sweep.py            # dry run
  python build_kmax_sweep.py --write
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
DAYS, SAFE = 13.1, 0.8
# (parent, tau_0, MEASURED 13.1 d front in m, first filenumber). The front is
# measured, not estimated, and bounds its own block -- see the domain note.
PARENTS = {10.36: (632960, 831.0, 633080),
           11.53: (632961, 932.0, 633090),
           12.71: (632962, 1082.0, 633100)}
CONTRASTS = [2.75, 6.0, 13.0, 30.0, 65.0, 145.0, 400.0]
TAIYI_NEAR, TAIYI_FAR = 1.1e-12, 4.0e-13


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
    ap.add_argument("--tau", nargs="+", type=float, default=[10.36],
                    help="which stress states to build; parents are fixed")
    a = ap.parse_args(argv)
    for t in a.tau:
        assert t in PARENTS, f"no parent registered for tau_0 = {t}"
    allrows, allmade = [], []
    for TAU in a.tau:
        PARENT, R_MAX_PARENT, FIRST = PARENTS[TAU]
        _build(TAU, PARENT, R_MAX_PARENT, FIRST, a.write, allrows, allmade)
    if allmade:
        print(f"\nall {len(allmade)} decks verified across "
              f"{len(a.tau)} stress state(s). Submit from {IN}:")
        for n, c, kx, mn in allmade:
            d = dict(read_deck(IN / f"res{n}.in"))
            print(f"  sbatch march26_submit_hbi_git_scratch.sh -i res{n}.in "
                  f"-w {d['injection_file'].strip(chr(34))} -p {mn}")
    elif not a.write:
        print("\ndry run -- nothing written. re-run with --write")


def _build(TAU, PARENT, R_MAX_PARENT, FIRST, write, allrows, allmade):
    pairs = read_deck(IN / f"res{PARENT}.in")
    base = dict(pairs)
    ds = ff(base["ds"]) * 1000.0
    IMAX = int(base["imax"]); half = IMAX * ds / 2
    kx0, km = ff(base["kpmax"]), ff(base["kpmin"])
    eta, phi, beta = ff(base["eta"]), ff(base["phi"]), ff(base["beta"])
    sig, f0, mu = ff(base["sigmainit"]), ff(base["f0"]), ff(base["muinit"])
    den = eta * phi * beta
    import re
    m0 = re.search(r"disc(\d+)", base["parameter_file"])
    disc = float(m0.group(1))

    assert abs(f0 - 0.60) < 1e-9, f"parent f0 is {f0}"
    assert abs(mu * sig - TAU) < 0.02, \
        f"res{PARENT}: tau_0 is {mu*sig:.2f}, expected {TAU}"
    assert abs(ff(base["kp"]) - km) < 1e-20, "parent kp != kpmin"

    tag = " (MEASURED)" if abs(TAU - 10.36) < 0.02 else ""
    print(f"\nparent res{PARENT}.in: tau_0 {mu*sig:.2f} MPa{tag}, "
          f"disc {disc:.0f} m, ds {ds:.0f} m")
    print(f"  kpmax {kx0:.2e}, kpmin {km:.1e}  ->  contrast {kx0/km:.0f}x")
    print(f"  D_far = kpmin/(eta phi beta) = {km/den:.4f} m2/s, FIXED "
          f"(kpmin held)")
    print(f"  Taiyi: k_near {TAIYI_NEAR:.1e} (contrast "
          f"{TAIYI_NEAR/TAIYI_FAR:.2f}x), k_far {TAIYI_FAR:.1e}\n")
    print(f"{'run':>7} {'contrast':>9} {'kpmax':>10} {'D_disc':>8} "
          f"{'vs now':>7} {'R bound':>8} {'L/half':>7}  note")

    rows, bad = [], 0
    for i, c in enumerate(CONTRASTS):
        kx = c * km
        D = kx / den
        # CONCLUSIONS.md: the front holds 0.98-1.03x across a 100x change in
        # kpmax, so it is weakly sensitive. sqrt(kpmax ratio) is used as a
        # deliberately PESSIMISTIC bound for the domain check -- it would be
        # the scaling if the front were set entirely by the disc diffusivity,
        # which it is not, since most of its travel is through kpmin rock.
        R = R_MAX_PARENT * max(1.0, math.sqrt(kx / kx0))
        ratio = (R + 398.0) / half
        bad += ratio >= SAFE
        note = ("Taiyi's contrast" if abs(c - TAIYI_NEAR/TAIYI_FAR) < 0.1
                else "= Taiyi's k_far" if abs(kx - TAIYI_FAR) < 1e-15 else "")
        rows.append((FIRST + i, c, kx, D, R, ratio))
        print(f"{FIRST+i:>7} {c:>8.2f}x {kx:>10.2e} {D:>8.3f} "
              f"{D/(kx0/den):>6.2f}x {R:>7.0f}m {ratio:>7.2f}  {note}"
              + ("" if ratio < SAFE else "  <-- BREACH"))
    print(f"\n  present value {kx0:.1e} = {kx0/km:.0f}x sits between "
          f"{CONTRASTS[-2]:.0f}x and {CONTRASTS[-1]:.0f}x, and already exists "
          f"as res{PARENT}.in")
    allrows.extend(rows)
    if bad:
        sys.exit(f"\n{bad} deck(s) breach L/half {SAFE}. Nothing written.")
    if not write:
        return

    c_ = (IMAX - 1) // 2
    rr = np.hypot(*(np.mgrid[0:IMAX, 0:IMAX] - c_)) * ds
    made = []
    for n, c, kx, D, R, ratio in rows:
        mn = f"perm_2zone_{IMAX}_ds{int(ds)}_disc{int(disc)}_kmax{kx:.2e}.txt"
        if not (IN / mn).exists():
            kp = np.where(rr.ravel() <= disc + 1e-9, kx, km)
            (IN / mn).write_text("kp\n" + "\n".join(f"{v:.6e}" for v in kp) + "\n")
            print(f"wrote {mn}")
        hdr = [
            f"! CYCLE 2 ROUND 5: kpmax SWEEP at the MEASURED tau_0.",
            f"! This run: kpmax {kx:.3e} m^2, contrast kpmax/kpmin = {c:.2f}x.",
            "!",
            f"! kpmax = {kx0:.1e} has been the value in EVERY run of this project,",
            "! ~90 of them, and it is one of the parameters with no justification",
            "! recorded anywhere in the repository. CONCLUSIONS.md records that the",
            "! front is insensitive to it across a 100x range while the wellhead",
            "! plateau is not -- so it is constrained by PRESSURE ALONE, and the",
            "! pressure calibration has since been corrected twice (the Holl datum,",
            "! and the segmented pipe friction). It needs re-testing rather than",
            "! inheriting.",
            "!",
            "! THE RANGE IS SET BY THE COMPARISON WITH WANG & DUNHAM. With kpmin",
            f"! held at {km:.1e}, the contrast spans {CONTRASTS[0]:.2f}x -- Taiyi's own,",
            f"! his zone being {TAIYI_NEAR:.1e} within 150 m and {TAIYI_FAR:.1e} beyond,",
            f"! i.e. nearly uniform -- to {CONTRASTS[-1]:.0f}x, the ratio of his",
            "! far-field value to our kpmin. Note where that puts the two models:",
            f"! OUR PRESENT kpmax {kx0:.1e} IS BELOW HIS FAR-FIELD {TAIYI_FAR:.1e}.",
            "! His reservoir is a nearly uniform, far more permeable fault zone;",
            "! ours is low-permeability rock with a small sharp conductive patch.",
            "! This sweep asks whether the data can tell them apart.",
            "!",
            "! SEVEN VALUES, LOGARITHMICALLY SPACED, because the range is a factor",
            "! of 145 and linear spacing would put five of seven in the top decade.",
            "!",
            f"! kpmin AND kp STAY AT {km:.1e} by instruction, so the far field is",
            f"! untouched: D_far = {km/den:.4f} m^2/s in every run, while",
            f"! D_disc = kpmax/(eta phi beta) = {D:.3f} m^2/s here against",
            f"! {kx0/den:.3f} at the present value, a factor {D/(kx0/den):.2f}.",
            "!",
            "! TWO KEYS DIFFER from the parent, kpmax and parameter_file, because",
            "! the MAP encodes kpmax inside the disc. A deck whose kpmax disagrees",
            "! with its map maximum is a silent inconsistency and presubmit.py",
            "! checks for it, so the pair is one physical change written twice.",
            f"! The disc radius is unchanged at {disc:.0f} m.",
            "!",
            f"! tau_0 = {mu*sig:.2f} MPa{tag}. dp_crit is unchanged at",
            f"! sigmabar_0(1 - muinit/f0) = {sig*(1-mu/f0):.2f} MPa, so any change in",
            "! slip comes from the pressure field and not from the failure threshold.",
            "!",
            "! WHAT TO WATCH FOR: at low contrast the disc is barely more permeable",
            "! than the background, injectivity collapses and dp should rise",
            "! steeply. That may drive MORE slip rather than less, dp_crit being",
            "! fixed. The low-contrast runs are the ones that may end on nstep",
            "! rather than on tmax.",
            "!",
            f"! DOMAIN: bounded by the parent's MEASURED {R_MAX_PARENT:.0f} m at",
            f"! {DAYS} d scaled by sqrt(kpmax ratio) = {max(1.0, math.sqrt(kx/kx0)):.2f},",
            f"! giving {R:.0f} m. That is deliberately pessimistic -- CONCLUSIONS.md",
            "! measured the front as nearly insensitive to kpmax, because most of",
            f"! its travel is through kpmin rock. L/half {ratio:.2f}.",
        ]
        ov = {"kpmax": f"{kx:.3e}", "parameter_file": f'"{mn}"'}
        lines = [f"filenumber {n}"] + hdr
        for k, v in pairs:
            if k == "filenumber":
                continue
            lines.append(f"{k} {ov.get(k, v)}")
        (IN / f"res{n}.in").write_text("\n".join(lines) + "\n")
        made.append((n, c, kx, mn))

    print("  verifying")
    nbad = 0
    ref = dict(pairs)
    for n, c, kx, mn in made:
        d = dict(read_deck(IN / f"res{n}.in"))
        diff = sorted(k for k in set(ref) | set(d)
                      if ref.get(k) != d.get(k) and k != "filenumber")
        ok = diff == ["kpmax", "parameter_file"]
        mp = np.loadtxt(IN / mn, skiprows=1)
        agree = (abs(mp.max() - ff(d["kpmax"])) < 1e-20
                 and abs(mp.min() - ff(d["kpmin"])) < 1e-20)
        rad = math.sqrt(int((mp == mp.max()).sum()) / math.pi) * ds
        radok = abs(rad - disc) <= 1.5 * ds
        nbad += not (ok and agree and radok)
        print(f"  {n}  {c:>6.2f}x  changed {diff}  map max {mp.max():.3e} "
              f"min {mp.min():.1e} r {rad:.0f} m  "
              + ("OK" if ok and agree and radok else "<-- WRONG"))
    if nbad:
        sys.exit(f"\n{nbad} problem(s). Nothing submitted.")
    allmade.extend(made)


if __name__ == "__main__":
    main()
