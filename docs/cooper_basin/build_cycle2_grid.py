#!/usr/bin/env python3
"""CYCLE 2, Round 3: fill the (tau_0 x disc radius) grid at arm 2. Disc <= 400 m.

WHAT EXISTED BEFORE THIS. A CROSS, not a grid. tau_0 was swept at disc 300 m
only (632960-632964) and the disc was swept at tau_0 = 10.36 MPa only
(632993 at 450 m, 632994 at 600 m). Seven of fifteen cells, and every interior
cell empty:

              300    450    600
    10.36     RUN    RUN    RUN
    11.53     RUN     .      .
    12.71     RUN     .      .
    13.86     RUN     .      .
    15.00     RUN     .      .

WHY THE INTERIOR MATTERS. The two levers do different things, measured:

    disc 300 -> 450 -> 600 at tau_0 10.36:
        lambda  187.8 -> 154.4 -> 129.1   (0.89 -> 0.73 -> 0.61 x observed)
        dp      +55.3 -> +32.6 ->  +2.0 %

    tau_0 10.36 -> 15.00 at disc 300:
        lambda  187.8 -> 411.4            (0.89 -> 1.96 x observed)
        dp      +55.3 ->  -17.5 %

So tau_0 raises lambda and lowers dp, while the disc lowers BOTH. Neither alone
can put lambda at 1.00x with dp at zero -- tau_0 13.86 gets dp to -0.2% but
overshoots lambda to 1.49x -- but a combination can, because the two move
lambda in opposite directions per unit of dp reduction. A first-order model
fitted to the cross above, lambda = lambda_300(tau) x f(disc) and
dp = dp_300(tau) + g(disc), predicts the dp = 0 locus as:

    disc 300   tau_0 13.85   lambda 1.49x
    disc 400   tau_0 12.9    lambda 1.10x     <- the best this grid can reach
    disc 450   tau_0 12.44   lambda 0.96x
    disc 600   tau_0 10.49   lambda 0.62x

DISC IS CAPPED AT 400 m BY INSTRUCTION, and the cost of that cap is on the
record above: the predicted best joint fit inside the cap is lambda 1.10x at
dp = 0, against 0.96x if 450 m were allowed. 632993 and 632994 already exist at
450 and 600 m and are not being deleted, but nothing new goes past 400.

THE GRID: disc 350 and 400 m x the five existing tau_0 values, 10 runs, filling
a 5 x 3 grid with the 300 m column already done. Every run differs from its
neighbour in the row by muinit alone and from its neighbour in the column by
parameter_file alone.

DOMAIN. The disc REDUCES the front rather than growing it -- bigger disc, lower
dp, less overpressure, less slip, smaller front, measured as 831 -> 821 -> 811 m
at tau_0 10.36 for disc 300 -> 450 -> 600. So the worst cell for the domain is
the one already run: tau_0 15.00 at disc 300 (632964), which reached 1713 m,
L/half = 0.70. Enlarging the disc can only reduce that. An earlier version of
this reasoning had the sign backwards and refused a 900 m disc for a boundary
breach it would not have caused; the estimate here is anchored on measured
fronts instead.

Usage:
  python build_cycle2_grid.py            # dry run
  python build_cycle2_grid.py --write
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
PARENT = 632960                       # arm 2, tau_0 10.36, disc 300
DAYS, SAFE = 13.1, 0.8
DISC_MAX = 400.0
DISCS = [350.0, 400.0]
MUINIT = [0.3700, 0.4120, 0.4540, 0.4950, 0.5360]
FIRST = 633000
# measured, from the existing cross -- used only for the domain check
R_MAX_300 = {0.3700: 831.0, 0.4120: 932.0, 0.4540: 1082.0,
             0.4950: 1312.0, 0.5360: 1713.0}


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
    sig, f0 = ff(base["sigmainit"]), ff(base["f0"])
    D_far = ff(base["kp"]) / (ff(base["eta"]) * ff(base["phi"]) * ff(base["beta"]))
    Ld = math.sqrt(4 * D_far * DAYS * 86400)

    assert abs(f0 - 0.60) < 1e-9, f"parent f0 is {f0}, expected 0.60 (FIXED)"
    assert max(DISCS) <= DISC_MAX, \
        f"disc {max(DISCS)} exceeds the {DISC_MAX:.0f} m cap"
    assert "skin" not in base, "parent sets skin; the grid assumes it absent"

    print(f"parent res{PARENT}.in: arm 2, ds {ds:.0f} m, imax {IMAX}, disc 300 m")
    print(f"  disc cap {DISC_MAX:.0f} m; sweeping {DISCS} x {len(MUINIT)} tau_0")
    print(f"  far-field L({DAYS} d) = {Ld:.0f} m, half-domain {half:.0f} m, "
          f"usable front <= {SAFE*half-Ld:.0f} m\n")
    print(f"{'run':>7} {'disc':>6} {'muinit':>7} {'tau_0':>7} {'R pred':>7} "
          f"{'L/half':>7}")

    rows, bad, n = [], 0, FIRST
    for dr in DISCS:
        for mu in MUINIT:
            tau = mu * sig
            # the disc reduces the front, so the disc-300 measurement is an
            # UPPER bound for any larger disc at the same tau_0
            R = R_MAX_300[mu]
            ratio = (R + Ld) / half
            bad += ratio >= SAFE
            rows.append((n, dr, mu, tau, R, ratio))
            print(f"{n:>7} {dr:>5.0f}m {mu:>7.4f} {tau:>6.2f}M {R:>6.0f}m "
                  f"{ratio:>7.2f}" + ("" if ratio < SAFE else "  <-- BREACH"))
            n += 1
    if bad:
        sys.exit(f"\n{bad} deck(s) breach L/half = {SAFE}. Nothing written.")
    print(f"\n  all {len(rows)} within L/half < {SAFE}; bounds are the measured "
          f"disc-300 fronts,\n  which over-estimate because a larger disc "
          f"SHRINKS the front")

    if not a.write:
        print("\ndry run -- nothing written. re-run with --write")
        return

    c = (IMAX - 1) // 2
    rr = np.hypot(*(np.mgrid[0:IMAX, 0:IMAX] - c)) * ds
    maps = {}
    for dr in DISCS:
        mn = f"perm_2zone_{IMAX}_ds{int(ds)}_disc{int(dr)}_kmax{kx:.1e}.txt"
        if not (IN / mn).exists():
            kp = np.where(rr.ravel() <= dr + 1e-9, kx, km)
            (IN / mn).write_text("kp\n" + "\n".join(f"{v:.6e}" for v in kp) + "\n")
            print(f"wrote {mn}  ({int((kp == kx).sum())} disc cells)")
        maps[dr] = mn

    made = []
    for n, dr, mu, tau, R, ratio in rows:
        hdr = [
            f"! CYCLE 2 ROUND 3: THE (tau_0 x disc) GRID AT ARM 2.",
            f"! This cell: disc {dr:.0f} m, muinit {mu:.4f}, tau_0 {tau:.2f} MPa.",
            "!",
            "! What existed before was a CROSS, not a grid: tau_0 was swept at",
            "! disc 300 m only (632960-4) and the disc at tau_0 10.36 only (632993",
            "! at 450 m, 632994 at 600 m). Seven of fifteen cells; every interior",
            "! cell empty. This fills the interior at disc <= 400 m.",
            "!",
            "! WHY THE INTERIOR. Measured, the two levers differ:",
            "!   disc 300->450->600 at tau_0 10.36: lambda 187.8->154.4->129.1,",
            "!                                     dp +55.3->+32.6->+2.0%",
            "!   tau_0 10.36->15.00 at disc 300:    lambda 187.8->411.4,",
            "!                                     dp +55.3->-17.5%",
            "! tau_0 raises lambda and lowers dp; the disc lowers BOTH. So neither",
            "! alone reaches lambda 1.00x with dp 0 -- tau_0 13.86 gets dp to -0.2%",
            "! but lambda to 1.49x -- while a combination can.",
            "!",
            "! PREDICTION ON RECORD, from a first-order fit to the cross",
            "! (lambda = lambda_300(tau) x f(disc), dp = dp_300(tau) + g(disc)),",
            "! for the dp = 0 locus:",
            "!     disc 300  tau_0 13.85  lambda 1.49x",
            "!     disc 400  tau_0 12.9   lambda 1.10x   <- best inside the cap",
            "!     disc 450  tau_0 12.44  lambda 0.96x",
            "! DISC IS CAPPED AT 400 m BY INSTRUCTION and the cost is stated: the",
            "! best joint fit reachable here is lambda ~1.10x at dp ~0, against",
            "! 0.96x if 450 m were allowed. If the grid confirms the model, that",
            "! cap is the thing to revisit -- not another parameter.",
            "!",
            f"! ONE KEY from each neighbour: muinit alone along the row,",
            f"! parameter_file alone down the column. Everything else is",
            f"! res{PARENT}.in's, including permev T, kL 1e-3, kT 1e15 (healing",
            f"! effectively disabled), kpmax {kx:.1e} inside the disc and",
            f"! kpmin {km:.1e} outside, f0 {f0:.2f} FIXED, sigmainit {sig:.2f} FIXED.",
            "!",
            f"! tau_0 = muinit*sigmainit = {mu:.4f} x {sig:.2f} = {tau:.2f} MPa;",
            f"! dp_crit = sigmabar_0(1 - muinit/f0) = {sig*(1-mu/f0):.2f} MPa.",
            "!",
            f"! DOMAIN: the disc SHRINKS the front (831->821->811 m at tau_0 10.36",
            f"! for disc 300->450->600), so the measured disc-300 front at this",
            f"! tau_0, {R:.0f} m, is an UPPER bound here. Plus {Ld:.0f} m far-field",
            f"! against {half:.0f} m half-domain gives L/half {ratio:.2f}, inside",
            f"! the {SAFE} limit. Worst cell in the whole grid is one already run,",
            "! 632964 at tau_0 15.00 disc 300, which reached 1713 m (L/half 0.70).",
        ]
        ov = {"muinit": f"{mu:.4f}", "parameter_file": f'"{maps[dr]}"'}
        lines = [f"filenumber {n}"] + hdr
        for k, v in pairs:
            if k == "filenumber":
                continue
            lines.append(f"{k} {ov.get(k, v)}")
        (IN / f"res{n}.in").write_text("\n".join(lines) + "\n")
        made.append((n, dr, mu))

    print("\nverifying each deck against the parent")
    nbad = 0
    ref = dict(pairs)
    for n, dr, mu in made:
        d = dict(read_deck(IN / f"res{n}.in"))
        diff = sorted(k for k in set(ref) | set(d)
                      if ref.get(k) != d.get(k) and k != "filenumber")
        # parameter_file always differs; muinit differs only where this cell's
        # tau_0 is not the parent's own. Demanding both unconditionally flagged
        # the two tau_0 = 10.36 cells, which are correct.
        want = ["parameter_file"] if abs(mu - ff(ref["muinit"])) < 1e-9 \
            else ["muinit", "parameter_file"]
        ok = diff == want
        nbad += not ok
        mp = np.loadtxt(IN / d["parameter_file"].strip('"'), skiprows=1)
        rad = math.sqrt(int((mp == mp.max()).sum()) / math.pi) * ds
        radok = abs(rad - dr) <= 1.5 * ds and abs(mp.max() - kx) < 1e-20
        nbad += not radok
        print(f"  {n}  disc {dr:.0f} m  changed {diff}  "
              f"map r_eff {rad:.0f} m  "
              + ("OK" if ok and radok else f"<-- expected {want}"))
    if nbad:
        sys.exit(f"\n{nbad} problem(s). Nothing submitted.")
    print(f"\nall {len(made)} verified. Submit from {IN}:")
    for n, dr, mu in made:
        d = dict(read_deck(IN / f"res{n}.in"))
        inj, perm = d["injection_file"].strip('"'), d["parameter_file"].strip('"')
        assert (IN / inj).exists() and (IN / inj).stat().st_size > 0
        print(f"  sbatch march26_submit_hbi_git_scratch.sh -i res{n}.in "
              f"-w {inj} -p {perm}")


if __name__ == "__main__":
    main()
