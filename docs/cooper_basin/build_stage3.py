#!/usr/bin/env python3
"""Stage 3: sweep tau_0 on a COMMON grid at both sigmabar_0.

The question is "how understressed can the fault be and still match both
targets?" That is a question about tau_0, not about muinit, so tau_0 is the swept
axis and muinit is derived: muinit = tau_0 / sigmainit.

Consequence, agreed deliberately: at sigmainit 27.99, reaching Taiyi's
tau_0 = 15.0 needs muinit = 0.536, above the 0.50 the physical-ranges table
calls realistic. That is accepted because 0.536 is exactly Taiyi's OWN implied
friction -- his setup_model.m sets taux_as_0 = 15e6 and his stress state gives
sigmabar_0 = 27.9892, so 15.0/27.9892 = 0.5359. Using it is matching him, not
departing from the ranges arbitrarily. The alternative (a common muinit grid)
was rejected because the sigmainit 27.99 column would top out at tau_0 13.99 and
never reach his value.

GRID: tau_0 = 11.0 to 15.0 MPa in 0.5 MPa steps, 9 values, at both sigmainit,
at both parents. 36 runs, run as two blocks of 18 via --start.

fluid and permev are FIXED at B (eta 1.27e-4) and T (enhancement on), because
those are the model under test: fluid A is the Taiyi-faithful baseline and
permev F is the enhancement control, and both are already covered at
muinit 0.37 by Stage 1 and Stage 2. Carrying them into the tau_0 sweep would be
sweeping the baselines rather than the model. parent stays an axis because on the
runs available it points two ways -- 1808 gives the better front (0.73 vs 0.66)
and the far worse wellhead (+310% vs +72%).

  sigmainit 30.00 -> muinit 0.3667 ... 0.5000
  sigmainit 27.99 -> muinit 0.3930 ... 0.5359

tau_0 = 11.0 sits BELOW the sigmainit 30.0 failure floor of 11.45 MPa and ABOVE
the 27.99 floor of 10.24, so it doubles as an empirical check that the floor is
real: the 30.0 twin should not slip and the 27.99 twin should.

The base cell must be chosen from Stage 2 results -- whichever (parent, fluid,
permev) combination came closest. That is why this script REQUIRES --base and
does not guess.

Usage:
    python build_stage3.py --base 632812            # dry run, prints the grid
    python build_stage3.py --base 632812 --write
"""
import argparse
import sys
from pathlib import Path

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")

TAU0_GRID = [11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0]
SIGMAS = ["30.0", "27.99"]
DP_OBS = 10.92          # peak measured downhole overpressure, 0-5 d, MPa
START = 632830


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
    ap.add_argument("--base", type=int, required=True,
                    help="the Stage 2 run whose (parent, fluid, permev) cell won")
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--start", type=int, default=START,
                    help="first filenumber; two parents need two blocks")
    a = ap.parse_args()

    src = IN / f"res{a.base}.in"
    if not src.exists():
        sys.exit(f"no such deck: {src}")
    pairs = read_deck(src)
    d = dict(pairs)
    f0 = float(d["f0"])

    print(f"base cell: res{a.base}.in")
    print(f"  eta {d['eta']}  beta {d['beta']}  phi {d['phi']}  permev {d['permev']}")
    print(f"  perm {d.get('parameter_file', 'uniform kp ' + d.get('kp',''))}")
    print(f"  kpmax {d.get('kpmax')}  kpmin {d.get('kpmin')}  f0 {f0}\n")
    print(f"{'run':>7s} {'sigma':>7s} {'tau0':>6s} {'muinit':>8s} {'dp_crit':>8s} "
          f"{'vs 10.92':>9s} {'can slip?':>10s}")
    print("-" * 64)

    n = a.start
    made = []
    for sig in SIGMAS:
        s = float(sig)
        for tau0 in TAU0_GRID:
            mu = tau0 / s
            dpc = s * (1 - mu / f0)
            ok = dpc <= DP_OBS
            print(f"{n:>7d} {sig:>7s} {tau0:>6.1f} {mu:>8.4f} {dpc:>8.2f} "
                  f"{dpc - DP_OBS:>+9.2f} {('yes' if ok else 'NO'):>10s}")
            if a.write:
                changes = {"filenumber": str(n), "sigmainit": sig,
                           "muinit": f"{mu:.4f}"}
                hdr = [
                    "! STAGE 3 -- tau_0 sweep on a COMMON grid at both sigmabar_0.",
                    f"! Cell: tau_0 = {tau0:.1f} MPa at sigmainit {sig}"
                    f"  ->  muinit = tau_0/sigmainit = {mu:.4f}",
                    f"! Built from res{a.base}.in; only sigmainit and muinit differ.",
                    "!",
                    "! tau_0 is the swept axis, not muinit, because the question is how",
                    "! understressed the fault can be -- a statement about tau_0. A common",
                    "! tau_0 grid also makes the two sigmainit columns comparable at equal",
                    "! shear stress.",
                    "!",
                    f"! dp_crit = sigmainit(1 - muinit/f0) = {dpc:.2f} MPa against a peak",
                    f"! MEASURED downhole overpressure of {DP_OBS:.2f} MPa over 0-5 d, so this",
                    "! cell " + ("CAN fail at the observed pressure."
                                 if ok else
                                 "CANNOT fail without over-pressurising."),
                ]
                if abs(mu - 0.5359) < 0.01:
                    hdr += [
                        "!",
                        "! NOTE muinit 0.536 is above the 0.50 the physical-ranges table",
                        "! calls realistic. It is used deliberately: it is Taiyi's OWN",
                        "! implied friction (his taux_as_0 = 15e6 with sigmabar_0 27.9892",
                        "! gives 15.0/27.9892 = 0.5359), so this cell matches him rather",
                        "! than departing from the ranges arbitrarily.",
                    ]
                lines = [f"filenumber {n}"] + hdr
                for k, v in pairs:
                    if k == "filenumber":
                        continue
                    lines.append(f"{k} {changes.get(k, v)}")
                (IN / f"res{n}.in").write_text("\n".join(lines) + "\n")
                made.append(n)
            n += 1

    if not a.write:
        print("\n(dry run -- pass --write to create the decks)")
        return
    print(f"\nwrote {len(made)} decks: {made[0]}-{made[-1]}")
    print("\nverifying only sigmainit and muinit changed:")
    base = dict(read_deck(src))
    bad = 0
    for m in made:
        new = dict(read_deck(IN / f"res{m}.in"))
        diff = {k for k in set(new) | set(base) if new.get(k) != base.get(k)}
        extra = diff - {"filenumber", "sigmainit", "muinit"}
        if extra:
            bad += 1
            print(f"  res{m}.in  UNEXPECTED: {sorted(extra)}")
    print(f"  {len(made)-bad}/{len(made)} clean")
    sys.exit(0 if bad == 0 else 1)


if __name__ == "__main__":
    main()
