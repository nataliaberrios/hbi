#!/usr/bin/env python3
"""Stage 10: raise tau_0 slightly at the kpmax where the wellhead is already in band.

Stage 9 established the lever structure, and it is not what I first described.
Raising kpmax lowers the near-well PEAK as 1/kpmax, as predicted -- but it lowers
the PLATEAU too, and slip amplitude tracks the plateau's excess over
dp_crit = sigma_0 - tau_0/f:

     run     kpmax   plateau   peak above   dp(r_w)   plateau-dp_crit   peak slip
  632875   2.5e-13     11.90        7.05     18.96            +1.18      6.03 cm
  632896   2.5e-12     10.66        0.70     11.36            -0.06      1.27
  632897   2.5e-11      9.78        0.07      9.85            -0.94      0.45
  632898   2.5e-10      8.63        0.01      8.64            -2.09      0.03

with dp_crit = 10.72 MPa at muinit 0.370. So kpmax alone cannot deliver both: it
buys the wellhead by draining the very pressure that drives slip.

TWO FACTS THAT MAKE A WINDOW EXIST ANYWAY.

1. The wellhead bias has a FLOOR near +10.7%, reached by 632898 even at
   dp(r_w) = 8.64 MPa. It is not set by the near-well pressure -- it is the
   shut-in periods, which HBI cannot follow because it has no wellbore
   bleed-off. So pushing kpmax past 2.5e-11 costs slip and buys nothing.

2. At kpmax 2.5e-11 the wellhead is ALREADY in band, +12.2%. The plateau there
   is 9.78 MPa, fixed by kpmax. The only missing ingredient is slip, and slip
   needs plateau > dp_crit.

So lower dp_crit rather than raise the plateau. dp_crit = sigma_0 - tau_0/f, and
Stage 3 measured that tau_0 moves slip while leaving the wellhead flat to +/-0.1%
across 36 runs. Sweeping muinit at fixed kpmax 2.5e-11 should therefore raise
slip without disturbing the wellhead match:

    muinit   tau_0    dp_crit   plateau - dp_crit
     0.370   10.36      10.72       -0.94   (632897, measured: 0.45 cm)
     0.385   10.78      10.03       -0.25
     0.397   11.11       9.47       +0.31
     0.410   11.48       8.86       +0.92

632875 reached 6.03 cm at an excess of +1.18, so +0.31 should land near the
observed 2.81 cm at 5 d. tau_0 = 11.11 MPa is still 26% below Wang & Dunham's
15.0, so this stays an understressed fault.

WHAT WOULD FALSIFY IT. If the wellhead moves with muinit here, the Stage-3
independence does not survive at high kpmax and the two targets are coupled
after all. If slip rises but the FRONT grows past the band, the front was not
volume-controlled and Stage 9's agreement was luck.

Only muinit changes. The map, kpmax, phi grading and everything else are copied
from res632897.in unchanged.

Usage:  python build_stage10.py [--write]
"""
import argparse
import sys
from pathlib import Path

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
PARENT = 632897
MUINIT = [("632900", 0.385), ("632901", 0.397), ("632902", 0.410)]
PLATEAU = 9.78          # MPa, measured for 632897 at 5 d


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
    sig, f0 = float(base["sigmainit"]), float(base["f0"])
    mu0 = float(base["muinit"])
    print(f"parent res{PARENT}.in   sigmainit {sig}  f0 {f0}  muinit {mu0}")
    print(f"  kpmax {base['kpmax']}  map {base['parameter_file']}")
    print(f"  measured plateau at 5 d: {PLATEAU:.2f} MPa\n")
    print(f"{'run':>8s} {'muinit':>7s} {'tau_0':>7s} {'dp_crit':>8s} "
          f"{'excess':>7s}   (632875 gave 6.03 cm at +1.18)")
    print(f"{PARENT:>8d} {mu0:>7.3f} {mu0*sig:>7.2f} {sig-mu0*sig/f0:>8.2f} "
          f"{PLATEAU-(sig-mu0*sig/f0):>+7.2f}   measured 0.45 cm")
    made = []
    for tag, m in MUINIT:
        dpc = sig - m * sig / f0
        print(f"{tag:>8s} {m:>7.3f} {m*sig:>7.2f} {dpc:>8.2f} "
              f"{PLATEAU-dpc:>+7.2f}")
        if not a.write:
            continue
        hdr = [
            "! STAGE 10 -- raise tau_0 slightly where the wellhead is already in band.",
            f"! Built from res{PARENT}.in; ONLY muinit changes.",
            "!",
            "! Stage 9 showed kpmax lowers the near-well peak as 1/kpmax but lowers the",
            "! PLATEAU with it, and slip tracks plateau - dp_crit:",
            "!   632875 kpmax 2.5e-13  plateau 11.90  excess +1.18  ->  6.03 cm",
            "!   632896 kpmax 2.5e-12  plateau 10.66  excess -0.06  ->  1.27 cm",
            "!   632897 kpmax 2.5e-11  plateau  9.78  excess -0.94  ->  0.45 cm",
            "!   632898 kpmax 2.5e-10  plateau  8.63  excess -2.09  ->  0.03 cm",
            "! So kpmax alone cannot give both: it drains the pressure that drives slip.",
            "!",
            "! But the wellhead bias FLOORS near +10.7% -- 632898 reaches it at",
            "! dp(r_w) = 8.64 MPa -- because the residual is the shut-in periods, which",
            "! HBI cannot follow without wellbore bleed-off, not the near-well pressure.",
            "! At kpmax 2.5e-11 the wellhead is already in band at +12.2%, so the only",
            "! missing ingredient is slip, and slip needs plateau > dp_crit.",
            "!",
            f"! muinit {mu0:.3f} -> {m:.3f}, so tau_0 {mu0*sig:.2f} -> {m*sig:.2f} MPa and",
            f"! dp_crit {sig-mu0*sig/f0:.2f} -> {dpc:.2f}, giving an excess of",
            f"! {PLATEAU-dpc:+.2f} MPa against the measured plateau of {PLATEAU:.2f}.",
            f"! 632875 reached 6.03 cm at an excess of +1.18, so this should give",
            f"! roughly {6.03*max(PLATEAU-dpc,0)/1.18:.1f} cm against the observed 2.81 cm at 5 d.",
            "!",
            f"! tau_0 = {m*sig:.2f} MPa is still {100*(1-m*sig/15.0):.0f}% below Wang &",
            "! Dunham's 15.0, so this remains an understressed fault.",
            "!",
            "! Stage 3 measured the wellhead flat to +/-0.1% across a 36-run tau_0",
            "! sweep. If it moves here, that independence does not survive at high",
            "! kpmax. If slip rises but the front leaves its band, the front was not",
            "! volume-controlled and Stage 9's agreement was luck.",
        ]
        lines = [f"filenumber {tag}"] + hdr
        for k, v in pairs:
            if k == "filenumber":
                continue
            lines.append(f"muinit {m:.4f}" if k == "muinit" else f"{k} {v}")
        (IN / f"res{tag}.in").write_text("\n".join(lines) + "\n")
        made.append((tag, m))

    if not a.write:
        print("\ndry run -- nothing written. re-run with --write")
        return
    print("\nverifying only muinit changed")
    bad = 0
    for tag, m in made:
        d = dict(read_deck(IN / f"res{tag}.in"))
        diff = {k for k in set(d) | set(base) if d.get(k) != base.get(k)}
        ok = (not (diff - {"filenumber", "muinit"})
              and abs(float(d["muinit"]) - m) < 1e-9)
        bad += not ok
        print(f"  res{tag}.in muinit {m:.3f}: changed {sorted(diff)}  "
              + ("OK" if ok else "PROBLEM"))
    print(f"\n  {len(made)-bad}/{len(made)} clean")
    sys.exit(0 if bad == 0 else 1)


if __name__ == "__main__":
    main()
