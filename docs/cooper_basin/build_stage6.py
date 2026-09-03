#!/usr/bin/env python3
"""Stage 6: permeability enhancement ON, on the Taiyi configuration.

WHY THIS PAIR AND NOT ANOTHER POINT IN THE GRID.

Of 68 runs, 27 match the wellhead to within 15%. TWENTY-FIVE of those produce no
slip whatsoever -- peak slip / dc is exactly 0. Only two match the wellhead AND
slip, and they are the two Taiyi-parameter runs:

    632881   wellhead  +8.2%   front 0.22   pk/dc  68   permev F
    632880   wellhead +11.1%   front 0.14   pk/dc 253   permev F

Exactly one run matches the front:

    632875   front 1.03        wellhead +80.0%          permev T

So there are two possible bases, and they are not equally good. 632875 needs its
wellhead overshoot cut 5.3x, and the wellhead is a direct measurement. 632881
already matches the wellhead and already slips; its front needs to grow 4.5x.

And 632880/632881 have permev **F**. The hypothesis this whole project rests on
-- that permeability enhancement, which Wang & Dunham did not have, plus a
nonuniform initial permeability, reaches the match at lower tau_0 -- has never
been tested on the one configuration that matches the wellhead. That is the gap
these two runs close.

CONSISTENCY WITH THE GOVERNING RULE. kpmax is the ceiling, kp = kpmin is the
floor and the initial value; when the initial permeability is nonuniform the
near-well disc IS kpmax and the background IS kpmin. perm_taiyi_601_ds20.txt
already satisfies this: verified 177 cells at 1.1e-12 (equivalent radius 150 m)
and 361024 at 4e-13, which are the paper's Table 1 values, and the parent decks
already carry kp 4e-13. So enabling enhancement adds kpmax 1.1e-12 and
kpmin 4e-13 and invents nothing. Asserted below rather than assumed.

DOMAIN IS SAFE. Enhancement raises k only where slip occurs, and slip reaches at
most ~80 m in the parents. The far field stays at kpmin = 4e-13, so
D = k/(eta*phi*beta) stays 4.49 m^2/s and the 5 d diffusion length stays 2.79 km
against a 6.01 km half-domain -- 46%, the same as the parents. This is why the
pair is run at ds 20 m and not the 5 m used by the other 66: at ds 5 m the
601-cell grid is only 1.5 km half-width and 2.79 km would run off the edge.

WHAT WOULD FALSIFY THE HYPOTHESIS. If the front stays near 0.22 with enhancement
on, then enhancement is not the missing ingredient at a wellhead-matching
pressure, and the 4.5x has to come from somewhere else. If the wellhead climbs
out of band as the front grows, then enhancement is on the same trade-off line as
every other lever and does not span the plane either.

Usage:  python build_stage6.py            # dry run, prints what it would write
        python build_stage6.py --write
"""
import argparse
import sys
from pathlib import Path

import numpy as np

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")

PAIRS = [(632884, 632880), (632885, 632881)]
MAP = "perm_taiyi_601_ds20.txt"


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

    # Bounds come from the map, never from a literal typed here.
    k = np.loadtxt(IN / MAP, skiprows=1)
    kpmax, kpmin = k.max(), k.min()
    n_near = int((k == kpmax).sum())
    print(f"{MAP}: {k.size} cells, kpmax {kpmax:.4e} x{n_near}, "
          f"kpmin {kpmin:.4e} x{int((k == kpmin).sum())}")

    made = []
    for new, parent in PAIRS:
        pairs = read_deck(IN / f"res{parent}.in")
        d = dict(pairs)

        # Refuse to build on a parent that is not what this script assumes.
        assert d["permev"].upper() == "F", \
            f"res{parent}.in already has permev {d['permev']}"
        assert d["parameter_file"].strip('"') == MAP, \
            f"res{parent}.in uses {d['parameter_file']}, not {MAP}"
        assert abs(float(d["kp"].replace("d", "e")) - kpmin) / kpmin < 1e-9, \
            f"res{parent}.in kp {d['kp']} != map min {kpmin:.4e}"
        assert "kpmax" not in d and "kpmin" not in d, \
            f"res{parent}.in already declares kpmax/kpmin"

        dc = d["dc"]
        print(f"\nres{new}.in  <- res{parent}.in   dc {dc}   "
              f"permev F -> T,  +kpmax {kpmax:.3e}  +kpmin {kpmin:.3e}")
        if not a.write:
            continue

        hdr = [
            "! STAGE 6 -- permeability enhancement ON, on the Taiyi configuration.",
            f"! Built from res{parent}.in; permev F -> T plus the two bounds. Nothing else.",
            "!",
            "! WHY. Of 68 runs, 27 match the wellhead within 15% and 25 of those",
            "! produce NO slip at all (pk/dc exactly 0). Only two match the wellhead",
            "! and also slip -- 632881 (+8.2%, front 0.22) and 632880 (+11.1%, front",
            "! 0.14) -- and both have permev F. Exactly one run matches the front,",
            "! 632875 at 1.03, but its wellhead is +80%. So the wellhead-matching",
            "! base needs its front to grow 4.5x, while the front-matching base needs",
            "! its overshoot cut 5.3x against a direct measurement. This pair tests",
            "! enhancement -- the project's central hypothesis -- on the base that",
            "! already matches the wellhead, which has never been done.",
            "!",
            "! BOUNDS. kpmax is the ceiling, kp = kpmin the floor and initial value;",
            "! the near-well disc IS kpmax and the background IS kpmin. Verified from",
            f"! {MAP}: {n_near} cells at {kpmax:.4e} (equivalent radius 150 m),",
            f"! the rest at {kpmin:.4e}, which are the paper's Table 1 values. The",
            f"! parent already carries kp {kpmin:.0e}, so nothing is invented here.",
            "!",
            "! DOMAIN. Enhancement raises k only where slip occurs (<= ~80 m in the",
            "! parents), so the far field stays at kpmin and D stays 4.49 m^2/s. The",
            "! 5 d diffusion length is 2.79 km against a 6.01 km half-domain (46%),",
            "! unchanged from the parent. ds is 20 m rather than the 5 m used by the",
            "! other 66 runs precisely for this reason: at ds 5 m the half-domain is",
            "! only 1.5 km and 2.79 km would run off the edge.",
            "!",
            f"! dc {dc} -- inherited from the parent, so this run is directly",
            f"! comparable to res{parent}.in and differs from it in permeability",
            "! evolution alone.",
        ]
        lines = [f"filenumber {new}"] + hdr
        for key, val in pairs:
            if key == "filenumber":
                continue
            if key == "permev":
                lines.append("permev T")
                lines.append(f"kpmax {kpmax:.4e}")
                lines.append(f"kpmin {kpmin:.4e}")
                continue
            lines.append(f"{key} {val}")
        (IN / f"res{new}.in").write_text("\n".join(lines) + "\n")
        made.append((new, parent))

    if not a.write:
        print("\ndry run -- nothing written. re-run with --write")
        return

    print("\nverifying each deck differs from its parent ONLY in the intended keys")
    bad = 0
    for new, parent in made:
        pd = dict(read_deck(IN / f"res{parent}.in"))
        nd = dict(read_deck(IN / f"res{new}.in"))
        diff = {k for k in set(pd) | set(nd) if pd.get(k) != nd.get(k)}
        extra = diff - {"filenumber", "permev", "kpmax", "kpmin"}
        kx = abs(float(nd["kpmax"]) - kpmax) / kpmax < 1e-9
        kn = abs(float(nd["kpmin"]) - kpmin) / kpmin < 1e-9
        ok = (not extra) and kx and kn and nd["permev"] == "T"
        bad += not ok
        print(f"  res{new}.in  changed {sorted(diff)}  "
              f"kpmax==map.max {kx}  kpmin==map.min {kn}  "
              f"{'OK' if ok else 'PROBLEM ' + str(sorted(extra))}")
    print(f"\n  {len(made)-bad}/{len(made)} clean")
    sys.exit(0 if bad == 0 else 1)


if __name__ == "__main__":
    main()
