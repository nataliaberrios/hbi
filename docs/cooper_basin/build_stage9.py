#!/usr/bin/env python3
"""Stage 9: raise kpmax on the one configuration that matches the front.

THE ARGUMENT. 632875 is the only run of 70 whose slip front matches
(lambda 1.03), but its wellhead is +80%. 28 runs match the wellhead, and 24 of
those produce essentially no slip. Nothing matches both. Every attempt so far
moved a lever that shifts BOTH targets, so they traded along a line.

kpmax does not. Measured, from 632812 vs 632818, which differ mainly in kpmax:

                    kpmax     peak above plateau     R      wellhead
    632812        2.5e-13           6.37 MPa       235 m     +70.7%
    632818        5.0e-14          31.90 MPa       250 m     +309.1%

The front moves 6% while the wellhead moves 4.4x. The peak above the plateau
scales as 1/kpmax, exactly as deltaP = q*eta/(4*pi*kpmax) predicts -- kpmax fell
5x and the peak rose 5.0x. So:

  * the FRONT is set by volume balance into the sealed disc, which depends on
    phi*beta and dp_crit but NOT on kpmax
  * the WELLHEAD is the plateau (pinned near dp_crit = sigma_0 - tau_0/f by the
    slip-enhancement feedback) PLUS a logarithmic peak that scales as 1/kpmax

Every run to date sat at kpmax ~2.5e-13, so the peak was always 6-7 MPa on top
of the plateau and always blew the wellhead out. Raising kpmax should collapse
the peak without touching the front.

PREDICTION, from 632875's own numbers (plateau 12.15 MPa, peak above it
6.81 MPa, Delta-p(r_w) 18.96 MPa against a measured 10.92):

     kpmax      peak above plateau     dp(r_w)     wellhead bias
    2.5e-13           6.81 MPa         18.96        +80%   (632875, measured)
    2.5e-12           0.68             12.83        ~+4%
    2.5e-11           0.068            12.22        ~+3%
    2.5e-10           0.0068           12.16        ~+3%

with the front staying at lambda ~1.03 throughout. The floor is the plateau
itself, 12.15 against 10.92, so the wellhead cannot come below about +3% by this
route alone -- but +3% is inside the +/-15% band, and the front is already
inside its band. That would be both, at muinit 0.370: understressed.

Three points rather than one, because the value of this is the TREND. If the
peak falls as 1/kpmax and the front holds, the mechanism is confirmed whatever
the absolute numbers do. If the front moves, the volume-balance picture is wrong
and one run at one kpmax would not have told us which.

CONSISTENCY. The rule is that the near-well disc IS kpmax and the background IS
kp = kpmin. So each deck gets its OWN map with the disc raised to match its
kpmax; kp and kpmin stay at 1e-15. The phi grading (G3: 0.005 near, 0.020 far)
and every other value are copied from 632875 unchanged.

Usage:  python build_stage9.py [--write]
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
PARENT = 632875
KPMAX = [("632896", 2.5e-12), ("632897", 2.5e-11), ("632898", 2.5e-10)]


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
    src_map = base["parameter_file"].strip('"')
    names = (IN / src_map).read_text().split("\n", 1)[0].split()
    arr = np.loadtxt(IN / src_map, skiprows=1)
    ikp, iph = names.index("kp"), names.index("phi")
    kp_old, phi = arr[:, ikp].copy(), arr[:, iph]
    old_max, old_min = kp_old.max(), kp_old.min()
    disc = kp_old == old_max
    eta = float(base["eta"].replace("d", "e"))
    beta = float(base["beta"])

    print(f"parent res{PARENT}.in   map {src_map}")
    print(f"  disc {disc.sum()} cells at {old_max:.1e}, background "
          f"{(~disc).sum()} at {old_min:.1e}")
    print(f"  phi {phi[disc][0]:.3f} near / {phi[~disc][0]:.3f} far  (G3)")
    print()
    print(f"{'run':>8s} {'kpmax':>10s} {'range':>9s} {'D disc':>10s} "
          f"{'D far':>9s} {'gamma':>8s} {'peak pred':>10s}")
    ds_km = float(base["ds"])
    den = math.log(0.2 * ds_km * 1e3 / float(base["rw"].replace("d", "e")))
    sw = float(base["Sw_fwid"].replace("d", "e"))
    made = []
    for tag, kx in KPMAX:
        D_disc = kx / (eta * phi[disc][0] * beta)
        D_far = old_min / (eta * phi[~disc][0] * beta)
        T = 2 * math.pi * kx / eta / den
        g = (sw / 60.0) / ((sw / 60.0) + T)
        # peak above plateau scales as 1/kpmax from 632875's measured 6.81 MPa
        print(f"{tag:>8s} {kx:>10.1e} {kx/old_min:>8.0f}x {D_disc:>10.1f} "
              f"{D_far:>9.2e} {g:>8.5f} {6.81*old_max/kx:>9.3f} MPa")
        if not a.write:
            continue
        newmap = f"permphi_2zone_601_ds5_kmax{kx:.1e}_G3.txt"
        kp_new = np.where(disc, kx, old_min)
        out = [" ".join(names)]
        out += [f"{k:.6e} {p:.6e}" for k, p in zip(kp_new, phi)]
        (IN / newmap).write_text("\n".join(out) + "\n")
        hdr = [
            "! STAGE 9 -- raise kpmax on the one configuration that matches the front.",
            f"! Built from res{PARENT}.in; only kpmax and the map's disc value change.",
            "!",
            "! 632875 is the only run of 70 whose front matches (lambda 1.03), but its",
            "! wellhead is +80%. 28 runs match the wellhead and 24 of those barely slip.",
            "! Nothing matches both, because every lever tried so far moves both targets.",
            "!",
            "! kpmax does not. From 632812 vs 632818, which differ mainly in kpmax:",
            "!   632812  kpmax 2.5e-13   peak above plateau  6.37 MPa   R 235 m   +70.7%",
            "!   632818  kpmax 5.0e-14   peak above plateau 31.90 MPa   R 250 m   +309%",
            "! The front moves 6% while the wellhead moves 4.4x, and the peak scales as",
            "! 1/kpmax exactly as deltaP = q*eta/(4*pi*kpmax) predicts. The front is set",
            "! by volume balance into the sealed disc (phi*beta and dp_crit, NOT kpmax);",
            "! the wellhead is a plateau pinned near dp_crit PLUS a peak going as 1/kpmax.",
            "!",
            f"! kpmax {old_max:.1e} -> {kx:.1e} ({kx/old_max:.0f}x). Predicted peak above",
            f"! plateau {6.81*old_max/kx:.3f} MPa against 632875's 6.81, so dp(r_w)",
            f"! {12.15 + 6.81*old_max/kx:.2f} MPa against a measured 10.92, and a wellhead",
            "! bias of roughly +3 to +4%. The front should be unchanged at ~1.03.",
            "!",
            "! CONSISTENCY: the near-well disc IS kpmax and the background IS kp = kpmin.",
            f"! {newmap} therefore has its {disc.sum()}-cell disc raised to {kx:.1e} while",
            f"! kp and kpmin stay at {old_min:.1e}. The G3 phi grading and every other",
            "! value are copied from the parent unchanged.",
            "!",
            f"! Domain: the far field stays at kpmin, D = {D_far:.2e} m^2/s, so the 5 d",
            f"! diffusion length is {math.sqrt(4*D_far*5*86400):.0f} m against a 1502 m",
            "! half-domain. Inside the disc D is large but the disc is bounded by the",
            "! unslipped kpmin region.",
            "!",
            f"! gamma = {g:.5f} at h = 60 s, against 632875's 0.14 -- a high kpmax makes",
            "! the Peaceman well WELL damped, so this should integrate more easily than",
            "! the parent, not less.",
        ]
        lines = [f"filenumber {tag}"] + hdr
        for k, v in pairs:
            if k == "filenumber":
                continue
            if k == "kpmax":
                lines.append(f"kpmax {kx:.1e}")
            elif k == "parameter_file":
                lines.append(f'parameter_file "{newmap}"')
            else:
                lines.append(f"{k} {v}")
        (IN / f"res{tag}.in").write_text("\n".join(lines) + "\n")
        made.append((tag, kx, newmap))

    if not a.write:
        print("\ndry run -- nothing written. re-run with --write")
        return

    print("\nverifying: map disc == kpmax, background == kpmin == kp, phi untouched,")
    print("and the deck differs from the parent only in filenumber/kpmax/parameter_file")
    bad = 0
    for tag, kx, newmap in made:
        d = dict(read_deck(IN / f"res{tag}.in"))
        v = np.loadtxt(IN / newmap, skiprows=1)
        diff = {k for k in set(d) | set(base) if d.get(k) != base.get(k)}
        ck = {
            "only 3 keys changed": not (diff - {"filenumber", "kpmax",
                                                "parameter_file"}),
            "kpmax == map max": abs(v[:, ikp].max() - kx) / kx < 1e-6,
            "kpmin == map min": abs(v[:, ikp].min() - old_min) / old_min < 1e-6,
            "kp == kpmin": d["kp"] == d["kpmin"],
            "phi identical": np.allclose(v[:, iph], phi),
            "disc size same": int((v[:, ikp] == v[:, ikp].max()).sum()) == int(disc.sum()),
        }
        ok = all(ck.values())
        bad += not ok
        print(f"  res{tag}.in kpmax {kx:.1e}: "
              + ("OK" if ok else "PROBLEM " + str([k for k, x in ck.items() if not x])))
    print(f"\n  {len(made)-bad}/{len(made)} clean")
    sys.exit(0 if bad == 0 else 1)


if __name__ == "__main__":
    main()
