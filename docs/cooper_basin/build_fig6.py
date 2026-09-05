#!/usr/bin/env python3
"""Build the two runs Fig 6 needs: 807's physics on current code, and 632510 to 17 d.

Fig 6 has three panels -- observed (top), Job 807 (middle), 632510 (bottom). The
bottom panel's cusp is fixed by the existing post-fix run in fix_632510/, but
that stops at 16.37 d, and the middle panel's run is unreproducible as written.
This builds both.

632894 -- 807's PHYSICS on the current code.
  res807.in itself cannot be rerun. Three reasons:
    * `shear_mod 24` is not a recognised key. read_inputfile has only
      case('rigid') (main_LH.f90:2695) -- the key was renamed at some point --
      so 807 silently ran at the default rigid = 32.04 GPa (m_const.f90:4)
      rather than the 24 it asked for. Now written as `rigid 24`.
    * `Sw_fwid` is absent, so it defaulted to 1.0 (main_LH.f90:180) against the
      7.4e-9 in use. 807 predates the wellbore model. Now set explicitly.
    * it uses oct_clean_300.txt, which is in the OLDER injection format:
        nwell / iwell jwell npoint / rates / times
      against what input_well now reads (m_diffusion.f90:49):
        nwell / npoint / times / iwell jwell / rates
      The current reader would take "300 300 499" and use 300 as npoint, then
      read the RATE line as times. Silently wrong, not an error.

  INJECTION: june_clean.txt, not a converted oct_clean_300.txt. They are the
  same November 2012 record sampled differently -- identical 31.39 d span,
  identical 8.7443e-3 peak, total volume per width 6085.9 vs 6090.4, a 0.07%
  difference -- but they differ pointwise by up to 9.8e-4 (11% of peak), and the
  well index differs (300,300 vs 301,301). Using june_clean for both simulated
  panels makes them directly comparable and matches the bottom panel, which is
  what the figure is for. Converting the old file would preserve 807's exact
  sampling at the cost of that comparability.

  Sw_fwid = 7.4e-9, NOT the 7.4e-7 used for the constant-rate runs. That larger
  value was needed because those inject 4.2e-3 from t = 0, and at kp = 1e-15 the
  Peaceman well is stiff (gamma = 0.977, pressure ramping at q/Sw_fwid).
  june_clean starts at 3.3e-4, a 13x gentler ramp, and 632510 has already run at
  kp = 1e-15 with Sw_fwid 7.4e-9 on this same injection without trouble. 7.4e-9
  also keeps both panels consistent.

632895 -- 632510 extended to 17 d.
  The post-fix run used nstep 40000 and stopped at 16.37 d having exhausted it,
  not at tmax. Deck is otherwise identical, including Sw_fwid 7.4e-9, so the new
  run continues the same physics rather than changing it -- important, since the
  16.37 d panel is already made and this should extend it, not replace it with
  something different.

Both: imax/jmax 601 (807 used 600, which has no centre cell), tmax 17 d to match
the figure, dtout 0.0002 for finer curves, nstep 400000.

Domain is safe for both -- the far field stays at kpmin = 1e-15, so
D = kp/(eta*phi*beta) is 4.99e-3 and 5.61e-3 m^2/s and the 17 d diffusion length
is 171 and 181 m against a 1502 m half-domain.

Usage:  python build_fig6.py [--write]
"""
import argparse
import math
import sys
from pathlib import Path

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
DAYS = 17.0
JOBS = {632894: dict(parent="807",
                     note="Job 807's physics on the current code"),
        632895: dict(parent="632510",
                     note="632510 post-fix, extended past nstep to reach 17 d")}


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
    tmax_yr = DAYS / 365.0

    common = {
        "imax": "601", "jmax": "601",
        "tmax": f"{tmax_yr:.8f}",
        "dtout": "0.0002",
        "nstep": "400000",
        "injectionfromfile": "T",
        "injection_file": '"june_clean.txt"',
        "rigid": "24",          # replaces 807's unrecognised shear_mod
        "Sw_fwid": "7.4e-9",
        "rw": "8.9e-2",
    }
    drop = {"shear_mod", "injection", "qinj"}

    made = []
    for new, spec in JOBS.items():
        src = IN / f"res{spec['parent']}.in"
        pairs = read_deck(src)
        base = dict(pairs)
        eta = float(base["eta"].replace("d", "e"))
        phi, beta = float(base["phi"]), float(base["beta"].replace("d", "e"))
        kp = float(base["kp"].replace("d", "e"))
        D = kp / (eta * phi * beta)
        L = math.sqrt(4 * D * DAYS * 86400)
        half = 0.005 * 601 / 2 * 1000
        assert L / half < 0.8, f"domain: L/half = {L/half:.2f}"
        print(f"res{new}.in  <- res{spec['parent']}.in   {spec['note']}")
        print(f"   kp {kp:.0e}  kpmax {base['kpmax']}  "
              f"parameterfromfile {base['parameterfromfile']}  "
              f"D_far {D:.3e} m^2/s  L({DAYS:.0f}d) {L:.0f} m / {half:.0f} m "
              f"= {L/half:.2f}")
        if not a.write:
            continue

        hdr = [
            f"! Fig 6 remake -- {spec['note']}.",
            f"! Built from res{spec['parent']}.in.",
            "!",
        ]
        if spec["parent"] == "807":
            hdr += [
                "! res807.in cannot be rerun as written:",
                "!  * shear_mod 24 is not a recognised key -- read_inputfile has only",
                "!    case('rigid') (main_LH.f90:2695), the key having been renamed --",
                "!    so 807 silently ran at the default rigid = 32.04 GPa",
                "!    (m_const.f90:4), not the 24 it asked for. Now rigid 24.",
                "!  * Sw_fwid was absent and defaults to 1.0 (main_LH.f90:180)",
                "!    against the 7.4e-9 in use. 807 predates the wellbore model.",
                "!  * it used oct_clean_300.txt, in the OLDER injection format",
                "!    (nwell / iwell jwell npoint / rates / times) against what",
                "!    input_well now reads (nwell / npoint / times / iwell jwell /",
                "!    rates). The current reader would take 300 as npoint and then",
                "!    read the rate line as times -- silently wrong, not an error.",
                "!",
                "! INJECTION is june_clean.txt rather than a converted",
                "! oct_clean_300.txt. Same November 2012 record sampled differently:",
                "! identical 31.39 d span and 8.7443e-3 peak, volume per width",
                "! 6085.9 vs 6090.4 (0.07%), but pointwise differences up to 9.8e-4,",
                "! 11% of peak, and a different well index (300,300 vs 301,301).",
                "! june_clean makes this panel directly comparable to the 632510",
                "! panel, which is the point of the figure.",
                "!",
                "! Sw_fwid 7.4e-9, NOT the 7.4e-7 used for the constant-rate runs.",
                "! That was needed because those inject 4.2e-3 from t=0 and at",
                "! kp = 1e-15 the Peaceman well is stiff. june_clean starts at",
                "! 3.3e-4, a 13x gentler ramp, and 632510 already runs at kp = 1e-15",
                "! with Sw_fwid 7.4e-9 on this injection without trouble.",
            ]
        else:
            hdr += [
                "! The post-fix run of 632510 (in fix_632510/) stopped at 16.37 d",
                "! having exhausted nstep 40000, NOT at tmax. Only nstep, tmax,",
                "! dtout and the filenumber change here, so this extends the same",
                "! physics rather than replacing it -- the 16.37 d panel is already",
                "! made from that run and this should continue it.",
            ]
        hdr += [
            "!",
            f"! Domain: far field stays at kpmin, D = {D:.3e} m^2/s, so the",
            f"! {DAYS:.0f} d diffusion length is {L:.0f} m against a {half:.0f} m",
            f"! half-domain ({100*L/half:.0f}%).",
            "!",
            "! imax/jmax 601 for an exact centre cell (807 used 600, which has none).",
        ]
        lines = [f"filenumber {new}"] + hdr
        seen = set()
        for k, v in pairs:
            if k == "filenumber" or k in drop:
                continue
            lines.append(f"{k} {common.get(k, v)}")
            seen.add(k)
        for k, v in common.items():
            if k not in seen:
                lines.append(f"{k} {v}")
        (IN / f"res{new}.in").write_text("\n".join(lines) + "\n")
        made.append((new, spec["parent"]))

    if not a.write:
        print("\ndry run -- nothing written. re-run with --write")
        return

    print("\nverifying")
    bad = 0
    PHYS = ("a", "b", "dc", "f0", "muinit", "sigmainit", "kp", "kpmin", "kpmax",
            "kL", "kT", "phi", "beta", "eta", "ds", "problem", "permev",
            "limitsigma", "minsig", "parameterfromfile", "backslip", "vpl",
            "velinit", "velmin", "pfinit", "pressurediffusion",
            "bcl", "bcr", "bct", "bcb", "pbcl", "pbcr", "pbct", "pbcb")
    for new, parent in made:
        b = dict(read_deck(IN / f"res{parent}.in"))
        d = dict(read_deck(IN / f"res{new}.in"))
        ck = {
            "physics == parent": all(d.get(k) == b.get(k) for k in PHYS
                                     if k in b),
            "rigid 24": d.get("rigid") == "24",
            "Sw_fwid 7.4e-9": d.get("Sw_fwid") == "7.4e-9",
            "no shear_mod": "shear_mod" not in d,
            # strip quotes: the deck stores it as "june_clean.txt" with them
            "june_clean": d.get("injection_file", "").strip('"') == "june_clean.txt",
            "601": d.get("imax") == "601" and d.get("jmax") == "601",
            "17 d": abs(float(d["tmax"]) * 365 - DAYS) < 1e-4,
            "nstep raised": int(d["nstep"]) > int(b.get("nstep", 0)),
        }
        ok = all(ck.values())
        bad += not ok
        print(f"  res{new}.in: " + ("OK" if ok else "PROBLEM "
                                    + str([k for k, v in ck.items() if not v])))
    print(f"\n  {len(made)-bad}/{len(made)} clean")
    sys.exit(0 if bad == 0 else 1)


if __name__ == "__main__":
    main()
