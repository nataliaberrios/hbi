#!/usr/bin/env python3
"""Stage 2: two-zone initial permeability, consistent bounds, full 16-cell slice.

The slice is parent x sigmabar_0 x fluid combination x permev = 2x2x2x2 = 16
cells. Four already exist and are NOT rebuilt:

    632522  1807, sigma 30.0, fluid A, permev T
    632523  1807, sigma 30.0, fluid A, permev F
    632524  1808, sigma 30.0, fluid A, permev T
    632525  1808, sigma 30.0, fluid A, permev F

This builds the remaining 12, from those four as parents.

Axes:
  parent   1807 -> kpmax 2.5e-13, map perm_2zone_601_ds5_kmax2.5e-13.txt (250x)
           1808 -> kpmax 5e-14,   map perm_2zone_601_ds5_kmax5e-14.txt   (50x)
           In both, the map's near-well disc IS kpmax and its background IS
           kp = kpmin = 1e-15. That consistency is the whole point: across 112
           archived runs, the 35 obeying it had 29 reach dc, the 77 that did not
           had 2.
  sigma    30.0 (res1807/res1808's round number) or 27.99 (derived from the
           stress data and verified against setup_model.m). dp_crit =
           sigma*(1 - muinit/f0) is 11.50 vs 10.73 MPa, against a peak measured
           overpressure of 10.92 -- so sigma decides whether the fault can fail
           at the observed pressure at all.
  fluid    A = the parent's own eta and beta.
           B = eta 1.27e-4 with beta x 7.007874, holding D = k/(eta*phi*beta)
           fixed. NOT a null change: T = 2*pi*k/eta (m_diffusion.f90:669)
           depends on eta alone, so the wellbore coupling gamma shifts.
  permev   enhancement on or off. Off is a real alternative model, not a
           degenerate case -- with a two-zone initial field the near-well
           conduit already exists, so permev F asks whether slip needs to
           extend it at all.

tmax is 5 d here, against 8 d for the four existing runs. Both cover the 0-5 d
comparison window, which is what every score uses.

Usage:  python build_stage2.py [--check]
"""
import argparse
import sys
from pathlib import Path

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")

ETA_OLD, ETA_NEW = 8.9e-4, 1.27e-4
RATIO = ETA_OLD / ETA_NEW               # 7.007874015748031
TMAX_5D = "0.01369863"
DTOUT = "0.0002"

# parent run to copy from -> (which res180x it descends from, its beta)
BASE = {1807: (632522, "2.25e-8"), 1808: (632524, "2.005e-8")}
SIGMAS = ["30.0", "27.99"]
FLUIDS = ["A", "B"]
PERMEV = ["T", "F"]

# the four cells that already exist and must not be rebuilt
EXISTING = {
    (1807, "30.0", "A", "T"): 632522, (1807, "30.0", "A", "F"): 632523,
    (1808, "30.0", "A", "T"): 632524, (1808, "30.0", "A", "F"): 632525,
}

RUNS = []
_n = 632810
for parent in (1807, 1808):
    for sig in SIGMAS:
        for fl in FLUIDS:
            for pev in PERMEV:
                if (parent, sig, fl, pev) in EXISTING:
                    continue
                RUNS.append((_n, parent, sig, fl, pev))
                _n += 1


def read_deck(p):
    out = []
    for line in Path(p).read_text().splitlines():
        if line.startswith("!") or not line.strip():
            continue
        w = line.split()
        if len(w) >= 2:
            out.append((w[0], " ".join(w[1:])))
    return out


def build(num, parent, sig, fl, pev):
    src, beta0 = BASE[parent]
    pairs = read_deck(IN / f"res{src}.in")
    d = dict(pairs)
    eta = "0.89e-3" if fl == "A" else "1.27e-4"
    beta = beta0 if fl == "A" else f"{float(beta0) * RATIO:.4e}"
    mu, f0, phi = float(d["muinit"]), float(d["f0"]), float(d["phi"])
    dpc = float(sig) * (1 - mu / f0)
    kpmax, kpmin = d["kpmax"], d["kpmin"]

    changes = {"filenumber": str(num), "sigmainit": sig, "eta": eta,
               "beta": beta, "permev": pev, "tmax": TMAX_5D, "dtout": DTOUT}
    hdr = [
        f"! STAGE 2 -- two-zone initial permeability, consistent bounds.",
        f"! Cell: parent {parent} | sigmainit {sig} | fluid {fl} | permev {pev}",
        f"! Built from res{src}.in, which is the (30.0, A, T/F) cell of this slice.",
        "!",
        f"! MAP: {d['parameter_file']}",
        f"!   near-well 150 m disc = kpmax = {kpmax}",
        f"!   background           = kp = kpmin = {kpmin}",
        "!   The map's two values ARE the evolution bounds -- nothing else appears.",
        "!   Of 112 archived runs with enhancement and a nonuniform map, the 35",
        "!   obeying this had 29 reach dc; the 77 that did not had 2.",
        "!",
        f"! permev {pev}: "
        + ("slip drives the background up toward kpmax, extending the conduit."
           if pev == "T" else
           "enhancement OFF -- the conduit is the initial condition only. With a"),
        "!" if pev == "T" else
        "!   two-zone field this is a real alternative model, not a degenerate case:",
        "!" if pev == "T" else
        "!   it asks whether slip needs to extend the conduit at all.",
        "!",
        f"! fluid {fl}: eta {eta}, beta {beta}"
        + ("  (the parent's own values)" if fl == "A" else
           f"  (beta x {RATIO:.6f} so D is unchanged from A)"),
    ]
    if fl == "B":
        hdr += [
            "!   Holding D fixed is NOT a null change: the Peaceman well index",
            "!   T = 2*pi*k/eta (m_diffusion.f90:669) depends on eta alone, so the",
            "!   wellbore coupling gamma = (Sw_fwid/h)/((Sw_fwid/h)+T) shifts.",
            f"!   Also phi*beta = {phi*float(beta):.3e}, ~3.9x above what phi<=0.02 and",
            "!   beta<=2e-8 allow. Deliberate, for diffusivity consistency.",
        ]
    hdr += [
        "!",
        f"! sigmainit {sig} MPa"
        + ("  (res1807/1808's round number, NOT measurement-derived)" if sig == "30.0"
           else "  (DERIVED: sv 100, sH 160, dip 10 deg, p_pore 73.82 -> 27.99)"),
        f"! tau0 = {mu*float(sig):.2f} MPa (Taiyi 15.0)   "
        f"dp_crit = {dpc:.2f} MPa",
        f"! Peak MEASURED downhole overpressure in 5 d is 10.92 MPa, so this cell "
        + ("CAN" if dpc <= 10.92 else "CANNOT"),
        "! fail at the observed pressure"
        + ("." if dpc <= 10.92 else " -- it can only slip by over-pressurising."),
        "!",
        f"! tmax {TMAX_5D} yr = 5.00 d. res{src}.in ran 8 d; both cover the 0-5 d",
        "! comparison window that every score uses.",
    ]
    lines = [f"filenumber {num}"] + hdr
    for k, v in pairs:
        if k == "filenumber":
            continue
        lines.append(f"{k} {changes.get(k, v)}")
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    print(f"Stage 2 slice: 2 parents x 2 sigma x 2 fluid x 2 permev = 16 cells")
    print(f"  {len(EXISTING)} already exist: "
          + ", ".join(str(v) for v in EXISTING.values()))
    print(f"  {len(RUNS)} to build\n")
    print(f"{'run':>7s} {'parent':>7s} {'sigma':>7s} {'fluid':>6s} {'permev':>7s} "
          f"{'eta':>9s} {'beta':>11s} {'tau0':>6s} {'dp_crit':>8s} {'can slip?':>10s}")
    print("-" * 96)
    for num, parent, sig, fl, pev in RUNS:
        src, beta0 = BASE[parent]
        d = dict(read_deck(IN / f"res{src}.in"))
        eta = "0.89e-3" if fl == "A" else "1.27e-4"
        beta = beta0 if fl == "A" else f"{float(beta0)*RATIO:.4e}"
        mu = float(d["muinit"])
        dpc = float(sig) * (1 - mu / float(d["f0"]))
        print(f"{num:>7d} {parent:>7d} {sig:>7s} {fl:>6s} {pev:>7s} {eta:>9s} "
              f"{float(beta):>11.4e} {mu*float(sig):>6.2f} {dpc:>8.2f} "
              f"{('yes' if dpc <= 10.92 else 'NO'):>10s}")
        if not a.check:
            (IN / f"res{num}.in").write_text(build(num, parent, sig, fl, pev))
    if a.check:
        print("\n--check: nothing written")
        return

    print(f"\nwrote {len(RUNS)} decks into {IN}")
    print("\nverifying each against the res6325xx cell it was built from:")
    ok = True
    for num, parent, sig, fl, pev in RUNS:
        src, _ = BASE[parent]
        new = dict(read_deck(IN / f"res{num}.in"))
        old = dict(read_deck(IN / f"res{src}.in"))
        diff = {k for k in set(new) | set(old) if new.get(k) != old.get(k)}
        allowed = {"filenumber", "sigmainit", "eta", "beta", "permev",
                   "tmax", "dtout"}
        extra = diff - allowed
        # the map and the bounds must be untouched -- that is the invariant
        for k in ("parameter_file", "kpmax", "kpmin", "kp", "muinit", "f0",
                  "imax", "jmax", "ds", "injection_file", "rw", "Sw_fwid"):
            if new.get(k) != old.get(k):
                extra.add(k)
        if extra:
            ok = False
        print(f"  res{num}.in vs res{src}.in: {sorted(diff)}"
              + (f"   UNEXPECTED: {sorted(extra)}" if extra else "   OK"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
