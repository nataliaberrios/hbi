#!/usr/bin/env python3
"""Stage 1: rerun 1807/1808 physics on fixed HBI, uniform initial permeability.

Builds eight decks from res1807.in and res1808.in, changing ONLY what has to
change, and writes them into the shared grid_search_inputs directory so they can
be diffed against every other deck.

Three axes: parent (1807 / 1808) x fluid combination (A / B) x sigmabar_0
(30.0 / 27.99). Every run therefore has a sigmainit twin differing in that one
value and nothing else.

Two combinations per parent:
  A  the parent's own fluid and storage (eta 0.89e-3, beta as-is)
  B  the correct reservoir viscosity eta 1.27e-4 with beta scaled by
     0.89e-3 / 1.27e-4 = 7.007874 so that D = k/(eta*phi*beta) is UNCHANGED

Note what B does and does not hold fixed. eta and beta do not enter the solver
only through their product:

    cdiff = kpG/(eta*beta*phi)      m_diffusion.f90:443-445   invariant
    str   = beta*phi                m_diffusion.f90:444       beta alone, x7.008
    T     = 2*pi*kpG/eta/(log(0.2*ds/rw)+skin)
                                    m_diffusion.f90:669       eta alone,  x7.008
    TT    = T/(str*ds^2)                                      invariant
    gamma = (Sw_fwid/h)/((Sw_fwid/h)+T)
                                    m_diffusion.f90:671       NOT invariant

So B preserves the diffusion coefficient exactly but changes the Peaceman
wellbore coupling: gamma 0.977 -> 0.858 at the initial kp, and 0.145 -> 0.024
once enhancement drives kp to kpmax. A-vs-B therefore measures the well-coupling
consequence of the viscosity correction; it is not a null test.

Deliberate deviations from the parent decks, and why:
  injection_file  april_clean_fixed2.txt -> june_clean.txt
                  The pressure and seismicity data being matched are the June
                  2012 stage. Without this the run is not comparable to the
                  target at all. This means these are NOT reproductions of
                  1807/1808 -- same physics, different stimulation.
  imax/jmax       600 -> 601
                  An odd grid puts a cell exactly at the centre, so the injector
                  is the single well-defined cell (301,301). On 600 it sits half
                  a cell (2.5 m) off centre.
  tmax            0.084 yr (30.66 d) -> 0.01369863 yr (5.00 d)
                  Both wellhead spikes (43.15 MPa at 0.95 d, 44.73 MPa at 1.58 d)
                  are inside the first flowing block, which ends at 1.576 d; the
                  well is then shut until 3.56 d. 5 d covers both spikes, the
                  shut-in recovery, and 1.4 d of cycle 2.
  dtout           0.001 yr (8.8 h) -> 0.0002 yr (1.75 h)
                  0.001 would put only ~14 frames in 5 d, too coarse to resolve
                  spikes 0.6 d apart.

Everything else is inherited verbatim, including rw 8.9e-2 and Sw_fwid 7.4e-9 --
which are ALREADY in 1807/1808 and are not new parameters, contrary to
expectation. Sw_fwid especially must never be dropped: its in-code default is
1.0 (main_LH.f90:180), eight orders of magnitude away.

Usage:  python build_stage1.py           # write the decks
        python build_stage1.py --check    # verify only, write nothing
"""
import argparse
import sys
from pathlib import Path

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")

ETA_OLD, ETA_NEW = 8.9e-4, 1.27e-4
RATIO = ETA_OLD / ETA_NEW               # 7.007874015748031

TMAX_5D = "0.01369863"                  # 5.00 d at HBI's 365 d year
DTOUT = "0.0002"
GRID = "601"
INJ = "june_clean.txt"

# sigmabar_0 is an AXIS, not a fixed choice. 30.0 is what res1807/res1808 use --
# a round number, not a measurement. 27.99 is derived from the stress data:
# sigma_n = sv*cos^2(dip) + sH*sin^2(dip) = 100*cos^2(10) + 160*sin^2(10) = 101.81,
# minus p_pore 73.82 -> 27.99 MPa. Since dp_crit = sigmabar_0*(1 - muinit/f0)
# scales with sigmabar_0, the choice decides whether the fault can fail at the
# measured overpressure at all: at muinit 0.37, dp_crit is 11.50 MPa at 30.0
# (above the 10.92 MPa measured peak, so it CANNOT fail) but 10.73 MPa at 27.99
# (below it, so it can). Every run therefore gets a twin differing ONLY in
# sigmainit.
SIGMAS = [("30.0", "res1807/res1808's own value, not measurement-derived"),
          ("27.99", "derived from sv 100, sH 160, dip 10 deg, p_pore 73.82")]

# (filenumber, parent, label, eta, beta) -- sigma twin appended below
_BASE = [
    (1807, "A", "0.89e-3", "2.25e-8"),
    (1807, "B", "1.27e-4", f"{2.25e-8 * RATIO:.4e}"),
    (1808, "A", "0.89e-3", "2.005e-8"),
    (1808, "B", "1.27e-4", f"{2.005e-8 * RATIO:.4e}"),
]
RUNS = []
_n = 632800
for _parent, _label, _eta, _beta in _BASE:
    for _sig, _ in SIGMAS:
        RUNS.append((_n, _parent, _label, _eta, _beta, _sig))
        _n += 1


def read_deck(p):
    """Ordered (key, value) pairs, skipping comments and blanks."""
    out = []
    for line in Path(p).read_text().splitlines():
        if line.startswith("!") or not line.strip():
            continue
        w = line.split()
        if len(w) >= 2:
            out.append((w[0], " ".join(w[1:])))
    return out


def build(num, parent, label, eta, beta, sig):
    pairs = read_deck(IN / f"res{parent}.in")
    d = dict(pairs)
    kpmax, kpmin, kp = d["kpmax"], d["kpmin"], d["kp"]
    phi = float(d["phi"])
    d_far = float(kpmin) / (float(eta.replace("d", "e")) * phi * float(beta))

    changes = {
        "filenumber": str(num),
        "imax": GRID, "jmax": GRID,
        "tmax": TMAX_5D, "dtout": DTOUT,
        "injection_file": f'"{INJ}"',
        "eta": eta, "beta": beta, "sigmainit": sig,
    }
    hdr = [
        f"! STAGE 1 -- res{parent}.in physics rerun on FIXED HBI, uniform initial perm.",
        f"! Combination {label}: "
        + ("the parent's own fluid and storage."
           if label == "A" else
           f"correct reservoir viscosity eta {eta} with beta scaled by "
           f"{RATIO:.6f}"),
        "!" if label == "A" else
        f"! so that D = k/(eta*phi*beta) is UNCHANGED from combination A.",
        "!",
        f"! permev T with UNIFORM kp = kpmin = {kp}; slip grows it toward",
        f"! kpmax = {kpmax}. No permeability map -- parameterfromfile F. The",
        "! near-well high-permeability zone is therefore created DYNAMICALLY by slip,",
        "! which is what distinguishes Stage 1 from Stage 2.",
        "!",
        f"! phi {d['phi']} x beta {beta} ;  D(far, at kpmin) = {d_far:.4e} m^2/s",
        f"! sigmainit {sig} MPa"
        + ("  (res1807/1808's own value; a round number, NOT measurement-derived)"
           if sig == "30.0" else
           "  (DERIVED from sv 100, sH 160, dip 10 deg, p_pore 73.82)"),
        f"! tau0 = muinit x sigmainit = {float(d['muinit'])*float(sig):.2f} MPa"
        f"  (Taiyi 15.0, so understressed)",
        f"! dp_crit = sigmainit(1 - muinit/f0) = "
        f"{float(sig)*(1-float(d['muinit'])/float(d['f0'])):.2f} MPa; the peak MEASURED"
        " downhole",
        "! overpressure in the first 5 d is 10.92 MPa, so this run "
        + ("CANNOT fail without over-pressurising."
           if float(sig)*(1-float(d['muinit'])/float(d['f0'])) > 10.92
           else "CAN fail at the observed pressure."),
        "! Its sigmainit twin differs in sigmainit and nothing else.",
        "!",
        "! DELIBERATE deviations from the parent deck:",
        f"!   injection_file -> {INJ}. The data being matched are the June 2012",
        "!     stage; the parent used april_clean_fixed2.txt. This is NOT a",
        "!     reproduction of the parent -- same physics, different stimulation.",
        f"!   imax/jmax 600 -> {GRID}, so the injector is the exact centre cell (301,301).",
        f"!   tmax -> {TMAX_5D} yr = 5.00 d. Both wellhead spikes (0.95 d, 1.58 d) are",
        "!     inside the first flowing block, which ends at 1.576 d.",
        f"!   dtout -> {DTOUT} yr (1.75 h); 0.001 gave only ~14 frames in 5 d.",
        "!",
        "! rw 8.9e-2 and Sw_fwid 7.4e-9 are inherited unchanged -- they were ALREADY",
        "! in the parent deck and are not new parameters.",
    ]
    if label == "B":
        hdr += [
            "!",
            "! CAVEAT: holding D fixed does NOT reproduce combination A. The Peaceman",
            "! well index T = 2*pi*k/eta (m_diffusion.f90:669) depends on eta ALONE, so",
            "! gamma = (Sw_fwid/h)/((Sw_fwid/h)+T) (m_diffusion.f90:671) is not preserved:",
            "! 0.977 -> 0.858 at the initial kp, and 0.145 -> 0.024 once kp reaches kpmax.",
            f"! Also phi*beta = {phi*float(beta):.3e} is ~3.9x above the maximum the stated",
            "! physical ranges allow (phi<=0.02, beta<=2e-8). Deliberate, for diffusivity",
            "! consistency; see docs/figs/cooper_basin_calibration.",
        ]

    # filenumber MUST stay on line 1: the harness reads it with awk 'NR==1 {print $2}'
    lines = [f"filenumber {num}"] + hdr
    for k, v in pairs:
        if k == "filenumber":
            continue
        lines.append(f"{k} {changes.get(k, v)}")
    return "\n".join(lines) + "\n", changes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    print(f"eta ratio {ETA_OLD:.3e} / {ETA_NEW:.3e} = {RATIO:.9f}\n")
    print(f"{'run':>7s} {'parent':>7s} {'cmb':>4s} {'sigma':>7s} {'eta':>9s} "
          f"{'beta':>11s} {'phi*beta':>10s} {'tau0':>8s} {'dp_crit':>8s} "
          f"{'vs 10.92':>9s} {'D_far m2/s':>11s}")
    print("-" * 108)
    for num, parent, label, eta, beta, sig in RUNS:
        text, _ = build(num, parent, label, eta, beta, sig)
        d = dict(read_deck(IN / f"res{parent}.in"))
        phi, mu, f0 = float(d["phi"]), float(d["muinit"]), float(d["f0"])
        dfar = float(d["kpmin"]) / (float(eta.replace("d", "e")) * phi * float(beta))
        dpc = float(sig) * (1 - mu / f0)
        print(f"{num:>7d} {parent:>7d} {label:>4s} {sig:>7s} {eta:>9s} "
              f"{float(beta):>11.4e} {phi*float(beta):>10.3e} {mu*float(sig):>8.2f} "
              f"{dpc:>8.2f} {'CANNOT' if dpc > 10.92 else 'can fail':>9s} {dfar:>11.4e}")
        if not a.check:
            (IN / f"res{num}.in").write_text(text)
    if a.check:
        print("\n--check: nothing written")
        return
    print(f"\nwrote {len(RUNS)} decks into {IN}")

    # every deck must differ from its parent ONLY in the intended keys
    print("\nverifying each deck against its parent (non-comment lines only):")
    ok = True
    for num, parent, label, eta, beta, sig in RUNS:
        new = dict(read_deck(IN / f"res{num}.in"))
        old = dict(read_deck(IN / f"res{parent}.in"))
        diff = {k for k in set(new) | set(old) if new.get(k) != old.get(k)}
        expect = {"filenumber", "imax", "jmax", "tmax", "dtout",
                  "injection_file", "eta", "beta", "sigmainit"}
        if sig == old.get("sigmainit"):
            diff.discard("sigmainit")
            expect.discard("sigmainit")
        if label == "A":
            expect -= {"eta"}          # A keeps the parent's eta
            if new["eta"] == old["eta"]:
                diff.discard("eta")
        extra = diff - expect
        missing = {"filenumber", "imax", "jmax", "tmax", "dtout",
                   "injection_file"} - diff
        status = "OK" if not extra and not missing else "PROBLEM"
        if status != "OK":
            ok = False
        print(f"  res{num}.in vs res{parent}.in: changed {sorted(diff)}  -> {status}")
        if extra:
            print(f"     UNEXPECTED changes: {sorted(extra)}")
        if missing:
            print(f"     MISSING expected changes: {sorted(missing)}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
