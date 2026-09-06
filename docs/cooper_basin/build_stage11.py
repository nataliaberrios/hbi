#!/usr/bin/env python3
"""Stage 11: 18 d extensions of the five plotted runs, plus the a-b sweep.

TWO SETS, built together because they answer independent questions.

SET A -- 18 d versions of the five runs in compare_story5/slip_poster_story5.
The observed slip file has 3,5,...,17 d, but every scored run stops at 5 d, so
the poster-style figure carries two curves against the original's eight, and
both of those nearly coincide (the observed 3 d and 5 d columns are
byte-identical, and simulated slip stalls after the first flowing period --
632901 goes 0.487 -> 0.488 cm). 18 d clears 17 with margin, which matters: the
last run to target 17 d ended at 16.999999 d and a naive tolerance dropped its
most important curve.

It also tests something the 5 d window cannot. lambda_obs falls from 0.1866 over
0-5 d to 0.1412 over 0-8 d, so if the front stops tracking sqrt(t) past 5 d the
match is a 5 d artifact. That is worth more than the figure.

  632910 <- 632901   the current best: front 1.06x, wellhead peak +1.6%
  632911 <- 632896   2.6x more slip, wellhead +19.1%
  632912 <- 632875   front matches, wellhead +80%
  632913 <- 632812   matches the observed AMPLITUDE to 3%, front 0.62x
  632914 <- 632881   Wang & Dunham's own parameters -- 14 d, NOT 18, see below

DOMAIN, checked before building rather than after. For the four ds = 5 m runs the
far field stays at kpmin = 1e-15, so D = 4.994e-3 m^2/s and L(18 d) = 176 m
against a 1502 m half-domain, L/half = 0.12. Ample.

632881 is different and is the one case that fails. It is Wang & Dunham's
config: permev F at kp = 4e-13 with ds = 20 m, so D = 4.494 m^2/s -- three
orders above the others -- and L/half reaches 0.80 at 14.9 d:

     14 d  L 4663 m  L/half 0.78  OK
     15 d  L 4827 m  L/half 0.80  too big
     18 d  L 5288 m  L/half 0.88  too big

Reaching 18 d safely would need imax >= 700, 1.4x the cells, which changes the
discretisation and makes it a different run rather than a longer one. So 632914
is 14 d. That still covers observed times 3,5,7,9,11,13 -- six of eight -- and
the plotting script trims and reports rather than extrapolating.

SET B -- the a-b sweep, off 632901. Every pressure route to slip amplitude is
now closed. Stage 9: amplitude crosses the observed value near kpmax 8e-13,
where interpolating between 632875 (+80%) and 632896 (+19.1%) puts the wellhead
at +35-40%, outside the band. Stage 10: raising tau_0 does nothing because the
plateau tracks dp_crit down with it. Friction is the only untouched knob --
a = 0.015 and b = 0.012 in all 75 scored runs.

The fault never reaches failure; it creeps at
v = vref*exp[(tau/sigmabar - f0)/(a-b)], so a-b sets how far below f0 it keeps
creeping before it stalls. Measured on 632901 that window is d(tau/sigmabar) =
0.0278, i.e. 9.3 e-folds of velocity, and the 9.3 is set by ln(velocity ratio)
rather than by a-b -- so slip should be roughly LINEAR in a-b, needing about
0.017 for the observed 2.81 cm. a-b stays POSITIVE: the fault remains
velocity-strengthening, which also means no nucleation-size constraint against
ds = 5 m, since h* only exists for b > a.

  632915  a 0.015  b 0.005  a-b 0.010   isolates a-b with a unchanged
  632916  a 0.022  b 0.005  a-b 0.017   targets the observed amplitude
  632917  a 0.022  b 0.012  a-b 0.010   CONTROL: same a-b as 632915, different a

632915 vs 632917 is the run that earns the parameterisation. If they agree, only
the difference matters. If they diverge, `a` acts on its own through the direct
effect and the sweep has to become two-dimensional.

Set B is also 18 d, not 5 d. A longer run is scored on the same 0-5 d window
(WINDOW_D in score_grid.py, which already handles the 8 d runs 632522-525), so
there is no comparability cost, and a run that turns out well needs no rerun for
the figure. RISK, stated in advance: faster creep may stiffen the integrator.
632902 already died at 4.51 d from stiffness at muinit 0.410 -- below the 4.75 d
usability gate. These keep 632901's muinit 0.397, which completed, but a stall
before 4.75 d would waste the run. nstep is 400000 either way, so a stall is the
only failure mode, not an nstep exhaustion like the first 632510 rerun.

Usage:  python build_stage11.py [--write]
"""
import argparse
import math
import sys
from pathlib import Path

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
# (new, parent, days, extra overrides)
JOBS = [
    (632910, 632901, 18.0, {}),
    (632911, 632896, 18.0, {}),
    (632912, 632875, 18.0, {}),
    (632913, 632812, 18.0, {}),
    (632914, 632881, 14.0, {}),          # domain-limited, see the docstring
    (632915, 632901, 18.0, {"a": "0.015", "b": "0.005"}),
    (632916, 632901, 18.0, {"a": "0.022", "b": "0.005"}),
    (632917, 632901, 18.0, {"a": "0.022", "b": "0.012"}),
]
NSTEP = "400000"
SAFE = 0.8


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


def extent(base, days, lam_km_sqrtd=None):
    """(L_est, half, note) -- how far pressure actually gets, not just D_far.

    THE CHECK THIS REPLACES WAS WRONG FOR ENHANCEMENT RUNS. It used only
    L = sqrt(4*D_far*t) with D_far from kp = kpmin, and reported L/half = 0.12
    for the ds = 5 m 18 d runs. That is the right question for a run whose
    permeability never changes, but with permev T the high-k zone is CARRIED
    OUTWARD by slip -- 632901's map disc is 150 m while its pressure plateau and
    slip front both reach 380 m at 5 d. The extent is therefore the enhanced
    edge PLUS diffusion beyond it, and the honest numbers for Stage 11 are
    0.48-0.67, not 0.12.

    For permev T, the enhanced edge is extrapolated from the run's own measured
    0-5 d front fit, R[km] = lam*sqrt(t[d]).
    """
    D = ff(base["kp"]) / (ff(base["eta"]) * ff(base["phi"]) * ff(base["beta"]))
    L_diff = math.sqrt(4 * D * days * 86400)
    half = int(base["imax"]) * ff(base["ds"]) * 1000 / 2
    on = str(base.get("permev", "F")).upper().startswith("T")
    if on and lam_km_sqrtd:
        R = lam_km_sqrtd * math.sqrt(days) * 1000.0
        return R + L_diff, half, f"enhanced edge {R:.0f}m + diff {L_diff:.0f}m"
    return L_diff, half, f"fixed perm, D_far {D:.3e}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()

    print(f"{'new':>7s} {'parent':>7s} {'days':>5s} {'D_far':>10s} {'L':>7s} "
          f"{'half':>6s} {'L/half':>7s} {'a-b':>7s} {'change':>22s}")
    rows, bad = [], 0
    for new, par, days, extra in JOBS:
        pairs = read_deck(IN / f"res{par}.in")
        b = dict(pairs)
        D = ff(b["kp"]) / (ff(b["eta"]) * ff(b["phi"]) * ff(b["beta"]))
        L = math.sqrt(4 * D * days * 86400)
        half = int(b["imax"]) * ff(b["ds"]) * 1000 / 2
        amb = ff(extra.get("a", b["a"])) - ff(extra.get("b", b["b"]))
        ok = L / half < SAFE
        bad += not ok
        chg = ", ".join(f"{k} {v}" for k, v in extra.items()) or "duration only"
        print(f"{new:>7d} {par:>7d} {days:>5.0f} {D:>10.3e} {L:>6.0f}m "
              f"{half:>5.0f}m {L/half:>7.2f} {amb:>+7.3f} {chg:>22s}"
              + ("" if ok else "   <-- DOMAIN FAIL"))
        rows.append((new, par, days, extra, pairs, b))
    if bad:
        sys.exit(f"\n{bad} deck(s) would exceed L/half = {SAFE}. Not writing.")
    print(f"\nall {len(rows)} within L/half < {SAFE}")

    if not a.write:
        print("\ndry run -- nothing written. re-run with --write")
        return

    made = []
    for new, par, days, extra, pairs, base in rows:
        D = ff(base["kp"]) / (ff(base["eta"]) * ff(base["phi"]) * ff(base["beta"]))
        L = math.sqrt(4 * D * days * 86400)
        half = int(base["imax"]) * ff(base["ds"]) * 1000 / 2
        amb = ff(extra.get("a", base["a"])) - ff(extra.get("b", base["b"]))
        hdr = [f"! STAGE 11 -- built from res{par}.in.", "!"]
        if extra:
            hdr += [
                "! a-b SWEEP. Every pressure route to slip amplitude is closed:",
                "! Stage 9 put the amplitude crossing near kpmax 8e-13, where the",
                "! wellhead interpolates to +35-40% (outside the band), and Stage 10",
                "! showed tau_0 does nothing because the plateau tracks dp_crit down",
                "! with it. Friction is the only untouched knob -- a 0.015 / b 0.012",
                "! in all 75 scored runs.",
                "!",
                "! The fault never reaches failure; it creeps at",
                "! v = vref*exp[(tau/sigmabar - f0)/(a-b)], so a-b sets how far below",
                "! f0 it keeps creeping before it stalls. On 632901 that window is",
                "! d(tau/sigmabar) = 0.0278, i.e. 9.3 e-folds, and the 9.3 is set by",
                "! ln(velocity ratio) not by a-b -- so slip should be roughly LINEAR",
                "! in a-b, needing ~0.017 for the observed 2.81 cm at 5 d.",
                "!",
                f"! a {base['a']} -> {extra.get('a', base['a'])}, "
                f"b {base['b']} -> {extra.get('b', base['b'])}, so a-b "
                f"{ff(base['a'])-ff(base['b']):+.3f} -> {amb:+.3f}.",
                "! a-b stays POSITIVE: the fault remains velocity-strengthening, so",
                "! there is no nucleation-size constraint against ds = 5 m (h* exists",
                "! only for b > a).",
                "!",
                "! 632915 (a 0.015, b 0.005) and 632917 (a 0.022, b 0.012) share",
                "! a-b = 0.010 with different a. If they agree, only the difference",
                "! matters; if not, a acts on its own through the direct effect.",
                "!",
            ]
        else:
            hdr += [
                "! DURATION ONLY -- physics identical to the parent, so this extends",
                "! that run rather than replacing it.",
                "!",
                "! The observed slip file has 3,5,...,17 d but the parent stops at 5 d,",
                "! so the poster-style figure carries two curves against the",
                "! original's eight -- and both nearly coincide, the observed 3 d and",
                "! 5 d columns being byte-identical and simulated slip stalling after",
                "! the first flowing period (632901: 0.487 -> 0.488 cm).",
                "!",
                "! It also tests what the 5 d window cannot: lambda_obs falls from",
                "! 0.1866 over 0-5 d to 0.1412 over 0-8 d, so if the front stops",
                "! tracking sqrt(t) past 5 d the match is a 5 d artifact.",
                "!",
            ]
        if days != 18.0:
            hdr += [
                f"! {days:.0f} d, NOT 18. This is Wang & Dunham's config -- permev F at",
                f"! kp = {ff(base['kp']):.0e} with ds = {ff(base['ds'])*1000:.0f} m, so",
                f"! D = {D:.3f} m^2/s, three orders above the ds = 5 m runs, and",
                "! L/half reaches 0.80 at 14.9 d (18 d would be 0.88). Reaching 18 d",
                "! safely needs imax >= 700, 1.4x the cells, which changes the",
                "! discretisation and makes it a different run rather than a longer",
                f"! one. {days:.0f} d still covers observed times 3,5,7,9,11,13.",
                "!",
            ]
        hdr += [
            f"! Domain: D_far = {D:.3e} m^2/s, L({days:.0f} d) = {L:.0f} m against a",
            f"! {half:.0f} m half-domain, L/half = {L/half:.2f}.",
            f"! nstep {NSTEP} so a stall is the only failure mode, not exhaustion.",
        ]
        ov = {"tmax": f"{days/365.0:.8f}", "nstep": NSTEP, **extra}
        lines = [f"filenumber {new}"] + hdr
        for k, v in pairs:
            if k == "filenumber":
                continue
            lines.append(f"{k} {ov.get(k, v)}")
        (IN / f"res{new}.in").write_text("\n".join(lines) + "\n")
        made.append((new, par, days, extra))

    print("\nverifying each deck differs from its parent ONLY as intended")
    nbad = 0
    for new, par, days, extra in made:
        b = dict(read_deck(IN / f"res{par}.in"))
        d = dict(read_deck(IN / f"res{new}.in"))
        allowed = {"filenumber", "tmax", "nstep"} | set(extra)
        diff = {k for k in set(d) | set(b) if d.get(k) != b.get(k)}
        ck = {
            "only allowed keys changed": not (diff - allowed),
            "tmax correct": abs(ff(d["tmax"]) * 365 - days) < 1e-4,
            "nstep raised": int(d["nstep"]) >= int(b.get("nstep", 0)),
            "injection unchanged": d.get("injection_file") == b.get("injection_file"),
            "map unchanged": d.get("parameter_file") == b.get("parameter_file"),
            "kpmax unchanged": d.get("kpmax") == b.get("kpmax"),
            "muinit unchanged": d.get("muinit") == b.get("muinit"),
        }
        for k, v in extra.items():
            ck[f"{k} = {v}"] = d.get(k) == v
        ok = all(ck.values())
        nbad += not ok
        print(f"  res{new}.in ({days:.0f} d): "
              + ("OK" if ok else "PROBLEM " + str([k for k, x in ck.items() if not x]))
              + f"   changed {sorted(diff)}")
    print(f"\n  {len(made)-nbad}/{len(made)} clean")
    sys.exit(0 if nbad == 0 else 1)


if __name__ == "__main__":
    main()
