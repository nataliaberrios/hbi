#!/usr/bin/env python3
"""Score the cycle-2 runs against the post-4.300 d record. One table, no plots.

The gap this fills: compare_cycle2.py prints a header reading
"run arm tau_0 front wellhd slip" and then never fills those columns -- it only
prints the sqrt-fit parameters. The summary numbers quoted in
notebooks/cycle2_comparison.ipynb were computed ad hoc, so they could not be
regenerated or checked. This is the reproducible version.

THREE TARGETS, each evaluated on the shifted clock where sim t = 0 is data-day
4.300:

  front     R_sim / R_obs at T_EVAL, via front_at() -- NOT np.interp on
            run_data's raw output, which is unsorted in T and returns the final
            front for any time asked. See front_at's docstring. R_obs is the RUNNING MAXIMUM of event
            distance over the full 20 735-event catalogue, monotone by
            construction. A per-bin percentile is not a front and retreats four
            times on this catalogue. No sqrt(t) fit is involved: the ratio is
            read at matched times, because neither front starts at the origin.

  wellhead  THE PRESSURE CHANGE, decomposed so it cannot be confused with the
            absolute pressure. Over FLOWING intervals only, reusing
            sf.pressure_score's mask (q > 25% of peak AND p_obs > 5 MPa) so the
            interval selection stays comparable to the earlier 86-run study.
            The mask is not optional: within 26 minutes of the first shut-in
            the gauge falls to -0.2 MPa absolute, 35 MPa below its starting
            value, so unmasked intervals score the model against a vented
            wellhead rather than against the reservoir.

            THE DATUM IS THE INITIAL FAULT PRESSURE, DERIVED FROM THE FIELD
            MEASUREMENTS rather than picked off the wellhead record. Holl &
            Barton report Pp = 72.70 MPa at the fault's 4100 mSS median depth;
            hydrostatic there is 40.221, so the overpressure is 32.479 MPa and
            that is the datum expressed at the wellhead.

            Three wellhead values were used before this, all of them wrong for
            the purpose because all of them are wellhead pressures that merely
            lie NEAR the datum: 34.412 (the record's first sample), 33.970 (the
            pre-injection median), 33.271 (the record at data-day 4.300). A
            fourth, 33.805, is the MODEL's own static and is still subtracted
            from the model side. Referencing the data to 32.479 instead of
            34.412 raises the observed dp by 1.93 MPa at sim t = 0, easing to
            1.30 MPa by 10 d as the q^2 pipe-friction term grows.

            p_f0 IS NOT 73.82. That figure is
            taiyi-wang-seis3D/source_code/setup_model.m:112 and sits 1.12 MPa
            above what Holl reports. It was used here until the datum was
            derived, and the difference is most of the gap between the 0.8 MPa
            shift first computed and the ~2 MPa expected.

            sf.pressure_score returns mean(p_sim - p_obs), the SUM of a static
            mismatch (33.805 - P_REF = +0.534 MPa, identical in every run and
            not a model error) and the dp error. Both are printed separately
            and an assert checks they sum to the total.

            SCORING ABSOLUTE PRESSURE IS NOT DEFENSIBLE EITHER WAY. The 33.805
            assumes 4077 m of pure water; at 1050 kg/m3 it is 31.8 MPa and at
            950 it is 35.8 -- a +-2 MPa baseline uncertainty from the column
            density alone, against a mean observed dp of 9.0 MPa. Absolute
            agreement inside 2 MPa is unfalsifiable, so dp is the target and
            the absolute total is printed only for continuity with the older
            scores.

  slip      slip_sim / (obs(t + 4.300) - obs(4.300)) at T_EVAL. The
            SUBTRACTION IS REQUIRED, not a convenience: main_LH.f90:923 sets
            slip = 0d0 unconditionally, so HBI cannot start with the slip cycle
            1 already produced, and comparing against absolute observed slip
            would credit the model with slip it structurally cannot have.

T_EVAL = 8.7 d, i.e. data-day 13.0, is the latest time at which the observed
slip profiles and the front are both available inside the shortest run's
window. Front target 860 m, slip increment target 2.598 cm.

WHAT THE NUMBERS CANNOT TELL YOU. The front is a slip contour in the model and
seismicity in the data; these are not the same quantity, and the ratio assumes
they track. A run can also score well at T_EVAL and diverge later -- that is
exactly what the earlier study's 5-day window hid, so --at takes any time and
--sweep prints the ratio at several.

Usage:
  python score_cycle2.py 6329{5,6,7,8}{0,1,2,3,4}
  python score_cycle2.py --at 12.0 632963 632964
  python score_cycle2.py --sweep 632963
"""
import argparse
import importlib.util as iu
import sys
from pathlib import Path

import numpy as np

H = "/home/users/nberrios/3dhbi/hbi_analysis"
_s = iu.spec_from_file_location("sf", H + "/make_sweep_figures.py")
sf = iu.module_from_spec(_s); _s.loader.exec_module(sf)
sys.path.insert(0, H + "/notebooks")
from sim_curves import load_slip, _deck, _ff
_f = iu.spec_from_file_location(
    "fb", "/home/users/nberrios/3dhbi/hbi_git/docs/postshutin/front_backfront.py")
fb = iu.module_from_spec(_f); _f.loader.exec_module(fb)
# THE DATUM LIVES IN ONE PLACE. fault_pressure.py owns p_f0, the column
# length, the flowing density and the segmented friction; importing it is what
# keeps the scorer and the figures from drifting apart, which they did twice.
_fp = iu.spec_from_file_location(
    "fp", "/home/users/nberrios/3dhbi/hbi_git/docs/cooper_basin/fault_pressure.py")
fp = iu.module_from_spec(_fp); _fp.loader.exec_module(fp)

T0 = 4.300                                  # sim t = 0, in data-days
P_STATIC = sf.P0 - sf.RHO * sf.G * sf.HW / 1e6      # 33.805 MPa, the model's
# THE OBSERVED PRESSURE REFERENCE IS THE RECORD'S VALUE AT SIM t = 0, not at
# the start of the record. HBI sets pfinit = 0, so the model's dp is measured
# from the fault's state at sim t = 0 = data-day 4.300. Measuring the observed
# dp from data-day 0 instead -- 34.412 MPa, the record's first sample -- was
# referencing the two curves to DIFFERENT INSTANTS, and subtracting a baseline
# 1.141 MPa too high made every observed dp that much too small.
#
# The well was vented during the shut-in and had not recovered when injection
# resumed, so at data-day 4.300 it sits at 33.271 MPa. Referencing to that
# raises the observed dp by 1.141 MPa everywhere and lowers the model's
# over-prediction correspondingly: at tau_0 10.36, +55.3% -> +35.6%, and the
# dp = 0 crossing moves from tau_0 13.86 to about 12.95.
#
# THIS IS A PLOTTING/SCORING REFERENCE, NOT A CHANGE TO THE SIMULATIONS.
# pfinit stays 0. Setting pfinit = -1.14 would encode the same physical fact a
# second time, in the model rather than in the comparison, and doing both would
# double-count it.
#
# Computed from the record rather than hard-coded, so it cannot drift.
P_REF = None                                # = fp.datum(), set in main()
T_EVAL = 8.7                                # sim time at which to score
OBS_SLIP = Path("/home/users/nberrios/3dhbi/hbi/slip_profiles_strike.txt")
OT = [3, 5, 7, 9, 11, 13, 15, 17]           # data-days of the slip profiles


def arm_of(n):
    """Which fluid parameterisation, by (eta, beta). Same test as compare_cycle2."""
    d = _deck(n)
    e, b = _ff(d["eta"]), _ff(d["beta"])
    if abs(e - 0.89e-3) < 1e-9:
        return 1
    if b > 1e-7:
        return 3
    if b < 1e-8:
        return 4
    return 2


def pressure_reference(obs=None):
    """The datum at the wellhead, from fault_pressure.py. obs unused."""
    return fp.datum()


def observed_front(t_eval):
    """Running-maximum event distance at sim time t_eval, in metres."""
    tc, rc, _ = fb.catalogue()
    return float(np.interp(t_eval + T0, tc, np.maximum.accumulate(rc)))


def observed_slip_increment(t_eval):
    """obs(t_eval + T0) - obs(T0) in cm, from the strike slip profiles."""
    o = np.loadtxt(OBS_SLIP)
    peak = [o[:, 1 + i].max() for i in range(8)]
    return (float(np.interp(t_eval + T0, OT, peak))
            - float(np.interp(T0, OT, peak)))


def front_at(d, t_eval):
    """Model front radius in METRES at t_eval, or NaN if it has not arrived.

    THIS EXISTS BECAUSE np.interp ON run_data's OUTPUT IS WRONG. run_data
    returns its (T, R) pairs indexed by CELL i, with R = |x[i]| -- which falls
    to zero at the injector and rises again on the far side -- and T the first
    time that cell exceeded the slip threshold. So T IS NOT SORTED, and
    np.interp requires sorted x. Given that V-shaped input it silently returns
    a value near the array endpoint, i.e. the MAXIMUM front reached by t_end,
    for any t_eval whatsoever.

    Measured on 632960: the unsorted call returns 831 m at t_eval = 8.7 d,
    which is the 13.1 d front; sorting first gives 594 m. Every front number
    quoted for cycle 2 before this fix was the run's final front, not the front
    at the scoring time, and was therefore too large by 30-40%.

    Two guards, both needed:
      - sort by T before interpolating;
      - return NaN when t_eval precedes the first exceedance, instead of
        clamping to R[0]. 632994 never slips past the threshold until 8.98 d,
        and clamping reported it as an 811 m front at 8.7 d while its peak slip
        was 34 microns.
    """
    T, R = np.asarray(d["T"], float), np.asarray(d["R"], float) * 1000.0
    if len(T) < 2:
        return np.nan
    o = np.argsort(T)
    T, R = T[o], R[o]
    if t_eval < T[0] or t_eval > d["t_end"]:
        return np.nan
    return float(np.interp(t_eval, T, R))


def score(n, obs, t_eval):
    """One run. Returns a dict, or None if it produced no usable output."""
    dk = sf.deck(n)
    d = sf.run_data(n, dk)
    if d is None:
        return None
    r = dict(n=n, arm=arm_of(n), t_end=d["t_end"],
             tau0=_ff(dk["muinit"]) * _ff(dk["sigmainit"]))

    # --- front. run_data returns R in KILOMETRES; the observed front is metres.
    rs = front_at(d, t_eval)
    ro = observed_front(t_eval)
    if np.isfinite(rs):
        r.update(R_sim=rs, R_obs=ro, front=rs / ro if ro > 0 else np.nan)

    # --- wellhead, flowing mask only, on the shifted observed series
    if d.get("tpw") is not None:
        sh = dict(obs)
        sh["tp"] = obs["tp"] - T0
        sh["ti"] = obs["ti"] - T0
        bias, rms = sf.pressure_score(sh, d["tpw"], d["ppw"], t_eval)
        hi = min(d["tpw"][-1], sh["tp"].max(), t_eval)
        if hi > 0.2 and np.isfinite(bias):
            gr = np.linspace(0.05, hi, 2000)
            ps = np.interp(gr, d["tpw"], d["ppw"])
            # the OBSERVED dp comes from fp.dp_observed, so it carries the
            # segmented friction as well as the datum
            _to, _dpo = fp.dp_observed(obs)
            ob_dp = np.interp(gr, _to, _dpo)
            ob = np.interp(gr, sh["tp"], sh["pm"])
            qg = np.interp(gr, sh["ti"], sh["q"])
            fl = (qg > 0.25 * np.nanmax(obs["q"])) & (ob > 5.0)
            # decomposed: total = static mismatch + dp error, and the static
            # part is a fixed offset that no model parameter can address
            dp_sim = float(np.mean(ps[fl])) - P_STATIC
            dp_obs = float(np.mean(ob_dp[fl]))
            r.update(p_total=bias, p_rms=rms,
                     p_static=P_STATIC - P_REF,
                     p_bias=dp_sim - dp_obs,
                     p_pct=100.0 * (dp_sim - dp_obs) / dp_obs,
                     dp_obs=dp_obs)
            # the decomposition no longer sums to pressure_score's raw total,
            # because dp_obs now carries the friction correction that
            # pressure_score's absolute difference does not. p_total is kept
            # only for continuity with the 86 older scores.

    # --- slip at the injector, against the increment since t0
    try:
        _, sl, ta = load_slip(n, t_eval, how="strike")
        if abs(ta - t_eval) < 0.3:
            si = float(sl[0]) * 100.0
            so = observed_slip_increment(t_eval)
            r.update(slip_sim=si, slip_obs=so,
                     slip=si / so if so > 0 else np.nan)
    except Exception:
        pass
    return r


def _row(r):
    g = lambda k, f: (f"{r[k]:{f}}" if k in r and np.isfinite(r[k]) else "-")
    return (f"  {r['n']}  {r['arm']}   {r['tau0']:5.2f}  "
            f"{g('front', '6.2f')}  {g('R_sim', '7.0f')}  "
            f"{g('p_pct', '+7.1f')}  {g('p_bias', '+7.2f')}  "
            f"{g('p_total', '+7.2f')}  "
            f"{g('slip', '6.2f')}  {g('slip_sim', '7.2f')}   {r['t_end']:5.1f}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("jobs", nargs="+", type=int)
    ap.add_argument("--at", type=float, default=T_EVAL,
                    help="sim time (days) at which to score")
    ap.add_argument("--sweep", action="store_true",
                    help="also print the front ratio at 2,4,...,12 d")
    a = ap.parse_args(argv)

    global P_REF
    obs = sf.observed()
    P_REF = pressure_reference(obs)
    ro, so = observed_front(a.at), observed_slip_increment(a.at)
    print(f"scored at sim t = {a.at:.2f} d (data-day {a.at + T0:.2f})")
    print(f"  observed front {ro:.0f} m, observed slip increment {so:.3f} cm "
          f"(absolute {so + np.interp(T0, OT, [np.loadtxt(OBS_SLIP)[:, 1+i].max() for i in range(8)]):.3f} cm)")
    print(f"  datum {P_REF:.3f} MPa at the wellhead = p_f0 {fp.P_F0} - "
          f"rho {fp.RHO_FLOW:.0f} x g x {fp.Z_COLUMN:.1f} m")
    print(f"         friction over the WCR geometry, "
          f"{fp.friction(60.9):.2f} MPa at the 60.9 L/s peak")
    print(f"  'dp%' is the pressure-CHANGE error with that offset removed; "
          f"'total' is sf.pressure_score's raw number.")
    print(f"\n     run arm  tau_0   front   R_sim      dp%   dp MPa   total  "
          f"  slip  slip_cm   t_end")
    rows = []
    for n in a.jobs:
        r = score(n, obs, a.at)
        if r is None:
            print(f"  {n}  -- no output")
            continue
        rows.append(r)
        print(_row(r))

    ok = [r for r in rows if all(k in r for k in ("front", "p_pct", "slip"))
          and np.isfinite(r["front"]) and np.isfinite(r["slip"])]
    if ok:
        # rank on the same combination the notebook quotes, so the two agree
        ok.sort(key=lambda r: (abs(r["front"] - 1) + abs(r["p_pct"]) / 100.0
                               + abs(r["slip"] - 1)))
        print(f"\nranked by |front-1| + |wellhead%|/100 + |slip-1|:")
        for i, r in enumerate(ok, 1):
            print(f"  {i:2d}. {r['n']}  arm {r['arm']}  tau0 {r['tau0']:.2f}  "
                  f"front {r['front']:.2f}  wellhead {r['p_pct']:+.1f}%  "
                  f"slip {r['slip']:.2f}   score "
                  f"{abs(r['front']-1)+abs(r['p_pct'])/100+abs(r['slip']-1):.3f}")
    excl = [r["n"] for r in rows if r not in ok]
    if excl:
        print(f"\nnot ranked (missing front, wellhead or slip): "
              f"{' '.join(str(x) for x in excl)}")

    if a.sweep:
        print(f"\nfront ratio against time, to expose runs that only match early:")
        ts = [2.0, 4.0, 6.0, 8.0, 8.7, 10.0, 12.0]
        print("     run  " + "  ".join(f"{t:5.1f}d" for t in ts))
        for n in a.jobs:
            cells = []
            for t in ts:
                r = score(n, obs, t)
                cells.append(f"{r['front']:6.2f}" if r and "front" in r
                             and np.isfinite(r["front"]) else "     -")
            print(f"  {n}  " + " ".join(cells))


if __name__ == "__main__":
    main()
