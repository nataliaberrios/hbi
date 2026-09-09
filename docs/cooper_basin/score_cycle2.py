#!/usr/bin/env python3
"""Score the cycle-2 runs against the post-4.300 d record. One table, no plots.

The gap this fills: compare_cycle2.py prints a header reading
"run arm tau_0 front wellhd slip" and then never fills those columns -- it only
prints the sqrt-fit parameters. The summary numbers quoted in
notebooks/cycle2_comparison.ipynb were computed ad hoc, so they could not be
regenerated or checked. This is the reproducible version.

THREE TARGETS, each evaluated on the shifted clock where sim t = 0 is data-day
4.300:

  front     R_sim / R_obs at T_EVAL. R_obs is the RUNNING MAXIMUM of event
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

            THREE BASELINES EXIST AND MIXING THEM IS THE TRAP:

                33.805  the MODEL's static wellhead, P0 - rho g H with
                        P0 = 73.8 MPa, rho = 1000, H = 4077 m
                        (make_sweep_figures.py:55, applied in run_data)
                34.412  the observed record's first sample -- Wang & Dunham's
                        reference, and what the slide figure uses
                33.970  the observed pre-injection median over t < 0.501 d

            sf.pressure_score returns mean(p_sim - p_obs), which is the SUM of
            a static mismatch (33.805 - P_REF = -0.607 MPa, identical in every
            run, not a model error) and the dp error. Quoting that total and
            then dividing it by dp measured from 34.412 mixes the model's
            baseline with the observed one and UNDERSTATES the dp error by
            0.607 MPa -- about 8 percentage points, enough to move the
            best-fitting tau_0 from 13.86 to 13.3. Both parts are printed
            separately, and an assert checks they sum to the total.

            SCORING ABSOLUTE PRESSURE IS NOT DEFENSIBLE ANYWAY. The 33.805
            assumes 4077 m of pure water. At 1050 kg/m3 it is 31.8 MPa, at 950
            it is 35.8 -- a +-2 MPa baseline uncertainty from the column
            density alone, against a mean observed dp of 7.9 MPa. Absolute
            agreement inside 2 MPa is therefore unfalsifiable. dp is the
            target; the absolute total is printed only for continuity with the
            older scores.

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

T0 = 4.300                                  # sim t = 0, in data-days
P_STATIC = sf.P0 - sf.RHO * sf.G * sf.HW / 1e6      # 33.805 MPa, the model's
P_REF = 34.412                              # observed first sample; see docstring
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


def score(n, obs, t_eval):
    """One run. Returns a dict, or None if it produced no usable output."""
    dk = sf.deck(n)
    d = sf.run_data(n, dk)
    if d is None:
        return None
    r = dict(n=n, arm=arm_of(n), t_end=d["t_end"],
             tau0=_ff(dk["muinit"]) * _ff(dk["sigmainit"]))

    # --- front. run_data returns R in KILOMETRES; the observed front is metres.
    if d["t_end"] >= t_eval and len(d["T"]) > 2:
        rs = float(np.interp(t_eval, d["T"], d["R"])) * 1000.0
        ro = observed_front(t_eval)
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
            ob = np.interp(gr, sh["tp"], sh["pm"])
            qg = np.interp(gr, sh["ti"], sh["q"])
            fl = (qg > 0.25 * np.nanmax(obs["q"])) & (ob > 5.0)
            # decomposed: total = static mismatch + dp error, and the static
            # part is a fixed offset that no model parameter can address
            dp_sim = float(np.mean(ps[fl])) - P_STATIC
            dp_obs = float(np.mean(ob[fl])) - P_REF
            r.update(p_total=bias, p_rms=rms,
                     p_static=P_STATIC - P_REF,
                     p_bias=dp_sim - dp_obs,
                     p_pct=100.0 * (dp_sim - dp_obs) / dp_obs,
                     dp_obs=dp_obs)
            assert abs(bias - (r["p_static"] + r["p_bias"])) < 1e-6, \
                "decomposition does not sum to pressure_score's total"

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

    obs = sf.observed()
    ro, so = observed_front(a.at), observed_slip_increment(a.at)
    print(f"scored at sim t = {a.at:.2f} d (data-day {a.at + T0:.2f})")
    print(f"  observed front {ro:.0f} m, observed slip increment {so:.3f} cm "
          f"(absolute {so + np.interp(T0, OT, [np.loadtxt(OBS_SLIP)[:, 1+i].max() for i in range(8)]):.3f} cm)")
    print(f"  wellhead: model static {P_STATIC:.3f} MPa vs observed reference "
          f"{P_REF:.3f}, a fixed {P_STATIC - P_REF:+.3f} MPa offset in EVERY run.")
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
