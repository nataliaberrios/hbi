#!/usr/bin/env python3
"""Overlay several HBI runs against the Cooper Basin data: wellhead pressure, R-T, R-V.

The notebook already produces R-T and R-V one job at a time. This does the thing
it cannot: put several runs on the SAME axes so a parameter change can be read
off directly, and fit every lambda over a COMMON window so the comparison is real.

All observed-data handling and all front/fit machinery is copied VERBATIM from
cooper_basin_validation_stage_june15.ipynb (cells 10, 14, 74, 75, 76, 77, 80, 81)
so the numbers are directly comparable to the per-job figures.

*** WHY THE COMMON WINDOW MATTERS ***
fit_sqrt_front least-squares fits over every point handed to it, so lambda depends
on how far the run got. Measured on 632510: the SAME run gives 0.1847 fit over
30.66 d and 0.1730 over 16.37 d -- 6.5% from the window alone, bigger than most
effects being chased. Runs of different duration therefore CANNOT have their
published lambdas compared. Everything here is refit on the shared span, and the
per-run full-span value is printed alongside so the size of the effect is visible.

Usage:
    python compare_runs_pressure_RT_RV.py 632522 632523 632524 632525
    python compare_runs_pressure_RT_RV.py --baseline 1808 632522 632524
"""
import argparse
import os
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from scipy.io import loadmat

try:
    from scipy.integrate import cumulative_trapezoid as cumtrapz
except ImportError:
    from scipy.integrate import cumtrapz

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path("/home/users/nberrios/3dhbi/hbi_analysis")
OUTDIR = HERE / "figures" / "cooper_basin_calibration"
SCRATCH = Path("/scratch/users/nberrios/3dhbi/output")
INPUT_DIR = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")

# ---- absolute-pressure constants, verbatim from notebook cell 14 --------------
RHO, G, HW, DW, FD = 1000.0, 9.81, 4077.0, 0.178, 0.015
P0_ABS_MPA = 73.8
FAULT_WIDTH_M = 6.0

DC = 1e-4

# validated with scripts/validate_palette.js --mode light (4 slots, all PASS,
# worst adjacent CVD dE 17.2 deutan). Assigned in fixed order, never cycled.
# Four hand-picked colourblind-safe hues for small comparisons. For more runs
# than that, extend with viridis rather than letting zip() truncate: `zip(runs,
# SERIES)` silently dropped 6 of a 10-run request, and the figure looked
# complete while missing most of what was asked for.
_SERIES4 = ["#D55E00", "#0072BD", "#009E73", "#8E44AD"]


def series_colors(n):
    """n distinct colours. Keeps the 4-colour palette exactly when n <= 4."""
    if n <= len(_SERIES4):
        return _SERIES4[:n]
    import matplotlib.pyplot as _plt
    return [_plt.cm.viridis(v) for v in __import__("numpy").linspace(0, 0.92, n)]


SERIES = _SERIES4
MEAS, OBSFRONT = "#6b6b66", "#a8071a"
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"

plt.rcParams.update({
    "font.size": 10.5, "axes.titlesize": 12, "axes.labelsize": 10.5,
    "axes.edgecolor": MUTED, "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
})


def style(ax):
    ax.grid(True, color=GRID, lw=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return ax


# ---- verbatim from notebook cell 75 ------------------------------------------
def calculate_seismicity_front_percentiles(event_times, event_distances, bin_size,
                                           lower_percentile, upper_percentile):
    sorted_indices = np.argsort(event_times)
    sorted_times = event_times[sorted_indices]
    sorted_distances = event_distances[sorted_indices]
    front_times, front_distances = [], []
    for start_idx in range(0, len(sorted_times), bin_size):
        end_idx = start_idx + bin_size
        bin_times = sorted_times[start_idx:end_idx]
        bin_distances = sorted_distances[start_idx:end_idx]
        if len(bin_times) == 0:
            continue
        lower_dist = np.percentile(bin_distances, lower_percentile)
        upper_dist = np.percentile(bin_distances, upper_percentile)
        mask = (bin_distances >= lower_dist) & (bin_distances <= upper_dist)
        front_times.extend(bin_times[mask])
        front_distances.extend(bin_distances[mask])
    return np.array(front_times), np.array(front_distances)


def dist_km(lat_vals, lon_vals, lat0, lon0):
    lat_diff = (lat_vals - lat0) * 111.0
    lon_diff = (lon_vals - lon0) * 111.0 * np.cos(np.radians(lat0))
    return np.sqrt(lat_diff**2 + lon_diff**2)


# ---- verbatim from notebook cell 76 ------------------------------------------
def extract_slip_cross_section(slip, nt, imax, jmax, axis, ds_km):
    if axis == 'downdip':
        cs_index = int(jmax / 2)
        slip_cross_section = np.zeros((imax, nt))
        for t in range(nt):
            grid = slip[:, t].reshape(imax, jmax)
            slip_cross_section[:, t] = grid[:, cs_index]
        x_coords = np.linspace(-ds_km * imax / 2, ds_km * imax / 2, imax)
        n_points = imax
    else:
        cs_index = int(imax / 2)
        slip_cross_section = np.zeros((jmax, nt))
        for t in range(nt):
            grid = slip[:, t].reshape(imax, jmax)
            slip_cross_section[:, t] = grid[cs_index, :]
        x_coords = np.linspace(-ds_km * jmax / 2, ds_km * jmax / 2, jmax)
        n_points = jmax
    return slip_cross_section, x_coords, n_points, cs_index


def calculate_slip_front(slip_cross_section, x_coords, time_years, Dc, n_points):
    cross_spatial, cross_temporal = [], []
    for i in range(n_points):
        above_dc = slip_cross_section[i, :] > Dc
        if np.any(above_dc):
            first_cross = np.argmax(above_dc)
            cross_spatial.append(abs(x_coords[i]))
            cross_temporal.append(time_years[first_cross])
    return np.array(cross_spatial), np.array(cross_temporal) * 365.0


# ---- verbatim from notebook cell 81 ------------------------------------------
def fit_sqrt_front(x, R):
    if len(x) < 2:
        return np.nan, None
    b = np.sqrt(x)
    lam = np.sum(b * R) / np.sum(b**2)
    return lam, lam * b


def read_deck(num):
    p = INPUT_DIR / f"res{num}.in"
    d = {}
    for line in p.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2 and not line.startswith("!"):
            d[parts[0]] = parts[1].strip('"')
    return d


def load_observed():
    """Notebook cells 10, 74, 77 -- everything that does not depend on a sim."""
    inj = loadmat(HERE / "Cooper_Basin_HAB_4_Injection_Rate.mat")["d"][0, 0]
    press = loadmat(HERE / "Cooper_Basin_HAB_4_Wellhead_Pressure.mat")["d"][0, 0]
    time_inj = inj["Date"].squeeze()
    rate = inj["Injection_rate"].squeeze() * (1000 / 60)      # m^3/min -> L/s
    time_inj_days = time_inj - time_inj[0]
    time_p = press["Date"].squeeze()
    pressure = press["Wellhead_pressure"].squeeze()
    time_p_days = time_p - time_p[0]

    mat = loadmat(HERE / "Cooper_Basin_Catalog_HAB_4.mat", squeeze_me=True,
                  struct_as_record=False)
    cat = {e.field: e.val for e in mat["Catalog"]}
    time = cat["Time"].astype(float)
    lat, lon, ml = cat["Lat"], cat["Long"], cat["ML"]
    dts = np.array([datetime.fromordinal(int(t)) + timedelta(days=t % 1)
                    - timedelta(days=366) for t in time])
    during = (dts >= datetime(2012, 11, 13)) & (dts <= datetime(2012, 11, 30))
    after = dts > datetime(2012, 11, 30)
    inj_lat, inj_lon = -27.8115, 140.7596

    d_during = dist_km(lat[during], lon[during], inj_lat, inj_lon)
    d_after = dist_km(lat[after], lon[after], inj_lat, inj_lon)
    t_during = time[during] - time[during][0]
    t_after = time[after] - time[during][0]

    cumvol = cumtrapz(rate, time_inj_days * 86400.0, initial=0.0)   # litres

    origin = np.median(d_during[:10])
    ft, fd = calculate_seismicity_front_percentiles(
        t_during, d_during - origin, bin_size=100,
        lower_percentile=90, upper_percentile=95)
    lam_obs, _ = fit_sqrt_front(ft, fd)
    print(f"observed seismicity-front lambda (full span): {lam_obs:.4f}")

    return dict(time_inj_days=time_inj_days, rate=rate, cumvol=cumvol,
                time_p_days=time_p_days, pressure=pressure,
                t_during=t_during, d_during=d_during, ml_during=ml[during],
                t_after=t_after, d_after=d_after, ml_after=ml[after],
                origin=origin, ft=ft, fd=fd, lam_obs=lam_obs)


def volume_ML(obs, t_days):
    return np.interp(t_days, obs["time_inj_days"], obs["cumvol"] / 1e6)


def pipe_friction_MPa(obs, t_days, deck):
    """Notebook cell 14: turbulent Darcy-Weisbach term from the HBI schedule."""
    inj = (INPUT_DIR / deck["injection_file"]).read_text().split("\n")
    tq = np.array(inj[2].split(), float)
    q = np.array(inj[4].split(), float) * FAULT_WIDTH_M          # m^2/s -> m^3/s
    qi = np.interp(t_days * 86400.0, tq, q)
    return FD * (8.0 * HW * RHO * qi**2) / (np.pi**2 * DW**5) / 1e6


def to_wellhead_abs(obs, t_days, dp_MPa, deck):
    """dp at depth -> absolute wellhead pressure, exactly as notebook cell 14."""
    return (P0_ABS_MPA + dp_MPa) - RHO * G * HW / 1e6 + pipe_friction_MPa(obs, t_days, deck)


def load_run(num, obs):
    base = SCRATCH / str(num)
    deck = read_deck(num)
    imax, jmax = int(deck["imax"]), int(deck["jmax"])
    ds_km = float(deck["ds"])
    ncell = imax * jmax
    r = dict(num=num, deck=deck, imax=imax, jmax=jmax, ds_km=ds_km)

    tfile = base / f"time{num}.dat"
    if not tfile.exists():
        raise FileNotFoundError(f"no output yet for {num} ({tfile})")
    t_days = np.loadtxt(tfile)[:, 1] / 86400.0

    # --- pressure: pw if the run wrote it, else pf at the injector cell
    inj = (INPUT_DIR / deck["injection_file"]).read_text().split("\n")
    i_inj, j_inj = (int(v) for v in inj[3].split())
    well_idx = (i_inj - 1) * jmax + (j_inj - 1)

    pwf = base / f"pw{num}.dat"
    if pwf.exists() and os.path.getsize(pwf) > 0:
        a = np.loadtxt(pwf)
        r["t_pw"], r["dp_pw"] = a[:, 1] / 86400.0, a[:, 2]
    else:
        r["t_pw"] = r["dp_pw"] = None

    pff = base / f"pf{num}.dat"
    nt_pf = min(os.path.getsize(pff) // (8 * ncell), len(t_days))
    pf = np.memmap(pff, np.float64, "r", shape=(nt_pf, ncell))
    r["t_pf"], r["dp_pf"] = t_days[:nt_pf], np.asarray(pf[:, well_idx])

    # --- slip fronts
    sf = base / f"slip{num}.dat"
    nt = min(os.path.getsize(sf) // (8 * ncell), len(t_days))
    slip = np.memmap(sf, np.float64, "r", shape=(nt, ncell)).T
    ty = t_days[:nt] / 365.0
    r["t_end"] = t_days[nt - 1]
    r["fronts"] = {}
    for axis in ("downdip", "perpendicular"):
        cs, x, n, _ = extract_slip_cross_section(slip, nt, imax, jmax, axis, ds_km)
        R, T = calculate_slip_front(cs, x, ty, DC, n)
        r["fronts"][axis] = dict(R=R, T=T, V=volume_ML(obs, T))
    print(f"  {num}: {nt} frames, 0-{r['t_end']:.3f} d, "
          f"pw {'yes' if r['t_pw'] is not None else 'NO (using pf at injector)'}, "
          + (f"permev T, kpmax {deck.get('kpmax','unset')}"
             if deck.get('permev','F').upper().startswith('T')
             else "permev F (fixed perm)"))
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jobs", nargs="+", type=int)
    ap.add_argument("--baseline", type=int, default=None,
                    help="older run to show as the pf-at-injector comparison")
    ap.add_argument("--tag", default="newruns")
    args = ap.parse_args()

    obs = load_observed()
    runs, skipped = [], []
    for n in args.jobs:
        try:
            runs.append(load_run(n, obs))
        except Exception as e:
            skipped.append((n, str(e)))
    base = None
    if args.baseline is not None:
        try:
            base = load_run(args.baseline, obs)
        except Exception as e:
            skipped.append((args.baseline, str(e)))
    if skipped:
        for n, e in skipped:
            print(f"  [skip] {n}: {e}")
    if not runs:
        sys.exit("no runs with output yet")

    # the shared span for every lambda
    t_cut = min(r["t_end"] for r in runs + ([base] if base else []))
    v_cut = float(volume_ML(obs, np.array([t_cut])))
    print(f"\ncommon fitting window: 0-{t_cut:.3f} d  (= {v_cut:.2f} ML)\n")

    def lam(r, axis, key):
        f = r["fronts"][axis]
        m = f["T"] <= t_cut
        l_cut, _ = fit_sqrt_front(f[key][m], f["R"][m])
        l_full, _ = fit_sqrt_front(f[key], f["R"])
        return l_cut, l_full

    print(f"{'job':>8s} {'permev':>7s} {'kpmax':>9s} "
          f"{'lamT_dd':>8s} {'lamT_pp':>8s} {'lamV_dd':>8s} {'lamV_pp':>8s} {'end(d)':>7s}")
    print("-" * 78)
    for r in runs + ([base] if base else []):
        ltd, _ = lam(r, "downdip", "T"); ltp, _ = lam(r, "perpendicular", "T")
        lvd, _ = lam(r, "downdip", "V"); lvp, _ = lam(r, "perpendicular", "V")
        print(f"{r['num']:>8d} {r['deck'].get('permev','?'):>7s} "
              f"{r['deck'].get('kpmax','-'):>9s} {ltd:8.4f} {ltp:8.4f} "
              f"{lvd:8.4f} {lvp:8.4f} {r['t_end']:7.2f}")
    print(f"{'observed':>8s} {'':>7s} {'':>9s} {obs['lam_obs']:8.4f}")

    OUTDIR.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------- 1. wellhead pressure
    fig, ax = plt.subplots(figsize=(10, 5.2), constrained_layout=True)
    ax.plot(obs["time_p_days"], obs["pressure"], lw=1.2, color=MEAS, alpha=0.85,
            label="Measured wellhead (Habanero 4)")
    if base is not None:
        p = to_wellhead_abs(obs, base["t_pf"], base["dp_pf"], base["deck"])
        ax.plot(base["t_pf"], p, lw=2.0, color=INK, ls=":",
                label=f"{base['num']} (old): $p_f$ at injector cell")
    for r, c in zip(runs, series_colors(len(runs))):
        if r["t_pw"] is not None:
            p = to_wellhead_abs(obs, r["t_pw"], r["dp_pw"], r["deck"])
            ax.plot(r["t_pw"], p, lw=1.8, color=c,
                    label=(f"{r['num']}: $p_w$  "
                           + (f"permev T, kpmax {r['deck'].get('kpmax','unset')}"
                              if r['deck'].get('permev','F').upper().startswith('T')
                              else "permev F (fixed perm)")))
    ax.axvline(t_cut, color=MUTED, lw=1, ls="--")
    style(ax).set(xlabel="days since injection began",
                  ylabel="absolute wellhead pressure (MPa)",
                  xlim=(0, max(t_cut * 1.05, 1)),
                  title="Wellhead pressure: new $p_w$ vs the old $p_f$-at-injector proxy")
    ax.legend(frameon=False, fontsize=9, labelcolor=INK)
    for ext in ("png", "pdf"):
        fig.savefig(OUTDIR / f"pressure_compare_{args.tag}.{ext}", dpi=300,
                    bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote pressure_compare_{args.tag}.png/.pdf")

    # ------------------------------------------------------------- 2/3. R-T, R-V
    for key, xlabel, fname, obsx in (
            ("T", "Time since injection start (days)", "RT", obs["t_during"]),
            ("V", "Cumulative injected volume (ML)", "RV", volume_ML(obs, obs["t_during"]))):
        fig, ax = plt.subplots(figsize=(9.5, 5.6), constrained_layout=True)
        big = obs["ml_during"] >= 3.0
        ax.scatter(obsx[~big], obs["d_during"][~big], s=4, alpha=0.22, color="gray",
                   edgecolors="none", label="Observed events during injection")
        ax.scatter(obsx[big], obs["d_during"][big], s=70, marker="*", color="#1f4e9c",
                   edgecolors="k", linewidths=0.5, label="Observed M$\\geq$3")
        fx = obs["ft"] if key == "T" else volume_ML(obs, obs["ft"])
        ax.scatter(fx, obs["fd"] + obs["origin"], s=13, color=OBSFRONT,
                   label="Observed seismicity front")
        o = np.argsort(fx)
        lo, _ = fit_sqrt_front(fx, obs["fd"])
        ax.plot(fx[o], (lo * np.sqrt(fx) + obs["origin"])[o], color=OBSFRONT, lw=2,
                label=rf"Observed fit: $R={lo:.3f}\sqrt{{{key}}}$")
        _rs = runs + ([base] if base else [])
        for r, c in zip(_rs, series_colors(len(runs)) + [INK]):
            f = r["fronts"]["downdip"]
            m = f["T"] <= t_cut
            l_cut, _ = fit_sqrt_front(f[key][m], f["R"][m])
            o = np.argsort(f[key][m])
            ax.plot(f[key][m][o], (l_cut * np.sqrt(f[key][m]))[o], lw=2, color=c,
                    ls=":" if (base is not None and r is base) else "-",
                    label=rf"{r['num']} downdip: $R={l_cut:.3f}\sqrt{{{key}}}$")
        style(ax).set(xlabel=xlabel, ylabel="Distance from injector (km)",
                      title=f"{fname}  —  every fit over the common window "
                            f"0–{t_cut:.1f} d ({v_cut:.1f} ML)")
        ax.legend(frameon=False, fontsize=8.5, labelcolor=INK, loc="upper left")
        for ext in ("png", "pdf"):
            fig.savefig(OUTDIR / f"{fname}_compare_{args.tag}.{ext}", dpi=300,
                        bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {fname}_compare_{args.tag}.png/.pdf")


if __name__ == "__main__":
    main()
