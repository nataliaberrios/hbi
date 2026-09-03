#!/usr/bin/env python3
"""Per-run figures, written into figures/<jobnum>/ alongside that run's params.txt.

Each job gets its own folder containing its own three figures plus params.txt, so a
single run can be looked at without untangling a multi-run overlay:

    figures/<job>/pressure_<job>.png   wellhead pressure vs the measurement
    figures/<job>/RT_<job>.png         slip front vs time, with the observed front
    figures/<job>/RV_<job>.png         slip front vs injected volume
    figures/<job>/params.txt           written by make_run_report.py

Front and observed-data machinery is copied verbatim from
cooper_basin_validation_stage_june15.ipynb (cells 10, 14, 74-77, 80, 81).

Two things this gets right that the overlay figures did not:
  - the front threshold is FIXED at FRONT_THR = 1e-4 m for every run. An earlier
    version of this file used each run's own dc and justified it on the grounds that
    dc sets both the rate-and-state length scale and the front definition. That
    double role is precisely the reason NOT to use it: both effects have the same
    sign, so a lower dc lowers the contour level (more cells counted as slipped)
    AND shrinks the nucleation size (slip propagates further), and a run scored at
    its own dc cannot separate the measurement artifact from the physics.
    Per-run dc is the right normalisation for "did this patch weaken at all" --
    that question is answered separately by peak_slip/dc in score_grid.py. It is
    the wrong normalisation for "how far did the front get relative to the data",
    which is a cross-run comparison against a single observed lambda and therefore
    needs one threshold for everybody. 64 of the 66 grid decks have dc = 1e-4, so
    this leaves their scores unchanged and makes the dc = 1.53e-5 Taiyi runs
    comparable rather than flattered.
  - the pressure metric is reported on FLOWING periods. During shut-ins the surface
    gauge reads a closed or bled well, which p_wh = p_downhole - rho g H + p_pipe
    cannot represent because it assumes flow. Including those windows inflated the
    ratio from 1.00 to 1.80.
  - kpmax/kL are shown ONLY for permev T runs. They are meaningless when
    permeability is fixed, and printing "kpmax ?" for a permev F deck is noise.

Usage:  python make_run_figures.py 632721 632720 ...
"""
import argparse
import glob
import os
import sys
from datetime import datetime as dt, timedelta as td
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

H = Path("/home/users/nberrios/3dhbi/hbi_analysis")
FIG = H / "figures"
IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
RHO, G, HW, DW, FD, P0, W = 1000.0, 9.81, 4077.0, 0.178, 0.015, 73.8, 6.0
FRONT_THR = 1e-4         # fixed slip threshold for the front, all runs -- see docstring

# validated with scripts/validate_palette.js --mode light: all six checks PASS
SIM, OBSC, MEAS = "#0072BD", "#a8071a", "#6b6b66"
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"

plt.rcParams.update({
    "font.size": 10.5, "axes.titlesize": 11.5, "axes.labelsize": 10.5,
    "axes.edgecolor": MUTED, "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
})


def style(ax):
    ax.grid(True, color=GRID, lw=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return ax


def ffloat(x):
    try:
        return float(str(x).replace("d", "e").replace("D", "e"))
    except (TypeError, ValueError):
        return None


def deck(n):
    d = {}
    for ln in open(IN / f"res{n}.in"):
        w = ln.split()
        if len(w) >= 2 and not ln.startswith("!"):
            d[w[0]] = w[1].strip('"')
    return d


def label_for(n, dk):
    """Only quote parameters that actually act in this run."""
    pb = (ffloat(dk.get("phi")) or 0) * (ffloat(dk.get("beta")) or 0)
    bits = [f"$\\mu_0$={ffloat(dk['muinit']):.2f}",
            f"$\\bar\\sigma_0$={ffloat(dk['sigmainit']):.2f} MPa",
            f"$\\phi\\beta$={pb:.1e}",
            f"dc={ffloat(dk['dc']):.2e}"]
    if dk.get("permev", "F").upper().startswith("T"):
        bits.append(f"permev T, kpmax={dk.get('kpmax', 'unset')}, kL={dk.get('kL', '?')}")
    else:
        bits.append("permev F (fixed perm)")
    return f"{n}:  " + ",  ".join(bits)


def observed():
    inj = loadmat(H / "Cooper_Basin_HAB_4_Injection_Rate.mat")["d"][0, 0]
    pr = loadmat(H / "Cooper_Basin_HAB_4_Wellhead_Pressure.mat")["d"][0, 0]
    ti = inj["Date"].squeeze(); ti = ti - ti[0]
    q = inj["Injection_rate"].squeeze() * (1000 / 60)
    m = np.isfinite(ti) & np.isfinite(q); ti, q = ti[m], q[m]
    o = np.argsort(ti); ti, q = ti[o], q[o]
    tp = pr["Date"].squeeze(); tp = tp - tp[0]
    pm = pr["Wellhead_pressure"].squeeze()
    m = np.isfinite(tp) & np.isfinite(pm); tp, pm = tp[m], pm[m]
    o = np.argsort(tp); tp, pm = tp[o], pm[o]

    mat = loadmat(H / "Cooper_Basin_Catalog_HAB_4.mat", squeeze_me=True,
                  struct_as_record=False)
    cat = {e.field: e.val for e in mat["Catalog"]}
    tt = cat["Time"].astype(float)
    dts = np.array([dt.fromordinal(int(x)) + td(days=x % 1) - td(days=366) for x in tt])
    dur = (dts >= dt(2012, 11, 13)) & (dts <= dt(2012, 11, 30))
    la, lo = -27.8115, 140.7596
    dd = np.sqrt(((cat["Lat"][dur] - la) * 111.0) ** 2
                 + ((cat["Long"][dur] - lo) * 111.0 * np.cos(np.radians(la))) ** 2)
    tdur = tt[dur] - tt[dur][0]
    org = np.median(dd[:10])
    o = np.argsort(tdur); st, sd = tdur[o], (dd - org)[o]
    ft, fd = [], []
    for i in range(0, len(st), 100):
        bt, bd = st[i:i + 100], sd[i:i + 100]
        if not len(bt):
            continue
        l_, h_ = np.percentile(bd, 90), np.percentile(bd, 95)
        k = (bd >= l_) & (bd <= h_)
        ft.extend(bt[k]); fd.extend(bd[k])
    return dict(ti=ti, q=q, tp=tp, pm=pm, ev_t=tdur, ev_d=dd, org=org,
                ft=np.array(ft), fd=np.array(fd),
                cumvol=cumtrapz(q, ti * 86400.0, initial=0.0))


def fit(x, R):
    if len(x) < 2:
        return np.nan
    b = np.sqrt(x)
    return float(np.sum(b * R) / np.sum(b ** 2))


def run_data(n, dk, thr=None):
    g = (glob.glob(f"/scratch/users/nberrios/3dhbi/runs/*/output/slip{n}.dat")
         + glob.glob(f"/scratch/users/nberrios/3dhbi/output/{n}/slip{n}.dat"))
    if not g:
        return None
    p = g[0]
    IM, JM = int(dk["imax"]), int(dk["jmax"])
    ds, NC = ffloat(dk["ds"]), IM * JM
    t = np.atleast_2d(np.loadtxt(p.replace("slip", "time")))[:, 1] / 86400.0
    nt = min(os.path.getsize(p) // (8 * NC), len(t))
    if nt < 5:
        return None
    thr = FRONT_THR if thr is None else thr
    s = np.memmap(p, np.float64, "r", shape=(nt, NC))
    out = {}
    for axis in ("downdip", "perpendicular"):
        if axis == "downdip":
            cs = np.array([np.asarray(s[k]).reshape(IM, JM)[:, int(JM / 2)]
                           for k in range(nt)]).T
            x = np.linspace(-ds * IM / 2, ds * IM / 2, IM); npt = IM
        else:
            cs = np.array([np.asarray(s[k]).reshape(IM, JM)[int(IM / 2), :]
                           for k in range(nt)]).T
            x = np.linspace(-ds * JM / 2, ds * JM / 2, JM); npt = JM
        R, T = [], []
        for i in range(npt):
            a = cs[i, :] > thr
            if a.any():
                R.append(abs(x[i])); T.append(t[np.argmax(a)])
        out[axis] = (np.array(R), np.array(T))
    pw = p.replace("slip", "pw")
    tp_, pp_ = None, None
    if os.path.exists(pw) and os.path.getsize(pw):
        a = np.atleast_2d(np.loadtxt(pw))
        f = (IN / dk["injection_file"]).read_text().split("\n")
        tq = np.array(f[2].split(), float); qq = np.array(f[4].split(), float) * W
        tp_ = a[:, 1] / 86400.0
        qi = np.interp(tp_ * 86400.0, tq, qq)
        pp_ = (P0 + a[:, 2]) - RHO * G * HW / 1e6 \
            + FD * (8.0 * HW * RHO * qi ** 2) / (np.pi ** 2 * DW ** 5) / 1e6
    return dict(fronts=out, t_end=t[nt - 1], thr=thr, tpw=tp_, ppw=pp_)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jobs", nargs="+", type=int)
    a = ap.parse_args()
    obs = observed()
    for n in a.jobs:
        if not (IN / f"res{n}.in").exists():
            print(f"  {n}: no deck"); continue
        dk = deck(n)
        d = run_data(n, dk)
        if d is None:
            print(f"  {n}: no usable output yet"); continue
        folder = FIG / str(n); folder.mkdir(parents=True, exist_ok=True)
        lab = label_for(n, dk)
        tcut = d["t_end"]

        # ---- pressure
        if d["tpw"] is not None:
            fig, ax = plt.subplots(figsize=(9.5, 4.8), constrained_layout=True)
            ax.plot(obs["tp"], obs["pm"], lw=1.2, color=MEAS, label="Measured wellhead")
            ax.plot(d["tpw"], d["ppw"], lw=1.8, color=SIM, label=f"{n} simulated $p_w$")
            hi = min(d["tpw"][-1], obs["tp"].max())
            gr = np.linspace(0.05, hi, 2000)
            ps = np.interp(gr, d["tpw"], d["ppw"])
            ob = np.interp(gr, obs["tp"], obs["pm"])
            qg = np.interp(gr, obs["ti"], obs["q"])
            fl = (qg > 0.25 * np.nanmax(obs["q"])) & (ob > 5.0)
            r = np.mean(ps[fl] / ob[fl]) if fl.sum() > 30 else np.nan
            rms = np.sqrt(np.mean((ps[fl] - ob[fl]) ** 2)) if fl.sum() > 30 else np.nan
            ax.set(xlabel="days since injection began",
                   ylabel="absolute wellhead pressure (MPa)", xlim=(0, hi * 1.02))
            # The measured wellhead record is the JUNE 2012 stage. A run driven by
            # a different injection file (the 1807/1808 family uses April) cannot
            # be scored against it -- say so on the figure instead of printing a
            # ratio that compares two different stimulations.
            inj = dk.get("injection_file", "?")
            if inj != "june_clean.txt":
                ax.set_title(f"{lab}\nNOT COMPARABLE: this run is driven by {inj}, "
                             f"the measured curve is the June 2012 stage",
                             fontsize=10, color=OBSC)
            else:
                ax.set_title(f"{lab}\nflowing-period ratio {r:.2f} "
                             f"({100*(r-1):+.1f}%), RMS {rms:.2f} MPa", fontsize=10)
            style(ax).legend(frameon=False, fontsize=9, labelcolor=INK)
            for e in ("png", "pdf"):
                fig.savefig(folder / f"pressure_{n}.{e}", dpi=200, bbox_inches="tight")
            plt.close(fig)

        # ---- R-T and R-V
        for key, xlab, fname in (("T", "days since injection began", "RT"),
                                 ("V", "cumulative injected volume (ML)", "RV")):
            fig, ax = plt.subplots(figsize=(9.5, 5.2), constrained_layout=True)
            ex = obs["ev_t"] if key == "T" else np.interp(
                obs["ev_t"], obs["ti"], obs["cumvol"] / 1e6)
            ax.scatter(ex, obs["ev_d"], s=4, alpha=0.20, color="gray",
                       edgecolors="none", label="observed events")
            fx = obs["ft"] if key == "T" else np.interp(
                obs["ft"], obs["ti"], obs["cumvol"] / 1e6)
            ax.scatter(fx, obs["fd"] + obs["org"], s=12, color=OBSC,
                       label="observed seismicity front")
            mo = fx <= (tcut if key == "T" else np.interp(tcut, obs["ti"], obs["cumvol"] / 1e6))
            lo_ = fit(fx[mo], obs["fd"][mo])
            o = np.argsort(fx[mo])
            ax.plot(fx[mo][o], (lo_ * np.sqrt(fx[mo]) + obs["org"])[o], color=OBSC, lw=2,
                    label=rf"observed fit $\lambda$={lo_:.3f}")
            R, T = d["fronts"]["downdip"]
            xx = T if key == "T" else np.interp(T, obs["ti"], obs["cumvol"] / 1e6)
            m = T <= tcut
            if m.sum() >= 2:
                ls_ = fit(xx[m], R[m])
                o = np.argsort(xx[m])
                ax.scatter(xx[m], R[m], s=12, color=SIM, alpha=0.8,
                           label="simulated slip front (downdip)")
                ax.plot(xx[m][o], (ls_ * np.sqrt(xx[m]))[o], color=SIM, lw=2,
                        label=rf"simulated fit $\lambda$={ls_:.3f}  "
                              rf"({ls_/lo_:.2f}$\times$ observed)")
            ax.set(xlabel=xlab, ylabel="distance from injector (km)")
            ax.set_title(f"{lab}\nfront threshold = {d['thr']:.2e} m (fixed, all runs); "
                         f"fits on 0–{tcut:.2f} d", fontsize=10)
            style(ax).legend(frameon=False, fontsize=8.5, labelcolor=INK, loc="upper left")
            for e in ("png", "pdf"):
                fig.savefig(folder / f"{fname}_{n}.{e}", dpi=200, bbox_inches="tight")
            plt.close(fig)
        print(f"  {n}: wrote {folder}/  (pressure_{n}, RT_{n}, RV_{n})")


if __name__ == "__main__":
    main()
