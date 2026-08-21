#!/usr/bin/env python3
"""Live progress across all Cooper Basin runs: pressure match, front fit, figures.

Works on runs that are still going. The harness only rsyncs output into
/scratch/.../output/<num> when a run FINISHES, so until then everything lives in
that job's rundir. This discovers the mapping by reading the deck the harness
copied into each rundir, rather than needing a hand-maintained jobid table.

Produces, into hbi_analysis/figures/:
    comparisons/progress_pressure.png   every run's wellhead pressure vs measured
    comparisons/progress_front.png      lambda against the strength/storage params
    comparisons/progress_table.txt      the numbers
    <num>/params.txt                    refreshed per-run summary

Metrics reported, and why these ones:
  - mean(sim/measured) on pressure, NOT peak-vs-peak. Peak comparison flattered an
    earlier result badly: 632522 looks 1.9x on peaks and is 4.8x on the mean.
  - lambda refit on a window SHARED by every run being compared, because
    fit_sqrt_front fits whatever span it is given and the coefficient moves 6.5%
    between a 16 d and a 31 d fit of the same run.

Usage:  python progress_report.py            # everything it can find
        python progress_report.py 632538 632539
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
FIG = HERE / "figures"
IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
SCRATCH = Path("/scratch/users/nberrios/3dhbi/output")
OAK = Path("/oak/stanford/groups/edunham/nberrios/3doutput")
RUNS = Path("/scratch/users/nberrios/3dhbi/runs")

RHO, G, HW, DW, FD, P0, W = 1000.0, 9.81, 4077.0, 0.178, 0.015, 73.8, 6.0
DC = 1e-4
DAYS_PER_YEAR = 365.0
LAMBDA_OBS = 0.171          # published, full observed span
LAMBDA_OBS_WIN = 0.171      # refit per common window at runtime
TCUT_OVERRIDE = None        # --tcut, to fit on a longer span than the shortest run

INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
MEAS, OBS = "#6b6b66", "#a8071a"


def ffloat(x):
    try:
        return float(str(x).replace("d", "e").replace("D", "e"))
    except (ValueError, TypeError):
        return None


def read_deck(p):
    d = {}
    for line in Path(p).read_text().splitlines():
        if line.startswith("!"):
            continue
        w = line.split()
        if len(w) >= 2:
            d[w[0]] = w[1].strip('"')
    return d


def discover():
    """filenumber -> (data dir, deck path). Rundirs win only if not yet finished."""
    found = {}
    for base in (SCRATCH, OAK):
        if not base.exists():
            continue
        for d in base.iterdir():
            if d.is_dir() and d.name.isdigit() and (d / f"time{d.name}.dat").exists():
                found[int(d.name)] = (d, IN / f"res{d.name}.in")
    for rd in RUNS.glob("*/res*.in"):
        try:
            num = int(read_deck(rd).get("filenumber", -1))
        except ValueError:
            continue
        t = rd.parent / "output" / f"time{num}.dat"
        if num > 0 and num not in found and t.exists() and os.path.getsize(t):
            found[num] = (rd.parent / "output", rd)
    return dict(sorted(found.items()))


def observed_front(obs):
    """The observed seismicity front, verbatim from notebook cells 74/75/77."""
    d, t = obs["d_during"], obs["t_during"]
    origin = np.median(d[:10])
    ev_t, ev_d = t, d - origin
    o = np.argsort(ev_t); st, sd = ev_t[o], ev_d[o]
    ft, fdi = [], []
    for i in range(0, len(st), 100):
        bt, bd = st[i:i+100], sd[i:i+100]
        if not len(bt):
            continue
        lo, hi = np.percentile(bd, 90), np.percentile(bd, 95)
        k = (bd >= lo) & (bd <= hi)
        ft.extend(bt[k]); fdi.extend(bd[k])
    return np.array(ft), np.array(fdi), origin


def observed():
    inj = loadmat(HERE / "Cooper_Basin_HAB_4_Injection_Rate.mat")["d"][0, 0]
    pr = loadmat(HERE / "Cooper_Basin_HAB_4_Wellhead_Pressure.mat")["d"][0, 0]
    ti = inj["Date"].squeeze(); rate = inj["Injection_rate"].squeeze() * (1000 / 60)
    ti = ti - ti[0]
    tp = pr["Date"].squeeze(); tp = tp - tp[0]; pm = pr["Wellhead_pressure"].squeeze()
    g = np.isfinite(tp) & np.isfinite(pm) & (pm > 0)
    tp, pm = tp[g], pm[g]
    o = np.argsort(tp); tp, pm = tp[o], pm[o]
    from datetime import datetime as _dt, timedelta as _td
    mat = loadmat(HERE / "Cooper_Basin_Catalog_HAB_4.mat", squeeze_me=True,
                  struct_as_record=False)
    cat = {e.field: e.val for e in mat["Catalog"]}
    tt = cat["Time"].astype(float)
    dts = np.array([_dt.fromordinal(int(x)) + _td(days=x % 1) - _td(days=366) for x in tt])
    dur = (dts >= _dt(2012, 11, 13)) & (dts <= _dt(2012, 11, 30))
    la, lo0 = -27.8115, 140.7596
    dd = np.sqrt(((cat["Lat"][dur]-la)*111.0)**2
                 + ((cat["Long"][dur]-lo0)*111.0*np.cos(np.radians(la)))**2)
    return dict(ti=ti, rate=rate, tp=tp, pm=pm,
                cumvol=cumtrapz(rate, ti * 86400.0, initial=0.0),
                t_during=tt[dur]-tt[dur][0], d_during=dd)


def inj_rate(deck, t_days):
    inj = (IN / deck["injection_file"]).read_text().split("\n")
    tq = np.array(inj[2].split(), float)
    qq = np.array(inj[4].split(), float) * W
    return np.interp(t_days * 86400.0, tq, qq)


def wellhead(deck, t_days, dp):
    q = inj_rate(deck, t_days)
    pipe = FD * (8.0 * HW * RHO * q**2) / (np.pi**2 * DW**5) / 1e6
    return (P0 + dp) - RHO * G * HW / 1e6 + pipe


def pressure_scores(deck, t, p, obs):
    """Wellhead agreement in MPa, split by whether the well is flowing.

    NOT mean(sim/measured). The measured wellhead falls to ~0 during shut-ins, so
    that ratio divides by a near-zero number: the same run reads 2.17 on one
    sampling and 0.32 on another, and goes negative if the pm>0 filter is dropped.
    Bias and RMS in MPa are division-free and mean the same thing everywhere.

    The split matters physically. HBI has no wellbore storage or bleed-off path,
    so when injection stops the model cannot vent the well -- simulated wellhead
    stays near its flowing value while the real one drops ~20 MPa. That is a
    structural limitation of the model, identical for every parameter set, so it
    cannot discriminate between them; only the flowing periods can.
    """
    hi = min(t[-1], obs["tp"].max())
    if hi <= 0.2:
        return {}
    grid = np.linspace(0.05, hi, 2000)          # uniform, sampling-independent
    ps = np.interp(grid, t, p)
    om = np.interp(grid, obs["tp"], obs["pm"])
    q = inj_rate(deck, grid)
    flow = q > 1e-6
    d = ps - om
    out = dict(span=hi, frac_flow=float(flow.mean()),
               rms_all=float(np.sqrt(np.mean(d**2))), bias_all=float(d.mean()))
    if flow.sum() >= 5:
        df = d[flow]
        out.update(rms_flow=float(np.sqrt(np.mean(df**2))), bias_flow=float(df.mean()),
                   sim_flow=float(ps[flow].mean()), meas_flow=float(om[flow].mean()))
    if (~flow).sum() >= 5:
        out.update(bias_shut=float(d[~flow].mean()))
    return out


def slip_front(base, num, deck, axis="downdip"):
    imax, jmax = int(deck["imax"]), int(deck["jmax"])
    ds_km = ffloat(deck["ds"]); nc = imax * jmax
    sf = base / f"slip{num}.dat"
    t = np.loadtxt(base / f"time{num}.dat")
    t = np.atleast_2d(t)[:, 1] / 86400.0
    nt = min(os.path.getsize(sf) // (8 * nc), len(t))
    if nt < 5:
        return None
    slip = np.memmap(sf, np.float64, "r", shape=(nt, nc)).T
    if axis == "downdip":
        idx, n, coords = int(jmax / 2), imax, np.linspace(-ds_km*imax/2, ds_km*imax/2, imax)
        cs = np.array([slip[:, k].reshape(imax, jmax)[:, idx] for k in range(nt)]).T
    else:
        idx, n, coords = int(imax / 2), jmax, np.linspace(-ds_km*jmax/2, ds_km*jmax/2, jmax)
        cs = np.array([slip[:, k].reshape(imax, jmax)[idx, :] for k in range(nt)]).T
    R, T = [], []
    for i in range(n):
        a = cs[i, :] > DC
        if a.any():
            R.append(abs(coords[i])); T.append(t[np.argmax(a)])
    return np.array(R), np.array(T), t[nt - 1]


def fit_sqrt(x, R):
    if len(x) < 2:
        return np.nan
    b = np.sqrt(x)
    return float(np.sum(b * R) / np.sum(b**2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jobs", nargs="*", type=int)
    ap.add_argument("--tcut", type=float, default=None,
                    help="front-fit window in days; default = shortest run's reach")
    a = ap.parse_args()
    if a.tcut:
        globals()["TCUT_OVERRIDE"] = a.tcut
    obs = observed()
    found = discover()
    if a.jobs:
        found = {k: v for k, v in found.items() if k in a.jobs}
    if not found:
        sys.exit("nothing found")

    rows = []
    for num, (base, deckp) in found.items():
        if not Path(deckp).exists():
            continue
        deck = read_deck(deckp)
        try:
            fr = slip_front(base, num, deck)
        except Exception as e:
            print(f"  {num}: front failed ({e})"); fr = None
        pw = base / f"pw{num}.dat"
        r = dict(num=num, deck=deck, base=base)
        if pw.exists() and os.path.getsize(pw):
            arr = np.atleast_2d(np.loadtxt(pw))
            r["t"] = arr[:, 1] / 86400.0
            r["p"] = wellhead(deck, r["t"], arr[:, 2])
        if fr:
            r["R"], r["T"], r["t_end"] = fr
        mu, sig, f0 = (ffloat(deck.get(k)) for k in ("muinit", "sigmainit", "f0"))
        phi, beta = ffloat(deck.get("phi")), ffloat(deck.get("beta"))
        r["dp_crit"] = sig - mu * sig / f0 if None not in (mu, sig, f0) else None
        r["phibeta"] = phi * beta if None not in (phi, beta) else None
        r["mu"] = mu
        rows.append(r)

    have_front = [r for r in rows if "T" in r]
    t_cut = min(r["t_end"] for r in have_front) if have_front else 0.0
    # A run launched minutes ago would otherwise drag the shared window down to
    # its own reach and make every lambda in the table a fit over a fraction of a
    # day -- 632721 reads 1.15 on 7.8 d and 0.37 on 0.87 d, purely from the window.
    if TCUT_OVERRIDE:
        t_cut = TCUT_OVERRIDE
    limiter = [r["num"] for r in have_front if abs(r["t_end"] - t_cut) < 1e-6]
    longest = max((r["t_end"] for r in have_front), default=0.0)

    # observed seismicity front, refit on the SAME window as the sims
    obs_ft, obs_fd, obs_origin = observed_front(obs)
    m_obs = obs_ft <= t_cut
    lam_obs_win = fit_sqrt(obs_ft[m_obs], obs_fd[m_obs]) if m_obs.sum() > 2 else LAMBDA_OBS
    globals()["LAMBDA_OBS_WIN"] = lam_obs_win

    FIG.mkdir(parents=True, exist_ok=True)
    (FIG / "comparisons").mkdir(exist_ok=True)
    L = [f"progress as of {datetime.now():%Y-%m-%d %H:%M}",
         f"common front-fit window: 0-{t_cut:.2f} d "
         f"(lambda is window-dependent, so every fit uses this span)"]
    if longest > 1.5 * t_cut and limiter:
        L += [f"  WARNING: the window is set by run(s) {limiter} at {t_cut:.2f} d, "
              f"while the longest run reaches {longest:.2f} d.",
              f"  Every lambda below is therefore a short-window fit. Pass "
              f"--tcut to fit on a longer span and drop the short run(s)."]
    L += ["",
          f"{'run':>8s} {'muinit':>7s} {'dp_crit':>8s} {'phi*beta':>9s} {'permev':>7s} "
          f"{'kpmax':>9s} {'reached':>8s} {'lam_win':>8s} {'lam_own':>8s} "
          f"{'pRMS_fl':>8s} {'pbias_fl':>9s} {'pbias_shut':>11s}",
          "-" * 120]
    for r in sorted(rows, key=lambda x: x["num"]):
        lam = fit_sqrt(r["T"][r["T"] <= t_cut], r["R"][r["T"] <= t_cut]) if "T" in r else np.nan
        # lam on the run's OWN full span too, so the shared-window choice hides nothing
        lam_own = fit_sqrt(r["T"], r["R"]) if "T" in r else np.nan
        sc = pressure_scores(r["deck"], r["t"], r["p"], obs) if "t" in r else {}
        r["scores"] = sc
        L.append(f"{r['num']:>8d} {r['mu'] if r['mu'] else float('nan'):>7.2f} "
                 f"{r['dp_crit'] if r['dp_crit'] else float('nan'):>8.2f} "
                 f"{r['phibeta'] if r['phibeta'] else float('nan'):>9.2e} "
                 f"{r['deck'].get('permev','?'):>7s} {r['deck'].get('kpmax','-'):>9s} "
                 f"{r.get('t_end', float('nan')):>8.2f} {lam:>8.4f} {lam_own:>8.4f} "
                 f"{sc.get('rms_flow', float('nan')):>8.2f} "
                 f"{sc.get('bias_flow', float('nan')):>+9.2f} "
                 f"{sc.get('bias_shut', float('nan')):>+11.2f}")
        r["lam"] = lam
    L += ["",
          f"observed lambda on this window = {LAMBDA_OBS_WIN:.4f}",
          "pressure columns are MPa, sim minus measured, on a uniform time grid:",
          "  pRMS_fl / pbias_fl  while the well is FLOWING -- the only periods that",
          "                      discriminate between parameter sets",
          "  pbias_shut          during SHUT-INS. Large and positive for every run:",
          "                      HBI has no wellbore bleed-off, so it cannot follow a",
          "                      shut-in. Structural, not a parameter choice.",
          "mean(sim/measured) is deliberately not reported -- measured wellhead falls",
          "to ~0 during shut-ins, so the ratio is a division by near-zero."]
    txt = "\n".join(L)
    (FIG / "comparisons" / "progress_table.txt").write_text(txt + "\n")
    print(txt)

    # ---- pressure figure
    fig, ax = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
    ax.plot(obs["tp"], obs["pm"], lw=1.3, color=MEAS, label="Measured wellhead")
    cm = plt.get_cmap("viridis")
    withp = [r for r in rows if "t" in r]
    for r, c in zip(withp, [cm(x) for x in np.linspace(0.08, 0.92, max(len(withp), 1))]):
        ax.plot(r["t"], r["p"], lw=1.4, color=c,
                label=f"{r['num']} mu{r['mu']:.2f} pb{r['phibeta']:.0e}")
    ax.set(xlabel="days since injection began", ylabel="absolute wellhead pressure (MPa)",
           title="Wellhead pressure, all runs (some still in progress)")
    ax.grid(True, color=GRID, lw=0.6, alpha=0.7); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=7.5, ncol=2, labelcolor=INK)
    for e in ("png", "pdf"):
        fig.savefig(FIG / "comparisons" / f"progress_pressure.{e}", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ---- lambda vs the two parameters that the scaling says control it
    good = [r for r in rows if np.isfinite(r.get("lam", np.nan)) and r["dp_crit"] and r["phibeta"]]
    if good:
        fig, ax = plt.subplots(figsize=(8, 5.2), constrained_layout=True)
        x = np.array([r["phibeta"] * r["dp_crit"] for r in good])
        y = np.array([r["lam"] for r in good])
        ax.axhline(LAMBDA_OBS_WIN, color=OBS, lw=2, ls="--")
        ax.annotate(f"observed {LAMBDA_OBS_WIN:.3f} (this window)", xy=(0.98, LAMBDA_OBS_WIN),
                    xycoords=("axes fraction", "data"), xytext=(0, 5),
                    textcoords="offset points", ha="right", fontsize=9, color=OBS)
        sc = ax.scatter(x, y, c=[r["mu"] for r in good], s=90, cmap="plasma",
                        edgecolors=INK, linewidths=0.6, zorder=3)
        for r, xx, yy in zip(good, x, y):
            ax.annotate(str(r["num"]), xy=(xx, yy), xytext=(6, 4),
                        textcoords="offset points", fontsize=7.5, color=INK)
        # volume balance: lambda ~ (phi beta dp_crit)^-1/2
        xs = np.logspace(np.log10(x.min()*0.7), np.log10(x.max()*1.4), 50)
        ref = y[np.argmin(np.abs(x - np.median(x)))] * np.sqrt(np.median(x) / xs)
        ax.plot(xs, ref, color=MUTED, lw=1, ls=":",
                label=r"volume balance: $\lambda\propto(\phi\beta\,\Delta p_{crit})^{-1/2}$")
        ax.set_xscale("log")
        ax.set(xlabel=r"$\phi\beta\,\Delta p_{crit}$   (Pa$^{-1}$ x MPa)",
               ylabel=r"front coefficient $\lambda$")
        ax.set_title(f"Front vs storage x strength margin, fits over 0-{t_cut:.1f} d")
        plt.colorbar(sc, ax=ax, label="muinit")
        ax.grid(True, color=GRID, lw=0.6, alpha=0.7); ax.set_axisbelow(True)
        for s in ("top", "right"): ax.spines[s].set_visible(False)
        ax.legend(frameon=False, fontsize=8.5, labelcolor=INK)
        for e in ("png", "pdf"):
            fig.savefig(FIG / "comparisons" / f"progress_front.{e}", dpi=200, bbox_inches="tight")
        plt.close(fig)
    # ---- THE JOINT CONSTRAINT: both targets on one axes -----------------------
    # Matching pressure alone is not the goal. A run is only a candidate if it sits
    # near (1,1): pressure ratio 1 AND lambda/lambda_obs 1. Plotting them separately
    # hides runs that win one target by sacrificing the other.
    joint = []
    for r in rows:
        if not np.isfinite(r.get("lam", np.nan)) or "t" not in r:
            continue
        hi = min(r["t"][-1], obs["tp"].max())
        if hi <= 0.2:
            continue
        sc = r.get("scores") or pressure_scores(r["deck"], r["t"], r["p"], obs)
        if "sim_flow" not in sc:
            continue
        # ratio of MEANS on flowing periods. Mean-of-ratios over the whole span
        # divides by a measured wellhead that falls to ~0 in shut-ins.
        joint.append((r, sc["sim_flow"] / sc["meas_flow"], r["lam"] / LAMBDA_OBS_WIN))
    if joint:
        fig, ax = plt.subplots(figsize=(8.4, 6.4), constrained_layout=True)
        ax.axhline(1.0, color=OBS, lw=1.2, ls="--")
        ax.axvline(1.0, color=OBS, lw=1.2, ls="--")
        ax.plot([1], [1], marker="*", ms=22, color=OBS, zorder=5)
        ax.annotate("TARGET\nboth matched", xy=(1, 1), xytext=(10, -26),
                    textcoords="offset points", fontsize=9.5, color=OBS, weight="bold")
        for r, pr, lr in joint:
            # agreed tolerance: +/-15% on BOTH targets (Taiyi does not match
            # the wellhead data exactly either); +/-10% flagged as strong
            ok = abs(pr - 1) <= 0.15 and abs(lr - 1) <= 0.15
            strong = abs(pr - 1) <= 0.10 and abs(lr - 1) <= 0.10
            ax.scatter([pr], [lr], s=110 if ok else 70,
                       color=("#009E73" if strong else "#B8860B") if ok else "#0072BD",
                       edgecolors=INK, linewidths=0.7, zorder=4)
            ax.annotate(f"{r['num']}\n$\\mu$={r['mu']:.2f}", xy=(pr, lr),
                        xytext=(7, 5), textcoords="offset points", fontsize=7.5, color=INK)
        ax.set(xlabel="pressure while flowing:  mean(sim) / mean(measured)      1.0 = matched",
               ylabel=r"front:  $\lambda\,/\,\lambda_{obs}$      1.0 = matched",
               title="Both constraints at once -- a candidate must be near the star")
        ax.grid(True, color=GRID, lw=0.6, alpha=0.7); ax.set_axisbelow(True)
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
        for e in ("png", "pdf"):
            fig.savefig(FIG / "comparisons" / f"progress_joint.{e}", dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"\njoint-constraint candidates (both within tolerance):")
        for r, pr, lr in sorted(joint, key=lambda z: abs(z[1]-1)+abs(z[2]-1)):
            ok = abs(pr-1) <= 0.15 and abs(lr-1) <= 0.15
            flag = ("  <== STRONG (both within 10%)" if abs(pr-1) <= 0.10 and abs(lr-1) <= 0.10
                    else "  <== PASS (both within 15%)" if ok else "")
            print(f"  {r['num']}  mu {r['mu']:.2f}  p_flow {pr:6.2f}  lam/lam_obs {lr:5.2f}{flag}")
    print(f"\nwrote {FIG}/comparisons/progress_*")


if __name__ == "__main__":
    main()
