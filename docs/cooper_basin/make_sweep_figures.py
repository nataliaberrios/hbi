#!/usr/bin/env python3
"""Meeting figures: one panel row per parameter sweep, showing how the fit responds.

Each sweep varies ONE parameter with everything else held fixed, and gets one
figure with the three panels that matter:

    pressure  |  R-T (front vs time)  |  R-V (front vs injected volume)

so the effect of that parameter on both targets is read off a single row. Nothing
new is plotted here -- these are the same three plots as make_run_figures.py, just
overlaid across a sweep and organised so the argument is visible.

The front and observed-data machinery is copied verbatim from
cooper_basin_validation_stage_june15.ipynb (cells 10, 14, 74-77, 80, 81).

Two conventions that keep the comparison honest:

  - lambda is fit on ONE window shared by every run in the sweep, and the observed
    front is refit on that same window. lambda depends on the window: the same run
    gives 0.1847 over 30.7 d and 0.1730 over 16.4 d, so comparing a run fit over
    8 d against one fit over 17 d compares windows, not physics. Each figure states
    its window, and any run too short to reach it is dropped and named.

  - pressure is scored on FLOWING periods, in MPa. HBI has no wellbore bleed-off,
    so during a shut-in the simulated wellhead stays near its flowing value while
    the measurement drops ~20 MPa. That gap is identical for every parameter set,
    so it cannot discriminate between them -- and mean(sim/measured) over the whole
    span divides by a measured pressure that falls to ~0, which made the same run
    read 2.17 on one sampling and 0.32 on another.

Usage:  python make_sweep_figures.py              # all sweeps
        python make_sweep_figures.py muinit       # one sweep by name
"""
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
OUT = H / "figures" / "sweeps"
IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
RHO, G, HW, DW, FD, P0, W = 1000.0, 9.81, 4077.0, 0.178, 0.015, 73.8, 6.0
TAU0_TAIYI = 15.0        # Wang & Dunham's initial shear stress; "understressed" is below this

# validated with scripts/validate_palette.js --mode light: all six checks PASS
OBSC, MEAS = "#a8071a", "#6b6b66"
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
PASS_LO, PASS_HI = 0.85, 1.15        # +/-15%, the agreed tolerance

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "axes.edgecolor": MUTED, "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "legend.fontsize": 8.2,
})

# ---------------------------------------------------------------- sweep definitions
# Every sweep holds everything fixed but the named parameter. The "fixed" string is
# printed on the figure so a reader never has to trust that claim -- and params.txt
# in each run's folder has the full deck.
SWEEPS = [
    dict(key="muinit",
         title=r"Effect of initial friction $\mu_0$  (fault strength)",
         fixed=r"fixed: $\bar\sigma_0$=30.0 MPa, $\phi\beta$=1e-10, permev F, "
               r"Taiyi perm map, dc=1e-4",
         param=lambda d: f"$\\mu_0$={ffloat(d['muinit']):.2f}"
                         f"  ($\\tau_0$={ffloat(d['muinit'])*ffloat(d['sigmainit']):.1f} MPa)",
         runs=[632532, 632534, 632536, 632507],
         note="mu_0 is the ONLY difference. 632507 at mu=0.50 is Taiyi-equivalent "
              "stress (tau_0 15.0 MPa); the rest are understressed."),
    dict(key="muinit_permevT",
         title=r"Effect of $\mu_0$ with permeability enhancement ON",
         fixed=r"fixed: $\bar\sigma_0$=27.99 MPa, $\phi\beta$=1e-10, permev T, "
               r"kpmax=2e-11, kpmin=4e-13, kL=1e-5, Taiyi perm map, dc=1e-4",
         param=lambda d: f"$\\mu_0$={ffloat(d['muinit']):.2f}"
                         f"  ($\\tau_0$={ffloat(d['muinit'])*ffloat(d['sigmainit']):.1f} MPa)",
         runs=[632568, 632580, 632592],
         note="The enhancement-ON counterpart to sweep_muinit. Note 632533/632535/"
              "632537 look like the natural T twins of that figure (same sigma, same "
              "phi*beta) but leave kpmax AND kpmin unset, so permeability decays "
              "uncontrolled -- k reached 1.1e-16 in 632537. They are excluded."),
    dict(key="muinit_permevT_kmax1e-10",
         title=r"Effect of $\mu_0$ with enhancement ON at the highest defensible "
               r"kpmax",
         fixed=r"fixed: $\bar\sigma_0$=27.99 MPa, $\phi\beta$=1e-10, permev T, "
               r"kpmax=1e-10 ($\approx$100 D), kpmin=4e-13, kL=1e-5, Taiyi perm map",
         param=lambda d: f"$\\mu_0$={ffloat(d['muinit']):.2f}"
                         f"  ($\\tau_0$={ffloat(d['muinit'])*ffloat(d['sigmainit']):.1f} MPa)",
         runs=[632570, 632582, 632594],
         note="Same mu_0 sweep at kpmax 1e-10, the upper edge of a defensible "
              "fault-zone conduit, to check the mu_0 response is not specific to "
              "one enhancement ceiling."),
    dict(key="storage",
         title=r"Effect of storage $\phi\beta$ at understressed $\mu_0$=0.46",
         fixed=r"fixed: $\mu_0$=0.46, $\bar\sigma_0$=27.99 MPa ($\tau_0$=12.9 MPa, "
               r"understressed), permev F, Taiyi perm map",
         param=lambda d: f"$\\phi\\beta$={ffloat(d['phi'])*ffloat(d['beta']):.1e}",
         runs=[632700, 632702, 632704, 632720, 632721],
         note="Storage is the lever that moves the front: smaller phi*beta spreads "
              "pressure faster, so the front runs further on the same injected volume."),
    dict(key="storage_mu050",
         title=r"Effect of storage $\phi\beta$ at $\mu_0$=0.50 (Taiyi-equivalent stress)",
         fixed=r"fixed: $\mu_0$=0.50, $\bar\sigma_0$=27.99 MPa, permev F, Taiyi perm map",
         param=lambda d: f"$\\phi\\beta$={ffloat(d['phi'])*ffloat(d['beta']):.1e}",
         runs=[632706, 632708, 632710, 632722],
         note="Same storage sweep at a stronger fault, to separate the effect of "
              "storage from the effect of strength."),
    dict(key="enhancement",
         title="Effect of permeability enhancement (permev T vs F), matched pairs",
         fixed=r"fixed within each pair: $\mu_0$=0.46, $\bar\sigma_0$=27.99 MPa, "
               r"Taiyi perm map;  permev T uses kpmax=2e-11, kL=1e-5",
         param=lambda d: (f"$\\phi\\beta$={ffloat(d['phi'])*ffloat(d['beta']):.1e}, "
                          + ("permev T" if d.get("permev", "F").upper().startswith("T")
                             else "permev F")),
         runs=[632700, 632701, 632702, 632703, 632704, 632705],
         pairs=True,
         note="Enhancement is a NEGATIVE feedback here: the slip-formed high-k "
              "channel bleeds pressure away from the well, lowering peak p_f and "
              "shrinking the front rather than extending it."),
    dict(key="enhancement_at_match",
         title="Enhancement turned ON at the one parameter point that matched "
               "(632721)",
         fixed=r"fixed: $\mu_0$=0.46, $\bar\sigma_0$=27.99 MPa ($\tau_0$=12.88 MPa, "
               r"understressed), $\phi\beta$=5e-12, dc=1e-4, Taiyi perm map, "
               r"kpmin=4e-13.  Enhancement is the ONLY difference.",
         param=lambda d: ("permev F (fixed perm)"
                          if not d.get("permev", "F").upper().startswith("T")
                          else f"permev T, kpmax={d.get('kpmax','unset')}, "
                               f"kL={d.get('kL','?')}"),
         runs=[632721, 632752, 632751, 632750],
         note="The decisive test. Enhancement IMPROVES the pressure match "
              "(632752 reaches -0.8%) and destroys the front in the same move: "
              "the slip-opened high-k channel bleeds pressure away from the well, "
              "so effective normal stress never drops far enough to keep the fault "
              "failing. No enhancement setting keeps both targets."),
    dict(key="kpmax",
         title="Effect of enhanced permeability ceiling kpmax",
         fixed=r"fixed: $\mu_0$=0.37, $\bar\sigma_0$=27.99 MPa, $\phi\beta$=1e-10, "
               r"permev T, kL=1e-3, Taiyi perm map",
         param=lambda d: f"kpmax={d.get('kpmax','unset')}",
         runs=[632551, 632548, 632549, 632550],
         note="Raising the ceiling two orders of magnitude lowers the pressure the "
              "well can hold, which is the same negative feedback seen in the pairs."),
]


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


def fit(x, R):
    if len(x) < 2:
        return np.nan
    b = np.sqrt(x)
    return float(np.sum(b * R) / np.sum(b ** 2))


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


def run_data(n, dk):
    """Front (downdip) and wellhead pressure, using this run's own dc as threshold."""
    g = (glob.glob(f"/scratch/users/nberrios/3dhbi/runs/*/output/slip{n}.dat")
         + glob.glob(f"/scratch/users/nberrios/3dhbi/output/{n}/slip{n}.dat")
         + glob.glob(f"/oak/stanford/groups/edunham/nberrios/3doutput/{n}/slip{n}.dat"))
    if not g:
        return None
    p = g[0]
    IM, JM = int(dk["imax"]), int(dk["jmax"])
    ds, NC = ffloat(dk["ds"]), IM * JM
    t = np.atleast_2d(np.loadtxt(p.replace("slip", "time")))[:, 1] / 86400.0
    nt = min(os.path.getsize(p) // (8 * NC), len(t))
    if nt < 5:
        return None
    thr = ffloat(dk["dc"])
    s = np.memmap(p, np.float64, "r", shape=(nt, NC))
    cs = np.array([np.asarray(s[k]).reshape(IM, JM)[:, int(JM / 2)]
                   for k in range(nt)]).T
    x = np.linspace(-ds * IM / 2, ds * IM / 2, IM)
    R, T = [], []
    for i in range(IM):
        a = cs[i, :] > thr
        if a.any():
            R.append(abs(x[i])); T.append(t[np.argmax(a)])
    pwf = p.replace("slip", "pw")
    tp_ = pp_ = None
    if os.path.exists(pwf) and os.path.getsize(pwf):
        a = np.atleast_2d(np.loadtxt(pwf))
        f = (IN / dk["injection_file"]).read_text().split("\n")
        tq = np.array(f[2].split(), float); qq = np.array(f[4].split(), float) * W
        tp_ = a[:, 1] / 86400.0
        qi = np.interp(tp_ * 86400.0, tq, qq)
        pp_ = (P0 + a[:, 2]) - RHO * G * HW / 1e6 \
            + FD * (8.0 * HW * RHO * qi ** 2) / (np.pi ** 2 * DW ** 5) / 1e6
    return dict(R=np.array(R), T=np.array(T), t_end=t[nt - 1], thr=thr,
                tpw=tp_, ppw=pp_)


def pressure_score(obs, tpw, ppw, tcut):
    """Bias and RMS in MPa on flowing periods only -- see the module docstring."""
    hi = min(tpw[-1], obs["tp"].max(), tcut)
    if hi <= 0.2:
        return np.nan, np.nan
    gr = np.linspace(0.05, hi, 2000)
    ps = np.interp(gr, tpw, ppw)
    ob = np.interp(gr, obs["tp"], obs["pm"])
    qg = np.interp(gr, obs["ti"], obs["q"])
    fl = (qg > 0.25 * np.nanmax(obs["q"])) & (ob > 5.0)
    if fl.sum() < 30:
        return np.nan, np.nan
    d = ps[fl] - ob[fl]
    return float(d.mean()), float(np.sqrt(np.mean(d ** 2)))


def make(sw, obs):
    runs = []
    for n in sw["runs"]:
        if not (IN / f"res{n}.in").exists():
            continue
        dk = deck(n)
        d = run_data(n, dk)
        if d is not None:
            runs.append((n, dk, d))
    if not runs:
        print(f"  {sw['key']}: no runs with output, skipping")
        return None

    # one shared window for the whole sweep, and say which runs it excluded
    reaches = sorted((d["t_end"], n) for n, _, d in runs)
    tcut = max(r[0] for r in reaches)
    keep = [(n, dk, d) for n, dk, d in runs if d["t_end"] >= 0.95 * tcut]
    if len(keep) < 2:                       # too aggressive; fall back to the shortest
        tcut = min(r[0] for r in reaches)
        keep = runs
    dropped = [n for n, _, d in runs if d["t_end"] < 0.95 * tcut]

    vcut = float(np.interp(tcut, obs["ti"], obs["cumvol"] / 1e6))
    mo = obs["ft"] <= tcut
    lam_obs = fit(obs["ft"][mo], obs["fd"][mo])

    cmap = plt.get_cmap("viridis")
    cols = [cmap(x) for x in np.linspace(0.06, 0.86, max(len(keep), 1))]
    if sw.get("pairs"):     # colour by pair, dashed for the permev T member
        cmap = plt.get_cmap("tab10")
        cols = [cmap((i // 2) % 10) for i in range(len(keep))]

    fig, axes = plt.subplots(1, 3, figsize=(18.5, 6.0), constrained_layout=True)
    lines = []

    # ---------------- panel 1: wellhead pressure
    ax = axes[0]
    ax.plot(obs["tp"], obs["pm"], lw=1.1, color=MEAS, zorder=1,
            label="measured wellhead")
    for (n, dk, d), c in zip(keep, cols):
        if d["tpw"] is None:
            continue
        ls = "--" if (sw.get("pairs")
                      and dk.get("permev", "F").upper().startswith("T")) else "-"
        m = d["tpw"] <= tcut
        bias, rms = pressure_score(obs, d["tpw"], d["ppw"], tcut)
        ax.plot(d["tpw"][m], d["ppw"][m], lw=1.9, ls=ls, color=c, zorder=3,
                label=f"{sw['param'](dk)}\n  bias {bias:+.1f}, RMS {rms:.1f} MPa")
    style(ax).set(xlabel="days since injection began",
                  ylabel="absolute wellhead pressure (MPa)", xlim=(0, tcut * 1.02),
                  title="Wellhead pressure  (scored while flowing)")
    ax.legend(frameon=False, labelcolor=INK, loc="lower right")

    # ---------------- panels 2 and 3: R-T and R-V
    for ax, key, xlab, ttl in ((axes[1], "T", "days since injection began",
                                "Slip front vs time  (R–T)"),
                               (axes[2], "V", "cumulative injected volume (ML)",
                                "Slip front vs injected volume  (R–V)")):
        ex = obs["ev_t"] if key == "T" else np.interp(obs["ev_t"], obs["ti"],
                                                      obs["cumvol"] / 1e6)
        ax.scatter(ex, obs["ev_d"], s=3.5, alpha=0.18, color="gray",
                   edgecolors="none", zorder=1, label="observed events")
        fx = obs["ft"] if key == "T" else np.interp(obs["ft"], obs["ti"],
                                                   obs["cumvol"] / 1e6)
        lim = tcut if key == "T" else vcut
        mm = fx <= lim
        ax.scatter(fx[mm], obs["fd"][mm] + obs["org"], s=10, color=OBSC, zorder=3,
                   label="observed seismicity front")
        o = np.argsort(fx[mm])
        ax.plot(fx[mm][o], (lam_obs * np.sqrt(fx[mm]) + obs["org"])[o], color=OBSC,
                lw=2.4, zorder=4, label=rf"observed fit $\lambda$={lam_obs:.3f}")
        for (n, dk, d), c in zip(keep, cols):
            xx = d["T"] if key == "T" else np.interp(d["T"], obs["ti"],
                                                     obs["cumvol"] / 1e6)
            m = d["T"] <= tcut
            if m.sum() < 2:
                # the fault never slipped anywhere: pressure never reached
                # dp_crit = sigmainit*(1 - muinit/f0). Say so on the figure rather
                # than dropping the run and leaving a silent gap in the sweep.
                ax.plot([], [], lw=2.0, color=c,
                        label=f"{n}  {sw['param'](dk)}\n  NO SLIP anywhere")
                continue
            ls_ = fit(xx[m], d["R"][m])
            ls = "--" if (sw.get("pairs")
                          and dk.get("permev", "F").upper().startswith("T")) else "-"
            o = np.argsort(xx[m])
            ax.scatter(xx[m], d["R"][m], s=6, color=c, alpha=0.45, zorder=2,
                       edgecolors="none")
            ax.plot(xx[m][o], (ls_ * np.sqrt(xx[m]))[o], lw=2.0, ls=ls, color=c,
                    zorder=4,
                    label=f"{n}  {sw['param'](dk)}\n  "
                          rf"$\lambda$={ls_:.3f} ({ls_/lam_obs:.2f}$\times$ obs)")
        ax.set(xlabel=xlab, ylabel="distance from injector (km)",
               xlim=(0, lim * 1.02))
        style(ax).set_title(ttl)
        ax.legend(frameon=False, labelcolor=INK, loc="upper left")

    sub = (f"{sw['fixed']}\nall fits on 0–{tcut:.2f} d, observed front refit on the "
           f"same window; front threshold = each run's own dc")
    if dropped:
        sub += (f"\nexcluded (has not reached {tcut:.2f} d yet): "
                + ", ".join(str(x) for x in dropped))
    fig.suptitle(sw["title"] + "\n" + sub, fontsize=12.5)

    OUT.mkdir(parents=True, exist_ok=True)
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"sweep_{sw['key']}.{e}", dpi=180, bbox_inches="tight")
    plt.close(fig)

    # ---------------- the numbers, so the figure can be checked
    L = [sw["title"], "=" * 78, sw["fixed"].replace("$", "").replace("\\", ""),
         f"window 0-{tcut:.2f} d, observed lambda on this window = {lam_obs:.4f}", ""]
    if dropped:
        L += [f"excluded, too short: {dropped}", ""]
    L += [f"{'run':>8s} {'param':>26s} {'lambda':>8s} {'/obs':>6s} "
          f"{'p bias':>8s} {'p RMS':>7s} {'tau0':>6s} {'reached':>8s}", "-" * 86]
    for n, dk, d in keep:
        xx = d["T"]; m = xx <= tcut
        ls_ = fit(xx[m], d["R"][m]) if m.sum() >= 2 else np.nan
        bias, rms = (pressure_score(obs, d["tpw"], d["ppw"], tcut)
                     if d["tpw"] is not None else (np.nan, np.nan))
        tau0 = (ffloat(dk.get("muinit")) or 0) * (ffloat(dk.get("sigmainit")) or 0)
        plain = sw["param"](dk).replace("$", "").replace("\\", "")
        lam_s = "NO SLIP" if not np.isfinite(ls_) else f"{ls_:.4f}"
        rat_s = "  -   " if not np.isfinite(ls_) else f"{ls_/lam_obs:.2f}"
        L.append(f"{n:>8d} {plain:>26s} {lam_s:>8s} {rat_s:>6s} "
                 f"{bias:>+8.2f} {rms:>7.2f} {tau0:>6.2f} {d['t_end']:>8.2f}")
    L += ["", "NOTE: " + sw["note"]]
    (OUT / f"sweep_{sw['key']}.txt").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    print(f"  -> {OUT}/sweep_{sw['key']}.png\n")
    return tcut


def main():
    want = sys.argv[1:]
    obs = observed()
    OUT.mkdir(parents=True, exist_ok=True)
    for sw in SWEEPS:
        if want and sw["key"] not in want:
            continue
        make(sw, obs)


if __name__ == "__main__":
    main()
