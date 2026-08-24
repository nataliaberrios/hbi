#!/usr/bin/env python3
"""Comparison figures for the grid: one 3-panel row per sweep, plus two summaries.

Per sweep: pressure | R-T | R-V, the same three panels as the per-run figures,
overlaid across the sweep so the effect of one axis is read off a single row.

Conventions, all inherited from score_grid.py so the figures and the tables agree:
  - all fits on 0-5 d, observed front REFIT on that window (lambda_obs 0.1866).
    The same runs score 0.97-0.99 over 0-8 d, so the window must be stated.
  - pressure scored on FLOWING periods only, in MPa and as a percentage. HBI has
    no wellbore bleed-off, so it cannot follow the 1.58-3.56 d shut-in at all.
  - runs whose peak slip never reaches dc are labelled NO SLIP on the figure
    rather than dropped, because that is a result about dp_crit.

Two summaries:
  grid_tau0_summary   lambda/lambda_obs against tau_0, one line per
                      (map, sigmainit), with the +/-15% band. The Stage 3 story
                      in one panel.
  grid_plane          all 56 runs on the front-vs-wellhead plane, coloured by
                      fluid combination. Shows the trade directly: fluid A has
                      the front and not the pressure, fluid B the reverse.

Usage:  python make_grid_sweeps.py
"""
import importlib.util as iu
import json
import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

H = Path("/home/users/nberrios/3dhbi/hbi_analysis")
OUT = Path("/home/users/nberrios/3dhbi/hbi_git/docs/figs/cooper_grid/sweeps")
IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
WINDOW = 5.0

_s = iu.spec_from_file_location("sf", str(H / "make_sweep_figures.py"))
sf = iu.module_from_spec(_s)
_s.loader.exec_module(sf)

# same palette as every other figure in this project
OBSC, MEAS = "#a8071a", "#6b6b66"
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
OK = "#009E73"
plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "axes.edgecolor": MUTED, "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "legend.fontsize": 8.2,
})

MAP1807 = "perm_2zone_601_ds5_kmax2.5e-13.txt"
MAP1808 = "perm_2zone_601_ds5_kmax5e-14.txt"

SWEEPS = [
    dict(key="fluid_uniform",
         title="Fluid A vs B, uniform initial permeability",
         fixed=r"fixed: $\mu_0$=0.37, permev T, uniform kp = kpmin = 1e-15, "
               r"parameterfromfile F"
               "\nA = $\\eta$ 8.9e-4 with the parent's $\\beta$;  "
               "B = $\\eta$ 1.27e-4 with $\\beta\\times$7.007874 so D is unchanged",
         runs=[632800, 632802, 632804, 632806],
         dash=lambda d: sf.ffloat(d["eta"]) < 5e-4,
         note="A and B have IDENTICAL D by construction, so any difference here is "
              "str = beta*phi and the wellbore coupling gamma -- the two quantities "
              "the D-matching does not preserve."),
    dict(key="fluid_twozone",
         title="Fluid A vs B, two-zone initial permeability",
         fixed=r"fixed: $\mu_0$=0.37, permev T, near-well disc = kpmax, "
               r"background = kp = kpmin = 1e-15"
               "\ndashed = fluid B ($\\eta$ 1.27e-4, D-matched)",
         runs=[632812, 632810, 632818, 632816],
         dash=lambda d: sf.ffloat(d["eta"]) < 5e-4,
         note="Fluid B on the 250x map has NO front at all while its wellhead sits "
              "at -0.5% (RMS 1.4 MPa), the best pressure match this project has "
              "produced."),
    dict(key="sigma_twins",
         title=r"$\bar\sigma_0$ twins: 30.0 vs 27.99 MPa, uniform perm",
         fixed=r"fixed: $\mu_0$=0.37, permev T, uniform kp; each pair differs in "
               r"$\bar\sigma_0$ and nothing else"
               "\ndashed = 27.99 MPa (derived from the stress data)",
         runs=[632800, 632801, 632804, 632805],
         dash=lambda d: abs(sf.ffloat(d["sigmainit"]) - 27.99) < 0.1,
         note="dp_crit is 11.50 MPa at 30.0 and 10.73 at 27.99, against a peak "
              "measured overpressure of 10.92 -- so only the 27.99 twin can fail at "
              "the observed pressure."),
    dict(key="permev_twozone",
         title="Enhancement on vs off, two-zone initial permeability",
         fixed=r"fixed: $\mu_0$=0.37, $\bar\sigma_0$=27.99 MPa, fluid A; "
               r"dashed = permev F"
               "\nwith a two-zone field the conduit already exists, so permev F asks "
               "whether slip needs to extend it",
         runs=[632812, 632813, 632818, 632819],
         dash=lambda d: not d.get("permev", "F").upper().startswith("T"),
         note="Enhancement buys about 17% on the front and leaves the pressure "
              "essentially untouched."),
    dict(key="tau0_1807",
         title=r"$\tau_0$ sweep, 250$\times$ map, $\bar\sigma_0$=27.99 MPa",
         fixed=r"fixed: fluid B ($\eta$ 1.27e-4), permev T, "
               r"near-well 2.5e-13 / background 1e-15"
               "\n$\\mu_0$ derived as $\\tau_0/\\bar\\sigma_0$, 0.393 to 0.536",
         runs=list(range(632839, 632848)),
         note="NO SLIP at every tau_0 including Taiyi's own 15.0 MPa. The pressure "
              "is identical to three figures across the whole sweep, because with no "
              "slip the permeability never evolves and the hydraulics decouple from "
              "friction entirely."),
    dict(key="tau0_1808",
         title=r"$\tau_0$ sweep, 50$\times$ map, $\bar\sigma_0$=27.99 MPa",
         fixed=r"fixed: fluid B ($\eta$ 1.27e-4), permev T, "
               r"near-well 5e-14 / background 1e-15"
               "\n$\\mu_0$ derived as $\\tau_0/\\bar\\sigma_0$, 0.393 to 0.536",
         runs=list(range(632859, 632868)),
         note="lambda/lambda_obs climbs 0.03 -> 0.34 across tau_0 11.0 -> 15.0. "
              "Monotonic but weak: extrapolating to 1.0 needs tau_0 well above "
              "Taiyi's 15.0, i.e. a MORE stressed fault."),
]


def style(ax):
    ax.grid(True, color=GRID, lw=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return ax


def label(n, dk):
    tau0 = sf.ffloat(dk["muinit"]) * sf.ffloat(dk["sigmainit"])
    fl = "B" if sf.ffloat(dk["eta"]) < 5e-4 else "A"
    mp = dk.get("parameter_file", "uniform")
    con = ("uniform" if mp == "uniform"
           else "250x" if "2.5e-13" in mp else "50x")
    return (f"{n}  $\\tau_0$={tau0:.1f}, $\\bar\\sigma_0$={sf.ffloat(dk['sigmainit']):.2f}"
            f", fluid {fl}, {con}, permev {dk.get('permev','F')}")


def load(n):
    dk = sf.deck(n)
    d = sf.run_data(n, dk)
    if d is None:
        return None
    # peak slip decides whether a front exists at all
    for base in (f"/scratch/users/nberrios/3dhbi/output/{n}",):
        sp = f"{base}/slip{n}.dat"
        if os.path.exists(sp):
            IM, JM = int(dk["imax"]), int(dk["jmax"])
            NC = IM * JM
            nt = os.path.getsize(sp) // (8 * NC)
            with open(sp, "rb") as f:
                f.seek((nt - 1) * NC * 8)
                d["peak"] = float(np.frombuffer(f.read(NC * 8), np.float64).max())
    d["dk"] = dk
    return d


def pressure_score(obs, d):
    hi = min(d["tpw"][-1], obs["tp"].max(), WINDOW)
    gr = np.linspace(0.05, hi, 2000)
    ps = np.interp(gr, d["tpw"], d["ppw"])
    ob = np.interp(gr, obs["tp"], obs["pm"])
    qg = np.interp(gr, obs["ti"], obs["q"])
    fl = (qg > 0.25 * np.nanmax(obs["q"])) & (ob > 5.0)
    if fl.sum() < 30:
        return np.nan, np.nan
    dd = ps[fl] - ob[fl]
    return 100 * (ps[fl].mean() / ob[fl].mean() - 1), float(np.sqrt(np.mean(dd ** 2)))


def make_sweep(sw, obs, lam_obs):
    runs = [(n, d) for n in sw["runs"] if (d := load(n)) is not None]
    if len(runs) < 2:
        print(f"  {sw['key']}: <2 runs, skipped")
        return
    vcut = float(np.interp(WINDOW, obs["ti"], obs["cumvol"] / 1e6))
    cmap = plt.get_cmap("viridis")
    cols = [cmap(x) for x in np.linspace(0.06, 0.88, len(runs))]
    if "dash" in sw:                     # pairs: colour by pair, dash the B member
        cmap = plt.get_cmap("tab10")
        cols = [cmap((i // 2) % 10) for i in range(len(runs))]

    fig, axes = plt.subplots(1, 3, figsize=(18.5, 6.0), constrained_layout=True)

    ax = axes[0]
    ax.plot(obs["tp"], obs["pm"], lw=1.1, color=MEAS, zorder=1,
            label="measured wellhead")
    for (n, d), c in zip(runs, cols):
        if d["tpw"] is None:
            continue
        ls = "--" if ("dash" in sw and sw["dash"](d["dk"])) else "-"
        m = d["tpw"] <= WINDOW
        pct, rms = pressure_score(obs, d)
        ax.plot(d["tpw"][m], d["ppw"][m], lw=1.9, ls=ls, color=c, zorder=3,
                label=f"{n}  {pct:+.1f}%, RMS {rms:.1f} MPa")
    style(ax).set(xlabel="days since injection began",
                  ylabel="absolute wellhead pressure (MPa)",
                  xlim=(0, WINDOW * 1.02),
                  title="Wellhead pressure  (scored while flowing)")
    ax.legend(frameon=False, labelcolor=INK, loc="upper left")

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
        lim = WINDOW if key == "T" else vcut
        mm = fx <= lim
        ax.scatter(fx[mm], obs["fd"][mm] + obs["org"], s=10, color=OBSC, zorder=3,
                   label="observed seismicity front")
        o = np.argsort(fx[mm])
        ax.plot(fx[mm][o], (lam_obs * np.sqrt(fx[mm]) + obs["org"])[o], color=OBSC,
                lw=2.4, zorder=4, label=rf"observed fit $\lambda$={lam_obs:.3f}")
        for (n, d), c in zip(runs, cols):
            ls = "--" if ("dash" in sw and sw["dash"](d["dk"])) else "-"
            xx = d["T"] if key == "T" else np.interp(d["T"], obs["ti"],
                                                     obs["cumvol"] / 1e6)
            m = d["T"] <= WINDOW
            if m.sum() < 2 or d.get("peak", 9) < d["thr"]:
                ax.plot([], [], lw=2.0, ls=ls, color=c,
                        label=f"{n}  NO SLIP (peak slip < dc)")
                continue
            ls_ = sf.fit(xx[m], d["R"][m])
            o = np.argsort(xx[m])
            ax.scatter(xx[m], d["R"][m], s=6, color=c, alpha=0.45, zorder=2,
                       edgecolors="none")
            ax.plot(xx[m][o], (ls_ * np.sqrt(xx[m]))[o], lw=2.0, ls=ls, color=c,
                    zorder=4,
                    label=rf"{n}  $\lambda$={ls_:.3f} ({ls_/lam_obs:.2f}$\times$)")
        ax.set(xlabel=xlab, ylabel="distance from injector (km)", xlim=(0, lim * 1.02))
        style(ax).set_title(ttl)
        ax.legend(frameon=False, labelcolor=INK, loc="upper left")

    fig.suptitle(f"{sw['title']}\n{sw['fixed']}\nall fits on 0–{WINDOW:.0f} d, "
                 f"observed front refit on the same window", fontsize=12)
    OUT.mkdir(parents=True, exist_ok=True)
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"sweep_{sw['key']}.{e}", dpi=170, bbox_inches="tight")
    plt.close(fig)
    (OUT / f"sweep_{sw['key']}.txt").write_text(
        sw["title"] + "\n" + "=" * 76 + "\n"
        + sw["fixed"].replace("$", "").replace("\\", "") + "\n\n"
        + f"window 0-{WINDOW:.0f} d, lambda_obs {lam_obs:.4f}\n\n"
        + "\n".join(
            f"  {n}  lam/obs "
            + ("NO SLIP" if d.get("peak", 9) < d["thr"]
               else f"{sf.fit(d['T'][d['T']<=WINDOW], d['R'][d['T']<=WINDOW])/lam_obs:.3f}")
            + f"   pressure {pressure_score(obs,d)[0]:+.1f}%"
            for n, d in runs)
        + "\n\nNOTE: " + sw["note"] + "\n")
    print(f"  wrote sweeps/sweep_{sw['key']}.png")


def summaries(obs, lam_obs):
    scores = json.loads((OUT.parent / "grid_scores.json").read_text())
    by = {r["n"]: r for r in scores}

    # ---- lambda vs tau_0, one line per (map, sigma)
    groups = [
        ("250x map, $\\bar\\sigma_0$=30.0", list(range(632830, 632839)), "#0072BD", "-"),
        ("250x map, $\\bar\\sigma_0$=27.99", list(range(632839, 632848)), "#0072BD", "--"),
        ("50x map, $\\bar\\sigma_0$=30.0", list(range(632850, 632859)), "#D55E00", "-"),
        ("50x map, $\\bar\\sigma_0$=27.99", list(range(632859, 632868)), "#D55E00", "--"),
    ]
    fig, ax = plt.subplots(figsize=(9.4, 6.2), constrained_layout=True)
    ax.axhspan(0.85, 1.15, color=OK, alpha=0.12, zorder=0)
    ax.axhline(1.0, color=OBSC, lw=1.6, ls="--", zorder=1)
    ax.annotate("observed front", xy=(0.99, 1.0), xycoords=("axes fraction", "data"),
                xytext=(0, 6), textcoords="offset points", ha="right",
                fontsize=9.5, color=OBSC)
    for lab, runs, c, ls in groups:
        t, y = [], []
        for n in runs:
            r = by.get(n)
            if not r:
                continue
            t.append(r["tau0"])
            y.append(r.get("lam_ratio", 0.0) if r.get("peak_over_dc", 0) >= 1 else 0.0)
        if t:
            ax.plot(t, y, marker="o", lw=2, ms=6, color=c, ls=ls, label=lab)
    ax.axvline(15.0, color=MUTED, lw=1.1, ls=":")
    ax.annotate("Taiyi $\\tau_0$=15.0", xy=(15.0, 0.02),
                xycoords=("data", "axes fraction"), xytext=(-6, 0),
                textcoords="offset points", ha="right", fontsize=9, color=MUTED)
    ax.set(xlabel=r"initial shear stress $\tau_0$ (MPa)",
           ylabel=r"$\lambda\,/\,\lambda_{observed}$      0 = no slip at all",
           ylim=(-0.03, 1.25))
    ax.set_title("Stage 3: the front against shear stress, fluid B (correct viscosity)\n"
                 "the 250x map never slips; the 50x map reaches 0.34 at Taiyi's own "
                 r"$\tau_0$", pad=10)
    style(ax).legend(frameon=False, labelcolor=INK, loc="upper left")
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"grid_tau0_summary.{e}", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("  wrote sweeps/grid_tau0_summary.png")

    # ---- front vs wellhead plane, all runs
    fig, ax = plt.subplots(figsize=(10.0, 7.0), constrained_layout=True)
    ax.axhspan(0.85, 1.15, color=OK, alpha=0.10, zorder=0)
    ax.axvspan(-15, 15, color=OK, alpha=0.10, zorder=0)
    ax.axhline(1.0, color=OBSC, lw=1.4, ls="--", zorder=1)
    ax.axvline(0.0, color=OBSC, lw=1.4, ls="--", zorder=1)
    ax.plot([0], [1], marker="*", ms=26, color=OBSC, zorder=8)
    ax.annotate("TARGET", xy=(0, 1), xytext=(13, -6), textcoords="offset points",
                fontsize=11, color=OBSC, weight="bold")
    seen = set()
    for r in scores:
        if "p_pct" not in r:
            continue
        fl = "B" if r["eta"] < 5e-4 else "A"
        c = "#D55E00" if fl == "A" else "#0072BD"
        mk = "o" if r["perm"] == "uniform" else ("D" if "2.5e-13" in str(r["perm"])
                                                 else "s")
        y = r.get("lam_ratio", 0.0) if r.get("peak_over_dc", 0) >= 1 else 0.0
        ax.scatter([r["p_pct"]], [y], marker=mk, s=80, color=c, alpha=0.85,
                   edgecolors=INK, linewidths=0.5, zorder=5)
        seen.add((fl, mk))
    ax.set_xscale("symlog", linthresh=20)
    ax.set(xlabel="wellhead pressure error while flowing (%)      0 = matched",
           ylabel=r"$\lambda\,/\,\lambda_{observed}$      0 = no slip",
           ylim=(-0.05, 1.3))
    ax.axvline(100, color=MUTED, lw=1, ls=":")
    ax.annotate("wellhead 2x measured", xy=(100, 0.42), xycoords=("data", "data"),
                xytext=(5, 0), textcoords="offset points", fontsize=8.5,
                color=MUTED, rotation=90, va="center")
    ax.set_title("Every grid run on both targets at once\n"
                 "fluid A (orange) has the front and not the pressure; "
                 "fluid B (blue) the reverse", pad=10)
    style(ax)
    ax.legend(handles=[
        Line2D([], [], marker="*", ls="", ms=16, color=OBSC, label="target"),
        Line2D([], [], marker="o", ls="", ms=9, color="#D55E00",
               label=r"fluid A, $\eta$ 8.9e-4 (Taiyi)"),
        Line2D([], [], marker="o", ls="", ms=9, color="#0072BD",
               label=r"fluid B, $\eta$ 1.27e-4 (corrected)"),
        Line2D([], [], marker="o", ls="", ms=9, color=MUTED, label="uniform perm"),
        Line2D([], [], marker="D", ls="", ms=9, color=MUTED, label=r"two-zone 250$\times$"),
        Line2D([], [], marker="s", ls="", ms=9, color=MUTED, label=r"two-zone 50$\times$"),
    ], frameon=False, labelcolor=INK, loc="upper right")
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"grid_plane.{e}", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("  wrote sweeps/grid_plane.png")


def main():
    obs = sf.observed()
    mo = obs["ft"] <= WINDOW
    lam_obs = sf.fit(obs["ft"][mo], obs["fd"][mo])
    print(f"window 0-{WINDOW:.0f} d, lambda_obs {lam_obs:.4f}\n")
    for sw in SWEEPS:
        make_sweep(sw, obs, lam_obs)
    summaries(obs, lam_obs)


if __name__ == "__main__":
    main()
