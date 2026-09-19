#!/usr/bin/env python3
"""PAPER FIGURE: the seismicity front is blind to the permeability that sets
total slip.

Seven runs at tau_0 = 11.53 MPa, disc 300 m, with kpmax swept over a 31-fold
range of contrast kpmax/kpmin (13x to 400x). ONLY kpmax moves -- stress state,
disc radius, fluid, injection history and mesh are identical -- so nothing but
the permeability can explain any difference between the curves.

    over a 31-fold range of permeability contrast
      lambda/lambda_obs   spans 1.01x        (1.02 1.03 1.03 1.03 1.03 1.02)
      R/R_obs at 8.7 d    spans 1.03x        (0.86 to 0.89)
      TOTAL SLIP          spans 12.15x       (4.06x down to 0.33x)
      wellhead dp         spans 533 points   (+527% down to -6%)

The seismicity front is constant to ONE PERCENT -- inside the scatter of the
observed front itself -- while total slip varies by an order of magnitude. The
front does not constrain total slip, and the pressure history that separates
these models is not a refinement but a necessity.

PANEL (d) IS THE POINT, AND IT HAS ONE AXIS. All three observables are plotted
as model/observed, so a perfect match is the dashed line at 1 and the three
curves can share a log axis. lambda lies ON that line across the whole sweep;
slip crosses it near 65x; pressure crosses it near 290x. That the crossings sit
at DIFFERENT contrasts is the paper's open question, and a reader should see it
without being told.

An earlier version put dp on a twin axis in percent, which needed a second
coloured spine and coloured ticks and still did not let the reader compare the
three. Expressing dp as a ratio removes the twin axis entirely.

DESIGN. Journal register, not slide register: 190 mm AGU two-column width, bare
bold panel letters and NO descriptive titles (the caption does that work), no
grid, no in-panel annotations, top and right spines off, one legend per panel at
most. Observed data is the same dark red in every panel and always the thickest
line, so it is identified once and read everywhere.

DEFINITIONS ARE IMPORTED, NOT REIMPLEMENTED. The seismicity front is
compare_arms.seismicity_front (events sorted by time, binned 100 at a time, the
90th-95th percentile band of each bin's distance, initial cloud radius
subtracted); lambda is compare_arms.lam, least squares through the origin; the
pressure datum is fault_pressure.datum(). This figure cannot disagree with the
scorer about what any of those mean.

TWO CURVES LEAVE PANEL (b). At 13x and 30x contrast the model reaches 73 and
37 MPa; the axis stops at 28 so the matched cases stay legible. That belongs in
the caption rather than in a rescale that would flatten the plateau the figure
exists to show.

Usage:
  python fig_tradeoff.py                 # tau_0 = 11.53, the paper's row
  python fig_tradeoff.py --tau 10.36
"""
import argparse
import glob
import importlib.util as iu
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.cm import ScalarMappable

H = "/home/users/nberrios/3dhbi/hbi_analysis"
_s = iu.spec_from_file_location("sf", H + "/make_sweep_figures.py")
sf = iu.module_from_spec(_s); _s.loader.exec_module(sf)
_a = iu.spec_from_file_location(
    "ca", "/home/users/nberrios/3dhbi/hbi_git/docs/cooper_basin/compare_arms.py")
ca = iu.module_from_spec(_a); _a.loader.exec_module(ca)
_f = iu.spec_from_file_location(
    "fp", "/home/users/nberrios/3dhbi/hbi_git/docs/cooper_basin/fault_pressure.py")
fp = iu.module_from_spec(_f); _f.loader.exec_module(fp)
sys.path.insert(0, H + "/notebooks")
from sim_curves import load_slip

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
OUTD = Path("/scratch/users/nberrios/3dhbi/output")
FIG = Path(H) / "figures" / "cycle2"
OBS_SLIP = Path("/home/users/nberrios/3dhbi/hbi/slip_profiles_strike.txt")
OT = [3, 5, 7, 9, 11, 13, 15, 17]
T0, TMAX, T_EVAL = 4.300, 13.2, 8.7
DP_TOP = 28.0
INK, OBSC = "#141413", "#a8071a"
CB_TICKS = [13, 30, 65, 145, 250, 400]

plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 9, "legend.fontsize": 8,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "axes.edgecolor": INK, "axes.linewidth": 0.7,
    "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": INK, "ytick.color": INK,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 3.0, "ytick.major.size": 3.0,
    "axes.grid": False, "legend.frameon": False,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def panel(ax, letter):
    """Bold letter above the axes, outside the data area. No title."""
    ax.text(-0.17, 1.02, f"({letter})", transform=ax.transAxes,
            fontsize=10, fontweight="bold", va="bottom", ha="left")


def members(tau):
    """[(contrast, run)] for the kpmax sweep at this tau_0, completed only."""
    out = []
    for f in sorted(glob.glob(str(IN / "res6329[0-9][0-9].in"))
                    + glob.glob(str(IN / "res6330[0-9][0-9].in"))
                    + glob.glob(str(IN / "res6331[0-9][0-9].in"))):
        n = int(re.search(r"res(\d+)", f).group(1))
        if not (OUTD / str(n) / f"slip{n}.dat").exists():
            continue
        d = sf.deck(n)
        try:
            if abs(sf.ffloat(d["muinit"]) * sf.ffloat(d["sigmainit"]) - tau) > 0.02:
                continue
            if "d4300" not in d["injection_file"]: continue
            if d.get("permev", "T").strip() != "T": continue
            m = re.search(r"disc(\d+)", d["parameter_file"])
            if not m or int(m.group(1)) != 300: continue
            # THE ARM-2 FLUID, EXACTLY. Without beta and Sw_fwid the storage
            # sweep contaminates the row: 632981 (phi varied) entered at
            # lambda 1.86x and 632971 at 0.08x with zero slip, and 633069-71
            # (Sw_fwid varied) duplicated 632961's 250x. One stray run at a
            # repeated contrast made panel (d) non-monotonic and inflated the
            # quoted spreads to lambda 23x and slip 44848x.
            if abs(sf.ffloat(d["eta"]) - 1.27e-4) > 1e-12: continue
            if abs(sf.ffloat(d["phi"]) - 0.01) > 1e-12: continue
            if abs(sf.ffloat(d["beta"]) - 2.25e-8) > 1e-12: continue
            if abs(sf.ffloat(d["Sw_fwid"]) - 7.4e-9) > 1e-13: continue
            if "skin" in d: continue
            c = sf.ffloat(d["kpmax"]) / sf.ffloat(d["kpmin"])
        except KeyError:
            continue
        out.append((round(c), n))
    seen = {}
    for c, n in sorted(out):
        seen.setdefault(c, n)          # one run per contrast
    return sorted(seen.items())


def slip_series(n):
    """(t_days, slip_cm) at the injector, on the run's own output cadence."""
    t, v = [], []
    for t_ in np.linspace(0.2, 13.0, 40):
        try:
            _, sl, ta = load_slip(n, t_, how="strike")
            if abs(ta - t_) < 0.3:
                t.append(t_); v.append(sl[0] * 100)
        except Exception:
            pass
    return np.array(t), np.array(v)


def crossing(cs, y):
    """Contrast at which the model/observed ratio passes 1, log-interpolated."""
    lc, ly = np.log(np.asarray(cs, float)), np.log(np.asarray(y, float))
    o = np.argsort(ly)
    if not (ly.min() < 0 < ly.max()):
        return np.nan
    return float(np.exp(np.interp(0.0, ly[o], lc[o])))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tau", type=float, default=11.53)
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)

    runs = members(a.tau)
    if len(runs) < 3:
        sys.exit(f"only {len(runs)} completed runs at tau_0 = {a.tau}")

    obs = sf.observed()
    tcat, rcat, _ = ca.fb.catalogue()
    k = (tcat - T0) > 0
    te, re_ = tcat[k] - T0, rcat[k]
    org = float(np.median(re_[np.argsort(te)][:10]))
    ft, fd = ca.seismicity_front(te, re_ - org)
    LOBS = ca.lam(ft, fd)
    R_OBS = float(np.interp(T_EVAL, ft, np.maximum.accumulate(fd + org)))
    to, dpo = fp.dp_observed(obs)
    o = np.loadtxt(OBS_SLIP)
    peak = [o[:, 1 + i].max() for i in range(8)]
    s0 = float(np.interp(T0, OT, peak))
    S_OBS = float(np.interp(T_EVAL + T0, OT, peak)) - s0
    P_STATIC = sf.P0 - sf.RHO * sf.G * sf.HW / 1e6

    fig, ax = plt.subplots(2, 2, figsize=(7.48, 5.9), dpi=400,
                           constrained_layout=True)
    (ar, ap_), (asl, am) = ax
    cs = [c for c, _ in runs]
    norm = LogNorm(vmin=min(cs) / 1.35, vmax=max(cs) * 1.35)
    cmap = plt.cm.viridis
    tf = np.linspace(0.0, TMAX, 300)

    ar.scatter(ft, fd + org, s=3.5, color=OBSC, alpha=0.28, lw=0, zorder=1)
    ar.plot(tf, LOBS * np.sqrt(tf) + org, "-", lw=2.0, color=OBSC, zorder=5,
            label="observed")
    ap_.plot(to, dpo, lw=1.0, color=OBSC, alpha=0.9, zorder=5)
    tg = np.linspace(0.2, 13.1, 80)
    asl.plot(tg, np.interp(tg + T0, OT, peak) - s0, lw=2.0, color=OBSC, zorder=5)

    L, DPR, SL = [], [], []
    print(f"tau_0 = {a.tau} MPa, disc 300 m, {len(runs)} runs")
    print(f"observed: lambda {LOBS:.1f} m/sqrt(d), R({T_EVAL} d) {R_OBS:.0f} m, "
          f"slip {S_OBS:.3f} cm, datum {fp.datum():.3f} MPa")
    print(f"{'run':>7} {'contrast':>9} {'lam/obs':>8} {'R/R_obs':>8} "
          f"{'dp %':>8} {'dp/obs':>7} {'slip':>6}")
    for c, n in runs:
        col = cmap(norm(c))
        d = sf.run_data(n, sf.deck(n))
        oo = np.argsort(d["T"])
        T, R = np.asarray(d["T"])[oo], np.asarray(d["R"])[oo] * 1000.0
        lam = ca.lam(T, R)
        Rev = float(np.interp(T_EVAL, T, R))
        ar.plot(T, R, lw=1.1, color=col, zorder=3)

        dpm = d["ppw"] - P_STATIC
        ap_.plot(d["tpw"], np.where(dpm > DP_TOP, np.nan, dpm),
                 lw=1.1, color=col, zorder=3)
        gr = np.linspace(0.05, min(d["tpw"][-1], to.max(), T_EVAL), 2000)
        ps = np.interp(gr, d["tpw"], dpm)
        ob = np.interp(gr, to, dpo)
        qg = np.interp(gr, obs["ti"] - T0, obs["q"])
        fl = (qg > 0.25 * np.nanmax(obs["q"])) \
            & (np.interp(gr, obs["tp"] - T0, obs["pm"]) > 5.0)
        # dp AS A RATIO, so panel (d) needs no twin axis
        dpr = float(np.mean(ps[fl])) / float(np.mean(ob[fl]))

        st, sv = slip_series(n)
        asl.plot(st, sv, lw=1.1, color=col, zorder=3)
        sr = float(np.interp(T_EVAL, st, sv)) / S_OBS if len(st) else np.nan

        L.append(lam / LOBS); DPR.append(dpr); SL.append(sr)
        print(f"{n:>7} {c:>8.0f}x {lam/LOBS:>8.2f} {Rev/R_OBS:>8.2f} "
              f"{100*(dpr-1):>+8.1f} {dpr:>7.2f} {sr:>6.2f}")

    # --- panel (d): one axis, three observables, a match is the line at 1
    am.axhline(1.0, color=INK, lw=0.8, ls=(0, (4, 3)), zorder=1)
    mk = dict(ms=4.0, mew=0.6, mec=INK, lw=1.5, zorder=3)
    am.plot(cs, L, "o-", color="#1d4ed8", label="seismicity front", **mk)
    am.plot(cs, SL, "s-", color="#D55E00", label="total slip", **mk)
    am.plot(cs, DPR, "^-", color="#0e7a63", label="wellhead pressure", **mk)
    am.set_xscale("log"); am.set_yscale("log")
    am.set(xlabel=r"Permeability contrast, $k_{p,\max}/k_{p,\min}$",
           ylabel="Model / observed", ylim=(0.2, 9))
    am.set_xticks(CB_TICKS)
    am.set_xticklabels([str(t) for t in CB_TICKS])
    am.minorticks_off()
    am.set_yticks([0.25, 0.5, 1, 2, 4, 8])
    am.set_yticklabels(["0.25", "0.5", "1", "2", "4", "8"])
    am.legend(loc="upper right", handlelength=1.8, borderaxespad=0.2,
              labelspacing=0.35)

    ar.set(xlabel="Time since injection resumed (d)",
           ylabel="Front radius (m)", xlim=(0, TMAX), ylim=(0, 1150))
    ar.legend(loc="lower right", handlelength=1.6)
    ap_.set(xlabel="Time since injection resumed (d)",
            ylabel=r"Pressure change from $p_{f0}$ (MPa)",
            xlim=(0, TMAX), ylim=(0, DP_TOP))
    asl.set(xlabel="Time since injection resumed (d)",
            ylabel="Slip at the injector (cm)", xlim=(0, TMAX), ylim=(0, 24))
    for x, l in zip((ar, ap_, asl, am), "abcd"):
        panel(x, l)

    sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    # SLIM AND FULL HEIGHT. A short colorbar centred on the figure reads as a
    # floating fifth element; spanning the full height makes it furniture.
    cb = fig.colorbar(sm, ax=ax.ravel().tolist(), location="right",
                      shrink=0.94, aspect=46, pad=0.012, ticks=CB_TICKS)
    cb.ax.set_yticklabels([str(t) for t in CB_TICKS])
    cb.ax.minorticks_off()
    cb.ax.tick_params(length=2.5, width=0.6)
    cb.outline.set_linewidth(0.6)
    cb.set_label(r"$k_{p,\max}/k_{p,\min}$", fontsize=8.5, labelpad=4)

    print(f"\nover a {max(cs)/min(cs):.0f}-fold contrast range: "
          f"lambda spans {max(L)/min(L):.2f}x, slip spans {max(SL)/min(SL):.2f}x, "
          f"dp spans {max(DPR)/min(DPR):.2f}x")
    print(f"model/observed = 1 crossings:  slip {crossing(cs, SL):.0f}x, "
          f"pressure {crossing(cs, DPR):.0f}x, front never "
          f"(flat at {np.mean(L):.2f})")

    stem = a.out or f"fig_tradeoff_tau{int(round(a.tau*100))}"
    FIG.mkdir(parents=True, exist_ok=True)
    for e_ in ("png", "pdf"):
        fig.savefig(FIG / f"{stem}.{e_}")
    plt.close(fig)
    print(f"wrote {FIG}/{stem}.png")


if __name__ == "__main__":
    main()
