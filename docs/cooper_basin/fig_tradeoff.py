#!/usr/bin/env python3
"""PAPER FIGURE: the seismicity front is blind to the permeability that sets
total slip.

This is the figure the paper is built on. Everything else either sets it up
(the model description, the analytical solution) or follows from it (the joint
constraint, the slip partitioning).

ONE ROW OF THE GRID. tau_0 = 11.53 MPa, disc 300 m, and kpmax swept over a
31-fold range of contrast kpmax/kpmin (13x to 400x). Only kpmax moves; the
stress state, the disc, the fluid, the injection history and the mesh are
identical across the seven runs, so nothing but the permeability can explain any
difference between them.

WHAT THE THREE PANELS SAY

  (a)  the seven model fronts lie on top of each other and on the data
  (b)  the seven wellhead pressure histories fan across the whole axis
  (c)  lambda/lambda_obs is flat while total slip falls monotonically

The point is the CONTRAST between (a) and (b), so they share a time axis and are
drawn at the same scale. A reader should be able to see that the curves which
are indistinguishable in (a) are obviously distinguishable in (b).

Panel (c) is the quantitative version, and the numbers are the result:

    over a 31-fold range of permeability contrast,
      lambda/lambda_obs spans 1.01x        (1.02, 1.03, 1.03, 1.03, 1.03, 1.02)
      R/R_obs at 8.7 d   spans 1.03x       (0.86 to 0.89)
      TOTAL SLIP         spans 12.15x      (4.06x down to 0.33x)
      wellhead dp        spans 533 points  (+527% down to -6%)

So the seismicity front is constant to ONE PERCENT -- far inside the scatter of
the observed front itself -- while total slip varies by an order of magnitude.
The front does not constrain total slip at all, and the pressure history that
separates these models is not a refinement but a necessity.

THE TWO TARGETS DO NOT AGREE, and the figure shows where. Total slip matches the
catalogue at contrast 65x, where dp is +89%; dp matches at 250-320x, where slip
is 0.36-0.40x. Either the crack-model inversion of the catalogue overestimates
slip by ~2.5x, or the model is missing slip. That is the paper's open question
and panel (c) is where a reader should be able to see it.

THE STRONGEST FORM OF THIS IS NOT IN THE SWEEP. Comparing 633001 with 633120,
where permeability EVOLUTION is off and the background is 1.76x higher, the
front moves 14% (lambda 0.99x -> 0.85x) and total slip moves a factor of TEN
(0.42x -> 4.06x). That pair is not plotted here because 633120 changes the model
structure rather than one parameter, so it does not belong on a one-parameter
axis; it is quoted in the caption and shown in static_perm_compare.

DEFINITIONS ARE IMPORTED, NOT REIMPLEMENTED. The seismicity front is
compare_arms.seismicity_front (events sorted by time, binned 100 at a time,
the 90th-95th percentile band of each bin's distance, initial cloud radius
subtracted); lambda is compare_arms.lam, least squares through the origin; the
pressure datum is fault_pressure.datum(). A figure from this script cannot
disagree with the scorer about what any of those mean.

TWO FRONT NUMBERS EXIST AND BOTH ARE PRINTED. lambda/lambda_obs is the slope of
the sqrt(t) fit; R/R_obs at 8.7 d is the radius. They differ (0.99 vs 0.88 for
633001) because the model grows diffusively at the right rate from a cloud that
starts too small. Quoting one as though it were the other has misled this
project twice.

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
INK, MUTED, OBSC = "#1a1a19", "#6b6b66", "#a8071a"

plt.rcParams.update({"font.size": 9.5, "axes.titlesize": 10.5,
                     "axes.labelsize": 10, "axes.edgecolor": MUTED,
                     "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED,
                     "xtick.labelcolor": INK, "ytick.labelcolor": INK,
                     "legend.fontsize": 8, "legend.frameon": False})


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
            # sweep contaminates the row: 632981 (phi varied) came in at
            # lambda 1.86x and 632971 at 0.08x with zero slip, and 633069-71
            # (Sw_fwid varied) duplicated 632961's contrast. One stray run at a
            # repeated contrast makes panel (c) non-monotonic and inflates the
            # quoted spreads by orders of magnitude.
            if abs(sf.ffloat(d["eta"]) - 1.27e-4) > 1e-12: continue
            if abs(sf.ffloat(d["phi"]) - 0.01) > 1e-12: continue
            if abs(sf.ffloat(d["beta"]) - 2.25e-8) > 1e-12: continue
            if abs(sf.ffloat(d["Sw_fwid"]) - 7.4e-9) > 1e-13: continue
            if "skin" in d: continue
            c = sf.ffloat(d["kpmax"]) / sf.ffloat(d["kpmin"])
        except KeyError:
            continue
        out.append((round(c), n))
    # ONE RUN PER CONTRAST. Duplicates would plot two points at the same x.
    seen = {}
    for c, n in sorted(out):
        seen.setdefault(c, n)
    return sorted(seen.items())


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

    fig = plt.figure(figsize=(13.2, 4.3), dpi=300)
    # A DEDICATED COLORBAR COLUMN. Attaching the colorbar to [ar, ap_] with
    # fraction/pad put its label on top of panel (c)'s y-label.
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 0.045, 1], wspace=0.40,
                          left=0.05, right=0.97, bottom=0.135, top=0.90)
    ar = fig.add_subplot(gs[0, 0])
    ap_ = fig.add_subplot(gs[0, 1])
    cax = fig.add_subplot(gs[0, 2])
    am = fig.add_subplot(gs[0, 3])

    cs = [c for c, _ in runs]
    norm = LogNorm(vmin=min(cs) * 0.8, vmax=max(cs) * 1.25)
    cmap = plt.cm.viridis
    tf = np.linspace(0.0, TMAX, 300)

    ar.scatter(ft, fd + org, s=7, color=OBSC, alpha=0.45, lw=0, zorder=1,
               label="observed")
    ar.plot(tf, LOBS * np.sqrt(tf) + org, "-", lw=2.2, color=OBSC, zorder=4,
            label=r"$\lambda_{\rm obs}$ = " f"{LOBS:.0f} m/" r"$\sqrt{\rm d}$")
    ap_.plot(to, dpo, lw=1.2, color=OBSC, alpha=0.85, zorder=5,
             label="observed")

    L, DP, SL = [], [], []
    print(f"tau_0 = {a.tau} MPa, disc 300 m, {len(runs)} runs")
    print(f"observed: lambda {LOBS:.1f} m/sqrt(d), R({T_EVAL} d) {R_OBS:.0f} m, "
          f"slip {S_OBS:.3f} cm, datum {fp.datum():.3f} MPa")
    print(f"{'run':>7} {'contrast':>9} {'lam/obs':>8} {'R/R_obs':>8} "
          f"{'dp %':>7} {'slip':>6}")
    for c, n in runs:
        col = cmap(norm(c))
        d = sf.run_data(n, sf.deck(n))
        oo = np.argsort(d["T"])
        T, R = np.asarray(d["T"])[oo], np.asarray(d["R"])[oo] * 1000.0
        lam = ca.lam(T, R)
        Rev = float(np.interp(T_EVAL, T, R))
        ar.plot(T, R, lw=1.5, color=col, alpha=0.95, zorder=3)

        dpm = d["ppw"] - P_STATIC
        ap_.plot(d["tpw"], dpm, lw=1.5, color=col, alpha=0.95, zorder=3)
        gr = np.linspace(0.05, min(d["tpw"][-1], to.max(), T_EVAL), 2000)
        ps = np.interp(gr, d["tpw"], dpm)
        ob = np.interp(gr, to, dpo)
        qg = np.interp(gr, obs["ti"] - T0, obs["q"])
        fl = (qg > 0.25 * np.nanmax(obs["q"])) \
            & (np.interp(gr, obs["tp"] - T0, obs["pm"]) > 5.0)
        e = 100.0 * float(np.mean(ps[fl] - ob[fl])) / float(np.mean(ob[fl]))

        sv, st = [], []
        for t_ in np.linspace(0.2, 13.0, 40):
            try:
                _, sl, ta = load_slip(n, t_, how="strike")
                if abs(ta - t_) < 0.3:
                    st.append(t_); sv.append(sl[0] * 100)
            except Exception:
                pass
        sr = float(np.interp(T_EVAL, st, sv)) / S_OBS if st else np.nan

        L.append(lam / LOBS); DP.append(e); SL.append(sr)
        print(f"{n:>7} {c:>8.0f}x {lam/LOBS:>8.2f} {Rev/R_OBS:>8.2f} "
              f"{e:>+7.1f} {sr:>6.2f}")

    # --- panel (c): the two curves that matter, on one log axis
    am.axhline(1.0, color=MUTED, lw=0.9, ls=":", zorder=1)
    am.plot(cs, L, "o-", color="#1d4ed8", lw=1.8, ms=5, zorder=3,
            label=r"$\lambda/\lambda_{\rm obs}$  (seismicity front)")
    am.plot(cs, SL, "s-", color="#D55E00", lw=1.8, ms=5, zorder=3,
            label=r"$\delta/\delta_{\rm obs}$  (total slip)")
    am.set_xscale("log"); am.set_yscale("log")
    am.set(xlabel=r"Permeability contrast  $k_{p,\max}/k_{p,\min}$",
           ylabel="Model / observed")
    am.set_title("(c)  The front is flat; the slip is not")
    am.legend(loc="upper left")
    at = am.twinx()
    at.plot(cs, DP, "^--", color="#009E73", lw=1.5, ms=5, zorder=2)
    at.axhline(0.0, color="#009E73", lw=0.7, ls=":", alpha=0.6)
    at.set_ylabel("Wellhead pressure error (%)", color="#009E73")
    at.tick_params(axis="y", colors="#009E73", labelcolor="#009E73")
    at.spines["right"].set_color("#009E73")

    ar.set(xlabel=f"Days since data-day {T0}", ylabel="Front radius (m)",
           xlim=(0, TMAX), ylim=(0, 1200))
    ar.set_title("(a)  Seismicity front")
    ar.legend(loc="lower right")
    ap_.set(xlabel=f"Days since data-day {T0}",
            ylabel=r"Pressure change from $p_{f0}$ (MPa)",
            xlim=(0, TMAX), ylim=(0, 28))
    ap_.set_title("(b)  Wellhead pressure")
    ap_.legend(loc="upper left")

    sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    # LABEL ABOVE THE BAR, TICKS ON ITS LEFT. With the label on the bar's
    # right it ran into panel (c)'s y-label even with a dedicated column, since
    # a colorbar label extends outside its own axes.
    fig.colorbar(sm, cax=cax)
    cax.yaxis.set_ticks_position("left")
    cax.set_title(r"$k_{p,\max}/k_{p,\min}$", fontsize=8.5, pad=7)

    spread_l = max(L) / min(L)
    spread_s = max(SL) / min(SL)
    print(f"\nover a {max(cs)/min(cs):.0f}-fold contrast range: "
          f"lambda spans {spread_l:.2f}x, total slip spans {spread_s:.2f}x, "
          f"dp spans {max(DP)-min(DP):.0f} points")

    stem = a.out or f"fig_tradeoff_tau{int(round(a.tau*100))}"
    FIG.mkdir(parents=True, exist_ok=True)
    for e_ in ("png", "pdf"):
        fig.savefig(FIG / f"{stem}.{e_}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {FIG}/{stem}.png")


if __name__ == "__main__":
    main()
