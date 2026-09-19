#!/usr/bin/env python3
"""The figure that prompted the datum correction, shown on both datums.

arm2_compare.py isolates tau_0 with everything else held fixed -- one fluid,
one disc radius, one kpmax -- and that is why a systematic pressure offset was
visible in it. Its panel (d) was the figure that prompted the question "what if
the pressure change were referenced to the initial fault pressure?". The later
figures (compare_disc, and the kpmax sweep) add levers that can absorb such an
offset, so it stops being visible in them.

This shows the same five runs, 632960-632964, with panel (d) drawn twice:

  LEFT   the original convention -- dp measured from the first sample of the
         wellhead record, 34.412 MPa at data-day 0
  RIGHT  the datum derived from the field measurements -- Holl & Barton's
         p_f0 = 72.70 MPa at -4100 mAHD, carried to the wellhead over the
         4182.5 m column at the WCR's measured flowing density of 1001 kg/m3,
         with the pipe friction over the real segmented geometry

Nothing about the SIMULATIONS differs between the two panels. The model curves
are identical; only the measured curve moves, by +2.78 MPa at sim t = 0 easing
to +1.53 by day 10 as the q^2 friction term grows. See docs/pressure_datum.tex
for the derivation and docs/figs/cycle2/DATUM.md for the arithmetic.

WHAT THE PAIR SHOWS. On the left the measured plateau sits at ~9.0 MPa and
every model curve overshoots it; the best tau_0 on pressure alone is 13.86,
where the front is 1.49x too fast. On the right the measured plateau is at
~11.8 MPa, the tau_0 11.53 curve lies on it, and that same tau_0 gives
lambda = 1.03x. The two targets stop disagreeing about the stress -- which is
the whole content of the correction, and it is a statement about how the DATA
were reduced, not about the model.

Usage:  python datum_before_after.py
"""
import importlib.util as iu
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = "/home/users/nberrios/3dhbi/hbi_analysis"
_s = iu.spec_from_file_location("sf", H + "/make_sweep_figures.py")
sf = iu.module_from_spec(_s); _s.loader.exec_module(sf)
_a = iu.spec_from_file_location(
    "ca", "/home/users/nberrios/3dhbi/hbi_git/docs/cooper_basin/compare_arms.py")
ca = iu.module_from_spec(_a); _a.loader.exec_module(ca)
_f = iu.spec_from_file_location(
    "fp", "/home/users/nberrios/3dhbi/hbi_git/docs/cooper_basin/fault_pressure.py")
fp = iu.module_from_spec(_f); _f.loader.exec_module(fp)

RUNS = list(range(632960, 632965))
OLD_REF = 34.412          # the record's first sample, at data-day 0
T0, TMAX = 4.300, 13.2
FIG = Path(H) / "figures" / "cycle2"
INK, MUTED = "#1a1a19", "#6b6b66"

plt.rcParams.update({"font.size": 11, "axes.titlesize": 12,
                     "axes.labelsize": 11.5, "axes.edgecolor": MUTED,
                     "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED,
                     "xtick.labelcolor": INK, "ytick.labelcolor": INK,
                     "legend.fontsize": 9})


def main(argv=None):
    obs = sf.observed()
    tcat, rcat, _ = ca.fb.catalogue()
    te = tcat - T0
    k = te > 0
    org = float(np.median(rcat[k][np.argsort(te[k])][:10]))
    ft, fd = ca.seismicity_front(te[k], rcat[k] - org)
    LOBS = ca.lam(ft, fd)

    t_old = obs["tp"] - T0
    dp_old = obs["pm"] - OLD_REF
    t_new, dp_new = fp.dp_observed(obs)
    P_STATIC = sf.P0 - sf.RHO * sf.G * sf.HW / 1e6

    fig, (aL, aR) = plt.subplots(1, 2, figsize=(14.0, 5.8), dpi=200,
                                 sharey=True, constrained_layout=True)
    cols = plt.cm.viridis(np.linspace(0.05, 0.85, len(RUNS)))

    aL.plot(t_old, dp_old, lw=1.4, color=MUTED, alpha=0.9, label="measured")
    aR.plot(t_new, dp_new, lw=1.4, color=MUTED, alpha=0.9, label="measured")

    print(f"observed lambda {LOBS:.1f} m/sqrt(d)")
    print(f"{'run':>7} {'tau_0':>6} {'lam/obs':>8} {'dp OLD %':>9} "
          f"{'dp NEW %':>9}")
    for c, n in zip(cols, RUNS):
        dk = sf.deck(n)
        d = sf.run_data(n, dk)
        if d is None or d.get("tpw") is None:
            continue
        tau = sf.ffloat(dk["muinit"]) * sf.ffloat(dk["sigmainit"])
        o = np.argsort(d["T"])
        L = ca.lam(np.asarray(d["T"])[o], np.asarray(d["R"])[o] * 1000.0)
        dpm = d["ppw"] - P_STATIC          # identical on both panels
        err = {}
        for tag, (to_, dpo_) in (("old", (t_old, dp_old)),
                                 ("new", (t_new, dp_new))):
            gr = np.linspace(0.05, min(d["tpw"][-1], to_.max(), 8.7), 2000)
            ps = np.interp(gr, d["tpw"], dpm)
            ob = np.interp(gr, to_, dpo_)
            qg = np.interp(gr, obs["ti"] - T0, obs["q"])
            fl = (qg > 0.25 * np.nanmax(obs["q"])) \
                & (np.interp(gr, obs["tp"] - T0, obs["pm"]) > 5.0)
            err[tag] = 100.0 * float(np.mean(ps[fl] - ob[fl])) \
                / float(np.mean(ob[fl]))
        lab = r"$\tau_0$ = " f"{tau:.2f} MPa" r"  |  $\lambda$ = " \
              f"{L/LOBS:.2f}" r"$\times$"
        aL.plot(d["tpw"], dpm, lw=1.8, color=c,
                label=lab + f",  {err['old']:+.0f}%")
        aR.plot(d["tpw"], dpm, lw=1.8, color=c,
                label=lab + f",  {err['new']:+.0f}%")
        print(f"{n:>7} {tau:>6.2f} {L/LOBS:>8.2f} {err['old']:>+9.1f} "
              f"{err['new']:>+9.1f}")

    for ax, ttl, sub in (
            (aL, "(a)  Original convention",
             f"$\\Delta p$ from the record's first sample, {OLD_REF:.3f} MPa"),
            (aR, "(b)  Datum from the field measurements",
             f"$\\Delta p$ from $p_{{f0}}$ = {fp.P_F0:.2f} MPa "
             f"($\\equiv$ {fp.datum():.3f} MPa at the wellhead)")):
        ax.axhline(0, color=INK, lw=1.0, ls="--", alpha=0.6)
        ax.set(xlabel="Days since injection resumed (data-day 4.300)",
               xlim=(0, TMAX), ylim=(-3, 22))
        ax.set_title(f"{ttl}\n{sub}", fontsize=11)
        ax.legend(loc="lower right", framealpha=0.95)
    aL.set_ylabel("Pressure change (MPa)")
    fig.suptitle("The same five simulations, two pressure references. "
                 "Model curves are identical; only the measured curve moves.",
                 fontsize=12.5)

    FIG.mkdir(parents=True, exist_ok=True)
    for e in ("png", "pdf"):
        fig.savefig(FIG / f"datum_before_after.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {FIG}/datum_before_after.png")


if __name__ == "__main__":
    main()
