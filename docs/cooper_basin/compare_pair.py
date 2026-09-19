#!/usr/bin/env python3
"""Any set of runs on the data-day-4.300 clock, side by side. Four panels.

compare_arms.py fixes the fluid and colours by tau_0; compare_disc.py fixes
tau_0 and colours by disc radius. Both DISCOVER their members from the decks,
which is right for a sweep and useless for "show me these two". This takes run
numbers on the command line and does nothing else, so a specific pair can be
looked at without inventing a sweep to contain them.

Helpers are imported, not reimplemented: seismicity_front and lam from
compare_arms (the project's own percentile-band definition and the least-squares
fit through the origin), dp_observed from fault_pressure (the Holl datum,
31.629 MPa at the wellhead). A figure from this script therefore cannot disagree
with the sweep figures about what the front is or where zero pressure sits.

EVERY RUN MUST START AT DATA-DAY 4.300. That is asserted from the deck's
injection_file rather than assumed -- 632995 starts at 3.5575 and plotting it
here would misalign it by 0.74 d, which is exactly the error that made the
trickle look significant. Off-clock runs go to compare_trickle.py.

FOUR PANELS
  (a) front radius vs time, with lambda fitted per run
  (b) front radius vs cumulative injected volume
  (c) wellhead pressure change on the Holl datum, with the flowing-mask error
  (d) slip at the injector against the observed increment since t0

THE TWO FRONT NUMBERS DISAGREE AND BOTH ARE REPORTED. lambda/lambda_obs is the
slope of the sqrt(t) fit; R_sim/R_obs at 8.7 d is the radius itself. 633001 is
0.99x on the first and 0.78x on the second, which is not a contradiction: the
model grows diffusively at the right rate from a cloud that starts too small.
Quoting only one of them has misled this project twice.

Usage:
  python compare_pair.py 633001 633110
  python compare_pair.py 633001 633110 632961 --out my_name
"""
import argparse
import importlib.util as iu
import re
import sys
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
sys.path.insert(0, H + "/notebooks")
from sim_curves import load_slip

FIG = Path(H) / "figures" / "cycle2"
OBS_SLIP = Path("/home/users/nberrios/3dhbi/hbi/slip_profiles_strike.txt")
OT = [3, 5, 7, 9, 11, 13, 15, 17]
T0, TMAX, T_EVAL = 4.300, 13.2, 8.7
INK, MUTED, OBSC = "#1a1a19", "#6b6b66", "#a8071a"
COLS = ("#1d4ed8", "#009E73", "#D55E00", "#8E44AD", "#0072B2")

plt.rcParams.update({"font.size": 10.5, "axes.titlesize": 11.5,
                     "axes.labelsize": 11, "axes.edgecolor": MUTED,
                     "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED,
                     "xtick.labelcolor": INK, "ytick.labelcolor": INK,
                     "legend.fontsize": 9})


def describe(dk):
    """The keys that actually distinguish cycle-2 runs, for the legend."""
    tau = sf.ffloat(dk["muinit"]) * sf.ffloat(dk["sigmainit"])
    m = re.search(r"disc(\d+)", dk.get("parameter_file", ""))
    kmax = sf.ffloat(dk["kpmax"])
    kmin = sf.ffloat(dk.get("kpmin", "1e-15"))
    return (tau, int(m.group(1)) if m else None, kmax,
            kmax / kmin if kmin > 0 else np.nan)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+", type=int)
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)

    obs = sf.observed()
    tcat, rcat, _ = ca.fb.catalogue()
    k = (tcat - T0) > 0
    te, re_ = tcat[k] - T0, rcat[k]
    org = float(np.median(re_[np.argsort(te)][:10]))
    ft, fd = ca.seismicity_front(te, re_ - org)
    LOBS = ca.lam(ft, fd)
    tvol, vol = ca.volume_since_t0(obs)
    Vf = np.interp(ft, tvol, vol)
    LOBSV = ca.lam(Vf, fd)
    to, dpo = fp.dp_observed(obs)
    o = np.loadtxt(OBS_SLIP)
    peak = [o[:, 1 + i].max() for i in range(8)]
    s0 = float(np.interp(T0, OT, peak))
    P_STATIC = sf.P0 - sf.RHO * sf.G * sf.HW / 1e6
    R_OBS = float(np.interp(T_EVAL, ft, np.maximum.accumulate(fd + org)))

    fig, ax = plt.subplots(2, 2, figsize=(15.0, 9.6), dpi=200,
                           constrained_layout=True)
    (art, arv), (adp, asl) = ax
    tf = np.linspace(0.0, TMAX, 300)
    vf = np.linspace(0.0, vol.max(), 300)

    art.scatter(ft, fd + org, s=13, color=OBSC, alpha=0.75, lw=0,
                label=f"observed seismicity front ({len(ft)} points)")
    art.plot(tf, LOBS * np.sqrt(tf) + org, "-", lw=2.4, color=OBSC,
             label=r"   $\lambda$ = " f"{LOBS:.1f} m/" r"$\sqrt{d}$")
    arv.scatter(Vf, fd + org, s=13, color=OBSC, alpha=0.75, lw=0,
                label="observed seismicity front")
    arv.plot(vf, LOBSV * np.sqrt(vf) + org, "-", lw=2.4, color=OBSC,
             label=r"   $\lambda_V$ = " f"{LOBSV:.1f}")
    adp.plot(to, dpo, lw=1.4, color=MUTED, alpha=0.9,
             label=r"measured $\Delta p$ (Holl datum, "
                   f"{fp.datum():.2f} MPa)")
    tg = np.linspace(0.2, 13.1, 60)
    asl.plot(tg, np.interp(tg + T0, OT, peak) - s0, lw=2.4, color=OBSC,
             label="observed increment since $t_0$")

    print(f"observed lambda {LOBS:.1f} m/sqrt(d), lambda_V {LOBSV:.1f}, "
          f"R({T_EVAL} d) {R_OBS:.0f} m; datum {fp.datum():.3f} MPa")
    print(f"{'run':>7} {'tau_0':>6} {'disc':>5} {'kpmax':>9} {'ratio':>7} "
          f"{'lam/obs':>8} {'R/R_obs':>8} {'dp %':>7} {'slip':>6}")
    for c, n in zip(COLS, a.runs):
        dk = sf.deck(n)
        inj = dk.get("injection_file", "")
        if "d4300" not in inj:
            raise SystemExit(
                f"{n} uses {inj.strip()}, so its clock does not start at "
                f"data-day {T0}. Use compare_trickle.py for off-clock runs.")
        d = sf.run_data(n, dk)
        if d is None or len(d["T"]) < 3:
            print(f"{n:>7}  no usable output"); continue
        tau, disc, kmax, ratio = describe(dk)

        oo = np.argsort(d["T"])
        T, R = np.asarray(d["T"])[oo], np.asarray(d["R"])[oo] * 1000.0
        Vs = np.interp(T, tvol, vol)
        L, LV = ca.lam(T, R), ca.lam(Vs, R)
        Rev = float(np.interp(T_EVAL, T, R)) if T[0] <= T_EVAL <= T[-1] else np.nan

        gr = np.linspace(0.05, min(d["tpw"][-1], to.max(), T_EVAL), 2000)
        ps = np.interp(gr, d["tpw"], d["ppw"]) - P_STATIC
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
        sr = (np.interp(T_EVAL, st, sv)
              / (float(np.interp(T_EVAL + T0, OT, peak)) - s0)) if st else np.nan

        head = (f"{n}  |  " r"$\tau_0$ " f"{tau:.2f} MPa,  disc {disc} m,  "
                r"$k_{p,max}$ " f"{kmax:.2e}")
        art.scatter(T, R, s=9, color=c, alpha=0.7, lw=0,
                    label=head + r"   |   $\lambda$ = "
                          f"{L:.1f} ({L/LOBS:.2f}$\\times$),  "
                          f"R/R$_{{obs}}$ {Rev/R_OBS:.2f}$\\times$")
        art.plot(tf, L * np.sqrt(tf), "-", lw=1.6, color=c, alpha=0.9)
        arv.scatter(Vs, R, s=9, color=c, alpha=0.7, lw=0,
                    label=f"{n}   |   {LV/LOBSV:.2f}$\\times$")
        arv.plot(vf, LV * np.sqrt(vf), "-", lw=1.6, color=c, alpha=0.9)
        adp.plot(d["tpw"], d["ppw"] - P_STATIC, lw=1.9, color=c,
                 label=head + f"   |   {e:+.1f}%")
        asl.plot(st, sv, lw=1.9, color=c,
                 label=f"{n}   |   {sr:.2f}$\\times$")
        print(f"{n:>7} {tau:>6.2f} {disc:>5} {kmax:>9.2e} {ratio:>6.0f}x "
              f"{L/LOBS:>8.2f} {Rev/R_OBS:>8.2f} {e:>+7.1f} {sr:>6.2f}")

    art.set(xlabel=f"Days since data-day {T0}", ylabel="Front radius (m)",
            xlim=(0, TMAX), ylim=(0, 1600))
    art.set_title("(a)  Front radius vs time")
    art.legend(loc="upper left", framealpha=0.93)
    arv.set(xlabel="Cumulative injected volume since $t_0$ (ML)",
            ylabel="Front radius (m)", ylim=(0, 1600))
    arv.set_title("(b)  Front radius vs injected volume")
    arv.legend(loc="upper left", framealpha=0.93)
    adp.axhline(0, color=INK, lw=1.0, ls="--", alpha=0.5)
    adp.set(xlabel=f"Days since data-day {T0}",
            ylabel="Pressure change from $p_{f0}$ (MPa)",
            xlim=(0, TMAX), ylim=(-2, 22))
    adp.set_title("(c)  Wellhead pressure change")
    adp.legend(loc="lower right", framealpha=0.93)
    asl.set(xlabel=f"Days since data-day {T0}",
            ylabel="Slip at the injector (cm)", xlim=(0, TMAX))
    asl.set_title("(d)  Cumulative slip at the injector")
    asl.legend(loc="upper left", framealpha=0.93)

    stem = a.out or ("compare_" + "_".join(str(r) for r in a.runs))
    FIG.mkdir(parents=True, exist_ok=True)
    for e_ in ("png", "pdf"):
        fig.savefig(FIG / f"{stem}.{e_}", bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {FIG}/{stem}.png")


if __name__ == "__main__":
    main()
