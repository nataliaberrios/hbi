#!/usr/bin/env python3
"""One figure per STRESS STATE, coloured by disc radius. The transpose of
compare_arms.py.

compare_arms.py fixes the fluid and colours by tau_0. This fixes tau_0 and
colours by the initial high-permeability disc radius, 150 to 450 m. Same four
panels, same definitions, same datum -- the helpers are imported from
compare_arms rather than reimplemented, so the seismicity front, the lambda fit
and the pressure datum cannot differ between the two views.

WHY BOTH VIEWS ARE NEEDED. The ranking shows two surfaces crossing: lambda/obs
= 1 lies along the tau_0 = 11.53 row, while dp = 0 runs diagonally from
(10.36, 450 m) through (11.53, 400 m) to (12.71, 150 m). compare_arms shows
the first; this shows the second. Read at tau_0 = 11.53 the two coincide, which
is the result -- and it is a crossing of two families, not one lucky run.

Also visible here and nowhere else: THE DISC'S GRIP ON dp WEAKENS AS tau_0
RISES. At tau_0 10.36 the disc spans 29 points of dp (+31 to +2%); at tau_0
15.00 only 3 points (-36 to -39%). More slip means more permeability
enhancement, which swamps the initial disc, so the disc is a lever only at low
stress. On the tau_0 15.00 figure the seven curves in panel (d) lie on top of
each other; on the 10.36 figure they fan across the whole axis.

Runs are DISCOVERED from the decks, not listed: arm-2 fluid (eta 1.27e-4,
beta 2.25e-8, phi 0.01), no skin, the d4300 injection file, and an output
directory that exists. So the phi, Sw_fwid and trickle runs are excluded
automatically rather than by hand.

Usage:
  python compare_disc.py                 # all five stress states
  python compare_disc.py --tau 11.53
"""
import argparse
import glob
import importlib.util as iu
import re
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

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
OUTD = Path("/scratch/users/nberrios/3dhbi/output")
FIG = Path(H) / "figures" / "cycle2"
T0, TMAX = 4.300, 13.2
INK, MUTED, OBSC = "#1a1a19", "#6b6b66", "#a8071a"

plt.rcParams.update({"font.size": 10.5, "axes.titlesize": 11.5,
                     "axes.labelsize": 11, "axes.edgecolor": MUTED,
                     "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED,
                     "xtick.labelcolor": INK, "ytick.labelcolor": INK,
                     "legend.fontsize": 8.5})


def discover():
    """{tau_0: [(disc, run), ...]} for the comparable arm-2 runs."""
    out = {}
    for f in sorted(glob.glob(str(IN / "res6329[0-9][0-9].in"))
                    + glob.glob(str(IN / "res6330[0-9][0-9].in"))):
        n = int(re.search(r"res(\d+)", f).group(1))
        if not (OUTD / str(n)).is_dir():
            continue
        d = sf.deck(n)
        try:
            if abs(sf.ffloat(d["eta"]) - 1.27e-4) > 1e-12: continue
            if abs(sf.ffloat(d["beta"]) - 2.25e-8) > 1e-12: continue
            if abs(sf.ffloat(d["phi"]) - 0.01) > 1e-12: continue
            if abs(sf.ffloat(d["Sw_fwid"]) - 7.4e-9) > 1e-13: continue
            if "skin" in d: continue
            if "d4300" not in d["injection_file"]: continue
            m = re.search(r"disc(\d+)", d["parameter_file"])
            if not m: continue
        except KeyError:
            continue
        tau = round(sf.ffloat(d["muinit"]) * sf.ffloat(d["sigmainit"]), 2)
        out.setdefault(tau, []).append((int(m.group(1)), n))
    return {k: sorted(v) for k, v in sorted(out.items())}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--tau", nargs="+", type=float)
    a = ap.parse_args(argv)

    groups = discover()
    if a.tau:
        groups = {k: v for k, v in groups.items()
                  if any(abs(k - t) < 0.02 for t in a.tau)}
    obs = sf.observed()
    tcat, rcat, _ = ca.fb.catalogue()
    kc = (tcat - T0) > 0
    te, re_ = tcat[kc] - T0, rcat[kc]
    org = float(np.median(re_[np.argsort(te)][:10]))
    ft, fd = ca.seismicity_front(te, re_ - org)
    LOBS = ca.lam(ft, fd)
    tvol, vol = ca.volume_since_t0(obs)
    Vf = np.interp(ft, tvol, vol)
    LOBSV = ca.lam(Vf, fd)
    to, dpo = fp.dp_observed(obs)
    tf = np.linspace(0.0, TMAX, 300)
    vf = np.linspace(0.0, vol.max(), 300)
    print(f"observed lambda {LOBS:.1f} m/sqrt(d), lambda_V {LOBSV:.1f}; "
          f"datum {fp.datum():.3f} MPa")

    for tau, members in groups.items():
        cols = plt.cm.plasma(np.linspace(0.05, 0.85, len(members)))
        fig, ax = plt.subplots(2, 2, figsize=(15.0, 9.6), dpi=200,
                               constrained_layout=True)
        (art, arv), (asz, adp) = ax

        art.scatter(ft, fd + org, s=13, color=OBSC, alpha=0.75, lw=0,
                    label=f"observed seismicity front ({len(ft)} points)")
        art.plot(tf, LOBS*np.sqrt(tf) + org, "-", lw=2.4, color=OBSC,
                 label=r"   $\lambda$ = " f"{LOBS:.1f} m/" r"$\sqrt{d}$")
        arv.scatter(Vf, fd + org, s=13, color=OBSC, alpha=0.75, lw=0,
                    label="observed seismicity front")
        arv.plot(vf, LOBSV*np.sqrt(vf) + org, "-", lw=2.4, color=OBSC,
                 label=r"   $\lambda_V$ = " f"{LOBSV:.1f}")
        asz.scatter(te, re_, s=3.0, alpha=0.18, color=MUTED, lw=0,
                    label=f"{int(kc.sum())} events after $t_0$")
        adp.plot(to, dpo, lw=1.4, color=MUTED, alpha=0.9,
                 label=r"measured $\Delta p$")

        print(f"\ntau_0 = {tau:.2f} MPa")
        print(f"  {'disc':>5} {'run':>7} {'lam/obs':>8} {'dp %':>7}")
        for c, (disc, n) in zip(cols, members):
            d = sf.run_data(n, sf.deck(n))
            if d is None or len(d["T"]) < 3:
                continue
            o = np.argsort(d["T"])
            T, R = np.asarray(d["T"])[o], np.asarray(d["R"])[o]*1000.0
            Vs = np.interp(T, tvol, vol)
            L, LV = ca.lam(T, R), ca.lam(Vs, R)
            gr = np.linspace(0.05, min(d["tpw"][-1], to.max(), 8.7), 2000)
            ps = np.interp(gr, d["tpw"], d["ppw"]) \
                - (sf.P0 - sf.RHO*sf.G*sf.HW/1e6)
            ob = np.interp(gr, to, dpo)
            qg = np.interp(gr, obs["ti"]-T0, obs["q"])
            fl = (qg > 0.25*np.nanmax(obs["q"])) \
                & (np.interp(gr, obs["tp"]-T0, obs["pm"]) > 5.0)
            e = 100.0*float(np.mean(ps[fl]-ob[fl]))/float(np.mean(ob[fl]))
            lab = f"disc {disc} m  ({n})"
            art.scatter(T, R, s=9, color=c, alpha=0.7, lw=0,
                        label=lab + r"   |   $\lambda$ = "
                              f"{L:.1f}  ({L/LOBS:.2f}$\\times$)")
            art.plot(tf, L*np.sqrt(tf), "-", lw=1.4, color=c, alpha=0.9)
            arv.scatter(Vs, R, s=9, color=c, alpha=0.7, lw=0,
                        label=lab + f"   |   {LV/LOBSV:.2f}$\\times$")
            arv.plot(vf, LV*np.sqrt(vf), "-", lw=1.4, color=c, alpha=0.9)
            asz.plot(T, R, lw=1.9, color=c, label=lab)
            adp.plot(d["tpw"], d["ppw"] - (sf.P0 - sf.RHO*sf.G*sf.HW/1e6),
                     lw=1.7, color=c, label=lab + f"   |   {e:+.1f}%")
            print(f"  {disc:>5} {n:>7} {L/LOBS:>8.2f} {e:>+7.1f}")

        art.set(xlabel="Days since injection resumed (data-day 4.300)",
                ylabel="Front radius (m)", xlim=(0, TMAX), ylim=(0, 3000))
        art.set_title("(a)  Front radius vs time")
        art.legend(loc="upper left", framealpha=0.93)
        arv.set(xlabel="Cumulative injected volume since $t_0$ (ML)",
                ylabel="Front radius (m)", ylim=(0, 3000))
        arv.set_title("(b)  Front radius vs injected volume")
        arv.legend(loc="upper left", framealpha=0.93)
        asz.set(xlabel="Days since injection resumed",
                ylabel="Distance from injection point (m)",
                xlim=(0, TMAX), ylim=(0, 3000))
        asz.set_title("(c)  Seismicity, with the same model fronts")
        asz.legend(loc="upper left", framealpha=0.93)
        adp.set(xlabel="Days since injection resumed",
                ylabel="Pressure change from $p_{f0}$ (MPa)",
                xlim=(0, TMAX), ylim=(-2, 25))
        adp.set_title("(d)  Pressure change from the initial fault pressure")
        adp.legend(loc="upper left", framealpha=0.93, ncol=2)
        fig.suptitle(r"$\tau_0$ = " f"{tau:.2f} MPa "
                     r"($\mu_0$ = " f"{tau/27.99:.4f}),  arm 2,  "
                     "coloured by initial disc radius", fontsize=13)

        stem = f"stress{int(round(tau*100))}_compare"
        FIG.mkdir(parents=True, exist_ok=True)
        for e_ in ("png", "pdf"):
            fig.savefig(FIG / f"{stem}.{e_}", bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {FIG}/{stem}.png")


if __name__ == "__main__":
    main()
