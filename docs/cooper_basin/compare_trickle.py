#!/usr/bin/env python3
"""With and without the 17.8 h trickle before the rate step-up.

Two runs, identical except for where their clock starts and therefore how much
of the record they are given:

    632960   t0 = data-day 4.300, june_clean_from_d4300.txt, tmax 13.10 d
    632995   t0 = data-day 3.558, june_clean_from_d3558.txt, tmax 13.84 d

3.558 d is where the first shut-in actually ENDS -- flow resumes there at
2.0-2.7 L/s and the wellhead responds, going 32.50 to 33.90 MPa between 3.50
and 3.56 d. 4.300 d is where the rate STEPS UP to 7.9 L/s. The project's decks
have called 4.300 "injection resumption" throughout, which is wrong: it is the
step-up, and 632960 therefore discards 17.8 h and 0.1417 ML of real injection.
The same mistake sits at the record's start, where T_ON = 0.501 is the step-up
and flow actually begins at 0.010 d.

THE CLOCKS DIFFER BY 0.7425 d, so nothing can be compared until both are put on
one axis. Everything here is plotted against DAYS SINCE DATA-DAY 4.300, which
means 632995 begins at -0.7425 and 632960 at 0. That choice keeps 632960 -- the
convention all ~100 cycle-2 runs use -- starting at zero, so the figure reads as
"what the trickle adds" rather than as a different experiment.

The observed series are on the same axis, so the measured curves are identical
between the two comparisons; only the model curves shift.

FOUR PANELS
  (a) injection rate as each deck sees it, which is the only input difference
  (b) seismicity front, with lambda fitted on the project's definition
  (c) wellhead pressure change on the Holl datum
  (d) slip at the injector against the observed increment

Usage:  python compare_trickle.py
"""
import importlib.util as iu
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

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
FIG = Path(H) / "figures" / "cycle2"
OBS_SLIP = Path("/home/users/nberrios/3dhbi/hbi/slip_profiles_strike.txt")
OT = [3, 5, 7, 9, 11, 13, 15, 17]
TREF = 4.300                      # the common axis origin, in data-days
# (run, its own t0 in data-days, label, colour, linewidth, alpha)
#
# THE TWO CURVES VERY NEARLY COINCIDE, which is the result -- so the first case
# is drawn thick and translucent and the second thin and solid on top of it.
# Equal linewidths simply hid 632960 under 632995 everywhere after t = 0 and
# the figure read as though only one run had been plotted.
CASES = [(632960, 4.300, "632960 — no trickle, $t_0$ = 4.300 d",
          "#1d4ed8", 4.2, 0.45),
         (632995, 3.5575, "632995 — trickle included, $t_0$ = 3.5575 d",
          "#009E73", 1.6, 1.0)]
INK, MUTED, OBSC = "#1a1a19", "#6b6b66", "#a8071a"
W = sf.W

plt.rcParams.update({"font.size": 10.5, "axes.titlesize": 11.5,
                     "axes.labelsize": 11, "axes.edgecolor": MUTED,
                     "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED,
                     "xtick.labelcolor": INK, "ytick.labelcolor": INK,
                     "legend.fontsize": 9})


def deck_rate(name):
    """(t_days on the file's own clock, q in L/s) from an injection file."""
    f = (IN / name).read_text().split("\n")
    return (np.array(f[2].split(), float) / 86400.0,
            np.array(f[4].split(), float) * W * 1000.0)


def main(argv=None):
    obs = sf.observed()
    tcat, rcat, _ = ca.fb.catalogue()
    te = tcat - TREF
    k = te > -1.0
    org = float(np.median(rcat[k][np.argsort(te[k])][:10]))
    ft, fd = ca.seismicity_front(te[k][te[k] > 0], (rcat[k] - org)[te[k] > 0])
    LOBS = ca.lam(ft, fd)
    to, dpo = fp.dp_observed(obs)          # already on the 4.300 clock
    o = np.loadtxt(OBS_SLIP)
    peak = [o[:, 1 + i].max() for i in range(8)]
    P_STATIC = sf.P0 - sf.RHO * sf.G * sf.HW / 1e6

    fig, ax = plt.subplots(2, 2, figsize=(14.5, 9.0), dpi=200,
                           constrained_layout=True)
    (aq, ar), (ap, asl) = ax
    tf = np.linspace(-0.9, 13.2, 400)
    tpos = tf[tf > 0]

    ar.scatter(ft, fd + org, s=12, color=OBSC, alpha=0.7, lw=0,
               label=f"observed front ({len(ft)} pts)")
    ar.plot(tpos, LOBS * np.sqrt(tpos) + org, "-", lw=2.2, color=OBSC,
            label=r"   $\lambda$ = " f"{LOBS:.1f}")
    ap.plot(to, dpo, lw=1.3, color=MUTED, alpha=0.9, label="measured")
    tg = np.linspace(0.2, 13.1, 60)
    asl.plot(tg, np.interp(tg + TREF, OT, peak)
             - float(np.interp(TREF, OT, peak)), lw=2.2, color=OBSC,
             label="observed increment")

    print(f"observed lambda {LOBS:.1f} m/sqrt(d); axis origin data-day {TREF}")
    print(f"{'run':>7} {'t0':>7} {'shift':>7} {'lam/obs':>8} {'dp %':>7} "
          f"{'slip':>6}")
    for n, t0, lab, col, lw, al in CASES:
        shift = t0 - TREF              # negative for the earlier start
        dk = sf.deck(n)
        d = sf.run_data(n, dk)
        q_t, q_v = deck_rate(dk["injection_file"].strip('"'))
        aq.plot(q_t + shift, q_v, lw=lw, color=col, alpha=al, label=lab)

        # lambda is FITTED ON THE RUN'S OWN CLOCK -- R = lam*sqrt(t) has an
        # origin, and for 632995 that origin is 3.5575, not 4.300. Fitting it on
        # the shifted axis would put sqrt(t) through the wrong zero and inflate
        # lambda. Only the DRAWING is shifted, via sqrt(t_common - shift).
        oo = np.argsort(d["T"])
        T_own = np.asarray(d["T"])[oo]
        R = np.asarray(d["R"])[oo] * 1000.0
        L = ca.lam(T_own, R)
        ar.scatter(T_own + shift, R, s=8 if al == 1.0 else 26, color=col,
                   alpha=0.7 if al == 1.0 else 0.35, lw=0,
                   label=lab.split(" — ")[0]
                   + r"  |  $\lambda$ = " f"{L:.1f} ({L/LOBS:.2f}" r"$\times$)")
        tdraw = tf[tf > shift]
        ar.plot(tdraw, L * np.sqrt(tdraw - shift), "-", lw=lw, color=col,
                alpha=al * 0.9)

        dpm = d["ppw"] - P_STATIC
        tpw = d["tpw"] + shift
        gr = np.linspace(0.05, min(tpw[-1], to.max(), 8.7), 2000)
        ps = np.interp(gr, tpw, dpm)
        ob = np.interp(gr, to, dpo)
        qg = np.interp(gr, obs["ti"] - TREF, obs["q"])
        fl = (qg > 0.25 * np.nanmax(obs["q"])) \
            & (np.interp(gr, obs["tp"] - TREF, obs["pm"]) > 5.0)
        e = 100.0 * float(np.mean(ps[fl] - ob[fl])) / float(np.mean(ob[fl]))
        ap.plot(tpw, dpm, lw=lw, color=col, alpha=al,
                label=lab.split(" — ")[0] + f"  |  {e:+.1f}%")

        sv, st = [], []
        for t_ in np.linspace(0.2, 13.0, 40):
            try:
                _, sl, ta = load_slip(n, t_, how="strike")
                if abs(ta - t_) < 0.3:
                    st.append(t_ + shift); sv.append(sl[0] * 100)
            except Exception:
                pass
        asl.plot(st, sv, lw=lw, color=col, alpha=al,
                 label=lab.split(" — ")[0])
        spk = np.interp(8.7, st, sv) if st else np.nan   # st is already shifted
        print(f"{n:>7} {t0:>7.4f} {shift:>+7.4f} {L/LOBS:>8.2f} {e:>+7.1f} "
              f"{spk/2.598:>6.2f}")

    aq.set(ylabel="Injection rate (L/s)", xlim=(-1.0, 13.2), ylim=(0, 55))
    aq.set_title("(a)  Injection rate as each deck sees it")
    aq.legend(loc="lower right", framealpha=0.95)
    ar.set(ylabel="Front radius (m)", xlim=(-1.0, 13.2), ylim=(0, 1400))
    ar.set_title("(b)  Seismicity front")
    ar.legend(loc="upper left", framealpha=0.95)
    ap.set(xlabel=f"Days since data-day {TREF}",
           ylabel="Pressure change from $p_{f0}$ (MPa)",
           xlim=(-1.0, 13.2), ylim=(-2, 24))
    ap.set_title("(c)  Wellhead pressure change")
    ap.legend(loc="lower right", framealpha=0.95)
    asl.set(xlabel=f"Days since data-day {TREF}",
            ylabel="Slip at the injector (cm)", xlim=(-1.0, 13.2))
    asl.set_title("(d)  Cumulative slip at the injector")
    asl.legend(loc="upper left", framealpha=0.95)
    fig.suptitle("The 17.8 h of 2.0–2.7 L/s flow before the rate step-up. "
                 "Both runs on one axis; 632995 starts 0.74 d earlier.",
                 fontsize=12.5)

    FIG.mkdir(parents=True, exist_ok=True)
    for e_ in ("png", "pdf"):
        fig.savefig(FIG / f"trickle_compare.{e_}", bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {FIG}/trickle_compare.png")


if __name__ == "__main__":
    main()
