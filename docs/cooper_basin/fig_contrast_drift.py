#!/usr/bin/env python3
"""Why the late-time pressure cannot be fixed with permeability.

Three targets, one permeability field, and they are mutually exclusive. This
plots the measurement that shows it, from runs that already existed.

THE SOFTENING METRIC. For a linear well dp/q is the same on every sustained
flow window. Measured on the 4-8 d plateau (q 26.5 L/s) and the 10.2-12.5 d
late window (q 37.8), the record gives 439 -> 431 kPa per L/s, a DRIFT OF -2%:
Habanero 4's injectivity does not change over the 13 days. Any run whose drift
is strongly negative is softening as it goes, and will miss the day-9 step no
matter how well it fits the plateau.

WHAT THE kpmax SWEEP MEASURES. Across contrast kpmax/kpmin from 3x to 400x at
two stress states, the drift is monotonic in contrast and nothing else:

    3x   +2%      65x   -9%
    6x   +1%     145x  -16%
   13x   -4%     400x  -21%
   30x   -8%

Contrast is the amount of room kp has to grow. At 3-6x it is already near
kpmax, so almost no enhancement happens and the well stays linear -- matching
the observed -2%. At 400x the enhanced zone grows through the whole run and the
well softens 21%.

AND THE FRONT IS FINE THERE TOO. 633091 at 6x gives lambda = 0.99x AND drift
+1%: the seismicity front and the well's linearity are both matched, together,
at low contrast. That has not happened in any other run.

SO WHY NOT USE IT. Because the LEVEL is wrong by an order of magnitude -- 633091
sits +1051% on the plateau. With kpmin pinned at 1e-15 through this whole
project, contrast and absolute permeability are the SAME KNOB, so low contrast
necessarily means low permeability and absurd pressure.

AND FREEING kpmin DOES NOT RESCUE IT. Scale every permeability by s: dp goes as
1/s and the front as sqrt(s). Fixing 633091's 11.5x pressure excess needs
s = 11.5, which multiplies the front by 3.4 -- so lambda goes 0.99x to ~3.4x.
The front wants low permeability (D ~ 0.06 m2/s, kp ~ 1.8e-15) and the well's
injectivity wants high (kp ~ 2.5e-13). That 140x gap is what FORCES a large
contrast, and a large contrast is what makes the well soften.

THE CONCLUSION IS STRUCTURAL, not a tuning failure. Front, injectivity level
and injectivity linearity are three constraints on one scalar field, and in a
two-zone model with slip-driven enhancement they cannot be satisfied at once.
The late-time pressure therefore needs something that is not permeability.

633120 is the independent check: evolution off, contrast 142x, drift +27%. With
no growth allowed the well stiffens instead, so the drift is set by how much
enhancement happens, and it brackets the observed -2% between 633001 (-20%) and
633120 (+27%).

Usage:  python fig_contrast_drift.py
"""
import importlib.util as iu
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
H = "/home/users/nberrios/3dhbi/hbi_analysis"
_s = iu.spec_from_file_location("sf", H + "/make_sweep_figures.py")
sf = iu.module_from_spec(_s); _s.loader.exec_module(sf)
_f = iu.spec_from_file_location("fp", HERE / "fault_pressure.py")
fp = iu.module_from_spec(_f); _f.loader.exec_module(fp)
_a = iu.spec_from_file_location("ca", HERE / "compare_arms.py")
ca = iu.module_from_spec(_a); _a.loader.exec_module(ca)

T0 = 4.300
PLAT, LATE = (4.0, 8.0), (10.2, 12.5)
FAM = {11.53: (633090, 633091, 633092, 633093, 633094, 633095, 633096),
       12.71: (633100, 633101, 633102, 633103, 633104, 633105, 633106)}
EXTRA = {633001: ("633001", "#1d4ed8"), 633120: ("633120, permev F", "#009E73")}
FIG = Path(H) / "figures" / "cycle2"
INK, MUTED, OBSC = "#1a1a19", "#6b6b66", "#a8071a"

plt.rcParams.update({"font.size": 10.5, "axes.titlesize": 11.5,
                     "axes.labelsize": 11, "axes.edgecolor": MUTED,
                     "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED,
                     "xtick.labelcolor": INK, "ytick.labelcolor": INK,
                     "legend.fontsize": 9})


def avg(t, y, w):
    m = (np.asarray(t) >= w[0]) & (np.asarray(t) <= w[1])
    return float(np.mean(np.asarray(y)[m])) if m.any() else np.nan


def main(argv=None):
    obs = sf.observed()
    to, dpo = fp.dp_observed(obs)
    ti, q = obs["ti"] - T0, obs["q"]
    PS = sf.P0 - sf.RHO * sf.G * sf.HW / 1e6
    tc, rc, _ = ca.fb.catalogue()
    k = (tc - T0) > 0
    te, re_ = tc[k] - T0, rc[k]
    org = float(np.median(re_[np.argsort(te)][:10]))
    ft, fd = ca.seismicity_front(te, re_ - org)
    LOBS = ca.lam(ft, fd)

    qp, ql = avg(ti, q, PLAT), avg(ti, q, LATE)
    op, ol = avg(to, dpo, PLAT), avg(to, dpo, LATE)
    drift_obs = 100.0 * ((ol / ql) - (op / qp)) / (op / qp)
    print(f"observed dp/q {1000*op/qp:.0f} -> {1000*ol/ql:.0f} kPa per L/s "
          f"({drift_obs:+.0f}%)")

    def measure(n):
        dk = sf.deck(n)
        d = sf.run_data(n, dk)
        if d is None or d.get("tpw") is None or len(d["T"]) < 3:
            return None
        dpm = d["ppw"] - PS
        mp, ml = avg(d["tpw"], dpm, PLAT), avg(d["tpw"], dpm, LATE)
        o = np.argsort(d["T"])
        L = ca.lam(np.asarray(d["T"])[o], np.asarray(d["R"])[o] * 1000.0)
        return dict(c=sf.ffloat(dk["kpmax"]) / sf.ffloat(dk["kpmin"]),
                    drift=100.0 * ((ml / ql) - (mp / qp)) / (mp / qp),
                    lev=100.0 * (mp - op) / op, lam=L / LOBS)

    fig, (aL, aR) = plt.subplots(1, 2, figsize=(13.4, 5.4), dpi=200,
                                 constrained_layout=True)
    cols = {11.53: "#8E44AD", 12.71: "#D55E00"}
    for tau, runs in FAM.items():
        m = [x for x in (measure(n) for n in runs) if x]
        m.sort(key=lambda x: x["c"])
        c = [x["c"] for x in m]
        aL.plot(c, [x["drift"] for x in m], "o-", lw=1.8, ms=5,
                color=cols[tau], label=r"$\tau_0$ = " f"{tau:.2f} MPa")
        aR.plot(c, [abs(x["lev"]) for x in m], "o-", lw=1.8, ms=5,
                color=cols[tau], label=r"$\tau_0$ = " f"{tau:.2f} MPa")
        for x in m:
            print(f"  tau {tau}  contrast {x['c']:6.0f}x  drift {x['drift']:+6.1f}%"
                  f"  level {x['lev']:+9.1f}%  lambda {x['lam']:.2f}")
    for n, (lab, col) in EXTRA.items():
        x = measure(n)
        if x:
            aL.plot([x["c"]], [x["drift"]], "D", ms=8, color=col, label=lab)
            aR.plot([x["c"]], [abs(x["lev"])], "D", ms=8, color=col, label=lab)
            print(f"  {lab}: contrast {x['c']:.0f}x  drift {x['drift']:+.1f}%"
                  f"  level {x['lev']:+.1f}%  lambda {x['lam']:.2f}")

    aL.axhline(drift_obs, color=OBSC, lw=2.2,
               label=f"observed, {drift_obs:+.0f}%")
    aL.axhspan(drift_obs - 3, drift_obs + 3, color=OBSC, alpha=0.12, lw=0)
    aL.set(xscale="log", xlabel=r"permeability contrast $k_{p,max}/k_{p,min}$",
           ylabel="drift in $\\Delta p/q$, plateau $\\to$ late (%)")
    aL.set_title("(a)  The well's linearity is set by the contrast\n"
                 "negative = the model softens as it goes", fontsize=11)
    aL.legend(loc="lower left", framealpha=0.95)
    aR.axhline(10, color=MUTED, lw=1.2, ls="--", label="10% error")
    aR.set(xscale="log", yscale="log",
           xlabel=r"permeability contrast $k_{p,max}/k_{p,min}$",
           ylabel="|plateau $\\Delta p$ error| (%)")
    aR.set_title("(b)  But the pressure LEVEL needs a large contrast\n"
                 "with $k_{p,min}$ pinned at 1e-15 these are one knob",
                 fontsize=11)
    aR.legend(loc="upper right", framealpha=0.95)
    fig.suptitle("Linearity wants a small contrast, level wants a large one. "
                 "One permeability field cannot do both.", fontsize=12.5)

    FIG.mkdir(parents=True, exist_ok=True)
    for e in ("png", "pdf"):
        fig.savefig(FIG / f"contrast_drift.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {FIG}/contrast_drift.png")


if __name__ == "__main__":
    main()
