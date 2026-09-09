#!/usr/bin/env python3
"""Does the deck's injection file reproduce the measured flow rate?

Every result in cycle 2 rests on this file being a faithful shift of the
measured record. Nothing has checked it since it was written, and the check is
worth having on its own terms: if the model is fed the wrong rate history then
the front and the wellhead dp are both wrong for a reason no parameter sweep
can reach.

WHAT IS COMPARED

  measured   Cooper_Basin_HAB_4_Injection_Rate.mat, the field record, read
             exactly as make_sweep_figures.observed() does -- rate x 1000/60 to
             go from the file's units to L/s -- then shifted by -4.300 d so it
             shares the simulation clock.

  deck       june_clean_from_d4300.txt, the file HBI actually opens. Its
             structure is nwell / npoint / times / i j / rates, read at
             m_diffusion.f90:58-66.

UNITS, and this is where a comparison like this usually goes wrong. The deck
holds times in SECONDS and rates in m^2/s -- a rate per unit fault width, not a
volumetric rate. Converting to L/s needs the width:

    q[L/s] = q[file] * W * 1000,    W = 6.0 m

W is not in the deck. It is hard-coded at make_sweep_figures.py:55 and applied
at :284, so the deck file and the field record cannot be compared without it,
and a wrong W would rescale the whole curve while leaving its shape perfect.
The peak checks out: file max 0.007994 m^2/s -> 47.96 L/s against a measured
peak of 47.96 L/s after t0.

WHAT THE THREE PANELS TEST

  (a) rate against time. Tests the SHAPE -- whether the 285-point piecewise
      representation follows a record sampled ~1.3 million times.
  (b) cumulative volume. Tests the INTEGRAL, which is what actually drives the
      pressure field; small rate errors that cancel do not matter here, and a
      persistent bias does.
  (c) residual, deck minus measured. Tests WHERE the disagreement sits rather
      than only how large it is -- a few large spikes at rate steps are
      harmless, a systematic offset over a plateau is not.

THE FILE EXTENDS PAST THE DATA and that is deliberate: it runs to 27.09 d while
the record ends 13.155 d after t0. Exactly one point lies beyond, (27.09 d, 0),
padding so HBI has a bracketing value; injection in the file stops at 12.95 d.
So the simulation sees no injection the data does not have.

Usage:  python compare_injection.py
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

IN = Path("/home/groups/edunham/nberrios/3dhbi/examples/grid_search_inputs")
DECK_Q = "june_clean_from_d4300.txt"
T0 = 4.300
W = sf.W                                  # 6.0 m, from make_sweep_figures.py:55
OUT = Path(H) / "figures" / "cycle2"
INK, MUTED, GRID = "#1a1a19", "#6b6b66", "#d8d8d4"
MEAS, DECK = "#6b6b66", "#1d4ed8"

plt.rcParams.update({"font.size": 10.5, "axes.titlesize": 11.5,
                     "axes.labelsize": 11, "axes.edgecolor": MUTED,
                     "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED,
                     "xtick.labelcolor": INK, "ytick.labelcolor": INK,
                     "legend.fontsize": 9})


def deck_rate(name=DECK_Q):
    """(t_days, q_L_per_s) from the injection file HBI opens."""
    f = (IN / name).read_text().split("\n")
    nwell, npoint = int(f[0]), int(f[1])
    t = np.array(f[2].split(), float)
    q = np.array(f[4].split(), float)
    assert len(t) == npoint and len(q) == npoint, \
        f"{name}: header says npoint = {npoint}, got {len(t)} times, {len(q)} rates"
    assert nwell == 1, f"{name}: nwell = {nwell}, this reader assumes 1"
    return t / 86400.0, q * W * 1000.0


def cumulative(t_days, q_Ls):
    """Trapezoidal cumulative volume in ML. q is L/s and t is days."""
    return np.concatenate([[0.0], np.cumsum(
        np.diff(t_days) * 86400.0 * 0.5 * (q_Ls[1:] + q_Ls[:-1]))]) / 1e6


def main(argv=None):
    obs = sf.observed()
    to = obs["ti"] - T0
    k = to >= 0
    to, qo = to[k], obs["q"][k]
    td, qd = deck_rate()

    Vo, Vd = cumulative(to, qo), cumulative(td, qd)
    hi = min(to.max(), td.max())
    gr = np.linspace(0.0, hi, 20000)
    ro = np.interp(gr, to, qo)
    rd = np.interp(gr, td, qd)
    res = rd - ro

    print(f"deck file {DECK_Q}: {len(td)} points, 0 to {td.max():.4f} d")
    print(f"measured record after t0: {len(to)} samples, 0 to {to.max():.4f} d")
    print(f"  W = {W} m applied to the deck's m^2/s; peaks "
          f"{qd.max():.3f} (deck) vs {qo.max():.3f} L/s (measured)")
    print(f"\nVOLUME to {hi:.3f} d")
    vo = float(np.interp(hi, to, Vo)); vd = float(np.interp(hi, td, Vd))
    print(f"  measured {vo:.4f} ML, deck {vd:.4f} ML, "
          f"difference {100*(vd-vo)/vo:+.4f}%")
    print(f"\nRATE RESIDUAL, deck - measured, on {len(gr)} points")
    print(f"  mean {res.mean():+.4f} L/s, rms {np.sqrt((res**2).mean()):.4f}, "
          f"max |.| {np.abs(res).max():.3f} at {gr[np.argmax(np.abs(res))]:.3f} d")
    for thr in (1.0, 2.0, 5.0):
        f_ = 100.0 * (np.abs(res) > thr).mean()
        print(f"  |residual| > {thr:.0f} L/s on {f_:.2f}% of the window")
    # The total volume error is two unrelated things and quoting only the
    # total hides which one matters. Split at the scoring time.
    print(f"\nWHERE THE VOLUME EXCESS ACCRUES")
    for t in (4.0, 8.0, 8.7, 12.0, 12.8, 12.9, hi):
        a = float(np.interp(t, to, Vo)); b = float(np.interp(t, td, Vd))
        print(f"  {t:6.3f} d   measured {a:8.4f} ML   deck {b:8.4f}   "
              f"{1000*(b-a):+8.1f} kL   {100*(b-a)/a:+6.3f}%")
    nz, nzd = np.where(qo > 0.05)[0], np.where(qd > 0.05)[0]
    print(f"  the jump after 12.8 d is the SHUT-IN TAIL: the measured record "
          f"stops abruptly at\n  {to[nz[-1]]:.4f} d, while the deck's last "
          f"nonzero point is {qd[nzd[-1]]:.2f} L/s at {td[nzd[-1]]:.4f} d "
          f"followed by\n  ({td[nzd[-1]+1]:.4f} d, 0), so HBI ramps down "
          f"linearly over {(td[nzd[-1]+1]-td[nzd[-1]])*24:.2f} h and injects\n"
          f"  ~133 kL that the record does not. It lands AFTER the 8.7 d "
          f"scoring time, where the\n  error is only +0.38%.")
    print(f"\ninjection in the deck stops at "
          f"{td[np.max(np.where(qd > 0))]:.3f} d; last point ({td[-1]:.2f} d, "
          f"{qd[-1]:.1f} L/s) is padding")

    # THREE PANELS, ONE QUESTION EACH, and no twin axes. The previous version
    # put the volume DIFFERENCE on a twin axis inside the volume panel, where a
    # red curve crossed two overlapping curves and read as a third dataset in
    # the same space; and it gave a whole panel to the rate residual, which is
    # a field of spikes at rate steps and looks like noise rather than a
    # result. The difference now gets its own panel, in PERCENT, which is the
    # number anyone actually wants.
    fig, ax = plt.subplots(3, 1, figsize=(12.0, 10.0), dpi=200, sharex=True,
                           gridspec_kw=dict(height_ratios=[1.5, 1.1, 1.1],
                                            hspace=0.28))
    aq, av, ar = ax
    SCORE_T = 8.7

    # (a) the measured record raw and pale, the deck bold over it. The grey
    # hash is the record's own momentary dropouts, ~1.1 M samples of it; the
    # point of the panel is that the blue line follows the BODY of the grey.
    aq.plot(to, qo, lw=1.0, color=MEAS, alpha=0.40,
            label=f"measured, {len(to):,} samples".replace(",", " "))
    aq.plot(td, qd, lw=2.0, color=DECK,
            label=f"{DECK_Q}, {len(td)} points")
    aq.set(ylabel="Injection rate (L/s)", ylim=(0, 55))
    aq.set_title("(a)  Injection rate")
    aq.legend(loc="lower center", framealpha=0.95, ncol=2)

    # (b) the integral. The two curves lying on top of each other IS the
    # result, so nothing is added to make them look different.
    av.plot(to, Vo, lw=3.0, color=MEAS, alpha=0.55, label="measured")
    av.plot(td, Vd, lw=1.5, color=DECK, label=DECK_Q)
    av.set(ylabel="Cumulative volume (ML)")
    av.set_title("(b)  Cumulative volume")
    av.legend(loc="upper left", framealpha=0.95)

    # (c) how wrong, as a percentage of the measured volume to that time.
    pct = 100.0 * (np.interp(gr, td, Vd) - np.interp(gr, to, Vo)) \
        / np.maximum(np.interp(gr, to, Vo), 1e-9)
    ar.axhline(0, color=MUTED, lw=1.0)
    ar.plot(gr, pct, lw=1.8, color="#a8071a")
    ar.axvline(SCORE_T, color=INK, lw=1.4, ls="--")
    # NO IN-PANEL TEXT. The dashed vertical marks the scoring time; what it
    # is worth there, and what the step at 12.85 d is, are printed by main()
    # and belong in whatever the figure is placed into.
    p_at = float(np.interp(SCORE_T, gr, pct))
    ar.set(xlabel="Days since injection resumed (data-day 4.300)",
           ylabel="Volume error (%)", xlim=(0, hi), ylim=(-0.35, 1.25))
    ar.set_title("(c)  Volume error, file relative to measured")

    OUT.mkdir(parents=True, exist_ok=True)
    for e in ("png", "pdf"):
        fig.savefig(OUT / f"injection_check.{e}", bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {OUT}/injection_check.png")


if __name__ == "__main__":
    main()
